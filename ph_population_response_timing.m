function ph_population_response_timing(population,trials,keys)
warning('off','MATLAB:catenate:DimensionMismatch');

keys.PSTH_binwidth=keys.ON.PSTH_binwidth;
keys.gaussian_kernel=keys.ON.PSTH_binwidth;
keys.kernel_type=keys.ON.kernel_type;
keys.n_consecutive_bins_significant=1; %%!

%% tuning table preparation and grouping
[TT,idx,group_values,unique_group_values]=ph_readout_tuning_table(keys);
if isempty(unique_group_values)
    return;
end
complete_unit_list={population.unit_ID}';
[unit_valid,TM]=ismember(complete_unit_list,TT(:,idx.unitID));
population=population(unit_valid);
complete_unit_list={population.unit_ID}';
population_group=group_values(TM(unit_valid));
[trials.accepted]=deal(true);
[UC, CM]=ph_get_condition_matrix(trials,keys);

%% Convert to ipsi/contra, Baseline subtraction, normalization, re-ordering, and gaussian RFs
[~, ~, condition]=ph_condition_normalization(population,trials,keys,UC,CM);

%condition_matrix            = combvec(CM',UC.hemifield)';
conditions_out              = combvec(UC.effector,CM')';
conditions_hf               = combvec(UC.hemifield,conditions_out')';

%% condition comparison
comparisons=keys.ON.comparisons_per_effector;
condition_parameters_comparissons = [{'hemifield'} {'effector'} keys.condition_parameters];

% comparisons_per_effector is misleading, because effector comparison is possible as well
for t=1:size(condition,1)
    typ=UC.type(t); %UC.type(mod(t-1,numel(UC.type))+1);
    current=[condition(t,:).per_hemifield];
    current_unit=vertcat(current.unit);
    current_window=vertcat(current.window);
    for comp=1:numel(comparisons)
        
        %% condition selection
        cM1=true(size(conditions_hf));
        cM2=true(size(conditions_hf));
        for k=1:numel(condition_parameters_comparissons)
            if isfield(comparisons(comp),condition_parameters_comparissons{k}) &&...
                    ~ isempty(comparisons(comp).(condition_parameters_comparissons{k}))
                cM1(~ismember(conditions_hf(:,k),comparisons(comp).(condition_parameters_comparissons{k}){1}),k)=false;
                cM2(~ismember(conditions_hf(:,k),comparisons(comp).(condition_parameters_comparissons{k}){2}),k)=false;
            end
        end
        c1=find(all(cM1,2));
        c2=find(all(cM2,2));
        
        keys=ph_get_epoch_keys(keys,typ,UC.effector,sum(UC.type_effector(:,1)==typ)>1);%% does it make sense to distinguish by effector?
        wo=find(ismember(keys.PSTH_WINDOWS(:,1),comparisons(comp).order_onset{1})):size(keys.PSTH_WINDOWS,1);
        
        tuning_direction=NaN(size(complete_unit_list));
        tuning_onset=NaN(size(complete_unit_list));
        for u=1:numel(complete_unit_list)
            %% baseline definition (for epoch tuning so far...)
            epoch_FRs=vertcat(current_unit(unique([c1; c2]),u).epoch_FRs);
            baseline=epoch_FRs(:,ismember(keys.EPOCHS(:,1),comparisons(comp).baseline_epoch));
            onset_found=0;
            for wn=1:numel(current(1).window)
                direction=ph_compare_by_bin_by_trial(current_window(:,wn),c1,c2,baseline,u,keys); %% has additional bins in the begininning and end
                direction=direction(keys.n_consecutive_bins_significant:end-(keys.n_consecutive_bins_significant-1)); % return to original window size
                sigbins(t).comparison(comp).window(wn).bins(u,:)=direction;
                
                %% tuning_onset (once per unit, (certain time before/after), until whenever), also defines order onset
                if wn >= wo(1)  %% could be better...
                    if wn==wo(1)
                        n_bins_tooearly=round((comparisons(comp).order_onset{2}-keys.PSTH_WINDOWS{wn,3})/keys.PSTH_binwidth);
                    else
                        n_bins_tooearly=0;
                    end
                    onset=find(~isnan(direction(n_bins_tooearly+1:end)),1);
                    if ~isempty(onset) && onset==1 && n_bins_tooearly>0 % if tuning was already there in the first bin, go backwards !!
                        onset=find(isnan([direction(n_bins_tooearly:-1:1) NaN]),1)*-1+2; %% had +2 here...
                    end
                    onset=onset+n_bins_tooearly;
                end
                if wn < wo(1) || isempty(onset) || onset_found % only 1 onset per unit, assuming earlier windows come earlier in time!!
                    sigbins(t).comparison(comp).window(wn).tuning_onset(u)=NaN;
                    sigbins(t).comparison(comp).window(wn).tuning_direction(u)=NaN;
                else
                    tuning_direction(u)=direction(onset);
                    sigbins(t).comparison(comp).window(wn).tuning_direction(u)=tuning_direction(u);
                    onset= onset+round(keys.PSTH_WINDOWS{wn,3}/keys.PSTH_binwidth);
                    sigbins(t).comparison(comp).window(wn).tuning_onset(u)=onset;
                    onset_found=1;
                    tuning_onset(u)=onset+1000*wn;
                end
            end
        end
        
        %% unit_order dependent on tuning onset
        % Think about ordering by effect size!
        [~, unit_order]=sort(tuning_onset);
        sigbins(t).comparison(comp).order_by_onset=unit_order;
        
        %% always produce other orderings
        c1=tuning_direction(unit_order)==-1;
        c2=tuning_direction(unit_order)==1;
        c3=~c1 & ~c2;
        unit_order=[unit_order(c1); unit_order(c2); unit_order(c3)];
        
        sigbins(t).comparison(comp).order_by_condition_onset=unit_order;
    end
end

%% plots
for t=1:numel(sigbins)
    typ=UC.type(mod(t-1,numel(UC.type))+1);
    for g=1:numel(unique_group_values)
        
        %% here check if N_total makes sense
        g_idx=ismember(population_group,unique_group_values{g});
        N_total=sum(g_idx);
        % N_total=sum(ismember(population_group,unique_group_values{g}));
        current=sigbins(t);
        [fig_title,filename]=ph_figure_titles(keys);
        
        %% save metadata
        unit_IDs=complete_unit_list(g_idx);
        comparison=current.comparison;
        save([keys.basepath_to_save, keys.project_version, filesep, 'population_meta_data', filesep, 'response timing', filesep, filename], 'keys','TT','unit_IDs','comparison');
        
        
        %% tuning per bin plot
        %% different ordering at once!
        order_onsets={'order_by_condition_onset','order_by_onset'};
        for k=1:numel(order_onsets)
            
            plot_1_title            = [fig_title   order_onsets{k} ' per bin'];
            %plot_1_title            = [fig_title  makestringifnumber(unique_group_values{g}) ' per bin'];
            %function makestringifnumber
            PSTH_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',plot_1_title);
            column=1;
            keys=ph_get_epoch_keys(keys,typ,UC.effector,sum(UC.type_effector(:,1)==typ)>1);%% does it make sense to distinguish by effector?
            n=1;
            for comp=1:numel(comparisons)
                wo=find(ismember(keys.PSTH_WINDOWS(:,1),comparisons(comp).order_onset{1}));
                subplot(numel(comparisons),1,(comp-1)+column)
                hold on
                col=comparisons(comp).colors;
                unit_order=current(n).comparison(comp).(order_onsets{k});
                state_shift=0;
                for w=1:size(keys.PSTH_WINDOWS,1)
                    t_before_state=keys.PSTH_WINDOWS{w,3};
                    t_after_state=keys.PSTH_WINDOWS{w,4};
                    bins=t_before_state:keys.PSTH_binwidth:t_after_state;
                    bins=bins+state_shift-t_before_state;
                    to_plot1=current(n).comparison(comp).window(w).bins(unit_order(g_idx),:);
                    to_plot1(isnan(to_plot1))=0;
                    to_plot2=to_plot1;
                    to_plot1(to_plot1==-1)=0;
                    to_plot2(to_plot2==1)=0;
                    to_plot2(to_plot2==-1)=1;
                    to_plot_white=to_plot1==0&to_plot2==0;
                    clear C
                    C(:,:,1) = to_plot1*col(1,1)+to_plot2*col(2,1)+to_plot_white;
                    C(:,:,2) = to_plot1*col(1,2)+to_plot2*col(2,2)+to_plot_white;
                    C(:,:,3) = to_plot1*col(1,3)+to_plot2*col(2,3)+to_plot_white;
                    image([state_shift state_shift+t_after_state-t_before_state],[0 size(to_plot1,1)],C);
                    if w==wo(1)
                        %make line for onset
                        x1=bins(1)+keys.ON.comparisons_per_effector(comp).order_onset{2}-t_before_state;
                        onset_handle(comp)=plot([x1 x1],[0 size(to_plot1,1)],'color',[218 165 32]/255,'linewidth',2);
                    end
                    state_shift=state_shift+t_after_state-t_before_state+0.1;
                end
                title([comparisons(comp).title ' Aligned to: ' comparisons(comp).order_onset{1} ', BL: ' comparisons(comp).baseline_epoch]);
                y_lim(comp)=size(to_plot1,1);
            end
            
            
            % subplot appearance, and tuning lines
            ylimmax=max(max(y_lim));
            ylimmin=0;
            y_lim=[ylimmin ylimmax];
            column=1;
            for comp=1:numel(comparisons)
                hold on
                subplot(numel(comparisons),1,(comp-1)+column)
                
                %% choices? hands? errm here its only relevant for event and epoch onsets
                tr=[trials.type]==typ & ismember([trials.effector],UC.effector) & [trials.completed]==1;
                ph_PSTH_background(trials(tr),y_lim,y_lim,y_lim,keys,1)
                % set ydata
                uistack(onset_handle(comp), 'top');
            end
            ph_title_and_save(PSTH_summary_handle,  [filename ' ' order_onsets{k} ' PSTHs'],plot_1_title,keys)
        end
        
        %% n tuned cells plot
        plot_3_title            = [fig_title  ' n tuned cells'];
        PSTH_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',plot_3_title);
        column=1;
        keys=ph_get_epoch_keys(keys,typ,UC.effector,sum(UC.type_effector(:,1)==typ)>1);%% does it make sense to distinguish by effector?
        
        clear x_lim y_lim sp
        for comp=1:numel(comparisons)
            subplot(numel(comparisons),1,(comp-1)+column)
            hold on
            col=comparisons(comp).colors;
            state_shift=0;
            
            for w=1:size(keys.PSTH_WINDOWS,1)
                concatinated_after_onset=current.comparison(comp).window(w).bins(g_idx,:);
                n_tuned_cells=[sum(concatinated_after_onset==1,1); sum(concatinated_after_onset==-1,1)];
                
                
                
                t_before_state=keys.PSTH_WINDOWS{w,3};
                t_after_state=keys.PSTH_WINDOWS{w,4};
                bins=t_before_state:keys.PSTH_binwidth:t_after_state;
                bins=bins+state_shift-t_before_state;
                to_plot1=n_tuned_cells(1,:)/N_total*100;
                to_plot2=n_tuned_cells(2,:)/N_total*100;
                plot(bins,to_plot1,'color',col(1,:),'linewidth',1);
                plot(bins,to_plot2,'color',col(2,:),'linewidth',1);
                state_shift=state_shift+t_after_state-t_before_state+0.1;
                
                
                maxbin1=nanmean(find(to_plot1==max(to_plot1)))-1;   maxbin2=nanmean(find(to_plot2==max(to_plot2)))-1;
                bins_before=t_before_state/keys.PSTH_binwidth;
                max1=maxbin1+ bins_before; max2=maxbin2 + bins_before;
                
                x1=bins(1)+maxbin1*keys.PSTH_binwidth;
                plot([x1 x1],[0 max(to_plot1)],'color',col(1,:),'linewidth',1);
                texttoplot=['Max=' num2str(max(to_plot1)) '/' b2t(max1,keys)];
                text(x1,max(to_plot1),texttoplot,'color',col(1,:));
                
                x2=bins(1)+maxbin2*keys.PSTH_binwidth;
                plot([x2 x2],[0 max(to_plot2)],'color',col(2,:),'linewidth',1);
                texttoplot=['Max=' num2str(max(to_plot2)) '/' b2t(max2,keys)];
                text(x2,max(to_plot2),texttoplot,'color',col(2,:));
            end
            title([comparisons(comp).title ' Aligned to: ' comparisons(comp).order_onset{1} ', BL: ' comparisons(comp).baseline_epoch]);
            y_lim(comp)=diff(get(gca,'ylim'));
        end
        
        %% subplot appearance, and tuning lines
        if keys.ON.link_y_lim  % key for autscale
            y_lim(1:numel(comparisons))=max(max(y_lim));
        end
        column=1;
        for comp=1:numel(comparisons)
            hold on
            subplot(numel(comparisons),1,(comp-1)+column)
            lims=[0 y_lim(comp)];
            %% choices? hands? errm here its only relevant for event and epoch onsets
            tr=[trials.type]==typ & ismember([trials.effector],UC.effector) & [trials.completed]==1;
            ph_PSTH_background(trials(tr),lims,lims,lims,keys,1)
        end
        
        ph_title_and_save(PSTH_summary_handle,  [filename ' fraction tuned'],plot_3_title,keys)
        
        
        %% tuning onset plot cumsum
        plot_2_title            = [fig_title  ' tuning onset cumsum'];
        PSTH_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',plot_2_title);
        column=1;
        keys=ph_get_epoch_keys(keys,typ,UC.effector,sum(UC.type_effector(:,1)==typ)>1);%% does it make sense to distinguish by effector?
        
        n=1;% find(UC.effector==eff);
        clear x_lim sp
        for comp=1:numel(comparisons)
            
            wo=find(ismember(keys.PSTH_WINDOWS(:,1),comparisons(comp).order_onset{1}));
            
            subplot(numel(comparisons),1,(comp-1)+column)
            hold on
            col=comparisons(comp).colors;
            %             unit_order=current(n).comparison(comp).order_by_onset;
            %             N_total=numel(unit_order); % for getting fraction
            state_shift=0;
            to_plot1_before=0;
            to_plot2_before=0;
            for w=1:size(keys.PSTH_WINDOWS,1)
                tuning_direction=current(n).comparison(comp).window(w).tuning_direction(g_idx);
                tuning_onset=current(n).comparison(comp).window(w).tuning_onset(g_idx);
                
                
                t_before_state=keys.PSTH_WINDOWS{w,3};
                t_after_state=keys.PSTH_WINDOWS{w,4};
                bins=t_before_state:keys.PSTH_binwidth:t_after_state;
                bins=bins+state_shift-t_before_state;
                onset1=tuning_onset(tuning_direction==1);
                onset2=tuning_onset(tuning_direction==-1);
                to_plot1=hist(onset1-t_before_state/keys.PSTH_binwidth,1:numel(bins))/N_total*100;
                to_plot2=hist(onset2-t_before_state/keys.PSTH_binwidth,1:numel(bins))/N_total*100;
                to_plot11=to_plot1;to_plot11(1)=to_plot11(1)+to_plot1_before;
                to_plot22=to_plot2;to_plot22(1)=+to_plot22(1)+to_plot2_before;
                to_plot1_before=sum(to_plot11);
                to_plot2_before=sum(to_plot22);
                plot(bins,cumsum(to_plot11),'color',col(1,:),'linewidth',1);
                plot(bins,cumsum(to_plot22),'color',col(2,:),'linewidth',1);
                state_shift=state_shift+t_after_state-t_before_state+0.1;
                if w==wo(1)
                    maxbin1=nanmean(find(to_plot1==max(to_plot1)))-1;   maxbin2=nanmean(find(to_plot2==max(to_plot2)))-1;
                    bins_before=t_before_state/keys.PSTH_binwidth;
                    max1=maxbin1+ bins_before; max2=maxbin2 + bins_before;
                    mean1=nanmean(onset1)-1; mean2=nanmean(onset2)-1;
                    median1=nanmedian(onset1)-1; median2=nanmedian(onset2)-1;
                    std1=nanstd(to_plot1); std2=nanstd(to_plot2);
                    
                    x1=bins(1)+maxbin1*keys.PSTH_binwidth;
                    plot([x1 x1],[0 max(to_plot1)],'color',col(1,:),'linewidth',1);
                    texttoplot=['Max=' num2str(max(to_plot1)) '/' b2t(max1,keys)  ', u=' b2t(mean1,keys)  ' std=' b2t(std1,keys) ' med=' b2t(median1,keys)];
                    text(x1,max(to_plot1),texttoplot,'color',col(1,:));
                    
                    x2=bins(1)+maxbin2*keys.PSTH_binwidth;
                    plot([x2 x2],[0 max(to_plot2)],'color',col(2,:),'linewidth',1);
                    texttoplot=['Max=' num2str(max(to_plot2)) '/' b2t(max2,keys)  ', u=' b2t(mean2,keys)  ' std=' b2t(std2,keys) ' med=' b2t(median2,keys)];
                    text(x2,max(to_plot2),texttoplot,'color',col(2,:));
                    
                    %make line for onset
                    x1=bins(1)+keys.ON.comparisons_per_effector(comp).order_onset{2}-t_before_state;
                    onset_handle(comp)=plot([x1 x1],[0 max([to_plot1 to_plot2])],'color',[218 165 32]/255,'linewidth',2);
                    
                end
            end
            title([comparisons(comp).title ' Aligned to: ' comparisons(comp).order_onset{1} ', BL: ' comparisons(comp).baseline_epoch]);
            y_lim(comp)=diff(get(gca,'ylim'));
        end
        
        %% subplot appearance, and tuning lines
        if keys.ON.link_y_lim % key for autscale
            y_lim(1:numel(comparisons))=max(max(y_lim));
        end
        column=1;
        for comp=1:numel(comparisons)
            hold on
            subplot(numel(comparisons),1,(comp-1)+column)
            lims=[0 y_lim(comp)];
            %% choices? hands? errm here its only relevant for event and epoch onsets
            tr=[trials.type]==typ & ismember([trials.effector],UC.effector) & [trials.completed]==1;
            % set ydata
            ph_PSTH_background(trials(tr),lims,lims,lims,keys,1)
            uistack(onset_handle(comp), 'top');
        end
        
        ph_title_and_save(PSTH_summary_handle,  [filename ' tuning onset cumsum'],plot_2_title,keys)
        
        
        
        %% tuning onset plot
        plot_2_title            = [fig_title  ' tuning onset'];
        PSTH_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',plot_2_title);
        column=1;
        keys=ph_get_epoch_keys(keys,typ,UC.effector,sum(UC.type_effector(:,1)==typ)>1);%% does it make sense to distinguish by effector?
        
        n=1;% find(UC.effector==eff);
        clear x_lim sp
        for comp=1:numel(comparisons)
            
            wo=find(ismember(keys.PSTH_WINDOWS(:,1),comparisons(comp).order_onset{1}));
            
            subplot(numel(comparisons),1,(comp-1)+column)
            hold on
            col=comparisons(comp).colors;
            state_shift=0;
            for w=1:size(keys.PSTH_WINDOWS,1)
                tuning_direction=current(n).comparison(comp).window(w).tuning_direction(g_idx);
                tuning_onset=current(n).comparison(comp).window(w).tuning_onset(g_idx);
                
                
                t_before_state=keys.PSTH_WINDOWS{w,3};
                t_after_state=keys.PSTH_WINDOWS{w,4};
                bins=t_before_state:keys.PSTH_binwidth:t_after_state;
                bins=bins+state_shift-t_before_state;
                onset1=tuning_onset(tuning_direction==1);
                onset2=tuning_onset(tuning_direction==-1);
                to_plot1=hist(onset1-t_before_state/keys.PSTH_binwidth,1:numel(bins))/N_total*100;
                to_plot2=hist(onset2-t_before_state/keys.PSTH_binwidth,1:numel(bins))/N_total*100;
                plot(bins,to_plot1,'color',col(1,:),'linewidth',1);
                plot(bins,to_plot2,'color',col(2,:),'linewidth',1);
                state_shift=state_shift+t_after_state-t_before_state+0.1;
                if w==wo(1)
                    maxbin1=nanmean(find(to_plot1==max(to_plot1)))-1;   maxbin2=nanmean(find(to_plot2==max(to_plot2)))-1;
                    bins_before=t_before_state/keys.PSTH_binwidth;
                    max1=maxbin1+ bins_before; max2=maxbin2 + bins_before;
                    mean1=nanmean(onset1)-1; mean2=nanmean(onset2)-1;
                    median1=nanmedian(onset1)-1; median2=nanmedian(onset2)-1;
                    std1=nanstd(to_plot1); std2=nanstd(to_plot2);
                    
                    x1=bins(1)+maxbin1*keys.PSTH_binwidth;
                    plot([x1 x1],[0 max(to_plot1)],'color',col(1,:),'linewidth',1);
                    texttoplot=['Max=' num2str(max(to_plot1)) '/' b2t(max1,keys)  ', u=' b2t(mean1,keys)  ' std=' b2t(std1,keys) ' med=' b2t(median1,keys)];
                    text(x1,max(to_plot1),texttoplot,'color',col(1,:));
                    
                    x2=bins(1)+maxbin2*keys.PSTH_binwidth;
                    plot([x2 x2],[0 max(to_plot2)],'color',col(2,:),'linewidth',1);
                    texttoplot=['Max=' num2str(max(to_plot2)) '/' b2t(max2,keys)  ', u=' b2t(mean2,keys)  ' std=' b2t(std2,keys) ' med=' b2t(median2,keys)];
                    text(x2,max(to_plot2),texttoplot,'color',col(2,:));
                    
                    %make line for onset
                    x1=bins(1)+keys.ON.comparisons_per_effector(comp).order_onset{2}-t_before_state;
                    onset_handle(comp)=plot([x1 x1],[0 max([to_plot1 to_plot2])],'color',[218 165 32]/255,'linewidth',2);
                    
                end
            end
            title([comparisons(comp).title ' Aligned to: ' comparisons(comp).order_onset{1} ', BL: ' comparisons(comp).baseline_epoch]);
            y_lim(comp)=diff(get(gca,'ylim'));
        end
        
        %% subplot appearance, and tuning lines
        if keys.ON.link_y_lim % key for autscale
            y_lim(1:numel(comparisons))=max(max(y_lim));
        end
        column=1;
        for comp=1:numel(comparisons)
            hold on
            subplot(numel(comparisons),1,(comp-1)+column)
            lims=[0 y_lim(comp)];
            %% choices? hands? errm here its only relevant for event and epoch onsets
            tr=[trials.type]==typ & ismember([trials.effector],UC.effector) & [trials.completed]==1;
            % set ydata
            ph_PSTH_background(trials(tr),lims,lims,lims,keys,1)
            uistack(onset_handle(comp), 'top');
        end
        
        ph_title_and_save(PSTH_summary_handle,  [filename ' tuning onset'],plot_2_title,keys)
    end
end
end

function hn=ph_compare_by_bin_by_trial(in,c1,c2,baseline,u,keys)
C1=[];
C2=[];
for cc=c1(:)'
    C1=[C1;vertcat(in(cc).unit(u).average_spike_density)];
end
if all(c1==c2) %% epoch comparison (only possibly desired comparison if conditions are exactly identical)
    C2=baseline*ones(1,size(C1,2));
    ho=do_stats(C1,C2,keys,1);
else
    for cc=c2(:)'
        C2=[C2;vertcat(in(cc).unit(u).average_spike_density)];
    end
    ho=do_stats(C1,C2,keys,0);
end

% not sure why there are nan values though (STS memory saccades)
ho(isnan(ho))=0;

%% keep only if n_consecutive_bins are significant
starts=find(diff([0 ho])==1);
ends=find(diff([ho 0])==-1);
hn=NaN(size(ho));
for x=1:numel(starts)
    if ends(x) >= starts(x)+keys.n_consecutive_bins_significant-1
        hn(starts(x):ends(x))=1;
    end
end

%% negative for condition 2 larger
hn(hn==1 & nanmean(C2)>nanmean(C1))=-1;
end

function texttoplot=b2t(mean1,keys)
texttoplot=num2str(round((mean1*keys.PSTH_binwidth)*1000));
end

function h=do_stats(Amat,Bmat,keys,paired)
%h = ttest2(A,B);
%ho=ttest2(A,B,0.05,'both','equal',1);
% with baseline per trial we can always do paired stats!

if keys.ON.permutation_tests
    
    n_permutations=10000;
    [clusts, p_values, ~, ~] = permutest( Amat', Bmat', paired, 0.05, 1000, true);
    indexes_sig=[clusts{p_values<0.05}]; %% p_values<??  --> what is p_crit?
    h=false(1,size(Amat,2));
    h(indexes_sig)=true;
    
else
    for k=1:size(Amat,2)
        A=Amat(:,k);
        B=Bmat(:,k);
        if paired
            if any(~isnan(A)&~isnan(B))
                h(k) = ttest(A,B);
            else
                h(k)=0;
            end
        else
            if any(~isnan(A)) && any (~isnan(B))
                h(k) = ttest2(A,B);
            else
                h(k)=0;
            end
        end
        
        %         if paired
        %             if any(~isnan(A)&~isnan(B))
        %                 [~, h(k)] = signrank(A,B);
        %             else
        %                 h(k)=0;
        %             end
        %         else
        %             if any(~isnan(A)) && any (~isnan(B))
        %                 [~, h(k)] = ranksum(A,B);
        %             else
        %                 h(k)=0;
        %             end
        %         end
    end
end
end