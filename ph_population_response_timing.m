function ph_population_response_timing(population,modified_keys)
warning('off','MATLAB:catenate:DimensionMismatch');
%%% !!!!!!!! make sure epochs (baseline, cue, .... make sense for all effectors)

for fn=fieldnames(modified_keys)'
    keys.(fn{:})=modified_keys.(fn{:});
end
keys.n_consecutive_bins_significant=1; %%!
keys.PSTH_binwidth=keys.ON.PSTH_binwidth;
keys.gaussian_kernel=keys.ON.gaussian_kernel;
keys.kernel_type                        ='gaussian'; %'box'; % we take box kernel here!!

%% tuning table preparation and grouping
keys.normalization_field='ON';
[tuning_per_unit_table]                 = ph_load_extended_tuning_table(keys);
[tuning_per_unit_table, Sel_for_title]  = ph_reduce_tuning_table(tuning_per_unit_table,keys);
idx_group_parameter=DAG_find_column_index(tuning_per_unit_table,keys.ON.group_parameter);
idx_unitID=DAG_find_column_index(tuning_per_unit_table,'unit_ID');
group_values=tuning_per_unit_table(:,idx_group_parameter);
group_values=cellfun(@num2str, group_values, 'UniformOutput', false);
cell_in_any_group=[false; ~ismember(group_values(2:end),keys.ON.group_excluded)];
unique_group_values=unique(group_values(cell_in_any_group));
if isempty(unique_group_values)
    disp('no relevant groups found');
    return;
end
tuning_per_unit_table=tuning_per_unit_table(cell_in_any_group,:);
group_values=group_values(cell_in_any_group);
complete_unit_list={population.unit_ID}';
population=population(ismember(complete_unit_list,tuning_per_unit_table(:,idx_unitID)));

all_trialz=[population.trial];
[UC, CM]=ph_get_condition_matrix(all_trialz,keys);

% reduce trials to only valid ... obsolete?
unit_valid=true(size(population));
for u=1:numel(population)
    poptr=population(u).trial;
    valid=ismember([poptr.effector],UC.effector) & ismember([poptr.type],UC.type);% & ismember([poptr.choice],UC.choice) & ismember([poptr.reach_hand],UC.reach_hand);
    population(u).trial=population(u).trial(valid);
    if sum(valid)==0
        unit_valid(u)=false;
    end
end
population=population(unit_valid);
complete_unit_list={population.unit_ID}';
unit_valid=ismember(tuning_per_unit_table(:,idx_unitID),complete_unit_list);
group_values=group_values(unit_valid);
tuning_per_unit_table=tuning_per_unit_table(unit_valid,:);

%% Convert to ipsi/contra, Baseline subtraction, normalization, re-ordering, and gaussian RFs
[~, ~, condition]=ph_condition_normalization(population,keys,UC,CM);

%condition_matrix            = combvec(CM',UC.hemifield)';
conditions_out              = combvec(UC.effector,CM')';
conditions_hf               = combvec(UC.hemifield,conditions_out')';

%% condition comparison
comparisons=keys.ON.comparisons_per_effector;
condition_parameters_comparissons = [{'hemifield'} {'effector'} keys.condition_parameters];

% comparisons_per_effector is misleading, because effector comparison is possible as well
for g=1:numel(unique_group_values)
    for t=1:size(condition,1)
        typ=UC.type(t); %UC.type(mod(t-1,numel(UC.type))+1);
        current=[condition(t,:).per_hemifield];
        current_unit=vertcat(current.unit);
        current_window=vertcat(current.window);
        %     why do we need this at all????
        %     idx_existing=DAG_find_column_index(tuning_per_unit_table,['existing_' 'in_AH_' type_effector_short{tye} '_' keys.arrangement(1:3)]);
        %     tya_existing{t}= [true; ismember(vertcat(tuning_per_unit_table{2:end,idx_existing}),true)];
        unitidx=ismember(complete_unit_list,tuning_per_unit_table(ismember(group_values,unique_group_values(g)),idx_unitID));
        units=find(all(unitidx,2));
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
            
            for u=1:numel(units)
                %% baseline definition (for epoch tuning so far...)
                epoch_averages=vertcat(current_unit(unique([c1; c2]),units(u)).epoch_FRs);
                baseline=epoch_averages(:,ismember(keys.EPOCHS(:,1),comparisons(comp).baseline_epoch));
                onset_found=0;
                for wn=1:numel(current(1).window)
                    temp_sigbins=ph_compare_by_bin_by_trial(current_window(:,wn),c1,c2,baseline,units(u),keys);
                    sigbins(t).per_group(g).comparison(comp).window(wn).bins(u,:)=temp_sigbins(keys.n_consecutive_bins_significant+1:end-keys.n_consecutive_bins_significant); % return to original window size
                    
                    %% tuning_onset (once per unit, (certain time before/after), until whenever), also defines order onset
                    if wn >= wo(1)  %% could be better...
                        direction=temp_sigbins(1:end-keys.n_consecutive_bins_significant); % keep the first few bins (consecutive)
                        if wn==wo(1)
                            n_bins_disregarded_beginning    =keys.n_consecutive_bins_significant+round((comparisons(comp).order_onset{2}-keys.PSTH_WINDOWS{wo(1),3})/keys.PSTH_binwidth);
                        else
                            n_bins_disregarded_beginning=keys.n_consecutive_bins_significant;
                        end
                        onset=find(~isnan(direction(n_bins_disregarded_beginning+1:end)),1);
                        if ~isempty(onset) && onset==1 && n_bins_disregarded_beginning>0 % if tuning was already there in the first bin, go backwards !!
                            onset=find(isnan([direction(n_bins_disregarded_beginning:-1:1) NaN]),1)*-1+2;
                        end
                        if isempty(onset) || onset_found % only 1 onset per unit, assuming earlier windows come earlier in time!!
                            sigbins(t).per_group(g).comparison(comp).window(wn).tuning_onset(u)=NaN;
                            sigbins(t).per_group(g).comparison(comp).window(wn).tuning_direction(u)=NaN;
                        else
                            onset= onset+n_bins_disregarded_beginning;
                            sigbins(t).per_group(g).comparison(comp).window(wn).tuning_direction(u)=direction(onset);
                            onset= onset-keys.n_consecutive_bins_significant+round(keys.PSTH_WINDOWS{wo(1),3}/keys.PSTH_binwidth);
                            sigbins(t).per_group(g).comparison(comp).window(wn).tuning_onset(u)=onset;
                            onset_found=1;
                        end
                    end
                end
            end
            
            %% N_tuned
            for wn=1:numel(current(1).window)
                concatinated_after_onset=[sigbins(t).per_group(g).comparison(comp).window(wn).bins];
                sigbins(t).per_group(g).comparison(comp).window(wn).n_tuned_cells=[sum(concatinated_after_onset==1,1); sum(concatinated_after_onset==-1,1)];
            end
            
            %% unit_order dependent on tuning onset
            tuning_onset_all_windows=vertcat(sigbins(t).per_group(g).comparison(comp).window(wo).tuning_onset);
            tuning_onset=tuning_onset_all_windows(1,:);
            for ww=numel(wo):-1:1
                to_replace=~isnan(tuning_onset_all_windows(ww,:));
                tuning_onset(to_replace)=tuning_onset_all_windows(ww,to_replace)+ww*1000; %% assuming no window has more than 1000 bins
            end
            [~, sigbins(t).per_group(g).comparison(comp).unit_order]=sort(tuning_onset);
        end
    end
end


%% plots
for t=1:numel(sigbins)
    typ=UC.type(mod(t-1,numel(UC.type))+1);
    for g=1:numel(unique_group_values)
        N_total=sum(ismember(group_values,unique_group_values{g}));
        current=[sigbins(t).per_group(g)];
        
     %% ideally, the names are taken from conditions as well, not hardcoded tt.hands and tt.choices !!
        
        condition_title=ph_get_condition_title(keys);
        fig_title=sprintf('Selection: %s %s %s %s %s, %s = %s N=%s ',...
            [Sel_for_title{:}],keys.monkey,[keys.conditions_to_plot{:}],condition_title,keys.arrangement,keys.ON.group_parameter,unique_group_values{g},mat2str(double(N_total)));
        filename=sprintf('%s %s %s = %s %s %s %s',...
            [Sel_for_title{:}],keys.monkey,keys.ON.group_parameter,unique_group_values{g},[keys.conditions_to_plot{:}],keys.arrangement(1:3),condition_title);
        
        %% save metadata
        unitidx=ismember(complete_unit_list,tuning_per_unit_table(ismember(group_values,unique_group_values(g)),idx_unitID));
        unit_IDs=complete_unit_list(all(unitidx,2));
        comparison=current.comparison;
        save([keys.basepath_to_save, keys.project_version, filesep, 'population_meta_data', filesep, 'response timing', filesep, keys.ON.comparisons_title ' ' filename], 'keys','tuning_per_unit_table','unit_IDs','comparison');
        plot_1_title            = [keys.ON.comparisons_title ' ' fig_title  ' per bin'];
        
        %% tuning per bin plot
        PSTH_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',plot_1_title);
        column=1;
        keys=ph_get_epoch_keys(keys,typ,UC.effector,sum(UC.type_effector(:,1)==typ)>1);%% does it make sense to distinguish by effector?
        n=1;
        for comp=1:numel(comparisons)
            wo=find(ismember(keys.PSTH_WINDOWS(:,1),comparisons(comp).order_onset{1}));
            subplot(numel(comparisons),1,(comp-1)+column)
            hold on
            col=comparisons(comp).colors;
            unit_order=current(n).comparison(comp).unit_order;
            state_shift=0;
            for w=1:size(keys.PSTH_WINDOWS,1)
                t_before_state=keys.PSTH_WINDOWS{w,3};
                t_after_state=keys.PSTH_WINDOWS{w,4};
                bins=t_before_state:keys.PSTH_binwidth:t_after_state;
                bins=bins+state_shift-t_before_state;
                to_plot1=current(n).comparison(comp).window(w).bins(unit_order,:);
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
            tr=[all_trialz.type]==typ & ismember([all_trialz.effector],UC.effector) & [all_trialz.completed]==1;
            ph_PSTH_background(all_trialz(tr),y_lim,y_lim,y_lim,keys,1)
            % set ydata
            uistack(onset_handle(comp), 'top');
        end
        ph_title_and_save(PSTH_summary_handle,  [keys.ON.comparisons_title ' ' filename ' PSTHs'],plot_1_title,keys)
        
        
        %% n tuned cells plot
        plot_3_title            = [keys.ON.comparisons_title ' ' fig_title  ' n tuned cells'];
        PSTH_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',plot_3_title);
        column=1;
        keys=ph_get_epoch_keys(keys,typ,UC.effector,sum(UC.type_effector(:,1)==typ)>1);%% does it make sense to distinguish by effector?
        n=1;
        clear x_lim y_lim sp
        for comp=1:numel(comparisons)
            subplot(numel(comparisons),1,(comp-1)+column)
            hold on
            col=comparisons(comp).colors;
            unit_order=current(n).comparison(comp).unit_order;
            state_shift=0;
            for w=1:size(keys.PSTH_WINDOWS,1)
                t_before_state=keys.PSTH_WINDOWS{w,3};
                t_after_state=keys.PSTH_WINDOWS{w,4};
                bins=t_before_state:keys.PSTH_binwidth:t_after_state;
                bins=bins+state_shift-t_before_state;
                to_plot1=sigbins(t).per_group(g).comparison(comp).window(w).n_tuned_cells(1,:)/N_total*100;
                to_plot2=sigbins(t).per_group(g).comparison(comp).window(w).n_tuned_cells(2,:)/N_total*100;
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
            tr=[all_trialz.type]==typ & ismember([all_trialz.effector],UC.effector) & [all_trialz.completed]==1;
            ph_PSTH_background(all_trialz(tr),lims,lims,lims,keys,1)
        end
        
        ph_title_and_save(PSTH_summary_handle,  [keys.ON.comparisons_title ' ' filename ' fraction tuned'],plot_3_title,keys)
        
        
        %% tuning onset plot
        plot_2_title            = [keys.ON.comparisons_title ' ' fig_title  ' tuning onset'];
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
            unit_order=current(n).comparison(comp).unit_order;
            N_total=numel(unit_order); % for getting fraction
            state_shift=0;
            for w=1:size(keys.PSTH_WINDOWS,1)
                tuning_direction=current(n).comparison(comp).window(w).tuning_direction;
                t_before_state=keys.PSTH_WINDOWS{w,3};
                t_after_state=keys.PSTH_WINDOWS{w,4};
                bins=t_before_state:keys.PSTH_binwidth:t_after_state;
                bins=bins+state_shift-t_before_state;
                onset1=sigbins(t).per_group(g).comparison(comp).window(w).tuning_onset(tuning_direction==1);
                onset2=sigbins(t).per_group(g).comparison(comp).window(w).tuning_onset(tuning_direction==-1);
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
            tr=[all_trialz.type]==typ & ismember([all_trialz.effector],UC.effector) & [all_trialz.completed]==1;
            % set ydata
            ph_PSTH_background(all_trialz(tr),lims,lims,lims,keys,1)
            uistack(onset_handle(comp), 'top');
        end
        
        ph_title_and_save(PSTH_summary_handle,  [keys.ON.comparisons_title ' ' filename ' tuning onset'],plot_2_title,keys)
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
    if ends(x)+1 > starts(x)+keys.n_consecutive_bins_significant
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
    confidence_criterion=0.95;
    Dav=abs(nanmean(Amat)-nanmean(Bmat));
    P_all=[Amat;Bmat];
    
    na=size(Amat,1);
    nb=size(Bmat,1);    
    N=na+nb;
    n_permutations=1000;
    n_for_sig=n_permutations*confidence_criterion;
    
    if paired
        for p=1:n_permutations
            to_switch=randi([0 1], na,1)==1;
            Ap=Amat;
            Bp=Bmat;
            Ap(to_switch,:)=Bmat(to_switch,:);
            Bp(to_switch,:)=Amat(to_switch,:);
            Pav=abs(nanmean(Ap)-nanmean(Bp));
            Pwrong(p,:)=Dav>Pav;
        end
        h=sum(Pwrong)>=n_for_sig;
    else
        for p=1:n_permutations
            p_idx=randsample(1:N,N);
            Pav=abs(nanmean(P_all(p_idx(1:na),:))-nanmean(P_all(p_idx(na+1:N),:)));
            Pwrong(p,:)=Dav>Pav;
        end
        h=sum(Pwrong)>=n_for_sig;
    end
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