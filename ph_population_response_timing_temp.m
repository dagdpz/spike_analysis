function ph_population_response_timing_temp(population,modified_keys)
%load('W:\Projects\Pulv_microstim_behavior\ephys\ephys_analysis_v5_July2016_coordinates\monkeys_Combined_mem.mat')
keys.n_consecutive_bins_significant=5;
warning('off','MATLAB:catenate:DimensionMismatch');

%%% !!!!!!!! make sure epochs (baseline, cue, .... make sense for all effectors)

Left_hemisphere_targets={'dPulv_l','pdSTS_L','FST_L','pTPO_L','MIP','MIP_L','unknown'};
Right_hemisphere_targets={'dPulv_r','MIP_R'};

for fn=fieldnames(modified_keys)'
    keys.(fn{:})=modified_keys.(fn{:});
end

legend_labels={'NH IS IN' 'NH IS CH' 'IH IS IN' 'IH IS CH' 'CH IS IN' 'CH IS CH' 'NH CS IN' 'NH CS CH' 'IH CS IN' 'IH CS CH' 'CH CS IN' 'CH CS CH' };
cols=keys.colors;
keys.line_colors=[cols.NH_IS_IN;cols.NH_IS_CH;cols.IH_IS_IN;cols.IH_IS_CH;cols.CH_IS_IN;cols.CH_IS_CH;
    cols.NH_CS_IN;cols.NH_CS_CH;cols.IH_CS_IN;cols.IH_CS_CH;cols.CH_CS_IN;cols.CH_CS_CH;]/255;
n_col=size(keys.line_colors,1);

%% tuning table preparation and grouping
[tuning_per_unit_table]                 = ph_load_extended_tuning_table(keys);
[tuning_per_unit_table, Sel_for_title]  = ph_reduce_tuning_table(tuning_per_unit_table,keys);
% if keys.FR_subtract_baseline
%     Sel_for_title =[Sel_for_title,{'base';'=';keys.epoch_BL;', '}];
% end
idx_unitID=DAG_find_column_index(tuning_per_unit_table,'unit_ID');
idx_group_parameter=DAG_find_column_index(tuning_per_unit_table,keys.ON.group_parameter);
group_values=tuning_per_unit_table(:,idx_group_parameter);
group_values=cellfun(@num2str, group_values, 'UniformOutput', false);
cell_in_any_group=[false; ~ismember(group_values(2:end),keys.ON.group_excluded)];
unique_group_values=unique(group_values(cell_in_any_group));
if isempty(unique_group_values)
    disp('no relevant groups found');
    return;
end
complete_unit_list={population.unit_ID}';

%% define conditions to look at
all_trialz=[population.trial];

%condition_parameters  ={'type','effector','reach_hand','choice'};
%condition_parameters  ={'reach_hand','choice','hemifield'};
condition_parameters  ={'reach_hand','choice','perturbation'};
per_trial.types       =[all_trialz.type];
per_trial.effectors   =[all_trialz.effector];
per_trial.hands       =[all_trialz.reach_hand];
per_trial.choice      =[all_trialz.choice];
per_trial.perturbation=[all_trialz.perturbation];


%hemifields=unique([whatisthis.trial.hemifield]);
%u_hemifields=[-1,0,1]; % why does this have to be hardcoded? ---> Because case not defined yet, case defines positions !!
u_hemifields=[-1,1]; % why does this have to be hardcoded? And what do we do with vertical targets?
u_types     =unique(per_trial.types);
u_effectors =unique(per_trial.effectors);
u_hands     =unique(per_trial.hands);
u_choice    =unique(per_trial.choice);
u_perturbation    =unique(per_trial.perturbation);
u_perturbation=u_perturbation(~isnan(u_perturbation));

%% limit conditions key?
u_hands     =u_hands(ismember(u_hands,keys.tt.hands));
u_choice    =u_choice(ismember(u_choice,keys.tt.choices));

all_type_effectors      = combvec(u_types,u_effectors)';
type_effectors =[];

% redifine type_effectors to include only relevant
for t=1:size(all_type_effectors,1)
    typ=all_type_effectors(t,1);
    eff=all_type_effectors(t,2);
    [~, type_effector_short{t}]=get_type_effector_name(typ,eff);
    if ~ismember(type_effector_short{t},keys.conditions_to_plot) %|| sum(tr_con)<1
        continue;
    end
    type_effectors=[type_effectors; all_type_effectors(t,:)];
end
type_effector_short(~ismember(type_effector_short,keys.conditions_to_plot))=[];
u_types     =unique(type_effectors(:,1))';
u_effectors =unique(type_effectors(:,2))';

condition_matrix    = combvec(u_hands,u_choice, u_perturbation,u_hemifields)';
conditions_out      = combvec(u_effectors,u_hands,u_choice, u_perturbation)';
conditions_hf       = combvec(u_hemifields,conditions_out')';


%% Convert to ipsi/contra, Baseline subtraction, normalization, re-ordering, and gaussian RFs
tr_con=ismember([all_trialz.completed],keys.cal.completed);
[whatisthis]=ph_arrange_positions_and_plots(all_trialz(tr_con),keys);
positions=unique(vertcat(whatisthis.trial.position),'rows');
clear whatisthis

%% type? (+effector?)
for tye=1:size(type_effectors,1) %typ=unique(per_trial.types)
    typ=type_effectors(tye,1);
    eff=type_effectors(tye,2);
    conditions_eff=conditions_out(conditions_out(:,1)==eff,2:end);
    t=find(u_types==typ);
    
    %% here we need to work on...!!! 'in_NH_' NOT IDEAL AT ALL
    % why do we need this at all????
    idx_existing=DAG_find_column_index(tuning_per_unit_table,['existing_' 'in_AH_' type_effector_short{tye} '_' keys.arrangement(1:3)]);
    tya_existing{t}= [true; ismember(vertcat(tuning_per_unit_table{2:end,idx_existing}),true)];
    keys=ph_get_epoch_keys(keys,typ,eff,sum(type_effectors(:,1)==typ)>1);
    
    for g=1:numel(unique_group_values)
        unitidx=ismember(complete_unit_list,tuning_per_unit_table(ismember(group_values,unique_group_values(g))&tya_existing{t},idx_unitID));
        units=find(all(unitidx,2))';
        for u=units
            tr_con=ismember([population(u).trial.completed],keys.cal.completed);
            [pop]=ph_arrange_positions_and_plots(population(u).trial(tr_con),keys);
            pop=ph_LR_to_CI(pop,population(u).target);
            
            %% average PSTH per unit
            for ce=1:size(conditions_eff,1)
                c=find(ismember(conditions_out,[eff conditions_eff(ce,:)],'rows'));
                clear trpar
                
                for par=1:numel(condition_parameters)
                    fn=condition_parameters{par};
                    trpar(par,:)=[pop.trial.(fn)]==conditions_eff(ce,par); %condition_matrix_wohf?
                end
                trpar(end+1,:)=[pop.trial.type]==typ & [pop.trial.effector]==eff;
                tr_con=all(trpar,1);
                per_epoch=vertcat(pop.trial(tr_con).epoch);
                
                % why hemifields independently in this piece of code? cause
                % meybe we want to add positions...
                for f=1:numel(u_hemifields) %hemifield (as compared to positions!
                    trtemp=[pop.trial.hemifield]==u_hemifields(f);
                    trials_for_SD=find(tr_con & trtemp);
                    
                    % epoch averages for specific conditions
                    if any(trtemp(tr_con))
                        condition(t,c).per_hemifield(f).unit(u).epoch_averages=...
                            nanmean(reshape([per_epoch(trtemp(tr_con),:).FR],size(per_epoch(trtemp(tr_con),:))),1);
                    else
                        condition(t,c).per_hemifield(f).unit(u).epoch_averages=NaN(size(vertcat(pop.trial(1).epoch)));
                    end
                    
                    for wn=1:size(keys.PSTH_WINDOWS,1)
                        % calculate significance for a few bins surrounding
                        temp_window=keys.PSTH_WINDOWS(wn,:);
                        keys.PSTH_WINDOWS{wn,3}=keys.PSTH_WINDOWS{wn,3}-keys.n_consecutive_bins_significant*keys.PSTH_binwidth;
                        keys.PSTH_WINDOWS{wn,4}=keys.PSTH_WINDOWS{wn,4}+keys.n_consecutive_bins_significant*keys.PSTH_binwidth;
                        for tr=1:numel(trials_for_SD)
                            condition(t,c).per_hemifield(f).window(wn).unit(u).average_spike_density(tr,:)=...
                                ph_spike_density(pop.trial(trials_for_SD(tr)),wn,keys,0,1);
                        end
                        if isempty(trials_for_SD)
                            condition(t,c).per_hemifield(f).window(wn).unit(u).average_spike_density(1,:)=...
                                NaN(size(ph_spike_density(pop.trial(1),wn,keys,0,1)));
                        end
                        keys.PSTH_WINDOWS(wn,:)=temp_window;
                    end
                end
            end
        end
    end
end
%end

%% condition comparison
comparisons_per_effector=keys.ON.comparisons_per_effector;
condition_parameters_comparissons = [{'hemifield'} {'effector'} condition_parameters];

% effector comparison only possible if several_effectors
for t=1:size(condition,1)
    typ=u_types(t); %u_types(mod(t-1,numel(u_types))+1);
    current=[condition(t,:).per_hemifield];
    current_unit=vertcat(current.unit);
    current_window=vertcat(current.window);
    unitidx=ismember(complete_unit_list,tuning_per_unit_table(ismember(group_values,unique_group_values(g)) & tya_existing{t},idx_unitID));
    units=find(all(unitidx,2));
    for comp=1:numel(comparisons_per_effector)
        cM1=true(size(conditions_hf));
        cM2=true(size(conditions_hf));
        for k=1:numel(condition_parameters_comparissons)
            if isfield(comparisons_per_effector,condition_parameters_comparissons{k}) &&...
                    ~ isempty(comparisons_per_effector(comp).(condition_parameters_comparissons{k}))
                cM1(~ismember(conditions_hf(:,k),comparisons_per_effector(comp).(condition_parameters_comparissons{k}){1}),k)=false;
                cM2(~ismember(conditions_hf(:,k),comparisons_per_effector(comp).(condition_parameters_comparissons{k}){2}),k)=false;
            end
        end
        n=1;
        c1=find(all(cM1,2));
        c2=find(all(cM2,2));
        for u=1:numel(units)
            %% baseline definition (for epoch tuning so far... CHECK dimensions of epoch_averages!!)
            epoch_averages=vertcat(current_unit(unique([c1; c2]),units(u)).epoch_averages);
            baseline=nanmean(epoch_averages(:,ismember(keys.EPOCHS(:,1),comparisons_per_effector(comp).baseline_epoch)));
            for wn=1:numel(current(1).window)
                temp_sigbins=ph_compare_by_bin_by_trial(current_window(:,wn),c1,c2,baseline,units(u),keys.n_consecutive_bins_significant);
                sigbins(t).per_effector(n).comparison(comp).window(wn).bins(u,:)=temp_sigbins(keys.n_consecutive_bins_significant+1:end-keys.n_consecutive_bins_significant);
            end
        end
        
        %% re-order dependent on tuning onset
        keys=ph_get_epoch_keys(keys,typ,u_effectors,sum(type_effectors(:,1)==typ)>1);%% does it make sense to distinguish by effector?
        wo=find(ismember(keys.PSTH_WINDOWS(:,1),comparisons_per_effector(comp).order_onset{1})):size(keys.PSTH_WINDOWS,1);
        n_bins_disregarded_beginning    =round((comparisons_per_effector(comp).order_onset{2}-keys.PSTH_WINDOWS{wo(1),3})/keys.PSTH_binwidth);
        n_bins_disregarded_end          =round((comparisons_per_effector(comp).order_onset{3}-keys.PSTH_WINDOWS{wo(1),3})/keys.PSTH_binwidth);
        n_bins_regarded                 =round((comparisons_per_effector(comp).order_onset{3}-comparisons_per_effector(comp).order_onset{2})/keys.PSTH_binwidth+1);
        concatinated_after_onset=[sigbins(t).per_effector(n).comparison(comp).window(wo).bins];
        concatinated_after_onset(:,end+1)=0;
        clear tuning_onset
        % i feel there should be a smarter way then just looping again
        for u=1:numel(units)
            tuning_onset(u)=find(~isnan([concatinated_after_onset(u,n_bins_disregarded_beginning+1:end)]),1); % in bins?
            if tuning_onset(u)==1 && n_bins_disregarded_beginning>0%% if tuning was already there in the first bin, go backwards !!
                tuning_onset(u)=find(isnan([concatinated_after_onset(u,n_bins_disregarded_beginning:-1:1) NaN]),1)*-1+2; % in bins?
            end
            tuning_direction(u)=concatinated_after_onset(u,tuning_onset(u)+n_bins_disregarded_beginning);
        end
        concatinated_after_onset(:,1:n_bins_disregarded_beginning)=[];
        tuning_direction(tuning_onset>n_bins_disregarded_end-n_bins_disregarded_beginning)=NaN;
        [~, sigbins(t).per_effector(n).comparison(comp).unit_order]=sort(tuning_onset);
        sigbins(t).per_effector(n).comparison(comp).tuning_onset=tuning_onset;%+comparisons_per_effector(comp).order_onset{2};
        sigbins(t).per_effector(n).comparison(comp).tuning_direction=tuning_direction;
        sigbins(t).per_effector(n).comparison(comp).n_tuned_cells=[sum(concatinated_after_onset(:,1:n_bins_regarded)==1,1); sum(concatinated_after_onset(:,1:n_bins_regarded)==-1,1)];
    end
end

%% plots
for t=1:numel(sigbins)
    typ=u_types(mod(t-1,numel(u_types))+1);
    current=[sigbins(t).per_effector];
    
    fig_title=sprintf('Selection: %s %s %s hnd %s ch %s %s grouped by %s',...
        [Sel_for_title{:}],keys.monkey,[keys.conditions_to_plot{:}],mat2str(u_hands),mat2str(double(u_choice)),keys.arrangement,keys.ON.group_parameter);
    filename=sprintf('%s %s %s %s %s hnd %s ch %s',...
        [Sel_for_title{:}],keys.monkey,keys.ON.group_parameter,[keys.conditions_to_plot{:}],keys.arrangement(1:3),mat2str(u_hands),mat2str(double(u_choice)));
    
    plot_1_title            = [keys.ON.comparisons_title ' ' fig_title  ' per bin'];
    
    %% tuning per bin plot
    PSTH_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_1_title);
    % for g=1:numel(unique_group_values)
    %         unitidx=ismember(complete_unit_list,tuning_per_unit_table(ismember(group_values,unique_group_values(g)) & tya_existing{t},idx_unitID));
    %         group_units=find(all(unitidx,2))';
    %
    %         if true %%reducing to only units that have each condition
    %             current_window=vertcat(current.window);
    %             current_units=vertcat(current_window(:,1).unit);
    %             group_units=intersect(group_units,find(all(arrayfun(@(x) ~isempty(x.average_spike_density),current_units),1)));
    %         end
    column=1;
    keys=ph_get_epoch_keys(keys,typ,u_effectors,sum(type_effectors(:,1)==typ)>1);%% does it make sense to distinguish by effector?
    n=1;
    for comp=1:numel(comparisons_per_effector)
        subplot(numel(comparisons_per_effector),1,(comp-1)+column)
        hold on
        col=comparisons_per_effector(comp).colors;
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
            state_shift=state_shift+t_after_state-t_before_state+0.1;
        end
        %                 legend_line_handles=[legend_line_handles errorbarhandle.mainLine];
        title(comparisons_per_effector(comp).title);
    end
    
    %title(sprintf('%s = %s, N (NH/CH/IH) =%s ',keys.ON.group_parameter,unique_group_values{g},num2str(n_units)),'interpreter','none');
    y_lim(n,comp)=size(to_plot1,1);
    
    %% subplot appearance, and tuning lines
    ylimmax=max(max(y_lim));
    ylimmin=0;
    y_lim=[ylimmin ylimmax];
    %     if ~isempty(keys.population_ylim)
    %         y_lim= keys.population_ylim;
    %     end
    column=1;
    for comp=1:numel(comparisons_per_effector)
        hold on
        subplot(numel(comparisons_per_effector),1,(comp-1)+column)
        
        %% choices? hands? errm here its only relevant for event and epoch onsets
        tr=[all_trialz.type]==typ & ismember([all_trialz.effector],u_effectors) & [all_trialz.completed]==1;
        %set(gca,'xtick',get(gca,'xtick')*keys.PSTH_binwidth);
        ph_PSTH_background(all_trialz(tr),y_lim,y_lim,y_lim,keys,1)
        
    end
    ph_title_and_save(PSTH_summary_handle,  [keys.ON.comparisons_title ' ' filename ' PSTHs'],plot_1_title,keys)
    
    
    %% tuning onset plot
    plot_2_title            = [keys.ON.comparisons_title ' ' fig_title  ' tuning onset'];
    PSTH_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_2_title);
    % for g=1:numel(unique_group_values)
    %         unitidx=ismember(complete_unit_list,tuning_per_unit_table(ismember(group_values,unique_group_values(g)) & tya_existing{t},idx_unitID));
    %         group_units=find(all(unitidx,2))';
    %
    %         if true %%reducing to only units that have each condition
    %             current_window=vertcat(current.window);
    %             current_units=vertcat(current_window(:,1).unit);
    %             group_units=intersect(group_units,find(all(arrayfun(@(x) ~isempty(x.average_spike_density),current_units),1)));
    %         end
    column=1;
    keys=ph_get_epoch_keys(keys,typ,u_effectors,sum(type_effectors(:,1)==typ)>1);%% does it make sense to distinguish by effector?
    
    n=1;% find(u_effectors==eff);
    clear x_lim sp
    for comp=1:numel(comparisons_per_effector)
        sp(comp)=subplot(numel(comparisons_per_effector),1,(comp-1)+column);
        hold on
        col=comparisons_per_effector(comp).colors;
        
        unit_order=current(n).comparison(comp).unit_order;
        N_total=numel(unit_order);
        tuning_direction=current(n).comparison(comp).tuning_direction;
        to_plot1=(sigbins(t).per_effector(n).comparison(comp).tuning_onset(tuning_direction==1)-1)*keys.PSTH_binwidth*1000;
        to_plot2=(sigbins(t).per_effector(n).comparison(comp).tuning_onset(tuning_direction==-1)-1)*keys.PSTH_binwidth*1000;
        to_plot1(to_plot1<0)=[];to_plot2(to_plot2<0)=[]; % remove tuning onsets before respective threshold
        
        maxbin=max([to_plot1 to_plot2])/1000+keys.ON.comparisons_per_effector(comp).order_onset{2};
        bins=(keys.ON.comparisons_per_effector(comp).order_onset{2}:keys.PSTH_binwidth:maxbin)*1000;%keys.ON.comparisons_per_effector(comp).order_onset{3})*1000;
        mean1=nanmean(to_plot1); mean2=nanmean(to_plot2);
        median1=nanmedian(to_plot1); median2=nanmedian(to_plot2);
        std1=std(to_plot1); std2=std(to_plot2);
        to_plot1=hist(to_plot1,(bins-bins(1)))/N_total*100;
        to_plot2=hist(to_plot2,(bins-bins(1)))/N_total*100;
        %to_plot1(1)=0; to_plot1(2)=0; % removing
        
        plot(bins,to_plot1,'color',col(1,:),'linewidth',1);
        text(bins(2),max([to_plot1,to_plot2])*2/3,num2str(round(sum(to_plot1))),'color',col(1,:));
        text(bins(2),max([to_plot1,to_plot2])/3,[num2str(round(mean1)) ' ' num2str(round(median1)) ' std=' num2str(round(std1))] ,'color',col(1,:));
        plot(bins,to_plot2,'color',col(2,:),'linewidth',1);
        text(bins(4),max([to_plot1,to_plot2])*2/3,num2str(round(sum(to_plot2))),'color',col(2,:));
        text(bins(4),max([to_plot1,to_plot2])/3,[num2str(round(mean2)) ' ' num2str(round(median2)) ' std=' num2str(round(std2))] ,'color',col(2,:));
        x_lim(comp,:)=get(gca,'xlim');
        title(comparisons_per_effector(comp).title);
    end
    
    % setting same x axis for all subplots
    max_x_diff=max(diff(x_lim,1,2));
    for comp=1:numel(comparisons_per_effector)
        subplot(sp(comp));
        set(gca,'xlim',[x_lim(comp,1) x_lim(comp,1)+max_x_diff]);
    end
    
    ph_title_and_save(PSTH_summary_handle,  [keys.ON.comparisons_title ' ' filename ' tuning onset'],plot_2_title,keys)
    %             legend(legend_line_handles,legend_labels(legend_label_indexes));
    
    %% n tuned cells plot
    plot_3_title            = [keys.ON.comparisons_title ' ' fig_title  ' n tuned cells'];
    PSTH_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_3_title);
    % for g=1:numel(unique_group_values)
    %         unitidx=ismember(complete_unit_list,tuning_per_unit_table(ismember(group_values,unique_group_values(g)) & tya_existing{t},idx_unitID));
    %         group_units=find(all(unitidx,2))';
    %
    %         if true %%reducing to only units that have each condition
    %             current_window=vertcat(current.window);
    %             current_units=vertcat(current_window(:,1).unit);
    %             group_units=intersect(group_units,find(all(arrayfun(@(x) ~isempty(x.average_spike_density),current_units),1)));
    %         end
    column=1;
    keys=ph_get_epoch_keys(keys,typ,u_effectors,sum(type_effectors(:,1)==typ)>1);%% does it make sense to distinguish by effector?
    n=1;
    clear x_lim sp
    for comp=1:numel(comparisons_per_effector)
        
        sp(comp)=subplot(numel(comparisons_per_effector),numel(u_effectors),(comp-1)*numel(u_effectors)+column);
        hold on
        col=comparisons_per_effector(comp).colors;
        unit_order=current(n).comparison(comp).unit_order;
        N_total=numel(unit_order); % for getting fraction
        
        to_plot1=sigbins(t).per_effector(n).comparison(comp).n_tuned_cells(1,:)/N_total*100;
        to_plot2=sigbins(t).per_effector(n).comparison(comp).n_tuned_cells(2,:)/N_total*100;
        bins=(keys.ON.comparisons_per_effector(comp).order_onset{2}:keys.PSTH_binwidth:keys.ON.comparisons_per_effector(comp).order_onset{3})*1000;%keys.ON.comparisons_per_effector(comp).order_onset{3})*1000;
        
        max1=bins(to_plot1==max(to_plot1));
        max2=bins(to_plot2==max(to_plot2));
        plot(bins,to_plot1,'color',col(1,:),'linewidth',1);
        plot(bins,to_plot2,'color',col(2,:),'linewidth',1);
        y_lim=get(gca,'ylim');
        if ~isempty(max1)
            max1=nanmean(max1);
            line([max1 max1],[0 y_lim(2)],'color',col(1,:),'linewidth',1);
            text(max1-(bins(end)-bins(1))/100,y_lim(2)/2,num2str(max1),'rotation',90);
        end
        if ~isempty(max2)
            max2=nanmean(max2);
            line([max2 max2],[0 y_lim(2)],'color',col(2,:),'linewidth',1);
            text(max2-(bins(end)-bins(1))/100,y_lim(2)/2,num2str(max2),'rotation',90);
        end
        
        %                 legend_line_handles=[legend_line_handles errorbarhandle.mainLine];
        x_lim(comp,:)=get(gca,'xlim');
        title(comparisons_per_effector(comp).title);
    end
    
    % setting same x axis for all subplots
    max_x_diff=max(diff(x_lim,1,2));
    for comp=1:numel(comparisons_per_effector)
        subplot(sp(comp));
        set(gca,'xlim',[x_lim(comp,1) x_lim(comp,1)+max_x_diff]);
    end
    
    ph_title_and_save(PSTH_summary_handle,  [keys.ON.comparisons_title ' ' filename ' fraction tuned'],plot_3_title,keys)
    %             legend(legend_line_handles,legend_labels(legend_label_indexes));
    
end
end

function hn=ph_compare_by_bin_by_trial(in,c1,c2,baseline,u,n_consecutive_bins)
C1=[];
C2=[];
for cc=c1(:)'
    C1=[C1;vertcat(in(cc).unit(u).average_spike_density)];
end
if all(c1==c2) %% epoch comparison (only possibly desired comparison if conditions are exactly identical)
    ho=ttest(C1,baseline,0.05,'both',1);
    C2=ones(1,size(C1,2))*baseline;
else
    for cc=c2(:)'
        C2=[C2;vertcat(in(cc).unit(u).average_spike_density)];
    end
    ho=ttest2(C1,C2,0.05,'both','equal',1);
end

% not sure why there are nan values though (STS memory saccades)
ho(isnan(ho))=0;

%% keep only if n_consecutive_bins are significant
starts=find(diff([0 ho])==1);
ends=find(diff([ho 0])==-1);
hn=NaN(size(ho));
for x=1:numel(starts)
    if ends(x)+1 > starts(x)+n_consecutive_bins
        hn(starts(x):ends(x))=1;
    end
end

%% negative for condition 2 larger
hn(hn==1 & nanmean(C2)>nanmean(C1))=-1;
end