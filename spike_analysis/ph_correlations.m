function ph_correlations(population,modified_keys)
%load('W:\Projects\Pulv_microstim_behavior\ephys\ephys_analysis_v5_July2016_coordinates\monkeys_Combined_mem.mat')
global MA_STATES
warning('off','MATLAB:catenate:DimensionMismatch');

%%% !!!!!!!! make sure epochs (baseline, cue, .... make sense for all effectors)

keys.saccade_states={'PreS','PeriS','PostS'};
keys.reach_states={'PreR','PeriR','PostR'};

Left_hemisphere_targets={'dPulv_l','pdSTS_L','FST_L','pTPO_L','MIP','MIP_L','unknown'};
Right_hemisphere_targets={'dPulv_r','MIP_R'};

for fn=fieldnames(modified_keys)'
    keys.(fn{:})=modified_keys.(fn{:});
end
keys.path_to_save=[keys.drive filesep keys.basepath_to_save filesep keys.project_version filesep 'correlation_analysis' filesep];
if ~exist(keys.path_to_save,'dir')
    mkdir([keys.drive filesep keys.basepath_to_save filesep keys.project_version ], 'correlation_analysis');
end


%% tuning table preparation and grouping
keys.tt.tasktypes={};

[tuning_per_unit_table]                 = ph_load_extended_tuning_table(keys);
[tuning_per_unit_table, Sel_for_title]  = ph_reduce_tuning_table(tuning_per_unit_table,keys);
if keys.FR_subtract_baseline
    Sel_for_title =[Sel_for_title,{'base';'=';keys.epoch_BL;', '}];
end
idx_unitID=find_column_index(tuning_per_unit_table,'unit_ID');

% idx_group_parameter=find_column_index(tuning_per_unit_table,keys.population_group_parameter);
% idx_unitID=find_column_index(tuning_per_unit_table,'unit_ID');
% group_values=tuning_per_unit_table(:,idx_group_parameter);
% group_values=cellfun(@num2str, group_values, 'UniformOutput', false);
% cell_in_any_group=[false; ~ismember(group_values(2:end),keys.population_group_excluded)];
% unique_group_values=unique(group_values(cell_in_any_group));
% if isempty(unique_group_values)
%     disp('no relevant groups found');
%     return;
% end
% complete_unit_list={population.unit_ID}';

%% define conditions to look at
all_trialz=[population.trial];

%condition_parameters  ={'type','effector','reach_hand','choice'};
%condition_parameters  ={'reach_hand','choice','hemifield'};
%condition_parameters  ={'reach_hand','choice'};
per_trial.types       =[all_trialz.type];
per_trial.effectors   =[all_trialz.effector];
per_trial.hands       =[all_trialz.reach_hand];
per_trial.choice      =[all_trialz.choice];


%hemifields=unique([whatisthis.trial.hemifield]);
%u_hemifields=[-1,1]; % why does this have to be hardcoded? And what do we do with vertical targets?
u_types     =unique(per_trial.types);
u_effectors =unique(per_trial.effectors);
u_hands     =unique(per_trial.hands);
u_choice    =unique(per_trial.choice);
keys.conditions_to_plot={'Ddre'};
hand_labels_LR={'NH','LH','RH'};
condition_parameters_pairs  ={'reach_hand','choice','pos'};

%% limit conditions key?
u_hands     =u_hands(ismember(u_hands,keys.limit_conditions.hands));
u_choice    =u_choice(ismember(u_choice,keys.limit_conditions.choice));

all_type_effectors      = combvec(u_types,u_effectors)';
type_effectors =[];


% redifine type_effectors to include only relevant
for t=1:size(all_type_effectors,1)
    typ=all_type_effectors(t,1);
    eff=all_type_effectors(t,2);
    [~, type_effector_short]=get_type_effector_name(typ,eff);
    if ~ismember(type_effector_short,keys.conditions_to_plot) %|| sum(tr_con)<1
        continue;
    end
    type_effectors=[type_effectors; all_type_effectors(t,:)];
end
u_types     =unique(type_effectors(:,1));
u_effectors =unique(type_effectors(:,2));
%[population.session]=deal(1);
%% Convert to ipsi/contra, Baseline subtraction, normalization, re-ordering, and gaussian RFs
for a=1:numel(keys.position_and_plotting_arrangements)
    keys.arrangement=keys.position_and_plotting_arrangements{a}; % important for ph_arrange_positions_and_plots
    
    tr_con=[all_trialz.success]==1;
    [whatisthis]=ph_arrange_positions_and_plots(all_trialz(tr_con),keys);
    positions=unique(vertcat(whatisthis.trial.position),'rows');
    
    clear whatisthis
    
    %% type? (+effector?)
    for t=1:size(type_effectors,1) %careful, because type_effectors is limited by keys.conditions_to_plot
        typ=type_effectors(t,1);
        eff=type_effectors(t,2);
        get_expected_MA_states(typ,eff,sum(type_effectors(:,1)==typ)>1); %% does it make sense to distinguish by effector?
        keys.ALL_EPOCHS=keys.EPOCHS_PER_TYPE{typ};
        keys.EPOCHS=keys.ALL_EPOCHS(ismember([keys.ALL_EPOCHS{:,2}],MA_STATES.all_states),:);
        for hnd=1:numel(u_hands)
            
        condition_matrix_pairs    = combvec(u_hands(hnd),u_choice,[positions(:,1)+1i*positions(:,2)]')'; %% to take only instructed !
        [~, type_effector_short]=get_type_effector_name(typ,eff);
        tyc=hnd+(t-1)*numel(u_hands)+(a-1)*(numel(u_types)*numel(u_hands));
        tyc_label{tyc}=[type_effector_short hand_labels_LR{u_hands(hnd)+1} keys.arrangement];
      
        for u1=1:numel(population)
            tr_u1=[population(u1).trial.success]==1 & [population(u1).trial.type]==typ & [population(u1).trial.effector]==eff;
            
                pop_corr(u1,tyc).session=population(u1).session;
                pop_corr(u1,tyc).target=population(u1).target;
                pop_corr(u1,tyc).unit_ID=population(u1).unit_ID;
            if sum(tr_u1)==0
               continue; 
            end
            [pop1]=ph_arrange_positions_and_plots(population(u1).trial(tr_u1),keys);
                trpos=vertcat(pop1.trial.position);
                AAA=num2cell([trpos(:,1)+1i*trpos(:,2)].');
                [pop1.trial.pos]=AAA{:};
            n_trials_u1=0;
            for c=1:size(condition_matrix_pairs,1)
                clear trpar1
                for par=1:numel(condition_parameters_pairs)
                    fn=condition_parameters_pairs{par};
                    trpar1(par,:)=[pop1.trial.(fn)]==condition_matrix_pairs(c,par);
                end
                per_epoch1=vertcat(pop1.trial(all(trpar1,1)).epoch);
                n_trials_c=size(per_epoch1,1);
                if n_trials_c>0
                FR_per_epoch1=reshape([per_epoch1.FR],n_trials_c,size(per_epoch1,2));
                pop_corr(u1,tyc).FR(n_trials_u1+1:n_trials_u1+n_trials_c,:)=...
                    bsxfun(@minus,FR_per_epoch1,nanmean(FR_per_epoch1,1));
                
                pop_corr(u1,tyc).block(n_trials_u1+1:n_trials_u1+n_trials_c)=[pop1.trial(all(trpar1,1)).block]';
                
                %% adding tuning information
                for e=1:size(keys.EPOCHS,1)
                    if ismember(pop_corr(u1,tyc).target,Right_hemisphere_targets)
                       hand_labels_CI={'NH','CH','IH'}; 
                    else
                        hand_labels_CI={'NH','IH','CH'};                        
                    end
                    %tuning_title=['in_' hand_labels_CI{hnd} '_' keys.EPOCHS{e,1} '_position_' type_effector_short '_' keys.arrangement(1:3)];
                    tuning_title=['in_' hand_labels_CI{hnd} '_PreR_position_' type_effector_short '_' keys.arrangement(1:3)];
                    
                    row_idx=find(ismember(tuning_per_unit_table(:,idx_unitID),population(u1).unit_ID));
                    column_idx=find_column_index(tuning_per_unit_table,tuning_title);
                    if ~isempty(row_idx) && ~isempty(column_idx)
                        pop_corr(u1,tyc).tuning{e}=tuning_per_unit_table{row_idx,column_idx};
                    else
                        pop_corr(u1,tyc).tuning{e}='-';
                    end
                end
                
                n_trials_u1=n_trials_u1+n_trials_c;
                end
            end
        end
        
        
        %% loop through conditions - only from here on it is different to population analysis
        n_pairs.ipsi=0;
        n_pairs.cont=0;
        for u1=1:size(pop_corr,1)
            
            
            same_area_units=find(ismember({pop_corr(:,tyc).target},pop_corr(u1,tyc).target));
            same_session_units=find([pop_corr(:,tyc).session]==pop_corr(u1,tyc).session);
            same_session_units=same_session_units(same_session_units>u1); %% so we dont take same pair twice
            u1_blocks=unique([pop_corr(u1,tyc).block]);
            for u2=same_session_units
                
                u12_blocks=intersect(u1_blocks,unique([pop_corr(u2,tyc).block]));
                if ~isempty(u12_blocks)
                    if any(same_area_units==u2)
                        pairstruct='ipsi';
                    elseif true
                        pairstruct='cont';
                    end
                    n_pairs.(pairstruct)=n_pairs.(pairstruct)+1;
                    Pairs(tyc).(pairstruct)(n_pairs.(pairstruct)).unit1=pop_corr(u1,tyc).unit_ID;
                    Pairs(tyc).(pairstruct)(n_pairs.(pairstruct)).unit2=pop_corr(u2,tyc).unit_ID;
                    
                            per_epoch1=pop_corr(u1,tyc).FR(ismember([pop_corr(u1,tyc).block],u12_blocks),:);
                            per_epoch2=pop_corr(u2,tyc).FR(ismember([pop_corr(u2,tyc).block],u12_blocks),:);
                     for e=1:size(keys.EPOCHS,1)
                       % for c=1:size(condition_matrix_pairs,1)
%                             trpar1(end+1,:)=ismember([pop1.trial.block],u12_blocks);
%                             trpar2(end+1,:)=ismember([pop2.trial.block],u12_blocks);
                            
                            if ~isempty(per_epoch1) && ~isempty(per_epoch2)
                                [RHO,PVAL]=corr(per_epoch1(:,e),per_epoch2(:,e));
                                Pairs(tyc).(pairstruct)(n_pairs.(pairstruct)).epoch(e)=RHO;
                                Pairs(tyc).(pairstruct)(n_pairs.(pairstruct)).signi(e)=double(PVAL<0.05)*sign(RHO);
                                
                            else
                                Pairs(tyc).(pairstruct)(n_pairs.(pairstruct)).epoch(e)=NaN;
                                Pairs(tyc).(pairstruct)(n_pairs.(pairstruct)).signi(e)=NaN;
                            end
                            if strcmp(pop_corr(u1,tyc).tuning{e},'-') || strcmp(pop_corr(u2,tyc).tuning{e},'-') ||...
                               strcmp(pop_corr(u1,tyc).tuning{e},'')  || strcmp(pop_corr(u2,tyc).tuning{e},'');
                               Pairs(tyc).(pairstruct)(n_pairs.(pairstruct)).tuning(e)=0;
                            elseif strcmp(pop_corr(u1,tyc).tuning{e},pop_corr(u2,tyc).tuning{e});
                               Pairs(tyc).(pairstruct)(n_pairs.(pairstruct)).tuning(e)=1;
                            else
                               Pairs(tyc).(pairstruct)(n_pairs.(pairstruct)).tuning(e)=-1;                                
                            end
                            %end
                     end
                end
            end
        end
        end
    end
end


fig_title=sprintf('Selection: %s %s %s hnd %s ch %s %s',...
    [Sel_for_title{:}],keys.monkey,[keys.conditions_to_plot{:}],mat2str(u_hands),mat2str(double(u_choice)),keys.position_and_plotting_arrangements{a});
filename=sprintf('%s %s %s %s hnd %s ch %s c%d',...
    [Sel_for_title{:}],keys.monkey,keys.population_group_parameter,[keys.conditions_to_plot{:}],mat2str(u_hands),mat2str(double(u_choice)),a);
plot_1_title            = [fig_title  ' Correlation'];

Correlation_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_1_title);

pair_col_struct={'ipsi','ipsi','cont','cont'};
pair_col_tuning=[1,-1,-1,1]; %% 'same' tuning (f.e. ipsi space) for units in opposite hemispheres is actually opposite
pair_col_titles={'same hemisphere same tuning','same hemisphere opposite tuning','opposite hemisphere same tuning','opposite hemisphere opposite tuning'};
all_epochs_to_show={'Fhol','Cue','EDel','Del','PreS','PeriS','PreR','PeriR'};

% keys.saccade_states={'PreS','PeriS','PostS'};
% keys.reach_states={'PreR','PeriR','PostR'};

clear N_units N_positive N_negative mean_corr sterr_corr


for tyc=1:numel(Pairs)
    tyc_epochs_to_show=all_epochs_to_show;
    if any(strfind(tyc_label{tyc},'Ddre'));
      tyc_epochs_to_show(ismember(tyc_epochs_to_show,keys.saccade_states))=[];
    elseif any(strfind(tyc_label{tyc},'Ddsa'));
      tyc_epochs_to_show(ismember(tyc_epochs_to_show,keys.reach_states))=[];        
    end
    
    epochs_to_show=find(ismember(keys.EPOCHS(:,1),tyc_epochs_to_show))';
    
    for c=1:numel(pair_col_struct)
        if ~isfield(Pairs(tyc),pair_col_struct{c})
            continue
        end
            
        current_correlations=vertcat(Pairs(tyc).(pair_col_struct{c}).epoch);
        current_significances=vertcat(Pairs(tyc).(pair_col_struct{c}).signi);
        current_tuning=vertcat(Pairs(tyc).(pair_col_struct{c}).tuning)==pair_col_tuning(c);
        asterisks=repmat({''},numel(epochs_to_show),1);
        for es=1:numel(epochs_to_show)
            e=epochs_to_show(es);
            mean_signi(es)=ttest(current_correlations(current_tuning(:,e),e));
            if mean_signi(es)==1
               asterisks{es}='*'; 
            end
            mean_corr(es)=nanmean(current_correlations(current_tuning(:,e),e));
            sterr_corr(es)=sterr(current_correlations(current_tuning(:,e),e));
            N_units(tyc,c,es)=sum(current_tuning(:,e));
            N_positive(tyc,c,es)=sum(current_significances(current_tuning(:,e),e)==1);
            N_negative(tyc,c,es)=sum(current_significances(current_tuning(:,e),e)==-1);
        end
    csp(tyc,c)=subplot(numel(Pairs),numel(pair_col_struct),(tyc-1)*numel(pair_col_struct)+c);
    hold on
    %bar(mean_corr);
    errorbar_constant_barsize_working(1:numel(mean_corr),mean_corr,sterr_corr,sterr_corr,0.2,'o')    
    text(1:numel(epochs_to_show),mean_corr+sterr_corr,asterisks,'color','k')
    set(gca,'xtick',1:numel(mean_corr),'xticklabel',keys.EPOCHS(epochs_to_show,1))
    if tyc==1;
    title(pair_col_titles{c})
    end
    if c==1
    ylabel([tyc_label(tyc); 'correlation coefficient']);
    end
    current_ylim=get(gca,'ylim');
    y_min(tyc,c)=current_ylim(1);
    y_max(tyc,c)=current_ylim(2);
    end
end
% ylim overwritten... (add key)
% y_lim=[min(min([y_min y_max])) max(max([y_min y_max]))];
% y_lim(2)=y_lim(2)+diff(y_lim)*0.2;
y_lim=[-0.1 0.25];
for tyc=1:size(csp,1)
    for c=1:size(csp,2)
        subplot(csp(tyc,c))
        set(gca,'ylim',y_lim);
        current_U=N_units(tyc,c,:);
        current_P=N_positive(tyc,c,:);
        current_N=N_negative(tyc,c,:);
        text(1:size(N_units,3),repmat(y_lim(2)-diff(y_lim)*0.05,size(N_units,3),1),num2str(current_U(:)),'color','k')
        text(1:size(N_units,3),repmat(y_lim(2)-diff(y_lim)*0.1,size(N_units,3),1),num2str(current_N(:)),'color','b')
        text(1:size(N_units,3),repmat(y_lim(2)-diff(y_lim)*0.15,size(N_units,3),1),num2str(current_P(:)),'color','r')
    end
end
 title_and_save(Correlation_summary_handle,  [filename ' PSTHs'],plot_1_title,keys)
end


function title_and_save(figure_handle,filename,plot_title,keys)
mtit(figure_handle,  plot_title, 'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none');
stampit;
if keys.plot.export
    export_fig([keys.path_to_save, filename], '-pdf','-transparent') % pdf by run
    close all
end
end
