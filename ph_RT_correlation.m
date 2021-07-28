function ph_RT_correlation(population,modified_keys)
warning('off','MATLAB:catenate:DimensionMismatch');%
% %%% !!!!!!!! make sure epochs (baseline, cue, .... make sense for all effectors)

for fn=fieldnames(modified_keys)'
    keys.(fn{:})=modified_keys.(fn{:});
end

legend_labels_hf={'NH IS IN' 'NH IS CH' 'IH IS IN' 'IH IS CH' 'CH IS IN' 'CH IS CH' ...
    'NH VS IN' 'NH VS CH' 'IH VS IN' 'IH VS CH' 'CH VS IN' 'CH VS CH' ...
    'NH CS IN' 'NH CS CH' 'IH CS IN' 'IH CS CH' 'CH CS IN' 'CH CS CH' ...
    'NH IS IN P' 'NH IS CH P' 'IH IS IN P' 'IH IS CH P' 'CH IS IN P' 'CH IS CH P'...
    'NH VS IN P' 'NH VS CH P' 'IH VS IN P' 'IH VS CH P' 'CH VS IN P' 'CH VS CH P'...
    'NH CS IN P' 'NH CS CH P' 'IH CS IN P' 'IH CS CH P' 'CH CS IN P' 'CH CS CH P'};
legend_labels_pref={
    'NH NP IN' 'NH NP CH' 'IH NP IN' 'IH NP CH' 'CH NP IN' 'CH NP CH' ...
    'NH PF IN' 'NH PF CH' 'IH PF IN' 'IH PF CH' 'CH PF IN' 'CH PF CH' ...
    'NH NP IN P' 'NH NP CH P' 'IH NP IN P' 'IH NP CH P' 'CH NP IN P' 'CH NP CH P'...
    'NH PF IN P' 'NH PF CH P' 'IH PF IN P' 'IH PF CH P' 'CH PF IN P' 'CH PF CH P'};
legend_labels_pos={
    'NH IN' 'NH CH' 'IH IN' 'IH CH' 'CH IN' 'CH CH' ...
    'NH IN P' 'NH CH P' 'IH IN P' 'IH CH P' 'CH IN P' 'CH CH P'};

cols=keys.colors;

%% there are just too many colors once we include vertical targets, so for now we just keep use the same ones again...
keys.line_colors=[[cols.NH_IS_IN;cols.NH_IS_CH;cols.IH_IS_IN;cols.IH_IS_CH;cols.CH_IS_IN;cols.CH_IS_CH;...
    cols.NH_VS_IN;cols.NH_VS_CH;cols.IH_VS_IN;cols.IH_VS_CH;cols.CH_VS_IN;cols.CH_VS_CH;...
    cols.NH_CS_IN;cols.NH_CS_CH;cols.IH_CS_IN;cols.IH_CS_CH;cols.CH_CS_IN;cols.CH_CS_CH;]/255;...
    [cols.NH_IS_IN;cols.NH_IS_CH;cols.IH_IS_IN;cols.IH_IS_CH;cols.CH_IS_IN;cols.CH_IS_CH;...
    cols.NH_VS_IN;cols.NH_VS_CH;cols.IH_VS_IN;cols.IH_VS_CH;cols.CH_VS_IN;cols.CH_VS_CH;...
    cols.NH_CS_IN;cols.NH_CS_CH;cols.IH_CS_IN;cols.IH_CS_CH;cols.CH_CS_IN;cols.CH_CS_CH;]/510]; %%temporary for inactivation
keys.pref_colors=[[cols.NH_IS_IN;cols.NH_IS_CH;cols.IH_IS_IN;cols.IH_IS_CH;cols.CH_IS_IN;cols.CH_IS_CH;...
    cols.NH_CS_IN;cols.NH_CS_CH;cols.IH_CS_IN;cols.IH_CS_CH;cols.CH_CS_IN;cols.CH_CS_CH;]/255;...
    [cols.NH_IS_IN;cols.NH_IS_CH;cols.IH_IS_IN;cols.IH_IS_CH;cols.CH_IS_IN;cols.CH_IS_CH;...
    cols.NH_CS_IN;cols.NH_CS_CH;cols.IH_CS_IN;cols.IH_CS_CH;cols.CH_CS_IN;cols.CH_CS_CH;]/510]; %%temporary for inactivation
keys.pos_colors=[[cols.NH_IN;cols.NH_CH;cols.IH_IN;cols.IH_CH;cols.CH_IN;cols.CH_CH]/255;...
    [cols.NH_IN;cols.NH_CH;cols.IH_IN;cols.IH_CH;cols.CH_IN;cols.CH_CH]/510]; %%temporary for inactivation

%% tuning table preparation and grouping
[tuning_per_unit_table]                 = ph_load_extended_tuning_table(keys);
[tuning_per_unit_table, Sel_for_title]  = ph_reduce_tuning_table(tuning_per_unit_table,keys);
idx_group_parameter=DAG_find_column_index(tuning_per_unit_table,keys.RT.group_parameter);
idx_unitID=DAG_find_column_index(tuning_per_unit_table,'unit_ID');
idx_RF_frame=DAG_find_column_index(tuning_per_unit_table,keys.RT.RF_frame_parameter);
group_values=tuning_per_unit_table(:,idx_group_parameter);
group_values=cellfun(@num2str, group_values, 'UniformOutput', false);
cell_in_any_group=[false; ~ismember(group_values(2:end),keys.RT.group_excluded)];
unique_group_values=unique(group_values(cell_in_any_group));
if isempty(unique_group_values)
    disp('no relevant groups found');
    return;
end
tuning_per_unit_table=tuning_per_unit_table(cell_in_any_group,:);
group_values=group_values(cell_in_any_group);
complete_unit_list={population.unit_ID}';
population=population(ismember(complete_unit_list,tuning_per_unit_table(:,idx_unitID)));
%complete_unit_list={population.unit_ID}';


all_trialz=[population.trial];
per_trial.types       =[all_trialz.type];
per_trial.effectors   =[all_trialz.effector];
u_con.type     =unique(per_trial.types);
u_con.effector =unique(per_trial.effectors);
all_type_effectors      = combvec(u_con.type,u_con.effector)';
type_effectors =[];

% redifine type_effectors to include only relevant
for t=1:size(all_type_effectors,1)
    typ=all_type_effectors(t,1);
    eff=all_type_effectors(t,2);
    [~, type_effector_short{t}]=MPA_get_type_effector_name(typ,eff);
    if ~ismember(type_effector_short{t},keys.conditions_to_plot) %|| sum(tr_con)<1
        continue;
    end
    type_effectors=[type_effectors; all_type_effectors(t,:)];
end
type_effector_short(~ismember(type_effector_short,keys.conditions_to_plot))=[];
u_con.type     =unique(type_effectors(:,1))';
u_con.effector =unique(type_effectors(:,2))';


%% define conditions to look at
all_trialz=[population.trial];
per_trial.types       =[all_trialz.type];
per_trial.effectors   =[all_trialz.effector];

tr_con=ismember([all_trialz.completed],keys.cal.completed);
[whatisthis]=ph_arrange_positions_and_plots(keys,all_trialz(tr_con));

condition_parameters  ={'reach_hand','choice','perturbation'};
per_trial.types       =[all_trialz.type];
per_trial.effectors   =[all_trialz.effector];
per_trial.hands       =[all_trialz.reach_hand];
per_trial.choice      =[all_trialz.choice];
per_trial.perturbation=[all_trialz.perturbation];
per_trial.hemifield   =[whatisthis.trial.hemifield];
per_trial.perturbation(ismember(per_trial.perturbation, keys.cal.perturbation_groups{1}))=0;
per_trial.perturbation(ismember(per_trial.perturbation, keys.cal.perturbation_groups{2}))=1;

u_con.hemifield=unique(per_trial.hemifield); %[-1,0,1]; % why does this have to be hardcoded? ---> Because case not defined yet, case defines positions !!

u_con.reach_hand     =unique(per_trial.hands);
u_con.choice    =unique(per_trial.choice);
u_con.perturbation    =unique(per_trial.perturbation);
u_con.perturbation=u_con.perturbation(~isnan(u_con.perturbation));

%% limit conditions key?
if ~any(keys.tt.hands==0) % cause hands 0 is any hand
    u_con.reach_hand     =u_con.reach_hand(ismember(u_con.reach_hand,keys.tt.hands));
end
u_con.choice    =u_con.choice(ismember(u_con.choice,keys.tt.choices));

% reduce trials to only valid
unit_valid=true(size(population));
for u=1:numel(population)
    poptr=population(u).trial;
    valid=ismember([poptr.effector],u_con.effector) & ismember([poptr.type],u_con.type) & ismember([poptr.choice],u_con.choice) & ismember([poptr.reach_hand],u_con.reach_hand);
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

%% defining set of conditions dynmically! --> use this as input for ph_condition_normalization as well
u_condition_definitions={'effector','reach_hand','choice','perturbation'};
for c=1:numel(u_condition_definitions)
    u_condition_parameters{c} = u_con.(u_condition_definitions{c});
end

u_condition_hf_definitions={'hemifield','effector','reach_hand','choice','perturbation'}; %rename
conditions_out            = combvec(u_condition_parameters{:})';
%condition_matrix            = combvec(u_con.reach_hand,u_con.choice, u_con.perturbation,u_con.hemifield)';
%conditions_out              = combvec(u_con.effector,u_con.reach_hand,u_con.choice, u_con.perturbation)';
conditions_hf               = combvec(u_con.hemifield,conditions_out')';
conditions_hf_complete      = combvec(u_con.hemifield,conditions_out')';
conditions_pref             = combvec([0 1],conditions_out')';

if any(u_con.reach_hand~=0) && any(u_con.perturbation==1) %splitting to all 4 hand space conditions if hands are involved
    [~,~,columns_hf] = unique(conditions_hf(:,[1,3]),'rows');
    [~,~,columns_pref] = unique(conditions_pref(:,[1,3]),'rows');
else
    columns_hf          = ones(size(conditions_hf,1),1);
    columns_pref        = ones(size(conditions_pref,1),1);
end




%% finding positions and fixations
positions=unique(vertcat(whatisthis.trial.position),'rows');
keys.normalization_field='RT';
[~, ~, condition_per_trial, ~]=ph_condition_normalization(population,keys);


typ=4; %% not for all types yet
condition_per_hf=[condition_per_trial.per_hemifield]; %check order!!!
ep=find(ismember(keys.EPOCHS_PER_TYPE{typ}(:,1),keys.RT.epoch_RT));

for c=1:numel(condition_per_hf)
    for u=1:numel(condition_per_hf(c).unit)
        [coef, pval] =corr(condition_per_hf(c).unit(u).epoch_FRs(:,ep),condition_per_hf(c).unit(u).sac_lat');
        condition_per_hf(c).unit(u).sac_corr=coef;
        condition_per_hf(c).unit(u).sac_p=pval;
        condition_per_hf(c).unit(u).sac_sig=pval<0.05;
        
        [coef, pval] =corr(condition_per_hf(c).unit(u).epoch_FRs(:,ep),condition_per_hf(c).unit(u).rea_lat');
        condition_per_hf(c).unit(u).rea_corr=coef;
        condition_per_hf(c).unit(u).rea_p=pval;
        condition_per_hf(c).unit(u).rea_sig=pval<0.05;
    end
end

%% plots



%find conditions for different effectors
con_comb=find(ismember(conditions_hf(:,ismember(u_condition_hf_definitions,'effector')),6)); %% could also be [2,6] in the future
con_sac=find(ismember(conditions_hf(:,ismember(u_condition_hf_definitions,'effector')),3)); %% could also be [0,3] in the future
con_rea=find(ismember(conditions_hf(:,ismember(u_condition_hf_definitions,'effector')),4)); %% could also be [1,4] in the future

n_conditions=numel(con_comb);

%% summary of correlation values dissociated only
plot_title              = ['RT correlations with FR in ' keys.RT.epoch_RT ', one effector'];
figure_handle              = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);

for c=1:n_conditions
    label_index=(conditions_hf(con_comb(c),1)+1)*6+conditions_hf(con_comb(c),3)*2+1;
    subplot(ceil(sqrt(n_conditions)),ceil(sqrt(n_conditions)),c)
    
    hold on
    plot([-1 1],[-1 1],'k');
    plot([-1 1],[1 -1],'k');
    corr_1=[condition_per_hf(con_sac(c)).unit.sac_corr];
    corr_2=[condition_per_hf(con_rea(c)).unit.rea_corr];        
    sig_1=[condition_per_hf(con_sac(c)).unit.sac_sig];
    sig_2=[condition_per_hf(con_rea(c)).unit.rea_sig];
    sig_c=sig_1+2*sig_2;
    
    for sig=0:3
        u=sig_c==sig;
        switch sig
            case 0
                col=[0.5 0.5 0.5];
            case 1
                col='r';
            case 2
                col='g';
            case 3
                col=[0 0 0];
        end
        scatter(corr_1(u),corr_2(u),10,col,'filled');
    end
    n0=sum(sig_c==0);
    n1=sum(sig_c==1);
    n2=sum(sig_c==2);
    n3=sum(sig_c==3);
    text(-0.9,-0.9,['N=' num2str(n0)],'color',[0.5 0.5 0.5]);
    text(-0.9,-0.8,['N=' num2str(n1)],'color','r');
    text(-0.9,-0.7,['N=' num2str(n2)],'color','g');
    text(-0.9,-0.6,['N=' num2str(n3)],'color',[0 0 0]);
    title(legend_labels_hf{label_index});
    xlim([-1.001 1.001]);
    ylim([-1.001 1.001]);
    axis square
    xlabel('Saccade');
    ylabel('Reach');
end
ph_title_and_save(figure_handle,plot_title,plot_title,keys)

%% summary of correlation values dissociated only
plot_title              = ['RT correlations with FR in ' keys.RT.epoch_RT ', both effectors'];
figure_handle              = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);

for c=1:n_conditions
    label_index=(conditions_hf(con_comb(c),1)+1)*6+conditions_hf(con_comb(c),3)*2+1;
    subplot(ceil(sqrt(n_conditions)),ceil(sqrt(n_conditions)),c)
    
    hold on
    plot([-1 1],[-1 1],'k');
    plot([-1 1],[1 -1],'k');
    corr_1=[condition_per_hf(con_comb(c)).unit.sac_corr];
    corr_2=[condition_per_hf(con_comb(c)).unit.rea_corr];        
    sig_1=[condition_per_hf(con_comb(c)).unit.sac_sig];
    sig_2=[condition_per_hf(con_comb(c)).unit.rea_sig];
    sig_c=sig_1+2*sig_2;
    
    for sig=0:3
        u=sig_c==sig;
        switch sig
            case 0
                col=[0.5 0.5 0.5];
            case 1
                col='r';
            case 2
                col='g';
            case 3
                col=[0 0 0];
        end
        scatter(corr_1(u),corr_2(u),10,col,'filled');
    end
    n0=sum(sig_c==0);
    n1=sum(sig_c==1);
    n2=sum(sig_c==2);
    n3=sum(sig_c==3);
    text(-0.9,-0.9,['N=' num2str(n0)],'color',[0.5 0.5 0.5]);
    text(-0.9,-0.8,['N=' num2str(n1)],'color','r');
    text(-0.9,-0.7,['N=' num2str(n2)],'color','g');
    text(-0.9,-0.6,['N=' num2str(n3)],'color',[0 0 0]);
    title(legend_labels_hf{label_index});
    xlim([-1.001 1.001]);
    ylim([-1.001 1.001]);
    axis square
    xlabel('Saccade');
    ylabel('Reach');
end
ph_title_and_save(figure_handle,plot_title,plot_title,keys)

%% summary of correlation values Saccades
plot_title              = ['RT correlations with FR in ' keys.RT.epoch_RT ', Saccades'];
figure_handle              = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);

for c=1:n_conditions
    label_index=(conditions_hf(con_comb(c),1)+1)*6+conditions_hf(con_comb(c),3)*2+1;
    subplot(ceil(sqrt(n_conditions)),ceil(sqrt(n_conditions)),c)
    
    hold on
    plot([-1 1],[-1 1],'k');
    plot([-1 1],[1 -1],'k');
    corr_1=[condition_per_hf(con_sac(c)).unit.sac_corr];
    corr_2=[condition_per_hf(con_comb(c)).unit.sac_corr];        
    sig_1=[condition_per_hf(con_sac(c)).unit.sac_sig];
    sig_2=[condition_per_hf(con_comb(c)).unit.sac_sig];
    sig_c=sig_1+2*sig_2;
    
    for sig=0:3
        u=sig_c==sig;
        switch sig
            case 0
                col=[0.5 0.5 0.5];
            case 1
                col='r';
            case 2
                col='b';
            case 3
                col=[0 0 0];
        end
        scatter(corr_1(u),corr_2(u),10,col,'filled');
    end
    n0=sum(sig_c==0);
    n1=sum(sig_c==1);
    n2=sum(sig_c==2);
    n3=sum(sig_c==3);
    text(-0.9,-0.9,['N=' num2str(n0)],'color',[0.5 0.5 0.5]);
    text(-0.9,-0.8,['N=' num2str(n1)],'color','r');
    text(-0.9,-0.7,['N=' num2str(n2)],'color','b');
    text(-0.9,-0.6,['N=' num2str(n3)],'color',[0 0 0]);
    title(legend_labels_hf{label_index});
    xlim([-1.001 1.001]);
    ylim([-1.001 1.001]);
    axis square
    xlabel('Saccade only');
    ylabel('Combined');
end
ph_title_and_save(figure_handle,plot_title,plot_title,keys)

%% summary of correlation values Reaches
plot_title              = ['RT correlations with FR in ' keys.RT.epoch_RT ', Reaches'];
figure_handle              = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);

for c=1:n_conditions
    label_index=(conditions_hf(con_comb(c),1)+1)*6+conditions_hf(con_comb(c),3)*2+1;
    subplot(ceil(sqrt(n_conditions)),ceil(sqrt(n_conditions)),c)
    
    hold on
    plot([-1 1],[-1 1],'k');
    plot([-1 1],[1 -1],'k');
    corr_1=[condition_per_hf(con_rea(c)).unit.rea_corr];
    corr_2=[condition_per_hf(con_comb(c)).unit.rea_corr];        
    sig_1=[condition_per_hf(con_rea(c)).unit.rea_sig];
    sig_2=[condition_per_hf(con_comb(c)).unit.rea_sig];
    sig_c=sig_1+2*sig_2;
    
    for sig=0:3
        u=sig_c==sig;
        switch sig
            case 0
                col=[0.5 0.5 0.5];
            case 1
                col='g';
            case 2
                col='b';
            case 3
                col=[0 0 0];
        end
        scatter(corr_1(u),corr_2(u),10,col,'filled');
    end
    n0=sum(sig_c==0);
    n1=sum(sig_c==1);
    n2=sum(sig_c==2);
    n3=sum(sig_c==3);
    text(-0.9,-0.9,['N=' num2str(n0)],'color',[0.5 0.5 0.5]);
    text(-0.9,-0.8,['N=' num2str(n1)],'color','g');
    text(-0.9,-0.7,['N=' num2str(n2)],'color','b');
    text(-0.9,-0.6,['N=' num2str(n3)],'color',[0 0 0]);
    title(legend_labels_hf{label_index});
    xlim([-1.001 1.001]);
    ylim([-1.001 1.001]);
    axis square
    xlabel('Reaches only');
    ylabel('Combined');
end
ph_title_and_save(figure_handle,plot_title,plot_title,keys)

