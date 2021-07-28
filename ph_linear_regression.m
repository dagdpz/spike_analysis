function ph_linear_regression(population,modified_keys)
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
idx_group_parameter=DAG_find_column_index(tuning_per_unit_table,keys.RE.group_parameter);
idx_unitID=DAG_find_column_index(tuning_per_unit_table,'unit_ID');
idx_RF_frame=DAG_find_column_index(tuning_per_unit_table,keys.RE.RF_frame_parameter);
group_values=tuning_per_unit_table(:,idx_group_parameter);
group_values=cellfun(@num2str, group_values, 'UniformOutput', false);
cell_in_any_group=[false; ~ismember(group_values(2:end),keys.RE.group_excluded)];
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
keys.normalization_field='RE';
keys.WINDOWS_PER_TYPE=keys.RE.WINDOWS_PER_TYPE;
if true
[~, condition,~,pref_valid]=ph_condition_normalization(population,keys);

save([keys.path_to_save, filesep, 'bootstrapped_normalized'],'keys','condition','complete_unit_list');
else
load([keys.path_to_save, filesep, 'bootstrapped_normalized'],'keys','condition','complete_unit_list');
end


%% regression condition matrix
regressors=keys.RE.regressors;
solution=keys.RE.solution;


% for g=1:numel(unique_group_values)
%     for t=1:size(condition,1)
%     end
% end

for r=1:numel(regressors)
    for c=1:numel(u_condition_hf_definitions)
        con4reg_tmp(:,c) = ismember(conditions_hf(:,c),regressors(r).(u_condition_hf_definitions{c}));
    end
    con4reg(:,r)=all(con4reg_tmp,2);
end
for c=1:numel(u_condition_hf_definitions)
    con4reg_tmp(:,c) =  ismember(conditions_hf(:,c),solution.(u_condition_hf_definitions{c}));
end

con4reg(:,r+1)=all(con4reg_tmp,2);

condition_per_hf=[condition.per_hemifield]; %check order!!!


%% for movement period, calculate difference of saccade and reach reaction times (per condition!),
%% cut off the respective amount of bins from reaches and combined,
%% and align saccades"

%find conditions for different effectors
con_comb=find(ismember(conditions_hf(:,ismember(u_condition_hf_definitions,'effector')),6)); %% could also be [2,6] in the future
wR=find(ismember([keys.WINDOWS_PER_TYPE{typ}(:,1)],'Reach')); %% typ!!
wS=find(ismember([keys.WINDOWS_PER_TYPE{typ}(:,1)],'Saccade')); %% typ!!
wM=find(ismember([keys.WINDOWS_PER_TYPE{typ}(:,1)],'Movement')); %% typ!!
con_sac=find(ismember(conditions_hf(:,ismember(u_condition_hf_definitions,'effector')),3)); %% could also be [0,3] in the future
con_rea=find(ismember(conditions_hf(:,ismember(u_condition_hf_definitions,'effector')),4)); %% could also be [1,4] in the future


tmpstrct=vertcat(condition_per_hf.window);
tmpstrct=tmpstrct(sub2ind(size(tmpstrct),con_comb,repmat(wM,size(con_comb))));
tmpstrct=[tmpstrct.unit];
sac_rt_comb=[tmpstrct.sac_lat];
rea_rt_comb=[tmpstrct.rea_lat];
[~,idx]=max(abs(rea_rt_comb-sac_rt_comb));
max_shift=rea_rt_comb(idx)-sac_rt_comb(idx);
max_shift_in_bin=round(double(abs(max_shift/keys.PSTH_binwidth)));
n_conditions=sum(con4reg(:,1));

for u=1:numel(condition_per_hf(c).window(wM).unit)
    for g=1:n_conditions
        reg_con_tmp=find(con4reg(:,end));
        c=reg_con_tmp(g);
        sac_RT=condition_per_hf(c).window(wM).unit(u).sac_lat;
        rea_RT=condition_per_hf(c).window(wM).unit(u).rea_lat;
        shift=round((rea_RT-sac_RT)/keys.PSTH_binwidth);
        for r=1:size(con4reg,2)
            reg_con_tmp=find(con4reg(:,r));
            c=reg_con_tmp(g);
            if ismember(c,con_sac)
                %% remove bins one
                if max_shift>0
                    condition_per_hf(c).window(wM).unit(u).average_spike_density = condition_per_hf(c).window(wM).unit(u).average_spike_density(shift+1:end-max_shift_in_bin+shift);
%                     condition_per_hf(c).window(w).unit(u).SD_for_corr           = condition_per_hf(c).window(wS).unit(u).average_spike_density(shift+1:end-max_shift_in_bin+shift);
%                     condition_per_hf(c).window(w).unit(u).bootstrapped          = condition_per_hf(c).window(wS).unit(u).bootstrapped(:,shift+1:end-max_shift_in_bin+shift);
                elseif max_shift<0
                    condition_per_hf(c).window(wM).unit(u).average_spike_density = condition_per_hf(c).window(wM).unit(u).average_spike_density(max_shift_in_bin+shift+1:end+shift);
%                     condition_per_hf(c).window(w).unit(u).SD_for_corr           = condition_per_hf(c).window(wS).unit(u).average_spike_density(max_shift_in_bin+shift+1:end+shift);
%                     condition_per_hf(c).window(w).unit(u).bootstrapped          = condition_per_hf(c).window(wS).unit(u).bootstrapped(:,max_shift_in_bin+shift+1:end+shift); %% this is assuming shift is negative for this condition/unit as well
                end
            elseif ismember(c,con_rea)
                if max_shift>0
                    condition_per_hf(c).window(wM).unit(u).average_spike_density = condition_per_hf(c).window(wM).unit(u).average_spike_density(1:end-max_shift_in_bin);
%                     condition_per_hf(c).window(w).unit(u).SD_for_corr           = condition_per_hf(c).window(wR).unit(u).average_spike_density(1:end-max_shift_in_bin);
%                     condition_per_hf(c).window(w).unit(u).bootstrapped          = condition_per_hf(c).window(wR).unit(u).bootstrapped(:,1:end-max_shift_in_bin);
                elseif max_shift<0
                    condition_per_hf(c).window(wM).unit(u).average_spike_density = condition_per_hf(c).window(wM).unit(u).average_spike_density(max_shift_in_bin+1:end);
%                     condition_per_hf(c).window(w).unit(u).SD_for_corr           = condition_per_hf(c).window(wR).unit(u).average_spike_density(max_shift_in_bin+1:end);
%                     condition_per_hf(c).window(w).unit(u).bootstrapped          = condition_per_hf(c).window(wR).unit(u).bootstrapped(:,max_shift_in_bin+1:end);
                end
            else 
                if max_shift>0
                    condition_per_hf(c).window(wM).unit(u).average_spike_density = condition_per_hf(c).window(wM).unit(u).average_spike_density(1:end-max_shift_in_bin);
%                     condition_per_hf(c).window(w).unit(u).SD_for_corr           = condition_per_hf(c).window(w).unit(u).average_spike_density(1:end-max_shift_in_bin);
%                     condition_per_hf(c).window(w).unit(u).bootstrapped          = condition_per_hf(c).window(w).unit(u).bootstrapped(:,1:end-max_shift_in_bin);
                elseif max_shift<0
                    condition_per_hf(c).window(wM).unit(u).average_spike_density = condition_per_hf(c).window(wM).unit(u).average_spike_density(max_shift_in_bin+1:end);
%                     condition_per_hf(c).window(w).unit(u).SD_for_corr           = condition_per_hf(c).window(w).unit(u).average_spike_density(max_shift_in_bin+1:end);
%                     condition_per_hf(c).window(w).unit(u).bootstrapped          = condition_per_hf(c).window(w).unit(u).bootstrapped(:,max_shift_in_bin+1:end);
                end
            end
        end
    end
end



%% revert by unit ordering
for c=1:numel(condition_per_hf)
    for w=1:numel(condition_per_hf(c).window)
        for u=1:numel(condition_per_hf(c).window(w).unit)
            PSTHs_per_unit(u).window(w).condition(c)= condition_per_hf(c).window(w).unit(u);
        end
    end
end

%% do the actual regression - some mismatch here??
for u=1:numel(PSTHs_per_unit)
    for w=1:numel(PSTHs_per_unit(u).window)
        clear reg4fit reg4fitSEM reg_perbinSEM reg_tmp sol_tmp reg_perbin reg_bootstrapped sol_bootstrapped corrR1_bt betas
        for r=1:numel(regressors)
            if w==wM && all(ismember(find(con4reg(:,r)),con_sac))
                reg_tmp(:,r)= [PSTHs_per_unit(u).window(wS).condition(con4reg(:,r)).average_spike_density]';
                reg_bootstrapped{r}= [PSTHs_per_unit(u).window(wS).condition(con4reg(:,r)).bootstrapped]';
            elseif w==wM && all(ismember(find(con4reg(:,r)),con_rea))
                reg_tmp(:,r)= [PSTHs_per_unit(u).window(wR).condition(con4reg(:,r)).average_spike_density]';
                reg_bootstrapped{r}= [PSTHs_per_unit(u).window(wR).condition(con4reg(:,r)).bootstrapped]';
            else
                reg_tmp(:,r)= [PSTHs_per_unit(u).window(w).condition(con4reg(:,r)).average_spike_density]';
                reg_bootstrapped{r}= [PSTHs_per_unit(u).window(w).condition(con4reg(:,r)).bootstrapped]';
            end
            
            reg4fit(:,r)= [PSTHs_per_unit(u).window(w).condition(con4reg(:,r)).average_spike_density]';
            reg4fitSEM(:,r)= [PSTHs_per_unit(u).window(w).condition(con4reg(:,r)).SEM_spike_density]';
            reg_perbin{r}= vertcat(PSTHs_per_unit(u).window(w).condition(con4reg(:,r)).average_spike_density);
            reg_perbinSEM{r}= vertcat(PSTHs_per_unit(u).window(w).condition(con4reg(:,r)).SEM_spike_density);
        end
        if w==wM && all(ismember(find(con4reg(:,1)),con_sac))
            sol_bootstrapped{1}=[PSTHs_per_unit(u).window(wS).condition(con4reg(:,end)).bootstrapped]';
            sol_tmp(:,1)       =[PSTHs_per_unit(u).window(wS).condition(con4reg(:,end)).average_spike_density]';
            sol_bootstrapped{2}=[PSTHs_per_unit(u).window(wR).condition(con4reg(:,end)).bootstrapped]';
            sol_tmp(:,2)       =[PSTHs_per_unit(u).window(wR).condition(con4reg(:,end)).average_spike_density]';
        elseif w==wM && all(ismember(find(con4reg(:,1)),con_rea))
            sol_bootstrapped{1}=[PSTHs_per_unit(u).window(wR).condition(con4reg(:,end)).bootstrapped]';
            sol_tmp(:,1)       =[PSTHs_per_unit(u).window(wR).condition(con4reg(:,end)).average_spike_density]';
            sol_bootstrapped{2}=[PSTHs_per_unit(u).window(wS).condition(con4reg(:,end)).bootstrapped]';
            sol_tmp(:,2)       =[PSTHs_per_unit(u).window(wS).condition(con4reg(:,end)).average_spike_density]';
        else
            sol_bootstrapped{1}=[PSTHs_per_unit(u).window(w).condition(con4reg(:,end)).bootstrapped]';
            sol_tmp(:,1)       =[PSTHs_per_unit(u).window(w).condition(con4reg(:,end)).average_spike_density]';
            sol_bootstrapped{2}=[PSTHs_per_unit(u).window(w).condition(con4reg(:,end)).bootstrapped]';
            sol_tmp(:,2)       =[PSTHs_per_unit(u).window(w).condition(con4reg(:,end)).average_spike_density]';
        end
        sol4fit=[PSTHs_per_unit(u).window(w).condition(con4reg(:,end)).average_spike_density]';
        sol_perbin=vertcat(PSTHs_per_unit(u).window(w).condition(con4reg(:,end)).average_spike_density);
        %         %% remove average reg_tmp from both reg_tmp and sol_temp - this is a bs idea, cause it leads to super-correlation
        %         sol_tmp=[PSTHs_per_unit(u).window(w).condition(con4reg(:,end)).average_spike_density]'-mean(reg_tmp,2);
        %         reg_tmp=reg_tmp-repmat(mean(reg_tmp,2),[1,numel(regressors)]);
        [corrR,corrP]=corr(reg4fit);
        [corrR1,corrP1]=corr([reg_tmp(:,1) sol_tmp(:,1)]);
        [corrR2,corrP2]=corr([reg_tmp(:,2) sol_tmp(:,2)]);
        [corrR3,corrP3]=corr([reg_tmp(:,1)+reg_tmp(:,2) sol_tmp(:,2)]);
        
        per_window(w).corr1(u)=corrR1(1,2);
        per_window(w).corr2(u)=corrR2(1,2);
        per_window(w).corr3(u)=corrR3(1,2);
        for bt=1:size(sol_bootstrapped{1},2)
            corrR1_tmp=corr([reg_bootstrapped{1}(:,bt) sol_bootstrapped{1}(:,bt)]);
            corrR2_tmp=corr([reg_bootstrapped{2}(:,bt) sol_bootstrapped{2}(:,bt)]);
            corrR3_tmp=corr([reg_bootstrapped{1}(:,bt)+reg_bootstrapped{2}(:,bt) sol_bootstrapped{2}(:,bt)]); %% what to do for movement window here??
            corrR1_bt(bt)=corrR1_tmp(1,2);
            corrR2_bt(bt)=corrR2_tmp(1,2);
            corrR3_bt(bt)=corrR3_tmp(1,2);
        end
        per_window(w).corrCI1(u,:)=prctile(corrR1_bt,[2.5 97.5]);
        per_window(w).corrCI2(u,:)=prctile(corrR2_bt,[2.5 97.5]);
        per_window(w).corrCI3(u,:)=prctile(corrR3_bt,[2.5 97.5]);
        
        if per_window(w).corr1(u)>per_window(w).corrCI2(u,2) && per_window(w).corr2(u)<per_window(w).corrCI1(u,1) && per_window(w).corrCI1(u,1)>0 %% reg1 larger
            per_window(w).corrCIsig(u)=1;
        elseif per_window(w).corr2(u)>per_window(w).corrCI1(u,2) && per_window(w).corr1(u)<per_window(w).corrCI2(u,1) && per_window(w).corrCI2(u,1)>0 %% reg2 larger
            per_window(w).corrCIsig(u)=2;
        elseif per_window(w).corrCI1(u,1)>0 && per_window(w).corrCI2(u,1)>0 %% both significant, but none significantly larger
            per_window(w).corrCIsig(u)=3;
        else
            per_window(w).corrCIsig(u)=0;
        end
        
        if per_window(w).corr1(u)>per_window(w).corrCI3(u,2) && per_window(w).corr3(u)<per_window(w).corrCI1(u,1) && per_window(w).corrCI1(u,1)>0 %% reg1 larger
            per_window(w).corrCIsig1(u)=1;
        elseif per_window(w).corr3(u)>per_window(w).corrCI1(u,2) && per_window(w).corr1(u)<per_window(w).corrCI3(u,1) && per_window(w).corrCI3(u,1)>0 %% sol larger
            per_window(w).corrCIsig1(u)=2;
        elseif per_window(w).corrCI1(u,1)>0 && per_window(w).corrCI3(u,1)>0 %% both significant, but none significantly larger
            per_window(w).corrCIsig1(u)=3;
        else
            per_window(w).corrCIsig1(u)=0;
        end
        
        
        if per_window(w).corr3(u)>per_window(w).corrCI2(u,2) && per_window(w).corr2(u)<per_window(w).corrCI3(u,1) && per_window(w).corrCI3(u,1)>0 %% sol larger
            per_window(w).corrCIsig2(u)=1;
        elseif per_window(w).corr2(u)>per_window(w).corrCI3(u,2) && per_window(w).corr3(u)<per_window(w).corrCI2(u,1) && per_window(w).corrCI2(u,1)>0 %% reg2 larger
            per_window(w).corrCIsig2(u)=2;
        elseif per_window(w).corrCI3(u,1)>0 && per_window(w).corrCI2(u,1)>0 %% both significant, but none significantly larger
            per_window(w).corrCIsig2(u)=3;
        else
            per_window(w).corrCIsig2(u)=0;
        end
        
        for b=1:size(sol_perbin,2)
            if all(~isnan(reg_perbin{1}(:,b))) && all(~isnan(reg_perbin{2}(:,b)))
                [betas,~,~,in] = stepwisefit([reg_perbin{1}(:,b) reg_perbin{2}(:,b)] ,sol_perbin(:,b),'display','off');
%                 
%                 %bl=min(reg_perbin{1}(:,b),reg_perbin{2}(:,b));
%                 bl=(reg_perbin{1}(:,b)+reg_perbin{2}(:,b))/2;
%                 x=reg_perbin{1}(:,b)-bl;
%                 y=reg_perbin{2}(:,b)-bl;
%                 z=sol_perbin(:,b)-bl;
% %                 
% %                 %  %b1       %b2    % intercept
% %                 LB=[0         0     -inf];
% %                 X0=[0         0        0];
% %                 UB=[inf       inf    inf];          
%                 %  %b1       %b2    % intercept
%                 LB=[0         0        0];
%                 X0=[1         1        0];
%                 UB=[2         2      inf];           
%                 fitT=fittype('ph_fit_linear(x,y,b1,b2,IC)','dependent',{'z'},'independent',{'x','y'},'coefficients',{'b1','b2','IC'});
%                 fitopts=fitoptions('method','NonlinearLeastSquares','Lower',LB,'Upper',UB,'StartPoint',X0); % 'Weights' !
%                 
%                 [fitobj, Goodness] =fit([x,y],z,fitT,fitopts);
%                 ci = confint(fitobj);
%                 in=all(ci>0);
%                 betas(1)=fitobj.b1;
%                 betas(2)=fitobj.b2;
                %[~,~,~,in] = stepwisefit([betas(1)*x betas(2)*y] ,z,'display','off','maxiter',1);
                per_window(w).betas1(u,b)=betas(1);
                per_window(w).betas2(u,b)=betas(2);
                per_window(w).sigbins1(u,b)=in(1);
                per_window(w).sigbins2(u,b)=in(2);
            else
                per_window(w).betas1(u,b)=NaN;
                per_window(w).betas2(u,b)=NaN;
                per_window(w).sigbins1(u,b)=0;
                per_window(w).sigbins2(u,b)=0;
            end
        end
        
        per_unit(u).window(w).name= keys.WINDOWS_PER_TYPE{typ}(w,1);
        
        
            similarity_criterion=1/10;
idx=abs(diff(reg4fit,1,2))>similarity_criterion*mean(mean(reg4fit)); %% works only for 2 regressors so far     
        
        if any(any(isnan(reg4fit))) || any(isnan(sol4fit)) || sum(idx)<10
            per_unit(u).window(w).betas=zeros(size(regressors)+1)';
            per_unit(u).window(w).intercept=0;
            per_unit(u).window(w).pval=ones(size(regressors)+1)';
            per_unit(u).window(w).in=[false false false];
            per_unit(u).window(w).SS=NaN;
            per_unit(u).window(w).SSres=NaN;
            per_unit(u).window(w).corrR=corrR;
            per_unit(u).window(w).corrP=corrP;
        else
            
       
            [per_unit(u).window(w).betas,SE,PVAL,in,stats,nextstep,history] = stepwisefit([reg4fit(idx,:) reg4fit(idx,1).*reg4fit(idx,2)],sol4fit(idx),'display','off');
            
            
            
            
%             
%                 bl=(reg4fit(:,1)+reg4fit(:,2))/2;
%                 x=[reg4fit(:,1)-bl]';
%                 y=[reg4fit(:,2)-bl]';
%                 z=[sol4fit-bl]';
% %                 
% %                 %  %b1       %b2    % intercept
% %                 LB=[0         0     -inf];
% %                 X0=[0         0        0];
% %                 UB=[inf       inf    inf];           
% 
%                 
%                 %  %b1       %b2    % intercept
%                 LB=[0         0        0];
%                 X0=[1         1        0];
%                 UB=[2         2      inf];     
%                 fitT=fittype('ph_fit_linear(x,y,b1,b2,IC)','dependent',{'z'},'independent',{'x','y'},'coefficients',{'b1','b2','IC'});
%                 fitopts=fitoptions('method','NonlinearLeastSquares','Lower',LB,'Upper',UB,'StartPoint',X0);
%                 
%                 [fitobj, Goodness] =fit([x',y'],z',fitT,fitopts);
%                 ci = confint(fitobj);
%                 in=all(ci>0);
%                 
%                 betas(1)=fitobj.b1;
%                 betas(2)=fitobj.b2;
%             
%                 [~,SE,PVAL,~,stats,nextstep,history] = stepwisefit([betas(1)*x betas(2)*y] ,z,'display','off','maxiter',1);
            
            per_unit(u).window(w).betas=betas;
            per_unit(u).window(w).pval=PVAL;        %% needs to change
            per_unit(u).window(w).in=in;
            per_unit(u).window(w).SS=stats.SStotal;
            per_unit(u).window(w).SSres=stats.SSresid;  %% not sure if this is the right thing, needed??
            per_unit(u).window(w).corrR=corrR;
            per_unit(u).window(w).corrP=corrP;
            
%             
%             per_unit(u).window(w).intercept=stats.intercept;
%             per_unit(u).window(w).pval=PVAL;
%             per_unit(u).window(w).in=in;
%             per_unit(u).window(w).SS=stats.SStotal;
%             per_unit(u).window(w).SSres=stats.SSresid;
%             per_unit(u).window(w).corrR=corrR;
%             per_unit(u).window(w).corrP=corrP;
        end
    end
end

%% save keys, beta values, complete?
ph_append_to_anova_table(keys,'regression',complete_unit_list,per_unit);

save([keys.path_to_save, filesep, 'regression data'],'keys','per_unit','per_window','complete_unit_list','condition_per_hf');



t=1; % not working across types
g=1; % not working yet for specific groups
typ=4; % not working across types
keys=ph_get_epoch_keys(keys,typ,u_con.effector,sum(type_effectors(:,1)==typ)>1);

wM=find(ismember([keys.WINDOWS_PER_TYPE{typ}(:,1)],'Movement')); %% typ!!
%% adjust window
if max_shift>0
    keys.PSTH_WINDOWS{wM,4}=keys.PSTH_WINDOWS{wM,4}-max_shift_in_bin*keys.PSTH_binwidth;
else
    keys.PSTH_WINDOWS{wM,3}=keys.PSTH_WINDOWS{wM,3}-max_shift_in_bin*keys.PSTH_binwidth;
end


n_windows=numel(per_unit(u).window);

%% plot per bin summary
figure
sph(1)=subplot(2,1,1);
hold on
state_shift=0;
for w=1:n_windows
    clear to_plot
    t_before_state=keys.PSTH_WINDOWS{w,3};
    t_after_state=keys.PSTH_WINDOWS{w,4};
    bins=t_before_state:keys.PSTH_binwidth:t_after_state;
    bins=bins+state_shift-t_before_state;
    
    to_plot(1,:)=sum(per_window(w).sigbins1,1);
    to_plot(2,:)=sum(per_window(w).sigbins2,1);
    to_plot(3,:)=sum(per_window(w).sigbins1 & per_window(w).sigbins2,1);
    plot(bins,to_plot(1,:),'color',regressors(1).color);
    plot(bins,to_plot(2,:),'color',regressors(2).color);
    plot(bins,to_plot(3,:),'color',solution.color);
    
    %
    %     for r=1:size(con4reg,2)-1
    %         per_window(w).betas1(u,b)=betas(1);
    %         per_window(w).betas2(u,b)=betas(2);
    %         per_window(w).sigbins2(u,b)=in(2);
    %
    %         PSTHs_per_unit(u).window(w).condition(reg_con_tmp(c)).average_spike_density*beta(r)+in/2;
    %
    %     end
    %     plot(bins,sum(to_plot),'color','k');
    %     r=r+1;
    %     reg_con_tmp=find(con4reg(:,r));
    %     plot(bins,PSTHs_per_unit(u).window(w).condition(reg_con_tmp(c)).average_spike_density,'color',solution.color)
    state_shift=state_shift+t_after_state-t_before_state+0.1;
    %title(Betas);
end
y_lim(1,:)=get(gca,'ylim');


%% subplot appearance, and tuning lines
for spn=1:numel(sph)
    subplot(sph(spn));
    hold on
    
    %% completed? choices? hands?
    tr=[all_trialz.type]==typ & ismember([all_trialz.completed],keys.cal.completed) &...
        ismember([all_trialz.completed],keys.cal.completed) & ismember([all_trialz.reach_hand],u_con.reach_hand) & ismember([all_trialz.choice],u_con.choice);
    ph_PSTH_background(all_trialz(tr),y_lim(spn,:),y_lim(spn,:),y_lim(spn,:),keys,keys.RE.fontsize_factor)
    
end

plot_title              = ['regression per bin summary'];
ph_title_and_save(gcf,  plot_title,plot_title,keys)

%% summary with beta values
figure;
clear corrP
all_windows=vertcat(per_unit.window);
for w=1:n_windows
    subplot(ceil(sqrt(n_windows)),ceil(sqrt(n_windows)),w)
    hold on
    current_win=all_windows(:,w);
    PVAL=[current_win.pval];    
    betas=[current_win.betas];
    for k=1:numel(current_win)
    corrP(k)=current_win(k).corrP(1,2)<0.05;        
    end
   
    sig1=PVAL(1,:)<0.05;
    sig2=PVAL(2,:)<0.05;
    scatter(betas(1,(sig1|sig2)&corrP),betas(2,(sig1|sig2)&corrP),15,'b');
    scatter(betas(1,~sig1&~sig2&corrP),betas(2,~sig1&~sig2&corrP),15,[0.5 0.5 0.5]);
    scatter(betas(1,~sig1&~sig2&~corrP),betas(2,~sig1&~sig2&~corrP),15,'k');
    scatter(betas(1,sig1&~sig2& ~corrP),betas(2,sig1&~sig2& ~corrP),20,regressors(1).color,'<','filled');
    scatter(betas(1,~sig1&sig2& ~corrP),betas(2,~sig1&sig2& ~corrP),20,regressors(2).color,'^','filled');
    scatter(betas(1,sig1&sig2& ~corrP),betas(2,sig1&sig2& ~corrP),15,'k','filled');
    tr=[' \color[rgb]{0 0 0} ' num2str(sum(sig1&sig2&~corrP)) ' \color[rgb]{' num2str(regressors(1).color) '} ' num2str(sum(sig1&~sig2&~corrP)) ...
        ' \color[rgb]{' num2str(regressors(2).color) '} ' num2str(sum(~sig1&sig2&~corrP)) ' \color[rgb]{0.7 0.7 0.7} ' num2str(sum(~sig1&~sig2&~corrP))...
        ' \color[rgb]{0.4 0.4 0.4} ' num2str(sum(~sig1&~sig2&corrP)) ' \color[rgb]{0 0 1} ' num2str(sum((sig1|sig2)&corrP))];
    title([keys.PSTH_WINDOWS{w,1} ': ' tr],'interpreter','tex'); %%
    %title (keys.PSTH_keys.PSTH_WINDOWS{w,1})
    axis square
    axis equal
end
plot_title              = ['Beta values summary '];
ph_title_and_save(gcf,  [plot_title ', linear regresion'],plot_title,keys)


%% summary of correlation values
plot_title              = ['effector correlations'];
figure_handle              = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
for w=1:n_windows
    subplot(ceil(sqrt(n_windows)),ceil(sqrt(n_windows)),w)
    
    hold on
    plot([-1 1],[-1 1],'k');
    plot([0 0],[-1 0],'k:');
    plot([-1 0],[0 0],'k:');
    
    for sig=0:3
        units=find([per_window(w).corrCIsig]==sig);
        for u=units
            switch per_window(w).corrCIsig(u)
                case 0
                    col=[0.5 0.5 0.5];
                case 1
                    col=regressors(1).color;
                case 2
                    col=regressors(2).color;
                case 3
                    col=[0 0 0];
            end
            if per_window(w).corrCIsig(u)>0 && per_window(w).corrCIsig(u)~=3
                plot(per_window(w).corrCI1(u,:),[per_window(w).corr2(u) per_window(w).corr2(u)],'Color',col,'linewidth',0.001);
                plot([per_window(w).corr1(u) per_window(w).corr1(u)],per_window(w).corrCI2(u,:),'Color',col,'linewidth',0.001);
            end
            scatter(per_window(w).corr1(u),per_window(w).corr2(u),10,col,'filled');
        end
    end
    n0=sum([per_window(w).corrCIsig]==0);
    n1=sum([per_window(w).corrCIsig]==1);
    n2=sum([per_window(w).corrCIsig]==2);
    n3=sum([per_window(w).corrCIsig]==3);
    text(-0.9,-0.9,['N=' num2str(n0)],'color',[0.5 0.5 0.5]);
    text(-0.9,-0.8,['N=' num2str(n1)],'color',regressors(1).color);
    text(-0.9,-0.7,['N=' num2str(n2)],'color',regressors(2).color);
    text(-0.9,-0.6,['N=' num2str(n3)],'color',[0 0 0]);
    %     tr=[' \color[rgb]{0 0 0} ' num2str(sum(sig1&sig2)) ' \color[rgb]{' num2str(regressors(1).color) '} ' num2str(sum(sig1&~sig2)) ' \color[rgb]{' num2str(regressors(2).color) '} ' num2str(sum(~sig1&sig2)) ' \color[rgb]{0.5 0.5 0.5} ' num2str(sum(~sig1&~sig2))];
    %     title([keys.PSTH_WINDOWS{w,1} ': ' tr],'interpreter','tex'); %%
    title(keys.WINDOWS_PER_TYPE{typ}(w,1))
    xlim([-1.001 1.001]);
    ylim([-1.001 1.001]);
    axis square
end
ph_title_and_save(figure_handle,plot_title,plot_title,keys)

%% summary of correlation values (saccades versus sum)
plot_title              = ['effector correlations sum versus saccades only'];
figure_handle              = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
for w=1:n_windows
    subplot(ceil(sqrt(n_windows)),ceil(sqrt(n_windows)),w)
    
    hold on
    plot([-1 1],[-1 1],'k');
    plot([0 0],[-1 0],'k:');
    plot([-1 0],[0 0],'k:');
    
    for sig=0:3
        units=find([per_window(w).corrCIsig1]==sig);
        for u=units
            switch per_window(w).corrCIsig1(u)
                case 0
                    col=[0.5 0.5 0.5];
                case 1
                    col=regressors(1).color;
                case 2
                    col=solution.color;
                case 3
                    col=[0 0 0];
            end
            if per_window(w).corrCIsig1(u)>0 && per_window(w).corrCIsig1(u)~=3
                plot(per_window(w).corrCI1(u,:),[per_window(w).corr3(u) per_window(w).corr3(u)],'Color',col,'linewidth',0.001);
                plot([per_window(w).corr1(u) per_window(w).corr1(u)],per_window(w).corrCI3(u,:),'Color',col,'linewidth',0.001);
            end
            scatter(per_window(w).corr1(u),per_window(w).corr3(u),10,col,'filled');
        end
    end
    n0=sum([per_window(w).corrCIsig1]==0);
    n1=sum([per_window(w).corrCIsig1]==1);
    n2=sum([per_window(w).corrCIsig1]==2);
    n3=sum([per_window(w).corrCIsig1]==3);
    text(-0.9,-0.9,['N=' num2str(n0)],'color',[0.5 0.5 0.5]);
    text(-0.9,-0.8,['N=' num2str(n1)],'color',regressors(1).color);
    text(-0.9,-0.7,['N=' num2str(n2)],'color',solution(1).color);
    text(-0.9,-0.6,['N=' num2str(n3)],'color',[0 0 0]);
    %     tr=[' \color[rgb]{0 0 0} ' num2str(sum(sig1&sig2)) ' \color[rgb]{' num2str(regressors(1).color) '} ' num2str(sum(sig1&~sig2)) ' \color[rgb]{' num2str(regressors(2).color) '} ' num2str(sum(~sig1&sig2)) ' \color[rgb]{0.5 0.5 0.5} ' num2str(sum(~sig1&~sig2))];
    %     title([keys.PSTH_WINDOWS{w,1} ': ' tr],'interpreter','tex'); %%
    title(keys.WINDOWS_PER_TYPE{typ}(w,1))
    xlim([-1.001 1.001]);
    ylim([-1.001 1.001]);
    axis square
end
ph_title_and_save(figure_handle,plot_title,plot_title,keys)

%% summary of correlation values (reaches versus sum)
plot_title              = ['effector correlations sum versus reaches only'];
figure_handle              = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
for w=1:n_windows
    subplot(ceil(sqrt(n_windows)),ceil(sqrt(n_windows)),w)
    
    hold on
    plot([-1 1],[-1 1],'k');
    plot([0 0],[-1 0],'k:');
    plot([-1 0],[0 0],'k:');
    
    for sig=0:3
        units=find([per_window(w).corrCIsig2]==sig);
        for u=units
            switch per_window(w).corrCIsig2(u)
                case 0
                    col=[0.5 0.5 0.5];
                case 1
                    col=solution(1).color;
                case 2
                    col=regressors(2).color;
                case 3
                    col=[0 0 0];
            end
            if per_window(w).corrCIsig2(u)>0 && per_window(w).corrCIsig2(u)~=3
                plot(per_window(w).corrCI3(u,:),[per_window(w).corr2(u) per_window(w).corr2(u)],'Color',col,'linewidth',0.001);
                plot([per_window(w).corr3(u) per_window(w).corr3(u)],per_window(w).corrCI2(u,:),'Color',col,'linewidth',0.001);
            end
            scatter(per_window(w).corr3(u),per_window(w).corr2(u),10,col,'filled');
        end
    end
    n0=sum([per_window(w).corrCIsig2]==0);
    n1=sum([per_window(w).corrCIsig2]==1);
    n2=sum([per_window(w).corrCIsig2]==2);
    n3=sum([per_window(w).corrCIsig2]==3);
    text(-0.9,-0.9,['N=' num2str(n0)],'color',[0.5 0.5 0.5]);
    text(-0.9,-0.8,['N=' num2str(n1)],'color',solution(1).color);
    text(-0.9,-0.7,['N=' num2str(n2)],'color',regressors(2).color);
    text(-0.9,-0.6,['N=' num2str(n3)],'color',[0 0 0]);
    %     tr=[' \color[rgb]{0 0 0} ' num2str(sum(sig1&sig2)) ' \color[rgb]{' num2str(regressors(1).color) '} ' num2str(sum(sig1&~sig2)) ' \color[rgb]{' num2str(regressors(2).color) '} ' num2str(sum(~sig1&sig2)) ' \color[rgb]{0.5 0.5 0.5} ' num2str(sum(~sig1&~sig2))];
    %     title([keys.PSTH_WINDOWS{w,1} ': ' tr],'interpreter','tex'); %%
    title(keys.WINDOWS_PER_TYPE{typ}(w,1))
    xlim([-1.001 1.001]);
    ylim([-1.001 1.001]);
    axis square
end
ph_title_and_save(figure_handle,plot_title,plot_title,keys)


%% plot examples
for u=1:numel(per_unit)
    clear to_plot
    figure
    n_windows=numel(per_unit(u).window);
    n_conditions=sum(con4reg(:,1));
    for c=1:n_conditions
        sph(c)=subplot(n_conditions,1,c);
        hold on
        state_shift=0;
        Betas='';
        for w=1:n_windows
            clear to_plot
            Betas=[Betas ' w' num2str(w) ': '];
            t_before_state=keys.PSTH_WINDOWS{w,3};
            t_after_state=keys.PSTH_WINDOWS{w,4};
            bins=t_before_state:keys.PSTH_binwidth:t_after_state;
            bins=bins+state_shift-t_before_state;
            for r=1:size(con4reg,2)-1
                reg_con_tmp=find(con4reg(:,r));
                beta(r)=per_unit(u).window(w).betas(r);
                in=per_unit(u).window(w).intercept;
                Betas=[Betas ', ' num2str(beta(r))];
                to_plot(r,:)=PSTHs_per_unit(u).window(w).condition(reg_con_tmp(c)).average_spike_density; %*beta(r)+in/2;
                plot(bins,to_plot(r,:),'color',regressors(r).color);
            end
            %mean_to_plot=sum(to_plot);
            plot(bins,sum(to_plot),'color','k');
            %             for r=1:size(con4reg,2)-1
            %                 to_plot_temp=(to_plot(r,:)-mean_to_plot)*beta(r)+mean_to_plot+in/numel(regressors); %% splitting intercept to all regressors
            %                 plot(bins,to_plot_temp,'color',regressors(r).color);
            %             end
            %            plot(bins,sum((to_plot-repmat(mean_to_plot,size(to_plot,1),1)).*repmat(beta',1,size(to_plot,2)))+mean_to_plot+in,'color','k');
            r=r+1;
            reg_con_tmp=find(con4reg(:,r));
            plot(bins,PSTHs_per_unit(u).window(w).condition(reg_con_tmp(c)).average_spike_density,'color',solution.color)
            state_shift=state_shift+t_after_state-t_before_state+0.1;
            title(Betas);
        end
        y_lim(c,:)=get(gca,'ylim');
    end
    
    %% subplot appearance, and tuning lines
    %         spf(sph==0)=[];
    sph(sph==0)=[];
    ylimmax=max(max(y_lim));
    ylimmin=min(min(y_lim));
    y_lim=[ylimmin ylimmax];
    %         if ~isempty(keys.RE.y_lim)
    %             y_lim= keys.RE.y_lim;
    %         end
    %
    for spn=1:numel(sph)
        
        %ef=spf(spn);
        %set(0, 'CurrentFigure', PSTH_summary_handle(ef));
        subplot(sph(spn));
        hold on
        
        %% completed? choices? hands?
        tr=[all_trialz.type]==typ & ismember([all_trialz.completed],keys.cal.completed) &...
            ismember([all_trialz.completed],keys.cal.completed) & ismember([all_trialz.reach_hand],u_con.reach_hand) & ismember([all_trialz.choice],u_con.choice);
        ph_PSTH_background(all_trialz(tr),y_lim,y_lim,y_lim,keys,keys.RE.fontsize_factor)
        
    end
    %             plot_title              = [fig_title plot_title_part ', ' type_effector_short ];
    %             ph_title_and_save(PSTH_summary_handle(ef),  [filename plot_title_part ', ' type_effector_short],plot_title,keys)
    %         for eff=u_con.effector
    %             ef=find(u_con.effector==eff);
    %             [~, type_effector_short] = MPA_get_type_effector_name(typ,eff);
    %         end
    
    plot_title              = ['unit ' complete_unit_list{u}];
    ph_title_and_save(gcf,  [plot_title ', linear regresion'],plot_title,keys)
end