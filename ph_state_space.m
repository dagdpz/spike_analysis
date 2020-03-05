function ph_state_space(population,modified_keys)

for fn=fieldnames(modified_keys)'
    keys.(fn{:})=modified_keys.(fn{:});
end
keys.path_to_save=[keys.basepath_to_save keys.project_version filesep 'state_space_analysis' filesep];
if ~exist(keys.path_to_save,'dir')
    mkdir([keys.basepath_to_save keys.project_version ], 'state_space_analysis');
end

legend_labels_hf={'NH IS IN' 'NH IS CH' 'IH IS IN' 'IH IS CH' 'CH IS IN' 'CH IS CH' ...
    'NH VS IN' 'NH VS CH' 'IH VS IN' 'IH VS CH' 'CH VS IN' 'CH VS CH' ...
    'NH CS IN' 'NH CS CH' 'IH CS IN' 'IH CS CH' 'CH CS IN' 'CH CS CH' ...
    'NH IS IN P' 'NH IS CH P' 'IH IS IN P' 'IH IS CH P' 'CH IS IN P' 'CH IS CH P'...
    'NH VS IN P' 'NH VS CH P' 'IH VS IN P' 'IH VS CH P' 'CH VS IN P' 'CH VS CH P'...
    'NH CS IN P' 'NH CS CH P' 'IH CS IN P' 'IH CS CH P' 'CH CS IN P' 'CH CS CH P'};
legend_labels_pref={'NH NP IN' 'NH NP CH' 'IH NP IN' 'IH NP CH' 'CH NP IN' 'CH NP CH' ...
    'NH PF IN' 'NH PF CH' 'IH PF IN' 'IH PF CH' 'CH PF IN' 'CH PF CH' ...
    'NH NP IN P' 'NH NP CH P' 'IH NP IN P' 'IH NP CH P' 'CH NP IN P' 'CH NP CH P'  ...
    'NH PF IN P' 'NH PF CH P' 'IH PF IN P' 'IH PF CH P' 'CH PF IN P' 'CH PF CH P'};

epoch_labels ={1 2 3 4 5 6 8 62};
epoch_labels(2,:) = {'INI' 'Fix_acq' 'Fix_hol' 'Tar_acq' 'Tar_hol' 'Cue' 'Del' 'Reach_ini'};

cols=keys.colors;
%% there are just too many colors once we include vertical targets, so for now we just keep use the same ones again...
keys.line_colors=[[cols.NH_IS_IN;cols.NH_IS_CH;cols.IH_IS_IN;cols.IH_IS_CH;cols.CH_IS_IN;cols.CH_IS_CH;...
    cols.NH_VS_IN;cols.NH_VS_CH;cols.IH_VS_IN;cols.IH_VS_CH;cols.CH_VS_IN;cols.CH_VS_CH;...
    cols.NH_CS_IN;cols.NH_CS_CH;cols.IH_CS_IN;cols.IH_CS_CH;cols.CH_CS_IN;cols.CH_CS_CH;]/255;...
    [cols.NH_IS_IN;cols.NH_IS_CH;cols.IH_IS_IN;cols.IH_IS_CH;cols.CH_IS_IN;cols.CH_IS_CH;...
    cols.NH_VS_IN;cols.NH_VS_CH;cols.IH_VS_IN;cols.IH_VS_CH;cols.CH_VS_IN;cols.CH_VS_CH;...
    cols.NH_CS_IN;cols.NH_CS_CH;cols.IH_CS_IN;cols.IH_CS_CH;cols.CH_CS_IN;cols.CH_CS_CH;]/400]; %%temporary for inactivation /610
keys.pref_colors=[[cols.NH_IS_IN;cols.NH_IS_CH;cols.IH_IS_IN;cols.IH_IS_CH;cols.CH_IS_IN;cols.CH_IS_CH;...
    cols.NH_CS_IN;cols.NH_CS_CH;cols.IH_CS_IN;cols.IH_CS_CH;cols.CH_CS_IN;cols.CH_CS_CH;]/255;...
    [cols.NH_IS_IN;cols.NH_IS_CH;cols.IH_IS_IN;cols.IH_IS_CH;cols.CH_IS_IN;cols.CH_IS_CH;...
    cols.NH_CS_IN;cols.NH_CS_CH;cols.IH_CS_IN;cols.IH_CS_CH;cols.CH_CS_IN;cols.CH_CS_CH;]/610]; %%temporary for inactivation



% n_col=size(keys.line_colors,1);

%% tuning table preparation and grouping
[tuning_per_unit_table]                 = ph_load_extended_tuning_table(keys);
[tuning_per_unit_table, Sel_for_title]  = ph_reduce_tuning_table(tuning_per_unit_table,keys);
if keys.ST.FR_subtract_baseline
    Sel_for_title =[Sel_for_title,{'base';'=';keys.ST.epoch_BL;', '}];
end
idx_group_parameter=DAG_find_column_index(tuning_per_unit_table,keys.ST.group_parameter);
idx_unitID=DAG_find_column_index(tuning_per_unit_table,'unit_ID');
group_values=tuning_per_unit_table(:,idx_group_parameter);
group_values=cellfun(@num2str, group_values, 'UniformOutput', false);
cell_in_any_group=[false; ~ismember(group_values(2:end),keys.ST.group_excluded)];
unique_group_values=unique(group_values(cell_in_any_group));
if isempty(unique_group_values)
    disp('no relevant groups found');
    return;
end
complete_unit_list={population.unit_ID}';

%% convert lateral to lesional in case of left hemisphere

%invert labels and colors for left hemisphere in order to have
%contralateral to inactivation
if any(strfind(Sel_for_title{3,1},'_L'))
    
    legend_labels_hf={'NH CS IN' 'NH CS CH' 'CH CS IN' 'CH CS CH' 'IH CS IN' 'IH CS CH' ...
        'NH VS IN' 'NH VS CH' 'CH VS IN' 'CH VS CH' 'IH VS IN' 'IH VS CH' ...
        'NH IS IN' 'NH IS CH' 'CH IS IN' 'CH IS CH' 'IH IS IN' 'IH IS CH' ...
        'NH CS IN P' 'NH CS CH P' 'CH CS IN P' 'CH CS CH P' 'IH CS IN P' 'IH CS CH P'...
        'NH VS IN P' 'NH VS CH P' 'CH VS IN P' 'CH VS CH P' 'IH VS IN P' 'IH VS CH P'...
        'NH IS IN P' 'NH IS CH P' 'CH IS IN P' 'CH IS CH P' 'IH IS IN P' 'IH IS CH P'};
    
    keys.line_colors=[[cols.NH_CS_IN;cols.NH_CS_CH;cols.CH_CS_IN;cols.CH_CS_CH;cols.IH_CS_IN;cols.IH_CS_CH;...
        cols.NH_VS_IN;cols.NH_VS_CH;cols.CH_VS_IN;cols.CH_VS_CH;cols.IH_VS_IN;cols.IH_VS_CH;...
        cols.NH_IS_IN;cols.NH_IS_CH;cols.CH_IS_IN;cols.CH_IS_CH;cols.IH_IS_IN;cols.IH_IS_CH;]/255;...
        [cols.NH_CS_IN;cols.NH_CS_CH;cols.CH_CS_IN;cols.CH_CS_CH;cols.IH_CS_IN;cols.IH_CS_CH;...
        cols.NH_VS_IN;cols.NH_VS_CH;cols.CH_VS_IN;cols.CH_VS_CH;cols.IH_VS_IN;cols.IH_VS_CH;...
        cols.NH_IS_IN;cols.NH_IS_CH;cols.CH_IS_IN;cols.CH_IS_CH;cols.IH_IS_IN;cols.IH_IS_CH;]/400]; %%temporary for inactivation /610
end

%% markers for different monkeys and colors for different grid holes
idx_grid_x=DAG_find_column_index(tuning_per_unit_table,'grid_x');
idx_grid_y=DAG_find_column_index(tuning_per_unit_table,'grid_y');
idx_grid_z=DAG_find_column_index(tuning_per_unit_table,'electrode_depth');

monkey_markers={};
for m=1:numel(keys.batching.monkeys)
    if isfield(keys,keys.batching.monkeys{m})
        monkey_markers=[monkey_markers keys.(keys.batching.monkeys{m}).marker];
    else
        monkey_markers=[monkey_markers 'o'];
    end
end
gridhole_colors={[70 0 0;255 0 0],[0 70 0;0 255 0],[0 70 70; 0 255 255],[70 70 0; 255 255 0]};


monkey_per_cell=cellfun(@(x) x(1:3),tuning_per_unit_table(:,idx_unitID),'uniformoutput',false);
[unique_monkeys]=unique(monkey_per_cell(2:end));
for m=1:numel(unique_monkeys)
    idx_m=cell_in_any_group&ismember(monkey_per_cell,unique_monkeys(m));
    [unique_grid2_values{m}]=unique(cell2mat([tuning_per_unit_table(idx_m,idx_grid_x), tuning_per_unit_table(idx_m,idx_grid_y)]),'rows');
    
    if ~isempty(unique_grid2_values{m})
        unique_gridx_values{m}=[unique_grid2_values{m}(:,1)];
        unique_gridy_values{m}=[unique_grid2_values{m}(:,2)];
        
        %electrode_depth_range{m}=[min([tuning_per_unit_table{idx_m,idx_grid_z}]) max([tuning_per_unit_table{idx_m,idx_grid_z}])];
        for gh=1:numel(unique_gridx_values{m})
            idx_g=[false; ismember([tuning_per_unit_table{2:end,idx_grid_x}],unique_gridx_values{m}(gh))' & ismember([tuning_per_unit_table{2:end,idx_grid_y}],unique_gridy_values{m}(gh) )'];
            electrode_depth_range{m}{gh}=[min([tuning_per_unit_table{idx_m & idx_g,idx_grid_z}]) max([tuning_per_unit_table{idx_m & idx_g,idx_grid_z}])];
        end
    end
end

%% gaussian fit settings
fitsettings.sd_max_x=12;
fitsettings.sd_x_min_ratio=0.125;
fitsettings.sd_max_y=fitsettings.sd_max_x;
fitsettings.sd_xy_min_ratio=0.25;
fitsettings.sd_xy_max_ratio=1;
fitsettings.sd_y_min_ratio=fitsettings.sd_x_min_ratio;
fitsettings.xout=[-30:30]; % range for gaussian fit
fitsettings.yout=[-15:15];
fitsettings.range_factor=1;

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
per_trial.perturbation(ismember(per_trial.perturbation, keys.cal.perturbation_groups{1}))=0;
per_trial.perturbation(ismember(per_trial.perturbation, keys.cal.perturbation_groups{2}))=1;

%hemifields=unique([whatisthis.trial.hemifield]);
u_hemifields=[-1,0,1]; % why does this have to be hardcoded? ---> Because case not defined yet, case defines positions !!

u_types     =unique(per_trial.types);
u_effectors =unique(per_trial.effectors);
u_hands     =unique(per_trial.hands);
u_choice    =unique(per_trial.choice);
u_perturbation    =unique(per_trial.perturbation);
u_perturbation=u_perturbation(~isnan(u_perturbation));

if any(u_perturbation==1) % temporary, better solution
    u_hemifields=[-1,1]; % why does this have to be hardcoded? And what do we do with vertical targets?
end
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
u_types     =unique(type_effectors(:,1));
u_effectors =unique(type_effectors(:,2));

condition_matrix    = combvec(u_hands,u_choice, u_perturbation,u_hemifields)';
conditions_out      = combvec(u_effectors,u_hands,u_choice, u_perturbation)';


conditions_hf       = combvec(u_hemifields,conditions_out')';
conditions_hf_complete       = combvec(u_hemifields,conditions_out')';
conditions_pref     = combvec([1 0],conditions_out')';

[~,~,columns_hf] = unique(conditions_hf(:,[1,3,4]),'rows');
[~,~,columns_pref] = unique(conditions_pref(:,[1,3,4]),'rows');

if any(u_perturbation==1) % temporary, better solution
    if strcmp(keys.arrangement,'hands_inactivation_in_ch')
        [~,~,columns_hf] = unique(conditions_hf(:,[1,3]),'rows');
    end
    condition_matrix(condition_matrix(:,1)==1 & condition_matrix(:,2)==1 & condition_matrix(:,4)==1,:)=[];
    condition_matrix(condition_matrix(:,1)==2 & condition_matrix(:,2)==1 & condition_matrix(:,4)==-1,:)=[];
    conditions_hf(conditions_hf(:,1)==1 & conditions_hf(:,3)==1 & conditions_hf(:,4)==1,:)=[];
    conditions_hf(conditions_hf(:,1)==-1 & conditions_hf(:,3)==2 & conditions_hf(:,4)==1,:)=[];
else
    
    columns_hf          = ones(size(conditions_hf,1),1);
    columns_pref        = ones(size(conditions_pref,1),1);
end




%% Convert to ipsi/contra, Baseline subtraction, normalization, re-ordering, and gaussian RFs
% for a=1:numel(keys.position_and_plotting_arrangements)
%     keys.arrangement=keys.position_and_plotting_arrangements{a}; % important for ph_arrange_positions_and_plots
a=1;
tr_con=ismember([all_trialz.completed],keys.cal.completed);

[whatisthis]=ph_arrange_positions_and_plots(all_trialz(tr_con),keys);
positions=unique(vertcat(whatisthis.trial.position),'rows');
fixations_temp=unique(vertcat(whatisthis.trial.fixation),'rows');
fix_temp_idx=true(size(fixations_temp,1),1);
for x=1:size(fixations_temp,1)-1
    if any(all(abs(bsxfun(@minus,fixations_temp(x+1:end,:),fixations_temp(x,:)))<4,2)) %% precision....
        fix_temp_idx(x)=false;
    end
end
fixations=fixations_temp(fix_temp_idx,:);
%        per_arrangement(a).positions=positions;
%     per_arrangement(a).fixations=fixations;
clear whatisthis

%% type? (+effector?)

for t=1:size(type_effectors,1) %typ=unique(per_trial.types)
    typ=type_effectors(t,1);
    eff=type_effectors(t,2);
    
    conditions_eff=conditions_out(conditions_out(:,1)==eff,2:end);
    tya=find(u_types==typ)+a*(numel(u_types))-1;
    idx_existing=DAG_find_column_index(tuning_per_unit_table,['existing_' 'in_AH_' type_effector_short{t} '_' keys.arrangement(1:3)]);
    tya_existing{tya}= [true; ismember(vertcat(tuning_per_unit_table{2:end,idx_existing}),true)];
    
    keys=ph_get_epoch_keys(keys,typ,eff,sum(type_effectors(:,1)==typ)>1);
    if ~isfield(keys.ST,'epoch_GB')
        keys.ST.epoch_GB=[keys.ANOVAS_PER_TYPE(typ).epoch{strcmp(keys.ANOVAS_PER_TYPE(typ).epoch(:,2),keys.ST.epoch_RF),1}]; %% gaussian baseline depends on the epoch for resposne field
    end
    if isempty(keys.ST.epoch_GB); keys.ST.epoch_GB={''}; end
    normalization_epoch     =find(ismember(keys.EPOCHS(:,1),keys.ST.epoch_for_normalization));
    receptive_field_epoch   =find(ismember(keys.EPOCHS(:,1),keys.ST.epoch_RF));
    baseline_epoch          =find(ismember(keys.EPOCHS(:,1),keys.ST.epoch_BL));
    preference_epoch        =find(ismember(keys.EPOCHS(:,1),keys.ST.epoch_PF));
    gaussian_bl_epoch       =find(ismember(keys.EPOCHS(:,1),keys.ST.epoch_GB));
    
    pref_valid(tya,:)=zeros(size(complete_unit_list,1),1);
    for g=1:numel(unique_group_values)
        unitidx=ismember(complete_unit_list,tuning_per_unit_table(ismember(group_values,unique_group_values(g))&tya_existing{tya},idx_unitID));
        units=find(all(unitidx,2))';
        for u=units
            tr_con=ismember([population(u).trial.completed],keys.cal.completed) & [population(u).trial.accepted];
            [pop]=ph_arrange_positions_and_plots(population(u).trial(tr_con),keys);
            
            [pop.trial(ismember([pop.trial.perturbation], keys.cal.perturbation_groups{1})).perturbation]=deal(0);
            [pop.trial(ismember([pop.trial.perturbation], keys.cal.perturbation_groups{2})).perturbation]=deal(1);
            pop=ph_LR_to_CI(pop,population(u).target);
            
            %% normalization factor
            for c=1:size(condition_matrix,1)
                clear trpar
                switch keys.ST.normalization
                    case 'by_condition'
                        condition_parameters_tmp=[condition_parameters {'hemifield'}];
                        for par=1:numel(condition_parameters_tmp)
                            fn=condition_parameters_tmp{par};
                            trpar(par,:)=[pop.trial.(fn)]==condition_matrix(c,par);
                        end
                        trpar(end+1,:)=[pop.trial.type]==typ & [pop.trial.effector]==eff;
                    case 'by_perturbation'
                        trpar=[pop.trial.type]==typ & [pop.trial.effector]==eff & [pop.trial.perturbation]==condition_matrix(c,strcmp(condition_parameters,'perturbation'));
                    case 'by_effector'
                        trpar=[pop.trial.type]==typ & [pop.trial.effector]==eff;
                    case 'by_type'
                        trpar=[pop.trial.type]==typ;
                    case 'by_all_trials'
                        trpar=true(size([pop.trial]));
                    case 'none' %% really doesnt need an extra rule actusally?? Not here at least..
                        trpar=[pop.trial.type]==typ;
                    case 'z_score' %% really doesnt need an extra rule actusally?? Not here at least..
                        trpar=[pop.trial.type]==typ & [pop.trial.effector]==eff;
                end
                tr_con=all(trpar,1);
                per_epoch=vertcat(pop.trial(tr_con).epoch);
                gaussian_baseline=vertcat(per_epoch(:,gaussian_bl_epoch).FR);
                if isempty(gaussian_baseline)
                    gaussian_baseline=zeros(size(per_epoch,1),1);
                end
                if keys.ST.FR_subtract_baseline
                    baseline(c)=nanmean(vertcat(per_epoch(:,baseline_epoch).FR));
                else
                    baseline(c)=0;
                end
                
                %what to do if per_epoch is empty or
                if strcmp(keys.ST.normalization,'none') || isempty(normalization_epoch) %normalization_epoch is empty? %% check if this works, for no multiplicative normalization
                    norm_factor(u,c)= 1;
                elseif isempty(per_epoch)
                    norm_factor(u,c)=NaN;
                elseif ~isempty(normalization_epoch)
                    norm_factor(u,c)=nanmean([per_epoch(:,normalization_epoch).FR NaN]);
                end
            end
            
            %% z-scoring
            if strcmp(keys.ST.normalization,'z_score')
                u_blocks=unique([pop.trial.block]);
                AT_tmp={};
                prev_block_end=0;
                trial_onset_times=zeros(size(pop.trial));
                for b=1:numel(u_blocks)
                    b_idx=[pop.trial.block]==u_blocks(b);
                    AT_tmp=[AT_tmp arrayfun(@(x) single(x.arrival_times + x.trial_onset_time + prev_block_end),pop.trial(b_idx),'uniformoutput',false)];
                    last_trial_idx=find(b_idx,1,'last');
                    trial_onset_times(b_idx)=[pop.trial(b_idx).trial_onset_time]+prev_block_end;
                    prev_block_end=prev_block_end+pop.trial(last_trial_idx).trial_onset_time+pop.trial(last_trial_idx).states_onset(pop.trial(last_trial_idx).states==98);
                    
                end
                AT_tmp=unique(vertcat(AT_tmp{:}));
                PSTH_ms =ceil(double(trial_onset_times(1))/keys.PSTH_binwidth)*keys.PSTH_binwidth:0.001:...
                    floor(double(trial_onset_times(end)+pop.trial(end).states_onset(pop.trial(end).states==98))/keys.PSTH_binwidth)*keys.PSTH_binwidth;
                SD_ms= conv(hist(AT_tmp,PSTH_ms),normpdf(-5*keys.gaussian_kernel:0.001:5*keys.gaussian_kernel,0,keys.gaussian_kernel),'same');
                
                t_idx=false(size(PSTH_ms));
                for t_temp=find([pop.trial.type]==typ & [pop.trial.effector]==eff)
                    min_temp=ceil(trial_onset_times(t_temp)/keys.PSTH_binwidth)*keys.PSTH_binwidth;
                    max_temp=floor((trial_onset_times(t_temp)+double(pop.trial(t_temp).states_onset(pop.trial(t_temp).states==50)))/keys.PSTH_binwidth)*keys.PSTH_binwidth;
                    t_idx(PSTH_ms+0.0002>=min_temp & PSTH_ms+0.0002<max_temp)=true;
                end
                SD=nanmean(reshape(SD_ms(t_idx),keys.PSTH_binwidth/0.001,sum(t_idx)/keys.PSTH_binwidth*0.001),1);
                norm_factor(u,:)=deal(std(SD));
                baseline(:)=deal(nanmean(SD));
            end
            
            norm_factor(u,:)=deal(max([norm_factor(u,:) 0])); % now we always normalize to maximum condition, 0 makes sure some value is there..
            norm_factor(norm_factor==0)=1;
            
            
            
            
            
            %% average PSTH per unit
            data_to_PCA = struct();
            for w=1:size(conditions_eff,1)
                clear trpar
                %condition trial index
                % condition(c).parameters=condition_matrix(c,:);
                
                %c=find(all(bsxfun(@minus,conditions_out,[eff,condition_matrix(w,1),condition_matrix(w,2),condition_matrix(w,3)])==0,2));
                c=w;
                
                for par=1:numel(condition_parameters)
                    fn=condition_parameters{par};
                    trpar(par,:)=[pop.trial.(fn)]==conditions_eff(w,par); %condition_matrix_wohf?
                    %trpar(par,:)=arrayfun(@(x) any([x.(fn)]==condition_matrix(c,par)), pop.trial);
                end
                trpar(end+1,:)=[pop.trial.type]==typ & [pop.trial.effector]==eff;
                tr_con=all(trpar,1);
                per_epoch=vertcat(pop.trial(tr_con).epoch);
                
                for wn=1:size(keys.PSTH_WINDOWS,1)
                    for f=1:numel(u_hemifields) %hemifield
                        trtemp=[pop.trial.hemifield]==u_hemifields(f);
                        %n=(f-1)*size(condition_matrix,1)/numel(u_hemifields);
                        %n=(f-1)*size(condition_matrix,1)/numel(u_hemifields);
                        n=max([1,find(ismember(condition_matrix(:,1:3),conditions_eff(w,:),'rows') & condition_matrix(:,end)==u_hemifields(f))]);
                        average_FR_per_unit.unit(u).condition(c).hemifield(f).windows(wn).average_spike_density = ph_spike_density(pop.trial(tr_con & trtemp),wn,keys,baseline(c),norm_factor(u,n));
                        %                         condition(tya,c).per_hemifield(f).window(wn).unit(u).average_spike_density=...
                        %                             ph_spike_density(pop.trial(tr_con & trtemp),wn,keys,baseline(c),norm_factor(u,n));
                        %                         condition(tya,c).per_hemifield(f).effector=eff;
                        
                    end
                    
                    %                     for p=1:size(positions,1)
                    %                         f=sign(positions(p,1)); %what happens if position is neither left nor right??
                    %                         %n=(find(u_hemifields==f)-1)*size(condition_matrix,1)/numel(u_hemifields);
                    %                         n=max([1,find(ismember(condition_matrix(:,1:3),conditions_eff(w,:),'rows') & condition_matrix(:,end)==f)]);
                    %
                    %
                    %                         %     baseline_epoch          =find(ismember(keys.EPOCHS(:,1),keys.ST.epoch_BL));
                    %                         %     preference_epoch        =find(ismember(keys.EPOCHS(:,1),keys.ST.epoch_PF));
                    %
                    %                         trpos=all(abs(bsxfun(@minus,vertcat(pop.trial.position),positions(p,:)))<1.5,2)';
                    %
                    %
                    %                         condition(tya,c).per_position(p).window(wn).unit(u).average_spike_density=...
                    %                             ph_spike_density(pop.trial(tr_con & trpos),wn,keys,baseline(c),norm_factor(u,n));
                    %                         condition(tya,c).per_position(p).position=positions(p,:);
                    %                         condition(tya,c).per_position(p).fixation=1;
                    %                         condition(tya,c).per_position(p).effector=eff;
                    %
                    %
                    %
                    %
                    %
                    %
                    %                     end
                end
                
            end
        end
    end
end
%% prepare data for PCA


%first keep units indexes (which unit go through the PCA and delete empty
%units
for idx = 1:size(average_FR_per_unit.unit,2)
    empty_units(idx) = ~isempty(average_FR_per_unit.unit(idx).condition);
end
average_FR_per_unit.unit = average_FR_per_unit.unit(logical(empty_units)); %remove empty units from the structure

%calculate windows size and condition size (to know what is what in the PCA
%input)
condition_index =0;
for ws = 1:size(keys.PSTH_WINDOWS,1)
    windows_index(1,ws) = size(average_FR_per_unit.unit(1).condition(1).hemifield(1).windows(ws).average_spike_density,2);
    condition_index = condition_index + windows_index(1,ws);
end

% check now that there are value for every conditions/hemifields/windows
% and remove units if not the case
units_to_remove = ones(size(average_FR_per_unit.unit,2),1);
for un = 1:size(average_FR_per_unit.unit,2)
    for cn = 1:size(average_FR_per_unit.unit(un).condition,2)
        for hm = 1:size(average_FR_per_unit.unit(un).condition(cn).hemifield,2)
            for wn = size(keys.PSTH_WINDOWS,1)
                empty_fields.unit(un).condition(cn).hemifield(hm).windows(wn) = ~isempty(average_FR_per_unit.unit(un).condition(cn).hemifield(hm).windows(wn).average_spike_density);
                if units_to_remove(un) == 0;
                    continue
                elseif ~any(empty_fields.unit(un).condition(cn).hemifield(hm).windows)
                    units_to_remove(un) = 0;
                else
                    units_to_remove(un) = 1;
                end
            end
        end
    end
end
average_FR_per_unit.unit = average_FR_per_unit.unit(logical(units_to_remove)); %remove empty units from the structure








%need to concatinate per unit for PCA (here windows first)
%[1st window 2nd window]
for un = 1:size(average_FR_per_unit.unit,2)
    for cn = 1:size(average_FR_per_unit.unit(un).condition,2)
        for hm = 1:size(average_FR_per_unit.unit(un).condition(cn).hemifield,2)
            concat_windows.unit(un).condition(cn).hemifield(hm).average_spike_density = horzcat(average_FR_per_unit.unit(un).condition(cn).hemifield(hm).windows(1).average_spike_density,...
                average_FR_per_unit.unit(un).condition(cn).hemifield(hm).windows(2).average_spike_density);
        end
    end
end




%need to concatinate per unit for PCA (here hand\space condition)
%[XHIS XHCS] (not flexible, hemifield hardcoded)
for un = 1:size(average_FR_per_unit.unit,2)
    for cn = 1:size(average_FR_per_unit.unit(un).condition,2)
        concat_hand_space.unit(un).condition(cn).average_spike_density = horzcat(concat_windows.unit(un).condition(cn).hemifield(1).average_spike_density,...
            concat_windows.unit(un).condition(cn).hemifield(2).average_spike_density);
    end
end

%need to concatinate per unit for PCA (last step, all conditions)
%[IHISCT IHCSCT CHISCT CHCSCT IHISPT IHCSPT CHISPT CHCSPT]

data_to_PCA=[];
data_to_PCA_CT=[];
data_to_PCA_PT=[];

for un = 1:size(average_FR_per_unit.unit,2)
    if keys.ST.combine_exp_conditions %calculate PCA based on all conditions (CT and PT) or calculate CT and PT separately
        data_to_PCA(un,:) = horzcat(concat_hand_space.unit(un).condition(1).average_spike_density,concat_hand_space.unit(un).condition(2).average_spike_density,...
            concat_hand_space.unit(un).condition(3).average_spike_density,concat_hand_space.unit(un).condition(4).average_spike_density);
    else
        data_to_PCA_CT(un,:) = horzcat(concat_hand_space.unit(un).condition(1).average_spike_density(1:condition_index),concat_hand_space.unit(un).condition(2).average_spike_density(1:condition_index),...
            concat_hand_space.unit(un).condition(3).average_spike_density(1:condition_index),concat_hand_space.unit(un).condition(4).average_spike_density(1:condition_index));
        data_to_PCA_PT(un,:) = horzcat(concat_hand_space.unit(un).condition(1).average_spike_density(condition_index+1:2*condition_index),concat_hand_space.unit(un).condition(2).average_spike_density(condition_index+1:2*condition_index),...
            concat_hand_space.unit(un).condition(3).average_spike_density(condition_index+1:2*condition_index),concat_hand_space.unit(un).condition(4).average_spike_density(condition_index+1:2*condition_index));
    end
end

%% prepare data for DataHigh
% data_to_dataHigh = struct();
% conditions_plot=conditions_hf_complete;
% it_cond = 0;
% for cond = 1:size(conditions_plot,1)
%  col=(conditions_plot(cond,1)+1)*6 + (conditions_plot(cond,3))*2 + (conditions_plot(cond,4)+1) + (conditions_plot(cond,5)*18);   %ritrieve color informatiom
%  data_to_dataHigh(cond).data = data_to_PCA(:,it_cond*condition_index+1:cond*condition_index);
%  data_to_dataHigh(cond).traj = 'traj';
%  data_to_dataHigh(cond).epochStarts = 1;
%  data_to_dataHigh(cond).epochColors = keys.line_colors(col,:);
%  data_to_dataHigh(cond).condition = char(legend_labels_hf(col));
%  it_cond = it_cond + 1;
% % plot3(coeff(it_cond*condition_index+1:cond*condition_index,1),coeff(it_cond*condition_index+1:cond*condition_index,2),coeff(it_cond*condition_index+1:cond*condition_index,3),'color',keys.line_colors(col,:),'LineWidth',1.5)
%  end
%
% % DataHigh(data_to_dataHigh,'DimReduce')
% file_to_save = [char(Sel_for_title(3,:)) '_to_dataHigh_' keys.ST.normalization];
% save([keys.drive filesep keys.basepath_to_save filesep keys.project_version filesep 'state_space_analysis' filesep file_to_save ],'data_to_dataHigh')


%%

%% run the PCA
if keys.ST.combine_exp_conditions
    [coeff,score,latent,~,explained] = pca(data_to_PCA);
    PCA_out =  pca(data_to_PCA);
    file_to_save = [char(Sel_for_title(3,:)) '_PCA_out_' keys.ST.normalization];
    save([keys.drive filesep keys.basepath_to_save filesep keys.project_version filesep 'state_space_analysis' filesep file_to_save ],'PCA_out')
    %calculate euclidian distances between hand/space conditions and experimental
    %conditions (HARDCODED!!!)
    euclidian_distances = struct();
    it_cond = 0;
    conditions_plot=conditions_hf_complete;
    for cond = 1:size(conditions_plot,1)
        for time_point = 1:condition_index %[IHISCT IHCSCT CHISCT CHCSCT IHISPT IHCSPT CHISPT CHCSPT]
            %control
            euclidian_distances.Control.IH(time_point) = calculate_3D_euc_dist(coeff(time_point,1:3),coeff(condition_index + time_point,1:3));
            euclidian_distances.Control.CH(time_point) = calculate_3D_euc_dist(coeff(2*condition_index + time_point,1:3),coeff(3*condition_index + time_point,1:3));
            euclidian_distances.Control.IS(time_point)  = calculate_3D_euc_dist(coeff(time_point,1:3),coeff(2*condition_index + time_point,1:3));
            euclidian_distances.Control.CS(time_point)  = calculate_3D_euc_dist(coeff(condition_index + time_point,1:3),coeff(3*condition_index + time_point,1:3));
            %perturbation
            euclidian_distances.Perturbation.IH(time_point) = calculate_3D_euc_dist(coeff(4*condition_index + time_point,1:3),coeff(5*condition_index + time_point,1:3));
            euclidian_distances.Perturbation.CH(time_point) = calculate_3D_euc_dist(coeff(6*condition_index + time_point,1:3),coeff(7*condition_index + time_point,1:3));
            euclidian_distances.Perturbation.IS(time_point)  = calculate_3D_euc_dist(coeff(4*condition_index + time_point,1:3),coeff(6*condition_index + time_point,1:3));
            euclidian_distances.Perturbation.CS(time_point)  = calculate_3D_euc_dist(coeff(5*condition_index + time_point,1:3),coeff(7*condition_index + time_point,1:3));
            %perturbation - control for same condition
            euclidian_distances.Differences.IHIS(time_point) = calculate_3D_euc_dist(coeff(time_point,1:3),coeff(4*condition_index + time_point,1:3));
            euclidian_distances.Differences.IHCS(time_point) = calculate_3D_euc_dist(coeff(condition_index + time_point,1:3),coeff(5*condition_index + time_point,1:3));
            euclidian_distances.Differences.CHIS(time_point) = calculate_3D_euc_dist(coeff(2*condition_index + time_point,1:3),coeff(6*condition_index + time_point,1:3));
            euclidian_distances.Differences.CHCS(time_point) = calculate_3D_euc_dist(coeff(3*condition_index + time_point,1:3),coeff(7*condition_index + time_point,1:3));
        end
    end
    
    
else
    
    PCA_out_CT =  pca(data_to_PCA_CT);
    PCA_out_PT =  pca(data_to_PCA_PT);
    file_to_save_CT = [char(Sel_for_title(3,:)) '_PCA_out_CT' keys.ST.normalization];
    save([keys.drive filesep keys.basepath_to_save filesep keys.project_version filesep 'state_space_analysis' filesep file_to_save_CT ],'PCA_out_CT')
    file_to_save_PT = [char(Sel_for_title(3,:)) '_PCA_out_PT' keys.ST.normalization];
    save([keys.drive filesep keys.basepath_to_save filesep keys.project_version filesep 'state_space_analysis' filesep file_to_save_CT ],'PCA_out_PT')
end


%% main plot
%general plot settings
%plot all trajectories on same plot
if keys.ST.combine_exp_conditions
    plot_title_part = ' neural_state_space';
    fig_title=sprintf('Selection: %s %s %s hnd %s ch %s %s grouped by %s normalized %s in %s',...
        [Sel_for_title{:}],keys.monkey,[keys.conditions_to_plot{:}],mat2str(u_hands),mat2str(double(u_choice)),keys.arrangement,keys.ST.group_parameter, keys.ST.normalization,keys.ST.epoch_for_normalization);
    filename=sprintf('%s %s %s %s %s hnd %s ch %s N_%s %s',...
        [Sel_for_title{:}],keys.monkey,keys.ST.group_parameter,[keys.conditions_to_plot{:}],keys.arrangement(1:3),mat2str(u_hands),mat2str(double(u_choice)),keys.ST.normalization,keys.ST.epoch_for_normalization);
    
    keys=ph_get_epoch_keys(keys,typ,eff,sum(type_effectors(:,1)==typ)>1);
    [~, type_effector_short] = get_type_effector_name(typ,eff);
    plot_title              = [fig_title type_effector_short plot_title_part];
    state_space_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
    conditions_plot=conditions_hf_complete;
    
    it_cond = 0;
    for cond = 1:size(conditions_plot,1)
        col=(conditions_plot(cond,1)+1)*6 + (conditions_plot(cond,3))*2 + (conditions_plot(cond,4)+1) + (conditions_plot(cond,5)*18);   %ritrieve color informatiom
        plot3(smooth(coeff(it_cond*condition_index+1:cond*condition_index,1)',0.5),smooth(coeff(it_cond*condition_index+1:cond*condition_index,2)',0.5),smooth(coeff(it_cond*condition_index+1:cond*condition_index,3)',0.5),'color',keys.line_colors(col,:),'LineWidth',1.5)
        %         plot3(smooth(coeff(it_cond*condition_index+1:cond*condition_index,1)'),smooth(coeff(it_cond*condition_index+1:cond*condition_index,2)'),smooth(coeff(it_cond*condition_index+1:cond*condition_index,3)'),'color',keys.line_colors(col,:),'LineWidth',1.5)
        %plot3(coeff(it_cond*condition_index+1:cond*condition_index,1)',coeff(it_cond*condition_index+1:cond*condition_index,2)',coeff(it_cond*condition_index+1:cond*condition_index,3)','color',keys.line_colors(col,:),'LineWidth',1.5)
        hold on
        % scatter3(coeff(it_cond*condition_index+1,1),coeff(it_cond*condition_index+1,2),coeff(it_cond*condition_index+1,3),'filled','k')
        % scatter3(coeff(cond*windows_index(1)+1,1),coeff(cond*windows_index(1)+1,2),coeff(cond*windows_index(1)+1,3),'filled','r')
        % scatter3(coeff(cond*condition_index,1),coeff(cond*condition_index,2),coeff(cond*condition_index,3),'filled','b')
        it_cond = it_cond +1;
        legend_labels{cond} = char(legend_labels_hf(col));
    end
    legend(legend_labels)
    %estetique
    xlabel('PC 1')
    ylabel('PC 2')
    zlabel('PC 3')
    set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],'XTick',[],'YTick',[],'ZTick',[])
    title_and_save(state_space_summary_handle,  [filename plot_title_part ', ' type_effector_short],plot_title,keys)
    %% plot latent variable separately
    plot_title_part = ' latent_variables';
    fig_title=sprintf('Selection: %s %s %s hnd %s ch %s %s grouped by %s normalized %s in %s',...
        [Sel_for_title{:}],keys.monkey,[keys.conditions_to_plot{:}],mat2str(u_hands),mat2str(double(u_choice)),keys.arrangement,keys.ST.group_parameter, keys.ST.normalization,keys.ST.epoch_for_normalization);
    filename=sprintf('%s %s %s %s %s hnd %s ch %s N_%s %s',...
        [Sel_for_title{:}],keys.monkey,keys.ST.group_parameter,[keys.conditions_to_plot{:}],keys.arrangement(1:3),mat2str(u_hands),mat2str(double(u_choice)),keys.ST.normalization,keys.ST.epoch_for_normalization);
    
    keys=ph_get_epoch_keys(keys,typ,eff,sum(type_effectors(:,1)==typ)>1);
    [~, type_effector_short] = get_type_effector_name(typ,eff);
    plot_title              = [fig_title type_effector_short plot_title_part];
    state_space_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
    conditions_plot=conditions_hf_complete;
    
    for n_latent = 1:3
        it_cond = 0;
        subplot(3,1,n_latent)
        for cond = 1:size(conditions_plot,1)
            col=(conditions_plot(cond,1)+1)*6 + (conditions_plot(cond,3))*2 + (conditions_plot(cond,4)+1) + (conditions_plot(cond,5)*18);   %ritrieve color informatiom
            plot(smooth(coeff(it_cond*condition_index+1:cond*condition_index,n_latent)',0.8),'color',keys.line_colors(col,:),'LineWidth',1.5)
            %              plot(coeff(it_cond*condition_index+1:cond*condition_index,n_latent)','color',keys.line_colors(col,:),'LineWidth',1.5)
            hold on
            it_cond = it_cond +1;
            legend_labels{cond} = char(legend_labels_hf(col));
            legend(legend_labels);
            title(['latent variable ' num2str(n_latent) ' (' num2str(explained(n_latent,1)) '%)']);
            set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],'XTick',[],'YTick',[],'ZTick',[])
        end
    end
    
    
    title_and_save(state_space_summary_handle,  [filename plot_title_part ', ' type_effector_short],plot_title,keys)
    %% plot variance explained
    plot_title_part = ' variance_explained';
    fig_title=sprintf('Selection: %s %s %s hnd %s ch %s %s grouped by %s normalized %s in %s',...
        [Sel_for_title{:}],keys.monkey,[keys.conditions_to_plot{:}],mat2str(u_hands),mat2str(double(u_choice)),keys.arrangement,keys.ST.group_parameter, keys.ST.normalization,keys.ST.epoch_for_normalization);
    filename=sprintf('%s %s %s %s %s hnd %s ch %s N_%s %s',...
        [Sel_for_title{:}],keys.monkey,keys.ST.group_parameter,[keys.conditions_to_plot{:}],keys.arrangement(1:3),mat2str(u_hands),mat2str(double(u_choice)),keys.ST.normalization,keys.ST.epoch_for_normalization);
    
    keys=ph_get_epoch_keys(keys,typ,eff,sum(type_effectors(:,1)==typ)>1);
    [~, type_effector_short] = get_type_effector_name(typ,eff);
    plot_title              = [fig_title type_effector_short plot_title_part];
    state_space_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
    
    
    for n_exp = 1:length(explained)
        cumulative_var_explained(n_exp) = sum(explained(1:n_exp));
    end
    
    plot([1:length(explained)],cumulative_var_explained,'LineWidth',1.5)
    ylabel('variance explained (%)')
    xlabel('principal components')
    
    title_and_save(state_space_summary_handle,  [filename plot_title_part ', ' type_effector_short],plot_title,keys)
    %% plot euclidian distances
    plot_title_part = ' euclidian_distances';
    fig_title=sprintf('Selection: %s %s %s hnd %s ch %s %s grouped by %s normalized %s in %s',...
        [Sel_for_title{:}],keys.monkey,[keys.conditions_to_plot{:}],mat2str(u_hands),mat2str(double(u_choice)),keys.arrangement,keys.ST.group_parameter, keys.ST.normalization,keys.ST.epoch_for_normalization);
    filename=sprintf('%s %s %s %s %s hnd %s ch %s N_%s %s',...
        [Sel_for_title{:}],keys.monkey,keys.ST.group_parameter,[keys.conditions_to_plot{:}],keys.arrangement(1:3),mat2str(u_hands),mat2str(double(u_choice)),keys.ST.normalization,keys.ST.epoch_for_normalization);
    keys=ph_get_epoch_keys(keys,typ,eff,sum(type_effectors(:,1)==typ)>1);
    [~, type_effector_short] = get_type_effector_name(typ,eff);
    plot_title              = [fig_title type_effector_short plot_title_part];
    state_space_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
    
    conditions_plot = [0 4 1 0 0;0 4 2 0 0;-1 4 0 0 0;1 4 0 0 0;0 4 1 0 1;0 4 2 0 1;-1 4 0 0 1;1 4 0 0 1;conditions_hf_complete(5:8,:)]; %hardcoded, such as the pca matrix is
    exp_cond_names = fieldnames(euclidian_distances);
    it_cond = 1;
    eff = keys.cal.types;
    
    for mw = 1:size(keys.WINDOWS_PER_TYPE{1,eff},1)
        epoch_line(mw,1) = modified_keys.WINDOWS_PER_TYPE{1,eff}(mw,2);
        epoch_line(mw,2) = modified_keys.WINDOWS_PER_TYPE{1,eff}(mw,3);
        epoch_index = find([epoch_labels{1,:}] == epoch_line{mw,1});
        epoch_line(mw,3) = epoch_labels(2,epoch_index);
    end
    
    for nstr = 1:size(exp_cond_names,1)
        subplot(size(exp_cond_names,1),1,nstr)
        hold on
        hnd_spa_names = fieldnames(euclidian_distances.(exp_cond_names{nstr}));
        for ncond = 1:numel(hnd_spa_names)
            col=(conditions_plot(it_cond,1)+1)*6 + (conditions_plot(it_cond,3))*2 + (conditions_plot(it_cond,4)+1) + (conditions_plot(it_cond,5)*18);   %ritrieve color informatiom
            plot(smooth(euclidian_distances.(exp_cond_names{nstr}).(hnd_spa_names{ncond}),0.8),'color',keys.line_colors(col,:),'LineWidth',1.5);
%             plot(euclidian_distances.(exp_cond_names{nstr}).(hnd_spa_names{ncond}),'color',keys.line_colors(col,:),'LineWidth',1.5);
            it_cond = it_cond + 1;
            legend_labels{ncond} = char(legend_labels_hf(col));
        end
        
        legend(legend_labels)
        allYLim(nstr,:) = get(gca, 'YLim');
        allAxes(nstr) = gca;
        title(char(exp_cond_names{nstr}))
        xlabel('time bin')
        ylabel('euclidian distance')
        
        for n_ep = 1:size(epoch_line,1)
            ind_ep = 1;
            if n_ep == 1
                line([abs(epoch_line{n_ep,2})*100 abs(epoch_line{n_ep,2})*100],[0 1],'Color','k');
                text(abs(epoch_line{n_ep,2})*100, max(allYLim(:)),epoch_line{n_ep,3});
                
            else
                line([(windows_index(ind_ep) + abs(epoch_line{n_ep,2})*100) (windows_index(ind_ep) + abs(epoch_line{n_ep,2})*100)],[0 1],'Color','k');
                text((windows_index(ind_ep) + abs(epoch_line{n_ep,2})*100), max(allYLim(:)),epoch_line{n_ep,3});
                ind_ep = ind_ep + 1;
            end
        end
    end
    
    set(allAxes, 'YLim', [min(allYLim(:)), max(allYLim(:))]);
    title_and_save(state_space_summary_handle,  [filename plot_title_part ', ' type_effector_short],plot_title,keys)
    
    %% latent variable
    
else
    %%
    %calculate and plot trajectories separately for control and
    %experimental condition
    condition_name = {'data_to_PCA_CT','data_to_PCA_PT'};
    for cd = size(condition_name)
        [coeff,score,latent,~,explained] = pca(eval(char(condition_name(cd))));
        plot_title_part = ' neural_state_space';
        fig_title=sprintf('Selection: %s %s %s hnd %s ch %s %s grouped by %s normalized %s in %s',...
            [Sel_for_title{:}],keys.monkey,[keys.conditions_to_plot{:}],mat2str(u_hands),mat2str(double(u_choice)),keys.arrangement,keys.ST.group_parameter, keys.ST.normalization,keys.ST.epoch_for_normalization);
        filename=sprintf('%s %s %s %s %s hnd %s ch %s N_%s %s %s',...
            [Sel_for_title{:}],keys.monkey,keys.ST.group_parameter,[keys.conditions_to_plot{:}],keys.arrangement(1:3),mat2str(u_hands),mat2str(double(u_choice)),keys.ST.normalization,keys.ST.epoch_for_normalization,char(condition_name(cd)));
        
        keys=ph_get_epoch_keys(keys,typ,eff,sum(type_effectors(:,1)==typ)>1);
        [~, type_effector_short] = get_type_effector_name(typ,eff);
        plot_title              = [fig_title type_effector_short plot_title_part];
        state_space_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
        if strcmp(char(condition_name(cd)),'data_to_PCA_CT')
            conditions_plot=conditions_hf_complete(1:4,:);
        else
            conditions_plot=conditions_hf_complete(5:8,:);
        end
        
        %loop to plot
        
        
        it_cond = 0;
        for cond = 1:size(conditions_plot,1)
            
            col=(conditions_plot(cond,1)+1)*6 + (conditions_plot(cond,3))*2 + (conditions_plot(cond,4)+1) + (conditions_plot(cond,5)*18);   %ritrieve color informatiom
            %             if cd == 1
            plot3(smooth(coeff(it_cond*condition_index+1:cond*condition_index,1)',0.3),smooth(coeff(it_cond*condition_index+1:cond*condition_index,2)',0.3),smooth(coeff(it_cond*condition_index+1:cond*condition_index,3)',0.3),'color',keys.line_colors(col,:),'LineWidth',1.5)
            %             else
            %                 plot3(smooth(coeff(it_cond*condition_index+1:cond*condition_index,1)',0.3),smooth(coeff(it_cond*condition_index+1:cond*condition_index,2)',0.3),smooth(coeff(it_cond*condition_index+1:cond*condition_index,3)',0.3),'color',keys.line_colors(col,:)/255,'LineWidth',1.5)
            %             end
            hold on
            it_cond = it_cond +1;
            legend_labels{cond} = char(legend_labels_hf(col));
        end
        legend(legend_labels)
        %estetique
        xlabel('PC 1')
        ylabel('PC 2')
        zlabel('PC 3')
        set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],'XTick',[],'YTick',[],'ZTick',[])
        title_and_save(state_space_summary_handle,  [filename plot_title_part ', ' type_effector_short],plot_title,keys)
    end
end

%%
    function title_and_save(figure_handle,filename,plot_title,keys)
        mtit(figure_handle,  plot_title, 'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none');
        stampit;
        ax = gca;
        ax.SortMethod='ChildOrder';
        texLabels = findall(figure_handle, 'type','text', 'FontWeight','bold');
        symbolIdx = ~cellfun('isempty',strfind({texLabels.String},'\'));
        set(texLabels(symbolIdx), 'FontWeight','normal');
        
        wanted_size=[50 30];
        set(figure_handle, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
        if keys.plot.export
            export_fig([keys.path_to_save, filename], '-pdf','-transparent'); % pdf by run
            close all
        end
    end

    function [dist] = calculate_3D_euc_dist(PA,PB)
        dist = sqrt((PB(1)-PA(1))^2 + (PB(2)-PA(2))^2 + (PB(3)-PA(3))^2);
    end

end



