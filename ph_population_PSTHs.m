function ph_population_PSTHs(population,modified_keys)
warning('off','MATLAB:catenate:DimensionMismatch');%
% %%% !!!!!!!! make sure epochs (baseline, cue, .... make sense for all effectors)

for fn=fieldnames(modified_keys)'
    keys.(fn{:})=modified_keys.(fn{:});
end

% legend_labels_hf={'NH IS IN' 'NH IS CH' 'IH IS IN' 'IH IS CH' 'CH IS IN' 'CH IS CH' ...
%     'NH VS IN' 'NH VS CH' 'IH VS IN' 'IH VS CH' 'CH VS IN' 'CH VS CH' ...
%     'NH CS IN' 'NH CS CH' 'IH CS IN' 'IH CS CH' 'CH CS IN' 'CH CS CH' ...
%     'NH IS IN P' 'NH IS CH P' 'IH IS IN P' 'IH IS CH P' 'CH IS IN P' 'CH IS CH P'...
%     'NH VS IN P' 'NH VS CH P' 'IH VS IN P' 'IH VS CH P' 'CH VS IN P' 'CH VS CH P'...
%     'NH CS IN P' 'NH CS CH P' 'IH CS IN P' 'IH CS CH P' 'CH CS IN P' 'CH CS CH P'};

legend_labels_hf={'NH CS IN' 'NH CS CH' 'CH CS IN' 'CH CS CH' 'IH CS IN' 'IH CS CH' ...
    'NH VS IN' 'NH VS CH' 'CH VS IN' 'CH VS CH' 'IH VS IN' 'IH VS CH' ...
    'NH IS IN' 'NH IS IH' 'CH IS IN' 'CH IS CH' 'IH IS IN' 'IH IS CH' ...
    'NH CS IN P' 'NH CS CH P' 'CH CS IN P' 'CH CS CH P' 'IH CS IN P' 'IH CS CH P'...
    'NH VS IN P' 'NH VS CH P' 'CH VS IN P' 'CH VS CH P' 'IH VS IN P' 'IH VS CH P'...
    'NH IS IN P' 'NH IS CH P' 'CH IS IN P' 'CH IS CH P' 'IH IS IN P' 'IH IS CH P'};

legend_labels_pref={'NH NP IN' 'NH NP CH' 'IH NP IN' 'IH NP CH' 'CH NP IN' 'CH NP CH' ...
    'NH PF IN' 'NH PF CH' 'IH PF IN' 'IH PF CH' 'CH PF IN' 'CH PF CH' ...
    'NH NP IN P' 'NH NP CH P' 'IH NP IN P' 'IH NP CH P' 'CH NP IN P' 'CH NP CH P'  ...
    'NH PF IN P' 'NH PF CH P' 'IH PF IN P' 'IH PF CH P' 'CH PF IN P' 'CH PF CH P'};

cols=keys.colors;

%% there are just too many colors once we include vertical targets, so for now we just keep use the same ones again...
keys.line_colors=[[cols.NH_IS_IN;cols.NH_IS_CH;cols.IH_IS_IN;cols.IH_IS_CH;cols.CH_IS_IN;cols.CH_IS_CH;...
    cols.NH_VS_IN;cols.NH_VS_CH;cols.IH_VS_IN;cols.IH_VS_CH;cols.CH_VS_IN;cols.CH_VS_CH;...
    cols.NH_CS_IN;cols.NH_CS_CH;cols.IH_CS_IN;cols.IH_CS_CH;cols.CH_CS_IN;cols.CH_CS_CH;]/255;...
    [cols.NH_IS_IN;cols.NH_IS_CH;cols.IH_IS_IN;cols.IH_IS_CH;cols.CH_IS_IN;cols.CH_IS_CH;...
    cols.NH_VS_IN;cols.NH_VS_CH;cols.IH_VS_IN;cols.IH_VS_CH;cols.CH_VS_IN;cols.CH_VS_CH;...
    cols.NH_CS_IN;cols.NH_CS_CH;cols.IH_CS_IN;cols.IH_CS_CH;cols.CH_CS_IN;cols.CH_CS_CH;]/610]; %%temporary for inactivation
keys.pref_colors=[[cols.NH_IS_IN;cols.NH_IS_CH;cols.IH_IS_IN;cols.IH_IS_CH;cols.CH_IS_IN;cols.CH_IS_CH;...
    cols.NH_CS_IN;cols.NH_CS_CH;cols.IH_CS_IN;cols.IH_CS_CH;cols.CH_CS_IN;cols.CH_CS_CH;]/255;...
    [cols.NH_IS_IN;cols.NH_IS_CH;cols.IH_IS_IN;cols.IH_IS_CH;cols.CH_IS_IN;cols.CH_IS_CH;...
    cols.NH_CS_IN;cols.NH_CS_CH;cols.IH_CS_IN;cols.IH_CS_CH;cols.CH_CS_IN;cols.CH_CS_CH;]/610]; %%temporary for inactivation

%% tuning table preparation and grouping
[tuning_per_unit_table]                 = ph_load_extended_tuning_table(keys);
[tuning_per_unit_table, Sel_for_title]  = ph_reduce_tuning_table(tuning_per_unit_table,keys);
if keys.PO.FR_subtract_baseline
    Sel_for_title =[Sel_for_title,{'base';'=';keys.PO.epoch_BL;', '}];
    keys.FR_subtract_baseline=0;
    population=ph_epochs(population,keys); %to undo baseline subtraction from FR per epochs, so that we can subtract in a meaningful way for population PSTH
end
idx_group_parameter=DAG_find_column_index(tuning_per_unit_table,keys.PO.group_parameter);
idx_unitID=DAG_find_column_index(tuning_per_unit_table,'unit_ID');
idx_RF_frame=DAG_find_column_index(tuning_per_unit_table,keys.PO.RF_frame_parameter);
group_values=tuning_per_unit_table(:,idx_group_parameter);
group_values=cellfun(@num2str, group_values, 'UniformOutput', false);
cell_in_any_group=[false; ~ismember(group_values(2:end),keys.PO.group_excluded)];
unique_group_values=unique(group_values(cell_in_any_group));
if isempty(unique_group_values)
    disp('no relevant groups found');
    return;
end
complete_unit_list={population.unit_ID}';

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
fitsettings.sd_x_min_ratio=0.125;%0.125;
fitsettings.sd_max_y=fitsettings.sd_max_x;
fitsettings.sd_xy_min_ratio=0.25;
fitsettings.sd_xy_max_ratio=1;
fitsettings.sd_y_min_ratio=fitsettings.sd_x_min_ratio;
fitsettings.xout=[-30:30]; % range for gaussian fit
fitsettings.yout=[-15:15];
fitsettings.range_factor=1;
fitsettings.fittypes=keys.PO.fittypes;

%% define conditions to look at
all_trialz=[population.trial];
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

u_hemifields=unique(per_trial.hemifield); %[-1,0,1]; % why does this have to be hardcoded? ---> Because case not defined yet, case defines positions !!

u_types     =unique(per_trial.types);
u_effectors =unique(per_trial.effectors);
u_hands     =unique(per_trial.hands);
u_choice    =unique(per_trial.choice);
u_perturbation    =unique(per_trial.perturbation);
u_perturbation=u_perturbation(~isnan(u_perturbation));

%% limit conditions key?
if ~any(keys.tt.hands==0) % cause hands 0 is any hand
u_hands     =u_hands(ismember(u_hands,keys.tt.hands));
end
u_choice    =u_choice(ismember(u_choice,keys.tt.choices));

all_type_effectors      = combvec(u_types,u_effectors)';
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
u_types     =unique(type_effectors(:,1))';
u_effectors =unique(type_effectors(:,2))';

condition_matrix            = combvec(u_hands,u_choice, u_perturbation,u_hemifields)';
conditions_out              = combvec(u_effectors,u_hands,u_choice, u_perturbation)';
conditions_hf               = combvec(u_hemifields,conditions_out')';
conditions_hf_complete      = combvec(u_hemifields,conditions_out')';
conditions_pref             = combvec([1 0],conditions_out')';

if any(u_hands~=0) && any(u_perturbation==1) %splitting to all 4 hand space conditions if hands are involved
    [~,~,columns_hf] = unique(conditions_hf(:,[1,3]),'rows');
    [~,~,columns_pref] = unique(conditions_pref(:,[1,3]),'rows');
else
    columns_hf          = ones(size(conditions_hf,1),1);
    columns_pref        = ones(size(conditions_pref,1),1);
end
%
% if any(u_perturbation==1) % temporary, better solution
%     if strcmp(keys.arrangement,'hands_inactivation_in_ch')
%         [~,~,columns_hf] = unique(conditions_hf(:,[1,3]),'rows');
%     end
%     condition_matrix(condition_matrix(:,1)==1 & condition_matrix(:,2)==1 & condition_matrix(:,4)==1,:)=[];
%     condition_matrix(condition_matrix(:,1)==2 & condition_matrix(:,2)==1 & condition_matrix(:,4)==-1,:)=[];
%     conditions_hf(conditions_hf(:,1)==1 & conditions_hf(:,3)==1 & conditions_hf(:,4)==1,:)=[];
%     conditions_hf(conditions_hf(:,1)==-1 & conditions_hf(:,3)==2 & conditions_hf(:,4)==1,:)=[];
% else
%     columns_hf          = ones(size(conditions_hf,1),1);
%     columns_pref        = ones(size(conditions_pref,1),1);
% end




%% finding positions and fixations
positions=unique(vertcat(whatisthis.trial.position),'rows');

%adjust fixations... should be done inside ph_arrange_positions_and_plots
%!?
fixations_temp=unique(vertcat(whatisthis.trial.fixation),'rows');
fix_temp_idx=true(size(fixations_temp,1),1);
for x=1:size(fixations_temp,1)-1
    if any(all(abs(bsxfun(@minus,fixations_temp(x+1:end,:),fixations_temp(x,:)))<4,2)) %% precision....
        fix_temp_idx(x)=false;
    end
end
fixations=fixations_temp(fix_temp_idx,:);
clear whatisthis

%% Convert to ipsi/contra, Baseline subtraction, normalization, re-ordering, and gaussian RFs
for tye=1:size(type_effectors,1)
    typ=type_effectors(tye,1);
    eff=type_effectors(tye,2);
    conditions_eff=conditions_out(conditions_out(:,1)==eff,2:end);
    t=find(u_types==typ);
    e=find(u_effectors==eff);
    keys=ph_get_epoch_keys(keys,typ,eff,sum(type_effectors(:,1)==typ)>1);
    
    % find epochs... here is still a problem, especially next line
    if ~isfield(keys.PO,'epoch_GB')
        keys.PO.epoch_GB=[keys.ANOVAS_PER_TYPE(typ).epoch{strcmp(keys.ANOVAS_PER_TYPE(typ).epoch(:,2),keys.PO.epoch_RF),1}]; %% gaussian baseline depends on the epoch for resposne field
    end
    if isempty(keys.PO.epoch_GB); keys.PO.epoch_GB={''}; end
    DN_epoch     =find(ismember(keys.EPOCHS(:,1),keys.PO.epoch_for_normalization));
    RF_epoch   =find(ismember(keys.EPOCHS(:,1),keys.PO.epoch_RF));
    BL_epoch          =find(ismember(keys.EPOCHS(:,1),keys.PO.epoch_BL));
    PF_epoch        =find(ismember(keys.EPOCHS(:,1),keys.PO.epoch_PF));
    gaussian_bl_epoch       =find(ismember(keys.EPOCHS(:,1),keys.PO.epoch_GB));
    
    pref_valid(t,:)=zeros(size(complete_unit_list,1),1);
    norm_factor=NaN(size(complete_unit_list,1),size(condition_matrix,1));
    for g=1:numel(unique_group_values)
        unitidx=ismember(complete_unit_list,tuning_per_unit_table(ismember(group_values,unique_group_values(g)),idx_unitID));
        units=find(all(unitidx,2))';
        for u=units
            tr_con=ismember([population(u).trial.completed],keys.cal.completed) & [population(u).trial.accepted];
            [pop]=ph_arrange_positions_and_plots(keys,population(u).trial(tr_con),population(u));
            
            [pop.trial(ismember([pop.trial.perturbation], keys.cal.perturbation_groups{1})).perturbation]=deal(0);
            [pop.trial(ismember([pop.trial.perturbation], keys.cal.perturbation_groups{2})).perturbation]=deal(1);
            pop=ph_LR_to_CI(keys,pop);
            per_epoch=vertcat(pop.trial.epoch);
            if keys.PO.FR_subtract_baseline
                baseline=repmat(vertcat(per_epoch(:,BL_epoch).FR),1,size(condition_matrix,1));
            else
                baseline=zeros(numel(pop.trial),size(condition_matrix,1));
            end
            %% normalization factor
            for n=1:size(condition_matrix,1)
                clear trpar
                switch keys.PO.normalization
                    case 'by_condition'
                        condition_parameters_tmp=[condition_parameters {'hemifield'}];
                        for par=1:numel(condition_parameters_tmp)
                            fn=condition_parameters_tmp{par};
                            trpar(par,:)=[pop.trial.(fn)]==condition_matrix(n,par);
                        end
                        trpar(end+1,:)=[pop.trial.type]==typ & [pop.trial.effector]==eff;
                    case 'by_perturbation'
                        trpar=[pop.trial.type]==typ & [pop.trial.effector]==eff & [pop.trial.perturbation]==condition_matrix(n,strcmp(condition_parameters,'perturbation'));
                    case 'by_effector'
                        trpar=[pop.trial.type]==typ & [pop.trial.effector]==eff;
                    case 'by_type'
                        trpar=[pop.trial.type]==typ;
                    case 'by_all_trials'
                        trpar=true(size([pop.trial]));
                    case 'none' %% really doesnt need an extra rule actually?? Not here at least..
                        trpar=[pop.trial.type]==typ;
                    case 'z_score' %% really doesnt need an extra rule actually?? Not here at least..
                        trpar=[pop.trial.type]==typ & [pop.trial.effector]==eff;
                    case 'percent_change' %% really doesnt need an extra rule actually?? Not here at least..
                        trpar=[pop.trial.type]==typ & [pop.trial.effector]==eff;
                end
                tr_con=all(trpar,1);
                if ~any(tr_con) %
                    continue
                end
                per_epoch=vertcat(pop.trial(tr_con).epoch);
                if keys.PO.FR_subtract_baseline && ~keys.PO.baseline_per_trial
                    baseline(:,n)=deal(nanmean(vertcat(per_epoch(:,BL_epoch).FR)));
                end
                
                %what to do if per_epoch is empty or
                if strcmp(keys.PO.normalization,'none') || isempty(DN_epoch) %DN_epoch is empty? %% check if this works, for no multiplicative normalization
                    norm_factor(u,n)= 1;
%                 elseif isempty(per_epoch)
%                     norm_factor(u,n)=NaN;
                elseif strcmp(keys.PO.normalization,'percent_change')
                    norm_factor(u,n)=deal(nanmean([per_epoch(:,BL_epoch).FR per_epoch(:,DN_epoch).FR NaN]))/100;
                    baseline(:,n)=deal(nanmean([per_epoch(:,BL_epoch).FR per_epoch(:,DN_epoch).FR NaN]));
                elseif ~isempty(DN_epoch)
                    norm_factor(u,n)=nanmean([per_epoch(:,DN_epoch).FR NaN]);
                end
            end
            
            if all(isnan(norm_factor(u,:))) &&  all(all(isnan(baseline)))
                continue
            end
            
            %% z-scoring (way too complicated implementation)
            if strcmp(keys.PO.normalization,'z_score')
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
            %% correct normalization factors if they are too low
            if any(norm_factor(u,:)<1)
                baseline=baseline+1-nanmean(norm_factor(u,:));
                norm_factor(u,:)=deal(1);
            end
            
            
            %% preferred and unpreferred location (taken from instructed trials only! (?) Not necessary though
            tr_con=[pop.trial.type]==typ & [pop.trial.effector]==eff;
            FR_for_pref=NaN(size(positions,1),1);
            per_epoch=vertcat(pop.trial(tr_con).epoch);
            for p=1:size(positions,1)
                tr_hemi=all(abs(bsxfun(@minus,vertcat(pop.trial.position),positions(p,:)))<1.5,2) & [pop.trial.choice]'==0;
                FR_for_pref(p)=nanmean([per_epoch((tr_hemi(tr_con)),PF_epoch).FR]);
            end
            if ischar(unique_group_values{g}) && strcmp(unique_group_values{g},'su') %very specific rule, be aware of this one!
                [~,pref_idx]=min(FR_for_pref);
            else
                [~,pref_idx]=max(FR_for_pref);
            end
            [~,unpref_idx]=min(abs(nanmean([positions(:,1)+positions(pref_idx,1) positions(:,2)-positions(pref_idx,2)],2)));
            pref_valid(t,u)=true;
            for ch=u_choice
                pref_valid(t,u)=pref_valid(t,u) && pref_idx~=unpref_idx && ...
                    sum(all(abs(bsxfun(@minus,vertcat(pop.trial(tr_con).position),positions(pref_idx,:)))<1.5,2)   & [pop.trial(tr_con).choice]'==ch) >=keys.cal.min_trials_per_condition && ...
                    sum(all(abs(bsxfun(@minus,vertcat(pop.trial(tr_con).position),positions(unpref_idx,:)))<1.5,2) & [pop.trial(tr_con).choice]'==ch) >=keys.cal.min_trials_per_condition;
            end
            
            %% average PSTH per unit
            for ce=1:size(conditions_eff,1)
                %c=ce;
                c=find(ismember(conditions_out,[eff conditions_eff(ce,:)],'rows')); %% c is conditions, n is including hemifields -1,0,1
                clear trpar
                
                for par=1:numel(condition_parameters)
                    fn=condition_parameters{par};
                    trpar(par,:)=[pop.trial.(fn)]==conditions_eff(ce,par); %condition_matrix_wohf?
                end
                trpar(end+1,:)=[pop.trial.type]==typ & [pop.trial.effector]==eff;
                tr_con=all(trpar,1);
                per_epoch=vertcat(pop.trial(tr_con).epoch);
                
                for w=1:size(keys.PSTH_WINDOWS,1)
                    for f=1:numel(u_hemifields) %hemifield
                        tr_hemi=[pop.trial.hemifield]==u_hemifields(f); 
                        ix = tr_con & tr_hemi;
                        n=max([1,find(ismember(condition_matrix(:,1:3),conditions_eff(ce,:),'rows') & condition_matrix(:,end)==u_hemifields(f))]);
                        condition(t,c).per_hemifield(f).window(w).unit(u).average_spike_density=...
                            ph_spike_density(pop.trial(ix),w,keys,baseline(ix,n),norm_factor(u,n));
                        condition(t,c).per_hemifield(f).effector=eff;
                        for x=1:size(fixations,1)
                            tr_fix=all(abs(bsxfun(@minus,vertcat(pop.trial.fixation),fixations(x,:)))<4,2)'; %4 is precision --> add to keys?
                            ix = tr_con & tr_hemi & tr_fix;
                            condition(t,c).per_hf_fixation(f,x).window(w).unit(u).average_spike_density=...
                                ph_spike_density(pop.trial(ix),w,keys,baseline(ix,n),norm_factor(u,n));
                            condition(t,c).per_hf_fixation(f,x).fixation=fixations(x,:);
                            condition(t,c).per_hf_fixation(f,x).effector=eff;
                            condition(t,c).per_hf_fixation(f,x).position=[u_hemifields(f) 0];
                            condition(t,c).per_hf_fixation(f,x).sign.unit(u)=1;
                        end
                    end
                    
                    for p=1:size(positions,1)
                        f=sign(positions(p,1)); %what happens if position is neither left nor right??
                        n=max([1,find(ismember(condition_matrix(:,1:3),conditions_eff(ce,:),'rows') & condition_matrix(:,end)==f)]);
                        tr_pos=all(abs(bsxfun(@minus,vertcat(pop.trial.position),positions(p,:)))<1.5,2)';
                        ix = tr_con & tr_pos;
                        
                        condition(t,c).per_position(p).window(w).unit(u).average_spike_density= ph_spike_density(pop.trial(ix),w,keys,baseline(ix,n),norm_factor(u,n));
                        condition(t,c).per_position(p).position=positions(p,:);
                        condition(t,c).per_position(p).fixation=1;
                        condition(t,c).per_position(p).effector=eff;
                        % problem here if tr_con is empty (which it can be,
                        % because we are not requiring every condition to
                        % be valid for all units
                        ix = tr_pos(tr_con);
                        if any(ix) && ttest([per_epoch(ix,PF_epoch).FR], [per_epoch(ix,BL_epoch).FR])==1
                            condition(t,c).per_position(p).sign.unit(u)=sign(mean([per_epoch(ix,PF_epoch).FR]- [per_epoch(ix,BL_epoch).FR]));
                        else
                            condition(t,c).per_position(p).sign.unit(u)=0;
                        end
                        %% supposedly matches with target position precision...
                        for x=1:size(fixations,1)
                            tr_fix=all(abs(bsxfun(@minus,vertcat(pop.trial.fixation),fixations(x,:)))<4,2)';
                            ix = tr_con & tr_pos & tr_fix;
                            
                            % within condition normalization here...
                            % norm_fact=nanmean([per_epoch(trpos(tr_con) & trfix(tr_con),BL_epoch).FR per_epoch(trpos(tr_con) & trfix(tr_con),DN_epoch).FR NaN])/100;
                            % base_line=nanmean([per_epoch(trpos(tr_con) & trfix(tr_con),BL_epoch).FR per_epoch(trpos(tr_con) & trfix(tr_con),DN_epoch).FR NaN]);
                            % condition(t,c).per_position_fixation(p,x).window(w).unit(u).average_spike_density=...
                            % ph_spike_density(pop.trial(tr_con & trpos & trfix),w,keys,base_line,norm_fact);
                            
                            condition(t,c).per_position_fixation(p,x).window(w).unit(u).average_spike_density=...
                                ph_spike_density(pop.trial(ix),w,keys,baseline(ix,n),norm_factor(u,n));
                            condition(t,c).per_position_fixation(p,x).position=positions(p,:);
                            condition(t,c).per_position_fixation(p,x).fixation=fixations(x,:);
                            condition(t,c).per_position_fixation(p,x).effector=eff;
                            condition(t,c).per_position_fixation(p,x).sign.unit(u)=0;
                            % problem here if tr_con is empty (which it can be,
                            % because we are not requiring every condition to
                            % be valid for all units
                            ix = tr_pos(tr_con) & tr_fix(tr_con);
                            if any(ix) && ttest([per_epoch(ix,PF_epoch).FR], [per_epoch(ix,BL_epoch).FR])==1
                                condition(t,c).per_position_fixation(p,x).sign.unit(u)=sign(mean([per_epoch(ix,PF_epoch).FR]- [per_epoch(ix,BL_epoch).FR]));
                            end
                        end
                        ix = tr_con & tr_pos;
                        if p==pref_idx                            
                            condition(t,c).per_preference(1).window(w).unit(u).average_spike_density=...
                                ph_spike_density(pop.trial(ix),w,keys,baseline(ix,n),norm_factor(u,n));
                            condition(t,c).per_preference(1).effector=eff;
                        end
                        if p==unpref_idx
                            condition(t,c).per_preference(2).window(w).unit(u).average_spike_density=...
                                ph_spike_density(pop.trial(ix),w,keys,baseline(ix,n),norm_factor(u,n));
                            condition(t,c).per_preference(2).effector=eff;
                        end
                    end
                end
                
                %% gaussian response fields (needs cleanup)
                if ~keys.PO.plot_RF
                    continue;
                end
                
                zin_per_pos=NaN(size(positions,1),1);
                
                gaussian_baseline=vertcat(per_epoch(:,gaussian_bl_epoch).FR);
                if isempty(gaussian_baseline)
                    gaussian_baseline=zeros(size(per_epoch,1),1);
                end
                %% supposedly matches with target position precision...
                for p=1:size(positions,1) %for FR plot only
                    tr=all(abs(bsxfun(@minus,vertcat(pop.trial(tr_con).position),positions(p,:)))<1.5,2);
                    if sum(tr)==0; continue; end
                    zin_per_pos(p)=nanmean(vertcat(per_epoch(tr,RF_epoch).FR)-gaussian_baseline(tr));
                end
                gaussian_positions=vertcat(pop.trial(tr_con).position);
                zin=vertcat(per_epoch(:,RF_epoch).FR);
                
                %%double because single in Linus/Curius Pulvinar gaze ?
                %nonanidx=~isnan(zin_per_pos);
                %FR_tmp=struct('FR',num2cell(zin_per_pos(nonanidx)),'x',num2cell(positions(nonanidx,1)),'y',num2cell(positions(nonanidx,2)));
                FR_tmp=struct('FR',num2cell(zin_per_pos),'x',num2cell(positions(:,1)),'y',num2cell(positions(:,2)));
                RF_tmp=ph_fit_target_positions_2D(gaussian_positions(:,1),gaussian_positions(:,2),zin,gaussian_baseline,fitsettings);
                
                condition(t,c).fitting.unit(u).parameters=RF_tmp;
                condition(t,c).fitting.unit(u).positions =FR_tmp;
                
            end
        end
    end
end

%% condition comparison
comparisons_per_effector(1).reach_hand{1}=0;
comparisons_per_effector(1).reach_hand{2}=0;
comparisons_per_effector(1).hemifield{1}=[-1];
comparisons_per_effector(1).hemifield{2}=[1];
comparisons_per_effector(1).choice{1}=0;
comparisons_per_effector(1).choice{2}=0;
comparisons_per_effector(1).color=[1 0 0];
condition_parameters_comparissons = [{'hemifield'} {'effector'} condition_parameters];

% effector comparison only possible if several_effectors
for t=1:size(condition,1)
    current=[condition(t,:).per_hemifield];
    current_window=vertcat(current.window);
    for w=1:numel(current(1).window) %% epochs might be different for different effectors.......... !!!!====???
        for g=1:numel(unique_group_values)
            unitidx=ismember(complete_unit_list,tuning_per_unit_table(ismember(group_values,unique_group_values(g)),idx_unitID));
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
                n=0;
                for eff=u_effectors
                    n=n+1;
                    cM1(conditions_hf(:,2)~=eff,2)=false;
                    cM2(conditions_hf(:,2)~=eff,2)=false;
                    c1=find(all(cM1,2));
                    c2=find(all(cM2,2));
                    sigbins(t).group(g).per_effector(n).window(w).bins(comp,:)=ph_compare_by_bin(current_window(:,w),c1,c2,units);
                end
            end
        end
    end
end

%% plots

logscale_127=(1:127)'*255/127;
logscale_255=(log(255)-log(255:-1:129)')*255/(log(255)-log(127));
if ~isempty(gaussian_bl_epoch)
    RF_colormap=[logscale_127 logscale_127 ones(127,1)*255; 255 255 255; ones(127,1)*255 flipud(logscale_127) flipud(logscale_127)]/256;
else
    RF_colormap=[255*ones(1,255);255:-1:1;255:-1:1]'/255;
end

for t=1:size(condition,1)
    typ=u_types(mod(t-1,numel(u_types))+1);
    
    fig_title=sprintf('%s %s %s hnd %s ch %s %s normalized %s in %s grouped by %s ',...
        keys.monkey,[keys.conditions_to_plot{:}],keys.arrangement,     mat2str(keys.tt.hands),mat2str(double(keys.tt.choices)),[Sel_for_title{:}],keys.PO.normalization,keys.PO.epoch_for_normalization,keys.PO.group_parameter);
    if strcmp(keys.PO.normalization,'percent_change')
        filename=sprintf('%s %s %s hnd %s ch %s %s N_prct %s to %s %s ',...
            keys.monkey,[keys.conditions_to_plot{:}],keys.arrangement(1:3),mat2str(keys.tt.hands),mat2str(double(keys.tt.choices)),[Sel_for_title{:}],keys.PO.epoch_BL,keys.PO.epoch_for_normalization,keys.PO.group_parameter);
    else
        filename=sprintf('%s %s %s hnd %s ch %s %s N_%s %s %s ',...
            keys.monkey,[keys.conditions_to_plot{:}],keys.arrangement(1:3),mat2str(keys.tt.hands),mat2str(double(keys.tt.choices)),[Sel_for_title{:}],keys.PO.normalization,keys.PO.epoch_for_normalization,keys.PO.group_parameter);
    end
    
    %% PSTH plot
    current=[condition(t,:).per_hemifield];
    conditions_PSTH=conditions_hf_complete;
    legend_labels=legend_labels_hf;
    plot_title_part=' PSTHs';
    units_valid=ones(size(complete_unit_list,1),1);
    column_indexes=columns_hf;
    plot_PSTH_no_empties
    
    %% PSTH preferred and unpreferred plot
    current=[condition(t,:).per_preference];
    conditions_PSTH=conditions_pref;
    legend_labels=legend_labels_pref;
    plot_title_part=[' pref in ' keys.PO.epoch_PF ' PSTHs'];
    units_valid=pref_valid(t,:)';
    column_indexes=columns_pref;
    if ~any(u_perturbation==1)
        plot_PSTH_no_empties
    end
    
    if   keys.PO.plot_per_position
        unique_group_values_tmp=unique_group_values;
        for gt=1:numel(unique_group_values_tmp)
            unique_group_values=unique_group_values_tmp(gt);
            %% PSTH per position plot
            current=[condition(t,:).per_position_fixation];
            current=current(:);
            conditions_PSTH=conditions_pref; %%???
            legend_labels={'-15', '0', '+15'};
            plot_title_part=['=' unique_group_values{1} ' PSTHs per position'];
            units_valid=ones(size(complete_unit_list,1),1);
            column_indexes=columns_pref;
            plot_PSTH_no_empties
            
            %% PSTH per position plot
            current=[condition(t,:).per_position_fixation];
            current=current(:);
            conditions_PSTH=conditions_pref; %%???
            legend_labels={'-15', '0', '+15'};
            plot_title_part=['=' unique_group_values{1} ' PSTHs per position F ensu'];
            units_valid=ones(size(complete_unit_list,1),1);
            column_indexes=columns_pref;
            plot_PSTH_no_empties
            
            current=[condition(t,:).per_position];
            current=current(:);
            conditions_PSTH=conditions_pref; %%???
            legend_labels={'-15', '0', '+15'};
            plot_title_part=['=' unique_group_values{1} ' PSTHs per position ensu'];
            units_valid=ones(size(complete_unit_list,1),1);
            column_indexes=columns_pref;
            plot_PSTH_no_empties
            
            current=[condition(t,:).per_hf_fixation];
            current=current(:);
            conditions_PSTH=conditions_pref; %%???
            legend_labels={'-15', '0', '+15'};
            plot_title_part=['=' unique_group_values{1} ' PSTHs per position hf'];
            units_valid=ones(size(complete_unit_list,1),1);
            column_indexes=columns_pref;
            plot_PSTH_no_empties
        end
        unique_group_values=unique_group_values_tmp;
    end
    
    %% RF and FR plots
    if ~keys.PO.plot_RF
        continue;
    end
    angles=[0:pi/100:(2+1/100)*pi];
    circle_x=cos(angles);
    circle_y=sin(angles);
    for g=1:numel(unique_group_values)
        unitidx=ismember(complete_unit_list,tuning_per_unit_table(ismember(group_values,unique_group_values(g)),idx_unitID));
        group_units=find(all(unitidx,2))';
        if isempty(keys.PO.RF_columns) || isempty(keys.PO.RF_rows)
            RF_columns=ceil(sqrt(numel(group_units)+1));
            RF_rows=ceil(sqrt(numel(group_units)+1));
        else
            RF_columns=keys.PO.RF_columns;
            RF_rows=keys.PO.RF_rows;
        end
        [~,complete_list_tuning_table_idx]= ismember(complete_unit_list,tuning_per_unit_table(:,idx_unitID));
        RF_frame_entries=tuning_per_unit_table(complete_list_tuning_table_idx(group_units),idx_RF_frame);
        for c=1:size(condition,2)
            per_unit=[condition(t,c).fitting.unit(group_units)];
            RFparameters            =[per_unit.parameters];
            pos=[per_unit.positions];
            bestfits={RFparameters.bestfit};
            secondbestfits={RFparameters.secondbestfit};
            
            [~,~,idx_fittype]=unique({RFparameters.bestfit});
            idx_fittype(ismember(bestfits,'none'))=100;
            idx_within_fittype=inf(size(idx_fittype));
            
            fittypes=keys.PO.fittypes;
            fittypes_all=[fittypes,'none'];
            RFsizes=zeros(size(RFparameters));
            for f=1:numel(fittypes_all)
                fittype=fittypes_all{f};
                fitidx.(fittype)=ismember(bestfits,fittype);
                fitidx2nd.(fittype)=ismember(secondbestfits,fittype);
                
                Parameters_temp=[RFparameters(fitidx.(fittype)).(fittype)];
                if isempty(Parameters_temp)
                    continue;
                end
                switch fittype
                    case {'sigmoidal','linear'}
                        [~, SC] =sort([Parameters_temp.phi]);
                    case {'gaussian1'}
                        RFsizes(fitidx.(fittype))                =2*2*sqrt([Parameters_temp.sx].*[Parameters_temp.sy]);
                        [~, SC] =sort(RFsizes(fitidx.(fittype)));
                        SC      =fliplr(SC);
                        SC=SC+1000*[Parameters_temp.zmax]>0;
                    case {'gaussian2','gaussian15'}
                        Zmax                     =[Parameters_temp.zmax1; Parameters_temp.zmax2];
                        RFsizes2D                  =2*[2*sqrt([Parameters_temp.sx1].*[Parameters_temp.sy1]);2*sqrt([Parameters_temp.sx2].*[Parameters_temp.sy2])];
                        [~,zone_idx]                =max(RFsizes2D,[],1);  %%
                        zone_idx=sub2ind(size(RFsizes2D),zone_idx,1:size(RFsizes2D,2));
                        RFsizes(fitidx.(fittype))                =RFsizes2D(zone_idx);
                        [~, SC] =sort(RFsizes(fitidx.(fittype)));
                        SC      =fliplr(SC);
                        SC=SC+1000*Zmax(zone_idx)>0;
                end
                if ~isempty(idx_RF_frame)
                    [~,SC]=ismember(RF_frame_entries(fitidx.(fittype)),keys.PO.RF_frame_entries);
                    switch fittype
                        case {'gaussian1'}
                            SC=SC+1000*([Parameters_temp.zmax]'>0);
                    end
                    
                end
                idx_within_fittype(fitidx.(fittype))=SC;
                
            end
            [~, sort_by_size_index] =sort(RFsizes);
            sort_by_size_index      =fliplr(sort_by_size_index);
            RF_sorting_matrix      =[idx_fittype;idx_within_fittype]';
            [~, RF_sort_index] = sortrows(RF_sorting_matrix);
            [~, RF_sort_index] = sort(RF_sort_index);
            
            %normalize firing rates?
            for u=1:size(pos,2)
                FR=[pos(:,u).FR];
                if ~isempty(gaussian_bl_epoch)
                    [FRmax(u), maxposition(u)]=max([FR(FR>0) 0]);
                    FRmin=min([FR(FR<0) 0]);
                    FRmax(u)=max([abs(FRmax(u)) abs(FRmin)]);
                    
                    FR255=num2cell(round((FR+FRmax(u))/2/FRmax(u)*254)+1);
                else
                    [FRmax(u), maxposition(u)]=max(FR);
                    FR255=num2cell(round(FR/FRmax(u)*254)+1);
                end
                [pos(:,u).FR255_GAU]=deal(FR255{:});
                [pos(:,u).FR255]=deal(FR255{:});
            end
            
            
            %% RF plot
            plot_title_part        = ['=' unique_group_values{g} ' con' num2str(c) ' RF in ' keys.PO.epoch_RF ' BL ' keys.PO.epoch_GB];
            f_handle              = figure('units','normalized','outerposition',[0 0 1 1],'name',[fig_title plot_title_part]);
            colormap(RF_colormap);
            
            
            for u=1:numel(group_units)
                n=RF_sort_index(u);
                subplot(RF_rows,RF_columns,n);
                %title(population(group_units(u)).unit_ID,'interpreter','none');
                RF=RFparameters(u);
                %
                if ~isempty(gaussian_bl_epoch)
                    Zout=round(RF.Zout/FRmax(u)/2*254 + 127)+1;
                else
                    Zout=round(RF.Zout/FRmax(u)*254)+1;
                end
                Zout(Zout>255)=255;
                Zout(Zout<1)=1;
                
                
                kk=cat(3,Zout);
                image(fitsettings.xout,-fitsettings.yout,rot90(nanmean(kk,3)))
                caxis([1 255]);
                
                %% two ellipses
                hold on
                BF=RF.(RF.bestfit);
                switch RF.bestfit
                    case 'gaussian1'
                        center=[BF.xmax BF.ymax];
                        ellipse_r=2*BF.sx*BF.sy./sqrt(BF.sy.^2*cos(angles - BF.phi).^2+BF.sx.^2*sin(angles  - BF.phi).^2);
                        ellipse_x = circle_x.*ellipse_r; ellipse_y = circle_y.*ellipse_r;
                        line(ellipse_x+center(1),ellipse_y+center(2),'linewidth',1,'color','k');
                    case {'gaussian2','gaussian15'}
                        center=[BF.xmax1 BF.ymax1];
                        ellipse_r=2*BF.sx1*BF.sy1./sqrt(BF.sy1.^2*cos(angles - BF.phi1).^2+BF.sx1.^2*sin(angles  - BF.phi1).^2);
                        ellipse_x = circle_x.*ellipse_r; ellipse_y = circle_y.*ellipse_r;
                        line(ellipse_x+center(1),ellipse_y+center(2),'linewidth',1,'color','k');
                        center=[BF.xmax2 BF.ymax2];
                        ellipse_r=2*BF.sx2*BF.sy2./sqrt(BF.sy2.^2*cos(angles - BF.phi2).^2+BF.sx2.^2*sin(angles  - BF.phi2).^2);
                        ellipse_x = circle_x.*ellipse_r; ellipse_y = circle_y.*ellipse_r;
                        line(ellipse_x+center(1),ellipse_y+center(2),'linewidth',1,'color','k');
                end
                
                if ~isempty(idx_RF_frame)
                    range_x=max(fitsettings.xout)-min(fitsettings.xout);
                    range_y=max(fitsettings.yout)-min(fitsettings.yout);
                    rh=rectangle('Position',[min(fitsettings.xout)+range_x/100 min(fitsettings.yout)+range_y/100  range_x*98/100 range_y*98/100]);
                    set(rh,'edgecolor',keys.PO.RF_frame_colors{ismember(keys.PO.RF_frame_entries,RF_frame_entries(u))}/256);
                    title([population(group_units(u)).unit_ID ' R2= ' num2str(round(RF.R2_adjusted*100)/100) ' ' RF_frame_entries{u}],'fontsize',3,'interpreter','none');
                else
                    title([population(group_units(u)).unit_ID ' R2= ' num2str(round(RF.R2_adjusted*100)/100) ],'fontsize',3,'interpreter','none');
                end
                axis equal
                set(gca,'Ydir','normal','Xtick',[],'Ytick',[],'xlim',[min(fitsettings.xout) max(fitsettings.xout)],'ylim',[min(fitsettings.yout) max(fitsettings.yout)]);
            end
            subplot(RF_rows,RF_columns,numel(group_units)+1);
            colorbar;
            ph_title_and_save(f_handle,  [filename plot_title_part],[fig_title plot_title_part],keys);
            
            
            %% fittype R2 values histogram
            plot_title_part       = ['=' unique_group_values{g} ' con' num2str(c) ' RF R2 ' 'in ' keys.PO.epoch_RF];
            f_handle              = figure('units','normalized','outerposition',[0 0 1 1],'name',[fig_title plot_title_part]);
            R2adjusted=[];
            for f=1:numel(fittypes)
                fittype=fittypes{f};
                fittemp=[RFparameters.(fittype)];
                R2adjusted_temp=hist([fittemp.R2],bins);
                R2adjusted=[R2adjusted;R2adjusted_temp];
                
            end
            bar(bins,R2adjusted','stacked');
            legend(fittypes);
            ph_title_and_save(f_handle,  [filename plot_title_part],[fig_title plot_title_part],keys);
            
            %% fittype R2 values histogram
            plot_title_part       = ['=' unique_group_values{g} ' con' num2str(c) ' RF R2 adjusted ' 'in ' keys.PO.epoch_RF];
            f_handle              = figure('units','normalized','outerposition',[0 0 1 1],'name',[fig_title plot_title_part]);
            bins=-1:0.05:1;
            R2adjusted=[];
            for f=1:numel(fittypes)
                fittype=fittypes{f};
                fittemp=[RFparameters.(fittype)];
                R2adjusted_temp=hist([fittemp.R2_adjusted],bins);
                R2adjusted=[R2adjusted;R2adjusted_temp];
                
            end
            bar(bins,R2adjusted','stacked');
            legend(fittypes);
            ph_title_and_save(f_handle,  [filename plot_title_part],[fig_title plot_title_part],keys);
            
            
            %% fittype R2 values histogram
            plot_title_part       = ['=' unique_group_values{g} ' con' num2str(c) ' RF R2 adjusted win ' 'in ' keys.PO.epoch_RF];
            f_handle              = figure('units','normalized','outerposition',[0 0 1 1],'name',[fig_title plot_title_part]);
            bins=-0.05:0.05:2;
            cols=colormap(jet(numel(fittypes)));
            for f=1:numel(fittypes)
                fittype=fittypes{f};
                subplot(numel(fittypes),1,f)
                hold on;
                title([fittype, ', N= ' num2str(sum(fitidx.(fittype)))]);
                R2adjusted=[];
                fittemp=[RFparameters.(fittype)];
                fittemp=fittemp((fitidx.(fittype)));
                for f2=1:numel(fittypes)
                    fittype2=fittypes{f2};
                    fittemp2=[RFparameters.(fittype2)];
                    fittemp2=fittemp2((fitidx.(fittype)));
                    tmp_idx=fitidx2nd.(fittype2);
                    tmp_idx=tmp_idx(fitidx.(fittype));
                    R2_differences=[fittemp(tmp_idx).R2_adjusted]-[fittemp2(tmp_idx).R2_adjusted];
                    R2adjusted_temp=hist(R2_differences,bins);
                    if f==f2
                        R2adjusted_temp=zeros(size(R2adjusted_temp));
                    end
                    R2adjusted=[R2adjusted;R2adjusted_temp];
                    mean_temp(f2)=mean(R2_differences);
                    median_temp(f2)=median(R2_differences);
                    sem_temp(f2)=sterr(R2_differences);
                end
                bar(bins,R2adjusted','stacked');
                y_lim=get(gca,'ylim');
                for f2=1:numel(fittypes)
                    plot([mean_temp(f2) mean_temp(f2)],y_lim,'color',cols(f2,:));
                    plot([mean_temp(f2)+sem_temp(f2) mean_temp(f2)+sem_temp(f2)],y_lim,':','color',cols(f2,:));
                    plot([mean_temp(f2)-sem_temp(f2) mean_temp(f2)-sem_temp(f2)],y_lim,':','color',cols(f2,:));
                    text(mean_temp(f2),y_lim(1),['mean: ' num2str(round(mean_temp(f2)*100)/100) ' + ' num2str(round(sem_temp(f2)*100)/100)],'rotation',90,'color',cols(f2,:));
                end
                legend(fittypes);
            end
            ph_title_and_save(f_handle,  [filename plot_title_part],[fig_title plot_title_part],keys);
            
            
            %% FR plot
            plot_title_part       = ['=' unique_group_values{g} ' con' num2str(c) ' FR ' 'in ' keys.PO.epoch_RF];
            f_handle              = figure('units','normalized','outerposition',[0 0 1 1],'name',[fig_title plot_title_part]);
            colormap(RF_colormap);
            for u=1:numel(group_units)
                n=RF_sort_index(u);
                FR255=[pos(:,u).FR255];
                subplot(RF_rows,RF_columns,n);
                title(population(group_units(u)).unit_ID,'fontsize',3,'interpreter','none');
                hold on;
                scatter([pos(:,u).x],[pos(:,u).y],25,FR255,'filled');
                caxis([1 255]);
                axis equal
                sp_position=get(gca,'Position');%sp_position(3)=sp_position(3)*1.4;sp_position(4)=sp_position(4)*1.4;
                set(gca,'Position',sp_position,'Ydir','normal','Xtick',[],'Ytick',[],'xlim',[min(fitsettings.xout) max(fitsettings.xout)],'ylim',[min(fitsettings.yout) max(fitsettings.yout)]);
            end
            subplot(RF_rows,RF_columns,numel(group_units)+1);
            colorbar
            ph_title_and_save(f_handle,  [filename plot_title_part],[fig_title plot_title_part],keys);
            
            %% FR peak histogram plot
            plot_title_part        = ['=' unique_group_values{g} ' con' num2str(c) ' FR in ' keys.PO.epoch_RF ' histogram' ];
            f_handle              = figure('units','normalized','outerposition',[0 0 1 1],'name',[fig_title plot_title_part]);
            hold on;
            colormap(RF_colormap);
            caxis([1 255]);
            for p=1:size(pos,1)
                N_per_pos(p)=sum(maxposition==p);
            end
            for p=1:size(pos,1)
                FRmean=round(N_per_pos(p)/max(N_per_pos)*254)+1;
                if ~isnan(FRmean)
                    plot(pos(p,1).x,pos(p,1).y,'o','markerfacecolor',RF_colormap(FRmean,:),'markeredgecolor','none','markersize',5);
                    text(pos(p,1).x,pos(p,1).y,num2str(N_per_pos(p)))
                end
            end
            axis equal
            colorbar;
            set(gca,'Xtick',[],'Ytick',[],'xlim',[min(fitsettings.xout) max(fitsettings.xout)],'ylim',[min(fitsettings.yout) max(fitsettings.yout)]);
            ph_title_and_save(f_handle,  [filename plot_title_part],[fig_title plot_title_part],keys);
            
            
            %% FR summary plot
            plot_title_part        = ['=' unique_group_values{g} ' con' num2str(c) ' FR in ' keys.PO.epoch_RF ' average' ];
            f_handle              = figure('units','normalized','outerposition',[0 0 1 1],'name',[fig_title plot_title_part]);
            hold on;
            colormap(RF_colormap);
            caxis([1 255]);
            for p=1:size(pos,1)
                FRmean(p)=nanmean([pos(p,:).FR]);
            end
            for p=1:size(pos,1)
                if ~isempty(gaussian_bl_epoch)
                    FRmeancolidx=round(FRmean(p)/FRmax(u)/2*254 + 127)+1;
                else
                    FRmeancolidx=round(FRmean(p)/FRmax(u)*254)+1;
                end
                if ~isnan(FRmean(p))
                    plot(pos(p,1).x,pos(p,1).y,'o','markerfacecolor',RF_colormap(FRmeancolidx,:),'markeredgecolor','none','markersize',5);
                    text(pos(p,1).x,pos(p,1).y,num2str(round(FRmean(p)*100)/100));
                end
            end
            axis equal
            colorbar;
            set(gca,'Xtick',[],'Ytick',[],'xlim',[min(fitsettings.xout) max(fitsettings.xout)],'ylim',[min(fitsettings.yout) max(fitsettings.yout)]);
            ph_title_and_save(f_handle,  [filename plot_title_part],[fig_title plot_title_part],keys);
            
            %% RF centers
            plot_title_part        = ['=' unique_group_values{g} ' con' num2str(c) ' RF  in ' keys.PO.epoch_RF  ' summary'];
            f_handle              = figure('units','normalized','outerposition',[0 0 1 1],'name',[fig_title plot_title_part]);
            %sorting by RF size (to plot small ones on top of large ones)
            for u=group_units(sort_by_size_index)
                if RFsizes(group_units==u)==0
                    continue
                end
                hold on;
                Parameters_temp=[RFparameters(group_units==u).(RFparameters(group_units==u).bestfit)];
                for ellipsn=1:2
                    switch RFparameters(group_units==u).bestfit
                        case 'gaussian1'
                            if ellipsn==2
                                continue
                            end
                            sx=Parameters_temp.sx;
                            sy=Parameters_temp.sy;
                            xmax=Parameters_temp.xmax;
                            ymax=Parameters_temp.ymax;
                            zmax=Parameters_temp.zmax;
                        case {'gaussian2','gaussian15'}
                            sx=[Parameters_temp.(['sx' num2str(ellipsn)])];
                            sy=[Parameters_temp.(['sy' num2str(ellipsn)])];
                            xmax=Parameters_temp.(['xmax' num2str(ellipsn)]);
                            ymax=Parameters_temp.(['ymax' num2str(ellipsn)]);
                            zmax=Parameters_temp.(['zmax' num2str(ellipsn)]);
                    end
                    RFsize              =2*sqrt(sx.*sy);
                    center              =[xmax ymax];
                    Allmonkeys={'Linus','Curius','Cornelius'};
                    current_monkey=Allmonkeys{cellfun(@(x) any(strfind(x,population(u).unit_ID(1:3))),Allmonkeys)};
                    
                    if sign(zmax)==1
                        col='r';
                    elseif sign(zmax)==-1
                        col='b';
                    else
                        col='k';
                    end
                    monkey_marker=keys.(current_monkey).marker;
                    RF_size_factor=0.1;
                    Radius=RF_size_factor*RFsize;
                    if strcmp(monkey_marker,'o')
                        ellipse_x = circle_x.*Radius;
                        ellipse_y = circle_y.*Radius;
                        line(ellipse_x+center(1),ellipse_y+center(2),'color',col,'linewidth',4);
                    elseif strcmp(monkey_marker,'s')
                        square_x = [-1,-1,1,1,-1].*Radius*sqrt(pi)/2;
                        square_y = [-1,1,1,-1,-1].*Radius*sqrt(pi)/2;
                        line(square_x+center(1),square_y+center(2),'color',col,'linewidth',4);
                    end
                end
            end
            
            axis equal
            set(gca,'Ydir','normal','xlim',[min(fitsettings.xout) max(fitsettings.xout)],'ylim',[min(fitsettings.yout) max(fitsettings.yout)]);
            ph_title_and_save(f_handle,  [filename plot_title_part],[fig_title plot_title_part],keys);
            
            
            %% RF sizes
            plot_title_part         = ['=' unique_group_values{g} ' con' num2str(c) ' gaussian RF  sizes in ' keys.PO.epoch_RF];
            f_handle              = figure('units','normalized','outerposition',[0 0 1 1],'name',[fig_title plot_title_part]);
            %sorting by RF size (to plot small ones on top of large ones)
            RFsizes_valid=RFsizes(RFsizes~=0);
            minsiz=floor(min(RFsizes_valid));
            maxsiz=ceil(max(RFsizes_valid));
            sizebins=minsiz:(maxsiz-minsiz)/20:maxsiz;
            FRhist=hist(RFsizes_valid,sizebins);
            hold on
            if any(FRhist~=0)
                bar(sizebins,FRhist);
                text(sizebins(2), max(FRhist)-max(FRhist/10),['u=' num2str(nanmean(RFsizes_valid)) ', med=' num2str(nanmedian(RFsizes_valid)) ', std=' num2str(nanstd(RFsizes_valid))]);
            end
            ph_title_and_save(f_handle,  [filename plot_title_part],[fig_title plot_title_part],keys);
            save([keys.path_to_save filename plot_title_part '.mat'],'RFsizes','RFparameters');
        end
    end
end

    function plot_PSTH_no_empties
        for eff=u_effectors %% one figure for ech effector, somehting probably does not work here, because (!!) suplots in each figure
            ef=find(u_effectors==eff);
            keys=ph_get_epoch_keys(keys,typ,eff,sum(type_effectors(:,1)==typ)>1);
            [~, type_effector_short] = MPA_get_type_effector_name(typ,eff);
            plot_title              = [fig_title plot_title_part ', ' type_effector_short ];
            PSTH_summary_handle(ef)     = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            for g=1:numel(unique_group_values)
                %%reducing to only units that have at least one condition
                unitidx=ismember(complete_unit_list,tuning_per_unit_table(ismember(group_values,unique_group_values(g)),idx_unitID));
                current_window=vertcat(current.window);
                current_units=vertcat(current_window(:,1).unit);
               
                empty_conditions_and_units=arrayfun(@(x) isempty(x.average_spike_density) || all(isnan( x.average_spike_density)),current_units) ; %nans not needed hopefully 
                empty_conditions_and_units(:,end+1:numel(unitidx))=true; % from the last valid to the end
                condition_empty=all(empty_conditions_and_units,2);
                
                if any(strfind(plot_title_part,'per position'))
                    [positions, ~,pos_sub_idx]=unique(vertcat(current.position),'rows');
                    [~, ~,fix_sub_idx]=unique(vertcat(current.fixation),'rows');
                    [subplot_pos, columns_to_loop, rows_to_loop]= DAG_dynamic_positions({positions});
                    conditions_to_loop=find([current.effector]'==eff & ~condition_empty);
                    sp_per_effector=max(subplot_pos);
                else
                    conditions_to_loop=find([current.effector]'==eff & ~condition_empty);
                    [~,columns_to_loop]=ismember(column_indexes,unique(column_indexes(conditions_to_loop)));
                    sp_per_effector=max(columns_to_loop)*numel(unique_group_values);
                end
                n_units_title_part=repmat({''},sp_per_effector*numel(u_effectors),1);
                spl=zeros(sp_per_effector*numel(u_effectors),1);
                
                legend_label_indexes=[];
                legend_line_handles=[];
                
                ensu_units{ef,1}=[];ensu_units{ef,2}=[];ensu_units{ef,3}=[];
                for c=conditions_to_loop(:)'
                    group_condition_units=find(all(unitidx,2) & units_valid & ~empty_conditions_and_units(c,:)')';
                    
                    %% this seems rather complicated for retrieving color information, but i
                    if strcmp(plot_title_part,' PSTHs')
                        column=columns_to_loop(c);
                        spn=(g-1)*max(columns_to_loop)+column+(ef-1)*sp_per_effector;
                        sph(spn)=subplot(numel(unique_group_values),max(columns_to_loop),spn-(ef-1)*sp_per_effector);
                        col=(conditions_PSTH(c,1)+1)*6 + (conditions_PSTH(c,3))*2 + (conditions_PSTH(c,4)+1) + (conditions_PSTH(c,5)*18);
                        current_color=keys.line_colors(col,:);
                        spl(spn)=spl(spn)+1;
                    elseif any(strfind(plot_title_part,'pref'))
                        column=columns_to_loop(c);
                        spn=(g-1)*max(columns_to_loop)+column+(ef-1)*sp_per_effector;
                        sph(spn)=subplot(numel(unique_group_values),max(columns_to_loop),spn-(ef-1)*sp_per_effector);
                        col=(conditions_PSTH(c,1))*6 + (conditions_PSTH(c,3))*2 + (conditions_PSTH(c,4)+1) + (conditions_PSTH(c,5)*12);
                        current_color=keys.pref_colors(col,:);
                        spl(spn)=spl(spn)+1;
                    elseif any(strfind(plot_title_part,'per position'))
                        %% still need to fix colors if different conditions are present, AND fix overwriting of different groups.... (how about legends?)
                        column=c; %irrelevant
                        spn=subplot_pos(pos_sub_idx(c))+(ef-1)*sp_per_effector;
                        sph(spn)=subplot(rows_to_loop,max(columns_to_loop),spn-(ef-1)*sp_per_effector);
                        col=fix_sub_idx(c);
                        current_color=keys.colors.fix_offset(col,:);
                        signs=[current(c).sign.unit];
                        spl(spn)=spl(spn)+1;
                    end
                    spf(spn)=ef;
                    %specific additional loop for enhancement and suppression
                    ensu_loops=1;
                    if any(strfind(plot_title_part,'ensu'))
                        ensu_loops=2;
                        ensucolors={[0 0 1]/spl(spn),[1 0 0]/spl(spn)};
                    end
                    
                    for ensu=1:ensu_loops
                        if any(strfind(plot_title_part,'ensu'))
                            current_color=ensucolors{ensu};
                            units=intersect(group_condition_units,find(signs==(ensu-1.5)*2));
                            if any(ismember(find(signs==1),group_condition_units)) && any(ismember(find(signs==-1),group_condition_units))
                                ensu_units{ef,3}=union(ensu_units{ef,3},units);
                                ensu_units{ef,ensu}=union(ensu_units{ef,ensu},units);
                            end
                        else
                            units=group_condition_units;
                        end
                        hold on
                        legend_label_indexes=[legend_label_indexes col];
                        n_units_title_part{spn}=[n_units_title_part{spn} '\color[rgb]{' num2str(current_color) '}' num2str(numel(units)) ' '];
                        if  numel(units)==0
                            continue;
                        end
                        state_shift=0;
                        for w=1:size(keys.PSTH_WINDOWS,1)
                            t_before_state=keys.PSTH_WINDOWS{w,3};
                            t_after_state=keys.PSTH_WINDOWS{w,4};
                            bins=t_before_state:keys.PSTH_binwidth:t_after_state;
                            bins=bins+state_shift-t_before_state;
                            props={'color',current_color,'linewidth',1};
                            errorbarhandle=shadedErrorBar(bins,nanmean(vertcat(current(c).window(w).unit(units).average_spike_density),1),...
                                sterr(vertcat(current(c).window(w).unit(units).average_spike_density),1),props,1); %% STERR!!!!
                            state_shift=state_shift+t_after_state-t_before_state+0.1;
                        end
                        legend_line_handles=[legend_line_handles errorbarhandle.mainLine];
                        title(sprintf('%s = %s, N =%s ',keys.PO.group_parameter,unique_group_values{g},n_units_title_part{spn}),'interpreter','tex'); %%
                    end
                end
                y_lim(spn,:)=get(gca,'ylim');
            end
            
            if keys.plot.population_PSTH_legends
                legend(legend_line_handles,legend_labels(legend_label_indexes),'location','southwest');
            end
        end
        
        %% subplot appearance, and tuning lines
        spf(sph==0)=[];
        sph(sph==0)=[];
        ylimmax=max(max(y_lim));
        ylimmin=min(min(y_lim));
        y_lim=[ylimmin ylimmax];
        if ~isempty(keys.PO.y_lim)
            y_lim= keys.PO.y_lim;
        end
        
        for spn=1:numel(sph)
            
                    ef=spf(spn);
                    set(0, 'CurrentFigure', PSTH_summary_handle(ef));
            subplot(sph(spn));
            hold on
            
            %% completed? choices? hands?
            %                 tr=[all_trialz.type]==typ & [all_trialz.effector]==eff & ismember([all_trialz.completed],keys.cal.completed) &...
            %                     ismember([all_trialz.completed],keys.cal.completed) & ismember([all_trialz.reach_hand],keys.tt.hands) & ismember([all_trialz.choice],u_choice);
            tr=[all_trialz.type]==typ & ismember([all_trialz.completed],keys.cal.completed) &...
                ismember([all_trialz.completed],keys.cal.completed) & ismember([all_trialz.reach_hand],u_hands) & ismember([all_trialz.choice],u_choice);
            ph_PSTH_background(all_trialz(tr),y_lim,y_lim,y_lim,keys,keys.PO.fontsize_factor)
            
        end
        for eff=u_effectors
            ef=find(u_effectors==eff);
            [~, type_effector_short] = MPA_get_type_effector_name(typ,eff);
            plot_title              = [fig_title plot_title_part ', ' type_effector_short ];
            ph_title_and_save(PSTH_summary_handle(ef),  [filename plot_title_part ', ' type_effector_short],plot_title,keys)
        end
        
        %% ensu summary figure
        if any(strfind(plot_title_part,'ensu'))
            
            for eff=u_effectors
                ef=find(u_effectors==eff);
                [~, type_effector_short] = MPA_get_type_effector_name(typ,eff);
                plot_title              = [fig_title plot_title_part ', ' type_effector_short ];
                PSTH_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
                for g=1:numel(unique_group_values)
                    unitidx=ismember(complete_unit_list,tuning_per_unit_table(ismember(group_values,unique_group_values(g)),idx_unitID));
                    group_units=find(all(unitidx,2) & units_valid)';
                    %%reducing to only units that have at least one condition
                    current_window=vertcat(current.window);
                    current_units=vertcat(current_window(:,1).unit);
                    empty_conditions_and_units=arrayfun(@(x) isempty(x.average_spike_density),current_units);
                    empty_conditions_and_units(:,end+1:numel(unitidx))=true;
                    condition_empty=all(empty_conditions_and_units,2);
                    
                    if any(strfind(plot_title_part,'per position'))
                        conditions_to_loop=find([current.effector]'==eff & ~condition_empty);
                    else
                        conditions_to_loop=find([current.effector]'==eff & ~condition_empty);
                    end
                    
                    posneg_idx=0;
                    for c=conditions_to_loop(:)'
                        signs=[current(c).sign.unit];
                        units_pos=intersect(group_units,find(signs==1));
                        units_neg=intersect(group_units,find(signs==-1));
                        if  numel(units_pos)==0 && numel(units_pos)==0 %% not sure if this is enough
                            continue;
                        end
                        posneg_idx=posneg_idx+1;
                        for w=1:size(keys.PSTH_WINDOWS,1)
                            if ~isempty(units_pos) && ~isempty(units_neg)
                                positive(w).psth(posneg_idx,:)=nanmedian(vertcat(current(c).window(w).unit(units_pos).average_spike_density),1);
                                negative(w).psth(posneg_idx,:)=nanmedian(vertcat(current(c).window(w).unit(units_neg).average_spike_density),1);
                                posneg(w).psth(posneg_idx,:)=positive(w).psth(posneg_idx,:)-negative(w).psth(posneg_idx,:);
                            end
                        end
                    end
                    
                    subplot(1,numel(unique_group_values),g);
                    hold on
                    
                    state_shift=0;
                    for w=1:size(keys.PSTH_WINDOWS,1)
                        t_before_state=keys.PSTH_WINDOWS{w,3};
                        t_after_state=keys.PSTH_WINDOWS{w,4};
                        bins=t_before_state:keys.PSTH_binwidth:t_after_state;
                        bins=bins+state_shift-t_before_state;
                        props={'color','r','linewidth',1};
                        if exist('positive','var') && exist('negative','var')
                            errorbarhandle=shadedErrorBar(bins,nanmean(vertcat(positive(w).psth),1),sterr(vertcat(positive(w).psth),1),props,1);
                            props={'color','b','linewidth',1};
                            errorbarhandle=shadedErrorBar(bins,nanmean(vertcat(negative(w).psth),1),sterr(vertcat(negative(w).psth),1),props,1);
                            props={'color','g','linewidth',1};
                            errorbarhandle=shadedErrorBar(bins,nanmean(vertcat(posneg(w).psth),1),sterr(vertcat(posneg(w).psth),1),props,1);
                        end
                        state_shift=state_shift+t_after_state-t_before_state+0.1;
                    end
                    title(['N contributing = ' num2str(numel(ensu_units{ef,2})) ' POS, ' num2str(numel(ensu_units{ef,1})) ' NEG, ' num2str(numel(ensu_units{ef,3})) ' POS-NEG']);
                    %% this part should eventually go in extra loop for same dimensions in each subplot
                    y_lim=get(gca,'ylim');
                    tr=[all_trialz.type]==typ & [all_trialz.effector]==eff & ismember([all_trialz.completed],keys.cal.completed) &...
                        ismember([all_trialz.completed],keys.cal.completed) & ismember([all_trialz.reach_hand],u_hands) & ismember([all_trialz.choice],u_choice);
                    ph_PSTH_background(all_trialz(tr),y_lim,y_lim,y_lim,keys,1)
                end
                ph_title_and_save(PSTH_summary_handle,  [filename plot_title_part ' summary, ' type_effector_short ],plot_title,keys)
            end
            
        end
    end

end

function h=ph_compare_by_bin(in,c1,c2,units)

% to compare only the same cells (!)
c12=[c1(:);c2(:)];
for cc=c12'
    units_c=find(~cellfun(@isempty,{in(cc).unit.average_spike_density}));
    units=intersect(units,units_c);
end

if ~isempty(units)
    C1=[];
    C2=[];
    for cc=c1(:)'
        C1=[C1;vertcat(in(cc).unit(units).average_spike_density)];
    end
    for cc=c2(:)'
        C2=[C2;vertcat(in(cc).unit(units).average_spike_density)];
    end
    h=ttest2(C1,C2,0.05,'both','equal',1); h(h==0)=NaN;
else
    h=NaN; %% size = number of bins ?
end
end

% 
% function plot_PSTH
% for eff=u_effectors %% one figure for ech effector
%     keys=ph_get_epoch_keys(keys,typ,eff,sum(type_effectors(:,1)==typ)>1);
%     [~, type_effector_short] = MPA_get_type_effector_name(typ,eff);
%     plot_title              = [fig_title type_effector_short plot_title_part];
%     PSTH_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
%     for g=1:numel(unique_group_values)
%         unitidx=ismember(complete_unit_list,tuning_per_unit_table(ismember(group_values,unique_group_values(g)) & tya_existing{t},idx_unitID));
%         group_units=find(all(unitidx,2) & units_valid)';
%         %%reducing to only units that have each condition
%         current_window=vertcat(current.window);
%         current_units=vertcat(current_window(:,1).unit);
%         empty_conditions_and_units=arrayfun(@(x) isempty(x.average_spike_density),current_units);
%         condition_empty=all(empty_conditions_and_units,2);
%         % this is to make sure cells are present in all conditions...
%         
%         if any(strfind(plot_title_part,'per position'))
%             
%             
%             [positions, ~,pos_sub_idx]=unique(vertcat(current.position),'rows');
%             [~, ~,fix_sub_idx]=unique(vertcat(current.fixation),'rows');
%             [subplot_pos, columns_to_loop, rows_to_loop]= DAG_dynamic_positions({positions});
%             conditions_to_loop=find([current.effector]'==eff & ~condition_empty);
%         else
%             
%             conditions_to_loop=find([current.effector]'==eff & ~condition_empty);
%             %conditions_to_loop=find(conditions_PSTH(:,2)==eff & ~condition_empty);
%             [~,columns_to_loop]=ismember(column_indexes,unique(column_indexes(conditions_to_loop)));
%             group_units=intersect(group_units,find(all(~empty_conditions_and_units(~condition_empty,:),1)));
%             
%             
%         end
%         
%         n_units=[];
%         legend_label_indexes=[];
%         legend_line_handles=[];
%         for c=conditions_to_loop(:)'
%             
%             %% this seems rather complicated for retrieving color information
%             if strcmp(plot_title_part,' PSTHs')
%                 column=columns_to_loop(c);
%                 spn=(g-1)*max(columns_to_loop)+column;
%                 sph(spn)=subplot(numel(unique_group_values),max(columns_to_loop),spn);
%                 col=(conditions_PSTH(c,1)+1)*6 + (conditions_PSTH(c,3))*2 + (conditions_PSTH(c,4)+1) + (conditions_PSTH(c,5)*18);
%                 current_color=keys.line_colors(col,:);
%             elseif any(strfind(plot_title_part,'pref'))
%                 column=columns_to_loop(c);
%                 spn=(g-1)*max(columns_to_loop)+column;
%                 sph(spn)=subplot(numel(unique_group_values),max(columns_to_loop),spn);
%                 col=(conditions_PSTH(c,1))*6 + (conditions_PSTH(c,3))*2 + (conditions_PSTH(c,4)+1) + (conditions_PSTH(c,5)*12);
%                 current_color=keys.pref_colors(col,:);
%             elseif any(strfind(plot_title_part,'per position'))
%                 %% still need to fix colors if different conditions are present, AND fix overwriting of different groups....
%                 column=c;
%                 spn=subplot_pos(pos_sub_idx(c));
%                 sph(spn)=subplot(rows_to_loop,max(columns_to_loop),spn);
%                 col=fix_sub_idx(c);
%                 current_color=keys.colors.fix_offset(col,:);
%             end
%             
%             hold on
%             legend_label_indexes=[legend_label_indexes col];
%             
%             units=intersect(group_units,find(arrayfun(@(x) ~isempty(x.average_spike_density),current(c).window(1).unit))) ; %why here again?
%             n_units=[n_units numel(units)];
%             if  numel(units)==0
%                 continue;
%             end
%             state_shift=0;
%             for w=1:size(keys.PSTH_WINDOWS,1)
%                 t_before_state=keys.PSTH_WINDOWS{w,3};
%                 t_after_state=keys.PSTH_WINDOWS{w,4};
%                 bins=t_before_state:keys.PSTH_binwidth:t_after_state;
%                 bins=bins+state_shift-t_before_state;
%                 props={'color',current_color,'linewidth',1};
%                 errorbarhandle=shadedErrorBar(bins,nanmean(vertcat(current(c).window(w).unit(units).average_spike_density),1),...
%                     sterr(vertcat(current(c).window(w).unit(units).average_spike_density),1),props,1); %% STERR!!!!
%                 state_shift=state_shift+t_after_state-t_before_state+0.1;
%             end
%             legend_line_handles=[legend_line_handles errorbarhandle.mainLine];
%         end
%         title(sprintf('%s = %s, N (NH/CH/IH) =%s ',keys.PO.group_parameter,unique_group_values{g},num2str(n_units)),'interpreter','none');
%         y_lim(spn,:)=get(gca,'ylim');
%     end
%     
%     if keys.plot.population_PSTH_legends
%         legend(legend_line_handles,legend_labels(legend_label_indexes),'location','southwest');
%     end
% end
% 
% %% subplot appearance, and tuning lines
% sph(sph==0)=[];
% ylimmax=max(max(y_lim));
% ylimmin=min(min(y_lim));
% y_lim=[ylimmin ylimmax];
% if ~isempty(keys.PO.y_lim)
%     y_lim= keys.PO.y_lim;
% end
% for eff=u_effectors
%     [~, type_effector_short] = MPA_get_type_effector_name(typ,eff);
%     plot_title              = [fig_title type_effector_short plot_title_part];
%     %                 unitidx=ismember(complete_unit_list,tuning_per_unit_table(ismember(group_values,unique_group_values(g)) & tya_existing{t},idx_unitID));
%     %                 group_units=find(all(unitidx,2) & units_valid)';
%     %                 %%reducing to only units that have each condition
%     %                 current_window=vertcat(current.window);
%     %                 current_units=vertcat(current_window(:,1).unit);
%     %                 empty_conditions_and_units=arrayfun(@(x) isempty(x.average_spike_density),current_units);
%     %                 condition_empty=all(empty_conditions_and_units,2);
%     %                 conditions_to_loop=find(conditions_PSTH(:,2)==eff & ~condition_empty);
%     %                 [~,columns_to_loop]=ismember(column_indexes,unique(column_indexes(conditions_to_loop)));
%     for spn=1:numel(sph)
%         
%         subplot(sph(spn));
%         hold on
%         
%         %% completed? choices? hands?
%         tr=[all_trialz.type]==typ & [all_trialz.effector]==eff & ismember([all_trialz.completed],keys.cal.completed) &...
%             ismember([all_trialz.completed],keys.cal.completed) & ismember([all_trialz.reach_hand],u_hands) & ismember([all_trialz.choice],u_choice);
%         ph_PSTH_background(all_trialz(tr),y_lim,y_lim,y_lim,keys,1)
%         
%         %                 %% peak times
%         %                                 Contra_handles=[];
%         %                                 Contra_indexes=[];
%         %                                 Ipsi_handles=[];
%         %                                 Ipsi_indexes=[];
%         %                                 for l=lines
%         %                                     Contra=vertcat(conditions(c).case(a).group(g).contra.line(l).unit.per_state);
%         %                                     Ipsi=vertcat(conditions(c).case(a).group(g).ipsi.line(l).unit.per_state);
%         %                                     mean_c=(nanmean(vertcat(Contra(:,s).peak_bin),1)-1)*keys.PSTH_binwidth;
%         %                                     mean_i=(nanmean(vertcat(Ipsi(:,s).peak_bin),1)-1)*keys.PSTH_binwidth;
%         %                                     sterr_c=(sterr(vertcat(Contra(:,s).peak_bin),1))*keys.PSTH_binwidth;
%         %                                     sterr_i=(sterr(vertcat(Ipsi(:,s).peak_bin),1))*keys.PSTH_binwidth;
%         %                                     if ~isempty(mean_c)
%         %                                     CL=line([mean_c mean_c]+state_shift,y_lim,'color',keys.hnd_choice_colors_C(l,:));
%         %                                     line([mean_c+sterr_c mean_c+sterr_c]+state_shift,y_lim,'color',keys.hnd_choice_colors_C(l,:),'linestyle',':');
%         %                                     line([mean_c-sterr_c mean_c-sterr_c]+state_shift,y_lim,'color',keys.hnd_choice_colors_C(l,:),'linestyle',':');
%         %                                     Contra_handles=[Contra_handles CL];
%         %                                     Contra_indexes=[Contra_indexes l];
%         %                                     end
%         %                                     if ~isempty(mean_i)
%         %                                     IL=line([mean_i mean_i]+state_shift,y_lim,'color',keys.hnd_choice_colors_I(l,:));
%         %                                     line([mean_i+sterr_i mean_i+sterr_i]+state_shift,y_lim,'color',keys.hnd_choice_colors_I(l,:),'linestyle',':');
%         %                                     line([mean_i-sterr_i mean_i-sterr_i]+state_shift,y_lim,'color',keys.hnd_choice_colors_I(l,:),'linestyle',':');
%         %                                     Ipsi_handles=[Ipsi_handles IL];
%         %                                     Ipsi_indexes=[Ipsi_indexes l];
%         %                                     end
%         %                                 end
%         
%         
%         
%         %                     %% space tuning
%         %                     for l=lines
%         %                         l_height=y_lim(2)-diff(y_lim)*(l*0.01+0.1);
%         %                         if ~isempty(conditions(c).case(a).group(g).space_tuning.line(l).per_state(s).h)
%         %                             plot(bins,conditions(c).case(a).group(g).space_tuning.line(l).per_state(s).h*l_height,...
%         %                                 'linewidth',2,'color',keys.hnd_choice_colors(l,:))
%         %                         end
%         %                     end
%         %
%         %                     %% hand tuning line
%         %                     for ic=0:1
%         %                         CH_idx=find(choices_hands(2,:)==1 & choices_hands(1,:)==ic );
%         %                         IH_idx=find(choices_hands(2,:)==2 & choices_hands(1,:)==ic );
%         %                         if ~ismember(CH_idx,lines) || ~ismember(IH_idx,lines)
%         %                             continue
%         %                         end
%         %                         for side=1:2
%         %                             l_height=y_lim(2)-diff(y_lim)*((2*side+ic+5)*0.01+0.1);
%         %                             if ~isempty(conditions(c).case(a).group(g).hand_tuning(ic+1).(CIFN{side}).per_state(s).h)
%         %                                 plot(bins,conditions(c).case(a).group(g).hand_tuning(ic+1).(CIFN{side}).per_state(s).h*l_height,...
%         %                                     'linewidth',2,'color',keys.space_colors(side,:))
%         %                             end
%         %                         end
%         %                     end
%         %                     state_shift=state_shift+t_after_state-t_before_state+0.1;
%     end
%     ph_title_and_save(PSTH_summary_handle,  [filename plot_title_part],plot_title,keys)
% end
% end


%                 %% peak times
%                                 Contra_handles=[];
%                                 Contra_indexes=[];
%                                 Ipsi_handles=[];
%                                 Ipsi_indexes=[];
%                                 for l=lines
%                                     Contra=vertcat(conditions(c).case(a).group(g).contra.line(l).unit.per_state);
%                                     Ipsi=vertcat(conditions(c).case(a).group(g).ipsi.line(l).unit.per_state);
%                                     mean_c=(nanmean(vertcat(Contra(:,s).peak_bin),1)-1)*keys.PSTH_binwidth;
%                                     mean_i=(nanmean(vertcat(Ipsi(:,s).peak_bin),1)-1)*keys.PSTH_binwidth;
%                                     sterr_c=(sterr(vertcat(Contra(:,s).peak_bin),1))*keys.PSTH_binwidth;
%                                     sterr_i=(sterr(vertcat(Ipsi(:,s).peak_bin),1))*keys.PSTH_binwidth;
%                                     if ~isempty(mean_c)
%                                     CL=line([mean_c mean_c]+state_shift,y_lim,'color',keys.hnd_choice_colors_C(l,:));
%                                     line([mean_c+sterr_c mean_c+sterr_c]+state_shift,y_lim,'color',keys.hnd_choice_colors_C(l,:),'linestyle',':');
%                                     line([mean_c-sterr_c mean_c-sterr_c]+state_shift,y_lim,'color',keys.hnd_choice_colors_C(l,:),'linestyle',':');
%                                     Contra_handles=[Contra_handles CL];
%                                     Contra_indexes=[Contra_indexes l];
%                                     end
%                                     if ~isempty(mean_i)
%                                     IL=line([mean_i mean_i]+state_shift,y_lim,'color',keys.hnd_choice_colors_I(l,:));
%                                     line([mean_i+sterr_i mean_i+sterr_i]+state_shift,y_lim,'color',keys.hnd_choice_colors_I(l,:),'linestyle',':');
%                                     line([mean_i-sterr_i mean_i-sterr_i]+state_shift,y_lim,'color',keys.hnd_choice_colors_I(l,:),'linestyle',':');
%                                     Ipsi_handles=[Ipsi_handles IL];
%                                     Ipsi_indexes=[Ipsi_indexes l];
%                                     end
%                                 end



%                     %% space tuning
%                     for l=lines
%                         l_height=y_lim(2)-diff(y_lim)*(l*0.01+0.1);
%                         if ~isempty(conditions(c).case(a).group(g).space_tuning.line(l).per_state(s).h)
%                             plot(bins,conditions(c).case(a).group(g).space_tuning.line(l).per_state(s).h*l_height,...
%                                 'linewidth',2,'color',keys.hnd_choice_colors(l,:))
%                         end
%                     end
%
%                     %% hand tuning line
%                     for ic=0:1
%                         CH_idx=find(choices_hands(2,:)==1 & choices_hands(1,:)==ic );
%                         IH_idx=find(choices_hands(2,:)==2 & choices_hands(1,:)==ic );
%                         if ~ismember(CH_idx,lines) || ~ismember(IH_idx,lines)
%                             continue
%                         end
%                         for side=1:2
%                             l_height=y_lim(2)-diff(y_lim)*((2*side+ic+5)*0.01+0.1);
%                             if ~isempty(conditions(c).case(a).group(g).hand_tuning(ic+1).(CIFN{side}).per_state(s).h)
%                                 plot(bins,conditions(c).case(a).group(g).hand_tuning(ic+1).(CIFN{side}).per_state(s).h*l_height,...
%                                     'linewidth',2,'color',keys.space_colors(side,:))
%                             end
%                         end
%                     end
%                     state_shift=state_shift+t_after_state-t_before_state+0.1;
