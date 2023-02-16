function [pop]=ph_arrange_positions_and_plots(keys,o,pop_in)
% this nasty piece of code defines positions, single cell plotting appearence (batching in figures/subplots/lines and color), as well as minimum number of trials per condition,
% all in one defined by keys.arrangement
% con_for_figure:           which trials to place on the same figure
% con_for_column:           which trials to place in the same subplot column
% con_for_row:              which trials to place in the same subplot row
% con_for_line:             which trials to batch for one PSTH line
% con_for_trial_crit:       which batching to use for calculating minimum trials per condition (apart from positions and conditions which are already split in figures (con_for_figure)
% fig_title:                part of the figure name indicating by which parameter they are different
% val_for_figure:           The values for those parameters (particularly relevant to create unique PDF names)
% sub_title:                subplot title name for per position plots indicating by what (position) parameter subplots are selected
% val_for_sub_assignment:   The values of this (position) parameter per subplot
% position_indexes:         Unique position indexes used for further analysis (very important to solve precision issues for small changes in position, batching close postiions together)
% val_for_pos_assignment:   The position values associated to the respective indexes
% hemifield_indexes:        Unique hemifield indexes used for further analysis (very important if hemifield is defined relative to fixation for example)
% fixation_per_trial:       Fixation values used for further processing
% out_analysis.PSTH_perpos_colors:  colors used in per position plots
% out_analysis.PSTH_summary_colors: colors used in type summary plots
% out_analysis.line_labels: lables for the different lines in the summary plot legend
%
%% taking over per unit informaiton from original pop if applicable
if nargin>2
    pop=rmfield(pop_in,'trial');
end

%% DEFINITION OF CONDITION INDICES TO PLOT, CURRENTLY TARGET LOCATION, FIXATION LOCATION, MOVEMENT VECTORS, CHOICE, HANDS
[~, displacement_types] = center_displacement_working(o,keys);
all_val=displacement_types(:,[1,2,5,6]);
fix_val=displacement_types(:,1:2);
mov_val=displacement_types(:,3:4);
tar_val=displacement_types(:,5:6);
cue_val=displacement_types(:,7:8);
stm_val=displacement_types(:,9:10);

all_idx=displacement_types(:,11);
fix_idx=displacement_types(:,12);
mov_idx=displacement_types(:,13);
tar_idx=displacement_types(:,14);
cue_idx=displacement_types(:,15);
stm_idx=displacement_types(:,16);
non_idx=ones(size(displacement_types,1),1);

[~,u_all_idx_idx]=unique(all_idx);
[~,u_fix_idx_idx]=unique(fix_idx);
[~,u_tar_idx_idx]=unique(tar_idx);
[~,u_cue_idx_idx]=unique(cue_idx);
[~,u_mov_idx_idx]=unique(mov_idx);
[~,u_stm_idx_idx]=unique(stm_idx);

choices                 =[o.choice]';
hands                   =[o.reach_hand]';
cueshape                =[o.cue_shape]';
success                 =[o.success]';
effectors               =[o.effector]';
perturbations_orig      =[o.perturbation]';
difficulty              = [o.difficulty];
stimuli_in_2hemifields  = [o.stimuli_in_2hemifields];
n_distractors           = [o.n_distractors];
n_nondistractors        = [o.n_nondistractors];
StimulusType            = [o.stimulustype];

%% perturbation (painful, because its either in groups - for analysis, or by block (actuallz original perturbation value from table) - for single cell plotting
perturbations=zeros(size(perturbations_orig));
unique_perturbations_temp=unique(perturbations_orig);
for PT=1:numel(unique_perturbations_temp) %% this doesnt work properly if there are only 2/3 sets though
    perturbations(perturbations_orig==unique_perturbations_temp(PT))=PT;
end
perturbation_group=nan(size(perturbations_orig));
for PT=1:numel(keys.cal.perturbation_groups)
    perturbation_group(ismember(perturbations_orig,keys.cal.perturbation_groups{PT}))  =PT;
end
blocks_temp   =[o.block]';
blocks=zeros(size(blocks_temp));
blo_values_o=unique(blocks_temp);
for PT=1:numel(blo_values_o) %% this doesnt work properly if there are only 2/3 sets though
    blocks(blocks_temp==blo_values_o(PT))=PT;
end

%% indexes for later assignments
[cho_values,~,cho_idx]      =unique(choices);
[hnd_values,~,hnd_idx]      =unique(hands);
[shp_values,~,shp_idx]      =unique(cueshape);
[suc_values,~,suc_idx]      =unique(success);
[ptb_values,~,ptb_idx]      =unique(perturbations);
[blo_values,~,blo_idx]      =unique(blocks);
[eff_values,~,eff_idx]      =unique(effectors);
[diff_values,~,diff_idx]    =unique(difficulty);
%[diff_values,~,diff_idx]    =unique(difficulty+3*success'); % i believe
%this separation is not needed any more
[Styp_values,~,Styp_idx]    =unique(StimulusType');
[hnd_cho_values,~,hch_idx]  =unique([hands choices],'rows');
[hnd_ptb_values,~,hpt_idx]  =unique([hands perturbations],'rows');
[hnd_blo_values,~,hbl_idx]  =unique([hands blocks],'rows');
[hnd_eff_values,~,hef_idx]  =unique([hands effectors],'rows');
[cho_ptb_values,~,chp_idx]  =unique([choices perturbations],'rows');

%% defaults
con_for_figure          = non_idx;
con_for_column          = non_idx;
con_for_row             = eff_idx;
fig_title               = '';
sub_title               = 'movement vector ';
val_for_figure          = {[]};
val_for_sub_assignment  = mov_val(u_mov_idx_idx,:);
val_for_pos_assignment  = mov_val(u_mov_idx_idx,:);
position_indexes        = mov_idx;
fixation_indexes        = fix_idx;
hemifield_indexes       = (real([o.tar_pos]' - [o.fix_pos]')>0)+1;
fixation_per_trial      = fix_val;


%% specifics
switch keys.arrangement
    
    case 'StimulusType_Difficulty_Position'
        fig_title               = 'StimType';
        con_for_figure          = Styp_idx;
        val_for_figure          = num2cell(Styp_values);
        val_for_sub_assignment  = tar_val(u_tar_idx_idx,:);
        val_for_pos_assignment  = tar_val(u_tar_idx_idx,:);
        position_indexes        = tar_idx;
        sub_title               = 'stimulus position';
        
    case 'StimulusType_Difficulty_Position_Successful'
        [diff_values,~,diff_idx]        = unique(difficulty);
        fig_title               = 'StimuType_Diff_Pos_Suc';
        con_for_figure          = Styp_idx;
        val_for_figure          = num2cell(Styp_values);
        position_indexes        = stm_idx;
        val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
        val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
        sub_title               = 'stimulus position';
        
    case 'DoubleSameTargets_Position'
        [diff_values,~,diff_idx]      =unique(difficulty+3*success');
        fig_title               = 'StimulusType';
        con_for_figure          = Styp_idx;
        val_for_figure          = num2cell(Styp_values);
        position_indexes        = stm_idx;
        val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
        val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
        sub_title               = 'stimulus position';
    case 'StimType_Diff_Pos_ErVsCor'
        difficulty =  difficulty+3*success';
        %              difficulty =  difficulty(success ==1)';
        %              idx_difficulty = find(difficulty <3);
        %              difficulty(idx_difficulty)     =[];
        %              stm_idx(idx_difficulty)        =[];
        %              Styp_idx(idx_difficulty)       =[];
        
        [diff_values,~,diff_idx]        = unique(difficulty);
        fig_title               = 'StimType_Diff_Pos_ErVsCor';
        con_for_figure          = Styp_idx;
        val_for_figure          = num2cell(Styp_values);                
        position_indexes        = stm_idx;
        val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
        val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
        sub_title               = 'stimulus position';
        
    case 'StimulusType_Difficulty_Position' %% this is a guess that this one is supposed to be a different one     
        [diff_values,~,diff_idx]        = unique(difficulty);
        fig_title               = 'StimulusType_Difficulty_Position';
        con_for_figure          = Styp_idx;
        val_for_figure          = num2cell(Styp_values);
        val_for_sub_assignment  = tar_val(u_tar_idx_idx,:);
        val_for_pos_assignment  = tar_val(u_tar_idx_idx,:);
        position_indexes        = tar_idx;
        
        sub_title               = 'stimulus position';
        
    case 'StimulusType_Difficulty_Position_Successful'
        [diff_values,~,diff_idx]        = unique(difficulty);
        fig_title               = 'StimuType_Diff_Pos_Suc';
        con_for_figure          = Styp_idx;
        val_for_figure          = num2cell(Styp_values);
        position_indexes        = stm_idx;
        val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
        val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
        sub_title               = 'stimulus position';
        
    case 'DoubleSameTargets_Position'
        [diff_values,~,diff_idx]      =unique(difficulty+3*success');
        fig_title               = 'StimulusType';
        con_for_figure          = Styp_idx;
        val_for_figure          = num2cell(Styp_values);
        position_indexes        = stm_idx;
        val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
        val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
        sub_title               = 'stimulus position';
        % con_for_figure          = suc_idx;
        
    case 'StimType_Diff_Pos_ErVsCor'
        difficulty =  difficulty+3*success';
        %              difficulty =  difficulty(success ==1)';
        %              idx_difficulty = find(difficulty <3);
        %              difficulty(idx_difficulty)     =[];
        %              stm_idx(idx_difficulty)        =[];
        %              Styp_idx(idx_difficulty)       =[];
        
    case 'StimTyp_Diff_Pos_Suc'
        [diff_values,~,diff_idx]        = unique(difficulty);        
        fig_title               = 'StimTyp_Diff_Pos_Suc';
        con_for_figure          = Styp_idx;
        val_for_figure          = num2cell(Styp_values);        
        sub_title               = 'stimulus position';        
        if length(diff_values) == 3
            fig_title               = 'StimTyp_Diff_Pos_Suc';            
            position_indexes        = stm_idx;
            val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
            val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
        end
        
    case 'Sgl_Diff_Pos_Suc_SaccadeEpoch'
        [diff_values,~,diff_idx]        = unique(difficulty);
        fig_title               = 'StimTyp_Diff_Pos_Suc';
        con_for_figure          = Styp_idx;
        val_for_figure          = num2cell(Styp_values);
        val_for_sub_assignment  = tar_val(u_tar_idx_idx,:);
        val_for_pos_assignment  = tar_val(u_tar_idx_idx,:);
        position_indexes        = tar_idx;       
        sub_title               = 'stimulus position';
        
    case 'DisTar_Diff_Pos_Suc'
        [diff_values,~,diff_idx]        = unique(difficulty);
        sub_title               = 'stimulus position';
        val_for_sub_assignment  = tar_val(u_tar_idx_idx,:);
        val_for_pos_assignment  = tar_val(u_tar_idx_idx,:);
        position_indexes        = tar_idx;
        
    case 'SpatialCompetition_Targets'
        % left vs Right is defined where the saccade was made to
        sub_title               = 'stimulus position';
        val_for_sub_assignment  = tar_val(u_tar_idx_idx,:);
        val_for_pos_assignment  = tar_val(u_tar_idx_idx,:);
        position_indexes        = tar_idx;
        
    case 'SpatialCompetition_Distractor'
        [diff_values,~,diff_idx]        = unique(difficulty);
        fig_title               = 'DifficultyLevel';
        con_for_figure          = diff_idx';
        val_for_figure          = num2cell(diff_values);
        sub_title               = 'stimulus position';
        position_indexes        = stm_idx;
        val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
        val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
        
    case 'StimType_Diff_Pos_ErVsCor'
        difficulty =  difficulty+3*success';
        fig_title               = 'StimType_Diff_Pos_ErVsCor';
        con_for_figure          = Styp_idx;
        val_for_figure          = num2cell(Styp_values);
        position_indexes        = stm_idx;
        val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
        val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
        sub_title               = 'stimulus position';
        
    case 'StimulusType_Difficulty_Position_ErrorVsCorrect_lastVers'
        fig_title               = 'StimulusType_Difficulty_Position';
        con_for_figure          = Styp_idx';
        val_for_figure          = num2cell(Styp_values);
        position_indexes        = stm_idx;
        val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
        val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
        sub_title               = 'stimulus position';
        
    case 'success_in_cue'
        fig_title               = '';
        val_for_figure          = {shp_values};
        sub_title               = 'cue position';
        val_for_sub_assignment  = cue_val(u_cue_idx_idx,:);
        val_for_pos_assignment  = cue_val(u_cue_idx_idx,:);
        position_indexes        = cue_idx;
        
    case 'cue_position'
        con_for_figure          = suc_idx;
        fig_title               = 'success ';
        sub_title               = 'cue position';
        val_for_figure          = num2cell(suc_values);
        val_for_sub_assignment  = cue_val(u_cue_idx_idx,:);
        val_for_pos_assignment  = cue_val(u_cue_idx_idx,:);
        position_indexes        = cue_idx;
        
    case 'hand_choices'
        fig_title               = '';
        val_for_figure          = {hnd_cho_values};
        
    case 'options'
        con_for_figure          = hnd_idx;
        fig_title               = 'hand ';
        val_for_figure          = num2cell(hnd_values);
        
    case 'hands_inactivation'
        con_for_figure          = cho_idx;
        fig_title               = 'choice ';
        val_for_figure          = num2cell(cho_values);
        [~,~,con_for_column]    = unique([hands hemifield_indexes],'rows');
        
    case 'hands_inactivation_in_ch'
        con_for_figure          = eff_idx;
        fig_title               = 'choice_instructed_comp ';
        val_for_figure          = num2cell(eff_values);
        [~,~,con_for_column]    = unique([hands hemifield_indexes],'rows');
     
    case 'hands_in_ch'
        con_for_figure          = eff_idx;
        fig_title               = 'choice_instructed_comp ';
        val_for_figure          = num2cell(eff_values);
        [~,~,con_for_column]    = unique([hands hemifield_indexes],'rows');
  
    case 'hands'
        con_for_figure          = cho_idx;
        con_for_row             = ones(size(con_for_figure));
        fig_title               = 'choice ';
        sub_title               = 'movement vector ';
        val_for_figure          = num2cell(cho_values);
        con_for_column          = eff_idx;
        %Precision=5;
        
    case 'fixation'
        con_for_row             = eff_idx;
        sub_title               = 'fixation at ';
        val_for_sub_assignment  = fix_val(u_fix_idx_idx,:);
        val_for_pos_assignment  = fix_val(u_fix_idx_idx,:);
        position_indexes        = fix_idx;
        
    case 'movement vectors'
        sub_title               = 'movement vector ';
        
    case 'target location by origin'
        sub_title               = 'target position ';
        val_for_sub_assignment  = tar_val(u_tar_idx_idx,:);
        val_for_pos_assignment  = tar_val(u_tar_idx_idx,:);
        position_indexes        = tar_idx;
end

%% subplot positions
[subplot_pos, pop.columns, pop.rows]= DAG_dynamic_positions({val_for_sub_assignment});

%% some preallocations
pop.figure_title_part      =fig_title;
for n=1:numel(val_for_figure)
    temp     =[arrayfun(@num2str, val_for_figure{n}, 'unif', 0)'; repmat({' '},size(val_for_figure{n},1),1)'];
    pop.figure_title_value{n,1}     =[temp{:}];
end

%% Assigning each trial
pop.trial=o;
for t=1:numel(o)
    pop.trial(t).figure         =con_for_figure(t);
    pop.trial(t).column         =con_for_column(t);
    pop.trial(t).row            =con_for_row(t);
    pop.trial(t).fixation       =fixation_per_trial(t,:);
    pop.trial(t).title_part     =sub_title;
    pop.trial(t).subplot_pos    =subplot_pos(position_indexes(t));
    pop.trial(t).pos_index      =position_indexes(t);
    pop.trial(t).fix_index      =fixation_indexes(t);
    pop.trial(t).position       =val_for_pos_assignment(position_indexes(t),:);
    pop.trial(t).hemifield      =-1*(pop.trial(t).position(1)<0)+1*(pop.trial(t).position(1)>0);
end
end

function [s_c, displacement_types] = center_displacement_working(trial,keys)
Prec_fix=keys.cal.precision_fix;
Prec_pos=keys.cal.precision_tar;

movement            =NaN(size(trial'));
fixation            =NaN(size(trial'));
target              =NaN(size(trial'));
cuepos              =NaN(size(trial'));
stmpos              =NaN(size(trial'));

s_a=unique_positions([trial.fix_pos],Prec_fix);
s_b=unique_positions([trial.tar_pos] - [trial.fix_pos],Prec_pos);
s_c=unique_positions([trial.tar_pos],Prec_fix);
s_d=unique_positions([trial.cue_pos],Prec_pos);
if isfield(trial,'stm_pos')
    s_e=unique_positions([trial.stm_pos],Prec_pos);
else
    s_e=s_d;
end
% for t=1:numel(trial)
%     for k=1:numel(s_a)
%         if abs(trial(t).fix_pos - s_a(k)) < Prec_fix
%             fixation(t)=s_a(k);
%         end
%     end
%     for k=1:numel(s_b)
%         if abs(trial(t).tar_pos - trial(t).fix_pos - s_b(k)) < Prec_pos
%             movement_direction(t)=s_b(k);
%         end
%     end
%     for k=1:numel(s_c)
%         if abs(trial(t).tar_pos - s_c(k)) < Prec_fix
%             target(t)=s_c(k);
%         end
%     end
%     for k=1:numel(s_d)
%         if abs(trial(t).cue_pos - s_d(k)) < Prec_pos
%             cuepos(t)=s_d(k);
%         end
%     end
%     for k=1:numel(s_e)
%         if abs(trial(t).stm_pos - s_e(k)) < Prec_pos
%             stmpos(t)=s_e(k);
%         end
%     end
% end

for k=1:numel(s_a)
    t=abs([trial.fix_pos] - s_a(k)) < Prec_fix;
    fixation(t)=s_a(k);
end
for k=1:numel(s_b)
    t=abs([trial.tar_pos] - [trial.fix_pos] - s_b(k)) < Prec_pos;
    movement(t)=s_b(k);
end
for k=1:numel(s_c)
    t= abs([trial.tar_pos] - s_c(k)) < Prec_fix;
    target(t)=s_c(k);
end
for k=1:numel(s_d)
    t=abs([trial.cue_pos] - s_d(k)) < Prec_pos;
    cuepos(t)=s_d(k);
end
for k=1:numel(s_e)
    t=abs([trial(t).stm_pos] - s_e(k)) < Prec_pos;
    stmpos(t)=s_e(k);
end


[~,~,unique_condition]      =unique([real(fixation),imag(fixation),real(target),imag(target)],'rows');
[~,~,fixation_location]     =unique([real(fixation),imag(fixation)],'rows');
[~,~,movement_direction]    =unique([real(movement),imag(movement)],'rows');
[~,~,target_location]       =unique([real(target),imag(target)],'rows');
[~,~,cue_location]          =unique([real(cuepos),imag(cuepos)],'rows');
[~,~,stimulus_location]     =unique([real(stmpos),imag(stmpos)],'rows');

fix_y=imag(nanmean(fixation));
fixation=fixation-1i*fix_y;
target=target-1i*fix_y;

displacement_types=[real(fixation) imag(fixation) real(movement) imag(movement) real(target) imag(target) real(cuepos) imag(cuepos), real(stmpos) imag(stmpos),...
    unique_condition fixation_location movement_direction target_location, cue_location, stimulus_location]; 

if numel(trial)==0
    displacement_types=NaN(1,16);
end
end

function unique_target_positions=unique_positions(all_target_positions,Precision)
target_positions    =   unique(all_target_positions(~isnan(all_target_positions)));
n_targets           =   numel(target_positions);
for t=1:n_targets
    target_positions(abs(target_positions-target_positions(t))<Precision)=target_positions(t);
end
unique_target_positions=unique(target_positions);
end
