function [o,info]=ph_arrange_positions_and_plots(keys,o)
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
val_for_figure          = mov_val;
val_for_sub_assignment  = mov_val;
val_for_pos_assignment  = mov_val;
position_indexes        = mov_idx;
fixation_indexes        = fix_idx;
hemifield_indexes       = (real([o.tar_pos]' - [o.fix_pos]')>0)+1;
fixation_per_trial      = fix_val;


%% specifics
switch keys.arrangement
    
    case 'success_in_cue'
        fig_title               = '';
        val_for_figure          = cueshape;
        sub_title               = 'cue position';
        val_for_sub_assignment  = cue_val;
        val_for_pos_assignment  = cue_val;
        position_indexes        = cue_idx;
        
    case 'cue_position'
        con_for_figure          = suc_idx;
        fig_title               = 'success ';
        sub_title               = 'cue position';
        val_for_figure          = success;
        val_for_sub_assignment  = cue_val;
        val_for_pos_assignment  = cue_val;
        position_indexes        = cue_idx;
        
    case 'hand_choices'
        fig_title               = '';
        val_for_figure          = [hands choices];
        
    case 'options'
        con_for_figure          = hnd_idx;
        fig_title               = 'hand ';
        val_for_figure          = hands;
        
    case 'hands_inactivation'
        con_for_figure          = cho_idx;
        fig_title               = 'choice ';
        val_for_figure          = choices;
        [~,~,con_for_column]    = unique([hands hemifield_indexes],'rows');
        
    case 'hands_inactivation_in_ch'
        con_for_figure          = eff_idx;
        fig_title               = 'choice_instructed_comp ';
        val_for_figure          = effectors;
        [~,~,con_for_column]    = unique([hands hemifield_indexes],'rows');
        
    case 'hands_in_ch'
        con_for_figure          = eff_idx;
        fig_title               = 'choice_instructed_comp ';
        val_for_figure          = effectors;
        [~,~,con_for_column]    = unique([hands hemifield_indexes],'rows');
        
    case 'hands'
        con_for_figure          = cho_idx;
        con_for_row             = ones(size(con_for_figure));
        fig_title               = 'choice ';
        sub_title               = 'movement vector ';
        val_for_figure          = choices;
        con_for_column          = eff_idx;
        %Precision=5;
        
    case 'fixation'
        con_for_row             = eff_idx;
        sub_title               = 'fixation at ';
        val_for_sub_assignment  = fix_val;
        val_for_pos_assignment  = fix_val;
        position_indexes        = fix_idx;
        
    case 'movement vectors'
        sub_title               = 'movement vector ';
        
    case 'target location by origin'
        sub_title               = 'target position ';
        val_for_sub_assignment  = tar_val;
        val_for_pos_assignment  = tar_val;
        position_indexes        = tar_idx;
end

%% subplot positions fro singel cell plotting -> individually for each session (and tasktype ?)
[sessions,~,session_index]=unique([o.session]);
for s=sessions
    idx_s=session_index==s;
    typs=unique([o(idx_s).type]);
    for t=typs
        idx=idx_s & [o.type]'==t;
        u_val=unique(val_for_sub_assignment(idx,:),'rows');
        [~,~,p_idx]=unique(position_indexes(idx));
        [subplot_pos, columns, rows]= DAG_dynamic_positions({u_val});
        [o(idx).columns]     =deal(columns);
        [o(idx).rows]        =deal(rows);
        subplot_pos_cell     =num2cell(subplot_pos(p_idx));
        [o(idx).subplot_pos]   =deal(subplot_pos_cell{:});
    end
end
[o.figure_title_part]        =deal(fig_title);


%% figure title value seems to be quite annoying...
for n=1:numel(o)
    o(n).figure_title_value=num2str(val_for_figure(n,:));
end

%% Assigning each trial
position            =num2cell(val_for_pos_assignment,2);
fixation            =num2cell(fixation_per_trial,2);
hemifield           =num2cell(sign(val_for_pos_assignment(:,1)));

con_for_figure      =num2cell(con_for_figure);
con_for_column      =num2cell(con_for_column);
con_for_row         =num2cell(con_for_row);
position_indexes    =num2cell(position_indexes);
fixation_indexes    =num2cell(fixation_indexes);

[o.title_part]     =deal(sub_title);
[o.figure]         =deal(con_for_figure{:});
[o.column]         =deal(con_for_column{:});
[o.row]            =deal(con_for_row{:});
[o.pos_index]      =deal(position_indexes{:});
[o.fix_index]      =deal(fixation_indexes{:});

[o.fixation]       =deal(fixation{:});
[o.position]       =deal(position{:});
[o.hemifield]      =deal(hemifield{:});
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
