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
keys.precision.fixation=4;
keys.precision.position=2;
[~, displacement_types] = center_displacement_working(o,keys);
all_val=displacement_types(:,1:4);
fix_val=displacement_types(:,1:2);
tar_val=displacement_types(:,3:4);
cue_val=displacement_types(:,5:6);
mov_val=displacement_types(:,3:4) - displacement_types(:,1:2);
stm_val=displacement_types(:,7:8);

all_idx=displacement_types(:,9);
fix_idx=displacement_types(:,10);
mov_idx=displacement_types(:,11);
tar_idx=displacement_types(:,12);
cue_idx=displacement_types(:,13);
stm_idx=displacement_types(:,14);
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
%
%% KK STUFF
%% spatial distractor task
difficulty = [o.difficulty];
stimuli_in_2hemifields = [o.stimuli_in_2hemifields];
n_distractors = [o.n_distractors];
n_nondistractors = [o.n_nondistractors];
StimulusType = [o.stimulustype];


%% KK distractor task - define new index
%1. Nb of stimuli -> 2 or 3
Settings.TaskParameter.target_color_dim                    = [128,0,0];
Settings.TaskParameter.fixation_color_dim                  = [60,60,60];
Settings.TaskParameter.analyze_distr_colors                = 'any';

num_targets_per_trial   = cellfun(@(x) sum(ismember(x,Settings.TaskParameter.target_color_dim,'rows')),{o.col_dim},'UniformOutput',1);
num_distr_per_trial     = cellfun(@(x) sum(~ismember(x,Settings.TaskParameter.fixation_color_dim,'rows') & ~ismember(x,Settings.TaskParameter.target_color_dim,'rows')),{o.col_dim},'UniformOutput',1);
StimulusType            = sum([num_targets_per_trial; num_distr_per_trial] );

for i = 1 : length(StimulusType)
    if num_targets_per_trial(i) == 0 && num_distr_per_trial(i) == 1 ||  num_targets_per_trial(i) == 1 && num_distr_per_trial(i) == 0
        StimulusType(i) = 1;
    elseif num_targets_per_trial(i) == 2 || num_distr_per_trial(i) == 2
        StimulusType(i) = 2;
    elseif num_targets_per_trial(i) == 1 && num_distr_per_trial(i) == 1
        StimulusType(i) = 3;
    end
end


trial_cond.single_targ = find(num_targets_per_trial == 1 & num_distr_per_trial == 0); % target only trials
trial_cond.single_distr = find(num_targets_per_trial == 0 & num_distr_per_trial == 1 ); % distr only trials

%trial_cond.targ_targ_2HF = find(num_targets_per_trial == 2 & x_pos1HF_stimuli' == 3); % target-target trials

% distractor_color
distr_color_per_trial = cellfun(@(x) x(~ismember(x,Settings.TaskParameter.fixation_color_dim,'rows') & ~ismember(x,Settings.TaskParameter.target_color_dim,'rows'),:),{o.col_dim},'UniformOutput',0);


if isequal(Settings.TaskParameter.analyze_distr_colors,'any')
    num_distr_colors2ana = 5; %Color of distractors is hard coded
else
    num_distr_colors2ana = size(Settings.analyze_distr_colors,1);
end
distr_colors_all = [NaN,NaN,NaN];
idx = 1;
distr_ID = 0;
distractors = nan(num_distr_colors2ana, size([o.n],2));
for trial_no = 1:size([o.n],2)
    curr_distr_colors = distr_color_per_trial{trial_no};
    if size(distr_color_per_trial{trial_no},1) > 1
        if ~isequal(curr_distr_colors(1,:),curr_distr_colors(2,:));
            warning(char(strcat(sprintf('There were two distractors with different colors in trial %d in',trial_no),{' '},file_to_load,{'. The one with index 1 was used for analysis!'})));
        end
        curr_distr_colors = curr_distr_colors(1,:);
    end
    % This will overwrite logicals in "distractors" with the
    % last distractor color per trial in "distr_color_per_trial"
    
    if ~ismember(curr_distr_colors,distr_colors_all,'rows') % "distr_colors_all" is the sequence of appearance of the distractor colors in one run
        distr_colors_all(idx,:) = curr_distr_colors;
        idx = idx + 1;
        distr_ID = distr_ID + 1;
        distractors(distr_ID,trial_no) = trial_no; % "distractors": rows represent distractor colors in order of appearance as in "distr_colors_all"
    else
        [~,distr_ID_member] = ismember(curr_distr_colors,distr_colors_all,'rows');
        distractors(distr_ID_member,trial_no) = trial_no;
    end
end;
distr_colors2ana = distr_colors_all;
num_distr_colors2ana = size(distr_colors2ana,1);

% calculate G/R ratio of RGB values
G_R_ratio = nan(num_distr_colors2ana,1);
for col = 1:num_distr_colors2ana
    G_R_ratio(col,1) = distr_colors2ana(col,2) / distr_colors2ana(col,1);
end

% sort data according to G/R ratio: from yellow (small ratio) to red (high ratio)
[G_R_ratio_sorted,ind]     = sort(G_R_ratio);
distr_colors2ana_sorted    = distr_colors2ana(ind,:);
distractor_color           = distractors(ind,:); % indices for each distractor color(rows) in each trial (columns)


%% Each row contains one Difficulty level (3 = only target, 1 = Diff to 2 = Easy)
AllDifficulties= distractor_color;difficulty= [];

for i_diff = 1: length(ind)
    AllDifficulties(i_diff, ~isnan(distractor_color( i_diff ,:))) = i_diff ;
    AllDifficulties(i_diff, isnan(distractor_color( i_diff ,:))) = 0;  % distr only trials
end
AllDifficulties(i_diff +1, sum(AllDifficulties,1 ) == 0) = i_diff +1; %all only target trials

% if length(G_R_ratio_sorted) > 2
% AllDifficulties(4,(AllDifficulties(4,:) == 4)) = 3;
% AllDifficulties(3,:) = 0;
% AllDifficulties(3,:) = AllDifficulties(4,:);
% AllDifficulties(4,:)= [];
% end
difficulty = sum(AllDifficulties, 1);

target_selected      =[o.target_selected]'; %Which target selected

% 
% disp(['Nb of sgl stimuli ', num2str(numel(difficulty(StimulusType == 1)))])
% disp(['Nb of sgl stimuli - target ', num2str(sum(sum([difficulty == 3 ;  StimulusType == 1])== 2))])
% disp(['Nb of sgl stimuli easy distr', num2str(sum(sum([difficulty == 2 ;  StimulusType == 1])== 2))])
% disp(['Nb of sgl stimuli diff distr', num2str(sum(sum([difficulty == 1 ;  StimulusType == 1])== 2))])
% 
% disp(['Nb of dbl same stimuli ', num2str(numel(difficulty(StimulusType == 2)))])
% disp(['Nb of dbl same stimuli - target ', num2str(sum(sum([ difficulty == 3 ;  StimulusType == 2])== 2))])
% disp(['Nb of dbl same stimuli - target ', num2str(sum(sum([target_selected' == 1 ; difficulty == 3 ;  StimulusType == 2])== 3))])
% disp(['Nb of dbl same stimuli - target ', num2str(sum(sum([target_selected' == 2 ;difficulty == 3 ;  StimulusType == 2])== 3))])
% disp(['Nb of dbl same stimuli - target ', num2str(sum(sum([target_selected' == 3 ;difficulty == 3 ;  StimulusType == 2])== 3))])
% 
% disp(['Nb of dbl same stimuli - easy distr ', num2str(sum(sum([difficulty == 2 ;  StimulusType == 2])== 2))])
% disp(['Nb of dbl same stimuli - diff distr ', num2str(sum(sum([difficulty == 1 ;  StimulusType == 2])== 2))])
% 
% disp(['Nb of dbl same stimuli ', num2str(numel(difficulty(StimulusType == 2)))])
% 
% 
% disp(['Nb of dbl diff stimuli ', num2str(numel(difficulty(StimulusType == 3)))])


% red_targets_per_trial = cellfun(@(x) find(ismember(x,Settings.TaskParameter.target_color_dim,'rows')),{o.col_dim},'UniformOutput',0);
% distr_per_trial = cellfun(@(x) find(~ismember(x,Settings.TaskParameter.fixation_color_dim,'rows') & ~ismember(x,Settings.TaskParameter.target_color_dim,'rows')),{o.col_dim},'UniformOutput',0);
% fix_per_trial = cellfun(@(x) find(ismember(x,Settings.TaskParameter.fixation_color_dim,'rows')),{o.col_dim},'UniformOutput',0);
%
%
% red_targets_per_trial_empty = cellfun(@(x) isempty(x),red_targets_per_trial,'UniformOutput',1);
% two_red_targets_per_trial = cellfun(@(x) isequal(size(x,1),2),red_targets_per_trial,'UniformOutput',1);
% two_distr_per_trial = cellfun(@(x) isequal(size(x,1),2),distr_per_trial,'UniformOutput',1);
%
% correct_target_per_trial = red_targets_per_trial;
% correct_target_per_trial(red_targets_per_trial_empty) = fix_per_trial(red_targets_per_trial_empty);
%
% red_targ_or_distr_per_trial = red_targets_per_trial; %sgl stimuli
% red_targ_or_distr_per_trial(red_targets_per_trial_empty) = distr_per_trial(red_targets_per_trial_empty);
% non_fix_stimuli_per_trial = red_targ_or_distr_per_trial; % also includes targ_targ and distr_distr trials
% red_targ_or_distr_per_trial(two_repd_targets_per_trial) = {[]};
% red_targ_or_distr_per_trial(two_distr_per_trial) = {[]}; % contains indices of the target (if present and only if it is a single target or a target-distractor trial)
%






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
[hnd_cho_values,~,hch_idx]  =unique([hands choices],'rows');
[hnd_ptb_values,~,hpt_idx]  =unique([hands perturbations],'rows');
[hnd_blo_values,~,hbl_idx]  =unique([hands blocks],'rows');
[hnd_eff_values,~,hef_idx]  =unique([hands effectors],'rows');
[cho_ptb_values,~,chp_idx]  =unique([choices perturbations],'rows');

%[diff_values,~,diff_idx]      =unique(difficulty+3*success');
[Styp_values,~,Styp_idx]      =unique(StimulusType');




%% labels
hand_labels         ={'NH','NH CH','LH','LH CH','RH','RH CH'};
choice_labels       ={'IN','CH','','','',''};
hand_choice_labels  ={'NH IN','NH CH','LH IN','LH CH','RH IN','RH CH'};
hand_ptb_labels     =cellstr([repmat('NH ',numel(blo_values_o),1) num2str(blo_values_o); repmat('LH ',numel(blo_values_o),1) num2str(blo_values_o); repmat('RH ',numel(blo_values_o),1) num2str(blo_values_o)])';
hand_eff_labels     =cellstr([repmat('NH ',numel(eff_values),1) num2str(eff_values); repmat('LH ',numel(eff_values),1) num2str(eff_values); repmat('RH ',numel(eff_values),1) num2str(eff_values)])';
fix_labels          =num2cell(num2str(unique(fix_idx)))';
cue_shape_labels ={};
for n_shp=1:numel(shp_values)
    cue_shape_labels{n_shp}=sprintf('%.1f',shp_values(n_shp));
end
%% KK stuff
[diff_values,~,diff_idx]        = unique(difficulty);
[diff_values,~,diff_idx]      =unique(difficulty+3*success');
[Styp_values,~,Styp_idx]      =unique(StimulusType');

%
%
% %% labels
% hand_labels         ={'NH','NH CH','LH','LH CH','RH','RH CH'};
% choice_labels       ={'IN','CH','','','',''};
% hand_choice_labels  ={'NH IN','NH CH','LH IN','LH CH','RH IN','RH CH'};
% hand_ptb_labels     =cellstr([repmat('NH ',numel(blo_values_o),1) num2str(blo_values_o); repmat('LH ',numel(blo_values_o),1) num2str(blo_values_o); repmat('RH ',numel(blo_values_o),1) num2str(blo_values_o)])';
% hand_eff_labels     =cellstr([repmat('NH ',numel(eff_values),1) num2str(eff_values); repmat('LH ',numel(eff_values),1) num2str(eff_values); repmat('RH ',numel(eff_values),1) num2str(eff_values)])';
% fix_labels          =num2cell(num2str(unique(fix_idx)))';
% cue_shape_labels ={};
% for n_shp=1:numel(shp_values)
%     cue_shape_labels{n_shp}=sprintf('%.1f',shp_values(n_shp));
% end

% %% color assignment
% cols=keys.colors;
% keys.hnd_choice_colors_L=[cols.NH_LS_IN;cols.NH_LS_CH;cols.LH_LS_IN;cols.LH_LS_CH;cols.RH_LS_IN;cols.RH_LS_CH]/255;
% keys.hnd_choice_colors_R=[cols.NH_RS_IN;cols.NH_RS_CH;cols.LH_RS_IN;cols.LH_RS_CH;cols.RH_RS_IN;cols.RH_RS_CH]/255;
% keys.hnd_choice_colors  =[cols.NH_IN;cols.NH_CH;cols.LH_IN;cols.LH_CH;cols.RH_IN;cols.RH_CH]/255;
% b_pt=(max(blo_values):-1*max(blo_values)/2/(numel(blo_values)-1):max(blo_values)/2)'/max(blo_values);
% b_pt(isnan(b_pt))=1; %% for the case of onlý one block... -.-
% keys.hnd_ptb_colors_L=[b_pt*cols.NH_LS_IN;b_pt*cols.LH_LS_IN;b_pt*cols.RH_LS_IN]/255;
% keys.hnd_ptb_colors_R=[b_pt*cols.NH_RS_IN;b_pt*cols.LH_RS_IN;b_pt*cols.RH_RS_IN]/255;
% keys.hnd_ptb_colors  =[b_pt*cols.NH_IN   ;b_pt*cols.LH_IN   ;b_pt*cols.RH_IN   ]/255;
%
% % hnd_set_values(hnd_set_values(:,2)==min(hnd_set_values(:,2)),2)=1;
% % hnd_set_values(hnd_set_values(:,2)==max(hnd_set_values(:,2)),2)=2;
% hand_choice_color_combination=[0 0; 0 1; 1 0; 1 1; 2 0; 2 1];
% hand_blo_color_combination   =fliplr(combvec(unique(blocks)',[0 1 2])');
% hand_eff_color_combination   =fliplr(combvec(unique(effectors)',[0 1 2])');
% %hand_pt_color_combination   =[0 1; 0 2; 0 3;1 1; 1 2; 1 3; 2 1; 2 2; 2 3];
% hand_color_combination=[0; NaN; 1; NaN; 2; NaN];
% choice_color_combination=[0; 1; NaN; NaN; NaN; NaN];
%
% b_ef=(numel(eff_values):-1*numel(eff_values)/2/(numel(eff_values)-1):numel(eff_values)/2)'/numel(eff_values);
% b_ef(isnan(b_ef))=1; %% for the case of onlý one block... -.-
% keys.hnd_eff_colors_L=[b_ef*cols.NH_LS_IN;b_ef*cols.LH_LS_IN;b_ef*cols.RH_LS_IN]/255;
% keys.hnd_eff_colors_R=[b_ef*cols.NH_RS_IN;b_ef*cols.LH_RS_IN;b_ef*cols.RH_RS_IN]/255;
% keys.hnd_eff_colors  =[b_ef*cols.NH_IN   ;b_ef*cols.LH_IN   ;b_ef*cols.RH_IN   ]/255;

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
fixation_indexes            = fix_idx;
hemifield_indexes       = (real([o.tar_pos]' - [o.fix_pos]')>0)+1;
fixation_per_trial      = fix_val;


%% Position for the spatial distractor task

%% specifics
switch keys.arrangement
    
    
    case 'StimulusType_Difficulty_Position'
        difficulty =  difficulty+3*success';
        %              difficulty =  difficulty(success ==1)';
        %              idx_difficulty = find(difficulty <3);
        %              difficulty(idx_difficulty)     =[];
        %              stm_idx(idx_difficulty)        =[];
        %              Styp_idx(idx_difficulty)       =[];
        
        % if length(G_R_ratio_sorted) > 2
        %     difficulty(difficulty == 0) = [];
        % end
        [diff_values,~,diff_idx]        = unique(difficulty);
        
        
        
        fig_title               = 'StimulusType_Difficulty_Position';
        con_for_figure          = Styp_idx;
        val_for_figure          = num2cell(Styp_values);
        
        con_for_line            = diff_idx';  %%% there was no ' !?
        if length(diff_values) == 3
            pop.line_labels        =   {'Diff', 'Easy','Tar'};
            
        else
            pop.line_labels        =   {'Diff','Diff2', 'Easy','Tar'};
        end
        
        position_indexes        = stm_idx;
        val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
        val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
        sub_title               = 'stimulus position';
        if length(diff_values) == 3
            col_left = autumn(3);
            col_right = winter(3);
            pop.PSTH_perpos_colors =   [[col_left(2,:);col_left(3,:);col_left(1,:)];[col_right(2,:);col_right(3,:);col_right(1,:)]] ;
            pop.PSTH_summary_colors=   [[col_left(2,:);col_left(3,:);col_left(1,:)];[col_right(2,:);col_right(3,:);col_right(1,:)]] ;
        else
            pop.PSTH_perpos_colors =   [autumn(4);winter(4)] ;
            pop.PSTH_summary_colors=   [autumn(4);winter(4)] ;
        end
        hemifield_indexes       = (real(stm_val(stm_idx))>0)+1;
        
        %
        % %     disp(['Nb of sgl stimuli - Diff R ', num2str(sum(sum([difficulty == 1 ;  StimulusType == 1 ;  hemifield_indexes' == 1])== 3))])
        % %     disp(['Nb of sgl stimuli - Diff L ', num2str(sum(sum([difficulty == 1 ;  StimulusType == 1 ;  hemifield_indexes' == 2])== 3))])
        % %
        % %     disp(['Nb of sgl stimuli - Diff R ', num2str(sum(sum([difficulty == 2 ;  StimulusType == 1 ;  hemifield_indexes' == 1])== 3))])
        % %     disp(['Nb of sgl stimuli - Diff L ', num2str(sum(sum([difficulty == 2 ;  StimulusType == 1 ;  hemifield_indexes' == 2])== 3))])
        % %
        % %     disp(['Nb of sgl stimuli - target R ', num2str(sum(sum([difficulty == 3 ;  StimulusType == 1 ;  hemifield_indexes' == 1])== 3))])
        % %     disp(['Nb of sgl stimuli - target L ', num2str(sum(sum([difficulty == 3 ;  StimulusType == 1 ;  hemifield_indexes' == 2])== 3))])
        %
        %      con_for_trial_crit      = con_for_line;
        %
        %
        %                 case 'StimulusType_Difficulty_Position_Successful'
        %               difficulty =  difficulty(success ==1);
        %               stm_idx  =  stm_idx(success ==1);
        %               Styp_idx =  Styp_idx(success ==1);
        %               stm_val  =  stm_val(success ==1,:);
        %               [~,u_stm_idx_idx]=unique(stm_idx);
        %
        %             [diff_values,~,diff_idx]        =unique(difficulty);
        %
        %      fig_title               = 'StimulusType_Difficulty_Position_Successful';
        %      con_for_figure          = Styp_idx;
        %      val_for_figure          = num2cell(Styp_values);
        %
        %      con_for_line            = diff_idx';  %%% there was no ' !?
        %      pop.line_labels        =   {'Diff','Easy','Tar'};
        %      position_indexes        = stm_idx;
        %      val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
        %      val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
        %      sub_title               = 'stimulus position';
        %      pop.PSTH_perpos_colors =   [autumn(3);winter(3)] ;
        %      pop.PSTH_summary_colors=   [autumn(3);winter(3)] ;
        %     hemifield_indexes       = (real(stm_val(stm_idx))>0)+1;
        %
        %    % (real([o.cue_pos]' - [o.tar_pos]')>0)+1;
        %      StimulusType  =  StimulusType(success ==1);
        %
        %     disp(['Nb of sgl stimuli - Diff R ', num2str(sum(sum([difficulty == 1 ;  StimulusType == 1 ;  hemifield_indexes' == 1])== 3))])
        %     disp(['Nb of sgl stimuli - Diff L ', num2str(sum(sum([difficulty == 1 ;  StimulusType == 1 ;  hemifield_indexes' == 2])== 3))])
        %
        %     disp(['Nb of sgl stimuli - Diff R ', num2str(sum(sum([difficulty == 2 ;  StimulusType == 1 ;  hemifield_indexes' == 1])== 3))])
        %     disp(['Nb of sgl stimuli - Diff L ', num2str(sum(sum([difficulty == 2 ;  StimulusType == 1 ;  hemifield_indexes' == 2])== 3))])
        %
        %     disp(['Nb of sgl stimuli - target R ', num2str(sum(sum([difficulty == 3 ;  StimulusType == 1 ;  hemifield_indexes' == 1])== 3))])
        %     disp(['Nb of sgl stimuli - target L ', num2str(sum(sum([difficulty == 3 ;  StimulusType == 1 ;  hemifield_indexes' == 2])== 3))])
        %
        %      con_for_trial_crit      = con_for_line;
        %
        %                case 'DoubleSameTargets_Position'
        %    [diff_values,~,diff_idx]      =unique(difficulty+3*success');
        %
        %      fig_title               = 'StimulusType';
        %      con_for_figure          = Styp_idx;
        %      val_for_figure          = num2cell(Styp_values);
        %
        %      con_for_line            = diff_idx';
        %      pop.line_labels        =   {'EDiff','EEasy','ETar','CDiff','CEasy','CTar'};
        %
        %      position_indexes        = stm_idx;
        %      val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
        %      val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
        %      sub_title               = 'stimulus position';
        %      pop.PSTH_perpos_colors =   [summer(6);winter(6)] ;
        %      pop.PSTH_summary_colors=   [summer(6);winter(6)] ;
        %     hemifield_indexes       = (real(stm_val(stm_idx))>0)+1;
        %      con_for_trial_crit      = con_for_line;
        %     % con_for_figure          = suc_idx;
        %
        %         case 'StimulusType_Difficulty_Position_ErrorVsCorrect'
        %                 difficulty =  difficulty+3*success';
        % %              difficulty =  difficulty(success ==1)';
        % %              idx_difficulty = find(difficulty <3);
        % %              difficulty(idx_difficulty)     =[];
        % %              stm_idx(idx_difficulty)        =[];
        % %              Styp_idx(idx_difficulty)       =[];
        %
        %      [diff_values,~,diff_idx]        = unique(difficulty);
        %      fig_title               = 'StimulusType_Difficulty_Position_ErrorVsCorrect';
        %      con_for_figure          = Styp_idx;
        %      val_for_figure          = num2cell(Styp_values);
        %
        %      con_for_line            = diff_idx';  %%% there was no ' !?
        %      pop.line_labels        =   {'EDiff','EEasy','ETar','CDiff','CEasy','CTar'};
        %
        %      position_indexes        = stm_idx;
        %      val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
        %      val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
        %      sub_title               = 'stimulus position';
        %
        %      pop.PSTH_perpos_colors =   [summer(6);winter(6)] ;
        %      pop.PSTH_summary_colors=   [summer(6);winter(6)] ;
        %     hemifield_indexes       = (real(stm_val(stm_idx))>0)+1;
        %      con_for_trial_crit      = con_for_line;
        %
        %         case 'SglStimuli_Difficulty_Position_ErrorVsCorrect'
        %      val_for_figure          = {Styp_values};
        %      con_for_line            = diff_idx';
        %      position_indexes        = stm_idx;
        %      val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
        %      val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
        %      sub_title               = 'stimulus position';
        %      pop.PSTH_perpos_colors =   [autumn(5);winter(5)] ;
        %      pop.PSTH_summary_colors=   [autumn(5);winter(5)] ;
        %     hemifield_indexes       = (real(stm_val(stm_idx))>0)+1;
        %      con_for_trial_crit      = con_for_line;
        %
        %            case 'StimulusType_Difficulty_Position_ErrorVsCorrect_lastVers'
        %    [diff_values,~,diff_idx]      =unique(difficulty+3*success');
        %
        %      fig_title               = 'StimulusType_Difficulty_Position';
        %      con_for_figure          = Styp_idx';
        %      val_for_figure          = num2cell(Styp_values);
        %
        %      con_for_line            = diff_idx';
        %      pop.line_labels        =   {'EDiff','EEasy','ETar','CDiff','CEasy','CTar'};
        %
        %      position_indexes        = stm_idx;
        %      val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
        %      val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
        %      sub_title               = 'stimulus position';
        %      pop.PSTH_perpos_colors =   [summer(6);winter(6)] ;
        %      pop.PSTH_summary_colors=   [summer(6);winter(6)] ;
        %     hemifield_indexes       = (real(stm_val(stm_idx))>0)+1;
        %      con_for_trial_crit      = con_for_line;
        %
        %
        %     case 'success_in_cue'
        %         con_for_line            = suc_idx;
        %
        % %% specifics
            
            case 'StimulusType_Difficulty_Position'
                fig_title               = 'StimType';
                con_for_figure          = Styp_idx;
                val_for_figure          = num2cell(Styp_values);
                
                con_for_line            = diff_idx';  %%% there was no ' !?
                %         if length(diff_values) == 3
                %             pop.line_labels        =   {'Tar', 'Easy','Diff'}; %{'Diff', 'Easy','Tar'};
                %
                %         else
                %             pop.line_labels        =   {'Diff','Diff2', 'Easy','Tar'};
                %         end
                
                %position_indexes        = stm_idx;
                %val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
                %val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
                
                %      position_indexes = mov_idx;
                %      val_for_sub_assignment  = tar_val(u_mov_idx_idx,:);
                %      val_for_pos_assignment  = mov_val(u_mov_idx_idx,:);
                
                val_for_sub_assignment  = tar_val(u_tar_idx_idx,:);
                val_for_pos_assignment  = tar_val(u_tar_idx_idx,:);
                position_indexes        = tar_idx;
                
                sub_title               = 'stimulus position';
                %         if length(diff_values) == 3
                %             col_left = autumn(3);
                %             col_right = winter(3);
                %             tar_purple = [1 0 1 ];
                % %
                % %             pop.PSTH_perpos_colors =   [[col_left(2,:);col_left(3,:);col_left(1,:)];[col_right(1,:);col_right(2,:);tar_purple]] ;
                % %             pop.PSTH_summary_colors=   [[col_left(2,:);col_left(3,:);col_left(1,:)];[col_right(1,:);col_right(2,:);tar_purple]] ;
                %             pop.PSTH_perpos_colors =   [autumn(4);winter(4)] ;
                %             pop.PSTH_summary_colors=   [autumn(4);winter(4)] ;
                %         else
                %             pop.PSTH_perpos_colors =   [autumn(4);winter(4)] ;
                %             pop.PSTH_summary_colors=   [autumn(4);winter(4)] ;
                %         end
                %             pop.PSTH_perpos_colors =   [autumn(6);winter(6)] ;
                %             pop.PSTH_summary_colors=   [autumn(6);winter(6)] ;
                hemifield_indexes       = (real(stm_val(stm_idx))>0)+1;
                con_for_trial_crit      = con_for_line;
                
            case 'StimulusType_Difficulty_Position_Successful'
                [diff_values,~,diff_idx]        = unique(difficulty);
                
                
                
                fig_title               = 'StimuType_Diff_Pos_Suc';
                con_for_figure          = Styp_idx;
                val_for_figure          = num2cell(Styp_values);
                
                con_for_line            = diff_idx';  %%% there was no ' !?
                if length(diff_values) == 3
                    pop.line_labels        =   {'Diff', 'Easy','Tar'};
                    
                else
                    pop.line_labels        =   {'Diff','Diff2', 'Easy','Tar'};
                end
                
                position_indexes        = stm_idx;
                val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
                val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
                sub_title               = 'stimulus position';
                
                col_left = autumn(3);
                col_right = winter(3);
                tar_purple = [1 0 1 ];
                
                pop.PSTH_perpos_colors =   [[col_left(2,:);col_left(3,:);col_left(1,:)];[col_right(1,:);col_right(2,:);tar_purple]] ;
                pop.PSTH_summary_colors=   [[col_left(2,:);col_left(3,:);col_left(1,:)];[col_right(1,:);col_right(2,:);tar_purple]] ;
                
                hemifield_indexes       = (real(stm_val(stm_idx))>0)+1;
                con_for_trial_crit      = con_for_line;
                
            case 'DoubleSameTargets_Position'
                [diff_values,~,diff_idx]      =unique(difficulty+3*success');
                
                fig_title               = 'StimulusType';
                con_for_figure          = Styp_idx;
                val_for_figure          = num2cell(Styp_values);
                
                con_for_line            = diff_idx';
                pop.line_labels        =   {'EDiff','EEasy','ETar','CDiff','CEasy','CTar'};
                
                position_indexes        = stm_idx;
                val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
                val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
                sub_title               = 'stimulus position';
                pop.PSTH_perpos_colors =   [summer(6);winter(6)] ;
                pop.PSTH_summary_colors=   [summer(6);winter(6)] ;
                hemifield_indexes       = (real(stm_val(stm_idx))>0)+1;
                con_for_trial_crit      = con_for_line;
                % con_for_figure          = suc_idx;
                
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
                
                con_for_line            = diff_idx';
                pop.line_labels        =   {'EDiff','EEasy','ETar','CDiff','CEasy','CTar'};
                
                position_indexes        = stm_idx;
                val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
                val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
                sub_title               = 'stimulus position';
                
                
                col_left = autumn(10);
                col_right = winter(6);
                tar_pink   = [0.4235    0.2510    0.3922];
                tar_purple = [1 0 1 ];
                Diff_R_easy = [0.2000    1.0000    0.8000 ];
                Diff_L_easy = [0.6000    0.8000    1.0000 ];
                
                
                pop.PSTH_perpos_colors =   [[col_left(4,:);col_left(7,:);col_left(2,:)  ;col_left(5,:); col_left(10,:)  ;col_left(1,:)];    [col_right(1,:);Diff_L_easy;tar_pink; col_right(4,:);Diff_R_easy;tar_purple]] ;
                pop.PSTH_summary_colors=   [[col_left(4,:);col_left(7,:);col_left(2,:)  ;col_left(5,:); col_left(10,:)  ;col_left(1,:)];    [col_right(1,:);Diff_L_easy;tar_pink; col_right(4,:);Diff_R_easy;tar_purple]] ;
                % uisetcolor([0.6 0.8 1])
                % pop.PSTH_perpos_colors =   [autumn(6);winter(6)] ;
                % pop.PSTH_summary_colors=   [autumn(6);winter(6)] ;
                hemifield_indexes       = (real(stm_val(stm_idx))>0)+1;
                con_for_trial_crit      = con_for_line;
                
                [diff_values,~,diff_idx]        = unique(difficulty);
                
                fig_title               = 'StimulusType_Difficulty_Position';
                con_for_figure          = Styp_idx;
                val_for_figure          = num2cell(Styp_values);
                con_for_line            = diff_idx';  %%% there was no ' !?
                if length(diff_values) == 3
                    pop.line_labels        =   {'Tar', 'Easy','Diff'}; %{'Diff', 'Easy','Tar'};
                else
                    pop.line_labels        =   {'Tar', 'Easy','Diff','Diff2',};
                end
                
                %position_indexes        = stm_idx;
                %val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
                %val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
                
                %      position_indexes = mov_idx;
                %      val_for_sub_assignment  = tar_val(u_mov_idx_idx,:);
                %      val_for_pos_assignment  = mov_val(u_mov_idx_idx,:);
                
                val_for_sub_assignment  = tar_val(u_tar_idx_idx,:);
                val_for_pos_assignment  = tar_val(u_tar_idx_idx,:);
                position_indexes        = tar_idx;
                
                sub_title               = 'stimulus position';
                if length(diff_values) == 3
                    col_left      = autumn(6);
                    col_right     = winter(3);
                    tar_purple    = [0.5    0.2510    0.3922];
                    col_fix       = gray(6);
                    
                    pop.PSTH_perpos_colors =   [[col_left(1,:);col_left(6,:);col_left(3,:)];[tar_purple; col_right(2,:);col_right(1,:)]] ;
                    pop.PSTH_summary_colors=   [[col_left(1,:);col_left(6,:);col_left(3,:)];[tar_purple ; col_right(2,:);col_right(1,:)]] ;
                    
                    pop.PSTH_perpos_colors =   [[col_left(1,:);col_left(6,:);col_left(3,:)];[col_fix(1,:);col_fix(4,:);col_fix(3,:)];[tar_purple; col_right(2,:);col_right(1,:)]] ;
                    pop.PSTH_summary_colors=   [[col_left(1,:);col_left(6,:);col_left(3,:)];[col_fix(1,:);col_fix(4,:);col_fix(3,:)];[tar_purple ; col_right(2,:);col_right(1,:)]] ;
                else
                    pop.PSTH_perpos_colors =   [autumn(4);winter(4)] ;
                    pop.PSTH_summary_colors=   [autumn(4);winter(4)] ;
                end
                hemifield_indexes       = (real(stm_val(stm_idx))>0)+1;
                con_for_trial_crit      = con_for_line;
                
            case 'StimulusType_Difficulty_Position_Successful'
                [diff_values,~,diff_idx]        = unique(difficulty);
                
                
                
                fig_title               = 'StimuType_Diff_Pos_Suc';
                con_for_figure          = Styp_idx;
                val_for_figure          = num2cell(Styp_values);
                
                con_for_line            = diff_idx';  %%% there was no ' !?
                if length(diff_values) == 3
                    pop.line_labels        =   {'Tar', 'Easy','Diff'}; %{'Diff', 'Easy','Tar'};
                else
                    pop.line_labels        =   {'Tar', 'Easy','Diff','Diff2',};
                end
                
                position_indexes        = stm_idx;
                val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
                val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
                sub_title               = 'stimulus position';
                
                if length(diff_values) == 3
                    col_left      = autumn(6);
                    col_right     = winter(3);
                    tar_purple    = [0.5    0.2510    0.3922];
                    col_fix       = gray(6);
                    
                    pop.PSTH_perpos_colors =   [[col_left(1,:);col_left(6,:);col_left(3,:)];[tar_purple; col_right(2,:);col_right(1,:)]] ;
                    pop.PSTH_summary_colors=   [[col_left(1,:);col_left(6,:);col_left(3,:)];[tar_purple ; col_right(2,:);col_right(1,:)]] ;
                    
                    %pop.PSTH_perpos_colors =   [[col_left(1,:);col_left(6,:);col_left(3,:)];[col_fix(1,:);col_fix(4,:);col_fix(3,:)];[tar_purple; col_right(2,:);col_right(1,:)]] ;
                    %pop.PSTH_summary_colors=   [[col_left(1,:);col_left(6,:);col_left(3,:)];[col_fix(1,:);col_fix(4,:);col_fix(3,:)];[tar_purple ; col_right(2,:);col_right(1,:)]] ;
                else
                    pop.PSTH_perpos_colors =   [autumn(4);winter(4)] ;
                    pop.PSTH_summary_colors=   [autumn(4);winter(4)] ;
                end
                
                hemifield_indexes       = (real(stm_val(stm_idx))>0)+1;
                con_for_trial_crit      = con_for_line;
                
            case 'DoubleSameTargets_Position'
                [diff_values,~,diff_idx]      =unique(difficulty+3*success');
                
                fig_title               = 'StimulusType';
                con_for_figure          = Styp_idx;
                val_for_figure          = num2cell(Styp_values);
                
                con_for_line            = diff_idx';
                pop.line_labels        =   {'EDiff','EEasy','ETar','CDiff','CEasy','CTar'};
                
                position_indexes        = stm_idx;
                val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
                val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
                sub_title               = 'stimulus position';
                pop.PSTH_perpos_colors =   [summer(6);winter(6)] ;
                pop.PSTH_summary_colors=   [summer(6);winter(6)] ;
                hemifield_indexes       = (real(stm_val(stm_idx))>0)+1;
                con_for_trial_crit      = con_for_line;
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
                
                con_for_line            = diff_idx';  %%% there was no ' !?
                if length(diff_values) == 3
                    pop.line_labels        =   {'Tar', 'Easy','Diff'}; %{'Diff', 'Easy','Tar'};
                else
                    pop.line_labels        =   {'Tar', 'Easy','Diff','Diff2',};
                end
                
                
                sub_title               = 'stimulus position';
                
                %         pop.line_labels        =   {'Tar', 'Easy'}; %{'Diff', 'Easy','Tar'};
                %          col_left      = autumn(6);
                %             col_right     = winter(3);
                %             tar_purple    = [0.5    0.2510    0.3922];
                %             col_fix       = gray(6);
                %
                %             pop.PSTH_perpos_colors =   [[col_left(1,:);col_left(6,:)];[tar_purple; col_right(2,:)]] ;
                %             pop.PSTH_summary_colors=   [[col_left(1,:);col_left(6,:)];[tar_purple ; col_right(2,:)]] ;
                %                  position_indexes        = stm_idx;
                %                  val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
                %                  val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
                %
                
                if length(diff_values) == 3
                    fig_title               = 'StimTyp_Diff_Pos_Suc';
                    
                    col_left      = autumn(6);
                    col_right     = winter(3);
                    tar_purple    = [0.5    0.2510    0.3922];
                    col_fix       = gray(6);
                    
                    pop.PSTH_perpos_colors =   [[col_left(1,:);col_left(6,:);col_left(3,:)];[tar_purple; col_right(2,:);col_right(1,:)]] ;
                    pop.PSTH_summary_colors=   [[col_left(1,:);col_left(6,:);col_left(3,:)];[tar_purple ; col_right(2,:);col_right(1,:)]] ;
                    position_indexes        = stm_idx;
                    val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
                    val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
                    
                    %             pop.PSTH_perpos_colors =   [[col_left(1,:);col_left(6,:);col_left(3,:)];[col_fix(1,:);col_fix(4,:);col_fix(3,:)];[tar_purple; col_right(2,:);col_right(1,:)]] ;
                    %             pop.PSTH_summary_colors=   [[col_left(1,:);col_left(6,:);col_left(3,:)];[col_fix(1,:);col_fix(4,:);col_fix(3,:)];[tar_purple ; col_right(2,:);col_right(1,:)]] ;
                    %             val_for_sub_assignment  = tar_val(u_tar_idx_idx,:);
                    %             val_for_pos_assignment  = tar_val(u_tar_idx_idx,:);
                    %             position_indexes        = tar_idx;
                    
                else
                    pop.PSTH_perpos_colors =   [autumn(4);winter(4)] ;
                    pop.PSTH_summary_colors=   [autumn(4);winter(4)] ;
                end
                
                hemifield_indexes       = (real(stm_val(stm_idx))>0)+1;
                con_for_trial_crit      = con_for_line;
                
            case 'Sgl_Diff_Pos_Suc_SaccadeEpoch'
                [diff_values,~,diff_idx]        = unique(difficulty);
                fig_title               = 'StimTyp_Diff_Pos_Suc';
                con_for_figure          = Styp_idx;
                val_for_figure          = num2cell(Styp_values);
                
                val_for_sub_assignment  = tar_val(u_tar_idx_idx,:);
                val_for_pos_assignment  = tar_val(u_tar_idx_idx,:);
                position_indexes        = tar_idx;
                
                con_for_line            = diff_idx';  %%%
                pop.line_labels        =   {'Tar', 'Easy','Diff'};
                sub_title               = 'stimulus position';
                pop.PSTH_perpos_colors =   [autumn(3);winter(3)] ;
                pop.PSTH_summary_colors=   [autumn(3);winter(3)] ;
                hemifield_indexes       = (real(tar_val(stm_idx))>0)+1;
                con_for_trial_crit      = con_for_line;
            case 'DisTar_Diff_Pos_Suc'
                [diff_values,~,diff_idx]        = unique(difficulty);
                
                con_for_line           = diff_idx';
                pop.line_labels        =   {'T-Deasy', 'T-Ddiff'};
                pop.line_labels        =   {'T-Deasy', 'T-Ddiff'};
                
                sub_title               = 'stimulus position';
                
                col_left      = autumn(6);
                col_right     = winter(3);
                pop.PSTH_perpos_colors =   [[col_left(6,:);col_left(3,:)];[col_right(2,:);col_right(1,:)]] ;
                pop.PSTH_summary_colors=   [[col_left(6,:);col_left(3,:)];[col_right(2,:);col_right(1,:)]] ;
                
                %hemifield_indexes       = (real(stm_val(stm_idx))>0)+1;
                hemifield_indexes       = (real(tar_val(stm_idx))>0)+1;
                
                con_for_trial_crit      = con_for_line;
                
                val_for_sub_assignment  = tar_val(u_tar_idx_idx,:);
                val_for_pos_assignment  = tar_val(u_tar_idx_idx,:);
                position_indexes        = tar_idx;
                
                
            case 'SpatialCompetition_Targets'
                % left vs Right is defined where the saccade was made to
                con_for_line           = Styp_idx;
                pop.line_labels        =   {'SglT','DblT-2Hf','DblT-1HF'};
                pop.line_labels        =   {'SglT','DblT-2Hf','DblT-1HF'};
                
                sub_title               = 'stimulus position';
                
                col_left      = autumn(6);
                col_right     = winter(3);
                tar_purple    = [0.5    0.2510    0.3922];
                col_fix       = gray(6);
                
                %
                pop.PSTH_perpos_colors =   [[col_left(1,:);col_left(3,:);[0.6350, 0.0780, 0.1840] 	];    [tar_purple; col_right(2,:);col_right(1,:)]] ;
                pop.PSTH_summary_colors=   [[col_left(1,:);col_left(3,:);[0.6350, 0.0780, 0.1840] 	];    [tar_purple ; col_right(2,:);col_right(1,:) ]] ;
                
                hemifield_indexes       = (real(stm_val(stm_idx))>0)+1;
                con_for_trial_crit      = con_for_line;
                
                val_for_sub_assignment  = tar_val(u_tar_idx_idx,:);
                val_for_pos_assignment  = tar_val(u_tar_idx_idx,:);
                position_indexes        = tar_idx;
                %unique(stm_val(find(StimulusType == 4),:), 'rows' )
                %unique(stm_val(find(StimulusType == 2),:), 'rows' )
                
            case 'SpatialCompetition_Distractor'
                [diff_values,~,diff_idx]        = unique(difficulty);
                
                fig_title               = 'DifficultyLevel';
                con_for_figure          = diff_idx';
                val_for_figure          = num2cell(diff_values);
                
                con_for_line           = Styp_idx;
                pop.line_labels        =   {'SglD','DblD-2Hf','DblD-1HF'};
                
                sub_title               = 'stimulus position';
                pop.PSTH_perpos_colors =   [spring(3);winter(3)] ;
                pop.PSTH_summary_colors=   [spring(3);winter(3)] ;
                
                hemifield_indexes       = (real(stm_val(stm_idx))>0)+1;
                con_for_trial_crit      = con_for_line;
                
                position_indexes        = stm_idx;
                val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
                val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
                
            case 'StimType_Diff_Pos_ErVsCor'
                difficulty =  difficulty+3*success';
                
                [diff_values,~,diff_idx]        = unique(difficulty);
                fig_title               = 'StimType_Diff_Pos_ErVsCor';
                con_for_figure          = Styp_idx;
                val_for_figure          = num2cell(Styp_values);
                
                con_for_line            = diff_idx';
                pop.line_labels        =   {'EDiff','EEasy','ETar','CDiff','CEasy','CTar'};
                pop.line_labels        =   {'ETar','EEasy','EDiff','CTar','CEasy','CDiff'};
                
                position_indexes        = stm_idx;
                val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
                val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
                sub_title               = 'stimulus position';
                
                
                col_left = autumn(10);
                col_right = winter(6);
                tar_pink   = [0.4235    0.2510    0.3922];
                tar_purple = [1 0 1 ];
                Diff_R_easy = [0.2000    1.0000    0.8000 ];
                Diff_L_easy = [0.6000    0.8000    1.0000 ];
                
                pop.PSTH_perpos_colors =   [[col_left(4,:);col_left(7,:);col_left(2,:)  ;col_left(5,:); col_left(10,:)  ;col_left(1,:)];    [col_right(1,:);Diff_L_easy;tar_pink; col_right(4,:);Diff_R_easy;tar_purple]] ;
                pop.PSTH_summary_colors=   [[col_left(4,:);col_left(7,:);col_left(2,:)  ;col_left(5,:); col_left(10,:)  ;col_left(1,:)];    [col_right(1,:);Diff_L_easy;tar_pink; col_right(4,:);Diff_R_easy;tar_purple]] ;
                pop.PSTH_perpos_colors =   [[col_left(2,:);col_left(7,:);col_left(4,:)  ;col_left(1,:); col_left(10,:)  ;col_left(5,:)];    [tar_pink; Diff_L_easy; col_right(1,:);tar_purple; Diff_R_easy; col_right(4,:)]] ;
                pop.PSTH_summary_colors=   [[col_left(2,:);col_left(7,:);col_left(4,:)  ;col_left(1,:); col_left(10,:)  ;col_left(5,:)];    [tar_pink; Diff_L_easy; col_right(1,:);tar_purple; Diff_R_easy; col_right(4,:)]] ;
                % uisetcolor([0.6 0.8 1])
                % pop.PSTH_perpos_colors =   [autumn(6);winter(6)] ;
                % pop.PSTH_summary_colors=   [autumn(6);winter(6)] ;
                hemifield_indexes       = (real(stm_val(stm_idx))>0)+1;
                con_for_trial_crit      = con_for_line;
                
            case 'StimulusType_Difficulty_Position_ErrorVsCorrect_lastVers'
                [diff_values,~,diff_idx]      =unique(difficulty+3*success');
                
                fig_title               = 'StimulusType_Difficulty_Position';
                con_for_figure          = Styp_idx';
                val_for_figure          = num2cell(Styp_values);
                
                con_for_line            = diff_idx';
                pop.line_labels        =   {'EDiff','EEasy','ETar','CDiff','CEasy','CTar'};
                
                position_indexes        = stm_idx;
                val_for_sub_assignment  = stm_val(u_stm_idx_idx,:);
                val_for_pos_assignment  = stm_val(u_stm_idx_idx,:);
                sub_title               = 'stimulus position';
                pop.PSTH_perpos_colors =   [summer(6);winter(6)] ;
                pop.PSTH_summary_colors=   [summer(6);winter(6)] ;
                hemifield_indexes       = (real(stm_val(stm_idx))>0)+1;
                con_for_trial_crit      = con_for_line;
                
                
            case 'success_in_cue'
                con_for_trial_crit      = suc_idx;
                fig_title               = '';
                val_for_figure          = {shp_values};
                sub_title               = 'cue position';
                val_for_sub_assignment  = cue_val(u_cue_idx_idx,:);
                val_for_pos_assignment  = cue_val(u_cue_idx_idx,:);
                position_indexes        = cue_idx;
                hemifield_indexes       = (real([o.cue_pos]' - [o.fix_pos]')>0)+1;
                
            case 'cue_position'
                con_for_figure          = suc_idx;
                con_for_trial_crit      = shp_idx;
                fig_title               = 'success ';
                sub_title               = 'cue position';
                val_for_figure          = num2cell(suc_values);
                val_for_sub_assignment  = cue_val(u_cue_idx_idx,:);
                val_for_pos_assignment  = cue_val(u_cue_idx_idx,:);
                position_indexes        = cue_idx;
                hemifield_indexes       = (real([o.cue_pos]' - [o.fix_pos]')>0)+1;
                
            case 'hand_choices'
                con_for_trial_crit      = hch_idx;
                fig_title               = '';
                val_for_figure          = {hnd_cho_values};
                
            case 'options'
                con_for_figure          = hnd_idx;
                con_for_trial_crit      = cho_idx;
                fig_title               = 'hand ';
                val_for_figure          = num2cell(hnd_values);
                %
            case 'hands_inactivation'
                con_for_figure          = cho_idx;
                con_for_line            = hpt_idx;
                con_for_trial_crit      = perturbation_group;
                fig_title               = 'choice ';
                val_for_figure          = num2cell(cho_values);
                [~,~,con_for_column]    = unique([hands hemifield_indexes],'rows');
                %         color_idx=ismember(hand_blo_color_combination(:,1),hnd_values) & ismember(hand_blo_color_combination(:,2),blo_values);
                %         pop.PSTH_perpos_colors     = keys.hnd_ptb_colors(ismember(hand_blo_color_combination,hnd_blo_values,'rows'),:);
                %         pop.PSTH_summary_colors    = [keys.hnd_ptb_colors_L(color_idx,:); keys.hnd_ptb_colors_R(color_idx,:)] ;
                %         pop.line_labels            = hand_ptb_labels(ismember(hand_blo_color_combination,hnd_blo_values,'rows'));
                
            case 'hands_inactivation_in_ch'
                con_for_figure          = eff_idx;
                con_for_trial_crit      = chp_idx;
                fig_title               = 'choice_instructed_comp ';
                val_for_figure          = num2cell(eff_values);
                [~,~,con_for_column]    = unique([hands hemifield_indexes],'rows');
                
            case 'hands'
                con_for_figure          = cho_idx;
                con_for_row             = ones(size(con_for_figure));
                con_for_trial_crit      = hef_idx;
                fig_title               = 'choice ';
                sub_title               = 'movement vector ';
                val_for_figure          = num2cell(cho_values);
                con_for_column          = eff_idx;
                %Precision=5;
                
            case 'fixation'
                con_for_row             = eff_idx;
                con_for_trial_crit      = non_idx;
                sub_title               = 'fixation at ';
                val_for_sub_assignment  = fix_val(u_fix_idx_idx,:);
                val_for_pos_assignment  = fix_val(u_fix_idx_idx,:);
                position_indexes        = fix_idx;
                hemifield_indexes       = (real([o.fix_pos])'>0)+1;
                
            case 'movement vectors'
                con_for_trial_crit      = fix_idx;
                sub_title               = 'movement vector ';
                
            case 'target location by origin'
                con_for_trial_crit      = fix_idx;
                sub_title               = 'target position ';
                val_for_sub_assignment  = tar_val(u_tar_idx_idx,:);
                val_for_pos_assignment  = tar_val(u_tar_idx_idx,:);
                position_indexes        = tar_idx;
                hemifield_indexes       = (real([o.tar_pos]')>0)+1;
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
        
        %% trial criterion (very temporary) - needs to be separate for ch and IN, also for hands!...
        pop.position_combinations  = [con_for_figure con_for_trial_crit position_indexes];
        pop.hemifield_combinations = [con_for_figure con_for_trial_crit hemifield_indexes];
        
end

    function [s_c, displacement_types] = center_displacement_working(trial,keys)
        Prec_fix=keys.precision.fixation;
        Prec_pos=keys.precision.position;
        
        movement_direction  =NaN(size(trial'));
        fixation            =NaN(size(trial'));
        target              =NaN(size(trial'));
        cuepos              =NaN(size(trial'));
        stmpos              =NaN(size(trial'));
        
        
        for t=1:numel(trial)
            allpos=trial(t).all_tar_pos;
            correct_tar=trial(t).correct_targets;
            correct_pos=allpos(correct_tar);
            fix_pos = [trial.fix_pos];
            if correct_pos == fix_pos(1) %% distractor only or double_distractor
                correct_pos=allpos(allpos~=0);
            end
            trial(t).stm_pos=correct_pos(1);
        end
        
        s_a=unique_positions([trial.fix_pos],Precision);
        s_b=unique_positions([trial.tar_pos] - [trial.fix_pos],Precision);
        s_c=unique_positions([trial.tar_pos],Precision);
        s_d=unique_positions([trial.cue_pos],Precision);
        s_e=unique_positions([trial.stm_pos],Precision);
        
        
        s_a=unique_positions([trial.fix_pos],Prec_fix);
        s_b=unique_positions([trial.tar_pos] - [trial.fix_pos],Prec_pos);
        s_c=unique_positions([trial.tar_pos],Prec_fix);
        s_d=unique_positions([trial.cue_pos],Prec_pos);
        if isfield(trial,'stm_pos')
            s_e=unique_positions([trial.stm_pos],Prec_pos);
        else
            s_e=s_d;
        end
        for t=1:numel(trial)
            for k=1:numel(s_a)
                if abs(trial(t).fix_pos - s_a(k)) < Prec_fix
                    fixation(t)=s_a(k);
                end
            end
            for k=1:numel(s_b)
                if abs(trial(t).tar_pos - trial(t).fix_pos - s_b(k)) < Prec_pos
                    movement_direction(t)=k;
                end
            end
            for k=1:numel(s_c)
                if abs(trial(t).tar_pos - s_c(k)) < Prec_fix
                    target(t)=s_c(k);
                end
            end
            for k=1:numel(s_d)
                if abs(trial(t).cue_pos - s_d(k)) < Prec_pos
                    cuepos(t)=s_d(k);
                end
            end
            for k=1:numel(s_e)
                if abs(trial(t).stm_pos - s_e(k)) < Precision
                    stmpos(t)=s_e(k);
                end
            end
            for k=1:numel(s_e)
                if abs(trial(t).stm_pos - s_e(k)) < Prec_pos
                    stmpos(t)=s_e(k);
                end
            end
        end
        
        [~,~,unique_condition]      =unique([real(fixation),imag(fixation),real(target),imag(target)],'rows');
        [~,~,fixation_location]     =unique([real(fixation),imag(fixation)],'rows');
        [~,~,movement_direction]    =unique([real(movement_direction),imag(movement_direction)],'rows');
        [~,~,target_location]       =unique([real(target),imag(target)],'rows');
        [~,~,cue_location]          =unique([real(cuepos),imag(cuepos)],'rows');
        [~,~,stimulus_location]     =unique([real(stmpos),imag(stmpos)],'rows');
        
        fix_y=imag(nanmean(fixation));
        fixation=fixation-1i*fix_y;
        target=target-1i*fix_y;
        
        displacement_types=[real(fixation) imag(fixation) real(target) imag(target) real(cuepos) imag(cuepos), real(stmpos) imag(stmpos),...
            unique_condition fixation_location movement_direction target_location, cue_location, stimulus_location]; %
        
        if numel(trial)==0
            displacement_types=NaN(1,14);
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
