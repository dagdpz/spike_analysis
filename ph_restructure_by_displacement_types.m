function [out_analysis]=ph_restructure_by_displacement_types(o,displacement_types,keys)

%% DEFINITION OF CONDITION INDICES TO PLOT, CURRENTLY TARGET LOCATION, FIXATION LOCATION, MOVEMENT VECTORS, CHOICE, HANDS
                        
all_val=displacement_types(:,1:4);
fix_val=displacement_types(:,1:2);
tar_val=displacement_types(:,3:4);
mov_val=displacement_types(:,3:4) - displacement_types(:,1:2);

all_idx=displacement_types(:,5);
fix_idx=displacement_types(:,6);
mov_idx=displacement_types(:,7);
tar_idx=displacement_types(:,8);
non_idx=ones(size(displacement_types,1),1);

[~,u_all_idx_idx]=unique(all_idx);
[~,u_fix_idx_idx]=unique(fix_idx);
[~,u_mov_idx_idx]=unique(mov_idx);
[~,u_tar_idx_idx]=unique(tar_idx);

choices=[o.choice]';
hands=[o.reach_hand]';
%hands(isnan(hands))=0;
[cho_values,~,cho_idx]=unique(choices);
[hnd_values,~,hnd_idx]=unique(hands);
[hnd_cho_values,~,hch_idx]=unique([hands choices],'rows');
%hnd_cho_values_inverted=hnd_cho_values';
hand_choice_color_combination=[0 0; 0 1; 1 0; 1 1; 2 0; 2 1];
hand_color_combination=[0; NaN; 1; NaN; 2; NaN];
choice_color_combination=[0; 1; NaN; NaN; NaN; NaN];
hand_labels={'NH','NH CH','LH','LH CH','RH','RH CH'};
choice_labels={'IN','CH','','','',''};
hand_choice_labels={'NH IN','NH CH','LH IN','LH CH','RH IN','RH CH'};
fix_labels={'1','2','3','4','5','6','7','8','9'};

switch keys.plot_type 
    case 'hand_choices'
        con_for_figure  = non_idx;
        %con_for_subplot = mov_idx;
        con_for_line    = hch_idx;
        val_for_figure  = {hnd_cho_values};
        %val_for_subplot = mov_val(u_mov_idx_idx,:);
        val_for_sub_assignment = mov_val(u_mov_idx_idx,:);
        val_for_pos_assignment = mov_val(u_mov_idx_idx,:);
        fig_title       = '';
        sub_title       = 'movement vector ';   
        position_indexes = mov_idx;
        %values_per_position = mov_val(u_mov_idx_idx,:);
        position_per_trial = mov_val;
        fixation_per_trial = fix_val;
        color_idx=ismember(hand_choice_color_combination(:,1),hnd_values) & ismember(hand_choice_color_combination(:,2),cho_values);
        out_analysis.PSTH_colors=keys.hnd_choice_colors(ismember(hand_choice_color_combination,hnd_cho_values,'rows'),:);
        out_analysis.PSTH_summary_colors=[keys.hnd_choice_colors_L(color_idx,:); keys.hnd_choice_colors_R(color_idx,:)] ;
        out_analysis.line_labels=hand_choice_labels(ismember(hand_choice_color_combination,hnd_cho_values,'rows'));
    case 'options'
        con_for_figure  = hnd_idx;
        %con_for_subplot = mov_idx;
        con_for_line    = cho_idx;
        val_for_figure  = num2cell(hnd_values);
        %val_for_subplot = mov_val(u_mov_idx_idx,:);
        val_for_sub_assignment = mov_val(u_mov_idx_idx,:);
        val_for_pos_assignment = mov_val(u_mov_idx_idx,:);
        fig_title       = 'hand ';
        sub_title       = 'movement vector ';   
        position_indexes = mov_idx;
        %values_per_position = mov_val(u_mov_idx_idx,:);
        position_per_trial = mov_val;
        fixation_per_trial = fix_val;
        color_idx=ismember(choice_color_combination,cho_values,'rows');
        out_analysis.PSTH_colors=keys.hnd_choice_colors(color_idx,:);
        out_analysis.PSTH_summary_colors=[keys.hnd_choice_colors_L(color_idx,:); keys.hnd_choice_colors_R(color_idx,:)] ;
        out_analysis.line_labels=choice_labels(ismember(choice_color_combination,cho_values));
        
    case 'hands'
        con_for_figure  = cho_idx;
        %con_for_subplot = mov_idx;
        con_for_line    = hnd_idx;
        val_for_figure  = num2cell(cho_values);
        val_for_sub_assignment = mov_val(u_mov_idx_idx,:);
        val_for_pos_assignment = mov_val(u_mov_idx_idx,:);
        fig_title       = 'choice ';
        sub_title       = 'movement vector ';   
        position_indexes = mov_idx;
        %position_per_trial = mov_val(ismember(con_for_subplot,u_mov_idx),:);
        position_per_trial = mov_val;
        fixation_per_trial = fix_val;
        color_idx=ismember(hand_color_combination,hnd_values,'rows');
        out_analysis.PSTH_colors=keys.hnd_choice_colors(color_idx,:);
        out_analysis.PSTH_summary_colors=[keys.hnd_choice_colors_L(color_idx,:); keys.hnd_choice_colors_R(color_idx,:)] ;
        out_analysis.line_labels=hand_labels(ismember(hand_color_combination,hnd_values));
        
    case 'fixation'
        con_for_figure  = non_idx;
        %con_for_subplot = fix_idx;
        con_for_line    = non_idx;
        val_for_figure  = {NaN};
        val_for_sub_assignment = fix_val(u_fix_idx_idx,:);
        val_for_pos_assignment = fix_val(u_fix_idx_idx,:);
        fig_title       = '';
        sub_title       = 'fixation at ';     
        position_indexes = fix_idx;
        %position_per_trial = fix_val(ismember(con_for_subplot,u_fix_idx),:);
        position_per_trial = fix_val;
        fixation_per_trial = fix_val;
        color_idx=ismember(choice_color_combination,cho_values,'rows');
        out_analysis.PSTH_colors=keys.offset_colors;
        out_analysis.PSTH_summary_colors=[keys.hnd_choice_colors_L(color_idx,:); keys.hnd_choice_colors_R(color_idx,:)] ;
        out_analysis.line_labels=fix_labels;
    case 'movement vectors'
        con_for_figure  = non_idx;
        %con_for_subplot = mov_idx;
        con_for_line    = fix_idx;
        val_for_figure  = {NaN};
        val_for_sub_assignment = mov_val(u_mov_idx_idx,:);
        val_for_pos_assignment = mov_val(u_mov_idx_idx,:);
        fig_title       = '';
        sub_title       = 'movement vector ';   
        position_indexes = mov_idx;
        %position_per_trial = mov_val(ismember(con_for_subplot,u_mov_idx),:);
        position_per_trial = mov_val;
        fixation_per_trial = fix_val;
        color_idx=ismember(choice_color_combination,cho_values,'rows');
        out_analysis.PSTH_colors=keys.offset_colors;
        out_analysis.PSTH_summary_colors=[keys.hnd_choice_colors_L(color_idx,:); keys.hnd_choice_colors_R(color_idx,:)] ;
        out_analysis.line_labels=fix_labels;
        
    case 'target location by origin'
        con_for_figure  = non_idx;
        %con_for_subplot = tar_idx;
        con_for_line    = fix_idx;
        val_for_figure  = {NaN};
        val_for_sub_assignment = tar_val(u_tar_idx_idx,:);
        val_for_pos_assignment = tar_val(u_tar_idx_idx,:);
        fig_title       = '';
        sub_title       = 'target position ';   
        position_indexes = tar_idx;
        %position_per_trial = tar_val(ismember(con_for_subplot,u_tar_idx),:);
        position_per_trial = tar_val;
        fixation_per_trial = fix_val;
        color_idx=ismember(choice_color_combination,cho_values,'rows');
        out_analysis.PSTH_colors=keys.offset_colors;
        out_analysis.PSTH_summary_colors=[keys.hnd_choice_colors_L(color_idx,:); keys.hnd_choice_colors_R(color_idx,:)] ;
        out_analysis.line_labels=fix_labels;
        
    case 'all conditions'
        con_for_figure  = fix_idx;
        %con_for_subplot = mov_idx;
        con_for_line    = non_idx;
        val_for_figure  = num2cell(fix_val(u_fix_idx_idx,:));
        val_for_sub_assignment = mov_val(u_mov_idx_idx,:);
        val_for_pos_assignment = all_val(u_all_idx_idx,:);
        fig_title       = 'fixation at ';
        sub_title       = 'movement vector ';   
        position_indexes = mov_idx; %all_idx;
        %position_per_trial = mov_val(ismember(con_for_subplot,u_mov_idx),:);
        position_per_trial = mov_val;
        fixation_per_trial = fix_val;
        color_idx=ismember(choice_color_combination,cho_values,'rows');
        out_analysis.PSTH_colors=keys.offset_colors;
        out_analysis.PSTH_summary_colors=[keys.hnd_choice_colors_L(color_idx,:); keys.hnd_choice_colors_R(color_idx,:)] ;
        out_analysis.line_labels=fix_labels;
end

%% subplot positions
[subplot_pos, out_analysis.columns, out_analysis.rows]= dynamic_positions({val_for_sub_assignment});

L_idx=position_per_trial(:,1)<0;
R_idx=position_per_trial(:,1)>0;
D_idx=position_per_trial(:,2)<0;
U_idx=position_per_trial(:,2)>0;

%% some preallocations
n_trials_P=zeros(numel(val_for_sub_assignment),1); % for positions
n_trials_L=0;
n_trials_R=0;
n_trials_U=0;
n_trials_D=0;

out_analysis.left                   =struct();
out_analysis.right                  =struct();
out_analysis.up                     =struct();
out_analysis.down                   =struct();
out_analysis.position           	=struct();
out_analysis.figure_title_part      =fig_title;
for n=1:numel(val_for_figure)
    temp     =[arrayfun(@num2str, val_for_figure{n}, 'unif', 0)'; repmat({' '},size(val_for_figure{n},1),1)'];
    out_analysis.figure_title_value{n,1}     =[temp{:}];
end
%% Assigning each trial
for t=1:numel(o)
    fig  = con_for_figure(t);
    %sub  = subplot_pos(con_for_subplot(t));
    pos  = subplot_pos(position_indexes(t));
    lin  = con_for_line(t);
    hand = hands(t);
    choice = choices(t);
    o(t).line=lin;
    o(t).figure=fig;
    o(t).hand=hand;
    o(t).choice=choice; 
    [o(t).spikes_per_state.line]=deal(lin);
    [o(t).spikes_per_state.figure]=deal(fig);   
    [o(t).spikes_per_state.hand]=deal(hand);
    [o(t).spikes_per_state.choice]=deal(choice);   
    [o(t).spikes_per_state.fixation]=deal(fixation_per_trial(t,:));  
    
    n_trials_P(pos)                                      =n_trials_P(pos)+1;
    out_analysis.position(pos).trial(n_trials_P(pos))    =o(t);
    out_analysis.position(pos).position                  =val_for_pos_assignment(position_indexes(t),:);
    out_analysis.position(pos).title_part                =sub_title;
    
    if L_idx(t)
        n_trials_L=n_trials_L+1;
        out_analysis.left.trial(n_trials_L)=o(t);
    elseif R_idx(t)
        n_trials_R=n_trials_R+1;
        out_analysis.right.trial(n_trials_R)=o(t);
    end
    
    if U_idx(t)
        n_trials_U=n_trials_U+1;
        out_analysis.up.trial(n_trials_U)=o(t);
    elseif D_idx(t)
        n_trials_D=n_trials_D+1;
        out_analysis.down.trial(n_trials_D)=o(t);
    end
end

end

function [subplot_pos, n_columns, n_rows]= dynamic_positions(U_POS)
sin_U_POS=round(vertcat(U_POS{:})*100)/100;
U_POS_adjusted=sin_U_POS;
[heights,~,d]=unique(U_POS_adjusted(:,2));
prev_n_pos_at_same_height=[];
n_pos_at_same_height=hist(d,1:max(d));
while max(n_pos_at_same_height)*(numel(n_pos_at_same_height)-1)>size(sin_U_POS,1)
    if numel(prev_n_pos_at_same_height)== numel(n_pos_at_same_height) && all(prev_n_pos_at_same_height==n_pos_at_same_height)
        break;
    end
    heights(n_pos_at_same_height>=floor(sqrt(size(sin_U_POS,1))))=[];
    heights_neg=heights(heights<=0);
    heights_pos=heights(heights>=0);
    heights_0=sum(~any(heights==0));
    height_to_adjust=[diff(heights_neg)==min(diff(heights_neg)); false(heights_0); diff(heights_pos)==min(diff(heights_pos))];
    for h=find(height_to_adjust)'
        to_adjust=[heights(h),heights(h+1)];
        [~,abs_idx]=sort(abs(to_adjust));
        U_POS_adjusted(U_POS_adjusted(:,2)==to_adjust(abs_idx(2)),2)=to_adjust(abs_idx(1));
    end
    [heights,~,d]=unique(U_POS_adjusted(:,2));
    prev_n_pos_at_same_height=n_pos_at_same_height;
    n_pos_at_same_height=hist(d,1:max(d));
end

n_columns=max(n_pos_at_same_height);
n_rows=numel(n_pos_at_same_height);

potential_rows_to_adjust_to=find(n_pos_at_same_height==max(n_pos_at_same_height));
row_d_to_mid=abs(potential_rows_to_adjust_to-ceil(n_rows/2));%adjusting to the most middle row!
row_to_adjust_to=potential_rows_to_adjust_to(row_d_to_mid==min(row_d_to_mid));
x_values_to_adjust_to=sort(U_POS_adjusted(U_POS_adjusted(:,2)==heights(row_to_adjust_to(1)),1));

for r=1:n_rows
    x_values_in_this_row=U_POS_adjusted(U_POS_adjusted(:,2)==heights(r),1);
    for c=1:n_pos_at_same_height(r)
        current_pos=U_POS_adjusted(:,2)==heights(r) & U_POS_adjusted(:,1)==x_values_in_this_row(c);
        x_diff=abs(x_values_to_adjust_to-U_POS_adjusted(current_pos,1));
        closest_x=x_values_to_adjust_to(x_diff==min(x_diff));
        U_POS_adjusted(current_pos,1)=closest_x(1);
    end    
end

unique_x_adjusted=unique(U_POS_adjusted(:,1));
unique_y_adjusted=flipud(unique(U_POS_adjusted(:,2)));
subplot_pos=NaN(size(U_POS_adjusted,1),1);

subplot_counter=0;
for r=1:n_rows
    for c=1:n_columns
        subplot_counter=subplot_counter+1;
        current_position_index=U_POS_adjusted(:,1)==unique_x_adjusted(c) & U_POS_adjusted(:,2)==unique_y_adjusted(r);
        if any(current_position_index)
            subplot_pos(current_position_index)=subplot_counter;
        end
    end
end
end

