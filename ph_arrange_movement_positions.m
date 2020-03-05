function [out_analysis]=ph_arrange_movement_positions(o,keys)
o=[o.movement];
out_analysis.movement=o;
% 
% keys.movement_angle_binwidth=60;
% keys.movement_amplitude_binwidth=10;
% keys.movement_plot_type='pre_movement';


%angles_for_title=0:keys.movement_angle_binwidth:359;
angles_for_title=mod((0:keys.movement_angle_binwidth:359)+180,360);
binwidth=keys.movement_angle_binwidth*pi/180;
angular_movement_bins=-pi:binwidth:pi;
amplitudal_movement_bins=[0:keys.movement_amplitude_binwidth:20 Inf];

lower_ang_bin_th=angular_movement_bins(1:end-1) +binwidth/2;
upper_ang_bin_th=angular_movement_bins(2:end)   +binwidth/2;
lower_amp_bin_th=amplitudal_movement_bins(1:end-1);
upper_amp_bin_th=amplitudal_movement_bins(2:end);

endpos_ang  =angle([o.endpos])  +binwidth/2;
strpos_ang  =angle([o.startpos])+binwidth/2;
vector_ang  =angle([o.vector])  +binwidth/2;

endpos_amp  =abs([o.endpos]);
strpos_amp  =abs([o.startpos]);
vector_amp  =abs([o.vector]);


switch keys.movement_arrangement
    case 'revealed_or_not'
%         con_for_figure  = shp_idx;
%         con_for_line    = suc_idx;
        val_for_figure  = num2cell(1);
        val_for_sub_assignment = [real(exp(1i*angular_movement_bins(1:end-1)).') imag(exp(1i*angular_movement_bins(1:end-1)).')];
%        val_for_pos_assignment = cue_val(u_cue_idx_idx,:);
%        fig_title       = 'shape '; 
%        position_indexes = cue_idx;
%        fixation_per_trial = fix_val;
        %color_idx=ismember(hand_choice_color_combination(:,1),hnd_values) & ismember(hand_choice_color_combination(:,2),cho_values);
        out_analysis.PSTH_colors=jet(3);
        out_analysis.PSTH_summary_colors=[1 0 0; 1 1 0; 0 1 0; 0 1 1];
        out_analysis.figure_title_part = 'shape ';
        out_analysis.line_labels={'none','at_start','at_end'};
        
        
        out_analysis.sub_title       = 'movement direction ';  
        revealed_fieldname_part='crossed';
        subplot_fieldname='vector_ang_index';
        line_fieldname='target_revealed_at_start_or_end';
        FR_fieldname='FR_peri';
    case 'pre'
%         con_for_figure  = shp_idx;
%         con_for_line    = suc_idx;
        val_for_figure  = num2cell(1);
        val_for_sub_assignment = [real(exp(1i*angular_movement_bins(1:end-1)).') imag(exp(1i*angular_movement_bins(1:end-1)).')];
%        val_for_pos_assignment = cue_val(u_cue_idx_idx,:);
%        fig_title       = 'shape '; 
%        position_indexes = cue_idx;
%        fixation_per_trial = fix_val;
        %color_idx=ismember(hand_choice_color_combination(:,1),hnd_values) & ismember(hand_choice_color_combination(:,2),cho_values);
        out_analysis.PSTH_colors=jet(numel(amplitudal_movement_bins)-1);
        out_analysis.PSTH_summary_colors=[1 0 0; 1 1 0; 0 1 0; 0 1 1];
        out_analysis.figure_title_part = 'shape ';
        out_analysis.line_labels=cellstr(num2str(upper_amp_bin_th(:)));
        
        
        out_analysis.sub_title       = 'movement direction ';  
        revealed_fieldname_part='at_start';
        subplot_fieldname='vector_ang_index';
        line_fieldname='vector_amp_index';
        FR_fieldname='FR_pre';
    case 'peri'
%         con_for_figure  = shp_idx;
%         con_for_line    = suc_idx;
        val_for_figure  = num2cell(1);
        val_for_sub_assignment = [real(exp(1i*angular_movement_bins(1:end-1)).') imag(exp(1i*angular_movement_bins(1:end-1)).')];
%        val_for_pos_assignment = cue_val(u_cue_idx_idx,:);
%        fig_title       = 'shape '; 
%        position_indexes = cue_idx;
%        fixation_per_trial = fix_val;
        %color_idx=ismember(hand_choice_color_combination(:,1),hnd_values) & ismember(hand_choice_color_combination(:,2),cho_values);
        out_analysis.PSTH_colors=jet(numel(amplitudal_movement_bins)-1);
        out_analysis.PSTH_summary_colors=[1 0 0; 1 1 0; 0 1 0; 0 1 1];
        out_analysis.figure_title_part = 'shape ';
        out_analysis.line_labels=cellstr(num2str(upper_amp_bin_th(:)));
        
        
        out_analysis.sub_title       = 'movement direction ';  
        revealed_fieldname_part='crossed';
        subplot_fieldname='vector_ang_index';
        line_fieldname='vector_amp_index';
        FR_fieldname='FR_peri';
        
    case 'post'
%         con_for_figure  = shp_idx;
%         con_for_line    = suc_idx;
        val_for_figure  = num2cell(1);
        val_for_sub_assignment = [real(exp(1i*angular_movement_bins(1:end-1)).') imag(exp(1i*angular_movement_bins(1:end-1)).')];
%        val_for_pos_assignment = cue_val(u_cue_idx_idx,:);
%        fig_title       = 'shape '; 
%        position_indexes = cue_idx;
%        fixation_per_trial = fix_val;
        %color_idx=ismember(hand_choice_color_combination(:,1),hnd_values) & ismember(hand_choice_color_combination(:,2),cho_values);
        out_analysis.PSTH_colors=jet(numel(amplitudal_movement_bins)-1);
        out_analysis.PSTH_summary_colors=[1 0 0; 1 1 0; 0 1 0; 0 1 1];
        out_analysis.figure_title_part = 'shape ';
        out_analysis.line_labels=cellstr(num2str(upper_amp_bin_th(:)));
        
        
        out_analysis.sub_title       = 'movement direction ';  
        revealed_fieldname_part='at_end';
        subplot_fieldname='vector_ang_index';
        line_fieldname='vector_amp_index';
        FR_fieldname='FR_post';
end

%% subplot positions
[subplot_pos, out_analysis.columns, out_analysis.rows]= DAG_dynamic_positions({val_for_sub_assignment});
[~,temp_sub_index]=sort(subplot_pos);
out_analysis.subplot_titles=angles_for_title(temp_sub_index);
for m= 1:numel(o) 
    out_analysis.movement(m).endpos_ang_index            =find(lower_ang_bin_th<=endpos_ang(m) & upper_ang_bin_th>endpos_ang(m));
    out_analysis.movement(m).stapos_ang_index            =find(lower_ang_bin_th<=strpos_ang(m) & upper_ang_bin_th>strpos_ang(m));
    out_analysis.movement(m).vector_ang_index            =find(lower_ang_bin_th<=vector_ang(m) & upper_ang_bin_th>vector_ang(m));
    
    out_analysis.movement(m).endpos_amp_index            =find(lower_amp_bin_th<=endpos_amp(m) & upper_amp_bin_th>endpos_amp(m));
    out_analysis.movement(m).stapos_amp_index            =find(lower_amp_bin_th<=strpos_amp(m) & upper_amp_bin_th>strpos_amp(m));
    out_analysis.movement(m).vector_amp_index            =find(lower_amp_bin_th<=vector_amp(m) & upper_amp_bin_th>vector_amp(m));
        
    
    %o(m).shape_crossed=o(m).shape_at_start; %%% !!!!!!!!!!!!!!
    out_analysis.movement(m).target_revealed            =o(m).(['tar_' revealed_fieldname_part]);
    out_analysis.movement(m).shape_revealed             =o(m).(['shape_' revealed_fieldname_part]);
    out_analysis.movement(m).no_target_revealed         =(o(m).tar_at_start==0 && o(m).tar_at_end==0 && o(m).tar_crossed==0);
    %out_analysis.movement(m).no_target_revealed         =~(o(m).tar_at_start || o(m).tar_at_end || o(m).tar_crossed);
    
    
    if out_analysis.movement(m).no_target_revealed
    out_analysis.movement(m).target_revealed_at_start_or_end=0;
    elseif  o(m).tar_at_start==1 && o(m).tar_at_end==0
    out_analysis.movement(m).target_revealed_at_start_or_end=1;
    elseif  o(m).tar_at_start==0 && o(m).tar_at_end==1
    out_analysis.movement(m).target_revealed_at_start_or_end=2;
    else
    out_analysis.movement(m).target_revealed_at_start_or_end=NaN;
        
    end
    
    out_analysis.movement(m).subplot_pos                =subplot_pos(out_analysis.movement(m).(subplot_fieldname));
    out_analysis.movement(m).subplot_value              =angular_movement_bins(out_analysis.movement(m).(subplot_fieldname));
    
    out_analysis.movement(m).line                       =out_analysis.movement(m).(line_fieldname);
    out_analysis.movement(m).FR                         =out_analysis.movement(m).(FR_fieldname);
    if isempty(out_analysis.movement(m).subplot_pos)
        out_analysis.movement(m).subplot_pos=NaN;
    end
    if isempty(out_analysis.movement(m).line)
        out_analysis.movement(m).line=NaN;
    end
            
    out_analysis.movement(m).figure                     =1;    
    % which field to put for FR ...?
end

%o.figure_title_part
% 
% %% some preallocations
% out_analysis.figure_title_part      =fig_title;
for n=1:numel(val_for_figure)
    temp     =[arrayfun(@num2str, val_for_figure{n}, 'unif', 0)'; repmat({' '},size(val_for_figure{n},1),1)'];
    out_analysis.figure_title_value{n,1}     =[temp{:}];
end
% %% Assigning each trial
% out_analysis.trial=o;
% for t=1:numel(o)
%     out_analysis.trial(t).line           =con_for_line(t);
%     out_analysis.trial(t).figure         =con_for_figure(t);
%     out_analysis.trial(t).hand           =hands(t);
%     out_analysis.trial(t).choice         =choices(t);
%     out_analysis.trial(t).fixation       =fixation_per_trial(t,:);
%     out_analysis.trial(t).title_part     =sub_title;
%     out_analysis.trial(t).subplot_pos    =subplot_pos(position_indexes(t));
%     out_analysis.trial(t).position       =val_for_pos_assignment(position_indexes(t),:);
%     out_analysis.trial(t).hemifield      =-1*(out_analysis.trial(t).position(1)<0)+1*(out_analysis.trial(t).position(1)>0);
% end

end

function [subplot_pos, n_columns, n_rows]= DAG_dynamic_positions(U_POS)
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

% function [s_c, displacement_types] = center_displacement_working(trial)
% Precision=2;
% 
% movement_direction  =NaN(size(trial'));
% fixation            =NaN(size(trial'));
% target              =NaN(size(trial'));
% cuepos              =NaN(size(trial'));
% 
% s_a=unique_positions([trial.fix_pos],Precision);
% s_b=unique_positions([trial.tar_pos] - [trial.fix_pos],Precision);
% s_c=unique_positions([trial.tar_pos],Precision);
% s_d=unique_positions([trial.cue_pos],Precision);
% 
% for t=1:numel(trial)
%     for k=1:numel(s_a)
%         if abs(trial(t).fix_pos - s_a(k)) < Precision
%             fixation(t)=s_a(k);
%         end
%     end
%     for k=1:numel(s_b)
%         if abs(trial(t).tar_pos - trial(t).fix_pos - s_b(k)) < Precision
%             movement_direction(t)=k;
%         end
%     end
%     for k=1:numel(s_c)
%         if abs(trial(t).tar_pos - s_c(k)) < Precision
%             target(t)=s_c(k);
%         end
%     end    
%     for k=1:numel(s_d)
%         if abs(trial(t).cue_pos - s_d(k)) < Precision
%             cuepos(t)=s_d(k);
%         end
%     end
% end
% 
% [~,~,unique_condition]      =unique([real(fixation),imag(fixation),real(target),imag(target)],'rows');
% [~,~,fixation_location]     =unique([real(fixation),imag(fixation)],'rows');
% [~,~,movement_direction]    =unique([real(movement_direction),imag(movement_direction)],'rows');
% [~,~,target_location]       =unique([real(target),imag(target)],'rows');
% [~,~,cue_location]          =unique([real(cuepos),imag(cuepos)],'rows');
% 
% fix_y=imag(nanmean(fixation));
% fixation=fixation-1i*fix_y;
% target=target-1i*fix_y;
% 
% displacement_types=[real(fixation) imag(fixation) real(target) imag(target) real(cuepos) imag(cuepos), ...
%                     unique_condition fixation_location movement_direction target_location, cue_location];
% end

% function unique_target_positions=unique_positions(all_target_positions,Precision)
% target_positions    =   unique(all_target_positions(~isnan(all_target_positions)));
% n_targets           =   numel(target_positions);
% for t=1:n_targets
%     target_positions(abs(target_positions-target_positions(t))<Precision)=target_positions(t);
% end
% unique_target_positions=unique(target_positions);
% end
