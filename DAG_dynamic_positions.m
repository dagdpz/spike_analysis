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
    heights(n_pos_at_same_height>=floor(size(sin_U_POS,1)/(numel(n_pos_at_same_height)-1)))=[];
    heights_neg=heights(heights<=0);
    heights_pos=heights(heights>=0);

    height_to_adjust=diff(heights_pos)==min(diff(heights_pos));
    for h=find(height_to_adjust)'
        to_adjust=[heights_pos(h),heights_pos(h+1)];
        [~,abs_idx]=sort(abs(to_adjust));
        U_POS_adjusted(U_POS_adjusted(:,2)==to_adjust(abs_idx(2)),2)=to_adjust(abs_idx(1));
    end

    height_to_adjust=diff(heights_neg)==min(diff(heights_neg));
    for h=find(height_to_adjust)'
        to_adjust=[heights_neg(h),heights_neg(h+1)];
        [~,abs_idx]=sort(abs(to_adjust));
        U_POS_adjusted(U_POS_adjusted(:,2)==to_adjust(abs_idx(2)),2)=to_adjust(abs_idx(1));
    end
    
    [heights,~,d]=unique(U_POS_adjusted(:,2));
    prev_n_pos_at_same_height=n_pos_at_same_height;
    n_pos_at_same_height=hist(d,1:max(d));
    if numel(heights_pos)<=1 && numel(heights_neg)<=1
        break;
    end
end

n_columns=max(n_pos_at_same_height);
n_rows=numel(n_pos_at_same_height);

potential_rows_to_adjust_to=find(n_pos_at_same_height==max(n_pos_at_same_height));
row_d_to_mid=abs(potential_rows_to_adjust_to-ceil(n_rows/2));%adjusting to the most middle row!
row_to_adjust_to=potential_rows_to_adjust_to(row_d_to_mid==min(row_d_to_mid));
x_values_to_adjust_to=sort(U_POS_adjusted(U_POS_adjusted(:,2)==heights(row_to_adjust_to(1)),1));

U_POS_final=U_POS_adjusted;
for r=1:n_rows
    x_values_in_this_row=U_POS_adjusted(U_POS_adjusted(:,2)==heights(r),1);
    for c=1:n_pos_at_same_height(r)
        current_pos=U_POS_adjusted(:,2)==heights(r) & U_POS_adjusted(:,1)==x_values_in_this_row(c);
        x_diff=abs(x_values_to_adjust_to-U_POS_adjusted(current_pos,1));
        closest_x=x_values_to_adjust_to(x_diff==min(x_diff));
        U_POS_final(current_pos,1)=closest_x(1);
    end    
end

unique_x_adjusted=unique(U_POS_final(:,1));
unique_y_adjusted=flipud(unique(U_POS_final(:,2)));
subplot_pos=NaN(size(U_POS_final,1),1);

subplot_counter=0;
for r=1:n_rows
    for c=1:n_columns
        subplot_counter=subplot_counter+1;
        current_position_index=U_POS_final(:,1)==unique_x_adjusted(c) & U_POS_final(:,2)==unique_y_adjusted(r);
        if any(current_position_index)
            subplot_pos(current_position_index)=subplot_counter;
        end
    end
end
end