function [TT,idx,group_values,unique_group_values]=ph_readout_tuning_table(keys)

PF=keys.normalization_field;

keys.(PF).FR_subtract_baseline=~strcmp(keys.(PF).epoch_BL,'none'); %% this one we should not need here any more !?

TT=keys.tuning_table;
idx.group_parameter=DAG_find_column_index(TT,keys.(PF).group_parameter);
idx.unitID=DAG_find_column_index(TT,'unit_ID');
idx.RF_frame=DAG_find_column_index(TT,keys.(PF).RF_frame_parameter);


group_values=TT(:,idx.group_parameter);
group_values=cellfun(@num2str, group_values, 'UniformOutput', false);
cell_in_any_group=[false; ~ismember(group_values(2:end),keys.(PF).group_excluded)];

TT=TT(cell_in_any_group,:);
group_values=group_values(cell_in_any_group);
unique_group_values=unique(group_values);
if isempty(unique_group_values)
    disp(['group ' keys.(PF).group_parameter ' empty or not found']);
%    return;
end
end