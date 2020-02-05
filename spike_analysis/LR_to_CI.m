function tuning_per_unit_table=LR_to_CI(tuning_per_unit_table)
idx_target=find_column_index(tuning_per_unit_table,'target');
Left_hemisphere_targets={'dPulv_l','pdSTS_L','FST_L','pTPO_L','MIP','MIP_L','unknown'};
Right_hemisphere_targets={'dPulv_r'};

R_hem=ismember(tuning_per_unit_table(:,idx_target),Right_hemisphere_targets);
L_hem=ismember(tuning_per_unit_table(:,idx_target),Left_hemisphere_targets);

O_entries={'LS','RS','LH','RH'};
R_entries={'CS','IS','CH','IH'};
L_entries={'IS','CS','IH','CH'};

tuning_table_temp=tuning_per_unit_table(R_hem,:);
for e=1:numel(O_entries)
    tuning_table_temp(cellfun(@(x) strcmp(x,O_entries{e}),tuning_table_temp))=R_entries(e);
end
tuning_per_unit_table(R_hem,:)=tuning_table_temp;

tuning_table_temp=tuning_per_unit_table(L_hem,:);
for e=1:numel(O_entries)
    tuning_table_temp(cellfun(@(x) strcmp(x,O_entries{e}),tuning_table_temp))=L_entries(e);
end
tuning_per_unit_table(L_hem,:)=tuning_table_temp;

% inverting effect sizes for right hemisphere, so R-L becomes C-I
ES_columns=cellfun(@(x) any(strfind(x,'_ES_')),tuning_per_unit_table(1,:));
tuning_per_unit_table(R_hem,ES_columns)= cellfun(@(x) x*-1,tuning_per_unit_table(R_hem,ES_columns),'uniformoutput',false);

IX_columns=cellfun(@(x) any(strfind(x,'_IX_')),tuning_per_unit_table(1,:));
tuning_per_unit_table(R_hem,IX_columns)= cellfun(@(x) x*-1,tuning_per_unit_table(R_hem,IX_columns),'uniformoutput',false);

LH_columns=cellfun(@(x) any(strfind(x,'_LH_')),tuning_per_unit_table(1,:));
RH_columns=cellfun(@(x) any(strfind(x,'_RH_')),tuning_per_unit_table(1,:));
Lhem_LH=tuning_per_unit_table(L_hem,LH_columns);
Lhem_RH=tuning_per_unit_table(L_hem,RH_columns);

%here should be a loop to avoid crazyness
tuning_per_unit_table(L_hem,LH_columns)=Lhem_RH;
tuning_per_unit_table(L_hem,RH_columns)=Lhem_LH;

for c=find(LH_columns)
    table_title=tuning_per_unit_table{1,c};
    position_in_title_to_change=strfind(table_title,'_LH_')+1;
    tuning_per_unit_table{1,c}(position_in_title_to_change)='C';
end
for c=find(RH_columns)
    table_title=tuning_per_unit_table{1,c};
    position_in_title_to_change=strfind(table_title,'_RH_')+1;
    tuning_per_unit_table{1,c}(position_in_title_to_change)='I';
end
% for in_or_ch={'in','ch'}
%     columns_L=find(~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),[in_or_ch{:} '_L_'])));
%     %columns_R=find(~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),[in_or_ch{:} '_R_'])));
%     for c=1:numel(columns_L)
%         idx_L=columns_L(c);
%         idx_R=find(~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),[tuning_per_unit_table{1,idx_L}(1:3) 'R' tuning_per_unit_table{1,idx_L}(5:end)])));
%         tuning_per_unit_table{1,end+1}=[tuning_per_unit_table{1,idx_L}(1:3) 'C' tuning_per_unit_table{1,idx_L}(5:end)];
%         tuning_per_unit_table(R_hem,end)=tuning_per_unit_table(R_hem,idx_L);
%         tuning_per_unit_table(L_hem,end)=tuning_per_unit_table(L_hem,idx_R);
%         tuning_per_unit_table{1,end+1}=[tuning_per_unit_table{1,idx_L}(1:3) 'I' tuning_per_unit_table{1,idx_L}(5:end)];
%         tuning_per_unit_table(R_hem,end)=tuning_per_unit_table(R_hem,idx_R);
%         tuning_per_unit_table(L_hem,end)=tuning_per_unit_table(L_hem,idx_L);
%     end
% end



end