function ph_add_Subregions_to_table(folder,filename,Subregions)

load([folder filesep filename '.mat']);
% Subregions{1}{1}.monkey='Fla';
% Subregions{1}{1}.target='dPulv_l';
% Subregions{1}{1}.grid_x=-11;
% Subregions{1}{1}.grid_y=1;
% Subregions{1}{1}.z_min=45;
% Subregions{1}{1}.z_max=49;
%
% Subregions{1}{2}.monkey='Lin';
% Subregions{1}{2}.target='dPulv_l';
% Subregions{1}{2}.grid_x=-5;
% Subregions{1}{2}.grid_y=4;
% Subregions{1}{2}.z_min=20;
% Subregions{1}{2}.z_max=49;
%
% Subregions{2}{1}.monkey='Fla';
% Subregions{2}{1}.target='dPulv_l';
% Subregions{2}{1}.grid_x=-11;
% Subregions{2}{1}.grid_y=1;
% Subregions{2}{1}.z_min=40;
% Subregions{2}{1}.z_max=45;

idx_ID          =find_column_index(tuning_per_unit_table,'unit_ID');
idx_target      =find_column_index(tuning_per_unit_table,'target');
idx_x           =find_column_index(tuning_per_unit_table,'grid_x');
idx_y           =find_column_index(tuning_per_unit_table,'grid_y');
idx_z           =find_column_index(tuning_per_unit_table,'electrode_depth');
idx_subregion   =find_column_index(tuning_per_unit_table,'Subregion');
if isempty(idx_subregion)
    idx_subregion   =size(tuning_per_unit_table,2)+1;
end
tuning_per_unit_table{1,idx_subregion}='Subregion';
for r=1:numel(Subregions)
    for h=1:numel(Subregions{r})
        i1=~cellfun(@isempty,strfind(tuning_per_unit_table(2:end,idx_ID),Subregions{r}{h}.monkey));
        i2=~cellfun(@isempty,strfind(tuning_per_unit_table(2:end,idx_target),Subregions{r}{h}.target));
        i3=vertcat(tuning_per_unit_table{2:end,idx_x})==Subregions{r}{h}.grid_x;
        i4=vertcat(tuning_per_unit_table{2:end,idx_y})==Subregions{r}{h}.grid_y;
        i5=vertcat(tuning_per_unit_table{2:end,idx_z})>Subregions{r}{h}.z_min & vertcat(tuning_per_unit_table{2:end,idx_z})<=Subregions{r}{h}.z_max;
        index=[false; i1&i2&i3&i4&i5];
        tuning_per_unit_table(index,idx_subregion)={r};
    end
end
save([folder filesep filename '.mat'],'tuning_per_unit_table');
end