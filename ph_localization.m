function ph_localization(keys)
% xyz: grid_hole_x, grid_hole_y, z depth mm from "top of the brain" or "from chamber top"

FN=fieldnames(keys.LO);
for f=1:numel(FN)
keys.(FN{f})=keys.LO.(FN{f});
end
TT=keys.tuning_table;

co = {'b','m','r'}; % visual,visuomotor, motor
switch keys.significance_to_plot
    case 'ungrouped'
        task_column=DAG_find_column_index(TT,'ungrouped');
        sign_column=DAG_find_column_index(TT,'ungrouped');
        co = {'b','m','r'}; % visual,visuomotor, motor
    case 'memory'
        task_column=DAG_find_column_index(TT,'in_epoch_main_Msac_fix');
        sign_column=DAG_find_column_index(TT,'in_NH_Cue_position_Msac_opt');
        co = {'b','m','r'}; % visual,visuomotor, motor
    case 'direct'
        task_column=DAG_find_column_index(TT,'in_epoch_main_Msac_fix');
        sign_column=DAG_find_column_index(TT,'in_NH_Cue_position_Vsac_opt');
        co = {'b','m','r'}; % visual,visuomotor, motor
    case 'Ddre'
        task_column=DAG_find_column_index(TT,'in_epoch_main_Msac_fix');
        sign_column=DAG_find_column_index(TT,'in_NH_Cue_position_Ddre_han');
        co = {'b','m','r'}; % visual,visuomotor, motor
    case 'Ddre_han'
        task_column=DAG_find_column_index(TT,'in_epoch_main_Ddre_han');
        sign_column=DAG_find_column_index(TT,'in_AH_Cue_position_Ddre_han');
        co = {'b','m','r'}; % visual,visuomotor, motor
    case 'visuomotor'
        task_column=DAG_find_column_index(TT,'in_epoch_main_Msac_opt');
        sign_column=[DAG_find_column_index(TT,'visual_Msac_opt'),...
            DAG_find_column_index(TT,'visuomotor_Msac_opt'),...
            DAG_find_column_index(TT,'motor_Msac_opt')];
        co = {'b','m','r'}; % visual,visuomotor, motor
    case 'gaze'
        task_column=DAG_find_column_index(TT,'in_epoch_main_Msac_fix');
        coloumn_temp=DAG_find_column_index(TT,'in_AH_Fhol_gaze_modulation_x_Msac_fix');
        coloumn_temp2=DAG_find_column_index(TT,'in_AH_Fhol_gaze_pref_x_Msac_fix');
        coloumn_temp3=DAG_find_column_index(TT,'in_AH_Fhol_position_Msac_fix');
        temp_table=[TT(:,coloumn_temp) TT(:,coloumn_temp2) TT(:,coloumn_temp3)];
        for k=1:size(temp_table,1)
            gaze_column{k}=[temp_table{k,:}];
        end
        gaze_column=gaze_column';
        n_columns=size(TT,2);
        TT(:,n_columns+1)=num2cell(strcmp(gaze_column,'monotonousCStrue'));
        TT(:,n_columns+2)=num2cell(strcmp(gaze_column,'monotonousIStrue'));
        TT(:,n_columns+3)=num2cell(strcmp(gaze_column,'nonmonotonousCEtrue'));
        TT(:,n_columns+4)=num2cell(strcmp(gaze_column,'nonmonotonousPEtrue'));
        sign_column=[n_columns+1,n_columns+2,n_columns+3,n_columns+4];
        co = {'m','r','b','g'}; % mono contra, mono ipsi, nonmonocentral nonmonoperipheral
end

clear penetration_date xyz target notes

target              = TT(:,DAG_find_column_index(TT,'target'));
penetration_date    = TT(:,DAG_find_column_index(TT,'unit_ID'));
task_type           = TT(:,task_column);
for c=1:numel(sign_column)
    significant_c{c}       = TT(:,sign_column(c));
end

row_index           = cellfun(@(x) ~isempty(strfind(x,keys.target_area)),target) & cellfun(@(x) ~isempty(strfind(x,keys.monkey)),penetration_date) & ~cellfun(@isempty,task_type); %% index by target location and monkey initials

idx_x=DAG_find_column_index(TT,'grid_x');
idx_y=DAG_find_column_index(TT,'grid_y');
idx_z=DAG_find_column_index(TT,'electrode_depth');
xyz                 = cell2mat(TT(row_index,[idx_x idx_y idx_z]));
xyz(:,3)            = -xyz(:,3);
xyz_nojitter        = xyz;

switch keys.saggital_or_coronal
    case 'coronal'
        xyz(:,1)            = xyz(:,1) + (rand(size(xyz(:,1)))-0.5)*1.5; % jitter
    case 'sagittal'
        xyz(:,2)            = xyz(:,2) + (rand(size(xyz(:,2)))-0.5)*1.5; % jitter
end
penetration_date    = penetration_date(row_index);
for c=1:numel(sign_column)
    significant(:,c) = cellfun(@(x) (~isa(x,'char')&&x==1)||(~isempty(x)&&isa(x,'char')&&~strcmp(x,'-')&&~strcmp(x,'false')),significant_c{c}(row_index));
end


keys.significant=significant;
keys.xyz = xyz;
keys.xyz_nojitter =xyz_nojitter;
keys.penetration_date = penetration_date;

CL_plot_electrode_localization(keys,keys.significance_to_plot,co,0,keys.saggital_or_coronal)
h =  findobj('type','figure');
for n = 1:length(h);
    export_fig(h(n),[keys.path_to_save  keys.significance_to_plot '_'  keys.monkey '_'  keys.target_area '_' num2str(n)], '-pdf','-transparent')
    close(h(n));
end

end