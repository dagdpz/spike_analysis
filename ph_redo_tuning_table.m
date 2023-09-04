function ph_redo_tuning_table(population,trials,keys)

if keys.TU.redo_statistics
    %% need to make sure we have uncorrected, complete R/L table loaded
    keys.anova_table_file=[keys.basepath_to_save keys.project_version filesep 'tuning_table_combined.mat'];
    if exist(keys.anova_table_file,'file')
        load(keys.anova_table_file)
        keys.tuning_table=tuning_per_unit_table; 
    else
        keys.tuning_table={'unit_ID'};
    end
    tuning_per_unit_table=ph_ANOVAS(population,trials,keys); % main function
    keys.tuning_table_foldername=[keys.basepath_to_save keys.project_version filesep];
    keys.tuning_table_filename='tuning_table_combined';
    keys.tuning_table=tuning_per_unit_table;
    save([keys.tuning_table_foldername keys.tuning_table_filename],'tuning_per_unit_table');
    ph_format_tuning_table(tuning_per_unit_table,keys);
    excel_table=ph_format_excel_tuning_table(keys.tuning_table,keys.conditions_to_plot,keys);
    xlswrite([keys.tuning_table_foldername 'tuning_table_' keys.target '_' keys.TU.unique_title '.xlsx'],excel_table);
else
    excel_table=ph_format_excel_tuning_table(keys.tuning_table,keys.conditions_to_plot(:)',keys);
    xlswrite([keys.basepath_to_save keys.project_version filesep 'tuning_table_' keys.target '_' keys.TU.unique_title '.xlsx'],excel_table);
end
end