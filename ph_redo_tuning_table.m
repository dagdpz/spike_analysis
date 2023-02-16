function ph_redo_tuning_table(population,keys)



if keys.TU.redo_statistics
    %% LS has this line, MP kept it?:
    population = ph_accept_trials_per_unit(population,keys);
    population = ph_assign_perturbation_group(keys,population); %% MP removed this line (hmm... what do we actually store as perutrbation ??)
    population = ph_epochs(population,keys);
    
    
    %% probably need to make sure we do not have contra-ipsi table loaded (?)
    keys.anova_table_file=[keys.basepath_to_save keys.project_version filesep 'tuning_table_combined.mat'];
    if exist(keys.anova_table_file,'file')
        load(keys.anova_table_file)
        keys.tuning_table=tuning_per_unit_table; %% load tuning table
    else
        keys.tuning_table={'unit_ID'};
    end
    %tuning_per_unit_table=ph_ANOVAS(population,keys); % main function
    tuning_per_unit_table=ph_ANOVAS_tmp(population,keys); % main function
    
    keys.tuning_table_foldername=[keys.basepath_to_save keys.project_version filesep];
    keys.tuning_table_filename='tuning_table_combined';
    keys.tuning_table=tuning_per_unit_table;
    save([keys.tuning_table_foldername keys.tuning_table_filename],'tuning_per_unit_table');
    ph_format_tuning_table(tuning_per_unit_table,keys);
    excel_table=ph_format_excel_tuning_table(keys.tuning_table,keys.conditions_to_plot,keys);
    xlswrite([keys.tuning_table_foldername 'tuning_table_' keys.TU.unique_title],excel_table); %not sure why this one would not work...
    
end
    
    excel_table=ph_format_excel_tuning_table(keys.tuning_table,keys.conditions_to_plot,keys);
    xlswrite([keys.basepath_to_save keys.project_version filesep 'tuning_table_' keys.TU.unique_title],excel_table);


end