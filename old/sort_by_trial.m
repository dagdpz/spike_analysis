function Trial = sort_by_trial(data_per_run)
fieldnames_to_remove={'time_axis','x_eye','y_eye','x_hnd','y_hnd','unit','channel','TDT_CAP1','TDT_ECG1','TDT_ECG4','TDT_POX1'};
    
Trial=[data_per_run.trial];
Trial=rmfield(Trial,fieldnames_to_remove(ismember(fieldnames_to_remove,fieldnames(Trial))));
end