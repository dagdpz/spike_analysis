function ph_initiation(projects,versions,replace_version_specific_settings)
if nargin<1
    projects={'Pulv_eye_hand'};%'PPC_pulv_eye_hand'; %'Pulv_microstim_behavior'; %'Pulv_eye_hand'; % 'STS_memory_saccades'};
end
% if nargin<3
%     replace_version_specific_settings=1;
% end
keys=struct;
for p=1:numel(projects)
    project=projects{p};
    keys=ph_general_settings(project,keys);
    project_specific_settings=[keys.db_folder project filesep 'phys_settings.m'];
    run(project_specific_settings)
    if nargin>=2
        keys.pdf_folders=versions;
    end
    
    for f=1:numel(keys.pdf_folders) % running multiple versions of the same project at once !
        if ~isempty(keys.pdf_folders{f})
            keys.pdf_folder=keys.pdf_folders{f};
        end
        if ~exist([keys.db_folder project filesep keys.pdf_folder filesep],'dir')
            mkdir([keys.db_folder project], keys.pdf_folder);
        end
        keys.version_specific_settings=[keys.db_folder project filesep keys.pdf_folder filesep 'phys_settings.m'];
        
        %if replace_version_specific_settings
            delete(keys.version_specific_settings);
            copyfile(project_specific_settings,keys.version_specific_settings);
%         else
%             run(keys.version_specific_settings)
%         end
        keys.additional_settings=[keys.db_folder project filesep keys.pdf_folder filesep 'additional_phys_settings.m'];
        if exist(keys.additional_settings,'file')
            run(keys.additional_settings)
        end
        if ~exist([keys.drive filesep keys.basepath_to_save filesep keys.pdf_folder filesep],'dir')
            mkdir([keys.drive filesep keys.basepath_to_save], keys.pdf_folder);
        end
        user=getUserName;
        dropboxpath=['C:\Users\' user '\Dropbox'];
        for m=1:numel(keys.monkeys)
            monkey=keys.monkeys{m};
            keys.sorted_neurons_filename =keys.(monkey).sorted_neurons_filename;
            keys.date=keys.(monkey).date;
            keys.monkey                     =[monkey '_phys'];
            keys.sorted_neurons_foldername  =[dropboxpath '\DAG\phys\' monkey '_phys_dpz'];
            keys.tuning_table_foldername    =[keys.drive filesep keys.basepath_to_save filesep keys.pdf_folder];
            keys.population_foldername      =[keys.drive filesep keys.basepath_to_save filesep keys.pdf_folder];
            keys.population_filename        =['population_' monkey];
            keys.tuning_table_filename      =['tuning_table_combined'];
            keys.filelist_formatted                         = keys.(monkey).filelist_formatted;
            [~, keys_out]                   =ph_multi_summary(keys); %% error here for flaff?
            save([keys.drive filesep keys.basepath_to_save filesep keys.pdf_folder filesep 'keys_' monkey],'keys_out');
        end
    end
end


ph_initiate_analysis(projects)
ph_initiate_cell_count(projects)
ph_initiate_population_analysis(projects)
end
