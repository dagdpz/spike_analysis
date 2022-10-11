function ph_initiation(project,versions,keep_version_settings)
% Input format: ph_initiation('Pulv_eye_hand',{'20161212'},1)
% 'Pulv_eye_hand' ... PROJECT FOLDER to save the output on Y:\Projects and
% to get the settings on Dropbox: \Dropbox\DAG\DAG_toolbox\spike_analysis

% ph_initiation is a all in one function: It processes and plots single units (ph_session_processing),
% copies the relevant sorted_neuron excel tables in the project/version folder, saves the used keys per monkey and population data per session in the same folder,
% creates the combined tuning table using ph_ANOVAS, formats it into contra and ipsi format (tuning_table_combined_CI) and creates an excel table with the same information for easier browsing (using ph_format_tuning_table),
% then performs population analysis (ph_initiate_cell_count,ph_initiate_population_analysis, ph_initiate_scatter)
%
% IMPORTANT: all processing options are defined in four different files, and there is certain hierarchy to those files:
% 1) ph_general_settings serves merely as a placeholder, so that keys will never be undefined
% 2) ph_project_settings (in dag_toolbox\spike_analysis\project) for keys that are constant across different analysis versions for this project
% 3) ph_project_version_settings (in dag_toolbox\spike_anlysis\project\version) for keys that are changing in different versions of the  same project
% (4) ph_additional_settings (in dag_toolbox\spike_anlysis\project\version) main usage is to select specific blocks (and its autmatically created when using phys_GUI to create PSTHs)
% Now since the priority is 4>3>2>1, if you are changing keys, make sure to do it in ph_project_version_settings. Once you are happy with your
% settings, you can copy that file to replace ph_project_settings to serve as a basis for a potential second version of your project
% When no ph_project_version_settings are available, a copy of ph_project_settings will be created in the respective subfolder serving as future ph_project_version_settings
% set keep_version_settings (input) to 0 whenever you want to ignore and overwrite ph_project_version_settings
%
% because quite often you might want to change only population analysis related keys,
% you can also run ph_initiate_cell_count,ph_initiate_population_analysis,
% ph_initiate_scatter as well as the ANOVAS (ph_initiate_analysis) independently.

if nargin<3
    keep_version_settings=1;
end
keys=struct;
keys=ph_general_settings(project,keys);
project_specific_settings=[keys.db_folder 'ph_project_settings.m'];
run(project_specific_settings)
if nargin>=2
    keys.project_versions=versions;
end
for f=1:numel(keys.project_versions) % running multiple versions of the same project at once !

    if ~isempty(keys.project_versions{f})
        keys.project_version=keys.project_versions{f};
    end
    if ~exist([keys.db_folder keys.project_version filesep],'dir')
        mkdir([keys.db_folder ], keys.project_version);
    end
    if ~exist([keys.basepath_to_save keys.project_version filesep],'dir')
        mkdir([keys.basepath_to_save ], keys.project_version);
    end
    seed_filename=[keys.basepath_to_save keys.project_version filesep 'seed.mat']; %might not be defined yet
    if exist(seed_filename,'file');
        load(seed_filename);
        rng(seed);
    else
        seed=rng;
        save(seed_filename,'seed');
    end
    keys.version_specific_settings=[keys.db_folder keys.project_version filesep 'ph_project_version_settings.m'];
    if keep_version_settings && exist(keys.version_specific_settings,'file')
        run(keys.version_specific_settings)
        keys.project_versions=versions;
        keys.project_version=keys.project_versions{f};
    else
        delete(keys.version_specific_settings);
        copyfile(project_specific_settings,keys.version_specific_settings);
    end
    keys.additional_settings=[keys.db_folder keys.project_version filesep 'ph_additional_settings.m'];
    if exist(keys.additional_settings,'file')
        run(keys.additional_settings)
    end
    if ~exist([keys.basepath_to_save keys.project_version filesep],'dir')
        mkdir(keys.basepath_to_save, keys.project_version);
    end
    user=getUserName;
    dropboxpath=['C:\Users\' user '\Dropbox'];
    for m=1:numel(keys.batching.monkeys)
        monkey=keys.batching.monkeys{m};
        keys.sorted_neurons_filename =keys.(monkey).sorted_neurons_filename;
        keys.sorted_neurons_sheetname=keys.(monkey).sorted_neurons_sheetname;
        keys.date=keys.(monkey).date;
        keys.monkey                     =[monkey '_phys'];
        keys.sorted_neurons_foldername  =[dropboxpath filesep 'DAG' filesep 'phys' filesep monkey '_phys_dpz'];
        keys.tuning_table_foldername    =[keys.basepath_to_save keys.project_version filesep];
        keys.tuning_table_filename      =['tuning_table_combined'];
        keys.population_foldername      =[keys.basepath_to_save keys.project_version filesep];
        keys.population_filename        =['population_' monkey];
        keys.sites_filename             =['sites_' monkey];
        keys.filelist_formatted         =keys.(monkey).filelist_formatted;
        if exist([keys.tuning_table_foldername keys.tuning_table_filename '.mat'],'file')
            load([keys.tuning_table_foldername keys.tuning_table_filename '.mat']);
            keys.tuning_per_unit_table=tuning_per_unit_table;
        else
            keys.tuning_per_unit_table= {'unit_ID'};
        end
        keys_out                        =ph_session_processing(keys); %% error here for flaff?
        tuning_per_unit_table           =keys_out.tuning_per_unit_table;
         save([keys.basepath_to_save keys.project_version filesep 'keys_' monkey],'keys_out');
         save([keys.tuning_table_foldername keys.tuning_table_filename],'tuning_per_unit_table');
    end
     ph_format_tuning_table(tuning_per_unit_table,keys);
end
 %ph_get_filelist(project,versions);
 ph_initiate_population_analysis(project,versions);
end
