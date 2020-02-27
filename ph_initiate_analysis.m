function ph_initiate_analysis(varargin)
keys=struct;
project=varargin{1};
keys=ph_general_settings(project,keys);
project_specific_settings=[keys.db_folder project filesep 'ph_project_settings.m'];
run(project_specific_settings);

if nargin>1
    keys.project_versions=varargin{2};
end
for f=1:numel(keys.project_versions) % running multiple versions of the same project at once !
if nargin>1
    keys.project_versions=varargin{2};
end
    if ~isempty(keys.project_versions{f})
        keys.project_version=keys.project_versions{f};
    end
    keys.version_specific_settings=[keys.db_folder project filesep keys.project_version filesep 'ph_project_version_settings.m'];
    run(keys.version_specific_settings);
    keys.tuning_table_foldername=[keys.drive '\Projects\' project '\ephys\' keys.project_version filesep];
    keys.tuning_table_filename='tuning_table_combined';
    population=ph_load_population([keys.drive filesep keys.basepath_to_save filesep keys.project_version],['population_']);
    if ~isempty(population); population(arrayfun(@(x) isempty(x.unit_ID),population))=[]; end;
    if ~isempty(population)
        if exist([keys.tuning_table_foldername filesep keys.tuning_table_filename '.mat'],'file')
            load([keys.tuning_table_foldername filesep keys.tuning_table_filename '.mat']);
            keys.tuning_per_unit_table=tuning_per_unit_table;
        else
            keys.tuning_per_unit_table= {'unit_ID'};
        end
        clear tuning_per_unit_table
        population = ph_accept_trials_per_unit(population);   
        population = ph_assign_perturbation_group(keys,population);       
        population = ph_epochs(population,keys);
        tuning_per_unit_table=ph_ANOVAS(population,keys); % main function
        save([keys.tuning_table_foldername filesep keys.tuning_table_filename],'tuning_per_unit_table');
        ph_format_tuning_table(tuning_per_unit_table,keys);
    end
end
end
