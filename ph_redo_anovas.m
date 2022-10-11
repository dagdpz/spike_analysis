function ph_redo_anovas(varargin)
% 3rd input for redo_ANOVAS
keys=struct;
project=varargin{1};
keys=ph_general_settings(project,keys);
project_specific_settings=[keys.db_folder  'ph_project_settings.m'];
run(project_specific_settings);
redo_anovas=1;

if nargin>1
    keys.project_versions=varargin{2};
end
if nargin>2
    redo_anovas=varargin{3};
end
for f=1:numel(keys.project_versions) % running multiple versions of the same project at once !
    if ~isempty(keys.project_versions{f})
        keys.project_version=keys.project_versions{f};
    end
    keys.version_specific_settings=[keys.db_folder keys.project_version filesep 'ph_project_version_settings.m'];
    run(keys.version_specific_settings);
    keys.tuning_table_foldername=[keys.basepath_to_save keys.project_version filesep];
    keys.tuning_table_filename='tuning_table_combined';
    
    if exist([keys.tuning_table_foldername filesep keys.tuning_table_filename '.mat'],'file')
        load([keys.tuning_table_foldername filesep keys.tuning_table_filename '.mat']);
        keys.tuning_per_unit_table=tuning_per_unit_table;
    else
        keys.tuning_per_unit_table= {'unit_ID'};
    end
    
    if redo_anovas
        population=ph_load_population([keys.basepath_to_save keys.project_version],['population_']);
        if ~isempty(population); population(arrayfun(@(x) isempty(x.unit_ID),population))=[]; end;
        if ~isempty(population)
            seed_filename=[keys.basepath_to_save keys.project_version filesep 'seed.mat'];
            if exist(seed_filename,'file');
                load(seed_filename);
                rng(seed);
            else
                seed=rng;
                save(seed_filename,'seed');
            end
            
            clear tuning_per_unit_table
            %% LS has this line, MP kept it?:
            population = ph_accept_trials_per_unit(population,keys); 
            population = ph_assign_perturbation_group(keys,population); %% MP removed this line (hmm... what do we actually store as perutrbation ??)
            
            population = ph_epochs(population,keys);
            tuning_per_unit_table=ph_ANOVAS(population,keys); % main function
            save([keys.tuning_table_foldername keys.tuning_table_filename],'tuning_per_unit_table');
        end
    end
    
    ph_format_tuning_table(tuning_per_unit_table,keys);
end
end
