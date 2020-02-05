function ph_redo_LR_to_CI(varargin)
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
    
    % keys.anova_table_file=[keys.drive '\Projects\' project '\ephys\' keys.project_version '\tuning_table_combined_CI.mat'];
    keys.tuning_table_foldername=[keys.drive '\Projects\' project '\ephys\' keys.project_version filesep];
    keys.tuning_table_filename='tuning_table_combined';
    load([keys.tuning_table_foldername filesep keys.tuning_table_filename '.mat']);
    keys.tuning_per_unit_table=tuning_per_unit_table;
    
    ph_format_tuning_table(tuning_per_unit_table,keys);
end
end
