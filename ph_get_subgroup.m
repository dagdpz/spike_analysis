function ph_get_subgroup(anova_table_file)



keys=struct;
keys=ph_general_settings(project,keys);
project_specific_settings=[keys.db_folder project filesep 'ph_project_settings.m'];
run(project_specific_settings)
















keys.anova_table_file=anova_table_file;

keys.monkey='';
keys.labels.choices=
keys.labels.handsIC=
keys.labels.perturbations=

keys.tt.tasktypes
keys.tt.hands
keys.tt.choices
keys.tt.perturbations

keys.tt.IC_for_criterion
keys.cal.min_trials_per_condition
keys.tt.trial_criterion_in
keys.tt.trial_criterion_ch



keys.tt.selection
keys.tt.unselect


keys.tt.combine_tuning_properties

tuning_per_unit_table=ph_load_extended_tuning_table(keys);
ph_reduce_tuning_table(tuning_per_unit_table,keys);







end