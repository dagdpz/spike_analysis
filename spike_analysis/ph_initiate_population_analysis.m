function ph_initiate_population_analysis(varargin)
if nargin==0
    projects={'PPC_pulv_eye_hand','Pulv_eye_hand','PPC_pulv_eye_hand','STS_memory_saccades','Pulv_microstim_behavior','Pulv_eye_gaze_position'};
else
    projects=varargin{1};
end
keys=struct;
keys.instructed_choice='in';
for p=1:numel(projects)
    project=projects{p};
    keys=ph_general_settings(project,keys);
    project_specific_settings=[keys.db_folder filesep project filesep 'phys_settings.m'];
    run(project_specific_settings);
    for f=1:numel(keys.pdf_folders) % running multiple versions of the same project at once !
        if ~isempty(keys.pdf_folders{f})
            keys.pdf_folder=keys.pdf_folders{f};
        end
        pdf_folder=keys.pdf_folder;
        keys.version_specific_settings=[keys.db_folder project filesep keys.pdf_folder filesep 'phys_settings.m'];
        run(keys.version_specific_settings);
        keys.pdf_folder=pdf_folder;
        if keys.combine_monkeys
            keys.monkeys={''};
        end
        
        for m=1:numel(keys.monkeys)
            keys.monkey=keys.monkeys{m};
            keys.anova_table_file=[keys.drive '\Projects\' project '\ephys\' keys.pdf_folder '\tuning_table_combined_CI.mat'];
            population=load_population([keys.drive filesep keys.basepath_to_save filesep keys.pdf_folder],['population_' keys.monkey]);
            for t=1:numel(keys.targets)
                target=keys.targets{t};
                for cc=1:numel(keys.pop)
                    
                    %% keys assignment (!)                   
                    
                    keys.FR_subtract_baseline               = assign_if_not_empty(keys.pop(cc).FR_subtract_baseline,0);
                    keys.epoch_BL                           = assign_if_not_empty(keys.pop(cc).epoch_BL,'INI'); 
                    keys.population_normalization           = assign_if_not_empty(keys.pop(cc).normalization,'none');
                    keys.epoch_for_normalization            = assign_if_not_empty(keys.pop(cc).epoch_for_normalization,'');
                    keys.plot_RF                            = assign_if_not_empty(keys.pop(cc).plot_RF,0);
                    keys.epoch_RF     = assign_if_not_empty(keys.pop(cc).epoch_RF,'');
                    keys.population_group_excluded          = assign_if_not_empty(keys.pop(cc).group_excluded,{'','-'});
                    keys.population_ylim                    = assign_if_not_empty(keys.pop(cc).y_lim,[]);
                    keys.tt.only_both_hands                 = assign_if_not_empty(keys.pop(cc).only_both_hands,0);
                    keys.tt.only_overlapping_tasktypes      = assign_if_not_empty(keys.pop(cc).only_overlapping_tasktypes,0);
                    
                    keys.tt.unselect                        = keys.pop(cc).unselect;
                    keys.tt.combine_tuning_properties       = keys.pop(cc).combine_tuning_properties;
                    keys.population_group_parameter         = keys.pop(cc).group_parameter;
                    conditions_to_plot                      = keys.pop(cc).conditions_to_plot;
                    population_selection                    = {'target',target};
                    population_selection                    = [population_selection ; keys.pop(cc).Selection];
                    for subregion=1:keys.n_Subregions
                        if keys.Subregions_separately
                            keys.tt.selection               = [population_selection ; {'Subregion', subregion}];
                        else
                            keys.tt.selection               = population_selection;
                        end
                        for condition_to_plot=conditions_to_plot
                            %%temp
                            c=1;
                             keys.case={keys.case_summaries{c}(1:3)};
                            keys.conditions_to_plot=condition_to_plot;
                            keys.tt.tasktypes=strcat(condition_to_plot, ['_' keys.case{1}]);
                            ph_population_analysis_contra_ipsi(population,keys);
                        end
                    end
                end
            end
        end
    end
end
end

function output=assign_if_not_empty(input,default)
    if ~isempty(input)
        output                       = input;
    else
        output                       = default;        
    end
end