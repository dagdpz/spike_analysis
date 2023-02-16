function ph_run_state_space(varargin)
%ph_run_population_analysis('PPC_pulv_eye_hand',{'LIP_dPul_inj_working'})

population_analysis_to_perform={'sta'}; %,'pop'};

keys=struct;
project=varargin{1};
keys=ph_general_settings(project,keys);
project_specific_settings=[keys.db_folder filesep project filesep 'ph_project_settings.m'];
run(project_specific_settings);

if nargin>1
    keys.project_versions=varargin{2};
end
for f=1:numel(keys.project_versions) % running multiple versions of the same project at once !
    if ~isempty(keys.project_versions{f})
        keys.project_version=keys.project_versions{f};
    end
    version_folder=keys.project_version;
    keys.version_specific_settings=[keys.db_folder project filesep keys.project_version filesep 'ph_project_version_settings.m'];
    run(keys.version_specific_settings);
    keys.project_version=version_folder;
    keys.monkeys=keys.batching.monkeys;
    if keys.batching.combine_monkeys
        keys.batching.monkeys={''};
    end
    
    for m=1:numel(keys.batching.monkeys)
        keys.monkey=keys.batching.monkeys{m};
        keys.anova_table_file=[keys.basepath_to_save keys.project_version filesep 'tuning_table_combined_CI.mat'];
        population=ph_load_population([keys.basepath_to_save keys.project_version],['population_' keys.monkey]);
        population=ph_assign_perturbation_group(keys,population);
        population=ph_epochs(population,keys);
        
        for t=1:numel(keys.batching.targets)
            target=keys.batching.targets{t};
            
            for subregion=1:keys.batching.n_Subregions
                %population_selection={};
                if keys.batching.Subregions_separately
                    population_selection            = {'target',target;'Subregion', subregion};
                else
                    population_selection            = {'target',target};
                end
                keys.arrangement=keys.position_and_plotting_arrangements{1};
                loop_through_plots(population,keys,population_selection,population_analysis_to_perform)
                
            end
        end
    end
end
end

function loop_through_plots(population,keys_in,population_selection,population_analysis_to_perform)
for P=population_analysis_to_perform
    ana=P{:};
    AN=upper(ana(1:2));
    
    
    for cc=1:numel(keys_in.(ana))
        keys=keys_in;
        
        % presets
        keys.(AN).logarithmic_scale=0;
        keys.(AN).absolutes=0;
        keys.(AN).VMI='';
        keys.(AN).hist_column='';
        keys.(AN).color_option='monkeys_by_color';
        
        
        keys.(AN).epoch_BL='INI';
        keys.(AN).epoch_PF='Cue';
        %keys.(AN).epoch_NM='Cue'; %%%epoch_for_normalization ... divisive!
        keys.(AN).epoch_for_normalization='INI'; %%%epoch_for_normalization ... divisive!
        keys.(AN).epoch_RF='';
        
        
        keys.(AN).FR_subtract_baseline=0;
        keys.(AN).normalization='none'; %%%normalization
        keys.(AN).plot_RF=0; %%%normalization
        keys.(AN).epoch_RF=''; %%%normalization
        keys.(AN).y_lim=[]; %%%normalization
        keys.(AN).group_excluded={'','-'}; %%%normalization
        
        
        keys.(AN).IC_to_plot='in';
        keys.(AN).plot_per_position=0;
        
        %            keys.comparisons_per_effector           = keys.ons(cc).comparisons_per_effector;
        %            keys.comparisons_title                  = keys.ons(cc).comparisons_title;
        
        %conditions_to_plot                      = keys.ons(cc).conditions_to_plot;
        
        
        %% defaults
        
        keys.tt.choices                         = 0;
        keys.tt.perturbations                   = 0;
%         keys.tt.hands                           = [0 1 2];
        keys.tt.hands                           = 0;
        keys.tt.IC_for_criterion                = 'in';
        keys.tt.unselect                        = {};
        keys.tt.combine_tuning_properties       = {};
        
        % plot specific keys
        FN=fieldnames(keys_in.(ana)(cc));
        FN=FN(~ismember(FN,'tt'));
        for fn=FN'
            if ~isempty(keys_in.(ana)(cc).(fn{:}))
                keys.(AN).(fn{:})=keys_in.(ana)(cc).(fn{:});
            end
        end
        if isstruct(keys_in.(ana)(cc).tt)
            for fn=fieldnames(keys_in.(ana)(cc).tt)'
                keys.tt.(fn{:})=keys_in.(ana)(cc).tt.(fn{:});
            end
        end
      
        if isfield(keys_in.(ana)(cc),'tt') && isfield(keys_in.(ana)(cc).tt,'selection')
            keys.tt.selection                       = [population_selection ; keys_in.(ana)(cc).tt.selection];
        else
            keys.tt.selection = population_selection;
        end
        
        if isfield(keys.tt,'tasktypes')
            trial_criterion_independently_for_conditions_to_plot=0;
        else
            trial_criterion_independently_for_conditions_to_plot=1;
        end
        
        %% conditions to plot irrelevant for scatter, but still there should be something to be able to loop through
        if ~isfield(keys_in.(ana)(cc),'conditions_to_plot') || numel(keys_in.(ana)(cc).conditions_to_plot)==0;
            keys_in.(ana)(cc).conditions_to_plot={'ololol'};
        end
        
        for condition_to_plot=keys_in.(ana)(cc).conditions_to_plot
            keys.conditions_to_plot=condition_to_plot;
            keys.(AN).tasktypes=condition_to_plot;
            if trial_criterion_independently_for_conditions_to_plot
                keys.tt.tasktypes=strcat(condition_to_plot, ['_' keys_in.arrangement(1:3)]);
            end
            
            seed_filename=[keys.basepath_to_save keys.project_version filesep 'seed.mat'];
            if exist(seed_filename,'file');
                load(seed_filename);
                rng(seed);
            else
                seed=rng;
                save(seed_filename,'seed');
            end
            switch ana
                case 'ons'
                    ph_population_response_timing(population,keys);
                case 'pop'
                    ph_population_PSTHs(population,keys);
                case 'sct'
                    %     keys.monkey             = '';
                    %     keys.batching.combine_monkeys=1;
                    ph_scatter(keys);
                case 'ccs'
                    keys.(AN).epochs=keys.(AN).epochs.(condition_to_plot{:});
                    %     keys.monkey             = '';
                    %     keys.batching.combine_monkeys=1;
                    for pie=[0,1]
                        for percent=[0,1]
                            keys.CC.plot_as_pie=pie;
                            keys.CC.percent=percent;
                            ph_anova_cell_count(keys);
                        end
                    end
                case 'gaz'
                    ph_gaze(population,keys);
                case'sta'
                    ph_state_space(population,keys);
            end
        end
    end
    
end
end



