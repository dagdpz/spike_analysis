function ph_initiate_population_analysis(varargin)
%ph_initiate_population_analysis('PPC_pulv_eye_hand',{'LIP_dPul_inj_working'},{'pop','ons','sct','ccs'})
population_analysis_to_perform={'pop','ons','sct','ccs','gaz','ref','gfl','hst','beh'}; 

if nargin>2
    population_analysis_to_perform=varargin{3};
end
keys=struct;
project=varargin{1};
keys=ph_general_settings(project,keys);
project_specific_settings=[keys.db_folder 'ph_project_settings.m'];
run(project_specific_settings);

if nargin>1
    keys.project_versions=varargin{2};
end
for f=1:numel(keys.project_versions) % running multiple versions of the same project at once !
    if ~isempty(keys.project_versions{f})
        keys.project_version=keys.project_versions{f};
    end
    version_folder=keys.project_version;
    keys.version_specific_settings=[keys.db_folder keys.project_version filesep 'ph_project_version_settings.m'];
    run(keys.version_specific_settings);
    keys.project_version=version_folder;
    keys.monkeys=keys.batching.monkeys;
    keys.project=project;
    population_analysis_to_perform=population_analysis_to_perform(ismember(population_analysis_to_perform,fieldnames(keys)));
    if keys.batching.combine_monkeys
        keys.batching.monkeys={''};
    end    
    keys=ph_make_subfolder('version_log',keys);
    keys=ph_make_subfolder('population_meta_data',keys);
    save([keys.basepath_to_save keys.project_version filesep 'version_log' filesep datestr(clock,'YYYYmmDD-HHMM')],'keys');
    for m=1:numel(keys.batching.monkeys)
        keys.monkey=keys.batching.monkeys{m};
        keys.anova_table_file=[keys.basepath_to_save keys.project_version filesep 'tuning_table_combined_CI.mat'];
        if any(ismember(population_analysis_to_perform,{'ons','pop','gaz','ref','gfl','clf','hst','reg','rtc','prf','ndt'}))
            population=ph_load_population([keys.basepath_to_save keys.project_version],['population_' keys.monkey]);
            population=ph_assign_perturbation_group(keys,population);
            population=ph_epochs(population,keys);
        else
            population=[];
        end
        for t=1:numel(keys.batching.targets)
            target=keys.batching.targets{t};            
            for subregion=1:keys.batching.n_Subregions
                if keys.batching.Subregions_separately
                    population_selection            = {'target',target;'Subregion', subregion};
                else
                    population_selection            = {'target',target};
                end
                    loop_through_plots(population,keys,population_selection,population_analysis_to_perform);            
            end
        end
    end
end
end

function loop_through_plots(population_full,keys_in,population_selection,population_analysis_to_perform)
for P=population_analysis_to_perform
    ana=P{:};
    AN=upper(ana(1:2));
    for cc=1:numel(keys_in.(ana))
        keys=keys_in;        
        keys.(AN).position_and_plotting_arrangements=keys_in.position_and_plotting_arrangements;
        
        %% DEFAULTS
        
        %% PSTH keys and binsize

        keys.(AN).PSTH_binwidth=keys.PSTH_binwidth;
        keys.(AN).gaussian_kernel=keys.gaussian_kernel;
        
        %grouping keys
        keys.(AN).group_parameter='ungrouped';
        keys.(AN).group_excluded={'','-'}; 
        
        % presets
        keys.(AN).logarithmic_scale=0;
        keys.(AN).absolutes=0;                  %%%??????????????????????
        keys.(AN).VMI='';
        keys.(AN).hist_column='';
        keys.(AN).color_option='monkeys_by_color';      
        keys.(AN).epoch_BL='none';
        keys.(AN).epoch_PF='none';
        keys.(AN).epoch_DN='none';
        keys.(AN).epoch_RF='none';
        keys.(AN).fittypes={'sigmoidal','linear','gaussian1'}; %,'gaussian2' %,'linear'
        keys.(AN).baseline_per_trial=0;
        keys.(AN).normalization='none'; 
        keys.(AN).unpref_def='horizontally_opposite'; 
        keys.(AN).permutation_tests=1; 
        keys.(AN).plot_as_pie=1; %%%%%%%%%%
        keys.(AN).percent=0;
        
        % plotting keys
        %keys.(AN).plot_RF=0; 
        keys.(AN).plot_per_position=0;
        keys.(AN).plot_posneg=0;
        keys.(AN).y_lim=[]; 
        keys.(AN).link_y_lim=1; 
        keys.(AN).IC_to_plot='in';
        keys.(AN).fontsize_factor=1.5;
        % keys.(AN).fontsize_factor=4;  %% MP: this should be in the settings file
        keys.(AN).split_SUs={''};
        keys.(AN).RF_frame_colors                 	= {};
        keys.(AN).RF_frame_entries                 	= {};
        keys.(AN).RF_frame_parameter                = '';
        keys.(AN).RF_columns                        = [];
        keys.(AN).RF_rows                           = [];
        
        %% key asignment
        FN=fieldnames(keys_in.(ana)(cc));
        FN=FN(~ismember(FN,{'tt','cal'}));
        for fn=FN'
            if ~isempty(keys_in.(ana)(cc).(fn{:}))
                keys.(AN).(fn{:})=keys_in.(ana)(cc).(fn{:});
                if strcmp(fn{:},'condition_parameters')
                keys.(fn{:})=keys_in.(ana)(cc).(fn{:});
                end
            end
        end
        
        keys.(AN).FR_subtract_baseline=~strcmp(keys.(AN).epoch_BL,'none');
        if isfield(keys_in.(ana)(cc),'tt') && isstruct(keys_in.(ana)(cc).tt)
            for fn=fieldnames(keys_in.(ana)(cc).tt)'
                keys.tt.(fn{:})=keys_in.(ana)(cc).tt.(fn{:});
            end
        end
        % new part for condition limiting !1
        if isfield(keys_in.(ana)(cc),'cal') && isstruct(keys_in.(ana)(cc).cal)
            for fn=fieldnames(keys_in.(ana)(cc).cal)'
                keys.cal.(fn{:})=keys_in.(ana)(cc).cal.(fn{:});
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
        
        population = ph_accept_trials_per_unit(population_full,keys);
        %population=population_full;
        
        for a=1:numel(keys.(AN).position_and_plotting_arrangements)            
            keys.arrangement= keys.(AN).position_and_plotting_arrangements{a};  
            for condition_to_plot=keys_in.(ana)(cc).conditions_to_plot
                keys.conditions_to_plot=condition_to_plot;
                keys.(AN).tasktypes=condition_to_plot;
                if trial_criterion_independently_for_conditions_to_plot
                    keys.tt.tasktypes=strcat(condition_to_plot, ['_' keys.arrangement(1:3)]);
                    %keys.tt.tasktypes=strcat(condition_to_plot{:}, ['_' keys.arrangement(1:3)]);
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
                    case 'rtc'                        
                        keys=ph_make_subfolder('reaction_time_correlation',keys);
                        ph_RT_correlation(population,keys);
                    case 'reg'                        
                        keys=ph_make_subfolder(['population_meta_data' filesep 'regression'],keys);
                        keys=ph_make_subfolder('regression',keys);
                        ph_linear_regression(population,keys);
                    case 'ons'                        
                        keys=ph_make_subfolder(['population_meta_data' filesep 'response timing'],keys);
                        keys=ph_make_subfolder('response timing',keys);
                        ph_population_response_timing(population,keys);
                    case 'beh'
                        keys=ph_make_subfolder(['population_meta_data' filesep 'behavior'],keys);
                        keys=ph_make_subfolder('behavior',keys);
                        ph_ephys_behavior(keys);
                    case 'pop'
                        keys=ph_make_subfolder(['population_meta_data' filesep 'population_analysis'],keys);
                        keys=ph_make_subfolder('population_analysis',keys);
                        ph_population_PSTHs(population,keys);
                    case 'sct'
                        keys=ph_make_subfolder(['population_meta_data' filesep 'scatter'],keys);
                        keys=ph_make_subfolder('scatter',keys);
                        ph_scatter(keys);
                    case 'ccs'
                        temp_epochs=keys.(AN).epochs;
                        keys.(AN).epochs=keys.(AN).epochs.(condition_to_plot{1});
                        %% missing storing stuff
                        ph_anova_cell_count(keys);
                        keys.(AN).epochs=temp_epochs;
                    case 'gaz'
                        keys=ph_make_subfolder(['population_meta_data' filesep 'gaze_analysis'],keys);
                        keys=ph_make_subfolder('gaze_analysis',keys);
                        ph_gaze(population,keys);
                    case 'ref'
                        keys=ph_make_subfolder(['population_meta_data' filesep 'scatter'],keys);
                        keys=ph_make_subfolder('scatter',keys);
                        ph_reference_frame(population,keys); 
                    case 'prf'
                        keys=ph_make_subfolder(['population_meta_data' filesep 'target_preference'],keys);
                        keys=ph_make_subfolder('target_preference',keys);
                        ph_target_preference(population,keys); 
                    case 'gfl'
                        keys=ph_make_subfolder(['population_meta_data' filesep 'tuning_curves'],keys);
                        keys=ph_make_subfolder('tuning_curves',keys);
                        ph_gainfields(population,keys); 
                    case 'clf'
                        keys=ph_make_subfolder(['population_meta_data' filesep 'classification'],keys);
                        keys=ph_make_subfolder('classification',keys);
                        ph_classification(population,keys); 
                    case 'hst'
                        keys=ph_make_subfolder(['population_meta_data' filesep 'hand_space'],keys);
                        keys=ph_make_subfolder('hand_space',keys);
                        temp_epochs=keys.(AN).epochs;
                        keys.(AN).epochs=keys.(AN).epochs.(condition_to_plot{:});
                        ph_hand_space_tuning_vector(population,keys); 
                        keys.(AN).epochs=temp_epochs;
                    case 'cpy'
                        keys=ph_make_subfolder(keys.CP.foldername,keys);
                        ph_copy_single_units(keys); 
                    case 'ndt'
                        keys=ph_make_subfolder(['population_meta_data' filesep 'decoding'],keys);
                        keys=ph_make_subfolder('decoding',keys);
                        ph_decoding(population,keys); 
                end
            end
        end
    end    
end
end

function keys=ph_make_subfolder(foldername,keys)
keys.path_to_save=[keys.basepath_to_save keys.project_version filesep foldername filesep];
if ~exist(keys.path_to_save,'dir')
    mkdir([keys.basepath_to_save keys.project_version ], foldername);
end
end