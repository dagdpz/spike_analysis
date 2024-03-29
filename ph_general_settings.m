function keys=ph_general_settings(project,keys)



%% general settings (multi-summary PSTH)
keys.PSTH_binwidth                      =0.01;                      % resolution of PSTH's (in seconds)
keys.gaussian_kernel                    =0.02;                      % std for the convolution to derive spie density (in seconds)
keys.kernel_type                        ='gaussian';                % could also be 'box'
keys.FR_at_peak                         =0;                         % currently not used
keys.position_and_plotting_arrangements ={'hands'};                 % defines position batching and which conditions go into different figures/lines
keys.condition_parameters={'choice','reach_hand','perturbation'};
keys.contra_ipsi_relative_to='target';

%% labels for conditions to plot - related to colors fieldnames!
keys.labels.handsIC                 ={'AH','IH','CH'};  %% AH!??
keys.labels.perturbation            ={'','PT','PT2','PT3','PT4','PT5','PT6','PT7','PT8'};
keys.labels.reach_hand              ={'AH','IH','CH'};
keys.labels.reach_handLR            ={'AH','LH','RH'}; 
keys.labels.choice                  ={'in','ch'};
keys.labels.stimulustype            ={'SS','TT','TD'};
keys.labels.difficulty              ={'TA','D1','D2'};
keys.labels.success                 ={'ER','SU'};
keys.labels.hemifield               ={'IS','VS','CS'};
keys.labels.fix_index               ={'IF','MF','CF'};
keys.labels.preferred               ={'NP','PF'};
keys.labels.stimuli_in_2hemifields  ={'1H','2H'};

[keys.colors,keys.linestyles]=ph_default_color_settings;

%% labels for tuning table significance entries
keys.TTconditions.hands             ={'AH','LH','RH'};
keys.TTconditions.choice            ={'in','ch'};
keys.TTconditions.hemifield         ={'LS','RS'};
keys.TTconditions.perturbation      ={'','_PT'};
keys.TTlabels.UD                    ={'DN','-','UP'};
keys.TTlabels.CR                    ={'UC','-','CR'};
keys.TTlabels.choices               ={'in','-','ch'}; %% why choice(S!)
keys.TTlabels.true                  ={'YE','NO','YE'}; %% 
keys.TTlabels.SglTar_Suc            ={'LS','-','RS'};
keys.TTlabels.Difficulty_Easy       ={'Ta','-','eD'}; %higher FR for T , higher FR for D
keys.TTlabels.Difficulty_Diff       ={'Ta','-','dD'}; %higher FR for T , higher FR for D
keys.TTlabels.SpatialComp_1HFTar    ={'ST','-','1T'}; %higher FR for T , higher FR for D
keys.TTlabels.SpatialComp_2HFTar    ={'ST','-','2T'}; %higher FR for single T , higher FR for double target
keys.TTlabels.epoch                 ={'su','-','en','bi'};
keys.TTlabels.hemifield             ={'LS','-','RS'}; 
keys.TTlabels.hands                 ={'LH','-','RH'};
keys.TTlabels.perturbation          ={'SU','-','EN'};

%% Batching per figure ! subregion keys..?
keys.batching.monkeys                    ={'Curius','Linus'};       % monkeys on the project
keys.batching.combine_monkeys            =0;                        % for population analysis, monkeys can be separately or combined
keys.batching.targets                    ={'dPulv_r','dPulv_l'};    % to combine both, just put {'dPulv'}
keys.batching.Subregions_separately      =0;                        % subregions can be processed independently (if only defining target is too crude
keys.batching.Subregions{1}{1}           =struct('monkey',{''},'target',{''},'grid_x',{NaN},'grid_y',{NaN},'z_min',{NaN},'z_max',{NaN}); % subregion definitions by  z range for each grid hole
keys.batching.n_Subregions               =numel(keys.batching.Subregions);

%% criterions to exclude trials and units
keys.cal.process_spikes                  =1;      % you can choose not to run spikes at all
keys.cal.process_sites                   =1;      % you can choose not to run lfp sites at all (saving processing time)
keys.cal.process_by_block                =1;      % you can choose not to run by block (body signals f.e.) at all (saving processing time)
keys.cal.MA_selection                   ={'display',0,'keep_raw_data',1,'saccade_definition',4,'reach_1st_pos',1,'correlation_conditions',{}};                        % if you want to run MA with specific settings
keys.cal.units_from_sorting_table       =1;                         % exclude units that are not in the sorting table (and therefore apply stability/single/SNT ratings)
keys.cal.datasets                       =[];
keys.cal.completed                      =1;                         % problematic, because of where and how it is used. so far, keep it 1
keys.cal.effectors                      =[0,1,2,3,4,5,6];           % excluding trials with non-matching effectors
keys.cal.types                          =[1,2,3,4,5,6];             % excluding trials with non-matching types
keys.cal.reach_hand                     =[0,1,2];                   % excluding trials with non-matching reach_hand
keys.cal.perturbation                   =[0,1];                   % excluding trials with non-matching reach_hand

keys.cal.choice                         =[0,1];                     % excluding trials with non-matching chocie
keys.cal.stablity                       =[0,1];                     % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the table
keys.cal.single_rating                  =[1,2,3];                   % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the tabl

keys.cal.FR                             =[0,inf];       % min and max value accepted    
keys.cal.n_spikes                       =[0,inf];       % min and max value accepted

keys.cal.automatic_stablity             =0;                         % using automatic stability assessment
keys.cal.automatic_SNR                  =0;                         % using automatic SNR assessment
keys.cal.SNR_rating                     =[1,2,3,4];                 % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the tabl
%keys.cal.min_trials_per_condition       =5;                         % minimum trials per conditon (look at ph_arrange_positions to see how conditions are defined)
keys.cal.min_spikes_per_unit            =50;                        % excluding units that have in total less spikes (workaround for sortcode assignment bug) - to be removed
keys.cal.perturbation_groups            ={0,[2,3,4,5,6,7,8]};       % which perturbation values from excel table will be assigned to control and perturbation for comparisons and population analysis
keys.cal.remove_trials_without_spikes=1;
keys.cal.remove_trials_with_outlying_FR=1;


keys.cal.min_trials_in                  =5;                   % minimum number of trials instructed
keys.cal.min_trials_ch                  =5;                   % minimum number of trials choice

%% precision for merging positions (and fixations)
keys.cal.precision_fix=4;
keys.cal.precision_tar=2;

%% ANOVA settings
keys.AN.test_types='parametric'; %% as opposed to 'nonparametric' and 'permutations'
keys.AN.n_permutations                 = 10000; % in case of permutations test_type
keys.AN.n_permutations_randanova2      = 1000;
keys.AN.n_permutations_randanova1      = 1000;
keys.AN.check_normality=0; %% very time consuming

%%normalization for ANOVAs
keys.AN.normalization='none';
keys.AN.epoch_DN='none';
keys.AN.epoch_RF='none';
keys.AN.epoch_BL='none';
keys.AN.epoch_PF='none';
keys.AN.baseline_per_trial=0;

% multicomparison - correct p-values on the fly to update significance table entries
keys.AN.multicomparison='none';
keys.AN.factors={'epoch','hemifield','hands','SxH','SglTar_Suc'};
keys.AN.factors_per_hand={'epoch','position','fixation','PxF','distance','angle','DxA','prefH', 'prefP','positionx','positiony','_positionxy','gaze_modulation_x','gaze_modulation_y','gaze_pref_x','gaze_pref_y'};
keys.AN.factors_per_hemifield={'Difficulty_Easy', 'Difficulty_Diff', 'SpatialComp_2HFTar', 'SpatialComp_1HFTar'};
keys.AN.factors_per_handhemifield={'epoch','prefH', 'prefP', 'SpatialComp_2HFTar', 'SpatialComp_1HFTar'}; %% basically perturbation and effector comparison


%% folders and filenames
keys.filelist_as_blocks     =0;
keys.drive                  =DAG_get_server_IP;
keys.basepath_to_save       =[keys.drive 'Projects' filesep project filesep 'ephys' filesep];
spike_analysis_location     =which('ph_initiation');
keys.db_folder              =[spike_analysis_location(1:strfind(spike_analysis_location,['spike_analysis' filesep 'ph_initiation'])-1) 'Settings' filesep project filesep 'spike_analysis' filesep];
%% this folder defines where to take settings from

keys.All_monkeys={'Flaffus','Linus','Curius','Tesla','Cornelius','Magnus','TDT_brain','Bacchus'};
for m=1:numel(keys.All_monkeys)
    keys.(keys.All_monkeys{m}).sorted_neurons_filename    =[keys.All_monkeys{m}(1:3) '_sorted_neurons'];
    keys.(keys.All_monkeys{m}).sorted_neurons_sheetname    ='final_sorting';
    keys.(keys.All_monkeys{m}).filelist_formatted         ={};
end

keys.Flaffus.color    ='r';
keys.Linus.color      =[0 0 255]/255;
keys.Curius.color     =[255 0 0]/255;
keys.Tesla.color      ='y';
keys.Cornelius.color  ='m';
keys.Bacchus.color    ='b';

keys.Flaffus.marker    ='<';
keys.Linus.marker      ='s';
keys.Curius.marker     ='o';
keys.Tesla.marker      ='d';
keys.Cornelius.marker  ='v';
keys.Bacchus.marker    ='x';

%% plotting settings
keys.plot.anova_tables                  ='off';     % display anova result tables for each unit, please keep off, it will get very messy if you turn this on
keys.plot.single_trials                 =0;         % not used yet
%keys.plot.single_cells                  =1;         % perform single cell plotting
keys.plot.waveforms                     =1;         % plot the waveform summary plots per session
keys.plot.polars_on_extra_figure        =0;         % recommended if there are too many conditions and the single cell heatmap plots are too small/crowded
keys.plot.eye_hand_traces               =1;         % Incredible performance booster if turned off
keys.plot.average_PSTH_line             =0;         % One PSTH line on top that represents the average of all others
%keys.plot.average_heat_maps             =0;
keys.plot.export                        =1;         % save plots as pdfs, you typically want this
keys.plot.events                        =1:100;     % select events that should be plotted (vertical lines) on all PSTH like plots
keys.plot.population_PSTH_legends       =1;         % if population legends should be plotted or not
keys.plot.cell_count_legends            =1;         % if population legends should be plotted or not
keys.plot.scatter_legends               =1;         % if population legends should be plotted or not
keys.plot.rotate_epoch_labels           =0;
keys.plot.rotate_time_labels            =0;



%% tuning table readout options (excluding particular subsets)
keys.tt.combine_tuning_properties   ={'place_name_here'}; %% additional table entry from combining columns
keys.tt.replace_tuning              ={}; %% additional table entry from combining columns
keys.tt.perturbation                = 0;
keys.tt.choice                      = 0;
keys.tt.reach_hand                  = 0;
keys.tt.IC_for_criterion            = 'in';
keys.tt.trial_criterion_in          ='per_position';
keys.tt.trial_criterion_ch          ='per_hemifield';
keys.tt.selection                   ={};
keys.tt.unselect                    ={};
keys.tt.selected_list               ={};
keys.tt.unselected_list             ={};
keys.tt.type_effectors              ={'Msac'};

%% population
keys.sct=struct([]);
keys.ccs=struct([]);
keys.ons=struct([]);
keys.pop=struct([]);
keys.hst=struct([]);

%% EPOCH SETTINGS 
% For each task type seperately, analysis epochs are defined

%% Fixation only type
keys.EPOCHS_PER_TYPE{1}={...
    'INI',      2,  -0.4,	-0.1,   'INI';...
    'Facq',     2,	0.05,	0.15,   'INI';...
    'Fhol',     3,	0.2,	0.7,      'INI';...
    };
keys.WINDOWS_PER_TYPE{1}={...
    'Aquisition',   3,	-0.5,	0.4;...
    'Fixation',     3,	0.5,   1;...
    };
keys.ANOVAS_PER_TYPE(1).epoch={'INI' 'Facq';...
    'Facq' 'Fhol'};

keys.ANOVAS_PER_TYPE(1).hemifield          ={'Facq','Fhol'}';
keys.ANOVAS_PER_TYPE(1).positions          ={'Facq','Fhol'}';
keys.ANOVAS_PER_TYPE(1).hands              ={'Facq','Fhol'}';
keys.ANOVAS_PER_TYPE(1).SxH                ={'Facq','Fhol'}';
keys.ANOVAS_PER_TYPE(1).main               ={'Facq','Fhol'}';

%% Visually guided type
keys.EPOCHS_PER_TYPE{2}={...
    'INI',      2,	-0.4,	-0.1,   'INI';...
    'Facq',     3,	0.05,	0.15,   'INI';...
    'Fhol',     4,	-0.3,	0,      'INI';...
    'Cue',      4,	0.05,	0.15,   'INI';...
    'PreS',     60,	-0.1,	-0.01,  'INI';...
    'PeriS',    60,	-0.01,	0.05,   'INI';...
    'Tacq',     5,	0,      0.15,   'INI';...
    'Thol',     20,	-0.3,	0,      'INI';...
    };
keys.ANOVAS_PER_TYPE(2).epoch={'INI' 'Facq';...
    'INI' 'Fhol';...
    'Fhol' 'Cue';...
    'Fhol' 'PreS';...
    'Fhol' 'PeriS';...
    'Fhol' 'Tacq';...
    'Fhol' 'Thol'};
keys.WINDOWS_PER_TYPE{2}={...
    'Initiation',   2,	-0.5,	0;...
    'Visual onset', 4,	-0.8,   0.17;...
    'Saccade',      60,	-0.01,  0.22;...
    'T hold',       20,	-0.3,   0.1;...
    };
keys.ANOVAS_PER_TYPE(2).hemifield            ={'Cue','PreS','PeriS','Tacq','Thol'}';
keys.ANOVAS_PER_TYPE(2).positions          ={'Cue','PreS','PeriS','Tacq','Thol'}';
keys.ANOVAS_PER_TYPE(2).hands              ={'Facq','Fhol','Cue','PreS','PeriS','Tacq','Thol'}';
keys.ANOVAS_PER_TYPE(2).SxH                ={'Cue','PreS','PeriS','Tacq','Thol'}';
keys.ANOVAS_PER_TYPE(2).main               ={'Facq','Fhol','Cue','PreS','PeriS','Tacq','Thol'}';

%% Memory type
keys.EPOCHS_PER_TYPE{3}={...
    'INI',      2,	-0.4,	-0.1,   'INI';...
    'Facq',     3,	0.05,	0.15,   'INI';...
    'Fhol',     6,	-0.3,	0,      'INI';...
    'Cue',      6,	0.05,   0.15,   'INI';...
    'MemE',     7, 	0,      0.2,    'INI';...
    'MemL',     9,	-0.3, 	0,      'INI';...
    'PreS',     60,	-0.1, 	-0.01,  'INI';...
    'PeriS',	60,	-0.01, 	0.05,   'INI';...
    'TIhol',	10,	0,      0.1,    'INI';...
    'Thol',     20,	-0.3,   0,      'INI';...
    };
keys.WINDOWS_PER_TYPE{3}={...
    'Initiation',   2,	-0.5,	0;...
    'Cue',          6,	-0.8,   0.78;...
    'Saccade',      60,	-0.7,   0.05;...
    'T hold',       20,	-0.6,   0;...
    };
keys.ANOVAS_PER_TYPE(3).epoch={'INI' 'Facq';...
    'INI' 'Fhol';...
    'Fhol' 'Cue';...
    'Fhol' 'MemE';...
    'Fhol' 'MemL';...
    'Fhol' 'PreS';...
    'Fhol' 'PeriS';...
    'Fhol' 'TIhol';...
    'Fhol' 'Thol';...
    };
keys.ANOVAS_PER_TYPE(3).hemifield            ={'Cue','MemE','MemL','PreS','PeriS','TIhol','Thol'}';
keys.ANOVAS_PER_TYPE(3).positions             ={'Cue','MemE','MemL','PreS','PeriS','TIhol','Thol'}';
keys.ANOVAS_PER_TYPE(3).hands              ={'Facq','Fhol','Cue','MemE','MemL','PreS','PeriS','TIhol','Thol'}';
keys.ANOVAS_PER_TYPE(3).SxH                ={'Cue','MemE','MemL','PreS','PeriS','TIhol','Thol'}';
keys.ANOVAS_PER_TYPE(3).main              ={'Facq','Fhol','Cue','MemE','MemL','PreS','PeriS','TIhol','Thol'}';

%% Delay type
keys.EPOCHS_PER_TYPE{4}={...
    'INI',      2,	-0.4,	-0.1,   'INI';...
    'Facq',     3,	-0.4,	-0.1,   'INI';...
    'Fhol',     6,	-0.3,	0,      'INI';...
    'Cue',      6,	0.05,   0.15,   'INI';...
    'EDel',     4, 	-0.6,   -0.3,   'INI';...
    'Del',      4, 	-0.3,   0,      'INI';...
    'PreS',     60,	-0.1, 	-0.01,  'INI';...
    'PeriS',	60,	-0.01, 	0.05,   'INI';...
    'PostS',	61,	0.05,   0.2,    'INI';...
    'PreR',     62,	-0.3, 	-0.01,  'INI';...
    'PeriR',	62,	-0.05, 	0.15,   'INI';...
    'PostR',	63,	0.05,   0.2,    'INI';...
    'Thol',     20,	-0.3,   0,      'INI';...
    };
keys.WINDOWS_PER_TYPE{4}={...
    'Initiation',   2,	-0.5,	0.4;...
    'Fixation',     3,	-1.2,   0.17;...
    'Delay Period', 6,	-0.33,  1.35;...
    'Saccade',      60,	-0.1,   0.6;...
    'Reach',        62,	-0.35,  0.7;...
    };
keys.ANOVAS_PER_TYPE(4).epoch={'INI' 'Facq';...
    'INI' 'Fhol';...
    'Fhol' 'Cue';...
    'Fhol' 'EDel';...
    'EDel' 'Del';...
    'EDel' 'PreS';...
    'EDel' 'PeriS';...
    'EDel' 'PostS';...
    'EDel' 'PreR';...
    'EDel' 'PeriR';...
    'EDel' 'PostR';...
    'Fhol' 'Thol';...
    };
keys.ANOVAS_PER_TYPE(4).hemifield            ={'INI','Facq','Fhol','Cue','EDel','Del','PreS','PeriS','PostS','PreR','PeriR','PostR','Thol'}';
keys.ANOVAS_PER_TYPE(4).positions          ={'INI','Facq','Fhol','Cue','EDel','Del','PreS','PeriS','PostS','PreR','PeriR','PostR','Thol'}';
keys.ANOVAS_PER_TYPE(4).hands              ={'INI','Facq','Fhol','Cue','EDel','Del','PreS','PeriS','PostS','PreR','PeriR','PostR','Thol'}';
keys.ANOVAS_PER_TYPE(4).SxH                ={'INI','Facq','Fhol','Cue','EDel','Del','PreS','PeriS','PostS','PreR','PeriR','PostR','Thol'}';
keys.ANOVAS_PER_TYPE(4).main                ={'INI','Facq','Fhol','Cue','EDel','Del','PreS','PeriS','PostS','PreR','PeriR','PostR','Thol'}';

%% M2S
keys.EPOCHS_PER_TYPE{5}={...
    'INI',      2,      -0.4,	-0.1,   'INI';...
    'Facq',     3,      0.05,	0.15,   'INI';...
    'Fhol',     6,      -0.3,	0,      'INI';...
    'Cue',      6,      0.05,   0.2,   'INI';...
    'MemE',     7,      0,      0.2,    'INI';...
    'PreS',     60,	    -0.1, 	-0.01,  'INI';...
    'PeriS',	60,     -0.01, 	0.05,   'INI';...
    'Exp',      14,     -0.8, 	-0.2,   'INI';...
    'Thol',     90,     -0.95,  -0.8,   'INI';...
    };
keys.WINDOWS_PER_TYPE{5}={...
    'INI',          2,	-0.5,	0.3     ;...
    'Sample',       6,	-0.4,   2       ;...
    'Selection',    90,	-2,     0       ;...
    };
keys.ANOVAS_PER_TYPE(5).epoch={'INI' 'Facq';...
    'Facq' 'Fhol';...
    'Fhol' 'Cue';...
    'Fhol' 'MemE';...
    'Fhol' 'PreS';...
    'Fhol' 'PeriS';...
    'Fhol' 'Thol';...
    };
keys.ANOVAS_PER_TYPE(5).hemifield            ={'Cue','MemE','PreS','PeriS','Thol'}';
keys.ANOVAS_PER_TYPE(5).positions          ={'Cue','MemE','PreS','PeriS','Thol'}';
keys.ANOVAS_PER_TYPE(5).hands              ={'Facq','Fhol','Cue','MemE','PreS','PeriS','Thol'}';
keys.ANOVAS_PER_TYPE(5).SxH                ={'Cue','MemE','PreS','PeriS','Thol'}';
keys.ANOVAS_PER_TYPE(5).main              ={'Facq','Fhol','Cue','MemE','PreS','PeriS','Thol'}';

%% M2S
keys.EPOCHS_PER_TYPE{6}={...
    'INI',      2,      -0.4,	-0.1,   'INI';...
    'Facq',     3,      0.05,	0.15,   'INI';...
    'Fhol',     6,      -0.3,	0,      'INI';...
    'Cue',      6,      0.05,   0.2,   'INI';...
    'MemE',     7,      0,      0.2,    'INI';...
    'PreS',     60,	    -0.1, 	-0.01,  'INI';...
    'PeriS',	60,     -0.01, 	0.05,   'INI';...
    'Exp',      14,     -0.8, 	-0.2,   'INI';...
    'Thol',     90,     -0.95,  -0.8,   'INI';...
    };
keys.WINDOWS_PER_TYPE{6}={...
    'INI',          2,	-0.5,	0.3     ;...
    'Sample',       6,	-0.4,   2       ;...
    'Selection',    90,	-2,     0       ;...
    };
keys.ANOVAS_PER_TYPE(6).epoch={'INI' 'Facq';...
    'Facq' 'Fhol';...
    'Fhol' 'Cue';...
    'Fhol' 'MemE';...
    'Fhol' 'PeriS';...
    'Fhol' 'Exp';...
    'Fhol' 'Thol';...
    };
keys.ANOVAS_PER_TYPE(6).hemifield            ={'Cue','MemE','PeriS','Exp','Thol'}';
keys.ANOVAS_PER_TYPE(6).positions          ={'Cue','MemE','PeriS','Exp','Thol'}';
keys.ANOVAS_PER_TYPE(6).hands              ={'Facq','Fhol','Cue','MemE','PeriS','Exp','Thol'}';
keys.ANOVAS_PER_TYPE(6).SxH                ={'Cue','MemE','PreS','Exp','Thol'}';
keys.ANOVAS_PER_TYPE(6).main              ={'Facq','Fhol','Cue','MemE','PeriS','Exp','Thol'}';

%% multicomparison
keys.AN.multicomp_epochs=keys.ANOVAS_PER_TYPE;
keys.AN.multicomparison='none';

%% general population settings
  % multicomparison
  keys.population_defaults.separate_multicomparison=0;
  
  %grouping keys
  keys.population_defaults.group_parameter='ungrouped';
  keys.population_defaults.group_excluded={'','-'};
  
  % presets
  keys.population_defaults.logarithmic_scale=0;
  keys.population_defaults.absolutes=0;                  %%%??????????????????????
  keys.population_defaults.permutation_tests=1; % this refers to a specific population analysis, no?
  keys.population_defaults.VMI='';
  keys.population_defaults.hist_column='';
  keys.population_defaults.color_option='monkeys_by_color';
  keys.population_defaults.epoch_BL='none';
  keys.population_defaults.epoch_PF='none';
  keys.population_defaults.epoch_DN='none';
  keys.population_defaults.epoch_RF='none';
  keys.population_defaults.fittypes={'sigmoidal','linear','gaussian1'}; %,'gaussian2' %,'linear'
  keys.population_defaults.baseline_per_trial=0;
  keys.population_defaults.normalization='none';
  keys.population_defaults.unpref_def='horizontally_opposite';
  keys.population_defaults.plot_as_pie=1; % cell counts
  
  % plotting keys
  keys.population_defaults.unique_title='';
  keys.population_defaults.plot_per_position=0;
  keys.population_defaults.plot_posneg=0;
  keys.population_defaults.y_lim=[];
  keys.population_defaults.link_y_lim=1;
  keys.population_defaults.IC_to_plot='in';
  keys.population_defaults.fontsize_factor=1.5;
  keys.population_defaults.link_y_lim                        = 1;
  % keys.population_defaults.fontsize_factor=4;  %% MP: this should be in the settings file
  %keys.population_defaults.split_SUs={''};            %% ??
  keys.population_defaults.RF_frame_colors                 	= {};
  keys.population_defaults.RF_frame_entries                 	= {};
  keys.population_defaults.RF_frame_parameter                = '';
  keys.population_defaults.RF_columns                        = [];
  keys.population_defaults.RF_rows                           = [];
  

  keys.population_defaults.redo_statistics                   = 0; %% maybe even put 1 here


%% Single cell plot settings
% ylimits
%keys.plot.FR_max_for_ylim               =50;
keys.population_defaults.trials_max_for_ylim           =50;
keys.population_defaults.excentricity_max_for_ylim     =30;
keys.population_defaults.eyetrace_factor               =0.5;
keys.population_defaults.hndtrace_factor               =0.5;
keys.population_defaults.line_labelling                ='left/right';

% ANOVA labels for single unit plots
keys.population_defaults.anova_main    ={'E','in_epoch_main','','S','in_hemifield_main','','C','ch_hemifield_main','','H','in_hands_main','','ExS','in_ExS','','ExH','in_ExH','','SxH','in_SxH',''};
keys.population_defaults.anova_effector={'E','in_epoch_main','','S','in_hemifield_main','','C','ch_hemifield_main','','H','in_hands_main','','ExS','in_ExS','','ExH','in_ExH','','SxH','in_SxH',''};
keys.population_defaults.anova_epoch1  ={'E','in_AH','epoch','S','in','hemifield','C','ch','hemifield','H','in','hands','SxH','in','SxH'};
keys.population_defaults.anova_epoch2  ={'LL','in_LH_LS','PT','RL','in_LH_RS','PT','LR','in_RH_LS','PT','RR','in_RH_RS','PT'};


keys.population_defaults.PSTH_perpos_width          =0.5;
keys.population_defaults.raster_width               =0.1;
keys.population_defaults.PSTH_summary_width         =1;
  
  

