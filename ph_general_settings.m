function keys=ph_general_settings(project,keys)
%% condition_identifiers
keys.labels.handsLR={'AH','LH','RH'};
keys.labels.handsIC={'AH','IH','CH'};  %% AH!??
keys.labels.perturbations={'','_PT'};  %% AH!??
keys.labels.choices={'in','ch'};

keys.FR_subtract_baseline=0;

%% general settings (multi-summary PSTH)

keys.PSTH_binwidth                      =0.01;                      % resolution of PSTH's (in seconds)
keys.gaussian_kernel                    =0.02;                      % std for the convolution to derive spie density (in seconds)
keys.subtract_baseline_for_anovas       =0;                         % subtracting the FR in baseline period (see below) for anovas
keys.FR_at_peak                         =0;                         % currently not used
keys.position_and_plotting_arrangements ={'hands'};                 % defines position batching and which conditions go into different figures/lines

%% Batching per figure ! subregion keys..?
keys.batching.monkeys                    ={'Curius','Linus'};       % monkeys on the project
keys.batching.combine_monkeys            =0;                        % for population analysis, monkeys can be separately or combined
keys.batching.targets                    ={'dPulv_r','dPulv_l'};    % to combine both, just put {'dPulv'}
keys.batching.Subregions_separately      =0;                        % subregions can be processed independently (if only defining target is too crude
keys.batching.Subregions{1}{1}           =struct('monkey',{''},'target',{''},'grid_x',{NaN},'grid_y',{NaN},'z_min',{NaN},'z_max',{NaN}); % subregion definitions by  z range for each grid hole
keys.batching.n_Subregions               =numel(keys.batching.Subregions);

%% criterions to exclude trials and units
keys.cal.MA_selection                   ={'display',0,'keep_raw_data',1,'saccade_definition',4,'reach_1st_pos',1,'correlation_conditions',{}};                        % if you want to run MA with specific settings
keys.cal.units_from_sorting_table       =1;                         % exclude units that are not in the sorting table (and therefore apply stability/single/SNT ratings)
keys.cal.datasets                       =[];
keys.cal.completed                      =1;                         % problematic, because of where and how it is used. so far, keep it 1
keys.cal.effectors                      =[0,1,2,3,4,5,6];           % excluding trials with non-matching effectors
keys.cal.types                          =[1,2,3,4,5,6];             % excluding trials with non-matching types
keys.cal.reach_hand                     =[0,1,2];                   % excluding trials with non-matching reach_hand
keys.cal.choice                         =[0,1];                     % excluding trials with non-matching chocie
keys.cal.stablity                       =[0,1];                     % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the table
keys.cal.single_rating                  =[1,2,3];                   % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the tabl
keys.cal.SNR_rating                     =[1,2,3,4];                 % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the tabl
keys.cal.min_trials_per_condition       =5;                         % minimum trials per conditon (look at ph_arrange_positions to see how conditions are defined)
keys.cal.min_spikes_per_unit            =50;                        % excluding units that have in total less spikes (workaround for sortcode assignment bug) - to be removed
keys.cal.perturbation_groups            ={0,[2,3,4,5,6,7,8]};       % which perturbation values from excel table will be assigned to control and perturbation for comparisons and population analysis

%% folders and filenames
keys.filelist_as_blocks     =0;
keys.drive                  =DAG_get_server_IP;
keys.basepath_to_save       =[keys.drive 'Projects' filesep project filesep 'ephys' filesep];
keys.project_versions       ={''};
spike_analysis_location     =which('ph_initiation');
keys.db_folder              =[spike_analysis_location(1:strfind(spike_analysis_location,['spike_analysis' filesep 'ph_initiation'])-1) 'Settings' filesep 'spike_analysis' filesep];

keys.Flaffus.sorted_neurons_filename    ='Fla_sorted_neurons';
keys.Linus.sorted_neurons_filename      ='Lin_sorted_neurons';
keys.Curius.sorted_neurons_filename     ='Cur_sorted_neurons';
keys.Tesla.sorted_neurons_filename      ='Tes_sorted_neurons';
keys.Cornelius.sorted_neurons_filename  ='Cor_sorted_neurons';
keys.Magnus.sorted_neurons_filename     ='Mag_sorted_neurons';
keys.TDT_brain.sorted_neurons_filename  ='TDT_sorted_neurons';
keys.Bacchus.sorted_neurons_filename    ='Bac_sorted_neurons';

keys.Flaffus.sorted_neurons_sheetname    ='final_sorting';
keys.Linus.sorted_neurons_sheetname      ='final_sorting';
keys.Curius.sorted_neurons_sheetname     ='final_sorting';
keys.Tesla.sorted_neurons_sheetname      ='final_sorting';
keys.Cornelius.sorted_neurons_sheetname  ='final_sorting';
keys.Magnus.sorted_neurons_sheetname     ='final_sorting';
keys.TDT_brain.sorted_neurons_sheetname  ='final_sorting';
keys.Bacchus.sorted_neurons_sheetname    ='final_sorting';

keys.Flaffus.filelist_formatted         ={};
keys.Linus.filelist_formatted           ={};
keys.Curius.filelist_formatted          ={};
keys.Tesla.filelist_formatted           ={};
keys.Cornelius.filelist_formatted       ={};
keys.Magnus.filelist_formatted          ={};
keys.TDT_brain.filelist_formatted       ={};
keys.Bacchus.filelist_formatted          ={};

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
keys.plot.single_cells                  =1;         % perform single cell plotting
keys.plot.waveforms                     =1;         % plot the waveform summary plots per session
keys.plot.polars_on_extra_figure        =0;         % recommended if there are too many conditions and the single cell heatmap plots are too small/crowded
keys.plot.eye_hand_traces               =1;         % Incredible performance booster if turned off
keys.plot.average_PSTH_line             =0;         % One PSTH line on top that represents the average of all others
keys.plot.average_heat_maps             =0;
keys.plot.export                        =1;         % save plots as pdfs, you typically want this
keys.plot.events                        =1:100;     % select events that should be plotted (vertical lines) on all PSTH like plots
keys.plot.population_PSTH_legends       =1;         % if population legends should be plotted or not
keys.plot.cell_count_legends            =1;         % if population legends should be plotted or not
keys.plot.scatter_legends               =1;         % if population legends should be plotted or not

% ylimits
%keys.plot.FR_max_for_ylim               =50;
keys.plot.trials_max_for_ylim           =50;
keys.plot.excentricity_max_for_ylim     =30;
keys.plot.eyetrace_factor               =0.5;
keys.plot.hndtrace_factor               =0.5;

%% are these used?
keys.width.PSTH_perpos          =0.5;
keys.width.raster               =0.1;
keys.width.PSTH_summary         =1;

%% colors
% traces
keys.colors.eye_ver         =[0.8 0 0];
keys.colors.eye_hor         =[1 0 0];
keys.colors.rhd_ver         =[0 0.8 0];
keys.colors.rhd_hor         =[0 1 0];
keys.colors.lhd_ver         =[0 0 0.8];
keys.colors.lhd_hor         =[0 0 1];

% cell count colors
keys.colors.NO_AN   =[255 255 255];
keys.colors.NO_TU   =[128 128 128];
keys.colors.EP_EN   =[239  65  54];
keys.colors.EP_BI   =[106 189  69];
keys.colors.EP_SU   =[0   174 239];
keys.colors.CR      =[255   0   0];
keys.colors.UC      =[127   0   0];

keys.colors.NH_IN   =[0 255 0];
keys.colors.NH_CH   =[0 128 0];
keys.colors.IH_IN   =[64 0 255];
keys.colors.IH_CH   =[32 0 128];
keys.colors.CH_IN   =[255 255 0];
keys.colors.CH_CH   =[128 128 0];
keys.colors.IS_IN   =[0 255 255];
keys.colors.IS_CH   =[0 128 128];
keys.colors.CS_IN   =[255 0 64];
keys.colors.CS_CH   =[128 0 32];

% single cell PSTH colors per position
keys.colors.LH_IN   =[64 0 255];
keys.colors.LH_CH   =[32 0 128];
keys.colors.RH_IN   =[255 255 0];
keys.colors.RH_CH   =[128 128 0];

% single cell PSTH colors per hemifield
keys.colors.NH_RS_IN=[255 0 64];
keys.colors.NH_RS_CH=[128 0 32];
keys.colors.LH_RS_IN=[255 0 255];
keys.colors.LH_RS_CH=[128 0 128];
keys.colors.RH_RS_IN=[255 128 0];
keys.colors.RH_RS_CH=[128 64 0];
keys.colors.NH_LS_IN=[0 255 255];
keys.colors.NH_LS_CH=[0 128 128];
keys.colors.LH_LS_IN=[0 128 255];
keys.colors.LH_LS_CH=[0 64 128];
keys.colors.RH_LS_IN=[0 255 0];
keys.colors.RH_LS_CH=[0 128 0];

% population contra ipsi PSTH colors
keys.colors.NH_CS_IN=[255 0 64];
keys.colors.NH_CS_CH=[128 0 32];
keys.colors.IH_CS_IN=[255 0 255];
keys.colors.IH_CS_CH=[128 0 128];
keys.colors.CH_CS_IN=[255 128 0];
keys.colors.CH_CS_CH=[128 64 0];
keys.colors.NH_IS_IN=[0 255 255];
keys.colors.NH_IS_CH=[0 128 128];
keys.colors.IH_IS_IN=[0 128 255];
keys.colors.IH_IS_CH=[0 64 128];
keys.colors.CH_IS_IN=[0 255 0];
keys.colors.CH_IS_CH=[0 128 0];
%% these need to be fixed somehow
keys.colors.NH_VS_IN=[150 150 150];
keys.colors.NH_VS_CH=[80 80 80];
keys.colors.IH_VS_IN=[64 0 255];
keys.colors.IH_VS_CH=[32 0 128];
keys.colors.CH_VS_IN=[255 255 0];
keys.colors.CH_VS_CH=[128 128 0];

% population effector colors

keys.colors.EF_SA   =[0 255  0];
keys.colors.EF_RE   =[0 255  0];
keys.colors.EF_FG   =[0 0  255];

%% overlapping tt and cc

% keys.cc.instructed_choice           ={'in','ch'}; %to loop through them
keys.colors.per_monkey          =[0 1 0; 1 0 0];
keys.colors.fix_offset          =[236 32 38; 16 159 218; 247 148 36]/255;
% keys.tt.instructed_choice           ='in';


%% tuning table readout options (excluding particular subsets)

keys.tt.combine_tuning_properties   ={'place_name_here'}; %% additional table entry from combining columns
keys.tt.epoch_criterion             ='none';
keys.tt.space_criterion             ='none';
keys.tt.hands_criterion             ='none';
keys.tt.SXH_criterion               ='none';
keys.tt.perturbations               = 0;
keys.tt.choices                     = 0;
keys.tt.hands                       = 0;
keys.tt.IC_for_criterion            = 'in';
keys.tt.trial_criterion_in          ='per_position';
keys.tt.trial_criterion_ch          ='per_hemifield';
keys.tt.selection                   ={};
keys.tt.unselect                    ={};
    keys.tt.type_effectors={'Msac'};

%% population
keys.sct=struct([]);
keys.ccs=struct([]);
keys.ons=struct([]);
keys.pop=struct([]);
keys.hst=struct([]);

%% EPOCH SETTINGS %%
%% For each task type seperately, analysis epochs are defined
%%

%% Fixation only type
keys.EPOCHS_PER_TYPE{1}={...
    'INI',      2,  -0.4,	-0.1,   'INI';...
    'Facq',     3,	0.05,	0.15,   'INI';...
    'Fhol',     4,	-0.3,	0,      'INI';...
    'PreS',     60,	-0.1,	-0.01,  'INI';...
    'PeriS',    60,	-0.01,	0.05,   'INI';...
    'PreR',     62,	-0.3, 	-0.01,  'INI';...
    'PeriR',	62,	-0.05, 	0.15,   'INI';...
    'Thol',     20,	-0.3,	0,      'INI';...
    };
keys.WINDOWS_PER_TYPE{1}={...
    'Aquisition',   3,	-0.5,	0.4;...
    'Fixation',     4,	-1.2,   0.17;...
    'Saccade',      60,	-0.1,   0.6;...
    'Reach',        62,	-0.35,  0.7;...
    };
keys.ANOVAS_PER_TYPE(1).epoch={'INI' 'Facq';...
    'Facq' 'Fhol';...
    'Fhol' 'PreS';...
    'Fhol' 'PeriS';...
    'Fhol' 'PreR';...
    'Fhol' 'PeriR';...
    'Fhol' 'Thol'};

keys.ANOVAS_PER_TYPE(1).spaceLR            ={'Facq','Fhol','PreS','PeriS','PreR','PeriR','Thol'}';
keys.ANOVAS_PER_TYPE(1).positions          ={'Facq','Fhol','PreS','PeriS','PreR','PeriR','Thol'}';
keys.ANOVAS_PER_TYPE(1).hands              ={'Facq','Fhol','PreS','PeriS','PreR','PeriR','Thol'}';
keys.ANOVAS_PER_TYPE(1).SxH                ={'Facq','Fhol','PreS','PeriS','PreR','PeriR','Thol'}';
keys.ANOVAS_PER_TYPE(1).main               ={'Facq','Fhol','PreS','PeriS','PreR','PeriR','Thol'}';

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
keys.ANOVAS_PER_TYPE(2).spaceLR            ={'Cue','PreS','PeriS','Tacq','Thol'}';
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
keys.ANOVAS_PER_TYPE(3).spaceLR            ={'Cue','MemE','MemL','PreS','PeriS','TIhol','Thol'}';
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
keys.ANOVAS_PER_TYPE(4).spaceLR            ={'INI','Facq','Fhol','Cue','EDel','Del','PreS','PeriS','PostS','PreR','PeriR','PostR','Thol'}';
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
keys.ANOVAS_PER_TYPE(5).spaceLR            ={'Cue','MemE','PreS','PeriS','Thol'}';
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
keys.ANOVAS_PER_TYPE(6).spaceLR            ={'Cue','MemE','PeriS','Exp','Thol'}';
keys.ANOVAS_PER_TYPE(6).positions          ={'Cue','MemE','PeriS','Exp','Thol'}';
keys.ANOVAS_PER_TYPE(6).hands              ={'Facq','Fhol','Cue','MemE','PeriS','Exp','Thol'}';
keys.ANOVAS_PER_TYPE(6).SxH                ={'Cue','MemE','PreS','Exp','Thol'}';
keys.ANOVAS_PER_TYPE(6).main              ={'Facq','Fhol','Cue','MemE','PeriS','Exp','Thol'}';

