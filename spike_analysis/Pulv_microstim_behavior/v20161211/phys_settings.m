keys.pdf_folders={''};
%keys.pdf_folder='v20161207b';
keys.pdf_folder='v20161211';
keys.filelist_formatted={};

keys.effectors_on_same_figure   =0;
keys.extend_population          =0;
keys.plot_FR_separately         =0;

%% to check carefully
keys.combine_monkeys            =1;
%keys.case_summaries         ={'options'};
keys.case_summaries         ={'hands'};
keys.FR_subtract_baseline   =0;

%% computation settings
keys.effectors                  =[0];
keys.reach_hand                 =[0];
keys.types                      =[2,3];

keys.targets                ={'dPulv_r'};%,'dPulv_l'};
keys.monkeys                ={'Curius','Linus'};
%keys.Curius.date            ='[20150515 20150814]';
keys.Curius.date            ='[20150617 20150617]';
keys.Linus.date             ='[20150508 20151030]';


%% cell count settings
keys.cell_count_factors         ={'epoch','space'};
keys.cc.conditions_to_plot      ={'Msac','Vsac'};


%% for cell counts in general, but could also be used for PSTH!?
keys.tt.epoch_criterion             ='SxE or epoch only';
keys.tt.space_criterion             ='interaction or space only';

%% epochs
keys.EPOCHS_PER_TYPE{2}={...
    'INI',      2,	-0.5,	0,      -0.4,	-0.1,   'INI';...
    'Facq',     3,	-0.4,   0.17,	0.05,	0.15,   'INI';...
    'Fhol',     4,	-0.33,	0,      -0.3,	0,      'INI';...
    'Cue',      4,	-0.8,   0.17,	0.05,	0.15,   'INI';...
    'PreS',     60,	-0.1,   -0.01,  -0.1,	-0.01,  'INI';...
    'PeriS',    60,	-0.01,	0.22,	-0.01,	0.05,   'INI';...
    'Tacq',     5,	0,      0.2,	0,      0.15,   'INI';...
    'Thol',     20,	-0.3,	0.1,	-0.3,	0,      'INI';...
    };
keys.epoch_comparisons{2}={'INI' 'Facq';...
    'INI' 'Fhol';...
    'Fhol' 'Cue';...
    'Fhol' 'PreS';...
    'Fhol' 'PeriS';...
    'Fhol' 'Tacq';...
    'Fhol' 'Thol'};

keys.epochs_for_multicomparison{2}          ={'INI','Facq','Fhol','Cue','PreS','PeriS','Tacq','Thol'}';
keys.epochs_spaceLR_multicomp{2}            ={'Cue','PreS','PeriS','Tacq','Thol'}';
keys.epochs_choice_multicomp{2}             ={'Cue','PreS','PeriS','Tacq','Thol'}';
keys.epochs_hands_multicomp{2}              ={'Facq','Fhol','Cue','PreS','PeriS','Tacq','Thol'}';
keys.epochs_SxH_multicomp{2}                ={'Cue','PreS','PeriS','Tacq','Thol'}';
keys.epochs_to_plot_PSTH{2}                 ={'Cue','PeriS','Thol'}';

keys.EPOCHS_PER_TYPE{3}={...
    'INI',      2,	-0.5,	0,      -0.4,	-0.1,   'INI';...
    'Facq',     3,	-0.4,   0.17,	0.05,	0.15,   'INI';...
    'Fhol',     6,	-0.33,	0,      -0.3,	0,      'INI';...
    'Cue',      6,	-0.8,   0.78, 	0.05,   0.15,   'INI';...
    'MemE',     7, 	0,      0.5,    0,      0.2,    'INI';...
    'MemL',     9,	-0.5,   0,      -0.3, 	0,      'INI';...
    'PreS',     60,	-0.2,   0,      -0.1, 	-0.01,  'INI';...
    'PeriS',	60,	-0.7,   0.05,	-0.01, 	0.05,   'INI';...
    'TIhol',	10,	0,      0.2,    0,      0.1,    'INI';...
    'Thol',     20,	-0.6,   0,      -0.3,   0,      'INI';...
    };

keys.epoch_comparisons{3}={'INI' 'Facq';...
    'INI' 'Fhol';...
    'Fhol' 'Cue';...
    'Fhol' 'MemE';...
    'Fhol' 'MemL';...
    'Fhol' 'PreS';...
    'Fhol' 'PeriS';...
    'Fhol' 'TIhol';...
    'Fhol' 'Thol';...
    };

keys.epochs_for_multicomparison{3}          ={'INI','Facq','Fhol','Cue','MemE','MemL','PreS','PeriS','TIhol','Thol'}';
keys.epochs_spaceLR_multicomp{3}            ={'Cue','MemE','MemL','PreS','PeriS','TIhol','Thol'}';
keys.epochs_choice_multicomp{3}             ={'Cue','MemE','MemL','PreS','PeriS','TIhol','Thol'}';
keys.epochs_hands_multicomp{3}              ={'Facq','Fhol','Cue','MemE','MemL','PreS','PeriS','TIhol','Thol'}';
keys.epochs_SxH_multicomp{3}                ={'Cue','MemE','MemL','PreS','PeriS','TIhol','Thol'}';
keys.epochs_to_plot_PSTH{3}                 ={'INI','Cue','PeriS','Thol'}';


%% population PSTH settings
cc=0;

% cc=cc+1;% 1 ungrouped raw
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'ungrouped';
% keys.pop(cc).conditions_to_plot      = {'Vsac','Msac'};
% keys.pop(cc).ylim                    = [-0.5 3.5];

cc=cc+1;% 1 ungrouped baseline subtraction
keys.pop(cc).normalization           = 'none';
keys.pop(cc).group_parameter         = 'ungrouped';
keys.pop(cc).conditions_to_plot      = {'Vsac'};
keys.pop(cc).FR_subtract_baseline    = 1;
keys.pop(cc).ylim                    = [-0.5 3.5];
cc=cc+1;% 1 ungrouped baseline subtraction
keys.pop(cc).normalization           = 'none';
keys.pop(cc).group_parameter         = 'ungrouped';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).FR_subtract_baseline    = 1;
keys.pop(cc).ylim                    = [-0.5 3.5];
cc=cc+1;% 1 ungrouped baseline subtraction
keys.pop(cc).normalization           = 'none';
keys.pop(cc).group_parameter         = 'ungrouped';
keys.pop(cc).conditions_to_plot      = {'Vsac'};
keys.pop(cc).FR_subtract_baseline    = 0;
keys.pop(cc).ylim                    = [-0.5 3.5];
cc=cc+1;% 1 ungrouped baseline subtraction
keys.pop(cc).normalization           = 'none';
keys.pop(cc).group_parameter         = 'ungrouped';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).FR_subtract_baseline    = 0;
keys.pop(cc).ylim                    = [-0.5 3.5];


cc=cc+1;% 2 Cue spacial tuning
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_Cue_spaceLR_Msac_opt';
keys.pop(cc).conditions_to_plot         = {'Msac'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).combine_tuning_properties  = {'S_ExS','in_spaceLR_main_Msac_opt','in_ExS_Msac_opt'};
keys.pop(cc).unselect                   = {'S_ExS','-0'};
keys.pop(cc).y_lim                      = [-4 21];

cc=cc+1;% 3 Peri-saccadic enhancement7suppression
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_NH_PeriS_epoch_Msac_opt';
keys.pop(cc).epoch_RF                   = 'PeriS';
keys.pop(cc).conditions_to_plot         = {'Msac'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).group_excluded             ={'','-','bi'};
keys.pop(cc).combine_tuning_properties  ={'E_ExS','in_epoch_main_Msac_opt','in_ExS_Msac_opt'};
keys.pop(cc).unselect                   ={'E_ExS','00'};
keys.pop(cc).y_lim                      = [-4 8.5];

cc=cc+1;% 4 Target hold spacial tuning
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_Thol_spaceLR_Msac_opt';
keys.pop(cc).epoch_RF                   = 'Thol';
keys.pop(cc).conditions_to_plot         = {'Msac'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).combine_tuning_properties  ={'S_ExS','in_spaceLR_main_Msac_opt','in_ExS_Msac_opt'};
keys.pop(cc).unselect                   ={'S_ExS','-0'};
keys.pop(cc).y_lim                      = [-4 8.5];

cc=cc+1;% 5 receptive fields
keys.pop(cc).normalization              = 'by_effector';
keys.pop(cc).group_parameter            = 'in_NH_Cue_position_Msac_opt';
keys.pop(cc).epoch_RF                  	= 'Cue';
keys.pop(cc).epoch_for_normalization  	= 'Cue';
keys.pop(cc).conditions_to_plot       	= {'Msac'};
keys.pop(cc).plot_RF                 	= 1;
keys.pop(cc).y_lim                      = [0 1.6];

% cc=cc+1;% 2 Cue spacial tuning
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'in_Cue_spaceLR_Vsac_opt';
% keys.pop(cc).conditions_to_plot         = {'Vsac'};
% keys.pop(cc).FR_subtract_baseline       = 1;
% keys.pop(cc).combine_tuning_properties  = {'S_ExS','in_spaceLR_main_Vsac_opt','in_ExS_Vsac_opt'};
% keys.pop(cc).unselect                   = {'S_ExS','-0'};

% cc=cc+1;% 3 Peri-saccadic enhancement7suppression
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'in_NH_PeriS_epoch_Vsac_opt';
% keys.pop(cc).epoch_RF                   = 'PeriS';
% keys.pop(cc).conditions_to_plot         = {'Vsac'};
% keys.pop(cc).FR_subtract_baseline       = 1;
% keys.pop(cc).group_excluded             ={'','-','bi'};
% keys.pop(cc).combine_tuning_properties  ={'E_ExS','in_epoch_main_Vsac_opt','in_ExS_Vsac_opt'};
% keys.pop(cc).unselect                   ={'E_ExS','00'};
% 
% cc=cc+1;% 4 Target hold spacial tuning
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'in_Thol_spaceLR_Vsac_opt';
% keys.pop(cc).epoch_RF                   = 'Thol';
% keys.pop(cc).conditions_to_plot         = {'Vsac'};
% keys.pop(cc).FR_subtract_baseline       = 1;
% keys.pop(cc).combine_tuning_properties	={'S_ExS','in_spaceLR_main_Vsac_opt','in_ExS_Vsac_opt'};
% keys.pop(cc).unselect                   ={'S_ExS','-0'};


% 
% cc=cc+1;% 214
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'in_NH_Cue_position_Vsac_opt';
% keys.pop(cc).epoch_RF                   = 'Cue';
% keys.pop(cc).conditions_to_plot     	= {'Vsac'};
% keys.pop(cc).plot_RF                 	= 1;
% keys.pop(cc).FR_subtract_baseline    	= 1;
% 
% cc=cc+1;% 215
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'in_NH_Cue_position_Msac_opt';
% keys.pop(cc).epoch_RF                   = 'Cue';
% keys.pop(cc).conditions_to_plot         = {'Msac'};
% keys.pop(cc).plot_RF                    = 1;
% keys.pop(cc).FR_subtract_baseline       = 1;

% cc=cc+1;% 228
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'in_NH_Cue_position_Vsac_opt';
% keys.pop(cc).epoch_RF                   = 'Cue';
% keys.pop(cc).conditions_to_plot      	= {'Vsac'};
% keys.pop(cc).plot_RF                  	= 1;
% keys.pop(cc).FR_subtract_baseline     	= 1;
% keys.pop(cc).combine_tuning_properties 	={'S_ExS','in_spaceLR_main_Vsac_opt','in_ExS_Vsac_opt'};
% keys.pop(cc).unselect                   ={'S_ExS','-0'};

% cc=cc+1;% 229
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'in_NH_Cue_position_Msac_opt';
% keys.pop(cc).epoch_RF                   = 'Cue';
% keys.pop(cc).conditions_to_plot         = {'Msac'};
% keys.pop(cc).plot_RF                    = 1;
% keys.pop(cc).FR_subtract_baseline   	= 1;
% keys.pop(cc).combine_tuning_properties 	={'S_ExS','in_spaceLR_main_Msac_opt','in_ExS_Msac_opt'};
% keys.pop(cc).unselect                   ={'S_ExS','-0'};
% 
% cc=cc+1;% 2160
% keys.pop(cc).normalization              = 'by_effector';
% keys.pop(cc).group_parameter            = 'in_NH_Cue_position_Msac_opt';
% keys.pop(cc).epoch_RF                   = 'Cue';
% keys.pop(cc).conditions_to_plot         = {'Msac'};
% keys.pop(cc).plot_RF                  	= 1;
% keys.pop(cc).FR_subtract_baseline    	= 1;
% keys.pop(cc).y_lim                      = [0 1.6];
% 
% cc=cc+1;% 217
% keys.pop(cc).normalization              = 'by_effector';
% keys.pop(cc).group_parameter            = 'in_Cue_spaceLR_Msac_opt';
% keys.pop(cc).epoch_RF                   = 'Cue';
% keys.pop(cc).epoch_for_normalization   	= 'Cue';
% keys.pop(cc).conditions_to_plot       	= {'Msac'};
% keys.pop(cc).plot_RF                  	= 1;
% keys.pop(cc).combine_tuning_properties 	={'S_ExS','in_spaceLR_main_Msac_opt','in_ExS_Msac_opt'};
% keys.pop(cc).unselect                   ={'S_ExS','-0'};
% 
% cc=cc+1;% 219
% keys.pop(cc).normalization              = 'by_effector';
% keys.pop(cc).group_parameter            = 'in_Cue_spaceLR_Vsac_opt';
% keys.pop(cc).epoch_RF                   = 'Cue';
% keys.pop(cc).epoch_for_normalization   	= 'Cue';
% keys.pop(cc).conditions_to_plot        	= {'Vsac'};
% keys.pop(cc).plot_RF                  	= 1;
% keys.pop(cc).combine_tuning_properties 	={'S_ExS','in_spaceLR_main_Vsac_opt','in_ExS_Vsac_opt'};
% keys.pop(cc).unselect                   ={'S_ExS','-0'};


% cc=cc+1;% 3 %memory excluding non task responsive
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'ungrouped';
% keys.pop(cc).epoch_RF                   = 'Cue';
% keys.pop(cc).epoch_for_normalization  	= 'Cue';
% keys.pop(cc).conditions_to_plot       	= {'Msac'};
% keys.pop(cc).combine_tuning_properties	={'E_S_ExS','in_epoch_main_Msac_opt','in_spaceLR_main_Msac_opt','in_ExS_Msac_opt'};
% keys.pop(cc).unselect                   ={'E_S_ExS','0-0'};
% 
% cc=cc+1;% 4 %direct excluding non task responsive
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'ungrouped';
% keys.pop(cc).epoch_RF                   = 'Cue';
% keys.pop(cc).epoch_for_normalization  	= 'Cue';
% keys.pop(cc).conditions_to_plot       	= {'Vsac'};
% keys.pop(cc).combine_tuning_properties 	={'E_S_ExS','in_epoch_main_Vsac_opt','in_spaceLR_main_Vsac_opt','in_ExS_Vsac_opt'};
% keys.pop(cc).unselect                   ={'E_S_ExS','0-0'};
% 
% cc=cc+1;% 5 memory cue tuning
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'in_Cue_spaceLR_Msac_opt';
% keys.pop(cc).conditions_to_plot      	= {'Msac'};
% keys.pop(cc).FR_subtract_baseline      	= 1;
% 
% cc=cc+1;% 6 direct cue tuning
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'in_Cue_spaceLR_Vsac_opt';
% keys.pop(cc).conditions_to_plot         = {'Vsac'};
% keys.pop(cc).FR_subtract_baseline      	= 1;







