keys.filelist_formatted={};

%% to check carefully
keys.position_and_plotting_arrangements         ={'options'};

%% computation settings
keys.cal.effectors                  =[0];
keys.cal.reach_hand                 =[0];
keys.cal.types                      =[2,3];


% keys.cal.types                          =[1,2,3,4,5,6];             % excluding trials with non-matching types
% keys.cal.reach_hand                     =[0,1,2];                   % excluding trials with non-matching reach_hand
% keys.cal.choice                         =[0,1];                     % excluding trials with non-matching chocie
% keys.cal.stablity                       =[0,1];                     % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the table
% keys.cal.single_rating                  =[1,2,3];                   % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the tabl
% keys.cal.SNR_rating                     =[1,2,3,4];                 % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the tabl
keys.cal.min_trials_per_condition       =60;                         % was actually at least 60 in total...  minimum trials per conditon (look at ph_arrange_positions to see how conditions are defined)

keys.tt.trial_criterion_in          ='total';

%% batching
keys.batching.combine_monkeys       =1;
keys.batching.targets               ={'dPulv_r'};%,'dPulv_l'};
keys.batching.monkeys               ={'Curius','Linus'};
keys.Curius.date                    ='[20150515 20150814]';
%keys.Curius.date                   ='[20150617 20150617]';
keys.Linus.date                     ='[20150508 20151030]';

%% cell count settings
keys.cc.factors                 ={'epoch','space'};
keys.cc.conditions_to_plot      ={'Msac','Vsac'};

%% for cell counts in general, but could also be used for PSTH!?
keys.tt.epoch_criterion             ='SxE or epoch only';
keys.tt.space_criterion             ='interaction or space only';

%% epochs
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
    'Visual onset', 4,	-0.8,   0.17;...
    'Saccade',      60,	-0.01,  0.8;...
    'T hold',       20,	-0.3,   0.1;...
    };

keys.ANOVAS_PER_TYPE(2).spaceLR            ={'Cue','PreS','PeriS','Tacq','Thol'}';
keys.ANOVAS_PER_TYPE(2).positions          ={'Cue','PreS','PeriS','Tacq','Thol'}';
keys.ANOVAS_PER_TYPE(2).hands              ={'Facq','Fhol','Cue','PreS','PeriS','Tacq','Thol'}';
keys.ANOVAS_PER_TYPE(2).SxH                ={'Cue','PreS','PeriS','Tacq','Thol'}';

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
    'Cue',          6,	-0.8,   0.78;...
    'Saccade',      60,	-0.7,   0.8;...
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

%% population PSTH settings
cc=0;

% cc=cc+1;% enhancement/suppression, POS-NEG
% keys.pop(cc).epoch_PF                = 'Thol';
% keys.pop(cc).epoch_RF                = 'Thol';
% keys.pop(cc).epoch_BL                = 'Fhol';
% keys.pop(cc).epoch_GB                = 'none';
% keys.pop(cc).epoch_for_normalization = 'Thol'; % for percent change, use this isntead!
% keys.pop(cc).normalization           = 'percent_change';
% keys.pop(cc).group_parameter         = 'ungrouped';
% keys.pop(cc).conditions_to_plot      = {'Vsac'};
% keys.pop(cc).position_and_plotting_arrangements      = {'options'};
% keys.pop(cc).FR_subtract_baseline    = 0;
% keys.pop(cc).ylim                    = [];
% keys.pop(cc).plot_per_position       = 1;
% keys.pop(cc).tt.choices              = 0;
% keys.pop(cc).tt.hands                = 0;
% keys.pop(cc).plot_RF                 = 0;
% 
% cc=cc+1;% enhancement/suppression, POS-NEG
% keys.pop(cc).epoch_PF                = 'Thol';
% keys.pop(cc).epoch_RF                = 'Thol';
% keys.pop(cc).epoch_BL                = 'Fhol';
% keys.pop(cc).epoch_GB                = 'none';
% keys.pop(cc).epoch_for_normalization = 'Thol'; % for percent change, use this isntead!
% keys.pop(cc).normalization           = 'percent_change';
% keys.pop(cc).group_parameter         = 'ungrouped';
% keys.pop(cc).conditions_to_plot      = {'Msac'};
% keys.pop(cc).position_and_plotting_arrangements      = {'options'};
% keys.pop(cc).FR_subtract_baseline    = 0;
% keys.pop(cc).ylim                    = [];
% keys.pop(cc).plot_per_position       = 1;
% keys.pop(cc).tt.choices              = 0;
% keys.pop(cc).tt.hands                = 0;
% keys.pop(cc).plot_RF                 = 0;

% cc=cc+1;% 1 ungrouped raw
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'ungrouped';
% keys.pop(cc).conditions_to_plot      = {'Vsac','Msac'};
% keys.pop(cc).ylim                    = [-0.5 3.5];
% 
% cc=cc+1;% 1 ungrouped baseline subtraction
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'ungrouped';
% keys.pop(cc).conditions_to_plot      = {'Vsac'};
% keys.pop(cc).FR_subtract_baseline    = 1;
% keys.pop(cc).ylim                    = [-0.5 3.5];
cc=cc+1;% 1 ungrouped baseline subtraction
keys.pop(cc).normalization           = 'none';
keys.pop(cc).group_parameter         = 'ungrouped';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).FR_subtract_baseline    = 1;
keys.pop(cc).ylim                    = [-0.5 3.5];
% cc=cc+1;% 1 ungrouped baseline subtraction
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'ungrouped';
% keys.pop(cc).conditions_to_plot      = {'Vsac'};
% keys.pop(cc).FR_subtract_baseline    = 0;
% keys.pop(cc).ylim                    = [-0.5 3.5];
% cc=cc+1;% 1 ungrouped baseline subtraction
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'ungrouped';
% keys.pop(cc).conditions_to_plot      = {'Msac'};
% keys.pop(cc).FR_subtract_baseline    = 0;
% keys.pop(cc).ylim                    = [-0.5 3.5];
% cc=cc+1;% 2 Cue spacial tuning
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'in_Cue_spaceLR_Msac_opt';
% keys.pop(cc).conditions_to_plot         = {'Msac'};
% keys.pop(cc).FR_subtract_baseline       = 1;
% keys.pop(cc).combine_tuning_properties  = {'S_ExS','in_spaceLR_main_Msac_opt','in_ExS_Msac_opt'};
% keys.pop(cc).unselect                   = {'S_ExS','-0'};
% keys.pop(cc).y_lim                      = [-4 21];
% 
% cc=cc+1;% 3 Peri-saccadic enhancement7suppression
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'in_NH_PeriS_epoch_Msac_opt';
% keys.pop(cc).epoch_RF                   = 'PeriS';
% keys.pop(cc).conditions_to_plot         = {'Msac'};
% keys.pop(cc).FR_subtract_baseline       = 1;
% keys.pop(cc).group_excluded             ={'','-','bi'};
% keys.pop(cc).combine_tuning_properties  ={'E_ExS','in_epoch_main_Msac_opt','in_ExS_Msac_opt'};
% keys.pop(cc).unselect                   ={'E_ExS','00'};
% keys.pop(cc).y_lim                      = [-4 8.5];
% 
% cc=cc+1;% 4 Target hold spacial tuning
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'in_Thol_spaceLR_Msac_opt';
% keys.pop(cc).epoch_RF                   = 'Thol';
% keys.pop(cc).conditions_to_plot         = {'Msac'};
% keys.pop(cc).FR_subtract_baseline       = 1;
% keys.pop(cc).combine_tuning_properties  ={'S_ExS','in_spaceLR_main_Msac_opt','in_ExS_Msac_opt'};
% keys.pop(cc).unselect                   ={'S_ExS','-0'};
% keys.pop(cc).y_lim                      = [-4 8.5];
% 
% cc=cc+1;% 5 receptive fields
% keys.pop(cc).normalization              = 'by_effector';
% keys.pop(cc).group_parameter            = 'in_NH_Cue_position_Msac_opt';
% keys.pop(cc).epoch_RF                  	= 'Cue';
% keys.pop(cc).epoch_for_normalization  	= 'Cue';
% keys.pop(cc).conditions_to_plot       	= {'Msac'};
% keys.pop(cc).FR_subtract_baseline       = 1;
% keys.pop(cc).plot_RF                 	= 1;
% keys.pop(cc).y_lim                      = [-0.3 1];
% 
% 
% 
