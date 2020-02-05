keys.project_versions={''};
keys.project_version='test';
keys.filelist_formatted={};
keys.filelist_as_blocks=0;

%% to check carefully
keys.position_and_plotting_arrangements             ={'hand_choices'};
keys.cal.units_from_sorting_table       =1;                         % exclude units that are not in the sorting table

%% computation settings
keys.cal.effectors                  =[0,1,6];
keys.cal.reach_hand                 =[0,1,2];
keys.cal.types                      =[1,4];
% keys.cal.stablity                       =[0,1];                     % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the table
% keys.cal.single_rating                  =[1,2,3];                   % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the tabl
% keys.cal.SNR_rating                     =[1,2,3,4];                   % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the tabl

%% batching
keys.batching.combine_monkeys           =0;
keys.batching.targets                   ={'MIP'};
keys.batching.monkeys                   ={'Magnus'};
keys.Magnus.date                        ='[20180809 20180809]';

%% cell count settings
keys.cc.factors                         ={'epoch','space','hand'};
keys.cc.conditions_to_plot              ={'Dcfr','Ddre','Ddsa'};

%% epochs
keys.EPOCHS_PER_TYPE{4}={...
    'INI',      2,	-0.4,	-0.1,   'INI';...
    'Facq',     3,	-0.4,	-0.1,   'INI';...
    'Fhol',     6,	-0.3,	0,      'INI';...
    'Cue',      6,	0.05,   0.15,   'INI';...
    'PreR',     62,	-0.3, 	-0.01,  'INI';...
    'PeriR',	62,	-0.05, 	0.15,   'INI';...
    'PostR',	63,	0.05,   0.2,    'INI';...
    'Thol',     20,	-0.3,   0,      'INI';...
    };

keys.WINDOWS_PER_TYPE{4}={...
    'Fixation',     3,	-0.7,   0.3;...
    'Delay Period', 6,	-0.33,  1;...
    };  

keys.ANOVAS_PER_TYPE(4).epoch={'INI' 'Facq';...
    'INI' 'Fhol';...
    'Fhol' 'Cue';...
    'Fhol' 'PreR';...
    'Fhol' 'PeriR';...
    'Fhol' 'PostR';...
    'Fhol' 'Thol';...
    };
keys.ANOVAS_PER_TYPE(4).spaceLR            ={'INI','Facq','Fhol','Cue','PreR','PeriR','PostR','Thol'}';
keys.ANOVAS_PER_TYPE(4).positions          ={'INI','Facq','Fhol','Cue','PreR','PeriR','PostR','Thol'}';
keys.ANOVAS_PER_TYPE(4).hands              ={'INI','Facq','Fhol','Cue','PreR','PeriR','PostR','Thol'}';
keys.ANOVAS_PER_TYPE(4).SxH                ={'INI','Facq','Fhol','Cue','PreR','PeriR','PostR','Thol'}';
keys.ANOVAS_PER_TYPE(4).main                ={'INI','Facq','Fhol','Cue','PreR','PeriR','PostR','Thol'}';

%% population PSTH settings
cc=0;
% 1
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'ungrouped';
keys.pop(cc).conditions_to_plot         = {'Ddre'}; %%this wont work any more
keys.pop(cc).FR_subtract_baseline       = 1;

% 1
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'ungrouped';
keys.pop(cc).conditions_to_plot         = {'Ddsa'}; %%this wont work any more
keys.pop(cc).FR_subtract_baseline       = 1;

% 1
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'ungrouped';
keys.pop(cc).conditions_to_plot         = {'Dcfr'}; %%this wont work any more
keys.pop(cc).FR_subtract_baseline       = 1;

% 2
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_PreR_spaceLR_Ddre_han';
keys.pop(cc).conditions_to_plot     	= {'Ddre'};
keys.pop(cc).FR_subtract_baseline       = 1;

% 3
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_PeriR_spaceLR_Ddre_han';
keys.pop(cc).conditions_to_plot     	= {'Ddre'};
keys.pop(cc).FR_subtract_baseline       = 1;

% 4
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_PostR_spaceLR_Ddre_han';
keys.pop(cc).conditions_to_plot     	= {'Ddre'};
keys.pop(cc).FR_subtract_baseline       = 1;


% 5
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_PreR_hands_Ddre_han';
keys.pop(cc).conditions_to_plot      	= {'Ddre'};
keys.pop(cc).FR_subtract_baseline       = 1;

% 6
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_PeriR_hands_Ddre_han';
keys.pop(cc).conditions_to_plot      	= {'Ddre'};
keys.pop(cc).FR_subtract_baseline       = 1;

% 7
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_PostR_hands_Ddre_han';
keys.pop(cc).conditions_to_plot         = {'Ddre'};
keys.pop(cc).FR_subtract_baseline       = 1;


% 8
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_PreS_spaceLR_Ddsa_han';
keys.pop(cc).conditions_to_plot         = {'Ddsa'};
keys.pop(cc).FR_subtract_baseline       = 1;

% 9
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_PeriS_spaceLR_Ddsa_han';
keys.pop(cc).conditions_to_plot     	= {'Ddsa'};
keys.pop(cc).FR_subtract_baseline       = 1;

% 10
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_PostS_spaceLR_Ddsa_han';
keys.pop(cc).conditions_to_plot     	= {'Ddsa'};
keys.pop(cc).FR_subtract_baseline       = 1;


% 11
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_PreS_hands_Ddsa_han';
keys.pop(cc).conditions_to_plot      	= {'Ddsa'};
keys.pop(cc).FR_subtract_baseline       = 1;

% 12
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_PeriS_hands_Ddsa_han';
keys.pop(cc).conditions_to_plot     	= {'Ddsa'};
keys.pop(cc).FR_subtract_baseline       = 1;

% 13
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_PostS_hands_Ddsa_han';
keys.pop(cc).conditions_to_plot       	= {'Ddsa'};
keys.pop(cc).FR_subtract_baseline       = 1;

