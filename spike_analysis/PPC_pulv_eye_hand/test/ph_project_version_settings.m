keys.project_versions={''};
keys.project_version='test';
keys.filelist_formatted={};

%% to check carefully
%keys.position_and_plotting_arrangements             ={'hands'};
keys.position_and_plotting_arrangements             ={'hands_inactivation'};

%% computation settings
keys.cal.effectors                  =[3,4,6];
keys.cal.reach_hand                 =[1,2];
keys.cal.types                      =[4];
keys.cal.units_from_sorting_table   =1;

%% batching
keys.batching.combine_monkeys           =0;
keys.batching.targets                   ={'MIP'};
keys.batching.monkeys                   ={'Linus','Tesla'};
keys.Tesla.date                         ='[20160217 20180101]';
%keys.Linus.date                         ='[20161103 20180101]';
keys.Linus.date                         ='[20170622 20180101]';

%% cell count settings
keys.cc.factors                         ={'epoch','space','hand'};
keys.cc.conditions_to_plot              ={'Dcfr','Ddre','Ddsa'};

%% epochs
keys.EPOCHS_PER_TYPE{4}={...
    'INI',      2,	-0.4,	-0.1,   'INI';...
    'Facq',     3,	-0.4,	-0.1,   'INI';...
    'Fhol',     6,	-0.3,	0,      'INI';...
    'Cue',      6,	0.05,   0.15,   'INI';...
    'EDel',     8, 	0.3,    0.6,   'INI';...
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
    'Fixation',     3,	-1.2,   0.17;...
    'Delay Period', 6,	-0.33,  1.35;...
    'Saccade',      60,	-0.1,   0.6;...
    'Reach',        62,	-0.35,  0.7;...
    };  

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

