keys.pdf_folder='v20161123';
keys.targets                ={'dPulv'};
keys.monkeys                ={'Linus','Curius'};
keys.Linus.date             ='[20151119 20151216]';%keys.date='[20151118 20151216]';??
keys.Curius.date            ='[20151204 20160115]';
keys.datasets               ={'Msac'};
%keys.case_summaries         ={'fixation','movement vectors','target location by origin'};
keys.case_summaries         ={'movement vectors'}; %,'movement vectors'
keys.cell_count_factors     ={'epoch','space'};
keys.FR_subtract_baseline   =0; %needs to be removed soon

keys.effectors                  =[0];
keys.reach_hand                 =[0];
keys.types                      =[3];
keys.effectors_on_same_figure   =0;
keys.plot_FR_separately         =0;


%% renamed Cue Ton
% keys.EPOCHS_PER_TYPE{2}={...
%     'INI',      2,	-0.5,	0,      -0.4,	-0.1,   'INI';...
%     'Facq',     3,	-0.4,   0.17,	0.05,	0.15,   'INI';...
%     'Fhol',     4,	-0.33,	0,      -0.3,	0,      'INI';...
%     'Cue',      4,	0,      0.12,	0.05,	0.12,   'INI';...
%     'PreS',     60,	-0.1,   0,      -0.1,	-0.01,  'INI';...
%     'PeriS',    60,	-0.01,	0.05,	-0.01,	0.05,   'INI';...
%     'Tacq',     5,	-0.05,	0.17,	0,      0.15,   'INI';...
%     'Thol',     20,	-0.33,	0.1,	-0.3,	0,      'INI';...
%     };


%
% keys.EPOCHS_PER_TYPE{2}={...
%     'INI',      2,	-0.5,	0,      -0.4,	-0.1,   'INI';...
%     'Facq',     3,	-0.4,   0.17,	0.05,	0.15,   'INI';...
%     'Fhol',     4,	-0.33,	0,      -0.3,	0,      'INI';...
%     'Ton',      4,	0,      0.12,	0.05,	0.12,   'INI';...
%     'PreS',     60,	-0.1,   0,      -0.1,	-0.01,  'INI';...
%     'PeriS',    60,	-0.01,	0.05,	-0.01,	0.05,   'INI';...
%     'Tacq',     5,	-0.05,	0.17,	0,      0.15,   'INI';...
%     'Thol',     20,	-0.33,	0.1,	-0.3,	0,      'INI';...
%     };
% %keys.epochs_for_PSTH{2}={'INI','Facq','Fhol','Ton','PreS','PeriS','Tacq','Thol'}';
% keys.epoch_comparisons{2}={'INI' 'Facq';...
%     'Facq' 'Fhol';...
%     'Fhol' 'Ton';...
%     'Fhol' 'PreS';...
%     'Fhol' 'PeriS';...
%     'Fhol' 'Tacq';...
%     'Fhol' 'Thol'};
%
% keys.epochs_for_multicomparison{2}          ={'INI','Facq','Fhol','Ton','PreS','PeriS','Tacq','Thol'}';
% keys.epochs_spaceLR_multicomp{2}            ={'Ton','PreS','PeriS','Tacq','Thol'}';
% keys.epochs_choice_multicomp{2}             ={'Ton','PreS','PeriS','Tacq','Thol'}';
% keys.epochs_hands_multicomp{2}              ={'Facq','Fhol','Ton','PreS','PeriS','Tacq','Thol'}';
% keys.epochs_SxH_multicomp{2}   ={'Ton','PreS','PeriS','Tacq','Thol'}';

% case 500
%                             keys.population_normalization           = 'by_effector';
%                             keys.population_group_parameter         = 'in_Cue_spaceLR_Msac_fix';
%                             keys.epoch_for_receptive_field_plot     = 'Cue';
%                             keys.epoch_for_normalization            = 'Cue';
%                             keys.conditions_to_plot                 = {'Msac'};
%                             keys.plot_RF                            = 1;
%                         case 501
%                             keys.population_normalization           = 'by_effector';
%                             keys.population_group_parameter         = 'in_Cue_spaceLR_Msac_mov';
%                             keys.epoch_for_receptive_field_plot     = 'Cue';
%                             keys.epoch_for_normalization            = 'Cue';
%                             keys.conditions_to_plot                 = {'Msac'};
%                             keys.plot_RF                            = 1;
%                         case 502
%                             keys.population_normalization           = 'by_effector';
%                             keys.population_group_parameter         = 'in_Cue_spaceLR_Msac_tar';
%                             keys.epoch_for_receptive_field_plot     = 'Cue';
%                             keys.epoch_for_normalization            = 'Cue';
%                             keys.conditions_to_plot                 = {'Msac'};
%                             keys.plot_RF                            = 1;
%                         case 503
%                             keys.population_normalization           = 'by_effector';
%                             keys.population_group_parameter         = 'in_PeriS_spaceLR_Msac_fix';
%                             keys.epoch_for_receptive_field_plot     = 'PeriS';
%                             keys.epoch_for_normalization            = 'PeriS';
%                             keys.conditions_to_plot                 = {'Msac'};
%                             keys.plot_RF                            = 0;
%                         case 504
%                             keys.population_normalization           = 'by_effector';
%                             keys.population_group_parameter         = 'in_PeriS_spaceLR_Msac_mov';
%                             keys.epoch_for_receptive_field_plot     = 'PeriS';
%                             keys.epoch_for_normalization            = 'PeriS';
%                             keys.conditions_to_plot                 = {'Msac'};
%                             keys.plot_RF                            = 0;
%                         case 505
%                             keys.population_normalization           = 'by_effector';
%                             keys.population_group_parameter         = 'in_PeriS_spaceLR_Msac_tar';
%                             keys.epoch_for_receptive_field_plot     = 'PeriS';
%                             keys.epoch_for_normalization            = 'PeriS';
%                             keys.conditions_to_plot                 = {'Msac'};
%                             keys.plot_RF                            = 0;
%                         case 506
%                             keys.population_normalization           = 'by_effector';
%                             keys.population_group_parameter         = 'in_TIhol_spaceLR_Msac_fix';
%                             keys.epoch_for_receptive_field_plot     = 'TIhol';
%                             keys.epoch_for_normalization            = 'TIhol';
%                             keys.conditions_to_plot                 = {'Msac'};
%                             keys.plot_RF                            = 0;
%                         case 507
%                             keys.population_normalization           = 'by_effector';
%                             keys.population_group_parameter         = 'in_TIhol_spaceLR_Msac_mov';
%                             keys.epoch_for_receptive_field_plot     = 'TIhol';
%                             keys.epoch_for_normalization            = 'TIhol';
%                             keys.conditions_to_plot                 = {'Msac'};
%                             keys.plot_RF                            = 0;
%                         case 508
%                             keys.population_normalization           = 'by_effector';
%                             keys.population_group_parameter         = 'in_TIhol_spaceLR_Msac_tar';
%                             keys.epoch_for_receptive_field_plot     = 'TIhol';
%                             keys.epoch_for_normalization            = 'TIhol';
%                             keys.conditions_to_plot                 = {'Msac'};
%                             keys.plot_RF                            = 0;
                            
                            