keys.pdf_folders                ={''};
keys.pdf_folder                 ='testtest';
keys.filelist_formatted         ={};

keys.effectors_on_same_figure   =0;
keys.extend_population          =0;
keys.plot_FR_separately         =0;
keys.contra_color               ='orange';

%% to check carefully
keys.combine_monkeys            =1;
keys.case_summaries         ={'options'};
keys.FR_subtract_baseline   =0;

%% computation settings
keys.effectors                  =[0];
keys.reach_hand                 =[0];
keys.types                      =[3];

keys.targets                ={'pdSTS_L','FST_L'};
keys.monkeys                ={'Cornelius'};
keys.Cornelius.date         ='[20160506 20160814]';


%% cell count settings
keys.cell_count_factors         ={'epoch','space'};
keys.cc.conditions_to_plot      ={'Msac'};


%% for cell counts in general, but could also be used for PSTH!?
keys.tt.epoch_criterion             ='SxE or epoch only';
keys.tt.space_criterion             ='interaction or space only';

%% population PSTH settings
cc=0;

cc=cc+1;% 1 % memory raw
keys.pop(cc).normalization                      = 'none';
keys.pop(cc).group_parameter                    = 'ungrouped';
keys.pop(cc).conditions_to_plot                 = {'Msac'};

cc=cc+1;% 2 memory baseline subtracted
keys.pop(cc).normalization                      = 'none';
keys.pop(cc).group_parameter                    = 'ungrouped';
keys.pop(cc).conditions_to_plot                 = {'Msac'};
keys.pop(cc).FR_subtract_baseline               = 1;

cc=cc+1;% 3
keys.pop(cc).normalization                      = 'none';
keys.pop(cc).group_parameter                    = 'in_Cue_spaceLR_Msac_opt';
keys.pop(cc).conditions_to_plot                 = {'Msac'};

cc=cc+1;% 4
keys.pop(cc).normalization                      = 'by_effector';
keys.pop(cc).group_parameter                    = 'in_Cue_spaceLR_Msac_opt';
keys.pop(cc).epoch_for_normalization            = 'Cue';
keys.pop(cc).conditions_to_plot                 = {'Msac'};

cc=cc+1;% 5
keys.pop(cc).normalization                      = 'none';
keys.pop(cc).group_parameter                    = 'in_MemE_spaceLR_Msac_opt';
keys.pop(cc).epoch_for_receptive_field_plot     = 'MemE';
keys.pop(cc).epoch_for_normalization            = 'MemE';
keys.pop(cc).conditions_to_plot                 = {'Msac'};

cc=cc+1;% 6
keys.pop(cc).normalization                      = 'by_effector';
keys.pop(cc).group_parameter                    = 'in_MemE_spaceLR_Msac_opt';
keys.pop(cc).epoch_for_receptive_field_plot     = 'MemE';
keys.pop(cc).epoch_for_normalization            = 'MemE';
keys.pop(cc).conditions_to_plot                 = {'Msac'};

cc=cc+1;% 7
keys.pop(cc).normalization                      = 'none';
keys.pop(cc).group_parameter                    = 'in_MemL_spaceLR_Msac_opt';
keys.pop(cc).epoch_for_receptive_field_plot     = 'MemL';
keys.pop(cc).epoch_for_normalization            = 'MemL';
keys.pop(cc).conditions_to_plot                 = {'Msac'};

cc=cc+1;% 8
keys.pop(cc).normalization                      = 'by_effector';
keys.pop(cc).group_parameter                    = 'in_MemL_spaceLR_Msac_opt';
keys.pop(cc).epoch_for_receptive_field_plot     = 'MemL';
keys.pop(cc).epoch_for_normalization            = 'MemL';
keys.pop(cc).conditions_to_plot                 = {'Msac'};

cc=cc+1;% 9
keys.pop(cc).normalization                      = 'none';
keys.pop(cc).group_parameter                    = 'in_PeriS_epoch_Msac_opt';
keys.pop(cc).epoch_for_receptive_field_plot     = 'PeriS';
keys.pop(cc).epoch_for_normalization            = 'PeriS';
keys.pop(cc).conditions_to_plot                 = {'Msac'};

cc=cc+1;% 10
keys.pop(cc).normalization                      = 'by_effector';
keys.pop(cc).group_parameter                    = 'in_PeriS_epoch_Msac_opt';
keys.pop(cc).epoch_for_receptive_field_plot     = 'PeriS';
keys.pop(cc).epoch_for_normalization            = 'PeriS';
keys.pop(cc).conditions_to_plot                 = {'Msac'};

cc=cc+1;% 11
keys.pop(cc).normalization                      = 'by_effector';
keys.pop(cc).group_parameter                    = 'in_PeriS_spaceLR_Msac_opt';
keys.pop(cc).epoch_for_receptive_field_plot     = 'PeriS';
keys.pop(cc).epoch_for_normalization            = 'PeriS';
keys.pop(cc).conditions_to_plot                 = {'Msac'};

cc=cc+1;% 12
keys.pop(cc).normalization                      = 'none';
keys.pop(cc).group_parameter                    = 'in_TIhol_spaceLR_Msac_opt';
keys.pop(cc).epoch_for_receptive_field_plot     = 'TIhol';
keys.pop(cc).epoch_for_normalization            = 'TIhol';
keys.pop(cc).conditions_to_plot                 = {'Msac'};

cc=cc+1;% 13
keys.pop(cc).normalization                      = 'by_effector';
keys.pop(cc).group_parameter                    = 'in_TIhol_spaceLR_Msac_opt';
keys.pop(cc).epoch_for_normalization            = 'TIhol';
keys.pop(cc).conditions_to_plot                 = {'Msac'};

cc=cc+1;% 14
keys.pop(cc).normalization                      = 'none';
keys.pop(cc).group_parameter                    = 'in_Thol_spaceLR_Msac_opt';
keys.pop(cc).conditions_to_plot                 = {'Msac'};

cc=cc+1;% 15
keys.pop(cc).normalization                      = 'by_effector';
keys.pop(cc).group_parameter                    = 'in_Thol_spaceLR_Msac_opt';
keys.pop(cc).epoch_for_normalization            = 'Thol';
keys.pop(cc).conditions_to_plot                 = {'Msac'};

%% response fields
cc=cc+1;% 16
keys.pop(cc).normalization                      = 'by_effector';
keys.pop(cc).group_parameter                    = 'in_NH_Cue_position_Msac_opt';
keys.pop(cc).epoch_for_receptive_field_plot     = 'Cue';
keys.pop(cc).epoch_for_normalization            = 'Cue';
keys.pop(cc).conditions_to_plot                 = {'Msac'};
keys.pop(cc).plot_RF                            = 1;

cc=cc+1;% 17
keys.pop(cc).normalization                      = 'by_effector';
keys.pop(cc).group_parameter                    = 'in_NH_PeriS_position_Msac_opt';
keys.pop(cc).epoch_for_receptive_field_plot     = 'PeriS';
keys.pop(cc).epoch_for_normalization            = 'PeriS';
keys.pop(cc).conditions_to_plot                 = {'Msac'};
keys.pop(cc).plot_RF                            = 1;

