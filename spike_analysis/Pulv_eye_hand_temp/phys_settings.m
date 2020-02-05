%keys.pdf_folder='v20161108_INI_extended_examples';
%keys.pdf_folders={'v20161208interleaved','v20161208blocked','v20161208smalltargets','v20161208largetargets'};
keys.pdf_folders={'v20161208blocked','v20161208smalltargets','v20161208largetargets','v20161208interleaved'};
% keys.pdf_folders={'v20161207b'};
% keys.pdf_folder='v20161207b';
keys.filelist_formatted={};
%keys.pdf_folders={'L_bha','L_iha','L_sst','L_lst'};
%keys.pdf_folders={'L_sst','L_lst','L_all'};

keys.effectors_on_same_figure   =1;
keys.extend_population          =1;
keys.plot_FR_separately         =1;

%% to check carefully
keys.combine_monkeys            =0;
keys.case_summaries         ={'hands'};
keys.FR_subtract_baseline   =0;

%% computation settings
keys.epochs_to_plot_PSTH{4}                 ={'Facq','Cue','PeriS','PeriR'}';

keys.effectors                  =[3,4,6];
keys.reach_hand                 =[1,2];
keys.types                      =[4];
keys.targets                ={'dPulv'};
keys.monkeys                ={'Linus'};%,'Flaffus'};%,'Flaffus'};%'Flaffus',
keys.Flaffus.date           ='[20160203 20161206]';
keys.Linus.date             ='[20160203 20160606]';

keys.cal.stablity       =[1];
keys.cal.single_rating  =[1,2];

%% cell count settings
keys.cell_count_factors                 ={'epoch','space','hand'};
keys.cc.conditions_to_plot              ={'Dcfr','Ddre','Ddsa'};
keys.cc.only_both_hands                 =1; 
keys.tt.SXH_criterion               ='SXH SXE HXE SHE';

%% population PSTH settings
cc=0;
% 1
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'ungrouped';
keys.pop(cc).conditions_to_plot         = {'Ddre','Ddsa','Dcfr'}; %%this wont work any more
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).only_both_hands            = 1;
% 2
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_PreR_spaceLR_Ddre_han';
keys.pop(cc).conditions_to_plot     	= {'Ddre'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).only_both_hands            = 1;
% 3
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_PeriR_spaceLR_Ddre_han';
keys.pop(cc).conditions_to_plot     	= {'Ddre'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).only_both_hands            = 1;
% 4
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_PostR_spaceLR_Ddre_han';
keys.pop(cc).conditions_to_plot     	= {'Ddre'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).only_both_hands            = 1;

% 5
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_PreR_hands_Ddre_han';
keys.pop(cc).conditions_to_plot      	= {'Ddre'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).only_both_hands            = 1;
% 6
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_PeriR_hands_Ddre_han';
keys.pop(cc).conditions_to_plot      	= {'Ddre'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).only_both_hands            = 1;
% 7
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_PostR_hands_Ddre_han';
keys.pop(cc).conditions_to_plot         = {'Ddre'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).only_both_hands            = 1;

% 8
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_PreS_spaceLR_Ddsa_han';
keys.pop(cc).conditions_to_plot         = {'Ddsa'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).only_both_hands            = 1;
% 9
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_PeriS_spaceLR_Ddsa_han';
keys.pop(cc).conditions_to_plot     	= {'Ddsa'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).only_both_hands            = 1;
% 10
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_PostS_spaceLR_Ddsa_han';
keys.pop(cc).conditions_to_plot     	= {'Ddsa'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).only_both_hands            = 1;

% 11
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_PreS_hands_Ddsa_han';
keys.pop(cc).conditions_to_plot      	= {'Ddsa'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).only_both_hands            = 1;
% 12
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_PeriS_hands_Ddsa_han';
keys.pop(cc).conditions_to_plot     	= {'Ddsa'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).only_both_hands            = 1;
% 13
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_PostS_hands_Ddsa_han';
keys.pop(cc).conditions_to_plot       	= {'Ddsa'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).only_both_hands            = 1;

% 
% 
% %% overlapp   in enhancement/supression
% % 64
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_PreR_epoch_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% keys.pop(cc).Selection                          = {'in_NH_PreS_epoch_Ddsa_han','su'};
% % 65
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_PeriR_epoch_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% keys.pop(cc).Selection                          = {'in_NH_PeriS_epoch_Ddsa_han','su'};
% % 66
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_PeriR_epoch_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% keys.pop(cc).Selection                          = {'in_NH_PeriS_epoch_Ddsa_han','en'};
% % 67
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_PreR_epoch_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% keys.pop(cc).Selection                          = {'in_NH_PreS_epoch_Ddsa_han','en'};
% 
% %% Cue/Del period and response period NH
% % 68
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_Cue_epoch_Ddsa_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddsa'};
% % 69
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_Del_epoch_Ddsa_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddsa'};
% % 70
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_PreS_epoch_Ddsa_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddsa'};
% % 71
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_PeriS_epoch_Ddsa_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddsa'};
% % 72
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_PostS_epoch_Ddsa_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddsa'};
% 
% 
% % 73
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_Cue_epoch_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% % 74
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_Del_epoch_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% % 75
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_PreR_epoch_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% % 76
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_PeriR_epoch_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% % 78
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_PostR_epoch_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% 
% % 79
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_Cue_epoch_Dcfr_han';
% keys.pop(cc).conditions_to_plot                 = {'Dcfr'};
% % 80
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_EDel_epoch_Dcfr_han';
% keys.pop(cc).conditions_to_plot                 = {'Dcfr'};
% % 81
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_PreR_epoch_Dcfr_han';
% keys.pop(cc).conditions_to_plot                 = {'Dcfr'};
% % 82
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_PeriR_epoch_Dcfr_han';
% keys.pop(cc).conditions_to_plot                 = {'Dcfr'};
% % 83
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_PostR_epoch_Dcfr_han';
% keys.pop(cc).conditions_to_plot                 = {'Dcfr'};
% 
% %% same, contra hand
% % 84
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_CH_Cue_epoch_Ddsa_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddsa'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 85
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_CH_Del_epoch_Ddsa_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddsa'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 86
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_CH_PreS_epoch_Ddsa_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddsa'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 87
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_CH_PeriS_epoch_Ddsa_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddsa'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 88
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_CH_PostS_epoch_Ddsa_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddsa'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 89
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_CH_Cue_epoch_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 90
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_CH_Del_epoch_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 91
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_CH_PreR_epoch_Ddre_han';
% keys.pop(cc).epoch_RF     = 'Cue';
% keys.pop(cc).epoch_for_normalization            = 'PreR';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% % 92
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_CH_PeriR_epoch_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 93
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_CH_PostR_epoch_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 94
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_CH_Cue_epoch_Dcfr_han';
% keys.pop(cc).conditions_to_plot                 = {'Dcfr'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 95
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_CH_EDel_epoch_Dcfr_han';
% keys.pop(cc).conditions_to_plot                 = {'Dcfr'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 96
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_CH_PreR_epoch_Dcfr_han';
% keys.pop(cc).conditions_to_plot                 = {'Dcfr'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 97
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_CH_PeriR_epoch_Dcfr_han';
% keys.pop(cc).conditions_to_plot                 = {'Dcfr'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 98
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_CH_PostR_epoch_Dcfr_han';
% keys.pop(cc).conditions_to_plot                 = {'Dcfr'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% %% and ipsi hand
% % 99
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_IH_Cue_epoch_Ddsa_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddsa'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 100
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_IH_Del_epoch_Ddsa_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddsa'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 101
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_IH_PreS_epoch_Ddsa_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddsa'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 102
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_IH_PeriS_epoch_Ddsa_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddsa'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 103
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_IH_PostS_epoch_Ddsa_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddsa'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 104
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_IH_Cue_epoch_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 105
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_IH_Del_epoch_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 106
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_IH_PreR_epoch_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 107
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_IH_PeriR_epoch_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 108
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_IH_PostR_epoch_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 109
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_IH_Cue_epoch_Dcfr_han';
% keys.pop(cc).conditions_to_plot                 = {'Dcfr'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 110
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_IH_EDel_epoch_Dcfr_han';
% keys.pop(cc).conditions_to_plot                 = {'Dcfr'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 111
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_IH_PreR_epoch_Dcfr_han';
% keys.pop(cc).conditions_to_plot                 = {'Dcfr'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 112
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_IH_PeriR_epoch_Dcfr_han';
% keys.pop(cc).conditions_to_plot                 = {'Dcfr'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% % 113
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_IH_PostR_epoch_Dcfr_han';
% keys.pop(cc).conditions_to_plot                 = {'Dcfr'};
% keys.pop(cc).group_excluded          = {'','-','bi'};
% 
% 
% %% response period overlapping dataset
% % 114
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_PreR_epoch_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% keys.pop(cc).Selection                          = {'existing_Ddre_han','1';'existing_Ddsa_han','1'};
% % 115
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_PeriR_epoch_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% keys.pop(cc).Selection                          = {'existing_Ddre_han','1';'existing_Ddsa_han','1'};
% % 116
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_PostR_epoch_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% keys.pop(cc).Selection                          = {'existing_Ddre_han','1';'existing_Ddsa_han','1'};
% % 117
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_PreS_epoch_Ddsa_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddsa'};
% keys.pop(cc).Selection                          = {'existing_Ddre_han','1';'existing_Ddsa_han','1'};
% % 118
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_PeriS_epoch_Ddsa_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddsa'};
% keys.pop(cc).Selection                          = {'existing_Ddre_han','1';'existing_Ddsa_han','1'};
% % 119
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_PostS_epoch_Ddsa_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddsa'};
% keys.pop(cc).Selection                          = {'existing_Ddre_han','1';'existing_Ddsa_han','1'};
% 
% 
% %% epoch in saccades
% % 120
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_Cue_epoch_Ddsa_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddsa'};
% keys.pop(cc).Selection                          = {'existing_Ddre_han','1';'existing_Ddsa_han','1'};
% % 121
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_EDel_epoch_Ddsa_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddsa'};
% keys.pop(cc).Selection                          = {'existing_Ddre_han','1';'existing_Ddsa_han','1'};
% % 122
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_Del_epoch_Ddsa_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddsa'};
% keys.pop(cc).Selection                          = {'existing_Ddre_han','1';'existing_Ddsa_han','1'};
% 
% % 123
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_Cue_epoch_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% keys.pop(cc).Selection                          = {'existing_Ddre_han','1';'existing_Ddsa_han','1'};
% % 124
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_EDel_epoch_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% keys.pop(cc).Selection                          = {'existing_Ddre_han','1';'existing_Ddsa_han','1'};
% % 125
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_Del_epoch_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% keys.pop(cc).Selection                          = {'existing_Ddre_han','1';'existing_Ddsa_han','1'};
% 
% % 126
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_Cue_epoch_Dcfr_han';
% keys.pop(cc).conditions_to_plot                 = {'Dcfr'};
% keys.pop(cc).Selection                          = {'existing_Ddre_han','1';'existing_Ddsa_han','1'};
% % 127
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_EDel_epoch_Dcfr_han';
% keys.pop(cc).conditions_to_plot                 = {'Dcfr'};
% keys.pop(cc).Selection                          = {'existing_Ddre_han','1';'existing_Ddsa_han','1'};
% % 128
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_NH_Del_epoch_Dcfr_han';
% keys.pop(cc).conditions_to_plot                 = {'Dcfr'};
% keys.pop(cc).Selection                          = {'existing_Ddre_han','1';'existing_Ddsa_han','1'};
% 
% % 129
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_Del_hands_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% % 130
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_Del_hands_Dcfr_han';
% keys.pop(cc).conditions_to_plot                 = {'Dcfr'};
% % 131
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_INI_hands_Dcfr_han';
% keys.pop(cc).conditions_to_plot                 = {'Dcfr'};
% % 132
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_INI_hands_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% % 133
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_INI_hands_Ddsa_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddsa'};
% % 134
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_Facq_hands_Dcfr_han';
% keys.pop(cc).conditions_to_plot                 = {'Dcfr'};
% % 135
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_Facq_hands_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% % 136
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_Facq_hands_Ddsa_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddsa'};
% 
% 
% % 137
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_INI_hands_Dcfr_han';
% keys.pop(cc).conditions_to_plot                 = {'Dcfr'};
% keys.pop(cc).FR_subtract_baseline               = 1;
% % 138
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_INI_hands_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% keys.pop(cc).FR_subtract_baseline               = 1;
% % 139
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_INI_hands_Ddsa_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddsa'};
% keys.pop(cc).FR_subtract_baseline               = 1;
% % 140
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_Facq_hands_Dcfr_han';
% keys.pop(cc).conditions_to_plot                 = {'Dcfr'};
% keys.pop(cc).FR_subtract_baseline               = 1;
% % 141
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_Facq_hands_Ddre_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddre'};
% keys.pop(cc).FR_subtract_baseline               = 1;
% % 142
% cc=cc+1;
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'in_Facq_hands_Ddsa_han';
% keys.pop(cc).conditions_to_plot                 = {'Ddsa'};
% keys.pop(cc).FR_subtract_baseline               = 1;
% 