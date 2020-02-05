%keys.project_versions={'v20161218blocked','v20161218interleaved'}; 
keys.project_versions={''}; 
keys.project_version='test_20170327'; 
 
%% to check carefully 
keys.position_and_plotting_arrangements         ={'hands'}; 
 
%% computation settings 
keys.cal.stablity       =[1]; 
keys.cal.single_rating  =[1,2]; 
keys.cal.effectors      =[3,4,6]; 
keys.cal.reach_hand     =[1,2]; 
keys.cal.types          =[4]; 
 
%% batching 
keys.batching.combine_monkeys       =0; 
keys.batching.targets               ={'dPulv'}; 
keys.batching.monkeys               ={'Flaffus'}; 
keys.Flaffus.date                   ='[20160617 20160617]'; 
%keys.Flaffus.date                  ='[20160203 20161206]'; 
keys.Linus.date                     ='[20160203 20160606]'; 
 
%% cell count settings 
keys.cc.factors                     ={'epoch','space','hand'}; 
keys.cc.conditions_to_plot          ={'Dcfr','Ddre','Ddsa'}; 
keys.tt.SXH_criterion               ='SXH SXE HXE SHE'; 
 
%% population PSTH settings 
cc=0; 
% 1 
cc=cc+1; 
keys.pop(cc).normalization              = 'none'; 
keys.pop(cc).group_parameter            = 'ungrouped'; 
keys.pop(cc).conditions_to_plot         = {'Ddre','Ddsa','Dcfr'}; %%this wont work any more 
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
 
 
 
 
 
 
 
 
 
