 keys.project_versions={''};
keys.project_version='singlecells';
keys.filelist_formatted={};

%% to check carefully
keys.position_and_plotting_arrangements             ={'hands_inactivation'};
keys.plot.single_cells                  =1;         % perform single cell plotting
keys.plot.polars_on_extra_figure        =1;
keys.plot.events                        =[2:10, 60, 61,62,63,20];

%% computation settings
keys.cal.units_from_sorting_table       =0;                         % exclude units that are not in the sorting table
keys.cal.effectors                  =[3,4,6];
keys.cal.reach_hand                 =[1,2];
keys.cal.types                      =[4];
keys.cal.stablity                       =[0,1];                     % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the table
keys.cal.single_rating                  =[1,2,3];                   % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the tabl
keys.cal.SNR_rating                     =[1,2,3,4];                   % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the tabl
keys.cal.min_trials_per_condition       =5;                        % about to be used !!
%keys.cal.min_spikes_per_unit            =50;                        % excluding units that have in total less spikes (workaround for sortcode assignment bug)


%% batching
keys.batching.combine_monkeys           =0;
keys.batching.targets                   ={'MIP_L','MIP_R','LIP_L','LIP_R'};
%keys.batching.targets                   ={'_L','_R'};
keys.batching.monkeys                   ={'Linus'};%,'Tesla'};
keys.Tesla.date                         ='[20160217 20180101]';
keys.Tesla.date                         ='[20161103 20180101]';

%% cell count settings
keys.cc.plot_types          ={'hands_inactivation'};%{'per_epoch','per_task','space_x_hand','hands_inactivation'};
keys.cc.factors             ={'space','hand','epoch'}; % only for per_epoch and per_task 
keys.cc.conditions_to_plot              ={'Ddre','Ddsa'};%'Dcfr',

keys.cc.hands=[1 2];
keys.cc.choices=[0];

keys.cc.epochsE.Ddre        ={'Fhol','Cue','EDel','Del','PreS','PeriS','PostS','PreR','PeriR','PostR','Thol'}';
keys.cc.epochsS.Ddre        ={'Cue','Del','PeriS','PeriR','Thol'}';
keys.cc.epochsE.Dcfr        ={'Fhol','Cue','EDel','Del','PreS','PeriS','PostS','PreR','PeriR','PostR','Thol'}';
keys.cc.epochsS.Dcfr        ={'Cue','Del','PeriS','PeriR','Thol'}';
keys.cc.epochsE.Ddsa        ={'Fhol','Cue','EDel','Del','PreS','PeriS','PostS','PreR','PeriR','PostR','Thol'}';
keys.cc.epochsS.Ddsa        ={'Cue','Del','PeriS','PeriR','Thol'}';


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

%% tuning onset plots
% cc=0;
% 
% cc=cc+1;% 1 ungrouped baseline subtraction
% ce=0;
% keys.ons(cc).comparisons_title       = 'Inactivation effects ipsi hand contra space';
% keys.ons(cc).group_parameter         = 'ungrouped';
% keys.ons(cc).conditions_to_plot      = {'Ddre'};
% keys.ons(cc).choices=0; %for cell exclusion
% keys.ons(cc).hands=[1 2]; %for cell exclusion
% ce=ce+1;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=1;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=1;
% keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[1];
% keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
% keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).choice{2}=0;
% keys.ons(cc).comparisons_per_effector(ce).perturbation{1}=[0];
% keys.ons(cc).comparisons_per_effector(ce).perturbation{2}=[1];
% keys.ons(cc).comparisons_per_effector(ce).order_onset={'Delay Period', -0.33, 1.35};
% keys.ons(cc).comparisons_per_effector(ce).colors=[[0 0 0]; [1 0 0]];
% keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
% keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Cue';
% ce=ce+1;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=1;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=1;
% keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[1];
% keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
% keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).choice{2}=0;
% keys.ons(cc).comparisons_per_effector(ce).perturbation{1}=[0];
% keys.ons(cc).comparisons_per_effector(ce).perturbation{2}=[1];
% keys.ons(cc).comparisons_per_effector(ce).order_onset={'Reach', -0.35, 0.7};
% keys.ons(cc).comparisons_per_effector(ce).colors=[[0 0 0]; [1 0 0]];
% keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='MemL';
% keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Reach';

%% population PSTH settings
cc=0;
% 1
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'ungrouped';
keys.pop(cc).conditions_to_plot         = {'Ddre'}; %%this wont work any more
keys.pop(cc).choices                    = [0];
keys.pop(cc).hands                      = [1 2];
keys.pop(cc).FR_subtract_baseline               = 1;
% 1
cc=cc+1;
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'ungrouped';
keys.pop(cc).conditions_to_plot         = {'Ddsa'}; %%this wont work any more
keys.pop(cc).choices                    = [0];
keys.pop(cc).hands                      = [1 2];
keys.pop(cc).FR_subtract_baseline               = 1;
cc=cc+1;
keys.pop(cc).normalization              = 'by_effector';
keys.pop(cc).group_parameter            = 'ungrouped';
keys.pop(cc).conditions_to_plot         = {'Ddre'}; %%this wont work any more
keys.pop(cc).epoch_for_normalization  	= 'INI';
keys.pop(cc).choices                    = [0];
keys.pop(cc).hands                      = [1 2];
keys.pop(cc).FR_subtract_baseline               = 0;
% 1
cc=cc+1;
keys.pop(cc).normalization              = 'by_effector';
keys.pop(cc).group_parameter            = 'ungrouped';
keys.pop(cc).conditions_to_plot         = {'Ddsa'}; %%this wont work any more
keys.pop(cc).epoch_for_normalization  	= 'INI';
keys.pop(cc).choices                    = [0];
keys.pop(cc).hands                      = [1 2];
keys.pop(cc).FR_subtract_baseline               = 0;
% % 2
% cc=cc+1;
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'in_PreR_spaceLR_Ddre_han';
% keys.pop(cc).conditions_to_plot     	= {'Ddre'};
% keys.pop(cc).FR_subtract_baseline       = 1;
% 
% % 3
% cc=cc+1;
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'in_PeriR_spaceLR_Ddre_han';
% keys.pop(cc).conditions_to_plot     	= {'Ddre'};
% keys.pop(cc).FR_subtract_baseline       = 1;
% 
% % 4
% cc=cc+1;
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'in_PostR_spaceLR_Ddre_han';
% keys.pop(cc).conditions_to_plot     	= {'Ddre'};
% keys.pop(cc).FR_subtract_baseline       = 1;
% 
% 
% % 5
% cc=cc+1;
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'in_PreR_hands_Ddre_han';
% keys.pop(cc).conditions_to_plot      	= {'Ddre'};
% keys.pop(cc).FR_subtract_baseline       = 1;
% 
% % 6
% cc=cc+1;
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'in_PeriR_hands_Ddre_han';
% keys.pop(cc).conditions_to_plot      	= {'Ddre'};
% keys.pop(cc).FR_subtract_baseline       = 1;
% 
% % 7
% cc=cc+1;
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'in_PostR_hands_Ddre_han';
% keys.pop(cc).conditions_to_plot         = {'Ddre'};
% keys.pop(cc).FR_subtract_baseline       = 1;
% 
% 
% % 8
% cc=cc+1;
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'in_PreS_spaceLR_Ddsa_han';
% keys.pop(cc).conditions_to_plot         = {'Ddsa'};
% keys.pop(cc).FR_subtract_baseline       = 1;
% 
% % 9
% cc=cc+1;
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'in_PeriS_spaceLR_Ddsa_han';
% keys.pop(cc).conditions_to_plot     	= {'Ddsa'};
% keys.pop(cc).FR_subtract_baseline       = 1;
% 
% % 10
% cc=cc+1;
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'in_PostS_spaceLR_Ddsa_han';
% keys.pop(cc).conditions_to_plot     	= {'Ddsa'};
% keys.pop(cc).FR_subtract_baseline       = 1;
% 
% 
% % 11
% cc=cc+1;
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'in_PreS_hands_Ddsa_han';
% keys.pop(cc).conditions_to_plot      	= {'Ddsa'};
% keys.pop(cc).FR_subtract_baseline       = 1;
% 
% % 12
% cc=cc+1;
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'in_PeriS_hands_Ddsa_han';
% keys.pop(cc).conditions_to_plot     	= {'Ddsa'};
% keys.pop(cc).FR_subtract_baseline       = 1;
% 
% % 13
% cc=cc+1;
% keys.pop(cc).normalization              = 'none';
% keys.pop(cc).group_parameter            = 'in_PostS_hands_Ddsa_han';
% keys.pop(cc).conditions_to_plot       	= {'Ddsa'};
% keys.pop(cc).FR_subtract_baseline       = 1;
% 
% 
%  
