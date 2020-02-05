keys.project_versions={''};
keys.project_version='Thesis';
keys.filelist_formatted={};

%% to check carefully
keys.position_and_plotting_arrangements         ={'options'};

%% computation settings
keys.cal.effectors                  =[0];
keys.cal.reach_hand                 =[0];
keys.cal.types                      =[3];

keys.cal.stablity                       =[1];                     % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the table
keys.cal.single_rating                  =[1,2,3];                   % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the tabl
keys.cal.SNR_rating                     =[1,2,3];                   % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the tabl


%% batching
keys.batching.combine_monkeys       =1;
keys.batching.targets               ={'dPulv'};
keys.batching.monkeys               ={'Linus'};%
keys.Curius.date                    ='[20150515 20150814]';
%keys.Curius.date                   ='[20150617 20150617]';
keys.Linus.date                     ='[20150508 20151030]';

%% cell count settings
keys.cc.hands               =[0];
keys.cc.choices             =[0];
keys.cc.factors                 ={'epoch','space'};
keys.cc.conditions_to_plot      ={'Msac'};
keys.cc.epochsE.Msac        ={'INI', 'Fhol','Cue','Cue2','MemE','MemL','PreS','PeriS','TIhol','Thol'}';
keys.cc.epochsS.Msac        ={'Cue','Cue2','MemE','MemL','PreS','PeriS','TIhol'}';

%% for cell counts in general, but could also be used for PSTH!?
keys.tt.epoch_criterion             ='SxE or epoch only';
keys.tt.space_criterion             ='interaction or space only';

%% epochs

keys.EPOCHS_PER_TYPE{3}={...
    'INI',      2,	-0.4,	-0.1,   'INI';...
    'Fhol',     6,	-0.3,	0,      'INI';...
    'Cue',      6,	0.07,   0.17,   'INI';...
    'MemL',     9,	-0.3, 	-0,     'INI';... 
    'PreS',     60,	-0.1, 	-0.01,  'INI';...
    'PeriS',	60,	-0.01,  0.05,   'INI';...
    'TIhol',	5,	-0.15,  0,      'INI';...
    'Tons',     5,	0,      0.1,    'INI';...
    'Thol',     5,	0.2,    0.5,    'INI';...
    };

keys.WINDOWS_PER_TYPE{3}={...
    'Initiation',   2,	-0.5,	0;...
    'Cue',          6,	-0.8,   0.78;...
    'Saccade',      60,	-0.7,   0.2;...
    'Hold',         5,	0,      0.5;...
    };

keys.ANOVAS_PER_TYPE(3).epoch={'INI' 'Fhol';...
    'Fhol' 'Cue';...
    'Fhol' 'MemL';...
    'MemL' 'PreS';...
    'MemL' 'PeriS';...
    'MemL' 'TIhol';...
    'MemL' 'Tons';...
    'MemL' 'Thol';...
    };

keys.ANOVAS_PER_TYPE(3).spaceLR            ={'Cue','MemL','PreS','Pre2','PeriS','Peri2','TIhol','Tons','Thol'}';
keys.ANOVAS_PER_TYPE(3).positions          ={'Cue','MemL','PreS','Pre2','PeriS','Peri2','TIhol','Tons','Thol'}';
keys.ANOVAS_PER_TYPE(3).hands              ={'Fhol','Cue','MemL','PreS','Pre2','PeriS','Peri2','TIhol','Tons','Thol'}';
keys.ANOVAS_PER_TYPE(3).SxH                ={'Cue','MemL','PreS','Pre2','PeriS','Peri2','TIhol','Tons','Thol'}';
keys.ANOVAS_PER_TYPE(3).main              ={'INI','Fhol','Cue','MemL','PreS','PeriS','Tons','Thol'}';



% %% tuning onset settings
% cc=0;
% 
% cc=cc+1;% 'Msac epoch tuning';
% ce=0;
% keys.ons(cc).comparisons_title       = 'Msac epoch tuning';
% keys.ons(cc).group_parameter         = 'ungrouped';
% keys.ons(cc).conditions_to_plot      = {'Msac'};
% 
% keys.ons(cc).choices=[0]; %for cell exclusion
% ce=ce+1;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
% keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1 1];
% keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[-1 1];
% keys.ons(cc).comparisons_per_effector(ce).choice{1}=[0];
% keys.ons(cc).comparisons_per_effector(ce).choice{2}=[0];
% keys.ons(cc).comparisons_per_effector(ce).order_onset={'Cue', 0, 0.28};
% keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.EP_EN/255; keys.colors.EP_SU/255];
% keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='Fhol';
% keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Cue';
% ce=ce+1;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
% keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1 1];
% keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[-1 1];
% keys.ons(cc).comparisons_per_effector(ce).choice{1}=[0];
% keys.ons(cc).comparisons_per_effector(ce).choice{2}=[0];
% keys.ons(cc).comparisons_per_effector(ce).order_onset={'Saccade', -0.20, 0};
% keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.EP_EN/255; keys.colors.EP_SU/255];
% keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='MemL';
% keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Saccade';
% ce=ce+1;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
% keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1 1];
% keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[-1 1];
% keys.ons(cc).comparisons_per_effector(ce).choice{1}=[0];
% keys.ons(cc).comparisons_per_effector(ce).choice{2}=[0];
% keys.ons(cc).comparisons_per_effector(ce).order_onset={'Hold', 0, 0.5};
% keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.EP_EN/255; keys.colors.EP_SU/255];
% keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='MemL';
% keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Hold';
% 
% cc=cc+1;% 1 ungrouped baseline subtraction
% ce=0;
% keys.ons(cc).comparisons_title       = 'Msac instructed space tuning';
% keys.ons(cc).group_parameter         = 'ungrouped';
% keys.ons(cc).conditions_to_plot      = {'Msac'};
% keys.ons(cc).choices=0; %for cell exclusion
% ce=ce+1;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
% keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1];
% keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
% keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).choice{2}=0;
% keys.ons(cc).comparisons_per_effector(ce).order_onset={'Cue', 0, 0.28};
% keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IS_IN/255; keys.colors.NH_CS_IN/255];
% keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
% keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Cue';
% ce=ce+1;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
% keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1];
% keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
% keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).choice{2}=0;
% keys.ons(cc).comparisons_per_effector(ce).order_onset={'Saccade', -0.20, 0};
% keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IS_IN/255; keys.colors.NH_CS_IN/255];
% keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
% keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Saccade';
% ce=ce+1;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
% keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1];
% keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
% keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).choice{2}=0;
% keys.ons(cc).comparisons_per_effector(ce).order_onset={'Hold', 0, 0.5};
% keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IS_IN/255; keys.colors.NH_CS_IN/255];
% keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='MemL';
% keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Hold';
% 
% cc=cc+1;% 1 ungrouped baseline subtraction
% ce=0;
% keys.ons(cc).comparisons_title       = 'Msac choice tuning';
% keys.ons(cc).group_parameter         = 'ungrouped';
% keys.ons(cc).conditions_to_plot      = {'Msac'};
% 
% keys.ons(cc).choices=[0 1]; %for cell exclusion
% ce=ce+1;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
% keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1];
% keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[-1];
% keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).choice{2}=1;
% keys.ons(cc).comparisons_per_effector(ce).order_onset={'Cue', 0, 0.28};
% keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IN/255; keys.colors.NH_CH/255];
% keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
% keys.ons(cc).comparisons_per_effector(ce).title='Ipsilateral';
% ce=ce+1;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
% keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[1];
% keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
% keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).choice{2}=1;
% keys.ons(cc).comparisons_per_effector(ce).order_onset={'Cue', 0, 0.28};
% keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IN/255; keys.colors.NH_CH/255];
% keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
% keys.ons(cc).comparisons_per_effector(ce).title='Contralateral';
% 
% 
% 
% cc=cc+1;% 1 ungrouped baseline subtraction
% ce=0;
% keys.ons(cc).comparisons_title       = 'Vsac instructed space tuning';
% keys.ons(cc).group_parameter         = 'ungrouped';
% keys.ons(cc).conditions_to_plot      = {'Vsac'};
% 
% keys.ons(cc).choices=0; %for cell exclusion
% ce=ce+1;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
% keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1];
% keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
% keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).choice{2}=0;
% keys.ons(cc).comparisons_per_effector(ce).order_onset={'Visual onset', 0, 0.28};
% keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IS_IN/255; keys.colors.NH_CS_IN/255];
% keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
% keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Visual onset';
% ce=ce+1;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
% keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1];
% keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
% keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).choice{2}=0;
% keys.ons(cc).comparisons_per_effector(ce).order_onset={'Saccade', -0.2, 0.2};
% keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IS_IN/255; keys.colors.NH_CS_IN/255];
% keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
% keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Saccade';
% 
% cc=cc+1;% 1 ungrouped baseline subtraction
% ce=0;
% keys.ons(cc).comparisons_title       = 'Vsac choice tuning';
% keys.ons(cc).group_parameter         = 'ungrouped';
% keys.ons(cc).conditions_to_plot      = {'Vsac'};
% keys.ons(cc).choices=[0 1]; %for cell exclusion
% ce=ce+1;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
% keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1];
% keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[-1];
% keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).choice{2}=1;
% keys.ons(cc).comparisons_per_effector(ce).order_onset={'Visual onset', 0, 0.28};
% keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IN/255; keys.colors.NH_CH/255];
% keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
% keys.ons(cc).comparisons_per_effector(ce).title='Ipsilateral';
% ce=ce+1;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
% keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[1];
% keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
% keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).choice{2}=1;
% keys.ons(cc).comparisons_per_effector(ce).order_onset={'Visual onset', 0, 0.28};
% keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IN/255; keys.colors.NH_CH/255];
% keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
% keys.ons(cc).comparisons_per_effector(ce).title='Contralateral';
% 
% cc=cc+1;% 1 ungrouped baseline subtraction
% ce=0;
% keys.ons(cc).comparisons_title       = 'Vsac epoch space tuning';
% keys.ons(cc).group_parameter         = 'ungrouped';
% keys.ons(cc).conditions_to_plot      = {'Vsac'};
% keys.ons(cc).choices=[0]; %for cell exclusion
% ce=ce+1;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
% keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1 1];
% keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[-1 1];
% keys.ons(cc).comparisons_per_effector(ce).choice{1}=[0];
% keys.ons(cc).comparisons_per_effector(ce).choice{2}=[0];
% keys.ons(cc).comparisons_per_effector(ce).order_onset={'Visual onset', 0, 0.28};
% keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.EP_EN/255; keys.colors.EP_SU/255];
% keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='Fhol';
% keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Visual onset';
% ce=ce+1;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
% keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1 1];
% keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[-1 1];
% keys.ons(cc).comparisons_per_effector(ce).choice{1}=[0];
% keys.ons(cc).comparisons_per_effector(ce).choice{2}=[0];
% keys.ons(cc).comparisons_per_effector(ce).order_onset={'Saccade', -0.2, 0.2};
% keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.EP_EN/255; keys.colors.EP_SU/255];
% keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='Fhol';
% keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Saccade';


% %% population PSTH settings
% cc=0;
% cc=cc+1;% 1 ungrouped baseline subtraction
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'ungrouped';
% keys.pop(cc).conditions_to_plot      = {'Vsac'};
% keys.pop(cc).FR_subtract_baseline    = 1;
% keys.pop(cc).ylim                    = [-0.5 3.5];
% 
% cc=cc+1;% 1 ungrouped baseline subtraction
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'ungrouped';
% keys.pop(cc).conditions_to_plot      = {'Msac'};
% keys.pop(cc).FR_subtract_baseline    = 1;
% keys.pop(cc).ylim                    = [-0.5 3.5];
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
% cc=cc+1;% 1 ungrouped baseline subtraction
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'ungrouped';
% keys.pop(cc).conditions_to_plot      = {'Msac'};
% keys.pop(cc).FR_subtract_baseline    = 1;
% keys.pop(cc).ylim                    = [-0.5 3.5];
% cc=cc+1;% 1 ungrouped baseline subtraction
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'ungrouped';
% keys.pop(cc).conditions_to_plot      = {'Vsac'};
% keys.pop(cc).FR_subtract_baseline    = 1;
% keys.pop(cc).ylim                    = [-0.5 3.5];
% cc=cc+1;% 1 ungrouped baseline subtraction
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).group_parameter         = 'ungrouped';
% keys.pop(cc).conditions_to_plot      = {'Msac'};
% keys.pop(cc).FR_subtract_baseline    = 1;
% keys.pop(cc).ylim                    = [-0.5 3.5];
% 
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
% keys.pop(cc).plot_RF                 	= 1;
% keys.pop(cc).y_lim                      = [0 1.6];
