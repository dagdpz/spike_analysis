keys.project_versions={''};
keys.project_version='20180222';
keys.filelist_formatted={};

keys.plot.single_cells                  =0;         % perform single cell plotting

%% to check carefully
keys.position_and_plotting_arrangements         ={'options'};

%% computation settings
keys.cal.effectors                  =[0];
keys.cal.reach_hand                 =[0];
keys.cal.types                      =[2,3];
keys.cal.min_trials_per_condition       =4;                        % about to be used !!


keys.cal.stablity                       =[1];                     % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the table
keys.cal.single_rating                  =[1,2,3];                   % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the tabl
keys.cal.SNR_rating                     =[1,2,3];                   % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the tabl


%% batching
keys.batching.combine_monkeys       =1;
keys.batching.targets               ={'dPulv'};
keys.batching.monkeys               ={'Curius','Linus'};
keys.Curius.date                    ='[20150515 20150814]';
keys.Linus.date                     ='[20150508 20151030]';
keys.Linus.color      =[0 0 255]/255;
keys.Curius.color     =[255 0 0]/255;

%% for cell counts in general, but could also be used for PSTH!?
keys.tt.epoch_criterion             ='none';
keys.tt.position_criterion          ='none';
keys.tt.space_criterion             ='interaction or space only'; %?????
keys.tt.space_criterion             ='none';

%% epochs
keys.EPOCHS_PER_TYPE{2}={...
    'INI',      2,	-0.4,	-0.1,   'INI';...
    'Fhol',     4,	-0.3,	0,      'INI';...
    'Cue',      4,	0.07,	0.17,   'INI';...
    'CueG',     4,	0.04,	0.14,   'INI';...
    'PreS',     60,	-0.1,	-0.01,  'INI';...
    'PeriS',    60,	-0.01,  0.05,   'INI';...
    'Tacq',     5,	0,      0.1,    'INI';...
    'TholG',     5,	0,      0.15,    'INI';...
    'Thol',     20,	-0.3,	0,      'INI';...
    };

keys.ANOVAS_PER_TYPE(2).epoch={'INI' 'Facq';...
    'INI' 'Fhol';...
    'Fhol' 'Cue';...
    'Fhol' 'CueG';...
    'Fhol' 'PreS';...
    'Fhol' 'PeriS';...
    'Fhol' 'Tacq';...
    'Fhol' 'TholG';...
    'Fhol' 'Thol'};

keys.WINDOWS_PER_TYPE{2}={...
    'Visual onset', 4,	-0.8,   0.17;...
    'Saccade',      60,	-0.01,  0.22;...
    'Go',           4,	-0.2,   0.3;...
    'T hold',       20,	-0.3,   0.1;...
    };
% keys.WINDOWS_PER_TYPE{2}={...
%     'Saccade',      60,	-0.8,  0.6;...
%     };

keys.ANOVAS_PER_TYPE(2).spaceLR            ={'Fhol','Cue','CueG','PreS','PeriS','Tacq','Thol','TholG'}';
keys.ANOVAS_PER_TYPE(2).positions          ={'Fhol','Cue','CueG','PreS','PeriS','Tacq','Thol','TholG'}';
keys.ANOVAS_PER_TYPE(2).hands              ={'Fhol','Cue','CueG','PreS','PeriS','Tacq','Thol','TholG'}';
keys.ANOVAS_PER_TYPE(2).SxH                ={'Fhol','Cue','CueG','PreS','PeriS','Tacq','Thol','TholG'}';
keys.ANOVAS_PER_TYPE(2).main              ={'INI','Fhol','Cue','PreS','PeriS','Tacq','Thol','TholG'}';

keys.EPOCHS_PER_TYPE{3}={...
    'INI',      2,	-0.4,	-0.1,   'INI';...
    'Fhol',     6,	-0.3,	0,      'INI';...
    'Cue',      6,	0.07,   0.17,   'INI';...
    'CueG',     6,	0.04,   0.14,   'INI';...
    'MemE',     7,	0.07,   0.17,   'INI';... 
    'MemL',     9,	-0.3, 	-0,     'INI';... 
    'PreS',     60,	-0.1, 	-0.01,  'INI';...
    'Pre2',     60,	-0.1, 	-0.01,  'INI';...
    'PreG',     9,	0.04,   0.14,   'INI';...
    'PeriS',	60,	-0.01,  0.05,   'INI';...
    'Peri2',	60,	-0.01,  0.05,   'INI';...
    'TIhol',	10,	0,      0.15,   'INI';...
    'TIholG',	10,	0,      0.15,   'INI';...
    'Tons',     4,	0.02,   0.12,    'INI';...
    'Thol',     5,	0.2,    0.5,    'INI';...
    };

keys.WINDOWS_PER_TYPE{3}={...
    'Cue',          6,	-0.8,   0.78;...
    'Saccade',      60,	-0.7,   0.8;...
    'Go',           9,	-0.2,   0.3;...
    'Target',       4,	-0.2,   0.7;...
    };

keys.ANOVAS_PER_TYPE(3).epoch={'INI' 'Fhol';...
    'Fhol' 'Cue';...
    'Fhol' 'CueG';...
    'Fhol' 'MemE';...
    'Fhol' 'MemL';...
    'MemL' 'PreS';...
    'Fhol' 'Pre2';...
    'Fhol' 'PreG';...
    'MemL' 'PeriS';...
    'Fhol' 'Peri2';...
    'MemL' 'TIhol';...
    'Fhol' 'TIholG';...
    'MemL' 'Tons';...
    'MemL' 'Thol';...
    };

keys.ANOVAS_PER_TYPE(3).spaceLR            ={'Fhol','Cue','CueG','MemE','MemL','PreS','Pre2','PreG','PeriS','Peri2','TIhol','TIholG','Tons','Thol'}';
keys.ANOVAS_PER_TYPE(3).positions          ={'Fhol','Cue','CueG','MemE','MemL','PreS','Pre2','PreG','PeriS','Peri2','TIhol','TIholG','Tons','Thol'}';
keys.ANOVAS_PER_TYPE(3).hands              ={'Fhol','Cue','CueG','MemE','MemL','PreS','Pre2','PreG','PeriS','Peri2','TIhol','TIholG','Tons','Thol'}';
keys.ANOVAS_PER_TYPE(3).SxH                ={'Fhol','Cue','CueG','MemE','MemL','PreS','Pre2','PreG','PeriS','Peri2','TIhol','TIholG','Tons','Thol'}';
keys.ANOVAS_PER_TYPE(3).main              ={'INI','Fhol','Cue','MemE','MemL','PreS','PeriS','TIhol','Tons','TIhol','TIholG','Thol'}'; %% 'TIhol' wasnt there??


%% cell count settings
cc=0;


cc=cc+1;
keys.ccs(cc).IC_to_plot             ='ch';
keys.ccs(cc).tt.choices             =[0,1];
keys.ccs(cc).tt.hands               =0;
keys.ccs(cc).factor                 ='space_position';
keys.ccs(cc).conditions_to_plot     ={'Msac'};
keys.ccs(cc).plot_type              ='per_epoch';
keys.ccs(cc).epochs.Msac            ={'INI', 'Fhol','Cue','MemE','MemL','PreS','PeriS','TIhol','Tons','Thol'}';

cc=cc+1;
keys.ccs(cc).tt.choices             =0;
keys.ccs(cc).tt.hands             =0;
keys.ccs(cc).factor                 ='space_position';
keys.ccs(cc).conditions_to_plot      ={'Msac'};
keys.ccs(cc).plot_type              ='per_epoch';
keys.ccs(cc).epochs.Msac        ={'INI', 'Fhol','Cue','MemE','MemL','PreS','PeriS','TIhol','Tons','Thol'}';

cc=cc+1;
keys.ccs(cc).tt.choices             =0;
keys.ccs(cc).tt.hands             =0;
keys.ccs(cc).factor                 ='position_space';
keys.ccs(cc).conditions_to_plot      ={'Msac'};
keys.ccs(cc).plot_type              ='per_epoch';
keys.ccs(cc).epochs.Msac        ={'INI', 'Fhol','Cue','MemE','MemL','PreS','PeriS','TIhol','Tons','Thol'}';

cc=cc+1;
keys.ccs(cc).tt.choices             =0;
keys.ccs(cc).tt.hands             =0;
keys.ccs(cc).factor                 ='epoch';
keys.ccs(cc).conditions_to_plot      ={'Msac'};
keys.ccs(cc).plot_type              ='visuomotor';
keys.ccs(cc).epochs.Msac        ={'INI', 'Fhol','Cue','MemE','MemL','PreS','PeriS','TIhol','Tons','Thol'}';

cc=cc+1;
keys.ccs(cc).tt.choices             =0;
keys.ccs(cc).tt.hands             =0;
keys.ccs(cc).factor                 ='epoch';
keys.ccs(cc).conditions_to_plot      ={'Msac'};
keys.ccs(cc).plot_type              ='per_epoch';
keys.ccs(cc).epochs.Msac        ={'INI', 'Fhol','Cue','MemE','MemL','PreS','PeriS','TIhol','Tons','Thol'}';

%% tuning onset settings
cc=0;

cc=cc+1;% 'Msac epoch tuning';
ce=0;
keys.ons(cc).comparisons_title       = 'Msac epoch tuning';
keys.ons(cc).group_parameter         = 'ungrouped';
keys.ons(cc).conditions_to_plot      = {'Msac'};
keys.ons(cc).tt.choices=0; %for cell exclusion
keys.ons(cc).tt.hands=0; %for cell exclusion
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1 1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[-1 1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=[0];
keys.ons(cc).comparisons_per_effector(ce).choice{2}=[0];
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Cue', 0, 0.3};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.EP_EN/255; keys.colors.EP_SU/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='Fhol';
keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Cue';
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1 1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[-1 1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=[0];
keys.ons(cc).comparisons_per_effector(ce).choice{2}=[0];
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Saccade', -0.2, 0.7};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.EP_EN/255; keys.colors.EP_SU/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='MemL';
keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Saccade';
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1 1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[-1 1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=[0];
keys.ons(cc).comparisons_per_effector(ce).choice{2}=[0];
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Target', -0.2, 0.7};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.EP_EN/255; keys.colors.EP_SU/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='MemL';
keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Target';


cc=cc+1;% 'Msac space tuning';
ce=0;
keys.ons(cc).comparisons_title       = 'Msac instructed space tuning';
keys.ons(cc).group_parameter         = 'ungrouped';
keys.ons(cc).conditions_to_plot      = {'Msac'};
keys.ons(cc).tt.choices=0; %for cell exclusion
keys.ons(cc).tt.hands=0; %for cell exclusion
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
keys.ons(cc).comparisons_per_effector(ce).choice{2}=0;
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Cue', 0, 0.3};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IS_IN/255; keys.colors.NH_CS_IN/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Cue';
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
keys.ons(cc).comparisons_per_effector(ce).choice{2}=0;
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Saccade', -0.2, 0.7};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IS_IN/255; keys.colors.NH_CS_IN/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='MemL';
keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Saccade';
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
keys.ons(cc).comparisons_per_effector(ce).choice{2}=0;
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Target', -0.2, 0.7};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IS_IN/255; keys.colors.NH_CS_IN/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='MemL';
keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Target';

% 
% 
% cc=cc+1;% 'Msac choice tuning - choice vs instructed';
% ce=0;
% keys.ons(cc).comparisons_title       = 'Msac choice tuning';
% keys.ons(cc).group_parameter         = 'ungrouped';
% keys.ons(cc).conditions_to_plot      = {'Msac'};
% keys.ons(cc).tt.hands=0; %for cell exclusion
% keys.ons(cc).tt.choices=[0 1]; %for cell exclusion
% ce=ce+1;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
% keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1];
% keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[-1];
% keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).choice{2}=1;
% keys.ons(cc).comparisons_per_effector(ce).order_onset={'Cue', 0, 0.78};
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
% keys.ons(cc).comparisons_per_effector(ce).order_onset={'Cue', 0, 0.78};
% keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IN/255; keys.colors.NH_CH/255];
% keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
% keys.ons(cc).comparisons_per_effector(ce).title='Contralateral';
% 
% cc=cc+1;% msac choice preference
% ce=0;
% keys.ons(cc).comparisons_title       = 'Msac choice preference';
% keys.ons(cc).group_parameter         = 'ungrouped';
% keys.ons(cc).conditions_to_plot      = {'Msac'};
% keys.ons(cc).tt.hands=0; %for cell exclusion
% keys.ons(cc).tt.choices=[1]; %for cell exclusion
% ce=ce+1;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
% keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
% keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1];
% keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
% keys.ons(cc).comparisons_per_effector(ce).choice{1}=1;
% keys.ons(cc).comparisons_per_effector(ce).choice{2}=1;
% keys.ons(cc).comparisons_per_effector(ce).order_onset={'Cue', 0, 0.78};
% keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.IS_CH/255; keys.colors.CS_CH/255];
% keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
% keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Cue';


%% population PSTH settings
cc=0;
 
cc=cc+1;% enhancement/suppression, POS-NEG
keys.pop(cc).epoch_PF                = 'Thol';
keys.pop(cc).epoch_RF                = 'Thol';
keys.pop(cc).epoch_BL                = 'Fhol';
keys.pop(cc).epoch_GB                = 'none';
keys.pop(cc).epoch_for_normalization = 'Thol'; % for percent change, use this isntead!
keys.pop(cc).normalization           = 'percent_change';
keys.pop(cc).group_parameter         = 'in_AH_Thol_position_Msac_opt';
keys.pop(cc).conditions_to_plot      = {'Vsac'};
keys.pop(cc).position_and_plotting_arrangements      = {'options'};
keys.pop(cc).FR_subtract_baseline    = 0;
keys.pop(cc).ylim                    = [];
keys.pop(cc).plot_per_position       = 1;
keys.pop(cc).tt.choices              = 0;
keys.pop(cc).tt.hands                = 0;
keys.pop(cc).plot_RF                 = 0;

cc=cc+1;% enhancement/suppression, POS-NEG
keys.pop(cc).epoch_PF                = 'Thol';
keys.pop(cc).epoch_RF                = 'Thol';
keys.pop(cc).epoch_BL                = 'Fhol';
keys.pop(cc).epoch_GB                = 'none';
keys.pop(cc).epoch_for_normalization = 'Thol'; % for percent change, use this isntead!
keys.pop(cc).normalization           = 'percent_change';
keys.pop(cc).group_parameter         = 'in_AH_Thol_position_Msac_opt';
keys.pop(cc).group_parameter         = 'ungrouped';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).position_and_plotting_arrangements      = {'options'};
keys.pop(cc).FR_subtract_baseline    = 0;
keys.pop(cc).ylim                    = [];
keys.pop(cc).plot_per_position       = 1;
keys.pop(cc).tt.choices              = 0;
keys.pop(cc).tt.hands                = 0;
keys.pop(cc).plot_RF                 = 0;

% cc=cc+1;% instructed, divisive normalization in Fhol
% keys.pop(cc).group_parameter         = 'ungrouped';
% keys.pop(cc).group_excluded          ={};
% keys.pop(cc).conditions_to_plot      = {'Msac'};
% keys.pop(cc).tt.choices                 = [0];
% keys.pop(cc).tt.hands                   = [0];
% keys.pop(cc).normalization           = 'by_effector';
% keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% keys.pop(cc).epoch_PF                = 'Cue';
% keys.pop(cc).FR_subtract_baseline    = 0;
% 
% 
% cc=cc+1;% instructed, subtracting INI
% keys.pop(cc).group_parameter         = 'ungrouped';
% keys.pop(cc).group_excluded          ={};
% keys.pop(cc).conditions_to_plot      = {'Msac'};
% keys.pop(cc).tt.choices                 = [0];
% keys.pop(cc).tt.hands                   = [0];
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% keys.pop(cc).epoch_PF                = 'Cue';
% keys.pop(cc).epoch_BL                = 'INI';
% keys.pop(cc).FR_subtract_baseline    = 1;
% 
% %% Response fields
% % 
% cc=cc+1;% cue response fields 
% keys.pop(cc).normalization              = 'by_effector';
% keys.pop(cc).group_parameter            = 'in_AH_Cue_position_Msac_opt';
% keys.pop(cc).group_excluded             = {'false'};
% keys.pop(cc).epoch_PF                  	= 'Cue';
% keys.pop(cc).epoch_RF                  	= 'Cue';
% keys.pop(cc).epoch_BL                  	= 'Fhol';
% keys.pop(cc).epoch_for_normalization  	= 'Fhol';
% keys.pop(cc).conditions_to_plot       	= {'Msac'};
% keys.pop(cc).tt.choices                    = [0];
% keys.pop(cc).tt.hands                      = [0];
% keys.pop(cc).plot_RF                 	= 1;
% 
% cc=cc+1;% cue preference for enhanced only
% keys.pop(cc).tt.selection             ={'in_AH_Cue_epoch_Msac_opt','en'};
% keys.pop(cc).group_parameter         = 'ungrouped';
% keys.pop(cc).normalization              = 'by_effector';
% keys.pop(cc).epoch_PF                  	= 'Cue';
% keys.pop(cc).epoch_RF                  	= 'Cue';
% keys.pop(cc).epoch_BL                  	= 'Fhol';
% keys.pop(cc).epoch_for_normalization  	= 'Fhol';
% keys.pop(cc).conditions_to_plot       	= {'Msac'};
% keys.pop(cc).tt.choices                    = [0];
% keys.pop(cc).tt.hands                      = [0];
% keys.pop(cc).plot_RF                 	= 1;
% 
% cc=cc+1;% PostS response fields 
% keys.pop(cc).normalization              = 'by_effector';
% keys.pop(cc).group_parameter            = 'in_AH_TIhol_position_Msac_opt';
% keys.pop(cc).group_excluded             = {'false'};
% keys.pop(cc).epoch_PF                  	= 'TIhol';
% keys.pop(cc).epoch_RF                  	= 'TIhol';
% keys.pop(cc).epoch_BL                  	= 'MemL';
% keys.pop(cc).epoch_for_normalization  	= 'MemL';
% keys.pop(cc).conditions_to_plot       	= {'Msac'};
% keys.pop(cc).tt.choices                    = [0];
% keys.pop(cc).tt.hands                      = [0];
% keys.pop(cc).plot_RF                 	= 1;
% 
% %% categories
% 
% % cc=cc+1;% visual response?
% % keys.pop(cc).group_parameter         = 'visual_Msac_opt';
% % keys.pop(cc).group_excluded          ={'0'};
% % keys.pop(cc).conditions_to_plot      = {'Msac'};
% % keys.pop(cc).tt.choices                 = [0];
% % keys.pop(cc).tt.hands                   = [0];
% % keys.pop(cc).normalization           = 'by_effector';
% % keys.pop(cc).epoch_PF                = 'Cue';
% % keys.pop(cc).epoch_RF                = 'Cue';
% % keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% % keys.pop(cc).FR_subtract_baseline    = 0;
% % keys.pop(cc).plot_RF                 = 1;
% % 
% % cc=cc+1;% motor output only?
% % keys.pop(cc).group_parameter         = 'motor_Msac_opt';
% % keys.pop(cc).group_excluded          ={'0'};
% % keys.pop(cc).conditions_to_plot      = {'Msac'};
% % keys.pop(cc).tt.choices                 = [0];
% % keys.pop(cc).tt.hands                   = [0];
% % keys.pop(cc).normalization           = 'by_effector';
% % keys.pop(cc).epoch_PF                = 'TIhol';
% % keys.pop(cc).epoch_RF                = 'TIhol';
% % keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% % keys.pop(cc).FR_subtract_baseline    = 0;
% % keys.pop(cc).plot_RF                 = 1;
% % 
% % cc=cc+1;% visuomotor cells???
% % keys.pop(cc).group_parameter         = 'visuomotor_Msac_opt';
% % keys.pop(cc).group_excluded          ={'0'};
% % keys.pop(cc).conditions_to_plot      = {'Msac'};
% % keys.pop(cc).tt.choices                 = [0];
% % keys.pop(cc).tt.hands                   = [0];
% % keys.pop(cc).normalization           = 'by_effector';
% % keys.pop(cc).epoch_PF                = 'Cue';
% % keys.pop(cc).epoch_RF                = 'Cue';
% % keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% % keys.pop(cc).FR_subtract_baseline    = 0;
% % keys.pop(cc).plot_RF                 = 1;
% % 
% % cc=cc+1;% visuomotor cells???
% % keys.pop(cc).group_parameter         = 'visuomotor_Msac_opt';
% % keys.pop(cc).group_excluded          ={'0'};
% % keys.pop(cc).conditions_to_plot      = {'Msac'};
% % keys.pop(cc).tt.choices                 = [0];
% % keys.pop(cc).tt.hands                   = [0];
% % keys.pop(cc).normalization           = 'by_effector';
% % keys.pop(cc).epoch_PF                = 'TIhol';
% % keys.pop(cc).epoch_RF                = 'TIhol';
% % keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% % keys.pop(cc).FR_subtract_baseline    = 0;
% % keys.pop(cc).plot_RF                 = 1;
% 
% 
% cc=cc+1;% visual response?
% keys.pop(cc).tt.selection             ={'visual_Msac_opt','1'};
% keys.pop(cc).group_parameter         = 'in_AH_Cue_epoch_Msac_opt';
% keys.pop(cc).group_excluded          = {'','-','bi'}; 
% keys.pop(cc).conditions_to_plot      = {'Msac'};
% keys.pop(cc).tt.choices              = [0];
% keys.pop(cc).tt.hands                = [0];
% keys.pop(cc).normalization           = 'by_effector';
% keys.pop(cc).epoch_PF                = 'Cue';
% keys.pop(cc).epoch_RF                = 'Cue';
% keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% keys.pop(cc).FR_subtract_baseline    = 0;
% keys.pop(cc).plot_RF                 = 0;
% keys.pop(cc).y_lim                   = [0 4];
% 
% cc=cc+1;% visual visuomotor response?
% keys.pop(cc).tt.selection             ={'visuomotor_Msac_opt','1'};
% keys.pop(cc).group_parameter         = 'in_AH_Cue_epoch_Msac_opt';
% keys.pop(cc).group_excluded          = {'','-','bi'}; 
% keys.pop(cc).conditions_to_plot      = {'Msac'};
% keys.pop(cc).tt.choices              = [0];
% keys.pop(cc).tt.hands                = [0];
% keys.pop(cc).normalization           = 'by_effector';
% keys.pop(cc).epoch_PF                = 'Cue';
% keys.pop(cc).epoch_RF                = 'Cue';
% keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% keys.pop(cc).FR_subtract_baseline    = 0;
% keys.pop(cc).plot_RF                 = 0;
% keys.pop(cc).y_lim                   = [0 4];
% 
% cc=cc+1;% motor visuomotor response?
% keys.pop(cc).tt.selection             ={'visuomotor_Msac_opt','1'};
% keys.pop(cc).group_parameter         = 'in_AH_TIhol_epoch_Msac_opt';
% keys.pop(cc).group_excluded          = {'','-','bi'}; 
% keys.pop(cc).conditions_to_plot      = {'Msac'};
% keys.pop(cc).tt.choices              = [0];
% keys.pop(cc).tt.hands                = [0];
% keys.pop(cc).normalization           = 'by_effector';
% keys.pop(cc).epoch_PF                = 'TIhol';
% keys.pop(cc).epoch_RF                = 'TIhol';
% keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% keys.pop(cc).FR_subtract_baseline    = 0;
% keys.pop(cc).plot_RF                 = 0;
% keys.pop(cc).y_lim                   = [0 4];
% 
% cc=cc+1;% motor response?
% keys.pop(cc).tt.selection             ={'motor_Msac_opt','1'};
% keys.pop(cc).group_parameter         = 'in_AH_TIhol_epoch_Msac_opt';
% keys.pop(cc).group_excluded          = {'','-','bi'}; 
% keys.pop(cc).conditions_to_plot      = {'Msac'};
% keys.pop(cc).tt.choices              = [0];
% keys.pop(cc).tt.hands                = [0];
% keys.pop(cc).normalization           = 'by_effector';
% keys.pop(cc).epoch_PF                = 'TIhol';
% keys.pop(cc).epoch_RF                = 'TIhol';
% keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% keys.pop(cc).FR_subtract_baseline    = 0;
% keys.pop(cc).plot_RF                 = 0;
% keys.pop(cc).y_lim                   = [0 4];
% % 
% % %% categories premotor2
% % cc=cc+1;% visual response?
% % keys.pop(cc).group_parameter         = 'visual_pre2_Msac_opt';
% % keys.pop(cc).group_excluded          ={'0'};
% % keys.pop(cc).conditions_to_plot      = {'Msac'};
% % keys.pop(cc).tt.choices                 = [0];
% % keys.pop(cc).tt.hands                   = [0];
% % keys.pop(cc).normalization           = 'by_effector';
% % keys.pop(cc).epoch_PF                = 'Cue';
% % keys.pop(cc).epoch_RF                = 'Cue';
% % keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% % keys.pop(cc).FR_subtract_baseline    = 0;
% % keys.pop(cc).plot_RF                 = 1;
% % 
% % cc=cc+1;% motor output only?
% % keys.pop(cc).group_parameter         = 'motor_pre2_Msac_opt';
% % keys.pop(cc).group_excluded          ={'0'};
% % keys.pop(cc).conditions_to_plot      = {'Msac'};
% % keys.pop(cc).tt.choices                 = [0];
% % keys.pop(cc).tt.hands                   = [0];
% % keys.pop(cc).normalization           = 'by_effector';
% % keys.pop(cc).epoch_PF                = 'Pre2';
% % keys.pop(cc).epoch_RF                = 'Pre2';
% % keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% % keys.pop(cc).FR_subtract_baseline    = 0;
% % keys.pop(cc).plot_RF                 = 1;
% % 
% % cc=cc+1;% visuomotor cells???
% % keys.pop(cc).group_parameter         = 'visuomotor_pre2_Msac_opt';
% % keys.pop(cc).group_excluded          ={'0'};
% % keys.pop(cc).conditions_to_plot      = {'Msac'};
% % keys.pop(cc).tt.choices                 = [0];
% % keys.pop(cc).tt.hands                   = [0];
% % keys.pop(cc).normalization           = 'by_effector';
% % keys.pop(cc).epoch_PF                = 'Cue';
% % keys.pop(cc).epoch_RF                = 'Cue';
% % keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% % keys.pop(cc).FR_subtract_baseline    = 0;
% % keys.pop(cc).plot_RF                 = 1;
% % 
% % cc=cc+1;% visuomotor cells???
% % keys.pop(cc).group_parameter         = 'visuomotor_pre2_Msac_opt';
% % keys.pop(cc).group_excluded          ={'0'};
% % keys.pop(cc).conditions_to_plot      = {'Msac'};
% % keys.pop(cc).tt.choices                 = [0];
% % keys.pop(cc).tt.hands                   = [0];
% % keys.pop(cc).normalization           = 'by_effector';
% % keys.pop(cc).epoch_PF                = 'Pre2';
% % keys.pop(cc).epoch_RF                = 'Pre2';
% % keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% % keys.pop(cc).FR_subtract_baseline    = 0;
% % keys.pop(cc).plot_RF                 = 1;
% 
% %% categories premotor2
% 
% cc=cc+1;% visual response?
% keys.pop(cc).tt.selection             ={'visual_pre2_Msac_opt','1'};
% keys.pop(cc).group_parameter         = 'in_AH_Cue_epoch_Msac_opt';
% keys.pop(cc).group_excluded          = {'','-','bi'}; 
% keys.pop(cc).conditions_to_plot      = {'Msac'};
% keys.pop(cc).tt.choices              = [0];
% keys.pop(cc).tt.hands                = [0];
% keys.pop(cc).normalization           = 'by_effector';
% keys.pop(cc).epoch_PF                = 'Cue';
% keys.pop(cc).epoch_RF                = 'Cue';
% keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% keys.pop(cc).FR_subtract_baseline    = 0;
% keys.pop(cc).plot_RF                 = 0;
% keys.pop(cc).y_lim                   = [0 4];
% 
% cc=cc+1;% visual visuomotor response?
% keys.pop(cc).tt.selection             ={'visuomotor_pre2_Msac_opt','1'};
% keys.pop(cc).group_parameter         = 'in_AH_Cue_epoch_Msac_opt';
% keys.pop(cc).group_excluded          = {'','-','bi'}; 
% keys.pop(cc).conditions_to_plot      = {'Msac'};
% keys.pop(cc).tt.choices              = [0];
% keys.pop(cc).tt.hands                = [0];
% keys.pop(cc).normalization           = 'by_effector';
% keys.pop(cc).epoch_PF                = 'Cue';
% keys.pop(cc).epoch_RF                = 'Cue';
% keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% keys.pop(cc).FR_subtract_baseline    = 0;
% keys.pop(cc).plot_RF                 = 0;
% keys.pop(cc).y_lim                   = [0 4];
% 
% cc=cc+1;% motor visuomotor response?
% keys.pop(cc).tt.selection             ={'visuomotor_pre2_Msac_opt','1'};
% keys.pop(cc).group_parameter         = 'in_AH_Pre2_epoch_Msac_opt';
% keys.pop(cc).group_excluded          = {'','-','bi'}; 
% keys.pop(cc).conditions_to_plot      = {'Msac'};
% keys.pop(cc).tt.choices              = [0];
% keys.pop(cc).tt.hands                = [0];
% keys.pop(cc).normalization           = 'by_effector';
% keys.pop(cc).epoch_PF                = 'Pre2';
% keys.pop(cc).epoch_RF                = 'Pre2';
% keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% keys.pop(cc).FR_subtract_baseline    = 0;
% keys.pop(cc).plot_RF                 = 0;
% keys.pop(cc).y_lim                   = [0 4];
% 
% cc=cc+1;% motor response?
% keys.pop(cc).tt.selection             ={'motor_pre2_Msac_opt','1'};
% keys.pop(cc).group_parameter         = 'in_AH_Pre2_epoch_Msac_opt';
% keys.pop(cc).group_excluded          = {'','-','bi'}; 
% keys.pop(cc).conditions_to_plot      = {'Msac'};
% keys.pop(cc).tt.choices              = [0];
% keys.pop(cc).tt.hands                = [0];
% keys.pop(cc).normalization           = 'by_effector';
% keys.pop(cc).epoch_PF                = 'Pre2';
% keys.pop(cc).epoch_RF                = 'Pre2';
% keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% keys.pop(cc).FR_subtract_baseline    = 0;
% keys.pop(cc).plot_RF                 = 0;
% keys.pop(cc).y_lim                   = [0 4];
% 
% %% epoch & space tuning in categories
% % cc=cc+1;% visual cells???
% % keys.pop(cc).tt.combine_tuning_properties  ={'comb_cue','in_AH_Cue_epoch_Msac_opt','in_Cue_spaceLR_Msac_opt'};
% % keys.pop(cc).group_parameter         = 'comb_cue';
% % keys.pop(cc).group_excluded          = {''};
% % keys.pop(cc).conditions_to_plot      = {'Msac'};
% % keys.pop(cc).tt.choices                 = [0];
% % keys.pop(cc).tt.hands                   = [0];
% % keys.pop(cc).normalization           = 'by_effector';
% % keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% % keys.pop(cc).FR_subtract_baseline    = 0;
% % keys.pop(cc).tt.selection                  ={'visual_Msac_opt','1'};
% % cc=cc+1;% visuomotor cells???
% % keys.pop(cc).tt.combine_tuning_properties  ={'comb_cue','in_AH_Cue_epoch_Msac_opt','in_Cue_spaceLR_Msac_opt'};
% % keys.pop(cc).group_parameter         = 'comb_cue';
% % keys.pop(cc).conditions_to_plot      = {'Msac'};
% % keys.pop(cc).tt.choices                 = [0];
% % keys.pop(cc).tt.hands                   = [0];
% % keys.pop(cc).normalization           = 'by_effector';
% % keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% % keys.pop(cc).FR_subtract_baseline    = 0;
% % keys.pop(cc).tt.selection                  ={'visuomotor_Msac_opt','1'};
% % cc=cc+1;% visuomotor cells???
% % keys.pop(cc).tt.combine_tuning_properties  ={'comb_sac','in_AH_TIhol_epoch_Msac_opt','in_TIhol_spaceLR_Msac_opt'};
% % keys.pop(cc).group_parameter         = 'comb_sac';
% % keys.pop(cc).group_excluded          = {''};
% % keys.pop(cc).conditions_to_plot      = {'Msac'};
% % keys.pop(cc).tt.choices                 = [0];
% % keys.pop(cc).tt.hands                   = [0];
% % keys.pop(cc).normalization           = 'by_effector';
% % keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% % keys.pop(cc).FR_subtract_baseline    = 0;
% % keys.pop(cc).tt.selection                  ={'visuomotor_Msac_opt','1'};
% % cc=cc+1;% motor cells???
% % keys.pop(cc).tt.combine_tuning_properties  ={'comb_sac','in_AH_TIhol_epoch_Msac_opt','in_TIhol_spaceLR_Msac_opt'};
% % keys.pop(cc).group_parameter         = 'comb_sac';
% % keys.pop(cc).group_excluded          = {''};
% % keys.pop(cc).conditions_to_plot      = {'Msac'};
% % keys.pop(cc).tt.choices                 = [0];
% % keys.pop(cc).tt.hands                   = [0];
% % keys.pop(cc).normalization           = 'by_effector';
% % keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% % keys.pop(cc).FR_subtract_baseline    = 0;
% % keys.pop(cc).tt.selection                  ={'motor_Msac_opt',1};
% % 
% % 
% % 
% % %% epoch & space tuning in categories pre
% % cc=cc+1;% visual cells???
% % keys.pop(cc).tt.combine_tuning_properties  ={'comb_cue','in_AH_Cue_epoch_Msac_opt','in_Cue_spaceLR_Msac_opt'};
% % keys.pop(cc).group_parameter         = 'comb_cue';
% % keys.pop(cc).group_excluded          = {''};
% % keys.pop(cc).conditions_to_plot      = {'Msac'};
% % keys.pop(cc).tt.choices                 = [0];
% % keys.pop(cc).tt.hands                   = [0];
% % keys.pop(cc).normalization           = 'by_effector';
% % keys.pop(cc).epoch_for_normalization = 'Fhol';
% % keys.pop(cc).FR_subtract_baseline    = 0;
% % keys.pop(cc).tt.selection                  ={'visual_pre2_Msac_opt','1'};
% % cc=cc+1;% visuomotor cells???
% % keys.pop(cc).tt.combine_tuning_properties  ={'comb_cue','in_AH_Cue_epoch_Msac_opt','in_Cue_spaceLR_Msac_opt'};
% % keys.pop(cc).group_parameter         = 'comb_cue';
% % keys.pop(cc).conditions_to_plot      = {'Msac'};
% % keys.pop(cc).tt.choices                 = [0];
% % keys.pop(cc).tt.hands                   = [0];
% % keys.pop(cc).normalization           = 'by_effector';
% % keys.pop(cc).epoch_for_normalization = 'Fhol'; 
% % keys.pop(cc).FR_subtract_baseline    = 0;
% % keys.pop(cc).tt.selection                  ={'visuomotor_pre2_Msac_opt','1'};
% % cc=cc+1;% visuomotor cells???
% % keys.pop(cc).tt.combine_tuning_properties  ={'comb_sac','in_AH_Pre2_epoch_Msac_opt','in_Pre2_spaceLR_Msac_opt'};
% % keys.pop(cc).group_parameter         = 'comb_sac';
% % keys.pop(cc).group_excluded          = {''};
% % keys.pop(cc).conditions_to_plot      = {'Msac'};
% % keys.pop(cc).tt.choices                 = [0];
% % keys.pop(cc).tt.hands                   = [0];
% % keys.pop(cc).normalization           = 'by_effector';
% % keys.pop(cc).epoch_for_normalization = 'Fhol'; 
% % keys.pop(cc).FR_subtract_baseline    = 0;
% % keys.pop(cc).tt.selection                  ={'visuomotor_pre2_Msac_opt','1'};
% % cc=cc+1;% motor cells???
% % keys.pop(cc).tt.combine_tuning_properties  ={'comb_sac','in_AH_Pre2_epoch_Msac_opt','in_Pre2_spaceLR_Msac_opt'};
% % keys.pop(cc).group_parameter         = 'comb_sac';
% % keys.pop(cc).group_excluded          = {''};
% % keys.pop(cc).conditions_to_plot      = {'Msac'};
% % keys.pop(cc).tt.choices                 = [0];
% % keys.pop(cc).tt.hands                   = [0];
% % keys.pop(cc).normalization           = 'by_effector';
% % keys.pop(cc).epoch_for_normalization = 'Fhol'; 
% % keys.pop(cc).FR_subtract_baseline    = 0;
% % keys.pop(cc).tt.selection                  ={'motor_pre2_Msac_opt',1};
% 
% % 
% 
% %% gaze cells (?)
% cc=cc+1;% gaze cells???
% keys.pop(cc).group_parameter         = 'in_Thol_spaceLR_Msac_opt';
% keys.pop(cc).conditions_to_plot      = {'Msac'};
% keys.pop(cc).tt.choices                 = [0];
% keys.pop(cc).tt.hands                   = [0];
% keys.pop(cc).normalization           = 'by_effector';
% keys.pop(cc).epoch_for_normalization = 'Fhol';
% keys.pop(cc).FR_subtract_baseline    = 0;
% 
% 
% %% problematic thingy (Pre- and Peri- saccadic enhance/suppression and space preference)
% % 
% % cc=cc+1;% 3 Ramping in 2 steps
% % keys.pop(cc).epoch_PF                = 'PeriS';
% % keys.pop(cc).group_parameter            = 'Ramp_pre';
% % keys.pop(cc).conditions_to_plot         = {'Msac'};
% % keys.pop(cc).group_excluded             ={'-bi','bi-','--','-su','-en','en-','su-','subi','enbi','bisu','bien'};
% % keys.pop(cc).tt.combine_tuning_properties  ={'Ramp_pre','in_AH_MemL_epoch_Msac_opt','in_AH_PreS_epoch_Msac_opt'};
% % keys.pop(cc).tt.choices                 = [0];
% % keys.pop(cc).tt.hands                   = [0];
% % keys.pop(cc).normalization           = 'by_effector';
% % keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% % keys.pop(cc).FR_subtract_baseline    = 0;
% % 
% cc=cc+1;% TIhol position but not HF preference
% keys.pop(cc).epoch_PF                = 'TIhol';
% keys.pop(cc).epoch_RF                = 'TIhol';
% keys.pop(cc).group_parameter            = 'TIhol_nonHF_but_position';
% keys.pop(cc).conditions_to_plot         = {'Msac'};
% keys.pop(cc).tt.selection                  ={'TIhol_nonHF_but_position','-true'};
% keys.pop(cc).tt.combine_tuning_properties  ={'TIhol_nonHF_but_position','in_TIhol_spaceLR_Msac_opt','in_AH_TIhol_position_Msac_opt'};
% keys.pop(cc).tt.choices                 = [0];
% keys.pop(cc).tt.hands                   = [0];
% keys.pop(cc).normalization           = 'by_effector';
% keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% keys.pop(cc).FR_subtract_baseline    = 0;
% keys.pop(cc).plot_RF                 	= 1;
% % 
% % cc=cc+1;% 3 Ramping in 2 steps
% % keys.pop(cc).epoch_PF                = 'PeriS';
% % keys.pop(cc).group_parameter            = 'Ramp_peri';
% % keys.pop(cc).conditions_to_plot         = {'Msac'};
% % keys.pop(cc).tt.combine_tuning_properties  ={'Ramp_peri','in_AH_MemL_epoch_Msac_opt','in_AH_PeriS_epoch_Msac_opt'};
% % keys.pop(cc).tt.choices                 = [0];
% % keys.pop(cc).tt.hands                   = [0];
% % keys.pop(cc).normalization           = 'by_effector';
% % keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% % keys.pop(cc).FR_subtract_baseline    = 0;
% % keys.pop(cc).tt.selection              = {'Ramp_peri','enen'};
% % 
% % cc=cc+1;% 3 Ramping in 2 steps
% % keys.pop(cc).epoch_PF                = 'PeriS';
% % keys.pop(cc).group_parameter            = 'Ramp_peri';
% % keys.pop(cc).conditions_to_plot         = {'Msac'};
% % %keys.pop(cc).group_excluded             ={'-bi','bi-','--','-su','-en','en-','su-','subi','enbi','bisu','bien'};
% % keys.pop(cc).tt.combine_tuning_properties  ={'Ramp_peri','in_AH_MemL_epoch_Msac_opt','in_AH_PeriS_epoch_Msac_opt'};
% % keys.pop(cc).tt.choices                 = [0];
% % keys.pop(cc).tt.hands                   = [0];
% % keys.pop(cc).normalization           = 'by_effector';
% % keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
% % keys.pop(cc).FR_subtract_baseline    = 0;
% % keys.pop(cc).tt.selection              = {'Ramp_peri','susu'};
% % 
% % cc=cc+1;% Pre-saccadic enhancement/suppression
% % keys.pop(cc).group_parameter         = 'in_AH_PreS_epoch_Msac_opt';
% % keys.pop(cc).epoch_RF                = 'PreS';
% % keys.pop(cc).epoch_PF                = 'PreS';
% % keys.pop(cc).conditions_to_plot      = {'Msac'};
% % keys.pop(cc).FR_subtract_baseline    = 1;
% % keys.pop(cc).group_excluded          ={'','-','bi'};
% % keys.pop(cc).tt.choices                 = [0];
% % keys.pop(cc).tt.hands                   = [0];
% % keys.pop(cc).normalization           = 'by_effector';
% % keys.pop(cc).epoch_for_normalization = 'Fhol';
% % keys.pop(cc).FR_subtract_baseline    = 0;
% % cc=cc+1;% 3 Pre-saccadic enhancement/suppression (relative to Fhol)... to be removed, this is to reconstruct the plots that we had previously
% % keys.pop(cc).group_parameter         = 'in_AH_Pre2_epoch_Msac_opt';
% % keys.pop(cc).epoch_RF                = 'Pre2';
% % keys.pop(cc).epoch_PF                = 'Pre2';
% % keys.pop(cc).conditions_to_plot      = {'Msac'};
% % keys.pop(cc).group_excluded          ={'','-','bi'};
% % keys.pop(cc).tt.choices                 = [0];
% % keys.pop(cc).tt.hands                   = [0];
% % keys.pop(cc).normalization           = 'by_effector';
% % keys.pop(cc).epoch_for_normalization = 'Fhol';
% % keys.pop(cc).FR_subtract_baseline    = 0;
% % 
% cc=cc+1;% Peri-saccadic enhancement/suppression
% keys.pop(cc).group_parameter         = 'in_AH_PeriS_epoch_Msac_opt';
% keys.pop(cc).epoch_RF                = 'PeriS';
% keys.pop(cc).epoch_PF                = 'PeriS';
% keys.pop(cc).conditions_to_plot      = {'Msac'};
% keys.pop(cc).FR_subtract_baseline    = 1;
% keys.pop(cc).group_excluded          ={'','-','bi'};
% keys.pop(cc).tt.choices                 = [0];
% keys.pop(cc).tt.hands                   = [0];
% keys.pop(cc).normalization           = 'by_effector';
% keys.pop(cc).epoch_for_normalization = 'Fhol';
% keys.pop(cc).FR_subtract_baseline    = 0;
% cc=cc+1;% 3 Peri-saccadic enhancement/suppression (relative to Fhol)... to be removed, this is to reconstruct the plots that we had previously
% keys.pop(cc).group_parameter         = 'in_AH_Peri2_epoch_Msac_opt';
% keys.pop(cc).epoch_RF                = 'Peri2';
% keys.pop(cc).epoch_PF                = 'Peri2';
% keys.pop(cc).conditions_to_plot      = {'Msac'};
% keys.pop(cc).group_excluded          ={'','-','bi'};
% keys.pop(cc).tt.choices                 = [0];
% keys.pop(cc).tt.hands                   = [0];
% keys.pop(cc).normalization           = 'by_effector';
% keys.pop(cc).epoch_for_normalization = 'Fhol';
% keys.pop(cc).FR_subtract_baseline    = 0;
% 
% cc=cc+1;% 3 Peri-saccadic enhancement/suppression (relative to Fhol)... to be removed, this is to reconstruct the plots that we had previously
% keys.pop(cc).group_parameter         = 'in_AH_Peri2_epoch_Msac_opt';
% keys.pop(cc).epoch_RF                = 'Peri2';
% keys.pop(cc).epoch_PF                = 'Peri2';
% keys.pop(cc).conditions_to_plot      = {'Msac'};
% keys.pop(cc).group_excluded          ={'','-','bi'};
% keys.pop(cc).tt.choices                 = [0];
% keys.pop(cc).tt.hands                   = [0];
% keys.pop(cc).normalization           = 'none';
% keys.pop(cc).epoch_for_normalization = 'Fhol';
% keys.pop(cc).FR_subtract_baseline    = 1;
% 
% cc=cc+1;% memory space preference
% keys.pop(cc).group_parameter         = 'in_AH_MemL_position_Msac_opt';
% keys.pop(cc).epoch_PF                = 'MemL';
% keys.pop(cc).conditions_to_plot      = {'Msac'};
% keys.pop(cc).tt.choices                 = [0 1];
% keys.pop(cc).tt.hands                   = [0];
% keys.pop(cc).normalization           = 'by_effector';
% keys.pop(cc).epoch_for_normalization = 'Fhol';
% keys.pop(cc).FR_subtract_baseline    = 0;
% cc=cc+1;% memory space preference
% keys.pop(cc).group_parameter         = 'in_MemL_spaceLR_Msac_opt';
% keys.pop(cc).epoch_PF                = 'MemL';
% keys.pop(cc).conditions_to_plot      = {'Msac'};
% keys.pop(cc).tt.choices                 = [0 1];
% keys.pop(cc).tt.hands                   = [0];
% keys.pop(cc).normalization           = 'by_effector';
% keys.pop(cc).epoch_for_normalization = 'Fhol';
% keys.pop(cc).FR_subtract_baseline    = 0;
% % cc=cc+1;% memory enhancement/suppression
% % keys.pop(cc).group_parameter         = 'in_AH_MemL_epoch_Msac_opt';
% % keys.pop(cc).epoch_PF                = 'MemL';
% % keys.pop(cc).conditions_to_plot      = {'Msac'};
% % keys.pop(cc).group_excluded          ={'','-','bi'};
% % keys.pop(cc).tt.choices                 = [0 1];
% % keys.pop(cc).tt.hands                   = [0];
% % keys.pop(cc).normalization           = 'by_effector';
% % keys.pop(cc).epoch_for_normalization = 'Fhol';
% % keys.pop(cc).FR_subtract_baseline    = 0;
% cc=cc+1;% memory space preference
% keys.pop(cc).group_parameter         = 'in_MemL_spaceLR_Msac_opt';
% keys.pop(cc).epoch_PF                = 'MemL';
% keys.pop(cc).conditions_to_plot      = {'Msac'};
% keys.pop(cc).tt.choices                 = [0];
% keys.pop(cc).tt.hands                   = [0];
% keys.pop(cc).normalization           = 'by_effector';
% keys.pop(cc).epoch_for_normalization = 'Fhol';
% keys.pop(cc).FR_subtract_baseline    = 0;
% cc=cc+1;% memory enhancement/suppression
% keys.pop(cc).group_parameter         = 'in_AH_MemL_epoch_Msac_opt';
% keys.pop(cc).epoch_PF                = 'MemL';
% keys.pop(cc).conditions_to_plot      = {'Msac'};
% keys.pop(cc).group_excluded          ={'','-','bi'};
% keys.pop(cc).tt.choices                 = [0];
% keys.pop(cc).tt.hands                   = [0];
% keys.pop(cc).normalization           = 'by_effector';
% keys.pop(cc).epoch_for_normalization = 'Fhol';
% keys.pop(cc).FR_subtract_baseline    = 0;

% %% Scatter keys
cs=0;  
% % to test
% cs=cs+1;
% keys.sct(cs).tt.tasktypes={'Msac','Vsac'};
% keys.sct(cs).X='in_AH_IS_CueG_epoch_DF_Msac_opt';
% keys.sct(cs).Y='in_AH_IS_PreG_epoch_DF_Msac_opt';
% keys.sct(cs).X_sig='in_AH_CueG_epoch_Msac_opt';
% keys.sct(cs).Y_sig='in_AH_PreG_epoch_Msac_opt';
% keys.sct(cs).color_option='ENSU_as_color';
% keys.sct(cs).categories={'visual_Msac_opt';'visuomotor_Msac_opt';'motor_Msac_opt'};

%% visually guided versus memory
% 

%% uNCOMMENT

% supra or sublinear
% visually guided versus memory guided DF CONTRA CUE
cs=cs+1;  
keys.sct(cs).tt.tasktypes={'Msac_opt','Vsac_opt'};
keys.sct(cs).X='in_AH_CS_CueG_epoch_DF_Vsac_opt';
keys.sct(cs).Y='PreCueSum_CS_Msac_opt';
keys.sct(cs).X_sig='in_AH_CS_CueG_epoch_Vsac_opt';
keys.sct(cs).Y_sig='in_AH_CS_CueG_epoch_Vsac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).categories={'visual_Msac_opt';'visuomotor_Msac_opt';'motor_Msac_opt'};

% visually guided versus memory guided DF IPSI CUE
cs=cs+1;  
keys.sct(cs).tt.tasktypes={'Msac_opt','Vsac_opt'};
keys.sct(cs).X='in_AH_IS_CueG_epoch_DF_Vsac_opt';
keys.sct(cs).Y='PreCueSum_IS_Msac_opt';
keys.sct(cs).X_sig='in_AH_IS_CueG_epoch_Vsac_opt';
keys.sct(cs).Y_sig='in_AH_IS_CueG_epoch_Vsac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).categories={'visual_Msac_opt';'visuomotor_Msac_opt';'motor_Msac_opt'};

% visually guided versus memory guided FR CONTRA CUE
cs=cs+1;  
keys.sct(cs).tt.tasktypes={'Msac_opt','Vsac_opt'};
keys.sct(cs).X='in_AH_CS_CueG_epoch_FR_Vsac_opt';
keys.sct(cs).Y='PreCueMean_CS_Msac_opt';
keys.sct(cs).X_sig='in_AH_CS_CueG_epoch_Vsac_opt';
keys.sct(cs).Y_sig='in_AH_CS_CueG_epoch_Vsac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).categories={'visual_Msac_opt';'visuomotor_Msac_opt';'motor_Msac_opt'};

% visually guided versus memory guided FR IPSI CUE
cs=cs+1;  
keys.sct(cs).tt.tasktypes={'Msac_opt','Vsac_opt'};
keys.sct(cs).X='in_AH_IS_CueG_epoch_FR_Vsac_opt';
keys.sct(cs).Y='PreCueMean_IS_Msac_opt';
keys.sct(cs).X_sig='in_AH_IS_CueG_epoch_Vsac_opt';
keys.sct(cs).Y_sig='in_AH_IS_CueG_epoch_Vsac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).categories={'visual_Msac_opt';'visuomotor_Msac_opt';'motor_Msac_opt'};

%%% UNBTIL HERE



% 
% 
% % visually guided versus memory guided FR CONTRA CUE
% cs=cs+1;  
% keys.sct(cs).tt.tasktypes={'Msac_opt','Vsac_opt'};
% keys.sct(cs).X='in_AH_CS_CueG_epoch_DF_Vsac_opt';
% keys.sct(cs).Y='in_AH_CS_CueG_epoch_DF_Msac_opt';
% keys.sct(cs).X_sig='in_AH_CS_CueG_epoch_Vsac_opt';
% keys.sct(cs).Y_sig='in_AH_CS_CueG_epoch_Msac_opt';
% keys.sct(cs).color_option='ENSU_as_color';
% keys.sct(cs).categories={'visual_Msac_opt';'visuomotor_Msac_opt';'motor_Msac_opt'};
% 
% % visually guided versus memory guided FR IPSI CUE
% cs=cs+1;  
% keys.sct(cs).tt.tasktypes={'Msac_opt','Vsac_opt'};
% keys.sct(cs).X='in_AH_IS_CueG_epoch_DF_Vsac_opt';
% keys.sct(cs).Y='in_AH_IS_CueG_epoch_DF_Msac_opt';
% keys.sct(cs).X_sig='in_AH_IS_CueG_epoch_Vsac_opt';
% keys.sct(cs).Y_sig='in_AH_IS_CueG_epoch_Msac_opt';
% keys.sct(cs).color_option='ENSU_as_color';
% keys.sct(cs).categories={'visual_Msac_opt';'visuomotor_Msac_opt';'motor_Msac_opt'};
% 
% % visually guided versus memory guided FR CONTRA PreS
% cs=cs+1;
% keys.sct(cs).tt.tasktypes={'Msac_opt','Vsac_opt'};
% keys.sct(cs).X='in_AH_CS_CueG_epoch_DF_Vsac_opt';
% keys.sct(cs).Y='in_AH_CS_PreG_epoch_DF_Msac_opt';
% keys.sct(cs).X_sig='in_AH_CS_CueG_epoch_Vsac_opt';
% keys.sct(cs).Y_sig='in_AH_CS_PreG_epoch_Msac_opt';
% keys.sct(cs).color_option='ENSU_as_color';
% keys.sct(cs).categories={'visual_Msac_opt';'visuomotor_Msac_opt';'motor_Msac_opt'};
% % visually guided versus memory guided FR IPSI PreS
% cs=cs+1;
% keys.sct(cs).tt.tasktypes={'Msac_opt','Vsac_opt'};
% keys.sct(cs).X='in_AH_IS_CueG_epoch_DF_Vsac_opt';
% keys.sct(cs).Y='in_AH_IS_PreG_epoch_DF_Msac_opt';
% keys.sct(cs).X_sig='in_AH_IS_CueG_epoch_Vsac_opt';
% keys.sct(cs).Y_sig='in_AH_IS_PreG_epoch_Msac_opt';
% keys.sct(cs).color_option='ENSU_as_color';
% keys.sct(cs).categories={'visual_Msac_opt';'visuomotor_Msac_opt';'motor_Msac_opt'};
% 
% % visually guided versus memory guided FR CONTRA Thol
% cs=cs+1;
% keys.sct(cs).tt.tasktypes={'Msac_opt','Vsac_opt'};
% keys.sct(cs).X='in_AH_CS_Thol_epoch_DF_Vsac_opt';
% keys.sct(cs).Y='in_AH_CS_TIhol_epoch_DF_Msac_opt';
% keys.sct(cs).X_sig='in_AH_CS_Thol_epoch_Vsac_opt';
% keys.sct(cs).Y_sig='in_AH_CS_TIhol_epoch_Msac_opt';
% keys.sct(cs).color_option='ENSU_as_color';
% keys.sct(cs).categories={'visual_Msac_opt';'visuomotor_Msac_opt';'motor_Msac_opt'};
% % visually guided versus memory guided FR IPSI Thol
% cs=cs+1;
% keys.sct(cs).tt.tasktypes={'Msac_opt','Vsac_opt'};
% keys.sct(cs).X='in_AH_IS_Thol_epoch_DF_Vsac_opt';
% keys.sct(cs).Y='in_AH_IS_TIhol_epoch_DF_Msac_opt';
% keys.sct(cs).X_sig='in_AH_IS_Thol_epoch_Vsac_opt';
% keys.sct(cs).Y_sig='in_AH_IS_TIhol_epoch_Msac_opt';
% keys.sct(cs).color_option='ENSU_as_color';
% keys.sct(cs).categories={'visual_Msac_opt';'visuomotor_Msac_opt';'motor_Msac_opt'};

%% comparison choice and instructed tuning

% preferred position

% Cue
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_Cue_prefP_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_Cue_prefP_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_Cue_position_Msac_opt';
keys.sct(cs).Y_sig='ch_AH_Cue_position_Msac_opt';
keys.sct(cs).color_option='monkeys_by_marker';
keys.sct(cs).tt.choices=[0,1];

% MemE
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_MemE_prefP_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_MemE_prefP_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_MemE_position_Msac_opt';
keys.sct(cs).Y_sig='ch_AH_MemE_position_Msac_opt';
keys.sct(cs).color_option='monkeys_by_marker';
keys.sct(cs).tt.choices=[0,1];

% MemL
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_MemL_prefP_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_MemL_prefP_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_MemL_position_Msac_opt';
keys.sct(cs).Y_sig='ch_AH_MemL_position_Msac_opt';
keys.sct(cs).color_option='monkeys_by_marker';
keys.sct(cs).tt.choices=[0,1];

% PreS
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_PreS_prefP_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_PreS_prefP_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_PreS_position_Msac_opt';
keys.sct(cs).Y_sig='ch_AH_PreS_position_Msac_opt';
keys.sct(cs).color_option='monkeys_by_marker';
keys.sct(cs).tt.choices=[0,1];

% PeriS
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_PeriS_prefP_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_PeriS_prefP_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_PeriS_position_Msac_opt';
keys.sct(cs).Y_sig='ch_AH_PeriS_position_Msac_opt';
keys.sct(cs).color_option='monkeys_by_marker';
keys.sct(cs).tt.choices=[0,1];

% TIhol
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_TIhol_prefP_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_TIhol_prefP_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_TIhol_position_Msac_opt';
keys.sct(cs).Y_sig='ch_AH_TIhol_position_Msac_opt';
keys.sct(cs).color_option='monkeys_by_marker';
keys.sct(cs).tt.choices=[0,1];

%% index

% Cue
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_Cue_spaceCI_IX_Msac_opt';
keys.sct(cs).Y='ch_Cue_spaceCI_IX_Msac_opt';
keys.sct(cs).X_sig='in_Cue_spaceLR_Msac_opt';
keys.sct(cs).Y_sig='ch_Cue_spaceLR_Msac_opt';
keys.sct(cs).color_option='monkeys_by_marker';
keys.sct(cs).tt.choices=[0,1];

% MemE
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_MemE_spaceCI_IX_Msac_opt';
keys.sct(cs).Y='ch_MemE_spaceCI_IX_Msac_opt';
keys.sct(cs).X_sig='in_MemE_spaceLR_Msac_opt';
keys.sct(cs).Y_sig='ch_MemE_spaceLR_Msac_opt';
keys.sct(cs).color_option='monkeys_by_marker';
keys.sct(cs).tt.choices=[0,1];

% MemL
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_MemL_spaceCI_IX_Msac_opt';
keys.sct(cs).Y='ch_MemL_spaceCI_IX_Msac_opt';
keys.sct(cs).X_sig='in_MemL_spaceLR_Msac_opt';
keys.sct(cs).Y_sig='ch_MemL_spaceLR_Msac_opt';
keys.sct(cs).color_option='monkeys_by_marker';
keys.sct(cs).tt.choices=[0,1];

% PreS
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_PreS_spaceCI_IX_Msac_opt';
keys.sct(cs).Y='ch_PreS_spaceCI_IX_Msac_opt';
keys.sct(cs).X_sig='in_PreS_spaceLR_Msac_opt';
keys.sct(cs).Y_sig='ch_PreS_spaceLR_Msac_opt';
keys.sct(cs).color_option='monkeys_by_marker';
keys.sct(cs).tt.choices=[0,1];

% PeriS
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_PeriS_spaceCI_IX_Msac_opt';
keys.sct(cs).Y='ch_PeriS_spaceCI_IX_Msac_opt';
keys.sct(cs).X_sig='in_PeriS_spaceLR_Msac_opt';
keys.sct(cs).Y_sig='ch_PeriS_spaceLR_Msac_opt';
keys.sct(cs).color_option='monkeys_by_marker';
keys.sct(cs).tt.choices=[0,1];

% TIhol
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_TIhol_spaceCI_IX_Msac_opt';
keys.sct(cs).Y='ch_TIhol_spaceCI_IX_Msac_opt';
keys.sct(cs).X_sig='in_TIhol_spaceLR_Msac_opt';
keys.sct(cs).Y_sig='ch_TIhol_spaceLR_Msac_opt';
keys.sct(cs).color_option='monkeys_by_marker';
keys.sct(cs).tt.choices=[0,1];

%% FR differences

% Cue
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_Cue_spaceLR_DF_Msac_opt';
keys.sct(cs).Y='ch_Cue_spaceLR_DF_Msac_opt';
keys.sct(cs).X_sig='in_Cue_spaceLR_Msac_opt';
keys.sct(cs).Y_sig='ch_Cue_spaceLR_Msac_opt';
keys.sct(cs).color_option='monkeys_by_marker';
keys.sct(cs).tt.choices=[0,1];

% MemE
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_MemE_spaceLR_DF_Msac_opt';
keys.sct(cs).Y='ch_MemE_spaceLR_DF_Msac_opt';
keys.sct(cs).X_sig='in_MemE_spaceLR_Msac_opt';
keys.sct(cs).Y_sig='ch_MemE_spaceLR_Msac_opt';
keys.sct(cs).color_option='monkeys_by_marker';
keys.sct(cs).tt.choices=[0,1];

% MemL
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_MemL_spaceLR_DF_Msac_opt';
keys.sct(cs).Y='ch_MemL_spaceLR_DF_Msac_opt';
keys.sct(cs).X_sig='in_MemL_spaceLR_Msac_opt';
keys.sct(cs).Y_sig='ch_MemL_spaceLR_Msac_opt';
keys.sct(cs).color_option='monkeys_by_marker';
keys.sct(cs).tt.choices=[0,1];

% PreS
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_PreS_spaceLR_DF_Msac_opt';
keys.sct(cs).Y='ch_PreS_spaceLR_DF_Msac_opt';
keys.sct(cs).X_sig='in_PreS_spaceLR_Msac_opt';
keys.sct(cs).Y_sig='ch_PreS_spaceLR_Msac_opt';
keys.sct(cs).color_option='monkeys_by_marker';
keys.sct(cs).tt.choices=[0,1];

% PeriS
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_PeriS_spaceLR_DF_Msac_opt';
keys.sct(cs).Y='ch_PeriS_spaceLR_DF_Msac_opt';
keys.sct(cs).X_sig='in_PeriS_spaceLR_Msac_opt';
keys.sct(cs).Y_sig='ch_PeriS_spaceLR_Msac_opt';
keys.sct(cs).color_option='monkeys_by_marker';
keys.sct(cs).tt.choices=[0,1];

% TIhol
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_TIhol_spaceLR_DF_Msac_opt';
keys.sct(cs).Y='ch_TIhol_spaceLR_DF_Msac_opt';
keys.sct(cs).X_sig='in_TIhol_spaceLR_Msac_opt';
keys.sct(cs).Y_sig='ch_TIhol_spaceLR_Msac_opt';
keys.sct(cs).color_option='monkeys_by_marker';
keys.sct(cs).tt.choices=[0,1];

%% Choice


%% Cue
% PrefIN CHOICE HEMI Cue pref vs unpreferred ALL
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_Cue_prefHO_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_Cue_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_Cue_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_Cue_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefIN CHOICE vs INSTRUCTED HEMI Cue ALL
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_Cue_prefHI_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_Cue_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_Cue_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_Cue_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefCH CHOICE HEMI Cue pref vs unpreferred ALL
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_Cue_prefHOch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_Cue_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_Cue_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_Cue_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefCH CHOICE vs INSTRUCTED HEMI Cue ALL
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_Cue_prefHIch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_Cue_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_Cue_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_Cue_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefIN CHOICE HEMI Cue pref vs unpreferred TUNED
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_Cue_prefHO_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_Cue_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_Cue_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_Cue_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_Cue_spaceLR_Msac_opt','-'};


% PrefIN CHOICE vs INSTRUCTED HEMI Cue TUNED
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_Cue_prefHI_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_Cue_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_Cue_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_Cue_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_Cue_spaceLR_Msac_opt','-'};


% PrefCH CHOICE HEMI Cue pref vs unpreferred TUNED
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_Cue_prefHOch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_Cue_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_Cue_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_Cue_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_Cue_spaceLR_Msac_opt','-'};


% PrefCH CHOICE vs INSTRUCTED HEMI Cue TUNED
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_Cue_prefHIch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_Cue_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_Cue_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_Cue_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_Cue_spaceLR_Msac_opt','-'};



%% MemE
% PrefIN CHOICE HEMI MemE pref vs unpreferred ALL
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_MemE_prefHO_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_MemE_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_MemE_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_MemE_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefIN CHOICE vs INSTRUCTED HEMI MemE ALL
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_MemE_prefHI_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_MemE_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_MemE_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_MemE_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefCH CHOICE HEMI MemE pref vs unpreferred ALL
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_MemE_prefHOch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_MemE_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_MemE_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_MemE_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefCH CHOICE vs INSTRUCTED HEMI MemE ALL
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_MemE_prefHIch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_MemE_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_MemE_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_MemE_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefIN CHOICE HEMI MemE pref vs unpreferred TUNED
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_MemE_prefHO_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_MemE_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_MemE_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_MemE_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_MemE_spaceLR_Msac_opt','-'};


% PrefIN CHOICE vs INSTRUCTED HEMI MemE TUNED
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_MemE_prefHI_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_MemE_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_MemE_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_MemE_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_MemE_spaceLR_Msac_opt','-'};


% PrefCH CHOICE HEMI MemE pref vs unpreferred TUNED
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_MemE_prefHOch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_MemE_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_MemE_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_MemE_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_MemE_spaceLR_Msac_opt','-'};


% PrefCH CHOICE vs INSTRUCTED HEMI MemE TUNED
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_MemE_prefHIch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_MemE_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_MemE_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_MemE_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_MemE_spaceLR_Msac_opt','-'};


%% MemL
% PrefIN CHOICE HEMI MemL pref vs unpreferred ALL
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_MemL_prefHO_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_MemL_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_MemL_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_MemL_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefIN CHOICE vs INSTRUCTED HEMI MemL ALL
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_MemL_prefHI_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_MemL_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_MemL_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_MemL_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefCH CHOICE HEMI MemL pref vs unpreferred ALL
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_MemL_prefHOch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_MemL_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_MemL_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_MemL_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefCH CHOICE vs INSTRUCTED HEMI MemL ALL
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_MemL_prefHIch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_MemL_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_MemL_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_MemL_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefIN CHOICE HEMI MemL pref vs unpreferred TUNED
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_MemL_prefHO_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_MemL_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_MemL_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_MemL_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_MemL_spaceLR_Msac_opt','-'};


% PrefIN CHOICE vs INSTRUCTED HEMI MemL TUNED
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_MemL_prefHI_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_MemL_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_MemL_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_MemL_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_MemL_spaceLR_Msac_opt','-'};


% PrefCH CHOICE HEMI MemL pref vs unpreferred TUNED
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_MemL_prefHOch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_MemL_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_MemL_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_MemL_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_MemL_spaceLR_Msac_opt','-'};


% PrefCH CHOICE vs INSTRUCTED HEMI MemL TUNED
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_MemL_prefHIch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_MemL_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_MemL_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_MemL_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_MemL_spaceLR_Msac_opt','-'};


%% PreS
% PrefIN CHOICE HEMI PreS pref vs unpreferred ALL
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_PreS_prefHO_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_PreS_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_PreS_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_PreS_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefIN CHOICE vs INSTRUCTED HEMI PreS ALL
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_PreS_prefHI_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_PreS_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_PreS_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_PreS_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefCH CHOICE HEMI PreS pref vs unpreferred ALL
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_PreS_prefHOch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_PreS_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_PreS_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_PreS_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefCH CHOICE vs INSTRUCTED HEMI PreS ALL
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_PreS_prefHIch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_PreS_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_PreS_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_PreS_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefIN CHOICE HEMI PreS pref vs unpreferred TUNED
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_PreS_prefHO_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_PreS_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_PreS_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_PreS_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_PreS_spaceLR_Msac_opt','-'};


% PrefIN CHOICE vs INSTRUCTED HEMI PreS TUNED
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_PreS_prefHI_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_PreS_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_PreS_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_PreS_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_PreS_spaceLR_Msac_opt','-'};


% PrefCH CHOICE HEMI PreS pref vs unpreferred TUNED
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_PreS_prefHOch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_PreS_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_PreS_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_PreS_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_PreS_spaceLR_Msac_opt','-'};


% PrefCH CHOICE vs INSTRUCTED HEMI PreS TUNED
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_PreS_prefHIch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_PreS_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_PreS_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_PreS_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_PreS_spaceLR_Msac_opt','-'};


%% Pre2
% PrefIN CHOICE HEMI Pre2 pref vs unpreferred ALL
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_Pre2_prefHO_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_Pre2_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_Pre2_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_Pre2_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefIN CHOICE vs INSTRUCTED HEMI Pre2 ALL
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_Pre2_prefHI_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_Pre2_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_Pre2_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_Pre2_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefCH CHOICE HEMI Pre2 pref vs unpreferred ALL
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_Pre2_prefHOch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_Pre2_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_Pre2_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_Pre2_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefCH CHOICE vs INSTRUCTED HEMI Pre2 ALL
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_Pre2_prefHIch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_Pre2_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_Pre2_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_Pre2_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefIN CHOICE HEMI Pre2 pref vs unpreferred TUNED
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_Pre2_prefHO_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_Pre2_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_Pre2_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_Pre2_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_Pre2_spaceLR_Msac_opt','-'};


% PrefIN CHOICE vs INSTRUCTED HEMI Pre2 TUNED
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_Pre2_prefHI_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_Pre2_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_Pre2_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_Pre2_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_Pre2_spaceLR_Msac_opt','-'};


% PrefCH CHOICE HEMI Pre2 pref vs unpreferred TUNED
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_Pre2_prefHOch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_Pre2_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_Pre2_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_Pre2_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_Pre2_spaceLR_Msac_opt','-'};


% PrefCH CHOICE vs INSTRUCTED HEMI Pre2 TUNED
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_Pre2_prefHIch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_Pre2_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_Pre2_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_Pre2_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_Pre2_spaceLR_Msac_opt','-'};


%% PeriS
% PrefIN CHOICE HEMI PeriS pref vs unpreferred ALL
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_PeriS_prefHO_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_PeriS_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_PeriS_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_PeriS_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefIN CHOICE vs INSTRUCTED HEMI PeriS ALL
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_PeriS_prefHI_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_PeriS_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_PeriS_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_PeriS_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefCH CHOICE HEMI PeriS pref vs unpreferred ALL
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_PeriS_prefHOch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_PeriS_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_PeriS_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_PeriS_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefCH CHOICE vs INSTRUCTED HEMI PeriS ALL
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_PeriS_prefHIch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_PeriS_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_PeriS_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_PeriS_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefIN CHOICE HEMI PeriS pref vs unpreferred TUNED
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_PeriS_prefHO_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_PeriS_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_PeriS_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_PeriS_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_PeriS_spaceLR_Msac_opt','-'};


% PrefIN CHOICE vs INSTRUCTED HEMI PeriS TUNED
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_PeriS_prefHI_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_PeriS_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_PeriS_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_PeriS_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_PeriS_spaceLR_Msac_opt','-'};


% PrefCH CHOICE HEMI PeriS pref vs unpreferred TUNED
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_PeriS_prefHOch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_PeriS_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_PeriS_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_PeriS_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_PeriS_spaceLR_Msac_opt','-'};


% PrefCH CHOICE vs INSTRUCTED HEMI PeriS TUNED
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_PeriS_prefHIch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_PeriS_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_PeriS_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_PeriS_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_PeriS_spaceLR_Msac_opt','-'};

%% Peri2
% PrefIN CHOICE HEMI Peri2 pref vs unpreferred ALL
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_Peri2_prefHO_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_Peri2_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_Peri2_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_Peri2_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefIN CHOICE vs INSTRUCTED HEMI Peri2 ALL
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_Peri2_prefHI_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_Peri2_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_Peri2_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_Peri2_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefCH CHOICE HEMI Peri2 pref vs unpreferred ALL
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_Peri2_prefHOch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_Peri2_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_Peri2_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_Peri2_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefCH CHOICE vs INSTRUCTED HEMI Peri2 ALL
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_Peri2_prefHIch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_Peri2_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_Peri2_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_Peri2_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefIN CHOICE HEMI Peri2 pref vs unpreferred TUNED
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_Peri2_prefHO_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_Peri2_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_Peri2_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_Peri2_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_Peri2_spaceLR_Msac_opt','-'};


% PrefIN CHOICE vs INSTRUCTED HEMI Peri2 TUNED
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_Peri2_prefHI_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_Peri2_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_Peri2_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_Peri2_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_Peri2_spaceLR_Msac_opt','-'};


% PrefCH CHOICE HEMI Peri2 pref vs unpreferred TUNED
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_Peri2_prefHOch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_Peri2_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_Peri2_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_Peri2_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_Peri2_spaceLR_Msac_opt','-'};


% PrefCH CHOICE vs INSTRUCTED HEMI Peri2 TUNED
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_Peri2_prefHIch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_Peri2_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_Peri2_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_Peri2_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_Peri2_spaceLR_Msac_opt','-'};

%% TIhol
% PrefIN CHOICE HEMI TIhol pref vs unpreferred ALL
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_TIhol_prefHO_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_TIhol_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_TIhol_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_TIhol_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefIN CHOICE vs INSTRUCTED HEMI TIhol ALL
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_TIhol_prefHI_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_TIhol_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_TIhol_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_TIhol_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefCH CHOICE HEMI TIhol pref vs unpreferred ALL
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_TIhol_prefHOch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_TIhol_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_TIhol_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_TIhol_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefCH CHOICE vs INSTRUCTED HEMI TIhol ALL
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_TIhol_prefHIch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_TIhol_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_TIhol_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_TIhol_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];


% PrefIN CHOICE HEMI TIhol pref vs unpreferred TUNED
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_TIhol_prefHO_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_TIhol_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_TIhol_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_TIhol_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_TIhol_spaceLR_Msac_opt','-'};


% PrefIN CHOICE vs INSTRUCTED HEMI TIhol TUNED
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_TIhol_prefHI_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_TIhol_prefHI_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_TIhol_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_TIhol_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_TIhol_spaceLR_Msac_opt','-'};


% PrefCH CHOICE HEMI TIhol pref vs unpreferred TUNED
cs=cs+1; 
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='ch_AH_TIhol_prefHOch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_TIhol_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_TIhol_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_TIhol_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_TIhol_spaceLR_Msac_opt','-'};


% PrefCH CHOICE vs INSTRUCTED HEMI TIhol TUNED
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_TIhol_prefHIch_DF_Msac_opt';
keys.sct(cs).Y='ch_AH_TIhol_prefHIch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_TIhol_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_TIhol_epoch_Msac_opt';
keys.sct(cs).color_option='ENSU_as_color';
keys.sct(cs).absolutes=1;
keys.sct(cs).tt.choices=[0,1];
keys.sct(cs).tt.unselect={'in_TIhol_spaceLR_Msac_opt','-'};

% VMI plots
% VMI ipsi versus contra POST
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='VMI_postEN_IS_Msac_opt';
keys.sct(cs).Y='VMI_postEN_CS_Msac_opt';
keys.sct(cs).X_sig='VMI_postEN_Msac_opt';
keys.sct(cs).Y_sig='VMI_postEN_Msac_opt';
keys.sct(cs).hist_column='VMI_postEN_Msac_opt';
keys.sct(cs).categories={'visual_Msac_opt';'visuomotor_Msac_opt';'motor_Msac_opt'};
keys.sct(cs).color_option='FR_as_color';

% VMI ipsi vs contra PERI
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='VMI_periEN_IS_Msac_opt';
keys.sct(cs).Y='VMI_periEN_CS_Msac_opt';
keys.sct(cs).X_sig='VMI_periEN_Msac_opt';
keys.sct(cs).Y_sig='VMI_periEN_Msac_opt';
keys.sct(cs).hist_column='VMI_periEN_Msac_opt';
keys.sct(cs).categories={'visual_Msac_opt';'visuomotor_Msac_opt';'motor_Msac_opt'};
keys.sct(cs).color_option='FR_as_color';

% IPSI peri2 vs cue enhancement, VMI as colors
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_IS_Cue_epoch_DF_Msac_opt';
keys.sct(cs).Y='in_AH_IS_Peri2_epoch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_IS_Cue_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_IS_Peri2_epoch_Msac_opt';
keys.sct(cs).VMI='VMI_periEN_IS_Msac_opt';
keys.sct(cs).hist_column='VMI_periEN_Msac_opt';
keys.sct(cs).categories={'visual_peri_Msac_opt';'visuomotor_peri_Msac_opt';'motor_peri_Msac_opt'};
keys.sct(cs).color_option='VMI_as_color';
keys.sct(cs).logarithmic_scale=1;

% CONTRA peri2 vs cue enhancement, VMI as colors
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_CS_Cue_epoch_DF_Msac_opt';
keys.sct(cs).Y='in_AH_CS_Peri2_epoch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_CS_Cue_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_CS_Peri2_epoch_Msac_opt';
keys.sct(cs).VMI='VMI_periEN_CS_Msac_opt';
keys.sct(cs).hist_column='VMI_periEN_CS_Msac_opt';
keys.sct(cs).categories={'visual_peri_Msac_opt';'visuomotor_peri_Msac_opt';'motor_peri_Msac_opt'};
keys.sct(cs).color_option='VMI_as_color';
keys.sct(cs).logarithmic_scale=1;

% IPSI periS vs cue enhancement, VMI as colors
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_IS_Cue_epoch_DF_Msac_opt';
keys.sct(cs).Y='in_AH_IS_PeriS_epoch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_IS_Cue_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_IS_PeriS_epoch_Msac_opt';
keys.sct(cs).VMI='VMI_periEN_IS_Msac_opt';
keys.sct(cs).hist_column='VMI_periEN_IS_Msac_opt';
keys.sct(cs).categories={'visual_peri_Msac_opt';'visuomotor_peri_Msac_opt';'motor_peri_Msac_opt'};
keys.sct(cs).color_option='VMI_as_color';
keys.sct(cs).logarithmic_scale=1;

% CONTRA periS vs cue enhancement, VMI as colors
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_CS_Cue_epoch_DF_Msac_opt';
keys.sct(cs).Y='in_AH_CS_PeriS_epoch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_CS_Cue_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_CS_PeriS_epoch_Msac_opt';
keys.sct(cs).VMI='VMI_periEN_CS_Msac_opt';
keys.sct(cs).hist_column='VMI_periEN_CS_Msac_opt';
keys.sct(cs).categories={'visual_peri_Msac_opt';'visuomotor_peri_Msac_opt';'motor_peri_Msac_opt'};
keys.sct(cs).color_option='VMI_as_color';
keys.sct(cs).logarithmic_scale=1;

% IPSI post vs cue enhancement, VMI as colors
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_IS_Cue_epoch_DF_Msac_opt';
keys.sct(cs).Y='in_AH_IS_TIhol_epoch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_IS_Cue_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_IS_TIhol_epoch_Msac_opt';
keys.sct(cs).VMI='VMI_postEN_IS_Msac_opt';
keys.sct(cs).hist_column='VMI_postEN_IS_Msac_opt';
keys.sct(cs).categories={'visual_Msac_opt';'visuomotor_Msac_opt';'motor_Msac_opt'};
keys.sct(cs).color_option='VMI_as_color';
keys.sct(cs).logarithmic_scale=1;

% CONTRA post vs cue enhancement, VMI as colors
cs=cs+1;
keys.sct(cs).tt.tasktypes={'Msac_opt'};
keys.sct(cs).X='in_AH_CS_Cue_epoch_DF_Msac_opt';
keys.sct(cs).Y='in_AH_CS_TIhol_epoch_DF_Msac_opt';
keys.sct(cs).X_sig='in_AH_CS_Cue_epoch_Msac_opt';
keys.sct(cs).Y_sig='in_AH_CS_TIhol_epoch_Msac_opt';
keys.sct(cs).VMI='VMI_postEN_CS_Msac_opt';
keys.sct(cs).hist_column='VMI_postEN_CS_Msac_opt';
keys.sct(cs).categories={'visual_Msac_opt';'visuomotor_Msac_opt';'motor_Msac_opt'};
keys.sct(cs).color_option='VMI_as_color';
keys.sct(cs).logarithmic_scale=1;




