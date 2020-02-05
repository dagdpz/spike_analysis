keys.project_versions={''};
keys.project_version='20171214';
keys.filelist_formatted={};

keys.plot.single_cells                  =1;         % perform single cell plotting

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
keys.batching.monkeys               ={'Linus','Curius'};
keys.Curius.date                    ='[20150515 20150814]';
keys.Linus.date                     ='[20150508 20151030]';
keys.Linus.color      =[0 0 255]/255;
keys.Curius.color     =[255 0 0]/255;

%% cell count settings
keys.cc.choices             =[0];
keys.cc.factors                 ={'epoch','space'};
keys.cc.conditions_to_plot      ={'Msac','Vsac'};
keys.cc.plot_types              ={'visuomotor','per_epoch','per_task','space_and_epoch'};
keys.cc.epochsE.Msac        ={'INI', 'Fhol','Cue','MemE','MemL','PreS','PeriS','TIhol','Tons','Thol'}';
keys.cc.epochsS.Msac        ={'Cue','MemL','PreS','PeriS','TIhol','Tons','Thol'}';

%% for cell counts in general, but could also be used for PSTH!?
keys.tt.epoch_criterion             ='SxE or epoch only';
keys.tt.space_criterion             ='interaction or space only';

%% epochs
keys.EPOCHS_PER_TYPE{2}={...
    'INI',      2,	-0.4,	-0.1,   'INI';...
    'Fhol',     4,	-0.3,	0,      'INI';...
    'Cue',      4,	0.07,	0.17,   'INI';...
    'PreS',     60,	-0.1,	-0.01,  'INI';...
    'PeriS',    60,	-0.01,  0.05,   'INI';...
    'Tacq',     5,	0,      0.1,    'INI';...
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
    'Initiation',   2,	-0.5,	0;...
    'Visual onset', 4,	-0.8,   0.17;...
    'Saccade',      60,	-0.01,  0.22;...
    'T hold',       20,	-0.3,   0.1;...
    };

keys.ANOVAS_PER_TYPE(2).spaceLR            ={'Cue','PreS','PeriS','Tacq','Thol'}';
keys.ANOVAS_PER_TYPE(2).positions          ={'Cue','PreS','PeriS','Tacq','Thol'}';
keys.ANOVAS_PER_TYPE(2).hands              ={'Fhol','Cue','PreS','PeriS','Tacq','Thol'}';
keys.ANOVAS_PER_TYPE(2).SxH                ={'Cue','PreS','PeriS','Tacq','Thol'}';
keys.ANOVAS_PER_TYPE(2).main              ={'INI','Fhol','Cue','PreS','PeriS','Tacq','Thol'}';

keys.EPOCHS_PER_TYPE{3}={...
    'INI',      2,	-0.4,	-0.1,   'INI';...
    'Fhol',     6,	-0.3,	0,      'INI';...
    'Cue',      6,	0.07,   0.17,   'INI';...
    'MemE',     7,	0.07,   0.17,   'INI';... 
    'MemL',     9,	-0.3, 	-0,     'INI';... 
    'PreS',     60,	-0.1, 	-0.01,  'INI';...
    'Pre2',     60,	-0.1, 	-0.01,  'INI';...
    'PeriS',	60,	-0.01,  0.05,   'INI';...
    'Peri2',	60,	-0.01,  0.05,   'INI';...
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
    'Fhol' 'MemE';...
    'Fhol' 'MemL';...
    'MemL' 'PreS';...
    'Fhol' 'Pre2';...
    'MemL' 'PeriS';...
    'Fhol' 'Peri2';...
    'MemL' 'TIhol';...
    'MemL' 'Tons';...
    'MemL' 'Thol';...
    };

keys.ANOVAS_PER_TYPE(3).spaceLR            ={'Cue','MemE','MemL','PreS','Pre2','PeriS','Peri2','TIhol','Tons','Thol'}';
keys.ANOVAS_PER_TYPE(3).positions          ={'Cue','MemE','MemL','PreS','Pre2','PeriS','Peri2','TIhol','Tons','Thol'}';
keys.ANOVAS_PER_TYPE(3).hands              ={'Fhol','Cue','MemE','MemL','PreS','Pre2','PeriS','Peri2','TIhol','Tons','Thol'}';
keys.ANOVAS_PER_TYPE(3).SxH                ={'Cue','MemE','MemL','PreS','Pre2','PeriS','Peri2','TIhol','Tons','Thol'}';
keys.ANOVAS_PER_TYPE(3).main              ={'INI','Fhol','Cue','MemE','MemL','PreS','PeriS','Tons','Thol'}';

%% tuning onset settings
cc=0;

cc=cc+1;% 'Msac epoch tuning';
ce=0;
keys.ons(cc).comparisons_title       = 'Msac epoch tuning';
keys.ons(cc).group_parameter         = 'ungrouped';
keys.ons(cc).conditions_to_plot      = {'Msac'};
keys.ons(cc).choices=0; %for cell exclusion
keys.ons(cc).hands=0; %for cell exclusion
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1 1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[-1 1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=[0];
keys.ons(cc).comparisons_per_effector(ce).choice{2}=[0];
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Cue', 0, 0.28};
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
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Saccade', -0.2, 0.2};
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
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Hold', 0, 0.5};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.EP_EN/255; keys.colors.EP_SU/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='MemL';
keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Hold';

cc=cc+1;% 1 ungrouped baseline subtraction
ce=0;
keys.ons(cc).comparisons_title       = 'Msac instructed space tuning';
keys.ons(cc).group_parameter         = 'ungrouped';
keys.ons(cc).conditions_to_plot      = {'Msac'};
keys.ons(cc).choices=0; %for cell exclusion
keys.ons(cc).hands=0; %for cell exclusion
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
keys.ons(cc).comparisons_per_effector(ce).choice{2}=0;
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Cue', 0, 0.28};
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
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Saccade', -0.2, 0.2};
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
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Hold', 0, 0.5};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IS_IN/255; keys.colors.NH_CS_IN/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='MemL';
keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Hold';

cc=cc+1;% 1 ungrouped baseline subtraction
ce=0;
keys.ons(cc).comparisons_title       = 'Msac choice tuning';
keys.ons(cc).group_parameter         = 'ungrouped';
keys.ons(cc).conditions_to_plot      = {'Msac'};
keys.ons(cc).hands=0; %for cell exclusion
keys.ons(cc).choices=[0 1]; %for cell exclusion
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[-1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
keys.ons(cc).comparisons_per_effector(ce).choice{2}=1;
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Cue', 0, 0.28};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IN/255; keys.colors.NH_CH/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
keys.ons(cc).comparisons_per_effector(ce).title='Ipsilateral';
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
keys.ons(cc).comparisons_per_effector(ce).choice{2}=1;
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Cue', 0, 0.28};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IN/255; keys.colors.NH_CH/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
keys.ons(cc).comparisons_per_effector(ce).title='Contralateral';
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[-1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
keys.ons(cc).comparisons_per_effector(ce).choice{2}=1;
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Cue', 0.28, 0.78};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IN/255; keys.colors.NH_CH/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
keys.ons(cc).comparisons_per_effector(ce).title='Ipsilateral';
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
keys.ons(cc).comparisons_per_effector(ce).choice{2}=1;
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Cue', 0.28, 0.78};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IN/255; keys.colors.NH_CH/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
keys.ons(cc).comparisons_per_effector(ce).title='Contralateral';


cc=cc+1;% msac choice preference
ce=0;
keys.ons(cc).comparisons_title       = 'Msac choice preference';
keys.ons(cc).group_parameter         = 'ungrouped';
keys.ons(cc).conditions_to_plot      = {'Msac'};
keys.ons(cc).hands=0; %for cell exclusion
keys.ons(cc).choices=[1]; %for cell exclusion
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=1;
keys.ons(cc).comparisons_per_effector(ce).choice{2}=1;
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Cue', 0, 0.28};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.IS_CH/255; keys.colors.CS_CH/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Cue';
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=1;
keys.ons(cc).comparisons_per_effector(ce).choice{2}=1;
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Cue', 0.28, 0.78};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.IS_CH/255; keys.colors.CS_CH/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Cue';

cc=cc+1;% 1 ungrouped baseline subtraction
ce=0;
keys.ons(cc).comparisons_title       = 'Vsac instructed space tuning';
keys.ons(cc).group_parameter         = 'ungrouped';
keys.ons(cc).conditions_to_plot      = {'Vsac'};
keys.ons(cc).choices=0; %for cell exclusion
keys.ons(cc).hands=0; %for cell exclusion
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
keys.ons(cc).comparisons_per_effector(ce).choice{2}=0;
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Visual onset', 0, 0.28};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IS_IN/255; keys.colors.NH_CS_IN/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Visual onset';
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
keys.ons(cc).comparisons_per_effector(ce).choice{2}=0;
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Saccade', -0.2, 0.2};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IS_IN/255; keys.colors.NH_CS_IN/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Saccade';

cc=cc+1;% 1 ungrouped baseline subtraction
ce=0;
keys.ons(cc).comparisons_title       = 'Vsac choice tuning';
keys.ons(cc).group_parameter         = 'ungrouped';
keys.ons(cc).conditions_to_plot      = {'Vsac'};
keys.ons(cc).choices=[0 1]; %for cell exclusion
keys.ons(cc).hands=0; %for cell exclusion
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[-1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
keys.ons(cc).comparisons_per_effector(ce).choice{2}=1;
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Visual onset', 0, 0.28};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IN/255; keys.colors.NH_CH/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
keys.ons(cc).comparisons_per_effector(ce).title='Ipsilateral';
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
keys.ons(cc).comparisons_per_effector(ce).choice{2}=1;
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Visual onset', 0, 0.28};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IN/255; keys.colors.NH_CH/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
keys.ons(cc).comparisons_per_effector(ce).title='Contralateral';

cc=cc+1;% choice preference
ce=0;
keys.ons(cc).comparisons_title       = 'Vsac choice preference';
keys.ons(cc).group_parameter         = 'ungrouped';
keys.ons(cc).conditions_to_plot      = {'Vsac'};
keys.ons(cc).choices=[1]; %for cell exclusion
keys.ons(cc).hands=0; %for cell exclusion
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=1;
keys.ons(cc).comparisons_per_effector(ce).choice{2}=1;
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Visual onset', 0, 0.28};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.IS_CH/255; keys.colors.CS_CH/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
keys.ons(cc).comparisons_per_effector(ce).title='Aligned to cue';

cc=cc+1;% 1 ungrouped baseline subtraction
ce=0;
keys.ons(cc).comparisons_title       = 'Vsac epoch space tuning';
keys.ons(cc).group_parameter         = 'ungrouped';
keys.ons(cc).conditions_to_plot      = {'Vsac'};
keys.ons(cc).choices=[0]; %for cell exclusion
keys.ons(cc).hands=0; %for cell exclusion
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1 1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[-1 1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=[0];
keys.ons(cc).comparisons_per_effector(ce).choice{2}=[0];
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Visual onset', 0, 0.28};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.EP_EN/255; keys.colors.EP_SU/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='Fhol';
keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Visual onset';
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1 1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[-1 1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=[0];
keys.ons(cc).comparisons_per_effector(ce).choice{2}=[0];
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Saccade', -0.2, 0.2};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.EP_EN/255; keys.colors.EP_SU/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='Fhol';
keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Saccade';
 
%% population PSTH settings
cc=0;

% All
cc=cc+1;% 1 ungrouped baseline subtraction
keys.pop(cc).normalization           = 'none';
keys.pop(cc).group_parameter         = 'ungrouped';
keys.pop(cc).conditions_to_plot      = {'Vsac'};
keys.pop(cc).FR_subtract_baseline    = 1;
keys.pop(cc).choices                 = [0 1];
keys.pop(cc).hands                   = [0];
cc=cc+1;% 2 ungrouped baseline subtraction
keys.pop(cc).normalization           = 'none';
keys.pop(cc).group_parameter         = 'ungrouped';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).FR_subtract_baseline    = 1;
keys.pop(cc).choices                 = [0 1];
keys.pop(cc).hands                   = [0];
cc=cc+1;% 3 ungrouped raw
keys.pop(cc).normalization           = 'none';
keys.pop(cc).group_parameter         = 'ungrouped';
keys.pop(cc).conditions_to_plot      = {'Vsac','Msac'};
keys.pop(cc).choices                 = [0 1];
keys.pop(cc).hands                   = [0];

% cathegories
cc=cc+1;% fixation cells?
keys.pop(cc).normalization           = 'none';
keys.pop(cc).group_parameter         = 'fixation_only_Msac_opt';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).FR_subtract_baseline    = 1;
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
cc=cc+1;% fixation cells?
keys.pop(cc).normalization           = 'none';
keys.pop(cc).group_parameter         = 'fixation_and_sac_suppression_Msac_opt';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).FR_subtract_baseline    = 1;
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
cc=cc+1;% visual response?
keys.pop(cc).normalization           = 'none';
keys.pop(cc).group_parameter         = 'visual_Msac_opt';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).FR_subtract_baseline    = 1;
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
cc=cc+1;% motor output only?
keys.pop(cc).normalization           = 'none';
keys.pop(cc).group_parameter         = 'motor_Msac_opt';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).FR_subtract_baseline    = 1;
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
cc=cc+1;% visuomotor cells???
keys.pop(cc).normalization           = 'none';
keys.pop(cc).group_parameter         = 'visuomotor_Msac_opt';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).FR_subtract_baseline    = 1;
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];


cc=cc+1;% 2 Cue spacial tuning
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_Cue_spaceLR_Msac_opt';
keys.pop(cc).conditions_to_plot         = {'Msac'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
cc=cc+1;% 2 Cue spacial tuning
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_NH_Cue_epoch_Msac_opt';
keys.pop(cc).conditions_to_plot         = {'Msac'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
cc=cc+1;% gaze cells?
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_TIhol_spaceLR_Msac_opt';
keys.pop(cc).epoch_RF                   = 'TIhol';
keys.pop(cc).conditions_to_plot         = {'Msac'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
cc=cc+1;% gaze cells?
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_NH_TIhol_epoch_Msac_opt';
keys.pop(cc).epoch_RF                   = 'TIhol';
keys.pop(cc).conditions_to_plot         = {'Msac'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
cc=cc+1;% gaze cells?
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_Thol_spaceLR_Msac_opt';
keys.pop(cc).epoch_RF                   = 'Thol';
keys.pop(cc).conditions_to_plot         = {'Msac'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
cc=cc+1;% gaze cells?
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_NH_Thol_epoch_Msac_opt';
keys.pop(cc).epoch_RF                   = 'Thol';
keys.pop(cc).conditions_to_plot         = {'Msac'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
cc=cc+1;% 3 Target onset space tuning
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_Tons_spaceLR_Msac_opt';
keys.pop(cc).epoch_RF                   = 'Tons';
keys.pop(cc).conditions_to_plot         = {'Msac'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
cc=cc+1;% 3 Target onset enhancement
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_NH_Tons_epoch_Msac_opt';
keys.pop(cc).epoch_RF                   = 'Tons';
keys.pop(cc).conditions_to_plot         = {'Msac'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];

% ramping ...
cc=cc+1;% ramping cells???
keys.pop(cc).normalization           = 'none';
keys.pop(cc).group_parameter         = 'in_NH_MemL_epoch_Msac_opt';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).FR_subtract_baseline    = 1;
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
cc=cc+1;% ramping cells???
keys.pop(cc).normalization           = 'none';
keys.pop(cc).group_parameter         = 'in_MemL_spaceLR_Msac_opt';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).FR_subtract_baseline    = 1;
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
cc=cc+1;% 3 Ramping in 2 steps
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'Ramp';
keys.pop(cc).epoch_RF                   = 'Peri2';
keys.pop(cc).conditions_to_plot         = {'Msac'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).group_excluded             ={'-bi','bi-','--','-su','-en','en-','su-'};
keys.pop(cc).combine_tuning_properties  ={'Ramp','in_NH_MemL_epoch_Msac_opt','in_NH_PreS_epoch_Msac_opt'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
cc=cc+1;% Pre-saccadic enhancement7suppression
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_NH_PreS_epoch_Msac_opt';
keys.pop(cc).epoch_RF                   = 'PreS';
keys.pop(cc).conditions_to_plot         = {'Msac'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).group_excluded             ={'','-','bi'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
cc=cc+1;% Peri-saccadic enhancement7suppression
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_NH_PeriS_epoch_Msac_opt';
keys.pop(cc).epoch_RF                   = 'PeriS';
keys.pop(cc).conditions_to_plot         = {'Msac'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).group_excluded             ={'','-','bi'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];

cc=cc+1;% 3 Pre-saccadic enhancement/suppression... to be removed, this is to reconstruct the plots that we had previously
keys.pop(cc).normalization              = 'none';
keys.pop(cc).group_parameter            = 'in_NH_Peri2_epoch_Msac_opt';
keys.pop(cc).epoch_RF                   = 'Peri2';
keys.pop(cc).conditions_to_plot         = {'Msac'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).group_excluded             ={'','-','bi'};
keys.pop(cc).combine_tuning_properties  ={'E_ExS','in_epoch_main_Msac_opt','in_ExS_Msac_opt'};
keys.pop(cc).unselect                   ={'E_ExS','00'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];




% cathegories
cc=cc+1;% fixation cells?
keys.pop(cc).group_parameter         = 'fixation_only_Msac_opt';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_condition';
keys.pop(cc).epoch_for_normalization = 'Fhol';
keys.pop(cc).FR_subtract_baseline    = 0;
cc=cc+1;% fixation cells?
keys.pop(cc).group_parameter         = 'ixation_and_sac_suppression_Msac_opt';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_condition';
keys.pop(cc).epoch_for_normalization = 'Fhol';
keys.pop(cc).FR_subtract_baseline    = 0;
cc=cc+1;% visual response?
keys.pop(cc).group_parameter         = 'visual_Msac_opt';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_condition';
keys.pop(cc).epoch_for_normalization = 'Cue';
keys.pop(cc).FR_subtract_baseline    = 0;
cc=cc+1;% motor output only?
keys.pop(cc).group_parameter         = 'motor_Msac_opt';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_condition';
keys.pop(cc).epoch_for_normalization = 'PeriS';
keys.pop(cc).FR_subtract_baseline    = 0;
cc=cc+1;% visuomotor cells???
keys.pop(cc).group_parameter         = 'visuomotor_Msac_opt';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_condition';
keys.pop(cc).epoch_for_normalization = 'Cue';
keys.pop(cc).FR_subtract_baseline    = 0;


cc=cc+1;% 2 Cue spacial tuning
keys.pop(cc).group_parameter         = 'in_Cue_spaceLR_Msac_opt';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_condition';
keys.pop(cc).epoch_for_normalization = 'Cue';
keys.pop(cc).FR_subtract_baseline    = 0;
cc=cc+1;% 2 Cue spacial tuning
keys.pop(cc).group_parameter         = 'in_NH_Cue_epoch_Msac_opt';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_condition';
keys.pop(cc).epoch_for_normalization = 'Cue';
keys.pop(cc).FR_subtract_baseline    = 0;
cc=cc+1;% gaze cells?
keys.pop(cc).group_parameter         = 'in_TIhol_spaceLR_Msac_opt';
keys.pop(cc).epoch_RF                = 'TIhol';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_condition';
keys.pop(cc).epoch_for_normalization = 'TIhol';
keys.pop(cc).FR_subtract_baseline    = 0;
cc=cc+1;% gaze cells?
keys.pop(cc).group_parameter         = 'in_NH_TIhol_epoch_Msac_opt';
keys.pop(cc).epoch_RF                = 'TIhol';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_condition';
keys.pop(cc).epoch_for_normalization = 'TIhol';
keys.pop(cc).FR_subtract_baseline    = 0;
cc=cc+1;% gaze cells?
keys.pop(cc).group_parameter         = 'in_Thol_spaceLR_Msac_opt';
keys.pop(cc).epoch_RF                = 'Thol';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_condition';
keys.pop(cc).epoch_for_normalization = 'Thol';
keys.pop(cc).FR_subtract_baseline    = 0;
cc=cc+1;% gaze cells?
keys.pop(cc).group_parameter         = 'in_NH_Thol_epoch_Msac_opt';
keys.pop(cc).epoch_RF                = 'Thol';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_condition';
keys.pop(cc).epoch_for_normalization = 'Thol';
keys.pop(cc).FR_subtract_baseline    = 0;
cc=cc+1;% gaze cells?
keys.pop(cc).group_parameter         = 'in_Tons_spaceLR_Msac_opt';
keys.pop(cc).epoch_RF                = 'Tons';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_condition';
keys.pop(cc).epoch_for_normalization = 'Tons';
keys.pop(cc).FR_subtract_baseline    = 0;
cc=cc+1;% 3 Target onset
keys.pop(cc).group_parameter         = 'in_NH_Tons_epoch_Msac_opt';
keys.pop(cc).epoch_RF                = 'Tons';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_condition';
keys.pop(cc).epoch_for_normalization = 'Tons';
keys.pop(cc).FR_subtract_baseline    = 0;

% ramping ...
cc=cc+1;% ramping cells???
keys.pop(cc).group_parameter         = 'in_NH_MemL_epoch_Msac_opt';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_condition';
keys.pop(cc).epoch_for_normalization = 'MemL';
keys.pop(cc).FR_subtract_baseline    = 0;
cc=cc+1;% ramping cells???
keys.pop(cc).group_parameter         = 'in_MemL_spaceLR_Msac_opt';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_condition';
keys.pop(cc).epoch_for_normalization = 'MemL';
keys.pop(cc).FR_subtract_baseline    = 0;
cc=cc+1;% 3 Ramping in 2 steps
keys.pop(cc).group_parameter            = 'Ramp';
keys.pop(cc).epoch_RF                   = 'Peri2';
keys.pop(cc).conditions_to_plot         = {'Msac'};
keys.pop(cc).group_excluded             ={'-bi','bi-','--','-su','-en','en-','su-'};
keys.pop(cc).combine_tuning_properties  ={'Ramp','in_NH_MemL_epoch_Msac_opt','in_NH_PreS_epoch_Msac_opt'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_condition';
keys.pop(cc).epoch_for_normalization = 'Peri2';
keys.pop(cc).FR_subtract_baseline    = 0;
cc=cc+1;% Pre-saccadic enhancement7suppression
keys.pop(cc).group_parameter         = 'in_NH_PreS_epoch_Msac_opt';
keys.pop(cc).epoch_RF                = 'PreS';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).group_excluded          ={'','-','bi'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_condition';
keys.pop(cc).epoch_for_normalization = 'PreS';
keys.pop(cc).FR_subtract_baseline    = 0;
cc=cc+1;% Peri-saccadic enhancement7suppression
keys.pop(cc).group_parameter            = 'in_NH_PeriS_epoch_Msac_opt';
keys.pop(cc).epoch_RF                   = 'PeriS';
keys.pop(cc).conditions_to_plot         = {'Msac'};
keys.pop(cc).FR_subtract_baseline       = 1;
keys.pop(cc).group_excluded             ={'','-','bi'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_condition';
keys.pop(cc).epoch_for_normalization = 'PeriS';
keys.pop(cc).FR_subtract_baseline    = 0;

cc=cc+1;% 3 Pre-saccadic enhancement/suppression... to be removed, this is to reconstruct the plots that we had previously
keys.pop(cc).group_parameter            = 'in_NH_Peri2_epoch_Msac_opt';
keys.pop(cc).epoch_RF                   = 'Peri2';
keys.pop(cc).conditions_to_plot         = {'Msac'};
keys.pop(cc).group_excluded             ={'','-','bi'};
keys.pop(cc).combine_tuning_properties  ={'E_ExS','in_epoch_main_Msac_opt','in_ExS_Msac_opt'};
keys.pop(cc).unselect                   ={'E_ExS','00'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_condition';
keys.pop(cc).epoch_for_normalization = 'Peri2';
keys.pop(cc).FR_subtract_baseline    = 0;


cc=cc+1;% 1 ungrouped baseline subtraction
keys.pop(cc).normalization           = 'by_effector';
keys.pop(cc).group_parameter         = 'ungrouped';
keys.pop(cc).conditions_to_plot      = {'Vsac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).FR_subtract_baseline    = 0;
keys.pop(cc).epoch_for_normalization = 'Cue';

cc=cc+1;% 1 ungrouped baseline subtraction
keys.pop(cc).normalization           = 'by_effector';
keys.pop(cc).group_parameter         = 'ungrouped';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).FR_subtract_baseline    = 0;
keys.pop(cc).epoch_for_normalization = 'Cue';

