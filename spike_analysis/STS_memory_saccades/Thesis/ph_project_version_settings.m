keys.project_versions={''};
keys.project_version='Thesis';
keys.filelist_formatted={};

keys.plot.single_cells                  =0;         % perform single cell plotting

%% to check carefully
keys.position_and_plotting_arrangements         ={'options'};
keys.cal.units_from_sorting_table       =1;                         % exclude units that are not in the sorting table

%% computation settings
keys.cal.effectors                  =[0];
keys.cal.reach_hand                 =[0];
keys.cal.types                      =[3];

keys.cal.stablity                       =[1];                     % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the table
keys.cal.single_rating                  =[1,2,3];                   % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the tabl
keys.cal.SNR_rating                     =[1,2,3];                   % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the tabl


%% batching
keys.batching.combine_monkeys       =1;
keys.batching.targets               ={'ST'};
keys.batching.monkeys               ={'Cornelius'};
keys.Cornelius.date                 ='[20150515 20170814]';

%% cell count settings
keys.cc.factors                 ={'epoch','space'};
keys.cc.conditions_to_plot      ={'Msac'};
keys.cc.plot_types              ={'space_and_epoch','per_task','visuomotor'}; 
keys.cc.epochsE.Msac        ={'INI', 'Fhol','Cue','MemL','PreS','PeriS','TIhol','Tons','Thol'}';
keys.cc.epochsS.Msac        ={'Cue','MemL','PreS','PeriS','TIhol','Tons','Thol'}';
keys.cc.choices             =[0,1];

%% for cell counts in general, but could also be used for PSTH!?
keys.tt.epoch_criterion             ='SxE or epoch only';
keys.tt.space_criterion             ='interaction or space only';

%% epochs
keys.EPOCHS_PER_TYPE{3}={...
    'INI',      2,	-0.4,	-0.1,   'INI';...
    'Fhol',     6,	-0.3,	0,      'INI';...
    'Cue',      6,	0.04,   0.14,   'INI';...
    'MemE',     7,	0.04,   0.14,   'INI';... 
    'MemL',     9,	-0.3, 	-0,     'INI';... 
    'PreS',     60,	-0.1, 	-0.01,  'INI';...
    'Pre2',     60,	-0.1, 	-0.01,  'INI';...
    'PeriS',	60,	-0.01,  0.05,   'INI';...
    'Peri2',	60,	-0.01,  0.05,   'INI';...
    'TIhol',	10,	0,      0.15,   'INI';...
    'Tons',     5,	0,      0.1,    'INI';...
    'Thol',     5,	0.2,    0.5,    'INI';...
    };

keys.WINDOWS_PER_TYPE{3}={...
    'Cue',          6,	-0.55,   0.58;...
    'Saccade',      60,	-0.6,   0.9;...
    'Target',      5,	-0.2,   0.45;...
    };

keys.ANOVAS_PER_TYPE(3).epoch={'Fhol' 'Cue';...
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
keys.ANOVAS_PER_TYPE(3).main              ={'Fhol','Cue','MemE','MemL','PreS','PeriS','Tons','Thol'}';

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
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Cue', 0, 0.78};
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
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Target', -0.2, 0.45};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.EP_EN/255; keys.colors.EP_SU/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='MemL';
keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Target';


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
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Cue', 0, 0.78};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IS_IN/255; keys.colors.NH_CS_IN/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='Fhol';
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
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Target', -0.2, 0.45};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IS_IN/255; keys.colors.NH_CS_IN/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='MemL';
keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Target';


%% population PSTH settings
cc=0;
cc=cc+1;% 3 ungrouped raw
keys.pop(cc).normalization           = 'none';
keys.pop(cc).group_parameter         = 'ungrouped';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
cc=cc+1;% 3 ungrouped raw
keys.pop(cc).normalization           = 'none';
keys.pop(cc).group_parameter         = 'ungrouped';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0 1];
keys.pop(cc).hands                   = [0];
cc=cc+1;% 3 ungrouped raw
keys.pop(cc).normalization           = 'by_effector';
keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
keys.pop(cc).group_parameter         = 'ungrouped';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];


%% categories
cc=cc+1;% visual response?
keys.pop(cc).group_parameter         = 'visual_Msac_opt';
keys.pop(cc).group_excluded          ={'0'};
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_effector';
keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
keys.pop(cc).FR_subtract_baseline    = 0;
cc=cc+1;% motor output only?
keys.pop(cc).group_parameter         = 'motor_Msac_opt';
keys.pop(cc).group_excluded          ={'0'};
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_effector';
keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
keys.pop(cc).FR_subtract_baseline    = 0;
cc=cc+1;% visuomotor cells???
keys.pop(cc).group_parameter         = 'visuomotor_Msac_opt';
keys.pop(cc).group_excluded          ={'0'};
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_effector';
keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
keys.pop(cc).FR_subtract_baseline    = 0;
cc=cc+1;% visual response?
keys.pop(cc).group_parameter         = 'visual_Msac_opt';
keys.pop(cc).group_excluded          ={'0'};
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'none';
keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
keys.pop(cc).FR_subtract_baseline    = 0;
cc=cc+1;% motor output only?
keys.pop(cc).group_parameter         = 'motor_Msac_opt';
keys.pop(cc).group_excluded          ={'0'};
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'none';
keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
keys.pop(cc).FR_subtract_baseline    = 0;
cc=cc+1;% visuomotor cells???
keys.pop(cc).group_parameter         = 'visuomotor_Msac_opt';
keys.pop(cc).group_excluded          ={'0'};
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'none';
keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
keys.pop(cc).FR_subtract_baseline    = 0;

%% space tuning in categories
cc=cc+1;% visual cells???
keys.pop(cc).combine_tuning_properties  ={'comb_cue','in_NH_Cue_epoch_Msac_opt','in_Cue_spaceLR_Msac_opt'};
keys.pop(cc).group_parameter         = 'comb_cue';
keys.pop(cc).group_excluded          = {''};
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_effector';
keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
keys.pop(cc).FR_subtract_baseline    = 0;
keys.pop(cc).Selection                  ={'visual_Msac_opt','1'};
cc=cc+1;% visuomotor cells???
keys.pop(cc).combine_tuning_properties  ={'comb_cue','in_NH_Cue_epoch_Msac_opt','in_Cue_spaceLR_Msac_opt'};
keys.pop(cc).group_parameter         = 'comb_cue';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_effector';
keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
keys.pop(cc).FR_subtract_baseline    = 0;
keys.pop(cc).Selection                  ={'visuomotor_Msac_opt','1'};
cc=cc+1;% visuomotor cells???
keys.pop(cc).combine_tuning_properties  ={'comb_sac','in_NH_TIhol_epoch_Msac_opt','in_TIhol_spaceLR_Msac_opt'};
keys.pop(cc).group_parameter         = 'comb_sac';
keys.pop(cc).group_excluded          = {''};
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_effector';
keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
keys.pop(cc).FR_subtract_baseline    = 0;
keys.pop(cc).Selection                  ={'visuomotor_Msac_opt','1'};
cc=cc+1;% motor cells???
keys.pop(cc).combine_tuning_properties  ={'comb_sac','in_NH_TIhol_epoch_Msac_opt','in_TIhol_spaceLR_Msac_opt'};
keys.pop(cc).group_parameter         = 'comb_sac';
keys.pop(cc).group_excluded          = {''};
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_effector';
keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
keys.pop(cc).FR_subtract_baseline    = 0;
keys.pop(cc).Selection                  ={'motor_Msac_opt',1};

%% enhancement/suppression in categories
cc=cc+1;% visual cells???
keys.pop(cc).group_parameter         = 'in_NH_Cue_epoch_Msac_opt';
keys.pop(cc).group_excluded          = {''};
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_effector';
keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
keys.pop(cc).FR_subtract_baseline    = 0;
keys.pop(cc).Selection                  ={'visual_Msac_opt','1'};
cc=cc+1;% visuomotor cells???
keys.pop(cc).group_parameter         = 'comb_cue_sac';
keys.pop(cc).combine_tuning_properties  ={'comb_cue_sac','in_NH_Cue_epoch_Msac_opt','in_NH_TIhol_epoch_Msac_opt'};
keys.pop(cc).group_excluded          = {''};
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_effector';
keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
keys.pop(cc).FR_subtract_baseline    = 0;
keys.pop(cc).Selection                  ={'visuomotor_Msac_opt','1'};
keys.pop(cc).unselect                  ={'comb','enbi';'comb','bien';'comb','subi';'comb','bisu'};
cc=cc+1;% motor cells???
keys.pop(cc).group_parameter         = 'in_NH_TIhol_epoch_Msac_opt';
keys.pop(cc).group_excluded          = {''};
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_effector';
keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
keys.pop(cc).FR_subtract_baseline    = 0;
keys.pop(cc).Selection                  ={'motor_Msac_opt','1'};

cc=cc+1;% visual cells???
keys.pop(cc).group_parameter         = 'in_NH_Cue_epoch_Msac_opt';
keys.pop(cc).group_excluded          = {''};
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'none';
keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
keys.pop(cc).FR_subtract_baseline    = 0;
keys.pop(cc).Selection                  ={'visual_Msac_opt','1'};
cc=cc+1;% visuomotor cells???
keys.pop(cc).group_parameter         = 'comb';
keys.pop(cc).combine_tuning_properties  ={'comb','in_NH_Cue_epoch_Msac_opt','in_NH_TIhol_epoch_Msac_opt'};
keys.pop(cc).group_excluded          = {''};
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'none';
keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
keys.pop(cc).FR_subtract_baseline    = 0;
keys.pop(cc).Selection                  ={'visuomotor_Msac_opt','1'};
keys.pop(cc).unselect                  ={'comb','enbi';'comb','bien';'comb','subi';'comb','bisu'};
cc=cc+1;% motor cells???
keys.pop(cc).group_parameter         = 'in_NH_TIhol_epoch_Msac_opt';
keys.pop(cc).group_excluded          = {''};
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'none';
keys.pop(cc).epoch_for_normalization = 'Fhol'; %%!!!
keys.pop(cc).FR_subtract_baseline    = 0;
keys.pop(cc).Selection                  ={'motor_Msac_opt','1'};

%% choice and memory activity
cc=cc+1;% ramping cells???
keys.pop(cc).group_parameter         = 'in_MemL_spaceLR_Msac_opt';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_effector';
keys.pop(cc).epoch_for_normalization = 'Fhol';
keys.pop(cc).FR_subtract_baseline    = 0;
cc=cc+1;% ramping cells???
keys.pop(cc).group_parameter         = 'in_NH_MemL_epoch_Msac_opt';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_effector';
keys.pop(cc).epoch_for_normalization = 'Fhol';
keys.pop(cc).FR_subtract_baseline    = 0;


cc=cc+1;% ramping cells???
keys.pop(cc).group_parameter         = 'in_MemL_spaceLR_Msac_opt';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0 1];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_effector';
keys.pop(cc).epoch_for_normalization = 'Fhol';
keys.pop(cc).FR_subtract_baseline    = 0;
cc=cc+1;% ramping cells???
keys.pop(cc).group_parameter         = 'in_NH_MemL_epoch_Msac_opt';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0 1];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_effector';
keys.pop(cc).epoch_for_normalization = 'Fhol';
keys.pop(cc).FR_subtract_baseline    = 0;

cc=cc+1;% gaze cells???
keys.pop(cc).group_parameter         = 'in_Thol_spaceLR_Msac_opt';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_effector';
keys.pop(cc).epoch_for_normalization = 'Fhol';
keys.pop(cc).FR_subtract_baseline    = 0;
cc=cc+1;% gaze cells???
keys.pop(cc).group_parameter         = 'in_NH_Thol_epoch_Msac_opt';
keys.pop(cc).conditions_to_plot      = {'Msac'};
keys.pop(cc).choices                 = [0];
keys.pop(cc).hands                   = [0];
keys.pop(cc).normalization           = 'by_effector';
keys.pop(cc).epoch_for_normalization = 'Fhol';
keys.pop(cc).FR_subtract_baseline    = 0;



