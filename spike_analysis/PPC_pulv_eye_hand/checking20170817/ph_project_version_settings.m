keys.project_versions={''};
keys.project_version='checking20170817';
keys.filelist_formatted={};

%% what to plot
keys.plot.single_cells =1;
keys.plot.waveforms=0;
%% to check carefully
%keys.position_and_plotting_arrangements             ={'hands'};
keys.position_and_plotting_arrangements             ={'hands_inactivation'};

%% computation settings
keys.cal.datasets                   =[41];
keys.cal.effectors                  =[3,4,6];
keys.cal.reach_hand                 =[1,2];
keys.cal.types                      =[4];
keys.cal.units_from_sorting_table   =1;

%% batching
keys.batching.combine_monkeys           =1;
keys.batching.targets                   ={'MIP_L','MIP_R'};
keys.batching.monkeys                   ={'Linus'};
% keys.Tesla.date                         ='[20160217 20180101]';
%keys.Linus.date                         ='[20161103 20180101]';
keys.Linus.date                         ='[20170714 20170817]';
keys.cal.perturbation_groups            ={0,[2,3,4,5,6,7,8]}; 

%% cell count settings
keys.cc.plot_types                      ={'hands_inactivation'};
keys.cc.factors                         ={'epoch'};%{'epoch','space','hand'};
keys.cc.conditions_to_plot              ={'Dcfr','Ddre','Ddsa'};
keys.plot.polars_on_extra_figure        =0;

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
    'Reach',        62,	-0.35,  0.7;...
    };  




%% tuning table selection
% Most important parameter to understand how cell seleciton works in the population analysis, is the trial criterion per condition:
% units with less than this amount of trials for one of the conditions will be excluded. 
keys.cal.min_trials_per_condition       =5;
% This criterion is applied every time the tuning table is loaded via ph_load_extended_tuning_table and ph_reduce_tuning_table. 
% The tricky part here is to define the conditions in which we ask for a specific number of trials. 
% Each condition is defined by tasktype(contains type,effector,and
% position_arrangement),position (depends on arrangement),hand used and choice/instruced trial (MISSING HERE: PERTURBATION)
keys.tt.tasktypes                   ={'Ddre_han'}; % typically only one tasktype defines
keys.tt.hands                       =[1 2];
keys.tt.choices                     =[0 1]; %IMPORTANT: and also not really perfect, for choice trials trial criterion is applied by hemifield, not by position.
% Each unique combination of the above parameters has to contain at least keys.cal.min_trials_per_condition trials, if not the cell is excluded in ph_reduce_tuning_table
keys.tt.selection                   ={'target','MIP'};                         % easy to use if there is a parameter in the tuning table for which you want your cells to have the same value
%                                       'in_NH_TIhol_position_Msac_opt','true'};  % each row in the cell arryáy will be used to exclude cells that don't have the specifie characteristic
keys.tt.unselect                    ={}; % easy to use if there is a parameter in the tuning table for which your cells shouldn't have a specific value
keys.tt.combine_tuning_properties   ={'hand_tuning','in_IH_Facq_epoch_Ddre_han','in_CH_Facq_epoch_Ddre_han'}; % results in enen/ensu/suen/en-/-en/su-/-su/--
% ph_load_extended_tuning_table will create an additional column from combining existing columns: 
% 1st entry is the name of the new column (to refer to it), and the following entries specify the columns which should be combined.
% this additional column can be used for grouping, specific selection and/or unselection ph_load_extended_tuning_table also creates new columns which can be used
% for multicomparison correction (sort of) by excluding cells that did not show certain combination of ANOVA effects (this is ONLY used in cell_counts so far
keys.tt.epoch_criterion             ='none'; % only relevant for cell counts
keys.tt.space_criterion             ='none';
keys.tt.hands_criterion             ='none';
keys.tt.SXH_criterion               ='none';


%% tuning onset keys

cc=0;

cc=cc+1;% 'Msac choice tuning - choice vs instructed';
ce=0;
keys.ons(cc).comparisons_title       = 'Ddre perturbation effects';
keys.ons(cc).group_parameter         = 'ungrouped';
keys.ons(cc).conditions_to_plot      = {'Ddre'};
keys.ons(cc).hands=[1 2]; %for cell exclusion
keys.ons(cc).choices=[0]; %for cell exclusion
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=1;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[-1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
keys.ons(cc).comparisons_per_effector(ce).choice{2}=0;
keys.ons(cc).comparisons_per_effector(ce).perturbation{1}=0;
keys.ons(cc).comparisons_per_effector(ce).perturbation{2}=1;
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Fixation', -1,  0.17};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.IH_IS_IN/255; keys.colors.IH_IS_CH/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
keys.ons(cc).comparisons_per_effector(ce).title='ipsi hand Ipsilateral';
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=1;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
keys.ons(cc).comparisons_per_effector(ce).choice{2}=0;
keys.ons(cc).comparisons_per_effector(ce).perturbation{1}=0;
keys.ons(cc).comparisons_per_effector(ce).perturbation{2}=1;
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Fixation', -1,  0.17};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.IH_CS_IN/255; keys.colors.IH_CS_CH/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
keys.ons(cc).comparisons_per_effector(ce).title='Ipsi hand Contralateral';
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=2;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=2;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[-1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
keys.ons(cc).comparisons_per_effector(ce).choice{2}=0;
keys.ons(cc).comparisons_per_effector(ce).perturbation{1}=0;
keys.ons(cc).comparisons_per_effector(ce).perturbation{2}=1;
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Fixation', -1,  0.17};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.CH_IS_IN/255; keys.colors.CH_IS_CH/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
keys.ons(cc).comparisons_per_effector(ce).title='Contra hand Ipsilateral';
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=2;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=2;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
keys.ons(cc).comparisons_per_effector(ce).choice{2}=0;
keys.ons(cc).comparisons_per_effector(ce).perturbation{1}=0;
keys.ons(cc).comparisons_per_effector(ce).perturbation{2}=1;
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Fixation', -1,  0.17};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.CH_CS_IN/255; keys.colors.CH_CS_CH/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';
keys.ons(cc).comparisons_per_effector(ce).title='Contra hand Contralateral';


%% population PSTH settings
cc=0;
% 1
cc=cc+1;
keys.pop(cc).combine_tuning_properties  ={'hand_tuning','in_IH_Facq_epoch_Ddre_han','in_CH_Facq_epoch_Ddre_han'};
keys.pop(cc).group_parameter            = 'ungrouped'; %hand_tuning
keys.pop(cc).normalization              = 'none';
keys.pop(cc).epoch_for_normalization    = 'Cue'; %%!!!
keys.pop(cc).conditions_to_plot         = {'Ddre'}; 
keys.pop(cc).FR_subtract_baseline       = 0;
keys.pop(cc).hands                      = [1 2];
keys.pop(cc).perturbation               = [0 1];
keys.pop(cc).choices                    = 0;
keys.pop(cc).Selection                  ={};
keys.pop(cc).unselect                   ={};
keys.pop(cc).group_excluded             = {};%{'susu','ensu','suen','-su','su-','--'};


%% scatter keys
cs=0;
%% reaches inactivation stuff

% example for plotting FR control versus perturbation
cs=cs+1;
keys.sct(cs).X='in_IH_IS_Fhol_CT_FR_Ddre_han';
keys.sct(cs).Y='in_IH_IS_Fhol_PT_FR_Ddre_han';
keys.sct(cs).X_sig='in_IH_IS_Fhol_PT_Ddre_han'; % O > ^ o (filled) 
keys.sct(cs).Y_sig='in_IH_IS_Fhol_PT_Ddre_han';
keys.sct(cs).tt.tasktypes={'Ddre_han'};
keys.sct(cs).tt.IC_for_criterion='in';
keys.sct(cs).tt.hands=[1 2];
keys.sct(cs).tt.choices=[0];


cs=cs+1;
keys.sct(cs).X='in_IH_CS_Fhol_CT_FR_Ddre_han';
keys.sct(cs).Y='in_IH_CS_Fhol_PT_FR_Ddre_han';
keys.sct(cs).X_sig='in_IH_CS_Fhol_PT_Ddre_han';
keys.sct(cs).Y_sig='in_IH_CS_Fhol_PT_Ddre_han';
keys.sct(cs).tt.tasktypes={'Ddre_han'};
keys.sct(cs).tt.IC_for_criterion='in';
keys.sct(cs).tt.hands=[1 2];
keys.sct(cs).tt.choices=[0];

cs=cs+1;
keys.sct(cs).X='in_CH_CS_Fhol_CT_FR_Ddre_han';
keys.sct(cs).Y='in_CH_CS_Fhol_PT_FR_Ddre_han';
keys.sct(cs).X_sig='in_CH_CS_Fhol_PT_Ddre_han';
keys.sct(cs).Y_sig='in_CH_CS_Fhol_PT_Ddre_han';
keys.sct(cs).tt.tasktypes={'Ddre_han'};
keys.sct(cs).tt.IC_for_criterion='in';
keys.sct(cs).tt.hands=[1 2];
keys.sct(cs).tt.choices=[0];

cs=cs+1;
keys.sct(cs).X='in_CH_IS_Fhol_CT_FR_Ddre_han';
keys.sct(cs).Y='in_CH_IS_Fhol_PT_FR_Ddre_han';
keys.sct(cs).X_sig='in_CH_IS_Fhol_PT_Ddre_han';
keys.sct(cs).Y_sig='in_CH_IS_Fhol_PT_Ddre_han';
keys.sct(cs).tt.tasktypes={'Ddre_han'};
keys.sct(cs).tt.IC_for_criterion='in';
keys.sct(cs).tt.hands=[1 2];
keys.sct(cs).tt.choices=[0];

keys.sct(cs).X='in_IH_IS_Cue_CT_FR_Ddre_han';
keys.sct(cs).Y='in_IH_IS_Cue_PT_FR_Ddre_han';
keys.sct(cs).X_sig='in_IH_IS_Cue_PT_Ddre_han'; % O > ^ o (filled) 
keys.sct(cs).Y_sig='in_IH_IS_Cue_PT_Ddre_han';
keys.sct(cs).tt.tasktypes={'Ddre_han'};
keys.sct(cs).tt.IC_for_criterion='in';
keys.sct(cs).tt.hands=[1 2];
keys.sct(cs).tt.choices=[0];


cs=cs+1;
keys.sct(cs).X='in_IH_CS_Cue_CT_FR_Ddre_han';
keys.sct(cs).Y='in_IH_CS_Cue_PT_FR_Ddre_han';
keys.sct(cs).X_sig='in_IH_CS_Cue_PT_Ddre_han';
keys.sct(cs).Y_sig='in_IH_CS_Cue_PT_Ddre_han';
keys.sct(cs).tt.tasktypes={'Ddre_han'};
keys.sct(cs).tt.IC_for_criterion='in';
keys.sct(cs).tt.hands=[1 2];
keys.sct(cs).tt.choices=[0];

cs=cs+1;
keys.sct(cs).X='in_CH_CS_Cue_CT_FR_Ddre_han';
keys.sct(cs).Y='in_CH_CS_Cue_PT_FR_Ddre_han';
keys.sct(cs).X_sig='in_CH_CS_Cue_PT_Ddre_han';
keys.sct(cs).Y_sig='in_CH_CS_Cue_PT_Ddre_han';
keys.sct(cs).tt.tasktypes={'Ddre_han'};
keys.sct(cs).tt.IC_for_criterion='in';
keys.sct(cs).tt.hands=[1 2];
keys.sct(cs).tt.choices=[0];

cs=cs+1;
keys.sct(cs).X='in_CH_IS_Cue_CT_FR_Ddre_han';
keys.sct(cs).Y='in_CH_IS_Cue_PT_FR_Ddre_han';
keys.sct(cs).X_sig='in_CH_IS_Cue_PT_Ddre_han';
keys.sct(cs).Y_sig='in_CH_IS_Cue_PT_Ddre_han';
keys.sct(cs).tt.tasktypes={'Ddre_han'};
keys.sct(cs).tt.IC_for_criterion='in';
keys.sct(cs).tt.hands=[1 2];
keys.sct(cs).tt.choices=[0];

keys.sct(cs).X='in_IH_IS_PreR_CT_FR_Ddre_han';
keys.sct(cs).Y='in_IH_IS_PreR_PT_FR_Ddre_han';
keys.sct(cs).X_sig='in_IH_IS_PreR_PT_Ddre_han'; % O > ^ o (filled) 
keys.sct(cs).Y_sig='in_IH_IS_PreR_PT_Ddre_han';
keys.sct(cs).tt.tasktypes={'Ddre_han'};
keys.sct(cs).tt.IC_for_criterion='in';
keys.sct(cs).tt.hands=[1 2];
keys.sct(cs).tt.choices=[0];


cs=cs+1;
keys.sct(cs).X='in_IH_CS_PreR_CT_FR_Ddre_han';
keys.sct(cs).Y='in_IH_CS_PreR_PT_FR_Ddre_han';
keys.sct(cs).X_sig='in_IH_CS_PreR_PT_Ddre_han';
keys.sct(cs).Y_sig='in_IH_CS_PreR_PT_Ddre_han';
keys.sct(cs).tt.tasktypes={'Ddre_han'};
keys.sct(cs).tt.IC_for_criterion='in';
keys.sct(cs).tt.hands=[1 2];
keys.sct(cs).tt.choices=[0];

cs=cs+1;
keys.sct(cs).X='in_CH_CS_PreR_CT_FR_Ddre_han';
keys.sct(cs).Y='in_CH_CS_PreR_PT_FR_Ddre_han';
keys.sct(cs).X_sig='in_CH_CS_PreR_PT_Ddre_han';
keys.sct(cs).Y_sig='in_CH_CS_PreR_PT_Ddre_han';
keys.sct(cs).tt.tasktypes={'Ddre_han'};
keys.sct(cs).tt.IC_for_criterion='in';
keys.sct(cs).tt.hands=[1 2];
keys.sct(cs).tt.choices=[0];

cs=cs+1;
keys.sct(cs).X='in_CH_IS_PreR_CT_FR_Ddre_han';
keys.sct(cs).Y='in_CH_IS_PreR_PT_FR_Ddre_han';
keys.sct(cs).X_sig='in_CH_IS_PreR_PT_Ddre_han';
keys.sct(cs).Y_sig='in_CH_IS_PreR_PT_Ddre_han';
keys.sct(cs).tt.tasktypes={'Ddre_han'};
keys.sct(cs).tt.IC_for_criterion='in';
keys.sct(cs).tt.hands=[1 2];
keys.sct(cs).tt.choices=[0];

% cs=cs+1;
% keys.sct(cs).X='in_IH_IS_Cue_PT_ES_Ddre_han';
% keys.sct(cs).Y='in_IH_CS_Cue_PT_ES_Ddre_han';
% keys.sct(cs).X_sig='in_IH_IS_Cue_PT_Ddre_han';
% keys.sct(cs).Y_sig='in_IH_CS_Cue_PT_Ddre_han';
% keys.sct(cs).tt.tasktypes={'Ddre_han'};
% keys.sct(cs).tt.IC_for_criterion='in';
% keys.sct(cs).tt.hands=[1 2];
% keys.sct(cs).tt.choices=[0];

% cs=cs+1;
% keys.sct(cs).X='in_CH_IS_Cue_PT_ES_Ddre_han';
% keys.sct(cs).Y='in_CH_CS_Cue_PT_ES_Ddre_han';
% keys.sct(cs).X_sig='in_CH_IS_Cue_PT_Ddre_han';
% keys.sct(cs).Y_sig='in_CH_CS_Cue_PT_Ddre_han';
% keys.sct(cs).tt.tasktypes={'Ddre_han'};
% keys.sct(cs).tt.IC_for_criterion='in';
% keys.sct(cs).tt.hands=[1 2];
% keys.sct(cs).tt.choices=[0];
% 
% cs=cs+1;
% keys.sct(cs).X='in_IH_IS_Del_PT_ES_Ddre_han';
% keys.sct(cs).Y='in_IH_CS_Del_PT_ES_Ddre_han';
% keys.sct(cs).X_sig='in_IH_IS_Del_PT_Ddre_han';
% keys.sct(cs).Y_sig='in_IH_CS_Del_PT_Ddre_han';
% keys.sct(cs).tt.tasktypes={'Ddre_han'};
% keys.sct(cs).tt.IC_for_criterion='in';
% keys.sct(cs).tt.hands=[1 2];
% keys.sct(cs).tt.choices=[0];
% 
% cs=cs+1;
% keys.sct(cs).X='in_CH_IS_Del_PT_ES_Ddre_han';
% keys.sct(cs).Y='in_CH_CS_Del_PT_ES_Ddre_han';
% keys.sct(cs).X_sig='in_CH_IS_Del_PT_Ddre_han';
% keys.sct(cs).Y_sig='in_CH_CS_Del_PT_Ddre_han';
% keys.sct(cs).tt.tasktypes={'Ddre_han'};
% keys.sct(cs).tt.IC_for_criterion='in';
% keys.sct(cs).tt.hands=[1 2];
% keys.sct(cs).tt.choices=[0];
% 
% cs=cs+1;
% keys.sct(cs).X='in_IH_IS_PreR_PT_ES_Ddre_han';
% keys.sct(cs).Y='in_IH_CS_PreR_PT_ES_Ddre_han';
% keys.sct(cs).X_sig='in_IH_IS_PreR_PT_Ddre_han';
% keys.sct(cs).Y_sig='in_IH_CS_PreR_PT_Ddre_han';
% keys.sct(cs).tt.tasktypes={'Ddre_han'};
% keys.sct(cs).tt.IC_for_criterion='in';
% keys.sct(cs).tt.hands=[1 2];
% keys.sct(cs).tt.choices=[0];
% 
% cs=cs+1;
% keys.sct(cs).X='in_CH_IS_PreR_PT_ES_Ddre_han';
% keys.sct(cs).Y='in_CH_CS_PreR_PT_ES_Ddre_han';
% keys.sct(cs).X_sig='in_CH_IS_PreR_PT_Ddre_han';
% keys.sct(cs).Y_sig='in_CH_CS_PreR_PT_Ddre_han';
% keys.sct(cs).tt.tasktypes={'Ddre_han'};
% keys.sct(cs).tt.IC_for_criterion='in';
% keys.sct(cs).tt.hands=[1 2];
% keys.sct(cs).tt.choices=[0];
% 
% cs=cs+1;
% keys.sct(cs).X='in_IH_IS_PeriR_PT_ES_Ddre_han';
% keys.sct(cs).Y='in_IH_CS_PeriR_PT_ES_Ddre_han';
% keys.sct(cs).X_sig='in_IH_IS_PeriR_PT_Ddre_han';
% keys.sct(cs).Y_sig='in_IH_CS_PeriR_PT_Ddre_han';
% keys.sct(cs).tt.tasktypes={'Ddre_han'};
% keys.sct(cs).tt.IC_for_criterion='in';
% keys.sct(cs).tt.hands=[1 2];
% keys.sct(cs).tt.choices=[0];
% 
% cs=cs+1;
% keys.sct(cs).X='in_CH_IS_PeriR_PT_ES_Ddre_han';
% keys.sct(cs).Y='in_CH_CS_PeriR_PT_ES_Ddre_han';
% keys.sct(cs).X_sig='in_CH_IS_PeriR_PT_Ddre_han';
% keys.sct(cs).Y_sig='in_CH_CS_PeriR_PT_Ddre_han';
% keys.sct(cs).tt.tasktypes={'Ddre_han'};
% keys.sct(cs).tt.IC_for_criterion='in';
% keys.sct(cs).tt.hands=[1 2];
% keys.sct(cs).tt.choices=[0];

%% saccades inactivation stuff
% cs=cs+1;
% keys.sct(cs).X='in_IH_IS_Fhol_PT_ES_Ddsa_han';
% keys.sct(cs).Y='in_IH_CS_Fhol_PT_ES_Ddsa_han';
% keys.sct(cs).X_sig='in_IH_IS_Fhol_PT_Ddsa_han';
% keys.sct(cs).Y_sig='in_IH_CS_Fhol_PT_Ddsa_han';
% keys.sct(cs).tt.tasktypes={'Ddre_han'};
% keys.sct(cs).tt.IC_for_criterion='in';
% keys.sct(cs).tt.hands=[1 2];
% keys.sct(cs).tt.choices=[0];
% 
% cs=cs+1;
% keys.sct(cs).X='in_CH_IS_Fhol_PT_ES_Ddsa_han';
% keys.sct(cs).Y='in_CH_CS_Fhol_PT_ES_Ddsa_han';
% keys.sct(cs).X_sig='in_CH_IS_Fhol_PT_Ddsa_han';
% keys.sct(cs).Y_sig='in_CH_CS_Fhol_PT_Ddsa_han';
% keys.sct(cs).tt.tasktypes={'Ddre_han'};
% keys.sct(cs).tt.IC_for_criterion='in';
% keys.sct(cs).tt.hands=[1 2];
% keys.sct(cs).tt.choices=[0];
% 
% cs=cs+1;
% keys.sct(cs).X='in_IH_IS_Cue_PT_ES_Ddsa_han';
% keys.sct(cs).Y='in_IH_CS_Cue_PT_ES_Ddsa_han';
% keys.sct(cs).X_sig='in_IH_IS_Cue_PT_Ddsa_han';
% keys.sct(cs).Y_sig='in_IH_CS_Cue_PT_Ddsa_han';
% keys.sct(cs).tt.tasktypes={'Ddre_han'};
% keys.sct(cs).tt.IC_for_criterion='in';
% keys.sct(cs).tt.hands=[1 2];
% keys.sct(cs).tt.choices=[0];
% 
% cs=cs+1;
% keys.sct(cs).X='in_CH_IS_Cue_PT_ES_Ddsa_han';
% keys.sct(cs).Y='in_CH_CS_Cue_PT_ES_Ddsa_han';
% keys.sct(cs).X_sig='in_CH_IS_Cue_PT_Ddsa_han';
% keys.sct(cs).Y_sig='in_CH_CS_Cue_PT_Ddsa_han';
% keys.sct(cs).tt.tasktypes={'Ddre_han'};
% keys.sct(cs).tt.IC_for_criterion='in';
% keys.sct(cs).tt.hands=[1 2];
% keys.sct(cs).tt.choices=[0];
% 
% cs=cs+1;
% keys.sct(cs).X='in_IH_IS_Del_PT_ES_Ddsa_han';
% keys.sct(cs).Y='in_IH_CS_Del_PT_ES_Ddsa_han';
% keys.sct(cs).X_sig='in_IH_IS_Del_PT_Ddsa_han';
% keys.sct(cs).Y_sig='in_IH_CS_Del_PT_Ddsa_han';
% keys.sct(cs).tt.tasktypes={'Ddre_han'};
% keys.sct(cs).tt.IC_for_criterion='in';
% keys.sct(cs).tt.hands=[1 2];
% keys.sct(cs).tt.choices=[0];
% 
% cs=cs+1;
% keys.sct(cs).X='in_CH_IS_Del_PT_ES_Ddsa_han';
% keys.sct(cs).Y='in_CH_CS_Del_PT_ES_Ddsa_han';
% keys.sct(cs).X_sig='in_CH_IS_Del_PT_Ddsa_han';
% keys.sct(cs).Y_sig='in_CH_CS_Del_PT_Ddsa_han';
% keys.sct(cs).tt.tasktypes={'Ddre_han'};
% keys.sct(cs).tt.IC_for_criterion='in';
% keys.sct(cs).tt.hands=[1 2];
% keys.sct(cs).tt.choices=[0];
% 
% cs=cs+1;
% keys.sct(cs).X='in_PreS_IS_IH_PT_ES_Ddsa_han';
% keys.sct(cs).Y='in_PreS_CS_IH_PT_ES_Ddsa_han';
% keys.sct(cs).X_sig='in_PreS_IS_IH_PT_Ddsa_han';
% keys.sct(cs).Y_sig='in_PreS_CS_IH_PT_Ddsa_han';
% keys.sct(cs).tt.tasktypes={'Ddre_han'};
% keys.sct(cs).tt.IC_for_criterion='in';
% keys.sct(cs).tt.hands=[1 2];
% keys.sct(cs).tt.choices=[0];
% 
% cs=cs+1;
% keys.sct(cs).X='in_PreS_IS_CH_PT_ES_Ddsa_han';
% keys.sct(cs).Y='in_PreS_CS_CH_PT_ES_Ddsa_han';
% keys.sct(cs).X_sig='in_PreS_IS_CH_PT_Ddsa_han';
% keys.sct(cs).Y_sig='in_PreS_CS_CH_PT_Ddsa_han';
% keys.sct(cs).tt.tasktypes={'Ddre_han'};
% keys.sct(cs).tt.IC_for_criterion='in';
% keys.sct(cs).tt.hands=[1 2];
% keys.sct(cs).tt.choices=[0];
% 
% cs=cs+1;
% keys.sct(cs).X='in_PeriS_IS_IH_PT_ES_Ddsa_han';
% keys.sct(cs).Y='in_PeriS_CS_IH_PT_ES_Ddsa_han';
% keys.sct(cs).X_sig='in_PeriS_IS_IH_PT_Ddsa_han';
% keys.sct(cs).Y_sig='in_PeriS_CS_IH_PT_Ddsa_han';
% keys.sct(cs).tt.tasktypes={'Ddre_han'};
% keys.sct(cs).tt.IC_for_criterion='in';
% keys.sct(cs).tt.hands=[1 2];
% keys.sct(cs).tt.choices=[0];
% 
% cs=cs+1;
% keys.sct(cs).X='in_PeriS_IS_CH_PT_ES_Ddsa_han';
% keys.sct(cs).Y='in_PeriS_CS_CH_PT_ES_Ddsa_han';
% keys.sct(cs).X_sig='in_PeriS_IS_CH_PT_Ddsa_han';
% keys.sct(cs).Y_sig='in_PeriS_CS_CH_PT_Ddsa_han';
% keys.sct(cs).tt.tasktypes={'Ddre_han'};
% keys.sct(cs).tt.IC_for_criterion='in';
% keys.sct(cs).tt.hands=[1 2];
% keys.sct(cs).tt.choices=[0];
% 

