keys.pdf_folders={''};
keys.pdf_folder='testtesttest';
keys.filelist_formatted={};

keys.effectors_on_same_figure   =1;
keys.extend_population          =1;
keys.plot_FR_separately         =1;
keys.contra_color               ='orange';

%% to check carefully
keys.combine_monkeys            =0;
keys.case_summaries             ={'hands'};
keys.FR_subtract_baseline       =0;

%% computation settings
keys.epochs_to_plot_PSTH{4}                 ={'Facq','Cue','PeriS','PeriR'}';

keys.effectors                  =[0,1,2,3,4,6];
keys.reach_hand                 =[1,2];
keys.types                      =[4];
keys.targets                    ={'MIP'};
keys.monkeys                    ={'Tesla'};
keys.Tesla.date                 ='[20160217 20170101]';


%% cell count settings
keys.cell_count_factors         ={'epoch','space','hand'};
keys.cc.conditions_to_plot      ={'Dcfr','Ddre','Ddsa'};
keys.cc.only_both_hands         =1;


%% epochs
keys.EPOCHS_PER_TYPE{4}={...
    'INI',      2,	-0.5,	0,      -0.3,	-0,     'INI';...
    'Facq',     3,	-0.6,   0,      -0.35,	0,      'INI';...
    'Fhol',     6,	-0.6,	0,      -0.45,	-0.05,  'INI';...
    'Cue',      6,	0,      0.4, 	0.05,   0.35,   'INI';...
    'Del',      4, 	-1,     0,      -0.35,  -0.05,  'INI';...
    'PreS',     60,	-0.22,  0,      -0.22,  -0.02,  'INI';...
    'PeriS',	60,	-0.1,   0.05,	-0.03,  0.05,   'INI';...
    'PostS',	61,	0,      0.22,	0.02,   0.22,   'INI';...
    'PreR',     62,	-0.22,  0,      -0.22, 	-0.02,  'INI';...
    'PeriR',	62,	-0.05,  0.2,	-0.03, 	0.15,   'INI';...
    'PostR',	63,	0,      0.22,	0.02,   0.22,   'INI';...
    'Thol',     5,	0.1,   0.5,     0.15,   0.45,   'INI';...
    };

%% subregions
keys.Subregions_separately=0;

keys.Subregions{1}{1}.monkey='Tes';
keys.Subregions{1}{1}.target='MIP';
keys.Subregions{1}{1}.grid_x=4;
keys.Subregions{1}{1}.grid_y=0;
keys.Subregions{1}{1}.z_min=0;
keys.Subregions{1}{1}.z_max=100;

keys.Subregions{1}{2}.monkey='Tes';
keys.Subregions{1}{2}.target='MIP';
keys.Subregions{1}{2}.grid_x=5;
keys.Subregions{1}{2}.grid_y=0;
keys.Subregions{1}{2}.z_min=0;
keys.Subregions{1}{2}.z_max=100;

keys.Subregions{1}{3}.monkey='Tes';
keys.Subregions{1}{3}.target='MIP';
keys.Subregions{1}{3}.grid_x=6;
keys.Subregions{1}{3}.grid_y=-2;
keys.Subregions{1}{3}.z_min=0;
keys.Subregions{1}{3}.z_max=100;

keys.Subregions{2}{1}.monkey='Tes';
keys.Subregions{2}{1}.target='MIP';
keys.Subregions{2}{1}.grid_x=7;
keys.Subregions{2}{1}.grid_y=-4;
keys.Subregions{2}{1}.z_min=0;
keys.Subregions{2}{1}.z_max=100;

keys.Subregions{2}{2}.monkey='Tes';
keys.Subregions{2}{2}.target='MIP';
keys.Subregions{2}{2}.grid_x=7;
keys.Subregions{2}{2}.grid_y=-8;
keys.Subregions{2}{2}.z_min=0;
keys.Subregions{2}{2}.z_max=100;

keys.n_Subregions=numel(keys.Subregions);

%% population PSTH settings
cc=0;
cc=cc+1;%  1
keys.pop(cc).normalization                     = 'none';
keys.pop(cc).group_parameter                   = 'ungrouped';
keys.pop(cc).epoch_for_normalization           = 'INI';
keys.pop(cc).conditions_to_plot                = {'Ddre','Ddsa','Dcfr'};

cc=cc+1;%  2
keys.pop(cc).normalization                     = 'none';
keys.pop(cc).group_parameter                   = 'untuned_Cue_response';
keys.pop(cc).conditions_to_plot                = {'Ddre','Ddsa','Dcfr'};
keys.pop(cc).combine_tuning_properties         ={'untuned_Cue_response','in_Cue_spaceLR_Ddre_han','in_Cue_epoch_Ddre_han'};
keys.pop(cc).group_excluded                    ={'','-','--','IS-','ISen','ISsu','ISbi','CS-','CSen','CSsu','CSbi','-bi'};

cc=cc+1;%  3
keys.pop(cc).normalization                     = 'by_effector';
keys.pop(cc).group_parameter                   = 'untuned_Cue_response';
keys.pop(cc).epoch_for_normalization           = 'Cue';
keys.pop(cc).conditions_to_plot                = {'Ddre','Ddsa','Dcfr'};
keys.pop(cc).combine_tuning_properties         ={'untuned_Cue_response','in_Cue_spaceLR_Ddre_han','in_Cue_epoch_Ddre_han'};
keys.pop(cc).group_excluded                    ={'','-','--','IS-','ISen','ISsu','ISbi','CS-','CSen','CSsu','CSbi','-bi'};

cc=cc+1;%  4
keys.pop(cc).normalization                     = 'none';
keys.pop(cc).group_parameter                   = 'all_tasks_exist';
keys.pop(cc).conditions_to_plot                = {'Ddre','Ddsa','Dcfr'};
keys.pop(cc).combine_tuning_properties         ={'all_tasks_exist','existing_Ddre_han','existing_Ddsa_han','existing_Dcfr_han'};
keys.pop(cc).group_excluded                    ={'000','100','101','110','010','011','001'};

cc=cc+1;%  5
keys.pop(cc).normalization                     = 'none';
keys.pop(cc).group_parameter                   = 'Ddre_Ddsa_exist';
keys.pop(cc).conditions_to_plot                = {'Ddre','Ddsa'};
keys.pop(cc).combine_tuning_properties         ={'Ddre_Ddsa_exist','existing_Ddre_han','existing_Ddsa_han'};
keys.pop(cc).group_excluded                    ={'00','10','01'};

cc=cc+1;%  6
keys.pop(cc).normalization                     = 'none';
keys.pop(cc).group_parameter                   = 'Ddre_Dcfr_exist';
keys.pop(cc).conditions_to_plot                = {'Ddre','Dcfr'};
keys.pop(cc).combine_tuning_properties         ={'Ddre_Dcfr_exist','existing_Ddre_han','existing_Dcfr_han'};
keys.pop(cc).group_excluded                    ={'00','10','01'};

cc=cc+1;%  7
keys.pop(cc).normalization                     = 'none';
keys.pop(cc).group_parameter                   = 'Ddsa_Dcfr_exist';
keys.pop(cc).conditions_to_plot                = {'Ddsa','Dcfr'};
keys.pop(cc).combine_tuning_properties         ={'Ddsa_Dcfr_exist','existing_Ddsa_han','existing_Dcfr_han'};
keys.pop(cc).group_excluded                    ={'00','10','01'};
