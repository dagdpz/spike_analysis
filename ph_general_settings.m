function keys=ph_general_settings(project,keys)
%% condition_identifiers

% labels for legends and colors

keys.condition_parameters={'choice','reach_hand','perturbation'};
keys.contra_ipsi_relative_to='target';

keys.labels.handsIC                 ={'AH','IH','CH'};  %% AH!??
keys.labels.perturbations           ={'','_PT'};
keys.labels.reach_hand              ={'AH','IH','CH'};
keys.labels.reach_handLR            ={'AH','LH','RH'}; 
keys.labels.choice                  ={'in','ch'};
keys.labels.perturbation            ={'','PT'};
keys.labels.stimulustype            ={'SS','TT','TD'};
keys.labels.difficulty              ={'TA','D1','D2'};
keys.labels.success                 ={'ER','SU'};
keys.labels.hemifield               ={'IS','VS','CS'};
keys.labels.fix_index               ={'IF','MF','CF'};
keys.labels.preferred               ={'NP','PF'};
keys.labels.stimuli_in_2hemifields  ={'1H','2H'};

%% labels for tuning table entries
keys.TTconditions.hands             ={'AH','LH','RH'};
keys.TTconditions.choice            ={'in','ch'};
keys.TTconditions.hemifield         ={'LS','RS'};
keys.TTlabels.UD                    ={'DN','-','UP'};
keys.TTlabels.CR                    ={'UC','-','CR'};
keys.TTlabels.choices               ={'in','-','ch'};
keys.TTlabels.true                  ={false,true}; %% this should go i think... should be a logical 0 or 1
keys.TTlabels.SglTar_Suc            ={'LS','-','RS'};
keys.TTlabels.Difficulty_Easy       ={'Ta','-','eD'}; %higher FR for T , higher FR for D
keys.TTlabels.Difficulty_Diff       ={'Ta','-','dD'}; %higher FR for T , higher FR for D
keys.TTlabels.SpatialComp_1HFTar    ={'ST','-','1T'}; %higher FR for T , higher FR for D
keys.TTlabels.SpatialComp_2HFTar    ={'ST','-','2T'}; %higher FR for single T , higher FR for double target
keys.TTlabels.epoch                 ={'su','-','en','bi'};
keys.TTlabels.hemifield             ={'LS','-','RS'}; %% call this hemifield
keys.TTlabels.hands                 ={'LH','-','RH'};
keys.TTlabels.PT                    ={'SU','-','EN'};

keys.cal.precision_fix=4;
keys.cal.precision_tar=2;
keys.cal.remove_trials_without_spikes=1;

%% general settings (multi-summary PSTH)
keys.PSTH_binwidth                      =0.01;                      % resolution of PSTH's (in seconds)
keys.gaussian_kernel                    =0.02;                      % std for the convolution to derive spie density (in seconds)
keys.kernel_type                        ='gaussian';
keys.FR_at_peak                         =0;                         % currently not used
keys.position_and_plotting_arrangements ={'hands'};                 % defines position batching and which conditions go into different figures/lines

%% Batching per figure ! subregion keys..?
keys.batching.monkeys                    ={'Curius','Linus'};       % monkeys on the project
keys.batching.combine_monkeys            =0;                        % for population analysis, monkeys can be separately or combined
keys.batching.targets                    ={'dPulv_r','dPulv_l'};    % to combine both, just put {'dPulv'}
keys.batching.Subregions_separately      =0;                        % subregions can be processed independently (if only defining target is too crude
keys.batching.Subregions{1}{1}           =struct('monkey',{''},'target',{''},'grid_x',{NaN},'grid_y',{NaN},'z_min',{NaN},'z_max',{NaN}); % subregion definitions by  z range for each grid hole
keys.batching.n_Subregions               =numel(keys.batching.Subregions);

%% criterions to exclude trials and units
keys.cal.process_spikes                  =1;      % you can choose not to run spikes at all
keys.cal.process_sites                   =1;      % you can choose not to run lfp sites at all (saving processing time)
keys.cal.process_by_block                =1;      % you can choose not to run by block (body signals f.e.) at all (saving processing time)
keys.cal.MA_selection                   ={'display',0,'keep_raw_data',1,'saccade_definition',4,'reach_1st_pos',1,'correlation_conditions',{}};                        % if you want to run MA with specific settings
keys.cal.units_from_sorting_table       =1;                         % exclude units that are not in the sorting table (and therefore apply stability/single/SNT ratings)
keys.cal.datasets                       =[];
keys.cal.completed                      =1;                         % problematic, because of where and how it is used. so far, keep it 1
keys.cal.effectors                      =[0,1,2,3,4,5,6];           % excluding trials with non-matching effectors
keys.cal.types                          =[1,2,3,4,5,6];             % excluding trials with non-matching types
keys.cal.reach_hand                     =[0,1,2];                   % excluding trials with non-matching reach_hand
keys.cal.choice                         =[0,1];                     % excluding trials with non-matching chocie
keys.cal.stablity                       =[0,1];                     % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the table
keys.cal.single_rating                  =[1,2,3];                   % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the tabl
keys.cal.automatic_stablity             =0;                         % using automatic stability assessment
keys.cal.SNR_rating                     =[1,2,3,4];                 % not assigning sorting table information if criterion is not met. Therefore only excludes when taking only units in the tabl
keys.cal.min_trials_per_condition       =5;                         % minimum trials per conditon (look at ph_arrange_positions to see how conditions are defined)
%keys.cal.min_spikes_per_unit            =50;                        % excluding units that have in total less spikes (workaround for sortcode assignment bug) - to be removed
keys.cal.perturbation_groups            ={0,[2,3,4,5,6,7,8]};       % which perturbation values from excel table will be assigned to control and perturbation for comparisons and population analysis

%% ANOVA normalization
keys.AN.normalization='none';
keys.AN.epoch_for_normalization='INI';
keys.AN.epoch_RF='INI';
keys.AN.epoch_BL='INI';
keys.AN.epoch_PF='INI';
keys.AN.epoch_GB='INI';
keys.AN.baseline_per_trial=0;
keys.AN.FR_subtract_baseline=0;
keys.AN.test_types='parametric'; %% as opposed to 'nonparametric'

%% folders and filenames
keys.filelist_as_blocks     =0;
keys.drive                  =DAG_get_server_IP;
keys.basepath_to_save       =[keys.drive 'Projects' filesep project filesep 'ephys' filesep];
spike_analysis_location     =which('ph_initiation');
keys.db_folder              =[spike_analysis_location(1:strfind(spike_analysis_location,['spike_analysis' filesep 'ph_initiation'])-1) 'Settings' filesep project filesep 'spike_analysis' filesep];
%% this folder defines where to take settings from

keys.All_monkeys={'Flaffus','Linus','Curius','Tesla','Cornelius','Magnus','TDT_brain','Bacchus'};
for m=1:numel(keys.All_monkeys)
    keys.(keys.All_monkeys{m}).sorted_neurons_filename    =[keys.All_monkeys{m}(1:3) '_sorted_neurons'];
    keys.(keys.All_monkeys{m}).sorted_neurons_sheetname    ='final_sorting';
    keys.(keys.All_monkeys{m}).filelist_formatted         ={};
end

keys.Flaffus.color    ='r';
keys.Linus.color      =[0 0 255]/255;
keys.Curius.color     =[255 0 0]/255;
keys.Tesla.color      ='y';
keys.Cornelius.color  ='m';
keys.Bacchus.color    ='b';

keys.Flaffus.marker    ='<';
keys.Linus.marker      ='s';
keys.Curius.marker     ='o';
keys.Tesla.marker      ='d';
keys.Cornelius.marker  ='v';
keys.Bacchus.marker    ='x';

%% plotting settings
keys.plot.anova_tables                  ='off';     % display anova result tables for each unit, please keep off, it will get very messy if you turn this on
keys.plot.single_trials                 =0;         % not used yet
keys.plot.single_cells                  =1;         % perform single cell plotting
keys.plot.waveforms                     =1;         % plot the waveform summary plots per session
keys.plot.polars_on_extra_figure        =0;         % recommended if there are too many conditions and the single cell heatmap plots are too small/crowded
keys.plot.eye_hand_traces               =1;         % Incredible performance booster if turned off
keys.plot.average_PSTH_line             =0;         % One PSTH line on top that represents the average of all others
%keys.plot.average_heat_maps             =0;
keys.plot.export                        =1;         % save plots as pdfs, you typically want this
keys.plot.events                        =1:100;     % select events that should be plotted (vertical lines) on all PSTH like plots
keys.plot.population_PSTH_legends       =1;         % if population legends should be plotted or not
keys.plot.cell_count_legends            =1;         % if population legends should be plotted or not
keys.plot.scatter_legends               =1;         % if population legends should be plotted or not

% ylimits
%keys.plot.FR_max_for_ylim               =50;
keys.plot.trials_max_for_ylim           =50;
keys.plot.excentricity_max_for_ylim     =30;
keys.plot.eyetrace_factor               =0.5;
keys.plot.hndtrace_factor               =0.5;

% ANOVA labels for single unit plots
keys.plot.anova_main    ={'E','in_epoch_main','','S','in_hemifield_main','','C','ch_hemifield_main','','H','in_hands_main','','ExS','in_ExS','','ExH','in_ExH','','SxH','in_SxH',''};
keys.plot.anova_effector={'E','in_epoch_main','','S','in_hemifield_main','','C','ch_hemifield_main','','H','in_hands_main','','ExS','in_ExS','','ExH','in_ExH','','SxH','in_SxH',''};
keys.plot.anova_epoch1  ={'E','in_AH','epoch','S','in','hemifield','C','ch','hemifield','H','in','hands','SxH','in','SxH'};
keys.plot.anova_epoch2  ={'LL','in_LH_LS','PT','RL','in_LH_RS','PT','LR','in_RH_LS','PT','RR','in_RH_RS','PT'};


keys.plot.PSTH_perpos_width          =0.5;
keys.plot.raster_width               =0.1;
keys.plot.PSTH_summary_width         =1;


%% colors & legends
% keys.condition_parameters={'reach_hand','choice','perturbation'};
% keys.labels.reach_hand={'NH','IH','CH'};
% keys.labels.hemifield={'IS','VS','CS'};
% keys.labels.fix_index={'IF','MF','CF'};
% keys.labels.preferred={'NP','PF'};
% keys.labels.choice={'IN','CH'};
% keys.labels.perturbation={'','PT'};
% keys.labels.stimuli_in_2hemifields={'1H','2H'};
% 
% labels.UD={'DN','-','UP'};
% labels.CR={'UC','-','CR'};
% labels.choices={'in','-','ch'};
% labels.true={'false','true'};
% labels.SglTar_Suc={'LS','-','RS'};
% labels.Difficulty_Easy={'Ta','-','eD'}; %higher FR for T , higher FR for D
% labels.Difficulty_Diff={'Ta','-','dD'}; %higher FR for T , higher FR for D
% labels.SpatialComp_1HFTar={'ST','-','1T'}; %higher FR for T , higher FR for D
% labels.SpatialComp_2HFTar={'ST','-','2T'}; %higher FR for single T , higher FR for double target
% labels.epoch={'su','-','en','bi'};
% labels.spaceLR={'LS','-','RS'};
% labels.hands={'LH','-','RH'};
% labels.PT={'SU','-','EN'};


% traces
keys.colors.eye_ver         =[0.8 0 0];
keys.colors.eye_hor         =[1 0 0];
keys.colors.rhd_ver         =[0 0.8 0];
keys.colors.rhd_hor         =[0 1 0];
keys.colors.lhd_ver         =[0 0 0.8];
keys.colors.lhd_hor         =[0 0 1];


keys.colors.AV   =[0 0 0];

keys.colors.IF   =[236 32 38];
keys.colors.MF   =[16 159 218];
keys.colors.CF   =[247 148 36];

keys.colors.IF_CS   =[236 32 38];
keys.colors.MF_CS   =[16 159 218];
keys.colors.CF_CS   =[247 148 36];

keys.colors.IF_VS   =[236 32 38]*2/3;
keys.colors.MF_VS   =[16 159 218]*2/3;
keys.colors.CF_VS   =[247 148 36]*2/3;

keys.colors.IF_IS   =[236 32 38]*1/3;
keys.colors.MF_IS   =[16 159 218]*1/3;
keys.colors.CF_IS   =[247 148 36]*1/3;

% cell count colors
keys.colors.NO_AN   =[255 255 255];
keys.colors.NO_TU   =[128 128 128];
keys.colors.EP_EN   =[255  65  0];
keys.colors.EP_BI   =[106 189  69];
keys.colors.EP_SU   =[0   65 255];
keys.colors.CR      =[255   0   0];
keys.colors.UC      =[127   0   0];


% single cell PSTH colors per position - i think those are wrong and should
% be reversed in name
keys.colors.in_LH   =[64 0 255];
keys.colors.ch_LH   =[32 0 128];
keys.colors.in_RH   =[255 255 0];
keys.colors.ch_RH   =[128 128 0];

keys.colors.in_AH   =[0 255 0];
keys.colors.ch_AH   =[0 128 0];
keys.colors.in_IH   =[64 0 255];
keys.colors.ch_IH   =[32 0 128];
keys.colors.in_CH   =[255 255 0];
keys.colors.ch_CH   =[128 128 0];

keys.colors.in_IS   =[0 255 255];
keys.colors.ch_IS   =[0 128 128];
keys.colors.in_CS   =[255 0 64];
keys.colors.ch_CS   =[128 0 32];

% 
% % single cell PSTH colors per hemifield
% keys.colors.NH_RS_IN=[255 0 64];
% keys.colors.NH_RS_CH=[128 0 32];
% keys.colors.LH_RS_IN=[255 0 255];
% keys.colors.LH_RS_CH=[128 0 128];
% keys.colors.RH_RS_IN=[255 128 0];
% keys.colors.RH_RS_CH=[128 64 0];
% keys.colors.NH_LS_IN=[0 255 255];
% keys.colors.NH_LS_CH=[0 128 128];
% keys.colors.LH_LS_IN=[0 128 255];
% keys.colors.LH_LS_CH=[0 64 128];
% keys.colors.RH_LS_IN=[0 255 0];
% keys.colors.RH_LS_CH=[0 128 0];
% 
% % population contra ipsi PSTH colors
% % keys.colors.NH_CS_IN=[255 0 64];
% % 
% % 
% %% MP changed colors for choice because he added dotted lines for that
% keys.colors.NH_CS_IN=[255 102 0];
% keys.colors.NH_CS_IN_P=[204 51 0];
% % keys.colors.NH_CS_CH=[128 0 32];
% % keys.colors.NH_CS_CH=[255 0 64];
% keys.colors.NH_CS_CH=[255 102 0];
% keys.colors.NH_CS_CH_P=[204 51 0];
% keys.colors.IH_CS_IN=[255 0 255];
% % keys.colors.IH_CS_CH=[128 0 128];
% keys.colors.IH_CS_CH=[255 0 255];
% keys.colors.CH_CS_IN=[255 128 0];
% % keys.colors.CH_CS_CH=[128 64 0];
% keys.colors.CH_CS_CH=[255 128 0];
% % keys.colors.NH_IS_IN=[0 255 255];
% keys.colors.NH_IS_IN=[0 119 255];
% keys.colors.NH_IS_IN_P=[0 70 141];
% % keys.colors.NH_IS_CH=[0 128 128];
% % keys.colors.NH_IS_CH=[0 255 255];
% keys.colors.NH_IS_CH=[0 119 255];
% keys.colors.NH_IS_CH_P=[0 70 141];
% keys.colors.IH_IS_IN=[0 128 255];
% % keys.colors.IH_IS_CH=[0 64 128];
% keys.colors.IH_IS_CH=[0 128 255];
% keys.colors.CH_IS_IN=[0 255 0];
% 
% keys.colors.CH_IS_CH=[0 128 0];
% %keys.colors.CH_IS_CH=[0 255 0];
% 
% %% these need to be fixed somehow
% keys.colors.NH_VS_IN=[150 150 150];
% keys.colors.NH_VS_CH=[80 80 80];
% keys.colors.IH_VS_IN=[64 0 255];
% keys.colors.IH_VS_CH=[32 0 128];
% keys.colors.CH_VS_IN=[255 255 0];
% keys.colors.CH_VS_CH=[128 128 0];


% population contra ipsi PSTH colors -- these are correct!
keys.colors.in_AH_CS=[255 0 64];
keys.colors.ch_AH_CS=[128 0 32];
keys.colors.in_IH_CS=[255 0 255];
keys.colors.ch_IH_CS=[128 0 128];
keys.colors.in_CH_CS=[255 128 0];
keys.colors.ch_CH_CS=[128 64 0];
keys.colors.in_AH_IS=[0 255 255];
keys.colors.ch_AH_IS=[0 128 128];
keys.colors.in_IH_IS=[0 128 255];
keys.colors.ch_IH_IS=[0 64 128];
keys.colors.in_CH_IS=[0 255 0];
keys.colors.ch_CH_IS=[0 128 0];


%% these need to be fixed somehow
keys.colors.in_AH_VS=[150 150 150];
keys.colors.ch_AH_VS=[80 80 80];
keys.colors.in_IH_VS=[64 0 255];
keys.colors.ch_IH_VS=[32 0 128];
keys.colors.in_CH_VS=[255 255 0];
keys.colors.ch_CH_VS=[128 128 0];
% preferred

keys.colors.in_AH_PF=[255 0 64];
keys.colors.ch_AH_PF=[128 0 32];
keys.colors.in_IH_PF=[255 0 255];
keys.colors.ch_IH_PF=[128 0 128];
keys.colors.in_CH_PF=[255 128 0];
keys.colors.ch_CH_PF=[128 64 0];
keys.colors.in_AH_NP=[0 255 255];
keys.colors.ch_AH_NP=[0 128 128];
keys.colors.in_IH_NP=[0 128 255];
keys.colors.ch_IH_NP=[0 64 128];
keys.colors.in_CH_NP=[0 255 0];
keys.colors.ch_CH_NP=[0 128 0];

temp_colors_I=autumn(18)*255;
temp_colors_C=winter(18)*255;
temp_colors_V=spring(18)*255;

% per position colors

keys.colors.SS_TA_ER=temp_colors_C(1,:);
keys.colors.TT_TA_ER=temp_colors_C(2,:);
keys.colors.TD_TA_ER=temp_colors_C(3,:);
keys.colors.SS_D1_ER=temp_colors_C(4,:);
keys.colors.TT_D1_ER=temp_colors_C(5,:);
keys.colors.TD_D1_ER=temp_colors_C(6,:);
keys.colors.SS_D2_ER=temp_colors_C(7,:);
keys.colors.TT_D2_ER=temp_colors_C(8,:);
keys.colors.TD_D2_ER=temp_colors_C(9,:);

keys.colors.SS_TA_SU=temp_colors_I(10,:);
keys.colors.TT_TA_SU=temp_colors_I(11,:);
keys.colors.TD_TA_SU=temp_colors_I(12,:);
keys.colors.SS_D1_SU=temp_colors_I(13,:);
keys.colors.TT_D1_SU=temp_colors_I(14,:);
keys.colors.TD_D1_SU=temp_colors_I(15,:);
keys.colors.SS_D2_SU=temp_colors_I(16,:);
keys.colors.TT_D2_SU=temp_colors_I(17,:);
keys.colors.TD_D2_SU=temp_colors_I(18,:);

% per Position & hemifield
keys.colors.SS_TA_1H_ER=temp_colors_C(1,:);
keys.colors.TT_TA_1H_ER=temp_colors_C(2,:);
keys.colors.TD_TA_1H_ER=temp_colors_C(3,:);
keys.colors.SS_D1_1H_ER=temp_colors_C(4,:);
keys.colors.TT_D1_1H_ER=temp_colors_C(5,:);
keys.colors.TD_D1_1H_ER=temp_colors_C(6,:);
keys.colors.SS_D2_1H_ER=temp_colors_C(7,:);
keys.colors.TT_D2_1H_ER=temp_colors_C(8,:);
keys.colors.TD_D2_1H_ER=temp_colors_C(9,:);

keys.colors.SS_TA_1H_SU=temp_colors_I(10,:);
keys.colors.TT_TA_1H_SU=temp_colors_I(11,:);
keys.colors.TD_TA_1H_SU=temp_colors_I(12,:);
keys.colors.SS_D1_1H_SU=temp_colors_I(13,:);
keys.colors.TT_D1_1H_SU=temp_colors_I(14,:);
keys.colors.TD_D1_1H_SU=temp_colors_I(15,:);
keys.colors.SS_D2_1H_SU=temp_colors_I(16,:);
keys.colors.TT_D2_1H_SU=temp_colors_I(17,:);
keys.colors.TD_D2_1H_SU=temp_colors_I(18,:);

keys.colors.SS_TA_2H_ER=temp_colors_C(1,:);
keys.colors.TT_TA_2H_ER=temp_colors_C(2,:);
keys.colors.TD_TA_2H_ER=temp_colors_C(3,:);
keys.colors.SS_D1_2H_ER=temp_colors_C(4,:);
keys.colors.TT_D1_2H_ER=temp_colors_C(5,:);
keys.colors.TD_D1_2H_ER=temp_colors_C(6,:);
keys.colors.SS_D2_2H_ER=temp_colors_C(7,:);
keys.colors.TT_D2_2H_ER=temp_colors_C(8,:);
keys.colors.TD_D2_2H_ER=temp_colors_C(9,:);

keys.colors.SS_TA_2H_SU=temp_colors_I(10,:);
keys.colors.TT_TA_2H_SU=temp_colors_I(11,:);
keys.colors.TD_TA_2H_SU=temp_colors_I(12,:);
keys.colors.SS_D1_2H_SU=temp_colors_I(13,:);
keys.colors.TT_D1_2H_SU=temp_colors_I(14,:);
keys.colors.TD_D1_2H_SU=temp_colors_I(15,:);
keys.colors.SS_D2_2H_SU=temp_colors_I(16,:);
keys.colors.TT_D2_2H_SU=temp_colors_I(17,:);
keys.colors.TD_D2_2H_SU=temp_colors_I(18,:);

%% Contra vs IPSI
keys.colors.SS_TA_ER_IS=temp_colors_I(1,:);
keys.colors.TT_TA_ER_IS=temp_colors_I(2,:);
keys.colors.TD_TA_ER_IS=temp_colors_I(3,:);
keys.colors.SS_D1_ER_IS=temp_colors_I(4,:);
keys.colors.TT_D1_ER_IS=temp_colors_I(5,:);
keys.colors.TD_D1_ER_IS=temp_colors_I(6,:);
keys.colors.SS_D2_ER_IS=temp_colors_I(7,:);
keys.colors.TT_D2_ER_IS=temp_colors_I(8,:);
keys.colors.TD_D2_ER_IS=temp_colors_I(9,:);

keys.colors.SS_TA_ER_VS=temp_colors_V(1,:);
keys.colors.TT_TA_ER_VS=temp_colors_V(2,:);
keys.colors.TD_TA_ER_VS=temp_colors_V(3,:);
keys.colors.SS_D1_ER_VS=temp_colors_V(4,:);
keys.colors.TT_D1_ER_VS=temp_colors_V(5,:);
keys.colors.TD_D1_ER_VS=temp_colors_V(6,:);
keys.colors.SS_D2_ER_VS=temp_colors_V(7,:);
keys.colors.TT_D2_ER_VS=temp_colors_V(8,:);
keys.colors.TD_D2_ER_VS=temp_colors_V(9,:);

keys.colors.SS_TA_ER_CS=temp_colors_C(1,:);
keys.colors.TT_TA_ER_CS=temp_colors_C(2,:);
keys.colors.TD_TA_ER_CS=temp_colors_C(3,:);
keys.colors.SS_D1_ER_CS=temp_colors_C(4,:);
keys.colors.TT_D1_ER_CS=temp_colors_C(5,:);
keys.colors.TD_D1_ER_CS=temp_colors_C(6,:);
keys.colors.SS_D2_ER_CS=temp_colors_C(7,:);
keys.colors.TT_D2_ER_CS=temp_colors_C(8,:);
keys.colors.TD_D2_ER_CS=temp_colors_C(9,:);


keys.colors.SS_TA_SU_IS=temp_colors_I(10,:);
keys.colors.TT_TA_SU_IS=temp_colors_I(11,:);
keys.colors.TD_TA_SU_IS=temp_colors_I(12,:);
keys.colors.SS_D1_SU_IS=temp_colors_I(13,:);
keys.colors.TT_D1_SU_IS=temp_colors_I(14,:);
keys.colors.TD_D1_SU_IS=temp_colors_I(15,:);
keys.colors.SS_D2_SU_IS=temp_colors_I(16,:);
keys.colors.TT_D2_SU_IS=temp_colors_I(17,:);
keys.colors.TD_D2_SU_IS=temp_colors_I(18,:);

keys.colors.SS_TA_SU_VS=temp_colors_V(10,:);
keys.colors.TT_TA_SU_VS=temp_colors_V(11,:);
keys.colors.TD_TA_SU_VS=temp_colors_V(12,:);
keys.colors.SS_D1_SU_VS=temp_colors_V(13,:);
keys.colors.TT_D1_SU_VS=temp_colors_V(14,:);
keys.colors.TD_D1_SU_VS=temp_colors_V(15,:);
keys.colors.SS_D2_SU_VS=temp_colors_V(16,:);
keys.colors.TT_D2_SU_VS=temp_colors_V(17,:);
keys.colors.TD_D2_SU_VS=temp_colors_V(18,:);

col_left      = autumn(6);

keys.colors.SS_TA_SU_CS= col_left(1,:);
keys.colors.TT_TA_SU_CS= col_left(1,:);
keys.colors.TD_TA_SU_CS= col_left(1,:);
keys.colors.SS_D1_SU_CS=col_left(6,:);
keys.colors.TT_D1_SU_CS=col_left(6,:);
keys.colors.TD_D1_SU_CS=col_left(6,:);
keys.colors.SS_D2_SU_CS= col_left(3,:);
keys.colors.TT_D2_SU_CS=col_left(3,:);
keys.colors.TD_D2_SU_CS=col_left(3,:);

% keys.colors.SS_TA_SU_CS=temp_colors_C(10,:);
% keys.colors.TT_TA_SU_CS=temp_colors_C(11,:);
% keys.colors.TD_TA_SU_CS=temp_colors_C(12,:);
% keys.colors.SS_D1_SU_CS=temp_colors_C(13,:);
% keys.colors.TT_D1_SU_CS=temp_colors_C(14,:);
% keys.colors.TD_D1_SU_CS=temp_colors_C(15,:);
% keys.colors.SS_D2_SU_CS=temp_colors_C(16,:);
% keys.colors.TT_D2_SU_CS=temp_colors_C(17,:);
% keys.colors.TD_D2_SU_CS=temp_colors_C(18,:);

%% + hemifield
keys.colors.SS_TA_1H_ER_IS=temp_colors_I(1,:);
keys.colors.TT_TA_1H_ER_IS=temp_colors_I(2,:);
keys.colors.TD_TA_1H_ER_IS=temp_colors_I(3,:);
keys.colors.SS_D1_1H_ER_IS=temp_colors_I(4,:);
keys.colors.TT_D1_1H_ER_IS=temp_colors_I(5,:);
keys.colors.TD_D1_1H_ER_IS=temp_colors_I(6,:);
keys.colors.SS_D2_1H_ER_IS=temp_colors_I(7,:);
keys.colors.TT_D2_1H_ER_IS=temp_colors_I(8,:);
keys.colors.TD_D2_1H_ER_IS=temp_colors_I(9,:);

keys.colors.SS_TA_1H_ER_VS=temp_colors_V(1,:);
keys.colors.TT_TA_1H_ER_VS=temp_colors_V(2,:);
keys.colors.TD_TA_1H_ER_VS=temp_colors_V(3,:);
keys.colors.SS_D1_1H_ER_VS=temp_colors_V(4,:);
keys.colors.TT_D1_1H_ER_VS=temp_colors_V(5,:);
keys.colors.TD_D1_1H_ER_VS=temp_colors_V(6,:);
keys.colors.SS_D2_1H_ER_VS=temp_colors_V(7,:);
keys.colors.TT_D2_1H_ER_VS=temp_colors_V(8,:);
keys.colors.TD_D2_1H_ER_VS=temp_colors_V(9,:);

keys.colors.SS_TA_1H_ER_CS=temp_colors_C(1,:);
keys.colors.TT_TA_1H_ER_CS=temp_colors_C(2,:);
keys.colors.TD_TA_1H_ER_CS=temp_colors_C(3,:);
keys.colors.SS_D1_1H_ER_CS=temp_colors_C(4,:);
keys.colors.TT_D1_1H_ER_CS=temp_colors_C(5,:);
keys.colors.TD_D1_1H_ER_CS=temp_colors_C(6,:);
keys.colors.SS_D2_1H_ER_CS=temp_colors_C(7,:);
keys.colors.TT_D2_1H_ER_CS=temp_colors_C(8,:);
keys.colors.TD_D2_1H_ER_CS=temp_colors_C(9,:);

keys.colors.SS_TA_2H_ER_IS=temp_colors_I(1,:);
keys.colors.TT_TA_2H_ER_IS=temp_colors_I(2,:);
keys.colors.TD_TA_2H_ER_IS=temp_colors_I(3,:);
keys.colors.SS_D1_2H_ER_IS=temp_colors_I(4,:);
keys.colors.TT_D1_2H_ER_IS=temp_colors_I(5,:);
keys.colors.TD_D1_2H_ER_IS=temp_colors_I(6,:);
keys.colors.SS_D2_2H_ER_IS=temp_colors_I(7,:);
keys.colors.TT_D2_2H_ER_IS=temp_colors_I(8,:);
keys.colors.TD_D2_2H_ER_IS=temp_colors_I(9,:);

keys.colors.SS_TA_2H_ER_VS=temp_colors_V(1,:);
keys.colors.TT_TA_2H_ER_VS=temp_colors_V(2,:);
keys.colors.TD_TA_2H_ER_VS=temp_colors_V(3,:);
keys.colors.SS_D1_2H_ER_VS=temp_colors_V(4,:);
keys.colors.TT_D1_2H_ER_VS=temp_colors_V(5,:);
keys.colors.TD_D1_2H_ER_VS=temp_colors_V(6,:);
keys.colors.SS_D2_2H_ER_VS=temp_colors_V(7,:);
keys.colors.TT_D2_2H_ER_VS=temp_colors_V(8,:);
keys.colors.TD_D2_2H_ER_VS=temp_colors_V(9,:);

keys.colors.SS_TA_2H_ER_CS=temp_colors_C(1,:);
keys.colors.TT_TA_2H_ER_CS=temp_colors_C(2,:);
keys.colors.TD_TA_2H_ER_CS=temp_colors_C(3,:);
keys.colors.SS_D1_2H_ER_CS=temp_colors_C(4,:);
keys.colors.TT_D1_2H_ER_CS=temp_colors_C(5,:);
keys.colors.TD_D1_2H_ER_CS=temp_colors_C(6,:);
keys.colors.SS_D2_2H_ER_CS=temp_colors_C(7,:);
keys.colors.TT_D2_2H_ER_CS=temp_colors_C(8,:);
keys.colors.TD_D2_2H_ER_CS=temp_colors_C(9,:);


keys.colors.SS_TA_1H_SU_IS=temp_colors_I(10,:);
keys.colors.TT_TA_1H_SU_IS=temp_colors_I(11,:);
keys.colors.TD_TA_1H_SU_IS=temp_colors_I(12,:);
keys.colors.SS_D1_1H_SU_IS=temp_colors_I(13,:);
keys.colors.TT_D1_1H_SU_IS=temp_colors_I(14,:);
keys.colors.TD_D1_1H_SU_IS=temp_colors_I(15,:);
keys.colors.SS_D2_1H_SU_IS=temp_colors_I(16,:);
keys.colors.TT_D2_1H_SU_IS=temp_colors_I(17,:);
keys.colors.TD_D2_1H_SU_IS=temp_colors_I(18,:);

keys.colors.SS_TA_1H_SU_VS=temp_colors_V(10,:);
keys.colors.TT_TA_1H_SU_VS=temp_colors_V(11,:);
keys.colors.TD_TA_1H_SU_VS=temp_colors_V(12,:);
keys.colors.SS_D1_1H_SU_VS=temp_colors_V(13,:);
keys.colors.TT_D1_1H_SU_VS=temp_colors_V(14,:);
keys.colors.TD_D1_1H_SU_VS=temp_colors_V(15,:);
keys.colors.SS_D2_1H_SU_VS=temp_colors_V(16,:);
keys.colors.TT_D2_1H_SU_VS=temp_colors_V(17,:);
keys.colors.TD_D2_1H_SU_VS=temp_colors_V(18,:);

keys.colors.SS_TA_1H_SU_CS=temp_colors_C(10,:);
keys.colors.TT_TA_1H_SU_CS=temp_colors_C(11,:);
keys.colors.TD_TA_1H_SU_CS=temp_colors_C(12,:);
keys.colors.SS_D1_1H_SU_CS=temp_colors_C(13,:);
keys.colors.TT_D1_1H_SU_CS=temp_colors_C(14,:);
keys.colors.TD_D1_1H_SU_CS=temp_colors_C(15,:);
keys.colors.SS_D2_1H_SU_CS=temp_colors_C(16,:);
keys.colors.TT_D2_1H_SU_CS=temp_colors_C(17,:);
keys.colors.TD_D2_1H_SU_CS=temp_colors_C(18,:);

keys.colors.SS_TA_2H_SU_IS=temp_colors_I(10,:);
keys.colors.TT_TA_2H_SU_IS=temp_colors_I(11,:);
keys.colors.TD_TA_2H_SU_IS=temp_colors_I(12,:);
keys.colors.SS_D1_2H_SU_IS=temp_colors_I(13,:);
keys.colors.TT_D1_2H_SU_IS=temp_colors_I(14,:);
keys.colors.TD_D1_2H_SU_IS=temp_colors_I(15,:);
keys.colors.SS_D2_2H_SU_IS=temp_colors_I(16,:);
keys.colors.TT_D2_2H_SU_IS=temp_colors_I(17,:);
keys.colors.TD_D2_2H_SU_IS=temp_colors_I(18,:);

keys.colors.SS_TA_2H_SU_VS=temp_colors_V(10,:);
keys.colors.TT_TA_2H_SU_VS=temp_colors_V(11,:);
keys.colors.TD_TA_2H_SU_VS=temp_colors_V(12,:);
keys.colors.SS_D1_2H_SU_VS=temp_colors_V(13,:);
keys.colors.TT_D1_2H_SU_VS=temp_colors_V(14,:);
keys.colors.TD_D1_2H_SU_VS=temp_colors_V(15,:);
keys.colors.SS_D2_2H_SU_VS=temp_colors_V(16,:);
keys.colors.TT_D2_2H_SU_VS=temp_colors_V(17,:);
keys.colors.TD_D2_2H_SU_VS=temp_colors_V(18,:);

keys.colors.SS_TA_2H_SU_CS=temp_colors_C(10,:);
keys.colors.TT_TA_2H_SU_CS=temp_colors_C(11,:);
keys.colors.TD_TA_2H_SU_CS=temp_colors_C(12,:);
keys.colors.SS_D1_2H_SU_CS=temp_colors_C(13,:);
keys.colors.TT_D1_2H_SU_CS=temp_colors_C(14,:);
keys.colors.TD_D1_2H_SU_CS=temp_colors_C(15,:);
keys.colors.SS_D2_2H_SU_CS=temp_colors_C(16,:);
keys.colors.TT_D2_2H_SU_CS=temp_colors_C(17,:);
keys.colors.TD_D2_2H_SU_CS=temp_colors_C(18,:);

%%
color_fieldnames=fieldnames(keys.colors);
for fn=1:numel(color_fieldnames)
    switch color_fieldnames{fn}(end-1:end)
        case 'CS'
            keys.colors.([color_fieldnames{fn}(1:end-2) 'PF'])=keys.colors.(color_fieldnames{fn});
        case 'IS'
            keys.colors.([color_fieldnames{fn}(1:end-2) 'NP'])=keys.colors.(color_fieldnames{fn});
    end
end

% population effector colors
keys.colors.EF_SA   =[0 255  0];
keys.colors.EF_RE   =[0 255  0];
keys.colors.EF_FG   =[0 0  255];

%% overlapping tt and cc
keys.colors.per_monkey          =[0 1 0; 1 0 0];


%% tuning table readout options (excluding particular subsets)
keys.tt.combine_tuning_properties   ={'place_name_here'}; %% additional table entry from combining columns
keys.tt.perturbations               = 0;
keys.tt.choices                     = 0;
keys.tt.hands                       = 0;
keys.tt.IC_for_criterion            = 'in';
keys.tt.trial_criterion_in          ='per_position';
keys.tt.trial_criterion_ch          ='per_hemifield';
keys.tt.selection                   ={};
keys.tt.unselect                    ={};
keys.tt.selected_list               ={};
keys.tt.unselected_list             ={};
keys.tt.type_effectors              ={'Msac'};

%% population
keys.sct=struct([]);
keys.ccs=struct([]);
keys.ons=struct([]);
keys.pop=struct([]);
keys.hst=struct([]);

%% EPOCH SETTINGS 
% For each task type seperately, analysis epochs are defined

%% Fixation only type
keys.EPOCHS_PER_TYPE{1}={...
    'INI',      2,  -0.4,	-0.1,   'INI';...
    'Facq',     3,	0.05,	0.15,   'INI';...
    'Fhol',     4,	-0.3,	0,      'INI';...
    'PreS',     60,	-0.1,	-0.01,  'INI';...
    'PeriS',    60,	-0.01,	0.05,   'INI';...
    'PreR',     62,	-0.3, 	-0.01,  'INI';...
    'PeriR',	62,	-0.05, 	0.15,   'INI';...
    'Thol',     20,	-0.3,	0,      'INI';...
    };
keys.WINDOWS_PER_TYPE{1}={...
    'Aquisition',   3,	-0.5,	0.4;...
    'Fixation',     4,	-1.2,   0.17;...
    'Saccade',      60,	-0.1,   0.6;...
    'Reach',        62,	-0.35,  0.7;...
    };
keys.ANOVAS_PER_TYPE(1).epoch={'INI' 'Facq';...
    'Facq' 'Fhol';...
    'Fhol' 'PreS';...
    'Fhol' 'PeriS';...
    'Fhol' 'PreR';...
    'Fhol' 'PeriR';...
    'Fhol' 'Thol'};

keys.ANOVAS_PER_TYPE(1).hemifield            ={'Facq','Fhol','PreS','PeriS','PreR','PeriR','Thol'}';
keys.ANOVAS_PER_TYPE(1).positions          ={'Facq','Fhol','PreS','PeriS','PreR','PeriR','Thol'}';
keys.ANOVAS_PER_TYPE(1).hands              ={'Facq','Fhol','PreS','PeriS','PreR','PeriR','Thol'}';
keys.ANOVAS_PER_TYPE(1).SxH                ={'Facq','Fhol','PreS','PeriS','PreR','PeriR','Thol'}';
keys.ANOVAS_PER_TYPE(1).main               ={'Facq','Fhol','PreS','PeriS','PreR','PeriR','Thol'}';

%% Visually guided type
keys.EPOCHS_PER_TYPE{2}={...
    'INI',      2,	-0.4,	-0.1,   'INI';...
    'Facq',     3,	0.05,	0.15,   'INI';...
    'Fhol',     4,	-0.3,	0,      'INI';...
    'Cue',      4,	0.05,	0.15,   'INI';...
    'PreS',     60,	-0.1,	-0.01,  'INI';...
    'PeriS',    60,	-0.01,	0.05,   'INI';...
    'Tacq',     5,	0,      0.15,   'INI';...
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
keys.ANOVAS_PER_TYPE(2).hemifield            ={'Cue','PreS','PeriS','Tacq','Thol'}';
keys.ANOVAS_PER_TYPE(2).positions          ={'Cue','PreS','PeriS','Tacq','Thol'}';
keys.ANOVAS_PER_TYPE(2).hands              ={'Facq','Fhol','Cue','PreS','PeriS','Tacq','Thol'}';
keys.ANOVAS_PER_TYPE(2).SxH                ={'Cue','PreS','PeriS','Tacq','Thol'}';
keys.ANOVAS_PER_TYPE(2).main               ={'Facq','Fhol','Cue','PreS','PeriS','Tacq','Thol'}';

%% Memory type
keys.EPOCHS_PER_TYPE{3}={...
    'INI',      2,	-0.4,	-0.1,   'INI';...
    'Facq',     3,	0.05,	0.15,   'INI';...
    'Fhol',     6,	-0.3,	0,      'INI';...
    'Cue',      6,	0.05,   0.15,   'INI';...
    'MemE',     7, 	0,      0.2,    'INI';...
    'MemL',     9,	-0.3, 	0,      'INI';...
    'PreS',     60,	-0.1, 	-0.01,  'INI';...
    'PeriS',	60,	-0.01, 	0.05,   'INI';...
    'TIhol',	10,	0,      0.1,    'INI';...
    'Thol',     20,	-0.3,   0,      'INI';...
    };
keys.WINDOWS_PER_TYPE{3}={...
    'Initiation',   2,	-0.5,	0;...
    'Cue',          6,	-0.8,   0.78;...
    'Saccade',      60,	-0.7,   0.05;...
    'T hold',       20,	-0.6,   0;...
    };
keys.ANOVAS_PER_TYPE(3).epoch={'INI' 'Facq';...
    'INI' 'Fhol';...
    'Fhol' 'Cue';...
    'Fhol' 'MemE';...
    'Fhol' 'MemL';...
    'Fhol' 'PreS';...
    'Fhol' 'PeriS';...
    'Fhol' 'TIhol';...
    'Fhol' 'Thol';...
    };
keys.ANOVAS_PER_TYPE(3).hemifield            ={'Cue','MemE','MemL','PreS','PeriS','TIhol','Thol'}';
keys.ANOVAS_PER_TYPE(3).positions             ={'Cue','MemE','MemL','PreS','PeriS','TIhol','Thol'}';
keys.ANOVAS_PER_TYPE(3).hands              ={'Facq','Fhol','Cue','MemE','MemL','PreS','PeriS','TIhol','Thol'}';
keys.ANOVAS_PER_TYPE(3).SxH                ={'Cue','MemE','MemL','PreS','PeriS','TIhol','Thol'}';
keys.ANOVAS_PER_TYPE(3).main              ={'Facq','Fhol','Cue','MemE','MemL','PreS','PeriS','TIhol','Thol'}';

%% Delay type
keys.EPOCHS_PER_TYPE{4}={...
    'INI',      2,	-0.4,	-0.1,   'INI';...
    'Facq',     3,	-0.4,	-0.1,   'INI';...
    'Fhol',     6,	-0.3,	0,      'INI';...
    'Cue',      6,	0.05,   0.15,   'INI';...
    'EDel',     4, 	-0.6,   -0.3,   'INI';...
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
    'Initiation',   2,	-0.5,	0.4;...
    'Fixation',     3,	-1.2,   0.17;...
    'Delay Period', 6,	-0.33,  1.35;...
    'Saccade',      60,	-0.1,   0.6;...
    'Reach',        62,	-0.35,  0.7;...
    };
keys.ANOVAS_PER_TYPE(4).epoch={'INI' 'Facq';...
    'INI' 'Fhol';...
    'Fhol' 'Cue';...
    'Fhol' 'EDel';...
    'EDel' 'Del';...
    'EDel' 'PreS';...
    'EDel' 'PeriS';...
    'EDel' 'PostS';...
    'EDel' 'PreR';...
    'EDel' 'PeriR';...
    'EDel' 'PostR';...
    'Fhol' 'Thol';...
    };
keys.ANOVAS_PER_TYPE(4).hemifield            ={'INI','Facq','Fhol','Cue','EDel','Del','PreS','PeriS','PostS','PreR','PeriR','PostR','Thol'}';
keys.ANOVAS_PER_TYPE(4).positions          ={'INI','Facq','Fhol','Cue','EDel','Del','PreS','PeriS','PostS','PreR','PeriR','PostR','Thol'}';
keys.ANOVAS_PER_TYPE(4).hands              ={'INI','Facq','Fhol','Cue','EDel','Del','PreS','PeriS','PostS','PreR','PeriR','PostR','Thol'}';
keys.ANOVAS_PER_TYPE(4).SxH                ={'INI','Facq','Fhol','Cue','EDel','Del','PreS','PeriS','PostS','PreR','PeriR','PostR','Thol'}';
keys.ANOVAS_PER_TYPE(4).main                ={'INI','Facq','Fhol','Cue','EDel','Del','PreS','PeriS','PostS','PreR','PeriR','PostR','Thol'}';

%% M2S
keys.EPOCHS_PER_TYPE{5}={...
    'INI',      2,      -0.4,	-0.1,   'INI';...
    'Facq',     3,      0.05,	0.15,   'INI';...
    'Fhol',     6,      -0.3,	0,      'INI';...
    'Cue',      6,      0.05,   0.2,   'INI';...
    'MemE',     7,      0,      0.2,    'INI';...
    'PreS',     60,	    -0.1, 	-0.01,  'INI';...
    'PeriS',	60,     -0.01, 	0.05,   'INI';...
    'Exp',      14,     -0.8, 	-0.2,   'INI';...
    'Thol',     90,     -0.95,  -0.8,   'INI';...
    };
keys.WINDOWS_PER_TYPE{5}={...
    'INI',          2,	-0.5,	0.3     ;...
    'Sample',       6,	-0.4,   2       ;...
    'Selection',    90,	-2,     0       ;...
    };
keys.ANOVAS_PER_TYPE(5).epoch={'INI' 'Facq';...
    'Facq' 'Fhol';...
    'Fhol' 'Cue';...
    'Fhol' 'MemE';...
    'Fhol' 'PreS';...
    'Fhol' 'PeriS';...
    'Fhol' 'Thol';...
    };
keys.ANOVAS_PER_TYPE(5).hemifield            ={'Cue','MemE','PreS','PeriS','Thol'}';
keys.ANOVAS_PER_TYPE(5).positions          ={'Cue','MemE','PreS','PeriS','Thol'}';
keys.ANOVAS_PER_TYPE(5).hands              ={'Facq','Fhol','Cue','MemE','PreS','PeriS','Thol'}';
keys.ANOVAS_PER_TYPE(5).SxH                ={'Cue','MemE','PreS','PeriS','Thol'}';
keys.ANOVAS_PER_TYPE(5).main              ={'Facq','Fhol','Cue','MemE','PreS','PeriS','Thol'}';

%% M2S
keys.EPOCHS_PER_TYPE{6}={...
    'INI',      2,      -0.4,	-0.1,   'INI';...
    'Facq',     3,      0.05,	0.15,   'INI';...
    'Fhol',     6,      -0.3,	0,      'INI';...
    'Cue',      6,      0.05,   0.2,   'INI';...
    'MemE',     7,      0,      0.2,    'INI';...
    'PreS',     60,	    -0.1, 	-0.01,  'INI';...
    'PeriS',	60,     -0.01, 	0.05,   'INI';...
    'Exp',      14,     -0.8, 	-0.2,   'INI';...
    'Thol',     90,     -0.95,  -0.8,   'INI';...
    };
keys.WINDOWS_PER_TYPE{6}={...
    'INI',          2,	-0.5,	0.3     ;...
    'Sample',       6,	-0.4,   2       ;...
    'Selection',    90,	-2,     0       ;...
    };
keys.ANOVAS_PER_TYPE(6).epoch={'INI' 'Facq';...
    'Facq' 'Fhol';...
    'Fhol' 'Cue';...
    'Fhol' 'MemE';...
    'Fhol' 'PeriS';...
    'Fhol' 'Exp';...
    'Fhol' 'Thol';...
    };
keys.ANOVAS_PER_TYPE(6).hemifield            ={'Cue','MemE','PeriS','Exp','Thol'}';
keys.ANOVAS_PER_TYPE(6).positions          ={'Cue','MemE','PeriS','Exp','Thol'}';
keys.ANOVAS_PER_TYPE(6).hands              ={'Facq','Fhol','Cue','MemE','PeriS','Exp','Thol'}';
keys.ANOVAS_PER_TYPE(6).SxH                ={'Cue','MemE','PreS','Exp','Thol'}';
keys.ANOVAS_PER_TYPE(6).main              ={'Facq','Fhol','Cue','MemE','PeriS','Exp','Thol'}';

