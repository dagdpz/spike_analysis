%% Linus ephys [20160203 20160506] datasets

clc
clear batch,
close all
warning('off','all')
% EH_compare_groups(batch,testing)

global GLO
GLO.accuracy_as_absolute            =   1; %% 1 meaning it computes the averages of x and y first and then creates the euclidean
GLO.delete_last                     =   0;
GLO.fontsize_titles_big             =   16;
GLO.fontsize_small                  =   8;
GLO.fontsize_ticks                  =   12;
GLO.fontsize_labels                 =   12;
GLO.linewidth                       =   2;
GLO.plot_raw_endpoints              =   1; %1 means 1 point per trial, 0 means average across trials
GLO.calculate_statististics         =   1;
GLO.parametric_testing              =   0;
GLO.plot_statististics              =   1;
GLO.plot_it                         =   1;
GLO.create_pdf                      =   1;
GLO.append_pdfs                     =   0;
GLO.parent_folder                   =   '';
GLO.folder_to_save                  =   'Y:\Projects\Pulv_eye_gaze_position\behavior\Linus';
GLO.type_of_free_gaze               =   '6';
GLO.one_subject                     =   0;
GLO.trial_by_trial                  =   0; % for statistics
GLO.CDF                             =   0;
GLO.text_in_plot                    =   0;
GLO.same_day                        =   0;
GLO.testing_patient                 =   0;
GLO.instructed_only                 =   0;
GLO.choice_only                     =   0;
GLO.only_significant                =   1; % for sigstar
GLO.only_success_for_accuracy       =   0;
GLO.only_between_group_comparissons =   0;
GLO.point_per_batch                 =   1; %0 average across session , 1 % 1 point per run
GLO.summary                         =   [-1]; %which plots
GLO.target_locations_in_raw         =   0;
GLO.saccade_in_raw                  =   0;
GLO.modify_positions                =   0;
GLO.euclideans_reach                =   [-15, 15];
GLO.trial_numbers                   =   0;
GLO.keep_raw_output                 =   0;
GLO.hits_in_plot                    =   0;
GLO.min_hits                        =   0; %or 1 for 50 hits min
GLO.eccentricity                    =   'ALL'; %%new

% steady.passfilter                   =   {'saccades','lat',0.08, 0.5;'reaches','lat',0.1, 10};

%steady.passfilter                   =   {'saccades','lat',0, 3;'reaches','lat',0.2, 10;}; % add
steady.passfilter                   =   {'saccades','lat',0, 3;'reaches','lat',0, 10;}; % add

% steady.passfilter                   =   {'saccades','lat',0.01, 10;'reaches','lat',0.01, 10;}; %% !!!!!! remove  


GLO.remove_outliers                 =   0;
% 1 ERROR BARS
% 2 HISTOGRAMS
% 3 ACCURACY
% 4 CORRELATIONS
% 5 HAND RELEASE
% 6 RAW XY 2D 
% 7X X VS TIME
% 7Y Y VS TIME
% 8 FILELIST
% 9 ERRORS
% 10 CH - IN

% % 1: saccade, that ended up closest to the target and was big enough
% % 2: biggest saccade that ended up close enough
% % 3: last saccade in the state
% % 4: first saccade in the state
% % 5: the first that is bigger than 'sac_min_amp'
steady.reach_1st_pos                =   1; %%new used for all but
steady.reach_1st_pos_in             =   0; %%new
steady.reach_pos_at_state_change    =   0; %%new

steady.remove_outliers              =   GLO.remove_outliers;
steady.sac_ini_t                  = 200;
steady.sac_end_t                  = 50;
%steady.min_dur                    = 0.03;
steady.eyetracker_sample_rate     = 220; % Hertz
steady.correct_offset             = 1;
steady.saccade_definition         = 4;
steady.smoothing_samples          = 15; %downsampling does a diff to find relevant changing samples, then it inerpolates to 1 ms between such samples and then it smoothens, the smoothing is set such that it grabs one sample before and one after plus a bit more but not reaching 2 samples 4.54 ms each sample of the eyetracker
steady.sac_min_amp                = 2;
steady.keep_raw_data              = 1;
steady.correlation_mode           = 'pearson';
steady.display                    = 0;
steady.downsampling               = 1;
steady.saccade_definition         = 4;
        steady.tar_range_x                  =   [NaN;NaN];
        steady.tar_range_y                  =   [NaN;NaN];

load('Y:\Projects\PPC_pulv_eye_hand\ephys\MIP_dPul_inj_working_post_sfn\behaviour_filelist.mat');

filelist_formatted_control=filelist_formatted.Lin_MIP_R_PT0_Ddre_han;
filelist_formatted_inactivation=filelist_formatted.Lin_MIP_R_PT1_Ddre_han;

% Linus dPul inactivation MIP recordings datasets
subject_ID{1}='Control';
group{1}                        = repmat({'Linus_phys'},size(filelist_formatted_control,1),1);
dates_subject_in{1}             = filelist_formatted_control(:,1);
batching{1}.runs                = filelist_formatted_control(:,2);  % either empty or specific runs specified
batching_type{1}                = 1; % 1 run by run, 2 session by session, 3 group by group
batching{1}.range_of_dates      = 0;

subject_ID{2}='Experimental';
group{2}                        = repmat({'Linus_phys'},size(filelist_formatted_inactivation,1),1);
dates_subject_in{2}                = filelist_formatted_inactivation(:,1);
batching{2}.runs                = filelist_formatted_inactivation(:,2);  % either empty or specific runs specified
batching_type{2}                = 1; % 1 run by run, 2 session by session, 3 group by group
batching{2}.range_of_dates      = 0;


%% over all
GLO.folder_to_save                  = 'Y:\Projects\PPC_pulv_eye_hand\behavior\beh_analysis\Linus_dPul_inj_MIP_per_run';
run beh_run_analysis




