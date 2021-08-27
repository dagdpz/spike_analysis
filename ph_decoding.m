function ph_decoding(population,keys)

wn=2;
type=4;
effector=4;
choice=0;

keys.PSTH_binwidth=0.001;
keys.gaussian_kernel=0.001;
keys.kernel_type='box';

keys.PSTH_WINDOWS=keys.WINDOWS_PER_TYPE{type};


path_to_save=[keys.basepath_to_save keys.project_version filesep 'decoding'];
raster_data_directory_name=[path_to_save filesep 'rasters'];
if ~exist(raster_data_directory_name,'dir')
    mkdir(path_to_save, 'rasters');
end


for u=1:numel(population)
    clear raster_data raster_labels raster_site_info
    unit=population(u);
    tr = [unit.trial.completed]==1 &  [unit.trial.effector]==effector &  [unit.trial.choice]==choice ;
    unit.trial=unit.trial(tr);
    
    [unit]=ph_arrange_positions_and_plots(keys,unit.trial,unit);
    unit=ph_LR_to_CI(keys,unit);
    
    for t=1:numel(unit.trial)
        raster_data(t,:)=ph_spike_density(unit.trial(t),wn,keys,0,1);
        raster_labels.handspace{t}=num2str([unit.trial(t).reach_hand unit.trial(t).hemifield]);
    end
    
    raster_site_info.session_ID=unit.session;
    raster_site_info.recording_channel= unit.channel;
    raster_site_info.unit=[unit.block_unit{:}];
    raster_site_info.alignment_event_time=keys.PSTH_WINDOWS{wn,3}*-1/0.001;
    
    save([raster_data_directory_name filesep unit.unit_ID '.mat'],'raster_data','raster_labels','raster_site_info');
end

%%  4.  Bin the data
save_prefix_name = 'Binned_data';
bin_width = 30; 
step_size = 10;  
binned_data_file_name = create_binned_data_from_raster_data(raster_data_directory_name, save_prefix_name, bin_width, step_size);
 

%%  5.  Calculate how many times each stimulus has been shown to each neuron
load(binned_data_file_name);  % load the binned data
specific_binned_labels_names = 'handspace';
for i = 0:60
    inds_of_sites_with_at_least_k_repeats = find_sites_with_k_label_repetitions(binned_labels.handspace, i);
    num_sites_with_k_repeats(i + 1) = length(inds_of_sites_with_at_least_k_repeats);
end


%%  Begin the decoding analysis  %%
%%  6.  Create a datasource object
% we will use object identity labels to decode which object was shown (disregarding the position of the object)

% use 20 cross-validation splits (which means that 19 examples of each object are used for training and 1 example of each object is used for testing)
num_cv_splits = 20; 

% create the basic datasource object
ds = basic_DS(binned_data_file_name, specific_binned_labels_names,  num_cv_splits);



% other useful options:

% if using the Poison Naive Bayes classifier, load the data as spike counts by setting the load_data_as_spike_counts flag to 1
%ds = basic_DS(binned_data_file_name, specific_binned_labels_names,  num_cv_splits, 1);

% can have multiple repetitions of each label in each cross-validation split (which is a faster way to run the code that uses most of the data)
%ds.num_times_to_repeat_each_label_per_cv_split = 2;

 % optionally can specify particular sites to use
%ds.sites_to_use = find_sites_with_k_label_repetitions(the_labels_to_use, num_cv_splits);  

% can do the decoding on a subset of labels
%ds.label_names_to_use =  {'kiwi', 'flower', 'guitar', 'hand'};




%%   7.  Create a feature preprocessor object

% create a feature preprocess that z-score normalizes each feature
the_feature_preprocessors{1} = zscore_normalize_FP;  


% other useful options:   

% can include a feature-selection features preprocessor to only use the top k most selective neurons
% fp = select_or_exclude_top_k_features_FP;
% fp.num_features_to_use = 25;   % use only the 25 most selective neurons as determined by a univariate one-way ANOVA
% the_feature_preprocessors{2} = fp;




%%  8.  Create a classifier object 

% select a classifier
the_classifier = max_correlation_coefficient_CL;


% other useful options:   

% use a poisson naive bayes classifier (note: the data needs to be loaded as spike counts to use this classifier)
%the_classifier = poisson_naive_bayes_CL;  

% use a support vector machine (see the documentation for all the optional parameters for this classifier)
%the_classifier = libsvm_CL;


%%  9.  create the cross-validator 


the_cross_validator = standard_resample_CV(ds, the_classifier, the_feature_preprocessors);  

the_cross_validator.num_resample_runs = 2;  % usually more than 2 resample runs are used to get more accurate results, but to save time we are using a small number here


% other useful options:   

% can greatly speed up the run-time of the analysis by not creating a full TCT matrix (i.e., only trainging and testing the classifier on the same time bin)
% the_cross_validator.test_only_at_training_times = 1;  




%%  10.  Run the decoding analysis   

% if calling the code from a script, one can log the code so that one can recreate the results 
%log_code_obj = log_code_object;
%log_code_obj.log_current_file; 


% run the decoding analysis 
DECODING_RESULTS = the_cross_validator.run_cv_decoding; 



%%  11.  Save the results

% save the results
save_file_name = [path_to_save 'decoding_results'];
save(save_file_name, 'DECODING_RESULTS'); 

% if logged the code that was run using a log_code_object, save the code
%LOGGED_CODE = log_code_obj.return_logged_code_structure;
%save(save_file_name, '-v7.3', 'DECODING_RESULTS', 'LOGGED_CODE'); 



%%  12.  Plot the basic results


% which results should be plotted (only have one result to plot here)
result_names{1} = save_file_name;

% create an object to plot the results
plot_obj = plot_standard_results_object(result_names);

% put a line at the time when the stimulus was shown
plot_obj.significant_event_times = 0;   


% optional argument, can plot different types of results
%plot_obj.result_type_to_plot = 2;  % for example, setting this to 2 plots the normalized rank results


plot_obj.plot_results;   % actually plot the results





%%  13.  Plot the TCT matrix

plot_obj = plot_standard_results_TCT_object(save_file_name);

plot_obj.significant_event_times = 0;   % the time when the stimulus was shown


% optional parameters when displaying the TCT movie
%plot_obj.movie_time_period_titles.title_start_times = [-500 0];
%plot_obj.movie_time_period_titles.title_names = {'Fixation Period', 'Stimulus Period'}

plot_obj.plot_results;  % plot the TCT matrix and a movie showing if information is coded by a dynamic population code
end
