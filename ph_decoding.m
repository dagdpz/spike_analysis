function ph_decoding(population,keys)
%%


wn=2;
type=4;
effector=4;
choice=0;

keys.PSTH_binwidth=0.001;
keys.gaussian_kernel=0.001;
keys.kernel_type='box';


%% define conditions to look at
all_trialz=[population.trial];
per_trial.types       =[all_trialz.type];
per_trial.effectors   =[all_trialz.effector];

tr_con=ismember([all_trialz.completed],keys.cal.completed);
[whatisthis]=ph_arrange_positions_and_plots(keys,all_trialz(tr_con));

per_trial.types       =[all_trialz.type];
per_trial.effectors   =[all_trialz.effector];
per_trial.hands       =[all_trialz.reach_hand];
per_trial.choice      =[all_trialz.choice];
per_trial.perturbation=[all_trialz.perturbation];
per_trial.hemifield   =[whatisthis.trial.hemifield];
per_trial.perturbation(ismember(per_trial.perturbation, keys.cal.perturbation_groups{1}))=0;
per_trial.perturbation(ismember(per_trial.perturbation, keys.cal.perturbation_groups{2}))=1;

uq.hemifield=unique(per_trial.hemifield); %[-1,0,1]; % why does this have to be hardcoded? ---> Because case not defined yet, case defines positions !!

uq.reach_hand     =unique(per_trial.hands);
uq.choice    =unique(per_trial.choice);
uq.effector    =unique(per_trial.effectors);
uq.perturbation    =unique(per_trial.perturbation);
uq.perturbation=uq.perturbation(~isnan(uq.perturbation));
uq.handspace    =combvec(uq.reach_hand,uq.hemifield);


condition_parameters           = {'effector','reach_hand','hemifield','choice','perturbation'};
conditions_matrix              = combvec(uq.effector,uq.reach_hand,uq.hemifield,uq.choice, uq.perturbation)';


keys.PSTH_WINDOWS=keys.WINDOWS_PER_TYPE{type};

for wn= 3:size(keys.PSTH_WINDOWS,1)
    window_name=keys.PSTH_WINDOWS{wn,1};
    
    path_to_save=[keys.basepath_to_save keys.project_version filesep 'decoding'];
    raster_data_directory_name=[path_to_save filesep window_name '_rasters' filesep];
    if ~exist(raster_data_directory_name,'dir')
        mkdir(path_to_save, [window_name '_rasters']);
    end
    
    
    for u=1:numel(population)
        clear raster_data raster_labels raster_site_info
        unit=population(u);
        tr = [unit.trial.completed]==1 ;%&  [unit.trial.effector]==effector &  [unit.trial.choice]==choice ;
        unit.trial=unit.trial(tr);
        
        [unit]=ph_arrange_positions_and_plots(keys,unit.trial,unit);
        unit=ph_LR_to_CI(keys,unit);
        
        for t=1:numel(unit.trial)
            raster_data(t,:)=ph_spike_density(unit.trial(t),wn,keys,0,1);
            raster_labels.handspace{t}=num2str([unit.trial(t).reach_hand unit.trial(t).hemifield]);
            raster_labels.reach_hand{t}=num2str(unit.trial(t).reach_hand);
            raster_labels.hemifield{t}=num2str(unit.trial(t).hemifield);
            raster_labels.effector{t}=num2str(unit.trial(t).effector);
            raster_labels.choice{t}=num2str(unit.trial(t).choice);
            raster_labels.perturbation{t}=num2str(unit.trial(t).perturbation);
        end
        
        raster_site_info.session_ID=unit.session;
        raster_site_info.recording_channel= unit.channel;
        raster_site_info.unit=[unit.block_unit{:}];
        raster_site_info.alignment_event_time=keys.PSTH_WINDOWS{wn,3}*-1/0.001;
        
        save([raster_data_directory_name unit.unit_ID '.mat'],'raster_data','raster_labels','raster_site_info');
    end
    
    %%  4.  Bin the data
    save_prefix_name = [path_to_save filesep window_name '_binned'];
    bin_width = 30;
    step_size = 10;
    binned_data_file_name = create_binned_data_from_raster_data(raster_data_directory_name, save_prefix_name, bin_width, step_size);
    
    
    
    
    %%  5.  Loop through parameters and subsets for each parameter
    %%% here add removing certain trials dependent on what parameters we are
    %%% decoding
    decoding_parameters={'handspace','reach_hand','hemifield','effector'};
    
    load(binned_data_file_name);  % load the binned data
    all_binned_labels=binned_labels;
    all_binned_data=binned_data;
    FN_labels=fieldnames(all_binned_labels);
    
    for par=decoding_parameters
        
        parameter=par{:};
        par_idx=find(~ismember(condition_parameters,parameter));
        switch parameter
            case 'handspace'
                par_idx=find(~ismember(condition_parameters,{'reach_hand','hemifield'}));
        end
        CM=unique(conditions_matrix(:,par_idx),'rows');
        CP=condition_parameters(par_idx);
        
        for subset=1:size(CM,1)
            clear binned_data binned_labels
            
            subset_name=[parameter '_'];
            for c=1:numel(CP)
                subset_name=[subset_name CP{c}(1:2) num2str(CM(subset,c))];
            end
            for u= 1:numel(all_binned_data)
                subset_index=true(size(all_binned_data{u},1),1);
                subset_index(any(isnan(all_binned_data{u}),2))=false;
                for c=1:numel(CP)
                    subset_index(~ismember(all_binned_labels.(CP{c}){u},num2str(CM(subset,c))))=false;
                end
                binned_data{u}=all_binned_data{u}(subset_index,:);
                for fn=1:numel(FN_labels)
                    binned_labels.(FN_labels{fn}){u}=all_binned_labels.(FN_labels{fn}){u}(subset_index);
                end
                %% make sure only units including a) data for the given subset and b) all possible values of the decoded parameter are taken!
                unit_valid(u)=~isempty(binned_data{u}) && numel(unique(binned_labels.(parameter){u}))==size(uq.(parameter),2);
            end
            % remove empty units
            binned_data=binned_data(unit_valid);
            for fn=1:numel(FN_labels)
                binned_labels.(FN_labels{fn})=binned_labels.(FN_labels{fn})(unit_valid);
            end
            if isempty(binned_data)
                continue;
            end
            
            
            %% save file !!
            subset_data_filename = [path_to_save filesep window_name '_' num2str(bin_width) 'ms_bins_' num2str(step_size) 'ms_sampled' '_' subset_name];
            save(subset_data_filename, 'binned_data', 'binned_labels', 'binned_site_info');
            
            %% Calculate how many times each stimulus has been shown to each neuron ==> THIS IS NOT USED AFTERWARDS...??
            for i = 0:60
                inds_of_sites_with_at_least_k_repeats = find_sites_with_k_label_repetitions(binned_labels.(parameter), i);
                num_sites_with_k_repeats(i + 1) = length(inds_of_sites_with_at_least_k_repeats);
            end
            
            
            %%  Begin the decoding analysis
            %%  6.  Create a datasource object
            % we will use object identity labels to decode which object was shown (disregarding the position of the object)
            
            % use 20 cross-validation splits (which means that 19 examples of each object are used for training and 1 example of each object is used for testing)
            num_cv_splits = 20;
            
            % create the basic datasource object
            ds = basic_DS(subset_data_filename, parameter,  num_cv_splits);
            
            ds.sites_to_use = find_sites_with_k_label_repetitions(binned_labels.(parameter),num_cv_splits);
            if isempty(ds.sites_to_use)
                continue;
            end
            
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
            results_file_name = [path_to_save filesep window_name  '_decoding_results_' subset_name];
            save(results_file_name, 'DECODING_RESULTS');
            
            % if logged the code that was run using a log_code_object, save the code
            %LOGGED_CODE = log_code_obj.return_logged_code_structure;
            %save(save_file_name, '-v7.3', 'DECODING_RESULTS', 'LOGGED_CODE');
            
            
            
            %%  12.  Plot the basic results
            
            
            % which results should be plotted (only have one result to plot here)
            result_names{1} = results_file_name;
            
            % create an object to plot the results
            clear plot_obj
            plot_obj = plot_standard_results_object(result_names);
            
            % put a line at the time when the stimulus was shown
            plot_obj.significant_event_times = 0;
            
            
            % optional argument, can plot different types of results
            %plot_obj.result_type_to_plot = 2;  % for example, setting this to 2 plots the normalized rank results
            
            
            plot_obj.plot_results;   % actually plot the results
            
            filename=[window_name '_' subset_name '_decoding_accuracy'];
            plot_title=[window_name '_' subset_name];
            ph_title_and_save(gcf,filename,plot_title,keys)
            
            %% this part i dont know what it actually is... ? ==> guessing its taking decoding weights based on each bin and plots the entire timecourse for those weights
            % %%  13.  Plot the TCT matrix
            %
            % plot_obj = plot_standard_results_TCT_object(save_file_name);
            %
            % plot_obj.significant_event_times = 0;   % the time when the stimulus was shown
            %
            %
            % % optional parameters when displaying the TCT movie
            % %plot_obj.movie_time_period_titles.title_start_times = [-500 0];
            % %plot_obj.movie_time_period_titles.title_names = {'Fixation Period', 'Stimulus Period'}
            %
            % plot_obj.plot_results;  % plot the TCT matrix and a movie showing if information is coded by a dynamic population code
        end
    end
end
end
