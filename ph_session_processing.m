function  keys = ph_session_processing(modified_keys)
% ph_session_processing is the first step in the spike analysis pipeline.
% It is performed on all data from one monkey, the analyzed data is defined by a range of dates (keys.monkey.dates) or by specific runs (typically found in ph_additional_settings and has priority),
% datasets (to be implemented), and various selections of trial type, effector, reach_hand, choice (defined in keys.cal)
% All analysis is performed per session, so the main loop is through sessions, because the same unit may appear several times in the same session.
% The first step is reading in the data per block and applying formatting and crude behavioral analysis (for reaction times mainly) using monkeypsych_analyze.
% Then, in ph_run_state_alignment_per_trial several things happen:
% 1) events are defined (i.e. state changes and saccade/reach onset/offset). state changes and behavioral events are treated in the same way for all future analysis
% 2) Spikes from the previous trial are appended (to allow analysis of firing rates before state 2 (which is time 0) later on, e.g. INI_TRI baseline)
% 3) nonmatching trial types are excluded (defined in keys.cal)
% 4) entries from the sorted_neurons excel table are added per UNIT AND RUN/BLOCK (including unit_ID, recording location and stability/single/SNR rankings)
% 5) If a unit is not found in the excel table (by session, run/block, channel, and unit (a,b,c,...), a generic unit_ID will be given containing the information about session/run/block/channel/unit.
%     This is relevant to allow analysis without filling in the excel table and for resorting purposes based on tuning (i.e. to decide about what to combine in the same_cells file)
% After doing so for all runs in the session, the data is restructured in order to combine data from the same unit across different runs/blocks
% Note that because stability/single/SNR rankings can differ for the same unit in different blocks, the selection criteria regarding rankings are applied here already, in the following way:
% the sorted_neurons table entries for units that dont meet the criteria in A GIVEN RUN/BLOCK is not taken over in ph_run_state_alignment_per_trial,
% so that these runs/blocks will not be combined with other runs/blocks with the same unit (because it now has a generic unit_ID)
% Spike waveforms are plotted independently for units that meet the criteria (i.e. are found in the sorted_neurons table) and units which don't.
% Units which were not found in the sorted_neurons table are excluded (and therefore stability/single/SNR criteria are applied) only if keys.cal.units_from_sorting_table is set to 1,
% so set keys.cal.units_from_sorting_table to 1 once the sorting is completed.
% Finally, single_cell plots (ph_plot_unit_per_condition) are created (if keys.plot.single_cells is selected)
% and a reduced version (removing hand/eye traces and individual spike shapes) of the population data is stored PER SESSION.
% Because we want to also have an overview of the significant responses in the per unit plots, ANOVAS and post-hoc comparisons (ph_ANOVAS) are applied before single cell plotting.
% For the epoch wise ANOVAs, we first need to calculate firing rates in the defined epochs (done in ph_epochs),
% and mark trials which exceed the confidence intervals of firing rates (ph_accept_trials_per_unit) to not take them into consideration.
% The latter marks trials as not accepted when firing rates accross the entire trial exceed confidence intervals (across trials) on a logarithmic scale.
% This means in particular trials with no spike at all will not be taken into account further (which allows excluding trials in plexon sort code by cutting off in time!).
% ph_ANOVAS will create a table containing all analysis results for each unit. This table is stored in the keys, and updated across sessions,
% such that in the end we have one table containing all ANOVA results.
% Because that is a very convenient format (one row per unit), other information like minimum number of trials per conditition and firing rates/modulation per epoch is also stored in this table.
% It should be noted here, that due to the importance of these table entries for population analysis, they can also be recalculated independently (using ph_initiate_analyis).
% IMPORTANT: because positions might be defined differently dependent on condition arrangement (assigned in ph_arrange_positions_and_subplots), ANOVAS are calculated indepently for each condition arrangement .
% As the name of ph_arrange_positions_and_subplots already says, this function not only defines positions, but also what is plotted in which figure/subplot for single cell plots as well as colors!

%% KEY DEFINITION
%% overwriting keys with input
for fn=fieldnames(modified_keys)'
    keys.(fn{:})=modified_keys.(fn{:});
end

%% loading sorted neuron excel table
keys = get_sorting_table(keys);

%% folder creation
keys.path_to_save=[keys.basepath_to_save keys.project_version filesep];
if ~exist(keys.path_to_save,'dir')
    mkdir([keys.basepath_to_save], keys.project_version);
end
if ~exist([keys.path_to_save 'single_cells'],'dir')
    mkdir(keys.path_to_save, 'single_cells');
end
if ~exist([keys.path_to_save 'spike_shapes'],'dir')
    mkdir(keys.path_to_save, 'spike_shapes');
end

%% Get filelist
data_path                                       = [keys.drive filesep 'Data' filesep];
if isempty(keys.filelist_formatted)
    disp('get_filelist_from_folder is taking folder and dates from input')
    dates                                           = str2num(keys.date);
    [filelist_complete, filelist_formatted, filelist_session] = get_filelist_from_folder_run_block(keys,[data_path keys.monkey '_combined_monkeypsych_TDT'], dates);
    disp('Will analyze:')
else
    filelist_formatted=cell(0,3);
    filecounter=1;
    for ses=1:size(keys.filelist_formatted,1)
        n_files=numel(keys.filelist_formatted{ses,2});
        session_as_string=num2str(keys.filelist_formatted{ses,1});
        filelist_formatted(filecounter:filecounter+n_files-1,1)=repmat({session_as_string},n_files,1);
        blockorrun=run2block([keys.drive, 'Data', filesep keys.monkey '_combined_monkeypsych_TDT' filesep session_as_string],keys.filelist_formatted{ses,2},keys.filelist_as_blocks);
        if keys.filelist_as_blocks
            blocks=num2cell(keys.filelist_formatted{ses,2});
            runs=num2cell(blockorrun);
            bor='blocks';
        else
            runs=num2cell(keys.filelist_formatted{ses,2});
            blocks=num2cell(blockorrun);
            bor='runs';
        end
        filelist_formatted(filecounter:filecounter+n_files-1,2)=blocks;
        filelist_formatted(filecounter:filecounter+n_files-1,3)=runs;
        filecounter=filecounter+n_files;
    end
    to_remove=[filelist_formatted{:,2}]==0 | [filelist_formatted{:,2}]==0;
    filelist_formatted(to_remove,:)=[];
    
    filelist_session=filelist_formatted;
    for k=1:size(filelist_formatted,1)
        filelist_complete{k,1}=[keys.drive 'Data' filesep keys.monkey '_combined_monkeypsych_TDT' filesep filelist_formatted{k,1} filesep...
            keys.monkey(1:3) 'combined' filelist_formatted{k,1}(1:4) '-' filelist_formatted{k,1}(5:6) '-' filelist_formatted{k,1}(7:8) ...
            '_' sprintf('%02d',filelist_formatted{k,3}) '_block_' sprintf('%02d',filelist_formatted{k,2}) '.mat'];
    end
    disp(['Analyzing indivdually selected ' bor ' from sessions:'])
end
disp(filelist_complete)

%% MAIN LOOP
sessions=unique(filelist_session(:,1));
for current_date = sessions(:)'
    %% Preparing data for the full session
    keys.date       = current_date{1};
    session_indexes =ismember(filelist_session(:,1),current_date);
    blocks          =filelist_formatted(session_indexes,2);
    runs            =filelist_formatted(session_indexes,3);
    files           =filelist_complete(session_indexes);
    %% blockwise processing (actually it's filewise...!)
    for run_idx = 1:size(runs,1)
        keys.block=blocks{run_idx}; keys.run = runs{run_idx};
        Selection=keys.cal.MA_selection;
        MA_out=monkeypsych_analyze_working(files(run_idx),Selection); %% loading data including analyzed behaviour and physiology
        data_per_run(run_idx)=ph_run_state_alignment_per_trial(MA_out{1},keys); %
    end
    
    %% Sorting by unit/site/block ID & plotting waveform/spikesorting overview
    
    if keys.cal.process_sites
        sites = sort_by_site_ID(data_per_run);
    end
    if keys.cal.process_spikes
        pop_resorted = sort_by_unit_ID(data_per_run);
        traces_per_block = sort_traces_by_block(data_per_run);
    end
    if keys.cal.process_by_block
        by_block = sort_by_block(data_per_run);
    end
    clear data_per_run;
    
    if keys.cal.process_spikes
        pop_resorted(ismember({pop_resorted.unit_ID},{'no unit'}))=[];              %% remove units that were not taken over from excel sheet (??)
        pop_resorted = ph_accept_trials_per_unit(pop_resorted,keys);                %% add field accepted for each trial per unit
        
        % plot the cells not meeting criteria - not needed any more (?)
        idx_Neuron_ID=DAG_find_column_index(keys.sorting_table,'Neuron_ID');
        all_unit_IDs=keys.sorting_table_units(:,idx_Neuron_ID);
        if keys.plot.waveforms && any(~ismember({pop_resorted.unit_ID},all_unit_IDs))
            keys.path_to_save=[keys.basepath_to_save keys.project_version filesep 'spike_shapes' filesep];
            plot_sorted_waveforms(pop_resorted((~ismember({pop_resorted.unit_ID},all_unit_IDs))),keys,'units not found in sorting table');
            plot_across_time(pop_resorted((~ismember({pop_resorted.unit_ID},all_unit_IDs))),keys,'units not found in sorting table',ch_start_end,'FR');
        end
        
        % Excluding cells that are not found in the sortign table... i.e. pop_resorted.unit_ID wont be assigned
        if keys.cal.units_from_sorting_table
            pop_resorted(~ismember({pop_resorted.unit_ID},all_unit_IDs))=[];
        end
        
        % plot the cells used for further processing
        if keys.plot.waveforms && ~isempty(pop_resorted)
            keys.path_to_save=[keys.basepath_to_save keys.project_version filesep 'spike_shapes' filesep];
            %% go by channel 
            [channels csortidx]=sort([pop_resorted.channel]);
            channel_diffs=[find([0 diff(channels)]) numel(channels)];
            start_idx=0;
            plot_cross_correlation(pop_resorted,keys,'sorted units correlation');
            while start_idx < numel(channels)
                end_idx=channel_diffs(find(floor((channel_diffs-start_idx)/37)<1,1,'last'));
                if isempty(end_idx) % this is the case if theres more than 36 units in one channel!
                    end_idx=channel_diffs(1);
                end
                to_plot=csortidx(start_idx+1:end_idx);
                ch_start_end=['ch_' num2str(channels(start_idx+1)) '-' num2str(channels(end_idx))];
                plot_sorted_waveforms(pop_resorted(to_plot),keys,['sorted units, ' ch_start_end]);
                plot_sorted_ISI(pop_resorted(to_plot),keys,['sorted units ISI, ' ch_start_end]);
                plot_across_time(pop_resorted(to_plot),keys,'sorted units',ch_start_end,'FR');
                plot_across_time(pop_resorted(to_plot),keys,'sorted units',ch_start_end,'SNR');
                plot_across_time(pop_resorted(to_plot),keys,'sorted units',ch_start_end,'amp');
                plot_across_time(pop_resorted(to_plot),keys,'sorted units',ch_start_end,'noise');
                %plot_SNR_across_time(pop_resorted(to_plot),keys,['sorted units SNR over time, ' ch_start_end]);
                start_idx=find(channels==channels(end_idx),1,'last');
            end
        end
        
        if ~isempty(pop_resorted)
            [pop_resorted.monkey]=deal(keys.monkey);
            
            %% Save population mat file per session and output
            population=ph_reduce_population(pop_resorted);
            save([keys.population_foldername filesep keys.population_filename '_' current_date{1} '.mat'],'population');
            %% Save population mat file per session and output 
            save([keys.population_foldername filesep keys.traces_filename '_' current_date{1} '.mat'],'traces_per_block');
        end
    end
    
    if keys.cal.process_sites
        idx_Site_ID=DAG_find_column_index(keys.sorting_table,'Site_ID');
        all_site_IDs=keys.sorting_table_sites(:,idx_Site_ID);
        % Excluding sites that dont match criterions... i.e. sites.unit_ID wont be assigned
        if keys.cal.units_from_sorting_table
            sites(~ismember({sites.site_ID},all_site_IDs))=[];
        end
        if ~isempty(sites)      %% Save by site mat file per session
            save([keys.population_foldername filesep keys.sites_filename '_' current_date{1} '.mat'],'sites');
        end
    end
    if keys.cal.process_by_block && ~isempty(by_block)     %% Save by block mat file per session
        save([keys.population_foldername filesep 'by_block_' keys.monkey(1:end-5) '_' current_date{1} '.mat'],'by_block');
    end
end

end

%% data resorting functions
%% date block trial N ??

function pop_resorted = sort_by_site_ID(o_t)
pop_resorted=struct('site_ID',{});
site_index=0;
fields_to_remove={'unit','channel','x_eye','y_eye','x_hnd','y_hnd','time_axis','TDT_LFPx'...
    'TDT_CAP1','TDT_CAP1_SR','TDT_CAP1_tStart',...
    'TDT_ECG1','TDT_ECG1_SR','TDT_ECG1_tStart',...
    'TDT_ECG4','TDT_ECG4_SR','TDT_ECG4_tStart',...
    'TDT_POX1','TDT_POX1_SR','TDT_POX1_tStart'};

for b=1:size(o_t,2)
    for t=1:size(o_t(b).trial,2)
        n_chans_s = 0;
        if isfield(o_t(b).trial(t),'TDT_LFPx')
            n_chans_s =size(o_t(b).trial(t).TDT_LFPx,1);  % is this nonempty channels only?
        end
        for c=1:n_chans_s
            current_site_already_processed=ismember({pop_resorted.site_ID},{o_t(b).trial(t).channel(c).site_ID});
            
            %% take over channel information to trial structure directly
            o_t(b).trial(t).dataset                         =o_t(b).trial(t).channel(c).dataset;
            o_t(b).trial(t).perturbation                    =o_t(b).trial(t).channel(c).perturbation;
            o_t(b).trial(t).LFP                             =o_t(b).trial(t).TDT_LFPx(c,:);
            
            if any(current_site_already_processed)
                %% append to existing site
                s=find(current_site_already_processed);
                n_trial(s)                   =n_trial(s)+1;
            else
                %% create new site
                site_index=site_index+1;
                s=site_index;
                n_trial(s)=1;
                pop_resorted(s).site_ID          =o_t(b).trial(t).channel(c).site_ID;
                pop_resorted(s).target           =o_t(b).trial(t).channel(c).target;
                pop_resorted(s).perturbation_site=o_t(b).trial(t).channel(c).perturbation_site;
                pop_resorted(s).grid_x           =o_t(b).trial(t).channel(c).grid_x;
                pop_resorted(s).grid_y           =o_t(b).trial(t).channel(c).grid_y;
                pop_resorted(s).electrode_depth  =o_t(b).trial(t).channel(c).electrode_depth;
            end
            trial_fieldnames_to_remove=fields_to_remove(ismember(fields_to_remove,fieldnames(o_t(b).trial(t))));
            tmp=rmfield(o_t(b).trial(t),trial_fieldnames_to_remove);
            pop_resorted(s).trial(n_trial(s))=tmp;
        end
    end
end

end

function pop_resorted = sort_by_unit_ID_old(o_t)
pop_resorted=struct('unit_ID',{});
unit_index=0;
fields_to_remove={'unit','channel','TDT_CAP1','TDT_POX1','TDT_ECG1','TDT_ECG4','TDT_LFPx'};

for b=1:size(o_t,2)
    for t=1:size(o_t(b).trial,2)
        for c=1:size(o_t(b).trial(t).unit,1)
            for u=1:size(o_t(b).trial(t).unit,2)
                current_unit_already_processed=ismember({pop_resorted.unit_ID},{o_t(b).trial(t).unit(c,u).unit_ID});
                
                o_t(b).trial(t).waveforms                 =o_t(b).trial(t).unit(c,u).waveforms;
                o_t(b).trial(t).stability_rating          =o_t(b).trial(t).unit(c,u).stability_rating;
                o_t(b).trial(t).arrival_times             =o_t(b).trial(t).unit(c,u).arrival_times;
                o_t(b).trial(t).dataset                   =o_t(b).trial(t).unit(c,u).dataset;
                o_t(b).trial(t).perturbation              =o_t(b).trial(t).unit(c,u).perturbation;
                o_t(b).trial(t).FR_average                =o_t(b).trial(t).unit(c,u).FR_average;
                o_t(b).trial(t).stability_rating          =o_t(b).trial(t).unit(c,u).stability_rating;
                o_t(b).trial(t).SNR_rating                =o_t(b).trial(t).unit(c,u).SNR_rating;
                n_spikes                                  =size(o_t(b).trial(t).unit(c,u).waveforms,1);
                av_waveform=mean(o_t(b).trial(t).waveforms,1);
                std_waveform=std(o_t(b).trial(t).waveforms,1);
                o_t(b).trial(t).waveform_amplitude=max(av_waveform)-min(av_waveform);
                o_t(b).trial(t).waveform_std=nanmean(std_waveform);
                o_t(b).trial(t).trialwise_SNR=(max(av_waveform)-min(av_waveform))/nanmean(std_waveform);

                if any(current_unit_already_processed)
                    %% append to existing unit
                    U=find(current_unit_already_processed);
                    pop_resorted(U).n_spikes         =pop_resorted(U).n_spikes+n_spikes;
                    pop_resorted(U).SNR_rating       =[pop_resorted(U).SNR_rating o_t(b).trial(t).unit(c,u).SNR_rating*n_spikes];
                    pop_resorted(U).Single_rating    =[pop_resorted(U).Single_rating o_t(b).trial(t).unit(c,u).Single_rating*n_spikes];
                    pop_resorted(U).stability_rating =[pop_resorted(U).stability_rating o_t(b).trial(t).unit(c,u).stability_rating*n_spikes];
                    if ~ismember(num2str(o_t(b).block),[pop_resorted(U).block_unit{1,:}])
                        pop_resorted(U).block_unit       =[pop_resorted(U).block_unit,{num2str(o_t(b).block);char(96+u);' '}];
                    end
                    n_trial(U)                   =n_trial(U)+1;
                else
                    %% create new unit
                    unit_index=unit_index+1;
                    U=unit_index;
                    n_trial(U)=1;
                    pop_resorted(U).n_spikes      =n_spikes;
                    pop_resorted(U).channel          =c;
                    pop_resorted(U).block_unit       ={num2str(o_t(b).block);char(96+u);' '};
                    pop_resorted(U).unit_ID          =o_t(b).trial(t).unit(c,u).unit_ID;
                    pop_resorted(U).SNR_rating       =o_t(b).trial(t).unit(c,u).SNR_rating*n_spikes;
                    pop_resorted(U).Single_rating    =o_t(b).trial(t).unit(c,u).Single_rating*n_spikes;
                    pop_resorted(U).stability_rating =o_t(b).trial(t).unit(c,u).stability_rating*n_spikes;
                    pop_resorted(U).site_ID          =o_t(b).trial(t).unit(c,u).site_ID;
                    pop_resorted(U).target           =o_t(b).trial(t).unit(c,u).target;
                    pop_resorted(U).perturbation_site=o_t(b).trial(t).unit(c,u).perturbation_site;
                    pop_resorted(U).grid_x           =o_t(b).trial(t).unit(c,u).grid_x;
                    pop_resorted(U).grid_y           =o_t(b).trial(t).unit(c,u).grid_y;
                    pop_resorted(U).electrode_depth  =o_t(b).trial(t).unit(c,u).electrode_depth;
                end
                trial_fieldnames_to_remove=fields_to_remove(ismember(fields_to_remove,fieldnames(o_t(b).trial(t))));
                tmp=rmfield(o_t(b).trial(t),trial_fieldnames_to_remove);
                pop_resorted(U).trial(n_trial(U))=tmp;
            end
        end
    end
end

for u=1:numel(pop_resorted)
    pop_resorted(u).SNR_rating       = nansum(pop_resorted(u).SNR_rating)/pop_resorted(u).n_spikes;
    pop_resorted(u).Single_rating    = nansum(pop_resorted(u).Single_rating)/pop_resorted(u).n_spikes;
    pop_resorted(u).stability_rating = nansum(pop_resorted(u).stability_rating)/pop_resorted(u).n_spikes;
    
    pop_resorted(u).waveform_average = nanmean(vertcat(pop_resorted(u).trial.waveforms),1);
    pop_resorted(u).waveform_std     = nanstd(vertcat(pop_resorted(u).trial.waveforms),0,1);
    %% compute waveform_width (should work both for positive and negative spikes)
    resampling_factor=10;
    resampled_wf=resample(double(pop_resorted(u).waveform_average),resampling_factor,1);
    wf_minmax=[min(resampled_wf) max(resampled_wf)];
    [~,peaklocation]=max(abs(resampled_wf));
    peaksign=sign(resampled_wf(peaklocation));
    %t1
    t1=find(resampled_wf(1:peaklocation)*peaksign<(diff(wf_minmax)/2),1,'last');
    t2=peaklocation-1+find(resampled_wf(peaklocation:end)*peaksign<(max(abs(wf_minmax))-diff(wf_minmax)/2),1,'first');
    if ~isempty(wf_minmax) && all(~isnan(wf_minmax))
        %pop_resorted(u).waveform_width = sum(abs(pop_resorted(u).waveform_average-pop_resorted(u).waveform_average(end)))/24414.0625/diff(wf_minmax); %% sampling rate hardcoded here
        pop_resorted(u).waveform_width = (t2-t1)/24414.0625/resampling_factor; %% sampling rate hardcoded here
    else
        pop_resorted(u).waveform_width = -1; 
    end
    %% compute Signal-to-noice ratio for one unit 
    pop_resorted(u).waveform_amplitude = diff(wf_minmax);  
    %pop_resorted(u).quantSNR =  pop_resorted(u).waveform_amplitude /nanmean(pop_resorted(u).waveform_std); 
end
end

function pop_resorted = sort_by_unit_ID(o_t)
pop_resorted=struct('unit_ID',{},'trial',{},'waveforms',{},'SNR_rating',{},'Single_rating',{},'stability_rating',{},'block_unit',{});
fields_to_remove={'unit','channel','TDT_CAP1','TDT_POX1','TDT_ECG1','TDT_ECG4','TDT_LFPx'};
for b=1:size(o_t,2)
    for u=1:numel(o_t(b).trial(1).unit)
        [c, un]                     = ind2sub(size(o_t(b).trial(1).unit),u);
        
        current_unit                = arrayfun(@(x) x.unit(u),o_t(b).trial);
        ID=current_unit.unit_ID;
        if strcmp(ID,'no unit')
           continue; 
        end
        U=find(ismember({pop_resorted.unit_ID},ID),1);
        if isempty(U)
            U=numel(pop_resorted)+1;
            pop_resorted(U).n_spikes         =0;
            pop_resorted(U).unit_ID          =ID;
            pop_resorted(U).channel     = c;
            pop_resorted(U).site_ID          =current_unit(1).site_ID;
            pop_resorted(U).target           =current_unit(1).target;
            pop_resorted(U).perturbation_site=current_unit(1).perturbation_site;
            pop_resorted(U).grid_x           =current_unit(1).grid_x;
            pop_resorted(U).grid_y           =current_unit(1).grid_y;
            pop_resorted(U).electrode_depth  =current_unit(1).electrode_depth;
        end
        pop_resorted(U).block_unit  = [pop_resorted(U).block_unit,{num2str(o_t(b).block);char(96+un);' '}];
        n_spikes=size(vertcat(current_unit.waveforms),1);
        pop_resorted(U).n_spikes         =pop_resorted(U).n_spikes + n_spikes;
        pop_resorted(U).SNR_rating       =[pop_resorted(U).SNR_rating       current_unit(1).SNR_rating*n_spikes];
        pop_resorted(U).Single_rating    =[pop_resorted(U).Single_rating    current_unit(1).Single_rating*n_spikes];
        pop_resorted(U).stability_rating =[pop_resorted(U).stability_rating current_unit(1).stability_rating*n_spikes];      
        
        trial_fieldnames_to_remove=fields_to_remove(ismember(fields_to_remove,fieldnames(o_t(b).trial)));
        tmp=rmfield(o_t(b).trial,trial_fieldnames_to_remove);
        [tmp.arrival_times]=current_unit.arrival_times;
        [tmp.waveforms]=current_unit.waveforms;
        [tmp.dataset]=current_unit.dataset;
        [tmp.perturbation]=current_unit.perturbation;
        [tmp.FR_average]=current_unit.FR_average;
        [tmp.stability_rating]=current_unit.stability_rating;
        [tmp.SNR_rating]=current_unit.SNR_rating;
        
        amp_wf=arrayfun(@(x) max(mean(x.waveforms,1)-min(mean(x.waveforms,1))),current_unit,'uniformoutput',false);
        std_wf=arrayfun(@(x) nanmean(std(x.waveforms,1)),current_unit,'uniformoutput',false);
        snr_wf=cellfun(@(x,y) x/y,amp_wf,std_wf,'uniformoutput',false);
        [tmp.waveform_amplitude]=deal(amp_wf{:});
        [tmp.waveform_std]      =deal(std_wf{:});
        [tmp.trialwise_SNR]     =deal(snr_wf{:});
        pop_resorted(U).trial   =[pop_resorted(U).trial tmp];
    end
end

    
for u=1:numel(pop_resorted)
    pop_resorted(u).SNR_rating       = nansum(pop_resorted(u).SNR_rating)/pop_resorted(u).n_spikes;
    pop_resorted(u).Single_rating    = nansum(pop_resorted(u).Single_rating)/pop_resorted(u).n_spikes;
    pop_resorted(u).stability_rating = nansum(pop_resorted(u).stability_rating)/pop_resorted(u).n_spikes;    
    pop_resorted(u).waveform_average = nanmean(vertcat(pop_resorted(u).trial.waveforms),1);
    pop_resorted(u).waveform_std     = nanstd(vertcat(pop_resorted(u).trial.waveforms),0,1);
    
    %% compute waveform_width (should work both for positive and negative spikes)
    resampling_factor=10;
    resampled_wf=resample(double(pop_resorted(u).waveform_average),resampling_factor,1);
    wf_minmax=[min(resampled_wf) max(resampled_wf)];
    [~,peaklocation]=max(abs(resampled_wf));
    peaksign=sign(resampled_wf(peaklocation));
    %t1
    t1=find(resampled_wf(1:peaklocation)*peaksign<(diff(wf_minmax)/2),1,'last');
    t2=peaklocation-1+find(resampled_wf(peaklocation:end)*peaksign<(max(abs(wf_minmax))-diff(wf_minmax)/2),1,'first');
    if ~isempty(wf_minmax) && all(~isnan(wf_minmax))
        %pop_resorted(u).waveform_width = sum(abs(pop_resorted(u).waveform_average-pop_resorted(u).waveform_average(end)))/24414.0625/diff(wf_minmax); %% sampling rate hardcoded here
        pop_resorted(u).waveform_width = (t2-t1)/24414.0625/resampling_factor; %% sampling rate hardcoded here
    else
        pop_resorted(u).waveform_width = -1;
    end
    pop_resorted(u).waveform_amplitude = diff(wf_minmax);
end
end


function pop_resorted = sort_by_block(o_t)
fields_to_remove={'unit','channel','x_eye','y_eye','x_hnd','y_hnd','time_axis',...
    'TDT_LFPx','TDT_LFPx_SR','TDT_LFPx_tStart'};
% need to fix same block issues.. (-_-)
block=0;
for b=1:size(o_t,2)
    if block~=o_t(b).block % new block
      t_start=0;
    end
    block=o_t(b).block;
    for t=1:size(o_t(b).trial,2)
        %o_t(b).trial(t).block=b;
        o_t(b).trial(t).perturbation                    =o_t(b).trial(t).channel(1).perturbation; %%??
        trial_fieldnames_to_remove=fields_to_remove(ismember(fields_to_remove,fieldnames(o_t(b).trial(t))));
        tmp=rmfield(o_t(b).trial(t),trial_fieldnames_to_remove);
        pop_resorted(block).trial(t_start+t)=tmp;
    end
    t_start=t_start+t;
end
end

function pop_resorted = sort_traces_by_block(o_t)
fields_to_remove={'unit','channel','TDT_LFPx','TDT_LFPx_SR','TDT_LFPx_tStart','TDT_CAP1','TDT_POX1','TDT_ECG1','TDT_ECG4'};
for b=1:size(o_t,2)
    for t=1:size(o_t(b).trial,2)
        trial_fieldnames_to_remove=fields_to_remove(ismember(fields_to_remove,fieldnames(o_t(b).trial(t))));
        tmp=rmfield(o_t(b).trial(t),trial_fieldnames_to_remove);
        pop_resorted(b).trial(t)=tmp;
    end
end
end

function pop=ph_reduce_population(pop_in)
pop=pop_in;
fields_to_remove={'waveforms','x_eye','y_eye','x_hnd','y_hnd','time_axis'};
for u=1:numel(pop)
    pop(u).trial=rmfield(pop(u).trial,fields_to_remove);
end
end


%% data path functions

function [filelist_complete, filelist_formatted, filelist_session] = get_filelist_from_folder_run_block(keys,main_folder,dates)
dir_main=dir(main_folder); % dir
session_folders=[];
ctr=1;
for k=1: length(dir_main)
    X=str2double(dir_main(k).name);
    if X<=  dates(2) && X >=  dates(1)
        session_folders{ctr}= dir_main(k).name;
        ctr=ctr+1;
    end
end

r=1; % run or row counter
for s = 1:length(session_folders)
    session_folder = [main_folder filesep session_folders{s}]; % session of interest
    d_individual_day_folder=dir(session_folder); % dir
    session_files={d_individual_day_folder.name}'; % files inside session folders
    session_files=session_files(cellfun(@(x) length(x) > 4 && strcmp(x(end-3:end),'.mat'),session_files));
    for f = 1:length(session_files) % start looping within the session folder
        filelist_complete{r,1}=[session_folder filesep session_files{f}];
        filesepindx=findstr(filelist_complete{r,1}, filesep);
        filelist_formatted(r,:)= {filelist_complete{r,1}(1:filesepindx(end)-1) str2double(filelist_complete{r,1}(end-5:end-4)) str2double(filelist_complete{r,1}(end-14:end-13))};
        filelist_session(r,:)= {filelist_complete{r,1}(filesepindx(end-1)+1:filesepindx(end)-1) filelist_complete{r,1}(end-5:end-4) filelist_complete{r,1}(end-14:end-13)};
        filelist_internal(r,:)= [str2double(filelist_complete{r,1}(filesepindx(end-1)+1:filesepindx(end)-1)) str2double(filelist_complete{r,1}(end-5:end-4)) str2double(filelist_complete{r,1}(end-14:end-13))];
        r=r+1;
    end
end

if ~isempty(keys.cal.datasets)
    idx_dataset=DAG_find_column_index(keys.sorting_table,'Set');
    idx_session=DAG_find_column_index(keys.sorting_table,'Date');
    idx_run=DAG_find_column_index(keys.sorting_table,'Run');
    idx_block=DAG_find_column_index(keys.sorting_table,'Block');
    rows=[false; ismember(cell2mat(keys.sorting_table(2:end,idx_dataset)),keys.cal.datasets)];
    filelist_fitting_datasets=cell2mat(keys.sorting_table(rows,[idx_session  idx_block  idx_run]));
    index_in_dataset=ismember(filelist_internal,filelist_fitting_datasets,'rows');
    filelist_complete=filelist_complete(index_in_dataset,:);
    filelist_formatted=filelist_formatted(index_in_dataset,:);
    filelist_session=filelist_session(index_in_dataset,:);
end
end

function output_blocks_or_runs=run2block(base_path,given_blocks_or_runs,reverse)
folder_content=dir([base_path filesep '*.mat']);
complete_filelist=vertcat(folder_content.name);
filelist_blocks=str2num(complete_filelist(:,end-5:end-4));
filelist_runs=str2num(complete_filelist(:,end-14:end-13));

if strcmp(reverse,'reverse') || reverse==1
    [~,idx] =ismember(given_blocks_or_runs,filelist_blocks);
    output_blocks_or_runs(idx~=0)=filelist_runs(idx(idx~=0));
else
    [~,idx] =ismember(given_blocks_or_runs,filelist_runs);
    output_blocks_or_runs(idx~=0)=filelist_blocks(idx(idx~=0));
end
output_blocks_or_runs(idx==0)=0;
end

function keys = get_sorting_table(keys)
keys.sorting_table_units = {};
keys.sorting_table_sites = {};
keys.sorting_table       = {};
temp_xlsx=dir(fullfile(keys.sorted_neurons_foldername,[keys.sorted_neurons_filename '.xls*']));
if ~isempty(temp_xlsx)
    excel_sheet= [keys.sorted_neurons_foldername filesep temp_xlsx(1).name];
    [~,~,sorting_table] = xlsread(excel_sheet,keys.sorted_neurons_sheetname); %sorting_table.NUMM doesnt work properly, because whole columns were missing
    
    keys.sorting_table_units = sorting_table;
    keys.sorting_table_sites = sorting_table;
    keys.sorting_table       = sorting_table;
    usable_index=DAG_find_column_index(sorting_table,'Usable');
    to_exclude_s=~ismember([sorting_table{2:end,usable_index}]',1); % think about other site criterias and maybe we want an option to include not usable?
%     keys.sorting_table_units([false;to_exclude_u],:) = [];
     keys.sorting_table_sites([false;to_exclude_s],:) = [];
end
xlswrite([keys.tuning_table_foldername filesep keys.sorted_neurons_filename],sorting_table);
end


%% plot functions

function plot_across_time(o,keys,which_units,ch_start_end,whattoplot)
title_part=[which_units ' ' whattoplot ' over time, ' ch_start_end];
fig_title=sprintf('%s, session %s, %s',keys.monkey,keys.date,title_part);
FR_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',fig_title);

n_columns_rows=ceil(numel(o)^(1/2));
alltrialz=[o.trial];
firstbin=min([alltrialz.run_onset_time]);
[lastbin,lasttrial_idx]=max([alltrialz.run_onset_time]+[alltrialz.trial_onset_time]);
lastbin=lastbin+max(alltrialz(lasttrial_idx).states_onset);

units=1:numel(o);
for u=units
    subplot(n_columns_rows,n_columns_rows,u);
    hold on;
    binsize=60;
    AT=[];
    WF=[];
    for t=1:numel(o(u).trial)
        ATt=o(u).trial(t).arrival_times;
        WFt=o(u).trial(t).waveforms;
        AT=vertcat(AT,ATt(ATt>0 & ATt<o(u).trial(t).states_onset(end-1))+o(u).trial(t).trial_onset_time+o(u).trial(t).run_onset_time-firstbin);
        WF=vertcat(WF,WFt);%(ATt>0 & ATt<o(u).trial(t).states_onset(end-1),:));
    end
    %bins=(o(u).trial(1).trial_onset_time):binsize:(o(u).trial(end).run_onset_time-o(u).trial(1).run_onset_time+o(u).trial(end).trial_onset_time+max(o(u).trial(end).states_onset));
    bins=0:binsize:(lastbin-firstbin);
    
    if ismember(whattoplot,{'SNR','amp','noise'})
        snr=NaN(size(bins));
        amp=NaN(size(bins));
        noi=NaN(size(bins));
        
        for b=1:numel(bins)
            idx=AT>bins(b)-binsize/2 & AT<bins(b)+binsize/2;
            meanwf=mean(WF(idx,:),1);
            noi(b)=mean(std(WF(idx,:),1));
            amp(b)=abs(max(meanwf)-min(meanwf));
            snr(b)=amp(b)/noi(b);
        end
        
    end
    switch whattoplot
        case 'FR'
            toplot=hist(AT,bins')/binsize;
            toplot_per_trial=[o(u).trial.FR_average];
            toplot_per_trial(isnan(toplot_per_trial))=0;
        case 'SNR'
            toplot=snr;
            toplot_per_trial=[o(u).trial.SNR_rating];
        case 'amp'
            toplot=amp;
            toplot_per_trial=zeros(numel(o(u).trial),1);
        case 'noise'
            toplot=noi;
            toplot_per_trial=zeros(numel(o(u).trial),1);
    end
    plot(bins,toplot);
    y_lim=ylim(gca);
    trial_blocks=[o(u).trial.block];
    trial_stability=[o(u).trial.stability_rating];
    unique_blocks=unique(trial_blocks);
    for b=unique_blocks
        tr_idx=trial_blocks==b & ~isnan(trial_stability);
        if sum(tr_idx)<2; continue; end;            % it can happen that an entire block is not accepted if FR changed drastically
        %FR_std=double(nanstd(FR_smoothed(tr_idx)));
        block_mean=double(nanmean(toplot_per_trial(tr_idx)));
        start_block=o(u).trial(find(tr_idx,1,'first')).run_onset_time-firstbin+o(u).trial(find(tr_idx,1,'first')).trial_onset_time;
        end_block=start_block+o(u).trial(find(tr_idx,1,'last')).trial_onset_time-o(u).trial(find(tr_idx,1,'first')).trial_onset_time;
        fanoish_factor=trial_stability(tr_idx);fanoish_factor=fanoish_factor(1);
        if fanoish_factor > 5 %% replace with keys
            col='g';
        elseif fanoish_factor> 2.5
            col='b';
        else
            col='r';
        end
        plot([start_block end_block],[block_mean block_mean],col,'linewidth',2)
        plot([start_block start_block],[0 block_mean],col,'linewidth',2)
        plot([end_block end_block],[0 block_mean],col,'linewidth',2)
        if strcmp(whattoplot,'FR')
            text(double(start_block+(end_block-start_block)/2), diff(y_lim)/2,sprintf('%0.1f',fanoish_factor),'HorizontalAlignment', 'Center')
        end
        %plot()
    end
    unit_title={sprintf('%s ',o(u).unit_ID),...
        sprintf(['ch/De: %d/%.2f b&u: %s' ],o(u).channel,o(u).electrode_depth,[o(u).block_unit{:}])}; %MP add number of spikes
    title(unit_title,'interpreter','none');
end
ph_title_and_save(FR_summary_handle,fig_title,fig_title,keys)
end

function plot_cross_correlation(population,keys,title_part)
chs=unique([population.channel]);
n_rows=ceil(sqrt(numel(chs)));
n_columns=ceil(sqrt(numel(chs)));
FR_summary_handle=figure;

fig_title=sprintf('%s, session %s, %s',keys.monkey,keys.date,title_part);
for c=1:numel(chs)
    CRL=[];
    ch=chs(c);
    subplot(n_rows,n_columns,c)
    pop=population([population.channel]==ch);
    unit_ticks=zeros(size(pop));
    for u1=1:numel(pop)
        u2=0;
        unit_ticks(u1)=str2double(pop(u1).unit_ID(end-1:end));
        while u2<u1
            u2=u2+1;
            u1_runs=unique([pop(u1).trial.run]);
            u2_runs=unique([pop(u2).trial.run]);
            u1_blocks=unique([pop(u1).trial.block]);
            u2_blocks=unique([pop(u2).trial.block]);
            r=intersect(u1_runs,u2_runs);
            b=intersect(u1_blocks,u2_blocks);
            if ~isempty(r)
                r1=ismember([pop(u1).trial.run],r) & ismember([pop(u1).trial.block],b);
                r2=ismember([pop(u2).trial.run],r) & ismember([pop(u2).trial.block],b);
                coeff=corr([[pop(u1).trial(r1).FR_average]',[pop(u2).trial(r2).FR_average]']);
                CRL(u1,u2)=coeff(1,2);
            else
                CRL(u1,u2)=0;
            end
        end
        
    end
    CRL=CRL+CRL';
    CRL(CRL>1)=1;
    
    X=1:numel(unit_ticks);
    imagesc(X,X,CRL,[0 1])
    set(gca,'ytick',X);
    set(gca,'yticklabel',unit_ticks);
    set(gca,'xtick',X);
    set(gca,'xticklabel',unit_ticks);
    xlabel('unit ID');
    ylabel('unit ID');
    title(['channel ' num2str(ch)]);
    colormap('jet')
    cb = colorbar;
    set(get(cb,'title'),'string', 'corr trial FR', 'fontsize',8);
end
ph_title_and_save(FR_summary_handle,fig_title,fig_title,keys)
close(gcf);

end

function plot_sorted_waveforms(o,keys,title_part)
fig_title=sprintf('%s, session %s, %s',keys.monkey,keys.date,title_part);
WF_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',fig_title);
for n_unit=1:numel(o)
    n_columns_rows=ceil(numel(o)^(1/2));
    subplot(n_columns_rows,n_columns_rows,n_unit)
    block_trials=[find([true diff([o(n_unit).trial(1:end-1).block])~=0]) numel(o(n_unit).trial)];
    wf_per_block=[];
    for b=1:numel(block_trials)-1
        meanblockwf=nanmean(cat(1,o(n_unit).trial(block_trials(b):block_trials(b+1)).waveforms),1);
        if ~isempty(meanblockwf);  wf_per_block(b,:)=meanblockwf; end
    end
    all_spikes_wf = cat(1,o(n_unit).trial.waveforms);
    n_all_spikes_wf = size(all_spikes_wf,1);
    
    % REDUCING WFS DISPLAYED
    n_sel_spike_wf = length(1:50:n_all_spikes_wf);
    if n_sel_spike_wf>0
        set(gca,'ColorOrder',jet(n_sel_spike_wf)),hold on
        plot(all_spikes_wf(1:50:n_sel_spike_wf*50,:)');
    end
    plot(wf_per_block','-k','linewidth',2);
    unit_title={sprintf('%s SN/Si/St: %.1f/%.1f/%.1f',o(n_unit).unit_ID,o(n_unit).SNR_rating,o(n_unit).Single_rating,o(n_unit).stability_rating),...
        sprintf(['spk: %d ch/De: %d/%.2f b&u: %s'],n_all_spikes_wf, o(n_unit).channel,o(n_unit).electrode_depth,[o(n_unit).block_unit{:}])}; %MP add number of spikes
        %sprintf('SNR/Wam/Std: %%d/%d/%d',round(o(n_unit).quantSNR),round(o(n_unit).waveform_amplitude),round(nanmean(o(n_unit).waveform_std))),... %%KK stuff
    title(unit_title,'interpreter','none','fontsize',6)
    if  max(max(all_spikes_wf(1:50:n_sel_spike_wf*50,:))) > min(min(all_spikes_wf(1:50:n_sel_spike_wf*50,:))) % not sure what this bug is about... Lin 20160303 
    set(gca,'xtick',[],'xcolor',[1 1 1],'FontSize',6,'ytick',...
        [min(min(all_spikes_wf(1:50:n_sel_spike_wf*50,:)')) max(max(all_spikes_wf(1:50:n_sel_spike_wf*50,:)'))]);          %MP remove X-axis keep Y-axis to have scale and show max/min values on Y axis
    end
end
ph_title_and_save(WF_summary_handle,fig_title,fig_title,keys)
end

function plot_sorted_ISI(o,keys,title_part)
fig_title=sprintf('%s, session %s, %s',keys.monkey,keys.date,title_part);
ISI_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',fig_title);

x_bins=logspace(-3,0,30);
x_bins=horzcat(0,x_bins);
%ISI = struct();
for n_unit=1:numel(o)
    n_columns_rows=ceil(numel(o)^(1/2));
    subplot(n_columns_rows,n_columns_rows,n_unit)
    %ISI(n_unit).trial(1).isi(1)=NaN;
    AT=NaN;
    for t=1:numel(o(n_unit).trial)
        AT=[AT o(n_unit).trial(t).arrival_times'+o(n_unit).trial(t).trial_onset_time];
        %         for n_spike_in_trial = 1:numel(o(n_unit).trial(t).arrival_times)
        %             if t==1 && n_spike_in_trial==1
        %                 continue;
        %             elseif n_spike_in_trial==1
        %                ISI(n_unit).trial(t).isi(n_spike_in_trial)=o(n_unit).trial(t).arrival_times(n_spike_in_trial) - o(n_unit).trial(t-1).arrival_times(end) + ...
        %                                                           o(n_unit).trial(t).trial_onset_time - o(n_unit).trial(t-1).trial_onset_time ;
        %             else
        %             ISI(n_unit).trial(t).isi(n_spike_in_trial)=o(n_unit).trial(t).arrival_times(n_spike_in_trial) - o(n_unit).trial(t).arrival_times(n_spike_in_trial-1);
        %             end
        %         end
    end
    AT=unique(AT); % due to ovrelapping end and beginning of trial, spikes can be counted twice
    all_ISI = diff(AT); %cat(2,ISI(n_unit).trial.isi);
    
    if ~isempty(all_ISI)
        hist_values = histc(all_ISI,x_bins);
        perc_first_bin = (hist_values(1)/sum(hist_values))*100;
        %         histogram(all_ISI,x_bins,'FaceColor','r','EdgeColor','w');
        %         hist(all_ISI,x_bins)
        %         bar(log10(x_bins(2:end)),hist_values(2:end));
        x_bins_bar=logspace(-3,0,31);
        %         bar(log10(x_bins(2:end-1)), hist_values(2:end-1))
        bar(log10(x_bins_bar), hist_values)
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor','r','EdgeColor','w');
        %         x_bins_ticks = [10^-3 10^-2 10^-1 0.5];
        x_bins_ticks = [-3:1];
        x_bins_ticks_label = [10^-3 10^-2 10^-1 0.5];
        set(gca,'ylim',[0 max([1 hist_values])],'xtick',x_bins_ticks,'box','off');%max(hist(all_ISI,x_bins)),'ycolor',[1 1 1],'xscale','log'
        set(gca,'XTickLabel',sprintf('%2.0e|',x_bins_ticks_label));
        %         h2 = findobj(gca,'Type','line');
        %         set(h2,'Marker','none');
        %         xlim([0,0.7]);
        text(-0.45,max(hist_values),[num2str(perc_first_bin,'%4.1f') '%<1ms'],'FontSize',8)
        
    end
    unit_title={sprintf('%s SN/Si/St: %d/%d/%d',o(n_unit).unit_ID,o(n_unit).SNR_rating,o(n_unit).Single_rating,o(n_unit).stability_rating)...
        sprintf(['spk: %d ch/De: %d/%.2f b: ' num2str(unique([o(n_unit).trial.block]))],numel(all_ISI), o(n_unit).channel,o(n_unit).electrode_depth)}; %MP add number of spikes
    title(unit_title,'interpreter','none','fontsize',10)
    
    
end
ph_title_and_save(ISI_summary_handle,fig_title,fig_title,keys)
end

