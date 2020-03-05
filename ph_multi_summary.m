function  [population_internal, keys] = ph_multi_summary(modified_keys)

global MA_STATES
%% KEY DEFINITION

%% Path keys
keys.drive                  ='X:';
keys.monkey                 ='Aureus_phys';
keys.date                   =datestr(date,'yyyymmdd');
keys.filelist_formatted     ={};
% keys.filelist_formatted={'20160203',[4,5,7];...
%                          '20160205',[3,4,5];...
%                          '20160205',[2];...
%                          '20160210',[8]}
keys.filelist_as_blocks     =0;
keys.run                    =[];
keys.block                  =[];

keys.tuning_table_filename      ='';
keys.sorted_neurons_filename    ='';
keys.population_filename        ='population';
keys.tuning_table_foldername    ='';
keys.sorted_neurons_foldername  ='';
keys.population_foldername      ='';

keys.basepath_to_save       ='Projects\Pulv_microstim_behavior\ephys';
keys.pdf_folder             ='ephys_analysis_v3_June2016';

%% General processing
%run 'epoch_definitions'; %% getting all epoch related keys from here
keys.exclude_channels       =[];
keys.MA_selection           ={};
keys.effectors              =[0,1,2,3,4,5,6];
keys.types                  =[1,2,3,4,5,6];
keys.reach_hand             =[0,1,2];

keys.save_preprocessed_data=0;
keys.case_summaries={'hands'};
keys.process_only_units_in_the_table=0;
keys.minimum_trials_per_condition=1;
keys.min_n_spikes_per_unit=50;
keys.FR_at_peak=0;
keys.FR_subtract_baseline=0;

%% Pdf management
keys.create_pdfs            =0;
keys.append_pdfs            =0;
keys.delete_previous_pdfs   =0;

%% Plotting
keys.effectors_on_same_figure   =1;
keys.plot_trials                =0;    % not used yet
keys.plot_vertical_positons_PSTH=0; %not used yet
keys.plot_average_line      =0;
keys.plot_FR_separately     =1;
keys.average_heat_maps      =0;

keys.plot                   =1;
keys.plot_waveforms         =1;
keys.plot_eye_hand_traces   =1; %% Incredible performance booster if turned off
keys.summary                =[1 2];
keys.PSTH_binwidth          =0.01;
keys.unsmoothed_histogram   =0;
keys.uname                  ={'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j'};
keys.display_anova_tables   ='off';

%% ylimits
keys.max_firing_rate_expected   =50;
keys.max_n_trials_expected      =50;
keys.max_excentricity_expected  =30;
keys.eye_offset_multiplier      =[110 0.5];
keys.hnd_offset_multiplier      =[140 0.5];

%% linewidth and colors

keys.l_width_his_s1         =0.5;
keys.l_width_ras_s1         =0.1;
keys.l_width_his_s2         =1;
keys.eye_ver_color          =[0.8 0 0];
keys.eye_hor_color          =[1 0 0];
keys.rhnd_ver_color         =[0 0.8 0];
keys.rhnd_hor_color         =[0 1 0];
keys.lhnd_ver_color         =[0 0 0.8];
keys.lhnd_hor_color         =[0 0 1];
keys.offset_colors          =[0 0 1; 0 1 1; 0 1 0];
% keys.hnd_choice_colors      =[1 0 0.5; 1 0.5 0; 0 0 1; 0 0 0.4; 0 1 0; 0 0.5 0]; %NH IN, NH CH, LH IN, LH CH, RH IN, RH CH;
% keys.hnd_choice_colors_L    =[1 0 0.5; 0.6 0 0.3; 0 0 1; 0 0 0.4; 1 0.7 0; 0.5 0.5 0]; %NH IN, NH CH, LH IN, LH CH, RH IN, RH CH;
% keys.hnd_choice_colors_R    =[0.9 0.6 0.25; 0.36 0.24 0.1; 0.8 0.5 1; 0.4 0 0.6; 0 1 0; 0 0.5 0]; %NH IN, NH CH, LH IN, LH CH, RH IN, RH CH;
% keys.color_average_line     =[0 0 0];

keys.hnd_choice_colors      =[1 0 0; 0.7 0 0; 0 0 1; 0 0 0.4; 0 1 0; 0 0.5 0]; %NH IN, NH CH, LH IN, LH CH, RH IN, RH CH;
keys.hnd_choice_colors_L    =[255 0 178; 127 0 89; 171 0 252; 85 0 126; 145 143 56; 77 76 28]/255; %NH IN, NH CH, LH IN, LH CH, RH IN, RH CH;
keys.hnd_choice_colors_R    =[255 153 20; 127 77 10; 125 130 255; 62 65 127; 222 220 0; 110 110 0]/255; %NH IN, NH CH, LH IN, LH CH, RH IN, RH CH;
keys.color_average_line     =[0 0 0];

%% overwriting keys with input
for fn=fieldnames(modified_keys)'
    keys.(fn{:})=modified_keys.(fn{:});
end

%% loading sorted neuron excel table
keys.xlsx_table = get_excel_table(keys);

%% folder creation
keys.path_to_save=[keys.drive filesep keys.basepath_to_save filesep keys.pdf_folder filesep];
if ~exist(keys.path_to_save,'dir')
    mkdir([keys.drive filesep keys.basepath_to_save], keys.pdf_folder);
end
if ~exist([keys.path_to_save filesep 'single_cell_examples'],'dir')
    mkdir(keys.path_to_save, 'single_cell_examples');
end
if ~exist([keys.path_to_save filesep 'spike_shapes'],'dir')
    mkdir(keys.path_to_save, 'spike_shapes');
end

%% Get filelist
data_path                                       = [keys.drive '\Data\'];
if isempty(keys.filelist_formatted)%(isempty(keys.block) && isempty(keys.run))
    dates                                           = str2num(keys.date);
    if size(dates,2) ~= 2
        dates                                           = repmat(dates,1,2);
    end
    
    [filelist_complete, filelist_formatted, filelist_session] = get_filelist_from_folder_run_block([data_path keys.monkey '_combined_monkeypsych_TDT'], dates);
    disp('Will analyze:')
    disp(filelist_complete)
    
    if keys.delete_previous_pdfs
        pdf_dir_to_delete   = unique(filelist_formatted(:,1));
        for idx_del         = 1:numel(pdf_dir_to_delete)
            delete([pdf_dir_to_delete{idx_del} filesep '*.pdf'])
        end
    end
else
    
    filelist_formatted=cell(0,3);
    filecounter=1;
    for ses=1:size(keys.filelist_formatted,1)
        n_files=numel(keys.filelist_formatted{ses,2});
        session_as_string=num2str(keys.filelist_formatted{ses,1});
        filelist_formatted(filecounter:filecounter+n_files-1,1)=repmat({session_as_string},n_files,1);
        blockorrun=run2block([keys.drive, filesep, 'Data', filesep keys.monkey '_combined_monkeypsych_TDT' filesep session_as_string],keys.filelist_formatted{ses,2},keys.filelist_as_blocks);
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
    disp(['Analyzing indivdually selected ' bor ' from sessions:'])
    disp(['    session      block   run'])
    disp(filelist_formatted)
    %     else
    %     dates                                           = str2num(keys.date);
    %
    %     if size(dates,2) ~= 2
    %         dates                                           = repmat(dates,1,2);
    %     end
    %     [filelist_complete_temp, filelist_formatted_temp, filelist_session_temp] = get_filelist_from_folder_run_block([data_path keys.monkey '_combined_monkeypsych_TDT'], dates);
    %     if ~isempty(keys.block) && isempty(keys.run)
    %         for i=1:size(filelist_formatted_temp,1)
    %             if filelist_formatted_temp{i,2} == keys.block
    %                 keys.run=filelist_formatted_temp{i,3};
    %             end
    %         end
    %     end
    %     if ~isempty(keys.run) && isempty(keys.block)
    %         for i=1:size(filelist_formatted_temp,1)
    %             if filelist_formatted_temp{i,3} == keys.run
    %                 keys.block=filelist_formatted_temp{i,2};
    %             end
    %         end
    %     end
    %     disp('Will analyze:')
    %     disp([keys.monkey ' ' keys.date ' run ' sprintf('%02d',keys.run) ' block ' sprintf('%02d',keys.block)])
    %     filelist_session={keys.date, num2str(keys.block), num2str(keys.run)};
    %     filelist_formatted={keys.date, keys.block, keys.run};
end

%% MAIN LOOP
sessions=unique(filelist_session(:,1));
% n_population=0;
n_existing_files=0;
for current_date = sessions(:)'
    
    n_population=0;
    o_a=[];
    
    %% Preparing data for the full session
    session_indexes =ismember(filelist_session(:,1),current_date);
    blocks          =filelist_formatted(session_indexes,2);
    runs            =filelist_formatted(session_indexes,3);
    
    %% blockwise processing
    for blo_idx = 1:size(blocks,1)
        keys.block=blocks{blo_idx};
        keys.date =current_date{1};    keys.run = runs{blo_idx};    keys.block = blocks{blo_idx};
        Selection=[{'display',0,'keep_raw_data',1,'saccade_definition',4,'reach_1st_pos',1} keys.MA_selection];
        keys.runpath2=[keys.drive, filesep, 'Data', filesep keys.monkey '_combined_monkeypsych_TDT' filesep keys.date];
        MA_out=monkeypsych_analyze_working({keys.runpath2,keys.block},Selection); %% loading data including analyzed behaviour and physiology
        o(blo_idx)=ph_run_state_alignment_per_trial(MA_out{1},keys); %
    end
    
    %% Sorting by unit ID & saving preprocessed output
    o_resorted = sort_by_unit_ID(o,keys);
    clear o;
    
    o_resorted(ismember({o_resorted.unit_ID},{'no unit'}))=[];
    idx_Neuron_ID=find_column_index(keys.xlsx_table,'Neuron_ID');
    all_unit_IDs=keys.xlsx_table(:,idx_Neuron_ID);
    
    if keys.plot_waveforms && any(~ismember({o_resorted.unit_ID},all_unit_IDs))
        plot_sorted_by_unit_ID(o_resorted((~ismember({o_resorted.unit_ID},all_unit_IDs)) | ([o_resorted.n_waveforms]<keys.min_n_spikes_per_unit)),'skipped units',keys);
    end
    o_resorted([o_resorted.n_waveforms]<keys.min_n_spikes_per_unit)=[];
    
    % Excluding units without ID or not
    if keys.process_only_units_in_the_table
        o_resorted(~ismember({o_resorted.unit_ID},all_unit_IDs))=[];
    end
    
    if keys.plot_waveforms
        plot_sorted_by_unit_ID(o_resorted,'analyzed units',keys);
    end
    
    %% Restructuring by conditions and displacement types! + excluding not selected effectors/types/hands ?
    for unit=1:numel(o_resorted)
        effectors=[o_resorted(unit).trial.effector];
        types=[o_resorted(unit).trial.type];
        hands=[o_resorted(unit).trial.reach_hand];
        effector_loop_idx=unique(effectors);
        type_loop_idx=unique(types);
        effector_loop_idx=effector_loop_idx(ismember(effector_loop_idx,keys.effectors));
        type_loop_idx=type_loop_idx(ismember(type_loop_idx,keys.types));
        for summary = 1:numel(keys.case_summaries);
            keys.plot_type = keys.case_summaries{summary};
            condition_index=0;
            for eff=effector_loop_idx
                for typ=type_loop_idx
                    disp(['unit: ' o_resorted(unit).unit_ID ' Current case: ', keys.case_summaries{summary}, ' ', get_type_effector_name(typ,eff)])
                    tr_index= effectors==eff & types == typ & ismember(hands,keys.reach_hand);
                    if keys.list_successful_only
                        tr_index = tr_index & [o_resorted(unit).trial.success]==1;
                    end
                    if sum(tr_index & [o_resorted(unit).trial.choice]==0) < keys.minimum_trials_per_condition
                        continue;
                    end
                    condition_index=condition_index+1;
                    o_c.trial=o_resorted(unit).trial(tr_index);
                    hands_on=any([o_resorted(unit).trial.reach_hand]>0);
                    
                    get_expected_MA_states(typ,eff,keys.effectors_on_same_figure);
                    keys.ALL_EPOCHS=keys.EPOCHS_PER_TYPE{typ};
                    keys.EPOCHS=keys.ALL_EPOCHS(ismember([keys.ALL_EPOCHS{:,2}],MA_STATES.all_states),:);
                    if keys.FR_at_peak %%%% to check !!!
                        [o_c.trial,keys]=ph_FR_at_peak(o_c.trial,keys);
                    end
                    if keys.FR_subtract_baseline %%%% to check !!!
                        [o_c.trial,keys]=ph_FR_subtract_baseline(o_c.trial,keys);
                    end
                    
                    [~, displacement_types] = center_displacement_working(o_c.trial);
                    
                    %% ph_restructure_by_displacement_types! source of most potential problems
                    o_temp=ph_restructure_by_displacement_types(o_c.trial,displacement_types,keys); % keys output for n columns
                    FN={'position','left','right','up','down'};
                    [o_a(unit).condition(condition_index).case(summary)]=add_structure_per_states(o_temp,FN,keys);
                    o_a(unit).condition(condition_index).effector=eff;
                    o_a(unit).condition(condition_index).type=typ;
                    o_a(unit).condition(condition_index).hands_on=hands_on;
                    o_a(unit).unit_ID=o_resorted(unit).unit_ID;
                    [o_a(unit).unit_ID o_a(unit).stability_rating o_a(unit).SNR_rating o_a(unit).Single_rating o_a(unit).target o_a(unit).grid_x o_a(unit).grid_y o_a(unit).electrode_depth o_a(unit).channel o_a(unit).block_unit] =...
                        deal(o_resorted(unit).unit_ID, o_resorted(unit).stability_rating, o_resorted(unit).SNR_rating, o_resorted(unit).Single_rating, o_resorted(unit).target, o_resorted(unit).grid_x, o_resorted(unit).grid_y, o_resorted(unit).electrode_depth, o_resorted(unit).channel, o_resorted(unit).block_unit);
                end
            end
        end
    end
    clear o_resorted
    
    %% statistics and table output
    if ~isempty(o_a); o_a(arrayfun(@(x) isempty(x.unit_ID),o_a))=[]; end;
    if ~isempty(o_a)
        if exist([keys.tuning_table_foldername filesep keys.tuning_table_filename '.mat'],'file')
            load([keys.tuning_table_foldername filesep keys.tuning_table_filename '.mat']);
            keys.tuning_per_unit_table=tuning_per_unit_table;
        else
            keys.tuning_per_unit_table= {'unit_ID'};
        end
        clear tuning_per_unit_table
        tuning_per_unit_table=ph_analyze_one_session(o_a,keys); % main function
        save([keys.tuning_table_foldername filesep keys.tuning_table_filename '.mat'],'tuning_per_unit_table');
    end
    population_internal=struct([]);
    if ~isempty(o_a)
        population_internal=ph_population(o_a);
        %population_internal(n_population+1:n_population+numel(o_a))=ph_population(o_a);
        n_population=numel(population_internal);
    end
    
    %% Plotting
    if keys.plot==1
        for unit=1:numel(o_a)
            [keys.unit_ID keys.stability_rating keys.SNR_rating keys.Single_rating keys.target keys.grid_x keys.grid_y keys.electrode_depth keys.channel keys.block_unit] =...
                deal(o_a(unit).unit_ID, o_a(unit).stability_rating, o_a(unit).SNR_rating, o_a(unit).Single_rating, o_a(unit).target, o_a(unit).grid_x, o_a(unit).grid_y, o_a(unit).electrode_depth, o_a(unit).channel, o_a(unit).block_unit);
            if keys.effectors_on_same_figure;
                [unique_types, ~, type_index]=unique([o_a(unit).condition.type]);
            else
                unique_types=1:numel(o_a(unit).condition);
                type_index= unique_types;
            end
            for condition=1:numel(unique_types)
                condition_index=type_index==condition;
                typ=unique([o_a(unit).condition(condition_index).type]);
                eff=[o_a(unit).condition(condition_index).effector];
                hands_on=any([o_a(unit).condition(condition_index).hands_on]);
                get_expected_MA_states(typ,eff,keys.effectors_on_same_figure);
                keys.ALL_EPOCHS=keys.EPOCHS_PER_TYPE{typ};
                keys.EPOCHS=keys.ALL_EPOCHS(ismember([keys.ALL_EPOCHS{:,2}],MA_STATES.all_states),:);
                [keys.effector, keys.type, keys.hands_on] = deal(eff, typ, hands_on);
                keys.current_unit_tuning= [tuning_per_unit_table(1,:); tuning_per_unit_table(ismember(tuning_per_unit_table(:,1), o_a(unit).unit_ID),:)];
                to_plot_temp=o_a(unit).condition(condition_index);
                to_plot_temp=vertcat(to_plot_temp.case);
                for summary = 1:size(to_plot_temp,2);
                    keys.plot_type=keys.case_summaries{summary}; %% not correct either
                    disp(['unit: ' o_a(unit).unit_ID ' Current case: ', keys.plot_type]);
                    %                     disp(['unit: ' o_a(unit).unit_ID ' Current case: ', keys.plot_type, ' ', get_type_effector_name(typ,eff)]);
                    ph_plot_unit_per_condition([to_plot_temp(:,summary)],keys);
                    get_expected_MA_states(typ,eff,keys.effectors_on_same_figure);
                end
            end
        end
    end
    
    
    %% Save population mat file
    if ~isempty(population_internal)
        n_files=ceil(numel(population_internal)/100);
        for f=1:n_files
            if f==n_files
                population=population_internal((f-1)*100+1:end);
            else
                population=population_internal((f-1)*100+1:f*100);
            end
            save([keys.population_foldername filesep keys.population_filename '_' num2str(f+n_existing_files,'%02u') '.mat'],'population');
        end
        n_existing_files=n_existing_files+n_files;
    end
    
end

end

function o_resorted = sort_by_unit_ID(o_t,keys)
clear o_resorted
o_resorted.unit_ID='';
unit_index=0;
for b=1:size(o_t,2)
    for t=1:size(o_t(b).trial,2)
        o_t(b).trial(t).block=b;
        for c=1:size(o_t(b).trial(t).unit,1)
            for u=1:size(o_t(b).trial(t).unit,2)
                current_unit_already_processed=ismember({o_resorted.unit_ID},{o_t(b).trial(t).unit(c,u).unit_ID});
                if any(current_unit_already_processed)
                    % append to existing unit
                    unit_index_processed=find(current_unit_already_processed);
                    if ~ismember(num2str(o_t(b).block),[o_resorted(unit_index_processed).block_unit{1,:}])
                        o_resorted(unit_index_processed).block_unit     =[o_resorted(unit_index_processed).block_unit,{num2str(o_t(b).block);keys.uname{u};' '}];
                        o_resorted(unit_index_processed).SNR_rating       =[o_resorted(unit_index_processed).SNR_rating o_t(b).trial(t).unit(c,u).SNR_rating];
                        o_resorted(unit_index_processed).Single_rating    =[o_resorted(unit_index_processed).Single_rating o_t(b).trial(t).unit(c,u).Single_rating];
                        o_resorted(unit_index_processed).stability_rating =[o_resorted(unit_index_processed).stability_rating o_t(b).trial(t).unit(c,u).stability_rating];
                    end
                    n_trial(unit_index_processed)                   =n_trial(unit_index_processed)+1;
                    o_t(b).trial(t).waveforms                       =o_t(b).trial(t).unit(c,u).waveforms;
                    o_t(b).trial(t).spikes_per_state                =o_t(b).trial(t).unit(c,u).spikes_per_state;
                    tmp=rmfield(o_t(b).trial(t),'unit');
                    o_resorted(unit_index_processed).trial(n_trial(unit_index_processed))=tmp;
                    o_resorted(unit_index_processed).n_waveforms =o_resorted(unit_index_processed).n_waveforms + size(o_t(b).trial(t).unit(c,u).waveforms,1);
                else
                    %create new unit
                    unit_index=unit_index+1;
                    n_trial(unit_index)=1;
                    o_resorted(unit_index).unit_ID          =o_t(b).trial(t).unit(c,u).unit_ID;
                    o_resorted(unit_index).channel          =c;
                    o_resorted(unit_index).SNR_rating       =o_t(b).trial(t).unit(c,u).SNR_rating;
                    o_resorted(unit_index).Single_rating    =o_t(b).trial(t).unit(c,u).Single_rating;
                    o_resorted(unit_index).stability_rating =o_t(b).trial(t).unit(c,u).stability_rating;
                    o_resorted(unit_index).target           =o_t(b).trial(t).unit(c,u).target;
                    o_resorted(unit_index).grid_x           =o_t(b).trial(t).unit(c,u).grid_x;
                    o_resorted(unit_index).grid_y           =o_t(b).trial(t).unit(c,u).grid_y;
                    o_resorted(unit_index).electrode_depth  =o_t(b).trial(t).unit(c,u).electrode_depth;
                    o_resorted(unit_index).block_unit       ={num2str(o_t(b).block);keys.uname{u};' '};
                    o_t(b).trial(t).waveforms               =o_t(b).trial(t).unit(c,u).waveforms;
                    o_t(b).trial(t).spikes_per_state        =o_t(b).trial(t).unit(c,u).spikes_per_state;
                    tmp=rmfield(o_t(b).trial(t),'unit');
                    o_resorted(unit_index).trial(n_trial(unit_index))=tmp;
                    o_resorted(unit_index).n_waveforms      =size(o_t(b).trial(t).unit(c,u).waveforms,1);
                end
            end
        end
    end
end

for u=1:numel(o_resorted)
    o_resorted(u).SNR_rating       = round(mean(o_resorted(u).SNR_rating)*10)/10;
    o_resorted(u).Single_rating    = round(mean(o_resorted(u).Single_rating)*10)/10;
    o_resorted(u).stability_rating = round(mean(o_resorted(u).stability_rating)*10)/10;
end
end

function plot_sorted_by_unit_ID(o_resorted,title_part,keys)
fig_title=sprintf('%s, session %s, %s',keys.monkey,keys.date,title_part);
WF_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',fig_title);

for n_unit=1:numel(o_resorted)
    n_columns_rows=ceil(numel(o_resorted)^(1/2));
    subplot(n_columns_rows,n_columns_rows,n_unit)
    block_trials=find([true diff([o_resorted(n_unit).trial(1:end-1).block])~=0 true]);
    wf_per_block=[];
    for b=1:numel(block_trials)-1
        meanblockwf=mean(cat(1,o_resorted(n_unit).trial(block_trials(b):block_trials(b+1)).waveforms));
        if ~isempty(meanblockwf);  wf_per_block(b,:)=meanblockwf; end
    end
    all_spikes_wf = cat(1,o_resorted(n_unit).trial.waveforms);
    n_all_spikes_wf = size(all_spikes_wf,1);
    
    % REDUCING WFS DISPLAYED
    n_sel_spike_wf = length(1:50:n_all_spikes_wf);
    if n_sel_spike_wf>0
        set(gca,'ColorOrder',jet(n_sel_spike_wf)),hold on
        plot(all_spikes_wf(1:50:n_sel_spike_wf*50,:)');
    end
    plot(wf_per_block','-k','linewidth',2);
    unit_tile=sprintf('%s ch/De: %d/%.2f SN/Si/St: %d/%d/%d',o_resorted(n_unit).unit_ID,o_resorted(n_unit).channel,...
        o_resorted(n_unit).electrode_depth,o_resorted(n_unit).SNR_rating,o_resorted(n_unit).Single_rating,o_resorted(n_unit).stability_rating);
    title(unit_tile,'interpreter','none','fontsize',8)
    axis off
end
title_and_save(WF_summary_handle,fig_title,keys)
end

function [filelist_complete, filelist_formatted, filelist_session] = get_filelist_from_folder_run_block(varargin)

if numel(varargin)==0
    folder_with_session_days=transfer_parameters.folders.TDT_data;
    dates=transfer_parameters.dates;
    disp('get_filelist_from_folder is taking folder and dates from transfer_parameters')
else
    folder_with_session_days=varargin{1};
    dates=varargin{2};
    disp('get_filelist_from_folder is taking folder and dates from input')
end

dir_folder_with_session_days=dir(folder_with_session_days); % dir
if all(isnan(dates)) || all(isempty(dates)) || all(dates==0)
    all_files_in_base_path=keep_only_numeric_cell_entries({dir_folder_with_session_days.name}); % all files from the main monkey folder
else
    all_files_in_base_path=[];
    ctr=1;
    for k=1: length(dir_folder_with_session_days)
        X=str2double(dir_folder_with_session_days(k).name);
        if X==dates(1) ||  ( X<=  dates(2) && X >  dates(1))
            all_files_in_base_path{ctr}= dir_folder_with_session_days(k).name;
            ctr=ctr+1;
        end
    end
    
end

i_run=1;
for in_folders = 1:length(all_files_in_base_path)
    individual_day_folder = [folder_with_session_days filesep all_files_in_base_path{in_folders}]; % session of interest
    d_individual_day_folder=dir(individual_day_folder); % dir
    files_inside_session_folder={d_individual_day_folder.name}'; % files inside session folders
    for number_of_files = 1:length(files_inside_session_folder) % start looping within the session folder
        if length(files_inside_session_folder{number_of_files}) > 4 && strcmp(files_inside_session_folder{number_of_files}(end-3:end),'.mat')
            filelist_complete(i_run,:)=[individual_day_folder filesep files_inside_session_folder{number_of_files}];
            filesepindx=findstr(filelist_complete(i_run,:), filesep);
            filelist_formatted(i_run,:)= {filelist_complete(i_run,1:filesepindx(end)-1) str2double(filelist_complete(i_run,end-5:end-4)) str2double(filelist_complete(i_run,end-14:end-13))};
            filelist_session(i_run,:)= {filelist_complete(i_run,filesepindx(end-1)+1:filesepindx(end)-1) filelist_complete(i_run,end-5:end-4) filelist_complete(i_run,end-14:end-13)};
            i_run=i_run+1;
        end
    end
    
end
end

function xlsx_table = get_excel_table(keys)
xlsx_table={};
temp_xlsx=dir(fullfile(keys.sorted_neurons_foldername,[keys.sorted_neurons_filename '*.xls*']));
if ~isempty(temp_xlsx)
    excel_sheet= [keys.sorted_neurons_foldername filesep temp_xlsx.name];
    [~,~,xlsx_table] = xlsread(excel_sheet); %xlsx_table.NUMM doesnt work properly, because whole columns were missing
    
    stability_index=find_column_index(xlsx_table,'Stability rank');
    single_index=find_column_index(xlsx_table,'Single rank');
    to_exclude=~ismember([xlsx_table{2:end,stability_index}]',keys.cal.stablity) | ~ismember([xlsx_table{2:end,single_index}]',keys.cal.single_rating);
    
    %xlsx_table.NUMM(to_exclude,:)=[];
    xlsx_table([false;to_exclude],:)=[];
    %xlsx_table([false;to_exclude],:)=[];
    
    
    %xlsx_table.masternum_orig=[NaN(1,size(xlsx_table.NUMM,2));xlsx_table.NUMM];
end
xlswrite([keys.tuning_table_foldername filesep keys.sorted_neurons_filename],xlsx_table);
end

function title_and_save(figure_handle,plot_title,keys)
mtit(figure_handle,  plot_title, 'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none');
stampit;
if keys.create_pdfs
    switch keys.append_pdfs
        case 0
            export_fig([keys.path_to_save 'spike_shapes' filesep plot_title], '-pdf','-transparent') % pdf by run
        case 1
            export_fig([keys.path_to_save 'spike_shapes' filesep '_appended'], '-pdf', '-append','-transparent') % pdf by run
    end
    close all
end
end

function out_analysis_extended=add_structure_per_states(out_analysis,FN,keys)
out_analysis_extended=out_analysis;
for fn=1:numel(FN)
    
    if ~isfield(out_analysis.(FN{fn}),'trial')
        out_analysis_extended.(FN{fn}).trial=struct('arrival_times',{});%critical!?
    end
    for idx=1:numel(out_analysis.(FN{fn}))
        for sta=1:size(keys.EPOCHS,1)
            n=0;
            out_analysis_extended.(FN{fn})(idx).per_state(sta).state=keys.EPOCHS{sta,1}; %%%%%!!!!!!!!!!!!!!!!!!!!!!
            if ~isfield(out_analysis.(FN{fn})(idx),'trial') || numel(out_analysis.(FN{fn})(idx).trial)==0
                out_analysis_extended.(FN{fn})(idx).per_state(sta).trial=struct('arrival_times',{});%critical!?
                continue;
            end
            for tri=1:numel(out_analysis.(FN{fn})(idx).trial)
                if numel(out_analysis.(FN{fn})(idx).trial(tri).spikes_per_state)>=sta && ~isempty(out_analysis.(FN{fn})(idx).trial(tri).spikes_per_state(sta).state)
                    n=n+1;
                    out_analysis_extended.(FN{fn})(idx).per_state(sta).trial(n)=out_analysis.(FN{fn})(idx).trial(tri).spikes_per_state(sta);
                end
            end
        end
    end
end
end

function [s_c, displacement_types] = center_displacement_working(trial)
Precision=2;

movement_direction  =NaN(size(trial'));
fixation            =NaN(size(trial'));
target              =NaN(size(trial'));

s_a=unique_positions([trial.fix_pos],Precision);
s_b=unique_positions([trial.tar_pos] - [trial.fix_pos],Precision);
s_c=unique_positions([trial.tar_pos],Precision);

for t=1:numel(trial)
    for k=1:numel(s_a)
        if abs(trial(t).fix_pos - s_a(k)) < Precision
            fixation(t)=s_a(k);
        end
    end
    for k=1:numel(s_b)
        if abs(trial(t).tar_pos - trial(t).fix_pos - s_b(k)) < Precision
            movement_direction(t)=k;
        end
    end
    for k=1:numel(s_c)
        if abs(trial(t).tar_pos - s_c(k)) < Precision
            target(t)=s_c(k);
        end
    end
end

[~,~,unique_condition]      =unique([real(fixation),imag(fixation),real(target),imag(target)],'rows');
[~,~,fixation_location]     =unique([real(fixation),imag(fixation)],'rows');
[~,~,movement_direction]    =unique([real(movement_direction),imag(movement_direction)],'rows');
[~,~,target_location]       =unique([real(target),imag(target)],'rows');

fix_y=imag(nanmean(fixation));
fixation=fixation-1i*fix_y;
target=target-1i*fix_y;

displacement_types=[real(fixation) imag(fixation) real(target) imag(target) unique_condition fixation_location movement_direction target_location];
end

function unique_target_positions=unique_positions(all_target_positions,Precision)
target_positions    =   unique(all_target_positions(~isnan(all_target_positions)));
n_targets           =   numel(target_positions);
for t=1:n_targets
    target_positions(abs(target_positions-target_positions(t))<Precision)=target_positions(t);
end
unique_target_positions=unique(target_positions);
end

function [oo,keys]=ph_FR_subtract_baseline(o,keys)
oo=o;
for t=1:numel(o)
    for s=1:numel(o(t).spikes_per_state)
        b=ismember(keys.EPOCHS(:,1),keys.EPOCHS(s,7));
        oo(t).spikes_per_state(s).FR=o(t).spikes_per_state(s).FR - o(t).spikes_per_state(b).FR;
    end
end
end

function [oo,keys]=ph_FR_at_peak(o,keys)
oo=o;
for t=1:numel(o)
    for s=1:size(keys.EPOCHS,1)
        restructured_by_states(s).trial(t)=o(t).spikes_per_state(s);
        restructured_by_states(s).state=o(t).spikes_per_state(s).state;
    end
end

for s=1:numel(restructured_by_states)
    state=restructured_by_states(s).state;
    state_idx=ismember(keys.EPOCHS(:,1),state);
    t_before_state=keys.EPOCHS{state_idx,3};
    t_after_state=keys.EPOCHS{state_idx,4};
    bins=t_before_state+keys.PSTH_binwidth/2:keys.PSTH_binwidth:t_after_state-keys.PSTH_binwidth/2;
    histo=hist(vertcat(restructured_by_states(s).trial.arrival_times),bins);
    histo=smooth(histo-histo(1),1,1.5);
    t_max=find(abs(histo)==max(abs(histo)));
    peak_location(s)=bins(t_max(1));
    keys.EPOCHS{state_idx,3}=peak_location(s)-keys.PSTH_binwidth/2;
    keys.EPOCHS{state_idx,4}=peak_location(s)+keys.PSTH_binwidth/2;
end
for t=1:numel(o)
    for s=1:numel(o(t).spikes_per_state)
        at=o(t).spikes_per_state(s).arrival_times;
        arrival_times=at(at>peak_location(s)-keys.PSTH_binwidth & at<peak_location(s)+keys.PSTH_binwidth);
        oo(t).spikes_per_state(s).FR=numel(arrival_times)/keys.PSTH_binwidth;
    end
end
end

function pop=ph_population(o_a)

FN={'left','right','position'};
pop=rmfield(o_a,{'channel','block_unit'});
fields_to_remove1={'PSTH_colors','PSTH_summary_colors','line_labels','up','down','figure_title_part','figure_title_value'};
fields_to_remove2={'trial'};
fields_to_remove3={'state','arrival_times','x_eye','y_eye','x_hnd','y_hnd','time_axis','line','figure'};

for u=1:numel(pop)
    for c=1:numel(pop(u).condition)
        pop(u).condition(c).case=rmfield(pop(u).condition(c).case,fields_to_remove1);
        for a=1:numel(pop(u).condition(c).case)
            for fn=1:numel(FN)
                pop(u).condition(c).case(a).(FN{fn})=rmfield(pop(u).condition(c).case(a).(FN{fn}),fields_to_remove2);
                for p=1:numel(pop(u).condition(c).case(a).(FN{fn}))
                    for s=1:numel(pop(u).condition(c).case(a).(FN{fn})(p).per_state)
                        % Cornelius 20160506, unit 1 condition 1 case 1 position 2
                        if isfield(pop(u).condition(c).case(a).(FN{fn})(p).per_state(s).trial,'x_eye')
                            pop(u).condition(c).case(a).(FN{fn})(p).per_state(s).trial=...
                                rmfield(pop(u).condition(c).case(a).(FN{fn})(p).per_state(s).trial,fields_to_remove3);
                        end
                    end
                end
            end
        end
    end
end
end