function o=ph_run_state_alignment_per_trial(MA_out,keys)
global MA_STATES
%% Defining maximum number of units
n_trials=numel(MA_out.physiology);
tr_in=[MA_out.physiology];
n_chans_u = max(arrayfun(@(x) size(x.spike_arrival_times,1),tr_in));  % is this nonempty channels only?
n_chans_s = 1;
if isfield(tr_in,'TDT_LFPx')
    nonempty=find(arrayfun(@(x) ~isempty(x.TDT_LFPx),tr_in));
    n_chans_s = size(tr_in(nonempty(1)).TDT_LFPx,1);  % first nonempty trial, is this nonempty channels only?
end
%n_chans=max([n_chans_u n_chans_s]);
n_units=0;
%% unitsperchannelmatrix JUST FOR DEBUGGING PURPOSES
unitsperchannelmatrix=false(n_chans_u,1);
for t=1:n_trials
    if ~isempty(tr_in(t).spike_arrival_times)
        existent_in_this_trial=~cellfun(@isempty,tr_in(t).spike_arrival_times);
        n_this_trial=size(existent_in_this_trial,2);
        n_accumulated=size(unitsperchannelmatrix,2);
        unitsperchannelmatrix=[unitsperchannelmatrix | existent_in_this_trial(:,1:n_accumulated) existent_in_this_trial(:,n_accumulated+1:n_this_trial) ];
        n_units = max([n_units size(tr_in(t).spike_arrival_times,2)]);
    end
end
for c=1:n_chans_s
    from_excel_per_channel(c) = get_sorted_neuron_out_xlsx_from_excel(keys.sorting_table_sites, keys, 1, c, 1,0);
end
for c=1:n_chans_u
    for u=1:n_units
        from_excel_per_unit(c,u) = get_sorted_neuron_out_xlsx_from_excel(keys.sorting_table_units, keys, u, c, unitsperchannelmatrix(c,u),1);
    end
end

%% Preallocation !?
NaNpar=num2cell(NaN(1,n_trials));
trial=struct('type',NaNpar,'effector',NaNpar,'reach_hand',NaNpar,'choice',NaNpar,'success',NaNpar,'fix_pos',NaNpar,'tar_pos',NaNpar,'trial_onset_time',NaNpar,...
    'time_axis',NaNpar,'x_eye',NaNpar,'y_eye',NaNpar,'x_hnd',NaNpar,'y_hnd',NaNpar,'states_onset',NaNpar); %add all fields of trial(?)
invalid_trials=[];

%% main loop
for t=1:n_trials
    trial_states=MA_out.states(t).TDT_states;
    trial_states_onset=MA_out.states(t).TDT_state_onsets;
    %% what to do with last trial? is the last state 1?
    %% unit Lin_20150521_02 debug
    %%W:\Data\Linus_phys_combined_monkeypsych_TDT\20150521\Lincombined2015-05-21_03_block_01.mat
    if isempty(trial_states_onset) || (t==n_trials && trial_states(end)~=1) || (t==n_trials-1 && trial_states(end)~=1) || (trial_states(end)~=90 && trial_states(end)~=99 && trial_states(end)~=1)%% Curius 20150603 block 4 trial 166
        %&& (t==n_trials-1 || t==n_trials) % last trial bug that should get fixed in TDT_trial_struct_working  line 30
        invalid_trials=[invalid_trials t];
        continue;
    end
    trial(t).date          = str2num(keys.date);
    trial(t).block         = keys.block;
    trial(t).run           = MA_out.selected(t).run;
    trial(t).n             = MA_out.selected(t).trials;
    trial(t).type          = MA_out.task(t).type;
    trial(t).effector      = MA_out.task(t).effector;
    trial(t).reach_hand    = MA_out.task(t).reach_hand;
    if isnan(trial(t).reach_hand); trial(t).reach_hand=0; end;
    trial(t).choice        = MA_out.binary(t).choice;
    trial(t).success       = MA_out.binary(t).success;
    trial(t).completed     = MA_out.binary(t).completed;
    
    
    %% Saccade or reach positions dependend on effector
    MA_out.reaches(t).ini_all=MA_out.reaches(t).ini;
    MA_out.reaches(t).end_all=MA_out.reaches(t).ini+MA_out.reaches(t).dur;
    MA_out.reaches(t).endpos_all=MA_out.reaches(t).endpos;
    MA_out.reaches(t).startpos_all=MA_out.reaches(t).startpos;
    if ismember(trial(t).effector,[0,3])
        Movement=MA_out.saccades(t);
    else
        Movement=MA_out.reaches(t);
    end
    trial(t).sac_off     = MA_out.saccades(t).endpos-MA_out.saccades(t).tar_pos;
    trial(t).sac_lat     = MA_out.saccades(t).lat;
    trial(t).rea_off     = MA_out.reaches(t).endpos-MA_out.reaches(t).tar_pos;
    trial(t).rea_lat     = MA_out.reaches(t).lat;
    
    trial(t).fix_pos     = Movement.fix_pos;
    trial(t).tar_pos     = Movement.tar_pos;
    trial(t).cue_pos     = Movement.cue_pos;
    trial(t).cue_shape   = Movement.all_convexities(1);
    
    % for distractor task    
    trial(t).all_tar_pos        = Movement.all_tar_pos;
    trial(t).col_dim            = Movement.col_dim;
    trial(t).col_bri            = Movement.col_bri;
    trial(t).target_selected    = Movement.target_selected;
    
    sac_ini = MA_out.states(t).start_obs + MA_out.saccades(t).lat;
    sac_end = sac_ini + MA_out.saccades(t).dur;
    rea_ini = MA_out.states(t).start_obs + MA_out.reaches(t).lat;
    rea_end = rea_ini + MA_out.reaches(t).dur;
    tri_end = MA_out.states(t).start_end;
    
    MPA_get_expected_states(trial(t).type,trial(t).effector,0);    %% set to 1 to allow later processing
    movement_states       = [MA_STATES.SAC_INI MA_STATES.SAC_END MA_STATES.REA_INI MA_STATES.REA_END MA_STATES.TRI_END];
    movement_onsets       = [sac_ini sac_end rea_ini rea_end tri_end];
    movement_onsets       = movement_onsets(ismember(movement_states,MA_STATES.all_states));
    movement_states       = movement_states(ismember(movement_states,MA_STATES.all_states));
    
    trial_states            =[trial_states(1:end-1) movement_states 98]; %MA_STATES.ITI ?
    trial_states_onset      =[trial_states_onset(1:end-1) movement_onsets trial_states_onset(end)];
    [~,tr_state_idx]        =unique(trial_states,'last');
    trial(t).states         =trial_states(sort(tr_state_idx));
    trial(t).states_onset   =trial_states_onset(sort(tr_state_idx));
    
    %% excluding unwanted trials from the beginning... this is actually sort of problematic
    if ~ismember(trial(t).type,keys.cal.types) || ~ismember(trial(t).effector,keys.cal.effectors) || ~ismember(trial(t).reach_hand,keys.cal.reach_hand)
        invalid_trials=[invalid_trials t];
        continue;
    end
    
    %     %% for several movements per trial
    %     mov_idx=~isnan(Movement.ini_all) & Movement.ini_all>=MA_out.states(t).start_obs;
    %     trial(t).movement_onsets        =Movement.ini_all(mov_idx);
    %     trial(t).movement_ends          =Movement.end_all(mov_idx);
    %     trial(t).movement_endpos        =Movement.endpos_all(mov_idx);
    %     trial(t).movement_startpos      =Movement.startpos_all(mov_idx);
    %
    %     state_indexes_start=arrayfun(@(x) find(MA_out.states(t).MP_states_onset>=x,1,'first')-1, trial(t).movement_onsets);
    %     state_indexes_end=arrayfun(@(x) find(MA_out.states(t).MP_states_onset<=x,1,'last'), trial(t).movement_ends);
    %
    %     trial(t).movement_tar_at_start    =MA_out.states(t).MP_states(state_indexes_start)==MA_out.states(t).state_1ao;
    %     trial(t).movement_tar_at_end      =MA_out.states(t).MP_states(state_indexes_end)==MA_out.states(t).state_1ao;
    %     trial(t).movement_tar_crossed     =state_indexes_end - state_indexes_start>2; %% works so far, because saccade can only be crossing one target if there are two state transitions within the saccade;
    %
    %     tar_ins=arrayfun(@(x) sum(MA_out.states(t).MP_states(1:x)==MA_out.states(t).state_1ao), state_indexes_start);
    %     tar_ins(tar_ins==0)=1;
    %     trial(t).movement_shape_at_start    = Movement.all_convexities(Movement.targets_inspected(tar_ins));
    %     trial(t).movement_shape_crossed    = Movement.all_convexities(Movement.targets_inspected(tar_ins));
    %     tar_ins=arrayfun(@(x) sum(MA_out.states(t).MP_states(1:x)==MA_out.states(t).state_1ao), state_indexes_end);
    %     tar_ins(tar_ins==0)=1;
    %     trial(t).movement_shape_at_end    = Movement.all_convexities(Movement.targets_inspected(tar_ins));
    %
    %     trial(t).movement_shape_at_start(~trial(t).movement_tar_at_start)=NaN;
    %     trial(t).movement_shape_crossed(~trial(t).movement_tar_crossed)=NaN;
    %     trial(t).movement_shape_at_end(~trial(t).movement_tar_at_end)=NaN;
    
    
    %% adding previous trial
    trial(t).run_onset_time              =MA_out.states(t).run_onset_time;
    trial(t).trial_onset_time        =MA_out.states(t).trial_onset_time;
    
    if t>1 % && keys.add_previous_trial_spikes     shift_in_seconds=1;
        trial(t).time_axis = [MA_out.raw(t-1).time_axis-(trial(t).trial_onset_time -MA_out.states(t-1).trial_onset_time) MA_out.raw(t).time_axis];
        trial(t).x_eye = [MA_out.raw(t-1).x_eye MA_out.raw(t).x_eye];
        trial(t).y_eye = [MA_out.raw(t-1).y_eye MA_out.raw(t).y_eye];
        trial(t).x_hnd = [MA_out.raw(t-1).x_hnd MA_out.raw(t).x_hnd];
        trial(t).y_hnd = [MA_out.raw(t-1).y_hnd MA_out.raw(t).y_hnd];
    else
        trial(t).time_axis = MA_out.raw(t).time_axis;
        trial(t).x_eye = MA_out.raw(t).x_eye;
        trial(t).y_eye = MA_out.raw(t).y_eye;
        trial(t).x_hnd = MA_out.raw(t).x_hnd;
        trial(t).y_hnd = MA_out.raw(t).y_hnd;
    end
    %% resample for speed !!!
    trace_idx=1:MA_out.keys.calc_val.i_sample_rate*keys.PSTH_binwidth:numel(trial(t).time_axis);
    trial(t).time_axis  = trial(t).time_axis(trace_idx);
    trial(t).x_eye      = trial(t).x_eye(trace_idx);
    trial(t).y_eye      = trial(t).y_eye(trace_idx);
    trial(t).x_hnd      = trial(t).x_hnd(trace_idx);
    trial(t).y_hnd      = trial(t).y_hnd(trace_idx);
    
    
    %% LFPs & other streams??
    shift_in_seconds=1;
    stream_fieldnames=fieldnames(tr_in);
    stream_fieldnames=stream_fieldnames(~ismember(stream_fieldnames,{'spike_arrival_times','spike_waveforms','streams_tStart'}));
    stream_fieldnames=stream_fieldnames(~cellfun(@(x) any(strfind(x,'_SR')),stream_fieldnames))';
    for FN=stream_fieldnames
        trial(t).(FN{:})=tr_in(t).(FN{:});
        trial(t).([FN{:} '_SR'])=tr_in(t).([FN{:} '_SR']);
        shift_n_samples=round(shift_in_seconds*trial(t).([FN{:} '_SR']));
        to_next_trial(t).(FN{:})=trial(t).(FN{:})(:,end-shift_n_samples:end);
        trial(t).(FN{:})(:,end-shift_n_samples:end)=[]; % cut off end of current trial (delta t = shift_in_seconds)            
        if t>1 %% append end of previous trial to the current one (delta t = shift_in_seconds)
            trial(t).(FN{:})=[to_next_trial(t-1).(FN{:}) trial(t).(FN{:})];
            trial(t).([FN{:} '_tStart'])=  tr_in(t).streams_tStart-shift_n_samples/trial(t).([FN{:} '_SR']); %% this might have to depend on sampling rate ideally...
        else
            n_samples_to_delete=round((shift_in_seconds*-1-tr_in(t).streams_tStart)*trial(t).([FN{:} '_SR']));
            if n_samples_to_delete>0
                trial(t).(FN{:})(:,1:n_samples_to_delete)=[];
                trial(t).([FN{:} '_tStart'])=tr_in(t).streams_tStart+n_samples_to_delete/trial(t).([FN{:} '_SR']);
            else
                trial(t).([FN{:} '_tStart'])=tr_in(t).streams_tStart;
            end
        end
    end
    %trial(t).streams_tStart=tr_in(t).streams_tStart;
    
    %% unspecific excel table data %% does this mean there are
    for c=1:n_chans_s,
        trial(t).channel(c).grid_x               =from_excel_per_channel(c).x{1} ;
        trial(t).channel(c).grid_y               =from_excel_per_channel(c).y{1} ;
        trial(t).channel(c).electrode_depth      =from_excel_per_channel(c).electrode_depth{1} ;
        trial(t).channel(c).target               =from_excel_per_channel(c).target{1} ;
        trial(t).channel(c).dataset              =from_excel_per_channel(c).dataset{1} ;
        trial(t).channel(c).perturbation         =from_excel_per_channel(c).perturbation{1} ;
        trial(t).channel(c).perturbation_site    =from_excel_per_channel(c).perturbation_site{1} ;
        trial(t).channel(c).site_ID              =from_excel_per_channel(c).site{1} ;        
    end
    
    %% SPIKES
    for c=1:n_chans_u,
        for u=1:n_units
            cat_wf=cat(3,tr_in.spike_waveforms);
            cat_wf_tr=vertcat(cat_wf(c,u,:));
            if any(any(~isnan(vertcat(cat_wf_tr{:}))))
                trial(t).unit(c,u).unit_ID=from_excel_per_unit(c,u).unique_neuron{1} ;
            else
                trial(t).unit(c,u).unit_ID='no unit';
            end
            trial(t).unit(c,u).site_ID              =from_excel_per_unit(c,u).site{1} ;
            trial(t).unit(c,u).SNR_rating           =from_excel_per_unit(c,u).SNR_rating{1} ;
            trial(t).unit(c,u).Single_rating        =from_excel_per_unit(c,u).Single_rating{1} ;
            trial(t).unit(c,u).stability_rating     =from_excel_per_unit(c,u).stability_rating{1} ;
            trial(t).unit(c,u).grid_x               =from_excel_per_unit(c,u).x{1} ;
            trial(t).unit(c,u).grid_y               =from_excel_per_unit(c,u).y{1} ;
            trial(t).unit(c,u).electrode_depth      =from_excel_per_unit(c,u).electrode_depth{1} ;
            trial(t).unit(c,u).target               =from_excel_per_unit(c,u).target{1} ;
            trial(t).unit(c,u).dataset              =from_excel_per_unit(c,u).dataset{1} ;
            trial(t).unit(c,u).perturbation         =from_excel_per_unit(c,u).perturbation{1} ;
            trial(t).unit(c,u).perturbation_site    =from_excel_per_unit(c,u).perturbation_site{1} ;
            
            % Waveforms (cutting off ITI ??)
            wf_idx          = tr_in(t).spike_arrival_times{c,u}>0 & tr_in(t).spike_arrival_times{c,u}<MA_out.states(t).TDT_state_onsets(end-1);
            trial(t).unit(c,u).waveforms            =tr_in(t).spike_waveforms{c,u}(wf_idx,:);
            
            %% ADDING PREVIOUS SPIKES TO CURRENT TRIAL (shift_in_seconds before trial onset)
            if t>1
                prev_trial_spike_arrival_times = tr_in(t-1).spike_arrival_times{c,u} -(trial(t).trial_onset_time -MA_out.states(t-1).trial_onset_time);
                tr_in(t).spike_arrival_times{c,u} = [prev_trial_spike_arrival_times(prev_trial_spike_arrival_times>-shift_in_seconds); tr_in(t).spike_arrival_times{c,u}];
            end
            trial(t).unit(c,u).arrival_times= tr_in(t).spike_arrival_times{c,u};
        end
    end
    
end
trial(invalid_trials)=[];
o.trial=trial;
o.block=keys.block;
end

function from_excel = get_sorted_neuron_out_xlsx_from_excel(xlsx_table, keys, unit, channel, unitexistsindata, per_unit)

idx_Date            =DAG_find_column_index(xlsx_table,'Date');
idx_Run             =DAG_find_column_index(xlsx_table,'Run');
idx_Block           =DAG_find_column_index(xlsx_table,'Block');
idx_Channel         =DAG_find_column_index(xlsx_table,'Chan');
idx_Unit            =DAG_find_column_index(xlsx_table,'Unit');
idx_Neuron_ID       =DAG_find_column_index(xlsx_table,'Neuron_ID');
idx_Site_ID         =DAG_find_column_index(xlsx_table,'Site_ID');
idx_Target          =DAG_find_column_index(xlsx_table,'Target');
idx_Hemisphere      =DAG_find_column_index(xlsx_table,'Hemisphere');
idx_Set             =DAG_find_column_index(xlsx_table,'Set');
idx_Perturbation    =DAG_find_column_index(xlsx_table,'Perturbation');
idx_Perturbation_site    =DAG_find_column_index(xlsx_table,'Perturbation_site');
idx_SNR             =DAG_find_column_index(xlsx_table,'SNR_rank');
idx_Single          =DAG_find_column_index(xlsx_table,'Single_rank');
idx_x               =DAG_find_column_index(xlsx_table,'x');
idx_y               =DAG_find_column_index(xlsx_table,'y');
idx_Electrode_depth =DAG_find_column_index(xlsx_table,'Aimed_electrode_depth');
idx_Stability       =DAG_find_column_index(xlsx_table,'Stability_rank');


r_Date  = DAG_find_row_index(xlsx_table(:,idx_Date),str2double(keys.date));
r_Block = DAG_find_row_index(xlsx_table(:,idx_Block),keys.block);
r_Run   = DAG_find_row_index(xlsx_table(:,idx_Run),keys.run);
r_Unit  = DAG_find_row_index(xlsx_table(:,idx_Unit),char(96+unit));
r_Chan  = DAG_find_row_index(xlsx_table(:,idx_Channel),channel);

unit_identifier=[keys.date '_block_' num2str(keys.block) '_run_' num2str(keys.run) '_ch_' num2str(channel) '_u_' char(96+unit)];
site_identifier=[keys.date '_block_' num2str(keys.block) '_run_' num2str(keys.run) '_ch_' num2str(channel)];
if per_unit
    row                         =find(r_Date & r_Block & r_Run & r_Chan & r_Unit );
    if numel(row)>1
        fprintf(2,'several table entries for %s, assuming 1st matching table entry \n', unit_identifier);
    end
else
    row                         =find(r_Date & r_Block & r_Run & r_Chan);
end
if isempty(row)
    from_excel.unique_neuron           = {unit_identifier};
    from_excel.site                    = {site_identifier};
    from_excel.SNR_rating              = {-1};
    from_excel.Single_rating           = {-1};
    from_excel.stability_rating        = {-1};
    from_excel.x                       = {-100};
    from_excel.y                       = {-100};
    from_excel.electrode_depth         = {-1};
    from_excel.target                  = {'unknown'};
    from_excel.dataset                 = {0};
    from_excel.perturbation            = {0};
    from_excel.perturbation_site       = {'NA'};
    if unitexistsindata
        fprintf(2, 'no matching sorting for %s \n', unit_identifier);
    end
else
    row=row(1);
    from_excel.unique_neuron           =xlsx_table(row,idx_Neuron_ID);
    from_excel.site                    =xlsx_table(row,idx_Site_ID);
    from_excel.SNR_rating              =xlsx_table(row,idx_SNR);
    from_excel.Single_rating           =xlsx_table(row,idx_Single);
    from_excel.stability_rating        =xlsx_table(row,idx_Stability);
    from_excel.x                       =xlsx_table(row,idx_x);
    from_excel.y                       =xlsx_table(row,idx_y);
    from_excel.electrode_depth         =xlsx_table(row,idx_Electrode_depth);
    from_excel.target                  =xlsx_table(row,idx_Target);
    from_excel.dataset                 =xlsx_table(row,idx_Set);
    from_excel.perturbation            =xlsx_table(row,idx_Perturbation);
    from_excel.perturbation_site       =xlsx_table(row,idx_Perturbation_site);
    % defining target including hemisphere (in case it wasnt used for target definition in the sorting table already)
    if ~isempty(idx_Hemisphere) && isempty(strfind(from_excel.target{:},'_')) && numel(xlsx_table{row,idx_Hemisphere})>0
        from_excel.target = {[from_excel.target{:} '_' upper(xlsx_table{row,idx_Hemisphere}(1))]};
    end
end
end

