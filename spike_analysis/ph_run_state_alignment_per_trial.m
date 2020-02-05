function o=ph_run_state_alignment_per_trial(MA_out,keys)
global MA_STATES

%% Defining maximum number of units
n_trials=numel(MA_out.physiology);
tr_in=[MA_out.physiology];
n_chans = size(tr_in(1).spike_arrival_times,1);
n_units=0;
for t=1:n_trials
    if ~isempty(tr_in(t).spike_arrival_times)
        n_units = max([n_units size(tr_in(t).spike_arrival_times,2)]);
    end
end

for c=1:n_chans,
    for u=1:n_units
        from_excel(c,u) = get_sorted_neuron_out_xlsx_from_excel(keys.xlsx_table, keys, u, c, 0);
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
    
    trial(t).type          =MA_out.task(t).type;
    trial(t).effector      =MA_out.task(t).effector;
    trial(t).reach_hand    =MA_out.task(t).reach_hand;
    if isnan(trial(t).reach_hand); trial(t).reach_hand=0; end;
    trial(t).choice        =MA_out.binary(t).choice;
    trial(t).success       =MA_out.binary(t).success;
    
    %% Saccade or reach positions dependend on effector
    if ismember(trial(t).effector,[0,3])
        trial(t).fix_pos        =MA_out.saccades(t).fix_pos;
        trial(t).tar_pos        =MA_out.saccades(t).tar_pos;
    else
        trial(t).fix_pos        =MA_out.reaches(t).fix_pos;
        trial(t).tar_pos        =MA_out.reaches(t).tar_pos;
    end
    
    %% adding previous trial
    trial(t).trial_onset_time        =MA_out.states(t).trial_onset_time;
    if t>1 % && keys.add_previous_trial_spikes
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
    
    get_expected_MA_states(trial(t).type,trial(t).effector,keys.effectors_on_same_figure);    %% set to 1 to allow later processing(?)
    keys.ALL_EPOCHS=keys.EPOCHS_PER_TYPE{trial(t).type};
    keys.EPOCHS=keys.ALL_EPOCHS(ismember([keys.ALL_EPOCHS{:,2}],MA_STATES.all_states),:);
    states_to_process_trial=[keys.EPOCHS{:,2}];    
    
    sac_ini = MA_out.states(t).start_obs + MA_out.saccades(t).lat;
    sac_end = sac_ini + MA_out.saccades(t).dur;
    rea_ini = MA_out.states(t).start_obs + MA_out.reaches(t).lat;
    rea_end = rea_ini + MA_out.reaches(t).dur;
    
    trial_states=[trial_states(1:end-1) MA_STATES.SAC_INI MA_STATES.SAC_END MA_STATES.REA_INI MA_STATES.REA_END 90]; %MA_STATES.ITI ?
    trial_states_onset=[trial_states_onset(1:end-1) sac_ini sac_end rea_ini rea_end trial_states_onset(end) NaN];
    trial(t).states_onset=trial_states_onset(1:end-1);
    
    for c=1:n_chans,
        for u=1:n_units
            cat_wf=cat(3,MA_out.physiology.spike_waveforms);
            cat_wf_tr=vertcat(cat_wf(c,u,:));
            if any(any(~isnan(vertcat(cat_wf_tr{:}))))
                trial(t).unit(c,u).unit_ID=from_excel(c,u).unique_neuron{1} ;
            else
                trial(t).unit(c,u).unit_ID='no unit';
            end
            trial(t).unit(c,u).SNR_rating           =from_excel(c,u).SNR_rating{1} ;
            trial(t).unit(c,u).Single_rating        =from_excel(c,u).Single_rating{1} ;
            trial(t).unit(c,u).stability_rating     =from_excel(c,u).stability_rating{1} ;
            trial(t).unit(c,u).grid_x               =from_excel(c,u).x{1} ;
            trial(t).unit(c,u).grid_y               =from_excel(c,u).y{1} ;
            trial(t).unit(c,u).electrode_depth      =from_excel(c,u).electrode_depth{1} ;
            trial(t).unit(c,u).target               =from_excel(c,u).target{1} ;
            
            % Waveforms (cutting off ITI)
            wf_idx=MA_out.physiology(t).spike_arrival_times{c,u}>0 & MA_out.physiology(t).spike_arrival_times{c,u}<MA_out.states(t).TDT_state_onsets(end-1);
            trial(t).unit(c,u).waveforms            =MA_out.physiology(t).spike_waveforms{c,u}(wf_idx,:);
            
            %% ADDING PREVIOUS SPIKES TO CURRENT TRIAL
            if t>1 % && keys.add_previous_trial_spikes
                prev_trial_spike_arrival_times = MA_out.physiology(t-1).spike_arrival_times{c,u} - trial(t-1).states_onset(end) + trial(t).states_onset(1);
                MA_out.physiology(t).spike_arrival_times{c,u} = [prev_trial_spike_arrival_times; MA_out.physiology(t).spike_arrival_times{c,u}];
                PSTH_x_axis=trial(t-1).states_onset(1)- trial(t-1).states_onset(end) + trial(t).states_onset(1):keys.PSTH_binwidth:trial(t).states_onset(end)+0.2;%% 200 additional ms
            else
                PSTH_x_axis=-2:keys.PSTH_binwidth:trial(t).states_onset(end);
            end
            trial(t).unit(c,u).trial_spike_density=double(proper_PSTH(MA_out.physiology(t).spike_arrival_times{c,u},PSTH_x_axis)/keys.PSTH_binwidth);
            [~,onset_indexes_to_process]=ismember(states_to_process_trial,trial_states);
            onset_indexes_to_process(onset_indexes_to_process==0)=numel(trial_states_onset); %% so they can be assgined to NAns later
            
            
            for s=1:numel(states_to_process_trial);
                current_state=states_to_process_trial(s);
                trial(t).unit(c,u).spikes_per_state(s).state=keys.EPOCHS{s,1};
                if ismember(current_state,trial_states) %% why not defining at least states?
                    state_idx=trial_states==current_state;
                    trial(t).unit(c,u).spikes_per_state(s).state_onset=trial_states_onset(state_idx);
                    trial(t).unit(c,u).spikes_per_state(s).state_onsets=trial_states_onset(onset_indexes_to_process)-trial_states_onset(state_idx);
                    
                    %% RASTER and PSTH
                    arrival_times       =MA_out.physiology(t).spike_arrival_times{c,u}-trial_states_onset(state_idx);
                    t_before_state      =keys.EPOCHS{s,3};
                    t_after_state       =keys.EPOCHS{s,4};
%                     arrival_times_PSTH  =arrival_times(arrival_times>t_before_state-keys.PSTH_binwidth & arrival_times<t_after_state+keys.PSTH_binwidth);
                    arrival_times_PSTH  =arrival_times(arrival_times>t_before_state & arrival_times<t_after_state);
                    trial(t).unit(c,u).spikes_per_state(s).arrival_times=arrival_times_PSTH;
                    if mod((PSTH_x_axis(1)-trial_states_onset(state_idx))*(1/keys.PSTH_binwidth)-0.5,1)<0.0001 %%LInus 20151106 block 3
                        trial_states_onset(state_idx)=trial_states_onset(state_idx)-0.0003;
                    elseif mod((PSTH_x_axis(1)-trial_states_onset(state_idx))*(1/keys.PSTH_binwidth)-0.5,1)>0.9999 %%LInus 20151106 block 3
                        trial_states_onset(state_idx)=trial_states_onset(state_idx)+0.0003;                        
                    end
                    PSTH_x_shifted=round(double(PSTH_x_axis-trial_states_onset(state_idx))/keys.PSTH_binwidth)*keys.PSTH_binwidth;
                    PSTH_x_axis_indexes=PSTH_x_shifted-t_before_state>-0.0001 & PSTH_x_shifted-t_after_state<0.0001;
                    %PSTH_x_axis_indexes=PSTH_x_shifted-t_before_state>-0 & PSTH_x_shifted-t_after_state<0; %important here to reproduce microstim ephys stuff
                    trial(t).unit(c,u).spikes_per_state(s).spike_density=trial(t).unit(c,u).trial_spike_density(PSTH_x_axis_indexes);
                                       
                    %% eye and hand
                    state_time=trial(t).time_axis-trial_states_onset(state_idx);
                    time_indexes=state_time>=t_before_state & state_time<=t_after_state;
                    trial(t).unit(c,u).spikes_per_state(s).x_eye=trial(t).x_eye(time_indexes);
                    trial(t).unit(c,u).spikes_per_state(s).y_eye=trial(t).y_eye(time_indexes);
                    trial(t).unit(c,u).spikes_per_state(s).x_hnd=trial(t).x_hnd(time_indexes);
                    trial(t).unit(c,u).spikes_per_state(s).y_hnd=trial(t).y_hnd(time_indexes);
                    trial(t).unit(c,u).spikes_per_state(s).time_axis=state_time(time_indexes);
                    
                    %% FR
                    t_before_state=keys.EPOCHS{s,5};
                    t_after_state=keys.EPOCHS{s,6};
%                    PSTH_x_axis_indexes=PSTH_x_shifted-t_before_state>-0.0001 & PSTH_x_shifted-t_after_state<0.0001;
                    trial(t).unit(c,u).spikes_per_state(s).FR=sum(arrival_times>=t_before_state & arrival_times<=t_after_state)/(t_after_state-t_before_state);
%                     PSTH_x_axis_indexes=PSTH_x_shifted>=single(t_before_state) & PSTH_x_shifted<=single(t_after_state); %important here to reproduce microstim ephys stuff
%                      trial(t).unit(c,u).spikes_per_state(s).FR=sum(trial(t).unit(c,u).trial_spike_density(PSTH_x_axis_indexes))/sum(PSTH_x_axis_indexes);
% %                     
                 
                else %create dummy
                    trial(t).unit(c,u).spikes_per_state(s)=struct('state',NaN,'state_onset',NaN,'state_onsets',NaN,'arrival_times',NaN,'spike_density',NaN,'x_eye',NaN,'y_eye',NaN,'x_hnd',NaN,'y_hnd',NaN,'time_axis',NaN,'FR',NaN);
                end
            end
        end
    end
    
end
trial(invalid_trials)=[];
o.trial=trial;
o.block=keys.block;
end

function from_excel = get_sorted_neuron_out_xlsx_from_excel(xlsx_table, keys, unit, channel, verbose)

idx_Date=find_column_index(xlsx_table,'Date');
idx_Run=find_column_index(xlsx_table,'Run');
idx_Block=find_column_index(xlsx_table,'Block');
idx_Channel=find_column_index(xlsx_table,'Chan');
idx_Unit=find_column_index(xlsx_table,'Unit');
idx_Neuron_ID=find_column_index(xlsx_table,'Neuron_ID');
idx_Target=find_column_index(xlsx_table,'Target');
idx_SNR=find_column_index(xlsx_table,'SNR rank');
idx_Single=find_column_index(xlsx_table,'Single rank');
idx_x=find_column_index(xlsx_table,'x');
idx_y=find_column_index(xlsx_table,'y');
idx_Electrode_depth=find_column_index(xlsx_table,'Aimed electrode_depth');
idx_Stability=find_column_index(xlsx_table,'Stability rank');

r_Date  = find_row_index_working(xlsx_table(:,idx_Date),str2double(keys.date));
r_Block = find_row_index_working(xlsx_table(:,idx_Block),keys.block);
r_Run   = find_row_index_working(xlsx_table(:,idx_Run),keys.run);
r_Unit  = find_row_index_working(xlsx_table(:,idx_Unit),keys.uname(unit));
r_Chan  = find_row_index_working(xlsx_table(:,idx_Channel),channel);

neuron_idx                         =find(r_Date & r_Block & r_Run & r_Unit & r_Chan);
if isempty(neuron_idx)
    from_excel.unique_neuron           = {[keys.date '_block_' num2str(keys.block) '_run_' num2str(keys.run) '_ch_' num2str(channel) '_u_' keys.uname{unit}]};
    from_excel.SNR_rating              = {-1};
    from_excel.Single_rating           = {-1};
    from_excel.x                       = {-100};
    from_excel.y                       = {-100};
    from_excel.electrode_depth         = {-1};
    from_excel.stability_rating        = {-1};
    from_excel.target                  = {'unknown'};
    if verbose
        disp(['[' 8 'NO MATCHING SORTING EXISTS!!!' ']' 8])
    end
else
    from_excel.unique_neuron           =xlsx_table(neuron_idx,idx_Neuron_ID);
    from_excel.SNR_rating              =xlsx_table(neuron_idx,idx_SNR);
    from_excel.Single_rating           =xlsx_table(neuron_idx,idx_Single);
    from_excel.x                       =xlsx_table(neuron_idx,idx_x);
    from_excel.y                       =xlsx_table(neuron_idx,idx_y);
    from_excel.electrode_depth         =xlsx_table(neuron_idx,idx_Electrode_depth);
    from_excel.stability_rating        =xlsx_table(neuron_idx,idx_Stability);
    from_excel.target                  =xlsx_table(neuron_idx,idx_Target);
end
end

function PSTH=proper_PSTH(spike_arrival_times,t,s)
bin=t(2)-t(1);
if nargin<3
    s=2*bin; %%2*bin
end
PSTH=zeros(size(t));
spike_arrival_times=spike_arrival_times(spike_arrival_times>=t(1) & spike_arrival_times<=t(end));
for spk=1:numel(spike_arrival_times)
    gaussian=exp(-(t-spike_arrival_times(spk)).^2/(2*s^2));
    gaussian=gaussian/sum(gaussian);
    PSTH=PSTH+gaussian;
end
end
