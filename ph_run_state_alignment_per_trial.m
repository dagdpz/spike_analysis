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
n_units=0;
%% unitsperchannelmatrix 
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

tix=get_indexes(keys.sorting_table_sites);
for c=1:n_chans_s
    from_excel_per_channel(c) = get_recording_info(keys.sorting_table_sites, tix,keys, 1, c, 0,0);
end
tix=get_indexes(keys.sorting_table_units);
cat_wf=cat(3,tr_in.spike_waveforms);
for c=1:n_chans_u
    for u=1:n_units
        from_excel_per_unit(c,u) = get_recording_info(keys.sorting_table_units, tix,keys, u, c, unitsperchannelmatrix(c,u),1);
        cat_wf_tr=vertcat(cat_wf(c,u,:));
        if ~any(any(~isnan(vertcat(cat_wf_tr{:}))))
            from_excel_per_unit(c,u).unit_ID='no unit';
        end
    end
end


%% Preallocation !?
NaNpar=num2cell(NaN(1,n_trials));
trial=struct('type',NaNpar,'effector',NaNpar,'reach_hand',NaNpar,'choice',NaNpar,'success',NaNpar,'fix_pos',NaNpar,'tar_pos',NaNpar,'trial_onset_time',NaNpar,...
     'time_axis',NaNpar,'x_eye',NaNpar,'y_eye',NaNpar,'x_hnd',NaNpar,'y_hnd',NaNpar,'states_onset',NaNpar); %add all fields of trial(?)
trials_wo_phys=[];
trials_wo_cond=[];

[trial.unit]=deal(from_excel_per_unit);
[trial.channel]=deal(from_excel_per_channel);
    
%% LFPs & other streams??
shift_in_seconds=1;
stream_fieldnames=fieldnames(tr_in);
stream_fieldnames=stream_fieldnames(~ismember(stream_fieldnames,{'spike_arrival_times','spike_waveforms','streams_tStart'}));
stream_fieldnames=stream_fieldnames(~cellfun(@(x) any(strfind(x,'_SR')) || any(strfind(x,'_t0_from_rec_start') ),stream_fieldnames))';

for FN=stream_fieldnames
    sr={tr_in.([FN{:} '_SR'])};
    sr(cellfun(@isempty,sr))={0};
    sr=cell2mat(sr);
    shift_n_samples=round(shift_in_seconds*sr)';
    % adding last second of previous trial to the beginning of the next trial
    tempstruct=[{tr_in(1).(FN{:})}; arrayfun(@(x,y,z) [x.(FN{:})(:,end-z+1:end) y.(FN{:})(:,1:end-z)],tr_in(1:end-1),tr_in(2:end),shift_n_samples(1:end-1),'UniformOutput',false)];
    % shorten first trial (remove stuff way before task)
    n_samples_to_delete=round((shift_in_seconds*-1-tr_in(1).streams_tStart)*sr(1));
    if n_samples_to_delete>-1
        trial(1).(FN{:})(:,1:n_samples_to_delete)=[];
        tstart={tr_in(1).streams_tStart+n_samples_to_delete/sr(1)};
    else
        tstart={tr_in(1).streams_tStart};
    end
    tstart=[tstart; arrayfun(@(x,y) x.streams_tStart-y,tr_in(2:end),shift_n_samples(2:end)./sr(2:end)','UniformOutput',false)];    
    [trial.([FN{:} '_tStart'])]=  deal(tstart{:});
    [trial.(FN{:})]=deal(tempstruct{:});
    [trial.([FN{:} '_SR'])]=tr_in.([FN{:} '_SR']);
    [trial.([FN{:} '_t0_from_rec_start'])]=tr_in.([FN{:} '_t0_from_rec_start']);
end
%% eye and hand traces
tempstruct=[{MA_out.raw(1).time_axis}; arrayfun(@(x,y,z,a) [x.time_axis-(y.trial_onset_time-z.trial_onset_time) a.time_axis],...
             MA_out.raw(1:end-1),MA_out.states(2:end),MA_out.states(1:end-1),MA_out.raw(2:end),'UniformOutput',false)];
         
trace_idx=cellfun(@(x) 1:MA_out.keys.calc_val.i_sample_rate*keys.PSTH_binwidth:numel(x),tempstruct,'UniformOutput',false);
tempstruct=cellfun(@(x,y) x(y),tempstruct,trace_idx,'UniformOutput',false);      
[trial.time_axis] = deal(tempstruct{:});

trace_fieldnames={'x_eye','y_eye','x_hnd','y_hnd'};
for FN=trace_fieldnames
    fn=FN{:};
    tempstruct=[{MA_out.raw(1).(fn)}; arrayfun(@(x,y,z) [x.(fn) y.(fn)],MA_out.raw(1:end-1),MA_out.raw(2:end),'UniformOutput',false)];    
    %% resampling    
    tempstruct=cellfun(@(x,y) x(y),tempstruct,trace_idx,'UniformOutput',false);      
    [trial.(fn)] = deal(tempstruct{:});
end
%% spikes

for t=1:numel(tr_in)
    
    t1=min(MA_out.states(t).TDT_state_onsets);
    t2=max(MA_out.states(t).TDT_state_onsets(1:end-1));
    if ~isempty(tr_in(t).spike_waveforms)
    unit_wf=cell2struct(tr_in(t).spike_waveforms,'waveforms',3);
    [trial(t).unit.waveforms]=unit_wf.waveforms;
    if t>1
        %% add previous trial's spikes to the beginning
        prev_t_correction=-(MA_out.states(t).trial_onset_time -MA_out.states(t-1).trial_onset_time);
        AA=cellfun(@(x,y) [x(x+prev_t_correction>-shift_in_seconds)+prev_t_correction;y],tr_in(t-1).spike_arrival_times,tr_in(t).spike_arrival_times,'Uniformoutput',false);
        unit_at=cell2struct(AA,'arrival_times',3);
    else
        unit_at=cell2struct(tr_in(t).spike_arrival_times,'arrival_times',3);
    end
    [trial(t).unit.arrival_times]=unit_at.arrival_times;
    
    AA=arrayfun(@(x) sum(x.arrival_times>t1 & x.arrival_times<t2)/(t2-t1),trial(t).unit,'Uniformoutput',false);
    [trial(t).unit.FR_average]=AA{:};
    else
        
    end
end

%% main loop
for t=1:n_trials
    trial_states=MA_out.states(t).TDT_states;
    trial_states_onset=MA_out.states(t).TDT_state_onsets;
    %% what to do with last trial? is the last state 1?
    %% unit Lin_20150521_02 debug
    %%W:\Data\Linus_phys_combined_monkeypsych_TDT\20150521\Lincombined2015-05-21_03_block_01.mat
    if isempty(trial_states_onset) ||...
            (t==n_trials && trial_states(end)~=1) ||...
            (t==n_trials-1 && trial_states(end)~=1) ||...
            (trial_states(end)~=90 && trial_states(end)~=99 && trial_states(end)~=1)%% Curius 20150603 block 4 trial 166
        %&& (t==n_trials-1 || t==n_trials) % last trial bug that should get fixed in TDT_trial_struct_working  line 30
        % if trial_states(end) is 99, this can be fixed here! no need to  discard those trials
        trials_wo_phys=[trials_wo_phys t];
        continue;
    end
    trial(t).date               = str2num(keys.date);
    trial(t).block              = keys.block;
    trial(t).run                = keys.run;%MA_out.selected(t).run; % cause MA run is wrong, it actually corresponds to block!
    trial(t).n                  = MA_out.selected(t).trials;
    trial(t).type               = MA_out.task(t).type;
    trial(t).effector           = MA_out.task(t).effector;
    trial(t).reach_hand         = MA_out.task(t).reach_hand;
    trial(t).correct_targets    = MA_out.task(t).correct_targets;
    
    trial(t).n_nondistractors          = MA_out.task(t).n_nondistractors;
    trial(t).n_distractors             = MA_out.task(t).n_distractors;
    trial(t).difficulty                = MA_out.task(t).difficulty;
    trial(t).stimuli_in_2hemifields    = MA_out.task(t).stimuli_in_2hemifields;
    
    if trial(t).n_nondistractors == 0 && trial(t).n_distractors == 2 ||  trial(t).n_nondistractors == 1 && trial(t).n_distractors == 1
        trial(t).stimulustype = 1; %% single stimulus
    elseif trial(t).n_nondistractors == 2 || trial(t).n_distractors == 3
        trial(t).stimulustype = 2; %% TT / DD
    elseif trial(t).n_nondistractors == 1 && trial(t).n_distractors == 2
        trial(t).stimulustype = 3; %% target distractor
    else
        trial(t).stimulustype = 1;
    end
    
    if isnan(trial(t).reach_hand); trial(t).reach_hand=0; end;
    trial(t).choice        = MA_out.binary(t).choice;
    trial(t).success       = MA_out.binary(t).success;
    trial(t).completed     = MA_out.binary(t).completed;
    
    
    %% Saccade or reach positions dependend on effector
    MA_out.reaches(t).ini_all       = MA_out.reaches(t).ini;
    MA_out.reaches(t).end_all       = MA_out.reaches(t).ini+MA_out.reaches(t).dur;
    MA_out.reaches(t).endpos_all    = MA_out.reaches(t).endpos;
    MA_out.reaches(t).startpos_all  = MA_out.reaches(t).startpos;
    trial(t).sac_off                = MA_out.saccades(t).endpos-MA_out.saccades(t).tar_pos;
    trial(t).sac_lat                = MA_out.saccades(t).lat;
    trial(t).rea_off                = MA_out.reaches(t).endpos-MA_out.reaches(t).tar_pos;
    trial(t).rea_lat                = MA_out.reaches(t).lat;
    
    sac_ini = MA_out.states(t).start_obs + MA_out.saccades(t).lat;
    sac_end = sac_ini + MA_out.saccades(t).dur;
    rea_ini = MA_out.states(t).start_obs + MA_out.reaches(t).lat;
    rea_end = rea_ini + MA_out.reaches(t).dur;
    tri_end = MA_out.states(t).start_end;
    
    if ismember(trial(t).effector,[0,3])
        Movement=MA_out.saccades(t);
        mov_ini=sac_ini;
        mov_end=sac_end;
    else
        Movement=MA_out.reaches(t);
        mov_ini=rea_ini;
        mov_end=rea_end;
    end
    
    trial(t).fix_pos     = Movement.fix_pos;
    trial(t).tar_pos     = Movement.tar_pos;
    trial(t).cue_pos     = Movement.cue_pos;
    trial(t).cue_shape   = Movement.all_convexities(1);
    
    % for distractor task
    trial(t).all_tar_pos        = Movement.all_tar_pos;
    trial(t).col_dim            = Movement.col_dim;
    trial(t).col_bri            = Movement.col_bri;
    trial(t).target_selected    = Movement.target_selected;
    
    if  numel(trial(t).all_tar_pos) == 1 %% KK NOTLÖSUNG
        correct_pos = trial(t).all_tar_pos;
    else
        correct_pos=trial(t).all_tar_pos(trial(t).correct_targets);
    end
    if correct_pos == trial(t).fix_pos(1) %% distractor only or double_distractor
        correct_pos=trial(t).all_tar_pos(trial(t).all_tar_pos~=0);
    end
    trial(t).stm_pos=correct_pos(1);
    
    MPA_get_expected_states(trial(t).type,trial(t).effector,0);    %% set to 1 to allow later processing
    movement_states       = [MA_STATES.SAC_INI MA_STATES.SAC_END MA_STATES.REA_INI MA_STATES.REA_END MA_STATES.MOV_INI MA_STATES.MOV_END MA_STATES.TRI_END];
    movement_onsets       = [sac_ini sac_end rea_ini rea_end mov_ini mov_end tri_end];
    movement_onsets       = movement_onsets(ismember(movement_states,MA_STATES.all_states));
    movement_states       = movement_states(ismember(movement_states,MA_STATES.all_states));
    
    trial_states            =[trial_states(1:end-1) movement_states MA_STATES.ITI_END];
    trial_states_onset      =[trial_states_onset(1:end-1) movement_onsets trial_states_onset(end)];
    [~,tr_state_idx]        =unique(trial_states,'last');
    trial(t).states         =trial_states(sort(tr_state_idx));
    trial(t).states_onset   =trial_states_onset(sort(tr_state_idx));
    
    %% adding previous trial
    trial(t).run_onset_time          =MA_out.states(t).run_onset_time;
    trial(t).trial_onset_time        =MA_out.states(t).trial_onset_time;
    
%     if t>1 % && keys.add_previous_trial_spikes     shift_in_seconds=1;
%         trial(t).time_axis = [MA_out.raw(t-1).time_axis-(trial(t).trial_onset_time -MA_out.states(t-1).trial_onset_time) MA_out.raw(t).time_axis];
%         trial(t).x_eye = [MA_out.raw(t-1).x_eye MA_out.raw(t).x_eye];
%         trial(t).y_eye = [MA_out.raw(t-1).y_eye MA_out.raw(t).y_eye];
%         trial(t).x_hnd = [MA_out.raw(t-1).x_hnd MA_out.raw(t).x_hnd];
%         trial(t).y_hnd = [MA_out.raw(t-1).y_hnd MA_out.raw(t).y_hnd];
%     else
%         trial(t).time_axis = MA_out.raw(t).time_axis;
%         trial(t).x_eye = MA_out.raw(t).x_eye;
%         trial(t).y_eye = MA_out.raw(t).y_eye;
%         trial(t).x_hnd = MA_out.raw(t).x_hnd;
%         trial(t).y_hnd = MA_out.raw(t).y_hnd;
%     end
%     %% resample traces for speed !!!
%     trace_idx=1:MA_out.keys.calc_val.i_sample_rate*keys.PSTH_binwidth:numel(trial(t).time_axis);
%     trial(t).time_axis  = trial(t).time_axis(trace_idx);
%     trial(t).x_eye      = trial(t).x_eye(trace_idx);
%     trial(t).y_eye      = trial(t).y_eye(trace_idx);
%     trial(t).x_hnd      = trial(t).x_hnd(trace_idx);
%     trial(t).y_hnd      = trial(t).y_hnd(trace_idx);
    
    
%     %% LFPs & other streams??
%     shift_in_seconds=1;
%     stream_fieldnames=fieldnames(tr_in);
%     stream_fieldnames=stream_fieldnames(~ismember(stream_fieldnames,{'spike_arrival_times','spike_waveforms','streams_tStart'}));
%     stream_fieldnames=stream_fieldnames(~cellfun(@(x) any(strfind(x,'_SR')) || any(strfind(x,'_t0_from_rec_start') ),stream_fieldnames))';
%     for FN=stream_fieldnames
%         trial(t).(FN{:})=tr_in(t).(FN{:});
%         trial(t).([FN{:} '_SR'])=tr_in(t).([FN{:} '_SR']);
%         trial(t).([FN{:} '_t0_from_rec_start'])=tr_in(t).([FN{:} '_t0_from_rec_start']);
%         shift_n_samples=round(shift_in_seconds*trial(t).([FN{:} '_SR']));
%         to_next_trial(t).(FN{:})=trial(t).(FN{:})(:,end-shift_n_samples:end);
%         trial(t).(FN{:})(:,end-shift_n_samples:end)=[]; % cut off end of current trial (delta t = shift_in_seconds)
%         if t>1 %% append end of previous trial to the current one (delta t = shift_in_seconds)
%             trial(t).(FN{:})=[to_next_trial(t-1).(FN{:}) trial(t).(FN{:})];
%             trial(t).([FN{:} '_tStart'])=  tr_in(t).streams_tStart-shift_n_samples/trial(t).([FN{:} '_SR']); %% this might have to depend on sampling rate ideally...
%         else
%             n_samples_to_delete=round((shift_in_seconds*-1-tr_in(t).streams_tStart)*trial(t).([FN{:} '_SR']));
%             if n_samples_to_delete>0
%                 trial(t).(FN{:})(:,1:n_samples_to_delete)=[];
%                 trial(t).([FN{:} '_tStart'])=tr_in(t).streams_tStart+n_samples_to_delete/trial(t).([FN{:} '_SR']);
%             else
%                 trial(t).([FN{:} '_tStart'])=tr_in(t).streams_tStart;
%             end
%         end
%     end
%     %trial(t).streams_tStart=tr_in(t).streams_tStart;
%     
%     %% unspecific excel table data
%     for c=1:n_chans_s,
%         trial(t).channel(c).grid_x               =from_excel_per_channel(c).x{1} ;
%         trial(t).channel(c).grid_y               =from_excel_per_channel(c).y{1} ;
%         trial(t).channel(c).electrode_depth      =from_excel_per_channel(c).electrode_depth{1} ;
%         trial(t).channel(c).target               =from_excel_per_channel(c).target{1} ;
%         trial(t).channel(c).dataset              =from_excel_per_channel(c).dataset{1} ;
%         trial(t).channel(c).perturbation         =from_excel_per_channel(c).perturbation{1} ;
%         trial(t).channel(c).perturbation_site    =from_excel_per_channel(c).perturbation_site{1} ;
%         trial(t).channel(c).site_ID              =from_excel_per_channel(c).site{1} ;
%     end
%     
%     %% SPIKES
%     for c=1:n_chans_u,
%         for u=1:n_units
%             cat_wf=cat(3,tr_in.spike_waveforms);
%             cat_wf_tr=vertcat(cat_wf(c,u,:));
%             if any(any(~isnan(vertcat(cat_wf_tr{:}))))
%                 trial(t).unit(c,u).unit_ID=from_excel_per_unit(c,u).unit_ID ;
%             else
%                 trial(t).unit(c,u).unit_ID='no unit';
%             end
%             trial(t).unit(c,u).site_ID              =from_excel_per_unit(c,u).site{1} ;
%             trial(t).unit(c,u).SNR_rating           =from_excel_per_unit(c,u).SNR_rating{1} ;
%             trial(t).unit(c,u).Single_rating        =from_excel_per_unit(c,u).Single_rating{1} ;
%             trial(t).unit(c,u).stability_rating     =from_excel_per_unit(c,u).stability_rating{1} ;
%             trial(t).unit(c,u).grid_x               =from_excel_per_unit(c,u).x{1} ;
%             trial(t).unit(c,u).grid_y               =from_excel_per_unit(c,u).y{1} ;
%             trial(t).unit(c,u).electrode_depth      =from_excel_per_unit(c,u).electrode_depth{1} ;
%             trial(t).unit(c,u).target               =from_excel_per_unit(c,u).target{1} ;
%             trial(t).unit(c,u).dataset              =from_excel_per_unit(c,u).dataset{1} ;
%             trial(t).unit(c,u).perturbation         =from_excel_per_unit(c,u).perturbation{1} ;
%             trial(t).unit(c,u).perturbation_site    =from_excel_per_unit(c,u).perturbation_site{1} ;
%             
%             % Waveforms (cutting off ITI ??)
%             t1=min(trial(t).states_onset);
%             t2=max(trial(t).states_onset(1:end-1));
%             %wf_idx                                = tr_in(t).spike_arrival_times{c,u}>t1 & tr_in(t).spike_arrival_times{c,u}<t2;
%             trial(t).unit(c,u).waveforms          = tr_in(t).spike_waveforms{c,u};%(wf_idx,:);
%             
%             %% ADDING PREVIOUS SPIKES TO CURRENT TRIAL (shift_in_seconds before trial onset)
%             if t>1 && ~ismember(t-1,trials_wo_phys)
%                 prev_trial_spike_arrival_times = tr_in(t-1).spike_arrival_times{c,u} -(trial(t).trial_onset_time -MA_out.states(t-1).trial_onset_time);
%                 tr_in(t).spike_arrival_times{c,u} = [prev_trial_spike_arrival_times(prev_trial_spike_arrival_times>-shift_in_seconds); tr_in(t).spike_arrival_times{c,u}];
%             end
%             AT= tr_in(t).spike_arrival_times{c,u};
%             trial(t).unit(c,u).arrival_times=AT;
%             trial(t).unit(c,u).FR_average=sum(AT>t1 & AT<t2)/(t2-t1);
%         end
%     end
end

invalid_trials=sort([trials_wo_phys trials_wo_cond]); % differentiation between phys not present and condition mismatches (which are NOT excluded at this stage any more !)
trial(invalid_trials)=[];

%% automatic stability (dependent on fano factor of Frs per trial
if ~isempty(trial) && (keys.cal.automatic_stablity || keys.cal.automatic_SNR)
    units_cat=cat(3,trial.unit);
    for c=1:n_chans_u,
        for u=1:n_units
            
            % stability
            FRs_cat=[units_cat(c,u,:).FR_average];
            %FRs_cat=FRs_cat(FRs_cat~=0); %% .... hmmmmmmmhmmmmm
            %           unit_mean=nanmean(FRs_cat);
            %           unit_std=nanstd(FRs_cat);
            %           confidence_interval=3*unit_std;
            cutoff=FRs_cat<0.5; % | FRs_cat<unit_mean-confidence_interval;
            cut_diff=diff([true cutoff]);
            first_valid=1;
            last_valid=numel(FRs_cat);
            if cutoff(1)
                first_valid=find(cut_diff,1,'first');
            end
            if cutoff(end)
                last_valid=find(cut_diff,1,'last')-1;
            end
            if all(cutoff)
                last_valid=1;
                first_valid=numel(FRs_cat);
            end
            
            %FRs_cat=FRs_cat(FRs_cat>(log(1+exp(unit_mean-confidence_interval))) & FRs_cat<(log(1+exp(unit_mean+confidence_interval))));
            FRs_cat=FRs_cat(first_valid:last_valid);
            FRs_cat=smooth(FRs_cat,10);
            stability=nanmean(FRs_cat)/nanstd(FRs_cat);
            
            % SNR
            WFs_cat=vertcat(units_cat(c,u,first_valid:last_valid).waveforms);
            waveform_average=mean(WFs_cat,1);
            waveform_std=std(WFs_cat,0,1);
            waveform_amplitude=max(waveform_average)-min(waveform_average);
            snr=waveform_amplitude/mean(waveform_std); % redefine "noise" based on broadband (?)
            
            % single-unit'ness as it was defined by Kim et al., 2009,
            % J.Neuro:
            % "Units with <1% of ISIs <3 ms were classified as single 
            % units. All others were classified as multiunits."
            ISI_cat = arrayfun(@(x) diff(x.arrival_times)', units_cat(c,u,:), 'UniformOutput', false);
            ISI_cat = [ISI_cat{:}];
            singleunitness = sum(ISI_cat < 0.003)/length(ISI_cat); % fraction of ISIs shorter than 3ms
            
            for t=1:numel(trial) % another loop? not so cool
                if keys.cal.automatic_stablity
                    if t>=first_valid && t<=last_valid
                        trial(t).unit(c,u).stability_rating=single(stability);
                    else
                        trial(t).unit(c,u).stability_rating=single(NaN);
                    end
                end
                if keys.cal.automatic_SNR
                    if ~isempty(snr)
                        trial(t).unit(c,u).SNR_rating=single(snr);
                    else
                        trial(t).unit(c,u).SNR_rating=single(NaN);
                    end
                end
                if keys.cal.automatic_singleunitness
                    if ~isempty(singleunitness)
                        trial(t).unit(c,u).Single_rating=single(singleunitness);
                    else
                        trial(t).unit(c,u).Single_rating=single(NaN);
                    end
                end
            end
        end
    end
end

o.trial=trial;
o.block=keys.block;
end

function tix=get_indexes(xlsx_table)

tix.Date                =DAG_find_column_index(xlsx_table,'Date');
tix.Run                 =DAG_find_column_index(xlsx_table,'Run');
tix.Block               =DAG_find_column_index(xlsx_table,'Block');
tix.Channel             =DAG_find_column_index(xlsx_table,'Chan');
tix.Unit                =DAG_find_column_index(xlsx_table,'Unit');
tix.Neuron_ID           =DAG_find_column_index(xlsx_table,'Neuron_ID');
tix.Site_ID             =DAG_find_column_index(xlsx_table,'Site_ID');
tix.Target              =DAG_find_column_index(xlsx_table,'Target');
tix.Hemisphere          =DAG_find_column_index(xlsx_table,'Hemisphere');
tix.Set                 =DAG_find_column_index(xlsx_table,'Set');
tix.Perturbation        =DAG_find_column_index(xlsx_table,'Perturbation');
tix.Perturbation_site   =DAG_find_column_index(xlsx_table,'Perturbation_site');
tix.SNR                 =DAG_find_column_index(xlsx_table,'SNR_rank');
tix.Single              =DAG_find_column_index(xlsx_table,'Single_rank');
tix.x                   =DAG_find_column_index(xlsx_table,'x');
tix.y                   =DAG_find_column_index(xlsx_table,'y');
tix.Electrode_depth     =DAG_find_column_index(xlsx_table,'Aimed_electrode_depth');
tix.Stability           =DAG_find_column_index(xlsx_table,'Stability_rank');
end

function from_excel = get_recording_info(xlsx_table, tix, keys, unit, channel, unitexistsindata, per_unit)

r_Date  = [0 xlsx_table{2:end,tix.Date}]'   ==str2double(keys.date);
r_Block = [0 xlsx_table{2:end,tix.Block}]'  ==keys.block;
r_Run   = [0 xlsx_table{2:end,tix.Run}]'    ==keys.run;
r_Chan  = [0 xlsx_table{2:end,tix.Channel}]'==channel;
r_Unit  = ismember(xlsx_table(:,tix.Unit),char(96+unit));

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
    from_excel.unit_ID                 = unit_identifier;
    from_excel.site_ID                 = site_identifier;
    from_excel.SNR_rating              = -1;
    from_excel.Single_rating           = -1;
    from_excel.stability_rating        = -1;
    from_excel.grid_x                  = -100;
    from_excel.grid_y                  = -100;
    from_excel.electrode_depth         = -1;
    from_excel.target                  = 'unknown';
    from_excel.dataset                 = 0;
    from_excel.perturbation            = 0;
    from_excel.perturbation_site       = 'NA';
    if unitexistsindata
        fprintf(2, 'no matching sorting for %s \n', unit_identifier);
    end
else
    row=row(1);
    from_excel.unit_ID                 =xlsx_table{row,tix.Neuron_ID};
    from_excel.site_ID                 =xlsx_table{row,tix.Site_ID};
    from_excel.SNR_rating              =xlsx_table{row,tix.SNR};
    from_excel.Single_rating           =xlsx_table{row,tix.Single};
    from_excel.stability_rating        =xlsx_table{row,tix.Stability};
    from_excel.grid_x                  =xlsx_table{row,tix.x};
    from_excel.grid_y                  =xlsx_table{row,tix.y};
    from_excel.electrode_depth         =xlsx_table{row,tix.Electrode_depth};
    from_excel.target                  =xlsx_table{row,tix.Target};
    from_excel.dataset                 =xlsx_table{row,tix.Set};
    from_excel.perturbation            =xlsx_table{row,tix.Perturbation};
    from_excel.perturbation_site       =xlsx_table{row,tix.Perturbation_site};
    % defining target including hemisphere (in case it wasnt used for target definition in the sorting table already)
    if ~isempty(tix.Hemisphere) && isempty(strfind(from_excel.target,'_')) && numel(xlsx_table{row,tix.Hemisphere})>0
        from_excel.target = [from_excel.target '_' upper(xlsx_table{row,tix.Hemisphere}(1))];
    end
end
end

