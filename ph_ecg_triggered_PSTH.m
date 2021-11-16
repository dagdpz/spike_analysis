% Todo - block conditions
% several units

load('Y:\Projects\Pulv_distractor_spatial_choice\Data\Bacchus\ECG\20211001\20211001_ecg.mat')
load('Y:\Projects\Pulv_distractor_spatial_choice\ephys\StimulusType_Difficulty_Position_LS\population_Bacchus_20211001.mat')

ECG_event=-1;
keys.PSTH_WINDOWS={'ECG',ECG_event,-0.5,0.5};
keys.PSTH_binwidth=0.01;
keys.kernel_type='gaussian';
keys.gaussian_kernel=0.02;

for u=1:numel(population)
    pop=population(u);
    trcell=num2cell(pop.trial);
    arrival_times=cellfun(@(x) [x.arrival_times]+x.trial_onset_time,trcell,'Uniformoutput',false);
    
    for o=1:numel(out)
        block=out(o).nrblock;
        tr=[pop.trial.block]==block;
        AT=vertcat(arrival_times{tr});
        RPEAK_ts=[out(o).Rpeak_t];
        %% reduce RPEAK_ts potentially ? (f.e.: longer than recorded ephys, inter-trial-interval?)
        
        
        %% now the tricky part: sort by ECG peaks ...
        for t=1:numel(RPEAK_ts)
            unit(u).block(o).trial(t).states=ECG_event;
            unit(u).block(o).trial(t).states_onset=0;
            AT_temp=AT-RPEAK_ts(t);
            unit(u).block(o).trial(t).arrival_times=AT_temp(AT_temp>keys.PSTH_WINDOWS{1,3}-0.2 & AT_temp<keys.PSTH_WINDOWS{1,4}+0.2);
        end
        
        trial=unit(u).block(o).trial;
        [SD  bins SD_VAR SD_SEM]=ph_spike_density(trial,1,keys,zeros(size(trial)),ones(size(trial)));
        
        plot((keys.PSTH_WINDOWS{1,3}:keys.PSTH_binwidth:keys.PSTH_WINDOWS{1,4})*1000,SD);
        
    end
    
    
    
end


