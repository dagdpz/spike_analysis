function pop_resorted = sort_by_unit_ID_temp(o_t)
fields_to_remove={'unit','channel','TDT_CAP1','TDT_POX1','TDT_ECG1','TDT_ECG4','TDT_LFPx'};
fields_to_keep={'block','run','n','dataset','perturbation','FR_average','stability_rating','SNR_rating'};
ttt=[fields_to_keep; cell(size(fields_to_keep))];
pop_resorted=struct('unit_ID',{},'trial',{},'waveforms',{},'unit_SNR_rating',{},'Single_rating',{},'unit_stability_rating',{},'block_unit',{},ttt{:});
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
        
        trial_fields=fieldnames(o_t(b).trial);
        %trial_fieldnames_to_remove=fields_to_remove(ismember(fields_to_remove,fieldnames(o_t(b).trial)));
        %trial_fields_to_remove=trial_fields(~ismember(trial_fields,fields_to_keep));
        %tmp=rmfield(o_t(b).trial,trial_fields_to_remove);
        tmp=rmfield(o_t(b).trial,trial_fields);
        %tmp=struct(current_unit);
        
        [tmp.arrival_times]=current_unit.arrival_times;
        [tmp.waveforms]=current_unit.waveforms;
        
        for fn=fields_to_keep
           pop_resorted(U).(fn{:})= [pop_resorted(U).(fn{:}) current_unit.dataset];
        end
%         
%         [tmp.dataset]=current_unit.dataset;
%         [tmp.perturbation]=current_unit.perturbation;
%         [tmp.FR_average]=current_unit.FR_average;
%         [tmp.stability_rating]=current_unit.stability_rating;
%         [tmp.SNR_rating]=current_unit.SNR_rating;
%         
%         amp_wf=arrayfun(@(x) max(mean(x.waveforms,1)-min(mean(x.waveforms,1))),current_unit,'uniformoutput',false);
%         std_wf=arrayfun(@(x) nanmean(std(x.waveforms,1)),current_unit,'uniformoutput',false);
%         snr_wf=cellfun(@(x,y) x/y,amp_wf,std_wf,'uniformoutput',false);
%         [tmp.waveform_amplitude]=deal(amp_wf{:});
%         [tmp.waveform_std]      =deal(std_wf{:});
%         [tmp.trialwise_SNR]     =deal(snr_wf{:});
        
        
        amp_wf=arrayfun(@(x) max(mean(x.waveforms,1)-min(mean(x.waveforms,1))),current_unit,'uniformoutput',false);
        std_wf=arrayfun(@(x) nanmean(std(x.waveforms,1)),current_unit,'uniformoutput',false);
        snr_wf=cellfun(@(x,y) x/y,amp_wf,std_wf,'uniformoutput',false);
%         [pop_resorted(U).waveform_amplitude]=deal(amp_wf{:});
%         [pop_resorted(U).waveform_std]      =deal(std_wf{:});
%         [pop_resorted(U).trialwise_SNR]     =deal(snr_wf{:});
        
        
        pop_resorted(U).waveform_amplitude=[amp_wf{:}];
        pop_resorted(U).waveform_std      =[std_wf{:}];
        pop_resorted(U).trialwise_SNR     =[snr_wf{:}];
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
