function [SD  bins]=ph_spike_density(trial,wn,keys,baseline,norm_factor)
%keys.gaussian_kernel=0.02;
sta=keys.PSTH_WINDOWS{wn,2};
t_before_state=keys.PSTH_WINDOWS{wn,3};
t_after_state=keys.PSTH_WINDOWS{wn,4};
bins=t_before_state:keys.PSTH_binwidth:t_after_state;

if ~any([trial(:).states]==sta)
    SD=NaN(1,numel(bins));
end

arrival_times=[];
n_trials=0;
PSTH_ms =t_before_state-1:0.001:t_after_state+1;
t_idx=PSTH_ms+0.0001>=t_before_state & PSTH_ms<=t_after_state + keys.PSTH_binwidth -0.0001;
gaussian=normpdf(-5*keys.gaussian_kernel:0.001:5*keys.gaussian_kernel,0,keys.gaussian_kernel);

%% how to treat trials that dont contain the state to which we align?
% for t=1:numel(trial)
%     idx=trial(t).states==sta;
%     if any(idx)
%         arrival_times=[arrival_times; trial(t).arrival_times-trial(t).states_onset(idx)];
%         n_trials=n_trials+1;
%     end
% end
% SD_ms= conv(hist(arrival_times,PSTH_ms),normpdf(-5*keys.gaussian_kernel:0.001:5*keys.gaussian_kernel,0,keys.gaussian_kernel),'same')/max([n_trials 1]);
% t_idx=PSTH_ms+0.0001>=t_before_state & PSTH_ms<=t_after_state + keys.PSTH_binwidth -0.0001;
% SD=nanmean(reshape(SD_ms(t_idx),keys.PSTH_binwidth/0.001,sum(t_idx)/keys.PSTH_binwidth*0.001),1);
% if n_trials==0
%     SD=[];
% else
%     SD=(SD-baseline)/norm_factor;
% end



%try out zscore only on windows of analysis
% if strcmp(keys.PO.normalization,'z_score')
%     for t=1:numel(trial)
%         idx_sta=trial(t).states==sta;
%         if any(idx_sta)
%             arrival_times=trial(t).arrival_times-trial(t).states_onset(idx_sta);
%             n_trials=n_trials+1;
%             SD_ms(n_trials,:)= conv(hist(arrival_times,PSTH_ms),gaussian,'same');
%             
%         end
%     end
%     
%     if n_trials==0
%         SD=NaN(1,sum(t_idx)/keys.PSTH_binwidth*0.001);  
%     else
%         baseline =  nanmean(reshape(SD_ms,[1 size(SD_ms,1)*size(SD_ms,2)]));
%         SD_ms = SD_ms - baseline;
%         SD=sum(SD_ms,1)/n_trials;
%         SD=nanmean(reshape(SD(t_idx),keys.PSTH_binwidth/0.001,sum(t_idx)/keys.PSTH_binwidth*0.001),1);
%         norm_std = nanstd(reshape(SD_ms,[1 size(SD_ms,1)*size(SD_ms,2)]));;
%         SD=SD/norm_std;
%     end
%       
%     
% else
    for t=1:numel(trial)
        idx_sta=trial(t).states==sta;
        if any(idx_sta)
            arrival_times=trial(t).arrival_times-trial(t).states_onset(idx_sta);
            n_trials=n_trials+1;
            SD_ms(n_trials,:)= conv(hist(arrival_times,PSTH_ms),gaussian,'same')-baseline(t);
        end
    end
    
    if n_trials==0
        SD=NaN(1,sum(t_idx)/keys.PSTH_binwidth*0.001);
    else
        SD=sum(SD_ms,1)/n_trials;
        SD=nanmean(reshape(SD(t_idx),keys.PSTH_binwidth/0.001,sum(t_idx)/keys.PSTH_binwidth*0.001),1);
%         SD=SD/norm_factor;
        SD=(SD-norm_factor)/norm_factor;
    end
% end
%
% SD=single(resample(SD_ms(t_idx),1,keys.PSTH_binwidth/0.001));
%
% SD=histc(SD_ms,1:keys.PSTH_binwidth/0.001:sum(t_idx));

end