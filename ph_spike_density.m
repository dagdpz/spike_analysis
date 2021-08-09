function [SD  bins SD_VAR]=ph_spike_density(trial,wn,keys,baseline,norm_factor)
sta=keys.PSTH_WINDOWS{wn,2};
t_before_state=keys.PSTH_WINDOWS{wn,3};
t_after_state=keys.PSTH_WINDOWS{wn,4};
bins=t_before_state:keys.PSTH_binwidth:t_after_state;

if ~any([trial(:).states]==sta)
    SD=NaN(1,numel(bins));
end
%
% n_trials=0;
% PSTH_ms =t_before_state-1:0.001:t_after_state+1;
% t_idx=PSTH_ms+0.0001>=t_before_state & PSTH_ms<=t_after_state + keys.PSTH_binwidth -0.0001;
% switch keys.kernel_type
%     case 'gaussian'
%         Kernel=normpdf(-5*keys.gaussian_kernel:0.001:5*keys.gaussian_kernel,0,keys.gaussian_kernel);
%     case 'box'
%         n_bins=round(2*keys.gaussian_kernel/0.001);
%         Kernel=ones(1,n_bins)/n_bins*1000; %%*1000 cause a one full spike in one 1ms bin means 1000sp/s locally
% end
%
% for t=1:numel(trial)
%     idx_sta=trial(t).states==sta;
%     if any(idx_sta)
%         arrival_times=trial(t).arrival_times-trial(t).states_onset(idx_sta);
%         n_trials=n_trials+1;
%         SD_ms(n_trials,:)= (conv(hist(arrival_times,PSTH_ms),Kernel,'same')-baseline(t))/norm_factor(t);
%     end
% end
%
% if n_trials==0
%     SD=NaN(1,sum(t_idx)/keys.PSTH_binwidth*0.001);
%     SD_VAR=NaN(1,sum(t_idx)/keys.PSTH_binwidth*0.001);
% else
%     SD=sum(SD_ms,1)/n_trials;
%     %SD_SEM=sterr(SD_ms,1);
%     SD_VAR=var(SD_ms,0,1);
%     SD_VAR=nanmean(reshape(SD_VAR(t_idx),keys.PSTH_binwidth/0.001,sum(t_idx)/keys.PSTH_binwidth*0.001),1);
%     SD=nanmean(reshape(SD(t_idx),keys.PSTH_binwidth/0.001,sum(t_idx)/keys.PSTH_binwidth*0.001),1);
% end


%
%
% n_trials=0;
% PSTH_ms =t_before_state-1:keys.PSTH_binwidth:t_after_state+1;
% t_idx=PSTH_ms+0.0001>=t_before_state & PSTH_ms<=t_after_state + keys.PSTH_binwidth -0.0001;
% switch keys.kernel_type
%     case 'gaussian'
%         Kernel=normpdf(-5*keys.gaussian_kernel:keys.PSTH_binwidth:5*keys.gaussian_kernel,0,keys.gaussian_kernel);
%     case 'box'
%         n_bins=round(2*keys.gaussian_kernel/keys.PSTH_binwidth);
%         Kernel=ones(1,n_bins)/n_bins/keys.PSTH_binwidth; %%*1000 cause a one full spike in one 1ms bin means 1000sp/s locally
% end
%
% for t=1:numel(trial)
%     idx_sta=trial(t).states==sta;
%     if any(idx_sta)
%         arrival_times=trial(t).arrival_times-trial(t).states_onset(idx_sta);
%         n_trials=n_trials+1;
%         SD_ms(n_trials,:)= (conv(hist(arrival_times,PSTH_ms),Kernel,'same')-baseline(t))/norm_factor(t);
%     end
% end
%
% if n_trials==0
%     SD=NaN(1,sum(t_idx));
%     SD_VAR=NaN(1,sum(t_idx));
% else
%     SD=sum(SD_ms,1)/n_trials;
%     %SD_SEM=sterr(SD_ms,1);
%     SD_VAR=var(SD_ms,0,1);
% %     SD_VAR=nanmean(reshape(SD_VAR(t_idx),keys.PSTH_binwidth/0.001,sum(t_idx)/keys.PSTH_binwidth*0.001),1);
% %     SD=nanmean(reshape(SD(t_idx),keys.PSTH_binwidth/0.001,sum(t_idx)/keys.PSTH_binwidth*0.001),1);
% end


SD=NaN(1,numel(bins));
SD_VAR=NaN(1,numel(bins));

n_trials=0;
PSTH_ms =t_before_state-1:0.001:t_after_state+1;
t_idx=PSTH_ms+0.0001>=t_before_state & PSTH_ms<=t_after_state + keys.PSTH_binwidth -0.0001;

switch keys.kernel_type
    case 'gaussian'
        Kernel=normpdf(-5*keys.gaussian_kernel:0.001:5*keys.gaussian_kernel,0,keys.gaussian_kernel);
    case 'box'
        n_bins=round(2*keys.gaussian_kernel/0.001);
        Kernel=ones(1,n_bins)/n_bins*1000; %%*1000 cause a one full spike in one 1ms bin means 1000sp/s locally
end

if nargout>2 || (numel(baseline>=2) && (any(diff(baseline)) || any(diff(norm_factor))))
    SD_ms=zeros(numel(trial),numel(PSTH_ms));
    for t=1:numel(trial)
        idx_sta=trial(t).states==sta;
        if any(idx_sta)
            arrival_times=trial(t).arrival_times-trial(t).states_onset(idx_sta);
            n_trials=n_trials+1;
            SD_ms(n_trials,:)= (conv(hist(arrival_times,PSTH_ms),Kernel,'same')-baseline(t))/norm_factor(t);
        end
    end
    SD_ms=SD_ms(1:n_trials,:);
    if n_trials==0
        SD=NaN(1,sum(t_idx)/keys.PSTH_binwidth*0.001);
        SD_VAR=NaN(1,sum(t_idx)/keys.PSTH_binwidth*0.001);
    else
        SD=sum(SD_ms,1)/n_trials;
        %SD_SEM=sterr(SD_ms,1);
        SD_VAR=var(SD_ms,0,1);
        SD_VAR=nanmean(reshape(SD_VAR(t_idx),keys.PSTH_binwidth/0.001,sum(t_idx)/keys.PSTH_binwidth*0.001),1);
        SD=nanmean(reshape(SD(t_idx),keys.PSTH_binwidth/0.001,sum(t_idx)/keys.PSTH_binwidth*0.001),1);
    end
else
    arrival_times=[];
    for t=1:numel(trial)
        idx_sta=trial(t).states==sta;
        if any(idx_sta)
            arrival_times=[arrival_times; trial(t).arrival_times-trial(t).states_onset(idx_sta)];
            n_trials=n_trials+1;
        end
    end
    if n_trials==0
        SD=NaN(1,sum(t_idx)/keys.PSTH_binwidth*0.001);
    else
        SD_ms(1,:)= conv(hist(arrival_times,PSTH_ms),Kernel,'same');
        SD=(SD_ms/n_trials-baseline(1))/norm_factor(1);
        SD=nanmean(reshape(SD(t_idx),keys.PSTH_binwidth/0.001,sum(t_idx)/keys.PSTH_binwidth*0.001),1);
    end
end
aa=1;
end