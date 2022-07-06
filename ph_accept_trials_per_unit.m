function pop_resorted = ph_accept_trials_per_unit(pop_resorted,keys)
for u=1:numel(pop_resorted)
    FRs=[pop_resorted(u).trial.FR_average];
    FRsT=FRs;
    if keys.cal.remove_trials_without_spikes
        FRsT=FRsT(FRsT~=0); % this has to be commented for pulv_oculomotor!!
    end
    unit_mean=double(nanmedian(FRsT));
    unit_std=double(nanstd(FRsT));
    confidence_interval=3*unit_std; %poisson?
    
    % here is the actual criterion
    accepted=num2cell(FRs>(log(1+exp(unit_mean-confidence_interval))) & FRs<(log(1+exp(unit_mean+confidence_interval))));
    [pop_resorted(u).trial.accepted]=deal(accepted{:});
    %[pop_resorted(u).trial.accepted]=deal(true);
end
end