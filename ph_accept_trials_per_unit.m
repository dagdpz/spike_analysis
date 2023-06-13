function pop_resorted = ph_accept_trials_per_unit(pop_resorted,keys)
for u=1:numel(pop_resorted)
    FRs=[pop_resorted(u).trial.FR_average];
    correct_task=ismember([pop_resorted(u).trial.effector],keys.cal.effectors) ...
               & ismember([pop_resorted(u).trial.type],keys.cal.types) ...
               & ismember([pop_resorted(u).trial.completed],keys.cal.completed);
    for c=1:numel(keys.condition_parameters)
        par=keys.condition_parameters{c};
        if ~all(isnan(keys.cal.(par))) && ~isempty(keys.cal.(par))
            correct_task=correct_task & ismember([pop_resorted(u).trial.(par)],keys.cal.(par));
        end
    end
    FRsT=FRs(correct_task);
    if keys.cal.remove_trials_without_spikes % pulv_oculomotor!!
        FRsT=FRsT(FRsT~=0); 
    end
    unit_mean=double(nanmedian(FRsT));
    unit_std=double(nanstd(FRsT));
    confidence_interval=3*unit_std; %poisson?
    
    % here is the actual outliar criterion
    if keys.cal.remove_trials_with_outlying_FR
        accepted=num2cell(FRs>(log(1+exp(unit_mean-confidence_interval))) & FRs<(log(1+exp(unit_mean+confidence_interval))) & correct_task);
    else
        accepted=num2cell(correct_task);
    end
    [pop_resorted(u).trial.accepted]=deal(accepted{:});
end
end