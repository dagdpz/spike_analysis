function pop_resorted = ph_accept_trials_per_unit(pop_resorted,trials,keys)
for u=1:numel(pop_resorted)
    p=pop_resorted(u);
    if isfield(p,'accepted')
        p=rmfield(p,'accepted');
    end
    pt=ph_get_unit_trials(p,trials);
    %    FRs=[pop_resorted(u).trial.FR_average];
    correct_task=ismember([pt.effector],keys.cal.effectors) ...
        & ismember([pt.type],keys.cal.types) ...
        & ismember([pt.completed],keys.cal.completed) ...
        & ~isnan([p.stability_rating]);
    for c=1:numel(keys.condition_parameters)
        par=keys.condition_parameters{c};
        if ~all(isnan(keys.cal.(par))) && ~isempty(keys.cal.(par))
            correct_task=correct_task & ismember([pt.(par)],keys.cal.(par));
        end
    end
    pop_resorted(u).accepted=correct_task;
    
end
end