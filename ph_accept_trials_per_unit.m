function pop_resorted = ph_accept_trials_per_unit(pop_resorted,trials,keys)
for u=1:numel(pop_resorted)
    p=pop_resorted(u);
    if isfield(p,'accepted')
        p=rmfield(p,'accepted');
    end
    pt=ph_get_unit_trials(p,trials);
    %    FRs=[pop_resorted(u).trial.FR_average];
    
    aborted_before_desired_state=false(size(pt));
    if numel(keys.cal.only_aborted_after_state)>1
        for r=1:size(keys.cal.only_aborted_after_state,1)
            aborted_before_desired_state_this_type=[pt.type]==keys.cal.only_aborted_after_state(r,1) & ...
              ~arrayfun(@(x) any(ismember(keys.cal.only_aborted_after_state(r,2),x.states)),pt);
              aborted_before_desired_state=aborted_before_desired_state | aborted_before_desired_state_this_type;
        end
    end
    
    correct_task=~aborted_before_desired_state...
        & ismember([pt.effector],keys.cal.effectors) ...
        & ismember([pt.type],keys.cal.types) ...
        & ismember([pt.completed],keys.cal.completed) ...
        & ~isnan([p.stability_rating]) ...
        & [p.exclusion_reason]==0 ;  %% stability including stability across blocks
    for c=1:numel(keys.condition_parameters)
        par=keys.condition_parameters{c};
        if ~all(isnan(keys.cal.(par))) && ~isempty(keys.cal.(par))
            correct_task=correct_task & ismember([pt.(par)],keys.cal.(par));
        end
    end
    pop_resorted(u).accepted=correct_task;
    
end
end