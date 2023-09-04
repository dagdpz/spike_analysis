function t=ph_get_unit_trials(u,trials)
m=ismember({trials.monkey},u.monkey) & ismember([trials.date],u.date);
t=trials(m);
trial_id_u=[u.block' u.run' u.n'];
trial_id_t=[[t.block]' [t.run]' [t.n]'];
u_trials=ismember(trial_id_t,trial_id_u,'rows');
t=t(u_trials);

u_trials=ismember(trial_id_u,trial_id_t,'rows');

%% add more?
fields_to_take_over={'accepted'};
for fn=fields_to_take_over
    if isfield(u,fn{:})
    current_field=num2cell(u.(fn{:})(u_trials));
    [t.(fn{:})]=deal(current_field{:});
    end
end

end