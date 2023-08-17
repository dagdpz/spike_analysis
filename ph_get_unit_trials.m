function pt=ph_get_unit_trials(p,trials)
m=ismember({trials.monkey},p.monkey) & ismember([trials.date],p.date);
pt=trials(m);
trial_id_u=[p.block' p.run' p.n'];
trial_id_t=[[pt.block]' [pt.run]' [pt.n]'];
u_trials=ismember(trial_id_t,trial_id_u,'rows');
pt=pt(u_trials);

u_trials=ismember(trial_id_u,trial_id_t,'rows');

fields_to_take_over={'dataset','perturbation','accepted'};
for fn=fields_to_take_over
    current_field=num2cell(p.(fn{:})(u_trials));
    [pt.(fn{:})]=deal(current_field{:});
end
end