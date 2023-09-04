function trials=ph_assign_perturbation_group(keys,trials)
nanidx=true(size(trials));
for g=1:numel(keys.cal.perturbation_groups)
    idx=ismember([trials.perturbation],keys.cal.perturbation_groups{g});
    [trials(idx).perturbation]=deal(g-1);
    nanidx=nanidx&~idx;
end
[trials(nanidx).perturbation]=deal(NaN);   

end