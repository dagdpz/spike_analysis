function population=ph_assign_perturbation_group(keys,population)
for u=1:numel(population)
    if ~isfield(population(u).trial,'perturbation') % quickfix to load microstim data
    [population(u).trial.perturbation]=deal(0);
       continue 
    end
    nanidx=true(size([population(u).trial.perturbation]));
    for g=1:numel(keys.cal.perturbation_groups)
        idx=ismember([population(u).trial.perturbation],keys.cal.perturbation_groups{g});
        [population(u).trial(idx).perturbation]=deal(g-1);
        nanidx=nanidx&~idx;
    end
    %     idx0=ismember([population(u).trial.perturbation],keys.cal.perturbation_groups{1});
    %     idx1=ismember([population(u).trial.perturbation],keys.cal.perturbation_groups{2});
    %     [population(u).trial(idx0).perturbation]=deal(0);
    %     [population(u).trial(idx1).perturbation]=deal(1);
    [population(u).trial(nanidx).perturbation]=deal(NaN);
end
end