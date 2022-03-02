function pop=ph_target_reassign(keys,pop)

col_index_ref = find(strcmp(pop(1,:),keys.contra_ipsi_relative_to));
col_index_target = find(strcmp(pop(1,:),'target'));
for i = 2:size(pop,1)
    ref = char(pop(i,col_index_ref));
    target = char(pop(i,col_index_target));
    if strcmpi(ref(end-1:end),'_L') & strcmp(keys.contra_ipsi_relative_to,'perturbation_site')
        
        if strcmpi(target(end-1:end),'_L')
            pop(i,col_index_target) = {[target(1:end-1) 'R']};
        elseif  strcmpi(target(end-1:end),'_R')
            pop(i,col_index_target) = {[target(1:end-1) 'L']};
        end
    end
end
end