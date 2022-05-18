function TT=ph_target_reassign(keys,TT)

col_index_ref = find(strcmp(TT(1,:),keys.contra_ipsi_relative_to));
col_index_target = find(strcmp(TT(1,:),'target'));
for i = 2:size(TT,1)
    ref = char(TT(i,col_index_ref));
    target = char(TT(i,col_index_target));
    if strcmp(keys.contra_ipsi_relative_to,'perturbation_site')
        target_right=strcmpi(target(end-1:end),'_R') || strcmpi(target(end-1:end),'_r');
        ref_right=strcmpi(ref(end-1:end),'_R') || strcmpi(ref(end-1:end),'_r');
        
        if (target_right && ref_right) || (~target_right && ~ref_right)
            TT(i,col_index_target) = {[target(1:end-1) 'L']}; %% left becomes ipsi!
        elseif (target_right && ~ref_right) || (~target_right && ref_right)
            TT(i,col_index_target) = {[target(1:end-1) 'R']}; %% Right becomes contra!
        end
    end
end
end