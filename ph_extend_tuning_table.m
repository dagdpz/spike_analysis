function ph_extend_tuning_table(keys)

%% Reading in table 

load(keys.anova_table_file);
tuning_per_unit_table=LR_to_CI(tuning_per_unit_table);

%% Reducing table according to target
% tar_idx=strcmp(tuning_per_unit_table(1,:),'target');
% to_remove=false(size(tuning_per_unit_table,1),1);
% for k=2:size(tuning_per_unit_table,1)
%     if ~isempty(keys.target)
%         to_remove(k)=isempty(strfind(tuning_per_unit_table{k,tar_idx},keys.target)) ;
%     end
% end

%% create space OR interaction column
for d=1:numel(keys.cc.tasktypes)
    keys.dataset=keys.cc.tasktypes{d};
    idx_space   =find_column_index(tuning_per_unit_table,[keys.instructed_choice '_spaceLR_main_'   keys.dataset '_' keys.case ] );
    idx_epoch   =find_column_index(tuning_per_unit_table,[keys.instructed_choice '_epoch_main_'     keys.dataset '_' keys.case ] );
    idx_hands   =find_column_index(tuning_per_unit_table,[keys.instructed_choice '_hands_main_'     keys.dataset '_' keys.case ] );
    idx_SxH     =find_column_index(tuning_per_unit_table,[keys.instructed_choice '_SxH_'            keys.dataset '_' keys.case] );
    idx_ExS     =find_column_index(tuning_per_unit_table,[keys.instructed_choice '_ExS_'            keys.dataset '_' keys.case] );
    idx_ExH     =find_column_index(tuning_per_unit_table,[keys.instructed_choice '_ExH_'            keys.dataset '_' keys.case] );
    idx_ESH     =find_column_index(tuning_per_unit_table,[keys.instructed_choice '_ExSxH_'          keys.dataset '_' keys.case] );
    
    switch keys.cc.space_criterion %% consider only neurons that show 1) SxE or space main effect; 2) SxE; 3) take all (no criterion)?
        case 'interaction or space only'
            space_or_interaction=([tuning_per_unit_table{2:end,idx_ExS}]==1)' | ~ismember(tuning_per_unit_table(2:end,idx_space),'-');
        case 'interaction only'
            space_or_interaction=([tuning_per_unit_table{2:end,idx_ExS}]==1)';
        case 'none'
            space_or_interaction=true(size(tuning_per_unit_table,1)-1,1);
    end
    %%% goofy
    switch keys.cc.epoch_criterion %% consider only neurons that show 1) SxE or space main effect; 2) SxE; 3) take all (no criterion)?
        case 'SxE or epoch only'
            epoch_or_interaction=([tuning_per_unit_table{2:end,idx_ExS}]==1)' | ([tuning_per_unit_table{2:end,idx_epoch}]==1)';
        case 'HxE or epoch only'
            epoch_or_interaction=([tuning_per_unit_table{2:end,idx_ExH}]==1)' | ([tuning_per_unit_table{2:end,idx_epoch}]==1)';
        case 'SHE or epoch only'
            epoch_or_interaction=([tuning_per_unit_table{2:end,idx_ESH}]==1)' | ([tuning_per_unit_table{2:end,idx_epoch}]==1)';
        case 'SxE only'
            epoch_or_interaction=([tuning_per_unit_table{2:end,idx_ExS}]==1)';
        case 'HxE only'
            epoch_or_interaction=([tuning_per_unit_table{2:end,idx_ExH}]==1)';
        case 'SHE only'
            epoch_or_interaction=([tuning_per_unit_table{2:end,idx_ESH}]==1)';
        case 'none'
            epoch_or_interaction=true(size(tuning_per_unit_table,1)-1,1);
    end
    switch keys.cc.hands_criterion %% consider only neurons that show 1) HxE or hand main effect; 2) HxE; 3) take all (no criterion)?
        case 'interaction or hands only'
            hands_or_interaction=([tuning_per_unit_table{2:end,idx_ExH}]==1)' | ~ismember(tuning_per_unit_table(2:end,idx_hands),'-');
        case 'interaction only'
            hands_or_interaction=([tuning_per_unit_table{2:end,idx_ExH}]==1)';
        case 'none'
            hands_or_interaction=true(size(tuning_per_unit_table,1)-1,1);
    end
    switch keys.cc.SXH_criterion %% consider only neurons that show 1) SxHxE or SxH; 2) SxHxE; 3) take all (no criterion)?
        case 'interaction or SXH only'
            SXH_or_interaction=([tuning_per_unit_table{2:end,idx_ESH}]==1)' | ~ismember(tuning_per_unit_table(2:end,idx_SxH),'-');
        case 'interaction only'
            SXH_or_interaction=([tuning_per_unit_table{2:end,idx_ESH}]==1)';
        case 'none'
            SXH_or_interaction=true(size(tuning_per_unit_table,1)-1,1);
    end
    n_column=size(tuning_per_unit_table,2);
    
    n_column=n_column+1;
    tuning_per_unit_table{1,n_column}=['space_or_interaction_' keys.dataset];
    tuning_per_unit_table(2:end,n_column)=num2cell(space_or_interaction);
    %tuning_per_unit_table(2:end,n_column)=cellstr(num2str(NaN(size(space_or_interaction))));
    
    n_column=n_column+1;
    tuning_per_unit_table{1,n_column}=['epoch_or_interaction_' keys.dataset];
    tuning_per_unit_table(2:end,n_column)=num2cell(epoch_or_interaction);
    %tuning_per_unit_table(2:end,n_column)=cellstr(num2str(NaN(size(epoch_or_interaction))));
    
    n_column=n_column+1;
    tuning_per_unit_table{1,n_column}=['hands_or_interaction_' keys.dataset];
    tuning_per_unit_table(2:end,n_column)=num2cell(hands_or_interaction);
    %tuning_per_unit_table(2:end,n_column)=cellstr(num2str(NaN(size(hands_or_interaction))));
    
    n_column=n_column+1;
    tuning_per_unit_table{1,n_column}=['SXH_or_interaction_' keys.dataset];
    tuning_per_unit_table(2:end,n_column)=num2cell(SXH_or_interaction);
    %tuning_per_unit_table(2:end,n_column)=cellstr(num2str(NaN(size(SXH_or_interaction))));
    
%     %% Reducing table according to keys.cc.only_both_hands and Selection
%     for k=2:size(tuning_per_unit_table,1)
%         if keys.cc.only_both_hands && any(idx_hands)
%             to_remove(k)=to_remove(k) || strcmp(tuning_per_unit_table(k,idx_hands),'');
%         end
%     end
%     
%     for sel=1:size(keys.cc.Selection,1)
%         selidx=find_column_index(tuning_per_unit_table,keys.cc.Selection{sel,1});
%         for k=2:size(tuning_per_unit_table,1)
%             if ischar(keys.cc.Selection{sel,2})
%                 to_remove(k)=to_remove(k) || isempty(strfind(tuning_per_unit_table{k,selidx},keys.cc.Selection{sel,2}));
%             elseif tuning_per_unit_table{k,selidx} ~= keys.cc.Selection{sel,2} %% find a way of not only selecting one specific value, but a range
%                 to_remove(k)=true;
%             end
%         end
%     end
%     
%     %% reducing to the neurons that have this dataset !!
%     to_remove(cellfun(@(x) isempty(x),tuning_per_unit_table(:,idx_epoch)))=true;
end

keys.xlsx_table=tuning_per_unit_table(~to_remove,:);
keys.space_or_interaction=space_or_interaction(~to_remove(2:end));
keys.epoch_or_interaction=epoch_or_interaction(~to_remove(2:end));
keys.hands_or_interaction=hands_or_interaction(~to_remove(2:end));
keys.SXH_or_interaction=SXH_or_interaction(~to_remove(2:end));
end


% function xlsx_table = get_mat_table(keys)
% xlsx_table={};
% xlsx_table.RAWW=tuning_per_unit_table;
% xlsx_table.STRR=tuning_per_unit_table;
% xlsx_table.STRR(cellfun(@isempty,xlsx_table.STRR))={''};
% xlsx_table.RAWW(cellfun(@isempty,xlsx_table.STRR))={NaN}; %%??
% end