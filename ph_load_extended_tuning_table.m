function tuning_per_unit_table=ph_load_extended_tuning_table(keys)

%% Reading in table 

tuning_per_unit_table={};
load(keys.anova_table_file);

tuning_per_unit_table(:,end+1)=num2cell(repmat('Y',size(tuning_per_unit_table,1),1));
tuning_per_unit_table(1,end)={'ungrouped'};
taskcaseexistingindex=find(~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),'in_epoch_main'))); %not ideal
taskcases={};
for idx=taskcaseexistingindex
    task_existing_column=num2cell(~cellfun(@isempty,tuning_per_unit_table(:,idx)));
    task_existing_column{1,1}=['existing' tuning_per_unit_table{1,idx}(14:end)];
    tuning_per_unit_table(:,end+1)=task_existing_column;
    taskcases = [taskcases, {tuning_per_unit_table{1,idx}(15:end)}];
end
%% to be combined with previous 
taskcaseexistingindex=find(~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),'ch_epoch_main'))); %not ideal
%taskcases={};
for idx=taskcaseexistingindex
    task_existing_column=num2cell(~cellfun(@isempty,tuning_per_unit_table(:,idx)));
    task_existing_column{1,1}=['existing_' tuning_per_unit_table{1,idx}(1:2) tuning_per_unit_table{1,idx}(14:end)];
    tuning_per_unit_table(:,end+1)=task_existing_column;
    %taskcases = [taskcases, {tuning_per_unit_table{1,idx}(15:end)}];
end

combined_column=cell(size(tuning_per_unit_table,1),1);
for props=2:size(keys.tt.combine_tuning_properties,2)
    column_index=find_column_index(tuning_per_unit_table,keys.tt.combine_tuning_properties{1,props});
    if isempty(column_index)
        disp(['column not found: ' keys.tt.combine_tuning_properties{sel,1}]);
        continue;
    end
    combined_column=cellfun(@(x,y) cat(2,x,num2str(y)),combined_column,tuning_per_unit_table(:,column_index),'uniformoutput',false);
end
if ~isempty(keys.tt.combine_tuning_properties)
combined_column{1,1}=keys.tt.combine_tuning_properties{1};
else
combined_column{1,1}='combined_effects';    
end
tuning_per_unit_table(:,end+1)=combined_column;

%% create space OR interaction column
for t=1:numel(taskcases)
    taskcase=taskcases{t};
    idx_space   =find_column_index(tuning_per_unit_table,[keys.instructed_choice '_spaceLR_main_'   taskcase] );
    idx_epoch   =find_column_index(tuning_per_unit_table,[keys.instructed_choice '_epoch_main_'     taskcase] );
    idx_hands   =find_column_index(tuning_per_unit_table,[keys.instructed_choice '_hands_main_'     taskcase] );
    idx_SxH     =find_column_index(tuning_per_unit_table,[keys.instructed_choice '_SxH_'            taskcase] );
    idx_ExS     =find_column_index(tuning_per_unit_table,[keys.instructed_choice '_ExS_'            taskcase] );
    idx_ExH     =find_column_index(tuning_per_unit_table,[keys.instructed_choice '_ExH_'            taskcase] );
    idx_ESH     =find_column_index(tuning_per_unit_table,[keys.instructed_choice '_ExSxH_'          taskcase] );
    
    
if ~isempty(idx_epoch); tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx_epoch))],idx_epoch)  ={NaN}; end;
if ~isempty(idx_space); tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx_space))],idx_space)  ={''}; end;
if ~isempty(idx_hands); tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx_hands))],idx_hands)  ={''}; end;
if ~isempty(idx_SxH);   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx_SxH))],idx_SxH)      ={''}; end;
if ~isempty(idx_ExS);   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx_ExS))],idx_ExS)      ={NaN}; end;
if ~isempty(idx_ExH);   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx_ExH))],idx_ExH)      ={NaN}; end;
if ~isempty(idx_ESH);   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx_ESH))],idx_ESH)      ={NaN}; end;
    
    
    switch keys.tt.space_criterion %% consider only neurons that show 1) SxE or space main effect; 2) SxE; 3) take all (no criterion)?
        case 'interaction or space only'
            space_or_interaction=([tuning_per_unit_table{2:end,idx_ExS}]==1)' | ~ismember(tuning_per_unit_table(2:end,idx_space),'-');
        case 'interaction only'
            space_or_interaction=([tuning_per_unit_table{2:end,idx_ExS}]==1)';
        case 'none'
            space_or_interaction=true(size(tuning_per_unit_table,1)-1,1);
    end
    %%% goofy
    switch keys.tt.epoch_criterion %% consider only neurons that show 1) SxE or space main effect; 2) SxE; 3) take all (no criterion)?
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
    switch keys.tt.hands_criterion %% consider only neurons that show 1) HxE or hand main effect; 2) HxE; 3) take all (no criterion)?
        case 'interaction or hands only'
            hands_or_interaction=([tuning_per_unit_table{2:end,idx_ExH}]==1)' | ~ismember(tuning_per_unit_table(2:end,idx_hands),'-');
        case 'interaction only'
            hands_or_interaction=([tuning_per_unit_table{2:end,idx_ExH}]==1)';
        case 'none'
            hands_or_interaction=true(size(tuning_per_unit_table,1)-1,1);
    end
    switch keys.tt.SXH_criterion %% consider only neurons that show 1) SxHxE or SxH; 2) SxHxE; 3) take all (no criterion)?
         case 'SXH SXE HXE SHE'
             SXH_or_interaction=([tuning_per_unit_table{2:end,idx_ESH}]==1)' | ~ismember(tuning_per_unit_table(2:end,idx_SxH),'-')...
                               |([tuning_per_unit_table{2:end,idx_ExH}]==1)' | ([tuning_per_unit_table{2:end,idx_ExS}]==1)';
        case 'interaction or SXH only'
            SXH_or_interaction=([tuning_per_unit_table{2:end,idx_ESH}]==1)' | ~ismember(tuning_per_unit_table(2:end,idx_SxH),'-');
        case 'interaction only'
            SXH_or_interaction=([tuning_per_unit_table{2:end,idx_ESH}]==1)';
        case 'none'
            SXH_or_interaction=true(size(tuning_per_unit_table,1)-1,1);
    end
    n_column=size(tuning_per_unit_table,2);
    
    n_column=n_column+1;
    tuning_per_unit_table{1,n_column}=['space_or_interaction_' taskcase];
    tuning_per_unit_table(2:end,n_column)=num2cell(space_or_interaction);
    
    n_column=n_column+1;
    tuning_per_unit_table{1,n_column}=['epoch_or_interaction_' taskcase];
    tuning_per_unit_table(2:end,n_column)=num2cell(epoch_or_interaction);
    
    n_column=n_column+1;
    tuning_per_unit_table{1,n_column}=['hands_or_interaction_' taskcase];
    tuning_per_unit_table(2:end,n_column)=num2cell(hands_or_interaction);
    
    n_column=n_column+1;
    tuning_per_unit_table{1,n_column}=['SXH_or_interaction_' taskcase];
    tuning_per_unit_table(2:end,n_column)=num2cell(SXH_or_interaction);
end
end