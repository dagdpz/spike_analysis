function [tuning_per_unit_table Sel_for_title]=ph_reduce_tuning_table(tuning_per_unit_table,keys)

%% monkey
idx_unitID=find_column_index(tuning_per_unit_table,'unit_ID');
if numel(keys.monkey)>3
    tuning_per_unit_table=tuning_per_unit_table([true; ~cellfun(@isempty,strfind(tuning_per_unit_table(2:end,idx_unitID),keys.monkey(1:3)))],:);
end

%% overlapping tasktypes_hands_choices
tasktypes_hands_choices=combvec((1:numel(keys.tt.tasktypes)),keys.tt.hands,keys.tt.choices,keys.tt.perturbations); % what is going on here???
cell_idx_tuning_table=true(size(tuning_per_unit_table,1)-1,1);
for t=1:size(tasktypes_hands_choices,2) %% crashes for not existing monkey
    tasktype=tasktypes_hands_choices(1,t);
    hand    =tasktypes_hands_choices(2,t)+1;
    choice  =tasktypes_hands_choices(3,t)+1;
    perturbation  =tasktypes_hands_choices(4,t)+1;
    column_title=['existing_' keys.labels.choices{choice} '_' keys.labels.handsIC{hand} keys.labels.perturbations{perturbation} '_' keys.tt.tasktypes{tasktype}];
    column_index=find_column_index(tuning_per_unit_table,column_title);
    if ~isempty(column_index)
        cell_idx_tuning_table=cell_idx_tuning_table & cell2mat(tuning_per_unit_table(2:end,column_index));
    else
        disp('combination of tasktype,hands,and choices not existing')
    end
end
tuning_per_unit_table=tuning_per_unit_table([true;cell_idx_tuning_table],:);


%% selection (pick only specific entries)
if size(tuning_per_unit_table,1)>1
    cell_idx_tuning_table=true(size(tuning_per_unit_table,1)-1,1);
    for sel=1:size(keys.tt.selection,1)
        column_index=find_column_index(tuning_per_unit_table,keys.tt.selection{sel,1});
        if isempty(column_index)
            disp(['column not found: ' keys.tt.selection{sel,1}]);
            continue;
        end
        if ischar(keys.tt.selection{sel,2})
            tuning_per_unit_table(2:end,column_index)=cellfun(@num2str,tuning_per_unit_table(2:end,column_index),'UniformOutput',0); %% crazy cellstring bug
            cell_idx_tuning_table=cell_idx_tuning_table & ~cellfun(@isempty,strfind(tuning_per_unit_table(2:end,column_index),keys.tt.selection{sel,2}));
            
        else %% is number!?
            cell_idx_tuning_table=cell_idx_tuning_table & cell2mat(tuning_per_unit_table(2:end,column_index))==keys.tt.selection{sel,2};
        end
    end
    tuning_per_unit_table=tuning_per_unit_table([true;cell_idx_tuning_table],:);
end

%% unselect (remove only specific entries)
if size(tuning_per_unit_table,1)>1
    cell_idx_tuning_table=true(size(tuning_per_unit_table,1)-1,1);
    for sel=1:size(keys.tt.unselect,1)
        column_index=find_column_index(tuning_per_unit_table,keys.tt.unselect{sel,1});
        if isempty(column_index)
            disp(['column not found: ' keys.tt.unselect{sel,1}]);
            continue;
        end
        if ischar(keys.tt.unselect{sel,2})
            cell_idx_tuning_table=cell_idx_tuning_table & cellfun(@(x) isempty(strfind(x,keys.tt.unselect{sel,2})),tuning_per_unit_table(2:end,column_index));
        else %% is number!?
            cell_idx_tuning_table=cell_idx_tuning_table & cell2mat(tuning_per_unit_table(2:end,column_index))~=keys.tt.unselect{sel,2};
        end
    end
    tuning_per_unit_table=tuning_per_unit_table([true;cell_idx_tuning_table],:);
end

%% Selection title (making sure files are named differently)
if ~isempty(keys.tt.selection)
    Sel_for_title=[keys.tt.selection(:,1) repmat({' = '},size(keys.tt.selection(:,1)))...
        cellfun(@num2str,keys.tt.selection(:,2),'uniformoutput',false)  repmat({', '},size(keys.tt.selection(:,1)))]';
else
    Sel_for_title={'';'';'';''};
end;
if ~isempty(keys.tt.unselect)
    Sel_for_title=[Sel_for_title, [keys.tt.unselect(:,1) repmat({' ~= '},size(keys.tt.unselect(:,1)))...
        cellfun(@num2str,keys.tt.unselect(:,2),'uniformoutput',false)  repmat({', '},size(keys.tt.unselect(:,1)))]'];
end;

% if keys.FR_subtract_baseline
%     Sel_for_title =[Sel_for_title,{'base';'=';'INI';', '}];
% end


% if ~isempty(keys.tt.selection)
%     Sel_for_title=[keys.tt.selection(:,1) repmat({' = '},size(keys.tt.selection(:,1)))...
%         cellfun(@num2str,keys.tt.selection(:,2),'uniformoutput',false)  repmat({', '},size(keys.tt.selection(:,1)))]';
% else
%     Sel_for_title={'';'';'';''};
% end;

% selected=cellfun(@num2str,keys.tt.selection(:,2),'Uniformoutput',false);
% if ~isempty(keys.tt.selection)
%     Selection=strcat(keys.tt.selection(:,1),'= ', selected, ',')';
% else
%     Selection={};
% end


end