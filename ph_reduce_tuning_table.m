function [TT Sel_for_title]=ph_reduce_tuning_table(TT,keys)
%% monkey
idx_unitID=DAG_find_column_index(TT,'unit_ID');
if numel(keys.monkey)>3
    TT=TT([true; ~cellfun(@isempty,strfind(TT(2:end,idx_unitID),keys.monkey(1:3)))],:);
end

%% overlapping tasktypes_hands_choices - trial criterion exclusion
tasktypes_hands_choices=combvec((1:numel(keys.tt.tasktypes)),keys.tt.hands,keys.tt.choices,keys.tt.perturbations);
row_index=true(size(TT,1)-1,1);
for t=1:size(tasktypes_hands_choices,2) %% crashes for not existing monkey
    tasktype=tasktypes_hands_choices(1,t);
    hand    =tasktypes_hands_choices(2,t)+1;
    choice  =tasktypes_hands_choices(3,t)+1;
    perturbation  =tasktypes_hands_choices(4,t)+1;
    column_title=['existing_' keys.labels.choices{choice} '_' keys.labels.handsIC{hand} keys.labels.perturbations{perturbation} '_' keys.tt.tasktypes{tasktype}];
    column_index=DAG_find_column_index(TT,column_title);
    if ~isempty(column_index)
        row_index=row_index & cell2mat(TT(2:end,column_index));
    else
        disp(['combination of tasktype ' keys.tt.tasktypes{tasktype} ',hands ' keys.labels.handsIC{hand} 'and choices ' keys.labels.choices{choice} 'not existing'])
        if ~strfind(keys.tt.tasktypes{tasktype},'_')
            disp('tasktype needs to contain arrangement as well to work properly, f.e.: Ddre_han')
        end
    end
end
TT=TT([true;row_index],:);

%% selection (pick only specific entries)
if size(TT,1)>1
    row_index=true(size(TT,1)-1,1);
    for sel=1:size(keys.tt.selection,1)
        column_index=DAG_find_column_index(TT,keys.tt.selection{sel,1});
        if isempty(column_index)
            disp(['column not found: ' keys.tt.selection{sel,1}]);
            continue;
        end
        if ischar(keys.tt.selection{sel,2})
            TT(2:end,column_index)=cellfun(@num2str,TT(2:end,column_index),'UniformOutput',0); %% crazy cellstring bug
            row_index=row_index & ~cellfun(@isempty,strfind(TT(2:end,column_index),keys.tt.selection{sel,2}));
            
        else %% is number!?
            row_index=row_index & cell2mat(TT(2:end,column_index))==keys.tt.selection{sel,2};
        end
    end
    TT=TT([true;row_index],:);
end

%% unselect (remove only specific entries)
if size(TT,1)>1
    row_index=true(size(TT,1)-1,1);
    for sel=1:size(keys.tt.unselect,1)
        column_index=DAG_find_column_index(TT,keys.tt.unselect{sel,1});
        if isempty(column_index)
            disp(['column not found: ' keys.tt.unselect{sel,1}]);
            continue;
        end
        if ischar(keys.tt.unselect{sel,2})
            row_index=row_index & cellfun(@(x) isempty(strfind(x,keys.tt.unselect{sel,2})),TT(2:end,column_index));
        else %% is number!?
            row_index=row_index & cell2mat(TT(2:end,column_index))~=keys.tt.unselect{sel,2};
        end
    end
    TT=TT([true;row_index],:);
    if ~isempty(keys.tt.unselected_list)
        load([keys.tt.unselected_list{:} filesep 'list']);
        TT=TT([true; ~ismember(TT(2:end,idx_unitID),example_list)],:);
    end
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

end