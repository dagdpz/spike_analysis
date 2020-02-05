function [tuning_per_unit_table Sel_for_title]=ph_reduce_tuning_table(tuning_per_unit_table,keys)

%% monkey
idx_unitID=find_column_index(tuning_per_unit_table,'unit_ID');
if numel(keys.monkey)>3
    tuning_per_unit_table=tuning_per_unit_table([true; ~cellfun(@isempty,strfind(tuning_per_unit_table(2:end,idx_unitID),keys.monkey(1:3)))],:);
end

%% selection (pick only specific entries)
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

%% unselect (remove only specific entries)
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


%% overlapping task types
if keys.tt.only_overlapping_tasktypes
    cell_idx_tuning_table=true(size(tuning_per_unit_table,1)-1,1);
    for d=1:numel(keys.tt.tasktypes) %% crashes for not existing monkey
        if isfield(keys,'instructed_choice') && strcmp(keys.instructed_choice,'ch')
        column_index=find_column_index(tuning_per_unit_table,['existing_ch_' keys.tt.tasktypes{d}]);
        else
        column_index=find_column_index(tuning_per_unit_table,['existing_' keys.tt.tasktypes{d}]);
            
        end
        cell_idx_tuning_table=cell_idx_tuning_table & cell2mat(tuning_per_unit_table(2:end,column_index));
    end
    tuning_per_unit_table=tuning_per_unit_table([true;cell_idx_tuning_table],:);
end

%% only two hands
if keys.tt.only_both_hands
    cell_idx_tuning_table=true(size(tuning_per_unit_table,1)-1,1);
    for d=1:numel(keys.tt.tasktypes)
        column_index=find_column_index(tuning_per_unit_table,[keys.instructed_choice '_hands_main_' keys.tt.tasktypes{d}]);
        tuning_per_unit_table(2:end,column_index)=cellfun(@num2str,tuning_per_unit_table(2:end,column_index),'UniformOutput',0);
        cell_idx_tuning_table=cell_idx_tuning_table & ~strcmp(tuning_per_unit_table(2:end,column_index),'');
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