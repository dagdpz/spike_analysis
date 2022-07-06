function [TT Sel_for_title]=ph_reduce_tuning_table(TT,keys)
%% monkey
idx_unitID=DAG_find_column_index(TT,'unit_ID');
if numel(keys.monkey)>3
    TT=TT([true; ~cellfun(@isempty,strfind(TT(2:end,idx_unitID),keys.monkey(1:3)))],:);
end

%% new stuff to make everything consistent within each other
for c=1:numel(keys.condition_parameters)
    %% hands --> reach_hand; choices --> choice; perturbations --> perturbation;
    CM_cell{c}=keys.tt.(keys.condition_parameters{c});
end
CM=combvec(CM_cell{:})';
for t=1:numel( keys.tt.tasktypes)
    tasktype=keys.tt.tasktypes{t};
    if ~strfind(tasktype,'_')
        disp('keys.tt.tasktypes needs to contain arrangement as well to work properly, f.e.: Ddre_han')
    end
    for r=1:size(CM,1)
        label='';
        for c=1:size(CM,2)
            add_to_label_index=0;
            switch keys.condition_parameters{c}
                case {'choice','reach_hand','perturbation','success','difficulty', 'stimuli_in_2hemifields'}
                    add_to_label_index=1;
            end
            label_index=CM(r,c)+add_to_label_index;
            if ~isnan(label_index)
                to_add=keys.labels.(keys.condition_parameters{c}){label_index};
                if ~isempty(to_add)
                    label=[label '_' to_add];
                end
            else
                label=[label '*']; %% this is essentially the new part for allowing to not consider a certain parameter (?)
            end
        end
        if strcmp(label(1),'*') % not nice, but okay for now. if first is a star, dont remove that one
            labels{r+(t-1)*size(CM,1)}=[label '_' tasktype];
        else
            labels{r+(t-1)*size(CM,1)}=[label(2:end) '_' tasktype];
        end
    end
end

%labels=unique(labels); %? if you do this, conditions dont match any more!
%we should be able to simply remove condition (?) BUT we want either or in * cases

%% trial criterion exclusion looking for labels !
row_index=true(size(TT,1)-1,1);
TT_titles=TT(1,:);
for l=1:numel(labels) %% crashes for not existing monkey
    column_title=['existing_' labels{l}];
    starpos=[0 strfind(column_title,'*') numel(column_title)+1];
    column_index=true(size(TT_titles));
    for s=1:numel(starpos)-1 % this is the loop looking for name parts... (NOT IDEAL if conditions are named similarly)
        column_index=column_index & cellfun(@(x) any(strfind(x,column_title(starpos(s)+1:starpos(s+1)-1))),TT_titles);
    end
    if sum(column_index)>0
        row_index=row_index & any(cell2mat(TT(2:end,column_index)),2);
    else
        disp([keys.condition_parameters{:} labels{l} 'not existing'])
    end
end
TT=TT([true;row_index],:);
TT=ph_target_reassign(keys,TT);

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