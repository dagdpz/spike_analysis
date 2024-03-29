function TT=ph_load_tuning_table(keys)
%% Reading in table
load(keys.anova_table_file);
TT=tuning_per_unit_table;

% fix blanks in column titles
for k=1:size(TT,2)
    blanks=strfind(TT{1,k},' ');
    if any(blanks)
        TT{1,k}(blanks) ='_';
    end
end

% add column for no grouping
TT(:,end+1)=num2cell(repmat('Y',size(TT,1),1));
TT(1,end)={'ungrouped'};

%% add trial criterion column (exisiting_tasktype)
% this part is linked to several definitions
% First of all, regardless of conditions, there is always going to be a
% separation by instructed and choice trials
% the minimum number of trials as well as the spatial conditions to have
% this minimum number of trials are independently definable for in/ch
% to check for spatial conditions available or to add new ones, see
% ph_get_minimum_trials
for criterium={'in','ch'}
    crit=criterium{:};
    switch keys.tt.(['trial_criterion_' crit])
        case 'per_hemifield_and_perturbation'
            strtofind='trials_per_hemifield'; %% ??? - there will be a problem with perturbation stuff !!
            disp(['per_hemifield_and_perturbation not supported as trial criterion' ]);
    end
    
    strtofind=['trials_' keys.tt.(['trial_criterion_' crit])];
    taskcaseexistingindex=find(~cellfun(@isempty,strfind(TT(1,:),strtofind)));
    for I=taskcaseexistingindex
        if ~(any(strfind(TT{1,I},[crit '_'])==1) || any(strfind(TT{1,I},['_' crit '_'])))
            continue;
        end
        task_existing_column=num2cell(cellfun(@(x) ~isempty(x)&&~ischar(x)&&x>=keys.cal.(['min_trials_' crit]),TT(:,I)));
        strpos=strfind(TT{1,I},strtofind);
        task_existing_column{1,1}=['existing_' TT{1,I}([1:strpos-1,strpos+numel(strtofind)+1:end])];
        TT(:,end+1)=task_existing_column;
    end
end
end