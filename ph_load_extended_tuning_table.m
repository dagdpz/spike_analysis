function tuning_per_unit_table=ph_load_extended_tuning_table(keys)

%% Reading in table

tuning_per_unit_table={};
load(keys.anova_table_file);

tuning_per_unit_table(:,end+1)=num2cell(repmat('Y',size(tuning_per_unit_table,1),1));
tuning_per_unit_table(1,end)={'ungrouped'};
taskcaseexistingindex=find(~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),'_trials_per_condition')));
taskcaseexistingindex_hf=find(~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),'_trials_per_hemifield')));
taskcaseexistingindex_chf=find(~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),'_trials_per_congruent_hand_hemifield')));
taskcaseexistingindex_tot=find(~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),'_trials_total')));
if ~isempty(taskcaseexistingindex)
    taskcases={};
    %% Trial criterion per position
    for I=taskcaseexistingindex
        if ~((any(strfind(tuning_per_unit_table{1,I},'in_')==1) && strcmp(keys.tt.trial_criterion_in,'per_position')) ||...
                (any(strfind(tuning_per_unit_table{1,I},'ch_')==1) && strcmp(keys.tt.trial_criterion_ch,'per_position')) )
            continue;
        end
        task_existing_column=num2cell(cellfun(@(x) ~isempty(x)&&~ischar(x)&&x>=keys.cal.min_trials_per_condition,tuning_per_unit_table(:,I)));
        strpos=strfind(tuning_per_unit_table{1,I},'_trials_per_condition');
        task_existing_column{1,1}=['existing_' tuning_per_unit_table{1,I}([1:strpos,strpos+22:end])];
        tuning_per_unit_table(:,end+1)=task_existing_column;
        taskcases = [taskcases, {tuning_per_unit_table{1,I}(end-7:end)}];
    end
    %% trial criterion per hemifield (especially useful for choices)
    for I=taskcaseexistingindex_hf
        if ~((any(strfind(tuning_per_unit_table{1,I},'in_')==1) && strcmp(keys.tt.trial_criterion_in,'per_hemifield')) ||...
                (any(strfind(tuning_per_unit_table{1,I},'ch_')==1) && strcmp(keys.tt.trial_criterion_ch,'per_hemifield')))
            continue;
        end
        task_existing_column=num2cell(cellfun(@(x) ~isempty(x)&&~ischar(x)&&x>=keys.cal.min_trials_per_condition,tuning_per_unit_table(:,I)));
        strpos=strfind(tuning_per_unit_table{1,I},'_trials_per_hemifield');
        task_existing_column{1,1}=['existing_' tuning_per_unit_table{1,I}([1:strpos,strpos+22:end])];
        tuning_per_unit_table(:,end+1)=task_existing_column;
        taskcases = [taskcases, {tuning_per_unit_table{1,I}(end-7:end)}];
    end
    %%  To analyse only congruent hand\space condition
    for I=taskcaseexistingindex_chf
        if ~((any(strfind(tuning_per_unit_table{1,I},'in_')==1) && strcmp(keys.tt.trial_criterion_in,'per_congruent_hand_hemifield')) ||...
                (any(strfind(tuning_per_unit_table{1,I},'ch_')==1) && strcmp(keys.tt.trial_criterion_ch,'per_congruent_hand_hemifield')))
            continue;
        end
        task_existing_column=num2cell(cellfun(@(x) ~isempty(x)&&~ischar(x)&&x>=keys.cal.min_trials_per_condition,tuning_per_unit_table(:,I)));
        strpos=strfind(tuning_per_unit_table{1,I},'_trials_per_congruent_hand_hemifield');
        task_existing_column{1,1}=['existing_' tuning_per_unit_table{1,I}([1:strpos,strpos+37:end])];
        tuning_per_unit_table(:,end+1)=task_existing_column;
        taskcases = [taskcases, {tuning_per_unit_table{1,I}(end-7:end)}];
    end
    %%  General minimum number of trials criterion (across all hand/hemifield/position conditions)
    for I=taskcaseexistingindex_tot
        if ~((any(strfind(tuning_per_unit_table{1,I},'in_')==1) && strcmp(keys.tt.trial_criterion_in,'total')) ||...
                (any(strfind(tuning_per_unit_table{1,I},'ch_')==1) && strcmp(keys.tt.trial_criterion_ch,'total')))
            continue;
        end
        task_existing_column=num2cell(cellfun(@(x) ~isempty(x)&&~ischar(x)&&x>=keys.cal.min_trials_per_condition,tuning_per_unit_table(:,I)));
        strpos=strfind(tuning_per_unit_table{1,I},'_trials_total');
        task_existing_column{1,1}=['existing_' tuning_per_unit_table{1,I}([1:strpos,strpos+14:end])];
        tuning_per_unit_table(:,end+1)=task_existing_column;
        taskcases = [taskcases, {tuning_per_unit_table{1,I}(end-7:end)}];
    end
else % this is for previous version (backwards compatibility) ... not needed any more??
    taskcaseexistingindex=find(~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),'in_epoch_main'))); %not ideal
    taskcases={};
    for I=taskcaseexistingindex
        task_existing_column=num2cell(~cellfun(@isempty,tuning_per_unit_table(:,I)));
        task_existing_column{1,1}=['existing' tuning_per_unit_table{1,I}(14:end)];
        tuning_per_unit_table(:,end+1)=task_existing_column;
        taskcases = [taskcases, {tuning_per_unit_table{1,I}(15:end)}];
    end
    %% to be combined with previous
    taskcaseexistingindex=find(~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),'ch_epoch_main'))); %not ideal
    %taskcases={};
    for I=taskcaseexistingindex
        task_existing_column=num2cell(~cellfun(@isempty,tuning_per_unit_table(:,I)));
        task_existing_column{1,1}=['existing' tuning_per_unit_table{1,I}(1:3) tuning_per_unit_table{1,I}(14:end)];
        tuning_per_unit_table(:,end+1)=task_existing_column;
        %taskcases = [taskcases, {tuning_per_unit_table{1,I}(15:end)}];
    end
    
end


combined_column=cell(size(tuning_per_unit_table,1),1);
for props=2:size(keys.tt.combine_tuning_properties,2)
    column_index=DAG_find_column_index(tuning_per_unit_table,keys.tt.combine_tuning_properties{1,props});
    if isempty(column_index)
        disp(['column not found: ' keys.tt.combine_tuning_properties{1,props}]);
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
taskcases=unique(taskcases);
epochs={'Fhol','Cue','Cue2','EDel','Del','MemE','MemL','PreS','PreR','PeriS','Pre2','Peri2','PreG','CueG','TIhol','THol'};
for t=1:numel(taskcases)
    taskcase=taskcases{t};
    idx.space   =DAG_find_column_index(tuning_per_unit_table,[keys.tt.IC_for_criterion '_spaceLR_main_'   taskcase] );
    idx.epoch   =DAG_find_column_index(tuning_per_unit_table,[keys.tt.IC_for_criterion '_epoch_main_'     taskcase] );
    idx.hands   =DAG_find_column_index(tuning_per_unit_table,[keys.tt.IC_for_criterion '_hands_main_'     taskcase] );
    idx.SxH     =DAG_find_column_index(tuning_per_unit_table,[keys.tt.IC_for_criterion '_SxH_'            taskcase] );
    idx.ExS     =DAG_find_column_index(tuning_per_unit_table,[keys.tt.IC_for_criterion '_ExS_'            taskcase] );
    idx.ExH     =DAG_find_column_index(tuning_per_unit_table,[keys.tt.IC_for_criterion '_ExH_'            taskcase] );
    idx.ESH     =DAG_find_column_index(tuning_per_unit_table,[keys.tt.IC_for_criterion '_ExSxH_'          taskcase] );
    
    
    for e=1:numel(epochs)
        idx.(epochs{e})                     =DAG_find_column_index(tuning_per_unit_table,[keys.tt.IC_for_criterion '_AH_' epochs{e} '_epoch_'   taskcase] );
        idx.([epochs{e} '_space_DF'])       =DAG_find_column_index(tuning_per_unit_table,[keys.tt.IC_for_criterion epochs{e} '_spaceLR_DF_'   taskcase] );
        idx.([epochs{e} '_IS_FR'])          =DAG_find_column_index(tuning_per_unit_table,[keys.tt.IC_for_criterion '_AH_IS_' epochs{e} '_epoch_FR_'  taskcase] );
        idx.([epochs{e} '_CS_FR'])          =DAG_find_column_index(tuning_per_unit_table,[keys.tt.IC_for_criterion '_AH_CS_' epochs{e} '_epoch_FR_'  taskcase] );
        idx.([epochs{e} '_IS_EN'])          =DAG_find_column_index(tuning_per_unit_table,[keys.tt.IC_for_criterion '_AH_IS_' epochs{e} '_epoch_DF_'  taskcase] );
        idx.([epochs{e} '_CS_EN'])          =DAG_find_column_index(tuning_per_unit_table,[keys.tt.IC_for_criterion '_AH_CS_' epochs{e} '_epoch_DF_'  taskcase] );
        idx.([epochs{e} '_in_space_DF'])    =DAG_find_column_index(tuning_per_unit_table,['in_' epochs{e} '_spaceLR_DF_'   taskcase] );
        idx.([epochs{e} '_ch_space_DF'])    =DAG_find_column_index(tuning_per_unit_table,['ch_' epochs{e} '_spaceLR_DF_'   taskcase] );
        idx.([epochs{e} '_in_IS_FR'])       =DAG_find_column_index(tuning_per_unit_table,['in_AH_IS_' epochs{e} '_epoch_FR_'  taskcase] );
        idx.([epochs{e} '_in_CS_FR'])       =DAG_find_column_index(tuning_per_unit_table,['in_AH_CS_' epochs{e} '_epoch_FR_'  taskcase] );
        idx.([epochs{e} '_ch_IS_FR'])       =DAG_find_column_index(tuning_per_unit_table,['ch_AH_IS_' epochs{e} '_epoch_FR_'  taskcase] );
        idx.([epochs{e} '_ch_CS_FR'])       =DAG_find_column_index(tuning_per_unit_table,['ch_AH_CS_' epochs{e} '_epoch_FR_'  taskcase] );
        
        
        idx.([epochs{e} '_in_IH_space'])       =DAG_find_column_index(tuning_per_unit_table,['in_IH_' epochs{e} '_spaceLR_'  taskcase] );
        idx.([epochs{e} '_in_CH_space'])       =DAG_find_column_index(tuning_per_unit_table,['in_CH_' epochs{e} '_spaceLR_'  taskcase] );
        idx.([epochs{e} '_in_IS_hands'])       =DAG_find_column_index(tuning_per_unit_table,['in_IS_' epochs{e} '_hands_'  taskcase] );
        idx.([epochs{e} '_in_CS_hands'])       =DAG_find_column_index(tuning_per_unit_table,['in_CS_' epochs{e} '_hands_'  taskcase] );
        
    end
    
    
    if ~isempty(idx.epoch); tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.epoch))],idx.epoch)  ={NaN}; end;
    if ~isempty(idx.space); tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.space))],idx.space)  ={''}; end;
    if ~isempty(idx.hands); tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.hands))],idx.hands)  ={''}; end;
    if ~isempty(idx.SxH);   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.SxH))],idx.SxH)      ={''}; end;
    if ~isempty(idx.ExS);   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.ExS))],idx.ExS)      ={NaN}; end;
    if ~isempty(idx.ExH);   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.ExH))],idx.ExH)      ={NaN}; end;
    if ~isempty(idx.ESH);   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.ESH))],idx.ESH)      ={NaN}; end;
    
    
    for e=1:numel(epochs)
        
        if ~isempty(idx.(epochs{e}));   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.(epochs{e})))],idx.(epochs{e}))      ={''}; end;
        if ~isempty(idx.([epochs{e} '_space_DF']));   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.([epochs{e} '_space_DF'])))],idx.([epochs{e} '_space_DF']))={single(NaN)}; end;
        if ~isempty(idx.([epochs{e} '_IS_FR']));   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.([epochs{e} '_IS_FR'])))],idx.([epochs{e} '_IS_FR']))={single(NaN)}; end;
        if ~isempty(idx.([epochs{e} '_CS_FR']));   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.([epochs{e} '_CS_FR'])))],idx.([epochs{e} '_CS_FR']))={single(NaN)}; end;
        if ~isempty(idx.([epochs{e} '_IS_EN']));   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.([epochs{e} '_IS_EN'])))],idx.([epochs{e} '_IS_EN']))={single(NaN)}; end;
        if ~isempty(idx.([epochs{e} '_CS_EN']));   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.([epochs{e} '_CS_EN'])))],idx.([epochs{e} '_CS_EN']))={single(NaN)}; end;
        if ~isempty(idx.([epochs{e} '_in_space_DF']));   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.([epochs{e} '_in_space_DF'])))],idx.([epochs{e} '_in_space_DF']))={single(NaN)}; end;
        if ~isempty(idx.([epochs{e} '_ch_space_DF']));   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.([epochs{e} '_ch_space_DF'])))],idx.([epochs{e} '_ch_space_DF']))={single(NaN)}; end;
        if ~isempty(idx.([epochs{e} '_in_IS_FR']));   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.([epochs{e} '_in_IS_FR'])))],idx.([epochs{e} '_in_IS_FR']))={single(NaN)}; end;
        if ~isempty(idx.([epochs{e} '_in_CS_FR']));   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.([epochs{e} '_in_CS_FR'])))],idx.([epochs{e} '_in_CS_FR']))={single(NaN)}; end;
        if ~isempty(idx.([epochs{e} '_ch_IS_FR']));   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.([epochs{e} '_ch_IS_FR'])))],idx.([epochs{e} '_ch_IS_FR']))={single(NaN)}; end;
        if ~isempty(idx.([epochs{e} '_ch_CS_FR']));   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.([epochs{e} '_ch_CS_FR'])))],idx.([epochs{e} '_ch_CS_FR']))={single(NaN)}; end;
        
        
        if ~isempty(idx.([epochs{e} '_in_IH_space']));   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.([epochs{e} '_in_IH_space'])))],idx.([epochs{e} '_in_IH_space']))      ={'-'}; end;
        if ~isempty(idx.([epochs{e} '_in_CH_space']));   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.([epochs{e} '_in_CH_space'])))],idx.([epochs{e} '_in_CH_space']))      ={'-'}; end;
        if ~isempty(idx.([epochs{e} '_in_IS_hands']));   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.([epochs{e} '_in_IS_hands'])))],idx.([epochs{e} '_in_IS_hands']))      ={'-'}; end;
        if ~isempty(idx.([epochs{e} '_in_CS_hands']));   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.([epochs{e} '_in_CS_hands'])))],idx.([epochs{e} '_in_CS_hands']))      ={'-'}; end;
    end
    
    %
    % if ~isempty(idx.Peri2_IS_EN);   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.Peri2_IS_EN))],idx.Peri2_IS_EN)  ={single(NaN)}; end;
    % if ~isempty(idx.Peri2_CS_EN);   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.Peri2_CS_EN))],idx.Peri2_CS_EN)  ={single(NaN)}; end;
    % if ~isempty(idx.PeriS_IS_EN);   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.PeriS_IS_EN))],idx.PeriS_IS_EN)  ={single(NaN)}; end;
    % if ~isempty(idx.PeriS_CS_EN);   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.PeriS_CS_EN))],idx.PeriS_CS_EN)  ={single(NaN)}; end;
    % if ~isempty(idx.TIhol_IS_EN);   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.TIhol_IS_EN))],idx.TIhol_IS_EN)  ={single(NaN)}; end;
    % if ~isempty(idx.TIhol_CS_EN);   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.TIhol_CS_EN))],idx.TIhol_CS_EN)  ={single(NaN)}; end;
    % if ~isempty(idx.Cue_IS_EN);     tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.Cue_IS_EN))],idx.Cue_IS_EN)      ={single(NaN)}; end;
    % if ~isempty(idx.Cue_CS_EN);     tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.Cue_CS_EN))],idx.Cue_CS_EN)      ={single(NaN)}; end;
    % if ~isempty(idx.CueG_IS_EN);     tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.CueG_IS_EN))],idx.CueG_IS_EN)      ={single(NaN)}; end;
    % if ~isempty(idx.CueG_CS_EN);     tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.CueG_CS_EN))],idx.CueG_CS_EN)      ={single(NaN)}; end;
    % if ~isempty(idx.PreG_IS_EN);     tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.PreG_IS_EN))],idx.PreG_IS_EN)      ={single(NaN)}; end;
    % if ~isempty(idx.PreG_CS_EN);     tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.PreG_CS_EN))],idx.PreG_CS_EN)      ={single(NaN)}; end;
    % if ~isempty(idx.CueG_IS_FR);     tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.CueG_IS_FR))],idx.CueG_IS_FR)      ={single(NaN)}; end;
    % if ~isempty(idx.CueG_CS_FR);     tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.CueG_CS_FR))],idx.CueG_CS_FR)      ={single(NaN)}; end;
    % if ~isempty(idx.PreG_IS_FR);     tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.PreG_IS_FR))],idx.PreG_IS_FR)      ={single(NaN)}; end;
    % if ~isempty(idx.PreG_CS_FR);     tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.PreG_CS_FR))],idx.PreG_CS_FR)      ={single(NaN)}; end;
    % if ~isempty(idx.Peri2_IS_FR);   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.Peri2_IS_FR))],idx.Peri2_IS_FR)  ={single(NaN)}; end;
    % if ~isempty(idx.Peri2_CS_FR);   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.Peri2_CS_FR))],idx.Peri2_CS_FR)  ={single(NaN)}; end;
    % if ~isempty(idx.PeriS_IS_FR);   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.PeriS_IS_FR))],idx.PeriS_IS_FR)  ={single(NaN)}; end;
    % if ~isempty(idx.PeriS_CS_FR);   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.PeriS_CS_FR))],idx.PeriS_CS_FR)  ={single(NaN)}; end;
    % if ~isempty(idx.TIhol_IS_FR);   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.TIhol_IS_FR))],idx.TIhol_IS_FR)  ={single(NaN)}; end;
    % if ~isempty(idx.TIhol_CS_FR);   tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.TIhol_CS_FR))],idx.TIhol_CS_FR)  ={single(NaN)}; end;
    % if ~isempty(idx.Cue_IS_FR);     tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.Cue_IS_FR))],idx.Cue_IS_FR)      ={single(NaN)}; end;
    % if ~isempty(idx.Cue_CS_FR);     tuning_per_unit_table([false; cellfun(@(x) isempty(x) ,tuning_per_unit_table(2:end,idx.Cue_CS_FR))],idx.Cue_CS_FR)      ={single(NaN)}; end;
    
    switch keys.tt.space_criterion %% consider only neurons that show 1) SxE or space main effect; 2) SxE; 3) take all (no criterion)?
        case 'interaction or space only'
            space_or_interaction=([tuning_per_unit_table{2:end,idx.ExS}]==1)' | ~ismember(tuning_per_unit_table(2:end,idx.space),'-');
        case 'interaction only'
            space_or_interaction=([tuning_per_unit_table{2:end,idx.ExS}]==1)';
        case 'none'
            space_or_interaction=true(size(tuning_per_unit_table,1)-1,1);
    end
    %%% goofy
    switch keys.tt.epoch_criterion %% consider only neurons that show 1) SxE or space main effect; 2) SxE; 3) take all (no criterion)?
        case 'SxE or epoch only'
            epoch_or_interaction=([tuning_per_unit_table{2:end,idx.ExS}]==1)' | ([tuning_per_unit_table{2:end,idx.epoch}]==1)';
        case 'HxE or epoch only'
            epoch_or_interaction=([tuning_per_unit_table{2:end,idx.ExH}]==1)' | ([tuning_per_unit_table{2:end,idx.epoch}]==1)';
        case 'SHE or epoch only'
            epoch_or_interaction=([tuning_per_unit_table{2:end,idx.ESH}]==1)' | ([tuning_per_unit_table{2:end,idx.epoch}]==1)';
        case 'SxE only'
            epoch_or_interaction=([tuning_per_unit_table{2:end,idx.ExS}]==1)';
        case 'HxE only'
            epoch_or_interaction=([tuning_per_unit_table{2:end,idx.ExH}]==1)';
        case 'SHE only'
            epoch_or_interaction=([tuning_per_unit_table{2:end,idx.ESH}]==1)';
        case 'none'
            epoch_or_interaction=true(size(tuning_per_unit_table,1)-1,1);
    end
    switch keys.tt.hands_criterion %% consider only neurons that show 1) HxE or hand main effect; 2) HxE; 3) take all (no criterion)?
        case 'interaction or hands only'
            hands_or_interaction=([tuning_per_unit_table{2:end,idx.ExH}]==1)' | ~ismember(tuning_per_unit_table(2:end,idx.hands),'-');
        case 'interaction only'
            hands_or_interaction=([tuning_per_unit_table{2:end,idx.ExH}]==1)';
        case 'none'
            hands_or_interaction=true(size(tuning_per_unit_table,1)-1,1);
    end
    switch keys.tt.SXH_criterion %% consider only neurons that show 1) SxHxE or SxH; 2) SxHxE; 3) take all (no criterion)?
        case 'SXH SXE HXE SHE'
            SXH_or_interaction=([tuning_per_unit_table{2:end,idx.ESH}]==1)' | ~ismember(tuning_per_unit_table(2:end,idx.SxH),'-')...
                |([tuning_per_unit_table{2:end,idx.ExH}]==1)' | ([tuning_per_unit_table{2:end,idx.ExS}]==1)';
        case 'interaction or SXH only'
            SXH_or_interaction=([tuning_per_unit_table{2:end,idx.ESH}]==1)' | ~ismember(tuning_per_unit_table(2:end,idx.SxH),'-');
        case 'interaction only'
            SXH_or_interaction=([tuning_per_unit_table{2:end,idx.ESH}]==1)';
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
    
    %% adding columns for visual, motor, visuomotor, and fixation cells
    if any(idx.Cue) && any(idx.TIhol)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['visual_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell(ismember(tuning_per_unit_table(2:end,idx.Cue),{'en','su','bi'}) & ~ismember(tuning_per_unit_table(2:end,idx.TIhol),{'en','su','bi'}));
    end
    
    if any(idx.Cue) && any(idx.TIhol)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['motor_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell(~ismember(tuning_per_unit_table(2:end,idx.Cue),{'en','su','bi'}) & ismember(tuning_per_unit_table(2:end,idx.TIhol),{'en','su','bi'}));
    end
    
    if any(idx.Cue) && any(idx.TIhol)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['visuomotor_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell(ismember(tuning_per_unit_table(2:end,idx.Cue),{'en','su','bi'}) & ismember(tuning_per_unit_table(2:end,idx.TIhol),{'en','su','bi'}));
    end
    
    if any(idx.Cue) && any(idx.PreS)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['visual_pre_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell(ismember(tuning_per_unit_table(2:end,idx.Cue),{'en','su','bi'}) & ~ismember(tuning_per_unit_table(2:end,idx.PreS),{'en','su','bi'}));
    end
    
    if any(idx.Cue) && any(idx.PreS)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['motor_pre_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell(~ismember(tuning_per_unit_table(2:end,idx.Cue),{'en','su','bi'}) & ismember(tuning_per_unit_table(2:end,idx.PreS),{'en','su','bi'}));
    end
    
    if any(idx.Cue) && any(idx.PreS)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['visuomotor_pre_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell(ismember(tuning_per_unit_table(2:end,idx.Cue),{'en','su','bi'}) & ismember(tuning_per_unit_table(2:end,idx.PreS),{'en','su','bi'}));
    end
    
    
    if any(idx.Cue) && any(idx.Pre2)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['visual_pre2_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell(ismember(tuning_per_unit_table(2:end,idx.Cue),{'en','su','bi'}) & ~ismember(tuning_per_unit_table(2:end,idx.Pre2),{'en','su','bi'}));
    end
    
    if any(idx.Cue) && any(idx.Pre2)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['motor_pre2_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell(~ismember(tuning_per_unit_table(2:end,idx.Cue),{'en','su','bi'}) & ismember(tuning_per_unit_table(2:end,idx.Pre2),{'en','su','bi'}));
    end
    
    if any(idx.Cue) && any(idx.Pre2)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['visuomotor_pre2_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell(ismember(tuning_per_unit_table(2:end,idx.Cue),{'en','su','bi'}) & ismember(tuning_per_unit_table(2:end,idx.Pre2),{'en','su','bi'}));
    end
    
    
    if any(idx.Cue) && any(idx.TIhol)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['notclassified_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell(~ismember(tuning_per_unit_table(2:end,idx.Cue),{'en','su','bi'}) & ~ismember(tuning_per_unit_table(2:end,idx.TIhol),{'en','su','bi'}));
    end
    %
    %     if any(idx.mem)
    %     n_column=n_column+1;
    %     tuning_per_unit_table{1,n_column}=['delay_' taskcase];
    %     tuning_per_unit_table(2:end,n_column)=num2cell(ismember(tuning_per_unit_table(2:end,idx.mem),'en') );
    %     end
    %
    if any(idx.Fhol) && any(idx.PeriS) && any(idx.Cue) && any(idx.TIhol)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['fixation_only_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell(ismember(tuning_per_unit_table(2:end,idx.Fhol),'en') & ismember(tuning_per_unit_table(2:end,idx.PeriS),'-') & ~ismember(tuning_per_unit_table(2:end,idx.Cue),'en') & ~ismember(tuning_per_unit_table(2:end,idx.TIhol),'en'));
    end
    
    if any(idx.Fhol) && any(idx.PeriS) && any(idx.Cue) && any(idx.TIhol)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['fixation_and_sac_suppression_' taskcase];
        %tuning_per_unit_table(2:end,n_column)=num2cell(ismember(tuning_per_unit_table(2:end,idx.Fhol),'en') & ismember(tuning_per_unit_table(2:end,idx.PeriS),'su') & ismember(tuning_per_unit_table(2:end,idx.TIhol),'en'));
        tuning_per_unit_table(2:end,n_column)=num2cell(ismember(tuning_per_unit_table(2:end,idx.Fhol),'en') & ismember(tuning_per_unit_table(2:end,idx.PeriS),'su') & ~ismember(tuning_per_unit_table(2:end,idx.Cue),'en') & ~ismember(tuning_per_unit_table(2:end,idx.TIhol),'en'));
    end
    if any(idx.Fhol) && any(idx.PeriS) && any(idx.Cue) && any(idx.TIhol)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['Sac_supression_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell(ismember(tuning_per_unit_table(2:end,idx.Fhol),'-') & ismember(tuning_per_unit_table(2:end,idx.PeriS),'su') & ~ismember(tuning_per_unit_table(2:end,idx.Cue),'en') & ~ismember(tuning_per_unit_table(2:end,idx.TIhol),'en'));
    end
    
    
    %     if any(idx.Cue) && any(idx.PeriS)
    %     n_column=n_column+1;
    %     tuning_per_unit_table{1,n_column}=['visual_peri_' taskcase];
    %     tuning_per_unit_table(2:end,n_column)=num2cell(ismember(tuning_per_unit_table(2:end,idx.Cue),'en') & ~ismember(tuning_per_unit_table(2:end,idx.PeriS),'en'));
    %      end
    %
    %      if any(idx.Cue) && any(idx.PeriS)
    %     n_column=n_column+1;
    %     tuning_per_unit_table{1,n_column}=['motor_peri_' taskcase];
    %     tuning_per_unit_table(2:end,n_column)=num2cell(~ismember(tuning_per_unit_table(2:end,idx.Cue),'en') & ismember(tuning_per_unit_table(2:end,idx.PeriS),'en'));
    %      end
    %
    %    if any(idx.Cue) && any(idx.PeriS)
    %     n_column=n_column+1;
    %     tuning_per_unit_table{1,n_column}=['visuomotor_peri_' taskcase];
    %     tuning_per_unit_table(2:end,n_column)=num2cell(ismember(tuning_per_unit_table(2:end,idx.Cue),'en') & ismember(tuning_per_unit_table(2:end,idx.PeriS),'en'));
    %    end
    
    if any(idx.Cue) && any(idx.PeriS)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['visual_peri_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell(ismember(tuning_per_unit_table(2:end,idx.Cue),{'en','su','bi'}) & ~ismember(tuning_per_unit_table(2:end,idx.PeriS),{'en','su','bi'}));
    end
    
    if any(idx.Cue) && any(idx.PeriS)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['motor_peri_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell(~ismember(tuning_per_unit_table(2:end,idx.Cue),{'en','su','bi'}) & ismember(tuning_per_unit_table(2:end,idx.PeriS),{'en','su','bi'}));
    end
    
    if any(idx.Cue) && any(idx.PeriS)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['visuomotor_peri_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell(ismember(tuning_per_unit_table(2:end,idx.Cue),{'en','su','bi'}) & ismember(tuning_per_unit_table(2:end,idx.PeriS),{'en','su','bi'}));
    end
    %
    %
    %     if any(idx.mem)
    %     n_column=n_column+1;
    %     tuning_per_unit_table{1,n_column}=['delay_' taskcase];
    %     tuning_per_unit_table(2:end,n_column)=num2cell(ismember(tuning_per_unit_table(2:end,idx.mem),'en') );
    %     end
    %
    %     if any(idx.Fhol) && any(idx.PeriS) && any(idx.Cue)
    %     n_column=n_column+1;
    %     tuning_per_unit_table{1,n_column}=['fixation_only_' taskcase];
    %     tuning_per_unit_table(2:end,n_column)=num2cell(ismember(tuning_per_unit_table(2:end,idx.Fhol),'en') & ismember(tuning_per_unit_table(2:end,idx.PeriS),'-') & ismember(tuning_per_unit_table(2:end,idx.Cue),'-'));
    %     end
    %
    %     if any(idx.Fhol) && any(idx.PeriS)
    %     n_column=n_column+1;
    %     tuning_per_unit_table{1,n_column}=['fixation_and_sac_suppression_' taskcase];
    %     tuning_per_unit_table(2:end,n_column)=num2cell(ismember(tuning_per_unit_table(2:end,idx.Fhol),'en') & ismember(tuning_per_unit_table(2:end,idx.PeriS),'su'));
    %     end
    %     if any(idx.Fhol) && any(idx.PeriS)
    %     n_column=n_column+1;
    %     tuning_per_unit_table{1,n_column}=['Sac_supression_' taskcase];
    %     tuning_per_unit_table(2:end,n_column)=num2cell(ismember(tuning_per_unit_table(2:end,idx.Fhol),'-') & ismember(tuning_per_unit_table(2:end,idx.PeriS),'su'));
    %     end
    
    if any(idx.PeriS_IS_FR) && any(idx.Cue_IS_FR)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['VMI_peri_IS_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell((cell2mat(tuning_per_unit_table(2:end,idx.PeriS_IS_FR)) - cell2mat(tuning_per_unit_table(2:end,idx.Cue_IS_FR)))./...
            (cell2mat(tuning_per_unit_table(2:end,idx.PeriS_IS_FR)) + cell2mat(tuning_per_unit_table(2:end,idx.Cue_IS_FR))));
    end
    
    if any(idx.PeriS_CS_FR) && any(idx.Cue_CS_FR)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['VMI_peri_CS_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell((cell2mat(tuning_per_unit_table(2:end,idx.PeriS_CS_FR)) - cell2mat(tuning_per_unit_table(2:end,idx.Cue_CS_FR)))./...
            (cell2mat(tuning_per_unit_table(2:end,idx.PeriS_CS_FR)) + cell2mat(tuning_per_unit_table(2:end,idx.Cue_CS_FR))));
    end
    
    
    if any(idx.PeriS_IS_EN) && any(idx.Cue_IS_EN)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['VMI_periEN_IS_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell((abs(cell2mat(tuning_per_unit_table(2:end,idx.PeriS_IS_EN))) - abs(cell2mat(tuning_per_unit_table(2:end,idx.Cue_IS_EN))))./...
            (abs(cell2mat(tuning_per_unit_table(2:end,idx.PeriS_IS_EN))) + abs(cell2mat(tuning_per_unit_table(2:end,idx.Cue_IS_EN)))));
    end
    if any(idx.PeriS_CS_EN) && any(idx.Cue_CS_EN)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['VMI_periEN_CS_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell((abs(cell2mat(tuning_per_unit_table(2:end,idx.PeriS_CS_EN))) - abs(cell2mat(tuning_per_unit_table(2:end,idx.Cue_CS_EN))))./...
            (abs(cell2mat(tuning_per_unit_table(2:end,idx.PeriS_CS_EN))) + abs(cell2mat(tuning_per_unit_table(2:end,idx.Cue_CS_EN)))));
    end
    
    if any(idx.PeriS_IS_EN) && any(idx.Cue_IS_EN) && any(idx.PeriS_CS_EN) && any(idx.Cue_CS_EN)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['VMI_periEN_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell((abs(cell2mat(tuning_per_unit_table(2:end,idx.PeriS_IS_EN))) + abs(cell2mat(tuning_per_unit_table(2:end,idx.PeriS_CS_EN))) - ...
            abs(cell2mat(tuning_per_unit_table(2:end,idx.Cue_IS_EN))) - abs(cell2mat(tuning_per_unit_table(2:end,idx.Cue_CS_EN))))./...
            (abs(cell2mat(tuning_per_unit_table(2:end,idx.PeriS_IS_EN))) + abs(cell2mat(tuning_per_unit_table(2:end,idx.PeriS_CS_EN))) + ...
            abs(cell2mat(tuning_per_unit_table(2:end,idx.Cue_IS_EN))) + abs(cell2mat(tuning_per_unit_table(2:end,idx.Cue_CS_EN)))));
    end
    
    if any(idx.TIhol_IS_FR) && any(idx.Cue_IS_FR)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['VMI_post_IS_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell((cell2mat(tuning_per_unit_table(2:end,idx.TIhol_IS_FR)) - cell2mat(tuning_per_unit_table(2:end,idx.Cue_IS_FR)))./...
            (cell2mat(tuning_per_unit_table(2:end,idx.TIhol_IS_FR)) + cell2mat(tuning_per_unit_table(2:end,idx.Cue_IS_FR))));
    end
    if any(idx.TIhol_CS_FR) && any(idx.Cue_CS_FR)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['VMI_post_CS_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell((cell2mat(tuning_per_unit_table(2:end,idx.TIhol_CS_FR)) - cell2mat(tuning_per_unit_table(2:end,idx.Cue_CS_FR)))./...
            (cell2mat(tuning_per_unit_table(2:end,idx.TIhol_CS_FR)) + cell2mat(tuning_per_unit_table(2:end,idx.Cue_CS_FR))));
    end
    
    
    
    
    if any(idx.TIhol_IS_EN) && any(idx.Cue_IS_EN)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['VMI_postEN_IS_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell((abs(cell2mat(tuning_per_unit_table(2:end,idx.TIhol_IS_EN))) - abs(cell2mat(tuning_per_unit_table(2:end,idx.Cue_IS_EN))))./...
            (abs(cell2mat(tuning_per_unit_table(2:end,idx.TIhol_IS_EN))) + abs(cell2mat(tuning_per_unit_table(2:end,idx.Cue_IS_EN)))));
    end
    if any(idx.TIhol_CS_EN) && any(idx.Cue_CS_EN)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['VMI_postEN_CS_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell((abs(cell2mat(tuning_per_unit_table(2:end,idx.TIhol_CS_EN))) - abs(cell2mat(tuning_per_unit_table(2:end,idx.Cue_CS_EN))))./...
            (abs(cell2mat(tuning_per_unit_table(2:end,idx.TIhol_CS_EN))) + abs(cell2mat(tuning_per_unit_table(2:end,idx.Cue_CS_EN)))));
    end
    
    if any(idx.TIhol_IS_EN) && any(idx.Cue_IS_EN) && any(idx.TIhol_CS_EN) && any(idx.Cue_CS_EN)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['VMI_postEN_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell((abs(cell2mat(tuning_per_unit_table(2:end,idx.TIhol_IS_EN))) + abs(cell2mat(tuning_per_unit_table(2:end,idx.TIhol_CS_EN))) - ...
            abs(cell2mat(tuning_per_unit_table(2:end,idx.Cue_IS_EN))) - abs(cell2mat(tuning_per_unit_table(2:end,idx.Cue_CS_EN))))./...
            (abs(cell2mat(tuning_per_unit_table(2:end,idx.TIhol_IS_EN))) + abs(cell2mat(tuning_per_unit_table(2:end,idx.TIhol_CS_EN))) + ...
            abs(cell2mat(tuning_per_unit_table(2:end,idx.Cue_IS_EN))) + abs(cell2mat(tuning_per_unit_table(2:end,idx.Cue_CS_EN)))));
        %     tuning_per_unit_table(2:end,n_column)=num2cell((abs(cell2mat(tuning_per_unit_table(2:end,idx.TIhol_IS_EN))) + abs(cell2mat(tuning_per_unit_table(2:end,idx.Cue_IS_EN))) + ...
        %                                                     abs(cell2mat(tuning_per_unit_table(2:end,idx.TIhol_CS_EN))) + abs(cell2mat(tuning_per_unit_table(2:end,idx.Cue_CS_EN))))/4);
    end
    
    
    if any(idx.CueG_IS_EN) && any(idx.PreG_IS_EN)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['PreCueSum_IS_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell(cell2mat(tuning_per_unit_table(2:end,idx.CueG_IS_EN)) + cell2mat(tuning_per_unit_table(2:end,idx.PreG_IS_EN)));
    end
    
    if any(idx.CueG_CS_EN) && any(idx.PreG_CS_EN)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['PreCueSum_CS_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell(cell2mat(tuning_per_unit_table(2:end,idx.CueG_CS_EN)) + cell2mat(tuning_per_unit_table(2:end,idx.PreG_CS_EN)));
    end
    
    
    if any(idx.CueG_IS_FR) && any(idx.PreG_IS_FR)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['PreCueMean_IS_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell((cell2mat(tuning_per_unit_table(2:end,idx.CueG_IS_FR)) + cell2mat(tuning_per_unit_table(2:end,idx.PreG_IS_FR)))/2);
    end
    
    if any(idx.CueG_CS_FR) && any(idx.PreG_CS_FR)
        n_column=n_column+1;
        tuning_per_unit_table{1,n_column}=['PreCueMean_CS_' taskcase];
        tuning_per_unit_table(2:end,n_column)=num2cell((cell2mat(tuning_per_unit_table(2:end,idx.CueG_CS_FR)) + cell2mat(tuning_per_unit_table(2:end,idx.PreG_CS_FR)))/2); %% /2 !! for mean
    end
    
    
    for e=1:numel(epochs)
        if any(idx.([epochs{e} '_in_space_DF'])) && any(idx.([epochs{e} '_in_IS_FR'])) && any(idx.([epochs{e} '_in_CS_FR'])) && any(idx.([epochs{e} '_ch_IS_FR'])) && any(idx.([epochs{e} '_ch_CS_FR']))
            n_column=n_column+1;
            tuning_per_unit_table{1,n_column}=['in_' epochs{e} '_spaceCI_IX_' taskcase];
            tuning_per_unit_table(2:end,n_column)=num2cell(cell2mat(tuning_per_unit_table(2:end,idx.([epochs{e} '_in_space_DF'])))./...
                (cell2mat(tuning_per_unit_table(2:end,idx.([epochs{e} '_in_IS_FR']))) + cell2mat(tuning_per_unit_table(2:end,idx.([epochs{e} '_in_CS_FR']))) +...
                cell2mat(tuning_per_unit_table(2:end,idx.([epochs{e} '_ch_IS_FR']))) + cell2mat(tuning_per_unit_table(2:end,idx.([epochs{e} '_ch_CS_FR']))))); %% /2 !! for mean
        end
    end
    
    
    for e=1:numel(epochs)
        
        if any(idx.([epochs{e} '_ch_space_DF'])) && any(idx.([epochs{e} '_in_IS_FR'])) && any(idx.([epochs{e} '_in_CS_FR'])) && any(idx.([epochs{e} '_ch_IS_FR'])) && any(idx.([epochs{e} '_ch_CS_FR']))
            
            n_column=n_column+1;
            tuning_per_unit_table{1,n_column}=['ch_' epochs{e} '_spaceCI_IX_' taskcase];
            tuning_per_unit_table(2:end,n_column)=num2cell(cell2mat(tuning_per_unit_table(2:end,idx.([epochs{e} '_ch_space_DF'])))./...
                (cell2mat(tuning_per_unit_table(2:end,idx.([epochs{e} '_in_IS_FR']))) + cell2mat(tuning_per_unit_table(2:end,idx.([epochs{e} '_in_CS_FR']))) +...
                cell2mat(tuning_per_unit_table(2:end,idx.([epochs{e} '_ch_IS_FR']))) + cell2mat(tuning_per_unit_table(2:end,idx.([epochs{e} '_ch_CS_FR']))))); %% /2 !! for mean
        end
        
        
        
        if any(idx.([epochs{e} '_in_IH_space'])) && any(idx.([epochs{e} '_in_CH_space'])) && any(idx.([epochs{e} '_in_IS_hands'])) && any(idx.([epochs{e} '_in_CS_hands']))
            
            
            IH_CS=strcmp(tuning_per_unit_table(:,idx.([epochs{e} '_in_IH_space'])),'CS');
            IH_IS=strcmp(tuning_per_unit_table(:,idx.([epochs{e} '_in_IH_space'])),'IS');
            CH_CS=strcmp(tuning_per_unit_table(:,idx.([epochs{e} '_in_CH_space'])),'CS');
            CH_IS=strcmp(tuning_per_unit_table(:,idx.([epochs{e} '_in_CH_space'])),'IS');
            IS_CH=strcmp(tuning_per_unit_table(:,idx.([epochs{e} '_in_IS_hands'])),'CH');
            IS_IH=strcmp(tuning_per_unit_table(:,idx.([epochs{e} '_in_IS_hands'])),'IH');
            CS_CH=strcmp(tuning_per_unit_table(:,idx.([epochs{e} '_in_CS_hands'])),'CH');
            CS_IH=strcmp(tuning_per_unit_table(:,idx.([epochs{e} '_in_CS_hands'])),'IH');
            incongruent_space=  (IH_CS & CH_IS) | (IH_IS & CH_CS);
            incongruent_hands=  (IS_CH & CS_IH) | (IS_IH & CS_CH);
            
            CS= (IH_CS | CH_CS) & ~incongruent_space;
            IS= (IH_IS | CH_IS) & ~incongruent_space;
            CH= (IS_CH | CS_CH) & ~incongruent_hands;
            IH= (IS_IH | CS_IH) & ~incongruent_hands;
            
            
            n_column=n_column+1;
            tuning_per_unit_table(CS,n_column)={'CS'};
            tuning_per_unit_table(IS,n_column)={'IS'};
            tuning_per_unit_table(~(CS|IS),n_column)={'-'};
            tuning_per_unit_table(incongruent_space,n_column)={'incongruent'};
            tuning_per_unit_table{1,n_column}=['in_' epochs{e} '_space_perhand_' taskcase];
            
            
            n_column=n_column+1;
            tuning_per_unit_table(CH,n_column)={'CH'};
            tuning_per_unit_table(IH,n_column)={'IH'};
            tuning_per_unit_table(~(CH|IH),n_column)={'-'};
            tuning_per_unit_table(incongruent_hands,n_column)={'incongruent'};
            tuning_per_unit_table{1,n_column}=['in_' epochs{e} '_hands_perspace_' taskcase];
        end
    end
end

%% across tasks combinations
idx.preS_spaceperhand=DAG_find_column_index(tuning_per_unit_table,'in_PreS_space_perhand_Ddsa_han');
idx.preR_spaceperhand=DAG_find_column_index(tuning_per_unit_table,'in_PreR_space_perhand_Ddre_han');
if any(idx.preS_spaceperhand) && any(idx.preR_spaceperhand)
    
    R_CS=strcmp(tuning_per_unit_table(:,idx.preR_spaceperhand),'CS');
    R_IS=strcmp(tuning_per_unit_table(:,idx.preR_spaceperhand),'IS');
    S_CS=strcmp(tuning_per_unit_table(:,idx.preS_spaceperhand),'CS');
    S_IS=strcmp(tuning_per_unit_table(:,idx.preS_spaceperhand),'IS');
    incongruent_space=  (R_CS & S_IS) | (R_IS & S_CS);
    CS= (R_CS | S_CS) & ~incongruent_space;
    IS= (R_IS | S_IS) & ~incongruent_space;
    
    
    n_column=n_column+1;
    tuning_per_unit_table(CS | IS,n_column)={'tuned'};
    %tuning_per_unit_table(IS,n_column)={'IS'};
    tuning_per_unit_table(~(CS|IS),n_column)={'-'};
    tuning_per_unit_table(incongruent_space,n_column)={'incongruent'};
    tuning_per_unit_table{1,n_column}='in_Pre_space_perhand_Ddre_or_Ddsa';
end

idx.preS_handperspace=DAG_find_column_index(tuning_per_unit_table,'in_PreS_hands_perspace_Ddsa_han');
idx.preR_handperspace=DAG_find_column_index(tuning_per_unit_table,'in_PreR_hands_perspace_Ddre_han');
if any(idx.preS_handperspace) && any(idx.preR_handperspace)
    
    R_CH=strcmp(tuning_per_unit_table(:,idx.preR_handperspace),'CH');
    R_IH=strcmp(tuning_per_unit_table(:,idx.preR_handperspace),'IH');
    S_CH=strcmp(tuning_per_unit_table(:,idx.preS_handperspace),'CH');
    S_IH=strcmp(tuning_per_unit_table(:,idx.preS_handperspace),'IH');
    incongruent_space=  (R_CH & S_IH) | (R_IH & S_CH);
    CH= (R_CH | S_CH) & ~incongruent_space;
    IH= (R_IH | S_IH) & ~incongruent_space;
    
    
    n_column=n_column+1;
    tuning_per_unit_table(CH | IH,n_column)={'tuned'};
    %tuning_per_unit_table(IS,n_column)={'IS'};
    tuning_per_unit_table(~(CH|IH),n_column)={'-'};
    tuning_per_unit_table(incongruent_space,n_column)={'incongruent'};
    tuning_per_unit_table{1,n_column}='in_Pre_hands_perspace_Ddre_or_Ddsa';
end

idx.CueS_spaceperhand=DAG_find_column_index(tuning_per_unit_table,'in_Cue_space_perhand_Ddsa_han');
idx.CueR_spaceperhand=DAG_find_column_index(tuning_per_unit_table,'in_Cue_space_perhand_Ddre_han');
if any(idx.preS_spaceperhand) && any(idx.preR_spaceperhand)
    
    R_CS=strcmp(tuning_per_unit_table(:,idx.CueR_spaceperhand),'CS');
    R_IS=strcmp(tuning_per_unit_table(:,idx.CueR_spaceperhand),'IS');
    S_CS=strcmp(tuning_per_unit_table(:,idx.CueS_spaceperhand),'CS');
    S_IS=strcmp(tuning_per_unit_table(:,idx.CueS_spaceperhand),'IS');
    incongruent_space=  (R_CS & S_IS) | (R_IS & S_CS);
    CS= (R_CS | S_CS) & ~incongruent_space;
    IS= (R_IS | S_IS) & ~incongruent_space;
    
    
    n_column=n_column+1;
    tuning_per_unit_table(CS | IS,n_column)={'tuned'};
    %tuning_per_unit_table(IS,n_column)={'IS'};
    tuning_per_unit_table(~(CS|IS),n_column)={'-'};
    tuning_per_unit_table(incongruent_space,n_column)={'incongruent'};
    tuning_per_unit_table{1,n_column}='in_Cue_space_perhand_Ddre_or_Ddsa';
end

end