function TT=ph_load_extended_tuning_table(keys)

%% Reading in table

TT={};
load(keys.anova_table_file);
TT=tuning_per_unit_table;

% fix blanks in column titles
for k=1:size(TT,2)
    blanks=strfind(TT{1,k},' ');
    if any(blanks)
        TT{1,k}(blanks) ='_';
    end
end

TT(:,end+1)=num2cell(repmat('Y',size(TT,1),1));
TT(1,end)={'ungrouped'};
taskcaseexistingindex=find(~cellfun(@isempty,strfind(TT(1,:),'_trials_per_condition')));
taskcaseexistingindex_hf=find(~cellfun(@isempty,strfind(TT(1,:),'_trials_per_hemifield')));
taskcaseexistingindex_chf=find(~cellfun(@isempty,strfind(TT(1,:),'_trials_per_congruent_hand_hemifield')));
taskcaseexistingindex_tot=find(~cellfun(@isempty,strfind(TT(1,:),'_trials_total')));
if ~isempty(taskcaseexistingindex)
    taskcases={};
    %% Trial criterion per position
    for I=taskcaseexistingindex
        if ~((any(strfind(TT{1,I},'in_')==1) && strcmp(keys.tt.trial_criterion_in,'per_position')) ||...
                (any(strfind(TT{1,I},'ch_')==1) && strcmp(keys.tt.trial_criterion_ch,'per_position')) )
            continue;
        end
        task_existing_column=num2cell(cellfun(@(x) ~isempty(x)&&~ischar(x)&&x>=keys.cal.min_trials_per_condition,TT(:,I)));
        strpos=strfind(TT{1,I},'_trials_per_condition');
        task_existing_column{1,1}=['existing_' TT{1,I}([1:strpos,strpos+22:end])];
        TT(:,end+1)=task_existing_column;
        taskcases = [taskcases, {TT{1,I}(end-7:end)}];
    end
    %% trial criterion per hemifield (especially useful for choices)
    for I=taskcaseexistingindex_hf
        if ~((any(strfind(TT{1,I},'in_')==1) && strcmp(keys.tt.trial_criterion_in,'per_hemifield')) ||...
                (any(strfind(TT{1,I},'ch_')==1) && strcmp(keys.tt.trial_criterion_ch,'per_hemifield')))
            continue;
        end
        task_existing_column=num2cell(cellfun(@(x) ~isempty(x)&&~ischar(x)&&x>=keys.cal.min_trials_per_condition,TT(:,I)));
        strpos=strfind(TT{1,I},'_trials_per_hemifield');
        task_existing_column{1,1}=['existing_' TT{1,I}([1:strpos,strpos+22:end])];
        TT(:,end+1)=task_existing_column;
        taskcases = [taskcases, {TT{1,I}(end-7:end)}];
    end
    %%  To analyse only congruent hand\space condition
    for I=taskcaseexistingindex_chf
        if ~((any(strfind(TT{1,I},'in_')==1) && strcmp(keys.tt.trial_criterion_in,'per_congruent_hand_hemifield')) ||...
                (any(strfind(TT{1,I},'ch_')==1) && strcmp(keys.tt.trial_criterion_ch,'per_congruent_hand_hemifield')))
            continue;
        end
        task_existing_column=num2cell(cellfun(@(x) ~isempty(x)&&~ischar(x)&&x>=keys.cal.min_trials_per_condition,TT(:,I)));
        strpos=strfind(TT{1,I},'_trials_per_congruent_hand_hemifield');
        task_existing_column{1,1}=['existing_' TT{1,I}([1:strpos,strpos+37:end])];
        TT(:,end+1)=task_existing_column;
        taskcases = [taskcases, {TT{1,I}(end-7:end)}];
    end
    
    %% trial criterion per hemifield separately for perturbation group (include if any condition is fullfilled, not ideal, but works for effector = 0)
    if strcmp(keys.tt.trial_criterion_in,'per_hemifield_and_perturbation')
        counter = 0;
        for I=taskcaseexistingindex_hf
            if ~((any(strfind(TT{1,I},'in_')==1) && strcmp(keys.tt.trial_criterion_in,'per_hemifield_and_perturbation')) ||...
                    (any(strfind(TT{1,I},'ch_')==1) && strcmp(keys.tt.trial_criterion_ch,'per_hemifield_and_perturbation')))
                continue;
            end
            counter = counter + 1;
            task_existing_column=num2cell(cellfun(@(x) ~isempty(x)&&~ischar(x)&&x>=keys.cal.min_trials_per_condition,TT(:,I)));
            strpos=strfind(TT{1,I},'_trials_per_hemifield');
            task_existing_column{1,1}=['existing_' TT{1,I}([1:strpos,strpos+22:end])];
            TT(:,end+1)=task_existing_column;
            taskcases = [taskcases, {TT{1,I}(end-7:end)}];
        end
        %this nasty line put 1 to all to existing column if at least 1 is
        %equal to 1, for each unit separately
        TT(2:end,end-(counter-1):end) = repmat(num2cell(any(cellfun(@(x) x == 1, TT(2:end,end-(counter-1):end)),2)),...
            1, size(TT(2:end,end-(counter-1):end),2));
        
    end
    
    %%  General minimum number of trials criterion (across all hand/hemifield/position conditions)
    for I=taskcaseexistingindex_tot
        if ~((any(strfind(TT{1,I},'in_')==1) && strcmp(keys.tt.trial_criterion_in,'total')) ||...
                (any(strfind(TT{1,I},'ch_')==1) && strcmp(keys.tt.trial_criterion_ch,'total')))
            continue;
        end
        task_existing_column=num2cell(cellfun(@(x) ~isempty(x)&&~ischar(x)&&x>=keys.cal.min_trials_per_condition,TT(:,I)));
        strpos=strfind(TT{1,I},'_trials_total');
        task_existing_column{1,1}=['existing_' TT{1,I}([1:strpos,strpos+14:end])];
        TT(:,end+1)=task_existing_column;
        taskcases = [taskcases, {TT{1,I}(end-7:end)}];
    end
else % this is for previous version (backwards compatibility) ... not needed any more??
    taskcaseexistingindex=find(~cellfun(@isempty,strfind(TT(1,:),'in_epoch_main'))); %not ideal
    taskcases={};
    for I=taskcaseexistingindex
        task_existing_column=num2cell(~cellfun(@isempty,TT(:,I)));
        task_existing_column{1,1}=['existing' TT{1,I}(14:end)];
        TT(:,end+1)=task_existing_column;
        taskcases = [taskcases, {TT{1,I}(15:end)}];
    end
    %% to be combined with previous
    taskcaseexistingindex=find(~cellfun(@isempty,strfind(TT(1,:),'ch_epoch_main'))); %not ideal
    %taskcases={};
    for I=taskcaseexistingindex
        task_existing_column=num2cell(~cellfun(@isempty,TT(:,I)));
        task_existing_column{1,1}=['existing' TT{1,I}(1:3) TT{1,I}(14:end)];
        TT(:,end+1)=task_existing_column;
        %taskcases = [taskcases, {TT{1,I}(15:end)}];
    end
    
end



%% create space OR interaction column
taskcases=unique(taskcases);
epochs={'Fhol','Cue','Cue2','EDel','Del','MemE','MemL','PreS','PreR','PeriS','PeriR','Pre2','Peri2','PreG','CueG','TIhol','THol'};
for t=1:numel(taskcases)
    taskcase=taskcases{t};
    idx.space   =DAG_find_column_index(TT,[keys.tt.IC_for_criterion '_spaceLR_main_'   taskcase] );
    idx.epoch   =DAG_find_column_index(TT,[keys.tt.IC_for_criterion '_epoch_main_'     taskcase] );
    idx.hands   =DAG_find_column_index(TT,[keys.tt.IC_for_criterion '_hands_main_'     taskcase] );
    idx.SxH     =DAG_find_column_index(TT,[keys.tt.IC_for_criterion '_SxH_'            taskcase] );
    idx.ExS     =DAG_find_column_index(TT,[keys.tt.IC_for_criterion '_ExS_'            taskcase] );
    idx.ExH     =DAG_find_column_index(TT,[keys.tt.IC_for_criterion '_ExH_'            taskcase] );
    idx.ESH     =DAG_find_column_index(TT,[keys.tt.IC_for_criterion '_ExSxH_'          taskcase] );
    
    
    for e=1:numel(epochs)
        idx.(epochs{e})                     =DAG_find_column_index(TT,[keys.tt.IC_for_criterion '_AH_' epochs{e} '_epoch_'   taskcase] );
        idx.([epochs{e} '_space_DF'])       =DAG_find_column_index(TT,[keys.tt.IC_for_criterion epochs{e} '_spaceLR_DF_'   taskcase] );
        idx.([epochs{e} '_IS_FR'])          =DAG_find_column_index(TT,[keys.tt.IC_for_criterion '_AH_IS_' epochs{e} '_epoch_FR_'  taskcase] );
        idx.([epochs{e} '_CS_FR'])          =DAG_find_column_index(TT,[keys.tt.IC_for_criterion '_AH_CS_' epochs{e} '_epoch_FR_'  taskcase] );
        idx.([epochs{e} '_IS_EN'])          =DAG_find_column_index(TT,[keys.tt.IC_for_criterion '_AH_IS_' epochs{e} '_epoch_DF_'  taskcase] );
        idx.([epochs{e} '_CS_EN'])          =DAG_find_column_index(TT,[keys.tt.IC_for_criterion '_AH_CS_' epochs{e} '_epoch_DF_'  taskcase] );
        idx.([epochs{e} '_in_space_DF'])    =DAG_find_column_index(TT,['in_' epochs{e} '_spaceLR_DF_'   taskcase] );
        idx.([epochs{e} '_ch_space_DF'])    =DAG_find_column_index(TT,['ch_' epochs{e} '_spaceLR_DF_'   taskcase] );
        idx.([epochs{e} '_in_IS_FR'])       =DAG_find_column_index(TT,['in_AH_IS_' epochs{e} '_epoch_FR_'  taskcase] );
        idx.([epochs{e} '_in_CS_FR'])       =DAG_find_column_index(TT,['in_AH_CS_' epochs{e} '_epoch_FR_'  taskcase] );
        idx.([epochs{e} '_ch_IS_FR'])       =DAG_find_column_index(TT,['ch_AH_IS_' epochs{e} '_epoch_FR_'  taskcase] );
        idx.([epochs{e} '_ch_CS_FR'])       =DAG_find_column_index(TT,['ch_AH_CS_' epochs{e} '_epoch_FR_'  taskcase] );
        
        
        idx.([epochs{e} '_in_IH_space'])       =DAG_find_column_index(TT,['in_IH_' epochs{e} '_spaceLR_'  taskcase] );
        idx.([epochs{e} '_in_CH_space'])       =DAG_find_column_index(TT,['in_CH_' epochs{e} '_spaceLR_'  taskcase] );
        idx.([epochs{e} '_in_IS_hands'])       =DAG_find_column_index(TT,['in_IS_' epochs{e} '_hands_'  taskcase] );
        idx.([epochs{e} '_in_CS_hands'])       =DAG_find_column_index(TT,['in_CS_' epochs{e} '_hands_'  taskcase] );
        
    end
    
    TT=cleanup_empties(TT,idx,'epoch','double');
    TT=cleanup_empties(TT,idx,'space','char');
    TT=cleanup_empties(TT,idx,'hands','char');
    TT=cleanup_empties(TT,idx,'SxH','char');
    TT=cleanup_empties(TT,idx,'ExS','double');
    TT=cleanup_empties(TT,idx,'ExH','double');
    TT=cleanup_empties(TT,idx,'ESH','double');
    %
    %     if ~isempty(idx.epoch); TT([false; cellfun(@(x) isempty(x) ,TT(2:end,idx.epoch))],idx.epoch)  ={NaN}; end;
    %     if ~isempty(idx.space); TT([false; cellfun(@(x) isempty(x) ,TT(2:end,idx.space))],idx.space)  ={''}; end;
    %     if ~isempty(idx.hands); TT([false; cellfun(@(x) isempty(x) ,TT(2:end,idx.hands))],idx.hands)  ={''}; end;
    %     if ~isempty(idx.SxH);   TT([false; cellfun(@(x) isempty(x) ,TT(2:end,idx.SxH))],idx.SxH)      ={''}; end;
    %     if ~isempty(idx.ExS);   TT([false; cellfun(@(x) isempty(x) ,TT(2:end,idx.ExS))],idx.ExS)      ={NaN}; end;
    %     if ~isempty(idx.ExH);   TT([false; cellfun(@(x) isempty(x) ,TT(2:end,idx.ExH))],idx.ExH)      ={NaN}; end;
    %     if ~isempty(idx.ESH);   TT([false; cellfun(@(x) isempty(x) ,TT(2:end,idx.ESH))],idx.ESH)      ={NaN}; end;
    
    for e=1:numel(epochs)
        TT=cleanup_empties(TT,idx,epochs{e});
        TT=cleanup_empties(TT,idx,[epochs{e} '_space_DF']);
        TT=cleanup_empties(TT,idx,[epochs{e} '_IS_FR']);
        TT=cleanup_empties(TT,idx,[epochs{e} '_CS_FR']);
        TT=cleanup_empties(TT,idx,[epochs{e} '_IS_EN']);
        TT=cleanup_empties(TT,idx,[epochs{e} '_CS_EN']);
        TT=cleanup_empties(TT,idx,[epochs{e} '_in_space_DF']);
        TT=cleanup_empties(TT,idx,[epochs{e} '_ch_space_DF']);
        TT=cleanup_empties(TT,idx,[epochs{e} '_in_IS_FR']);
        TT=cleanup_empties(TT,idx,[epochs{e} '_in_CS_FR']);
        TT=cleanup_empties(TT,idx,[epochs{e} '_ch_IS_FR']);
        TT=cleanup_empties(TT,idx,[epochs{e} '_ch_CS_FR']);
        TT=cleanup_empties(TT,idx,[epochs{e} '_in_IH_space']);
        TT=cleanup_empties(TT,idx,[epochs{e} '_in_CH_space']);
        TT=cleanup_empties(TT,idx,[epochs{e} '_in_IS_hands']);
        TT=cleanup_empties(TT,idx,[epochs{e} '_in_CS_hands']);
    end
    
    
    switch keys.tt.space_criterion %% consider only neurons that show 1) SxE or space main effect; 2) SxE; 3) take all (no criterion)?
        case 'interaction or space only'
            space_or_interaction=([TT{2:end,idx.ExS}]==1)' | ~ismember(TT(2:end,idx.space),'-');
        case 'interaction only'
            space_or_interaction=([TT{2:end,idx.ExS}]==1)';
        case 'none'
            space_or_interaction=true(size(TT,1)-1,1);
    end
    
    %%% goofy
    switch keys.tt.epoch_criterion %% consider only neurons that show 1) SxE or space main effect; 2) SxE; 3) take all (no criterion)?
        case 'SxE or epoch only'
            epoch_or_interaction=([TT{2:end,idx.ExS}]==1)' | ([TT{2:end,idx.epoch}]==1)';
        case 'HxE or epoch only'
            epoch_or_interaction=([TT{2:end,idx.ExH}]==1)' | ([TT{2:end,idx.epoch}]==1)';
        case 'SHE or epoch only'
            epoch_or_interaction=([TT{2:end,idx.ESH}]==1)' | ([TT{2:end,idx.epoch}]==1)';
        case 'SxE only'
            epoch_or_interaction=([TT{2:end,idx.ExS}]==1)';
        case 'HxE only'
            epoch_or_interaction=([TT{2:end,idx.ExH}]==1)';
        case 'SHE only'
            epoch_or_interaction=([TT{2:end,idx.ESH}]==1)';
        case 'none'
            epoch_or_interaction=true(size(TT,1)-1,1);
    end
    switch keys.tt.hands_criterion %% consider only neurons that show 1) HxE or hand main effect; 2) HxE; 3) take all (no criterion)?
        case 'interaction or hands only'
            hands_or_interaction=([TT{2:end,idx.ExH}]==1)' | ~ismember(TT(2:end,idx.hands),'-');
        case 'interaction only'
            hands_or_interaction=([TT{2:end,idx.ExH}]==1)';
        case 'none'
            hands_or_interaction=true(size(TT,1)-1,1);
    end
    switch keys.tt.SXH_criterion %% consider only neurons that show 1) SxHxE or SxH; 2) SxHxE; 3) take all (no criterion)?
        case 'SXH SXE HXE SHE'
            SXH_or_interaction=([TT{2:end,idx.ESH}]==1)' | ~ismember(TT(2:end,idx.SxH),'-')...
                |([TT{2:end,idx.ExH}]==1)' | ([TT{2:end,idx.ExS}]==1)';
        case 'interaction or SXH only'
            SXH_or_interaction=([TT{2:end,idx.ESH}]==1)' | ~ismember(TT(2:end,idx.SxH),'-');
        case 'interaction only'
            SXH_or_interaction=([TT{2:end,idx.ESH}]==1)';
        case 'none'
            SXH_or_interaction=true(size(TT,1)-1,1);
    end
    n_column=size(TT,2);
    
    n_column=n_column+1;
    TT{1,n_column}=['space_or_interaction_' taskcase];
    TT(2:end,n_column)=num2cell(space_or_interaction);
    
    n_column=n_column+1;
    TT{1,n_column}=['epoch_or_interaction_' taskcase];
    TT(2:end,n_column)=num2cell(epoch_or_interaction);
    
    n_column=n_column+1;
    TT{1,n_column}=['hands_or_interaction_' taskcase];
    TT(2:end,n_column)=num2cell(hands_or_interaction);
    
    n_column=n_column+1;
    TT{1,n_column}=['SXH_or_interaction_' taskcase];
    TT(2:end,n_column)=num2cell(SXH_or_interaction);
    
    %% adding columns for visual, motor, visuomotor, and fixation cells
    if any(idx.Fhol) && any(idx.PeriS) && any(idx.Cue) && any(idx.TIhol)
        TT{1,n_column+1}=['fixation_only_' taskcase];
        TT{1,n_column+2}=['fixation_and_sac_suppression_' taskcase];
        TT{1,n_column+3}=['Sac_supression_' taskcase];
        TT(2:end,n_column+1)=num2cell(ismember(TT(2:end,idx.Fhol),'en') & ismember(TT(2:end,idx.PeriS),'-')  & ~ismember(TT(2:end,idx.Cue),'en') & ~ismember(TT(2:end,idx.TIhol),'en'));
        TT(2:end,n_column+2)=num2cell(ismember(TT(2:end,idx.Fhol),'en') & ismember(TT(2:end,idx.PeriS),'su') & ~ismember(TT(2:end,idx.Cue),'en') & ~ismember(TT(2:end,idx.TIhol),'en'));
        TT(2:end,n_column+3)=num2cell(ismember(TT(2:end,idx.Fhol),'-')  & ismember(TT(2:end,idx.PeriS),'su') & ~ismember(TT(2:end,idx.Cue),'en') & ~ismember(TT(2:end,idx.TIhol),'en'));
        n_column=n_column+3;
    end
    
    TT=get_VM(TT,idx.Cue,idx.TIhol,'first',['visual_' taskcase]);
    TT=get_VM(TT,idx.Cue,idx.TIhol,'second',['motor_' taskcase]);
    TT=get_VM(TT,idx.Cue,idx.TIhol,'both',['visuomotor_' taskcase]);
    TT=get_VM(TT,idx.Cue,idx.TIhol,'none',['notclassified_' taskcase]);
    
    TT=get_VM(TT,idx.Cue,idx.PreS,'first',['visual_pre_' taskcase]);
    TT=get_VM(TT,idx.Cue,idx.PreS,'second',['motor_pre_' taskcase]);
    TT=get_VM(TT,idx.Cue,idx.PreS,'both',['visuomotor_pre_' taskcase]);
    
    TT=get_VM(TT,idx.Cue,idx.Pre2,'first',['visual_pre2_' taskcase]);
    TT=get_VM(TT,idx.Cue,idx.Pre2,'second',['motor_pre2_' taskcase]);
    TT=get_VM(TT,idx.Cue,idx.Pre2,'both',['visuomotor_pre2_' taskcase]);
    
    TT=get_VM(TT,idx.Cue,idx.PeriS,'first',['visual_peri_' taskcase]);
    TT=get_VM(TT,idx.Cue,idx.PeriS,'second',['motor_peri_' taskcase]);
    TT=get_VM(TT,idx.Cue,idx.PeriS,'both',['visuomotor_peri_' taskcase]);
    
    
    TT=get_VMI(TT,idx.PeriS_IS_FR,idx.Cue_IS_FR,['VMI_peri_IS_' taskcase],'signed');
    TT=get_VMI(TT,idx.PeriS_CS_FR,idx.Cue_CS_FR,['VMI_peri_CS_' taskcase],'signed');
    TT=get_VMI(TT,idx.PeriS_IS_EN,idx.Cue_IS_EN,['VMI_periEN_IS_' taskcase],'absolute');
    TT=get_VMI(TT,idx.PeriS_CS_EN,idx.Cue_CS_EN,['VMI_periEN_CS_' taskcase],'absolute');
    
    TT=get_VMI(TT,idx.TIhol_IS_FR,idx.Cue_IS_FR,['VMI_post_IS_' taskcase],'signed');
    TT=get_VMI(TT,idx.TIhol_CS_FR,idx.Cue_CS_FR,['VMI_post_CS_' taskcase],'signed');
    TT=get_VMI(TT,idx.TIhol_IS_EN,idx.Cue_IS_EN,['VMI_postEN_IS_' taskcase],'absolute');
    TT=get_VMI(TT,idx.TIhol_CS_EN,idx.Cue_CS_EN,['VMI_postEN_CS_' taskcase],'absolute');
    
    %% MP
    
%     TT=get_VM(TT,idx.Cue,idx.Del,'first',['visual_only_' taskcase]); %% doesnt use suppression though !
%     TT=get_VM(TT,idx.Cue,idx.Del,'both',['visuomotor_' taskcase]); %% doesnt use suppression though ! --> does use it for this one Oo
%     TT=get_VM(TT,idx.Cue,idx.Del,'second',['motor_only_' taskcase]); %% doesnt use suppression though !
%     TT=get_VM(TT,idx.Cue,idx.Del,'first',['visual_only_' taskcase]); %% doesnt use suppression though !
%     TT=get_VM(TT,idx.Cue,idx.Del,'first',['visual_only_' taskcase]); %% doesnt use suppression though !
    
    
    n_column=size(TT,2);
    if any(idx.Cue) % && any(idx.TIhol)
        n_column=n_column+1;
        TT{1,n_column}=['visual_en_' taskcase];
        TT(2:end,n_column)=num2cell(ismember(TT(2:end,idx.Cue),{'en','bi'}));% & ~ismember(TT(2:end,idx.TIhol),{'en','su','bi'}));
    end
    
    if  any(idx.Del)
        n_column=n_column+1;
        TT{1,n_column}=['motor_en_' taskcase];
        TT(2:end,n_column)=num2cell(ismember(TT(2:end,idx.Del),{'en','bi'})); %& ismember(TT(2:end,idx.Del),{'en','su','bi'}));
    end
    
    
    %     if any(idx.Cue)  && any(idx.Del)
    %         n_column=n_column+1;
    %         tuning_per_unit_table{1,n_column}=['visual_only_' taskcase];
    %         tuning_per_unit_table(2:end,n_column)=num2cell(ismember(tuning_per_unit_table(2:end,idx.Cue),{'en','bi'}) & ~ismember(tuning_per_unit_table(2:end,idx.Del),{'en','bi'}));
    %     end
    %
    %     if any(idx.Cue)  && any(idx.Del)
    %         n_column=n_column+1;
    %         tuning_per_unit_table{1,n_column}=['visuomotor_' taskcase];
    %         tuning_per_unit_table(2:end,n_column)=num2cell(ismember(tuning_per_unit_table(2:end,idx.Cue),{'en','bi'}) & ismember(tuning_per_unit_table(2:end,idx.Del),{'en','bi'}));
    %     end
    %
    %     if any(idx.Cue) && any(idx.Del)
    %         n_column=n_column+1;
    %         tuning_per_unit_table{1,n_column}=['visuomotor_' taskcase];
    %         tuning_per_unit_table(2:end,n_column)=num2cell(ismember(tuning_per_unit_table(2:end,idx.Cue),{'en','su','bi'}) & ismember(tuning_per_unit_table(2:end,idx.Del),{'en','su','bi'}));
    %     end
    %
    %     if any(idx.Cue)  && any(idx.Del)
    %         n_column=n_column+1;
    %         tuning_per_unit_table{1,n_column}=['motor_only_' taskcase];
    %         tuning_per_unit_table(2:end,n_column)=num2cell(~ismember(tuning_per_unit_table(2:end,idx.Cue),{'en','bi'}) & ismember(tuning_per_unit_table(2:end,idx.Del),{'en','bi'}));
    %     end
    %
    
    %% MP end
    
    
    %% why this one was there for MP Mods?
    %     if any(idx.PeriS_IS_FR) && any(idx.Cue_IS_FR)
    %         n_column=n_column+1;
    %         tuning_per_unit_table{1,n_column}=['VMI_peri_IS_' taskcase];
    %         tuning_per_unit_table(2:end,n_column)=num2cell((cell2mat(tuning_per_unit_table(2:end,idx.PeriS_IS_FR)) - cell2mat(tuning_per_unit_table(2:end,idx.Cue_IS_FR)))./...
    %             (cell2mat(tuning_per_unit_table(2:end,idx.PeriS_IS_FR)) + cell2mat(tuning_per_unit_table(2:end,idx.Cue_IS_FR))));
    %     end
    
    %     if any(idx.Cue) && any(idx.TIhol)
    %         n_column=n_column+1;
    %         tuning_per_unit_table{1,n_column}=['notclassified_' taskcase];
    %         tuning_per_unit_table(2:end,n_column)=num2cell(~ismember(tuning_per_unit_table(2:end,idx.Cue),{'en','su','bi'}) & ~ismember(tuning_per_unit_table(2:end,idx.TIhol),{'en','su','bi'}));
    %     end
    %
    %     if any(idx.mem)
    %     n_column=n_column+1;
    %     tuning_per_unit_table{1,n_column}=['delay_' taskcase];
    %     tuning_per_unit_table(2:end,n_column)=num2cell(ismember(tuning_per_unit_table(2:end,idx.mem),'en') );
    %     end
    %
    %
    
    %     if any(idx.Cue) && any(idx.PeriS)
    %         n_column=n_column+1;
    %         tuning_per_unit_table{1,n_column}=['visual_peri_' taskcase];
    %         tuning_per_unit_table(2:end,n_column)=num2cell(ismember(tuning_per_unit_table(2:end,idx.Cue),{'en','su','bi'}) & ~ismember(tuning_per_unit_table(2:end,idx.PeriS),{'en','su','bi'}));
    %     end
    %
    %     if any(idx.Cue) && any(idx.PeriS)
    %         n_column=n_column+1;
    %         tuning_per_unit_table{1,n_column}=['motor_peri_' taskcase];
    %         tuning_per_unit_table(2:end,n_column)=num2cell(~ismember(tuning_per_unit_table(2:end,idx.Cue),{'en','su','bi'}) & ismember(tuning_per_unit_table(2:end,idx.PeriS),{'en','su','bi'}));
    %     end
    %
    %     if any(idx.Cue) && any(idx.PeriS)
    %         n_column=n_column+1;
    %         tuning_per_unit_table{1,n_column}=['visuomotor_peri_' taskcase];
    %         tuning_per_unit_table(2:end,n_column)=num2cell(ismember(tuning_per_unit_table(2:end,idx.Cue),{'en','su','bi'}) & ismember(tuning_per_unit_table(2:end,idx.PeriS),{'en','su','bi'}));
    %     end
    %
    %
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
    
    
    
    
    
    n_column=size(TT,2);
    if any(idx.PeriS_IS_EN) && any(idx.Cue_IS_EN) && any(idx.PeriS_CS_EN) && any(idx.Cue_CS_EN)
        n_column=n_column+1;
        TT{1,n_column}=['VMI_periEN_' taskcase];
        TT(2:end,n_column)=num2cell((abs(cell2mat(TT(2:end,idx.PeriS_IS_EN))) + abs(cell2mat(TT(2:end,idx.PeriS_CS_EN))) - ...
            abs(cell2mat(TT(2:end,idx.Cue_IS_EN))) - abs(cell2mat(TT(2:end,idx.Cue_CS_EN))))./...
            (abs(cell2mat(TT(2:end,idx.PeriS_IS_EN))) + abs(cell2mat(TT(2:end,idx.PeriS_CS_EN))) + ...
            abs(cell2mat(TT(2:end,idx.Cue_IS_EN))) + abs(cell2mat(TT(2:end,idx.Cue_CS_EN)))));
    end
    
    
    if any(idx.TIhol_IS_EN) && any(idx.Cue_IS_EN) && any(idx.TIhol_CS_EN) && any(idx.Cue_CS_EN)
        n_column=n_column+1;
        TT{1,n_column}=['VMI_postEN_' taskcase];
        TT(2:end,n_column)=num2cell((abs(cell2mat(TT(2:end,idx.TIhol_IS_EN))) + abs(cell2mat(TT(2:end,idx.TIhol_CS_EN))) - ...
            abs(cell2mat(TT(2:end,idx.Cue_IS_EN))) - abs(cell2mat(TT(2:end,idx.Cue_CS_EN))))./...
            (abs(cell2mat(TT(2:end,idx.TIhol_IS_EN))) + abs(cell2mat(TT(2:end,idx.TIhol_CS_EN))) + ...
            abs(cell2mat(TT(2:end,idx.Cue_IS_EN))) + abs(cell2mat(TT(2:end,idx.Cue_CS_EN)))));
    end
    
    
    if any(idx.CueG_IS_EN) && any(idx.PreG_IS_EN)
        n_column=n_column+1;
        TT{1,n_column}=['PreCueSum_IS_' taskcase];
        TT(2:end,n_column)=num2cell(cell2mat(TT(2:end,idx.CueG_IS_EN)) + cell2mat(TT(2:end,idx.PreG_IS_EN)));
    end
    
    if any(idx.CueG_CS_EN) && any(idx.PreG_CS_EN)
        n_column=n_column+1;
        TT{1,n_column}=['PreCueSum_CS_' taskcase];
        TT(2:end,n_column)=num2cell(cell2mat(TT(2:end,idx.CueG_CS_EN)) + cell2mat(TT(2:end,idx.PreG_CS_EN)));
    end
    
    
    if any(idx.CueG_IS_FR) && any(idx.PreG_IS_FR)
        n_column=n_column+1;
        TT{1,n_column}=['PreCueMean_IS_' taskcase];
        TT(2:end,n_column)=num2cell((cell2mat(TT(2:end,idx.CueG_IS_FR)) + cell2mat(TT(2:end,idx.PreG_IS_FR)))/2);
    end
    
    if any(idx.CueG_CS_FR) && any(idx.PreG_CS_FR)
        n_column=n_column+1;
        TT{1,n_column}=['PreCueMean_CS_' taskcase];
        TT(2:end,n_column)=num2cell((cell2mat(TT(2:end,idx.CueG_CS_FR)) + cell2mat(TT(2:end,idx.PreG_CS_FR)))/2); %% /2 !! for mean
    end
    
    
    
    for e=1:numel(epochs)
        if any(idx.([epochs{e} '_in_space_DF'])) && any(idx.([epochs{e} '_in_IS_FR'])) && any(idx.([epochs{e} '_in_CS_FR'])) && any(idx.([epochs{e} '_ch_IS_FR'])) && any(idx.([epochs{e} '_ch_CS_FR']))
            n_column=n_column+1;
            TT{1,n_column}=['in_' epochs{e} '_spaceCI_IX_' taskcase];
            TT(2:end,n_column)=num2cell(cell2mat(TT(2:end,idx.([epochs{e} '_in_space_DF'])))./...
                (cell2mat(TT(2:end,idx.([epochs{e} '_in_IS_FR']))) + cell2mat(TT(2:end,idx.([epochs{e} '_in_CS_FR']))) +...
                cell2mat(TT(2:end,idx.([epochs{e} '_ch_IS_FR']))) + cell2mat(TT(2:end,idx.([epochs{e} '_ch_CS_FR']))))); %% /2 !! for mean
        end
    end
    
    
    for e=1:numel(epochs)
        
        if any(idx.([epochs{e} '_ch_space_DF'])) && any(idx.([epochs{e} '_in_IS_FR'])) && any(idx.([epochs{e} '_in_CS_FR'])) && any(idx.([epochs{e} '_ch_IS_FR'])) && any(idx.([epochs{e} '_ch_CS_FR']))
            
            n_column=n_column+1;
            TT{1,n_column}=['ch_' epochs{e} '_spaceCI_IX_' taskcase];
            TT(2:end,n_column)=num2cell(cell2mat(TT(2:end,idx.([epochs{e} '_ch_space_DF'])))./...
                (cell2mat(TT(2:end,idx.([epochs{e} '_in_IS_FR']))) + cell2mat(TT(2:end,idx.([epochs{e} '_in_CS_FR']))) +...
                cell2mat(TT(2:end,idx.([epochs{e} '_ch_IS_FR']))) + cell2mat(TT(2:end,idx.([epochs{e} '_ch_CS_FR']))))); %% /2 !! for mean
        end
        
        
        if any(idx.([epochs{e} '_in_IH_space'])) && any(idx.([epochs{e} '_in_CH_space'])) && any(idx.([epochs{e} '_in_IS_hands'])) && any(idx.([epochs{e} '_in_CS_hands']))
            
            IH_CS=strcmp(TT(:,idx.([epochs{e} '_in_IH_space'])),'CS');
            IH_IS=strcmp(TT(:,idx.([epochs{e} '_in_IH_space'])),'IS');
            CH_CS=strcmp(TT(:,idx.([epochs{e} '_in_CH_space'])),'CS');
            CH_IS=strcmp(TT(:,idx.([epochs{e} '_in_CH_space'])),'IS');
            IS_CH=strcmp(TT(:,idx.([epochs{e} '_in_IS_hands'])),'CH');
            IS_IH=strcmp(TT(:,idx.([epochs{e} '_in_IS_hands'])),'IH');
            CS_CH=strcmp(TT(:,idx.([epochs{e} '_in_CS_hands'])),'CH');
            CS_IH=strcmp(TT(:,idx.([epochs{e} '_in_CS_hands'])),'IH');
            incongruent_space=  (IH_CS & CH_IS) | (IH_IS & CH_CS);
            incongruent_hands=  (IS_CH & CS_IH) | (IS_IH & CS_CH);
            
            CS= (IH_CS | CH_CS) & ~incongruent_space;
            IS= (IH_IS | CH_IS) & ~incongruent_space;
            CH= (IS_CH | CS_CH) & ~incongruent_hands;
            IH= (IS_IH | CS_IH) & ~incongruent_hands;
            
            n_column=n_column+1;
            TT(CS,n_column)={'CS'};
            TT(IS,n_column)={'IS'};
            TT(~(CS|IS),n_column)={'-'};
            TT(incongruent_space,n_column)={'incongruent'};
            TT{1,n_column}=['in_' epochs{e} '_space_perhand_' taskcase];
            
            n_column=n_column+1;
            TT(CH,n_column)={'CH'};
            TT(IH,n_column)={'IH'};
            TT(~(CH|IH),n_column)={'-'};
            TT(incongruent_hands,n_column)={'incongruent'};
            TT{1,n_column}=['in_' epochs{e} '_hands_perspace_' taskcase];
        end
    end
end

%% across tasks combinations
idx.CueS_spaceperhand=DAG_find_column_index(TT,'in_Cue_space_perhand_Ddsa_han');
idx.CueR_spaceperhand=DAG_find_column_index(TT,'in_Cue_space_perhand_Ddre_han');
idx.preS_spaceperhand=DAG_find_column_index(TT,'in_PreS_space_perhand_Ddsa_han');
idx.preR_spaceperhand=DAG_find_column_index(TT,'in_PreR_space_perhand_Ddre_han');
idx.preS_handperspace=DAG_find_column_index(TT,'in_PreS_hands_perspace_Ddsa_han');
idx.preR_handperspace=DAG_find_column_index(TT,'in_PreR_hands_perspace_Ddre_han');

TT=get_across_task(TT,idx.CueS_spaceperhand,idx.CueR_spaceperhand,'space','in_Cue_space_perhand_Ddre_or_Ddsa');
TT=get_across_task(TT,idx.preS_spaceperhand,idx.preR_spaceperhand,'space','in_Pre_space_perhand_Ddre_or_Ddsa');
TT=get_across_task(TT,idx.preS_handperspace,idx.preR_handperspace,'hands','in_Pre_hands_perspace_Ddre_or_Ddsa');
%
% if any(idx.preS_spaceperhand) && any(idx.preR_spaceperhand)
%     R_CS=strcmp(TT(:,idx.preR_spaceperhand),'CS');
%     R_IS=strcmp(TT(:,idx.preR_spaceperhand),'IS');
%     S_CS=strcmp(TT(:,idx.preS_spaceperhand),'CS');
%     S_IS=strcmp(TT(:,idx.preS_spaceperhand),'IS');
%     incongruent_space=  (R_CS & S_IS) | (R_IS & S_CS);
%     CS= (R_CS | S_CS) & ~incongruent_space;
%     IS= (R_IS | S_IS) & ~incongruent_space;
%     n_column=n_column+1;
%     TT(CS | IS,n_column)={'tuned'};
%     %TT(IS,n_column)={'IS'};
%     TT(~(CS|IS),n_column)={'-'};
%     TT(incongruent_space,n_column)={'incongruent'};
%     TT{1,n_column}='in_Pre_space_perhand_Ddre_or_Ddsa';
% end
%
% idx.preS_handperspace=DAG_find_column_index(TT,'in_PreS_hands_perspace_Ddsa_han');
% idx.preR_handperspace=DAG_find_column_index(TT,'in_PreR_hands_perspace_Ddre_han');
% if any(idx.preS_handperspace) && any(idx.preR_handperspace)
%     R_CH=strcmp(TT(:,idx.preR_handperspace),'CH');
%     R_IH=strcmp(TT(:,idx.preR_handperspace),'IH');
%     S_CH=strcmp(TT(:,idx.preS_handperspace),'CH');
%     S_IH=strcmp(TT(:,idx.preS_handperspace),'IH');
%     incongruent_space=  (R_CH & S_IH) | (R_IH & S_CH);
%     CH= (R_CH | S_CH) & ~incongruent_space;
%     IH= (R_IH | S_IH) & ~incongruent_space;
%     n_column=n_column+1;
%     TT(CH | IH,n_column)={'tuned'};
%     %TT(IS,n_column)={'IS'};
%     TT(~(CH|IH),n_column)={'-'};
%     TT(incongruent_space,n_column)={'incongruent'};
%     TT{1,n_column}='in_Pre_hands_perspace_Ddre_or_Ddsa';
% end
%
% idx.CueS_spaceperhand=DAG_find_column_index(TT,'in_Cue_space_perhand_Ddsa_han');
% idx.CueR_spaceperhand=DAG_find_column_index(TT,'in_Cue_space_perhand_Ddre_han');
% if any(idx.preS_spaceperhand) && any(idx.preR_spaceperhand)
%     R_CS=strcmp(TT(:,idx.CueR_spaceperhand),'CS');
%     R_IS=strcmp(TT(:,idx.CueR_spaceperhand),'IS');
%     S_CS=strcmp(TT(:,idx.CueS_spaceperhand),'CS');
%     S_IS=strcmp(TT(:,idx.CueS_spaceperhand),'IS');
%     incongruent_space=  (R_CS & S_IS) | (R_IS & S_CS);
%     CS= (R_CS | S_CS) & ~incongruent_space;
%     IS= (R_IS | S_IS) & ~incongruent_space;
%     n_column=n_column+1;
%     TT(CS | IS,n_column)={'tuned'};
%     %TT(IS,n_column)={'IS'};
%     TT(~(CS|IS),n_column)={'-'};
%     TT(incongruent_space,n_column)={'incongruent'};
%     TT{1,n_column}='in_Cue_space_perhand_Ddre_or_Ddsa';
% end


combined_column=cell(size(TT,1),1);
for props=2:size(keys.tt.combine_tuning_properties,2)
    column_index=DAG_find_column_index(TT,keys.tt.combine_tuning_properties{1,props});
    if isempty(column_index)
        disp(['column not found: ' keys.tt.combine_tuning_properties{1,props}]);
        continue;
    end
    combined_column=cellfun(@(x,y) cat(2,x,num2str(y)),combined_column,TT(:,column_index),'uniformoutput',false);
end
if ~isempty(keys.tt.combine_tuning_properties)
    combined_column{1,1}=keys.tt.combine_tuning_properties{1};
else
    combined_column{1,1}='combined_effects';
end
TT(:,end+1)=combined_column;

end

function TT=cleanup_empties(TT,idx,field_name,info_class)
if ~isempty(idx.(field_name));
    temp_info=TT([false; cellfun(@(x) ~isempty(x) ,TT(2:end,idx.(field_name)))],idx.(field_name));
    if nargin <4
        info_class=class(temp_info{1});
    end
    switch info_class
        case 'double'
            TT([false; cellfun(@(x) isempty(x) ,TT(2:end,idx.(field_name)))],idx.(field_name))={NaN};
        case 'single'
            TT([false; cellfun(@(x) isempty(x) ,TT(2:end,idx.(field_name)))],idx.(field_name))={single(NaN)};
        case 'char'
            TT([false; cellfun(@(x) isempty(x) ,TT(2:end,idx.(field_name)))],idx.(field_name))={'-'};
        otherwise
            TT([false; cellfun(@(x) isempty(x) ,TT(2:end,idx.(field_name)))],idx.(field_name))={NaN};
    end
end;
end

function TT=get_VM(TT,idx1,idx2,modulation,field_name)
if any(idx1) && any(idx2)
    n_column=size(TT,2)+1;
    idx1_modulated=ismember(TT(2:end,idx1),{'en','su','bi'});
    idx2_modulated=ismember(TT(2:end,idx2),{'en','su','bi'});
    switch modulation
        case 'first'
            to_add=num2cell( idx1_modulated &~idx2_modulated);
        case 'second'
            to_add=num2cell(~idx1_modulated & idx2_modulated);
        case 'both'
            to_add=num2cell( idx1_modulated & idx2_modulated);
        case 'none'
            to_add=num2cell(~idx1_modulated &~idx2_modulated);
    end
    TT{1,n_column}=field_name;
    TT(2:end,n_column)= to_add;
end
end

function TT=get_VMI(TT,idx1,idx2,field_name,mode)
if any(idx1) && any(idx2)
    n_column=size(TT,2)+1;
    idx1_EN=cell2mat(TT(2:end,idx1));
    idx2_EN=cell2mat(TT(2:end,idx2));
    switch mode
        case 'absolute'
            to_add=num2cell((abs(idx1_EN) - abs(idx2_EN))./(abs(idx1_EN) + abs(idx2_EN)));
        case 'signed'
            to_add=num2cell((idx1_EN - idx2_EN)./(idx1_EN + cidx2_EN));
    end
    TT{1,n_column}=field_name;
    TT(2:end,n_column)= to_add;
end
end

function TT=get_across_task(TT,idx1,idx2,mode,field_name)

n_column=size(TT,2)+1;
switch mode
    case 'hands'
        entries={'CH','IH'};
    case 'space'
        entries={'CS','IS'};
end

if any(idx1) && any(idx1)
    IX1C=strcmp(TT(:,idx1),entries{1});
    IX1I=strcmp(TT(:,idx1),entries{2});
    IX2C=strcmp(TT(:,idx2),entries{1});
    IX2I=strcmp(TT(:,idx2),entries{2});
    incongruent=  (IX1C & IX2I) | (IX1I & IX2C);
    Contra= (IX1C | IX2C) & ~incongruent;
    Ipsi= (IX1I | IX2I) & ~incongruent;
    TT(Contra | Ipsi,n_column)={'tuned'};
    TT(~(Contra|Ipsi),n_column)={'-'};
    TT(incongruent,n_column)={'incongruent'};
    TT{1,n_column}=field_name;
end
end