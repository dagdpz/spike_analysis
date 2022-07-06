function TT=ph_load_extended_tuning_table(keys)
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

TT(:,end+1)=num2cell(repmat('Y',size(TT,1),1));
TT(1,end)={'ungrouped'};

if ~all(cellfun(@isempty,strfind(TT(1,:),'_trials_per_position'))) %~isempty(taskcaseexistingindex)
    taskcases={};
    for criterium={'in','ch'}
        crit=criterium{:};
        switch keys.tt.(['trial_criterion_' crit])
            case 'per_hemifield_and_perturbation'
                strtofind='trials_per_hemifield'; %% ??? - there will be a problem with perturbation stuff !!
                disp(['per_hemifield_and_perturbation not supported as trial criterion' ]);
        end
        
        strtofind=['_trials_' keys.tt.(['trial_criterion_' crit])];
        taskcaseexistingindex=find(~cellfun(@isempty,strfind(TT(1,:),strtofind)));
        for I=taskcaseexistingindex
            if ~(any(strfind(TT{1,I},[crit '_'])==1) || any(strfind(TT{1,I},['_' crit '_'])))
                continue;
            end
            task_existing_column=num2cell(cellfun(@(x) ~isempty(x)&&~ischar(x)&&x>=keys.cal.min_trials_per_condition,TT(:,I)));
            strpos=strfind(TT{1,I},strtofind);
            task_existing_column{1,1}=['existing_' TT{1,I}([1:strpos-1,strpos+numel(strtofind)+1:end])];
            TT(:,end+1)=task_existing_column;
            taskcases = [taskcases, {TT{1,I}(end-7:end)}];
        end
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
    idx.space   =DAG_find_column_index(TT,[keys.tt.IC_for_criterion '_hemifield_main_'   taskcase] );
    idx.epoch   =DAG_find_column_index(TT,[keys.tt.IC_for_criterion '_epoch_main_'     taskcase] );
    idx.hands   =DAG_find_column_index(TT,[keys.tt.IC_for_criterion '_hands_main_'     taskcase] );
    idx.SxH     =DAG_find_column_index(TT,[keys.tt.IC_for_criterion '_SxH_'            taskcase] );
    idx.ExS     =DAG_find_column_index(TT,[keys.tt.IC_for_criterion '_ExS_'            taskcase] );
    idx.ExH     =DAG_find_column_index(TT,[keys.tt.IC_for_criterion '_ExH_'            taskcase] );
    idx.ESH     =DAG_find_column_index(TT,[keys.tt.IC_for_criterion '_ExSxH_'          taskcase] );
    
    for e=1:numel(epochs)
        idx.(epochs{e})                     =DAG_find_column_index(TT,[keys.tt.IC_for_criterion '_AH_' epochs{e} '_epoch_'   taskcase] );
        idx.([epochs{e} '_space_DF'])       =DAG_find_column_index(TT,[keys.tt.IC_for_criterion '_' epochs{e} '_hemifield_DF_'   taskcase] );
        idx.([epochs{e} '_IS_FR'])          =DAG_find_column_index(TT,[keys.tt.IC_for_criterion '_AH_IS_' epochs{e} '_epoch_FR_'  taskcase] );
        idx.([epochs{e} '_CS_FR'])          =DAG_find_column_index(TT,[keys.tt.IC_for_criterion '_AH_CS_' epochs{e} '_epoch_FR_'  taskcase] );
        idx.([epochs{e} '_IS_EN'])          =DAG_find_column_index(TT,[keys.tt.IC_for_criterion '_AH_IS_' epochs{e} '_epoch_DF_'  taskcase] );
        idx.([epochs{e} '_CS_EN'])          =DAG_find_column_index(TT,[keys.tt.IC_for_criterion '_AH_CS_' epochs{e} '_epoch_DF_'  taskcase] );
        idx.([epochs{e} '_in_space_DF'])    =DAG_find_column_index(TT,['in_' epochs{e} '_hemifield_DF_'   taskcase] );
        idx.([epochs{e} '_ch_space_DF'])    =DAG_find_column_index(TT,['ch_' epochs{e} '_hemifield_DF_'   taskcase] );
        idx.([epochs{e} '_in_IS_FR'])       =DAG_find_column_index(TT,['in_AH_IS_' epochs{e} '_epoch_FR_'  taskcase] );
        idx.([epochs{e} '_in_CS_FR'])       =DAG_find_column_index(TT,['in_AH_CS_' epochs{e} '_epoch_FR_'  taskcase] );
        idx.([epochs{e} '_ch_IS_FR'])       =DAG_find_column_index(TT,['ch_AH_IS_' epochs{e} '_epoch_FR_'  taskcase] );
        idx.([epochs{e} '_ch_CS_FR'])       =DAG_find_column_index(TT,['ch_AH_CS_' epochs{e} '_epoch_FR_'  taskcase] );
        
        idx.([epochs{e} '_in_IH_space'])       =DAG_find_column_index(TT,['in_IH_' epochs{e} '_hemifield_'  taskcase] );
        idx.([epochs{e} '_in_CH_space'])       =DAG_find_column_index(TT,['in_CH_' epochs{e} '_hemifield_'  taskcase] );
        idx.([epochs{e} '_in_IS_hands'])       =DAG_find_column_index(TT,['in_IS_' epochs{e} '_hands_'  taskcase] );
        idx.([epochs{e} '_in_CS_hands'])       =DAG_find_column_index(TT,['in_CS_' epochs{e} '_hands_'  taskcase] );
    end
    
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
    %% MP end
    
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
            to_add=num2cell((idx1_EN - idx2_EN)./(idx1_EN + idx2_EN));
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