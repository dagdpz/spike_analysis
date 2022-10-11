function ph_format_tuning_table(tuning_per_unit_table,keys)
keys.tuning_table_foldername    =[keys.basepath_to_save keys.project_version];
tuning_per_unit_table           =ph_add_Subregions_to_table(tuning_per_unit_table,keys.batching.Subregions);
tuning_per_unit_table           =convert_LR_to_CI(tuning_per_unit_table,keys);
save([keys.tuning_table_foldername filesep 'tuning_table_combined_CI'],'tuning_per_unit_table')

tasktypes={};
tasktypes_index=1;
for effector=keys.cal.effectors
    for type=keys.cal.types
        [~, tasktypes{tasktypes_index}]=MPA_get_type_effector_name(type,effector);
        tasktypes_index=tasktypes_index+1;
    end
end
excel_table=format_excel_tuning_table(tuning_per_unit_table,tasktypes,keys);
xlswrite([keys.tuning_table_foldername filesep keys.tuning_table_filename],excel_table);
end

function TT=convert_LR_to_CI(TT,keys)
idx_target=DAG_find_column_index(TT,keys.contra_ipsi_relative_to);
R_hem=cellfun(@(x) strcmpi(x(end-1:end),'_R'),TT(:,idx_target)); % row index for right hemisphere
L_hem=cellfun(@(x) strcmpi(x(end-1:end),'_L'),TT(:,idx_target)); % row index for left hemisphere

%% Replacing left/right hemifield and hand preference by ipsi/contra
O_entries={'LS','RS','LH','RH'}; % Original entries
R_entries={'CS','IS','CH','IH'}; % What they become for right hemisphere
L_entries={'IS','CS','IH','CH'}; % What they become for left hemisphere
tuning_table_temp=TT(R_hem,:);   % reduced table with only units form right hemisphere
for e=1:numel(O_entries)         % simply finding each table entry that refers to space (LS, RS, LH, RH), and replace it
    tuning_table_temp(cellfun(@(x) strcmp(x,O_entries{e}),tuning_table_temp))=R_entries(e);
end
TT(R_hem,:)=tuning_table_temp;

tuning_table_temp=TT(L_hem,:);
for e=1:numel(O_entries)
    tuning_table_temp(cellfun(@(x) strcmp(x,O_entries{e}),tuning_table_temp))=L_entries(e);
end
TT(L_hem,:)=tuning_table_temp;

%% inverting hemifield/hand differences (DF) and normalized indexes (IX) for right hemisphere, so R-L becomes C-I (ONLY INVERTING RIGHT HEMISPHERE HERE)
DF_columns=cellfun(@(x) any(strfind(x,'_hemifield_DF')) || any(strfind(x,'_hands_DF')),TT(1,:));
TT(R_hem,DF_columns)= cellfun(@(x) x*-1,TT(R_hem,DF_columns),'uniformoutput',false);

IX_columns=cellfun(@(x) any(strfind(x,'_hemifield_IX')) || any(strfind(x,'_hands_IX')),TT(1,:));
TT(R_hem,IX_columns)= cellfun(@(x) x*-1,TT(R_hem,IX_columns),'uniformoutput',false);

IX_columns=cellfun(@(x) any(strfind(x,'_hemifield_PV')) || any(strfind(x,'_hands_PV')),TT(1,:));
TT(R_hem,IX_columns)= cellfun(@(x) x*-1,TT(R_hem,IX_columns),'uniformoutput',false);

%% completing tuning table columns to have both hands for each respective parameter
TT=create_missing_mirrored_columns(TT,'LH_','RH_');
%% but we also need to complete for space now :)
TT=create_missing_mirrored_columns(TT,'LS_','RS_');

%% This next part is a bit confusing. We replaced all the entries so far ...
% BUT the column headers are still in LH/RH LS/RS notation
% and the tricky part here is, that we have to shuffle these columns,
% because for 2 units from a different hemisphere, f.e. the value for CH
% comes from a different original column (LH or RH)
% The fact that we have both hand and space to modify doesnt
% make it more complicated, we can address them independently in two steps (which
% are pretty much the same)
% We are only swapping columns for one hemisphere (the left, but this is irrelevant, as long as column headers are named accordingly) !!

%% replacing LH and RH with IH and CH
LH_columns=find_string_in_the_beginning_or_separated_by_underscore(TT,'LH_');
RH_columns=find_string_in_the_beginning_or_separated_by_underscore(TT,'RH_');
Lhem_LH=TT(L_hem,LH_columns);                                    % reduced table Left hemisphere left  hand (temporary stored information)
Lhem_RH=TT(L_hem,RH_columns);                                    % reduced table Left hemisphere right hand (temporary stored information)
for c=LH_columns                                                 % looping through all LH columns to replace left hemisphere rows with respective RH entries
    table_title=TT{1,c};
    position_in_title_to_change=strfind(table_title,'LH_');
    table_title(position_in_title_to_change)='R';                % we replaced LH with RH, so this is now the header of the column we want to replace left hemisphere entries
    counterhand_column=DAG_find_column_index(TT,table_title);    % position of the respective RH column in the full table
    TT(L_hem,c)=Lhem_RH(:,RH_columns==counterhand_column);       % here we replace all left hemisphere rows for the current (LH) column with the respective temporarily stored RH column entries
    TT{1,c}(position_in_title_to_change)='C';                    % now we can rename the column, since we swapped RH entries to LH for left hemisphere, it's now contra hand (for both hemispheres!)
end
for c=RH_columns                                                 % this is the same as above, but now for RH columns
    table_title=TT{1,c};
    position_in_title_to_change=strfind(table_title,'RH_');
    table_title(position_in_title_to_change)='C';                % here's the only difference: there is no LH in the column headers any more, becuase we replaced it with CH
    counterhand_column=DAG_find_column_index(TT,table_title);
    TT(L_hem,c)=Lhem_LH(:,LH_columns==counterhand_column);
    TT{1,c}(position_in_title_to_change)='I';
end

%% replacing LS and RS with IS and CS - now the way tuning table is created now, we are not avoiding missing equivalents any more necessarily
LS_columns=find_string_in_the_beginning_or_separated_by_underscore(TT,'LS_');
RS_columns=find_string_in_the_beginning_or_separated_by_underscore(TT,'RS_');
Lhem_LS=TT(L_hem,LS_columns);
Lhem_RS=TT(L_hem,RS_columns);
for c=LS_columns
    table_title=TT{1,c};
    position_in_title_to_change=strfind(table_title,'LS_');
    table_title(position_in_title_to_change)='R';
    counterhand_column=DAG_find_column_index(TT,table_title);
    TT(L_hem,c)=Lhem_RS(:,RS_columns==counterhand_column);
    TT{1,c}(position_in_title_to_change)='C';
end
for c=RS_columns
    table_title=TT{1,c};
    position_in_title_to_change=strfind(table_title,'RS_');
    table_title(position_in_title_to_change)='C';
    counterhand_column=DAG_find_column_index(TT,table_title);
    TT(L_hem,c)=Lhem_LS(:,LS_columns==counterhand_column);
    TT{1,c}(position_in_title_to_change)='I';
end

%% find columns that have empties and replace those empties with appropriate entries
temp_info=[false(1,size(TT,2)); cellfun(@(x) isempty(x) ,TT(2:end,:))];
columns_to_cleanup=find(any(temp_info,1));
for c=columns_to_cleanup
    temp_info=TT([false; cellfun(@(x) ~isempty(x) ,TT(2:end,c))],c);
    if ~isempty(temp_info)
        info_class=class(temp_info{1});
        TT=cleanup_empties(TT,c,info_class);
    else
        disp(['all entries empty for ' TT{1,c}]);
    end
end
end

function idx=find_string_in_the_beginning_or_separated_by_underscore(TT,str1)
idx=find(cellfun(@(x) any(strfind(x,['_' str1])) || (any(strfind(x,str1)) && strfind(x,str1)==1),TT(1,:)));
end

function TT=cleanup_empties(TT,idx,info_class)
switch info_class
    case 'double'
        TT([false; cellfun(@(x) isempty(x) ,TT(2:end,idx))],idx)={NaN};
    case 'single'
        TT([false; cellfun(@(x) isempty(x) ,TT(2:end,idx))],idx)={single(NaN)};
    case 'char'
        TT([false; cellfun(@(x) isempty(x) ,TT(2:end,idx))],idx)={'-'};
    case 'logical'
        TT([false; cellfun(@(x) isempty(x) ,TT(2:end,idx))],idx)={false};
    otherwise
        disp('unknown table entry format, going for double NaN as default')
        TT([false; cellfun(@(x) isempty(x) ,TT(2:end,idx))],idx)={NaN};
end
end

function TT=create_missing_mirrored_columns(TT,str1,str2)
%  entries will be empty, but at least the column names are there - we need this for the next step
idx_ips=find_string_in_the_beginning_or_separated_by_underscore(TT,str1);
idx_con=find_string_in_the_beginning_or_separated_by_underscore(TT,str2);
labels_con=TT(1,idx_con);
labels_ips=TT(1,idx_ips);
labels_ips_as_con=cellfun(@(x) strrep(x,str1,str2),labels_ips,'uniformoutput',false);
labels_con_as_ips=cellfun(@(x) strrep(x,str2,str1),labels_con,'uniformoutput',false);
idx_ips_missing=~ismember(labels_con,labels_ips_as_con);
idx_con_missing=~ismember(labels_ips_as_con,labels_con);

labels_ips_missing=labels_con_as_ips(idx_ips_missing);
labels_con_missing=labels_ips_as_con(idx_con_missing);

% now create 2 tables for replacement that have the correct format of
% emptyness (single/double NaN, string);
TT_ips_to_con=TT(2:end,idx_ips(idx_con_missing));
TT_con=[labels_con_missing; cell(size(TT_ips_to_con))];
for c=1:sum(idx_con_missing)
    temp_info=TT_ips_to_con([false; cellfun(@(x) ~isempty(x) ,TT_ips_to_con(2:end,c))],c);
    info_class=class(temp_info{1});
    TT_con=cleanup_empties(TT_con,c,info_class);
end

TT_con_to_ips=TT(2:end,idx_con(idx_ips_missing));
TT_ips=[labels_ips_missing; cell(size(TT_con_to_ips))];
for c=1:sum(idx_ips_missing)
    temp_info=TT_con_to_ips([false; cellfun(@(x) ~isempty(x) ,TT_con_to_ips(2:end,c))],c);
    info_class=class(temp_info{1});
    TT_ips=cleanup_empties(TT_ips,c,info_class);
end
TT=[TT TT_ips TT_con];
end

function TT=ph_add_Subregions_to_table(TT,Subregions)
idx_ID          =DAG_find_column_index(TT,'unit_ID');
idx_target      =DAG_find_column_index(TT,'target');
idx_x           =DAG_find_column_index(TT,'grid_x');
idx_y           =DAG_find_column_index(TT,'grid_y');
idx_z           =DAG_find_column_index(TT,'electrode_depth');
idx_subregion   =DAG_find_column_index(TT,'Subregion');
if isempty(idx_subregion)
    idx_subregion   =size(TT,2)+1;
end
TT{1,idx_subregion}='Subregion';
for r=1:numel(Subregions)
    for h=1:numel(Subregions{r})
        i1=~cellfun(@isempty,strfind(TT(2:end,idx_ID),Subregions{r}{h}.monkey));
        i2=~cellfun(@isempty,strfind(TT(2:end,idx_target),Subregions{r}{h}.target));
        i3=vertcat(TT{2:end,idx_x})==Subregions{r}{h}.grid_x;
        i4=vertcat(TT{2:end,idx_y})==Subregions{r}{h}.grid_y;
        i5=vertcat(TT{2:end,idx_z})>Subregions{r}{h}.z_min & vertcat(TT{2:end,idx_z})<=Subregions{r}{h}.z_max;
        index=[false; i1&i2&i3&i4&i5];
        TT(index,idx_subregion)={r};
    end
end
end

function completed_table=format_excel_tuning_table(TT,tasks,keys)
cases=keys.position_and_plotting_arrangements;
N_columns_unchanged=12; 
in_or_ch={'in','ch'};
hands={'AH','IH','CH'}; 
sides={'IS','CS'}; 

keys.AN.general_factors={'epoch_main','hemifield_main','hands_main','ExS','ExH','SxH','ExSxH'}; %factor for ANOVA: epoch*hand*space
keys.AN.general_factors_per_hand={'position_main'};
keys.AN.general_factors_per_hand_space={'PT_main','ExP','epoch_main'};
keys.AN.factors={'epoch','hemifield','hands','SxH','SglTar_Suc'};
%keys.AN.factors_per_hand={'epoch','position','fixation','PxF','RF_shift','gaze_modulation_x','gaze_modulation_y','RF_choice1','RF_choice2'};
keys.AN.factors_per_hand={'epoch','position','fixation','PxF','distance','angle','DxA','prefH', 'prefP','positionx','positiony','_positionxy','gaze_modulation_x','gaze_modulation_y','gaze_pref_x','gaze_pref_y'};
keys.AN.factors_per_hemifield={'Difficulty_Easy', 'Difficulty_Diff', 'SpatialComp_2HFTar', 'SpatialComp_1HFTar'};
keys.AN.factors_per_handhemifield={'epoch','prefH', 'prefP', 'SpatialComp_2HFTar', 'SpatialComp_1HFTar'};

potential_affixes={'','_DF','_IX','_FR','_PV','_NM','_SC'}; %%?????
factors=create_all_string_combinations(keys.AN.factors,potential_affixes);
factors_per_hand=create_all_string_combinations(keys.AN.factors_per_hand,potential_affixes);
factors_per_hemifield=create_all_string_combinations(keys.AN.factors_per_hemifield,potential_affixes);


idx_subregion   =DAG_find_column_index(TT,'Subregion');
temp_table=[TT(:,1:N_columns_unchanged) TT(:,idx_subregion)];

n_table=0;
for t=1:numel(tasks)
    current_task=tasks{t};
    temp_table2=[temp_table vertcat({'task'},repmat({current_task},size(TT,1)-1,1))];
    for c=1:numel(cases)
        current_case=cases{c}(1:3);
        idx_N_trials=~cellfun(@isempty,strfind(TT(1,:),['trials_' keys.tt.trial_criterion_in '_' current_task '_' current_case])) & ~cellfun(@isempty,strfind(TT(1,:),'in_'));
        underscore_idx=cellfun(@(x) strfind(x,'_'),TT(1,idx_N_trials),'UniformOutput',false);
        N_trial_titles=cellfun(@(x,y) ['N_trials_' x(1:y(end-4)-1)],TT(1,idx_N_trials),underscore_idx,'UniformOutput',false);
        temp_table_IN=vertcat(N_trial_titles,TT(2:end,idx_N_trials));
        idx_N_trials=~cellfun(@isempty,strfind(TT(1,:),['trials_' keys.tt.trial_criterion_ch '_' current_task '_' current_case])) & ~cellfun(@isempty,strfind(TT(1,:),'ch_'));
        underscore_idx=cellfun(@(x) strfind(x,'_'),TT(1,idx_N_trials),'UniformOutput',false);
        N_trial_titles=cellfun(@(x,y) ['N_trials_' x(1:y(end-4)-1)],TT(1,idx_N_trials),underscore_idx,'UniformOutput',false);
        temp_table_CH=vertcat(N_trial_titles,TT(2:end,idx_N_trials));
        temp_table3=[temp_table2 vertcat({'case'},repmat({current_case},size(TT,1)-1,1)) temp_table_IN temp_table_CH];
        unique_epoch_title_indxes=(~cellfun(@isempty,strfind(TT(1,:),'in_')) | ~cellfun(@isempty,strfind(TT(1,:),'ch_'))) &...
            ~cellfun(@isempty,strfind(TT(1,:),'_epoch_')) &...
            ~cellfun(@isempty,strfind(TT(1,:),'_AH_')) &...
            cellfun(@isempty,strfind(TT(1,:),'_main_')) & ~cellfun(@isempty,strfind(TT(1,:),tasks{t})) &...
            cellfun(@isempty,strfind(TT(1,:),'_DF_'))   & cellfun(@isempty,strfind(TT(1,:),'_IX_')) &...
            cellfun(@isempty,strfind(TT(1,:),'_EN_'))   & cellfun(@isempty,strfind(TT(1,:),'_SC_')) &...
            cellfun(@isempty,strfind(TT(1,:),'_CS_'))   & cellfun(@isempty,strfind(TT(1,:),'_IS_')) &...
            cellfun(@isempty,strfind(TT(1,:),'_bilateral')) &...
            ~cellfun(@isempty,strfind(TT(1,:),current_case));
        epochs=TT(1,unique_epoch_title_indxes);
        epochs=cellfun(@(x) x(7:strfind(x,'_epoch_')-1),epochs,'UniformOutput',false);
        epochs=unique(epochs);
        
        for ic=1:numel(in_or_ch)
            current_INORCH=in_or_ch{ic};
            temp_table4=[temp_table3 vertcat({'in_ch'},repmat({current_INORCH},size(TT,1)-1,1))];            
            temp_table5=temp_table4;
            
            % general factors
            for m=1:numel(keys.AN.general_factors)
                current_factor=keys.AN.general_factors{m};
                idx=DAG_find_column_index(TT,[current_INORCH '_' current_factor '_' current_task '_' current_case]);
                if isempty(idx); continue; end;
                temp_table5=[temp_table5 vertcat({[current_factor '_general']},TT(2:end,idx))];
            end            
            temp_table6=temp_table5;
            
            % per hand & per handspace factors
            for h=1:numel(hands)
                for m=1:numel(keys.AN.general_factors_per_hand)
                    current_factor=[hands{h} keys.AN.general_factors_per_hand{m}];
                    idx=DAG_find_column_index(TT,[current_INORCH '_' current_factor '_' current_task '_' current_case]);
                    if isempty(idx); continue; end;
                    temp_table6=[temp_table6 vertcat({[current_factor '_general']},TT(2:end,idx))];
                end
                for s=1:numel(sides)
                    for f=1:numel(keys.AN.general_factors_per_hand_space)
                        idx=DAG_find_column_index(TT,[current_INORCH '_' hands{h} '_' sides{s} '_' keys.AN.general_factors_per_hand_space{f} '_' current_task '_' current_case]);
                        if isempty(idx); continue; end;
                        temp_table6=[temp_table6 vertcat({[hands{h} '_' sides{s} '_' keys.AN.general_factors_per_hand_space{f} '_tuning']},TT(2:end,idx))];
                    end
                end
            end
            
            % per epoch factors
            for e=1:numel(epochs)
                n_table=n_table+1;
                epoch=epochs{e};
                temp_table7=[temp_table6 vertcat({'epoch'},repmat({epoch},size(TT,1)-1,1))];
                for f=1:numel(factors)
                    current_factor=factors{f};
                    idx=DAG_find_column_index(TT,[current_INORCH '_' epoch '_' current_factor '_' current_task '_' current_case]);
                    if isempty(idx); continue; end;
                    temp_table7=[temp_table7 vertcat({[current_factor  '_tuning']},TT(2:end,idx))];
                end
                for s=1:numel(sides)
                    for f=1:numel(factors_per_hemifield)
                        idx=DAG_find_column_index(TT,[current_INORCH '_' sides{s} '_' epoch '_' factors_per_hemifield{f} '_' current_task '_' current_case]);
                        if isempty(idx); continue; end;
                        temp_table7=[temp_table7 vertcat({[sides{s} '_' factors_per_hemifield{f} '_tuning']},TT(2:end,idx))];
                    end
                end
                for h=1:numel(hands)
                    for f=1:numel(factors_per_hand)
                        idx=DAG_find_column_index(TT,[current_INORCH '_' hands{h} '_' epoch '_' factors_per_hand{f} '_' current_task '_' current_case]);
                        if isempty(idx); continue; end;
                        temp_table7=[temp_table7 vertcat({[hands{h} '_' factors_per_hand{f} '_tuning']},TT(2:end,idx))];
                    end
                end
                final_cell_table{n_table}=temp_table7;
            end
        end
    end
end

completed_table=final_cell_table{1}(1,:);
for n=1:numel(final_cell_table)
    table_for_updating=final_cell_table{n};
    n_rows=size(completed_table,1);
    n_columns=size(completed_table,2);
    for c=1:size(table_for_updating,2)
        idx=DAG_find_column_index(completed_table,table_for_updating{1,c});
        if isempty(idx)
            completed_table(1,n_columns+1) = table_for_updating(1,c);
            completed_table(n_rows+1:n_rows+size(table_for_updating,1)-1,n_columns+1) = table_for_updating(2:end,c);
            n_columns=n_columns+1;
        else
            completed_table(n_rows+1:n_rows+size(table_for_updating,1)-1,idx)=table_for_updating(2:end,c);
        end
    end
end
end

function C=create_all_string_combinations(A,B)
A=A(:);
C={};
for x=1:numel(B)
    C=[C;strcat(A,B(x))];
end
C=C';
end