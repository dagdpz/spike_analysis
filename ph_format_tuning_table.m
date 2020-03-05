function ph_format_tuning_table(tuning_per_unit_table,keys)
keys.tuning_table_foldername    =[keys.basepath_to_save keys.project_version];
tuning_per_unit_table           =ph_add_Subregions_to_table(tuning_per_unit_table,keys.batching.Subregions);
tuning_per_unit_table           =convert_LR_to_CI(tuning_per_unit_table);
save([keys.tuning_table_foldername filesep 'tuning_table_combined_CI'],'tuning_per_unit_table')

tasktypes={};
tasktypes_index=1;
for effector=keys.cal.effectors
    for type=keys.cal.types
        [~, tasktypes{tasktypes_index}]=get_type_effector_name(type,effector);
        tasktypes_index=tasktypes_index+1;
    end
end
excel_table=format_excel_tuning_table(tuning_per_unit_table,tasktypes,keys.position_and_plotting_arrangements);
xlswrite([keys.tuning_table_foldername filesep keys.tuning_table_filename],excel_table);
end

function tuning_per_unit_table=convert_LR_to_CI(tuning_per_unit_table)
idx_target=DAG_find_column_index(tuning_per_unit_table,'target');
R_hem=cellfun(@(x) strcmpi(x(end-1:end),'_R'),tuning_per_unit_table(:,idx_target));
L_hem=cellfun(@(x) strcmpi(x(end-1:end),'_L'),tuning_per_unit_table(:,idx_target));

%% Replacing left/right hemifield and hand preference by ipsi/contra
O_entries={'LS','RS','LH','RH'};
R_entries={'CS','IS','CH','IH'};
L_entries={'IS','CS','IH','CH'};
tuning_table_temp=tuning_per_unit_table(R_hem,:);
for e=1:numel(O_entries)
    tuning_table_temp(cellfun(@(x) strcmp(x,O_entries{e}),tuning_table_temp))=R_entries(e);
end
tuning_per_unit_table(R_hem,:)=tuning_table_temp;

tuning_table_temp=tuning_per_unit_table(L_hem,:);
for e=1:numel(O_entries)
    tuning_table_temp(cellfun(@(x) strcmp(x,O_entries{e}),tuning_table_temp))=L_entries(e);
end
tuning_per_unit_table(L_hem,:)=tuning_table_temp;

%% inverting effect sizes (ES) and index (IX) for right hemisphere, so R-L becomes C-I
ES_columns=cellfun(@(x) any(strfind(x,'_ES_')),tuning_per_unit_table(1,:));
tuning_per_unit_table(R_hem,ES_columns)= cellfun(@(x) x*-1,tuning_per_unit_table(R_hem,ES_columns),'uniformoutput',false);

DF_columns=cellfun(@(x) any(strfind(x,'_spaceLR_DF')) || any(strfind(x,'_hands_DF')),tuning_per_unit_table(1,:));
tuning_per_unit_table(R_hem,DF_columns)= cellfun(@(x) x*-1,tuning_per_unit_table(R_hem,DF_columns),'uniformoutput',false);

IX_columns=cellfun(@(x) any(strfind(x,'_spaceLR_IX')) || any(strfind(x,'_hands_IX')),tuning_per_unit_table(1,:));
tuning_per_unit_table(R_hem,IX_columns)= cellfun(@(x) x*-1,tuning_per_unit_table(R_hem,IX_columns),'uniformoutput',false);

%% completing tuning table columns to have both hands for each respective parameter
LH_columns=find(cellfun(@(x) any(strfind(x,'_LH_')),tuning_per_unit_table(1,:)));
RH_columns=find(cellfun(@(x) any(strfind(x,'_RH_')),tuning_per_unit_table(1,:)));
LH_column_names_as_RH=cellfun(@(x) strrep(x,'_LH_','_RH_'),tuning_per_unit_table(1,LH_columns),'uniformoutput',false);
RH_column_names_as_RH=tuning_per_unit_table(1,RH_columns);
LH_missing=RH_column_names_as_RH(~ismember(RH_column_names_as_RH,LH_column_names_as_RH));
RH_missing=LH_column_names_as_RH(~ismember(LH_column_names_as_RH,RH_column_names_as_RH));
LH_missing=cellfun(@(x) strrep(x,'_RH_','_LH_'),LH_missing,'uniformoutput',false);
if ~isempty(RH_missing)
    tuning_per_unit_table(1,end+1:end+numel(RH_missing))=RH_missing;
end
if ~isempty(LH_missing)
    tuning_per_unit_table(1,end+1:end+numel(LH_missing))=LH_missing;
end

%% replacing LH and RH with IH and CH
LH_columns=find(cellfun(@(x) any(strfind(x,'_LH_')),tuning_per_unit_table(1,:)));
RH_columns=find(cellfun(@(x) any(strfind(x,'_RH_')),tuning_per_unit_table(1,:)));
Lhem_LH=tuning_per_unit_table(L_hem,LH_columns);
Lhem_RH=tuning_per_unit_table(L_hem,RH_columns);
for c=LH_columns
    table_title=tuning_per_unit_table{1,c};
    position_in_title_to_change=strfind(table_title,'_LH_')+1;
    table_title(position_in_title_to_change)='R';
    counterhand_column=DAG_find_column_index(tuning_per_unit_table,table_title);
    tuning_per_unit_table(L_hem,c)=Lhem_RH(:,RH_columns==counterhand_column);
    tuning_per_unit_table{1,c}(position_in_title_to_change)='C';
end
for c=RH_columns
    table_title=tuning_per_unit_table{1,c};
    position_in_title_to_change=strfind(table_title,'_RH_')+1;
    table_title(position_in_title_to_change)='C';
    counterhand_column=DAG_find_column_index(tuning_per_unit_table,table_title);
    tuning_per_unit_table(L_hem,c)=Lhem_LH(:,LH_columns==counterhand_column);
    tuning_per_unit_table{1,c}(position_in_title_to_change)='I';
end

%% replacing LS and RS with IS and CS
LS_columns=find(cellfun(@(x) any(strfind(x,'_LS_')),tuning_per_unit_table(1,:)));
RS_columns=find(cellfun(@(x) any(strfind(x,'_RS_')),tuning_per_unit_table(1,:)));
Lhem_LH=tuning_per_unit_table(L_hem,LS_columns);
Lhem_RH=tuning_per_unit_table(L_hem,RS_columns);
for c=LS_columns
    table_title=tuning_per_unit_table{1,c};
    position_in_title_to_change=strfind(table_title,'_LS_')+1;
    table_title(position_in_title_to_change)='R';
    counterhand_column=DAG_find_column_index(tuning_per_unit_table,table_title);
    tuning_per_unit_table(L_hem,c)=Lhem_RH(:,RS_columns==counterhand_column);
    tuning_per_unit_table{1,c}(position_in_title_to_change)='C';
end
for c=RS_columns
    table_title=tuning_per_unit_table{1,c};
    position_in_title_to_change=strfind(table_title,'_RS_')+1;
    table_title(position_in_title_to_change)='C';
    counterhand_column=DAG_find_column_index(tuning_per_unit_table,table_title);
    tuning_per_unit_table(L_hem,c)=Lhem_LH(:,LS_columns==counterhand_column);
    tuning_per_unit_table{1,c}(position_in_title_to_change)='I';
end
end

function tuning_per_unit_table=ph_add_Subregions_to_table(tuning_per_unit_table,Subregions)
idx_ID          =DAG_find_column_index(tuning_per_unit_table,'unit_ID');
idx_target      =DAG_find_column_index(tuning_per_unit_table,'target');
idx_x           =DAG_find_column_index(tuning_per_unit_table,'grid_x');
idx_y           =DAG_find_column_index(tuning_per_unit_table,'grid_y');
idx_z           =DAG_find_column_index(tuning_per_unit_table,'electrode_depth');
idx_subregion   =DAG_find_column_index(tuning_per_unit_table,'Subregion');
if isempty(idx_subregion)
    idx_subregion   =size(tuning_per_unit_table,2)+1;
end
tuning_per_unit_table{1,idx_subregion}='Subregion';
for r=1:numel(Subregions)
    for h=1:numel(Subregions{r})
        i1=~cellfun(@isempty,strfind(tuning_per_unit_table(2:end,idx_ID),Subregions{r}{h}.monkey));
        i2=~cellfun(@isempty,strfind(tuning_per_unit_table(2:end,idx_target),Subregions{r}{h}.target));
        i3=vertcat(tuning_per_unit_table{2:end,idx_x})==Subregions{r}{h}.grid_x;
        i4=vertcat(tuning_per_unit_table{2:end,idx_y})==Subregions{r}{h}.grid_y;
        i5=vertcat(tuning_per_unit_table{2:end,idx_z})>Subregions{r}{h}.z_min & vertcat(tuning_per_unit_table{2:end,idx_z})<=Subregions{r}{h}.z_max;
        index=[false; i1&i2&i3&i4&i5];
        tuning_per_unit_table(index,idx_subregion)={r};
    end
end
end

function completed_table=format_excel_tuning_table(tuning_per_unit_table,tasks,cases)
N_columns_unchanged=11;
in_or_ch={'in','ch'};
hands={'AH','IH','CH'}; %%!!!
sides={'IS','CS'}; %%!!!
general_factors_per_hand={'position_main'};
%factors_per_hand={'epoch','epoch_ES','epoch_IX','position','position_ES','position_IX'}; %position independently for both hands somehow
factors_per_hand={'epoch','epoch_ES','epoch_IX','position','position_PV','fixation','fixation_PV','fixation','fixation_PV','PxF','PxF_PV','RF_shift','gaze_modulation_x','gaze_modulation_y','RF_choice1','RF_choice2'}; %position independently for both hands somehow
general_factors={'epoch_main','spaceLR_main','hands_main','ExS','ExH','SxH','ExSxH'};
factors={'epoch','spaceLR','spaceLR_ES','spaceLR_IX','hands','hands_ES','hands_IX','SxH','SxH_ES','SxH_IX'}; %position independently for both hands somehow
general_factors_per_hand_space={'PT_main','ExP','epoch_main'}; %position independently for both hands somehow

idx_subregion   =DAG_find_column_index(tuning_per_unit_table,'Subregion');
temp_table=[tuning_per_unit_table(:,1:N_columns_unchanged) tuning_per_unit_table(:,idx_subregion)];
        
n_table=0;
for t=1:numel(tasks)
    current_task=tasks{t};
    temp_table2=[temp_table vertcat({'task'},repmat({current_task},size(tuning_per_unit_table,1)-1,1))];
    for c=1:numel(cases)
        current_case=cases{c}(1:3);
        idx_N_trials=~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),['_trials_per_condition_' current_task '_' current_case]));
        underscore_idx=cellfun(@(x) strfind(x,'_'),tuning_per_unit_table(1,idx_N_trials),'UniformOutput',false);
        N_trial_titles=cellfun(@(x,y) ['N_trials_' x(1:y(end-4)-1)],tuning_per_unit_table(1,idx_N_trials),underscore_idx,'UniformOutput',false);
        temp_table3=[temp_table2 vertcat({'case'},repmat({current_case},size(tuning_per_unit_table,1)-1,1)) vertcat(N_trial_titles,tuning_per_unit_table(2:end,idx_N_trials))];
        
        unique_epoch_title_indxes=~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),'in_'))    & ~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),'_epoch_')) &...
            ~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),'_AH_')) &...
            cellfun(@isempty,strfind(tuning_per_unit_table(1,:),'_main_')) & ~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),tasks{t})) &...
            cellfun(@isempty,strfind(tuning_per_unit_table(1,:),'_DF_'))   & cellfun(@isempty,strfind(tuning_per_unit_table(1,:),'_IX_')) &...
            cellfun(@isempty,strfind(tuning_per_unit_table(1,:),'_EN_'))   & cellfun(@isempty,strfind(tuning_per_unit_table(1,:),'_SC_')) &...
            cellfun(@isempty,strfind(tuning_per_unit_table(1,:),'_CS_'))   & cellfun(@isempty,strfind(tuning_per_unit_table(1,:),'_IS_')) &...
            cellfun(@isempty,strfind(tuning_per_unit_table(1,:),'_bilateral')) &...
            ~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),current_case));
        epochs_tmp=tuning_per_unit_table(1,unique_epoch_title_indxes);
        epoch_string_end=cell2mat(strfind(epochs_tmp,'_epoch_'))-1;
        
        for ic=1:numel(in_or_ch)
            current_INORCH=in_or_ch{ic};
            temp_table4=[temp_table3 vertcat({'in_ch'},repmat({current_INORCH},size(tuning_per_unit_table,1)-1,1))];
            
            temp_table5=temp_table4;
            for m=1:numel(general_factors)
                current_factor=general_factors{m};
                idx=DAG_find_column_index(tuning_per_unit_table,[current_INORCH '_' current_factor '_' current_task '_' current_case]);
                if isempty(idx); continue; end;
                temp_table5=[temp_table5 vertcat({[current_factor '_general']},tuning_per_unit_table(2:end,idx))];
            end
            
            temp_table6=temp_table5;
            % per hand factors
            for h=1:numel(hands)
                for m=1:numel(general_factors_per_hand)
                    current_factor=[hands{h} general_factors_per_hand{m}];
                    idx=DAG_find_column_index(tuning_per_unit_table,[current_INORCH '_' current_factor '_' current_task '_' current_case]);
                    if isempty(idx); continue; end;
                    temp_table6=[temp_table6 vertcat({[current_factor '_general']},tuning_per_unit_table(2:end,idx))];
                end
                for s=1:numel(sides)
                    for f=1:numel(general_factors_per_hand_space)
                        idx=DAG_find_column_index(tuning_per_unit_table,[current_INORCH '_' hands{h} '_' sides{s} '_' general_factors_per_hand_space{f} '_' current_task '_' current_case]);
                        if isempty(idx); continue; end;
                        temp_table6=[temp_table6 vertcat({[hands{h} '_' sides{s} '_' general_factors_per_hand_space{f} '_tuning']},tuning_per_unit_table(2:end,idx))];
                    end
                end
            end
            
            for e=1:numel(epochs_tmp)
                n_table=n_table+1;
                epoch=epochs_tmp{e}(7:epoch_string_end(e));
                temp_table7=[temp_table6 vertcat({'epoch'},repmat({epoch},size(tuning_per_unit_table,1)-1,1))];
                for f=1:numel(factors)
                    current_factor=factors{f};
                    idx=DAG_find_column_index(tuning_per_unit_table,[current_INORCH '_' epoch '_' current_factor '_' current_task '_' current_case]);
                    if isempty(idx); continue; end;
                    temp_table7=[temp_table7 vertcat({[current_factor  '_tuning']},tuning_per_unit_table(2:end,idx))];
                end
                for h=1:numel(hands)
                    for f=1:numel(factors_per_hand)
                        idx=DAG_find_column_index(tuning_per_unit_table,[current_INORCH '_' hands{h} '_' epoch '_' factors_per_hand{f} '_' current_task '_' current_case]);
                        if isempty(idx); continue; end;
                        temp_table7=[temp_table7 vertcat({[hands{h} '_' factors_per_hand{f} '_tuning']},tuning_per_unit_table(2:end,idx))];
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
