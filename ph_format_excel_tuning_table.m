function completed_table=ph_format_excel_tuning_table(TT,tasks,keys)
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
        underscore_idx=cellfun(@(x) strfind(x,'trials'),TT(1,idx_N_trials),'UniformOutput',false);
        N_trial_titles=cellfun(@(x,y) ['N_trials_' x(1:y-2)],TT(1,idx_N_trials),underscore_idx,'UniformOutput',false);
        temp_table_IN=vertcat(N_trial_titles,TT(2:end,idx_N_trials));
        idx_N_trials=~cellfun(@isempty,strfind(TT(1,:),['trials_' keys.tt.trial_criterion_ch '_' current_task '_' current_case])) & ~cellfun(@isempty,strfind(TT(1,:),'ch_'));
        underscore_idx=cellfun(@(x) strfind(x,'trials'),TT(1,idx_N_trials),'UniformOutput',false);
        N_trial_titles=cellfun(@(x,y) ['N_trials_' x(1:y-2)],TT(1,idx_N_trials),underscore_idx,'UniformOutput',false);
        temp_table_CH=vertcat(N_trial_titles,TT(2:end,idx_N_trials));
        temp_table3=[temp_table2 vertcat({'case'},repmat({current_case},size(TT,1)-1,1)) temp_table_IN temp_table_CH];
        unique_epoch_title_indxes=(~cellfun(@isempty,strfind(TT(1,:),'in_')) | ~cellfun(@isempty,strfind(TT(1,:),'ch_'))) &...
            ~cellfun(@isempty,strfind(TT(1,:),'_epoch_')) &...
            cellfun(@isempty,strfind(TT(1,:),'_AH_')) &... %% not there any more!
            cellfun(@isempty,strfind(TT(1,:),'_IH_')) & cellfun(@isempty,strfind(TT(1,:),'_CH_')) &... %% not there any more!
            cellfun(@isempty,strfind(TT(1,:),'_main_')) & ~cellfun(@isempty,strfind(TT(1,:),tasks{t})) &...
            cellfun(@isempty,strfind(TT(1,:),'_DF_'))   & cellfun(@isempty,strfind(TT(1,:),'_IX_')) &...
            cellfun(@isempty,strfind(TT(1,:),'_EN_'))   & cellfun(@isempty,strfind(TT(1,:),'_SC_')) &...
            cellfun(@isempty,strfind(TT(1,:),'_CS_'))   & cellfun(@isempty,strfind(TT(1,:),'_IS_')) &...
            cellfun(@isempty,strfind(TT(1,:),'_bilateral')) &...
            ~cellfun(@isempty,strfind(TT(1,:),current_case));
        epochs=TT(1,unique_epoch_title_indxes);
        epochs=cellfun(@(x) x(4:strfind(x,'_epoch_')-1),epochs,'UniformOutput',false);
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
                    current_factor=[hands{h} '_' keys.AN.general_factors_per_hand{m}];
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