function format_tuning_table(folder,filename,tasks,cases)
%format_tuning_table('W:\Projects\Pulv_microstim_behavior\ephys\ephys_analysis_v7_July2016_coordinates','tuning_table_mem_combined')
%
% folder='W:\Projects\Pulv_microstim_behavior\ephys\ephys_analysis_v7_July2016_coordinates';
% filename='tuning_table_mem_combined';

load([folder filesep filename '.mat']);



sides={'LS','RS'};
% !!!!
N_columns_unchanged=8;
% tasks={'Msac'};
% cases={'opt'};

in_or_ch={'in','ch'};
hands={'NH','LH','RH'}; %%!!!
general_factors_per_hand={'position_main'};
%factors_per_hand={'epoch','epoch_ES','epoch_IX','position','position_ES','position_IX'}; %position independently for both hands somehow
factors_per_hand={'monotonous','epoch','epoch_ES','epoch_IX','position','position_PV','fixation','fixation_PV','fixation','fixation_PV','PxF','PxF_PV','RF_choice1','RF_choice2'}; %position independently for both hands somehow



general_factors={'epoch_main','spaceLR_main','hands_main','ExS','ExH','SxH','ExSxH'};
factors={'spaceLR','spaceLR_ES','spaceLR_IX','hands','hands_ES','hands_IX','SxH','SxH_ES','SxH_IX'}; %position independently for both hands somehow

%N_cells=size(tuning_per_unit_table,1)-1;
temp_table=tuning_per_unit_table(:,1:N_columns_unchanged);
idx_subregion   =find_column_index(tuning_per_unit_table,'Subregion');
temp_table=[temp_table tuning_per_unit_table(:,idx_subregion)];
n_table=0;
for t=1:numel(tasks)
    current_task=tasks{t};
    
    temp_table2=[temp_table vertcat({'task'},repmat({current_task},size(tuning_per_unit_table,1)-1,1))];
    for c=1:numel(cases)
        current_case=cases{c}(1:3);
        temp_table3=[temp_table2 vertcat({'case'},repmat({current_case},size(tuning_per_unit_table,1)-1,1))];
        
        unique_epoch_title_indxes=~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),'in_'))    & ~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),'_epoch_')) &...
            ~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),'_NH_')) &...
            cellfun(@isempty,strfind(tuning_per_unit_table(1,:),'_main_')) & ~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),tasks{t})) &...
            cellfun(@isempty,strfind(tuning_per_unit_table(1,:),'_ES_'))   & cellfun(@isempty,strfind(tuning_per_unit_table(1,:),'_IX_')) &...
            ~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),current_case));
        epochs_tmp=tuning_per_unit_table(1,unique_epoch_title_indxes);
        epoch_string_end=cell2mat(strfind(epochs_tmp,'_epoch_'))-1;
        
        for ic=1:numel(in_or_ch)
            current_INORCH=in_or_ch{ic};
            temp_table4=[temp_table3 vertcat({'in_ch'},repmat({current_INORCH},size(tuning_per_unit_table,1)-1,1))];
            
            temp_table5=temp_table4;
            for m=1:numel(general_factors)
                current_factor=general_factors{m};
                idx=find_column_index(tuning_per_unit_table,[current_INORCH '_' current_factor '_' current_task '_' current_case]);
                if isempty(idx); continue; end;
                temp_table5=[temp_table5 vertcat({[current_factor '_general']},tuning_per_unit_table(2:end,idx))];
            end
            
            temp_table51=temp_table5;
            % per hand factors
            for h=1:numel(hands)
                %             temp_table51=[temp_table51 vertcat({'in_ch'},repmat({current_INORCH},size(tuning_per_unit_table,1)-1,1))];
                %             temp_table52=temp_table51;
                for m=1:numel(general_factors_per_hand)
                    current_factor=[hands{h} general_factors_per_hand{m}];
                    idx=find_column_index(tuning_per_unit_table,[current_INORCH '_' current_factor '_' current_task '_' current_case]);
                    if isempty(idx); continue; end;
                    temp_table51=[temp_table51 vertcat({[current_factor '_general']},tuning_per_unit_table(2:end,idx))];
                end
            end
            
            for e=1:numel(epochs_tmp)
                n_table=n_table+1;
                
                epoch=epochs_tmp{e}(7:epoch_string_end(e));
                temp_table6=[temp_table51 vertcat({'epoch'},repmat({epoch},size(tuning_per_unit_table,1)-1,1))];
                for f=1:numel(factors)
                    current_factor=factors{f};
                    %                     if strcmp(current_factor,'epoch')
                    % %                         idx=find_column_index(tuning_per_unit_table,[current_INORCH '_L_' epoch '_' current_task '_' current_case]);
                    % %                         if isempty(idx); continue; end;
                    % %                         temp_table6=[temp_table6 vertcat({[current_factor '_L']},tuning_per_unit_table(2:end,idx))];
                    % %                         idx=find_column_index(tuning_per_unit_table,[current_INORCH '_R_' epoch '_' current_task '_' current_case]);
                    % %                         if isempty(idx); continue; end;
                    % %                         temp_table6=[temp_table6 vertcat({[current_factor '_R']},tuning_per_unit_table(2:end,idx))];
                    %                         idx=find_column_index(tuning_per_unit_table,[current_INORCH '_' current_factor '_' epoch '_' current_task '_' current_case]);
                    %                         if isempty(idx); continue; end;
                    %                         temp_table6=[temp_table6 vertcat({current_factor},tuning_per_unit_table(2:end,idx))];
                    %                     else
                    idx=find_column_index(tuning_per_unit_table,[current_INORCH '_' epoch '_' current_factor '_' current_task '_' current_case]);
                    if isempty(idx); continue; end;
                    temp_table6=[temp_table6 vertcat({[current_factor  '_tuning']},tuning_per_unit_table(2:end,idx))];
                    %                    end
                end
                
                for h=1:numel(hands)
                    for f=1:numel(factors_per_hand)
                        idx=find_column_index(tuning_per_unit_table,[current_INORCH '_' hands{h} '_' epoch '_' factors_per_hand{f} '_' current_task '_' current_case]);
                        if isempty(idx); continue; end;
                        temp_table6=[temp_table6 vertcat({[hands{h} '_' factors_per_hand{f} '_tuning']},tuning_per_unit_table(2:end,idx))];
                        %
                    end
                end
                
                final_cell_table{n_table}=temp_table6;
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
        idx=find_column_index(completed_table,table_for_updating{1,c});
        if isempty(idx)
            completed_table(1,n_columns+1) = table_for_updating(1,c);
            completed_table(n_rows+1:n_rows+size(table_for_updating,1)-1,n_columns+1) = table_for_updating(2:end,c);
            n_columns=n_columns+1;
        else
            completed_table(n_rows+1:n_rows+size(table_for_updating,1)-1,idx)=table_for_updating(2:end,c);
        end
    end
end


xlswrite([folder filesep filename],completed_table);
end
