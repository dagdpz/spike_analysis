function ph_append_to_anova_table(keys,analysis_type,varargin)
switch analysis_type
    case 'regression'
                unit_IDs=varargin{1};
        U=varargin{2};
        clear table;
        table{1,1}='unit_ID';
        for u=1:numel(U)
            W=U(u).window;
            table{u+1,1}=unit_IDs{u};
            for w=1:numel(W)
                C=W(w);
                wn=C.name{:};
                if u==1
                    table(1,w+1:w+9)={[wn '_b1'],[wn '_b2'],[wn '_p1'],[wn '_p2'],[wn '_in1'],[wn '_in2'],[wn '_SS'],[wn '_corrP'],[wn '_corrR']};
                end
                table(u+1,w+1:w+9)={C.betas(1),C.betas(2),C.pval(1),C.pval(2),C.in(1),C.in(2),C.SS,C.corrP(1,2),C.corrR(1,2)}; %% SE(!)???
            end
        end
end

load([keys.basepath_to_save keys.project_version filesep 'tuning_table_combined_CI.mat']);
table_for_updating=table(1,:);
rows_to_update=[];
for u=1:numel(unit_IDs)
    ix=find(ismember(tuning_per_unit_table(:,1),unit_IDs{u}));
    table_for_updating(ix,:)=table(u+1,:);
    rows_to_update=[rows_to_update ix];  
end
tuning_per_unit_table = DAG_update_mastertable_cell(tuning_per_unit_table,table_for_updating,rows_to_update);
save([keys.basepath_to_save keys.project_version filesep 'tuning_table_combined_CI.mat'],'tuning_per_unit_table');

end