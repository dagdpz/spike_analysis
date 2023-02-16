function ph_append_to_anova_table(keys,analysis_type,varargin)
        table{1,1}='unit_ID';
switch analysis_type
    case 'regression'
        unit_IDs=varargin{1};
        U=varargin{2};
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
    case 'RFs'
        U=varargin{2};
        typ_eff=varargin{3};
        label=varargin{4};
        for u=1:numel(U)   
        unit_IDs=varargin{1};       
            table{u+1,1}=unit_IDs{u};
            RFstruct=U(u).(U(u).bestfit);  
            if strcmp(U(u).bestfit,'gaussian15')
               if RFstruct.zmax1>0
               RFstruct.sx=RFstruct.sx1; 
               RFstruct.sy=RFstruct.sy1; 
               RFstruct.phi=RFstruct.phi1; 
               RFstruct.xmax=RFstruct.xmax1; 
               RFstruct.ymax=RFstruct.ymax1; 
               RFstruct.zmax=RFstruct.zmax1; 
               else
               RFstruct.sx=RFstruct.sx2; 
               RFstruct.sy=RFstruct.sy2; 
               RFstruct.phi=RFstruct.phi2; 
               RFstruct.xmax=RFstruct.xmax2; 
               RFstruct.ymax=RFstruct.ymax2; 
               RFstruct.zmax=RFstruct.zmax2; 
               end
            end
            if u==1
                table(1,2:3)={[label '_' keys.RF.epoch_RF '_center_x_' typ_eff '_' keys.arrangement(1:3)],[label '_' keys.RF.epoch_RF '_size_x_' typ_eff '_' keys.arrangement(1:3)]}; %% add RF epoch
            end
            if any(strcmp(fieldnames(RFstruct),'sx'))
                RFstruct.xsize=max(abs(RFstruct.sx*cos(RFstruct.phi)),abs(RFstruct.sy*sin(RFstruct.phi)));
                
                if RFstruct.zmax>0 %% exclude suppressed ones here
                table(u+1,2:3)={RFstruct.xmax,RFstruct.xsize};
                else
                    table(u+1,2:3)={[],[]};
                end
                
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