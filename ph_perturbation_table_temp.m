function ph_perturbation_table_Copy(project,versions)
% This function aim to compare the effect of perturbation PER UNIT accross
% different task condition and epochs
keys=struct;
keys=ph_general_settings(project,keys);

keys.tt.IC_for_criterion            = 'in';
%github test
%test bruno

project_specific_settings=[keys.db_folder filesep project filesep 'ph_project_settings.m'];
run(project_specific_settings);

keys.project_versions=versions;
subplot_rows = 2;

for f=1:numel(keys.project_versions)
    
   keys.project_versions=versions;
%    if ~isempty(keys.project_versions{f})
%         keys.project_version=keys.project_versions{f};
%    end
    version_folder=char(keys.project_versions);
    keys.version_specific_settings=[keys.db_folder project filesep version_folder filesep 'ph_project_version_settings.m'];
    run(keys.version_specific_settings);
    keys.project_version=version_folder;
    
    
    keys.tt.choices             = 0;
    keys.tt.hands               = [1 2];
    keys.tt.perturbations       = [0 1];
    keys.tt.tasktypes           = {'Ddre_han'};
    
    for m=1:numel(keys.batching.monkeys)
        keys.monkey=keys.batching.monkeys{m};
        for t=1:numel(keys.batching.targets)
            target=keys.batching.targets{t};
            
            keys.tt.selection                   ={'target',target};
            tasktype=keys.tt.tasktypes{1};
            keys.anova_table_file= ['Y:\Projects\' project '\ephys\' version_folder '\tuning_table_combined_CI.mat'];
            [tuning_per_unit_table]                 = ph_load_extended_tuning_table(keys);
            [tuning_per_unit_table, Sel_for_title]  = ph_reduce_tuning_table(tuning_per_unit_table,keys);
            
            
            for kk=1:size(tuning_per_unit_table,2)
                idx.(tuning_per_unit_table{1,kk})=kk;
            end
            
            
            if any(strfind(keys.project_version,'MIP'))
                epochs={'Facq','Fhol','Cue','Del','PreR','PeriR'};
                
            elseif  any(strfind(keys.project_version,'LIP'))
                epochs={'Facq','Fhol','Cue','Del','PreS','PeriS'};
            end
            
            figure
            for e=1:numel(epochs)
                
                epoch=epochs{e};
                
                subplot(subplot_rows,numel(epochs)/2,e);
               
                HSC = [1 2 3 4 5]; %Four Hand-Space conditons +1
                
                for unit = 1:(size(tuning_per_unit_table,1)-1)
                    
                   if any(strfind(target,'_L'))
                        
                        eff_IH_IS{unit} = [tuning_per_unit_table{unit+1, idx.(['in_CH_CS_' epoch '_PT_' tasktype])}];
                        eff_IH_CS{unit} = [tuning_per_unit_table{unit+1, idx.(['in_CH_IS_' epoch '_PT_' tasktype])}];
                        eff_CH_IS{unit} = [tuning_per_unit_table{unit+1, idx.(['in_IH_CS_' epoch '_PT_' tasktype])}];
                        eff_CH_CS{unit} = [tuning_per_unit_table{unit+1, idx.(['in_IH_IS_' epoch '_PT_' tasktype])}];
                        
                        
                       eff_AH_AS_L=vertcat(eff_IH_IS,eff_IH_CS,eff_CH_IS,eff_CH_CS)';
                       eff_L= strrep(eff_AH_AS_L,{'EN'}, '1');
                       eff_L= strrep(eff_L,{'SU'}, '2');
                       eff_L= strrep(eff_L,{'-'}, '0');
                       %add an extra row and column of zeros so we Map to units,ie. 5x(no. of units + 1)
                       Df_L= cellfun(@str2double,[eff_L repmat({'0'},size(eff_L,1),1)]);
                       Df_L=vertcat(Df_L, zeros(1,size(eff_L,2)+1));
                       units_L=[1:size(vertcat(tuning_per_unit_table(2:size(tuning_per_unit_table,1),1) ,'-'))];
%                        
                     

                   else
                        
                        eff_IH_IS{unit} = [tuning_per_unit_table{unit+1, idx.(['in_IH_IS_' epoch '_PT_' tasktype])}];
                        eff_IH_CS{unit} = [tuning_per_unit_table{unit+1, idx.(['in_IH_CS_' epoch '_PT_' tasktype])}];
                        eff_CH_IS{unit} = [tuning_per_unit_table{unit+1, idx.(['in_CH_IS_' epoch '_PT_' tasktype])}];
                        eff_CH_CS{unit} = [tuning_per_unit_table{unit+1, idx.(['in_CH_CS_' epoch '_PT_' tasktype])}];
                        
                        
                       eff_AH_AS_R=vertcat(eff_IH_IS,eff_IH_CS,eff_CH_IS,eff_CH_CS)';
                       eff_R= strrep(eff_AH_AS_R,{'EN'}, '1');
                       eff_R= strrep(eff_R,{'SU'}, '2');
                       eff_R= strrep(eff_R,{'-'}, '0');
                     %add an extra row and column of zeros so we Map to units,ie. 5x(no. of units + 1)
                       Df_R= cellfun(@str2double,[eff_R repmat({'0'},size(eff_R,1),1)]);
                       Df_R=vertcat(Df_R, zeros(1,size(eff_R,2)+1));
                       units_R=[1:size(vertcat(tuning_per_unit_table(2:size(tuning_per_unit_table,1),1) ,'-'))];
                         
                        
                   end
                   

% %              
    
                end
                
                unit_IDs=tuning_per_unit_table(2:end, idx.unit_ID);
                if any(strfind(target,'_L'))
                pcolor(HSC,units_L,Df_L);
                else
                pcolor(HSC,units_R,Df_R);
                end
                x=HSC;
                set(gca,'XTick',[1.5,2.5,3.5,4.5])
                set(gca,'XTickLabel',{'IH-IS','IH-CS','CH-IS','CH-CS'})
                map = [1 1 1;  1 0 0; 0 0 0];%colors
                colormap(map)   
                title(epoch);
                 
            end

        end
    end
end
            
% %                    xlabel('Hand Space Conditions')
% % 
% %                    ylabel('Units')  

disp('done');


end

