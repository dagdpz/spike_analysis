function ph_perturbation_table(project,versions)
% This function aim to compare the effect of perturbation PER UNIT accross
% different task condition and epochs
keys=struct;
keys=ph_general_settings(project,keys);

keys.tt.IC_for_criterion            = 'in';


project_specific_settings=[keys.db_folder filesep project filesep 'ph_project_settings.m'];
run(project_specific_settings);

keys.project_versions=versions;
subplot_rows = 2;

for f=1:numel(keys.project_versions)
    
   keys.project_versions=versions;
   if ~isempty(keys.project_versions{f})
        keys.project_version=keys.project_versions{f};
   end
    version_folder=keys.project_version;
    keys.version_specific_settings=[keys.db_folder project filesep keys.project_version filesep 'ph_project_version_settings.m'];
    run(keys.version_specific_settings);
    keys.project_version=version_folder;
    
    
    keys.tt.choices             = 0;
    keys.tt.hands               = [1 2];
    keys.tt.perturbations       = [0 1];
%     keys.tt.tasktypes           = {4};
    keys.tt.tasktypes           = {'Ddre_han'};
    
    for m=1:numel(keys.batching.monkeys)
        keys.monkey=keys.batching.monkeys{m};
        for t=1:numel(keys.batching.targets)
            target=keys.batching.targets{t};
            
            keys.tt.selection                   ={'target',target};
            tasktype=keys.tt.tasktypes{1};
            keys.anova_table_file=horzcat(['Y:\Projects\' project '\ephys\' version_folder '\tuning_table_combined_CI.mat']);
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
            
            
            for e=1:numel(epochs)
                
                epoch=epochs{e};
                
                subplot(subplot_rows,numel(epochs)/2,e);
                title(epoch);
                
                for unit = 1:(size(tuning_per_unit_table,1)-1)
                    
                    if any(strfind(target,'_L'))
                        
                        eff_IH_IS{unit} = [tuning_per_unit_table{unit+1, idx.(['in_CH_CS_' epoch '_PT_' tasktype])}];
                        eff_IH_CS{unit} = [tuning_per_unit_table{unit+1, idx.(['in_CH_IS_' epoch '_PT_' tasktype])}];
                        eff_CH_IS{unit} = [tuning_per_unit_table{unit+1, idx.(['in_IH_CS_' epoch '_PT_' tasktype])}];
                        eff_CH_CS{unit} = [tuning_per_unit_table{unit+1, idx.(['in_IH_IS_' epoch '_PT_' tasktype])}];
                   
                    else
                        
                        eff_IH_IS{unit} = [tuning_per_unit_table{unit+1, idx.(['in_IH_IS_' epoch '_PT_' tasktype])}];
                        eff_IH_CS{unit} = [tuning_per_unit_table{unit+1, idx.(['in_IH_CS_' epoch '_PT_' tasktype])}];
                        eff_CH_IS{unit} = [tuning_per_unit_table{unit+1, idx.(['in_CH_IS_' epoch '_PT_' tasktype])}];
                        eff_CH_CS{unit} = [tuning_per_unit_table{unit+1, idx.(['in_CH_CS_' epoch '_PT_' tasktype])}];
                        
                    end
                    
                end
                unit_IDs=tuning_per_unit_table(2:end, idx.unit_ID);
            end
            
        end
    end
end


disp('done');


end

