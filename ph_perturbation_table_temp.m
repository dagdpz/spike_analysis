function ph_perturbation_table_temp(project,versions)
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
        %avg_unit_epoch_L=[];
        %avg_unit_epoch_R=[];
        for t=1:numel(keys.batching.targets)
            target=keys.batching.targets{t};
            
            keys.tt.selection                   ={'target',target};
            tasktype=keys.tt.tasktypes{1};
            keys.anova_table_file= ['Y:\Projects\' project '\ephys\' version_folder '\tuning_table_combined_CI.mat'];
            [tuning_per_unit_table]                 = ph_load_extended_tuning_table(keys);
            [tuning_per_unit_table, Sel_for_title]  = ph_reduce_tuning_table(tuning_per_unit_table,keys);
            
            if ~exist([keys.basepath_to_save filesep keys.project_version filesep 'perturbation_table'],'dir')
                mkdir([keys.basepath_to_save filesep keys.project_version],'perturbation_table');
            end
            
            
            for kk=1:size(tuning_per_unit_table,2)
                idx.(tuning_per_unit_table{1,kk})=kk;
            end
            
            
            if any(strfind(keys.project_version,'MIP'))
                epochs={'INI', 'Facq','Fhol','Cue','EDel','Del','PreR','PeriR','PostR','Thol'};
                
            elseif  any(strfind(keys.project_version,'LIP'))
                epochs={'INI', 'Facq','Fhol','Cue','EDel','Del','PreR','PeriR','PostR','Thol'};
            end
            
            figure
            for e=1:numel(epochs)
                
                epoch=epochs{e};
                
                subplot(subplot_rows,numel(epochs)/2,e);
                
                HSC = [1 2 3 4 5 6 7]; %Four Hand-Space conditons +1
                
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
%                         Df_L=vertcat(Df_L, zeros(1,size(eff_L,2)+1));
                        units_L=[1:size(vertcat(tuning_per_unit_table(2:size(tuning_per_unit_table,1)- 1,1) ,'-'))];
                        
                        %create average response per unit per epoch L side(BG)
                        for counter_average_L = 1 : size( Df_L,1)
                            if sum((Df_L(counter_average_L, :))==1)>=2 && sum((Df_L(counter_average_L, :))==2)<2 %at least 2 red but no more than 1 black
                                avg_unit_response_L (counter_average_L,:)= [1 0];
                            elseif sum((Df_L(counter_average_L, :))==2)>=2 && sum((Df_L(counter_average_L, :))==1)<2 % at least 2 black but no more than 1 red
                                avg_unit_response_L (counter_average_L,:)= [2 0];
                            elseif sum((Df_L(counter_average_L, :))==0)>=2 && ((sum((Df_L(counter_average_L, :))==1))+(sum((Df_L(counter_average_L, :))==2)))<2 %at least 2 white but not 1 red and 1 black
                                avg_unit_response_L (counter_average_L, :)= [0 0];
                            else 
                                avg_unit_response_L (counter_average_L,:)= [3 0]; %inconsistent cases
                            end
                        end
                        
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
%                         Df_R=vertcat(Df_R, zeros(1,size(eff_R,2)+1));
                        units_R=[1:size(vertcat(tuning_per_unit_table(2:size(tuning_per_unit_table,1) - 1,1) ,'-'))];
                        
                        %create average response per unit per epoch R side(BG)
                        avg_unit_response_R = [];
                        for counter_average_R = 1 : size( Df_R, 1)
                            if sum((Df_R(counter_average_R, :))==1)>=2 && sum((Df_R(counter_average_R, :))==2)<2 %at least 2 red but no more than 1 black
                                avg_unit_response_R (counter_average_R,:)= [1 0];
                            elseif sum((Df_R(counter_average_R, :))==2)>=2 && sum((Df_R(counter_average_R, :))==1)<2 % at least 2 black but no more than 1 red
                                avg_unit_response_R (counter_average_R,:)= [2 0];
                            elseif sum((Df_R(counter_average_R, :))==0)>=2 && ((sum((Df_R(counter_average_R, :))==1))+(sum((Df_R(counter_average_R, :))==2)))<2 %at least 2 white but not 1 red and 1 black
                                avg_unit_response_R (counter_average_R,:)= [0 0];
                            else 
                                avg_unit_response_R (counter_average_R,:)= [3 0]; %inconsistent cases
                            end
                        end
                    end

                end
                
                unit_IDs=tuning_per_unit_table(2:end, idx.unit_ID);

                if any(strfind(target,'_L'))
                    to_plot_L = horzcat(Df_L,avg_unit_response_L);
                    pcolor(HSC,units_L,to_plot_L);
                    
                    %create table for left hemisphere with unit_ID and
                    %average response per epoch (BG)
                    %unit_IDs_L = [unit_IDs];
                    %avg_unit_epoch_L = horzcat(avg_unit_epoch_L, avg_unit_response_L(:,1));
                    avg_unit_epoch_MIP_L.unit_ID = [unit_IDs];
                    avg_unit_epoch_MIP_L.(epochs{e}) = avg_unit_response_L(:,1); 
                    %avg_unit_epoch_L_table = table(unit_IDs_L,  avg_unit_epoch_L); 
                else
                    to_plot_R = horzcat(Df_R,avg_unit_response_R);
                    pcolor(HSC,units_R,to_plot_R);

                    %create table for right hemisphere with unit_ID and
                    %average response per epoch (BG)
                    %avg_unit_epoch_R = horzcat(avg_unit_epoch_R, avg_unit_response_R(:,1));
                    %unit_IDs_R = [unit_IDs];
                    %avg_unit_epoch_R_table = table(unit_IDs_R,  avg_unit_epoch_R); 
                    avg_unit_epoch_MIP_R.unit_ID = [unit_IDs];
                    avg_unit_epoch_MIP_R.(epochs{e}) = avg_unit_response_R(:,1);

                end
                x=HSC;
                set(gca,'XTick',[1.5,2.5,3.5,4.5, 5.5, 6.5])
                set(gca,'XTickLabel',{'IH-IS','IH-CS','CH-IS','CH-CS', '', 'Avg Response'})
                map = [1 1 1;  1 0 0; 0 0 0; 0 1 0];%colors
                cb = colorbar;
                cb.Label.String = 'Effect of Inactivation';
                caxis([0 3]);
                %colorbar('LimitsMode', 'manual', 'Limits', [0 3], 'Ticks',[0,1,2,3],'TickLabels',{'No Effect','Enhance','Supression','Inconsistent'})
                %set (hh, 'ylim', [0 3]);
                colormap(map)
                title(epoch);
               
            end
            figure_title=[keys.monkey '_' target '_' tasktype];
            filename=[keys.basepath_to_save filesep keys.project_version filesep 'perturbation_table' filesep figure_title];
            mtit(figure_title, 'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none');
            export_fig(filename, '-pdf','-transparent') % pdf by run
            
        end
        %avg_unit_epoch_L_table = splitvars(avg_unit_epoch_L_table, 'avg_unit_epoch_L', 'NewVariableNames',{'INI', 'Facq','Fhol','Cue','EDel','Del','PreR','PeriR','PostR','Thol'});
        %avg_unit_epoch_R_table = splitvars(avg_unit_epoch_R_table, 'avg_unit_epoch_R', 'NewVariableNames',{'INI', 'Facq','Fhol','Cue','EDel','Del','PreR','PeriR','PostR','Thol'});
        %avg_unit_epoch_MIP_L = table2struct (avg_unit_epoch_L_table);
        %avg_unit_epoch_MIP_R = table2struct (avg_unit_epoch_R_table);     
        struct_path =  [keys.basepath_to_save filesep keys.project_version filesep 'perturbation_table'];
        save([struct_path filesep keys.monkey '_MIP_L_' tasktype '_' 'avg_units_struct.mat'], 'avg_unit_epoch_MIP_L');
        save([struct_path filesep keys.monkey '_MIP_R_' tasktype '_' 'avg_units_struct.mat'], 'avg_unit_epoch_MIP_R');

    end

end

% %                    xlabel('Hand Space Conditions')
% %
% %                    ylabel('Units')
close all

disp('done');


end

