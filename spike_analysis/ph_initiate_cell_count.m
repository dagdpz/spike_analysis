function ph_initiate_cell_count(varargin)
%projects={'PPC_pulv_eye_hand','Pulv_eye_gaze_position','Pulv_eye_hand','PPC_pulv_eye_hand','STS_memory_saccades','Pulv_microstim_behavior'};
if nargin==0
    projects={'PPC_pulv_eye_hand','Pulv_eye_gaze_position','Pulv_eye_hand','PPC_pulv_eye_hand','STS_memory_saccades','Pulv_microstim_behavior'};
else
    projects=varargin{1};
end

Selection_a={};%Selection={'stability_rating',1};

in_or_ch={'in'};
%in_or_ch={'ch'};

keys=struct();

for p=1:numel(projects)
    project=projects{p};
    keys=ph_general_settings(project,keys);
    project_specific_settings=[keys.db_folder filesep project filesep 'phys_settings.m'];
    run(project_specific_settings);
    for f=1:numel(keys.pdf_folders) % running multiple versions of the same project at once !
        if ~isempty(keys.pdf_folders{f})
            keys.pdf_folder=keys.pdf_folders{f};
        end
        if nargin>1
            keys.pdf_folder=varargin{2};
        end
        pdf_folder=keys.pdf_folder;
        keys.version_specific_settings=[keys.db_folder project filesep keys.pdf_folder filesep 'phys_settings.m'];
        run(keys.version_specific_settings);
        keys.pdf_folder=pdf_folder;
        %% move from here
        
        keys.tt.only_overlapping_tasktypes=1; %%
        % keys.tt.only_both_hands        =1; %% ?
        
        keys.anova_table_file=[keys.drive '\' keys.basepath_to_save '\' keys.pdf_folder '\tuning_table_combined_CI.mat'];
        
        if keys.combine_monkeys
            keys.monkeys={''};
        end
        for summary={'space_and_epoch'} %%{'per_epoch','per_task','space_x_hand','space_and_epoch'} %
            keys.cc.plot_type=summary;
            for ic=1:numel(in_or_ch)
                keys.instructed_choice=in_or_ch{ic};
                for m=1:numel(keys.monkeys)
                    keys.monkey=keys.monkeys{m};
                    %                     if ~isempty(keys.monkey)
                    %                         Selection_m=[Selection_a; {'unit_ID', keys.monkey(1:3)}];
                    %                     else
                    %                         Selection_m=Selection_a;
                    %                     end
                    for c=1:numel(keys.case_summaries) %%% need to rename this one....
                        keys.case={keys.case_summaries{c}(1:3)};
                        for t=1:numel(keys.targets) %% target= very critical !!!
                            Selection_t=[Selection_a; {'target',keys.targets{t}}];
                            %keys.target=keys.targets{t};
                            
                            % this stuff should probably be inside the cell
                            % count function itself (tasktype selection)
                            
                            if ismember(keys.cc.plot_type,{'per_epoch'})
                                n_tasks=1;
                            else
                                n_tasks=numel(keys.cc.conditions_to_plot);
                            end
                            if ismember(keys.cc.plot_type,{'per_epoch','per_task'})
                                factor_loop=keys.cell_count_factors;
                            else
                                factor_loop={''};
                                n_tasks=numel(keys.cc.conditions_to_plot);
                            end
                            for factor_cell=factor_loop
                                keys.cc.factors=factor_cell{:};
                                if ismember(keys.cc.factors,{'choice1','choice2'})
                                    n_tasks=numel(keys.cc.conditions_to_plot);
                                        keys.tt.only_overlapping_tasktypes=1;
                                end
                                for k=1:n_tasks
                                    if ismember(keys.cc.plot_type,{'per_epoch'}) && ~ismember(keys.cc.factors,{'choice1','choice2'})% %ismember(keys.cc.plot_type,{'per_epoch','per_task'}) || ismember(keys.cc.factors,{'choice1','choice2'})%% to simplify
                                        tasktypes=keys.cc.conditions_to_plot;
                                        keys.tt.only_overlapping_tasktypes=1;
                                    else
                                        tasktypes=keys.cc.conditions_to_plot(k);
                                    end
                                    keys.cc.epochs={};
                                    for d=1:numel(tasktypes)
                                        if strcmp(keys.cc.plot_type,'per_task')
                                            keys.cc.epochs=[keys.cc.epochs; keys.cc.epochsS.(tasktypes{d})(~ismember(keys.cc.epochsS.(tasktypes{d}),keys.cc.epochs))];
                                        else
                                            keys.cc.epochs=[keys.cc.epochs; keys.cc.epochsE.(tasktypes{d})(~ismember(keys.cc.epochsE.(tasktypes{d}),keys.cc.epochs))];
                                        end
                                    end
                                    for subregion=1:keys.n_Subregions
                                        if keys.Subregions_separately
                                            keys.tt.selection=[Selection_t; {'Subregion', subregion}];
                                        else
                                            keys.tt.selection=Selection_t;
                                        end
                                        for plot_as_pie=1%:1
                                            for percent=0:1
                                                keys.cc.percent=percent;
                                                keys.cc.plot_as_pie=plot_as_pie;
                                                keys.cc.tasktypes=tasktypes;
                                                keys.tt.tasktypes=strcat(tasktypes, ['_' keys.case{1}]);
                                                keys.tt.only_both_hands             =keys.cc.only_both_hands;
                                                
%                                                 %  for fixation x target anovas
%                                                 keys.case{1}='mov';
%                                                 keys.case{2}='tar';
                                                ph_anova_cell_count(keys); %% main function
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

