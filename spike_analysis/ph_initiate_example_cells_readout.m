function example_list=ph_initiate_example_cells_readout(project,versions)
%projects={'PPC_pulv_eye_hand','Pulv_eye_gaze_position','Pulv_eye_hand','PPC_pulv_eye_hand','STS_memory_saccades','Pulv_microstim_behavior'};
% if nargin==0
%     projects={'PPC_pulv_eye_hand','Pulv_eye_gaze_position','Pulv_eye_hand','PPC_pulv_eye_hand','STS_memory_saccades','Pulv_microstim_behavior'};
% else
%     projects=varargin{1};
% end
Selection_a={};%Selection={'stability_rating',1};
keys=struct();

%
% example_cell_definitions={'visual_enhanced',    'visual',1,     'in_NH_Cue_epoch',      'en', 'in_Cue_spaceLR',     'CS';...
%                           'visual_suppressed',  'visual',1,     'in_NH_Cue_epoch',      'su', 'in_Cue_spaceLR',     'IS';...
%                           'motor_enhanced',     'motor',1,      'in_NH_TIhol_epoch',    'en', 'in_TIhol_spaceLR',   '-';...
%                           'motor_suppressed',   'motor',1,      'in_NH_TIhol_epoch',    'su', 'in_TIhol_spaceLR',   '-';...
%                           'vismot_enhanced',    'visuomotor',1, 'in_NH_TIhol_epoch',    'en', 'in_TIhol_spaceLR',   'CS';...
%                           'vismot_suppressed',  'visuomotor',1, 'in_NH_TIhol_epoch',    'su', 'in_TIhol_spaceLR',   '-';...
%                           'mem_enhanced',       'in_NH_PeriS_epoch','en',          'in_NH_MemL_epoch',     'en', 'in_MemL_spaceLR',    'CS';...
%                           'mem_suppressed',     '','',          'in_NH_MemL_epoch',     'su', 'in_MemL_spaceLR',    '-';...
%                           'gaze_contra',        '','',          'in_NH_Thol_epoch',     'en', 'in_Thol_spaceLR',    'CS';...
%                           'gaze_ipsi',          '','',          'in_NH_Thol_epoch',     'su', 'in_Thol_spaceLR',    'CS'};
%
%
% 
example_cell_definitions={'retinotopic',            'in_AH_Cue_fixation_Msac_mov','false','in_AH_Cue_position_Msac_mov','true', 'in_AH_Cue_PxF_Msac_mov','false';...
    'gaze',                   'in_AH_Cue_fixation_Msac_mov','true', 'in_AH_Cue_position_Msac_mov','false','in_AH_Cue_PxF_Msac_mov','false';...
    'gaze_retinotopic_main',  'in_AH_Cue_fixation_Msac_mov','true', 'in_AH_Cue_position_Msac_mov','true', 'in_AH_Cue_PxF_Msac_mov','false';...
    'interaction',            'in_AH_Cue_PxF_Msac_mov','true',     'in_AH_Cue_PxF_Msac_mov','true',     'in_AH_Cue_PxF_Msac_mov','true';...
    };

keys.tt.choices=0;
keys.tt.tasktypes={'Msac_mov'};

% 
% example_cell_definitions={'fixation_gaze_effect',            'in_AH_Fhol_position_Msac_fix','true'   };
% keys.tt.choices=0;
% keys.tt.tasktypes={'Msac_fix'};

% 
% example_cell_definitions={'contra_but_no_position_effect','in_AH_Cue_position_Msac_mov','false','in_Cue_spaceLR_Msac_mov',     'CS';...
%                           'ipsi_but_no_position_effect','in_AH_Cue_position_Msac_mov','false','in_Cue_spaceLR_Msac_mov',     'IS'  }; 
                      
%                       
% example_cell_definitions={'fixThol_contracontra','in_Fhol_spaceLR_Msac_fix','CS','in_Thol_spaceLR_Msac_tar',     'CS';...
%                           'fixThol_ipsiipsi','in_Fhol_spaceLR_Msac_fix','IS','in_Thol_spaceLR_Msac_tar',     'IS';...
%                           'fixThol_contraipsi','in_Fhol_spaceLR_Msac_fix','CS','in_Thol_spaceLR_Msac_tar',     'IS';...
%                           'fixThol_ipsicontra','in_Fhol_spaceLR_Msac_fix','IS','in_Thol_spaceLR_Msac_tar',     'CS'}; 


%example_cell_definitions={'good_examples','in_AH_Cue_epoch_Msac_opt','en', 'in_Cue_spaceLR_Msac_opt','CS','in_AH_TIhol_epoch_Msac_opt','en', 'in_TIhol_spaceLR_Msac_opt','CS','in_AH_Thol_epoch_Msac_opt','en', 'in_Thol_spaceLR_Msac_opt',  'CS' };
% 
% example_cell_definitions={'hemi_wo_ANOVA_CS','in_Cue_spaceLR_Msac_opt','CS','in_AH_Cue_position_Msac_opt','false';...
%                           'hemi_wo_ANOVA_IS','in_Cue_spaceLR_Msac_opt','IS','in_AH_Cue_position_Msac_opt','false'};
% 
% keys.tt.choices=0;
% keys.tt.tasktypes={'Msac_opt'};
% 
% example_cell_definitions={'Initial_final_gaze_interaction','in_AH_Thol_PxF_Msac_tar','true','in_AH_Thol_position_Msac_tar','true';...
%                           'PreSGazetuning','in_AH_PreS_position_Msac_tar','true','in_AH_PreS_position_Msac_tar','true';...
%                           }; 
% keys.tt.choices=0;
% keys.tt.tasktypes={'Msac_mov'};

% example_cell_definitions={'Post_reti_only','in_AH_TIhol_position_Msac_mov','true','in_AH_TIhol_fixation_Msac_mov','false','in_AH_TIhol_PxF_Msac_mov','false';...
%                           'Post_gaze_only','in_AH_TIhol_position_Msac_mov','false','in_AH_TIhol_fixation_Msac_mov','true','in_AH_TIhol_PxF_Msac_mov','false';...
%                           'Post_both_main','in_AH_TIhol_position_Msac_mov','true','in_AH_TIhol_fixation_Msac_mov','true','in_AH_TIhol_PxF_Msac_mov','false';...
%                           'Post_reti_interaction','in_AH_TIhol_position_Msac_mov','true','in_AH_TIhol_fixation_Msac_mov','false','in_AH_TIhol_PxF_Msac_mov','true';...
%                           'Post_gaze_interaction','in_AH_TIhol_position_Msac_mov','false','in_AH_TIhol_fixation_Msac_mov','true','in_AH_TIhol_PxF_Msac_mov','true';...
%                           'Post_interaction_only','in_AH_TIhol_position_Msac_mov','true','in_AH_TIhol_fixation_Msac_mov','true','in_AH_TIhol_PxF_Msac_mov','true';...
%                           'Thol_reti_only','in_AH_Thol_position_Msac_mov','true','in_AH_Thol_fixation_Msac_mov','false','in_AH_Thol_PxF_Msac_mov','false';...
%                           'Thol_gaze_only','in_AH_Thol_position_Msac_mov','false','in_AH_Thol_fixation_Msac_mov','true','in_AH_Thol_PxF_Msac_mov','false';...
%                           'Thol_both_main','in_AH_Thol_position_Msac_mov','true','in_AH_Thol_fixation_Msac_mov','true','in_AH_Thol_PxF_Msac_mov','false';...
%                           'Thol_reti_interaction','in_AH_Thol_position_Msac_mov','true','in_AH_Thol_fixation_Msac_mov','false','in_AH_Thol_PxF_Msac_mov','true';...
%                           'Thol_gaze_interaction','in_AH_Thol_position_Msac_mov','false','in_AH_Thol_fixation_Msac_mov','true','in_AH_Thol_PxF_Msac_mov','true';...
%                           'Thol_interaction_only','in_AH_Thol_position_Msac_mov','true','in_AH_Thol_fixation_Msac_mov','true','in_AH_Thol_PxF_Msac_mov','true';...
%                           }; 
% keys.tt.choices=0;
% keys.tt.tasktypes={'Msac_mov'};
% 
% example_cell_definitions={'instructed contra',  'in_Cue_spaceLR_Msac_opt','CS','ch_Cue_spaceLR_Msac_opt','-';...
%     'instructed ipsi',  'in_Cue_spaceLR_Msac_opt','IS','ch_Cue_spaceLR_Msac_opt','-';...
%     'choice contra',  'in_Cue_spaceLR_Msac_opt','-','ch_Cue_spaceLR_Msac_opt','CS';...
%     'choice ipsi',  'in_Cue_spaceLR_Msac_opt','-','ch_Cue_spaceLR_Msac_opt','IS';...
%     'both contra',  'in_Cue_spaceLR_Msac_opt','CS','ch_Cue_spaceLR_Msac_opt','CS';...
%     'both ipsi',  'in_Cue_spaceLR_Msac_opt','IS','ch_Cue_spaceLR_Msac_opt','IS';...
%     'no cue tuning',  'in_Cue_spaceLR_Msac_opt','-','ch_Cue_spaceLR_Msac_opt','-';...
%     };
% keys.tt.choices=[0,1];
% keys.tt.tasktypes={'Msac_opt'};

keys.tt.IC_for_criterion='in';
keys.tt.hands=0;
mutually_exclusive=0;



% for p=1:numel(projects)
keys=ph_general_settings(project,keys);
project_specific_settings=[keys.db_folder filesep project filesep 'ph_project_settings.m'];
run(project_specific_settings);
for f=1:numel(versions) % running multiple versions of the same project at once !
    %         if ~isempty(keys.project_versions{f})
    %             keys.project_version=keys.project_versions{f};
    %         end
    %         if nargin>1
    keys.project_version=versions{f};
    % end
    version_folder=keys.project_version;
    keys.version_specific_settings=[keys.db_folder project filesep keys.project_version filesep 'ph_project_version_settings.m'];
    run(keys.version_specific_settings);
    keys.project_version=version_folder;
    keys.anova_table_file=[keys.drive '\' keys.basepath_to_save '\' keys.project_version '\tuning_table_combined_CI.mat'];
    
    if keys.batching.combine_monkeys
        keys.batching.monkeys={''};
    end
    
    for m=1:numel(keys.batching.monkeys)
        keys.monkey=keys.batching.monkeys{m};
        for c=1:numel(keys.position_and_plotting_arrangements)
            keys.case={keys.position_and_plotting_arrangements{c}(1:3)};
            for t=1:numel(keys.batching.targets) %% target= very critical !!!
                Selection_t=[Selection_a; {'target',keys.batching.targets{t}}];
                %  n_tasks=numel(keys.cc.conditions_to_plot); %% cc??
                %
                %                     for k=1:n_tasks
                %                         tasktypes=keys.cc.conditions_to_plot(k); %% cc??
                for subregion=1:keys.batching.n_Subregions
                    if keys.batching.Subregions_separately
                        keys.tt.selection=[Selection_t; {'Subregion', subregion}];
                    else
                        keys.tt.selection=Selection_t;
                    end
                    
                    %% keys.cc.instructed_choice for criterion???
                    % for choices=keys.cc.choices
                    %for hands=keys.cc.hands
                    
                    
                    % CASES ???
                    %=strcat(tasktypes, ['_' keys.case{1}]);
                    
                    
                    [tuning_per_unit_table]=ph_load_extended_tuning_table(keys);
                    [tuning_per_unit_table, keys.selection_title]=ph_reduce_tuning_table(tuning_per_unit_table,keys);
                    tuning_per_unit_table(cellfun(@(x) isempty(x) & ~islogical(x),tuning_per_unit_table))={''};
                    idx_ID=find_column_index(tuning_per_unit_table,'unit_ID');
                    for ex=1:size(example_cell_definitions,1)
                        row_idx=true(size(tuning_per_unit_table,1),1);
                        
                        for par=2:2:size(example_cell_definitions,2)
                            entry=example_cell_definitions{ex,par+1};
                            column_index=find_column_index(tuning_per_unit_table,[example_cell_definitions{ex,par}]);
                            
                            
                            if ischar(entry) && any(column_index)
                                tuning_per_unit_table(2:end,column_index)=cellfun(@num2str,tuning_per_unit_table(2:end,column_index),'UniformOutput',0); %% crazy cellstring bug
                                row_idx=row_idx & [false; ~cellfun(@isempty,strfind(tuning_per_unit_table(2:end,column_index),entry))];
                                
                            elseif any(column_index) %% is number!?
                                row_idx=row_idx & [false; cell2mat(tuning_per_unit_table(2:end,column_index))==entry];
                            end
                            
                        end
                        example_list(ex).unit_IDs=tuning_per_unit_table(row_idx,idx_ID);
                        
                        example_list(ex).definition=example_cell_definitions{ex,1};
                    end
                    if mutually_exclusive
                        example_indexes=1:numel(example_list);
                        for ex=example_indexes
                            example_list_exclusive(ex).unit_IDs=example_list(ex).unit_IDs(~ismember(example_list(ex).unit_IDs,vertcat(example_list(example_indexes~=ex).unit_IDs)));
                            example_list_exclusive(ex).definition=example_list(ex).definition;
                        end
                        example_list=example_list_exclusive;
                    end
                end
            end
        end
    end
end
%end
end
