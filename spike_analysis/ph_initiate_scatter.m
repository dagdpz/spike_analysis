
keys.plot                               =0;
keys.plot_waveforms                     =0;
keys.delete_previous_pdfs               =0;
keys.create_pdfs                        =1;
keys.append_pdfs                        =0;
keys.list_successful_only               =1;
keys.minimum_trials_per_condition       =60;
keys.min_n_spikes_per_unit              =50;
keys.process_only_units_in_the_table    =1;

keys.monkey_colors          =[0 1 0; 1 0 0];

%project='Pulv_microstim_behavior';%'PPC_pulv_eye_hand'; %'Pulv_microstim_behavior'; %'Pulv_eye_hand'; % 'STS_memory_saccades'
%projects={'Pulv_eye_gaze_position'};
projects={'Pulv_microstim_behavior'};
%keys=ph_settings(project,keys);
for p=1:numel(projects)
    project=projects{p};
    keys=ph_general_settings(project,keys);
    project_specific_settings=[keys.db_folder project filesep 'phys_settings.m'];
    run(project_specific_settings)
%     if nargin>=2
%         keys.pdf_folders=versions;
%     end
    
    for f=1:numel(keys.pdf_folders) % running multiple versions of the same project at once !
        if ~isempty(keys.pdf_folders{f})
            keys.pdf_folder=keys.pdf_folders{f};
        end
        pdf_folder=keys.pdf_folder;
        keys.version_specific_settings=[keys.db_folder project filesep keys.pdf_folder filesep 'phys_settings.m'];
        run(keys.version_specific_settings)
        keys.pdf_folder=pdf_folder;
        
        keys.SX.instructed_choice   ='in';
        keys.SX.epoch               ='Cue';
        keys.SX.parameter           ='spaceLR';
        keys.SX.dataset             ='Msac';
        keys.SX.case                ='opt';
        
        keys.SY.instructed_choice   ='in';
        keys.SY.epoch               ='Cue';
        keys.SY.parameter           ='spaceLR';
        keys.SY.dataset             ='Vsac';
        keys.SY.case                ='opt';

             keys.tt.tasktypes={'Msac_opt','Vsac_opt'};         
%         
%         keys.SX.instructed_choice   ='in_NH';
%         keys.SX.epoch               ='Cue';
%         keys.SX.parameter           ='fixation_PV';
%         keys.SX.dataset             ='Msac';
%         keys.SX.case                ='mov';
%         
%         keys.SY.instructed_choice   ='in_NH';
%         keys.SY.epoch               ='Cue';
%         keys.SY.parameter           ='position_PV';
%         keys.SY.dataset             ='Msac';
%         keys.SY.case                ='mov';
%         
%             %%?
%             keys.tt.tasktypes={'Msac_mov'};
            
            
        keys.table_full_path=[keys.drive '\' keys.basepath_to_save '\' keys.pdf_folder '\tuning_table_combined.mat'];
        keys.folder_to_save=[keys.drive filesep keys.basepath_to_save filesep keys.pdf_folder];
        keys.subfolder_prefix='scatter';
            keys.anova_table_file=[keys.drive '\Projects\' project '\ephys\' keys.pdf_folder '\tuning_table_combined_CI.mat'];
        
            
        for t=1:numel(keys.targets)
           % keys.target=[target{:}];
             keys.tt.selection={'target',keys.targets{t}};
            for comparison=1:6
                switch comparison
                    case 1
                        keys.SX.epoch               ='Cue';
                        keys.SY.epoch               ='Cue';
                    case 2
                        keys.SX.epoch               ='PreS';
                        keys.SY.epoch               ='PreS';
                    case 3
                        keys.SX.epoch               ='PeriS';
                        keys.SY.epoch               ='PeriS';
                    case 4
                        keys.SX.epoch               ='Thol';
                        keys.SY.epoch               ='Thol';
                    case 5
                        keys.SX.epoch               ='TIhol';
                        keys.SY.epoch               ='Tacq';
                end
%                 irrelevant
%              keys.tt.epoch_criterion='none';
%              keys.tt.space_criterion='none';
              keys.tt.only_overlapping_tasktypes=0;
                keys.monkey             = '';
                ph_scatter(keys);
            end
        end
    end
end