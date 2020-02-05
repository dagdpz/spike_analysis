function filelist=ph_get_filelist(varargin)
project=varargin{1};

keys=struct;
keys.combine_monkeys=0;
keys=ph_general_settings(project,keys);
project_specific_settings=[keys.db_folder project filesep 'ph_project_settings.m'];
run(project_specific_settings)

keys.project_versions=varargin{2};
for f=1:numel(keys.project_versions) % running multiple versions of the same project at once !
    keys.project_version=keys.project_versions{f};
    version_folder=keys.project_version;
    keys.version_specific_settings=[keys.db_folder project filesep keys.project_version filesep 'ph_project_version_settings.m'];
    run(keys.version_specific_settings);
    keys.project_version=version_folder;
    
    keys.folder_to_save=[keys.drive filesep keys.basepath_to_save filesep keys.project_version];
    keys.subfolder_prefix='scatter';
    
    filelist=struct;
    filelist_formatted=struct;
    selection=keys.tt.selection;
    for a=1:numel(keys.position_and_plotting_arrangements)
        current_case=keys.position_and_plotting_arrangements{a};
        keys.tt.tasktypes=strcat(keys.tt.type_effectors,'_', current_case(1:3)); %% conditions_to_plot is not ideal here
        for m=1:numel(keys.batching.monkeys)
            keys.monkey=keys.batching.monkeys{m};
            keys.anova_table_file=[keys.drive '\Projects\' project '\ephys\' keys.project_version '\tuning_table_combined_CI.mat'];
            if nargin<3
                population=ph_load_population([keys.drive filesep keys.basepath_to_save filesep keys.project_version],['population_' keys.monkey]);
            else
                population=varargin{3};
            end
            
            for t=1:numel(keys.batching.targets)
                target=keys.batching.targets{t};
                keys.tt.selection                    = [selection; {'target',target}];
                [tuning_table]=ph_load_extended_tuning_table(keys);
                [tuning_table, keys.selection_title]=ph_reduce_tuning_table(tuning_table,keys);
                [~,units_valid,~]=intersect({population.unit_ID},tuning_table(:,1));
                pop=population(units_valid);
                pop_trials=[pop.trial];
                for pt=1:numel(keys.cal.perturbation_groups)
                    pop_trials_pt=pop_trials(ismember([pop_trials.perturbation],keys.cal.perturbation_groups{pt}));
                    if isempty(pop_trials_pt)
                       continue; 
                    end
                    blocks=unique([[pop_trials_pt.date]',[pop_trials_pt.block]'],'rows');
                    sessions=unique(blocks(:,1));
                    for s=1:numel(sessions)
                        session=sessions(s);
                        idx=blocks(:,1)==session;
                        base_path=['Y:\Data\' keys.monkey '_phys_combined_monkeypsych_TDT\' num2str(session)];
                        runs=run2block(base_path,blocks(idx,2),1);
                        blocks(idx,3)=runs;
                        filelist_formatted.([keys.monkey(1:3) '_' target '_PT' num2str(pt-1) '_' keys.tt.tasktypes{:}]){s,1}=session;
                        filelist_formatted.([keys.monkey(1:3) '_' target '_PT' num2str(pt-1) '_' keys.tt.tasktypes{:}]){s,2}=runs;
                    end
                    runs0=arrayfun(@(x) num2str(x), blocks(:,1),'UniformOutput',false);
                    runs1=arrayfun(@(x) num2str(floor(x/10000)), blocks(:,1),'UniformOutput',false);
                    runs2=arrayfun(@(x) sprintf('%02d',mod(x,100)), blocks(:,1),'UniformOutput',false);
                    runs3=arrayfun(@(x) sprintf('%02d',mod(floor(x/100),100)), blocks(:,1),'UniformOutput',false);
                    runs4=arrayfun(@(x) sprintf('%02d',x), blocks(:,3),'UniformOutput',false);
                    
                    
                    current_filelist=strcat(['Y:\Data\' keys.monkey '_phys' filesep], runs0, [keys.monkey(1:3)], runs1, '-', runs2, '-', runs3, '_', runs4,'.mat');
                    
                    filelist.([keys.monkey(1:3) '_' target '_PT' num2str(pt-1) '_' keys.tt.tasktypes{:}])=current_filelist;
                end
            end
        end
        
    end
    save([keys.folder_to_save filesep 'behaviour_filelist.mat'],'filelist','filelist_formatted');
end

end


function output_blocks_or_runs=run2block(base_path,given_blocks_or_runs,reverse)
folder_content=dir([base_path filesep '*.mat']);
complete_filelist=vertcat(folder_content.name);
filelist_blocks=str2num(complete_filelist(:,end-5:end-4));
filelist_runs=str2num(complete_filelist(:,end-14:end-13));

if strcmp(reverse,'reverse') || reverse==1
    [~,idx] =ismember(given_blocks_or_runs,filelist_blocks);
    output_blocks_or_runs(idx~=0)=filelist_runs(idx(idx~=0));
else
    [~,idx] =ismember(given_blocks_or_runs,filelist_runs);
    output_blocks_or_runs(idx~=0)=filelist_blocks(idx(idx~=0));
end
output_blocks_or_runs(idx==0)=0;
end
