function population=ph_torn_MUs(population,modified_keys)

for fn=fieldnames(modified_keys)'
    keys.(fn{:})=modified_keys.(fn{:});
end

%% tuning table preparation and grouping
[tuning_per_unit_table]                 = ph_load_extended_tuning_table(keys);
[tuning_per_unit_table, Sel_for_title]  = ph_reduce_tuning_table(tuning_per_unit_table,keys);
if keys.TO.FR_subtract_baseline
    Sel_for_title =[Sel_for_title,{'base';'=';keys.TO.epoch_BL;', '}];
end
idx_group_parameter=DAG_find_column_index(tuning_per_unit_table,keys.TO.group_parameter);
group_values=tuning_per_unit_table(:,idx_group_parameter);
group_values=cellfun(@num2str, group_values, 'UniformOutput', false);
cell_in_any_group=[false; ~ismember(group_values(2:end),keys.TO.group_excluded)];
unique_group_values=unique(group_values(cell_in_any_group));
if isempty(unique_group_values)
    disp('no relevant groups found');
    return;
end
complete_unit_list={population.unit_ID}';
%
idx_unitID=DAG_find_column_index(tuning_per_unit_table,'unit_ID');


%% som eundefined stuff
dropboxpath=DAG_get_Dropbox_path;
keys.batching.monkeys={'Linus','Curius'};

keys.SNR_single=[1,2];
keys.SiR_single=1;
keys.StR_single=1;

sorting_table={};
for m=1:numel(keys.batching.monkeys)
    monkey=keys.batching.monkeys{m};
    keys.sorted_neurons_filename =keys.(monkey).sorted_neurons_filename;
    keys.sorted_neurons_sheetname=keys.(monkey).sorted_neurons_sheetname;
    keys.date=keys.(monkey).date;
    keys.monkey                     =[monkey '_phys'];
    keys.sorted_neurons_foldername  =[dropboxpath filesep 'DAG' filesep 'phys' filesep monkey '_phys_dpz'];
    keys.tuning_table_foldername    =[keys.drive filesep keys.basepath_to_save filesep keys.project_version];
    keys.tuning_table_filename      =['tuning_table_combined'];
    keys.population_foldername      =[keys.drive filesep keys.basepath_to_save filesep keys.project_version];
    keys.population_filename        =['population_' monkey];
    keys.sites_filename             =['sites_' monkey];
    keys.filelist_formatted         =keys.(monkey).filelist_formatted;
    
    
    temp_xlsx=dir(fullfile(keys.sorted_neurons_foldername,[keys.sorted_neurons_filename '*.xls*']));
    if ~isempty(temp_xlsx)
        excel_sheet= [keys.sorted_neurons_foldername filesep temp_xlsx(1).name];
        [~,~,ST] = xlsread(excel_sheet,keys.sorted_neurons_sheetname); 
        sorting_table(1,:)=ST(1,1:25); %%% 1:25 because size is not identical
        sorting_table=[sorting_table;ST(2:end,1:25)];
    end
    
end

ID_ix=DAG_find_column_index(sorting_table,'Neuron_ID');
block_ix=DAG_find_column_index(sorting_table,'Block');
SNR_ix=DAG_find_column_index(sorting_table,'SNR_rank');
Single_ix=DAG_find_column_index(sorting_table,'Single_rank');
Stability_ix=DAG_find_column_index(sorting_table,'Stability_rank');

for g=1:numel(unique_group_values)
    clear suindexes
    unitidx=ismember(complete_unit_list,tuning_per_unit_table(ismember(group_values,unique_group_values(g)),idx_unitID));
    units=find(all(unitidx,2))';
    whatsthis=[population(units).trial];
    alltypes=[whatsthis.type];
    u_types=unique(alltypes);
    
    for typ=u_types
        n=0;
        for u=units
            n=n+1;
            ID=complete_unit_list{u};
            trials=[population(u).trial];
            trials=trials([trials.type]==typ);
            if size(trials,2)>0
                blocks=unique([trials.block]);
                rows=[false; ismember(sorting_table(2:end,ID_ix),ID) & ismember(vertcat(sorting_table{2:end,block_ix}),blocks)];
                SNR=[sorting_table{rows,SNR_ix}];
                SiR=[sorting_table{rows,Single_ix}];
                StR=[sorting_table{rows,Stability_ix}];
            else
                SNR=NaN;
                SiR=NaN;
                StR=NaN;
            end
            
            unit_torn=  (any(ismember(SNR,keys.SNR_single)) & any(~ismember(SNR,keys.SNR_single))) |...
                (any(ismember(SiR,keys.SiR_single)) & any(~ismember(SiR,keys.SiR_single))) |...
                (any(ismember(StR,keys.StR_single)) & any(~ismember(StR,keys.StR_single)));
            
            
            per_type(typ).rankings(n).SNR=SNR;
            per_type(typ).rankings(n).SiR=SiR;
            per_type(typ).rankings(n).StR=StR;
            per_type(typ).rankings(n).unit_torn=unit_torn;
            
        end
    end
end






for g=1:numel(unique_group_values)
    clear suindexes
    unitidx=ismember(complete_unit_list,tuning_per_unit_table(ismember(group_values,unique_group_values(g)),idx_unitID));
    units=find(all(unitidx,2))';
    onerunonly=arrayfun(@(x) numel(x.block_unit)==3,population(units));
    IDS={population(units).unit_ID};
    LL=IDS(~onerunonly);
    LL_index=ismember(IDS,LL);
    singlerating=[population(units).Single_rating];
    stabilityrating=[population(units).stability_rating];
    
    suindexes(1,:)=singlerating==1 & stabilityrating==1;
    suindexes(2,:)=singlerating==2 | (singlerating==1 & stabilityrating<1) & ~ LL_index;
    suindexes(3,:)=singlerating>2;
    LL_index=LL_index & ~suindexes(3,:) & ~suindexes(1,:);
    
    notenoughtrials_forsingle_idx=ismember(IDS,{'Cur_20151204_22','Cur_20151204_24'}); % excluded from singleunitanalysis
    
    suindexes(1,:)= suindexes(1,:) | (LL_index &~ notenoughtrials_forsingle_idx);
    suindexes(2,:)=suindexes(2,:) | notenoughtrials_forsingle_idx;
    
    
    
    
    FR=[population(units).FR];
    bins=2.5:5:ceil(max(FR)*10)/10;
    
    %% figure 1 SU/SU+/MU + six units,
    
    plot_title=[keys.tt.tasktypes{:} ' FR SU_SU+_MU with torn units'];
    filename=plot_title;
    figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
    histoindex= suindexes ;
    histoindex(end+1,:)=true([1 size(histoindex,2)]);
    
    for h=1:size(histoindex,1)
        histo(:,h)=hist(FR(histoindex(h,:)),bins);
        means(h)=nanmean(FR(histoindex(h,:)));
        medians(h)=nanmedian(FR(histoindex(h,:)));
        stds(h)=nanstd(FR(histoindex(h,:)));
    end
    %bar(bins,histo,'stacked');
    cols=colormap;
    for h=1:size(histoindex,1)
        subplot(size(histoindex,1),1,h)
        bar(bins,histo(:,h));
        y_lim=get(gca,'ylim');
        x_lim=get(gca,'xlim');
        texttoplot=['M:' num2str(means(h)) ' +' num2str(stds(h)) ', med:'  num2str(medians(h)) ];
        text(x_lim(1)+diff(x_lim)/20,y_lim(2)-diff(y_lim)/40*h,texttoplot,'color',cols(h,:));
    end
    legend({'SU','SU+','MU'})
    xlabel('Firing rate')
    ylabel('N units');
    ph_title_and_save(figure_handle,filename,plot_title,keys)
    
    %% figure 2 SU + six units/MU,
    
    suindexes(3,:)=suindexes(2,:) | suindexes(3,:);
    suindexes(2,:)=LL_index &~ notenoughtrials_forsingle_idx;
    suindexes(1,:)=singlerating==1 & stabilityrating==1;
    
    plot_title=[keys.tt.tasktypes{:}  ' FR SU_MU with torn units'];
    filename=plot_title;
    figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
    histoindex= suindexes ;
    histoindex(end+1,:)=true([1 size(histoindex,2)]);
    
    for h=1:size(histoindex,1)
        histo(:,h)=hist(FR(histoindex(h,:)),bins);
        means(h)=nanmean(FR(histoindex(h,:)));
        medians(h)=nanmedian(FR(histoindex(h,:)));
        stds(h)=nanstd(FR(histoindex(h,:)));
    end
    bar(bins,histo(:,end-1),'stacked');
    cols=colormap;
    y_lim=get(gca,'ylim');
    x_lim=get(gca,'xlim');
    for h=1:size(histoindex,1)
        texttoplot=['M:' num2str(means(h)) ' +' num2str(stds(h)) ', med:'  num2str(medians(h)) ];
        text(x_lim(1)+diff(x_lim)/20,y_lim(2)-diff(y_lim)/40*h,texttoplot,'color',cols(h,:));
    end
    legend({'SU','torn','MU'})
    xlabel('Firing rate')
    ylabel('N units');
    ph_title_and_save(figure_handle,filename,plot_title,keys)
    
    
    %% figure 3 SU/MU with gain field units,
    
    e=3;
    EPOCHS=keys.EPOCHS_PER_TYPE{3};
    clear sugainindex
    
    popuation_total=population;
    load('Y:\Projects\Pulv_eye_gaze_position\ephys\20190309\tuning_curves\population_FR_vector.mat')
    RF_per_epoch=vertcat(population.RF_per_epoch); % cue
    RF_per_epoch_cue=RF_per_epoch(units,e);
    
    positionanovaindex=DAG_find_column_index(tuning_per_unit_table,['in_AH_' EPOCHS{e,1} '_position_' keys.tt.tasktypes{:}]);
    gazeanovaindex=DAG_find_column_index(tuning_per_unit_table,['in_AH_' EPOCHS{e,1} '_fixation_' keys.tt.tasktypes{:}]);
    interactionanovaindex=DAG_find_column_index(tuning_per_unit_table,['in_AH_' EPOCHS{e,1} '_PxF_' keys.tt.tasktypes{:}]);
    unit_IDs_over_all_position_effect=tuning_per_unit_table(strcmp(['false'; tuning_per_unit_table(2:end,positionanovaindex)],'true'),1);
    unit_IDs_gazedependence=tuning_per_unit_table(strcmp(['false'; tuning_per_unit_table(2:end,gazeanovaindex)],'true') | strcmp(['false'; tuning_per_unit_table(2:end,interactionanovaindex)],'true'),1);
    
    ANOVA_position_effect=ismember(IDS,unit_IDs_over_all_position_effect);
    gaze_dependence=ismember(IDS,unit_IDs_gazedependence)';
    %popuation_total=population;
    sugainindex(1,:)=(suindexes(1,:)|suindexes(2,:))&~(ANOVA_position_effect&[RF_per_epoch_cue.significant_gain]);
    sugainindex(2,:)=(suindexes(1,:)|suindexes(2,:))&ANOVA_position_effect&[RF_per_epoch_cue.significant_gain];
    %     sugainindex(3,:)=suindexes(2,:)&~(ANOVA_position_effect&[RF_per_epoch_cue.significant_gain]);
    %     sugainindex(4,:)=suindexes(2,:)&ANOVA_position_effect&[RF_per_epoch_cue.significant_gain];
    %     sugainindex(5,:)=suindexes(3,:)&~(ANOVA_position_effect&[RF_per_epoch_cue.significant_gain]);
    %     sugainindex(6,:)=suindexes(3,:)&ANOVA_position_effect&[RF_per_epoch_cue.significant_gain];
    sugainindex(3,:)=suindexes(3,:)&~(ANOVA_position_effect&[RF_per_epoch_cue.significant_gain]);
    sugainindex(4,:)=suindexes(3,:)&ANOVA_position_effect&[RF_per_epoch_cue.significant_gain];
    
    
    plot_title=[keys.tt.tasktypes{:} ' FR SU_MU with gain field units'];
    filename=plot_title;
    figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
    histoindex=sugainindex;
    histoindex(end+1,:)=true([1 size(histoindex,2)]);
    hold on
    
    cols=[0 0 128; 0 0 255; 128 0 0; 255 0 0]/255;
    colormap(cols);
    for h=1:size(histoindex,1)
        if h==1 || h==3
            histo(:,h)=hist(FR(histoindex(h,:) | histoindex(h+1,:)),bins);
            means(h)=nanmean(FR(histoindex(h,:) | histoindex(h+1,:)));
            medians(h)=nanmedian(FR(histoindex(h,:) | histoindex(h+1,:)));
            stds(h)=nanstd(FR(histoindex(h,:) | histoindex(h+1,:)));
        else
            histo(:,h)=hist(FR(histoindex(h,:)),bins);
            means(h)=nanmean(FR(histoindex(h,:)));
            medians(h)=nanmedian(FR(histoindex(h,:)));
            stds(h)=nanstd(FR(histoindex(h,:)));
        end
        Ns(h)=sum(histoindex(h,:));
    end
    bar(bins,histo(:,1:end-1),'stacked');
    cols=[cols; 0 0 0];
    y_lim=get(gca,'ylim');
    x_lim=get(gca,'xlim');
    for h=1:size(histoindex,1)
        texttoplot=['M:' num2str(means(h)) ' +' num2str(stds(h)) ', med:'  num2str(medians(h)) 'N: ' num2str(Ns(h))];
        text(x_lim(1)+diff(x_lim)/20,y_lim(2)-diff(y_lim)/40*h,texttoplot,'color',cols(h,:));
    end
    
    %legend({'SU','SU+gain','six','six+gain','MU','MU+gain'})
    legend({'SU','SU+gain','MU','MU+gain'})
    xlabel('Firing rate')
    ylabel('N units');
    ph_title_and_save(figure_handle,filename,plot_title,keys)
    
    %% figure 4 number of trials
    
    N_trials_indx=DAG_find_column_index(tuning_per_unit_table,['in_AH_' 'trials_per_condition_' keys.tt.tasktypes{:} ]);
    N_per_condition=[tuning_per_unit_table{2:end,N_trials_indx}];
    n_total=arrayfun(@(x) sum([x.trial.accepted] & [x.trial.success]),population);
    
    plot_title=[keys.tt.tasktypes{:} ' number of trials'];
    filename=plot_title;
    figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
    
    subplot(3,1,1)
    bins=0:1:max(N_per_condition);
    hist(N_per_condition,bins)
    y_lim=get(gca,'ylim');
    x_lim=get(gca,'xlim');
    
    texttoplot=['M:' num2str(mean(N_per_condition)) ' +' num2str(std(N_per_condition)) ', med:'  num2str(median(N_per_condition)) ', range:'  num2str(min(N_per_condition)) ' to ' num2str(max(N_per_condition))];
    text(x_lim(1)+diff(x_lim)/20,y_lim(2)-diff(y_lim)/40,texttoplot);
    xlabel('MMinimum trials per condition')
    ylabel('N units');
    
    subplot(3,1,2)
    bins=0:5:ceil(max(n_total(units))*10)/10;
    hist(n_total(units),bins)
    y_lim=get(gca,'ylim');
    x_lim=get(gca,'xlim');
    texttoplot=['M:' num2str(mean(n_total(units))) ' +' num2str(std(n_total(units))) ', med:'  num2str(median(n_total(units))) ', range:'  num2str(min(n_total(units))) ' to ' num2str(max(n_total(units)))];
    text(x_lim(1)+diff(x_lim)/20,y_lim(2)-diff(y_lim)/40,texttoplot);
    xlabel('Total trials per condition')
    ylabel('N units');
    
    
    n_total=n_total/24;
    subplot(3,1,3)
    bins=0:1:ceil(max(n_total(units)));
    hist(n_total(units),bins)
    y_lim=get(gca,'ylim');
    x_lim=get(gca,'xlim');
    texttoplot=['M:' num2str(mean(n_total(units))) ' +' num2str(std(n_total(units))) ', med:'  num2str(median(n_total(units))) ', range:'  num2str(min(n_total(units))) ' to ' num2str(max(n_total(units)))];
    text(x_lim(1)+diff(x_lim)/20,y_lim(2)-diff(y_lim)/40,texttoplot);
    xlabel('Average trials per condition')
    ylabel('N units');
    ph_title_and_save(figure_handle,filename,plot_title,keys)
end


end