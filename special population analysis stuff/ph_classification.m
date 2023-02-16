function population=ph_classification(population,modified_keys)

for fn=fieldnames(modified_keys)'
    keys.(fn{:})=modified_keys.(fn{:});
end

[TT,idx,group_values,unique_group_values]=ph_readout_tuning_table(keys);
if isempty(unique_group_values)
    disp('no relevant groups found');
    return;
end
complete_unit_list={population.unit_ID}';

for g=1:numel(unique_group_values)
    clear suindexes
    unitidx=ismember(complete_unit_list,TT(ismember(group_values,unique_group_values(g)),idx.unitID));
    units=find(all(unitidx,2))';
        
    singlerating=[population(units).Single_rating];
    stabilityrating=[population(units).stability_rating]; % all of them are supposed to be stable
    snrrating=[population(units).SNR_rating];
    
    suindexes(1,:)=singlerating==1 & snrrating<=2;
    suindexes(2,:)=(singlerating>1 | snrrating>=2);

    FR=[population(units).FR];
    bins=2.5:5:ceil(max(FR)*10)/10;

    %% figure 1 ONLY single ranking ,
    
    plot_title=[keys.tt.tasktypes{:}  ' FR dependent on Su ranking'];
    filename=plot_title;
    figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
    histoindex= suindexes ;
    colormap(jet(size(histoindex,1)));
    histoindex(end+1,:)=true([1 size(histoindex,2)]);
    for h=1:size(histoindex,1)
        histo(:,h)=hist(FR(histoindex(h,:)),bins);
        means(h)=nanmean(FR(histoindex(h,:)));
        medians(h)=nanmedian(FR(histoindex(h,:)));
        stds(h)=nanstd(FR(histoindex(h,:)));
    end
    bar(bins,histo(:,1:end-1),'stacked');
    cols=[jet(size(suindexes,1)); 0 0 0];
    y_lim=get(gca,'ylim');
    x_lim=get(gca,'xlim');
    for h=1:size(histoindex,1)
        texttoplot=['Mean:' num2str(means(h)) ' +' num2str(stds(h)) ', med:'  num2str(medians(h)) ];
        text(x_lim(2)-diff(x_lim)/2,y_lim(2)-diff(y_lim)/40*h,texttoplot,'color',cols(h,:));
    end
    legend({'SU','MU'})
    xlabel('Firing rate')
    ylabel('N units');
    ph_title_and_save(figure_handle,filename,plot_title,keys)

    
    %% figure 2 SU/SU+/MU thiss is same as 1 - i guess the idea was to split furter...??
    
    plot_title=[keys.tt.tasktypes{:} ' FR SU_MU'];
    filename=plot_title;
    figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
    histoindex= suindexes ;
    cols=[jet(size(histoindex,1)); 0 0 0];
    histoindex(end+1,:)=true([1 size(histoindex,2)]);
    
    for h=1:size(histoindex,1)
        histo(:,h)=hist(FR(histoindex(h,:)),bins);
        means(h)=nanmean(FR(histoindex(h,:)));
        medians(h)=nanmedian(FR(histoindex(h,:)));
        stds(h)=nanstd(FR(histoindex(h,:)));
    end
    %bar(bins,histo,'stacked');
    for h=1:size(histoindex,1)
        subplot(size(histoindex,1),1,h)
        bar(bins,histo(:,h));
        y_lim=get(gca,'ylim');
        x_lim=get(gca,'xlim');
        texttoplot=['M:' num2str(means(h)) ' +' num2str(stds(h)) ', med:'  num2str(medians(h)) ];
        text(x_lim(1)+diff(x_lim)/20,y_lim(2)-diff(y_lim)/40*h,texttoplot,'color',cols(h,:));
    end
    legend({'SU','MU'})
    xlabel('Firing rate')
    ylabel('N units');
    ph_title_and_save(figure_handle,filename,plot_title,keys)
    
    %% figure 3 SU/MU with gain field units,
    if false ; %exist()
        e=3;
        EPOCHS=keys.EPOCHS_PER_TYPE{3};
        clear sugainindex
        
        popuation_total=population;
        load('Y:\Projects\Pulv_eye_gaze_position\ephys\20190309\tuning_curves\population_FR_vector.mat')
        RF_per_epoch=vertcat(population.RF_per_epoch); % cue
        RF_per_epoch_cue=RF_per_epoch(units,e);
        
        positionanovaindex=DAG_find_column_index(TT,['in_AH_' EPOCHS{e,1} '_position_' keys.tt.tasktypes{:}]);
        gazeanovaindex=DAG_find_column_index(TT,['in_AH_' EPOCHS{e,1} '_fixation_' keys.tt.tasktypes{:}]);
        interactionanovaindex=DAG_find_column_index(TT,['in_AH_' EPOCHS{e,1} '_PxF_' keys.tt.tasktypes{:}]);
        unit_IDs_over_all_position_effect=TT(strcmp(['false'; TT(2:end,positionanovaindex)],'true'),1);
        unit_IDs_gazedependence=TT(strcmp(['false'; TT(2:end,gazeanovaindex)],'true') | strcmp(['false'; TT(2:end,interactionanovaindex)],'true'),1);
        
        ANOVA_position_effect=ismember(IDS,unit_IDs_over_all_position_effect);
        gaze_dependence=ismember(IDS,unit_IDs_gazedependence)';
        %popuation_total=population;
        sugainindex(1,:)=(suindexes(1,:)|suindexes(2,:))&~(ANOVA_position_effect&[RF_per_epoch_cue.significant_gain]);
        sugainindex(2,:)=(suindexes(1,:)|suindexes(2,:))&ANOVA_position_effect&[RF_per_epoch_cue.significant_gain];
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
    end
    
    %% figure 4 number of trials
    
     for criterium={'in','ch'}
        crit=criterium{:};
        switch keys.tt.(['trial_criterion_' crit])
            case 'per_hemifield_and_perturbation'
                strtofind='trials_per_hemifield'; %% ??? - there will be a problem with perturbation stuff !!
                disp(['per_hemifield_and_perturbation not supported as trial criterion' ]);
        end
        
        strtofind=['trials_' keys.tt.(['trial_criterion_' crit])]; %% why this needed to be a cell here but not in load_tuning_table
        N_trial_idx.(crit)=find(~cellfun(@isempty,strfind(keys.tuning_table(1,:),strtofind)));
        N_titles=keys.tuning_table(1,N_trial_idx.(crit));
        tmp_idx=find(~cellfun(@isempty,strfind(N_titles,crit)));
         N_trial_idx.(crit)= N_trial_idx.(crit)(tmp_idx);
     end
    
    instructed_N_titles=keys.tuning_table(1,N_trial_idx.in);
        tmp_idx=find(~cellfun(@isempty,strfind(instructed_N_titles,keys.tt.tasktypes{:})));
    N_per_condition=[TT{:,N_trial_idx.in(tmp_idx)}];
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
    xlabel('Minimum trials per condition')
    ylabel('N units');
    title('Minimum N trials per condition');
    
    subplot(3,1,2)
    bins=0:5:ceil(max(n_total(units))*10)/10;
    hist(n_total(units),bins)
    y_lim=get(gca,'ylim');
    x_lim=get(gca,'xlim');
    texttoplot=['M:' num2str(mean(n_total(units))) ' +' num2str(std(n_total(units))) ', med:'  num2str(median(n_total(units))) ', range:'  num2str(min(n_total(units))) ' to ' num2str(max(n_total(units)))];
    text(x_lim(1)+diff(x_lim)/20,y_lim(2)-diff(y_lim)/40,texttoplot);
    xlabel('Total trials per condition')
    ylabel('N units');
    title('total N trials per condition');
    
    
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
    title('Minimum N trials per condition');
    ph_title_and_save(figure_handle,filename,plot_title,keys)
end


end