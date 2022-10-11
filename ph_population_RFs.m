function ph_population_RFs(population,modified_keys)
warning('off','MATLAB:catenate:DimensionMismatch');

for fn=fieldnames(modified_keys)'
    keys.(fn{:})=modified_keys.(fn{:});
end

keys.PO.FR_subtract_baseline=~strcmp(keys.PO.epoch_BL,'none'); %% this one we should not need here any more !?

%% gaussian fit settings
fitsettings.sd_max_x=12; 
fitsettings.sd_x_min_ratio=0.125;
%fitsettings.sd_max_y=8; %fitsettings.sd_max_x; %8?
fitsettings.sd_xy_min_ratio=0.25;
%fitsettings.sd_xy_max_ratio=1;
%fitsettings.sd_y_min_ratio=fitsettings.sd_x_min_ratio;
fitsettings.xout=[-30:30]; 
fitsettings.yout=[-15:15]; 
%fitsettings.range_factor=1;
fitsettings.fittypes=keys.PO.fittypes;
fitsettings.baseline_subtracted=keys.PO.FR_subtract_baseline;
keys.PO.fitsettings=fitsettings;

%% tuning table preparation and grouping
[tuning_per_unit_table]                 = ph_load_extended_tuning_table(keys);
[tuning_per_unit_table, Sel_for_title]  = ph_reduce_tuning_table(tuning_per_unit_table,keys);
idx_group_parameter=DAG_find_column_index(tuning_per_unit_table,keys.PO.group_parameter);
idx_unitID=DAG_find_column_index(tuning_per_unit_table,'unit_ID');
idx_RF_frame=DAG_find_column_index(tuning_per_unit_table,keys.PO.RF_frame_parameter);
group_values=tuning_per_unit_table(:,idx_group_parameter);
group_values=cellfun(@num2str, group_values, 'UniformOutput', false);
cell_in_any_group=[false; ~ismember(group_values(2:end),keys.PO.group_excluded)];
unique_group_values=unique(group_values(cell_in_any_group));
if isempty(unique_group_values)
    disp('no relevant groups found');
    return;
end
tuning_per_unit_table=tuning_per_unit_table(cell_in_any_group,:);
group_values=group_values(cell_in_any_group);
complete_unit_list={population.unit_ID}';
population=population(ismember(complete_unit_list,tuning_per_unit_table(:,idx_unitID)));

%population_all_tr = ph_accept_trials_per_unit(population,keys);
all_trialz=[population.trial];
%all_trialz2=[population_all_tr.trial];

[UC, CM, labels]=ph_get_condition_matrix(all_trialz,keys);

%% fix labels --> and with labels colors!
legend_labels_hem={};
legend_labels_prf={};
legend_labels_pos={};
for h=UC.hemifield % append hemifield labels, careful with the order!
    legend_labels_hem=[legend_labels_hem; strcat(labels,['_' keys.labels.hemifield{h+2}])];
end
for h=1:2          % append preference labels, careful with the order!
    legend_labels_prf=[legend_labels_prf; strcat(labels,['_' keys.labels.preferred{h}])];
end
for p=1:size(UC.position,1)
    legend_labels_pos=[legend_labels_pos; labels];
end
% repeat for all effectors
legend_labels_hem=repmat(reshape(legend_labels_hem,numel(legend_labels_hem),1),numel(UC.effector),1);
legend_labels_prf=repmat(reshape(legend_labels_prf,numel(legend_labels_prf),1),numel(UC.effector),1);
legend_labels_pos=repmat(reshape(legend_labels_pos,numel(legend_labels_pos),1),numel(UC.effector),1);

complete_unit_list={population.unit_ID}';
unit_valid=ismember(tuning_per_unit_table(:,idx_unitID),complete_unit_list);
group_values=group_values(unit_valid);
tuning_per_unit_table=tuning_per_unit_table(unit_valid,:);

conditions_out              = combvec(UC.effector,CM')';
conditions_hf               = combvec(UC.hemifield,conditions_out')';
conditions_pref             = combvec([0 1],conditions_out')';

if isfield(UC,'reach_hand') && any(UC.reach_hand~=0) && any(UC.perturbation==1) %splitting to all 4 hand space conditions if hands are involved
    [~,~,columns_hf] = unique(conditions_hf(:,[1,3]),'rows');
    [~,~,columns_pref] = unique(conditions_pref(:,[1,3]),'rows');
else
    columns_hf          = ones(size(conditions_hf,1),1);
    columns_pref        = ones(size(conditions_pref,1),1);
end


%% normalization
keys.normalization_field='RF';
[~, condition,~,pref_valid]=ph_condition_normalization(population,keys,UC,CM);

%% plots

logscale_127=(1:127)'*255/127;
if keys.PO.FR_subtract_baseline
    RF_colormap=[logscale_127 logscale_127 ones(127,1)*255; 255 255 255; ones(127,1)*255 flipud(logscale_127) flipud(logscale_127)]/256;
else
    RF_colormap=[255*ones(1,255);255:-1:1;255:-1:1]'/255;
end

for t=1:size(condition,1)
    typ=UC.type(mod(t-1,numel(UC.type))+1);
    if strcmp(keys.PO.normalization,'percent_change')
        normalization_label=sprintf('N_prct %s to %s',keys.PO.epoch_BL,keys.PO.epoch_DN);
    elseif ~keys.PO.FR_subtract_baseline
        normalization_label=sprintf('N_%s in %s',keys.PO.normalization,keys.PO.epoch_DN);
    elseif keys.PO.baseline_per_trial
        normalization_label=sprintf('N_%s in %s, - %s per trial',keys.PO.normalization,keys.PO.epoch_DN,keys.PO.epoch_BL);
    else
        normalization_label=sprintf('N_%s (X-%s)by%s',keys.PO.normalization,keys.PO.epoch_BL,keys.PO.epoch_DN);
    end
     %% ideally, the names are taken from conditions as well, not hardcoded tt.hands and tt.choices !!
    condition_title=ph_get_condition_title(keys);
    fig_title=sprintf('%s %s %s %s %s %s grouped by %s ',...
        keys.monkey,[keys.conditions_to_plot{:}],keys.arrangement,condition_title,[Sel_for_title{:}],normalization_label,keys.PO.group_parameter);
    filename=sprintf('%s %s %s %s %s %s %s ',...
        keys.monkey,[keys.conditions_to_plot{:}],keys.arrangement(1:3),condition_title,[Sel_for_title{:}],normalization_label,keys.PO.group_parameter);
    
    %save metadata
    unit_IDs=complete_unit_list;
    save([keys.basepath_to_save, keys.project_version, filesep, 'population_meta_data', filesep, 'RF_analysis', filesep, filename], 'keys','tuning_per_unit_table','unit_IDs');
    
    %% plots  
    unique_group_values_tmp=unique_group_values;
    for gt=1:numel(unique_group_values_tmp)
        unique_group_values=unique_group_values_tmp(gt);
        
        unitidx=ismember(complete_unit_list,tuning_per_unit_table(ismember(group_values,unique_group_values_tmp(gt)),idx_unitID));
        group_units=find(all(unitidx,2))';
            
        for c=1:size(condition,2)
            clear FRmax maxposition
            %% typ we know already, this is only for naming!
            eff=conditions_out(c,1); % this here be wrong
            
            per_unit=[condition(t,c).fitting.unit(group_units)];
            pos=[per_unit.positions];
            
            %normalize firing rates?
            for u=1:size(pos,2)
                FR=[pos(:,u).FR];
                if keys.PO.FR_subtract_baseline
                    [FRmax(u), maxposition(u)]=max([FR(FR>0) 0]);
                    FRmin=min([FR(FR<0) 0]);
                    FRmax(u)=max([abs(FRmax(u)) abs(FRmin)]);
                    FR255=num2cell(round((FR+FRmax(u))/2/FRmax(u)*254)+1);
                else
                    [FRmax(u), maxposition(u)]=max(FR);
                    FR255=num2cell(round(FR/FRmax(u)*254)+1);
                end
                [pos(:,u).FR255_GAU]=deal(FR255{:});
                [pos(:,u).FR255]=deal(FR255{:});
            end
            
            cond(c).FRmax=FRmax;
            cond(c).pos=pos;
            cond(c).maxposition=maxposition;
            
            conditiontitle{c}=['T' num2str(typ) 'E' num2str(eff) ' ' labels{mod(c-1,size(conditions_out,1)/numel(UC.effector))+1}];
        end
        
        
        if   keys.PO.plot_per_position            
            %% FR summary plot
            plot_title_part        = ['=' unique_group_values_tmp{1} ' FR in ' keys.PO.epoch_RF ' average' ];
            f_handle              = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',[fig_title plot_title_part]);
            for c=1:size(condition,2)
                subplot(ceil(sqrt(size(condition,2))),ceil(sqrt(size(condition,2))),c)
                title(strrep(conditiontitle{c},'_',' '));
                hold on;
                colormap(RF_colormap);
                caxis([1 255]);
                pos=cond(c).pos;
                FRmax=cond(c).FRmax;
                for p=1:size(pos,1)
                    FRmean(p)=nanmean([pos(p,:).FR]);
                end
                for p=1:size(pos,1)
                    if ~keys.PO.FR_subtract_baseline %% these 2 lines were reversed?
                        FRmeancolidx=round(FRmean(p)/mean(FRmax)*254)+1;
                    else
                        FRmeancolidx=round(FRmean(p)/mean(FRmax)/2*254 + 127)+1;
                    end
                    if ~isnan(FRmean(p))
                        plot(pos(p,1).x,pos(p,1).y,'o','markerfacecolor',RF_colormap(FRmeancolidx,:),'markeredgecolor','none','markersize',5);
                        text(pos(p,1).x,pos(p,1).y,num2str(round(FRmean(p)*100)/100));
                    end
                end
                axis equal
                colorbar;
                set(gca,'Xtick',[],'Ytick',[],'xlim',[min(fitsettings.xout) max(fitsettings.xout)],'ylim',[min(fitsettings.yout) max(fitsettings.yout)]);
            end
            ph_title_and_save(f_handle,  [filename plot_title_part],[fig_title plot_title_part],keys);
            
            %% FR peak histogram plot
            plot_title_part        = ['=' unique_group_values_tmp{1} ' FR in ' keys.PO.epoch_RF ' histogram' ];
            f_handle              = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',[fig_title plot_title_part]);
            for c=1:size(condition,2)
                subplot(ceil(sqrt(size(condition,2))),ceil(sqrt(size(condition,2))),c)
                title(strrep(conditiontitle{c},'_',' '));
                hold on;
                colormap(RF_colormap);
                caxis([1 255]);
                pos=cond(c).pos;
                maxposition=cond(c).maxposition;
                for p=1:size(pos,1)
                    N_per_pos(p)=sum(maxposition==p);
                end
                for p=1:size(pos,1)
                    FRmean=round(N_per_pos(p)/max(N_per_pos)*254)+1;
                    if ~isnan(FRmean)
                        plot(pos(p,1).x,pos(p,1).y,'o','markerfacecolor',RF_colormap(FRmean,:),'markeredgecolor','none','markersize',5);
                        text(pos(p,1).x,pos(p,1).y,num2str(N_per_pos(p)))
                    end
                end
                axis equal
                colorbar;
                set(gca,'Xtick',[],'Ytick',[],'xlim',[min(fitsettings.xout) max(fitsettings.xout)],'ylim',[min(fitsettings.yout) max(fitsettings.yout)]);
            end
            ph_title_and_save(f_handle,  [filename plot_title_part],[fig_title plot_title_part],keys);
        end
        
        
        %% RF and FR plots
        if ~strcmp(keys.PO.epoch_RF,'none')
            angles=[0:pi/100:(2+1/100)*pi];
            circle_x=cos(angles);
            circle_y=sin(angles);
            g=1;
            % for g=1:numel(unique_group_values)
            unitidx=ismember(complete_unit_list,tuning_per_unit_table(ismember(group_values,unique_group_values(g)),idx_unitID));
            group_units=find(all(unitidx,2))';
            if isempty(keys.PO.RF_columns) || isempty(keys.PO.RF_rows)
                RF_columns=ceil(sqrt(numel(group_units)+1));
                RF_rows=ceil(sqrt(numel(group_units)+1));
            else
                RF_columns=keys.PO.RF_columns;
                RF_rows=keys.PO.RF_rows;
            end
            [~,complete_list_tuning_table_idx]= ismember(complete_unit_list,tuning_per_unit_table(:,idx_unitID));
            RF_frame_entries=tuning_per_unit_table(complete_list_tuning_table_idx(group_units),idx_RF_frame);
            for c=1:size(condition,2)
                per_unit=[condition(t,c).fitting.unit(group_units)];
                RFparameters            =[per_unit.parameters];
                %pos=[per_unit.positions];
                                
                FRmax=cond(c).FRmax;
                pos= cond(c).pos;
                maxposition=cond(c).maxposition;
                
                bestfits={RFparameters.bestfit};
                secondbestfits={RFparameters.secondbestfit};
                
                [~,~,idx_fittype]=unique({RFparameters.bestfit});
                idx_fittype(ismember(bestfits,'none'))=100;
                idx_within_fittype=inf(size(idx_fittype));
                
                fittypes=keys.PO.fittypes;
                fittypes_all=[fittypes,'none'];
                RFsizes=zeros(size(RFparameters));
                for f=1:numel(fittypes_all)
                    fittype=fittypes_all{f};
                    fitidx.(fittype)=ismember(bestfits,fittype);
                    fitidx2nd.(fittype)=ismember(secondbestfits,fittype);
                    
                    Parameters_temp=[RFparameters(fitidx.(fittype)).(fittype)];
                    if isempty(Parameters_temp)
                        continue;
                    end
                    switch fittype
                        case {'none'}
                            SC=1:numel(Parameters_temp);
                        case {'sigmoidal','linear'}
                            [~, SC] =sort([Parameters_temp.phi]);
                        case {'gaussian1'}
                            RFsizes(fitidx.(fittype))                =2*2*sqrt([Parameters_temp.sx].*[Parameters_temp.sy]);
                            [~, SC] =sort(RFsizes(fitidx.(fittype)));
                            SC      =fliplr(SC);
                            SC=SC+1000*[Parameters_temp.zmax]>0;
                        case {'gaussian2','gaussian15'}
                            Zmax                     =[Parameters_temp.zmax1; Parameters_temp.zmax2];
                            RFsizes2D                  =2*[2*sqrt([Parameters_temp.sx1].*[Parameters_temp.sy1]);2*sqrt([Parameters_temp.sx2].*[Parameters_temp.sy2])];
                            [~,zone_idx]                =max(RFsizes2D,[],1);  %%
                            zone_idx=sub2ind(size(RFsizes2D),zone_idx,1:size(RFsizes2D,2));
                            RFsizes(fitidx.(fittype))                =RFsizes2D(zone_idx);
                            [~, SC] =sort(RFsizes(fitidx.(fittype)));
                            SC      =fliplr(SC);
                            SC=SC+1000*Zmax(zone_idx)>0;
                    end
                    if ~isempty(idx_RF_frame)
                        [~,SC]=ismember(RF_frame_entries(fitidx.(fittype)),keys.PO.RF_frame_entries);
                        switch fittype
                            case {'gaussian1'}
                                SC=SC+1000*([Parameters_temp.zmax]'>0);
                        end
                    end
                    idx_within_fittype(fitidx.(fittype))=SC;
                end
                [~, sort_by_size_index] =sort(RFsizes);
                sort_by_size_index      =fliplr(sort_by_size_index);
                RF_sorting_matrix       =[idx_fittype;idx_within_fittype]';
                [~, RF_sort_index]      = sortrows(RF_sorting_matrix);
                [~, RF_sort_index]      = sort(RF_sort_index);
                
                
                %% RF plot
                plot_title_part        = ['=' unique_group_values{g} ' con' num2str(c) ' RF in ' keys.PO.epoch_RF];
                f_handle              = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',[fig_title plot_title_part]);
                colormap(RF_colormap);                
                for u=1:numel(group_units)
                    n=RF_sort_index(u);
                    subplot(RF_rows,RF_columns,n);
                    RF=RFparameters(u);
                    if keys.PO.FR_subtract_baseline
                        Zout=round(RF.Zout/FRmax(u)/2*254 + 127)+1;
                    else
                        Zout=round(RF.Zout/FRmax(u)*254)+1;
                    end
                    Zout(Zout>255)=255;
                    Zout(Zout<1)=1;
                    image(fitsettings.xout,-fitsettings.yout,rot90(nanmean(cat(3,Zout),3)))
                    caxis([1 255]);
                    
                    % two ellipses
                    hold on
                    BF=RF.(RF.bestfit);
                    switch RF.bestfit
                        case 'gaussian1'
                            center=[BF.xmax BF.ymax];
                            ellipse_r=2*BF.sx*BF.sy./sqrt(BF.sy.^2*cos(angles - BF.phi).^2+BF.sx.^2*sin(angles  - BF.phi).^2);
                            ellipse_x = circle_x.*ellipse_r; ellipse_y = circle_y.*ellipse_r;
                            line(ellipse_x+center(1),ellipse_y+center(2),'linewidth',1,'color','k');
                        case {'gaussian2','gaussian15'}
                            center=[BF.xmax1 BF.ymax1];
                            ellipse_r=2*BF.sx1*BF.sy1./sqrt(BF.sy1.^2*cos(angles - BF.phi1).^2+BF.sx1.^2*sin(angles  - BF.phi1).^2);
                            ellipse_x = circle_x.*ellipse_r; ellipse_y = circle_y.*ellipse_r;
                            line(ellipse_x+center(1),ellipse_y+center(2),'linewidth',1,'color','k');
                            center=[BF.xmax2 BF.ymax2];
                            ellipse_r=2*BF.sx2*BF.sy2./sqrt(BF.sy2.^2*cos(angles - BF.phi2).^2+BF.sx2.^2*sin(angles  - BF.phi2).^2);
                            ellipse_x = circle_x.*ellipse_r; ellipse_y = circle_y.*ellipse_r;
                            line(ellipse_x+center(1),ellipse_y+center(2),'linewidth',1,'color','k');
                    end
                    
                    if ~isempty(idx_RF_frame)
                        range_x=max(fitsettings.xout)-min(fitsettings.xout);
                        range_y=max(fitsettings.yout)-min(fitsettings.yout);
                        rh=rectangle('Position',[min(fitsettings.xout)+range_x/100 min(fitsettings.yout)+range_y/100  range_x*98/100 range_y*98/100]);
                        set(rh,'edgecolor',keys.PO.RF_frame_colors{ismember(keys.PO.RF_frame_entries,RF_frame_entries(u))}/256);
                        title([population(group_units(u)).unit_ID ' R2= ' num2str(round(RF.R2_adjusted*100)/100) ' ' RF_frame_entries{u}],'fontsize',3,'interpreter','none');
                    else
                        title([population(group_units(u)).unit_ID ' R2= ' num2str(round(RF.R2_adjusted*100)/100) ],'fontsize',3,'interpreter','none');
                    end
                    axis equal
                    set(gca,'Ydir','normal','Xtick',[],'Ytick',[],'xlim',[min(fitsettings.xout) max(fitsettings.xout)],'ylim',[min(fitsettings.yout) max(fitsettings.yout)]);
                end
                subplot(RF_rows,RF_columns,numel(group_units)+1);
                colorbar;
                ph_title_and_save(f_handle,  [filename plot_title_part],[fig_title plot_title_part],keys);
                
                %% fittype R2 values histogram
                plot_title_part       = ['=' unique_group_values{g} ' con' num2str(c) ' RF R2 ' 'in ' keys.PO.epoch_RF];
                f_handle              = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',[fig_title plot_title_part]);
                R2adjusted=[];
                for f=1:numel(fittypes)
                    fittype=fittypes{f};
                    fittemp=[RFparameters.(fittype)];
                    R2adjusted_temp=hist([fittemp.R2],bins);
                    R2adjusted=[R2adjusted;R2adjusted_temp];
                end
                bar(bins,R2adjusted','stacked');
                legend(fittypes);
                ph_title_and_save(f_handle,  [filename plot_title_part],[fig_title plot_title_part],keys);
                
                %% fittype R2 values histogram
                plot_title_part       = ['=' unique_group_values{g} ' con' num2str(c) ' RF R2 adjusted ' 'in ' keys.PO.epoch_RF];
                f_handle              = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',[fig_title plot_title_part]);
                bins=-1:0.05:1;
                R2adjusted=[];
                for f=1:numel(fittypes)
                    fittype=fittypes{f};
                    fittemp=[RFparameters.(fittype)];
                    R2adjusted_temp=hist([fittemp.R2_adjusted],bins);
                    R2adjusted=[R2adjusted;R2adjusted_temp];
                end
                bar(bins,R2adjusted','stacked');
                legend(fittypes);
                ph_title_and_save(f_handle,  [filename plot_title_part],[fig_title plot_title_part],keys);
                
                
                %% fittype R2 values histogram - indicating winner
                plot_title_part       = ['=' unique_group_values{g} ' con' num2str(c) ' RF R2 adjusted win ' 'in ' keys.PO.epoch_RF];
                f_handle              = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',[fig_title plot_title_part]);
                bins=-0.05:0.05:2;
                cols=colormap(jet(numel(fittypes)));
                for f=1:numel(fittypes)
                    fittype=fittypes{f};
                    subplot(numel(fittypes),1,f)
                    hold on;
                    title([fittype, ', N= ' num2str(sum(fitidx.(fittype)))]);
                    R2adjusted=[];
                    fittemp=[RFparameters.(fittype)];
                    fittemp=fittemp((fitidx.(fittype)));
                    for f2=1:numel(fittypes)
                        fittype2=fittypes{f2};
                        fittemp2=[RFparameters.(fittype2)];
                        fittemp2=fittemp2((fitidx.(fittype)));
                        tmp_idx=fitidx2nd.(fittype2);
                        tmp_idx=tmp_idx(fitidx.(fittype));
                        R2_differences=[fittemp(tmp_idx).R2_adjusted]-[fittemp2(tmp_idx).R2_adjusted];
                        R2adjusted_temp=hist(R2_differences,bins);
                        if f==f2
                            R2adjusted_temp=zeros(size(R2adjusted_temp));
                        end
                        R2adjusted=[R2adjusted;R2adjusted_temp];
                        mean_temp(f2)=mean(R2_differences);
                        median_temp(f2)=median(R2_differences);
                        sem_temp(f2)=sterr(R2_differences);
                    end
                    bar(bins,R2adjusted','stacked');
                    y_lim=get(gca,'ylim');
                    for f2=1:numel(fittypes)
                        plot([mean_temp(f2) mean_temp(f2)],y_lim,'color',cols(f2,:));
                        plot([mean_temp(f2)+sem_temp(f2) mean_temp(f2)+sem_temp(f2)],y_lim,':','color',cols(f2,:));
                        plot([mean_temp(f2)-sem_temp(f2) mean_temp(f2)-sem_temp(f2)],y_lim,':','color',cols(f2,:));
                        text(mean_temp(f2),y_lim(1),['mean: ' num2str(round(mean_temp(f2)*100)/100) ' + ' num2str(round(sem_temp(f2)*100)/100)],'rotation',90,'color',cols(f2,:));
                    end
                    legend(fittypes);
                end
                ph_title_and_save(f_handle,  [filename plot_title_part],[fig_title plot_title_part],keys);
                
                
                %% FR plot
                plot_title_part       = ['=' unique_group_values{g} ' con' num2str(c) ' FR ' 'in ' keys.PO.epoch_RF];
                f_handle              = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',[fig_title plot_title_part]);
                colormap(RF_colormap);
                for u=1:numel(group_units)
                    n=RF_sort_index(u);
                    FR255=[pos(:,u).FR255];
                    subplot(RF_rows,RF_columns,n);
                    title(population(group_units(u)).unit_ID,'fontsize',3,'interpreter','none');
                    hold on;
                    scatter([pos(:,u).x],[pos(:,u).y],25,FR255,'filled');
                    caxis([1 255]);
                    axis equal
                    sp_position=get(gca,'Position');%sp_position(3)=sp_position(3)*1.4;sp_position(4)=sp_position(4)*1.4;
                    set(gca,'Position',sp_position,'Ydir','normal','Xtick',[],'Ytick',[],'xlim',[min(fitsettings.xout) max(fitsettings.xout)],'ylim',[min(fitsettings.yout) max(fitsettings.yout)]);
                end
                subplot(RF_rows,RF_columns,numel(group_units)+1);
                colorbar
                ph_title_and_save(f_handle,  [filename plot_title_part],[fig_title plot_title_part],keys);
                
                %% RF centers
                plot_title_part        = ['=' unique_group_values{g} ' con' num2str(c) ' RF  in ' keys.PO.epoch_RF  ' summary'];
                f_handle              = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',[fig_title plot_title_part]);
                %sorting by RF size (to plot small ones on top of large ones)
                for u=group_units(sort_by_size_index)
                    if RFsizes(group_units==u)==0
                        continue
                    end
                    hold on;
                    Parameters_temp=[RFparameters(group_units==u).(RFparameters(group_units==u).bestfit)];
                    for ellipsn=1:2
                        switch RFparameters(group_units==u).bestfit
                            case 'gaussian1'
                                if ellipsn==2
                                    continue
                                end
                                sx=Parameters_temp.sx;
                                sy=Parameters_temp.sy;
                                xmax=Parameters_temp.xmax;
                                ymax=Parameters_temp.ymax;
                                zmax=Parameters_temp.zmax;
                            case {'gaussian2','gaussian15'}
                                sx=[Parameters_temp.(['sx' num2str(ellipsn)])];
                                sy=[Parameters_temp.(['sy' num2str(ellipsn)])];
                                xmax=Parameters_temp.(['xmax' num2str(ellipsn)]);
                                ymax=Parameters_temp.(['ymax' num2str(ellipsn)]);
                                zmax=Parameters_temp.(['zmax' num2str(ellipsn)]);
                        end
                        RFsize              =2*sqrt(sx.*sy);
                        center              =[xmax ymax];
                        Allmonkeys={'Linus','Curius','Cornelius'};
                        current_monkey=Allmonkeys{cellfun(@(x) any(strfind(x,population(u).unit_ID(1:3))),Allmonkeys)};
                        
                        if sign(zmax)==1
                            col='r';
                        elseif sign(zmax)==-1
                            col='b';
                        else
                            col='k';
                        end
                        monkey_marker=keys.(current_monkey).marker;
                        RF_size_factor=0.1;
                        Radius=RF_size_factor*RFsize;
                        if strcmp(monkey_marker,'o')
                            ellipse_x = circle_x.*Radius;
                            ellipse_y = circle_y.*Radius;
                            line(ellipse_x+center(1),ellipse_y+center(2),'color',col,'linewidth',4);
                        elseif strcmp(monkey_marker,'s')
                            square_x = [-1,-1,1,1,-1].*Radius*sqrt(pi)/2;
                            square_y = [-1,1,1,-1,-1].*Radius*sqrt(pi)/2;
                            line(square_x+center(1),square_y+center(2),'color',col,'linewidth',4);
                        end
                    end
                end
                axis equal
                set(gca,'Ydir','normal','xlim',[min(fitsettings.xout) max(fitsettings.xout)],'ylim',[min(fitsettings.yout) max(fitsettings.yout)]);
                ph_title_and_save(f_handle,  [filename plot_title_part],[fig_title plot_title_part],keys);
                
                %% RF sizes
                plot_title_part         = ['=' unique_group_values{g} ' con' num2str(c) ' RF sizes in ' keys.PO.epoch_RF];
                f_handle              = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',[fig_title plot_title_part]);
                %sorting by RF size (to plot small ones on top of large ones)
                RFsizes_valid=RFsizes(RFsizes~=0);
                minsiz=floor(min(RFsizes_valid));
                maxsiz=ceil(max(RFsizes_valid));
                sizebins=minsiz:(maxsiz-minsiz)/20:maxsiz;
                FRhist=hist(RFsizes_valid,sizebins);
                hold on
                if any(FRhist~=0)
                    bar(sizebins,FRhist);
                    text(sizebins(2), max(FRhist)-max(FRhist/10),['u=' num2str(nanmean(RFsizes_valid)) ', med=' num2str(nanmedian(RFsizes_valid)) ', std=' num2str(nanstd(RFsizes_valid))]);
                end
                ph_title_and_save(f_handle,  [filename plot_title_part],[fig_title plot_title_part],keys);
                
                %% these two vary from condition to condition and should be saved independently for each condition
                
                [~, type_effector_short] = MPA_get_type_effector_name(typ,eff);
                ph_append_to_anova_table(keys,'RFs',complete_unit_list(group_units),RFparameters,type_effector_short,labels{c});
                
                save([keys.path_to_save filename plot_title_part '.mat'],'RFsizes','RFparameters');
            end
            % end
            
        end
    end
    unique_group_values=unique_group_values_tmp;
end

end

