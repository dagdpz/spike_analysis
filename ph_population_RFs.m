function ph_population_RFs(population,modified_keys)
warning('off','MATLAB:catenate:DimensionMismatch');

for fn=fieldnames(modified_keys)'
    keys.(fn{:})=modified_keys.(fn{:});
end

keys.RF.FR_subtract_baseline=~strcmp(keys.RF.epoch_BL,'none'); %% this one we should not need here any more !?

%% gaussian fit settings
fitsettings.sd_max_x=12;
fitsettings.sd_x_min_ratio=0.125;
fitsettings.sd_xy_min_ratio=0.25;
fitsettings.xout=[-30:30];
fitsettings.yout=[-15:15];
fitsettings.fittypes=keys.RF.fittypes;
fitsettings.baseline_subtracted=keys.RF.FR_subtract_baseline;
keys.RF.fitsettings=fitsettings;

%% tuning table preparation and grouping
% TT=keys.tuning_table;
% 
% idx.group_parameter=DAG_find_column_index(TT,keys.RF.group_parameter);
% idx.unitID=DAG_find_column_index(TT,'unit_ID');
% idx.RF_frame=DAG_find_column_index(TT,keys.RF.RF_frame_parameter);
% group_values=TT(:,idx.group_parameter);
% group_values=cellfun(@num2str, group_values, 'UniformOutput', false);
% cell_in_any_group=[false; ~ismember(group_values(2:end),keys.RF.group_excluded)];
% unique_group_values=unique(group_values(cell_in_any_group));
% if isempty(unique_group_values)
%     disp('no relevant groups found');
%     return;
% end
% TT=TT(cell_in_any_group,:);
% group_values=group_values(cell_in_any_group);

[TT,idx,group_values,unique_group_values]=ph_readout_tuning_table(keys);
if isempty(unique_group_values)
    return;
end

complete_unit_list={population.unit_ID}';
[unit_valid,TM]=ismember(complete_unit_list,TT(:,idx.unitID));
population_RF_frames=TT(:,idx.RF_frame);
population=population(unit_valid);
complete_unit_list={population.unit_ID}';
population_group=group_values(TM(unit_valid));   
population_RF_frames=population_RF_frames(TM(unit_valid));   


% complete_unit_list={population.unit_ID}';
% population=population(ismember(complete_unit_list,TT(:,idx.unitID)));

all_trialz=[population.trial];
[UC, CM, labels]=ph_get_condition_matrix(all_trialz,keys);
% 
% unit_valid=ismember(TT(:,idx.unitID),complete_unit_list);
% group_values=group_values(unit_valid);
% group_values=group_values(unit_valid);
% TT=TT(unit_valid,:);

conditions_out              = combvec(UC.effector,CM')';

%% normalization
keys.normalization_field='RF';
[~, condition]=ph_condition_normalization(population,keys,UC,CM);

%% plots

logscale_127=(1:127)'*255/127;
if keys.RF.FR_subtract_baseline
    RF_colormap=[logscale_127 logscale_127 ones(127,1)*255; 255 255 255; ones(127,1)*255 flipud(logscale_127) flipud(logscale_127)]/256;
else
    RF_colormap=[255*ones(1,255);255:-1:1;255:-1:1]'/255;
end

for t=1:size(condition,1)
    typ=UC.type(mod(t-1,numel(UC.type))+1);
    [fig_title,filename]=ph_figure_titles(keys);
    %save metadata
    unit_IDs=complete_unit_list;
    save([keys.basepath_to_save, keys.project_version, filesep, 'population_meta_data', filesep, 'response fields', filesep, filename], 'keys','TT','unit_IDs');
    
    %% plots
    unique_group_values_tmp=unique_group_values;
    for gt=1:numel(unique_group_values_tmp)
        unique_group_values=unique_group_values_tmp(gt);
        group_units=find(ismember(population_group,unique_group_values));
        RF_frame_entries=population_RF_frames(group_units);
%         unitidx=ismember(complete_unit_list,TT(ismember(group_values,unique_group_values_tmp(gt)),idx.unitID));
%         group_units=find(all(unitidx,2))';
%         
        for c=1:size(condition,2)
            clear FRmax maxposition
            %% typ we know already, this is only for naming!
            eff=conditions_out(c,1); % this here be wrong
            
            
            per_unit=[condition(t,c).fitting.unit(group_units)];
            pos=[per_unit.positions];
            %normalize firing rates?
            for u=1:size(pos,2)
                FR=[pos(:,u).FR];
                if keys.RF.FR_subtract_baseline
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
        
        %% FR summary plot
        plot_title_part        = ['=' unique_group_values_tmp{1} ' FR in ' keys.RF.epoch_RF ' average' ];
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
                if ~keys.RF.FR_subtract_baseline %% these 2 lines were reversed?
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
        plot_title_part        = ['=' unique_group_values_tmp{1} ' FR in ' keys.RF.epoch_RF ' histogram' ];
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
        
        
        %% RF and FR plots
        if strcmp(keys.RF.epoch_RF,'none')
            disp('No valid response field epoch selected');
            continue;
        end    
        angles=[0:pi/100:(2+1/100)*pi];
        circle_x=cos(angles);
        circle_y=sin(angles);
        
%         unique_group_values=unique_group_values_tmp(gt);
%         g_idx=ismember(population_group,unique_group_values(1));
%         g=1;
%         unitidx=ismember(complete_unit_list,TT(ismember(group_values,unique_group_values(g)),idx.unitID));
%         group_units=find(all(unitidx,2))';
        if isempty(keys.RF.RF_columns) || isempty(keys.RF.RF_rows)
            RF_columns=ceil(sqrt(numel(group_units)+1));
            RF_rows=ceil(sqrt(numel(group_units)+1));
        else
            RF_columns=keys.RF.RF_columns;
            RF_rows=keys.RF.RF_rows;
        end
%         [~,complete_list_tuning_table_idx]= ismember(complete_unit_list,TT(:,idx.unitID));
%         RF_frame_entries=TT(complete_list_tuning_table_idx(g_idx),idx.RF_frame);
        
        
        
        
        for c=1:size(condition,2)
            per_unit=[condition(t,c).fitting.unit(group_units)];
            RFparameters            =[per_unit.parameters];
            FRmax=cond(c).FRmax;
            pos= cond(c).pos;
            maxposition=cond(c).maxposition;
            bestfits={RFparameters.bestfit};
            secondbestfits={RFparameters.secondbestfit};
            [~,~,idx_fittype]=unique({RFparameters.bestfit});
            idx_fittype(ismember(bestfits,'none'))=100;
            idx_within_fittype=inf(size(idx_fittype));
            fittypes=keys.RF.fittypes;
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
                if ~isempty(idx.RF_frame)
                    [~,SC]=ismember(RF_frame_entries(fitidx.(fittype)),keys.RF.RF_frame_entries);
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
            plot_title_part        = ['=' unique_group_values{g} ' con' num2str(c) ' RF in ' keys.RF.epoch_RF];
            f_handle              = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',[fig_title plot_title_part]);
            colormap(RF_colormap);
            for u=1:numel(group_units)
                n=RF_sort_index(u);
                subplot(RF_rows,RF_columns,n);
                RF=RFparameters(u);
                if keys.RF.FR_subtract_baseline
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
                
                if ~isempty(idx.RF_frame)
                    range_x=max(fitsettings.xout)-min(fitsettings.xout);
                    range_y=max(fitsettings.yout)-min(fitsettings.yout);
                    rh=rectangle('Position',[min(fitsettings.xout)+range_x/100 min(fitsettings.yout)+range_y/100  range_x*98/100 range_y*98/100]);
                    set(rh,'edgecolor',keys.RF.RF_frame_colors{ismember(keys.RF.RF_frame_entries,RF_frame_entries(u))}/256);
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
            plot_title_part       = ['=' unique_group_values{g} ' con' num2str(c) ' RF R2 ' 'in ' keys.RF.epoch_RF];
            f_handle              = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',[fig_title plot_title_part]);
            bins=-1:0.05:1;
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
            plot_title_part       = ['=' unique_group_values{g} ' con' num2str(c) ' RF R2 adjusted ' 'in ' keys.RF.epoch_RF];
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
            plot_title_part       = ['=' unique_group_values{g} ' con' num2str(c) ' RF R2 adjusted win ' 'in ' keys.RF.epoch_RF];
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
            plot_title_part       = ['=' unique_group_values{g} ' con' num2str(c) ' FR ' 'in ' keys.RF.epoch_RF];
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
            plot_title_part        = ['=' unique_group_values{g} ' con' num2str(c) ' RF  in ' keys.RF.epoch_RF  ' summary'];
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
            plot_title_part         = ['=' unique_group_values{g} ' con' num2str(c) ' RF sizes in ' keys.RF.epoch_RF];
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
        
    end
    unique_group_values=unique_group_values_tmp;
end

end

