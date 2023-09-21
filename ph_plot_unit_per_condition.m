function ph_plot_unit_per_condition(population,trials,keys)
%tuning_per_unit_table=keys.tuning_per_unit_table;
%% colors and ylim preparation
eye_offset  =   -keys.UN.trials_max_for_ylim-keys.UN.excentricity_max_for_ylim*keys.UN.eyetrace_factor;
hnd_offset  =   eye_offset-keys.UN.excentricity_max_for_ylim*keys.UN.eyetrace_factor*2;
eye_col_v   =   keys.colors.eye_ver;
eye_col_h   =   keys.colors.eye_hor;

for u=1:numel(population)
    
    pop=population(u);
    UT=ph_get_unit_trials(pop,trials);
    
    types       =[UT.type];
    effectors   =[UT.effector];
    hands       =[UT.reach_hand];
    current_unit_tuning= [keys.tuning_table(1,:); keys.tuning_table(ismember(keys.tuning_table(:,1), pop.unit_ID),:)];
    if size(current_unit_tuning,1)<2
        continue;
    end
    for type=unique(types)
        
        %% selection of epochs to plot
        effectors_effector_loop=unique(effectors(ismember(effectors,keys.cal.effectors)));
        [keys all_sta]=ph_get_epoch_keys(keys,type,effectors_effector_loop,1);
        n_states=numel(all_sta);
        
        te_idx= types == type & ismember(effectors,effectors_effector_loop);% & ismember(hands,keys.reach_hand);
        eff_for_typ=effectors(te_idx);
        effectors_effector_loop=unique(eff_for_typ);
        
        if sum(te_idx)==0; % this line wont work as it is
            fprintf('%s has no trials for effector %.0f type %.0f hands %s',population(u).unit_ID,effector,type,mat2str(keys.cal.reach_hand));
            continue;
        end
        
        T=UT;
        if strcmp(keys.UN.line_labelling,'contra/ipsi')
            T=ph_LR_to_CI(keys,pop,T); % Convert to ipsi/contra? not sure if necessary here
        end
        nrows=T(1).rows;
        ncolumns=T(1).columns;
        
        %o=o.trial;
        if any(hands(te_idx)>0)
            y_lim_PSTH=hnd_offset-keys.UN.excentricity_max_for_ylim*keys.UN.eyetrace_factor;
        else
            y_lim_PSTH=eye_offset-keys.UN.excentricity_max_for_ylim*keys.UN.eyetrace_factor;
        end
        unique_figures=unique([T.figure]);
        for fig=unique_figures
            fig_idx     =[T.figure]==fig;
            title_part  =unique({T(fig_idx).figure_title_part});
            title_value =unique({T(fig_idx).figure_title_value});
            title_part  =[title_part{:}];
            title_value =[title_value{:}];
            
            [UC, CM, labels]=ph_get_condition_matrix(T(fig_idx),keys);
            
            
            %% defining trials that go into each line fpr current figure
            line_idx=true(size(CM,1),size(fig_idx,2));
            for L=1:size(CM,1)
                for con=1:size(CM,2)
                    line_idx(L,:)=line_idx(L,:) & [T.(keys.condition_parameters{con})]==CM(L,con);
                end
            end
            if keys.plot.average_PSTH_line
                line_idx(end+1,:)=true(1,size(fig_idx,2));
                labels{end+1}='AV';
            end
            n_lines=size(line_idx,1);
            
            
            %% here basically decide about labels dependent on L/R or C/I plotting
            if strcmp(keys.UN.line_labelling,'contra/ipsi')
                legend_labels=labels;
                side_labels={'I','C'};
                side_labels_col={'IS','CS'};
                if any(UC.hemifield==0) %&& fig==unique_figures(1)% vertical targets, KK uncomment, why figure and whz colors?
                    side_labels={'I','V','C'};
                    side_labels_col={'IS','VS','CS'};
                end
            else
                legend_labels=strrep(labels,'IH','LH');
                legend_labels=strrep(legend_labels,'CH','RH');
                side_labels={'L','R'};
                side_labels_col={'IS','CS'};
                if any(UC.hemifield==0) %&& fig==unique_figures(1)% vertical targets, KK uncomment, why figure and whz colors?
                    side_labels={'L','V','R'};
                    side_labels_col={'IS','VS','CS'};
                end
            end
            labels_hem={};
            for h=UC.hemifield % append hemifield labels, careful with the order!
                labels_hem=[labels_hem; strcat(labels,['_' keys.labels.hemifield{h+2}])];
            end
            labels_hem=labels_hem';
            labels_hem=labels_hem(:);
            
            
            effectors_on_figure=effectors_effector_loop(ismember(effectors_effector_loop,eff_for_typ(fig_idx(te_idx))));
            for e=1:numel(effectors_on_figure)
                clear hh hs hb
                effector=effectors_on_figure(e);
                te_idx= te_idx & effectors==effector;
                [type_effector_full, type_effector_short type_string] = MPA_get_type_effector_name(type,effector);
                
                %% ANOVA results to print
                anova_title_part=get_anova_results(keys,current_unit_tuning,type_effector_short,keys.UN.anova_main,'');
                fig_title=sprintf('%s, %s, %s, %s %s %s',population(u).unit_ID, population(u).target, type_effector_full, keys.arrangement, title_part, title_value);
                fig_title_part=sprintf(', Stability %.1f, SNR %.1f, Single %.1f Grid: %d/%d Depth %.2f ANOVA %s Channel: %d Blocks&Units: %s ', ...
                    population(u).avg_stability, population(u).avg_SNR, population(u).avg_single_rating, population(u).grid_x, population(u).grid_y,...
                    population(u).electrode_depth, anova_title_part, population(u).channel, [population(u).block_unit{:}]);
                
                %% Per position PSTH figure
                plot_1_title            = [fig_title  ' PSTHs'];
                PSTH_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',[plot_1_title fig_title_part]);
                subplot_indizes=unique([T(fig_idx & te_idx).subplot_pos]);
                y_max =1;
                yL_min=0;
                for sub=subplot_indizes
                    subplot(nrows,ncolumns,sub,'XlimMode','Manual','YlimMode','Manual','XtickMode','Manual')
                    hold on
                    subplot_struct=T([T.subplot_pos]==sub & te_idx);
                    title([subplot_struct(1).title_part ' ' num2str(round(subplot_struct(1).position))]);
                    state_seperator=zeros(n_lines,1);
                    state_seperator_max=0.2;
                    for w=1:size(keys.PSTH_WINDOWS,1)
                        sta=keys.PSTH_WINDOWS{w,2};
                        t_before_state=keys.PSTH_WINDOWS{w,3};
                        t_after_state=keys.PSTH_WINDOWS{w,4};
                        histo=[];
                        raster_y=0;
                        for lin=1:n_lines
                            col=keys.colors.(labels{lin})/255;
                            ix=line_idx(lin,:) & te_idx & [T.subplot_pos]==sub & [T.figure]==fig;
                            line_struct=T(ix);
                            ATs=pop.trial(ix);
                            n_trials=numel(line_struct);
                            state_shift=state_seperator(lin)-t_before_state;
                            state_seperator(lin)=state_shift + t_after_state + 0.1;
                            if n_trials==0 || ~any([line_struct.states]==sta)
                                continue;
                            end
                            
                            %% PSTH
                            [histo(lin,:), bins, ~, SEM]=ph_spike_density(ATs,line_struct,w,keys,zeros(numel(line_struct),1),ones(numel(line_struct),1));
                            bins=bins+state_shift;
                            lineProps={'color',col,'linewidth',keys.UN.PSTH_perpos_width};
                            shadedErrorBar(bins,histo(lin,:),SEM,lineProps,1);
                            
                            %% skipping average raster
                            if keys.plot.average_PSTH_line && lin==n_lines
                                continue;
                            end
                            
                            %% Raster and eye hand traces  + Number of trials line and text
                            line([state_shift+t_before_state state_shift+t_after_state],[raster_y raster_y],'Color',col,'LineWidth',0.5); hold on;
                            raster=[NaN; NaN];
                            for t=1:numel(line_struct)
                                if line_struct(t).reach_hand==1
                                    hnd_col_v=keys.colors.lhd_ver;
                                    hnd_col_h=keys.colors.lhd_hor;
                                elseif line_struct(t).reach_hand==2
                                    hnd_col_v=keys.colors.rhd_ver;
                                    hnd_col_h=keys.colors.rhd_hor;
                                else
                                    hnd_col_v=[1 1 1];
                                    hnd_col_h=[1 1 1];
                                end
                                trial_state_onset=line_struct(t).states_onset(line_struct(t).states==sta);
                                at_idx=ATs(t).arrival_times-trial_state_onset>=t_before_state &...
                                    ATs(t).arrival_times-trial_state_onset<=t_after_state;
                                %if ~strcmp(keys.UN.line_labelling,'contra/ipsi')
                                time_axis=line_struct(t).time_axis-trial_state_onset+state_shift;
                                t_idx=line_struct(t).time_axis-trial_state_onset>=t_before_state &...
                                    line_struct(t).time_axis-trial_state_onset<=t_after_state;
                                if keys.plot.eye_hand_traces
                                    line(time_axis(t_idx),line_struct(t).x_eye(t_idx)*keys.UN.eyetrace_factor + eye_offset,'color',eye_col_h);
                                    line(time_axis(t_idx),line_struct(t).y_eye(t_idx)*keys.UN.eyetrace_factor + eye_offset,'color',eye_col_v);
                                    line(time_axis(t_idx),line_struct(t).x_hnd(t_idx)*keys.UN.hndtrace_factor + hnd_offset,'color',hnd_col_h);
                                    line(time_axis(t_idx),line_struct(t).y_hnd(t_idx)*keys.UN.hndtrace_factor + hnd_offset,'color',hnd_col_v);
                                end
                                %end
                                AT=[ATs(t).arrival_times(at_idx)]'+state_shift-trial_state_onset;
                                raster=[raster [AT;repmat(raster_y-t*0.1,size(AT))]];
                                %ig_make_raster(AT,raster_y-t,1,0,'Color',col,'LineWidth',keys.UN.raster_width);
                            end
                            scatter(raster(1,:),raster(2,:),10,col,'s','filled');
                            line([state_shift+t_before_state state_shift+t_after_state],[raster_y-t raster_y-t],'Color',col,'LineWidth',0.5);
                            
                            if w==1
                                not_accepted=find(~[line_struct.accepted])*-1+raster_y;
                                plot((state_shift+t_before_state)*ones(size(not_accepted)),not_accepted,'marker','s','markeredgecolor','k','markerfacecolor','k','linestyle','none');
                                text(state_shift+t_before_state,raster_y,['N=' num2str(t)],'Color',col,'fontsize',5,'verticalalignment','top');
                            end
                            raster_y=raster_y-n_trials;
                            
                        end
                        state_seperator_max=max([state_seperator_max;state_seperator]);
                        y_max=max([y_max max(histo)]);
                        yL_min=min([yL_min raster_y]);
                    end
                end
                y_lim=[y_lim_PSTH,y_max];
                for sub=subplot_indizes
                    subplot(nrows,ncolumns,sub,'XlimMode','Manual','YlimMode','Manual','XtickMode','Manual')
                    ph_PSTH_background(T(te_idx),y_lim,[yL_min y_lim(2)],[0 y_lim(2)],keys,0.5)
                    set(gca,'ylim',y_lim,'xlim',[0 state_seperator_max-0.1]);
                end
                ph_title_and_save(PSTH_summary_handle,plot_1_title,[plot_1_title fig_title_part],keys);
                
                %% FR figure
                plot_title        =[fig_title ' FR'];
                FR_summary_handle =figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',[plot_title fig_title_part]);
                tr_idx=fig_idx & te_idx & [T.accepted]==1;
                if sum(tr_idx)==0 % might be the case if all trials are excluded cause of no spikes
                    continue;
                end
                subplot_indizes=1:max([T.subplot_pos]);
                clear FR_heat;
                
                for sub=subplot_indizes
                    for sta=all_sta
                        for lin=1:n_lines
                            FR_idx = line_idx(lin,:) & [T.subplot_pos]==sub & [T.figure]==fig;
                            FR_idx = FR_idx([T.type]==type);
                            if sum(FR_idx)==0
                                FR_heat(sta,lin).pos(sub)    =NaN;
                            else
                                FR_heat(sta,lin).pos(sub)    =nanmean(pop.epochs_per_type{type}(FR_idx,sta));
                            end
                        end
                    end
                end
                
                % Plotting heat maps
                for sta=all_sta
                    for lin=1:n_lines
                        state_label =keys.EPOCHS{sta,1};
                        hh(sta+n_states*(lin-1))=subplot_assignment(keys,'FR',n_states,n_lines,sta,lin,e,0,effectors_on_figure,1);
                        plot_firing_rate_heatmap(FR_heat(sta,lin).pos,ncolumns,nrows,1:max(subplot_indizes));
                        title([strrep(legend_labels{lin},'_',' ') ' ' state_label]);
                    end
                    % one single colorbar
                    if sta==n_states && exist('hh','var')
                        hp4 = get(hh(end),'Position');
                        hc = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.015  hp4(4)*n_lines]);
                        set(get(hc,'title'),'String','Sp/s');
                    end
                end
                ig_set_caxis_equal_lim(hh);
                
                %% Plotting polars
                if keys.plot.polars_on_extra_figure
                    ph_title_and_save(FR_summary_handle,plot_title,[plot_title fig_title_part],keys);
                    plot_title        =[fig_title ' polars'];
                    FR_summary_handle   =figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',[plot_title fig_title_part]);
                end
                clear FR_means FR_sterr;
                positions               =vertcat(T(tr_idx).position);
                position_angles         =angle(positions(:,1)+1i*positions(:,2));
                position_norm           =abs(positions(:,1)+1i*positions(:,2));
                position_norm(position_norm==0)=1;
                norm_vector=positions./[position_norm position_norm];
                [unique_angles]=unique(position_angles);
                
                for sta=all_sta
                    sta_FR=pop.epochs_per_type{type}(tr_idx([T.type]==type),sta);
                    norm_factor=max(sta_FR);
                    sta_vector=[sta_FR.*norm_vector(:,1) sta_FR.*norm_vector(:,2)]/norm_factor;
                    for lin=1:n_lines
                        lin_idx=line_idx(lin,tr_idx);
                        if sum(lin_idx)==0
                            FR_means(sta,lin).vector=[NaN NaN];
                        else
                            FR_means(sta,lin).vector=nanmean(sta_vector(lin_idx,:),1);
                        end
                        for ang=1:numel(unique_angles)
                            states_trials       =pop.epochs_per_type{type}(position_angles==unique_angles(ang) & lin_idx',sta);
                            if isempty(states_trials) %unfortunate debugging for empty positions !
                                FR_means(sta,lin).pos(ang)    =NaN;
                                FR_sterr(sta,lin).pos(ang)    =NaN;
                            else
                                FR_means(sta,lin).pos(ang)    =nanmean(states_trials)/norm_factor;
                                FR_sterr(sta,lin).pos(ang)    =sterr(states_trials)/norm_factor;
                            end
                        end
                    end
                end
                
                FR_mean_max=max([FR_means.pos]+[FR_sterr.pos]);
                for sta=all_sta
                    for lin=1:n_lines
                        subplot_assignment(keys,'polar',n_states,n_lines,sta,lin,e,0,effectors_on_figure,1);
                        line_spec.color=keys.colors.(labels{lin})/255;
                        line_spec.linewidth=1;
                        line_spec.visible='on';
                        ph_polar_plot(FR_means(sta,lin).pos,FR_sterr(sta,lin).pos,unique_angles',line_spec,FR_mean_max)
                        tmp_factor=FR_mean_max/norm(FR_means(sta,lin).vector);
                        line([0,FR_means(sta,lin).vector(1)]*tmp_factor,[0,FR_means(sta,lin).vector(2)]*tmp_factor,'linewidth',2,'color',line_spec.color);
                    end
                end
                ph_title_and_save(FR_summary_handle,plot_title,[plot_title fig_title_part],keys);
            end
            
            %% Summary plot - average in summary???
            fig_title=sprintf('%s, %s, %s type, %s %s %s' ,population(u).unit_ID, population(u).target, type_string, keys.arrangement, title_part, title_value);
            plot_title        =[fig_title ' effector summary'];
            FR_summary_handle   =figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',[plot_title fig_title_part]);
            
            unique_columns=unique([T.column]);
            unique_rows=unique([T.row]);
            
            
            hs_ylim=0;
            hb_ylim=[0 0];
            hs=zeros(numel(unique_rows),numel(unique_columns));
            for r=1:numel(unique_rows)
                row=unique_rows(r);
                effectors_in_row=unique([T(fig_idx & [T.row]==row).effector]);
                for c=unique_columns
                    hs(r,c)=subplot_assignment(keys,'PSTH',n_states,n_lines,1,lin,r,c,unique_rows,unique_columns);
                    hold on
                    state_seperator=0;
                    for w=1:size(keys.PSTH_WINDOWS,1)
                        t_before_state=keys.PSTH_WINDOWS{w,3};
                        t_after_state=keys.PSTH_WINDOWS{w,4};
                        state_shift             =state_seperator-t_before_state;
                        bins                    =t_before_state:keys.PSTH_binwidth:t_after_state;
                        bins                    =bins+state_shift;
                        for h=1:numel(UC.hemifield)
                            for lin=1:size(CM,1)
                                col=keys.colors.([labels{lin} '_' side_labels_col{h}])/255;
                                tr_idx=fig_idx & te_idx & line_idx(lin,:) &[T.row]==row & [T.column]==c & [T.accepted]==1 & [T.hemifield]==UC.hemifield(h);
                                if sum(tr_idx)==0
                                    continue;
                                end
                                [PSTH,~,~,SEM]=ph_spike_density(pop.trial(tr_idx),T(tr_idx),w,keys,zeros(sum(tr_idx),1),ones(sum(tr_idx),1));
                                lineProps={'color',col,'linewidth',keys.UN.PSTH_summary_width};
                                shadedErrorBar(bins,PSTH,SEM,lineProps,1);
                                hs_ylim=max(hs_ylim,max(PSTH));
                            end
                        end
                        state_seperator=state_shift + t_after_state + 0.1;
                    end
                end
                
                %% Bar plots
                clear FR_summ FR_1quarters FR_2quarters FR_3quarters FR_labels
                
                for sta=all_sta
                    for er=1:numel(effectors_in_row)
                        effector=effectors_in_row(er);
                        e=find(effectors_on_figure==effector);
                        R=numel(effectors_on_figure)*(r-1)+e;
                        [type_effector_full, type_effector_short] = MPA_get_type_effector_name(type,effector);
                        con=0;
                        for h=1:numel(UC.hemifield)
                            for lin=1:size(CM,1)
                                con=con+1;
                                tr_idx=fig_idx & line_idx(lin,:) &[T.effector]==effector & [T.row]==row & [T.accepted]==1 & [T.hemifield]==UC.hemifield(h);
                                FR_labels{R,sta,con}=strrep([side_labels{h} ' ' legend_labels{lin}],'_',' ');
                                [FR_summ{R,sta,con},FR_1quarters(R,sta,con),FR_2quarters(R,sta,con),FR_3quarters(R,sta,con)]=deal(NaN);
                                if sum(tr_idx)==0
                                    continue;
                                end
                                %per_epoch=vertcat(o(tr_index).epoch);
                                FR_summ{R,sta,con}=pop.epochs_per_type{type}(tr_idx([T.type]==type),sta);
                                FR_1quarters(R,sta,con)=quantile(FR_summ{R,sta,con},0.25);
                                FR_2quarters(R,sta,con)=quantile(FR_summ{R,sta,con},0.5);
                                FR_3quarters(R,sta,con)=quantile(FR_summ{R,sta,con},0.75);
                                
                                hb_ylim=[min(hb_ylim(1),max([FR_1quarters(R,sta,con) nanmean(FR_summ{R,sta,con})])) ...
                                    max(hb_ylim(2),max([FR_3quarters(R,sta,con) nanmean(FR_summ{R,sta,con})]))];
                            end
                        end
                        
                        hb(sta+numel(all_sta)*(e-1))=subplot_assignment(keys,'Bars',n_states,n_lines,sta,lin,r,c,unique_rows,unique_columns,er,numel(effectors_in_row));
                        %% get title
                        state_label =keys.EPOCHS{sta,1};
                        anova_title_part1=get_anova_results(keys,current_unit_tuning,type_effector_short,keys.UN.anova_epoch1,['_' state_label '_']);
                        anova_title_part2=get_anova_results(keys,current_unit_tuning,type_effector_short,keys.UN.anova_epoch2,['_' state_label '_']);
                        title({state_label;anova_title_part1;anova_title_part2},'fontsize',8);
                        hold on
                        for b=1:size(FR_summ,3)
                            bm=nanmean(FR_summ{R,sta,b});
                            ci=[sterr(FR_summ{R,sta,b})*1.96 0];
                            ci=ci(1);
                            col=keys.colors.(labels_hem{b})/255;
                            bar(b,bm,1,'FaceColor',col);
                            eb=errorbar_constant_barsize_working(b, bm, ci, ci, 0.3);
                            patch([b-0.1 b+0.1 b+0.1 b-0.1],[FR_1quarters(R,sta,b) FR_1quarters(R,sta,b) FR_3quarters(R,sta,b) FR_3quarters(R,sta,b)],'k','EdgeColor','none');
                            scatter(b,FR_2quarters(R,sta,b),'o','markerfacecolor','w','markeredgecolor','k');
                            set(eb,'color','k');
                        end
                        set(gca,'xlim',[0.5 b+0.5],'Xtick',1:b,'Xticklabel',FR_labels(R,sta,:),'fontsize',8);
                        rotateXLabels(gca, 45);
                        
                        %% vertical titles per row
                        if sta==1
                            anova_title_part=get_anova_results(keys,current_unit_tuning,type_effector_short,keys.UN.anova_effector,'');
                            ytitle={type_effector_full; anova_title_part};
                            text(-2,0,ytitle,'rotation',90,'interpreter','none');
                            yticks=get(gca,'Ytick');
                            set(gca,'Yticklabel',repmat({''},size(yticks)),'fontsize',8);
                        end
                    end
                end
            end
            
            %% adding spikes and equalizing summary PSTH axes, adding background and text
            y_lim=[0, hs_ylim];
            spike_length=y_lim(2)*0.002/(1+sum(fig_idx & [T.accepted])*0.002);
            for r=1:numel(unique_rows)
                row=unique_rows(r);
                for c=unique_columns
                    subplot(hs(r,c))
                    hold on
                    
                    %% RASTER
                    state_seperator=0;
                    for w=1:size(keys.PSTH_WINDOWS,1)
                        sta=keys.PSTH_WINDOWS{w,2};
                        t_before_state=keys.PSTH_WINDOWS{w,3};
                        t_after_state=keys.PSTH_WINDOWS{w,4};
                        state_shift             =state_seperator-t_before_state;
                        raster_y=0;
                        for h=1:numel(UC.hemifield)
                            for lin=1:size(CM,1)
                                col=keys.colors.([labels{lin} '_' side_labels_col{h}])/255;
                                tr_idx=fig_idx & line_idx(lin,:) & te_idx & [T.row]==row & [T.column]==c & [T.accepted]==1 & [T.hemifield]==UC.hemifield(h);
                                n_trials     =sum(tr_idx);
                                line_struct=T(tr_idx);
                                if n_trials==0 || ~any([line_struct.states]==sta)
                                    raster_y=raster_y-n_trials*spike_length;
                                    continue;
                                end
                                
                                ATs=pop.trial(tr_idx);
                                line([state_shift+t_before_state state_shift+t_after_state],[raster_y raster_y],'Color',col,'LineWidth',0.5); hold on;
                                raster=[NaN; NaN];
                                for t=1:numel(line_struct)
                                    trial_state_onset=line_struct(t).states_onset(line_struct(t).states==sta);
                                    at_idx=ATs(t).arrival_times-trial_state_onset>=t_before_state & ATs(t).arrival_times-trial_state_onset<=t_after_state;
                                    AT=[ATs(t).arrival_times(at_idx)]'+state_shift-trial_state_onset;
                                    raster=[raster [AT;repmat(raster_y-t*spike_length,size(AT))]];
                                    %ig_make_raster(AT,raster_y-t*spike_length,spike_length,0,'Color',col,'LineWidth',keys.UN.raster_width);
                                end
                                scatter(raster(1,:),raster(2,:),10,col,'s','filled');
                                line([state_shift+t_before_state state_shift+t_after_state],[raster_y-t*spike_length raster_y-t*spike_length],'Color',col,'LineWidth',0.5);
                                if w==1
                                    text(state_shift+t_before_state,raster_y,['N=' num2str(t)],'Color',col,'fontsize',5,'verticalalignment','top');
                                end
                                raster_y=raster_y-n_trials*spike_length;
                            end
                        end
                        state_seperator=state_shift + t_after_state + 0.1;
                        y_lim(1)=min([y_lim(1), raster_y]);
                    end
                end
            end
            for r=1:numel(unique_rows)
                for c=unique_columns
                    subplot(hs(r,c))
                    ph_PSTH_background(T(fig_idx & te_idx),y_lim,y_lim,[0 y_lim(2)],keys,1)
                end
            end
            
            %% axes limits for bar plots
            if diff(hb_ylim)==0
                hb_ylim(2)=hb_ylim(1)+1;
            end
            if exist('hb','var')
                for hb_idx=find(hb~=0)
                    subplot(hb(hb_idx))
                    set(gca,'ylim',[hb_ylim]);
                end
            end
            ph_title_and_save(FR_summary_handle,plot_title,[plot_title fig_title_part],keys);
        end
    end
    
end

end

function plot_firing_rate_heatmap(R,n_columns,n_rows,subplot_pos)
X =NaN(n_columns,n_rows);
X(subplot_pos(~isnan(subplot_pos)))=R;
X=X';
map = gray;
map = flipud(map);
colormap(map);
X = [[X nan(size(X,1),1)] ; nan(1,size(X,2)+1)];
pcolor(X); set(gca,'Ydir','reverse');
axis equal; axis off;
end

function  anova_title_part=get_anova_results(keys,current_unit_tuning,type_effector_short,to_look_for,state_label)
anova_title_part='';
for ix=1:3:numel(to_look_for)
    label=[to_look_for{ix+1} state_label to_look_for{ix+2} '_' type_effector_short '_' keys.arrangement(1:3)];
    idx= DAG_find_column_index(current_unit_tuning,label);
    if ~isempty(idx)
        to_add=current_unit_tuning{2,idx};
    else
        to_add='-';
    end
    if ~isstr(to_add)
        to_add=num2str(to_add);
    end
    anova_title_part = [anova_title_part sprintf('%s:%s',to_look_for{ix},to_add)];
end
end

function ig_make_raster(t,inipos,len,dis,varargin)
% make_raster	- draws raster (short vertical line) plot for spike train
%--------------------------------------------------------------------------------
% Input(s): 	t	- spike arrival times
%			  t should be vector (one row) or 2D matrix (several rows)
%		inipos	- initial vertical position (optional, default 1)
%		len	- raster line length (optional, default 1)
%		dis	- distance between rows (optional, default 0)
%		varargin- Line properties (as in plot) (optional)
% Output(s):	h	- handler to raster
% Usage:	make_raster(t,inipos,len,dis,varargin);
%
% Last modified 23.11.02
% Copyright (c) 2002 Igor Kagan
% kigor@tx.technion.ac.il
% http://igoresha.virtualave.net
%--------------------------------------------------------------------------------
if nargin < 2,
    inipos = 1;
    len = 1;
    dis = 0;
elseif nargin < 3,
    len = 1;
    dis = 0;
elseif nargin < 4,
    dis = 0;
end

[rows, cols] = size(t);
if min(size(t))==1, t = t(:)'; end


for r = 1:rows,
    
    x = zeros(3*cols,1);
    x(1:3:3*cols,1) = t(r,:)';
    x(2:3:3*cols,1) = t(r,:)';
    x(3:3:3*cols,1) = NaN;
    
    y = zeros(3*cols,1);
    y1 = inipos + (r-1)*ones(cols,1) + (r-1)*dis;
    y2 = inipos + (r-1+len)*ones(cols,1);
    
    y(1:3:3*cols,:) = y1;
    y(2:3:3*cols,:) = y2;
    y(3:3:3*cols,:) = NaN;
    
    line(x,y,varargin{:});
    %h = [h hh];
end

end

function handle=subplot_assignment(keys,subplot_case,n_states,n_lines,sta,lin,row,column,rows,columns,effector,effectors_per_row)
if keys.plot.polars_on_extra_figure
    n_rows_FR           =n_lines;
    n_rows_PSTH         =2*numel(rows);
    n_columns_PSTH      =numel(columns);
    FR_row_shift        =n_states*(lin-1);
    polar_row_shift     =0;
else
    n_rows_FR           =2*n_lines;
    n_rows_PSTH         =2*numel(rows);
    n_columns_PSTH      =numel(columns);
    FR_row_shift        =2*n_states*(lin-1);
    polar_row_shift     =n_states;
end
switch subplot_case
    case 'PSTH'
        handle=subplot(n_rows_PSTH,n_columns_PSTH,2*n_columns_PSTH*(row-1)+column);
    case 'Bars'
        handle=subplot(n_rows_PSTH*effectors_per_row,n_states,sta+n_states*(effectors_per_row*(2*row-1) + effector-1));
        pos=get(handle,'Position'); rescale=80/100; posnew=[pos(1)+pos(3)*(1-rescale) pos(2)+pos(4)*(1-rescale) pos(3)*rescale pos(4)*rescale];
        set(handle,'Position',posnew);
    case 'FR'
        handle=subplot(n_rows_FR,n_states,sta+FR_row_shift);
    case 'polar'
        handle=subplot(n_rows_FR,n_states,sta+FR_row_shift+polar_row_shift);
end
end

function ig_set_caxis_equal_lim(h)
% ig_set_axis_equal_lim	- set CAXIS equal (min max of all)
%--------------------------------------------------------------------------------
% Input(s): 	h - set of handlers
%               MODE - 'Xlim' or 'Ylim'
% Output(s):	none
% Usage:	set_axes_equal_lim(h,'Ylim')

%
% Last modified 11.12.02
% Copyright (c) 2002 Igor Kagan
% kigor@tx.technion.ac.il
% http://igoresha.virtualave.net
%--------------------------------------------------------------------------------

if nargin < 1,
    h = ig_get_figure_axes;
end

if length(h)<2,
    return;
end

Lim = get(h,'Clim');

a=0;
isreallyaxes=false(size(h));
for i=1:size(Lim,1);
    if isnumeric(Lim{i})
        a=a+1;
        isreallyaxes(i)=true;
        L(a,:) = Lim{i};
    end
end
set(h(isreallyaxes),'Clim',[min(min(L)) max(max(L))]);
end



