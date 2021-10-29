function ph_plot_unit_per_condition2(population,keys)
tuning_per_unit_table=keys.tuning_per_unit_table;

keys.path_to_save=[keys.basepath_to_save keys.project_version filesep 'single_cell_examples' filesep];
%% colors and ylim preparation
eye_offset=-keys.plot.trials_max_for_ylim-keys.plot.excentricity_max_for_ylim*keys.plot.eyetrace_factor;
hnd_offset=eye_offset-keys.plot.excentricity_max_for_ylim*keys.plot.eyetrace_factor*2;
eye_col_v=keys.colors.eye_ver;
eye_col_h=keys.colors.eye_hor;

for unit=1:numel(population)
    types       =[population(unit).trial.type];
    effectors   =[population(unit).trial.effector];
    hands       =[population(unit).trial.reach_hand];
    
    % copy of unit identifiers to keys... redundant?
    [keys.unit_ID keys.stability_rating keys.SNR_rating keys.Single_rating keys.target keys.grid_x keys.grid_y keys.electrode_depth keys.channel keys.block_unit] =...
        deal(population(unit).unit_ID, population(unit).stability_rating, population(unit).SNR_rating, population(unit).Single_rating, population(unit).target,...
        population(unit).grid_x, population(unit).grid_y, population(unit).electrode_depth, population(unit).channel, population(unit).block_unit);
    
    current_unit_tuning= [tuning_per_unit_table(1,:); tuning_per_unit_table(ismember(tuning_per_unit_table(:,1), population(unit).unit_ID),:)];
    for a=1:numel(keys.position_and_plotting_arrangements)
        keys.arrangement=keys.position_and_plotting_arrangements{a};
        for type=unique(types)
            
            %% selection of epochs to plot
            effectors_effector_loop=unique(effectors);
            [keys all_sta]=ph_get_epoch_keys(keys,type,effectors_effector_loop,1);
            n_states=numel(all_sta);
            
            o_index= types == type & ismember(effectors,effectors_effector_loop);% & ismember(hands,keys.reach_hand);
            if keys.cal.completed
                o_index = o_index & [population(unit).trial.completed]==1;  %% success/completed
            end
            eff_for_typ=effectors(o_index);
            effectors_effector_loop=unique(eff_for_typ);
            
            if sum(o_index)==0;
                disp(sprintf('%s has no trials for effector %.0f type %.0f hands %s completed= %.0f',population(unit).unit_ID,effector,type,mat2str(keys.cal.reach_hand),keys.cal.completed));
                continue;
            end
            o=ph_arrange_positions_and_plots(keys,population(unit).trial(o_index),population(unit));
            T=o.trial;
            if any(hands(o_index)>0)
                y_lim_PSTH=hnd_offset-keys.plot.excentricity_max_for_ylim*keys.plot.eyetrace_factor;
            else
                y_lim_PSTH=eye_offset-keys.plot.excentricity_max_for_ylim*keys.plot.eyetrace_factor;
            end
            unique_figures=unique([T.figure]);
            for fig=unique_figures
                title_value=vertcat(o.figure_title_value{fig});
                title_part=o.figure_title_part;
                fig_idx=[T.figure]==fig;
                
                %% new part!
                [UC, CM, labels]=ph_get_condition_matrix(T(fig_idx),keys);
                %% for general use
                for L=1:size(CM,1)
                    line_idx(L,:)=true(size(fig_idx));
                    for con=1:size(CM,2)
                        line_idx(L,:)=line_idx(L,:) & [T.(keys.condition_parameters{con})]==CM(L,con);
                    end
                end
                n_lines=size(CM,1);
                %unique_lines=1:size(CM,1); %% only needed now to specify average plotting
                legend_labels_hem={};
                for h=UC.hemifield % append hemifield labels, careful with the order!
                    legend_labels_hem=[legend_labels_hem; strcat(labels,['_' keys.labels.hemifield{h+2}])];
                end
                legend_labels_hem=reshape(legend_labels_hem,numel(legend_labels_hem),1);
                side_labels={'L','R'};
                side_labels_col={'IS','CS'};
                %                 if keys.plot.average_heat_maps
                %                     unique_lines(end)=[];
                %                 end
                if any(UC.hemifield==0) %&& fig==unique_figures(1)% vertical targets, KK uncomment, why figure and whz colors?
                    side_labels={'L','V','R'};
                    side_labels_col={'IS','VS','CS'};
                end
                %% IDEA!: take hemifields as conditions here already as well!!!
                % if ~keys.plot.average_PSTH_line
                %                                         % this just makes it simpler, looping several times but plotting the same all over
                %                                         tr_index= tr_index &  [T.line]==lin;
                %                                     end
                
                effectors_on_figure=effectors_effector_loop(ismember(effectors_effector_loop,eff_for_typ(fig_idx)));
                for e=1:numel(effectors_on_figure)
                    clear hh hs hb
                    effector=effectors_on_figure(e);
                    [type_effector_full, type_effector_short type_string] = MPA_get_type_effector_name(type,effector);
                    
                    %% ANOVA results to print
                    anova_title_part=get_anova_results(keys,current_unit_tuning,type_effector_short,keys.plot.anova_main,'');
                    
                    fig_title=sprintf('%s, %s, %s, %s %s %s',keys.unit_ID, keys.target, type_effector_full, keys.arrangement, title_part, title_value);
                    fig_title_part=sprintf(', Stability %d, SNR %d, Single %d Grid: %d/%d Depth %.2f ANOVA %s Channel: %d Blocks&Units: %s ', ...
                        keys.stability_rating, keys.SNR_rating, keys.Single_rating, keys.grid_x, keys.grid_y, keys.electrode_depth, anova_title_part, keys.channel, [keys.block_unit{:}]);
                    
                    %% Per position PSTH figure
                    plot_1_title            = [fig_title  ' PSTHs'];
                    PSTH_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',[plot_1_title fig_title_part]);
                    subplot_indizes=unique([T(fig_idx & [T.effector]==effector).subplot_pos]);
                    y_max =1;
                    yL_min=0;
                    for sub=subplot_indizes
                        subplot(o.rows,o.columns,sub,'XlimMode','Manual','YlimMode','Manual','XtickMode','Manual')
                        hold on
                        subplot_struct=T([T.subplot_pos]==sub & [T.effector]==effector);
                        title([subplot_struct(1).title_part num2str(round(subplot_struct(1).position))]);
                        state_seperator=zeros(n_lines,1);
                        state_seperator_max=0.2;
                        for w=1:size(keys.PSTH_WINDOWS,1)
                            sta=keys.PSTH_WINDOWS{w,2};
                            t_before_state=keys.PSTH_WINDOWS{w,3};
                            t_after_state=keys.PSTH_WINDOWS{w,4};
                            histo=[];
                            lines_for_average=[];
                            raster_y=0;
                            for lin=1:size(CM,1)
                                col=keys.colors.(labels{lin})/255;
                                line_struct=T(line_idx(lin,:) & [T.effector]==effector & [T.figure]==fig);
                                n_trials=numel(line_struct);
                                state_shift=state_seperator(lin)-t_before_state;
                                if n_trials==0 || ~any([line_struct.states]==sta)
                                    state_seperator(lin)=state_shift + t_after_state + 0.1;
                                    continue;
                                else
                                    lines_for_average =[lines_for_average; lin];
                                end
                                
                                %% Raster and eye hand traces  + Number of trials line and text
                                line([state_shift+t_before_state state_shift+t_after_state],[raster_y raster_y],'Color',col,'LineWidth',0.5); hold on;
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
                                    time_axis=line_struct(t).time_axis-trial_state_onset+state_shift;
                                    t_idx=line_struct(t).time_axis-trial_state_onset>=t_before_state &...
                                        line_struct(t).time_axis-trial_state_onset<=t_after_state;
                                    at_idx=line_struct(t).arrival_times-trial_state_onset>=t_before_state &...
                                        line_struct(t).arrival_times-trial_state_onset<=t_after_state;
                                    if keys.plot.eye_hand_traces
                                        line(time_axis(t_idx),line_struct(t).x_eye(t_idx)*keys.plot.eyetrace_factor + eye_offset ,'color',eye_col_h);
                                        line(time_axis(t_idx),line_struct(t).y_eye(t_idx)*keys.plot.eyetrace_factor + eye_offset,'color',eye_col_v);
                                        line(time_axis(t_idx),line_struct(t).x_hnd(t_idx)*keys.plot.hndtrace_factor + hnd_offset,'color',hnd_col_h);
                                        line(time_axis(t_idx),line_struct(t).y_hnd(t_idx)*keys.plot.hndtrace_factor + hnd_offset,'color',hnd_col_v);
                                    end
                                    ig_make_raster([line_struct(t).arrival_times(at_idx)]'+state_shift-trial_state_onset,raster_y-t,1,0,'Color',col,'LineWidth',keys.width.raster);
                                end
                                line([state_shift+t_before_state state_shift+t_after_state],[raster_y-t raster_y-t],'Color',col,'LineWidth',0.5);
                                
                                if w==1
                                    not_accepted=find(~[line_struct.accepted])*-1+raster_y;
                                    plot((state_shift+t_before_state)*ones(size(not_accepted)),not_accepted,'marker','s','markeredgecolor','k','markerfacecolor','k','linestyle','none');
                                    text(state_shift+t_before_state,raster_y,['N=' num2str(t)],'Color',col,'fontsize',5,'verticalalignment','top');
                                end
                                raster_y=raster_y-n_trials;
                                
                                % PSTH
                                state_seperator(lin)=state_shift + t_after_state + 0.1;
                                [histo(lin,:), bins, ~, SEM]=ph_spike_density(line_struct,w,keys,zeros(numel(line_struct),1),ones(numel(line_struct),1));
                                bins=bins+state_shift;
                                
                                lineProps={'color',col,'linewidth',keys.width.PSTH_perpos};
                                shadedErrorBar(bins,histo(lin,:),SEM,lineProps,1);
                            end
                            
                            if keys.plot.average_PSTH_line
                                line(bins,nanmean(histo(lines_for_average,:),1),'color','k','LineWidth',keys.width.PSTH_perpos);          %PSTH
                            end
                            state_seperator_max=max([state_seperator_max;state_seperator]);
                            y_max=max([y_max max(histo)]);
                            yL_min=min([yL_min raster_y]);
                        end
                    end
                    y_lim=[y_lim_PSTH,y_max];
                    for sub=subplot_indizes
                        subplot(o.rows,o.columns,sub,'XlimMode','Manual','YlimMode','Manual','XtickMode','Manual')
                        ph_PSTH_background(T([T.effector]==effector),y_lim,[yL_min y_lim(2)],[0 y_lim(2)],keys,0.5)
                        set(gca,'ylim',y_lim,'xlim',[0 state_seperator_max-0.1]);
                    end
                    ph_title_and_save(PSTH_summary_handle,plot_1_title,[plot_1_title fig_title_part],keys);
                    
                    %% FR figure
                    plot_title        =[fig_title ' FR'];
                    FR_summary_handle   =figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',[plot_title fig_title_part]);
                    tr_index=fig_idx & [T.effector]==effector & [T.accepted]==1;
                    if sum(tr_index)==0 % might be the case if all trials are excluded cause of no spikes
                        continue;
                    end
                    % Computing mean firing rates per position
                    %                     if keys.plot.average_heat_maps
                    %                         o.line_labels=[o.line_labels {'Av'}];
                    %                         unique_lines=[unique_lines max(unique_lines)+1];
                    %                         o.PSTH_perpos_colors=[o.PSTH_perpos_colors; 0 0 0];
                    %                     end
                    subplot_indizes=1:max([T.subplot_pos]);
                    clear FR_heat;
                    
                    for sub=subplot_indizes
                        for sta=all_sta
                            for lin=1:size(CM,1)
                                lin_idx = line_idx(lin,:) & [T.subplot_pos]==sub & [T.figure]==fig;
                                if sum(lin_idx)==0
                                    FR_heat(sta,lin).pos(sub)    =NaN;
                                else
                                    line_struct=T(lin_idx);
                                    states_trials       =vertcat(line_struct.epoch);
                                    FR_heat(sta,lin).pos(sub)    =nanmean([states_trials(:,sta).FR]);
                                end
                            end
                            lin_idx=tr_index & [T.subplot_pos]==sub & [T.figure]==fig;
                            if keys.plot.average_PSTH_line && sum(lin_idx)>0
                                states_trials=   vertcat(T(lin_idx).epoch);
                                FR_heat(sta,lin+1).pos(sub)=nanmean([states_trials(:,sta).FR]);
                            elseif keys.plot.average_PSTH_line
                                FR_heat(sta,lin+1).pos(sub)=NaN;
                            end
                        end
                    end
                    
                    % Plotting heat maps
                    for sta=all_sta
                        for lin=1:size(CM,1)
                            state_label =keys.EPOCHS{sta,1};
                            hh(sta+n_states*(lin-1))=subplot_assignment(keys,'FR',n_states,n_lines,sta,lin,e,0,effectors_on_figure,1);
                            plot_firing_rate_heatmap(FR_heat(sta,lin).pos,o.columns,o.rows,1:max(subplot_indizes));
                            title([strrep(labels{lin},'_',' ') ' ' state_label]);
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
                    polar_struct=T(tr_index);
                    positions=vertcat(polar_struct.position);
                    position_angles=angle(positions(:,1)+1i*positions(:,2));
                    position_normfactor=abs(positions(:,1)+1i*positions(:,2));
                    position_normfactor(position_normfactor==0)=1;
                    norm_vector=positions./[position_normfactor position_normfactor];
                    [unique_angles]=unique(position_angles);
                    
                    states_trials_sp       =vertcat(polar_struct.epoch);
                    for sta=all_sta
                        sta_FR=[states_trials_sp(:,sta).FR]';
                        norm_factor=max(sta_FR);
                        sta_vector=[sta_FR.*norm_vector(:,1) sta_FR.*norm_vector(:,2)]/norm_factor;
                        for lin=1:n_lines
                            lin_idx=line_idx(lin,tr_index);
                            
                            FR_means(sta,lin).vector=nanmean(sta_vector(lin_idx,:),1);
                            for ang=1:numel(unique_angles)
                                if keys.plot.average_PSTH_line && lin==n_lines
                                    line_struct=polar_struct(position_angles==unique_angles(ang));
                                else
                                    line_struct=polar_struct(position_angles==unique_angles(ang) & lin_idx');
                                end
                                states_trials       =vertcat(line_struct.epoch);
                                if isempty(states_trials)% || isempty(subplot_struct.position) %unfortunate debugging for empty positions !
                                    FR_means(sta,lin).pos(ang)    =NaN;
                                    FR_sterr(sta,lin).pos(ang)    =NaN;
                                else
                                    FR_means(sta,lin).pos(ang)    =nanmean([states_trials(:,sta).FR])/norm_factor;
                                    FR_sterr(sta,lin).pos(ang)    =sterr([states_trials(:,sta).FR])/norm_factor;
                                end
                            end
                        end
                    end
                    
                    FR_mean_max=max([FR_means.pos]+[FR_sterr.pos]);
                    for sta=all_sta
                        for lin=1:size(CM,1)
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
                
                %% Summary plot
                fig_title=sprintf('%s, %s, %s type, %s %s %s' ,keys.unit_ID, keys.target, type_string, keys.arrangement, title_part, title_value);
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
                    %% PSTH summary R/L and U/D subplot
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
                                    tr_index=fig_idx & line_idx(lin,:) &[T.row]==row & [T.column]==c & [T.accepted]==1 & [T.hemifield]==UC.hemifield(h);
                                    %                                     if ~keys.plot.average_PSTH_line
                                    %                                         % this just makes it simpler, looping several times but plotting the same all over
                                    %                                         tr_index= tr_index &  [T.line]==lin;
                                    %                                     end
                                    if sum(tr_index)==0
                                        continue;
                                    end
                                    [PSTH,~,~,SEM]=ph_spike_density(T(tr_index),w,keys,zeros(sum(tr_index),1),ones(sum(tr_index),1));
                                    lineProps={'color',col,'linewidth',keys.width.PSTH_summary};
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
                                    tr_index=fig_idx & line_idx(lin,:) &[T.effector]==effector & [T.row]==row & [T.accepted]==1 & [T.hemifield]==UC.hemifield(h);
                                    FR_labels{R,sta,con}=''; [FR_summ{R,sta,con},FR_1quarters(R,sta,con),FR_2quarters(R,sta,con),FR_3quarters(R,sta,con)]=deal(NaN);
                                    if sum(tr_index)==0
                                        continue;
                                    end
                                    per_epoch=vertcat(T(tr_index).epoch);
                                    FR_summ{R,sta,con}=[per_epoch(:,sta).FR];
                                    FR_1quarters(R,sta,con)=quantile(FR_summ{R,sta,con},0.25);
                                    FR_2quarters(R,sta,con)=quantile(FR_summ{R,sta,con},0.5);
                                    FR_3quarters(R,sta,con)=quantile(FR_summ{R,sta,con},0.75);
                                    FR_labels{R,sta,con}=strrep([side_labels{h} ' ' labels{lin}],'_',' ');
                                    
                                    hb_ylim=[min(hb_ylim(1),max([FR_1quarters(R,sta,con) nanmean(FR_summ{R,sta,con})])) ...
                                        max(hb_ylim(2),max([FR_3quarters(R,sta,con) nanmean(FR_summ{R,sta,con})]))];
                                end
                            end
                            
                            hb(sta+numel(all_sta)*(e-1))=subplot_assignment(keys,'Bars',n_states,n_lines,sta,lin,r,c,unique_rows,unique_columns,er,numel(effectors_in_row));
                            %% get title
                            state_label =keys.EPOCHS{sta,1};
                            anova_title_part1=get_anova_results(keys,current_unit_tuning,type_effector_short,keys.plot.anova_epoch1,['_' state_label '_']);
                            anova_title_part2=get_anova_results(keys,current_unit_tuning,type_effector_short,keys.plot.anova_epoch2,['_' state_label '_']);
                            title({state_label;anova_title_part1;anova_title_part2},'fontsize',8);
                            hold on
                            for b=1:size(FR_summ,3)
                                bm=nanmean(FR_summ{R,sta,b});
                                ci=[sterr(FR_summ{R,sta,b})*1.96 0];
                                ci=ci(1);
                                col=keys.colors.(legend_labels_hem{b})/255;
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
                                anova_title_part=get_anova_results(keys,current_unit_tuning,type_effector_short,keys.plot.anova_effector,'');
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
                                for lin=1:n_lines
                                    col=keys.colors.([labels{lin} '_' side_labels_col{h}])/255;
                                    tr_index=fig_idx & line_idx(lin,:) & [T.row]==row & [T.column]==c & [T.accepted]==1 & [T.hemifield]==UC.hemifield(h);
                                    %                                     if ~keys.plot.average_PSTH_line
                                    %                                         % this just makes it simpler, looping several times but plotting the same all over
                                    %                                         tr_index= tr_index &  [T.line]==lin;
                                    %                                     end
                                    n_trials     =sum(tr_index);
                                    line_struct=T(tr_index);
                                    spike_length=y_lim(2)*0.002/(1+sum(fig_idx & [T.accepted])*0.002);
                                    if n_trials==0 || ~any([line_struct.states]==sta)
                                        raster_y=raster_y-n_trials*spike_length;
                                        continue;
                                    end
                                    
                                    %% Raster
                                    line([state_shift+t_before_state state_shift+t_after_state],[raster_y raster_y],'Color',col,'LineWidth',0.5); hold on;
                                    for t=1:numel(line_struct)
                                        trial_state_onset=line_struct(t).states_onset(line_struct(t).states==sta);
                                        at_idx=line_struct(t).arrival_times-trial_state_onset>=t_before_state &...
                                            line_struct(t).arrival_times-trial_state_onset<=t_after_state;
                                        ig_make_raster([line_struct(t).arrival_times(at_idx)]'+state_shift-trial_state_onset,raster_y-t*spike_length,spike_length,0,'Color',col,'LineWidth',keys.width.raster);
                                    end
                                    
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
                        ph_PSTH_background(T,y_lim,y_lim,[0 y_lim(2)],keys,1)
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

% function title_and_save(figure_handle,filename,plot_title,keys)
% mtit(figure_handle,  plot_title, 'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 8,'Interpreter', 'none');
% stampit;
% if keys.plot.export
%     wanted_size=[50 30];
%     set(figure_handle, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
%     export_fig([keys.tuning_table_foldername filesep 'single_cell_examples' filesep filename], '-pdf','-transparent') % pdf by run
%     close all
% end
% end

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


                
%                 %% per movement plot
%                 keys.movement_angle_binwidth=60;
%                 keys.movement_amplitude_binwidth=10;
%                 keys.movement_plot_type='revealed_or_not';
%                 keys.movement_plot_types={'pre','peri','post'};
%                 
%                 binwidth=keys.movement_angle_binwidth*pi/180;
%                 angular_movement_bins=-pi:binwidth:pi;
%                 
%                 for e=1:numel(effectors_on_figure)
%                     effector=effectors_on_figure(e);
%                     o_e=o;
%                     o_e.trial=o_e.trial(of_index & [T.effector]==effector);
%                     [type_effector_full, type_effector_short] = MPA_get_type_effector_name(type,effector);
%                     om=ph_arrange_movement_positions(population(unit).trial(oe_index),keys); %why not using o_e => this would separate into error and success
%                     PSTH_perpos_colors=om.PSTH_perpos_colors;
%                     
%                     plot_4_title            = [fig_title  ' movement PSTHs '];
%                     PSTH_movement_handle     = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',[plot_4_title]);
%                     fig_title_part=[type_effector_short];
%                     clear sp
%                     for sub= unique([om.movement(~isnan([om.movement.subplot_pos])).subplot_pos])
%                         sp(sub)=subplot(om.rows,om.columns,sub);
%                         os=om.movement([om.movement.subplot_pos]==sub);
%                         title([om.sub_title num2str(om.subplot_titles(sub))]);
%                         hold on
%                         raster_y=20;
%                         unique_lines=unique([os.line]);
%                         unique_lines(isnan(unique_lines))=[];
%                         for l=1:numel(unique_lines)
%                             lin=unique_lines(l);
%                             line_idx=[os.line]==lin & arrayfun(@(x) sum(~isnan(x.spike_density))>0,os); %% unfortunate
%                             to_plot=vertcat(os(line_idx).spike_density);
%                             raster_to_plot={os(line_idx).arrival_times};
%                             t_b=-0.2;
%                             t_a= 0.2;
%                             PSTH_x_axis=t_b:keys.PSTH_binwidth:t_a;
%                             for t=1:numel(raster_to_plot)
%                                 ig_make_raster(raster_to_plot{t}',raster_y+t*0.1,1,0,'Color',PSTH_perpos_colors(l,:),'LineWidth',1)
%                             end
%                             raster_y=raster_y+t*0.1;
%                             line([t_b t_a],[raster_y+1 raster_y+1],'Color',PSTH_perpos_colors(l,:),'LineWidth',0.5);
%                             text(t_b,raster_y,['N=' num2str(t)],'Color',PSTH_perpos_colors(l,:),'fontsize',5);
%                             line(PSTH_x_axis,nanmean(to_plot,1),'color',PSTH_perpos_colors(l,:));
%                         end
%                         y_lim_m(sub,:)=get(gca,'ylim');
%                     end
%                     legend(om.line_labels)
%                     y_lim_m_max=max(max(y_lim_m));
%                     y_lim_m_min=min(min(y_lim_m));
%                     for sub= unique([om.movement(~isnan([om.movement.subplot_pos])).subplot_pos])
%                         subplot(sp(sub));
%                         rectangle('Position',[-0.15 y_lim_m_min 0.14 0.05*diff([y_lim_m_min y_lim_m_max])],'EdgeColor','none','FaceColor',[0.9 0.9 0.9]) %frame for indicating state separation
%                         rectangle('Position',[-0.01 y_lim_m_min 0.06 0.05*diff([y_lim_m_min y_lim_m_max])],'EdgeColor','none','FaceColor',[0.6 0.6 0.6]) %frame for indicating state separation
%                         rectangle('Position',[ 0.05 y_lim_m_min 0.1  0.05*diff([y_lim_m_min y_lim_m_max])],'EdgeColor','none','FaceColor',[0.3 0.3 0.3]) %frame for indicating state separation
%                         
%                         line([0 0],[y_lim_m_min y_lim_m_max],'color','k');
%                         set(gca,'ylim',[y_lim_m_min y_lim_m_max]);
%                     end
%                     ph_title_and_save(PSTH_movement_handle,plot_4_title,[plot_4_title fig_title_part],keys);
%                     
%                     
%                     plot_5_title            = [fig_title  ' movement polar '];
%                     fig_title_part=[type_effector_short];
%                     polar_movement_handle     = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',[plot_5_title]);
%                     keys.movement_plot_types={'pre','peri','post'};
%                     row_fields={'no_target_revealed','target_revealed'};
%                     n_columns=numel(keys.movement_plot_types);
%                     % get unique shapes
%                     oo=[population(unit).trial(oe_index)];
%                     oo=[oo.movement];
%                     unique_shapes=unique([[oo.shape_at_start] [oo.shape_at_end]]);
%                     unique_shapes(isnan(unique_shapes))=[];
%                     %clear PSTH_perpos_colors
%                     % calculate mean and std of FR per angle
%                     FR_rev_mean=NaN(numel(keys.movement_plot_types),1);
%                     FR_rev_sem=NaN(numel(keys.movement_plot_types),1);
%                     for n_column=1:numel(keys.movement_plot_types)
%                         epoch=keys.movement_plot_types{n_column};
%                         keys.movement_plot_type=epoch;
%                         om=ph_arrange_movement_positions(population(unit).trial(oe_index),keys);
%                         %why not using o_e => this would separate into error and success
%                         PSTH_colors_cell{n_column}=om.PSTH_perpos_colors;
%                         for n_row=1:numel(row_fields)
%                             row_field=row_fields{n_row};
%                             subplot(numel(row_fields)+1,n_columns,n_column+(n_row-1)*n_columns);
%                             or=om.movement([om.movement.(row_field)]==1);
%                             unique_lines=unique([or.line]);
%                             %for lin=unique_lines
%                             for pos=1:numel(angular_movement_bins)
%                                 pos_m=pos;
%                                 if pos==numel(angular_movement_bins)
%                                     pos_m=1;
%                                 end
%                                 %m_index=[or.subplot_pos]==pos_m;% & [or.line]==lin;
%                                 m_index=[or.subplot_value]==angular_movement_bins(pos_m);% & [or.line]==lin;
%                                 FR_mean(l,pos,n_column,n_row)=nanmean([or(m_index).FR]);
%                                 FR_sem(l,pos,n_column,n_row)=sterr([or(m_index).FR]);
%                             end
%                             %end
%                         end
%                         for s=1:numel(unique_shapes)
%                             shape=unique_shapes(s);
%                             s_index=[or.shape_revealed]==shape;
%                             FR_rev_mean(n_column,s)=nanmean([or(s_index).FR]);
%                             FR_rev_sem(n_column,s)=sterr([or(s_index).FR]);
%                         end
%                     end
%                     FR_mean_max=max(FR_mean(:));
%                     FR_revealed_max=max(FR_rev_mean(:)+FR_rev_sem(:));
%                     for n_column=1:numel(keys.movement_plot_types)
%                         epoch=keys.movement_plot_types{n_column};
%                         for n_row=1:numel(row_fields)
%                             row_field=row_fields{n_row};
%                             subplot(numel(row_fields)+1,n_columns,n_column+(n_row-1)*n_columns);
%                             line_spec.visible='off';
%                             polar(0,FR_mean_max,line_spec);%,'color',PSTH_perpos_colors(lin,:))
%                             hold on
%                             %unique_lines=unique([or.line]);
%                             %for lin=unique_lines
%                             line_spec.color=PSTH_colors_cell{n_column}(l,:);
%                             line_spec.linewidth=2;
%                             line_spec.visible='on';
%                             ph_polar_plot(FR_mean(l,:,n_column,n_row),FR_sem(l,:,n_column,n_row),angular_movement_bins,line_spec)
%                             %polar(angular_movement_bins,FR_mean(lin,:,n_column,n_row),line_spec);%,'color',PSTH_perpos_colors(lin,:))
%                             hold on
%                             %end
%                             if n_row==1
%                                 title(epoch,'interpreter','none') %om.sub_title ' '
%                             end
%                             if n_column==1
%                                 ylabel(row_field,'interpreter','none')
%                             end
%                         end
%                         
%                         subplot(numel(row_fields)+1,n_columns,n_column+numel(row_fields)*n_columns);
%                         hold on
%                         bar(FR_rev_mean(n_column,:),'FaceColor','g');
%                         erb=errorbar_constant_barsize_working([1:size(FR_rev_mean,2)]', FR_rev_mean(n_column,:), FR_rev_sem(n_column,:), FR_rev_sem(n_column,:), 0.3, 'o');
%                         set(erb,'color','k');
%                         set(gca,'xticklabel',unique_shapes,'xlim',[0.5 size(FR_rev_mean,2)+0.5],'ylim',[0 max([FR_revealed_max 1])]);
%                         xlabel('shape','interpreter','none')
%                         if n_column==1
%                             ylabel('FR');
%                         end
%                     end
%                     ph_title_and_save(polar_movement_handle,plot_5_title,[plot_5_title fig_title_part],keys);
%                 end


