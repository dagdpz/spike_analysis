function ph_plot_unit_per_condition(population,keys)
tuning_per_unit_table=keys.tuning_per_unit_table;

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
            
            if any(hands(o_index)>0)
                y_lim_PSTH=hnd_offset-keys.plot.excentricity_max_for_ylim*keys.plot.eyetrace_factor;
            else
                y_lim_PSTH=eye_offset-keys.plot.excentricity_max_for_ylim*keys.plot.eyetrace_factor;
            end
            unique_figures=unique([o.trial.figure]);
            for fig=unique_figures
                title_value=vertcat(o.figure_title_value{fig});
                title_part=o.figure_title_part;
                of_index=[o.trial.figure]==fig;
                effectors_on_figure=effectors_effector_loop(ismember(effectors_effector_loop,eff_for_typ(of_index)));
                unique_lines=unique([o.trial(of_index).line]);
                for e=1:numel(effectors_on_figure)
                    clear hh hs hb
                    effector=effectors_on_figure(e);
                    [type_effector_full, type_effector_short type_string] = MPA_get_type_effector_name(type,effector);
                    
                    %% ANOVA results to print
                    anova_labels.epoch  = ['in_epoch_main_' type_effector_short '_' keys.arrangement(1:3)];
                    anova_labels.space  = ['in_spaceLR_main_'  type_effector_short '_' keys.arrangement(1:3)];
                    anova_labels.choice = ['ch_spaceLR_main_'  type_effector_short '_' keys.arrangement(1:3)];
                    anova_labels.hands  = ['in_hands_main_' type_effector_short '_' keys.arrangement(1:3)];
                    anova_labels.ExS    = ['in_ExS_' type_effector_short '_' keys.arrangement(1:3)];
                    anova_labels.ExH    = ['in_ExH_' type_effector_short '_' keys.arrangement(1:3)];
                    anova_labels.SxH    = ['in_SxH_' type_effector_short '_' keys.arrangement(1:3)];
                    anova_labels.SxH    = ['in_SxH_' type_effector_short '_' keys.arrangement(1:3)];
                    for FN=fieldnames(anova_labels)'
                        idx= DAG_find_column_index(current_unit_tuning,anova_labels.(FN{:}));
                        if ~isempty(idx)
                            anova_results(e).(FN{:})=current_unit_tuning{2,idx};
                        else
                            anova_results(e).(FN{:})='-';
                        end
                        if ~isstr(anova_results(e).(FN{:}))
                            anova_results(e).(FN{:})=num2str(anova_results(e).(FN{:}));
                        end
                    end
                    
                    fig_title=sprintf('%s, %s, %s, %s %s %s' ,...
                        keys.unit_ID, keys.target, type_effector_full, keys.arrangement, title_part, title_value);
                    fig_title_part=sprintf(', Stability %d, SNR %d, Single %d Grid: %d/%d Depth %.2f ANOVA E:%s S:%s C:%s H:%s ExS:%s ExH:%s SxH:%s Channel: %d Blocks&Units: %s ', ...
                        keys.stability_rating, keys.SNR_rating, keys.Single_rating, keys.grid_x, keys.grid_y, keys.electrode_depth, ...
                        anova_results(e).epoch, anova_results(e).space, anova_results(e).choice, anova_results(e).hands,...
                        anova_results(e).ExS, anova_results(e).ExH, anova_results(e).SxH, keys.channel, [keys.block_unit{:}]);
                    
                    %% Per position PSTH figure
                    plot_1_title            = [fig_title  ' PSTHs'];
                    %PSTH_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',[plot_1_title fig_title_part]);
                    PSTH_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',[plot_1_title fig_title_part]);
                    subplot_indizes=unique([o.trial(of_index & [o.trial.effector]==effector).subplot_pos]);
                    PSTH_perpos_colors=o.PSTH_perpos_colors;
                    y_max =1;
                    yL_min=0;
                    for sub=subplot_indizes
                        subplot(o.rows,o.columns,sub,'XlimMode','Manual','YlimMode','Manual','XtickMode','Manual')
                        hold on
                        subplot_struct=o.trial([o.trial.subplot_pos]==sub & [o.trial.effector]==effector);
                        title([subplot_struct(1).title_part num2str(round(subplot_struct(1).position))]);
                        state_seperator=zeros(max(unique_lines),1);
                        state_seperator_max=0.2;
                        
                        for w=1:size(keys.PSTH_WINDOWS,1)
                            sta=keys.PSTH_WINDOWS{w,2};
                            t_before_state=keys.PSTH_WINDOWS{w,3};
                            t_after_state=keys.PSTH_WINDOWS{w,4};
                            histo=[];
                            line_counter=0;
                            lines_for_average=[];
                            current_raster_y=0;
                            
                            for lin=unique_lines
                                line_counter=line_counter+1;
                                col=mod(lin-1,size(PSTH_perpos_colors,1))+1;
                                line_struct=subplot_struct([subplot_struct.line]==lin & [subplot_struct.figure]==fig);
                                n_trials=numel(line_struct);
                                state_shift=state_seperator(lin)-t_before_state;
                                if n_trials==0 || ~any([line_struct.states]==sta)
                                    state_seperator(lin)=state_shift + t_after_state + 0.1;
                                    continue;
                                else
                                    lines_for_average =[lines_for_average; lin];
                                end
                                
                                %% Raster and eye hand traces  + Number of trials line and text
                                line([state_shift+t_before_state state_shift+t_after_state],[current_raster_y current_raster_y],'Color',PSTH_perpos_colors(col,:),'LineWidth',0.5); hold on;
                                for t=1:numel(line_struct)
                                    if line_struct(t).hand==1
                                        hnd_col_v=keys.colors.lhd_ver;
                                        hnd_col_h=keys.colors.lhd_hor;
                                    elseif line_struct(t).hand==2
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
                                    ig_make_raster([line_struct(t).arrival_times(at_idx)]'+state_shift-trial_state_onset,current_raster_y-t,1,0,'Color',PSTH_perpos_colors(col,:),'LineWidth',keys.width.raster);
                                end
                                line([state_shift+t_before_state state_shift+t_after_state],[current_raster_y-t current_raster_y-t],'Color',PSTH_perpos_colors(col,:),'LineWidth',0.5);
                                
                                if w==1
                                    not_accepted=find(~[line_struct.accepted])*-1+current_raster_y;
                                    plot((state_shift+t_before_state)*ones(size(not_accepted)),not_accepted,'marker','s','markeredgecolor','k','markerfacecolor','k','linestyle','none');
                                    text(state_shift+t_before_state,current_raster_y,['N=' num2str(t)],'Color',PSTH_perpos_colors(col,:),'fontsize',5,'verticalalignment','top');
                                end
                                current_raster_y=current_raster_y-n_trials;
                                
                                % PSTH
                                state_seperator(lin)=state_shift + t_after_state + 0.1;
                                [histo(lin,:), bins, ~, SEM]=ph_spike_density(line_struct,w,keys,zeros(numel(line_struct),1),ones(numel(line_struct),1));
                                bins=bins+state_shift;
                                
                                lineProps={'color',PSTH_perpos_colors(col,:),'linewidth',keys.width.PSTH_perpos};
                                shadedErrorBar(bins,histo(lin,:),SEM,lineProps,1);
                                %line(bins,histo(lin,:),'color',PSTH_perpos_colors(col,:),'LineWidth',keys.width.PSTH_perpos);          %PSTH
                            end
                            if keys.plot.average_PSTH_line
                                line(bins,nanmean(histo(lines_for_average,:),1),'color','k','LineWidth',keys.width.PSTH_perpos);          %PSTH
                            end
                            state_seperator_max=max([state_seperator_max;state_seperator]);
                            y_max=max([y_max max(histo)]);
                            yL_min=min([yL_min current_raster_y]);
                        end
                    end
                    y_lim=[y_lim_PSTH,y_max];
                    for sub=subplot_indizes
                        subplot(o.rows,o.columns,sub,'XlimMode','Manual','YlimMode','Manual','XtickMode','Manual')
                        ph_PSTH_background(subplot_struct,y_lim,[yL_min y_lim(2)],[0 y_lim(2)],keys,0.5)
                        set(gca,'ylim',y_lim,'xlim',[0 state_seperator_max-0.1]);
                    end
                    title_and_save(PSTH_summary_handle,plot_1_title,[plot_1_title fig_title_part],keys);
                    
                    %% FR figure
                    plot_title        =[fig_title ' FR'];
                    FR_summary_handle   =figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',[plot_title fig_title_part]);
                    tr_index=of_index & [o.trial.effector]==effector & [o.trial.accepted]==1;
                    if sum(tr_index)==0 % might be the case if all trials are excluded cause of no spikes
                        continue;
                    end
                    % Computing mean firing rates per position
                    if keys.plot.average_heat_maps
                        o.line_labels=[o.line_labels {'Av'}];
                        unique_lines=[unique_lines max(unique_lines)+1];
                        o.PSTH_perpos_colors=[o.PSTH_perpos_colors; 0 0 0];
                    end
                    subplot_indizes=1:max([o.trial.subplot_pos]);                    
                    clear FR_heat;
                    
                    for sub=subplot_indizes
                        subplot_struct=o.trial(tr_index & [o.trial.subplot_pos]==sub);
                        for sta=all_sta
                            for lin=unique_lines
                                if keys.plot.average_PSTH_line && lin==max(unique_lines)
                                    line_struct=subplot_struct([subplot_struct.figure]==fig);
                                else
                                    line_struct=subplot_struct([subplot_struct.line]==lin & [subplot_struct.figure]==fig);
                                end
                                states_trials       =vertcat(line_struct.epoch);
                                if isempty(states_trials)% || isempty(subplot_struct.position) %unfortunate debugging for empty positions !
                                    FR_heat(sta,lin).pos(sub)    =NaN;
                                else
                                    FR_heat(sta,lin).pos(sub)    =nanmean([states_trials(:,sta).FR]);
                                end
                            end
                        end
                    end
                    
                    % Plotting heat maps
                    for sta=all_sta
                        for lin=unique_lines
                            state_label =keys.EPOCHS{sta,1};
                            line_label  =o.line_labels{lin};
                            hh(sta+n_states*(lin-1))=subplot_assignment(keys,'FR',n_states,max(unique_lines),sta,lin,e,0,effectors_on_figure,1);
                            plot_firing_rate_heatmap(FR_heat(sta,lin).pos,o.columns,o.rows,1:max(subplot_indizes));
                            title([line_label ' ' state_label]);
                        end
                        % one single colorbar
                        if sta==n_states && exist('hh','var')
                            hp4 = get(hh(end),'Position');
                            hc = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.015  hp4(4)*numel(unique_lines)]);
                            set(get(hc,'title'),'String','Sp/s');
                        end
                    end
                    ig_set_caxis_equal_lim(hh);
                    
                    %% Plotting polars
                    if keys.plot.polars_on_extra_figure
                        title_and_save(FR_summary_handle,plot_title,[plot_title fig_title_part],keys);
                        plot_title        =[fig_title ' polars'];
                        FR_summary_handle   =figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',[plot_title fig_title_part]);
                    end
                    clear FR_means FR_sterr;
                    polar_struct=o.trial(tr_index);
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
                        for lin=unique_lines
                            lin_idx=[polar_struct.line]==lin;
                            FR_means(sta,lin).vector=nanmean(sta_vector(lin_idx,:),1);
                            for ang=1:numel(unique_angles)
                                if keys.plot.average_PSTH_line && lin==max(unique_lines)
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
                    
                    FR_mean_max=max([FR_means.pos]+[FR_sterr.pos]); %+ [FR_sterr.pos]);
                    for sta=all_sta
                        for lin=unique_lines
                            subplot_assignment(keys,'polar',n_states,max(unique_lines),sta,lin,e,0,effectors_on_figure,1);
                            line_spec.color=o.PSTH_perpos_colors(lin,:);
                            line_spec.linewidth=1;
                            line_spec.visible='on';
                            ph_polar_plot(FR_means(sta,lin).pos,FR_sterr(sta,lin).pos,unique_angles',line_spec,FR_mean_max)
                            tmp_factor=FR_mean_max/norm(FR_means(sta,lin).vector);
                            line([0,FR_means(sta,lin).vector(1)]*tmp_factor,[0,FR_means(sta,lin).vector(2)]*tmp_factor,'linewidth',2,'color',line_spec.color);
                        end
                    end
                    title_and_save(FR_summary_handle,plot_title,[plot_title fig_title_part],keys);
                end
                
                %% Summary plot
                fig_title=sprintf('%s, %s, %s type, %s %s %s' ,keys.unit_ID, keys.target, type_string, keys.arrangement, title_part, title_value);
                plot_title        =[fig_title ' effector summary'];
                FR_summary_handle   =figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',[plot_title fig_title_part]);   
                hemifields=[o.trial.hemifield];
                unique_hemifields=unique(hemifields);
                unique_columns=unique([o.trial.column]);
                unique_rows=unique([o.trial.row]);
                side_labels={'L','R'}; 
                if keys.plot.average_heat_maps
                    unique_lines(end)=[];
                end
                if any(unique_hemifields==0) && fig==unique_figures(1)% vertical targets
                    side_labels={'L','V','R'};
                    n_colors=size(o.PSTH_summary_colors,1)/2;
                    V_col=linspace(0,0.5,n_colors)'*[1 1 1];
                    o.PSTH_summary_colors=[o.PSTH_summary_colors(1:n_colors,:);V_col;o.PSTH_summary_colors(n_colors+1:end,:)];
                end
                
                hs_ylim=0;
                hb_ylim=[0 0];
                hs=zeros(numel(unique_rows),numel(unique_columns));
                for r=1:numel(unique_rows)
                    row=unique_rows(r);
                    effectors_in_row=unique([o.trial(of_index & [o.trial.row]==row).effector]);
                    %% PSTH summary R/L and U/D subplot
                    for c=unique_columns
                        hs(r,c)=subplot_assignment(keys,'PSTH',n_states,numel(unique_lines),1,lin,r,c,unique_rows,unique_columns);
                        hold on
                        state_seperator=0;
                        for w=1:size(keys.PSTH_WINDOWS,1)
                            t_before_state=keys.PSTH_WINDOWS{w,3};
                            t_after_state=keys.PSTH_WINDOWS{w,4};
                            con_counter=0;
                            state_shift             =state_seperator-t_before_state;
                            bins                    =t_before_state:keys.PSTH_binwidth:t_after_state;
                            bins                    =bins+state_shift;
                            for h=1:numel(unique_hemifields)
                                for lin=1:max(unique_lines)
                                    con_counter=con_counter+1;
                                    tr_index=of_index & [o.trial.row]==row & [o.trial.column]==c & [o.trial.accepted]==1 & hemifields==unique_hemifields(h);
                                    if ~keys.plot.average_PSTH_line
                                        % this just makes it simpler, looping several times but plotting the same all over
                                        tr_index= tr_index &  [o.trial.line]==lin;
                                    end
                                    if sum(tr_index)==0
                                        continue;
                                    end
                                    [PSTH,~,~,SEM]=ph_spike_density(o.trial(tr_index),w,keys,zeros(sum(tr_index),1),ones(sum(tr_index),1));
                                    
                                    
                                    lineProps={'color',o.PSTH_summary_colors(con_counter,:),'linewidth',keys.width.PSTH_summary};
                                    shadedErrorBar(bins,PSTH,SEM,lineProps,1);
                                    %line(bins,PSTH,'color',o.PSTH_summary_colors(con_counter,:),'LineWidth',keys.width.PSTH_summary);
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
                            con_counter=0;
                            for h=1:numel(unique_hemifields)
                                for lin=1:max(unique_lines)
                                    con_counter=con_counter+1;
                                    tr_index=of_index & [o.trial.effector]==effector & [o.trial.row]==row & [o.trial.accepted]==1 & hemifields==unique_hemifields(h) & [o.trial.line]==lin;
                                    FR_labels{R,sta,con_counter}='';
                                    [FR_summ{R,sta,con_counter},FR_1quarters(R,sta,con_counter),FR_2quarters(R,sta,con_counter),FR_3quarters(R,sta,con_counter)]=deal(NaN);
                                    if sum(tr_index)==0
                                        continue;
                                    end
                                    per_epoch=vertcat(o.trial(tr_index).epoch);
                                    FR_summ{R,sta,con_counter}=[per_epoch(:,sta).FR];
                                    FR_1quarters(R,sta,con_counter)=quantile(FR_summ{R,sta,con_counter},0.25);
                                    FR_2quarters(R,sta,con_counter)=quantile(FR_summ{R,sta,con_counter},0.5);
                                    FR_3quarters(R,sta,con_counter)=quantile(FR_summ{R,sta,con_counter},0.75);
                                    FR_labels{R,sta,con_counter}=[side_labels{h} o.line_labels{lin}]; % get_label([FN{:}],hnd,ch);
                                    
                                    hb_ylim=[min(hb_ylim(1),max([FR_1quarters(R,sta,con_counter) nanmean(FR_summ{R,sta,con_counter})])) ...
                                        max(hb_ylim(2),max([FR_3quarters(R,sta,con_counter) nanmean(FR_summ{R,sta,con_counter})]))];
                                end
                            end
                            
                            state_label =keys.EPOCHS{sta,1};
                            hb(sta+numel(all_sta)*(e-1))=subplot_assignment(keys,'Bars',n_states,numel(unique_lines),sta,lin,r,c,unique_rows,unique_columns,er,numel(effectors_in_row));
                            tuning_labels.space     = ['in_' state_label '_spaceLR_' type_effector_short '_' keys.arrangement(1:3)];
                            tuning_labels.hands     = ['in_' state_label '_hands_' type_effector_short '_' keys.arrangement(1:3)];
                            tuning_labels.SXH       = ['in_' state_label '_SxH_' type_effector_short '_' keys.arrangement(1:3)];
                            tuning_labels.epoch     = ['in_AH_' state_label '_epoch_' type_effector_short '_' keys.arrangement(1:3)];
                            tuning_labels.choice    = ['ch_' state_label '_spaceLR_' type_effector_short '_' keys.arrangement(1:3)];
                            tuning_labels.LL_PT    =  ['in_LH_LS_' state_label '_PT_' type_effector_short '_' keys.arrangement(1:3)];
                            tuning_labels.RL_PT    =  ['in_LH_RS_' state_label '_PT_' type_effector_short '_' keys.arrangement(1:3)];
                            tuning_labels.LR_PT    =  ['in_RH_LS_' state_label '_PT_' type_effector_short '_' keys.arrangement(1:3)];
                            tuning_labels.RR_PT    =  ['in_RH_RS_' state_label '_PT_' type_effector_short '_' keys.arrangement(1:3)];
                            for FN=fieldnames(tuning_labels)'
                                idx= DAG_find_column_index(current_unit_tuning,tuning_labels.(FN{:}));
                                if ~isempty(idx)
                                    tuning.(FN{:})=current_unit_tuning{2,idx};
                                else
                                    tuning.(FN{:})='-';
                                end
                                if ~isstr(tuning.(FN{:}))
                                    tuning.(FN{:})=num2str(tuning.(FN{:}));
                                end
                            end
                            
                            title({state_label;sprintf('S %s E %s C %s H %s/%s',tuning.space,tuning.epoch,tuning.choice,tuning.hands,tuning.SXH);sprintf('PT: %s/%s/%s/%s',tuning.LL_PT,tuning.RL_PT,tuning.LR_PT,tuning.RR_PT)},'fontsize',8);
                            hold on
                            for b=1:size(FR_summ,3)
                                bm=nanmean(FR_summ{R,sta,b});
                                ci=[sterr(FR_summ{R,sta,b})*1.96 0];
                                ci=ci(1);
                                bar(b,bm,1,'FaceColor',o.PSTH_summary_colors(b,:));
                                eb=errorbar_constant_barsize_working(b, bm, ci, ci, 0.3);
                                patch([b-0.1 b+0.1 b+0.1 b-0.1],[FR_1quarters(R,sta,b) FR_1quarters(R,sta,b) FR_3quarters(R,sta,b) FR_3quarters(R,sta,b)],'k','EdgeColor','none');
                                scatter(b,FR_2quarters(R,sta,b),'o','markerfacecolor','w','markeredgecolor','k');
                                set(eb,'color','k');
                            end
                            set(gca,'xlim',[0.5 b+0.5],'Xtick',1:b,'Xticklabel',FR_labels(R,sta,:),'fontsize',8);
                            rotateXLabels(gca, 45);
                            
                            % vertical titles
                            if sta==1
                            ytitle={type_effector_full; sprintf('E:%s ExS:%s ExH:%s SxH:%s',...
                                anova_results(e).epoch, anova_results(e).ExS, anova_results(e).ExH, anova_results(e).SxH)};
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
                            con_counter=0;
                            current_raster_y=0;
                            for h=1:numel(unique_hemifields)
                                for lin=1:max(unique_lines)
                                    con_counter=con_counter+1;
                                    tr_index=of_index & [o.trial.row]==row & [o.trial.column]==c & [o.trial.accepted]==1 & hemifields==unique_hemifields(h);
                                    if ~keys.plot.average_PSTH_line
                                        % this just makes it simpler, looping several times but plotting the same all over
                                        tr_index= tr_index &  [o.trial.line]==lin;
                                    end
                                    n_trials     =sum(tr_index);
                                    line_struct=o.trial(tr_index);
                                    spike_length=y_lim(2)*0.002/(1+sum(of_index & [o.trial.accepted])*0.002);
                                    if n_trials==0 || ~any([line_struct.states]==sta)
                                    current_raster_y=current_raster_y-n_trials*spike_length;
                                        continue;
                                    end
                                    
                                    %% Raster
                                    line([state_shift+t_before_state state_shift+t_after_state],[current_raster_y current_raster_y],'Color',o.PSTH_summary_colors(con_counter,:),'LineWidth',0.5); hold on;
                                    for t=1:numel(line_struct)
                                        trial_state_onset=line_struct(t).states_onset(line_struct(t).states==sta);
                                        at_idx=line_struct(t).arrival_times-trial_state_onset>=t_before_state &...
                                            line_struct(t).arrival_times-trial_state_onset<=t_after_state;
                                        ig_make_raster([line_struct(t).arrival_times(at_idx)]'+state_shift-trial_state_onset,current_raster_y-t*spike_length,spike_length,0,'Color',o.PSTH_summary_colors(con_counter,:),'LineWidth',keys.width.raster);
                                    end
                                    
                                    line([state_shift+t_before_state state_shift+t_after_state],[current_raster_y-t*spike_length current_raster_y-t*spike_length],'Color',o.PSTH_summary_colors(con_counter,:),'LineWidth',0.5);
                                    if w==1
                                        text(state_shift+t_before_state,current_raster_y,['N=' num2str(t)],'Color',o.PSTH_summary_colors(con_counter,:),'fontsize',5,'verticalalignment','top');
                                    end
                                    current_raster_y=current_raster_y-n_trials*spike_length;
                                end
                            end
                            state_seperator=state_shift + t_after_state + 0.1;
                            y_lim(1)=min([y_lim(1), current_raster_y]);
                        end
                    end
                end
                for r=1:numel(unique_rows)
                    for c=unique_columns
                        subplot(hs(r,c))
                        ph_PSTH_background(o.trial,y_lim,y_lim,[0 y_lim(2)],keys,1)
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
                title_and_save(FR_summary_handle,plot_title,[plot_title fig_title_part],keys);
                
                %% here we stop plotting temporarily - the rest is movement related and doesnt work if movement directions are not assigned
                continue;
                
                %% per movement plot
                keys.movement_angle_binwidth=60;
                keys.movement_amplitude_binwidth=10;
                keys.movement_plot_type='revealed_or_not';
                keys.movement_plot_types={'pre','peri','post'};
                
                binwidth=keys.movement_angle_binwidth*pi/180;
                angular_movement_bins=-pi:binwidth:pi;
                
                for e=1:numel(effectors_on_figure)
                    effector=effectors_on_figure(e);
                    o_e=o;
                    o_e.trial=o_e.trial(of_index & [o.trial.effector]==effector);
                    [type_effector_full, type_effector_short] = MPA_get_type_effector_name(type,effector);
                    om=ph_arrange_movement_positions(population(unit).trial(oe_index),keys); %why not using o_e => this would separate into error and success
                    PSTH_perpos_colors=om.PSTH_perpos_colors;
                    
                    plot_4_title            = [fig_title  ' movement PSTHs '];
                    PSTH_movement_handle     = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',[plot_4_title]);
                    fig_title_part=[type_effector_short];
                    clear sp
                    for sub= unique([om.movement(~isnan([om.movement.subplot_pos])).subplot_pos])
                        sp(sub)=subplot(om.rows,om.columns,sub);
                        os=om.movement([om.movement.subplot_pos]==sub);
                        title([om.sub_title num2str(om.subplot_titles(sub))]);
                        hold on
                        current_raster_y=20;
                        unique_lines=unique([os.line]);
                        unique_lines(isnan(unique_lines))=[];
                        for l=1:numel(unique_lines)
                            lin=unique_lines(l);
                            line_idx=[os.line]==lin & arrayfun(@(x) sum(~isnan(x.spike_density))>0,os); %% unfortunate
                            to_plot=vertcat(os(line_idx).spike_density);
                            raster_to_plot={os(line_idx).arrival_times};
                            t_b=-0.2;
                            t_a= 0.2;
                            PSTH_x_axis=t_b:keys.PSTH_binwidth:t_a;
                            for t=1:numel(raster_to_plot)
                                ig_make_raster(raster_to_plot{t}',current_raster_y+t*0.1,1,0,'Color',PSTH_perpos_colors(l,:),'LineWidth',1)
                            end
                            current_raster_y=current_raster_y+t*0.1;
                            line([t_b t_a],[current_raster_y+1 current_raster_y+1],'Color',PSTH_perpos_colors(l,:),'LineWidth',0.5);
                            text(t_b,current_raster_y,['N=' num2str(t)],'Color',PSTH_perpos_colors(l,:),'fontsize',5);
                            line(PSTH_x_axis,nanmean(to_plot,1),'color',PSTH_perpos_colors(l,:));
                        end
                        y_lim_m(sub,:)=get(gca,'ylim');
                    end
                    legend(om.line_labels)
                    y_lim_m_max=max(max(y_lim_m));
                    y_lim_m_min=min(min(y_lim_m));
                    for sub= unique([om.movement(~isnan([om.movement.subplot_pos])).subplot_pos])
                        subplot(sp(sub));
                        rectangle('Position',[-0.15 y_lim_m_min 0.14 0.05*diff([y_lim_m_min y_lim_m_max])],'EdgeColor','none','FaceColor',[0.9 0.9 0.9]) %frame for indicating state separation
                        rectangle('Position',[-0.01 y_lim_m_min 0.06 0.05*diff([y_lim_m_min y_lim_m_max])],'EdgeColor','none','FaceColor',[0.6 0.6 0.6]) %frame for indicating state separation
                        rectangle('Position',[ 0.05 y_lim_m_min 0.1  0.05*diff([y_lim_m_min y_lim_m_max])],'EdgeColor','none','FaceColor',[0.3 0.3 0.3]) %frame for indicating state separation
                        
                        line([0 0],[y_lim_m_min y_lim_m_max],'color','k');
                        set(gca,'ylim',[y_lim_m_min y_lim_m_max]);
                    end
                    title_and_save(PSTH_movement_handle,plot_4_title,[plot_4_title fig_title_part],keys);
                    
                    
                    plot_5_title            = [fig_title  ' movement polar '];
                    fig_title_part=[type_effector_short];
                    polar_movement_handle     = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',[plot_5_title]);
                    keys.movement_plot_types={'pre','peri','post'};
                    row_fields={'no_target_revealed','target_revealed'};
                    n_columns=numel(keys.movement_plot_types);
                    % get unique shapes
                    oo=[population(unit).trial(oe_index)];
                    oo=[oo.movement];
                    unique_shapes=unique([[oo.shape_at_start] [oo.shape_at_end]]);
                    unique_shapes(isnan(unique_shapes))=[];
                    %clear PSTH_perpos_colors
                    % calculate mean and std of FR per angle
                    FR_rev_mean=NaN(numel(keys.movement_plot_types),1);
                    FR_rev_sem=NaN(numel(keys.movement_plot_types),1);
                    for n_column=1:numel(keys.movement_plot_types)
                        epoch=keys.movement_plot_types{n_column};
                        keys.movement_plot_type=epoch;
                        om=ph_arrange_movement_positions(population(unit).trial(oe_index),keys);
                        %why not using o_e => this would separate into error and success
                        PSTH_colors_cell{n_column}=om.PSTH_perpos_colors;
                        for n_row=1:numel(row_fields)
                            row_field=row_fields{n_row};
                            subplot(numel(row_fields)+1,n_columns,n_column+(n_row-1)*n_columns);
                            or=om.movement([om.movement.(row_field)]==1);
                            unique_lines=unique([or.line]);
                            %for lin=unique_lines
                            for pos=1:numel(angular_movement_bins)
                                pos_m=pos;
                                if pos==numel(angular_movement_bins)
                                    pos_m=1;
                                end
                                %m_index=[or.subplot_pos]==pos_m;% & [or.line]==lin;
                                m_index=[or.subplot_value]==angular_movement_bins(pos_m);% & [or.line]==lin;
                                FR_mean(l,pos,n_column,n_row)=nanmean([or(m_index).FR]);
                                FR_sem(l,pos,n_column,n_row)=sterr([or(m_index).FR]);
                            end
                            %end
                        end
                        for s=1:numel(unique_shapes)
                            shape=unique_shapes(s);
                            s_index=[or.shape_revealed]==shape;
                            FR_rev_mean(n_column,s)=nanmean([or(s_index).FR]);
                            FR_rev_sem(n_column,s)=sterr([or(s_index).FR]);
                        end
                    end
                    FR_mean_max=max(FR_mean(:));
                    FR_revealed_max=max(FR_rev_mean(:)+FR_rev_sem(:));
                    for n_column=1:numel(keys.movement_plot_types)
                        epoch=keys.movement_plot_types{n_column};
                        for n_row=1:numel(row_fields)
                            row_field=row_fields{n_row};
                            subplot(numel(row_fields)+1,n_columns,n_column+(n_row-1)*n_columns);
                            line_spec.visible='off';
                            polar(0,FR_mean_max,line_spec);%,'color',PSTH_perpos_colors(lin,:))
                            hold on
                            %unique_lines=unique([or.line]);
                            %for lin=unique_lines
                            line_spec.color=PSTH_colors_cell{n_column}(l,:);
                            line_spec.linewidth=2;
                            line_spec.visible='on';
                            ph_polar_plot(FR_mean(l,:,n_column,n_row),FR_sem(l,:,n_column,n_row),angular_movement_bins,line_spec)
                            %polar(angular_movement_bins,FR_mean(lin,:,n_column,n_row),line_spec);%,'color',PSTH_perpos_colors(lin,:))
                            hold on
                            %end
                            if n_row==1
                                title(epoch,'interpreter','none') %om.sub_title ' '
                            end
                            if n_column==1
                                ylabel(row_field,'interpreter','none')
                            end
                        end
                        
                        subplot(numel(row_fields)+1,n_columns,n_column+numel(row_fields)*n_columns);
                        hold on
                        bar(FR_rev_mean(n_column,:),'FaceColor','g');
                        erb=errorbar_constant_barsize_working([1:size(FR_rev_mean,2)]', FR_rev_mean(n_column,:), FR_rev_sem(n_column,:), FR_rev_sem(n_column,:), 0.3, 'o');
                        set(erb,'color','k')
                        %colors
                        set(gca,'xticklabel',unique_shapes,'xlim',[0.5 size(FR_rev_mean,2)+0.5],'ylim',[0 max([FR_revealed_max 1])]);
                        xlabel('shape','interpreter','none') %om.sub_title ' '
                        
                        if n_column==1
                            ylabel('FR');
                        end
                    end
                    title_and_save(polar_movement_handle,plot_5_title,[plot_5_title fig_title_part],keys);
                end
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

function title_and_save(figure_handle,filename,plot_title,keys)
mtit(figure_handle,  plot_title, 'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 8,'Interpreter', 'none');
stampit;
if keys.plot.export
    wanted_size=[50 30];
    set(figure_handle, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
    
    
    export_fig([keys.tuning_table_foldername filesep 'single_cell_examples' filesep filename], '-pdf','-transparent') % pdf by run
    close all
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


