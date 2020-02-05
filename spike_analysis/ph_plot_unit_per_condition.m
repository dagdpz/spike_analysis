function ph_plot_unit_per_condition(o,keys)
if isempty(o)
    return;
end
%keys.runpath_to_save=[keys.basepath_to_save keys.date filesep];
global MA_STATES
keys.saccade_states={'PreS','PeriS','PostS'};
keys.reach_states={'PreR','PeriR','PostR'};
%% colors and ylim preparation
PSTH_colors=o.PSTH_colors;
eye_offset=keys.max_firing_rate_expected+keys.max_n_trials_expected+keys.max_excentricity_expected*keys.eye_offset_multiplier(2);
hnd_offset=eye_offset+keys.max_excentricity_expected*keys.eye_offset_multiplier(2)*3;
eye_col_v=keys.eye_ver_color;
eye_col_h=keys.eye_hor_color;

if keys.hands_on && keys.plot_eye_hand_traces
    y_lim_PSTH=hnd_offset+keys.max_excentricity_expected*keys.eye_offset_multiplier(2);
elseif keys.plot_eye_hand_traces
    y_lim_PSTH=eye_offset+keys.max_excentricity_expected*keys.eye_offset_multiplier(2);
else
    y_lim_PSTH=keys.max_firing_rate_expected+keys.max_n_trials_expected;
end

%% selection of states to plot
epoch_to_plot_PSTH=ismember(keys.EPOCHS(:,1),keys.epochs_to_plot_PSTH{keys.type})';
states_to_plot=[keys.EPOCHS{:,2}];
[~,state_order_indexes]=ismember(states_to_plot,MA_STATES.all_states);
[sta_index_in_order,all_sta]=sort(state_order_indexes);
sta_in_order=all_sta(sta_index_in_order>0 & epoch_to_plot_PSTH(all_sta));
n_states=numel(all_sta);
tmp=[o(1).position.trial];
unique_figures=unique([tmp.figure]);
get_expected_MA_states(keys.type,keys.effector,keys.effectors_on_same_figure);

% 
%     get_expected_MA_states(conditions(c).type,conditions(c).effector,keys.effectors_on_same_figure);
%     %get_expected_MA_states(conditions(c).type,conditions(c).effector);
%     keys.EPOCHS=keys.EPOCHS_PER_TYPE{conditions(c).type};
%     %keys.EPOCHS= keys.EPOCHS(ismember(keys.EPOCHS(:,1),keys.epochs_to_plot_PSTH{conditions(c).type}),:);
%     epoch_to_plot_PSTH=ismember(keys.EPOCHS(:,1),keys.epochs_to_plot_PSTH{conditions(c).type})';
%     states_to_plot=[keys.EPOCHS{:,2}];
%     [~,state_order_indexes]=ismember(states_to_plot,MA_STATES.all_states);
%     [sta_index_in_order,sta_in_order]=sort(state_order_indexes);
%     sta_in_order=sta_in_order(sta_index_in_order>0 & epoch_to_plot_PSTH(sta_in_order));


for fig=unique_figures
    for e=1:numel(o) %% looping through effectors (potentially)
        if fig > numel(o(e).figure_title_value)
            continue;
        end
        title_value=vertcat(o(e).figure_title_value{fig});
        title_part=o(e).figure_title_part;
        [type_effector_full, type_effector_short] = get_type_effector_name(keys.type,keys.effector(e));
        
        anova_labels.epoch  = ['in_epoch_main_' type_effector_short '_' keys.plot_type(1:3)];
        anova_labels.space  = ['in_spaceLR_main_'  type_effector_short '_' keys.plot_type(1:3)];
        anova_labels.choice = ['ch_spaceLR_main_'  type_effector_short '_' keys.plot_type(1:3)];
        anova_labels.hands  = ['in_hands_main_' type_effector_short '_' keys.plot_type(1:3)];
        anova_labels.ExS    = ['in_ExS_' type_effector_short '_' keys.plot_type(1:3)];
        anova_labels.ExH    = ['in_ExH_' type_effector_short '_' keys.plot_type(1:3)];
        anova_labels.SxH    = ['in_SxH_' type_effector_short '_' keys.plot_type(1:3)];
        for FN=fieldnames(anova_labels)'
            idx= find_column_index(keys.current_unit_tuning,anova_labels.(FN{:}));
            if ~isempty(idx)
                anova_results(e).(FN{:})=keys.current_unit_tuning{2,idx};
            else
                anova_results(e).(FN{:})='-';
            end
            if ~isstr(anova_results(e).(FN{:}))
                anova_results(e).(FN{:})=num2str(anova_results(e).(FN{:}));
            end
        end
        
        
        fig_title=sprintf('%s, %s, %s, %s %s %s' ,...
            keys.unit_ID, keys.target, type_effector_full, keys.plot_type, title_part, title_value);
        fig_title_part=sprintf(', Stability %d, SNR %d, Single %d Grid: %d/%d Depth %.2f ANOVA E:%s S:%s C:%s H:%s ExS:%s ExH:%s SxE:%s Channel: %d Blocks&Units: %s ', ...
            keys.stability_rating, keys.SNR_rating, keys.Single_rating, keys.grid_x, keys.grid_y, keys.electrode_depth, ...
            anova_results(e).epoch, anova_results(e).space, anova_results(e).choice, anova_results(e).hands,...
            anova_results(e).ExS, anova_results(e).ExH, anova_results(e).SxH, keys.channel, [keys.block_unit{:}]);
        
        
        
        %% PSTH figure
        if any(ismember(keys.summary,[-1,1]))
            plot_1_title            = [fig_title  ' PSTHs'];
            PSTH_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',[plot_1_title fig_title_part]);
            background_plotted=zeros(numel(o(e).position),numel(states_to_plot));
            %set(gca,'XlimMode','Manual','YlimMode','Manual','XtickMode','Manual');
            for sub=1:numel(o(e).position)
                subplot(o(e).rows,o(e).columns,sub,'XlimMode','Manual','YlimMode','Manual','XtickMode','Manual')
                hold on
                subplot_struct=o(e).position(sub);
                y_lim=[0,y_lim_PSTH];
                state_seperator_max=0.2;
                if isempty(subplot_struct.trial);
                    continue;
                end;
                
                title([subplot_struct.title_part num2str(round(subplot_struct.position))]);
                unique_lines=unique([subplot_struct.trial.line]);
                state_seperator=zeros(max(unique_lines),1);
                
                for sta_sp=1:numel(sta_in_order)
                    sta=sta_in_order(sta_sp);
                    % epoch definitions
                    state_label=keys.EPOCHS{sta,1};
                    t_before_state=keys.EPOCHS{sta,3};
                    t_after_state=keys.EPOCHS{sta,4};
                    patch_x=[keys.EPOCHS{sta,5} keys.EPOCHS{sta,6}-keys.EPOCHS{sta,5}];
                    state_onsets_tmp=[];
                    clear histo
                    line_counter=0;
                    lines_for_average=[];
                    
                    for lin=unique_lines
                        
                        line_counter=line_counter+1;
                        col=mod(lin-1,size(PSTH_colors,1))+1;
                        per_state=subplot_struct.per_state(sta);
                        per_state.trial=per_state.trial([per_state.trial.line]==lin & [per_state.trial.figure]==fig);
                        n_trials=numel(per_state.trial);
                        if n_trials==0
                            continue;
                        end
                        
                        %% y definition (for rasters) and plotting shaded areas below everything else
                        state_shift=state_seperator(lin)-t_before_state;
                        if ~background_plotted(sub,sta)
                            current_raster_y=keys.max_firing_rate_expected;
                            rectangle('Position',[patch_x(1)+state_shift y_lim(1) patch_x(2) diff(y_lim)],'EdgeColor','none','FaceColor',[0.9 0.9 0.9]) %frame for indicating state separation
                            rectangle('Position',[state_seperator(lin) y_lim(1) t_after_state-t_before_state diff(y_lim)]) %frame for indicating state separation
                            if state_shift>=state_seperator(lin) && state_shift<= state_seperator(lin)+t_after_state-t_before_state
                                line([state_shift state_shift],y_lim,'color',[0.7 0.7 0.7]);
                            end
                            text(state_seperator(lin),(diff(y_lim)*-0.04)/2,state_label,'fontsize',7);                               %label for state
                            background_plotted(sub,sta)=1;
                        end
                        
                        %% Raster and eye hand traces  + Number of trials line and text
                        line([state_shift+t_before_state state_shift+t_after_state],[current_raster_y+1 current_raster_y+1],'Color',PSTH_colors(col,:),'LineWidth',0.5); hold on;
                        
                        for t=1:numel(per_state.trial)
                            if per_state.trial(t).hand==1
                                hnd_col_v=keys.lhnd_ver_color;
                                hnd_col_h=keys.lhnd_hor_color;
                            elseif per_state.trial(t).hand==2
                                hnd_col_v=keys.rhnd_ver_color;
                                hnd_col_h=keys.rhnd_hor_color;
                            else
                                hnd_col_v=[1 1 1];
                                hnd_col_h=[1 1 1];
                            end
                            time_axis=per_state.trial(t).time_axis+state_shift;
                            if keys.plot_eye_hand_traces
                                line(time_axis,per_state.trial(t).x_eye*keys.eye_offset_multiplier(2) + eye_offset ,'color',eye_col_h);
                                line(time_axis,per_state.trial(t).y_eye*keys.eye_offset_multiplier(2) + eye_offset,'color',eye_col_v);
                                line(time_axis,per_state.trial(t).x_hnd*keys.hnd_offset_multiplier(2) + hnd_offset,'color',hnd_col_h);
                                line(time_axis,per_state.trial(t).y_hnd*keys.hnd_offset_multiplier(2) + hnd_offset,'color',hnd_col_v);
                            end
                            ig_make_raster([per_state.trial(t).arrival_times]'+state_shift,current_raster_y+t,1,0,'Color',PSTH_colors(col,:),'LineWidth',keys.l_width_ras_s1);
                        end
                        
                        line([state_shift+t_before_state state_shift+t_after_state],[current_raster_y+t+1 current_raster_y+t+1],'Color',PSTH_colors(col,:),'LineWidth',0.5);
                        if sta==sta_in_order(1)
                            text(state_shift+t_before_state,current_raster_y,['N=' num2str(t)],'Color',PSTH_colors(col,:),'fontsize',5);
                        end
                        
                        current_raster_y=current_raster_y+n_trials;
                        
                        %% PSTH
                        state_seperator(lin)=state_shift + t_after_state + 0.1;
                        state_onsets_tmp=[state_onsets_tmp per_state.trial.state_onset]; %% THIS IS INCORRECT!!!
                        state_shifts(sta_sp)=state_shift;
                        bins=t_before_state:keys.PSTH_binwidth:t_after_state;
                        bins=bins+state_shift;
                        if ~all(isempty(vertcat(per_state.trial.spike_density)))
                            lines_for_average=[lines_for_average lin];
                            histo(lin,:)=sum(vertcat(per_state.trial.spike_density),1)/n_trials;
                            line(bins,histo(lin,:),'color',PSTH_colors(col,:),'LineWidth',keys.l_width_his_s1);          %PSTH
                        end
                        
                    end
                    
                    if keys.plot_average_line
                        line(bins,nanmean(histo(lines_for_average,:)),'color',keys.color_average_line,'LineWidth',keys.l_width_his_s1);          %PSTH
                    end
                    state_seperator_max=max([state_seperator_max;state_seperator]);
                    text(state_shift+t_before_state,diff(y_lim)*-0.1,sprintf('%.0f',t_before_state*1000),'fontsize',5);
                    text(state_shift+t_after_state,diff(y_lim)*-0.1,sprintf('%.0f',t_after_state*1000),'Horizontalalignment','right','fontsize',5);
                    state_onsets(sta_sp)=nanmean(state_onsets_tmp); %% THIS IS INCORRECT!!!
                end
                
                set(gca,'ylim',y_lim,'xlim',[0 state_seperator_max-0.1]);
                [u_state_onsets,u_state_onset_indexes,~]=unique(state_onsets);
                text(state_shifts(u_state_onset_indexes),diff(y_lim)*-0.15*ones(size(u_state_onset_indexes)),num2str(round(u_state_onsets'*1000)),'fontsize',6,'horizontalAlignment','center','color',[0.7 0.7 0.7]);
                remove_axis('x')
            end
            %end
            title_and_save(PSTH_summary_handle,plot_1_title,[plot_1_title fig_title_part],keys);
        end
        
        
        %% FR figure
        if any(ismember(keys.summary,[-1,2]))
            %for e=1:numel(o)
            plot_2_title        =[fig_title ' FR'];
            FR_summary_handle   =figure('units','normalized','outerposition',[0 0 1 1],'name',[plot_2_title fig_title_part]);
            
            %% Computing mean firing rates per position
            postrials=[o(e).position.trial];
            if keys.average_heat_maps
                unique_lines=1;
            else
                unique_lines=unique([postrials.line]);
            end
            clear FR;
            for sub=1:numel(o(e).position)
                subplot_struct  =o(e).position(sub);
                for sta=all_sta
                    for lin=unique_lines
                        if isempty(subplot_struct.trial) || isempty(subplot_struct.position) %unfortunate debugging for empty positions !
                            FR(sta,lin).pos(sub)    =NaN;
                        else
                            per_state=subplot_struct.per_state(sta);
                            if keys.plot_average_line
                                per_state.trial=per_state.trial([per_state.trial.figure]==fig);
                            else
                                per_state.trial=per_state.trial([per_state.trial.line]==lin & [per_state.trial.figure]==fig);
                            end
                            states_trials       =[per_state.trial];
                            FR(sta,lin).pos(sub)    =nanmean([states_trials.FR]);
                        end
                    end
                end
            end
            
%             y_lim       =[0,max([FR.pos])];
%             max_duration=max([keys.EPOCHS{all_sta,4}]-[keys.EPOCHS{all_sta,3}]);
%             
            %% Plotting Firing rates
            for sta_sp=1:numel(all_sta)
                sta=all_sta(sta_sp);
                hs_ylim(sta)=0;
                hb_ylim=[0 0];
                unique_lines=1:size(FR,2);
                for lin=unique_lines
                    state_label =keys.EPOCHS{sta,1};
                    line_label  =o(e).line_labels{lin};
                    hh(sta_sp+n_states*(lin-1))=subplot_assignment('FR',n_states,numel(unique_lines),sta_sp,keys,lin,e);
                    plot_firing_rate_heatmap(FR(sta,lin).pos,o(e).columns,o(e).rows,1:numel(o(e).position));
                    title([line_label ' ' state_label]);
                end
                
                % one single colorbar
                if sta==n_states
                    hp4 = get(hh(end),'Position');
                    hc = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.015  hp4(4)*numel(unique_lines)]);
                    set(get(hc,'title'),'String','Sp/s');
                end
            end
            
            ig_set_caxis_equal_lim(hh);
            if keys.plot_FR_separately; title_and_save(FR_summary_handle,plot_2_title,[plot_2_title fig_title_part],keys); end
            
        end
    end
    
    if any(ismember(keys.summary,[-1,2]))
        if keys.effectors_on_same_figure; fig_title=sprintf('%s, %s, effector summary, %s %s %s' ,keys.unit_ID, keys.target, keys.plot_type, title_part, title_value); end
        plot_3_title        =[fig_title ' summary'];
        if keys.plot_FR_separately; FR_summary_handle   =figure('units','normalized','outerposition',[0 0 1 1],'name',[plot_3_title fig_title_part]); end
        
        for e=1:numel(o)
            
            %% Summary plot preparations
            fields_to_plot_PSTH={'left','right'};
            unique_hands=[];
            unique_choices=[];
            for fn=1:numel(fields_to_plot_PSTH)
                FN=fields_to_plot_PSTH{fn};
                fig_idx=[o(e).(FN).trial.figure]==fig;
                unique_hands=unique([unique_hands unique([o(e).(FN).trial(fig_idx).hand])]);
                unique_choices=unique([unique_choices unique([o(e).(FN).trial(fig_idx).choice])]);
            end
            %con_col=ones(size(o.PSTH_summary_colors,1));
            
            
            
            state_seperator=0;
            for sta_sp=1:numel(sta_in_order)
                %% PSTH summary R/L and U/D subplot
                hs(e)=subplot_assignment('PSTH',n_states,numel(unique_lines),1,keys,lin,e);
                hold on
                sta=sta_in_order(sta_sp);
                hold on
                t_before_state  =keys.EPOCHS{sta,3};
                t_after_state   =keys.EPOCHS{sta,4};
                
                con_counter=0;
                state_shift             =state_seperator-t_before_state;
                state_shifts(sta_sp)    =state_shift;
                bins                    =t_before_state:keys.PSTH_binwidth:t_after_state;
                bins                    =bins+state_shift;
                state_onsets(sta_sp)=0;
                for fn=1:numel(fields_to_plot_PSTH)
                    FN=fields_to_plot_PSTH{fn};
                    if ~isfield(o(e).(FN),'per_state') %unfortunate debugging for empty positions !!!!
                        con_counter=con_counter+numel(unique_hands)*numel(unique_choices);
                        %col_counter=numel(o.PSTH_summary_colors)/numel(fields_to_plot_PSTH);
                        continue;
                    end
                    for hnd=unique_hands
                        for ch=unique_choices
                            con_counter=con_counter+1;
                            
                            tr_index=true(size([o(e).(FN).per_state(sta).trial]));
                            tr_index= tr_index &  [o(e).(FN).per_state(sta).trial.figure]==fig;
                            if numel(unique_choices)>1
                                tr_index=tr_index &  [o(e).(FN).per_state(sta).trial.choice]==ch ;
                            end
                            if numel(unique_hands)>1
                                tr_index=tr_index & [o(e).(FN).per_state(sta).trial.hand]==hnd;
                            end
                            n_trials     =sum(tr_index);
                            if n_trials==0
                                continue;
                            end
                            current_trials=o(e).(FN).per_state(sta).trial(tr_index);
                            PSTH=sum(vertcat(current_trials.spike_density),1)/n_trials;
                            if isempty(PSTH)
                                continue;
                            end
                            line(bins,PSTH,'color',o(e).PSTH_summary_colors(con_counter,:),'LineWidth',keys.l_width_his_s2);
                            
                            state_onsets(sta_sp)=max([state_onsets(sta_sp), nanmean([current_trials.state_onset])]); %% THIS IS INCORRECT!!!
                            hs_ylim(sta)=max(hs_ylim(sta),max(PSTH));
%                             hb_ylim=[   min(hb_ylim(1),max([FR_1quarters(sta,con_counter) mean(FR_summ{sta,con_counter})])) ...
%                                 max(hb_ylim(2),max([FR_3quarters(sta,con_counter) mean(FR_summ{sta,con_counter})]))];
                        end
                    end
                end
                %set(gca,'Layer','top','xlim',[t_before_state t_before_state+max_duration],'fontsize',8);
                state_seperator=state_shift + t_after_state + 0.1;
            end
            
            %% Bar plots
            for sta_sp=1:numel(all_sta)
                %% PSTH summary bar plots
                sta=all_sta(sta_sp);
                con_counter=0;
                %FR_grouping_vector{sta}=[];
                for fn=1:numel(fields_to_plot_PSTH)
                    FN=fields_to_plot_PSTH{fn};
                    if ~isfield(o(e).(FN),'per_state') %unfortunate debugging for empty positions !!!!
                        con_counter=con_counter+numel(unique_hands)*numel(unique_choices);
                        continue;
                    end
                    
                    mean_state_onsets(e,sta_sp).(FN)=nanmean(vertcat(o(e).(FN).per_state(sta).trial([o(e).(FN).per_state(sta).trial.figure]==fig).state_onsets),1);
                    for hnd=unique_hands
                        for ch=unique_choices
                            con_counter=con_counter+1;
                            
                            FR_labels{sta,con_counter}='';
                            [FR_summ{sta,con_counter},FR_1quarters(sta,con_counter),FR_2quarters(sta,con_counter),FR_3quarters(sta,con_counter)]=deal(NaN);
                            
                            tr_index=true(size([o(e).(FN).per_state(sta).trial]));
                            tr_index= tr_index &  [o(e).(FN).per_state(sta).trial.figure]==fig;
                            if numel(unique_choices)>1
                                tr_index=tr_index &  [o(e).(FN).per_state(sta).trial.choice]==ch ;
                            end
                            if numel(unique_hands)>1
                                tr_index=tr_index & [o(e).(FN).per_state(sta).trial.hand]==hnd;
                            end
                            n_trials     =sum(tr_index);
                            if n_trials==0
                                continue;
                            end
                            current_trials=o(e).(FN).per_state(sta).trial(tr_index);
                            %% for bar plots
                            FR_summ{sta,con_counter}=[current_trials.FR];
                            FR_1quarters(sta,con_counter)=quantile(FR_summ{sta,con_counter},0.25);
                            FR_2quarters(sta,con_counter)=quantile(FR_summ{sta,con_counter},0.5);
                            FR_3quarters(sta,con_counter)=quantile(FR_summ{sta,con_counter},0.75);
                            FR_labels{sta,con_counter}=get_label(FN,hnd,ch);
                            
                            hb_ylim=[   min(hb_ylim(1),max([FR_1quarters(sta,con_counter) mean(FR_summ{sta,con_counter})])) ...
                                max(hb_ylim(2),max([FR_3quarters(sta,con_counter) mean(FR_summ{sta,con_counter})]))];
                        end
                    end
                end
                %set(gca,'Layer','top','xlim',[t_before_state t_before_state+max_duration],'fontsize',8);
                state_seperator=state_shift + t_after_state + 0.1;
            end
            
            for sta_sp=1:numel(all_sta)
                sta=all_sta(sta_sp);
                state_label =keys.EPOCHS{sta,1};
                hb(sta+numel(all_sta)*(e-1))=subplot_assignment('Bars',n_states,numel(unique_lines),sta_sp,keys,lin,e);
                [type_effector_full, type_effector_short] = get_type_effector_name(keys.type,keys.effector(e));
                tuning_labels.space     = ['in_' state_label '_spaceLR_' type_effector_short '_' keys.plot_type(1:3)];
                tuning_labels.hands     = ['in_' state_label '_hands_' type_effector_short '_' keys.plot_type(1:3)];
                tuning_labels.SXH       = ['in_' state_label '_SxH_' type_effector_short '_' keys.plot_type(1:3)];
                tuning_labels.epoch     = ['in_NH_' state_label '_epoch_' type_effector_short '_' keys.plot_type(1:3)];
                %tuning_labels.epochR   = ['in_R_' state_label '_' type_effector_short '_' keys.plot_type(1:3)];
                tuning_labels.choice    = ['ch_' state_label '_spaceLR_' type_effector_short '_' keys.plot_type(1:3)];
                for FN=fieldnames(tuning_labels)'
                    idx= find_column_index(keys.current_unit_tuning,tuning_labels.(FN{:}));
                    if ~isempty(idx)
                        tuning.(FN{:})=keys.current_unit_tuning{2,idx};
                    else
                        tuning.(FN{:})='-';
                    end
                    if ~isstr(tuning.(FN{:}))
                        tuning.(FN{:})=num2str(tuning.(FN{:}));
                    end
                end
                
                title({state_label;sprintf('S %s E %s C %s H %s/%s',tuning.space,tuning.epoch,tuning.choice,tuning.hands,tuning.SXH)},'fontsize',8);
                hold on
                for b=1:size(FR_summ,2)
                    bm(sta,b)=nanmean(FR_summ{sta,b});
                    ci=[sem(FR_summ{sta,b})*1.96 0];
                    ci=ci(1);
                    bar(b,bm(sta,b),1,'FaceColor',o(e).PSTH_summary_colors(b,:));
                    eb=errorbar_constant_barsize_working(b, bm(sta,b), ci, ci, 0.3);
                    patch([b-0.1 b+0.1 b+0.1 b-0.1],[FR_1quarters(sta,b) FR_1quarters(sta,b) FR_3quarters(sta,b) FR_3quarters(sta,b)],'k','EdgeColor','none');
                    %rectangle('Position',[b-0.1 FR_1quarters(sta,b) 0.2 FR_3quarters(sta,b)-FR_1quarters(sta,b)],'FaceColor','k','EdgeColor','none');
                    scatter(b,FR_2quarters(sta,b),'o','markerfacecolor','w','markeredgecolor','k');
                    set(eb,'color','k');
                end
                set(gca,'xlim',[0.5 b+0.5],'Xtick',1:b,'Xticklabel',FR_labels(sta,:),'fontsize',8);
                rotateXLabels(gca, 90);
            end
            
            
            %% vertical titles
            
            subplot(hb(all_sta(1)+numel(all_sta)*(e-1)));
            ytitle={type_effector_full; sprintf('E:%s S:%s C:%s H:%s ExS:%s ExH:%s SxH:%s',...
                anova_results(e).epoch, anova_results(e).space, anova_results(e).choice, anova_results(e).hands,...
                anova_results(e).ExS, anova_results(e).ExH, anova_results(e).SxH)};
            text(-2,0,ytitle,'rotation',90,'interpreter','none');
        end
        
        %% mean state onsets (combining the left and right subfields)
        
           
        state_onsets_to_plot=cell(numel(o),numel(all_sta));
        for e=1:size(state_onsets_to_plot,1)
            for s=1:size(state_onsets_to_plot,2)
                state_onsets_to_plot{e,s}=mean([mean_state_onsets(e,s).left; mean_state_onsets(e,s).right],1);
            end
        end
        
        %% equalizing summary PSTH axes, adding background and text
        y_lim=[0, max(max(hs_ylim))];
        for e=1:numel(hs)
            subplot(hs(e))
            state_seperator         =0;
            state_seperator_max     =0;
            for sta_sp=1:numel(sta_in_order)
                sta=sta_in_order(sta_sp);
                t_before_state  =keys.EPOCHS{sta,3};
                t_after_state   =keys.EPOCHS{sta,4};
                state_shift     =state_seperator-t_before_state;
                
                
                %% background
                for sta_sta=1:numel(all_sta)
                    sta_t=all_sta(sta_sta);
                    if ~ismember(keys.EPOCHS{sta,2},MA_STATES.all_states) || ~ismember(keys.EPOCHS{sta_t,2},MA_STATES.all_states)
                       continue; 
                    end
                    relative_state_onset=state_onsets_to_plot{e,sta};
                    state_shift_relative=state_shift+relative_state_onset(sta_t);
                    state_label     =keys.EPOCHS{sta_t,1};
                    t_b             =keys.EPOCHS{sta_t,5};
                    t_a             =keys.EPOCHS{sta_t,6};
                    if state_shift_relative>=state_seperator-0.0001 && state_shift_relative <= state_seperator+t_after_state-t_before_state+0.0001
                        if ismember(state_label,keys.saccade_states)
                            li=line([state_shift_relative state_shift_relative],y_lim,'color',[1 0 0]);
                            ptch=rectangle('Position',[t_b+state_shift_relative y_lim(1) t_a-t_b diff(y_lim)*0.05],'EdgeColor','none','FaceColor',[1 0 0]); %frame for indicating FR window
                        elseif ismember(state_label,keys.reach_states)
                            li=line([state_shift_relative state_shift_relative],y_lim,'color',[0 1 0]);
                            ptch=rectangle('Position',[t_b+state_shift_relative y_lim(1) t_a-t_b diff(y_lim)*0.05],'EdgeColor','none','FaceColor',[0 1 0]); %frame for indicating FR window
                        else
                            li=line([state_shift_relative state_shift_relative],y_lim,'color',[0.8 0.8 0.8]);
                            ptch=rectangle('Position',[t_b+state_shift_relative y_lim(1) t_a-t_b diff(y_lim)*0.05],'EdgeColor','none','FaceColor',[0.8 0.8 0.8]); %frame for indicating FR window
                        end
                        
                        uistack(ptch, 'bottom');
                        text(double(state_shift_relative+t_b),diff(y_lim)*-0.03,state_label,'fontsize',12);  %label for state
                    end
                end
                 
                rectangle('Position',[state_seperator y_lim(1) t_after_state-t_before_state diff(y_lim)]) %frame for indicating state separation
                 
                 
                state_seperator=state_shift + t_after_state + 0.1;
                state_seperator_max=max([state_seperator_max;state_seperator]);
                text(state_shift+t_before_state,diff(y_lim)*-0.07,sprintf('%.0f',t_before_state*1000),'fontsize',8);
                text(state_shift+t_after_state,diff(y_lim)*-0.07,sprintf('%.0f',t_after_state*1000),'Horizontalalignment','right','fontsize',8);
            end
            set(gca,'ylim',y_lim,'xlim',[0 state_seperator_max-0.1]);
            [u_state_onsets,u_state_onset_indexes,~]=unique(state_onsets);
            text(state_shifts(u_state_onset_indexes),diff(y_lim)*-0.1*ones(size(u_state_onset_indexes)),num2str(round(u_state_onsets'*1000)),'fontsize',8,'horizontalAlignment','center','color',[0.7 0.7 0.7]);
            remove_axis('x')
        end
        for hb_idx=find(hb~=0)
            subplot(hb(hb_idx))
            set(gca,'ylim',[hb_ylim]);
        end
        
        title_and_save(FR_summary_handle,plot_3_title,[plot_3_title fig_title_part],keys);
    end
    
end
end

function label=get_label(FN,hnd,ch)
switch hnd
    case 0
        hnd_str='';
    case 1
        hnd_str='LH';
    case 2
        hnd_str='RH';
end
switch ch
    case 0
        cho_str='I';
    case 1
        cho_str='C';
end

label=[upper(FN(1)) hnd_str cho_str];
end

function plot_firing_rate_heatmap(R,n_columns,n_rows,subplot_pos)
X =NaN(n_columns,n_rows);
X(subplot_pos(~isnan(subplot_pos)))=R;
X=X';
map = jet;
colormap(map);
X = [[X nan(size(X,1),1)] ; nan(1,size(X,2)+1)];
pcolor(X); set(gca,'Ydir','reverse');
axis equal; axis off;
end

function title_and_save(figure_handle,filename,plot_title,keys)
mtit(figure_handle,  plot_title, 'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none');
stampit;
if keys.create_pdfs
    switch keys.append_pdfs
        case 0
            export_fig([keys.path_to_save 'single_cell_examples' filesep filename], '-pdf','-transparent') % pdf by run
        case 1
            export_fig([keys.path_to_save 'single_cell_examples' filesep filename '_appended'], '-pdf', '-append','-transparent') % pdf by run
    end
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

function handle=subplot_assignment(subplot_case,n_states,n_lines,sta_sp,keys,lin,effector)
if keys.plot_FR_separately
    n_rows_FR           =n_lines;
    n_rows_PSTH         =2*numel(keys.effector);
    FR_row_shift        =n_states*(lin-1);
else
    n_rows_FR           =n_lines + 2;
    n_rows_PSTH         =n_lines + 2*numel(keys.effector);
    FR_row_shift        =n_states*(lin+1);
end
switch subplot_case
    case 'PSTH'
        handle=subplot(n_rows_PSTH,1,2*effector-1);
    case 'Bars'
        handle=subplot(n_rows_PSTH,n_states,sta_sp+n_states*(2*effector-1));
    case 'FR'
        handle=subplot(n_rows_FR,n_states,sta_sp+FR_row_shift);
end
end
