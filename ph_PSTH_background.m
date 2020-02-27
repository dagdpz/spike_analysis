function ph_PSTH_background(trial,y_lim,y_limframe,y_limPSTH,keys,fontsize_factor)
if nargin<5
    fontsize_factor=1;
end
if diff(y_lim)==0
    y_lim(2)=y_lim(1)+1;
end
if diff(y_limframe)==0
    y_limframe(2)=y_limframe(1)+1;
end
if diff(y_limPSTH)==0
    y_limPSTH(2)=y_limPSTH(1)+1;
end
state_seperator         =0;
state_seperator_max     =0;
for w=1:size(keys.PSTH_WINDOWS,1)
    sta             =keys.PSTH_WINDOWS{w,2};
    t_before_state  =keys.PSTH_WINDOWS{w,3};
    t_after_state   =keys.PSTH_WINDOWS{w,4};
    window_label    =keys.PSTH_WINDOWS{w,1};
    state_shift     =state_seperator-t_before_state;
    
    %% background
    [state_names,absolute_state_onsets,relative_state_onset,relative_epochs,epoch_names,states]=ph_state_onsets(trial,sta,keys);
    states_in_window=relative_state_onset>=t_before_state & relative_state_onset<=t_after_state;
    epochs_in_window=relative_epochs(:,1)>=t_before_state & relative_epochs(:,2)<=t_after_state ;
    epochs_in_window(~ismember(epoch_names,keys.ANOVAS.main))=false;
    states_in_window(~ismember(states,keys.plot.events))=false;
    
    % plot only desired events and epochs
    state_positions=relative_state_onset(states_in_window)+state_shift;
    if ~isempty(state_positions)
        line([state_positions; state_positions],y_limframe,'color',[0.8 0.8 0.8],'linewidth',1);
        %here we loop, just so not everything is one text object later, which is sort of a bug in export_fig
        temp_state_onsets=absolute_state_onsets(states_in_window);
        for p=1:numel(state_positions)
            text(double(state_positions(p)),diff(y_limframe)*-0.07+y_limframe(1),num2str(round(temp_state_onsets(p)*1000)),...
                'fontsize',8*fontsize_factor,'horizontalalignment','center','color',[p p p]/255); %this is a trick: becuase they are plotted in different colors, they will be diffrent elements in adobe
        end
        text(double(state_positions'),ones(size(state_positions'))*(y_limframe(2)-diff(y_limframe)*0.15),state_names(states_in_window),...
            'fontsize',6*fontsize_factor,'interpreter','none','rotation',90,'verticalalignment','top');  %label for state
        line([state_shift state_shift],y_limframe,'color',[1 0 1],'linewidth',1);
    end
    text(state_seperator+(t_after_state-t_before_state)/2,y_limPSTH(2)-diff(y_limPSTH)*0.05,window_label,...
        'color',[1 0 1],'fontsize',6*fontsize_factor,'horizontalalignment','center','verticalalignment','top');
    for ep=find(epochs_in_window')
        ptch=rectangle('Position',[relative_epochs(ep,1)+state_shift y_limPSTH(1) diff(relative_epochs(ep,:)) diff(y_limPSTH)*0.1],'EdgeColor','none','FaceColor',[0.8 0.8 0.8]); %frame for indicating FR window
        text(double(state_shift+relative_epochs(ep,1)+diff(relative_epochs(ep,:))/2),y_limPSTH(1)+diff(y_limPSTH)*0.05,epoch_names(ep),'color',[ep ep ep]/255,'verticalalignment','middle','horizontalalignment','center','fontsize',6*fontsize_factor);
        uistack(ptch, 'bottom');
    end
    line([state_seperator state_shift+t_after_state],[y_limframe(1) y_limframe(1)],'linewidth',1,'color','k');
    state_seperator=state_shift + t_after_state + 0.1;
    state_seperator_max=max([state_seperator_max;state_seperator]);
end
line([0 0],y_limframe,'linewidth',1,'color','k');
set(gca,'ylim',y_lim,'xlim',[0 state_seperator_max-0.1]);
y_tick      =get(gca,'ytick');
y_tick      =y_tick(y_tick>=y_limPSTH(1));
for k=1:numel(y_tick)
    text(-(state_seperator_max-0.1)/400, y_tick(k),num2str(y_tick(k)),'fontsize',8*fontsize_factor,'horizontalalignment','right');
    line([0 (state_seperator_max-0.1)/400],[y_tick(k) y_tick(k)],'color','k');
end
remove_axis('xy')
end