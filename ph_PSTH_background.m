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

%% this part here is quite annoying, but speeds up calculation of event lines later:
% trial_onsets will be a [number of trials] x [number of states] matrix,
unique_states=unique([trial.states]);
trial_onsets=NaN(numel(trial),numel(unique_states));
for t=1:numel(trial)
    [~,si]=ismember(unique_states,trial(t).states);
    [sisorted,ssi]=sort(si);
    ssi=ssi(sisorted~=0);
    trial_onsets(t,ssi)=trial(t).states_onset;
end

for w=1:size(keys.PSTH_WINDOWS,1)
    sta             =keys.PSTH_WINDOWS{w,2};
    t_before_state  =keys.PSTH_WINDOWS{w,3};
    t_after_state   =keys.PSTH_WINDOWS{w,4};
    window_label    =keys.PSTH_WINDOWS{w,1};
    state_shift     =state_seperator-t_before_state;
    
    %% background
    [state_names,absolute_state_onsets,relative_state_onset,relative_epochs,epoch_names,states]=ph_state_onsets(trial,trial_onsets,sta,keys);
    %[state_names,absolute_state_onsets,relative_state_onset,relative_epochs,epoch_names,states]=ph_state_onsets(trial,sta,keys);
    
    states_in_window=relative_state_onset>=t_before_state & relative_state_onset<=t_after_state;
    epochs_in_window=relative_epochs(:,1)>=t_before_state & relative_epochs(:,2)<=t_after_state ;
    epochs_in_window(~ismember(epoch_names,keys.ANOVAS.main))=false;
    states_in_window(~ismember(states,keys.plot.events))=false;
    % plot only desired events and epochs
    state_positions=relative_state_onset(states_in_window)+state_shift;
    [state_positions, state_positions_idx]=sort(state_positions);
    state_names_in_window=state_names(states_in_window);
    state_names_in_window=state_names_in_window(state_positions_idx);
    x_limframe=get(gca,'xlim');
    if ~isempty(state_positions)
        %         if w == 2
        %             line([state_positions([1 4]); state_positions([1 4])],y_limframe,'color',[0.8 0.8 0.8],'linewidth',1,'LineStyle', '--');
        %         end
        %here we loop, just so not everything is one text object later, which is sort of a bug in export_fig
        temp_state_onsets=absolute_state_onsets(states_in_window);
        temp_state_onsets=temp_state_onsets(state_positions_idx);
        for p=1:numel(state_positions)
            if keys.plot.rotate_time_labels
                text(double(state_positions(p)),diff(y_limframe)*-0.02+y_limframe(1),num2str(round(temp_state_onsets(p)*1000)),...
                    'rotation',45,'fontsize',8*fontsize_factor,'horizontalalignment','center','verticalalignment','top'); %this is a trick: becuase they are plotted in different colors, they will be different elements in adobe
            else
                text(double(state_positions(p)),diff(y_limframe)*-0.02+y_limframe(1),num2str(round(temp_state_onsets(p)*1000)),...
                    'fontsize',8*fontsize_factor,'horizontalalignment','center','verticalalignment','top'); %this is a trick: becuase they are plotted in different colors, they will be different elements in adobe
            end
            if p>1 && state_positions(p)-state_positions(p-1)<diff(x_limframe)/100
                current_y=prev_y-diff(y_limframe)*0.15;
            else
                current_y=y_limframe(2)-diff(y_limframe)*0.15;
            end
            text(double(state_positions(p)),current_y,state_names_in_window(p),'fontsize',8*fontsize_factor,'interpreter','none','rotation',90,'verticalalignment','top');  %label for state
            prev_y=current_y;
        end
        line([unique(state_positions); unique(state_positions)],y_limframe,'color',[0.8 0.8 0.8],'linewidth',1,'LineStyle', '--');
        line([state_shift state_shift],y_limframe,'color',[0 0 0],'linewidth',1.5,'LineStyle', '--');
        
        %         text(double(state_positions'),ones(size(state_positions'))*(y_limframe(2)-diff(y_limframe)*0.15),state_names(states_in_window),...
        %             'fontsize',6*fontsize_factor,'interpreter','none','rotation',90,'verticalalignment','top');  %label for state
        
        %% MP: WHAT IS THIS?
        %        text(state_shift-0.075,diff(y_limframe)*-0.07+y_limframe(1),num2str(state_shift),'fontsize',10*fontsize_factor)
    end
    
    %% this here is LS version
    text(state_seperator+(t_after_state-t_before_state)/2,y_limPSTH(2)-diff(y_limPSTH)*0.05,window_label,...
        'color',[0 0 0],'fontsize',11*fontsize_factor,'horizontalalignment','center','verticalalignment','top');
    %     %% MP
    %     if w == 1
    %         text(state_shift+0.1,y_limPSTH(2)-diff(y_limPSTH)*0.05,window_label,...
    %             'color',[0 0 0],'fontsize',11*fontsize_factor,'horizontalalignment','center','verticalalignment','top');
    %     else
    %         text(state_shift-0.15,y_limPSTH(2)-diff(y_limPSTH)*0.05,window_label,...
    %             'color',[0 0 0],'fontsize',11*fontsize_factor,'horizontalalignment','center','verticalalignment','top');
    %     end
    for ep=find(epochs_in_window')
        %         ptch=rectangle('Position',[relative_epochs(ep,1)+state_shift y_limPSTH(1) diff(relative_epochs(ep,:)) diff(y_limPSTH)*0.1],'EdgeColor','none','FaceColor',[0.8 0.8 0.8]); %frame for indicating FR window
        %         text(double(state_shift+relative_epochs(ep,1)+diff(relative_epochs(ep,:))/2),y_limPSTH(1)+diff(y_limPSTH)*0.01,epoch_names(ep),'color',[ep ep ep]/255,'verticalalignment','middle','horizontalalignment','left','rotation',90,'fontsize',6*fontsize_factor);
        %% try MP here: 45° epoch labels
        ptch=rectangle('Position',[relative_epochs(ep,1)+state_shift y_limPSTH(1) diff(relative_epochs(ep,:)) diff(y_limPSTH)],'EdgeColor','none','FaceColor',[0.8 0.8 0.8]); %frame for indicating FR window
        if keys.plot.rotate_epoch_labels
            text(double(state_shift+relative_epochs(ep,1)+diff(relative_epochs(ep,:))/2),y_limPSTH(1)+diff(y_limPSTH)*0.02,epoch_names(ep),'color',[0 0 0]/255,'verticalalignment','middle','horizontalalignment','center','fontsize',8*fontsize_factor,'rotation',45);
        else
            text(double(state_shift+relative_epochs(ep,1)+diff(relative_epochs(ep,:))/2),y_limPSTH(1)+diff(y_limPSTH)*0.02,epoch_names(ep),'color',[0 0 0]/255,'verticalalignment','middle','horizontalalignment','center','fontsize',8*fontsize_factor);
            
        end
        
        uistack(ptch, 'bottom');
    end
    line([state_seperator state_shift+t_after_state],[y_limframe(1) y_limframe(1)],'linewidth',1,'color','k');
    state_seperator=state_shift + t_after_state + 0.1;
    state_seperator_max=max([state_seperator_max;state_seperator]);
end
%% time scaling
w1end=keys.PSTH_WINDOWS{1,4}-keys.PSTH_WINDOWS{1,3};
line([w1end-0.05 w1end+0.15],[y_limframe(1)+y_limframe(2)*0.05 y_limframe(1)+y_limframe(2)*0.05],'linewidth',1,'color','k')
text(w1end+0.05,(y_limframe(1)+y_limframe(2)*0.05),'200 ms','horizontalalignment','center','verticalalignment','bottom','fontsize',4*fontsize_factor)
if keys.plot.rotate_time_labels
    text(0,diff(y_limframe)*-0.02+y_limframe(1),'0','fontsize',8*fontsize_factor,'rotation',45)
else
    text(0,diff(y_limframe)*-0.02+y_limframe(1),'0','fontsize',8*fontsize_factor)
end
line([0 0],y_limframe,'linewidth',1,'color','k');
set(gca,'ylim',single(y_lim),'xlim',[0 state_seperator_max-0.1]);
y_tick      =get(gca,'ytick');
y_tick      =y_tick(y_tick>=y_limPSTH(1));
y_tick(abs(y_tick)<5*eps)=0;
for k=1:numel(y_tick)
    text(-(state_seperator_max-0.1)/400, y_tick(k),num2str(y_tick(k)),'fontsize',8*fontsize_factor,'horizontalalignment','right');
    line([0 (state_seperator_max-0.1)/400],[y_tick(k) y_tick(k)],'color','k');
end
% if strcmp(keys.PO.normalization,'none') & ~keys.PO.FR_subtract_baseline
%     text(-0.2, (round(numel(y_tick)/2)-1) ,'Firing rate (spk/s)','fontsize',8*fontsize_factor,'rotation',90)
% elseif keys.PO.FR_subtract_baseline
%  text(-0.2, (round(numel(y_tick)/2)-1) ,'Normalized firing rate (spk/s)','fontsize',8*fontsize_factor,'rotation',90)
% else
% text(-0.2, (round(numel(y_tick)/2)-1) ,'Normalized firing rate','fontsize',8*fontsize_factor,'rotation',90)
% end
%  text(state_shift-1,diff(y_limframe)*-0.09+y_limframe(1),'Time(s)','fontsize',8*fontsize_factor)
remove_axis('xy')
end