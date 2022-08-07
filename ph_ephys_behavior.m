function ph_ephys_behavior(modified_keys) %project,version,dataset,epoch,basefolder)
% project='Pulv_eye_gaze_position';
% version='paper';
% dataset='_dPulv_PT0_Msac_mov';


for fn=fieldnames(modified_keys)'
    keys.(fn{:})=modified_keys.(fn{:});
end
recalibrated=1;
corrective=1;
epoch={'Thol',5,	0.2,    0.5};
databasefolder='C:\Users\lschneider\Desktop';


target=keys.tt.selection{1,2};

project=keys.project;
version=keys.project_version;
dataset=['_' target '_PT' num2str(keys.tt.perturbations) '_' keys.conditions_to_plot{:} '_' keys.arrangement(1:3)]; %'_dPulv_PT0_Msac_opt';
datapath=[databasefolder filesep project dataset];
load(['Y:\Projects\' project '\ephys\' version '\behaviour_filelist.mat']);

%keys.colors.fix_offset          =[236 32 38; 16 159 218; 247 148 36]/255;




if corrective==1
    prefix='corr ';
    %  keys.path_to_save=['Y:\Projects\' project '\behavior\correctiverecalibrated'];
    keys.cal.MA_selection                   ={'correct_offset',0,'success',1,'display',0,'summary',0,'keep_raw_data',1,'saccade_definition',4,'reach_1st_pos',1,'correlation_conditions',{},'nsacc_max', 20, 'sac_ini_t', 20,  'sac_end_t',  10, 'sac_min_dur', 0, 'sac_min_dur', 0, 'which_saccade', 'corrective'};                        % if you want to run MA with specific settings

else
    prefix='';
    %  keys.path_to_save=['Y:\Projects\' project '\behavior\correctivenotrecalibrated'];
    keys.cal.MA_selection                   ={'correct_offset',0,'success',1,'display',0,'keep_raw_data',1,'saccade_definition',4,'reach_1st_pos',1,'correlation_conditions',{}};                        % if you want to run MA with specific settings
    
end
if recalibrated==1
    prefix=[prefix 'recal '];
    % keys.path_to_save=['Y:\Projects\' project '\behavior\correctiverecalibrated'];
else
    prefix=[prefix 'rawMP '];
    % keys.path_to_save=['Y:\Projects\' project '\behavior\correctivenotrecalibrated'];
end
keys.plot.export=1;


%% modify filelist to create batches for each session
for monkeys=[{'Both'} keys.monkeys]
    monkey=monkeys{:};
    if strcmp(monkey,'Both')
        monkey='Lin';
        taskfield=[monkey dataset];
        mainfolder=filelist.(taskfield){1};
        dashstr=strfind(mainfolder,'\');
        mainfolder=mainfolder(1:dashstr(end));
        filelist_in=filelist_formatted.(taskfield);
        filelist_in(:,1)=cellstr(num2str(vertcat(filelist_in{:,1})));
        filelist_in(:,1)=strcat(mainfolder, filelist_in(:,1));
        
        if recalibrated==1
            aa=dir([datapath filesep monkey filesep]);
            sessions={aa([aa.isdir]).name};
            sessions=sessions(3:end);
            for s=1:numel(sessions)
                filelist_in{s,1}=[datapath filesep monkey filesep sessions{s}];
            end
        end
        
        monkey='Cur';
        taskfield=[monkey dataset];
        mainfolder=filelist.(taskfield){1};
        dashstr=strfind(mainfolder,'\');
        mainfolder=mainfolder(1:dashstr(end));
        filelist_in2=filelist_formatted.(taskfield);
        filelist_in2(:,1)=cellstr(num2str(vertcat(filelist_in2{:,1})));
        filelist_in2(:,1)=strcat(mainfolder, filelist_in2(:,1));
        
        
        if recalibrated==1
            aa=dir([datapath filesep monkey filesep]);
            sessions={aa([aa.isdir]).name};
            sessions=sessions(3:end);
            for s=1:numel(sessions)
                filelist_in2{s,1}=[datapath filesep monkey filesep sessions{s}];
            end
        end
        
        filelist_in=vertcat(filelist_in, filelist_in2);
        monkey='Both';
    else
        taskfield=[monkey dataset];
        mainfolder=filelist.(taskfield){1};
        dashstr=strfind(mainfolder,'\');
        mainfolder=mainfolder(1:dashstr(end));
        filelist_in=filelist_formatted.(taskfield);
        filelist_in(:,1)=cellstr(num2str(vertcat(filelist_in{:,1})));
        filelist_in(:,1)=strcat(mainfolder, filelist_in(:,1));
        
        
        if recalibrated==1
            aa=dir([datapath filesep monkey filesep]);
            sessions={aa([aa.isdir]).name};
            sessions=sessions(3:end);
            for s=1:numel(sessions)
                filelist_in{s,1}=[datapath filesep monkey filesep sessions{s}];
            end
        end
    end
    %end
    MA_input={};
    for c=1:size(filelist_in,1)
        MA_input=[MA_input, {{filelist_in{c,1}, filelist_in{c,2}}},{keys.cal.MA_selection}];
    end
    %if ~exist(['Y:\Projects\' project '\ephys\' version '\behaviour_' dataset '.mat'],'file')
    data=monkeypsych_analyze_working(MA_input{:});
    %save(['Y:\Projects\' project '\ephys\' version '\behaviour_' dataset '.mat'],'data');
    %     else
    %         load(['Y:\Projects\' project '\ephys\' version '\behaviour_' dataset '.mat']);
    %     end
    
    condition_parameters  ={'reach_hand','choice','perturbation'};
    
    %     per_trial.perturbation=[all_trialz.perturbation];
    %     per_trial.perturbation(ismember(per_trial.perturbation, keys.cal.perturbation_groups{1}))=0;
    %     per_trial.perturbation(ismember(per_trial.perturbation, keys.cal.perturbation_groups{2}))=1;
    %     u_perturbation    =unique(per_trial.perturbation);
    %     u_perturbation=u_perturbation(~isnan(u_perturbation));
    
    u_hemifields=[-1,0,1]; % why does this have to be hardcoded? ---> Because case not defined yet, case defines positions !!
    
    
    
    all_data=[data{:}];
    all_task=vertcat(all_data.task);
    all_binary=vertcat(all_data.binary);
    %all_reaches=vertcat(all_data.reaches);
    u_types     =unique([all_task.type]);
    u_effectors =unique([all_task.effector]);
    all_hands       =[all_task.reach_hand];
    all_hands(isnan(all_hands))=0;
    all_choices       =[all_binary.choice];
    u_hands     =unique(all_hands);
    u_choice    =unique(all_choices);
    u_hands     =u_hands(ismember(u_hands,keys.tt.hands));
    u_choice    =u_choice(ismember(u_choice,keys.tt.choices));
    
    all_type_effectors      = combvec(u_types,u_effectors)';
    
    condition_matrix            = combvec(u_hands,u_choice)';%, u_perturbation,u_hemifields)';
    %conditions_out              = combvec(u_effectors,u_hands,u_choice, u_perturbation)';
    % conditions_hf               = combvec(u_hemifields,conditions_out')';
    % conditions_hf_complete      = combvec(u_hemifields,conditions_out')';
    % conditions_pref             = combvec([1 0],conditions_out')';
    type_effectors =[];
    
    for t=1:size(all_type_effectors,1)
        typ=all_type_effectors(t,1);
        eff=all_type_effectors(t,2);
        [~, type_effector_short{t}]=MPA_get_type_effector_name(typ,eff);
        if ~ismember(type_effector_short{t},keys.conditions_to_plot) %|| sum(tr_con)<1
            type_effector_short{t}=[];
            continue;
        end
        type_effectors=[type_effectors; all_type_effectors(t,:)];
    end
    type_effector_short(cellfun(@isempty,type_effector_short))=[];
    
    
    for te=1:size(type_effectors,1)
        type=type_effectors(te,1);
        effector=type_effectors(te,2);
        current_effector='saccades';
        all_saccades=vertcat(all_data.(current_effector));
        positions=unique(round([all_saccades.tar_pos] - [all_saccades.fix_pos]));
        positions(isnan(positions))=[];
        fixations=unique(real([all_saccades.fix_pos])); % will need to use precision here
        distances=abs(positions);
        unique_distances=unique(distances);
        for c= 1:size(condition_matrix,1) %% not sure if dom here is 1
            hand=condition_matrix(c,1);
            choice=condition_matrix(c,2);
            condition_title=['hn' num2str(hand) 'ch' num2str(choice)];
            
            clear per_pos per_hemi per_session
            for s=1:numel(data)
                saccades=[data{s}.saccades];
                states=[data{s}.states];
                raw=[data{s}.raw];
                task=[data{s}.task];
                reach_hand=[task.reach_hand];reach_hand(isnan(reach_hand))=0;
                binary=[data{s}.binary];
                te_idx=[task.type]==type & [task.effector]==effector & [binary.choice]==choice & reach_hand==hand;
                saccades=saccades(te_idx);
                states=states(te_idx);
                raw=raw(te_idx);
                per_session.RTs(s)=nanmean([saccades.lat])*1000;
                per_session.dur(s)=nanmean([saccades.dur])*1000;
                per_session.amp(s)=nanmean(abs([saccades.endpos]-[saccades.startpos]));
                per_session.vel(s)=nanmean([saccades.velocity]);
                per_session.pre(s)=nanmean([saccades.endpos]-[saccades.fix_pos]);
                per_session.percent_sac(s)=sum(~isnan([saccades.lat]))/sum(te_idx)*100;
                per_session.Ntr(s)=sum(te_idx);
                per_session.Ntr_l(s)=sum(real([saccades.tar_pos]-[saccades.fix_pos])<0);
                
                fix_pos=[saccades.fix_pos];
                for h=1:2
                    lr=(h-1.5)*2;
                    tr_idx=sign(real([saccades.tar_pos] - fix_pos))==lr;
                    
                    
                    tr_idx_int=find(tr_idx);
                    if isempty(tr_idx_int)
                        continue;
                    end
                    av_pos=[];
                    acq_dur=[];
                    for t=1:numel(tr_idx_int)
                        idx=tr_idx_int(t);
                        on=states(idx).MP_states_onset([states(idx).MP_states]==epoch{2});
                        time=raw(idx).time_axis >= on+epoch{3}  & raw(idx).time_axis <= on+epoch{4}; %>=?
                        av_pos(t)=median(raw(idx).x_eye(time))+1i*median(raw(idx).y_eye(time));
                        acq_dur(t)=1;
                        %  acq_dur(t)=states(idx).MP_states_onset([states(idx).MP_states]==10)-...
                        %  states(idx).MP_states_onset([states(idx).MP_states]==9);
                    end
                    
                    per_hemi(h).sac_end_all{s}=[saccades(tr_idx).endpos]-fix_pos(tr_idx)+real(fix_pos(tr_idx));
                    per_hemi(h).eye_in_thol_all{s}=av_pos-fix_pos(tr_idx)+real(fix_pos(tr_idx));
                    
                    per_hemi(h).sac_off(s)=nanmean([saccades(tr_idx).accuracy_xy]);
                    per_hemi(h).RTs(s)=nanmean([saccades(tr_idx).lat])*1000;
                    per_hemi(h).dur(s)=nanmean([saccades(tr_idx).dur])*1000;
                    per_hemi(h).vel(s)=nanmean([saccades(tr_idx).velocity]);
                    per_hemi(h).acq_dur(s)=nanmean(acq_dur-[saccades(tr_idx).lat])*1000;
                    per_hemi(h).tar_rad(s)=nanmean([saccades(tr_idx).tar_rad]);
                    per_hemi(h).sac_end(s)=nanmean([saccades(tr_idx).endpos]-fix_pos(tr_idx)+real(fix_pos(tr_idx)));
                    per_hemi(h).sac_prec(s)=nanstd(real([saccades(tr_idx).endpos]-fix_pos(tr_idx)));
                    per_hemi(h).eye_in_thol(s)=nanmean(av_pos-fix_pos(tr_idx)+real(fix_pos(tr_idx)));
                    per_hemi(h).amp(s)=nanmean(abs([saccades(tr_idx).endpos]-[saccades(tr_idx).startpos]));
                    per_hemi(h).finalamp(s)=nanmean(abs(av_pos-fix_pos(tr_idx)));
                    per_hemi(h).final_prec(s)=nanstd(av_pos-fix_pos(tr_idx));
                    
                    per_hemi(h).meanRTs= nanmean([per_hemi(h).RTs]);
                    per_hemi(h).meandur= nanmean([per_hemi(h).dur]);
                    per_hemi(h).meanvel= nanmean([per_hemi(h).vel]);
                    per_hemi(h).meanamp= nanmean([per_hemi(h).amp]);
                    per_hemi(h).meanacq_dur= nanmean([per_hemi(h).acq_dur]);
                    per_hemi(h).meanfinalamp= nanmean([per_hemi(h).finalamp]);
                    %
                    %                         per_hemi(h).fixation(s)=f;
                    %                         per_hemi(h).distance(s)=abs(positions(p));
                    %                         per_hemi(h).position(s)=p;
                    per_hemi(h).percent_sac(s)=sum(~isnan([saccades(tr_idx).lat]))/sum(tr_idx)*100;
                    per_hemi(h).Ntr(s)=sum(tr_idx);
                    per_hemi(h).percent_trials(s)=sum(tr_idx)/per_session.Ntr(s)*100;
                end
                
                for f=1:numel(fixations)
                    for p=1:numel(positions)
                        tr_idx=(round([saccades.tar_pos] - fix_pos))==positions(p) & real([saccades.fix_pos])==fixations(f);
                        
                        
                        tr_idx_int=find(tr_idx);
                        
                        %preallocate...
                        
                        %                         per_pos(f,p).sac_end_all{s}=[saccades(tr_idx).endpos]-fix_pos(tr_idx)+fixations(f);
                        %                         per_pos(f,p).eye_in_thol_all{s}=av_pos-fix_pos(tr_idx)+fixations(f);
                        %
                        per_pos(f,p).sac_off(s)=NaN;
                        per_pos(f,p).RTs(s)=NaN;
                        per_pos(f,p).dur(s)=NaN;
                        per_pos(f,p).vel(s)=NaN;
                        per_pos(f,p).acq_dur(s)=NaN;
                        per_pos(f,p).tar_rad(s)=NaN;
                        per_pos(f,p).sac_end(s)=NaN;
                        per_pos(f,p).sac_prec(s)=NaN;
                        per_pos(f,p).eye_in_thol(s)=NaN;
                        per_pos(f,p).amp(s)=NaN;
                        per_pos(f,p).finalamp(s)=NaN;
                        per_pos(f,p).final_prec(s)=NaN;
                        
                        per_pos(f,p).fixation(s)=NaN;
                        per_pos(f,p).distance(s)=NaN;
                        per_pos(f,p).position(s)=NaN;
                        per_pos(f,p).per_session.Ntr(s)=NaN;
                        per_pos(f,p).percent_trials(s)=NaN;
                        
                        if isempty(tr_idx_int)
                            continue;
                        end
                        av_pos=[];
                        acq_dur=[];
                        for t=1:numel(tr_idx_int)
                            idx=tr_idx_int(t);
                            on=states(idx).MP_states_onset([states(idx).MP_states]==epoch{2});
                            time=raw(idx).time_axis >= on+epoch{3}  & raw(idx).time_axis <= on+epoch{4}; %>=?
                            av_pos(t)=median(raw(idx).x_eye(time))+1i*median(raw(idx).y_eye(time));
                            acq_dur(t)=1;
                            %  acq_dur(t)=states(idx).MP_states_onset([states(idx).MP_states]==10)-...
                            %  states(idx).MP_states_onset([states(idx).MP_states]==9);
                        end
                        
                        per_pos(f,p).sac_end_all{s}=[saccades(tr_idx).endpos]-fix_pos(tr_idx)+fixations(f);
                        per_pos(f,p).eye_in_thol_all{s}=av_pos-fix_pos(tr_idx)+fixations(f);
                        
                        per_pos(f,p).sac_off(s)=nanmean([saccades(tr_idx).accuracy_xy]);
                        per_pos(f,p).RTs(s)=nanmean([saccades(tr_idx).lat])*1000;
                        per_pos(f,p).dur(s)=nanmean([saccades(tr_idx).dur])*1000;
                        per_pos(f,p).vel(s)=nanmean([saccades(tr_idx).velocity]);
                        per_pos(f,p).acq_dur(s)=nanmean(acq_dur-[saccades(tr_idx).lat])*1000;
                        per_pos(f,p).tar_rad(s)=nanmean([saccades(tr_idx).tar_rad]);
                        per_pos(f,p).sac_end(s)=nanmean([saccades(tr_idx).endpos]-fix_pos(tr_idx)+fixations(f));
                        per_pos(f,p).sac_prec(s)=nanstd(real([saccades(tr_idx).endpos]-fix_pos(tr_idx)));
                        per_pos(f,p).eye_in_thol(s)=nanmean(av_pos-fix_pos(tr_idx)+fixations(f));
                        per_pos(f,p).amp(s)=nanmean(abs([saccades(tr_idx).endpos]-[saccades(tr_idx).startpos]));
                        per_pos(f,p).finalamp(s)=nanmean(abs(av_pos-fix_pos(tr_idx)));
                        per_pos(f,p).final_prec(s)=nanstd(av_pos-fix_pos(tr_idx));
                        
                        per_pos(f,p).meanRTs= nanmean([per_pos(f,p).RTs]);
                        per_pos(f,p).meandur= nanmean([per_pos(f,p).dur]);
                        per_pos(f,p).meanvel= nanmean([per_pos(f,p).vel]);
                        per_pos(f,p).meanamp= nanmean([per_pos(f,p).amp]);
                        per_pos(f,p).meanacq_dur= nanmean([per_pos(f,p).acq_dur]);
                        per_pos(f,p).meanfinalamp= nanmean([per_pos(f,p).finalamp]);
                        
                        per_pos(f,p).fixation(s)=f;
                        per_pos(f,p).distance(s)=abs(positions(p));
                        per_pos(f,p).position(s)=p;
                        per_pos(f,p).per_session.Ntr(s)=sum(tr_idx);
                        per_pos(f,p).percent_trials(s)=sum(tr_idx)/per_session.Ntr(s)*100;
                        per_pos(f,p).percent_sac(s)=sum(~isnan([saccades(tr_idx).lat]))/sum(tr_idx)*100;
                    end
                end
            end
            
            %% plotting
            
            
            %% trials per position
            plot_title=[prefix ' ' monkey ' ' type_effector_short{te} ' ' condition_title ' percent saccades detected'];
            filename=plot_title;
            figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            bins=0:10:100;
            histo=hist(per_session.percent_sac,bins);
            bar(0:10:100,histo);
            [~, p_lr]=ttest(per_hemi(1).percent_sac,per_hemi(2).percent_sac);
            
            title([num2str(round(nanmean(per_session.percent_sac))) ' + ' num2str(round(sterr(per_session.percent_sac))) ' total; ' ...
                   num2str(round(nanmean(per_hemi(1).percent_sac))) ' + ' num2str(round(sterr(per_hemi(1).percent_sac))) ' / '...
                   num2str(round(nanmean(per_hemi(2).percent_sac))) ' + ' num2str(round(sterr(per_hemi(2).percent_sac)))' ' L/R, p=' num2str(round(p_lr*1000)/1000)]);
            ph_title_and_save(figure_handle,filename,plot_title,keys)
            
            
            %% trials per position
            plot_title=[prefix ' ' monkey ' ' type_effector_short{te} ' ' condition_title ' percent trials per position'];
            filename=plot_title;
            figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            Colors=cool(round(max(max([per_pos.percent_trials])))+1);
            colormap(Colors);
            angles=0:pi/180:2*pi;
            
            hold on
            for f=1:numel(fixations)
                
                for p=1:numel(positions)
                    
                    r=nanmean([per_pos(f,p).percent_trials])/5;
                    rs=r-sterr([per_pos(f,p).percent_trials])/5;
                    rl=r+sterr([per_pos(f,p).percent_trials])/5;
                    x0=positions(p)+fixations(f);
                    %scatter(real(x0),imag(x0),'+','markeredgecolor','k');
                    patch(real(x0)+r*sin(angles), imag(x0)+r*cos(angles),Colors(round(nanmean([per_pos(f,p).percent_trials]))+1,:));
                    plot(real(x0)+rs*sin(angles), imag(x0)+rs*cos(angles),'color','k');%Colors(round(nanmean([per_pos(f,p).percent_trials]))+1,:));
                    plot(real(x0)+rl*sin(angles), imag(x0)+rl*cos(angles),'color','k');%Colors(round(nanmean([per_pos(f,p).percent_trials]))+1,:));
                    text(real(x0),imag(x0),num2str(nanmean([per_pos(f,p).percent_trials])),'HorizontalAlignment','center');
                end
            end
            [~, p_lr]=ttest(per_hemi(1).percent_trials,per_hemi(2).percent_trials);
            title([num2str(round(nanmean(per_hemi(1).percent_trials))) ' + ' num2str(round(sterr(per_hemi(1).percent_trials))) ' / '...
                   num2str(round(nanmean(per_hemi(2).percent_trials))) ' + ' num2str(round(sterr(per_hemi(2).percent_trials)))' ' L/R, p=' num2str(round(p_lr*1000)/1000)]);
            axis equal
            colorbar;
            ph_title_and_save(figure_handle,filename,plot_title,keys)
            
            %% accuracy
            plot_title=[prefix ' ' monkey ' ' type_effector_short{te} ' ' condition_title ' accuracy'];
            filename=plot_title;
            figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            Colors=keys.colors.fix_offset;
            angles=0:pi/180:2*pi;
            
            hold on
            for f=1:numel(fixations)
                for p=1:numel(positions)
                    
                    r=5;
                    x0=positions(p)+fixations(f);
                    scatter(real(x0),imag(x0),'+','markeredgecolor','k');
                    plot(real(x0)+r*sin(angles), imag(x0)+r*cos(angles),'k');
                    plot([(nanmean([per_pos(f,p).sac_end])) , nanmean([per_pos(f,p).eye_in_thol])],'color',Colors(f,:));
                    scatter(real(nanmean([per_pos(f,p).sac_end])),imag(nanmean([per_pos(f,p).sac_end])),30,Colors(f,:),'filled');
                end
            end
            axis equal
            ph_title_and_save(figure_handle,filename,plot_title,keys)
            
            
            %% accuracy per session
            plot_title=[prefix ' ' monkey ' ' type_effector_short{te} ' ' condition_title ' accuracy per session'];
            filename=plot_title;
            figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            Colors=keys.colors.fix_offset;
            angles=0:pi/180:2*pi;
            
            hold on
            for f=1:numel(fixations)
                for p=1:numel(positions)
                    r=5;
                    x0=positions(p)+fixations(f);
                    scatter(real(x0),imag(x0),'+','markeredgecolor','k');
                    plot(real(x0)+r*sin(angles), imag(x0)+r*cos(angles),'k');
                    
                    %plot([(nanmean(per_pos(f,p).sac_end)) , nanmean(per_pos(f,p).eye_in_thol)],'color',Colors(f,:));
                    scatter(real([per_pos(f,p).sac_end]),imag([per_pos(f,p).sac_end]),30,Colors(f,:),'filled');
                    nanmean(per_pos(f,p).dur);
                end
            end
            axis equal
            ph_title_and_save(figure_handle,filename,plot_title,keys)
            
            %% accuracy per session
            plot_title=[prefix ' ' monkey ' ' type_effector_short{te} ' ' condition_title ' final accuracy per session'];
            filename=plot_title;
            figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            Colors=keys.colors.fix_offset;
            angles=0:pi/180:2*pi;
            
            hold on
            for f=1:numel(fixations)
                for p=1:numel(positions)
                    r=5;
                    x0=positions(p)+fixations(f);
                    scatter(real(x0),imag(x0),'+','markeredgecolor','k');
                    plot(real(x0)+r*sin(angles), imag(x0)+r*cos(angles),'k');
                    
                    %plot([(nanmean(per_pos(f,p).sac_end)) , nanmean(per_pos(f,p).eye_in_thol)],'color',Colors(f,:));
                    plot([[per_pos(f,p).sac_end]' , [per_pos(f,p).eye_in_thol]']','color',Colors(f,:));
                    scatter(real([per_pos(f,p).eye_in_thol]),imag([per_pos(f,p).eye_in_thol]),15,Colors(f,:),'filled');
                    %nanmean(per_pos(f,p).dur);
                end
            end
            axis equal
            ph_title_and_save(figure_handle,filename,plot_title,keys)
            
            
            %% accuracy per session
            Colors=keys.colors.fix_offset;
            angles=0:pi/180:2*pi;
            for s=1:numel(data)
                plot_title=[prefix ' ' monkey ' ' type_effector_short{te} ' ' condition_title ' final accuracy, session ' num2str(s)];
                filename=plot_title;
                figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
                hold on
                for f=1:numel(fixations)
                    for p=1:numel(positions)
                        if s <=numel(per_pos(f,p).sac_end_all)
                            r=5;
                            x0=positions(p)+fixations(f);
                            scatter(real(x0),imag(x0),'+','markeredgecolor','k');
                            plot(real(x0)+r*sin(angles), imag(x0)+r*cos(angles),'k');
                            
                            %plot([(nanmean(per_pos(f,p).sac_end)) , nanmean(per_pos(f,p).eye_in_thol)],'color',Colors(f,:));
                            plot([[per_pos(f,p).sac_end_all{s}]' , [per_pos(f,p).eye_in_thol_all{s}]']','color',Colors(f,:));
                            scatter(real([per_pos(f,p).eye_in_thol_all{s}]),imag([per_pos(f,p).eye_in_thol_all{s}]),15,Colors(f,:),'filled');
                            %nanmean(per_pos(f,p).dur);
                        end
                    end
                end
                axis equal
                ph_title_and_save(figure_handle,filename,plot_title,keys)
            end
            
            %% precision per session
            plot_title=[prefix ' ' monkey ' ' type_effector_short{te} ' ' condition_title ' precision per session'];
            filename=plot_title;
            figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            Colors=keys.colors.fix_offset;
            
            angles=0:pi/180:2*pi;
            
            hold on
            for f=1:numel(fixations)
                for p=1:numel(positions)
                    r=5;
                    x0=positions(p)+fixations(f);
                    scatter(real(x0),imag(x0),'+','markeredgecolor','k');
                    plot(real(x0)+r*sin(angles), imag(x0)+r*cos(angles),'k');
                    
                    x0=nanmean([per_pos(f,p).sac_end]);
                    r=nanmean([per_pos(f,p).sac_prec]);
                    plot(real(x0)+r*sin(angles), imag(x0)+r*cos(angles),'color',Colors(f,:));
                    drawring([real(x0) imag(x0)],r-sterr([per_pos(f,p).sac_prec]), r+sterr([per_pos(f,p).sac_prec]),angles,Colors(f,:));
                    %scatter(real([per_pos(f,p).sac_prec]),imag([per_pos(f,p).sac_prec]),30,Colors(f,:),'filled');
                    %nanmean(per_pos(f,p).dur);
                end
            end
            axis equal
            ph_title_and_save(figure_handle,filename,plot_title,keys)
            
            
            %% precision per session
            plot_title=[prefix ' ' monkey ' ' type_effector_short{te} ' ' condition_title ' precision of session means'];
            filename=plot_title;
            figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            Colors=keys.colors.fix_offset;
            angles=0:pi/180:2*pi;
            
            hold on
            for f=1:numel(fixations)
                for p=1:numel(positions)
                    r=5;
                    x0=positions(p)+fixations(f);
                    scatter(real(x0),imag(x0),'+','markeredgecolor','k');
                    plot(real(x0)+r*sin(angles), imag(x0)+r*cos(angles),'k');
                    
                    x0=nanmean([per_pos(f,p).sac_end]);
                    r=nanstd([per_pos(f,p).sac_end]);
                    plot(real(x0)+r*sin(angles) , imag(x0)+r*cos(angles),'color',Colors(f,:));
                    %scatter(real([per_pos(f,p).sac_prec]),imag([per_pos(f,p).sac_prec]),30,Colors(f,:),'filled');
                    %nanmean(per_pos(f,p).dur);
                end
            end
            axis equal
            ph_title_and_save(figure_handle,filename,plot_title,keys)
            
            %% final precision per session
            plot_title=[prefix ' ' monkey ' ' type_effector_short{te} ' ' condition_title ' final precision per session'];
            filename=plot_title;
            figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            Colors=keys.colors.fix_offset;
            angles=0:pi/180:2*pi;
            hold on
            for f=1:numel(fixations)
                for p=1:numel(positions)
                    r=5;
                    x0=positions(p)+fixations(f);
                    scatter(real(x0),imag(x0),'+','markeredgecolor','k');
                    plot(real(x0)+r*sin(angles), imag(x0)+r*cos(angles),'k');
                    x0=nanmean([per_pos(f,p).eye_in_thol]);
                    r=nanmean([per_pos(f,p).final_prec]);
                    plot(real(x0)+r*sin(angles) , imag(x0)+r*cos(angles),'color',Colors(f,:));
                    drawring([real(x0) imag(x0)],r-sterr([per_pos(f,p).final_prec]), r+sterr([per_pos(f,p).final_prec]),angles,Colors(f,:));
                end
            end
            axis equal
            ph_title_and_save(figure_handle,filename,plot_title,keys)
            
            
            %% final precision per session
            plot_title=[prefix ' ' monkey ' ' type_effector_short{te} ' ' condition_title ' final precision of session means'];
            filename=plot_title;
            figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            Colors=keys.colors.fix_offset;
            angles=0:pi/180:2*pi;
            
            hold on
            for f=1:numel(fixations)
                for p=1:numel(positions)
                    r=5;
                    x0=positions(p)+fixations(f);
                    scatter(real(x0),imag(x0),'+','markeredgecolor','k');
                    plot(real(x0)+r*sin(angles), imag(x0)+r*cos(angles),'k');
                    x0=nanmean([per_pos(f,p).eye_in_thol]);
                    r=nanstd([per_pos(f,p).eye_in_thol]);
                    plot(real(x0)+r*sin(angles) , imag(x0)+r*cos(angles),'color',Colors(f,:));
                end
            end
            axis equal
            ph_title_and_save(figure_handle,filename,plot_title,keys)
            
            %% duration
            
            plot_title=[prefix ' ' monkey ' ' type_effector_short{te} ' ' condition_title ' duration'];
            filename=plot_title;
            figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            
            bins=min(floor([per_pos.meandur])):1:max(ceil([per_pos.meandur]));
            clear histo
            for d=1:numel(unique_distances)
                histo(:,d)=hist([per_pos(:,distances==unique_distances(d)).meandur],bins);
            end
            bar(bins,histo,'stacked');
            x_lim=get(gca,'xlim');
            y_lim=get(gca,'ylim');
            
            for f=1:numel(fixations)
                for d=1:numel(unique_distances)
                    persession=nanmean(vertcat(per_pos(f,distances==unique_distances(d)).dur),1);
                    text_to_plot=[num2str(nanmean(persession)) '+' num2str(sterr(persession))];
                    text(x_lim(1)+diff(x_lim)/5*f,y_lim(2)-diff(y_lim)/20*d,text_to_plot);
                end
            end
            
            pd = anovan([per_pos.dur]',[[per_pos.fixation]' [per_pos.distance]'],'model','full','display','off');
            pp = anovan([per_pos.dur]',[[per_pos.fixation]' [per_pos.position]'],'model','full','display','off');
            
            text_to_plot=['p anova FxD: ' num2str(pd(1)) ', ' num2str(pd(2)) ', ' num2str(pd(3)) ];
            text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*(d+1),text_to_plot);
            text_to_plot=['p anova FxP: ' num2str(pp(1)) ', ' num2str(pp(2)) ', ' num2str(pp(3)) ];
            text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*(d+2),text_to_plot);
            
            xlabel('saccade duration [ms]');
            ylabel('N conditions');
            title(['across targets: ' num2str(mean([per_pos(:).meandur])) ' + ' num2str(sterr([per_pos(:).meandur]))...
                'across sessions: ' num2str(mean(per_session.dur)) ' + ' num2str(sterr(per_session.dur)) ...
                'across targets, end as tihol: ' num2str(mean([per_pos(:).meanacq_dur])) ' + ' num2str(sterr([per_pos(:).meanacq_dur])) ]);
            
            legend(cellstr(num2str(unique_distances')));
            ph_title_and_save(figure_handle,filename,plot_title,keys)
            
            %% latencies
            
            plot_title=[prefix ' ' monkey ' ' type_effector_short{te} ' ' condition_title ' latency'];
            filename=plot_title;
            figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            
            subplot(2,1,1)
            clear histo
            bins=min(floor([per_hemi.RTs])):1:max(ceil([per_hemi.RTs]));
            for d=1:numel(per_hemi)
                histo(:,d)=hist([per_hemi(d).RTs],bins);
            end
            colormap([0 0 1;0 1 0])
            plot(bins,histo);
            [~, p_lr]=ttest(per_hemi(1).RTs,per_hemi(2).RTs);
            title([num2str(round(nanmean(per_hemi(1).RTs))) ' + ' num2str(round(sterr(per_hemi(1).RTs))) ' / '...
                   num2str(round(nanmean(per_hemi(2).RTs))) ' + ' num2str(round(sterr(per_hemi(2).RTs)))' ' L/R, p=' num2str(round(p_lr*1000)/1000)]);
            
            

            subplot(2,1,2)
            bins=min(floor([per_pos.meanRTs])):1:max(ceil([per_pos.meanRTs]));
            clear histo
            for d=1:numel(unique_distances)
                histo(:,d)=hist([per_pos(:,distances==unique_distances(d)).meanRTs],bins);
            end
            bar(bins,histo,'stacked');
            x_lim=get(gca,'xlim');
            y_lim=get(gca,'ylim');

            for f=1:numel(fixations)
                for d=1:numel(unique_distances)
                    persession=nanmean(vertcat(per_pos(f,distances==unique_distances(d)).RTs),1);
                    text_to_plot=[num2str(nanmean(persession)) '+' num2str(sterr(persession))];
                    text(x_lim(1)+diff(x_lim)/5*f,y_lim(2)-diff(y_lim)/20*d,text_to_plot);
                end
            end
            
            pd = anovan([per_pos.RTs]',[[per_pos.fixation]' [per_pos.distance]'],'model','full','display','off');
            pp = anovan([per_pos.RTs]',[[per_pos.fixation]' [per_pos.position]'],'model','full','display','off');
            

            text_to_plot=['p anova FxD: ' num2str(pd(1)) ', ' num2str(pd(2)) ', ' num2str(pd(3)) ];
            text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*(d+1),text_to_plot);
            text_to_plot=['p anova FxP: ' num2str(pp(1)) ', ' num2str(pp(2)) ', ' num2str(pp(3)) ];
            text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*(d+2),text_to_plot);
            
            xlabel('saccade reaction time [ms]');
            ylabel('N conditions');
            title(['across targets: ' num2str(mean([per_pos(:).meanRTs])) ' + ' num2str(sterr([per_pos(:).meanRTs]))...
                'across sessions: ' num2str(mean(per_session.RTs)) ' + ' num2str(sterr(per_session.RTs))]);
            
            legend(cellstr(num2str(unique_distances')));
            ph_title_and_save(figure_handle,filename,plot_title,keys)
            
            %% amplitude
            
            plot_title=[prefix ' ' monkey ' ' type_effector_short{te} ' ' condition_title ' amplitude'];
            filename=plot_title;
            figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            
            bins=min(floor([per_pos.meanamp])):1:max(ceil([per_pos.meanamp]));
            clear histo
            for d=1:numel(unique_distances)
                histo(:,d)=hist([per_pos(:,distances==unique_distances(d)).meanamp],bins);
            end
            bar(bins,histo,'stacked');
            x_lim=get(gca,'xlim');
            y_lim=get(gca,'ylim');

            
            for f=1:numel(fixations)
                for d=1:numel(unique_distances)
                    persession=nanmean(vertcat(per_pos(f,distances==unique_distances(d)).amp),1);
                    text_to_plot=[num2str(nanmean(persession)) '+' num2str(sterr(persession))];
                    text(x_lim(1)+diff(x_lim)/5*f,y_lim(2)-diff(y_lim)/20*d,text_to_plot);
                end
            end
            
            pd = anovan([per_pos.amp]',[[per_pos.fixation]' [per_pos.distance]'],'model','full','display','off');
            pp = anovan([per_pos.amp]',[[per_pos.fixation]' [per_pos.position]'],'model','full','display','off');
            
            text_to_plot=['p anova FxD: ' num2str(pd(1)) ', ' num2str(pd(2)) ', ' num2str(pd(3)) ];
            text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*(d+1),text_to_plot);
            text_to_plot=['p anova FxP: ' num2str(pp(1)) ', ' num2str(pp(2)) ', ' num2str(pp(3)) ];
            text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*(d+2),text_to_plot);
            
            xlabel('saccade amplitude [°]');
            ylabel('N conditions');
            title(['across targets: ' num2str(mean([per_pos(:).meanamp])) ' + ' num2str(sterr([per_pos(:).meanamp]))...
                'across sessions: ' num2str(mean(per_session.amp)) ' + ' num2str(sterr(per_session.amp))  ]);
            
            legend(cellstr(num2str(unique_distances')));
            ph_title_and_save(figure_handle,filename,plot_title,keys)
            
            %% final amplitude
            
            plot_title=[prefix ' ' monkey ' ' type_effector_short{te} ' ' condition_title ' final amplitude'];
            filename=plot_title;
            figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            
            bins=min(floor([per_pos.meanfinalamp])):1:max(ceil([per_pos.meanfinalamp]));
            clear histo
            for d=1:numel(unique_distances)
                histo(:,d)=hist([per_pos(:,distances==unique_distances(d)).meanfinalamp],bins);
            end
            bar(bins,histo,'stacked');
            x_lim=get(gca,'xlim');
            y_lim=get(gca,'ylim');

            
            for f=1:numel(fixations)
                for d=1:numel(unique_distances)
                    persession=nanmean(vertcat(per_pos(f,distances==unique_distances(d)).finalamp),1);
                    text_to_plot=[num2str(nanmean(persession)) '+' num2str(sterr(persession))];
                    text(x_lim(1)+diff(x_lim)/5*f,y_lim(2)-diff(y_lim)/20*d,text_to_plot);
                end
            end
            
            pd = anovan([per_pos.finalamp]',[[per_pos.fixation]' [per_pos.distance]'],'model','full','display','off');
            pp = anovan([per_pos.finalamp]',[[per_pos.fixation]' [per_pos.position]'],'model','full','display','off');
            
            text_to_plot=['p anova FxD: ' num2str(pd(1)) ', ' num2str(pd(2)) ', ' num2str(pd(3)) ];
            text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*(d+1),text_to_plot);
            text_to_plot=['p anova FxP: ' num2str(pp(1)) ', ' num2str(pp(2)) ', ' num2str(pp(3)) ];
            text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*(d+2),text_to_plot);
            
            xlabel('final distance [°]');
            ylabel('N conditions');
            title(['across targets: ' num2str(mean([per_pos(:).meanfinalamp])) ' + ' num2str(sterr([per_pos(:).meanfinalamp])) ]);
            
            legend(cellstr(num2str(unique_distances')));
            ph_title_and_save(figure_handle,filename,plot_title,keys)
            
            
            
            %% velocity
            
            plot_title=[prefix ' ' monkey ' ' type_effector_short{te} ' ' condition_title ' velocity'];
            filename=plot_title;
            figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            
            bins=min(floor([per_pos.meanvel]/10)*10):10:max(ceil([per_pos.meanvel]/10)*10);
            clear histo
            for d=1:numel(unique_distances)
                histo(:,d)=hist([per_pos(:,distances==unique_distances(d)).meanvel],bins);
            end
            bar(bins,histo,'stacked');
            x_lim=get(gca,'xlim');
            y_lim=get(gca,'ylim');

            
            for f=1:numel(fixations)
                for d=1:numel(unique_distances)
                    persession=nanmean(vertcat(per_pos(f,distances==unique_distances(d)).vel),1);
                    text_to_plot=[num2str(nanmean(persession)) '+' num2str(sterr(persession))];
                    text(x_lim(1)+diff(x_lim)/5*f,y_lim(2)-diff(y_lim)/20*d,text_to_plot);
                end
            end
            
            pd = anovan([per_pos.vel]',[[per_pos.fixation]' [per_pos.distance]'],'model','full','display','off');
            pp = anovan([per_pos.vel]',[[per_pos.fixation]' [per_pos.position]'],'model','full','display','off');
            
            text_to_plot=['p anova FxD: ' num2str(pd(1)) ', ' num2str(pd(2)) ', ' num2str(pd(3)) ];
            text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*(d+1),text_to_plot);
            text_to_plot=['p anova FxP: ' num2str(pp(1)) ', ' num2str(pp(2)) ', ' num2str(pp(3)) ];
            text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*(d+2),text_to_plot);
            
            xlabel('velocity [°/s]');
            ylabel('N conditions');
            title(['across targets: ' num2str(mean([per_pos(:).meanvel])) ' + ' num2str(sterr([per_pos(:).meanvel]))...
                'across sessions: ' num2str(mean(per_session.vel)) ' + ' num2str(sterr(per_session.vel))  ]);
            
            legend(cellstr(num2str(unique_distances')));
            ph_title_and_save(figure_handle,filename,plot_title,keys)
            
            
            %% precision
            
            plot_title=[prefix ' ' monkey ' ' type_effector_short{te} ' ' condition_title ' precision'];
            filename=plot_title;
            figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            
            bins=min(floor([per_pos.meanvel]/10)*10):10:max(ceil([per_pos.meanvel]/10)*10);
            clear histo
            for d=1:numel(unique_distances)
                histo(:,d)=hist([per_pos(:,distances==unique_distances(d)).meanvel],bins);
            end
            bar(bins,histo,'stacked');
            x_lim=get(gca,'xlim');
            y_lim=get(gca,'ylim');

            
            pd = anovan([per_pos.sac_prec]',[[per_pos.fixation]' [per_pos.distance]'],'model','full','display','off');
            pp = anovan([per_pos.sac_prec]',[[per_pos.fixation]' [per_pos.position]'],'model','full','display','off');
            
            text_to_plot=['p anova FxD: ' num2str(pd(1)) ', ' num2str(pd(2)) ', ' num2str(pd(3)) ];
            text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*(d+1),text_to_plot);
            text_to_plot=['p anova FxP: ' num2str(pp(1)) ', ' num2str(pp(2)) ', ' num2str(pp(3)) ];
            text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*(d+2),text_to_plot);
            
            xlabel('precision [°]');
            ylabel('N conditions');
            title(['across targets: ' num2str(mean([per_pos(:).sac_prec])) ' + ' num2str(sterr([per_pos(:).sac_prec]))]);
            
            legend(cellstr(num2str(unique_distances')));
            ph_title_and_save(figure_handle,filename,plot_title,keys)
        end
    end
end
end

function hArrow = drawArrow(p0,p1,color)
% drawArrow(p0,p1)
% Draws a simple arrow in 2D, from p0 to p1.
% from: https://de.mathworks.com/matlabcentral/fileexchange/55181-drawarrow by Matthew Kelly
%
% INPUTS:
%   p0 = [x0; y0] = position of the tail
%   p1 = [x1; y1] = position of the tip
%   color = arrow color. Optional: default is black
%       --> can be 'r','g','b','c','m','y','w', 'k' or a 1x3 color vector
%
% OUTPUTS:
%   hArrow = handle to the patch object representing the arrow
%
% Defaults:
if nargin == 2
    color = 'k';
end
% Parameters:
W1 = 0.08;   % half width of the arrow head, normalized by length of arrow
W2 = 0.014;  % half width of the arrow shaft
L1 = 0.18;   % Length of the arrow head, normalized by length of arrow
L2 = 0.13;  % Length of the arrow inset
% Unpack the tail and tip of the arrow
x0 = p0(1);
y0 = p0(2);
x1 = p1(1);
y1 = p1(2);
% Start by drawing an arrow from 0 to 1 on the x-axis
P = [...
    0, (1-L2), (1-L1), 1, (1-L1), (1-L2), 0;
    W2,    W2,     W1, 0,    -W1,    -W2, -W2];
% Scale,rotate, shift and plot:
dx = x1-x0;
dy = y1-y0;
Length = sqrt(dx*dx + dy*dy);
Angle = atan2(-dy,dx);
P = P.*repmat([Length;1],1,7);   %Scale Length*
P = [cos(Angle), sin(Angle); -sin(Angle), cos(Angle)]*P;  %Rotate
P = p0(:)*ones(1,7) + P;  %Shift
% Plot!
hArrow = patch(P(1,:), P(2,:),color);  axis equal;
set(hArrow,'EdgeColor',color);
end

function drawring(center,rin, rout,angles,color)
xin = center(1) + rin*cos(angles);
xout = center(1) + rout*cos(angles);
yin = center(2) + rin*sin(angles);
yout = center(2) + rout*sin(angles);
% Make patch
hp = patch([xout,xin],[yout,yin],color,'linestyle','none');
end
