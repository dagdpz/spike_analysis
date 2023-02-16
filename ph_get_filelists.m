function filelist=ph_get_filelists(population,keys)

[TT,idx,group_values,unique_group_values]=ph_readout_tuning_table(keys);
if isempty(unique_group_values)
    return;
end

complete_unit_list={population.unit_ID}';
[unit_valid,TM]=ismember(complete_unit_list,TT(:,idx.unitID));
population=population(unit_valid);
% complete_unit_list={population.unit_ID}';
% population_group=group_values(TM(unit_valid));
all_trialz=[population.trial];
[UC, CM, labels]=ph_get_condition_matrix(all_trialz,keys);
CP_out                      = [{'type','effector'},keys.condition_parameters];
conditions_out              = combvec(UC.type,UC.effector,CM')';

pop_monkeys={population.monkey};
%pop_target={population.target};

%keys.target=keys.batching.targets{1}; %% to_remove...

%% these here should be inputs
keys.colors.fix_offset=[keys.colors.IF; keys.colors.MF; keys.colors.IF_CS];
corrective=1;
recalibrated=0;
toplot=1;
current_effector='saccades';
epoch={'Thol',5,	0.2,    0.5}; %% for final precision etc.
datapath_recalibrated=''; %% data path for recalibrated data

%%% actually the idea of splitting into datasets is kind of obsolete now?

%% create filelist from population files
for m=1:numel(keys.monkeys)
    keys.monkey=keys.monkeys{m};
    pop_m=population(ismember(pop_monkeys,[keys.monkey '_phys']));
    pop_trials=[pop_m.trial];
    for c=1:size(conditions_out,1)
        eff=conditions_out(c,strcmp(CP_out,'effector'));
        typ=conditions_out(c,strcmp(CP_out,'type'));
        clear trpar
        for par=1:numel(CP_out)
            fn=CP_out{par};
            trpar(par,:)=[pop_trials.(fn)]==conditions_out(c,par);
            keys.tt.(fn)=conditions_out(c,par);
        end
        c_trials=pop_trials(all(trpar,1));
        [~,~,condition_title]=ph_figure_titles(keys);
        condition_title=strrep(condition_title,' ','_');
        if isempty(c_trials)
            continue;
        end
        prefix=[keys.monkey(1:3) '_T' num2str(typ) 'E' num2str(eff) '_' keys.arrangement '_' condition_title];
        if corrective==1
            prefix=[prefix '_corr'];
        else
            prefix=[prefix '_unco'];
        end
        if recalibrated==1
            prefix=[prefix '_recal'];
            base_path=datapath_recalibrated;
        else
            prefix=[prefix '_rawMP'];
            %base_path=['Y:\Data\' keys.monkey '_phys_combined_monkeypsych_TDT\'];
            base_path='Y:\Data'; %\' keys.monkey '_phys_combined_monkeypsych_TDT\'];
        end
        
        blocks=unique([[c_trials.date]',[c_trials.block]'],'rows');
        sessions=unique(blocks(:,1));
        for s=1:numel(sessions)
            session=sessions(s);
            idx=blocks(:,1)==session;
            runs=run2block([base_path filesep keys.monkey '_phys_combined_monkeypsych_TDT' filesep num2str(session)],blocks(idx,2),1);
            blocks(idx,3)=runs;
            filelist_formatted.(prefix){s,1}=session;
            filelist_formatted.(prefix){s,2}=runs;
        end
        runs0=arrayfun(@(x) num2str(x), blocks(:,1),'UniformOutput',false);
        runs1=arrayfun(@(x) num2str(floor(x/10000)), blocks(:,1),'UniformOutput',false);
        runs2=arrayfun(@(x) sprintf('%02d',mod(floor(x/100),100)), blocks(:,1),'UniformOutput',false);
        runs3=arrayfun(@(x) sprintf('%02d',mod(x,100)), blocks(:,1),'UniformOutput',false);
        runs4=arrayfun(@(x) sprintf('%02d',x), blocks(:,3),'UniformOutput',false);
        
        %% here we can make a case for, uhm, recalibrated data?
        current_filelist=strcat([base_path filesep keys.monkey filesep], runs0, keys.monkey(1:3), runs1, '-', runs2, '-', runs3, '_', runs4,'.mat');
        filelist.(prefix)=current_filelist;
    end
end
if exist([keys.path_to_save filesep 'behaviour_filelist.mat'],'file'); % make sure to save filelist without intervening with current plotting
    filelist_tmp=filelist;
    filelist_formatted_tmp=filelist_formatted;
    load([keys.path_to_save filesep 'behaviour_filelist.mat'])
    FN=fieldnames(filelist_tmp);
    for f=1:numel(FN)
        filelist.(FN{f})=filelist_tmp.(FN{f});
        filelist_formatted.(FN{f})=filelist_formatted_tmp.(FN{f});
    end
    save([keys.path_to_save filesep 'behaviour_filelist.mat'],'filelist','filelist_formatted');
    filelist=filelist_tmp;
    filelist_formatted=filelist_formatted_tmp;
    clear filelist_formatted_tmp filelist_tmp
else
    
    save([keys.path_to_save filesep 'behaviour_filelist.mat'],'filelist','filelist_formatted');
end


if corrective==1
    keys.cal.MA_selection                   ={'correct_offset',0,'success',1,'display',0,'summary',0,'keep_raw_data',1,'saccade_definition',11,'reach_1st_pos',1,'correlation_conditions',{},'nsacc_max', 20, 'sac_ini_t', 20,  'sac_end_t',  10, 'sac_min_dur', 0, 'sac_min_dur', 0};                        % if you want to run MA with specific settings
else
    keys.cal.MA_selection                   ={'correct_offset',0,'success',1,'display',0,'keep_raw_data',1,'saccade_definition',4,'reach_1st_pos',1,'correlation_conditions',{}};                        % if you want to run MA with specific settings
end


FN=fieldnames(filelist);
for f=1:numel(FN)
    FN_nomonkey{f}=FN{f}(5:end);
end
datasets=unique(FN_nomonkey);
allmonkeys=[{'All'} keys.monkeys];
for D=1:numel(datasets)
    dataset=datasets{D};
    for m=1:numel(allmonkeys)
        monkey=allmonkeys{m};
        if strcmp(monkey,'All')
            filelist_in={};
            for M=1:numel(keys.monkeys)
                monkey=keys.monkeys{M};
                filelist_in=vertcat(filelist_in,get_filelist_per_session(monkey,dataset,filelist,filelist_formatted,recalibrated,datapath_recalibrated));
            end
        else
            filelist_in=get_filelist_per_session(monkey,dataset,filelist,filelist_formatted,recalibrated,datapath_recalibrated);
        end
        MA_input={};
        for f=1:size(filelist_in,1)
            MA_input=[MA_input, {{filelist_in{f,1}, filelist_in{f,2}}},{keys.cal.MA_selection}];
        end
        data=monkeypsych_analyze_working(MA_input{:});
                
        all_data=[data{:}];
        all_task=vertcat(all_data.task);
        all_binary=vertcat(all_data.binary);
        [all_task.choice]=all_binary.choice;
        hemifields=num2cell([all_binary.eyetar_r]-[all_binary.eyetar_l]+[all_binary.hndtar_r]-[all_binary.hndtar_l]);
        [all_task.hemifield]=hemifields{:};
        u_hemifields=unique([all_task.hemifield]);
        
        
        for c= 1:size(conditions_out,1) %% not sure if dimension here is 1
            all_saccades=vertcat(all_data.(current_effector));
            positions=unique(round([all_saccades.tar_pos] - [all_saccades.fix_pos]));
            positions(isnan(positions))=[];
            fixations=unique(real([all_saccades.fix_pos])); % will need to use precision here
            distances=abs(positions);
            unique_distances=unique(distances);
            
            
            clear per_pos per_hemi per_session
            for s=1:numel(data)
                saccades=[data{s}.(current_effector)];
                states=[data{s}.states];
                raw=[data{s}.raw];
                task=[data{s}.task];
                binary=[data{s}.binary];
                [task.choice]=binary.choice;
                hemifields=num2cell([binary.eyetar_r]-[binary.eyetar_l]+[binary.hndtar_r]-[binary.hndtar_l]);
                [task.hemifield]=hemifields{:};
                                
                clear trpar
                for par=1:numel(CP_out)
                    fn=CP_out{par};
                    trpar(par,:)=[task.(fn)]==conditions_out(c,par);
                end
                trs=all(trpar,1);
                
                saccades=saccades(trs);
                states=states(trs);
                raw=raw(trs);
                hemifields=[task(trs).hemifield];
                per_session.RTs(s)=nanmean([saccades.lat])*1000;
                per_session.dur(s)=nanmean([saccades.dur])*1000;
                per_session.amp(s)=nanmean(abs([saccades.endpos]-[saccades.startpos]));
                per_session.vel(s)=nanmean([saccades.velocity]);
                per_session.pre(s)=nanmean([saccades.endpos]-[saccades.fix_pos]);
                per_session.percent_sac(s)=sum(~isnan([saccades.lat]))/sum(trs)*100;
                per_session.Ntr(s)=sum(trs);
                per_session.Ntr_l(s)=sum(real([saccades.tar_pos]-[saccades.fix_pos])<0);
                
                fix_pos=[saccades.fix_pos];
                for h=1:numel(u_hemifields)
                    tr_idx=hemifields==u_hemifields(h);
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
            
            %% store here
            mon(D,m).condition(c).per_session=per_session;
            mon(D,m).condition(c).per_hemi=per_hemi;
            mon(D,m).condition(c).per_pos=per_pos;
            
            
        end
    end
end

%% plotting
if toplot
for D=1:numel(datasets)
    dataset=datasets{D};
    for m=1:numel(allmonkeys)
        monkey=allmonkeys{m};
        for c= 1:size(conditions_out,1) %% not sure if dimension here is 1     

            
            per_session=mon(D,m).condition(c).per_session;
            per_hemi=mon(D,m).condition(c).per_hemi;
            per_pos=mon(D,m).condition(c).per_pos;
            
            eff=conditions_out(c,strcmp(CP_out,'effector'));
            typ=conditions_out(c,strcmp(CP_out,'type'));            
            for par=1:numel(CP_out)
                fn=CP_out{par};
                keys.tt.(fn)=conditions_out(c,par);
            end
            [fig_title,filename,condition_title]=ph_figure_titles(keys);
            condition_title=[monkey(1:3) '_T' num2str(typ) 'E' num2str(eff) '_' keys.arrangement '_' strrep(condition_title,' ','_')];            
            if corrective
                condition_title=['Corr ' condition_title];
            end
            if recalibrated
                condition_title=['Recal ' condition_title];
            end
            
            %% trials per position
            plot_title=[condition_title ' percent saccades detected'];
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
            plot_title=[condition_title ' percent trials per position'];
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
            plot_title=[condition_title ' accuracy'];
            filename=plot_title;
            figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            Colors=keys.colors.fix_offset/256;
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
            plot_title=[condition_title ' accuracy per session'];
            filename=plot_title;
            figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            Colors=keys.colors.fix_offset/256;
            angles=0:pi/180:2*pi;
            
            hold on
            for f=1:numel(fixations)
                for p=1:numel(positions)
                    r=5;
                    x0=positions(p)+fixations(f);
                    scatter(real(x0),imag(x0),'+','markeredgecolor','k');
                    plot(real(x0)+r*sin(angles), imag(x0)+r*cos(angles),'k');
                    scatter(real([per_pos(f,p).sac_end]),imag([per_pos(f,p).sac_end]),30,Colors(f,:),'filled');
                    nanmean(per_pos(f,p).dur);
                end
            end
            axis equal
            ph_title_and_save(figure_handle,filename,plot_title,keys)
            
            %% accuracy per session
            plot_title=[condition_title ' final accuracy per session'];
            filename=plot_title;
            figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            Colors=keys.colors.fix_offset/256;
            angles=0:pi/180:2*pi;
                        
            hold on
            for f=1:numel(fixations)
                for p=1:numel(positions)
                    r=5;
                    x0=positions(p)+fixations(f);
                    scatter(real(x0),imag(x0),'+','markeredgecolor','k');
                    plot(real(x0)+r*sin(angles), imag(x0)+r*cos(angles),'k');
                    plot([[per_pos(f,p).sac_end]' , [per_pos(f,p).eye_in_thol]']','color',Colors(f,:));
                    scatter(real([per_pos(f,p).eye_in_thol]),imag([per_pos(f,p).eye_in_thol]),15,Colors(f,:),'filled');
                end
            end
            axis equal
            ph_title_and_save(figure_handle,filename,plot_title,keys)
            
            
            %% accuracy per session
            Colors=keys.colors.fix_offset/256;
            angles=0:pi/180:2*pi;
            for s=1:numel(data)
                plot_title=[condition_title ' final accuracy, session ' num2str(s)];
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
                            plot([[per_pos(f,p).sac_end_all{s}]' , [per_pos(f,p).eye_in_thol_all{s}]']','color',Colors(f,:));
                            scatter(real([per_pos(f,p).eye_in_thol_all{s}]),imag([per_pos(f,p).eye_in_thol_all{s}]),15,Colors(f,:),'filled');
                        end
                    end
                end
                axis equal
                ph_title_and_save(figure_handle,filename,plot_title,keys)
            end
            
            %% precision per session
            plot_title=[condition_title ' precision per session'];
            filename=plot_title;
            figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            Colors=keys.colors.fix_offset/256;
            
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
                end
            end
            axis equal
            ph_title_and_save(figure_handle,filename,plot_title,keys)
            
            
            %% precision per session
            plot_title=[condition_title ' precision of session means'];
            filename=plot_title;
            figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            Colors=keys.colors.fix_offset/256;
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
                end
            end
            axis equal
            ph_title_and_save(figure_handle,filename,plot_title,keys)
            
            %% final precision per session
            plot_title=[condition_title ' final precision per session'];
            filename=plot_title;
            figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            Colors=keys.colors.fix_offset/256;
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
            plot_title=[condition_title ' final precision of session means'];
            filename=plot_title;
            figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            Colors=keys.colors.fix_offset/256;
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
            
            plot_title=[condition_title ' duration'];
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
            
            plot_title=[condition_title ' latency'];
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
                num2str(round(nanmean(per_hemi(2).RTs))) ' + ' num2str(round(sterr(per_hemi(2).RTs))) ' L/R, p=' num2str(round(p_lr*1000)/1000)]);
            
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
            
            plot_title=[condition_title ' amplitude'];
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
            
            plot_title=[condition_title ' final amplitude'];
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
            
            plot_title=[condition_title ' velocity'];
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
            
            plot_title=[condition_title ' precision'];
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

aaa=1;
%  the idea would be to shorten plotting part substantially
%     function plot_for_all_positions
% 
%             for f=1:numel(fixations)
%                 for p=1:numel(positions)
%                     r=5;
%                     x0=positions(p)+fixations(f);
%                     scatter(real(x0),imag(x0),'+','markeredgecolor','k');
%                     plot(real(x0)+r*sin(angles), imag(x0)+r*cos(angles),'k');
%                     x0=nanmean([per_pos(f,p).(stuct_to_plot)]);
%                     r=nanstd([per_pos(f,p).(stuct_to_plot)]);
%                     plot(real(x0)+r*sin(angles) , imag(x0)+r*cos(angles),'color',Colors(f,:));
%                 end
%             end
%     end

%% oculomotor readout 
% 
% %curius
% CHL=mon(2,2).condition(2).per_hemi(1).RTs;
% CHR=mon(2,2).condition(2).per_hemi(2).RTs;
% INL=mon(1,2).condition(1).per_hemi(1).RTs;
% INR=mon(1,2).condition(1).per_hemi(2).RTs;
% [~,p_cur_L]=ttest(INL,CHL)
% [~,p_cur_R]=ttest(INR,CHR)
% 
% %curius
% CHL=mon(2,3).condition(2).per_hemi(1).RTs;
% CHR=mon(2,3).condition(2).per_hemi(2).RTs;
% INL=mon(1,3).condition(1).per_hemi(1).RTs;
% INR=mon(1,3).condition(1).per_hemi(2).RTs;
% [~,p_lin_L]=ttest(INL,CHL)
% [~,p_lin_R]=ttest(INR,CHR)
end
%end
%end

function filelist_in=get_filelist_per_session(monkey,dataset,filelist,filelist_formatted,recalibrated,datapath_recalibrated)
taskfield=[monkey(1:3) '_' dataset];
mainfolder=filelist.(taskfield){1};
dashstr=strfind(mainfolder,'\');
mainfolder=mainfolder(1:dashstr(end));
filelist_in=filelist_formatted.(taskfield);
filelist_in(:,1)=cellstr(num2str(vertcat(filelist_in{:,1})));
filelist_in(:,1)=strcat(mainfolder, filelist_in(:,1));

if recalibrated==1
    aa=dir([datapath_recalibrated filesep monkey filesep]);
    sessions={aa([aa.isdir]).name};
    sessions=sessions(3:end);
    for s=1:numel(sessions)
        filelist_in{s,1}=[datapath_recalibrated filesep monkey filesep sessions{s}];
    end
end
end


function output_blocks_or_runs=run2block(base_path,given_blocks_or_runs,reverse)
folder_content=dir([base_path filesep '*.mat']);
complete_filelist=vertcat(folder_content.name);
filelist_blocks=str2num(complete_filelist(:,end-5:end-4));
filelist_runs=str2num(complete_filelist(:,end-14:end-13));

if strcmp(reverse,'reverse') || reverse==1
    [~,idx] =ismember(given_blocks_or_runs,filelist_blocks);
    output_blocks_or_runs(idx~=0)=filelist_runs(idx(idx~=0));
else
    [~,idx] =ismember(given_blocks_or_runs,filelist_runs);
    output_blocks_or_runs(idx~=0)=filelist_blocks(idx(idx~=0));
end
output_blocks_or_runs(idx==0)=0;
end



function drawring(center,rin, rout,angles,color)
xin = center(1) + rin*cos(angles);
xout = center(1) + rout*cos(angles);
yin = center(2) + rin*sin(angles);
yout = center(2) + rout*sin(angles);
% Make patch
hp = patch([xout,xin],[yout,yin],color,'linestyle','none');
end
