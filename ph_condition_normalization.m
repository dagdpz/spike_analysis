function [pop_out, condition, condition_per_trial, pref_valid]=ph_condition_normalization(population,keys)
%% assume arrangement was chosen already
switch keys.normalization_field
    case 'PO'
condition_per_trial=struct;
    case 'ON'
condition=struct;
end
K=keys.(keys.normalization_field);
%% define conditions to look at
all_trialz=[population.trial];
%% reduce all_trials by conditions_to_plot for population
%% baseline is now per unit (and not per trial any more)
%% type_effector_short?

tr_con=ismember([all_trialz.completed],keys.cal.completed);
[whatisthis]=ph_arrange_positions_and_plots(keys,all_trialz(tr_con));
%% finding positions and fixations
positions=unique(vertcat(whatisthis.trial.position),'rows');

per_trial.types       =[all_trialz.type];
per_trial.effectors   =[all_trialz.effector];
per_trial.hands       =[all_trialz.reach_hand];
per_trial.choice      =[all_trialz.choice];
per_trial.perturbation=[all_trialz.perturbation];
per_trial.hemifield   =[whatisthis.trial.hemifield];
per_trial.perturbation(ismember(per_trial.perturbation, keys.cal.perturbation_groups{1}))=0;
per_trial.perturbation(ismember(per_trial.perturbation, keys.cal.perturbation_groups{2}))=1;

u_hemifields    =unique(per_trial.hemifield);
u_types         =unique(per_trial.types);
u_effectors     =unique(per_trial.effectors);
u_hands         =unique(per_trial.hands);
u_choice        =unique(per_trial.choice);
u_perturbation  =unique(per_trial.perturbation);
u_perturbation  =u_perturbation(~isnan(u_perturbation));


type_effectors      = combvec(u_types,u_effectors)';

u_types     =unique(type_effectors(:,1))';
u_effectors =unique(type_effectors(:,2))';

CP_out  ={'effector','reach_hand','choice','perturbation'};
conditions_out              = combvec(u_effectors,u_hands,u_choice, u_perturbation)';

%adjust fixations... should be done inside ph_arrange_positions_and_plots
%!?
fixations_temp=unique(vertcat(whatisthis.trial.fixation),'rows');
fix_temp_idx=true(size(fixations_temp,1),1);
for x=1:size(fixations_temp,1)-1
    if any(all(abs(bsxfun(@minus,fixations_temp(x+1:end,:),fixations_temp(x,:)))<4,2)) %% precision....
        fix_temp_idx(x)=false;
    end
end
fixations=fixations_temp(fix_temp_idx,:);
clear whatisthis

%% Convert to ipsi/contra, Baseline subtraction, normalization, re-ordering, and gaussian RFs

%% change order:

%% define normalization matrix

CM      =1;
CP      ={'completed'};
switch K.normalization
    case 'by_position'
        CM=combvec(u_effectors,u_hands,u_choice, u_perturbation,1:size(positions,1));
        CP  ={'effector','reach_hand','choice','perturbation','p'}; %% p would need to be added,
    case 'by_condition' % aka 'by_hemifield'
        CM=combvec(u_effectors,u_hands,u_choice, u_perturbation,u_hemifields);
        CP  ={'effector','reach_hand','choice','perturbation','hemifield'};
    case 'by_perturbation'
        CM=combvec(u_effectors, u_perturbation);
        CP  ={'effector','perturbation'};
    case 'by_hand'
        CM=combvec(u_effectors,u_hands);
        CP  ={'effector','reach_hand'};
    case {'by_effector','percent_change','z_score'}
        CM=combvec(u_effectors);
        CP  ={'effector'};
    case {'by_type','none'} %this is standard
end

        pref_valid=false(numel(u_types),numel(population));
for u=1:numel(population)
    clear trcon tr
    tr_con=ismember([population(u).trial.completed],keys.cal.completed) & [population(u).trial.accepted];
    [pop]=ph_arrange_positions_and_plots(keys,population(u).trial(tr_con),population(u));
    
    [pop.trial(ismember([pop.trial.perturbation], keys.cal.perturbation_groups{1})).perturbation]=deal(0);
    [pop.trial(ismember([pop.trial.perturbation], keys.cal.perturbation_groups{2})).perturbation]=deal(1);
    pop_out(u)=pop;
    
    %% assign position index?
        
    %% normalization factor and baseline (per trial!)
    norm_factor=ones(size(pop.trial)); %% instead of uxn, make it txcH, same for baseline here?, could still be by trial (for plotting purposes only?)
    baseline=zeros(size(pop.trial)); %% instead of uxn, make it txcH, same for baseline here?, could still be by trial (for plotting purposes only?)
       
    for t=1:numel(u_types)
        typ=u_types(t);
        eff=type_effectors(type_effectors(:,1)==typ,2);
        keys=ph_get_epoch_keys(keys,typ,eff(1),numel(eff)>1); %taking eff(1) is admittedlzy sloppy here
        
        trtyp=[pop.trial.type]==typ;
        per_epoch=vertcat(pop.trial(trtyp).epoch); %only current type!?
        present_epochs={per_epoch(1,:).state};
        
        DN_epoch    =find(ismember(present_epochs,K.epoch_for_normalization));
        RF_epoch    =find(ismember(present_epochs,K.epoch_RF));
        BL_epoch    =find(ismember(present_epochs,K.epoch_BL));
        PF_epoch    =find(ismember(present_epochs,K.epoch_PF));
        gaussian_bl_epoch       =find(ismember(present_epochs,K.epoch_GB));
                        
        trpar   =trtyp;
        if strcmp(K.normalization,'by_all_trials')
            trpar=true(size([pop.trial]));
        end      
        
        %% normalization factor (condition wise)
        for n=1:size(CM,2)
            %clear trpar
            for p=1:size(CM,1)
                fn=CP{p};
                trcon(p,:)=[pop.trial.(fn)]==CM(p,n);
            end
            tr(:,n)=all(trcon,1) & trpar; 
            if ~any(tr(:,n)) %
                continue
            end
            per_epoch=vertcat(pop.trial(tr(:,n)).epoch);
            if strcmp(K.normalization,'percent_change') %not sure if this really percent change, the scale should be logarithmic?
                norm_factor(tr(:,n))=nanmean([per_epoch(:,BL_epoch).FR per_epoch(:,DN_epoch).FR NaN])/100;
                baseline(tr(:,n))   =nanmean([per_epoch(:,BL_epoch).FR per_epoch(:,DN_epoch).FR NaN]);
            else
                baseline(tr(:,n))   =nanmean(vertcat(per_epoch(:,BL_epoch).FR));
                norm_factor(tr(:,n))=nanmean(vertcat(per_epoch(:,DN_epoch).FR));
            %norm_factor=nanmean([per_epoch(:,DN_epoch).FR NaN]);
            end
            
        end
        
        %% per trial baseline/normalization
        if K.baseline_per_trial            
            per_epoch=vertcat(pop.trial(trtyp).epoch);
            baseline=[per_epoch(:,BL_epoch).FR];
        end
        %% validate
        if strcmp(K.normalization,'none') || isempty(DN_epoch) %DN_epoch is empty? %% check if this works, for no multiplicative normalization
            norm_factor(:)= 1;
        end        
        if ~K.FR_subtract_baseline  || isempty(BL_epoch)
            baseline(:)=0;
        end                
        %% not sure how (and why) to implement this 
%         if all(isnan(norm_factor(t,:))) &&  all(all(isnan(baseline)))
%             continue
%         end
        
        %% z-scoring (way too complicated implementation)
        if strcmp(K.normalization,'z_score')
            u_blocks=unique([pop.trial.block]);
            AT_tmp={};
            prev_block_end=0;
            trial_onset_times=zeros(size(pop.trial));
            for b=1:numel(u_blocks)
                b_idx=[pop.trial.block]==u_blocks(b);
                AT_tmp=[AT_tmp arrayfun(@(x) single(x.arrival_times + x.trial_onset_time + prev_block_end),pop.trial(b_idx),'uniformoutput',false)];
                last_trial_idx=find(b_idx,1,'last');
                trial_onset_times(b_idx)=[pop.trial(b_idx).trial_onset_time]+prev_block_end;
                prev_block_end=prev_block_end+pop.trial(last_trial_idx).trial_onset_time+pop.trial(last_trial_idx).states_onset(pop.trial(last_trial_idx).states==98);
            end
            AT_tmp=unique(vertcat(AT_tmp{:}));
            PSTH_ms =ceil(double(trial_onset_times(1))/keys.PSTH_binwidth)*keys.PSTH_binwidth:0.001:...
                floor(double(trial_onset_times(end)+pop.trial(end).states_onset(pop.trial(end).states==98))/keys.PSTH_binwidth)*keys.PSTH_binwidth;
            SD_ms= conv(hist(AT_tmp,PSTH_ms),normpdf(-5*keys.gaussian_kernel:0.001:5*keys.gaussian_kernel,0,keys.gaussian_kernel),'same');
            
            t_idx=false(size(PSTH_ms));
            for t_temp=find([pop.trial.type]==typ & [pop.trial.effector]==eff)
                min_temp=ceil(trial_onset_times(t_temp)/keys.PSTH_binwidth)*keys.PSTH_binwidth;
                max_temp=floor((trial_onset_times(t_temp)+double(pop.trial(t_temp).states_onset(pop.trial(t_temp).states==50)))/keys.PSTH_binwidth)*keys.PSTH_binwidth;
                t_idx(PSTH_ms+0.0002>=min_temp & PSTH_ms+0.0002<max_temp)=true;
            end
            SD=nanmean(reshape(SD_ms(t_idx),keys.PSTH_binwidth/0.001,sum(t_idx)/keys.PSTH_binwidth*0.001),1);
            norm_factor(:)=deal(std(SD));
            baseline(:)=deal(nanmean(SD));
        end
        
        %% this line here, what to do with it?
        %norm_factor(t,:)=deal(max([norm_factor(t,:) 0])); % now we always normalize to maximum condition, 0 makes sure some value is there..
        %% correct normalization factors if they are too low  --> still normalize differently in different trials, or in general add 1 ? 
        if any(norm_factor(:)<1)
            %baseline=baseline+1-nanmean(norm_factor(t,:)); %not sure what this was meant for?
            norm_factor(:)=deal(1);
        end
        
        clear per_epoch_FR
        tr=trtyp;
        per_epoch=vertcat(pop.trial(tr).epoch)';
        for e=1:size(per_epoch,1) %% looping may not be elegant, but at least confusions are avoided
            per_epoch_FR(e,:)=([per_epoch(e,:).FR]-baseline(tr))./norm_factor(tr); %not sure if dimension is correct here!!!!!!!!!!!!
        end
        per_epoch_FR=num2cell(per_epoch_FR);
        [per_epoch.FR]=per_epoch_FR{:};
        for trl=find(tr)
            [pop_out(u).trial(trl).epoch.FR]=per_epoch_FR{:,trl};
        end
    end
    
    
  
        if strcmp(keys.normalization_field,'AN')
            continue
        end
        
        
        %% PSTH calculation
    
    pop=ph_LR_to_CI(keys,pop_out(u)); % Convert to ipsi/contra? not sure if necessary here 
    for t=1:numel(u_types)
        typ=u_types(t);
        eff=type_effectors(type_effectors(:,1)==typ,2);
        keys=ph_get_epoch_keys(keys,typ,eff(1),numel(eff)>1); %taking eff(1) is admittedlzy sloppy here   
        trtyp=[pop.trial.type]==typ;
        
        
        
        %% preferred and unpreferred location (taken from instructed trials only! (?) Not necessary though
        % So far, prefernce is estimated per type
        
        FR_for_pref=NaN(size(positions,1),1);
        per_epoch=vertcat(pop.trial(trtyp).epoch);
        for p=1:size(positions,1)
            tr_hemi=all(abs(bsxfun(@minus,vertcat(pop.trial.position),positions(p,:)))<1.5,2) & [pop.trial.choice]'==0;
            FR_for_pref(p)=nanmean([per_epoch((tr_hemi(trtyp)),PF_epoch).FR])-nanmean(baseline(tr_hemi & trtyp'));
        end
%         if ischar(unique_group_values{g}) && strcmp(unique_group_values{g},'su') %very specific rule, be aware of this one!
%             [~,pref_idx]=min(FR_for_pref);
%         else
            [~,pref_idx]=max(FR_for_pref);
%         end
        [~,unpref_idx]=min(abs(nanmean([positions(:,1)+positions(pref_idx,1) positions(:,2)-positions(pref_idx,2)],2)));
         pref_valid(t,u)=true;
        for ch=u_choice
            pref_valid(t,u)=pref_valid(t,u) && pref_idx~=unpref_idx && ...
                sum(all(abs(bsxfun(@minus,vertcat(pop.trial(trtyp).position),positions(pref_idx,:)))<1.5,2)   & [pop.trial(trtyp).choice]'==ch) >=keys.cal.min_trials_per_condition && ...
                sum(all(abs(bsxfun(@minus,vertcat(pop.trial(trtyp).position),positions(unpref_idx,:)))<1.5,2) & [pop.trial(trtyp).choice]'==ch) >=keys.cal.min_trials_per_condition;
        end
        
   % for e=1:numel(u_effectors)
%         eff=u_effectors(e);
%         typ=type_effectors(tye,1);
%         eff=type_effectors(tye,2);
%         t=find(u_types==typ);
%         e=find(u_effectors==eff);
%         keys=ph_get_epoch_keys(keys,typ,eff,sum(type_effectors(:,1)==typ)>1);
        
     %   conditions_eff=conditions_out(conditions_out(:,1)==eff,2:end);
    
      
        
        %% average PSTH per unit
        for c=1:size(conditions_out,1)
            eff=conditions_out(1,strcmp(CP_out,'effector'));
            %c=ce;
            %c=find(ismember(conditions_out,[eff conditions_eff(ce,:)],'rows')); %% c is conditions, n is including hemifields -1,0,1
            clear trpar
            
            for par=1:numel(CP_out)
                fn=CP_out{par};
                trpar(par,:)=[pop.trial.(fn)]==conditions_out(c,par); %condition_matrix_wohf?
            end
            trpar(end+1,:)=trtyp;
            tr_con=all(trpar,1);
            per_epoch=vertcat(pop.trial(tr_con).epoch); % take already normalized values???
            
            for w=1:size(keys.PSTH_WINDOWS,1)
                for f=1:numel(u_hemifields) %hemifield
                    tr_hemi=[pop.trial.hemifield]==u_hemifields(f);
                    switch keys.normalization_field
                        case 'PO'
                            ix = tr_con & tr_hemi;
                            %n=max([1,find(ismember(CM(:,1:4),conditions_out(c,:),'rows') & CM(:,end)==u_hemifields(f))]);
                            condition(t,c).per_hemifield(f).window(w).unit(u).average_spike_density=...
                                ph_spike_density(pop.trial(ix),w,keys,baseline(ix),norm_factor(ix));
                            condition(t,c).per_hemifield(f).effector=eff;
                            for x=1:size(fixations,1)
                                tr_fix=all(abs(bsxfun(@minus,vertcat(pop.trial.fixation),fixations(x,:)))<4,2)'; %4 is precision --> add to keys?
                                ix = tr_con & tr_hemi & tr_fix;
                                condition(t,c).per_hf_fixation(f,x).window(w).unit(u).average_spike_density=...
                                    ph_spike_density(pop.trial(ix),w,keys,baseline(ix),norm_factor(ix));
                                condition(t,c).per_hf_fixation(f,x).fixation=fixations(x,:);
                                condition(t,c).per_hf_fixation(f,x).effector=eff;
                                condition(t,c).per_hf_fixation(f,x).position=[u_hemifields(f) 0];
                                condition(t,c).per_hf_fixation(f,x).sign.unit(u)=1;
                            end
                            
                        case 'ON'
                            %% alternative per trial
                            trials_for_SD=find(tr_con & tr_hemi);
                            temp_window=keys.PSTH_WINDOWS(w,:);
                            keys.PSTH_WINDOWS{w,3}=keys.PSTH_WINDOWS{w,3}-keys.n_consecutive_bins_significant*keys.PSTH_binwidth;
                            keys.PSTH_WINDOWS{w,4}=keys.PSTH_WINDOWS{w,4}+keys.n_consecutive_bins_significant*keys.PSTH_binwidth;
                            for tr=1:numel(trials_for_SD)
                                tr_tmp=trials_for_SD(tr);
                                condition_per_trial(t,c).per_hemifield(f).window(w).unit(u).average_spike_density(tr,:)=...
                                    ph_spike_density(pop.trial(tr_tmp),w,keys,baseline(tr_tmp),norm_factor(tr_tmp));
                            condition_per_trial(t,c).per_hemifield(f).unit(u).epoch_averages=...
                                nanmean(reshape([per_epoch(tr_hemi(tr_con),:).FR],size(per_epoch(tr_hemi(tr_con),:))),1);
                            end
                            if isempty(trials_for_SD)
                                condition_per_trial(t,c).per_hemifield(f).window(w).unit(u).average_spike_density(1,:)=...
                                    NaN(size(ph_spike_density(pop.trial(1),w,keys,baseline(1),norm_factor(1))));
                            condition_per_trial(t,c).per_hemifield(f).unit(u).epoch_averages=NaN(1,size(keys.EPOCHS,1));
                            end
                            keys.PSTH_WINDOWS(w,:)=temp_window;
                    end
                end
                
                
                if strcmp(keys.normalization_field,'PO') %% && plot_position??
                    for p=1:size(positions,1)
                        tr_pos=all(abs(bsxfun(@minus,vertcat(pop.trial.position),positions(p,:)))<1.5,2)';
                        ix = tr_con & tr_pos;
                        
                        condition(t,c).per_position(p).window(w).unit(u).average_spike_density= ph_spike_density(pop.trial(ix),w,keys,baseline(ix),norm_factor(ix));
                        condition(t,c).per_position(p).position=positions(p,:);
                        condition(t,c).per_position(p).fixation=1;
                        condition(t,c).per_position(p).effector=eff;
                        % problem here if tr_con is empty (which it can be,
                        % because we are not requiring every condition to
                        % be valid for all units
                        ix = tr_pos(tr_con);
                        if any(ix) && ttest([per_epoch(ix,PF_epoch).FR], [per_epoch(ix,BL_epoch).FR])==1
                            condition(t,c).per_position(p).sign.unit(u)=sign(mean([per_epoch(ix,PF_epoch).FR]- [per_epoch(ix,BL_epoch).FR]));
                        else
                            condition(t,c).per_position(p).sign.unit(u)=0;
                        end
                        %% supposedly matches with target position precision...
                        for x=1:size(fixations,1)
                            tr_fix=all(abs(bsxfun(@minus,vertcat(pop.trial.fixation),fixations(x,:)))<4,2)';
                            ix = tr_con & tr_pos & tr_fix;
                            
                            % within condition normalization here...
                            % norm_fact=nanmean([per_epoch(trpos(tr_con) & trfix(tr_con),BL_epoch).FR per_epoch(trpos(tr_con) & trfix(tr_con),DN_epoch).FR NaN])/100;
                            % base_line=nanmean([per_epoch(trpos(tr_con) & trfix(tr_con),BL_epoch).FR per_epoch(trpos(tr_con) & trfix(tr_con),DN_epoch).FR NaN]);
                            % condition(t,c).per_position_fixation(p,x).window(w).unit(u).average_spike_density=...
                            % ph_spike_density(pop.trial(tr_con & trpos & trfix),w,keys,base_line,norm_fact);
                            
                            condition(t,c).per_position_fixation(p,x).window(w).unit(u).average_spike_density=...
                                ph_spike_density(pop.trial(ix),w,keys,baseline(ix),norm_factor(ix));
                            condition(t,c).per_position_fixation(p,x).position=positions(p,:);
                            condition(t,c).per_position_fixation(p,x).fixation=fixations(x,:);
                            condition(t,c).per_position_fixation(p,x).effector=eff;
                            condition(t,c).per_position_fixation(p,x).sign.unit(u)=0;
                            % problem here if tr_con is empty (which it can be,
                            % because we are not requiring every condition to
                            % be valid for all units
                            ix = tr_pos(tr_con) & tr_fix(tr_con);
                            if any(ix) && ttest([per_epoch(ix,PF_epoch).FR], [per_epoch(ix,BL_epoch).FR])==1
                                condition(t,c).per_position_fixation(p,x).sign.unit(u)=sign(mean([per_epoch(ix,PF_epoch).FR]- [per_epoch(ix,BL_epoch).FR]));
                            end
                        end
                        ix = tr_con & tr_pos;
                        if p==unpref_idx
                            condition(t,c).per_preference(1).window(w).unit(u).average_spike_density=...
                                ph_spike_density(pop.trial(ix),w,keys,baseline(ix),norm_factor(ix));
                            condition(t,c).per_preference(1).effector=eff;
                        end
                        if p==pref_idx
                            condition(t,c).per_preference(2).window(w).unit(u).average_spike_density=...
                                ph_spike_density(pop.trial(ix),w,keys,baseline(ix),norm_factor(ix));
                            condition(t,c).per_preference(2).effector=eff;
                        end
                    end
                end
            end
            
            %% gaussian response fields (needs cleanup)
            if ~K.plot_RF
                continue;
            end
            
            %% gaussian fit settings
            fitsettings.sd_max_x=12;
            fitsettings.sd_x_min_ratio=0.125;%0.125;
            fitsettings.sd_max_y=fitsettings.sd_max_x;
            fitsettings.sd_xy_min_ratio=0.25;
            fitsettings.sd_xy_max_ratio=1;
            fitsettings.sd_y_min_ratio=fitsettings.sd_x_min_ratio;
            fitsettings.xout=[-30:30]; % range for gaussian fit
            fitsettings.yout=[-15:15];
            fitsettings.range_factor=1;
            fitsettings.fittypes=K.fittypes;
            
            
            zin_per_pos=NaN(size(positions,1),1);
            
            gaussian_baseline=vertcat(per_epoch(:,gaussian_bl_epoch).FR);
            if isempty(gaussian_baseline)
                gaussian_baseline=zeros(size(per_epoch,1),1);
            end
            %% supposedly matches with target position precision...
            for p=1:size(positions,1) %for FR plot only
                tr=all(abs(bsxfun(@minus,vertcat(pop.trial(tr_con).position),positions(p,:)))<1.5,2);
                if sum(tr)==0; continue; end
                zin_per_pos(p)=nanmean(vertcat(per_epoch(tr,RF_epoch).FR)-gaussian_baseline(tr));
            end
            gaussian_positions=vertcat(pop.trial(tr_con).position);
            zin=vertcat(per_epoch(:,RF_epoch).FR);
            
            %%double because single in Linus/Curius Pulvinar gaze ?
            %nonanidx=~isnan(zin_per_pos);
            %FR_tmp=struct('FR',num2cell(zin_per_pos(nonanidx)),'x',num2cell(positions(nonanidx,1)),'y',num2cell(positions(nonanidx,2)));
            FR_tmp=struct('FR',num2cell(zin_per_pos),'x',num2cell(positions(:,1)),'y',num2cell(positions(:,2)));
            RF_tmp=ph_fit_target_positions_2D(gaussian_positions(:,1),gaussian_positions(:,2),zin,gaussian_baseline,fitsettings);
            
            condition(t,c).fitting.unit(u).parameters=RF_tmp;
            condition(t,c).fitting.unit(u).positions =FR_tmp;
            
        end
    %end
    %     end
    end
end
end
