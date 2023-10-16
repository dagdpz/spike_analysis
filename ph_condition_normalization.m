function [pop_out, condition, condition_per_trial, pref_valid]=ph_condition_normalization(population,trials,keys,UC,CM)
%% assume arrangement was chosen already!
%% position precision relevant for this part here as well!

%% placeholders for non-used outputs
switch keys.normalization_field
    case {'PO','RE','PR','RF'}
        condition_per_trial=struct;
    case {'ON','RT'}
        condition=struct;
end
K=keys.(keys.normalization_field);


K.FR_subtract_baseline=~strcmp(K.epoch_BL,'none');
CP_out                      = [{'effector'},keys.condition_parameters];
conditions_out              = combvec(UC.effector,CM')';

%% define normalization matrix
CM      ={};
CP      ={'completed'};
switch K.normalization
    case 'by_position'
        %CM=combvec(UC.effector,UC.reach_hand,UC.choice, UC.perturbation,1:size(UC.position,1));
        CP  =[{'effector'},keys.condition_parameters,{'pos_index'}];
    case 'by_condition' % aka 'by_hemifield'
        CP  =[{'effector'},keys.condition_parameters,{'hemifield'}];
    case 'by_condition_DisTask' %
        CP  ={'stimulustype','difficulty','choice','success','hemifield'};
    case 'by_perturbation'
        CP  ={'effector','perturbation'};
    case 'by_hand'
        CP  ={'effector','reach_hand'};
    case {'by_effector','percent_change','z_score'}
        CP  ={'effector'};
    case {'by_type','none'} %this is standard
end
for par=1:numel(CP)
    CM=[CM UC.(CP{par})];
end
CM=combvec(CM{:});


%% Convert to ipsi/contra, Baseline subtraction, normalization, re-ordering, and gaussian RFs
pref_valid=false(numel(UC.type),numel(population));
for u=1:numel(population)
    clear trcon tr
    pop=population(u);
    ut=ph_get_unit_trials(pop,trials);
    %ut_con=ismember([ut.completed],keys.cal.completed) & [ut.accepted];
    ut_acc=[ut.accepted];
    alltypes=[ut.type];
    ut=ut(ut_acc);
    pop.trial=pop.trial(ut_acc);
    % why store?
    pop_out(u)=pop;
    
    %% normalization factor and baseline (per trial!)
    norm_factor=ones(size(ut));
    baseline=zeros(size(ut));
    PF_epoch_FR=NaN(size(ut));
    
    for t=1:numel(UC.type)
        typ=UC.type(t);
        eff=UC.type_effector(UC.type_effector(:,1)==typ,2);
        keys=ph_get_epoch_keys(keys,typ,eff(1),numel(eff)>1); %taking eff(1) is admittedly sloppy here
        
        present_epochs=keys.EPOCHS_PER_TYPE{typ}(:,1); % now we pretend all epochs are there...
        DN_epoch{t}    =ismember(present_epochs,K.epoch_DN);
        RF_epoch{t}    =ismember(present_epochs,K.epoch_RF);
        BL_epoch{t}    =ismember(present_epochs,K.epoch_BL);
        PF_epoch{t}    =ismember(present_epochs,K.epoch_PF);
        
        
        ut_typ=[ut.type]==typ;
        pop.epochs_per_type{typ}=pop.epochs_per_type{typ}(ut_acc(alltypes==typ),:);
        
        %trpar   =ut_con(trtyp);
        
        %per_epoch=pop.epochs_per_type{typ}(trpar,:); %only current type!
        
        if strcmp(K.normalization,'by_all_trials') %% this one does not work right now
            trpar=true(size([pop.trial]));
        end
        
        
        
        %% normalization factor (condition wise)
        for n=1:size(CM,2)
            for p=1:size(CM,1)
                fn=CP{p};
                trcon(p,:)=[ut.(fn)]==CM(p,n);
            end
            tr(:,n)=all(trcon,1) & ut_typ;
            if ~any(tr(:,n)) %
                continue
            end
            per_epoch=pop.epochs_per_type{typ}(tr(ut_typ,n),:);
            if strcmp(K.normalization,'percent_change') %not sure if this really percent change, the scale should be logarithmic?
                norm_factor(tr(:,n))=nanmean([per_epoch(:,BL_epoch{t}) per_epoch(:,DN_epoch{t}) NaN])/100;
                baseline(tr(:,n))   =nanmean([per_epoch(:,BL_epoch{t}) per_epoch(:,DN_epoch{t}) NaN]);
            else
                if any(BL_epoch{t})
                    baseline(tr(:,n))   =nanmean(vertcat(per_epoch(:,BL_epoch{t})));
                end
                if any(DN_epoch{t})
                    norm_factor(tr(:,n))=nanmean(vertcat(per_epoch(:,DN_epoch{t})));
                end
            end
            if any(PF_epoch{t})
                PF_epoch_FR(tr(:,n))   =per_epoch(:,PF_epoch{t});
            end
        end
        
        %% per trial baseline/normalization! should divisive be by trial as well?
        if K.baseline_per_trial
            per_epoch=pop.epochs_per_type{typ}(ut_acc(ut_typ),:);
            baseline(ut_typ)=per_epoch(:,BL_epoch{t});
        end
        %% validate
        if strcmp(K.normalization,'none') || sum(DN_epoch{t})==0
            norm_factor(:)= 1;
        end
        if (~K.FR_subtract_baseline  || sum(BL_epoch{t})==0) && ~strcmp(K.normalization,'percent_change') % 20211102 percent_change has subtract baseline included basically
            baseline(:)=0;
        end
        
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
        
        %% correct normalization factors if they are too low  --> still normalize differently in different trials, or in general add 1 ?
        if any(norm_factor(:)<1)
            baseline=baseline+1-nanmean(norm_factor); %when we use both baseline and norm_factor, cant simply only saturate norm_factor
            norm_factor(:)=deal(1);
        end
        
        clear per_epoch_FR
        per_epoch=pop.epochs_per_type{typ}(ut_acc(ut_typ),:);
        per_epoch_FR=(per_epoch-repmat(double(baseline(ut_acc(ut_typ))'),1,size(per_epoch,2)))./repmat(norm_factor(ut_acc(ut_typ))',1,size(per_epoch,2));
        pop_out(u).epochs_per_type{typ}=per_epoch_FR;
    end
    
    % skip PSTH computation per condition if we are only interested in FR epoch normalization
    if strcmp(keys.normalization_field,'AN') || strcmp(keys.normalization_field,'TU')
        continue
    end
    
    %% PSTH calculation
    ut=ph_LR_to_CI(keys,pop,ut); % Convert to ipsi/contra? not sure if necessary here
    for t=1:numel(UC.type)
        typ=UC.type(t);
        eff=UC.type_effector(UC.type_effector(:,1)==typ,2);
        keys=ph_get_epoch_keys(keys,typ,eff(1),numel(eff)>1); %taking eff(1) is admittedly sloppy here
        ut_typ=[ut.type]==typ;
        
        %% preferred and unpreferred location (taken from instructed trials only! (?) Not necessary though
        if any(ut_typ)
            FR_for_pref=NaN(size(UC.position,1),1);
            for p=1:size(UC.position,1)
                ut_hemi=all(abs(bsxfun(@minus,vertcat(ut.position),UC.position(p,:)))<keys.cal.precision_tar,2) & [ut.choice]'==0; %% tr_pos?
                FR_for_pref(p)=nanmean((PF_epoch_FR(ut_hemi & ut_typ')-baseline(ut_hemi & ut_typ'))./norm_factor(ut_hemi & ut_typ'));
            end
            [~,pref_idx]=max(FR_for_pref);
            switch K.unpref_def
                case 'minimum'
                    [~,unpref_idx]=min(FR_for_pref);
                case 'horizontally_opposite'
                    [~,unpref_idx]=min(abs(nanmean([UC.position(:,1)+UC.position(pref_idx,1) UC.position(:,2)-UC.position(pref_idx,2)],2)));
                case 'diagonally_opposite'
                    [~,unpref_idx]=min(abs(nanmean([UC.position(:,1)+UC.position(pref_idx,1) UC.position(:,2)+UC.position(pref_idx,2)],2)));
            end
            
            pref_valid(t,u)=true;
            for ch=UC.choice
                pref_valid(t,u)=pref_valid(t,u) && pref_idx~=unpref_idx && ...
                    sum(all(abs(bsxfun(@minus,vertcat(ut(ut_typ).position),UC.position(pref_idx,:)))<keys.cal.precision_tar,2)   & [ut(ut_typ).choice]'==ch) >=keys.cal.min_trials_per_condition && ...
                    sum(all(abs(bsxfun(@minus,vertcat(ut(ut_typ).position),UC.position(unpref_idx,:)))<keys.cal.precision_tar,2) & [ut(ut_typ).choice]'==ch) >=keys.cal.min_trials_per_condition;
            end
        end
        
        %% average PSTH per unit
        for c=1:size(conditions_out,1)
            eff=conditions_out(c,strcmp(CP_out,'effector'));
            clear trpar
            
            for par=1:numel(CP_out)
                fn=CP_out{par};
                trpar(par,:)=[ut.(fn)]==conditions_out(c,par);
            end
            trpar(end+1,:)=ut_typ; %% comment this line to see where it leads??
            ut_con=all(trpar,1);
            %per_epoch=get_per_epoch(pop.trial(tr_con));% take already normalized values
            per_epoch=pop.epochs_per_type{typ}(ut_con(ut_typ),:);
            for ep=1:size(per_epoch,2)
                for p=1:size(UC.position,1)
                    ut_pos=all(abs(bsxfun(@minus,vertcat(ut.position),UC.position(p,:)))<keys.cal.precision_tar,2) & ut_con' & [ut.choice]'==0;
                    condition(t,c).per_position(p).epoch(ep).unit(u).FR=nanmean(per_epoch(ut_pos(ut_typ),ep))-nanmean(baseline(ut_pos & ut_typ'));
                end
            end
            if strcmp(keys.normalization_field,'PR')
                continue
            end
            
            for w=1:size(keys.PSTH_WINDOWS,1)
                for f=1:numel(UC.hemifield) %hemifield
                    ut_hemi=[ut.hemifield]==UC.hemifield(f);
                    switch keys.normalization_field
                        case {'PO','RE'}
                            ix = ut_con & ut_hemi & ut_typ;
                            if strcmp(keys.normalization_field,'RE')
                                [condition(t,c).per_hemifield(f).window(w).unit(u).average_spike_density,~,...
                                    condition(t,c).per_hemifield(f).window(w).unit(u).VAR_spike_density]=...
                                    ph_spike_density(pop.trial(ix),ut(ix),w,keys,baseline(ix),norm_factor(ix));
                            else
                                [condition(t,c).per_hemifield(f).window(w).unit(u).average_spike_density]=...
                                    ph_spike_density(pop.trial(ix),ut(ix),w,keys,baseline(ix),norm_factor(ix));
                            end
                            condition(t,c).per_hemifield(f).effector=eff;
                            condition(t,c).per_hemifield(f).window(w).unit(u).sac_lat=nanmean([ut(ix).sac_lat]);
                            condition(t,c).per_hemifield(f).window(w).unit(u).rea_lat=nanmean([ut(ix).rea_lat]);
                            condition(t,c).per_hemifield(f).window(w).unit(u).sac_sem=sterr([ut(ix).sac_lat]);
                            condition(t,c).per_hemifield(f).window(w).unit(u).rea_sem=sterr([ut(ix).rea_lat]);
                            
                            %% not sure if this makes much sense or if its even correct like this
                            ix=ix & ~isnan(PF_epoch_FR)&~isnan(baseline);
                            if sum(ix)>1 && any(ttest(PF_epoch_FR(ix), baseline(ix))==1) %% what what what why ttest here
                                condition(t,c).per_hemifield(f).sign.unit(u)=sign(mean(PF_epoch_FR(ix) - baseline(ix)));
                            else
                                condition(t,c).per_hemifield(f).sign.unit(u)=0;
                            end
                            
                        case 'ON'
                            %% alternative per trial
                            ix = ut_con & ut_hemi & ut_typ;
                            trials_for_SD=find(ix);
                            ep_hemi=find(ut_con(ut_typ) & ut_hemi(ut_typ));
                            temp_window=keys.PSTH_WINDOWS(w,:);
                            keys.PSTH_WINDOWS{w,3}=keys.PSTH_WINDOWS{w,3}-(keys.n_consecutive_bins_significant-1)*keys.PSTH_binwidth; % we need a few more so that first bins are not unsignificant by definition
                            keys.PSTH_WINDOWS{w,4}=keys.PSTH_WINDOWS{w,4}+(keys.n_consecutive_bins_significant-1)*keys.PSTH_binwidth;
                            condition_per_trial(t,c).per_hemifield(f).unit(u).epoch_FRs=...
                                reshape(per_epoch(ep_hemi,:),size(per_epoch(ep_hemi,:)))./repmat(norm_factor(ix)',1,size(per_epoch,2));
                            for tr=1:numel(trials_for_SD)
                                tr_tmp=trials_for_SD(tr);
                                condition_per_trial(t,c).per_hemifield(f).window(w).unit(u).average_spike_density(tr,:)=...
                                    ph_spike_density(pop.trial(tr_tmp),ut(tr_tmp),w,keys,baseline(tr_tmp),norm_factor(tr_tmp));
                            end
                            if sum(ix)==0
                                condition_per_trial(t,c).per_hemifield(f).window(w).unit(u).average_spike_density(1,:)=...
                                    NaN(size(ph_spike_density(pop.trial(1),ut(1),w,keys,baseline(1),norm_factor(1))));
                                condition_per_trial(t,c).per_hemifield(f).unit(u).epoch_FRs=NaN(1,size(keys.EPOCHS,1));
                            end
                            keys.PSTH_WINDOWS(w,:)=temp_window;
                            
                        case 'RT'
                            trials_for_SD=find(ut_acc & ut_hemi);
                            condition_per_trial(t,c).per_hemifield(f).unit(u).epoch_FRs=...
                                reshape([per_epoch(ut_hemi(ut_acc),:).FR],size(per_epoch(ut_hemi(ut_acc),:)))./repmat(norm_factor(trials_for_SD)',1,size(per_epoch,2));
                            condition_per_trial(t,c).per_hemifield(f).unit(u).sac_lat=[pop.trial(trials_for_SD).sac_lat];
                            condition_per_trial(t,c).per_hemifield(f).unit(u).rea_lat=[pop.trial(trials_for_SD).rea_lat];
                    end
                    
                    if strcmp(keys.normalization_field,'RE') %% bootstrapping PSTHs(!?)
                        n_iterations=100;
                        for it=1:n_iterations
                            ix = find(ut_acc & ut_hemi);
                            ix=randsample(ix,round(numel(ix)/10*8));
                            [condition(t,c).per_hemifield(f).window(w).unit(u).bootstrapped(it,:)]=...
                                ph_spike_density(pop.trial(ix),ut(ix),w,keys,baseline(ix),norm_factor(ix));
                        end
                    end
                end
                
                if strcmp(keys.normalization_field,'PO') %% && plot_position??
                    for p=1:size(UC.position,1)
                        ix=all(abs(bsxfun(@minus,vertcat(ut.position),UC.position(p,:)))<keys.cal.precision_tar,2)';
                        
                        condition(t,c).per_position(p).window(w).unit(u).average_spike_density= ph_spike_density(pop.trial(ix),ut(ix),w,keys,baseline(ix),norm_factor(ix));
                        condition(t,c).per_position(p).position=UC.position(p,:);
                        condition(t,c).per_position(p).fixation=1;
                        condition(t,c).per_position(p).effector=eff;
                        
                        if p==unpref_idx
                            condition(t,c).per_preference(1).window(w).unit(u).average_spike_density=...
                                ph_spike_density(pop.trial(ix),ut(ix),w,keys,baseline(ix),norm_factor(ix));
                            condition(t,c).per_preference(1).effector=eff;
                        end
                        if p==pref_idx
                            condition(t,c).per_preference(2).window(w).unit(u).average_spike_density=...
                                ph_spike_density(pop.trial(ix),ut(ix),w,keys,baseline(ix),norm_factor(ix));
                            condition(t,c).per_preference(2).effector=eff;
                        end
                        
                        ix=ix &  ~isnan(PF_epoch_FR)&~isnan(baseline);
                        % problem here if tr_con is empty (which it can be, because we are not requiring every condition to be valid for all units
                        if sum(ix)>1 && any(ttest(PF_epoch_FR(ix), baseline(ix))==1)
                            condition(t,c).per_position(p).sign.unit(u)=sign(mean(PF_epoch_FR(ix) - baseline(ix)));
                        else
                            condition(t,c).per_position(p).sign.unit(u)=0;
                        end
                    end
                end
            end
            
            zin_per_pos=NaN(size(UC.position,1),1);
            
            %% Always compute firing rates per position if theres at least one trial for htis condition (?)
            if sum(ut_acc) == 0
                continue
            end
            
            for p=1:size(UC.position,1) %for FR plot only
                tr=all(abs(bsxfun(@minus,vertcat(ut.position),UC.position(p,:)))<keys.cal.precision_tar,2);
                tr=tr(ut_typ) & ut_con(ut_typ)';
                if sum(tr)==0; continue; end
                zin_per_pos(p)=nanmean(vertcat(per_epoch(tr,RF_epoch{t})));
            end
            RF_positions=vertcat(ut(ut_typ).position);
            zin=vertcat(per_epoch(:,RF_epoch{t}));
            FR_tmp=struct('FR',num2cell(zin_per_pos),'x',num2cell(UC.position(:,1)),'y',num2cell(UC.position(:,2)));
            condition(t,c).fitting.unit(u).positions =FR_tmp;
            
            %% gaussian response fields (needs cleanup)
            if strcmp(keys.normalization_field,'RF') && ~strcmp(K.epoch_RF,'none')
                fitsettings=K.fitsettings;
                RF_tmp=ph_fit_target_positions_2D(RF_positions(:,1),RF_positions(:,2),zin,fitsettings);
                condition(t,c).fitting.unit(u).parameters=RF_tmp;
            end
            
        end
    end
end
end

