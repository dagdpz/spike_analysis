function tuning_per_unit_table=ph_ANOVAS(population,keys)
tuning_per_unit_table=keys.tuning_per_unit_table;
for unit=1:numel(population)
    anova_struct_current_unit=struct();
    if sum([population(unit).trial.accepted]==1 & [population(unit).trial.completed]==1)>0 % temporary? for pulv_oculomotor_dataset
        for a=1:numel(keys.position_and_plotting_arrangements) % arrangement defines poaitions, therefore also hemifield (which is part of conditions)
            keys.arrangement=keys.position_and_plotting_arrangements{a};
            %o=ph_arrange_positions_and_plots(keys,population(unit).trial(tr_considered)); % arrangement
            keys.labels.reach_hand=keys.labels.reach_handLR ;
            [UC, CM, lables]=ph_get_condition_matrix(population(unit).trial,keys);
            keys.normalization_field='AN';
            o=ph_condition_normalization(population(unit),keys,UC,CM); %condition wise normalization (also reduces conditions!???)
            
            for type=UC.type
                %% check carefully multicomp epochs !!
                keys.main_multicomp                     =keys.ANOVAS_PER_TYPE(type).main;
                keys.epoch_multicomp                    =keys.ANOVAS_PER_TYPE(type).epoch;
                keys.epoch_hemifield_multicomp          =keys.ANOVAS_PER_TYPE(type).hemifield;
                keys.epoch_hands_multicomp              =keys.ANOVAS_PER_TYPE(type).hands;
                keys.epoch_SxH_multicomp                =keys.ANOVAS_PER_TYPE(type).SxH;
                keys.epoch_position_multicomp           =keys.ANOVAS_PER_TYPE(type).positions;
                
                keys.anova_epochs=keys.ANOVAS_PER_TYPE(type); %% this here should eventually replace having to repeat the above
                for effector=UC.effector
                    keys=ph_get_epoch_keys(keys,type,effector,sum(UC.type_effector(:,1)==type)>1);
                    [~, condition_fieldname_part]=MPA_get_type_effector_name(type,effector);
                    tr_index= [o.trial.effector]==effector & [o.trial.type]==type ;
                    if sum(tr_index)==0;
                        fprintf('no trials for effector %.0f type %.0f \n',effector,type);continue;
                    end
                    o_e=o;
                    o_e.trial=o.trial(tr_index);
                    
                    trial_criterion=ph_get_minimum_trials(keys,o_e,CM,UC,lables);
                    [FR,epochs,idx,u_pos,u_fix]=ph_get_anova_inputs(o_e,keys);
                    [anova_struct]=n_way_anova(keys,epochs,FR,idx,u_pos,u_fix);
                    anova_struct_current_unit.([condition_fieldname_part '_'  keys.arrangement(1:3)])=anova_struct; clear anova_struct
                    FN_crit=fieldnames(trial_criterion);
                    for f=1:numel(FN_crit)
                        anova_struct_current_unit.([condition_fieldname_part '_'  keys.arrangement(1:3)]).(FN_crit{f})= trial_criterion.(FN_crit{f});
                    end
                end
                
                %% effector comparison !
                effectors_for_this_type=UC.type_effector(UC.type_effector(:,1)==type,2)';
                effectors_to_compare=combvec(effectors_for_this_type,effectors_for_this_type);
                effectors_to_compare=effectors_to_compare(:,effectors_to_compare(1,:)<effectors_to_compare(2,:));
                for comp=1:size(effectors_to_compare,2)
                    comp_eff=effectors_to_compare(:,comp);
                    keys=ph_get_epoch_keys(keys,type,comp_eff,1);
                    [~, condition_fieldname_part1]=MPA_get_type_effector_name(type,comp_eff(1));
                    [~, condition_fieldname_part2]=MPA_get_type_effector_name(type,comp_eff(2));
                    condition_fieldname_part=[condition_fieldname_part1 '_vs_' condition_fieldname_part2];
                    tr_index= ismember([o.trial.effector],comp_eff) & [o.trial.type] == type & ismember([o.trial.perturbation],keys.cal.perturbation_groups{1});
                    if sum(tr_index)==0;
                        fprintf('no trials for effectors %s type %.0f perturbation %s \n',mat2str(comp_eff'),type,mat2str(keys.cal.perturbation_groups{1}));continue;
                    end
                    o_e=o;
                    o_e.trial=o.trial(tr_index);
                    keys.TTlabels.effector={condition_fieldname_part1 '-' condition_fieldname_part2};
                    [FR,epochs,idx,u_pos,u_fix]=ph_get_anova_inputs(o_e,keys);
                    [anova_struct]=effector_comparison_anova(keys,epochs,FR,idx,u_pos,u_fix);
                    anova_struct_current_unit.([condition_fieldname_part '_'  keys.arrangement(1:3)])=anova_struct; clear anova_struct
                end
            end
        end
    end
    rows_to_update=find(ismember(tuning_per_unit_table(:,DAG_find_column_index(tuning_per_unit_table,'unit_ID')),{population(unit).unit_ID}));
    if isempty(rows_to_update)
        rows_to_update=size(tuning_per_unit_table,1)+1;
    end
    clear unit_table;
    inital_fieldnames={'unit_ID','monkey','target','perturbation_site','grid_x','grid_y','electrode_depth','FR','stability_rating','SNR_rating','Single_rating','waveform_width'};
    %inital_fieldnames={'unit_ID','monkey','target','grid_x','grid_y','electrode_depth','FR','stability_rating','SNR_rating','Single_rating','waveform_width'};
    unit_table(1,1:numel(inital_fieldnames))=inital_fieldnames;
    for fn=1:numel(inital_fieldnames)
        unit_table{2,fn}=population(unit).(inital_fieldnames{fn});
    end
    title_counter=numel(inital_fieldnames);
    FN=fieldnames(anova_struct_current_unit);
    for fn=1:numel(FN)
        FNsub=fieldnames(anova_struct_current_unit.(FN{fn}));
        for fnsub=1:numel(FNsub)
            title_counter=title_counter+1;
            unit_table{1,title_counter}=[FNsub{fnsub} '_' FN{fn}];
            unit_table{2,title_counter}=anova_struct_current_unit.(FN{fn}).(FNsub{fnsub});
        end
    end
%     tic
%     tuning_per_unit_table=DAG_update_mastertable_cell(tuning_per_unit_table,unit_table,rows_to_update);
%     toc

    
    titles_update=unit_table(1,:);
    titles_original_table=tuning_per_unit_table(1,:);
    [~, title_positions]=ismember(titles_update,titles_original_table);
    for n_update_column=1:numel(titles_update)
        XXX=title_positions(n_update_column);
        if XXX==0
            n_column=size(tuning_per_unit_table,2)+1;
            tuning_per_unit_table(1,n_column)=titles_update(n_update_column);
        else
            n_column=XXX;
        end
        tuning_per_unit_table(rows_to_update,n_column)=unit_table(2,n_update_column);
    end
end
end

function anova_struct=n_way_anova(keys,epochs,FR,idx,Positions,Fixations)
labels=keys.TTlabels;
labels.control_test={};
conditions=keys.TTconditions;
IN=conditions.choice{1};
CH=conditions.choice{2};
conditions.choice=conditions.choice([sum(idx.in)>0 sum(idx.ch)>0]);

[~, ~, epoch_idx]=unique(epochs);
idx_fieldnames=fieldnames(idx)';

% epochs for epoch comparison
epoch_multicomp=keys.epoch_multicomp;
idx_ep=ismember(epoch_multicomp(:,1),epochs') &  ismember(epoch_multicomp(:,2),epochs');
epoch_multicomp=epoch_multicomp(idx_ep,:);
%% Independently for Instructed and Choice !!
for ch=1:numel(conditions.choice)
    INCH=conditions.choice{ch};
    labels.INCH=conditions.choice{ch};
    %% perturbation here is defined as always the lowest accessible for control... (???)
    tr_ch=idx.(INCH) & idx.PT==0;
    for fn=idx_fieldnames
        idx.(['tr_' fn{:}])=      idx.(fn{:}) & repmat(tr_ch,1,size(idx.(fn{:}),2)) ;
    end
    
    %% number of trials for each condition to decide about type of (main) anova
    n_lr=sum(idx.tr_LS)>0 && sum(idx.tr_RS)>0;
    n_2hands=sum(idx.tr_LH_LS)>0 && sum(idx.tr_LH_RS)>0 && sum(idx.tr_RH_LS)>0 && sum(idx.tr_RH_RS)>0;
    if n_2hands % three-way anova epoch, space, hand
        keys.anova_varnames={'epoch' 'hemifield' 'hands'};
        tr=idx.tr_LH_LS | idx.tr_LH_RS | idx.tr_RH_LS | idx.tr_RH_RS;
        Par=[epoch_idx, idx.tr_RS, idx.tr_RH];
    elseif n_lr % two-way anova epoch, space
        keys.anova_varnames={'epoch' 'hemifield'};
        tr=idx.tr_LS | idx.tr_RS;
        Par=[epoch_idx, idx.tr_RS];
    else
        disp(['Not fitting dataset... no left or right ' INCH ' present?'])
        anova_struct.([INCH '_epoch_main'])=[];
        continue;
    end
        
    %% main effects and general interactions
    tr_main=idx.(INCH) & ismember(epochs,keys.main_multicomp) & idx.PT==0; %& ismember(epochs,epochs_for_multicomparison);
    [anova_out.p,anova_out.table,anova_out.stats,anova_out.terms] = anovan(FR(tr_main),Par(tr_main,:),'model','full','varnames',keys.anova_varnames,'display',keys.plot.anova_tables);
    anova_struct.([INCH '_epoch_main'])=anova_out.p(1)<0.05; %main effect on epoch!
    handindexes=[1,2,3];
    handindexes=handindexes(any(idx.tr_hands));

    %% Epoch X space/hand tuning (regardless of the other)
    for k=2:numel(keys.anova_varnames)
        multicomp_epochs=keys.(['epoch_' keys.anova_varnames{k} '_multicomp']);
        %multicomp_epochs=multicomp_epochs(ismember(multicomp_epochs,epochs')); %% to think about: add &  ismember(keys.epoch_multicomp(:,2),keys.EPOCHS(:,1));
        label=labels.(keys.anova_varnames{k});
        idx1m=tr_main & Par(:,k)==0 & ismember(epochs,multicomp_epochs);
        idx2m=tr_main & Par(:,k)==1 & ismember(epochs,multicomp_epochs);
        labelindex=(anova_out.p(k)<0.05) *sign(nanmean(FR(idx2m))-nanmean(FR(idx1m)))+2; labelindex(isnan(labelindex))=2;
        anova_struct.([INCH '_' keys.anova_varnames{k} '_main'])=label{labelindex}; %main effect with direction!
        anova_struct.([INCH '_Ex' upper(keys.anova_varnames{k}(1))])=anova_out.p(k+numel(keys.anova_varnames)-1)<0.05; %TRUE/FALSE labels??
    end

    %% space x hand anovas per epoch
    if ismember('hemifield',keys.anova_varnames) && ismember('hands',keys.anova_varnames)
        k=6;
        tr_SXH=tr_main & ismember(epochs,keys.epoch_SxH_multicomp);
        labelindexcr=(anova_out.p(k)<0.05) *sign(nanmean(FR(tr_SXH & idx.tr_CR))-nanmean(FR(tr_SXH & idx.tr_UC)))+2; labelindexcr(isnan(labelindexcr))=2;% 3 crossed>uncrossed ,1 crossed<uncrossed
        anova_struct.([INCH '_SxH'])=labels.CR{labelindexcr};
        anova_struct.([INCH '_ExSxH'])=anova_out.p(7)<0.05;
               
        for s=keys.epoch_SxH_multicomp'
            idxS= tr & ismember(epochs,s);
            if any(~isnan(FR(idxS))) && any(idx.tr_LS(idxS)) && any(idx.tr_RS(idxS)) && any(idx.tr_LH(idxS)) && any(idx.tr_RH(idxS))
                [anova_outs.p,anova_outs.table,anova_outs.stats,anova_outs.terms] = anovan(FR(idxS),double([idx.tr_RS(idxS),idx.tr_RH(idxS)]),'model','full','varnames',{'space';'hands'},'display',keys.plot.anova_tables);
                
                h=anova_outs.p(1)<0.05;
                DF=nanmean(FR(idx.tr_RS & idxS))-nanmean(FR(idx.tr_LS & idxS));
                labelindexsp=h*sign(DF)+2; labelindexsp(isnan(labelindexsp))=2;
                anova_struct.([INCH '_' s{:} '_SxH_hemifield'])= labels.hemifield{labelindexsp};
                anova_struct.([INCH '_' s{:} '_SxH_hemifield_PV'])= single(anova_outs.p(1)*sign(DF));
                anova_struct.([INCH '_' s{:} '_SxH_hemifield_DF'])= DF;
                anova_struct.([INCH '_' s{:} '_SxH_hemifield_IX'])= DF/(nanmean(FR(idx.tr_RS & idxS))+nanmean(FR(idx.tr_LS & idxS)));
                
                h=anova_outs.p(2)<0.05;
                DF=nanmean(FR(idx.tr_RH & idxS))-nanmean(FR(idx.tr_LH & idxS));
                labelindexha=h*sign(DF)+2; labelindexha(isnan(labelindexha))=2;
                anova_struct.([INCH '_' s{:} '_SxH_hands'])= labels.hands{labelindexha};
                anova_struct.([INCH '_' s{:} '_SxH_hands_PV'])= single(anova_outs.p(2)*sign(DF));
                anova_struct.([INCH '_' s{:} '_SxH_hands_DF'])= DF;
                anova_struct.([INCH '_' s{:} '_SxH_hands_IX'])= DF/(nanmean(FR(idx.tr_RH & idxS))+nanmean(FR(idx.tr_LH & idxS)));
                
                h=anova_outs.p(3)<0.05;
                DF=nanmean(FR(idx.tr_CR & idxS))-nanmean(FR(idx.tr_UC & idxS));
                labelindexcr=h*sign(DF)+2; labelindexcr(isnan(labelindexcr))=2;
                anova_struct.([INCH '_' s{:} '_SxH'])=labels.CR{labelindexcr};
                anova_struct.([INCH '_' s{:} '_SxH_PV'])= single(anova_outs.p(3)*sign(DF));
                anova_struct.([INCH '_' s{:} '_SxH_DF'])=DF;
                anova_struct.([INCH '_' s{:} '_SxH_IX'])=DF/(nanmean(FR(idx.tr_CR & idxS))+nanmean(FR(idx.tr_UC & idxS)));
                
            end
        end
    end    
    
    %% epoch tuning across all trials    
    labels.conditions={''};
    idx.conditions=idx.tr_all;
    anova_struct=do_epoch_stats(keys,anova_struct,FR,epochs,idx,labels,'epoch_multicomp','epoch');
    
    %% epoch tuning per hand & space
    labels.conditions={'AH_LS','AH_RS','LH_LS','LH_RS','RH_LS','RH_RS'};
    idx.conditions=[idx.tr_LS idx.tr_RS idx.tr_LH_LS idx.tr_LH_RS idx.tr_RH_LS idx.tr_RH_RS];
    anova_struct=do_epoch_stats(keys,anova_struct,FR,epochs,idx,labels,'epoch_multicomp','epoch');
    
    %% epoch tuning per hand (taking tuning for both HFs into account, 'bi' if opposing direction for both HFs)
    
    for row=1:size(epoch_multicomp,1)
        s=epoch_multicomp(row,2);
        b=epoch_multicomp(row,1);
        for hn=handindexes
            LHRH=conditions.hands{hn};
            % select label dependend on both sides - make sure they are both present?
            combined_tuning=[anova_struct.([INCH '_' LHRH '_LS_' s{:} '_epoch'])    anova_struct.([INCH '_' LHRH '_RS_' s{:} '_epoch'])];
            combined_DF=    [anova_struct.([INCH '_' LHRH '_LS_' s{:} '_epoch_DF']) anova_struct.([INCH '_' LHRH '_RS_' s{:} '_epoch_DF'])];
            switch combined_tuning
                case {'-su','su-','susu'};   labelindex=1;
                case {'--'};                 labelindex=2;
                case {'-en','en-','enen'};   labelindex=3;
                case {'suen','ensu'};        labelindex=4;
            end
            anova_struct.([INCH '_' LHRH  '_' s{:} '_epoch'])=labels.epoch{labelindex};
            [~,DFidx]=max(abs(combined_DF)); %this does not necessarily make sense to be honest
            DF=combined_DF(DFidx);
            anova_struct.([INCH '_' LHRH  '_' s{:} '_epoch_DF'])=DF; %%% to think about !!!!
            idxbl=tr & ismember(epochs,b) & idx.sides(:,DFidx) & idx.tr_hands(:,hn); % take respective baseline index for side!
            anova_struct.([INCH '_' LHRH  '_' s{:} '_epoch_SC'])=DF/nanmean(FR(idxbl)); %%% to think about !!!! now it corresponds to % signal change
            labelindex=-1*strcmp(combined_tuning,'susu')+strcmp(combined_tuning,'enen')+2;
            anova_struct.([INCH '_' LHRH  '_' s{:} '_epoch_bilateral'])=labels.epoch{labelindex};
        end
    end
      
    %% hand and space tuning per epoch (across all trials)
    labels.conditions={''};
    idx.conditions=idx.tr_all;
    anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,'epoch_hemifield_multicomp','tr_LS','tr_RS','hemifield');
    anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,'epoch_hands_multicomp','tr_LH','tr_RH','hands');
    
    %% hemifield tuning per hand
    labels.conditions=conditions.hands;
    idx.conditions=idx.tr_hands;
    anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,'epoch_SxH_multicomp','tr_LS','tr_RS','hemifield');
            
    %% hand tuning per hemifield
    labels.conditions=conditions.hemifield;
    idx.conditions=idx.sides;
    anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,'epoch_SxH_multicomp','tr_LH','tr_RH','hands');
           
    %%  Difficulty stuff    
    labels.conditions=conditions.hemifield;
    idx.conditions=idx.tr_sides;
    anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,'epoch_hemifield_multicomp','tr_Diff0','tr_Diff1','Difficulty_Easy');
    anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,'epoch_hemifield_multicomp','tr_Diff0','tr_Diff2','Difficulty_Diff');
    anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,'epoch_hemifield_multicomp','tr_nonDistr1','tr_TT1HF','SpatialComp_1HFTar');
    anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,'epoch_hemifield_multicomp','tr_nonDistr1','tr_TT2HF','SpatialComp_2HFTar');
        
    %%  single trials per space - this did not do what it's supposed to do !
    labels.conditions={''};
    idx.conditions=idx.tr_all;
    anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,'epoch_hemifield_multicomp','tr_SglL','tr_SglR','SglTar_Suc');
    
    %% old version to compare
    %     multicomp_epochs=keys.(['epoch_hemifield_multicomp']);
    %     multicomp_epochs=multicomp_epochs(ismember(multicomp_epochs,epochs'));
    %     label={'SC','-','SI'}; %higher FR for CS , higher FR for IS
    %     idx1= idx.nonDistr1 &  idx.suc1 ;
    %     idx2= idx.nonDistr1 &  idx.suc1 ;
    %     % idx1= sum([(idx.nonDistr1 &  idx.suc1) , (idx.Distr2 &  idx.suc0)],2) ;
    %     % idx2= sum([(idx.nonDistr1 &  idx.suc1) , (idx.Distr2 &  idx.suc0)],2) ;
    %
    %     %t-test
    %     for sideindex=1:2 %left * right
    %         for s=multicomp_epochs(:)'
    %             idxS=ismember(epochs,s)  & idx.tr_sides(:,sideindex);
    %             h=do_stats(FR(idx1 & idxS),FR(idx2 & idxS),keys,0); %not paired
    %             DF=nanmean(FR(idx2 & idxS))-nanmean(FR(idx1 & idxS));
    %             labelindex=h*sign(DF)+2; labelindex(isnan(labelindex))=2;
    %             anova_struct.([INCH '_' LSRSnamepart{sideindex} '_' s{:} '_' 'SglTar_Suc' ])=label{labelindex};
    %             anova_struct.([INCH '_' LSRSnamepart{sideindex} '_' s{:} '_' 'SglTar_Suc' '_DF'])=DF;
    %             anova_struct.([INCH '_' LSRSnamepart{sideindex} '_' s{:} '_' 'SglTar_Suc' '_IX'])=DF/(nanmean(FR(idx2 & idxS))+ nanmean(FR(idx1 & idxS)));
    %         end
    %     end
    
    %% Position anova (per hand) + comparison CH vs IN
    epochs_for_position_comparison=keys.epoch_position_multicomp(ismember(keys.epoch_position_multicomp,keys.EPOCHS(:,1)));
    multicomp_epochs=epochs_for_position_comparison(ismember(epochs_for_position_comparison,epochs));
    
    u_hnd                       =find(any(idx.tr_hands,1));
    for hn=u_hnd
        LHRH=conditions.hands{hn};
        tr=idx.(INCH) & idx.tr_hands(:,hn) & ismember(epochs,epochs_for_position_comparison);
        if size(Fixations,1)==1
            varnames={'State','position'};
            [anova_out.p,anova_out.table,anova_out.stats,anova_out.terms] =...
                anovan(FR(tr),{epochs(tr),idx.pos(tr)},'model','full','varnames',varnames,'display',keys.plot.anova_tables);
        else
            varnames={'State','position','fixation'};
            [anova_out.p,anova_out.table,anova_out.stats,anova_out.terms] =...
                anovan(FR(tr),{epochs(tr),idx.pos(tr),idx.fix(tr)},'model','full','varnames',varnames,'display',keys.plot.anova_tables);
            anova_struct.([INCH '_' LHRH '_fixation_main'])=anova_out.p(3)<0.05; %main effect on epoch!
        end
        anova_struct.([INCH '_' LHRH '_epoch_main'])=anova_out.p(1)<0.05; %main effect on epoch!
        anova_struct.([INCH '_' LHRH '_position_main'])=anova_out.p(2)<0.05; %main effect on epoch!
        
        
        for s=multicomp_epochs(:)'
            tr=idx.(INCH) & idx.tr_hands(:,hn) & ismember(epochs,s);  %this is not ideal, because tr is used in a different way here
            b=epoch_multicomp(ismember(epoch_multicomp(:,2),s),1);
            
            if any(~isnan(FR(tr)))
                [anova_outs.p] = anovan(FR(tr),[idx.pos(tr),idx.fix(tr)],'model','full','varnames',{'Positions','Fixations'},'display',keys.plot.anova_tables);
                if any(isnan(anova_outs.p)) %this is the case for target position by location in pulvinar eye gaze project AND for fixation only
                    % take only positions which have combinations of different fixations to compute interaction
                    pos_fix=[idx.pos(tr) idx.fix(tr)];
                    u_pos_fix=unique(pos_fix,'rows');
                    u_pos_fix1_valid=find(hist(u_pos_fix(:,1),unique(idx.pos(tr)))>1);
                    u_pos_fix2_valid=find(hist(u_pos_fix(:,2),unique(idx.fix(tr)))>1);
                    u_pos_fix=u_pos_fix(ismember(u_pos_fix(:,1),u_pos_fix1_valid') & ismember(u_pos_fix(:,2),u_pos_fix2_valid'),:);
                    pos_fix_valid=ismember(pos_fix,u_pos_fix,'rows');
                    FR_temp=FR(tr);
                    if any(pos_fix_valid)
                        [ptemp] = anovan(FR_temp(pos_fix_valid),pos_fix(pos_fix_valid,:),'model','full','varnames',{'Positions','Fixations'},'display',keys.plot.anova_tables);
                    else
                        ptemp(2)=NaN;
                        ptemp(3)=NaN;
                    end
                    
                    [anova_outs.p] = anovan(FR(tr),idx.pos(tr),'model','full','varnames',{'Positions'},'display',keys.plot.anova_tables);
                    anova_outs.p(2)=NaN;
                    anova_outs.p(3)=ptemp(2)>0.05 || ptemp(3)>0.05; % this is not ideal, interaction here only significant if there is a main effect as well....
                end
                
                h1=anova_outs.p(1)<0.05;
                h2=anova_outs.p(2)<0.05;
                h3=anova_outs.p(3)<0.05;
                
                anova_struct.([INCH '_' LHRH '_' s{:} '_position'])     =labels.true{h1+1};
                anova_struct.([INCH '_' LHRH '_' s{:} '_fixation'])     =labels.true{h2+1};
                anova_struct.([INCH '_' LHRH '_' s{:} '_PxF'])          =labels.true{h3+1};
                anova_struct.([INCH '_' LHRH '_' s{:} '_position_PV'])  =single(anova_outs.p(1));
                anova_struct.([INCH '_' LHRH '_' s{:} '_fixation_PV'])  =single(anova_outs.p(2));
                anova_struct.([INCH '_' LHRH '_' s{:} '_PxF_PV'])       =single(anova_outs.p(3));
                
                
                pecc = anovan(FR(tr),[idx.ecc(tr),idx.ang(tr)],'model','full','varnames',{'Eccentricity','Angle'},'display',keys.plot.anova_tables);
                
                h1=pecc(1)<0.05;
                h2=pecc(2)<0.05;
                h3=pecc(3)<0.05;
                anova_struct.([INCH '_' LHRH '_' s{:} '_distance'])     =labels.true{h1+1};
                anova_struct.([INCH '_' LHRH '_' s{:} '_angle'])        =labels.true{h2+1};
                anova_struct.([INCH '_' LHRH '_' s{:} '_DxA'])          =labels.true{h3+1};
                anova_struct.([INCH '_' LHRH '_' s{:} '_distance_PV'])  =single(pecc(1));
                anova_struct.([INCH '_' LHRH '_' s{:} '_angle_PV'])     =single(pecc(2));
                anova_struct.([INCH '_' LHRH '_' s{:} '_DxA_PV'])       =single(pecc(3));
                
                % x and y separately - didnt add p_values here yet!
                pxy = anovan(FR(tr),[idx.pos_x(tr),idx.pos_y(tr)],'model','full','varnames',{'x','y'},'display',keys.plot.anova_tables);
                h1=pxy(1)<0.05;
                h2=pxy(2)<0.05;
                h3=pxy(3)<0.05;
                anova_struct.([INCH '_' LHRH '_' s{:} '_positionx'])=labels.true{h1+1};
                anova_struct.([INCH '_' LHRH '_' s{:} '_positiony'])=labels.true{h2+1};
                anova_struct.([INCH '_' LHRH '_' s{:} '_positionxy'])=labels.true{h3+1};
                for xory={'x','y'}
                    xorytag=xory{:};
                    switch xorytag
                        case 'x'
                            labelxy=labels.hemifield;
                            hxy=h1;
                        case 'y'
                            labelxy=labels.UD;
                            hxy=h2;
                    end
                    if hxy
                        clear meanFR_per_position ttestperpos directionality_per_position
                        u_xy_pos=unique(idx.(['pos_' xorytag])(tr));
                        for f=u_xy_pos'
                            meanFR_per_position(f)=nanmean(FR(tr & idx.(['pos_' xorytag])==f));
                            ttestperpos(f)=do_stats(FR(tr & idx.(['pos_' xorytag])==f),FR(tr & idx.(['pos_' xorytag])==(mod(f,numel(u_xy_pos))+1)),keys,0);
                            directionality_per_position(f)=sign(nanmean(FR(tr & idx.(['pos_' xorytag])==(mod(f,numel(u_xy_pos))+1)))-nanmean(FR(tr & idx.(['pos_' xorytag])==f)));
                        end
                        [~,pref_idx]=max(meanFR_per_position);
                        directionality_per_position(end)=directionality_per_position(end)*-1; %need to invert, because we are doing the full circle i.e. comparing last to first
                        if ttestperpos(end)==1 && all(diff(directionality_per_position(ttestperpos==1))==0) % new definition of monotonic
                            anova_struct.([INCH '_' LHRH '_' s{:} '_gaze_modulation_' xorytag])='monotonous';
                            anova_struct.([INCH '_' LHRH '_' s{:} '_gaze_pref_' xorytag])=labelxy{(directionality_per_position(end)==1)*2+1};
                        else
                            anova_struct.([INCH '_' LHRH '_' s{:} '_gaze_modulation_' xorytag])='nonmonotonous';
                            if pref_idx==1 || pref_idx==numel(u_xy_pos)
                                anova_struct.([INCH '_' LHRH '_' s{:} '_gaze_pref_' xorytag])='PE';
                            else
                                anova_struct.([INCH '_' LHRH '_' s{:} '_gaze_pref_' xorytag])='CE';
                            end
                        end
                    else % not even nonmonotoneous
                        anova_struct.([INCH '_' LHRH '_' s{:} '_gaze_modulation_' xorytag])='-';
                        anova_struct.([INCH '_' LHRH '_' s{:} '_gaze_pref_' xorytag])='-';
                    end
                end
                
                %% choice part for position with strongest response
                if ch==2  && ~isempty(anova_struct.in_epoch_main) && isfield(anova_struct,['in_' LHRH '_' s{:} '_epoch_DF'])%any(idx.in) %% there has to be choice and instructed
                    clear Average_FR_per_position_IN Average_FR_per_position_CH
                    for p=1:max(idx.pos)
                        Average_FR_per_position_IN(p)=nanmean(FR(idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.pos==p));
                        Average_FR_per_position_CH(p)=nanmean(FR(idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.pos==p));
                    end
                    % preference (highest or lowest) defined by over all
                    % increase or decrease in this epoch
                    if false % anova_struct.(['in_' LHRH '_' s{:} '_epoch_DF']) <0
                        invertsign=-1;
                        takemin=1;
                    else
                        invertsign=1;
                        takemin=0;
                    end
                    
                    if takemin
                        [~,RF_position_index_IN]=min(Average_FR_per_position_IN);
                        [~,RF_position_index_CH]=min(Average_FR_per_position_CH);
                    else
                        [~,RF_position_index_IN]=max(Average_FR_per_position_IN);
                        [~,RF_position_index_CH]=max(Average_FR_per_position_CH);
                    end

                    % bootstrapping choice preference: for every iteration,
                    % 50% of trials to each hemifield are taken to estimate
                    % preferred hemifield, and the difference is computed
                    % based on the remaining 50%


                        %% preference based on choice trials (??)
%                         RF_CH_RS_IDX=idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.RS;
%                         RF_CH_LS_IDX=idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.LS;           
%                         RF_in_hemifield_index_CH =invertsign*sign(nanmean(FR(RF_CH_RS_IDX)) - nanmean(FR(RF_CH_LS_IDX))) +2;
%                         RF_out_hemifield_index_CH=invertsign*sign(nanmean(FR(RF_CH_RS_IDX)) - nanmean(FR(RF_CH_LS_IDX)))*-1 +2;

%                         RF_in_hemifield_index_CH(isnan(RF_in_hemifield_index_CH))=2;
%                         RF_out_hemifield_index_CH(isnan(RF_out_hemifield_index_CH))=2;
%                         idx_IN_prefHI_ch  =idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_in_hemifield_index_CH);
%                         idx_IN_prefHO_ch  =idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_out_hemifield_index_CH);
%                         idx_CH_prefHI_ch  =idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_in_hemifield_index_CH);
%                         idx_CH_prefHO_ch  =idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_out_hemifield_index_CH);
                        
                    
                    n_boots=100;
                    RF_IN_RS_IDX=idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.RS;
                    RF_IN_LS_IDX=idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.LS;
                    RF_IN_RS_IDX_bl=idx.in & idx.hands(:,hn) & ismember(epochs,b) & idx.RS;
                    RF_IN_LS_IDX_bl=idx.in & idx.hands(:,hn) & ismember(epochs,b) & idx.LS;
                    FR_LS_in=FR(RF_IN_LS_IDX);
                    FR_RS_in=FR(RF_IN_RS_IDX);
                    FR_LS_in_bl=FR(RF_IN_LS_IDX_bl);
                    FR_RS_in_bl=FR(RF_IN_RS_IDX_bl);
                    [~,ix_L]=sort(rand(n_boots,numel(FR_LS_in)),2);
                    [~,ix_R]=sort(rand(n_boots,numel(FR_RS_in)),2);
                        
                    for boots=1:n_boots
                        FR_R_L=nanmean(FR_RS_in(ix_R(boots,1:round(size(ix_R,2)/2)))) - nanmean(FR_LS_in(ix_L(boots,1:round(size(ix_L,2)/2))));
                        RF_in_hemifield_index_IN =invertsign*sign(FR_R_L) +2;
                        RF_out_hemifield_index_IN=invertsign*sign(FR_R_L)*-1 +2;
                        RF_in_hemifield_index_IN(isnan(RF_in_hemifield_index_IN))=2;
                        RF_out_hemifield_index_IN(isnan(RF_out_hemifield_index_IN))=2;
                        
                        FR_R_other50=nanmean(FR_RS_in(ix_R(boots,round(size(ix_R,2)/2)+1:end)));
                        FR_L_other50=nanmean(FR_LS_in(ix_L(boots,round(size(ix_L,2)/2)+1:end)));
                        FR_R_other50_bl=nanmean(FR_RS_in_bl(ix_R(boots,round(size(ix_R,2)/2)+1:end)));
                        FR_L_other50_bl=nanmean(FR_LS_in_bl(ix_L(boots,round(size(ix_L,2)/2)+1:end)));
                        
                        if RF_in_hemifield_index_IN ==1 %% left preferred                            
                            bootstrapped.in_prefHI_FR(boots)=FR_L_other50;
                            bootstrapped.in_prefHO_FR(boots)=FR_R_other50;  
                            bootstrapped.in_prefHI_FR_bl(boots)=FR_L_other50_bl;
                            bootstrapped.in_prefHO_FR_bl(boots)=FR_R_other50_bl;  
                        elseif RF_in_hemifield_index_IN ==3 %% right preferred                            
                            bootstrapped.in_prefHI_FR(boots)=FR_R_other50;
                            bootstrapped.in_prefHO_FR(boots)=FR_L_other50;  
                            bootstrapped.in_prefHI_FR_bl(boots)=FR_R_other50_bl;
                            bootstrapped.in_prefHO_FR_bl(boots)=FR_L_other50_bl;  
                        else %% no preference                            
                            bootstrapped.in_prefHI_FR(boots)=NaN;
                            bootstrapped.in_prefHO_FR(boots)=NaN;        
                            bootstrapped.in_prefHI_FR_bl(boots)=NaN;
                            bootstrapped.in_prefHO_FR_bl(boots)=NaN;                        
                        end
                        %                         idx_IN_prefHI_in  =idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_in_hemifield_index_IN)  & RF_TEST;
                        %                         idx_IN_prefHO_in  =idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_out_hemifield_index_IN) & RF_TEST;
                        idx_CH_prefHI_in  =idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_in_hemifield_index_IN);
                        idx_CH_prefHO_in  =idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_out_hemifield_index_IN);
                        idx_CH_prefHI_in_bl  =idx.ch & idx.hands(:,hn) & ismember(epochs,b) & idx.LR(:,RF_in_hemifield_index_IN);
                        idx_CH_prefHO_in_bl  =idx.ch & idx.hands(:,hn) & ismember(epochs,b) & idx.LR(:,RF_out_hemifield_index_IN);
                        
                        bootstrapped.ch_prefHI_FR(boots)=nanmean(FR(idx_CH_prefHI_in));
                        bootstrapped.ch_prefHO_FR(boots)=nanmean(FR(idx_CH_prefHO_in));     
                        bootstrapped.ch_prefHI_FR_bl(boots)=nanmean(FR(idx_CH_prefHI_in_bl));
                        bootstrapped.ch_prefHO_FR_bl(boots)=nanmean(FR(idx_CH_prefHO_in_bl));                      
                    end

%                     for boots=1:100
%                         %% preference based on instructed trials
%                         RF_IN_RS_IDX=idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.RS;
%                         RF_IN_LS_IDX=idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.LS;
%                         RF_IN_RS_IDX(randsample(find(RF_IN_RS_IDX),round(sum(RF_IN_RS_IDX)/2)))=false;  % taking only 50% of trials for preference estimation
%                         RF_IN_LS_IDX(randsample(find(RF_IN_LS_IDX),round(sum(RF_IN_LS_IDX)/2)))=false;
%                         RF_TEST=~RF_IN_RS_IDX & ~RF_IN_LS_IDX;
%                         RF_in_hemifield_index_IN =invertsign*sign(nanmean(FR(RF_IN_RS_IDX)) - nanmean(FR(RF_IN_LS_IDX))) +2;
%                         RF_out_hemifield_index_IN=invertsign*sign(nanmean(FR(RF_IN_RS_IDX)) - nanmean(FR(RF_IN_LS_IDX)))*-1 +2;
%                         
%                         %% preference based on choice trials (??)
%                         RF_CH_RS_IDX=idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.RS;
%                         RF_CH_LS_IDX=idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.LS;           
%                         RF_in_hemifield_index_CH =invertsign*sign(nanmean(FR(RF_CH_RS_IDX)) - nanmean(FR(RF_CH_LS_IDX))) +2;
%                         RF_out_hemifield_index_CH=invertsign*sign(nanmean(FR(RF_CH_RS_IDX)) - nanmean(FR(RF_CH_LS_IDX)))*-1 +2;
%                         
%                         RF_in_hemifield_index_IN(isnan(RF_in_hemifield_index_IN))=2;
%                         RF_out_hemifield_index_IN(isnan(RF_out_hemifield_index_IN))=2;
%                         RF_in_hemifield_index_CH(isnan(RF_in_hemifield_index_CH))=2;
%                         RF_out_hemifield_index_CH(isnan(RF_out_hemifield_index_CH))=2;
%                         
%                         %% full hemifield
%                         idx_IN_prefHI_in  =idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_in_hemifield_index_IN)  & RF_TEST;
%                         idx_IN_prefHO_in  =idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_out_hemifield_index_IN) & RF_TEST;
%                         idx_CH_prefHI_in  =idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_in_hemifield_index_IN)  & RF_TEST;
%                         idx_CH_prefHO_in  =idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_out_hemifield_index_IN) & RF_TEST;
%                         
%                         idx_IN_prefHI_ch  =idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_in_hemifield_index_CH);
%                         idx_IN_prefHO_ch  =idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_out_hemifield_index_CH);
%                         idx_CH_prefHI_ch  =idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_in_hemifield_index_CH);
%                         idx_CH_prefHO_ch  =idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_out_hemifield_index_CH);
%                         
%                         bootstrapped.in_prefHI_FR(boots)=nanmean(FR(idx_IN_prefHI_in));
%                         bootstrapped.ch_prefHI_FR(boots)=nanmean(FR(idx_CH_prefHI_in));
%                         bootstrapped.in_prefHO_FR(boots)=nanmean(FR(idx_IN_prefHO_in));
%                         bootstrapped.ch_prefHO_FR(boots)=nanmean(FR(idx_CH_prefHO_in));
%                     end
                    
                    %% baselines per hemifield???
%                     idxb_IN_prefHI_in  =idx.in & idx.hands(:,hn) & ismember(epochs,b) & idx.LR(:,RF_in_hemifield_index_IN);
%                     idxb_IN_prefHO_in  =idx.in & idx.hands(:,hn) & ismember(epochs,b) & idx.LR(:,RF_out_hemifield_index_IN);
%                     idxb_CH_prefHI_in  =idx.ch & idx.hands(:,hn) & ismember(epochs,b) & idx.LR(:,RF_in_hemifield_index_IN);
%                     idxb_CH_prefHO_in  =idx.ch & idx.hands(:,hn) & ismember(epochs,b) & idx.LR(:,RF_out_hemifield_index_IN);
                    
%                     idxb_IN_prefHI_ch  =idx.in & idx.hands(:,hn) & ismember(epochs,b) & idx.LR(:,RF_in_hemifield_index_CH);
%                     idxb_IN_prefHO_ch  =idx.in & idx.hands(:,hn) & ismember(epochs,b) & idx.LR(:,RF_out_hemifield_index_CH);
%                     idxb_CH_prefHI_ch  =idx.ch & idx.hands(:,hn) & ismember(epochs,b) & idx.LR(:,RF_in_hemifield_index_CH);
%                     idxb_CH_prefHO_ch  =idx.ch & idx.hands(:,hn) & ismember(epochs,b) & idx.LR(:,RF_out_hemifield_index_CH);
                    
                    %if numel(FR_LS_in)>=keys.cal.min_trials_per_condition && numel(FR_RS_in)>=keys.cal.min_trials_per_condition && sum(idx_CH_prefHI_in)>=keys.cal.min_trials_per_condition    
                    if min([numel(FR_LS_in),numel(FR_RS_in),sum(idx_CH_prefHO_in),sum(idx_CH_prefHI_in)])>=keys.cal.min_trials_per_condition 
                        p= sum([bootstrapped.ch_prefHI_FR]-[bootstrapped.in_prefHI_FR]<0)/n_boots;
                        p(p>0.5)=p-1;p=p*2;
                        %h= prctile([bootstrapped.ch_prefHI_FR]-[bootstrapped.in_prefHI_FR],2.5) >0 | prctile([bootstrapped.ch_prefHI_FR]-[bootstrapped.in_prefHI_FR],97.5) <0;
                        h= p<0.05;
                        DF=nanmean([bootstrapped.ch_prefHI_FR]-[bootstrapped.in_prefHI_FR]);
                        labelindex=h*sign(DF)+2; labelindex(isnan(labelindex))=2;
                        anova_struct.([IN '_' LHRH '_' s{:} '_prefH'])       =labels.choices{labelindex};
                        anova_struct.([CH '_' LHRH '_' s{:} '_prefH'])       =labels.choices{labelindex};
                        anova_struct.([IN '_' LHRH '_' s{:} '_prefH_PV'])    =p;
                        anova_struct.([IN '_' LHRH '_' s{:} '_prefH_FR'])    =nanmean([bootstrapped.in_prefHI_FR]);
                        anova_struct.([CH '_' LHRH '_' s{:} '_prefH_FR'])    =nanmean([bootstrapped.ch_prefHI_FR]);
                        anova_struct.([IN '_' LHRH '_' s{:} '_prefH_DF'])    =nanmean([bootstrapped.in_prefHI_FR]-[bootstrapped.in_prefHO_FR]);
                        anova_struct.([CH '_' LHRH '_' s{:} '_prefH_DF'])    =nanmean([bootstrapped.ch_prefHI_FR]-[bootstrapped.ch_prefHO_FR]);
                        anova_struct.([IN '_' LHRH '_' s{:} '_prefH_IX'])    =nanmean([bootstrapped.in_prefHI_FR]-[bootstrapped.in_prefHO_FR])/nanmean([bootstrapped.in_prefHI_FR,bootstrapped.in_prefHO_FR]);
                        anova_struct.([CH '_' LHRH '_' s{:} '_prefH_IX'])    =nanmean([bootstrapped.ch_prefHI_FR]-[bootstrapped.ch_prefHO_FR])/nanmean([bootstrapped.ch_prefHI_FR,bootstrapped.ch_prefHO_FR]);
                        
                        
                        %% these here had to be adjusted to bootstrapping method
                        anova_struct.([CH '_' LHRH '_' s{:} '_prefHI_FR'])   =nanmean(bootstrapped.ch_prefHI_FR);
                        anova_struct.([CH '_' LHRH '_' s{:} '_prefHO_FR'])   =nanmean(bootstrapped.ch_prefHO_FR);
                        anova_struct.([CH '_' LHRH '_' s{:} '_prefHI_DF'])   =nanmean(bootstrapped.ch_prefHI_FR)-nanmean(bootstrapped.ch_prefHI_FR_bl);
                        anova_struct.([CH '_' LHRH '_' s{:} '_prefHO_DF'])   =nanmean(bootstrapped.ch_prefHO_FR)-nanmean(bootstrapped.ch_prefHO_FR_bl);
                        
                        anova_struct.([IN '_' LHRH '_' s{:} '_prefHI_FR'])   =nanmean(bootstrapped.in_prefHI_FR);
                        anova_struct.([IN '_' LHRH '_' s{:} '_prefHO_FR'])   =nanmean(bootstrapped.in_prefHO_FR);
                        anova_struct.([IN '_' LHRH '_' s{:} '_prefHI_DF'])   =nanmean(bootstrapped.in_prefHI_FR)-nanmean(bootstrapped.in_prefHI_FR_bl);
                        anova_struct.([IN '_' LHRH '_' s{:} '_prefHO_DF'])   =nanmean(bootstrapped.in_prefHO_FR)-nanmean(bootstrapped.in_prefHO_FR_bl);
                        
                        % based on preferred choice
%                         anova_struct.([CH '_' LHRH '_' s{:}  '_prefHIch_FR'])   =nanmean(FR(idx_CH_prefHI_ch));
%                         anova_struct.([CH '_' LHRH '_' s{:}  '_prefHOch_FR'])   =nanmean(FR(idx_CH_prefHO_ch));
%                         anova_struct.([CH '_' LHRH '_' s{:}  '_prefHIch_DF'])   =nanmean(FR(idx_CH_prefHI_ch))-nanmean(FR(idxb_CH_prefHI_ch));
%                         anova_struct.([CH '_' LHRH '_' s{:}  '_prefHOch_DF'])   =nanmean(FR(idx_CH_prefHO_ch))-nanmean(FR(idxb_CH_prefHO_ch));
% 
%                         anova_struct.([IN '_' LHRH '_' s{:}  '_prefHIch_FR'])   =nanmean(FR(idx_IN_prefHI_ch));
%                         anova_struct.([IN '_' LHRH '_' s{:}  '_prefHOch_FR'])   =nanmean(FR(idx_IN_prefHO_ch));
%                         anova_struct.([IN '_' LHRH '_' s{:}  '_prefHIch_DF'])   =nanmean(FR(idx_IN_prefHI_ch))-nanmean(FR(idxb_IN_prefHI_ch));
%                         anova_struct.([IN '_' LHRH '_' s{:}  '_prefHOch_DF'])   =nanmean(FR(idx_IN_prefHO_ch))-nanmean(FR(idxb_IN_prefHO_ch));
                        
                    end
%                     %% based on preferred choice...
%                     if sum(idx_IN_prefHI_ch)>=keys.cal.min_trials_per_condition && sum(idx_CH_prefHI_ch)>=keys.cal.min_trials_per_condition
%                         [h,p,n]=do_stats(FR(idx_IN_prefHI_ch),FR(idx_CH_prefHI_ch),keys,0);
%                         DF=nanmean(FR(idx_CH_prefHI_ch))-nanmean(FR(idx_IN_prefHI_ch));
%                         labelindex=h*sign(DF)+2; labelindex(isnan(labelindex))=2;
%                         anova_struct.([IN '_' LHRH '_' s{:} '_prefHch'])       =labels.choices{labelindex};
%                         anova_struct.([CH '_' LHRH '_' s{:} '_prefHch'])       =labels.choices{labelindex};
%                         anova_struct.([IN '_' LHRH '_' s{:} '_prefHch_PV'])    =p*sign(DF);
%                         anova_struct.([CH '_' LHRH '_' s{:} '_prefHch_PV'])    =p*sign(DF);
%                         anova_struct.([IN '_' LHRH '_' s{:} '_prefHch_FR'])    =nanmean(FR(idx_IN_prefHI_ch));
%                         anova_struct.([CH '_' LHRH '_' s{:} '_prefHch_FR'])    =nanmean(FR(idx_CH_prefHI_ch));
%                         anova_struct=add_normality_results(anova_struct,[IN '_' LHRH '_' s{:} '_prefHch'],keys,n);
%                         anova_struct=add_normality_results(anova_struct,[CH '_' LHRH '_' s{:} '_prefHch'],keys,n);
%                     end

%                     if numel(FR_LS_in)>=keys.cal.min_trials_per_condition && numel(FR_RS_in)>=keys.cal.min_trials_per_condition && sum(idx_CH_prefHO_in)>=keys.cal.min_trials_per_condition
%                         anova_struct.([CH '_' LHRH '_' s{:} '_prefHI_FR'])   =nanmean(FR(idx_CH_prefHI_in));
%                         anova_struct.([CH '_' LHRH '_' s{:} '_prefHO_FR'])   =nanmean(FR(idx_CH_prefHO_in));
%                         anova_struct.([CH '_' LHRH '_' s{:} '_prefHI_DF'])   =nanmean(FR(idx_CH_prefHI_in))-nanmean(FR(idxb_CH_prefHI_in));
%                         anova_struct.([CH '_' LHRH '_' s{:} '_prefHO_DF'])   =nanmean(FR(idx_CH_prefHO_in))-nanmean(FR(idxb_CH_prefHO_in));
%                     end
%                     if sum(idx_IN_prefHI_in)>=keys.cal.min_trials_per_condition && sum(idx_IN_prefHO_in)>=keys.cal.min_trials_per_condition
%                         anova_struct.([IN '_' LHRH '_' s{:} '_prefHI_FR'])   =nanmean(FR(idx_IN_prefHI_in));
%                         anova_struct.([IN '_' LHRH '_' s{:} '_prefHO_FR'])   =nanmean(FR(idx_IN_prefHO_in));
%                         anova_struct.([IN '_' LHRH '_' s{:} '_prefHI_DF'])   =nanmean(FR(idx_IN_prefHI_in))-nanmean(FR(idxb_IN_prefHI_in));
%                         anova_struct.([IN '_' LHRH '_' s{:} '_prefHO_DF'])   =nanmean(FR(idx_IN_prefHO_in))-nanmean(FR(idxb_IN_prefHO_in));
%                     end
%                     
%                     if sum(idx_IN_prefHI_ch)>=keys.cal.min_trials_per_condition && sum(idx_CH_prefHI_ch)>=keys.cal.min_trials_per_condition
%                         [h,p,n]=do_stats(FR(idx_IN_prefHI_ch),FR(idx_CH_prefHI_ch),keys,0);
%                         DF=nanmean(FR(idx_CH_prefHI_ch))-nanmean(FR(idx_IN_prefHI_ch));
%                         labelindex=h*sign(DF)+2; labelindex(isnan(labelindex))=2;
%                         anova_struct.([IN '_' LHRH '_' s{:} '_prefHch'])       =labels.choices{labelindex};
%                         anova_struct.([CH '_' LHRH '_' s{:} '_prefHch'])       =labels.choices{labelindex};
%                         anova_struct.([IN '_' LHRH '_' s{:} '_prefHch_PV'])    =p*sign(DF);
%                         anova_struct.([CH '_' LHRH '_' s{:} '_prefHch_PV'])    =p*sign(DF);
%                         anova_struct.([IN '_' LHRH '_' s{:} '_prefHch_FR'])    =nanmean(FR(idx_IN_prefHI_ch));
%                         anova_struct.([CH '_' LHRH '_' s{:} '_prefHch_FR'])    =nanmean(FR(idx_CH_prefHI_ch));
%                         anova_struct=add_normality_results(anova_struct,[IN '_' LHRH '_' s{:} '_prefHch'],keys,n);
%                         anova_struct=add_normality_results(anova_struct,[CH '_' LHRH '_' s{:} '_prefHch'],keys,n);
%                     end
%                     if sum(idx_CH_prefHI_ch)>=keys.cal.min_trials_per_condition && sum(idx_CH_prefHO_ch)>=keys.cal.min_trials_per_condition
%                         anova_struct.([CH '_' LHRH '_' s{:}  '_prefHIch_FR'])   =nanmean(FR(idx_CH_prefHI_ch));
%                         anova_struct.([CH '_' LHRH '_' s{:}  '_prefHOch_FR'])   =nanmean(FR(idx_CH_prefHO_ch));
%                         anova_struct.([CH '_' LHRH '_' s{:}  '_prefHIch_DF'])   =nanmean(FR(idx_CH_prefHI_ch))-nanmean(FR(idxb_CH_prefHI_ch));
%                         anova_struct.([CH '_' LHRH '_' s{:}  '_prefHOch_DF'])   =nanmean(FR(idx_CH_prefHO_ch))-nanmean(FR(idxb_CH_prefHO_ch));
%                     end
%                     if sum(idx_IN_prefHI_ch)>=keys.cal.min_trials_per_condition && sum(idx_IN_prefHO_ch)>=keys.cal.min_trials_per_condition
%                         anova_struct.([IN '_' LHRH '_' s{:}  '_prefHIch_FR'])   =nanmean(FR(idx_IN_prefHI_ch));
%                         anova_struct.([IN '_' LHRH '_' s{:}  '_prefHOch_FR'])   =nanmean(FR(idx_IN_prefHO_ch));
%                         anova_struct.([IN '_' LHRH '_' s{:}  '_prefHIch_DF'])   =nanmean(FR(idx_IN_prefHI_ch))-nanmean(FR(idxb_IN_prefHI_ch));
%                         anova_struct.([IN '_' LHRH '_' s{:}  '_prefHOch_DF'])   =nanmean(FR(idx_IN_prefHO_ch))-nanmean(FR(idxb_IN_prefHO_ch));
%                     end
                    
                    
                    % preferred location - this one should be bootstrapped too i suppose
                    idx_IN_prefPI_in  =idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.pos==RF_position_index_IN;
                    idx_IN_prefPO_in  =idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.opp==RF_position_index_IN;
                    idx_CH_prefPI_in  =idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.pos==RF_position_index_IN;
                    idx_CH_prefPO_in  =idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.opp==RF_position_index_IN;
                    
                    idxb_IN_prefPI_in  =idx.in & idx.hands(:,hn) & ismember(epochs,b) & idx.pos==RF_position_index_IN;
                    idxb_IN_prefPO_in  =idx.in & idx.hands(:,hn) & ismember(epochs,b) & idx.opp==RF_position_index_IN;
                    idxb_CH_prefPI_in  =idx.ch & idx.hands(:,hn) & ismember(epochs,b) & idx.pos==RF_position_index_IN;
                    idxb_CH_prefPO_in  =idx.ch & idx.hands(:,hn) & ismember(epochs,b) & idx.opp==RF_position_index_IN;
                    
                    if sum(idx_IN_prefPI_in)>=keys.cal.min_trials_per_condition && sum(idx_CH_prefPI_in)>=keys.cal.min_trials_per_condition
                        [h,p,n]=do_stats(FR(idx_IN_prefPI_in),FR(idx_CH_prefPI_in),keys,0);
                        DF=nanmean(FR(idx_CH_prefPI_in))-nanmean(FR(idx_IN_prefPI_in));
                        labelindex=h*sign(DF)+2; labelindex(isnan(labelindex))=2;
                        anova_struct.([IN '_' LHRH '_' s{:}  '_prefP'])      =labels.choices{labelindex};
                        anova_struct.([CH '_' LHRH '_' s{:}  '_prefP'])      =labels.choices{labelindex};
                        anova_struct.([IN '_' LHRH '_' s{:}  '_prefP_PV'])   =p*sign(DF);
                        anova_struct.([CH '_' LHRH '_' s{:}  '_prefP_PV'])   =p*sign(DF);
                        anova_struct=add_normality_results(anova_struct,[IN '_' LHRH '_' s{:}  '_prefP'],keys,n);
                        anova_struct=add_normality_results(anova_struct,[CH '_' LHRH '_' s{:}  '_prefP'],keys,n);
                    end
                    
                    if sum(idx_IN_prefPI_in)>=keys.cal.min_trials_per_condition && sum(idx_CH_prefPI_in)>=keys.cal.min_trials_per_condition &&...
                            sum(idx_IN_prefPO_in)>=keys.cal.min_trials_per_condition && sum(idx_CH_prefPO_in)>=keys.cal.min_trials_per_condition
                        anova_struct.([IN '_' LHRH '_' s{:}  '_prefP_DF'])   =nanmean(FR(idx_IN_prefPI_in))-nanmean(FR(idx_IN_prefPO_in));
                        anova_struct.([CH '_' LHRH '_' s{:}  '_prefP_DF'])   =nanmean(FR(idx_CH_prefPI_in))-nanmean(FR(idx_CH_prefPO_in));
                    end
                    if sum(idx_IN_prefPI_in)>=keys.cal.min_trials_per_condition
                        anova_struct.([IN '_' LHRH '_' s{:}  '_prefPI_FR'])  =nanmean(FR(idx_IN_prefPI_in));
                        anova_struct.([IN '_' LHRH '_' s{:}  '_prefPI_DF'])  =nanmean(FR(idx_IN_prefPI_in))-nanmean(FR(idxb_IN_prefPI_in));
                    end
                    if sum(idx_IN_prefPO_in)>=keys.cal.min_trials_per_condition
                        anova_struct.([IN '_' LHRH '_' s{:}  '_prefPO_FR'])  =nanmean(FR(idx_IN_prefPO_in));
                        anova_struct.([IN '_' LHRH '_' s{:}  '_prefPO_DF'])  =nanmean(FR(idx_IN_prefPO_in))-nanmean(FR(idxb_IN_prefPO_in));
                    end
                    if sum(idx_CH_prefPI_in)>=keys.cal.min_trials_per_condition
                        anova_struct.([CH '_' LHRH '_' s{:}  '_prefPI_FR'])  =nanmean(FR(idx_CH_prefPI_in));
                        anova_struct.([CH '_' LHRH '_' s{:}  '_prefPI_DF'])  =nanmean(FR(idx_CH_prefPI_in))-nanmean(FR(idxb_CH_prefPI_in));
                    end
                    if sum(idx_CH_prefPO_in)>=keys.cal.min_trials_per_condition
                        anova_struct.([CH '_' LHRH '_' s{:}  '_prefPO_FR'])  =nanmean(FR(idx_CH_prefPO_in));
                        anova_struct.([CH '_' LHRH '_' s{:}  '_prefPO_DF'])  =nanmean(FR(idx_CH_prefPO_in))-nanmean(FR(idxb_CH_prefPO_in));
                    end
                    
                    % preferred location
                    idx_IN_prefPI_ch  =idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.pos==RF_position_index_CH;
                    idx_IN_prefPO_ch  =idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.opp==RF_position_index_CH;
                    idx_CH_prefPI_ch  =idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.pos==RF_position_index_CH;
                    idx_CH_prefPO_ch  =idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.opp==RF_position_index_CH;
                    
                    idxb_IN_prefPI_ch  =idx.in & idx.hands(:,hn) & ismember(epochs,b) & idx.pos==RF_position_index_CH;
                    idxb_IN_prefPO_ch  =idx.in & idx.hands(:,hn) & ismember(epochs,b) & idx.opp==RF_position_index_CH;
                    idxb_CH_prefPI_ch  =idx.ch & idx.hands(:,hn) & ismember(epochs,b) & idx.pos==RF_position_index_CH;
                    idxb_CH_prefPO_ch  =idx.ch & idx.hands(:,hn) & ismember(epochs,b) & idx.opp==RF_position_index_CH;
                    
                    if sum(idx_IN_prefPI_ch)>=keys.cal.min_trials_per_condition && sum(idx_CH_prefPI_ch)>=keys.cal.min_trials_per_condition
                        [h,p,n]=do_stats(FR(idx_IN_prefPI_ch),FR(idx_CH_prefPI_ch),keys,0);
                        DF=nanmean(FR(idx_CH_prefPI_ch))-nanmean(FR(idx_IN_prefPI_ch));
                        labelindex=h*sign(DF)+2; labelindex(isnan(labelindex))=2;
                        anova_struct.([IN '_' LHRH '_' s{:}  '_prefPch'])      =labels.choices{labelindex};
                        anova_struct.([CH '_' LHRH '_' s{:}  '_prefPch'])      =labels.choices{labelindex};
                        anova_struct.([IN '_' LHRH '_' s{:}  '_prefPch_PV'])   =p*sign(DF);
                        anova_struct.([CH '_' LHRH '_' s{:}  '_prefPch_PV'])   =p*sign(DF);
                        anova_struct=add_normality_results(anova_struct,[IN '_' LHRH '_' s{:}  '_prefPch'],keys,n);
                        anova_struct=add_normality_results(anova_struct,[CH '_' LHRH '_' s{:}  '_prefPch'],keys,n);
                    end
                    
                    if sum(idx_IN_prefPI_ch)>=keys.cal.min_trials_per_condition && sum(idx_CH_prefPI_ch)>=keys.cal.min_trials_per_condition &&...
                            sum(idx_IN_prefPO_ch)>=keys.cal.min_trials_per_condition && sum(idx_CH_prefPO_ch)>=keys.cal.min_trials_per_condition
                        anova_struct.([IN '_' LHRH '_' s{:}  '_prefPch_DF'])   =nanmean(FR(idx_IN_prefPI_ch))-nanmean(FR(idx_IN_prefPO_ch));
                        anova_struct.([CH '_' LHRH '_' s{:}  '_prefPch_DF'])   =nanmean(FR(idx_CH_prefPI_ch))-nanmean(FR(idx_CH_prefPO_ch));
                    end
                    if sum(idx_IN_prefPI_ch)>=keys.cal.min_trials_per_condition
                        anova_struct.([IN '_' LHRH '_' s{:}  '_prefPIch_FR'])   =nanmean(FR(idx_IN_prefPI_ch));
                        anova_struct.([IN '_' LHRH '_' s{:}  '_prefPIch_DF'])   =nanmean(FR(idx_IN_prefPI_ch))-nanmean(FR(idxb_IN_prefPI_ch));
                    end
                    if sum(idx_IN_prefPO_ch)>=keys.cal.min_trials_per_condition
                        anova_struct.([IN '_' LHRH '_' s{:}  '_prefPOch_FR'])   =nanmean(FR(idx_IN_prefPO_ch));
                        anova_struct.([IN '_' LHRH '_' s{:}  '_prefPOch_DF'])   =nanmean(FR(idx_IN_prefPO_ch))-nanmean(FR(idxb_IN_prefPO_ch));
                    end
                    if sum(idx_CH_prefPI_ch)>=keys.cal.min_trials_per_condition
                        anova_struct.([CH '_' LHRH '_' s{:}  '_prefPIch_FR'])   =nanmean(FR(idx_CH_prefPI_ch));
                        anova_struct.([CH '_' LHRH '_' s{:}  '_prefPIch_DF'])   =nanmean(FR(idx_CH_prefPI_ch))-nanmean(FR(idxb_CH_prefPI_ch));
                    end
                    if sum(idx_CH_prefPO_ch)>=keys.cal.min_trials_per_condition
                        anova_struct.([CH '_' LHRH '_' s{:}  '_prefPOch_FR'])   =nanmean(FR(idx_CH_prefPO_ch));
                        anova_struct.([CH '_' LHRH '_' s{:}  '_prefPOch_DF'])   =nanmean(FR(idx_CH_prefPO_ch))-nanmean(FR(idxb_CH_prefPO_ch));
                    end
                end
            else
                anova_struct.([INCH '_' LHRH '_' s{:} '_position'])=labels.hemifield{2};
            end
        end
    end
    
    %% perturbation ANOVA (!) 
    % temporarly recompute this stuff    :((
    tr=idx.(INCH) & ismember(epochs,epoch_multicomp) & (idx.PT==0 | idx.PT==1); %& ismember(epochs,epochs_for_multicomparison);
    tr_main=idx.(INCH) & ismember(epochs,keys.main_multicomp) & (idx.PT==0 | idx.PT==1); %& ismember(epochs,epochs_for_multicomparison);
    for fn=idx_fieldnames
        idx.(['tr_' fn{:}])=      idx.(fn{:}) & repmat(tr,1,size(idx.(fn{:}),2)) ;
    end
    idx.tr_CT=idx.tr_PT==0;
    
    %
    %     epoch_multicomp=keys.epoch_multicomp(idx_ep,:);
    %     multicomp_epochs=keys.epoch_SxH_multicomp;
    %     hand_space_multicomp=multicomp_epochs(ismember(multicomp_epochs,epochs'));
    
    perturbations_to_compare=unique(idx.tr_PT);
    if numel(perturbations_to_compare)>1
        labels.control_test={'EP','BL'};
        labels.conditions={'LH_LS','LH_RS','RH_LS','RH_RS','LS','RS'};
        idx.conditions=[idx.tr_LH_LS idx.tr_LH_RS idx.tr_RH_LS idx.tr_RH_RS idx.tr_LS idx.tr_RS];
        anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,'epoch_SxH_multicomp','tr_CT','tr_PT','PT');
        
        for c=1:numel(labels.conditions)
            idx_con=idx.(['tr_' labels.conditions{c}]);
            CON=labels.conditions{c};
            if any(idx.tr_PT(tr_main & idx_con)==0) && any(idx.tr_PT(tr_main & idx_con)==1)
                [anova_out.p,anova_out.table,anova_out.stats,anova_out.terms] = anovan(FR(tr_main & idx_con),[epoch_idx(tr_main & idx_con) idx.tr_PT(tr_main & idx_con)],'model','full','varnames',{'epoch', 'perturbation'},'display',keys.plot.anova_tables);
                anova_struct.([INCH '_' CON '_epoch_main'])    =anova_out.p(1)<0.05; %main effect on epoch!
                anova_struct.([INCH '_' CON '_PT_main'])       =anova_out.p(2)<0.05; %main effect on epoch!
                anova_struct.([INCH '_' CON '_ExP'])           =anova_out.p(3)<0.05; %main effect on epoch!
            end
            
            %% this part here i would basically like to remove in total...
            % a) you can get PT vs. CT test by using baseline correct AN
            % normalization
            % b) I'm not sure what is done with PT_epoch_DF afterwards
            for row=1:size(epoch_multicomp,1)
                s=epoch_multicomp(row,2);
                b=epoch_multicomp(row,1);
                idxS   = ismember(epochs,s);
                idxB   = ismember(epochs,b);
                idx1   = idx.tr_PT==0     & idx_con & idxS;
                idx1b  = idx.tr_PT==0     & idx_con & idxB;
                idx2   = idx.tr_PT==1     & idx_con & idxS;
                idx2b  = idx.tr_PT==1     & idx_con & idxB;
                
                if sum(idx1)>0 && sum(idx2)>0
                    [h,n]=do_stats(FR(idx1)-FR(idx1b),FR(idx2)-FR(idx2b),keys,0);
                    DFPT = nanmean(FR(idx2)-FR(idx2b));
                    DFCT = nanmean(FR(idx1)-FR(idx1b));
                    DF=DFPT-DFCT;
                else
                    h=false; n=NaN;
                    DFPT=single(NaN);
                    DFCT=single(NaN);
                    DF=single(NaN);
                end
                
                labelindexpt=h*sign(DF)+2; labelindexpt(isnan(labelindexpt))=2;
                anova_struct.([INCH '_' CON '_' s{:} '_PTbl']) = labels.PT{labelindexpt}; %
                anova_struct.([INCH '_' CON '_' s{:} '_PTbl_DF']) = DF; %
                anova_struct.([INCH '_' CON '_' s{:} '_PT_epoch_DF']) = DFPT; %
                anova_struct.([INCH '_' CON '_' s{:} '_CT_epoch_DF']) = DFCT; %
                prefix=[INCH '_' CON '_' s{:}];
                anova_struct=add_normality_results(anova_struct,prefix,keys,n);
            end
            
            
        end
    end
end
end

function anova_struct=effector_comparison_anova(keys,epochs,FR,idx,Positions,Fixations)
anova_struct=struct();
conditions=keys.TTconditions;
conditions.choice=conditions.choice([sum(idx.in)>0 sum(idx.ch)>0]);
labels=keys.TTlabels;

%% Independently for Instructed and Choice !!
idx_fieldnames=fieldnames(idx)';
for ch=1:numel(conditions.choice)
    
    %% this part should go outside the loop! ... NAH
    tr=idx.(conditions.choice{ch}); 
    labels.INCH=conditions.choice{ch};
    for fn=idx_fieldnames
        idx.(['tr_' fn{:}])=      idx.(fn{:}) & repmat(tr,1,size(idx.(fn{:}),2));
    end
    
    labels.control_test={};
    
    labels.conditions={''};
    idx.conditions=idx.tr_all;
    % epoch_SxH_multicomp needs to be replaced by something meaningful
    anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,'epoch_SxH_multicomp','tr_E1','tr_E2','effector');
        
    labels.conditions={'LH_LS','LH_RS','RH_LS','RH_RS','LS','RS'};
    idx.conditions=[idx.tr_LH_LS idx.tr_LH_RS idx.tr_RH_LS idx.tr_RH_RS idx.tr_LS idx.tr_RS];
    anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,'epoch_SxH_multicomp','tr_E1','tr_E2','effector');
    
%     %% 3-way ANOVA per epoch...
%     
%     multicomp_epochs=keys.epoch_SxH_multicomp; %% needs to be replaced by something meaningful
%     for s=multicomp_epochs(:)'
%         idxS= tr & ismember(epochs,s);
%         if any(~isnan(FR(idxS))) && any(idx.tr_LS(idxS)) && any(idx.tr_RS(idxS)) && any(idx.tr_LH(idxS)) && any(idx.tr_RH(idxS))
%             [anova_outs.p,anova_outs.table,anova_outs.stats,anova_outs.terms] = anovan(FR(idxS),double([idx.tr_RS(idxS),idx.tr_RH(idxS),idx.tr_E2(idxS)]),'model','full','varnames',{'space';'hands';'effector'},'display',keys.plot.anova_tables);
%         end
%     end
end
end

%% stats per epoch and condition
function anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,epochs_FN,idx1_FN,idx2_FN,fieldnamepart)
multicomp_epochs=keys.(epochs_FN);
multicomp_epochs=multicomp_epochs(ismember(multicomp_epochs,epochs'));
label=labels.(fieldnamepart);
idx1= idx.(idx1_FN);
idx2= idx.(idx2_FN);
INCH= labels.INCH;
COND=labels.conditions;

SEP='_';
for  c=1:size(idx.conditions,2)
    if isempty(COND{c})
        SEP='';
    end
    for s=multicomp_epochs(:)'
        idxS=ismember(epochs,s)  & idx.conditions(:,c);
        if any(idx1 & idxS) && any(idx2 & idxS)
            [h,p,n]=do_stats(FR(idx1 & idxS),FR(idx2 & idxS),keys,0);
            DF=nanmean(FR(idx2 & idxS))-nanmean(FR(idx1 & idxS));
            labelindex=h*sign(DF)+2; labelindex(isnan(labelindex))=2;
            prefix=[INCH '_' COND{c} SEP s{:} '_' fieldnamepart ];
            anova_struct.(prefix)=label{labelindex};
            anova_struct.([prefix '_DF'])=DF;
            anova_struct.([prefix '_PV'])=p;
            anova_struct.([prefix '_IX'])=DF/(nanmean(FR(idx2 & idxS))+ nanmean(FR(idx1 & idxS)));
            anova_struct=add_normality_results(anova_struct,prefix,keys,n);
            if ~isempty(labels.control_test)                
                anova_struct.([INCH '_' COND{c} SEP s{:} '_' labels.control_test{1} '_FR']) = nanmean(FR(idx1));
                anova_struct.([INCH '_' COND{c} SEP s{:} '_' labels.control_test{2} '_FR']) = nanmean(FR(idx2));
            end
        end
    end
end
end

%% similarly, make epoch effects
function anova_struct=do_epoch_stats(keys,anova_struct,FR,epochs,idx,labels,epochs_FN,fieldnamepart)
%% epoch tuning per hand
multicomp_epochs=keys.(epochs_FN);
idx_ep=ismember(multicomp_epochs(:,1),epochs') &  ismember(multicomp_epochs(:,2),epochs');
multicomp_epochs=multicomp_epochs(idx_ep,:);
label=labels.(fieldnamepart);
INCH= labels.INCH;
COND=labels.conditions;
for row=1:size(multicomp_epochs,1)
    s=multicomp_epochs(row,2);
    b=multicomp_epochs(row,1);
    for  c=1:size(idx.conditions,2)
        idx1=ismember(epochs,b) & idx.conditions(:,c); 
        idx2=ismember(epochs,s) & idx.conditions(:,c);
        if sum(idx1)>1 && sum(idx2)>1
            [h,p,n]=do_stats(FR(idx1),FR(idx2),keys,0); %%unpaired ?? 
        else
            h=false; n=NaN;
        end
        DF=(nanmean(FR(idx2))-nanmean(FR(idx1)));
        labelindex=h*sign(DF)+2; labelindex(isnan(labelindex))=2;
        if isempty(COND{c})
            prefix=[INCH '_' s{:} '_' fieldnamepart];
        else
            prefix=[INCH '_' COND{c} '_' s{:} '_' fieldnamepart];
        end
        anova_struct.([prefix '_FR'])=nanmean(FR(idx2));
        anova_struct.([prefix '_DF'])= DF;
        anova_struct.([prefix '_PV']) = p*sign(DF);
        anova_struct.(prefix)        = label{labelindex};
        anova_struct=add_normality_results(anova_struct,prefix,keys,n);
    end
end
end

%% embedded test selection
function [h,p,n]=do_stats(A,B,keys,paired)
switch keys.AN.test_types
    case 'parametric'
        if paired
            if any(~isnan(A)&~isnan(B))
                h = ttest(A,B);
            else
                h=0;
            end
        else
            if any(~isnan(A)) && any (~isnan(B))
                h = ttest2(A,B);
            else
                h=0;
            end
        end
    case 'nonparametric'
        if paired
            if any(~isnan(A)&~isnan(B))
                [p, h] = signrank(A,B);
            else
                h=0;
            end
        else
            if any(~isnan(A)) && any (~isnan(B))
                [p, h] = ranksum(A,B);
            else
                h=0;
            end
        end
    case 'permutations'
%     [p] = permtest(A,B,1000,'conservative');    
%     h=p<0.05;
    [h,p]=DAG_permutation_test(A,B,10000);
end
p=single(p);
if keys.AN.check_normality && numel(A)>4 && numel(B)>4
    %[H,pValue,KSstatistic,criticalValue] = lillietest(x,alpha,distr,mctol)
    n=lillietest(A)==0 && lillietest(B)==0;
else
    n=1;
end
end

%% normality subfunction

function anova_struct=add_normality_results(anova_struct,prefix,keys,n)
if keys.AN.check_normality
    if isnan(n)
        entry='NA';
    elseif n==1
        entry='YE';
    elseif n==0
        entry='NO';
    end
    anova_struct.([prefix '_NM'])=entry;
end
end

%% this function is almost ready to be implemented
function [oo,keys]=ph_FR_at_peak(o,keys)
oo=o;
for t=1:numel(o)
    for s=1:size(keys.EPOCHS,1)
        restructured_by_states(s).trial(t)=o(t).spikes_per_state(s);
        restructured_by_states(s).state=o(t).spikes_per_state(s).state;
    end
end

for s=1:numel(restructured_by_states)
    state=restructured_by_states(s).state;
    state_idx=ismember(keys.EPOCHS(:,1),state);
    t_before_state=keys.EPOCHS{state_idx,3};
    t_after_state=keys.EPOCHS{state_idx,4};
    bins=t_before_state+keys.PSTH_binwidth/2:keys.PSTH_binwidth:t_after_state-keys.PSTH_binwidth/2;
    histo=hist(vertcat(restructured_by_states(s).trial.arrival_times),bins);
    histo=smooth(histo-histo(1),1,1.5);
    t_max=find(abs(histo)==max(abs(histo)));
    peak_location(s)=bins(t_max(1));
    keys.EPOCHS{state_idx,3}=peak_location(s)-keys.PSTH_binwidth/2;
    keys.EPOCHS{state_idx,4}=peak_location(s)+keys.PSTH_binwidth/2;
end
for t=1:numel(o)
    for s=1:numel(o(t).spikes_per_state)
        at=o(t).spikes_per_state(s).arrival_times;
        arrival_times=at(at>peak_location(s)-keys.PSTH_binwidth & at<peak_location(s)+keys.PSTH_binwidth);
        oo(t).spikes_per_state(s).FR=numel(arrival_times)/keys.PSTH_binwidth;
    end
end
end

