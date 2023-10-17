function TT=ph_ANOVAS(population,trials,keys)
TT=keys.tuning_table;
for unit=1:numel(population)
    display(['anovas for: ' population(unit).unit_ID]); 
    anova_struct_current_unit=struct();
    pop=population(unit);
    pt=ph_get_unit_trials(pop,trials);
    valid_trials=[pt.accepted]==1 & [pt.completed]==1;
    if sum(valid_trials)>0 % temporary? for pulv_oculomotor_dataset
        pop.trial=pop.trial(valid_trials);        
        pt=pt(valid_trials);
        fn={'perturbation','accepted','block','n','run'};
        for f=fn
            pop.(f{:})=[pt.(f{:})];
        end
        keys.labels.reach_hand=keys.labels.reach_handLR ;
        [UC, CM, lables]=ph_get_condition_matrix(pt,keys);%% brakes if not all tasktypes are defined
        keys.normalization_field='AN';
        o=ph_condition_normalization(pop,pt,keys,UC,CM); %condition wise normalization (also reduces conditions!???)
        for type=UC.type
            keys.anova_epochs=keys.ANOVAS_PER_TYPE(type); %% this here should eventually replace having to repeat the above
            
            %% within effectors comparison
            for effector=UC.effector
                keys=ph_get_epoch_keys(keys,type,effector,sum(UC.type_effector(:,1)==type)>1);
                [~, condition_fieldname_part]=MPA_get_type_effector_name(type,effector);
                tr_index= [pt.effector]==effector & [pt.type]==type ;
                if sum(tr_index)==0;
                    fprintf('no trials for effector %.0f type %.0f \n',effector,type);continue;
                end
                o_e=o;
                o_e.trial=o.trial(tr_index);
                o_e.EPOCH=o.epochs_per_type{type};
                
                [UCt, CMt, lables]=ph_get_condition_matrix(pt(tr_index),keys); % this is to make sure we are only taking present conditions for trial criterion
                trial_criterion=ph_get_minimum_trials(keys,pt(tr_index),CMt,UCt,lables);
                [FR,epochs,idx,u_pos,u_fix]=ph_get_anova_inputs(o_e,pt(tr_index),keys);
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
        %         end
    end
    rows_to_update=find(ismember(TT(:,DAG_find_column_index(TT,'unit_ID')),{population(unit).unit_ID}));
    if isempty(rows_to_update)
        rows_to_update=size(TT,1)+1;
    end
    clear unit_table;
    inital_fieldnames={'unit_ID','monkey','target','perturbation_site','grid_x','grid_y','electrode_depth','FR','n_spikes','avg_stability','avg_SNR','avg_single_rating','waveform_width'};
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
    titles_original_table=TT(1,:);
    [~, title_positions]=ismember(titles_update,titles_original_table);
    for n_update_column=1:numel(titles_update)
        XXX=title_positions(n_update_column);
        if XXX==0
            n_column=size(TT,2)+1;
            TT(1,n_column)=titles_update(n_update_column);
        else
            n_column=XXX;
        end
        TT(rows_to_update,n_column)=unit_table(2,n_update_column);
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
conditions.perturbation=conditions.perturbation([sum(idx.PT==0)>0 sum(idx.PT==1)>0]);
pt_values=[0,1];
pt_values=pt_values([sum(idx.PT==0)>0 sum(idx.PT==1)>0]);

[~, ~, epoch_idx]=unique(epochs);
idx_fieldnames=fieldnames(idx)';

% epochs for epoch comparison
epoch_multicomp=keys.anova_epochs.epoch;
idx_ep=ismember(epoch_multicomp(:,1),epochs') &  ismember(epoch_multicomp(:,2),epochs');
epoch_multicomp=epoch_multicomp(idx_ep,:);
anova_struct=struct();
%% Independently for Instructed and Choice !!
for ch=1:numel(conditions.choice)
    INCH=conditions.choice{ch};
    labels.INCH=conditions.choice{ch};
    for pt=1:numel(pt_values)
        pt_value=pt_values(pt);
        PT=conditions.perturbation{pt};
        labels.PT=conditions.perturbation{pt};
        %% perturbation here is defined as always the lowest accessible for control... (???)
        tr_ch=idx.(INCH) & idx.PT==pt_value;
        for fn=idx_fieldnames
            idx.(['tr_' fn{:}])=      idx.(fn{:}) & repmat(tr_ch,1,size(idx.(fn{:}),2)) ;
            idx.(['pt_' fn{:}])=      idx.(fn{:}) & repmat(idx.PT==1,1,size(idx.(fn{:}),2)) ;
        end
        
        %% number of trials for each condition to decide about type of (main) anova
        n_lr=sum(idx.tr_LS)>0 && sum(idx.tr_RS)>0;
        n_2hands=sum(idx.tr_LH_LS)>0 && sum(idx.tr_LH_RS)>0 && sum(idx.tr_RH_LS)>0 && sum(idx.tr_RH_RS)>0;
        if n_2hands % three-way anova epoch, space, hand
            %anova_varnames={'epoch' 'hemifield' 'hands'};
            tr=idx.tr_LH_LS | idx.tr_LH_RS | idx.tr_RH_LS | idx.tr_RH_RS;
            %Par=[epoch_idx, idx.tr_RS, idx.tr_RH];
        elseif n_lr % two-way anova epoch, space
            %anova_varnames={'epoch' 'hemifield'};
            tr=idx.tr_LS | idx.tr_RS;
            %Par=[epoch_idx, idx.tr_RS];
        else
            disp(['Not fitting dataset... no left or right ' INCH PT ' present?'])
            %anova_struct.([INCH PT '_epoch_main'])=[];
            continue;
        end
        
        %% this part is basically useless (?)
        %
        %         %% main effects and general interactions
        %         tr_main=idx.(INCH) & ismember(epochs,keys.anova_epochs.main);
        %         pval = ph_anova1(FR(tr_main),Par(tr_main,:),keys);
        %         h=pval<0.05; %main effect on epoch!
        %         anova_struct.([INCH PT '_epoch_main'])=labels.true{h+2};
        %         handindexes=[1,2,3];
        %         handindexes=handindexes(any(idx.tr_hands));
        %
        %         %% Epoch X space/hand tuning (regardless of the other)
        %         for k=2:numel(anova_varnames)
        %             multicomp_epochs=keys.anova_epochs.(anova_varnames{k});
        %             %multicomp_epochs=multicomp_epochs(ismember(multicomp_epochs,epochs')); %% to think about: add &  ismember(keys.epoch_multicomp(:,2),keys.EPOCHS(:,1));
        %             label=labels.(anova_varnames{k});
        %             idx1m=tr_main & Par(:,k)==0 & ismember(epochs,multicomp_epochs);
        %             idx2m=tr_main & Par(:,k)==1 & ismember(epochs,multicomp_epochs);
        %             labelindex=(pval(k)<0.05) *sign(nanmean(FR(idx2m))-nanmean(FR(idx1m)))+2; labelindex(isnan(labelindex))=2;
        %             anova_struct.([INCH PT '_' anova_varnames{k} '_main'])=label{labelindex}; %main effect with direction!
        %             hx=pval(k+numel(anova_varnames)-1)<0.05;
        %             anova_struct.([INCH PT '_Ex' upper(anova_varnames{k}(1))])=labels.true{hx+2};
        %         end
        
        %% space x hand anovas per epoch
        if n_2hands
            k=6;
            tr_SXH=tr_main & ismember(epochs,keys.anova_epochs.SxH);
            labelindexcr=(pval(k)<0.05) *sign(nanmean(FR(tr_SXH & idx.tr_CR))-nanmean(FR(tr_SXH & idx.tr_UC)))+2; labelindexcr(isnan(labelindexcr))=2;% 3 crossed>uncrossed ,1 crossed<uncrossed
            anova_struct.([INCH PT '_SxH'])=labels.CR{labelindexcr};
            anova_struct.([INCH PT '_ExSxH'])=pval(7)<0.05;
            
            for s=keys.anova_epochs.SxH'
                idxS= tr & ismember(epochs,s);
                if any(~isnan(FR(idxS))) && any(idx.tr_LS(idxS)) && any(idx.tr_RS(idxS)) && any(idx.tr_LH(idxS)) && any(idx.tr_RH(idxS))
                    pval = ph_anova2(FR(idxS),{idx.tr_RS(idxS),idx.tr_RH(idxS)},1,keys);
                    h=pval(1)<0.05;
                    DF=nanmean(FR(idx.tr_RS & idxS))-nanmean(FR(idx.tr_LS & idxS));
                    labelindexsp=h*sign(DF)+2; labelindexsp(isnan(labelindexsp))=2;
                    anova_struct.([INCH PT '_' s{:} '_SxH_hemifield'])     = labels.hemifield{labelindexsp};
                    anova_struct.([INCH PT '_' s{:} '_SxH_hemifield_PV'])  = single(pval(1)*sign(DF));
                    anova_struct.([INCH PT '_' s{:} '_SxH_hemifield_DF'])  = DF;
                    anova_struct.([INCH PT '_' s{:} '_SxH_hemifield_IX'])  = DF/(nanmean(FR(idx.tr_RS & idxS))+nanmean(FR(idx.tr_LS & idxS)));
                    
                    h=pval(2)<0.05;
                    DF=nanmean(FR(idx.tr_RH & idxS))-nanmean(FR(idx.tr_LH & idxS));
                    labelindexha=h*sign(DF)+2; labelindexha(isnan(labelindexha))=2;
                    anova_struct.([INCH PT '_' s{:} '_SxH_hands'])         = labels.hands{labelindexha};
                    anova_struct.([INCH PT '_' s{:} '_SxH_hands_PV'])      = single(pval(2)*sign(DF));
                    anova_struct.([INCH PT '_' s{:} '_SxH_hands_DF'])      = DF;
                    anova_struct.([INCH PT '_' s{:} '_SxH_hands_IX'])      = DF/(nanmean(FR(idx.tr_RH & idxS))+nanmean(FR(idx.tr_LH & idxS)));
                    
                    h=pval(3)<0.05;
                    DF=nanmean(FR(idx.tr_CR & idxS))-nanmean(FR(idx.tr_UC & idxS));
                    labelindexcr=h*sign(DF)+2; labelindexcr(isnan(labelindexcr))=2;
                    anova_struct.([INCH PT '_' s{:} '_SxH'])               = labels.CR{labelindexcr};
                    anova_struct.([INCH PT '_' s{:} '_SxH_PV'])            = single(pval(3)*sign(DF));
                    anova_struct.([INCH PT '_' s{:} '_SxH_DF'])            = DF;
                    anova_struct.([INCH PT '_' s{:} '_SxH_IX'])            = DF/(nanmean(FR(idx.tr_CR & idxS))+nanmean(FR(idx.tr_UC & idxS)));
                end
            end
        end
        
        %% epoch tuning across all trials
        labels.conditions={''};
        idx.conditions=idx.tr_all;
        anova_struct=do_epoch_stats(keys,anova_struct,FR,epochs,idx,labels,'epoch','epoch');
        
        %% epoch tuning per hand & space
        labels.conditions={'AH_LS','AH_RS','LH_LS','LH_RS','RH_LS','RH_RS'};
        idx.conditions=[idx.tr_LS idx.tr_RS idx.tr_LH_LS idx.tr_LH_RS idx.tr_RH_LS idx.tr_RH_RS];
        anova_struct=do_epoch_stats(keys,anova_struct,FR,epochs,idx,labels,'epoch','epoch');
        
        %     %% epoch tuning per hand OR space
        %     labels.conditions={'AH','LH','RH','LS','RS'};
        %     idx.conditions=[idx.tr_AH idx.tr_LH idx.tr_RH idx.tr_LS idx.tr_RS];
        %     anova_struct=do_epoch_stats(keys,anova_struct,FR,epochs,idx,labels,'epoch','epoch');
        
        %% hand and space tuning per epoch (across all trials)
        labels.conditions={''};
        idx.conditions=idx.tr_all;
        anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,'hemifield','tr_LS','tr_RS','hemifield');
        anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,'hands','tr_LH','tr_RH','hands');
        
        %% hemifield tuning per hand
        labels.conditions=conditions.hands;
        idx.conditions=idx.tr_hands;
        anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,'SxH','tr_LS','tr_RS','hemifield');
        
        %% hand tuning per hemifield
        labels.conditions=conditions.hemifield;
        idx.conditions=idx.sides;
        anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,'SxH','tr_LH','tr_RH','hands');
        
        %%  Difficulty stuff
        labels.conditions=conditions.hemifield;
        idx.conditions=idx.tr_sides;
        anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,'hemifield','tr_Diff0','tr_Diff1','Difficulty_Easy');
        anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,'hemifield','tr_Diff0','tr_Diff2','Difficulty_Diff');
        anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,'hemifield','tr_nonDistr1','tr_TT1HF','SpatialComp_1HFTar');
        anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,'hemifield','tr_nonDistr1','tr_TT2HF','SpatialComp_2HFTar');
        
        %%  single trials per space - this did not do what it's supposed to do !
        labels.conditions={''};
        idx.conditions=idx.tr_all;
        anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,'hemifield','tr_SglL','tr_SglR','SglTar_Suc');
        
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
        %             anova_struct.([INCH PT '_' LSRSnamepart{sideindex} '_' s{:} '_' 'SglTar_Suc' ])=label{labelindex};
        %             anova_struct.([INCH PT '_' LSRSnamepart{sideindex} '_' s{:} '_' 'SglTar_Suc' '_DF'])=DF;
        %             anova_struct.([INCH PT '_' LSRSnamepart{sideindex} '_' s{:} '_' 'SglTar_Suc' '_IX'])=DF/(nanmean(FR(idx2 & idxS))+ nanmean(FR(idx1 & idxS)));
        %         end
        %     end
        
        %% Position anova (per hand) + comparison CH vs IN
        epochs_for_position_comparison=keys.anova_epochs.positions(ismember(keys.anova_epochs.positions,keys.EPOCHS(:,1)));
        multicomp_epochs=epochs_for_position_comparison(ismember(epochs_for_position_comparison,epochs));
        
        u_hnd                       =find(any(idx.tr_hands,1));
        for hn=u_hnd
            
            LHRH=conditions.hands{hn};
            tr=idx.tr_hands(:,hn) & ismember(epochs,epochs_for_position_comparison);
            
            if ~any(~isnan(FR(tr)))
                continue
            end
            %             if size(Fixations,1)==1
            %                 pval = ph_anova2(FR(tr),{epochs(tr),idx.pos(tr)},1,keys);
            %
            %             else
            %                 varnames={'State','position','fixation'};
            %                 pval = ph_anova1(FR(tr),{epochs(tr),idx.pos(tr),idx.fix(tr)},keys);
            %                 h3=pval(3)<0.05;
            %                 anova_struct.([INCH PT '_' LHRH '_fixation_main'])=labels.true{h3+2}; %main effect on epoch!
            %             end
            %             h1=pval(1)<0.05;
            %             h2=pval(2)<0.05;
            %             anova_struct.([INCH PT '_' LHRH '_epoch_main'])=labels.true{h1+2}; %main effect on epoch!
            %             anova_struct.([INCH PT '_' LHRH '_position_main'])=labels.true{h2+2}; %main effect on epoch!
            %
            
            for s=multicomp_epochs(:)'
                tr=idx.tr_hands(:,hn) & ismember(epochs,s);  %this is not ideal, because tr is used in a different way here
                b=epoch_multicomp(ismember(epoch_multicomp(:,2),s),1);
                if ~any(~isnan(FR(tr)))
                    anova_struct.([INCH PT '_' LHRH '_' s{:} '_position'])=labels.hemifield{2};
                    continue
                end
                %pval = anovan(FR(tr),[idx.pos(tr),idx.fix(tr)],'model','full','varnames',{'Positions','Fixations'},'display',keys.plot.anova_tables);
                pos_fix=[idx.pos(tr) idx.fix(tr)];
                u_pos_fix=unique(pos_fix,'rows');
                u_pos_fix1_valid=find(hist(u_pos_fix(:,1),unique(idx.pos(tr)))>1);
                u_pos_fix2_valid=find(hist(u_pos_fix(:,2),unique(idx.fix(tr)))>1);
                u_pos_fix=u_pos_fix(ismember(u_pos_fix(:,1),u_pos_fix1_valid') & ismember(u_pos_fix(:,2),u_pos_fix2_valid'),:);
                pos_fix_valid=ismember(pos_fix,u_pos_fix,'rows'); % this checks for at least 2 positions for every fixation and at least 2 fixations for every positon
                
                if any(~pos_fix_valid) %this is the case for target position by location in pulvinar eye gaze project AND for fixation only
                    % take only positions which have combinations of different fixations to compute interaction
                    
                    FR_temp=FR(tr);
                    if any(pos_fix_valid)
                        ptemp = ph_anova2(FR_temp(pos_fix_valid),pos_fix(pos_fix_valid,:),1,keys);
                    else
                        ptemp(2)=NaN;
                        ptemp(3)=NaN;
                    end
                    
                    pval = ph_anova1(FR(tr),idx.pos(tr),keys);
                    pval(2)=NaN;
                    pval(3)=ptemp(2)>0.05 || ptemp(3)>0.05; % this is not ideal, interaction here only significant if there is a main effect as well....
                else
                    pval= ph_anova2(FR(tr),{idx.pos(tr),idx.fix(tr)},1,keys);
                end
                
                h1=pval(1)<0.05;
                h2=pval(2)<0.05;
                h3=pval(3)<0.05;
                
                anova_struct.([INCH PT '_' LHRH '_' s{:} '_position'])     =labels.true{h1+2};
                anova_struct.([INCH PT '_' LHRH '_' s{:} '_fixation'])     =labels.true{h2+2};
                anova_struct.([INCH PT '_' LHRH '_' s{:} '_PxF'])          =labels.true{h3+2};
                anova_struct.([INCH PT '_' LHRH '_' s{:} '_position_PV'])  =single(pval(1));
                anova_struct.([INCH PT '_' LHRH '_' s{:} '_fixation_PV'])  =single(pval(2));
                anova_struct.([INCH PT '_' LHRH '_' s{:} '_PxF_PV'])       =single(pval(3));
                
                pecc = ph_anova2(FR(tr),{idx.ecc(tr),idx.ang(tr)},1,keys);
                h1=pecc(1)<0.05;
                h2=pecc(2)<0.05;
                h3=pecc(3)<0.05;
                anova_struct.([INCH PT '_' LHRH '_' s{:} '_distance'])     =labels.true{h1+2};
                anova_struct.([INCH PT '_' LHRH '_' s{:} '_angle'])        =labels.true{h2+2};
                anova_struct.([INCH PT '_' LHRH '_' s{:} '_DxA'])          =labels.true{h3+2};
                anova_struct.([INCH PT '_' LHRH '_' s{:} '_distance_PV'])  =single(pecc(1));
                anova_struct.([INCH PT '_' LHRH '_' s{:} '_angle_PV'])     =single(pecc(2));
                anova_struct.([INCH PT '_' LHRH '_' s{:} '_DxA_PV'])       =single(pecc(3));
                
                % x and y separately - didnt add p_values here yet!
                pxy = ph_anova2(FR(tr),{idx.pos_x(tr),idx.pos_y(tr)},1,keys);
                h1=pxy(1)<0.05;
                h2=pxy(2)<0.05;
                h3=pxy(3)<0.05;
                anova_struct.([INCH PT '_' LHRH '_' s{:} '_positionx'])=labels.true{h1+2};
                anova_struct.([INCH PT '_' LHRH '_' s{:} '_positiony'])=labels.true{h2+2};
                anova_struct.([INCH PT '_' LHRH '_' s{:} '_positionxy'])=labels.true{h3+2};
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
                            anova_struct.([INCH PT '_' LHRH '_' s{:} '_gaze_modulation_' xorytag])='monotonous';
                            anova_struct.([INCH PT '_' LHRH '_' s{:} '_gaze_pref_' xorytag])=labelxy{(directionality_per_position(end)==1)*2+1};
                        else
                            anova_struct.([INCH PT '_' LHRH '_' s{:} '_gaze_modulation_' xorytag])='nonmonotonous';
                            if pref_idx==1 || pref_idx==numel(u_xy_pos)
                                anova_struct.([INCH PT '_' LHRH '_' s{:} '_gaze_pref_' xorytag])='PE';
                            else
                                anova_struct.([INCH PT '_' LHRH '_' s{:} '_gaze_pref_' xorytag])='CE';
                            end
                        end
                    else % not even nonmonotoneous
                        anova_struct.([INCH PT '_' LHRH '_' s{:} '_gaze_modulation_' xorytag])='-';
                        anova_struct.([INCH PT '_' LHRH '_' s{:} '_gaze_pref_' xorytag])='-';
                    end
                end
                
                %% choice part for position with strongest response
                if ch==2  && isfield(anova_struct,['in' PT '_' LHRH '_' s{:} '_epoch_DF']) %% there has to be choice and instructed
                    
                    % preference (highest or lowest) defined by over all
                    % increase or decrease in this epoch
                    if false % anova_struct.(['in_' LHRH '_' s{:} '_epoch_DF']) <0
                        invertsign=-1;
                        takemin=1;
                    else
                        invertsign=1;
                        takemin=0;
                    end
                    
                    % bootstrapping choice preference: for every iteration,
                    % 50% of trials to each hemifield are taken to estimate
                    % preferred hemifield, and the difference is computed
                    % based on the remaining 50%
                    
                    n_boots=100;
                    RF_IN_RS_IDX=idx.pt_in & idx.hands(:,hn) & ismember(epochs,s) & idx.RS;
                    RF_IN_LS_IDX=idx.pt_in & idx.hands(:,hn) & ismember(epochs,s) & idx.LS;
                    RF_IN_RS_IDX_bl=idx.pt_in & idx.hands(:,hn) & ismember(epochs,b) & idx.RS;
                    RF_IN_LS_IDX_bl=idx.pt_in & idx.hands(:,hn) & ismember(epochs,b) & idx.LS;
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
                        idx_CH_prefHI_in  =idx.pt_ch & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_in_hemifield_index_IN);
                        idx_CH_prefHO_in  =idx.pt_ch & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_out_hemifield_index_IN);
                        idx_CH_prefHI_in_bl  =idx.pt_ch & idx.hands(:,hn) & ismember(epochs,b) & idx.LR(:,RF_in_hemifield_index_IN);
                        idx_CH_prefHO_in_bl  =idx.pt_ch & idx.hands(:,hn) & ismember(epochs,b) & idx.LR(:,RF_out_hemifield_index_IN);
                        
                        bootstrapped.ch_prefHI_FR(boots)=nanmean(FR(idx_CH_prefHI_in));
                        bootstrapped.ch_prefHO_FR(boots)=nanmean(FR(idx_CH_prefHO_in));
                        bootstrapped.ch_prefHI_FR_bl(boots)=nanmean(FR(idx_CH_prefHI_in_bl));
                        bootstrapped.ch_prefHO_FR_bl(boots)=nanmean(FR(idx_CH_prefHO_in_bl));
                    end
                    
                    if min([numel(FR_LS_in),numel(FR_RS_in),sum(idx_CH_prefHO_in),sum(idx_CH_prefHI_in)])>=keys.cal.min_trials_pref
                        p= sum([bootstrapped.ch_prefHI_FR]-[bootstrapped.in_prefHI_FR]<0)/n_boots;
                        p(p>0.5)=p-1;p=p*2;
                        %h= prctile([bootstrapped.ch_prefHI_FR]-[bootstrapped.in_prefHI_FR],2.5) >0 | prctile([bootstrapped.ch_prefHI_FR]-[bootstrapped.in_prefHI_FR],97.5) <0;
                        h= abs(p)<0.05;
                        DF=nanmean([bootstrapped.ch_prefHI_FR]-[bootstrapped.in_prefHI_FR]);
                        labelindex=h*sign(DF)+2; labelindex(isnan(labelindex))=2;
                        anova_struct.([IN PT '_' LHRH '_' s{:} '_prefH'])       =labels.choices{labelindex};
                        anova_struct.([CH PT '_' LHRH '_' s{:} '_prefH'])       =labels.choices{labelindex};
                        anova_struct.([IN PT '_' LHRH '_' s{:} '_prefH_PV'])    =p;
                        anova_struct.([IN PT '_' LHRH '_' s{:} '_prefH_FR'])    =nanmean([bootstrapped.in_prefHI_FR]);
                        anova_struct.([CH PT '_' LHRH '_' s{:} '_prefH_FR'])    =nanmean([bootstrapped.ch_prefHI_FR]);
                        anova_struct.([IN PT '_' LHRH '_' s{:} '_prefH_DF'])    =nanmean([bootstrapped.in_prefHI_FR]-[bootstrapped.in_prefHO_FR]);
                        anova_struct.([CH PT '_' LHRH '_' s{:} '_prefH_DF'])    =nanmean([bootstrapped.ch_prefHI_FR]-[bootstrapped.ch_prefHO_FR]);
                        anova_struct.([IN PT '_' LHRH '_' s{:} '_prefH_IX'])    =nanmean([bootstrapped.in_prefHI_FR]-[bootstrapped.in_prefHO_FR])/nanmean([bootstrapped.in_prefHI_FR,bootstrapped.in_prefHO_FR]);
                        anova_struct.([CH PT '_' LHRH '_' s{:} '_prefH_IX'])    =nanmean([bootstrapped.ch_prefHI_FR]-[bootstrapped.ch_prefHO_FR])/nanmean([bootstrapped.ch_prefHI_FR,bootstrapped.ch_prefHO_FR]);
                        
                        
                        %% these here had to be adjusted to bootstrapping method
                        anova_struct.([CH PT  '_' LHRH '_' s{:} '_prefHI_FR'])   =nanmean(bootstrapped.ch_prefHI_FR);
                        anova_struct.([CH PT  '_' LHRH '_' s{:} '_prefHO_FR'])   =nanmean(bootstrapped.ch_prefHO_FR);
                        anova_struct.([CH PT  '_' LHRH '_' s{:} '_prefHI_DF'])   =nanmean(bootstrapped.ch_prefHI_FR)-nanmean(bootstrapped.ch_prefHI_FR_bl);
                        anova_struct.([CH PT  '_' LHRH '_' s{:} '_prefHO_DF'])   =nanmean(bootstrapped.ch_prefHO_FR)-nanmean(bootstrapped.ch_prefHO_FR_bl);
                        
                        anova_struct.([IN PT '_' LHRH '_' s{:} '_prefHI_FR'])   =nanmean(bootstrapped.in_prefHI_FR);
                        anova_struct.([IN PT '_' LHRH '_' s{:} '_prefHO_FR'])   =nanmean(bootstrapped.in_prefHO_FR);
                        anova_struct.([IN PT '_' LHRH '_' s{:} '_prefHI_DF'])   =nanmean(bootstrapped.in_prefHI_FR)-nanmean(bootstrapped.in_prefHI_FR_bl);
                        anova_struct.([IN PT '_' LHRH '_' s{:} '_prefHO_DF'])   =nanmean(bootstrapped.in_prefHO_FR)-nanmean(bootstrapped.in_prefHO_FR_bl);
                    end
                    
                    % preferred location - this one should be bootstrapped too i suppose
                    
                    clear Average_FR_per_position_IN Average_FR_per_position_CH
                    for p=1:max(idx.pos)
                        Average_FR_per_position_IN(p)=nanmean(FR(idx.pt_in & idx.hands(:,hn) & ismember(epochs,s) & idx.pos==p));
                        Average_FR_per_position_CH(p)=nanmean(FR(idx.pt_ch & idx.hands(:,hn) & ismember(epochs,s) & idx.pos==p));
                    end
                    if takemin
                        [~,RF_position_index_IN]=min(Average_FR_per_position_IN);
                        [~,RF_position_index_CH]=min(Average_FR_per_position_CH);
                    else
                        [~,RF_position_index_IN]=max(Average_FR_per_position_IN);
                        [~,RF_position_index_CH]=max(Average_FR_per_position_CH);
                    end
                    
                    idx_IN_prefPI_in  =idx.pt_in & idx.hands(:,hn) & ismember(epochs,s) & idx.pos==RF_position_index_IN;
                    idx_IN_prefPO_in  =idx.pt_in & idx.hands(:,hn) & ismember(epochs,s) & idx.opp==RF_position_index_IN;
                    idx_CH_prefPI_in  =idx.pt_ch & idx.hands(:,hn) & ismember(epochs,s) & idx.pos==RF_position_index_IN;
                    idx_CH_prefPO_in  =idx.pt_ch & idx.hands(:,hn) & ismember(epochs,s) & idx.opp==RF_position_index_IN;
                    
                    idxb_IN_prefPI_in  =idx.pt_in & idx.hands(:,hn) & ismember(epochs,b) & idx.pos==RF_position_index_IN;
                    idxb_IN_prefPO_in  =idx.pt_in & idx.hands(:,hn) & ismember(epochs,b) & idx.opp==RF_position_index_IN;
                    idxb_CH_prefPI_in  =idx.pt_ch & idx.hands(:,hn) & ismember(epochs,b) & idx.pos==RF_position_index_IN;
                    idxb_CH_prefPO_in  =idx.pt_ch & idx.hands(:,hn) & ismember(epochs,b) & idx.opp==RF_position_index_IN;
                    
                    if sum(idx_IN_prefPI_in)>=keys.cal.min_trials_pref && sum(idx_CH_prefPI_in)>=keys.cal.min_trials_pref
                        [h,p,n]=do_stats(FR(idx_IN_prefPI_in),FR(idx_CH_prefPI_in),keys,0);
                        DF=nanmean(FR(idx_CH_prefPI_in))-nanmean(FR(idx_IN_prefPI_in));
                        labelindex=h*sign(DF)+2; labelindex(isnan(labelindex))=2;
                        anova_struct.([IN PT '_' LHRH '_' s{:}  '_prefP'])      =labels.choices{labelindex};
                        anova_struct.([CH PT '_' LHRH '_' s{:}  '_prefP'])      =labels.choices{labelindex};
                        anova_struct.([IN PT '_' LHRH '_' s{:}  '_prefP_PV'])   =p*sign(DF);
                        anova_struct.([CH PT '_' LHRH '_' s{:}  '_prefP_PV'])   =p*sign(DF);
                        anova_struct=add_normality_results(anova_struct,[IN PT '_' LHRH '_' s{:}  '_prefP'],keys,n);
                        anova_struct=add_normality_results(anova_struct,[CH PT  '_' LHRH '_' s{:}  '_prefP'],keys,n);
                    end
                    
                    if sum(idx_IN_prefPI_in)>=keys.cal.min_trials_pref && sum(idx_CH_prefPI_in)>=keys.cal.min_trials_pref &&...
                            sum(idx_IN_prefPO_in)>=keys.cal.min_trials_pref && sum(idx_CH_prefPO_in)>=keys.cal.min_trials_pref
                        anova_struct.([IN PT '_' LHRH '_' s{:}  '_prefP_DF'])   =nanmean(FR(idx_IN_prefPI_in))-nanmean(FR(idx_IN_prefPO_in));
                        anova_struct.([CH PT  '_' LHRH '_' s{:}  '_prefP_DF'])   =nanmean(FR(idx_CH_prefPI_in))-nanmean(FR(idx_CH_prefPO_in));
                    end
                    if sum(idx_IN_prefPI_in)>=keys.cal.min_trials_pref
                        anova_struct.([IN PT '_' LHRH '_' s{:}  '_prefPI_FR'])  =nanmean(FR(idx_IN_prefPI_in));
                        anova_struct.([IN PT '_' LHRH '_' s{:}  '_prefPI_DF'])  =nanmean(FR(idx_IN_prefPI_in))-nanmean(FR(idxb_IN_prefPI_in));
                    end
                    if sum(idx_IN_prefPO_in)>=keys.cal.min_trials_pref
                        anova_struct.([IN PT '_' LHRH '_' s{:}  '_prefPO_FR'])  =nanmean(FR(idx_IN_prefPO_in));
                        anova_struct.([IN PT '_' LHRH '_' s{:}  '_prefPO_DF'])  =nanmean(FR(idx_IN_prefPO_in))-nanmean(FR(idxb_IN_prefPO_in));
                    end
                    if sum(idx_CH_prefPI_in)>=keys.cal.min_trials_pref
                        anova_struct.([CH PT  '_' LHRH '_' s{:}  '_prefPI_FR'])  =nanmean(FR(idx_CH_prefPI_in));
                        anova_struct.([CH PT  '_' LHRH '_' s{:}  '_prefPI_DF'])  =nanmean(FR(idx_CH_prefPI_in))-nanmean(FR(idxb_CH_prefPI_in));
                    end
                    if sum(idx_CH_prefPO_in)>=keys.cal.min_trials_pref
                        anova_struct.([CH PT  '_' LHRH '_' s{:}  '_prefPO_FR'])  =nanmean(FR(idx_CH_prefPO_in));
                        anova_struct.([CH PT  '_' LHRH '_' s{:}  '_prefPO_DF'])  =nanmean(FR(idx_CH_prefPO_in))-nanmean(FR(idxb_CH_prefPO_in));
                    end
                    
                    % preferred location
                    idx_IN_prefPI_ch  =idx.pt_in & idx.hands(:,hn) & ismember(epochs,s) & idx.pos==RF_position_index_CH;
                    idx_IN_prefPO_ch  =idx.pt_in & idx.hands(:,hn) & ismember(epochs,s) & idx.opp==RF_position_index_CH;
                    idx_CH_prefPI_ch  =idx.pt_ch & idx.hands(:,hn) & ismember(epochs,s) & idx.pos==RF_position_index_CH;
                    idx_CH_prefPO_ch  =idx.pt_ch & idx.hands(:,hn) & ismember(epochs,s) & idx.opp==RF_position_index_CH;
                    
                    idxb_IN_prefPI_ch  =idx.pt_in & idx.hands(:,hn) & ismember(epochs,b) & idx.pos==RF_position_index_CH;
                    idxb_IN_prefPO_ch  =idx.pt_in & idx.hands(:,hn) & ismember(epochs,b) & idx.opp==RF_position_index_CH;
                    idxb_CH_prefPI_ch  =idx.pt_ch & idx.hands(:,hn) & ismember(epochs,b) & idx.pos==RF_position_index_CH;
                    idxb_CH_prefPO_ch  =idx.pt_ch & idx.hands(:,hn) & ismember(epochs,b) & idx.opp==RF_position_index_CH;
                    
                    if sum(idx_IN_prefPI_ch)>=keys.cal.min_trials_pref && sum(idx_CH_prefPI_ch)>=keys.cal.min_trials_pref
                        [h,p,n]=do_stats(FR(idx_IN_prefPI_ch),FR(idx_CH_prefPI_ch),keys,0);
                        DF=nanmean(FR(idx_CH_prefPI_ch))-nanmean(FR(idx_IN_prefPI_ch));
                        labelindex=h*sign(DF)+2; labelindex(isnan(labelindex))=2;
                        anova_struct.([IN PT '_' LHRH '_' s{:}  '_prefPch'])      =labels.choices{labelindex};
                        anova_struct.([CH PT  '_' LHRH '_' s{:}  '_prefPch'])      =labels.choices{labelindex};
                        anova_struct.([IN PT '_' LHRH '_' s{:}  '_prefPch_PV'])   =p*sign(DF);
                        anova_struct.([CH PT  '_' LHRH '_' s{:}  '_prefPch_PV'])   =p*sign(DF);
                        anova_struct=add_normality_results(anova_struct,[IN PT '_' LHRH '_' s{:}  '_prefPch'],keys,n);
                        anova_struct=add_normality_results(anova_struct,[CH PT  '_' LHRH '_' s{:}  '_prefPch'],keys,n);
                    end
                    
                    if sum(idx_IN_prefPI_ch)>=keys.cal.min_trials_pref && sum(idx_CH_prefPI_ch)>=keys.cal.min_trials_pref &&...
                            sum(idx_IN_prefPO_ch)>=keys.cal.min_trials_pref && sum(idx_CH_prefPO_ch)>=keys.cal.min_trials_pref
                        anova_struct.([IN PT '_' LHRH '_' s{:}  '_prefPch_DF'])   =nanmean(FR(idx_IN_prefPI_ch))-nanmean(FR(idx_IN_prefPO_ch));
                        anova_struct.([CH PT  '_' LHRH '_' s{:}  '_prefPch_DF'])   =nanmean(FR(idx_CH_prefPI_ch))-nanmean(FR(idx_CH_prefPO_ch));
                    end
                    if sum(idx_IN_prefPI_ch)>=keys.cal.min_trials_pref
                        anova_struct.([IN PT '_' LHRH '_' s{:}  '_prefPIch_FR'])   =nanmean(FR(idx_IN_prefPI_ch));
                        anova_struct.([IN PT '_' LHRH '_' s{:}  '_prefPIch_DF'])   =nanmean(FR(idx_IN_prefPI_ch))-nanmean(FR(idxb_IN_prefPI_ch));
                    end
                    if sum(idx_IN_prefPO_ch)>=keys.cal.min_trials_pref
                        anova_struct.([IN PT '_' LHRH '_' s{:}  '_prefPOch_FR'])   =nanmean(FR(idx_IN_prefPO_ch));
                        anova_struct.([IN PT '_' LHRH '_' s{:}  '_prefPOch_DF'])   =nanmean(FR(idx_IN_prefPO_ch))-nanmean(FR(idxb_IN_prefPO_ch));
                    end
                    if sum(idx_CH_prefPI_ch)>=keys.cal.min_trials_pref
                        anova_struct.([CH PT  '_' LHRH '_' s{:}  '_prefPIch_FR'])   =nanmean(FR(idx_CH_prefPI_ch));
                        anova_struct.([CH PT  '_' LHRH '_' s{:}  '_prefPIch_DF'])   =nanmean(FR(idx_CH_prefPI_ch))-nanmean(FR(idxb_CH_prefPI_ch));
                    end
                    if sum(idx_CH_prefPO_ch)>=keys.cal.min_trials_pref
                        anova_struct.([CH PT  '_' LHRH '_' s{:}  '_prefPOch_FR'])   =nanmean(FR(idx_CH_prefPO_ch));
                        anova_struct.([CH PT  '_' LHRH '_' s{:}  '_prefPOch_DF'])   =nanmean(FR(idx_CH_prefPO_ch))-nanmean(FR(idxb_CH_prefPO_ch));
                    end
                end
            end
        end
        
        
    end
    
    %% perturbation ANOVA (!)
    % temporarly recompute this stuff    :((
    tr=idx.(INCH) & ismember(epochs,epoch_multicomp) & (idx.PT==0 | idx.PT==1); %& ismember(epochs,epochs_for_multicomparison);
    %tr_main=idx.(INCH) & ismember(epochs,keys.anova_epochs.main) & (idx.PT==0 | idx.PT==1); %& ismember(epochs,epochs_for_multicomparison);
    for fn=idx_fieldnames
        idx.(['tr_' fn{:}])=      idx.(fn{:}) & repmat(tr,1,size(idx.(fn{:}),2)) ;
    end
    idx.tr_CT=idx.tr_PT==0;
    
    perturbations_to_compare=unique(idx.tr_PT);
    if numel(perturbations_to_compare)>1
        %labels.control_test={'EP','BL'};
        labels.PT='';
        labels.conditions={'LH_LS','LH_RS','RH_LS','RH_RS','LH','RH','LS','RS'};
        idx.conditions=[idx.tr_LH_LS idx.tr_LH_RS idx.tr_RH_LS idx.tr_RH_RS idx.tr_LS idx.tr_RS];
        anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,'SxH','tr_CT','tr_PT','perturbation');
        
        %         for c=1:numel(labels.conditions)
        %             idx_con=idx.(['tr_' labels.conditions{c}]);
        %             CON=labels.conditions{c};
        %             if any(idx.tr_PT(tr_main & idx_con)==0) && any(idx.tr_PT(tr_main & idx_con)==1)
        %                 %[anova_out.p,anova_out.table,anova_out.stats,anova_out.terms] = anovan(FR(tr_main & idx_con),[epoch_idx(tr_main & idx_con) idx.tr_PT(tr_main & idx_con)],'model','full','varnames',{'epoch', 'perturbation'},'display',keys.plot.anova_tables);
        %                 m=FR(tr_main & idx_con);
        %                 group={epoch_idx(tr_main & idx_con) idx.tr_PT(tr_main & idx_con)};
        %                 pval = ph_anova2(m,group,1,keys);
        %                 h1=pval(1)<0.05;
        %                 h2=pval(2)<0.05;
        %                 h3=pval(3)<0.05;
        %                 anova_struct.([INCH '_' CON '_epoch_main'])    =labels.true{h1+2};%main effect on epoch!
        %                 anova_struct.([INCH '_' CON '_PT_main'])       =labels.true{h2+2}; %main effect on perturbation!
        %                 anova_struct.([INCH '_' CON '_ExP'])           =labels.true{h3+2}; %epoch perturbation interaction!
        %             end
        %         end
    end
end
end

function anova_struct=effector_comparison_anova(keys,epochs,FR,idx,Positions,Fixations)
anova_struct=struct();
conditions=keys.TTconditions;
conditions.choice=conditions.choice([sum(idx.in)>0 sum(idx.ch)>0]);
labels=keys.TTlabels;
conditions.perturbation=conditions.perturbation([sum(idx.PT==0)>0 sum(idx.PT==1)>0]);
pt_values=[0,1];
pt_values=pt_values([sum(idx.PT==0)>0 sum(idx.PT==1)>0]);

%% Independently for Instructed and Choice !!
idx_fieldnames=fieldnames(idx)';
for ch=1:numel(conditions.choice)
    
    for pt=1:numel(pt_values)
        pt_value=pt_values(pt);
        labels.PT=conditions.perturbation{pt};
        %% this part should go outside the loop! ... NAH
        tr=idx.(conditions.choice{ch}) &  idx.PT==pt_value;
        labels.INCH=conditions.choice{ch};
        for fn=idx_fieldnames
            idx.(['tr_' fn{:}])=      idx.(fn{:}) & repmat(tr,1,size(idx.(fn{:}),2));
        end
        
        labels.control_test={};
        labels.PT={''};
        labels.conditions={''};
        idx.conditions=idx.tr_all;
        % epoch_SxH_multicomp needs to be replaced by something meaningful
        anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,'SxH','tr_E1','tr_E2','effector');
        
        labels.conditions={'LH_LS','LH_RS','RH_LS','RH_RS','LS','RS'};
        idx.conditions=[idx.tr_LH_LS idx.tr_LH_RS idx.tr_RH_LS idx.tr_RH_RS idx.tr_LS idx.tr_RS];
        anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,'SxH','tr_E1','tr_E2','effector');
    end
end
end

%% stats per epoch and condition
function anova_struct=do_stats_per_epoch_and_condition(keys,anova_struct,FR,epochs,idx,labels,epochs_FN,idx1_FN,idx2_FN,fieldnamepart)
multicomp_epochs=keys.anova_epochs.(epochs_FN);
multicomp_epochs=multicomp_epochs(ismember(multicomp_epochs,epochs'));
label=labels.(fieldnamepart);
idx1= idx.(idx1_FN);
idx2= idx.(idx2_FN);
INCH= labels.INCH;
PT= labels.PT;
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
            prefix=[INCH PT '_' COND{c} SEP s{:} '_' fieldnamepart ];
            anova_struct.(prefix)=label{labelindex};
            anova_struct.([prefix '_DF'])=DF;
            anova_struct.([prefix '_PV'])=p*sign(DF);
            anova_struct.([prefix '_IX'])=DF/(nanmean(FR(idx2 & idxS))+ nanmean(FR(idx1 & idxS)));
            anova_struct=add_normality_results(anova_struct,prefix,keys,n);
            %             if ~isempty(labels.control_test)
            %                 anova_struct.([INCH '_' COND{c} SEP s{:} '_' labels.control_test{1} '_FR']) = nanmean(FR(idx1));
            %                 anova_struct.([INCH '_' COND{c} SEP s{:} '_' labels.control_test{2} '_FR']) = nanmean(FR(idx2));
            %             end
        end
    end
end
end

%% similarly, make epoch effects
function anova_struct=do_epoch_stats(keys,anova_struct,FR,epochs,idx,labels,epochs_FN,fieldnamepart)
%% epoch tuning per hand
multicomp_epochs=keys.anova_epochs.(epochs_FN);
idx_ep=ismember(multicomp_epochs(:,1),epochs') &  ismember(multicomp_epochs(:,2),epochs');
multicomp_epochs=multicomp_epochs(idx_ep,:);
label=labels.(fieldnamepart);
INCH= labels.INCH;
PT= labels.PT;
COND=labels.conditions;
for row=1:size(multicomp_epochs,1)
    s=multicomp_epochs(row,2);
    b=multicomp_epochs(row,1);
    for  c=1:size(idx.conditions,2)
        idx1=ismember(epochs,b) & idx.conditions(:,c);
        idx2=ismember(epochs,s) & idx.conditions(:,c);
        if sum(idx1)>1 && sum(idx2)>1
            [h,p,n]=do_stats(FR(idx1),FR(idx2),keys,1); %%unpaired ??
        else
            h=false; n=NaN; p=NaN;
        end
        DF=(nanmean(FR(idx2))-nanmean(FR(idx1)));
        labelindex=h*sign(DF)+2; labelindex(isnan(labelindex))=2;
        if isempty(COND{c})
            prefix=[INCH PT '_' s{:} '_' fieldnamepart];
        else
            prefix=[INCH PT '_' COND{c} '_' s{:} '_' fieldnamepart];
        end
        anova_struct.([prefix '_FR'])= nanmean(FR(idx2));
        anova_struct.([prefix '_DF'])= DF;
        anova_struct.([prefix '_PV'])= p*sign(DF);
        anova_struct.(prefix)        = label{labelindex};
        anova_struct=add_normality_results(anova_struct,prefix,keys,n);
    end
end
end

%% embedded test selection
function [h,p,n]=do_stats(A,B,keys,paired)
p=1;
switch keys.AN.test_types
    case 'parametric'
        if paired
            if sum(~isnan(A)&~isnan(B))>1
                [h,p] = ttest(A,B);
            else
                h=0;
            end
        else
            if sum(~isnan(A))>1 && sum(~isnan(B))>1
                [h,p] =  ttest2(A,B);
            else
                h=0;
            end
        end
    case 'nonparametric'
        if paired
            if sum(~isnan(A)&~isnan(B))>1
                [p, h] = signrank(A,B);
            else
                h=0;
            end
        else
            if sum(~isnan(A))>1 && sum(~isnan(B))>1
                [p, h] = ranksum(A,B);
            else
                h=0;
            end
        end
    case 'permutations'
        %     [p] = permtest(A,B,1000,'conservative');
        %     h=p<0.05;
        
        
        if paired %% how would paired here even look like?
            % B=B-A; A=zeros(size(A));
        else
        end
        
        
        [h,p]=DAG_permutation_test(A,B,keys.AN.n_permutations);
end
p=single(p);
if keys.AN.check_normality && numel(A)>4 && numel(B)>4
    %[H,pValue,KSstatistic,criticalValue] = lillietest(x,alpha,distr,mctol)
    n=lillietest(A)==0 && lillietest(B)==0;
else
    n=1;
end
end

function p=ph_anova2(FR,idx,interaction,keys)
switch keys.AN.test_types
    case 'parametric'
        p=anovan(FR,idx,'model','full','varnames',{'x','y'},'display',keys.plot.anova_tables);
    case 'nonparametric'
        p=anovan(FR,idx,'model','full','varnames',{'x','y'},'display',keys.plot.anova_tables);
    case 'permutations'
        p=randanova2(FR,idx,keys.AN.n_permutations_randanova2,interaction);
end
end

function p=ph_anova1(FR,idx,keys)
switch keys.AN.test_types
    case 'parametric'
        p=anovan(FR,idx,'model','full','varnames',{'x'},'display',keys.plot.anova_tables);
    case 'nonparametric'
        p=anovan(FR,idx,'model','full','varnames',{'x'},'display',keys.plot.anova_tables);
    case 'permutations'
        p=randanova1(FR,idx,keys.AN.n_permutations_randanova1);
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

