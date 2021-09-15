function tuning_per_unit_table=ph_ANOVAS(population,keys)
tuning_per_unit_table=keys.tuning_per_unit_table;
for unit=1:numel(population)
    anova_struct_current_unit=struct();
    unit_ID=population(unit).unit_ID;
    types           =[population(unit).trial.type];
    effectors       =[population(unit).trial.effector];
    type_effectors  =combvec(unique(types),unique(effectors))';
    control    =ismember([population(unit).trial.perturbation],keys.cal.perturbation_groups{1});
    
    %% accepted AND completed!!
    tr_considered=[population(unit).trial.accepted]==1;
    if keys.cal.completed
        tr_considered = tr_considered & [population(unit).trial.completed]==1;
    end
    if sum(tr_considered)>0
        for a=1:numel(keys.position_and_plotting_arrangements) % arrangement defines poaitions, therefore alsohemifield (which is part of conditions)
            keys.arrangement=keys.position_and_plotting_arrangements{a};
            o=ph_arrange_positions_and_plots(keys,population(unit).trial(tr_considered)); % arrangement            
            keys.normalization_field='AN';
            o=ph_condition_normalization(o,keys); %condition wise normalization (also reduces conditions!???)

            for type=unique(types)
                %% check carefully multicomp epochs !!
                keys.main_multicomp                     =keys.ANOVAS_PER_TYPE(type).main;
                keys.epoch_multicomp                    =keys.ANOVAS_PER_TYPE(type).epoch;
                keys.epoch_spaceLR_multicomp            =keys.ANOVAS_PER_TYPE(type).spaceLR;
                keys.epoch_hands_multicomp              =keys.ANOVAS_PER_TYPE(type).hands;
                keys.epoch_SxH_multicomp                =keys.ANOVAS_PER_TYPE(type).SxH;
                keys.epoch_position_multicomp           =keys.ANOVAS_PER_TYPE(type).positions;
                
                for effector=unique(effectors)
                    keys=ph_get_epoch_keys(keys,type,effector,sum(type_effectors(:,1)==type)>1);
                    [~, condition_fieldname_part]=MPA_get_type_effector_name(type,effector);
                    tr_index= [o.trial.effector]==effector & [o.trial.type]==type ;% & ismember(hands,keys.cal.reach_hand);
                                       
                    if sum(tr_index)==0;
                        disp(sprintf('no trials for effector %.0f type %.0f hands %s completed= %.0f',effector,type,mat2str(keys.cal.reach_hand),keys.cal.completed));
                        continue;
                    end
                    o_e=o;
                    o_e.trial=o.trial(tr_index);
                    
                    trial_criterion=get_minimum_trials(keys,o_e,tr_index);
                    [FR,epochs,idx,u_pos,u_fix]=ph_get_anova_inputs(o_e,keys);
                    [anova_struct]=n_way_anova(keys,epochs,FR,idx,u_pos,u_fix,unit,unit_ID);
                    anova_struct_current_unit.([condition_fieldname_part '_'  keys.arrangement(1:3)])=anova_struct; clear anova_struct
                    FN_crit=fieldnames(trial_criterion);
                    for f=1:numel(FN_crit)
                        anova_struct_current_unit.([condition_fieldname_part '_'  keys.arrangement(1:3)]).(FN_crit{f})= trial_criterion.(FN_crit{f});
                    end
                end
                
                %% effector comparison !
                effectors_for_this_type=type_effectors(type_effectors(:,1)==type,2)';
                effectors_to_compare=combvec(effectors_for_this_type,effectors_for_this_type);
                effectors_to_compare=effectors_to_compare(:,effectors_to_compare(1,:)<effectors_to_compare(2,:));
                for comp=1:size(effectors_to_compare,2)
                    comp_eff=effectors_to_compare(:,comp);
                    keys=ph_get_epoch_keys(keys,type,comp_eff,1);
                    [~, condition_fieldname_part1]=MPA_get_type_effector_name(type,comp_eff(1));
                    [~, condition_fieldname_part2]=MPA_get_type_effector_name(type,comp_eff(2));
                    condition_fieldname_part=[condition_fieldname_part1 '_vs_' condition_fieldname_part2];
                    tr_index= ismember(effectors,comp_eff) & types == type & control;
                    tr_index= tr_index(tr_considered); %& ismember(hands,keys.cal.reach_hand);
                    
                    if sum(tr_index)==0;
                        disp(sprintf('no trials for effectors %s type %.0f hands %s completed= %.0f',mat2str(comp_eff'),type,mat2str(keys.cal.reach_hand),keys.cal.completed));
                        continue;
                    end
                    
                    o_e=o;
                    o_e.trial=o.trial(tr_index);
                    keys.labels.eff={condition_fieldname_part1 '-' condition_fieldname_part2};
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
    %inital_fieldnames={'unit_ID','monkey','target','perturbation_site','grid_x','grid_y','electrode_depth','FR','stability_rating','SNR_rating','Single_rating','waveform_width'};
    inital_fieldnames={'unit_ID','monkey','target','grid_x','grid_y','electrode_depth','FR','stability_rating','SNR_rating','Single_rating','waveform_width'};
    unit_table(1,1:numel(inital_fieldnames))=inital_fieldnames;
    for fn=1:numel(inital_fieldnames)
        unit_table{rows_to_update,fn}=population(unit).(inital_fieldnames{fn});
    end
    title_counter=numel(inital_fieldnames);
    FN=fieldnames(anova_struct_current_unit);
    for fn=1:numel(FN)
        FNsub=fieldnames(anova_struct_current_unit.(FN{fn}));
        for fnsub=1:numel(FNsub)
            title_counter=title_counter+1;
            unit_table{1,title_counter}=[FNsub{fnsub} '_' FN{fn}];
            unit_table{rows_to_update,title_counter}=anova_struct_current_unit.(FN{fn}).(FNsub{fnsub});
        end
    end
    tuning_per_unit_table=DAG_update_mastertable_cell(tuning_per_unit_table,unit_table,rows_to_update);
end
end

function anova_struct=n_way_anova(keys,epochs,FR,idx,Positions,Fixations,u,unit_ID)
INCHnamepart=keys.labels.choices([sum(idx.in)>0 sum(idx.ch)>0]);
LHRHnamepart=keys.labels.handsLR;
LSRSnamepart={'LS','RS'};
labelen={'su','-','en','bi'};
labelsp={'LS','-','RS'};
labelud={'DN','-','UP'};
labelha={'LH','-','RH'};
labelcr={'UC','-','CR'};
labelIN={'in','-','ch'};
true_labels={'false','true'};

%% doesnt even belong here...
%epochs_for_multicomparison=keys.epoch_for_multicomparison(ismember(keys.epoch_for_multicomparison,keys.EPOCHS(:,1)));
idx_ep=ismember(keys.epoch_multicomp(:,1),keys.EPOCHS(:,1)) &  ismember(keys.epoch_multicomp(:,2),keys.EPOCHS(:,1));
epoch_multicomp=keys.epoch_multicomp(idx_ep,:);
idx_ep=ismember(keys.main_multicomp,keys.EPOCHS(:,1));
main_multicomp=keys.main_multicomp(idx_ep,:);

[~, ~, epoch_idx]=unique(epochs);
idx_fieldnames=fieldnames(idx)';

%% Independently for Instructed and Choice !!
for ch=1:numel(INCHnamepart)
    %% this part should go NOT? outside the loop! % perturbation here is defined as always the lowest accessible for control...;
    tr_ch=idx.(INCHnamepart{ch}) & ismember(epochs,epoch_multicomp) & idx.PT==0;
    tr_main=idx.(INCHnamepart{ch}) & ismember(epochs,main_multicomp) & idx.PT==0; %& ismember(epochs,epochs_for_multicomparison);
    for fn=idx_fieldnames
        idx.(['tr_' fn{:}])=      idx.(fn{:}) & repmat(tr_ch,1,size(idx.(fn{:}),2)) ;
    end
    
    %% number of trials for each condition to decide about type of (main) anova
    n_lr=sum(idx.tr_LS)>0 && sum(idx.tr_RS)>0;
    n_2hands=sum(idx.tr_LH_LS)>0 && sum(idx.tr_LH_RS)>0 && sum(idx.tr_RH_LS)>0 && sum(idx.tr_RH_RS)>0;
    if n_2hands % three-way anova epoch, space, hand
        keys.anova_varnames={'epoch' 'spaceLR' 'hands'};
        tr=tr_ch & (idx.tr_LH_LS | idx.tr_LH_RS | idx.tr_RH_LS | idx.tr_RH_RS);
        Par=[epoch_idx, idx.tr_RS, idx.tr_RH];
    elseif n_lr % two-way anova epoch, space
        keys.anova_varnames={'epoch' 'spaceLR'};
        tr=tr_ch & (idx.tr_LS | idx.tr_RS);
        Par=[epoch_idx, idx.tr_RS];
    else
        disp(['Not fitting dataset... no left or right ' INCHnamepart{ch} ' present?'])
        anova_struct.([INCHnamepart{ch} '_epoch_main'])=[];
        continue;
    end
    
    %% main effects and general interactions
    [anova_out.p,anova_out.table,anova_out.stats,anova_out.terms] = anovan(FR(tr_main),Par(tr_main,:),'model','full','varnames',keys.anova_varnames,'display',keys.plot.anova_tables);
    anova_struct.([INCHnamepart{ch} '_epoch_main'])=anova_out.p(1)<0.05; %main effect on epoch!
    handindexes=[1,2,3];
    handindexes=handindexes(any(idx.tr_hands));
    
    %% epoch tuning per hand
    multicomp_epochs=epoch_multicomp;
    for row=1:size(multicomp_epochs,1)
        s=multicomp_epochs(row,2);
        b=multicomp_epochs(row,1);
        % per hand
        for hn=handindexes
            % separately define enhancement and supression for each side (left/right)
            for sideindex=1:2
                idx1=tr & ismember(epochs,b) & idx.sides(:,sideindex) & idx.tr_hands(:,hn);
                idx2=tr & ismember(epochs,s) & idx.sides(:,sideindex) & idx.tr_hands(:,hn);
                if sum(idx1)>1 && sum(idx2)>1
                    h=do_stats(FR(idx1),FR(idx2),keys,0); %%unpaired ?? ttest versus baseline
                else
                    h=false;
                end
                DF(sideindex)=(nanmean(FR(idx2))-nanmean(FR(idx1)));
                labelindex=h*sign(DF(sideindex))+2; labelindex(isnan(labelindex))=2;
                anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' LSRSnamepart{sideindex} '_' s{:} '_epoch_FR'])=nanmean(FR(idx2));
                anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' LSRSnamepart{sideindex} '_' s{:} '_epoch_DF'])= DF(sideindex) ;
                anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' LSRSnamepart{sideindex} '_' s{:} '_epoch'])= labelen{labelindex};
                ensu(sideindex)=h*sign(DF(sideindex))+2;
            end
            % select label dependend on both sides
            ensu(isnan(ensu))=2;
            if any(ensu==1) && any(ensu==3)
                labelindex=4;
            elseif any(ensu~=2)
                labelindex=unique(ensu(ensu~=2));
            else
                labelindex=2;
            end
            [~,DFidx]=max(abs(DF));
            DF=DF(DFidx);
            idxbl=tr & ismember(epochs,b) & idx.sides(:,DFidx) & idx.tr_hands(:,hn); % take respective baseline index for side!
            idxep=tr & ismember(epochs,s) & idx.sides(:,DFidx) & idx.tr_hands(:,hn); % take respective epoch index for side!
            anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn}  '_' s{:} '_epoch'])=labelen{labelindex};
            anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn}  '_' s{:} '_epoch_DF'])=DF; %%% to think about !!!!
            anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn}  '_' s{:} '_epoch_SC'])=DF/nanmean(FR(idxbl)); %%% to think about !!!! now it corresponds to % signal change
            labelindex=-1*all(ensu==1)+all(ensu==3)+2;
            anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn}  '_' s{:} '_epoch_bilateral'])=labelen{labelindex};
            
        end
    end
    
    %% epoch tuning across all trials
    multicomp_epochs=epoch_multicomp;
    for row=1:size(multicomp_epochs,1)
        s=multicomp_epochs(row,2);
        b=multicomp_epochs(row,1);
        idxbl=tr_ch & ismember(epochs,b);
        idxep=tr_ch & ismember(epochs,s);        
        if sum(idxbl)>1 && sum(idxep)>1
            h=do_stats(FR(idxbl),FR(idxep),keys,1); %%paired ?? ttest versus baseline
        else
            h=false;
        end
        DF=(nanmean(FR(idxep))-nanmean(FR(idxbl)));
        labelindex=h*sign(DF)+2; labelindex(isnan(labelindex))=2;
        anova_struct.([INCHnamepart{ch} '_' s{:} '_epoch'])=labelen{labelindex};
        %             anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn}  '_' s{:} '_epoch_DF'])=DF; %%% to think about !!!!
        %             anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn}  '_' s{:} '_epoch_SC'])=DF/nanmean(FR(idxbl)); %%% to think about !!!! now it corresponds to % signal change
        
    end
    
    
    %% hand OR space tuning (gets overwritten in case of both)
    for k=2:numel(keys.anova_varnames)
        multicomp_epochs=keys.(['epoch_' keys.anova_varnames{k} '_multicomp']);
        multicomp_epochs=multicomp_epochs(ismember(multicomp_epochs,epochs')); %% to hink about: add &  ismember(keys.epoch_multicomp(:,2),keys.EPOCHS(:,1));
        switch keys.anova_varnames{k}
            case 'spaceLR'
                label={'LS','-','RS'};
            case 'hands'
                label={'LH','-','RH'};
        end
        
        idx1=tr & Par(:,k)==0;
        idx2=tr & Par(:,k)==1;
        idx1m=tr_main & Par(:,k)==0 & ismember(epochs,multicomp_epochs);
        idx2m=tr_main & Par(:,k)==1 & ismember(epochs,multicomp_epochs);
        
        labelindex=(anova_out.p(k)<0.05) *sign(nanmean(FR(idx2m))-nanmean(FR(idx1m)))+2; labelindex(isnan(labelindex))=2;
        anova_struct.([INCHnamepart{ch} '_' keys.anova_varnames{k} '_main'])=label{labelindex}; %main effect with direction!
        anova_struct.([INCHnamepart{ch} '_Ex' upper(keys.anova_varnames{k}(1))])=anova_out.p(k+numel(keys.anova_varnames)-1)<0.05; %main effect with direction!
        
        for s=multicomp_epochs(:)'
            idxS=ismember(epochs,s);
            h=do_stats(FR(idx1 & idxS),FR(idx2 & idxS),keys,0);
            DF=nanmean(FR(idx2 & idxS))-nanmean(FR(idx1 & idxS));
            labelindex=h*sign(DF)+2; labelindex(isnan(labelindex))=2;
            anova_struct.([INCHnamepart{ch} '_' s{:} '_' keys.anova_varnames{k}])=label{labelindex};
            anova_struct.([INCHnamepart{ch} '_' s{:} '_' keys.anova_varnames{k} '_DF'])=DF;
            anova_struct.([INCHnamepart{ch} '_' s{:} '_' keys.anova_varnames{k} '_IX'])=DF/(nanmean(FR(idx2 & idxS))+ nanmean(FR(idx1 & idxS)));
        end
    end
    
    %% space x hand anovas
    if ismember('spaceLR',keys.anova_varnames) && ismember('hands',keys.anova_varnames)
        k=6;
        multicomp_epochs=keys.epoch_SxH_multicomp;
        multicomp_epochs=multicomp_epochs(ismember(multicomp_epochs,epochs'));
        tr_SXH=tr_main & ismember(epochs,multicomp_epochs);
        labelindexcr=(anova_out.p(k)<0.05) *sign(nanmean(FR(tr_SXH & idx.tr_CR))-nanmean(FR(tr_SXH & idx.tr_UC)))+2; labelindexcr(isnan(labelindexcr))=2;% 3 crossed>uncrossed ,1 crossed<uncrossed
        anova_struct.([INCHnamepart{ch} '_SxH'])=labelcr{labelindexcr};
        anova_struct.([INCHnamepart{ch} '_ExSxH'])=anova_out.p(7)<0.05;
        
        for s=multicomp_epochs(:)'
            idxS= tr & ismember(epochs,s);
            if any(~isnan(FR(idxS))) && any(idx.tr_LS(idxS)) && any(idx.tr_RS(idxS)) && any(idx.tr_LH(idxS)) && any(idx.tr_RH(idxS))
                [anova_outs.p,anova_outs.table,anova_outs.stats,anova_outs.terms] = anovan(FR(idxS),double([idx.tr_RS(idxS),idx.tr_RH(idxS)]),'model','full','varnames',{'space';'hands'},'display',keys.plot.anova_tables);
                
                %% replacing hand and space tuning with anova per epoch !!
                h=anova_outs.p(1)<0.05;
                DF=nanmean(FR(idx.tr_RS & idxS))-nanmean(FR(idx.tr_LS & idxS));
                labelindexsp=h*sign(DF)+2; labelindexsp(isnan(labelindexsp))=2;
                anova_struct.([INCHnamepart{ch} '_' s{:} '_SxH_spaceLR'])= labelsp{labelindexsp}; %
                anova_struct.([INCHnamepart{ch} '_' s{:} '_SxH_spaceLR_DF'])= DF; %
                anova_struct.([INCHnamepart{ch} '_' s{:} '_SxH_spaceLR_IX'])= DF/(nanmean(FR(idx.tr_RS & idxS))+nanmean(FR(idx.tr_LS & idxS))); %
                
                h=anova_outs.p(2)<0.05;
                DF=nanmean(FR(idx.tr_RH & idxS))-nanmean(FR(idx.tr_LH & idxS));
                labelindexha=h*sign(DF)+2; labelindexha(isnan(labelindexha))=2;
                anova_struct.([INCHnamepart{ch} '_' s{:} '_SxH_hands'])= labelha{labelindexha};
                anova_struct.([INCHnamepart{ch} '_' s{:} '_SxH_hands_DF'])= DF;
                anova_struct.([INCHnamepart{ch} '_' s{:} '_SxH_hands_IX'])= DF/(nanmean(FR(idx.tr_RH & idxS))+nanmean(FR(idx.tr_LH & idxS)));
                
                h=anova_outs.p(3)<0.05;
                DF=nanmean(FR(idx.tr_CR & idxS))-nanmean(FR(idx.tr_UC & idxS));
                labelindexcr=h*sign(DF)+2; labelindexcr(isnan(labelindexcr))=2;
                anova_struct.([INCHnamepart{ch} '_' s{:} '_SxH'])=labelcr{labelindexcr};
                anova_struct.([INCHnamepart{ch} '_' s{:} '_SxH_DF'])=DF;
                anova_struct.([INCHnamepart{ch} '_' s{:} '_SxH_IX'])=DF/(nanmean(FR(idx.tr_CR & idxS))+nanmean(FR(idx.tr_UC & idxS)));
                
                %% space per hand tuning and hand per space tuning
                
                for hn=handindexes
                    idxL=idxS & idx.tr_LS & idx.tr_hands(:,hn);
                    idxR=idxS & idx.tr_RS & idx.tr_hands(:,hn);
                    h=do_stats(FR(idxL),FR(idxR),keys,0); 
                    DF=(nanmean(FR(idxR))-nanmean(FR(idxL)));
                    labelindex=h*sign(DF)+2; labelindex(isnan(labelindex))=2;
                    anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_spaceLR_DF'])= DF;
                    anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_spaceLR'])= labelsp{labelindex};
                end
                
                for sideindex=1:2
                    idxL=idxS & idx.sides(:,sideindex) & idx.tr_LH;
                    idxR=idxS & idx.sides(:,sideindex) & idx.tr_RH;
                    h=do_stats(FR(idxL),FR(idxR),keys,0); 
                    DF=(nanmean(FR(idxR))-nanmean(FR(idxL)));
                    labelindex=h*sign(DF)+2; labelindex(isnan(labelindex))=2;
                    anova_struct.([INCHnamepart{ch} '_' LSRSnamepart{sideindex} '_' s{:} '_hands_DF'])= DF;
                    anova_struct.([INCHnamepart{ch} '_' LSRSnamepart{sideindex} '_' s{:} '_hands'])= labelha{labelindex};
                end
                
                
            end
        end
    end
    
    %% Position anova (per hand) + comparison CH vs IN
    epochs_for_position_comparison=keys.epoch_position_multicomp(ismember(keys.epoch_position_multicomp,keys.EPOCHS(:,1)));
    multicomp_epochs=epochs_for_position_comparison(ismember(epochs_for_position_comparison,epochs));
    
    u_hnd                       =find(any(idx.tr_hands,1));
    for hn=u_hnd
        tr=idx.(INCHnamepart{ch}) & idx.tr_hands(:,hn) & ismember(epochs,epochs_for_position_comparison);
        if size(Fixations,1)==1
            varnames={'State','position'};
            [anova_out.p,anova_out.table,anova_out.stats,anova_out.terms] =...
                anovan(FR(tr),{epochs(tr),idx.pos(tr)},'model','full','varnames',varnames,'display',keys.plot.anova_tables);
        else
            varnames={'State','position','fixation'};
            [anova_out.p,anova_out.table,anova_out.stats,anova_out.terms] =...
                anovan(FR(tr),{epochs(tr),idx.pos(tr),idx.fix(tr)},'model','full','varnames',varnames,'display',keys.plot.anova_tables);
            anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_fixation_main'])=anova_out.p(3)<0.05; %main effect on epoch!
        end
        anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_epoch_main'])=anova_out.p(1)<0.05; %main effect on epoch!
        anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_position_main'])=anova_out.p(2)<0.05; %main effect on epoch!
        
        
        for s=multicomp_epochs(:)'
            tr=idx.(INCHnamepart{ch}) & idx.tr_hands(:,hn) & ismember(epochs,s);  %this is not ideal, because tr is used in a different way here
            b=epoch_multicomp(ismember(epoch_multicomp(:,2),s),1);
            
            if any(~isnan(FR(tr)))
                [anova_outs.p] = anovan(FR(tr),[idx.pos(tr),idx.fix(tr)],'model','full','varnames',{'Positions','Fixations'},'display',keys.plot.anova_tables);
                if any(isnan(anova_outs.p)) %this is the case for target position by location in pulvinar eye gaze project AND for fixation only
                    % take only positions which have combinations of
                    % different fixations to compute interaction
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
                
                anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_position'])=true_labels{h1+1};
                anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_fixation'])=true_labels{h2+1};
                anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_PxF'])=true_labels{h3+1};


                pecc = anovan(FR(tr),[idx.ecc(tr),idx.ang(tr)],'model','full','varnames',{'Eccentricity','Angle'},'display',keys.plot.anova_tables);
                
                h1=pecc(1)<0.05;
                h2=pecc(2)<0.05;
                h3=pecc(3)<0.05;
                anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_distance'])=true_labels{h1+1};
                anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_angle'])=true_labels{h2+1};
                anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_DxA'])=true_labels{h3+1};

                % x and y separately
                pxy = anovan(FR(tr),[idx.pos_x(tr),idx.pos_y(tr)],'model','full','varnames',{'x','y'},'display',keys.plot.anova_tables);
                h1=pxy(1)<0.05;
                h2=pxy(2)<0.05;
                h3=pxy(3)<0.05;
                anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_positionx'])=true_labels{h1+1};
                anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_positiony'])=true_labels{h2+1};
                anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_positionxy'])=true_labels{h3+1};
                for xory={'x','y'}
                    xorytag=xory{:};
                    switch xorytag
                        case 'x'
                            labelxy=labelsp;
                            hxy=h1;
                        case 'y'
                            labelxy=labelud;
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
                            anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_gaze_modulation_' xorytag])='monotonous';
                            anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_gaze_pref_' xorytag])=labelxy{(directionality_per_position(end)==1)*2+1};
                            
                        else
                            anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_gaze_modulation_' xorytag])='nonmonotonous';
                            if pref_idx==1 || pref_idx==numel(u_xy_pos)
                                anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_gaze_pref_' xorytag])='PE';
                            else
                                anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_gaze_pref_' xorytag])='CE';
                            end
                        end
                    else
                        %% not even nonmonotoneous
                        anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_gaze_modulation_' xorytag])='-';
                        anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_gaze_pref_' xorytag])='-';
                    end
                end
                
                %% choice part for position with strongest response
                if ch==2  && ~isempty(anova_struct.in_epoch_main) && isfield(anova_struct,['in_' LHRHnamepart{hn} '_' s{:} '_epoch_DF'])%any(idx.in) %% there has to be choice and instructed
                    clear Average_FR_per_position_IN Average_FR_per_position_CH
                    for p=1:max(idx.pos)
                        Average_FR_per_position_IN(p)=nanmean(FR(idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.pos==p));
                        Average_FR_per_position_CH(p)=nanmean(FR(idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.pos==p));
                    end
                    % preference (highest or lowest) defined by over all
                    % increase or decrease in this epoch
                    if false % anova_struct.(['in_' LHRHnamepart{hn} '_' s{:} '_epoch_DF']) <0
                        invertsign=-1;
                        takemin=1;
                    else
                        invertsign=1;
                        takemin=0;
                    end
                    
                    %% preference based on instructed trials
                    RF_IN_RS_IDX=idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.RS;
                    RF_IN_LS_IDX=idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.LS;
                    RF_IN_RS_IDX(randsample(find(RF_IN_RS_IDX),round(sum(RF_IN_RS_IDX)/2)))=false;  % taking only 50% of trials for preference estimation
                    RF_IN_LS_IDX(randsample(find(RF_IN_LS_IDX),round(sum(RF_IN_LS_IDX)/2)))=false;
                    RF_TEST=~RF_IN_RS_IDX & ~RF_IN_LS_IDX;
                    RF_in_hemifield_index_IN =invertsign*sign(nanmean(FR(RF_IN_RS_IDX)) - nanmean(FR(RF_IN_LS_IDX))) +2;
                    RF_out_hemifield_index_IN=invertsign*sign(nanmean(FR(RF_IN_RS_IDX)) - nanmean(FR(RF_IN_LS_IDX)))*-1 +2;
                    
                     %% preference based on choice trials (??)
                    RF_CH_RS_IDX=idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.RS;
                    RF_CH_LS_IDX=idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.LS;            % here was a bug i believe, taking instructed trials for left space???
                    RF_in_hemifield_index_CH =invertsign*sign(nanmean(FR(RF_CH_RS_IDX)) - nanmean(FR(RF_CH_LS_IDX))) +2;
                    RF_out_hemifield_index_CH=invertsign*sign(nanmean(FR(RF_CH_RS_IDX)) - nanmean(FR(RF_CH_LS_IDX)))*-1 +2;
                    
                    RF_in_hemifield_index_IN(isnan(RF_in_hemifield_index_IN))=2;
                    RF_out_hemifield_index_IN(isnan(RF_out_hemifield_index_IN))=2;
                    RF_in_hemifield_index_CH(isnan(RF_in_hemifield_index_CH))=2;
                    RF_out_hemifield_index_CH(isnan(RF_out_hemifield_index_CH))=2;
                    if takemin
                        [~,RF_position_index_IN]=min(Average_FR_per_position_IN);
                        [~,RF_position_index_CH]=min(Average_FR_per_position_CH);
                    else
                        [~,RF_position_index_IN]=max(Average_FR_per_position_IN);
                        [~,RF_position_index_CH]=max(Average_FR_per_position_CH);
                    end
                    
                    %% full hemifield
                    idx_IN_prefHI_in  =idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_in_hemifield_index_IN)  & RF_TEST;
                    idx_IN_prefHO_in  =idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_out_hemifield_index_IN) & RF_TEST;
                    idx_CH_prefHI_in  =idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_in_hemifield_index_IN)  & RF_TEST;
                    idx_CH_prefHO_in  =idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_out_hemifield_index_IN) & RF_TEST;
                    
                    idx_IN_prefHI_ch  =idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_in_hemifield_index_CH);
                    idx_IN_prefHO_ch  =idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_out_hemifield_index_CH);
                    idx_CH_prefHI_ch  =idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_in_hemifield_index_CH);
                    idx_CH_prefHO_ch  =idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.LR(:,RF_out_hemifield_index_CH);
                    
                    %% baselines per hemifield???
                    idxb_IN_prefHI_in  =idx.in & idx.hands(:,hn) & ismember(epochs,b) & idx.LR(:,RF_in_hemifield_index_IN);
                    idxb_IN_prefHO_in  =idx.in & idx.hands(:,hn) & ismember(epochs,b) & idx.LR(:,RF_out_hemifield_index_IN);
                    idxb_CH_prefHI_in  =idx.ch & idx.hands(:,hn) & ismember(epochs,b) & idx.LR(:,RF_in_hemifield_index_IN);
                    idxb_CH_prefHO_in  =idx.ch & idx.hands(:,hn) & ismember(epochs,b) & idx.LR(:,RF_out_hemifield_index_IN);
                    
                    idxb_IN_prefHI_ch  =idx.in & idx.hands(:,hn) & ismember(epochs,b) & idx.LR(:,RF_in_hemifield_index_CH);
                    idxb_IN_prefHO_ch  =idx.in & idx.hands(:,hn) & ismember(epochs,b) & idx.LR(:,RF_out_hemifield_index_CH);
                    idxb_CH_prefHI_ch  =idx.ch & idx.hands(:,hn) & ismember(epochs,b) & idx.LR(:,RF_in_hemifield_index_CH);
                    idxb_CH_prefHO_ch  =idx.ch & idx.hands(:,hn) & ismember(epochs,b) & idx.LR(:,RF_out_hemifield_index_CH);
                    
                    if sum(idx_IN_prefHI_in)>=keys.cal.min_trials_per_condition && sum(idx_CH_prefHI_in)>=keys.cal.min_trials_per_condition
                        h=do_stats(FR(idx_IN_prefHI_in),FR(idx_CH_prefHI_in),keys,0);
                        DF=nanmean(FR(idx_CH_prefHI_in))-nanmean(FR(idx_IN_prefHI_in));
                        labelindex=h*sign(DF)+2; labelindex(isnan(labelindex))=2;
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:} '_prefH'])       =labelIN{labelindex};
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:} '_prefH'])       =labelIN{labelindex};
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:} '_prefH_FR'])    =nanmean(FR(idx_IN_prefHI_in));
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:} '_prefH_FR'])    =nanmean(FR(idx_CH_prefHI_in));
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:} '_prefH_DF'])    =nanmean(FR(idx_IN_prefHI_in))-nanmean(FR(idx_IN_prefHO_in));
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:} '_prefH_DF'])    =nanmean(FR(idx_CH_prefHI_in))-nanmean(FR(idx_CH_prefHO_in));
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:} '_prefH_IX'])    =(nanmean(FR(idx_IN_prefHI_in))-nanmean(FR(idx_IN_prefHO_in)))/nanmean(FR(idx_IN_prefHI_in | idx_IN_prefHO_in));
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:} '_prefH_IX'])    =(nanmean(FR(idx_CH_prefHI_in))-nanmean(FR(idx_CH_prefHO_in)))/nanmean(FR(idx_CH_prefHI_in | idx_CH_prefHO_in));
                    end
                    if sum(idx_CH_prefHI_in)>=keys.cal.min_trials_per_condition && sum(idx_CH_prefHO_in)>=keys.cal.min_trials_per_condition
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:}  '_prefHI_FR'])   =nanmean(FR(idx_CH_prefHI_in));
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:}  '_prefHO_FR'])   =nanmean(FR(idx_CH_prefHO_in));
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:}  '_prefHI_DF'])    =nanmean(FR(idx_CH_prefHI_in))-nanmean(FR(idxb_CH_prefHI_in));
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:}  '_prefHO_DF'])    =nanmean(FR(idx_CH_prefHO_in))-nanmean(FR(idxb_CH_prefHO_in));
                    end
                    if sum(idx_IN_prefHI_in)>=keys.cal.min_trials_per_condition && sum(idx_IN_prefHO_in)>=keys.cal.min_trials_per_condition
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:}  '_prefHI_FR'])   =nanmean(FR(idx_IN_prefHI_in));
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:}  '_prefHO_FR'])   =nanmean(FR(idx_IN_prefHO_in));
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:}  '_prefHI_DF'])    =nanmean(FR(idx_IN_prefHI_in))-nanmean(FR(idxb_IN_prefHI_in));
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:}  '_prefHO_DF'])    =nanmean(FR(idx_IN_prefHO_in))-nanmean(FR(idxb_IN_prefHO_in));
                    end
                    
                    if sum(idx_IN_prefHI_ch)>=keys.cal.min_trials_per_condition && sum(idx_CH_prefHI_ch)>=keys.cal.min_trials_per_condition
                        h=do_stats(FR(idx_IN_prefHI_ch),FR(idx_CH_prefHI_ch),keys,0);
                        DF=nanmean(FR(idx_CH_prefHI_ch))-nanmean(FR(idx_IN_prefHI_ch));
                        labelindex=h*sign(DF)+2; labelindex(isnan(labelindex))=2;
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:} '_prefHch'])       =labelIN{labelindex};
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:} '_prefHch'])       =labelIN{labelindex};
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:} '_prefHch_FR'])    =nanmean(FR(idx_IN_prefHI_ch));
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:} '_prefHch_FR'])    =nanmean(FR(idx_CH_prefHI_ch));
                    end
                    if sum(idx_CH_prefHI_ch)>=keys.cal.min_trials_per_condition && sum(idx_CH_prefHO_ch)>=keys.cal.min_trials_per_condition
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:}  '_prefHIch_FR'])   =nanmean(FR(idx_CH_prefHI_ch));
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:}  '_prefHOch_FR'])   =nanmean(FR(idx_CH_prefHO_ch));
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:}  '_prefHIch_DF'])    =nanmean(FR(idx_CH_prefHI_ch))-nanmean(FR(idxb_CH_prefHI_ch));
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:}  '_prefHOch_DF'])    =nanmean(FR(idx_CH_prefHO_ch))-nanmean(FR(idxb_CH_prefHO_ch));
                    end
                    if sum(idx_IN_prefHI_ch)>=keys.cal.min_trials_per_condition && sum(idx_IN_prefHO_ch)>=keys.cal.min_trials_per_condition
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:}  '_prefHIch_FR'])   =nanmean(FR(idx_IN_prefHI_ch));
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:}  '_prefHOch_FR'])   =nanmean(FR(idx_IN_prefHO_ch));
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:}  '_prefHIch_DF'])    =nanmean(FR(idx_IN_prefHI_ch))-nanmean(FR(idxb_IN_prefHI_ch));
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:}  '_prefHOch_DF'])    =nanmean(FR(idx_IN_prefHO_ch))-nanmean(FR(idxb_IN_prefHO_ch));
                    end
                    % preferred location
                    idx_IN_prefPI_in  =idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.pos==RF_position_index_IN;
                    idx_IN_prefPO_in  =idx.in & idx.hands(:,hn) & ismember(epochs,s) & idx.opp==RF_position_index_IN;
                    idx_CH_prefPI_in  =idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.pos==RF_position_index_IN;
                    idx_CH_prefPO_in  =idx.ch & idx.hands(:,hn) & ismember(epochs,s) & idx.opp==RF_position_index_IN;
                    
                    idxb_IN_prefPI_in  =idx.in & idx.hands(:,hn) & ismember(epochs,b) & idx.pos==RF_position_index_IN;
                    idxb_IN_prefPO_in  =idx.in & idx.hands(:,hn) & ismember(epochs,b) & idx.opp==RF_position_index_IN;
                    idxb_CH_prefPI_in  =idx.ch & idx.hands(:,hn) & ismember(epochs,b) & idx.pos==RF_position_index_IN;
                    idxb_CH_prefPO_in  =idx.ch & idx.hands(:,hn) & ismember(epochs,b) & idx.opp==RF_position_index_IN;
                    
                    %% here, the trial criterion was awkward, because > instead of >=
                    if sum(idx_IN_prefPI_in)>=keys.cal.min_trials_per_condition && sum(idx_CH_prefPI_in)>=keys.cal.min_trials_per_condition
                        h=do_stats(FR(idx_IN_prefPI_in),FR(idx_CH_prefPI_in),keys,0);
                        DF=nanmean(FR(idx_CH_prefPI_in))-nanmean(FR(idx_IN_prefPI_in));
                        labelindex=h*sign(DF)+2; labelindex(isnan(labelindex))=2;
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:}  '_prefP'])      =labelIN{labelindex};
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:}  '_prefP'])      =labelIN{labelindex};
                    end
                    
                    if sum(idx_IN_prefPI_in)>=keys.cal.min_trials_per_condition && sum(idx_CH_prefPI_in)>=keys.cal.min_trials_per_condition &&...
                            sum(idx_IN_prefPO_in)>=keys.cal.min_trials_per_condition && sum(idx_CH_prefPO_in)>=keys.cal.min_trials_per_condition
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:}  '_prefP_DF'])   =nanmean(FR(idx_IN_prefPI_in))-nanmean(FR(idx_IN_prefPO_in));
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:}  '_prefP_DF'])   =nanmean(FR(idx_CH_prefPI_in))-nanmean(FR(idx_CH_prefPO_in));
                    end
                    if sum(idx_IN_prefPI_in)>=keys.cal.min_trials_per_condition
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:}  '_prefPI_FR'])   =nanmean(FR(idx_IN_prefPI_in));
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:}  '_prefPI_DF'])   =nanmean(FR(idx_IN_prefPI_in))-nanmean(FR(idxb_IN_prefPI_in));
                    end
                    if sum(idx_IN_prefPO_in)>=keys.cal.min_trials_per_condition
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:}  '_prefPO_FR'])   =nanmean(FR(idx_IN_prefPO_in));
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:}  '_prefPO_DF'])   =nanmean(FR(idx_IN_prefPO_in))-nanmean(FR(idxb_IN_prefPO_in));
                    end
                    if sum(idx_CH_prefPI_in)>=keys.cal.min_trials_per_condition
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:}  '_prefPI_FR'])   =nanmean(FR(idx_CH_prefPI_in));
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:}  '_prefPI_DF'])   =nanmean(FR(idx_CH_prefPI_in))-nanmean(FR(idxb_CH_prefPI_in));
                    end
                    if sum(idx_CH_prefPO_in)>=keys.cal.min_trials_per_condition
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:}  '_prefPO_FR'])  =nanmean(FR(idx_CH_prefPO_in));
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:}  '_prefPO_DF'])   =nanmean(FR(idx_CH_prefPO_in))-nanmean(FR(idxb_CH_prefPO_in));
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
                    %% here, the trial criterion was awkward, because > instead of >=
                    if sum(idx_IN_prefPI_ch)>=keys.cal.min_trials_per_condition && sum(idx_CH_prefPI_ch)>=keys.cal.min_trials_per_condition
                        h=do_stats(FR(idx_IN_prefPI_ch),FR(idx_CH_prefPI_ch),keys,0);
                        DF=nanmean(FR(idx_CH_prefPI_ch))-nanmean(FR(idx_IN_prefPI_ch));
                        labelindex=h*sign(DF)+2; labelindex(isnan(labelindex))=2;
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:}  '_prefPch'])      =labelIN{labelindex};
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:}  '_prefPch'])      =labelIN{labelindex};
                    end
                    
                    if sum(idx_IN_prefPI_ch)>=keys.cal.min_trials_per_condition && sum(idx_CH_prefPI_ch)>=keys.cal.min_trials_per_condition &&...
                            sum(idx_IN_prefPO_ch)>=keys.cal.min_trials_per_condition && sum(idx_CH_prefPO_ch)>=keys.cal.min_trials_per_condition
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:}  '_prefPch_DF'])   =nanmean(FR(idx_IN_prefPI_ch))-nanmean(FR(idx_IN_prefPO_ch));
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:}  '_prefPch_DF'])   =nanmean(FR(idx_CH_prefPI_ch))-nanmean(FR(idx_CH_prefPO_ch));
                    end
                    if sum(idx_IN_prefPI_ch)>=keys.cal.min_trials_per_condition
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:}  '_prefPIch_FR'])   =nanmean(FR(idx_IN_prefPI_ch));
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:}  '_prefPIch_DF'])   =nanmean(FR(idx_IN_prefPI_ch))-nanmean(FR(idxb_IN_prefPI_ch));
                    end
                    if sum(idx_IN_prefPO_ch)>=keys.cal.min_trials_per_condition
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:}  '_prefPOch_FR'])  =nanmean(FR(idx_IN_prefPO_ch));
                        anova_struct.([INCHnamepart{1} '_' LHRHnamepart{hn} '_' s{:}  '_prefPOch_DF'])   =nanmean(FR(idx_IN_prefPO_ch))-nanmean(FR(idxb_IN_prefPO_ch));
                    end
                    if sum(idx_CH_prefPI_ch)>=keys.cal.min_trials_per_condition
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:}  '_prefPIch_FR'])   =nanmean(FR(idx_CH_prefPI_ch));
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:}  '_prefPIch_DF'])   =nanmean(FR(idx_CH_prefPI_ch))-nanmean(FR(idxb_CH_prefPI_ch));
                    end
                    if sum(idx_CH_prefPO_ch)>=keys.cal.min_trials_per_condition
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:}  '_prefPOch_FR'])  =nanmean(FR(idx_CH_prefPO_ch));
                        anova_struct.([INCHnamepart{2} '_' LHRHnamepart{hn} '_' s{:}  '_prefPOch_DF'])   =nanmean(FR(idx_CH_prefPO_ch))-nanmean(FR(idxb_CH_prefPO_ch));
                    end
                end
            else
                anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_position'])=labelsp{2};
            end
        end
    end
    
    %% subset ANOVA (!) - inactivation for example...
    % temporarly recompute this stuff    :((
    tr=idx.(INCHnamepart{ch}) & ismember(epochs,epoch_multicomp) & (idx.PT==0 | idx.PT==1); %& ismember(epochs,epochs_for_multicomparison);
    tr_main=idx.(INCHnamepart{ch}) & ismember(epochs,main_multicomp) & (idx.PT==0 | idx.PT==1); %& ismember(epochs,epochs_for_multicomparison);
    for fn=idx_fieldnames
        idx.(['tr_' fn{:}])=      idx.(fn{:}) & repmat(tr,1,size(idx.(fn{:}),2)) ;
    end
    
    multicomp_epochs=keys.epoch_SxH_multicomp;
    hand_space_multicomp=multicomp_epochs(ismember(multicomp_epochs,epochs'));
    
    perturbations_to_compare=unique(idx.tr_PT);
    if numel(perturbations_to_compare)>1
        conditions_to_compare={'LH_LS','LH_RS','RH_LS','RH_RS'};
        labelpt={'SU','-','EN'};
        for c=1:numel(conditions_to_compare)
            idx_con=idx.(['tr_' conditions_to_compare{c}]);
            
            if any(idx.tr_PT(tr_main & idx_con)==0) && any(idx.tr_PT(tr_main & idx_con)==1)
                [anova_out.p,anova_out.table,anova_out.stats,anova_out.terms] = anovan(FR(tr_main & idx_con),[epoch_idx(tr_main & idx_con) idx.tr_PT(tr_main & idx_con)],'model','full','varnames',{'epoch', 'perturbation'},'display',keys.plot.anova_tables);
                anova_struct.([INCHnamepart{ch} '_' conditions_to_compare{c} '_epoch_main'])    =anova_out.p(1)<0.05; %main effect on epoch!
                anova_struct.([INCHnamepart{ch} '_' conditions_to_compare{c} '_PT_main'])       =anova_out.p(2)<0.05; %main effect on epoch!
                anova_struct.([INCHnamepart{ch} '_' conditions_to_compare{c} '_ExP'])           =anova_out.p(3)<0.05; %main effect on epoch!
            end
            
            
            for s=hand_space_multicomp(:)'
                idxS= ismember(epochs,s);
                idx1=idx.tr_PT==0   & idx_con & idxS;
                idx2=idx.tr_PT==1   & idx_con & idxS;
                
                if sum(idx1)>0 && sum(idx2)>0
                    h=do_stats(FR(idx1),FR(idx2),keys,0);
                    DF=(nanmean(FR(idx2))-nanmean(FR(idx1)));
                else
                    h=false;
                    DF=single(NaN);
                end
                labelindexpt=h*sign(DF)+2; labelindexpt(isnan(labelindexpt))=2;
                anova_struct.([INCHnamepart{ch} '_' conditions_to_compare{c} '_' s{:} '_PT']) = labelpt{labelindexpt}; %
                anova_struct.([INCHnamepart{ch} '_' conditions_to_compare{c} '_' s{:} '_PT_DF']) = DF; %
                anova_struct.([INCHnamepart{ch} '_' conditions_to_compare{c} '_' s{:} '_PT_IX']) = DF/(nanmean(FR(idx2))+nanmean(FR(idx1)));
                anova_struct.([INCHnamepart{ch} '_' conditions_to_compare{c} '_' s{:} '_CT_FR']) = nanmean(FR(idx1)); %
                anova_struct.([INCHnamepart{ch} '_' conditions_to_compare{c} '_' s{:} '_PT_FR']) = nanmean(FR(idx2)); %
                
            end
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
                    h=do_stats(FR(idx1)-FR(idx1b),FR(idx2)-FR(idx2b),keys,0);
                    DFPT = nanmean(FR(idx2)-FR(idx2b));
                    DFCT = nanmean(FR(idx1)-FR(idx1b));
                    DF=DFPT-DFCT;
                else
                    h=false;
                    DFPT=single(NaN);
                    DFCT=single(NaN);
                    DF=single(NaN);
                end
                
                labelindexpt=h*sign(DF)+2; labelindexpt(isnan(labelindexpt))=2;
                anova_struct.([INCHnamepart{ch} '_' conditions_to_compare{c} '_' s{:} '_PTbl']) = labelpt{labelindexpt}; %
                anova_struct.([INCHnamepart{ch} '_' conditions_to_compare{c} '_' s{:} '_PTbl_DF']) = DF; %
                anova_struct.([INCHnamepart{ch} '_' conditions_to_compare{c} '_' s{:} '_PT_epoch_DF']) = DFPT; %
                anova_struct.([INCHnamepart{ch} '_' conditions_to_compare{c} '_' s{:} '_CT_epoch_DF']) = DFCT; %
            end
            
            
        end
        %         labelindexcr=(anova_out.p(k)<0.05) *sign(nanmean(FR(tr & idx.tr_CR))-nanmean(FR(tr & idx.tr_UC)))+2; labelindexcr(isnan(labelindexcr))=2;% 3 crossed>uncrossed ,1 crossed<uncrossed
        %         anova_struct.([INCHnamepart{ch} '_SxH'])=labelcr{labelindexcr};
        %         anova_struct.([INCHnamepart{ch} '_ExSxH'])=anova_out.p(7)<0.05;
        %
        %         for s=multicomp_epochs(:)'
        %             idxS= tr & ismember(epochs,s);
        %             if any(~isnan(FR(idxS))) && any(idx.tr_LS(idxS)) && any(idx.tr_RS(idxS)) && any(idx.tr_LH(idxS)) && any(idx.tr_RH(idxS))
        %                 [anova_outs.p,anova_outs.table,anova_outs.stats,anova_outs.terms] = anovan(FR(idxS),double([idx.tr_RS(idxS),idx.tr_RH(idxS)]),'model','full','varnames',{'space';'hands'},'display',keys.plot.anova_tables);
        %
        %                 %% replacing hand and space tuning with anova per epoch !!
        %                 h=anova_outs.p(1)<0.05;
        %                 ES=nanmean(FR(idx.tr_RS & idxS))-nanmean(FR(idx.tr_LS & idxS));
        %                 labelindexsp=h*sign(ES)+2; labelindexsp(isnan(labelindexsp))=2;
        %                 anova_struct.([INCHnamepart{ch} '_' s{:} '_spaceLR'])= labelsp{labelindexsp}; %
        %                 anova_struct.([INCHnamepart{ch} '_' s{:} '_spaceLR_ES'])= ES; %
        %                 anova_struct.([INCHnamepart{ch} '_' s{:} '_spaceLR_IX'])= ES/(nanmean(FR(idx.tr_RS & idxS))+nanmean(FR(idx.tr_LS & idxS))); %
        %
        %                 h=anova_outs.p(2)<0.05;
        %                 ES=nanmean(FR(idx.tr_RH & idxS))-nanmean(FR(idx.tr_LH & idxS));
        %                 labelindexha=h*sign(ES)+2; labelindexha(isnan(labelindexha))=2;
        %                 anova_struct.([INCHnamepart{ch} '_' s{:} '_hands'])= labelha{labelindexha};
        %                 anova_struct.([INCHnamepart{ch} '_' s{:} '_hands_ES'])= ES;
        %                 anova_struct.([INCHnamepart{ch} '_' s{:} '_hands_IX'])= ES/(nanmean(FR(idx.tr_RH & idxS))+nanmean(FR(idx.tr_LH & idxS)));
        %
        %                 h=anova_outs.p(3)<0.05;
        %                 ES=nanmean(FR(idx.tr_CR & idxS))-nanmean(FR(idx.tr_UC & idxS));
        %                 labelindexcr=h*sign(ES)+2; labelindexcr(isnan(labelindexcr))=2;
        %                 anova_struct.([INCHnamepart{ch} '_' s{:} '_SxH'])=labelcr{labelindexcr};
        %                 anova_struct.([INCHnamepart{ch} '_' s{:} '_SxH_ES'])=ES;
        %                 anova_struct.([INCHnamepart{ch} '_' s{:} '_SxH_IX'])=ES/(nanmean(FR(idx.tr_CR & idxS))+nanmean(FR(idx.tr_UC & idxS)));
        %             end
        %         end
        
    end
end
end

function anova_struct=effector_comparison_anova(keys,epochs,FR,idx,Positions,Fixations)
anova_struct=struct();
INCHnamepart=keys.labels.choices([sum(idx.in)>0 sum(idx.ch)>0]);
LHRHnamepart=keys.labels.handsLR;
labelen={'su','-','en','bi'};
labelsp={'LS','-','RS'};
labelha={'LH','-','RH'};
labelcr={'UC','-','CR'};
labelef=keys.labels.eff;
labelIN={'IN','-','CH'};
true_labels={'false','true'};

conditions_to_compare={'LH_LS','LH_RS','RH_LS','RH_RS'};

idx_ep=ismember(keys.epoch_multicomp(:,1),keys.EPOCHS(:,1)) &  ismember(keys.epoch_multicomp(:,2),keys.EPOCHS(:,1));
epoch_multicomp=keys.epoch_multicomp(idx_ep,:);
idx_ep=ismember(keys.main_multicomp,keys.EPOCHS(:,1));
main_multicomp=keys.main_multicomp(idx_ep,:);

[~, ~, epoch_idx]=unique(epochs);
%% Independently for Instructed and Choice !!
idx_fieldnames=fieldnames(idx)';
for ch=1:numel(INCHnamepart)
    
    %% this part should go outside the loop!
    tr=idx.(INCHnamepart{ch}) & ismember(epochs,epoch_multicomp); %& ismember(epochs,epochs_for_multicomparison);
    for fn=idx_fieldnames
        idx.(['tr_' fn{:}])=      idx.(fn{:}) & repmat(tr,1,size(idx.(fn{:}),2));
    end
    
    %% 3-way ANOVA per epoch...
    
    multicomp_epochs=keys.epoch_SxH_multicomp; %% needs to be replaced by something meaningful
    for s=multicomp_epochs(:)'
        idxS= tr & ismember(epochs,s);
        if any(~isnan(FR(idxS))) && any(idx.tr_LS(idxS)) && any(idx.tr_RS(idxS)) && any(idx.tr_LH(idxS)) && any(idx.tr_RH(idxS))
            [anova_outs.p,anova_outs.table,anova_outs.stats,anova_outs.terms] = anovan(FR(idxS),double([idx.tr_RS(idxS),idx.tr_RH(idxS),idx.tr_E2(idxS)]),'model','full','varnames',{'space';'hands';'effector'},'display',keys.plot.anova_tables);
            
            %% replacing hand and space tuning with anova per epoch !!
            h=anova_outs.p(1)<0.05;
            ES=nanmean(FR(idx.tr_RS & idxS))-nanmean(FR(idx.tr_LS & idxS));
            labelindexsp=h*sign(ES)+2; labelindexsp(isnan(labelindexsp))=2;
            anova_struct.([INCHnamepart{ch} '_' s{:} '_spaceLR'])= labelsp{labelindexsp}; %
            anova_struct.([INCHnamepart{ch} '_' s{:} '_spaceLR_ES'])= ES; %
            anova_struct.([INCHnamepart{ch} '_' s{:} '_spaceLR_IX'])= ES/(nanmean(FR(idx.tr_RS & idxS))+nanmean(FR(idx.tr_LS & idxS))); %
            
            h=anova_outs.p(2)<0.05;
            ES=nanmean(FR(idx.tr_RH & idxS))-nanmean(FR(idx.tr_LH & idxS));
            labelindexha=h*sign(ES)+2; labelindexha(isnan(labelindexha))=2;
            anova_struct.([INCHnamepart{ch} '_' s{:} '_hands'])= labelha{labelindexha};
            anova_struct.([INCHnamepart{ch} '_' s{:} '_hands_ES'])= ES;
            anova_struct.([INCHnamepart{ch} '_' s{:} '_hands_IX'])= ES/(nanmean(FR(idx.tr_RH & idxS))+nanmean(FR(idx.tr_LH & idxS)));
            
            h=anova_outs.p(3)<0.05;
            ES=nanmean(FR(idx.tr_E2 & idxS))-nanmean(FR(idx.tr_E1 & idxS));
            labelindexef=h*sign(ES)+2; labelindexef(isnan(labelindexef))=2;
            anova_struct.([INCHnamepart{ch} '_' s{:} '_effector'])=labelef{labelindexef};
            anova_struct.([INCHnamepart{ch} '_' s{:} '_effector_ES'])=ES;
            anova_struct.([INCHnamepart{ch} '_' s{:} '_effector_IX'])=ES/(nanmean(FR(idx.tr_E2 & idxS))+nanmean(FR(idx.tr_E1 & idxS)));
            
            for c=1:numel(conditions_to_compare)
                idx1=idx.tr_E1   & idx.(['tr_' conditions_to_compare{c}]) & idxS;
                idx2=idx.tr_E2   & idx.(['tr_' conditions_to_compare{c}]) & idxS;
                
                if sum(idx1)>0 && sum(idx2)>0
                    h=do_stats(FR(idx1),FR(idx2),keys,0);
                    ES=(nanmean(FR(idx2))-nanmean(FR(idx1)))/(nanmean(FR(idx2))+nanmean(FR(idx1)));
                else
                    h=false;
                    ES=single(NaN);
                end
                labelindexef=h*sign(ES)+2; labelindexef(isnan(labelindexef))=2;
                anova_struct.([INCHnamepart{ch} '_' conditions_to_compare{c} '_' s{:} '_effector']) = labelef{labelindexef};
                anova_struct.([INCHnamepart{ch} '_' conditions_to_compare{c} '_' s{:} '_effector_ES']) = ES; %
            end
            
            
            %freaking interactions
            %
            %                 h=anova_outs.p(4)<0.05;
            %                 ES=nanmean(FR(idx.tr_CR & idxS))-nanmean(FR(idx.tr_UC & idxS));
            %                 labelindexcr=h*sign(ES)+2; labelindexcr(isnan(labelindexcr))=2;
            %                 anova_struct.([INCHnamepart{ch} '_' s{:} '_SxH'])=labelcr{labelindexcr};
            %                 anova_struct.([INCHnamepart{ch} '_' s{:} '_SxH_ES'])=ES;
            %                 anova_struct.([INCHnamepart{ch} '_' s{:} '_SxH_IX'])=ES/(nanmean(FR(idx.tr_CR & idxS))+nanmean(FR(idx.tr_UC & idxS)));
        end
    end
    
end
end

function trial_criterion=get_minimum_trials(keys,o,tr_index)
INCHnamepart={'in','ch'};
LHRHnamepart={'AH','NH','LH','RH'};
PTnamepart={'','PT_'};

o.position_combinations=o.position_combinations(tr_index,:);
o.hemifield_combinations=o.hemifield_combinations(tr_index,:);

hands           =[o.trial.reach_hand];
choices         =[o.trial.choice];
perturbations    =[o.trial.perturbation];

possible_combinations=combvec([-1 keys.cal.reach_hand],keys.cal.choice,1:numel(keys.cal.perturbation_groups))';
unique_position_indexes=unique(o.position_combinations(:,3));
unique_hemifield_indexes=unique(o.hemifield_combinations(:,3));
nonnanidx=~any(isnan(o.position_combinations),2)';

clear trial_criterion %n_trials_per_condition
for c=1:size(possible_combinations,1)
    
    if possible_combinations(c,1)==-1 %% any hand, relevant later
        idx=nonnanidx & choices==possible_combinations(c,2) & ismember(perturbations,keys.cal.perturbation_groups{possible_combinations(c,3)});
    else
        idx=nonnanidx & hands==possible_combinations(c,1) & choices==possible_combinations(c,2) & ismember(perturbations,keys.cal.perturbation_groups{possible_combinations(c,3)});
    end
    namepart=[INCHnamepart{possible_combinations(c,2)+1} '_' LHRHnamepart{possible_combinations(c,1)+2} '_' PTnamepart{possible_combinations(c,3)}];
    
    
    %% trials_total    
    trial_criterion.([namepart 'trials_total'])=sum(idx);
    
    %% trials_per_condition (should be renamed per_position) **    
    [~, ~, unique_combinations]= unique(o.position_combinations(idx,:),'rows');
    if all(ismember(unique_position_indexes,o.position_combinations(idx,3)))
        least_trials_in_position=hist(unique_combinations,1:max(unique_combinations));
    else % in case the postiion wasnt succesfully selected at all (especially imoprtant for choice!)
        least_trials_in_position=0;
    end
    trial_criterion.([namepart 'trials_per_condition'])=min(least_trials_in_position);
    %n_trials_per_condition(c)=min(least_trials_in_position);
    
    %% trials_per_hemfield    
    [~, ~, unique_hemifield]= unique(o.hemifield_combinations(idx,:),'rows');
    if all(ismember(unique_hemifield_indexes,o.hemifield_combinations(idx,3)))
        least_trials_in_hemifield=hist(unique_hemifield,1:max(unique_hemifield));
    else % in case the postiion wasnt succesfully selected at all (especially imoprtant for choice!)
        least_trials_in_hemifield=0;
    end
    trial_criterion.([namepart 'trials_per_hemifield'])=min(least_trials_in_hemifield);
    
    %% trials_per_congruent_hand_hemifeld    
    unique_hemifield = o.hemifield_combinations(idx,3);
    if possible_combinations(c,1)==1 %left hand
        unique_hemifield(unique_hemifield==2)=[];
    elseif possible_combinations(c,1)==2 %right hand
        unique_hemifield(unique_hemifield==1)=[];
    end
    congruent_trials=numel(unique_hemifield);
    trial_criterion.([namepart 'trials_per_congruent_hand_hemifield'])=congruent_trials;
end
end

function h=do_stats(A,B,keys,paired)
%h = ttest2(A,B);
if paired
    if any(~isnan(A)&~isnan(B))
        [~, h] = signrank(A,B);
    else
        h=0;
    end
else
    if any(~isnan(A)) && any (~isnan(B))
        [~, h] = ranksum(A,B);
    else
        h=0;
    end
end
end

%% these are some subfunctions almost ready to be implemented

function oo=ph_FR_subtract_baseline(o,keys)
oo=o;
for t=1:size(o,1)
    epochs={o(t,:).state};
    for e=1:size(o,2)
        s=ismember(keys.EPOCHS(:,1),o(t,e).state);
        b=ismember(epochs,keys.EPOCHS(s,5));
        if any(b)
            oo(t,e).FR=o(t,e).FR - o(t,b).FR;
        end
    end
end
end

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
