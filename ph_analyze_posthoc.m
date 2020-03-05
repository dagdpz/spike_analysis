function tuning_per_unit_table=ph_analyze_posthoc(o_a,keys)
global MA_STATES
%population_analysis_table{1,:}={};
tuning_per_unit_table=keys.tuning_per_unit_table;


for unit=1:numel(o_a)
    anova_struct_current_unit=struct();
    conditions=[o_a(unit).condition];
    for con=1:numel(conditions)
        effector=conditions(con).effector;
        type=conditions(con).type;
        
        keys.epoch_main_multicomp               =keys.epoch_comparisons{type};
        keys.epoch_for_multicomparison          =keys.epochs_for_multicomparison{type};
        keys.epoch_spaceLR_multicomp            =keys.epochs_spaceLR_multicomp{type};
        keys.epoch_choice_multicomp             =keys.epochs_choice_multicomp{type};
        keys.epoch_hands_multicomp              =keys.epochs_hands_multicomp{type};
        keys.epoch_SxH_multicomp                =keys.epochs_SxH_multicomp{type};
        
        get_expected_MA_states(type,effector,keys.effectors_on_same_figure);
        keys.ALL_EPOCHS=keys.EPOCHS_PER_TYPE{type};
        keys.EPOCHS=keys.ALL_EPOCHS(ismember([keys.ALL_EPOCHS{:,2}],MA_STATES.all_states),:);
        
        %condition_fieldname_part=['eff_' num2str(effector) '_typ_' num2str(type)];
        [~, condition_fieldname_part]=get_type_effector_name(type,effector);
        for cas=1:size(o_a(unit).condition(con).case,2)
            n_left=numel([conditions(con).case(cas).left.trial]);
            n_right=numel([conditions(con).case(cas).right.trial]);
            if n_left>0 && n_right>0
                % Position wise
                pos                         =[o_a(unit).condition(con).case(cas).position];
                pos_emp                     =arrayfun(@(x) isempty(x.position), pos);
                pos(pos_emp)                =[];
                pos_tmp                     =arrayfun(@(x) repmat(x.position(1)+1i*x.position(2),1,numel(x.per_state)), pos,'UniformOutput',false);
                pos_tmp                     =num2cell([pos_tmp{:}]);
                Pos_per_state               =[pos.per_state];
                [Pos_per_state.position]    =pos_tmp{:};
                Pos_trial                   =[Pos_per_state.trial];
                sta_tmp                     =arrayfun(@(x) repmat({x.state},1,numel(x.trial)), Pos_per_state,'UniformOutput',false);
                pos_pos                     =arrayfun(@(x) repmat(x.position,1,numel(x.trial)), Pos_per_state,'UniformOutput',false);
                States                      =[sta_tmp{:}]';
                pos_pos                     =[pos_pos{:}]';
                [u_pos,~,pos_pid]           =unique(pos_pos);
                %[~, ~, pos_sta]             =unique(States);
                pos_hnd                     =vertcat(Pos_trial.hand);
                pos_cho                     =vertcat(Pos_trial.choice);
                pos_FR                      =vertcat(Pos_trial.FR);
                
                % Left Right
                Left=[conditions(con).case(cas).left.per_state.trial];
                Right=[conditions(con).case(cas).right.per_state.trial];
                Left_IN=Left([Left.choice]==0);
                Left_CH=Left([Left.choice]==1);
                Left_IN_LH=Left([Left.choice]==0 & [Left.hand]==1);
                Left_IN_RH=Left([Left.choice]==0 & [Left.hand]==2);
                Left_CH_LH=Left([Left.choice]==1 & [Left.hand]==1);
                Left_CH_RH=Left([Left.choice]==1 & [Left.hand]==2);
                Right_IN=Right([Right.choice]==0);
                Right_CH=Right([Right.choice]==1);
                Right_IN_LH=Right([Right.choice]==0 & [Right.hand]==1);
                Right_IN_RH=Right([Right.choice]==0 & [Right.hand]==2);
                Right_CH_LH=Right([Right.choice]==1 & [Right.hand]==1);
                Right_CH_RH=Right([Right.choice]==1 & [Right.hand]==2);
                
                % number of trials for each condition to decide about type
                % of anova
                n_IN=numel(Left_IN)>0 && numel(Right_IN)>0;
                n_CH=numel(Left_CH)>0 && numel(Right_CH)>0;
                n_IN_2hands=numel(Left_IN_LH)>0 && numel(Right_IN_LH)>0 && numel(Left_IN_RH)>0 && numel(Right_IN_RH)>0;
                n_CH_2hands=numel(Left_CH_LH)>0 && numel(Right_CH_LH)>0 && numel(Left_CH_RH)>0 && numel(Right_CH_RH)>0;
                if n_IN_2hands && n_CH_2hands % four-way anova epoch, space, hand, choice
                    keys.anova_varnames={'epoch' 'spaceLR' 'hands' 'choice'};
                    %[anova_struct]=n_way_anova(keys, Left_IN_LH,Left_CH_LH,Left_IN_RH,Left_CH_RH,Right_IN_LH,Right_CH_LH,Right_IN_RH,Right_CH_RH);
                    [anova_struct]=n_way_anova(keys, Left_IN_LH,Left_IN_RH,Right_IN_LH,Right_IN_RH,Left_CH_LH,Left_CH_RH,Right_CH_LH,Right_CH_RH);
                    anova_struct_current_unit.([condition_fieldname_part '_'  keys.case_summaries{cas}(1:3)])=anova_struct; clear anova_struct
                elseif n_IN_2hands && ~n_CH_2hands % three-way anova epoch, space, hand
                    keys.anova_varnames={'epoch' 'spaceLR' 'hands'};
                    %[anova_struct]=n_way_anova(keys, Left_IN_LH,Left_IN_RH,Right_IN_LH,Right_IN_RH);
                    [anova_struct]=n_way_anova(keys, Left_IN_LH,Left_IN_RH,Right_IN_LH,Right_IN_RH);
                    anova_struct_current_unit.([condition_fieldname_part '_'  keys.case_summaries{cas}(1:3)])=anova_struct; clear anova_struct
                elseif ~n_IN_2hands && n_IN && n_CH % three-way anova epoch, space, choice
                    keys.anova_varnames={'epoch' 'spaceLR' 'choice'};
                    %[anova_struct]=n_way_anova(keys, Left_IN,Left_CH,Right_IN,Right_CH);
                    [anova_struct]=n_way_anova(keys, Left_IN,Right_IN,Left_CH,Right_CH);
                    anova_struct_current_unit.([condition_fieldname_part '_'  keys.case_summaries{cas}(1:3)])=anova_struct; clear anova_struct
                elseif n_IN % two-way anova epoch, space
                    keys.anova_varnames={'epoch' 'spaceLR'};
                    [anova_struct]=n_way_anova(keys,Left,Right);
                    anova_struct_current_unit.([condition_fieldname_part '_'  keys.case_summaries{cas}(1:3)])=anova_struct; clear anova_struct
                    
                else
                    disp('Not fitting dataset')
                end
                
                pos_anova_tmp=ph_pos_anova(keys,pos_FR,States,pos_cho,pos_pid,u_pos,pos_hnd);
                for FN=fieldnames(pos_anova_tmp)'
                    anova_struct_current_unit.([condition_fieldname_part '_'  keys.case_summaries{cas}(1:3)]).(FN{:})=pos_anova_tmp.(FN{:});
                end
                anova_struct_current_unit.([condition_fieldname_part '_'  keys.case_summaries{cas}(1:3)]).trialsINLR= min(numel(Left_IN),numel(Right_IN))/size(keys.EPOCHS,1);
            else
                disp('No trial for left or right side')
            end
            
            
        end
    end
    
    
    rows_to_update=find(ismember(tuning_per_unit_table(:,find_column_index(tuning_per_unit_table,'unit_ID')),{o_a(unit).unit_ID}));
    if isempty(rows_to_update)
        rows_to_update=size(tuning_per_unit_table,1)+1;
    end
    
    clear unit_table;
    inital_fieldnames={'unit_ID','target','grid_x','grid_y','electrode_depth','stability_rating','SNR_rating','Single_rating'};
    unit_table(1,1:numel(inital_fieldnames))=inital_fieldnames;
    for fn=1:numel(inital_fieldnames)
        unit_table{rows_to_update,fn}=o_a(unit).(inital_fieldnames{fn});
    end
    title_counter=numel(inital_fieldnames);
    
    FN=fieldnames(anova_struct_current_unit);
    for fn=1:numel(FN)
        FNsub=fieldnames(anova_struct_current_unit.(FN{fn}));
        for fnsub=1:numel(FNsub)
            title_counter=title_counter+1;
            %table_struct.([FN{fn} '_' FNsub{fnsub}])=anova_struct_current_unit.(FN{fn}).(FNsub{fnsub});
            unit_table{1,title_counter}=[FNsub{fnsub} '_' FN{fn}];
            unit_table{rows_to_update,title_counter}=anova_struct_current_unit.(FN{fn}).(FNsub{fnsub});
        end
    end
    tuning_per_unit_table=update_mastertable_cell(tuning_per_unit_table,unit_table,rows_to_update);
    
end

end

function [anova_struct]=n_way_anova(keys,varargin)

idx_ep=ismember(keys.epoch_main_multicomp(:,1),keys.EPOCHS(:,1)) &  ismember(keys.epoch_main_multicomp(:,2),keys.EPOCHS(:,1));
states_for_multicomparison=keys.epoch_for_multicomparison(ismember(keys.epoch_for_multicomparison,keys.EPOCHS(:,1)));
epoch_main_multicomp=keys.epoch_main_multicomp(idx_ep,:);


choice_varname=ismember(keys.anova_varnames,'choice');
varnames_wo_choice=keys.anova_varnames(~choice_varname);
if any(choice_varname)
    numargin=numel(varargin)/2;
    in_ch=[0,numel(varargin)/2];
    parnamepart={'in','ch'};
else
    numargin=numel(varargin);
    in_ch=0;
    parnamepart={'in'};
end

%% Independently for Instructed and Choice !!
for ch=1:numel(in_ch)
    for n=1:numargin
        state_indexes=ismember({varargin{n+in_ch(ch)}.state},states_for_multicomparison);
        %tmp=[varargin{n}(state_indexes).trial];
        FR_cell{n}=[varargin{n+in_ch(ch)}(state_indexes).FR];
        States_cell{n}={varargin{n+in_ch(ch)}(state_indexes).state};
        N{n}=numel(FR_cell{n});
    end
    FR=[FR_cell{:}]';
    States=[States_cell{:}]';
    [~, ~, Par{1}]=unique(States);
    
    %% main effects and general interactions
    if numargin==8; %four-way
        Par{2}=[zeros(N{1}+N{2}+N{3}+N{4},1); ones(N{5}+N{6}+N{7}+N{8},1)];
        Par{3}=[zeros(N{1}+N{2},1); ones(N{3}+N{4},1); zeros(N{5}+N{6},1); ones(N{7}+N{8},1)];
        Par{4}=[zeros(N{1},1); ones(N{2},1); zeros(N{3},1); ones(N{4},1); zeros(N{5},1); ones(N{6},1); zeros(N{7},1); ones(N{8},1) ];
        [anova_out.p,anova_out.table,anova_out.stats,anova_out.terms] = anovan(FR,[Par{1},Par{2},Par{3},Par{4}],'model','full','varnames',varnames_wo_choice,'display',keys.display_anova_tables);
    elseif numargin==4; %three-way
        Par{2}=[zeros(N{1}+N{2},1); ones(N{3}+N{4},1)];
        Par{3}=[zeros(N{1},1); ones(N{2},1); zeros(N{3},1); ones(N{4},1)];
        [anova_out.p,anova_out.table,anova_out.stats,anova_out.terms] = anovan(FR,[Par{1},Par{2},Par{3}],'model','full','varnames',varnames_wo_choice,'display',keys.display_anova_tables);
    else
        Par{2}=[zeros(N{1},1); ones(N{2},1)];
        [anova_out.p,anova_out.table,anova_out.stats,anova_out.terms] = anovan(FR,[Par{1},Par{2}],'model','full','varnames',varnames_wo_choice,'display',keys.display_anova_tables);
    end
    anova_struct.([parnamepart{ch} '_epoch_main'])=anova_out.p(1)<0.05; %main effect on epoch!
    
    
    %% epoch tuning regardless of hand
    S=find(ismember(varnames_wo_choice,'spaceLR'));
    if S
        idx.sides={Par{S}==0,Par{S}==1};
    else
        idx.sides={true(size(States)),true(size(States))};
    end
    
    H=find(ismember(varnames_wo_choice,'hands'));
    if H
        handindexes=[1,2,3];
        idx.hands={true(size(States)),Par{H}==0,Par{H}==1}; %% check 0 is left, 1 is right
    else
        handindexes=[1];
        idx.hands={true(size(States))};
    end
    
    label={'su','-','en','bi'};
    %sides={'L','R'};
    LHRHnamepart={'NH','LH','RH'};
    multicomp_states=epoch_main_multicomp;
    for row=1:size(multicomp_states,1)
        for handindex=handindexes
            for sideindex=1:2
                %side=sides(sideindex);
                idx1=ismember(States,multicomp_states(row,1)) & idx.sides{sideindex} & idx.hands{handindex};
                idx2=ismember(States,multicomp_states(row,2)) & idx.sides{sideindex} & idx.hands{handindex};
                if sum(idx1)>0 && sum(idx2)>0
                    [~,p]=ttest(FR(idx1),FR(idx2)); %%paired ttest !!??
                    h=p<0.05;
                else
                    h=false;
                end
                %h=p*numel(multicomp_states)<0.05;%% bonferroni multicomparison correction
                ES(sideindex)=(nanmean(FR(idx2))-nanmean(FR(idx1)));
                ensu(sideindex)=h*sign(ES(sideindex))+2;
            end
            ensu(isnan(ensu))=2;
            if any(ensu==1) && any(ensu==3)
                labelindex=4;
            elseif any(ensu~=2)
                labelindex=unique(ensu(ensu~=2));
            else
                labelindex=2;
            end
            [~,ESidx]=max(abs(ES));
            ES=ES(ESidx);
            anova_struct.([parnamepart{ch} '_' LHRHnamepart{handindex}  '_' multicomp_states{row,2} '_epoch'])=label{labelindex};
            anova_struct.([parnamepart{ch} '_' LHRHnamepart{handindex}  '_' multicomp_states{row,2} '_epoch_ES'])=ES; %%% to think about !!!!
            anova_struct.([parnamepart{ch} '_' LHRHnamepart{handindex}  '_' multicomp_states{row,2} '_epoch_IX'])=ES/nanmean(FR(idx2))+nanmean(FR(idx1)); %%% to think about !!!!
        end
    end
        
    
    %% hand and space tuning
    for k=2:numel(varnames_wo_choice)
        epoch_interaction_effect=anova_out.p(k+numel(varnames_wo_choice)-1)<0.05;
        %Par_k=Pa.(['r' num2str(k)]);
        multicomp_states=keys.(['epoch_' varnames_wo_choice{k} '_multicomp']);
        multicomp_states=multicomp_states(ismember(multicomp_states,States'));
        switch varnames_wo_choice{k}
            case 'spaceLR'
                label={'LS','-','RS'};
            case 'hands'
                label={'LH','-','RH'};
        end
        
        idx1=Par{k}==0;
        idx2=Par{k}==1;
        labelindex=(anova_out.p(k)<0.05) *sign(nanmean(FR(idx2))-nanmean(FR(idx1)))+2; labelindex(isnan(labelindex))=2;
        anova_struct.([parnamepart{ch} '_' varnames_wo_choice{k} '_main'])=label{labelindex}; %main effect with direction!
        anova_struct.([parnamepart{ch} '_Ex' upper(varnames_wo_choice{k}(1))])=epoch_interaction_effect; %main effect with direction!
        
        for s=multicomp_states(:)'
            idxS=ismember(States,s);
            [~,p]=ttest2(FR(idx1 & idxS),FR(idx2 & idxS));
            h=p<0.05;
            %h=p*numel(multicomp_states)<0.05;%% bonferroni multicomparison correction
            ES=nanmean(FR(idx2 & idxS))-nanmean(FR(idx1 & idxS));
            labelindex=h*sign(ES)+2; labelindex(isnan(labelindex))=2;
            anova_struct.([parnamepart{ch} '_' s{:} '_' varnames_wo_choice{k}])=label{labelindex};
            anova_struct.([parnamepart{ch} '_' s{:} '_' varnames_wo_choice{k} '_ES'])=ES;
            anova_struct.([parnamepart{ch} '_' s{:} '_' varnames_wo_choice{k} '_IX'])=ES/(nanmean(FR(idx2 & idxS))+ nanmean(FR(idx1 & idxS)));
        end
    end
    
    %% space hand anovas
    if ismember('spaceLR',varnames_wo_choice) && ismember('hands',varnames_wo_choice)
        k=6;
        multicomp_states=keys.epoch_SxH_multicomp;
        multicomp_states=multicomp_states(ismember(multicomp_states,States'));
        idx_LS=Par{2}==0;
        idx_RS=Par{2}==1;
        idx_LH=Par{3}==0;
        idx_RH=Par{3}==1;
        crossed   = (idx_LS & idx_RH) | (idx_RS & idx_LH);
        uncrossed = (idx_LS & idx_LH) | (idx_RS & idx_RH);
        labelcr={'UC','-','CR'};
        labelha={'LH','-','RH'};
        labelsp={'LS','-','RS'};
        labelindexcr=(anova_out.p(k)<0.05) *sign(nanmean(FR(crossed))-nanmean(FR(uncrossed)))+2; labelindexcr(isnan(labelindexcr))=2;
        anova_struct.([parnamepart{ch} '_SxH'])=labelcr{labelindexcr}; % 1 crossed>uncrossed ,-1 crossed<uncrossed
        anova_struct.([parnamepart{ch} '_ExSxH'])=anova_out.p(7)<0.05; % 1 crossed>uncrossed ,-1 crossed<uncrossed
        
        for s=multicomp_states(:)'
            idxS=ismember(States,s);
            if any(~isnan(FR(idxS))) && any(Par{2}(idxS)==1) && any(Par{2}(idxS)==0) && any(Par{3}(idxS)==1) && any(Par{3}(idxS)==0)
                [anova_outs.p,anova_outs.table,anova_outs.stats,anova_outs.terms] = anovan(FR(idxS),[Par{2}(idxS),Par{3}(idxS)],'model','full','varnames',{'space';'hands'},'display',keys.display_anova_tables);
         
                %% replacing hand and space tuning with anova per epoch !!
                h=anova_outs.p(1)<0.05;
                ES=nanmean(FR(Par{2}==1 & idxS))-nanmean(FR(Par{2}==0 & idxS));
                labelindexsp=h*sign(ES)+2; labelindexsp(isnan(labelindexsp))=2;
                anova_struct.([parnamepart{ch} '_' s{:} '_spaceLR'])= labelsp{labelindexsp}; %    
                anova_struct.([parnamepart{ch} '_' s{:} '_spaceLR_ES'])= ES; %               
                anova_struct.([parnamepart{ch} '_' s{:} '_spaceLR_IX'])= ES/(nanmean(FR(Par{2}==1 & idxS))+nanmean(FR(Par{2}==0 & idxS))); %              
                
                h=anova_outs.p(2)<0.05;
                ES=nanmean(FR(Par{3}==1 & idxS))-nanmean(FR(Par{3}==0 & idxS));
                labelindexha=h*sign(ES)+2; labelindexha(isnan(labelindexha))=2;
                anova_struct.([parnamepart{ch} '_' s{:} '_hands'])= labelha{labelindexha}; 
                anova_struct.([parnamepart{ch} '_' s{:} '_hands_ES'])= ES; 
                anova_struct.([parnamepart{ch} '_' s{:} '_hands_IX'])= ES/(nanmean(FR(Par{3}==1 & idxS))+nanmean(FR(Par{3}==0 & idxS))); 
                
                h=anova_outs.p(3)<0.05;
                ES=nanmean(FR(crossed & idxS))-nanmean(FR(uncrossed & idxS));
                labelindexcr=h*sign(ES)+2; labelindexcr(isnan(labelindexcr))=2;
                anova_struct.([parnamepart{ch} '_' s{:} '_SxH'])=labelcr{labelindexcr}; 
                anova_struct.([parnamepart{ch} '_' s{:} '_SxH_ES'])=ES; 
                anova_struct.([parnamepart{ch} '_' s{:} '_SxH_IX'])=ES/(nanmean(FR(crossed & idxS))+nanmean(FR(uncrossed & idxS))); 
            end
        end
    end
end
end

function anova_struct=ph_pos_anova(keys,FR,States,choices,Position_idx,Positions,hands)
choices=choices+1;
hands=hands+1;
[~, ~, StateIndex]          =unique(States);
u_cho=unique(choices)';
u_hnd=unique(hands)';
INCHnamepart={'in','ch'};
LHRHnamepart={'NH','LH','RH'};

varnames={'State','position'};
states_for_multicomparison=keys.epoch_SxH_multicomp(ismember(keys.epoch_SxH_multicomp,keys.EPOCHS(:,1)));

multicomp_states=states_for_multicomparison(ismember(states_for_multicomparison,States'));
labelsp={'LS','-','RS'};

%% Independently for Instructed and Choice !!
for ch=u_cho
    for hn=u_hnd
        idx=choices==ch & hands==hn;
        [anova_out.p,anova_out.table,anova_out.stats,anova_out.terms] =...
        anovan(FR(idx),[StateIndex(idx),Position_idx(idx)],'model','full','varnames',varnames,'display',keys.display_anova_tables);
        anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_epoch_main_effect'])=anova_out.p(1)<0.05; %main effect on epoch!
        anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_position_main_effect'])=anova_out.p(2)<0.05; %main effect on epoch!
        for s=multicomp_states(:)'
            idx=ismember(States,s) & choices==ch & hands==hn;
            if any(~isnan(FR(idx)))
%                clear FRmean
%                 for p=1:max(Position_idx(idx))
%                     FRmean(p)=mean(FR(idx & Position_idx==p));
%                 end
%                 [~,FRmaxpos]=max(FRmean);
                idx_L= ismember(Position_idx,find(real(Positions)<0));
                idx_R= ismember(Position_idx,find(real(Positions)>0));
                
                %mean(FRmean(real(Positions(FRmaxpos))<0))*-1 + mean(FRmean(real(Positions(FRmaxpos))>0))
                [anova_outs.p] = anovan(FR(idx),{Position_idx(idx)},'model','full','varnames',{'Positions'},'display',keys.display_anova_tables);
                h=anova_outs.p(1)<0.05;
                %labelindexsp=h*sign(real(Positions(FRmaxpos)))+2; labelindexsp(isnan(labelindexsp))=2;
                labelindexsp=h*sign(nanmean(FR(idx & idx_R)) - nanmean(FR(idx & idx_L)))+2; labelindexsp(isnan(labelindexsp))=2;
                anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_position'])=labelsp{labelindexsp};
            else
                anova_struct.([INCHnamepart{ch} '_' LHRHnamepart{hn} '_' s{:} '_position'])=labelsp{2};
            end
        end
    end
end
end
