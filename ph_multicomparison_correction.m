function tuning_per_unit_table=ph_multicomparison_correction(tuning_per_unit_table,keys)
%% these are the ANOVA results we are looking for
% general ANOVAS are not corrected
% keys.AN.general_factors={'epoch_main','hemifield_main','hands_main','ExS','ExH','SxH','ExSxH'}; %factor for ANOVA: epoch*hand*space
% keys.AN.general_factors_per_hand={'position_main'};
% keys.AN.general_factors_per_hand_space={'PT_main','ExP','epoch_main'};
keys.AN.factors={'epoch','hemifield','hands','SxH','SglTar_Suc'};
%keys.AN.factors_per_hand={'epoch','position','fixation','PxF','RF_shift','gaze_modulation_x','gaze_modulation_y','RF_choice1','RF_choice2'};
keys.AN.factors_per_hand={'epoch','position','fixation','PxF','distance','angle','DxA','prefH', 'prefP','positionx','positiony','_positionxy','gaze_modulation_x','gaze_modulation_y','gaze_pref_x','gaze_pref_y'};
keys.AN.factors_per_hemifield={'Difficulty_Easy', 'Difficulty_Diff', 'SpatialComp_2HFTar', 'SpatialComp_1HFTar'};
keys.AN.factors_per_handhemifield={'epoch','prefH', 'prefP', 'SpatialComp_2HFTar', 'SpatialComp_1HFTar'}; %% basically perturbation and effector comparison

in_or_ch={'in','ch'};


factors=unique([keys.AN.factors,keys.AN.factors_per_hand,keys.AN.factors_per_hemifield,keys.AN.factors_per_handhemifield]);
default_factor_for_epochs='epoch';


tasktypes={};
tasktypes_index=1;
for effector=keys.cal.effectors
    for type=keys.cal.types
        [~, tasktypes{tasktypes_index}]=MPA_get_type_effector_name(type,effector);
        types(tasktypes_index)=type;
        tasktypes_index=tasktypes_index+1;
    end
end
cases=keys.position_and_plotting_arrangements;


% hands={'LH';'RH'};
% sides={'LS';'RS'};
% handhemifelds={'LH_LS';'LH_RS';'RH_LS';'RH_RS';};
% all_hands={'AH';'LH';'RH'};


hands={'IH';'CH'};
sides={'IS';'CS'};
handhemifelds={'IH_IS';'IH_CS';'CH_IS';'CH_CS';};
all_hands={'AH';'IH';'CH'};
TTlabelFNs=fieldnames(keys.TTlabels);
for f=1:numel(TTlabelFNs)
    FN=TTlabelFNs{f};
    for fe=1:numel(keys.TTlabels.(FN))
        current_field=keys.TTlabels.(FN){fe};
        current_field=strrep(current_field,'L','I');
        current_field=strrep(current_field,'R','C');
        keys.TTlabels.(FN){fe}=current_field;
    end
end


for t=1:numel(tasktypes)
    current_task=tasktypes{t};
    keys.anova_epochs=keys.ANOVAS_PER_TYPE(types(t));
    factor_epochs_defined=fieldnames(keys.anova_epochs);
    for c=1:numel(cases)
        current_case=cases{c};
        affix=[current_task '_' current_case(1:3)];
        %% standard factors...
        for f=1:numel(factors)
            factor=factors{f};
            if ismember(factor,factor_epochs_defined)
                epochs=keys.anova_epochs.(factor);
            else
                epochs=keys.anova_epochs.(default_factor_for_epochs);
                disp(['no anova epochs for factor:' factor ', assuming ' default_factor_for_epochs]);
            end
            for ch=1:numel(in_or_ch)
                INCH=in_or_ch{ch};
                
                
                if ismember(factor,keys.AN.factors) %% regardless of hand/space
                    prefixes={INCH};
                    tuning_per_unit_table=multicomp_correction(tuning_per_unit_table,factor,epochs,prefixes,affix,keys);
                end
                if ismember(factor,keys.AN.factors_per_hand) %% per hand (careful! AN (any hand does not need to be corrected for number of hands!)
                    prefixes={[INCH '_AH']};
                    tuning_per_unit_table=multicomp_correction(tuning_per_unit_table,factor,epochs,prefixes,affix,keys);
                    prefixes=strcat(INCH,'_',hands);
                    tuning_per_unit_table=multicomp_correction(tuning_per_unit_table,factor,epochs,prefixes,affix,keys);
                end
                if ismember(factor,keys.AN.factors_per_hemifield) %% per space
                    prefixes=strcat(INCH,'_',sides);
                    tuning_per_unit_table=multicomp_correction(tuning_per_unit_table,factor,epochs,prefixes,affix,keys);
                end
                if ismember(factor,keys.AN.factors_per_handhemifield) %% per hand & hemifield...!??!?
                    prefixes=strcat(INCH,'_AH_',sides);
                    tuning_per_unit_table=multicomp_correction(tuning_per_unit_table,factor,epochs,prefixes,affix,keys);
                    prefixes=strcat(INCH,'_',handhemifelds);
                    tuning_per_unit_table=multicomp_correction(tuning_per_unit_table,factor,epochs,prefixes,affix,keys);
                end
                if ~ismember(factor,keys.AN.factors) && ~ismember(factor,keys.AN.factors_per_hand)  && ~ismember(factor,keys.AN.factors_per_hemifield)
                    disp(['factor unknown:' factor]);
                end
            end
        end
        
        
        %% Need to add part where we combine results for two hemifields eventually
        epoch_multicomp=keys.anova_epochs.epoch;
        
        %% epoch tuning per hand (taking tuning for both HFs into account, 'bi' if opposing direction for both HFs)
        
        all_titles=tuning_per_unit_table(1,:);
        for ch=1:numel(in_or_ch)
            INCH=in_or_ch{ch};
            for hn=1:numel(all_hands)
                LHRH=all_hands{hn};
                for row=1:size(epoch_multicomp,1) %% might be able to remove "epoch" loop, but not too sure
                    if size(epoch_multicomp,2)==1
                    s=epoch_multicomp(row,1);
                    else
                    s=epoch_multicomp(row,2);
                    end
                    current_titles_L=[INCH '_' LHRH '_' sides{1} '_' s{:} '_epoch_' affix];
                    current_titles_R=[INCH '_' LHRH '_' sides{2} '_' s{:} '_epoch_' affix];
                    current_titles_C=[INCH '_' LHRH '_' s{:} '_epoch_' affix];
                    current_titles_B=[INCH '_' LHRH '_' s{:} '_epoch_bilateral_' affix];
                    [~,idx_L]=ismember(current_titles_L,all_titles);
                    [~,idx_R]=ismember(current_titles_R,all_titles);
                    [~,idx_C]=ismember(current_titles_C,all_titles);
                    [~,idx_B]=ismember(current_titles_B,all_titles);
                    if ~all([idx_L idx_R idx_C idx_B])
                        continue;
                    end
                    table_entries_L=tuning_per_unit_table(2:end,idx_L);
                    table_entries_R=tuning_per_unit_table(2:end,idx_R);
                    
                    % select label dependend on both sides - make sure they are both present?
                    combined_tuning=strcat(table_entries_L,table_entries_R);
                    combined_tuning(cellfun(@isempty,combined_tuning))={'--'};
                    labelindex=2*ones(size(combined_tuning));
                    labelindex(ismember(combined_tuning,{'-su','su-','susu'}))=1;
                    labelindex(ismember(combined_tuning,{'--'}))=2;
                    labelindex(ismember(combined_tuning,{'-en','en-','enen'}))=3;
                    labelindex(ismember(combined_tuning,{'suen','ensu'}))=4;
                    
                    tuning_per_unit_table(2:end,idx_C)=keys.TTlabels.epoch(labelindex);
                    labelindex=2*ones(size(combined_tuning));
                    labelindex(ismember(combined_tuning,{'susu'}))=1;
                    labelindex(ismember(combined_tuning,{'enen'}))=3;
                    tuning_per_unit_table(2:end,idx_B)=keys.TTlabels.epoch(labelindex);
                end
            end
        end
    end
end

end

function tuning_per_unit_table=multicomp_correction(tuning_per_unit_table,factor,epochs,prefixes,affix,keys)
labels=keys.labels;
TTlabels=keys.TTlabels;
choice_factors={'prefH','prefP'};
if size (epochs,2)>1 %% enhancment/suppression epochs
    epochs=epochs(:,2);
end
if isfield(TTlabels,factor)
    sig_labels=TTlabels.(factor);
elseif ismember(factor,choice_factors)
    sig_labels=TTlabels.choices; %% why choice(S!)?
else
    disp(['no designated TTlabels for factor:' factor ', assuming *false/true*']);
    sig_labels=TTlabels.true;
end
p_titles={''};
h_titles={''};
for p=1:numel(prefixes)
    p_titles=[p_titles; strcat(prefixes{p},'_',epochs,'_',factor,'_PV_',affix)];
    h_titles=[h_titles; strcat(prefixes{p},'_',epochs,'_',factor,'_',affix)];
end

all_titles=tuning_per_unit_table(1,:);
[~,idx_p]=ismember(p_titles,all_titles);
[~,idx_h]=ismember(h_titles,all_titles);

notthere=idx_p==0 |idx_h==0;
idx_p(notthere)=[];
idx_h(notthere)=[];
message_sent=0;
if ~isempty(idx_p) && ~isempty(idx_h)
    %% remove loop eventually for improving performance? not necessary i guess
    for u=2:size(tuning_per_unit_table,1)
        empties=cellfun(@isempty,tuning_per_unit_table(u,idx_p));
        pvals=[tuning_per_unit_table{u,idx_p}];
        if ~isempty(pvals)
            if (any(isnan(pvals)) || any(empties)) && ~message_sent
                disp(['found NAN or empty values for factor:' factor ', excluding them from multicomparison correction']);
                message_sent=1;
            end
            idx_h_current=idx_h(~empties);
            idx_h_current(isnan(pvals))=[];
            pvals(isnan(pvals))=[];
            sig_multi=correct_pvalue(pvals,keys);
            tuning_per_unit_table(u,idx_h_current)=sig_labels(sig_multi); %% i think this is correct...
        end
    end
end


end

function  sig_multi=correct_pvalue(pval,keys)
pval_abs=abs(pval);
pval_sign=sign(pval);
alpha=0.05;
%       h=[p_multi<0 && p_multi>-0.05, p_multi<-0.05 || p_multi>0.05, p_multi>0 && p_multi<0.05];
keys.AN.multicomparison='FDR'    ;
switch keys.AN.multicomparison
    case 'bonferoni'
        new_sig=pval*numel(pval)<alpha;
    case 'FDR' % Benjamini-Hochberg FDR algorithm
        [new_sig2]=fdr_bh(pval_abs);
        p_crit=(1:numel(pval))*1/numel(pval)*alpha;
        [pval_sorted, sort_idx]=sort(pval_abs);
        new_sig=pval_sorted-p_crit<=0;
        if all(new_sig)
            N_sig=numel(pval);
        elseif ~any(new_sig) 
            N_sig=0;
        else
            N_sig=find(diff([new_sig 0])==-1,1,'last');
        end
        new_sig=[ones(1,N_sig) zeros(1,numel(pval)-N_sig)];
        new_sig(sort_idx)=new_sig;
        aaaa=1;
end
sig_multi=pval_sign.*new_sig+2;
end
