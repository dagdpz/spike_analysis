function keys=ph_tuning_table_correction(keys)



TT=keys.tuning_table;
TT=ph_multicomparison_correction(TT,keys);
[TT]=ph_extend_tuning_table(TT,keys);
[TT, keys.selection_title]=ph_reduce_tuning_table(TT,keys);
keys.tuning_table=TT;
end

%% multicomparison and subfunctions
function TT=ph_multicomparison_correction(TT,keys)
in_or_ch={'in','ch'};
hands={'IH';'CH'};
sides={'IS';'CS'};
handhemifelds={'IH_IS';'IH_CS';'CH_IS';'CH_CS';};
all_hands={'AH';'IH';'CH'};

factors=unique([keys.AN.factors,keys.AN.factors_per_hand,keys.AN.factors_per_hemifield,keys.AN.factors_per_handhemifield]);
default_factor_for_epochs='epoch';
to_display={'Multicomparison issues:'};


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
    keys.anova_epochs=keys.AN.multicomp_epochs(types(t));
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
                to_display=[to_display, {['no anova epochs for factor:' factor ', assuming ' default_factor_for_epochs]}];
            end
            for ch=1:numel(in_or_ch)
                INCH=in_or_ch{ch};
                if ismember(factor,keys.AN.factors) %% regardless of hand/hemifield
                    prefixes={INCH};
                    [TT,to_display]=multicomp_correction(TT,to_display,factor,epochs,prefixes,affix,keys);
                end
                if ismember(factor,keys.AN.factors_per_hand) %% per hand (careful! AN (any hand does not need to be corrected for number of hands!)
                    prefixes={[INCH '_AH']};
                    [TT,to_display]=multicomp_correction(TT,to_display,factor,epochs,prefixes,affix,keys);
                    prefixes=strcat(INCH,'_',hands);
                    [TT,to_display]=multicomp_correction(TT,to_display,factor,epochs,prefixes,affix,keys);
                end
                if ismember(factor,keys.AN.factors_per_hemifield) %% per hemifield
                    prefixes=strcat(INCH,'_',sides);
                    [TT,to_display]=multicomp_correction(TT,to_display,factor,epochs,prefixes,affix,keys);
                end
                if ismember(factor,keys.AN.factors_per_handhemifield) %% per hand & hemifield...!??!?
                    prefixes=strcat(INCH,'_AH_',sides);
                    [TT,to_display]=multicomp_correction(TT,to_display,factor,epochs,prefixes,affix,keys);
                    prefixes=strcat(INCH,'_',handhemifelds);
                    [TT,to_display]=multicomp_correction(TT,to_display,factor,epochs,prefixes,affix,keys);
                end
                if ~ismember(factor,keys.AN.factors) && ~ismember(factor,keys.AN.factors_per_hand)  && ~ismember(factor,keys.AN.factors_per_hemifield)
                    to_display=[to_display, {['factor unknown:' factor]}];
                end
            end
        end
        
        
        %% Need to add part where we combine results for two hemifields eventually
        %% epoch tuning per hand (taking tuning for both HFs into account, 'bi' if opposing direction for both HFs)
        epoch_multicomp=keys.anova_epochs.epoch;
        
        all_titles=TT(1,:);
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
                    current_titles_L =[INCH '_' LHRH '_' sides{1} '_' s{:} '_epoch_' affix];
                    current_titles_R =[INCH '_' LHRH '_' sides{2} '_' s{:} '_epoch_' affix];
                    current_titles_l =[INCH '_' LHRH '_' sides{1} '_' s{:} '_epoch_DF_' affix];
                    current_titles_r =[INCH '_' LHRH '_' sides{2} '_' s{:} '_epoch_DF_' affix];
                    current_titles_C =[INCH '_' LHRH '_' s{:} '_epoch_' affix];
                    %current_titles_B =[INCH '_' LHRH '_' s{:} '_epoch_bilateral_' affix]; %% what was this even supposed to be?
                    [~,idx_L]=ismember(current_titles_L,all_titles);
                    [~,idx_R]=ismember(current_titles_R,all_titles);
                    [~,idx_l]=ismember(current_titles_l,all_titles);
                    [~,idx_r]=ismember(current_titles_r,all_titles);
                    [~,idx_C]=ismember(current_titles_C,all_titles);
                    %[~,idx_B]=ismember(current_titles_B,all_titles);
                    if ~all([idx_L idx_R]) %%~all([idx_L idx_R idx_C idx_B])
                        continue;
                    end
                    if idx_C==0 %&& idx_B==0 %%%%isempty(idx_C) && isempty(idx_B) %% either both are there (legacy) or not there yet
                       idx_C=size(TT,2)+1;%idx_B=size(TT,2)+2;
                    end
                        
                    table_entries_L=TT(2:end,idx_L);
                    table_entries_R=TT(2:end,idx_R);
                    table_entries_l=[TT{2:end,idx_l}];
                    table_entries_r=[TT{2:end,idx_r}];
                    [~,max_diff_idx]=max(abs([table_entries_l;table_entries_r])); %% had forgotten absolute here... -_-
                    table_entry_LR=[table_entries_L,table_entries_R];
                    table_entry_D=table_entry_LR([1:numel(max_diff_idx)]+(max_diff_idx-1)*numel(max_diff_idx))';
                    [~,labelindex_D]=ismember(table_entry_D,{'su','-','en'});
                    
                    % select label dependend on both sides - make sure they are both present?
                    combined_tuning=strcat(table_entries_L,table_entries_R);
                    combined_tuning(cellfun(@isempty,combined_tuning))={'--'};
                    labelindex=2*ones(size(combined_tuning));
                    labelindex(ismember(combined_tuning,{'-su','su-','susu'}))=1;
                    labelindex(ismember(combined_tuning,{'--'}))=2;
                    labelindex(ismember(combined_tuning,{'-en','en-','enen'}))=3;
                    %labelindex(ismember(combined_tuning,{'suen','ensu'}))=4;
                    labelindex(ismember(combined_tuning,{'suen','ensu'}))=labelindex_D(ismember(combined_tuning,{'suen','ensu'}));
                    
                    TT{1,idx_C}=current_titles_C;
                    TT(2:end,idx_C)=keys.TTlabels.epoch(labelindex);
%                     labelindex=2*ones(size(combined_tuning));
%                     labelindex(ismember(combined_tuning,{'susu'}))=1;
%                     labelindex(ismember(combined_tuning,{'enen'}))=3;
%                     TT{1,idx_B}=current_titles_B;
%                     TT(2:end,idx_B)=keys.TTlabels.epoch(labelindex);
                    
                    if strcmp(LHRH,'AH')
                        
                        current_titles_E =[INCH '_' s{:} '_epoch_' affix];
                        current_titles_H =[INCH '_' s{:} '_hemifield_' affix];
                        [~,idx_E]=ismember(current_titles_E,all_titles);
                        [~,idx_H]=ismember(current_titles_H,all_titles);
                        table_entries_H=TT(2:end,idx_H);
                        table_entries_E=TT(2:end,idx_E);
                        to_replace=[false; ismember(table_entries_H,{'IS','CS'}) & ismember(table_entries_E,{'-'})];
                        TT(to_replace,idx_E)=TT(to_replace,idx_C);
                    end
                end
            end
        end
    end
end

% display errors and such
to_display=unique(to_display);
for td=1:numel(to_display)
    disp(to_display{td});
end


end

function [TT,to_display]=multicomp_correction(TT,to_display,factor,epochs,prefixes,affix,keys)
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
    to_display=[to_display, {['no designated TTlabels for factor:' factor ', assuming *false/true*']}];
    sig_labels=TTlabels.true;
end
p_titles={''};
h_titles={''};
for p=1:numel(prefixes)
    p_titles=[p_titles; strcat(prefixes{p},'_',epochs,'_',factor,'_PV_',affix)];
    h_titles=[h_titles; strcat(prefixes{p},'_',epochs,'_',factor,'_',affix)];
end

all_titles=TT(1,:);
[~,idx_p]=ismember(p_titles,all_titles);
[~,idx_h]=ismember(h_titles,all_titles);

notthere=idx_p==0 |idx_h==0;
idx_p(notthere)=[];
idx_h(notthere)=[];
message_sent=0;
if ~isempty(idx_p) && ~isempty(idx_h)
    %% remove loop eventually for improving performance? not necessary i guess
    for u=2:size(TT,1)
        empties=cellfun(@isempty,TT(u,idx_p));
        pvals=[TT{u,idx_p}];
        if ~isempty(pvals) && ~all(isnan(pvals))
            if (any(isnan(pvals)) || any(empties)) && ~message_sent                
                to_display=[to_display, {['found NAN or empty values for factor:' factor ', excluding them from multicomparison correction']}];
                message_sent=1;
            end
            idx_h_current=idx_h(~empties);
            idx_h_current(isnan(pvals))=[];
            pvals(isnan(pvals))=[];
            sig_multi=correct_pvalue(pvals,keys);
            TT(u,idx_h_current)=sig_labels(sig_multi); %% i think this is correct...
        end
    end
end
end

function  sig_multi=correct_pvalue(pval,keys)
pval_abs=abs(pval);
pval_sign=sign(pval);
alpha=0.05;
switch keys.AN.multicomparison
    case 'none'
        new_sig=abs(pval)<alpha;
    case 'bonferoni'
        new_sig=pval*numel(pval)<alpha;
    case 'FDR' % Benjamini-Hochberg FDR algorithm
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
end
sig_multi=pval_sign.*new_sig+2;
end

%% extended TT entries subfunctions

function TT=ph_extend_tuning_table(TT,keys)

epochs={'Fhol','Cue','Cue2','EDel','Del','MemE','MemL','PreS','PreR','PeriS','PeriR','Pre2','Peri2','PreG','CueG','TIhol','THol'};
INCH=keys.tt.IC_for_criterion;

for t=1:numel(keys.tt.tasktypes)
    taskcase=keys.tt.tasktypes{t};
    idx.hemifield   =DAG_find_column_index(TT,[INCH '_hemifield_main_'   taskcase] );
    idx.epoch   =DAG_find_column_index(TT,[INCH '_epoch_main_'     taskcase] );
    idx.hands   =DAG_find_column_index(TT,[INCH '_hands_main_'     taskcase] );
    idx.SxH     =DAG_find_column_index(TT,[INCH '_SxH_'            taskcase] );
    idx.ExS     =DAG_find_column_index(TT,[INCH '_ExS_'            taskcase] );
    idx.ExH     =DAG_find_column_index(TT,[INCH '_ExH_'            taskcase] );
    idx.ESH     =DAG_find_column_index(TT,[INCH '_ExSxH_'          taskcase] );
    
    for e=1:numel(epochs)
        idx.(epochs{e})                     =DAG_find_column_index(TT,[INCH '_' epochs{e} '_epoch_'   taskcase] );
        idx.([epochs{e} '_AH'])             =DAG_find_column_index(TT,[INCH '_AH_' epochs{e} '_epoch_'   taskcase] );
        idx.([epochs{e} '_hemifield_DF'])   =DAG_find_column_index(TT,[INCH '_' epochs{e} '_hemifield_DF_'   taskcase] );
        idx.([epochs{e} '_IS_FR'])          =DAG_find_column_index(TT,[INCH '_AH_IS_' epochs{e} '_epoch_FR_'  taskcase] );
        idx.([epochs{e} '_CS_FR'])          =DAG_find_column_index(TT,[INCH '_AH_CS_' epochs{e} '_epoch_FR_'  taskcase] );
        idx.([epochs{e} '_IS_EN'])          =DAG_find_column_index(TT,[INCH '_AH_IS_' epochs{e} '_epoch_DF_'  taskcase] );
        idx.([epochs{e} '_CS_EN'])          =DAG_find_column_index(TT,[INCH '_AH_CS_' epochs{e} '_epoch_DF_'  taskcase] );
        idx.([epochs{e} '_in_hemifield_DF'])=DAG_find_column_index(TT,['in_' epochs{e} '_hemifield_DF_'   taskcase] );
        idx.([epochs{e} '_ch_hemifield_DF'])=DAG_find_column_index(TT,['ch_' epochs{e} '_hemifield_DF_'   taskcase] );
        idx.([epochs{e} '_in_IS_FR'])       =DAG_find_column_index(TT,['in_AH_IS_' epochs{e} '_epoch_FR_'  taskcase] );
        idx.([epochs{e} '_in_CS_FR'])       =DAG_find_column_index(TT,['in_AH_CS_' epochs{e} '_epoch_FR_'  taskcase] );
        idx.([epochs{e} '_ch_IS_FR'])       =DAG_find_column_index(TT,['ch_AH_IS_' epochs{e} '_epoch_FR_'  taskcase] );
        idx.([epochs{e} '_ch_CS_FR'])       =DAG_find_column_index(TT,['ch_AH_CS_' epochs{e} '_epoch_FR_'  taskcase] );
        
        idx.([epochs{e} '_in_IH_hemifield'])=DAG_find_column_index(TT,['in_IH_' epochs{e} '_hemifield_'  taskcase] );
        idx.([epochs{e} '_in_CH_hemifield'])=DAG_find_column_index(TT,['in_CH_' epochs{e} '_hemifield_'  taskcase] );
        idx.([epochs{e} '_in_IS_hands'])    =DAG_find_column_index(TT,['in_IS_' epochs{e} '_hands_'  taskcase] );
        idx.([epochs{e} '_in_CS_hands'])    =DAG_find_column_index(TT,['in_CS_' epochs{e} '_hands_'  taskcase] );
    end
    
    %% adding columns for visual, motor, visuomotor, and fixation cells
    if any(idx.Fhol) && any(idx.PeriS) && any(idx.Cue) && any(idx.TIhol)
        nc=size(TT,2);
        TT{1,nc+1}=['fixation_only_' taskcase];
        TT{1,nc+2}=['fixation_and_sac_suppression_' taskcase];
        TT{1,nc+3}=['Sac_supression_' taskcase];
        TT(2:end,nc+1)=num2cell(ismember(TT(2:end,idx.Fhol),'en') & ismember(TT(2:end,idx.PeriS),'-')  & ~ismember(TT(2:end,idx.Cue),'en') & ~ismember(TT(2:end,idx.TIhol),'en'));
        TT(2:end,nc+2)=num2cell(ismember(TT(2:end,idx.Fhol),'en') & ismember(TT(2:end,idx.PeriS),'su') & ~ismember(TT(2:end,idx.Cue),'en') & ~ismember(TT(2:end,idx.TIhol),'en'));
        TT(2:end,nc+3)=num2cell(ismember(TT(2:end,idx.Fhol),'-')  & ismember(TT(2:end,idx.PeriS),'su') & ~ismember(TT(2:end,idx.Cue),'en') & ~ismember(TT(2:end,idx.TIhol),'en'));
    end
    
    TT=get_VM(TT,idx.Cue,idx.TIhol,'first',['visual_' taskcase]);
    TT=get_VM(TT,idx.Cue,idx.TIhol,'second',['motor_' taskcase]);
    TT=get_VM(TT,idx.Cue,idx.TIhol,'both',['visuomotor_' taskcase]);
    TT=get_VM(TT,idx.Cue,idx.TIhol,'none',['notclassified_' taskcase]);
    
    TT=get_VM(TT,idx.Cue,idx.PreS,'first',['visual_preS_' taskcase]);
    TT=get_VM(TT,idx.Cue,idx.PreS,'second',['motor_preS_' taskcase]);
    TT=get_VM(TT,idx.Cue,idx.PreS,'both',['visuomotor_preS_' taskcase]);
    
    TT=get_VM(TT,idx.Cue,idx.Pre2,'first',['visual_pre2_' taskcase]);
    TT=get_VM(TT,idx.Cue,idx.Pre2,'second',['motor_pre2_' taskcase]);
    TT=get_VM(TT,idx.Cue,idx.Pre2,'both',['visuomotor_pre2_' taskcase]);
    
    TT=get_VM(TT,idx.Cue,idx.PeriS,'first',['visual_periS_' taskcase]);
    TT=get_VM(TT,idx.Cue,idx.PeriS,'second',['motor_periS_' taskcase]);
    TT=get_VM(TT,idx.Cue,idx.PeriS,'both',['visuomotor_periS_' taskcase]);
    
    
    TT=get_VMI(TT,idx.PeriS_IS_FR,idx.Cue_IS_FR,['VMI_peri_IS_' taskcase],'signed');
    TT=get_VMI(TT,idx.PeriS_CS_FR,idx.Cue_CS_FR,['VMI_peri_CS_' taskcase],'signed');
    TT=get_VMI(TT,idx.PeriS_IS_EN,idx.Cue_IS_EN,['VMI_periEN_IS_' taskcase],'absolute');
    TT=get_VMI(TT,idx.PeriS_CS_EN,idx.Cue_CS_EN,['VMI_periEN_CS_' taskcase],'absolute');
    
    TT=get_VMI(TT,idx.TIhol_IS_FR,idx.Cue_IS_FR,['VMI_post_IS_' taskcase],'signed');
    TT=get_VMI(TT,idx.TIhol_CS_FR,idx.Cue_CS_FR,['VMI_post_CS_' taskcase],'signed');
    TT=get_VMI(TT,idx.TIhol_IS_EN,idx.Cue_IS_EN,['VMI_postEN_IS_' taskcase],'absolute');
    TT=get_VMI(TT,idx.TIhol_CS_EN,idx.Cue_CS_EN,['VMI_postEN_CS_' taskcase],'absolute');
    
    %% MP
    
%     TT=get_VM(TT,idx.Cue,idx.Del,'first',['visual_only_' taskcase]); %% doesnt use suppression though !
%     TT=get_VM(TT,idx.Cue,idx.Del,'both',['visuomotor_' taskcase]); %% doesnt use suppression though ! --> does use it for this one Oo
%     TT=get_VM(TT,idx.Cue,idx.Del,'second',['motor_only_' taskcase]); %% doesnt use suppression though !
%     TT=get_VM(TT,idx.Cue,idx.Del,'first',['visual_only_' taskcase]); %% doesnt use suppression though !
%     TT=get_VM(TT,idx.Cue,idx.Del,'first',['visual_only_' taskcase]); %% doesnt use suppression though !
    
    
    nc=size(TT,2);
    if any(idx.Cue) % && any(idx.TIhol)
        nc=nc+1;
        TT{1,nc}=['visual_en_' taskcase];
        TT(2:end,nc)=num2cell(ismember(TT(2:end,idx.Cue),{'en','bi'}));% & ~ismember(TT(2:end,idx.TIhol),{'en','su','bi'}));
    end
    
    if  any(idx.Del)
        nc=nc+1;
        TT{1,nc}=['motor_en_' taskcase];
        TT(2:end,nc)=num2cell(ismember(TT(2:end,idx.Del),{'en','bi'})); %& ismember(TT(2:end,idx.Del),{'en','su','bi'}));
    end
    %% MP end
    
    nc=size(TT,2);
    if any(idx.PeriS_IS_EN) && any(idx.Cue_IS_EN) && any(idx.PeriS_CS_EN) && any(idx.Cue_CS_EN)
        nc=nc+1;
        TT{1,nc}=['VMI_periEN_' taskcase];
        TT(2:end,nc)=num2cell((abs(cell2mat(TT(2:end,idx.PeriS_IS_EN))) + abs(cell2mat(TT(2:end,idx.PeriS_CS_EN))) - ...
            abs(cell2mat(TT(2:end,idx.Cue_IS_EN))) - abs(cell2mat(TT(2:end,idx.Cue_CS_EN))))./...
            (abs(cell2mat(TT(2:end,idx.PeriS_IS_EN))) + abs(cell2mat(TT(2:end,idx.PeriS_CS_EN))) + ...
            abs(cell2mat(TT(2:end,idx.Cue_IS_EN))) + abs(cell2mat(TT(2:end,idx.Cue_CS_EN)))));
    end
        
    if any(idx.TIhol_IS_EN) && any(idx.Cue_IS_EN) && any(idx.TIhol_CS_EN) && any(idx.Cue_CS_EN)
        nc=nc+1;
        TT{1,nc}=['VMI_postEN_' taskcase];
        TT(2:end,nc)=num2cell((abs(cell2mat(TT(2:end,idx.TIhol_IS_EN))) + abs(cell2mat(TT(2:end,idx.TIhol_CS_EN))) - ...
            abs(cell2mat(TT(2:end,idx.Cue_IS_EN))) - abs(cell2mat(TT(2:end,idx.Cue_CS_EN))))./...
            (abs(cell2mat(TT(2:end,idx.TIhol_IS_EN))) + abs(cell2mat(TT(2:end,idx.TIhol_CS_EN))) + ...
            abs(cell2mat(TT(2:end,idx.Cue_IS_EN))) + abs(cell2mat(TT(2:end,idx.Cue_CS_EN)))));
    end
        
    if any(idx.CueG_IS_EN) && any(idx.PreG_IS_EN)
        nc=nc+1;
        TT{1,nc}=['PreCueSum_IS_' taskcase];
        TT(2:end,nc)=num2cell(cell2mat(TT(2:end,idx.CueG_IS_EN)) + cell2mat(TT(2:end,idx.PreG_IS_EN)));
    end
    
    if any(idx.CueG_CS_EN) && any(idx.PreG_CS_EN)
        nc=nc+1;
        TT{1,nc}=['PreCueSum_CS_' taskcase];
        TT(2:end,nc)=num2cell(cell2mat(TT(2:end,idx.CueG_CS_EN)) + cell2mat(TT(2:end,idx.PreG_CS_EN)));
    end
        
    if any(idx.CueG_IS_FR) && any(idx.PreG_IS_FR)
        nc=nc+1;
        TT{1,nc}=['PreCueMean_IS_' taskcase];
        TT(2:end,nc)=num2cell((cell2mat(TT(2:end,idx.CueG_IS_FR)) + cell2mat(TT(2:end,idx.PreG_IS_FR)))/2);
    end
    
    if any(idx.CueG_CS_FR) && any(idx.PreG_CS_FR)
        nc=nc+1;
        TT{1,nc}=['PreCueMean_CS_' taskcase];
        TT(2:end,nc)=num2cell((cell2mat(TT(2:end,idx.CueG_CS_FR)) + cell2mat(TT(2:end,idx.PreG_CS_FR)))/2); %% /2 !! for mean
    end
    
    for e=1:numel(epochs)
        if any(idx.([epochs{e} '_in_hemifield_DF'])) && any(idx.([epochs{e} '_in_IS_FR'])) && any(idx.([epochs{e} '_in_CS_FR'])) && any(idx.([epochs{e} '_ch_IS_FR'])) && any(idx.([epochs{e} '_ch_CS_FR']))
            nc=nc+1;
            TT{1,nc}=['in_' epochs{e} '_hemifieldCI_IX_' taskcase];
            TT(2:end,nc)=num2cell(cell2mat(TT(2:end,idx.([epochs{e} '_in_hemifield_DF'])))./...
                (cell2mat(TT(2:end,idx.([epochs{e} '_in_IS_FR']))) + cell2mat(TT(2:end,idx.([epochs{e} '_in_CS_FR']))) +...
                cell2mat(TT(2:end,idx.([epochs{e} '_ch_IS_FR']))) + cell2mat(TT(2:end,idx.([epochs{e} '_ch_CS_FR']))))); %% /2 !! for mean
        end
    end
    
    for e=1:numel(epochs)
        if any(idx.([epochs{e} '_ch_hemifield_DF'])) && any(idx.([epochs{e} '_in_IS_FR'])) && any(idx.([epochs{e} '_in_CS_FR'])) && any(idx.([epochs{e} '_ch_IS_FR'])) && any(idx.([epochs{e} '_ch_CS_FR']))
            nc=nc+1;
            TT{1,nc}=['ch_' epochs{e} '_hemifieldCI_IX_' taskcase];
            TT(2:end,nc)=num2cell(cell2mat(TT(2:end,idx.([epochs{e} '_ch_hemifield_DF'])))./...
                (cell2mat(TT(2:end,idx.([epochs{e} '_in_IS_FR']))) + cell2mat(TT(2:end,idx.([epochs{e} '_in_CS_FR']))) +...
                cell2mat(TT(2:end,idx.([epochs{e} '_ch_IS_FR']))) + cell2mat(TT(2:end,idx.([epochs{e} '_ch_CS_FR']))))); %% /2 !! for mean
        end
        
        if any(idx.([epochs{e} '_in_IH_hemifield'])) && any(idx.([epochs{e} '_in_CH_hemifield'])) && any(idx.([epochs{e} '_in_IS_hands'])) && any(idx.([epochs{e} '_in_CS_hands']))
            
            IH_CS=strcmp(TT(:,idx.([epochs{e} '_in_IH_hemifield'])),'CS');
            IH_IS=strcmp(TT(:,idx.([epochs{e} '_in_IH_hemifield'])),'IS');
            CH_CS=strcmp(TT(:,idx.([epochs{e} '_in_CH_hemifield'])),'CS');
            CH_IS=strcmp(TT(:,idx.([epochs{e} '_in_CH_hemifield'])),'IS');
            IS_CH=strcmp(TT(:,idx.([epochs{e} '_in_IS_hands'])),'CH');
            IS_IH=strcmp(TT(:,idx.([epochs{e} '_in_IS_hands'])),'IH');
            CS_CH=strcmp(TT(:,idx.([epochs{e} '_in_CS_hands'])),'CH');
            CS_IH=strcmp(TT(:,idx.([epochs{e} '_in_CS_hands'])),'IH');
            incongruent_hemifield=  (IH_CS & CH_IS) | (IH_IS & CH_CS);
            incongruent_hands=  (IS_CH & CS_IH) | (IS_IH & CS_CH);
            
            CS= (IH_CS | CH_CS) & ~incongruent_hemifield;
            IS= (IH_IS | CH_IS) & ~incongruent_hemifield;
            CH= (IS_CH | CS_CH) & ~incongruent_hands;
            IH= (IS_IH | CS_IH) & ~incongruent_hands;
            
            nc=nc+1;
            TT(CS,nc)={'CS'};
            TT(IS,nc)={'IS'};
            TT(~(CS|IS),nc)={'-'};
            TT(incongruent_hemifield,nc)={'incongruent'};
            TT{1,nc}=['in_' epochs{e} '_hemifield_perhand_' taskcase];
            
            nc=nc+1;
            TT(CH,nc)={'CH'};
            TT(IH,nc)={'IH'};
            TT(~(CH|IH),nc)={'-'};
            TT(incongruent_hands,nc)={'incongruent'};
            TT{1,nc}=['in_' epochs{e} '_hands_perhemifield_' taskcase];
        end
    end
end

%% across tasks combinations
idx.CueS_hemifieldperhand=DAG_find_column_index(TT,'in_Cue_hemifield_perhand_Ddsa_han');
idx.CueR_hemifieldperhand=DAG_find_column_index(TT,'in_Cue_hemifield_perhand_Ddre_han');
idx.preS_hemifieldperhand=DAG_find_column_index(TT,'in_PreS_hemifield_perhand_Ddsa_han');
idx.preR_hemifieldperhand=DAG_find_column_index(TT,'in_PreR_hemifield_perhand_Ddre_han');
idx.preS_handperhemifield=DAG_find_column_index(TT,'in_PreS_hands_perhemifield_Ddsa_han');
idx.preR_handperhemifield=DAG_find_column_index(TT,'in_PreR_hands_perhemifield_Ddre_han');

TT=get_across_task(TT,idx.CueS_hemifieldperhand,idx.CueR_hemifieldperhand,'hemifield','in_Cue_hemifield_perhand_Ddre_or_Ddsa');
TT=get_across_task(TT,idx.preS_hemifieldperhand,idx.preR_hemifieldperhand,'hemifield','in_Pre_hemifield_perhand_Ddre_or_Ddsa');
TT=get_across_task(TT,idx.preS_handperhemifield,idx.preR_handperhemifield,'hands','in_Pre_hands_perhemifield_Ddre_or_Ddsa');

combined_column=cell(size(TT,1),1);
for props=2:size(keys.tt.combine_tuning_properties,2)
    column_index=DAG_find_column_index(TT,keys.tt.combine_tuning_properties{1,props});
    if isempty(column_index)
        disp(['column not found: ' keys.tt.combine_tuning_properties{1,props}]);
        continue;
    end
    combined_column=cellfun(@(x,y) cat(2,x,num2str(y)),combined_column,TT(:,column_index),'uniformoutput',false);
end
if ~isempty(keys.tt.combine_tuning_properties)
    combined_column{1,1}=keys.tt.combine_tuning_properties{1};
else
    combined_column{1,1}='combined_effects';
end
TT(:,end+1)=combined_column;


for props=1:size(keys.tt.replace_tuning,1)
    column_index=DAG_find_column_index(TT,keys.tt.replace_tuning{props,1});
    if isempty(column_index)
        disp(['column not found: ' keys.tt.replace_tuning{props,1}]);
        continue;
    end
    to_replace=cellfun(@(x) strcmp(num2str(x),num2str(keys.tt.replace_tuning{props,2})),TT(:,column_index));
    TT(to_replace,column_index)=repmat({num2str(keys.tt.replace_tuning{props,3})},sum(to_replace),1);
end
end

function TT=get_VM(TT,idx1,idx2,modulation,field_name)
if any(idx1) && any(idx2)
    nc=size(TT,2)+1;
    idx1_modulated=ismember(TT(2:end,idx1),{'en','su','bi'});
    idx2_modulated=ismember(TT(2:end,idx2),{'en','su','bi'});
    switch modulation
        case 'first'
            to_add=num2cell( idx1_modulated &~idx2_modulated);
        case 'second'
            to_add=num2cell(~idx1_modulated & idx2_modulated);
        case 'both'
            to_add=num2cell( idx1_modulated & idx2_modulated);
        case 'none'
            to_add=num2cell(~idx1_modulated &~idx2_modulated);
    end
    TT{1,nc}=field_name;
    TT(2:end,nc)= to_add;
end
end

function TT=get_VMI(TT,idx1,idx2,field_name,mode)
if any(idx1) && any(idx2)
    nc=size(TT,2)+1;
    idx1_EN=cell2mat(TT(2:end,idx1));
    idx2_EN=cell2mat(TT(2:end,idx2));
    switch mode
        case 'absolute'
            to_add=num2cell((abs(idx1_EN) - abs(idx2_EN))./(abs(idx1_EN) + abs(idx2_EN)));
        case 'signed'
            to_add=num2cell((idx1_EN - idx2_EN)./(idx1_EN + idx2_EN));
    end
    TT{1,nc}=field_name;
    TT(2:end,nc)= to_add;
end
end

function TT=get_across_task(TT,idx1,idx2,mode,field_name)
nc=size(TT,2)+1;
switch mode
    case 'hands'
        entries={'CH','IH'};
    case 'hemifield'
        entries={'CS','IS'};
end

if any(idx1) && any(idx1)
    IX1C=strcmp(TT(:,idx1),entries{1});
    IX1I=strcmp(TT(:,idx1),entries{2});
    IX2C=strcmp(TT(:,idx2),entries{1});
    IX2I=strcmp(TT(:,idx2),entries{2});
    incongruent=  (IX1C & IX2I) | (IX1I & IX2C);
    Contra= (IX1C | IX2C) & ~incongruent;
    Ipsi= (IX1I | IX2I) & ~incongruent;
    TT(Contra | Ipsi,nc)={'tuned'};
    TT(~(Contra|Ipsi),nc)={'-'};
    TT(incongruent,nc)={'incongruent'};
    TT{1,nc}=field_name;
end
end

%% reducing TT entries 

function [TT Sel_for_title]=ph_reduce_tuning_table(TT,keys)
%% monkey
idx_unitID=DAG_find_column_index(TT,'unit_ID');
if numel(keys.monkey)>3
    TT=TT([true; ~cellfun(@isempty,strfind(TT(2:end,idx_unitID),keys.monkey(1:3)))],:);
end

%% new stuff to make everything consistent within each other
for c=1:numel(keys.condition_parameters)
    %% hands --> reach_hand; choices --> choice; perturbations --> perturbation;
    CM_cell{c}=keys.tt.(keys.condition_parameters{c});
end
CM=combvec(CM_cell{:})';
labels={};
for t=1:numel(keys.tt.tasktypes)
    tasktype=keys.tt.tasktypes{t};
    if ~strfind(tasktype,'_')
        disp('keys.tt.tasktypes needs to contain arrangement as well to work properly, f.e.: Ddre_han')
    end
    for r=1:size(CM,1)
        label='';
        for c=1:size(CM,2)
            add_to_label_index=0;
            switch keys.condition_parameters{c}
                case {'choice','reach_hand','perturbation','success','difficulty', 'stimuli_in_2hemifields'}
                    add_to_label_index=1;
            end
            label_index=CM(r,c)+add_to_label_index;
            if ~isnan(label_index)
                to_add=keys.labels.(keys.condition_parameters{c}){label_index};
                if ~isempty(to_add)
                    label=[label '_' to_add];
                end
            else
                label=[label '*']; %% this is essentially the new part for allowing to not consider a certain parameter (?)
            end
        end
        if strcmp(label(1),'*') % not nice, but okay for now. if first is a star, dont remove that one
            labels{r+(t-1)*size(CM,1)}=[label '_' tasktype];
        else
            labels{r+(t-1)*size(CM,1)}=[label(2:end) '_' tasktype];
        end
    end
end

%labels=unique(labels); %? if you do this, conditions dont match any more!
%we should be able to simply remove condition (?) BUT we want either or in * cases

%% trial criterion exclusion looking for labels !
row_index=true(size(TT,1)-1,1);
TT_titles=TT(1,:);
for l=1:numel(labels) %% crashes for not existing monkey
    column_title=['existing_' labels{l}];
    starpos=[0 strfind(column_title,'*') numel(column_title)+1];
    column_index=true(size(TT_titles));
    for s=1:numel(starpos)-1 % this is the loop looking for name parts... (NOT IDEAL if conditions are named similarly)
        column_index=column_index & cellfun(@(x) any(strfind(x,column_title(starpos(s)+1:starpos(s+1)-1))),TT_titles);
    end
    if sum(column_index)>0
        row_index=row_index & any(cell2mat(TT(2:end,column_index)),2);
    else
        disp([keys.condition_parameters{:} labels{l} 'not existing'])
    end
end
TT=TT([true;row_index],:);
TT=ph_target_reassign(keys,TT);


%% selection (pick only specific entries)
if size(TT,1)>1
    row_index=true(size(TT,1)-1,1);
    
    % ratings
    selection_criteria={'stability_rating','SNR_rating','Single_rating'};
    for sel=1:numel(selection_criteria)
        criterion=selection_criteria{sel};
        crit_min=min(keys.tt.(criterion));
        crit_max=max(keys.tt.(criterion));
        column_index=DAG_find_column_index(TT,criterion);
        row_index=row_index & cell2mat(TT(2:end,column_index))>=crit_min & cell2mat(TT(2:end,column_index))<=crit_max;
    end
    
    % specific selection
    for sel=1:size(keys.tt.selection,1)
        column_index=DAG_find_column_index(TT,keys.tt.selection{sel,1});
        if isempty(column_index)
            disp(['column not found: ' keys.tt.selection{sel,1}]);
            continue;
        end
        if ischar(keys.tt.selection{sel,2})
            TT(2:end,column_index)=cellfun(@num2str,TT(2:end,column_index),'UniformOutput',0); %% crazy cellstring bug
            row_index=row_index & ~cellfun(@isempty,strfind(TT(2:end,column_index),keys.tt.selection{sel,2}));
            
        else %% is number!?
            row_index=row_index & cell2mat(TT(2:end,column_index))==keys.tt.selection{sel,2};
        end
    end
    TT=TT([true;row_index],:);
end

%% unselect (remove only specific entries)
if size(TT,1)>1
    row_index=true(size(TT,1)-1,1);
    for sel=1:size(keys.tt.unselect,1)
        column_index=DAG_find_column_index(TT,keys.tt.unselect{sel,1});
        if isempty(column_index)
            disp(['column not found: ' keys.tt.unselect{sel,1}]);
            continue;
        end
        if ischar(keys.tt.unselect{sel,2})
            row_index=row_index & cellfun(@(x) isempty(strfind(x,keys.tt.unselect{sel,2})),TT(2:end,column_index));
        else %% is number!?
            row_index=row_index & cell2mat(TT(2:end,column_index))~=keys.tt.unselect{sel,2};
        end
    end
    TT=TT([true;row_index],:);
    if ~isempty(keys.tt.unselected_list)
        load([keys.tt.unselected_list{:} filesep 'list']);
        TT=TT([true; ~ismember(TT(2:end,idx_unitID),example_list)],:);
    end
end

%% Selection title (making sure files are named differently)
if ~isempty(keys.tt.selection)
    Sel_for_title=[keys.tt.selection(:,1) repmat({' = '},size(keys.tt.selection(:,1)))...
        cellfun(@num2str,keys.tt.selection(:,2),'uniformoutput',false)  repmat({', '},size(keys.tt.selection(:,1)))]';
else
    Sel_for_title={'';'';'';''};
end;
if ~isempty(keys.tt.unselect)
    Sel_for_title=[Sel_for_title, [keys.tt.unselect(:,1) repmat({' ~= '},size(keys.tt.unselect(:,1)))...
        cellfun(@num2str,keys.tt.unselect(:,2),'uniformoutput',false)  repmat({', '},size(keys.tt.unselect(:,1)))]'];
end;

end

%% reassign target (to be able to have hemifields defined relative to injection site?)
function TT=ph_target_reassign(keys,TT)

col_index_ref = find(strcmp(TT(1,:),keys.contra_ipsi_relative_to));
col_index_target = find(strcmp(TT(1,:),'target'));
for i = 2:size(TT,1)
    ref = char(TT(i,col_index_ref));
    target = char(TT(i,col_index_target));
    if strcmp(keys.contra_ipsi_relative_to,'perturbation_site')
        target_right=strcmpi(target(end-1:end),'_R') || strcmpi(target(end-1:end),'_r');
        ref_right=strcmpi(ref(end-1:end),'_R') || strcmpi(ref(end-1:end),'_r');
        
        if (target_right && ref_right) || (~target_right && ~ref_right)
            TT(i,col_index_target) = {[target(1:end-1) 'L']}; %% left becomes ipsi!
        elseif (target_right && ~ref_right) || (~target_right && ref_right)
            TT(i,col_index_target) = {[target(1:end-1) 'R']}; %% Right becomes contra!
        end
    end
end
end