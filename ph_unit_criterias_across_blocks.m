function  pop_resorted = ph_unit_criterias_across_blocks(pop_resorted,trials,keys)
[trials.position]=deal(0);
[trials.fixation]=deal(0);
[trials.hemifield]=deal(0);
[trials.fix_index]=deal(0);
[trials.pos_index]=deal(0);
[trials.accepted]=deal(true);
[UC, CM, labels]=ph_get_condition_matrix(trials,keys);
CP  =[keys.condition_parameters];
CM=CM'; %% something is fishy here
for u=1:numel(pop_resorted)
    U=pop_resorted(u);
    N_spikes_per_trial=arrayfun(@(x) numel(x.arrival_times),U.trial);
    T=ph_get_unit_trials(U,trials);%
    acc=[T.accepted];
    
    FR=[T.FR_average];
    block=[T.block];
    tdur=arrayfun(@(x) x.states_onset(x.states==90)-x.states_onset(x.states==2),T);
    
    
    for t=1:numel(UC.type)
        typ=UC.type(t);
        [~, ~, typ_label, ~]=MPA_get_type_effector_name(typ,0);
        
        ut_typ=[T.type]==typ;
        trcon=true(size(ut_typ));
        for c=1:size(CM,2)
            for par=1:size(CM,1)
                fn=CP{par};
                trcon(par,:)=[T.(fn)]==CM(par,c);
            end
            tr=all(trcon,1) & ut_typ & acc;
            
            u_blocks=unique(block(tr));
            % resampling per block
            spikes_per_block=[];
            stability_per_block=[];
            FR_per_block=[];            
            FRrs=[];
            brs=[];
            for b=1:numel(u_blocks)
                bl=u_blocks(b);
                trb=tr & block==bl;
                FRb=FR(trb);
                tb=double(cumsum(tdur(trb)));
                if numel(tb)>=2
                    FRrsb = ph_resample_FRs(FRb,tb);
                else
                    FRrsb= FRb;
                end
                FRrs=[FRrs FRrsb];
                brs=[brs repmat(bl,size(FRrsb))];
                sta=[T(trb).stability_rating];
                stability_per_block(b)=unique(sta(~isnan(sta)));
                spikes_per_block(b)=sum(N_spikes_per_trial(trb));
                FR_per_block(b)=mean(FRrsb);
            end
                        
            valid_blocks=u_blocks;
            valid_sta=stability_per_block;
            valid_nsp=spikes_per_block;
            valid_FRs=FR_per_block;
            p=0;
            while p<0.001 && numel(valid_blocks)>1
                % ANOVA to find main effect of block
                p = anova1(FRrs,brs,'off');
                if p<0.001
                    if any(valid_nsp<max(valid_nsp)/10)
                        [~,ix]=min(valid_nsp);
                        exclusion_reason=3;
                    elseif numel(valid_blocks)>2
                        [~,ix]=max(abs(valid_FRs-mean(valid_FRs)));
                        exclusion_reason=4;
                    else
                        [~,ix]=min(valid_sta);
                        exclusion_reason=5;
                    end
                    
                    block_to_remove=valid_blocks(ix);
                    valid_blocks(ix)=[];
                    valid_nsp(ix)=[];
                    valid_sta(ix)=[];
                    valid_FRs(ix)=[];
                    ix=brs==block_to_remove;
                    FRrs(ix)=[];
                    brs(ix)=[];
                    pop_resorted(u).exclusion_reason(block==block_to_remove)=exclusion_reason;
                    %tr(~ismember(block,valid_blocks))=false;
                end
            end
            %% set stability to NaN and accepted to false
            to_set_false=ut_typ & all(trcon,1) & ~ismember(block,valid_blocks);
            pop_resorted(u).accepted(to_set_false)=false;
            %% calculate per condition
            for_average=ut_typ & all(trcon,1) & ismember(block,valid_blocks);
            
            stability= nanmean(FRrs)/std(FRrs); % Fano-factor: variance / mean
            WFs_cat=vertcat(U.trial(for_average).waveforms);
            amps=max(abs(WFs_cat),[],2);
            WF_rescaled=WFs_cat./repmat(amps,1,size(WFs_cat,2));
            snr=1/mean(std(WF_rescaled,0,1));
            pop_resorted(u).(['criteria_stability_' typ_label(1) '_' labels{c}])         = stability;
            pop_resorted(u).(['criteria_SNR_' typ_label(1) '_' labels{c}])               = snr;
            
        end
    end
    
    for_average=pop_resorted(u).accepted;
    if sum(for_average)>=2
        FRrs = ph_resample_FRs(FR(for_average),cumsum(tdur(for_average)));
    else
        FRrs=FR(for_average);
    end
    stability= nanmean(FRrs)/std(FRrs) ; % Fano-factor: variance / mean
    
    WFs_cat=vertcat(U.trial(for_average).waveforms);
    amps=max(abs(WFs_cat),[],2);
    WF_rescaled=WFs_cat./repmat(amps,1,size(WFs_cat,2));
    snr=1/mean(std(WF_rescaled,0,1));
    
    if keys.cal.automatic_stablity
        if ~isempty(stability)
            pop_resorted(u).avg_stability=single(stability);
        else
            pop_resorted(u).avg_stability=single(NaN);
        end
    end
    if keys.cal.automatic_SNR
        if ~isempty(snr)
            pop_resorted(u).avg_SNR=single(snr);
        else
            pop_resorted(u).avg_SNR=single(NaN);
        end
    end
end
end