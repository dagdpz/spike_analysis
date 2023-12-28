function  pop_resorted = ph_unit_criterias_across_blocks(pop_resorted,trials,keys)
[trials.position]=deal(0);
[trials.fixation]=deal(0);
[trials.hemifield]=deal(0);
[trials.fix_index]=deal(0);
[trials.pos_index]=deal(0);
[trials.accepted]=deal(true);
[UC, CM, labels]=ph_get_condition_matrix(trials,keys);
%CP  =[{'effector'},keys.condition_parameters,{'pos_index'}];
CP  =[keys.condition_parameters];
CM=CM'; %% something is fishy here
for u=1:numel(pop_resorted)
    U=pop_resorted(u);    
    T=ph_get_unit_trials(U,trials);%     
     acc=[T.accepted];
    
    FR=[T.FR_average];
    block=[T.block];
    
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
            FR_per_block=[];
            for b=1:numel(u_blocks)
                bl=u_blocks(b);
                FR_per_block(b)=mean(FR(block==bl));
            end
            
            valid_blocks=u_blocks;
            valid_FRs=FR_per_block;
            p=0;
            while p<0.05 && numel(unique(block(tr)))>1
                % ANOVA to find main effect of block
                %p = anova1(FR(tr),block(tr),'off');
                p = kruskalwallis(FR(tr),block(tr),'off');
                if p<0.05
                    % remove block with lowest FR
                    valid_blocks(valid_FRs==min(valid_FRs))=[];
                    valid_FRs(valid_FRs==min(valid_FRs))=[];
                    tr(~ismember(block,valid_blocks))=false;
                    
                end
            end
            %% set stability to NaN and accepted to false
            to_set_false=ut_typ & all(trcon,1) & ~ismember(block,valid_blocks);
            pop_resorted(u).stability_rating(to_set_false)=NaN;
            pop_resorted(u).accepted(to_set_false)=false;
            %% calculate per condition 
            for_average=ut_typ & all(trcon,1) & ismember(block,valid_blocks);
            pop_resorted(u).criteria.(['stability_' typ_label(1) '_' labels{c}])         =nanmean(pop_resorted(u).stability_rating(for_average)); %% recompute based on FR per trial (?)
            pop_resorted(u).criteria.(['SNR_' typ_label(1) '_' labels{c}])               =nanmean(pop_resorted(u).SNR_rating(for_average));      %% recompute based on waveforms per trial (?)
            %pop_resorted(u).(['single_rating_' typ_label(1) '_' labels{n}])     =nanmean(pop_resorted(u).(stability)(for_average));
        end
    end
    pop_resorted(u).avg_stability=nanmean(pop_resorted(u).stability_rating);
    pop_resorted(u).avg_SNR=nanmean(pop_resorted(u).SNR_rating);
end
end