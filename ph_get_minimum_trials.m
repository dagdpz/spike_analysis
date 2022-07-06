function trial_criterion=ph_get_minimum_trials(keys,o,CM,UC,labels)
%% assign position indexes...
UC.pos_idx=1:max([o.trial.pos_index]);

% %% override labels
% for r=1:size(CM,1)
%     label='';
%     for c=1:size(CM,2)
%         add_to_label_index=0;
%         switch keys.condition_parameters{c}
%             case {'choice','reach_hand','perturbation','success','difficulty', 'stimuli_in_2hemifields'}
%                 add_to_label_index=1;
%         end
%         label_index=CM(r,c)+add_to_label_index;
%         to_add=keys.labels.(keys.condition_parameters{c}){label_index};
%         if ~isempty(to_add)
%             label=[label '_' to_add];
%         end
%     end
%     labels{r}=label(2:end);
% end

%% gotta find if and where reach hand condition is present, because of _trials_per_congruent_hand_hemifield
[~,reach_hand_idx]=ismember(keys.condition_parameters,'reach_hand');
reach_hand_idx=find(reach_hand_idx);

for hf=1:numel(UC.hemifield)
    trhf(hf,:)=[o.trial.hemifield]==UC.hemifield(hf);
end

for p=1:numel(UC.pos_idx)
    trpos(p,:)=[o.trial.pos_index]==UC.pos_idx(p);
end

for c=1:size(CM,1)
    clear trpar 
    for par=1:numel(keys.condition_parameters)
        fn=keys.condition_parameters{par};
        trpar(par,:)=[o.trial.(fn)]==CM(c,par);
    end
    namepart=labels{c};
    
    for hf=1:numel(UC.hemifield)
        trials_in_hemifield(hf)=sum(all([trpar; trhf(hf,:)],1));
    end
    for p=1:numel(UC.pos_idx)
        trials_in_position(p)=sum(all([trpar; trpos(p,:)],1));
    end
    if isempty(reach_hand_idx)
        trials_congruent=trials_in_hemifield;
    else
        ReachH=CM(c,reach_hand_idx);
        switch ReachH
            case 0 %% not too sure what that means at this level... ?
                hf=[];
            case 1
                hf=find(UC.hemifield==-1);
            case 2
                hf=find(UC.hemifield==1);
        end
        trials_congruent=sum(all([trpar; trhf(hf,:)],1));
    end    
    trial_criterion.([namepart '_trials_total_'])=sum(all(trpar,1));
    trial_criterion.([namepart '_trials_per_position_'])=min(trials_in_position); 
    trial_criterion.([namepart '_trials_per_hemifield_'])=min(trials_in_hemifield);
    trial_criterion.([namepart '_trials_per_congruent_hand_hemifield_'])=min(trials_congruent);
end
end