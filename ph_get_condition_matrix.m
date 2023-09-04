function [UC, CM, labels]=ph_get_condition_matrix(trials,keys)
%% only consider accepted and completed trials
tr_con=[trials.accepted] & [trials.completed]; 
trials=trials(tr_con);

%% type and effectors
utype=unique([trials.type]); 
ueffector=unique([trials.effector]);
all_type_effectors      = combvec(utype,ueffector)';
type_effectors =[];
if isfield(keys,'conditions_to_plot')
    % redefine type_effectors to include only relevant
    for t=1:size(all_type_effectors,1)
        typ=all_type_effectors(t,1);
        eff=all_type_effectors(t,2);
        [~, type_effector_short{t}]=MPA_get_type_effector_name(typ,eff);
        if ~ismember(type_effector_short{t},keys.conditions_to_plot) %|| sum(tr_con)<1
            continue;
        end
        type_effectors=[type_effectors; all_type_effectors(t,:)];
    end
    type_effector_short(~ismember(type_effector_short,keys.conditions_to_plot))=[];
else
    type_effectors=all_type_effectors;
    type_effector_short={''};
end

UC.type     =unique(type_effectors(:,1))';
UC.effector =unique(type_effectors(:,2))';
UC.type_effector=type_effectors;
UC.type_effector_labels=type_effector_short;

UC.completed=unique([trials.completed]);
UC.perturbation         =unique([trials.perturbation]);
for c=1:numel(keys.condition_parameters)
    UC.(keys.condition_parameters{c})=unique([trials.(keys.condition_parameters{c})]);
end

%% postitions 
UC.position             =unique(vertcat(trials.position),'rows');
UC.hemifield            =unique([trials.hemifield]);
UC.fix_index            =unique([trials.fix_index]);
UC.pos_index            =unique([trials.pos_index]);
fixations_temp=unique(vertcat(trials.fixation),'rows');
fix_temp_idx=true(size(fixations_temp,1),1);
for x=1:size(fixations_temp,1)-1
    if any(all(abs(bsxfun(@minus,fixations_temp(x+1:end,:),fixations_temp(x,:)))<keys.cal.precision_fix,2))
        fix_temp_idx(x)=false;
    end
end
UC.fixation=fixations_temp(fix_temp_idx,:);


for c=1:numel(keys.condition_parameters)
    par=keys.condition_parameters{c};
    CM_cell{c}=UC.(par);
end
CM=combvec(CM_cell{:})';

%% NOW, since we know the conditions, we can also define legends and colors (legend defines fields in keys.color.)!
for r=1:size(CM,1)
    label='';
    for c=1:size(CM,2)
        add_to_label_index=0;
        switch keys.condition_parameters{c}
            case {'choice','reach_hand','perturbation','success','difficulty', 'stimuli_in_2hemifields'}
                add_to_label_index=1;
        end
        label_index=CM(r,c)+add_to_label_index;
        to_add=keys.labels.(keys.condition_parameters{c}){label_index};
        if ~isempty(to_add)
            label=[label '_' to_add];
        end
    end
    labels{r}=label(2:end);
end
