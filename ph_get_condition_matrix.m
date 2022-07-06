function [UC, CM, labels]=ph_get_condition_matrix(trials,keys)
%% type and effectors
all_type_effectors      = combvec(unique([trials.type]),unique([trials.effector]))';
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


%% define conditions to look at
tr_con=ismember([trials.completed],keys.cal.completed);
[whatisthis]=ph_arrange_positions_and_plots(keys,trials(tr_con));

for c=1:numel(keys.condition_parameters)
    UC.(keys.condition_parameters{c})=unique([trials.(keys.condition_parameters{c})]);
end

%% this part is not ideal
perturbation=[trials.perturbation];
perturbation(ismember(perturbation, keys.cal.perturbation_groups{1}))=0;
perturbation(ismember(perturbation, keys.cal.perturbation_groups{2}))=1;
perturbation(isnan(perturbation))=0;
UC.perturbation     =unique(perturbation);
UC.position             =unique(vertcat(whatisthis.trial.position),'rows');
UC.hemifield            =unique([whatisthis.trial.hemifield]);
UC.fix_index            =unique([whatisthis.trial.fix_index]);
fixations_temp=unique(vertcat(whatisthis.trial.fixation),'rows');
fix_temp_idx=true(size(fixations_temp,1),1);
for x=1:size(fixations_temp,1)-1
    if any(all(abs(bsxfun(@minus,fixations_temp(x+1:end,:),fixations_temp(x,:)))<keys.cal.precision_fix,2))
        fix_temp_idx(x)=false;
    end
end
UC.fixation=fixations_temp(fix_temp_idx,:);
%% limit conditions key?
if ~any(keys.tt.hands==0) % cause hands 0 is any hand
    UC.reach_hand     =UC.reach_hand(ismember(UC.reach_hand,keys.tt.hands));
end
UC.choice    =UC.choice(ismember(UC.choice,keys.tt.choices));
%% up to here


for c=1:numel(keys.condition_parameters)
    CM_cell{c}=UC.(keys.condition_parameters{c});
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
