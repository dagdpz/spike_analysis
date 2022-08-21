function [UC, CM, labels]=ph_get_condition_matrix(trials,keys)

%% define conditions to look at
tr_con=[trials.accepted] & [trials.completed]; 
trials=trials(tr_con);
%% type and effectors
%% limit type and effector to plot HERE! --> NOT IN RUN_STATE_ALIGNMENT ANY MORE
utype=unique([trials.type]); 
%utype=utype(ismember(utype,keys.cal.types)); %% dont like that this one is called types with an 's'
ueffector=unique([trials.effector]);
%ueffector=ueffector(ismember(ueffector,keys.cal.effectors)); %% dont like that this one is called effectors with an 's'
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



for c=1:numel(keys.condition_parameters)
    UC.(keys.condition_parameters{c})=unique([trials.(keys.condition_parameters{c})]);
end

%% this part is not ideal (the perturbation one!)
perturbation=[trials.perturbation];
perturbation(ismember(perturbation, keys.cal.perturbation_groups{1}))=0;
perturbation(ismember(perturbation, keys.cal.perturbation_groups{2}))=1;
perturbation(isnan(perturbation))=0;
UC.perturbation         =unique(perturbation);

%% postitions are defiend by arrangement!
[whatisthis]=ph_arrange_positions_and_plots(keys,trials);
UC.position             =unique(vertcat(whatisthis.trial.position),'rows');
UC.hemifield            =unique([whatisthis.trial.hemifield]);
UC.fix_index            =unique([whatisthis.trial.fix_index]);
UC.pos_index            =unique([whatisthis.trial.pos_index]);
fixations_temp=unique(vertcat(whatisthis.trial.fixation),'rows');
fix_temp_idx=true(size(fixations_temp,1),1);
for x=1:size(fixations_temp,1)-1
    if any(all(abs(bsxfun(@minus,fixations_temp(x+1:end,:),fixations_temp(x,:)))<keys.cal.precision_fix,2))
        fix_temp_idx(x)=false;
    end
end
UC.fixation=fixations_temp(fix_temp_idx,:);

% %if ~any(keys.tt.hands==0) % cause hands 0 is any hand ... 
% UC.reach_hand     =UC.reach_hand(ismember(UC.reach_hand,keys.tt.hands));
% %end
% UC.choice    =UC.choice(ismember(UC.choice,keys.tt.choices));
%% up to here


for c=1:numel(keys.condition_parameters)
    par=keys.condition_parameters{c};
%     if ~all(isnan(keys.cal.(par))) && ~isempty(keys.cal.(par))%% limit conditions key?
%         UC.(par)    =UC.(par)(ismember(UC.(par),keys.cal.(par)));
%     end
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
