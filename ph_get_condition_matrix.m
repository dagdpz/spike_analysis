function [UC, CM, labels]=ph_get_condition_matrix(trials,keys)
types       =[trials.type];
effectors   =[trials.effector];
all_type_effectors      = combvec(unique(types),unique(effectors))';
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

hands       =[trials.reach_hand];
choice      =[trials.choice];
perturbation=[trials.perturbation];
hemifield   =[whatisthis.trial.hemifield];

position=vertcat(whatisthis.trial.position);
fixations_temp=unique(vertcat(whatisthis.trial.fixation),'rows');
fix_temp_idx=true(size(fixations_temp,1),1);
for x=1:size(fixations_temp,1)-1
    if any(all(abs(bsxfun(@minus,fixations_temp(x+1:end,:),fixations_temp(x,:)))<4,2)) %% precision....
        fix_temp_idx(x)=false;
    end
end
UC.fixation=fixations_temp(fix_temp_idx,:);
perturbation(ismember(perturbation, keys.cal.perturbation_groups{1}))=0;
perturbation(ismember(perturbation, keys.cal.perturbation_groups{2}))=1;
perturbation(isnan(perturbation))=0;

UC.hemifield        =unique(hemifield); 
UC.position         =unique(position,'rows');
UC.reach_hand       =unique(hands);
UC.choice           =unique(choice);
UC.perturbation     =unique(perturbation);
% UC.perturbation     =UC.perturbation(~isnan(UC.perturbation));

%% limit conditions key?
if ~any(keys.tt.hands==0) % cause hands 0 is any hand
    UC.reach_hand     =UC.reach_hand(ismember(UC.reach_hand,keys.tt.hands));
end
UC.choice    =UC.choice(ismember(UC.choice,keys.tt.choices));
UC.reach_hand=UC.reach_hand;

%% KK stuff
UC.stimulustype         =unique([trials.stimulustype]);
UC.difficulty           =unique([trials.difficulty]);
UC.success              =unique([trials.success]);
UC.fix_index            =unique([whatisthis.trial.fix_index]);

%CM_cell={};
for c=1:numel(keys.condition_parameters)
    CM_cell{c}=UC.(keys.condition_parameters{c});
    %Conditons_per_trial(:,c)=[trials.(keys.condition_parameters{c})];
end
%CM=unique(Conditons_per_trial,'rows');
CM=combvec(CM_cell{:})';

%% NOW, since we know the conditions, we can also define legends and colors (legend defines fields in keys.color.)!
for r=1:size(CM,1)
    label='';
    for c=1:size(CM,2)
        add_to_label_index=0;
        switch keys.condition_parameters{c}
            case {'choice','reach_hand','perturbation','success','difficulty'}
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