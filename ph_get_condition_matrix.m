function [UC, CM]=ph_get_condition_matrix(trials,keys)
types       =[trials.type];
effectors   =[trials.effector];
all_type_effectors      = combvec(unique(types),unique(effectors))';
type_effectors =[];

% redifine type_effectors to include only relevant
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
UC.types     =unique(type_effectors(:,1))';
UC.effectors =unique(type_effectors(:,2))';


%% define conditions to look at
tr_con=ismember([trials.completed],keys.cal.completed);
[whatisthis]=ph_arrange_positions_and_plots(keys,trials(tr_con));

types       =[trials.type];
effectors   =[trials.effector];
hands       =[trials.reach_hand];
choice      =[trials.choice];
perturbation=[trials.perturbation];
hemifield   =[whatisthis.trial.hemifield];
perturbation(ismember(perturbation, keys.cal.perturbation_groups{1}))=0;
perturbation(ismember(perturbation, keys.cal.perturbation_groups{2}))=1;

UC.hemifields=unique(hemifield); %[-1,0,1]; % why does this have to be hardcoded? ---> Because case not defined yet, case defines positions !!

UC.hands     =unique(hands);
UC.choice    =unique(choice);
UC.perturbation    =unique(perturbation);
UC.perturbation=UC.perturbation(~isnan(UC.perturbation));

%% limit conditions key?
if ~any(keys.tt.hands==0) % cause hands 0 is any hand
    UC.hands     =UC.hands(ismember(UC.hands,keys.tt.hands));
end
UC.choice    =UC.choice(ismember(UC.choice,keys.tt.choices));



condition_parameters  ={'reach_hand','choice','perturbation'};

CM_cell={};
for c=1:numel(condition_parameters)
    CM_cell{c}=UC.(condition_parameters{c});
end


condition_matrix            = combvec(UC.hands,UC.choice, UC.perturbation,UC.hemifields)';
conditions_out              = combvec(UC.effectors,UC.hands,UC.choice, UC.perturbation)';