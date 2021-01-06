function [state_names,absolute_state_onsets,relative_state_onsets,relative_epochs,epoch_names,states]=ph_state_onsets(trial,relative_state,keys)
%% this is used for plotting
global MA_STATES
FN=fieldnames(MA_STATES);
field_idx=structfun(@(x) numel(x)==1,MA_STATES);
NAMES=lower(FN(field_idx));
VALUES=struct2cell(MA_STATES);
VALUES=cell2mat(VALUES(field_idx));
unique_states=unique([trial.states]);
epoch_states=[keys.EPOCHS{:,2}];
relative_epochs=[];
epoch_names={};

tr_matrix=true(numel(trial),numel(unique_states)); %this line is still like 1/3 of the plotting time (down from 86%) ...
for t=1:numel(trial)
    if any(trial(t).states==relative_state)
        tr_matrix(t,:)=ismember(unique_states,trial(t).states) ;
    else
        tr_matrix(t,:)=false(size(unique_states));
    end
    %=any(trial(t).states==sta) && any(trial(t).states==relative_state);
end

for s=1:numel(unique_states)
    sta=unique_states(s);
    tr=tr_matrix(:,s); % trial(arrayfun(@(x) any(x.states==sta) && any(x.states==relative_state),trial));
    onsets=[trial(tr).states_onset];
    sta_idx=[trial(tr).states]==sta;
    rel_idx=[trial(tr).states]==relative_state;
    relative_state_onsets(s)=nanmean([onsets(sta_idx)-onsets(rel_idx) NaN]);
    absolute_state_onsets(s)=nanmean([onsets(sta_idx) NaN]);
    state_names(s)=NAMES(find(VALUES==sta,1,'first'));
    states(s)=sta;
    epoch_idx=epoch_states==sta;
    if any(epoch_idx)
        relative_epochs=[relative_epochs;vertcat(keys.EPOCHS{epoch_idx,3})+relative_state_onsets(s) vertcat(keys.EPOCHS{epoch_idx,4})+relative_state_onsets(s)];
        epoch_names=[epoch_names; keys.EPOCHS(epoch_idx,1)];
    else
        relative_epochs=[relative_epochs;NaN NaN];
        epoch_names=[epoch_names; {''}];
    end
end
end