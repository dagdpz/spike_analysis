function [state_names,absolute_state_onsets,relative_state_onsets,relative_epochs,epoch_names,states]=ph_state_onsets_2(trial,trial_onsets,relative_state,keys)
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
for s=1:numel(unique_states)
    sta=unique_states(s);
    rel_idx=unique_states==relative_state;
    relative_state_onsets(s)=nanmean([trial_onsets(:,s)-trial_onsets(:,rel_idx); NaN]);
    absolute_state_onsets(s)=nanmean([trial_onsets(:,s); NaN]);
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