function [keys all_sta]=ph_get_epoch_keys(keys,typ,eff,by_effector)
global MA_STATES
        get_expected_MA_states(typ,eff,by_effector); %% does it make sense to distinguish by effector? Probably not any more!
        keys.ALL_EPOCHS=keys.EPOCHS_PER_TYPE{typ};
        keys.EPOCHS=keys.ALL_EPOCHS(ismember([keys.ALL_EPOCHS{:,2}],MA_STATES.all_states),:);
        keys.PSTH_WINDOWS=keys.WINDOWS_PER_TYPE{typ};
        keys.ANOVAS=keys.ANOVAS_PER_TYPE(typ);
        all_sta=find(ismember([keys.EPOCHS{:,2}],MA_STATES.all_states));        
        %get_expected_MA_states(typ,eff,0); %% does it make sense to distinguish by effector?
end