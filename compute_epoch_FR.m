function FR=compute_epoch_FR(trial,unit)
        si=ismember(trial.states,s);
        if any(si)
            current_state_onset=trial.states_onset(si);
            t_before_state=EPOCHS{s,3} + current_state_onset;
            t_after_state=EPOCHS{s,4}  + current_state_onset;
            FR=sum(unit.arrival_times>=t_before_state & unit.arrival_times<=t_after_state)/(t_after_state-t_before_state);
        else
            FR=NaN;
        end
    end