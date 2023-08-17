function population=ph_epochs(population,trials,keys)
for u=1:numel(population)
    p=population(u);
    pt=ph_get_unit_trials(p,trials);
    u_types=unique([trials.type]);
    
    for t=u_types
        EPOCHS=keys.EPOCHS_PER_TYPE{t};
        tr=[pt.type]==t;
        population(u).epochs_per_type{t}=NaN(sum(tr),size(EPOCHS,1));
        for e=1:size(EPOCHS,1)
            s=EPOCHS{e,2};
            %             si=ismember({p.trial(tr).states},s);
            current_t=pt(tr);
            current_u=p.trial(tr);
            
            isthere=find(arrayfun(@(x) any(x.states==s),current_t));
            whichstate=arrayfun(@(x) find(x.states==s),current_t(isthere));  
%             if isempty(whichstate)
%                 continue;
%             end
            state_onsets=arrayfun(@(x,y) (x.states_onset(y)),current_t(isthere),whichstate);
            
            isthere=isthere(~isnan(state_onsets));
            state_onsets=state_onsets(~isnan(state_onsets));
            current_u=current_u(isthere);
            from=state_onsets + EPOCHS{e,3};
            to=  state_onsets + EPOCHS{e,4};
            population(u).epochs_per_type{t}(isthere,e)=arrayfun(@(x,y,z) single(sum(x.arrival_times>=y & x.arrival_times<=z)/(z-y)),current_u,from,to);
        end
    end
    population(u).FR                =nanmean([population(u).FR_average]);
end
%
%
%     function FR=compute_epoch_FR(trial,unit)
%         si=ismember(trial.states,s);
%         if any(si) && ~isnan(trial.states_onset(si))
%             current_state_onset=trial.states_onset(si);
%             t_before_state=EPOCHS{e,3} + current_state_onset;
%             t_after_state=EPOCHS{e,4}  + current_state_onset;
%             FR=single(sum(unit.arrival_times>=t_before_state & unit.arrival_times<=t_after_state)/(t_after_state-t_before_state));
%         else
%             FR=single(NaN);
%         end
%     end
end