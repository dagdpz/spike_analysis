function population=ph_epochs(population,trials,keys)
% also reduces to only completed trials......

for u=1:numel(population)
    p=population(u);
    ut=ph_get_unit_trials(p,trials);
    u_types=unique([trials.type]);
    
    fields_to_shorten={'FR_average','stability_rating','SNR_rating','block','run','n','trial'};
    for f=1:numel(fields_to_shorten)
        population(u).(fields_to_shorten{f})= p.(fields_to_shorten{f})([ut.completed]);
    end
    for t=u_types
        EPOCHS=keys.EPOCHS_PER_TYPE{t};
        tr=[ut.type]==t & [ut.completed]; %% epochs only for completed trials............
        population(u).epochs_per_type{t}=NaN(sum(tr),size(EPOCHS,1));
        current_t=ut(tr);
        if ~isempty(current_t)
            for e=1:size(EPOCHS,1)
                s=EPOCHS{e,2};
                current_u=p.trial(tr);
                isthere         = find(arrayfun(@(x) any(x.states==s),current_t));
                whichstate      = arrayfun(@(x) find(x.states==s),current_t(isthere));
                state_onsets    = arrayfun(@(x,y) (x.states_onset(y)),current_t(isthere),whichstate);
                isthere         = isthere(~isnan(state_onsets));
                state_onsets    = state_onsets(~isnan(state_onsets));
                current_u       = current_u(isthere);
                from            = state_onsets + EPOCHS{e,3};
                to              = state_onsets + EPOCHS{e,4};
                population(u).epochs_per_type{t}(isthere,e)=arrayfun(@(x,y,z) single(sum(x.arrival_times>=y & x.arrival_times<=z)/(z-y)),current_u,from,to);
            end
        end
    end
    population(u).FR                =nanmean([population(u).FR_average]);
end
end