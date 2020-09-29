function population=ph_epochs(population,keys)
for u=1:numel(population)
    for t=1:numel(population(u).trial)
        population(u).trial(t).epoch=struct();
        population(u).trial(t).FR=[];
        trial=population(u).trial(t);
        keys.EPOCHS=keys.EPOCHS_PER_TYPE{trial.type};
        states_for_epochs=[keys.EPOCHS{:,2}];
        
        % for M2S states are not unique !!
        [~,tr_state_idx]                =unique(trial.states,'last');
        trial_states                    =trial.states(tr_state_idx);
        trial_states_onset              =trial.states_onset(tr_state_idx); 
        trial.FR                        =sum(trial.arrival_times>=0 & trial.arrival_times<=max(trial_states_onset(1:end-1)))/max(trial_states_onset(1:end-1));
                
        for s=1:numel(states_for_epochs);
            current_state=states_for_epochs(s);
            trial.epoch(s).state=keys.EPOCHS{s,1};
            if ismember(current_state,trial_states) 
                state_idx=trial_states==current_state;
                current_state_onset=trial_states_onset(state_idx);
                trial.epoch(s).state_onset=current_state_onset; %% needed?
                t_before_state=keys.EPOCHS{s,3} + current_state_onset;
                t_after_state=keys.EPOCHS{s,4}  + current_state_onset;
                
                % spike count based FR!
                trial.epoch(s).FR=sum(trial.arrival_times>=t_before_state & trial.arrival_times<=t_after_state)/(t_after_state-t_before_state);
               
            else % for states that dont exist in this trial, trial.epoch(s).state_onset and trial.epoch(s).FR will be NaNs
                trial.epoch(s).state_onset=NaN;
                trial.epoch(s).FR=NaN;
            end
        end
        population(u).trial(t)=trial;
        
%         mov_idxs=find(~isnan(Movement.ini_all) & ~isnan(Movement.end_all) & ...
%             Movement.ini_all>=MA_out.states(t).start_obs); %% & before trial end...?, minimum amplitude?
%         
%         for m=1:numel(mov_idxs)
%             mov_idx=mov_idxs(m);
%             ini_t=Movement.ini_all(mov_idx);
%             end_t=Movement.end_all(mov_idx);
%             trial(t).unit(c,u).movement(m).endpos       =Movement.endpos_all(mov_idx);
%             trial(t).unit(c,u).movement(m).startpos     =Movement.startpos_all(mov_idx);
%             trial(t).unit(c,u).movement(m).vector       =trial(t).unit(c,u).movement(m).endpos - trial(t).unit(c,u).movement(m).startpos;
%             
%             start_state_idx=find(diff([MA_out.states(t).MP_states_onset<=ini_t])~=0);
%             end_state_idx=find(diff([MA_out.states(t).MP_states_onset>=end_t])~=0);
%             state_indexes=start_state_idx:end_state_idx;
%             
%             trial(t).unit(c,u).movement(m).tar_at_start=MA_out.states(t).MP_states(state_indexes(1))==MA_out.states(t).state_1ao;
%             trial(t).unit(c,u).movement(m).tar_at_end=MA_out.states(t).MP_states(state_indexes(end))==MA_out.states(t).state_1ao;
%             trial(t).unit(c,u).movement(m).tar_crossed=numel(state_indexes)>2; %% works so far, because saccade can only be crossing one target if there are two state transitions within the saccade;
%             
%             trial(t).unit(c,u).movement(m).shape_at_start   = NaN;
%             trial(t).unit(c,u).movement(m).shape_at_end     = NaN;
%             trial(t).unit(c,u).movement(m).shape_crossed    = NaN;
%             if trial(t).unit(c,u).movement(m).tar_at_start
%                 tar_ins=sum(MA_out.states(t).MP_states(1:state_indexes(1))==MA_out.states(t).state_1ao);
%                 trial(t).unit(c,u).movement(m).shape_at_start=Movement.all_convexities(max([1,Movement.targets_inspected(tar_ins)]));
%             end
%             
%             if trial(t).unit(c,u).movement(m).tar_crossed
%                 tar_ins=sum(MA_out.states(t).MP_states(1:state_indexes(end))==MA_out.states(t).state_1ao);
%                 trial(t).unit(c,u).movement(m).shape_crossed=Movement.all_convexities(max([1,Movement.targets_inspected(tar_ins)]));
%             end
%             
%             if trial(t).unit(c,u).movement(m).tar_at_end
%                 tar_ins=sum(MA_out.states(t).MP_states(1:state_indexes(end))==MA_out.states(t).state_1ao);
%                 trial(t).unit(c,u).movement(m).shape_at_end=Movement.all_convexities(max([1,Movement.targets_inspected(tar_ins)]));
%             end
%             
%             t_b=-0.2;
%             t_a= 0.2;
%             arrival_times_mov       =MA_out.physiology(t).spike_arrival_times{c,u}-ini_t;
%             %arrival_times_mov       =arrival_times_mov(arrival_times_mov>=(t_b-6*keys.PSTH_binwidth) & arrival_times_mov<(t_a+6*keys.PSTH_binwidth));
%             %
%             PSTH_x_axis=t_b:keys.PSTH_binwidth:t_a;
%             trial(t).unit(c,u).movement(m).arrival_times=arrival_times_mov(arrival_times_mov>=t_b & arrival_times_mov<=t_a);%/keys.PSTH_binwidth);
%             trial(t).unit(c,u).movement(m).spike_density=double(proper_PSTH(arrival_times_mov,PSTH_x_axis));%/keys.PSTH_binwidth);
%             
%             %% FR pre, peri and post
%             t_b=-0.15;
%             t_a= -0.01;
%             trial(t).unit(c,u).movement(m).FR_pre=sum(arrival_times_mov>=t_b & arrival_times_mov<=t_a)/(t_a-t_b);
%             
%             t_b= -0.01;
%             t_a= 0.05;
%             trial(t).unit(c,u).movement(m).FR_peri=sum(arrival_times_mov>=t_b & arrival_times_mov<=t_a)/(t_a-t_b);
%             
%             arrival_times_mov       =MA_out.physiology(t).spike_arrival_times{c,u}-end_t;
%             t_b= 0.05;
%             t_a= 0.15;
%             trial(t).unit(c,u).movement(m).FR_post=sum(arrival_times_mov>=t_b & arrival_times_mov<=t_a)/(t_a-t_b);
%         end
%         
%         if isempty(mov_idxs) %create dummy
%             trial(t).unit(c,u).movement=struct('endpos',NaN,'startpos',NaN,'vector',NaN,'tar_at_start',NaN,'tar_at_end',NaN,'tar_crossed',NaN,...
%                 'shape_at_start',NaN,'shape_crossed',NaN,'shape_at_end',NaN,'arrival_times',NaN,'spike_density',NaN,'FR_pre',NaN,'FR_peri',NaN,'FR_post',NaN);
%         end
    end
        population(u).FR                =nanmean([population(u).trial.FR]);
        
        if keys.FR_subtract_baseline
            population(u)=ph_FR_subtract_baseline(population(u),keys);
        end
        if keys.cal.divide_baseline_for_ANOVA
            population(u)=ph_divide_by_baseline_FR(population(u),keys);
        end
end
end


function [unit_out]=ph_FR_subtract_baseline(unit_in,keys)
unit_out=unit_in;
for t=1:numel(unit_in.trial)
    for s=1:numel(unit_in.trial(t).epoch)
        b=ismember(keys.EPOCHS(:,1),keys.EPOCHS(s,5));
        unit_out.trial(t).epoch(s).FR=unit_in.trial(t).epoch(s).FR - unit_in.trial(t).epoch(b).FR;
    end
end
end

function [unit_out]=ph_divide_by_baseline_FR(unit_in,keys)
unit_out=unit_in;
for t=1:numel(unit_in.trial)
    for s=1:numel(unit_in.trial(t).epoch)
        b=ismember(keys.EPOCHS(:,1),keys.EPOCHS(s,5));
        unit_out.trial(t).epoch(s).FR=unit_in.trial(t).epoch(s).FR/(max([unit_in.trial(t).epoch(b).FR 1]));
    end
end
end