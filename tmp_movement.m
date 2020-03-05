
%         trial_PSTH_time =-1:0.001:trial(t).states_onset(end);
%         PSTH_windows=keys.EPOCHS(:,1:4);
          
%         for w=1:size(PSTH_windows,1);
%             current_state=PSTH_windows{w,2};
%             trial.windows(w).name=PSTH_windows{w,1};
%             if ismember(current_state,trial_states) %% why not defining at least states?
%                 state_idx=trial_states==current_state;
%                 current_state_onset=trial_states_onset(state_idx);
%                  
%                 t_before_state=keys.EPOCHS{w,3} + current_state_onset;
%                 t_after_state=keys.EPOCHS{w,4}  + current_state_onset;
%                 PSTH_time=trial_PSTH_time>t_before_state-0.0001 & trial_PSTH_time<t_after_state+0.0001;
%                 
%                 
%                 % eye and hand
%                 start_index=find(trial.time_axis>=t_before_state,1);
%                 end_index=find(trial.time_axis>t_after_state,1);
%                 time_indexes=start_index:keys.PSTH_binwidth/eyetracker_interpolated_sampling_rate:end_index-1;
%                 
%                 % RASTER and PSTH
%                 trial.windows(w).arrival_times=trial.arrival_times(trial.arrival_times>t_before_state & trial.arrival_times<t_after_state+keys.PSTH_binwidth)-current_state_onset;
%                if sum(PSTH_time)==0
%                    trial.windows(w).spike_density=single(NaN(floor(size(PSTH_time)*keys.PSTH_binwidth/0.001))); 
%                else
%                 trial.windows(w).spike_density=single(resample(double(trial.spike_density(PSTH_time)),1,keys.PSTH_binwidth/0.001));            
%                end
%                 trial.windows(w).x_eye       =trial.x_eye(time_indexes);
%                 trial.windows(w).y_eye       =trial.y_eye(time_indexes);
%                 trial.windows(w).x_hnd       =trial.x_hnd(time_indexes);
%                 trial.windows(w).y_hnd       =trial.y_hnd(time_indexes);
%                 trial.windows(w).time_axis   =trial.time_axis(time_indexes)-current_state_onset;
%                 else %create dummy
%                 trial.windows(w)=struct('name',NaN,'arrival_times',NaN,'spike_density',NaN,'x_eye',NaN,'y_eye',NaN,'x_hnd',NaN,'y_hnd',NaN,'time_axis',NaN);
%             
%             end
%         end
%         