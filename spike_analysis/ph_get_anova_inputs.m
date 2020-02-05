function [FR,epochs,idx,u_pos,u_fix]=ph_get_anova_inputs(o,keys)
            trial_pos=vertcat(o.trial.position);
            per_epoch=vertcat(o.trial.epoch);
            %per_epoch=vertcat(population(unit).trial(tr_index).epoch);
            
            if keys.subtract_baseline_for_anovas
                per_epoch=ph_FR_subtract_baseline(per_epoch,keys);
            end
            per_epoch=per_epoch(:);
            FR=[per_epoch.FR]';
            epochs={per_epoch.state}';
            n_trials=size(trial_pos,1);
            n_epochs=numel(per_epoch)/n_trials;
            
            idx.LS =        repmat(trial_pos(:,1)<0,n_epochs,1); %% the limit here was set to 1!!???
            idx.RS =        repmat(trial_pos(:,1)>0,n_epochs,1);
            idx.NH =        repmat([o.trial.hand]'==0,n_epochs,1);
            idx.LH =        repmat([o.trial.hand]'==1,n_epochs,1);
            idx.RH =        repmat([o.trial.hand]'==2,n_epochs,1);
            idx.AH =        true(size(idx.LH));
            idx.ch =        repmat([o.trial.choice]'==1,n_epochs,1);
            idx.in =        repmat([o.trial.choice]'==0,n_epochs,1);
            idx.E1 =        repmat([o.trial.effector]'==min([o.trial.effector]),n_epochs,1);
            idx.E2 =        repmat([o.trial.effector]'==max([o.trial.effector]),n_epochs,1);
            idx.SS =        repmat([o.trial.dataset]',n_epochs,1); %% subset!            
            idx.SS(idx.SS==max(idx.SS))=1;
            idx.SS(idx.SS==min(idx.SS))=0;
            idx.SS(isnan(idx.SS))=0;
            
            temp_perturbation=repmat([o.trial.perturbation]',n_epochs,1); %% inactivation f.e.
            idx.PT =        -1*ones(size(temp_perturbation));
            idx.PT(ismember(temp_perturbation,keys.cal.perturbation_groups{1})) =0;
            idx.PT(ismember(temp_perturbation,keys.cal.perturbation_groups{2})) =1;
            
        	idx.LH_LS=   idx.LH & idx.LS;
        	idx.LH_RS=   idx.LH & idx.RS;
            idx.RH_LS=   idx.RH & idx.LS;
            idx.RH_RS=   idx.RH & idx.RS;
            
            idx.CR =        (idx.LS & idx.RH) | (idx.RS & idx.LH);
            idx.UC =        (idx.LS & idx.LH) | (idx.RS & idx.RH);
            idx.sides =     [idx.LS idx.RS]; 
            idx.LR =        [idx.LS false(size(idx.LS)) idx.RS];
            idx.hands =     [idx.AH idx.LH idx.RH];
            
            
            temp_pos=repmat(vertcat(o.trial.position),n_epochs,1);
            temp_fix=repmat(vertcat(o.trial.fixation),n_epochs,1);
            [~,~,idx.pos_x]           =unique(temp_pos(:,1));  
            [~,~,idx.pos_y]           =unique(temp_pos(:,2));  
            
            [u_pos,~,idx.pos]         =unique(temp_pos,'rows');                      
            [u_fix,~,idx.fix]         =unique(temp_fix,'rows');
            %% revert position !!
            o_pos=[u_pos(:,1)*-1 u_pos(:,2)];
            
            [~,idx.opp]=ismember(round(temp_pos*10),round(o_pos*10),'rows');
            u_pos=u_pos(:,1)+1i*u_pos(:,2);
end