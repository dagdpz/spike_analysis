function [FR,epochs,idx,u_pos,u_fix]=ph_get_anova_inputs(o,keys)
            trial_pos=vertcat(o.trial.position);
            per_epoch=vertcat(o.trial.epoch);
            
            per_epoch=per_epoch(:);
            FR=[per_epoch.FR]';
            epochs={per_epoch.state}';
            n_trials=size(trial_pos,1);
            n_epochs=numel(per_epoch)/n_trials;

           
            idx.LS =        repmat(trial_pos(:,1)<0,n_epochs,1); %% the limit here was set to 1!!???
            idx.RS =        repmat(trial_pos(:,1)>0,n_epochs,1);
            idx.NH =        repmat([o.trial.reach_hand]'==0,n_epochs,1);
            idx.LH =        repmat([o.trial.reach_hand]'==1,n_epochs,1);
            idx.RH =        repmat([o.trial.reach_hand]'==2,n_epochs,1);
            idx.AH =        true(size(idx.LH));
            idx.ch =        repmat([o.trial.choice]'==1,n_epochs,1);
            idx.in =        repmat([o.trial.choice]'==0,n_epochs,1);
            idx.E1 =        repmat([o.trial.effector]'==min([o.trial.effector]),n_epochs,1);
            idx.E2 =        repmat([o.trial.effector]'==max([o.trial.effector]),n_epochs,1);
            idx.SS =        repmat([o.trial.dataset]',n_epochs,1); %% subset!            
            idx.SS(idx.SS==max(idx.SS))=1;
            idx.SS(idx.SS==min(idx.SS))=0;
            idx.SS(isnan(idx.SS))=0;
            idx.all=   true (size(idx.SS)); 
            
            idx.suc0 =       repmat([o.trial.success]'==0,n_epochs,1);
            idx.suc1 =       repmat([o.trial.success]'==1,n_epochs,1);% easy
            idx.Diff0 =      repmat([o.trial.difficulty]'==0,n_epochs,1);
            idx.Diff1 =      repmat([o.trial.difficulty]'==1,n_epochs,1);% easy
            idx.Diff2 =      repmat([o.trial.difficulty]'==2,n_epochs,1); % difficult
            idx.StimIn2HF0 = repmat([o.trial.stimuli_in_2hemifields]'==0,n_epochs,1);
            idx.StimIn2HF1 = repmat([o.trial.stimuli_in_2hemifields]'==1,n_epochs,1);
            idx.nonDistr0 =  repmat([o.trial.n_nondistractors]'==0,n_epochs,1); % 0 target
            idx.nonDistr1 =  repmat([o.trial.n_nondistractors]'==1,n_epochs,1); % 1 targets
            idx.nonDistr2 =  repmat([o.trial.n_nondistractors]'==2,n_epochs,1); % 2 targets
            idx.Distr1 =     repmat([o.trial.n_distractors]'==1,n_epochs,1); % fixation
            idx.Distr2 =     repmat([o.trial.n_distractors]'==2,n_epochs,1); % one distractor
            idx.Distr3 =     repmat([o.trial.n_distractors]'==3,n_epochs,1); % 2 distractor
                                   
            idx.TT1HF =     idx.StimIn2HF0 & idx.nonDistr2; % double targets 1HF
            idx.TT2HF =     idx.StimIn2HF1 & idx.nonDistr2; % double targets
            idx.SglL =      idx.nonDistr1 &  idx.suc1 & idx.LS ; %% single left side
            idx.SglR =      idx.nonDistr1 &  idx.suc1 & idx.RS ; %% single right side
            
            temp_perturbation=repmat([o.trial.perturbation]',n_epochs,1); %% inactivation f.e.

            idx.PT =        -1*ones(size(temp_perturbation));
            idx.PT(isnan(temp_perturbation)) =0; %%???? this can only be useful if the respective field in the excel table is empty
            idx.PT(temp_perturbation==0) =0;
            idx.PT(temp_perturbation==1) =1;
            
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
            temp_ecc=round((temp_pos(:,1).^2 + temp_pos(:,2).^2).^(1/2));
            temp_ang=round(angle(temp_pos(:,1)+1i*temp_pos(:,2))*180/pi)*pi/180;
            [~,~,idx.pos_x]           =unique(temp_pos(:,1));  
            [~,~,idx.pos_y]           =unique(temp_pos(:,2));  
            
            [~,~,idx.ecc]           =unique(temp_ecc);  
            [~,~,idx.ang]           =unique(temp_ang);  
            [u_pos,~,idx.pos]         =unique(temp_pos,'rows');                      
            [u_fix,~,idx.fix]         =unique(temp_fix,'rows');
            %% revert position !!
            o_pos=[u_pos(:,1)*-1 u_pos(:,2)];
            
            [~,idx.opp]=ismember(round(temp_pos*10),round(o_pos*10),'rows');
            u_pos=u_pos(:,1)+1i*u_pos(:,2);
end