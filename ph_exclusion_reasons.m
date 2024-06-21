function ph_exclusion_reasons(population,trials,keys)

keys.path_to_save=[keys.basepath_to_save keys.project_version];

targets={population.target};
excluded_because=zeros(size(population));
both_tasks_left=zeros(size(population));
dur_F=zeros(size(population));
dur_V=zeros(size(population));
dur_F_completed=zeros(size(population));
dur_V_completed=zeros(size(population));


for u=1:numel(population)
    
    p=population(u);
    ut=ph_get_unit_trials(p,trials);
    types=[ut.type];
    completed=[ut.completed]==1;
    durations=arrayfun(@(x) x.states_onset(x.states==90)- x.states_onset(x.states==2),ut);
    
    
    ER=p.exclusion_reason;
    if all(ER==1)
        excluded_because(u)=1;
    elseif all(ER==1 | ER==2)
        excluded_because(u)=2;
    elseif ~any(ER==0)
        excluded_because(u)=3;
    else
        excluded_because(u)=0;
    end
    
    %both_tasks_left(u)=isfield(p.criteria,'SNR_V') && isfield(p.criteria,'SNR_F') && ~isnan(p.criteria.SNR_V) && ~isnan(p.criteria.SNR_F);
    
%     totake=p.exclusion_reason==0 & types==1;
%     dur_F(u)=sum(durations(totake));
%     totake=p.exclusion_reason==0 & types==2;
%     dur_V(u)=sum(durations(totake));
    
    % entire recording duration (during trial)

    totake=types==1;
    dur_F(u)=sum(durations(totake));
    totake=types==2 ;
    dur_V(u)=sum(durations(totake));

    totake=types==1 & completed;
    dur_F_completed(u)=sum(durations(totake));
    totake=types==2 & completed;
    dur_V_completed(u)=sum(durations(totake));         
    
    
    % left after chopping
    totake=ER~=1 & types==1;
    dur_F_eN1(u)=sum(durations(totake));
    totake=ER~=1 & types==2;
    dur_V_eN1(u)=sum(durations(totake));
        
    % left after ANOVA
    totake=ER==0 & types==1;
    dur_F_e0(u)=sum(durations(totake));
    totake=ER==0 & types==2;
    dur_V_e0(u)=sum(durations(totake));


    % left after chopping completed
    totake=ER~=1 & types==1 & completed;
    dur_F_eN1_completed(u)=sum(durations(totake));
    totake=ER~=1 & types==2 & completed;
    dur_V_eN1_completed(u)=sum(durations(totake));
        
    % left after ANOVA completed
    totake=ER==0 & types==1 & completed;
    dur_F_e0_completed(u)=sum(durations(totake));
    totake=ER==0 & types==2 & completed;
    dur_V_e0_completed(u)=sum(durations(totake));
end



%% now make comprehensive table ?
D={'Target','total','above FR','stable blocks'};
Tar=unique(targets);
for t=1:numel(Tar)
    T=Tar{t};
    TX=ismember(targets,T);

    D{t+1,1}=T;
    D{t+1,2}=sum(TX);
    D{t+1,3}=sum(TX & ~(excluded_because==1));
    D{t+1,4}=sum(TX & ~(excluded_because==1 | excluded_because==2));
end
save([keys.path_to_save filesep 'Exclusion_reasons'],'D');




ratio_F_completed=dur_F_completed./dur_F;
ratio_V_completed=dur_V_completed./dur_V;

ratio_F_after_chop=dur_F_eN1./dur_F;
ratio_V_after_chop=dur_V_eN1./dur_V;

ratio_F_after_anova=dur_F_e0./dur_F;
ratio_V_after_anova=dur_V_e0./dur_V;

ratio_F_after_chop_completed=dur_F_eN1_completed./dur_F_completed;
ratio_V_after_chop_completed=dur_V_eN1_completed./dur_V_completed;

ratio_F_after_anova_completed=dur_F_e0_completed./dur_F_completed;
ratio_V_after_anova_completed=dur_V_e0_completed./dur_V_completed;

for t=1:numel(Tar)
    T=Tar{t};
    TX=ismember(targets,T);
    
    
    h=figure;
    hold on
    
    
    bins=0:0.1:1;
    hist_FC=hist(ratio_F_completed(TX),bins);
    hist_VC=hist(ratio_V_completed(TX),bins);
    
    hist_FAC=hist(ratio_F_after_chop(TX),bins);
    hist_VAC=hist(ratio_V_after_chop(TX),bins);
    
    hist_FAA=hist(ratio_F_after_anova(TX),bins);
    hist_VAA=hist(ratio_V_after_anova(TX),bins);  
    
    hist_FACC=hist(ratio_F_after_chop_completed(TX),bins);
    hist_VACC=hist(ratio_V_after_chop_completed(TX),bins);
    
    hist_FAAC=hist(ratio_F_after_anova_completed(TX),bins);
    hist_VAAC=hist(ratio_V_after_anova_completed(TX),bins);
    
    
    
    
    plot(bins,hist_FC/sum(hist_FC),'color','k');
    plot(bins,hist_VC/sum(hist_VC),'color',[0.5 0.5 0.5]);    
    
    plot(bins,hist_FAC/sum(hist_FAC),'color',[0  0.7 0.4]);
    plot(bins,hist_VAC/sum(hist_VAC),'color',[0.4 0.7 0 ]);    
    
    plot(bins,hist_FAA/sum(hist_FAA),'color',[0.2 0.6 0.7 ]);
    plot(bins,hist_VAA/sum(hist_VAA),'color',[0.7 0.6 0.2 ]);
    
    
    plot(bins,hist_FACC/sum(hist_FACC),'color',[0 0 0.6 ]);
    plot(bins,hist_VACC/sum(hist_VACC),'color',[0.6 0 0]);    
    
    plot(bins,hist_FAAC/sum(hist_FAAC),'color',[0 0 1]);
    plot(bins,hist_VAAC/sum(hist_VAAC),'color',[1 0 0]);
    
    
    legend({'Rest completed','Task completed','Rest after chop','Task after chop','Rest after ANOVA','Task after ANOVA','Rest completed after chop','Task completed  after chop','Rest completed after ANOVA','Task completed after ANOVA'});
    xlabel('fraction of total trial duration remaining');
    ylabel('fraction of units');
    title([T ', ' num2str(sum(TX)) ' units' ]);
    
    export_fig(h,[keys.path_to_save filesep T '_Chopped_duration'],'-pdf');
    
end
close all;


end