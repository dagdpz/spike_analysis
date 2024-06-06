function ph_exclusion_reasons(population,trials)
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
    completed=[ut.completed];
    durations=arrayfun(@(x) x.states_onset(x.states==90)- x.states_onset(x.states==2),ut);
    
    
    ER=unique(p.exclusion_reason);
    if all(ER==1)
        excluded_because(u)=1;
    elseif all(ER==1 | ER==2)
        excluded_because(u)=2;
    elseif ~any(ER==0)
        excluded_because(u)=3;
    else
        excluded_because(u)=0;
    end
    
    both_tasks_left(u)=isfield(p.criteria,'SNR_V') && isfield(p.criteria,'SNR_F') && ~isnan(p.criteria.SNR_V) && ~isnan(p.criteria.SNR_F);
    
    totake=p.exclusion_reason==0 & types==1;
    dur_F(u)=sum(durations(totake));
    totake=p.exclusion_reason==0 & types==2;
    dur_V(u)=sum(durations(totake));
    totake=p.exclusion_reason==0 & types==1 & completed;
    dur_F_completed(u)=sum(durations(totake));
    totake=p.exclusion_reason==0 & types==2 & completed;
    dur_V_completed(u)=sum(durations(totake));
end



%% now make comprehensive table ?

D={'Target','total','above FR','stable blocks','200s + rest','200s + task','200s in both','200s + rest compl','200s + task compl','200s in both compl'};
Tar=unique(targets);
for t=1:numel(Tar)
    T=Tar{t};
    TX=ismember(targets,T);

    D{t+1,1}=T;
    D{t+1,2}=sum(TX);
    D{t+1,3}=sum(TX & ~(excluded_because==1));
    D{t+1,4}=sum(TX & ~(excluded_because==1 | excluded_because==2));
    D{t+1,5}=sum(TX & excluded_because==0 & dur_F>200);
    D{t+1,6}=sum(TX & excluded_because==0 & dur_V>200);
    D{t+1,7}=sum(TX & excluded_because==0 & dur_F>200 & dur_V>200);
    D{t+1,8}=sum(TX & excluded_because==0 & dur_F_completed>200);
    D{t+1,9}=sum(TX & excluded_because==0 & dur_V_completed>200);
    D{t+1,10}=sum(TX & excluded_because==0 & dur_F_completed>200 & dur_V_completed>200);
end
save('Magnus_data','D')

end