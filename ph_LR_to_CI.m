function trials=ph_LR_to_CI(keys,pop,trials)
%% assigning contra and ipsi instead of left/right, dependent on the recorded target (which contains hemisphere information in the end)
% reference is left hemisphere, meaning that right becomes contra and left ipsi
% that is why hemifields, positions and hands need to be inverted for right hemisphere targets.
% for right recording sites LR->CI  (hemifield & positions...?)
% positive space and hand==2 become contra

target=pop.(keys.contra_ipsi_relative_to);
if strcmpi(target(end-1:end),'_R') || strcmpi(target(end-1:end),'_r')
    poptr=num2cell([trials.hemifield]*-1);
    [trials.hemifield]=deal(poptr{:});
    poptr=cellfun(@(x) x.*[-1,1],{trials.position},'uniformoutput',false);
    [trials.position]=deal(poptr{:});
    poptr=cellfun(@(x) x.*[-1,1],{trials.fixation},'uniformoutput',false);
    [trials.fixation]=deal(poptr{:});
    hand1=[trials.reach_hand]==2;
    hand2=[trials.reach_hand]==1;
    [trials(hand1).reach_hand]=deal(1);
    [trials(hand2).reach_hand]=deal(2);
    
    %% fixation index--> temporary solution
    if isfield(trials(1),'fix_index')
        fixindex=[trials.fix_index];
        [trials(fixindex==3).fix_index]=deal(1);
        [trials(fixindex==1).fix_index]=deal(3);
    end
end
end