function trials=ph_load_traces(folder,filename,trials)

allfiles=dir([folder filesep filename '*']);
allfiles=allfiles(~[allfiles.isdir]);

fields_to_take_over={'time_axis','x_eye','y_eye','x_hnd','y_hnd'};
for f=1:numel(allfiles)
    load([folder filesep allfiles(f).name])
    monkey=traces(1).monkey;
    
    %% this part only, because traces contains monkey_phys as monkey for some reason....
    us_idxs=strfind(monkey,'_');
    if ~isempty(us_idxs)
        monkey=monkey(1:  us_idxs(1)-1);
    end
    
    monkey_trials={trials.monkey};
    ismonkey=ismember(monkey_trials,monkey);
    ID_trials=[[trials.date]; [trials.block]; [trials.run]; [trials.n]];
    ID_trials=ID_trials(:,ismonkey);
    ID_traces=[[traces.date]; [traces.block]; [traces.run]; [traces.n]];
    
    ix_traces=ismember(ID_traces',ID_trials','rows');
    ix_trials=ismember(ID_trials',ID_traces','rows');

    for fn=1:numel(fields_to_take_over)
        [trials(ix_trials).(fields_to_take_over{fn})]=traces(ix_traces).(fields_to_take_over{fn});
    end
    
end

end