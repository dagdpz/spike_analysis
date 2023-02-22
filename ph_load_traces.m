function pop=ph_load_traces(folder,filename,pop)

allfiles=dir([folder filesep filename '*']);
allfiles=allfiles(~[allfiles.isdir]);
n_block=0;
for f=1:numel(allfiles)
    load([folder filesep allfiles(f).name])
    traces(n_block+1:n_block+numel(traces_per_block))=traces_per_block;
    n_block=numel(traces);
end

sessions=arrayfun(@(x) x.trial(1).date,traces);
blocks=arrayfun(@(x) x.trial(1).block,traces);

for u=1:numel(pop)
    for t=1:numel(pop(u).trial)
        
        %find block
        b=sessions==pop(u).trial(t).date & blocks==pop(u).trial(t).block;
        %find trial
        T=pop(u).trial(t).n;
        
        pop(u).trial(t).x_eye=traces(b).trial(T).x_eye;
        pop(u).trial(t).y_eye=traces(b).trial(T).y_eye;
        pop(u).trial(t).x_hnd=traces(b).trial(T).x_hnd;
        pop(u).trial(t).y_hnd=traces(b).trial(T).y_hnd;
        pop(u).trial(t).time_axis=traces(b).trial(T).time_axis;
    end
end

end