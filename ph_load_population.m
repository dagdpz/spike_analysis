function pop=ph_load_population(folder,filename)

allfiles=dir([folder filesep filename '*']);
allfiles=allfiles(~[allfiles.isdir]);

fields_to_remove={'TDT_CAP1_SR','TDT_CAP1_t0_from_rec_start','TDT_CAP1_tStart','TDT_ECG1_SR','TDT_ECG1_t0_from_rec_start',...
    'TDT_ECG1_tStart','TDT_ECG4_SR','TDT_ECG4_t0_from_rec_start','TDT_ECG4_tStart','TDT_LFPx_SR','TDT_LFPx_t0_from_rec_start',...
    'TDT_LFPx_tStart','TDT_POX1_SR','TDT_POX1_t0_from_rec_start','TDT_POX1_tStart'};
n_population=0;
for f=1:numel(allfiles)
    fname=allfiles(f).name;
    und_idx=strfind(fname,'_');
    load([folder filesep fname])
    population=eval(fname(1:und_idx(1)-1));
    %temporary for noise correlation
    [population.session]=deal(f);
    [population.date]=deal(str2double(fname(und_idx(end)+1:end-4)));
    [population.monkey]=deal(fname(und_idx(1)+1:und_idx(2)-1));
    fn_population=fieldnames(population);
    population=rmfield(population,fn_population(ismember(fn_population,fields_to_remove)));
    pop(n_population+1:n_population+numel(population))=population;
    n_population=numel(pop);
end

%temporary stuff to plot population for pulvinar reach
fields_to_remove={};
% for p=1:numel(pop)
%     if isfield(pop(p).trial,'stability_rating')
%         [pop(p).trial]=rmfield(pop(p).trial,'stability_rating');
%     end
% end

if isempty(allfiles)
    pop=struct();
end
end