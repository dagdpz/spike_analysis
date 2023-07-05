function pop=ph_load_population(folder,filename)

allfiles=dir([folder filesep filename '*']);
allfiles=allfiles(~[allfiles.isdir]);
n_population=0;
for f=1:numel(allfiles)
    load([folder filesep allfiles(f).name])
    %temporary for noise correlation
    [population.session]=deal(f);
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