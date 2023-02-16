function ph_copy_single_units(keys)

TT=keys.tuning_table;
TT(cellfun(@(x) isempty(x) & ~islogical(x),TT))={''}; %% is this needed?
idx_ID=DAG_find_column_index(TT,'unit_ID');
ph_copy_file_list_to_dir(TT(:,idx_ID),keys);
end

function ph_copy_file_list_to_dir(example_list,keys)
additional_pattern = '';
remove_original = 0;
verbose = 0;
    
project=keys.project;
version=keys.project_version;
fromdir=['Y:\Projects\' project '\ephys\' version '\single_cells'];
todir=['Y:\Projects\' project '\ephys\' version '\single_cells' '\' keys.CP.foldername];
% if ~exist(todir,'dir')
%  mkdir(todir)
% end
save([todir '\list'],'example_list')
k=1;
for ex=1:numel(example_list)
    name=example_list{ex};
    if strfind(name,'*'),
        pattern = name;
    else
        pattern = ['*' name '*' additional_pattern '*'];
    end
    [d] = dir([fromdir filesep pattern]);
    n = length(d);
    if n>0,
        list(k:k+n-1) = {d.name};
        k = k + n;
    end
    if verbose,
        fprintf('Processing list entry %d, found %d files corresponding to pattern %s\n',ex,n,pattern);
    end
    
end
k = k - 1;

success=1;
for f=1:k,
    if remove_original,
        [success,message] = movefile([fromdir filesep list{f}],todir);
    else
        [success,message] = copyfile([fromdir filesep list{f}],todir);
    end
    if ~success,
        disp(message);
        return;
    end
end

if success,
    fprintf('Copied %d files from %s to %s\n',k,fromdir,todir);
end
end