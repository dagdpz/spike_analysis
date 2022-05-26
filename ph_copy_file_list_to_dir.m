function ph_copy_file_list_to_dir(example_list,project,versions, folder, additional_pattern,remove_original,verbose)
% example_list=ph_initiate_example_cells_readout(project,versions);
%ph_copy_file_list_to_dir(example_list,'Pulv_eye_gaze_position',{'20190309_Cur20151217'}, 'movement');
%ph_initiate_example_cells_readout(project,versions)

%IG_COPY_FILE_LIST_TO_DIR		- copy or move files in filelist from fromdir to todir
% e.g. ig_copy_file_list_to_dir('Y:\Projects\Pulv_oculomotor\ephys\20180222\single_cell_examples','C:\Temp','C:\Temp\cells_perisac_ehn_42.txt',0,'memory',0,1);


if nargin < 5,
    additional_pattern = '';
end

if nargin < 6,
    remove_original = 0;
end

if nargin < 7,
    verbose = 0;
end

for v=1:numel(versions)
    version=versions{v};
    fromdir=['Y:\Projects\' project '\ephys\' version '\single_cell_examples'];
    %[success,message] = copyfile([fromdir filesep list{f}],todir);
    
    
    k = 1;
    todir=['Y:\Projects\' project '\ephys\' version '\included_units\' folder];
    
    for ex=1:numel(example_list)
        name=example_list{ex,1};
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
    
    
    if ~isdir(todir),
        [success,message] = mkdir(todir);
        if ~success,
            disp(message);
            return;
        end
    end
    
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



end