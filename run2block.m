function output_blocks_or_runs=run2block(base_path,given_blocks_or_runs,reverse)
folder_content=dir([base_path filesep '*.mat']);
complete_filelist=vertcat(folder_content.name);
filelist_blocks=str2num(complete_filelist(:,end-5:end-4));
filelist_runs=str2num(complete_filelist(:,end-14:end-13));

if strcmp(reverse,'reverse') || reverse==1
    [~,idx] =ismember(given_blocks_or_runs,filelist_blocks);
    output_blocks_or_runs(idx~=0)=filelist_runs(idx(idx~=0));
else
    [~,idx] =ismember(given_blocks_or_runs,filelist_runs);
    output_blocks_or_runs(idx~=0)=filelist_blocks(idx(idx~=0));
end
output_blocks_or_runs(idx==0)=0;
end