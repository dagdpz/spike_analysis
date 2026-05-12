function wc_file_name = ph_figure_out_waveclus_file_by_channel_and_blocks(chNum, blkNum, concatenation_info)
% ph_figure_out_waveclus_file_by_channel_and_blocks
% This function figures out the corresponding waveclus file that contains
% data for a unit recorded on the specified channel and blocks. It's
% helpful for figuring out the spike detection thresholds.
%
% USAGE:
%	wc_file_name = ph_figure_out_waveclus_file_by_channel_and_blocks(chNum, blkNum, session_info)
%
% INPUTS:
%       chNum - channel number taken from population.channel in
%       spike_analysis output
%
%       blkNum - unique block list taken from [population.block]
%
%       session_info - the info on the current session that contains data
%       paths (path to waveclus files is required here)
%       
% OUTPUTS:
%		wc_file_name - string containing the full path to the corresponding
%		waveclus file
%
% Author(s):	L.N. Vasileva, DAG, DPZ
% URL:		http://www.dpz.eu/dag
%
% Change log:
% 2023-11-06:	Created function
% 
%%%%%%%%%%%%%%%%%%%%%%%%%[DAG mfile header version 1]%%%%%%%%%%%%%%%%%%%%%%%%%

chNumStr  = ['00' num2str(chNum)];
chNumStr  = chNumStr(end-2:end);

%concatWCfile = dir([session_info.Input_WC 'concatenation_info.mat']);
concatInfo = load([concatenation_info 'concatenation_info.mat']);

matchBlocks = cellfun(@(x) any(ismember(x, blkNum)), concatInfo.whattofindwhere{chNum}); % check that block lists overlap at least partially
    
fileNums = cellfun(@(x) x{chNum}, concatInfo.wheretofindwhat(concatInfo.whattofindwhere{chNum}{matchBlocks}));

wcFileNum = unique(fileNums);

wc_file_name = [concatenation_info 'dataspikes_ch' chNumStr '_' num2str(wcFileNum) '.mat'];

end
