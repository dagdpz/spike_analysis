function export_fig(varargin)
matlabversion=datevec(version('-date'));
if matlabversion(1)>=2014
export_fig_2014(varargin{:});    
else
export_fig_2011(varargin{:});
end
end
