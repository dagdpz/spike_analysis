function varargout=shadedErrorBar(x,y,errBar,lineProps,transparent)
matlabversion=datevec(version('-date'));
if matlabversion(1)>=2014
    out=shadedErrorBar_2014(x,y,errBar,lineProps,transparent);
else
    out=shadedErrorBar_2011(x,y,errBar,lineProps,transparent);
end
if ~iscell(out)
varargout={out};
else
varargout=out;
end
end
