function [FRrs,trs] = ph_resample_FRs(FR,t)
if numel(t)~=0
    trs=t(1):20:(t(end)+20);
    for xx=1:numel(trs)-1
        FRrs(xx)=mean(FR(trs(xx)<=t & t<trs(xx+1)));
    end
else
    trs=NaN;
    FRrs=NaN;
end

end