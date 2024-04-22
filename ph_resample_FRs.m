function [FRrs,trs] = ph_resample_FRs(FR,t)


trs=t(1):20:(t(end)+20);
for xx=1:numel(trs)-1
    FRrs(xx)=mean(FR(trs(xx)<=t & t<trs(xx+1))); 
end
% FRrs=smooth(FRrs);
% FRrs=FRrs';





%trs=trs+5;
%trs(end)=[];



%[FRrs,trs] = resample(double(FR),double(t),1/10,'linear');
% FRrs(end+1)=mean(FR(t>=trs(end)));
% trs(end+1)=trs(2)-trs(1)+trs(end);







% 
% 
% fs1 = 10;             % Original sampling frequency in Hz
% t1 = 0:1/fs1:1;       % Time vector
% x = t1 + 100;         % Define a linear sequence
% xpad = [repmat(x(1), 1, 10), x, repmat(x(end), 1, 10)];
% tpad = [-1/fs1*10 : 1/fs1: 0-1/fs1, t1, 1+1/fs1:1/fs1:1+1/fs1*10];
% ypad = resample(xpad,3,2);  % Now resample it
% t2 = (0:(length(ypad)-1))*2/(3*fs1) - 1;  % New time vector
% plot(t1,x,'*',t2,ypad,'o',(-0.5:0.01:1.5),(-0.5:0.01:1.5)+100,':')



end