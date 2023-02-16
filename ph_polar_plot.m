function ph_polar_plot(means,sterrs,angles,line_spec,limit)
if numel(means)>=1
    means(end+1)=means(1);
    sterrs(end+1)=sterrs(1);
    angles(end+1)=angles(1);
    separations=find(diff(angles)>2*pi/6);
end
msterr=means-sterrs;
msterr(sign(msterr)==-1)=0;
psterr=means+sterrs;
line_spec.visible='off';
polar(0,limit,line_spec);

hold on
rescale_factor=round(get(gca,'ylim')./limit*100)/100;

sterrangles=[angles,fliplr(angles)];
sterrdata=[msterr,fliplr(psterr)];

h=polar(sterrangles,sterrdata*rescale_factor(2),line_spec);
patch( get(h,'XData'), get(h,'YData'), line_spec.color*0.5,'linestyle','none');

line_spec.visible='on';
polar(angles,means*rescale_factor(2),line_spec);

%rescaling text
th = findall(gca,'Type','text');
for i = 1:length(th),
    set(th(i),'FontSize',6)
    if any(strfind(get(th(i),'String'),'  '))
        set(th(i),'String',num2str(str2double(get(th(i),'String'))/rescale_factor(2)));
    end
end

end