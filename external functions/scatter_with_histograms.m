function scatter_with_histograms(SC)

%% missing: scale
SC.markersize=10;
markerratio=1.5;
x_values_scat=[SC.x_values{:}];
y_values_scat=[SC.y_values{:}];
hold on

%% scatter
for idx=1:numel(x_values_scat)
    if ismember(SC.markers{idx},{'>','<','^','v'})
        markersize=SC.markersize*markerratio;
    else
        markersize=SC.markersize;
        
    end
    if ~isempty(SC.filled{idx})
        %scatter(x_values_scat{idx}, y_values_scat{idx}, markersize, SC.colors{idx},SC.markers{idx},SC.filled{idx});
        scatter(x_values_scat{idx}, y_values_scat{idx}, markersize, SC.colors{idx},SC.markers{idx},'filled','MarkerEdgeColor', [0 0 0] );
    else
        %scatter(x_values_scat{idx}, y_values_scat{idx}, markersize, SC.colors{idx},SC.markers{idx});
        sc_temp=scatter(x_values_scat{idx}, y_values_scat{idx}, markersize, SC.colors{idx},SC.markers{idx});
        set(sc_temp,'markerfacecolor',[1 1 1]); % not sure if this works in all cases with different edge colors (empty)
    end
end

axes_limits=[min([get(gca,'xlim') get(gca,'ylim')]) max([get(gca,'xlim') get(gca,'ylim')])];
set(gca,'xlim',axes_limits,'ylim',axes_limits);
axis square;
plot(axes_limits,axes_limits,'k'); %diagonal

%% histograms
for idx=1:size(SC.x_values,2)
    temp=[SC.x_values{:,idx}];
    SC.x{idx}=[temp{:}];
end
temp=[SC.x_values{:}];
SC.x{idx+1}=[temp{:}];
% exception (actually default)
if size(SC.y_values,1)==1
    for idx=1:size(SC.y_values,2)
        temp=[SC.y_values{:,idx}];
        SC.y{idx}=[temp{:}];
    end
else
    
    for idx=1:size(SC.y_values,1)
        temp=[SC.y_values{idx,:}];
        SC.y{idx}=[temp{:}];
    end
    
end
temp=[SC.y_values{:}];
SC.y{idx+1}=[temp{:}];
for idx=1:numel(SC.d_values)
    if ~isempty(SC.d_values{idx}{1})
        [corr_Rs{idx},corr_Ps{idx}]=corr(SC.d_values{idx}{1}(:),SC.d_values{idx}{2}(:));
    else
        corr_Rs{idx}=0;corr_Ps{idx}=1;
    end
    SC.d{idx}=SC.d_values{idx}{1} - SC.d_values{idx}{2};
    
    %% paried tests
    if ~isempty(SC.d_values{idx}{1})
        [~,ppaired(idx)]=ttest2([SC.d_values{idx}{1}],[SC.d_values{idx}{2}]);
        ppairedabs(idx) =signrank(abs([SC.d_values{idx}{1}]),abs([SC.d_values{idx}{2}]));
    else
        
        ppaired(idx)=1;
        ppairedabs(idx)=1;
    end
end
temp=vertcat(SC.d_values{:});
[corr_Rs{idx+1},corr_Ps{idx+1}]=corr([temp{:,1}]',[temp{:,2}]');
SC.d{idx+1}=[temp{:,1}]-[temp{:,2}];
[~,ppaired(idx+1)]=ttest2([temp{:,1}],[temp{:,2}]);
ppairedabs(idx+1)=signrank(abs([temp{:,1}]),abs([temp{:,2}]));


N_bins=20;
minmax=[min(axes_limits)-max(axes_limits) max(axes_limits)-min(axes_limits)];
bins_xy=axes_limits(1):diff(axes_limits)/N_bins:axes_limits(2);
bins_d=minmax(1):diff(minmax)/N_bins/2:minmax(2);

max_N_for_scale=ceil(max([hist([SC.x{end}],bins_xy) hist([SC.y{end}],bins_xy) hist([SC.d{end}],bins_d)])/10)*10;
hist_scale=diff(axes_limits)/max_N_for_scale/2;

histograms={'x','y','d'};

for counter=1:numel(histograms)
    SC.hist_colors{counter}{end+1}='k';
    fname=histograms{counter};
    scale=1;
    bins=bins_xy;
    flip_hist=1;
    switch fname
        case 'x'
            offset= [0 -diff(axes_limits)/5 + axes_limits(1)];
            phi= 0;
            flip_hist=-1;
        case 'y'
            offset= [-diff(axes_limits)/5 + axes_limits(1) 0];
            phi= pi/2;
        case 'd'
            OS=axes_limits(1)+diff(axes_limits)/2;
            offset= [OS OS];
            phi= -pi/4;
            scale=abs(sin(phi));
            bins=bins_d;
    end
    
    wid = (bins(2)-bins(1))/2.5*scale;%1/numel(SC.y_values);
    
    R = [cos(phi),-sin(phi);sin(phi),cos(phi)];
    to_plot=SC.(fname);
    
    kk=zeros(numel(bins),4);
    for idx=1:numel(to_plot)
        
        means.(fname){idx}=nanmean([to_plot{idx}]);
        [~,ps.(fname){idx}]=ttest([to_plot{idx}]);
        Ns.(fname){idx}=numel(to_plot{idx});
        if idx<numel(to_plot) % because the last one is all
            h=hist([to_plot{idx}],bins);
            b=bins;
            h=h*flip_hist;
            b=b*scale;
            kk = [zeros(numel(h),1) h(:) h(:) zeros(numel(h),1)] +kk;
            v=kk(:,2)~=0;
            kk2 = reshape(kk(v,:)',[1,numel(kk(v,:))]);
            bb = sort([b(v)-wid b(v)-wid b(v)+wid b(v)+wid]);%+wid*(idx-1);
            
            coords = [bb(:),kk2(:)*hist_scale];
            coords = R * coords' ;
            
            p=patch(coords(1,:)+offset(1),coords(2,:)+offset(2),SC.hist_colors{counter}{idx},'clipping','off');
            uistack(p,'bottom');
        end
        xy_line(:,1)=offset'+R*[means.(fname){idx} 0]';
        xy_line(:,2)=offset'+R*[means.(fname){idx} diff(axes_limits)/2*flip_hist]';
        %
        line(xy_line(1,:),xy_line(2,:),'linestyle',':','color',SC.hist_colors{counter}{idx},'clipping','off')
    end
end

y_lim_out=[axes_limits(1)-diff(axes_limits)/5-hist_scale*max_N_for_scale axes_limits(2)];
x_lim_out=[axes_limits(1)-diff(axes_limits)/5-hist_scale*max_N_for_scale axes_limits(2)];

%outerpos=get(gca,'outerposition')
innerpos=get(gca,'Position');
outerpos=get(gca,'outerposition');
% pos1=(diff(x_lim_out)-diff(axes_limits))/diff(x_lim_out);
% pos2=(diff(y_lim_out)-diff(axes_limits))/diff(y_lim_out);
pos1=innerpos(1)+((axes_limits(1)-x_lim_out(1))/diff(x_lim_out))*innerpos(3);
pos2=innerpos(2)+((axes_limits(1)-y_lim_out(1))/diff(y_lim_out))*innerpos(4);
pos3=(diff(axes_limits))/diff(x_lim_out)*innerpos(3);
pos4=(diff(axes_limits))/diff(y_lim_out)*innerpos(4);


%% add text
fonstsize=6;
space_available=[axes_limits(1)-x_lim_out(1) axes_limits(1)-y_lim_out(1)];

rows_start=y_lim_out(1)+space_available(2)*0.8;
col_start=x_lim_out(1)+space_available(1)/7;
rows_sep=-1*space_available(2)/14;
rows_sep2=rows_start;
col_sep=space_available(1)/7;

for counter=1:numel(histograms)
    fname=histograms{counter};
    to_plot=means.(fname);
    text(col_start,rows_sep2            ,[upper(fname) ': N= '],'fontsize',fonstsize)
    text(col_start,rows_sep +rows_sep2  ,[upper(fname) ': u= '],'fontsize',fonstsize)
    text(col_start,2*rows_sep +rows_sep2,[upper(fname) ': p= '],'fontsize',fonstsize)
    for idx=1:numel(to_plot)
        text(col_start+col_sep*idx,rows_sep2            ,sprintf('%d',Ns.(fname){idx}),'color', SC.hist_colors{counter}{idx},'fontsize',fonstsize) ;
        text(col_start+col_sep*idx,rows_sep +rows_sep2  ,sprintf('%0.1f', means.(fname){idx}),'color', SC.hist_colors{counter}{idx},'fontsize',fonstsize) ;
        text(col_start+col_sep*idx,2*rows_sep +rows_sep2,sprintf('%0.3f',ps.(fname){idx}),'color', SC.hist_colors{counter}{idx},'fontsize',fonstsize) ;
    end
    rows_sep2=rows_sep2+3*rows_sep;
end

% paired statistics

text(col_start,rows_sep2          ,'Ttest2: ','fontsize',fonstsize)
text(col_start,rows_sep2+ rows_sep,'Nonpar: ','fontsize',fonstsize)
for idx=1:numel(ppaired)    
    text(col_start+col_sep*idx,rows_sep2            ,sprintf('%0.3f', ppaired(idx)),'color', SC.hist_colors{3}{idx},'fontsize',fonstsize) ;
    text(col_start+col_sep*idx,rows_sep2 + rows_sep ,sprintf('%0.3f', ppairedabs(idx)),'color', SC.hist_colors{3}{idx},'fontsize',fonstsize) ;
end

%% add correlations text
rows_sep=-1*space_available(2)/14;
col_sep=space_available(1)/7;
rows_start=y_lim_out(2)-space_available(2)/20;
col_start=x_lim_out(1)+space_available(1)/7;

text(col_start,rows_start            ,'R= ','fontsize',fonstsize)
text(col_start,rows_start +rows_sep  ,'p= ','fontsize',fonstsize)
for idx=1:numel(corr_Ps)
    text(col_start+col_sep*idx,rows_start            ,sprintf('%0.2f',corr_Rs{idx}),'color', SC.hist_colors{3}{idx},'fontsize',fonstsize) ;
    text(col_start+col_sep*idx,rows_start +rows_sep  ,sprintf('%0.3f', corr_Ps{idx}),'color', SC.hist_colors{3}{idx},'fontsize',fonstsize) ;
end

text(col_start,rows_start +2*rows_sep  ,['N scale= ' num2str(max_N_for_scale)],'fontsize',fonstsize)


set(gca,'Position',[pos1 pos2 pos3 pos4])

%coords = [coords;b(end),0;b(1),0];
end