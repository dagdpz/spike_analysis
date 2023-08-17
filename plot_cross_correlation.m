function plot_cross_correlation(population,keys,title_part)
chs=unique([population.channel]);
n_rows=ceil(sqrt(numel(chs)));
n_columns=ceil(sqrt(numel(chs)));
FR_summary_handle=figure;
fig_title=sprintf('%s, session %s, %s',keys.monkey,keys.date,title_part);
for c=1:numel(chs)
    CRL=[];
    ch=chs(c);
    subplot(n_rows,n_columns,c)
    pop=population([population.channel]==ch);
    unit_ticks=zeros(size(pop));
    for u1=1:numel(pop)
        u2=0;
        unit_ticks(u1)=str2double(pop(u1).unit_ID(end-1:end));
        while u2<u1
            u2=u2+1;
            u1_runs=unique([pop(u1).run]);
            u2_runs=unique([pop(u2).run]);
            u1_blocks=unique([pop(u1).block]);
            u2_blocks=unique([pop(u2).block]);
            r=intersect(u1_runs,u2_runs);
            b=intersect(u1_blocks,u2_blocks);
            if ~isempty(r)
                r1=ismember([pop(u1).run],r) & ismember([pop(u1).block],b);
                r2=ismember([pop(u2).run],r) & ismember([pop(u2).block],b);
                coeff=corr([[pop(u1).FR_average(r1)]',[pop(u2).FR_average(r2)]']);
                CRL(u1,u2)=coeff(1,2);
            else
                CRL(u1,u2)=0;
            end
        end
        
    end
    CRL=CRL+CRL';
    CRL(CRL>1)=1;
    
    X=1:numel(unit_ticks);
    imagesc(X,X,CRL,[0 1])
    set(gca,'ytick',X);
    set(gca,'yticklabel',unit_ticks);
    set(gca,'xtick',X);
    set(gca,'xticklabel',unit_ticks);
    xlabel('unit ID');
    ylabel('unit ID');
    title(['channel ' num2str(ch)]);
    colormap('jet')
    cb = colorbar;
    set(get(cb,'title'),'string', 'corr trial FR', 'fontsize',8);
end
ph_title_and_save(FR_summary_handle,fig_title,fig_title,keys)
close(gcf);

end
