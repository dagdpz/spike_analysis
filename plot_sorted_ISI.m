function plot_sorted_ISI(o,trials,keys,title_part)
fig_title=sprintf('%s, session %s, %s',keys.monkey,keys.date,title_part);
ISI_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',fig_title);

trial_IDs=[trials.block; trials.run; trials.n]';
x_bins=logspace(-3,0,30);
x_bins=horzcat(0,x_bins);
for u=1:numel(o)
    unit_trial_ID=[o(u).block; o(u).run; o(u).n]';
    trials_in_unit=trials(ismember(trial_IDs,unit_trial_ID,'rows'));
    n_columns_rows=ceil(numel(o)^(1/2));
    subplot(n_columns_rows,n_columns_rows,u)
    AT=NaN;
    for t=1:numel(o(u).trial)
        AT=[AT o(u).trial(t).arrival_times'+trials_in_unit(t).trial_onset_time];
    end
    AT=unique(AT); % due to ovrelapping end and beginning of trial, spikes can be counted twice
    all_ISI = diff(AT); %cat(2,ISI(n_unit).trial.isi);
    
    if ~isempty(all_ISI)
        hist_values = histc(all_ISI,x_bins);
        perc_first_bin = (hist_values(1)/sum(hist_values))*100;
        x_bins_bar=logspace(-3,0,31);
        bar(log10(x_bins_bar), hist_values)
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor','r','EdgeColor','w');
        x_bins_ticks = [-3:1];
        x_bins_ticks_label = [10^-3 10^-2 10^-1 0.5];
        set(gca,'ylim',[0 max([1 hist_values])],'xtick',x_bins_ticks,'box','off');%max(hist(all_ISI,x_bins)),'ycolor',[1 1 1],'xscale','log'
        set(gca,'XTickLabel',sprintf('%2.0e|',x_bins_ticks_label));
        text(-0.45,max(hist_values),[num2str(perc_first_bin,'%4.1f') '%<1ms'],'FontSize',8)
        
    end
    unit_title={sprintf('%s SN/Si/St: %d/%d/%d',o(u).unit_ID,o(u).avg_SNR,o(u).avg_single_rating,o(u).avg_stability)...
        sprintf(['spk: %d ch/De: %d/%.2f b: ' num2str(unique([o(u).block]))],numel(all_ISI), o(u).channel,o(u).electrode_depth)}; %MP add number of spikes
    title(unit_title,'interpreter','none','fontsize',10)
    
    
end
ph_title_and_save(ISI_summary_handle,fig_title,fig_title,keys)
end