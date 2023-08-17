function plot_across_time(o,trials,keys,which_units,ch_start_end,whattoplot)
title_part=[which_units ' ' whattoplot ' over time, ' ch_start_end];
fig_title=sprintf('%s, session %s, %s',keys.monkey,keys.date,title_part);
FR_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',fig_title);

trial_IDs=[trials.block; trials.run; trials.n]';

n_columns_rows=ceil(numel(o)^(1/2));
%alltrialz=[o.trial];
firstbin=min([trials.run_onset_time]);
[lastbin,lasttrial_idx]=max([trials.run_onset_time]+[trials.trial_onset_time]);
lastbin=lastbin+max(trials(lasttrial_idx).states_onset);

units=1:numel(o);
for u=units
    
    unit_trial_ID=[o(u).block; o(u).run; o(u).n]';
    trials_in_unit=trials(ismember(trial_IDs,unit_trial_ID,'rows'));
    
    subplot(n_columns_rows,n_columns_rows,u);
    hold on;
    binsize=60;
    AT=[];
    WF=[];
    for t=1:numel(o(u).trial)
        ATt=o(u).trial(t).arrival_times;
        WFt=o(u).trial(t).waveforms;
        AT=vertcat(AT,ATt(ATt>0 & ATt<trials_in_unit(t).states_onset(end-1))+trials_in_unit(t).trial_onset_time+trials_in_unit(t).run_onset_time-firstbin);
        WF=vertcat(WF,WFt);%(ATt>0 & ATt<o(u).trial(t).states_onset(end-1),:));
    end
    %bins=(o(u).trial(1).trial_onset_time):binsize:(o(u).trial(end).run_onset_time-o(u).trial(1).run_onset_time+o(u).trial(end).trial_onset_time+max(o(u).trial(end).states_onset));
    bins=0:binsize:(lastbin-firstbin);
    
    if ismember(whattoplot,{'SNR','amp','noise'})
        snr=NaN(size(bins));
        amp=NaN(size(bins));
        noi=NaN(size(bins));
        
        for b=1:numel(bins)
            idx=AT>bins(b)-binsize/2 & AT<bins(b)+binsize/2;
            meanwf=mean(WF(idx,:),1);
            noi(b)=mean(std(WF(idx,:),1));
            amp(b)=abs(max(meanwf)-min(meanwf));
            snr(b)=amp(b)/noi(b);
        end
        
    end
    switch whattoplot
        case 'FR'
            toplot=hist(AT,bins')/binsize;
            toplot_per_trial=[o(u).FR_average];
            toplot_per_trial(isnan(toplot_per_trial))=0;
        case 'SNR'
            toplot=snr;
            toplot_per_trial=[o(u).SNR_rating];
        case 'amp'
            toplot=amp;
            toplot_per_trial=zeros(numel(o(u).trial),1);
        case 'noise'
            toplot=noi;
            toplot_per_trial=zeros(numel(o(u).trial),1);
    end
    plot(bins,toplot);
    y_lim=ylim(gca);
    trial_blocks=[o(u).block];
    trial_stability=[o(u).stability_rating];
    unique_blocks=unique(trial_blocks);
    for b=unique_blocks
        tr_idx=trial_blocks==b & ~isnan(trial_stability);
        if sum(tr_idx)<2; continue; end;            % it can happen that an entire block is not accepted if FR changed drastically
        %FR_std=double(nanstd(FR_smoothed(tr_idx)));
        block_mean=double(nanmean(toplot_per_trial(tr_idx)));
        start_block=trials_in_unit(find(tr_idx,1,'first')).run_onset_time-firstbin+trials_in_unit(find(tr_idx,1,'first')).trial_onset_time;
        end_block=start_block+trials_in_unit(find(tr_idx,1,'last')).trial_onset_time-trials_in_unit(find(tr_idx,1,'first')).trial_onset_time;
        fanoish_factor=trial_stability(tr_idx);fanoish_factor=fanoish_factor(1);
        if fanoish_factor > 5 %% replace with keys
            col='g';
        elseif fanoish_factor> 2.5
            col='b';
        else
            col='r';
        end
        plot([start_block end_block],[block_mean block_mean],col,'linewidth',2)
        plot([start_block start_block],[0 block_mean],col,'linewidth',2)
        plot([end_block end_block],[0 block_mean],col,'linewidth',2)
        if strcmp(whattoplot,'FR')
            text(double(start_block+(end_block-start_block)/2), diff(y_lim)/2,sprintf('%0.1f',fanoish_factor),'HorizontalAlignment', 'Center')
        end
        %plot()
    end
    unit_title={sprintf('%s ',o(u).unit_ID),...
        sprintf(['ch/De: %d/%.2f b&u: %s' ],o(u).channel,o(u).electrode_depth,[o(u).block_unit{:}])}; %MP add number of spikes
    title(unit_title,'interpreter','none');
end
ph_title_and_save(FR_summary_handle,fig_title,fig_title,keys)
end
