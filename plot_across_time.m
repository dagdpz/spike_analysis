
function plot_across_time(o,trials,keys,which_units,ch_start_end,whattoplot)
title_part=[which_units ' ' whattoplot ' over time, ' ch_start_end];
fig_title=sprintf('%s, session %s, %s',keys.monkey,keys.date,title_part);
FR_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',fig_title);
trial_IDs=[trials.block; trials.run; trials.n]';
n_columns_rows=ceil(numel(o)^(1/2));
firstbin=min([trials.run_onset_time]);
[lastbin,lasttrial_idx]=max([trials.run_onset_time]+[trials.trial_onset_time]);
lastbin=lastbin+max(trials(lasttrial_idx).states_onset);

tdur=arrayfun(@(x) x.states_onset(x.states==90)-x.states_onset(x.states==2),trials);

units=1:numel(o);
for u=units
    trial_blocks=[o(u).block];
    unique_blocks=unique(trial_blocks);
    unit_trial_ID=[o(u).block; o(u).run; o(u).n]';
    UT=trials(ismember(trial_IDs,unit_trial_ID,'rows'));
    subplot(n_columns_rows,n_columns_rows,u);
    hold on;
    binsize=60;
    AT=[];
    WF=[];
    for t=1:numel(o(u).trial)
        ATt=o(u).trial(t).arrival_times;
        WFt=o(u).trial(t).waveforms;
        AT=vertcat(AT,ATt(ATt>0 & ATt<UT(t).states_onset(end-1))+UT(t).trial_onset_time+UT(t).run_onset_time-firstbin);
        WF=vertcat(WF,WFt);%(ATt>0 & ATt<o(u).trial(t).states_onset(end-1),:));
    end
    
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
    
    trial_stability=[o(u).stability_rating];
    exclusion_reason=[o(u).exclusion_reason];
    switch whattoplot
        case 'FR'
            %% resample per block
            bins_rs=[];
            FR_rs=[];
            FR_rs_c=[];
            blocks_c=[];
            for b=unique_blocks
                btru=trial_blocks==b;
                FRb=o(u).FR_average(btru);
                
                
                
                btr=[UT.block]==b & ismember([UT.completed],keys.cal.completed) & ~isnan(trial_stability);
                tb=cumsum(tdur(btr));
                %% add half the duration here ?
                tb=[0 tb(1:end-1)]+ [UT(btr).run_onset_time]  - firstbin;
                tb_withITI=cumsum([0 diff([UT(btr).trial_onset_time])]) + [UT(btr).run_onset_time] - firstbin ; 
                
                FRrsb = ph_resample_FRs(FRb,tb);
                FR_rs_c=[FR_rs_c FRrsb];
                blocks_c=[blocks_c repmat(b,size(FRrsb))];
                
                
                btr=[UT.block]==b;
                tb=cumsum(tdur(btr));
                %% add half the duration here ?
                tb=[0 tb(1:end-1)]+ [UT(btr).run_onset_time]  - firstbin;
                tb_withITI=cumsum([0 diff([UT(btr).trial_onset_time])]) + [UT(btr).run_onset_time] - firstbin ; 
                
                binrsb= ph_resample_FRs(tb_withITI,tb);
                FRrsb = ph_resample_FRs(FRb,tb);
                
                bins_rs=[bins_rs binrsb];
                FR_rs=[FR_rs FRrsb];
            end
            bins=bins_rs;
            toplot=FR_rs;
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
    plot(bins,toplot,'k','linewidth',0.1);
    y_lim=ylim(gca);
    for b=unique_blocks
        
        tr_ok=trial_blocks==b;
        if sum(tr_ok)<2; continue; end;            % it can happen that an entire block is not accepted if FR changed drastically
        if all([UT(tr_ok).type]==1)
            style=':';
        else
            style='-';
        end
        if any(tr_ok & ~isnan(trial_stability))
            tr_bad = tr_ok & isnan(trial_stability);
            tr_ok = tr_ok & ~isnan(trial_stability);
            block_mean=double(nanmean(FR_rs_c(blocks_c==b)));
        else
            block_mean=0;
            tr_bad = tr_ok;
        end
        if any(tr_bad)
            bad_starts=find(diff([false tr_bad])==1);
            bad_ends  =find(diff([tr_bad false])==-1);
            % first bit
            
            start_block=UT(bad_starts(1)).run_onset_time-firstbin+UT(bad_starts(1)).trial_onset_time;
            end_block=start_block+UT(bad_ends(1)).trial_onset_time-UT(bad_starts(1)).trial_onset_time;
            plot([start_block end_block],[0 0],'color',[0.5 0.5 0.5],'linestyle',style,'linewidth',1.5)
            
            % second part
            if numel(bad_starts) == 2                
                start_block=UT(bad_starts(2)).run_onset_time-firstbin+UT(bad_starts(2)).trial_onset_time;
                end_block=start_block+UT(bad_ends(2)).trial_onset_time-UT(bad_starts(2)).trial_onset_time;
                plot([start_block end_block],[0 0],'color',[0.5 0.5 0.5],'linestyle',style,'linewidth',1.5)
            elseif numel(bad_starts) >2
                disp('3 invalid intervals for this block ??');
            end
        end
        
        %FR_std=double(nanstd(FR_smoothed(tr_idx)));
        start_block=UT(find(tr_ok,1,'first')).run_onset_time-firstbin+UT(find(tr_ok,1,'first')).trial_onset_time;
        end_block=start_block+UT(find(tr_ok,1,'last')).trial_onset_time-UT(find(tr_ok,1,'first')).trial_onset_time;
        fanoish_factor=trial_stability(tr_ok);fanoish_factor=fanoish_factor(1);
        exclusion_code=unique(exclusion_reason(tr_ok));
        if numel(exclusion_code)>1
                disp('several reasons to exclude block ??');
        end
        switch exclusion_code
            case 0 %% not excluded
            col='g';   
            case 1 %% low FR
            col=[0.5 0.5 0.5];
            case 2 %% unstable in general
            col='r';
            case 3 %% low spike count
            col='m';
            case 4 %% most different FR 
            col='b';
             case 5 %% lower stability
            col='c';
        end
        
        
        
        plot([start_block end_block],[block_mean block_mean],'color',col,'linestyle',style,'linewidth',1.5)
        plot([start_block start_block],[0 block_mean],'color',col,'linestyle',style,'linewidth',1.5)
        plot([end_block end_block],[0 block_mean],'color',col,'linestyle',style,'linewidth',1.5)
        if strcmp(whattoplot,'FR')
            text(double(start_block+(end_block-start_block)/2), diff(y_lim)/4,sprintf('%0.1f',fanoish_factor),'fontsize',8,'HorizontalAlignment', 'Center')
        end
    end
    unit_title={sprintf('%s %.1f Hz ch/De: %d/%.2f ',o(u).unit_ID,nanmean(o(u).FR_average),o(u).channel,o(u).electrode_depth),...
        sprintf('b&u: %s',[o(u).block_unit{:}])}; %MP add number of spikes
    title(unit_title,'interpreter','none','fontsize',6);
    switch whattoplot
        case 'FR'
            set(gca,'xlim',[0,lastbin-firstbin]);
    end
end
ph_title_and_save(FR_summary_handle,fig_title,fig_title,keys)
end