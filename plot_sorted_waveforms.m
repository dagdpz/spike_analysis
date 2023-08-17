function plot_sorted_waveforms(o,keys,title_part)
fig_title=sprintf('%s, session %s, %s',keys.monkey,keys.date,title_part);
WF_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',fig_title);
for n_unit=1:numel(o)
    n_columns_rows=ceil(numel(o)^(1/2));
    subplot(n_columns_rows,n_columns_rows,n_unit)
    block_trials=[find([true diff([o(n_unit).block(1:end-1)])~=0]) numel(o(n_unit).trial)];
    wf_per_block=[];
    for b=1:numel(block_trials)-1
        meanblockwf=nanmean(cat(1,o(n_unit).trial(block_trials(b):block_trials(b+1)).waveforms),1);
        if ~isempty(meanblockwf);  wf_per_block(b,:)=meanblockwf; end
    end
    all_spikes_wf = cat(1,o(n_unit).trial.waveforms);
    n_all_spikes_wf = size(all_spikes_wf,1);
    
    % REDUCING WFS DISPLAYED
    n_sel_spike_wf = length(1:50:n_all_spikes_wf);
    if n_sel_spike_wf>0
        set(gca,'ColorOrder',jet(n_sel_spike_wf)),hold on
        plot(all_spikes_wf(1:50:n_sel_spike_wf*50,:)');
    end
    plot(wf_per_block','-k','linewidth',2);
    unit_title={sprintf('%s SN/Si/St: %.1f/%.1f/%.1f',o(n_unit).unit_ID,o(n_unit).avg_SNR,o(n_unit).avg_single_rating,o(n_unit).avg_stability),...
        sprintf(['spk: %d ch/De: %d/%.2f b&u: %s'],n_all_spikes_wf, o(n_unit).channel,o(n_unit).electrode_depth,[o(n_unit).block_unit{:}])}; %MP add number of spikes
        %sprintf('SNR/Wam/Std: %%d/%d/%d',round(o(n_unit).quantSNR),round(o(n_unit).waveform_amplitude),round(nanmean(o(n_unit).waveform_std))),... %%KK stuff
    title(unit_title,'interpreter','none','fontsize',6)
    if  max(max(all_spikes_wf(1:50:n_sel_spike_wf*50,:))) > min(min(all_spikes_wf(1:50:n_sel_spike_wf*50,:))) % not sure what this bug is about... Lin 20160303 
    set(gca,'xtick',[],'xcolor',[1 1 1],'FontSize',6,'ytick',...
        [min(min(all_spikes_wf(1:50:n_sel_spike_wf*50,:)')) max(max(all_spikes_wf(1:50:n_sel_spike_wf*50,:)'))]);          %MP remove X-axis keep Y-axis to have scale and show max/min values on Y axis
    end
end
ph_title_and_save(WF_summary_handle,fig_title,fig_title,keys)
end
