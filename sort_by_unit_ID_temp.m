function sites = sort_by_unit_ID_temp(o_t) %sort_by_site_ID(o_t)
sites=struct('site_ID',{});
site_index=0;
for b=1:size(o_t,2)    
    if isfield(o_t(b).trial(1),'TDT_LFPx')
        AA=vertcat(o_t(b).trial.channel);
        LFP=[o_t(b).trial.TDT_LFPx];
        LFP_samples=arrayfun(@(x) size(x.TDT_LFPx,2),o_t(b).trial);
        n_chans_s=size(AA,2);
        for c=1:n_chans_s
            current_site_already_processed=ismember({sites.site_ID},AA(1,c).site_ID);            
            if any(current_site_already_processed)
                %% append to existing site
                s=find(current_site_already_processed);
            else
                %% create new site
                site_index=site_index+1;
                s=site_index;
                sites(s).site_ID          =AA(1,c).site_ID;
                sites(s).target           =AA(1,c).target;
                sites(s).perturbation_site=AA(1,c).perturbation_site;
                sites(s).grid_x           =AA(1,c).grid_x;
                sites(s).grid_y           =AA(1,c).grid_y;
                sites(s).electrode_depth  =AA(1,c).electrode_depth;
                sites(s).LFP              =[];
                sites(s).LFP_samples      =[];
                sites(s).dataset          =[];
                sites(s).perturbation     =[];
                sites(s).block            =[];
                sites(s).run              =[];
                sites(s).n                =[];
                sites(s).LFP_tStart       =[];
                sites(s).LFP_t0           =[];
            end
            sites(s).LFP_samples   =[sites(s).LFP_samples LFP_samples];
            sites(s).dataset       =[sites(s).dataset      AA(:,c).dataset];
            sites(s).perturbation  =[sites(s).perturbation AA(:,c).perturbation];
            sites(s).block         =[sites(s).block o_t(b).trial.block];
            sites(s).run           =[sites(s).run   o_t(b).trial.run];
            sites(s).n             =[sites(s).n     o_t(b).trial.n];
            sites(s).LFP_tStart    =[sites(s).LFP_tStart o_t(b).trial.TDT_LFPx_tStart];
            sites(s).LFP_t0        =[sites(s).LFP_t0 o_t(b).trial.TDT_LFPx_t0_from_rec_start];
            sites(s).LFP           =[sites(s).LFP   LFP(c,:)];
        end
    end
end

end
