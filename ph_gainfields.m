function population=ph_gainfields(population,modified_keys)
echo off %% ??
warning('off','MATLAB:catenate:DimensionMismatch');
do_plot=1;
type=3;
markersize=20;
n_iterations_bootstrap=1000;
significance_measure='paired';%'paired'; %'unpaired';
method='FR_vector'; %'FR_vector'; 'Von_mises_fit'; %

for fn=fieldnames(modified_keys)'
    keys.(fn{:})=modified_keys.(fn{:});
end
%% tuning table preparation and grouping
[tuning_per_unit_table]                 = ph_load_extended_tuning_table(keys);
[tuning_per_unit_table, Sel_for_title]  = ph_reduce_tuning_table(tuning_per_unit_table,keys);
if keys.GF.FR_subtract_baseline
    Sel_for_title =[Sel_for_title,{'base';'=';keys.GF.epoch_BL;', '}];
end
idx_group_parameter=DAG_find_column_index(tuning_per_unit_table,keys.GF.group_parameter);
group_values=tuning_per_unit_table(:,idx_group_parameter);
group_values=cellfun(@num2str, group_values, 'UniformOutput', false);
cell_in_any_group=[false; ~ismember(group_values(2:end),keys.GF.group_excluded)];
unique_group_values=unique(group_values(cell_in_any_group));
if isempty(unique_group_values)
    disp('no relevant groups found');
    return;
end
complete_unit_list={population.unit_ID}';
idx_unitID=DAG_find_column_index(tuning_per_unit_table,'unit_ID');

%%  here we go
N_units=numel(population);
angles_to_plot=-pi:0.2:pi;
epochs_names=keys.EPOCHS_PER_TYPE{type}(:,1);
epochs_with_baselines=keys.ANOVAS_PER_TYPE(type).epoch;
epochs_names_considered=epochs_names(ismember(epochs_names,epochs_with_baselines(:,2)));

n_epochs=numel(epochs_names_considered);

pct1 = 100*0.05/2; % alpha 0.05
pct2 = 100-pct1;

F = @(a,x) a(1)+exp(cos(x-a(2))*a(3))*a(3)/pi;
opts = optimset('lsqcurvefit');
opts.Display='off';
for u=1:N_units % %48%%90:N_units
    plot_title=population(u).unit_ID;
    filename=population(u).unit_ID;
    tr=[population(u).trial.completed];
    [pop]=ph_arrange_positions_and_plots(population(u).trial(tr),keys);
    pop.trial=pop.trial([pop.trial.accepted]);
    pop=ph_epochs(pop,keys);
    Positions=vertcat(pop.trial.position);
    Positions=Positions(:,1)+1i*Positions(:,2);
    [unique_positions,~,temp_pos_idx]=unique(Positions);
    unique_pos_indexes=unique(temp_pos_idx);
    N_per_pos=hist(temp_pos_idx,unique_pos_indexes);
    angles=angle(Positions);%*exp(1i*-pi/2));
    unique_lines_effector=unique([pop.trial.line]);
    per_epoch=vertcat(pop.trial.epoch);
    Positions_normalized=Positions./abs(Positions);
    
    for current_epoch=1:n_epochs
        ep=find(ismember(epochs_names,epochs_names_considered(current_epoch)));
        ep_bl_name=epochs_with_baselines(ismember(epochs_with_baselines(:,2),epochs_names_considered(current_epoch)),1);
        ep_bl=ismember(epochs_names,ep_bl_name);
        
        %% bootstrapping tuning vector
        % what to do with nan FRs?
        FRs_no_baseline=vertcat(per_epoch(:,ep).FR);
        FRs=vertcat(per_epoch(:,ep).FR)-vertcat(per_epoch(:,ep_bl).FR);
        FR_vectors=FRs_no_baseline.*Positions_normalized;
        clear bootstraped
        
        %% von mises fit parameters
        idx=~isnan(FRs);
        LB=double([max([abs(FRs);1])*-1.3 -pi    1/pi]); %% Amplitude Direction Width
        UB=double([max([abs(FRs);1])*1.3  pi     10000]);%% Amplitude Direction Width
        X0=double([0                      0      4/pi]); %% Amplitude Direction Width
        
        %% significant RF fit?
        if min(N_per_pos)>=keys.cal.min_trials_per_condition
            mises_pars=lsqcurvefit(F,double(X0),double(angles(idx)),double(FRs(idx)),LB,UB,opts);
            [~,~,~,significant_RF]=stepwisefit(F(mises_pars,angles(idx)),double(FRs(idx)),'display','off');
        else
            significant_RF=false;
        end
        %unit_valid_counter=0;
        unit_valid=false;
        significant_fit=false(1,max(unique_lines_effector));
        p_anova=false(1,max(unique_lines_effector));
        mises_pars=NaN(max(unique_lines_effector),3);
        for lin=1:max(unique_lines_effector)
            clear bootstat idexes_per_position FR_vector_per_position
            tr=[pop.trial.line]==lin & ~isnan(FRs');
            mean_FR_angle=angle(nanmean(FR_vectors(tr)));
            N_per_pos_lin=hist(temp_pos_idx(tr),unique_pos_indexes);
            
            %% significant RF fit for current line?
            if min(N_per_pos_lin)>=keys.cal.min_trials_per_condition
                mises_pars(lin,:)=lsqcurvefit(F,double(X0),double(angles(tr)),double(FRs(tr)),LB,UB,opts);
                [~,~,~,inmodel]=stepwisefit(F(mises_pars(lin,:),angles(tr)),double(FRs(tr)),'display','off');
                significant_fit(lin)=inmodel;
            else
                significant_fit(lin)=false;
            end
            [~,~,pos_idx]=unique(Positions(tr));
            p_anova=anova1(FRs(tr),pos_idx,'off');
            
            %% preparing samples for bootstrap
            n_samples_per_bootstrap=NaN;
            for p=1:numel(unique_positions)
                idexes_per_position{p}=find(Positions==unique_positions(p) & tr');
                FR_vector_per_position(p)=nanmean(FR_vectors(idexes_per_position{p}));
                n_samples_per_bootstrap=min([n_samples_per_bootstrap, numel(idexes_per_position{p})]);
            end
            
            %% bootstrapping
            if n_samples_per_bootstrap>=keys.cal.min_trials_per_condition
                for k=1:n_iterations_bootstrap
                    idx_idx=cellfun(@(x) x(randsample(numel(x),10,true)),idexes_per_position,'Uniformoutput',false);
                    idx_idx=[idx_idx{:}];
                    bootstat(k,:)=[abs(nansum(FR_vectors(idx_idx(:)))/size(idx_idx,1)) angle(nansum(FR_vectors(idx_idx(:))/size(idx_idx,1))) pi/4];
                end
            else
                bootstat=single(NaN(n_iterations_bootstrap,4));
            end
            
            bootstat(:,2) = mod(bootstat(:,2)-mean_FR_angle+5*pi,2*pi)-pi+mean_FR_angle; % shift so CIs are around the mean
            bootstraped(lin).CI_amp   = [prctile(bootstat(:,1),pct1,1) nanmean(bootstat(:,1)) prctile(bootstat(:,1),pct2,1)];
            bootstraped(lin).CI_angle = [prctile(bootstat(:,2),pct1,1) nanmean(bootstat(:,2)) prctile(bootstat(:,2),pct2,1)];
            
            bootstraped(lin).FR_vectors  = bootstat(:,1).*exp(1i*bootstat(:,2));
            bootstraped(lin).mean_FR_vector  = nanmean(bootstraped(lin).FR_vectors);
            bootstraped(lin).original_mean_FR_vector  = nansum(FR_vector_per_position);
            % for plotting
            bootstraped(lin).mean_amp   = nanmean(bootstat(:,1));
            bootstraped(lin).mean_angle = nanmean(bootstat(:,2));
            bootstraped(lin).mean_kappa = nanmean(bootstat(:,3));
            % for CIs
            bootstraped(lin).raw_amp   = bootstat(:,1);
            bootstraped(lin).raw_angle = bootstat(:,2);
            
            %% annoying part to see if there is "tuning" in the first place
            projected_vectors=bootstat(:,1).*cos(bootstat(:,2)-angle(nansum(bootstat(:,1).*exp(1i*bootstat(:,2)))));
            if prctile(projected_vectors,100,1)>0 && prctile(projected_vectors,0,1)>0 %% all bootstrapped endpoints have to be on the same side for the unit to be considered having tuning
                unit_valid=true;
            end
        end
        
        %% SIGNIFICANCE with CIs
        switch significance_measure
            case 'unpaired'
                CI1=vertcat(bootstraped.CI_amp);
                Lower_bound=CI1(:,1);
                Mean_values=CI1(:,2);
                Upper_bound=CI1(:,3);
                %sign_diff = ~(bsxfun( @le, Lower_bound, Mean_values' ) & bsxfun( @ge, Upper_bound, Mean_values' ));
                sign_diff = ~(bsxfun( @le, Lower_bound, Upper_bound' ) & bsxfun( @ge, Upper_bound, Lower_bound' ));
                significant_gain=any(any(sign_diff));
                
                for c_s=1:numel(bootstraped)+1 %% shifting always one by 2 pi to find overlaps around 0
                    CI1=vertcat(bootstraped.CI_angle);
                    if c_s<numel(bootstraped)
                        CI1(c_s,:)= CI1(c_s,:)+2*pi;
                    end
                    Lower_bound=CI1(:,1);
                    Mean_values=CI1(:,2);
                    Upper_bound=CI1(:,3);
                    %sign_diff = ~(bsxfun( @le, Lower_bound, Mean_values' ) & bsxfun( @ge, Upper_bound, Mean_values' ));
                    sign_diff = ~(bsxfun( @le, Lower_bound, Upper_bound' ) & bsxfun( @ge, Upper_bound, Lower_bound' ));
                    %
                    significant_shift(c_s)=any(any(sign_diff));
                end
                significant_shift=all(significant_shift); % only true if CIs are not overlapping for any 2*pi shift
                
            case 'paired'
                %% significance with CIs of differences
                temp1=true(numel(bootstraped));
                [temp2 temp3]=ind2sub(size(temp1),find(temp1));
                comparison_matrix=[temp2(temp2<temp3) temp3(temp2<temp3)];
                for c_s=1:size(comparison_matrix,1) %% shifting always one by 2 pi to find overlaps around 0
                    tt=bootstraped(comparison_matrix(c_s,1)).raw_amp-bootstraped(comparison_matrix(c_s,2)).raw_amp;
                    CI1_diff(c_s,:)=[prctile(tt,pct1,1) prctile(tt,pct2,1)];
                    
                    tt=bootstraped(comparison_matrix(c_s,1)).raw_angle-bootstraped(comparison_matrix(c_s,2)).raw_angle;
                    %tt=[tt tt+pi tt-pi tt+2*pi tt-2*pi tt*3+pi tt-3*pi tt+4*pi tt-4*pi];
                    tt=[tt tt+2*pi tt-2*pi tt+4*pi tt-4*pi];% ti ti+2*pi ti-2+pi ti+4*pi ti-4*pi];
                    [~,totest_for_CI_indexes]=min(abs(tt),[],2);
                    tt=tt(sub2ind(size(tt),1:size(tt,1),totest_for_CI_indexes'))';
                    CI2_diff(c_s,:)=[prctile(tt,pct1,1) prctile(tt,pct2,1)];
                end
                significant_gain=any(CI1_diff(:,1)>0 | CI1_diff(:,2)<0);
                significant_shift=any(CI2_diff(:,1)>0 | CI2_diff(:,2)<0);
        end
        
        %% storing parameters
        population(u).RF_per_epoch(ep).significant_gain=significant_gain;
        population(u).RF_per_epoch(ep).significant_shift=significant_shift;
        population(u).RF_per_epoch(ep).significant_fit=significant_fit;
        population(u).RF_per_epoch(ep).p_anova=p_anova;
        population(u).RF_per_epoch(ep).significant_RF=significant_RF;
        population(u).RF_per_epoch(ep).unit_valid=unit_valid;
        population(u).RF_per_epoch(ep).mises_pars=mises_pars;
        population(u).RF_per_epoch(ep).gain_signal_change=(max([bootstraped.mean_amp])-min([bootstraped.mean_amp]))/max([bootstraped.mean_amp]);        
        RF_per_epoch(ep).bootstraped=bootstraped;
    end
    
    if do_plot
        linespec.visible='off';
        for current_epoch=1:n_epochs
            ep=find(ismember(epochs_names,epochs_names_considered(current_epoch)));
            ep_bl_name=epochs_with_baselines(ismember(epochs_with_baselines(:,2),epochs_names_considered(current_epoch)),1);
            ep_bl=ismember(epochs_names,ep_bl_name);
            
            %% bootstrapping tuning vector
            % what to do with nan FRs?
            figure_handle=figure;
            %sp(1,ep)=subplot(n_rows_columns,n_rows_columns,ep);
            bootstraped=RF_per_epoch(ep).bootstraped;
            all_bt=vertcat(bootstraped);
            max_FR_vector=max(abs(vertcat(all_bt.FR_vectors)));
            SIG=population(u).RF_per_epoch(ep);
            switch method
                case 'FR_vector'
                    %% tuning vector bootstrapping plot
                    title_part=[keys.EPOCHS_PER_TYPE{type}{ep,1} 'val_RF_gain_shift ' num2str([SIG.unit_valid SIG.significant_RF SIG.significant_gain SIG.significant_shift])];
                    filenamepart='FR_vector';
                    polar(0, max_FR_vector,linespec);
                    hold on
                    for lin=1:max(unique_lines_effector)
                        scatter(real([bootstraped(lin).FR_vectors]),imag([bootstraped(lin).FR_vectors]),5,pop.PSTH_perpos_colors(lin,:));
                        FR_vector=bootstraped(lin).mean_FR_vector;
                        original_FR_vector=bootstraped(lin).original_mean_FR_vector;
                        FR_Amp_CI=bootstraped(lin).CI_amp([1,end])*FR_vector/abs(FR_vector);
                        FR_Ang_CI=exp(1i*bootstraped(lin).CI_angle)*abs(FR_vector);
                        FR_Ang_CI_rad=exp(1i*[bootstraped(lin).CI_angle(1):pi/100:bootstraped(lin).CI_angle(end) bootstraped(lin).CI_angle(end)])*abs(FR_vector);
                        plot(real(original_FR_vector),imag(original_FR_vector),'h','markersize',markersize,'markeredgecolor','k','color',pop.PSTH_perpos_colors(lin,:),'Markerfacecolor',pop.PSTH_perpos_colors(lin,:));
                        plot(real(FR_vector),imag(FR_vector),'s','markersize',markersize,'markeredgecolor','k','color',pop.PSTH_perpos_colors(lin,:),'Markerfacecolor',pop.PSTH_perpos_colors(lin,:));
                        plot(real(FR_Amp_CI),imag(FR_Amp_CI),'d-','linewidth',3,'markersize',markersize,'markeredgecolor','k','color',pop.PSTH_perpos_colors(lin,:));
                        plot(real(FR_Ang_CI(1)),imag(FR_Ang_CI(1)),'d','linewidth',3,'markersize',markersize,'markeredgecolor','k','color',pop.PSTH_perpos_colors(lin,:));
                        plot(real(FR_Ang_CI(end)),imag(FR_Ang_CI(end)),'d','linewidth',3,'markersize',markersize,'markeredgecolor','k','color',pop.PSTH_perpos_colors(lin,:));
                        line([real(FR_Ang_CI_rad)],[imag(FR_Ang_CI_rad)],'linewidth',3,'color',pop.PSTH_perpos_colors(lin,:))
                    end
                    
                case 'Von_mises_fit'
                    %% von mises fitting plot
                    OR=num2str(SIG.significant_fit(1));
                    BL='no';
                    RE='no';
                    if numel(SIG.significant_fit)>1
                        BL=num2str(SIG.significant_fit(2));
                    end
                    if numel(SIG.significant_fit)>2
                        RE=num2str(SIG.significant_fit(3));
                    end
                    title_part=[keys.EPOCHS_PER_TYPE{type}{ep,1} ' tot_' num2str(SIG.significant_RF) '_or_' OR '_bl_' BL  '_re_' RE];
                    filenamepart='tuning_profile';
                    FRs=vertcat(per_epoch(:,ep).FR)-vertcat(per_epoch(:,ep_bl).FR);
                    max_FR=max([FRs;1]);
                    min_FR=min([FRs;0]);
                    minmaxdiff=max_FR-min_FR;
                    if numel(FRs)>10
                        unique_FRs=unique([FRs;0]);
                        mean_FR_step=nanmean(diff(unique_FRs))/2;
                    else
                        mean_FR_step=1;
                    end
                    if isempty(mean_FR_step)
                        mean_FR_step=1;
                    end
                    ylim([min_FR-minmaxdiff*0.1 max_FR+minmaxdiff*0.1]);
                    hold on
                    for lin=1:max(unique_lines_effector)
                        tr=[pop.trial.line]==lin;
                        current_trial_indexes=find(tr' & ~isnan(FRs));
                        current_FRs=FRs(current_trial_indexes);
                        angles_ep=angles(current_trial_indexes);
                        coefficients=population(u).RF_per_epoch(ep).mises_pars(lin,:); %[bootstraped(lin).mean_amp bootstraped(lin).mean_angle bootstraped(lin).mean_kappa];
                        scatter(angles_ep,current_FRs+rand(size(current_FRs))*mean_FR_step,10,pop.PSTH_perpos_colors(lin,:),'o','filled')
                        plot(angles_to_plot,F(coefficients,angles_to_plot),'color',pop.PSTH_perpos_colors(lin,:),'linewidth',3);
                    end
            end
            ph_title_and_save(figure_handle,[filename ' ' title_part ' ' filenamepart],[plot_title ' ' title_part filenamepart],keys);
        end
    end
end
save([keys.path_to_save 'population_' method],'population');
end