function ph_error_signal(varargin)

keys=struct;
project=varargin{1};
keys=ph_general_settings(project,keys);
project_specific_settings=[keys.db_folder filesep project filesep 'ph_project_settings.m'];
run(project_specific_settings);

if nargin>1
    keys.project_versions=varargin{2};
end
for f=1:numel(keys.project_versions) % running multiple versions of the same project at once !
    if ~isempty(keys.project_versions{f})
        keys.project_version=keys.project_versions{f};
    end
    version_folder=keys.project_version;
    keys.version_specific_settings=[keys.db_folder project filesep keys.project_version filesep 'ph_project_version_settings.m'];
    run(keys.version_specific_settings);
    keys.project_version=version_folder;
    keys.monkeys=keys.batching.monkeys;
    if keys.batching.combine_monkeys
        keys.batching.monkeys={''};
    end
    
    for m=1:numel(keys.batching.monkeys)
        keys.monkey=keys.batching.monkeys{m};
        keys.anova_table_file=[keys.basepath_to_save keys.project_version '\tuning_table_combined_CI.mat'];
        population=ph_load_population([keys.basepath_to_save keys.project_version],['population_' keys.monkey]);
        population=ph_assign_perturbation_group(keys,population);
        population=ph_epochs(population,keys);
        
        keys.tt.IC_for_criterion='in';
        keys.tt.tasktypes={'Msac_opt'};
        keys.tt.hands=0;
        keys.tt.choices=0;
        keys.tt.perturbations=0;
        
        keys.arrangement='opt';
        %
        %         % specific subset
        %         keys.tt.combine_tuning_properties  = {'fix_tar_mono','in_AH_Fhol_position_Msac_fix','in_AH_Thol_position_Msac_tar','in_AH_Fhol_gaze_modulation_Msac_fix'};
        %         keys.tt.selection                  = {'fix_tar_mono','truefalsenonmonotoneus'};
        
        tasktype='Msac_opt';
        for t=1:numel(keys.batching.targets)
            target=keys.batching.targets{t};
            
            for subregion=1:keys.batching.n_Subregions
                %population_selection={};
                if keys.batching.Subregions_separately
                    keys.tt.selection            = {'target',target;'Subregion', subregion};
                else
                    keys.tt.selection             = {'target',target};
                end
                
                
                %% from here make subfunction
                keys.path_to_save=[keys.drive filesep keys.basepath_to_save filesep keys.project_version filesep 'error_signals' filesep];
                if ~exist(keys.path_to_save,'dir')
                    mkdir([keys.drive filesep keys.basepath_to_save filesep keys.project_version ], 'error_signals');
                end
                
                %% load tuning table and reduce population
                [tuning_per_unit_table]                 = ph_load_extended_tuning_table(keys);
                [tuning_per_unit_table, Sel_for_title]  = ph_reduce_tuning_table(tuning_per_unit_table,keys);
                
                all_IDs={population.unit_ID};
                valid_IDs=ismember(all_IDs,tuning_per_unit_table(:,1));
                population=population(valid_IDs);
                
                all_trials=[population.trial];
                u_types     =unique([all_trials.type]);
                u_effectors =unique([all_trials.effector]);
                type_effectors      = combvec(u_types,u_effectors)';
                for te=1:size(type_effectors,1) %typ=unique(per_trial.types)
                    clear coef pval coef_res pval_res coef_dir pval_dir coef_res_dir pval_res_dir sig_idx
                    
                    typ=type_effectors(te,1);
                    eff=type_effectors(te,2);
                    [~, type_effector_short]=get_type_effector_name(typ,eff);
                    idx_existing=DAG_find_column_index(tuning_per_unit_table,['existing_' 'in_AH_' type_effector_short '_' keys.arrangement(1:3)]);
                    idx_unit_ID=DAG_find_column_index(tuning_per_unit_table,['unit_ID']);
                    units_existing= [true; ismember(vertcat(tuning_per_unit_table{2:end,idx_existing}),true)];
                    unit_IDs_te=tuning_per_unit_table(units_existing,idx_unit_ID);
                    
                    population_te=population(ismember({population.unit_ID},unit_IDs_te));
                    unit_IDs_te={population_te.unit_ID};
                    % here we go
                    for u=1:numel(population_te)
                        pop=ph_LR_to_CI(keys,population_te(u));
                        pop=ph_epochs(pop,keys);
                        trials=[pop.trial];
                        trials=trials([trials.success] & [trials.accepted] & [trials.effector]==eff & [trials.type]==typ);
                        
                        unique_targets=unique([trials.tar_pos]);
                        unique_targets(isnan(unique_targets))=[];
                        idx_valid=true(size(unique_targets));
                        for tar_temp=1:numel(unique_targets)
                            idx_valid=idx_valid & ~(1:numel(unique_targets)<tar_temp & abs(unique_targets-unique_targets(tar_temp))<2);
                        end
                        unique_targets=unique_targets(idx_valid);
                        per_epoch=vertcat(trials.epoch);
                        per_epoch_residuals=NaN(size(per_epoch));
                        for ep=1:size(per_epoch,2)
                            for tar=1:numel(unique_targets)
                                tar_idx=abs([trials.tar_pos]-unique_targets(tar))<2;
                                %off=abs([trials(tar_idx).sac_off]);
                                FR_residuals=[per_epoch(tar_idx,ep).FR]-nanmean([per_epoch(tar_idx,ep).FR]);
                                
                                % correct by mean offset
                                temp_offsets=[trials(tar_idx).sac_off]-mean([trials(tar_idx).sac_off]);
                                temptemp=num2cell(temp_offsets);
                                [trials(tar_idx).sac_off]=deal(temptemp{:});
                                per_epoch_residuals(tar_idx,ep)=FR_residuals;
                                
                                
                                temptemp=num2cell(abs(abs(temp_offsets).*cos(angle(temp_offsets)-angle([trials(tar_idx).tar_pos]))));
                                [trials(tar_idx).sac_off_dir]=deal(temptemp{:});
                                
                                
                                %                                 temptemp=num2cell([trials(tar_idx).sac_off]-mean([trials(tar_idx).sac_off]));
                                %                                 [trials(tar_idx).sac_off]=deal(temptemp{:});
                            end
                            
                            
                            
                            off=abs([trials.sac_off]);
                            FR=[per_epoch(:,ep).FR];
                            [coef(u,ep), pval(u,ep)] = corr(FR(:),off(:),'type','Spearman');
                            
                            
                            FR=[per_epoch_residuals(:,ep)];
                            [coef_res(u,ep), pval_res(u,ep)] = corr(FR(:),off(:),'type','Spearman');
                            
                            
                            off=abs([trials.sac_off_dir]);
                            FR=[per_epoch(:,ep).FR];
                            [coef_dir(u,ep), pval_dir(u,ep)] = corr(FR(:),off(:),'type','Spearman');
                            
                            
                            FR=[per_epoch_residuals(:,ep)];
                            [coef_res_dir(u,ep), pval_res_dir(u,ep)] = corr(FR(:),off(:),'type','Spearman');
                        end
                    end
                    epochnames=keys.EPOCHS_PER_TYPE{typ}(:,1);
                    for ep=1:size(coef,2)                        
                        sig_col=DAG_find_column_index(tuning_per_unit_table,['in_' epochnames{ep} '_spaceLR_' tasktype]);
                        sig_units=cellfun(@(x) strcmp(x,'CS') || strcmp(x,'IS'),tuning_per_unit_table(:,sig_col));
                        sig_idx(:,ep)=ismember(unit_IDs_te,tuning_per_unit_table(sig_units,1))';
                        
                        histogram_per_epoch(coef,pval,ep,sig_idx,Sel_for_title,type_effector_short,epochnames{ep},' absolute errors, residuals',keys)
                        histogram_per_epoch(coef_res,pval_res,ep,sig_idx,Sel_for_title,type_effector_short,epochnames{ep},' absolute errors, residuals of both FR ERR',keys)
                        histogram_per_epoch(coef_dir,pval_dir,ep,sig_idx,Sel_for_title,type_effector_short,epochnames{ep},' directional errors, residuals',keys)
                        histogram_per_epoch(coef_res_dir,pval_res_dir,ep,sig_idx,Sel_for_title,type_effector_short,epochnames{ep},' directional errors, residuals of both FR ERR',keys)  
                    end
                    
                    %% add significances here as well?
                    
                    histogram_summary(pval,sig_idx,Sel_for_title,type_effector_short,epochnames,' absolute residual ERR errors summary',keys)
                    histogram_summary(pval_res,sig_idx,Sel_for_title,type_effector_short,epochnames,' absolute residual FR errors summary',keys)
                    histogram_summary(pval_dir,sig_idx,Sel_for_title,type_effector_short,epochnames,' directional residual ERR errors summary',keys)
                    histogram_summary(pval_res_dir,sig_idx,Sel_for_title,type_effector_short,epochnames,' directional residual FRERR errors summary',keys)
                end
                
                
            end
        end
    end
end
end

function [N,P,R]=scatter_per_sig(varargin)
unit=varargin{1};
idx=varargin{2};
N=sum(idx);
if N>2
    [R, P] = corr([unit(idx).FR_fix_C]'-[unit(idx).FR_fix_I]',[unit(idx).FR_tar_C]'-[unit(idx).FR_tar_I]');
else
    R=0;
    P=0;
end
scatter([unit(idx).FR_fix_C]-[unit(idx).FR_fix_I],[unit(idx).FR_tar_C]-[unit(idx).FR_tar_I],varargin{3:end});

end

function pop=ph_LR_to_CI_raw(pop,target)
%% assigning contra and ipsi instead of left/right, dependent on the recorded target (which contains hemisphere information in the end)
% reference is left hemisphere, meaning that right becomes contra and left ipsi
% that is why hemifields, positions and hands need to be inverted for right hemisphere targets.
% for right recording sites LR->CI  (hemifield & positions...?)
% positive space and hand==2 become contra

% Left_hemisphere_targets={'dPulv_l','pdSTS_L','FST_L','pTPO_L','LIP_L','MIP_L','unknown'};
% Right_hemisphere_targets={'dPulv_r','PUL_r','PUL_R','MIP_R','LIP_R'};

if strcmpi(target(end-1:end),'_R')
    temp=num2cell(real([pop.trial.fix_pos])*-1+1i*imag([pop.trial.fix_pos]));
    [pop.trial.fix_pos]=deal(temp{:});
    temp=num2cell(real([pop.trial.tar_pos])*-1+1i*imag([pop.trial.tar_pos]));
    [pop.trial.tar_pos]=deal(temp{:});
    hand1=[pop.trial.reach_hand]==2;
    hand2=[pop.trial.reach_hand]==1;
    [pop.trial(hand1).reach_hand]=deal(1);
    [pop.trial(hand2).reach_hand]=deal(2);
end
end


function histogram_per_epoch(coef,pval,ep,sig_idx,Sel_for_title,type_effector_short,epochname,title,keys)

figure_handle=figure;
hold on
minmax=[floor(min(min(coef))*10)/10 ceil(max(max(coef))*10)/10];
bins=[minmax(1):0.1:minmax(2)];


histsigsig=hist(coef(pval(:,ep)<0.05 & sig_idx(:,ep),ep), bins);
histsig=hist(coef(pval(:,ep)<0.05 & ~sig_idx(:,ep),ep), bins) + histsigsig;
histnonsigcor=hist(coef(pval(:,ep)>=0.05 & sig_idx(:,ep),ep), bins)+histsig +histsigsig;
histnonsig=hist(coef(pval(:,ep)>=0.05 & ~sig_idx(:,ep),ep), bins)+histnonsigcor +histsig +histsigsig;
bar(bins,histnonsig,'facecolor','w');
bar(bins,histnonsigcor,'facecolor',[0.5 0.5 0.5]);
bar(bins,histsig,'facecolor','k');
bar(bins,histsigsig,'facecolor','r');
xlabel('R');
ylabel('N');


wanted_size=[50 30];
set(figure_handle, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
export_fig([ keys.path_to_save Sel_for_title{:} ' ' type_effector_short ' ' epochname title], '-pdf','-transparent') % pdf by run
close(gcf);
end


function histogram_summary(pval,sig_idx,Sel_for_title,type_effector_short,epochnames,title,keys)

figure_handle=figure;
hold on
sigsig=sum(pval<0.05  &  sig_idx,1);
sigcor=sum(pval<0.05  & ~sig_idx,1);
sigtun=sum(pval>=0.05 &  sig_idx,1);
nonsig=sum(pval>=0.05 & ~sig_idx,1);

bh=bar(nonsig+sigtun+sigcor+sigsig,'facecolor','w');
bh=bar(sigtun+sigcor+sigsig,'facecolor',[0.5 0.5 0.5]);
bh=bar(sigcor+sigsig,'facecolor','k');
bh=bar(sigsig,'facecolor','r');

x_ticks=1:numel(sigsig);
set(gca,'xtick',x_ticks,'xticklabel',epochnames);

wanted_size=[50 30];
set(figure_handle, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
export_fig([ keys.path_to_save Sel_for_title{:} ' ' type_effector_short title], '-pdf','-transparent') % pdf by run
close(gcf);
end
