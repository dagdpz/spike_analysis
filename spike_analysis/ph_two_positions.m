function ph_two_positions(varargin)

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
    
    
    keys.basepath_to_save=['Projects' filesep project filesep 'ephys' filesep];
    seed_filename=[keys.drive keys.basepath_to_save keys.project_version filesep 'seed.mat'];
    if exist(seed_filename,'file');
        load(seed_filename);
        rng(seed);
    else
        seed=rng;
        save(seed_filename,'seed');
    end
    
    for counter=1:2
        
        switch counter
            case 1
                keys.FR_or_idx='FR';
            case 2
                keys.FR_or_idx='index';
        end
        
        for m=1:numel(keys.batching.monkeys)
            keys.monkey=keys.batching.monkeys{m};
            keys.anova_table_file=[keys.drive '\Projects\' project '\ephys\' keys.project_version '\tuning_table_combined_CI.mat'];
            population=ph_load_population([keys.drive filesep keys.basepath_to_save filesep keys.project_version],['population_' keys.monkey]);
            population=ph_assign_perturbation_group(keys,population);
            population=ph_epochs(population,keys);
            
            keys.tt.IC_for_criterion='in';
            keys.tt.tasktypes={'Msac_tar'};
            keys.tt.hands=0;
            keys.tt.choices=0;
            keys.tt.perturbations=0;
            
            % specific subset
%             keys.tt.combine_tuning_properties  = {'fix_tar_mono','in_AH_Fhol_position_Msac_fix','in_AH_Thol_position_Msac_tar','in_AH_Fhol_gaze_modulation_Msac_fix'};
%             keys.tt.selection                  = {'fix_tar_mono','truefalsenonmonotoneus'};
            
            [tuning_per_unit_table]                 = ph_load_extended_tuning_table(keys);
            [tuning_per_unit_table, Sel_for_title]  = ph_reduce_tuning_table(tuning_per_unit_table,keys);
            
            all_IDs={population.unit_ID};
            valid_IDs=ismember(all_IDs,tuning_per_unit_table(:,1));
            population=population(valid_IDs);
            
            for t=1:numel(keys.batching.targets)
                target=keys.batching.targets{t};
                
                for subregion=1:keys.batching.n_Subregions
                    %population_selection={};
                    if keys.batching.Subregions_separately
                        population_selection            = {'target',target;'Subregion', subregion};
                    else
                        population_selection            = {'target',target};
                    end
                    
                    %% load tuning table and reduce population
                    
                    
                    % here we go
                    for u=1:numel(population)
                        pop=ph_LR_to_CI_raw(population(u),population(u).target);
                        
                        trials=[pop.trial];
                        trials=trials([trials.completed] & [trials.accepted]);
                        fix_idxI=abs(real([trials.fix_pos]+15))<2;
                        fix_idxS=abs(real([trials.fix_pos]))<2;
                        fix_idxC=abs(real([trials.fix_pos]-15))<2;
                        
                        
                        tar_idxI=abs([trials.tar_pos]-1i*imag([trials.fix_pos])+15)<2;
                        tar_idxS=abs([trials.tar_pos]-1i*imag([trials.fix_pos]))<2;
                        tar_idxC=abs([trials.tar_pos]-1i*imag([trials.fix_pos])-15)<2;
                        
                        
                        
%                         tar_idxI=abs(real([trials.tar_pos]-1i*imag([trials.fix_pos])+15))<2;
%                         tar_idxC=abs(real([trials.tar_pos]-1i*imag([trials.fix_pos])-15))<2;
                        
                        per_epoch=vertcat(trials.epoch);
                        fix_FR=[per_epoch(:,~cellfun(@isempty,strfind({per_epoch(1,:).state},'Fhol'))).FR];
                        tar_FR=[per_epoch(:,~cellfun(@isempty,strfind({per_epoch(1,:).state},'Thol'))).FR];
                        
                        unit(u).FR_fix_I=nanmean(fix_FR(fix_idxI));
                        unit(u).FR_fix_C=nanmean(fix_FR(fix_idxC));
                        unit(u).FR_tar_I=nanmean(tar_FR(tar_idxI));
                        unit(u).FR_tar_C=nanmean(tar_FR(tar_idxC));
                        
                        
                        if strcmp(keys.FR_or_idx,'index')
                            
                            tmpFRI1=unit(u).FR_fix_I;
                            tmpFRC1=unit(u).FR_fix_C;
                            tmpFRI2=unit(u).FR_tar_I;
                            tmpFRC2=unit(u).FR_tar_C;
                            unit(u).FR_fix_I=unit(u).FR_fix_I/(tmpFRI1 + tmpFRC1);
                            unit(u).FR_fix_C=unit(u).FR_fix_C/(tmpFRI1 + tmpFRC1);
                            unit(u).FR_tar_I=unit(u).FR_tar_I/(tmpFRI2 + tmpFRC2);
                            unit(u).FR_tar_C=unit(u).FR_tar_C/(tmpFRI2 + tmpFRC2);
                        end
%                         
% %% perfroming 1000 iterations of repicking the same number of trials for fixation
% for it=1:1000
%     tmp_fix_idx_I=find(fix_idxI);
%     tmp_fix_idx_C=find(fix_idxC);
%    p_fix(it)= ttest2(fix_FR(randsample(tmp_fix_idx_C,sum(tar_idxC))),fix_FR(randsample(tmp_fix_idx_I,sum(tar_idxI))));
% end
% 
%                         unit(u).fix_sig=median(p_fix);                        %ttest2(fix_FR(fix_idxC),fix_FR(fix_idxI));
%                         unit(u).tar_sig=ttest2(tar_FR(tar_idxC),tar_FR(tar_idxI));
                        
% 
%     tmp_fix_idx_I=find(fix_idxI);
%     tmp_fix_idx_C=find(fix_idxC);
%                     unit(u).fix_sig=ttest2(fix_FR(randsample(tmp_fix_idx_C,sum(tar_idxC))),fix_FR(randsample(tmp_fix_idx_I,sum(tar_idxI))));                       %ttest2(fix_FR(fix_idxC),fix_FR(fix_idxI));
%                         unit(u).tar_sig=ttest2(tar_FR(tar_idxC),tar_FR(tar_idxI));
any_fix=fix_idxI | fix_idxS | fix_idxC;
any_tar=tar_idxI | tar_idxS | tar_idxC;

p_temp= anova1(fix_FR(any_fix),[fix_idxI(any_fix) + 2*fix_idxS(any_fix) + 3*fix_idxC(any_fix)],'off');

unit(u).fix_sig=p_temp<0.05;
p_temp= anova1(tar_FR(any_tar),[tar_idxI(any_tar) + 2*tar_idxS(any_tar) + 3*tar_idxC(any_tar)],'off');
unit(u).tar_sig= p_temp<0.05;
% [p,anovatab,stats] = anova1(fix_FR,[];)
%                         
%                     unit(u).fix_sig=ttest2(fix_FR(randsample(tmp_fix_idx_C,sum(tar_idxC))),fix_FR(randsample(tmp_fix_idx_I,sum(tar_idxI))));                       %ttest2(fix_FR(fix_idxC),fix_FR(fix_idxI));
%                         unit(u).tar_sig=ttest2(tar_FR(tar_idxC),tar_FR(tar_idxI));

%                         [~,unit(u).fix_sig]=ranksum(fix_FR(fix_idxC),fix_FR(fix_idxI));
%                          [~,unit(u).tar_sig]=ranksum(tar_FR(tar_idxC),tar_FR(tar_idxI));
                    end
                    
                    figure_handle=figure;
                    hold on
                    c=0;
                    idx=[unit.fix_sig]==1 & [unit.tar_sig]==1; c=c+1;
                    [N(c),P(c),R(c)]=scatter_per_sig(unit,idx,'o','filled');
                    idx=[unit.fix_sig]~=1 & [unit.tar_sig]==1; c=c+1;
                    [N(c),P(c),R(c)]=scatter_per_sig(unit,idx,'^','filled');
                    idx=[unit.fix_sig]==1 & [unit.tar_sig]~=1; c=c+1;
                    [N(c),P(c),R(c)]=scatter_per_sig(unit,idx,'>','filled');
                    idx=[unit.fix_sig]~=1 & [unit.tar_sig]~=1; c=c+1;
                    [N(c),P(c),R(c)]=scatter_per_sig(unit,idx,'o');
                    idx=true(size(unit)); c=c+1;
                    [N(c),P(c),R(c)]=scatter_per_sig(unit,idx,'o');
                    ylabel(['Final gaze preference Contra-Ipsi ' keys.FR_or_idx])
                    xlabel(['Initial gaze preference Contra-Ipsi ' keys.FR_or_idx])
                    title('comparing two positions only')
                    
                    
                    axes_limits=[min([get(gca,'xlim') get(gca,'ylim')]) max([get(gca,'xlim') get(gca,'ylim')])];
                    set(gca,'xlim',axes_limits,'ylim',axes_limits);
                    axis square;
                    
                    y_lim=ylim;
                    x_lim=xlim;
                    
                    text(repmat(x_lim(1)+1*diff(x_lim)/10,1,6),y_lim(2)-[1,2,3,4,5,6].*diff(y_lim)/20,{'N';num2str(N')});
                    text(repmat(x_lim(1)+2*diff(x_lim)/10,1,6),y_lim(2)-[1,2,3,4,5,6].*diff(y_lim)/20,{'R';num2str(round(R'*100)/100)});
                    text(repmat(x_lim(1)+3*diff(x_lim)/10,1,6),y_lim(2)-[1,2,3,4,5,6].*diff(y_lim)/20,{'P';num2str(round(P'*1000)/1000)});
                    
                    
                    
                    wanted_size=[50 30];
                    set(figure_handle, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
                    
                    
                    export_fig([keys.drive filesep keys.basepath_to_save filesep keys.project_version filesep 'scatter' filesep 'comparing two positions only ' keys.FR_or_idx], '-pdf','-transparent') % pdf by run
                    close(gcf);
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
