function ph_central_gaze_only(varargin)

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
        keys.anova_table_file=[keys.basepath_to_save keys.project_version filesep 'tuning_table_combined_CI.mat'];
        population=ph_load_population([keys.basepath_to_save keys.project_version],['population_' keys.monkey]);
        population=ph_assign_perturbation_group(keys,population);
        population=ph_epochs(population,keys);
        
        keys.tt.IC_for_criterion='in';
        keys.tt.tasktypes={'Msac_mov'};
        keys.tt.selection={'target','dPul'};
        keys.tt.hands=0;
        keys.tt.choices=0;
        keys.tt.perturbations=0;
        
        
        for counter=1:4
            
            switch counter
                case 1
                    keys.epochs_central_gaze={'Cue','Cue'};
                    keys.arrangements_central_gaze={'fix','mov'};
                    keys.FR_or_idx_central_gaze='FR';
                case 2
                    keys.epochs_central_gaze={'Cue','Cue'};
                    keys.arrangements_central_gaze={'fix','mov'};
                    keys.FR_or_idx_central_gaze='index';
                case 3
                    keys.epochs_central_gaze={'Fhol','Cue'};
                    keys.arrangements_central_gaze={'fix','mov'};
                    keys.FR_or_idx_central_gaze='FR';
                case 4
                    keys.epochs_central_gaze={'Fhol','Cue'};
                    keys.arrangements_central_gaze={'fix','mov'};
                    keys.FR_or_idx_central_gaze='index';
                    
            end
                
        ep1=keys.epochs_central_gaze{1};
        ep2=keys.epochs_central_gaze{2};
        ar1=keys.arrangements_central_gaze{1};
        ar2=keys.arrangements_central_gaze{2};
% 
%         % specific subset
%         keys.tt.combine_tuning_properties  = {'fix_tar_mono','in_AH_Fhol_position_Msac_fix','in_AH_Thol_position_Msac_tar','in_AH_Fhol_gaze_modulation_Msac_fix'};
%         keys.tt.selection                  = {'fix_tar_mono','truefalsenonmonotoneus'};

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
                    central_gaze_idx=abs(real([trials.fix_pos]))<2;
                    contra_gaze_idx=real([trials.fix_pos])<-10;
                    ipsi_gaze_idx=real([trials.fix_pos])>10;
                    
                    idxI.fix=abs(real([trials.fix_pos]+15))<2;
                    idxC.fix=abs(real([trials.fix_pos]-15))<2;       
                    idxI.mov=real([trials.tar_pos]-[trials.fix_pos])<-2 & central_gaze_idx ;
                    idxC.mov=real([trials.tar_pos]-[trials.fix_pos])>2 & central_gaze_idx ;
                    idxI.tar=real([trials.tar_pos])<-2 & central_gaze_idx ;
                    idxC.tar=real([trials.tar_pos])>2  & central_gaze_idx ;
                    
                    
                    idxIc.mov=real([trials.tar_pos]-[trials.fix_pos])<-2 & contra_gaze_idx ;
                    idxCc.mov=real([trials.tar_pos]-[trials.fix_pos])>2 & contra_gaze_idx ;
                    idxIc.tar=real([trials.tar_pos])<-2 & contra_gaze_idx ;
                    idxCc.tar=real([trials.tar_pos])>2  & contra_gaze_idx ;                    
                    
                    idxIi.mov=real([trials.tar_pos]-[trials.fix_pos])<-2 & ipsi_gaze_idx ;
                    idxCi.mov=real([trials.tar_pos]-[trials.fix_pos])>2 & ipsi_gaze_idx ;
                    idxIi.tar=real([trials.tar_pos])<-2 & ipsi_gaze_idx ;
                    idxCi.tar=real([trials.tar_pos])>2  & ipsi_gaze_idx ;
                    
                    per_epoch=vertcat(trials.epoch);
                    ep1_FR=[per_epoch(:,~cellfun(@isempty,strfind({per_epoch(1,:).state},ep1))).FR];
                    ep2_FR=[per_epoch(:,~cellfun(@isempty,strfind({per_epoch(1,:).state},ep2))).FR];
                    
                    unit(u).FRI_1=nanmean(ep1_FR(idxI.(ar1)));
                    unit(u).FRC_1=nanmean(ep1_FR(idxC.(ar1)));
                    unit(u).FRI_2=nanmean(ep2_FR(idxI.(ar2)));
                    unit(u).FRC_2=nanmean(ep2_FR(idxC.(ar2)));
                    
                    
                    unit(u).FRIc_2=nanmean(ep2_FR(idxIc.(ar2)));
                    unit(u).FRIi_2=nanmean(ep2_FR(idxIi.(ar2)));
                    unit(u).FRCc_2=nanmean(ep2_FR(idxCc.(ar2)));
                    unit(u).FRCi_2=nanmean(ep2_FR(idxCi.(ar2)));
                    
                    if strcmp(keys.FR_or_idx_central_gaze,'index')
                        tmpFRI1=unit(u).FRI_1;
                        tmpFRC1=unit(u).FRC_1;
                        tmpFRI2=unit(u).FRI_2;
                        tmpFRC2=unit(u).FRC_2;
                        unit(u).FRI_1=unit(u).FRI_1/(tmpFRI1 + tmpFRC1);
                        unit(u).FRC_1=unit(u).FRC_1/(tmpFRI1 + tmpFRC1);
                        unit(u).FRI_2=unit(u).FRI_2/(tmpFRI2 + tmpFRC2);
                        unit(u).FRC_2=unit(u).FRC_2/(tmpFRI2 + tmpFRC2);
                        
                        
                        tmpFRIi2=unit(u).FRIi_2;
                        tmpFRCi2=unit(u).FRCi_2;
                        unit(u).FRIi_2=unit(u).FRIi_2/(tmpFRIi2 + tmpFRCi2);
                        unit(u).FRCi_2=unit(u).FRCi_2/(tmpFRIi2 + tmpFRCi2);
                        tmpFRIc2=unit(u).FRIc_2;
                        tmpFRCc2=unit(u).FRCc_2;
                        unit(u).FRIc_2=unit(u).FRIc_2/(tmpFRIc2 + tmpFRCc2);
                        unit(u).FRCc_2=unit(u).FRCc_2/(tmpFRIc2 + tmpFRCc2);
                    end
                    
                    unit(u).ep1_sig=ttest2(ep1_FR(idxC.(ar1)),ep1_FR(idxI.(ar1)));
                    unit(u).ep2_sig=ttest2(ep2_FR(idxC.(ar2)),ep2_FR(idxI.(ar2)));
                end
                
                figure_handle=figure;
                hold on
                RM_totest=abs([[unit.FRCi_2]-[unit.FRIi_2];[unit.FRC_2]-[unit.FRI_2];[unit.FRCc_2]-[unit.FRIc_2]])';                
                [p_RManova] = anova_rm(abs(RM_totest), 'off');
                
                cols={'r','b','g',[0.5 0.5 0.5],'k'};
                c=0;
                idx=[unit.ep1_sig]==1 & [unit.ep2_sig]==1; c=c+1;
                [N(c),P(c),R(c),pRank(c),D(c)]=scatter_per_sig(unit,idx,1,'o','filled',cols{c});
                idx=[unit.ep1_sig]~=1 & [unit.ep2_sig]==1; c=c+1;
                [N(c),P(c),R(c),pRank(c),D(c)]=scatter_per_sig(unit,idx,1,'o','filled',cols{c});
                idx=[unit.ep1_sig]==1 & [unit.ep2_sig]~=1; c=c+1;
                [N(c),P(c),R(c),pRank(c),D(c)]=scatter_per_sig(unit,idx,1,'o','filled',cols{c});
                idx=[unit.ep1_sig]~=1 & [unit.ep2_sig]~=1; c=c+1;
                [N(c),P(c),R(c),pRank(c),D(c)]=scatter_per_sig(unit,idx,1,'o','MarkerEdgeColor',cols{c});
                idx=true(size(unit)); c=c+1;
                [N(c),P(c),R(c),pRank(c),D(c)]=scatter_per_sig(unit,idx,0,'o');
                ylabel([ep2 ' ' ar2 ' preference Contra-Ipsi' keys.FR_or_idx_central_gaze])
                xlabel([ep1 ' ' ar1 ' preference Contra-Ipsi' keys.FR_or_idx_central_gaze])
                title(['Central gaze only, RM for all gaze positions' num2str(p_RManova(1))])
                
                
                axes_limits=[min([get(gca,'xlim') get(gca,'ylim')]) max([get(gca,'xlim') get(gca,'ylim')])];
                set(gca,'xlim',axes_limits,'ylim',axes_limits);
                axis square;
                
                y_lim=ylim;
                x_lim=xlim;
                
                x_separation=diff(x_lim)/10;
                y_separation=diff(y_lim)/20;
                
                text(x_lim(1)+x_separation,y_lim(2)-1*y_separation,'N =');
                text(x_lim(1)+x_separation,y_lim(2)-2*y_separation,'R =');
                text(x_lim(1)+x_separation,y_lim(2)-3*y_separation,'P =');
                text(x_lim(1)+x_separation,y_lim(2)-4*y_separation,'PRanksAbs =');
                
                D_labels={'x','-','y'};
                for cc=1:c
                    tt=text(x_lim(1)+(cc+1)*x_separation,y_lim(2)-1*y_separation,num2str(N(cc)));                      set(tt,'color',cols{cc});
                    tt=text(x_lim(1)+(cc+1)*x_separation,y_lim(2)-2*y_separation,num2str(round(R(cc)*100)/100));       set(tt,'color',cols{cc});
                    tt=text(x_lim(1)+(cc+1)*x_separation,y_lim(2)-3*y_separation,num2str(round(P(cc)*1000)/1000));     set(tt,'color',cols{cc});
                    tt=text(x_lim(1)+(cc+1)*x_separation,y_lim(2)-4*y_separation,num2str(round(pRank(cc)*1000)/1000)); set(tt,'color',cols{cc});
                    tt=text(x_lim(1)+(cc+1)*x_separation,y_lim(2)-5*y_separation,D_labels{D(cc)+2}); set(tt,'color',cols{cc});
                    
                end
                
                wanted_size=[50 30];
                set(figure_handle, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
                
                
                export_fig([keys.basepath_to_save keys.project_version filesep 'scatter' filesep 'Central gaze only ' ep1 '_' ar1 '_vs_' ep2 '_' ar2 ', ' keys.FR_or_idx_central_gaze], '-pdf','-transparent') % pdf by run
                close(gcf);
            end
        end
        end
    end
end
end

function [N,P,R,P2,D]=scatter_per_sig(varargin)
unit=varargin{1};
idx=varargin{2};
toplot=varargin{3};
N=sum(idx);
if N>2
[R, P] = corr([unit(idx).FRC_1]'-[unit(idx).FRI_1]',[unit(idx).FRC_2]'-[unit(idx).FRI_2]');
P2 = signrank(abs([unit(idx).FRC_1]'-[unit(idx).FRI_1]'),abs([unit(idx).FRC_2]'-[unit(idx).FRI_2]'));
D = sign(nanmean(abs([unit(idx).FRC_2]'-[unit(idx).FRI_2]')) - nanmean(abs([unit(idx).FRC_1]'-[unit(idx).FRI_1]')));
else
   R=0;
   P=1;
   P2=1;
   D=0;
end



if toplot
scatter([unit(idx).FRC_1]-[unit(idx).FRI_1],[unit(idx).FRC_2]-[unit(idx).FRI_2],varargin{4:end});
end
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
