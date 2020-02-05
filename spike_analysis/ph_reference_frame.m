function epoch=ph_reference_frame(population,modified_keys)

warning('off','MATLAB:catenate:DimensionMismatch');

n_iterations=1000;

for fn=fieldnames(modified_keys)'
    keys.(fn{:})=modified_keys.(fn{:});
end
% keys.path_to_save=[keys.drive filesep keys.basepath_to_save filesep keys.project_version filesep 'scatter' filesep];
% if ~exist(keys.path_to_save,'dir')
%     mkdir([keys.drive filesep keys.basepath_to_save filesep keys.project_version ], 'scatter');
% end

%% tuning table preparation and grouping
[tuning_per_unit_table]                 = ph_load_extended_tuning_table(keys);
[tuning_per_unit_table, Sel_for_title]  = ph_reduce_tuning_table(tuning_per_unit_table,keys);
if keys.RE.FR_subtract_baseline
    Sel_for_title =[Sel_for_title,{'base';'=';keys.RE.epoch_BL;', '}];
end
idx_group_parameter=find_column_index(tuning_per_unit_table,keys.RE.group_parameter);
group_values=tuning_per_unit_table(:,idx_group_parameter);
group_values=cellfun(@num2str, group_values, 'UniformOutput', false);
cell_in_any_group=[false; ~ismember(group_values(2:end),keys.RE.group_excluded)];
unique_group_values=unique(group_values(cell_in_any_group));
if isempty(unique_group_values)
    disp('no relevant groups found');
    return;
end
complete_unit_list={population.unit_ID}';
idx_unitID=find_column_index(tuning_per_unit_table,'unit_ID');

all_IDs={population.unit_ID};
valid_IDs=ismember(all_IDs,tuning_per_unit_table(:,1));
population=population(valid_IDs);

%% define conditions to look at
all_trialz=[population.trial];

%condition_parameters  ={'reach_hand','choice','perturbation'};
per_trial.types       =[all_trialz.type];
per_trial.effectors   =[all_trialz.effector];
per_trial.hands       =[all_trialz.reach_hand];
per_trial.choice      =[all_trialz.choice];
per_trial.perturbation=[all_trialz.perturbation];
per_trial.perturbation(ismember(per_trial.perturbation, keys.cal.perturbation_groups{1}))=0;
per_trial.perturbation(ismember(per_trial.perturbation, keys.cal.perturbation_groups{2}))=1;

u_hemifields=[-1,0,1]; % why does this have to be hardcoded? ---> Because case not defined yet, case defines positions !!

u_types     =unique(per_trial.types);
u_effectors =unique(per_trial.effectors);
u_hands     =unique(per_trial.hands);
u_choice    =unique(per_trial.choice);
u_perturbation    =unique(per_trial.perturbation);
u_perturbation=u_perturbation(~isnan(u_perturbation));

if any(u_perturbation==1) % temporary, better solution
    u_hemifields=[-1,1]; % why does this have to be hardcoded? And what do we do with vertical targets?
end
%% limit conditions key?
u_hands     =u_hands(ismember(u_hands,keys.tt.hands));
u_choice    =u_choice(ismember(u_choice,keys.tt.choices));

all_type_effectors      = combvec(u_types,u_effectors)';
type_effectors =[];

% redifine type_effectors to include only relevant
for t=1:size(all_type_effectors,1)
    typ=all_type_effectors(t,1);
    eff=all_type_effectors(t,2);
    [~, type_effector_short{t}]=get_type_effector_name(typ,eff);
    if ~ismember(type_effector_short{t},keys.conditions_to_plot) %|| sum(tr_con)<1
        continue;
    end
    type_effectors=[type_effectors; all_type_effectors(t,:)];
end
type_effector_short(~ismember(type_effector_short,keys.conditions_to_plot))=[];
u_types     =unique(type_effectors(:,1));
u_effectors =unique(type_effectors(:,2));


% here we go

for t=1:size(type_effectors,1) %still didnt add effectors properly
    typ=type_effectors(t,1);
    eff=type_effectors(t,2);
    
    keys=get_epoch_keys(keys,typ,eff,sum(type_effectors(:,1)==typ)>1);
    %filenamepart=[Sel_for_title{:} type_effector_short{t} 'reference frame correlations'];
    
    fig_title=sprintf('%s %s %s hnd %s ch %s %s normalized %s in %s grouped by %s ',...
        keys.monkey,[keys.conditions_to_plot{:}],keys.arrangement,     mat2str(u_hands),mat2str(double(u_choice)),[Sel_for_title{:}],keys.RE.normalization,keys.RE.epoch_for_normalization,keys.RE.group_parameter);
    filename=sprintf('%s %s %s hnd %s ch %s %s N_%s %s %s ',...
        keys.monkey,[keys.conditions_to_plot{:}],keys.arrangement(1:3),mat2str(u_hands),mat2str(double(u_choice)),[Sel_for_title{:}],keys.RE.normalization,keys.RE.epoch_for_normalization,keys.RE.group_parameter);
    
    for u=1:numel(population)
        keys.arrangement='movement vectors';
        trials=[population(u).trial];
        trials=trials([trials.completed] & [trials.accepted] & [trials.type]==typ);
        tar_positions=[trials.tar_pos];
        gaz_positions=[trials.fix_pos];
        mov_positions=tar_positions-gaz_positions;
        
        u_mov_pos=unique(mov_positions);
        u_tar_pos=unique(tar_positions);
        u_gaze_pos=unique(gaz_positions);
        
        EPOCHS=keys.EPOCHS;
        for ep=1:size(EPOCHS,1)
            clear cor_tar cor_mov Bootstrap FR_mov FR_tar
            for g=1:numel(u_gaze_pos)
                for p=1:numel(u_mov_pos)
                    EPS=vertcat(trials(mov_positions==u_mov_pos(p) & gaz_positions==u_gaze_pos(g)).epoch);
                    EPS_idx=1:size(EPS,1);
                    FR_mov(p,g)=nanmean([EPS(:,ep).FR]);
                    for it=1:n_iterations
                        Bootstrap(it).FR_mov(p,g)=nanmean([EPS(randsample(EPS_idx,round(numel(EPS_idx)/10*8)),ep).FR]);
                    end
                end
                for p=1:numel(u_tar_pos)
                    EPS=vertcat(trials(tar_positions==u_tar_pos(p) & gaz_positions==u_gaze_pos(g)).epoch);
                    EPS_idx=1:size(EPS,1);
                    if ~isempty(EPS)
                        FR_tar(p,g)=nanmean([EPS(:,ep).FR]);
                        for it=1:n_iterations
                            Bootstrap(it).FR_tar(p,g)=nanmean([EPS(randsample(EPS_idx,round(numel(EPS_idx)/10*8)),ep).FR]);
                        end
                    else
                        FR_tar(p,g)=NaN;
                        for it=1:n_iterations
                            Bootstrap(it).FR_tar(p,g)=NaN;
                        end
                    end
                end
            end
            
            
            for g=1:numel(u_gaze_pos)
                h=mod(g,numel(u_gaze_pos))+1;
                idx_valid=~isnan(FR_tar(:,g)) & ~isnan(FR_tar(:,h));
                FR_temp1=FR_tar(idx_valid,g);
                FR_temp2=FR_tar(idx_valid,h);
                cor_tar(g)=corr(FR_temp1,FR_temp2);
                
                idx_valid=~isnan(FR_mov(:,g)) & ~isnan(FR_mov(:,h));
                FR_temp1=FR_mov(idx_valid,g);
                FR_temp2=FR_mov(idx_valid,h);
                cor_mov(g)=corr(FR_temp1,FR_temp2);
            end
            
            epoch(ep).av_cor_mov(u)=nanmean(cor_mov);
            epoch(ep).av_cor_tar(u)=nanmean(cor_tar);
            
            
            
            for it=1:n_iterations
                for g=1:numel(u_gaze_pos)
                    h=mod(g,numel(u_gaze_pos))+1;
                    idx_valid=~isnan(Bootstrap(it).FR_tar(:,g)) & ~isnan(Bootstrap(it).FR_tar(:,h));
                    FR_temp1=Bootstrap(it).FR_tar(idx_valid,g);
                    FR_temp2=Bootstrap(it).FR_tar(idx_valid,h);
                    cor_tar(g)=corr(FR_temp1,FR_temp2);
                    
                    idx_valid=~isnan(Bootstrap(it).FR_mov(:,g)) & ~isnan(Bootstrap(it).FR_mov(:,h));
                    FR_temp1=Bootstrap(it).FR_mov(idx_valid,g);
                    FR_temp2=Bootstrap(it).FR_mov(idx_valid,h);
                    cor_mov(g)=corr(FR_temp1,FR_temp2);
                end
                
                Bootstrap(it).av_cor_mov=nanmean(cor_mov);
                Bootstrap(it).av_cor_tar=nanmean(cor_tar);
            end
            
            epoch(ep).BT_CI_mov(u,:)=prctile([Bootstrap.av_cor_mov],[2.5 97.5]);
            epoch(ep).BT_CI_tar(u,:)=prctile([Bootstrap.av_cor_tar],[2.5 97.5]);
            epoch(ep).BT_cor_mov(u)=nanmean([Bootstrap.av_cor_mov]);
            epoch(ep).BT_cor_tar(u)=nanmean([Bootstrap.av_cor_tar]);
            
            if epoch(ep).av_cor_mov(u)>epoch(ep).BT_CI_tar(u,2) && epoch(ep).av_cor_tar(u)<epoch(ep).BT_CI_mov(u,1) && epoch(ep).BT_CI_mov(u,1)>0 %% mov larger
                epoch(ep).BT_CI_sig(u)=1;
            elseif epoch(ep).av_cor_tar(u)>epoch(ep).BT_CI_mov(u,2) && epoch(ep).av_cor_mov(u)<epoch(ep).BT_CI_tar(u,1) && epoch(ep).BT_CI_tar(u,1)>0 %% tar larger
                epoch(ep).BT_CI_sig(u)=2;
            else
                epoch(ep).BT_CI_sig(u)=0;
            end
        end
    end
    
    for only_ANOVA=0:1        
        figure_handle              = figure('units','normalized','outerposition',[0 0 1 1],'name',[fig_title type_effector_short{t} 'reference frame correlations']);
        for ep=1:size(EPOCHS,1)
            epoc=EPOCHS{ep,1};
            col1=strfind(tuning_per_unit_table(:,find_column_index(tuning_per_unit_table,['in_AH_' epoc '_position_Msac_mov'])),'true');
            col2=strfind(tuning_per_unit_table(:,find_column_index(tuning_per_unit_table,['in_AH_' epoc '_PxF_Msac_mov'])),'true');
            col3=strfind(tuning_per_unit_table(:,find_column_index(tuning_per_unit_table,['in_AH_' epoc '_position_Msac_tar'])),'true');
            col4=strfind(tuning_per_unit_table(:,find_column_index(tuning_per_unit_table,['in_AH_' epoc '_PxF_Msac_tar'])),'true');
            ANOVA_units=tuning_per_unit_table(cellfun(@(w,x,y,z) any(w)||any(x)||any(y)||any(z),col1,col2,col3,col4),1);
            
            subplot(ceil(sqrt(size(EPOCHS,1))),ceil(sqrt(size(EPOCHS,1))),ep)
            hold on
            plot([-1 1],[-1 1],'k');
            plot([0 0],[-1 0],'k:');
            plot([-1 0],[0 0],'k:');
            for sig=0:2
                units=find([epoch(ep).BT_CI_sig]==sig);
                for u=units
                    switch epoch(ep).BT_CI_sig(u)
                        case 0
                            col=[0.5 0.5 0.5];
                        case 1
                            col='r';
                        case 2
                            col='g';
                    end
                    if ismember(population(u).unit_ID,ANOVA_units)
                        if epoch(ep).BT_CI_sig(u)>0
                            plot( epoch(ep).BT_CI_mov(u,:),[epoch(ep).av_cor_tar(u) epoch(ep).av_cor_tar(u)],'Color',col,'linewidth',0.001);
                            plot([epoch(ep).av_cor_mov(u) epoch(ep).av_cor_mov(u)],epoch(ep).BT_CI_tar(u,:),'Color',col,'linewidth',0.001);
                        end
                        scatter(epoch(ep).av_cor_mov(u),epoch(ep).av_cor_tar(u),10,col,'filled');
                    elseif ~only_ANOVA
                        if epoch(ep).BT_CI_sig(u)>0
                            plot( epoch(ep).BT_CI_mov(u,:),[epoch(ep).av_cor_tar(u) epoch(ep).av_cor_tar(u)],'Color',col,'linewidth',0.001);
                            plot([epoch(ep).av_cor_mov(u) epoch(ep).av_cor_mov(u)],epoch(ep).BT_CI_tar(u,:),'Color',col,'linewidth',0.001);
                        end
                        scatter(epoch(ep).av_cor_mov(u),epoch(ep).av_cor_tar(u),10,col);
                    end
                end
                UID={population(units).unit_ID};
                epoch(ep).unit_IDs{sig+1,1}=UID;
                epoch(ep).unit_IDs{sig+1,2}=UID(ismember(UID,ANOVA_units));
            end
            if only_ANOVA
                anova_member_index=ismember({population.unit_ID},ANOVA_units);
                filename_suffix='only_ANOVA';
            else
                anova_member_index=true(1,numel(population))  ;
                filename_suffix='all';
            end
            text(-0.8,0.8,['N=' num2str(sum([epoch(ep).BT_CI_sig]==1 & anova_member_index))],'color','r');
            text(-0.8,0.7,['N=' num2str(sum([epoch(ep).BT_CI_sig]==2 & anova_member_index))],'color','g');
            text(-0.8,0.6,['N=' num2str(sum([epoch(ep).BT_CI_sig]==0 & anova_member_index))],'color',[0.5 0.5 0.5]);
            %                   [~,p]=ttest([epoch(ep).av_cor_mov],[epoch(ep).av_cor_tar]);
            %                    if mean([epoch(ep).av_cor_mov]-[epoch(ep).av_cor_tar])<0
            %                        largerone='screen larger, ';
            %                    else
            %                        largerone='ret larger, ';
            %                    end
            [~,p]=ttest([epoch(ep).BT_cor_mov(anova_member_index)],[epoch(ep).BT_cor_tar(anova_member_index)]);
            if mean([epoch(ep).BT_cor_mov(anova_member_index)]-[epoch(ep).BT_cor_tar(anova_member_index)])<0
                largerone='screen larger, ';
            else
                largerone='ret larger, ';
            end
            
            
            xlabel('retinocentric');
            ylabel('screen centered');
            title([EPOCHS{ep,1} ', ' largerone 'p= ' num2str(p)]);
            xlim([-1.001 1.001]);
            ylim([-1.001 1.001]);
            axis square
        end
        
        ph_title_and_save(figure_handle,[filename ' ' filename_suffix ', ' type_effector_short{t}],[filename ' ' filename_suffix ', ' type_effector_short{t}],keys)
        
%         wanted_size=[50 30];
%         set(figure_handle, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
%         export_fig([keys.drive filesep keys.basepath_to_save filesep keys.project_version filesep 'scatter' filesep filename ' ' filename_suffix ', ' type_effector_short{t}], '-pdf','-transparent') % pdf by run
%         close(gcf);
    end
    
    save([keys.path_to_save filesep keys.monkey '_' type_effector_short{t} '_' keys.arrangement '_correlation coefficients reference frame'],'epoch');
end
end
%
% end
% end
% end
% end

