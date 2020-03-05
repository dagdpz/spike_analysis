function ph_gaze(population,modified_keys)
warning('off','MATLAB:catenate:DimensionMismatch');

for fn=fieldnames(modified_keys)'
    keys.(fn{:})=modified_keys.(fn{:});
end
% keys.path_to_save=[keys.drive filesep keys.basepath_to_save filesep keys.project_version filesep 'gaze_analysis' filesep];
% if ~exist(keys.path_to_save,'dir')
%     mkdir([keys.drive filesep keys.basepath_to_save filesep keys.project_version ], 'gaze_analysis');
% end

%% tuning table preparation and grouping
[tuning_per_unit_table]                 = ph_load_extended_tuning_table(keys);
[tuning_per_unit_table, Sel_for_title]  = ph_reduce_tuning_table(tuning_per_unit_table,keys);
if keys.GA.FR_subtract_baseline
    Sel_for_title =[Sel_for_title,{'base';'=';keys.GA.epoch_BL;', '}];
end
idx_group_parameter=DAG_find_column_index(tuning_per_unit_table,keys.GA.group_parameter);
group_values=tuning_per_unit_table(:,idx_group_parameter);
group_values=cellfun(@num2str, group_values, 'UniformOutput', false);
cell_in_any_group=[false; ~ismember(group_values(2:end),keys.GA.group_excluded)];
unique_group_values=unique(group_values(cell_in_any_group));
if isempty(unique_group_values)
    disp('no relevant groups found');
    return;
end
complete_unit_list={population.unit_ID}';
%
idx_unitID=DAG_find_column_index(tuning_per_unit_table,'unit_ID');

%% define conditions to look at
all_trialz=[population.trial];

condition_parameters  ={'reach_hand','choice'};
per_trial.types       =[all_trialz.type];
per_trial.effectors   =[all_trialz.effector];
per_trial.hands       =[all_trialz.reach_hand];
per_trial.choice      =[all_trialz.choice];

u_types     =unique(per_trial.types);
u_effectors =unique(per_trial.effectors);
u_hands     =unique(per_trial.hands);
u_choice    =unique(per_trial.choice);

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
u_types     =unique(type_effectors(:,1));
u_effectors =unique(type_effectors(:,2));

condition_matrix    = combvec(u_hands,u_choice)';


%% Convert to ipsi/contra, Baseline subtraction, normalization, re-ordering, and gaussian RFs
a=1;
%% find positions
tr_con=ismember([all_trialz.completed],keys.cal.completed);
[whatisthis]=ph_arrange_positions_and_plots(all_trialz(tr_con),keys);
allpos=vertcat(whatisthis.trial.position);

clear whatisthis

for xory=keys.GA.xory
    clear condition norm_factor
    xorytag=xory{:};
    
    switch xorytag
        case 'x'
            xy_index=1;
        case 'y'
            xy_index=2;
    end
    positions_temp=unique(allpos(:,xy_index));
    pos_temp_idx=true(size(positions_temp,1),1);
    for x=1:size(positions_temp,1)-1
        if any(all(abs(bsxfun(@minus,positions_temp(x+1:end,:),positions_temp(x,:)))<4,2)) %% precision....
            pos_temp_idx(x)=false;
        end
    end
    positions=positions_temp(pos_temp_idx,:);
    
    
    condition_matrix_p    = combvec(condition_matrix',positions')';
    %% type? (+effector?)
    for t=1:size(type_effectors,1) %typ=unique(per_trial.types)
        typ=type_effectors(t,1);
        eff=type_effectors(t,2);
        tya=find(u_types==typ)+a*(numel(u_types))-1;
        idx_existing=DAG_find_column_index(tuning_per_unit_table,['existing_' type_effector_short{t} '_' keys.arrangement(1:3)]);
        tya_existing{tya}= [true; ismember(vertcat(tuning_per_unit_table{2:end,idx_existing}),true)];
        
        
        keys=ph_get_epoch_keys(keys,typ,eff,sum(type_effectors(:,1)==typ)>1);
        normalization_epoch     =find(ismember(keys.EPOCHS(:,1),keys.GA.epoch_for_normalization));
        
        for g=1:numel(unique_group_values)
            unitidx=ismember(complete_unit_list,tuning_per_unit_table(ismember(group_values,unique_group_values(g))&tya_existing{tya},idx_unitID));
            units=find(all(unitidx,2))';
            for u=units
                tr_con=ismember([population(u).trial.completed],keys.cal.completed) & [population(u).trial.accepted];
                [pop]=ph_arrange_positions_and_plots(population(u).trial(tr_con),keys);
                pop=ph_LR_to_CI(pop,population(u).target);
                poppos=vertcat(pop.trial.position);
                %% normalization factor
                for c=1:size(condition_matrix_p,1) %% position condition matrix...!
                    clear trpar
                    switch keys.GA.normalization
                        case 'by_condition'
                            condition_parameters_tmp=[condition_parameters {'hemifield'}];
                            for par=1:numel(condition_parameters_tmp)
                                fn=condition_parameters_tmp{par};
                                trpar(par,:)=[pop.trial.(fn)]==condition_matrix(c,par);
                            end
                            trpar(end+1,:)=[pop.trial.type]==typ & [pop.trial.effector]==eff;
                            
                        case 'by_position'
                            condition_parameters_tmp=condition_parameters;
                            for par=1:numel(condition_parameters_tmp)
                                fn=condition_parameters_tmp{par};
                                trpar(par,:)=[pop.trial.(fn)]==condition_matrix_p(c,par);
                            end
                            
                            trpar(end+1,:)=all(abs(bsxfun(@minus,poppos(:,xy_index),condition_matrix_p(c,par+1)))<4,2);
                            %trpar(end+1,:)=ismember(vertcat(pop.trial.position),condition_matrix_p(c,par+1:par+2),'rows');
                            trpar(end+1,:)=[pop.trial.type]==typ & [pop.trial.effector]==eff & [pop.trial.accepted];
                            
                        case 'by_effector'
                            trpar=[pop.trial.type]==typ & [pop.trial.effector]==eff;
                        case 'by_type'
                            trpar=[pop.trial.type]==typ;
                        case 'by_all_trials'
                            trpar=true(size([pop.trial]));
                        case 'none'
                            norm_factor(u,c)= 1;
                            continue
                    end
                    
                    tr_con=all(trpar,1);
                    per_epoch=vertcat(pop.trial(tr_con).epoch);
                    
                    if strcmp(keys.GA.normalization,'none') || isempty(normalization_epoch)
                        norm_factor(u,c)= 1;
                    elseif isempty(per_epoch)
                        norm_factor(u,c)=NaN;
                    elseif ~isempty(normalization_epoch)
                        norm_factor(u,c)=nanmean([per_epoch(:,normalization_epoch).FR NaN]);
                    end
                end
                %not there yet, maximum??
                norm_factor(u,:)=deal(max([norm_factor(u,:) 0]));
                norm_factor(norm_factor==0)=1;
                
                %% average FR per unit
                for e=1%:size(condition_matrix,1)/2 ??
                    clear trpar
                    c=1;
                    %                    c=find(all(bsxfun(@minus,conditions_out,[eff,condition_matrix(e,1),condition_matrix(e,2)])==0,2));
                    %
                    %                     for par=1:numel(condition_parameters)
                    %                         fn=condition_parameters{par};
                    %                         trpar(par,:)=[pop.trial.(fn)]==condition_matrix(ep,par); %condition_matrix_wohf?
                    %                         %trpar(par,:)=arrayfun(@(x) any([x.(fn)]==condition_matrix(c,par)), pop.trial);
                    %                     end
                    %                     trpar(end+1,:)=[pop.trial.type]==typ & [pop.trial.effector]==eff;
                    %                     tr_con=all(trpar,1);
                    %                     per_epoch=vertcat(pop.trial(tr_con).epoch);
                    
                    
                    for p=1:size(positions,1)
                        %% supposedly matches with target position precision...
                        tr=all(abs(bsxfun(@minus,poppos(:,xy_index),positions(p,:)))<4,2) & [pop.trial.accepted]'; %% precision
                        
                        per_epoch=vertcat(pop.trial(tr).epoch);
                        for ep=1:size(keys.EPOCHS,1)
                            sta=keys.EPOCHS{ep,2};
                            condition(tya,c).epoch(ep).unit(u).per_position(p)=...
                                nanmean([per_epoch(:,ep).FR])/norm_factor(u,c);
                        end
                    end
                    
                    
                end
            end
        end
    end
    
    %% plots
    for t=1:size(condition,1)
        typ=u_types(mod(t-1,numel(u_types))+1);
        current=[condition(t,:)];
        %% effector missing in figure title
        fig_title=sprintf('%s %s %s hnd %s ch %s %s normalized %s in %s grouped by %s ',...
            keys.monkey,[keys.conditions_to_plot{:}],keys.arrangement,     mat2str(u_hands),mat2str(double(u_choice)),[Sel_for_title{:}],keys.GA.normalization,keys.GA.epoch_for_normalization,keys.GA.group_parameter);
        filename=sprintf('%s %s %s hnd %s ch %s %s N_%s %s %s ',...
            keys.monkey,[keys.conditions_to_plot{:}],keys.arrangement(1:3),mat2str(u_hands),mat2str(double(u_choice)),[Sel_for_title{:}],keys.GA.normalization,keys.GA.epoch_for_normalization,keys.GA.group_parameter);
        
        for eff=u_effectors
            
            plot_1_title            = [fig_title  ' FR'];
            FR_figure_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_1_title);
            keys=ph_get_epoch_keys(keys,typ,eff,sum(type_effectors(:,1)==typ)>1);
            clear n_max
            for g=1:numel(unique_group_values)
                unitidx=ismember(complete_unit_list,tuning_per_unit_table(ismember(group_values,unique_group_values(g)) & tya_existing{t},idx_unitID));
                group_units=find(all(unitidx,2))';
                
                for ep=1:size(keys.EPOCHS,1)
                    sp_handle(g,ep)=subplot(numel(unique_group_values),size(keys.EPOCHS,1),(g-1)*size(keys.EPOCHS,1)+ep);
                    AA=vertcat(current.epoch(ep).unit(group_units).per_position)';
                    x_axis=repmat(1:size(AA,1),size(AA,2),1)';
                    [~, max_pos]=max(AA,[],1);
                    if keys.GA.center_at_max
                        x_axis=x_axis-repmat(max_pos,size(x_axis,1),1)+(size(AA,1));
                    end
                    plot(x_axis,AA)
                    
                    x_ticks=min(min(x_axis)):max(max(x_axis));
                    n_per_pos=hist(x_axis(:),x_ticks);
                    n_max(ep,g).per_pos=hist(max_pos,1:size(AA,1));
                    title(keys.EPOCHS{ep,1})
                    x_lim=get(gca,'xlim');
                    if ep==1
                        text(x_lim(1)-diff(x_lim)/2,1,[unique_group_values(g); {['N=' num2str(numel(group_units))]}],'rotation',90);
                    end
                    set(gca,'ylim',[0 1],'xlim',[min(x_ticks) max([x_ticks,2])],'xtick',x_ticks,'xticklabel',n_per_pos); %% consider non normalized?
                end
            end
            ph_title_and_save(FR_figure_handle,  [filename ' FR'],plot_1_title,keys)
            
            % histogram
            plot_2_title            = [fig_title  ' histogram'];
            Histogram_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_2_title);
            for ep=1:size(keys.EPOCHS,1)
                subplot(1,size(keys.EPOCHS,1),ep);
                bar(vertcat(n_max(ep,:).per_pos)','stacked');
                if ep==1
                    legend(unique_group_values)
                end
                title(keys.EPOCHS{ep,1})
            end
            ph_title_and_save(Histogram_handle,  [filename ' histogram'],plot_2_title,keys)
        end
    end
    
end
end
% 
% function title_and_save(figure_handle,filename,plot_title,keys)
% mtit(figure_handle,  plot_title, 'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none');
% stampit;
% wanted_size=[50 30];
% set(figure_handle, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
% if keys.plot.export
%     export_fig([keys.path_to_save, filename], '-pdf','-transparent') % pdf by run
%     close all
% end
% end
