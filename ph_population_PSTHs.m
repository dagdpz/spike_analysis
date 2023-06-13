function ph_population_PSTHs(population,modified_keys)
warning('off','MATLAB:catenate:DimensionMismatch');

for fn=fieldnames(modified_keys)'
    keys.(fn{:})=modified_keys.(fn{:});
end

%% check if this is set
% keys.normalization_field='PO';

keys.PO.FR_subtract_baseline=~strcmp(keys.PO.epoch_BL,'none'); %% this one we should not need here any more !?
[TT,idx,group_values,unique_group_values]=ph_readout_tuning_table(keys);

if isempty(unique_group_values)
    return;
end

complete_unit_list={population.unit_ID}';
[unit_valid,TM]=ismember(complete_unit_list,TT(:,idx.unitID));
population=population(unit_valid);
complete_unit_list={population.unit_ID}';
population_group=group_values(TM(unit_valid));
all_trialz=[population.trial];
[UC, CM, labels]=ph_get_condition_matrix(all_trialz,keys);

%% fix labels --> and with labels colors!
legend_labels_hem={};
legend_labels_prf={};
legend_labels_pos={};
for h=UC.hemifield % append hemifield labels, careful with the order!
    legend_labels_hem=[legend_labels_hem; strcat(labels,['_' keys.labels.hemifield{h+2}])];
end
for h=1:2          % append preference labels, careful with the order!
    legend_labels_prf=[legend_labels_prf; strcat(labels,['_' keys.labels.preferred{h}])];
end
for p=1:size(UC.position,1)
    legend_labels_pos=[legend_labels_pos; labels];
end
% repeat for all effectors
legend_labels_hem=repmat(reshape(legend_labels_hem,numel(legend_labels_hem),1),numel(UC.effector),1);
legend_labels_prf=repmat(reshape(legend_labels_prf,numel(legend_labels_prf),1),numel(UC.effector),1);
legend_labels_pos=repmat(reshape(legend_labels_pos,numel(legend_labels_pos),1),numel(UC.effector),1);

conditions_out              = combvec(UC.effector,CM')';
conditions_hf               = combvec(UC.hemifield,conditions_out')';
conditions_pref             = combvec([0 1],conditions_out')';

if isfield(UC,'reach_hand') && any(UC.reach_hand~=0) && any(UC.perturbation==1) %splitting to all 4 hand space conditions if hands are involved
    [~,~,columns_hf] = unique(conditions_hf(:,[1,3]),'rows');
    [~,~,columns_pref] = unique(conditions_pref(:,[1,3]),'rows');
else
    columns_hf          = ones(size(conditions_hf,1),1);
    columns_pref        = ones(size(conditions_pref,1),1);
end

% if any(UC.perturbation==1) % temporary, better solution
%     if strcmp(keys.arrangement,'hands_inactivation_in_ch')
%         [~,~,columns_hf] = unique(conditions_hf(:,[1,3]),'rows');
%     end
%     condition_matrix(condition_matrix(:,1)==1 & condition_matrix(:,2)==1 & condition_matrix(:,4)==1,:)=[];
%     condition_matrix(condition_matrix(:,1)==2 & condition_matrix(:,2)==1 & condition_matrix(:,4)==-1,:)=[];
%     conditions_hf(conditions_hf(:,1)==1 & conditions_hf(:,3)==1 & conditions_hf(:,4)==1,:)=[];
%     conditions_hf(conditions_hf(:,1)==-1 & conditions_hf(:,3)==2 & conditions_hf(:,4)==1,:)=[];
% else
%     columns_hf          = ones(size(conditions_hf,1),1);
%     columns_pref        = ones(size(conditions_pref,1),1);
% end

%% normalization
[~, condition,~,pref_valid]=ph_condition_normalization(population,keys,UC,CM);

%% condition comparison (kinda there to be used, but not used currently ???)
comparisons_per_effector(1).reach_hand{1}=0;
comparisons_per_effector(1).reach_hand{2}=0;
comparisons_per_effector(1).hemifield{1}=[-1];
comparisons_per_effector(1).hemifield{2}=[1];
comparisons_per_effector(1).choice{1}=0;
comparisons_per_effector(1).choice{2}=0;
comparisons_per_effector(1).color=[1 0 0];
condition_parameters_comparissons = [{'hemifield'} {'effector'} keys.condition_parameters];

%% population significance ?
for t=1:size(condition,1)
    current=[condition(t,:).per_hemifield];
    current_window=vertcat(current.window);
    for w=1:numel(current(1).window) %% epochs might be different for different effectors.......... !!!!====???
        for g=1:numel(unique_group_values)
            %%KK
            %unitidx=ismember(complete_unit_list,TT(ismember(group_values,unique_group_values(g)),idx.unitID));
            units=find(ismember(population_group,unique_group_values(g)));
            
            for comp=1:numel(comparisons_per_effector)
                cM1=true(size(conditions_hf));
                cM2=true(size(conditions_hf));
                for k=1:numel(condition_parameters_comparissons)
                    if isfield(comparisons_per_effector,condition_parameters_comparissons{k}) &&...
                            ~ isempty(comparisons_per_effector(comp).(condition_parameters_comparissons{k}))
                        cM1(~ismember(conditions_hf(:,k),comparisons_per_effector(comp).(condition_parameters_comparissons{k}){1}),k)=false;
                        cM2(~ismember(conditions_hf(:,k),comparisons_per_effector(comp).(condition_parameters_comparissons{k}){2}),k)=false;
                    end
                end
                n=0;
                for eff=UC.effector
                    n=n+1;
                    cM1(conditions_hf(:,2)~=eff,2)=false;
                    cM2(conditions_hf(:,2)~=eff,2)=false;
                    c1=find(all(cM1,2));
                    c2=find(all(cM2,2));
                    sigbins(t).group(g).per_effector(n).window(w).bins(comp,:)=ph_compare_by_bin(current_window(:,w),c1,c2,units);
                end
            end
        end
    end
end

%% plots
for t=1:size(condition,1)
    typ=UC.type(mod(t-1,numel(UC.type))+1);
    %save metadata
    [fig_title,filename]=ph_figure_titles(keys);
    unit_IDs=complete_unit_list;
    save([keys.basepath_to_save, keys.project_version, filesep, 'population_meta_data', filesep, 'population_analysis', filesep, filename], 'keys','TT','unit_IDs');
    
    %% PSTH plot
    current=[condition(t,:).per_hemifield];
    legend_labels=legend_labels_hem;
    plot_title_part=' PSTHs';
    units_valid=ones(size(complete_unit_list,1),1);
    column_indexes=columns_hf;
    plot_PSTH
    
    %% PSTH preferred and unpreferred plot
    current=[condition(t,:).per_preference];
    legend_labels=legend_labels_prf;
    plot_title_part=[' pref in ' keys.PO.epoch_PF ' PSTHs'];
    units_valid=pref_valid(t,:)';
    column_indexes=columns_pref;
    if ~any(UC.perturbation==1)
        plot_PSTH
    end
    
    unique_group_values_tmp=unique_group_values;
    for gt=1:numel(unique_group_values_tmp)
        unique_group_values=unique_group_values_tmp(gt);
        if   keys.PO.plot_per_position
            %% PSTH per position plot
            current=[condition(t,:).per_position]; current=current(:);
            legend_labels=legend_labels_pos;
            plot_title_part=['=' unique_group_values{1} ' PSTHs per position'];
            units_valid=ones(size(complete_unit_list,1),1);
            column_indexes=columns_pref;
            plot_PSTH
            
            %             %% unfortunately, these are changed withing plot_PSTH, but should be the same as here (?)
            %             unitidx=ismember(complete_unit_list,TT(ismember(group_values,unique_group_values_tmp(1)),idx.unitID));
            %             group_units=find(all(unitidx,2))';
            
            %
            %             %% PSTH per position plot (by initial fixation - these here need to be fixed!)
            %             if true
            %
            %                 current=[condition(t,:).per_position]; current=current(:);
            %                 plot_title_part=['=' unique_group_values{1} ' PSTHs per position ensu ' keys.PO.epoch_PF '-' keys.PO.epoch_BL];
            %                 units_valid=ones(size(complete_unit_list,1),1);
            %                 column_indexes=columns_pref;
            %                 plot_PSTH
            %
            % %                 current=[condition(t,:).per_hemifield]; current=current(:);
            % %                 plot_title_part=['=' unique_group_values{1} ' PSTHs ensu ' keys.PO.epoch_PF '-' keys.PO.epoch_BL];
            % %                 units_valid=ones(size(complete_unit_list,1),1);
            % %                 column_indexes=columns_hf;
            % %                 plot_PSTH
            %
            %             end
            
            
        end
    end
    unique_group_values=unique_group_values_tmp;
end

    function plot_PSTH
        %% one figure for each effector, somehting probably does not work here, because (!!) suplots in each figure
        spf=[];
        for eff=UC.effector 
            ef=find(UC.effector==eff);
            keys=ph_get_epoch_keys(keys,typ,eff,sum(UC.type_effector(:,1)==typ)>1);
            [~, type_effector_short] = MPA_get_type_effector_name(typ,eff);
            plot_title              = [fig_title plot_title_part ', ' type_effector_short ];
            PSTH_summary_handle(ef)     = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',plot_title);
            eff_empty=1;
            for g=1:numel(unique_group_values)
                group_units=find(ismember(population_group,unique_group_values(g)));
                
                current_window=vertcat(current.window);
                current_units=vertcat(current_window(:,1).unit);
                
                empty_conditions_and_units=arrayfun(@(x) isempty(x.average_spike_density) || all(isnan( x.average_spike_density)),current_units) ; %nans not needed hopefully
                %empty_conditions_and_units=arrayfun(@(x) isempty(x.average_spike_density) ,current_units) ; %nans not needed hopefully, still there in MP mods
                empty_conditions_and_units(:,end+1:numel(group_units))=true; % from the last valid to the end
                condition_empty=all(empty_conditions_and_units,2);
                
                if all(condition_empty)
                   continue;
                else
                    eff_empty=0;
                end
                
                if any(strfind(plot_title_part,'per position'))
                    [positions, ~,pos_sub_idx]=unique(vertcat(current.position),'rows');
                    [~, ~,fix_sub_idx]=unique(vertcat(current.fixation),'rows');
                    [subplot_pos, columns_to_loop, rows_to_loop]= DAG_dynamic_positions({positions});
                    conditions_to_loop=find([current.effector]'==eff & ~condition_empty);
                    sp_per_effector=max(subplot_pos);
                else
                    conditions_to_loop=find([current.effector]'==eff & ~condition_empty);
                    [~,columns_to_loop]=ismember(column_indexes,unique(column_indexes(conditions_to_loop)));
                    sp_per_effector=max(columns_to_loop)*numel(unique_group_values);
                end
                n_units_title_part=repmat({''},sp_per_effector*numel(UC.effector),1);
                spl=zeros(sp_per_effector*numel(UC.effector),1);
                
                legend_label_indexes=[];
                legend_line_handles=[];
                
                ensu_units{ef,1}=[];ensu_units{ef,2}=[];ensu_units{ef,3}=[];
                for c=conditions_to_loop(:)'
                    group_condition_units=find(ismember(population_group,unique_group_values(g)) & units_valid & ~empty_conditions_and_units(c,:)')';
                    
                    %group_units=find(ismember(population_group,unique_group_values));
                    current_color=keys.colors.(legend_labels{c})/255;
                    current_line_style=keys.linestyles.(legend_labels{c});
                    column=c; %irrelevant ?
                    
                    if any(strfind(plot_title_part,'per position'))%% still need to fix colors if different conditions are present, AND fix overwriting of different groups.... (how about legends?)
                        spn=subplot_pos(pos_sub_idx(c))+(ef-1)*sp_per_effector;
                        sph(spn)=subplot(rows_to_loop,max(columns_to_loop),spn-(ef-1)*sp_per_effector);
                        signs=[current(c).sign.unit];
                    elseif any(strfind(plot_title_part,'PSTHs'))
                        column=columns_to_loop(c);
                        spn=(g-1)*max(columns_to_loop)+column+(ef-1)*sp_per_effector;
                        sph(spn)=subplot(numel(unique_group_values),max(columns_to_loop),spn-(ef-1)*sp_per_effector);
                    end
                    spl(spn)=spl(spn)+1;
                    spf(spn)=ef;
                    
                    %specific additional loop for enhancement and suppression
                    ensu_loops=1;
                    if any(strfind(plot_title_part,'ensu'))
                        ensu_loops=2;
                        ensucolors={[0 0 1]/spl(spn),[1 0 0]/spl(spn)};
                    end
                    
                    for ensu=1:ensu_loops
                        if any(strfind(plot_title_part,'ensu'))
                            current_color=ensucolors{ensu};
                            units=intersect(group_condition_units,find(signs==(ensu-1.5)*2));
                            if any(ismember(find(signs==1),group_condition_units)) && any(ismember(find(signs==-1),group_condition_units))
                                ensu_units{ef,3}=union(ensu_units{ef,3},units);
                                ensu_units{ef,ensu}=union(ensu_units{ef,ensu},units);
                            end
                        else
                            units=group_condition_units;
                        end
                        hold on
                        if  numel(units)==0
                            continue;
                        end
                        state_shift=0;
                        nonnans_present=0;
                        for w=1:size(keys.PSTH_WINDOWS,1)
                            t_before_state=keys.PSTH_WINDOWS{w,3};
                            t_after_state=keys.PSTH_WINDOWS{w,4};
                            bins=t_before_state:keys.PSTH_binwidth:t_after_state;
                            bins=bins+state_shift-t_before_state;
                            props={'color',current_color,'linewidth',3,'LineStyle',current_line_style};
                            errorbarhandle=shadedErrorBar(bins,nanmean(vertcat(current(c).window(w).unit(units).average_spike_density),1),...
                                sterr(vertcat(current(c).window(w).unit(units).average_spike_density),1),props,1); %% STERR!!!!
                            state_shift=state_shift+t_after_state-t_before_state+0.1;
                            nonnans_present=nonnans_present | any(any(~isnan(vertcat(current(c).window(w).unit(units).average_spike_density))));
                        end
                        if nonnans_present
                            legend_label_indexes=[legend_label_indexes c];
                            n_units_title_part{spn}=[n_units_title_part{spn} '\color[rgb]{' num2str(current_color) '}' num2str(numel(units)) ' '];
                            legend_line_handles=[legend_line_handles errorbarhandle.mainLine];
                        end
                        group_para=keys.PO.group_parameter; group_para(strfind(group_para,'_'))='-';
                        group_val=unique_group_values{g}; group_val(strfind(group_val,'_'))='-';
                        title(sprintf('%s = %s, N =%s ',group_para,unique_group_values{g},n_units_title_part{spn}),'interpreter','tex'); %%
                    end
                end
                y_lim(spn,:)=get(gca,'ylim');
            end
            
            if eff_empty
                continue;
            end
            if keys.plot.population_PSTH_legends
                LL_to_plot=cellfun(@(x) strrep(x,'_',' '),legend_labels,'uniformoutput',false);
                LL_to_plot=LL_to_plot(legend_label_indexes);
                [~,uix]=unique(LL_to_plot);
                legend(legend_line_handles(uix),LL_to_plot(uix),'location','southwest');
            end
        end
        if isempty(spf)
           return;
        end

        %% subplot appearance, and tuning lines
        spf(sph==0)=[];
        sph(sph==0)=[];
        ylimmax=max(max(y_lim));
        ylimmin=min(min(y_lim));
        y_lim_to_plot=[ylimmin ylimmax];
        if ~isempty(keys.PO.y_lim)
            y_lim_to_plot= keys.PO.y_lim;
        end
        
        for spn=1:numel(sph)
            ef=spf(spn);
            axes(sph(spn));
            if ~keys.PO.link_y_lim
                y_lim_to_plot=y_lim(spn,:);
            end
            %% alternative (MP ?)
            %             ef=spf(spn);
            %             set(0, 'CurrentFigure', PSTH_summary_handle(ef));
            %             subplot(sph(spn));
            
            hold on
            
            %% type? completed? choices? hands? (added effector here!!)
            tr=[all_trialz.type]==typ & [all_trialz.effector]==UC.effector(ef) & ismember([all_trialz.completed],keys.cal.completed); % & ismember([all_trialz.reach_hand],UC.reach_hand) & ismember([all_trialz.choice],UC.choice);
            ph_PSTH_background(all_trialz(tr),y_lim_to_plot,y_lim_to_plot,y_lim_to_plot,keys,keys.PO.fontsize_factor)
        end
        for eff=UC.effector
            ef=find(UC.effector==eff);
            set(0, 'CurrentFigure', PSTH_summary_handle(ef));
            axes(sph(find(spf==ef,1)));
            [~, type_effector_short] = MPA_get_type_effector_name(typ,eff);
            plot_title              = [fig_title plot_title_part ', ' type_effector_short ];
            ph_title_and_save(PSTH_summary_handle(ef),  [filename plot_title_part ', ' type_effector_short],plot_title,keys)
        end
        
        %% ensu summary figure
        if any(strfind(plot_title_part,'ensu'))
            for eff=UC.effector
                ef=find(UC.effector==eff);
                [~, type_effector_short] = MPA_get_type_effector_name(typ,eff);
                plot_title              = [fig_title plot_title_part ', ' type_effector_short ];
                PSTH_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'color','w','name',plot_title);
                for g=1:numel(unique_group_values)
                    group_units=find(ismember(population_group,unique_group_values(g)) & units_valid)';
                    
                    %%reducing to only units that have at least one condition
                    current_window=vertcat(current.window);
                    current_units=vertcat(current_window(:,1).unit);
                    empty_conditions_and_units=arrayfun(@(x) isempty(x.average_spike_density),current_units);
                    empty_conditions_and_units(:,end+1:numel(units_valid))=true;
                    condition_empty=all(empty_conditions_and_units,2);
                    conditions_to_loop=find([current.effector]'==eff & ~condition_empty);
                    
                    posneg_idx=0;
                    for c=conditions_to_loop(:)'
                        signs=[current(c).sign.unit];
                        units_pos=intersect(group_units,find(signs==1));
                        units_neg=intersect(group_units,find(signs==-1));
                        if  numel(units_pos)==0 && numel(units_pos)==0 %% not sure if this is enough
                            continue;
                        end
                        posneg_idx=posneg_idx+1;
                        for w=1:size(keys.PSTH_WINDOWS,1)
                            if ~isempty(units_pos) && ~isempty(units_neg)
                                positive(w).psth(posneg_idx,:)=nanmedian(vertcat(current(c).window(w).unit(units_pos).average_spike_density),1);
                                negative(w).psth(posneg_idx,:)=nanmedian(vertcat(current(c).window(w).unit(units_neg).average_spike_density),1);
                                posneg(w).psth(posneg_idx,:)=positive(w).psth(posneg_idx,:)-negative(w).psth(posneg_idx,:);
                            end
                        end
                    end
                    
                    sp(g)=subplot(1,numel(unique_group_values),g);
                    hold on
                    
                    state_shift=0;
                    for w=1:size(keys.PSTH_WINDOWS,1)
                        t_before_state=keys.PSTH_WINDOWS{w,3};
                        t_after_state=keys.PSTH_WINDOWS{w,4};
                        bins=t_before_state:keys.PSTH_binwidth:t_after_state;
                        bins=bins+state_shift-t_before_state;
                        props={'color','r','linewidth',1};
                        if exist('positive','var') && exist('negative','var')
                            errorbarhandle=shadedErrorBar(bins,nanmean(vertcat(positive(w).psth),1),sterr(vertcat(positive(w).psth),1),props,1);
                            props={'color','b','linewidth',1};
                            errorbarhandle=shadedErrorBar(bins,nanmean(vertcat(negative(w).psth),1),sterr(vertcat(negative(w).psth),1),props,1);
                            props={'color','g','linewidth',1};
                            errorbarhandle=shadedErrorBar(bins,nanmean(vertcat(posneg(w).psth),1),sterr(vertcat(posneg(w).psth),1),props,1);
                        end
                        state_shift=state_shift+t_after_state-t_before_state+0.1;
                    end
                    title(['N contributing = ' num2str(numel(ensu_units{ef,2})) ' POS, ' num2str(numel(ensu_units{ef,1})) ' NEG, ' num2str(numel(ensu_units{ef,3})) ' POS-NEG']);
                    %% this part should eventually go in extra loop for same dimensions in each subplot
                    %% also, can we speed up by plotting the same background for every subplot???
                    y_lim_to_plot=get(gca,'ylim');
                    tr=[all_trialz.type]==typ & [all_trialz.effector]==eff & ismember([all_trialz.completed],keys.cal.completed) &...
                        ismember([all_trialz.reach_hand],UC.reach_hand) & ismember([all_trialz.choice],UC.choice);
                    ph_PSTH_background(all_trialz(tr),y_lim_to_plot,y_lim_to_plot,y_lim_to_plot,keys,1)
                end
                ph_title_and_save(PSTH_summary_handle,  [filename plot_title_part ' summary, ' type_effector_short ],plot_title,keys)
            end
        end
    end
end

function h=ph_compare_by_bin(in,c1,c2,units)
% to compare only the same cells (!)
c12=[c1(:);c2(:)];
for cc=c12'
    units_c=find(~cellfun(@isempty,{in(cc).unit.average_spike_density}));
    units=intersect(units,units_c);
end

if ~isempty(units)
    C1=[];
    C2=[];
    for cc=c1(:)'
        C1=[C1;vertcat(in(cc).unit(units).average_spike_density)];
    end
    for cc=c2(:)'
        C2=[C2;vertcat(in(cc).unit(units).average_spike_density)];
    end
    h=ttest2(C1,C2,0.05,'both','equal',1); h(h==0)=NaN;
else
    h=NaN; %% size = number of bins ?
end
end
