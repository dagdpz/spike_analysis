function ph_population_analysis_contra_ipsi(population,modified_keys)
%load('W:\Projects\Pulv_microstim_behavior\ephys\ephys_analysis_v5_July2016_coordinates\monkeys_Combined_mem.mat')
global MA_STATES
warning('off','MATLAB:catenate:DimensionMismatch');
%run 'epoch_definitions';

%% Keys
% 
% keys.anova_table_file                   = 'W:\Projects\Pulv_eye_hand\ephys\ephys_analysis_14_July2016_coordinates_minus_baseline\monkeys_Combined_del.mat';
% keys.effectors_on_same_figure           = 1;
% keys.PSTH_binwidth                      = 0.01;
% keys.population_normalization           = 'by_effector';
% keys.tt.selection               ={'target','dPulv_l'};
% keys.tt.unselect                ={};
% keys.population_group_parameter         ='in_Cue_spaceLR_Ddre_han';
% keys.combine_tuning_properties          ={'place_name_here'};
% keys.population_group_excluded          ={'','-'};
% keys.epoch_RF     = 'Cue';
% keys.epoch_for_normalization            = 'Cue';
% keys.conditions_to_plot                 = {'Ddre'};
% keys.monkey                             = {'Combined'};
% 
% keys.drive                  ='W:';
% keys.basepath_to_save       ='Projects\Default\ephys';
% keys.pdf_folder             ='ephys_analysis_v3_June2016';
% keys.create_pdfs            =1;
% keys.append_pdfs            =0;
% keys.plot_RF                =1;
% keys.population_ylim        =[];

keys.saccade_states={'PreS','PeriS','PostS'};
keys.reach_states={'PreR','PeriR','PostR'};

for fn=fieldnames(modified_keys)'
    keys.(fn{:})=modified_keys.(fn{:});
end
keys.path_to_save=[keys.drive filesep keys.basepath_to_save filesep keys.pdf_folder filesep 'population_analysis' filesep];
if ~exist(keys.path_to_save,'dir')
    mkdir([keys.drive filesep keys.basepath_to_save filesep keys.pdf_folder ], 'population_analysis');
end


if strcmp(keys.contra_color,'pink')
    keys.hnd_choice_colors_C    =[255 0 178; 127 0 89; 171 0 252; 85 0 126; 145 143 56; 77 76 28]/255; %NH IN, NH CH, CH IN, CH CH, IH IN, IH CH;
    keys.hnd_choice_colors_I    =[255 153 20; 127 77 10; 125 130 255; 62 65 127; 222 220 0; 110 110 0]/255; %NH IN, NH CH, CH IN, CH CH, IH IN, IH CH;
    keys.hnd_choice_colors      =[1 0 0; 0.7 0 0; 0 0 1; 0 0 0.4; 0 1 0; 0 0.5 0]; %NH IN, NH CH, CH IN, CH CH, IH IN, IH CH;
    keys.space_colors            =[255 0 150;255 150 0]/255; %CS, IS, CH, IH;
else
    keys.hnd_choice_colors_I    =[255 0 178; 127 0 89; 145 143 56; 77 76 28; 171 0 252; 85 0 126]/255; %NH IN, NH CH, CH IN, CH CH, IH IN, IH CH;
    keys.hnd_choice_colors_C    =[255 153 20; 127 77 10; 222 220 0; 110 110 0; 125 130 255; 62 65 127]/255; %NH IN, NH CH, CH IN, CH CH, IH IN, IH CH;
    keys.hnd_choice_colors      =[1 0 0; 0.7 0 0; 0 1 0; 0 0.5 0; 0 0 1; 0 0 0.4]; %NH IN, NH CH, CH IN, CH CH, IH IN, IH CH;
    keys.space_colors            =[255 150 0;255 0 150]/255; %CS, IS, CH, IH;
end


%% tuning table preparation and grouping
[tuning_per_unit_table]                 = ph_load_extended_tuning_table(keys);
[tuning_per_unit_table, Sel_for_title]  = ph_reduce_tuning_table(tuning_per_unit_table,keys);
if keys.FR_subtract_baseline
    Sel_for_title =[Sel_for_title,{'base';'=';keys.epoch_BL;', '}];
end
idx_group_parameter=find_column_index(tuning_per_unit_table,keys.population_group_parameter);
group_values=tuning_per_unit_table(:,idx_group_parameter);
group_values=cellfun(@num2str, group_values, 'UniformOutput', false);
%unique_group_values=unique(group_values(2:end));
cell_in_any_group=[false; ~ismember(group_values(2:end),keys.population_group_excluded)];
unique_group_values=unique(group_values(cell_in_any_group));


%% markers for different monkeys and colors for different grid holes
idx_grid_x=find_column_index(tuning_per_unit_table,'grid_x');
idx_grid_y=find_column_index(tuning_per_unit_table,'grid_y');
idx_grid_z=find_column_index(tuning_per_unit_table,'electrode_depth');

monkey_markers={'o','s'};
monkey_linecolor={'y','m'};
gridhole_colors={[70 0 0;255 0 0],[0 70 0;0 255 0],[0 70 70; 0 255 255]};

idx_unitID=find_column_index(tuning_per_unit_table,'unit_ID');
monkey_per_cell=cellfun(@(x) x(1:3),tuning_per_unit_table(:,idx_unitID),'uniformoutput',false);
[unique_monkeys]=unique(monkey_per_cell(2:end));
for m=1:numel(unique_monkeys)
    idx_m=cell_in_any_group&ismember(monkey_per_cell,unique_monkeys(m));
    [unique_grid2_values{m}]=unique(cell2mat([tuning_per_unit_table(idx_m,idx_grid_x), tuning_per_unit_table(idx_m,idx_grid_y)]),'rows');
    
    if ~isempty(unique_grid2_values{m})
        unique_gridx_values{m}=[unique_grid2_values{m}(:,1)];
        unique_gridy_values{m}=[unique_grid2_values{m}(:,2)];
        
        %electrode_depth_range{m}=[min([tuning_per_unit_table{idx_m,idx_grid_z}]) max([tuning_per_unit_table{idx_m,idx_grid_z}])];
        for gh=1:numel(unique_gridx_values{m})
            idx_g=[false; ismember([tuning_per_unit_table{2:end,idx_grid_x}],unique_gridx_values{m}(gh))' & ismember([tuning_per_unit_table{2:end,idx_grid_y}],unique_gridy_values{m}(gh) )'];
            electrode_depth_range{m}{gh}=[min([tuning_per_unit_table{idx_m & idx_g,idx_grid_z}]) max([tuning_per_unit_table{idx_m & idx_g,idx_grid_z}])];
        end
    end
end


if isempty(unique_group_values)
    disp('no relevant groups found');
    return;
end

hands=[0,1,2];
choices=[0,1];
choices_hands=combvec(choices,hands);
contra_labels={'NH CS IN' 'NH CS CH' 'CH CS IN' 'CH CS CH' 'IH CS IN' 'IH CS CH'};
ipsi_labels={'NH IS IN' 'NH IS CH' 'CH IS IN' 'CH IS CH' 'IH IS IN' 'IH IS CH'};
%lines=find(choices_hands(1,:)==0); %% choices are automatically excluded here... :(
lines=1:size(choices_hands,2);

%% gaussian fit settings
fitsettings.sd_max_x=12;
% fitsettings.sd_max_y=12;
%fitsettings.sd_x_min_ratio=0.25;
%fitsettings.sd_y_min_ratio=0.25;

fitsettings.sd_x_min_ratio=0.25;
fitsettings.sd_max_y=fitsettings.sd_max_x;
fitsettings.sd_xy_min_ratio=0.25;
fitsettings.sd_xy_max_ratio=1;
fitsettings.sd_y_min_ratio=fitsettings.sd_x_min_ratio;
fitsettings.xout=[-30:30]; % range for gaussian fit
fitsettings.yout=[-15:15];
fitsettings.range_factor=1;

FNpop={'left','right','position'};
FNCI={'contra','ipsi','position'}; % counter side, important for normalization
CIFN={'contra','ipsi'};   %for hand tuning
complete_unit_list={population.unit_ID}';

units_per_condition=arrayfun(@(x,y) repmat(y,1,numel(x.condition)),population,1:numel(population),'UniformOutput',false);
units_per_condition=[units_per_condition{:}];

all_conditions=[population.condition];
all_effectors=vertcat(all_conditions.effector);
all_types=vertcat(all_conditions.type);
[unique_conditions, con_idx, condition_indexes]=unique([all_effectors all_types],'rows');
%n_conditions=size(unique_conditions,1);
conditions=all_conditions(con_idx);
for c=1:numel(conditions)
    [~, type_effector_short]=get_type_effector_name(conditions(c).type,conditions(c).effector);
    if ~ismember(type_effector_short,keys.conditions_to_plot)
        continue;
    end
    conditions(c).case=rmfield(conditions(c).case,FNpop);
    eff=conditions(c).effector;
    typ=conditions(c).type;
    keys.EPOCHS=keys.EPOCHS_PER_TYPE{typ};
    normalization_epoch     =find(ismember(keys.EPOCHS(:,1),keys.epoch_for_normalization));    
    receptive_field_epoch   =find(ismember(keys.EPOCHS(:,1),keys.epoch_RF));
    baseline_epoch          =find(ismember(keys.EPOCHS(:,1),keys.epoch_BL));
    for a=1:numel(conditions(c).case)
        for g=1:numel(unique_group_values)
            g_idx= find(ismember(complete_unit_list,tuning_per_unit_table(ismember(group_values,unique_group_values(g)),idx_unitID)));
            u_idx= ismember(units_per_condition,g_idx)' & condition_indexes==c;
            if ~any(u_idx)
                continue;
            end
            c_units=units_per_condition(u_idx);
            
            %% position to subplot assignment struggle
            con_cas_pos     =vertcat(all_conditions(u_idx).case);
            con_cas_pos     ={con_cas_pos(:,a).position};
            subplotidx_tmp  =cellfun(@(x) 1:numel(x),con_cas_pos,'uniformoutput',false);
            con_cas_pos     =[con_cas_pos{:}];
            subplotidx_tmp  =[subplotidx_tmp{:}];
            valid_pos       =arrayfun(@(x) ~isempty(x.position),con_cas_pos);
            con_cas_pos     =con_cas_pos(valid_pos); %% empty position in cornz??
            subplotidx_tmp  =num2cell(subplotidx_tmp(valid_pos));
            [con_cas_pos.subplot_index] =subplotidx_tmp{:};
            [u_con_cas_pos,upos_idx]    =unique(round(vertcat(con_cas_pos.position)*100),'rows');  %% soritng positions !!!???
            [~,posresortindex]          =sort([subplotidx_tmp{upos_idx}]);
            u_con_cas_pos               =u_con_cas_pos(posresortindex,:)/100;
            
            %% normalization factor preparation
            for u=1:numel(c_units)
                pu = c_units(u);
                unit_effectors = [population(pu).condition.effector];
                unit_types = [population(pu).condition.type];
                uc = unit_effectors ==eff & unit_types == typ;
                %u_tmp= rmfield(population(pu),'condition');
                uL_per_state=[population(pu).condition(uc).case(a).left.per_state];
                uR_per_state=[population(pu).condition(uc).case(a).right.per_state];
                lineL_per_state=vertcat(uL_per_state.trial);
                lineR_per_state=vertcat(uR_per_state.trial);
                uL_trials=[uL_per_state(1).trial];
                uR_trials=[uR_per_state(1).trial];
                factor_tmp=[];
                for n=1:numel(lines)
                    ll=lines(n);
                    tr_idx_per_lineL{n}= [uL_trials.choice]==choices_hands(1,ll) & [uL_trials.hand]==choices_hands(2,ll);
                    tr_idx_per_lineR{n}= [uR_trials.choice]==choices_hands(1,ll) & [uR_trials.hand]==choices_hands(2,ll);
                    switch keys.population_normalization
                        case 'by_effector'
                            tmp_SDL=[vertcat(lineL_per_state(normalization_epoch,tr_idx_per_lineL{n}).spike_density)];
                            tmp_SDR=[vertcat(lineR_per_state(normalization_epoch,tr_idx_per_lineR{n}).spike_density)];
                            factor_tmp=[factor_tmp mean(tmp_SDL,1)  mean(tmp_SDR,1)];
                        case 'none'
                            factor_tmp= 1;
                    end
                    
                    %                             switch keys.population_normalization
                    %                                 case 'by_condition'
                    %                                     factor= nanmean(nanmean([line_per_state(normalization_epoch,tr_idx).spike_density]));
                    %                                 case 'by_effector'
                    %                                     factor= max(nanmean([line_per_state(normalization_epoch,:).spike_density linec_per_state(normalization_epoch,:).spike_density]));
                    %                                     %factor= norm_factor(u);
                    %                                 case 'by_type'
                    %                                 case 'by_all_trials'
                    %                             end
                end
                norm_factor(u)=max([factor_tmp,NaN]);
            end
            
            %% main part
            for fn=1:numel(FNCI)
                clear fn_tmp u_tmp
                for l=lines
                    if strcmp(FNCI{fn},'position')
                        for p=1:size(u_con_cas_pos,1)
                            pos_units=con_cas_pos([con_cas_pos.subplot_index]==p);
                            fn_tmp(p).position=u_con_cas_pos(p,:);
                            for u=1:numel(pos_units)
                                u_per_state=[pos_units(u).per_state];
                                u_trials=[u_per_state(1).trial];
                                tr_idx= [u_trials.choice]==choices_hands(1,l) & [u_trials.hand]==choices_hands(2,l);
                                line_per_state=vertcat(u_per_state.trial);
                                u_tmp.position=pos_units(u).position;
                                u_tmp.per_state=rmfield(pos_units(u).per_state,'trial');
                                fn_tmp(p).line(l).unit(u)=u_tmp;
                                if keys.FR_subtract_baseline
                                    baseline=mean(vertcat(line_per_state(baseline_epoch,tr_idx).spike_density),2); %% put here the epoch subtracted instead of 1 ... !!!
                                else
                                    baseline=zeros(sum(tr_idx),1);
                                end
                                for s= 1:numel(fn_tmp(p).line(l).unit(u).per_state)
                                    fn_tmp(p).line(l).unit(u).per_state(s).trial=line_per_state(s,tr_idx);
                                    Spike_density=vertcat(fn_tmp(p).line(l).unit(u).per_state(s).trial.spike_density);
                                    correction=baseline(arrayfun(@(x) ~isempty(x.spike_density), line_per_state(s,tr_idx)),ones(1,size(Spike_density,2)));
                                    if isempty(correction); correction=0; end;
                                    fn_tmp(p).line(l).unit(u).per_state(s).average_spike_density=...
                                        mean((Spike_density-correction)/(norm_factor(u)),1);
                                    [~, fn_tmp(p).line(l).unit(u).per_state(s).peak_bin]=max(fn_tmp(p).line(l).unit(u).per_state(s).average_spike_density);
                                end
                                
                                %                                  fn_tmp(p).line(l).unit(u).per_state(s).average_spike_density=...
                                %                                         mean(vertcat(fn_tmp(p).line(l).unit(u).per_state(s).trial.spike_density)/norm_factor(u),1);
                                %                                     [~, fn_tmp(p).line(l).unit(u).per_state(s).peak_bin]=max(fn_tmp(p).line(l).unit(u).per_state(s).average_spike_density);
                                
                            end
                        end
                    else
                        for u=1:numel(c_units)
                            pu = c_units(u);
                            unit_effectors = [population(pu).condition.effector];
                            unit_types = [population(pu).condition.type];
                            uc = unit_effectors ==eff & unit_types == typ;
                            %u_tmp= rmfield(population(pu),'condition');
                            
                            hand_CI=0;
                            
                            Left_hemisphere_targets={'dPulv_l','pdSTS_L','FST_L','pTPO_L','MIP','MIP_L','unknown'};
                            Right_hemisphere_targets={'dPulv_r'};
                            if ismember(population(pu).target,Left_hemisphere_targets)
                                if strcmp(FNCI{fn},'contra')
                                    FNLR= 'right';
                                elseif strcmp(FNCI{fn},'ipsi')
                                    FNLR= 'left';
                                end
                                if choices_hands(2,l)==1
                                    hand_CI=2;
                                elseif choices_hands(2,l)==2
                                    hand_CI=1;
                                end
                            elseif ismember(population(pu).target,Right_hemisphere_targets)
                                if strcmp(FNCI{fn},'contra')
                                    FNLR= 'left';
                                elseif strcmp(FNCI{fn},'ipsi')
                                    FNLR= 'right';
                                end
                                if choices_hands(2,l)==1
                                    hand_CI=1;
                                elseif choices_hands(2,l)==2
                                    hand_CI=2;
                                end
                            end
                            
                            u_per_state=[population(pu).condition(uc).case(a).(FNLR).per_state];
                            u_trials=[u_per_state(1).trial];
                            tr_idx= [u_trials.choice]==choices_hands(1,l) & [u_trials.hand]==hand_CI;
                            line_per_state=vertcat(u_per_state.trial);
                            u_tmp=rmfield(population(pu),'condition');
                            u_tmp.per_state=rmfield(population(pu).condition(uc).case(a).(FNLR).per_state,'trial');
                            fn_tmp.line(l).unit(u)=u_tmp;
                            
                            if keys.FR_subtract_baseline
                                baseline=mean(vertcat(line_per_state(baseline_epoch,tr_idx).spike_density),2); %% put here the epoch subtracted instead of 1 ... !!!
                            else
                                baseline=zeros(sum(tr_idx),1);
                            end
                            for s= 1:numel(fn_tmp.line(l).unit(u).per_state)
                                fn_tmp.line(l).unit(u).per_state(s).trial=line_per_state(s,tr_idx);
                                Spike_density=vertcat(fn_tmp.line(l).unit(u).per_state(s).trial.spike_density);
                                correction=baseline(arrayfun(@(x) ~isempty(x.spike_density), line_per_state(s,tr_idx)),ones(1,size(Spike_density,2)));
                                if isempty(correction); correction=0; end;
                                fn_tmp.line(l).unit(u).per_state(s).average_spike_density=...
                                    mean((Spike_density-correction)/(norm_factor(u)),1);
                                [~, fn_tmp.line(l).unit(u).per_state(s).peak_bin]=max(fn_tmp.line(l).unit(u).per_state(s).average_spike_density);
                            end
                        end
                    end
                end
                conditions(c).case(a).group(g).(FNCI{fn})= fn_tmp;
            end
            %% Gaussian receptive field fit
            for l=lines %% ??? not even caring about hand... either select this properly or leave it !!!!!!!
                if ~keys.plot_RF
                    continue;
                end
                for u=1:numel(c_units)
                    pu = c_units(u);
                    unit_effectors = [population(pu).condition.effector];
                    unit_types = [population(pu).condition.type];
                    uc = unit_effectors ==eff & unit_types == typ;
                    u_per_state=vertcat(population(pu).condition(uc).case(a).('position').per_state);
                    un_con_cas_pos=vertcat(population(pu).condition(uc).case(a).('position').position);
                    for s=1:size(u_per_state,2)
                        xin=un_con_cas_pos(:,1);
                        yin=un_con_cas_pos(:,2);
                        pos_trials={u_per_state(:,s).trial};
                        
                        
                        %tr_idx= [pos_trials.choice]==choices_hands(1,l); %& [u_trials.hand]==hand_CI;    TODO
                        pos_trials=cellfun(@(x) x([x.choice]==choices_hands(1,l)),pos_trials,'UniformOutput',false); %%cornz strikes again
                        %pos_trials(cellfun(@isempty,pos_trials))=[]; %%cornz strikes again
                        zin=cellfun(@(x) nanmean(double([x.FR])),pos_trials); %%double because single in Linus/Curius Pulvinar gaze
                        FR_tmp=struct('FR',cellfun(@(x) [x.FR],pos_trials,'uniformoutput',false),'x',num2cell(xin'),'y',num2cell(yin'));
                        if ~any(isnan(zin)) && numel(zin)>=4 && numel(unique(zin))>1
                            [ftmp.Zout ftmp.sx ftmp.sy ftmp.phi ftmp.xmax ftmp.ymax ftmp.zmax]=fit_gaussian(xin,yin,zin,fitsettings);
                        else
                            [ftmp.Zout ftmp.sx ftmp.sy ftmp.phi ftmp.xmax ftmp.ymax ftmp.zmax]=deal(NaN);
                        end
                        conditions(c).case(a).group(g).gaussian_fit.line(l).per_state(s).unit(u).parameters=ftmp;
                        conditions(c).case(a).group(g).gaussian_fit.line(l).per_state(s).unit(u).positions=FR_tmp;
                        conditions(c).case(a).group(g).gaussian_fit.line(l).per_state(s).unit(u).ID=population(pu).unit_ID;
                        conditions(c).case(a).group(g).gaussian_fit.line(l).per_state(s).unit(u).monkey=population(pu).unit_ID(1:3);
                        conditions(c).case(a).group(g).gaussian_fit.line(l).per_state(s).unit(u).X=population(pu).grid_x;
                        conditions(c).case(a).group(g).gaussian_fit.line(l).per_state(s).unit(u).Y=population(pu).grid_y;
                        conditions(c).case(a).group(g).gaussian_fit.line(l).per_state(s).unit(u).Z=population(pu).electrode_depth;
                    end
                end
            end
            
            %% significance space tuning
            for l=lines
                Contra  =vertcat(conditions(c).case(a).group(g).contra.line(l).unit.per_state);
                Ipsi    =vertcat(conditions(c).case(a).group(g).ipsi.line(l).unit.per_state);
                
                for s= 1:size(Contra,2)
                    C_state=vertcat(Contra(:,s).average_spike_density);
                    I_state=vertcat(Ipsi(:,s).average_spike_density);
                    if ~isempty(C_state) && ~isempty(I_state)
                        %h=ttest2(C_state,I_state,0.05,'both','equal',1); h(h==0)=NaN; %unpaired!
                        h=ttest(C_state,I_state,0.05,'both',1); h(h==0)=NaN; %unpaired!
                    else
                        h=NaN(1, max([size(C_state,2) size(I_state,2)]));
                    end
                    conditions(c).case(a).group(g).space_tuning.line(l).per_state(s).h=h;
                end
            end
            
            %% significance hand tuning
            for ic=0:1
                CH_idx=find(choices_hands(2,:)==1 & choices_hands(1,:)==ic );
                IH_idx=find(choices_hands(2,:)==2 & choices_hands(1,:)==ic );
                if ~ismember(CH_idx,lines) || ~ismember(IH_idx,lines)
                    continue
                end
                for side=1:2
                    LH=vertcat(conditions(c).case(a).group(g).(CIFN{side}).line(CH_idx).unit.per_state);
                    RH=vertcat(conditions(c).case(a).group(g).(CIFN{side}).line(IH_idx).unit.per_state);
                    for s= 1:size(LH,2)
                        L_state=vertcat(LH(:,s).average_spike_density);
                        R_state=vertcat(RH(:,s).average_spike_density);
                        if ~isempty(L_state) && ~isempty(R_state)
                            h=ttest2(L_state,R_state,0.05,0,1,1); h(h==0)=NaN;
                        else
                            h=NaN(1, max([size(L_state,2) size(R_state,2)]));
                        end
                        conditions(c).case(a).group(g).hand_tuning(ic+1).(CIFN{side}).per_state(s).h=h;
                    end
                end
            end
            
        end
    end
end

%% plot
for c=1:numel(conditions)
    [type_effector type_effector_short]=get_type_effector_name(conditions(c).type,conditions(c).effector);
    if ~ismember(type_effector_short,keys.conditions_to_plot)
        continue;
    end
    get_expected_MA_states(conditions(c).type,conditions(c).effector,keys.effectors_on_same_figure);
    %get_expected_MA_states(conditions(c).type,conditions(c).effector);
    keys.EPOCHS=keys.EPOCHS_PER_TYPE{conditions(c).type};
    %keys.EPOCHS= keys.EPOCHS(ismember(keys.EPOCHS(:,1),keys.epochs_to_plot_PSTH{conditions(c).type}),:);
    epoch_to_plot_PSTH=ismember(keys.EPOCHS(:,1),keys.epochs_to_plot_PSTH{conditions(c).type})';
    states_to_plot=[keys.EPOCHS{:,2}];
    [~,state_order_indexes]=ismember(states_to_plot,MA_STATES.all_states);
    [sta_index_in_order,all_sta]=sort(state_order_indexes);
    sta_in_order=all_sta(sta_index_in_order>0 & epoch_to_plot_PSTH(all_sta));
    get_expected_MA_states(conditions(c).type,conditions(c).effector,0);
    
    receptive_field_epoch=find(ismember(keys.EPOCHS(:,1),keys.epoch_RF));
    %baseline_epoch=find(ismember(keys.EPOCHS(:,1),keys.epoch_BL));
    for a=1:numel(conditions(c).case)
        fig_title=sprintf('Selection: %s %s %s %s grouped by %s normalized %s in %s',...
            [Sel_for_title{:}],keys.monkey,type_effector,keys.case_summaries{a},keys.population_group_parameter, keys.population_normalization,keys.epoch_for_normalization);
        filename=sprintf('%s %s %s %s c%d N_%s %s',[Sel_for_title{:}],keys.monkey,keys.population_group_parameter,type_effector_short,a,keys.population_normalization,keys.epoch_for_normalization);
        plot_1_title            = [fig_title  ' PSTHs'];
        PSTH_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_1_title);
        
        state_onsets_to_plot=cell(numel(conditions(c).case(a).group),max(all_sta));
        %% PSTH plot
        for g=1:numel(conditions(c).case(a).group)
            if isempty(conditions(c).case(a).group(g).contra)
                continue;
            end
            subplot(numel(conditions(c).case(a).group),1,g)
            hold on
            n_units=[];
            %handtuningplotted=0;
            for l=lines
                Contra=vertcat(conditions(c).case(a).group(g).contra.line(l).unit.per_state);
                Ipsi=vertcat(conditions(c).case(a).group(g).ipsi.line(l).unit.per_state);
                n_units(l)=sum(arrayfun(@(x,y) ~isempty(x.per_state(1).average_spike_density) | ~isempty(y.per_state(1).average_spike_density),...
                    conditions(c).case(a).group(g).contra.line(l).unit,conditions(c).case(a).group(g).ipsi.line(l).unit));
                state_shift=0;
                for s=sta_in_order
                    state_label=keys.EPOCHS{s,1};
                    t_before_state=keys.EPOCHS{s,3};
                    t_after_state=keys.EPOCHS{s,4};
                    bins=t_before_state:keys.PSTH_binwidth:t_after_state;
                    bins=bins+state_shift-t_before_state;
                    if ~isempty(mean(vertcat(Contra(:,s).average_spike_density))) && ismember(keys.EPOCHS{s,2},MA_STATES.all_states)
                        props={'color',keys.hnd_choice_colors_C(l,:),'linewidth',1};
                        shadedErrorBar(bins,nanmean(vertcat(Contra(:,s).average_spike_density),1),...
                            sterr(vertcat(Contra(:,s).average_spike_density),1),props,1)
                    end
                    if ~isempty(mean(vertcat(Ipsi(:,s).average_spike_density)))  && ismember(keys.EPOCHS{s,2},MA_STATES.all_states)
                        props={'color',keys.hnd_choice_colors_I(l,:),'linewidth',1};
                        shadedErrorBar(bins,nanmean(vertcat(Ipsi(:,s).average_spike_density),1),...
                            sterr(vertcat(Ipsi(:,s).average_spike_density),1),props,1)
                    end
                    Tr=[Ipsi(:,s).trial  Contra(:,s).trial];
                    State_onsets{g}{l,s}=[Tr.state_onset];
                    state_shift=state_shift+t_after_state-t_before_state+0.1;
                end
                %handtuningplotted=1;
                
                %% mean state onsets (relative to current state)
                for s=all_sta
                    contra_trials=[Contra(:,s).trial]';
                    ipsi_trials=[Ipsi(:,s).trial]';
                    mean_state_onsets(l,s).contra=nanmean(vertcat(contra_trials.state_onsets),1);
                    mean_state_onsets(l,s).ipsi=nanmean(vertcat(ipsi_trials.state_onsets),1);
                end
            end
            title(sprintf('%s = %s, N (NH/CH/IH) =%s ',keys.population_group_parameter,unique_group_values{g},num2str(n_units(lines))),'interpreter','none');
            y_lim(g,:)=get(gca,'ylim');
            
            
        for s=1:size(state_onsets_to_plot,2)
            state_onsets_to_plot{g,s}=mean([vertcat(mean_state_onsets(:,s).contra); vertcat(mean_state_onsets(:,s).ipsi)],1);
        end
        end
        
        %% subplot appearance, and tuning lines
        ylimmax=max(max(y_lim));
        ylimmin=min(min(y_lim));
        y_lim=[ylimmin ylimmax];
        if ~isempty(keys.population_ylim)
            y_lim= keys.population_ylim;
        end
        for g=1:numel(conditions(c).case(a).group)
            if isempty(conditions(c).case(a).group(g).contra)
                continue;
            end
            subplot(numel(conditions(c).case(a).group),1,g);
            state_shift=0;
            for s=sta_in_order
                State_onset=mean([State_onsets{g}{:,s}]);
                t_before_state=keys.EPOCHS{s,3};
                t_after_state=keys.EPOCHS{s,4};
                bins=t_before_state:keys.PSTH_binwidth:t_after_state;
                bins=bins+state_shift-t_before_state;
                
                for sta_sta=1:numel(all_sta)
                    sta_t=all_sta(sta_sta);
                    if ~ismember(keys.EPOCHS{s,2},MA_STATES.all_states) || ~ismember(keys.EPOCHS{sta_t,2},MA_STATES.all_states)
                       continue; 
                    end
                    relative_state_onset=state_onsets_to_plot{g,s};
                    state_shift_relative=state_shift-t_before_state+relative_state_onset(sta_t);
                    state_label     =keys.EPOCHS{sta_t,1};
                    t_b             =keys.EPOCHS{sta_t,5};
                    t_a             =keys.EPOCHS{sta_t,6};
                    if state_shift_relative>=state_shift-0.0001 && state_shift_relative<= state_shift+t_after_state-t_before_state+0.0001
                        if ismember(state_label,keys.saccade_states)
                            li=line([state_shift_relative state_shift_relative],y_lim,'color',[1 0 0]);
                            ptch=rectangle('Position',[t_b+state_shift_relative y_lim(1) t_a-t_b diff(y_lim)*0.05],'EdgeColor','none','FaceColor',[1 0 0]); %frame for indicating FR window
                        elseif ismember(state_label,keys.reach_states)
                            li=line([state_shift_relative state_shift_relative],y_lim,'color',[0 1 0]);
                            ptch=rectangle('Position',[t_b+state_shift_relative y_lim(1) t_a-t_b diff(y_lim)*0.05],'EdgeColor','none','FaceColor',[0 1 0]); %frame for indicating FR window
                        else
                            li=line([state_shift_relative state_shift_relative],y_lim,'color',[0.8 0.8 0.8]);
                            ptch=rectangle('Position',[t_b+state_shift_relative y_lim(1) t_a-t_b diff(y_lim)*0.05],'EdgeColor','none','FaceColor',[0.8 0.8 0.8]); %frame for indicating FR window
                        end
                        uistack(ptch, 'bottom');
                        text(double(state_shift_relative+t_b),y_lim(2)+diff(y_lim)*-0.08,state_label,'fontsize',12);  %label for state
                    end
                end
                rectangle('Position',[state_shift y_lim(1) t_after_state-t_before_state diff(y_lim)],'linewidth',1) %frame for indicating state separation && normalization window
                
                Prev_state_onset=State_onset;
                
                %                 %% peak times
                %                                 Contra_handles=[];
                %                                 Contra_indexes=[];
                %                                 Ipsi_handles=[];
                %                                 Ipsi_indexes=[];
                %                                 for l=lines
                %                                     Contra=vertcat(conditions(c).case(a).group(g).contra.line(l).unit.per_state);
                %                                     Ipsi=vertcat(conditions(c).case(a).group(g).ipsi.line(l).unit.per_state);
                %                                     mean_c=(nanmean(vertcat(Contra(:,s).peak_bin),1)-1)*keys.PSTH_binwidth;
                %                                     mean_i=(nanmean(vertcat(Ipsi(:,s).peak_bin),1)-1)*keys.PSTH_binwidth;
                %                                     sterr_c=(sterr(vertcat(Contra(:,s).peak_bin),1))*keys.PSTH_binwidth;
                %                                     sterr_i=(sterr(vertcat(Ipsi(:,s).peak_bin),1))*keys.PSTH_binwidth;
                %                                     if ~isempty(mean_c)
                %                                     CL=line([mean_c mean_c]+state_shift,y_lim,'color',keys.hnd_choice_colors_C(l,:));
                %                                     line([mean_c+sterr_c mean_c+sterr_c]+state_shift,y_lim,'color',keys.hnd_choice_colors_C(l,:),'linestyle',':');
                %                                     line([mean_c-sterr_c mean_c-sterr_c]+state_shift,y_lim,'color',keys.hnd_choice_colors_C(l,:),'linestyle',':');
                %                                     Contra_handles=[Contra_handles CL];
                %                                     Contra_indexes=[Contra_indexes l];
                %                                     end
                %                                     if ~isempty(mean_i)
                %                                     IL=line([mean_i mean_i]+state_shift,y_lim,'color',keys.hnd_choice_colors_I(l,:));
                %                                     line([mean_i+sterr_i mean_i+sterr_i]+state_shift,y_lim,'color',keys.hnd_choice_colors_I(l,:),'linestyle',':');
                %                                     line([mean_i-sterr_i mean_i-sterr_i]+state_shift,y_lim,'color',keys.hnd_choice_colors_I(l,:),'linestyle',':');
                %                                     Ipsi_handles=[Ipsi_handles IL];
                %                                     Ipsi_indexes=[Ipsi_indexes l];
                %                                     end
                %                                 end
                
                %% space tuning
                for l=lines
                    l_height=y_lim(2)-diff(y_lim)*(l*0.01+0.1);
                    if ~isempty(conditions(c).case(a).group(g).space_tuning.line(l).per_state(s).h)
                        plot(bins,conditions(c).case(a).group(g).space_tuning.line(l).per_state(s).h*l_height,...
                            'linewidth',2,'color',keys.hnd_choice_colors(l,:))
                    end
                end
                
                %% hand tuning line
                for ic=0:1
                    CH_idx=find(choices_hands(2,:)==1 & choices_hands(1,:)==ic );
                    IH_idx=find(choices_hands(2,:)==2 & choices_hands(1,:)==ic );
                    if ~ismember(CH_idx,lines) || ~ismember(IH_idx,lines)
                        continue
                    end
                    for side=1:2
                        l_height=y_lim(2)-diff(y_lim)*((2*side+ic+5)*0.01+0.1);
                        if ~isempty(conditions(c).case(a).group(g).hand_tuning(ic+1).(CIFN{side}).per_state(s).h)
                            plot(bins,conditions(c).case(a).group(g).hand_tuning(ic+1).(CIFN{side}).per_state(s).h*l_height,...
                                'linewidth',2,'color',keys.space_colors(side,:))
                        end
                    end
                end
                state_shift=state_shift+t_after_state-t_before_state+0.1;
            end
            set(gca,'xlim',[0 state_shift-0.1],'ylim',y_lim);
            %legend([Contra_handles Ipsi_handles],[contra_labels(Contra_indexes),ipsi_labels(Ipsi_indexes)])
            remove_axis('x');
        end
        title_and_save(PSTH_summary_handle,  [filename ' PSTHs'],plot_1_title,keys)
        if true
            %% FR distribution plot (sanity check -.-)
            if strcmp(keys.population_normalization,'none')
                plot_AFR_title            = [fig_title  ' FR distribution'];
                Average_FR_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_AFR_title);
                for g=1:numel(conditions(c).case(a).group)
                    bins=0:3:100;
                    if isempty(conditions(c).case(a).group(g).contra)
                        continue;
                    end
                    subplot(numel(conditions(c).case(a).group),1,g)
                    hold on
                    for l=lines
                        Contra=vertcat(conditions(c).case(a).group(g).contra.line(l).unit.per_state);
                        Ipsi=vertcat(conditions(c).case(a).group(g).ipsi.line(l).unit.per_state);
                        FR_average=[];
                        uu=0;
                        for u=1:size(Contra,1)
                            if ~isempty([Contra(u,:).average_spike_density Ipsi(u,:).average_spike_density])
                                uu=uu+1;
                                FR_average(uu)=mean([Contra(u,:).average_spike_density Ipsi(u,:).average_spike_density]);
                            end
                        end
                        y_lim(g,:)=get(gca,'ylim');
                        text(min(bins),y_lim(g,1)+diff(y_lim(g,:))*l/10,sprintf('mean=%.2f,median=%.2f,std=%.2f',nanmean(FR_average),nanmedian(FR_average),nanstd(FR_average)),'interpreter','none');
                        plot(bins,hist(FR_average,bins),'color',keys.hnd_choice_colors_C(l,:),'linewidth',1)
                    end
                    %title(sprintf('%s = %s, N (NH/CH/IH) =%s ',keys.population_group_parameter,unique_group_values{g},num2str(n_units(lines))),'interpreter','none');
                    title(sprintf('%s = %s',keys.population_group_parameter,unique_group_values{g}),'interpreter','none');
                    
                    %
                end
                title_and_save(Average_FR_handle,  [filename ' FR distribution'],plot_AFR_title,keys)
            end
        end
        
        %% PSTH per position plot
        if ~keys.plot_RF
            continue;
        end
        
        for g=1:numel(conditions(c).case(a).group)
            plot_g_title            = [fig_title  ' PSTHs' sprintf('%s = %s', keys.population_group_parameter,unique_group_values{g})];
            PSTH_per_position_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_g_title);%,'interpreter','none')
            if isempty(conditions(c).case(a).group(g).position)
                continue;
            end
            y_lim_p=cell(0);
            n_units=[];
            for p=1:numel(conditions(c).case(a).group(g).position)
                subplot(conditions(c).case(a).rows,conditions(c).case(a).columns,p)
                hold on
                for l=lines
                    PSTH_p=vertcat(conditions(c).case(a).group(g).position(p).line(l).unit.per_state);
                    n_units(l)=sum(arrayfun(@(x) ~isempty(x.per_state(1).average_spike_density),...
                        conditions(c).case(a).group(g).position(p).line(l).unit));
                    state_shift=0;
                    for s=sta_in_order
                        t_before_state=keys.EPOCHS{s,3};
                        t_after_state=keys.EPOCHS{s,4};
                        bins=t_before_state:keys.PSTH_binwidth:t_after_state;
                        bins=bins+state_shift-t_before_state;
                        if ~isempty(mean(vertcat(PSTH_p(:,s).average_spike_density)))
                            props={'color',keys.hnd_choice_colors(l,:),'linewidth',1};
                            shadedErrorBar(bins,nanmean(vertcat(PSTH_p(:,s).average_spike_density),1),...
                                sterr(vertcat(PSTH_p(:,s).average_spike_density),1),props,1)
                        end
                        Tr=[PSTH_p(:,s).trial];
                        State_onsets{g}{l,s}=[Tr.state_onset];
                        state_shift=state_shift+t_after_state-t_before_state+0.1;
                    end
                end
                title(sprintf('Position  %.2f / %.2f N (NH/CH/IH) =%s ',conditions(c).case(a).group(g).position(p).position,num2str(n_units(lines))),'interpreter','none');
                y_lim_p{p}(g,:)=get(gca,'ylim');
            end
            
            
            %% subplot appearance, and tuning lines
            ylimmax=max(max([y_lim_p{:}]));
            ylimmin=min(min([y_lim_p{:}]));
            y_lim=[ylimmin ylimmax];
            if ~isempty(keys.population_ylim)
                y_lim= keys.population_ylim;
            end
            for  p=1:numel(conditions(c).case(a).group(g).position)
                %             if isempty(conditions(c).case(a).group(g).contra)
                %                 continue;
                %             end
                subplot(conditions(c).case(a).rows,conditions(c).case(a).columns,p)
                state_shift=0;
                Prev_state_onset=NaN;
                for s=sta_in_order
                    
                    State_onset=mean([State_onsets{g}{:,s}]);
                    state_label=keys.EPOCHS{s,1};
                    if strcmp(keys.epoch_for_normalization,state_label) && ~strcmp(keys.population_normalization,'none')
                        rect_linewidth=4;
                    else
                        rect_linewidth=1;
                    end
                    
                    t_before_state=keys.EPOCHS{s,3};
                    t_after_state=keys.EPOCHS{s,4};
                    bins=t_before_state:keys.PSTH_binwidth:t_after_state;
                    bins=bins+state_shift-t_before_state;
                    patch_x=[keys.EPOCHS{s,5}-keys.EPOCHS{s,3} keys.EPOCHS{s,6}-keys.EPOCHS{s,5}];
                    rectangle('Position',[patch_x(1)+state_shift y_lim(1) patch_x(2) diff(y_lim)/25],'EdgeColor','none','FaceColor',[0.8 0.8 0.8]) %frame for indicating state separation
                    rectangle('Position',[state_shift y_lim(1) t_after_state-t_before_state diff(y_lim)],'linewidth',rect_linewidth) %frame for indicating state separation && normalization window
                    
                    if t_before_state<=0 && t_after_state>=0
                        line([state_shift-t_before_state state_shift-t_before_state],y_lim,'color',[0.8 0.8 0.8]);
                        if ~isnan(State_onset) && ~(State_onset==Prev_state_onset)
                            text(state_shift-t_before_state,y_lim(1)+diff(y_lim)*-0.05,num2str(round(State_onset*100)/100),'fontsize',7,'color',[0.8 0.8 0.8]);
                        end
                    end
                    text(state_shift+0.02,y_lim(2)+diff(y_lim)*-0.08,state_label,'fontsize',10);
                    Prev_state_onset=State_onset;
                    
                    %    os=state_shift-t_before_state;
                    %% peak times
                    Contra_handles=[];
                    Contra_indexes=[];
                    Ipsi_handles=[];
                    Ipsi_indexes=[];
                    for l=lines
                        Contra=vertcat(conditions(c).case(a).group(g).contra.line(l).unit.per_state);
                        Ipsi=vertcat(conditions(c).case(a).group(g).ipsi.line(l).unit.per_state);
                        mean_c=(nanmean(vertcat(Contra(:,s).peak_bin),1)-1)*keys.PSTH_binwidth;
                        mean_i=(nanmean(vertcat(Ipsi(:,s).peak_bin),1)-1)*keys.PSTH_binwidth;
                        sterr_c=(sterr(vertcat(Contra(:,s).peak_bin),1))*keys.PSTH_binwidth;
                        sterr_i=(sterr(vertcat(Ipsi(:,s).peak_bin),1))*keys.PSTH_binwidth;
                        if ~isempty(mean_c)
                            CL=line([mean_c mean_c]+state_shift,y_lim,'color',keys.hnd_choice_colors_C(l,:));
                            line([mean_c+sterr_c mean_c+sterr_c]+state_shift,y_lim,'color',keys.hnd_choice_colors_C(l,:),'linestyle',':');
                            line([mean_c-sterr_c mean_c-sterr_c]+state_shift,y_lim,'color',keys.hnd_choice_colors_C(l,:),'linestyle',':');
                            Contra_handles=[Contra_handles CL];
                            Contra_indexes=[Contra_indexes l];
                        end
                        if ~isempty(mean_i)
                            IL=line([mean_i mean_i]+state_shift,y_lim,'color',keys.hnd_choice_colors_I(l,:));
                            line([mean_i+sterr_i mean_i+sterr_i]+state_shift,y_lim,'color',keys.hnd_choice_colors_I(l,:),'linestyle',':');
                            line([mean_i-sterr_i mean_i-sterr_i]+state_shift,y_lim,'color',keys.hnd_choice_colors_I(l,:),'linestyle',':');
                            Ipsi_handles=[Ipsi_handles IL];
                            Ipsi_indexes=[Ipsi_indexes l];
                        end
                    end
                    
                    %% space tuning
                    for l=lines
                        l_height=y_lim(2)-diff(y_lim)*(l*0.01+0.1);
                        if ~isempty(conditions(c).case(a).group(g).space_tuning.line(l).per_state(s).h)
                            plot(bins,conditions(c).case(a).group(g).space_tuning.line(l).per_state(s).h*l_height,...
                                'linewidth',2,'color',keys.hnd_choice_colors(l,:))
                        end
                    end
                    
                    %% hand tuning line
                    for ic=0:1
                        CH_idx=find(choices_hands(2,:)==1 & choices_hands(1,:)==ic );
                        IH_idx=find(choices_hands(2,:)==2 & choices_hands(1,:)==ic );
                        if ~ismember(CH_idx,lines) || ~ismember(IH_idx,lines)
                            continue
                        end
                        for side=1:2
                            l_height=y_lim(2)-diff(y_lim)*((2*side+ic+5)*0.01+0.1);
                            if ~isempty(conditions(c).case(a).group(g).hand_tuning(ic+1).(CIFN{side}).per_state(s).h)
                                plot(bins,conditions(c).case(a).group(g).hand_tuning(ic+1).(CIFN{side}).per_state(s).h*l_height,...
                                    'linewidth',2,'color',keys.space_colors(side,:))
                            end
                        end
                    end
                    state_shift=state_shift+t_after_state-t_before_state+0.1;
                end
                set(gca,'xlim',[0 state_shift-0.1],'ylim',y_lim);
                legend([Contra_handles Ipsi_handles],[contra_labels(Contra_indexes),ipsi_labels(Ipsi_indexes)])
                remove_axis('x');
            end
            title_and_save(PSTH_per_position_handle,  [filename ' PSTHs pos ' unique_group_values{g}],plot_g_title,keys)
            
        end
        
        %% receptive field plots
        
        angles=[0:pi/100:(2+1/100)*pi];
        circle_x=cos(angles);
        circle_y=sin(angles);
        for g=1:numel(conditions(c).case(a).group)
            if isempty(conditions(c).case(a).group(g).contra) % why only contra...?
                continue;
            end
            %fig_title=sprintf('Selection: %s Type: %d Effector: %d case %d',[Sel_for_title{:}],conditions(c).type,conditions(c).effector,a);
            plot_2_title            = [fig_title  ' RF ' unique_group_values{g}];
            Receptive_fields_handle = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_2_title);
            for l=lines
                for s=receptive_field_epoch
                    receptive_field_parameters=[conditions(c).case(a).group(g).gaussian_fit.line(l).per_state(s).unit.parameters];
                    Gaussian_fit={receptive_field_parameters.Zout};
                    Gaussian_fit_sx={receptive_field_parameters.sx};
                    Gaussian_fit_sy={receptive_field_parameters.sy};
                    Gaussian_fit_phi={receptive_field_parameters.phi};
                    Gaussian_fit_x={receptive_field_parameters.xmax};
                    Gaussian_fit_y={receptive_field_parameters.ymax};
                    %                     RF_size_x=[receptive_field_parameters.sx];
                    %                     RF_size_y=[receptive_field_parameters.sy];
                    for u=1:numel(Gaussian_fit)
                        subplot(ceil(sqrt(numel(Gaussian_fit))),ceil(sqrt(numel(Gaussian_fit))),u);
                        
                        center=[Gaussian_fit_x{u} Gaussian_fit_y{u}];
                        ellipse_r=2*Gaussian_fit_sx{u}*Gaussian_fit_sy{u}./sqrt(Gaussian_fit_sy{u}.^2*cos(angles + Gaussian_fit_phi{u}).^2+Gaussian_fit_sx{u}.^2*sin(angles  + Gaussian_fit_phi{u}).^2);
                        ellipse_x = circle_x.*ellipse_r;
                        ellipse_y = circle_y.*ellipse_r;
                        kk=cat(3,Gaussian_fit{u});
                        imagesc(fitsettings.xout,-fitsettings.yout,rot90(mean(kk,3)))
                        line(ellipse_x+center(1),ellipse_y+center(2),'linewidth',4,'color','k');
                        axis equal
                        set(gca,'Ydir','normal','Xtick',[],'Ytick',[],'xlim',[min(fitsettings.xout) max(fitsettings.xout)],'ylim',[min(fitsettings.yout) max(fitsettings.yout)]);
                        title(conditions(c).case(a).group(g).gaussian_fit.line(l).per_state(s).unit(u).ID,'interpreter','none')
                        text(0,0,num2str(4*sqrt(Gaussian_fit_sx{u}*Gaussian_fit_sy{u})));
                    end
                end
            end
            title_and_save(Receptive_fields_handle,  [filename ' RF ' unique_group_values{g}],plot_2_title,keys);
            
            %% Firing rate plot
            
            %fig_title=sprintf('Selection: %s Type: %d Effector: %d case %d',[Sel_for_title{:}],conditions(c).type,conditions(c).effector,a);
            plot_3_title            = [fig_title  ' FR' unique_group_values{g}];
            FR_summary_handle     = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_2_title);
            for l=lines
                for s=receptive_field_epoch
                    %subplot(numel(conditions(c).case(a).group),numel(sta_in_order),(g-1)*numel(sta_in_order) + s)
                    receptive_field_parameters=[conditions(c).case(a).group(g).gaussian_fit.line(l).per_state(s).unit.parameters];
                    %Gaussian_fit={receptive_field_parameters.Zout};
                    for u=1:numel(receptive_field_parameters)
                        subplot(ceil(sqrt(numel(receptive_field_parameters))),ceil(sqrt(numel(receptive_field_parameters))),u);
                        title(conditions(c).case(a).group(g).gaussian_fit.line(l).per_state(s).unit(u).ID,'interpreter','none')
                        hold on
                        positions=[conditions(c).case(a).group(g).gaussian_fit.line(l).per_state(s).unit(u).positions];
                        FR=arrayfun(@(x) nanmean(double([x.FR])),positions);
                        FR=FR-min(FR)+1;
                        FRmax=max(FR);
                        FR255=round(FR/FRmax*255);
                        cols=jet(255);
                        for p=1:numel(positions)
                            if ~isnan(FR255(p))
                                plot(positions(p).x,positions(p).y,'o','markerfacecolor',cols(FR255(p),:),'markersize',20);
                            end
                        end
                        axis equal
                        set(gca,'Xtick',[],'Ytick',[],'xlim',[min(fitsettings.xout) max(fitsettings.xout)],'ylim',[min(fitsettings.yout) max(fitsettings.yout)]);
                    end
                end
                
            end
            title_and_save(FR_summary_handle,  [filename ' FR ' unique_group_values{g}],plot_3_title,keys);
            
        end
        
        %% RF centers
        plot_2b_title            = [fig_title  ' RF centers'];
        RF_centers_handle = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_2b_title);
        
        maxRFsize=sqrt(fitsettings.sd_max_x^2 +  fitsettings.sd_max_y^2);
        minRFsize=sqrt((fitsettings.sd_max_x*fitsettings.sd_x_min_ratio)^2 +  (fitsettings.sd_max_y*fitsettings.sd_y_min_ratio)^2 );
        for g=1:numel(conditions(c).case(a).group)
            if isempty(conditions(c).case(a).group(g).gaussian_fit) % why only contra...?
                continue;
            end
            hold on
            for l=lines
                for s=receptive_field_epoch
                    RF_units=[conditions(c).case(a).group(g).gaussian_fit.line(l).per_state(s).unit];
                    RFparameters=[RF_units.parameters];
                    %RFsizes=sqrt([RFparameters.sx].^2+[RFparameters.sy].^2);
                    RFsizes=2*sqrt([RFparameters.sx].*[RFparameters.sy]);
                    [~, sort_by_size_index]=sort(RFsizes);
                    RF_units=RF_units(fliplr(sort_by_size_index));
                    RFsizes=RFsizes(fliplr(sort_by_size_index));
                    
                    RF_size_factor=0.1;
                    
                    for u=1:numel(RF_units)
                        M=ismember(unique_monkeys,RF_units(u).monkey);
                        colidx=ismember(unique_gridx_values{M},RF_units(u).X)&ismember(unique_gridy_values{M},RF_units(u).Y);
                        colfactor=(RF_units(u).Z-electrode_depth_range{M}{colidx}(1))/diff(electrode_depth_range{M}{colidx});
                        colfactor(isnan(colfactor))=1;
                        col=(gridhole_colors{colidx}(1,:)+ diff(gridhole_colors{colidx})*colfactor)/255;
                        
                        center=[RF_units(u).parameters.xmax RF_units(u).parameters.ymax];
                        Radius=RF_size_factor*RFsizes(u);
                        
                        if strcmp(monkey_markers{M},'o')
                            ellipse_x = circle_x.*Radius;
                            ellipse_y = circle_y.*Radius;
                            line(ellipse_x+center(1),ellipse_y+center(2),'color',col,'linewidth',4);
                            %line(ellipse_x+center(1),ellipse_y+center(2),'color',[1 1 1],'linewidth',1);
                        elseif strcmp(monkey_markers{M},'s')
                            square_x = [-1,-1,1,1,-1].*Radius*sqrt(pi)/2;
                            square_y = [-1,1,1,-1,-1].*Radius*sqrt(pi)/2;
                            line(square_x+center(1),square_y+center(2),'color',col,'linewidth',4);
                            %line(square_x+center(1),square_y+center(2),'color',[1 1 1],'linewidth',1);
                        end
                    end
                    axis equal
                    set(gca,'Ydir','normal','xlim',[min(fitsettings.xout) max(fitsettings.xout)],'ylim',[min(fitsettings.yout) max(fitsettings.yout)]);
                end
            end
            %title_and_save(RF_summary_handle,  [filename ' RF summary ' unique_group_values{g}],plot_2b_title,keys);
        end
        title_and_save(RF_centers_handle,  [filename ' RF centers '],plot_2b_title,keys);
        
        %% RF centers
        
        for M=1:numel(unique_monkeys)
            for  gd=1:numel(unique_gridx_values{M})
                title_part=sprintf('%s_%d_%d',unique_monkeys{M},unique_gridx_values{M}(gd),unique_gridy_values{M}(gd));
                plot_2b_title            = [fig_title  ' RF centers' title_part];
                RF_centers_handle = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_2b_title);
                
                
                maxRFsize=sqrt(fitsettings.sd_max_x^2 +  fitsettings.sd_max_y^2);
                minRFsize=sqrt((fitsettings.sd_max_x*fitsettings.sd_x_min_ratio)^2 +  (fitsettings.sd_max_y*fitsettings.sd_y_min_ratio)^2 );
                hold on
                for l=lines
                    for s=receptive_field_epoch
                        All_groups=[conditions(c).case(a).group];
                        All_fits=[All_groups.gaussian_fit];
                        if isempty(All_fits) % why only contra...?
                            continue;
                        end
                        Current_line=vertcat(All_fits.line);
                        Current_state=vertcat(Current_line(:,l).per_state);
                        
                        Current_units=Current_state(:,s);
                        RF_units=[Current_units.unit];
                        uidx_monkey=ismember(vertcat(RF_units.monkey),unique_monkeys);
                        uidx_X=ismember(vertcat(RF_units.X),unique_gridx_values{M}(gd));
                        uidx_Y=ismember(vertcat(RF_units.Y),unique_gridy_values{M}(gd));
                        RF_units=RF_units(uidx_monkey&uidx_X&uidx_Y);
                        if isempty(RF_units)
                            continue;
                        end
                        %% sort by Z
                        [~, sort_by_size_Z]=sort([RF_units.Z]);
                        RF_units=RF_units(sort_by_size_Z);
                        
                        for u=1:numel(RF_units)
                            %M=ismember(unique_monkeys,RF_units(u).monkey);
                            colfactor=(RF_units(u).Z-electrode_depth_range{M}{gd}(1))/diff(electrode_depth_range{M}{gd});
                            colfactor(isnan(colfactor))=1;
                            col=(gridhole_colors{gd}(1,:)+ diff(gridhole_colors{gd})*colfactor)/255;
                            
                            center=[RF_units(u).parameters.xmax RF_units(u).parameters.ymax];
                            Gaussian_fit_sx=[RF_units(u).parameters.sx];
                            Gaussian_fit_sy=[RF_units(u).parameters.sy];
                            Gaussian_fit_phi=[RF_units(u).parameters.phi];
                            
                            ellipse_r=Gaussian_fit_sx*Gaussian_fit_sy./sqrt(Gaussian_fit_sy.^2*cos(angles + Gaussian_fit_phi).^2+Gaussian_fit_sx.^2*sin(angles  + Gaussian_fit_phi).^2);
                            ellipse_x = circle_x.*ellipse_r;
                            ellipse_y = circle_y.*ellipse_r;
                            line(ellipse_x+center(1),ellipse_y+center(2),'linewidth',8,'color',col);
                        end
                        axis equal
                        set(gca,'Ydir','normal','xlim',[min(fitsettings.xout) max(fitsettings.xout)],'ylim',[min(fitsettings.yout) max(fitsettings.yout)]);
                    end
                end
                title_and_save(RF_centers_handle,  [filename ' RF centers ' title_part],plot_2b_title,keys);
            end
        end
        
        %% RF standard deviation
        %         plot_2b_title            = [fig_title  ' RF ellipses'];
        %         RF_ellipses_handle = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_2b_title);
        %         for g=1:numel(conditions(c).case(a).group)
        %             %         plot_2b_title            = [fig_title  ' RF summary' unique_group_values{g}];
        %             %         RF_summary_handle = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_2b_title);
        %             if isempty(conditions(c).case(a).group(g).gaussian_fit) % why only contra...?
        %                 continue;
        %             end
        %             hold on
        %             for l=lines
        %                 for s=receptive_field_epoch
        %                     receptive_field_parameters=[conditions(c).case(a).group(g).gaussian_fit.line(l).per_state(s).unit.parameters];
        %                     Gaussian_fit_sx={receptive_field_parameters.sx};
        %                     Gaussian_fit_sy={receptive_field_parameters.sy};
        %                     Gaussian_fit_phi={receptive_field_parameters.phi};
        %                     Gaussian_fit_x={receptive_field_parameters.xmax};
        %                     Gaussian_fit_y={receptive_field_parameters.ymax};
        %                     for u=1:numel(Gaussian_fit_sx)
        %                         center=[Gaussian_fit_x{u} Gaussian_fit_y{u}];
        %                         ellipse_r=Gaussian_fit_sx{u}*Gaussian_fit_sy{u}./sqrt(Gaussian_fit_sy{u}.^2*cos(angles + Gaussian_fit_phi{u}).^2+Gaussian_fit_sx{u}.^2*sin(angles  + Gaussian_fit_phi{u}).^2);
        %                         ellipse_x = circle_x.*ellipse_r;
        %                         ellipse_y = circle_y.*ellipse_r;
        %                         line(ellipse_x+center(1),ellipse_y+center(2),'linewidth',4,'color','k');
        %                         axis equal
        %                         set(gca,'Ydir','normal','xlim',[min(fitsettings.xout) max(fitsettings.xout)],'ylim',[min(fitsettings.yout) max(fitsettings.yout)]);
        %                     end
        %                 end
        %             end
        %             %title_and_save(RF_summary_handle,  [filename ' RF summary ' unique_group_values{g}],plot_2b_title,keys);
        %         end
        %         title_and_save(RF_ellipses_handle,  [filename ' RF ellipses '],plot_2b_title,keys);
        %
       
        %% RF size plots
        RF_size_stepsize=2;
        plot_4_title            = [fig_title  ' RF size ' unique_group_values{g}];
        Receptive_field_size_handle = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_4_title);
        for g=1:numel(conditions(c).case(a).group)
            if isempty(conditions(c).case(a).group(g).contra) % why only contra...?
                continue;
            end
            %fig_title=sprintf('Selection: %s Type: %d Effector: %d case %d',[Sel_for_title{:}],conditions(c).type,conditions(c).effector,a);
            for l=lines
                s=receptive_field_epoch;
                receptive_field_parameters=[conditions(c).case(a).group(g).gaussian_fit.line(l).per_state(s).unit.parameters];
                %Gaussian_fit={receptive_field_parameters.Zout};
                RF_size_x=4*[receptive_field_parameters.sx];
                RF_size_y=4*[receptive_field_parameters.sy];
                RF_size_a{g}=sqrt(RF_size_x.*RF_size_y);
                %RF_eccentricity=sqrt([receptive_field_parameters.xmax].^2 + [receptive_field_parameters.xmax].^2);
                RF_eccentricity{g}=[receptive_field_parameters.xmax];
                %end
            end
            RF_bins=0:RF_size_stepsize:(max([fitsettings.xout,fitsettings.yout])-min([fitsettings.xout,fitsettings.yout]));
            Ex_bins=min([fitsettings.xout]):RF_size_stepsize:max([fitsettings.xout]);
            
            subplot(numel(conditions(c).case(a).group),4,(g-1)*4+1)
            RF_x_hist=hist(RF_size_x,RF_bins);
            bar(RF_bins,RF_x_hist/sum(RF_x_hist));
            title(sprintf('%s = %s, RF size X, N=%d',keys.population_group_parameter,unique_group_values{g},sum(RF_x_hist)),'interpreter','none');
            set(gca,'xlim',([min(RF_bins)-RF_size_stepsize max(RF_bins)+RF_size_stepsize]),'ylim',[0 1])
            ylabel('N cells');
            xlabel('RF size X');
            text(min(RF_bins),0.8,sprintf('mean=%d.2f,median=%.2f,std=%.2f',nanmean(RF_size_x),nanmedian(RF_size_x),nanstd(RF_size_x)),'interpreter','none');
            
            subplot(numel(conditions(c).case(a).group),4,(g-1)*4+2)
            RF_y_hist=hist(RF_size_y,RF_bins);
            bar(RF_bins,RF_y_hist/sum(RF_y_hist));
            title(sprintf('%s = %s, RF size Y, N=%d',keys.population_group_parameter,unique_group_values{g},sum(RF_y_hist)),'interpreter','none');
            set(gca,'xlim',([min(RF_bins)-RF_size_stepsize max(RF_bins)+RF_size_stepsize]),'ylim',[0 1])
            ylabel('N cells');
            xlabel('RF size Y');
            text(min(RF_bins),0.8,sprintf('mean=%.2f,median=%.2f,std=%.2f',nanmean(RF_size_y),nanmedian(RF_size_y),nanstd(RF_size_y)),'interpreter','none');
            
            subplot(numel(conditions(c).case(a).group),4,(g-1)*4+3)
            RF_r_hist=hist(RF_size_a{g},RF_bins);
            bar(RF_bins,RF_r_hist/sum(RF_r_hist));
            title(sprintf('%s = %s, RF size euclidean, N=%d',keys.population_group_parameter,unique_group_values{g},sum(RF_r_hist)),'interpreter','none');
            set(gca,'xlim',([min(RF_bins)-RF_size_stepsize max(RF_bins)+RF_size_stepsize]),'ylim',[0 1])
            ylabel('N cells');
            xlabel('RF Radius');
            text(min(RF_bins),0.8,sprintf('mean=%.2f,median=%.2f,std=%.2f',nanmean(RF_size_a{g}),nanmedian(RF_size_a{g}),nanstd(RF_size_a{g})),'interpreter','none');
            text(min(RF_bins),0.6,sprintf('All groups: mean=%.2f,median=%.2f,std=%.2f',nanmean([RF_size_a{:}]),nanmedian([RF_size_a{:}]),nanstd([RF_size_a{:}])),'interpreter','none');
            
            subplot(numel(conditions(c).case(a).group),4,(g-1)*4+4)
            RF_ex_hist=hist(RF_eccentricity{g},Ex_bins);
            bar(Ex_bins,RF_ex_hist/sum(RF_ex_hist));
            title(sprintf('%s = %s, RF size excentricity, N=%d',keys.population_group_parameter,unique_group_values{g},sum(RF_ex_hist)),'interpreter','none');
            set(gca,'xlim',([min(Ex_bins)-RF_size_stepsize max(Ex_bins)+RF_size_stepsize]),'ylim',[0 1])
            ylabel('N cells');
            xlabel('excentricity');
            text(min(Ex_bins),0.8,sprintf('mean=%.2f,median=%.2f,std=%.2f',nanmean(RF_eccentricity{g}),nanmedian(RF_eccentricity{g}),nanstd(RF_eccentricity{g})),'interpreter','none');
            text(min(Ex_bins),0.6,sprintf('mean=%.2f,median=%.2f,std=%.2f',nanmean([RF_eccentricity{:}]),nanmedian([RF_eccentricity{:}]),nanstd([RF_eccentricity{:}])),'interpreter','none');
        end
        
        title_and_save(Receptive_field_size_handle,  [filename ' RF size' unique_group_values{g}],plot_4_title,keys);
        
    end
end
end

function [Zout sx sy phi xmax ymax zmax]=fit_gaussian(xin,yin,zin,fitsettings)
X=[xin(:),yin(:)];
xout=fitsettings.xout;
yout=fitsettings.yout;
Xout=combvec(xout,yout)';
Z=zin(:);
Z=Z-min(Z); %normalization

range_factor=fitsettings.range_factor;
x_range=(max(max(xin))-min(min(xin)))*range_factor;
y_range=(max(max(yin))-min(min(yin)))*range_factor;
sd_max_x=fitsettings.sd_max_x;
sd_max_y=fitsettings.sd_max_y;

% Fr weighted target position or Zmax as starting point
start_pos_xy=[mean(xin.*Z)/mean(Z) mean(yin.*Z)/mean(Z)];

%[~, maxidx]=max(Z);
%start_pos_xy=[xin(maxidx) yin(maxidx)];



zscaling=1;

% 
% F = @(a,x) a(4)*zscaling*max(Z)*...
%     exp(-((x(:,1)-a(5)*x_range)*cos(a(3))-(x(:,2)-a(6)*y_range)*sin(a(3))).^2/((2*a(1)*sd_max_x)^2)...
%         -((x(:,1)-a(5)*x_range)*sin(a(3))+(x(:,2)-a(6)*y_range)*cos(a(3))).^2/((2*a(2)*sd_max_y)^2));
% 
% LB=[fitsettings.sd_x_min_ratio fitsettings.sd_y_min_ratio -pi/4 0.5 min(min(xin))*range_factor/x_range min(min(yin))*range_factor/y_range];
% UB=[1 1 pi/4 1.5 max(max(xin))*range_factor/x_range max(max(yin))*range_factor/y_range]; %sigmax sigmay rotation peak x0 y0
% X0=[fitsettings.sd_x_min_ratio fitsettings.sd_y_min_ratio 0 1 start_pos_xy(1)*range_factor/x_range start_pos_xy(2)*range_factor/y_range];
% 
% 
% opts = optimset('lsqcurvefit');
% opts.Display='off';
% coefficients = lsqcurvefit(F,double(X0),double(X),double(Z),LB,UB,opts);
% Zout=F(coefficients,Xout);
% Zout=reshape(Zout,numel(xout),numel(yout));
% 
% sx=coefficients(1)*sd_max_x;
% sy=coefficients(2)*sd_max_y;
% phi=coefficients(3);
% zmax=coefficients(4)*zscaling*max(Z);
% xmax=coefficients(5)*x_range;
% ymax=coefficients(6)*y_range;

F = @(a,x) a(4)*zscaling*max(Z)*...
    exp(-((x(:,1)-a(5)*x_range)*cos(a(3))-(x(:,2)-a(6)*y_range)*sin(a(3))).^2/((2*a(1)*sd_max_x)^2)...
        -((x(:,1)-a(5)*x_range)*sin(a(3))+(x(:,2)-a(6)*y_range)*cos(a(3))).^2/((2*(a(1)*a(2)+(fitsettings.sd_x_min_ratio-a(1)*a(2))*(sign(fitsettings.sd_x_min_ratio-a(1)*a(2))+1)/2)*sd_max_x)^2));
%sigmax(ratio to sd_max_x)              ratio xy                            rotation peak   x0                                      y0
LB=[fitsettings.sd_x_min_ratio          fitsettings.sd_xy_min_ratio         -pi/2    0.5      min(min(xin))*range_factor/x_range     min(min(yin))*range_factor/y_range];
X0=[fitsettings.sd_x_min_ratio/2+1/2    fitsettings.sd_xy_min_ratio/2+1/2   0       1      start_pos_xy(1)*range_factor/x_range start_pos_xy(2)*range_factor/y_range];
UB=[1                                   1                                   pi/2    1.5      max(max(xin))*range_factor/x_range     max(max(yin))*range_factor/y_range]; 
% fitsettings.sd_xy_min_ratio=0.25;
% fitsettings.sd_xy_max_ratio=1;

opts = optimset('lsqcurvefit');
opts.Display='off';
% try
coefficients = lsqcurvefit(F,double(X0),double(X),double(Z),LB,UB,opts);
% catch eee
%     aaa=1;
% end
Zout=F(coefficients,Xout);
Zout=reshape(Zout,numel(xout),numel(yout));

sx=coefficients(1)*sd_max_x;
sy=sd_max_y*(coefficients(1)*coefficients(2) +... 
(fitsettings.sd_x_min_ratio-coefficients(1)*coefficients(2))*(sign(fitsettings.sd_x_min_ratio-coefficients(1)*coefficients(2))+1)/2);
phi=coefficients(3);
zmax=coefficients(4)*zscaling*max(Z);
xmax=coefficients(5)*x_range;
ymax=coefficients(6)*y_range;
end

function title_and_save(figure_handle,filename,plot_title,keys)
mtit(figure_handle,  plot_title, 'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none');
stampit;
if keys.create_pdfs
    switch keys.append_pdfs
        case 0
            export_fig([keys.path_to_save, filename], '-pdf','-transparent') % pdf by run
        case 1
            export_fig([keys.path_to_save,filename '_appended batches'], '-pdf', '-append','-transparent') % pdf by run
    end
    close all
end
end
