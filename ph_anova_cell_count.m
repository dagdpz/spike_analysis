function pie_data=ph_anova_cell_count(modified_keys)
for fn=fieldnames(modified_keys)'
    keys.(fn{:})=modified_keys.(fn{:});
end
[tuning_per_unit_table]=ph_load_extended_tuning_table(keys);
[tuning_per_unit_table, keys.selection_title]=ph_reduce_tuning_table(tuning_per_unit_table,keys);

%% here, make sure the epochs match currently checked epochs!
% 
% keys.ccs(cc).conditions_to_plot     ={'Msac'};
% keys.ccs(cc).plot_type              ='visuomotor';
% keys.ccs(cc).epochs.Msac            ={'INI', 'Fhol','Cue','MemE','MemL','PreS','PeriS','TIhol','Tons','Thol'}';

for t=1:numel(keys.ANOVAS_PER_TYPE)
    FNs=fieldnames(keys.ANOVAS_PER_TYPE);
    for f=1:numel(FNs)
        keys.ANOVAS_PER_TYPE(t).(FNs{f})=keys.CC.epochs;
    end
end
    tuning_per_unit_table=ph_multicomparison_correction(tuning_per_unit_table,keys);
%tuning_per_unit_table(cellfun(@(x) isempty(x) & ~islogical(x),tuning_per_unit_table))={''};

%% sort of obsolete
keys.tuning_table=tuning_per_unit_table;

% %% ANOVA criterions readout
% criterions={'space_or_interaction','epoch_or_interaction','hands_or_interaction','SXH_or_interaction'};
% for c=1:numel(criterions)
%     parameter_criterion_columns=~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),criterions{c}));
%     for t=1:numel(keys.CC.tasktypes)
%         task_criterion_columns(t,:)=~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),keys.CC.tasktypes{t}));
%     end
%     criterion_columns=parameter_criterion_columns & any(task_criterion_columns,1);
%     keys.(criterions{c})=any(cell2mat(tuning_per_unit_table(2:end,criterion_columns)),2);
% end

%% get legends and colors
keys=get_legends(keys);

%% Calculation in subfunctions
[pie_data matrix_data]=summary_anova_multilevel(keys);
n_cells=size(keys.tuning_table,1)-1;


%% PLOT: pie or bar

if keys.CC.plot_as_pie %% bar or pie
    keys.subfolder_to_save='cell_counts_as_pie';
else
    keys.subfolder_to_save='cell_counts_as_bar';
end
if keys.CC.percent
    keys.subfolder_to_save=[keys.subfolder_to_save '_percent'];
else
    keys.subfolder_to_save=[keys.subfolder_to_save '_absolute'];
end

figure('units','normalized','outerposition',[0 0 1 1],'color','w')
set(0, 'DefaultFigureRenderer', 'painters');
hold on
N_columns_rows=ceil(sqrt(size(pie_data,1)));
r_dec=1/size(pie_data,2);
for sub= 1:size(pie_data,1)
    subplot(N_columns_rows,N_columns_rows, sub);
    title(keys.all_titles{sub})
    hold on
    r=1;
    if keys.CC.plot_as_pie %% bar or pie
        levels_in_order=size(pie_data,2):-1:1;
    else
        levels_in_order=1:size(pie_data,2);
    end
    for level=levels_in_order
        if keys.CC.percent
            x=pie_data{sub,level}/max([n_cells 1])*100;
            ylim_bar=100;
        else
            x=pie_data{sub,level};
            ylim_bar=max([size(keys.tuning_table,1)-1 1]);
        end
        if strcmp(keys.CC.plot_type,'per_epoch') || strcmp(keys.CC.plot_type,'per_task') || strcmp(keys.CC.plot_type,'effector') || ...
                strcmp(keys.CC.plot_type,'visuomotor')
            current_colors=repmat(keys.all_colors{1},numel(x)/(size(keys.all_colors{1},1)),1);
        else
            current_colors=repmat(keys.all_colors{level},numel(x)/(size(keys.all_colors{level},1)),1);
        end
        x_as_labels=num2cell(x(:));
        if keys.CC.percent
            x_as_labels=cellfun(@(x) [num2str(round(x)) '%'],x_as_labels,'UniformOutput',false);
        else
            x_as_labels=cellfun(@(x) num2str(round(x)),x_as_labels,'UniformOutput',false);
        end
        if keys.CC.plot_as_pie  %% bar or pie
            hstepsize=2;
            piechart(level,sub).handle=DAG_pie_chart(x,r,current_colors,[0,0],(1-r_dec/2),x_as_labels);
            text(0,-r+r_dec/2,keys.x_labels{level});
            r=r-r_dec;
        else
            zero_dummie=zeros(size(pie_data,2)-level+1,numel(x));
            x_dummie=repmat(level:(size(pie_data,2)+1),numel(x),1)';
            if strcmp(keys.CC.plot_type,'per_epoch')  %|| strcmp(keys.CC.plot_type,'gaze')
                close(gcf)
                figure('units','normalized','outerposition',[0 0 1 1],'color','w')
                set(0, 'DefaultFigureRenderer', 'painters');
                colormap(current_colors);
                if keys.CC.percent
                    x=vertcat(pie_data{:})/max([n_cells 1])*100;
                    x_as_labels=cellfun(@(x) [num2str(round(x)) '%'],num2cell(x(:)),'UniformOutput',false);
                    %[num2str(round(x(:))) '%'];
                else
                    x=vertcat(pie_data{:});
                    x_as_labels=cellfun(@(x) num2str(round(x)),num2cell(x(:)),'UniformOutput',false);
                end
                cumsum_X=cumsum(x,2);
                cumsum_index=(ones(size(cumsum_X,2),1)*(1:size(cumsum_X,1)))';
                piechart(level,sub).handle=bar(x,'stacked');
                text(cumsum_index(:),cumsum_X(:)-x(:)/2,x_as_labels);
            else
                piechart(level,sub).handle=bar(x_dummie,[x;zero_dummie],1,'stacked');
                for k=1:numel(piechart(level,sub).handle)
                    set(piechart(level,sub).handle(k),'facecolor',current_colors(k,:))
                end
                cumsum_X=cumsum(x);
                cumsum_index=[cumsum_X(1)~=0, diff(cumsum_X)>0];
                text(repmat(level,sum(cumsum_index),1),cumsum_X(cumsum_index)-ylim_bar/25,x_as_labels(cumsum_index));
            end
        end
    end
    if ~keys.CC.plot_as_pie %% bar or pie
        ylabel('N units');
        hstepsize=1;
        if strcmp(keys.CC.plot_type,'per_epoch') %|| strcmp(keys.CC.plot_type,'gaze')
            set(gca,'xlim',[0.5,size(pie_data,1)+0.5],'xtick',[1:size(pie_data,1)],'xticklabel',keys.all_titles,'ylim',[0 ylim_bar]);
            break;
        else
            set(gca,'xlim',[0.5,size(pie_data,2)+0.5],'xtick',[1:size(pie_data,2)],'xticklabel',keys.x_labels,'ylim',[0 ylim_bar]);
        end
    end
end

%% legends
if keys.plot.cell_count_legends
    legend(piechart(1,1).handle(1+(0:keys.n1-1)*hstepsize),keys.legends{1},'location','best','fontsize',3);
    if keys.n2>0
        legend(piechart(2,2).handle(1+(0:keys.n2-1)*hstepsize),keys.legends{2},'location','best','fontsize',3);
    end
end
keys.title_part='complete';
title_and_save(keys);


%% plot similarity matrix
if ~(strcmp(keys.CC.plot_type,'per_epoch') || strcmp(keys.CC.plot_type,'per_task'))
    return;
end
figure('units','normalized','outerposition',[0 0 1 1],'color','w')
fields_to_plot= {'same','diff','none'};
for to_plot=1:numel(fields_to_plot)
    FN=fields_to_plot{to_plot};
    for sub= 1:numel(matrix_data)
        subplot(3,numel(matrix_data), sub+(to_plot-1)*numel(matrix_data));
        x_lim=size(matrix_data(sub).(FN),1);
        y_lim=size(matrix_data(sub).(FN),2);
        image(1:x_lim,1:y_lim,matrix_data(sub).(FN)./matrix_data(sub).base*62.5)
        title(keys.all_titles{sub})
        textpositions=combvec(1:x_lim,1:y_lim)';
        if keys.CC.percent
            Matrix_percent=round(matrix_data(sub).(FN)./repmat(diag(matrix_data(sub).same),1,size(matrix_data(sub).(FN),1))*100);
            text(textpositions(:,2),textpositions(:,1),num2str(Matrix_percent(:)),'horizontalalignment','center')
        else
            text(textpositions(:,2),textpositions(:,1),num2str(matrix_data(sub).(FN)(:)),'horizontalalignment','center')
        end
        set(gca,'xlim',[0.5 x_lim+0.5],'xtick',1:x_lim,'xticklabel',keys.x_labels,'ylim',[0.5 y_lim+0.5],'ytick',1:y_lim,'yticklabel',keys.x_labels);
        if sub==1
            ylabel(FN);
        end
        axis square
    end
end
keys.title_part='one_by_one';
title_and_save(keys);

function keys=get_legends(keys)
keys.SxH_mod_index=0;
cols=keys.colors;

%% color and legend assignment
epoch_legend={'en','bi','su','No tuning','No anova'};
epoch_colors=[cols.EP_EN; cols.EP_BI; cols.EP_SU; cols.NO_TU; cols.NO_AN]/255;
epochall_legend={'en','su','No tuning','No anova'};
epochall_colors=[cols.EP_EN; cols.EP_SU; cols.NO_TU; cols.NO_AN]/255;
choice_legend={'IN','CH','No tuning','No anova'};
choice_colors=[cols.in_IS; cols.ch_IS; cols.NO_TU; cols.NO_AN;]/255;
space_and_hand_legend={'IS','IHIS','IH','IHCS','CS','CHCS','CH','CHIS','incongruent','None'};
% space_and_hand_colors=[cols.in_IS; cols.in_IH_IS; cols.in_IH; cols.in_IH_CS;...
%     cols.in_CS; cols.in_CH_CS; cols.in_CH; cols.in_CH_IS; cols.NO_TU; cols.NO_AN]/255;


space_and_hand_colors=[cols.in_IS; cols.in_IH; cols.in_IH; cols.in_IH;...
    cols.in_CS; cols.in_CH; cols.in_CH; cols.in_CH; cols.NO_TU; cols.NO_AN]/255;

SH_as_enhancement_legend={'IS','IHIS','IH','IHCS','CS','CHCS','CH','CHIS','IHISsu','IHCSsu','CHCSsu','CHISsu','CR','UC','en','su','none','no anova'};
SH_as_enhancement_colors=[cols.in_IS; cols.in_IH_IS; cols.in_IH; cols.in_IH_CS;...
    cols.in_CS; cols.in_CH_CS; cols.in_CH; cols.in_CH_IS; cols.ch_IH_IS; cols.ch_IH_CS;cols.ch_CH_CS; cols.ch_CH_IS;...
    [0 122 222];[222 122 222];cols.EP_EN; cols.EP_SU;cols.NO_TU; cols.NO_AN]/255;
space_legend={'CS','IS','No tuning','No anova'};
space_colors=[cols.in_CS; cols.in_IS; cols.NO_TU; cols.NO_AN]/255;
hand_legend={'CH','IH','No tuning','No anova'};
hand_colors=[cols.in_CH; cols.in_IH; cols.NO_TU; cols.NO_AN]/255;
effector_legend={'EF1','EF2','No tuning','No anova'};
effector_colors=[cols.EP_EN; cols.EP_SU; cols.NO_TU; cols.NO_AN]/255;
inactivation_legend={'CS:EN','CS:EN & IS:EN','IS:EN','CS:SU & IS:EN','CS:SU','CS:SU & IS:SU','IS:SU','CS:EN & IS:SU', 'No tuning'};
inactivation_colors=[255 0 0; 255 77 10; 255 153 20; 255 255 0; 0 255 178; 0 178 166; 0 100 255;128 50 205; cols.NO_TU]/255;
%     keys.all_colors{1}   =[cols.in_CS; (cols.in_CS+cols.in_IS)/2; cols.in_IS; (255-cols.in_CS+cols.in_IS)/2;...
%         255-cols.in_CS; (510-cols.in_CS-cols.in_IS)/2; 255-cols.in_IS;(255+cols.in_CS-cols.in_IS)/2; cols.NO_TU]/255;
space_position_legend={'CS','IS','Position','No tuning'};
space_position_colors=[cols.in_CS; cols.in_IS; cols.EP_BI; cols.NO_TU]/255;
crossed_uncrossed_legend={'UC','CR','No tuning','No anova'};
crossed_uncrossed_colors=[255 0 0; 125 0 0; cols.NO_TU; cols.NO_AN]/255;
visuomotor_legend={'Vis','Vmot','Mot','No anova'};
visuomotor_colors=[cols.in_CS; cols.EP_BI; cols.in_IS; cols.NO_AN]/255;
fixation_x_position_comb_legend={'Position','Gaze','G+P','Any GxP','None'};
fixation_x_position_comb_colors=[247 148 30; 238 43 123; 159 31 99; 199 45 196; 255 255 255]/255;
fixation_x_position_CI_legend={'Position','Gaze','G+P','None'};
fixation_x_position_CI_colors=[247 148 30; 238 43 123; 159 31 99; 255 255 255]/255;
eccentricity_x_angle_legend={'All','Ecc + Ang','Ecc + ExA','Ecc','Ang + ExA','Ang','ExA','None'};
eccentricity_x_angle_colors=[222 220 0; 255 150 0; 255 0 0; 255 0 150; 171 0 252; 0 0 255; 125 130 255; 255 255 255]/255;
gaze_legend={'Tar','Tar+Gaze','Gaze Fhol','None'};
gaze_colors=[127 127 0; 0 255 0;0 127 127;255 255 255]/255;
gaze_interaction_colors=[0 255 0; 255 255 255]/255;
gaze_interaction_legend={'Interaction','None'}; %% CHECK
gaze_and_fixation_x_position_legend={'Position','Gaze','G+P','G+P+GxP','P+GxP','G+GxP','GxP','None'};
gaze_and_fixation_x_position_colors=[222 220 0; 255 150 0; 255 0 0; 255 0 150; 171 0 252; 0 0 255; 125 130 255; 255 255 255]/255;
fixation_x_position_legend={'Position','Gaze','G+P','G+P+GxP','P+GxP','G+GxP','GxP','None'};
fixation_x_position_colors=[222 220 0; 255 150 0; 255 0 0; 255 0 150; 171 0 252; 0 0 255; 125 130 255; 255 255 255]/255; % not sure if this is correct
gaze_fhol_color=[0 255 0; 255 255 255]/255;
gaze_fhol_legend={'None','Gaze Fhol'}; %% CHECK
bilateral_legend={'en','su','No tuning','No anova'}; %CHECK (really)
bilateral_colors=[cols.EP_EN; cols.EP_SU; cols.NO_AN; cols.NO_TU]/255;

keys.legends={};
factors{1}=keys.CC.factor;

keys.all_titles=keys.CC.epochs;
switch keys.CC.plot_type
    case 'per_epoch_2levels' %% new case here
        factors{2}= factors{1};
        factors{1}=keys.CC.first_level_factor;
        keys.x_labels=factors;
    case 'per_epoch'
        keys.all_titles=keys.CC.epochs;
        keys.x_labels=keys.CC.tasktypes;
        [factors{1:numel(keys.CC.tasktypes)}]=deal(keys.CC.factor);
    case 'per_task'
        keys.all_titles=keys.CC.tasktypes;
        keys.x_labels=keys.CC.epochs;
        [factors{1:numel(keys.CC.epochs)}]=deal(keys.CC.factor);
    case 'effector'
        keys.x_labels=keys.CC.epochs;
        factors{1}='effector';
    case 'hemi_and_epoch'
        keys.x_labels={'space','epoch'};
        factors{1}='epoch';
        factors{2}='space';
    case 'space_and_epoch'
        keys.x_labels={'space','epoch'};
        factors{1}='epoch';
        factors{2}='space_position';
    case 'space_and_bilateral' %???
        keys.x_labels={'space','bilateral'};
        factors{1}='space';
        factors{2}='bilateral';
    case 'space_and_hand'
        keys.x_labels={'space','hand'};
        factors{1}='space';
        factors{2}='hand';
    case 'hand_and_epoch' % reverse order
        keys.x_labels={'epoch','hand'};
        factors{1}='epoch';
        factors{2}='hand';
    case 'space_x_hand'
        keys.x_labels={'Main effects','Interactions'};
        factors{1}='space_and_hand';
        factors{2}='space_x_hand';
    case 'hands_inactivation'
        keys.x_labels={'Contra hand','Ipsi hand'};
        factors{1}='hands_inactivation';
        factors{2}='hands_inactivation';
    case 'fixation_x_position' %% something strange
        keys.x_labels={'Retinotopic','initial gaze'}; %%?
        factors{1}='fixation_x_position';
        factors{2}='gaze_fhol';
    case 'fixation_x_position_comb' %% something strange
        keys.x_labels={''}; %%?
        factors{1}='fixation_x_position_comb';
    case 'fixation_x_position_CI' %% something strange
        keys.x_labels={''}; %%?
        factors{1}='fixation_x_position_CI';
    case 'gaze_and_fixation_x_position'
        keys.x_labels={'Retinotopic','Object centered'};
        factors{1}='gaze_and_fixation_x_position';
        factors{2}='gaze';
    case 'gaze'
        keys.x_labels={'Retinotopic','Object centered'};
        factors{1}='gaze';
        factors{2}='gaze_interaction';
    case 'eccentricity_x_angle'
        keys.x_labels={''};
        factors{1}='eccentricity_x_angle';
    case 'visuomotor'
        keys.x_labels={''}; %%?
        keys.all_titles={'category'};
        factors{1}='visuomotor';
end

levels=1:numel(factors); % numel(legends)?
switch keys.CC.plot_type %its really only one level (of complexity) here
    case {'per_epoch','per_task'}
        levels=1;
end

keys.all_titles=keys.CC.epochs;
keys.factors=factors;
for l=levels
    factor=factors{l};
    switch factor
        case 'hand'
            keys.all_colors{l}   =hand_colors;
            keys.legends{l}=hand_legend;
        case 'space'
            keys.all_colors{l}   =space_colors;
            keys.legends{l}=space_legend;
        case {'space_and_hand','space_hand'}
            keys.all_colors{l}   =space_and_hand_colors;
            keys.legends{l}=space_and_hand_legend;
        case 'SH_as_enhancement'
            keys.all_colors{l}   =SH_as_enhancement_colors;
            keys.legends{l}=SH_as_enhancement_legend;
        case {'space_position','position_space'}
            keys.all_colors{l}   =space_position_colors;
            keys.legends{l}=space_position_legend;
        case {'choice1','choice2'}
            keys.all_colors{l}   =choice_colors;
            keys.legends{l}=choice_legend;
        case 'epoch'
            keys.all_colors{l}   =epoch_colors;
            keys.legends{l}=epoch_legend;
        case 'effector';
            keys.all_colors{l}   =effector_colors;
            keys.legends{l}=effector_legend;
        case 'epoch_all'
            keys.all_colors{l}   =epochall_colors;
            keys.legends{l}=epochall_legend;
        case 'fixation_x_position'
            keys.all_colors{l}   =fixation_x_position_colors;
            keys.legends{l}=fixation_x_position_legend;
        case 'visuomotor'
            keys.all_colors{l}   =visuomotor_colors;
            keys.legends{l}=visuomotor_legend;
            %         case 'SxH'
            %             keys.all_colors{l}   =SxH_colors;
            %             keys.legends{l}=SxH_legend;
        case 'space_x_hand'
            keys.all_colors{l}   =crossed_uncrossed_colors;
            keys.legends{l}=crossed_uncrossed_legend;
        case 'hands_inactivation'
            keys.all_colors{l}   =inactivation_colors;
            keys.legends{l}=inactivation_legend;
        case 'fixation_x_position_comb'
            keys.all_colors{l}   =fixation_x_position_comb_colors;
            keys.legends{l}=fixation_x_position_comb_legend;
        case 'fixation_x_position_CI'
            keys.all_colors{l}   =fixation_x_position_CI_colors;
            keys.legends{l}=fixation_x_position_CI_legend;
        case 'gaze_fhol';
            keys.all_colors{l}   =gaze_fhol_color;
            keys.legends{l}=gaze_fhol_legend;
        case 'eccentricity_x_angle'
            keys.all_colors{l}   =eccentricity_x_angle_colors;
            keys.legends{l}=eccentricity_x_angle_legend;
        case 'gaze_and_fixation_x_position';
            keys.all_colors{l}   =gaze_and_fixation_x_position_colors;
            keys.legends{l}=gaze_and_fixation_x_position_legend;
        case 'gaze';
            keys.all_colors{l}   =gaze_colors;
            keys.legends{l}=gaze_legend;
        case 'gaze_interaction';
            keys.all_colors{l}   =gaze_interaction_colors;
            keys.legends{l}=gaze_interaction_legend;
        case 'bilateral';
            keys.all_colors{l}   =bilateral_colors;
            keys.legends{l}=bilateral_legend;
    end
end

keys.n1=numel(keys.legends{1});
if size(factors)>1
    keys.n2=numel(keys.legends{2});
else
    keys.n2=0;
end

function [multilevel_data matrix_data] = summary_anova_multilevel(keys)
xlsx_table=keys.tuning_table;
for i=1:size(keys.tuning_table,2)
    param(i) = keys.tuning_table(1,i);
    idx.(param{i})= DAG_find_column_index(xlsx_table,param{i});
end

multilevel_data={};
matrix_data=[];
no_anova_marker=keys.n1;
no_tuning_marker=keys.n1-1;

%% different loop depending on plot_type
first_loop=1:numel(keys.CC.epochs);
second_loop=1:numel(keys.factors);

switch keys.CC.plot_type
    case 'per_task'
        first_loop=1:numel(keys.CC.tasktypes);
        second_loop=1:numel(keys.CC.epochs);
    case 'per_epoch'
        second_loop=1:numel(keys.CC.tasktypes);
end

for e=first_loop
    for L=second_loop
        
        epoch             =   keys.CC.epochs{e};
        factor            =   keys.factors{L};
        tasktype          =   keys.CC.tasktypes{1};
        
        IC_to_plot        =   keys.CC.IC_to_plot;
        arrangement       =   keys.arrangement(1:3);
        hand_labels{1}    =   '_';
        parameters        =   {};
        epochs            =   {};
        arrangements      =   {};
        
        switch keys.CC.plot_type
            case 'per_task'
                tasktype          =   keys.CC.tasktypes{e};
                epoch             =   keys.CC.epochs{L};
            case 'per_epoch'
                tasktype          =   keys.CC.tasktypes{L};
        end
        
        %% level dependent
        if L==1
            keys.n=keys.n1;
            if strcmp(keys.CC.plot_type,'per_epoch_2levels')
                epoch             =   keys.CC.first_level_epochs{1};
                factor            =   keys.CC.first_level_factor;
            end
        elseif L==2
            keys.n=keys.n2;
        end
        
        switch factor
            
            case 'epoch'
                hand_labels{1}='_AH_';
                parameters{1}='_epoch_';
            case 'epoch_all'
                parameters{1}='_epoch_';
            case 'bilateral'
                parameters{1}='_epoch_bilateral_';
                %                 %%% why this is so confusing
                %                 keys.CC.factor='space';
                %                 keys.CC.factor='epoch';
                hand_labels{1}=['_' keys.labels.handsIC{keys.tt.hands(1)+1} '_'];
            case 'space'
                parameters{1}= '_hemifield_';
            case 'space_and_hand'
                parameters{1}='_hemifield_';
                parameters{2}='_hands_';
                hand_labels{2}='_';
            case 'space_hand'
                hand_labels{1}='_CH_';
                hand_labels{2}='_IH_';
                hand_labels{3}='_CS_';
                hand_labels{4}='_IS_';
                parameters{1}='_hemifield_';
                parameters{2}='_hemifield_';
                parameters{3}='_hands_';
                parameters{4}='_hands_';
                
            case 'space_position'
                parameters{1}='_hemifield_';
                parameters{2}='_position_';
                hand_labels{2}=['_' keys.labels.handsIC{keys.tt.hands(1)+1} '_'];
            case 'position_space'
                parameters{1}='_hemifield_';
                parameters{2}='_position_';
                hand_labels{2}=['_' keys.labels.handsIC{keys.tt.hands(1)+1} '_'];
            case 'hand'
                parameters{1}='_hands_';
            case 'choice1'
                hand_labels{1}='_AH_';
                parameters{1}='_RF_choice1_';
            case 'effector'
                hand_labels{1}=keys.CC.factor;
                parameters{1}=['_effector_' keys.CC.tasktypes{1} '_vs_'];
                tasktype=keys.CC.tasktypes{2};
            case 'space_x_hand'
                parameters{1}='_SxH_';
            case 'SH_as_enhancement'
                hand_labels{1}='_CH_CS_';
                hand_labels{2}='_CH_IS_';
                hand_labels{3}='_IH_CS_';
                hand_labels{4}='_IH_IS_';
                parameters{1}='_epoch_';
                parameters{2}='_epoch_';
                parameters{3}='_epoch_';
                parameters{4}='_epoch_';
            case 'hands_inactivation'
                parameters{1}='_PT_';
                parameters{2}='_PT_';
                if level==1
                    hand_labels{1}='_CH_CS_';
                    hand_labels{2}='_CH_IS_';
                elseif level==2
                    hand_labels{1}='_IH_CS_';
                    hand_labels{2}='_IH_IS_';
                end
            case  'eccentricity_x_angle';
                hand_labels{1}='_AH_';
                hand_labels{2}='_AH_';
                hand_labels{3}='_AH_';
                parameters{1}='_distance_';
                parameters{2}='_angle_';
                parameters{3}='_DxA_';
            case 'fixation_x_position_comb'
                hand_labels{1}='_AH_';
                hand_labels{2}='_AH_';
                hand_labels{3}='_AH_';
                parameters{1}='_fixation_';
                parameters{2}='_position_';
                parameters{3}='_PxF_';
            case 'fixation_x_position_CI'
                parameters{1}='_hemifield_';
                parameters{2}='_hemifield_';
                hand_labels{2}='_';
                arrangements{1}='fix';
                arrangements{2}='mov';
            case 'fixation_x_position'
                hand_labels{1}='_AH_';
                hand_labels{2}='_AH_';
                hand_labels{3}='_AH_';
                parameters{1}='_fixation_';
                parameters{2}='_position_';
                parameters{3}='_PxF_';
            case 'gaze_fhol'
                %% ?? keys.CC.factor='fixation_x_position';
                epochs{1}='Fhol';
                parameters{1}='_position_';
                arrangements{1}='fix';
            case 'gaze'
                hand_labels{1}='_AH_';
                hand_labels{2}='_AH_';
                epochs{1}='Fhol';
                epochs{2}=epoch;
                arrangements{1}='fix';
                arrangements{2}='tar';
                parameters{1}='_position_';
                parameters{2}='_position_';
            case 'gaze_interaction'
                parameters{1}='_PxF_';
                arrangement='tar';
            case 'visuomotor'
                tuning_variables{1}= ['visual_'   tasktype '_' arrangement];
                tuning_variables{2}= ['motor_'   tasktype '_' arrangement];
                tuning_variables{3}= ['visuomotor_'   tasktype '_' arrangement];
        end
        
        for tv=1:numel(parameters)
            if ~isempty(epochs)
                epoch=epochs{tv};
            end
            if ~isempty(arrangements)
                arrangement=arrangements{tv};
            end
            tuning_variables{tv}= [IC_to_plot hand_labels{tv} epoch parameters{tv} tasktype '_' arrangement];
        end
        keys.CC.factor=factor;
        tuning{L}=read_table_column_detail(keys,xlsx_table,idx,tuning_variables);
    end
    % KK
    pie_tmp=[];
    pie_tmp=get_pie_multilevel(keys,pie_tmp,tuning);
    multilevel_data(e,:)=pie_tmp;
    
    for d1=second_loop
        basic_quantity=tuning{d1}~=no_tuning_marker & tuning{d1}~=no_anova_marker;
        for d2=second_loop
            matrix_data(e).same(d1,d2)=sum(tuning{d1}(basic_quantity)==tuning{d2}(basic_quantity));
            matrix_data(e).diff(d1,d2)=sum(tuning{d1}(basic_quantity)~=tuning{d2}(basic_quantity) & tuning{d2}(basic_quantity)~=no_tuning_marker);
            matrix_data(e).none(d1,d2)=sum(tuning{d2}(basic_quantity)==no_tuning_marker);
            matrix_data(e).base(d1,d2)=sum(basic_quantity);
        end
    end
end

function read_out=read_table_column_detail(keys,table,index,tuning_variables)
%% first check if entries are valid
valid=true;
for k=1:size(tuning_variables,2)
    valid=valid && isfield(index,tuning_variables{k}) && any(index.(tuning_variables{k}));
end

if ~ valid
    read_out=ones(1,size(table,1)-1)*keys.n;
    return;
end

switch keys.CC.factor
    case {'CH_CS','CH_IS','IH_CS','IH_IS'}
        EF1=ismember(table(2:end,index.(tuning_variables{1})),keys.CC.tasktypes{1});
        EF2=ismember(table(2:end,index.(tuning_variables{1})),keys.CC.tasktypes{2});
        NA=ismember(table(2:end,index.(tuning_variables{1})),'-');
        
        read_out(EF1)                =1;
        read_out(EF2)                =2;
        read_out(~EF1 & ~EF2  &  ~NA)=3;
        read_out(NA)                 =4;
        
    case 'space_hand'
        CH_IS=ismember(table(2:end,index.(tuning_variables{1})),'IS');
        CH_CS=ismember(table(2:end,index.(tuning_variables{1})),'CS');
        IH_IS=ismember(table(2:end,index.(tuning_variables{2})),'IS');
        IH_CS=ismember(table(2:end,index.(tuning_variables{2})),'CS');
        CS_IH=ismember(table(2:end,index.(tuning_variables{3})),'IH');
        CS_CH=ismember(table(2:end,index.(tuning_variables{3})),'CH');
        IS_IH=ismember(table(2:end,index.(tuning_variables{4})),'IH');
        IS_CH=ismember(table(2:end,index.(tuning_variables{4})),'CH');
        
        IS=CH_IS |IH_IS;
        CS=CH_CS |IH_CS;
        IH=IS_IH |CS_IH;
        CH=CS_CH |IS_CH;
        
        IHIS=(IH & IS);
        IHCS=(IH & CS);
        CHIS=(CH & IS);
        CHCS=(CH & CS);
        
        incongruent= (IH_IS & CH_CS) | (CH_IS & IH_CS) | (CS_CH & IS_IH) | (CS_IH & IS_CH);
        na=~CH_IS &~CH_CS &~IH_IS &~IH_CS &~CS_IH &~CS_CH &~IS_IH &~IS_CH;
        %rest=~IS &~CS &~IH &~CH &~IHIS &~IHCS &~CHIS &~CHCS&~na&~incongruent;
        
        read_out(IS)            =mod(0+keys.SxH_mod_index,8)+1;
        read_out(IHIS)          =mod(1+keys.SxH_mod_index,8)+1;
        read_out(IH)            =mod(2+keys.SxH_mod_index,8)+1;
        read_out(IHCS)          =mod(3+keys.SxH_mod_index,8)+1;
        read_out(CS)            =mod(4+keys.SxH_mod_index,8)+1;
        read_out(CHCS)          =mod(5+keys.SxH_mod_index,8)+1;
        read_out(CH)            =mod(6+keys.SxH_mod_index,8)+1;
        read_out(CHIS)          =mod(7+keys.SxH_mod_index,8)+1;
        read_out(incongruent)   =9;
        read_out(na)            =10;
        
    case 'space_and_hand'
        IStuning=ismember(table(2:end,index.(tuning_variables{1})),'IS');
        CStuning=ismember(table(2:end,index.(tuning_variables{1})),'CS');
        IHtuning=ismember(table(2:end,index.(tuning_variables{2})),'IH');
        CHtuning=ismember(table(2:end,index.(tuning_variables{2})),'CH');
        IHIS=(IHtuning & IStuning);
        IHCS=(IHtuning & CStuning);
        CHIS=(CHtuning & IStuning);
        CHCS=(CHtuning & CStuning);
        IH=IHtuning & ~IStuning & ~CStuning;
        CH=CHtuning & ~IStuning & ~CStuning;
        IS=IStuning & ~IHtuning & ~CHtuning;
        CS=CStuning & ~IHtuning & ~CHtuning;
        read_out(IS)            =mod(0+keys.SxH_mod_index,8)+1;
        read_out(IHIS)          =mod(1+keys.SxH_mod_index,8)+1;
        read_out(IH)            =mod(2+keys.SxH_mod_index,8)+1;
        read_out(IHCS)          =mod(3+keys.SxH_mod_index,8)+1;
        read_out(CS)            =mod(4+keys.SxH_mod_index,8)+1;
        read_out(CHCS)          =mod(5+keys.SxH_mod_index,8)+1;
        read_out(CH)            =mod(6+keys.SxH_mod_index,8)+1;
        read_out(CHIS)          =mod(7+keys.SxH_mod_index,8)+1;
        read_out(~IHtuning & ~CHtuning & ~IStuning & ~CStuning)=9;
        %read_out(na)            =10;
        
    case 'SH_as_enhancement'
        CS_CH_en=ismember(table(2:end,index.(tuning_variables{1})),'en');
        IS_CH_en=ismember(table(2:end,index.(tuning_variables{2})),'en');
        CS_IH_en=ismember(table(2:end,index.(tuning_variables{3})),'en');
        IS_IH_en=ismember(table(2:end,index.(tuning_variables{4})),'en');
        CS_CH_su=ismember(table(2:end,index.(tuning_variables{1})),'su');
        IS_CH_su=ismember(table(2:end,index.(tuning_variables{2})),'su');
        CS_IH_su=ismember(table(2:end,index.(tuning_variables{3})),'su');
        IS_IH_su=ismember(table(2:end,index.(tuning_variables{4})),'su');
        
        IHIS=( IS_IH_en &~IS_CH_en &~CS_IH_en &~CS_CH_en) | (~IS_IH_su & IS_CH_su & CS_IH_su & CS_CH_su);
        IHCS=(~IS_IH_en &~IS_CH_en & CS_IH_en &~CS_CH_en) | ( IS_IH_su & IS_CH_su &~CS_IH_su & CS_CH_su);
        CHIS=(~IS_IH_en & IS_CH_en &~CS_IH_en &~CS_CH_en) | ( IS_IH_su &~IS_CH_su & CS_IH_su & CS_CH_su);
        CHCS=(~IS_IH_en &~IS_CH_en &~CS_IH_en & CS_CH_en) | ( IS_IH_su & IS_CH_su & CS_IH_su &~CS_CH_su);
        
        
        IHISsu=(~IS_IH_en & IS_CH_en & CS_IH_en & CS_CH_en) | ( IS_IH_su &~IS_CH_su &~CS_IH_su &~CS_CH_su);
        IHCSsu=( IS_IH_en & IS_CH_en &~CS_IH_en & CS_CH_en) | (~IS_IH_su &~IS_CH_su & CS_IH_su &~CS_CH_su);
        CHISsu=( IS_IH_en &~IS_CH_en & CS_IH_en & CS_CH_en) | (~IS_IH_su & IS_CH_su &~CS_IH_su &~CS_CH_su);
        CHCSsu=( IS_IH_en & IS_CH_en & CS_IH_en &~CS_CH_en) | (~IS_IH_su &~IS_CH_su &~CS_IH_su & CS_CH_su);
        
        IS=( IS_IH_en & IS_CH_en &~CS_IH_en &~CS_CH_en) | (~IS_IH_su &~IS_CH_su & CS_IH_su & CS_CH_su);
        CS=(~IS_IH_en &~IS_CH_en & CS_IH_en & CS_CH_en) | ( IS_IH_su & IS_CH_su &~CS_IH_su &~CS_CH_su);
        CH=(~IS_IH_en & IS_CH_en &~CS_IH_en & CS_CH_en) | ( IS_IH_su &~IS_CH_su & CS_IH_su &~CS_CH_su);
        IH=( IS_IH_en &~IS_CH_en & CS_IH_en &~CS_CH_en) | (~IS_IH_su & IS_CH_su &~CS_IH_su & CS_CH_su);
        
        CR=( IS_CH_en & CS_IH_en &~IS_IH_en &~CS_CH_en) | (~IS_CH_su &~CS_IH_su & IS_IH_su & CS_CH_su);
        UN=(~IS_CH_en &~CS_IH_en & IS_IH_en & CS_CH_en) | ( IS_CH_su & CS_IH_su &~IS_IH_su &~CS_CH_su);
        
        EN= IS_CH_en & CS_IH_en & IS_IH_en & CS_CH_en;
        SU= IS_CH_su & CS_IH_su & IS_IH_su & CS_CH_su;
        
        na=~IS_IH_en &~IS_CH_en &~CS_IH_en &~CS_CH_en &~IS_IH_su &~IS_CH_su &~CS_IH_su &~CS_CH_su;
        rest=~IHIS & ~IHCS & ~CHIS & ~CHCS & ~IS & ~CS & ~CH & ~IH & ~IHISsu & ~IHCSsu & ~CHISsu & ~CHCSsu & ~EN & ~SU & ~CR & ~UN & ~na;
        
        read_out(IS)            =mod(0+keys.SxH_mod_index,8)+1;
        read_out(IHIS)          =mod(1+keys.SxH_mod_index,8)+1;
        read_out(IH)            =mod(2+keys.SxH_mod_index,8)+1;
        read_out(IHCS)          =mod(3+keys.SxH_mod_index,8)+1;
        read_out(CS)            =mod(4+keys.SxH_mod_index,8)+1;
        read_out(CHCS)          =mod(5+keys.SxH_mod_index,8)+1;
        read_out(CH)            =mod(6+keys.SxH_mod_index,8)+1;
        read_out(CHIS)          =mod(7+keys.SxH_mod_index,8)+1;
        read_out(IHISsu)          =9;
        read_out(IHCSsu)          =10;
        read_out(CHCSsu)          =11;
        read_out(CHISsu)          =12;
        read_out(CR)            =13;
        read_out(UN)            =14;
        read_out(EN)            =15;
        read_out(SU)            =16;
        read_out(rest)          =17;
        read_out(na)            =18;
        
    case 'hand'
        IH=ismember(table(2:end,index.(tuning_variables{1})),'IH') ;
        CH=ismember(table(2:end,index.(tuning_variables{1})),'CH') ;
        read_out(CH)                =1;
        read_out(IH)                =2;
        read_out(~IH & ~CH)   =3;
        
    case 'space'
        IS=ismember(table(2:end,index.(tuning_variables{1})),'IS');
        CS=ismember(table(2:end,index.(tuning_variables{1})),'CS');
        read_out(CS)                =1;
        read_out(IS)                =2;
        read_out(~IS & ~CS)   =3;
        
    case 'epoch'
        en=ismember(table(2:end,index.(tuning_variables{1})),'en');
        su=ismember(table(2:end,index.(tuning_variables{1})),'su');
        bi=ismember(table(2:end,index.(tuning_variables{1})),'bi');
        read_out(en)                =1;
        read_out(bi)                =2;
        read_out(su)                =3;
        read_out(~en & ~su & ~bi)   =4;
        
    case 'epoch_all'
        en=ismember(table(2:end,index.(tuning_variables{1})),'en');
        su=ismember(table(2:end,index.(tuning_variables{1})),'su');
        read_out(en)                =1;
        read_out(su)                =2;
        read_out(~en & ~su)         =3;
        
    case 'choice1'
        IN=ismember(table(2:end,index.(tuning_variables{1})),'IN');
        CH=ismember(table(2:end,index.(tuning_variables{1})),'CH');
        read_out(IN)                =1;
        read_out(CH)                =2;
        read_out(~IN & ~CH)         =3;
        
    case 'space_x_hand'
        CR=ismember(table(2:end,index.(tuning_variables{1})),'CR');
        UC=ismember(table(2:end,index.(tuning_variables{1})),'UC');
        read_out(UC)                =1;
        read_out(CR)                =2;
        read_out(~CR & ~UC)         =3;
        
    case 'fixation_x_position'
        fixtuning=ismember(table(2:end,index.(tuning_variables{1})),'YE');
        postuning=ismember(table(2:end,index.(tuning_variables{2})),'YE');
        fxptuning=ismember(table(2:end,index.(tuning_variables{3})),'YE');
        Pos= ~fixtuning  & postuning  & ~fxptuning;
        Fix=  fixtuning  & ~postuning & ~fxptuning;
        FnP=  fixtuning  &  postuning & ~fxptuning;
        All=  fixtuning  &  postuning &  fxptuning;
        PnX= ~fixtuning  &  postuning &  fxptuning;
        FnX=  fixtuning  & ~postuning &  fxptuning;
        FxP= ~fixtuning  & ~postuning &  fxptuning;
        Na = ~fixtuning  & ~postuning & ~fxptuning;
        
        read_out(Pos)          =1;
        read_out(Fix)          =2;
        read_out(FnP)          =3;
        read_out(All)          =4;
        read_out(PnX)          =5;
        read_out(FnX)          =6;
        read_out(FxP)          =7;
        read_out(Na)           =8;
        
    case 'eccentricity_x_angle'
        ecctuning=ismember(table(2:end,index.(tuning_variables{1})),'YE');
        angtuning=ismember(table(2:end,index.(tuning_variables{2})),'YE');
        exatuning=ismember(table(2:end,index.(tuning_variables{3})),'YE');
        All=  ecctuning  &  angtuning &  exatuning;
        EnA=  ecctuning  &  angtuning & ~exatuning;
        EnX=  ecctuning  & ~angtuning &  exatuning;
        Ecc=  ecctuning  & ~angtuning & ~exatuning;
        AnX= ~ecctuning  &  angtuning &  exatuning;
        Ang= ~ecctuning  &  angtuning & ~exatuning;
        ExA= ~ecctuning  & ~angtuning &  exatuning;
        Na=  ~ecctuning  & ~angtuning & ~exatuning;
        
        read_out(All)          =1;
        read_out(EnA)          =2;
        read_out(EnX)          =3;
        read_out(Ecc)          =4;
        read_out(AnX)          =5;
        read_out(Ang)          =6;
        read_out(ExA)          =7;
        read_out(Na)           =8;
        
    case 'fixation_x_position_CI'
        fixtuning=ismember(table(2:end,index.(tuning_variables{1})),'IS') | ismember(table(2:end,index.(tuning_variables{1})),'CS');
        postuning=ismember(table(2:end,index.(tuning_variables{2})),'IS') | ismember(table(2:end,index.(tuning_variables{2})),'CS');
        %CStuning=ismember(table(2:end,index.(indexfnS)),'CS');
        Pos= ~fixtuning  & postuning  ;
        Fix=  fixtuning  & ~postuning ;
        FnP=  fixtuning  &  postuning ;
        Na = ~fixtuning  & ~postuning ;
        
        read_out(Pos)          =1;
        read_out(Fix)          =2;
        read_out(FnP)          =3;
        read_out(Na)           =4;
        
    case 'gaze'
        fixtuning=ismember(table(2:end,index.(tuning_variables{1})),'YE');
        tartuning=ismember(table(2:end,index.(tuning_variables{2})),'YE');
        
        Tar= ~fixtuning  &  tartuning;
        Fix=  fixtuning  & ~tartuning;
        FnT=  fixtuning  &  tartuning;
        Na = ~fixtuning  & ~tartuning;
        read_out(Tar)          =1;
        read_out(FnT)          =2;
        read_out(Fix)          =3;
        read_out(Na)           =4;
        
    case 'gaze_interaction';
        interaction=ismember(table(2:end,index.(tuning_variables{1})),'YE');
        read_out(interaction)          =1;
        read_out(~interaction)         =2;
        
    case 'visuomotor'
        vis=cell2mat(table(2:end,index.(tuning_variables{1})));
        mot=cell2mat(table(2:end,index.(tuning_variables{2})));
        vmt=cell2mat(table(2:end,index.(tuning_variables{3})));
        read_out(vis)                =1;
        read_out(vmt)                =2;
        read_out(mot)                =3;
        read_out(~vis & ~mot & ~vmt) =4;
        
    case 'hands_inactivation'
        CS_EN=ismember(table(2:end,index.(tuning_variables{1})),'EN');
        CS_SU=ismember(table(2:end,index.(tuning_variables{1})),'SU');
        IS_EN=ismember(table(2:end,index.(tuning_variables{2})),'EN');
        IS_SU=ismember(table(2:end,index.(tuning_variables{2})),'SU');
        
        read_out(CS_EN & ~(IS_EN | IS_SU))              =1;
        read_out(CS_EN & IS_EN)                         =2;
        read_out(IS_EN & ~(CS_EN | CS_SU))              =3;
        read_out(CS_SU & IS_EN)                         =4;
        read_out(CS_SU & ~(IS_EN | IS_SU))              =5;
        read_out(CS_SU & IS_SU)                         =6;
        read_out(IS_SU & ~(CS_EN | CS_SU))              =7;
        read_out(CS_EN & IS_SU)                         =8;
        read_out(~CS_EN & ~CS_SU & ~IS_SU & ~IS_EN)     =9;
        
    case 'space_position'
        CS=ismember(table(2:end,index.(tuning_variables{1})),'CS');
        IS=ismember(table(2:end,index.(tuning_variables{1})),'IS');
        PO=ismember(table(2:end,index.(tuning_variables{2})),'YE');
        clear read_out
        read_out(CS)                =1;
        read_out(IS)                =2;
        read_out(~IS & ~CS & PO)    =3;
        read_out(~IS & ~CS & ~PO)   =4;
        
    case 'position_space'
        CS=ismember(table(2:end,index.(tuning_variables{1})),'CS');
        IS=ismember(table(2:end,index.(tuning_variables{1})),'IS');
        PO=ismember(table(2:end,index.(tuning_variables{2})),'YE');
        clear read_out
        read_out(CS & PO)           =1;
        read_out(IS & PO)           =2;
        read_out(~IS & ~CS & PO)    =3;
        read_out(~PO)               =4;
        
end

function pie=get_pie_multilevel(keys,pie,tuning,previous_level_index,level)
if nargin<5
    level=0;
end
level=level+1;
if nargin<4
    previous_level_index=true(size(tuning{level}));
end
if numel(pie)<level
    pie{level}=[];
end
n_conditions=keys.n1;
for c=1:n_conditions
    next_level_index = tuning{level}==c;
    pie{level}(end+1)= sum(previous_level_index(:) & next_level_index(:));
    if level<numel(tuning)
        pie=get_pie_multilevel(keys,pie,tuning,previous_level_index(:) & next_level_index(:),level);
    end
end

function title_and_save(keys)
mtit(gcf,[keys.monkey ' ' [keys.CC.tasktypes{:}] ' ' [keys.CC.plot_type] ' ' keys.CC.factor ' ' keys.CC.IC_to_plot ' ' keys.arrangement(1:3) ' ' keys.selection_title{:} ' N: ' num2str(size(keys.tuning_table,1)-1) ' ' keys.title_part], 'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none');

wanted_size=[50 30];
set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])

folder_to_save=[keys.basepath_to_save keys.project_version];
if ~exist([folder_to_save filesep keys.subfolder_to_save],'dir');
    mkdir(folder_to_save,keys.subfolder_to_save);
end

%% just to have different pdfs
keys.all_titles=keys.CC.epochs;
if strcmp(keys.CC.plot_type,'per_epoch_2levels')
    keys.CC.plot_type=[keys.CC.first_level_factor '_' keys.CC.first_level_epochs{1} ' and '];
end
condition_title=ph_get_condition_title(keys);

export_fig(gcf, [folder_to_save filesep keys.subfolder_to_save filesep keys.monkey ' ' [keys.CC.tasktypes{:}] ' '  keys.arrangement(1:3) ' ' keys.CC.IC_to_plot ...
    ' ' condition_title ' ' keys.selection_title{:} ' ' [keys.CC.plot_type]  ' ' keys.CC.factor ', ' keys.title_part], '-pdf','-transparent') % pdf by run
close(gcf);
