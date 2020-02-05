function pie_data=ph_anova_cell_count(modified_keys)


for fn=fieldnames(modified_keys)'
    keys.(fn{:})=modified_keys.(fn{:});
end
[tuning_per_unit_table]=ph_load_extended_tuning_table(keys);
[tuning_per_unit_table, keys.selection_title]=ph_reduce_tuning_table(tuning_per_unit_table,keys);
tuning_per_unit_table(cellfun(@(x) isempty(x) & ~islogical(x),tuning_per_unit_table))={''};

%% sort of obsolete
keys.tuning_table=tuning_per_unit_table;

%% ANOVA criterions readout
criterions={'space_or_interaction','epoch_or_interaction','hands_or_interaction','SXH_or_interaction'};
for c=1:numel(criterions)
    parameter_criterion_columns=~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),criterions{c}));
    for t=1:numel(keys.CC.tasktypes)
        task_criterion_columns(t,:)=~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),keys.CC.tasktypes{t}));
    end
    criterion_columns=parameter_criterion_columns & any(task_criterion_columns,1);
    keys.(criterions{c})=any(cell2mat(tuning_per_unit_table(2:end,criterion_columns)),2);
end


%% color and legend assignment
epochlegend={'en','bi','su','No tuning','No anova'};
epochalllegend={'en','su','No tuning','No anova'};
choicelegend={'IN','CH','No tuning','No anova'};
keys.SxH_mod_index=0;
S_xH_legend={'IS','IHIS','IH','IHCS','CS','CHCS','CH','CHIS'};
spacelegend={'CS','IS','No tuning','No anova'};
handlegend={'CH','IH','No tuning','No anova'};

inactivation_legend={'CS:EN','CS:EN & IS:EN','IS:EN','CS:SU & IS:EN','CS:SU','CS:SU & IS:SU','IS:SU','CS:EN & IS:SU', 'No tuning'};
space_position_legend={'CS','IS','Position','No tuning'};
crossed_uncrossed_legend={'UC','CR','No tuning','No anova'};
visuomotor_legend={'Vis','Vmot','Mot','No anova'};
fixation_x_position_comb_legend={'Position','Gaze','G+P','Any GxP','None'};
fixation_x_position_CI_legend={'Position','Gaze','G+P','None'};
fixation_x_position_legend={'Position','Gaze','G+P','G+P+GxP','P+GxP','G+GxP','GxP','None'};
eccentricity_x_angle_legend={'All','Ecc + Ang','Ecc + ExA','Ecc','Ang + ExA','Ang','ExA','None'};
gaze_legend={'Tar','Tar+Gaze','Gaze Fhol','None'};
gaze_and_fixation_x_position_legend={'Position','Gaze','G+P','G+P+GxP','P+GxP','G+GxP','GxP','None'};

legend1={};
legend2={};
cols=keys.colors;
if strcmp(keys.CC.factor,'hand')
    keys.all_colors{1}   =[cols.CH_IN; cols.IH_IN; cols.NO_TU; cols.NO_AN]/255;
    legend1=handlegend;
elseif strcmp(keys.CC.factor,'space')
    keys.all_colors{1}   =[cols.CS_IN; cols.IS_IN; cols.NO_TU; cols.NO_AN]/255;
    legend1=spacelegend;
elseif strcmp(keys.CC.factor,'space_position')
    keys.all_colors{1}   =[cols.CS_IN; cols.IS_IN; cols.EP_BI; cols.NO_TU]/255;
    legend1=space_position_legend;
elseif strcmp(keys.CC.factor,'position_space')
    keys.all_colors{1}   =[cols.CS_IN; cols.IS_IN; cols.EP_BI; cols.NO_TU]/255;
    legend1=space_position_legend;
elseif strcmp(keys.CC.factor,'choice1') || strcmp(keys.CC.factor,'choice2')
    keys.all_colors{1}   =[cols.IS_CH; cols.IS_IN;  cols.NO_TU; cols.NO_AN;]/255;
    legend1=choicelegend;
elseif strcmp(keys.CC.factor,'epoch') %% epoch colors...
    keys.all_colors{1}   =[cols.EP_EN; cols.EP_BI; cols.EP_SU; cols.NO_TU; cols.NO_AN]/255;
    legend1=epochlegend;
elseif strcmp(keys.CC.factor,'epoch_all') %% epoch colors...
    keys.all_colors{1}   =[cols.EP_EN; cols.EP_SU; cols.NO_TU; cols.NO_AN]/255;
    legend1=epochalllegend;
end

all_titles=keys.CC.epochs;
if strcmp(keys.CC.plot_type,'per_epoch')
    x_labels=keys.CC.tasktypes;
elseif strcmp(keys.CC.plot_type,'per_task')
    all_titles=keys.CC.tasktypes;
    x_labels=keys.CC.epochs;
elseif strcmp(keys.CC.plot_type,'hemi_and_epoch')
    x_labels={'space','epoch'};
    keys.all_colors{1}   =[cols.EP_EN; cols.EP_BI; cols.EP_SU; cols.NO_TU; cols.NO_AN]/255;
    keys.all_colors{2}   =[cols.CS_IN; cols.IS_IN; cols.NO_TU; cols.NO_AN]/255;
    legend1=epochlegend;
    legend2=spacelegend;

elseif strcmp(keys.CC.plot_type,'space_and_epoch')
    x_labels={'space','epoch'};
    keys.all_colors{1}   =[cols.EP_EN; cols.EP_BI; cols.EP_SU; cols.NO_TU; cols.NO_AN]/255;
    keys.all_colors{2}   =[cols.CS_IN; cols.IS_IN; cols.EP_BI; cols.NO_TU]/255;
    legend1=epochlegend;
    legend2=space_position_legend;
elseif strcmp(keys.CC.plot_type,'space_and_bilateral') %???
    x_labels={'space','bilateral'};
    keys.all_colors{1}   =[cols.CS_IN; cols.IS_IN; cols.NO_TU; cols.NO_AN]/255;
    keys.all_colors{2}   =[cols.EP_EN; cols.EP_SU; cols.NO_AN; cols.NO_TU]/255;
elseif strcmp(keys.CC.plot_type,'space_and_hand')
    x_labels={'space','hand'};
    keys.all_colors{1}   =[cols.CS_IN; cols.IS_IN; cols.NO_TU; cols.NO_AN]/255;
    keys.all_colors{2}   =[cols.CH_IN; cols.IH_IN; cols.NO_TU; cols.NO_AN]/255;
    legend1=spacelegend;
    legend2=handlegend;
elseif strcmp(keys.CC.plot_type,'hand_and_epoch') % reverse order
    x_labels={'epoch','hand'};
    keys.all_colors{1}   =[cols.EP_EN; cols.EP_BI; cols.EP_SU; cols.NO_TU; cols.NO_AN]/255;
    keys.all_colors{2}   =[cols.CH_IN; cols.IH_IN; cols.NO_TU; cols.NO_AN]/255;
    legend1=epochlegend;
    legend2=handlegend;
elseif strcmp(keys.CC.plot_type,'space_x_hand')
    x_labels={'Main effects','Interactions'};
    keys.all_colors{1}   =[cols.IS_IN; cols.IH_IS_IN; cols.IH_IN; cols.IH_CS_IN;...
        cols.CS_IN; cols.CH_CS_IN; cols.CH_IN; cols.CH_IS_IN; cols.NO_TU; cols.NO_AN]/255;
    keys.all_colors{2}   =[255 0 0; 125 0 0; cols.NO_TU; cols.NO_AN]/255;
    legend1=S_xH_legend;
    legend2=crossed_uncrossed_legend;
elseif strcmp(keys.CC.plot_type,'hands_inactivation')
    x_labels={'Contra hand','Ipsi hand'};
%     keys.all_colors{1}   =[cols.CS_IN; (cols.CS_IN+cols.IS_IN)/2; cols.IS_IN; (255-cols.CS_IN+cols.IS_IN)/2;...
%         255-cols.CS_IN; (510-cols.CS_IN-cols.IS_IN)/2; 255-cols.IS_IN;(255+cols.CS_IN-cols.IS_IN)/2; cols.NO_TU]/255;
%     keys.all_colors{2}   =[cols.CS_IN; (cols.CS_IN+cols.IS_IN)/2; cols.IS_IN; (255-cols.CS_IN+cols.IS_IN)/2;...
%         255-cols.CS_IN; (510-cols.CS_IN-cols.IS_IN)/2; 255-cols.IS_IN;(255+cols.CS_IN-cols.IS_IN)/2; cols.NO_TU]/255;
    keys.all_colors{1}   =[255 0 0; 255 77 10; 255 153 20; 255 255 0;...
        0 255 178; 0 178 166; 0 100 255;128 50 205; cols.NO_TU]/255;
    keys.all_colors{2}   =[255 0 0; 255 77 10; 255 153 20; 255 255 0;...
        0 255 178; 0 178 166; 0 100 255;128 50 205; cols.NO_TU]/255;
    legend1=inactivation_legend;
    legend2=inactivation_legend;
elseif strcmp(keys.CC.plot_type,'fixation_x_position_comb') %% something strange
    x_labels={''}; %%?
    keys.all_colors{1}   =[247 148 30; 238 43 123; 159 31 99; 199 45 196; 255 255 255]/255;
    legend1=fixation_x_position_comb_legend;
elseif strcmp(keys.CC.plot_type,'fixation_x_position_CI') %% something strange
    x_labels={''}; %%?
    keys.all_colors{1}   =[247 148 30; 238 43 123; 159 31 99; 255 255 255]/255;
    legend1=fixation_x_position_CI_legend;
elseif strcmp(keys.CC.plot_type,'fixation_x_position') %% something strange
    x_labels={'Retinotopic','initial gaze'}; %%?
    keys.all_colors{1}   =[222 220 0; 255 150 0; 255 0 0; 255 0 150; 171 0 252; 0 0 255; 125 130 255; 255 255 255]/255;
    keys.all_colors{2}   =[0 255 0; 255 255 255]/255;
    legend1=fixation_x_position_legend;
    legend2={'None','Gaze Fhol'}; %% CHECK 
elseif strcmp(keys.CC.plot_type,'eccentricity_x_angle') %% something strange
    x_labels={''}; %%?
    keys.all_colors{1}   =[222 220 0; 255 150 0; 255 0 0; 255 0 150; 171 0 252; 0 0 255; 125 130 255; 255 255 255]/255;
    legend1=eccentricity_x_angle_legend;
elseif strcmp(keys.CC.plot_type,'gaze_and_fixation_x_position') 
    x_labels={'Retinotopic','Object centered'};
    keys.all_colors{1}   =[22 220 0; 255 150 0; 255 0 0; 255 0 150; 171 0 252; 0 0 255; 125 130 255; 255 255 255]/255;
    keys.all_colors{2}   =[127 127 0; 0 255 0;0 127 127;255 255 255]/255;
    legend1=gaze_and_fixation_x_position_legend;
    legend2=gaze_legend;
elseif strcmp(keys.CC.plot_type,'gaze')
    x_labels={'Retinotopic','Object centered'}; 
    keys.all_colors{1}   =[222 220 0; 255 150 0; 255 0 0; 255 255 255]/255;
    keys.all_colors{2}   =[0 255 0; 255 255 255]/255;
    legend1=gaze_legend;
    legend2={'Interaction','None'}; %% CHECK 
elseif strcmp(keys.CC.plot_type,'visuomotor')
    x_labels={''}; %%?
    all_titles={'category'};
    keys.all_colors{1}   =[cols.CS_IN; cols.EP_BI; cols.IS_IN; cols.NO_AN]/255;
    legend1=visuomotor_legend;
end

keys.n1=numel(legend1);
keys.n2=numel(legend2);

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

figure('units','normalized','outerposition',[0 0 1 1])
set(0, 'DefaultFigureRenderer', 'painters');
hold on
N_columns_rows=ceil(sqrt(size(pie_data,1)));
r_dec=1/size(pie_data,2);
for sub= 1:size(pie_data,1)
    subplot(N_columns_rows,N_columns_rows, sub);
    title(all_titles{sub})
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
        if strcmp(keys.CC.plot_type,'per_epoch') || strcmp(keys.CC.plot_type,'per_task') || ...
                strcmp(keys.CC.plot_type,'visuomotor')
            current_colors=repmat(keys.all_colors{1},numel(x)/(size(keys.all_colors{1},1)),1);
        elseif strcmp(keys.CC.plot_type,'space_and_epoch') || strcmp(keys.CC.plot_type,'hemi_and_epoch') || strcmp(keys.CC.plot_type,'space_x_hand') ||...
                strcmp(keys.CC.plot_type,'space_and_hand') || strcmp(keys.CC.plot_type,'hand_and_epoch') ||...
                strcmp(keys.CC.plot_type,'fixation_x_position') || strcmp(keys.CC.plot_type,'gaze_and_fixation_x_position') ||...
                strcmp(keys.CC.plot_type,'fixation_x_position_comb') || strcmp(keys.CC.plot_type,'gaze') ||...
                strcmp(keys.CC.plot_type,'fixation_x_position_CI') ||  strcmp(keys.CC.plot_type,'eccentricity_x_angle')  ||...
                strcmp(keys.CC.plot_type,'space_and_bilateral') || strcmp(keys.CC.plot_type,'hands_inactivation')
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
            piechart(level,sub).handle=pie_chart(x,r,current_colors,[0,0],(1-r_dec/2),x_as_labels);
            text(0,-r+r_dec/2,x_labels{level});
            r=r-r_dec;
        else
            zero_dummie=zeros(size(pie_data,2)-level+1,numel(x));
            x_dummie=repmat(level:(size(pie_data,2)+1),numel(x),1)';
            if strcmp(keys.CC.plot_type,'per_epoch')  %|| strcmp(keys.CC.plot_type,'gaze')
                close(gcf)
                figure('units','normalized','outerposition',[0 0 1 1])
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
            set(gca,'xlim',[0.5,size(pie_data,1)+0.5],'xtick',[1:size(pie_data,1)],'xticklabel',all_titles,'ylim',[0 ylim_bar]);
            break;
        else
            set(gca,'xlim',[0.5,size(pie_data,2)+0.5],'xtick',[1:size(pie_data,2)],'xticklabel',x_labels,'ylim',[0 ylim_bar]);
        end
    end
end

%% legends
if keys.plot.cell_count_legends
    legend(piechart(1,1).handle(1+(0:keys.n1-1)*hstepsize),legend1,'location','best','fontsize',3);
if keys.n2>0
    legend(piechart(2,2).handle(1+(0:keys.n2-1)*hstepsize),legend2,'location','best','fontsize',3);
end
end
keys.title_part='complete';
title_and_save(keys);


%% plot similarity matrix
if ~(strcmp(keys.CC.plot_type,'per_epoch') || strcmp(keys.CC.plot_type,'per_task'))
    return;
end
figure('units','normalized','outerposition',[0 0 1 1])
fields_to_plot= {'same','diff','none'};
for to_plot=1:numel(fields_to_plot)
    FN=fields_to_plot{to_plot};
    for sub= 1:numel(matrix_data)
        subplot(3,numel(matrix_data), sub+(to_plot-1)*numel(matrix_data));
        x_lim=size(matrix_data(sub).(FN),1);
        y_lim=size(matrix_data(sub).(FN),2);
        image(1:x_lim,1:y_lim,matrix_data(sub).(FN)./matrix_data(sub).base*62.5)
        title(all_titles{sub})
        textpositions=combvec(1:x_lim,1:y_lim)';
        if keys.CC.percent
            Matrix_percent=round(matrix_data(sub).(FN)./repmat(diag(matrix_data(sub).same),1,size(matrix_data(sub).(FN),1))*100);
            text(textpositions(:,2),textpositions(:,1),num2str(Matrix_percent(:)),'horizontalalignment','center')
        else
            text(textpositions(:,2),textpositions(:,1),num2str(matrix_data(sub).(FN)(:)),'horizontalalignment','center')
        end
        set(gca,'xlim',[0.5 x_lim+0.5],'xtick',1:x_lim,'xticklabel',x_labels,'ylim',[0.5 y_lim+0.5],'ytick',1:y_lim,'yticklabel',x_labels);
        if sub==1
            ylabel(FN);
        end
        axis square
    end
end
keys.title_part='one_by_one';
title_and_save(keys);

function [multilevel_data matrix_data] = summary_anova_multilevel(keys)

xlsx_table=keys.tuning_table;

for i=1:size(keys.tuning_table,2)
    param(i) = keys.tuning_table(1,i);
    idx.(param{i})= find_column_index(xlsx_table,param{i});
end

multilevel_data={};
matrix_data=[];
no_anova_marker=keys.n1;
no_tuning_marker=keys.n1-1;

switch keys.CC.plot_type
    case 'per_epoch'
        keys.n=keys.n1;
        for e=1:numel(keys.CC.epochs)
            for d=1:numel(keys.CC.tasktypes)
                switch keys.CC.factor
                    case 'epoch'
                        tuning_variables{1}= [keys.CC.IC_to_plot '_AH_' keys.CC.epochs{e} '_epoch_'   keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
                    case 'epoch_all'
                        tuning_variables{1}= [keys.CC.IC_to_plot '_' keys.CC.epochs{e} '_epoch_'   keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
                    case 'space'
                        tuning_variables{1}= [keys.CC.IC_to_plot '_' keys.CC.epochs{e} '_spaceLR_' keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
                    case 'space_position'
                        hand=keys.labels.handsIC{keys.tt.hands(1)+1};
                        tuning_variables{1}= [keys.CC.IC_to_plot '_' keys.CC.epochs{e} '_spaceLR_' keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
                        tuning_variables{2}= [keys.CC.IC_to_plot '_' hand '_' keys.CC.epochs{e} '_position_' keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
                    case 'position_space'
                        hand=keys.labels.handsIC{keys.tt.hands(1)+1};
                        tuning_variables{1}= [keys.CC.IC_to_plot '_' keys.CC.epochs{e} '_spaceLR_' keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
                        tuning_variables{2}= [keys.CC.IC_to_plot '_' hand '_' keys.CC.epochs{e} '_position_' keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
                    case 'hand'
                        tuning_variables{1}= [keys.CC.IC_to_plot '_' keys.CC.epochs{e} '_hands_'   keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
                    case 'choice1'
                        tuning_variables{1}= [keys.CC.IC_to_plot '_AH_' keys.CC.epochs{e} '_RF_choice1_'     keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
                        
                end
                tuning{d}=read_table_column_detail(keys,xlsx_table,idx,tuning_variables);
            end
            pie_tmp=[];
            pie_tmp=get_pie_multilevel(keys,pie_tmp,tuning);
            multilevel_data(e,:)=pie_tmp;
            
            for d1=1:numel(keys.CC.tasktypes)
                basic_quantity=tuning{d1}~=no_tuning_marker & tuning{d1}~=no_anova_marker;
                for d2=1:numel(keys.CC.tasktypes)
                    matrix_data(e).same(d1,d2)=sum(tuning{d1}(basic_quantity)==tuning{d2}(basic_quantity));
                    matrix_data(e).diff(d1,d2)=sum(tuning{d1}(basic_quantity)~=tuning{d2}(basic_quantity) & tuning{d2}(basic_quantity)~=no_tuning_marker);
                    matrix_data(e).none(d1,d2)=sum(tuning{d2}(basic_quantity)==no_tuning_marker);
                    matrix_data(e).base(d1,d2)=sum(basic_quantity);
                end
            end
        end
        
    case 'per_task'
        keys.n=keys.n1;
        for d=1:numel(keys.CC.tasktypes)
            for e=1:numel(keys.CC.epochs)
                switch keys.CC.factor
                    case 'epoch'
                        tuning_variables{1}= [keys.CC.IC_to_plot '_AH_' keys.CC.epochs{e} '_epoch_'   keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
                    case 'epoch_all'
                        tuning_variables{1}= [keys.CC.IC_to_plot '_' keys.CC.epochs{e} '_epoch_'   keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
                    case 'space'
                        tuning_variables{1}= [keys.CC.IC_to_plot '_' keys.CC.epochs{e} '_spaceLR_' keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
                    case 'space_position'
                        hand=keys.labels.handsIC{keys.tt.hands(1)+1};
                        tuning_variables{1}= [keys.CC.IC_to_plot '_' keys.CC.epochs{e} '_spaceLR_' keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
                        tuning_variables{2}= [keys.CC.IC_to_plot '_' hand '_' keys.CC.epochs{e} '_position_' keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
                    case 'position_space'
                        hand=keys.labels.handsIC{keys.tt.hands(1)+1};
                        tuning_variables{1}= [keys.CC.IC_to_plot '_' keys.CC.epochs{e} '_spaceLR_' keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
                        tuning_variables{2}= [keys.CC.IC_to_plot '_' hand '_' keys.CC.epochs{e} '_position_' keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
                    case 'hand'
                        tuning_variables{1}= [keys.CC.IC_to_plot '_' keys.CC.epochs{e} '_hands_'   keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
                    case 'choice1'
                        tuning_variables{1}= [keys.CC.IC_to_plot '_AH_' keys.CC.epochs{e} '_RF_choice1_'     keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
                end
                tuning{e}=read_table_column_detail(keys,xlsx_table,idx,tuning_variables);
            end
            pie_tmp=[];
            pie_tmp=get_pie_multilevel(keys,pie_tmp,tuning);
            multilevel_data(d,:)=pie_tmp;
            
            for e1=1:numel(keys.CC.epochs)
                basic_quantity=tuning{e1}~=no_tuning_marker & tuning{e1}~=no_anova_marker;
                for e2=1:numel(keys.CC.epochs)
                    matrix_data(d).same(e1,e2)=sum(tuning{e1}(basic_quantity)==tuning{e2}(basic_quantity));
                    matrix_data(d).diff(e1,e2)=sum(tuning{e1}(basic_quantity)~=tuning{e2}(basic_quantity) & tuning{e2}(basic_quantity)~=no_tuning_marker);
                    matrix_data(d).none(e1,e2)=sum(tuning{e2}(basic_quantity)==no_tuning_marker);
                    matrix_data(d).base(e1,e2)=sum(basic_quantity);
                end
            end
        end
        
    case 'space_and_epoch'
        for e=1:numel(keys.CC.epochs)
            d=1;
            keys.CC.factor='epoch';
            keys.n=keys.n1;
            tuning_variables{1}= [keys.CC.IC_to_plot '_AH_' keys.CC.epochs{e} '_epoch_'   keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tuninge{d}=read_table_column_detail(keys,xlsx_table,idx,tuning_variables);
            
            keys.CC.factor='space_position';
            keys.n=keys.n2;
            hand=keys.labels.handsIC{keys.tt.hands(1)+1};
            tuning_variables{1}= [keys.CC.IC_to_plot '_' keys.CC.epochs{e} '_spaceLR_' keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tuning_variables{2}= [keys.CC.IC_to_plot '_' hand '_' keys.CC.epochs{e} '_position_' keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tunings{d}=read_table_column_detail(keys,xlsx_table,idx,tuning_variables);
            
            pie_tmp=[];
            pie_tmp=get_pie_two_levels(keys,pie_tmp,tuninge,tunings);
            multilevel_data(e,:)=pie_tmp;
        end

    case 'hemi_and_epoch'
        for e=1:numel(keys.CC.epochs)
            d=1;
            keys.CC.factor='epoch';
            keys.n=keys.n1;
            tuning_variables{1}= [keys.CC.IC_to_plot '_AH_' keys.CC.epochs{e} '_epoch_'   keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tuninge{d}=read_table_column_detail(keys,xlsx_table,idx,tuning_variables);
            
            keys.CC.factor='space';
            keys.n=keys.n2;
            tuning_variables{1}= [keys.CC.IC_to_plot '_' keys.CC.epochs{e} '_spaceLR_' keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tunings{d}=read_table_column_detail(keys,xlsx_table,idx,tuning_variables);
            
            pie_tmp=[];
            pie_tmp=get_pie_two_levels(keys,pie_tmp,tuninge,tunings);
            multilevel_data(e,:)=pie_tmp;
        end
        
    case 'space_and_bilateral' %% does that work??
        hand=keys.labels.handsIC{keys.tt.hands(1)+1};
        for e=1:numel(keys.CC.epochs)
            d=1;
            keys.CC.factor='space';
            keys.n=keys.n1;
            tuning_variables{1}= [keys.CC.IC_to_plot '_' hand '_' keys.CC.epochs{e} '_position_' keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tunings{d}=read_table_column_detail(keys,xlsx_table,idx,tuning_variables);
            
            keys.CC.factor='epoch';
            keys.n=keys.n2;
            tuning_variables{1}= [keys.CC.IC_to_plot '_' hand '_' keys.CC.epochs{e} '_epoch_bilateral_'   keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tuningb{d}=read_table_column_detail(keys,xlsx_table,idx,tuning_variables);
            
            pie_tmp=[];
            pie_tmp=get_pie_two_levels(keys,pie_tmp,tunings,tuningb);
            multilevel_data(e,:)=pie_tmp;
        end
        
    case 'hand_and_epoch'
        for e=1:numel(keys.CC.epochs)
            d=1;
            
            keys.CC.factor='epoch';
            keys.n=keys.n1;
            tuning_variables{1}= [keys.CC.IC_to_plot '_AH_' keys.CC.epochs{e} '_epoch_'   keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tuninge{d}=read_table_column_detail(keys,xlsx_table,idx,tuning_variables);            
            
            keys.CC.factor='hand';
            keys.n=keys.n2;
            tuning_variables{1}= [keys.CC.IC_to_plot '_' keys.CC.epochs{e} '_hands_'   keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tuningh{d}=read_table_column_detail(keys,xlsx_table,idx,tuning_variables);
            
            pie_tmp=[];
            pie_tmp=get_pie_two_levels(keys,pie_tmp,tuninge,tuningh);
            multilevel_data(e,:)=pie_tmp;
        end
        
    case 'space_and_hand'
        for e=1:numel(keys.CC.epochs)
            d=1;
            keys.CC.factor='space';
            keys.n=keys.n1;
            tuning_variables{1}= [keys.CC.IC_to_plot '_' keys.CC.epochs{e} '_spaceLR_' keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tunings{d}=read_table_column_detail(keys,xlsx_table,idx,tuning_variables);
            
            keys.CC.factor='hand';
            keys.n=keys.n2;
            tuning_variables{1}= [keys.CC.IC_to_plot '_' keys.CC.epochs{e} '_hands_'   keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tuningh{d}=read_table_column_detail(keys,xlsx_table,idx,tuning_variables);
            
            pie_tmp=[];
            pie_tmp=get_pie_two_levels(keys,pie_tmp,tunings,tuningh);
            multilevel_data(e,:)=pie_tmp;
        end
        
    case 'space_x_hand'
        for e=1:numel(keys.CC.epochs)
            d=1;
            keys.CC.factor='space_and_hand';
            keys.n=keys.n1;
            tuning_variables{1}= [keys.CC.IC_to_plot '_' keys.CC.epochs{e} '_spaceLR_' keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tuning_variables{2}= [keys.CC.IC_to_plot '_' keys.CC.epochs{e} '_hands_'   keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tuningh{d}=read_table_column_detail(keys,xlsx_table,idx,tuning_variables);
            
            keys.CC.factor='space_x_hand';
            keys.n=keys.n2;
            clear tuning_variables
            tuning_variables{1}= [keys.CC.IC_to_plot '_' keys.CC.epochs{e} '_SxH_'     keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tuningx{d}=read_table_column_detail(keys,xlsx_table,idx,tuning_variables);
            pie_tmp=[];
            pie_tmp=get_pie_two_levels(keys,pie_tmp,tuningh,tuningx);
            multilevel_data(e,:)=pie_tmp;
        end
        
    case 'hands_inactivation'
        for e=1:numel(keys.CC.epochs)
            d=1;
            keys.CC.factor='inactivation_per_hand';
            keys.n=keys.n1;
            tuning_variables{1}= [keys.CC.IC_to_plot '_CH_CS_' keys.CC.epochs{e} '_PT_' keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tuning_variables{2}= [keys.CC.IC_to_plot '_CH_IS_' keys.CC.epochs{e} '_PT_' keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tuningh{d}=read_table_column_detail(keys,xlsx_table,idx,tuning_variables);
            
            keys.CC.factor='inactivation_per_hand';
            keys.n=keys.n2;
            tuning_variables{1}= [keys.CC.IC_to_plot '_IH_CS_' keys.CC.epochs{e} '_PT_' keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tuning_variables{2}= [keys.CC.IC_to_plot '_IH_IS_' keys.CC.epochs{e} '_PT_' keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tuninge{d}=read_table_column_detail(keys,xlsx_table,idx,tuning_variables);
            pie_tmp=[];
            pie_tmp=get_pie_two_levels(keys,pie_tmp,tuningh,tuninge);
            multilevel_data(e,:)=pie_tmp;
        end
        
    case 'fixation_x_position_comb'
        for e=1:numel(keys.CC.epochs)
            d=1;
            keys.CC.factor='fixation_x_position_comb';
            keys.n=keys.n1;
            tuning_variables{1}= [keys.CC.IC_to_plot '_AH_' keys.CC.epochs{e} '_fixation_'   keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tuning_variables{2}= [keys.CC.IC_to_plot '_AH_' keys.CC.epochs{e} '_position_'   keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tuning_variables{3}= [keys.CC.IC_to_plot '_AH_' keys.CC.epochs{e} '_PxF_'        keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tuning{d}=read_table_column_detail(keys,xlsx_table,idx,tuning_variables);
            pie_tmp=[];
            pie_tmp=get_pie_multilevel(keys,pie_tmp,tuning);
            multilevel_data(e,:)=pie_tmp;
        end
     
    case 'fixation_x_position_CI'
        for e=1:numel(keys.CC.epochs)
            d=1;
            keys.CC.factor='fixation_x_position_CI';
            keys.n=keys.n1;
            tuning_variables{1}= [keys.CC.IC_to_plot '_' keys.CC.epochs{e} '_spaceLR_'   keys.CC.tasktypes{d} '_' 'fix' ]; % hardcoded here
            tuning_variables{2}= [keys.CC.IC_to_plot '_' keys.CC.epochs{e} '_spaceLR_'   keys.CC.tasktypes{d} '_' 'mov' ]; % hardcoded here
            tuning{d}=read_table_column_detail(keys,xlsx_table,idx,tuning_variables);
            pie_tmp=[];
            pie_tmp=get_pie_multilevel(keys,pie_tmp,tuning);
            multilevel_data(e,:)=pie_tmp;
        end
    
     
    case 'eccentricity_x_angle'
        for e=1:numel(keys.CC.epochs)
            d=1;
            keys.CC.factor='eccentricity_x_angle';
            keys.n=keys.n1;
            tuning_variables{1}= [keys.CC.IC_to_plot '_AH_' keys.CC.epochs{e} '_distance_'   keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ]; 
            tuning_variables{2}= [keys.CC.IC_to_plot '_AH_' keys.CC.epochs{e} '_angle_'   keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ]; 
            tuning_variables{3}= [keys.CC.IC_to_plot '_AH_' keys.CC.epochs{e} '_DxA_'   keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ]; 
            tuning{d}=read_table_column_detail(keys,xlsx_table,idx,tuning_variables);
            pie_tmp=[];
            pie_tmp=get_pie_multilevel(keys,pie_tmp,tuning);
            multilevel_data(e,:)=pie_tmp;
        end
    
    case 'gaze_and_fixation_x_position'
        for e=1:numel(keys.CC.epochs)
            %for d=1:numel(keys.CC.tasktypes)
            d=1;
            keys.CC.factor='fixation_x_position';
            keys.n=keys.n1;
            tuning_variables{1}= [keys.CC.IC_to_plot '_AH_' keys.CC.epochs{e} '_fixation_'   keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tuning_variables{2}= [keys.CC.IC_to_plot '_AH_' keys.CC.epochs{e} '_position_'   keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tuning_variables{3}= [keys.CC.IC_to_plot '_AH_' keys.CC.epochs{e} '_PxF_'        keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tuningm{d}=read_table_column_detail(keys,xlsx_table,idx,tuning_variables);
            
            keys.CC.factor='fixation_x_position';
            keys.n=keys.n2;
            tuning_variables{1}= [keys.CC.IC_to_plot '_AH_' 'Fhol' '_position_'   keys.CC.tasktypes{d} '_' 'fix' ];            % hardcoded here
            tuning_variables{2}= [keys.CC.IC_to_plot '_AH_' keys.CC.epochs{e} '_position_'   keys.CC.tasktypes{d} '_' 'tar' ]; % hardcoded here
            tuning_variables{3}= [keys.CC.IC_to_plot '_AH_' keys.CC.epochs{e} '_PxF_'        keys.CC.tasktypes{d} '_' 'tar' ]; % hardcoded here
            tuningt{d}=read_table_column_detail(keys,xlsx_table,idx,tuning_variables);
            %end
            pie_tmp=[];
            pie_tmp=get_pie_two_levels(keys,pie_tmp,tuningm,tuningt);
            multilevel_data(e,:)=pie_tmp;
        end
        
    case 'gaze'
            keys.n=keys.n1;
        for e=1:numel(keys.CC.epochs)
            %for d=1:numel(keys.CC.tasktypes)
                
            d=1;
            keys.CC.factor='gaze';
            keys.n=keys.n1;
            tuning_variables{1}= [keys.CC.IC_to_plot '_AH_' 'Fhol' '_position_'              keys.CC.tasktypes{d} '_' 'fix' ];            % hardcoded here
            tuning_variables{2}= [keys.CC.IC_to_plot '_AH_' keys.CC.epochs{e} '_position_'   keys.CC.tasktypes{d} '_' 'tar' ]; % hardcoded here
            tuningm{d}=read_table_column_detail(keys,xlsx_table,idx,tuning_variables);
            
            clear tuning_variables
            keys.CC.factor='gaze_interaction';
            keys.n=keys.n2; 
            tuning_variables{1}= [keys.CC.IC_to_plot '_AH_' keys.CC.epochs{e} '_PxF_'   keys.CC.tasktypes{d} '_' 'tar' ]; % hardcoded here
            tuningi{d}=read_table_column_detail(keys,xlsx_table,idx,tuning_variables);
            %end
            %end
            
            pie_tmp=[];
            pie_tmp=get_pie_two_levels(keys,pie_tmp,tuningm,tuningi);
            multilevel_data(e,:)=pie_tmp;
%             pie_tmp=[];
%             pie_tmp=get_pie_multilevel(keys,pie_tmp,tuning);
%             multilevel_data(e,:)=pie_tmp;
            
            %% ?? not sure this is correct here
            for d1=1:numel(tuningm)
                basic_quantity=tuningm{d1}~=no_tuning_marker & tuningm{d1}~=no_anova_marker;
                for d2=1:numel(tuningi)
                    matrix_data(e).same(d1,d2)=sum(tuningm{d1}(basic_quantity)==tuningm{d2}(basic_quantity));
                    matrix_data(e).diff(d1,d2)=sum(tuningm{d1}(basic_quantity)~=tuningi{d2}(basic_quantity) & tuningi{d2}(basic_quantity)~=no_tuning_marker);
                    matrix_data(e).none(d1,d2)=sum(tuningi{d2}(basic_quantity)==no_tuning_marker);
                    matrix_data(e).base(d1,d2)=sum(basic_quantity);
                end
            end
        end
        
    case 'visuomotor'
            keys.n=keys.n1;
        keys.CC.factor='visuomotor'; % overwrite...
        for d=1:numel(keys.CC.tasktypes)
            tuning_variables{1}= ['visual_'   keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tuning_variables{2}= ['motor_'   keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tuning_variables{3}= ['visuomotor_'   keys.CC.tasktypes{d} '_' keys.arrangement(1:3) ];
            tuning{1}=read_table_column_detail(keys,xlsx_table,idx,tuning_variables);
            
            pie_tmp=[];
            pie_tmp=get_pie_multilevel(keys,pie_tmp,tuning);
            multilevel_data(d,:)=pie_tmp;
            
            %doesnt make much sense for this plot type, we define it nevertheless not to have errors
            matrix_data(d).same=1;
            matrix_data(d).diff=1;
            matrix_data(d).none=1;
            matrix_data(d).base=1;
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
    case 'space_and_hand'
        IStuning=ismember(table(2:end,index.(tuning_variables{1})),'IS') & keys.SXH_or_interaction; %keys.space_or_interaction;
        CStuning=ismember(table(2:end,index.(tuning_variables{1})),'CS') & keys.SXH_or_interaction; %keys.space_or_interaction;
        IHtuning=ismember(table(2:end,index.(tuning_variables{2})),'IH') & keys.SXH_or_interaction; %keys.hands_or_interaction;
        CHtuning=ismember(table(2:end,index.(tuning_variables{2})),'CH') & keys.SXH_or_interaction; %keys.hands_or_interaction;
        IHIS=(IHtuning & IStuning);
        IHCS=(IHtuning & CStuning);
        CHIS=(CHtuning & IStuning);
        CHCS=(CHtuning & CStuning);
        IH=IHtuning & ~IStuning & ~CStuning;
        CH=CHtuning & ~IStuning & ~CStuning;
        IS=IStuning & ~IHtuning & ~CHtuning;
        CS=CStuning & ~IHtuning & ~CHtuning;
        %na=~keys.space_or_interaction & ~keys.hands_or_interaction;
        na=~keys.SXH_or_interaction;
        read_out(IS)            =mod(0+keys.SxH_mod_index,8)+1;
        read_out(IHIS)          =mod(1+keys.SxH_mod_index,8)+1;
        read_out(IH)            =mod(2+keys.SxH_mod_index,8)+1;
        read_out(IHCS)          =mod(3+keys.SxH_mod_index,8)+1;
        read_out(CS)            =mod(4+keys.SxH_mod_index,8)+1;
        read_out(CHCS)          =mod(5+keys.SxH_mod_index,8)+1;
        read_out(CH)            =mod(6+keys.SxH_mod_index,8)+1;
        read_out(CHIS)          =mod(7+keys.SxH_mod_index,8)+1;
        read_out(~IHtuning & ~CHtuning & ~IStuning & ~CStuning & ~na)=9;
        read_out(na)            =10;
        
    case 'hand'
        IH=ismember(table(2:end,index.(tuning_variables{1})),'IH') & keys.hands_or_interaction;
        CH=ismember(table(2:end,index.(tuning_variables{1})),'CH') & keys.hands_or_interaction;
        na=~keys.hands_or_interaction;
        clear read_out
        read_out(CH)                =1;
        read_out(IH)                =2;
        read_out(~IH & ~CH & ~na)   =3;
        read_out(na)                =4;
        
    case 'space'
        IS=ismember(table(2:end,index.(tuning_variables{1})),'IS') & keys.space_or_interaction;
        CS=ismember(table(2:end,index.(tuning_variables{1})),'CS') & keys.space_or_interaction;
        na=~keys.space_or_interaction;
        clear read_out
        read_out(CS)                =1;
        read_out(IS)                =2;
        read_out(~IS & ~CS & ~na)   =3;
        read_out(na)                =4;
        
    case 'epoch'
        en=ismember(table(2:end,index.(tuning_variables{1})),'en') & keys.epoch_or_interaction;
        su=ismember(table(2:end,index.(tuning_variables{1})),'su') & keys.epoch_or_interaction;
        bi=ismember(table(2:end,index.(tuning_variables{1})),'bi') & keys.epoch_or_interaction;
        na=~keys.epoch_or_interaction;
        clear read_out
        read_out(en)                =1;
        read_out(bi)                =2;
        read_out(su)                =3;
        read_out(na)                =5;
        read_out(~en & ~su & ~bi & ~na)   =4;
        
    case 'epoch_all'
        en=ismember(table(2:end,index.(tuning_variables{1})),'en') & keys.epoch_or_interaction;
        su=ismember(table(2:end,index.(tuning_variables{1})),'su') & keys.epoch_or_interaction;
        na=~keys.epoch_or_interaction;
        clear read_out
        read_out(en)                =1;
        read_out(su)                =2;
        read_out(na)                =4;
        read_out(~en & ~su & ~na)   =3;
        
    case 'choice1'
        IN=ismember(table(2:end,index.(tuning_variables{1})),'IN') & keys.space_or_interaction;
        CH=ismember(table(2:end,index.(tuning_variables{1})),'CH') & keys.space_or_interaction;
        na=~keys.space_or_interaction;
        clear read_out
        read_out(IN)                =1;
        read_out(CH)                =2;
        read_out(~IN & ~CH & ~na)   =3;
        read_out(na)                =4;
        
    case 'space_x_hand'
        CR=ismember(table(2:end,index.(tuning_variables{1})),'CR') & keys.SXH_or_interaction;
        UC=ismember(table(2:end,index.(tuning_variables{1})),'UC') & keys.SXH_or_interaction;
        na=~keys.SXH_or_interaction;
        clear read_out
        read_out(UC)                =1;
        read_out(CR)                =2;
        read_out(~CR & ~UC & ~na)   =3;
        read_out(na)                =4;
        
    case 'fixation_x_position'
        %indexfnS=fixaion, indexfnH=position, indexfnX=pxF
        fixtuning=ismember(table(2:end,index.(tuning_variables{1})),'true'); %& keys.hands_or_interaction; %% anova criterion for this one doesnt make much sense, does it??
        postuning=ismember(table(2:end,index.(tuning_variables{2})),'true'); %& keys.hands_or_interaction;
        fxptuning=ismember(table(2:end,index.(tuning_variables{3})),'true'); %& keys.space_or_interaction;
        %CStuning=ismember(table(2:end,index.(indexfnS)),'CS') & keys.space_or_interaction;
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


%     case 'fixation_x_position'
%         %indexfnS=fixaion, indexfnH=position, indexfnX=pxF
%         fixtuning=ismember(table(2:end,index.(tuning_variables{1})),'true'); %& keys.hands_or_interaction; %% anova criterion for this one doesnt make much sense, does it??
%         postuning=ismember(table(2:end,index.(tuning_variables{2})),'true'); %& keys.hands_or_interaction;
%         fxptuning=ismember(table(2:end,index.(tuning_variables{3})),'true'); %& keys.space_or_interaction;
%         %CStuning=ismember(table(2:end,index.(indexfnS)),'CS') & keys.space_or_interaction;
%         Pos= ~fixtuning  & postuning  & ~fxptuning;
%         Fix=  fixtuning  & ~postuning & ~fxptuning;
%         FnP=  fixtuning  &  postuning & ~fxptuning;
%         All=  fixtuning  &  postuning &  fxptuning;
%         PnX= ~fixtuning  &  postuning &  fxptuning;
%         FnX=  fixtuning  & ~postuning &  fxptuning;
%         FxP= ~fixtuning  & ~postuning &  fxptuning;
%         Na = ~fixtuning  & ~postuning & ~fxptuning;
%         
%         read_out(Pos)          =1;
%         read_out(Fix)          =2;
%         read_out(FnP)          =3;
%         read_out(All)          =4;
%         read_out(PnX)          =5;
%         read_out(FnX)          =6;
%         read_out(FxP)          =7;
%         read_out(Na)           =8;
        
    case 'eccentricity_x_angle'
        %indexfnS=fixaion, indexfnH=position, indexfnX=pxF
        ecctuning=ismember(table(2:end,index.(tuning_variables{1})),'true'); %& keys.hands_or_interaction; %% anova criterion for this one doesnt make much sense, does it??
        angtuning=ismember(table(2:end,index.(tuning_variables{2})),'true'); %& keys.hands_or_interaction;
        exatuning=ismember(table(2:end,index.(tuning_variables{3})),'true'); %& keys.space_or_interaction;
        %CStuning=ismember(table(2:end,index.(indexfnS)),'CS') & keys.space_or_interaction;
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
        %indexfnS=fixaion, indexfnH=position, indexfnX=pxF
                
        fixtuning=ismember(table(2:end,index.(tuning_variables{1})),'IS') | ismember(table(2:end,index.(tuning_variables{1})),'CS');
        postuning=ismember(table(2:end,index.(tuning_variables{2})),'IS') | ismember(table(2:end,index.(tuning_variables{2})),'CS');
        %CStuning=ismember(table(2:end,index.(indexfnS)),'CS') & keys.space_or_interaction;
        Pos= ~fixtuning  & postuning  ;
        Fix=  fixtuning  & ~postuning ;
        FnP=  fixtuning  &  postuning ;  
        Na = ~fixtuning  & ~postuning ;
        
        read_out(Pos)          =1;
        read_out(Fix)          =2;
        read_out(FnP)          =3;
        read_out(Na)           =4;
        
    case 'gaze'
        %indexfnS=fixaion, indexfnH=position, indexfnX=pxF
        fixtuning=ismember(table(2:end,index.(tuning_variables{1})),'true'); %& keys.hands_or_interaction; %% anova criterion for this one doesnt make much sense, does it??
        tartuning=ismember(table(2:end,index.(tuning_variables{2})),'true'); %& keys.hands_or_interaction;
        
        Tar= ~fixtuning  &  tartuning;
        Fix=  fixtuning  & ~tartuning;
        FnT=  fixtuning  &  tartuning;
        Na = ~fixtuning  & ~tartuning;
        
        read_out(Tar)          =1;
        read_out(FnT)          =2;
        read_out(Fix)          =3;
        read_out(Na)           =4;
    case 'gaze_interaction';
        interaction=ismember(table(2:end,index.(tuning_variables{1})),'true');
        
        read_out(interaction)          =1;
        read_out(~interaction)           =2;
        
    case 'visuomotor'
        vis=cell2mat(table(2:end,index.(tuning_variables{1})));
        mot=cell2mat(table(2:end,index.(tuning_variables{2})));
        vmt=cell2mat(table(2:end,index.(tuning_variables{3})));
        clear read_out
        read_out(vis)                =1;
        read_out(vmt)                =2;
        read_out(mot)                =3;
        read_out(~vis & ~mot & ~vmt)   =4;
        
    case 'inactivation_per_hand'
        CS_EN=ismember(table(2:end,index.(tuning_variables{1})),'EN') & keys.space_or_interaction;
        CS_SU=ismember(table(2:end,index.(tuning_variables{1})),'SU') & keys.space_or_interaction;
        IS_EN=ismember(table(2:end,index.(tuning_variables{2})),'EN') & keys.space_or_interaction;
        IS_SU=ismember(table(2:end,index.(tuning_variables{2})),'SU') & keys.space_or_interaction;
        
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

        CS=ismember(table(2:end,index.(tuning_variables{1})),'CS') ; % & keys.space_or_interaction;
        IS=ismember(table(2:end,index.(tuning_variables{1})),'IS') ; % & keys.space_or_interaction;
        PO=ismember(table(2:end,index.(tuning_variables{2})),'true') ; % & keys.space_or_interaction;
%         na=~keys.space_or_interaction; %% not ideal either!
        clear read_out
        read_out(CS)                =1;
        read_out(IS)                =2;
        read_out(~IS & ~CS & PO)    =3;
        read_out(~IS & ~CS & ~PO)   =4;

    case 'position_space'

        CS=ismember(table(2:end,index.(tuning_variables{1})),'CS') ; % & keys.space_or_interaction;
        IS=ismember(table(2:end,index.(tuning_variables{1})),'IS') ; % & keys.space_or_interaction;
        PO=ismember(table(2:end,index.(tuning_variables{2})),'true') ; % & keys.space_or_interaction;
%         na=~keys.space_or_interaction; %% not ideal either!
        clear read_out
        read_out(CS & PO)                =1;
        read_out(IS & PO)                =2;
        read_out(~IS & ~CS & PO)    =3;
        read_out(~PO)   =4;
        
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

function pie=get_pie_two_levels(keys,pie,tuning1,tuning2)
pie{1}=[];
pie{2}=[];
for n1=1:keys.n1
    pie{1}(end+1)=sum(tuning1{1}==n1);
    for n2=1:keys.n2
        pie{2}(end+1)=sum(tuning1{1}(:)==n1 & tuning2{1}(:)==n2);
    end
end

function title_and_save(keys)
criterions=[' criterions: S: ' keys.tt.space_criterion ' H: ' keys.tt.hands_criterion ' SxH: ' keys.tt.SXH_criterion ' E: ' keys.tt.epoch_criterion];
mtit(gcf,[keys.monkey ' ' [keys.CC.tasktypes{:}] ' ' [keys.CC.plot_type] ' ' keys.CC.factor ' ' keys.CC.IC_to_plot ' ' keys.arrangement(1:3) ' ' keys.selection_title{:} ' N: ' num2str(size(keys.tuning_table,1)-1) ' ' keys.title_part  ' ' criterions], 'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none');

wanted_size=[50 30];
set(gcf, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])

folder_to_save=[keys.drive filesep keys.basepath_to_save filesep keys.project_version];
if ~exist([folder_to_save filesep keys.subfolder_to_save],'dir');
    mkdir(folder_to_save,keys.subfolder_to_save);
end
export_fig(gcf, [folder_to_save filesep keys.subfolder_to_save filesep keys.monkey ' ' [keys.CC.tasktypes{:}] ' '  keys.arrangement(1:3) ' ' keys.CC.IC_to_plot ...
    ' hnd ' num2str(keys.tt.hands) ' ch ' num2str(keys.tt.choices) ' ' keys.selection_title{:} ' ' [keys.CC.plot_type]  ' ' keys.CC.factor ', ' keys.title_part], '-pdf','-transparent') % pdf by run
close(gcf);




