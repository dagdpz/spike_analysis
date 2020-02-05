function pie_data=ph_anova_cell_count(modified_keys)

keys.monkey             = 'Curius';
keys.case               = {'opt'};
keys.target             = 'dPulv_l';
keys.instructed_choice  = 'in';
keys.contra_color       = 'pink';
keys.all_colors         = {[255 0 178; 171 0 252; 0 0 255; 125 130 255; 255 153 20; 222 220 0; 0 255 0; 145 143 56; 255 0 0; 255 255 255]/255}; %% overwritten anyway, never input

keys.cc.tasktypes           = {'Msac'};
keys.cc.percent            = 0;
keys.cc.plot_as_pie        = 0;
keys.cc.epochs             = {''};
keys.cc.space_criterion    = 'none';
keys.cc.epoch_criterion    = 'none';
keys.cc.hands_criterion    = 'none';
keys.cc.SXH_criterion      = 'none';
keys.cc.plot_type          = 'per_epoch';
keys.cc.factors            = 'space';
keys.cc.only_both_hands    = 0;
keys.cc.Selection          = {};

for fn=fieldnames(modified_keys)'
    keys.(fn{:})=modified_keys.(fn{:});
end
[tuning_per_unit_table]=ph_load_extended_tuning_table(keys);
[tuning_per_unit_table, keys.selection_title]=ph_reduce_tuning_table(tuning_per_unit_table,keys);
tuning_per_unit_table(cellfun(@(x) isempty(x) & ~islogical(x),tuning_per_unit_table))={''};

%% sort of obsolete
keys.xlsx_table=tuning_per_unit_table;

%% ANOVA cirterions readout
criterions={'space_or_interaction','epoch_or_interaction','hands_or_interaction','SXH_or_interaction'};
for c=1:numel(criterions)
    parameter_criterion_columns=~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),criterions{c}));
    for t=1:numel(keys.cc.tasktypes)
        task_criterion_columns(t,:)=~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),keys.cc.tasktypes{t}));
    end
    criterion_columns=parameter_criterion_columns & any(task_criterion_columns,1);
    keys.(criterions{c})=any(cell2mat(tuning_per_unit_table(2:end,criterion_columns)),2);
end

%% keys.contra_color defines on which side (of pie plot) to put contralateral (colors depend on side!)
epochlegend={'en','bi','su','No tuning','No anova'};
choicelegend={'IN','CH','No tuning','No anova'};
if strcmp(keys.contra_color,'pink')
    keys.SxH_mod_index=4;%% Ahaa! This is for assigning correct colors in the Hand x Space plot 
    keys.contra_color_index=1;
    keys.ipsi_index=4;
    S_xH_legend={'CS','CHCS','CH','CHIS','IS','IHIS','IH','IHCS'};
    spacelegend={'CS','IS','No tuning','No anova'};
    handlegend={'CH','IH','No tuning','No anova'};
else
    keys.SxH_mod_index=0;
    keys.contra_color_index=4;
    keys.ipsi_index=1;
    S_xH_legend={'IS','IHIS','IH','IHCS','CS','CHCS','CH','CHIS'};
    spacelegend={'IS','CS','No tuning','No anova'};
    handlegend={'IH','CH','No tuning','No anova'};
end

if strcmp(keys.cc.factors,'hand')
    keys.all_colors{1}   =[0 0 255; 255 255 255; 128 128 128; 0 255 0]/255;
elseif strcmp(keys.cc.factors,'space')
    keys.all_colors{1}   =[255 0 178; 255 255 255; 128 128 128; 255 153 20]/255;
elseif strcmp(keys.cc.factors,'choice1') || strcmp(keys.cc.factors,'choice2')
    keys.all_colors{1}   =[255 0 178; 255 255 255; 128 128 128; 255 153 20]/255;
elseif strcmp(keys.cc.factors,'epoch') %% epoch colors...
    keys.all_colors{1}   =[255 0 0; 0 255 0; 0 0 255; 255 255 255; 128 128 128]/255;
end

all_titles=keys.cc.epochs;
if strcmp(keys.cc.plot_type,'per_epoch')
    x_labels=keys.cc.tasktypes;
elseif strcmp(keys.cc.plot_type,'per_task')
    all_titles=keys.cc.tasktypes;
    x_labels=keys.cc.epochs;
elseif strcmp(keys.cc.plot_type,'space_and_epoch')
    x_labels={'space','epoch'};
    keys.all_colors{1}   =[255 0 178; 255 255 255; 128 128 128; 255 153 20]/255;
    keys.all_colors{2}   =[255 0 0; 0 255 0; 0 0 255; 255 255 255; 128 128 128]/255;
elseif strcmp(keys.cc.plot_type,'space_and_hand')
    x_labels={'space','hand'};
    keys.all_colors{1}   =[255 0 178; 255 255 255; 128 128 128; 255 153 20]/255;
    keys.all_colors{2}   =[0 0 255; 255 255 255; 128 128 128; 0 255 0]/255;
elseif strcmp(keys.cc.plot_type,'hand_and_epoch')
    x_labels={'hand','epoch'};
    keys.all_colors{1}   =[0 0 255; 255 255 255; 128 128 128; 0 255 0]/255;
    keys.all_colors{2}   =[255 0 0; 0 255 0; 0 0 255; 255 255 255; 128 128 128]/255;
elseif strcmp(keys.cc.plot_type,'space_x_hand')
    x_labels={'Main effects','Interactions'};
    keys.all_colors{1}   =[255 0 150; 171 0 252; 0 0 255; 125 130 255; 255 150 0; 222 220 0; 0 255 0; 145 143 56; 128 128 128; 255 255 255]/255;
    keys.all_colors{2}   =[255 0 0; 255 255 255; 128 128 128; 125 0 0]/255;
elseif strcmp(keys.cc.plot_type,'fixation_x_position')
    x_labels={'Retinotopic','Object centered'}; %%?
    keys.all_colors{1}   =[255 0 0; 255 0 150; 171 0 252; 0 0 255; 125 130 255; 255 150 0; 222 220 0; 255 255 255]/255;
    keys.all_colors{2}   =[255 0 0; 255 0 150; 171 0 252; 0 0 255; 125 130 255; 255 150 0; 222 220 0; 255 255 255]/255;
end


%% Calculation in subfunctions
[pie_data matrix_data]=summary_anova_multilevel(keys);
n_cells=size(keys.xlsx_table,1)-1;


%% PLOT: pie or bar

if keys.cc.plot_as_pie %% bar or pie
    keys.subfolder_to_save='cell_counts_as_pie';
else
    keys.subfolder_to_save='cell_counts_as_bar';
end
if keys.cc.percent
    keys.subfolder_to_save=[keys.subfolder_to_save '_percent'];
else
    keys.subfolder_to_save=[keys.subfolder_to_save '_absolute'];
end

%if ~strcmp(keys.cc.plot_type,'per_task')
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
    if keys.cc.plot_as_pie %% bar or pie
        levels_in_order=size(pie_data,2):-1:1;
    else
        levels_in_order=1:size(pie_data,2);
    end
    for level=levels_in_order
        if keys.cc.percent
            x=pie_data{sub,level}/n_cells*100;
            ylim_bar=100;
        else
            x=pie_data{sub,level};
            ylim_bar=max([size(keys.xlsx_table,1)-1 1]);
        end
        if strcmp(keys.cc.plot_type,'per_epoch') || strcmp(keys.cc.plot_type,'per_task')
            current_colors=repmat(keys.all_colors{1},numel(x)/(size(keys.all_colors{1},1)),1);
        elseif strcmp(keys.cc.plot_type,'space_and_epoch') || strcmp(keys.cc.plot_type,'space_x_hand') ||...
              strcmp(keys.cc.plot_type,'space_and_hand') || strcmp(keys.cc.plot_type,'hand_and_epoch') ||...
              strcmp(keys.cc.plot_type,'fixation_x_position') 
            current_colors=repmat(keys.all_colors{level},numel(x)/(size(keys.all_colors{level},1)),1);
        end
        x_as_labels=num2cell(x(:));
        if keys.cc.percent
        x_as_labels=cellfun(@(x) [num2str(round(x)) '%'],x_as_labels,'UniformOutput',false);            
        else
        x_as_labels=cellfun(@(x) num2str(round(x)),x_as_labels,'UniformOutput',false);
        end
        if keys.cc.plot_as_pie  %% bar or pie
            hstepsize=2;
            piechart(level,sub).handle=pie_chart(x,r,current_colors,[0,0],(1-r_dec/2),x_as_labels);
            text(0,-r+r_dec/2,x_labels{level});
            r=r-r_dec;
        else
            zero_dummie=zeros(size(pie_data,2)-level+1,numel(x));
            x_dummie=repmat(level:(size(pie_data,2)+1),numel(x),1)';
            piechart(level,sub).handle=bar(x_dummie,[x;zero_dummie],1,'stacked');
            for k=1:numel(piechart(level,sub).handle)
                set(piechart(level,sub).handle(k),'facecolor',current_colors(k,:))
            end
            cumsum_X=cumsum(x);
            cumsum_index=[cumsum_X(1)~=0, diff(cumsum_X)>0];
            text(repmat(level,sum(cumsum_index),1),cumsum_X(cumsum_index)-ylim_bar/25,x_as_labels(cumsum_index));
        end
    end
    if ~keys.cc.plot_as_pie %% bar or pie
        set(gca,'xlim',[0.5,size(pie_data,2)+0.5],'xtick',[1:size(pie_data,2)],'xticklabel',x_labels,'ylim',[0 ylim_bar]);
        ylabel('N units');
        hstepsize=1;
    end
end

%% legends
if strcmp(keys.cc.plot_type,'per_epoch')
    if strcmp(keys.cc.factors,'space')
        legend(piechart(1,1).handle([1 1+3*hstepsize 1+2*hstepsize 1+hstepsize]),spacelegend,'location','best');
    elseif strcmp(keys.cc.factors,'hand')
        legend(piechart(1,1).handle([1 1+3*hstepsize 1+2*hstepsize 1+hstepsize]),handlegend,'location','best');
    elseif strcmp(keys.cc.factors,'choice1') || strcmp(keys.cc.factors,'choice2')
        legend(piechart(1,1).handle([1 1+3*hstepsize 1+2*hstepsize 1+hstepsize]),choicelegend,'location','best');
    elseif strcmp(keys.cc.factors,'epoch')
        legend(piechart(1,1).handle([1 1+hstepsize 1+2*hstepsize 1+4*hstepsize 1+3*hstepsize]),epochlegend,'location','best');
    else
       % legend(piechart(1,1).handle(1:hstepsize:end),{'LS','LHLS','LH','LHRS','RS','RHRS','RH','RHLS','IA','No tuning'},'location','best');
    end
elseif strcmp(keys.cc.plot_type,'per_task')
    if strcmp(keys.cc.factors,'space')
        legend(piechart(1,1).handle([1 1+3*hstepsize 1+2*hstepsize 1+hstepsize]),spacelegend,'location','best');
    elseif strcmp(keys.cc.factors,'hand')
        legend(piechart(1,1).handle([1 1+3*hstepsize 1+2*hstepsize 1+hstepsize]),handlegend,'location','best');
    elseif strcmp(keys.cc.factors,'choice1') || strcmp(keys.cc.factors,'choice2')
        legend(piechart(1,1).handle([1 1+3*hstepsize 1+2*hstepsize 1+hstepsize]),choicelegend,'location','best');
    elseif strcmp(keys.cc.factors,'epoch')
        legend(piechart(1,1).handle([1 1+hstepsize 1+2*hstepsize 1+4*hstepsize 1+3*hstepsize]),epochlegend,'location','best');
    else
       % legend(piechart(1,1).handle(1:hstepsize:end),{'LS','LHLS','LH','LHRS','RS','RHRS','RH','RHLS','IA','No tuning'},'location','best');
    end
elseif strcmp(keys.cc.plot_type,'space_x_hand')
    subplot(N_columns_rows,N_columns_rows, 1)
    %legend(piechart(1,1).handle(1:hstepsize:end),{'LS','LHLS','LH','LHRS','RS','RHRS','RH','RHLS','No tuning','No anova'},'location','best');
    legend(piechart(1,1).handle(1:hstepsize:end),[S_xH_legend {'No tuning','No anova'}],'location','best');
    subplot(N_columns_rows,N_columns_rows, 2)
    legend(piechart(2,2).handle([1 1+3*hstepsize 1+2*hstepsize 1+hstepsize]),{'CR','UC','No tuning','No anova'},'location','best');
elseif strcmp(keys.cc.plot_type,'space_and_epoch')
    subplot(N_columns_rows,N_columns_rows, 1)
    legend(piechart(1,1).handle([1 1+3*hstepsize 1+2*hstepsize 1+hstepsize]),spacelegend,'location','best');
    subplot(N_columns_rows,N_columns_rows, 2)
    legend(piechart(2,2).handle([1 1+hstepsize 1+2*hstepsize 1+4*hstepsize 1+3*hstepsize]),epochlegend,'location','best');
elseif strcmp(keys.cc.plot_type,'space_and_hand')
    subplot(N_columns_rows,N_columns_rows, 1)
    legend(piechart(1,1).handle([1 1+3*hstepsize 1+2*hstepsize 1+hstepsize]),spacelegend,'location','best');
    subplot(N_columns_rows,N_columns_rows, 2)
    legend(piechart(2,2).handle([1 1+4*hstepsize 1+3*hstepsize 1+2*hstepsize 1+hstepsize]),handlegend,'location','best');
elseif strcmp(keys.cc.plot_type,'hand_and_epoch')
    subplot(N_columns_rows,N_columns_rows, 1)
    legend(piechart(1,1).handle([1 1+3*hstepsize 1+2*hstepsize 1+hstepsize]),handlegend,'location','best');
    subplot(N_columns_rows,N_columns_rows, 2)
    legend(piechart(2,2).handle([1 1+hstepsize 1+2*hstepsize 1+4*hstepsize 1+3*hstepsize]),epochlegend,'location','best');
elseif strcmp(keys.cc.plot_type,'fixation_x_position')
    subplot(N_columns_rows,N_columns_rows, 1)
    legend(piechart(1,1).handle([1 1+1*hstepsize 1+2*hstepsize 1+3*hstepsize  1+4*hstepsize  1+5*hstepsize 1+6*hstepsize 1+7*hstepsize]),{'Fix','Mov','F+FxM','M+FxM','F+M','F+M+FxM','FxM','none'},'location','best');
    subplot(N_columns_rows,N_columns_rows, 2)
    legend(piechart(2,2).handle([1 1+1*hstepsize 1+2*hstepsize 1+3*hstepsize  1+4*hstepsize  1+5*hstepsize 1+6*hstepsize 1+7*hstepsize]),{'Fix','Pos','F+FxP','P+FxP','F+P','F+P+FxP','FxP','none'},'location','best');
end
keys.title_part='complete';
title_and_save(keys);
%end

%% plot similarity matrix
if ~(strcmp(keys.cc.plot_type,'per_epoch') || strcmp(keys.cc.plot_type,'per_task'))
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
        if keys.cc.percent
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

xlsx_table=keys.xlsx_table;

for i=1:size(keys.xlsx_table,2)
    param(i) = keys.xlsx_table(1,i);
    idx.(param{i})= find_column_index(xlsx_table,param{i});
end

multilevel_data={};
matrix_data=[];
if strcmp(keys.cc.factors,'hand')    || strcmp(keys.cc.factors,'space') ||...
   strcmp(keys.cc.factors,'choice1') || strcmp(keys.cc.factors,'choice2');
    no_tuning_marker=3;
    no_anova_marker=2;
elseif strcmp(keys.cc.factors,'epoch');
    no_tuning_marker=5;
    no_anova_marker=4;   
else
    no_tuning_marker=9;
    no_anova_marker=10;
end

if strcmp(keys.cc.plot_type,'per_epoch')
    for e=1:numel(keys.cc.epochs)
        for d=1:numel(keys.cc.tasktypes)
            epochtuning{d}= [keys.instructed_choice '_NH_' keys.cc.epochs{e} '_epoch_'   keys.cc.tasktypes{d} '_' keys.case{1} ]; 
            spacetuning{d}= [keys.instructed_choice '_' keys.cc.epochs{e} '_spaceLR_' keys.cc.tasktypes{d} '_' keys.case{1} ];
            handstuning{d}= [keys.instructed_choice '_' keys.cc.epochs{e} '_hands_'   keys.cc.tasktypes{d} '_' keys.case{1} ];
            S_x_Htuning{d}= [keys.instructed_choice '_' keys.cc.epochs{e} '_SxH_'     keys.cc.tasktypes{d} '_' keys.case{1} ];
            choictuning{d}= [keys.instructed_choice '_NH_' keys.cc.epochs{e} '_RF_' keys.cc.factors '_'     keys.cc.tasktypes{d} '_' keys.case{1} ];
            tuning{d}=read_table_column_detail(keys,xlsx_table,idx,epochtuning{d},spacetuning{d},handstuning{d},S_x_Htuning{d},choictuning{d});
        end
        pie_tmp=[];
        pie_tmp=get_pie_multilevel(keys,pie_tmp,tuning);
        multilevel_data(e,:)=pie_tmp;
        
        for d1=1:numel(keys.cc.tasktypes)
            basic_quantity=tuning{d1}~=no_tuning_marker & tuning{d1}~=no_anova_marker;
            for d2=1:numel(keys.cc.tasktypes)
                matrix_data(e).same(d1,d2)=sum(tuning{d1}(basic_quantity)==tuning{d2}(basic_quantity));
                matrix_data(e).diff(d1,d2)=sum(tuning{d1}(basic_quantity)~=tuning{d2}(basic_quantity) & tuning{d2}(basic_quantity)~=no_tuning_marker);
                matrix_data(e).none(d1,d2)=sum(tuning{d2}(basic_quantity)==no_tuning_marker);
                matrix_data(e).base(d1,d2)=sum(basic_quantity);
            end
        end
    end
elseif strcmp(keys.cc.plot_type,'per_task')
    for d=1:numel(keys.cc.tasktypes)
        for e=1:numel(keys.cc.epochs)
            epochtuning{e}= [keys.instructed_choice '_NH_' keys.cc.epochs{e} '_epoch_'   keys.cc.tasktypes{d} '_' keys.case{1} ]; 
            spacetuning{e}= [keys.instructed_choice '_' keys.cc.epochs{e} '_spaceLR_' keys.cc.tasktypes{d} '_' keys.case{1} ];
            handstuning{e}= [keys.instructed_choice '_' keys.cc.epochs{e} '_hands_'   keys.cc.tasktypes{d} '_' keys.case{1} ];
            S_x_Htuning{e}= [keys.instructed_choice '_' keys.cc.epochs{e} '_SxH_'     keys.cc.tasktypes{d} '_' keys.case{1} ];
            choictuning{e}= [keys.instructed_choice '_NH_' keys.cc.epochs{e} '_RF_' keys.cc.factors '_'    keys.cc.tasktypes{d} '_' keys.case{1} ];
            tuning{e}=read_table_column_detail(keys,xlsx_table,idx,epochtuning{e},spacetuning{e},handstuning{e},S_x_Htuning{e},choictuning{e});
        end
        pie_tmp=[];
        pie_tmp=get_pie_multilevel(keys,pie_tmp,tuning);
        multilevel_data(d,:)=pie_tmp;
        
        for e1=1:numel(keys.cc.epochs)
            basic_quantity=tuning{e1}~=no_tuning_marker & tuning{e1}~=no_anova_marker;
            for e2=1:numel(keys.cc.epochs)
                matrix_data(d).same(e1,e2)=sum(tuning{e1}(basic_quantity)==tuning{e2}(basic_quantity));
                matrix_data(d).diff(e1,e2)=sum(tuning{e1}(basic_quantity)~=tuning{e2}(basic_quantity) & tuning{e2}(basic_quantity)~=no_tuning_marker);
                matrix_data(d).none(e1,e2)=sum(tuning{e2}(basic_quantity)==no_tuning_marker);
                matrix_data(d).base(e1,e2)=sum(basic_quantity);
            end
        end
    end
elseif strcmp(keys.cc.plot_type,'space_and_epoch')
    for e=1:numel(keys.cc.epochs)
        %for d=1:numel(keys.cc.tasktypes)
        d=1;
        epochtuning{d}= [keys.instructed_choice '_NH_' keys.cc.epochs{e} '_epoch_'   keys.cc.tasktypes{d} '_' keys.case{1} ]; 
        spacetuning{d}= [keys.instructed_choice '_' keys.cc.epochs{e} '_spaceLR_' keys.cc.tasktypes{d} '_' keys.case{1} ];
        handstuning{d}= [keys.instructed_choice '_' keys.cc.epochs{e} '_hands_'   keys.cc.tasktypes{d} '_' keys.case{1} ];
        S_x_Htuning{d}= [keys.instructed_choice '_' keys.cc.epochs{e} '_SxH_'     keys.cc.tasktypes{d} '_' keys.case{1} ];
        choictuning{d}= [keys.instructed_choice '_NH_' keys.cc.epochs{e} '_RF_' keys.cc.factors '_'  keys.cc.tasktypes{d} '_' keys.case{1} ];
            
        keys.cc.factors='space';
        tunings{d}=read_table_column_detail(keys,xlsx_table,idx,epochtuning{d},spacetuning{d},handstuning{d},S_x_Htuning{d},choictuning{d});
        keys.cc.factors='epoch';
        tuninge{d}=read_table_column_detail(keys,xlsx_table,idx,epochtuning{d},spacetuning{d},handstuning{d},S_x_Htuning{d},choictuning{d});
        %end
        pie_tmp=[];
        pie_tmp=get_pie_two_levels(keys,pie_tmp,tunings,tuninge);
        multilevel_data(e,:)=pie_tmp;
    end
    
elseif strcmp(keys.cc.plot_type,'hand_and_epoch')
    for e=1:numel(keys.cc.epochs)
        %for d=1:numel(keys.cc.tasktypes)
        d=1;
        epochtuning{d}= [keys.instructed_choice '_NH_' keys.cc.epochs{e} '_epoch_'   keys.cc.tasktypes{d} '_' keys.case{1} ]; 
        spacetuning{d}= [keys.instructed_choice '_' keys.cc.epochs{e} '_spaceLR_' keys.cc.tasktypes{d} '_' keys.case{1} ];
        handstuning{d}= [keys.instructed_choice '_' keys.cc.epochs{e} '_hands_'   keys.cc.tasktypes{d} '_' keys.case{1} ];
        S_x_Htuning{d}= [keys.instructed_choice '_' keys.cc.epochs{e} '_SxH_'     keys.cc.tasktypes{d} '_' keys.case{1} ];
        choictuning{d}= [keys.instructed_choice '_NH_' keys.cc.epochs{e} '_RF_' keys.cc.factors '_' keys.cc.tasktypes{d} '_' keys.case{1} ];
        keys.cc.factors='hand';
        tuningh{d}=read_table_column_detail(keys,xlsx_table,idx,epochtuning{d},spacetuning{d},handstuning{d},S_x_Htuning{d},choictuning{d});
        keys.cc.factors='epoch';
        tuninge{d}=read_table_column_detail(keys,xlsx_table,idx,epochtuning{d},spacetuning{d},handstuning{d},S_x_Htuning{d},choictuning{d});
        %end
        pie_tmp=[];
        pie_tmp=get_pie_two_levels(keys,pie_tmp,tuningh,tuninge);
        multilevel_data(e,:)=pie_tmp;
    end
    
elseif strcmp(keys.cc.plot_type,'space_and_hand')
    for e=1:numel(keys.cc.epochs)
        %for d=1:numel(keys.cc.tasktypes)
        d=1;
        epochtuning{d}= [keys.instructed_choice '_NH_' keys.cc.epochs{e} '_epoch_'   keys.cc.tasktypes{d} '_' keys.case{1} ];
        spacetuning{d}= [keys.instructed_choice '_' keys.cc.epochs{e} '_spaceLR_' keys.cc.tasktypes{d} '_' keys.case{1} ];
        handstuning{d}= [keys.instructed_choice '_' keys.cc.epochs{e} '_hands_'   keys.cc.tasktypes{d} '_' keys.case{1} ];
        S_x_Htuning{d}= [keys.instructed_choice '_' keys.cc.epochs{e} '_SxH_'     keys.cc.tasktypes{d} '_' keys.case{1} ];
        choictuning{d}= [keys.instructed_choice '_NH_' keys.cc.epochs{e} '_RF_' keys.cc.factors '_' keys.cc.tasktypes{d} '_' keys.case{1} ];
        keys.cc.factors='space';
        tuningh{d}=read_table_column_detail(keys,xlsx_table,idx,epochtuning{d},spacetuning{d},handstuning{d},S_x_Htuning{d},choictuning{d});
        keys.cc.factors='hand';
        tuninge{d}=read_table_column_detail(keys,xlsx_table,idx,epochtuning{d},spacetuning{d},handstuning{d},S_x_Htuning{d},choictuning{d});
        %end
        pie_tmp=[];
        pie_tmp=get_pie_two_levels(keys,pie_tmp,tuningh,tuninge);
        multilevel_data(e,:)=pie_tmp;
    end
elseif strcmp(keys.cc.plot_type,'space_x_hand')
    for e=1:numel(keys.cc.epochs)
        %for d=1:numel(keys.cc.tasktypes)
        d=1;
        epochtuning{d}= [keys.instructed_choice '_NH_' keys.cc.epochs{e} '_epoch_'  keys.cc.tasktypes{d} '_' keys.case{1} ]; 
        spacetuning{d}= [keys.instructed_choice '_' keys.cc.epochs{e} '_spaceLR_' keys.cc.tasktypes{d} '_' keys.case{1} ];
        handstuning{d}= [keys.instructed_choice '_' keys.cc.epochs{e} '_hands_'   keys.cc.tasktypes{d} '_' keys.case{1} ];
        S_x_Htuning{d}= [keys.instructed_choice '_' keys.cc.epochs{e} '_SxH_'     keys.cc.tasktypes{d} '_' keys.case{1} ];
        choictuning{d}= [keys.instructed_choice '_NH_' keys.cc.epochs{e} '_RF_' keys.cc.factors '_'  keys.cc.tasktypes{d} '_' keys.case{1} ];
        keys.cc.factors='space_and_hand';
        tuningh{d}=read_table_column_detail(keys,xlsx_table,idx,epochtuning{d},spacetuning{d},handstuning{d},S_x_Htuning{d},choictuning{d});
        keys.cc.factors='space_x_hand';
        tuninge{d}=read_table_column_detail(keys,xlsx_table,idx,epochtuning{d},spacetuning{d},handstuning{d},S_x_Htuning{d},choictuning{d});
        %end
        pie_tmp=[];
        pie_tmp=get_pie_two_levels(keys,pie_tmp,tuningh,tuninge);
        multilevel_data(e,:)=pie_tmp;
    end
elseif strcmp(keys.cc.plot_type,'fixation_x_position')
    for e=1:numel(keys.cc.epochs)
        %for d=1:numel(keys.cc.tasktypes)
        d=1;
        epochtuning{d}= [keys.instructed_choice '_NH_' keys.cc.epochs{e} '_epoch_'      keys.cc.tasktypes{d} '_' keys.case{1} ];
        fixattuning{d}= [keys.instructed_choice '_NH_' keys.cc.epochs{e} '_fixation_'   keys.cc.tasktypes{d} '_' keys.case{1} ];
        posittuning{d}= [keys.instructed_choice '_NH_' keys.cc.epochs{e} '_position_'   keys.cc.tasktypes{d} '_' keys.case{1} ];
        F_x_Ptuning{d}= [keys.instructed_choice '_NH_' keys.cc.epochs{e} '_PxF_'        keys.cc.tasktypes{d} '_' keys.case{1} ];
        choictuning{d}= [keys.instructed_choice '_NH_' keys.cc.epochs{e} '_RF_' keys.cc.factors '_'  keys.cc.tasktypes{d} '_' keys.case{1} ];
        keys.cc.factors='fixation_x_position';
        tuningm{d}=read_table_column_detail(keys,xlsx_table,idx,epochtuning{d},fixattuning{d},posittuning{d},F_x_Ptuning{d},choictuning{d});
        epochtuning{d}= [keys.instructed_choice '_NH_' keys.cc.epochs{e} '_epoch_'      keys.cc.tasktypes{d} '_' keys.case{2} ];
        fixattuning{d}= [keys.instructed_choice '_NH_' keys.cc.epochs{e} '_fixation_'   keys.cc.tasktypes{d} '_' keys.case{2} ];
        posittuning{d}= [keys.instructed_choice '_NH_' keys.cc.epochs{e} '_position_'   keys.cc.tasktypes{d} '_' keys.case{2} ];
        F_x_Ptuning{d}= [keys.instructed_choice '_NH_' keys.cc.epochs{e} '_PxF_'        keys.cc.tasktypes{d} '_' keys.case{2} ];   
        choictuning{d}= [keys.instructed_choice '_NH_' keys.cc.epochs{e} '_RF_' keys.cc.factors '_'  keys.cc.tasktypes{d} '_' keys.case{1} ];
        keys.cc.factors='fixation_x_position';
        tuningt{d}=read_table_column_detail(keys,xlsx_table,idx,epochtuning{d},fixattuning{d},posittuning{d},F_x_Ptuning{d},choictuning{d});
        %end
        pie_tmp=[];
        pie_tmp=get_pie_two_levels(keys,pie_tmp,tuningm,tuningt);
        multilevel_data(e,:)=pie_tmp;
    end
end

function read_out=read_table_column_detail(keys,table,index,indexfnE,indexfnS,indexfnH,indexfnX,indexfnC)
% 
valid1=isfield(index,indexfnX) && isfield(index,indexfnH) && isfield(index,indexfnS) &&...
    any(index.(indexfnX))   && any(index.(indexfnH))   && any(index.(indexfnS)) && strcmp(keys.cc.factors,'space_and_hand');
valid2=isfield(index,indexfnH) && any(index.(indexfnH)) && strcmp(keys.cc.factors,'hand');
valid3=isfield(index,indexfnS) && any(index.(indexfnS)) && strcmp(keys.cc.factors,'space');
valid4=isfield(index,indexfnE) && any(index.(indexfnE)) && strcmp(keys.cc.factors,'epoch');
valid5=isfield(index,indexfnX) && any(index.(indexfnX)) && strcmp(keys.cc.factors,'space_x_hand');
valid6=isfield(index,indexfnX) && any(index.(indexfnX)) && strcmp(keys.cc.factors,'fixation_x_position'); %% workaround, indexfnS=fixaion, indexfnH=position, indexfnX=pxF,
valid7=isfield(index,indexfnC) && any(index.(indexfnC)) && (strcmp(keys.cc.factors,'choice1') || strcmp(keys.cc.factors,'choice2')); 

if valid1 %space x hand
    IHtuning=ismember(table(2:end,index.(indexfnH)),'IH') & keys.SXH_or_interaction; %keys.hands_or_interaction;
    CHtuning=ismember(table(2:end,index.(indexfnH)),'CH') & keys.SXH_or_interaction; %keys.hands_or_interaction;
    IStuning=ismember(table(2:end,index.(indexfnS)),'IS') & keys.SXH_or_interaction; %keys.space_or_interaction;
    CStuning=ismember(table(2:end,index.(indexfnS)),'CS') & keys.SXH_or_interaction; %keys.space_or_interaction;
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
elseif  valid2 %% hand
    IH=ismember(table(2:end,index.(indexfnH)),'IH') & keys.hands_or_interaction;
    CH=ismember(table(2:end,index.(indexfnH)),'CH') & keys.hands_or_interaction;
    na=~keys.hands_or_interaction;
    clear read_out
    read_out(IH)                =keys.ipsi_index;
    read_out(CH)                =keys.contra_color_index;
    read_out(na)                =2;
    read_out(~IH & ~CH & ~na)   =3;
    
elseif valid3 %% space
    IS=ismember(table(2:end,index.(indexfnS)),'IS') & keys.space_or_interaction;
    CS=ismember(table(2:end,index.(indexfnS)),'CS') & keys.space_or_interaction;
    na=~keys.space_or_interaction;
    clear read_out
    read_out(IS)                =keys.ipsi_index;
    read_out(CS)                =keys.contra_color_index;
    read_out(na)                =2;
    read_out(~IS & ~CS & ~na)   =3;
    
elseif valid4 %% epoch
    en=ismember(table(2:end,index.(indexfnE)),'en') & keys.epoch_or_interaction;
    su=ismember(table(2:end,index.(indexfnE)),'su') & keys.epoch_or_interaction;
    bi=ismember(table(2:end,index.(indexfnE)),'bi') & keys.epoch_or_interaction;
    na=~keys.epoch_or_interaction;
    clear read_out
    read_out(en)                =1;
    read_out(bi)                =2;
    read_out(su)                =3;
    read_out(na)                =4;
    read_out(~en & ~su & ~bi & ~na)   =5;
    
elseif valid5 %% space X hand
    CR=ismember(table(2:end,index.(indexfnX)),'CR') & keys.SXH_or_interaction;
    UC=ismember(table(2:end,index.(indexfnX)),'UC') & keys.SXH_or_interaction;
    na=~keys.SXH_or_interaction;
    clear read_out
    read_out(UC)                =4;
    read_out(CR)                =1;
    read_out(na)                =2;
    read_out(~CR & ~UC & ~na)   =3;
    
elseif valid6 %% position x fixation
    %indexfnS=fixaion, indexfnH=position, indexfnX=pxF
    fixtuning=ismember(table(2:end,index.(indexfnS)),'true'); %& keys.hands_or_interaction; %% anova criterion for this one doesnt make much sense, does it??
    postuning=ismember(table(2:end,index.(indexfnH)),'true'); %& keys.hands_or_interaction;
    fxptuning=ismember(table(2:end,index.(indexfnX)),'true'); %& keys.space_or_interaction;
    %CStuning=ismember(table(2:end,index.(indexfnS)),'CS') & keys.space_or_interaction;
    Fix     =  fixtuning & ~postuning & ~fxptuning;
    Pos     = ~fixtuning &  postuning & ~fxptuning;
    FixI    =  fixtuning & ~postuning &  fxptuning;
    PosI    = ~fixtuning &  postuning &  fxptuning;
    FP      =  fixtuning &  postuning & ~fxptuning;
    FPI     =  fixtuning &  postuning &  fxptuning;
    I       = ~fixtuning & ~postuning &  fxptuning;
    Na      = ~fixtuning & ~postuning & ~fxptuning;
    
    read_out(Fix)        =1;
    read_out(Pos)        =2;
    read_out(FixI)       =3;
    read_out(PosI)       =4;
    read_out(FP)         =5;
    read_out(FPI)        =6;
    read_out(I)          =7;
    read_out(Na)         =8;
elseif valid7 %choice
    IN=ismember(table(2:end,index.(indexfnC)),'IN') & keys.space_or_interaction;
    CH=ismember(table(2:end,index.(indexfnC)),'CH') & keys.space_or_interaction;
    na=~keys.space_or_interaction;
    clear read_out
    read_out(IN)                =1;
    read_out(CH)                =4;
    read_out(na)                =2;
    read_out(~IN & ~CH & ~na)   =3;
    
else
    if strcmp(keys.cc.factors,'hand') || strcmp(keys.cc.factors,'space') || strcmp(keys.cc.factors,'epoch')
        read_out=ones(1,size(table,1)-1)*2;
    else
        read_out=ones(1,size(table,1)-1)*10;
    end
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
n_conditions=10;
if strcmp(keys.cc.factors,'hand') || strcmp(keys.cc.factors,'space') ||...
   strcmp(keys.cc.factors,'choice1') || strcmp(keys.cc.factors,'choice2')
    n_conditions=4;
elseif strcmp(keys.cc.factors,'epoch')
    n_conditions=5;
end
for c=1:n_conditions
    next_level_index = tuning{level}==c;
    pie{level}(end+1)= sum(previous_level_index(:) & next_level_index(:));
    if level<numel(tuning)
        pie=get_pie_multilevel(keys,pie,tuning,previous_level_index(:) & next_level_index(:),level);
    end
end

function pie=get_pie_two_levels(keys,pie,tuning1,tuning2)
if ismember(keys.cc.plot_type,{'space_and_hand','hand_and_epoch'}) 
    nt1=4;
    nt2=4;    
elseif ismember(keys.cc.plot_type,{'space_and_epoch'})
    nt1=4;
    nt2=5;
elseif strcmp(keys.cc.plot_type,'space_x_hand')
    nt1=10;
    nt2=4;
elseif strcmp(keys.cc.plot_type,'fixation_x_position')
    nt1=8;
    nt2=8;
end
pie{1}=[];
pie{2}=[];
for n1=1:nt1
    pie{1}(end+1)=sum(tuning1{1}==n1);
    for n2=1:nt2
        pie{2}(end+1)=sum(tuning1{1}(:)==n1 & tuning2{1}(:)==n2);
    end
end

function title_and_save(keys)
% selected=cellfun(@num2str,keys.tt.selection(:,2),'Uniformoutput',false);
% if ~isempty(keys.tt.selection)
%     Selection=strcat(keys.tt.selection(:,1),'= ', selected, ',')';
% else
%     Selection={};
% end
criterions=[' criterions: S: ' keys.tt.space_criterion ' H: ' keys.tt.hands_criterion ' SxH: ' keys.tt.SXH_criterion ' E: ' keys.tt.epoch_criterion];
mtit(gcf,[keys.monkey ' ' [keys.cc.tasktypes{:}] ' ' [keys.cc.plot_type{:}] ' ' keys.cc.factors ' ' keys.instructed_choice ' ' keys.case{1} ' ' keys.selection_title{:} ' N: ' num2str(size(keys.xlsx_table,1)-1) ' ' keys.title_part  ' ' criterions], 'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none');

folder_to_save=[keys.drive filesep keys.basepath_to_save filesep keys.pdf_folder];
if ~exist([folder_to_save filesep keys.subfolder_to_save],'dir');
    mkdir(folder_to_save,keys.subfolder_to_save);
end
export_fig(gcf, [folder_to_save filesep keys.subfolder_to_save filesep keys.monkey ' ' [keys.cc.tasktypes{:}] ' ' keys.instructed_choice ' '...
    keys.selection_title{:} ' ' [keys.cc.plot_type{:}]  ' ' keys.cc.factors ', ' keys.title_part], '-pdf','-transparent') % pdf by run
close(gcf);




