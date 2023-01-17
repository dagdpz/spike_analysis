function out=ph_scatter(modified_keys)
out=struct;
fontsize=6;

for fn=fieldnames(modified_keys)'
    keys.(fn{:})=modified_keys.(fn{:});
end

TT=keys.tuning_table;
if size(TT,1)==1
    return;
end

% %% ANOVA cirterions readout
% criterions={'space_or_interaction','epoch_or_interaction','hands_or_interaction','SXH_or_interaction'};
% for c=1:numel(criterions)
%     parameter_criterion_columns=~cellfun(@isempty,strfind(TT(1,:),criterions{c}));
%     for t=1:numel(keys.tt.tasktypes)
%         task_criterion_columns(t,:)=~cellfun(@isempty,strfind(TT(1,:),keys.tt.tasktypes{t}));
%     end
%     criterion_columns=parameter_criterion_columns & any(task_criterion_columns,1);
%     keys.(criterions{c})=any(cell2mat(TT(2:end,criterion_columns)),2);
% end

%% row indexes
x_label=keys.SC.X;
y_label=keys.SC.Y;

c_x    =DAG_find_column_index(TT,keys.SC.X);
c_y    =DAG_find_column_index(TT,keys.SC.Y);
cs_x   =DAG_find_column_index(TT,keys.SC.X_sig);
cs_y   =DAG_find_column_index(TT,keys.SC.Y_sig);
c_VM   =DAG_find_column_index(TT,keys.SC.VMI);
c_h    =DAG_find_column_index(TT,keys.SC.hist_column);

rs_x   =[false; cellfun(@(x) ~strcmp(x,'-') && ~strcmp(x,'--') &&  ~strcmp(x,'false'),TT(2:end,cs_x))];
rs_y   =[false; cellfun(@(x) ~strcmp(x,'-') && ~strcmp(x,'--') &&  ~strcmp(x,'false'),TT(2:end,cs_y))];


row_valid=cellfun(@(x) ~isempty(x) && ~isnan(x),TT(2:end,c_x)) & cellfun(@(x) ~isempty(x) && ~isnan(x),TT(2:end,c_y));
row_valid=[false;row_valid];

%% monkey per unit readout
col_unit_ID=DAG_find_column_index(TT,'unit_ID');
monkeys=cellfun(@(x) x(1:3),TT(2:end,col_unit_ID),'uniformoutput',false);
unique_monkeys=unique(monkeys);
monkeys=['___'; monkeys];


if ~any(row_valid)
    return;
end

rs_xy       = rs_x &  rs_y & row_valid;
rs_xn       = rs_x & ~rs_y & row_valid;
rs_ny       =~rs_x &  rs_y & row_valid;
rs_nn       =~rs_x & ~rs_y & row_valid;

if ismember(keys.SC.color_option,{'ENSU_as_color'})
    
end
%

if ismember(keys.SC.color_option,{'monkeys_by_marker'})
    
    cat_colors={'m','b','r',[0.5 0.5 0.5]}; %% easier to understand significant in both is the mixed color
    row_nocat=row_valid;
    for idx=1:3
        switch idx
            case 1
                r_cat{idx}= rs_x & ~rs_y & row_valid;
            case 2
                r_cat{idx}= rs_x &  rs_y & row_valid;
            case 3
                r_cat{idx}=~rs_x & ~rs_y & row_valid;
        end
        row_nocat =row_nocat & ~r_cat{idx};
    end
    r_cat{idx+1}=row_nocat;
end


if ismember(keys.SC.color_option,{'FR_as_color','VMI_as_color','category_as_color'})
    
    cat_colors={'b','m','r',[0.5 0.5 0.5]}; %% based on VMI still
    row_nocat=row_valid;
    for idx=1:numel(keys.SC.categories)
        col_cat   =DAG_find_column_index(TT,keys.SC.categories{idx});
        r_cat{idx}=[false; vertcat(TT{2:end,col_cat})]&row_valid;
        row_nocat =row_nocat & ~r_cat{idx};
    end
    r_cat{idx+1}=row_nocat;
end


%% monkey markers and columns assignment
for m=1:numel(unique_monkeys)
    row_monkey{m}=ismember(monkeys,unique_monkeys{m})& row_valid;
end

%% figure 1
fig_title=ph_figure_titles(keys);
%fig_title=[keys.monkey ' ' keys.selection_title{:} y_label '__vs__' x_label];
plot_1_title            = [fig_title  y_label '__vs__' x_label];
FR_index_handle = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_1_title);

if keys.SC.logarithmic_scale%%% transform to logarithmic scale
    temp_x_nonempty=TT(row_valid,c_x);
    temp_y_nonempty=TT(row_valid,c_y);
    TT(row_valid,c_x)=num2cell(sign([TT{row_valid,c_x}]).*log(abs([TT{row_valid,c_x}])+1));
    TT(row_valid,c_y)=num2cell(sign([TT{row_valid,c_y}]).*log(abs([TT{row_valid,c_y}])+1));
end
if keys.SC.absolutes%%% transform to logarithmic scale
    temp_x_nonempty=TT(row_valid,c_x);
    temp_y_nonempty=TT(row_valid,c_y);
    TT(row_valid,c_x)=num2cell(abs([TT{row_valid,c_x}]));
    TT(row_valid,c_y)=num2cell(abs([TT{row_valid,c_y}]));
end

SC.markers={};
SC.filled={};
SC.colors={};
SC.hist_colors={};
SC.x_values={};
SC.y_values={};
switch keys.SC.color_option
    case 'monkeys_by_color'
        for idx=1:numel(unique_monkeys)
            idx1=2*idx-1;
            idx2=idx1+1;
            idxd=3*(idx-1);
            current_monkey=keys.monkeys{cellfun(@(x) any(strfind(x,unique_monkeys{idx})),keys.monkeys)};
            current_color=keys.(current_monkey).color;
            SC.markers=[SC.markers {'o', '>', '^', 'o'}];
            SC.filled=[SC.filled {'filled', 'filled', 'filled', ''}];
            SC.colors=[SC.colors current_color current_color current_color current_color];
            SC.hist_colors{1}{idx1}=current_color;
            SC.hist_colors{2}{idx1}=current_color;
            SC.hist_colors{1}{idx2}=current_color/2;
            SC.hist_colors{2}{idx2}=current_color/2;
            SC.hist_colors{3}(idxd+1:idxd+3)={current_color current_color/2 [0.5 0.5 0.5]};
            r_m=row_monkey{idx};
            
            SC.x_values{1,idx1}=[{[TT{rs_xy & r_m,c_x}]} {[TT{rs_xn & r_m,c_x}]}];
            SC.x_values{1,idx2}=[{[TT{rs_ny & r_m,c_x}]} {[TT{rs_nn & r_m,c_x}]}];
            SC.y_values{idx1,1}=[{[TT{rs_xy & r_m,c_y}]} {[TT{rs_xn & r_m,c_y}]}];
            SC.y_values{idx2,1}=[{[TT{rs_ny & r_m,c_y}]} {[TT{rs_nn & r_m,c_y}]}];
            
            SC.d_values{idxd+1}= [{[TT{rs_xy & r_m,c_x}]}  {[TT{rs_xy & r_m,c_y}]}];
            SC.d_values{idxd+2}= [{[TT{(rs_xn | rs_ny) & r_m,c_x}]}  {[TT{(rs_xn | rs_ny) & r_m,c_y}]}];
            SC.d_values{idxd+3}= [{[TT{rs_nn & r_m,c_x}]}  {[TT{rs_nn & r_m,c_y}]}];
        end
        
    case 'monkeys_by_marker'
        r_am=false(size(row_monkey{1}));
        for idx=1:numel(unique_monkeys)
            idx1=2*idx-1;
            idx2=idx1+1;
            
            current_monkey=keys.monkeys{cellfun(@(x) any(strfind(x,unique_monkeys{idx})),keys.monkeys)};
            current_marker=keys.(current_monkey).marker;
            SC.markers=[SC.markers {current_marker} {current_marker} {current_marker} {current_marker}];
            SC.filled=[SC.filled; {'filled', 'filled', 'filled', ''}];
            SC.colors=[SC.colors; cat_colors];
            r_m=row_monkey{idx};
            r_am=r_am | r_m;
            SC.x_values{1,1}{idx}=[TT{rs_xy & r_m,c_x}];
            SC.x_values{1,2}{idx}=[TT{rs_xn & r_m,c_x}];
            SC.x_values{1,3}{idx}=[TT{rs_ny & r_m,c_x}];
            SC.x_values{1,4}{idx}=[TT{rs_nn & r_m,c_x}];
            
            SC.y_values{1,1}{idx}=[TT{rs_xy & r_m,c_y}];
            SC.y_values{2,1}{idx}=[TT{rs_xn & r_m,c_y}];
            SC.y_values{3,1}{idx}=[TT{rs_ny & r_m,c_y}];
            SC.y_values{4,1}{idx}=[TT{rs_nn & r_m,c_y}];
        end
        SC.colors=SC.colors(:);
        SC.filled=SC.filled(:);
        SC.hist_colors{1}=cat_colors;
        SC.hist_colors{2}=cat_colors;
        SC.hist_colors{3}=cat_colors;
        SC.d_values{1}= [{[TT{rs_xy & r_am,c_x}]}  {[TT{rs_xy & r_am,c_y}]}];
        SC.d_values{2}= [{[TT{rs_xn & r_am,c_x}]}  {[TT{rs_xn & r_am,c_y}]}];
        SC.d_values{3}= [{[TT{rs_ny & r_am,c_x}]}  {[TT{rs_ny & r_am,c_y}]}];
        SC.d_values{4}= [{[TT{rs_nn & r_am,c_x}]}  {[TT{rs_nn & r_am,c_y}]}];
        
        
    case 'FR_as_color'
        colormap([(1:1:255)' zeros(255,1) (255:-1:1)']/255)
        FR_values=[TT{row_valid,cs_x}]; %%% only by X???
        maxFR=max(FR_values);
        FR_values=round(FR_values/maxFR*255);
        for idx=1:numel(r_cat)
            SC.markers=[SC.markers {'o', '>', '^', 'o'}];
            SC.filled=[SC.filled {'filled', 'filled', 'filled', ''}];
            SC.colors=[SC.colors {FR_values(rs_xy(row_valid) & r_cat{idx}(row_valid))'} {FR_values(rs_xn(row_valid) & r_cat{idx}(row_valid))'} {FR_values(rs_ny(row_valid) & r_cat{idx}(row_valid))'} {FR_values(rs_nn(row_valid) & r_cat{idx}(row_valid))'}];
            SC.hist_colors{1}{idx}=cat_colors{idx};
            SC.hist_colors{2}{idx}=cat_colors{idx};
            SC.hist_colors{3}{idx}=cat_colors{idx};
            SC.x_values{idx}=[{[TT{rs_xy & r_cat{idx},c_x}]} {[TT{rs_xn & r_cat{idx},c_x}]} {[TT{rs_ny & r_cat{idx},c_x}]} {[TT{rs_nn & r_cat{idx},c_x}]}];
            SC.y_values{idx}=[{[TT{rs_xy & r_cat{idx},c_y}]} {[TT{rs_xn & r_cat{idx},c_y}]} {[TT{rs_ny & r_cat{idx},c_y}]} {[TT{rs_nn & r_cat{idx},c_y}]}];
            SC.d_values{idx}=[{[SC.x_values{idx}{:}]} {[SC.y_values{idx}{:}]}];
        end
        
    case 'VMI_as_color'
        SC.markers={'o', '>', '^', 'o'};
        SC.filled={'filled', 'filled', 'filled', ''};
        VM_values=[TT{row_valid,c_VM}];
        VM_values=round((VM_values+1)/2*255);
        
        for idx=1:numel(r_cat)
            SC.markers=[SC.markers {'o', '>', '^', 'o'}];
            SC.filled=[SC.filled {'filled', 'filled', 'filled', ''}];
            SC.colors=[SC.colors {VM_values(rs_xy(row_valid) & r_cat{idx}(row_valid))'} {VM_values(rs_xn(row_valid) & r_cat{idx}(row_valid))'} {VM_values(rs_ny(row_valid) & r_cat{idx}(row_valid))'} {VM_values(rs_nn(row_valid) & r_cat{idx}(row_valid))'}];
            SC.hist_colors{1}{idx}=cat_colors{idx};
            SC.hist_colors{2}{idx}=cat_colors{idx};
            SC.hist_colors{3}{idx}=cat_colors{idx};
            SC.x_values{idx}=[{[TT{rs_xy & r_cat{idx},c_x}]} {[TT{rs_xn & r_cat{idx},c_x}]} {[TT{rs_ny & r_cat{idx},c_x}]} {[TT{rs_nn & r_cat{idx},c_x}]}];
            SC.y_values{idx}=[{[TT{rs_xy & r_cat{idx},c_y}]} {[TT{rs_xn & r_cat{idx},c_y}]} {[TT{rs_ny & r_cat{idx},c_y}]} {[TT{rs_nn & r_cat{idx},c_y}]}];
            SC.d_values{idx}=[{[SC.x_values{idx}{:}]} {[SC.y_values{idx}{:}]}];
        end
        colormap([(1:1:255)' zeros(255,1) (255:-1:1)']/255)
        
    case 'category_as_color'
        for idx=1:numel(r_cat)
            SC.markers=[SC.markers {'o', '>', '^', 'o'}];
            SC.filled=[SC.filled {'filled', 'filled', 'filled', ''}];
            SC.colors=[SC.colors cat_colors(idx) cat_colors(idx) cat_colors(idx) cat_colors(idx)];
            SC.hist_colors{1}{idx}=cat_colors{idx};
            SC.hist_colors{2}{idx}=cat_colors{idx};
            SC.hist_colors{3}{idx}=cat_colors{idx};
            SC.x_values{idx}=[{[TT{rs_xy & r_cat{idx},c_x}]} {[TT{rs_xn & r_cat{idx},c_x}]} {[TT{rs_ny & r_cat{idx},c_x}]} {[TT{rs_nn & r_cat{idx},c_x}]}];
            SC.y_values{idx}=[{[TT{rs_xy & r_cat{idx},c_y}]} {[TT{rs_xn & r_cat{idx},c_y}]} {[TT{rs_ny & r_cat{idx},c_y}]} {[TT{rs_nn & r_cat{idx},c_y}]}];
            SC.d_values{idx}=[{[SC.x_values{idx}{:}]} {[SC.y_values{idx}{:}]}];
        end
        
    case 'ENSU_as_color'
        ensu_x=TT(:,cs_x);ensu_x(~row_valid)={''};
        ensu_y=TT(:,cs_y);ensu_y(~row_valid)={''};
        rE_x=ismember(ensu_x,'en') & row_valid;
        rE_y=ismember(ensu_y,'en')  & row_valid;
        rS_x=ismember(ensu_x,'su')  & row_valid;
        rS_y=ismember(ensu_y,'su')  & row_valid;
        
        SC.markers={'o', 'o', '>', 'o', 'o', '>', '^', '^', 'o'};
        SC.colors={'r', 'm', 'r', 'c', 'b', 'b', 'r', 'b', [0.5 0.5 0.5]};
        SC.filled={'filled', 'filled', 'filled', 'filled', 'filled', 'filled', 'filled', 'filled', ''};
        
        SC.hist_colors{1}={'r', 'b', [0.5 0.5 0.5]};
        SC.hist_colors{2}={'r', 'b', [0.5 0.5 0.5]};
        SC.hist_colors{3}={'r', 'b', 'm', 'c', [0.5 0.5 0.5]};
        
        SC.x_values(:,1) =[{{[TT{rs_xy & rE_x & rE_y,c_x}]}} {{[TT{rs_xy & rE_x & rS_y,c_x}]}} {{[TT{rs_xn & rE_x,c_x}]}}]; % x red (y red/blue/black)
        SC.x_values(:,2) =[{{[TT{rs_xy & rS_x & rE_y,c_x}]}} {{[TT{rs_xy & rS_x & rS_y,c_x}]}} {{[TT{rs_xn & rS_x,c_x}]}}]; % x blue (y red/blue/black)
        SC.x_values(:,3) =[{{[TT{rs_ny & rE_y,c_x}]}}        {{[TT{rs_ny & rS_y,c_x}]}}        {{[TT{rs_nn,c_x}]}}]; % x black (y red/blue/black)
        SC.y_values(1,:)=[{{[TT{rs_xy & rE_x & rE_y,c_y}]}} {{[TT{rs_xy & rS_x & rE_y,c_y}]}} {{[TT{rs_ny & rE_y,c_y}]}}]; % y red (x red/blue/black)
        SC.y_values(2,:)=[{{[TT{rs_xy & rE_x & rS_y,c_y}]}} {{[TT{rs_xy & rS_x & rS_y,c_y}]}} {{[TT{rs_ny & rS_y,c_y}]}}]; % x blue (x red/blue/black)
        SC.y_values(3,:)=[{{[TT{rs_xn & rE_x,c_y}]}}        {{[TT{rs_xn & rS_x,c_y}]}}        {{[TT{rs_nn,c_y}]}}]; % x black (x red/blue/black)
        
        r_any_EN=rs_xy & rE_x & rE_y | rs_xn & rE_x | rs_ny & rE_y;
        r_any_SU=rs_xy & rS_x & rS_y | rs_xn & rS_x | rs_ny & rS_y;
        SC.d_values{1}= [{[TT{r_any_EN,c_x}]}  {[TT{r_any_EN,c_y}]}];
        SC.d_values{2}= [{[TT{r_any_SU,c_x}]}  {[TT{r_any_SU,c_y}]}];
        SC.d_values{3}= [{[TT{rs_xy & rE_x & rS_y,c_x}]}  {[TT{rs_xy & rE_x & rS_y,c_y}]}];
        SC.d_values{4}= [{[TT{rs_xy & rS_x & rE_y,c_x}]}  {[TT{rs_xy & rS_x & rE_y,c_y}]}];
        SC.d_values{5}= [{[TT{rs_nn,c_x}]}  {[TT{rs_nn,c_y}]}];
        
    case 'ENSU_all_separated'
        ensu_x=TT(:,cs_x);ensu_x(~row_valid)={''};
        ensu_y=TT(:,cs_y);ensu_y(~row_valid)={''};
        rE_x=ismember(ensu_x,'en') & row_valid;
        rE_y=ismember(ensu_y,'en')  & row_valid;
        rS_x=ismember(ensu_x,'su')  & row_valid;
        rS_y=ismember(ensu_y,'su')  & row_valid;
        
        SC.markers={'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o'};
        SC.colors={'r', 'c', [0.5 0 0], 'c', 'b', [0 0.5 0], 'm', 'g', [0.5 0.5 0.5]};
        SC.filled={'filled', 'filled', 'filled', 'filled', 'filled', 'filled', 'filled', 'filled', ''};
        
        SC.hist_colors{1}={'r', 'b', [0.5 0.5 0.5]};
        SC.hist_colors{2}={'r', 'b', [0.5 0.5 0.5]};
        %SC.hist_colors{3}={'r', 'b', 'm', 'c', [0.5 0.5 0.5]};
        SC.hist_colors{3}={'r', 'c', [0.5 0 0], 'c', 'b', [0 0.5 0], 'm', 'g', [0.5 0.5 0.5]};
        
        SC.x_values(:,1) =[{{[TT{rs_xy & rE_x & rE_y,c_x}]}} {{[TT{rs_xy & rE_x & rS_y,c_x}]}} {{[TT{rs_xn & rE_x,c_x}]}}]; % x red (y red/blue/black)
        SC.x_values(:,2) =[{{[TT{rs_xy & rS_x & rE_y,c_x}]}} {{[TT{rs_xy & rS_x & rS_y,c_x}]}} {{[TT{rs_xn & rS_x,c_x}]}}]; % x blue (y red/blue/black)
        SC.x_values(:,3) =[{{[TT{rs_ny & rE_y,c_x}]}}        {{[TT{rs_ny & rS_y,c_x}]}}        {{[TT{rs_nn,c_x}]}}]; % x black (y red/blue/black)
        SC.y_values(1,:)=[{{[TT{rs_xy & rE_x & rE_y,c_y}]}} {{[TT{rs_xy & rS_x & rE_y,c_y}]}} {{[TT{rs_ny & rE_y,c_y}]}}]; % y red (x red/blue/black)
        SC.y_values(2,:)=[{{[TT{rs_xy & rE_x & rS_y,c_y}]}} {{[TT{rs_xy & rS_x & rS_y,c_y}]}} {{[TT{rs_ny & rS_y,c_y}]}}]; % x blue (x red/blue/black)
        SC.y_values(3,:)=[{{[TT{rs_xn & rE_x,c_y}]}}        {{[TT{rs_xn & rS_x,c_y}]}}        {{[TT{rs_nn,c_y}]}}]; % x black (x red/blue/black)
        
        r_any_EN=rs_xy & rE_x & rE_y | rs_xn & rE_x | rs_ny & rE_y;
        r_any_SU=rs_xy & rS_x & rS_y | rs_xn & rS_x | rs_ny & rS_y;
        SC.d_values{1}= [{[TT{rs_xy & rE_x & rE_y,c_x}]}  {[TT{rs_xy & rE_x & rE_y,c_y}]}];
        SC.d_values{2}= [{[TT{rs_xy & rE_x & rS_y,c_x}]}  {[TT{rs_xy & rE_x & rS_y,c_y}]}];
        SC.d_values{3}= [{[TT{rs_xn & rE_x,c_x}]}  {[TT{rs_xn & rE_x,c_y}]}];
        SC.d_values{4}= [{[TT{rs_xy & rS_x & rE_y,c_x}]}  {[TT{rs_xy & rS_x & rE_y,c_y}]}];
        SC.d_values{5}= [{[TT{rs_xy & rS_x & rS_y,c_x}]}  {[TT{rs_xy & rS_x & rS_y,c_y}]}];
        SC.d_values{6}= [{[TT{rs_xn & rS_x,c_x}]}  {[TT{rs_xn & rS_x,c_y}]}];
        SC.d_values{7}= [{[TT{rs_ny & rE_y,c_x}]}  {[TT{rs_ny & rE_y,c_y}]}];
        SC.d_values{8}= [{[TT{rs_ny & rS_y,c_x}]}  {[TT{rs_ny & rS_y,c_y}]}];
        SC.d_values{9}= [{[TT{rs_nn,c_x}]}  {[TT{rs_nn,c_y}]}];
        
        
%         
%         
%         SC.d_values{2}= [{[TT{r_any_SU,c_x}]}  {[TT{r_any_SU,c_y}]}];
%         SC.d_values{3}= [{[TT{rs_xy & rE_x & rS_y,c_x}]}  {[TT{rs_xy & rE_x & rS_y,c_y}]}];
%         SC.d_values{4}= [{[TT{rs_xy & rS_x & rE_y,c_x}]}  {[TT{rs_xy & rS_x & rE_y,c_y}]}];
%         SC.d_values{5}= [{[TT{rs_nn,c_x}]}  {[TT{rs_nn,c_y}]}];
        
end

%FR_index_handle = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_1_title);

subplot(1,2,1)
DAG_scatter_with_histograms(SC)
xlabel(x_label,'interpreter','none');
ylabel(y_label,'interpreter','none');
if strcmp(keys.SC.color_option,'VMI_as_color')
    colorbar('yTick',[0,63,127,190,255],'yTickLabel',[-1,-0.5,0,0.5,1]); % doesnt work??
end

axes_limits=get(gca,'xlim'); %% xlim, but not ylim?
if keys.SC.logarithmic_scale %log scale strikes back
    realmax=exp(max(abs(axes_limits)))-1;
    log_axes_values=0:round(realmax/5):realmax;
    log_axes_ticks=log(log_axes_values+1);
    log_axes_values=[-1*fliplr(log_axes_values(2:end)) log_axes_values];
    log_axes_ticks=[-1*fliplr(log_axes_ticks(2:end)) log_axes_ticks];
    set(gca,'xtick',log_axes_ticks,'xticklabel',log_axes_values,'ytick',log_axes_ticks,'yticklabel',log_axes_values);
    % resetting table entries
    TT(row_valid,c_x)=   temp_x_nonempty;
    TT(row_valid,c_y)=   temp_y_nonempty;
end


if ~isempty(keys.SC.hist_column)
    %% add histogram for separately defined parameter
    
    subplot(1,2,2);
    hold on;
    switch keys.SC.color_option
        case 'monkeys_by_color'
        case 'FR_as_color'
            for idx=1:numel(r_cat)
                hist_values{idx}=[TT{r_cat{idx},c_h}];
                hist_colors{idx}=cat_colors{idx};
            end
        case 'VMI_as_color'
            for idx=1:numel(r_cat)
                hist_values{idx}=[TT{r_cat{idx},c_h}];
                hist_colors{idx}=cat_colors{idx};
            end
        case 'category_as_color'
            for idx=1:numel(r_cat)
                hist_values{idx}=[TT{r_cat{idx},c_h}];
                hist_colors{idx}=cat_colors{idx};
            end
    end
    minmax=[min([hist_values{:}]) max([hist_values{:}])];
    minmax=round(minmax);
    minmax=[-1*max(abs(minmax)) max(abs(minmax))];
    bins=minmax(1):diff(minmax)/20:minmax(2);
    hprev=zeros(size(bins));
    for idx=1:numel(hist_values)
        h=hist(hist_values{idx},bins)+hprev;
        mean_val{idx}=nanmean(hist_values{idx});
        [~,p_val{idx}]=ttest(hist_values{idx});
        Ns{idx}=numel(hist_values{idx});
        hb=bar(bins,h,'FaceColor',hist_colors{idx});
        uistack(hb,'bottom');
        hprev=h;
    end
    mean_val{idx+1}=nanmean([hist_values{:}]);
    [~,p_val{idx+1}]=ttest([hist_values{:}]);
    Ns{idx+1}=numel([hist_values{:}]);
    hist_colors{idx+1}='k';
    y_lim=get(gca,'ylim');
    x_lim=get(gca,'xlim');
    
    row_inc=-1*diff(y_lim)/40;
    col_inc=diff(x_lim)/10;
    row_start=y_lim(2)+row_inc;
    col_start=x_lim(1)+col_inc;
    text(col_start,row_start,          'N=','color','k');
    text(col_start,row_start+row_inc,  'M=','color','k');
    text(col_start,row_start+2*row_inc,'p=','color','k');
    for idx=1:numel(mean_val)
        text(col_start+col_inc*idx,row_start,          sprintf('%d',  Ns{idx}),       'color',hist_colors{idx});
        text(col_start+col_inc*idx,row_start+row_inc,  sprintf('%0.2f',mean_val{idx}),'color',hist_colors{idx});
        text(col_start+col_inc*idx,row_start+2*row_inc,sprintf('%0.3f',p_val{idx}),   'color',hist_colors{idx});
        line([mean_val{idx} mean_val{idx}],y_lim,'linestyle',':','color',hist_colors{idx})
    end
    
    ylabel('N units');
    xlabel(keys.SC.hist_column,'interpreter','none');
    
end

save([keys.basepath_to_save, keys.project_version, filesep, 'population_meta_data', filesep, 'scatter', filesep, plot_1_title], 'keys');
ph_title_and_save(FR_index_handle,plot_1_title,plot_1_title,keys);
