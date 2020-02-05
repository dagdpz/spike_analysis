function pie_data=ph_scatter(modified_keys)
% get_pie('drive','X','monkey_folder','Curius_ephys_analysis','dataset','memory')

% default_list = {                                                                              ...
%     'folder_to_save',                   'W:\Projects\Pulv_eye_hand\ephys\'                  , ...
%     'subfolder_prefix',                'cell_counts'                                         , ...
%     'monkey',                           'Curius'                                            , ...
%     'datasets',                         {'Msac'}                                            , ...
%     'case',                             'opt'                                               , ...
%     'target',                           'dPulv_l'                                           , ...
%     'instructed_choice',                'in'                                                , ...
%     'only_both_hands',                  0                                                   , ...
%     'Selection',                        {}                                                  , ...
%     'table_full_path',                  'W:\Projects\Pulv_eye_hand\ephys\ephys_analysis_v7_June2016\monkeys_Combined_del.mat'                                        , ...
%     'epochs',                           {''}                                                , ...
%     'space_criterion',                  'none'                                              , ...
%     'epoch_criterion',                  'none'                                              , ...
%     'hands_criterion',                  'none'                                              , ...
%     'SXH_criterion',                    'none'                                              , ...
%     'summary'                           'per_epoch'                                         , ...
%     'factors'                           'space'                                             , ...
%     'percent',                          0                                                   , ...
%     'create_pdf',                       0                                                   , ...
%     'all_colors',                       {[255 0 178; 171 0 252; 0 0 255; 125 130 255; 255 153 20; 222 220 0; 0 255 0; 145 143 56; 255 0 0; 255 255 255]/255}, ...
%     'contra',                           'left'                                              , ...
%     'plot_as_pie',                      0                                                   , ...
%     };


keys.monkey             = 'Curius';
keys.case               = 'opt';
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
tuning_per_unit_table(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),tuning_per_unit_table))={0};
%tuning_per_unit_table(cellfun(@(x) isempty(x) & ~islogical(x),tuning_per_unit_table))={''};

%% sort of obsolete
%keys.xlsx_table=tuning_per_unit_table;

%% ANOVA cirterions readout
criterions={'space_or_interaction','epoch_or_interaction','hands_or_interaction','SXH_or_interaction'};
for c=1:numel(criterions)
    parameter_criterion_columns=~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),criterions{c}));
    for t=1:numel(keys.tt.tasktypes)
        task_criterion_columns(t,:)=~cellfun(@isempty,strfind(tuning_per_unit_table(1,:),keys.tt.tasktypes{t}));
    end
    criterion_columns=parameter_criterion_columns & any(task_criterion_columns,1);
    keys.(criterions{c})=any(cell2mat(tuning_per_unit_table(2:end,criterion_columns)),2);
end

%% temporary, scatter stuff is to be found here:


keys.path_to_save =[keys.folder_to_save filesep keys.subfolder_prefix filesep];
if ~exist(keys.path_to_save,'dir');
    mkdir(keys.folder_to_save,keys.subfolder_prefix);
end

col_unit_ID=find_column_index(tuning_per_unit_table,'unit_ID');
monkeys=cellfun(@(x) x(1:3),tuning_per_unit_table(2:end,col_unit_ID),'uniformoutput',false);

unique_monkeys=unique(monkeys);

%% figure 1
% 
% x_label=[keys.SX.instructed_choice '_' keys.SX.epoch '_' keys.SX.parameter '_ES_' keys.SX.dataset '_' keys.SX.case ];
% y_label=[keys.SY.instructed_choice '_' keys.SY.epoch '_' keys.SY.parameter '_ES_' keys.SY.dataset '_' keys.SY.case ];
% 
% col_scatterx   =find_column_index(tuning_per_unit_table,[keys.SX.instructed_choice '_' keys.SX.epoch '_' keys.SX.parameter '_ES_' keys.SX.dataset '_' keys.SX.case ] );
% col_scattery   =find_column_index(tuning_per_unit_table,[keys.SY.instructed_choice '_' keys.SY.epoch '_' keys.SY.parameter '_ES_' keys.SY.dataset '_' keys.SY.case ] );
   

x_label=[keys.SX.instructed_choice '_' keys.SX.epoch '_' keys.SX.parameter '_ES_' keys.SX.dataset '_' keys.SX.case ];
y_label=[keys.SY.instructed_choice '_' keys.SY.epoch '_' keys.SY.parameter '_ES_' keys.SY.dataset '_' keys.SY.case ];

col_scatterx   =find_column_index(tuning_per_unit_table,[keys.SX.instructed_choice '_' keys.SX.epoch '_' keys.SX.parameter '_ES_' keys.SX.dataset '_' keys.SX.case ] );
col_scattery   =find_column_index(tuning_per_unit_table,[keys.SY.instructed_choice '_' keys.SY.epoch '_' keys.SY.parameter '_ES_' keys.SY.dataset '_' keys.SY.case ] );


col_signifix   =find_column_index(tuning_per_unit_table,[keys.SX.instructed_choice '_' keys.SX.epoch '_' keys.SX.parameter '_' keys.SX.dataset '_' keys.SX.case ] );
col_signifiy   =find_column_index(tuning_per_unit_table,[keys.SY.instructed_choice '_' keys.SY.epoch '_' keys.SY.parameter '_' keys.SY.dataset '_' keys.SY.case ] );

row_signifix   =[false; cellfun(@(x) ~strcmp(x,'-'),tuning_per_unit_table(2:end,col_signifix))];
row_signifiy   =[false; cellfun(@(x) ~strcmp(x,'-'),tuning_per_unit_table(2:end,col_signifiy))];

fig_title=[keys.selection_title{:} '_' y_label '__vs__' x_label];
plot_1_title            = [fig_title  ' FR'];
FR_index_handle = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_1_title);
%nonnanidx=~isnan([tuning_per_unit_table{2:end,col_scatterx}])'& ~isnan([tuning_per_unit_table{2:end,col_scattery}])';
nonemptyidx=~cellfun(@isempty,tuning_per_unit_table(2:end,col_scatterx)) & ~cellfun(@isempty,tuning_per_unit_table(2:end,col_scattery));

x_lim(1)=min([tuning_per_unit_table{[false;nonemptyidx],col_scatterx},tuning_per_unit_table{[false;nonemptyidx],col_scattery}]);
x_lim(2)=max([tuning_per_unit_table{[false;nonemptyidx],col_scatterx},tuning_per_unit_table{[false;nonemptyidx],col_scattery}]);
y_lim=x_lim;
hold on
line(x_lim,[0 0],'color','k','linestyle','-');
line([0 0],y_lim,'color','k','linestyle','-');
set(gca,'ylim',y_lim,'xlim',x_lim);
for m=1:numel(unique_monkeys)
    
row_monkey=ismember(monkeys,unique_monkeys{m})& nonemptyidx;
row_monkey_all=ismember(monkeys,unique_monkeys{m});
row_signifixy       = row_signifix &  row_signifiy &[false;row_monkey];
row_signifixn       = row_signifix & ~row_signifiy&[false;row_monkey];
row_signifiny       =~row_signifix &  row_signifiy&[false;row_monkey];
row_signifixn_all   = row_signifix & [false;row_monkey_all];
row_signifiny_all   = row_signifiy & [false;row_monkey_all];
row_signifinn       =~row_signifix & ~row_signifiy&[false;row_monkey];
row_signifany       = (row_signifix |  row_signifiy)&[false;row_monkey];

current_color=keys.monkey_colors(m,:);
scatter([tuning_per_unit_table{row_signifixy,col_scatterx}],[tuning_per_unit_table{row_signifixy,col_scattery}],100,current_color,'o','filled');
scatter([tuning_per_unit_table{row_signifixn,col_scatterx}],[tuning_per_unit_table{row_signifixn,col_scattery}],100,current_color,'>','filled');
scatter([tuning_per_unit_table{row_signifiny,col_scatterx}],[tuning_per_unit_table{row_signifiny,col_scattery}],100,current_color,'^','filled');
scatter([tuning_per_unit_table{row_signifinn,col_scatterx}],[tuning_per_unit_table{row_signifinn,col_scattery}],100,current_color,'o');

yvalues=[tuning_per_unit_table{[false;row_monkey],col_scattery}]';
xvalues=[tuning_per_unit_table{[false;row_monkey],col_scatterx}]';
yvalues_all=[tuning_per_unit_table{[false;row_monkey_all],col_scattery}]';
xvalues_all=[tuning_per_unit_table{[false;row_monkey_all],col_scatterx}]';
yvalues_sig=[tuning_per_unit_table{row_signifiy&[false;row_monkey],col_scattery}]';
xvalues_sig=[tuning_per_unit_table{row_signifix&[false;row_monkey],col_scatterx}]';
yvalues_sig_all=[tuning_per_unit_table{row_signifixn_all,col_scattery}]';
xvalues_sig_all=[tuning_per_unit_table{row_signifiny_all,col_scatterx}]';

[~,p_yvalues]=ttest(yvalues);
[~,p_xvalues]=ttest(xvalues);
[~,p_yvalues_all]=ttest(yvalues_all);
[~,p_xvalues_all]=ttest(xvalues_all);
[~,p_yvalues_sig]=ttest(yvalues_sig);
[~,p_xvalues_sig]=ttest(xvalues_sig);
[~,p_yvalues_sig_all]=ttest(yvalues_sig_all);
[~,p_xvalues_sig_all]=ttest(xvalues_sig_all);

Psf=polyfit([tuning_per_unit_table{row_signifany,col_scatterx}],[tuning_per_unit_table{row_signifany,col_scattery}],1);
sline=refline(Psf(1),Psf(2));set(sline,'Color',current_color,'LineWidth',2);
[Rs,Ps]=corr([tuning_per_unit_table{row_signifany,col_scatterx}]',[tuning_per_unit_table{row_signifany,col_scattery}]');
text(x_lim(1)+diff(x_lim)/20,y_lim(2)-(6*m-0)*diff(y_lim)/20,['SIG: R= ' num2str(Rs) ', P= ' num2str(Ps)],'color',current_color);

Pnf=polyfit(xvalues,yvalues,1);
sline=refline(Pnf(1),Pnf(2));set(sline,'Color',current_color,'linestyle',':','LineWidth',2);
[Ra,Pa]=corr(xvalues,yvalues);
text(x_lim(1)+diff(x_lim)/20,y_lim(2)-(6*m-1)*diff(y_lim)/20,['ALL: R= ' num2str(Ra) ', P= ' num2str(Pa)],'color',current_color);



meansstds           =sprintf('plotted: mean^ %.3f ,std^ %.3f, p^: %.3f, mean> %.3f ,std> %.3f, p>: %.3f',nanmean(yvalues),nanstd(yvalues),p_yvalues,nanmean(xvalues),nanstd(xvalues),p_xvalues);
meansstds_all       =sprintf('ALL: mean^ %.3f ,std^ %.3f, p^: %.3f, mean> %.3f ,std> %.3f, p>: %.3f',nanmean(yvalues_all),nanstd(yvalues_all),p_yvalues_all,nanmean(xvalues_all),nanstd(xvalues_all),p_xvalues_all);
meansstds_sig       =sprintf('sig: mean^ %.3f ,std^ %.3f, p^: %.3f, mean> %.3f ,std> %.3f, p>: %.3f',nanmean(yvalues_sig),nanstd(yvalues_sig),p_yvalues_sig,nanmean(xvalues_sig),nanstd(xvalues_sig),p_xvalues_sig);
meansstds_sig_all   =sprintf('SIG_ALL: mean^ %.3f ,std^ %.3f, p^: %.3f, mean> %.3f ,std> %.3f, p>: %.3f',nanmean(yvalues_sig_all),nanstd(yvalues_sig_all),p_yvalues_sig_all,nanmean(xvalues_sig_all),nanstd(xvalues_sig_all),p_xvalues_sig_all);

text(x_lim(1)+diff(x_lim)/20,y_lim(2)-(6*m-4)*diff(y_lim)/20,meansstds,'color',current_color,'interpreter','none');
text(x_lim(1)+diff(x_lim)/20,y_lim(2)-(6*m-5)*diff(y_lim)/20,meansstds_all,'color',current_color,'interpreter','none');
text(x_lim(1)+diff(x_lim)/20,y_lim(2)-(6*m-2)*diff(y_lim)/20,meansstds_sig,'color',current_color,'interpreter','none');
text(x_lim(1)+diff(x_lim)/20,y_lim(2)-(6*m-3)*diff(y_lim)/20,meansstds_sig_all,'color',current_color,'interpreter','none');
end
xlabel(x_label,'interpreter','none');
ylabel(y_label,'interpreter','none');
axis square;

title_and_save(FR_index_handle,plot_1_title,plot_1_title,keys);

%% figure 2

x_label=[keys.SX.instructed_choice '_' keys.SX.epoch '_' keys.SX.parameter '_IX_' keys.SX.dataset '_' keys.SX.case ];
y_label=[keys.SY.instructed_choice '_' keys.SY.epoch '_' keys.SY.parameter '_IX_' keys.SY.dataset '_' keys.SY.case ];

col_scatterx   =find_column_index(tuning_per_unit_table,[keys.SX.instructed_choice '_' keys.SX.epoch '_' keys.SX.parameter '_IX_' keys.SX.dataset '_' keys.SX.case ] );
col_scattery   =find_column_index(tuning_per_unit_table,[keys.SY.instructed_choice '_' keys.SY.epoch '_' keys.SY.parameter '_IX_' keys.SY.dataset '_' keys.SY.case ] );
   
% col_signifix   =find_column_index(tuning_per_unit_table,[keys.SX.instructed_choice '_' keys.SX.epoch '_' keys.SX.parameter '_' keys.SX.dataset '_' keys.SX.case ] );
% col_signifiy   =find_column_index(tuning_per_unit_table,[keys.SY.instructed_choice '_' keys.SY.epoch '_' keys.SY.parameter '_' keys.SY.dataset '_' keys.SY.case ] );
% 
% row_signifix   =[false; cellfun(@(x) ~strcmp(x,'-'),tuning_per_unit_table(2:end,col_signifix))];
% row_signifiy   =[false; cellfun(@(x) ~strcmp(x,'-'),tuning_per_unit_table(2:end,col_signifiy))];

plot_1_title            = [fig_title  ' FR_index'];
FR_index_handle = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_1_title);
x_lim=[-1,1];
y_lim=[-1,1];
hold on
line(x_lim,[0 0],'color','k','linestyle','-');
line([0 0],y_lim,'color','k','linestyle','-');
set(gca,'ylim',y_lim,'xlim',x_lim);
nonemptyidx=~cellfun(@isempty,tuning_per_unit_table(2:end,col_scatterx)) & ~cellfun(@isempty,tuning_per_unit_table(2:end,col_scattery));  
for m=1:numel(unique_monkeys)
  
row_monkey=ismember(monkeys,unique_monkeys{m})& nonemptyidx;
row_monkey_all=ismember(monkeys,unique_monkeys{m});
row_signifixy  = row_signifix &  row_signifiy &[false;row_monkey];
row_signifixn  = row_signifix & ~row_signifiy&[false;row_monkey];
row_signifixn_all   = row_signifix & [false;row_monkey_all];
row_signifiny_all   = row_signifiy & [false;row_monkey_all];
row_signifiny  =~row_signifix &  row_signifiy&[false;row_monkey];
row_signifinn  =~row_signifix & ~row_signifiy&[false;row_monkey];
row_signifany  = (row_signifix |  row_signifiy)&[false;row_monkey];


yvalues=[tuning_per_unit_table{[false;row_monkey],col_scattery}]';
xvalues=[tuning_per_unit_table{[false;row_monkey],col_scatterx}]';
yvalues_all=[tuning_per_unit_table{[false;row_monkey_all],col_scattery}]';
xvalues_all=[tuning_per_unit_table{[false;row_monkey_all],col_scatterx}]';
yvalues_sig=[tuning_per_unit_table{row_signifiy&[false;row_monkey],col_scattery}]';
xvalues_sig=[tuning_per_unit_table{row_signifix&[false;row_monkey],col_scatterx}]';
yvalues_sig_all=[tuning_per_unit_table{row_signifixn_all,col_scattery}]';
xvalues_sig_all=[tuning_per_unit_table{row_signifiny_all,col_scatterx}]';

[~,p_yvalues]=ttest(yvalues);
[~,p_xvalues]=ttest(xvalues);
[~,p_yvalues_all]=ttest(yvalues_all);
[~,p_xvalues_all]=ttest(xvalues_all);
[~,p_yvalues_sig]=ttest(yvalues_sig);
[~,p_xvalues_sig]=ttest(xvalues_sig);
[~,p_yvalues_sig_all]=ttest(yvalues_sig_all);
[~,p_xvalues_sig_all]=ttest(xvalues_sig_all);

current_color=keys.monkey_colors(m,:);
scatter([tuning_per_unit_table{row_signifixy,col_scatterx}],[tuning_per_unit_table{row_signifixy,col_scattery}],100,current_color,'o','filled');
scatter([tuning_per_unit_table{row_signifixn,col_scatterx}],[tuning_per_unit_table{row_signifixn,col_scattery}],100,current_color,'>','filled');
scatter([tuning_per_unit_table{row_signifiny,col_scatterx}],[tuning_per_unit_table{row_signifiny,col_scattery}],100,current_color,'^','filled');
scatter([tuning_per_unit_table{row_signifinn,col_scatterx}],[tuning_per_unit_table{row_signifinn,col_scattery}],100,current_color,'o');

Psf=polyfit([tuning_per_unit_table{row_signifany,col_scatterx}],[tuning_per_unit_table{row_signifany,col_scattery}],1);
sline=refline(Psf(1),Psf(2));set(sline,'Color',current_color,'LineWidth',2);
[Rs,Ps]=corr([tuning_per_unit_table{row_signifany,col_scatterx}]',[tuning_per_unit_table{row_signifany,col_scattery}]');
text(x_lim(1)+diff(x_lim)/20,y_lim(2)-(6*m)*diff(y_lim)/20,['SIG: R= ' num2str(Rs) ', P= ' num2str(Ps)],'color',current_color);

Pnf=polyfit(xvalues,yvalues,1);
sline=refline(Pnf(1),Pnf(2));set(sline,'Color',current_color,'linestyle',':','LineWidth',2);
[Ra,Pa]=corr(xvalues,yvalues);
text(x_lim(1)+diff(x_lim)/20,y_lim(2)-(6*m-1)*diff(y_lim)/20,['ALL: R= ' num2str(Ra) ', P= ' num2str(Pa)],'color',current_color);

meansstds           =sprintf('plotted: mean^ %.3f ,std^ %.3f, p^: %.3f, mean> %.3f ,std> %.3f, p>: %.3f',nanmean(yvalues),nanstd(yvalues),p_yvalues,nanmean(xvalues),nanstd(xvalues),p_xvalues);
meansstds_all       =sprintf('ALL: mean^ %.3f ,std^ %.3f, p^: %.3f, mean> %.3f ,std> %.3f, p>: %.3f',nanmean(yvalues_all),nanstd(yvalues_all),p_yvalues_all,nanmean(xvalues_all),nanstd(xvalues_all),p_xvalues_all);
meansstds_sig       =sprintf('sig: mean^ %.3f ,std^ %.3f, p^: %.3f, mean> %.3f ,std> %.3f, p>: %.3f',nanmean(yvalues_sig),nanstd(yvalues_sig),p_yvalues_sig,nanmean(xvalues_sig),nanstd(xvalues_sig),p_xvalues_sig);
meansstds_sig_all   =sprintf('SIG_ALL: mean^ %.3f ,std^ %.3f, p^: %.3f, mean> %.3f ,std> %.3f, p>: %.3f',nanmean(yvalues_sig_all),nanstd(yvalues_sig_all),p_yvalues_sig_all,nanmean(xvalues_sig_all),nanstd(xvalues_sig_all),p_xvalues_sig_all);

text(x_lim(1)+diff(x_lim)/20,y_lim(2)-(6*m-4)*diff(y_lim)/20,meansstds,'color',current_color,'interpreter','none');
text(x_lim(1)+diff(x_lim)/20,y_lim(2)-(6*m-5)*diff(y_lim)/20,meansstds_all,'color',current_color,'interpreter','none');
text(x_lim(1)+diff(x_lim)/20,y_lim(2)-(6*m-2)*diff(y_lim)/20,meansstds_sig,'color',current_color,'interpreter','none');
text(x_lim(1)+diff(x_lim)/20,y_lim(2)-(6*m-3)*diff(y_lim)/20,meansstds_sig_all,'color',current_color,'interpreter','none');
end




row_monkey=nonemptyidx;
row_signifixn_all   = row_signifix;
row_signifiny_all   = row_signifiy;
row_signifany  = (row_signifix |  row_signifiy)&[false;row_monkey];


yvalues=[tuning_per_unit_table{[false;row_monkey],col_scattery}]';
xvalues=[tuning_per_unit_table{[false;row_monkey],col_scatterx}]';
yvalues_all=[tuning_per_unit_table{2:end,col_scattery}]';
xvalues_all=[tuning_per_unit_table{2:end,col_scatterx}]';
yvalues_sig=[tuning_per_unit_table{ row_signifiy&[false;row_monkey],col_scattery}]';
xvalues_sig=[tuning_per_unit_table{ row_signifix&[false;row_monkey],col_scatterx}]';
yvalues_sig_all=[tuning_per_unit_table{row_signifixn_all,col_scattery}]';
xvalues_sig_all=[tuning_per_unit_table{row_signifiny_all,col_scatterx}]';

[~,p_yvalues]=ttest(yvalues);
[~,p_xvalues]=ttest(xvalues);
[~,p_yvalues_all]=ttest(yvalues_all);
[~,p_xvalues_all]=ttest(xvalues_all);
[~,p_yvalues_sig]=ttest(yvalues_sig);
[~,p_xvalues_sig]=ttest(xvalues_sig);
[~,p_yvalues_sig_all]=ttest(yvalues_sig_all);
[~,p_xvalues_sig_all]=ttest(xvalues_sig_all);
m=3;
[Ra,Pa]=corr(xvalues,yvalues);
[Rs,Ps]=corr([tuning_per_unit_table{row_signifany,col_scatterx}]',[tuning_per_unit_table{row_signifany,col_scattery}]');
text(x_lim(1)+diff(x_lim)/20,y_lim(2)-(6*m)*diff(y_lim)/20,['SIG: R= ' num2str(Rs) ', P= ' num2str(Ps)],'color','k');
text(x_lim(1)+diff(x_lim)/20,y_lim(2)-(6*m-1)*diff(y_lim)/20,['ALL: R= ' num2str(Ra) ', P= ' num2str(Pa)],'color','k');

meansstds           =sprintf('plotted: mean^ %.3f ,std^ %.3f, p^: %.3f, mean> %.3f ,std> %.3f, p>: %.3f',nanmean(yvalues),nanstd(yvalues),p_yvalues,nanmean(xvalues),nanstd(xvalues),p_xvalues);
meansstds_all       =sprintf('ALL: mean^ %.3f ,std^ %.3f, p^: %.3f, mean> %.3f ,std> %.3f, p>: %.3f',nanmean(yvalues_all),nanstd(yvalues_all),p_yvalues_all,nanmean(xvalues_all),nanstd(xvalues_all),p_xvalues_all);
meansstds_sig       =sprintf('sig: mean^ %.3f ,std^ %.3f, p^: %.3f, mean> %.3f ,std> %.3f, p>: %.3f',nanmean(yvalues_sig),nanstd(yvalues_sig),p_yvalues_sig,nanmean(xvalues_sig),nanstd(xvalues_sig),p_xvalues_sig);
meansstds_sig_all   =sprintf('SIG_ALL: mean^ %.3f ,std^ %.3f, p^: %.3f, mean> %.3f ,std> %.3f, p>: %.3f',nanmean(yvalues_sig_all),nanstd(yvalues_sig_all),p_yvalues_sig_all,nanmean(xvalues_sig_all),nanstd(xvalues_sig_all),p_xvalues_sig_all);

text(x_lim(1)+diff(x_lim)/20,y_lim(2)-(6*m-4)*diff(y_lim)/20,meansstds,'color','k','interpreter','none');
text(x_lim(1)+diff(x_lim)/20,y_lim(2)-(6*m-5)*diff(y_lim)/20,meansstds_all,'color','k','interpreter','none');
text(x_lim(1)+diff(x_lim)/20,y_lim(2)-(6*m-2)*diff(y_lim)/20,meansstds_sig,'color','k','interpreter','none');
text(x_lim(1)+diff(x_lim)/20,y_lim(2)-(6*m-3)*diff(y_lim)/20,meansstds_sig_all,'color','k','interpreter','none');

xlabel(x_label,'interpreter','none');
ylabel(y_label,'interpreter','none');
axis square;

title_and_save(FR_index_handle,plot_1_title,plot_1_title,keys);



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


% function xlsx_table = get_mat_table(keys)
% xlsx_table={};
% load(keys.table_full_path);
% tuning_per_unit_table=LR_to_CI(tuning_per_unit_table);
% tuning_per_unit_table=tuning_per_unit_table;
% tuning_per_unit_table=tuning_per_unit_table;
% tuning_per_unit_table(cellfun(@isempty,tuning_per_unit_table))={''};
% tuning_per_unit_table(cellfun(@isempty,tuning_per_unit_table))={NaN};
% 

