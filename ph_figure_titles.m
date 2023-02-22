function [fig_title,filename,condition_title]=ph_figure_titles(keys)

PF=keys.normalization_field;
Sel_for_title=keys.selection_title;

if strcmp(keys.(PF).normalization,'percent_change')
    normalization_label=sprintf('N_prct %s to %s',keys.(PF).epoch_BL,keys.(PF).epoch_DN);
elseif strcmp(keys.(PF).epoch_BL,'none')
    if strcmp(keys.(PF).epoch_DN,'none') || strcmp(keys.(PF).normalization,'none')
        normalization_label=sprintf('N_none');
    else
        normalization_label=sprintf('N_%s in %s',keys.(PF).normalization,keys.(PF).epoch_DN);
    end
elseif keys.(PF).baseline_per_trial
    normalization_label=sprintf('N_%s in %s, - %s per trial',keys.(PF).normalization,keys.(PF).epoch_DN,keys.(PF).epoch_BL);
elseif strcmp(keys.(PF).epoch_DN,'none') || strcmp(keys.(PF).normalization,'none')
    normalization_label=sprintf('N_%s (X-%s)',keys.(PF).normalization,keys.(PF).epoch_BL);
else
    normalization_label=sprintf('N_%s (X-%s)by%s',keys.(PF).normalization,keys.(PF).epoch_BL,keys.(PF).epoch_DN);
end
[condition_title_short, condition_title_full]=ph_get_condition_title(keys);

if isempty(keys.(PF).unique_title)
    fig_title=sprintf('%s %s %s %s %s%s grouped by %s ',...
        keys.monkey,[keys.conditions_to_plot{:}],keys.arrangement,condition_title_full,[Sel_for_title{:}],normalization_label,keys.(PF).group_parameter);
    filename=sprintf('%s %s %s %s %s%s %s ',...
        keys.monkey,[keys.conditions_to_plot{:}],keys.arrangement(1:3),condition_title_short,[Sel_for_title{:}],normalization_label,keys.(PF).group_parameter);
else
    fig_title=sprintf('%s %s %s %s %s%s %s grouped by %s ',...
        keys.(PF).unique_title, keys.monkey,[keys.conditions_to_plot{:}],keys.arrangement,condition_title_full,[Sel_for_title{:}],normalization_label,keys.(PF).group_parameter);
    filename=sprintf('%s %s %s %s %s%s %s %s ',...
        keys.(PF).unique_title, keys.monkey,[keys.conditions_to_plot{:}],keys.arrangement(1:3),condition_title_short,[Sel_for_title{:}],normalization_label,keys.(PF).group_parameter);
end
end

function [condition_title_short, condition_title_full]=ph_get_condition_title(keys)
condition_title_short='';
condition_title_full='';
for c=1:numel(keys.condition_parameters)
    condition_title_short=[condition_title_short sprintf([keys.condition_parameters{c}(1:3) '=%s, '],mat2str(keys.tt.(keys.condition_parameters{c})))];
    condition_title_full=[condition_title_full   sprintf([keys.condition_parameters{c} '=%s, '],mat2str(keys.tt.(keys.condition_parameters{c})))];
end
end