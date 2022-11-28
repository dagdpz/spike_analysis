function [fig_title,filename]=ph_figure_titles(keys)

PF=keys.normalization_field;
Sel_for_title=keys.selection_title;

    if strcmp(keys.(PF).normalization,'percent_change')
        normalization_label=sprintf('N_prct %s to %s',keys.(PF).epoch_BL,keys.(PF).epoch_DN);
    elseif ~keys.(PF).FR_subtract_baseline
        normalization_label=sprintf('N_%s in %s',keys.(PF).normalization,keys.(PF).epoch_DN);
    elseif keys.(PF).baseline_per_trial
        normalization_label=sprintf('N_%s in %s, - %s per trial',keys.(PF).normalization,keys.(PF).epoch_DN,keys.(PF).epoch_BL);
    else
        normalization_label=sprintf('N_%s (X-%s)by%s',keys.(PF).normalization,keys.(PF).epoch_BL,keys.(PF).epoch_DN);
    end
    condition_title=ph_get_condition_title(keys);
    
    if isempty(keys.(PF).unique_title)
    fig_title=sprintf('%s %s %s %s %s %s grouped by %s ',...
        keys.monkey,[keys.conditions_to_plot{:}],keys.arrangement,condition_title,[Sel_for_title{:}],normalization_label,keys.(PF).group_parameter);
    filename=sprintf('%s %s %s %s %s %s %s ',...
        keys.monkey,[keys.conditions_to_plot{:}],keys.arrangement(1:3),condition_title,[Sel_for_title{:}],normalization_label,keys.(PF).group_parameter);
    else
    fig_title=sprintf('%s %s %s %s %s %s %s grouped by %s ',...
        keys.(PF).unique_title, keys.monkey,[keys.conditions_to_plot{:}],keys.arrangement,condition_title,[Sel_for_title{:}],normalization_label,keys.(PF).group_parameter);
    filename=sprintf('%s %s %s %s %s %s %s %s ',...
        keys.(PF).unique_title, keys.monkey,[keys.conditions_to_plot{:}],keys.arrangement(1:3),condition_title,[Sel_for_title{:}],normalization_label,keys.(PF).group_parameter);
    end
end

function condition_title=ph_get_condition_title(keys)
condition_title='';
for c=1:numel(keys.condition_parameters)
    condition_title=[condition_title sprintf([keys.condition_parameters{c} ' %s'],mat2str(keys.tt.(keys.condition_parameters{c})))];
end
end