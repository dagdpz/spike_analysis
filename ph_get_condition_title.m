function condition_title=ph_get_condition_title(keys)
condition_title='';
for c=1:numel(keys.condition_parameters)
    condition_title=[condition_title sprintf([keys.condition_parameters{c} ' %s'],mat2str(keys.tt.(keys.condition_parameters{c})))];
end
end