function completed_table = update_mastertable_cell(original_table,table_for_updating,rows_to_update)
% This function updates all rows specified with raw_data_rows of a original_table (supposedly a cell, with titles in the first row)
% With the corresponding rows of table_for_updating (supposedly a cell, with titles in the first row). 
% dimensions of both tables don't need to fit.
% Order of columns does not matter, since this function compares the title
% names in order to update the correct column (So title names are relevant and should match)
% all rows_to_update should be larger than 1, first row should be titles!


titles_update=table_for_updating(1,:);
titles_original_table=original_table(1,:);
%all_rows=[NaN, 1:max(size(original_table(2:end,1),1),size(table_for_updating(2:end,1),1))]; % 
all_rows=1:max(size(original_table,1),size(table_for_updating,1)); % 
logidx_columns=ismember(all_rows,rows_to_update);

completed_table=original_table;
for n_update_column=1:numel(titles_update)
    titles_completed_table=completed_table(1,:);
    all_updated_columns=1:size(titles_completed_table,2);
    if ismember(titles_update{n_update_column},titles_original_table)
        logidx_completed_column=ismember(titles_completed_table,titles_update{n_update_column});
        n_completed_column=all_updated_columns(logidx_completed_column);
    else
        n_completed_column=all_updated_columns(end)+1;
    end
    completed_table(1,n_completed_column)=titles_update(n_update_column);
    if ~all(isempty([table_for_updating{logidx_columns,n_update_column}])) ||...
            ~ismember(titles_update{n_update_column},titles_original_table)
        % only overwrite if new entry is meaningful
    completed_table(logidx_columns,n_completed_column)=table_for_updating(logidx_columns,n_update_column);
    end
end
end