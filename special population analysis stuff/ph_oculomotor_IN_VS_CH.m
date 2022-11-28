load('Y:\Projects\Pulv_oculomotor\ephys\paper_remastered\tuning_table_combined_CI.mat')

ix_center_x=DAG_find_column_index(tuning_per_unit_table,'NH_IN_Cue_center_x_Msac_options');
ix_size_x=DAG_find_column_index(tuning_per_unit_table,'NH_IN_Cue_size_x_Msac_options');
ix_CH_IN=DAG_find_column_index(tuning_per_unit_table,'in_AH_Cue_prefH_Msac_opt');
ix_CH_FR=DAG_find_column_index(tuning_per_unit_table,'ch_AH_Cue_prefH_FR_Msac_opt');
ix_IN_FR=DAG_find_column_index(tuning_per_unit_table,'in_AH_Cue_prefH_FR_Msac_opt');

ix_ch_trial_criterion=DAG_find_column_index(tuning_per_unit_table,'ch_AH_trials_per_hemifield_Msac_opt');

for kk=1:2

row_CH=cellfun(@(x) strcmp(x,'ch'),tuning_per_unit_table(:,ix_CH_IN));
row_IN=cellfun(@(x) strcmp(x,'in'),tuning_per_unit_table(:,ix_CH_IN));
row_NN=cellfun(@(x) strcmp(x,'-'),tuning_per_unit_table(:,ix_CH_IN));

if kk==2
  
row_CH=cellfun(@(x,y) ~isempty(x) && ~isempty(y) && isnumeric(x) && isnumeric(y) && x>y,tuning_per_unit_table(:,ix_CH_FR),tuning_per_unit_table(:,ix_IN_FR));  
row_IN=cellfun(@(x,y) ~isempty(x) && ~isempty(y) && isnumeric(x) && isnumeric(y) && x<y,tuning_per_unit_table(:,ix_CH_FR),tuning_per_unit_table(:,ix_IN_FR));  
row_NN=cellfun(@(x,y) ~isempty(x) && ~isempty(y) && isnumeric(x) && isnumeric(y) && x==y,tuning_per_unit_table(:,ix_CH_FR),tuning_per_unit_table(:,ix_IN_FR));  
end

row_size_valid=cellfun(@(x) ~isempty(x) && ~strcmp(x,'-'),tuning_per_unit_table(:,ix_center_x));
row_size_valid2=cellfun(@(x) ~isempty(x) && isnumeric(x) && x>=4,tuning_per_unit_table(:,ix_ch_trial_criterion)); %% trial crierion for choices !!
row_size_valid3=row_size_valid&row_size_valid2;


sizes.CH=2*vertcat(tuning_per_unit_table{row_CH & row_size_valid3,ix_size_x});
sizes.IN=2*vertcat(tuning_per_unit_table{row_IN & row_size_valid3,ix_size_x});
sizes.ALL=2*vertcat(tuning_per_unit_table{(row_IN | row_CH | row_NN) & row_size_valid3,ix_size_x});

count.CH=numel(sizes.CH);
count.IN=numel(sizes.IN);
count.ALL=numel(sizes.ALL);

centers.CH=vertcat(tuning_per_unit_table{row_CH & row_size_valid3,ix_center_x});
centers.IN=vertcat(tuning_per_unit_table{row_IN & row_size_valid3,ix_center_x});
centers.ALL=vertcat(tuning_per_unit_table{(row_IN | row_CH | row_NN) & row_size_valid3,ix_center_x});

[sig.siz_CH, ps.siz_CH]=ttest2(sizes.CH,sizes.IN);
[sig.cen_CH, ps.cen_CH]=ttest2(centers.CH,centers.IN);
[sig.cen_abs_CH, ps.cen_abs_CH]=ttest2(abs(centers.CH),abs(centers.IN));

means_sizes.CH=mean(sizes.CH);
means_sizes.IN=mean(sizes.IN);
means_sizes.ALL=mean(sizes.ALL);
std_sizes.CH=std(sizes.CH);
std_sizes.IN=std(sizes.IN);
std_sizes.ALL=std(sizes.ALL);
means_center.CH=mean(centers.CH);
means_center.IN=mean(centers.IN);
means_center.ALL=mean(centers.ALL);
means_center_abs.CH=mean(abs(centers.CH));
means_center_abs.IN=mean(abs(centers.IN));
means_center_abs.ALL=mean(abs(centers.ALL));

count
std_sizes
means_sizes
ps
end