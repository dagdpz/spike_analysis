function [result_table cells_per_epoch]=ph_tuning_readout(project,version)
%OO=vertcat(population.RF_per_epoch);

type=3;
method='FR_vector'; %'FR_vector'; 'Von_mises_fit'; %

Epochs_to_show= {'Fhol','Cue','MemL','PreS','PeriS','TIhol','Tons'};
Epochs_full_names= {'Fix','Cue','MemL','Pre','Peri','Post','Tons'};

keys=struct;
keys=ph_general_settings(project,keys);
project_specific_settings=[keys.db_folder filesep project filesep 'ph_project_settings.m'];
run(project_specific_settings);
keys.project_versions=version;
keys.project_version=keys.project_versions{1};

version_folder=keys.project_version;
keys.version_specific_settings=[keys.db_folder project filesep keys.project_version filesep 'ph_project_version_settings.m'];
run(keys.version_specific_settings);
keys.project_version=version_folder;

keys.drive='Y:\';
keys.basepath_to_save=['Projects' filesep project filesep 'ephys' filesep];
keys.path_to_save=[keys.drive keys.basepath_to_save keys.project_version filesep 'tuning_curves' filesep];
keys.anova_table_file=[keys.drive keys.basepath_to_save keys.project_version filesep 'tuning_table_combined_CI.mat'];

% if ~exist(keys.path_to_save,'dir');
%     mkdir(keys.basepath_to_save,'tuning_curves');
% end

load([keys.path_to_save 'population_' method]);
keys.tt.IC_for_criterion='in';
keys.monkey='';
keys.tt.tasktypes={'Msac_mov'};
keys.tt.hands=0;
keys.tt.choices=0;
keys.tt.perturbations=0;
keys.tt.selection={'target','dPul'};
[tuning_per_unit_table]=ph_load_extended_tuning_table(keys);
[tuning_per_unit_table, keys.selection_title]=ph_reduce_tuning_table(tuning_per_unit_table,keys);

unit_IDs=tuning_per_unit_table(2:end,1);
unit_IDs_pop={population.unit_ID};
valid_IDs=ismember(unit_IDs_pop,unit_IDs);
unit_IDs_pop=unit_IDs_pop(valid_IDs);
N=sum(valid_IDs);
OO=vertcat(population(valid_IDs).RF_per_epoch);



EPOCHS=keys.EPOCHS_PER_TYPE{type};

EP=ismember(EPOCHS(:,1),Epochs_to_show);
for e=2:size(EPOCHS,1)
% N_1_dim(e)=sum([OO(:,e).d_tuning_vector_shift]==1);
% N_2_dim(e)=sum([OO(:,e).d_tuning_vector_shift]==2);
% for u=1:size(OO,1)
% size_p(e,u).x=size(OO(u,e).p_tuning_vector_per_line);
% end
unit_valid=[(OO(:,e).unit_valid)];
gain_signal_change=[(OO(:,e).gain_signal_change)];
reasonable_tuning=true(size(unit_valid')); %% placeholder
%reasonable_tuning=unit_valid'; 
%ANOVA_position_effect=arrayfun(@(x) any(x.p_anova*3<0.05),OO(:,e));
significant_gain=vertcat(OO(:,e).significant_gain);
significant_shift=vertcat(OO(:,e).significant_shift);
significant_RF=vertcat(OO(:,e).significant_RF);


%% careful here!! ANOVA_position_effect

positionanovaindex=DAG_find_column_index(tuning_per_unit_table,['in_AH_' EPOCHS{e,1} '_position_Msac_mov']);
gazeanovaindex=DAG_find_column_index(tuning_per_unit_table,['in_AH_' EPOCHS{e,1} '_fixation_Msac_mov']);
interactionanovaindex=DAG_find_column_index(tuning_per_unit_table,['in_AH_' EPOCHS{e,1} '_PxF_Msac_mov']);
unit_IDs_over_all_position_effect=tuning_per_unit_table(strcmp(['false'; tuning_per_unit_table(2:end,positionanovaindex)],'true'),1);
unit_IDs_gazedependence=tuning_per_unit_table(strcmp(['false'; tuning_per_unit_table(2:end,gazeanovaindex)],'true') | strcmp(['false'; tuning_per_unit_table(2:end,interactionanovaindex)],'true'),1);
%unit_IDs_gazedependence=tuning_per_unit_table(strcmp(['false'; tuning_per_unit_table(2:end,interactionanovaindex)],'true'),1);

ANOVA_position_effect=ismember(unit_IDs_pop,unit_IDs_over_all_position_effect)';
gaze_dependence=ismember(unit_IDs_pop,unit_IDs_gazedependence)';
gaze_dependence=true(size(unit_valid')); %% placeholder
%ANOVA_position_effect=reasonable_tuning;
% reasonable_tuning=ismember(unit_IDs_pop,unit_IDs_gazedependence)';
% gaze_dependence=true(size(unit_valid')); %% placeholder

only_gain=       significant_gain & ~significant_shift & ANOVA_position_effect &  reasonable_tuning;
only_shift=     ~significant_gain &  significant_shift & ANOVA_position_effect &  reasonable_tuning;
gain_and_shift=  significant_gain &  significant_shift & ANOVA_position_effect &  reasonable_tuning;
only_tuning=    ~significant_gain & ~significant_shift & ANOVA_position_effect &  reasonable_tuning;
ANOVA_true=     ~significant_gain & ~significant_shift & ANOVA_position_effect & ~reasonable_tuning;

gainfield(e)=sum(only_gain);
RF_shift(e)=sum(only_shift);
Gain_and_shift(e)=sum(gain_and_shift);
Tuning_vector_with_gaze(e)=sum(only_tuning & gaze_dependence);
Tuning_vector_wo_gaze(e)=sum(only_tuning & ~gaze_dependence);
ANOVA_with_gaze(e)=sum(ANOVA_true & gaze_dependence);
ANOVA_wo_gaze(e)=sum(ANOVA_true & ~gaze_dependence);

cells_per_epoch(e).epoch=EPOCHS{e,1};
cells_per_epoch(e).gain=unit_IDs_pop(only_gain);
cells_per_epoch(e).RF_shift=unit_IDs_pop(only_shift);
cells_per_epoch(e).Gain_and_shift=unit_IDs_pop(gain_and_shift);
cells_per_epoch(e).ANOVA_only=unit_IDs_pop(significant_gain~=1 & significant_shift~=1 & ANOVA_position_effect & reasonable_tuning);
% 
% N_response_fields(e)=sum(significant_RF==1);
% gainfields_with_receptive_field(e)=sum(significant_RF==1 & significant_gain==1);
% shifts_with_receptive_field(e)=sum(significant_RF==1 & significant_shift==1);
% gainfields_with_receptive_field_any(e)=sum(ANOVA_position_effect==1 & significant_gain==1);
% shifts_with_receptive_field_any(e)=sum(ANOVA_position_effect==1 & significant_shift==1);
% valid_units(e)=sum(unit_valid);

mean_gain_effect(e)=nanmean(gain_signal_change(only_gain | gain_and_shift));
std_gain_effect(e)=nanstd(gain_signal_change(only_gain | gain_and_shift));
end
result_table=[...
['Epoch';EPOCHS(:,1)],... 
['Gain';num2cell(gainfield')],...
['Gain and Shift';num2cell(Gain_and_shift')],...
['Shift';num2cell(RF_shift')],...
['tuning_vector_wo_gaze';num2cell(Tuning_vector_wo_gaze')],...
['ANOVA_true_wo_gaze';num2cell(ANOVA_wo_gaze')],...
['gain_signal_change_mean';num2cell(mean_gain_effect')],...
['gain_signal_change_std';num2cell(std_gain_effect')],...
];

save([keys.path_to_save filesep 'RF_gain_shift_table'],'result_table'); 
fig=figure('units','normalized','outerposition',[0 0 1 1],'name','Gain field summary');
colormap([1 0 0; 0 1 0; 0 0 1; 0.3 0.3 0.3; 0.5 0.5 0.5; 0.7 0.7 0.7; 1 1 1]);
bar(1:sum(EP),[gainfield(EP)' Gain_and_shift(EP)' RF_shift(EP)' Tuning_vector_with_gaze(EP)' Tuning_vector_wo_gaze(EP)' ANOVA_with_gaze(EP)' ANOVA_wo_gaze(EP)']/N*100,'stacked');

topplot={gainfield,Gain_and_shift,RF_shift,Tuning_vector_with_gaze,Tuning_vector_wo_gaze,ANOVA_with_gaze,ANOVA_wo_gaze};
posy=0;
for k=1:numel(topplot)
    curr=topplot{k};
    text(1:sum(EP),posy+curr(EP)/N*100/2,strcat(cellstr(num2str(curr(EP)')), ' (', cellstr(num2str(round(curr(EP)'/N*100))), '%)'),'HorizontalAlignment','center');
    posy=posy+curr(EP)/N*100;
end
% text(1:sum(EP),gainfield(EP)/N*100 + Gain_and_shift(EP)/N*100/2,strcat(cellstr(num2str(Gain_and_shift(EP)')), ' (', cellstr(num2str(round(Gain_and_shift(EP)'/N*100))), '%)'),'HorizontalAlignment','center')
% text(1:sum(EP),gainfield(EP)/N*100 + Gain_and_shift(EP)/N*100 + RF_shift(EP)/N*100/2,strcat(cellstr(num2str(RF_shift(EP)')), ' (', cellstr(num2str(round(RF_shift(EP)'/N*100))), '%)'),'HorizontalAlignment','center')
% text(1:sum(EP),gainfield(EP)/N*100 + Gain_and_shift(EP)/N*100 + RF_shift(EP)/N*100 + N_smooth_tuning_only(EP)/N*100/2,...
% strcat(cellstr(num2str(N_smooth_tuning_only(EP)')), ' (', cellstr(num2str(round(N_smooth_tuning_only(EP)'/N*100))), '%)'),'HorizontalAlignment','center')
% text(1:sum(EP),gainfield(EP)/N*100 + Gain_and_shift(EP)/N*100 + RF_shift(EP)/N*100 + N_smooth_tuning_only(EP)/N*100 + Only_response_field(EP)/N*100/2,...
% strcat(cellstr(num2str(Only_response_field(EP)')), ' (', cellstr(num2str(round(Only_response_field(EP)'/N*100))), '%)'),'HorizontalAlignment','center')

legend('Gain','Gain and shift', 'Shift', 'tuning vector + gaze dep','tuning vector ', 'retinocentric + gaze dep', 'retinocentric');
ylabel(['Number of cells [%], N=' num2str(N)])
set(gca,'xticklabel',Epochs_full_names);


wanted_size=[50 30];
set(fig, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size])
export_fig([keys.path_to_save 'Gain field summary'], '-pdf','-transparent') % pdf by run
close all
end

