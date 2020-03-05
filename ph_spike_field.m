
session_folder='Y:\Projects\Pulv_oculomotor\ephys\20180222\';
dir_session_folder=dir([session_folder '*.mat']);
names={dir_session_folder.name};
population_names=names(cellfun(@(x) any(strfind(x,'population_')),names));

% load('Y:\Projects\Pulv_oculomotor\ephys\20180222\sites_Curius_20150618.mat')
% load('Y:\Projects\Pulv_oculomotor\ephys\20180222\population_Curius_20150618.mat')
for sess=1:numel(population_names)
    clear population sites
    load([session_folder population_names{sess}]);
    load([session_folder 'sites' population_names{sess}(11:end)]);    
    
All_site_IDs={sites.site_ID};
for u=1:numel(population)
   site_ID=population(u).site_ID;
   s=find(ismember(All_site_IDs,site_ID));
   site=sites(s);
   s_tr_mat=[site.trial.block; site.trial.run; site.trial.n];
   u_tr_mat=[population(u).trial.block; population(u).trial.run; population(u).trial.n];
   
   s_idx=ismember(s_tr_mat',u_tr_mat','rows');
   u_idx=ismember(u_tr_mat',s_tr_mat','rows');
   
   sum(s_idx)==sum(u_idx)
end
end