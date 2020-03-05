function reach_tuning_inactivation(project,versions)
%reach_tuning_inactivation('PPC_pulv_eye_hand',{'MIP_dPul_Inj_working','LIP_dPul_Inj_working'})
% project='PPC_pulv_eye_hand';
% version='MIP_dPul_Inj_working';
%keys.tt.tasktypes                    ={'Ddre_han'};


%
% keys.tt.choices             = 0;
% keys.tt.hands               = [1 2];
% keys.tt.perturbations       = 0;
% keys.cal.min_trials_per_condition=5;
% keys.tt.combine_tuning_properties   ={}; %% additional table entry from combining columns
% keys.monkey                    ='';
%
%
%
% keys.labels.handsLR={'AH','LH','RH'};
% keys.labels.handsIC={'AH','IH','CH'};  %% AH!??
% keys.labels.perturbations={'','_PT'};  %% AH!??
% keys.labels.choices={'in','ch'};
% keys.tt.epoch_criterion             ='none';
% keys.tt.space_criterion             ='none';
% keys.tt.hands_criterion             ='none';
% keys.tt.SXH_criterion               ='none';
%keys.tt.IC_for_criterion            = 'in';
% keys.tt.unselect                    ={};

keys=struct;
keys=ph_general_settings(project,keys);



angles=0:pi/180:2*pi;
x_circle=cos(angles)*0.5;
y_circle=sin(angles)*0.5;

keys.tt.IC_for_criterion            = 'in';


project_specific_settings=[keys.db_folder filesep project filesep 'ph_project_settings.m'];
run(project_specific_settings);

keys.project_versions=versions;
analyzeTuningIndices = true;

if analyzeTuningIndices
	subplot_rows = 6;
else
	subplot_rows = 4;
end
for f=1:numel(keys.project_versions)
	
	keys.project_versions=versions;
	if ~isempty(keys.project_versions{f})
		keys.project_version=keys.project_versions{f};
	end
	version_folder=keys.project_version;
	keys.version_specific_settings=[keys.db_folder project filesep keys.project_version filesep 'ph_project_version_settings.m'];
	run(keys.version_specific_settings);
	keys.project_version=version_folder;
	
	
	keys.tt.choices             = 0;
	keys.tt.hands               = [1 2];
	keys.tt.perturbations       = [0 1];
	
	for m=1:numel(keys.batching.monkeys)
		keys.monkey=keys.batching.monkeys{m};
		for t=1:numel(keys.batching.targets)
			target=keys.batching.targets{t};
			
			keys.tt.selection                   ={'target',target};
			tasktype=keys.tt.tasktypes{1};
			keys.anova_table_file=['Y:' filesep 'Projects'  filesep  project filesep 'ephys' filesep version_folder filesep 'tuning_table_combined_CI.mat'];
			[tuning_per_unit_table]                 = ph_load_extended_tuning_table(keys);
			[tuning_per_unit_table, Sel_for_title]  = ph_reduce_tuning_table(tuning_per_unit_table,keys);
			
			
			for kk=1:size(tuning_per_unit_table,2)
				idx.(tuning_per_unit_table{1,kk})=kk;
			end
			
			
			if any(strfind(keys.project_version,'MIP'))
				epochs={'Facq','Fhol','Cue','Del','PreR','PeriR'};
				% epochs={'Del'}; % for debug
			elseif  any(strfind(keys.project_version,'LIP'))
				epochs={'Facq','Fhol','Cue','Del','PreS','PeriS'};
				% epochs={'Cue','Del','PeriS'}; % for debug
			end
			
			
			for normalize_by_max_FR=0:1
				
				% variant 1: one task, two conditions (e.g. pre- and post-injection), two epochs: "baseline" and "response",
				% four trial types (contralesional and ipsilesional hand and space)
				
				%close all;
				figure;
				max_sum=0;
				for e=1:numel(epochs)
					
					epoch=epochs{e};
					%subplot(ceil(sqrt(numel(epochs)+1)),ceil(sqrt(numel(epochs)+1)),e);
					subplot(subplot_rows,numel(epochs),2*numel(epochs)+e);
					title(epoch);
					
					
					if any(strfind(target,'_L'))
						% 1 - control
						A1_isch = [tuning_per_unit_table{2:end, idx.(['in_IH_CS_' epoch '_CT_FR_' tasktype])}];
						A1_csch = [tuning_per_unit_table{2:end, idx.(['in_IH_IS_' epoch '_CT_FR_' tasktype])}];
						A1_csih = [tuning_per_unit_table{2:end, idx.(['in_CH_IS_' epoch '_CT_FR_' tasktype])}];
						A1_isih = [tuning_per_unit_table{2:end, idx.(['in_CH_CS_' epoch '_CT_FR_' tasktype])}];
						
						% 2 - inactivation
						A2_isch = [tuning_per_unit_table{2:end, idx.(['in_IH_CS_' epoch '_PT_FR_' tasktype])}];
						A2_csch = [tuning_per_unit_table{2:end, idx.(['in_IH_IS_' epoch '_PT_FR_' tasktype])}];
						A2_csih = [tuning_per_unit_table{2:end, idx.(['in_CH_IS_' epoch '_PT_FR_' tasktype])}];
						A2_isih = [tuning_per_unit_table{2:end, idx.(['in_CH_CS_' epoch '_PT_FR_' tasktype])}];
						
					else
						% 1 - control
						A1_isch = [tuning_per_unit_table{2:end, idx.(['in_CH_IS_' epoch '_CT_FR_' tasktype])}];
						A1_csch = [tuning_per_unit_table{2:end, idx.(['in_CH_CS_' epoch '_CT_FR_' tasktype])}];
						A1_csih = [tuning_per_unit_table{2:end, idx.(['in_IH_CS_' epoch '_CT_FR_' tasktype])}];
						A1_isih = [tuning_per_unit_table{2:end, idx.(['in_IH_IS_' epoch '_CT_FR_' tasktype])}];
						
						% 2 - inactivation
						A2_isch = [tuning_per_unit_table{2:end, idx.(['in_CH_IS_' epoch '_PT_FR_' tasktype])}];
						A2_csch = [tuning_per_unit_table{2:end, idx.(['in_CH_CS_' epoch '_PT_FR_' tasktype])}];
						A2_csih = [tuning_per_unit_table{2:end, idx.(['in_IH_CS_' epoch '_PT_FR_' tasktype])}];
						A2_isih = [tuning_per_unit_table{2:end, idx.(['in_IH_IS_' epoch '_PT_FR_' tasktype])}];
					end
					
					unit_IDs=tuning_per_unit_table(2:end, idx.unit_ID);
					
					
					map = hsv(360);
					
					max_per_cell = max([A1_isch; A1_csch; A1_csih; A1_isih]);
					if normalize_by_max_FR % normalize by max rate in condtion 1 in each cell
						A1_isch = A1_isch./max_per_cell;
						A1_csch = A1_csch./max_per_cell;
						A1_csih = A1_csih./max_per_cell;
						A1_isih = A1_isih./max_per_cell;
						
						A2_isch = A2_isch./max_per_cell;
						A2_csch = A2_csch./max_per_cell;
						A2_csih = A2_csih./max_per_cell;
						A2_isih = A2_isih./max_per_cell;
						
					end
					
					%cell_valid=max_per_cell>20;
					%                         A1_isch = A1_isch(cell_valid);
					%                         A1_csch = A1_csch(cell_valid);
					%                         A1_csih = A1_csih(cell_valid);
					%                         A1_isih = A1_isih(cell_valid);
					%
					%                         A2_isch = A2_isch(cell_valid);
					%                         A2_csch = A2_csch(cell_valid);
					%                         A2_csih = A2_csih(cell_valid);
					%                         A2_isih = A2_isih(cell_valid);
					
					
					%% Resorting by angle
					A_mat=[A1_isch' A1_csch' A1_csih' A1_isih'];
					A2_mat=[A2_isch' A2_csch' A2_csih' A2_isih'];
					phi=pi/4;
					R_mat=[sin(phi)+1i*cos(phi), -sin(phi)+1i*cos(phi), -sin(phi)-1i*cos(phi), sin(phi)-1i*cos(phi)];
					I1_mat=R_mat*A_mat';
					I2_mat=R_mat*A2_mat';
					angles_sum=round(180/pi*angle(I1_mat))+180;
                    
           
					%                                     [~,sort_idx]=sort(abs(I2_mat-I1_mat));
					%
					%
					% %sort_idx(end-2:end)=[];
					%
					%                                     A1_isch=A1_isch(sort_idx);
					%                                     A1_csch=A1_csch(sort_idx);
					%                                     A1_csih=A1_csih(sort_idx);
					%                                     A1_isih=A1_isih(sort_idx);
					%
					%                                     A2_isch=A2_isch(sort_idx);
					%                                     A2_csch=A2_csch(sort_idx);
					%                                     A2_csih=A2_csih(sort_idx);
					%                                     A2_isih=A2_isih(sort_idx);
					%                     sorted_unit_ID=unit_IDs(sort_idx);
					
					%                     to_exclude=ismember(sorted_unit_ID,{'Lin_20171026_10','Lin_20171109_09','Lin_20171109_08'});
					%                     A1_isch(to_exclude)=[];
					%                     A1_csch(to_exclude)=[];
					%                     A1_csih(to_exclude)=[];
					%                     A1_isih(to_exclude)=[];
					%                     A2_isch(to_exclude)=[];
					%                     A2_csch(to_exclude)=[];
					%                     A2_isih(to_exclude)=[];
					%
					%
					%                     I1_mat(to_exclude)=[];
					%                     I2_mat(to_exclude)=[];
					
					TI1_space = [];	TI2_space = []; TI1_hand = []; TI2_hand = [];
					TI1_pref_space = [];	TI2_pref_space = []; TI1_pref_hand = []; TI2_pref_hand = [];
					
					for k=1:length(A1_isch),
						[cx1 cy1 cx2 cy2] = plot_diag_two_conditions([1 2 3 4],[A1_isch(k) A1_csch(k) A1_csih(k) A1_isih(k)],...
							[1 2 3 4],[A2_isch(k) A2_csch(k) A2_csih(k) A2_isih(k)],'Marker','.','Color',map(angles_sum(k),:)); hold on
						
                        
                        
						% text(cx1,cy1,unit_IDs{k},'FontSize',5,'Interpreter','none'); % for debug
						
							
						
						if analyzeTuningIndices % tuning indices
							TI1_space(k) = csi(nanmean([A1_csch(k) A1_csih(k)]),nanmean([A1_isch(k) A1_isih(k)]));
							TI2_space(k) = csi(nanmean([A2_csch(k) A2_csih(k)]),nanmean([A2_isch(k) A2_isih(k)]));
							TI1_hand(k) = csi(nanmean([A1_csch(k) A1_isch(k)]),nanmean([A1_isih(k) A1_csih(k)]));
							TI2_hand(k) = csi(nanmean([A2_csch(k) A2_isch(k)]),nanmean([A2_isih(k) A2_csih(k)]));
							
							% now pref. (in condition 1) vs nonpref. (in condition 1)
							if nanmean([A1_csch(k) A1_csih(k)])>nanmean([A1_isch(k) A1_isih(k)]),
								TI1_pref_space(k) = TI1_space(k);
								TI2_pref_space(k) = TI2_space(k);
							else
								TI1_pref_space(k) = csi(nanmean([A1_isch(k) A1_isih(k)]),nanmean([A1_csch(k) A1_csih(k)]));
								TI2_pref_space(k) = csi(nanmean([A2_isch(k) A2_isih(k)]),nanmean([A2_csch(k) A2_csih(k)]));
							end
							
							if nanmean([A1_csch(k) A1_isch(k)])>nanmean([A1_isih(k) A1_csih(k)]),
								TI1_pref_hand(k) = TI1_hand(k);
								TI2_pref_hand(k) = TI2_hand(k);
							else
								TI1_pref_hand(k) = csi(nanmean([A1_isih(k) A1_csih(k)]),nanmean([A1_csch(k) A1_isch(k)]));
								TI2_pref_hand(k) = csi(nanmean([A2_isih(k) A2_csih(k)]),nanmean([A2_csch(k) A2_isch(k)]));
							end
						end % of if tuning indices
						
% 						if strcmp(unit_IDs{k},'Lin_20170707_07'), % for debug
% 							disp('Lin_20170707_07');
% 						end
						
					end
					
					% max1 = max(max([A1_isch; A1_csch; A1_csih; A1_isih]));
					% max2 = max(max([A2_isch; A2_csch; A2_csih; A2_isih]));
					% maxv = max([max1 max2]);
					%
					
					xlim = get(gca,'xlim');
					ylim = get(gca,'ylim');
					maxv=max(abs([xlim ylim]));
					set(gca,'xlim',[-maxv maxv]);
					set(gca,'ylim',[-maxv maxv]);
					
					axis square
					ig_add_diagonal_line;
					ig_add_zero_lines;
					xlim = get(gca,'xlim');
					ylim = get(gca,'ylim');
					line([xlim(1) xlim(2)],[ylim(2) ylim(1)],'Color','k','LineStyle',':');
					line(x_circle,y_circle,'Color','k','LineStyle',':')
					line(2*x_circle,2*y_circle,'Color','k','LineStyle',':')
					
					
					subplot(subplot_rows,numel(epochs),e);
					title([epoch ' ' num2str(numel(A2_isch))]);
					for k=1:length(A1_isch),
						plot_diag_rectangles([1 2 3 4],[A1_isch(k) A1_csch(k) A1_csih(k) A1_isih(k)],'Marker','.','Color',map(angles_sum(k),:)); hold on
					end
					line(x_circle,y_circle,'Color','k','LineStyle',':')
					
					subplot(subplot_rows,numel(epochs),numel(epochs)+e);
					
					for k=1:length(A1_isch),
						plot_diag_rectangles([1 2 3 4],[A2_isch(k) A2_csch(k) A2_csih(k) A2_isih(k)],'Marker','.','Color',map(angles_sum(k),:)); hold on
					end
					
					% summary
					sps_handle(e)=subplot(subplot_rows,numel(epochs),3*numel(epochs)+e);
					x1_average=nanmean(real(I1_mat));
					y1_average=nanmean(imag(I1_mat));
					x2_average=nanmean(real(I2_mat));
					y2_average=nanmean(imag(I2_mat));
					hold on
					plot(x1_average,y1_average,'o','color',[0.5 0.5 0.5]);
					line(x_circle,y_circle,'Color','k','LineStyle',':')
					%plot(x2_average,y2_average,'o');
					drawArrow([x1_average y1_average],[x2_average y2_average],[0.5 0.5 0.5]);
					% axis equal
					% axis square
					
					max_sum=max([max_sum abs(x2_average) abs(y2_average) abs(x1_average) abs(y1_average)]);
					
                   
                
                
					if analyzeTuningIndices,
						
						% uncomment below to use pref/nonpref
						% TI1_space = TI1_pref_space; TI2_space = TI2_pref_space; TI1_hand = TI1_pref_hand; TI2_hand = TI2_pref_hand;
						
						subplot(subplot_rows,numel(epochs),4*numel(epochs)+e); hold on
						for k=1:length(A1_isch),
							plot(TI1_space(k),TI2_space(k),'o','Color',map(angles_sum(k),:));
							% text(TI1_space(k),TI2_space(k),unit_IDs{k},'FontSize',5,'Interpreter','none'); % for debug
						end
						[h_space,p_space] = ttest(TI1_space,TI2_space);
						title(sprintf('d %.2f p %.2f',nanmean(TI2_space-TI1_space),p_space));
						set(gca,'Xlim',[-1 1],'Ylim',[-1 1]);
						axis square
						ig_add_diagonal_line;
						ig_add_zero_lines;
						
						subplot(subplot_rows,numel(epochs),5*numel(epochs)+e); hold on
						for k=1:length(A1_isch),
							plot(TI1_hand(k),TI2_hand(k),'o','Color',map(angles_sum(k),:));
							% text(TI1_hand(k),TI2_hand(k),unit_IDs{k},'FontSize',5,'Interpreter','none'); % for debug
						end
						[h_hand,p_hand] = ttest(TI1_hand,TI2_hand);
						title(sprintf('d %.2f p %.2f',nanmean(TI2_hand-TI1_hand),p_hand));
						set(gca,'Xlim',[-1 1],'Ylim',[-1 1]);
						axis square
						ig_add_diagonal_line;
						ig_add_zero_lines;
                    end
                    
                 
					
				end
				
				%% scaling (axis limits) across figures
				
				max_sum=0.7;
				
				for e=1:numel(epochs)
					
					subplot(sps_handle(e))
					set(gca,'xlim',[-max_sum max_sum],'ylim',[-max_sum max_sum]);
					ig_add_diagonal_line;
					ig_add_zero_lines;
					xlim = get(gca,'xlim');
					ylim = get(gca,'ylim');
					line([xlim(1) xlim(2)],[ylim(2) ylim(1)],'Color','k','LineStyle',':');
				end
				
				% subplot(ceil(sqrt(numel(epochs)+1)),ceil(sqrt(numel(epochs)+1)),e+1);
				%             subplot(2,numel(epochs)+1,e+1);
				%             colorbar;
				
				if ~exist([keys.drive keys.basepath_to_save filesep keys.project_version filesep 'tuning_inactivation'],'dir')
					mkdir([keys.drive keys.basepath_to_save filesep keys.project_version],'tuning_inactivation');
				end
				
				
				figure_title=[keys.monkey '_' target '_' tasktype '_normalized_' num2str(normalize_by_max_FR)];
				filename=[keys.drive keys.basepath_to_save filesep keys.project_version filesep 'tuning_inactivation' filesep figure_title];
				mtit(figure_title, 'xoff', 0, 'yoff', 0.05, 'color', [0 0 0], 'fontsize', 12,'Interpreter', 'none');
				export_fig(filename, '-pdf','-transparent') % pdf by run
                
                
				
			end
		end
	end
end

disp('Done.')

function [cx1 cy1 cx2 cy2] = plot_diag_two_conditions(diag1, mag1, diag2, mag2, varargin)
% diag 1 -> 2 -> 3 - 4 counterclockwise
% IS-CH CS-CH CS-IH IS-IH

% if ~all(size(diag) == size(mag)),
% end

for d = 1:length(diag1),
	switch diag(d),
		case 1
			x1(d) = sqrt(mag1(d)^2/2);
			y1(d) = sqrt(mag1(d)^2/2);
			x2(d) = sqrt(mag2(d)^2/2);
			y2(d) = sqrt(mag2(d)^2/2);
		case 2
			x1(d) = -sqrt(mag1(d)^2/2);
			y1(d) = sqrt(mag1(d)^2/2);
			x2(d) = -sqrt(mag2(d)^2/2);
			y2(d) = sqrt(mag2(d)^2/2);
		case 3
			x1(d) = -sqrt(mag1(d)^2/2);
			y1(d) = -sqrt(mag1(d)^2/2);
			x2(d) = -sqrt(mag2(d)^2/2);
			y2(d) = -sqrt(mag2(d)^2/2);
		case 4
			x1(d) = sqrt(mag1(d)^2/2);
			y1(d) = -sqrt(mag1(d)^2/2);
			x2(d) = sqrt(mag2(d)^2/2);
			y2(d) = -sqrt(mag2(d)^2/2);
	end
end
[dummy,sort_idx] = sort(diag1);
x1 = x1(sort_idx);
y1 = y1(sort_idx);
x2 = x2(sort_idx);
y2 = y2(sort_idx);

hold on
%
% h1 = plot([x1 x1(1)],[y1 y1(1)],'k-',varargin{:}); hold on;
% h2 = plot([x2 x2(1)],[y2 y2(1)],'k-',varargin{:}); hold on;
%
% set(h2,'Color',get(h1,'Color')*0.66);

%% ?????
% [cx1,cy1] = ig_comass(x1,y1);
% [cx2,cy2] = ig_comass(x2,y2);


cx1= sum(x1);
cy1= sum(y1);
cx2= sum(x2);
cy2= sum(y2);

% if ~(cx1 == 0 && cy1 == 0),
%     % hf1 = feather(cx1,cy1,'r');
%     hf1 = plot(cx1,cy1,'o');
%     set(hf1,varargin{end-1:end});
% end

% if ~(cx2 == 0 && cy2 == 0),
% 	% hf2 = feather(cx2,cy2,'r');
% 	hf2 = plot(cx2,cy2,'o');
% 	set(hf2,varargin{end-1:end});
% 	set(hf2,'Color',get(hf1,'Color')*0.66);
% end

drawArrow([cx1 cy1],[cx2 cy2],varargin{end});

cx1 = double(cx1);
cy1 = double(cy1);
cx2 = double(cx2);
cy2 = double(cy2);



function h = plot_diag_rectangles(diag1, mag1, varargin)
% diag 1 -> 2 -> 3 - 4 counterclockwise
% IS-CH CS-CH CS-IH IS-IH

% if ~all(size(diag) == size(mag)),
% end

for d = 1:length(diag1),
	switch diag(d),
		case 1
			x1(d) = sqrt(mag1(d)^2/2);
			y1(d) = sqrt(mag1(d)^2/2);
		case 2
			x1(d) = -sqrt(mag1(d)^2/2);
			y1(d) = sqrt(mag1(d)^2/2);
		case 3
			x1(d) = -sqrt(mag1(d)^2/2);
			y1(d) = -sqrt(mag1(d)^2/2);
		case 4
			x1(d) = sqrt(mag1(d)^2/2);
			y1(d) = -sqrt(mag1(d)^2/2);
	end
end
[dummy,sort_idx] = sort(diag1);
x1 = x1(sort_idx);
y1 = y1(sort_idx);

hold on
h = plot([x1 x1(1)],[y1 y1(1)],'k-',varargin{:}); hold on;



function h = plot_diag(diag, mag, varargin)

% diag 1 -> 2 -> 3 - 4 counterclockwise
% IS-CH CS-CH CS-IH IS-IH

% if ~all(size(diag) == size(mag)),
% end

for d = 1:length(diag),
	switch diag(d),
		case 1
			x(d) = sqrt(mag(d)^2/2);
			y(d) = sqrt(mag(d)^2/2);
		case 2
			x(d) = -sqrt(mag(d)^2/2);
			y(d) = sqrt(mag(d)^2/2);
		case 3
			x(d) = -sqrt(mag(d)^2/2);
			y(d) = -sqrt(mag(d)^2/2);
		case 4
			x(d) = sqrt(mag(d)^2/2);
			y(d) = -sqrt(mag(d)^2/2);
	end
end
[dummy,sort_idx] = sort(diag);
x = x(sort_idx);
y = y(sort_idx);

h = plot([x x(1)],[y y(1)],'k-',varargin{:}); hold on;
[cx,cy] = ig_comass(x,y);

if ~(cx == 0 && cy == 0),
	h = feather(cx,cy,'r');
	set(h,varargin{end-1:end});
end








function hArrow = drawArrow(p0,p1,color)
% drawArrow(p0,p1)
% Draws a simple arrow in 2D, from p0 to p1.
% from: https://de.mathworks.com/matlabcentral/fileexchange/55181-drawarrow by Matthew Kelly
%
% INPUTS:
%   p0 = [x0; y0] = position of the tail
%   p1 = [x1; y1] = position of the tip
%   color = arrow color. Optional: default is black
%       --> can be 'r','g','b','c','m','y','w', 'k' or a 1x3 color vector
%
% OUTPUTS:
%   hArrow = handle to the patch object representing the arrow
%
% Defaults:
if nargin == 2
	color = 'k';
end
% Parameters:
W1 = 0.08;   % half width of the arrow head, normalized by length of arrow
W2 = 0.014;  % half width of the arrow shaft
L1 = 0.18;   % Length of the arrow head, normalized by length of arrow
L2 = 0.13;  % Length of the arrow inset
% Unpack the tail and tip of the arrow
x0 = p0(1);
y0 = p0(2);
x1 = p1(1);
y1 = p1(2);
% Start by drawing an arrow from 0 to 1 on the x-axis
P = [...
	0, (1-L2), (1-L1), 1, (1-L1), (1-L2), 0;
	W2,    W2,     W1, 0,    -W1,    -W2, -W2];
% Scale,rotate, shift and plot:
dx = x1-x0;
dy = y1-y0;
Length = sqrt(dx*dx + dy*dy);
Angle = atan2(-dy,dx);
P = P.*repmat([Length;1],1,7);   %Scale Length*
P = [cos(Angle), sin(Angle); -sin(Angle), cos(Angle)]*P;  %Rotate
P = p0(:)*ones(1,7) + P;  %Shift
% Plot!
hArrow = patch(P(1,:), P(2,:),color);  axis equal;
set(hArrow,'EdgeColor',color);


function CSI = csi(Contra,Ipsi)
CSI = (Contra-Ipsi)./(Contra+Ipsi);
% CSI = (Contra-Ipsi)./max( [abs(Contra) abs(Ipsi)],[],2 );
% CSI = (Contra-Ipsi);

