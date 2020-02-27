function ph_hand_space_tuning_vector(keys)
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
%
% keys=struct;
% keys=ph_general_settings(project,keys);



angles=0:pi/180:2*pi;
x_circle=cos(angles)*0.5;
y_circle=sin(angles)*0.5;
phi=pi/4;
R_mat=[sin(phi)+1i*cos(phi), -sin(phi)+1i*cos(phi), -sin(phi)-1i*cos(phi), sin(phi)-1i*cos(phi)];
%
% keys.tt.IC_for_criterion            = 'in';
%
%
% project_specific_settings=[keys.db_folder filesep project filesep 'ph_project_settings.m'];
% run(project_specific_settings);
%
% keys.project_versions=versions;
analyzeTuningIndices = false; %true;

if analyzeTuningIndices
    subplot_rows = 6;
else
    subplot_rows = 4;
end

[tuning_per_unit_table]                 = ph_load_extended_tuning_table(keys);
[tuning_per_unit_table, Sel_for_title]  = ph_reduce_tuning_table(tuning_per_unit_table,keys);


for kk=1:size(tuning_per_unit_table,2)
    idx.(tuning_per_unit_table{1,kk})=kk;
end

epochs=keys.HS.epochs;
% if any(strfind(keys.project_version,'MIP'))
%     epochs={'Facq','Fhol','Cue','Del','PreR','PeriR'};
%     % epochs={'Del'}; % for debug
% elseif  any(strfind(keys.project_version,'LIP'))
%     epochs={'Facq','Fhol','Cue','Del','PreS','PeriS'};
%     % epochs={'Cue','Del','PeriS'}; % for debug
% end


for normalize_by_max_FR=0:1
    for tt=1:numel(keys.HS.tasktypes)
        tasktype=[keys.HS.tasktypes{tt} '_' keys.arrangement(1:3)];
        % variant 1: one task, two conditions (e.g. pre- and post-injection), two epochs: "baseline" and "response",
        % four trial types (contralesional and ipsilesional hand and space)
        
        plot_title_part=['normalized_' num2str(normalize_by_max_FR)];
        fig_title=sprintf('%s %s %s hnd %s ch %s %s, %s,%s',...
            keys.monkey,[keys.conditions_to_plot{:}],keys.arrangement,     mat2str(keys.tt.hands),mat2str(double(keys.tt.choices)),[Sel_for_title{:}],plot_title_part,tasktype);
        f_handle=figure('units','normalized','outerposition',[0 0 1 1],'name',fig_title);
        max_sum=0;
        for e=1:numel(epochs)
            epoch=epochs{e};
            if ismember(1,keys.tt.perturbation)
                A1_isch = [tuning_per_unit_table{2:end, idx.(['in_IH_CS_' epoch '_CT_FR_' tasktype])}];
                A1_csch = [tuning_per_unit_table{2:end, idx.(['in_IH_IS_' epoch '_CT_FR_' tasktype])}];
                A1_csih = [tuning_per_unit_table{2:end, idx.(['in_CH_IS_' epoch '_CT_FR_' tasktype])}];
                A1_isih = [tuning_per_unit_table{2:end, idx.(['in_CH_CS_' epoch '_CT_FR_' tasktype])}];
                
                A2_isch = [tuning_per_unit_table{2:end, idx.(['in_IH_CS_' epoch '_PT_FR_' tasktype])}];
                A2_csch = [tuning_per_unit_table{2:end, idx.(['in_IH_IS_' epoch '_PT_FR_' tasktype])}];
                A2_csih = [tuning_per_unit_table{2:end, idx.(['in_CH_IS_' epoch '_PT_FR_' tasktype])}];
                A2_isih = [tuning_per_unit_table{2:end, idx.(['in_CH_CS_' epoch '_PT_FR_' tasktype])}];
            else
                A1_isch = [tuning_per_unit_table{2:end, idx.(['in_IH_CS_' epoch '_epoch_FR_' tasktype])}];
                A1_csch = [tuning_per_unit_table{2:end, idx.(['in_IH_IS_' epoch '_epoch_FR_' tasktype])}];
                A1_csih = [tuning_per_unit_table{2:end, idx.(['in_CH_IS_' epoch '_epoch_FR_' tasktype])}];
                A1_isih = [tuning_per_unit_table{2:end, idx.(['in_CH_CS_' epoch '_epoch_FR_' tasktype])}];
                
                A2_isch = A1_isch;
                A2_csch = A1_csch;
                A2_csih = A1_csih;
                A2_isih = A1_isih;
            end
            
            %% normalize
            if normalize_by_max_FR % normalize by max rate in condtion 1 (control)  in each cell
                max_per_cell = max([A1_isch; A1_csch; A1_csih; A1_isih]);
                A1_isch = A1_isch./max_per_cell;
                A1_csch = A1_csch./max_per_cell;
                A1_csih = A1_csih./max_per_cell;
                A1_isih = A1_isih./max_per_cell;
                A2_isch = A2_isch./max_per_cell;
                A2_csch = A2_csch./max_per_cell;
                A2_csih = A2_csih./max_per_cell;
                A2_isih = A2_isih./max_per_cell;
            end
            
            
            %% rotating in right position
            A1_mat=[A1_isch' A1_csch' A1_csih' A1_isih'];
            I1_mat=R_mat*A1_mat';
            angles_sum=round(180/pi*angle(I1_mat))+180; angles_sum(isnan(angles_sum))=1; angles_sum(angles_sum==0)=360;
            A2_mat=[A2_isch' A2_csch' A2_csih' A2_isih'];
            I2_mat=R_mat*A2_mat';
            
            %unit_IDs=tuning_per_unit_table(2:end, idx.unit_ID);
            map = hsv(360);
            
            %% plot tuning rectangle per unit
            subplot(subplot_rows,numel(epochs),e);
            title([epoch ' ' num2str(numel(A1_isch))]);
            for k=1:length(A1_isch),
                plot_diag_rectangles([1 2 3 4],[A1_isch(k) A1_csch(k) A1_csih(k) A1_isih(k)],'Marker','.','Color',map(angles_sum(k),:)); hold on
            end
            line(x_circle,y_circle,'Color','k','LineStyle',':')
            
            subplot(subplot_rows,numel(epochs),numel(epochs)+e);
            title([epoch ' ' num2str(numel(A2_isch))]);
            for k=1:length(A1_isch),
                plot_diag_rectangles([1 2 3 4],[A2_isch(k) A2_csch(k) A2_csih(k) A2_isih(k)],'Marker','.','Color',map(angles_sum(k),:)); hold on
            end
            
            %% plot tuning vector (or tuning vector disposition after inactivation) per unit
            
            sp1_handle(e)=subplot(subplot_rows,numel(epochs),2*numel(epochs)+e);
            scaling_max_1(e)=max(max(abs([A1_isch A1_csch A1_csih A1_isih A2_isch A2_csch A2_csih A2_isih])));
            title(epoch);
            
            for k=1:length(A1_isch),
                if ismember(1,keys.tt.perturbation)
                    plot_diag_two_conditions([1 2 3 4],[A1_isch(k) A1_csch(k) A1_csih(k) A1_isih(k)],[1 2 3 4],[A2_isch(k) A2_csch(k) A2_csih(k) A2_isih(k)],'Marker','.','Color',map(angles_sum(k),:)); hold on
                else
                    plot_diag_one_condition([1 2 3 4],[A1_isch(k) A1_csch(k) A1_csih(k) A1_isih(k)],'MarkerEdgeColor',map(angles_sum(k),:),'MarkerFaceColor','none','linestyle','none','Markersize',5,'Marker','o'); hold on
                end
            end
            
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
            
            
            %% summary across units
            sps_handle(e)=subplot(subplot_rows,numel(epochs),3*numel(epochs)+e);
            if ismember(1,keys.tt.perturbation)
                x1_average=nanmean(real(I1_mat));
                y1_average=nanmean(imag(I1_mat));
                x2_average=nanmean(real(I2_mat));
                y2_average=nanmean(imag(I2_mat));
            else
                
                x1_average=0;
                y1_average=0;
                x2_average=nanmean(real(I2_mat));
                y2_average=nanmean(imag(I2_mat));
            end
            hold on
            plot(x1_average,y1_average,'o','color',[0.5 0.5 0.5]);
            line(x_circle,y_circle,'Color','k','LineStyle',':')
            %plot(x2_average,y2_average,'o');
            drawArrow([x1_average y1_average],[x2_average y2_average],[0.5 0.5 0.5]);
            % axis equal
            % axis square
            
            max_sum=max([max_sum abs(x2_average) abs(y2_average) abs(x1_average) abs(y1_average)]);
            
            
            
            
            if analyzeTuningIndices % tuning indices
                TI1_space = [];	TI2_space = []; TI1_hand = []; TI2_hand = [];
                TI1_pref_space = [];	TI2_pref_space = []; TI1_pref_hand = []; TI2_pref_hand = [];
                
                for k=1:length(A1_isch),
                    
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
                end
                
                
                
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
        
        %max_sum=1; %0.7;
        for e=1:numel(epochs)
            subplot(sps_handle(e))
            set(gca,'xlim',[-max_sum max_sum],'ylim',[-max_sum max_sum]);
            ig_add_diagonal_line;
            ig_add_zero_lines;
            xlim = get(gca,'xlim');
            ylim = get(gca,'ylim');
            line([xlim(1) xlim(2)],[ylim(2) ylim(1)],'Color','k','LineStyle',':');
        end
        
        max_sum=max(scaling_max_1);
        for e=1:numel(epochs)
            set(sp1_handle(e),'xlim',[-max_sum max_sum],'ylim',[-max_sum max_sum]);
        end
        
        
        % subplot(ceil(sqrt(numel(epochs)+1)),ceil(sqrt(numel(epochs)+1)),e+1);
        %             subplot(2,numel(epochs)+1,e+1);
        %             colorbar;
        
        filename=[fig_title];
        ph_title_and_save(f_handle,  filename,fig_title,keys);
    end
end
%disp('Done.')

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



function plot_diag_one_condition(diag1, mag1, varargin)
for d = 1:length(diag1),
    x1(d) = sqrt(mag1(d)^2/2);
    y1(d) = sqrt(mag1(d)^2/2);
    switch diag(d),
        case 2
            x1(d) = -x1(d);
        case 3
            x1(d) = -x1(d);
            y1(d) = -y1(d);
        case 4
            y1(d) = -y1(d);
    end
end
[~,sort_idx] = sort(diag1);
x1 = x1(sort_idx);
y1 = y1(sort_idx);

hold on
cx1= sum(x1);
cy1= sum(y1);

plot(cx1,cy1,varargin{:});

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

