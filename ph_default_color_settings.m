function [colors,linestyles]=ph_default_color_settings
% traces
colors.eye_ver         =[0.8 0 0];
colors.eye_hor         =[1 0 0];
colors.rhd_ver         =[0 0.8 0];
colors.rhd_hor         =[0 1 0];
colors.lhd_ver         =[0 0 0.8];
colors.lhd_hor         =[0 0 1];

% average PSTH across conditions
colors.AV   =[0 0 0];

% colors for 3 fixation positions and contra/ipsi (gaze paper)
colors.IF   =[236 32 38];
colors.MF   =[16 159 218];
colors.CF   =[247 148 36];

colors.IF_CS   =[236 32 38];
colors.MF_CS   =[16 159 218];
colors.CF_CS   =[247 148 36];

colors.IF_VS   =[236 32 38]*2/3;
colors.MF_VS   =[16 159 218]*2/3;
colors.CF_VS   =[247 148 36]*2/3;

colors.IF_IS   =[236 32 38]*1/3;
colors.MF_IS   =[16 159 218]*1/3;
colors.CF_IS   =[247 148 36]*1/3;

% cell count colors
colors.NO_AN   =[255 255 255];
colors.NO_TU   =[128 128 128];
colors.EP_EN   =[255  65  0];
colors.EP_BI   =[106 189  69];
colors.EP_SU   =[0   65 255];
colors.CR      =[255   0   0];
colors.UC      =[127   0   0];

% single cell PSTH colors per position 
colors.in   =[0 255 0];
colors.ch   =[0 128 0];

colors.in_LH   =[64 0 255];
colors.ch_LH   =[32 0 128];
colors.in_RH   =[255 255 0];
colors.ch_RH   =[128 128 0];

colors.in_AH   =[0 255 0];
colors.ch_AH   =[0 128 0];
colors.in_IH   =[64 0 255];
colors.ch_IH   =[32 0 128];
colors.in_CH   =[255 255 0];
colors.ch_CH   =[128 128 0];

colors.in_IS   =[0 255 255];
colors.ch_IS   =[0 128 128];
colors.in_CS   =[255 0 64];
colors.ch_CS   =[128 0 32];
colors.in_VS   =[20 20 20];
colors.ch_VS   =[50 50 50];
 
%% MP changed colors for choice because he added dotted lines for that
colors.in_AH=[255 102 0];
colors.ch_AH=[255 102 0];
colors.in_AH_PT=[204 51 0];
colors.ch_AH_PT=[204 51 0];

colors.in_AH_CS=[255 102 0];
colors.ch_AH_CS=[255 102 0];
colors.in_AH_IS=[0 119 255];
colors.ch_AH_IS=[0 119 255];

colors.in_AH_PT_CS=[204 51 0];
colors.ch_AH_PT_CS=[204 51 0];
colors.in_AH_PT_IS=[0 70 141];
colors.ch_AH_PT_IS=[0 70 141];

% now same for ipsi hand 
colors.in_IH=[255 0 255];
colors.ch_IH=[255 0 255];
colors.in_IH_PT=[204 0 204];
colors.ch_IH_PT=[204 0 204];

colors.in_IH_CS=[255 0 255];
colors.ch_IH_CS=[255 0 255];
colors.in_IH_IS=[0 128 255];
colors.ch_IH_IS=[0 128 255];

colors.in_IH_PT_CS=[204 0 204];
colors.ch_IH_PT_CS=[204 0 204];
colors.in_IH_PT_IS=[0 64 204];
colors.ch_IH_PT_IS=[0 64 204];

% now same for contra hand 
colors.in_CH=[255 128 0];
colors.ch_CH=[255 128 0];
colors.in_CH_PT=[204 64 0];
colors.ch_CH_PT=[204 64 0];

colors.in_CH_CS=[255 128 0];
colors.ch_CH_CS=[255 128 0];
colors.in_CH_IS=[0 255 0];
colors.ch_CH_IS=[0 155 0];

colors.in_CH_PT_CS=[204 64 0];
colors.ch_CH_PT_CS=[204 64 0];
colors.in_CH_PT_IS=[0 204 0];
colors.ch_CH_PT_IS=[0 204 0];

% %% population contra ipsi and vertical PSTH colors -overwriting MP settings for now
% colors.in_AH_CS=[255 0 64];
% colors.ch_AH_CS=[128 0 32];
% colors.in_IH_CS=[255 0 255];
% colors.ch_IH_CS=[128 0 128];
% colors.in_CH_CS=[255 128 0];
% colors.ch_CH_CS=[128 64 0];
% 
% colors.in_AH_IS=[0 255 255];
% colors.ch_AH_IS=[0 128 128];
% colors.in_IH_IS=[0 128 255];
% colors.ch_IH_IS=[0 64 128];
% colors.in_CH_IS=[0 255 0];
% colors.ch_CH_IS=[0 128 0];

colors.in_AH_VS=[150 150 150];
colors.ch_AH_VS=[80 80 80];
colors.in_IH_VS=[64 0 255];
colors.ch_IH_VS=[32 0 128];
colors.in_CH_VS=[255 255 0];
colors.ch_CH_VS=[128 128 0];

%% distractor task
temp_colors_I=autumn(18)*255;
temp_colors_C=winter(18)*255;
temp_colors_V=spring(18)*255;

% per position colors
colors.SS_TA_ER=temp_colors_C(1,:);
colors.TT_TA_ER=temp_colors_C(2,:);
colors.TD_TA_ER=temp_colors_C(3,:);
colors.SS_D1_ER=temp_colors_C(4,:);
colors.TT_D1_ER=temp_colors_C(5,:);
colors.TD_D1_ER=temp_colors_C(6,:);
colors.SS_D2_ER=temp_colors_C(7,:);
colors.TT_D2_ER=temp_colors_C(8,:);
colors.TD_D2_ER=temp_colors_C(9,:);

colors.SS_TA_SU=temp_colors_I(10,:);
colors.TT_TA_SU=temp_colors_I(11,:);
colors.TD_TA_SU=temp_colors_I(12,:);
colors.SS_D1_SU=temp_colors_I(13,:);
colors.TT_D1_SU=temp_colors_I(14,:);
colors.TD_D1_SU=temp_colors_I(15,:);
colors.SS_D2_SU=temp_colors_I(16,:);
colors.TT_D2_SU=temp_colors_I(17,:);
colors.TD_D2_SU=temp_colors_I(18,:);

% per Position & hemifield
colors.SS_TA_1H_ER=temp_colors_C(1,:);
colors.TT_TA_1H_ER=temp_colors_C(2,:);
colors.TD_TA_1H_ER=temp_colors_C(3,:);
colors.SS_D1_1H_ER=temp_colors_C(4,:);
colors.TT_D1_1H_ER=temp_colors_C(5,:);
colors.TD_D1_1H_ER=temp_colors_C(6,:);
colors.SS_D2_1H_ER=temp_colors_C(7,:);
colors.TT_D2_1H_ER=temp_colors_C(8,:);
colors.TD_D2_1H_ER=temp_colors_C(9,:);

colors.SS_TA_1H_SU=temp_colors_I(10,:);
colors.TT_TA_1H_SU=temp_colors_I(11,:);
colors.TD_TA_1H_SU=temp_colors_I(12,:);
colors.SS_D1_1H_SU=temp_colors_I(13,:);
colors.TT_D1_1H_SU=temp_colors_I(14,:);
colors.TD_D1_1H_SU=temp_colors_I(15,:);
colors.SS_D2_1H_SU=temp_colors_I(16,:);
colors.TT_D2_1H_SU=temp_colors_I(17,:);
colors.TD_D2_1H_SU=temp_colors_I(18,:);

colors.SS_TA_2H_ER=temp_colors_C(1,:);
colors.TT_TA_2H_ER=temp_colors_C(2,:);
colors.TD_TA_2H_ER=temp_colors_C(3,:);
colors.SS_D1_2H_ER=temp_colors_C(4,:);
colors.TT_D1_2H_ER=temp_colors_C(5,:);
colors.TD_D1_2H_ER=temp_colors_C(6,:);
colors.SS_D2_2H_ER=temp_colors_C(7,:);
colors.TT_D2_2H_ER=temp_colors_C(8,:);
colors.TD_D2_2H_ER=temp_colors_C(9,:);

colors.SS_TA_2H_SU=temp_colors_I(10,:);
colors.TT_TA_2H_SU=temp_colors_I(11,:);
colors.TD_TA_2H_SU=temp_colors_I(12,:);
colors.SS_D1_2H_SU=temp_colors_I(13,:);
colors.TT_D1_2H_SU=temp_colors_I(14,:);
colors.TD_D1_2H_SU=temp_colors_I(15,:);
colors.SS_D2_2H_SU=temp_colors_I(16,:);
colors.TT_D2_2H_SU=temp_colors_I(17,:);
colors.TD_D2_2H_SU=temp_colors_I(18,:);

% Contra vs IPSI
colors.SS_TA_ER_IS=temp_colors_I(1,:);
colors.TT_TA_ER_IS=temp_colors_I(2,:);
colors.TD_TA_ER_IS=temp_colors_I(3,:);
colors.SS_D1_ER_IS=temp_colors_I(4,:);
colors.TT_D1_ER_IS=temp_colors_I(5,:);
colors.TD_D1_ER_IS=temp_colors_I(6,:);
colors.SS_D2_ER_IS=temp_colors_I(7,:);
colors.TT_D2_ER_IS=temp_colors_I(8,:);
colors.TD_D2_ER_IS=temp_colors_I(9,:);

colors.SS_TA_ER_VS=temp_colors_V(1,:);
colors.TT_TA_ER_VS=temp_colors_V(2,:);
colors.TD_TA_ER_VS=temp_colors_V(3,:);
colors.SS_D1_ER_VS=temp_colors_V(4,:);
colors.TT_D1_ER_VS=temp_colors_V(5,:);
colors.TD_D1_ER_VS=temp_colors_V(6,:);
colors.SS_D2_ER_VS=temp_colors_V(7,:);
colors.TT_D2_ER_VS=temp_colors_V(8,:);
colors.TD_D2_ER_VS=temp_colors_V(9,:);

colors.SS_TA_ER_CS=temp_colors_C(1,:);
colors.TT_TA_ER_CS=temp_colors_C(2,:);
colors.TD_TA_ER_CS=temp_colors_C(3,:);
colors.SS_D1_ER_CS=temp_colors_C(4,:);
colors.TT_D1_ER_CS=temp_colors_C(5,:);
colors.TD_D1_ER_CS=temp_colors_C(6,:);
colors.SS_D2_ER_CS=temp_colors_C(7,:);
colors.TT_D2_ER_CS=temp_colors_C(8,:);
colors.TD_D2_ER_CS=temp_colors_C(9,:);

colors.SS_TA_SU_IS=temp_colors_I(10,:);
colors.TT_TA_SU_IS=temp_colors_I(11,:);
colors.TD_TA_SU_IS=temp_colors_I(12,:);
colors.SS_D1_SU_IS=temp_colors_I(13,:);
colors.TT_D1_SU_IS=temp_colors_I(14,:);
colors.TD_D1_SU_IS=temp_colors_I(15,:);
colors.SS_D2_SU_IS=temp_colors_I(16,:);
colors.TT_D2_SU_IS=temp_colors_I(17,:);
colors.TD_D2_SU_IS=temp_colors_I(18,:);

colors.SS_TA_SU_VS=temp_colors_V(10,:);
colors.TT_TA_SU_VS=temp_colors_V(11,:);
colors.TD_TA_SU_VS=temp_colors_V(12,:);
colors.SS_D1_SU_VS=temp_colors_V(13,:);
colors.TT_D1_SU_VS=temp_colors_V(14,:);
colors.TD_D1_SU_VS=temp_colors_V(15,:);
colors.SS_D2_SU_VS=temp_colors_V(16,:);
colors.TT_D2_SU_VS=temp_colors_V(17,:);
colors.TD_D2_SU_VS=temp_colors_V(18,:);

%not sure why different color scheme here
col_left      = autumn(6);

colors.SS_TA_SU_CS= col_left(1,:);
colors.TT_TA_SU_CS= col_left(1,:);
colors.TD_TA_SU_CS= col_left(1,:);
colors.SS_D1_SU_CS=col_left(6,:);
colors.TT_D1_SU_CS=col_left(6,:);
colors.TD_D1_SU_CS=col_left(6,:);
colors.SS_D2_SU_CS= col_left(3,:);
colors.TT_D2_SU_CS=col_left(3,:);
colors.TD_D2_SU_CS=col_left(3,:);

% colors.SS_TA_SU_CS=temp_colors_C(10,:);
% colors.TT_TA_SU_CS=temp_colors_C(11,:);
% colors.TD_TA_SU_CS=temp_colors_C(12,:);
% colors.SS_D1_SU_CS=temp_colors_C(13,:);
% colors.TT_D1_SU_CS=temp_colors_C(14,:);
% colors.TD_D1_SU_CS=temp_colors_C(15,:);
% colors.SS_D2_SU_CS=temp_colors_C(16,:);
% colors.TT_D2_SU_CS=temp_colors_C(17,:);
% colors.TD_D2_SU_CS=temp_colors_C(18,:);

% additonally separate by distractor position: same(1H)/opposite(2H) hemifield
colors.SS_TA_1H_ER_IS=temp_colors_I(1,:);
colors.TT_TA_1H_ER_IS=temp_colors_I(2,:);
colors.TD_TA_1H_ER_IS=temp_colors_I(3,:);
colors.SS_D1_1H_ER_IS=temp_colors_I(4,:);
colors.TT_D1_1H_ER_IS=temp_colors_I(5,:);
colors.TD_D1_1H_ER_IS=temp_colors_I(6,:);
colors.SS_D2_1H_ER_IS=temp_colors_I(7,:);
colors.TT_D2_1H_ER_IS=temp_colors_I(8,:);
colors.TD_D2_1H_ER_IS=temp_colors_I(9,:);

colors.SS_TA_1H_ER_VS=temp_colors_V(1,:);
colors.TT_TA_1H_ER_VS=temp_colors_V(2,:);
colors.TD_TA_1H_ER_VS=temp_colors_V(3,:);
colors.SS_D1_1H_ER_VS=temp_colors_V(4,:);
colors.TT_D1_1H_ER_VS=temp_colors_V(5,:);
colors.TD_D1_1H_ER_VS=temp_colors_V(6,:);
colors.SS_D2_1H_ER_VS=temp_colors_V(7,:);
colors.TT_D2_1H_ER_VS=temp_colors_V(8,:);
colors.TD_D2_1H_ER_VS=temp_colors_V(9,:);

colors.SS_TA_1H_ER_CS=temp_colors_C(1,:);
colors.TT_TA_1H_ER_CS=temp_colors_C(2,:);
colors.TD_TA_1H_ER_CS=temp_colors_C(3,:);
colors.SS_D1_1H_ER_CS=temp_colors_C(4,:);
colors.TT_D1_1H_ER_CS=temp_colors_C(5,:);
colors.TD_D1_1H_ER_CS=temp_colors_C(6,:);
colors.SS_D2_1H_ER_CS=temp_colors_C(7,:);
colors.TT_D2_1H_ER_CS=temp_colors_C(8,:);
colors.TD_D2_1H_ER_CS=temp_colors_C(9,:);

colors.SS_TA_2H_ER_IS=temp_colors_I(1,:);
colors.TT_TA_2H_ER_IS=temp_colors_I(2,:);
colors.TD_TA_2H_ER_IS=temp_colors_I(3,:);
colors.SS_D1_2H_ER_IS=temp_colors_I(4,:);
colors.TT_D1_2H_ER_IS=temp_colors_I(5,:);
colors.TD_D1_2H_ER_IS=temp_colors_I(6,:);
colors.SS_D2_2H_ER_IS=temp_colors_I(7,:);
colors.TT_D2_2H_ER_IS=temp_colors_I(8,:);
colors.TD_D2_2H_ER_IS=temp_colors_I(9,:);

colors.SS_TA_2H_ER_VS=temp_colors_V(1,:);
colors.TT_TA_2H_ER_VS=temp_colors_V(2,:);
colors.TD_TA_2H_ER_VS=temp_colors_V(3,:);
colors.SS_D1_2H_ER_VS=temp_colors_V(4,:);
colors.TT_D1_2H_ER_VS=temp_colors_V(5,:);
colors.TD_D1_2H_ER_VS=temp_colors_V(6,:);
colors.SS_D2_2H_ER_VS=temp_colors_V(7,:);
colors.TT_D2_2H_ER_VS=temp_colors_V(8,:);
colors.TD_D2_2H_ER_VS=temp_colors_V(9,:);

colors.SS_TA_2H_ER_CS=temp_colors_C(1,:);
colors.TT_TA_2H_ER_CS=temp_colors_C(2,:);
colors.TD_TA_2H_ER_CS=temp_colors_C(3,:);
colors.SS_D1_2H_ER_CS=temp_colors_C(4,:);
colors.TT_D1_2H_ER_CS=temp_colors_C(5,:);
colors.TD_D1_2H_ER_CS=temp_colors_C(6,:);
colors.SS_D2_2H_ER_CS=temp_colors_C(7,:);
colors.TT_D2_2H_ER_CS=temp_colors_C(8,:);
colors.TD_D2_2H_ER_CS=temp_colors_C(9,:);


colors.SS_TA_1H_SU_IS=temp_colors_I(10,:);
colors.TT_TA_1H_SU_IS=temp_colors_I(11,:);
colors.TD_TA_1H_SU_IS=temp_colors_I(12,:);
colors.SS_D1_1H_SU_IS=temp_colors_I(13,:);
colors.TT_D1_1H_SU_IS=temp_colors_I(14,:);
colors.TD_D1_1H_SU_IS=temp_colors_I(15,:);
colors.SS_D2_1H_SU_IS=temp_colors_I(16,:);
colors.TT_D2_1H_SU_IS=temp_colors_I(17,:);
colors.TD_D2_1H_SU_IS=temp_colors_I(18,:);

colors.SS_TA_1H_SU_VS=temp_colors_V(10,:);
colors.TT_TA_1H_SU_VS=temp_colors_V(11,:);
colors.TD_TA_1H_SU_VS=temp_colors_V(12,:);
colors.SS_D1_1H_SU_VS=temp_colors_V(13,:);
colors.TT_D1_1H_SU_VS=temp_colors_V(14,:);
colors.TD_D1_1H_SU_VS=temp_colors_V(15,:);
colors.SS_D2_1H_SU_VS=temp_colors_V(16,:);
colors.TT_D2_1H_SU_VS=temp_colors_V(17,:);
colors.TD_D2_1H_SU_VS=temp_colors_V(18,:);

colors.SS_TA_1H_SU_CS=temp_colors_C(10,:);
colors.TT_TA_1H_SU_CS=temp_colors_C(11,:);
colors.TD_TA_1H_SU_CS=temp_colors_C(12,:);
colors.SS_D1_1H_SU_CS=temp_colors_C(13,:);
colors.TT_D1_1H_SU_CS=temp_colors_C(14,:);
colors.TD_D1_1H_SU_CS=temp_colors_C(15,:);
colors.SS_D2_1H_SU_CS=temp_colors_C(16,:);
colors.TT_D2_1H_SU_CS=temp_colors_C(17,:);
colors.TD_D2_1H_SU_CS=temp_colors_C(18,:);

colors.SS_TA_2H_SU_IS=temp_colors_I(10,:);
colors.TT_TA_2H_SU_IS=temp_colors_I(11,:);
colors.TD_TA_2H_SU_IS=temp_colors_I(12,:);
colors.SS_D1_2H_SU_IS=temp_colors_I(13,:);
colors.TT_D1_2H_SU_IS=temp_colors_I(14,:);
colors.TD_D1_2H_SU_IS=temp_colors_I(15,:);
colors.SS_D2_2H_SU_IS=temp_colors_I(16,:);
colors.TT_D2_2H_SU_IS=temp_colors_I(17,:);
colors.TD_D2_2H_SU_IS=temp_colors_I(18,:);

colors.SS_TA_2H_SU_VS=temp_colors_V(10,:);
colors.TT_TA_2H_SU_VS=temp_colors_V(11,:);
colors.TD_TA_2H_SU_VS=temp_colors_V(12,:);
colors.SS_D1_2H_SU_VS=temp_colors_V(13,:);
colors.TT_D1_2H_SU_VS=temp_colors_V(14,:);
colors.TD_D1_2H_SU_VS=temp_colors_V(15,:);
colors.SS_D2_2H_SU_VS=temp_colors_V(16,:);
colors.TT_D2_2H_SU_VS=temp_colors_V(17,:);
colors.TD_D2_2H_SU_VS=temp_colors_V(18,:);

colors.SS_TA_2H_SU_CS=temp_colors_C(10,:);
colors.TT_TA_2H_SU_CS=temp_colors_C(11,:);
colors.TD_TA_2H_SU_CS=temp_colors_C(12,:);
colors.SS_D1_2H_SU_CS=temp_colors_C(13,:);
colors.TT_D1_2H_SU_CS=temp_colors_C(14,:);
colors.TD_D1_2H_SU_CS=temp_colors_C(15,:);
colors.SS_D2_2H_SU_CS=temp_colors_C(16,:);
colors.TT_D2_2H_SU_CS=temp_colors_C(17,:);
colors.TD_D2_2H_SU_CS=temp_colors_C(18,:);

%% preferred and unpreferred same as contra/ipsi + add different perturbation block colors
color_fieldnames=fieldnames(colors);
for fn=1:numel(color_fieldnames)
    switch color_fieldnames{fn}(end-1:end)
        case 'CS'
            colors.([color_fieldnames{fn}(1:end-2) 'PF'])=colors.(color_fieldnames{fn});
        case 'IS'
            colors.([color_fieldnames{fn}(1:end-2) 'NP'])=colors.(color_fieldnames{fn});
    end
    pt_index=strfind(color_fieldnames{fn},'PT');
    if any(pt_index)
        colors.([color_fieldnames{fn}(1:pt_index+1) '2' color_fieldnames{fn}(pt_index+2:end)])=colors.(color_fieldnames{fn})*7/8;
        colors.([color_fieldnames{fn}(1:pt_index+1) '3' color_fieldnames{fn}(pt_index+2:end)])=colors.(color_fieldnames{fn})*6/8;
        colors.([color_fieldnames{fn}(1:pt_index+1) '4' color_fieldnames{fn}(pt_index+2:end)])=colors.(color_fieldnames{fn})*5/8;
        colors.([color_fieldnames{fn}(1:pt_index+1) '5' color_fieldnames{fn}(pt_index+2:end)])=colors.(color_fieldnames{fn})*4/8;
        colors.([color_fieldnames{fn}(1:pt_index+1) '6' color_fieldnames{fn}(pt_index+2:end)])=colors.(color_fieldnames{fn})*3/8;
        colors.([color_fieldnames{fn}(1:pt_index+1) '7' color_fieldnames{fn}(pt_index+2:end)])=colors.(color_fieldnames{fn})*2/8;
        colors.([color_fieldnames{fn}(1:pt_index+1) '8' color_fieldnames{fn}(pt_index+2:end)])=colors.(color_fieldnames{fn})*1/8;
    end
end

%% population effector colors
colors.EF_SA   =[0 255  0];
colors.EF_RE   =[0 255  0];
colors.EF_FG   =[0 0  255];

colors.in_PF   =colors.EP_EN;
colors.in_NP   =colors.EP_SU;

%% overlapping tt and cc (??)
colors.per_monkey          =[0 1 0; 1 0 0];

%% linestyles for choices: dotted lines!
color_fieldnames=fieldnames(colors);
for fn=1:numel(color_fieldnames)
    if strcmp(color_fieldnames{fn}(1:2),'ch')        
        linestyles.(color_fieldnames{fn})=':';    
    else
        linestyles.(color_fieldnames{fn})='-';
    end
end

end