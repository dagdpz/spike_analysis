%% general batching
keys.position_and_plotting_arrangements ={'hands'};                 % defines position batching and which conditions go into different figures/lines
keys.batching.monkeys                    ={'Curius','Linus'};       % monkeys on the project
keys.batching.combine_monkeys            =0;                        % for population analysis, monkeys can be separately or combined
keys.batching.targets                    ={'dPulv_r','dPulv_l'};    % to combine both hemispheres, just put {'dPulv'}
keys.batching.Subregions_separately      =0;                        % subregions can be processed independently (if only defining target is too crude)
keys.batching.Subregions{1}{1}           =struct('monkey',{''},'target',{''},'grid_x',{NaN},'grid_y',{NaN},'z_min',{NaN},'z_max',{NaN}); % subregion definitions by  z range for each grid hole
keys.batching.n_Subregions               =numel(keys.batching.Subregions);

keys.cal.perturbation_groups            ={0,[2,3,4,5,6,7,8]}; 

%% tuning table selection
% Most important parameter to understand how cell seleciton works in the population analysis, is the trial criterion per condition:
% units with less than this amount of trials for one of the conditions will be excluded. 
keys.cal.min_trials_per_condition       =5;
% This criterion is applied every time the tuning table is loaded via ph_load_extended_tuning_table and ph_reduce_tuning_table. 
% The tricky part here is to define the conditions in which we ask for a specific number of trials. 
% Each condition is defined by tasktype(contains type,effector,and
% position_arrangement),position (depends on arrangement),hand used and choice/instruced trial (MISSING HERE: PERTURBATION)
keys.tt.tasktypes                   ={'Ddre_han','Ddsa_han'}; % typically only one tasktype defines
keys.tt.hands                       =[1 2];
keys.tt.perturbations               =[0 1]; % from keys.cal.perturbation_groups: 0 is first group, 1 is second group and so on
keys.tt.choices                     =[0 1]; %IMPORTANT: and also not really perfect, for choice trials trial criterion is applied by hemifield, not by position.
% Each unique combination of the above parameters has to contain at least keys.cal.min_trials_per_condition trials, if not the cell is excluded in ph_reduce_tuning_table
keys.tt.selection                   ={'target','dPulv';                         % easy to use if there is a parameter in the tuning table for which you want your cells to have the same value
                                      'in_NH_TIhol_position_Msac_opt','true'};  % each row in the cell arryáy will be used to exclude cells that don't have the specifie characteristic
keys.tt.unselect                    ={'in_NH_TIhol_position_Msac_opt','false'}; % easy to use if there is a parameter in the tuning table for which your cells shouldn't have a specific value
keys.tt.combine_tuning_properties   ={'hand_tuning','in_IH_Facq_epoch_Ddre_han','in_CH_Facq_epoch_Ddre_han'}; % results in enen/ensu/suen/en-/-en/su-/-su/--
% ph_load_extended_tuning_table will create an additional column from combining existing columns: 
% 1st entry is the name of the new column (to refer to it), and the following entries specify the columns which should be combined.
% this additional column can be used for grouping, specific selection and/or unselection ph_load_extended_tuning_table also creates new columns which can be used
% for multicomparison correction (sort of) by excluding cells that did not show certain combination of ANOVA effects (this is ONLY used in cell_counts so far
keys.tt.epoch_criterion             ='none'; % only relevant for cell counts
keys.tt.space_criterion             ='none';
keys.tt.hands_criterion             ='none';
keys.tt.SXH_criterion               ='none';

%% tuning table format
%% NH ---> AH


% It is crucial to know how different tuning table (ANOVA results) entries are coded:


% !main ANOVA effects and interactions
%               example : 'in_epoch_main_Msac_opt'
% column name fragments : in/ch               _ epoch_main/spaceLR_main/hands_main/ExS/ExH/SxH/ExSxH              _ 'tasktype'            _ 'arrangement'
%               meaning : instructed/choice   _ parameter tested (see below)                                      _ short type & effector _ first three letters of position arrangement
%            Parameters : epoch_main/spaceLR_main/hands_main | ExS/ExH/SxH/ExSxH - interactions E(epoch)/S(hemifield)/H(hand) 
%       Possible values : true,false                         | true,false                                             

% !main ANOVA effects and interactions per hand
%               example : 'in_IH_epoch_main_Msac_opt' 
% column name fragments : in/ch               _ AH/IH/CH              _ epoch_main/position_main/fixation_main   _ 'tasktype'            _ 'arrangement'
%               meaning : instructed/choice   _ any/ipsi/contra hand  _ parameter tested (see below)             _ short type & effector _ first three letters of position arrangement
%            Parameters : epoch_main/spaceLR_main/hands_main | ExS/ExH/SxH/ExSxH - interactions E(epoch)/S(hemifield)/H(hand) 
%       Possible values : true,false                         | true,false                                             

% !per epoch tests
%               example : 'in_Thol_spaceLR_Msac_opt'   
% column name fragments : in/ch               _ 'epoch'                          _ spaceLR/hands/SxH             _ 'tasktype'            _ 'arrangement'
%               meaning : instructed/choice   _ name defined in EPOCHS_PER_TYPE  _ parameter tested (see below)  _ short type & effector _ first three letters of position arrangement
%            Parameters : spaceLR (hemifield preference) | hands(hand preference)     | SxH (interaction)
%       Possible values : IS/CS/- (ipsi/contra hemifield)| IH/CH/- (ipsi/contra hand) | CR/UC (crossed/uncrossed preferred)
% additional extensions : _DF/_IX                        | _DF/_IX                    | _DF/_IX  

% !per epoch tests per hand
%               example : 'in_AH_TIhol_epoch_Msac_opt' 
% column name fragments : in/ch               _ AH/IH/CH              _ 'epoch'                          _ epoch/position/prefH/gaze_modulation   _ 'tasktype'            _ 'arrangement'
%               meaning : instructed/choice   _ any/ipsi/contra hand  _ name defined in EPOCHS_PER_TYPE  _ parameter tested (see below)           _ short type & effector _ first three letters of position arrangement
%            Parameters : epoch (compared to respective baseline) | fixation/position/PxF (ANOVA main effects/interaction) | prefH/prefP (for preferred hemifield/position) | gaze_modulation
%       Possible values : /en/su/-                                | true,false                                             | in/ch/-  (intructed or choice)                 | monotoneus/nonmonotoneous
% additional extensions : _DF/_SC/bilateral                       |                                                        | _FR/I_FR/O_FR

% !per epoch tests per hand and hemifield
%               example : 'in_IH_CS_TIhol_epoch_Msac_opt' 
% column name fragments : in/ch               _ AH/IH/CH              _ IS/CS     _ 'epoch'                         _ epoch/PT/CT/PTbl             _ 'tasktype'            _ 'arrangement'
%               meaning : instructed/choice   _ any/ipsi/contra hand  _ hemifield _ name defined in EPOCHS_PER_TYPE _ parameter tested (see below) _ short type & effector _ first three letters of position arrangement
%            Parameters : epoch (compared to respective baseline) | PT (effect of perturbation) | PTbl (taking epoch baselines into account) | CT (control) 
%       Possible values : en/su/-                                 | EN/SU/-                     | EN/SU/-                                    |              
% additional extensions : _DF/_EN                                 |_FR/_DF                      | _DF                                        | _FR          

%            EXTENSIONS : | _FR                                 | _DF  (contra-ipsi, epoch-baseline)                      | _IX                     | I_FR                          | O_FR
%               meaning : average firing rate (accross trials)  | FR difference (A-B) between the two tested conditions   | index=(A-B)/(A+B)       | firing rate preferred (in)    | firing rate nonpreferred (out)

%% population PSTH settings 

cc=0;
% 1
cc=cc+1; % one increment of this index per plot                                        
% The following parameters will be assigned to keys.tt (please check in the previous section)
keys.pop(cc).tt.hands                      = [1 2];    % also used to exclude trials (not only units) to be plotted
keys.pop(cc).tt.perturbations              = [0 1];    % also used to exclude trials (not only units) to be plotted
keys.pop(cc).tt.choices                    = 0;        % also used to exclude trials (not only units) to be plotted
keys.pop(cc).tt.combine_tuning_properties  ={'hand_tuning','in_IH_Facq_epoch_Ddre_han','in_CH_Facq_epoch_Ddre_han'};
keys.pop(cc).tt.selection                  ={};
keys.pop(cc).tt.unselect                   ={};
keys.pop(cc).tt.tasktypes                  ={'Ddre_han'};


keys.pop(cc).group_parameter            = 'hand_tuning';        % defines the parameter from the tuning table by which population PSTHs will be grouped (rows)
keys.pop(cc).group_excluded             = {'','-','bi'};        % exclude units with these values for the selected parameter (group_parameter)
keys.pop(cc).epoch_PF                   = 'PreR';               % epoch in which preference defines target location for "pref" plots
keys.pop(cc).epoch_RF                   = 'PreR';               % epoch for which gaussian response fields will be plotted (if plot_RF ~ 0)
keys.pop(cc).epoch_for_normalization    = 'Fhol';               % epoch used for (divisive) normalization
keys.pop(cc).normalization              = 'by_effector';        % separate (divisive) normalization factor for trials grouped by effector; other options:
                                                                % 'by_condition','by_effector','by_type','by_all_trials','z_score','none'
keys.pop(cc).epoch_BL                   = 'INI';                % Epoch to subtract trial by trial (if FR_subtract_baseline ~ 0)
keys.pop(cc).FR_subtract_baseline       = 0;                    % trial by trial subtraction of the FR in the epoch defined by epoch_BL
keys.pop(cc).plot_RF                    = 1;                    % calculate and plot gaussian RFs
keys.pop(cc).y_lim                      = [0 100];              % y limts can be hardcoded to have the same scale across multiple plots

keys.pop(cc).conditions_to_plot         = {'Ddre'};             % although this parameter is designed for multiple conditions, only one entry here is recommended
keys.pop(cc).IC_for_criterion           = 'in';                 % used? for?



 
%% response timing keys (for figures saved in response timing folder formerly known as categorization)
cc=cc+1;% again, one increment of this index per plot
ce=0;                                        
% The following parameters will be assigned to keys.tt (please check in the previous section)
keys.ons(cc).tt.hands                      = [1 2];    
keys.ons(cc).tt.perturbations              = [0 1];    
keys.ons(cc).tt.choices                    = 0;        
keys.ons(cc).tt.combine_tuning_properties  ={'hand_tuning','in_IH_Facq_epoch_Ddre_han','in_CH_Facq_epoch_Ddre_han'};
keys.ons(cc).tt.selection                  ={};
keys.ons(cc).tt.unselect                   ={};
keys.ons(cc).tt.tasktypes                  ={'Msac_opt'};
% Same functionality as for population PSTHs
keys.ons(cc).group_parameter         = 'ungrouped';             % not sure if and how this works for response timing
keys.ons(cc).group_excluded          = 'ungrouped';
keys.ons(cc).conditions_to_plot      = {'Msac'};                % although this parameter is designed for multiple conditions, only one entry here is recommended
% because conditions are defined very specificly, a unique title needs to be defined
keys.ons(cc).comparisons_title       = 'Msac instructed space tuning';
ce=ce+1;                            % increment for each subplot
% conditions that are compared are defined separately for each subplot (NOT IDEAL, perturbation missing ?)
% IMPORTANT: if all conditions are identical, epoch modulation will be evaluated 
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1 1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[-1 1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
keys.ons(cc).comparisons_per_effector(ce).choice{2}=0;
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Cue', 0, 0.3};                                  % here we define the ordering of cells (by first significant bin in the time window around the event)
                                                                                                        % AND the window plotted in other figures (most cells tuned and tuning onset)
                                                                                                        % IMPORTANT: event here is defined by the alignment for the window corresponding to keys.WINDOWS_PER_TYPE, 
                                                                                                        % and the window can not exceed the respective window in keys.WINDOWS_PER_TYPE
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IS_IN/255; keys.colors.NH_CS_IN/255];  % colors for both conditions are defined here (first condition higher/second condition higher)
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='INI';                                         % epoch used as baseline for epoch tuning (enhancement/suppression)
keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Cue';                                       % subplot title
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
keys.ons(cc).comparisons_per_effector(ce).choice{2}=0;
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Saccade', -0.2, 0.7};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IS_IN/255; keys.colors.NH_CS_IN/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='MemL';
keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Saccade';
ce=ce+1;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{1}=0;
keys.ons(cc).comparisons_per_effector(ce).reach_hand{2}=0;
keys.ons(cc).comparisons_per_effector(ce).hemifield{1}=[-1];
keys.ons(cc).comparisons_per_effector(ce).hemifield{2}=[1];
keys.ons(cc).comparisons_per_effector(ce).choice{1}=0;
keys.ons(cc).comparisons_per_effector(ce).choice{2}=0;
keys.ons(cc).comparisons_per_effector(ce).order_onset={'Target', -0.2, 0.7};
keys.ons(cc).comparisons_per_effector(ce).colors=[keys.colors.NH_IS_IN/255; keys.colors.NH_CS_IN/255];
keys.ons(cc).comparisons_per_effector(ce).baseline_epoch='MemL';
keys.ons(cc).comparisons_per_effector(ce).title='Aligned to Target';

%% cell count keys
% overlapping tt and cc

cc=0;  % again, one increment of this index per plot

cc=cc+1;  % again, one increment of this index per plot
keys.ccs(cc).tt.hands               =[0 1 2]; 
keys.ccs(cc).tt.choices             =[0 1]; 
keys.sct(cc).tt.tasktypes           ={'Ddre_han'};
keys.ccs(cc).tt.plot_type           ='per_epoch'; %,'space_and_bilateral' 'space_x_hand' 'space_and_epoch' 'fixation_x_position' 'per_task'(try to avoid, because it slows down the processing)
keys.ccs(cc).tt.factor              ='space'; %  'hand' 'epoch' only for per_epoch and per_task relevant


keys.ccs(cc).conditions_to_plot  ={'Dcfr','Ddre','Ddsa'};
keys.ccs(cc).epochsE.Msac        ={'INI', 'Fhol','Cue','MemE','MemL','PreS','PeriS','TIhol','Thol'}';
keys.ccs(cc).epochsS.Msac        ={'Cue','MemE','MemL','PreS','PeriS','TIhol','Thol'}';
keys.ccs(cc).epochsE.Vsac        ={'INI', 'Facq','Cue','PreS','PeriS','Tacq','Thol'}';
keys.ccs(cc).epochsS.Vsac        ={'Cue','PreS','PeriS','Tacq','Thol'}';
keys.ccs(cc).epochsE.Ddre        ={'INI', 'Facq','Fhol','Cue','EDel','Del','PreS','PeriS','PostS','PreR','PeriR','PostR','Thol'}';
keys.ccs(cc).epochsS.Ddre        ={'Cue','Del','PeriS','PeriR','Thol'}';
keys.ccs(cc).epochsE.Dcfr        ={'INI', 'Facq','Fhol','Cue','EDel','Del','PreS','PeriS','PostS','PreR','PeriR','PostR','Thol'}';
keys.ccs(cc).epochsS.Dcfr        ={'Cue','Del','PeriS','PeriR','Thol'}';
keys.ccs(cc).epochsE.Ddsa        ={'INI', 'Facq','Fhol','Cue','EDel','Del','PreS','PeriS','PostS','PreR','PeriR','PostR','Thol'}';
keys.ccs(cc).epochsS.Ddsa        ={'Cue','Del','PeriS','PeriR','Thol'}';

%% scatter keys

cc=0;  % again, one increment of this index per plot

cc=cc+1;  
keys.sct(cc).tt.tasktypes={'Msac_opt','Vsac_opt'};
keys.sct(cc).X='in_NH_CS_CueG_epoch_EN_Vsac_opt';
keys.sct(cc).Y='in_NH_CS_CueG_epoch_EN_Msac_opt';
keys.sct(cc).X_sig='in_NH_CS_CueG_epoch_Vsac_opt';
keys.sct(cc).Y_sig='in_NH_CS_CueG_epoch_Msac_opt';
keys.sct(cc).color_option='category_as_color';
keys.sct(cc).categories={'visual_Msac_opt';'visuomotor_Msac_opt';'motor_Msac_opt'};



%% condition_identifiers - check if these keys are relevant still
keys.labels.handsLR={'NH','LH','RH'};
keys.labels.handsIC={'NH','IH','CH'};
keys.labels.choices={'in','ch'};


