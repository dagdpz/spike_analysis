function get_MA_STATES
global MA_STATES
MA_STATES.INI_TRI       = 1; % initialize trial
MA_STATES.FIX_ACQ       = 2; % fixation acquisition
MA_STATES.FIX_HOL       = 3; % fixation hold
MA_STATES.TAR_ACQ       = 4; % target acquisition
MA_STATES.TAR_HOL       = 5; % target hold
MA_STATES.CUE_ON        = 6; % cue on
MA_STATES.MEM_PER       = 7; % memory period
MA_STATES.DEL_PER       = 8; % delay period
MA_STATES.TAR_ACQ_INV   = 9; % target acquisition invisible
MA_STATES.TAR_HOL_INV   = 10; % target hold invisible
MA_STATES.MAT_ACQ       = 11; % target acquisition in sample to match
MA_STATES.MAT_HOL       = 12; % target acquisition in sample to match
MA_STATES.MAT_ACQ_MSK   = 13; % target acquisition in sample to match
MA_STATES.MAT_HOL_MSK   = 14; % target acquisition in sample to match
MA_STATES.SEN_RET       = 15; % return to sensors for poffenberger
MA_STATES.ABORT         = 19;
MA_STATES.SUCCESS       = 20;
MA_STATES.REWARD        = 21;
MA_STATES.ITI           = 50;
MA_STATES.SAC_INI       = 60;
MA_STATES.SAC_END       = 61;
MA_STATES.REA_INI       = 62;
MA_STATES.REA_END       = 63;
MA_STATES.TRI_END       = 90;
MA_STATES.ITI_END       = 98;
MA_STATES.CLOSE         = 99;
MA_STATES.SUCCESS_ABORT =-1;
end