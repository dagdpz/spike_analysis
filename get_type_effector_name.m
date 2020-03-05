function [type_effector type_effector_short typ eff]=get_type_effector_name(type,effector)
switch type
    case 1
        typ='Fixation';
        typs='F';
    case 2
        typ='Visually_guided';
        typs='V';
    case 3
        typ='Memory';
        typs='M';
    case 4
        typ='Delay';
        typs='D';
    case 5
        typ='M2S';
        typs='M2';
    case 6
        typ='S2S';
        typs='S2';
end

switch effector
    case 0
        eff='saccades';
        effs='sac';
    case 1
        eff='free_gaze_reaches';
        effs='fgr';
    case 2
        eff='joint';
        effs='joi';
    case 3
        eff='diss_saccades';
        effs='dsa';
    case 4
        eff='diss_reaches';
        effs='dre';
    case 5
        eff='sen_reaches';
        effs='sen';
    case 6
        eff='central_fix_reaches';
        effs='cfr';
end
type_effector=[typ '_' eff];
type_effector_short=[typs effs];
end