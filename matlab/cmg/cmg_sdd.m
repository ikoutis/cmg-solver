% cmg_sdd: wrapper for cmg_precondition
%
% supports backward combatibility with previous version
function pfun = cmg_sdd(A)

if (nargin>1)
    pfun = [];
    disp('You are invoking an older version of CMG. Check cmg_precondition.m')
    return
end

pfun = cmg_precondition(A);

