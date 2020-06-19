%% function ang_diff = angle2vectors(v1,v2)
%% MCB,NIOZ,17-09-05
%% calcs angle between 2 vectors
%% output in radians
%% NOTE that individual vectors are located in collums of v1 and v2
function ang_diff = angle2vectors(v1,v2)

ang_diff = acos(dot(v1,v2)./(sqrt(v1(1,:).^2+v1(2,:).^2).*sqrt(v2(1,:).^2+v2(2,:).^2)));
