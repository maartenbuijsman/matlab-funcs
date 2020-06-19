%% function ang_diff = angle2vectors_complex(v1,v2)
%% MCB,UCLA,19-03-2008
%% calcs angle between 2 vectors, which are complex
%% output ang_diff in radians
function ang_diff = angle2vectors_complex(v1,v2)

ang_diff = acos(real(dot(v1,v2))./(abs(v1).*abs(v2)));
