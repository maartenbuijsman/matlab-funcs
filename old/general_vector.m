%% general_vector
%% play around with gradients

clear all

r=[2 2 0];
f=[3 3 0];

dot(r,f)

fabs = sqrt(sum(f.^2));
rabs = sqrt(sum(r.^2));

alp = acos(dot(r,f)/(fabs*rabs))*180/pi

fabs*rabs

ds = sqrt(2^2+2^2+2^2)