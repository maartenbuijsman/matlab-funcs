function AA = area_grid(LON,LAT);

%% AA = area_grid(LON,LAT);
%  MCB, USM, 2019-18-07
%  computes area given rectangular matrices of LON and LAT
%  LON and LAT need to be equidistant in degrees: dy x dy
%  LON and LAT need to be MxN matrices with M=latitude and N=longitide
%  Output AA is the same size as LON

% why is there an anomaly close to 0, 180 or 360 degrees latitide???
dy = spheric_dist(LAT(2,1),LAT(1,1),LON(2,1),LON(1,1));
dx = spheric_dist(LAT(:,2),LAT(:,1),LON(:,2),LON(:,1));

aar = dx*dy;
AA = repmat(aar,[1 size(LON,2)]);    


return

figure
mypcolor(LON,LAT,AA); shading flat





%%  what is below shows that there is a singularity near 0 180 360 degrees long and 0 lat
dy  = LAT(2,1) - LAT(1,1);
lat = LAT(:,1); 
E = wgs84Ellipsoid;

aar = ones(size(LON,1),1);
lonc = 180;
for j=1:length(lat)
    lonp = [lonc-dy/2   lonc+dy/2   lonc+dy/2   lonc-dy/2   lonc-dy/2];
    %lonp = [       0           dy          dy          dy          0];
    latp = [lat(j)-dy/2 lat(j)-dy/2 lat(j)+dy/2 lat(j)+dy/2 lat(j)-dy/2];                 
    aar(j) = areaint(latp,lonp,E);
end
AA = repmat(aar,[1 size(LON,2)]);    


% figure
% plot(lonp,latp)
% 
% figure
% mypcolor(LON,LAT,AA); shading flat

% figure
% plot(lat,aar)
% holder

%% more variations of the same

dy  = LAT(2,1) - LAT(1,1);
E = wgs84Ellipsoid;

AA = ones(size(LON));
for i=1:size(LON,2)
    disp(num2str(i))
    for j=1:size(LON,1)
        lonc = LON(j,i);
        latc = LAT(j,i);
        lonp = [lonc-dy/2   lonc+dy/2   lonc+dy/2   lonc-dy/2   lonc-dy/2];
        latp = [latc-dy/2   latc-dy/2   latc+dy/2   latc+dy/2   latc-dy/2];
        AA(j,i) = areaint(latp,lonp,E);
    end
end

figure
mypcolor(LON,LAT,AA); shading flat
colorbar

figure
mypcolor(LON,LAT,AA-AA2); shading flat
colorbar


%% submitted to mathworks

dy  = 0.5;
E = wgs84Ellipsoid;

%lonc = -5:1:5;
lonc = [175:1:185]-dy/2;
latc = [-5:1:5]-dy/2;

AA = ones(length(latc),length(lonc));
for i=1:length(lonc)
    disp(num2str(i))
    for j=1:length(latc)
        lonp = [lonc(i)-dy/2   lonc(i)+dy/2   lonc(i)+dy/2  lonc(i)-dy/2  lonc(i)-dy/2];
        latp = [latc(j)-dy/2   latc(j)-dy/2   latc(j)+dy/2  latc(j)+dy/2  latc(j)-dy/2];
        AA(j,i) = areaint(latp,lonp,E);
    end
end

figure
pcolor(lonc,latc,AA); 
colorbar













