function [outVar] = var_ave_quick(LON,LAT,lon,lat,gdepth,ddeg,inVar,inlon,inlat,inDepth,inAAA);

%% [outVar] = var_ave_quick(LON,LAT,lon,lat,gdepth,ddeg,inVar,inlon,inlat,inDepth,inAAA);
%  MCB, USM, 2020-3-25
% call from main program

% resolution of ingrid
jmed  = median(1:size(inlon,1));
inRes = inlon(jmed,2)-inlon(jmed,1);

%number of buffer points used
nump = ceil(ddeg/inRes)*2;   

%% make the matrices larger before interpolating
%  for gridding purposes

% left side
npx1 = ceil((inlon(jmed,1) - lon(1))/inRes);

% right side
npx2 = ceil((lon(end) - inlon(jmed,end))/inRes);

inlon   = [  inlon(:,end-npx1+1:end)-360  inlon     inlon(:,1:npx2)+360];
inlat   = [  inlat(:,end-npx1+1:end)      inlat     inlat(:,1:npx2)];
inVar   = [  inVar(:,end-npx1+1:end)      inVar     inVar(:,1:npx2)]; 
inDepth = [inDepth(:,end-npx1+1:end)      inDepth inDepth(:,1:npx2)]; 
inAAA   = [  inAAA(:,end-npx1+1:end)      inAAA     inAAA(:,1:npx2)]; 

% check
inlon(jmed,[1:npx1+5])  
inlon(jmed,[end-npx2-5:end])

[ny,nx] = size(inlon);

disp(['dimensions ny and nx of expanded grid are ' num2str([ny,nx])]);


%% loop over the new grid - outer loop is over lon
outVar  = ones(size(LON))*NaN;

lon2 = LON(1,:);
lat2 = LAT(:,1)';

for i=1:length(lon)-1
    
    % output to screen every 10 step
    if rem(i,10)==0
        disp(['i=' num2str(i) ' of ' num2str(length(lon2))]);
    end
    
    % polygon of the lon coordinates of the grid
    lonp = [lon(i) lon(i+1) lon(i+1) lon(i)   lon(i)];
    
    % ilat is the index of plat
    % it is advanced at the end of the loop
    ilat = 1;
    
    % loop over lat for constant lon
    for j=1:length(lat)-1
        
        %if ~isnan(gflux(j,i))
        if gdepth(j,i)>0 %then ocean; this retains more of the data           
            % find plon(ilat,:) nearest to center of cell
            % cell center is ( lon(i) + lon(i+1) )/2
            % Iselx is index of longitude
            [d,Iselx] = min(abs(inlon(ilat,:) - (lon(i) + lon(i+1))/2));
            
            % find latitude corners of cell polygon
            latp = [lat(j) lat(j)   lat(j+1) lat(j+1) lat(j)];
            
            % speed up search
            % find plat(:,Iselx) nearest to center of cell
            % (lat(j) + lat(j+1))/2 is cell center
            % thus  Isely and Iselx are the indices on plon and plat closest to the cell center
            [d,Isely] = min(abs( inlat(:,Iselx) - (lat(j) + lat(j+1))/2 ));
            
            % select all indices in larger box marked by Isely and Iselx  +/- nump points
            II = max([Iselx-nump 1]):min([Iselx+nump nx]);
            JJ = max([Isely-nump 1]):min([Isely+nump ny]);
            %II,JJ, size(AAA);
            % find indices inside polygon, which marks the cell boundaries
            % isel is used for averaging
            isel = inpolygon(inlon(JJ,II),inlat(JJ,II),lonp,latp);
            
            % exlclude land and NaNs
            IN = find(isel==1 & inDepth(JJ,II)>0 & ~isnan(inVar(JJ,II)));
            numpa(j,i) = length(IN);
            
            % do averaging
            if ~isempty(IN)
                varTOTa2 = inVar(JJ,II);
                AAA2     = inAAA(JJ,II);
                % area averaging
                outVar(j,i) = sum(varTOTa2(IN).*AAA2(IN))/sum(AAA2(IN));
            end
        end
        
        %% DOES THIS NEED TO BE INSIDE  if gdepth(j,i)>0 ????
        % advance ilat
        ilat = min([JJ(end) ny]);
    end
end