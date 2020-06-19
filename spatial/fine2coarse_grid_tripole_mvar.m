function [outVar] = fine2coarse_grid_tripole_mvar(lon,lat,LON,ddeg,gdepth,inVar,inlon,inlat,inDepth,inAAA,latlim);

%% [outVar] = fine2coarse_grid_tripole_mvar(lon,lat,LON,ddeg,gdepth,inVar,inlon,inlat,inDepth,inAAA,latlim);
%  Maarten Buijsman, USM, 2019-6-6
%  average fine grid to coarse grid
%  do this for multiple 2D variables stacked into one 3D array inVar
%
%  lon and lat are vectors and grid faces of coarse grid
%  LON is the cell center of lon, but is an array; is used to declare the size of the output variables
%  gdepth is the depth of the coarse grid


% resolution of ingrid
inRes = inlon(500,2)-inlon(500,1);
nump = ceil(ddeg/inRes)*2;   %number of buffer points used
[ny,nx] = size(inlat);

nzz = size(inVar,3);
outVar = ones([size(LON) nzz])*NaN;
size(outVar)

for i=1:length(lon)-1
    
    % output to screen every 10 step
    %if rem(i,10)==0
        disp(['i=' num2str(i) ' of ' num2str(length(lon)-1)]);
    %end
    
    % polygon of the lon coordinates of the grid
    lonp = [lon(i) lon(i+1) lon(i+1) lon(i)   lon(i)];
    
    % ilat is the index of plat
    % it is advanced at the end of the loop
    ilat = 1;
    
    % only select north of latlim
    Iall = find(inlon > min(lonp)-ddeg & inlon < min(lonp)+ddeg & inlat > latlim-ddeg);
    
    % loop over lat for constant lon
    for j=1:length(lat)-1
        
        %if ~isnan(gflux(j,i))
        if gdepth(j,i)>0 %then ocean; this retains more of the data
            % find plon(ilat,:) nearest to center of cell
            % cell center is ( lon(i) + lon(i+1) )/2
            % Iselx is index of longitude
            if lat(j)< latlim
            
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
                IN = find(isel==1 & inDepth(JJ,II)>0 & ~isnan(inVar(JJ,II,1)));
                numpa(j,i) = length(IN);
                
                % do averaging
                if ~isempty(IN)

                    AAA2 = inAAA(JJ,II);
                    
                    for ii = 1:nzz
                        varTOTa2 = inVar(JJ,II,ii);
                        
                        % area averaging
                        outVar(j,i,ii) = sum(varTOTa2(IN).*AAA2(IN))/sum(AAA2(IN));
                    end
                end
                
                
            else % in tripole grid area above latlim
                
                % find latitude corners of cell polygon
                latp = [lat(j) lat(j)   lat(j+1) lat(j+1) lat(j)];

                % find indices inside polygon, which marks the cell boundaries
                % isel is used for averaging
                isel = inpolygon(inlon(Iall),inlat(Iall),lonp,latp);

                % exlclude land and NaNs
                varTOTa0 = inVar(:,:,1);
                IN = find(isel==1 & inDepth(Iall)>0 & ~isnan(varTOTa0(Iall)));
                numpa(j,i) = length(IN);
                
                % do averaging
                if ~isempty(IN)
                    AAA2 = inAAA(Iall);
                    
                    for ii = 1:nzz
                        varTOTa1 = inVar(:,:,ii);
                        varTOTa2 = varTOTa1(Iall);
                        
                        % area averaging
                        outVar(j,i,ii)  = sum(varTOTa2(IN).*AAA2(IN))/sum(AAA2(IN));
                    end

                end                
            
            end
            
%             % test
%             plot(lonp,latp,'r-')
%             hold on
%             plot(inlon(JJ,II),inlat(JJ,II),'k.')
%             hold off
%             pause
            
            % advance ilat
            % is this the problem????
            ilat = min([JJ(end) ny]);
            
        end
        
        %                 % advance ilat
        %                 ilat = min([JJ(end) ny]);
    end
end

% ------------------------------------------
