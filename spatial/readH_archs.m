function aout = readH_archs(fid,idm,jdm,varname);
% MCB, USM, 2024/4/6
% readH_archs reads variables varname at a certain layer LH in hycom
% assumes the file is already opened with identifier FID
% the variables are listed in the .b file at line LS
% idm and jdm are lengths of x and y axis
%
% list of variables ====================================
%
% montg1 
% srfhgt 
% u_btrop
% v_btrop
% u-vel. 
% v-vel. 
% temp   
% salin  

% grid --------------------------------------    
% plon
% plat
% qlon
% qlat
% ulon
% ulat
% vlon
% vlat
% pang
% pscx
% pscy
% qscx
% qscy
% uscx
% uscy
% vscx
% vscy
% cori
% pasp
%     
% topography --------------------------------------    
% depth



%% fields --------------------------------------
if     strcmp(varname,'montg1');   LS=1;  k=1;
elseif strcmp(varname,'srfhgt');   LS=2;  k=1;
%elseif strcmp(varname,'steric');   LS=3;  k=1;
elseif strcmp(varname,'u_btrop');  LS=3; k=1;
elseif strcmp(varname,'v_btrop');  LS=4; k=1;
elseif strcmp(varname,'u-vel');    LS=5; k=1;
elseif strcmp(varname,'v-vel');    LS=6; k=1;
elseif strcmp(varname,'temp');     LS=7; k=1;
elseif strcmp(varname,'salin');    LS=8; k=1;

   
% grid --------------------------------------    
elseif strcmp(varname,'plon');   LS=1;   k=1;
elseif strcmp(varname,'plat');   LS=2;   k=1;
elseif strcmp(varname,'qlon');   LS=3;   k=1;
elseif strcmp(varname,'qlat');   LS=4;   k=1;
elseif strcmp(varname,'ulon');   LS=5;   k=1;
elseif strcmp(varname,'ulat');   LS=6;   k=1;
elseif strcmp(varname,'vlon');   LS=7;   k=1;
elseif strcmp(varname,'vlat');   LS=8;   k=1;
elseif strcmp(varname,'pang');   LS=9;   k=1;
elseif strcmp(varname,'pscx');   LS=10;  k=1;
elseif strcmp(varname,'pscy');   LS=11;  k=1;
elseif strcmp(varname,'qscx');   LS=12;  k=1;
elseif strcmp(varname,'qscy');   LS=13;  k=1;
elseif strcmp(varname,'uscx');   LS=14;  k=1;
elseif strcmp(varname,'uscy');   LS=15;  k=1;
elseif strcmp(varname,'vscx');   LS=16;  k=1;
elseif strcmp(varname,'vscy');   LS=17;  k=1;
elseif strcmp(varname,'cori');   LS=18;  k=1;
elseif strcmp(varname,'pasp');   LS=19;  k=1;
    
% topography --------------------------------------    
elseif strcmp(varname,'depth');   LS=1;  k=1;

end 

%disp([num2str([LS k])])

%% extract data --------------------------------------
nt=idm*jdm;   
nta=floor((nt+4095)/4096)*4096;  %total length of line

% -- search from the beginning
% LS-1 is starting position
start_after_numbytes = ( (LS-1)+(k-1)*5 )*nta*4;
status = fseek(fid,start_after_numbytes,-1); % multiply times 4 bytes !
aout   = fread(fid,[1,    nt],'real*4');
aout   = reshape(aout,idm,jdm);

% --- so for now reshape works I need to transpose the matrix
aout=aout';
 
% --- remove land mask huge=2^100=1.2677*10^30 
aout(aout>1e30)=NaN; 

return

