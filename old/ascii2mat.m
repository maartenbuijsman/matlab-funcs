%% program ascii2mat.m
%% Maarten Buijsman, NIOZ, 02-11-2004
%% ascii2mat.m converts any ascii to mat files
%% uses dirin, namein, ext

clear all

%% for 2004 bottom
%year = 2004; period = [1:9 11 22:33]; 
%namein = 'bot_n'; ext = 'out'; dir = ['C:\data\TESO_data\getdata\p',int2str(year),'\bottom'];

%% for GPS_sc60nav
% year = 2004; period = [21:33]; 
% namein = 'sc60nav_n'; ext = 'out'; dir = ['C:\data\TESO_data\getdata\p',int2str(year),'\GPS_sc60nav'];

%% for GPS_sc60nav
% year = 2004; period = [31]; 
% namein = 'GPS_DH_n'; ext = 'out'; dir = ['C:\data\TESO_data\getdata\p',int2str(year),'\GPS'];
% %namein = 'GPS_TX_n'; ext = 'out'; dir = ['C:\data\TESO_data\getdata\p',int2str(year),'\GPS'];

%% for sbe21d
year = 2004; period = [2:9 11:26 30:33]; 
namein = 'sbe21_n'; ext = 'out'; dir = ['C:\data\TESO_data\getdata\p',int2str(year),'\sbe21d'];

for ii = 1:length(period);
    step = period(ii);
    
    %%-- load data    
    disp(['start loading ADCP ',int2str(step)])
    hh = load ([dir,'\',namein,int2str(step),'.',ext]);    
    disp('finished loading')
   
    eval(['save ',dir,'\',namein,int2str(step),'.mat hh'])
end



