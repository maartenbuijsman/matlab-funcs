%% function [X1,X2,Xh1,Xh2,fk,y1,y2,t,N,Ns,dt] = Xspectra_BenPier(t,y1,y2,Kb,rem_ave);
%% Maarten Buijsman, NIOZ, 29-05-06
%% calculates spectra for two variables according to Bendat and Piersol (1986)
%% INPUT: vectors t (time) and y1 and y2 (data), Kb number of blocks, and rem_ave to remove mean (not=0, yes=1)
%% OUTPUT: smoothed X1 and X2 and Hanning windowed Xh1 and Xh2 spectra as a function of frequency fk
%% and data modified data y1, y2, t, N, Ns, and dt 
%% NOTE; 1) rem_ave=0 and win=0; rem_ave=1 and win=0; rem_ave=1 and win=1 => NO  leftside lobe
%%       2) rem_ave=0 and win=1;                                          => YES leftside lobe
%% leftside lobe = very large X value in first frequency (why???)

function [X1,X2,Xh1,Xh2,fk,y1,y2,t,N,Ns,dt] = Xspectra_BenPier(t,y1,y2,Kb,rem_ave);

%% test data
% rem_ave = 1; Kb = 5;
% t = [0:1:12*40]; y1 = sin(2*pi/12*tg)+2; y2 = y1; %figure; plot(tg,xg)

%% length data and delta t
N = length(t);
dt = t(2)-t(1);

%% make power spectra using hamming windows and (50%) block averaging
%% 1) make sure that Ns is even
%Kb = 3;           % non-overlapping blocks
Ns = N/Kb;
Ndum = N; % start value
while rem(Ns,2)>0
    Ns = Ndum/Kb; % length block
    if rem(Ns,2)>0; Ndum = Ndum-1; end
end

%% 2) remove data points
if Ndum<N
    y1 = y1(1:Ndum); % shorten time series
    y2 = y2(1:Ndum); 
    t  = t(1:Ndum);
    N = Ndum;
end

%% 3) with 50% block avaraging determine the indices of the blocks 
%% loop over N/(Ns/2)-1 blocks
istart = 1; k1 = 0; Ib = [];
while istart-1+Ns<=N
    k1 = k1 + 1;
    Ib(k1,:) = [istart:istart-1+Ns]; %% cont. indices of data
    
    istart = istart+Ns/2;
end

%% 4) loop over blocks, and first determine the 'autospectral density function'
X1 = []; X2 = []; Xh1 = []; Xh2 = [];
for II=1:size(Ib,1)
    %% remove mean in each block
    if rem_ave == 1
        %% subtract simple mean
%         yd1 = y1(Ib(II,:)) - mean(y1(Ib(II,:)));
%         yd2 = y2(Ib(II,:)) - mean(y2(Ib(II,:)));
        %% or subtract linear fit
        [x_s,y_fit,r,p,cf] = line_fit(Ib(II,:),y1(Ib(II,:)),1); yd1 = y1(Ib(II,:))-(Ib(II,:)*cf(1)+cf(2));        
        [x_s,y_fit,r,p,cf] = line_fit(Ib(II,:),y2(Ib(II,:)),1); yd2 = y2(Ib(II,:))-(Ib(II,:)*cf(1)+cf(2));                
    else
        yd1 = y1(Ib(II,:));
        yd2 = y2(Ib(II,:));
    end
    %mean(yd1),mean(yd2)
    
    %% fourier spectrum y1
    [X1(II,:),Xh1(II,:),fk] = fft_BenPier(t(Ib(II,:)),yd1);
    
    if mean(yd1) ~= mean(yd2) %% to save time
        %% fourier spectrum y2
        [X2(II,:),Xh2(II,:),fk] = fft_BenPier(t(Ib(II,:)),yd2);
    else;  X2(II,:) = X1(II,:); Xh2(II,:) = Xh1(II,:);
    end
end
