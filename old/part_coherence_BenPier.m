%% function [gam2s,phis,Gyy2s,Gs,noise,fk,Nb] = part_coherence_BenPier(tgrid,xm,Kb,conf,rem_ave,win);
%% Maarten Buijsman, NIOZ, 07-06-06
%% partial coherence according to Bendat and Piersol (1986), examp. on p.419
%% performs size(xm,1)-1 steps, in which partial coherence of size(xm,1)-1 independent vars
%% with dependent var is calculated, in each step the coherence of the previous var on the current 
%% vars is subtracted, so that only the partial coherence gam2s with the dependent var remains
%% INPUT: time vector tgrid, matrix xm with indepentdent vars and dependent var in last line; 
%% NOTE that independent vars in xm should have declining energy,
%% number of blocks Kb, confidence level in [%], rem_ave=0 (keep trend) or 1 (remove trend),
%% win = 1 (yes) or 0 (no) Hanning window
%% OUTPUT: struct array gam2s containing part. coher. per frequecy fk of size xm with on 
%% 1st line ordinary and in last array multiple coherence, struct array phis with the phase in [rad]
%% size matches part. coher. arrays, struct array Gyy2s with autospectra of all ind. vars., 
%% struct array Gs of length(gam2s)^3 containing all cross-correlation spectra for all steps+1
%% vector noise is the remaining unexplained coherence per fk, number of 50% blocks really used Nb 
%% NOTE that positive phi means that dep.var.(e.g.Q) leads indep.var. (e.g.tau) 
%% and negative phi means that indep.var. leads dep.var. (see test in coherence_BenPier)

function [gam2s,phis,Gyy2s,Gs,noise,fk,Nb] = part_coherence_BenPier(tgrid,xm,Kb,conf,rem_ave,win);

%% make auto and cross-correlation matrices
%% store in structure
ms = size(xm,1);
for ii=1:ms; for hh=1:ms; for jj=1:ms; Gs(hh,jj,ii).sp = 0; end; end; end;
       
for hh=1:ms     %% vertical    
    disp(['row ',num2str(hh)]);
    for jj=1:ms %% horizontal
        if jj<hh %% do the mirror thing
            Gs(hh,jj,1).sp = conj(Gs(jj,hh,1).sp);                
        else
            [Gs(hh,jj,1).sp,fk,phase,Nb] = autocross_BenPier(tgrid,xm(hh,:),xm(jj,:),Kb,rem_ave,win);
        end
    end
end

%test = []; for hh=1:ms; for jj=1:ms; test(hh,jj) = Gs(hh,jj,1).sp(5); end; end

%% loop over steps in process
Ls = []; gam2s = []; Gyy2s = []; phis = [];
for ii=1:ms-1 
    %% fill Lm, use Gm as a blue print
    for jj=ii+1:ms
        Ls(ii,jj).sp = Gs(ii,jj,ii).sp./Gs(ii,ii,ii).sp;
    end
    
    %% fill Gm using Lm, this has to be done for every frequency
    for hh=ii+1:ms     %% vertical
        for jj=ii+1:ms %% horizontal
            Gs(hh,jj,ii+1).sp = Gs(hh,jj,ii).sp-Ls(ii,jj).sp.*Gs(hh,ii,ii).sp;
            if jj<hh %% do the mirror thing
                %% Gm(hh,jj,ii+1) = -Gm(jj,hh,ii+1);
                Gs(hh,jj,ii+1).sp = conj(Gs(jj,hh,ii+1).sp);                
            end
        end
    end
    
    %% ordinary and partial coherence functions
    gam2s(ii).sp = real((abs(Gs(ii,ms,ii).sp).^2) ./ ( Gs(ii,ii,ii).sp.*Gs(ms,ms,ii).sp ));
    
    %% get phases from the crossterms; NOTE that only the partial coherence has a phase
    phis(ii).sp = atan2(-imag(Gs(ii,ms,ii).sp),real(Gs(ii,ms,ii).sp));
    
    %% 2X check
    Gyy2s(ii).sp = real(gam2s(ii).sp.*Gs(ms,ms,ii).sp);
end
%% multiple coherence function
gam2s(ms).sp = 1-real(Gs(ms,ms,ms).sp./Gs(ms,ms,1).sp);
noise = real(Gs(ms,ms,ms).sp./Gs(ms,ms,1).sp);

%% check, last term is the noise
%% Gyy2  = Gs(ms,ms,1).sp; (OK, with complex Gm matrix)
Gyy2sum = 0; for hh=1:length(Gyy2s); Gyy2sum = Gyy2sum + Gyy2s(hh).sp; end
Gyy2 = Gyy2sum + real(Gs(ms,ms,ii+1).sp);
%figure; loglog(fk,Gs(ms,ms,1).sp,'b-',fk,Gyy2,'r--'); %% check

disp(['max. duration is ',num2str(1/fk(2)),' days; number of 50% blocks is ',num2str(Nb)]);

