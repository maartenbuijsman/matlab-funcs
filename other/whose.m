%% whose
%% MCB, 2010-03-16
%% lists total anmount of bytes

SSS = whos;
mem=0;
for i=1:length(SSS)
    mem=mem+SSS(i).bytes;
end
disp(['total bytes in memory is ' num2str(mem) ' B (' num2str(mem/1e6) ' MB)'])



