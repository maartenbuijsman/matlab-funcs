%% function findstr_infile
%% Maarten Buijsman, GFDL, 2011-4-29
%% searches file for a variable and extracts the numeric value
%% between identifyers id1 and id2

function  value = findstr_infile(filename,variable,id1,id2)


%% test ----------------------------------
%id1='=';id2=',';
%filename = [dirin 'data']
%variable = 'deltaT'; 
%variable = 'select_rStar=2';
%% test ----------------------------------

fid = fopen(filename);
value = [];
while feof(fid)~=1
    tline = fgetl(fid);
    Istr1 = findstr(variable,tline);
    if ~isempty(Istr1)
        I1 = findstr(id1,tline); 
        I2 = findstr(id2,tline);
        if isempty(I2); %if no comma go to the end of the line
            I2 = length(tline); 
            value = str2num(tline(I1+1:I2)); disp([variable ' = ' num2str(value)]);             
        else
            value = str2num(tline(I1+1:I2-1)); disp([variable ' = ' num2str(value)]); 
        end
    end
end
fclose(fid);

return