%% function  value = findstr_infile_txt(filename,variable,id1,id2,char)
%% Maarten Buijsman, GFDL, 2011-6-22
%% based on findstr_infile
%% searches file for a variable and extracts the numeric value or text
%% between identifyers id1 and id2
%% char = 'txt' ot 'num'

function  value = findstr_infile_txt(filename,variable,id1,id2,char)


%% test ----------------------------------
%id1='=';id2='';
%filename = [dirin 'data']
%variable = 'diffKhT'; 
%% test ----------------------------------

fid = fopen(filename);
value = NaN;
while feof(fid)~=1
    tline = fgetl(fid);
    Istr1 = findstr(variable,tline);
    if ~isempty(Istr1)
        I1 = findstr(id1,tline); 
        I2 = findstr(id2,tline);
        if isempty(findstr('#',tline)) %% omit commented out lines
            if isempty(I2); %if no comma go to the end of the line
                I2 = length(tline); 
                if strcmp(char,'num')
                    value = str2num(tline(I1+1:I2)); disp([variable ' = ' num2str(value)]);             
                elseif  strcmp(char,'txt')
                    value = tline(I1+1:I2); disp([variable ' = ' value]);
                end
            else
                if strcmp(char,'num')
                    value = str2num(tline(I1+1:I2-1)); disp([variable ' = ' num2str(value)]); 
                elseif  strcmp(char,'txt')
                    value = tline(I1+1:I2-1); disp([variable ' = ' value]);                 
                end
            end
        end
    end
end
fclose(fid);

return