function GGM = rd_gfc(gfcfile)
%{
    Read the.gfc file in the ICGEM standard format
Input:
    gfc file path
Retrun:
    GGM struct cnm-format
%}
arguments
    gfcfile
end

fid = fopen(gfcfile);
cont = true;

while(~feof(fid))
    s = fgetl(fid);
    if ~strncmpi(s, 'product_type', 12) && cont == true
        continue
    else
        cont = false;
    end

    if strncmpi(s, 'tide_system', 11)
        s_temp = isspace(s);
        s = s(find(s_temp > 0, 1) : end);
        s = strtrim(s);
        s_temp = isspace(s);
        if any(s_temp > 0)
            s = s(1 : find(s_temp == 1, 1));
        end
        tide_system = s;
    end

    if strncmpi(s, 'max_degree', 10)
        s_temp = isspace(s);
        s = s(find(s_temp > 0, 1) : end);
        s = strtrim(s);
        s_temp = isspace(s);
        if any(s_temp > 0)
            s = s(1 : find(s_temp == 1, 1));
        end
        max_degree = str2double(s);
    end

    if strncmpi(s, 'earth_gravity_constant', 22)
        s_temp = isspace(s);
        s = s(find(s_temp>0, 1) : end);
        s = strtrim(s);
        s_temp = isspace(s);
        if any(s_temp > 0)
            s = s(1 : find(s_temp == 1, 1));
        end
        GM = str2double(s);
    end

    if strncmpi(s, 'radius', 6)
        s_temp = isspace(s);
        s = s(find(s_temp > 0, 1) : end);
        s = strtrim(s);
        s_temp = isspace(s);
        if any(s_temp > 0)
            s = s(1 : find(s_temp == 1, 1));
        end
        R = str2double(s);
    end

    if strncmpi(s, 'norm', 4)
        s_temp = isspace(s);
        s = s(find(s_temp > 0, 1) : end);
        s = strtrim(s);  
        norm_type = s;
    end

    if strncmpi(s, 'key', 3)
        s_temp = textscan(s, repmat('%s', 1, 10));
        if isempty(s_temp{10})
            colm = 6;
        else
            colm = 7;
        end
    end

    if strncmpi(s, 'end_of_head', 11)
        break
    end
end
clear s s_temp
str = ['%s', repmat('%f', 1, colm)];
CS = textscan(fid, str);

if all(strcmp(CS{1}, 'gfc') == 0)
    error('Unsupported format of the GGM file.')
else
    ind = ismember(CS{1}, {'gfc', 'gfct'});
end
fclose(fid);

CS = CS(2:7);
CS = cell2mat(CS);
if exist('ind', 'var')
    CS = CS(ind, :);
end

[~, cols] = size(CS);
if cols < 4
    error('Unsupported format of the GGM file.')
end

GGM = struct();
GGM.R = R;
GGM.GM = GM;
GGM.max_degree = max_degree;
GGM.tide_system = tide_system;
GGM.norm_type = norm_type;
GGM.CS = CS;

end


