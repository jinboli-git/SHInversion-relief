function [datac, re] = gridclip(Data, dlon, dlat, clip)
% a fuction to clip the griddata
% input
% Data  is the 2-D data matrix
% dlon and dlat is the resultion of lon ,lat
% clip is the vector [upclip,downclip,leftclip,rigthclip] (degree)
% output
% data grid after clip

arguments
    Data
    dlon
    dlat
    clip = 0
end
if isscalar(clip)
    if clip==0
        datac = Data;
        return
    else
        clip = repmat(clip,4,1);
    end
end

if isempty(dlat)
    error('please input the dlon and dlat');
else
%     up = floor(clip(1)/dlat);
%     down = floor(clip(2)/dlat);
    up = round(clip(1)/dlat);
    down = round(clip(2)/dlat);
end

if isempty(dlon)
    left = 0;
    right = 0;
else
    left = round(clip(3)/dlon);
    right = round(clip(4)/dlon);
end

datac = Data(up+1:end-down,left+1:end-right);
re = [up, down, left, right];
end
