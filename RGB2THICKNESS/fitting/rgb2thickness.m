rgb = ones(256,256,256,'uint16');
file = load('rm_s.mat');
rgbmatrix = file.rgbmatrix;
for i = 1:256
    i
    for j = 1:256
        for k = 1:256
            rgb(i,j,k) = rgb2n2(i,j,k,rgbmatrix);
        end
    end
end
function [thick] = rgb2n2(y1,y2,y3,rgbmatrix)
    delta = (rgbmatrix(1,1)-y1)^2+(rgbmatrix(2,1)-y2)^2+(rgbmatrix(3,1)-y3)^2;
    res = 0;
    for i = 2:300
        delta0 = (rgbmatrix(1,i)-y1)^2+(rgbmatrix(2,i)-y2)^2+(rgbmatrix(3,i)-y3)^2;
        if delta0 < delta
            delta = delta0;
            res = i;
        end
    end
    thick = res;
end
save('colormap7.mat','rgb')