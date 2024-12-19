clc
clear
namelist = dir("thickness_all/*.png");
RGB = [];
for i = 1:length(namelist)
    name = [namelist(i).folder,'/',namelist(i).name]
    c = imread(name);
    R = c(:,:,1);
    G = c(:,:,2);
    B = c(:,:,3);
    thick_char = namelist(i).name(1:end-4);
    thick_str = string(thick_char);
    thick = double(thick_str);
    r = mean(R(:));
    g = mean(G(:));
    b = mean(B(:));
    RGB = [RGB;[thick,r,g,b]];
end
csvwrite("thickness_all/thickness2rgb.csv",RGB);
