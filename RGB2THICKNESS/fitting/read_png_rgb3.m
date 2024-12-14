a = csvread("thick2.csv");
namelist = dir("thickness2/*.png");
RGB = [];
for i = 1:length(namelist)
    name = [namelist(i).folder,'/',namelist(i).name]
    b = imread(name);
    R = b(:,:,1);
    G = b(:,:,2);
    B = b(:,:,3);
    r = mean(R(:));
    g = mean(G(:));
    b = mean(B(:));
    RGB = [RGB;[r,g,b]];
end
a = [a,RGB];
csvwrite("thickness2rgb2.csv",a);