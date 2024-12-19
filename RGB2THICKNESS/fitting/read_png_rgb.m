a = csvread("consistence2thickness.csv");
namelist = dir("thickness1/*.png");
RGB = [];
for i = 1:length(namelist)
    name = [namelist(i).folder,'/',namelist(i).name]
    c = imread(name);
    R = c(:,:,1);
    G = c(:,:,2);
    B = c(:,:,3);
    r = mean(R(:));
    g = mean(G(:));
    b = mean(B(:));
    RGB = [RGB;[r,g,b]];
end
a = [a,RGB];
csvwrite("thickness2rgb2.csv",a);
