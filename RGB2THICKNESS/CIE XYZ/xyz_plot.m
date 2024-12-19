clc
clear
x = 1000;       %薄膜厚度/nm
n = 1.428;      %薄膜折射率
n0 = 1;         %空气折射率
ng = 3.836;     %硅片折射率
i = 0;          %折射角
lambda = 360:830;
x = 1:x;
delta = 4*pi*n*x'*cos(i)./lambda;
R = (n^2*(n0-ng)^2*cos(delta/2).^2+(n^2-n0*ng)^2*sin(delta/2).^2)./...
       (n^2*(n0+ng)^2*cos(delta/2).^2+(n^2+n0*ng)^2*sin(delta/2).^2);
cie = colorMatchFcn('1931_full');
XYZ = R*cie(:,2:4);
xyz_norm = normalize(XYZ,2,'norm',1);
rgb = xyz2rgb(xyz_norm);
subplot(2,3,1)
plot(xyz_norm(:,1))
title('x')
subplot(2,3,2)
plot(xyz_norm(:,2))
title('y')
subplot(2,3,3)
plot(xyz_norm(:,3))
title('z')
subplot(2,3,4)
plot(rgb(:,1))
title('r')
subplot(2,3,5)
plot(rgb(:,2))
title('g')
subplot(2,3,6)
plot(rgb(:,3))
text(100,100,'100')
title('b')