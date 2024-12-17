clc
clear
x = 1000;
[l, xFcn, yFcn, zFcn] = colorMatchFcn('1931_full');
rgbmatrix = zeros(3,x);
for x = 1:x
    X = 0;
    Y = 0;
    Z = 0;
    for lambda = 360:830
        index = (l == lambda);
        R = reflect(lambda,x);
        X = X + R *xFcn(index);
        Y = Y + R * yFcn(index);
        Z = Z + R * zFcn(index);
    end
    S = sum([X,Y,Z]);
    X = X/S;
    Y = Y/S;
    Z = Z/S;
    XYZ = xyz2rgb([X,Y,Z],"WhitePoint","d55");
    rgbmatrix(1,x) = XYZ(1);
    rgbmatrix(2,x) = XYZ(2);
    rgbmatrix(3,x) = XYZ(3);
end
r = rgbmatrix(1,:);
g = rgbmatrix(2,:);
b = rgbmatrix(3,:);
C(:,:,1) = r.*ones(100,1);
C(:,:,2) = g.*ones(100,1);
C(:,:,3) = b.*ones(100,1);
image(C)