clear
clc
lambda = 0.1:0.001:1;
B = 1.458;
C = 0.00354;
n = B + C./(lambda.^2);
plot(lambda,n)