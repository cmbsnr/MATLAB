clc
clear
test = imread("test.png");
siz = size(test);
test_xyz = rgb2xyz(test);
r = ones(siz,"uint8");
r(:,:,1) = test(:,:,1);
g = ones(siz,"uint8");
g(:,:,2) = test(:,:,3);
b = ones(siz,"uint8");
b(:,:,3) = test(:,:,3);
subplot(4,4,1)
image(test)
subplot(4,4,2)
image(r)
title('r')
subplot(4,4,3)
image(g)
title('g')
subplot(4,4,4)
image(b)
title('_b')
subplot(4,4,5)
image(test(:,:,1))
title('r')
subplot(4,4,9)
image(test(:,:,2))
title('g')
subplot(4,4,13)
image(test(:,:,3))
title('b')