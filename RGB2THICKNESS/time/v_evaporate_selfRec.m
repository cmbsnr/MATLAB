clear
clc
load("rm_f_2000.mat");
% rgbmatrix(1,:) = rgbmatrix(1,:) * 0.85;
index = 0;     % 
lim   = 50;      % 扫描范围 ±20 nm
% 输入视频
videoFile = 'Sout.mp4';

% 创建 VideoReader 对象
v = VideoReader(videoFile);

% 预分配（提高性能）
numFrames = floor(v.Duration * v.FrameRate);

frames = cell(numFrames, 1);   % 存图像
timeList = zeros(numFrames, 1); % 存时间

i = 1;

% 逐帧读取
while hasFrame(v)
    frame = readFrame(v);      % 读取当前帧
    
    frames{i} = frame;         % 存入cell
    timeList(i) = v.CurrentTime; % 当前时间（单位：秒）
    
    i = i + 1;
end
[H, W, ~] = size(frame);

h = H / 16;
w = W / 16;

y1 = 160; % floor(H/2 - h/2) + 1;
y2 = 1840; % floor(H/2 + h/2);

x1 = 1800; % floor(W/2 - w/2) + 1;
x2 = 2040; % floor(W/2 + w/2);
% 如果实际帧数小于预估，截断
frames = frames(1:i-1);
timeList = timeList(1:i-1);

N = numel(frames);
thick = ones([N,1]);
rgb = ones([N,3]);
for i = N:-1:1
    frame = frames{i};
    t = timeList(i);
    frame_center = frame(y1:y2, x1:x2, :);
    rgb_mean = squeeze(mean(frame_center, [1 2]));

    % 边界保护
    left  = max(1, index);
    right = min(size(rgbmatrix,2), index + 2 * lim);
    
    % 截取子区间
    sub_rgb = rgbmatrix(:, left:right);
    
    % 计算距离
    diff = sub_rgb - rgb_mean;
    dist = sum(diff.^2, 1);
    
    % 找局部最优
    [~, local_idx] = min(dist);
    
    % 转换回全局索引
    idx = left + local_idx - 1;
    rgb(i,:) = rgb_mean;
    thick(i) = idx;   % 单位：nm
    index = idx;
end
plot(timeList,thick,'o')
xlabel('Time(s)')
ylabel('Film thickness(nm)')
title('Evaporate Viscosity')