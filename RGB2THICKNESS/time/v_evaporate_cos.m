clear
clc
load("rm_f.mat");
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

h = H / 8;
w = W / 8;

y1 = floor(H/2 - h/2) + 1;
y2 = floor(H/2 + h/2);

x1 = floor(W/2 - w/2) + 1;
x2 = floor(W/2 + w/2);
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
    rgb_mean = squeeze(mean(frame, [1 2]));
    cos_sim = (rgbmatrix' * rgb_mean) ./ ...
          (vecnorm(rgbmatrix)' * norm(rgb_mean));
    [~, idx] = max(cos_sim);
    rgb(i,:) = rgb_mean;
    thick(i) = idx;   % 单位：nm
    index = idx;
end
plot(thick)