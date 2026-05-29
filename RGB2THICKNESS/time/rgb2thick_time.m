clear
clc
lim = 50; %扫描范围
index = 10;
load("rm_f.mat");
rgbmatrix(1,:) = rgbmatrix(1,:) * 0.85;
% 输入视频
videoFile = 'Sout.mp4';

% 创建 VideoReader 对象
v = VideoReader(videoFile);
frameRate = v.FrameRate;
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

y1 = floor(H/2 - h/2) + 1;
y2 = floor(H/2 + h/2);

x1 = floor(W/2 - w/2) + 1;
x2 = floor(W/2 + w/2);
% 如果实际帧数小于预估，截断
frames = frames(1:i-1);
timeList = timeList(1:i-1);

N = numel(frames);
thick = cell([N,1]);
rgb = ones([N,3]);
thick_avg = ones([N,1]);

for k = N:-1:1
    frame = frames{k};
    t = timeList(k);
    test = frames{k};
    test = imrotate(test,-90);
    
    siz = size(test);
    mat = nan(siz(1), siz(2));
    for i = 1:siz(1)
        for j = 1:siz(2)
            r = double(test(i,j,1));
            g = double(test(i,j,2));
            b = double(test(i,j,3));
            mid = rect_mean_exclude_nan(mat,i,j,index);
            if k < N
                mid2 = rect_mean_exclude_nan(thick{k+1},i,j,index);
                mid = (mid + mid2)/2;
            end
            mat(i,j) = rgb2Thick(rgbmatrix,r,g,b,mid,lim);
        end
    end
    thick{k} = mat;
    thick_avg(k) = mean(mat(y1:y2, x1:x2),"all");
    thick_avg(k)
end
plot(thick_avg)

% 创建视频对象
v = VideoWriter('mesh_video.mp4', 'MPEG-4');
v.FrameRate = frameRate;   % 帧率（可调）
open(v);

figure;

for k = 1:N
    mesh(thick{k})
    zlim([0,800]);   % 固定坐标范围（重要）
    title(['Frame = ', num2str(k)]);
    
    drawnow;
    
    % ==== 捕获当前帧 ====
    frame = getframe(gcf);
    
    % ==== 写入视频 ====
    writeVideo(v, frame);
end

% 关闭视频
close(v);
% customCMap = uint8(rgbmatrix(:,min(mat,[],'all'):max(mat,[],'all'))');
% colormap(customCMap);
% colorbar
% mesh(mat);

function thick = rgb2Thick(rgbmatrix, y1, y2, y3, mid, lim)

    % 处理 mid
    if isnan(mid)
        mid = 1;   % 比原来更合理，避免 start=0
    end

    % 限定扫描范围
    start = max(1, mid - lim);
    endp  = min(size(rgbmatrix,2), mid + lim);

    % 取子区间
    sub = rgbmatrix(:, start:endp);   % 3×N

    % 目标向量
    target = [y1; y2; y3];

    % 向量化计算欧氏距离平方
    diff = sub - target;              % 自动广播 (3×N)
    dist = sum(diff.^2, 1);           % 1×N

    % 找最小值位置（局部索引）
    [~, local_idx] = min(dist);

    % 转换为全局索引
    thick = start + local_idx - 1;

end

% function [thick] = rgb2Thick(rgbmatrix,y1,y2,y3,mid,lim)
%     if isnan(mid)
%         mid = 0;
%     end
%     start = max(1, mid - lim);
%     endp = min(mid + lim, size(rgbmatrix,2));
%     delta = (rgbmatrix(1,start)-y1)^2+(rgbmatrix(2,start)-y2)^2+(rgbmatrix(3,start)-y3)^2;
%     res = 0;
%     for i = start+1:endp
%         delta0 = (rgbmatrix(1,i)-y1)^2+(rgbmatrix(2,i)-y2)^2+(rgbmatrix(3,i)-y3)^2;
%         if delta0 < delta
%             delta = delta0;
%             res = i;
%         end
%     end
%     thick = res;
%     if thick == 0
%         thick = mid;
%     end
% end

function [co] = rect_mean_exclude_nan(mat, x, y, i)
    if x == 1 && y == 1
        co = 0;
    else
        [m, n] = size(mat);  % 获取矩阵尺寸
    
        % 计算有效的矩形范围
        row_start = max(1, x - i);
        row_end   = min(m, x);
        col_start = max(1, y - i);
        col_end   = min(n, y);
    
        % 提取子矩阵
        submat = mat(row_start:row_end, col_start:col_end);
        co = int16(mean(submat(~isnan(submat))));
    end
end
