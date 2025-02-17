% 读取原始视频
videoFile = 'test.mp4';
v = VideoReader(videoFile);

% 定义起始时间和结束时间（单位：秒）
startTime = 29;
endTime = 39;

% 读取指定时间段的视频帧
frames = {}; % 用于存储帧
while hasFrame(v)
    frame = readFrame(v);
    if v.CurrentTime >= startTime && v.CurrentTime <= endTime
        frames{end+1} = frame; % 存储帧
    elseif v.CurrentTime > endTime
        break; % 超出范围，停止读取
    end
end

% 倒序帧
frames = flip(frames);

% 创建新视频文件
outputFile = 'output_reversed.mp4';
vOut = VideoWriter(outputFile, 'MPEG-4');
vOut.FrameRate = v.FrameRate; % 保持相同帧率
open(vOut);

% 写入倒序的视频帧
for i = 1:length(frames)
    writeVideo(vOut, frames{i});
end

% 关闭视频文件
close(vOut);

disp('视频处理完成，已保存为 output_reversed.mp4');