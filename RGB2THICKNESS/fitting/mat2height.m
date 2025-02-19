videoFile = 'output_reversed.mp4';
v = VideoReader(videoFile);

outputFile = 'height2.mp4';
vOut = VideoWriter(outputFile, 'MPEG-4');
vOut.FrameRate = v.FrameRate; % 保持相同帧率
open(vOut);

for i = 1:length(mats)
    mesh(mats{i}')
    zlim([0 500])
    frame = getframe(gcf);
    writeVideo(vOut, frame);
end


close(vOut);

disp('视频处理完成，已保存为 height.mp4');