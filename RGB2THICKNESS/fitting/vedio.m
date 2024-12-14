vidObj = VideoReader("test.mp4");
while hasFrame(vidObj)
    vidFrame = readFrame(vidObj);

end
