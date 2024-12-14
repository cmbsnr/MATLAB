clc
clear
close all
%%
Fs = 44100;      
recObj = audiorecorder(Fs,24,1);
% SampleRate = 44100;
% recObj.BitsPerSample = 24;
recDuration = 5;
disp("Begin speaking.")
recordblocking(recObj,recDuration);
disp("End of recording.")
y = getaudiodata(recObj);
%%
figure
t0 = (0:length(y)-1)/Fs;
subplot(211)
plot(t0,y)
data = y(10001:10000+Fs);
fftdata = fft(data);
mag = abs(fftdata)*2/Fs;
f = (0:length(data)-1)*Fs/length(data);
subplot(212)
plot(f,mag)

deviceWriter = audioDeviceWriter;
t = 0:1/Fs:10;
audioData = sin(2*pi*1000*t);
% sound(audioData,Fs)
%audiowrite('11.wav',[audioData;audioData]',Fs)

