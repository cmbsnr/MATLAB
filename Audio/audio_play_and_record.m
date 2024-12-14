Fs = 44100; %采样率
nbit = 24;%位深度
recObj = audiorecorder(Fs,nbit,1);%录音对象
recDuration = 5;%录音持续时间，单位s
f = 200;%固定频率，单位Hz
sf = 5000;%扫频
%发声
t = 0:1/Fs:recDuration;
audioData = sin(2*pi*f*t);
%sf_t = 0:1/Fs:recDuration;
sf_audioData = chirp(t,0,recDuration,sf);
sound(audioData,Fs,nbit)
%录音
disp("开始录音.")
recordblocking(recObj,recDuration);
disp("结束录音")
y = getaudiodata(recObj);
%傅里叶变换
%function fft
% figure
% t0 = (0:length(y)-1)/Fs;
% subplot(211)
% plot(t0,y)
% fftdata = fft(y);
% mag = abs(fftdata)*2/Fs;
% f = (0:length(y)-1)*Fs/length(y);
% subplot(212)
% plot(f,mag)

Y = fft(y);
L = length(y);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')