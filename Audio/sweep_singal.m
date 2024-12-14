% 设置参数
fs = 44100; % 采样频率，常用的值为44100Hz
t = 10; % 持续时间，单位为秒
f_start = 1000; % 起始频率，单位为Hz
f_end = 10000; % 结束频率，单位为Hz

% 生成时间向量
time = linspace(0, t, fs * t);

% 生成扫频信号
sweep_signal = chirp(time, f_start, t, f_end);

% 保存为wav格式
filename = 'sweep_signal.wav';
audiowrite(filename, sweep_signal, fs);

% 播放生成的音频
sound(sweep_signal, fs);

disp('扫频音频已生成并保存为 sweep_signal.wav');
