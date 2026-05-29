import numpy as np
import cv2
import scipy.io as sio
import matplotlib.pyplot as plt

# === 1. 读取 .mat 数据 ===
data = sio.loadmat('rm_f_2000.mat')
rgbmatrix = data['rgbmatrix']   # shape: (3, N)

# === 2. 参数初始化 ===
index = 0
lim = 50

videoFile = 'Sout.mp4'

# === 3. 视频读取 ===
cap = cv2.VideoCapture(videoFile)

frames = []
timeList = []

fps = cap.get(cv2.CAP_PROP_FPS)

frame_id = 0
while True:
    ret, frame = cap.read()
    if not ret:
        break
    
    # OpenCV是BGR → 转RGB
    frame = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
    
    frames.append(frame)
    timeList.append(frame_id / fps)
    
    frame_id += 1

cap.release()

frames = np.array(frames)
timeList = np.array(timeList)

N = len(frames)

# === 4. 区域参数 ===
H, W, _ = frames[0].shape

y1, y2 = 160, 1840
x1, x2 = 1800, 2040

# === 5. 结果初始化 ===
thick = np.ones(N)
rgb = np.ones((N, 3))

# === 6. 逐帧处理（倒序） ===
for i in range(N-1, -1, -1):
    frame = frames[i]
    
    # 中心区域
    frame_center = frame[y1:y2, x1:x2, :]
    
    # 平均RGB
    rgb_mean = frame_center.mean(axis=(0, 1))  # shape (3,)
    
    # === 搜索窗口 ===
    left = max(0, index)   # Python从0开始
    right = min(rgbmatrix.shape[1], index + 2 * lim)
    
    sub_rgb = rgbmatrix[:, left:right]
    
    # === 欧氏距离 ===
    diff = sub_rgb - rgb_mean.reshape(3, 1)
    dist = np.sum(diff**2, axis=0)
    
    # === 找最优 ===
    local_idx = np.argmin(dist)
    idx = left + local_idx
    
    rgb[i, :] = rgb_mean
    thick[i] = idx   # 单位 nm
    
    index = idx  # 更新搜索中心

# === 7. 绘图 ===
plt.figure()
plt.plot(timeList, thick, 'o')
plt.xlabel('Time (s)')
plt.ylabel('Film thickness (nm)')
plt.title('Evaporate Viscosity')
plt.show()