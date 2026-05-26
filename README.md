# MATLAB

MATLAB 科研代码仓库，涵盖**薄膜干涉测厚**和**分子动力学模拟**两大方向，以及音频处理、图论、数据平滑等工具模块。

## 模块说明

### `RGB2THICKNESS/` — 薄膜干涉颜色→厚度转换

利用薄膜干涉颜色反推薄膜厚度。

- **`CIE XYZ/`** — 色度学基础：CIE 标准色匹配函数（1931/1964/Judd-Vos/Stiles-Burch）、多层膜干涉反射率模拟（Fresnel 方程）、Cauchy 色散方程
- **`fitting/`** — 核心拟合管线。通过预计算的 RGB-厚度参考表（`rgbmatrix`，3×N 矩阵），用最近邻查找将图像每个像素的 RGB 值映射为厚度。主要数据文件：`rm_s.mat`、`rm_f.mat`、`rm_f_2000.mat`
- **`time/`** — 时间分辨蒸发分析。逐帧读取视频，追踪薄膜厚度随时间变化

### `Dump/` — LAMMPS 轨迹后处理

读取 LAMMPS dump 文件（`.lammpstrj`），计算系统能量、应力、密度分布、对距离、自相关函数，构建 1-2/1-3/1-4 近邻列表用于力场参数分配。

### `ord2data/` — 分子文件转换

将 `.ord` 分子结构文件转换为 LAMMPS data 格式（用于 PDMS 模拟）。自动生成键、角、二面角拓扑信息。

### `Audio/` — 音频处理

录音、回放、FFT 频谱分析、啁啾信号生成。

### `movavg/` — 移动平均平滑

读取 CSV 格式的粘附力数据，应用中心移动平均窗口平滑，绘制力-深度曲线。

### `graph/` — 图论

枚举无向图中所有指定长度的路径（禁止立即折返），用于分子拓扑分析。

## 环境依赖

- MATLAB（基础模块 + Image Processing Toolbox）
- 部分脚本使用 `VideoReader` / `VideoWriter` 处理视频

## 使用方式

所有脚本在 MATLAB 中直接运行。`.mlx` 文件为 MATLAB Live Script（交互式笔记本），推荐在 MATLAB 编辑器中打开。多数脚本以 `clear; clc;` 开头。

数据文件（`.mat`、`.mp4`、`.png`、`.jpg`）已加入 `.gitignore`，中间产物需重新生成。
