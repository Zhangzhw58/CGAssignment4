# 光线追踪
## 中山大学计算机图形学期末大作业

### 运行环境

Windows + Visual Studio 2022
c++ 实现，依赖库GLM2

### 项目结构

- src: 项目主要代码，实现了光线追踪的主要功能
- image: 在项目中的场景中用到的图片文件

### 环境配置

请按`readme.pdf`配置好运行环境

### 运行参数

`src/main.cpp`中:

- aspect_ratio = 3.0 / 2.0; // 图像长宽比
- gWidth = 800; // 图像宽度
- samples_per_pixel = 5000;  // 每像素点采样数，值越大图像效果越好，运行时间越长
- scene_choice = 0; // 场景选择，目前支持(0~8)，0为默认最终场景