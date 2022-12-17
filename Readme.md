You may open the video: https://www.bilibili.com/video/BV1od4y1K71J/?spm_id_from=333.999.0.0&vd_source=6968856c8aca1ad25b542920e96c4e43 for further illustration.

# Introduction

该软件面向 CT滤波反投影算法的初学者，希望通过参数的交互，让使用者加深对该算法的理解。 

编程及项目管理：Visual Studio 2019 

图像处理：openCV 4.4.0（C++） 

交互式界面设计：Qt 5.12.2（C++） 

包含数字积分法和解析法（shepplogan模型）投影、S-L滤波、反投影过程。 

用户交互性强，支持多个参数的修改。 

本地图片支持彩色输入（输出为灰度图），支持长方形图片输入，图像类型支持*.png,*.jpg,*.bmp。 

提供重建评价指标，便于重建效果评估。 

注：该软件仅展示了平行束滤波反投影算法（不包含扇形束、锥形束等） 

# FILES

```c++
├── CT_application_boxed.exe               	// a software that can run independently
├── Test_pictures/         					
│   └── shepplogan.png       				// shepplogan-Created by program
├── Files/									//Illustration Files
│   ├── Introduction_of_What_These_Codes_Mean.pdf 
│   └── How_to_Use_the_Software.pdf        
├── Codes/              					// all the codes
│   ├── 1_CT_codes (without GUI)/	           	
│   └──	2_CT_application (with GUI)/			
```

# Reference

1. Projection of SheppLogan Model parameters and analytic method

Avinash C. Kak, and Malcolm Slaney, “Principles of Computerized Tomographic Imaging”, IEEE Press, 1999. https://engineering.purdue.edu/~malcolm/pct/，Chapter 3. Algorithms for Reconstruction with Nondiffracting Sources. 

2. integral projection algorithm

Siddon, Robert L. "Fast calculation of the exact radiological path for a three‐dimensional CT array." Medical physics 12.2 (1985): 252-255.