#pragma once
#ifndef CT_APPLICATION_H
#define CT_APPLICATION_H

#include <QtWidgets/QMainWindow>
#include "ui_CT_application.h"

#include<iostream>
#include<opencv2/opencv.hpp>

#include "basic_display_functions.h"
#include "basic_functions.h"
#include "shepplogan.h"
#include "radon_parallel_line.h"
#include "filter_parallel_line.h"

class CT_application : public QMainWindow
{
    Q_OBJECT

public:
    CT_application(QWidget *parent = Q_NULLPTR);
    

private:
    Ui::CT_applicationClass ui; 
    cv::Mat original_image; //存储原始图像（作为最后比对的标准）
    cv::Mat image; //初始化图像
    int radon_flag = radonFlags::analytical; // 投影模式（解析法or数字积分，可选择）默认shepplogan解析法
    //int radon_flag = radonFlags::integral;
    int radon_num_degree = 180;//数字积分（本地图片）的投影角度
    //int num_circles_shepplogan = shepplogan_const::num_circles;
    int size_of_shepplogan = 256;//解析模型的图像尺寸（默认参数256，可输入）
    bool use_filter = 1;//使用滤波器（默认使用，即滤波反投影），若不使用，则为直接反投影
    bool show_reconstruction_process = 1;//延时显示重建过程（默认不显示，即快速显示结果）
    int reconstruction_speed = 5;//重建速度（1-5）
    bool first_time_flag_text = 1;//是否是第一次进入“改变内置的解析法的shepplogan图像尺寸（输入框）”的函数，避免软件开启的显示问题
    bool first_time_flag_slider = 1;

private:
    void displayMat_color(QLabel* label, cv::Mat image);
    void displayMat_gray(QLabel* label, cv::Mat image);

private slots:
    void display_img_click(void);//点击【显示图片】【按钮】
    void size_of_shepplogan_changed(int); //改变内置的解析法的shepplogan图像尺寸【输入框】
    void size_of_shepplogan_changed_slider(int value);//改变内置的解析法的shepplogan图像尺寸【滑动条】
    void num_degrees_changed(int value);//改变数字积分法的投影总角度【输入框】
    void num_degrees_changed_slider(int value);//改变改变数字积分法的投影总角度【滑动条】
    void speed_changed_slider(int value);//改变重建速度的【滑动条】
    void show_reconstruct_result_click(void);//点击【开始重建】按钮
    void use_filter_check(bool value);//修改【使用滤波器】的勾选状态（不勾选等同于“直接反投影”，勾选后使用S-L滤波）
    void show_process_check(bool value);//修改【显示重建过程】的勾选状态
    void comboBox_imgtype_changed(int value);//修改投影类型（解析法or数字积分）的【复选框】
};


#endif/*CT_APPLICATION_H*/