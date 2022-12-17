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
    cv::Mat original_image; //�洢ԭʼͼ����Ϊ���ȶԵı�׼��
    cv::Mat image; //��ʼ��ͼ��
    int radon_flag = radonFlags::analytical; // ͶӰģʽ��������or���ֻ��֣���ѡ��Ĭ��shepplogan������
    //int radon_flag = radonFlags::integral;
    int radon_num_degree = 180;//���ֻ��֣�����ͼƬ����ͶӰ�Ƕ�
    //int num_circles_shepplogan = shepplogan_const::num_circles;
    int size_of_shepplogan = 256;//����ģ�͵�ͼ��ߴ磨Ĭ�ϲ���256�������룩
    bool use_filter = 1;//ʹ���˲�����Ĭ��ʹ�ã����˲���ͶӰ��������ʹ�ã���Ϊֱ�ӷ�ͶӰ
    bool show_reconstruction_process = 1;//��ʱ��ʾ�ؽ����̣�Ĭ�ϲ���ʾ����������ʾ�����
    int reconstruction_speed = 5;//�ؽ��ٶȣ�1-5��
    bool first_time_flag_text = 1;//�Ƿ��ǵ�һ�ν��롰�ı����õĽ�������shepploganͼ��ߴ磨����򣩡��ĺ��������������������ʾ����
    bool first_time_flag_slider = 1;

private:
    void displayMat_color(QLabel* label, cv::Mat image);
    void displayMat_gray(QLabel* label, cv::Mat image);

private slots:
    void display_img_click(void);//�������ʾͼƬ������ť��
    void size_of_shepplogan_changed(int); //�ı����õĽ�������shepploganͼ��ߴ硾�����
    void size_of_shepplogan_changed_slider(int value);//�ı����õĽ�������shepploganͼ��ߴ硾��������
    void num_degrees_changed(int value);//�ı����ֻ��ַ���ͶӰ�ܽǶȡ������
    void num_degrees_changed_slider(int value);//�ı�ı����ֻ��ַ���ͶӰ�ܽǶȡ���������
    void speed_changed_slider(int value);//�ı��ؽ��ٶȵġ���������
    void show_reconstruct_result_click(void);//�������ʼ�ؽ�����ť
    void use_filter_check(bool value);//�޸ġ�ʹ���˲������Ĺ�ѡ״̬������ѡ��ͬ�ڡ�ֱ�ӷ�ͶӰ������ѡ��ʹ��S-L�˲���
    void show_process_check(bool value);//�޸ġ���ʾ�ؽ����̡��Ĺ�ѡ״̬
    void comboBox_imgtype_changed(int value);//�޸�ͶӰ���ͣ�������or���ֻ��֣��ġ���ѡ��
};


#endif/*CT_APPLICATION_H*/