#pragma once
#pragma once
#ifndef SHEPPLOGAN_H_INCLUDED
#define SHEPPLOGAN_H_INCLUDED

using namespace std;
//using namespace cv;

//#define  pi  3.14159265358979323846  //!< const value,pi

#include "CT_application.h"
enum shepplogan_const {
	num_circles = 10 // shepplogan中的椭圆个数
};

cv::Mat Create_shepplogan(int rows = 256);
cv::Mat Calc_shepplogan_radon(int rows = 256);

void Calc_shepplogan_radon_show(int rows, int index, cv::Mat& radon_Img);

#endif // !SHEPPLOGAN_H_INCLUDED