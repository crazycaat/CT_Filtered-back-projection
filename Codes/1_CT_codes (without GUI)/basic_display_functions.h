#pragma once
#ifndef BASIC_DISPLAY_FUNCTIONS_H_INCLUDED
#define BASIC_DISPLAY_FUNCTIONS_H_INCLUDED

using namespace std;
using namespace cv;

Mat Convert_to_show_image(Mat Img);
Mat Convert_to_show_normalize(Mat Img);
Mat Convert_to_show_equalizeHist(Mat inputImg);
Mat Convert_to_show_LogTrans(Mat inputImg, double gamma = 1.5);

void showImg(string windowname, Mat img);
void write_and_show_Img(Mat image, string image_name, string type_name);
int WriteData(string fileName, cv::Mat& matData);
int LoadData(string fileName, cv::Mat& matData, int matRows, int matCols, int matChns);

double Calculate_error_sqrt(Mat image1, Mat image2);
double Calculate_error_abs(Mat image1, Mat image2);

#endif // !BASIC_DISPLAY_FUNCTIONS_H_INCLUDED

