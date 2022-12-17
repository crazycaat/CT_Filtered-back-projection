#pragma once
#ifndef BASIC_DISPLAY_FUNCTIONS_H_INCLUDED
#define BASIC_DISPLAY_FUNCTIONS_H_INCLUDED

using namespace std;
//using namespace cv;

cv::Mat Convert_to_show_image(cv::Mat Img);
cv::Mat Convert_to_show_normalize(cv::Mat Img);
cv::Mat Convert_to_show_equalizeHist(cv::Mat inputImg);
cv::Mat Convert_to_show_LogTrans(cv::Mat inputImg, double gamma = 1.5);

void showImg(string windowname, cv::Mat img);
void write_and_show_Img(cv::Mat image, string image_name, string type_name);
int WriteData(string fileName, cv::Mat& matData);
int LoadData(string fileName, cv::Mat& matData, int matRows, int matCols, int matChns);

double Calculate_error_sqrt(cv::Mat image1, cv::Mat image2);
double Calculate_error_abs(cv::Mat image1, cv::Mat image2);


#endif // !BASIC_DISPLAY_FUNCTIONS_H_INCLUDED

