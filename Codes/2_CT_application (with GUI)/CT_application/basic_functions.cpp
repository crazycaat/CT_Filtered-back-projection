#include <iostream>
#include<opencv2/opencv.hpp>

#include "basic_functions.h"

//#include <opencv2/highgui/highgui_c.h> 

using namespace std;
using namespace cv;


float sumMat(Mat& inputImg)
{
	float sum = 0.0;
	int rowNumber = inputImg.rows;
	int colNumber = inputImg.cols * inputImg.channels();
	for (int i = 0; i < rowNumber; i++)
	{
		uchar* data = inputImg.ptr<uchar>(i);
		for (int j = 0; j < colNumber; j++)
		{
			sum = data[j] + sum;
		}
	}
	return	sum;
}


Mat Resize_img_to_Square(Mat img)
{
	int new_rows = max(img.rows, img.cols);
	Mat new_img = Mat::zeros(new_rows, new_rows, COLOR_BGR2BGRA);//uchar类型
	if (img.rows <=  img.cols)//列大于行
	{
		int padding_length = floor((img.cols- img.rows) / 2);
		for (int i = padding_length; i < padding_length + img.rows; ++i)
		{
			for (int j = 0; j < img.cols; ++j)
			{
				uchar tmp = img.at<uchar>(i - padding_length, j);
				new_img.at<uchar>(i, j) = tmp;
			}
		}
	}
	else//行大于列
	{
		int padding_length = floor((img.rows - img.cols) / 2);
		for (int i = 0; i < img.rows; ++i)
		{
			for (int j = padding_length; j < padding_length + img.cols; ++j)
			{
				uchar tmp = img.at<uchar>(i, j - padding_length);
				new_img.at<uchar>(i, j) = tmp;
			}
		}
	}

	return new_img;

}


int max(int a, int b)
{
	if (a > b) return a;
	else return b;
}

int min(int a, int b)
{
	if (a < b) return a;
	else return b;
}

double max(double a, double b)
{
	if (a > b) return a;
	else return b;
}

double min(double a, double b)
{
	if (a < b) return a;
	else return b;
}

double max_3(double a, double b, double c)
{
	double max = a;
	if (b > max) max = b;
	if (c > max) max = c;
	return max;
}

double min_3(double a, double b, double c)
{
	double min = a;
	if (b < min) min = b;
	if (c < min) min = c;
	return min;
}

double distance_2d(Point2d A, Point2d B)
{
	double distance;
	distance = sqrt(pow(A.x-B.x, 2) + pow(A.y - B.y, 2));
	return distance;
}

