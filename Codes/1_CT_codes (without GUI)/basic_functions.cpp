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



