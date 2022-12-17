#pragma once
#ifndef BASIC_FUNCTIONS_H_INCLUDED
#define BASIC_FUNCTIONS_H_INCLUDED

using namespace std;
using namespace cv;

float sumMat(Mat& inputImg);
int max(int a, int b);
int min(int a, int b);
double max(double a, double b);
double min(double a, double b);
double max_3(double a, double b, double c);
double min_3(double a, double b, double c);
double distance_2d(Point2d A, Point2d B);

#endif // !BASIC_FUNCTIONS_H_INCLUDED

