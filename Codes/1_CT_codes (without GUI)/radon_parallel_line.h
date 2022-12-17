#pragma once
#ifndef RADON_PARALLEL_LINE_H_INCLUDED
#define RADON_PARALLEL_LINE_H_INCLUDED

using namespace std;
using namespace cv;

#define  pi  3.14159265358979323846  //!< const value,pi



static int Set_sensors_length_parallel_line(const Mat inputImg);
static void Set_Transmittors_and_Sensors_parallel_line(const int theta_degree, const int sensor_length,
	Point2d* Trm, Point2d* Rcv);

static double Cal_distance_point_to_line(Point2d Point, Point2d A1, Point2d A2);
static double Get_grayscale_iradon(double distance, const Mat radon_Img,
	const int sensor_length, const int theta_degree);



Mat Radon_parallel_line(Mat image, int num_degrees=180);
Mat iRadon_parallel_line(Mat radon_Img);//ÏñËØÇý¶¯


#endif // !RADON_PARALLEL_LINE_H_INCLUDED


