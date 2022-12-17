#pragma once
#ifndef RADON_PARALLEL_LINE_H_INCLUDED
#define RADON_PARALLEL_LINE_H_INCLUDED

using namespace std;
//using namespace cv;

//#define  pi  3.14159265358979323846  //!< const value,pi

enum radonFlags {
	/** choose the method for radon */
	analytical = 1, // for shepplogan
	integral = 2 // for existing image
};

static int Set_sensors_length_parallel_line(const cv::Mat inputImg);
static void Set_Transmittors_and_Sensors_parallel_line(const int theta_degree, const int sensor_length,
	cv::Point2d* Trm, cv::Point2d* Rcv);

static double Cal_distance_point_to_line(cv::Point2d Point, cv::Point2d A1, cv::Point2d A2);
static double Get_grayscale_iradon(double distance, const cv::Mat radon_Img,
	const int sensor_length, const int theta_degree);



cv::Mat Radon_parallel_line(cv::Mat image, int num_degrees=180);
cv::Mat iRadon_parallel_line(cv::Mat radon_Img);

void Radon_parallel_line_show(cv::Mat image, int sensor_length, int current_degree, cv::Mat& radon_Img);
void iRadon_parallel_line_show(cv::Mat radon_Img, int current_degree, cv::Mat& output_Img);


#endif // !RADON_PARALLEL_LINE_H_INCLUDED


