#pragma once
#ifndef BASIC_RADON_FUNCTIONS_H_INCLUDED
#define BASIC_RADON_FUNCTIONS_H_INCLUDED

using namespace std;
using namespace cv;

#define  pi  3.14159265358979323846  //!< const value,pi

void Calculate_range_of_alpha(const int i, const Point2d* Trm, const Point2d* Rcv, const Mat image,
	bool& flag_intersection, double& min_alpha, double& max_alpha);
void Calculate_range_of_index_to_image(const int i, const Point2d* Trm,
	const Point2d* Rcv, const Mat image, const double min_alpha, const double max_alpha,
	int& intersection_n_min, int& intersection_n_max, int& intersection_m_min, int& intersection_m_max);
void Calculate_alpha_x(const int i, const Point2d* Trm, const Point2d* Rcv,
	const Mat image, const int intersection_n_min, const int sizeof_alpha_x, double* alpha_x);
void Calculate_alpha_y(const int i, const Point2d* Trm, const Point2d* Rcv,
	const Mat image, const int intersection_m_min, const int sizeof_alpha_y, double* alpha_y);
void Calculate_alpha(const int sizeof_alpha_x, const int sizeof_alpha,
	const double* alpha_x, const double* alpha_y, double* alpha);
double Calculate_projection_oneValue(const int i, const Point2d* Trm,
	const Point2d* Rcv, const Mat image, const int sizeof_alpha, const double* alpha);
void Project_copyto_2Darray(const int theta_degree, const double* projection_1d, Mat& radon_Img);

#endif // !BASIC_RADON_FUNCTIONS_H_INCLUDED
