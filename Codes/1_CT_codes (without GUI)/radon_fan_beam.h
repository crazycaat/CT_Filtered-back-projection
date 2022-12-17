#pragma once
#pragma once
#ifndef RADON_FAN_BEAM_H_INCLUDED
#define RADON_FAN_BEAM_H_INCLUDED

using namespace std;
using namespace cv;

#define  pi  3.14159265358979323846  //!< const value,pi

Mat Radon_fan_beam(Mat image, int num_rotate_degree = 360, int num_detector_degree = 60);
//Mat iRadon_fan_beam(Mat radon_Img);

static void Set_Transmittors_and_Sensors_fan_beam(const int current_rotate_degree, const int num_detector_degree,
	const double R_circle, Point2d* Trm, Point2d* Rcv);



#endif // !RADON_FAN_BEAM_H_INCLUDED



