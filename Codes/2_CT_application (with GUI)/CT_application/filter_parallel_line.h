#pragma once
#ifndef FILTER_PARALLEL_H_INCLUDED
#define FILTER_PARALLEL_H_INCLUDED

using namespace std;
//using namespace cv;

//#define  pi  3.14159265358979323846  //!< const value,pi

enum filterFlags {
	/** choose the filter for the radon image */
	non_filter = 0,
	RL_filter = 1,//Ìí¼ÓÂË²¨|w|£¨¾ØÐÎ´°£©
	SL_filter = 2//Ìí¼ÓÂË²¨|w|*sinc(w/2B)
};


cv::Mat Filter_radon_parallel_line(cv::Mat inputImg, int filter_flag = 0);

#endif // !FILTER_PARALLEL_H_INCLUDED