#include <iostream>
#include<opencv2/opencv.hpp>

#include "radon_parallel_line.h"
#include "basic_functions.h"

//#include <opencv2/highgui/highgui_c.h> 

using namespace std;
using namespace cv;


/** @brief calculate min_alpha and max_alpha
alpha:The ratio of the distance from the transmitter to the intersection
		to the distance from the transmitter to the receiver
Reference:Siddon, Robert L.
		"Fast calculation of the exact radiological path for a three‐dimensional CT array."
		Medical physics 12.2 (1985): 252-255.

@param	i: The serial number of the transmitter(or the receiver);
@param	Trm:record all the positions of transmittors;
@param	Rcv:record all the positions of receivers;
@param	image:the original input image,to calculate the range of x and y;
@return	flag_intersection: turns 0 if the ray has no intersection with the image;
		CAUTION!It will be modified by this function.
@return	min_alpha:the minimum value of the ratio alpha;
		CAUTION!It will be modified by this function.
@return	min_alpha:the maximum value of the ratio alpha;
		CAUTION!It will be modified by this function.

 */
void Calculate_range_of_alpha(const int i, const Point2d* Trm, const Point2d* Rcv, const Mat image,
	bool& flag_intersection, double& min_alpha, double& max_alpha)
{
	/*   确定图像所在的空间范围   */
	double img_x_min = -0.5 * image.cols;
	double img_x_max = img_x_min + image.cols;
	double img_y_min = -0.5 * image.rows;
	double img_y_max = img_y_min + image.rows;
	double min_alpha_x = 1.0, min_alpha_y = 1, max_alpha_x = 0.0, max_alpha_y = 0.0;

	if ((Rcv[i].x - Trm[i].x) != 0)
	{
		min_alpha_x = min((img_x_min - Trm[i].x) / (Rcv[i].x - Trm[i].x),
			(img_x_max - Trm[i].x) / (Rcv[i].x - Trm[i].x));
		max_alpha_x = max((img_x_min - Trm[i].x) / (Rcv[i].x - Trm[i].x),
			(img_x_max - Trm[i].x) / (Rcv[i].x - Trm[i].x));
	}
	if ((Rcv[i].y - Trm[i].y) != 0)
	{
		min_alpha_y = min((img_y_min - Trm[i].y) / (Rcv[i].y - Trm[i].y),
			(img_y_max - Trm[i].y) / (Rcv[i].y - Trm[i].y));
		max_alpha_y = max((img_y_min - Trm[i].y) / (Rcv[i].y - Trm[i].y),
			(img_y_max - Trm[i].y) / (Rcv[i].y - Trm[i].y));
	}
	//cout << "min_alpha_x = " << min_alpha_x << ", max_alpha_x = " << max_alpha_x << endl;
	//cout << "min_alpha_y = " << min_alpha_y << ", max_alpha_y = " << max_alpha_y << endl;
	min_alpha = max_3(min_alpha_x, min_alpha_y, 0.0);
	//cout << "min_alpha = " << min_alpha<<endl;
	max_alpha = min_3(max_alpha_x, max_alpha_y, 1.0);
	//cout << "min_alpha = " << min_alpha << ", max_alpha = " << max_alpha << endl;
	if (min_alpha >= max_alpha)//没有交点，或者仅有一个点的交点（而不是线段），投影为0
	{
		flag_intersection = 0;
		//cout << "NO INTERSECTION!theta = " << theta_degree << ", number(i) of transmittor = " << i << endl;
	}
}


/** @brief Calculate the image indices m and n of the part where the ray intersects the image.
@param	i: The serial number of the transmitter(or the receiver);
@param	Trm:record all the positions of transmittors;
@param	Rcv:record all the positions of receivers;
@param	image:the original input image,to calculate the range of x and y;
@param 	min_alpha:the minimum value of the ratio alpha;
@param	min_alpha:the maximum value of the ratio alpha;
@return	intersection_n_min:Minimum value of column indices where the ray intersects the image
		CAUTION!It will be modified by this function.
@return	intersection_n_max:Maximum value of column indices where the ray intersects the image
		CAUTION!It will be modified by this function.
@return	intersection_m_min:Minimum value of row indices where the ray intersects the image
		CAUTION!It will be modified by this function.
@return	intersection_m_max:Maximum value of row indices where the ray intersects the image
		CAUTION!It will be modified by this function.
 */
void Calculate_range_of_index_to_image(const int i, const Point2d* Trm,
	const Point2d* Rcv, const Mat image, const double min_alpha, const double max_alpha,
	int& intersection_n_min, int& intersection_n_max, int& intersection_m_min, int& intersection_m_max)
{
	/*   确定图像所在的空间范围   */
	double img_x_min = -0.5 * image.cols;
	double img_x_max = img_x_min + image.cols;
	double img_y_min = -0.5 * image.rows;
	double img_y_max = img_y_min + image.rows;

	if (Trm[i].x < Rcv[i].x)
	{
		intersection_n_min = (image.cols) - (img_x_max - min_alpha * (Rcv[i].x - Trm[i].x) - Trm[i].x);
		intersection_n_max = -1 + (Trm[i].x + max_alpha * (Rcv[i].x - Trm[i].x) - img_x_min);
	}
	else
	{
		intersection_n_min = (image.cols) - (img_x_max - max_alpha * (Rcv[i].x - Trm[i].x) - Trm[i].x);
		intersection_n_max = -1 + (Trm[i].x + min_alpha * (Rcv[i].x - Trm[i].x) - img_x_min);
	}
	//m范围（m为图像的行，在坐标系中对应y）
	//由于随x增大，n增大；而随y增大，m减小，所以公式需要做调整(与文献中的方向有差别）。
	if (Trm[i].y < Rcv[i].y)
	{
		intersection_m_min = (image.rows) - (Trm[i].y + max_alpha * (Rcv[i].y - Trm[i].y) - img_y_min);
		intersection_m_max = img_y_max - min_alpha * (Rcv[i].y - Trm[i].y) - Trm[i].y - 1;
	}
	else
	{
		intersection_m_min = (image.rows) - (Trm[i].y + min_alpha * (Rcv[i].y - Trm[i].y) - img_y_min);
		intersection_m_max = img_y_max - max_alpha * (Rcv[i].y - Trm[i].y) - Trm[i].y - 1;
	}
}

/** @brief Calculate all the the ratios alpha_x that intersect the line x=a.
@param	i: The serial number of the transmitter(or the receiver);
@param	Trm:record all the positions of transmittors;
@param	Rcv:record all the positions of receivers;
@param	image:the original input image,to calculate the range of x and y;
@param	intersection_n_min:Minimum value of column indices where the ray intersects the image
@param	sizeof_alpha_x:The size of array alpha_x ( = the number of insections with the line x=a)
@return alpha_x: all the the ratios alpha_x that intersect the line x=a.
		CAUTION!It will be modified by this function.
 */
void Calculate_alpha_x(const int i, const Point2d* Trm, const Point2d* Rcv,
	const Mat image, const int intersection_n_min, const int sizeof_alpha_x, double* alpha_x)
{
	for (int k = 0; k < sizeof_alpha_x; ++k)
	{
		int x_position = intersection_n_min + k - image.cols / 2;
		alpha_x[k] = (x_position - Trm[i].x) / (Rcv[i].x - Trm[i].x);
	}
}

/** @brief Calculate all the the ratios alpha_y that intersect the line y=a.
@param	i: The serial number of the transmitter(or the receiver);
@param	Trm:record all the positions of transmittors;
@param	Rcv:record all the positions of receivers;
@param	image:the original input image,to calculate the range of x and y;
@param	intersection_m_min:Minimum value of row indices where the ray intersects the image
@param	sizeof_alpha_y:The size of array alpha_y ( = the number of insections with the line y=a)
@return alpha_y: all the the ratios alpha_y that intersect the line y=a.
		CAUTION!It will be modified by this function.
 */
void Calculate_alpha_y(const int i, const Point2d* Trm, const Point2d* Rcv,
	const Mat image, const int intersection_m_min, const int sizeof_alpha_y, double* alpha_y)
{
	for (int k = 0; k < sizeof_alpha_y; ++k)
	{
		int y_position = image.rows / 2 - (intersection_m_min + k);
		alpha_y[k] = (y_position - Trm[i].y) / (Rcv[i].y - Trm[i].y);
	}

}

/** @brief Calculate all the the ratios alpha that intersect the line x=a or y=a.
@param	sizeof_alpha_x:The size of array alpha_x ( = the number of insections with the line x=a)
@param	sizeof_alpha:The size of array alpha
@param	alpha_x: all the the ratios alpha_x that intersect the line x=a.
@param	alpha_y: all the the ratios alpha_y that intersect the line y=a.
@return	alpha: all the the ratios alpha that intersect the line x=a or y=a.
		CAUTION!It will be modified by this function.
 */
void Calculate_alpha(const int sizeof_alpha_x, const int sizeof_alpha,
	const double* alpha_x, const double* alpha_y, double* alpha)
{
	for (int k = 0; k < sizeof_alpha_x; ++k)
	{
		alpha[k] = alpha_x[k];
	}
	for (int k = sizeof_alpha_x; k < sizeof_alpha; ++k)
	{
		alpha[k] = alpha_y[k - sizeof_alpha_x];
	}

	sort(alpha, alpha + sizeof_alpha);//升序排列
}

/** @brief Calculate the projection at a certain degree received by certain receiver
@param	i: The serial number of the transmitter(or the receiver);
@param	Trm:record all the positions of transmittors;
@param	Rcv:record all the positions of receivers;
@param	image:the original input image,to get the grayscale value;
@param	sizeof_alpha:The size of array alpha
@param	alpha: all the the ratios alpha that intersect the line x=a or y=a.
@return	projection:the projection at a certain degree received by certain receiver(i)
 */
double Calculate_projection_oneValue(const int i, const Point2d* Trm,
	const Point2d* Rcv, const Mat image, const int sizeof_alpha, const double* alpha)
{
	double total_TR_length = distance_2d(Trm[i], Rcv[i]);
	double projection = 0.0;

	for (int k = 1; k < sizeof_alpha; ++k)
	{
		double length_of_little_route = (alpha[k] - alpha[k - 1]) * total_TR_length;
		double mid_alpha = 0.5 * (alpha[k] + alpha[k - 1]);
		Point2d M_pixel;//通过的像素的中点坐标
		M_pixel.x = Trm[i].x + mid_alpha * (Rcv[i].x - Trm[i].x);
		M_pixel.y = Trm[i].y + mid_alpha * (Rcv[i].y - Trm[i].y);

		int m_pixel, n_pixel;//对应的像素的索引
		//Caution!!n为列，对应x；m为行，对应y！！！！（就是这里弄反了，重建图像旋转了90度还有镜像）
		n_pixel = floor(M_pixel.x + 0.5 * image.cols);
		n_pixel = max(0, n_pixel);
		n_pixel = min(n_pixel, image.rows - 1);

		m_pixel = floor(0.5 * image.rows - M_pixel.y);
		m_pixel = max(0, m_pixel);
		m_pixel = min(m_pixel, image.cols - 1);
		double grayscale = image.at<uchar>(m_pixel, n_pixel);
		projection += length_of_little_route * grayscale;
	}
	return projection;
}


/** @brief copy the 1d projection to the 2d Mat(radon_Img)
@param	theta_degree:the rotation degree of transmittor and receiver
			0-180,represent the angle between the ray(Trm to Rcv) and the positive direction of the x-axis
			also corresponds the colomns of the Mat radon_Img
@param	projection_1d:the projection corresponds to a certain angle and a series of receivers
@return	radon_Img:Output 64-bit 1-channel image.(refresh a column)
		Size:rows = number of receivers = sqrt(2) * (maximum side of input image)
			 cols = number of degrees( default = 180)
		CAUTION!It will be modified by this function.

 */
void Project_copyto_2Darray(const int theta_degree, const double* projection_1d, Mat& radon_Img)
{
	for (int m = 0; m < radon_Img.rows; ++m)
	{
		radon_Img.at<double>(m, theta_degree) = projection_1d[m];
	}
}

