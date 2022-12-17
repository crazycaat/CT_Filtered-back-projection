#include <iostream>
#include<opencv2/opencv.hpp>

#include "radon_parallel_line.h"
#include "basic_functions.h"
#include "basic_radon_functions.h"

//#include <opencv2/highgui/highgui_c.h> 

using namespace std;
using namespace cv;

/** @brief set the sensor length

@param   inputImg: Input the original image (8-bit 1-channel image).
@return  sensor_length
		 sensor_length = the length of transmittors = the length of the receivers
			= the array size of Trm = the array size of Rcv = sqrt(2) * (maximum side of the image)

 */
static int Set_sensors_length_parallel_line(const Mat inputImg)
{
	//投影的尺寸计算
	//以图像的最长边为标准建立探测器长度，且为了在倾斜时能够覆盖所有图像像素，探测器的最小长度为图像边长的根号2倍
	int img_longest_edge;
	img_longest_edge = max(inputImg.rows, inputImg.cols);
	//cout << "img_longest_edge = " << img_longest_edge << endl;
	int sensor_length = ceil(img_longest_edge * sqrt(2)) + 2;
	//cout << "The radius of the midpoint of sensor = the minimum length of the sensor = " << sensor_length << endl;
	return sensor_length;
}


/** @brief set the positions of transmittors and sensors, and store the positions in the arrays called Trm and Rcv.

@param theta_degree:the rotation degree of transmittor and receiver
			0-180,represent the angle between the ray(Trm to Rcv) and the positive direction of the x-axis
@param  sensor_length 
		 sensor_length = the length of transmittors = the length of the receivers 
		    = the array size of Trm = the array size of Rcv = sqrt(2) * (maximum side of the image)
@param	Trm: record all the positions of transmittors when the angle equals theta
		CAUTION!The content of Trm will be modified by this function
@param	Rcv: record all the positions of receivers when the angle equals theta
		CAUTION!The content of Rcv will be modified by this function

@note
   -  the function will directly modify the parameters - Trm and Rcv
 */
static void Set_Transmittors_and_Sensors_parallel_line(const int theta_degree, const int sensor_length,
	Point2d* Trm,Point2d* Rcv)
{
	float theta = float(theta_degree) * pi / 180;

	//发射器中点
	Point2d mid_Trm;
	mid_Trm.x = 0.5 * (double(sensor_length) - 2) * cos(theta + pi);//离图像稍微近一点，保证45度时全覆盖（取整带来的问题），前面对应已经+2
	mid_Trm.y = 0.5 * (double(sensor_length) - 2) * sin(theta + pi);
	cout << "theta = " << theta_degree << " processing......" << endl;
	//cout << "mid_Trm = " << mid_Trm << endl;
	//发射器位置确定，Trm[0]~Trm[sensor_length-1]共sensor_length个
	//Point2d* Trm;
	
	Trm[0].x = mid_Trm.x + 0.5 * sensor_length * cos(0.5 * pi + theta);
	Trm[0].y = mid_Trm.y + 0.5 * sensor_length * sin(0.5 * pi + theta);
	//cout << "Trm[0]" << Trm[0] << endl;
	for (int i = 1; i < sensor_length; ++i)
	{
		Trm[i].x = Trm[i - 1].x - 1 * cos(0.5 * pi + theta);
		Trm[i].y = Trm[i - 1].y - 1 * sin(0.5 * pi + theta);
	}
	//cout << "Trm[" << sensor_length - 1 << "] = " << Trm[sensor_length - 1] << endl;


	//接收器中点
	Point2d mid_Rcv;
	mid_Rcv.x = 0.5 * double(sensor_length) * cos(theta);
	mid_Rcv.y = 0.5 * double(sensor_length) * sin(theta);
	//cout << "mid_Rcv = " << mid_Rcv << endl;
	//接收器位置确定，Rcv[0]~Rcv[sensor_length-1]共sensor_length个
	//Point2d* Rcv;
	
	Rcv[0].x = mid_Rcv.x + 0.5 * sensor_length * cos(0.5 * pi + theta);
	Rcv[0].y = mid_Rcv.y + 0.5 * sensor_length * sin(0.5 * pi + theta);
	//cout << "Rcv[0]" << Rcv[0] << endl;
	for (int i = 1; i < sensor_length; ++i)
	{
		Rcv[i].x = Rcv[i - 1].x - 1 * cos(0.5 * pi + theta);
		Rcv[i].y = Rcv[i - 1].y - 1 * sin(0.5 * pi + theta);
	}
	//cout << "Rcv[" << sensor_length - 1 << "] = " << Rcv[sensor_length - 1] << endl;
}



/** @brief Parallel line radon

@param image Input 8-bit 1-channel image.
@return  radon_Img Output 64-bit 1-channel image.  
		Size:rows = number of receivers = sqrt(2) * (maximum side of input image)
			 cols = number of degrees( default = 180)

The function references the article:
Siddon, Robert L. "Fast calculation of the exact radiological path for a three‐dimensional CT array."
Medical physics 12.2 (1985): 252-255.

@note
   -  the output image(radon_Img) has a narrow effective grayscale with a large grayscale range.
      Therefore, simple process of radon_Img before showing and storing the result is needed.@ref function-Convert_to_show_radon_parallel_line
 */
Mat Radon_parallel_line(Mat image, int num_degrees)
{
	//第一步：投影
	//投影的尺寸计算
	//以图像的最长边为标准建立探测器长度，且为了在倾斜时能够覆盖所有图像像素，探测器的最小长度为图像边长的根号2倍
	int sensor_length = Set_sensors_length_parallel_line(image);
	Mat radon_Img = Mat::zeros(sensor_length, num_degrees, COLOR_BGR2GRAY);//存储投影后的图像矩阵，一列为一个角度
		
	//degree大循环，默认参数为循环180度，步长为1度
	for (int theta_degree = 0; theta_degree < num_degrees; ++theta_degree)
	//for (int theta_degree = 0; theta_degree < 1; theta_degree+=90)
	{
		/*
		构造射线和探测器
		由于输入图像已经转换为单通道，所以image.at<uchar>(i,j)中i为图像的行，j为图像的列
		设Transmittor(Trm)为发射源，receiver(Rcv)为接收器，Trm，Rcv均为大小为sensor_length的数组。
		0度表示从左向右的射线，90度表示从下向上的射线
		以图像中心为原点建立直角坐标系，则发射源和探测器的中点永远在圆x^2+y^2 = (1/4)sensor_length^2上
		theta为射线方向与x轴正方向的夹角
		（根据45度时探测器与图像中心距离至少为1/2图像斜角长度，探测器长度=图像斜角长度计算而得，中点的轨迹的半径=1/2探测器长度）
		*/
		Point2d* Trm;
		Point2d* Rcv;
		Trm = new Point2d[sensor_length];
		Rcv = new Point2d[sensor_length];
		Set_Transmittors_and_Sensors_parallel_line(theta_degree, sensor_length, Trm, Rcv);//构造射线和探测器的位置

		/*    正式开始投影！！    */
		double* projection_1d;
		projection_1d = new double[sensor_length]; //记录特定角度下的各个探测器接收到的一系列投影

		//i小循环，特定某个角度theta_degree下的的各个发射器（or探测器）
		for (int i = 0; i < sensor_length; ++i)
		//for (int i = 55; i < 56; ++i)
		{
			//根据Trm[i]和Rcv[i]的坐标，求出对应Rcv[i]接收到的投影值projection_1d[i]
			//初始化第i个投影值，以及第i个发射器和接收器对应的alpha比例的最大和最小值
			//alpha的含义参考文献：Siddon, Robert L. "Fast calculation of the exact radiological path for a three‐dimensional CT array." Medical physics 12.2 (1985): 252-255.
			projection_1d[i] = 0.0;
			
			bool flag_intersection = 1;
			
			//计算min_alpha 和max_alpha
			double min_alpha = 0.0, max_alpha = 1.0;
			Calculate_range_of_alpha(i, Trm, Rcv, image, flag_intersection, min_alpha, max_alpha);
			if (flag_intersection==0){continue;}//没有交点，或者仅有一个点的交点（而不是线段），投影为0
			
			//代码运行到这里，表示已经获得了min_alpha，max_alpha，并且确定这根射线与图像存在交点
			//求解射线与图像相交的m，n的范围(m,n为图像的位置索引）
			int intersection_m_min = 0, intersection_n_min = 0;
			int intersection_m_max = image.rows, intersection_n_max = image.cols;
			Calculate_range_of_index_to_image(i, Trm, Rcv, image, min_alpha, max_alpha,
				intersection_n_min, intersection_n_max, intersection_m_min, intersection_m_max);

			//到这里，获得了行m的实际穿过范围，和列n的实际穿过范围
			//接下来求解所有的alpha_x和所有的alpha_y
			//与列标号intersection_n_min，到列标号intersection_n_max的所有交点的alpha值
			//等价于与直线x=intersection_n_min-image.cols/2，到与直线x=intersection_n_max-image.cols/2的所有交点的alpha值
			int sizeof_alpha;
			double* alpha;
			if ((theta_degree%90) != 0)
			{
				//求解alpha_x
				double* alpha_x;
				int sizeof_alpha_x = intersection_n_max - intersection_n_min + 1;
				alpha_x = new double[sizeof_alpha_x];
				Calculate_alpha_x(i, Trm, Rcv, image, intersection_n_min, sizeof_alpha_x, alpha_x);

				//求解alpha_y
				double* alpha_y;
				int sizeof_alpha_y = intersection_m_max - intersection_m_min + 1;
				//cout << "sizeof_alpha_y = " << sizeof_alpha_y << endl;
				alpha_y = new double[sizeof_alpha_y];
				Calculate_alpha_y(i, Trm, Rcv, image, intersection_m_min, sizeof_alpha_y, alpha_y);

				//接下来需要将alpha_x和alpha_y融合成一个数组，并且升序排列，得到alpha的求解结果
				sizeof_alpha = sizeof_alpha_x + sizeof_alpha_y;
				alpha = new double[sizeof_alpha];
				Calculate_alpha(sizeof_alpha_x, sizeof_alpha, alpha_x, alpha_y, alpha);
			}
			else
			{
				if ((theta_degree % 180) == 0)//平行于x轴，没有alpha_y
				{
					//求解alpha_x
					double* alpha_x;
					int sizeof_alpha_x = intersection_n_max - intersection_n_min + 1;
					alpha_x = new double[sizeof_alpha_x];
					Calculate_alpha_x(i, Trm, Rcv, image, intersection_n_min, sizeof_alpha_x, alpha_x);

					//求解alpha
					sizeof_alpha = sizeof_alpha_x;
					alpha = alpha_x;
					sort(alpha, alpha + sizeof_alpha);
				}
				else//平行于y轴，没有alpha_x
				{
					//求解alpha_y
					double* alpha_y;
					int sizeof_alpha_y = intersection_m_max - intersection_m_min + 1;
					alpha_y = new double[sizeof_alpha_y];
					Calculate_alpha_y(i, Trm, Rcv, image, intersection_m_min, sizeof_alpha_y, alpha_y);

					sizeof_alpha = sizeof_alpha_y;
					alpha = alpha_y;
					sort(alpha, alpha + sizeof_alpha);
				}
			}
			
			
			//接下来需要取邻近的alpha值（已升序排列），将长度（alpha值差 乘 射线总长度）与对应灰度值相乘
			//对应的数组索引通过mid_alpha = 0.5*[alpha(k)+alpha(k-1)]计算
			projection_1d[i] = Calculate_projection_oneValue(i, Trm, Rcv, image, sizeof_alpha, alpha);			
		}
		//将一维投影写入radon矩阵的对应列
		Project_copyto_2Darray(theta_degree,projection_1d,radon_Img);
	}
	return radon_Img;

}

void Radon_parallel_line_show(Mat image,int sensor_length, int current_degree, Mat& radon_Img)
{
	//第一步：投影
	//投影的尺寸计算
	//以图像的最长边为标准建立探测器长度，且为了在倾斜时能够覆盖所有图像像素，探测器的最小长度为图像边长的根号2倍
	//int sensor_length = Set_sensors_length_parallel_line(image);
	//Mat radon_Img = Mat::zeros(sensor_length, num_degrees, COLOR_BGR2GRAY);//存储投影后的图像矩阵，一列为一个角度

	//degree大循环，默认参数为循环180度，步长为1度
	//for (int theta_degree = 0; theta_degree < num_degrees; ++theta_degree)
		//for (int theta_degree = 0; theta_degree < 1; theta_degree+=90)
	//{
		/*
		构造射线和探测器
		由于输入图像已经转换为单通道，所以image.at<uchar>(i,j)中i为图像的行，j为图像的列
		设Transmittor(Trm)为发射源，receiver(Rcv)为接收器，Trm，Rcv均为大小为sensor_length的数组。
		0度表示从左向右的射线，90度表示从下向上的射线
		以图像中心为原点建立直角坐标系，则发射源和探测器的中点永远在圆x^2+y^2 = (1/4)sensor_length^2上
		theta为射线方向与x轴正方向的夹角
		（根据45度时探测器与图像中心距离至少为1/2图像斜角长度，探测器长度=图像斜角长度计算而得，中点的轨迹的半径=1/2探测器长度）
		*/
		Point2d* Trm;
		Point2d* Rcv;
		Trm = new Point2d[sensor_length];
		Rcv = new Point2d[sensor_length];
		Set_Transmittors_and_Sensors_parallel_line(current_degree, sensor_length, Trm, Rcv);//构造射线和探测器的位置

		/*    正式开始投影！！    */
		double* projection_1d;
		projection_1d = new double[sensor_length]; //记录特定角度下的各个探测器接收到的一系列投影

		//i小循环，特定某个角度theta_degree下的的各个发射器（or探测器）
		for (int i = 0; i < sensor_length; ++i)
			//for (int i = 55; i < 56; ++i)
		{
			//根据Trm[i]和Rcv[i]的坐标，求出对应Rcv[i]接收到的投影值projection_1d[i]
			//初始化第i个投影值，以及第i个发射器和接收器对应的alpha比例的最大和最小值
			//alpha的含义参考文献：Siddon, Robert L. "Fast calculation of the exact radiological path for a three‐dimensional CT array." Medical physics 12.2 (1985): 252-255.
			projection_1d[i] = 0.0;

			bool flag_intersection = 1;

			//计算min_alpha 和max_alpha
			double min_alpha = 0.0, max_alpha = 1.0;
			Calculate_range_of_alpha(i, Trm, Rcv, image, flag_intersection, min_alpha, max_alpha);
			if (flag_intersection == 0) { continue; }//没有交点，或者仅有一个点的交点（而不是线段），投影为0

			//代码运行到这里，表示已经获得了min_alpha，max_alpha，并且确定这根射线与图像存在交点
			//求解射线与图像相交的m，n的范围(m,n为图像的位置索引）
			int intersection_m_min = 0, intersection_n_min = 0;
			int intersection_m_max = image.rows, intersection_n_max = image.cols;
			Calculate_range_of_index_to_image(i, Trm, Rcv, image, min_alpha, max_alpha,
				intersection_n_min, intersection_n_max, intersection_m_min, intersection_m_max);

			//到这里，获得了行m的实际穿过范围，和列n的实际穿过范围
			//接下来求解所有的alpha_x和所有的alpha_y
			//与列标号intersection_n_min，到列标号intersection_n_max的所有交点的alpha值
			//等价于与直线x=intersection_n_min-image.cols/2，到与直线x=intersection_n_max-image.cols/2的所有交点的alpha值
			int sizeof_alpha;
			double* alpha;
			if ((current_degree % 90) != 0)
			{
				//求解alpha_x
				double* alpha_x;
				int sizeof_alpha_x = intersection_n_max - intersection_n_min + 1;
				alpha_x = new double[sizeof_alpha_x];
				Calculate_alpha_x(i, Trm, Rcv, image, intersection_n_min, sizeof_alpha_x, alpha_x);

				//求解alpha_y
				double* alpha_y;
				int sizeof_alpha_y = intersection_m_max - intersection_m_min + 1;
				//cout << "sizeof_alpha_y = " << sizeof_alpha_y << endl;
				alpha_y = new double[sizeof_alpha_y];
				Calculate_alpha_y(i, Trm, Rcv, image, intersection_m_min, sizeof_alpha_y, alpha_y);

				//接下来需要将alpha_x和alpha_y融合成一个数组，并且升序排列，得到alpha的求解结果
				sizeof_alpha = sizeof_alpha_x + sizeof_alpha_y;
				alpha = new double[sizeof_alpha];
				Calculate_alpha(sizeof_alpha_x, sizeof_alpha, alpha_x, alpha_y, alpha);
			}
			else
			{
				if ((current_degree % 180) == 0)//平行于x轴，没有alpha_y
				{
					//求解alpha_x
					double* alpha_x;
					int sizeof_alpha_x = intersection_n_max - intersection_n_min + 1;
					alpha_x = new double[sizeof_alpha_x];
					Calculate_alpha_x(i, Trm, Rcv, image, intersection_n_min, sizeof_alpha_x, alpha_x);

					//求解alpha
					sizeof_alpha = sizeof_alpha_x;
					alpha = alpha_x;
					sort(alpha, alpha + sizeof_alpha);
				}
				else//平行于y轴，没有alpha_x
				{
					//求解alpha_y
					double* alpha_y;
					int sizeof_alpha_y = intersection_m_max - intersection_m_min + 1;
					alpha_y = new double[sizeof_alpha_y];
					Calculate_alpha_y(i, Trm, Rcv, image, intersection_m_min, sizeof_alpha_y, alpha_y);

					sizeof_alpha = sizeof_alpha_y;
					alpha = alpha_y;
					sort(alpha, alpha + sizeof_alpha);
				}
			}


			//接下来需要取邻近的alpha值（已升序排列），将长度（alpha值差 乘 射线总长度）与对应灰度值相乘
			//对应的数组索引通过mid_alpha = 0.5*[alpha(k)+alpha(k-1)]计算
			projection_1d[i] = Calculate_projection_oneValue(i, Trm, Rcv, image, sizeof_alpha, alpha);
		}
		//将一维投影写入radon矩阵的对应列
		Project_copyto_2Darray(current_degree, projection_1d, radon_Img);
	//}
}









/** @brief Calculate the distance from a point to a line

@param Point:The point
@param A1,A2:The points on the line,represent the line.
@return  distance = the distance from the point to the line.
 */
static double Cal_distance_point_to_line(Point2d Point,Point2d A1,Point2d A2)
{
	//已知点为A1（x1，y1），A2（x2，y2），则这两点连成的直线为：(y2-y1)x-(x2-x1)y+(x2-x1)y1-(y2-y1)x1 = 0
	double para_A =A2.y - A1.y;
	double para_B = -(A2.x - A1.x);
	double para_C = (A2.x - A1.x) * A1.y - (A2.y - A1.y) * A1.x;
	/*cout << "Trm[0] = " << A1 << endl;
	cout << "Rcv[0] = " << A2 << endl;
	cout << "A = " << para_A << ", B = " << para_B << ", C = " << para_C << endl;*/
	double distance = abs(para_A * Point.x + para_B * Point.y + para_C)
		/ sqrt(para_A * para_A + para_B * para_B);

	return distance;
}


/** @brief Calculate the back-projection gray value corresponding to the pixel point by interpolation,
according to the [distance] from the center of the pixel point 
to the connection line between the 0th emitter and the detector,
as well as the [projection data].
(Under certain degree)

@param Point:distance: the distance from the pixel to the line between the 0th emitter and the detector
@param radon_Img:projection data
@param sensor_length:restrict the maximum of the index when searching in the projection data
@param theta_degree:current theta(0-179)
@return  grayscale = the grayscale that will be added to the pixel at the angle of theta_degree.
 */
static double Get_grayscale_iradon(double distance, const Mat radon_Img,
	const int sensor_length,const int theta_degree)
{
	double grayscale = 0;
	int k = floor(distance);
	if (k < 0)
	{
		cout << "Error!distance out of range(<0)." << endl;
		grayscale = radon_Img.at<double>(0, theta_degree);
		return grayscale;
	}

	//由于取整的问题，会有极少数的越界，或者+1之后越界，直接用最后一个数据即可
	if (k + 1 > sensor_length - 1)
	{
		cout << "Error!distance out of range(>max)." << endl;
		grayscale = radon_Img.at<double>(sensor_length - 1, theta_degree);
		return grayscale;
	}
	/*cout << "k = " << k << endl;
	cout << "distance - double(k) = " << distance - double(k) << endl;*/

	//一次插值
	grayscale = (distance - double(k)) * radon_Img.at<double>(k, theta_degree)
		+ (double(k) + 1 - distance) * radon_Img.at<double>(k + 1, theta_degree);

	return grayscale;
}

Mat iRadon_parallel_line(Mat radon_Img)
{
	cout << endl << "iRadon_parallel_line processing..." << endl;
	int sensor_length = radon_Img.rows;
	int num_degrees = radon_Img.cols;
	int img_rows = ceil(double(sensor_length) / sqrt(2));
	int img_cols = img_rows;
	Mat output_Img = Mat::zeros(img_rows, img_rows, COLOR_BGR2GRAY);

	cout << "sensor_length = " << sensor_length << ",num_degrees = " << num_degrees << endl;
	cout << "img_rows = " << img_rows << endl;

	//像素驱动方法
	for (int theta_degree = 0; theta_degree < num_degrees; ++theta_degree)
	//for (int theta_degree = 3; theta_degree < 4; ++theta_degree)
	{
		//float theta = float(theta_degree) * pi / 180;

		/*	构造射线和探测器(具体参考投影过程）	*/
		Point2d* Trm;
		Point2d* Rcv;
		Trm = new Point2d[sensor_length];
		Rcv = new Point2d[sensor_length];
		Set_Transmittors_and_Sensors_parallel_line(theta_degree, sensor_length, Trm, Rcv);//构造射线和探测器的位置

		for (int i = 0; i < img_rows; ++i)
		//for (int i = 100; i < 101; ++i)
		{
			for (int j = 0; j < img_cols; ++j)
			//for (int j = 100; j < 101; ++j)
			{
				Point2d center_point;
				center_point.x = (double(j) - double(img_cols) / 2);
				center_point.y = (-double(i) + double(img_rows) / 2);

				//求该 像素的中心【点】center_point 到 Trm[0]与Rcv[0]连成的【直线】 的距离 
				//求出的距离基本对应Rcv的序号（Rcv各个探测器之间的排列距离为1个单位）
				double distance = Cal_distance_point_to_line(center_point, Trm[0], Rcv[0]);
				//cout << "distance = " << distance << endl;
				//根据点到直线的距离，推算出所需要的投影对应的Trm[k]的序号k，据此插值返回对应反投影的灰度值
				double grayscale = Get_grayscale_iradon(distance, radon_Img, sensor_length, theta_degree);
				//cout << "grayscale = " << grayscale << endl;
				output_Img.at<double>(i, j) += grayscale;
			}

		}
	}
	
	//除以反投影的次数，使像素值在正常范围内
	for (int i = 0; i < img_rows; ++i)
	{
		for (int j = 0; j < img_cols; ++j)
		{
			output_Img.at<double>(i, j) /= num_degrees;
		}
	}

	cout << "iRadon_parallel_line finished." << endl;
	return output_Img;
}


void iRadon_parallel_line_show(Mat radon_Img, int current_degree, Mat& output_Img)
{
	int sensor_length = radon_Img.rows;
	int num_degrees = radon_Img.cols;
	int img_rows = ceil(double(sensor_length) / sqrt(2));
	int img_cols = img_rows;

	//像素驱动方法
	int theta_degree = current_degree;
	/*	构造射线和探测器(具体参考投影过程）	*/
	Point2d* Trm;
	Point2d* Rcv;
	Trm = new Point2d[sensor_length];
	Rcv = new Point2d[sensor_length];
	Set_Transmittors_and_Sensors_parallel_line(theta_degree, sensor_length, Trm, Rcv);//构造射线和探测器的位置

	for (int i = 0; i < img_rows; ++i)
		//for (int i = 100; i < 101; ++i)
	{
		for (int j = 0; j < img_cols; ++j)
			//for (int j = 100; j < 101; ++j)
		{
			Point2d center_point;
			center_point.x = (double(j) - double(img_cols) / 2);
			center_point.y = (-double(i) + double(img_rows) / 2);

			//求该 像素的中心【点】center_point 到 Trm[0]与Rcv[0]连成的【直线】 的距离 
			//求出的距离基本对应Rcv的序号（Rcv各个探测器之间的排列距离为1个单位）
			double distance = Cal_distance_point_to_line(center_point, Trm[0], Rcv[0]);
			//cout << "distance = " << distance << endl;
			//根据点到直线的距离，推算出所需要的投影对应的Trm[k]的序号k，据此插值返回对应反投影的灰度值
			double grayscale = Get_grayscale_iradon(distance, radon_Img, sensor_length, theta_degree);
			//cout << "grayscale = " << grayscale << endl;
			output_Img.at<double>(i, j) += grayscale;
		}

	}
}


//void iRadon_parallel_line_show(Mat radon_Img, int current_degree, Mat& output_Img)
//{
//	cout << endl << "iRadon_parallel_line processing..." << endl;
//	int sensor_length = radon_Img.rows;
//	int num_degrees = radon_Img.cols;
//	int img_rows = ceil(double(sensor_length) / sqrt(2));
//	int img_cols = img_rows;
//	output_Img = Mat::zeros(img_rows, img_rows, COLOR_BGR2GRAY);
//
//	cout << "sensor_length = " << sensor_length << ",num_degrees = " << num_degrees << endl;
//	cout << "img_rows = " << img_rows << endl;
//
//	//像素驱动方法
//	//for (int theta_degree = 0; theta_degree < num_degrees; ++theta_degree)
//	for (int theta_degree = current_degree; theta_degree < current_degree + 1; ++theta_degree)
//		//for (int theta_degree = 3; theta_degree < 4; ++theta_degree)
//	{
//		//float theta = float(theta_degree) * pi / 180;
//
//		/*	构造射线和探测器(具体参考投影过程）	*/
//		Point2d* Trm;
//		Point2d* Rcv;
//		Trm = new Point2d[sensor_length];
//		Rcv = new Point2d[sensor_length];
//		Set_Transmittors_and_Sensors_parallel_line(theta_degree, sensor_length, Trm, Rcv);//构造射线和探测器的位置
//
//		for (int i = 0; i < img_rows; ++i)
//			//for (int i = 100; i < 101; ++i)
//		{
//			for (int j = 0; j < img_cols; ++j)
//				//for (int j = 100; j < 101; ++j)
//			{
//				Point2d center_point;
//				center_point.x = (double(j) - double(img_cols) / 2);
//				center_point.y = (-double(i) + double(img_rows) / 2);
//
//				//求该 像素的中心【点】center_point 到 Trm[0]与Rcv[0]连成的【直线】 的距离 
//				//求出的距离基本对应Rcv的序号（Rcv各个探测器之间的排列距离为1个单位）
//				double distance = Cal_distance_point_to_line(center_point, Trm[0], Rcv[0]);
//				//cout << "distance = " << distance << endl;
//				//根据点到直线的距离，推算出所需要的投影对应的Trm[k]的序号k，据此插值返回对应反投影的灰度值
//				double grayscale = Get_grayscale_iradon(distance, radon_Img, sensor_length, theta_degree);
//				//cout << "grayscale = " << grayscale << endl;
//				output_Img.at<double>(i, j) += grayscale;
//			}
//
//		}
//	}
//}
