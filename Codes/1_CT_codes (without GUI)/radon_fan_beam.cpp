#include <iostream>
#include<opencv2/opencv.hpp>

#include "radon_fan_beam.h"
#include "basic_functions.h"
#include "basic_radon_functions.h"

//#include <opencv2/highgui/highgui_c.h> 

using namespace std;
using namespace cv;




//自己写的扇形束投影函数
/** @brief Fan beam radon

@param image Input 8-bit 1-channel image.
@param num_rotate_degree : The number of degrees of the transmittors
@param num_detector_degree : The number of degrees of receivers(The maximum angle between the transmitter and different receivers)
@return  radon_Img Output 64-bit 1-channel image.
		Size:rows = number of receivers = num_detector_degree( default = 60)
			 cols = number of degrees = num_rotate_degree( default = 360)
 */
//Mat Radon_fan_beam(Mat image, int num_rotate_degree, int num_detector_degree)
//{
//	cout << "fan beam radon processing......" << endl;
//	Mat radon_Img = Mat::zeros(num_detector_degree, num_rotate_degree, COLOR_BGR2GRAY);
//	double R_circle = max(image.rows, image.cols) * sqrt(2) / 2 + 2;  //+2保证在离散近似下，圆仍然能在图像外围
//	
//	for (int current_rotate_degree = 0; current_rotate_degree < num_rotate_degree; ++current_rotate_degree)
//	//for (int current_rotate_degree = 180; current_rotate_degree < 181; ++current_rotate_degree)
//	{
//		/*
//		【发射器和接收器位置的确定】
//		current_rotate_degree是中心射束和x轴正方向的夹角，
//		Rcv以顺时针方向排列，原点到中心接收器的射线，即为中心射束（与x轴正方向夹角为current_rotate_degree）
//
//		Trm矩阵中各个元素相等（扇形束只有一个发射器），为了和平行束的后续计算兼容，选择了数组存储；
//		Rcv共有num_detector_degree个元素，射束顺时针旋转1°，Rcv的索引加一。
//		因此num_detector_degree也是 Trm与Rcv[0]连线 与 Trm与Rcv[num_detector_degree-1]连线 之间的夹角。
//
//		Trm 和 Rcv都在圆心为原点，半径为R_circle的圆上
//		
//		*/
//		Point2d* Trm;
//		Point2d* Rcv;
//		Trm = new Point2d[num_detector_degree];
//		Rcv = new Point2d[num_detector_degree];
//		Set_Transmittors_and_Sensors_fan_beam(current_rotate_degree, num_detector_degree, R_circle, Trm, Rcv);//构造射线和探测器的位置
//
//		/*    正式开始投影！！    */
//		double* projection_1d;
//		projection_1d = new double[num_detector_degree]; //记录发射器特定角度（位置）下的各个探测器接收到的一系列投影
//
//		//i小循环，特定某个角度current_rotate_degree（发射器固定）下的的各个探测器
//		for (int i = 0; i < num_detector_degree; ++i)
//		//for (int i = 29; i < 30; ++i)
//		{
//			//根据Trm[i]和Rcv[i]的坐标，求出对应Rcv[i]接收到的投影值projection_1d[i]
//			//初始化第i个投影值，以及第i个发射器和接收器对应的alpha比例的最大和最小值
//			//alpha的含义参考文献：Siddon, Robert L. "Fast calculation of the exact radiological path for a three‐dimensional CT array." Medical physics 12.2 (1985): 252-255.
//			projection_1d[i] = 0.0;
//
//			bool flag_intersection = 1;
//
//			//计算min_alpha 和max_alpha
//			double min_alpha = 0.0, max_alpha = 1.0;
//			Calculate_range_of_alpha(i, Trm, Rcv, image, flag_intersection, min_alpha, max_alpha);
//			if (flag_intersection == 0) { continue; }//没有交点，或者仅有一个点的交点（而不是线段），投影为0
//
//			//代码运行到这里，表示已经获得了min_alpha，max_alpha，并且确定这根射线与图像存在交点
//			//求解射线与图像相交的m，n的范围(m,n为图像的位置索引）
//			int intersection_m_min = 0, intersection_n_min = 0;
//			int intersection_m_max = image.rows, intersection_n_max = image.cols;
//			Calculate_range_of_index_to_image(i, Trm, Rcv, image, min_alpha, max_alpha,
//				intersection_n_min, intersection_n_max, intersection_m_min, intersection_m_max);
//
//			//到这里，获得了行m的实际穿过范围，和列n的实际穿过范围
//			//接下来求解所有的alpha_x和所有的alpha_y
//			//与列标号intersection_n_min，到列标号intersection_n_max的所有交点的alpha值
//			//等价于与直线x=intersection_n_min-image.cols/2，到与直线x=intersection_n_max-image.cols/2的所有交点的alpha值
//			int sizeof_alpha;
//			double* alpha;
//			int theta_degree_ray = current_rotate_degree + num_detector_degree - 2 * i - 1;
//			if ((theta_degree_ray % 90) != 0)
//			{
//				//求解alpha_x
//				double* alpha_x;
//				int sizeof_alpha_x = intersection_n_max - intersection_n_min + 1;
//				alpha_x = new double[sizeof_alpha_x];
//				Calculate_alpha_x(i, Trm, Rcv, image, intersection_n_min, sizeof_alpha_x, alpha_x);
//
//				//求解alpha_y
//				double* alpha_y;
//				int sizeof_alpha_y = intersection_m_max - intersection_m_min + 1;
//				//cout << "sizeof_alpha_y = " << sizeof_alpha_y << endl;
//				alpha_y = new double[sizeof_alpha_y];
//				Calculate_alpha_y(i, Trm, Rcv, image, intersection_m_min, sizeof_alpha_y, alpha_y);
//
//				//接下来需要将alpha_x和alpha_y融合成一个数组，并且升序排列，得到alpha的求解结果
//				sizeof_alpha = sizeof_alpha_x + sizeof_alpha_y;
//				alpha = new double[sizeof_alpha];
//				Calculate_alpha(sizeof_alpha_x, sizeof_alpha, alpha_x, alpha_y, alpha);
//			}
//			else
//			{
//				if ((theta_degree_ray % 180) == 0)//平行于x轴，没有alpha_y
//				{
//					//求解alpha_x
//					double* alpha_x;
//					int sizeof_alpha_x = intersection_n_max - intersection_n_min + 1;
//					alpha_x = new double[sizeof_alpha_x];
//					Calculate_alpha_x(i, Trm, Rcv, image, intersection_n_min, sizeof_alpha_x, alpha_x);
//
//					//求解alpha
//					sizeof_alpha = sizeof_alpha_x;
//					alpha = alpha_x;
//					sort(alpha, alpha + sizeof_alpha);
//				}
//				else//平行于y轴，没有alpha_x
//				{
//					//求解alpha_y
//					double* alpha_y;
//					int sizeof_alpha_y = intersection_m_max - intersection_m_min + 1;
//					alpha_y = new double[sizeof_alpha_y];
//					Calculate_alpha_y(i, Trm, Rcv, image, intersection_m_min, sizeof_alpha_y, alpha_y);
//
//					sizeof_alpha = sizeof_alpha_y;
//					alpha = alpha_y;
//					sort(alpha, alpha + sizeof_alpha);
//				}
//			}
//
//
//			//接下来需要取邻近的alpha值（已升序排列），将长度（alpha值差 乘 射线总长度）与对应灰度值相乘
//			//对应的数组索引通过mid_alpha = 0.5*[alpha(k)+alpha(k-1)]计算
//			projection_1d[i] = Calculate_projection_oneValue(i, Trm, Rcv, image, sizeof_alpha, alpha);
//		}
//		//将一维投影写入radon矩阵的对应列
//		Project_copyto_2Darray(current_rotate_degree, projection_1d, radon_Img);
//
//	}
//
//	return radon_Img;
//}

Mat Radon_fan_beam(Mat image, int num_rotate_degree, int num_detector_degree)
{
	cout << "fan beam radon processing......" << endl;
	Mat radon_Img = Mat::zeros(num_detector_degree, num_rotate_degree, COLOR_BGR2GRAY);
	double R_circle = max(image.rows, image.cols) * sqrt(2) / 2 + 2;  //+2保证在离散近似下，圆仍然能在图像外围

	for (int current_rotate_degree = 0; current_rotate_degree < num_rotate_degree; ++current_rotate_degree)
		//for (int current_rotate_degree = 180; current_rotate_degree < 181; ++current_rotate_degree)
	{
		/*
		【发射器和接收器位置的确定】
		current_rotate_degree是中心射束和x轴正方向的夹角，
		Rcv以顺时针方向排列，原点到中心接收器的射线，即为中心射束（与x轴正方向夹角为current_rotate_degree）

		Trm矩阵中各个元素相等（扇形束只有一个发射器），为了和平行束的后续计算兼容，选择了数组存储；
		Rcv共有num_detector_degree个元素，射束顺时针旋转1°，Rcv的索引加一。
		因此num_detector_degree也是 Trm与Rcv[0]连线 与 Trm与Rcv[num_detector_degree-1]连线 之间的夹角。

		Trm 和 Rcv都在圆心为原点，半径为R_circle的圆上

		*/
		Point2d* Trm;
		Point2d* Rcv;
		Trm = new Point2d[num_detector_degree];
		Rcv = new Point2d[num_detector_degree];
		Set_Transmittors_and_Sensors_fan_beam(current_rotate_degree, num_detector_degree, R_circle, Trm, Rcv);//构造射线和探测器的位置

		/*    正式开始投影！！    */
		double* projection_1d;
		projection_1d = new double[num_detector_degree]; //记录发射器特定角度（位置）下的各个探测器接收到的一系列投影

		//i小循环，特定某个角度current_rotate_degree（发射器固定）下的的各个探测器
		for (int i = 0; i < num_detector_degree; ++i)
			//for (int i = 29; i < 30; ++i)
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
			int theta_degree_ray = current_rotate_degree + num_detector_degree - 2 * i - 1;
			if ((theta_degree_ray % 90) != 0)
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
				if ((theta_degree_ray % 180) == 0)//平行于x轴，没有alpha_y
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
		Project_copyto_2Darray(current_rotate_degree, projection_1d, radon_Img);

	}

	return radon_Img;
}


static void Set_Transmittors_and_Sensors_fan_beam(const int current_rotate_degree, const int num_detector_degree,
	const double R_circle, Point2d* Trm, Point2d* Rcv)
{
	//传入数据均为角度值，需要转换为弧度制！！
	double theta_rad = double(current_rotate_degree) * pi / 180;

	//设置发射器位置（特定角度下，扇形束发射器的位置为一个点）
	for (int i = 0; i < num_detector_degree; ++i)
	{
		Trm[i].x = R_circle * cos(theta_rad + pi);
		Trm[i].y = R_circle * sin(theta_rad + pi);
	}
	//cout << "Trm = " << Trm[0] << endl;

	//设置接收器位置
	for (int i = 0; i < num_detector_degree; ++i)
	{
		double tmp = double(current_rotate_degree) + double(num_detector_degree) - 2 * double(i) - 1;
		tmp = tmp * pi / 180;
		Rcv[i].x = R_circle * cos(tmp);
		Rcv[i].y = R_circle * sin(tmp);
		//cout << "Rcv[" << i << "] = " << Rcv[i] << endl;
	}
}


