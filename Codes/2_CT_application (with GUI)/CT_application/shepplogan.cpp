#include <iostream>
#include<opencv2/opencv.hpp>

#include "shepplogan.h"

using namespace cv;
#define  pi  3.14159265358979323846  //!< const value,pi


/** @brief add ellipse to the image
@param   inputImg: The image to be modified
@param   center_point: The center of the ellipse
@param   major_axis: The major_axis of the ellipse
@param   minor_axis: The minor_axis of the ellipse
@param   rotation_degree: The rotation degree of the ellipse
@param   refractive_index: The refractive_index of the ellipse( influence the grayscale of the image )
@return  outputImg:The output image( inputImg + ellipse )

 */
static Mat Add_ellipse(Mat inputImg, Point2d center_point, double major_axis, double minor_axis,
	int rotation_degree, double refractive_index)
{
	int rows = inputImg.rows;
	int cols = inputImg.cols;//实际上从构造的源头来看，rows恒等于cols
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
		{
			double m = center_point.x * (rows / 2);//平移参数m
			double n = center_point.y * (rows / 2);//平移参数n
			double x = j - double(cols) / 2;//当前像素在xy坐标系中的位置-横坐标
			double y = -i + double(rows) / 2;//当前像素在xy坐标系中的位置-纵坐标
			double a = major_axis * (rows / 2);//椭圆参数a
			double b = minor_axis * (rows / 2);//椭圆参数b
			double alpha = rotation_degree * pi / 180;	//椭圆旋转角度alpha
			double grayscale = 255 * refractive_index;

			//参考：平移，旋转后的椭圆公式 https://blog.csdn.net/qq_41685265/article/details/104267256
			if (pow((x - m) * cos(alpha) + (y - n) * sin(alpha), 2) / (a * a)
				+ pow((m - x) * sin(alpha) + (y - n) * cos(alpha), 2) / (b * b) < 1)
			{
				inputImg.at<double>(i, j) += grayscale;//以灰度值为255为基准
			}
		}
	}
	return inputImg;
}


//解析法计算shepplogan图像（仅供图片展示用，或者对比积分法的时候作为读入图片用）
Mat Create_shepplogan(int rows)
{
	int cols = rows;
	Mat Img = Mat::zeros(rows, rows, COLOR_BGR2GRAY);
	
	//参数参考：Avinash C. Kak, and Malcolm Slaney, 
	//“Principles of Computerized Tomographic Imaging”, IEEE Press, 1999. chapter3
	//Algorithms for Reconstruction with Nondiffracting Sources.
	//注：为了使图像对比度更明显，对各椭圆的灰度值进行了修改
	Img = Add_ellipse(Img, Point2d(0, 0), 0.92, 0.69, 90, 2.0);//主白圈
	Img = Add_ellipse(Img, Point2d(0, -0.0184), 0.874, 0.6624, 90, -1.48);//主灰圈
	Img = Add_ellipse(Img, Point2d(0.22, 0), 0.31, 0.11, 72, -0.2);//右大圈
	Img = Add_ellipse(Img, Point2d(-0.22, 0), 0.41, 0.16, 108, -0.2);//左大圈

	Img = Add_ellipse(Img, Point2d(0, 0.35), 0.25, 0.21, 90, 0.25);//中上大白圈
	Img = Add_ellipse(Img, Point2d(0, 0.1), 0.046, 0.046, 0, 0.25);//中间上面的小圈
	Img = Add_ellipse(Img, Point2d(0, -0.1), 0.046, 0.046, 0, 0.25);//中间下面的小圈

	Img = Add_ellipse(Img, Point2d(-0.08, -0.605), 0.046, 0.023, 0, 0.25);//下面左边的小圈
	Img = Add_ellipse(Img, Point2d(0, -0.605), 0.023, 0.023, 0, 0.25);//下面中间的小圈
	Img = Add_ellipse(Img, Point2d(0.06, -0.605), 0.046, 0.023, 90, 0.25);//下面右边的小圈

	return Img;
}


/** @brief add radon result of ellipse to the radon image
@param   radonImg: The image to be modified
@param   rows: The rows (and cols) of the original image (input by the user)
@param   center_point: The center of the ellipse
@param   major_axis: The major_axis of the ellipse
@param   minor_axis: The minor_axis of the ellipse
@param   rotation_degree: The rotation degree of the ellipse
@param   refractive_index: The refractive_index of the ellipse( influence the grayscale of the image )
@return  outputImg:The output image( radonImg + ellipse_radon )

 */
static Mat Add_ellipse_radon(Mat radonImg, int rows, Point2d center_point, double major_axis, double minor_axis,
	int rotation_degree, double refractive_index)
{
	int radon_rows = radonImg.rows;
	int num_degrees = radonImg.cols;

	for (int theta_degree = 0; theta_degree < 180; ++theta_degree)
	{
		double m = center_point.x * (rows / 2);//平移参数m
		double n = center_point.y * (rows / 2);//平移参数n
		double a = major_axis * (rows / 2);//椭圆参数a
		double b = minor_axis * (rows / 2);//椭圆参数b
		double alpha = rotation_degree * pi / 180;//椭圆旋转角度alpha

		double theta = (theta_degree - double(90)) * pi / 180;//投影角度theta（公式中以向上为正方向，程序中所有建模以射线向右为正方向，相差了90°）

		//公式参考：Avinash C. Kak, and Malcolm Slaney, 
		//“Principles of Computerized Tomographic Imaging”, IEEE Press, 1999. Chapter3
		double a_theta_2 = pow(a, 2) * pow(cos(theta - alpha), 2)
			+ pow(b, 2) * pow(sin(theta - alpha), 2);
		double s = sqrt(m * m + n * n);
		double gamma = 0;
		if (n == 0)
		{
			if (m >= 0)
				gamma = 0;
			else
				gamma = pi;
		}
		else
		{
			if ((m < 0))//第二、三象限（否则，m，n只能影响gamma，一三和二四象限内部无法区分，书中没有体现这个细节...）
			{
				gamma = atan(float(n / m)) + pi;
			}
			else  //第一、四象限
			{
				gamma = atan(float(n / m));
			}
		}

		for (int i = 0; i < radon_rows; ++i)
		{
			double t0 = i - 0.5 * radon_rows;
			double t = t0 - s * cos(gamma - theta);
			if (t * t <= a_theta_2)
			{
				radonImg.at<double>(i, theta_degree) += ((2 * double(255) * refractive_index * a * b) / a_theta_2)
					* sqrt(a_theta_2 - pow(t, 2));
			}	
		}
	}
	
	return radonImg;
}


//解析法直接计算shepplogan的投影
Mat Calc_shepplogan_radon(int rows)
{
	int cols = rows;
	int num_degrees = 180;
	int rows_radon = rows * sqrt(2);
	Mat radon_Img = Mat::zeros(rows_radon, num_degrees, COLOR_BGR2GRAY);

	radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0, 0), 0.92, 0.69, 90, 2.0);//主白圈
	radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0, -0.0184), 0.874, 0.6624, 90, -1.48);//主灰圈
	radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0.22, 0), 0.31, 0.11, 72, -0.2);//右大圈
	radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(-0.22, 0), 0.41, 0.16, 108, -0.2);//左大圈

	radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0, 0.35), 0.25, 0.21, 90, 0.25);//中上大白圈
	radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0, 0.1), 0.046, 0.046, 0, 0.25);//中间上面的小圈
	radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0, -0.1), 0.046, 0.046, 0, 0.25);//中间下面的小圈

	radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(-0.08, -0.605), 0.046, 0.023, 0, 0.25);//下面左边的小圈
	radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0, -0.605), 0.023, 0.023, 0, 0.25);//下面中间的小圈
	radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0.06, -0.605), 0.046, 0.023, 90, 0.25);//下面右边的小圈

	return radon_Img;
}

void Calc_shepplogan_radon_show(int rows, int index, Mat& radon_Img)
{
	int num_degrees = 180;
	int rows_radon = rows * sqrt(2);
	
	switch(index)
	{
		case 0:radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0, 0), 0.92, 0.69, 90, 2.0);//主白圈
			break;
		case 1:radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0, -0.0184), 0.874, 0.6624, 90, -1.48);//主灰圈
			break;
		case 2:radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0.22, 0), 0.31, 0.11, 72, -0.2);//右大圈
			break;
		case 3:radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(-0.22, 0), 0.41, 0.16, 108, -0.2);//左大圈
			break;
		case 4:radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0, 0.35), 0.25, 0.21, 90, 0.25);//中上大白圈
			break;
		case 5:radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0, 0.1), 0.046, 0.046, 0, 0.25);//中间上面的小圈
			break;
		case 6:radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0, -0.1), 0.046, 0.046, 0, 0.25);//中间下面的小圈
			break;
		case 7:radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(-0.08, -0.605), 0.046, 0.023, 0, 0.25);//下面左边的小圈
			break;
		case 8:radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0, -0.605), 0.023, 0.023, 0, 0.25);//下面中间的小圈
			break;
		case 9:radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0.06, -0.605), 0.046, 0.023, 90, 0.25);//下面右边的小圈
			break;
		
	}
}
