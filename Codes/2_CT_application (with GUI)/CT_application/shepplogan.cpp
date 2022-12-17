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
	int cols = inputImg.cols;//ʵ���ϴӹ����Դͷ������rows�����cols
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
		{
			double m = center_point.x * (rows / 2);//ƽ�Ʋ���m
			double n = center_point.y * (rows / 2);//ƽ�Ʋ���n
			double x = j - double(cols) / 2;//��ǰ������xy����ϵ�е�λ��-������
			double y = -i + double(rows) / 2;//��ǰ������xy����ϵ�е�λ��-������
			double a = major_axis * (rows / 2);//��Բ����a
			double b = minor_axis * (rows / 2);//��Բ����b
			double alpha = rotation_degree * pi / 180;	//��Բ��ת�Ƕ�alpha
			double grayscale = 255 * refractive_index;

			//�ο���ƽ�ƣ���ת�����Բ��ʽ https://blog.csdn.net/qq_41685265/article/details/104267256
			if (pow((x - m) * cos(alpha) + (y - n) * sin(alpha), 2) / (a * a)
				+ pow((m - x) * sin(alpha) + (y - n) * cos(alpha), 2) / (b * b) < 1)
			{
				inputImg.at<double>(i, j) += grayscale;//�ԻҶ�ֵΪ255Ϊ��׼
			}
		}
	}
	return inputImg;
}


//����������shepploganͼ�񣨽���ͼƬչʾ�ã����߶ԱȻ��ַ���ʱ����Ϊ����ͼƬ�ã�
Mat Create_shepplogan(int rows)
{
	int cols = rows;
	Mat Img = Mat::zeros(rows, rows, COLOR_BGR2GRAY);
	
	//�����ο���Avinash C. Kak, and Malcolm Slaney, 
	//��Principles of Computerized Tomographic Imaging��, IEEE Press, 1999. chapter3
	//Algorithms for Reconstruction with Nondiffracting Sources.
	//ע��Ϊ��ʹͼ��Աȶȸ����ԣ��Ը���Բ�ĻҶ�ֵ�������޸�
	Img = Add_ellipse(Img, Point2d(0, 0), 0.92, 0.69, 90, 2.0);//����Ȧ
	Img = Add_ellipse(Img, Point2d(0, -0.0184), 0.874, 0.6624, 90, -1.48);//����Ȧ
	Img = Add_ellipse(Img, Point2d(0.22, 0), 0.31, 0.11, 72, -0.2);//�Ҵ�Ȧ
	Img = Add_ellipse(Img, Point2d(-0.22, 0), 0.41, 0.16, 108, -0.2);//���Ȧ

	Img = Add_ellipse(Img, Point2d(0, 0.35), 0.25, 0.21, 90, 0.25);//���ϴ��Ȧ
	Img = Add_ellipse(Img, Point2d(0, 0.1), 0.046, 0.046, 0, 0.25);//�м������СȦ
	Img = Add_ellipse(Img, Point2d(0, -0.1), 0.046, 0.046, 0, 0.25);//�м������СȦ

	Img = Add_ellipse(Img, Point2d(-0.08, -0.605), 0.046, 0.023, 0, 0.25);//������ߵ�СȦ
	Img = Add_ellipse(Img, Point2d(0, -0.605), 0.023, 0.023, 0, 0.25);//�����м��СȦ
	Img = Add_ellipse(Img, Point2d(0.06, -0.605), 0.046, 0.023, 90, 0.25);//�����ұߵ�СȦ

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
		double m = center_point.x * (rows / 2);//ƽ�Ʋ���m
		double n = center_point.y * (rows / 2);//ƽ�Ʋ���n
		double a = major_axis * (rows / 2);//��Բ����a
		double b = minor_axis * (rows / 2);//��Բ����b
		double alpha = rotation_degree * pi / 180;//��Բ��ת�Ƕ�alpha

		double theta = (theta_degree - double(90)) * pi / 180;//ͶӰ�Ƕ�theta����ʽ��������Ϊ�����򣬳��������н�ģ����������Ϊ�����������90�㣩

		//��ʽ�ο���Avinash C. Kak, and Malcolm Slaney, 
		//��Principles of Computerized Tomographic Imaging��, IEEE Press, 1999. Chapter3
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
			if ((m < 0))//�ڶ��������ޣ�����m��nֻ��Ӱ��gamma��һ���Ͷ��������ڲ��޷����֣�����û���������ϸ��...��
			{
				gamma = atan(float(n / m)) + pi;
			}
			else  //��һ��������
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


//������ֱ�Ӽ���shepplogan��ͶӰ
Mat Calc_shepplogan_radon(int rows)
{
	int cols = rows;
	int num_degrees = 180;
	int rows_radon = rows * sqrt(2);
	Mat radon_Img = Mat::zeros(rows_radon, num_degrees, COLOR_BGR2GRAY);

	radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0, 0), 0.92, 0.69, 90, 2.0);//����Ȧ
	radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0, -0.0184), 0.874, 0.6624, 90, -1.48);//����Ȧ
	radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0.22, 0), 0.31, 0.11, 72, -0.2);//�Ҵ�Ȧ
	radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(-0.22, 0), 0.41, 0.16, 108, -0.2);//���Ȧ

	radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0, 0.35), 0.25, 0.21, 90, 0.25);//���ϴ��Ȧ
	radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0, 0.1), 0.046, 0.046, 0, 0.25);//�м������СȦ
	radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0, -0.1), 0.046, 0.046, 0, 0.25);//�м������СȦ

	radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(-0.08, -0.605), 0.046, 0.023, 0, 0.25);//������ߵ�СȦ
	radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0, -0.605), 0.023, 0.023, 0, 0.25);//�����м��СȦ
	radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0.06, -0.605), 0.046, 0.023, 90, 0.25);//�����ұߵ�СȦ

	return radon_Img;
}

void Calc_shepplogan_radon_show(int rows, int index, Mat& radon_Img)
{
	int num_degrees = 180;
	int rows_radon = rows * sqrt(2);
	
	switch(index)
	{
		case 0:radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0, 0), 0.92, 0.69, 90, 2.0);//����Ȧ
			break;
		case 1:radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0, -0.0184), 0.874, 0.6624, 90, -1.48);//����Ȧ
			break;
		case 2:radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0.22, 0), 0.31, 0.11, 72, -0.2);//�Ҵ�Ȧ
			break;
		case 3:radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(-0.22, 0), 0.41, 0.16, 108, -0.2);//���Ȧ
			break;
		case 4:radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0, 0.35), 0.25, 0.21, 90, 0.25);//���ϴ��Ȧ
			break;
		case 5:radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0, 0.1), 0.046, 0.046, 0, 0.25);//�м������СȦ
			break;
		case 6:radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0, -0.1), 0.046, 0.046, 0, 0.25);//�м������СȦ
			break;
		case 7:radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(-0.08, -0.605), 0.046, 0.023, 0, 0.25);//������ߵ�СȦ
			break;
		case 8:radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0, -0.605), 0.023, 0.023, 0, 0.25);//�����м��СȦ
			break;
		case 9:radon_Img = Add_ellipse_radon(radon_Img, rows, Point2d(0.06, -0.605), 0.046, 0.023, 90, 0.25);//�����ұߵ�СȦ
			break;
		
	}
}
