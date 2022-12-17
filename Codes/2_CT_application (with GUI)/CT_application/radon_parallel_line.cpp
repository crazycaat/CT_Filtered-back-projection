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
	//ͶӰ�ĳߴ����
	//��ͼ������Ϊ��׼����̽�������ȣ���Ϊ������бʱ�ܹ���������ͼ�����أ�̽��������С����Ϊͼ��߳��ĸ���2��
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

	//�������е�
	Point2d mid_Trm;
	mid_Trm.x = 0.5 * (double(sensor_length) - 2) * cos(theta + pi);//��ͼ����΢��һ�㣬��֤45��ʱȫ���ǣ�ȡ�����������⣩��ǰ���Ӧ�Ѿ�+2
	mid_Trm.y = 0.5 * (double(sensor_length) - 2) * sin(theta + pi);
	cout << "theta = " << theta_degree << " processing......" << endl;
	//cout << "mid_Trm = " << mid_Trm << endl;
	//������λ��ȷ����Trm[0]~Trm[sensor_length-1]��sensor_length��
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


	//�������е�
	Point2d mid_Rcv;
	mid_Rcv.x = 0.5 * double(sensor_length) * cos(theta);
	mid_Rcv.y = 0.5 * double(sensor_length) * sin(theta);
	//cout << "mid_Rcv = " << mid_Rcv << endl;
	//������λ��ȷ����Rcv[0]~Rcv[sensor_length-1]��sensor_length��
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
Siddon, Robert L. "Fast calculation of the exact radiological path for a three�\dimensional CT array."
Medical physics 12.2 (1985): 252-255.

@note
   -  the output image(radon_Img) has a narrow effective grayscale with a large grayscale range.
      Therefore, simple process of radon_Img before showing and storing the result is needed.@ref function-Convert_to_show_radon_parallel_line
 */
Mat Radon_parallel_line(Mat image, int num_degrees)
{
	//��һ����ͶӰ
	//ͶӰ�ĳߴ����
	//��ͼ������Ϊ��׼����̽�������ȣ���Ϊ������бʱ�ܹ���������ͼ�����أ�̽��������С����Ϊͼ��߳��ĸ���2��
	int sensor_length = Set_sensors_length_parallel_line(image);
	Mat radon_Img = Mat::zeros(sensor_length, num_degrees, COLOR_BGR2GRAY);//�洢ͶӰ���ͼ�����һ��Ϊһ���Ƕ�
		
	//degree��ѭ����Ĭ�ϲ���Ϊѭ��180�ȣ�����Ϊ1��
	for (int theta_degree = 0; theta_degree < num_degrees; ++theta_degree)
	//for (int theta_degree = 0; theta_degree < 1; theta_degree+=90)
	{
		/*
		�������ߺ�̽����
		��������ͼ���Ѿ�ת��Ϊ��ͨ��������image.at<uchar>(i,j)��iΪͼ����У�jΪͼ�����
		��Transmittor(Trm)Ϊ����Դ��receiver(Rcv)Ϊ��������Trm��Rcv��Ϊ��СΪsensor_length�����顣
		0�ȱ�ʾ�������ҵ����ߣ�90�ȱ�ʾ�������ϵ�����
		��ͼ������Ϊԭ�㽨��ֱ������ϵ������Դ��̽�������е���Զ��Բx^2+y^2 = (1/4)sensor_length^2��
		thetaΪ���߷�����x��������ļн�
		������45��ʱ̽������ͼ�����ľ�������Ϊ1/2ͼ��б�ǳ��ȣ�̽��������=ͼ��б�ǳ��ȼ�����ã��е�Ĺ켣�İ뾶=1/2̽�������ȣ�
		*/
		Point2d* Trm;
		Point2d* Rcv;
		Trm = new Point2d[sensor_length];
		Rcv = new Point2d[sensor_length];
		Set_Transmittors_and_Sensors_parallel_line(theta_degree, sensor_length, Trm, Rcv);//�������ߺ�̽������λ��

		/*    ��ʽ��ʼͶӰ����    */
		double* projection_1d;
		projection_1d = new double[sensor_length]; //��¼�ض��Ƕ��µĸ���̽�������յ���һϵ��ͶӰ

		//iСѭ�����ض�ĳ���Ƕ�theta_degree�µĵĸ�����������or̽������
		for (int i = 0; i < sensor_length; ++i)
		//for (int i = 55; i < 56; ++i)
		{
			//����Trm[i]��Rcv[i]�����꣬�����ӦRcv[i]���յ���ͶӰֵprojection_1d[i]
			//��ʼ����i��ͶӰֵ���Լ���i���������ͽ�������Ӧ��alpha������������Сֵ
			//alpha�ĺ���ο����ף�Siddon, Robert L. "Fast calculation of the exact radiological path for a three�\dimensional CT array." Medical physics 12.2 (1985): 252-255.
			projection_1d[i] = 0.0;
			
			bool flag_intersection = 1;
			
			//����min_alpha ��max_alpha
			double min_alpha = 0.0, max_alpha = 1.0;
			Calculate_range_of_alpha(i, Trm, Rcv, image, flag_intersection, min_alpha, max_alpha);
			if (flag_intersection==0){continue;}//û�н��㣬���߽���һ����Ľ��㣨�������߶Σ���ͶӰΪ0
			
			//�������е������ʾ�Ѿ������min_alpha��max_alpha������ȷ�����������ͼ����ڽ���
			//���������ͼ���ཻ��m��n�ķ�Χ(m,nΪͼ���λ��������
			int intersection_m_min = 0, intersection_n_min = 0;
			int intersection_m_max = image.rows, intersection_n_max = image.cols;
			Calculate_range_of_index_to_image(i, Trm, Rcv, image, min_alpha, max_alpha,
				intersection_n_min, intersection_n_max, intersection_m_min, intersection_m_max);

			//������������m��ʵ�ʴ�����Χ������n��ʵ�ʴ�����Χ
			//������������е�alpha_x�����е�alpha_y
			//���б��intersection_n_min�����б��intersection_n_max�����н����alphaֵ
			//�ȼ�����ֱ��x=intersection_n_min-image.cols/2������ֱ��x=intersection_n_max-image.cols/2�����н����alphaֵ
			int sizeof_alpha;
			double* alpha;
			if ((theta_degree%90) != 0)
			{
				//���alpha_x
				double* alpha_x;
				int sizeof_alpha_x = intersection_n_max - intersection_n_min + 1;
				alpha_x = new double[sizeof_alpha_x];
				Calculate_alpha_x(i, Trm, Rcv, image, intersection_n_min, sizeof_alpha_x, alpha_x);

				//���alpha_y
				double* alpha_y;
				int sizeof_alpha_y = intersection_m_max - intersection_m_min + 1;
				//cout << "sizeof_alpha_y = " << sizeof_alpha_y << endl;
				alpha_y = new double[sizeof_alpha_y];
				Calculate_alpha_y(i, Trm, Rcv, image, intersection_m_min, sizeof_alpha_y, alpha_y);

				//��������Ҫ��alpha_x��alpha_y�ںϳ�һ�����飬�����������У��õ�alpha�������
				sizeof_alpha = sizeof_alpha_x + sizeof_alpha_y;
				alpha = new double[sizeof_alpha];
				Calculate_alpha(sizeof_alpha_x, sizeof_alpha, alpha_x, alpha_y, alpha);
			}
			else
			{
				if ((theta_degree % 180) == 0)//ƽ����x�ᣬû��alpha_y
				{
					//���alpha_x
					double* alpha_x;
					int sizeof_alpha_x = intersection_n_max - intersection_n_min + 1;
					alpha_x = new double[sizeof_alpha_x];
					Calculate_alpha_x(i, Trm, Rcv, image, intersection_n_min, sizeof_alpha_x, alpha_x);

					//���alpha
					sizeof_alpha = sizeof_alpha_x;
					alpha = alpha_x;
					sort(alpha, alpha + sizeof_alpha);
				}
				else//ƽ����y�ᣬû��alpha_x
				{
					//���alpha_y
					double* alpha_y;
					int sizeof_alpha_y = intersection_m_max - intersection_m_min + 1;
					alpha_y = new double[sizeof_alpha_y];
					Calculate_alpha_y(i, Trm, Rcv, image, intersection_m_min, sizeof_alpha_y, alpha_y);

					sizeof_alpha = sizeof_alpha_y;
					alpha = alpha_y;
					sort(alpha, alpha + sizeof_alpha);
				}
			}
			
			
			//��������Ҫȡ�ڽ���alphaֵ�����������У��������ȣ�alphaֵ�� �� �����ܳ��ȣ����Ӧ�Ҷ�ֵ���
			//��Ӧ����������ͨ��mid_alpha = 0.5*[alpha(k)+alpha(k-1)]����
			projection_1d[i] = Calculate_projection_oneValue(i, Trm, Rcv, image, sizeof_alpha, alpha);			
		}
		//��һάͶӰд��radon����Ķ�Ӧ��
		Project_copyto_2Darray(theta_degree,projection_1d,radon_Img);
	}
	return radon_Img;

}

void Radon_parallel_line_show(Mat image,int sensor_length, int current_degree, Mat& radon_Img)
{
	//��һ����ͶӰ
	//ͶӰ�ĳߴ����
	//��ͼ������Ϊ��׼����̽�������ȣ���Ϊ������бʱ�ܹ���������ͼ�����أ�̽��������С����Ϊͼ��߳��ĸ���2��
	//int sensor_length = Set_sensors_length_parallel_line(image);
	//Mat radon_Img = Mat::zeros(sensor_length, num_degrees, COLOR_BGR2GRAY);//�洢ͶӰ���ͼ�����һ��Ϊһ���Ƕ�

	//degree��ѭ����Ĭ�ϲ���Ϊѭ��180�ȣ�����Ϊ1��
	//for (int theta_degree = 0; theta_degree < num_degrees; ++theta_degree)
		//for (int theta_degree = 0; theta_degree < 1; theta_degree+=90)
	//{
		/*
		�������ߺ�̽����
		��������ͼ���Ѿ�ת��Ϊ��ͨ��������image.at<uchar>(i,j)��iΪͼ����У�jΪͼ�����
		��Transmittor(Trm)Ϊ����Դ��receiver(Rcv)Ϊ��������Trm��Rcv��Ϊ��СΪsensor_length�����顣
		0�ȱ�ʾ�������ҵ����ߣ�90�ȱ�ʾ�������ϵ�����
		��ͼ������Ϊԭ�㽨��ֱ������ϵ������Դ��̽�������е���Զ��Բx^2+y^2 = (1/4)sensor_length^2��
		thetaΪ���߷�����x��������ļн�
		������45��ʱ̽������ͼ�����ľ�������Ϊ1/2ͼ��б�ǳ��ȣ�̽��������=ͼ��б�ǳ��ȼ�����ã��е�Ĺ켣�İ뾶=1/2̽�������ȣ�
		*/
		Point2d* Trm;
		Point2d* Rcv;
		Trm = new Point2d[sensor_length];
		Rcv = new Point2d[sensor_length];
		Set_Transmittors_and_Sensors_parallel_line(current_degree, sensor_length, Trm, Rcv);//�������ߺ�̽������λ��

		/*    ��ʽ��ʼͶӰ����    */
		double* projection_1d;
		projection_1d = new double[sensor_length]; //��¼�ض��Ƕ��µĸ���̽�������յ���һϵ��ͶӰ

		//iСѭ�����ض�ĳ���Ƕ�theta_degree�µĵĸ�����������or̽������
		for (int i = 0; i < sensor_length; ++i)
			//for (int i = 55; i < 56; ++i)
		{
			//����Trm[i]��Rcv[i]�����꣬�����ӦRcv[i]���յ���ͶӰֵprojection_1d[i]
			//��ʼ����i��ͶӰֵ���Լ���i���������ͽ�������Ӧ��alpha������������Сֵ
			//alpha�ĺ���ο����ף�Siddon, Robert L. "Fast calculation of the exact radiological path for a three�\dimensional CT array." Medical physics 12.2 (1985): 252-255.
			projection_1d[i] = 0.0;

			bool flag_intersection = 1;

			//����min_alpha ��max_alpha
			double min_alpha = 0.0, max_alpha = 1.0;
			Calculate_range_of_alpha(i, Trm, Rcv, image, flag_intersection, min_alpha, max_alpha);
			if (flag_intersection == 0) { continue; }//û�н��㣬���߽���һ����Ľ��㣨�������߶Σ���ͶӰΪ0

			//�������е������ʾ�Ѿ������min_alpha��max_alpha������ȷ�����������ͼ����ڽ���
			//���������ͼ���ཻ��m��n�ķ�Χ(m,nΪͼ���λ��������
			int intersection_m_min = 0, intersection_n_min = 0;
			int intersection_m_max = image.rows, intersection_n_max = image.cols;
			Calculate_range_of_index_to_image(i, Trm, Rcv, image, min_alpha, max_alpha,
				intersection_n_min, intersection_n_max, intersection_m_min, intersection_m_max);

			//������������m��ʵ�ʴ�����Χ������n��ʵ�ʴ�����Χ
			//������������е�alpha_x�����е�alpha_y
			//���б��intersection_n_min�����б��intersection_n_max�����н����alphaֵ
			//�ȼ�����ֱ��x=intersection_n_min-image.cols/2������ֱ��x=intersection_n_max-image.cols/2�����н����alphaֵ
			int sizeof_alpha;
			double* alpha;
			if ((current_degree % 90) != 0)
			{
				//���alpha_x
				double* alpha_x;
				int sizeof_alpha_x = intersection_n_max - intersection_n_min + 1;
				alpha_x = new double[sizeof_alpha_x];
				Calculate_alpha_x(i, Trm, Rcv, image, intersection_n_min, sizeof_alpha_x, alpha_x);

				//���alpha_y
				double* alpha_y;
				int sizeof_alpha_y = intersection_m_max - intersection_m_min + 1;
				//cout << "sizeof_alpha_y = " << sizeof_alpha_y << endl;
				alpha_y = new double[sizeof_alpha_y];
				Calculate_alpha_y(i, Trm, Rcv, image, intersection_m_min, sizeof_alpha_y, alpha_y);

				//��������Ҫ��alpha_x��alpha_y�ںϳ�һ�����飬�����������У��õ�alpha�������
				sizeof_alpha = sizeof_alpha_x + sizeof_alpha_y;
				alpha = new double[sizeof_alpha];
				Calculate_alpha(sizeof_alpha_x, sizeof_alpha, alpha_x, alpha_y, alpha);
			}
			else
			{
				if ((current_degree % 180) == 0)//ƽ����x�ᣬû��alpha_y
				{
					//���alpha_x
					double* alpha_x;
					int sizeof_alpha_x = intersection_n_max - intersection_n_min + 1;
					alpha_x = new double[sizeof_alpha_x];
					Calculate_alpha_x(i, Trm, Rcv, image, intersection_n_min, sizeof_alpha_x, alpha_x);

					//���alpha
					sizeof_alpha = sizeof_alpha_x;
					alpha = alpha_x;
					sort(alpha, alpha + sizeof_alpha);
				}
				else//ƽ����y�ᣬû��alpha_x
				{
					//���alpha_y
					double* alpha_y;
					int sizeof_alpha_y = intersection_m_max - intersection_m_min + 1;
					alpha_y = new double[sizeof_alpha_y];
					Calculate_alpha_y(i, Trm, Rcv, image, intersection_m_min, sizeof_alpha_y, alpha_y);

					sizeof_alpha = sizeof_alpha_y;
					alpha = alpha_y;
					sort(alpha, alpha + sizeof_alpha);
				}
			}


			//��������Ҫȡ�ڽ���alphaֵ�����������У��������ȣ�alphaֵ�� �� �����ܳ��ȣ����Ӧ�Ҷ�ֵ���
			//��Ӧ����������ͨ��mid_alpha = 0.5*[alpha(k)+alpha(k-1)]����
			projection_1d[i] = Calculate_projection_oneValue(i, Trm, Rcv, image, sizeof_alpha, alpha);
		}
		//��һάͶӰд��radon����Ķ�Ӧ��
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
	//��֪��ΪA1��x1��y1����A2��x2��y2���������������ɵ�ֱ��Ϊ��(y2-y1)x-(x2-x1)y+(x2-x1)y1-(y2-y1)x1 = 0
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

	//����ȡ�������⣬���м�������Խ�磬����+1֮��Խ�磬ֱ�������һ�����ݼ���
	if (k + 1 > sensor_length - 1)
	{
		cout << "Error!distance out of range(>max)." << endl;
		grayscale = radon_Img.at<double>(sensor_length - 1, theta_degree);
		return grayscale;
	}
	/*cout << "k = " << k << endl;
	cout << "distance - double(k) = " << distance - double(k) << endl;*/

	//һ�β�ֵ
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

	//������������
	for (int theta_degree = 0; theta_degree < num_degrees; ++theta_degree)
	//for (int theta_degree = 3; theta_degree < 4; ++theta_degree)
	{
		//float theta = float(theta_degree) * pi / 180;

		/*	�������ߺ�̽����(����ο�ͶӰ���̣�	*/
		Point2d* Trm;
		Point2d* Rcv;
		Trm = new Point2d[sensor_length];
		Rcv = new Point2d[sensor_length];
		Set_Transmittors_and_Sensors_parallel_line(theta_degree, sensor_length, Trm, Rcv);//�������ߺ�̽������λ��

		for (int i = 0; i < img_rows; ++i)
		//for (int i = 100; i < 101; ++i)
		{
			for (int j = 0; j < img_cols; ++j)
			//for (int j = 100; j < 101; ++j)
			{
				Point2d center_point;
				center_point.x = (double(j) - double(img_cols) / 2);
				center_point.y = (-double(i) + double(img_rows) / 2);

				//��� ���ص����ġ��㡿center_point �� Trm[0]��Rcv[0]���ɵġ�ֱ�ߡ� �ľ��� 
				//����ľ��������ӦRcv����ţ�Rcv����̽����֮������о���Ϊ1����λ��
				double distance = Cal_distance_point_to_line(center_point, Trm[0], Rcv[0]);
				//cout << "distance = " << distance << endl;
				//���ݵ㵽ֱ�ߵľ��룬���������Ҫ��ͶӰ��Ӧ��Trm[k]�����k���ݴ˲�ֵ���ض�Ӧ��ͶӰ�ĻҶ�ֵ
				double grayscale = Get_grayscale_iradon(distance, radon_Img, sensor_length, theta_degree);
				//cout << "grayscale = " << grayscale << endl;
				output_Img.at<double>(i, j) += grayscale;
			}

		}
	}
	
	//���Է�ͶӰ�Ĵ�����ʹ����ֵ��������Χ��
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

	//������������
	int theta_degree = current_degree;
	/*	�������ߺ�̽����(����ο�ͶӰ���̣�	*/
	Point2d* Trm;
	Point2d* Rcv;
	Trm = new Point2d[sensor_length];
	Rcv = new Point2d[sensor_length];
	Set_Transmittors_and_Sensors_parallel_line(theta_degree, sensor_length, Trm, Rcv);//�������ߺ�̽������λ��

	for (int i = 0; i < img_rows; ++i)
		//for (int i = 100; i < 101; ++i)
	{
		for (int j = 0; j < img_cols; ++j)
			//for (int j = 100; j < 101; ++j)
		{
			Point2d center_point;
			center_point.x = (double(j) - double(img_cols) / 2);
			center_point.y = (-double(i) + double(img_rows) / 2);

			//��� ���ص����ġ��㡿center_point �� Trm[0]��Rcv[0]���ɵġ�ֱ�ߡ� �ľ��� 
			//����ľ��������ӦRcv����ţ�Rcv����̽����֮������о���Ϊ1����λ��
			double distance = Cal_distance_point_to_line(center_point, Trm[0], Rcv[0]);
			//cout << "distance = " << distance << endl;
			//���ݵ㵽ֱ�ߵľ��룬���������Ҫ��ͶӰ��Ӧ��Trm[k]�����k���ݴ˲�ֵ���ض�Ӧ��ͶӰ�ĻҶ�ֵ
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
//	//������������
//	//for (int theta_degree = 0; theta_degree < num_degrees; ++theta_degree)
//	for (int theta_degree = current_degree; theta_degree < current_degree + 1; ++theta_degree)
//		//for (int theta_degree = 3; theta_degree < 4; ++theta_degree)
//	{
//		//float theta = float(theta_degree) * pi / 180;
//
//		/*	�������ߺ�̽����(����ο�ͶӰ���̣�	*/
//		Point2d* Trm;
//		Point2d* Rcv;
//		Trm = new Point2d[sensor_length];
//		Rcv = new Point2d[sensor_length];
//		Set_Transmittors_and_Sensors_parallel_line(theta_degree, sensor_length, Trm, Rcv);//�������ߺ�̽������λ��
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
//				//��� ���ص����ġ��㡿center_point �� Trm[0]��Rcv[0]���ɵġ�ֱ�ߡ� �ľ��� 
//				//����ľ��������ӦRcv����ţ�Rcv����̽����֮������о���Ϊ1����λ��
//				double distance = Cal_distance_point_to_line(center_point, Trm[0], Rcv[0]);
//				//cout << "distance = " << distance << endl;
//				//���ݵ㵽ֱ�ߵľ��룬���������Ҫ��ͶӰ��Ӧ��Trm[k]�����k���ݴ˲�ֵ���ض�Ӧ��ͶӰ�ĻҶ�ֵ
//				double grayscale = Get_grayscale_iradon(distance, radon_Img, sensor_length, theta_degree);
//				//cout << "grayscale = " << grayscale << endl;
//				output_Img.at<double>(i, j) += grayscale;
//			}
//
//		}
//	}
//}
