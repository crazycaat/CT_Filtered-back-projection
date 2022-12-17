#include <iostream>
#include<opencv2/opencv.hpp>

#include "radon_fan_beam.h"
#include "basic_functions.h"
#include "basic_radon_functions.h"

//#include <opencv2/highgui/highgui_c.h> 

using namespace std;
using namespace cv;




//�Լ�д��������ͶӰ����
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
//	double R_circle = max(image.rows, image.cols) * sqrt(2) / 2 + 2;  //+2��֤����ɢ�����£�Բ��Ȼ����ͼ����Χ
//	
//	for (int current_rotate_degree = 0; current_rotate_degree < num_rotate_degree; ++current_rotate_degree)
//	//for (int current_rotate_degree = 180; current_rotate_degree < 181; ++current_rotate_degree)
//	{
//		/*
//		���������ͽ�����λ�õ�ȷ����
//		current_rotate_degree������������x��������ļнǣ�
//		Rcv��˳ʱ�뷽�����У�ԭ�㵽���Ľ����������ߣ���Ϊ������������x��������н�Ϊcurrent_rotate_degree��
//
//		Trm�����и���Ԫ����ȣ�������ֻ��һ������������Ϊ�˺�ƽ�����ĺ���������ݣ�ѡ��������洢��
//		Rcv����num_detector_degree��Ԫ�أ�����˳ʱ����ת1�㣬Rcv��������һ��
//		���num_detector_degreeҲ�� Trm��Rcv[0]���� �� Trm��Rcv[num_detector_degree-1]���� ֮��ļнǡ�
//
//		Trm �� Rcv����Բ��Ϊԭ�㣬�뾶ΪR_circle��Բ��
//		
//		*/
//		Point2d* Trm;
//		Point2d* Rcv;
//		Trm = new Point2d[num_detector_degree];
//		Rcv = new Point2d[num_detector_degree];
//		Set_Transmittors_and_Sensors_fan_beam(current_rotate_degree, num_detector_degree, R_circle, Trm, Rcv);//�������ߺ�̽������λ��
//
//		/*    ��ʽ��ʼͶӰ����    */
//		double* projection_1d;
//		projection_1d = new double[num_detector_degree]; //��¼�������ض��Ƕȣ�λ�ã��µĸ���̽�������յ���һϵ��ͶӰ
//
//		//iСѭ�����ض�ĳ���Ƕ�current_rotate_degree���������̶����µĵĸ���̽����
//		for (int i = 0; i < num_detector_degree; ++i)
//		//for (int i = 29; i < 30; ++i)
//		{
//			//����Trm[i]��Rcv[i]�����꣬�����ӦRcv[i]���յ���ͶӰֵprojection_1d[i]
//			//��ʼ����i��ͶӰֵ���Լ���i���������ͽ�������Ӧ��alpha������������Сֵ
//			//alpha�ĺ���ο����ף�Siddon, Robert L. "Fast calculation of the exact radiological path for a three�\dimensional CT array." Medical physics 12.2 (1985): 252-255.
//			projection_1d[i] = 0.0;
//
//			bool flag_intersection = 1;
//
//			//����min_alpha ��max_alpha
//			double min_alpha = 0.0, max_alpha = 1.0;
//			Calculate_range_of_alpha(i, Trm, Rcv, image, flag_intersection, min_alpha, max_alpha);
//			if (flag_intersection == 0) { continue; }//û�н��㣬���߽���һ����Ľ��㣨�������߶Σ���ͶӰΪ0
//
//			//�������е������ʾ�Ѿ������min_alpha��max_alpha������ȷ�����������ͼ����ڽ���
//			//���������ͼ���ཻ��m��n�ķ�Χ(m,nΪͼ���λ��������
//			int intersection_m_min = 0, intersection_n_min = 0;
//			int intersection_m_max = image.rows, intersection_n_max = image.cols;
//			Calculate_range_of_index_to_image(i, Trm, Rcv, image, min_alpha, max_alpha,
//				intersection_n_min, intersection_n_max, intersection_m_min, intersection_m_max);
//
//			//������������m��ʵ�ʴ�����Χ������n��ʵ�ʴ�����Χ
//			//������������е�alpha_x�����е�alpha_y
//			//���б��intersection_n_min�����б��intersection_n_max�����н����alphaֵ
//			//�ȼ�����ֱ��x=intersection_n_min-image.cols/2������ֱ��x=intersection_n_max-image.cols/2�����н����alphaֵ
//			int sizeof_alpha;
//			double* alpha;
//			int theta_degree_ray = current_rotate_degree + num_detector_degree - 2 * i - 1;
//			if ((theta_degree_ray % 90) != 0)
//			{
//				//���alpha_x
//				double* alpha_x;
//				int sizeof_alpha_x = intersection_n_max - intersection_n_min + 1;
//				alpha_x = new double[sizeof_alpha_x];
//				Calculate_alpha_x(i, Trm, Rcv, image, intersection_n_min, sizeof_alpha_x, alpha_x);
//
//				//���alpha_y
//				double* alpha_y;
//				int sizeof_alpha_y = intersection_m_max - intersection_m_min + 1;
//				//cout << "sizeof_alpha_y = " << sizeof_alpha_y << endl;
//				alpha_y = new double[sizeof_alpha_y];
//				Calculate_alpha_y(i, Trm, Rcv, image, intersection_m_min, sizeof_alpha_y, alpha_y);
//
//				//��������Ҫ��alpha_x��alpha_y�ںϳ�һ�����飬�����������У��õ�alpha�������
//				sizeof_alpha = sizeof_alpha_x + sizeof_alpha_y;
//				alpha = new double[sizeof_alpha];
//				Calculate_alpha(sizeof_alpha_x, sizeof_alpha, alpha_x, alpha_y, alpha);
//			}
//			else
//			{
//				if ((theta_degree_ray % 180) == 0)//ƽ����x�ᣬû��alpha_y
//				{
//					//���alpha_x
//					double* alpha_x;
//					int sizeof_alpha_x = intersection_n_max - intersection_n_min + 1;
//					alpha_x = new double[sizeof_alpha_x];
//					Calculate_alpha_x(i, Trm, Rcv, image, intersection_n_min, sizeof_alpha_x, alpha_x);
//
//					//���alpha
//					sizeof_alpha = sizeof_alpha_x;
//					alpha = alpha_x;
//					sort(alpha, alpha + sizeof_alpha);
//				}
//				else//ƽ����y�ᣬû��alpha_x
//				{
//					//���alpha_y
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
//			//��������Ҫȡ�ڽ���alphaֵ�����������У��������ȣ�alphaֵ�� �� �����ܳ��ȣ����Ӧ�Ҷ�ֵ���
//			//��Ӧ����������ͨ��mid_alpha = 0.5*[alpha(k)+alpha(k-1)]����
//			projection_1d[i] = Calculate_projection_oneValue(i, Trm, Rcv, image, sizeof_alpha, alpha);
//		}
//		//��һάͶӰд��radon����Ķ�Ӧ��
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
	double R_circle = max(image.rows, image.cols) * sqrt(2) / 2 + 2;  //+2��֤����ɢ�����£�Բ��Ȼ����ͼ����Χ

	for (int current_rotate_degree = 0; current_rotate_degree < num_rotate_degree; ++current_rotate_degree)
		//for (int current_rotate_degree = 180; current_rotate_degree < 181; ++current_rotate_degree)
	{
		/*
		���������ͽ�����λ�õ�ȷ����
		current_rotate_degree������������x��������ļнǣ�
		Rcv��˳ʱ�뷽�����У�ԭ�㵽���Ľ����������ߣ���Ϊ������������x��������н�Ϊcurrent_rotate_degree��

		Trm�����и���Ԫ����ȣ�������ֻ��һ������������Ϊ�˺�ƽ�����ĺ���������ݣ�ѡ��������洢��
		Rcv����num_detector_degree��Ԫ�أ�����˳ʱ����ת1�㣬Rcv��������һ��
		���num_detector_degreeҲ�� Trm��Rcv[0]���� �� Trm��Rcv[num_detector_degree-1]���� ֮��ļнǡ�

		Trm �� Rcv����Բ��Ϊԭ�㣬�뾶ΪR_circle��Բ��

		*/
		Point2d* Trm;
		Point2d* Rcv;
		Trm = new Point2d[num_detector_degree];
		Rcv = new Point2d[num_detector_degree];
		Set_Transmittors_and_Sensors_fan_beam(current_rotate_degree, num_detector_degree, R_circle, Trm, Rcv);//�������ߺ�̽������λ��

		/*    ��ʽ��ʼͶӰ����    */
		double* projection_1d;
		projection_1d = new double[num_detector_degree]; //��¼�������ض��Ƕȣ�λ�ã��µĸ���̽�������յ���һϵ��ͶӰ

		//iСѭ�����ض�ĳ���Ƕ�current_rotate_degree���������̶����µĵĸ���̽����
		for (int i = 0; i < num_detector_degree; ++i)
			//for (int i = 29; i < 30; ++i)
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
			int theta_degree_ray = current_rotate_degree + num_detector_degree - 2 * i - 1;
			if ((theta_degree_ray % 90) != 0)
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
				if ((theta_degree_ray % 180) == 0)//ƽ����x�ᣬû��alpha_y
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
		Project_copyto_2Darray(current_rotate_degree, projection_1d, radon_Img);

	}

	return radon_Img;
}


static void Set_Transmittors_and_Sensors_fan_beam(const int current_rotate_degree, const int num_detector_degree,
	const double R_circle, Point2d* Trm, Point2d* Rcv)
{
	//�������ݾ�Ϊ�Ƕ�ֵ����Ҫת��Ϊ�����ƣ���
	double theta_rad = double(current_rotate_degree) * pi / 180;

	//���÷�����λ�ã��ض��Ƕ��£���������������λ��Ϊһ���㣩
	for (int i = 0; i < num_detector_degree; ++i)
	{
		Trm[i].x = R_circle * cos(theta_rad + pi);
		Trm[i].y = R_circle * sin(theta_rad + pi);
	}
	//cout << "Trm = " << Trm[0] << endl;

	//���ý�����λ��
	for (int i = 0; i < num_detector_degree; ++i)
	{
		double tmp = double(current_rotate_degree) + double(num_detector_degree) - 2 * double(i) - 1;
		tmp = tmp * pi / 180;
		Rcv[i].x = R_circle * cos(tmp);
		Rcv[i].y = R_circle * sin(tmp);
		//cout << "Rcv[" << i << "] = " << Rcv[i] << endl;
	}
}


