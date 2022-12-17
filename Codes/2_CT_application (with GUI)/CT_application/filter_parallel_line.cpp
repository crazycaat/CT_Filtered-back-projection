#include <iostream>
#include<opencv2/opencv.hpp>
#include <math.h>

#include "filter_parallel_line.h"
#include "basic_functions.h"

#define  pi  3.14159265358979323846  //!< const value,pi
using namespace cv;

Mat Filter_radon_parallel_line(Mat inputImg,int filter_flag)
{
	if (filter_flag == non_filter)
	{
		//  do not need to apply any filters
		//  direct radon back
		return inputImg;
	}


	//Mat filtered_radon_Img = inputImg.clone();


	/*inputImg = Mat::zeros(7, 4, COLOR_BGR2GRAY);
	inputImg.at<double>(0, 0) = 1;
	inputImg.at<double>(0, 1) = 2;
	inputImg.at<double>(0, 2) = 3;
	inputImg.at<double>(0, 3) = 4;
	inputImg.at<double>(0, 0) = 1;
	inputImg.at<double>(1, 0) = 2;
	inputImg.at<double>(2, 0) = 3;
	inputImg.at<double>(3, 0) = 4;

	inputImg.at<double>(0, 1) = 5;
	inputImg.at<double>(1, 1) = 6;
	inputImg.at<double>(2, 1) = 7;
	inputImg.at<double>(3, 1) = 8;*/

	//inputImg.at<double>(0, 2) = 4;
	//inputImg.at<double>(1, 2) = 3;
	//inputImg.at<double>(2, 2) = 2;
	//inputImg.at<double>(3, 2) = 2;

	//inputImg.at<double>(0, 3) = 2;
	//inputImg.at<double>(1, 3) = 4;
	//inputImg.at<double>(2, 3) = 1;
	//inputImg.at<double>(3, 3) = 5;



	cout << "start filtering..." << endl;
	Mat filtered_radon_Img = inputImg.clone();
	//ѭ�������зֱ����dft�����ҽ�ʵ�����鲿�ֱ����Output_real_dft��Output_imaginary_dft�Ķ�Ӧ��
	for (int j = 0; j < inputImg.cols; ++j)
	//for (int j = 0; j < 1; ++j)
	{
		cout << "degree = " << j << endl;
		//1.��ȡ��------------------------------------------------------------------------------------
		Mat I = inputImg(cv::Range::all(), cv::Range(j, j + 1));

		//2.������һ�е�dft----------���ο�cv��dft�����Ĺٷ����Դ��롿-----------------------------------------------------------------
		Mat padded;                 //��0�������ͼ�����
		int m = getOptimalDFTSize(I.rows);
		int n = getOptimalDFTSize(I.cols);

		//�������ͼ��I���������Ϊpadded���Ϸ����󷽲�����䴦��
		copyMakeBorder(I, padded, 0, m - I.rows, 0, n - I.cols, BORDER_CONSTANT, Scalar::all(0));

		Mat planes[] = { Mat_<float>(padded), Mat::zeros(padded.size(),CV_32F) };
		Mat complexI;
		merge(planes, 2, complexI);     //��planes�ںϺϲ���һ����ͨ������complexI

		dft(complexI, complexI);        //���и���Ҷ�任

		//�����ֵ��ת���������߶�(logarithmic scale)
		//=> log(1 + sqrt(Re(DFT(I))^2 + Im(DFT(I))^2))
		split(complexI, planes);        //planes[0] = Re(DFT(I),planes[1] = Im(DFT(I))
										//��planes[0]Ϊʵ��,planes[1]Ϊ�鲿
		Mat magI;
		magnitude(planes[0], planes[1], magI);     //planes[0] = magnitude

		/*cout << "radon_1d = " << endl;
		cout << I << endl;
		cout << "dft result = " << endl;
		cout << complexI << endl;*/
		//cout << "real part = " << endl;
		//cout << planes[0] << endl;
		//cout << "imaginary part = " << endl;
		//cout << planes[1] << endl;
		//cout << "magI = " << endl;
		//cout << magI << endl;

		//;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		//3.����һ�н����˲���R-L�˲����� or S-L�˲�������---------------------------------------------------------
		//R-L�˲�����: H_{RL} = |w|  (when |w|<B)
		if (filter_flag == RL_filter)
		{
			int DFT_length = complexI.rows;
			//�˲���(|w|����ԭ����������м�)��1��2��3��2��1�����жԳ��ԣ����ҿ��Թ�һ������0.33,0.66,1,0.66,0.33
			for (int w = 0; w < floor(0.5 * DFT_length); ++w)
			{
				double filterFactor = double(w) / (DFT_length / 2);
				//cout << "w = " << w << ", filterFactor = " << filterFactor << endl;
				complexI(cv::Range(w, w + 1), cv::Range::all()) *= filterFactor;//��w��
				complexI(cv::Range(DFT_length - w - 1, DFT_length - w), cv::Range::all()) *= filterFactor;//��DFT_length-w�������жԳ��ԣ�filterFactor���
			}
		}

		//S-L�˲�����: H_{SL} = |w|*sinc(w/2B)  (when |w|<B)
		if (filter_flag == SL_filter)
		{
			
			int DFT_length = complexI.rows;
			double Band_width = 1 * DFT_length;//���ô�����СGibbs effect
			for (int w = 0; w < floor(0.5 * DFT_length); ++w)
			{
				double filterFactor = double(w);
				double tmp = pi * ((w - floor(0.5 * DFT_length)) / (2 * Band_width));
				if (tmp != 0)
				{
					filterFactor *= sin(tmp) / tmp;//�ٳ���sinc����
				}
				else
				{
					filterFactor *= 1;
				}

				complexI(cv::Range(w, w + 1), cv::Range::all()) *= filterFactor;//��w��
				complexI(cv::Range(DFT_length - w - 1, DFT_length - w), cv::Range::all()) *= filterFactor;//��DFT_length-w-1�������жԳ��ԣ�filterFactor���
			}
		}
	
		

		//4.����һ��idft���������˲����Mat
		Mat idft_1d_result;
		//merge(planes, 2, complexI);     //��planes�ںϺϲ���һ����ͨ������complexI
		idft(complexI, idft_1d_result, DFT_SCALE);
		split(idft_1d_result, planes);
		Mat tmp = planes[0];//ȡʵ����Ϊidft�Ľ��
		//magnitude(planes[0], planes[1], planes[0]);     //planes[0] = magnitude
		//tmp = planes[0];
		for (int i = 0; i < inputImg.rows; ++i)
		{
			filtered_radon_Img.at<double>(i, j) = tmp.at<float>(i, 0);
		}
	}

	/*cout << "filtered_radon_Img = " << endl;
	cout << filtered_radon_Img << endl;*/

	return filtered_radon_Img;
	
}
