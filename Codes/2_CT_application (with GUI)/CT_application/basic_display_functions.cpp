#include <iostream>
#include<opencv2/opencv.hpp>

#include <fstream> 
#include <iterator> 
#include <vector> 

#include "basic_display_functions.h"
#include "basic_functions.h"

using namespace std;
using namespace cv;

Mat Convert_to_show_image(Mat Img)
{
	Mat Img_display = Img.clone();
	//radon_Img.convertTo(radon_Img_display, COLOR_BGR2BGRA);
	//equalizeHist(radon_Img_display, radon_Img_display); //ֱ��ͼ���⻯������Ϊ��ͨ��8λ�Ҷ�ͼ��

	double max_grayscale_value = 0;
	for (int m = 0; m < Img_display.rows; ++m)
	{
		for (int n = 0; n < Img_display.cols; ++n)
		{
			if (Img_display.at<double>(m, n) > max_grayscale_value)
				max_grayscale_value = Img_display.at<double>(m, n);
		}
	}
	cout << endl;
	cout << "max_grayscale_value = " << max_grayscale_value << endl;
	cout << "mean_grayscale_value = " << sumMat(Img) / (Img.rows * Img.cols) << endl;

	//DrawHist(radon_Img);

	for (int n = 0; n < Img_display.cols; ++n)
	{
		for (int m = 0; m < Img_display.rows; ++m)
		{
			//radon_Img_display.at<double>(m, n) = ((radon_Img_display.at<double>(m, n) - 7000) / (max_grayscale_value - 7000)) * 255;
			Img_display.at<double>(m, n) = ((Img_display.at<double>(m, n) - max_grayscale_value / 3) / (max_grayscale_value - max_grayscale_value / 3)) * 255;
			if (Img_display.at<double>(m, n) < 0) Img_display.at<double>(m, n) = 0;
		}
	}
	//radon_Img.convertTo(radon_Img_display, COLOR_BGR2BGRA);
	//normalize(radon_Img_display, radon_Img_display, 0, 255, NORM_MINMAX);

	//����Ҷ�ֵ��Χ����֮���ͼ������ֵ
	//for (int n = 0; n < num_degrees; ++n)
	/*for (int n = 136; n < 145; ++n)
	{
		cout << "degree = " << n << endl;
		for (int m = 0; m < sensor_length; ++m)
		{
			cout << "sensor" << m << "(grayscale normalized) = " << int(radon_Img_display.at<double>(m, n)) << endl;
		}
	}*/
	return Img_display;

}

Mat Convert_to_show_normalize(Mat Img)
{
	Mat Img_display = Img.clone();
	double max_grayscale_value = 0;
	double min_grayscale_value = 1000000000;
	for (int m = 0; m < Img_display.rows; ++m)
	{
		for (int n = 0; n < Img_display.cols; ++n)
		{
			if (Img_display.at<double>(m, n) > max_grayscale_value)
				max_grayscale_value = Img_display.at<double>(m, n);
			if (Img_display.at<double>(m, n) < min_grayscale_value)
				min_grayscale_value = Img_display.at<double>(m, n);
		}
	}
	cout << endl;
	cout << "max_grayscale_value = " << max_grayscale_value << endl;
	cout << "min_grayscale_value = " << min_grayscale_value << endl;
	cout << "mean_grayscale_value = " << sumMat(Img) / (Img.rows * Img.cols) << endl;

	for (int n = 0; n < Img_display.cols; ++n)
	{
		for (int m = 0; m < Img_display.rows; ++m)
		{
			Img_display.at<double>(m, n) = ((Img_display.at<double>(m, n) - min_grayscale_value) / (max_grayscale_value - min_grayscale_value)) * 255;
			if (Img_display.at<double>(m, n) < 0) Img_display.at<double>(m, n) = 0;
		}
	}
	return Img_display;

}


Mat Convert_to_show_equalizeHist(Mat inputImg)
{
	Mat outputImg;
	inputImg.convertTo(inputImg, COLOR_BGR2BGRA);
	equalizeHist(inputImg, outputImg);
	return outputImg;
}

Mat Convert_to_show_LogTrans(Mat inputImg, double gamma)
{
	Mat outputImg = inputImg.clone();
	for (int i = 0; i < outputImg.rows; i++)
	{
		for (int j = 0; j < outputImg.cols; j++)
		{
			//outputImg.at<double>(i, j) = 6 * log((double)(outputImg.at<double>(i, j)) + 1);  //�����任 s=6*log(r+1)
			outputImg.at<double>(i, j) = pow((outputImg.at<double>(i, j)), gamma);
		}
	}
	return outputImg;
}


void showImg(string windowname, Mat img)
{
	namedWindow(windowname, WINDOW_GUI_EXPANDED);
	//resizeWindow(windowname, 256, 256);
	imshow(windowname, img);
}


/** @brief write and show the image(equalized)

@param   image: Input the image (8-bit 1-channel image).
@param   image_name: the original image name.(eg.shepplogan,circle,rec,line)
@param   type_name: the type of the image(eg.filtered_radon,direct_iradon,filtered_iradon)

 */
void write_and_show_Img(Mat image, string image_name, string type_name)
{
	//showImg(type_name + " result", image);
	//imwrite("./pics/" + image_name + "_output_" + type_name + ".png", image);
	cout << "the size of " + type_name + " = (" << image.rows << ", " << image.cols << ")" << endl;

	//���˲����ͶӰ����дΪtxt
	if (WriteData("./txt/" + image_name + "_" + type_name + ".txt", image) == 0)
	{
		cout << type_name << " write mat to txt SUCCESS!" << endl;
	}
	//����ʵ��ֱ�ӷ�ͶӰ��ͼƬ�ĻҶȷ�Χ������10000���ң�����ͼ��ֱ����ʾ�޷����ֳ��ҶȲ��죬���ƶ�ֵͼ��Ч����
	//�����Ҫ�ԻҶȷ�Χ���е�����ʾ
	//Mat Img_display = Convert_to_show_image(image);
	//showImg(type_name+" equalized", Img_display);
	//imwrite("./pics/" + image_name + "_output_" + type_name + "_equalized0.png", Img_display);

	//Mat Img_display2 = Convert_to_show_normalize(image);
	////showImg(type_name+" equalized", Img_display2);
	//imwrite("./pics/" + image_name + "_output_" + type_name + "_normalized.png", Img_display2);

	//Mat Img_display3 = Convert_to_show_equalizeHist(Img_display2);
	//imwrite("./pics/" + image_name + "_output_" + type_name + "_normalized_and_equalized.png", Img_display3);

	//Mat Img_display4 = Convert_to_show_equalizeHist(image);
	//imwrite("./pics/" + image_name + "_output_" + type_name + "_equalized.png", Img_display4);

	//Mat Img_display5 = Convert_to_show_LogTrans(Img_display2);
	//Img_display5 = Convert_to_show_normalize(Img_display5);
	//imwrite("./pics/" + image_name + "_output_" + type_name + "_log.png", Img_display5);

	Mat Img_display = Convert_to_show_normalize(image);
	imwrite("./pics/" + image_name + "_output_" + type_name + ".png", Img_display);

	/*Mat Img_display_log = Convert_to_show_LogTrans(Img_display);
	Img_display_log = Convert_to_show_normalize(Img_display_log);
	imwrite("./pics/" + image_name + "_output_" + type_name + "_log.png", Img_display_log);*/


}



/*----------------------------
 * ���� : �� cv::Mat ����д�뵽 .txt �ļ�
 *----------------------------
 * ���� : WriteData
 * ���� : public
 * ���� : -1�����ļ�ʧ�ܣ�0��д�����ݳɹ���1������Ϊ��
 *
 * ���� : fileName [in] �ļ���
 * ���� : matData [in] ��������
 */
int WriteData(string fileName, cv::Mat& matData)
{
	int retVal = 0;

	// ���ļ� 
	ofstream outFile(fileName.c_str(), ios_base::out); //���½��򸲸Ƿ�ʽд�� 
	if (!outFile.is_open())
	{
		cout << "���ļ�ʧ��" << endl;
		retVal = -1;
		return (retVal);
	}

	// �������Ƿ�Ϊ�� 
	if (matData.empty())
	{
		cout << "����Ϊ��" << endl;
		retVal = 1;
		return (retVal);
	}

	// д������ 
	for (int r = 0; r < matData.rows; r++)
	{
		for (int c = 0; c < matData.cols; c++)
		{
			//if ((c+1) % 256 == 0 && c != 0)	outFile << endl;//�ﵽ��txt�ı���ʽÿ�е�����ַ���
			double data = matData.at<double>(r, c); //��ȡ���ݣ�at<type> - type �Ǿ���Ԫ�صľ������ݸ�ʽ 
			outFile << data << "\t"; //ÿ�������� tab ���� 
		}
		outFile << endl; //���� 
	}

	return (retVal);
}


/*----------------------------
 * ���� : �� .txt �ļ��ж������ݣ����浽 cv::Mat ����
 * - Ĭ�ϰ� float ��ʽ�������ݣ�
 * - ���û��ָ��������С��к�ͨ������������ľ����ǵ�ͨ����N �� 1 �е�
 *----------------------------
 * ���� : LoadData
 * ���� : public
 * ���� : -1�����ļ�ʧ�ܣ�0�����趨�ľ��������ȡ���ݳɹ���1����Ĭ�ϵľ��������ȡ����
 *
 * ���� : fileName [in] �ļ���
 * ���� : matData [out] ��������
 * ���� : matRows [in] ����������Ĭ��Ϊ 0
 * ���� : matCols [in] ����������Ĭ��Ϊ 0
 * ���� : matChns [in] ����ͨ������Ĭ��Ϊ 0
 */
int LoadData(string fileName, cv::Mat& matData, int matRows, int matCols, int matChns)
{
	int retVal = 0;

	// ���ļ� 
	ifstream inFile(fileName.c_str(), ios_base::in);
	if (!inFile.is_open())
	{
		cout << "��ȡ�ļ�ʧ��" << endl;
		retVal = -1;
		return (retVal);
	}

	// �������� 
	istream_iterator<float> begin(inFile); //�� float ��ʽȡ�ļ�����������ʼָ�� 
	istream_iterator<float> end; //ȡ�ļ�������ֹλ�� 
	vector<float> inData(begin, end); //���ļ����ݱ����� std::vector �� 
	cv::Mat tmpMat = cv::Mat(inData); //�������� std::vector ת��Ϊ cv::Mat 

	// ����������д��� 
	//copy(vec.begin(),vec.end(),ostream_iterator<double>(cout,"\t")); 

	// ����趨�ľ���ߴ��ͨ���� 
	size_t dataLength = inData.size();
	//1.ͨ���� 
	if (matChns == 0)
	{
		matChns = 1;
	}
	//2.������ 
	if (matRows != 0 && matCols == 0)
	{
		matCols = dataLength / matChns / matRows;
	}
	else if (matCols != 0 && matRows == 0)
	{
		matRows = dataLength / matChns / matCols;
	}
	else if (matCols == 0 && matRows == 0)
	{
		matRows = dataLength / matChns;
		matCols = 1;
	}
	//3.�����ܳ��� 
	if (dataLength != (matRows * matCols * matChns))
	{
		cout << "��������ݳ��� ������ �趨�ľ���ߴ���ͨ����Ҫ�󣬽���Ĭ�Ϸ�ʽ�������" << endl;
		retVal = 1;
		matChns = 1;
		matRows = dataLength;
	}

	// ���ļ����ݱ������������ 
	matData = tmpMat.reshape(matChns, matRows).clone();

	return (retVal);
}

double Calculate_error_sqrt(Mat image1, Mat image2)
{
	double error = -1;
	image1.convertTo(image1, COLOR_BGR2GRAY);
	image2.convertTo(image2, COLOR_BGR2GRAY);
	image1 = Convert_to_show_normalize(image1);
	image2 = Convert_to_show_normalize(image2);
	if (image1.size == image2.size)//ԭͼ�ߴ���ؽ�ͼ�ߴ���ȫ��ͬ
	{
		//����ԭʼͼ��image1������ƽ��ֵ
		double sum1 = 0;
		for (int i = 0; i < image1.rows; ++i)
		{
			for (int j = 0; j < image1.cols; ++j)
			{
				sum1 += image1.at<double>(i, j);

			}
		}
		double average1 = sum1 / (double(image1.rows) * double(image1.cols));

		double sum_numerator = 0;
		double sum_dominator = 0;
		for (int i = 0; i < image1.rows; ++i)
		{
			for (int j = 0; j < image1.cols; ++j)
			{
				sum_numerator += pow(image1.at<double>(i, j) - image2.at<double>(i, j), 2);
				sum_dominator += pow(image1.at<double>(i, j) - average1, 2);
			}
		}
		error = sqrt(sum_numerator / sum_dominator);
	}
	else//�ؽ�ͼ�ߴ��ԭͼ�ߴ粻��ȫƥ��
	{
		cout << "Caution!size not paired!" << endl;
		cout << "The size of original image = " << image1.size << endl;
		cout << "The size of iradon image = " << image2.size << endl;

		//����ԭʼͼ��image1������ƽ��ֵ
		double sum1 = 0;
		for (int i = 0; i < image1.rows; ++i)
		{
			for (int j = 0; j < image1.cols; ++j)
			{
				sum1 += image1.at<double>(i, j);

			}
		}
		double average1 = sum1 / (double(image1.rows) * double(image1.cols));

		double sum_numerator = 0;
		double sum_dominator = 0;

		//�ߴ粻ƥ�䣬��Сͼ���ڴ�ͼ����Ƚϣ�����ȽϿ��ص����֣�
		for (int i = 0; i < min(image1.rows, image2.rows); ++i)
		{
			for (int j = 0; j < min(image1.cols, image2.cols); ++j)
			{
				int i1, i2, j1, j2;
				if (image1.rows < image2.rows)
				{
					i1 = i;
					i2 = i + floor((image2.rows - image1.rows) / 2);
				}
				else
				{
					i1 = i + floor((image1.rows - image2.rows) / 2);
					i2 = i;
				}

				if (image1.cols < image2.cols)
				{
					j1 = j;
					j2 = j + floor((image2.cols - image1.cols) / 2);
				}
				else
				{
					j1 = j + floor((image1.cols - image2.cols) / 2);
					j2 = j;
				}
				sum_numerator += pow(image1.at<double>(i1, j1) - image2.at<double>(i2, j2), 2);
				sum_dominator += pow(image1.at<double>(i1, j1) - average1, 2);
			}
		}
		error = sqrt(sum_numerator / sum_dominator);
	}
	return error;
}


double Calculate_error_abs(Mat image1, Mat image2)
{
	double error = -1;
	image1.convertTo(image1, COLOR_BGR2GRAY);
	image2.convertTo(image2, COLOR_BGR2GRAY);
	image1 = Convert_to_show_normalize(image1);
	image2 = Convert_to_show_normalize(image2);
	if (image1.size == image2.size)//ԭͼ�ߴ���ؽ�ͼ�ߴ���ȫ��ͬ
	{
		double sum_numerator = 0;
		double sum_dominator = 0;
		for (int i = 0; i < image1.rows; ++i)
		{
			for (int j = 0; j < image1.cols; ++j)
			{
				sum_numerator += abs(image1.at<double>(i, j) - image2.at<double>(i, j));
				sum_dominator += abs(image1.at<double>(i, j));
			}
		}
		error = sum_numerator / sum_dominator;
	}
	else//�ؽ�ͼ�ߴ��ԭͼ�ߴ粻��ȫƥ��
	{
		cout << "Caution!size not paired!" << endl;
		cout << "The size of original image = " << image1.size << endl;
		cout << "The size of iradon image = " << image2.size << endl;

		double sum_numerator = 0;
		double sum_dominator = 0;

		//�ߴ粻ƥ�䣬��Сͼ���ڴ�ͼ����Ƚϣ�����ȽϿ��ص����֣�
		for (int i = 0; i < min(image1.rows, image2.rows); ++i)
		{
			for (int j = 0; j < min(image1.cols, image2.cols); ++j)
			{
				int i1, i2, j1, j2;
				if (image1.rows < image2.rows)
				{
					i1 = i;
					i2 = i + floor((image2.rows - image1.rows) / 2);
				}
				else
				{
					i1 = i + floor((image1.rows - image2.rows) / 2);
					i2 = i;
				}

				if (image1.cols < image2.cols)
				{
					j1 = j;
					j2 = j + floor((image2.cols - image1.cols) / 2);
				}
				else
				{
					j1 = j + floor((image1.cols - image2.cols) / 2);
					j2 = j;
				}
				sum_numerator += abs(image1.at<double>(i1, j1) - image2.at<double>(i2, j2));
				sum_dominator += abs(image1.at<double>(i1, j1));
			}
		}
		error = sum_numerator / sum_dominator;
	}
	return error;
}