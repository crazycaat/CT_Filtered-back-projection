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
	//equalizeHist(radon_Img_display, radon_Img_display); //直方图均衡化，输入为单通道8位灰度图像

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

	//输出灰度值范围处理之后的图像像素值
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
			//outputImg.at<double>(i, j) = 6 * log((double)(outputImg.at<double>(i, j)) + 1);  //对数变换 s=6*log(r+1)
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

	//将滤波后的投影数据写为txt
	if (WriteData("./txt/" + image_name + "_" + type_name + ".txt", image) == 0)
	{
		cout << type_name << " write mat to txt SUCCESS!" << endl;
	}
	//由于实际直接反投影的图片的灰度范围集中在10000左右，导致图像直接显示无法体现出灰度差异，类似二值图像效果。
	//因此需要对灰度范围进行调整显示
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
 * 功能 : 将 cv::Mat 数据写入到 .txt 文件
 *----------------------------
 * 函数 : WriteData
 * 访问 : public
 * 返回 : -1：打开文件失败；0：写入数据成功；1：矩阵为空
 *
 * 参数 : fileName [in] 文件名
 * 参数 : matData [in] 矩阵数据
 */
int WriteData(string fileName, cv::Mat& matData)
{
	int retVal = 0;

	// 打开文件 
	ofstream outFile(fileName.c_str(), ios_base::out); //按新建或覆盖方式写入 
	if (!outFile.is_open())
	{
		cout << "打开文件失败" << endl;
		retVal = -1;
		return (retVal);
	}

	// 检查矩阵是否为空 
	if (matData.empty())
	{
		cout << "矩阵为空" << endl;
		retVal = 1;
		return (retVal);
	}

	// 写入数据 
	for (int r = 0; r < matData.rows; r++)
	{
		for (int c = 0; c < matData.cols; c++)
		{
			//if ((c+1) % 256 == 0 && c != 0)	outFile << endl;//达到了txt文本格式每行的最大字符数
			double data = matData.at<double>(r, c); //读取数据，at<type> - type 是矩阵元素的具体数据格式 
			outFile << data << "\t"; //每列数据用 tab 隔开 
		}
		outFile << endl; //换行 
	}

	return (retVal);
}


/*----------------------------
 * 功能 : 从 .txt 文件中读入数据，保存到 cv::Mat 矩阵
 * - 默认按 float 格式读入数据，
 * - 如果没有指定矩阵的行、列和通道数，则输出的矩阵是单通道、N 行 1 列的
 *----------------------------
 * 函数 : LoadData
 * 访问 : public
 * 返回 : -1：打开文件失败；0：按设定的矩阵参数读取数据成功；1：按默认的矩阵参数读取数据
 *
 * 参数 : fileName [in] 文件名
 * 参数 : matData [out] 矩阵数据
 * 参数 : matRows [in] 矩阵行数，默认为 0
 * 参数 : matCols [in] 矩阵列数，默认为 0
 * 参数 : matChns [in] 矩阵通道数，默认为 0
 */
int LoadData(string fileName, cv::Mat& matData, int matRows, int matCols, int matChns)
{
	int retVal = 0;

	// 打开文件 
	ifstream inFile(fileName.c_str(), ios_base::in);
	if (!inFile.is_open())
	{
		cout << "读取文件失败" << endl;
		retVal = -1;
		return (retVal);
	}

	// 载入数据 
	istream_iterator<float> begin(inFile); //按 float 格式取文件数据流的起始指针 
	istream_iterator<float> end; //取文件流的终止位置 
	vector<float> inData(begin, end); //将文件数据保存至 std::vector 中 
	cv::Mat tmpMat = cv::Mat(inData); //将数据由 std::vector 转换为 cv::Mat 

	// 输出到命令行窗口 
	//copy(vec.begin(),vec.end(),ostream_iterator<double>(cout,"\t")); 

	// 检查设定的矩阵尺寸和通道数 
	size_t dataLength = inData.size();
	//1.通道数 
	if (matChns == 0)
	{
		matChns = 1;
	}
	//2.行列数 
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
	//3.数据总长度 
	if (dataLength != (matRows * matCols * matChns))
	{
		cout << "读入的数据长度 不满足 设定的矩阵尺寸与通道数要求，将按默认方式输出矩阵！" << endl;
		retVal = 1;
		matChns = 1;
		matRows = dataLength;
	}

	// 将文件数据保存至输出矩阵 
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
	if (image1.size == image2.size)//原图尺寸和重建图尺寸完全相同
	{
		//计算原始图像image1的像素平均值
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
	else//重建图尺寸和原图尺寸不完全匹配
	{
		cout << "Caution!size not paired!" << endl;
		cout << "The size of original image = " << image1.size << endl;
		cout << "The size of iradon image = " << image2.size << endl;

		//计算原始图像image1的像素平均值
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

		//尺寸不匹配，把小图放在大图中央比较（或仅比较可重叠部分）
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
	if (image1.size == image2.size)//原图尺寸和重建图尺寸完全相同
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
	else//重建图尺寸和原图尺寸不完全匹配
	{
		cout << "Caution!size not paired!" << endl;
		cout << "The size of original image = " << image1.size << endl;
		cout << "The size of iradon image = " << image2.size << endl;

		double sum_numerator = 0;
		double sum_dominator = 0;

		//尺寸不匹配，把小图放在大图中央比较（或仅比较可重叠部分）
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