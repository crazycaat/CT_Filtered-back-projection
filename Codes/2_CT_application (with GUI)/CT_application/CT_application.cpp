#include "CT_application.h"
#include "stdafx.h"

using namespace std;
using namespace cv;

//���캯�������ø���ѡ�򡢹�ѡ�򡢻������������ȵĳ�ʼ״̬��
CT_application::CT_application(QWidget *parent)
    : QMainWindow(parent)
{
    //��ʼ��״̬
    ui.setupUi(this);
    ui.spinBox_size_of_shepplogan->setValue(size_of_shepplogan);
    ui.horizontalSlider_size_of_shepplogan->setValue(size_of_shepplogan);

    ui.spinBox_num_degrees->setValue(radon_num_degree);
    ui.horizontalSlider_num_degrees->setValue(radon_num_degree);
    
    ui.horizontalSlider_speed->setValue(reconstruction_speed);
    

    //��ʹ���˲������Ĺ�ѡ���ʼ״̬����
    if(use_filter)
        ui.checkBox_filter->setCheckState(Qt::Checked);
    else
        ui.checkBox_filter->setCheckState(Qt::Unchecked);

    //����ʾ�ؽ����̡��Ĺ�ѡ���ʼ״̬����
    if(show_reconstruction_process)
        ui.checkBox_process->setCheckState(Qt::Checked);
    else
        ui.checkBox_process->setCheckState(Qt::Unchecked);

    //����ѡ�򡿣�ͶӰ���ͣ��ĳ�ʼ����������ع��ܵ�����/����
    if (radon_flag == radonFlags::analytical)
    {
        ui.comboBox_img->setCurrentIndex(0);//���ø�ѡ��Ĭ��Ϊ����shepploganģ��
        //�������ֻ�����ع���
        ui.horizontalSlider_num_degrees->setEnabled(false);
        ui.spinBox_num_degrees->setEnabled(false);
        //���ý�����ͶӰ��ع���
        ui.horizontalSlider_size_of_shepplogan->setEnabled(true);
        ui.spinBox_size_of_shepplogan->setEnabled(true);
    }
    else
    {
        ui.comboBox_img->setCurrentIndex(1);//���ø�ѡ��Ĭ��Ϊ�ӱ���ѡ��ͼƬ�����ֻ���
        //�������ֻ�����ع���
        ui.horizontalSlider_num_degrees->setEnabled(true);
        ui.spinBox_num_degrees->setEnabled(true);
        //���ý�����ͶӰ��ع���
        ui.horizontalSlider_size_of_shepplogan->setEnabled(false);
        ui.spinBox_size_of_shepplogan->setEnabled(false);
    }

    

    image = Create_shepplogan(size_of_shepplogan); //��ʼ��ͼ��
    original_image = image.clone();

}



void CT_application::displayMat_color(QLabel* label, cv::Mat image)
{
    Mat rgb;
    QImage img;
    QImage imgScaled;
    if (image.channels() == 3)
    {
        //cvt Mat BGR 2 QImage RGB
        cvtColor(image, rgb, COLOR_BGR2RGB);
        img = QImage((const unsigned char*)(rgb.data),
            rgb.cols, rgb.rows,
            rgb.cols * rgb.channels(),
            QImage::Format_RGB888);
    }
    else
    {
        img = QImage((const unsigned char*)(image.data),
            image.cols, image.rows,
            image.cols * image.channels(),
            QImage::Format_Indexed8);
    }
    imgScaled = img.scaled(label->size(), Qt::KeepAspectRatio);
    // ui->original_img->setPixmap(QPixmap::fromImage(img));
    label->setPixmap(QPixmap::fromImage(imgScaled));
    label->resize(label->pixmap()->size());
}


void CT_application::displayMat_gray(QLabel* label, cv::Mat image)
{
    QImage img;
    QImage imgScaled;
    
    image = Convert_to_show_normalize(image);
    image.convertTo(image, COLOR_BGR2BGRA);
   
    img = QImage((const unsigned char*)(image.data),
        image.cols, image.rows,
        image.cols * image.channels(),
        QImage::Format_Grayscale8);
    
    imgScaled = img.scaled(label->size(), Qt::KeepAspectRatio);
    // ui->original_img->setPixmap(QPixmap::fromImage(img));
    label->setPixmap(QPixmap::fromImage(imgScaled));
    label->resize(label->pixmap()->size());
}

//�������ʾͼƬ����ť
void CT_application::display_img_click()
{
    if (radon_flag == radonFlags::integral)//��ʾѡ��ӱ����ϴ�ͼƬ����ʾ
    {
        QFileDialog fileDlg;
        //QString imgFile = fileDlg.getOpenFileName(this, QString::fromLocal8Bit("����ͼƬ"), ".", QString::fromLocal8Bit("ͼƬ��*jpg *png *bmp *dfx��"));
        QString imgFile = fileDlg.getOpenFileName(this, QString::fromLocal8Bit("����ͼƬ"), ".", tr("Image Files (*.png *.jpg *.bmp)"));
        QFileInfo fileInf(imgFile);
        if (!fileInf.exists())
        {
            QMessageBox::warning(this, QString::fromLocal8Bit("����"), QString::fromLocal8Bit("���棺��ѡ��ĵ��ļ������ڣ�"), QMessageBox::Ok);
            return;
        }
        else
        {
            std::string pathStr = imgFile.toLocal8Bit();
            image = imread(pathStr.c_str());
        }
        displayMat_color(ui.original_img, image);
        cvtColor(image, image, COLOR_BGR2GRAY);
        original_image = image.clone();

        //����ͼ��ߴ�Ϊ�����Σ���0��
        if (image.rows != image.cols)
        {
            image = Resize_img_to_Square(image);
        }
    }
    else//default:Ĭ��ʹ��shepploganģ��
    {
        image = Create_shepplogan(size_of_shepplogan);
        original_image = image.clone();
        displayMat_gray(ui.original_img, image);
    }
}

//�ı����õĽ�������shepploganͼ��ߴ硾�����
void CT_application::size_of_shepplogan_changed(int value)
{
    if (first_time_flag_text == 0 && radon_flag == radonFlags::analytical)
    {
        image = Create_shepplogan(size_of_shepplogan);
        original_image = image.clone();
        displayMat_gray(ui.original_img, image);
    }
    first_time_flag_text = 0;

    size_of_shepplogan = value;
    ui.horizontalSlider_size_of_shepplogan->setValue(value); 
}

//�ı����õĽ�������shepploganͼ��ߴ硾��������
void CT_application::size_of_shepplogan_changed_slider(int value)
{
    if (first_time_flag_slider == 0 && radon_flag == radonFlags::analytical)
    {
        image = Create_shepplogan(size_of_shepplogan);
        original_image = image.clone();
        displayMat_gray(ui.original_img, image);
    }
    first_time_flag_slider = 0;

    size_of_shepplogan = value;
    ui.spinBox_size_of_shepplogan->setValue(value);
}

//�ı����ֻ��ַ���ͶӰ�ܽǶȡ������
void CT_application::num_degrees_changed(int value)
{
    radon_num_degree = value;
    ui.horizontalSlider_num_degrees->setValue(value);
}

//�ı�ı����ֻ��ַ���ͶӰ�ܽǶȡ���������
void CT_application::num_degrees_changed_slider(int value)
{
    radon_num_degree = value;
    ui.spinBox_num_degrees->setValue(value);
   
}

//�ı��ؽ��ٶȵġ���������
void CT_application::speed_changed_slider(int value)
{
    reconstruction_speed = value;
}

//�޸ġ�ʹ���˲������Ĺ�ѡ״̬������ѡ��ͬ�ڡ�ֱ�ӷ�ͶӰ������ѡ��ʹ��S-L�˲���
void CT_application::use_filter_check(bool value)
{
    use_filter = value;
}

//�޸ġ���ʾ�ؽ����̡��Ĺ�ѡ״̬
void CT_application::show_process_check(bool value)
{
    show_reconstruction_process = value;
}

//�������ʼ�ؽ�����ť
void CT_application::show_reconstruct_result_click()
{
    ui.comboBox_img->setEnabled(false);//�ؽ�����ʱ����ֹ�޸�ͼ��ͶӰ���ͣ����ø�ѡ��
    ui.Button_showimg->setEnabled(false);//�ؽ�����ʱ����ֹ�����µ�ͼƬ
    ui.pushButton->setEnabled(false);//�ؽ�����ʱ����ֹ�ظ��������ʼ�ؽ���
    
    if (radon_flag == radonFlags::analytical)
    {
        displayMat_gray(ui.original_img, image);//������δ��ʾͼ������ʾͼ��
    }
    
    cv::Mat iradon_img;
    if (show_reconstruction_process == 0)
    {
        //ͶӰ
        cv::Mat radonImg;
        if (radon_flag == radonFlags::integral)//�����ͼƬ����Ҫ���ֻ���ͶӰ
        {
            radonImg = Radon_parallel_line(image, radon_num_degree);//�Լ�д��ƽ����ͶӰ����
            displayMat_gray(ui.radon_img, radonImg);
        }
        else//default ����������shepplogan��ͶӰ
        {
            radonImg = Calc_shepplogan_radon(size_of_shepplogan);
            displayMat_gray(ui.radon_img, radonImg);
        }

        //�˲�
        if (use_filter == 1)
        {
            radonImg = Filter_radon_parallel_line(radonImg, SL_filter);
            //radonImg = Filter_radon_parallel_line(radonImg, RL_filter);
            displayMat_gray(ui.filtered_radon, radonImg);
        }
        else
        {
            displayMat_gray(ui.filtered_radon, radonImg);
        }
        
        //��ͶӰ
        iradon_img = iRadon_parallel_line(radonImg);
        displayMat_gray(ui.iradon_img, iradon_img);
    }
    else
    {
        //ͶӰ
        cv::Mat radonImg;
        if (radon_flag == radonFlags::integral)//�����ͼƬ����Ҫ���ֻ���ͶӰ
        {
            int img_longest_edge = max(image.rows, image.cols);
            int sensor_length = ceil(img_longest_edge * sqrt(2)) + 2;
            radonImg = cv::Mat::zeros(sensor_length, radon_num_degree, COLOR_BGR2GRAY);
            for (int theta_degree = 0; theta_degree < radon_num_degree; ++theta_degree)
            {
                Radon_parallel_line_show(image, sensor_length, theta_degree, radonImg);//�Լ�д��ƽ����ͶӰ����
                displayMat_gray(ui.radon_img, radonImg);
                cv::waitKey(6 - reconstruction_speed);
            }
        }
        else//default ����������shepplogan��ͶӰ
        {
            int rows = size_of_shepplogan;
            int rows_radon = rows * sqrt(2);
            radonImg = Mat::zeros(rows_radon, 180, COLOR_BGR2GRAY);
            for (int i = 0; i < 10; ++i)
            {
                Calc_shepplogan_radon_show(size_of_shepplogan, i, radonImg);
                displayMat_gray(ui.radon_img, radonImg);
                cv::waitKey(100 * (6 - reconstruction_speed));
            }
        }

        //�˲�
        if (use_filter == 1)
        {
            radonImg = Filter_radon_parallel_line(radonImg, SL_filter);
            displayMat_gray(ui.filtered_radon, radonImg);
        }
        else
        {
            displayMat_gray(ui.filtered_radon, radonImg);
        }

        //��ͶӰ
        int img_rows = ceil(double(radonImg.rows) / sqrt(2));
        iradon_img = Mat::zeros(img_rows, img_rows, COLOR_BGR2GRAY);
        for (int theta_degree = 0; theta_degree < radonImg.cols; ++theta_degree)
        {
            iRadon_parallel_line_show(radonImg, theta_degree, iradon_img);
            displayMat_gray(ui.iradon_img, iradon_img); 
            cv::waitKey(6 - reconstruction_speed);
        }
        
    }

    double error_d = Calculate_error_sqrt(original_image, iradon_img);
    double error_r = Calculate_error_abs(original_image, iradon_img);

    ui.textBrowser_error_d->setText(QString::number(error_d));
    ui.textBrowser_error_r->setText(QString::number(error_r));

    ui.comboBox_img->setEnabled(true);//�ؽ��������������ø�ѡ�������޸�ͼ��ͶӰ����
    ui.Button_showimg->setEnabled(true);//�ؽ����������������µ�ͼƬ
    ui.pushButton->setEnabled(true);//�ؽ��������ָ�����ʼ�ؽ�����ť����

}

//�޸�ͶӰ���ͣ�������or���ֻ��֣��ġ���ѡ��
void CT_application::comboBox_imgtype_changed(int value)
{
    if (value == 0) //��ѡ��ѡ��Shepploganģ�͡�
    {
        //�������ֻ�����ع���
        ui.horizontalSlider_num_degrees->setEnabled(false);
        ui.spinBox_num_degrees->setEnabled(false);
        //���ý�����ͶӰ��ع���
        ui.horizontalSlider_size_of_shepplogan->setEnabled(true);
        ui.spinBox_size_of_shepplogan->setEnabled(true);

        radon_flag = radonFlags::analytical;  //����������ͶӰ
        image = Create_shepplogan(size_of_shepplogan);
        original_image = image.clone();
    }
    else  // ��ѡ��ѡ�񡰴ӱ���ѡ��ͼƬ��
    {
        //�������ֻ�����ع���
        ui.horizontalSlider_num_degrees->setEnabled(true);
        ui.spinBox_num_degrees->setEnabled(true);
        //���ý�����ͶӰ��ع���
        ui.horizontalSlider_size_of_shepplogan->setEnabled(false);
        ui.spinBox_size_of_shepplogan->setEnabled(false);

        radon_flag = radonFlags::integral;  //���ֻ��ַ�����ͶӰ
        QFileDialog fileDlg;
        //QString imgFile = fileDlg.getOpenFileName(this, QString::fromLocal8Bit("����ͼƬ"), ".", QString::fromLocal8Bit("ͼƬ��*jpg *png *bmp *dfx��"));
        QString imgFile = fileDlg.getOpenFileName(this, QString::fromLocal8Bit("����ͼƬ"), ".", tr("Image Files (*.png *.jpg *.bmp)"));
        QFileInfo fileInf(imgFile);
        if (!fileInf.exists())
        {
            QMessageBox::warning(this, QString::fromLocal8Bit("����"), QString::fromLocal8Bit("���棺��ѡ��ĵ��ļ������ڣ�"), QMessageBox::Ok);
            return;
        }
        else
        {
            std::string pathStr = imgFile.toLocal8Bit();
            image = imread(pathStr.c_str());
        }
        displayMat_color(ui.original_img, image);

        //ui.textBrowser->setText(QString::number(image.type()));
        cvtColor(image, image, COLOR_BGR2GRAY); //��3ͨ��CV_8U��Ϊ��ͨ��CV_8U(uchar����)
        original_image = image.clone();

        //����ͼ��ߴ�Ϊ�����Σ���0��
        if (image.rows != image.cols)
        {
            image = Resize_img_to_Square(image);
        }
    }
}
