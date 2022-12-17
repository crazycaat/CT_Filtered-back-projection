#pragma once
#pragma once
#ifndef SHEPPLOGAN_H_INCLUDED
#define SHEPPLOGAN_H_INCLUDED

using namespace std;
using namespace cv;

#define  pi  3.14159265358979323846  //!< const value,pi

Mat Create_shepplogan(int rows = 256);
Mat Calc_shepplogan_radon(int rows = 256);

#endif // !SHEPPLOGAN_H_INCLUDED