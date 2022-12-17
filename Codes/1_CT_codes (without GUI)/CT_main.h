#pragma once
#ifndef CT_MAIN_H_INCLUDED
#define CT_MAIN_H_INCLUDED

using namespace std;
using namespace cv;

enum radonFlags {
	/** choose the method for radon */
	analytical = 1, // for shepplogan
	integral = 2 // for existing image
};

enum radonTypes {
	parallel = 1,
	sector = 2
};

#include "basic_display_functions.h"
#include "radon_parallel_line.h"
#include "radon_fan_beam.h"
#include "filter_parallel_line.h"
#include "shepplogan.h"


#endif // !CT_MAIN_H_INCLUDED#pragma once
