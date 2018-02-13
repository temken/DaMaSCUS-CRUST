#ifndef __Density_Estimation_hpp_
#define __Density_Estimation_hpp_

#include <vector>
#include <string>

#include "General_Utilities.hpp"


//Kernels
	extern double Gaussian_Kernel(double x);

//KDE
	extern Interpolation Perform_KDE(std::vector<DataPoint> data,double xMin,double xMax,double bw = 0);

#endif