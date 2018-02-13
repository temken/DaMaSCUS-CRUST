#ifndef __RN_Generators_hpp_
#define __RN_Generators_hpp_

#include <random>
#include <functional>


extern double ProbabilitySample(std::mt19937& PRNG);
extern double MaxwellSample(std::mt19937& PRNG);
extern double ThetaSample(std::mt19937& PRNG);
extern double PhiSample(std::mt19937& PRNG);

extern double Rejection_Sampling(std::function<double(double)> PDF,double xMin,double xMax,double yMin,double yMax,std::mt19937& PRNG);

#endif