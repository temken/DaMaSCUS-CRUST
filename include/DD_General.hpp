//Contains general functions related with direct detection rates.
#ifndef __DD_General_hpp_
#define __DD_General_hpp_

#include <vector>
#include <string>

#include "Numerics_Functions.hpp"
#include "DM_Particle.hpp"

//1. Attenuation
	extern std::vector<double> Compute_Attenuation(unsigned long long int particles, const std::vector<DataPoint> &data,double vCutoff);

//2. Eta-Function: Integral of f(v)/v with lower integration limit v_Min
	extern double EtaFunction(double vMin,double vE);
	extern Interpolation EtaFunction(const std::vector<double>& attenuation,const std::vector<DataPoint> &speeddata,double vCutoff,double vE);

//3. Compute recoil spectrum
	extern Interpolation Compute_Spectrum(const DM_Particle &DM,const std::string& experiment);
	extern Interpolation Compute_Spectrum(const DM_Particle &DM,const std::string& experiment,const std::vector<DataPoint> &data,const std::vector<double>& attenuation);

//4. Likelihood
	extern double Likelihood(const DM_Particle& DM,const std::string& experiment,const std::function<double(double)>& spectrum,double &N);
	extern std::vector<double> Likelihood(const DM_Particle& DM,const std::string& experiment,const std::function<double(double)>& spectrum,const std::vector<double>& attenuation,const std::vector<DataPoint> &data,double vCutoff,std::vector<double> &NMC);

//5. Total number of events rate
	extern double N_Signals(const std::function<double(double)>& spectrum,double E1,double E2,double exposure);

//6. Standard constraints
	extern double Upper_Bound(const DM_Particle& DM,const std::string& experiment,double CL);
	extern std::vector<std::vector<double>> Limit_Curve(DM_Particle DM,const std::string& experiment,const std::string& simID,double mMin,double mMax,int Nm,double CL,int rank=0);

#endif