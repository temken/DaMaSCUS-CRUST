//Contains functions related with direct detection rates.

#ifndef __Direct_Detection_hpp_
#define __Direct_Detection_hpp_

#include <vector>
#include <string>

#include "General_Utilities.hpp"

//Velocity necessary for a DM particle of mass mChi to recoil a nucleus of mass number A with a given recoil energy.
	extern double vMinimal(double RecoilEnergy,double mChi,double A);
	extern double Cutoff_Speed(std::string experiment,double mDM,int rank=0);
//Maximum recoil energy a DM particle of mass mChi and speed v can cause on a nucleus of mass number A
	extern double ERMax(double v,double mChi,double A);
//Eta-Function: Integral of f(v)/v with lower integration limit v_Min
	extern double EtaFunction(double vMin,double vE);
//Direct detection recoil spectrum
	extern double dRdER(double ER,double X,double A,double mChi,double sigma,bool ldm);
	extern Interpolation dRdER(double Emin,double Emax,double X,double A,double mDM,double sigma,bool ldm,const std::vector<double> attenuation,const std::vector<DataPoint> &speeddata);
//Compute recoil spectrum
	extern Interpolation Compute_Spectrum(std::string experiment,double mDM,double sigma,bool ldm);
	extern Interpolation Compute_Spectrum(std::string experiment,double mDM,double sigma,bool ldm,const std::vector<DataPoint> &data,std::vector<double> attenuation);
//Attenuation
	extern std::vector<double> Compute_Attenuation(unsigned long long int particles, const std::vector<DataPoint> &data,double vCutoff=0.0);
//Likelihood
	extern double Poisson_Likelihood(double expectation_value,unsigned int observed_events);
	extern double Likelihood(std::string experiment,std::function<double(double)> spectrum,double &N);
	extern std::vector<double> Likelihood(std::string experiment,std::function<double(double)> spectrum,std::vector<double> attenuation,std::vector<double> &NMC);
//Total number of events rate
	extern double N_Signals(std::function<double(double)> spectrum,double E1,double E2,double exposure);
//Standard Upper bound
	extern double Upper_Bound(std::string experiment,double mDM,bool ldm,double CL);

#endif