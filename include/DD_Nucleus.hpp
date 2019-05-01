//Contains functions related with direct detection rates for DM-Nucleus experiments.
#ifndef __DD_Nucleus_hpp_
#define __DD_Nucleus_hpp_

#include "Numerics_Functions.hpp"
#include "DM_Particle.hpp"

//1. Velocity necessary for a DM particle of mass mChi to recoil a nucleus of mass number A with a given recoil energy.
	extern double vMinimal_N(double RecoilEnergy,double mChi,double A);
//2. Maximum recoil energy a DM particle of mass mChi and speed v can cause on a nucleus of mass number A
	extern double ERMax(double v,double mChi,double A);
//3. Direct detection nuclear recoil spectrum, both analytic and MC
	extern double dRdER(double ER,const DM_Particle& DM,double X,int Z,double A);
	extern Interpolation dRdER(const DM_Particle& DM,double Emin,double Emax,double X,int Z,double A,const std::vector<double>& attenuation,const std::vector<DataPoint> &speeddata);

#endif