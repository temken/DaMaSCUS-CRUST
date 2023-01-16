#include "DD_Nucleus.hpp"

#include <cmath>
#include <iostream>

#include "DD_General.hpp"
#include "Physics_Functions.hpp"

// 1. Velocity necessary for a DM particle of mass mChi to recoil a nucleus of mass number A with a given recoil energy.
double vMinimal_N(double RecoilEnergy, double mChi, double A)
{
	return sqrt(mNucleon * A * RecoilEnergy / 2.0 / pow(Reduced_Mass(mChi, A * mNucleon), 2.0));
}
// 2. Maximum recoil energy a DM particle of mass mChi and speed v can cause on a nucleus of mass number A
double ERMax(double v, double mChi, double A)
{
	return 2.0 * v * v * pow(Reduced_Mass(mChi, A * mNucleon), 2.0) / A / mNucleon;
}

// 3. Direct detection recoil spectrum
double dRdER(double ER, const DM_Particle& DM, double X, int Z, double A)
{
	double mA  = A * mNucleon;
	double nDM = rhoDM / DM.mass;
	double v   = 1.0;	// Cancels in the following expression
	return X / mA * nDM * (v * v * DM.dSdER(ER, Z, A, v)) * EtaFunction(vMinimal_N(ER, DM.mass, A), vEarth);
}
Interpolation dRdER(const DM_Particle& DM, double Emin, double Emax, double X, int Z, double A, const std::vector<double>& attenuation, const std::vector<DataPoint>& speeddata, double vCutoff)
{
	// 0. maximal energy
	double ERmax = ERMax((vesc + vEarth), DM.mass, A);
	ERmax		 = std::min(ERmax, Emax);	// only calculate points where dRdER !=0 for the interpolation
	// Check if the DM is able to cause recoils
	double vMin_tot = vMinimal_N(Emin, DM.mass, A);
	if(vMin_tot > (vesc + vEarth))
	{
		int points = 10;
		double dER = (Emax - Emin) / (points - 1.0);
		std::vector<std::vector<double>> interpol_list;
		for(int i = 0; i < points; i++)
		{
			double ER = Emin + i * dER;
			interpol_list.push_back(std::vector<double> {ER, 0.0});
		}
		return Interpolation(interpol_list);
	}

	// 1. Prefactor
	double mA		 = A * mNucleon;
	double prefactor = X / mA * rhoDM / DM.mass;
	// 2. Find the Kernel density estimate for the speed distribution function
	Interpolation etaMC = EtaFunction(attenuation, speeddata, vCutoff, vEarth);

	// 3. Create list to interpolate
	int points = 200;
	double dER = (ERmax - Emin) / (points - 1.0);
	std::vector<std::vector<double>> interpol_list;
	for(int i = 0; i < points; i++)
	{
		double ER	= Emin + i * dER;
		double vMin = vMinimal_N(ER, DM.mass, A);
		if(vMin >= (vesc + vEarth))
			interpol_list.push_back(std::vector<double> {ER, 0.0});
		else
		{
			double v  = 1.0;
			double dr = prefactor * (v * v * DM.dSdER(ER, Z, A, v)) * etaMC(vMin);
			interpol_list.push_back(std::vector<double> {ER, dr});
		}
	}
	// 3.4 If ERMax < Emax, fill the rest of the interval with zeros, so that we have a interpolation on the full interval.
	if(ERmax < Emax)
	{
		points = 5;
		dER	   = (Emax - ERmax) / (points - 1);
		for(int i = 1; i < points; i++)
		{
			double ER = ERmax + i * dER;
			interpol_list.push_back(std::vector<double> {ER, 0.0});
		}
	}
	// 4. Interpolate and include prefactor.
	Interpolation drder(interpol_list);

	return drder;
}