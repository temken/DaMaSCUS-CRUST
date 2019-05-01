//Contains functions related with direct detection rates for DM-Electron experiments.
#ifndef __DD_Electron_hpp_
#define __DD_Electron_hpp_

#include <string>
#include <vector>

#include "DM_Particle.hpp"
#include "Numerics_Functions.hpp"

//1. Import crystal/ionization form factor
	extern void Import_FormFactor(const std::string& target);

//2. Semiconductors

	//Minimal velocity
	extern double vMinimal_e(double q,double Ee,double mDM);

	//Recoil spectra
	extern double dRdEe(double Ee,const DM_Particle& DM,const std::string& target,double efficiency);
	extern Interpolation dRdEe(const DM_Particle& DM,const std::string& target,double efficiency,const std::vector<double>& attenuation,const std::vector<DataPoint> &speeddata,double vCutoff);
	//Total event rates

	extern double TotalRate(const DM_Particle& DM,int Qthreshold,double efficiency,const std::string& target);
	extern double TotalRate(const DM_Particle& DM,int Qthreshold,double efficiency,const std::string& target,const std::vector<double>& attenuation,const std::vector<DataPoint> &speeddata,double vCutoff);

//3. Atomic shell struct for liquid noble targets.
	struct Shell 
	{
		std::string name;
		double Ebinding;
		int neSecondary;
		std::vector< std::vector<double> > FormFactor;
		double logkMin,logkMax,logqMin,logqMax, dlogk,dlogq;
		int nk,nq;
			//Constructors
				Shell(std::string nm,std::string filepath,double Eb,int Ns,double lnkMin,double lnkMax,int Nk,double lnqMin,double lnqMax,int Nq);
	};
	extern std::vector<Shell> Shells;

//4. Liquid noble gas experiments

	//Recoil energy spectrum
	extern double dRdlog10Ee(double Ee,const DM_Particle& DM,const std::string& target,const Shell& shell);

	//Spectrum number of electrons
	extern double PDFne(const std::string& target,unsigned int ne,double Ee,const Shell& shell);
	extern double dRdne(unsigned int ne,const DM_Particle& DM,const std::string& target,const Shell& shell);
	extern std::vector<double> dRdne(const DM_Particle& DM,const std::string& target,const std::vector<Shell>& shells,const std::vector<double> &attenuation,const std::vector<DataPoint> &speeddata,double vCutoff);
	
	//Spectrum number of PE
	extern double dRdnPE(unsigned int nPE,const DM_Particle& DM,const std::string& target,double muPE,double sigPE,const std::vector<Shell>& shells);
	extern std::vector<double> dRdnPE(double muPE,double sigPE,const std::vector<double>& dRdn);
	extern double dNdnPE(unsigned int nPE,const DM_Particle& DM,const std::string& target,const std::string& detector,const std::vector<Shell>& shells);
	extern std::vector<double> dNdnPE(const std::string& detector,const std::vector<double>& dRdnpe);

#endif