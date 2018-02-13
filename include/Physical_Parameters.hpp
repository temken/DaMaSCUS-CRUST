#ifndef __Physical_Parameters_hpp_
#define __Physical_Parameters_hpp_

#include <Eigen/Geometry>
#include <vector>


using namespace std;


//Units:
	//Energy
	extern const double GeV,eV,keV,MeV,TeV;
	//Mass
	extern const double gram,kg;
	//Length
	extern const double cm,meter,km,fm,pb,parsec,kpc,Mpc;
	//Time
	extern const double sec,minute,hour,day,year;
	//Temperature
	extern const double Kelvin;
	//Others
	extern const double erg;

//Specific Parameters:
	//Masses
	extern const double mPlanck,GNewton,mProton,mElectron,mNucleon;
	//Geographic Parameters
	extern const double mEarth,rEarth,rhoEarth,rhoCrust;
	//Solar Parameters
	extern const double mSun,rSun,rhoSun;
	//Dark Matter Halo Parameters
	extern const double v0,vesc,vEarth,rhoDM;
	extern const double Nesc;
	extern const double ymax;
	//degree to rad
	extern const double deg;

//Unit Conversion
	extern double InUnits(double quantity, double dimension);
//Reduced Mass
	extern double Reduced_Mass(double m1,double m2);


//Nucleus Mass
	extern double NucleusMass(double A);

//Wimp Nucleus Cross-section:
	extern double FormFactor(double q,double A);
	extern double sigmaSI(double mX,double sigman0,int Z,double A);
	extern double Total_CrossSection(double mDM,double sigma,int Z,double A,double vDM,bool ldm);
	extern double Differential_CrossSection(double ER,double mDM,double sigma,int Z,double A,double vDM,bool ldm);
//Coordinate Systems
	Eigen::Vector3d SphericalCoordinates(double r,double theta,double phi);

//DM Speed distribution
	extern double SpeedDistribution(double v,double vEarth);


#endif







