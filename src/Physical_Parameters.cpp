#include "Physical_Parameters.hpp"

#include <cmath>
#include <iostream>
#include <functional>

#include "General_Utilities.hpp"

using namespace std::placeholders;

//Units:
	//Energy
	const double GeV=1.0;
	const double eV=1.0E-9*GeV;
	const double keV=1.0E-6*GeV;
	const double MeV=1.0E-3*GeV;
	const double TeV=1.0E3*GeV;
	//Mass
	const double gram=5.617977528089887E23*GeV;
	const double kg=1E3*gram;
	//Length
	const double cm=5.068E13/GeV;
	const double meter=100*cm;
	const double km=1000*meter;
	const double fm=1E-15*meter;
	const double pb=1E-36*pow(cm,2);
	const double parsec=3.0857E16*meter;
	const double kpc=1E3*parsec;
	const double Mpc=1E6*parsec;
	//Time
	const double sec=299792458*meter;
	const double minute=60*sec;
	const double hour=60*minute;
	const double day=24*hour;
	const double year=365*day;
	//Temperature
	const double Kelvin=8.62E-14*GeV;
	//Others
	const double erg=gram*pow(cm/sec,2);

//Specific Parameters:
	//Masses
	const double mPlanck= 1.2209E19*GeV;
	const double GNewton=pow(mPlanck,-2);
	const double mProton=0.938*GeV ;
	const double mElectron= 0.511*MeV;
	const double mNucleon = 0.932*GeV;
	//Geographic Parameters
	const double mEarth=5.972E24*kg;
	const double rEarth=6371*km;
	const double rhoEarth=5.51*gram*pow(cm,-3);
	const double rhoCrust=2.7*gram*pow(cm,-3);
	//Solar Parameters
	const double mSun=1.989E30*kg;
	const double rSun=6.957E8*meter;
	const double rhoSun=151*gram*pow(cm,-3);
	//Dark Matter Halo Parameters
	const double v0=220*km/sec;
	const double vesc=544*km/sec;
	const double vEarth = 230*km/sec;
	const double rhoDM=0.3*GeV*pow(cm,-3);
	const double Nesc=M_PI*v0*v0*(sqrt(M_PI)*v0*erf(vesc/v0)-2*vesc*exp(-vesc*vesc/v0/v0));
	const double ymax=4*M_PI/Nesc*v0*v0*exp(-1);
	//degree to rad
	const double deg=M_PI/180.0;



//Unit Conversion
	double InUnits(double quantity, double dimension)
	{
		return quantity/dimension;
	}


//Physical Functions

	//Nucleus Mass
		double NucleusMass(double A)
		{
			return A*mNucleon;
		}
	//Reduced Mass
		double Reduced_Mass(double m1,double m2)
		{
			return m1*m2/(m1+m2);
		}


//DM nucleus cross-section
	//Squared Helm form factor
		double FormFactor(double q,double A)
		{
			if(q==0.0) return 1.0;
			else
			{
				double a = 0.52*fm;
				double c = (1.23*pow(A,1.0/3.0)-0.6)*fm;
				double s = 0.9*fm;
				double rn = sqrt(c*c+7.0/3.0*pow(M_PI*a,2.0)-5*s*s);
				double qr = q*rn;
				return 3.0*(sin(qr)/pow(qr,3.0)-cos(qr)/pow(qr,2.0))*exp(-q*q*s*s/2.0);
			}
			
		}
	//Zero momentum transfer spin-independent cross-section
		double sigmaSI(double mX,double sigman0,int Z,double A)
		{
			return sigman0*pow(Reduced_Mass(mX,NucleusMass(A)),2.0)/pow(Reduced_Mass(mX,mNucleon),2.0)*pow(A,2.0);
		}
	//Differential cross-section
		double Differential_CrossSection(double ER,double mDM,double sigma,int Z,double A,double vDM,bool ldm)
		{
			double mA = A*mNucleon;
			double ERmax = 2.0*pow(Reduced_Mass(mDM,mA)*vDM,2.0)/mA;
			if(ldm) return sigmaSI(mDM,sigma,Z,A)/ERmax;
			else
			{
				double q = sqrt(2.0*mA*ER);
				return sigmaSI(mDM,sigma,Z,A)/ERmax*pow(FormFactor(q,A),2.0); 
			}
			
		}
		double Total_CrossSection(double mDM,double sigma,int Z,double A,double vDM,bool ldm)
		{
			if(ldm) return sigmaSI(mDM,sigma,Z,A);
			else
			{
				double mA = A*mNucleon;
				double ERmax = 2.0*pow(Reduced_Mass(mDM,mA)*vDM,2.0)/mA;
				//integrate the diff cross section
				auto dodER = std::bind(Differential_CrossSection,_1,mDM,sigma,Z,A,vDM,ldm);
				//Numerical integration
					double a = 0.0;
					double b = ERmax;
					double c = (a + b)/2;
					double h = b - a;                                                                  
		  			double fa = dodER(a);
		  			double fb = dodER(b);
		  			double fc = dodER(c);
		  			double S = (h/6)*(fa + 4*fc + fb);
		  			double epsilon = 1e-5*S;
		  			double integral =AdaptiveSimpsons(dodER,a,b,epsilon);
					return integral;
			}
		}

//Coordinate Systems
	Eigen::Vector3d SphericalCoordinates(double r,double theta,double phi)
	{
		Eigen::Vector3d v(r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta));
		return  v;
	}
//DM speed distribution
	double StepFunction(double x)
	{
		if(x>=0) 	return 1.0;
		else		return 0.0;
	}
	double SpeedDistribution(double v,double vEarth)
	{
		return M_PI*v*v0*v0/Nesc/vEarth*(2*exp(-(v*v+vEarth*vEarth)/v0/v0)*sinh(2*v*vEarth/v0/v0)+(exp(-pow(v+vEarth,2.0)/v0/v0)-exp(-vesc*vesc/v0/v0))*StepFunction(abs(v+vEarth)-vesc)-(exp(-pow(v-vEarth,2.0)/v0/v0)-exp(-vesc*vesc/v0/v0))*StepFunction(abs(v-vEarth)-vesc) );
	}



	