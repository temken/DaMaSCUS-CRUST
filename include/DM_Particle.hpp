//Contains functions regarding dark matter. Especially scattering cross sections.
#ifndef __DM_Particle_hpp_
#define __DM_Particle_hpp_

#include <string>

//1. DM-nucleus scattering cross sections:
	extern double FormFactor_N(double q,double A,bool ldm);
	class DM_Particle
	{
		private:
			std::string ZorA;
			double FormFactor_A(double q,int Z)const;
		public:
			bool ldm;
			std::string formfactor;
			double mMediator;
			bool screening;
			double mass;
			double sigma_n;
			double sigma_e;

			//Constructor
			DM_Particle();
			DM_Particle(double mDM,double sn,double se=0.0,bool light=true,std::string ff="Contact",std::string za="A",bool scr=false,double mMed=0.0);
			
			void Set_Mass(double m);
			void Set_Sigma_n(double s);
			void Set_Sigma_e(double s);
			
			double FormFactor(double q)const;
			double sigmaSI(int Z,double A)const;
			double dSdER(double ER,int Z,double A,double vDM) const;
			double dSdq2(double q,int Z,double A,double vDM)const;
			double Sigma_Tot(int Z,double A,double vDM) const;
	};

	//Return the relevant stopping cross section
	extern double Stopping_CrossSection(const std::string& detector,double sigma,double mDM);

//2. DM Speed distribution
	extern double SpeedDistribution(double v,double vEarth);
	extern double Average_Speed(double vEarth,double vMin=0.0);

#endif