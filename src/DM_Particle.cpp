#include "DM_Particle.hpp"

#include <iostream>
#include <cmath>
#include <functional>

#include "Numerics_Functions.hpp"
#include "Physics_Functions.hpp"

//6. DM-nucleus scattering cross sections:
	//Nuclear Helm form factor
		double FormFactor_N(double q,double A,bool ldm)
		{
			if(ldm || q==0.0) return 1.0;
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
		
	//3. DM particle struct
		DM_Particle::DM_Particle()
		{
			mass=0.0;
			sigma_n=0.0;
			sigma_e = 0.0;
			ldm = true;
			formfactor=" ";
			ZorA = " ";
			screening=false;
			mMediator=0.0;
		}

		DM_Particle::DM_Particle(double mDM,double sn,double se,bool light,std::string ff,std::string za,bool scr,double mMed)
		{
			mass=mDM;
			sigma_n=sn;
			sigma_e = se;
			ldm = light;
			formfactor=ff;
			ZorA = za;
			screening=scr;
			mMediator=mMed;
		}
		//Atomic form factor (charge screening)
		double DM_Particle::FormFactor_A(double q,int Z) const
		{
			if(!screening) return 1.0;
			else
			{
				double a2q2 = pow(q*Thomas_Fermi_Radius(Z),2.0);
				return a2q2/(1.0+a2q2);
			}
		}
		void DM_Particle::Set_Mass(double m)
		{
			mass = m;
		}
		void DM_Particle::Set_Sigma_n(double s)
		{
			sigma_n = s;
			sigma_e = pow(Reduced_Mass(mass,mElectron)/Reduced_Mass(mass,mProton),2.0)*sigma_n;
		}
		void DM_Particle::Set_Sigma_e(double s)
		{
			sigma_e = s;
			sigma_n = pow(Reduced_Mass(mass,mProton)/Reduced_Mass(mass,mElectron),2.0)*sigma_e;
		}
	//DM form factor
		double DM_Particle::FormFactor(double q) const 
		{
			//Contact interactions
				if(formfactor=="Contact") return 1.0;
			//General dark photon
				else if(formfactor=="General")	return (qRef*qRef+mMediator*mMediator)/(q*q+mMediator*mMediator);
			//Long range interaction
				else if (formfactor=="Long-Range")	return qRef*qRef/q/q;
			//Electric dipole interaction
				else if (formfactor=="Electric-Dipole")		return qRef/q;
			//Error
				else
				{
					std::cerr <<"Error in FormFactor(): Form factor "<<formfactor <<"not recognized."<<endl;
					std::exit(EXIT_FAILURE);
				}
		}
	//Zero momentum transfer spin-independent cross-section
		double DM_Particle::sigmaSI(int Z,double A) const
		{
			double X=-1.0;
			double mAux=-1.0;
			if(ZorA=="Z") 
			{
				X=Z;
				mAux = mProton;
			}
			else if (ZorA=="A") 
			{
				X=A;
				mAux = mNucleon;
			}
			else 
			{
				std::cerr <<"Error in sigmaSI: ZorA = " <<ZorA <<" not recognized."<<endl;
				std::exit(EXIT_FAILURE);
			}
			return sigma_n*pow(Reduced_Mass(mass,NucleusMass(A)),2.0)/pow(Reduced_Mass(mass,mAux),2.0)*pow(X,2.0);
		}
	//Differential cross-sections
		double DM_Particle::dSdER(double ER,int Z,double A,double vDM) const
		{
			double mA = A*mNucleon;
			double ERmax = 2.0*pow(Reduced_Mass(mass,mA)*vDM,2.0)/mA;
			double q = sqrt(2.0*mA*ER);	
			if(q==0&&(formfactor=="Electric-Dipole"||formfactor=="Long-Range")&&screening) q=1e-15*keV; //Screening cause 1/q to cancel
			return sigmaSI(Z,A) / ERmax * pow(FormFactor_N(q,A,ldm),2.0) * pow(FormFactor(q),2.0) * pow(FormFactor_A(q,Z),2.0); 
		}
		double DM_Particle::dSdq2(double q,int Z,double A,double vDM) const
		{
			double mA = A*mNucleon;
			double qMax2 = 4.0*pow(Reduced_Mass(mass,mA)*vDM,2.0);
			if(q==0&&(formfactor=="Electric-Dipole"||formfactor=="Long-Range")&&screening) q=1e-15*keV; //Screening cause 1/q to cancel
			return sigmaSI(Z,A) / qMax2 * pow(FormFactor_N(q,A,ldm),2.0) * pow(FormFactor(q),2.0)* pow(FormFactor_A(q,Z),2.0); 
		}
		double DM_Particle::Sigma_Tot(int Z,double A,double vDM) const
		{
			if(!ldm)
			{
				double mA = A*mNucleon;
				double ERmax = 2.0*pow(Reduced_Mass(mass,mA)*vDM,2.0)/mA;
				double ERmin = 0.0;
				//integrate the diff cross section
					auto dodER = std::bind(&DM_Particle::dSdER,this,std::placeholders::_1,Z,A,vDM);
					//Numerical integration
			  			double epsilon = Find_Epsilon(dodER,ERmin,ERmax,1e-5);
			  			double integral =Integrate(dodER,ERmin,ERmax,epsilon);
					return integral;
			}
			else
			{
				double result = sigmaSI(Z,A);
				
				double q2max = pow(2.0*Reduced_Mass(mass,A*mNucleon)*vDM,2.0);
				double a2=pow(Thomas_Fermi_Radius(Z),2.0);
				//General interaction
					if(formfactor == "General")
					{
						double m2 = pow(mMediator,2.0);
						if(screening)
						{
							result*=(pow(a2,2.0)*pow(m2 + qRef*qRef,2.0)*(((a2*m2-1.0)*q2max*(q2max + m2*(2.0 + a2*q2max)))/((m2 + q2max)*(1.0 + a2*q2max)) + 2.0*m2*log((m2 + q2max)/(m2 + a2*m2*q2max))))/(pow(a2*m2-1.0,3.0)*q2max);
						}
						else
						{
							result*= pow(qRef*qRef+m2,2.0)/m2/(m2+q2max);
						}
					}
				//Contact interaction
					else if(formfactor == "Contact" && screening)
					{
						result *=(1.0+1.0/(1+a2*q2max)-2.0/a2/q2max*log1p(a2*q2max));
					}
				//Electric dipole interaction
					else if(formfactor == "Electric-Dipole")
					{
						result *= pow(qRef,2.0)/q2max*(log1p(a2*q2max)-a2*q2max/(1+a2*q2max));
					}
				//Long range interaction
					else if(formfactor== "Long-Range")
					{
						result *= pow(a2,2.0)*pow(qRef,4.0)/(1+a2*q2max);
					}

				return result;
			}
		}

	//Find the stopping cross section
		double Stopping_CrossSection(const std::string& detector,double sigma,double mDM)
		{
			if(detector == "Semiconductor"||detector=="SENSEI-surface"||detector=="SENSEI"||detector=="SuperCDMS"||detector=="DAMIC-M"||detector == "XENON10e"||detector == "XENON100e"||detector =="DarkSide-50")return pow(Reduced_Mass(mDM,mProton)/Reduced_Mass(mDM,mElectron),2.0)*sigma;
			else return sigma;
		}

//7. DM speed distribution
	double SpeedDistribution(double v,double vEarth)
	{
		return M_PI*v*v0*v0/Nesc/vEarth*(2*exp(-(v*v+vEarth*vEarth)/v0/v0)*sinh(2*v*vEarth/v0/v0)+(exp(-pow(v+vEarth,2.0)/v0/v0)-exp(-vesc*vesc/v0/v0))*StepFunction(abs(v+vEarth)-vesc)-(exp(-pow(v-vEarth,2.0)/v0/v0)-exp(-vesc*vesc/v0/v0))*StepFunction(abs(v-vEarth)-vesc) );
	}
 	double Average_Speed(double vEarth,double vMin)
 	{
 		//1. integrand
	 		std::function<double(double)> integrand = [vEarth] (double v)
	 		{
	 			return v*SpeedDistribution(v,vEarth);
	 		};
	 	//2. integrate
	 		double vMean = Integrate(integrand,vMin,(vesc+vEarth),1e-6);
	 	//3. Renormalize
	 		if(vMin>0.0)
	 		{
	 			
	 			std::function<double(double)> integrand2 = [vEarth] (double v)
		 		{
		 			return SpeedDistribution(v,vEarth);
		 		};
		 		double norm=Integrate(integrand2,vMin,vesc+vEarth,1e-6);
		 		vMean/=norm;
	 		}
	 		return vMean;
 	}