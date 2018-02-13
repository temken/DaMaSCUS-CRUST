#include "Direct_Detection.hpp"

#include <iostream>
#include <fstream>
#include <functional>
#include <algorithm>
#include <numeric> //for std::accumulate

#include "Physical_Parameters.hpp"
#include "Density_Estimation.hpp"


using namespace std::placeholders;
//Statistics
	//Gaussian
		double Gauss_Distribution(double x,double mu,double sigma)
		{
			return 1.0/sqrt(2.0*M_PI)/sigma*exp(-pow((x-mu)/sigma,2.0)/2.0);
		}
	//Poisson statistics 
		double Factorial(unsigned int n)
		{
			if(n==0) return 1;
			else return n*Factorial(n-1);
		}
		double Poisson_Distribution(double expected_events, unsigned int events)
		{
			double result = pow(expected_events,events)/Factorial(events)*exp(-expected_events);
			if(std::isnan(result)) result = 0.0; //pow() then gives inf resulting in nan 
			return result; 
		}
		double Poisson_Likelihood(double expectation_value,unsigned int observed_events)
		{
			double llh = 0.0;
			for(unsigned int k = 0; k<=observed_events;k++) llh+=Poisson_Distribution(expectation_value,k);
			return llh;
		}
	//Maximum Gap a la Yellin
		double C0(double x,double mu)
		{
			if(x==mu) return 1.0-exp(-mu);
			else
			{
				int m = mu/x;
				double sum=0.0;
				for(int k=0;k<=m;k++) 
				{
					double term = pow(k*x-mu,k)/Factorial(k)*exp(-k*x)*(1.0+k/(mu-k*x));
					sum+= term ;
					if (fabs(term)<1e-20) break;
				}
				return sum;
			}
		}
		double C0_error(double x,double mu,double mu_error)
		{
			if(x==mu) return mu_error*exp(-mu);
			else
			{
				int m = mu/x;
				double sum=0.0;
				for(int k=0;k<=m;k++) 
				{
					double term = k*pow(k*x-mu,(-2.0+k))/Factorial(k)*exp(-k*x)*(1.0+k*(x-1.0)-mu);
					sum+= term ;
					if (fabs(term)<1e-20) break;
				}
				return fabs(sum*mu_error);
			}
		}

//Velocity necessary for a DM particle of mass mChi to recoil a nucleus of mass number A with a given recoil energy.
	double vMinimal(double RecoilEnergy,double mChi,double A)
	{
		return sqrt(mNucleon*A*RecoilEnergy/2.0/pow(Reduced_Mass(mChi,A*mNucleon),2.0));
	}
	double Cutoff_Speed(std::string experiment,double mDM,int rank)
	{
		double vcut=0;
		if(experiment=="DAMIC")
		{
			double Ethr = 0.55*keV;
			double A = 28.0;
			vcut = vMinimal(Ethr,mDM,A);
		}
		else if(experiment == "XENON1T")
		{
			double Ethr = 5*keV;
			double A = 131.0;
			vcut = vMinimal(Ethr,mDM,A);
		}
		else if(experiment =="CRESST-II")
		{
			double Ethr = 307*eV;
			double resolution = Ethr/5.0;
			double Emin = Ethr-3.0*resolution;
			double A[3] = {16.,40.,184.};
			std::vector<double> v = {vMinimal(Emin,mDM,A[0]),vMinimal(Emin,mDM,A[1]),vMinimal(Emin,mDM,A[2])};
			vcut = *std::min_element(v.begin(),v.end());
		}
		else if(experiment =="CRESST-surface")
		{
			double Ethr = 19.7*eV;
			double resolution = 3.74*eV;
			double Emin = Ethr-3.0*resolution;
			double A[3] = {16.,27.};
			std::vector<double> v = {vMinimal(Emin,mDM,A[0]),vMinimal(Emin,mDM,A[1])};
			vcut = *std::min_element(v.begin(),v.end());
		}
		else
		{
			//Error message
				cout <<"Error in Cutoff_Speed(): Experiment " <<experiment <<" not recognized."<<endl;
				vcut = 0;
		}

		if(vcut>(vEarth+vesc))
		{
			if(rank==0) cout <<"Warning: vCutoff>vEsc+vEarth for mDM = " <<mDM <<"GeV."<<endl;
			vcut = -1;
		}
		return vcut;
	}

//Maximum recoil energy a DM particle of mass mChi and speed v can cause on a nucleus of mass number A
	double ERMax(double v,double mChi,double A)
	{
		return 2.0*v*v*pow(Reduced_Mass(mChi,A*mNucleon),2.0)/A/mNucleon;
	}

//Eta-Function: Integral of f(v)/v with lower integration limit v_Min
	double EtaFunction(double vMin,double vE)
	{
		double xMin = vMin/v0;
		double xEsc = vesc/v0;
		double xE=vE/v0;
		if(xMin>xE+xEsc) return 0.0;
		else if(fabs(xMin-xE-xEsc)<1e-8) return 0.0;
		else if(xMin>abs(xE-xEsc)) return pow(M_PI,1.5)*v0*v0/2.0/Nesc/xE*(erf(xEsc)-erf(xMin-xE)-2.0/sqrt(M_PI)*(xE+xEsc-xMin)*exp(-xEsc*xEsc));
		else if (xEsc>xE) return pow(M_PI,1.5)*v0*v0/2.0/Nesc/xE*(erf(xMin+xE)-erf(xMin-xE)-4.0/sqrt(M_PI)*xE*exp(-xEsc*xEsc));
		else return 1.0/v0/xE;
	}

//Direct detection recoil spectrum
	double dRdER(double ER,double X,double A,double mChi,double sigma,bool ldm)
	{
		double output = X/2.0*rhoDM/mChi*sigma*pow(A/Reduced_Mass(mChi,mNucleon),2.0)*EtaFunction(vMinimal(ER,mChi,A),vEarth);
		if(!ldm)
		{
			double q = sqrt(2.0*A*mNucleon*ER);
			output *= pow(FormFactor(q,A),2.0);
		}
		return output;
	}

	Interpolation dRdER(double Emin,double Emax,double X,double A,double mDM,double sigma,bool ldm,const std::vector<double> attenuation,const std::vector<DataPoint> &speeddata)
	{
		//0. maximal energy
			double ERmax = ERMax((vesc+vEarth),mDM,A);
			ERmax = std::min(ERmax,Emax); //only calculate points where dRdER !=0 for the interpolation
		//1. Prefactor
			double mA = A * mNucleon;
			double prefactor = X/mA*rhoDM/mDM*attenuation[0];
		//2. Find the Kernel density estimate for the speed distribution function
			double vCutoff=vMinimal(Emin,mDM,A);
			Interpolation kde = Perform_KDE(speeddata,vCutoff,(vesc+vEarth));
			// kde.Save_Function("kdetest.txt",100);
		//3. Create list to interpolate
			int points = 200;
			double dER = (ERmax-Emin)/(points-1.0);
			std::vector<std::vector<double>> interpol_list;
			for(int i = 0;i<points;i++)
			{
				double ER = Emin+ i*dER;
				//3.1 Create integrand.
					std::function<double(double)> integrand = [ER,mDM,sigma,A,&kde,ldm](double v)
					{
						return v*kde(v)*Differential_CrossSection(ER,mDM,sigma,1,A,v,ldm);
					};
				//3.2 Integrate.
					double vMin = vMinimal(ER,mDM,A);
					double integral;
					if(vMin>=(vesc+vEarth))
					{
						integral = 0.0;
					}
					else
					{
				  		//Integrate
				  		//NOTE: Gives precision warnings regularly due to a too small choice of epsilon. Needs to be fixed at some point.
				  			double epsilon = Find_Epsilon(integrand,vMin,(vesc+vEarth),1.0e-3);
							integral = AdaptiveSimpsons(integrand,vMin,(vesc+vEarth),epsilon);
	
					}
				//3.3 Append to list 
					interpol_list.push_back(std::vector<double> {ER,integral});
			}
			//3.4 If ERMax < Emax, fill the rest of the interval with zeros, so that we have a interpolation on the full interval.
				if(ERmax<Emax)
				{
					points = 5;
					dER = (Emax-ERmax)/(points-1);
					for(int i = 1;i<points;i++)
					{
						double ER = ERmax + i*dER;
						interpol_list.push_back(std::vector<double> {ER,0.0});
					}
				}
		//4. Interpolate and include prefactor.
			Interpolation drder(interpol_list);
			drder.Multiply(prefactor);
		
		return drder;
	}

//Compute recoil spectrum
	Interpolation Compute_Spectrum(std::string experiment,double mDM,double sigma,bool ldm)
	{
		if(experiment == "DAMIC")
		{
			//1. Detector input
				// double efficiency = 0.75;
				double efficiency = 1.0;
				double A=28.0;
				double X=1.0;
				double Ethr = 0.55*keV;
				double Emax = 7.0*keV;
			//2. Tabulate and interpolate
				//2.1 Find maximum ER, i.e. the domain of the spectrum
					double ERmax = ERMax((vesc+vEarth),mDM,A);
					ERmax = std::min(ERmax,Emax);
				//2.2 Tabulate the spectrum
					int points = 250;
					double dER = (ERmax-Ethr)/(points-1.0);
					std::vector<std::vector<double>> interpol_list;
					for(int i = 0;i<points;i++)
					{
						double ER = Ethr + i*dER;
						double dR = efficiency*dRdER(ER,X,A,mDM,sigma,ldm);
						interpol_list.push_back(std::vector<double> {ER,dR});
					}
				//2.3 if the maximum value is not Emax we append one more point, to avoid problems for the integration later on.
					if(ERmax<Emax)interpol_list.push_back(std::vector<double> {Emax,0.0});
				//2.4 Interpolate and return
					Interpolation spectrum(interpol_list);
					return spectrum;
				
		}
		else if(experiment=="XENON1T")
		{
			//1. Detector input
				double efficiency = 0.82;
				double A=131.0;
				double X=1.0;
				double Ethr = 5*keV;
				double Emax = 40.0*keV;
			//2. Tabulate and interpolate
				//2.1 Find maximum ER, i.e. the domain of the spectrum
					double ERmax = ERMax((vesc+vEarth),mDM,A);
					ERmax = std::min(ERmax,Emax);
				//2.2 Tabulate the spectrum
					int points = 250;
					double dER = (ERmax-Ethr)/(points-1.0);
					std::vector<std::vector<double>> interpol_list;
					for(int i = 0;i<points;i++)
					{
						double ER = Ethr + i*dER;
						double dR = efficiency*dRdER(ER,X,A,mDM,sigma,ldm);
						interpol_list.push_back(std::vector<double> {ER,dR});
					}
				//2.3 if the maximum value is not Emax we append one more point, to avoid problems for the integration later on.
					if(ERmax<Emax)interpol_list.push_back(std::vector<double> {Emax,0.0});
				//2.4 Interpolate and return
					Interpolation spectrum(interpol_list);
					return spectrum;
		}
		else if(experiment=="CRESST-II")
		{
			//Detector Input
				double A[3] = {16.,40.,184.};
				double X[3] = {0.22,0.14,0.64};				
				double Ethr = 307*eV;
				double resolution = Ethr/5.0;
				double Emax = 40*keV;
			//1. Import efficiencies
				Interpolation eff_O("../detectors/CRESST-II/Lise_eff_AR_O.dat",keV);
				Interpolation eff_Ca("../detectors/CRESST-II/Lise_eff_AR_Ca.dat",keV);
				Interpolation eff_W("../detectors/CRESST-II/Lise_eff_AR_W.dat",keV);

			//2. Tabulate dR/dE
				//2.1 Find the maximum E for which we have contributions.
					std::vector<double> ERmax_aux = {ERMax((vesc+vEarth),mDM,A[0]),ERMax((vesc+vEarth),mDM,A[1]),ERMax((vesc+vEarth),mDM,A[2])};
					double ERmax = *std::max_element(ERmax_aux.begin(),ERmax_aux.end());
					if((ERmax+6.0*resolution)<Emax) Emax =  (ERmax+6.0*resolution);
				//2.2 Tabulate dRdE
					std::vector<std::vector<double>> table;
					int steps = 200;
					double dE = (Emax - Ethr)/(steps-1.0);

					for(int i =0; i<steps; i++)
					{
						double E = Ethr+ i*dE;
						//2.2.1 Find minimum and maximum ER contributing to dR/dE(E):
							std::vector<double> aux = {E-6.0*resolution,Ethr-3.0*resolution};
							double eMin = *std::max_element(aux.begin(),aux.end());
							double eMax = E+6.0*resolution;
						//2.2.2 Define integrand
							std::function<double(double)> integrand = [E,X,A,mDM,sigma,resolution,ldm,&eff_O,&eff_Ca,&eff_W] (double ER)
							{
								double sum = eff_O(E)*dRdER(ER,X[0],A[0],mDM,sigma,ldm)+eff_Ca(E)*dRdER(ER,X[1],A[1],mDM,sigma,ldm)+eff_W(E)*dRdER(ER,X[2],A[2],mDM,sigma,ldm);
								return Gauss_Distribution(E,ER,resolution)*sum;
							};

						//2.2.3 Integrate and append to list
							// double epsilon = Find_Epsilon(integrand,ERmin,ERmax,1e-5);
							double epsilon = 1e-3*(eMax-eMin)*integrand(eMin);
							double  dRdE = AdaptiveSimpsons(integrand, eMin,eMax,epsilon);
							table.push_back(std::vector<double> {E,dRdE});


					}

				//2.3 If we lowered Emax we have to put a final value, to not mess with the numerical integration later on.
				if(Emax<40*keV)
				{
					table.push_back(std::vector<double> {Emax+resolution,0.0});
					table.push_back(std::vector<double> {40*keV+6.0*resolution,0.0});
				} 
			//3. Interpolate and return spectrum
				Interpolation spectrum(table);
				return spectrum;
		}
		else if(experiment=="CRESST-surface")
		{
			//Detector Input
				double A[2] = {16.,27.};
				double X[2] = {0.47,0.53};				
				double Ethr = 19.7*eV;
				double resolution = 3.74*eV;
				double Emax = 600*eV;
			//1. Import efficiencies
				//not necessary here
			//2. Tabulate dR/dE
				//2.1 Find the maximum E for which we have contributions.
					std::vector<double> ERmax_aux = {ERMax((vesc+vEarth),mDM,A[0]),ERMax((vesc+vEarth),mDM,A[1])};
					double ERmax = *std::max_element(ERmax_aux.begin(),ERmax_aux.end());
					if((ERmax+6.0*resolution)<Emax) Emax =  (ERmax+6.0*resolution);
				//2.2 Tabulate dRdE
					std::vector<std::vector<double>> table;
					int steps = 200;
					double dE = (Emax - Ethr)/(steps-1.0);
					for(int i =0; i<steps; i++)
					{
						double E = Ethr+ i*dE;
						//2.2.1 Find minimum and maximum ER contributing to dR/dE(E):
							std::vector<double> aux = {E-6.0*resolution,Ethr-3.0*resolution};
							double eMin = *std::max_element(aux.begin(),aux.end());
							double eMax = E+6.0*resolution;
						//2.2.2 Define integrand
							std::function<double(double)> integrand = [E,X,A,mDM,sigma,resolution,ldm] (double ER)
							{
								double sum = dRdER(ER,X[0],A[0],mDM,sigma,ldm)+dRdER(ER,X[1],A[1],mDM,sigma,ldm);
								return Gauss_Distribution(E,ER,resolution)*sum;
							};
						//2.2.3 Integrate and append to list
							double epsilon = 1e-3*(eMax-eMin)*integrand(eMin);
							double  dRdE = AdaptiveSimpsons(integrand, eMin,eMax,epsilon);
							table.push_back(std::vector<double> {E,dRdE});

					}
				//2.3 If we lowered Emax we have to put a final value, to not mess with the numerical integration later on.
				if(Emax<600*eV)
				{
					table.push_back(std::vector<double> {Emax+resolution,0.0});
					table.push_back(std::vector<double> {600*eV+6.0*resolution,0.0});
				} 
			//3. Interpolate and return spectrum
				Interpolation spectrum(table);
				return spectrum;
		}
		else
		{
			//Error message
				cout <<"Error in Compute_Spectrum(): Experiment " <<experiment <<" not recognized."<<endl;
				Interpolation spectrum;
				return spectrum;
		}

	}
	Interpolation Compute_Spectrum(std::string experiment,double mDM,double sigma,bool ldm,const std::vector<DataPoint> &data,std::vector<double> attenuation)
	{
		if(experiment == "DAMIC")
		{
			//Detector input
				double efficiency = 1.0;
				double A=28.0;
				double X=1.0;
				double Ethr = 0.55*keV;
				double Emax = 7.0*keV;
			//spectrum:
				Interpolation spectrum = dRdER(Ethr,Emax,X,A,mDM,sigma,ldm,attenuation,data);
				spectrum.Multiply(efficiency);
				return spectrum;
				
		}
		else if(experiment=="XENON1T")
		{
			//Detector input
				double efficiency = 0.82;
				double A=131.0;
				double X=1.0;
				double Ethr = 5*keV; 
				double Emax = 40.0*keV;
			//spectrum:
				Interpolation spectrum = dRdER(Ethr,Emax,X,A,mDM,sigma,ldm,attenuation,data);
				spectrum.Multiply(efficiency);
				return spectrum;
		}
		else if(experiment=="CRESST-II")
		{
			//0. Detector Parameter
				double A[3] = {16.,40.,184.};
				double X[3] = {0.22,0.14,0.64};				
				double Ethr = 307*eV;
				double resolution = Ethr/5.0;
				double Emax = 40*keV;
			//1. Import efficiencies
				Interpolation eff_O("../detectors/CRESST-II/Lise_eff_AR_O.dat",keV);
				Interpolation eff_Ca("../detectors/CRESST-II/Lise_eff_AR_Ca.dat",keV);
				Interpolation eff_W("../detectors/CRESST-II/Lise_eff_AR_W.dat",keV);
			//2. Tabulate and interpolate dR/dE
					
				//2.1 Interpolate the 3 recoil spectra dR/dER
					double ERmin = Ethr-3.0*resolution;
					double ERmax = Emax + 6.0*resolution;
					Interpolation dRdER_O = dRdER( ERmin,ERmax, X[0],A[0], mDM,sigma,ldm,attenuation,data);
					Interpolation dRdER_Ca = dRdER( ERmin,ERmax, X[1],A[1], mDM,sigma,ldm,attenuation,data);
					Interpolation dRdER_W = dRdER( ERmin,ERmax, X[2],A[2], mDM,sigma,ldm,attenuation,data);

				//2.2 Find the minimum and maximum ER for which we have contributions.
					std::vector<double> ERmax_aux = {ERMax((vesc+vEarth),mDM,A[0]),ERMax((vesc+vEarth),mDM,A[1]),ERMax((vesc+vEarth),mDM,A[2])};
					ERmax = *std::max_element(ERmax_aux.begin(),ERmax_aux.end());
					if((ERmax+6.0*resolution)<Emax) Emax =  (ERmax+6.0*resolution);
				//2.3 Tabulate
					std::vector<std::vector<double>> table;
					int steps = 200;
					double dE = (Emax - Ethr)/(steps-1.0);
					for(int i =0; i<steps; i++)
					{
						double E = Ethr+ i*dE;
						//2.3.1 Find minimum and maximum ER contributing to dR/dE(E):
							std::vector<double> aux = {E-6.0*resolution,Ethr-3.0*resolution};
							double eMin = *std::max_element(aux.begin(),aux.end());
							double eMax = E+6.0*resolution;
						//2.3.2 Define integrand
							std::function<double(double)> integrand = [E,resolution,&dRdER_O,&dRdER_Ca,&dRdER_W,&eff_O,&eff_Ca,&eff_W] (double ER)
							{
								double sum = eff_O(E)*dRdER_O(ER)+eff_Ca(E)*dRdER_Ca(ER)+eff_W(E)*dRdER_W(ER);
								return Gauss_Distribution(E,ER,resolution)*sum;
							};

						//2.3.3 Integrate and append to list							
							double epsilon = 1e-3*(eMax-eMin)*(dRdER_O(eMin)+dRdER_Ca(eMin)+dRdER_W(eMin));
							double dRdE=0.0;
							if(epsilon>0) //epsilon can become negative for reasons of double precision
							{
								dRdE = AdaptiveSimpsons(integrand, eMin,eMax,epsilon);
							}
							table.push_back(std::vector<double> {E,dRdE});

					}
				//2.4 If we lowered Emax we have to put a final value, to not mess with the numerical integration later on.
					
					if(Emax<40*keV)
					{
						table.push_back(std::vector<double> {Emax+resolution,0.0});
						table.push_back(std::vector<double> {40*keV+6.0*resolution,0.0});
					} 
			//3. Interpolate and return spectrum
				Interpolation spectrum(table);
				return spectrum;
		}
		else if(experiment=="CRESST-surface")
		{
			//0. Detector Parameter
				double A[2] = {16.,27.};
				double X[2] = {0.47,0.53};				
				double Ethr = 19.7*eV;
				double resolution = 3.74*eV;
				double Emax = 600*eV;
			//1. Import efficiencies
				//not necessary here
			//2. Tabulate and interpolate dR/dE
					
				//2.1 Interpolate the 3 recoil spectra dR/dER
					double ERmin = Ethr-3.0*resolution;
					double ERmax = Emax + 6.0*resolution;
					Interpolation dRdER_O = dRdER( ERmin,ERmax, X[0],A[0], mDM,sigma,ldm,attenuation,data);
					Interpolation dRdER_Al = dRdER( ERmin,ERmax, X[1],A[1], mDM,sigma,ldm,attenuation,data);

				//2.2 Find the minimum and maximum ER for which we have contributions.
					std::vector<double> ERmax_aux = {ERMax((vesc+vEarth),mDM,A[0]),ERMax((vesc+vEarth),mDM,A[1])};
					ERmax = *std::max_element(ERmax_aux.begin(),ERmax_aux.end());
					if((ERmax+6.0*resolution)<Emax) Emax =  (ERmax+6.0*resolution);
				//2.3 Tabulate
					std::vector<std::vector<double>> table;
					int steps = 200;
					double dE = (Emax - Ethr)/(steps-1.0);
					for(int i =0; i<steps; i++)
					{
						double E = Ethr+ i*dE;
						//2.3.1 Find minimum and maximum ER contributing to dR/dE(E):
							std::vector<double> aux = {E-6.0*resolution,Ethr-3.0*resolution};
							double eMin = *std::max_element(aux.begin(),aux.end());
							double eMax = E+6.0*resolution;
						//2.3.2 Define integrand
							std::function<double(double)> integrand = [E,resolution,&dRdER_O,&dRdER_Al] (double ER)
							{
								double sum = dRdER_O(ER)+dRdER_Al(ER);
								return Gauss_Distribution(E,ER,resolution)*sum;
							};

						//2.3.3 Integrate and append to list
							double epsilon = 1e-3*(eMax-eMin)*(dRdER_O(eMin)+dRdER_Al(eMin));
							double dRdE=0.0;
							if(epsilon>0) //epsilon can become negative for reasons of double precision
							{
								dRdE = AdaptiveSimpsons(integrand, eMin,eMax,epsilon);
							}
							table.push_back(std::vector<double> {E,dRdE});

					}
				//2.4 If we lowered Emax we have to put a final value, to not mess with the numerical integration later on.
					
					if(Emax<600*eV)
					{
						table.push_back(std::vector<double> {Emax+resolution,0.0});
						table.push_back(std::vector<double> {600*eV+6.0*resolution,0.0});
					} 
			//3. Interpolate and return spectrum
				Interpolation spectrum(table);
				return spectrum;
		}
		else
		{
			//Error message
				cout <<"Error in Compute_Spectrum(): Experiment " <<experiment <<" not recognized."<<endl;
				Interpolation spectrum;
				return spectrum;
		}
		
	}

//Attenuation
	std::vector<double> Compute_Attenuation(unsigned long long int particles, const std::vector<DataPoint> &data,double vCutoff)
	{
		//1. We only simulate particles with v>vMin, hence we have to find out how many particles we neglected.
			double Simulated_Fraction=1.0;
			//Define integrand and epsilon
				auto integrand = std::bind(SpeedDistribution,_1,vEarth);
				double epsilon = 0.00001;
			//Integrate
				Simulated_Fraction = AdaptiveSimpsons(integrand,vCutoff,vEarth+vesc,epsilon);
			//True number of particles
				double Ntrue = 1.0*particles/Simulated_Fraction;
		//2. Attenuation is the sum of the weights
			std::vector<double> att;
			double Weight_Sum = 0.0;
			double Weight2_Sum = 0.0;
			for(unsigned int i=0;i<data.size();i++)
			{
				Weight_Sum+=data[i].weight;
				Weight2_Sum+=data[i].weight*data[i].weight;
			}
			//Compute attenuation
				att.push_back(Weight_Sum/Ntrue);
				att.push_back(sqrt(Weight2_Sum)/Ntrue);
		return att;
	}

//Likelihood
	double Likelihood(std::string experiment,std::function<double(double)> spectrum, double &N)
	{
		double llh=-1.0;
		//DAMIC
		if(experiment=="DAMIC")
		{
			double exposure = 0.107*kg*day;
			double Ethreshold = 0.55*keV;
			double Emax = 7*keV;
			int nEvents = 106;
			N = N_Signals(spectrum,Ethreshold,Emax,exposure);
			llh = Poisson_Likelihood(N,nEvents);
		}
		else if(experiment=="XENON1T")
		{
			double exposure = 34.2*day*1042*kg;
			double Ethreshold = 5*keV;
			double Emax = 40*keV;
			N = N_Signals(spectrum,Ethreshold,Emax,exposure);
			llh = Poisson_Likelihood(N,0);
		}
		else if(experiment=="CRESST-II")
		{
			//0. Detector parameters
				double Ethr = 307*eV;
				double Emax = 40*keV;
				double exposure = 52.15*kg*day;

			//1. Import and sort events
				std::vector<double> events = Read_List("../detectors/CRESST-II/Lise_AR.dat",keV);
				events.push_back(Ethr);
				events.push_back(Emax);
				std::sort(events.begin(),events.end());
			//2. Integrate spectrum inbetween the events
				std::vector<double> x;
				for(unsigned int i = 0;i < (events.size()-1); i++)
				{
					//2.1 The gap
					double E1 = events[i];
					double E2 = events[i+1];
					//2.2 Integrate the gap
						double epsilon = Find_Epsilon(spectrum,E1,E2,1e-2);
						double xGap = exposure*AdaptiveSimpsons(spectrum, E1,E2,epsilon);
					//2.3 Append the x value to the list
					x.push_back(xGap);
				}
			//3. Find xMax
				double xMax = *std::max_element(x.begin(),x.end());
			//4. Return likelihood.
				N = std::accumulate(x.begin(),x.end(),0.0);
				llh = 1.0 - C0(xMax,N);
		}
		else if(experiment=="CRESST-surface")
		{
			//0. Detector parameters
				double Ethr = 19.7*eV;
				double Emax = 600*eV;
				double exposure = 0.046*gram*day;

			//1. Import and sort events
				std::vector<double> events = Read_List("../detectors/CRESST-surface/data.txt",keV);
				events.push_back(Ethr);
				events.push_back(Emax);
				std::sort(events.begin(),events.end());
			//2. Integrate spectrum inbetween the events
				std::vector<double> x;
				for(unsigned int i = 0;i < (events.size()-1); i++)
				{
					//2.1 The gap
					double E1 = events[i];
					double E2 = events[i+1];
					//2.2 Integrate the gap
						double epsilon = Find_Epsilon(spectrum,E1,E2,1e-2);
						double xGap = exposure*AdaptiveSimpsons(spectrum, E1,E2,epsilon);
					//2.3 Append the x value to the list
					x.push_back(xGap);
				}
			//3. Find xMax
				double xMax = *std::max_element(x.begin(),x.end());
			//4. Return likelihood.
				N = std::accumulate(x.begin(),x.end(),0.0);
				llh = 1.0 - C0(xMax,N);
		}
		else
		{
			cout <<"Error in Likelihood(): Experiment " <<experiment <<" not recognized."<<endl;
		}
		return llh;
	}
	std::vector<double> Likelihood(std::string experiment,std::function<double(double)> spectrum,std::vector<double> attenuation,std::vector<double> &NMC)
	{
		std::vector<double> llh;
		//DAMIC
		if(experiment=="DAMIC")
		{
			double exposure = 0.107*kg*day;
			double Ethreshold = 0.55*keV;
			double Emax = 7*keV;
			double nEvents=106;
			double n = N_Signals(spectrum,Ethreshold,Emax,exposure);
			NMC.push_back(n);
			NMC.push_back(attenuation[1]/attenuation[0]*n);
			llh.push_back( Poisson_Likelihood(n,nEvents) );
			if(llh[0]==0.0) llh.push_back(0.0);
			else llh.push_back( NMC[1]*exp(-NMC[0]) * pow(NMC[0],nEvents)/Factorial(nEvents) );
		}
		else if(experiment=="XENON1T")
		{
			double exposure = 34.2*day*1042*kg;
			double Ethreshold = 5*keV;
			double Emax = 40*keV;
			double n = N_Signals(spectrum,Ethreshold,Emax,exposure);
			NMC.push_back(n);
			NMC.push_back(attenuation[1]/attenuation[0]*n);
			llh.push_back( Poisson_Likelihood(n,0) );
			llh.push_back( NMC[1]*exp(-NMC[0]) );
		}
		else if(experiment=="CRESST-II")
		{
			//0. Detector parameters
				double Ethr = 307*eV;
				double Emax = 40*keV;
				double exposure = 52.15*kg*day;

			//1. Import and sort events
				std::vector<double> events = Read_List("../detectors/CRESST-II/Lise_AR.dat",keV);
				events.push_back(Ethr);
				events.push_back(Emax);
				std::sort(events.begin(),events.end());
			//2. Integrate spectrum inbetween the events
				std::vector<double> x;
				for(unsigned int i = 0;i < (events.size()-1); i++)
				{
					//2.1 The gap
					double E1 = events[i];
					double E2 = events[i+1];
					//2.2 Integrate the gap
						double epsilon = Find_Epsilon(spectrum,E1,E2,1e-2);
						double xGap = exposure*AdaptiveSimpsons(spectrum, E1,E2,epsilon);
					//2.3 Append the x value to the list
					x.push_back(xGap);
				}
			//3. Find xMax
				double xMax = *std::max_element(x.begin(),x.end());
			//4. Return likelihood.
				double N = std::accumulate(x.begin(),x.end(),0.0);
				double error = attenuation[1]/attenuation[0]*N;
				NMC.push_back( N );
				NMC.push_back( error );
				llh.push_back(1.0 - C0(xMax,N));
				llh.push_back(C0_error(xMax,N,error) );
		}
		else if(experiment=="CRESST-surface")
		{
			//0. Detector parameters
				double Ethr = 19.7*eV;
				double Emax = 600*eV;
				double exposure = 0.046*gram*day;

			//1. Import and sort events
				std::vector<double> events = Read_List("../detectors/CRESST-surface/data.txt",keV);
				events.push_back(Ethr);
				events.push_back(Emax);
				std::sort(events.begin(),events.end());
			//2. Integrate spectrum inbetween the events
				std::vector<double> x;
				for(unsigned int i = 0;i < (events.size()-1); i++)
				{
					//2.1 The gap
					double E1 = events[i];
					double E2 = events[i+1];
					//2.2 Integrate the gap
						double epsilon = Find_Epsilon(spectrum,E1,E2,1e-2);
						double xGap = exposure*AdaptiveSimpsons(spectrum, E1,E2,epsilon);
					//2.3 Append the x value to the list
					x.push_back(xGap);
				}
			//3. Find xMax
				double xMax = *std::max_element(x.begin(),x.end());
			//4. Return likelihood.
				double N = std::accumulate(x.begin(),x.end(),0.0);
				double error = attenuation[1]/attenuation[0]*N;
				NMC.push_back( N );
				NMC.push_back( error );
				llh.push_back(1.0 - C0(xMax,N));
				llh.push_back(C0_error(xMax,N,error) );
		}
		else
		{
			cout <<"Error in Likelihood(): Experiment " <<experiment <<" not recognized."<<endl;
		}
		return llh;
	}

//Total number of events
	double N_Signals(std::function<double(double)> spectrum,double E1,double E2,double exposure)
	{
		double epsilon = Find_Epsilon(spectrum,E1,E2,1e-4);
		return exposure*AdaptiveSimpsons(spectrum, E1,E2,epsilon);
	}

//Standard Upper bound
	 double Upper_Bound(std::string experiment,double mDM,bool ldm,double CL)
	 {
	 	double sigma=-1.0;

		//DAMIC
		if(experiment=="DAMIC")
		{
			double N_Limit = 120.452;
			double exposure = 0.107*kg*day;
			double Ethreshold = 0.55*keV;
			double Emax = 7*keV;
			std::function <double(double)> spectrum = Compute_Spectrum(experiment,mDM,1.0*pb,ldm);
			double N = N_Signals(spectrum,Ethreshold,Emax,exposure);
			sigma = N_Limit/N*pb;
		}
		else if(experiment=="XENON1T")
		{
			double N_Limit = -log(1.0-CL);
			double exposure = 34.2*day*1042*kg;
			double Ethreshold = 5*keV;
			double Emax = 40*keV;
			std::function <double(double)> spectrum = Compute_Spectrum(experiment,mDM,1.0*pb,ldm);
			double N = N_Signals(spectrum,Ethreshold,Emax,exposure);
			sigma = N_Limit/N*pb;
		}
		else if(experiment == "CRESST-II")
		{
			//Find critical cross-section using bisection
				double n;
				double s1=1e-45*cm*cm;
				double s2=1e-30*cm*cm;
				double s3 = pow(10.0, log10(s1*s2)/2.0);
				double e = 0.01;
				double llh=0.0;
				while(fabs(llh-(1.0-CL))>e)
				{
					Interpolation spectrum = Compute_Spectrum(experiment,mDM,s3,ldm);
					llh = Likelihood(experiment,spectrum,n);
					if(llh>(1.0-CL))
					{
						s1 = s3;
					}
					else
					{
						s2 = s3;
					}
					s3 = pow(10.0, log10(s1*s2)/2.0);
				}
				sigma=s3;
		}
		else if(experiment == "CRESST-surface")
		{
			//Find critical cross-section using bisection
				double n;
				double s1=1e-36*cm*cm;
				double s2=1e-28*cm*cm;
				double s3 = pow(10.0, log10(s1*s2)/2.0);
				double e = 0.001;
				double llh=0.0;
				while(fabs(llh-(1.0-CL))>e)
				{
					Interpolation spectrum = Compute_Spectrum(experiment,mDM,s3,ldm);
					llh = Likelihood(experiment,spectrum,n);
					if(llh>(1.0-CL))
					{
						s1 = s3;
					}
					else
					{
						s2 = s3;
					}
					s3 = pow(10.0, log10(s1*s2)/2.0);
				}
				sigma=s3;
		}
		else
		{
			cout <<"Error in Likelihood(): Experiment " <<experiment <<" not recognized."<<endl;
		}
		return sigma;
	 }

