#include "DD_General.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <algorithm>
#include <chrono>
#include <numeric> //for std::accumulate

#include "Physics_Functions.hpp"
#include "Input_Parameters.hpp"
#include "DD_Nucleus.hpp"
#include "DD_Electron.hpp"

using namespace std::placeholders;

//1. Attenuation
	std::vector<double> Compute_Attenuation(unsigned long long int particles, const std::vector<DataPoint> &data,double vCutoff)
	{
		//1. We only simulate particles with v>vMin, hence we have to find out how many particles we neglected.
			//Define integrand and epsilon
				auto integrand = std::bind(SpeedDistribution,_1,vEarth);
				double epsilon = 0.00001;
			//Integrate
				double Simulated_Fraction = Integrate(integrand,vCutoff,vEarth+vesc,epsilon);
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

//2. Eta-Function: Integral of f(v)/v with lower integration limit v_Min
	//Analytic
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
	
	//MC based
	Interpolation EtaFunction(const std::vector<double>& attenuation,const std::vector<DataPoint> &speeddata,double vCutoff,double vE)
	{
		Interpolation kde = Perform_KDE(speeddata,vCutoff,(vesc+vEarth));
		kde.Multiply(attenuation[0]);
		int points =200;
		double dv=(vesc+vE-vCutoff)/(points-1.0);
		std::vector<std::vector<double>> interpol_list;
		for(int i=0;i<199;i++)
		{
			double vMin = vCutoff+i*dv;
			std::function<double(double)> integrand = [&kde](double v)
			{
				return kde(v)/v;
			};
			double epsilon=1e-2*(vesc+vE-vMin)*(integrand(vMin)+integrand((vMin+vEarth+vesc)/2.0))/2.0;
			double value=Integrate(integrand,vMin,vesc+vE,epsilon);
			interpol_list.push_back(std::vector<double> {vMin,value});
		}
		interpol_list.push_back(std::vector<double> {(vesc+vE),0.0});
		interpol_list.push_back(std::vector<double> {1.0,0.0});
		interpol_list.push_back(std::vector<double> {10.0,0.0});

		Interpolation eta(interpol_list);
		return eta;
	}

//3. Compute general recoil spectrum
	Interpolation Compute_Spectrum(const DM_Particle &DM,const std::string& experiment)
	{
		if(experiment == "DAMIC")
		{
			//1. Detector input
				// double efficiency = 0.75;
				double efficiency = 1.0;
				int Z = 14;
				double A=28.0;
				double X=1.0;
				double Ethr = 0.55*keV;
				double Emax = 7.0*keV;
			//2. Tabulate and interpolate
				//2.1 Find maximum ER, i.e. the domain of the spectrum
					double ERmax = ERMax((vesc+vEarth),DM.mass,A);
					ERmax = std::min(ERmax,Emax);
				//2.2 Tabulate the spectrum
					int points = 250;
					double dER = (ERmax-Ethr)/(points-1.0);
					std::vector<std::vector<double>> interpol_list;
					for(int i = 0;i<points;i++)
					{
						double ER = Ethr + i*dER;
						double dR = efficiency*dRdER(ER,DM,X,Z,A);
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
				int Z = 54;
				double A=131.0;
				double X=1.0;
				double Ethr = 5*keV;
				double Emax = 40.0*keV;
			//2. Tabulate and interpolate
				//2.1 Find maximum ER, i.e. the domain of the spectrum
					double ERmax = ERMax((vesc+vEarth),DM.mass,A);
					ERmax = std::min(ERmax,Emax);
				//2.2 Tabulate the spectrum
					int points = 250;
					double dER = (ERmax-Ethr)/(points-1.0);
					std::vector<std::vector<double>> interpol_list;
					for(int i = 0;i<points;i++)
					{
						double ER = Ethr + i*dER;
						double dR = efficiency*dRdER(ER,DM,X,Z,A);
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
				int Z[3]	= {8,20,74};
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
					std::vector<double> ERmax_aux = {ERMax((vesc+vEarth),DM.mass,A[0]),ERMax((vesc+vEarth),DM.mass,A[1]),ERMax((vesc+vEarth),DM.mass,A[2])};
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
							std::function<double(double)> integrand = [E,X,Z,A,DM,resolution,&eff_O,&eff_Ca,&eff_W] (double ER)
							{
								double sum = eff_O(E)*dRdER(ER,DM,X[0],Z[0],A[0])+eff_Ca(E)*dRdER(ER,DM,X[1],Z[1],A[1])+eff_W(E)*dRdER(ER,DM,X[2],Z[2],A[2]);
								return PDF_Gauss(E,ER,resolution)*sum;
							};

						//2.2.3 Integrate and append to list
							// double epsilon = Find_Epsilon(integrand,ERmin,ERmax,1e-5);
							double epsilon = 1e-3*(eMax-eMin)*integrand(eMin);
							double  dRdE = Integrate(integrand, eMin,eMax,epsilon);
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
				double Z[2] = {8,13};
				double A[2] = {16.,27.};
				double X[2] = {0.47,0.53};				
				double Ethr = 19.7*eV;
				double resolution = 3.74*eV;
				double Emax = 600*eV;
			//1. Import efficiencies
				//not necessary here
			//2. Tabulate dR/dE
				//2.1 Find the maximum E for which we have contributions.
					std::vector<double> ERmax_aux = {ERMax((vesc+vEarth),DM.mass,A[0]),ERMax((vesc+vEarth),DM.mass,A[1])};
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
							std::function<double(double)> integrand = [E,X,Z,A,DM,resolution] (double ER)
							{
								double sum = dRdER(ER,DM,X[0],Z[0],A[0])+dRdER(ER,DM,X[1],Z[1],A[1]);
								return PDF_Gauss(E,ER,resolution)*sum;
							};
						//2.2.3 Integrate and append to list
							double epsilon = 1e-3*(eMax-eMin)*integrand(eMin);
							double  dRdE = Integrate(integrand, eMin,eMax,epsilon);
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
		else if(experiment == "Semiconductor"||experiment=="SENSEI"||experiment == "SENSEI-surface"||experiment=="SuperCDMS"||experiment=="DAMIC-M")
		{
			//Tabulate and interpolate
				//1. Domain of the spectrum
					double Eemax = 50*eV;
					double Eemin = 0.1*eV;
				//2. Tabulate the spectrum
					int points = 250;
					double dEe = (Eemax-Eemin)/(points-1.0);
					std::vector<std::vector<double>> interpol_list;
					for(int i = 0;i<points;i++)
					{
						double Ee = Eemin + i*dEe;
						double dR = dRdEe(Ee,DM,DMe_target,DMe_efficiency);
						interpol_list.push_back(std::vector<double> {Ee,dR});
					}
				//3. Interpolate and return
					Interpolation spectrum(interpol_list);
					return spectrum;
		}
		else if(experiment=="XENON10e"||experiment=="XENON100e")
		{
			//Tabulate and interpolate
				//1. Tabulate the PE spectrum
					int points = 250;
					std::vector<std::vector<double>> interpol_list;
					for(int nPE = 1;nPE<=points;nPE++)
					{
						double dN = dNdnPE(nPE,DM,DMe_target,experiment,Shells);
						interpol_list.push_back(std::vector<double> {1.0*nPE,dN});
					}
				//2. Interpolate and return
					Interpolation spectrum(interpol_list);
					return spectrum;
		}
		else if(experiment=="DarkSide-50")
		{
			//1. Tabulate the ne spectrum
				std::vector<std::vector<double>> interpol_list;
				for(int ne=1;ne<=10;ne++)
				{
					double dR=0.0;
					for(unsigned int i=0;i<Shells.size();i++)
					{						
						dR+= dRdne(ne,DM,DMe_target,Shells[i]);
					} 
					interpol_list.push_back(std::vector<double> {1.0*ne,dR});
				}
			//2. Interpolate and return
				Interpolation spectrum(interpol_list);
				return spectrum;
		}
		else
		{
			//Error
				cerr <<"Error in Compute_Spectrum(): Experiment " <<experiment <<" not recognized."<<endl;
				std::exit(EXIT_FAILURE);
		}
	}
	Interpolation Compute_Spectrum(const DM_Particle &DM,const std::string& experiment,const std::vector<DataPoint> &data,const std::vector<double>& attenuation)
	{
		if(experiment == "DAMIC")
		{
			//Detector input
				double efficiency = 1.0;
				int Z = 14;
				double A=28.0;
				double X=1.0;
				double Ethr = 0.55*keV;
				double Emax = 7.0*keV;
			//spectrum:
				Interpolation spectrum = dRdER(DM,Ethr,Emax,X,Z,A,attenuation,data);
				spectrum.Multiply(efficiency);
				return spectrum;
				
		}
		else if(experiment=="XENON1T")
		{
			//Detector input
				double efficiency = 0.82;
				int Z = 54;
				double A=131.0;
				double X=1.0;
				double Ethr = 5*keV; 
				double Emax = 40.0*keV;
			//spectrum:
				Interpolation spectrum = dRdER(DM,Ethr,Emax,X,Z,A,attenuation,data);
				spectrum.Multiply(efficiency);
				return spectrum;
		}
		else if(experiment=="CRESST-II")
		{

			//0. Detector Parameter
				int Z[3]	= {8,20,74};
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
					Interpolation dRdER_O = dRdER(DM, ERmin,ERmax, X[0],Z[0],A[0], attenuation,data);
					Interpolation dRdER_Ca = dRdER( DM,ERmin,ERmax, X[1],Z[1],A[1], attenuation,data);
					Interpolation dRdER_W = dRdER(DM, ERmin,ERmax, X[2],Z[2],A[2],attenuation,data);

				//2.2 Find the minimum and maximum ER for which we have contributions.
					std::vector<double> ERmax_aux = {ERMax((vesc+vEarth),DM.mass,A[0]),ERMax((vesc+vEarth),DM.mass,A[1]),ERMax((vesc+vEarth),DM.mass,A[2])};
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
								return PDF_Gauss(E,ER,resolution)*sum;
							};

						//2.3.3 Integrate and append to list							
							double epsilon = 1e-3*(eMax-eMin)*(dRdER_O(eMin)+dRdER_Ca(eMin)+dRdER_W(eMin));
							double dRdE=0.0;
							if(epsilon>0) //epsilon can become negative for reasons of double precision
							{
								dRdE = Integrate(integrand, eMin,eMax,epsilon);
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
				int Z[2]	= {8,13};
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
					Interpolation dRdER_O = 	dRdER(DM,ERmin,ERmax,X[0],Z[0],A[0],attenuation,data);
					Interpolation dRdER_Al = 	dRdER(DM,ERmin,ERmax,X[1],Z[1],A[1],attenuation,data);

				//2.2 Find the minimum and maximum ER for which we have contributions.
					std::vector<double> ERmax_aux = {ERMax((vesc+vEarth),DM.mass,A[0]),ERMax((vesc+vEarth),DM.mass,A[1])};
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
								return PDF_Gauss(E,ER,resolution)*sum;
							};

						//2.3.3 Integrate and append to list
							double epsilon = 1e-3*(eMax-eMin)*(dRdER_O(eMin)+dRdER_Al(eMin));
							double dRdE=0.0;
							if(epsilon>0) //epsilon can become negative for reasons of double precision
							{
								dRdE = Integrate(integrand, eMin,eMax,epsilon);
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
		else if(experiment=="Semiconductor"||experiment=="SENSEI"||experiment=="SENSEI-surface"||experiment=="XENON10e"||experiment=="XENON100e"||experiment=="DarkSide-50"||experiment=="SuperCDMS"||experiment=="DAMIC-M")
		{
			//not necessary
			return Interpolation();
		}
		else
		{
			cerr <<"Error in Compute_Spectrum(): Experiment " <<experiment <<" not recognized."<<endl;
			std::exit(EXIT_FAILURE);
		}		
	}


//4. Likelihood
	double Likelihood(const DM_Particle& DM,const std::string& experiment,const std::function<double(double)>& spectrum, double &N)
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
			llh = CDF_Poisson(N,nEvents);
		}
		else if(experiment=="XENON1T")
		{
			double exposure = 34.2*day*1042*kg;
			double Ethreshold = 5*keV;
			double Emax = 40*keV;
			N = N_Signals(spectrum,Ethreshold,Emax,exposure);
			llh = CDF_Poisson(N,0);
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
						double xGap = exposure*Integrate(spectrum, E1,E2,epsilon);
					//2.3 Append the x value to the list
					x.push_back(xGap);
				}
			//3. Find xMax
				double xMax = *std::max_element(x.begin(),x.end());
			//4. Return likelihood.
				N = std::accumulate(x.begin(),x.end(),0.0);
				llh = 1.0 - MaxGap_C0(xMax,N);
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
						double xGap = exposure*Integrate(spectrum, E1,E2,epsilon);
					//2.3 Append the x value to the list
					x.push_back(xGap);
				}
			//3. Find xMax
				double xMax = *std::max_element(x.begin(),x.end());
			//4. Return likelihood.
				N = std::accumulate(x.begin(),x.end(),0.0);
				llh = 1.0 - MaxGap_C0(xMax,N);
		}
		else if(experiment=="Semiconductor")
		{
			N = DMe_exposure*TotalRate(DM,DMe_threshold,DMe_efficiency,DMe_target);
			llh = CDF_Poisson(N,DMe_events);

		}
		else if(experiment=="SENSEI"||experiment=="SENSEI-surface"||experiment=="SuperCDMS"||experiment=="DAMIC-M")
		{
			//Efficiency and number of events
				std::vector<double> efficiency;
				std::vector<int> BinData;
				if(experiment=="SENSEI")
				{
					efficiency={1.0,0.62,0.48};
					BinData = {8516,87,0};

				}
				else if(experiment=="SENSEI-surface")
				{
					efficiency={0.668,0.41,0.32,0.27,0.24};
					BinData={140302,4676,131,1,0};
				}
				else if(experiment=="SuperCDMS")
				{
					//all points within 2sigma, taken from fig. 3 of [arXiv:1804.10697]
					efficiency={0.88,0.91,0.91,0.91,0.91,0.91};
					BinData={53000, 400, 74, 18, 7, 14};
				}
				else if(experiment=="DAMIC-M")
				{
					efficiency={1.0,1.0,1.0,1.0,1.0,1.0};
					BinData={100000,0,0,0,0,0};
				}
			//Compute rates
				std::vector<double> rates;
				for(unsigned int bin=1;bin<=BinData.size()+1;bin++) rates.push_back(TotalRate(DM,bin,DMe_efficiency,DMe_target));
			//Compute number of events and likelihood per bin
				std::vector<double> events;
				std::vector<double> llhs;
				for(unsigned int bin=(DMe_threshold-1);bin<BinData.size();bin++) 
				{
					events.push_back(efficiency[bin]*DMe_exposure*(rates[bin]-rates[bin+1]));
					llhs.push_back(CDF_Poisson(events[bin],BinData[bin]));
					N+=events[bin];
				}
			//Find minimum value
				llh = *std::min_element(llhs.begin(),llhs.end());
		}
		else if(experiment=="XENON10e")
		{
			//Data is binned in 7 bins [14,41),[41,68),[68,95),[95,122),[122,149),[149,176),[176,203)
				std::vector<int> data={126, 60, 12, 3, 2, 0, 2};
				std::vector<unsigned int> bins={14,41,68,95,122,149,176,203};
				std::vector<double> llhs;
			//Compute events in each bin.
				for(unsigned int bin=0;bin<data.size();bin++)
				{
					double Nbin=0;
					for(unsigned int PE= bins[bin];PE<bins[bin+1];PE++)
					{
						Nbin+=spectrum(PE);
					}
					N+=Nbin;
					double l=CDF_Poisson(Nbin,data[bin]);
					llhs.push_back(l);

				}
			llh = *std::min_element(llhs.begin(),llhs.end());
		}
		else if(experiment=="XENON100e")
		{
			//Data is binned in 3 bins [80,90),[90,110),[110,130)
				std::vector<int> data={794, 1218, 924, 776, 669, 630, 528, 488, 433, 387};
				std::vector<unsigned int> bins={80, 90, 110,130,150,170,190,210,230,250,270};
				std::vector<double> llhs;		
			//Compute events in each bin.
				//We only use the first 3 bins.
				for(unsigned int bin=0;bin<3;bin++)
				{
					double Nbin=0;
					for(unsigned int PE= bins[bin];PE<bins[bin+1];PE++)
					{
						Nbin+=dNdnPE(PE,DM,DMe_target,experiment,Shells);
					}
					N+=Nbin;
					double l=CDF_Poisson(Nbin,data[bin]);
					llhs.push_back(l);
				}
			llh = *std::min_element(llhs.begin(),llhs.end());
		}
		else if(experiment=="DarkSide-50")
		{
			std::vector<int>data={118643, 219893, 6131, 673, 252, 227, 198, 199, 189, 247, 230, 261,249, 329, 336, 349, 351, 352, 384, 411, 405, 461, 460, 436, 500, 546, 538, 536, 556, 583, 573, 630, 603, 635, 639, 682, 736, 755, 804, 811, 809, 882, 934, 935, 871, 965, 946, 1072, 997, 1060};
			std::vector<double> llhs;
			//Threshold and exposure
				double exposure = 6786.0*kg*day;
				int nThreshold=3;
			//Compute likelihood for each bin.
				for(unsigned int ne=nThreshold;ne<data.size();ne++)
				{
					double Ntot=0.0;
					//Add all shells
					for(unsigned int shell=0;shell<Shells.size();shell++)
					{					
						Ntot+=exposure/2.0 * (dRdne(ne,DM,DMe_target,Shells[shell])+dRdne(ne+1,DM,DMe_target,Shells[shell]));
					
					}
					N+=Ntot;
					double l=CDF_Poisson(Ntot,data[ne-1]);
					llhs.push_back(l);
				}
			llh = *std::min_element(llhs.begin(),llhs.end());

		}
		else
		{
			std::cerr <<"Error in Likelihood(): Experiment " <<experiment <<" not recognized."<<endl;
			std::exit(EXIT_FAILURE);
		}
		return llh;
	}

	std::vector<double> Likelihood(const DM_Particle& DM,const std::string& experiment,const std::function<double(double)>& spectrum,const std::vector<double>& attenuation,const std::vector<DataPoint> &data,double vCutoff,std::vector<double> &NMC)
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
			llh.push_back( CDF_Poisson(n,nEvents) );
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
			llh.push_back( CDF_Poisson(n,0) );
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
						double xGap = exposure*Integrate(spectrum, E1,E2,epsilon);
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
				llh.push_back(1.0 - MaxGap_C0(xMax,N));
				llh.push_back(MaxGap_C0_error(xMax,N,error) );
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
						double xGap = exposure*Integrate(spectrum, E1,E2,epsilon);
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
				llh.push_back(1.0 - MaxGap_C0(xMax,N));
				llh.push_back(MaxGap_C0_error(xMax,N,error) );
		}
		else if(experiment=="Semiconductor")
		{
			double N = DMe_exposure*TotalRate(DM,DMe_threshold,DMe_efficiency,DMe_target,attenuation,data,vCutoff);
			double error = attenuation[1]/attenuation[0]*N;
			NMC.push_back( N );
			NMC.push_back( error );
			llh.push_back(CDF_Poisson(N,DMe_events));
			if(llh[0]==0.0) llh.push_back(0.0);
			else llh.push_back( NMC[1]*PDF_Poisson(NMC[0],DMe_events) );

		}
		else if(experiment=="SENSEI" || experiment=="SENSEI-surface"||experiment=="SuperCDMS"||experiment=="DAMIC-M")
		{
			//Efficiency and number of events per bin
				std::vector<double> efficiency;
				std::vector<int> BinData;
				if(experiment=="SENSEI")
				{
					efficiency={1.0,0.62,0.48};
					BinData = {8516,87,0};
				}
				else if(experiment=="SENSEI-surface")
				{
					efficiency={0.668,0.41,0.32,0.27,0.24};
					BinData={140302,4676,131,1,0};
				}
				else if(experiment=="SuperCDMS")
				{
					//all points within 2sigma, taken from fig 3 of 1804.10697
					efficiency={0.88,0.91,0.91,0.91,0.91,0.91};
					BinData={53000, 400, 74, 18, 7, 14};
				}
				else if(experiment=="DAMIC-M")
				{
					efficiency={1.0,1.0,1.0,1.0,1.0,1.0};
					BinData={100000,0,0,0,0,0};
				}
			//Compute spectrum
				Interpolation spectrum = dRdEe(DM,DMe_target,DMe_efficiency,attenuation,data,vCutoff);
			//Compute rates
				std::vector<double> rates;
				double Egap = 1.11*eV;
				double eps = 3.6*eV;
				for(unsigned int bin=DMe_threshold;bin<BinData.size()+1;bin++) 
				{
					double E1= eps*(bin-1.0)+Egap;
					double E2= eps*(bin)+Egap;
					double epsilon = 1.0e-5*(E2-E1)*spectrum(E1);
					double rate = Integrate(spectrum,E1,E2,epsilon);
					rates.push_back(rate);
				}
			//Compute number of events and likelihood per bin, and total number of events
				std::vector<double> events;
				std::vector<double> llhs;
				double N = 0;
				for(unsigned int bin=(DMe_threshold-1);bin<BinData.size();bin++) 
				{
					events.push_back(efficiency[bin]*DMe_exposure*rates[bin]);
					llhs.push_back(CDF_Poisson(events[bin],BinData[bin]));
					N+=events[bin];
				}
			//Total number of events
				double error = attenuation[1]/attenuation[0]*N;
				NMC.push_back( N );
				NMC.push_back( error );
			//Find minimum value
				llh.push_back( *std::min_element(llhs.begin(),llhs.end()) );
				llh.push_back(0.0);

		}
		else if(experiment=="XENON10e")
		{
			std::vector<double> dRdn=dRdne(DM,DMe_target,Shells,attenuation,data,vCutoff);
			double muPE = 27.0;
			double sigPE=6.7;
			std::vector<double> dRdnpe=dRdnPE(muPE,sigPE,dRdn);
			std::vector<double> dNdnpe=dNdnPE(experiment,dRdnpe);
			double N=0.0;
			//Data is binned in 7 bins [14,41),[41,68),[68,95),[95,122),[122,149),[149,176),[176,203)
				std::vector<int> events={126, 60, 12, 3, 2, 0, 2};
				std::vector<unsigned int> bins={14,41,68,95,122,149,176,203};
				std::vector<double> llhs;

			//Compute events in each bin.
				for(unsigned int bin=0;bin<events.size();bin++)
				{
					double Nbin=0;
					for(unsigned int PE= bins[bin];PE<bins[bin+1];PE++)
					{
						Nbin+=dNdnpe[PE-1];						
					}
					N+=Nbin;
					double l=CDF_Poisson(Nbin,events[bin]);
					llhs.push_back(l);

				}
			double error = attenuation[1]/attenuation[0]*N;
			double llh_min = *std::min_element(llhs.begin(),llhs.end());

			NMC.push_back( N );
			NMC.push_back( error );
			llh.push_back(llh_min);
			if(llh[0]==0.0) llh.push_back(0.0);
			else llh.push_back( NMC[1]*exp(-NMC[0]) );

		}
		else if(experiment=="XENON100e")
		{
			std::vector<double> dRdn=dRdne(DM,DMe_target,Shells,attenuation,data,vCutoff);
			double muPE = 19.7;
			double sigPE=6.2;
			std::vector<double> dRdnpe=dRdnPE(muPE,sigPE,dRdn);
			std::vector<double> dNdnpe=dNdnPE(experiment,dRdnpe);
			double N=0.0;
			//Data is binned in 3 bins [80,90),[90,110),[110,130)
				std::vector<int> events={794, 1218, 924, 776, 669, 630, 528, 488, 433, 387};
				std::vector<unsigned int> bins={80, 90, 110,130,150,170,190,210,230,250,270};
				std::vector<double> llhs;

			//Compute events in each bin for the first 3 bins.
				for(unsigned int bin=0;bin<3;bin++)
				{
					double Nbin=0;
					for(unsigned int PE= bins[bin];PE<bins[bin+1];PE++)
					{
						Nbin+=dNdnpe[PE-1];
					}
					N+=Nbin;
					double l=CDF_Poisson(Nbin,events[bin]);
					llhs.push_back(l);

				}
			double error = attenuation[1]/attenuation[0]*N;
			double llh_min = *std::min_element(llhs.begin(),llhs.end());

			NMC.push_back( N );
			NMC.push_back( error );
			llh.push_back(llh_min);
			if(llh[0]==0.0) llh.push_back(0.0);
			else llh.push_back( NMC[1]*exp(-NMC[0]) );
		}
		else if(experiment=="DarkSide-50")
		{
			//Experiment's parameter
				std::vector<int>events={118643, 219893, 6131, 673, 252, 227, 198, 199, 189, 247, 230, 261,249, 329, 336, 349, 351, 352, 384, 411, 405, 461, 460, 436, 500, 546, 538, 536, 556, 583, 573, 630, 603, 635, 639, 682, 736, 755, 804, 811, 809, 882, 934, 935, 871, 965, 946, 1072, 997, 1060};
				double exposure = 6786.0*kg*day;
				int nThreshold=3;
			//Spectrum
				std::vector<double> dRdn=dRdne(DM,DMe_target,Shells,attenuation,data,vCutoff);
			//Compute likelihood for each bin.
				std::vector<double> llhs;
				double Ntot=0.0;
				for(unsigned int ne=nThreshold;ne<events.size();ne++)
				{
					// double N=exposure * dRdn[ne-1];
					double N=exposure/2.0 * (dRdn[ne-1]+dRdn[ne]);
					Ntot+=N;
					double l=CDF_Poisson(N,events[ne-1]);
					llhs.push_back(l);
				}
			//Find the minimum likelihood.	
				double llh_min = *std::min_element(llhs.begin(),llhs.end());
				double error = attenuation[1]/attenuation[0]*Ntot;

			NMC.push_back( Ntot );
			NMC.push_back( error );
			llh.push_back(llh_min);
			if(llh[0]==0.0) llh.push_back(0.0);
			else llh.push_back( NMC[1]*exp(-NMC[0]) );
		}
		else
		{
			cerr <<"Error in Likelihood(): Experiment " <<experiment <<" not recognized."<<endl;
			std::exit(EXIT_FAILURE);
		}
		return llh;
	}

//5. Total number of events
	double N_Signals(const std::function<double(double)>& spectrum,double E1,double E2,double exposure)
	{
		double epsilon = Find_Epsilon(spectrum,E1,E2,1e-4);
		return exposure*Integrate(spectrum, E1,E2,epsilon);
	}

//6. Standard constraint
	 double Upper_Bound(const DM_Particle& DM,const std::string& experiment,double CL)
	{
	 	double sigma=-1.0;

		//DAMIC
		if(experiment=="DAMIC")
		{
			double N_Limit = Inv_CDF_Poisson(106,1.0-CL);
			double exposure = 0.107*kg*day;
			double Ethreshold = 0.55*keV;
			double Emax = 7*keV;
			DM_Particle dm=DM;
			dm.Set_Sigma_n(pb);
			std::function <double(double)> spectrum = Compute_Spectrum(dm,experiment);
			double N = N_Signals(spectrum,Ethreshold,Emax,exposure);
			sigma = N_Limit/N*pb;
		}
		else if(experiment=="XENON1T")
		{
			double N_Limit = -log(1.0-CL);
			double exposure = 34.2*day*1042*kg;
			double Ethreshold = 5*keV;
			double Emax = 40*keV;
			DM_Particle dm=DM;
			dm.Set_Sigma_n(pb);
			std::function <double(double)> spectrum = Compute_Spectrum(dm,experiment);
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
					DM_Particle dm=DM;
					dm.Set_Sigma_n(s3);
					Interpolation spectrum = Compute_Spectrum(dm,experiment);
					llh = Likelihood(dm,experiment,spectrum,n);
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
					DM_Particle dm=DM;
					dm.Set_Sigma_n(s3);
					Interpolation spectrum = Compute_Spectrum(dm,experiment);
					llh = Likelihood(dm,experiment,spectrum,n);
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
		else if(experiment == "Semiconductor")
		{
			double N_Limit =Inv_CDF_Poisson(DMe_events,1.0-CL);
			DM_Particle dm=DM;
			dm.Set_Sigma_e(pb);
			double N = DMe_exposure*TotalRate(dm,DMe_threshold,DMe_efficiency,DMe_target);
			sigma = N_Limit/N*pb;
		}
		else if(experiment =="SENSEI-surface"||experiment=="SENSEI"||experiment=="SuperCDMS"||experiment=="DAMIC-M")
		{
			//Efficiency and number of events per bin
				std::vector<double> efficiency;
				std::vector<int> BinData;
				if(experiment=="SENSEI")
				{
					efficiency={1.0,0.62,0.48};
					BinData = {8516,87,0};

				}
				else if(experiment=="SENSEI-surface")
				{
					efficiency={0.668,0.41,0.32,0.27,0.24};
					BinData={140302,4676,131,1,0};
				}
				else if(experiment=="SuperCDMS")
				{
					// all points within 2sigma, taken from fig 3 of 1804.10697
					efficiency={0.88,0.91,0.91,0.91,0.91,0.91};
					BinData={53000, 400, 74, 18, 7, 14};
				}
				else if(experiment=="DAMIC-M")
				{
					efficiency={1.0,1.0,1.0,1.0,1.0,1.0};
					BinData={100000,0,0,0,0,0};
				}
			//Compute rates for 1 pb
				std::vector<double> rates;
				DM_Particle dm=DM;
				dm.Set_Sigma_e(pb);
				for(unsigned int bin=1;bin<=(BinData.size()+1);bin++) 
				{					
					rates.push_back(TotalRate(dm,bin,DMe_efficiency,DMe_target));
				}
			//Compute number of events and limit per bin
				std::vector<double>sigmas;
				for(unsigned int bin=(DMe_threshold-1);bin<BinData.size();bin++)
				{
					double N_Limit=Inv_CDF_Poisson(BinData[bin],1.0-CL);
					double Nbin =efficiency[bin]*DMe_exposure*(rates[bin]-rates[bin+1]);
					double sig = N_Limit/Nbin*pb;
					sigmas.push_back(sig);
				}
			//Find minimum value
				sigma = *std::min_element(sigmas.begin(),sigmas.end());
		}
		else if(experiment == "XENON10e")
		{
			//Data is binned in 7 bins [14,41),[41,68),[68,95),[95,122),[122,149),[149,176),[176,203)
				std::vector<int> data={126, 60, 12, 3, 2, 0, 2};
				std::vector<unsigned int> bins={14,41,68,95,122,149,176,203};

				std::vector<double>sigmas;
			//Compute events in each bin with 1pb.
				DM_Particle dm=DM;
				dm.Set_Sigma_e(pb);
				for(unsigned int bin=0;bin<data.size();bin++)
				{
					double N_Limit=Inv_CDF_Poisson(data[bin],1.0-CL);
					double Nbin=0;
					for(unsigned int PE= bins[bin];PE<bins[bin+1];PE++)
					{
						Nbin+=dNdnPE(PE,dm,DMe_target,experiment,Shells);
					}
					
					sigmas.push_back(N_Limit/Nbin*pb);
				}
			sigma = *std::min_element(sigmas.begin(),sigmas.end());
		}
		else if(experiment == "XENON100e")
		{
			//Data is binned in 3 bins [80,90),[90,110),[110,130)
				std::vector<int> data={794, 1218, 924, 776, 669, 630, 528, 488, 433, 387};
				std::vector<unsigned int> bins={80, 90, 110,130,150,170,190,210,230,250,270};
				
				//Inverting the Poisson CDF is not straight forward (look at inverse regularized gamma functions at some point)
					std::vector<double> Nl;
					for(unsigned int j=0;j<data.size();j++) Nl.push_back(Inv_CDF_Poisson(data[j],1.0-CL));
				
				std::vector<double>sigmas;
			//Compute events in each bin with 1pb.
				//We only use the first 3 bins.
				DM_Particle dm=DM;
				dm.Set_Sigma_e(pb);
				for(unsigned int bin=0;bin<3;bin++)
				{
					double N_Limit=Nl[bin];
					double Nbin=0;
					for(unsigned int PE= bins[bin];PE<bins[bin+1];PE++)
					{
						Nbin+=dNdnPE(PE,dm,DMe_target,experiment,Shells);
					}
					sigmas.push_back(N_Limit/Nbin*pb);
				}
			sigma = *std::min_element(sigmas.begin(),sigmas.end());
		}
		else if (experiment=="DarkSide-50")
		{
			//Threshold and exposure
				double exposure = 6786.0*kg*day;
				int nThreshold=3;
			//Data
				std::vector<int>data={118643, 219893, 6131, 673, 252, 227, 198, 199, 189, 247, 230, 261,249, 329, 336, 349, 351, 352, 384, 411, 405, 461, 460, 436, 500, 546, 538, 536, 556, 583, 573, 630, 603, 635, 639, 682, 736, 755, 804, 811, 809, 882, 934, 935, 871, 965, 946, 1072, 997, 1060};
			//Inverting the Poisson CDF
				std::vector<double> Nl;
				for(unsigned int j=0;j<data.size();j++) Nl.push_back(Inv_CDF_Poisson(data[j],1.0-CL));

			std::vector<double>sigmas;
			//Compute events and limit in each bin for 1pb.
				DM_Particle dm=DM;
				dm.Set_Sigma_e(pb);
				for(unsigned int ne=nThreshold;ne<data.size();ne++)
				{
					double N_Limit = Nl[ne-1];
					double Nbin=0.0;
					//Add all shells
					for(unsigned int shell=0;shell<Shells.size();shell++)
					{					
						Nbin += 1.0/2.0*exposure*(dRdne(ne,dm,DMe_target,Shells[shell])+dRdne(ne+1,dm,DMe_target,Shells[shell]));
					}
					if(Nbin>0)
					{
						sigmas.push_back(N_Limit/Nbin*pb);
					}
					else
					{
						break;
					}
					
					
				}
			//Take the minimal value
				sigma = *std::min_element(sigmas.begin(),sigmas.end());
		}
		else
		{
			std::cerr <<"Error in Likelihood(): Experiment " <<experiment <<" not recognized."<<endl;
			std::exit(EXIT_FAILURE);
		}
		return sigma;
	}

	std::vector<std::vector<double>> Limit_Curve(DM_Particle DM,const std::string& experiment,const std::string& simID,double mMin,double mMax,int Nm,double CL,int rank)
	{
		if(rank == 0)cout <<"Compute analytic limits."<<endl;
		std::vector<std::vector<double>> output;
		//Time
			std::chrono::high_resolution_clock::time_point tS,tE;
		  	tS = std::chrono::high_resolution_clock::now();

		//Output file
			ofstream f;
			if(rank == 0) f.open("../results/"+simID+"/Limits_Analytic.txt");
		//Mass scan
			double dlogm = (Nm==1)? 0 : (log10(mMax)-log10(mMin))/(Nm-1);
			for(int i=0;i<Nm;i++)
			{
				double mDM = pow(10.0,log10(mMin)+i*dlogm);
				DM.Set_Mass(mDM);
				double sigma = Upper_Bound(DM,experiment,CL);
				output.push_back(std::vector<double> {mDM,sigma});
				if(rank==0)
				{
					f<<mDM <<"\t" <<InUnits(sigma,cm*cm)<<endl;
					cout <<"\r                                                     \r" 
					<<Round(mDM) <<" GeV\t" <<Round(InUnits(sigma,cm*cm))<<"cm^2\t("<<floor(100.0*i/Nm)<<"\%)"<<flush;
				}
			}
		//Finish
			tE= std::chrono::high_resolution_clock::now();
			double duration =1e-6*std::chrono::duration_cast<std::chrono::microseconds>( tE - tS ).count();
			if(rank==0)
			{
				cout <<"\rDone.  ("<<Round(duration)<<"sec)                                 "<<endl;
				f.close();
			}	
			return output;
	}
