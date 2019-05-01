#include "DD_Electron.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "Physics_Functions.hpp"
#include "Input_Parameters.hpp"
#include "DD_General.hpp"

//1. Import crystal/ionization form factor
	//Crystal form factor
	double FFcrystalSi[900][500];
	double FFcrystalGe[900][500];

	//Ionization form factors for different atomic shells
	std::vector<Shell> Shells;

	//Import tables
	void Import_FormFactor(const std::string& target)
	{
		//a) Semiconductor
			if(target=="Si"||target=="Ge")
			{
				//Prefactors:
				double dE = 0.1*eV;
				double wk = 2.0/137.0;
				double prefactor = (target=="Si")? 2.0*eV: 1.8*eV;
				
				//Import data:
				ifstream f;
				f.open("../detectors/Semiconductors/C."+target+"137.dat");
				for(int Ei=1;Ei<=500;Ei++)
				{
					for(int qi=1;qi<=900;qi++)
					{
						if(target=="Si")
						{
							f >> FFcrystalSi[qi-1][Ei-1];
							FFcrystalSi[qi-1][Ei-1]*=prefactor*qi/dE*wk/4.0;
						}
						else
						{
							f >> FFcrystalGe[qi-1][Ei-1]; 
							FFcrystalGe[qi-1][Ei-1]*=prefactor*qi/dE*wk/4.0;

						}
					}
				}
				f.close();
			}
		//b) Xenon
			else if(target=="Xe")
			{
				Shells.push_back( Shell ("5p","../detectors/FF_Xenon/FF_5p.txt",12.4*eV,0,-2.4,2.4,97,-1,4,101) );
				Shells.push_back( Shell ("5s","../detectors/FF_Xenon/FF_5s.txt",25.7*eV,0,-2.4,1.8,43,-1,4,51) );
				Shells.push_back( Shell ("4d","../detectors/FF_Xenon/FF_4d.txt",75.6*eV,4,-2.4,1.75,84,-1,4,101) );
				Shells.push_back( Shell ("4p","../detectors/FF_Xenon/FF_4p.txt",163.5*eV,6,-2.4,1.65,36,-1,4,101) );
				Shells.push_back( Shell ("4s","../detectors/FF_Xenon/FF_4s.txt",213.8*eV,3,-2.4,1.8,43,-1,4,51) );
			}
		//c) Argon
			else if(target=="Ar")
			{
				Shells.push_back( Shell ("3p32","../detectors/FF_Argon/FF_3p32.txt",	15.7*eV,0,-2.4,1.8,14,-1,4,51) );
				Shells.push_back( Shell ("3p12","../detectors/FF_Argon/FF_3p12.txt",	15.9*eV,0,-2.4,1.8,14,-1,4,51) );
				Shells.push_back( Shell ("3s"	,"../detectors/FF_Argon/FF_3s.txt",		29.3*eV,0,-2.4,1.8,43,-1,4,51) );
			}
		//d) Error
			else
			{
				cerr<<"Error in Import_Formfactor(): Target "<<target<<" not recognized."<<endl;
				std::exit(EXIT_FAILURE);
			}
	}


//2. Semiconductors
	//Minimal velocity
	double vMinimal_e(double q,double Ee,double mDM)
	{
		return (Ee/q+q/2.0/mDM);
	}

	//Recoil spectra
	double dRdEe(double Ee,const DM_Particle& DM,const std::string& target,double efficiency)
	{
		double dE = 0.1*eV;
		double dq = 0.02*aEM*mElectron;
		
		//Target specific parameters:
		double Mcell;
		if(target == "Si") 			Mcell = 2.0*28.08*mNucleon;
		else if (target == "Ge")	Mcell = 2.0*72.6*mNucleon;
		else
		{
			std::cerr <<"Error in dRdEe(): target " <<target <<" not recognized."<<endl;
			std::exit(EXIT_FAILURE);
		}
		double (&FFcrystal)[900][500] = (target=="Si") ? FFcrystalSi : FFcrystalGe;

		double prefactor = efficiency*rhoDM/DM.mass/Mcell*DM.sigma_e*aEM*pow(mElectron/Reduced_Mass(DM.mass,mElectron),2.0);
		
		//Integrate by summing over the tabulated form factors
		int Ei = floor(Ee/dE)-1;
		double sum=0.0;
		for(int qi=0;qi<900;qi++) sum+=1.0/(qi+1)/(qi+1)/dq*EtaFunction(vMinimal_e((qi+1)*dq,Ee,DM.mass),vEarth)*pow(DM.FormFactor((qi+1)*dq),2.0)*FFcrystal[qi][Ei];
		
		return prefactor*sum;
	}
	Interpolation dRdEe(const DM_Particle& DM,const std::string& target,double efficiency,const std::vector<double>& attenuation,const std::vector<DataPoint> &speeddata,double vCutoff)
	{
		double dE = 0.1*eV;
		double dq = 0.02*aEM*mElectron;
		//Target specific parameters
		double Mcell,Egap,eps;
		if(target == "Si")
		{
			Mcell = 2.0*28.08*mNucleon;
			Egap = 1.11*eV;
			eps = 3.6*eV;
		}			
		else if (target == "Ge")	
		{
			Mcell = 2.0*72.6*mNucleon;
			Egap = 0.67*eV;
			eps = 2.9*eV;
		}
		else
		{
			std::cerr <<"Error in dRdEe(): Target " <<target <<" not recognized."<<endl;
			std::exit(EXIT_FAILURE);
		}
		double (&FFcrystal)[900][500] = (target=="Si") ? FFcrystalSi : FFcrystalGe;

		//MC eta function
		Interpolation etaMC = EtaFunction(attenuation,speeddata,vCutoff,vEarth);

		double prefactor = efficiency*rhoDM/DM.mass/Mcell*DM.sigma_e*aEM*pow(mElectron/Reduced_Mass(DM.mass,mElectron),2.0);

		//Interpolation
		int points =500;
		std::vector<std::vector<double>> interpol_list;
		double Ethr= eps*(DMe_threshold-1)+Egap;
		for(int Ei=Ethr/dE;Ei<points;Ei++)
		{
			double sum=0.0;
			for(int qi=0;qi<900;qi++)
			{
				sum+=1.0/(qi+1)/(qi+1)/dq*etaMC(vMinimal_e((qi+1)*dq,(Ei+1)*dE,DM.mass))*pow(DM.FormFactor((qi+1)*dq),2.0)*FFcrystal[qi][Ei];
			}
			interpol_list.push_back(std::vector<double> {(Ei)*dE,prefactor*sum});
		}
		Interpolation spectrum(interpol_list);

		return spectrum;
	}

	//Total event rates
	double TotalRate(const DM_Particle& DM,int Qthreshold,double efficiency,const std::string& target)
	{
		double dE = 0.1*eV;
		double dq = 0.02*aEM*mElectron;

		//Target specific parameters
		double Mcell,Egap,eps;
		if(target == "Si")
		{
			Mcell = 2.0*28.08*mNucleon;
			Egap = 1.11*eV;
			eps = 3.6*eV;
		}			
		else if (target == "Ge")	
		{
			Mcell = 2.0*72.6*mNucleon;
			Egap = 0.67*eV;
			eps = 2.9*eV;
		}
		else
		{
			std::cerr <<"Error in TotalRate(): Target " <<target <<" not recognized."<<endl;
			std::exit(EXIT_FAILURE);
		}
		double Ethr= eps*(Qthreshold-1)+Egap;
		double (&FFcrystal)[900][500] = (target=="Si") ? FFcrystalSi : FFcrystalGe;

		//Integrating over the form factors by summing the tabulated values.
		double prefactor =efficiency*rhoDM/DM.mass/Mcell*DM.sigma_e*aEM*pow(mElectron/Reduced_Mass(DM.mass,mElectron),2.0);
		double sum=0.0;
		for(int Ei=(Ethr/dE);Ei<500;Ei++)
		{
			for(int qi=0;qi<900;qi++) sum+=dE/(qi+1)/(qi+1)/dq*EtaFunction(vMinimal_e((qi+1)*dq,(Ei+1)*dE,DM.mass),vEarth)*pow(DM.FormFactor((qi+1)*dq),2.0)*FFcrystal[qi][Ei];
		}

		return prefactor*sum;
	}
	double TotalRate(const DM_Particle& DM,int Qthreshold,double efficiency,const std::string& target,const std::vector<double>& attenuation,const std::vector<DataPoint> &speeddata,double vCutoff)
	{
		double dE = 0.1*eV;
		double dq = 0.02*aEM*mElectron;
		//Target specific parameters
		double Mcell,Egap,eps;
		if(target == "Si")
		{
			Mcell = 2.0*28.08*mNucleon;
			Egap = 1.11*eV;
			eps = 3.6*eV;
		}			
		else if (target == "Ge")	
		{
			Mcell = 2.0*72.6*mNucleon;
			Egap = 0.67*eV;
			eps = 2.9*eV;
		}
		else
		{
			std::cerr <<"Error in TotalRate(): Target " <<target <<" not recognized."<<endl;
			std::exit(EXIT_FAILURE);
		}

		//Energy threshold
		double Ethr= eps*(Qthreshold-1)+Egap;
		
		//Crystal form factor
		double (&FFcrystal)[900][500] = (target=="Si") ? FFcrystalSi : FFcrystalGe;
		
		//MC eta function
		Interpolation etaMC = EtaFunction(attenuation,speeddata,vCutoff,vEarth);
		
		//Integrating over the form factors by summing the tabulated values.
		double prefactor = efficiency*rhoDM/DM.mass/Mcell*DM.sigma_e*aEM*pow(mElectron/Reduced_Mass(DM.mass,mElectron),2.0);
		double sum=0.0;
		for(int Ei=(Ethr/dE);Ei<500;Ei++)
		{
			for(int qi=0;qi<900;qi++) sum+=dE/(qi+1)/(qi+1)/dq*etaMC(vMinimal_e((qi+1)*dq,(Ei+1)*dE,DM.mass))*pow(DM.FormFactor((qi+1)*dq),2.0)*FFcrystal[qi][Ei];
		}
		
		return prefactor*sum;
	}


//3. Atomic shell struct	
	//Constructor
	Shell::Shell(std::string nm,std::string filepath,double Eb,int Ns,double lnkMin,double lnkMax,int Nk,double lnqMin,double lnqMax,int Nq)
	{
		name = nm;
		Ebinding = Eb;
		neSecondary	=Ns;
		//Ionization form factor
			//Data range in k and q (logarithmically in units of aEM mElectron)
			logkMin=lnkMin;
			logkMax=lnkMax;
			logqMin=lnqMin;
			logqMax=lnqMax;
			nk=Nk;
			nq=Nq;
			dlogk=(logkMax-logkMin)/(nk-1);
			dlogq=(logqMax-logqMin)/(nq-1);

			//Import form factor table
			ifstream inputfile;
			inputfile.open(filepath);
			if (inputfile.good())
			{
				for(int ki=0;ki<Nk;ki++)
				{
					std::vector<double> qVector;
					for(int qi=0;qi<Nq;qi++)
					{
						double x;
		            	inputfile >>x;
						qVector.push_back(x);
					}
					FormFactor.push_back(qVector);
				}
		        inputfile.close();
	   		}
	   		else 
   			{
   				cerr <<"Error in Shell(): File " <<filepath <<" not found"<<endl;
   				std::exit(EXIT_FAILURE);
   			}
	}

//4. Liquid noble gas experiments
	//Recoil energy spectrum
	double dRdlog10Ee(double Ee,const DM_Particle& DM,const std::string& target,const Shell& shell)
	{
		double A = (target=="Xe")? 131.0 : 40.0;
		double mA = A*mNucleon;
		double qref = aEM*mElectron;
		double prefactor =log(10)* 1.0/mA*rhoDM/DM.mass*DM.sigma_e/8/pow(Reduced_Mass(DM.mass,mElectron),2.0);
		double sum=0.0;
		int ki =std::round ( (1.0/2.0*log(2.0*mElectron*Ee/qref/qref)-shell.logkMin)/shell.dlogk );
		if(ki>=shell.nk||ki<0) return 0.0;
		for(int qi=0;qi<shell.nq;qi++)
		{
			double q = qref*exp(shell.logqMin+(qi)*shell.dlogq);
			double vMin = (shell.Ebinding+Ee)/q+q/2.0/DM.mass;
			sum+=shell.dlogq*q*q*EtaFunction(vMin,vEarth)*pow(DM.FormFactor(q),2.0)*shell.FormFactor[ki][qi];
		}
		return prefactor*sum;
	}

	//Spectrum number of electrons
	double PDFne(const std::string& target,unsigned int ne,double Ee,const Shell& shell)
	{
		double fR=0.0;
		double NxNi = 0.2;
		double W = (target=="Xe")? 13.8*eV : 19.6*eV;
		double fe = (1.0-fR)/(1.0+NxNi);
		double neMax = shell.neSecondary+floor(Ee/W);
		return PDF_Binomial(neMax,fe,ne-1);
	}

	double dRdne(unsigned int ne,const DM_Particle& DM,const std::string& target,const Shell& shell)
	{
		double A = (target=="Xe")? 131.0 : 40.0;
		double mA = A*mNucleon;
		double qref = aEM*mElectron;
		double prefactor = 1.0/mA*rhoDM/DM.mass*DM.sigma_e/8/pow(Reduced_Mass(DM.mass,mElectron),2.0);
		double sum=0.0;
		for(int ki=0;ki<shell.nk;ki++)
		{
			double Ee = pow(qref*exp(shell.logkMin+(ki)*shell.dlogk),2.0)/2.0/mElectron;
			for(int qi=0;qi<shell.nq;qi++)
			{
				double q = qref*exp(shell.logqMin+(qi)*shell.dlogq);
				double vMin = (shell.Ebinding+Ee)/q+q/2.0/DM.mass;
				sum+=2.0*shell.dlogk*shell.dlogq*PDFne(target,ne,Ee,shell)*q*q*EtaFunction(vMin,vEarth)*pow(DM.FormFactor(q),2.0)*shell.FormFactor[ki][qi];
			}
		}
		return prefactor*sum;
	}

	std::vector<double> dRdne(const DM_Particle& DM,const std::string& target,const std::vector<Shell>& shells,const std::vector<double> &attenuation,const std::vector<DataPoint> &speeddata,double vCutoff)
	{
		std::vector<double> output;
		double A = (target=="Xe")? 131.0 : 40.0;
		double mA = A*mNucleon;
		double qref = aEM*mElectron;
		double prefactor = 1.0/mA*rhoDM/DM.mass*DM.sigma_e/8.0/pow(Reduced_Mass(DM.mass,mElectron),2.0);

		Interpolation etaMC = EtaFunction(attenuation,speeddata,vCutoff,vEarth);

		for(int ne=1;ne<=50;ne++)
		{
			double sum=0.0;
			for(unsigned int s=0;s<shells.size();s++)
			{
				for(int ki=0;ki<shells[s].nk;ki++)
				{
					double Ee = pow(qref*exp(shells[s].logkMin+(ki)*shells[s].dlogk),2.0)/2.0/mElectron;
					for(int qi=0;qi<shells[s].nq;qi++)
					{
						double q = qref*exp(shells[s].logqMin+(qi)*shells[s].dlogq);
						double vMin = (shells[s].Ebinding+Ee)/q+q/2.0/DM.mass;
						sum+=2.0*shells[s].dlogk*shells[s].dlogq*PDFne(target,ne,Ee,shells[s])*q*q*etaMC(vMin)*pow(DM.FormFactor(q),2.0)*shells[s].FormFactor[ki][qi];
					}

				}
			}
			output.push_back(prefactor*sum);
		}

		return output;
	}
	
	//Spectrum number of PE
	double dRdnPE(unsigned int nPE,const DM_Particle& DM,const std::string& target,double muPE,double sigPE,const std::vector<Shell>& shells)
	{
		double output=0.0;
		for(int ne=1;ne<=15;ne++)
		{
			for(unsigned int i=0;i<shells.size();i++)
			{
				output+= PDF_Gauss(nPE,muPE*ne,sqrt(ne)*sigPE)*dRdne(ne,DM,target,shells[i]);
			}	
		}
		return output;
	}
	
	std::vector<double> dRdnPE(double muPE,double sigPE,const std::vector<double>& dRdn)
	{
		std::vector<double> output;
		for(int nPE=1;nPE<250;nPE++)
		{
			double sum=0.0;
			for(unsigned int ne=0;ne<dRdn.size();ne++)
			{
				sum+= PDF_Gauss(nPE,muPE*(ne+1),sqrt(ne+1)*sigPE)*dRdn[ne];
			}
			output.push_back(sum);

		}
		return output;
	}

	std::vector<double> PE_Trigger_Eff;
	std::vector<double> PE_Acc_Eff;
	double dNdnPE(unsigned int nPE,const DM_Particle& DM,const std::string& target,const std::string& detector,const std::vector<Shell>& shells)
	{
		double exposure,muPE,sigPE,TotalEff;

		if(detector=="XENON10e")
		{
			exposure=15*kg*day;
			muPE=27.0;
			sigPE=6.7;
			//Efficiency
				//Import trigger efficiency (taken from [arXiv:1206.2644])
				if(PE_Trigger_Eff.size()==0) PE_Trigger_Eff = Read_List("../detectors/XENON10e/PE_Trigger_Efficiency.txt");
				//Flat acceptance efficiency
				double PE_Acc_Eff_flat=0.92;
				TotalEff = PE_Acc_Eff_flat*PE_Trigger_Eff[nPE-1];
		}
		else if(detector=="XENON100e")
		{
			exposure=30*kg*year;
			muPE=19.7;
			sigPE=6.2;
			//Efficiencies
			//Import trigger and acceptance efficiency (taken from [arXiv:1605.06262])
			if(PE_Trigger_Eff.size()==0)
			{
				PE_Trigger_Eff = Read_List("../detectors/XENON100e/PE_Trigger_Efficiency.txt");
				PE_Acc_Eff = Read_List("../detectors/XENON100e/PE_Acceptance_Efficiency.txt");
			}
			TotalEff=PE_Trigger_Eff[nPE]*PE_Acc_Eff[nPE];
		}
		else
		{
			cerr <<"Error in dNdnPE(): Detector "<<detector <<" not recognized."<<endl;
			std::exit(EXIT_FAILURE);
		}
		return TotalEff*exposure*dRdnPE(nPE,DM,target,muPE,sigPE,shells);
	}

	std::vector<double> dNdnPE(const std::string& detector,const std::vector<double>& dRdnpe)
	{
		std::vector<double> output;
		double exposure;

		// Interpolation efficiency();
		if(detector=="XENON10e")
		{
			exposure=15*kg*day;
			//Efficiencies
			//Import trigger efficiency (taken from [arXiv:1206.2644])
			std::vector<double> PE_Trigger_Eff = Read_List("../detectors/XENON10e/PE_Trigger_Efficiency.txt");
			//Flat acceptance efficiency
			double PE_Acc_Eff=0.92;
			for(unsigned int i=0;i<dRdnpe.size();i++)
			{
				output.push_back(PE_Acc_Eff*PE_Trigger_Eff[i]*exposure*dRdnpe[i]);
			}
		}
		else if(detector=="XENON100e")
		{
			exposure=30*kg*year;
			//Efficiencies
			//Import trigger and acceptenace efficiency (taken from [arXiv:1605.06262])
			std::vector<double> PE_Trigger_Eff = Read_List("../detectors/XENON100e/PE_Trigger_Efficiency.txt");
			std::vector<double> PE_Acc_Eff = Read_List("../detectors/XENON100e/PE_Acceptance_Efficiency.txt");
			
			for(unsigned int i=0;i<dRdnpe.size();i++)
			{
				output.push_back(PE_Trigger_Eff[i]*PE_Acc_Eff[i]*exposure*dRdnpe[i]);
			}
		}
		else
		{
			cerr <<"Error in dNdnPE(): Detector "<<detector <<" not recognized."<<endl;
			std::exit(EXIT_FAILURE);
		}
		return output;
	}