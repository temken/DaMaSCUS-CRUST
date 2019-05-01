#include "Trajectory_Simulation.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>

#include "Input_Parameters.hpp"
#include "DD_Nucleus.hpp"

using namespace std::placeholders;


//Cut off speed for a given experiment
	double Cutoff_Speed(const std::string& experiment,double mDM,int rank)
	{
		double vcut=0;
		if(experiment=="DAMIC")
		{
			double Ethr = 0.55*keV;
			double A = 28.0;
			vcut = vMinimal_N(Ethr,mDM,A);
		}
		else if(experiment == "XENON1T")
		{
			double Ethr = 5*keV;
			double A = 131.0;
			vcut = vMinimal_N(Ethr,mDM,A);
		}
		else if(experiment =="CRESST-II")
		{
			double Ethr = 307*eV;
			double resolution = Ethr/5.0;
			double Emin = Ethr-3.0*resolution;
			double A[3] = {16.,40.,184.};
			std::vector<double> v = {vMinimal_N(Emin,mDM,A[0]),vMinimal_N(Emin,mDM,A[1]),vMinimal_N(Emin,mDM,A[2])};
			vcut = *std::min_element(v.begin(),v.end());
		}
		else if(experiment =="CRESST-surface")
		{
			double Ethr = 19.7*eV;
			double resolution = 3.74*eV;
			double Emin = Ethr-3.0*resolution;
			double A[3] = {16.,27.};
			std::vector<double> v = {vMinimal_N(Emin,mDM,A[0]),vMinimal_N(Emin,mDM,A[1])};
			vcut = *std::min_element(v.begin(),v.end());
		}
		else if(experiment == "Semiconductor"||experiment=="SENSEI-surface"||experiment=="SENSEI"||experiment=="SuperCDMS"||experiment=="DAMIC-M")
		{
			double eps = (DMe_target=="Si") ? 3.6*eV : 2.9*eV;
			double Egap = (DMe_target=="Si") ? 1.11*eV : 0.67*eV;
			double Ethr = eps*(DMe_threshold-1.0)+Egap;
			vcut = sqrt(2.0*Ethr/mDM);
		}
		else if(experiment == "XENON10e" || experiment == "XENON100e")
		{
			vcut = sqrt(2*12.4*eV/mDM);
		}
		else if(experiment == "DarkSide-50")
		{
			vcut = sqrt(2*15.7*eV/mDM);
		}
		else
		{
			std::cerr <<"Error in Cutoff_Speed(): Experiment " <<experiment <<" not recognized."<<endl;
			std::exit(EXIT_FAILURE);
		}
		if(vcut>(vEarth+vesc))
		{
			if(rank==0) cout <<"Warning: vCutoff>vEsc+vEarth for mDM = " <<mDM <<"GeV."<<endl;
			vcut = -1;
		}
		return vcut;
	}
//Initial Condition Generator
	Event InitialCondition(double tIni,Vector3D xIni,std::mt19937& PRNG,double vCutoff)
	{	
		//Speed with rejection sampling
			double ymax = (vCutoff>350*km/sec)? SpeedDistribution(vCutoff,vEarth) : 1000.0;
			auto pdf = std::bind(SpeedDistribution,_1,vEarth);
			double v = Rejection_Sampling(pdf,vCutoff,(vesc+vEarth),ymax,PRNG);
		//Velocity Direction
			double cosalpha = sqrt(ProbabilitySample(PRNG));
			Vector3D eV(sqrt(1.0-cosalpha*cosalpha),0,-cosalpha);
			Vector3D vIni=v*eV;
		return Event (tIni,xIni,vIni);
	}

//Propagate a particle until it scatters, reaches the detector or gets back into space
	void Propagate(Event& event,double &weight,std::mt19937& PRNG,double logXi=-1.0)
	{
		//If this is the first iteration we set log(1-xi) and the weight.
			if(logXi<0.0)
			{
				double Xi = ProbabilitySample(PRNG);
				logXi =(-1.0)*log(1.0-Xi);
				weight=(1.0+IS_MFP)*pow(1.0-Xi,IS_MFP);

			}
		//Initial layer and potentially next layer
			int layer_old = Current_Layer(event);
			double vz= event.Velocity()[2];
			int next_layer = (vz>0)? layer_old-1 : layer_old+1;
		//Distance to next layer boundary
			double v = event.Speed();
			double t_exit;
			if(next_layer==0) 					t_exit = -event.Position()[2]/vz;
			else if(next_layer==Detector_Index)	t_exit = (-Detector_Depth-event.Position()[2])/vz;
			else if(vz>0)						t_exit = (-Layers[layer_old-1].depth-event.Position()[2])/vz;
			else								t_exit = (-Layers[next_layer-1].depth-event.Position()[2])/vz;
			double d = t_exit*v;
		//Mean free path
			double mfp = Layers[layer_old-1].MFP(v)*(1.0+IS_MFP);
		//a) The particle did not leave the layer
			if(logXi<d/mfp)
			{
				double L = logXi * mfp;
				double dt = L/v;
				event.IncreaseTime(dt);
				event.IncreasePosition(dt*event.Velocity());
			}
		//b) The particle reaches space or the detector
			else if(next_layer == 0 || next_layer == Detector_Index)
			{
				double L = d+10.0*meter;
				double dt = L/v;
				event.IncreaseTime(dt);
				event.IncreasePosition(dt*event.Velocity());
			}
		//c) The particle reaches the next layer
			else
			{
				double dt = (d+1.0e-6*meter)/v;
				event.IncreaseTime(dt);
				event.IncreasePosition(dt*event.Velocity());
				if(Current_Layer(event)==layer_old) cout <<"Transition from " <<layer_old <<" to " <<layer_old <<" at z=" <<event.Position()[2]/meter <<"m."<<endl;
				logXi-=d/mfp;
				Propagate(event,weight,PRNG,logXi);
			}
	}

//Perform the scattering.
	double PDF_ScatterAngle(double cosalpha,const DM_Particle& DM,const std::vector<double>& target,double vDM,double deltaIS)
	{
		double pdf;
		double Z = target[0];
		double A = target[1];

		double mu = Reduced_Mass(DM.mass,A*mNucleon);
		double q2max = 4.0*mu*mu*vDM*vDM ;
		
		if(!DM.ldm)
		{
			double q = sqrt(q2max/2.0*(1.0-cosalpha));
			//Screening is included automatically.
			double sigmaTot = DM.Sigma_Tot(Z,A,vDM);
			pdf= (q2max/sigmaTot*DM.dSdq2(q,Z,A,vDM))/2.0; 
		}
		else
		{
			pdf = 1.0/2.0;
			double a2=pow(Thomas_Fermi_Radius(Z),2.0);
			double x =a2*q2max;
			//General interaction
				if(DM.formfactor == "General")
				{
					double m2 = pow(DM.mMediator,2.0);
					if(DM.screening)
					{
						pdf*=(4.0*pow(1.0 - cosalpha,2.0)*(m2 + q2max)*(1.0 + x)*pow(m2*x-q2max,3))/(pow(2.0*m2 + q2max - cosalpha*q2max,2.0)*pow(2.0 + (1.0 - cosalpha)*x,2.0)*(-2.0*m2*q2max - pow(q2max,2.0) + pow(m2,2.0)*x*(2.0 + x) + 2.0*m2*(m2 + q2max)*(1.0 + x)*log((m2 + q2max)/(m2 + m2*x))));
					}
					else
					{
						pdf*=4.0*m2*(m2+q2max)/pow(2.0*m2+q2max-cosalpha*q2max,2.0);
					}
				}
			//Contact interaction
				else if(DM.formfactor == "Contact" && DM.screening)
				{
					pdf*=x*x*x/4.0*(1.0+x)/(x*(2.0+x)-2.0*(1.0+x)*log(1.0+x))*pow((1.0-cosalpha),2.0)/pow((1.0+x/2.0*(1.0-cosalpha)),2.0);

				}
			//Electric dipole interaction
				else if(DM.formfactor == "Electric-Dipole")
				{
					pdf*=x*x/2.0*(1.0+x)/((1.0+x)*log1p(x)-x)*(1.0-cosalpha)/pow((1.0+x/2.0*(1.0-cosalpha)),2.0);
				}
			//Long range interaction
				else if(DM.formfactor== "Long-Range")
				{
					pdf*=(1.0+x)/pow((1.0+x/2.0*(1.0-cosalpha)),2.0);
				}
			return pdf;
		}

		return pdf + deltaIS/2.0*cosalpha;
	}
	double Sample_ScatteringAngle(const DM_Particle& DM,const std::vector<double>& target,double vDM,double delta, double& weight,std::mt19937& PRNG)
	{
		weight = 1.0;
		double cosalpha=0.0;
		if(delta<1e-6 || PDF_ScatterAngle(-1.0,DM,target,vDM,delta)<0.0) delta = 0;
		
		//Heavy DM. Include Helm form factor
			if(!DM.ldm)
			{
				
				auto pdf = std::bind(PDF_ScatterAngle,_1,DM,target,vDM,delta);
				//TODO ymax has to be found for general pdfs, or a whole new sampling method is necessary.
				cosalpha = Rejection_Sampling(pdf,-1.0,1.0,2.0,PRNG);
				weight = PDF_ScatterAngle(cosalpha,DM,target,vDM,0.0)/PDF_ScatterAngle(cosalpha,DM,target,vDM,delta);
			}
		//Light DM
			else
			{
				//a) General interaction
					if(DM.formfactor=="General")
					{
						double Xi=ProbabilitySample(PRNG);
						double Z = target[0];
						double A = target[1];
						if(delta<1e-6)
						{
							if(DM.screening)
							{
								double mu = Reduced_Mass(DM.mass,A*mNucleon);
								double q2max = 4.0*mu*mu*vDM*vDM ;
								double a2=pow(Thomas_Fermi_Radius(Z),2.0);
								double m2 = DM.mMediator*DM.mMediator;
								std::function<double(double)> cdf = [Xi,m2,q2max,a2](double cosa)
								{
									return Xi - (-(((m2 + q2max)*(1.0 + a2*q2max)*((-1.0 + a2*m2)*(((1.0 + cosa)*q2max)/((1.0 + a2*q2max)*(-2.0 + a2*(-1.0 + cosa)*q2max)) + pow(m2,2.0)*(1.0/(m2 + q2max) - 2.0/(2.0*m2 + q2max - cosa*q2max))) + 2.0*m2*(log1p(a2*q2max) + log((2.0*m2 + q2max - cosa*q2max)/(m2 + q2max)) - log(2.0 + a2*q2max - a2*cosa*q2max))))/((-1.0 + a2*m2)*q2max*(q2max + m2*(2.0 + a2*q2max)) + 2.0*m2*(m2 + q2max)*(1.0 + a2*q2max)*log((m2 + q2max)/(m2 + a2*m2*q2max)))));
								};
								cosalpha=Find_Root(cdf,-1.0,1.0,1e-6);
								weight=1.0;
							}
							else
							{
								double m2 = DM.mMediator*DM.mMediator;
								double mu = Reduced_Mass(DM.mass,A*mNucleon);
								double q2max = 4.0*mu*mu*vDM*vDM ;
								cosalpha = (m2*(2.0*Xi-1.0)+q2max*Xi)/(m2+q2max*Xi);
								weight=1.0;	
							}
						}
						else
						{
							std::cerr <<"Error in Sample_ScatteringAngle(): Importance sampling for the scattering angle should not be used with the general DM form factor. Set is_angle to zero in the cfg file."<<endl;
							std::exit(EXIT_FAILURE);
						}
							
					}
				//b) Contact interaction
					else if(DM.formfactor=="Contact")
					{
						if(DM.screening)
						{
							double Xi=ProbabilitySample(PRNG);
							double Z = target[0];
							double A = target[1];
							double mu = Reduced_Mass(DM.mass,A*mNucleon);
							double q2max = 4.0*mu*mu*vDM*vDM ;
							double a2=pow(Thomas_Fermi_Radius(Z),2.0);
							double x =a2*q2max;
							std::function<double(double)> cdf = [Xi,x](double cosa)
							{
								return Xi-(((1.0 + x)*(x + cosa*x - 2.0/(1.0 + x) + 4.0/(2.0 + x - cosa*x) - 4.0*log(2.0*(1.0 + x)) + 4.0*log(2.0 + x - cosa*x)))/(2.0*(x*(2.0 + x) - 2.0*(1.0 + x)*log(1.0 + x))));
							};
							cosalpha=Find_Root(cdf,-1.0,1.0,1e-4);
							weight=1.0;
						}
						else
						{
							double Xi=ProbabilitySample(PRNG);
							if(delta < 1e-6) 
							{
								cosalpha= 2.0*Xi-1.0;
								weight=1.0;
							}
							else 
							{
								cosalpha = 1.0/delta * (-1.0+sqrt((1-delta)*(1-delta)+4*delta*Xi));
								weight = 1.0/(1.0+delta*cosalpha);
							}
						}
					}
				//c) Long range interaction
					else if(DM.formfactor=="Long-Range")
					{
						double Xi=ProbabilitySample(PRNG);
						double Z = target[0];
						double A = target[1];
						double mu = Reduced_Mass(DM.mass,A*mNucleon);
						double q2max = 4.0*mu*mu*vDM*vDM ;
						double a2=pow(Thomas_Fermi_Radius(Z),2.0);
						double x =a2*q2max;
						if(delta<1e-6)
						{
							cosalpha=(Xi*(x+2.0)-1.0)/(1.0+x*Xi);
							weight=1.0;
						}
						else
						{
							std::cerr <<"Error in Sample_ScatteringAngle(): Importance sampling for the scattering angle should not be used with the 1/q^2 DM form factor. Set is_angle to zero in the cfg file."<<endl;
							std::exit(EXIT_FAILURE);
						}
					}
				//d) Electric dipole interaction
					else if(DM.formfactor=="Electric-Dipole")
					{
						if(delta<1e-6)
						{
							double Xi=ProbabilitySample(PRNG);
							double Z = target[0];
							double A = target[1];
							double mu = Reduced_Mass(DM.mass,A*mNucleon);
							double q2max = 4.0*mu*mu*vDM*vDM ;
							double a2=pow(Thomas_Fermi_Radius(Z),2.0);
							double x =a2*q2max;
					
							std::function<double(double)> cdf = [Xi,x](double cosa)
							{
								return Xi+((1.0 + x)*(-(1.0/(1.0 + x)) + 2.0/(2.0 + x - cosa*x) + log((1.0 + 0.5*(1.0 - cosa)*x)/(1.0 + x)))/(-x + (1.0 + x)*log1p(x)));
							};
							cosalpha=Find_Root(cdf,-1.0,1.0,1e-6);
							weight=1.0;
						}
						else
						{
							std::cerr <<"Error in Sample_ScatteringAngle(): Importance sampling for the scattering angle should not be used with the 1/q DM form factor. Set is_angle to zero in the cfg file."<<endl;
							std::exit(EXIT_FAILURE);
						}
					}

			}
		return cosalpha;
	}

	void Scatter(const DM_Particle& DM,Event& event, double& weight,std::mt19937& PRNG)
	{
		Vector3D vIni = event.Velocity(), ev = vIni.normalized();
		double vDM = vIni.norm();
		//Find target
		int layer = Current_Layer(event);
		std::vector<double> target = Layers[layer-1].Sample_Target(vDM,PRNG);//prob[ScatterNucleus(prob,PRNG)];
		//Find cos alpha
			double cosalpha = Sample_ScatteringAngle(DM,target,vDM,IS_Angle,weight,PRNG);
		//Construction of n, the unit vector pointing into the direction of vfinal.
			double cosphi=2.0*ProbabilitySample(PRNG)-1.0;
			double sinphi = sqrt(1.0-cosphi*cosphi);
			double sinalpha = sqrt(1.0-cosalpha*cosalpha);
			double aux = sqrt(1.0-pow(ev[2],2.0));
			Vector3D n
				(
				cosalpha*ev[0]+(sinalpha*(-ev[0]*ev[2]*cosphi+ev[1]*sinphi))/aux,
				cosalpha*ev[1]+(sinalpha*(-ev[1]*ev[2]*cosphi-ev[0]*sinphi))/aux,
				cosalpha*ev[2]+aux*cosphi*sinalpha
				);
		//Find post-scattering velocity, see Landau Lifshitz.
			double mNucleus = target[1]*mNucleon;
			Vector3D vNew = mNucleus/(DM.mass+mNucleus)*vDM*n+DM.mass/(DM.mass+mNucleus)*vIni;
		//TEST:DEACTIVATE DEFLECTIONS
			// vNew = vNew.norm()*vIni.normalized();
		//Set new velocity
			event.SetVelocity(vNew);
	}

//Summarizing parameters of a trajectory
	Result::Result()
	{
		success = false;
		nScattering = 0;
		max_horizontal_distance=0.0;
		weight = 1.0;
	}
	Result::Result(double w)
	{
		success = false;
		nScattering = 0;
		max_horizontal_distance=0.0;
		weight = w;
	}
	void Result::Summary()
	{
		cout 	<<endl<<"Success\t\t\t" <<success<<endl
				<<"Final depth[m]\t\t" <<InUnits(final_event.Position()[2],meter)<<endl
				<<"Final speed[km/sec]\t" <<InUnits(final_event.Speed(),km/sec)<<endl
				<<"nScattering\t\t"<<nScattering<<endl
				<<"lMax[m]\t\t\t" <<max_horizontal_distance/meter<<endl
				<<"weight\t\t\t" <<weight <<endl<<endl;
	}
//Simulate one particle track:
	std::vector<Result> Simulate_Trajectory(DM_Particle& DM, Event IC, double vMin,std::mt19937& PRNG,Result result,std::string split_ID)
	{
		std::vector<Result> results;
		Event event = IC;
		bool simulate = true;
		//1. Check if initial conditions are above ground
			if(event.Position()[2]>0)
			{
				double vz = event.Velocity()[2];
				if(vz>0)
				{
					std::cerr<<"Warning: Invalid IC. Velocity points up." <<endl;
					simulate = false;
				}
				else
				{
					double d = event.Position()[2];
					double tEntry=-(d+1.0e-6*meter)/vz;
					event.IncreaseTime(tEntry);
					event.IncreasePosition(tEntry*event.Velocity());
				}
			}
		//2. Simulation
			double I_1=Importance(event.Position()[2]);
			while(simulate)
			{
				//Next position
					double Weight_MFP;
					Propagate(event,Weight_MFP,PRNG);
					result.weight*=Weight_MFP;

				//Check the maximum horizontal distance
					double z = event.Position()[2];
				//a) Particle got reflected
					if(z>0.0) 
					{
						simulate = false;
						result.final_event=event;
						results.push_back(result);
					}
				//b) Particle scatters
					else if(z>(-1.0)*Detector_Depth)
					{
						double I_2 =Importance(z);
						//Splitting
						if(I_2 > I_1)
						{
							//Number of splits:
							int splits;
							double nu = I_2/I_1;
							double Delta = nu - (int) nu;
							if(Delta<1e-4)  splits = nu;
							else
							{
								double xi = ProbabilitySample(PRNG);
								if(Delta <=xi)	splits = (int) nu;
								else			splits = (int) (nu+1.0);
							}
							//Splitted weight
							result.weight = 1.0*result.weight/splits;
							for(int i=0;i<splits;i++)
							{
								Event split_event = event;
								Result split_result = result;

								double Weight_Angle;
								Scatter(DM,split_event,Weight_Angle,PRNG);
								split_result.weight*=Weight_Angle;
								split_result.nScattering++;
								if(split_event.Speed()<vMin)  
								{
									simulate = false;
									split_result.final_event=split_event;
									results.push_back(split_result);
								}
								else
								{
									std::string split_ID_new = split_ID + std::to_string(i+1);
									simulate = false;
									std::vector<Result> split_results =  Simulate_Trajectory(DM,split_event,vMin,PRNG,split_result,split_ID_new);
									results.insert(results.end(),split_results.begin(),split_results.end());
								}
								
							}
						}
						//Russian Roulette
						else if(I_2 < I_1)
						{
							double pSurv = 1.0*I_2/I_1;
							if(Survive_Russian_Roulette(pSurv,PRNG))
							{
								result.weight = 1.0*result.weight/pSurv;
								double Weight_Angle;
								Scatter(DM,event,Weight_Angle,PRNG);
								result.weight*=Weight_Angle;
								result.nScattering++;
								//Check if we fall below the speed cutoff
								if(event.Speed()<vMin)  
								{
									simulate = false;
									result.final_event=event;
									results.push_back(result);
								}

							}
							else
							{
								simulate = false;
								result.final_event=event;
								results.push_back(result);

							}

						}
						else
						{
							double Weight_Angle;
							Scatter(DM,event,Weight_Angle,PRNG);
							result.weight*=Weight_Angle;
							result.nScattering++;
							//Check if we fall below the speed cutoff
							if(event.Speed()<vMin)  
							{
								simulate = false;
								result.final_event=event;
								results.push_back(result);
							}
						}
						I_1=I_2;	
					}
				//c) Particle reaches the detector
					else
					{
						if(event.Speed()>vMin)	result.success  = true;
						simulate = false;
						result.final_event=event;
						results.push_back(result);
					}
			}
		return results;
	}
