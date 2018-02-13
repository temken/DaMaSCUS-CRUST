#include "Trajectory_Simulation.hpp"

#include <iostream>

#include "General_Utilities.hpp"

using namespace std::placeholders;

//Initial Condition Generator
	Event InitialCondition(double tIni,Eigen::Vector3d& xIni,std::mt19937& PRNG,double vCutoff)
	{	
		//Speed
			bool success=false;
			double v,y;
			double ymax=1000.0;
			if(vCutoff>350*km/sec) ymax = SpeedDistribution(vCutoff,vEarth);
			while(success==false)
			{
				v=vCutoff+ProbabilitySample(PRNG)*(vesc+vEarth-vCutoff);
				y=ProbabilitySample(PRNG)*ymax;

				if(y<=SpeedDistribution(v,vEarth)&&v>vCutoff)
				{
					success=true;
				}
			}
		//Velocity Direction
			double cosalpha = sqrt(ProbabilitySample(PRNG));
			Eigen::Vector3d eV(sqrt(1.0-cosalpha*cosalpha),0,-cosalpha);
			Eigen::Vector3d vIni=v*eV;
		Event output(tIni,xIni,vIni);
		return output;
	}

//Propagate a particle until it scatters, reaches the detector or gets back into space
	void Free_Propagation(Event& event,double& weight,std::vector<std::vector<double>>& prob,std::mt19937& PRNG,double logXi)
	{
		//If this is the first iteration we set log(1-xi) and the weight.
			if(logXi<0)
			{
				double Xi = ProbabilitySample(PRNG);
				logXi =(-1.0)*log(1.0-Xi);
				weight=(1.0+IS_MFP)*pow(1.0-Xi,IS_MFP);

			}
		//Current Layer
			Eigen::Vector3d x0 = event.Position();
			Eigen::Vector3d v = event.Velocity();
			int layer_index = Current_Layer(x0);
			
		//Potential free travelling distance
			double mfp = Layers[layer_index-1].Get_MFP();
			double L = mfp*(1.0+IS_MFP)*logXi;
			Eigen::Vector3d x = x0 + L * v.normalized();
		//Check if we passed a layer boundary
			if(Current_Layer(x)==layer_index)
			{
				event.SetPosition(x);
				event.IncreaseTime(L/v.norm());
			}
		//Transition to next layer
			else
			{
				Event exitevent = Leave_Layer(event);
				if(Current_Layer(event)==Current_Layer(exitevent)) cout <<"Transition from " <<Current_Layer(event) <<" to " <<Current_Layer(exitevent) <<" at z=" <<exitevent.Position()[2]/meter <<"m."<<endl;
				layer_index = Current_Layer(exitevent);
				L = (x0-exitevent.Position()).norm();
				//Check if we reached the detector or space without scattering
					if(layer_index == 0 || layer_index == Detector_Index)
					{
						event = exitevent;

					}
					else
					{
						prob = Layers[layer_index - 1].Get_Probability();
						logXi-= L/mfp/(1.0+IS_MFP);
						if(logXi<0.0) 
						{
							cout <<"\nError in Free_Propagation(): logXi < 0 (" <<logXi <<")"<<endl;
						}
						event = exitevent;
						Free_Propagation(event,weight,prob,PRNG,logXi);
					}
			}
	}

//Find the scattering nucleus species:
	std::vector<double> ScatterNucleus(std::vector<std::vector<double>>& prob,std::mt19937& PRNG)
	{
		std::vector<double> output;
		//Random number
			double xi = ProbabilitySample(PRNG);
			double sum=0.0;
		//Which element in the element?
		for(unsigned int i=0;i<prob.size();i++)
		{
			sum+=prob[i][0];
			if(sum>xi)
			{
				output=prob[i];
				break;
			}
		}
		return output;
	}

//Perform the scattering.
	double PDF_ScatterAngle(double cosalpha,double A,double mDM,double vDM,double deltaIS)
	{
		double mA = A*mNucleon;
		double mu = Reduced_Mass(mDM,mA);
		double q = sqrt(2.0*mu*mu*vDM*vDM*(1.0-cosalpha));
		//Numerator
			double num = pow(FormFactor(q,A),2.0);
		//Denominator, integrate the form factor
			std::function<double(double)> ff = [A,mu,vDM](double x)
			{
				double Q = sqrt(2.0*mu*mu*vDM*vDM*(1.0-x));
				return pow(FormFactor(Q,A),2.0);
			};
			double epsilon =1e-6*ff(1.0);
			double den = AdaptiveSimpsons(ff,-1.0,1.0,epsilon);
		return num/den+deltaIS/2.0*cosalpha;

	}
	void Scatter(std::vector<std::vector<double>>& prob, Event& event, double& weight,double mDM,std::mt19937& PRNG,bool ldm)
	{
		Eigen::Vector3d n, vIni = event.Velocity(), ev = vIni.normalized();
		double A = ScatterNucleus(prob,PRNG)[2];
		double cosalpha=0.0,mNucleus = A*mNucleon;
		// cout <<mNucleus/mNucleon <<endl;
		//Find cos alpha
			if(ldm)
			{
				if(IS_Angle<1e-4)
				{
					cosalpha= 2.0*ProbabilitySample(PRNG)-1.0;
					weight =1.0;
				}
				else
				{
					cosalpha = 1.0/IS_Angle * (-1.0+sqrt((1-IS_Angle)*(1-IS_Angle)+4*IS_Angle*ProbabilitySample(PRNG)));
					weight = 1.0/(1.0+IS_Angle*cosalpha);
				}
			}
			else
			{
				double vDM =event.Velocity().norm();
				// double qMax2 = 4.0*pow(Reduced_Mass(mDM,mNucleus)*event.Velocity().norm(),2.0);
				if(IS_Angle<1e-4)
				{
					auto pdf = std::bind(PDF_ScatterAngle,_1,A,mDM,vDM,0.0);
					cosalpha = Rejection_Sampling(pdf,-1.0,1.0,0.0,pdf(1.0),PRNG);
					weight = 1.0;
				}
				else
				{
					if( PDF_ScatterAngle(-1.0,A,mDM,vDM,IS_Angle)<0.0)
					{
						auto pdf = std::bind(PDF_ScatterAngle,_1,A,mDM,vDM,0.0);
						cosalpha = Rejection_Sampling(pdf,-1.0,1.0,0.0,pdf(1.0),PRNG);
						weight = 1.0;
					}
					else
					{
						auto pdf = std::bind(PDF_ScatterAngle,_1,A,mDM,vDM,IS_Angle);
						cosalpha = Rejection_Sampling(pdf,-1.0,1.0,0.0,pdf(1.0),PRNG);
						weight = PDF_ScatterAngle(cosalpha,A,mDM,vDM,0.0)/PDF_ScatterAngle(cosalpha,A,mDM,vDM,IS_Angle);
					}
				}
			}
		//Construction of n, the unit vector pointing into the direction of vfinal.
			Eigen::Vector3d ep(ev(2),ev(2),-(ev(0)+ev(1)));
			ep.normalize();
			double alpha = acos(cosalpha);	
			if(alpha<M_PI/2) 	n=ev+tan(alpha)*ep;
			else 				n=-ev-tan(alpha)*ep;
			n.normalize();
			//We chose a particular ep, therefore we finally rotate n around ev by a random angle.
				double phi=PhiSample(PRNG);
				n=(n.dot(ev))*ev+cos(phi)*(ev.cross(n)).cross(ev)+sin(phi)*ev.cross(n);
		//Find post-scattering velocity, see Landau Lifshitz..
			Eigen::Vector3d vNew = mNucleus/(mDM+mNucleus)*vIni.norm()*n+mDM/(mDM+mNucleus)*vIni;
			event.SetVelocity(vNew);
	}

//Simulate one particle track:
	Trajectory ParticleTrack(double mDM,double sigma,Event& IniCondi, double vcut, std::mt19937& PRNG,bool ldm)
	{
		// We start outside the earth
			bool Underground=true;
		//Vector for the output:
			Trajectory output;
		//Initial Conditions, which are above the ground.
			output.New_Event(IniCondi);
		//Current Event
			Event Current_Event = IniCondi;
		//Nucleus probability:
			std::vector<std::vector<double>> NucleusProbability;
		//Point of Entry
			double tEntry=-Current_Event.Position()[2]/Current_Event.Velocity()[2];
			if(tEntry<0.0)
			{
				Underground=false;
				std::cout<<"Error: Invalid IC. Velocity points up." <<endl;
			}
			else
			{
				Current_Event.IncreaseTime(tEntry);
				Eigen::Vector3d dx = tEntry*Current_Event.Velocity();
				dx[2] = (-1.0)*Current_Event.Position()[2]; //for numerical precision reasons.
				Current_Event.IncreasePosition(dx);
				output.New_Event(Current_Event);
				NucleusProbability = Layers[ Current_Layer(Current_Event) - 1 ].Get_Probability();
			}
		//Underground Trajectory
			while(Underground)
			{
				if(!ldm) Compute_All_MFP(mDM,sigma,ldm,Current_Event.Velocity().norm());
				//New position
					double MFP_Weight;
					Free_Propagation(Current_Event,MFP_Weight,NucleusProbability,PRNG);
					output.Update_Weight(MFP_Weight);
				//Check if the particle is still underground
					if(Current_Layer(Current_Event) == 0)
					{
						Underground = false;
						//Point of exit:
							output.New_Event(Current_Event);
						//Final point above ground:
							
							double dt = 1*km/Current_Event.Velocity()[2];
							Eigen::Vector3d dx = dt*Current_Event.Velocity();
							Current_Event.IncreaseTime(dt);
							Current_Event.IncreasePosition(dx);
							output.New_Event(Current_Event);
					}
				//Check if the particle reached the detector
					else if (Current_Layer(Current_Event) == Detector_Index)
					{
						output.New_Event(Current_Event);
						Underground = false;
					}
				//Else scattering
					else
					{
						double Angle_Weight;
						Scatter(NucleusProbability,Current_Event,Angle_Weight,mDM,PRNG,ldm);
						output.Update_Weight(Angle_Weight);
						output.New_Event(Current_Event);
						//Check if we reached the cutoff speed.
							if(Current_Event.Velocity().norm()<vcut) 
							{
								Underground = false;
							}
					}
			}
		return output;
	}
//Simulate particles to gather velocity data
	std::vector<DataPoint> Simulate_Trajectories(int SampleSize,double mDM, double sigma,double vCutoff,double& lMax,long long unsigned int& ParticleCounter,std::mt19937& PRNG,bool ldm,int rank)
	{
		std::vector<DataPoint> OutputData;
		double ScatteringCounter=0.0;
		ParticleCounter=0;
		double Weight_Sum=0;
		double Weight2_Sum=0;
		int DataCounter=0;
		lMax=0;
		int k = 1; //for progress bar
		Compute_All_MFP(mDM,sigma,ldm);
		// Do the random walks and collect velocity data at detector depth.
			//Progress bar that disappears afterwards.
				if(rank==0)		cout <<"                                       |\r";
			while(DataCounter<SampleSize)
			{
				//Initial Conditions
					Eigen::Vector3d xIni(0.0,0.0,10.0*meter);
				//Simulate Trajectory
					Event IC = InitialCondition(0.0,xIni,PRNG,vCutoff); //WARNING: If vCutoff is put in here, vCutoff must be set in ComputeAttenuation in main()
					Trajectory trajectory=ParticleTrack(mDM,sigma,IC,vCutoff,PRNG,ldm);
					ParticleCounter++;
				//Check if the trajectory was successful.
					double speed=trajectory.TrajectoryEnd().Velocity().norm();
					if(trajectory.Trajectory_Type()==3&&speed>vCutoff)
					{
						//Check for the maximum horizontal distance
							// double l = trajectory.Maximum_Horizontal_Distance();
							// if(l>lMax) lMax=l;
						//Extract data at detection depth
							DataCounter++;
							//Progress bar
								if(rank==0&&1.0*DataCounter/SampleSize>k*0.025)
								{
									cout <<"="<<flush;
									k++;
								}
							double weight=trajectory.Get_Weight();
							Weight_Sum+=weight;
							Weight2_Sum+=weight*weight;
							ScatteringCounter+=weight*trajectory.NoOfScatterings();
							OutputData.push_back( DataPoint (speed,weight) );
					}

			}
		//Erase progress bar
			if(rank==0) cout <<"\r                                          \r";
		//Return data
			return OutputData;
	}