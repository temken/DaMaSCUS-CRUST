#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <functional>
#include <cmath>

#include "mpi.h"

#include "Physical_Parameters.hpp"
#include "General_Utilities.hpp"
#include "Trajectory_Class.hpp"
#include "Trajectory_Simulation.hpp"
#include "Direct_Detection.hpp"
#include "Density_Estimation.hpp"

using namespace std;
using namespace std::chrono;
using namespace std::placeholders;

int main(int argc, char *argv[])
{
//INITIALIZATION	
////////////////////////////////////////////////////////////
	//MPI Enviroment
		// Initialize the MPI environment
	    	MPI_Init(NULL, NULL);
	    // Get the number of processes
	    	int numprocs;
	    	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	    // Get the ID number of the process
	    	int myRank;
	    	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	//Starting time
		high_resolution_clock::time_point tStart,tEnd;
	  	tStart = high_resolution_clock::now();
	//Initialize Random Number Generator:
		std::random_device rd;
	  	std::mt19937 PRNG(rd());
	//Initial output
		if(myRank==0)
		{
			std::time_t start_time = std::chrono::system_clock::to_time_t(tStart);
			cout <<"\n##############################"<<endl
			<<"DaMaSCUS - CRUST "<<version<<endl<<endl
			<<"Starting Time: " <<std::ctime(&start_time)<<endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	//Read in configuration file, defining the layer structure and analysis parameter.
		Read_Config_File(argv[1],myRank);
	//Find the local samplesize for each MPI process
		int Local_SampleSize = ceil(1.0*SampleSize/numprocs);
		int Global_SampleSize = numprocs*Local_SampleSize;


	//Cross-section scan step width
		double sigmaMax=1e-24*cm*cm;
		double MaxSigmaSteps=25;
	//Output files
		string filename1="../results/"+ID+"/"+Detector+"_Results.txt";
		string filename2="../results/"+ID+"/"+Detector+"_Limits.txt";
		ofstream f1,f2;
		if(myRank==0)
		{
			f1.open(filename1.c_str());
			f2.open(filename2.c_str());
		}		
		

	//Mass Scan
	for(int i=0;i<Masses;i++)
	{
		//Set new mass
			double mDM=pow(10.0,log10(mMin)+i*dm);
		//The kinetic speed cutoff for DAMIC
			vCutoff=Cutoff_Speed(Detector,mDM,myRank);
			if (vCutoff < 0) continue;
		
		if(myRank==0) cout <<endl<<"mDM = " <<mDM <<" GeV\t vCut = " <<InUnits(vCutoff,km/sec) <<" km/sec"<<endl;
		//Compute standard constraint
			double sigma_bound;
			if(myRank==0) 
			{
				sigma_bound = Upper_Bound(Detector,mDM,LDM,CL);
				cout <<"Upper bound: "<<InUnits(sigma_bound,cm*cm)<<"cm^2"<<endl;
			}
			MPI_Bcast(&sigma_bound,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		
		//Starting value of sigma below the usual bound.
			double sigmaMin=pow(10.0,floor(log10(sigma_bound/cm/cm))-1.0)*cm*cm;
		//Set rough cross section step size.
			double dsigma = (log10(sigmaMax)-log10(sigmaMin))/(MaxSigmaSteps-1);
		//CS Scan
			std::vector< std::vector<double>> interpolationlist;
			double sigma = sigmaMin;
			double Last_sigma = 0.0;
			double Last_Ntot=1e-10;
			bool SignalDecrease = false;
			bool LLH_Increase = false;
			int Decrease_Counter =3;
			// bool Intersection_Method = false;
			while(sigma<sigmaMax)
			{
				//Compute the next the mean free path.
					Compute_All_MFP(mDM,sigma,LDM);
				//Simulation outputs
					double Local_lMax=0.0;
					long long unsigned int Local_ParticleCounter=0;
				//Simulation
					std::vector<DataPoint> Local_Data=Simulate_Trajectories(Local_SampleSize,mDM,sigma,vCutoff,Local_lMax,Local_ParticleCounter,PRNG,LDM,myRank);
					std::vector<DataPoint> Global_Data;

				//Send the data to the master process 0:
					if(myRank ==0) Global_Data.resize(Global_SampleSize);
					double Global_lMax=0.0;
					long long unsigned int Global_ParticleCounter=0;
					//Data
						MPI_Datatype  mpi_datapoint;
						MPI_Type_contiguous(2, MPI_DOUBLE, &mpi_datapoint);
						MPI_Type_commit(&mpi_datapoint);
						MPI_Gather(&Local_Data.front(),Local_SampleSize,mpi_datapoint,&Global_Data.front(),Local_SampleSize,mpi_datapoint,0,MPI_COMM_WORLD);
					//Counter
						MPI_Reduce(&Local_ParticleCounter,&Global_ParticleCounter,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
					//Maximum distance
						MPI_Reduce(&Local_lMax,&Global_lMax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

				//The master process computes the event rates
				if(myRank == 0)
				{
					//1. Analytic result
						double Ntot;
						std::function<double(double)> spectrum = Compute_Spectrum(Detector,mDM,sigma,LDM);

						double LLH = Likelihood(Detector,spectrum,Ntot);
					//2. MC result
						//2.1 Compute attenuation
							std::vector <double> attenuation = Compute_Attenuation(Global_ParticleCounter,Global_Data,vCutoff);//WARNING: If vCutoff is put in here, vCutoff must be set in InitialConditions in SimulateTrajectory() as well.
						//2.2 Likelihood and total number of events.
							std::vector<double> Ntot_MC;
							// cout <<"\t Spectrum" <<endl;
							Interpolation spectrum_MC = Compute_Spectrum(Detector,mDM,sigma,LDM,Global_Data,attenuation);
							// cout <<"\tLLH" <<endl;
							std::vector<double> LLH_MC = Likelihood(Detector,spectrum_MC,attenuation,Ntot_MC);
						//2.3 Create some output
							f1 <<mDM <<"\t" <<InUnits(sigma,cm*cm)<<"\t"<<LLH_MC[0] <<"\t"<<LLH_MC[1] <<"\t"<<LLH <<"\t"<<Ntot_MC[0]<<"\t"<<Ntot_MC[1]<<"\t"<<Ntot<<"\t"<<attenuation[0]<<"\t"<<attenuation[1]<<endl;
							cout <<mDM<<"\t" <<InUnits(sigma,cm*cm) <<"cm2\t N = "<<Ntot_MC[0]<<"+-"<<Ntot_MC[1]<<" ("<<Ntot<<")\t\t\tL="<<LLH_MC[0] <<"+-"<<LLH_MC[1] <<" ("<<LLH<<")"<<endl;
						//2.4 append likelihood to list
							if(LLH_MC[0]>0.0) interpolationlist.push_back(std::vector<double> {log10(sigma),log10(LLH_MC[0])});
							else interpolationlist.push_back(std::vector<double> {log10(sigma),-50.0});
					//3. Find the next cross-section.
						//Check if the signal starts to decrease.
							if(!SignalDecrease&&Ntot_MC[0]<Last_Ntot) 
							{
								cout <<"Note: The detection signal starts to decrease. Decrease CS step size." <<endl;
								SignalDecrease = true;
								dsigma/=10.0;
							}
							else if(Ntot_MC[0]/Last_Ntot<0.3 && Decrease_Counter>0)
							{
								cout <<"Note: Signal decreases fast. Decrease CS step size further." <<endl;
								dsigma/=2.0;
								Decrease_Counter--;
							}
						//Check if we went below N_Limit
							if(SignalDecrease && LLH_MC[0] >0.0 && !LLH_Increase )
							{
								cout <<"Likelihood starts growing again.Decrease CS step size."<<endl;
								dsigma/=2.0;
								LLH_Increase = true;
							}
							if(SignalDecrease&& LLH_MC[0]>(1.0-CL))
							{
								cout <<"Likelihood grows above 1-CL. Use interpolation to find critical value."<<endl;
								//Interpolate
									Interpolation logL(interpolationlist);
									
									
								//Create function for root finding
									std::function<double(double)> func = [&logL] (double logSigma)
									{
										return ( logL(logSigma)-log10(1.0-CL) );
									};
									ofstream ft;
									ft.open("interpolationlist.txt");
									for(unsigned int i = 0; i<interpolationlist.size(); i++)
									{
										ft <<interpolationlist[i][0]<<"\t"<<interpolationlist[i][1]<<"\t"<<func(interpolationlist[i][0])<<endl;
									}
									ft.close();
								//Find the critical cross-section via regula falsi
									double logSigma = RegulaFalsi(func,log10(Last_sigma),log10(sigma),1e-6);
									sigma = pow(10.0,logSigma);
								//Save limits
									f2 <<mDM <<"\t" <<InUnits(sigma_bound,cm*cm)<<"\t" <<InUnits(sigma,cm*cm)<<endl;
								//Terminal output
									cout <<"Found limit:\t"<<mDM<<"\t"<<InUnits(sigma,cm*cm)<<endl<<endl;
								
								//End the cross-section scan
									sigma=cm*cm;
							}
						//Otherwise find the next cross-section
							else
							{
								Last_sigma =sigma;
								sigma = pow(10.0,(log10(sigma)+dsigma));
								Last_Ntot = Ntot_MC[0];
							}
				}
				//Share the new cross-section among the peasants
					MPI_Bcast(&sigma,1,MPI_DOUBLE,0,MPI_COMM_WORLD);								
			}
	}	


	if(myRank==0)
	{
		f1.close();
		f2.close();
	}
	MPI_Barrier(MPI_COMM_WORLD);  	
				
	//Ending time and computing time
		tEnd = high_resolution_clock::now();
		double duration =1e-6*duration_cast<microseconds>( tEnd - tStart ).count();
		if(myRank==0)cout <<"\nProcessing Time: "<< duration <<"s"<<endl;
	// Finalize the MPI environment.
	    MPI_Finalize();
    	return 0;
}
