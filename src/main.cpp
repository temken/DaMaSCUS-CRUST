#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <functional>
#include <cmath>
#include <algorithm>
#include <sstream>

#include "mpi.h"

#include "Physics_Functions.hpp"
#include "Input_Parameters.hpp"
#include "Trajectory_Simulation.hpp"
#include "DD_General.hpp"

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
	    //Ring Communication
	    	int source = (myRank==0)? numprocs - 1: myRank-1;
			int destination = (myRank==numprocs-1)? 0 : myRank+1;
			int tag =0;
	    //MPI status and requests.
	    	MPI_Status   status;
   			MPI_Request	send_request;
		//Datapoint mpi_type
			MPI_Datatype  mpi_datapoint;
			MPI_Type_contiguous(2, MPI_DOUBLE, &mpi_datapoint);
			MPI_Type_commit(&mpi_datapoint);
	//Starting time
	  	auto tStart = steady_clock::now();
	//Initialize Random Number Generator:
		std::random_device rd;
	  	std::mt19937 PRNG(rd());
	
	//Initial output
		if(myRank==0)
		{
			std::time_t start_date = std::chrono::system_clock::to_time_t(system_clock::now());
			cout <<"\n##############################"<<endl
			<<"DaMaSCUS - CRUST "<<version<<endl<<endl
			<<"Starting Time: " <<std::ctime(&start_date)<<endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	//Read in configuration file, defining the layer structure and analysis parameter.
		Read_Config_File(argv[1],myRank);

	//OPTION 1: Only compute the underground speed distribution for a given mass, cross section, and speed cutoff:
	if(argc > 2)
	{
		if(argc == 3 || argc > 6)
		{
			cerr <<"Error in main(): For PDF-only run, DaMaSCUS-CRUST requires between 2 and 4 additional command line arguments: mDM[GeV], sigma[cm^2] (optional: vMin[km/sec], SampleSize)."<<endl;
			std::exit(EXIT_FAILURE);
		}
		double mDM_pdf = atof( argv[2]) * GeV;
		double sigma_pdf = atof(argv[3]) * cm*cm;
		double vMin_pdf = (argc==5)? atof( argv[4]) * km/sec : 0;
		if(argc == 6) SampleSize = atof(argv[5]);
		DM.Set_Mass(mDM_pdf);
		DM.Set_Sigma_n(sigma_pdf);
		double Analytic_vMean =Average_Speed(vEarth,vMin_pdf);
		Compute_All_MFP(DM,vMin_pdf);
		int GIS_Domains=1;
		if(GIS)
		{
			GIS_Domains = std::min(3,Importance_Domains(DM,Analytic_vMean,vMin_pdf,GIS_Kappa,Layers));
			Compute_Importance_Boundaries(DM,Analytic_vMean,GIS_Domains,GIS_Splits,Layers);
		}
		//Pre-simulation output.
		if(myRank==0)
		{
			cout<<"\nCompute the DM speed distribution for\n\tsigma_n = "<<Round(InUnits(sigma_pdf,cm*cm)) <<"cm^2"<<endl
			<<"\tmDM = "<<((DM.mass<GeV)? Round(DM.mass/MeV) : Round(DM.mass)) <<(((DM.mass<GeV)? " MeV" :" GeV")) <<endl
			<<"\tvMin = " <<vMin_pdf/km*sec <<" km/sec"<<endl;
			if(GIS) cout <<"\tGIS domains:\t"<<GIS_Domains<<endl;
			else cout <<endl;
		}
		//Simulation outputs
		double lMax =0.0;
		long long unsigned int nData=0;
		long long unsigned int nData_local=0;
		long long unsigned int nData_new=0;
		long long unsigned int nParticles=0;
		long long unsigned int nReflections=0;
		double nScatteringsAv=0.0;
		std::vector<DataPoint> data_local;

		//Get the MPI ring communication going, and create a progress bar
		tag=0;
		int k = 0;
		int BarLength=50;
		if(myRank==0)
		{
			MPI_Isend(&nData,1,MPI_UNSIGNED_LONG_LONG,destination,tag,MPI_COMM_WORLD,&send_request);
			for(int j=0;j<BarLength;j++) cout <<" ";
			cout <<"|\r";						
		}
		while(nData < SampleSize)
		{
			//Simulate a DM particle.
				Event IC = InitialCondition(0.0,Vector3D(0.0,0.0,10.0*meter),PRNG,vMin_pdf);
				std::vector<Result> results= Simulate_Trajectory(DM,IC,vMin_pdf,PRNG);
				nParticles++;
			//Check if the DM particle made it to the detectory
				for(unsigned int i=0;i<results.size();i++)
				{
					if(results[i].success)
					{
						nData_local++;
						nData_new++;
						//Check for maximum horizontal displacement.
						if(lMax<results[i].max_horizontal_distance) lMax=results[i].max_horizontal_distance;
						
						nScatteringsAv=((nData_local-1.0)*nScatteringsAv+ results[i].nScattering)/nData_local;
						//Save the new data point
						data_local.push_back( DataPoint (results[i].final_event.Speed(),results[i].weight) );
					}
					else if(results[i].final_event.Position()[2]>0.0) nReflections++;
				}

			//Check if MPI token arrived.
			int flag;
			MPI_Iprobe(source,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
			if(flag)
			{
				MPI_Recv(&nData,1,MPI_UNSIGNED_LONG_LONG,source,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
				
				if(nData<SampleSize && (nData+nData_new)>=SampleSize) tag = source+1;
				else if (nData>=SampleSize) tag = status.MPI_TAG;
				nData+=nData_new;
				nData_new=0;
				//Pass on the token, unless you are the very last process.
				if(tag != (myRank+1) ) MPI_Isend(&nData,1,MPI_UNSIGNED_LONG_LONG,destination,tag,MPI_COMM_WORLD,&send_request);
				//Progress bar
				while(myRank==0&&std::min(1.0*nData/SampleSize,1.0)>1.0*k/BarLength)
				{
					cout <<"="<<flush;
					k++;
				}
			}
			
		}
		//Delete progress bar
		if(myRank==0)
		{
			cout <<"\r";
			for(int j=0;j<=BarLength;j++) cout <<" ";
			cout <<"\r";
		}
		//MPI reductions
		//Global variables
			long long unsigned int Global_ParticleCounter;
			long long unsigned int Global_ReflectionCounter;
			double Global_lMax;
			double Global_nSc;
			long long unsigned int Global_SampleSize;

		//Total number of data points
			MPI_Allreduce(&nData,&Global_SampleSize,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,MPI_COMM_WORLD);
		//Gather data at the master process
			unsigned long long int nData_locals[numprocs];
			int recv_counts[numprocs];
			int recv_displs[numprocs];
			MPI_Gather(&nData_local,1,MPI_UNSIGNED_LONG_LONG,&nData_locals,1,MPI_UNSIGNED_LONG_LONG,0,MPI_COMM_WORLD);
			for(int i=0;i<numprocs;i++) 
			{
				recv_counts[i] = nData_locals[i];
				recv_displs[i]= (i==0)? 0 : recv_displs[i-1] + recv_counts[i-1];
			}
			std::vector<DataPoint> Global_Data;
			if(myRank ==0) Global_Data.resize(Global_SampleSize);
			MPI_Gatherv(&data_local.front(),nData_local,mpi_datapoint,&Global_Data.front(),recv_counts,recv_displs,mpi_datapoint,0,MPI_COMM_WORLD);
		//Counter
			MPI_Allreduce(&nParticles,&Global_ParticleCounter,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
			MPI_Reduce(&nReflections,&Global_ReflectionCounter,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
		//Maximum distance
			MPI_Reduce(&lMax,&Global_lMax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
		//Average number of scatterings
			MPI_Reduce(&nScatteringsAv,&Global_nSc,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			Global_nSc/=numprocs;

		MPI_Barrier(MPI_COMM_WORLD);  	
		// The master process computes the event rates
		if(myRank == 0)
		{
			//0. Post simulation output:
				std::vector<double> vMean = Weighted_Average(Global_Data);
				cout 	<<"\tParticles: " <<Global_ParticleCounter
						<<"\tnData: " <<Global_SampleSize
						<<"\tData/Refl/Stop[\%] = "
						<<Round(100.0*Global_SampleSize/Global_ParticleCounter)<<" / "
						<<Round(100.0*Global_ReflectionCounter/Global_ParticleCounter)<<" / "
						<<Round(100.0*(Global_ParticleCounter-Global_ReflectionCounter-Global_SampleSize)/Global_ParticleCounter)
						<<endl<<"\t<v>=("<<Round(InUnits(vMean[0],km/sec))<<"+-"<<Round(InUnits(vMean[1],km/sec),2)<<")km/sec\t(" <<Round(InUnits(Analytic_vMean,km/sec))<<" km/sec)"
						<<"\t<nSC>="<<Round(Global_nSc)
						<<"\tlMax="<<Round(InUnits(Global_lMax,km))<<"km"<<endl<<endl;
			//1. MC result
				//1.1 Compute attenuation
				std::vector <double> attenuation = Compute_Attenuation(Global_ParticleCounter,Global_Data,vMin_pdf);
				cout<<"\tAttenuation: " <<attenuation[0] <<"+-" <<attenuation[1] <<endl;
				//1.2 Compute the PDF via KDE
				Interpolation kde = Perform_KDE(Global_Data,vMin_pdf,(vesc+vEarth));
				kde.Multiply(attenuation[0]);
				//1.3 Save the PDF to a file
				std::ostringstream streamObj;
				streamObj <<"../results/"<<ID <<"/DM_Speed_PDF_mDM="<< mDM_pdf<<"_GeV_sigma="<<sigma_pdf/cm/cm <<"_cm2_vMin="<<vMin_pdf /km*sec<<"_kps.txt";
				std::string filepath = streamObj.str();
				kde.Save_Function(filepath,200);
		}
	}
	//OPTION 2: Run a full parameter scan for a given experiment and return the constraints.
	else
	{

		//Find analytic limit curve:
			std::vector<std::vector<double>> AnalyticConstraints = Limit_Curve(DM,Detector,ID,mMin,mMax,Masses,CL,myRank);
		
		//Relative stopping powers
			DM.Set_Mass(GeV);
			DM.Set_Sigma_n(pb);
			double vDM = 1.0e-3;
			Compare_StoppingPower(DM,vDM,Layers,myRank);
		
		//Output files
			string filename1="../results/"+ID+"/Results_Detection.txt";
			string filename2="../results/"+ID+"/Limits_MC.txt";
			string filename3="../results/"+ID+"/Results_Simulation.txt";
			ofstream f1,f2,f3;
			if(myRank==0)
			{
				f1.open(filename1.c_str());
				f2.open(filename2.c_str());
				f3.open(filename3.c_str());
			}

		MPI_Barrier(MPI_COMM_WORLD);		
		////////////////////////////////////////////////////////////

		//Mass Scan
		int Mass_Scan = 1;	//not a boolean since it has to be MPI-broadcasted later
		bool Constraints = false;
		for(int i=0;i<Masses;i++)
		{
			if(Mass_Scan!=1) break;
			//Starting time of mass step
				auto tMass1 = steady_clock::now();
			//Set new mass
				DM.Set_Mass( pow(10.0,log10(mMin)+i*dm) );
			//The kinetic speed cutoff
				double vCutoff=Cutoff_Speed(Detector,DM.mass,myRank);
				if (vCutoff < 0) continue;
				double Analytic_vMean =Average_Speed(vEarth,vCutoff);
			//Analytic and MC bound.
				double sigma_Low_Analytic = AnalyticConstraints[i][1];
				double sigma_Low_MC=0.0;
			//Start of the cross section scan
				double sigmaMin=pow(10.0,floor(log10(sigma_Low_Analytic/cm/cm))-1.0)*cm*cm;
				if(myRank==0)
				{
					cout 	<<endl<<(i+1)<<"/"<<Masses<<")\nmDM = " <<((DM.mass<GeV)? Round(DM.mass/MeV) : Round(DM.mass)) <<(((DM.mass<GeV)? " MeV" :" GeV"))<<"\tvCut = " <<Round(InUnits(vCutoff,km/sec)) <<" km/sec"<<endl
							<<"\tAnalytic bound: "<<Round(InUnits(sigma_Low_Analytic,cm*cm))<<"cm^2"<<endl;
				}
			//Set initial cross section step size.
				double dsigma=dSigma;
			//CS Scan
				std::vector< std::vector<double>> interpolationlist;
				double Last_sigma = 0.0;
				double Last_Ntot=1e-10;
				double Last_LLH=2.0;
				bool SignalDecrease = false;
				bool FoundConstraint = false;
				bool LLH_Increase = false;
				int Decrease_Counter =3;
				int step_counter=0;
				int CS_Scan=1; //not a boolean since it has to be MPI-broadcasted later
				int GIS_Domains=1;

				double sigma = sigmaMin;
				while(CS_Scan==1)
				{
					//Timing of a given sigma simulation
		  				auto tSigma1 = steady_clock::now();
						step_counter++;
					//Find nuclear cross section
						DM.Set_Sigma_n( Stopping_CrossSection(Detector,sigma,DM.mass) );
					//Compute the next the mean free path.
						Compute_All_MFP(DM,vCutoff);
					//Compute GIS domains if activated;
						if(GIS)
						{
							int GIS_Domains_new = Importance_Domains(DM,Analytic_vMean,vCutoff,GIS_Kappa,Layers);
							GIS_Domains = std::min(25,std::min(GIS_Domains_new,GIS_Domains+1));
							Compute_Importance_Boundaries(DM,Analytic_vMean,GIS_Domains,GIS_Splits,Layers);
						}
					//Pre-simulation output.
						if(myRank==0)
						{
								cout<<"\t"<<step_counter<<")\t"<<Round(InUnits(sigma,cm*cm)) <<"cm^2";
								if(sigma!=DM.sigma_n) cout <<"\t("<<Round(InUnits(DM.sigma_n,cm*cm))<<"cm^2)";
								cout<<"\t(mDM = "<<((DM.mass<GeV)? Round(DM.mass/MeV) : Round(DM.mass)) <<(((DM.mass<GeV)? " MeV)" :" GeV)"));
								if(GIS) cout <<"\tGIS domains:\t"<<GIS_Domains<<endl;
								else cout <<endl;
						}
					//Simulation outputs
						double lMax =0.0;
						long long unsigned int nData=0;
						long long unsigned int nData_local=0;
						long long unsigned int nData_new=0;
						long long unsigned int nParticles=0;
						long long unsigned int nReflections=0;
						double nScatteringsAv=0.0;
						std::vector<DataPoint> data_local;

						//Get the MPI ring communication going, and create a progress bar
						tag=0;
						int k = 0;
						int BarLength=50;
						if(myRank==0)
						{
							MPI_Isend(&nData,1,MPI_UNSIGNED_LONG_LONG,destination,tag,MPI_COMM_WORLD,&send_request);
							for(int j=0;j<BarLength;j++) cout <<" ";
							cout <<"|\r";						
						}
						while(nData < SampleSize)
						{
							//Simulate a DM particle.
								Event IC = InitialCondition(0.0,Vector3D(0.0,0.0,10.0*meter),PRNG,vCutoff);
								std::vector<Result> results= Simulate_Trajectory(DM,IC,vCutoff,PRNG);
								nParticles++;
							//Check if the DM particle made it to the detectory
								for(unsigned int i=0;i<results.size();i++)
								{
									if(results[i].success)
									{
										nData_local++;
										nData_new++;
										//Check for maximum horizontal displacement.
										if(lMax<results[i].max_horizontal_distance) lMax=results[i].max_horizontal_distance;
										
										nScatteringsAv=((nData_local-1.0)*nScatteringsAv+ results[i].nScattering)/nData_local;
										//Save the new data point
										data_local.push_back( DataPoint (results[i].final_event.Speed(),results[i].weight) );
									}
									else if(results[i].final_event.Position()[2]>0.0) nReflections++;
								}

							//Check if MPI token arrived.
							int flag;
							MPI_Iprobe(source,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
							if(flag)
							{
								MPI_Recv(&nData,1,MPI_UNSIGNED_LONG_LONG,source,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
								
								if(nData<SampleSize && (nData+nData_new)>=SampleSize) tag = source+1;
								else if (nData>=SampleSize) tag = status.MPI_TAG;
								nData+=nData_new;
								nData_new=0;
								//Pass on the token, unless you are the very last process.
								if(tag != (myRank+1) ) MPI_Isend(&nData,1,MPI_UNSIGNED_LONG_LONG,destination,tag,MPI_COMM_WORLD,&send_request);
								//Progress bar
								while(myRank==0&&std::min(1.0*nData/SampleSize,1.0)>1.0*k/BarLength)
								{
									cout <<"="<<flush;
									k++;
								}
							}
							
						}
						//Delete progress bar
						if(myRank==0)
						{
							cout <<"\r";
							for(int j=0;j<=BarLength;j++) cout <<" ";
							cout <<"\r";
						}
						//MPI reductions
						//Global variables
							long long unsigned int Global_ParticleCounter;
							long long unsigned int Global_ReflectionCounter;
							double Global_lMax;
							double Global_nSc;
							long long unsigned int Global_SampleSize;

						//Total number of data points
							MPI_Allreduce(&nData,&Global_SampleSize,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,MPI_COMM_WORLD);
						//Gather data at the master process
							unsigned long long int nData_locals[numprocs];
							int recv_counts[numprocs];
							int recv_displs[numprocs];
							MPI_Gather(&nData_local,1,MPI_UNSIGNED_LONG_LONG,&nData_locals,1,MPI_UNSIGNED_LONG_LONG,0,MPI_COMM_WORLD);
							for(int i=0;i<numprocs;i++) 
							{
								recv_counts[i] = nData_locals[i];
								recv_displs[i]= (i==0)? 0 : recv_displs[i-1] + recv_counts[i-1];
							}
							std::vector<DataPoint> Global_Data;
							if(myRank ==0) Global_Data.resize(Global_SampleSize);
							MPI_Gatherv(&data_local.front(),nData_local,mpi_datapoint,&Global_Data.front(),recv_counts,recv_displs,mpi_datapoint,0,MPI_COMM_WORLD);
						//Counter
							MPI_Allreduce(&nParticles,&Global_ParticleCounter,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
							MPI_Reduce(&nReflections,&Global_ReflectionCounter,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
						//Maximum distance
							MPI_Reduce(&lMax,&Global_lMax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
						//Average number of scatterings
							MPI_Reduce(&nScatteringsAv,&Global_nSc,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
							Global_nSc/=numprocs;

					MPI_Barrier(MPI_COMM_WORLD);  	
					// The master process computes the event rates
					if(myRank == 0)
					{
						//0. Post simulation output:
							std::vector<double> vMean = Weighted_Average(Global_Data);
							cout 	<<"\t\tParticles: " <<Global_ParticleCounter
									<<"\tnData: " <<Global_SampleSize
									<<"\tData/Refl/Stop[\%] = "
									<<Round(100.0*Global_SampleSize/Global_ParticleCounter)<<" / "
									<<Round(100.0*Global_ReflectionCounter/Global_ParticleCounter)<<" / "
									<<Round(100.0*(Global_ParticleCounter-Global_ReflectionCounter-Global_SampleSize)/Global_ParticleCounter)
									<<endl<<"\t\t<v>=("<<Round(InUnits(vMean[0],km/sec))<<"+-"<<Round(InUnits(vMean[1],km/sec),2)<<")km/sec\t(" <<Round(InUnits(Analytic_vMean,km/sec))<<" km/sec)"
									<<"\t<nSC>="<<Round(Global_nSc)
									<<"\tlMax="<<Round(InUnits(Global_lMax,km))<<"km"<<endl;
						//1. Analytic result
							double Ntot=0.0;
							std::function<double(double)> spectrum = Compute_Spectrum(DM,Detector);
							double LLH = Likelihood(DM,Detector,spectrum,Ntot);
						//2. MC result
							//2.1 Compute attenuation
								std::vector <double> attenuation = Compute_Attenuation(Global_ParticleCounter,Global_Data,vCutoff);
							//2.2 Likelihood and total number of events.
								std::vector<double> Ntot_MC;
								Interpolation spectrum_MC = Compute_Spectrum(DM,Detector,Global_Data,attenuation);
								std::vector<double> LLH_MC = Likelihood(DM,Detector,spectrum_MC,attenuation,Global_Data,vCutoff,Ntot_MC);
							//2.3 Create some output
								f1 <<DM.mass <<"\t" <<InUnits(sigma,cm*cm)<<"\t"<<LLH_MC[0] <<"\t"<<LLH_MC[1] <<"\t"<<LLH <<"\t"<<Ntot_MC[0]<<"\t"<<Ntot_MC[1]<<"\t"<<Ntot<<endl;
								f3 <<DM.mass <<"\t" <<InUnits(sigma,cm*cm)<<"\t"<<Global_nSc<<"\t"<<Global_ParticleCounter<<"\t"<<Global_ReflectionCounter<<"\t"<<(Global_ParticleCounter-Global_ReflectionCounter-Global_SampleSize)<<"\t"<<attenuation[0]<<"\t"<<attenuation[1]<<endl;
							//2.4 append likelihood to list
								if(LLH_MC[0]>0.0) interpolationlist.push_back(std::vector<double> {log10(sigma),log10(LLH_MC[0])});
								else interpolationlist.push_back(std::vector<double> {log10(sigma),-30.0});
							//2.5 Post-analysis console output
								auto tSigma2 = steady_clock::now();
								double SigmaDuration =1e-6*duration_cast<microseconds>( tSigma2 - tSigma1 ).count();
								cout <<"\t\tN = "<<Round(Ntot_MC[0],3)<<"+-"<<Round(Ntot_MC[1],2)<<"\t("<<Round(Ntot,3)<<")\tL="<<Round(LLH_MC[0],2) <<"+-"<<Round(LLH_MC[1],2) <<" ("<<Round(LLH,2)<<")\t("<<Round(SigmaDuration) <<"sec)"<<endl<<endl;
						//3. Find the next cross-section.
							//Check if the signal starts to decrease.
								if(step_counter==1 && LLH_MC[0]<(1-CL))
								{
									cout <<"Note: LLH_MC < (1-CL) on the first step. No constraint possible for mDM = " <<DM.mass<<" GeV." <<endl;
									CS_Scan=0;
								}
								else if(!FoundConstraint && LLH_MC[0]<(1-CL))
								{
									FoundConstraint=true;
									cout <<"\tNote: LLH_MC <(1-CL)."<<flush;
										//Interpolate
											Interpolation logL(interpolationlist);
										//Create function for root finding
											std::function<double(double)> func = [&logL] (double logSigma)
											{
												return ( logL(logSigma)-log10(1.0-CL) );
											};
										//Find the critical cross-section
											double logSigma = Find_Root(func,log10(Last_sigma),log10(sigma),1e-6);
											sigma_Low_MC = pow(10.0,logSigma);
										//Terminal output
											cout <<"\tUpper bound:\t"<<((DM.mass<GeV)? Round(DM.mass/MeV) : Round(DM.mass)) <<(((DM.mass<GeV)? " MeV\t" :" GeV\t"))<<Round(InUnits(sigma_Low_MC,cm*cm))<<" cm^2 ("<<Round( InUnits(sigma_Low_Analytic,cm*cm))<<" cm^2)"<<endl<<endl;
											if(sigma_Low_MC>=1.5*sigma_Low_Analytic) cout <<"\tNote: MC bound considerably higher."<<endl<<endl;

								}
								if(!SignalDecrease&&Ntot_MC[0]<Last_Ntot) 
								{
									cout <<"\tNote: The detection signal starts to decrease. Decrease CS step size." <<endl<<endl;
									SignalDecrease = true;
									dsigma/=5.0;
								}
								else if(Ntot_MC[0]/Last_Ntot<0.1 && Decrease_Counter>0)
								{
									cout <<"\tNote: Signal decreases fast. Decrease CS step size further." <<endl<<endl;
									dsigma/=2.0;
									Decrease_Counter--;
								}
							//Check if likelihood starts increasing again.
								else if(!LLH_Increase&&LLH_MC[0]/Last_LLH>1.2)
								{
									cout <<"\tNote: Likelihood starts growing again. Decrease CS step size."<<endl<<endl;
									dsigma/=2.0;
									LLH_Increase=true;
								}
							//Check if we are done for this mass
								if((SignalDecrease||LLH_Increase)&& LLH_MC[0]>(1.0-CL))
								{
									if(FoundConstraint)
									{
										cout <<"\tNote: Likelihood grows above 1-CL."<<endl;
										//Interpolate
											Interpolation logL(interpolationlist);
											
											
										//Create function for root finding
											std::function<double(double)> func = [&logL] (double logSigma)
											{
												return ( logL(logSigma)-log10(1.0-CL) );
											};
										//Find the critical cross-section via regula falsi
											double logSigma = Find_Root(func,log10(Last_sigma),log10(sigma),1e-6);
											sigma = pow(10.0,logSigma);
											double sigma_Low;
											if(sigma_Low_MC>=1.5*sigma_Low_Analytic) sigma_Low = sigma_Low_MC;
											else sigma_Low = sigma_Low_Analytic;
										//Save limits
											f2 <<DM.mass <<"\t" <<InUnits(sigma_Low,cm*cm)<<"\t" <<InUnits(sigma,cm*cm)<<endl;
										//Terminal output
											auto tMass2 = steady_clock::now();
											double MassDuration =1e-6*duration_cast<microseconds>( tMass2 - tMass1 ).count();
											cout <<"\tFound limit:\t"<<((DM.mass<GeV)? Round(DM.mass/MeV) : Round(DM.mass)) <<(((DM.mass<GeV)? " MeV\t" :" GeV\t"))<<Round(InUnits(sigma_Low,cm*cm))<<" cm^2 - "<<Round(InUnits(sigma,cm*cm)) <<" cm^2\t("<<Round(MassDuration)<<"sec)"<<endl<<endl;
										
										//End the cross-section scan
											CS_Scan=0;
											Constraints=true;
									}
									else
									{
										cout<<"Likelihood grows above 1-CL, but no constraint could be imposed for mDM = " <<DM.mass<<" GeV.";
										//End the cross-section scan
										CS_Scan=0;
										//End the mass scan, if the constraints ended.
										if(Constraints)
										{
											Mass_Scan=0;
											cout <<"\tMass scan is aborted."<<endl;
										}
										else cout <<endl;

									}
									
								}
							//Otherwise find the next cross-section
								else
								{
									Last_LLH = LLH_MC[0];
									Last_sigma =sigma;
									Last_Ntot = Ntot_MC[0];
									sigma = pow(10.0,(log10(sigma)+dsigma));
								}
					}
					//Broadcast the new cross section and parameter scan indicator
						MPI_Bcast(&sigma,1,MPI_DOUBLE,0,MPI_COMM_WORLD);								
						MPI_Bcast(&CS_Scan,1,MPI_INT,0,MPI_COMM_WORLD);	
						MPI_Bcast(&Mass_Scan,1,MPI_INT,0,MPI_COMM_WORLD);							
				}
		}	
		if(myRank==0)
		{
			f1.close();
			f2.close();
		}
		MPI_Barrier(MPI_COMM_WORLD);  	
	}
///////////////////////////////////////////////////////////////////
	//Ending time and computing time
		auto tEnd = steady_clock::now();
		double duration =1e-6*duration_cast<microseconds>( tEnd - tStart ).count();
		if(myRank==0)cout <<"\nProcessing Time: "<< duration <<"s\a"<<endl;
	// Finalize the MPI environment.
	    MPI_Finalize();
    
    return 0;
}
