//Contains the global input parameters, which are set by the configuration file.
#ifndef __Input_Parameters_hpp_
#define __Input_Parameters_hpp_

#include <vector>
#include <string>

#include "Physics_Functions.hpp"
#include "DM_Particle.hpp"
#include "Simulation_Essentials.hpp"

//1. Global variables and input parameters
	//Version
	extern std::string version;

	//Simulation ID
	extern std::string ID;

	//Statistical parameter
	extern unsigned int SampleSize;
	//Certainty Level
	extern double CL;
	//Importance sampling:
	extern double IS_Angle;	//Scattering Angle
	extern double IS_MFP;	//MFP
	//Importance splitting
	extern bool GIS;
	extern double GIS_Splits;
	extern double GIS_Kappa;
	//Parameter scan
	extern double mMin,mMax,dm;
	extern int Masses;
	extern double dSigma;

	//Layer structure
	extern bool Atmosphere ;
	extern double Altitude;
	extern int Atmo_Layers;
	extern std::vector<Layer> Layers;

	//Dark Matter Interaction
	extern DM_Particle DM;

	//Detector
	extern std::string Detector;

	//DM electron scattering experiments with semiconductors
	extern std::string DMe_target;
	extern double DMe_exposure;
	extern int DMe_threshold;
	extern double DMe_efficiency;
	extern unsigned int DMe_events;

	//Detector Location
	extern double Detector_Depth;
	extern int Detector_Index;

//2. Read in input parameters from configuration file.
	extern void Read_Config_File(const char *inputfile,int rank=0);

#endif