//Contains the actual trajectory simulation function.
#ifndef __Trajector_Simulation_hpp_
#define __Trajector_Simulation_hpp_

#include <functional>
#include <random>

#include "Physics_Functions.hpp"
#include "Simulation_Essentials.hpp"

//Cut off speed for a given experiment
	extern double Cutoff_Speed(const std::string& experiment,double mDM,int rank=0);

//Initial Condition Generator
	extern Event InitialCondition(double tIni,Vector3D xIni,std::mt19937& PRNG,double vCutoff);

//Simulate one particle track:
	struct Result
	{
		bool success;
		Event final_event;
		long unsigned int nScattering;
		double max_horizontal_distance;
		double weight;

		//Constructors
		Result();
		Result(double w);

		void Summary();
	};
	extern std::vector<Result> Simulate_Trajectory(DM_Particle& DM, Event IC, double vMin,std::mt19937& PRNG,Result result=Result(),std::string split_ID ="Contact");
	
#endif