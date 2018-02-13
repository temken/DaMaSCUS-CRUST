#ifndef __Trajector_Simulation_hpp_
#define __Trajector_Simulation_hpp_

#include <Eigen/Geometry>
#include <functional>
#include <random>

#include "Physical_Parameters.hpp"
#include "General_Utilities.hpp"

#include "Layer_Class.hpp"
#include "RN_Generators.hpp"
#include "Trajectory_Class.hpp"

//Initial Condition Generator
	extern Event InitialCondition(double tIni,Eigen::Vector3d& xIni,std::mt19937& PRNG,double vCutoff = 0.0);

//Propagate a particle until it scatters, reaches the detector or gets back into space
	extern void Free_Propagation(Event& event,double& weight,std::vector<std::vector<double>>& prob,std::mt19937& PRNG,double logXi = -1.0);

//Find the scattering nucleus species:
	extern std::vector<double> ScatterNucleus(std::vector<std::vector<double>>& prob,std::mt19937& PRNG);

//Perform the scattering.
	extern double PDF_ScatterAngle(double cosalpha,double A,double mDM,double vDM,double deltaIS);
	extern void Scatter(std::vector<std::vector<double>>& prob,Event& event, double& weight,double mDM,std::mt19937& PRNG);

//Simulate one particle track:
	extern Trajectory ParticleTrack(double mDM,double sigman0,Event& IniCondi, double vcut, std::mt19937& PRNG,bool ldm);

//Generate velocity data
	extern std::vector<DataPoint> Simulate_Trajectories(int SampleSize,double mDM, double sigma,double vCutoff,double& hMax,long long unsigned int& ParticleCounter,std::mt19937& PRNG,bool ldm,int rank = 0);


#endif