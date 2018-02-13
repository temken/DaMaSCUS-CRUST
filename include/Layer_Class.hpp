#ifndef __Layer_Class_hpp_
#define __Layer_Class_hpp_

#include <vector>
#include <string>

#include "Physical_Parameters.hpp"
#include "Trajectory_Class.hpp"


class Layer
{
	std::string name;
	unsigned int index;
	double density,thickness,depth;
	std::vector<std::vector<double>> composition;
	
	double mfp=-1.0;
	std::vector<std::vector<double>> probability;

	public:
		//Constructor
			Layer(std::string& ID,unsigned int ind, double rho, double thick,double deep,std::vector<std::vector<double>>& comp);
		//Functions
			//Set and changes variables:
				void Set_Name(std::string newname);
				void Set_Density(double newdensity);
				void Set_Thickness(double newthickness);
				void Set_Composition(std::vector<std::vector<double>>& newcomposition);
			//Return values
				std::string Get_Name();
				double Get_Density();
				double Get_Thickness();
				double Get_Depth();
				std::vector<std::vector<double>> Get_Composition();	
				double Get_MFP();
				std::vector<std::vector<double>> Get_Probability();
			//Further Functions
				void Compute_MFP(double mDM, double sigma,double speed,bool ldm);
				void Print_Summary(); 

};

extern void Compute_All_MFP(double mDM,double sigma,bool ldm, double speed = 300*km/sec);
extern int Current_Layer(Eigen::Vector3d& x);
extern int Current_Layer(Event& event);
extern Event Leave_Layer(Event& x);

#endif