#include "RN_Generators.hpp"

#include <cmath>
#include <iostream>

#include "Physical_Parameters.hpp"


//Random Number Generation:
	// Real random number which is used as a seed for the PRNG:
	std::random_device rd;
	//The PRNG:
  	std::mt19937 generator(rd());
  	//Propability:
  	std::uniform_real_distribution<double> distrXi(0,1);
  	//Uniform angles:
  	std::uniform_real_distribution<double> distrcos(-1,1);
  	std::uniform_real_distribution<double> distrphi(0.0,2*M_PI);
  	//For Rejection sampling, to find the norm of the maxwell distributed initial velocity
  	std::uniform_real_distribution<double> distrx(0.0,vesc);
	std::uniform_real_distribution<double> distry(0.0,1/Nesc);

//Random numbers between 0 and 1:
	double ProbabilitySample(std::mt19937& PRNG)
	{
		return distrXi(PRNG);
	}

//MaxwellSample after integration over the directions
	double MaxwellPDF(double v)
	{
		//return sqrt(2.0/M_PI)*pow(sigma,-3.0)*pow(x,2.0)*exp(-pow(x/sigma,2.0)/2);
		return 4*M_PI/Nesc*v*v*exp(-v*v/v0/v0);
	}
//Finding a sample from a maxwell distribution with sigma=v0/sqrt2. Method: Rejection sampling
	double MaxwellSample(std::mt19937& PRNG){
		bool success=false;
		double x,y;
		while(success==false){
			x=distrx(PRNG);
			y=distrXi(PRNG)*ymax;

			if(y<=MaxwellPDF(x))
			{
				success=true;
			}
		}
		return x;
	}

//Random Angles
	double ThetaSample(std::mt19937& PRNG){
	  	return acos(distrcos(PRNG));
	}
	double PhiSample(std::mt19937& PRNG){
	  	return distrphi(PRNG);
	}

//Rejection sampling
	double Rejection_Sampling(std::function<double(double)> PDF,double xMin,double xMax,double yMin,double yMax,std::mt19937& PRNG)
	{
		bool success=false;
		double x,y;
		while(!success)
		{
			x=xMin + distrXi(PRNG) * (xMax-xMin);
			y=yMin + distrXi(PRNG) * (yMax-yMin);
			double pdf = PDF(x);
			if (pdf < 0.0) 
			{
				cout <<"Error in Rejection_Sampling(): f("<<x<<")=" <<pdf<<endl;
			}
			else if(y<=PDF(x)) success=true;
		}
		return x;
	}


