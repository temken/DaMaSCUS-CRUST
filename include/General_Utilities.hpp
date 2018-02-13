#ifndef __General_Utilities_hpp_
#define __General_Utilities_hpp_

#include <Eigen/Geometry>
#include <string>
#include <vector>
#include <functional>

#include "mpi.h"

#include "Layer_Class.hpp"

//Global variables
	//Version
		extern std::string version;
	//Simulation ID
		extern std::string ID;
	//Statistical parameter
		extern int SampleSize;
		//Certainty Level
			extern double CL;
		//Importance sampling:
			extern double IS_Angle;	//Scattering Angle
			extern double IS_MFP;	//MFP
		//Parameter scan
			extern double mMin,mMax,dm;
			extern int Masses;
	//Layer structure
		extern std::vector<Layer> Layers;
	//Dark Matter
		extern bool LDM;
	//Detector
		extern std::string Detector;

	//Derived variables
		//Speed Cutoff
			extern double vCutoff;
		//Detector Location
			extern double Detector_Depth;
			extern int Detector_Index;

//Read in cfg file
		extern void Read_Config_File(const char *inputfile,int rank=0);

//Read in table into vector
		extern std::vector<double> Read_List(std::string filepath,double dimension=1.0);

//Work division
	extern std::vector<int> WorkDivision(int WorldSize,int MassSteps);

//Numerical Integration
	extern double Find_Epsilon(std::function<double(double)> func, double a,double b,double precision);
	extern double AdaptiveSimpsons(std::function<double(double)> func, double a,double b, double epsilon,int maxRecursionDepth=20);
	double auxAdaptiveSimpson(std::function<double(double)> func, double a,double b, double epsilon,double S,double fa,double fb,double fc,int bottom,bool &warning);
//Root finding
	extern double RegulaFalsi(std::function<double(double)> func,double xLeft, double xRight,double epsilon);
//Sign function
	extern double Sign(double arg);

//Interpolation
	//Interpolation of tabulated data using Steffen splines 
	//(http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1990A%26A...239..443S&amp;data_type=PDF_HIGH&amp;whole_paper=YES&amp;type=PRINTER&amp;filetype=.pdf)
	class Interpolation 
	{
	private:
		std::vector<std::vector<double>> TabulatedData;
		unsigned int N_Data;
		std::vector<double> xDomain;
		//Steffen coefficients
		std::vector<double> a,b,c,d;
		//Pre-factor
		double preFactor;
		//compute Steffen coefficients
		void Compute_Steffen_Coefficients(std::vector<std::vector<double>>& data, std::vector<double> &a,std::vector<double> &b,std::vector<double> &c,std::vector<double> &d);
	public:
		//Constructor from data or a data file
			Interpolation();
			Interpolation(std::vector<std::vector<double>>& data);
			Interpolation(std::string filename,double dim1=1.0,double dim2=1.0);
		//Return values
			std::vector<std::vector<double>> Return_Data();
			std::vector<double> Return_Domain();
			double Return_Prefactor();
			std::vector<std::vector<double>> Return_Coefficients();
		//Set values
			void Set_Prefactor(double factor);
		//Interpolation
			double Interpolate(double x);
			double operator ()(double x)
		    {
		        return Interpolate(x);
		    };
		//Multiply by a constant
		    void Multiply(double factor);
		//Save function in a file
		    void Save_Function(std::string filename,unsigned int points);
	};
	
//Struct for weighted data points
	struct DataPoint
	{
		double value;
		double weight;
		//Constructor
			DataPoint(double v,double w);
			DataPoint(double v);
			DataPoint();

	};
	//Operators
		bool operator <(const DataPoint& lhs,const DataPoint& rhs);
		bool operator >(const DataPoint& lhs,const DataPoint& rhs);
		bool operator ==(const DataPoint& lhs,const DataPoint& rhs);
		std::ostream& operator <<(std::ostream &output,const DataPoint& dp);

#endif