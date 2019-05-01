//Contains auxiliary functions and classes for the MC simulations.
#ifndef __Simulation_Essentials_hpp_
#define __Simulation_Essentials_hpp_

#include <fstream>
#include <vector>
#include <random>
#include <functional>

#include "DM_Particle.hpp"

//1. Basic vector class
	class Vector3D
	{
		private:
			double comp[3];
		public:
			// Constructors:
			Vector3D();
			Vector3D(double x, double y, double z);
			Vector3D(const Vector3D& v);
			
			// Operator overloadings:
			Vector3D operator+(Vector3D v);
	 		Vector3D operator-(Vector3D v);
			Vector3D operator/(double s);
	  		Vector3D operator=(Vector3D v);
	  		Vector3D& operator+= (const Vector3D& v);
	  		Vector3D& operator-=(const Vector3D& v);
	  		double& operator[](int i)       { return comp[i]; }
    		const double& operator[](int i) const { return comp[i]; }
	  		
	  		// Functions:
	  		double norm();
	  		void normalize();
	  		Vector3D normalized();
	  		double dot(const Vector3D& v);
	  		Vector3D cross(const Vector3D& v);
	};
	Vector3D operator*(const Vector3D &v,double s);
	Vector3D operator*(double s,const Vector3D &v);
	std::ostream& operator <<(std::ostream &output,const Vector3D& v);
	bool operator ==(const Vector3D& v1,const Vector3D& v2);

	//Coordinate Systems
	Vector3D SphericalCoordinates(double r,double theta,double phi);

//2. Event Class
	class Event
	{
		double time;
		Vector3D position;
		Vector3D velocity;
		public:
			//Constructors
			Event();
			Event(double t, Vector3D& x, Vector3D& v);
			Event(double t, double x, double y, double z, double vx, double vy, double vz);
			
			void SetTime(double t);
			void SetPosition(Vector3D& newposition);
			void SetPosition(double x,double y,double z);
			void SetVelocity(Vector3D& newvelocity);
			void SetVelocity(double vx,double vy,double vz);
			
			void IncreaseTime(double dt);
			void IncreasePosition(Vector3D dx);
			
			//Return Values
			double Time() const;
			Vector3D Position() const;
			Vector3D Velocity() const;
			Event kmsec();
			
			//Vector Norms
			double Speed();
			
			//Overloading the output operator <<
			friend std::ostream& operator<<(std::ostream &output,const Event &event);
	};

//3. Shielding layers
	class Layer
	{
		private:
			std::vector<std::vector<double>>composition;

			//DM mean free path
			bool tabulated;
			double mfp;
			std::vector<double> prob;
			//Number of values in table
			int N=10000;
			double vMin;
			std::vector<double> mfp_array;
			std::vector<std::vector<double>> prob_array;

		public:
			int index;
			std::string name;
			double density;
			double thickness;
			double depth;
			
			//Constructor
			Layer(std::string& ID,unsigned int ind, double rho, double thick,double deep,std::vector<std::vector<double>>& comp);
			
			//Mean free path			
			void Compute_MFP(const DM_Particle &DM,double vmin=0.0);
			double MFP(double vDM);
			
			//Stopping power
			double Stopping_Power(const DM_Particle& DM, double vDM) const;
			std::vector<double> Sample_Target(double vDM,std::mt19937& PRNG);
			void Print_Summary(); 
	};
	extern void Compute_All_MFP(const DM_Particle& DM, double vmin);
	extern int Current_Layer(Vector3D& x);
	extern int Current_Layer(Event& event);
	extern Event Leave_Layer(Event& x);
	extern void Compare_StoppingPower(const DM_Particle& DM,double vDM,const std::vector<Layer>& layers,int rank);
	//Atmosphere
	extern std::vector<Layer> Atmospheric_Layers(int layers,double altitude);

//4. Random number generation and sampling
	extern double ProbabilitySample(std::mt19937& PRNG);
	extern double ThetaSample(std::mt19937& PRNG);
	extern double PhiSample(std::mt19937& PRNG);
	extern double Rejection_Sampling(const std::function<double(double)>& PDF,double xMin,double xMax,double yMax,std::mt19937& PRNG);

//5. Functions for splitting/Russian roulette
	extern bool Survive_Russian_Roulette(double pSurv,std::mt19937& PRNG);
	extern int Importance_Domains(const DM_Particle& DM,double vDM,double vCutoff,double kappa,std::vector<Layer>& layers);
	extern void Compute_Importance_Boundaries(const DM_Particle& DM,double vDM,int N_layers,double splits,std::vector<Layer>& Layers);
	extern double Importance(double z);

#endif