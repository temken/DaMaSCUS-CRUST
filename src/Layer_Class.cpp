#include "Layer_Class.hpp"

#include <iostream>
#include <iomanip>

#include "General_Utilities.hpp"

//Layer Class
	//Constructors:
		Layer::Layer(std::string& ID,unsigned int ind, double rho, double thick,double deep,std::vector<std::vector<double>>& comp)
		{
			//ID:
				name=ID;
				index=ind;
			//Dimension and location:
				depth=deep;
				thickness=thick;
			//Content:
				density=rho;
				composition=comp;
		}
	//Set Functions
		void Layer::Set_Name(std::string newname)
		{
			name=newname;
		}
		void Layer::Set_Density(double newdensity)
		{
			density=newdensity;
		}
		void Layer::Set_Thickness(double newthickness)
		{
			thickness=newthickness;
		}
		void Layer::Set_Composition(std::vector<std::vector<double>>& newcomposition)
		{
			composition=newcomposition;
		}
	//Get Functions
		std::string Layer::Get_Name()
		{
			return name;
		}
		double Layer::Get_Density()
		{
			return density;
		}
		double Layer::Get_Thickness()
		{
			return thickness;
		}
		double Layer::Get_Depth()
		{
			return depth;
		}
		std::vector<std::vector<double>> Layer::Get_Composition()
		{
			return composition;
		}
		double Layer::Get_MFP()
		{
			return mfp;
		}
		std::vector<std::vector<double>> Layer::Get_Probability()
		{
			return probability;
		}
	//Further Functions

		//Compute MFP and scattering probabilities
			void Layer::Compute_MFP(double mDM, double sigma,double speed,bool ldm)
			{
				//Delete previous data
					probability.clear();
					mfp=0.0;
				for(unsigned int i = 0; i < composition.size();i++)
				{
					double f = composition[i][0];
					double Z = composition[i][1];
					double A = composition[i][2];
					double lambdai = density*f/NucleusMass(A)*Total_CrossSection(mDM,sigma,Z,A,speed,ldm);
					mfp+=lambdai;
					std::vector<double> prob {lambdai,Z,A};
					probability.push_back(prob);
				}
				//Normalize probabilities
					for(unsigned int i = 0; i < probability.size();i++)
					{
						probability[i][0]=probability[i][0]/mfp;
					}
				mfp=1.0/mfp;
			}

		void Layer::Print_Summary()
		{
			cout 	<<"\tLayer "<<index<<":\t\t"		<<name 										<<endl
					<<"\tDepth [m]:\t\t" <<InUnits(depth,meter)<<endl
					<<"\tThickness [m]:\t\t"	<<InUnits(thickness,meter)			<<endl
					<<"\tDensity [g cm^-3]:\t"		<<InUnits(density,gram/cm/cm/cm)	<<endl
					<<"\tComposition:";
			if(probability.size()==0) cout<<"\t\tf[\%]\tZ\t A"	<<endl;
			else cout<<"\t\tf[\%]\tZ\tA\tp[\%]" <<endl;
			for(unsigned int i = 0;i<composition.size();i++)
			{
				cout <<"\t\t\t\t"<<100*composition[i][0]<<"\t"<<composition[i][1] <<"\t"<<composition[i][2];
				if(probability.size()==0) cout <<endl;
				else cout <<"\t" <<100*probability[i][0]<<endl;
			}
			if(probability.size()!=0) cout <<"\tMFP[m]:\t\t\t" <<InUnits(mfp,meter)<<endl;
			cout <<endl;
		}

//Compute mfp for all existing layers:
	void Compute_All_MFP(double mDM,double sigma, bool ldm,double speed)
	{
		for(unsigned i = 0; i<Layers.size();i++) Layers[i].Compute_MFP(mDM,sigma,speed,ldm);
	}

//Return the current layer
	int Current_Layer(Eigen::Vector3d& x)
	{
		int output=-1;
		double z = x[2];
		if(z>0) output = 0;
		else if (z<=(-1.0)*Detector_Depth) output=Detector_Index;
		else
		{
			double depth=0.0;
			for(unsigned int i=0; i<Layers.size();i++)
			{
				depth-=Layers[i].Get_Thickness();
				if(z>depth)
				{
					output=i+1;
					break;
				}
			}
		}
		return output;
	}
	int Current_Layer(Event& event)
	{
		Eigen::Vector3d x = event.Position();
		int output=-1;
		double z = x[2];
		if(z>0) output = 0;
		else if (z<=(-1.0)*Detector_Depth) output=Detector_Index;
		else
		{
			double depth=0.0;
			for(unsigned int i=0; i<Layers.size();i++)
			{
				depth-=Layers[i].Get_Thickness();
				//Particle sits on the upper boundary flying upwards
				if(z==(-1.0)*Layers[i].Get_Depth() && event.Velocity()[2]>0)
				{
					output=i;
					break;
				}
				else if(z>depth)
				{
					output=i+1;
					break;
				}
			}
		}
		return output;
	}
//Event when particle leave the current layer
	Event Leave_Layer(Event& x)
	{
		Event output=x;
		//Current layer
			Eigen::Vector3d posi = x.Position();
			int layer = Current_Layer(x);
		//Find t_Exit
			double z =posi[2];
			double vz =x.Velocity()[2];
			double zExit=(-1.0)*Layers[layer-1].Get_Depth();
			if(vz<0) zExit -= Layers[layer-1].Get_Thickness();
			double tExit = (zExit-z)/vz;
		//Set the new time and position
			double tNew =x.Time()+tExit;
			Eigen::Vector3d xNew = x.Position()+tExit*x.Velocity();
			//push it a tiny bit into the new layer for numerical reasons.
				xNew+=(1.0e-7*meter)*x.Velocity().normalized();
			output.SetTime(tNew);
			output.SetPosition(xNew);

		return output;
	}

