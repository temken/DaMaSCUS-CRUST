#include "General_Utilities.hpp"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <cstdlib>
#include <libconfig.h++>
#include <algorithm>    // std::min_element, std::max_element

#include <sys/types.h> // required for stat.h
#include <sys/stat.h>	//required to create a folder

#include "Physical_Parameters.hpp"

using namespace libconfig;
//Version
	std::string version = "1.0";
//Input variables
	//Simulation ID
		std::string ID="default";
	//Statistical parameter
		int SampleSize=-1;
		//Certainty Level
			double CL=-1.0;
		//Importance sampling:
			double IS_Angle=0.0;	//Scattering Angle
			double IS_MFP=0.0;	//MFP
		//Parameter Scan
			double mMin=0.0;
			double mMax=0.0;
			double dm=0.0;
			int Masses=0;
	//Layer structure
		std::vector<Layer> Layers;
	//Dark Matter
		bool LDM=false;
	//Detector
		std::string Detector="default";

//Derived variables
	//Speed Cutoff
		double vCutoff=-1.0;
	//Detector Location
		double Detector_Depth=-1.0;
		int Detector_Index=-1;



//Copy cfg file
	void Copy_Config_File(const char *inputfile)
	{
		std::ifstream    inFile;
		std::ofstream outFile;
		inFile.open(inputfile);
    	outFile.open("../results/"+ID+"/"+ID+".cfg");
		outFile << inFile.rdbuf();
		inFile.close();
		outFile.close();
	}

//Read in cfg file
	void Read_Config_File(const char *inputfile,int rank)
	{
		std::string line="----------";
		Config cfg;
		// Read the file. If there is an error, report it and exit.
			try
			{		
				cfg.readFile(inputfile);
				if(rank==0) cout <<setw(10)<<"Config file:"<<std::setw(25)<<inputfile<<endl<<line<<endl;
			}
			catch(const FileIOException &fioex)
			{
				std::cerr << "I/O error while reading configuration file." << std::endl;
				exit(EXIT_FAILURE);
			}
			catch(const ParseException &pex)
			{
				std::cerr << "Configurate file parse error at " << pex.getFile() << ":" << pex.getLine() << " - " << pex.getError() << std::endl;
				exit(EXIT_FAILURE);
			}
		//Simulation ID
			try
			{
				ID = cfg.lookup("simID").c_str();
				if(rank==0) cout <<left<<setw(20)<<"\tSimulation ID:" <<ID <<endl;
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "No 'simID' setting in configuration file." << endl;
				exit(EXIT_FAILURE);
			}
		//Create folder for results
			if(rank==0) cout <<endl<<"\tCreating folder /results/"+ID <<"." <<endl;
			std::string sPath = "../results/"+ID;
			mode_t nMode = 0733; // UNIX style permissions
			int nError = 0;
			#if defined(_WIN32)
			  nError = _mkdir(sPath.c_str()); // can be used on Windows
			#else 
			  nError = mkdir(sPath.c_str(),nMode); // can be used on non-Windows
			#endif
			if (nError != 0) 
			{
				if(rank==0) cout <<"\tThe folder already exists, data will be overwritten."<<endl<<endl;
			}
			else cout <<endl;
		//Copy cfg file into this folder
			Copy_Config_File(inputfile);

		//Initial SampleSize
			try
			{
				SampleSize = cfg.lookup("samplesize");
				if(rank==0) cout <<left<<setw(20)<<"\tSample size:" <<SampleSize <<endl;
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "No 'samplesize' setting in configuration file." << endl;
				exit(EXIT_FAILURE);
			}
		//Certainty Level
			try
			{
				CL = cfg.lookup("cl");
				if(rank==0) cout <<left<<setw(20)<<"\tCertainty level:" <<CL <<endl;
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "No 'cl' setting in configuration file." << endl;
				exit(EXIT_FAILURE);
			}
		//Importance sampling parameter for the scattering angle:
			try
			{
				IS_Angle = cfg.lookup("is_angle");
				if(rank==0) cout <<left<<setw(20)<<"\tIS (angle):" <<IS_Angle <<endl;
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "No 'is_angle' setting in configuration file." << endl;
				exit(EXIT_FAILURE);
			}
		//Importance sampling parameter for the mean free path:
			try
			{
				IS_MFP = cfg.lookup("is_mfp");
				if(rank==0) cout <<left<<setw(20)<<"\tIS (MFP):" <<IS_MFP <<endl;
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "No 'is_mfp' setting in configuration file." << endl;
				exit(EXIT_FAILURE);
			}
		//De-/Activate light DM option
			try
			{
				LDM = cfg.lookup("LDM");
				if(rank==0) cout <<left<<setw(20)<<boolalpha<<"\tLight DM:" <<LDM <<endl;
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "No 'LDM' setting in configuration file." << endl;
				exit(EXIT_FAILURE);
			}
		//Experiment
			try
			{
				Detector = cfg.lookup("experiment").c_str();
				if(rank==0) cout <<left<<setw(20)<<"\tDetector:" <<Detector <<endl<<endl;
			}
			catch(const SettingNotFoundException &nfex)
			{
				cerr << "No 'experiment' setting in configuration file." << endl;
				exit(EXIT_FAILURE);
			}
		//Parameter scan
			if(rank==0) cout <<left<<setw(20)<<"\tParameter scan:" <<endl;
			//mMin
				try
				{
					mMin = cfg.lookup("mMin");
					if(rank==0) cout <<left<<setw(20)<<"\tmMin [GeV]:" <<mMin<<endl;
				}
				catch(const SettingNotFoundException &nfex)
				{
					cerr << "No 'mMin' setting in configuration file." << endl;
					exit(EXIT_FAILURE);
				}
			//mMax
				try
				{
					mMax = cfg.lookup("mMax");
					if(rank==0) cout <<left<<setw(20)<<"\tmMax [GeV]:" <<mMax <<endl;
				}
				catch(const SettingNotFoundException &nfex)
				{
					cerr << "No 'mMax' setting in configuration file." << endl;
					exit(EXIT_FAILURE);
				}
			//masses
				try
				{
					Masses = cfg.lookup("masses");
					if(rank==0) cout <<left<<setw(20)<<"\tMass steps:" <<Masses <<endl;
				}
				catch(const SettingNotFoundException &nfex)
				{
					cerr << "No 'masses' setting in configuration file." << endl;
					exit(EXIT_FAILURE);
				}
			//Number of steps
				if(mMin==mMax) Masses=1;
				if(Masses==1) 	dm=0;
				else 			dm = (log10(mMax)-log10(mMin))/(Masses-1.0);
		//Layer Structure
			double layer_depth=0.0;
			const Setting& root = cfg.getRoot();
			try
			{
				const Setting &layers = root["layers"];
				int count=layers.getLength();
				if(rank==0) cout <<endl<<left<<setw(20)<<"\tLayers:" <<count<<endl;
			 	if(rank==0) cout<<"\t"<<line<<endl
			 				<<"\tLayer 0:\tSpace"<<endl;
				for(int i=0;i<count;i++)
				{
					const Setting &layer =layers[i];
					std::string layer_name=layer.lookup("name").c_str();
					double layer_density=layer.lookup("density");
					double layer_thickness=layer.lookup("thickness");
					//Units
						layer_density*=gram/cm/cm/cm;
						layer_thickness*=meter;
					//Composition
						std::vector<std::vector<double>> layer_composition;
						int element_count=layer.lookup("composition").getLength();
						for(int j=0;j<element_count;j++)
						{
							std::vector<double> element;
							for(int k=0;k<3;k++) element.push_back(layer.lookup("composition")[j][k]);
							layer_composition.push_back(element);	
						}
					//Construct Layer
						unsigned int layer_index=i+1;
						Layer next_layer(layer_name,layer_index,layer_density,layer_thickness,layer_depth,layer_composition);
						if(rank==0) cout<<"\t"<<line<<endl;
						if(rank==0) next_layer.Print_Summary();
						Layers.push_back(next_layer);	
					//Increase depth by the layer's thickness
						layer_depth+=layer_thickness;
				}
				if(rank==0) cout<<"\t"<<line<<endl;
			}
			catch(const SettingNotFoundException &nfex)
			{
			    // Ignore.
			}
		//Derived Quantities
			Detector_Depth = layer_depth;
			Detector_Index = Layers.size()+1;
			if(rank==0) cout 	<<"\tLayer " <<Detector_Index<<":\t\tDetector" <<endl
								<<"\tType:\t\t\t" <<Detector <<endl
			 					<<"\tDepth:\t\t\t" <<InUnits(Detector_Depth,meter)<<"m" <<endl
			 					<<"\t"<<line<<endl
			 					<<line<<endl;
	}

//Read list from file
	std::vector<double> Read_List(std::string filepath,double dimension)
	{
		std::vector<double> list;
		ifstream inputfile;
		inputfile.open(filepath);
		if (inputfile.good())
		{

	        while (!inputfile.eof())
	        {
	        	double x;
	            inputfile >>x;
	            list.push_back(x*dimension);
	        }
	        // Close the file.
	        inputfile.close();
   		}
   		else
   		{
        	cout << "Error in Interpolation(" <<filepath<<"): File does not exist."<<endl;
    	}

		return list;
	}


//Work division
	std::vector<int> WorkDivision(int WorldSize,int MassSteps)
	{
		std::vector<int> output;
		
		int overlap = WorldSize*ceil(1.0*MassSteps/WorldSize)-MassSteps;
		for(int i =0 ; i<=WorldSize; i++)
		{
			if(i==0) 			output.push_back(0);
			else if(i<=overlap)	output.push_back(i*floor(1.0*MassSteps/WorldSize));
			else				output.push_back(output[i-1]+ceil(1.0*MassSteps/WorldSize));
		}
		return output;
	}	




//Numerical integration with adaptive Simpsons rule, adaptive error control based on absolute error
		double Find_Epsilon(std::function<double(double)> func, double a,double b,double precision)
		{
			double c = (a + b)/2;
			double h = b - a;                                                                  
  			double fa = func(a);
  			double fb = func(b);
  			double fc = func(c);                                                           
  			double S = (h/6)*(fa + 4*fc + fb);
  			double epsilon = precision*S;
  			return epsilon;
		}
		
		double AdaptiveSimpsons(std::function<double(double)> func, double a,double b, double epsilon,int maxRecursionDepth)
		{
			double c = (a + b)/2;
			double h = b - a;                                                                  
  			double fa = func(a);
  			double fb = func(b);
  			double fc = func(c);                                                           
  			double S = (h/6)*(fa + 4*fc + fb);
  			bool warning = false;
  			double result =   auxAdaptiveSimpson(func, a, b, fabs(epsilon), S, fa, fb, fc, maxRecursionDepth,warning); 
  			if(warning)   
			{
  				cout <<"Warning: Numerical integration on the interval ("<<a<<","<<b<<") did not converge to the desired precision." <<endl;
				cout <<"\tDesired precision: " <<fabs(epsilon) <<" Result: " <<result<<endl;;
			}
  			return result;
		}
		double auxAdaptiveSimpson(std::function<double(double)> func, double a,double b, double epsilon,double S,double fa,double fb,double fc,int bottom,bool &warning)
		{
			double c = (a+b)/2;
			double h = b-a;
			double d = (a+c)/2;
			double e = (b+c)/2;
			double fd=func(d);
			double fe=func(e);
			double Sleft = (h/12)*(fa + 4*fd + fc);                                                           
			double Sright = (h/12)*(fc + 4*fe + fb);                                                          
			double S2 = Sleft + Sright; 
			if (bottom <= 0 || abs(S2 - S) <= 15*epsilon)//15 due to error analysis 
			{
				if(bottom<=0&&abs(S2 - S) > 15*epsilon) warning=true;
				return S2 + (S2 - S)/15;  
			}                                         
			else
			{
				return auxAdaptiveSimpson(func, a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1,warning)+auxAdaptiveSimpson(func, c, b, epsilon/2, Sright, fc, fb, fe, bottom-1,warning); 
			}    
		}
//Root finding with regula falsi
	double RegulaFalsi(std::function<double(double)> func,double xLeft, double xRight,double epsilon)
	{
		double root;
		if(xLeft>xRight)
		{
			printf("Error in RegulaFalsi(): xLeft > xRight\n");
			root = 0.0;
		}
		else if(func(xLeft)*func(xRight)>=0.0)
		{
			printf("Error in RegulaFalsi(): f(xLeft)*f(xRight)>0.0\n");
			root = 0.0;
		}
		else
		{
			double x3=xLeft;
			while(abs(func(x3))>epsilon)
			{
				double x1 = xLeft;
				double x2 = xRight;
				double y1 = func(x1);
				double y2 = func(x2);
				x3 = (x1*y2-x2*y1)/(y2-y1);
				double y3 = func(x3);
				if(y3*y1>0.0)	xLeft = x3;
				else		 	xRight = x3;
				// cout <<y3<<endl;
			}
			root = x3;
		}
		return root;
	}
//Sign function
		double Sign(double arg)
		{
			if(arg>0.0) 		return 1.0;
			else if(arg==0.0)	return 0.0;
			else 				return -1.0;
		}

//Interpolation of tabulated data
	//Computation of the coefficients
		void Interpolation::Compute_Steffen_Coefficients(std::vector<std::vector<double>>& data, std::vector<double> &a,std::vector<double> &b,std::vector<double> &c,std::vector<double> &d)
		{
			unsigned int N= data.size();
			//Compute the Steffen coefficients for the interpolation
				//1. h and s.
					double h[N-1],s[N-1];
					for(unsigned int i=0;i<N-1;i++)
					{
						double x_i = data[i][0];
						double x_ip1 = data[i+1][0];
						double y_i = data[i][1];
						double y_ip1 = data[i+1][1];
						h[i] = x_ip1-x_i;
						s[i]= (y_ip1-y_i)/h[i];
					}
				//2. p and dy
					double dy[N],p[N];
					for(unsigned int i=0;i<N;i++)
					{
						//First point
							if(i==0)
							{
								p[i] = s[i]*(1.0+h[i]/(h[i]+h[i+1]))-s[i+1]*h[i]/(h[i]+h[i+1]);
								dy[i] = (Sign(p[i])+Sign(s[i]))*std::min(abs(s[i]),0.5*abs(p[i]));
							}
						//Last point
							else if(i==N-1)
							{
								p[i] = s[i-1]*(1.0+h[i-1]/(h[i-1]+h[i-2]))-s[i-2]*h[i-1]/(h[i-1]+h[i-2]);
								dy[i] = (Sign(p[i])+Sign(s[i-1]))*std::min(abs(s[i-1]),0.5*abs(p[i]));
							}
						//Points in the middle
							else
							{
								p[i] = (s[i-1]*h[i]+s[i]*h[i-1])/(h[i-1]+h[i]);
								dy[i] = (Sign(s[i-1])+Sign(s[i]))*std::min(abs(p[i])/2.0 ,std::min(abs(s[i]) ,abs(s[i-1]) ) );
							}
					}
				//3. a,b,c, and d
					for(unsigned int i=0;i<N-1;i++)
					{
						a.push_back((dy[i]+dy[i+1]-2.0*s[i])/pow(h[i],2.0));
						b.push_back((3.0*s[i]-2.0*dy[i]-dy[i+1])/h[i]);
						c.push_back(dy[i]);
						d.push_back(data[i][1]);
						if(std::isnan(a.back())||std::isnan(b.back())||std::isnan(c.back())||std::isnan(d.back())) 
						{
							cout <<"Warning: Steffen coefficients in interpolation are NAN."<<endl;
						}
					}
		}
	//Constructors

		Interpolation::Interpolation()
		{
			
			//Generate some data to interpolate zero
				std::vector<std::vector<double>> data;
				data.push_back(std::vector<double> {-1.0,0.0});
				data.push_back(std::vector<double> {0.0,0.0});
				data.push_back(std::vector<double> {+1.0,0.0});
			//Define members
				preFactor = 1.0;
				TabulatedData = data;
				N_Data = data.size();
				xDomain = {-1.0,1.0};
			//Compute coefficients
				Compute_Steffen_Coefficients(data,a,b,c,d);

		}
		Interpolation::Interpolation(std::string filename,double dim1,double dim2)
		{
			//1. Import data.
				ifstream inputfile;
				inputfile.open(filename);
				if (inputfile.good())
				{

			        while (!inputfile.eof())
			        {
			        	double x,y;
			            inputfile >>x;
			            inputfile >>y;
			            TabulatedData.push_back(std::vector<double> {x*dim1,y*dim2});
			        }
			        // Close the file.
			        inputfile.close();
		   		}
		   		else
		   		{
		        	cout << "Error in Interpolation(" <<filename<<"): File does not exist."<<endl;
		    	}
			//2. Sort data.
				std::sort(TabulatedData.begin(), TabulatedData.end());
			//3. Interpolate data.
				preFactor = 1.0;
				N_Data = TabulatedData.size();
				//3.1 Find the function's domain:
					std::vector <double> x;
					for(unsigned i = 0; i < TabulatedData.size(); i++) x.push_back(TabulatedData[i][0]);
					xDomain.push_back( *min_element(x.begin(),x.end()));
					xDomain.push_back( *max_element(x.begin(),x.end()));
				//3.2 Compute the Steffen coefficients for the interpolation
					Compute_Steffen_Coefficients(TabulatedData,a,b,c,d);

		}
		
		Interpolation::Interpolation(std::vector<std::vector<double>>& data)
		{
			preFactor = 1.0;
			TabulatedData = data;
			N_Data = TabulatedData.size();
			//Find the function's domain:
				std::vector <double> x;
				for(unsigned i = 0; i < TabulatedData.size(); i++) x.push_back(TabulatedData[i][0]);
				xDomain.push_back( *min_element(x.begin(),x.end()));
				xDomain.push_back( *max_element(x.begin(),x.end()));
			//Compute the Steffen coefficients for the interpolation
				Compute_Steffen_Coefficients(TabulatedData,a,b,c,d);
		}
	//Return values
		std::vector<std::vector<double>> Interpolation::Return_Data()
		{
			return TabulatedData;
		}
		std::vector<double> Interpolation::Return_Domain()
		{
			return xDomain;
		}
		double Interpolation::Return_Prefactor()
		{
			return preFactor;
		}
		std::vector<std::vector<double>> Interpolation::Return_Coefficients()
		{
			std::vector<std::vector<double>> output;
			for (unsigned i = 0 ; i < a.size() ; i++)
			{
				std::vector<double> aux {a[i],b[i],c[i],d[i]};
				output.push_back(aux);
			}
			return output;
		}
	//Set values
		void Interpolation::Set_Prefactor(double factor)
		{
			preFactor=factor;
		}
	//Interpolation
		double Interpolation::Interpolate(double x)
		{
			// if((x<xDomain[0]||x>xDomain[1])&&fabs(xDomain[1]-x)>1e-10&&fabs(xDomain[0]-x)>1e-10)
			if( ((xDomain[0]-x)>1e-10) || ((x-xDomain[1])>1e-10) )
			{
				// cout <<(x-xDomain[1])<<"\t" <<(x==xDomain[1]) <<endl;
				printf("\nWarning: In Interpolate() %e lies outside the domain [%e,%e].\n\n",x,xDomain[0],xDomain[1]);
				return 0;
			}
			else
			{
				//Steffen spline interpolation
				for(unsigned int i=0;i<N_Data-1;i++)
				{
					double x_i = TabulatedData[i][0];
					double x_ip1 = TabulatedData[i+1][0];
					if(x<=x_ip1) 
					{
						return preFactor*(a[i]*pow(x-x_i,3.0)+b[i]*pow(x-x_i,2.0)+c[i]*(x-x_i)+d[i]);
					}
				}

			}
			return 0;
		}
	//Multiply by a factor
		void Interpolation::Multiply(double factor)
		{
			preFactor*=factor;
		}
	//Save function
	    void Interpolation::Save_Function(std::string filename,unsigned int points)
	    {
	    	std::ofstream f;
	    	f.open(filename);
	    	double dx = (xDomain[1]-xDomain[0])/(points-1);
	    	for(unsigned int i = 0;i < points; i++)
	 		{
	 			double x = xDomain[0] + i * dx;
	 			f <<x <<"\t" <<Interpolate(x)<<std::endl;
	 		}
	    	f.close();
	    }

//weighted data struct
	//Constructors
	    DataPoint::DataPoint(double v, double w)
	    {
	    	value=v;
	    	weight=w;
	    }
	    DataPoint::DataPoint(double v)
	    {
	    	value=v;
	    	weight=1.0;
	    }
	    DataPoint::DataPoint()
	    {
	    	value=0.0;
	    	weight=1.0;
	    }
	//Operators
	    bool operator <(const DataPoint& lhs,const DataPoint& rhs)
	    {
	    	return (lhs.value<rhs.value);
	    }
		bool operator >(const DataPoint& lhs,const DataPoint& rhs)
		{
			return (lhs.value>rhs.value);
		}
		bool operator ==(const DataPoint& lhs,const DataPoint& rhs)
		{
			return (lhs.value==rhs.value);
		}
		std::ostream& operator <<(std::ostream &output,const DataPoint& dp)
		{
			output <<dp.value<<"\t" <<dp.weight;
			return output;
		}

