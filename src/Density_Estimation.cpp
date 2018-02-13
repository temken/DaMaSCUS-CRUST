#include "Density_Estimation.hpp"

#include <iostream>
#include <fstream>
#include <cmath>

//Kernels
	//Gauss
		double Gaussian_Kernel(double x)
		{
			return 1.0/sqrt(2.0*M_PI)*exp(-x*x/2.0);
		}
//KDE
	Interpolation Perform_KDE(std::vector<DataPoint> data,double xMin,double xMax,double bw)
	{
		//Data points
			unsigned int N_Data = data.size();
			double Weight_Sum=0.0;
		//Count the weights
			for(unsigned int i = 0 ; i<N_Data ; i++) Weight_Sum += data[i].weight;

		//1. Bandwidth selection:
		//If the bandwidth is not set manually, we find it here.
			if(bw==0)
			{
				//1.1 Average.
					double AverageSum=0.0;
					for(unsigned int i = 0 ; i<N_Data ; i++)
					{
						AverageSum+=data[i].weight*data[i].value;
					}
					double Average = AverageSum/Weight_Sum;
				//1.2 Standard deviation
					double Variance=0.0;
					for(unsigned int i = 0 ; i<N_Data ; i++)
					{
						Variance += data[i].weight * pow(data[i].value - Average,2.0) / Weight_Sum;
					}
				//1.3 Bandwidth with rule-of-thumb estimator
					bw = sqrt(Variance)*pow(4.0/3.0/N_Data,0.2);
			}
		//2. Perform and tabulate the KDE
			int points = 150;
			double dx = (xMax-xMin)/(points-1);
			std::vector<std::vector<double>> Interpol_List;
			for(int j =0;j<points;j++)
			{
				double x = xMin + j*dx;
				double kde=0.0;
				for(unsigned int i = 0;i<N_Data;i++)
				{
					//Vanilla Gauss KDE
						kde+= data[i].weight*Gaussian_Kernel((x-data[i].value)/bw);
					//Reflection method
						// double xRefl = 2.0*xMin-data[i];
						// kde+= Gaussian_Kernel((x-data[i])/bw)+Gaussian_Kernel((x-xRefl)/bw);
						// kde+=weights[i]*(Gaussian_Kernel((x-data[i])/bw)+Gaussian_Kernel((x-xRefl)/bw));
				}
				//Pseudo data. Cowling and Hall's three-point rules
					//Sort data:
						std::sort(data.begin(), data.end());
					// 	Generate pseudo data
						int N_PseudoData = N_Data/3.0;
						for(int j = 0; j < N_PseudoData ; j++) 
						{
							double xPseudo = 4.0*xMin-6.0*data[j].value+4.0*data[2*j].value-data[3*j].value;
							// double wPseudo = (6.0*data[j].weight+4.0*data[2*j].weight+data[3*j].weight)/9.0;
							double wPseudo = (data[j].weight+data[2*j].weight+data[3*j].weight)/3.0;
							kde+= wPseudo*Gaussian_Kernel((x-xPseudo)/bw);
						}
				kde/= bw*Weight_Sum;
				Interpol_List.push_back( std::vector<double> {x,kde});
			}
		//3. Interpolate the list.
			Interpolation result(Interpol_List);
		//4. Check normalization/ re-normalize.
			double norm = AdaptiveSimpsons(result,xMin,xMax,0.000001);
			result.Multiply(1.0/norm);
		//5. Output
			return result;
	}