//Disclaimer:
//Some of the function implemenations were made with the help of the 
//"Numerical Recipes 3rd Edition: The Art of Scientific Computing"
//by William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery

#include "Numerics_Functions.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>

using namespace std;

//1. Special functions
	//1.1 Simple functions
		int Sign(double arg)
		{
			if(arg>0.0) 		return 1;
			else if(arg==0.0)	return 0;
			else 				return -1;
		}
		double Sign(double x, double y)
		{
			if(Sign(x) == Sign(y)) return x;
			else return -1.0*x;
		}
		double StepFunction(double x)
		{
			if(x>=0) 	return 1.0;
			else		return 0.0;
		}
		//This function should only be used to 
		double Round(double N,unsigned int digits)
		{
			if(N==0) return 0;
			if(digits>5)
			{
				std::cerr <<"Error in Round(): Significant digits > 5."<<std::endl;
				std::exit(EXIT_FAILURE);
			}
			//Make the argument a positive number.
				double sign = Sign(N);
				N*=sign;
			//Cut of the decimal power
				double DecimalPower = floor( log10(N) );
				double prefactor = N*pow(10,-DecimalPower);
			//Cut off all unneeded decimal points, which might cause rounding error
				prefactor = (int) floor(0.5+1.0* pow(10,digits)*prefactor);
				prefactor/=10;
			//Round the prefactor with k significant digits
				prefactor = floor(0.5 +1.0*prefactor);
				prefactor/= pow(10,digits-1);
			//Final output
				return sign*prefactor*pow(10,DecimalPower);
		}

	//1.2 Gamma Functions from "Numerical Recipes"
		//Factorial
			std::vector<double> FactorialList= {1.0};
			double Factorial(unsigned int n)
			{
				if (n>170)
				{

					std::cerr <<"Error in Factorial: Overflow for " <<n <<"!."<<std::endl;
					std::exit(EXIT_FAILURE);
				}
				else if(n<FactorialList.size()) return FactorialList[n];
				else
				{
					while(FactorialList.size()<=n)
					{
						FactorialList.push_back( FactorialList.back()*FactorialList.size());
					}
					return FactorialList.back();
				}
			}
		//Binomial coefficients
			double BinomialCoefficient(int n,int k)
			{
				if(k<0 ||n<0 )
				{
					std::cerr <<"Warning in BinomialCoefficient(): negative arguments. Return 0." <<std::endl;
					std::exit(EXIT_FAILURE);
				}
				else if(n<k)	return 0;
				else if (n>170)
				{
					return floor(0.5+exp(GammaLn(n+1.0)-GammaLn(k+1.0)-GammaLn(n-k+1.0)));
				}
				else	return floor(0.5+Factorial(n)/Factorial(k)/Factorial(n-k));
			}
		//Logarithmic gamma function
			double GammaLn(double x)
			{
				double cof[14] = {57.1562356658629235,-59.5979603554754912, 14.1360979747417471,-0.491913816097620199,.339946499848118887e-4, .465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3, -.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3, .844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};
				if(x<=0)
				{
					std::cerr<<"Error in GammaLn(x): x<=0."<<std::endl;
					std::exit(EXIT_FAILURE);
				}
				double sum = 0.999999999999997092;
				double y = x;
				double tmp = x+671.0/128.0;
				tmp = (x+0.5)*log(tmp)-tmp;
				for(int j=0;j<14;j++)
				{
					sum+=cof[j]/++y;
				}
				return tmp+log(2.5066282746310005*sum/x);
			}
		//Q(x,a) via integration;
			double GammaQint(double x,double a)
			{
				//Compute P(x,a)
					double gammaP;
					double gln=GammaLn(a);
					//How far to integrate N sqrt(a) around the peak at a-1:
						double N = 10;
						double tPeak = a-1.0;
						double tMin = std::max(0.0, tPeak - N*sqrt(a));
						double tMax = tPeak + N*sqrt(a);
		 			//If x lies far away from the peak:
						if(x>tMax) gammaP = 1.0;
						else if (x<tMin ) gammaP = 0.0;
					//Numerical integration
						else
						{
							
							//integrand
								std::function<double(double)> integrand = [a,gln] (double t)
								{
									return exp(-gln-t+log(t)*(a-1.0));
								};
								if(x<tMin) tMin = 0.0;
							//Precision
								double eps = Find_Epsilon(integrand,tMin,x,1e-5);
							//Integrate
								gammaP = Integrate(integrand,tMin,x,eps);
						}
				
				return 1.0-gammaP;
			}
		//Series expansion of P(x,a)
			double GammaPser(double x,double a)
			{
				double eps = std::numeric_limits<double>::epsilon();
				double sum,del,ap;
				double gln=GammaLn(a);
				ap = a;
				del = sum = 1.0/a;
				while(fabs(del)>fabs(sum)*eps)
				{
					ap++;
					del*=x/ap;
					sum+=del;
				}
				return sum*exp(-x+a*log(x)-gln);
			}
		//Continued fraction representation of Q(x,a)
			double GammaQcf(double x,double a)
			{
				//Precision
					double eps = std::numeric_limits<double>::epsilon();
					double FPMIN = std::numeric_limits<double>::min()/eps;
				double del=0.0;
				double gln=GammaLn(a);
				double b = x+1.0-a;
				double c=1.0/FPMIN;
				double d = 1.0/b;
				double h=d;
				int i=1;
				while(fabs(del-1.0)>eps)
				{
					double an = -1.0*i*(i-a);
					b+=2.0;
					d=an*d+b;
					if(fabs(d)<FPMIN) d=FPMIN;
					c=b+an/c;
					if(fabs(c)<FPMIN) c =FPMIN;
					d=1.0/d;
					del=d*c;
					h*=del;
				}
				return exp(-x+a*log(x)-gln)*h;
			}
		//Final function using different methods for different parts of the domain
			double GammaQ(double x,double a)
			{
				double aMax=100.0;
				if(x<0.0||a<=0.0)
				{
					std::cerr <<"Error in GammaQ("<<x<<","<<a<<"): Invalid arguments."<<std::endl;
					std::exit(EXIT_FAILURE);
				}
				else if(x==0)return 1.0;
				else if (a>aMax) return GammaQint(x,a);
				else if (x<a+1.0) return 1.0-GammaPser(x,a);
				else return GammaQcf(x,a);
			}
			double GammaP(double x,double a)
			{
				return 1.0-GammaQ(x,a);
			}
			
		//Inverse incomplete gamma function. 
			//Solves P(x,a)=p for x.
				double Inv_GammaP(double p,double a)
				{
					
					//Check the arguments
						if(a<=0.0)
						{
							std::cerr <<"Error in Inv_GammaP(): a must be positive."<<std::endl;
							std::exit(EXIT_FAILURE);
						}
						if(p>=1.0) return std::max(100.0,a+100.*sqrt(a));
						if(p<=0.0) return 0.0;
					//Parameter
						double x;
						double gln=GammaLn(a);
						double a1=a-1.0;
						double lna1 = log(a1);
						double afac = exp(a1*(lna1-1.0)-gln);
					//Initial guess 1
						if(a>1.0)
						{
							double pp = (p<0.5)? p : 1.0-p;
							double t = sqrt(-2.0*log(pp));
							x = (2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t;
							if(p<0.5) x=-x;
							x = std::max(1.0e-3,a*pow(1.0-1.0/(9.*a)-x/(3.*sqrt(a)),3.0));

						}
					//Initial guess 2
						else
						{
							double t = 1.0 - a*(0.253+a*0.12);
					        if (p < t) x = pow(p/t,1./a);
					        else x = 1.-log(1.-(p-t)/(1.-t));
						}
					//Halley's method
						double EPS = 1.0e-8;
						for(int i=0;i<12;i++)
						{
							if(x<=0.0) return 0.0;
							double error = GammaP(x,a) - p;
							double t;
							if(a>1.0) t = afac*exp(-(x-a1)+a1*(log(x)-lna1));
							else t = exp(-x+a1*log(x)-gln);
							double u = error/t;
							x-= (t = u/(1.-0.5*std::min(1.,u*((a-1.)/x - 1))));
							if (x <= 0.) x = 0.5*(x + t);
							if (fabs(t) < EPS*x ) break;
						}
					return x;
				}
			//Solves Q(x,a)=p for x.
				double Inv_GammaQ(double q,double a)
				{
					return Inv_GammaP(1.0-q,a);
				}



	//1.3 Statistics
		double PDF_Gauss(double x,double mu,double sigma)
		{
			return 1.0/sqrt(2.0*M_PI)/sigma*exp(-pow((x-mu)/sigma,2.0)/2.0);
		}
		double PDF_Binomial(int mu,double p,int trials)
		{
			return BinomialCoefficient(mu,trials)*pow(p,trials)*pow(1.0-p,(mu-trials));
		}
		double PDF_Poisson(double expected_events, unsigned int events)
		{
			if(expected_events==0&&events==0) return 1.0;
			else if(expected_events==0&&events>0) return 0.0;
			else
			{
				double sum=events*log(expected_events) -expected_events;
				for(unsigned int i=2;i<=events;i++)sum-= log(i);
				return exp(sum);
			}
		}
		double CDF_Poisson(double expectation_value,unsigned int observed_events)
		{
			double gq = GammaQ(expectation_value,observed_events+1);
			if (gq>=0) return gq;
			else return 0.0;
		}
		double Inv_CDF_Poisson(unsigned int observed_events, double cdf)
		{
			if(observed_events==0) return (-1.0)*log(cdf);
			else return Inv_GammaQ(cdf,observed_events+1); 	
		}
		//Maximum Gap a la Yellin
			double MaxGap_C0(double x,double mu)
			{
				if(x==mu) return 1.0-exp(-mu);
				else
				{
					int m = mu/x;
					double sum=0.0;
					for(int k=0;k<=m;k++) 
					{
						double term = pow(k*x-mu,k)/Factorial(k)*exp(-k*x)*(1.0+k/(mu-k*x));
						sum+= term ;
						if (fabs(term)<1e-20) break;
					}
					return sum;
				}
			}
			double MaxGap_C0_error(double x,double mu,double mu_error)
			{
				if(x==mu) return mu_error*exp(-mu);
				else
				{
					int m = mu/x;
					double sum=0.0;
					for(int k=0;k<=m;k++) 
					{
						double term = k*pow(k*x-mu,(-2.0+k))/Factorial(k)*exp(-k*x)*(1.0+k*(x-1.0)-mu);
						sum+= term ;
						if (fabs(term)<1e-20) break;
					}
					return fabs(sum*mu_error);
				}
			}

//2. Numerical Integration via Adaptive Simpson Method
	//Function to return a reasonable precision.
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
	//Recursive functions for one-dimensional integration
		double Integrate(std::function<double(double)> func, double a,double b, double epsilon,int maxRecursionDepth)
		{
			int sign = +1;
			if(a==b) return 0.0;
			else if(a>b)
			{
				double aux = a;
				a = b;
				b=aux;
				sign = -1;
			}
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
					std::cout <<"Warning in Integrate(): Numerical integration on the interval ("<<a<<","<<b<<") did not converge to the desired precision." <<std::endl;
					std::cout <<"\tDesired precision: " <<fabs(epsilon) <<" Result: " <<result<<std::endl;
				}
			if(std::isnan(result)) std::cout <<"Warning in Integrate(): Result is nan."<<std::endl;
			else if(std::isinf(result)) std::cout <<"Warning in Integrate(): Result is inf."<<std::endl;
			return sign*result;
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

//3. Finding roots
	//Root finding with Ridder's method
	double Find_Root(std::function<double(double)>& func,double xLeft, double xRight,double xAccuracy)
	{
		const int Max_Iterations = 50;
		//1. Check if xLeft<xRight, otherwise swap.
			if(xLeft>xRight)
			{
				double temp = xLeft;
				xLeft = xRight;
				xRight = temp;
			}
		//2. Compute functions at boundary
			double fLeft = func(xLeft);
			double fRight = func(xRight);
		
		//3. Check if xLeft and xRight bracket a root or already yield a root. Also check for NaN's.
			if(std::isnan(fLeft)||std::isnan(fRight))
			{
				std::cerr <<"Error in Find_Root(): Function returns nan at the brackets."<<std::endl;
				std::exit(EXIT_FAILURE);
			}
			else if(fLeft*fRight>=0.0)
			{
				if(fLeft==0) return xLeft;
				else if(fRight==0) return xRight;
				else
				{
					std::cerr<<"Error in Find_Root(): f(xLeft)*f(xRight) = ("<<func(xLeft)<<")*("<<func(xRight) <<")>0.0"<<std::endl;
					std::exit(EXIT_FAILURE);
				}
			}
		//4. Ridder's method
			else
			{
				double x1 = xLeft;
				double x2 = xRight;
				double f1 = fLeft;
				double f2 = fRight;
				double result = -9.9e99;
				for(int i=0;i<Max_Iterations;i++)
				{
					//Mid point
						double x3 = (x1+x2)/2.0;

						double f3 = func(x3);
					//New point
						double x4 = x3 +(x3-x1) * Sign(f1-f2)*f3/sqrt(f3*f3-f1*f2);
					//Check if we found the root
						if(fabs(x4-result)<xAccuracy) return x4;
					//Prepare next iteration
						result = x4;
						double f4 = func(x4);
						if(f4==0.0) return result;
						//a) x3 and x4 bracket the root
						if(Sign(f3,f4) != f3)
						{
							x1 = x3;
							f1 = f3;
							x2 = x4;
							f2 = f4;
						}
						//b) x1 and x4 bracket the root
						else if(Sign(f1,f4) != f1)
						{
							x2 = x4;
							f2 = f4;

						}
						//c) x2 and x4 bracket the root
						else if(Sign(f2,f4) != f2)
						{
							x1 = x4;
							f1 = f4;
						}
						else
						{
							std::cerr <<"Error in Find_Root(). Ridder's method does not reach the root."<<std::endl;
							std::exit(EXIT_FAILURE);
						}
				}
				std::cout <<"Warning in Find_Root(): Iterations exceed the maximum. Final value f("<<result<<")="<<func(result)<<std::endl;
				return result;
			}
	}


//4. One-dimensional interpolation of tabulated data using Steffen splines 
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
								dy[i] = (Sign(p[i])+Sign(s[i]))*std::min(1.0*abs(s[i]),0.5*abs(p[i]));
							}
						//Last point
							else if(i==N-1)
							{
								p[i] = s[i-1]*(1.0+h[i-1]/(h[i-1]+h[i-2]))-s[i-2]*h[i-1]/(h[i-1]+h[i-2]);
								dy[i] = (Sign(p[i])+Sign(s[i-1]))*std::min(1.0*abs(s[i-1]),0.5*abs(p[i]));
							}
						//Points in the middle
							else
							{
								p[i] = (s[i-1]*h[i]+s[i]*h[i-1])/(h[i-1]+h[i]);
								dy[i] = (Sign(s[i-1])+Sign(s[i]))*std::min(1.0*abs(p[i])/2.0 ,std::min(1.0*abs(s[i]) ,1.0*abs(s[i-1]) ) );
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
							std::cout <<"Warning: Steffen coefficients in interpolation are NAN."<<std::endl;
						}
					}
		}
		unsigned int Interpolation::Bisection(double x,int jLeft,int jRight)
		{
			while((jRight-jLeft)>1)
			{
				int jm = (jRight+jLeft) >>1 ;
				if (x >= TabulatedData[jm][0])	jLeft=jm;
				else jRight=jm;
			}
			return jLeft;
		}
		unsigned int Interpolation::Hunt(double x)
		{

			// 1. Hunting phase returns jd,ju which bracket j
				int dj = 1;
				int jd;
				unsigned int ju;
				//Hunt up
				if(x>TabulatedData[jLast][0])
				{
					jd = jLast;
					ju = jd+dj;
					while(x>TabulatedData[ju][0])
					{
						jd = ju;
						ju += dj;
						//Check if we ran off the range:
						if(ju > N_Data-1)
						{
							ju = N_Data-1;
							break;
						}
						else
						{
							dj += dj;
						}
					}
				}
				//Hunt down
				else if (x<TabulatedData[jLast][0])
				{
					ju = jLast;
					jd = ju-dj;
					while(x<TabulatedData[jd][0])
					{
						ju = jd;
						jd -= dj;
						//Check if we ran off the range:
						if(jd < 0)
						{
							jd = 0;
							break;
						}
						else
						{
							dj += dj;
						}
					}
				}
				else
				{
					return jLast;
				}
			// 2. Bisection phase
				if((ju-jd)>1) jd = Bisection(x,jd,ju);
			return jd;
		}
		// Find j such that list[j]<x<list[j+1]
		unsigned int Interpolation::Locate(double x)
		{
			if( ((xDomain[0]-x)>0.0) || ((x-xDomain[1])>0.0) )
			{
				printf("\nError in Interpolation::Locate(): x = %e lies outside the domain [%e,%e].\n\n",x,xDomain[0],xDomain[1]);
				std::exit(EXIT_FAILURE);
			}
			else
			{
				//Use Bisection() or the Hunt method, depending of the last calls were correlated.
				unsigned int j = corr ? Hunt(x): Bisection(x,0,N_Data-1);
				//Check if the points are still correlated.
				corr=(fabs((j-jLast))<10);

				jLast = j;
				return j;
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
				jLast = 0;
				corr=false;
			//Compute coefficients
				Compute_Steffen_Coefficients(data,a,b,c,d);

		}
		Interpolation::Interpolation(const std::string& filename,double dim1,double dim2)
		{
			//1. Import data.
				std::ifstream inputfile;
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
		        	std::cerr << "Error in Interpolation(" <<filename<<"): File does not exist."<<std::endl;
		    		std::exit(EXIT_FAILURE);
		    	}
			//2. Sort data.
				std::sort(TabulatedData.begin(), TabulatedData.end());
			//3. Interpolate data.
				preFactor = 1.0;
				N_Data = TabulatedData.size();
				jLast = 0;
				corr=false;
				//3.1 Find the function's domain:
					std::vector <double> x;
					for(unsigned i = 0; i < TabulatedData.size(); i++) x.push_back(TabulatedData[i][0]);
					xDomain.push_back( *min_element(x.begin(),x.end()));
					xDomain.push_back( *max_element(x.begin(),x.end()));
				//3.2 Compute the Steffen coefficients for the interpolation
					Compute_Steffen_Coefficients(TabulatedData,a,b,c,d);

		}
		Interpolation::Interpolation(const std::vector<std::vector<double>>& data,double dim1,double dim2)
		{
			preFactor = 1.0;
			for(unsigned int i =0;i<data.size();i++)
			{
				std::vector<double> aux ={data[i][0]*dim1,data[i][1]*dim2};
				TabulatedData.push_back(aux);
			}
			// TabulatedData = data;
			N_Data = TabulatedData.size();
			jLast = 0;
			corr=false;
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
			int j = Locate(x);
			double x_j = TabulatedData[j][0];
			return preFactor*(a[j]*pow(x-x_j,3.0)+b[j]*pow(x-x_j,2.0)+c[j]*(x-x_j)+d[j]);
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

//5. Struct for weighted data points
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
	//Weighted mean
		std::vector<double> Weighted_Average(std::vector<DataPoint>& data)
		{
			double sum =0.0;
			double wsum=0.0;
			long unsigned int N=data.size();
			//1. Average
				for(unsigned int i=0;i<N;i++)
				{
					sum+=data[i].weight*data[i].value;
					wsum+=data[i].weight;
				}
				double Average = sum/wsum;
				double wAverage= wsum/N;
			//2. Standard error (Cochran)
				double sum1=0.0,sum2=0.0,sum3=0.0;
				for(unsigned int i=0;i<data.size();i++)
				{
					sum1+=pow((data[i].weight*data[i].value-Average*wAverage),2.0);
					sum2+=(data[i].weight-wAverage)*(data[i].weight*data[i].value-Average*wAverage);
					sum3+=pow(data[i].weight-wAverage,2.0);
				}
				double SE=N/(N-1.0)/wsum/wsum*(sum1-2.0*Average*sum2+pow(Average,2.0)*sum3);
			//3. Return result
				return std::vector<double> {Average,sqrt(SE)};
		}

//6. Kernel density estimation
	double Gaussian_Kernel(double x)
	{
		return 1.0/sqrt(2.0*M_PI)*exp(-x*x/2.0);
	}
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
			//Sort data:
				std::sort(data.begin(), data.end());
			// Pseudo data
				unsigned int N_PseudoData = N_Data/3.0;
		//2. Perform and tabulate the KDE
			int points = 100;
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
					//Cowling and Hall pseudo data method
						if(i<N_PseudoData)
						{
							double xPseudo = 4.0*xMin-6.0*data[i].value+4.0*data[2*i].value-data[3*i].value;
							// double wPseudo = (6.0*data[i].weight+4.0*data[2*i].weight+data[3*i].weight)/9.0;
							double wPseudo = (data[i].weight+data[2*i].weight+data[3*i].weight)/3.0;
							kde+= wPseudo*Gaussian_Kernel((x-xPseudo)/bw);
						}
					//Reflection method
						// double xRefl = 2.0*xMin-data[i];
						// kde+= Gaussian_Kernel((x-data[i])/bw)+Gaussian_Kernel((x-xRefl)/bw);
						// kde+=weights[i]*(Gaussian_Kernel((x-data[i])/bw)+Gaussian_Kernel((x-xRefl)/bw));
				}				
				kde/= bw*Weight_Sum;
				Interpol_List.push_back( std::vector<double> {x,kde});
			}
		//3. Interpolate the list.
			Interpolation result(Interpol_List);
		//4. Check normalization/ re-normalize.
			double norm = Integrate(result,xMin,xMax,0.0001);
			result.Multiply(1.0/norm);
		//5. Output
			return result;
	}

//7. Read in table into vector
	std::vector<double> Read_List(std::string filepath,double dimension)
	{
		std::vector<double> list;
		std::ifstream inputfile;
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
        	std::cerr << "Error in Read_List(" <<filepath<<"): File does not exist."<<std::endl;
        	std::exit(EXIT_FAILURE);
    	}

		return list;
	}