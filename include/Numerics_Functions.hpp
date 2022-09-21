// Contains numerical methods and special functions.
#ifndef __Numerics_Functions_hpp_
#define __Numerics_Functions_hpp_

#include <functional>
#include <vector>

// 1. Special functions
// 1.1 Simple functions
extern int Sign(double arg);
extern double Sign(double x, double y);	  // Returns x with the sign of y.
extern double StepFunction(double x);
extern double Round(double N, unsigned int digits = 3);

// 1.2 Gamma Functions
extern double Factorial(unsigned int n);
extern double BinomialCoefficient(int n, int k);
extern double GammaLn(double x);
extern double GammaQ(double x, double a);		// Incomplete Gamma Q(x,a)
extern double GammaP(double x, double a);		// Incomplete Gamma P(x,a)
extern double Inv_GammaP(double p, double a);	// Solves P(x,a)=p for x.
extern double Inv_GammaQ(double q, double a);	// Solves Q(x,a)=q for x.

// 1.3 Statistics
extern double PDF_Gauss(double x, double mu, double sigma);
extern double PDF_Binomial(int mu, double p, int trials);
extern double PDF_Poisson(double expected_events, unsigned int events);
extern double CDF_Poisson(double expectation_value, unsigned int observed_events);
extern double Inv_CDF_Poisson(unsigned int observed_events, double cdf);   // Solves the CDF = cdf for mu
extern double MaxGap_C0(double x, double mu);
extern double MaxGap_C0_error(double x, double mu, double mu_error);

// 2. Numerical Integration via Adaptive Simpson Method
extern double Find_Epsilon(std::function<double(double)> func, double a, double b, double precision);
extern double Integrate(std::function<double(double)> func, double a, double b, double epsilon, int maxRecursionDepth = 20);
double auxAdaptiveSimpson(std::function<double(double)> func, double a, double b, double epsilon, double S, double fa, double fb, double fc, int bottom, bool& warning);

// 3. Root finding
extern double Find_Root(std::function<double(double)>& func, double xLeft, double xRight, double epsilon);

// 4. One-dimensional interpolation of tabulated data using Steffen splines
class Interpolation
{
  private:
	std::vector<std::vector<double>> TabulatedData;
	unsigned int N_Data;
	std::vector<double> xDomain;
	// Steffen coefficients
	std::vector<double> a, b, c, d;
	// Pre-factor
	double preFactor;
	// compute Steffen coefficients
	void Compute_Steffen_Coefficients(std::vector<std::vector<double>>& data, std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d);
	// Locate j in array such that x[j]<x<x[j+1].
	unsigned int jLast;
	bool corr;	 // if successive calls are correlated, then the hunt method can be faster.
	unsigned int Bisection(double x, int jLeft, int jRight);
	unsigned int Hunt(double x);
	unsigned int Locate(double x);

  public:
	// Constructor from data or a data file
	Interpolation();
	Interpolation(const std::vector<std::vector<double>>& data, double dim1 = 1.0, double dim2 = 1.0);
	Interpolation(const std::string& filename, double dim1 = 1.0, double dim2 = 1.0);

	// Return values
	std::vector<std::vector<double>> Return_Data();
	std::vector<double> Return_Domain();
	double Return_Prefactor();
	std::vector<std::vector<double>> Return_Coefficients();

	// Set values
	void Set_Prefactor(double factor);

	// Interpolation
	double Interpolate(double x);
	double operator()(double x) { return Interpolate(x); };

	// Multiply by a constant
	void Multiply(double factor);

	// Save function in a file
	void Save_Function(std::string filename, unsigned int points);
};

// 5. Struct for weighted data points
struct DataPoint
{
	double value;
	double weight;

	// Constructor
	DataPoint(double v, double w);
	DataPoint(double v);
	DataPoint();
};
// Operators
bool operator<(const DataPoint& lhs, const DataPoint& rhs);
bool operator>(const DataPoint& lhs, const DataPoint& rhs);
bool operator==(const DataPoint& lhs, const DataPoint& rhs);
std::ostream& operator<<(std::ostream& output, const DataPoint& dp);

extern std::vector<double> Weighted_Average(std::vector<DataPoint>& data);

// 6. Kernel density estimation
// Kernels
extern double Gaussian_Kernel(double x);
// KDE
extern Interpolation Perform_KDE(std::vector<DataPoint> data, double xMin, double xMax, double bw = 0);

// 7. Read in table into vector
extern std::vector<double> Read_List(std::string filepath, double dimension = 1.0);

#endif