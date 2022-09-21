#include "Simulation_Essentials.hpp"

#include <algorithm>   // std::min_element, std::max_element
#include <cmath>
#include <iostream>

#include "Input_Parameters.hpp"
#include "Numerics_Functions.hpp"
#include "Physics_Functions.hpp"

// 1. Basic vector class
// Constructors
Vector3D::Vector3D()
{
	comp[0] = 0.0;
	comp[1] = 0.0;
	comp[2] = 0.0;
}
Vector3D::Vector3D(double x, double y, double z)
{
	comp[0] = x;
	comp[1] = y;
	comp[2] = z;
}
Vector3D::Vector3D(const Vector3D& v)
{
	comp[0] = v[0];
	comp[1] = v[1];
	comp[2] = v[2];
}
// Operator overloadings:
Vector3D Vector3D::operator+(Vector3D v)
{
	return Vector3D((comp[0] + v[0]), (comp[1] + v[1]), (comp[2] + v[2]));
}
Vector3D Vector3D::operator-(Vector3D v)
{
	return Vector3D((comp[0] - v[0]), (comp[1] - v[1]), (comp[2] - v[2]));
}
Vector3D operator*(const Vector3D& v, double s)
{
	return Vector3D(s * v[0], s * v[1], s * v[2]);
}
Vector3D operator*(double s, const Vector3D& v)
{
	return Vector3D(s * v[0], s * v[1], s * v[2]);
}
Vector3D Vector3D::operator/(double s)
{
	return Vector3D(comp[0] / s, comp[1] / s, comp[2] / s);
}
Vector3D Vector3D::operator=(Vector3D v)
{
	comp[0] = v[0];
	comp[1] = v[1];
	comp[2] = v[2];
	return *this;
}
Vector3D& Vector3D::operator+=(const Vector3D& v)
{
	comp[0] += v[0];
	comp[1] += v[1];
	comp[2] += v[2];
	return *this;
}
Vector3D& Vector3D::operator-=(const Vector3D& v)
{
	comp[0] -= v[0];
	comp[1] -= v[1];
	comp[2] -= v[2];
	return *this;
}
bool operator==(const Vector3D& v1, const Vector3D& v2)
{
	return ((v1[0] == v2[0]) && (v1[1] == v2[1]) && (v1[2] == v2[2]));
}
std::ostream& operator<<(std::ostream& output, const Vector3D& v)
{
	output << "(" << v[0] << " , " << v[1] << " , " << v[2] << ")";
	return output;
}
// Functions:
double Vector3D::norm()
{
	return sqrt(comp[0] * comp[0] + comp[1] * comp[1] + comp[2] * comp[2]);
}
void Vector3D::normalize()
{
	double l = norm();
	comp[0] /= l;
	comp[1] /= l;
	comp[2] /= l;
}
Vector3D Vector3D::normalized()
{
	double l = norm();
	return Vector3D(comp[0] / l, comp[1] / l, comp[2] / l);
}
double Vector3D::dot(const Vector3D& v)
{
	return comp[0] * v[0] + comp[1] * v[1] + comp[2] * v[2];
}
Vector3D Vector3D::cross(const Vector3D& v)
{
	return Vector3D((comp[1] * v[2] - comp[2] * v[1]), (comp[2] * v[0] - comp[0] * v[2]), (comp[0] * v[1] - comp[1] * v[0]));
}

// Coordinate Systems
Vector3D SphericalCoordinates(double r, double theta, double phi)
{
	Vector3D v(r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta));
	return v;
}

// 2. Event Class
// Constructors
Event::Event()
{
	time	 = 0.0;
	position = Vector3D(0.0, 0.0, 0.0);
	velocity = Vector3D(0.0, 0.0, 0.0);
}
Event::Event(double t, Vector3D& x, Vector3D& v)
{
	time	 = t;
	position = x;
	velocity = v;
}
Event::Event(double t, double x, double y, double z, double vx, double vy, double vz)
{
	time	 = t;
	position = Vector3D(x, y, z);
	velocity = Vector3D(vx, vy, vz);
}
// Set Values
void Event::SetTime(double t)
{
	time = t;
}
void Event::SetPosition(Vector3D& newposition)
{
	position = newposition;
}
void Event::SetVelocity(Vector3D& newvelocity)
{
	velocity = newvelocity;
}
void Event::SetPosition(double x, double y, double z)
{
	position = Vector3D(x, y, z);
}
void Event::SetVelocity(double vx, double vy, double vz)
{
	velocity = Vector3D(vx, vy, vz);
}
// Change values
void Event::IncreaseTime(double dt)
{
	time += dt;
}
void Event::IncreasePosition(Vector3D dx)
{
	position += dx;
}
// Return Values
double Event::Time() const
{
	return time;
}
Vector3D Event::Position() const
{
	return position;
}
Vector3D Event::Velocity() const
{
	return velocity;
}
// Returns the event in sec, km, km/sec.
Event Event::kmsec()
{
	double t   = time / sec;
	Vector3D x = 1 / km * position;
	Vector3D v = sec / km * velocity;
	return Event(t, x, v);
}
// Vector norms
double Event::Speed()
{
	return velocity.norm();
}
// Overload <<
std::ostream& operator<<(std::ostream& output, const Event& event)
{
	return output << "{"
				  << event.time
				  << ","
				  << event.position
				  << ","
				  << event.velocity
				  << "}";
}

// 3. Shielding layers
// Constructors:
Layer::Layer(std::string& ID, unsigned int ind, double rho, double thick, double deep, std::vector<std::vector<double>>& comp)
{
	// ID:
	name  = ID;
	index = ind;
	// Dimension and location:
	depth	  = deep;
	thickness = thick;
	// Content:
	density		= rho;
	composition = comp;
}
void Layer::Compute_MFP(const DM_Particle& DM, double vmin)
{
	// a) MFP does not depend on DM speed
	if(DM.ldm && !DM.screening && DM.formfactor == "Contact")
	{
		tabulated = false;
		prob.clear();
		mfp = 0.0;
		for(unsigned int i = 0; i < composition.size(); i++)
		{
			double f	   = composition[i][0];
			double Z	   = composition[i][1];
			double A	   = composition[i][2];
			double lambdai = density * f / NucleusMass(A) * DM.Sigma_Tot(Z, A, 1.0e-3);
			mfp += lambdai;
			prob.push_back(lambdai);
		}
		// Normalize probabilities
		for(unsigned int i = 0; i < prob.size(); i++)
		{
			prob[i] = prob[i] / mfp;
		}

		mfp = 1.0 / mfp;
	}

	// b) MFP does depend on DM speed
	else
	{
		tabulated = true;
		vMin	  = vmin;
		mfp_array.clear();
		prob_array.clear();

		double dv = (vesc + vEarth - vMin) / (N - 1);
		for(int i = 0; i < N; i++)
		{
			double vDM	   = vMin + i * dv;
			double mfp_aux = 0.0;
			std::vector<double> prob_aux;
			for(unsigned int j = 0; j < composition.size(); j++)
			{
				double f	   = composition[j][0];
				double Z	   = composition[j][1];
				double A	   = composition[j][2];
				double lambdai = density * f / NucleusMass(A) * DM.Sigma_Tot(Z, A, vDM);
				mfp_aux += lambdai;
				prob_aux.push_back(lambdai);
			}
			// Normalize probabilities
			for(unsigned int j = 0; j < prob_aux.size(); j++)
			{
				prob_aux[j] = prob_aux[j] / mfp_aux;
			}
			mfp_aux = 1.0 / mfp_aux;
			mfp_array.push_back(mfp_aux);
			prob_array.push_back(prob_aux);
		}
	}
}

double Layer::MFP(double vDM)
{
	if(tabulated)
	{

		int i = std::floor(1.0 * (vDM - vMin) / (vesc + vEarth - vMin) * N);
		return mfp_array[i];
	}
	else
	{
		return mfp;
	}
}
// Stopping power
double Layer::Stopping_Power(const DM_Particle& DM, double vDM) const
{
	double stopping_power = 0.0;
	for(unsigned int i = 0; i < composition.size(); i++)
	{
		double f	 = composition[i][0];
		double Z	 = composition[i][1];
		double A	 = composition[i][2];
		double mA	 = A * mNucleon;
		double nT	 = f * density / mA;
		double ERmax = 2.0 * pow(Reduced_Mass(DM.mass, mA) * vDM, 2.0) / mA;
		if(DM.ldm && DM.formfactor != "General")
		{
			double a = Thomas_Fermi_Radius(Z);
			if(DM.formfactor == "Contact")
			{
				if(!DM.screening)
					stopping_power += ERmax * nT * DM.sigmaSI(Z, A) / 2.0;
				else
				{
					stopping_power += nT * DM.sigmaSI(Z, A) / ERmax * (-1.0 + 2.0 * pow(a, 2.0) * ERmax * mA * (-2.0 + pow(a, 2.0) * ERmax * mA) + 1.0 / (1.0 + 2.0 * pow(a, 2.0) * ERmax * mA) + 3.0 * log1p(2.0 * pow(a, 2.0) * ERmax * mA)) / (4. * pow(a, 4.0) * pow(mA, 2));
				}
			}
			else if(DM.formfactor == "Electric-Dipole")
			{
				stopping_power += nT * DM.sigmaSI(Z, A) / ERmax * -(pow(qRef, 2.0) * ((-2.0 * pow(a, 2.0) * ERmax * mA * (1.0 + pow(a, 2.0) * ERmax * mA)) / (1.0 + 2.0 * pow(a, 2.0) * ERmax * mA) + log1p(2.0 * pow(a, 2.0) * ERmax * mA))) / (2.0 * pow(a, 2.0) * pow(mA, 2.0));
			}
			else if(DM.formfactor == "Long-Range")
			{
				stopping_power += nT * DM.sigmaSI(Z, A) / ERmax * (pow(qRef, 4.0) * (-1.0 + 1.0 / (1.0 + 2.0 * pow(a, 2.0) * ERmax * mA) + log1p(2.0 * pow(a, 2.0) * ERmax * mA))) / (4.0 * pow(mA, 2.0));
			}
		}
		else
		{
			std::function<double(double)> integrand = [DM, Z, A, vDM](double ER) {
				return ER * DM.dSdER(ER, Z, A, vDM);
			};
			double epsilon	= Find_Epsilon(integrand, 0.0, ERmax, 1.0e-5);
			double integral = Integrate(integrand, 0.0, ERmax, epsilon);
			stopping_power += nT * integral;
		}
	}
	return stopping_power;
}

std::vector<double> Layer::Sample_Target(double vDM, std::mt19937& PRNG)
{
	// Random number
	double xi  = ProbabilitySample(PRNG);
	double sum = 0.0;
	// Which element in the element?
	unsigned int i;
	for(i = 0; i < composition.size(); i++)
	{
		if(tabulated)
		{
			int j = std::floor(1.0 * (vDM - vMin) / (vesc + vEarth - vMin) * N);
			sum += prob_array[j][i];
		}
		else
		{
			sum += prob[i];
		}
		if(sum > xi)
		{
			break;
		}
	}
	return std::vector<double> {composition[i][1], composition[i][2]};
}

void Layer::Print_Summary()
{
	std::cout << "\tLayer " << index << ":\t\t" << name << std::endl
			  << "\tDepth [m]:\t\t" << InUnits(depth, meter) << std::endl
			  << "\tThickness [m]:\t\t" << InUnits(thickness, meter) << std::endl
			  << "\tDensity [g cm^-3]:\t" << InUnits(density, gram / cm / cm / cm) << std::endl
			  << "\tComposition:";
	std::cout << "\t\tf[\%]\tZ\t A" << std::endl;
	for(unsigned int i = 0; i < composition.size(); i++)
	{
		std::cout << "\t\t\t\t" << 100 * composition[i][0] << "\t" << composition[i][1] << "\t" << composition[i][2] << std::endl;
	}
	std::cout << std::endl;
}

// Compute mfp for all existing layers:
void Compute_All_MFP(const DM_Particle& DM, double vmin)
{
	for(unsigned i = 0; i < Layers.size(); i++)
		Layers[i].Compute_MFP(DM, vmin);
}

// Return the current layer
int Current_Layer(Vector3D& x)
{
	int output = -1;
	double z   = x[2];
	if(z > 0)
		output = 0;
	else if(z <= (-1.0) * Detector_Depth)
		output = Detector_Index;
	else
	{
		double depth = 0.0;
		for(unsigned int i = 0; i < Layers.size(); i++)
		{
			depth -= Layers[i].thickness;
			if(z > depth)
			{
				output = i + 1;
				break;
			}
		}
	}
	return output;
}
int Current_Layer(Event& event)
{
	Vector3D x = event.Position();
	int output = -1;
	double z   = x[2];
	if(z > 0)
		output = 0;
	else if(z <= (-1.0) * Detector_Depth)
		output = Detector_Index;
	else
	{
		double depth = 0.0;
		for(unsigned int i = 0; i < Layers.size(); i++)
		{
			depth -= Layers[i].thickness;
			// Particle sits on the upper boundary flying upwards
			if(z == (-1.0) * Layers[i].depth && event.Velocity()[2] > 0)
			{
				output = i;
				break;
			}
			else if(z > depth)
			{
				output = i + 1;
				break;
			}
		}
	}
	return output;
}
// Event when particle leaves the current layer
Event Leave_Layer(Event& x)
{
	Event output = x;
	// Current layer
	Vector3D posi = x.Position();
	int layer	  = Current_Layer(x);
	// Find t_Exit
	double z	 = posi[2];
	double vz	 = x.Velocity()[2];
	double zExit = (-1.0) * Layers[layer - 1].depth;
	if(vz < 0)
		zExit -= Layers[layer - 1].thickness;
	double tExit = (zExit - z) / vz;
	// Set the new time and position
	double tNew	  = x.Time() + tExit;
	Vector3D xNew = x.Position() + tExit * x.Velocity();
	// push it a tiny bit into the new layer for numerical reasons.
	xNew += (1.0e-7 * meter) * x.Velocity().normalized();
	output.SetTime(tNew);
	output.SetPosition(xNew);

	return output;
}

void Compare_StoppingPower(const DM_Particle& DM, double vDM, const std::vector<Layer>& layers, int rank)
{
	if(layers.size() > 1)
	{
		std::vector<double> layer_stopping;
		double layer_stopping_tot = 0;

		for(unsigned int i = 0; i < layers.size(); i++)
		{
			double s = layers[i].thickness * layers[i].Stopping_Power(DM, vDM);
			layer_stopping.push_back(s);
			layer_stopping_tot += s;
		}
		// Output
		if(rank == 0)
		{
			std::cout << std::endl
					  << "Relative stopping power of the layers." << std::endl;
			for(unsigned int i = 0; i < layers.size(); i++)
			{
				std::cout << "\t" << layers[i].name << ": " << Round(100.0 * layer_stopping[i] / layer_stopping_tot) << "\%" << std::endl;	 //<<"\t("<<exp(-Layers[i].Get_Thickness()/Layers[i].Get_MFP())<<")"<<std::endl;
			}
		}
	}
}
std::vector<Layer> Atmospheric_Layers(int layers, double altitude)
{
	std::vector<Layer> output;
	// Atmosphere composition
	std::vector<double> f = {0.756, 0.231, 0.013, 0.0002, 0.00001};
	std::vector<double> Z = {7.0, 8.0, 18.0, 6.0, 10.0};
	std::vector<double> A = {14.0, 16.0, 40.0, 12.0, 20.0};
	std::vector<std::vector<double>> composition;
	for(unsigned int i = 0; i < f.size(); i++)
		composition.push_back({f[i], Z[i], A[i]});
	// Density profile of the atmosphere
	std::vector<std::vector<double>> data = {{0, 1.225}, {2, 1.007}, {4, 0.8193}, {6, 0.6601}, {8, 0.5258}, {10, 0.4135}, {12, 0.3119}, {14, 0.2279}, {16, 0.1665}, {18, 0.1216}, {20, 0.08891}, {22, 0.06451}, {24, 0.04694}, {26, 0.03426}, {28, 0.02508}, {30, 0.01841}, {32, 0.01355}, {34, 0.009887}, {36, 0.007257}, {38, 0.005366}, {40, 0.003995}, {42, 0.002995}, {44, 0.002259}, {46, 0.001714}, {48, 0.001317}, {50, 0.001027}, {52, 0.0008055}, {54, 0.0006389}, {56, 0.0005044}, {58, 0.0003962}, {60, 0.0003096}, {62, 0.0002407}, {64, 0.000186}, {66, 0.0001429}, {68, 0.0001091}, {70, 0.00008281}, {72, 0.00006236}, {74, 0.00004637}, {76, 0.0000343}, {78, 0.00002523}, {80, 0.00001845}, {82, 0.00001341}, {84, 9.69e-6}, {86, 6.955e-6}};
	Interpolation density(data, km, kg / meter / meter / meter);
	// Integrated density
	double integral = Integrate(density, altitude, 86.0 * km, 1e-7);
	// Divide the atmosphere into sub-layers
	std::vector<double> boundaries = {altitude};
	std::vector<double> densities;
	std::vector<double> thickness;
	for(int i = 0; i < layers; i++)
	{
		std::function<double(double)> func = [&density, &boundaries, layers, integral](double y) {
			double eps = Find_Epsilon(density, boundaries.back(), y, 1e-6);
			return (Integrate(density, boundaries.back(), y, eps) - 1.0 * integral / layers);
		};
		boundaries.push_back((i < (layers - 1)) ? Find_Root(func, boundaries.back(), 86.0 * km, 1e-4) : 86.0 * km);
		thickness.push_back((boundaries[i + 1] - boundaries[i]));
		densities.push_back(1.0 * integral / layers / thickness.back());
		// Construct layer
		std::string ID = "Atmosphere " + std::to_string(layers - i) + "/" + std::to_string(layers);
		int index	   = layers - i;
		double depth   = 86.0 * km - boundaries.back();
		output.push_back(Layer(ID, index, densities.back(), thickness.back(), depth, composition));
	}
	std::reverse(output.begin(), output.end());
	return output;
}

// 4. Random number generation and sampling
//  Real random number which is used as a seed for the PRNG:
std::random_device rd;
// The PRNG:
std::mt19937 generator(rd());
// Propability:
std::uniform_real_distribution<double> distrXi(0, 1);
// Uniform angles:
std::uniform_real_distribution<double> distrcos(-1, 1);
std::uniform_real_distribution<double> distrphi(0.0, 2 * M_PI);
// For Rejection sampling, to find the norm of the maxwell distributed initial velocity
std::uniform_real_distribution<double> distrx(0.0, vesc);
std::uniform_real_distribution<double> distry(0.0, 1 / Nesc);

// Random numbers between 0 and 1:
double ProbabilitySample(std::mt19937& PRNG)
{
	return distrXi(PRNG);
}
// Random Angles
double ThetaSample(std::mt19937& PRNG)
{
	return acos(distrcos(PRNG));
}
double PhiSample(std::mt19937& PRNG)
{
	return distrphi(PRNG);
}

// Rejection sampling
double Rejection_Sampling(const std::function<double(double)>& PDF, double xMin, double xMax, double yMax, std::mt19937& PRNG)
{
	bool success = false;
	double x;
	int count = 0;
	while(!success)
	{
		count++;
		if(count % 1000 == 0)
		{
			cout << "Warning in Rejection_Sampling(): Very inefficient sampling with N=" << count << "." << endl;
		}
		x		   = xMin + distrXi(PRNG) * (xMax - xMin);
		double y   = distrXi(PRNG) * yMax;
		double pdf = PDF(x);
		if(pdf < 0.0)
		{
			std::cerr << "Error in Rejection_Sampling(): PDF is negative -> f(" << x << ")=" << pdf << endl;
			std::exit(EXIT_FAILURE);
		}
		else if(pdf > yMax)
		{
			std::cout << "Warning in Rejection_Sampling(): PDF>yMax, yMax is set to PDF." << endl;
			return Rejection_Sampling(PDF, xMin, xMax, pdf, PRNG);
		}
		else if(y <= pdf)
			success = true;
	}
	return x;
}

// 5. Functions for splitting/Russian roulette
bool Survive_Russian_Roulette(double pSurv, std::mt19937& PRNG)
{
	if(ProbabilitySample(PRNG) < pSurv)
		return true;
	else
		return false;
}
int Importance_Domains(const DM_Particle& DM, double vDM, double vCutoff, double kappa, std::vector<Layer>& layers)
{
	// Integrated stopping power
	double S = 0.0;
	for(unsigned int i = 0; i < layers.size(); i++)
	{
		S += layers[i].thickness * layers[i].Stopping_Power(DM, vDM);
	}
	// Energy loss from average to cutoff
	double dE		= 0.5 * DM.mass * (vDM * vDM - vCutoff * vCutoff);
	double Ndomains = ceil(kappa * S / dE);
	return (Ndomains < 1e4) ? Ndomains : 1e4;
}

std::vector<double> importance_boundaries;
std::vector<double> importances;
void Compute_Importance_Boundaries(const DM_Particle& DM, double vDM, int N_layers, double splits, std::vector<Layer>& layers)
{
	// Delete old values;
	importance_boundaries.clear();
	importances.clear();
	if(N_layers > 1 && splits > 1)
	{
		importance_boundaries.push_back(0.0);
		importances.push_back(1);
		// Compute stopping powers of layers
		std::vector<double> layer_stopping;
		double layer_stopping_tot = 0;

		for(unsigned int i = 0; i < layers.size(); i++)
		{
			double s = layers[i].thickness * layers[i].Stopping_Power(DM, vDM);
			layer_stopping.push_back(s);
			layer_stopping_tot += s;
		}

		// Divide shielding into equal-stopping layers
		double equal_stopping = 1.0 * layer_stopping_tot / N_layers;
		double z			  = 0.0;
		unsigned int l		  = 0;
		int n				  = 1;
		double aux			  = equal_stopping;
		while(n <= N_layers)
		{
			// residual stopping in this layer
			double dz = (layers[l].depth + layers[l].thickness + z);
			double ds = dz / layers[l].thickness * layer_stopping[l];
			// a) next importance boundary is in this layer
			if(ds >= aux)
			{
				z -= 1.0 * aux / ds * dz;
				importance_boundaries.push_back(z);
				importances.push_back(pow(splits, n));
				n++;
				aux = equal_stopping;
			}
			// b) no more importance boundaries in this layer
			else
			{
				aux -= ds;
				l++;
				if(l >= layers.size())
					break;
				z = -layers[l].depth;
			}
		}
	}
}
double Importance(double z)
{
	if(importances.empty())
		return 1;
	else if(z >= importance_boundaries.front())
		return 0;
	else if(z <= importance_boundaries.back())
		return importances.back();
	else
	{
		unsigned int i;
		for(i = 0; i < importance_boundaries.size(); i++)
		{
			if(z >= importance_boundaries[i])
				break;
		}
		return importances[i - 1];
	}
}