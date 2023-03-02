#include "Input_Parameters.hpp"

#include <fstream>
#include <iostream>
#include <libconfig.h++>
#include <sys/stat.h>	 //required to create a folder
#include <sys/types.h>	 // required for stat.h

#include "DD_Electron.hpp"

using namespace libconfig;

// 1. Global variables and input parameters
// Version
std::string version = "1.2.0";
// Input variables
// Simulation ID
std::string ID = "default";
// Statistical parameter
unsigned int SampleSize = 0;
// Certainty Level
double CL = -1.0;
// Importance sampling:
double IS_Angle = 0.0;	 // Scattering Angle
double IS_MFP	= 0.0;	 // MFP
// Importance splitting
bool GIS		  = false;
int GIS_max_layers = 25;
double GIS_Splits = 0.0;
double GIS_Kappa  = 0.0;

// Parameter Scan
double mMin	  = 0.0;
double mMax	  = 0.0;
double dm	  = 0.0;
int Masses	  = 0;
double dSigma = 0.0;
// Layer structure
bool Atmosphere = false;
double Altitude = 0.0;
int Atmo_Layers = 0;
std::vector<Layer> Layers;
// Dark Matter
DM_Particle DM;
// Detector
std::string Detector = "default";
// DM-Electron scattering with semiconductors
std::string DMe_target	= "default";
double DMe_exposure		= 0.0;
int DMe_threshold		= 0.0;
double DMe_efficiency	= 1.0;
unsigned int DMe_events = 0;
// Detector Location
double Detector_Depth = -1.0;
int Detector_Index	  = -1;

// 2. Read in the input parameters
void Copy_Config_File(const char* inputfile)
{
	std::ifstream inFile;
	std::ofstream outFile;
	inFile.open(inputfile);
	outFile.open("../results/" + ID + "/" + ID + ".cfg");
	outFile << inFile.rdbuf();
	inFile.close();
	outFile.close();
}

void Read_Config_File(const char* inputfile, int rank)
{
	std::string line = "----------";
	Config cfg;
	// Read the file. If there is an error, report it and exit.
	try
	{
		cfg.readFile(inputfile);
		if(rank == 0)
			cout << "Config file:\t\t" << inputfile << endl
				 << line << endl;
	}
	catch(const FileIOException& fioex)
	{
		std::cerr << "I/O error while reading configuration file." << std::endl;
		exit(EXIT_FAILURE);
	}
	catch(const ParseException& pex)
	{
		std::cerr << "Configurate file parse error at " << pex.getFile() << ":" << pex.getLine() << " - " << pex.getError() << std::endl;
		exit(EXIT_FAILURE);
	}
	// Simulation ID
	try
	{
		ID = cfg.lookup("simID").c_str();
		if(rank == 0)
			cout << "\tSimulation ID:\t\t" << ID << endl;
	}
	catch(const SettingNotFoundException& nfex)
	{
		cerr << "No 'simID' setting in configuration file." << endl;
		exit(EXIT_FAILURE);
	}
	// Create folder for results
	if(rank == 0)
	{
		cout << endl
			 << "\tCreating folder /results/" + ID << "." << endl;
		std::string sPath = "../results/" + ID;
		mode_t nMode	  = 0733;	// UNIX style permissions
		int nError		  = 0;
#if defined(_WIN32)
		nError = _mkdir(sPath.c_str());	  // can be used on Windows
#else
		nError = mkdir(sPath.c_str(), nMode);	// can be used on non-Windows
#endif
		if(nError != 0)
		{
			cout << "\tThe folder already exists, data will be overwritten." << endl;
		}
	}
	// Copy cfg file into this folder
	Copy_Config_File(inputfile);

	// Initial SampleSize
	if(rank == 0)
		cout << endl
			 << "\tStatistics:" << endl;

	try
	{
		SampleSize = cfg.lookup("samplesize");
		if(rank == 0)
			cout << "\t\tMinimal sample size:\t" << SampleSize << endl;
	}
	catch(const SettingNotFoundException& nfex)
	{
		cerr << "No 'samplesize' setting in configuration file." << endl;
		exit(EXIT_FAILURE);
	}
	// Certainty Level
	try
	{
		CL = cfg.lookup("cl");
		if(rank == 0)
			cout << "\t\tCertainty level:\t" << CL << endl;
	}
	catch(const SettingNotFoundException& nfex)
	{
		cerr << "No 'cl' setting in configuration file." << endl;
		exit(EXIT_FAILURE);
	}
	// Halo parameter
	if(rank == 0)
		cout << endl
			 << "\tDM halo parameter:" << endl;
	// DM energy density
	try
	{
		rhoDM = cfg.lookup("rhoDM");
		rhoDM *= GeV / cm / cm / cm;
		if(rank == 0)
			cout << "\t\trhoDM\t\t" << InUnits(rhoDM, GeV / cm / cm / cm) << " GeV/cm^3" << endl;
	}
	catch(const SettingNotFoundException& nfex)
	{
		cerr << "No 'rhoDM' setting in configuration file." << endl;
		exit(EXIT_FAILURE);
	}
	// Velocity dispersion
	try
	{
		v0 = cfg.lookup("v0");
		v0 *= km / sec;
		if(rank == 0)
			cout << "\t\tv0:\t\t" << InUnits(v0, km / sec) << " km/sec" << endl;
	}
	catch(const SettingNotFoundException& nfex)
	{
		cerr << "No 'v0' setting in configuration file." << endl;
		exit(EXIT_FAILURE);
	}
	// Earth velocity
	try
	{
		vEarth = cfg.lookup("vEarth");
		vEarth *= km / sec;
		if(rank == 0)
			cout << "\t\tvE:\t\t" << InUnits(vEarth, km / sec) << " km/sec" << endl;
	}
	catch(const SettingNotFoundException& nfex)
	{
		cerr << "No 'vEarth' setting in configuration file." << endl;
		exit(EXIT_FAILURE);
	}
	// Galactic escape velocity
	try
	{
		vesc = cfg.lookup("vEscape");
		vesc *= km / sec;
		if(rank == 0)
			cout << "\t\tvEscape:\t" << InUnits(vesc, km / sec) << " km/sec" << endl;
	}
	catch(const SettingNotFoundException& nfex)
	{
		cerr << "No 'vEscape' setting in configuration file." << endl;
		exit(EXIT_FAILURE);
	}

	// Variation reduction
	if(rank == 0)
		cout << endl
			 << "\tVariation reduction techniques:" << endl;
	// Importance sampling parameter for the scattering angle:
	try
	{
		IS_Angle = cfg.lookup("is_angle");
	}
	catch(const SettingNotFoundException& nfex)
	{
		cerr << "No 'is_angle' setting in configuration file." << endl;
		exit(EXIT_FAILURE);
	}
	// Importance sampling parameter for the mean free path:
	try
	{
		IS_MFP = cfg.lookup("is_mfp");
	}
	catch(const SettingNotFoundException& nfex)
	{
		cerr << "No 'is_mfp' setting in configuration file." << endl;
		exit(EXIT_FAILURE);
	}
	if(rank == 0 && (IS_Angle > 0.0 || IS_MFP > 0.0))
	{
		cout << "\t\tImportance Sampling\t[x]" << endl;
		cout << "\t\t\tIS (angle):\t" << IS_Angle << endl;
		cout << "\t\t\tIS (MFP):\t" << IS_MFP << endl;
	}
	else if(rank == 0)
		cout << "\t\tImportance Sampling\t[ ]" << endl;
	// Importance splitting layer:
	try
	{
		GIS = cfg.lookup("splitting");
	}
	catch(const SettingNotFoundException& nfex)
	{
		cerr << "No 'splitting' setting in configuration file." << endl;
		exit(EXIT_FAILURE);
	}
	// Importance splits at boundary:
	try
	{
		GIS_Splits = cfg.lookup("splits");
		if(GIS_Splits <= 1 && GIS)
		{
			cerr << "Error: Option 'splits' in configuration file has to be larger than 1 for active GIS." << endl;
			exit(EXIT_FAILURE);
		}
	}
	catch(const SettingNotFoundException& nfex)
	{
		cerr << "No 'splits' setting in configuration file." << endl;
		exit(EXIT_FAILURE);
	}
	// Coefficients
	try
	{
		GIS_Kappa = cfg.lookup("kappa");
	}
	catch(const SettingNotFoundException& nfex)
	{
		cerr << "No 'kappa' setting in configuration file." << endl;
		exit(EXIT_FAILURE);
	}
	// Max number of layers
	try
	{
		GIS_max_layers = cfg.lookup("maxLayers");
	}
	catch(const SettingNotFoundException& nfex)
	{
		// if maxLayers is not in the config, using the default value - 25
		GIS_max_layers = 25;
	}
	if(rank == 0 && GIS)
	{
		cout << "\t\tImportance Splitting\t[x]" << endl;
		cout << "\t\t\tSplits:\t\t" << GIS_Splits << endl;
		cout << "\t\t\tKappa:\t\t" << GIS_Kappa << endl;
		cout << "\t\t\tMaxLayers:\t\t" << GIS_max_layers << endl        
			 << endl;
	}
	else if(rank == 0)
		cout << "\t\tImportance Splitting\t[ ]" << endl
			 << endl;

	// Interaction parameter
	if(rank == 0)
		cout << "\tDM interactions:" << endl;
	// De-/Activate light DM option
	bool LDM;
	try
	{
		LDM = cfg.lookup("LDM");
		if(rank == 0)
			cout << "\t\tLight DM:\t\t" << ((LDM) ? "[x]" : "[ ]") << endl;
	}
	catch(const SettingNotFoundException& nfex)
	{
		cerr << "No 'LDM' setting in configuration file." << endl;
		exit(EXIT_FAILURE);
	}
	// Form factor F_DM
	std::string FF_DM;
	try
	{
		FF_DM = cfg.lookup("DM_FormFactor").c_str();
		if(rank == 0)
			cout << "\t\tForm factor F_DM:\t" << FF_DM << endl;
	}
	catch(const SettingNotFoundException& nfex)
	{
		cerr << "No 'DM_FormFactor' setting in configuration file." << endl;
		exit(EXIT_FAILURE);
	}
	// Coupling to Z or A
	std::string ZorA;
	try
	{
		ZorA = cfg.lookup("ZorA").c_str();
		if(rank == 0)
			cout << "\t\tCoupling to:\t\t" << ZorA << endl;
	}
	catch(const SettingNotFoundException& nfex)
	{
		cerr << "No 'ZorA' setting in configuration file." << endl;
		exit(EXIT_FAILURE);
	}
	// Screening
	bool Screening;
	if(FF_DM == "Electric-Dipole" || FF_DM == "Long-Range")
		Screening = true;
	else
	{
		try
		{
			Screening = cfg.lookup("Screening");
		}
		catch(const SettingNotFoundException& nfex)
		{
			cerr << "No 'Screening' setting in configuration file." << endl;
			exit(EXIT_FAILURE);
		}
	}
	if(rank == 0)
		cout << "\t\tCharge screening:\t" << ((Screening) ? "[x]" : "[ ]") << endl;
	// Mediator mass
	double mMediator;
	if(FF_DM == "General")
	{
		try
		{
			mMediator = cfg.lookup("mMediator");
			mMediator *= MeV;
			if(rank == 0)
				cout << "\tmMediator:" << mMediator << endl;
		}
		catch(const SettingNotFoundException& nfex)
		{
			cerr << "No 'mMediator' setting in configuration file." << endl;
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		mMediator = 0.0;
	}
	DM = DM_Particle(0.0, 0.0, 0.0, LDM, FF_DM, ZorA, Screening, mMediator);
	if(rank == 0)
		cout << endl;
	// Experiment
	if(rank == 0)
		cout << "\tExperiment:" << endl;
	try
	{
		Detector = cfg.lookup("experiment").c_str();
		if(rank == 0)
			cout << "\t\tDetector:\t" << Detector << endl;
	}
	catch(const SettingNotFoundException& nfex)
	{
		cerr << "No 'experiment' setting in configuration file." << endl;
		exit(EXIT_FAILURE);
	}
	// DM-Electron
	if(Detector == "Semiconductor")
	{

		try
		{
			DMe_target = cfg.lookup("target").c_str();
			if(DMe_target != "Si" && DMe_target != "Ge")
			{
				std::cerr << "Error in cfg file. Semiconductor target " << DMe_target << " not recognized." << endl;
				std::exit(EXIT_FAILURE);
			}
			// Import form factor
			else
			{
				Import_FormFactor(DMe_target);
			}
			DMe_threshold  = cfg.lookup("threshold");
			DMe_exposure   = cfg.lookup("exposure");
			DMe_efficiency = cfg.lookup("efficiency");
			DMe_events	   = cfg.lookup("events");
			if(rank == 0)
				cout << "\t\tTarget:\t\t" << DMe_target << endl
					 << "\t\tThreshold:\t" << DMe_threshold << endl
					 << "\t\tExposure[g yr]:\t" << DMe_exposure << endl
					 << "\t\tEfficiency:\t" << DMe_efficiency << endl
					 << "\t\tEvents:\t\t" << DMe_events << endl
					 << endl;
			DMe_exposure *= gram * year;
		}
		catch(const SettingNotFoundException& nfex)
		{
			cerr << "No 'Semiconductor' setting in configuration file." << endl;
			exit(EXIT_FAILURE);
		}
	}
	else if(Detector == "SENSEI" || Detector == "SENSEI-surface")
	{

		DMe_target = "Si";
		Import_FormFactor(DMe_target);	 // Import form factor
		DMe_threshold  = 1;
		DMe_exposure   = (Detector == "SENSEI-surface") ? 0.07 * gram * 456 * minute : 0.246 * gram * day;
		DMe_efficiency = 1.0;
	}
	else if(Detector == "SuperCDMS" || Detector == "DAMIC-M")
	{

		DMe_target = "Si";
		Import_FormFactor(DMe_target);	 // Import form factor
		DMe_threshold  = 1;
		DMe_exposure   = (Detector == "SuperCDMS") ? 0.487 * gram * day : 1.0 * kg * year;
		DMe_efficiency = (Detector == "SuperCDMS") ? 0.9545 : 1.0;
	}
	else if(Detector == "XENON10e" || Detector == "XENON100e")
	{
		// Import form factor
		DMe_target = "Xe";
		// Import_FFion();
		Import_FormFactor(DMe_target);
	}
	else if(Detector == "DarkSide-50")
	{
		// Import form factor
		DMe_target = "Ar";
		// Import_FFion();
		Import_FormFactor(DMe_target);
	}
	// Parameter scan
	if(rank == 0)
		cout << "\tParameter scan:" << endl;
	// mMin
	try
	{
		mMin = cfg.lookup("mMin");
		if(rank == 0)
			cout << "\t\tmMin [GeV]:\t" << mMin << endl;
	}
	catch(const SettingNotFoundException& nfex)
	{
		cerr << "No 'mMin' setting in configuration file." << endl;
		exit(EXIT_FAILURE);
	}
	// mMax
	try
	{
		mMax = cfg.lookup("mMax");
		if(rank == 0)
			cout << "\t\tmMax [GeV]:\t" << mMax << endl;
	}
	catch(const SettingNotFoundException& nfex)
	{
		cerr << "No 'mMax' setting in configuration file." << endl;
		exit(EXIT_FAILURE);
	}
	// masses
	try
	{
		Masses = cfg.lookup("masses");
		if(rank == 0)
			cout << "\t\tMass steps:\t" << Masses << endl;
	}
	catch(const SettingNotFoundException& nfex)
	{
		cerr << "No 'masses' setting in configuration file." << endl;
		exit(EXIT_FAILURE);
	}
	// masses
	try
	{
		dSigma = cfg.lookup("dSigma");
		if(rank == 0)
			cout << "\t\tCS stepsize:\t" << dSigma << endl;
	}
	catch(const SettingNotFoundException& nfex)
	{
		cerr << "No 'dSigma' setting in configuration file." << endl;
		exit(EXIT_FAILURE);
	}
	// Number of steps
	if(mMin == mMax)
		Masses = 1;
	if(Masses == 1)
		dm = 0;
	else
		dm = (log10(mMax) - log10(mMin)) / (Masses - 1.0);
	// Layer Structure
	if(rank == 0)
		cout << endl
			 << "\tShielding Layers:" << endl
			 << endl;
	// Atmosphere
	try
	{
		Atmosphere = cfg.lookup("atmosphere");
	}
	catch(const SettingNotFoundException& nfex)
	{
		cerr << "No 'atmosphere' setting in configuration file." << endl;
		exit(EXIT_FAILURE);
	}
	if(Atmosphere)
	{
		if(rank == 0)
			cout << "\tAtmosphere\t[x]" << endl;
		try
		{
			Atmo_Layers = cfg.lookup("atmo_layers");
			if(rank == 0)
				cout << "\t\tAtm. layers:\t" << Atmo_Layers << endl;
		}
		catch(const SettingNotFoundException& nfex)
		{
			cerr << "No 'atmo_layers' setting in configuration file." << endl;
			exit(EXIT_FAILURE);
		}
		try
		{
			Altitude = cfg.lookup("altitude");
			if(rank == 0)
				cout << "\t\tAltitude[m]:\t" << Altitude << endl;
			if(Altitude >= 86000.0 || Altitude < 0.0)
			{
				cerr << "Error: Altitude must be chosen between 0.0 and 86000.0." << endl;
				exit(EXIT_FAILURE);
			}
			Altitude *= meter;
		}
		catch(const SettingNotFoundException& nfex)
		{
			cerr << "No 'altitude' setting in configuration file." << endl;
			exit(EXIT_FAILURE);
		}
		Layers = Atmospheric_Layers(Atmo_Layers, Altitude);
	}
	else
	{
		if(rank == 0)
			cout << "\tAtmosphere\t[ ]" << endl;
	}
	// User defined layers
	double layer_depth	= (Layers.size() == 0) ? 0.0 : Layers.back().depth + Layers.back().thickness;
	const Setting& root = cfg.getRoot();
	try
	{
		const Setting& layers = root["layers"];
		int count			  = layers.getLength();
		if(rank == 0)
			cout << endl
				 << "\tLayers:\t\t\t" << count << endl;
		if(rank == 0)
			cout << "\t" << line << endl
				 << "\tLayer 0:\t\tSpace" << endl;
		for(int i = 0; i < count; i++)
		{
			const Setting& layer   = layers[i];
			std::string layer_name = layer.lookup("name").c_str();
			double layer_density   = layer.lookup("density");
			double layer_thickness = layer.lookup("thickness");
			if(layer_thickness == 0.0 || layer_density == 0.0)
				continue;
			// Units
			layer_density *= gram / cm / cm / cm;
			layer_thickness *= meter;
			// Composition
			std::vector<std::vector<double>> layer_composition;
			int element_count = layer.lookup("composition").getLength();
			for(int j = 0; j < element_count; j++)
			{
				std::vector<double> element;
				for(int k = 0; k < 3; k++)
					element.push_back(layer.lookup("composition")[j][k]);
				layer_composition.push_back(element);
			}
			// Construct Layer
			unsigned int layer_index = Layers.size() + 1;
			Layer next_layer(layer_name, layer_index, layer_density, layer_thickness, layer_depth, layer_composition);
			Layers.push_back(next_layer);

			// Increase depth by the layer's thickness
			layer_depth += layer_thickness;
		}
		for(unsigned int i = 0; i < Layers.size(); i++)
		{
			if(rank == 0)
			{
				cout << "\t" << line << endl;
				Layers[i].Print_Summary();
			}
		}
		if(rank == 0)
			cout << "\t" << line << endl;
	}
	catch(const SettingNotFoundException& nfex)
	{
		// Ignore.
	}
	// Derived Quantities
	// Halo
	Nesc = M_PI * v0 * v0 * (sqrt(M_PI) * v0 * erf(vesc / v0) - 2 * vesc * exp(-vesc * vesc / v0 / v0));
	// Detector
	Detector_Depth = layer_depth;
	Detector_Index = Layers.size() + 1;
	if(rank == 0)
		cout << "\tLayer " << Detector_Index << ":\t\tDetector" << endl
			 << "\tType:\t\t\t" << Detector << endl
			 << "\tDepth:\t\t\t" << InUnits(Detector_Depth, meter) << "m" << endl
			 << "\t" << line << endl
			 << line << endl;
}