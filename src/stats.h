#ifndef _STATS_h
#define _STATS_h
#include "main.h"

#define MAX_ANGLE 180
#define PPG  3 //params per gaussian. amplitude, sigma, position

using namespace FGraph;
using namespace boost;
using namespace std;
using namespace arma;



enum class ModelMatch {Exponential, PowerLaw, LogNorm};

//holds the different functions that define different statistics models
//along with whether or not they use a varying xmin, and the number of random values (and output values) are consumed when generating sample values from the distribution
template<typename T>
struct ModelHolder
{
	const std::function<T(const vector<double>&,const vector<double>::iterator &)> calculate_params_func; //generate the model parameters from the data and xmin
	const std::function<double(const double, const double, const T&)> cdf_func; //calculate the cumulative distribution function for a given x value and xmin value
	const std::function<double(const double, const double, const T&)> pdf_func; //calculate the probability distribution function for a given x value and xmin value
	const std::function<double(const double, const double, const T&)> convert_single_func; //convert a single uniform random number (0 to 1) to an element from the distribution (r, xmin, parameters)
	const std::function<vector<double>(const vector<double>&, const double, const T&)> convert_multi_func; //convert multiple uniform random numbers, to the same number of model samples (r values, xmin, parameters).
	const bool use_xmin; //whether or not this distriubtion model uses xmin
	const int  samples_per_call; //number of random values to give and receive when generating random values
};

//struct to hold the relevant information for a particular realisation of a statistics model
//including the statistics model data, and an name, for convenience
template<typename T>
struct DistStats{
	const ModelHolder<T> model; //the functions for the type of statistics model
	const std::string name; //a name for the model
	T params; //the parameters that define this particular instance of the statistics model
	double xmin; //the xmin value (will be the min of the original dataset, if the model doesn't specifically use one)
	int xmin_ind; //the index of the minimum x value (relative to the sorted original dataset)
	double ks; //best Kolmogorov-Smirnoff(v?) value
	double pvalue; //the P value for this instance of the model, as compared to the original dataset. Larger is "better", in this case. (ie, larger pvalue <-> higher prob that the distribution is compatible with the dataset)
	double sl; //sum of the log-likelihood values
	DistStats<T>(ModelHolder<T> mdl, std::string n) : model(mdl), name(n), params({0}), xmin(0), xmin_ind(0), ks(0), pvalue(0), sl(0) {};
};

struct PowerLawParams
{
	double c;
	double alpha;
	
	friend std::ostream& operator<<(std::ostream &os, const PowerLawParams &p) {return os << "c = " << p.c << ", alpha = " << p.alpha;}
};

struct ExponentialParams
{
	double lambda;
	friend std::ostream& operator<<(std::ostream &os, const ExponentialParams &p) {return os << "lambda = " << p.lambda;}
};

struct LogNormParams
{
	double mu;
	double sigma;
	friend std::ostream& operator<<(std::ostream &os, const LogNormParams &p) {return os << "mu = " << p.mu << ", sigma = " << p.sigma;}
};

struct WeibullParams
{
	double lambda, beta;
	friend std::ostream& operator<<(std::ostream &os, const WeibullParams &p) {return os << "lambda = " << p.lambda << ", beta = " << p.beta;}
};

typedef boost::variant<DistStats<PowerLawParams>, DistStats<ExponentialParams>, DistStats<LogNormParams>, DistStats<WeibullParams>> DistStatsType;

struct StatsModelData
{
	//parameters for each model
	//match with each model
	std::vector<DistStatsType> models; //the different models and their parameters
	std::vector<std::vector<std::tuple<double, double>>> comparisons; //the R and p values for the comparisons of each pair of models
	std::vector<int> best_matches; //a vector that lists all the equal best matches (as there may be more than one model that counts as a "best match")
	int best_match; //a single best match (for convenience)
};

class STATS
{
	public:
		friend class GEOMETRIE;

		STATS();
		~STATS()
		{}

	void CreateStats(VECTOR lines);
	
	template <typename T>
	void DEManalysis(Graph& G, double radius, string filename, T r);
	
	template <typename T>
	void AnalyseRaster(VECTOR lines, T R);


	double PointExtractor(point_type P, double radius, double Transform[8], double** raster);
	
	StatsModelData GetLengthDist(VECTOR lines);
	
	void DoBoxCount(VECTOR lines);


	double MinVarBuf(line_type L,  double GeoTransform[8], double** raster);
	void KDE_estimation_strikes(VECTOR lines);
	
	
	
	
	
	void KMCluster(bool input, int No);
};
#endif
