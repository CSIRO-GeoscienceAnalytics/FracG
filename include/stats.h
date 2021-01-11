#ifndef _STATS_h
#define _STATS_h
#include "../include/fracg.h"

#define MAX_ANGLE 180
#define PPG  3 //params per gaussian. amplitude, sigma, position

namespace FracG
{

    enum class ModelMatch {Exponential, PowerLaw, LogNorm};

    //holds the different functions that define different statistics models
    //along with whether or not they use a varying xmin, and the number of random values (and output values) are consumed when generating sample values from the distribution
    template<typename T>
    struct ModelHolder
    {
            const std::function<T(const std::vector<double>&,const std::vector<double>::iterator &)> calculate_params_func; //generate the model parameters from the data and xmin
            const std::function<double(const double, const double, const T&)> cdf_func; //calculate the cumulative distribution function for a given x value and xmin value
            const std::function<double(const double, const double, const T&)> pdf_func; //calculate the probability distribution function for a given x value and xmin value
            const std::function<double(const double, const double, const T&)> convert_single_func; //convert a single uniform random number (0 to 1) to an element from the distribution (r, xmin, parameters)
            const std::function<std::vector<double>(const std::vector<double>&, const double, const T&)> convert_multi_func; //convert multiple uniform random numbers, to the same number of model samples (r values, xmin, parameters).
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


    typedef std::tuple<double, double, double> gauss_param; //the parameters are amplitude, width/sigma, and finally angle
    typedef std::vector<gauss_param> gauss_params; //these represent a sum of gaussians


    //this class represents a probability distribution of angles/orientations
    //it has one or more gaussians, and potentially a uniform distribution
    class AngleDistribution
    {
    private:
        //evaluate a single Gaussian.
        //have this as separate function, to only have to correct the angle when necessary
        static double EvaluateSingleGaussian(const gauss_param &gauss, double corrected_angle)
        {
            const double amp   = std::get<0>(gauss);
                    const double size  = std::get<1>(gauss);
                    const double pos   = std::get<2>(gauss);
                    const double dist1 = corrected_angle - pos;
                    const double dist2 = dist1 - std::copysign(180, dist1); //two differences in angle, to accound for how the angle values wrap around at 180 degrees
                    const double z1    = dist1 / size;
                    const double z2    = dist2 / size;
                    const double amp_eff = amp / (std::sqrt(2*M_PI) * size);

            const double sum1 = amp_eff * std::exp(-0.5 * z1 * z1);
            const double sum2 = amp_eff * std::exp(-0.5 * z2 * z2);

            return sum1 + sum2;
        }


    public:

        bool with_uniform = false; //true iff we are using a uniform distribution component
        double uniform_prob = 0; //the probability of the uniform components
        gauss_params gaussians; //the set of gaussians

        AngleDistribution() {};
        AngleDistribution(const gauss_params &gausses) : gaussians(gausses) {};
        AngleDistribution(const double unif, const gauss_params gausses) : with_uniform(true), uniform_prob(unif), gaussians(gausses) {};

        static double CorrectAngle(const double angle)
        {
            return fmod(fmod(angle, MAX_ANGLE) + MAX_ANGLE, MAX_ANGLE); //angle is now between 0 and 180 degrees
        }

        //evaluate a single gaussian, at a particular angle
        static double EvaluateGaussian(const gauss_param &gauss, const double angle)
        {
            const double corrected_angle = CorrectAngle(angle);
            return EvaluateSingleGaussian(gauss, corrected_angle);
        }

        //get the pdf value of the distribution at a particular angle
        double EvaluateDistribution(const double angle)
        {
            const double corrected_angle = CorrectAngle(angle);
            double sum = 0;
            for (auto it = gaussians.begin(); it < gaussians.end(); it++)
            {
                sum += EvaluateSingleGaussian(*it, corrected_angle);
            }
            if (with_uniform) sum += uniform_prob / MAX_ANGLE;

            return sum;
        }
    };

    struct StatsModelData
    {
            //parameters for each model
            //match with each model
            std::vector<DistStatsType> models; //the different models and their parameters
            std::vector<std::vector<std::tuple<double, double>>> comparisons; //the R and p values for the comparisons of each pair of models
            std::vector<int> best_matches; //a vector that lists all the equal best matches (as there may be more than one model that counts as a "best match")
            int best_match; //a single best match (for convenience)
    };

	void CreateStats(VECTOR &lines, AngleDistribution &angle_dist);
	double PointExtractor(point_type P, double radius, double Transform[8], double** raster);
	StatsModelData GetLengthDist(VECTOR lines);
	void DoBoxCount(VECTOR lines);
	double MinVarBuf(line_type L,  double GeoTransform[8], double** raster);
	AngleDistribution KdeEstimationStrikes(VECTOR &lines, const double param_penalty = 2);
	int CheckGaussians(AngleDistribution &angle_dist, double angle);
	void ScanLine(VECTOR &lines, int nb_scanlines, AngleDistribution &angle_dist, double min_spaceing);
	void RasterStatistics(VECTOR lines, double dist, std::string raster_filename);
};
#endif
