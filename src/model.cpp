#include "model.h"

MODEL::MODEL ()

{

}

void MODEL::GradientDecent(FSTATS data)
{
 
	
}

void MODEL::Covariance()
{
running_stat_vec<vec> stats;
vec sample;

for(uword i=0; i<10000; ++i)
  {
  sample = randu<vec>(5);
  stats(sample);
  }

cout << "mean = " << endl << stats.mean() << endl;
cout << "var  = " << endl << stats.var()  << endl;
cout << "max  = " << endl << stats.max()  << endl;

//----------------------------------------------------------------------
running_stat_vec<rowvec> more_stats(true);

for(uword i=0; i<20; ++i)
  {
  sample = randu<rowvec>(3);
  
  sample(1) -= sample(0);
  sample(2) += sample(1);
  
  more_stats(sample);
  }

cout << "covariance matrix = " << endl;
cout << more_stats.cov() << endl;

rowvec sd = more_stats.stddev();

cout << "correlations = " << endl;
cout << more_stats.cov() / (sd.t() * sd);
}





void MODEL::GaussianMix()
{
// create synthetic data with 2 Gaussians
uword d = 5;       // dimensionality
uword N = 10000;   // number of vectors
mat data(d, N, fill::zeros);
vec mean0 = linspace<vec>(1,d,d);
vec mean1 = mean0 + 2;
uword i = 0;


  while(i < N)
  {
	if(i < N)  { data.col(i) = mean0 + randn<vec>(d); ++i; }
	if(i < N)  { data.col(i) = mean0 + randn<vec>(d); ++i; }
	if(i < N)  { data.col(i) = mean1 + randn<vec>(d); ++i; }
  }

// model the data as a diagonal GMM with 2 Gaussians
gmm_diag model;
bool status = model.learn(data, 2, maha_dist, random_subset, 10, 5, 1e-10, true);

if(status == false)
  {
	cout << "learning failed" << endl;
  }

model.means.print("means:");

double  scalar_likelihood = model.log_p( data.col(0)    );
rowvec     set_likelihood = model.log_p( data.cols(0,9) );

double overall_likelihood = model.avg_log_p(data);

uword   gaus_id  = model.assign( data.col(0),    eucl_dist );
urowvec gaus_ids = model.assign( data.cols(0,9), prob_dist );

urowvec hist1 = model.raw_hist (data, prob_dist);
 rowvec hist2 = model.norm_hist(data, eucl_dist);

model.save("my_model.gmm");
	
}
