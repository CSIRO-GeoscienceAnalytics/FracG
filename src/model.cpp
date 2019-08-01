#include "model.h"

MODEL::MODEL ()

{

}

void MODEL::GaussianMix(vec data, int g)
{
	gmm_diag model;

	bool status = model.learn(data, g, maha_dist, random_subset, 10, 5, 1e-10, true);
	if(status == false)
		cout << "learning failed" << endl;
}
