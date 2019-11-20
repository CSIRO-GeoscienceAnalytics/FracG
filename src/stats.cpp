#include "stats.h"
#include "geometrie.h"
#include "GeoRef.h"
#include "graph.h"
#include "model.h"

STATS::STATS ()

{

}

vec STATS::Bootstrapping(vec Data)
{
	srand ( time(NULL) ); 
	int maxNb = Data.size(), size_seg, seg;
	float new_CV, CV =  0.9 * stddev(Data) / mean(Data);
	vec sampled_data = Data;

	do {
		size_seg = 0.05 * maxNb ;
		seg = rand() % size_seg;
		for (int i = 0; i < seg; i++)
		{
			int p1 = rand() % maxNb ;
			int p2 = rand() % maxNb ;
			sampled_data(p1) = Data(p2);
		}
		new_CV =  stddev(sampled_data) / mean(sampled_data);
	} while (new_CV > CV);
	cout << CV << endl;
	return(sampled_data);
}

//struct to hold the quad tree data
struct QuadTreeNode
{
	box bbox; //bounding box for the square
	vector<box> child_bboxes; //bboxes for this node's children
	vector<std::unique_ptr<QuadTreeNode>> children;
	int generation;
	//how do I have this be null until i set them? if i use a vector, i cant guarantee that the nodes will line up with the bounding boxes? c++ is wierd
	friend std::ostream& operator<<(std::ostream &os, const QuadTreeNode &n) {
		const point_type &min = n.bbox.min_corner();
		const point_type &max = n.bbox.max_corner();
		const double min_x = min.x();
		const double max_x = max.x();
		const double min_y = min.y();
		const double max_y = max.y();
		return os << "QuadTreeNode gen " << n.generation << " (" << min_x << ", " << min_y << ") -> (" << max_x << ", " << max_y << ")";}
};

//Make a new QuadTreeNode from its bounding box
QuadTreeNode NewQuadTreeNode(const box &bbox, const int generation = -1)
{
	QuadTreeNode node;
	node.generation = generation;
	node.bbox = bbox;
	const point_type &min = bbox.min_corner();
	const point_type &max = bbox.max_corner();
	const double min_x = min.x();
	const double max_x = max.x();
	const double min_y = min.y();
	const double max_y = max.y();
	const double half_x = (min_x + max_x)/2;
	const double half_y = (min_y + max_y)/2;
	node.child_bboxes.push_back(box(point_type( min_x,  min_y), point_type(half_x, half_y)));
	node.child_bboxes.push_back(box(point_type(half_x,  min_y), point_type( max_x, half_y)));
	node.child_bboxes.push_back(box(point_type( min_x, half_y), point_type(half_x,  max_y)));
	node.child_bboxes.push_back(box(point_type(half_x, half_y), point_type( max_x,  max_y)));
	for (unsigned int counter = 0; counter < node.child_bboxes.size(); counter++) node.children.push_back(nullptr); //why was this < an <=? did I change it? why didn't it cause problems earlier?
	return std::move(node);
}

//Check a fault segment, and add it to this node and this node's children (where they overlap with the fault segment)
void QuadTreeAddFaultSegment(QuadTreeNode &node, const line_type &fault, int layers)
{
	if (node.children.size() != 4)
	{
		cout << "ERROR - node " << node << " has children length " << node.children.size() << " at layer " << layers << endl;
		return;
	}
	const box fault_bbox = geometry::return_envelope<box>(fault);
	for (unsigned int i = 0; i < node.children.size(); i++){
		const box child_bbox = node.child_bboxes[i];
		if (geometry::disjoint(child_bbox, fault_bbox) || geometry::disjoint(child_bbox, fault)) continue; //check the fault's bounding box first, because it is faster
		std::unique_ptr<QuadTreeNode> &child = node.children[i];
		if (child == nullptr) child = std::move(std::make_unique<QuadTreeNode>(NewQuadTreeNode(child_bbox, node.generation+1)));
		if (layers < 1) continue;
// 		#pragma omp task shared(node, fault, child)
		{
			geometry::model::multi_linestring<line_type> child_faults;
			geometry::intersection(fault, node.child_bboxes[i], child_faults);
			for (auto it = child_faults.begin(); it != child_faults.end(); it++) QuadTreeAddFaultSegment(*node.children[i], *it, layers-1);
		}
	}
// 	#pragma omp taskwait
}

//Once the QuadTree has been constructed, call this to count the boxes
void CountQuadTreeBoxes(QuadTreeNode &node, std::vector<std::tuple<double, long long int>> &counts, unsigned int index)
{
	if (index >= counts.size()) //this should never be off by more than one bexo, because the preceding elements of counts should have been set when checking the parent boxes
	{
		const double side_length = (node.bbox.max_corner().x() - node.bbox.min_corner().x()); //assuming that the boxes are square
		counts.push_back(make_tuple(side_length, 0));
	}
	std::get<1>(counts[index]) += 1;//count this box
	//now check its children
	for (auto it = node.children.begin(); it != node.children.end(); it++)
	{
		if (*it != nullptr) CountQuadTreeBoxes(**it, counts, index+1);
	}
}

void STATS::BoxCountQuadTree(const vector<line_type> &faults, vector<std::tuple<double, long long int>> &counts, const double start_scale, const int max_depth, point_type &offset)
{
	box bbox(faults.front().front(), faults.front().front()); //bounding box for the entire fault set. start with an arbitrary point that is known to be inside the bounding box of the faults
	vector<box> bboxes;
	bboxes.reserve(faults.size());//bounding boxes for each fault
	for (auto it = faults.begin(); it != faults.end(); it++)
	{
		box fault_bbox = geometry::return_envelope<box>(*it);
		geometry::expand(bbox, fault_bbox);
		bboxes.push_back(fault_bbox);//std::move()
	}
	const double xOff = offset.x(), yOff = offset.y();
	const double xStart = bbox.min_corner().x() + ((xOff == 0) ? 0 : std::fmod(xOff, start_scale) - start_scale);
	const double yStart = bbox.min_corner().y() + ((yOff == 0) ? 0 : std::fmod(yOff, start_scale) - start_scale);
	const double xEnd = bbox.max_corner().x(), yEnd = bbox.max_corner().y();
	const int nX = lrint(ceil((xEnd - xStart)/start_scale)) +1; //plus one here, so that nX and nY are exclusive upper bounds. (and also measn that they're the acutal number of elements)
	const int nY = lrint(ceil((yEnd - yStart)/start_scale)) +1; //ie, the indices are in the range [0, nX) and [0, nY)
	
// 	cout << "FINAL bbox = x(" << bbox.min_corner().x()<<", "<<bbox.max_corner().x()<<"), y("<<bbox.min_corner().y()<<", "<<bbox.max_corner().y()<<")"<<endl;
// 	cout << "x = " << xStart << ", " << xEnd << ", y = " << yStart << ", " << yEnd << endl;
	
	//setup the top-most layer of the Quad Tree
	vector<vector<box>> node_bboxes; //the bounding boxes for the top-most layer of QuadTreeNode s
	vector<vector<unique_ptr<QuadTreeNode>>> nodes; //the (nodes that are) the top-most layer of QuadTreeNode s
	node_bboxes.reserve(nX);
	nodes.reserve(nX);
	for (int ix = 0; ix < nX; ix++)
	{
		vector<box> bbox_row;
		bbox_row.reserve(nY);
		vector<unique_ptr<QuadTreeNode>> node_row;
		node_row.reserve(nY);
		for (int iy = 0; iy < nY; iy++)
		{
			bbox_row.push_back(box(point_type(xStart+ix*start_scale, yStart+iy*start_scale), point_type(xStart+(ix+1)*start_scale, yStart+(iy+1)*start_scale)));
			node_row.push_back(nullptr);
		}
		node_bboxes.push_back(std::move(bbox_row));
		nodes.push_back(std::move(node_row));
	}
	
	//now make the QuadTree by adding the faults
	for (int fi = 0; fi < (int)faults.size(); fi++) //fi = "fault index"
	{
		const box &fault_bbox  = bboxes[fi];
		const line_type &fault = faults[fi];
		//only check the quadtree nodes that overlap with the fault's bounding box
		/*const */int min_x = std::lrint(std::floor((fault_bbox.min_corner().x() - xStart)/start_scale));
		/*const */int max_x = std::lrint(std::ceil ((fault_bbox.max_corner().x() - xStart)/start_scale));
		/*const */int min_y = std::lrint(std::floor((fault_bbox.min_corner().y() - yStart)/start_scale));
		/*const */int max_y = std::lrint(std::ceil ((fault_bbox.max_corner().y() - yStart)/start_scale));
		if (min_x < 0)
		{
			cout << "WARNING: box counting minimum x boundary returned " << min_x << endl;
			min_x = 0;
		}
		if (max_x >= nX)
		{
			cout << "WARNING: box counting maximum x boundary returned " << max_x << " when using " << nX << " boxes, from position (" << fault_bbox.min_corner().x() << ", " << fault_bbox.max_corner().x() << "), (" << fault_bbox.min_corner().y() << ", " << fault_bbox.max_corner().y() << ")" << endl;
			max_x = nX-1;
		}
		if (min_y < 0)
		{
			cout << "WARNING: box counting minimum y boundary returned " << min_y << endl;
			min_y = 0;
		}
		if (max_y >= nY)
		{
			cout << "WARNING: box counting maximum y boundary returned " << max_y << " when using " << nY << "boxes" << endl;
			max_y = nY-1;
		}
		//this parallelisation can be better, each fault might only be looking at a small range of the area, not enough to use all threads/cores
		//
#pragma omp parallel for collapse(2)
		for (int xind = min_x; xind <= max_x; xind++)
		{
			for (int yind = min_y; yind <= max_y; yind++)
			{
				const box &node_bbox = node_bboxes[xind][yind];
				std::unique_ptr<QuadTreeNode> &node = nodes[xind][yind];
				//if there is no overlap, check the next candidate node
				//check the simpler envelope first, then the full geometry of the fault
				if (geometry::disjoint(fault_bbox, node_bbox) || geometry::disjoint(fault, node_bbox)) continue;
				//if the node doesn't already exist, make it
				if (node == nullptr) node = std::move(std::make_unique<QuadTreeNode>(NewQuadTreeNode(node_bbox, 0)));
				//check the degments of the fault that overlap this node, and add them to the quadtree
				geometry::model::multi_linestring<line_type> fault_segments;
				geometry::intersection(fault, node_bbox, fault_segments);
				for (auto it = fault_segments.begin(); it != fault_segments.end(); it++) QuadTreeAddFaultSegment(*node, *it, max_depth);
			}
		}
	}
	
	//the quad tree has been constucted, now count them
	for (auto itx = nodes.begin(); itx != nodes.end(); itx++)
		for (auto ity = itx->begin(); ity != itx->end(); ity++)
			if (*ity != nullptr) CountQuadTreeBoxes(**ity, counts, 0);
}

void STATS::DoBoxCount(vector<line_type> faults, std::ofstream& txtF)
{
	cout << "Box Counting (QuadTree Method)" << endl;
	std::clock_t startcputime = std::clock();
	auto t_start = std::chrono::high_resolution_clock::now();
	
	//use the length of the longest fault as the starting box size
	double longest = 0;
	double shortest = std::numeric_limits<double>::infinity();
	for (auto it = faults.begin(); it != faults.end(); it++)
	{
		const double length = geometry::length(*it);
		longest = max(longest, length);
		shortest = min(shortest, length);
	}
	cout << "The fault length range is " << shortest << "m - " << longest << "m, for a range of 2^" << log2(longest/shortest) << endl;
	const int max_depth = 6;
	vector<std::tuple<double, long long int>> counts;
	point_type offset(0,0);
	BoxCountQuadTree(faults, counts, longest, max_depth, offset);
	
	//for debug purposes
// 	cout << "Box counting results:" << endl;
// 	for (auto t = counts.begin(); t < counts.end(); t++)
// 	{
// 		cout << std::get<0>(*t) << ", " << std::get<1>(*t) << endl;
// 	}
	
	auto t_end = std::chrono::high_resolution_clock::now();
	cout << "Boxcounting (QuadTree Method) 2: " << (clock() - startcputime) / (double)CLOCKS_PER_SEC << " seconds [CPU Clock] \n"
		 <<  std::chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;
	txtF <<" frequency \t size" << endl;
	
//Least-Square-Fit------------------------------------------------------
	int NB = 0;
	vec X(counts.size(),fill::zeros);
	vec Y(counts.size(),fill::zeros);
	vec XX(counts.size(),fill::zeros);
	vec XY(counts.size(),fill::zeros);
	
	for (auto it = counts.begin(); it != counts.end(); it++)
	{
		txtF << std::get<1>(*it) << "\t" << std::get<0>(*it) << endl;
		X[NB] = log(get<0>(*it));
		Y[NB] = log(get<1>(*it));
		XX[NB] = log(get<0>(*it)) * log(get<0>(*it));
		XY[NB] = log(get<0>(*it)) *  log(get<1>(*it));
		NB++;
	}
	
	double SUM_yres = 0, SUM_res = 0, R2;
	const int N = counts.size();
	double xsum  = arma::accu(X );
	double ysum  = arma::accu(Y );
	double x2sum = arma::accu(XX);
	double xysum = arma::accu(XY);

	double a = (N * xysum-xsum*ysum)/(N*x2sum-xsum*xsum);		 // slope
	double b = (x2sum*ysum-xsum*xysum)/(x2sum*N-xsum*xsum);		//intercept
	double y_interc = arma::mean(Y) - a * arma::mean(X);
	
	double *y_fit = new double[N]; 							//fitted values of y

	for (int ii = 0; ii < N; ii++)
	{
		y_fit[ii] = a*X[ii] + b;						//y(fitted) at x 
		SUM_yres += pow(Y[ii] - y_interc -( a * X[ii]),2);
		SUM_res  += pow(Y[ii] - arma::mean(Y),2);
	}
	
	R2 = (SUM_res - SUM_yres) / SUM_res;
	
	txtF <<  "No \t x(S) \t log(y(observed)) \t log(y(fitted))" << endl;
	
	for (int ii=0; ii<N; ii++)
		txtF <<"" << ii+1 << "\t"<< X[ii] << "\t" <<Y[ii] << "\t" << y_fit[ii] << endl;    

	txtF << " LS-Fit: "<< "\t" << a  << "x + " << b << endl; 
	txtF << " R^2: "   << "\t" << R2 <<endl; 
	cout << "The linear fit line is of the form: "<<a<<"x + " << b << " ( R^2 = " << R2 << " )"  << endl;

	delete[] y_fit;
}


//caluclates the Kolmogorov-Smirnov test for a given vector of values, a function that gives the cdf for the distribution being tested, the parameters for the cdf model, and the minimum x value to iterate from
//the data in values *must* be sorted (so that tha CDF for the data values is calculated in the correct order)
template<typename T>
double KSTest(vector<double> &values, const ModelHolder<T> &model, const T &params, const vector<double>::iterator &xmin)
{
// 	//for debug purposes
// 	std::ostringstream fname;
// 	fname << "kstest_" << (xmin - values.begin()) << ".txt";
// 	ofstream outF(fname.str());
// 	outF << "#x data_cdf_before data_cdf_after distribution_cdf" << endl;
	
	const int n_tail = values.end() - xmin;
// 	for (auto it = values.begin(); it != values.end(); it++) if (*it >= xmin) n_tail++; //not needed, because the data is already sorted
	double max_diff = 0; //maximum difference between the data values and the model CDF
	//int counter = 0;
	#pragma omp parallel for reduction(max:max_diff)
	for (auto it = xmin; it < values.end(); it++)
	{
		const int position = it - xmin;
		const double dist_cdf_value = model.cdf_func(*it, *xmin, params);
		const double data_cdf_before =  position      / (double)n_tail;//(++counter)
		const double data_cdf_after  = (position + 1) / (double)n_tail; //measure the difference both before and after the step function of the empirical CDF rises
		max_diff = std::max(max_diff, std::abs(data_cdf_before - dist_cdf_value));
		max_diff = std::max(max_diff, std::abs(data_cdf_after  - dist_cdf_value));
// 		outF << *it << " " << data_cdf_before << " " << data_cdf_after << " " << dist_cdf_value << endl;
	}
	
	return max_diff;
}

//Calculate the best paramaters for the distribution, for the given dataset. Iterates over each element of the values to check if it is the best option for xmin
template<typename T>
void CalculateBestDistParams(vector<double> &values, const ModelHolder<T> &model, T &params, vector<double>::iterator &xmin, double &best_ks)
{
	best_ks = std::numeric_limits<double>::infinity();
	//if we don't use xmin, then only calculate the results once
	if (!model.use_xmin)
	{
		xmin = values.begin();
		params = model.calculate_params_func(values, xmin);
		best_ks = KSTest(values, model, params, xmin);
		return;
	}
	//if we do need to select a value for xmin, try each value as a candidate
	vector<double>::iterator xmax = values.end();
	while ((xmax > values.begin()) && (xmax == values.end() || *xmax >= values.back())) --xmax; //don't use xmin here, it isn't inisialised. this is the function that is supposed to initialise it.
	#pragma omp parallel for
	for (auto it = values.begin(); it <= xmax; it++)
	{
		//first check to see if we have any values that are strictly larger than our proposed minimum
		if (*it >= values.back()) continue; //they're sorted, so if this one doesn't have enything higher than it, none of the others will either
		T this_params = model.calculate_params_func(values, it);
		const double this_ks = KSTest(values, model, this_params, it);
		#pragma omp critical
		{
			if (this_ks < best_ks)
			{
				xmin = it;
				params = std::move(this_params);
				best_ks = this_ks;
			}
		}
	}
}

//get the vector of values that are below xmin
void GetLowerValues(const vector<double> &values, const double xmin, vector<double> lower, int &n_tail)
{
	lower.clear();
	n_tail = 0;
	for (auto it = values.begin(); it != values.end(); it++) if (*it < xmin) lower.push_back(*it);
	n_tail = (int) (values.size() - lower.size());
}

//calculate/fit parameters for pdf p(x) = c * (x-xmin)**-alpha
//cdf(x, c, alpha) = (x/x_min)**-(alpha + 1)
//TODO: need to calculate optimal xmin, not just use the lowst value
//returns <alpha, c>
//default xmin = min(data)
PowerLawParams PowerLawCalculateParams(const vector<double> &data, const vector<double>::iterator &xmin)
{
	const int N = data.end() - xmin;
	double sum = 0;
	for (auto it = xmin; it != data.end(); it++)
	{
		sum += std::log(*it / *xmin); //need a gap so / * doesn't become an opening comment character
	}
	const double alpha = 1 + N/sum;
	PowerLawParams params = {(alpha - 1) / *xmin, alpha};
	return std::move(params);
}

//generate the Cumulative Density Function for a particular x value
double PowerLawCDF(const double x, const double xmin, const PowerLawParams &params)
{
	return 1 - std::pow(x/xmin, 1 - params.alpha);
}

//generate the Probability Density Function for a particular x value
double PowerLawPDF(const double x, const double xmin, const PowerLawParams &params)
{
	return params.c * std::pow(x/xmin, -params.alpha);
}

//generate an element in the distribution
double PowerLawGenerateValue(const double r, const double xmin, const PowerLawParams &params)
{
	return xmin * std::pow(r, -1/(params.alpha - 1));//1 - formerly (1 - r) //changing between (1-r) and (r) shouldn't change things, because they are both uniformly distributed between 0 and 1
}

//A "standard" way of generating (and consuming) multiple values at a time, using a converter that takes a single value at a time
template<typename T>
vector<double> StandardGenerateMultiValue(const vector<double> &r, const double xmin, const T &params, std::function<double(const double, const double, const T &params)> convert_func)
{
	vector<double> samples;
	samples.reserve(r.size());
	for (auto it = r.begin(); it != r.end(); it++) samples.push_back(convert_func(*it, xmin, params));
	return std::move(samples);
}

//calculate exponential distribution parameters for pdf p(x,lambda) = (lambda) * exp(- lambda * x)
//cdf(x, lambda) = 1 - exp[- lambda * x]
//default xmin = 0
ExponentialParams ExponentialCalculateParams(const vector<double> &data, const vector<double>::iterator &xmin)
{
	double sum = 0;
	const int count = data.end() - xmin;
	for (auto it = xmin; it != data.end(); it++)
	{
		sum += (*it); // //I think there is where to put xmin?
	}
	double lambda = 1/(sum / count - *xmin);//this is (a bit) faster
	ExponentialParams params = {lambda};
	return std::move(params);
}

double ExponentialCDF(const double x, const double xmin, const ExponentialParams &params)
{
	return 1 - exp(-params.lambda * (x - xmin));
}

double ExponentialPDF(const double x, const double xmin, const ExponentialParams &params)
{
	return params.lambda * exp(-params.lambda * (x - xmin));
}

double ExponentialGenerateValue(const double r, const double xmin, const ExponentialParams &params)
{
	return xmin - (1/params.lambda)*log(r);//1 - 
}

//calculate LogNormal distribution parameters, for pdf p(x,mu,sigma) = (1/(x sigma sqrt(2 pi))) * exp[-(ln(x) - mu)**2/(2 sigma**2)]
//cdf(x, mu, sigma) = 1/2 + (1/2) * erf[(ln(x) - mu) / (sqrt(2) sigma)]
//returns mu, sigma
//default xmin = 
LogNormParams LogNormCalculateParams(const vector<double> &data, const vector<double>::iterator &xmin)
{
	double sum = 0;
	int count = data.end() - xmin;
	for (auto it = xmin; it != data.end(); it++)
	{
		sum += log(*it);
	}
	const double mu = sum / count;
	sum = 0;
	for (auto it = xmin; it != data.end(); it++)
	{
		const double val = log(*it) - mu;
		sum += val*val;
	}
	const double sigma = sqrt(sum/count);
	LogNormParams params = {mu, sigma};
	return std::move(params);
}

double LogNormCDF(const double x, const double xmin, const LogNormParams &params)
{
	return 0.5 + 0.5*std::erf((log(x) - params.mu)/(sqrt(2)*params.sigma));
}

double LogNormPDF(const double x, const double xmin, const LogNormParams &params)
{
	const double top = log(x) - params.mu;
	const double expterm = (-top*top)/(2*params.sigma*params.sigma);
	return exp(expterm)/(x*params.sigma*sqrt(2*math::constants::pi<double>()));
}

double LogNormGenerateValue(const double r, const double xmin, const LogNormParams &params)
{
	//the paper says that there is no form for generating a single value, but this is the inverse of the LogNorm CDF function
	//exp(params.mu - 2 * params.sigma * params.sigma * inverf(1 - 2*r)); //need to find an implementation of the inverse error function, first
	return -1; //invalid operation
}

vector<double> LogNormGenerateMultiValue(const vector<double> &r, const double xmin, const LogNormParams &params)
{
	//this doesn't use mu, but then again the power law one doesn't use c either.
	const int N = (int)r.size() / 2; //making this many pairs of values
	vector<double> samples;
	samples.reserve(2*N);
	for (int i = 0; i < N; i++)
	{
		const double r1 = r[2*i  ];
		const double r2 = r[2*i+1];
		const double rho = std::sqrt(-2*params.sigma*params.sigma * std::log(r1));//1 - 
		const double theta = 2*math::constants::pi<double>()*r2;
		samples.push_back(std::exp(rho*std::sin(theta)));
		samples.push_back(std::exp(rho*std::cos(theta)));
	}
	return std::move(samples);
}

WeibullParams WeibullCalculateParams(const vector<double> &data, const vector<double>::iterator &xmin)
{
	const int N = data.end() - xmin;
	double lambda = 0, beta = 0;
	double logsum = 0;
	for (auto it = xmin; it < data.end(); it++) logsum += log(*it);
	logsum = logsum / N;
	boost::uintmax_t max_iters = 1000;
// 	boost::math::tools::eps_tolerance<double> tol(0.001);
	//need to solve beta numerically
	const double result = boost::math::tools::newton_raphson_iterate(
		[&](double b)->std::pair<double, double>{
			double tsum = 0, bsum = 0;
			double dtsum = 0, dbsum = 0;
			for (auto it = xmin; it < data.end(); it++)
			{
				const double log_val = log(*it);
				const double pow_val = std::pow(*it, b);
				 tsum += pow_val * log_val;
				 bsum += pow_val;
				dtsum += pow_val * log_val * log_val;
				dbsum += pow_val * log_val;
			}
			const double lxm = log(*xmin);
			const double  subval =   N*pow(*xmin, b);
			const double dsubval = subval*lxm;
			const double  topv =  tsum -  subval * lxm;
			const double dtopv = dtsum - dsubval * lxm;
			const double  botv =  bsum -  subval;
			const double dbotv = dbsum - dsubval;
			const double func_val = topv / botv - logsum;
			const double diff_val = (botv*dtopv - topv*dbotv) / (botv*botv); //d/dx (t(x)/b(x)) = dt(x)/db(x) = (b(x) dt(x) - t(x) db(x)) / (b(x)**2)
			return make_pair(func_val, diff_val); //I think these calculations are correct, but I'm not entirely sure
		},
		//intial guess, min, max, number of digits precision, max iterations
		2.0, 0.001, 20.0, 5, max_iters);
// 	cout << "used " << max_iters << " iterations" << endl;
	beta = result;//(result.first + result.second) / 2;
	const double xbmin = std::pow(*xmin, beta);
	for (auto it = xmin; it < data.end(); it++) lambda += std::pow(*it, beta) - xbmin;
	lambda = lambda / N;
	return WeibullParams{lambda, beta};
}

double WeibullCDF(const double x, const double xmin, const WeibullParams &params)
{
	return 1 - exp(-std::pow((x/*-xmin*/)/params.lambda, params.beta));
}

double WeibullPDF(const double x, const double xmin, const WeibullParams &params)
{
	return (params.beta / params.lambda) * std::pow((x/*-xmin*/)/params.lambda, params.beta - 1) * exp(-std::pow((x/*-xmin*/)/params.lambda, params.beta));
}

double WeibullGenerateValue(const double r, const double xmin, const WeibullParams &params)
{
	return /*xmin +*/ params.lambda * std::pow(-log(1-r), 1/params.beta);
}

template<typename T>
vector<double> GetSampleSingles(const ModelHolder<T> &model, const T &params, const int n_samples, const vector<double>::iterator &xmin, random::mt19937 &gen, const std::vector<double> &values)
{
	vector<double> samples;
	samples.reserve(n_samples);
	const int max_low = std::max(0, (int) (xmin - values.begin() - 1)); //index of the final value in the data that is below xmin. Inclusive value.
	const double p_lower = (xmin - values.begin())/values.size(); //probability with which to use the below-xmin portion of the samples 
	std::uniform_real_distribution<double> uniform01(0.0, 1.0);
	std::uniform_int_distribution<int> lower_selector(0, max_low);
	for (int i = 0; i < n_samples; i++)
	{
		double value;
		if (uniform01(gen) < p_lower){
			//get sample from lower part
			value = values[lower_selector(gen)];
		} else {
			//get from distribution
			value = model.convert_single_func(uniform01(gen), *xmin, params);
		}
		samples.push_back(value);
	}
	// if (rng < p_lower) return sample(lower);
	//return convert_func(rng(), xmin, params);
	std::sort(samples.begin(), samples.end());
	return std::move(samples);
}

template<typename T>
vector<double> GetSampleMulti(const ModelHolder<T> &model, const T &params, const int n_samples, const vector<double>::iterator &xmin, random::mt19937 &gen, const std::vector<double> &values)
{
	const int per_call = model.samples_per_call;
	const int n_generations = std::lrint(std::ceil(n_samples/per_call));
	vector<double> samples;
	samples.reserve(n_generations*per_call);
	const int max_low = std::max(0, (int) (xmin - values.begin() - 1)); //index of the final value in the data that is below xmin. Inclusive value.
	const double p_lower = (xmin - values.begin())/values.size(); //probability with which to use the below-xmin portion of the samples 
	std::uniform_real_distribution<double> uniform01(0.0, 1.0);
	std::uniform_int_distribution<int> lower_selector(0, max_low);
	vector<double> random_values;
	random_values.reserve(per_call);
	for (int i = 0; i < n_generations; i++)
	{
		if (uniform01(gen) < p_lower){
			for (int j = 0; j < per_call; j++)
			{
				//get sample from lower part
				samples.push_back(values[lower_selector(gen)]);
			}
		} else {
			//get from distribution
			random_values.clear();
			for (int j = 0; j < per_call; j++) random_values.push_back(uniform01(gen));
			vector<double> model_values = model.convert_multi_func(random_values, *xmin, params);
			samples.insert(samples.end(), model_values.begin(), model_values.end());
		}
	}
	// if (rng < p_lower) return sample(lower);
	//return convert_func(rng(), xmin, params);
	std::sort(samples.begin(), samples.end());
	return std::move(samples);
}

//visitor struct to print out the cdf of a single model
//we have to print out the different models separately, because different xmin's mean different data CDF's
struct PrintCDFSingle : public boost::static_visitor<void>//(DistStats<T> &stats)
{
	public:
		const std::vector<double> &values;
		const std::string fname;
		PrintCDFSingle(const std::vector<double> &v, const std::string fn) : values(v), fname(fn) {};
		
		template<typename T>
		void operator()(DistStats<T> &stats)
		{
			const int NV = values.size(); //total number of data values
			const int N  = NV - stats.xmin_ind; //number of data values that are above xmin
			const double ND = N; //save having to convert this multiple times
			std::string model_fname = fname+"_"+stats.name+".txt";
			std::replace(model_fname.begin(), model_fname.end(), ' ', '_');
			std::transform(model_fname.begin(), model_fname.end(), model_fname.begin(), ::tolower);
			std::ofstream fout(model_fname);
			fout << "#" << stats.name << ": " << stats.params << ", xmin = " << stats.xmin << ", xmin_ind = " << stats.xmin_ind << " of " << NV << std::endl;
			fout << "#x_value data_cdf_low data_cdf_high model_cdf" << endl;
			
			for (int i = stats.xmin_ind; i < NV; i++)
			{
				const double data_low = (i     - stats.xmin_ind) / ND;
				const double data_hi =  (i + 1 - stats.xmin_ind) / ND;
				const double model_cdf = stats.model.cdf_func(values.at(i), stats.xmin, stats.params);
				fout << values[i] << " " << data_low << " " << data_hi << " " << model_cdf << std::endl;
			}
		}
};

//print out all the model CDF's, as compared to their relative data CDF's
void PrintCDFs(const StatsModelData dists, const vector<double> values, const std::string fname)
{
	PrintCDFSingle print_visitor(values, fname);
	for (auto s : dists.models) //is is a DistStats<whatever> struct
		boost::apply_visitor(print_visitor, s);
}

template<typename T1, typename T2>
std::tuple<double, double> log_likelihood_ratio(const ModelHolder<T1> &model1, const T1 &params1, const ModelHolder<T2> &model2, const T2 params2, const std::vector<double> &values, const std::vector<double>::const_iterator &xmin1, const std::vector<double>::const_iterator &xmin2)
{
	std::vector<double>::const_iterator xmin = std::max(xmin1, xmin2); //use the range of data values where both models are valid
	const int N = values.end() - xmin;
	double l1 = 0, l2 = 0; //sum of l1 and l2
	std::vector<double> l1s(N), l2s(N); //l_i(x) = log-likelihood of x = ln(p_i(x)) //init them with enough values
	//scale the (log) likelihood values so that the pdf from the mutual xmin to the final data point still integrates to one, so scale by 1/(1-cdf(mutual xmin))
	const double log_scale1 = -log(1 - model1.cdf_func(*xmin, *xmin1, params1)); //remember -log(x) == log(1/x)
	const double log_scale2 = -log(1 - model2.cdf_func(*xmin, *xmin2, params2));
	l1s.reserve(N);
	l2s.reserve(N);
	#pragma omp parallel for reduction(+:l1,l2)
	for (auto it = xmin; it < values.end(); it++)
	{
		const int index = it - xmin;
		double l1v = log(model1.pdf_func(*it, *xmin1, params1)) + log_scale1; //and remember log(a*b) == log(a) + log(b)
		l1s[index] = l1v;
		l1 += l1v;
		double l2v = log(model2.pdf_func(*it, *xmin2, params2)) + log_scale2;
		l2s[index] = l2v;
		l2 += l2v;
	}
	const double l1m = l1/N; //calculate the mean so we can calculate the standard deviation
	const double l2m = l2/N;
	const double R = l1 - l2; //the ratio that is returned. R > 0 -> model1 is a closer match, R<0 -> model 2 is a closer match, R == 0, theyre the same
	double sigma_sum = 0;
	#pragma omp parallel for reduction(+:sigma_sum)
	for (auto it = xmin; it < values.end(); it++)
	{
		const int index = it - xmin;
		double l = (l1s[index] - l2s[index]) - (l1m - l2m);
		sigma_sum += l*l;
	}
	const double sigma = std::sqrt(sigma_sum/N);
// 	cout << "comparison: R = " << R << ", sigma = " << sigma << ", sigma_sum = " << sigma_sum <<", N = " << N << ", erfc input = " << std::abs(R)/(std::sqrt(2*N)*sigma) << endl;
	const double p = (sigma_sum <= 0.0) ? 1 : erfc(std::abs(R)/(std::sqrt(2*N)*sigma)); //lower <-> more confidence that R value is meaningful
	return std::make_tuple(R, p);//std::move()
}

//generate a number of samples from a distribution, including the values that are below xmin
//split into single and multi calls, so as to avoid extra overhead in the (more common) single value case
template<typename T>
vector<double> GetSampleValues(const ModelHolder<T> &model, const T &params, const int n_samples, const vector<double>::iterator &xmin, random::mt19937 &gen, const std::vector<double> &values)
{
	if (model.samples_per_call == 1) return std::move(GetSampleSingles<T>(model, params, n_samples, xmin, gen, values));
	return std::move(GetSampleMulti<T>(model, params, n_samples, xmin, gen, values));
}

template<typename T>
double GetPValue(const ModelHolder<T> &model, const T &params, const vector<double>::iterator &xmin, const std::vector<double> &values, const double data_ks, std::vector<random::mt19937> &gen)
{
	const double err = 0.01; //accuracy to two significant figures
	const int ntrials = std::lrint(std::ceil(0.25/(err*err))); //2500
	int npass = 0;
	#pragma omp parallel for reduction(+:npass)
	for (int i = 0; i < ntrials; i++)
	{
		random::mt19937 &local_generator = gen[omp_get_thread_num()]; //each thread uses its own RNG
		vector<double> samples = GetSampleValues<T>(model, params, (int)values.size(), xmin, local_generator, values);
		vector<double>::iterator sample_xmin;
		T sample_params;
		double sample_ks;
		CalculateBestDistParams(samples, model, sample_params, sample_xmin, sample_ks);
		if (sample_ks > data_ks) npass++;
	}
	return npass / (double) ntrials;
}

//boost visitor object, which creates the specified statistics model, including calculating the ks value and the pvalue
struct GenerateStats : public boost::static_visitor<void>//(DistStats<T> &stats)
{
	public:
		std::vector<double> &values;
		std::vector<random::mt19937> &generators;
		GenerateStats(std::vector<double> &v, std::vector<random::mt19937> &gens) : values(v), generators(gens) {};
		
		template<typename T>
		void operator()(DistStats<T> &stats)/* const*/
		{
			std::vector<double>::iterator xmin;
			CalculateBestDistParams(values, stats.model, stats.params, xmin, stats.ks);
			stats.pvalue = GetPValue(stats.model, stats.params, xmin, values, stats.ks, generators);
			stats.xmin = *xmin;
			stats.xmin_ind = xmin - values.begin();
			cout << "Model " << stats.name << " with parameters " << stats.params << " has ks value " << stats.ks << ", pvalue = " << stats.pvalue << ", with xmin = " << stats.xmin << ", at index " << stats.xmin_ind << " of " << values.size() << endl;
		}
};

//compare the models against each other, in order to determine which one(s) are the best fit
void DecideBestModel(StatsModelData &model_data, const vector<double> &values, const double ACCEPT_THRESHOLD, const double DIFFERENCE_THRESHOLD)
{
	model_data.best_match = -1;
	model_data.best_matches.clear();
	bool have_model = false;
	vector<DistStatsType> &models = model_data.models; //convenience name
	const int N = model_data.models.size(); //Number of models
	//initialise comparisons as a NxN 2D vector, initialised with R=0 (no difference), p = 1 (no certainty that there is a meaningful difference even if R is large)//= std::vector<
	model_data.comparisons = std::vector<std::vector<std::tuple<double,double>>>(N, std::vector<std::tuple<double, double>>(N));
	
	//first, check to see if any models are acceptable.
	for (auto it = models.begin(); it != models.end(); it++){
		double pvalue = boost::apply_visitor([](auto& x)->double{return x.pvalue;}, *it);
		if (pvalue >= ACCEPT_THRESHOLD)
		{
			have_model = true;
			model_data.best_match = it - models.begin();
			cout << "Model " << boost::apply_visitor([](auto &x){return x.name;}, *it) << " has not been rejected" << endl;
		}
	}
	cout << endl;
	
	//if none are acceptable, stop now. there is no point in comparing them
	if (!have_model)
	{
		cout << "All models are below the acceptance threshold" << endl;
		return;
	}
	
	//do a pairwise comparison between all the acceptable models
	for (auto it1 = models.begin(); it1 != models.end(); it1++)
	{
		if (boost::apply_visitor([&](auto &x){return x.pvalue < ACCEPT_THRESHOLD;}, *it1)) continue; //if we don't accept this model, don't compare it to other models
		int ind1 = it1 - models.begin();
		std::string name1 = boost::apply_visitor([](auto &x){return x.name;}, *it1);
		//need const_iterator because values is a const vector now. (also need const_iterator's in log_likelihood_ratio())
		const std::vector<double>::const_iterator xmin1 = std::next(values.begin(), boost::apply_visitor([](auto &x){return x.xmin_ind;}, *it1)); //we can save this one to save a bit of time
		for (auto it2 = it1; it2 != models.end(); it2++)
		{
			if (it1 == it2) continue;
			if (boost::apply_visitor([&](auto &x){return x.pvalue < ACCEPT_THRESHOLD;}, *it2)) continue;
			int ind2 = it2 - models.begin();
			std::string name2 = boost::apply_visitor([](auto &x){return x.name;}, *it2);
			std::tuple<double, double> comparison_result = boost::apply_visitor(
				[&](auto &x1, auto &x2)->std::tuple<double, double>{return log_likelihood_ratio(x1.model, x1.params, x2.model, x2.params, values, xmin1, std::next(values.begin(), x2.xmin_ind));}
				, *it1, *it2);
			double R, p;
			std::tie(R, p) = comparison_result;
			model_data.comparisons[ind1][ind2] = std::move(comparison_result);
			if (p > DIFFERENCE_THRESHOLD)
			{
				cout << "Models " << name1 << " and " << name2 << " cannot be distinguished, with R = " << R << ", and p = " << p << endl;
				continue;
			}
			std::string best_name, worse_name;
			if (R >= 0)
			{
				best_name = name1;
				worse_name = name2;
			} else {
				best_name = name2;
				worse_name = name1;
			}
			cout << "Model " << best_name << " is a better fit that model " << worse_name << ", with R = " << R << ", p = " << p << endl;
// 			if (model_data.best_match == worse_ind) model_data.best_match = best_ind;
		}
	}
	
	//now make sure that we have the absolute best
	std::vector<bool> seen(N, false); //keep track of which ones we've seen
	bool changed = true;
	bool have_loop = false;
	while (changed)
	{
		changed = false;
		if (seen[model_data.best_match]) {have_loop = true; break;} //we've seen this best_match before, we have a loop
		seen[model_data.best_match] = true;
		for (int i = 0; i < N; i++)
		{
			if (i == model_data.best_match) continue;
			int i1 = model_data.best_match, i2 = i;
			bool bfirst = i1 < i2; //whether or not b(est_match) is before i in the index list. comparisons is ordered such that the lower-valued index is first
			//and R is >0 if the first model is a better match, R is < 0 if the later model is a better match
			if (!bfirst) {i1 = i; i2 = model_data.best_match;}
			double R, p;
			std::tie(R,p) = model_data.comparisons[i1][i2];
			if (p > DIFFERENCE_THRESHOLD) continue;
			if ((bfirst && R < 0) || (!bfirst && R > 0))
			{
				model_data.best_match = i;
				changed = true;
				break;
			}
		}
	}
	
	if (have_loop)
	{
		cout << "There is a loop in the relation of the best models" << endl;
		std::fill(seen.begin(), seen.end(), false);
		int next_best = model_data.best_match;
		int curr_best = model_data.best_match;
		while (!seen[next_best])
		{
			curr_best = next_best;
			seen[curr_best] = true;
			for (int i = 0; i < N; i++)
			{
				int i1 = curr_best, i2 = i;
				if (i1 > i2) {i1 = i; i2 = model_data.best_match;}
				if (i1 == i2) model_data.best_matches.push_back(i);
				else
				{
					double R, p;
					std::tie(R, p) = model_data.comparisons[i1][i2];
					if (p > DIFFERENCE_THRESHOLD) model_data.best_matches.push_back(i);
					else if (curr_best < i && R < 0) next_best = i;
					else if (curr_best > i && R > 0) next_best = i;
				}
			}
		}
		cout << "The loop is: ";
		for (auto it = model_data.best_matches.begin(); it < model_data.best_matches.end(); it++)
			cout << ((it == model_data.best_matches.begin()) ? "" : ", ") << boost::apply_visitor([](auto &x){return x.name;}, model_data.models[*it]);
		cout << endl;
		return;
	}
	
	//now see which models are equal to the selected best match
	for (int i = 0; i < N; i++)
	{
		double R, p;
		int i1 = model_data.best_match, i2 = i;
		if (i1 > i2) {i1 = i; i2 = model_data.best_match;}
		std::tie(R, p) = model_data.comparisons[i1][i2];
		if (i == model_data.best_match || p > DIFFERENCE_THRESHOLD) model_data.best_matches.push_back(i); //add the best match, and any other models that cannot be distinguished from the best match
	}
	//done
}

StatsModelData STATS::CompareStatsModels(bool print, std::vector<double> &values, std::ofstream& txtF)
{
	StatsModelData results;
	std::sort(values.begin(), values.end());
	cout << "The values are in the range " << values.front() << ", " << values.back() << endl;
// 	for (auto it = values.begin(); it != values.end(); it++)
// 	{
// 		if (it != values.begin()) cout << ", ";
// 		cout << *it;
// 	}
// 	cout << endl;
	random::random_device rd{};
	const int NTHREADS = std::max(1, omp_get_max_threads());
	if (print)
		cout << "Using " << NTHREADS << " threads" << endl;
	std::vector<random::mt19937> rngs(NTHREADS);
	for (auto it = rngs.begin(); it < rngs.end(); it++) *it = std::move(random::mt19937{rd}); //make one random number generator per omp thread
	ModelHolder<PowerLawParams> power_law_model = {PowerLawCalculateParams, PowerLawCDF, PowerLawPDF, PowerLawGenerateValue,
		[](const vector<double>&r, const double xmin, const PowerLawParams& params) -> vector<double> {return StandardGenerateMultiValue<PowerLawParams>(r, (const double)xmin, params, PowerLawGenerateValue);}, true, 1};
	ModelHolder<PowerLawParams> power_law_all_model = {PowerLawCalculateParams, PowerLawCDF, PowerLawPDF, PowerLawGenerateValue,
		[](const vector<double>&r, const double xmin, const PowerLawParams& params) -> vector<double> {return StandardGenerateMultiValue<PowerLawParams>(r, (const double)xmin, params, PowerLawGenerateValue);}, false, 1};
		
	ModelHolder<ExponentialParams> exponential_model = {ExponentialCalculateParams, ExponentialCDF, ExponentialPDF, ExponentialGenerateValue,
		[](const vector<double>&r, const double xmin, const ExponentialParams& params) -> vector<double> {return StandardGenerateMultiValue<ExponentialParams>(r, (const double)xmin, params, ExponentialGenerateValue);}, true, 1};
	ModelHolder<ExponentialParams> exponential_all_model = {ExponentialCalculateParams, ExponentialCDF, ExponentialPDF, ExponentialGenerateValue,
		[](const vector<double>&r, const double xmin, const ExponentialParams& params) -> vector<double> {return StandardGenerateMultiValue<ExponentialParams>(r, (const double)xmin, params, ExponentialGenerateValue);}, false, 1};
		
	ModelHolder<LogNormParams> lognorm_model = {LogNormCalculateParams, LogNormCDF, LogNormPDF, LogNormGenerateValue, LogNormGenerateMultiValue, true, 2};
	ModelHolder<LogNormParams> lognorm_all_model = {LogNormCalculateParams, LogNormCDF, LogNormPDF, LogNormGenerateValue, LogNormGenerateMultiValue, false, 2};
	
	/*
	ModelHolder<WeibullParams> weibull_model = {WeibullCalculateParams, WeibullCDF, WeibullPDF, WeibullGenerateValue,
		[](const vector<double>&r, const double xmin, const WeibullParams& params) -> vector<double> {return StandardGenerateMultiValue<WeibullParams>(r, (const double)xmin, params, WeibullGenerateValue);}, true, 1};
	ModelHolder<WeibullParams> weibull_all_model = {WeibullCalculateParams, WeibullCDF, WeibullPDF, WeibullGenerateValue,
		[](const vector<double>&r, const double xmin, const WeibullParams& params) -> vector<double> {return StandardGenerateMultiValue<WeibullParams>(r, (const double)xmin, params, WeibullGenerateValue);}, false, 1};
	*/
	
	vector<DistStatsType> &models = results.models; //convenience name
	
	models.push_back(DistStats<PowerLawParams>{power_law_model, std::string{"Power Law"}});
	models.push_back(DistStats<PowerLawParams>{power_law_all_model, std::string{"Power Law All"}});
	
	models.push_back(DistStats<ExponentialParams>{exponential_model, "Exponential"});
	models.push_back(DistStats<ExponentialParams>{exponential_all_model, "Exponential All"});
	
	models.push_back(DistStats<LogNormParams>{lognorm_model, "Log Normal"});
	models.push_back(DistStats<LogNormParams>{lognorm_all_model, "Log Normal All"});
	
	//weibull is slow, because you need to iteratively solve for the parameters
// 	models.push_back(DistStats<WeibullParams>{weibull_model, "Weibull"});
// 	models.push_back(DistStats<WeibullParams>{weibull_all_model, "Weibull All"});

	GenerateStats generate_visitor{values, rngs};
	
	for (auto it = models.begin(); it != models.end(); it++)
	{
		 boost::apply_visitor(generate_visitor, *it);
	}
	
	cout << endl;
	
	PrintCDFs(results, values, "model_cdf"); //print the cdf's for testing purposes
	
	const double ACCEPT_THRESHOLD = 0.1; //models are accepted if their pvalue is *ABOVE* this threshold
	
	const double DIFFERENCE_THRESHOLD = 0.1; //comparisons are considered meaningful if the pvalue is *BELOW* this threshold
	DecideBestModel(results, values, ACCEPT_THRESHOLD, DIFFERENCE_THRESHOLD); //compare the models against each other
	
	if (results.best_matches.size() == 1)
	{
		const int best = results.best_match;
		if (best >= 0)
		{
			if (print)
				cout << "The best model is " << "\t" << boost::apply_visitor([](auto& x)->std::string{return x.name;}, results.models[best]) << endl << endl;
			txtF << "The best model is " << "\t" << boost::apply_visitor([](auto& x)->std::string{return x.name;}, results.models[best]) << endl;
			
		}
	} else if (results.best_matches.size() > 1) 
	{
		if (print)
			cout << "The " << results.best_matches.size() << "\t" << " best matches are: ";
		txtF << "The " << results.best_matches.size() << " best matches are: " << "\t";
		for (auto it = results.best_matches.begin(); it < results.best_matches.end(); it++)
		{
			if (print)
				cout << (it == results.best_matches.begin() ? "" : ", ") << boost::apply_visitor([](auto& x)->std::string{return x.name;}, results.models[*it]);
			txtF << (it == results.best_matches.begin() ? "" : ", ") << boost::apply_visitor([](auto& x)->std::string{return x.name;}, results.models[*it]);
		}
		cout << endl;
		txtF << endl;
	}
	return results;
}

double STATS::ScanLineDensity(vector<line_type> faults,  vector< std::pair<double, double>>Maximas)
{
/***********************************************************************
* Berkowitz, B. (1995). Analysis of fracture network connectivity 
* using percolation theory. 
* Mathematical Geology, 27(4), 467-483.
***********************************************************************/
	GEO georef;
	point_type rn_p;
	point_type p1, p2;
	line_type ScanLine1, ScanLine2;
	int line;// intersec1, intersec2; //commenting out these unsused variables
// 	double rX1, rX2, rY1, rY2;
	vector<pair<int, double>> TotalIntersec;
	polygon_type pl = georef.BoundingBox(GeoTransform, 10);
	vector<double> directions(100);
	random::mt19937 gen;
	random::uniform_int_distribution<> distX(GeoTransform[0], GeoTransform[0] + GeoTransform[1]*GeoTransform[6]);
	random::uniform_int_distribution<> distY(GeoTransform[3], GeoTransform[3] + GeoTransform[5]*GeoTransform[7]);
	
	for (unsigned int i =0; i < Maximas.size();i++)
	{
		for (int j =0; j < (Maximas[i].second * 100);j++)
			directions.push_back(Maximas[i].first);
	}
	

	line = 0;
	do{

		bool accept = false;
		do{
			//rn_p.set<0>(distX(gen));
			//rn_p.set<1>(distY(gen));	
			//if (geometry::within(rn_p,pl))
				accept = true;
		}while (!accept);
		
		
		

		BOOST_FOREACH(line_type F, faults)
		{

		}
		
		line = 10;
	}while (line < 10);
	

	return(0);

 }
 
double STATS::PointExtractor(point_type P, double radius, double Transform[8], double** raster)
{
	GEOMETRIE geom;
	box AOI;
	BUFFER envelop;
	point_type point;
	int maxX, minX, maxY, minY;
	double M = 0;
	vector<double> D;

	envelop = geom.DefinePointBuffer(P, radius);
	geometry::envelope(envelop, AOI);
	maxX = (geometry::get<geometry::max_corner, 0>(AOI) - Transform[0]) / Transform[1];
	minX = (geometry::get<geometry::min_corner, 0>(AOI) - Transform[0]) / Transform[1]; 
	maxY = abs((geometry::get<geometry::max_corner, 1>(AOI) - Transform[3]) / Transform[5]);
	minY = abs((geometry::get<geometry::min_corner, 1>(AOI) - Transform[3]) / Transform[5]);
	
	for (int x = minX; x < maxX; x++)
	{
		point.set<0>(Transform[0] + x * Transform[1]);	
		for (int y = maxY; y < minY; y++)
		{
			point.set<1>(Transform[3] + y * Transform[5]);
			
			if (geometry::within(point, envelop))
				D.push_back(raster[x][y]);
		}
	}
	
	if (D.size() > 0)
	{
		vec DATA(D.size(),fill::zeros);
		for(int i = 0; i < (int) D.size(); i++)
			 DATA(i) = D.at(i);
		M = arma::mean(DATA);
	}

	return(M);
}
 

 
 void STATS::AnalyseRasterFaults(vector<READ> input, vector<line_type> faults)
 {
	GEO georef;
	GEOMETRIE geom;
	BUFFER envelop;
	box AOI;
	point_type point;
	vector<BUFFER> envelopes;
	polygon_type pl;
	box bx;
	int fault_no = 1;
	string  Result, ResultsF, ResultsG, exten;
	double strike, length, sinus;
	 
	ResultsF =  "FaultData_";

	for(std::vector<READ>::size_type i = 0; i != input.size(); i++) 
	{ 
		pl = georef.BoundingBox(input[i].transform, 1);
		cout << "Analyse Raster '"<< input[i].name << "' for entire fault set" << endl;
// 		double len = (abs(input[i].transform[1]) + abs(input[i].transform[5])) / 2;
		Result = ResultsF + input[i].name + ".csv";
		ofstream txtF (Result);
		txtF << "File: " << input[i].name << endl;
			
		//create simple statisics of entire raster file
		vec rA_data(input[i].transform[6] * input[i].transform[7], fill::zeros);
		int I = 0;
		for (int x = 0; x < input[i].transform[6]; x++)
			for (int y = 0; y < input[i].transform[7]; y++)
			{
				rA_data(I) = input[i].values[x][y];
				I++;
			}
		txtF << "Mean: "	<< arma::mean(rA_data) << "\t";
		txtF << "Median: "	<< arma::median(rA_data) << "\t";
		txtF << "Stdev: "	<< arma::stddev(rA_data) << "\t";
		txtF << "Var: " 	<< arma::var(rA_data) << "\t";
		txtF << "Range: "	<< arma::range(rA_data) <<" \n" << endl;

		txtF << "No. \t Strike \t Length \t Sinuosity \t Centrod Value \t Gradient at Centre \t Mean Value \t  Mean Cross Gradient \t  Total Parallel Gradient \t Mean Parallel Gradient" << endl;
//geometric properties  of  trace --------------------------------------
		
		BOOST_FOREACH(line_type F, faults)
		{
			strike = (float)(atan2(F.front().x() - F.back().x(), F.front().y() - F.back().y())) 
					* (180 / math::constants::pi<double>());
			if (strike  < 0) 
				strike  += 180;
					
			length = (long double) geometry::length(F) / 1000;
			sinus  =  geometry::length(F) / geometry::distance(F.front(), F.back() );
				
			txtF << fault_no  << "\t" << strike << "\t" << length << "\t" << sinus << "\t";
				
//check if centroid is 
			boost::geometry::centroid(F, point);
//get value and gradient at centre of lineamnet-------------------------
			if (geometry::within(point, pl))
			{
				txtF << georef.getValue(point, input[i].transform,  input[i].values) << "\t";
				txtF << georef.CentreGradient(F, input[i].transform, input[i].values) << "\t";
			}
			else
				txtF << "- \t - \t";
				
//check if fault is within raster file----------------------------------
			if(geometry::within(F.front(), pl) && geometry::within(F.back(),pl))
			{
				txtF << georef.LineExtractor(F, input[i].transform, input[i].values) << "\t";
				
				if(georef.getElevation(F.front(), input[i].name) > georef.getElevation(F.back(), input[i].name))
					txtF << (georef.getElevation(F.front(), input[i].name) - georef.getElevation(F.back(), input[i].name))/geometry::length(F) << "\t";
				else
					txtF << (georef.getElevation(F.back(), input[i].name) - georef.getElevation(F.front(), input[i].name))/geometry::length(F)  << "\t";
				
				txtF << georef.CrossGradient(F, input[i].transform, input[i].values);
				
				GEOMETRIE::SegPoints	 <geometry::model::referring_segment<point_type>> functor2;
				
				functor2 = geometry::for_each_segment(F, functor2);

				vec MEAN_GRAD (functor2.Points.size(),fill::zeros);
				

				int ii = 0;
				for (vector<std::tuple<point_type, point_type, point_type>>::
				const_iterator I = functor2.Points.begin(); I != functor2.Points.end(); ++I) 
				{
					double grad1 = 0, grad2 = 0;
					int no = 0;
					if (geometry::within(get<0>(*I), pl) && geometry::within(get<1>(*I), pl))
					{
						if(georef.getElevation(get<0>(*I), input[i].name) > georef.getElevation(get<1>(*I), input[i].name))
							grad1 = georef.getElevation(get<0>(*I), input[i].name) - georef.getElevation(get<1>(*I), input[i].name);
						else
							grad1 = georef.getElevation(get<1>(*I), input[i].name) - georef.getElevation(get<0>(*I), input[i].name);
						no++;
					}
					if (geometry::within(get<1>(*I), pl) && geometry::within(get<2>(*I), pl))
					{
						if(georef.getElevation(get<0>(*I), input[i].name) > georef.getElevation(get<1>(*I), input[i].name))
							grad2 = georef.getElevation(get<1>(*I), input[i].name) - georef.getElevation(get<2>(*I), input[i].name);
						else
							grad2 = georef.getElevation(get<2>(*I), input[i].name) - georef.getElevation(get<1>(*I), input[i].name);
						no++;
					}
					if (no != 0)
					{
						MEAN_GRAD[ii] = (grad1+grad2)/no;
						ii++;
					}
				}
				txtF << arma::mean(MEAN_GRAD);
			}
			else 
				txtF << "- \t - \t - \t - \t" ;
		txtF << endl;
		fault_no++;
		}
		cout << "  done" << endl;
	}
 }
 
 void STATS::AnalyseRasterGraph(vector<READ> input, vector<line_type> faults, point_type source, point_type target, double Dist)
 {
	GEO geo;
	GRAPH G;

	map_vertex_type map;
	
	for(std::vector<READ>::size_type i = 0; i != input.size(); i++) 
	{
		Graph graph;
		
		//string g_name ="GRAPH_"+input[i].name+".shp";
		string g_name ="GRAPH.shp";
		const char *G_n = g_name.c_str();
		map.clear();
// 		double Dist = (input[i].transform[1]+input[i].transform[5]) / 4; //distance threshold for considering separated objects to be at the same locations //use the user-given value instead
		cout << "Using distance threshold of " << Dist << endl;
		
		ofstream txtF; 
		txtF.open ("Graph_"+input[i].name+".csv"); 
		
		cout << " Analyse Raster '"<< input[i].name <<"' for Graph" << endl;
		cout << "  Building new Graph for raster" << endl;
		G.ReadVEC4raster(graph, map, faults);// //convert the faults into a graph
		G.SplitFaults(graph, map, Dist); //split the faults in the graph into fault segments, according to the intersections of the faults
		G.RemoveSpurs(graph, map, Dist); //remove any spurs from the graph network
		cout << " Assigning weights to vertices and edges" << endl;
		G.GraphAnalysis(graph, faults, txtF, 10);
// 		G.ShortPath(graph, map, source, target, 500); //for debug purposes, to see if this gets the same path as using the other graph-reading functions
		geo.AssignValuesGraph(graph, input[i].transform, input[i].values);
		
		geo.AssignValuesAll(graph, input[i].name);
		
		geo.WriteSHP(graph, G_n); //write out the graph data as a file
		geo.WriteTxt(graph, input[i].name); //write out the graph data as a text file
		DGraph flow_graph = geo.MakeDirectedGraph(graph); //don't know if we need these yet//, input[i].transform, input[i].values
		geo.setup_maximum_flow(flow_graph);
		double max_flow = geo.maximum_flow(flow_graph, source, target);
		cout << "The maximum flow from (" << source.x() << ", " << source.y() << ") to (" << target.x() << ", " << target.y() << ") is " << max_flow << endl;
		cout << "  done" << endl;
	}
 }
 
void STATS::AddData(vector<line_type> faults, point_type source, point_type target, double Dist)
{
	GEO georef;
	GEOMETRIE geom;
	
	double old_GeoTransform[8];
	bool prompt;
	char choice;
	string name, rasFile;
	double** ADD_RASTER;	
	vector<READ>  input;

	cout<<"Additional raster files? y/n ";
	cin>>choice;

	if (choice == 'y')
	{
	std::copy(std::begin(GeoTransform), std::end(GeoTransform), std::begin(old_GeoTransform));
		do{
			prompt = false;
			cout <<"raster file ";
			cin >> name;
			rasFile = (string) name;
			georef.GetRasterProperties(rasFile, ADD_RASTER);
			
			input.push_back(READ(name, ADD_RASTER, GeoTransform));
	
			while (!prompt)
			{
				cout<<"Additional raster files? y/n ";
				cin>>choice;
				cout << endl;
			
				if (choice == 'y' || choice == 'n')	
					prompt = true;
			}
		} while(choice == 'y');
		
	AnalyseRasterFaults(input, faults);
	AnalyseRasterGraph(input, faults, source, target, Dist);
	std::copy(std::begin(old_GeoTransform), std::end(old_GeoTransform), std::begin(GeoTransform));
	}
}

FSTATS STATS::KMCluster(bool input, int No, FSTATS faultStats)
{
	bool status, prompt;
	char choice;
	int clust;
	vector <double> 	data1;
	vector <double> 	data2;
	FSTATS result;
	mat data(2, faultStats.size(),fill::zeros );
	mat means;

	for ( typename FSTATS::const_iterator it = faultStats.begin(); it != faultStats.end(); it++)
	{	 
		data1.push_back(it->first); 
		data2.push_back(it->second); 
	}
	assert (data1.size() == data2.size());
	for (int i = 0; i < (int) data1.size(); i++)
	{
			data(0,i) =  data1.at(i);
			data(1,i) =  data2.at(i);
	}

	if (input)
	{
		do {
			prompt = false;
			cout <<"Enter number of clusters ";
			std::cin >> clust;

			status = kmeans(means, data, clust, random_subset, 15, true);
		
			if(!status)
				cout << "clustering failed" << endl;
			means.print("means:");

			while (!prompt)
			{
				cout<<"Clustering ok? y/n ";
				cin>>choice;
				cout << endl;
				if (choice == 'y' || choice == 'n')	
					prompt = true;	
			}
		} while(choice == 'n');
		
		int g = 1;
		data1.clear();
		data2.clear();
		for(const auto& val : means)
		{
			if (g%2)
				data1.push_back(val);
			else 
				data2.push_back(val);
			g++;
		}
		assert (data1.size() == data2.size());
		for (int i = 0; i < (int) data1.size(); i++)
			result.push_back(make_pair(data1.at(i), data2.at(i)));
	}
	
	else if (!input)
	{
		status = kmeans(means, data, No, random_subset, 15, false);
	
		if(status == false)
			cout << "clustering failed" << endl;
		
		int g = 1;
		data1.clear();
		data2.clear();
		for(const auto& val : means)
		{
			if (g%2)
				data1.push_back(val);
			else 
				data2.push_back(val);
			g++;
		}
		assert (data1.size() == data2.size());
		for (int i = 0; i < (int) data1.size(); i++)
			result.push_back(make_pair(data1.at(i), data2.at(i)));
	}
	return result;
}

//gaussion kernel for kde
float k_gauss( double val)
{
 const double p = 1.0 / sqrt( 2.0 * M_PI);
 float result = 0.5 * (val*val);
 result = p * std::exp(- result);
 return result;
}

//simple KDE
vector< std::pair<double, double> > kde(arma::vec val, int len, float b)
 {
 vector< std::pair<double, double> > density( len );
 
 arma::vec val2(len,fill::zeros);
 int i = 0;
 for (auto& v : val)
	 val2[i] = v + 180;
	 
 vec val3 = join_cols<vec>(val, val2);
 len = val.size();
 
 const double p = 1.0 / (b * len );
	for(int r = 0; r < len; r++)
	{
		double sum = 0;
		for(int i = 0; i < len; i++)
		{
			double dist = abs(val[r] - val[i]);
			sum +=  k_gauss( dist / b );
		}
	density[r] = std::make_pair( val[r], p*sum );
	}
  return density;
 }
 
 void MA_filter(vector< std::pair<double,double>>&KDE, int window)
 {
	unsigned int wn = (int) (window - 1) / 2;
	vector<float>KDE_filtered(KDE.size());

	for (unsigned int i = 0; i < KDE.size(); i++) 
	{
		double av = 0;
		if ((i - wn) < 0)
		{
			for (unsigned int ii = KDE.size(); ii > (KDE.size() - wn); ii--)
				av += KDE[ii].second;
		}
		else
		{
			for (unsigned int ii = i; ii > (i-wn); ii--)
				av += KDE[ii].second;
		}

		if (i + wn > KDE.size())
		{
			for (unsigned int ii = 0; ii < wn; ii++)
				av += KDE[ii].second;
		}
		else
		{
			for (unsigned int iii = i; iii < (i+wn); iii++)
				av += KDE[iii].second;
		}
		KDE_filtered[i] = av / wn;
	}
	for (int i =0; i < KDE.size();i++)
	{
		KDE[i].second = KDE_filtered[i];
	}
 }
 
 
 void interpolate(vector< std::pair<double, double>>&KDE)
{
// <boost/math/interpolators/barycentric_rational.hpp>
	std::map<double, double> r;
	for (auto e : KDE)
	{
		r[e.first] = e.second;
	}
	auto x_range = boost::adaptors::keys(r);
    auto y_range = boost::adaptors::values(r);
    math::barycentric_rational<double> interpolant(x_range.begin(), x_range.end(), y_range.begin());

	KDE.clear();
     for (int i = 0; i < 180; i++)
		KDE.push_back(make_pair(i, interpolant(i)));
}

void STATS::KDE_estimation_strikes(vector<line_type> lineaments, ofstream& txtF)
 {
	int index = 0 ;
	vec ANGLE(lineaments.size(),fill::zeros);
	BOOST_FOREACH(line_type F,  lineaments)
	{ 
		double strike = (float)(atan2(F.front().x() - F.back().x(), F.front().y() - F.back().y())) 
							* (180 / math::constants::pi<double>());
		if (strike  < 0) 
			strike  += 180;
		ANGLE(index) = strike;
		index++;
	}
	 
//KDE estimation of orientation to obtain number of lineament sets------
 sort(ANGLE.begin(), ANGLE.end());
 vector< std::pair<double, double>> GAUSS =  kde( ANGLE, ANGLE.size(), 10);
 vector< std::pair<double, double>> Maximas;
 MA_filter(GAUSS, 5);
 interpolate(GAUSS);

	for (unsigned int i = 0; i < GAUSS.size(); i++)
	{
		if (i == 0)
			if (GAUSS[i].second > GAUSS[GAUSS.size()-1].second && 
				GAUSS[i].second > GAUSS[i+1].second)
					Maximas.push_back(make_pair(GAUSS[i].first, GAUSS[i].second));

		if (i > 0 && i+1 < GAUSS.size())
			if (GAUSS[i].second > GAUSS[i-1].second &&
				GAUSS[i].second > GAUSS[i+1].second)
					Maximas.push_back(make_pair(GAUSS[i].first, GAUSS[i].second));
		
		
		if (i == GAUSS.size())
			if (GAUSS[i].second > GAUSS[i-1].second && 
				GAUSS[i].second > GAUSS[0].second)
					Maximas.push_back(make_pair(GAUSS[i].first, GAUSS[i].second));
	}

	txtF << "Fault sets:" << "\t" << Maximas.size() << endl;
	txtF << "Angle \t Density \t Smoothed Density" << endl;
	for(unsigned int i =0; i < GAUSS.size(); i++)
		txtF << GAUSS[i].first << "\t " << GAUSS[i].second  <<  endl;   
	txtF << endl;
		
	txtF << "Generation \t Angle \t Probability" << endl;
	int G = 1;
	for(auto e : Maximas)
	{
		txtF << G << "\t" << e.first << "\t" << e.second << endl;
		G++;
	}	
 }


//create statisics from input file--------------------------------------
void STATS::CreateStats(ofstream& txtF, vector<line_type> faults)
{
	MODEL model;
	size_t index = 0;
	
	int FaultNumber = faults.size();
	float mean_angle, med_angle, stddev_angle, var_angle, range_angle, strike, sinuosity;
	float mean_length, med_length, stddev_length, var_length, range_length, length;
	float mean_sin, med_sin, stddev_sin, var_sin, range_sin;
	point_type centre;
// 	double bin_size; //unused
	//						index  wkt   strike  length  sinu   cdf_stri cdf_len cdf_sin nearf  distance to nearest fault
	std::vector<std::tuple<int, string, double, double, double, double, double, double, int, double, int, double>> faultInfo;
	vector<p_index> points;
	vector<pl_index> points2;
	
	GEO georef;
	FSTATS len_dir, sin_dir,len_sin, result;
	line_type multiL;
	
	std::vector<double> fault_lengths;
	fault_lengths.reserve(faults.size());

	std::vector< std::pair< int, vector<float>>> DirecDist;
	std::vector<float> Transfer;

	vec ANGLE(FaultNumber,fill::zeros);
	vec LENGTH(FaultNumber,fill::zeros);
	vec SINUOSITY(FaultNumber,fill::zeros);
	
	BOOST_FOREACH(line_type F, faults)
	{ 
		strike = (float)(atan2(F.front().x() - F.back().x(), F.front().y() - F.back().y())) 
							* (180 / math::constants::pi<double>());
		if (strike  < 0) 
			strike  += 180;
		length = (long double) geometry::length(F) / 1000;
		sinuosity =   geometry::length(F) / geometry::distance(F.front(), F.back());

		ANGLE(index) = strike;
		LENGTH(index) = length;	
		SINUOSITY(index) = sinuosity;
		
		len_dir.push_back(make_pair(length, strike));
		len_sin.push_back(make_pair(sinuosity, length));
		sin_dir.push_back(make_pair(sinuosity, strike));
		
		fault_lengths.push_back(length);

		int i = 0;
		//string coord = "LINESTRING(";
		string coord = "";
		BOOST_FOREACH(point_type p, F)
		{
			if (i < F.size())
			{
				geometry::append(p,multiL);
				coord += to_string(p.x()) + " ";
				coord += to_string(p.y()) + ", ";
			}
			else 
			{
				coord += to_string(p.x()) + " ";
				coord += to_string(p.y()) + "";
				
			}
			i++;
		}

		geometry::centroid(F, centre);
		
		points.push_back(std::make_pair(centre, index));
		points2.push_back(std::make_tuple(centre, index, geometry::length(F)));
		
		faultInfo.push_back(make_tuple(index, coord, strike, length, sinuosity, 0, 0, 0, 0, 0, 0, 0));
		index++;
	}

	//Orientation-----------------------------------------------------------
	mean_angle    = mean(ANGLE);
	med_angle     = median(ANGLE);
	stddev_angle  = stddev(ANGLE);
	var_angle     = var(ANGLE);
	range_angle   = arma::range(ANGLE);

	//Length----------------------------------------------------------------
	mean_length   = mean(LENGTH);
	med_length    = median(LENGTH);
	stddev_length = stddev(LENGTH);
	var_length    = var(LENGTH);
	range_length  = arma::range(LENGTH);

	//Sinuosity-------------------------------------------------------------
	mean_sin      = mean(SINUOSITY);
	med_sin       = median(SINUOSITY);
	stddev_sin    = stddev(SINUOSITY);
	var_sin       = var(SINUOSITY);
	range_sin     = arma::range(SINUOSITY);

	//correlation between geometric parameters------------------------------

	mat a = conv_to<mat>::from(LENGTH);
	mat b = conv_to<mat>::from(ANGLE);

	mat cor_length_strike = cor(LENGTH, ANGLE);
	mat cor_length_sinus  = cor(LENGTH, SINUOSITY);
	mat cor_strike_sinus  = cor(ANGLE, SINUOSITY);

	//covariance between geometric parameters-------------------------------
	mat cov_length_strike = cov(LENGTH, ANGLE);
	mat cov_length_sinus  = cov(LENGTH, SINUOSITY);
	mat cov_strike_sinus  = cov(ANGLE, SINUOSITY);

	//fault centre to fault centre distance---------------------------------
	vector<p_index> closest;
	georef.Point_Tree(points, closest);
	vec distance(points.size(),fill::zeros);
	
	//fault centre to fault centre distance of larger fault----------------
	/***********************************************************************
	* Bour, O., & Davy, P. (1999). 
	* Clustering and size distributions of fault patterns: Theory and measurements. 
	* Geophysical Research Letters, 26(13), 2001-2004.
	***********************************************************************/
	vector<pl_index> closest2;
	georef.Point_Tree2(points2, closest2, floor(LENGTH.max()*1000));
	vec distance2(points.size(),fill::zeros);

	size_t curPos = 0;
	for  (std::tuple<int, string, double, double, double, double, double, double, int, double, int, double> &it : faultInfo)
	{
		get<5>(it) = normcdf((double)ANGLE(curPos), (double)mean_angle, (double)stddev_angle);
		get<6>(it) = normcdf((double)LENGTH(curPos), (double)mean_length, (double)stddev_length);
		get<7>(it) = normcdf((double)SINUOSITY(curPos), (double)mean_sin, (double)stddev_sin);
		get<8>(it) = closest[curPos].second;
		get<9>(it) = geometry::distance(points[curPos].first, closest[curPos].first);
		get<10>(it) = get<1>(closest2[curPos]);
		get<11>(it) = geometry::distance(points[curPos].first, get<0>(closest2[curPos]));
		distance[curPos] = geometry::distance(points[curPos].first, closest[curPos].first);
		distance2[curPos] = geometry::distance(points[curPos].first, get<0>(closest2[curPos]));
		curPos++;
	}
	
	txtF.open ("Fault_Statistics.csv", ios::out | ios::app); 
	txtF << "index  \t  wkt  \t  strike  \t  length  \t sinuosity  \t cdf_strike  \t cdf_length \t cdf_sinuosity \t nearest fault \t distance to nearest fault \t neares larger fault \t distance to nearest larger fault" << endl;
	
	for  (std::tuple<int, string, double, double, double, double, double, double, int, double, int, double> it : faultInfo)
	{
		txtF << get<0>(it) << "\t"<< get<1>(it) << "\t"<< get<2>(it) << "\t"<< get<3>(it) << "\t"
			 << get<4>(it) << "\t"<< get<5>(it) << "\t"<< get<6>(it) << "\t"<< get<7>(it) << "\t"
			 << get<8>(it) << "\t"<< get<9>(it) << "\t"<< get<10>(it)<< "\t"<< get<11>(it) << endl;
	}

	txtF << "Mean:     "<< "\t" << " " << "\t" << mean_angle <<"\t" << mean_length << "\t" << mean_sin << "\t" << "\t" << "\t" << "\t" << "\t" << mean(distance)   	<< "\t" << "\t" << mean(distance2)   	<< "\n"
		<< "Median:   "<< "\t" << " " << "\t" << med_angle  <<"\t" << med_length  << "\t" << med_sin  << "\t" << "\t" << "\t" << "\t" << "\t" << median(distance) 	<< "\t" << "\t" << median(distance2) 	<< "\n"
		<< "Stddev:   "<< "\t" << " " << "\t" <<stddev_angle<<"\t" <<stddev_length<< "\t" <<stddev_sin<< "\t" << "\t" << "\t" << "\t" << "\t" << stddev(distance) 	<< "\t" << "\t" << stddev(distance2) 	<< "\n"
		<< "Variance: "<< "\t" << " " << "\t" << var_angle  <<"\t" << var_length  << "\t" << var_sin  << "\t" << "\t" << "\t" << "\t" << "\t" << var(distance)		<< "\t" << "\t" << var(distance2)	 	<< "\n"
		<< "Range:    "<< "\t" << " " << "\t" << range_angle<< "\t"<< range_length<< "\t" << range_sin<< "\t" << "\t" << "\t" << "\t" << "\t" << arma::range(distance)<< "\t" << "\t" <<arma::range(distance2)<< "\n"
		<< endl;;
		
		
//KDE estimation of orientation to obtain number of lineament sets------
 sort(ANGLE.begin(), ANGLE.end());
 vector< std::pair<double, double>> GAUSS =  kde( ANGLE, ANGLE.size(), 10);
 vector< std::pair<double, double>> Maximas;
	
 MA_filter(GAUSS, 5);
 //interpolate(GAUSS);

	for (unsigned int i = 0; i < GAUSS.size(); i++)
	{
		if (i == 0)
			if (GAUSS[i].second > GAUSS[GAUSS.size()-1].second && 
				GAUSS[i].second > GAUSS[i+1].second)
					Maximas.push_back(make_pair(GAUSS[i].first, GAUSS[i].second));

		if (i != 0 && i != GAUSS.size())
			if (GAUSS[i].second > GAUSS[i-1].second &&
				GAUSS[i].second > GAUSS[i+1].second)
					Maximas.push_back(make_pair(GAUSS[i].first, GAUSS[i].second));
		
		
		if (i == GAUSS.size())
			if (GAUSS[i].second > GAUSS[i-1].second && 
				GAUSS[i].second > GAUSS[0].second)
					Maximas.push_back(make_pair(GAUSS[i].first, GAUSS[i].second));
	}

	cout << "Estimated Lineament sets: "<< Maximas.size() << endl;


	//model.GaussianMix(GAUSS, Maximas.size());

/*
	unsigned int m = 0;
	vector< std::pair<double, double> >::iterator vs = GAUSS.begin();
	
	for (unsigned int i = 0; i < GAUSS.size(); i++)
	{
		//cout << GAUSS[i].first << " " << Maximas[m].first << endl; //this takes up a lot of screen space
		if (GAUSS[i].first == Maximas[m].first)
		{
			auto result = *std::min_element(vs, GAUSS.begin() + i, [](const auto& lhs, const auto& rhs) 
			{
				return lhs.second < rhs.second;
			});
			cout << result.first << " " << result.second << endl;
			vs = GAUSS.begin() + i;
			m++;
			
			if (m > Maximas.size())
				m = 0;
		}
	}
*/
	//ScanLineDensity(faults, Maximas);

	txtF << "Fault sets:" << "\t" << Maximas.size() << endl;
	txtF << "Angle \t Density \t Smoothed Density" << endl;
	for(unsigned int i =0; i < GAUSS.size(); i++)
		txtF << GAUSS[i].first << "\t " << GAUSS[i].second  <<  endl;   
	txtF << endl;
		
	txtF << "Generation \t Angle \t Probability" << endl;
	int G = 1;
	for(auto e : Maximas)
	{
		txtF << G << "\t" << e.first << "\t" << e.second << endl;
		G++;
	}	

	txtF << "K-MEANS CLUSTERING " << endl;
	//K-MEANS Clustering----------------------------------------------------
	txtF <<" k-means length and strike (" << Maximas.size() << ") sets: " << endl;
	result = KMCluster(false, Maximas.size(), len_dir);
	for ( typename FSTATS::const_iterator it = result.begin(); it != result.end(); it++)
		txtF << it->first << "\t" << it->second << endl;
	result.clear();
	txtF << endl; 
	
	txtF <<" k-means sinuosity and strike (" << Maximas.size() << ") sets" << endl;
	result = KMCluster(false, Maximas.size(), sin_dir);
	for ( typename FSTATS::const_iterator it = result.begin(); it != result.end(); it++)
		txtF << it->first << "\t" << it->second << endl;
	result.clear();
	txtF << endl; 
	
	txtF <<" k-means sinuosity and length (" << Maximas.size() << ") sets" << endl;;
	result = KMCluster(false, Maximas.size(), len_sin);
	for ( typename FSTATS::const_iterator it = result.begin(); it != result.end(); it++)
		txtF << it->first << "\t" << it->second << endl;
	result.clear();
	
	txtF << "BOX-COUNTING " << endl;
	DoBoxCount(faults, txtF);
	txtF << "KS-FITTING (LENGTH DISTRIBUTION) " << endl;

	CompareStatsModels(true, fault_lengths, txtF);
}

void STATS::DEManalysis(Graph& G, double radius, string filename, double** RASTER)
{
	std::clock_t startcputime = std::clock();
	GEO ref;
	GEOMETRIE geom;
	BUFFER Radius;
	box AOI;
	vertex_type U, u;
	point_type point;
	double m_front = 0, m_back = 0;
	vector<double>data;
	vector<double> slope;
	double Sl;

	ofstream Points("Points.txt");

	string f = filename.substr(filename.find("/")+1); 
	filename =  f.substr(0, f.find(".")).c_str();
	string output =  (string)"DEM-" + filename.c_str() + "_vertexData.csv";
	string output2 = (string)"DEM-" + filename.c_str() + "_EdgeData.csv";
	ofstream txtF (output);
	ofstream txtF2 (output2);
	txtF << "Vertex Data" << endl; 
	txtF << "Type \t Component \t Value[m]" << endl;
	
	txtF2 << "Edge Data" << endl; 
	txtF2 << "Edge No. \t Edge Type \t Component \t Mean Value1[m] \t Stdev1 \t Mean Value2[m] \t Stdev2 \t Slope" << endl;
	
	cout << "DEM ANALYSIS " << endl; 
	//Global vertex analysis-----------------------------------------------
	if (txtF.is_open())  
	{ 
		for (auto vd : boost::make_iterator_range(vertices(G))) 
		{
			if (!G[vd].Enode)
			{
				Radius = geom.DefinePointBuffer(G[vd].location, radius);
				geometry::envelope(Radius, AOI);
				int maxX = (geometry::get<geometry::max_corner, 0>(AOI) - GeoTransform[0]) / GeoTransform[1];
				int minX = (geometry::get<geometry::min_corner, 0>(AOI) - GeoTransform[0]) / GeoTransform[1]; 
				int maxY = abs((geometry::get<geometry::max_corner, 1>(AOI) - GeoTransform[3]) / GeoTransform[5]);
				int minY = abs((geometry::get<geometry::min_corner, 1>(AOI) - GeoTransform[3]) / GeoTransform[5]);
				
				minX = std::max(minX, 0);
				maxX = std::min(maxX, (int)GeoTransform[6]);
				minY = std::max(minY, 0);
				maxY = std::min(maxY, (int)GeoTransform[7]);
				
				for (int x = minX; x < maxX; x++)
				{
					point.set<0>(GeoTransform[0] + x * GeoTransform[1]);
					for (int y = maxY; y < minY; y++)
					{
						point.set<1>(GeoTransform[3] + y * GeoTransform[5]);
					
						if (geometry::within(point, Radius))
						{
							if (degree(vd, G) == 1)
								txtF << "I" <<"\t"<< G[vd].component << "\t" << RASTER[x][y] << endl;
							if (degree(vd, G) == 3)
								txtF << "Y" <<"\t"<< G[vd].component << "\t" << RASTER[x][y] << endl;
							if (degree(vd, G) == 4)
								txtF << "X" <<"\t"<< G[vd].component << "\t" << RASTER[x][y] << endl;
						}
						geometry::clear(point);
					}
				}
			}
		}
		txtF.close();
	}
	cout << " analysed vertices " << endl;
	//IndividualFaults------------------------------------------------------
	if (txtF2.is_open())  
	{ 
		int Eg = 1;
		for (auto ve : boost::make_iterator_range(edges(G))) 
		{
			U = source(ve, G);
			u = target(ve, G);
			if (!G[U].Enode && !G[u].Enode)
			{
				GEOMETRIE::Perpencicular <geometry::model::referring_segment<point_type>> functor;
				functor = geometry::for_each_segment(G[ve].trace, functor);

				int segm = 0;
				BOOST_FOREACH(line_type cross, functor.all)
				{
					bool have_front = false, have_back = false;
					txtF2 << Eg << "\t" << G[ve].BranchType << "\t"<< G[ve].component << "\t";
					point_type f = cross.front();
					point_type b = cross.back();
//front-----------------------------------------------------------------
					Radius = geom.DefinePointBuffer(f, radius);
					geometry::envelope(Radius, AOI);
					
					int maxX = (geometry::get<geometry::max_corner, 0>(AOI) - GeoTransform[0]) / GeoTransform[1];
					int minX = (geometry::get<geometry::min_corner, 0>(AOI) - GeoTransform[0]) / GeoTransform[1]; 
					int maxY = abs((geometry::get<geometry::max_corner, 1>(AOI) - GeoTransform[3]) / GeoTransform[5]);
					int minY = abs((geometry::get<geometry::min_corner, 1>(AOI) - GeoTransform[3]) / GeoTransform[5]);
					
					//cap the values to what is available in the raster file
					minX = std::max(minX, 0);
					maxX = std::min(maxX, (int)GeoTransform[6]);
					minY = std::max(minY, 0);
					maxY = std::min(maxY, (int)GeoTransform[7]);

					for (int x = minX; x < maxX; x++)
					{
						point.set<0>(GeoTransform[0] + x * GeoTransform[1]);
						
						for (int y = maxY; y < minY; y++)
						{
							point.set<1>(GeoTransform[3] + y * GeoTransform[5]);

							if (geometry::within(point, Radius))
							{
								Points <<  GeoTransform[0]+ x *  GeoTransform[1] << "\t" <<  GeoTransform[3]+ y *  GeoTransform[5] <<  "\t" << 0 << "\t"<< RASTER[x][y] <<  endl;;
								data.push_back(RASTER[x][y]);
							}
							geometry::clear(point);
						}
					}
					
					if (data.size() !=0)
					{
						have_front = true;
						vec DATA(data.size(),fill::zeros);
						for(int i =0; i < (int) data.size(); i++)
							DATA(i) = data.at(i);
						m_front = arma::mean(DATA);
						txtF2 << m_front << "\t" << arma::stddev(DATA) << "\t";
						data.clear();
					}
						
	//back------------------------------------------------------------------
					Radius = geom.DefinePointBuffer(b, radius);
					geometry::envelope(Radius, AOI);
					maxX = (geometry::get<geometry::max_corner, 0>(AOI) - GeoTransform[0]) / GeoTransform[1];
					minX = (geometry::get<geometry::min_corner, 0>(AOI) - GeoTransform[0]) / GeoTransform[1]; 
					maxY = abs((geometry::get<geometry::max_corner, 1>(AOI) - GeoTransform[3]) / GeoTransform[5]);
					minY = abs((geometry::get<geometry::min_corner, 1>(AOI) - GeoTransform[3]) / GeoTransform[5]);
					
					minX = std::max(minX, 0);
					maxX = std::min(maxX, (int)GeoTransform[6]);
					minY = std::max(minY, 0);
					maxY = std::min(maxY, (int)GeoTransform[7]);

					for (int x = minX; x < maxX; x++)
					{
						point.set<0>(GeoTransform[0] + x * GeoTransform[1]);
						for (int y = maxY; y < minY; y++)
						{
							point.set<1>(GeoTransform[3] + y * GeoTransform[5]);

							if (geometry::within(point, Radius))
							{
								Points <<  GeoTransform[0]+ x *  GeoTransform[1] << "\t" <<  GeoTransform[3]+ y *  GeoTransform[5] << "\t" << 1 << "\t"<< RASTER[x][y] << endl;
								data.push_back(RASTER[x][y]);								
							}
							
							geometry::clear(point);
						}
					}
					if (data.size() !=0)
					{
						have_back = true;
						vec DATA = conv_to<vec>::from(data);
						m_back = arma::mean(DATA);
						txtF2 << m_back << "\t" << arma::stddev(DATA) <<"\t";
						data.clear();
					}
					
	//slope-----------------------------------------------------------------
					if (have_front && have_back)
					{
						if (m_front > m_back) //NOTE: Uli: what happens if m_front == m_back?
							Sl = (m_front - m_back) / geometry::distance(f,b);
					
						else if (m_front < m_back)
							Sl = (m_back - m_front) / geometry::distance(f,b);
						txtF2 << Sl;
						segm++;
						slope.push_back(Sl);
					}
					txtF2 << endl;
				}
			}
			if (slope.size() > 0)
			{	
				vec DATA = conv_to<vec>::from(slope);
				G[ve].CrossGrad = arma::mean(DATA);
			}
			Eg++;
		}
		txtF2.close();
	}
	cout << " analysed edges " << endl;
	//----------------------------------------------------------------------
	cout << "DEM Analysis: " << (clock() - startcputime) / (double)CLOCKS_PER_SEC << " seconds [CPU Clock] \n" << endl;
}
