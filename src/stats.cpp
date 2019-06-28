#include "stats.h"
#include "geometrie.h"
#include "GeoRef.h"

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
	//how do I have this be null until i set them? if i use a vector, i cant guarantee that the nodes will line up with the bounding boxes? c++ is wierd
};

//Make a new QuadTreeNode from its bounding box
QuadTreeNode NewQuadTreeNode(const box &bbox)
{
	QuadTreeNode node;
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
	for (unsigned int counter = 0; counter <= node.child_bboxes.size(); counter++) node.children.push_back(nullptr);
	return std::move(node);
}

//Check a fault segment, and add it to this node and this node's children (where they overlap with the fault segment)
void QuadTreeAddFaultSegment(QuadTreeNode &node, line_type &fault, int layers)
{
	const box fault_bbox = geometry::return_envelope<box>(fault);
	for (unsigned int i = 0; i < node.children.size(); i++){
		const box &child_bbox = node.child_bboxes[i];
		if (geometry::disjoint(child_bbox, fault_bbox) || geometry::disjoint(child_bbox, fault)) continue; //check the fault's bounding box first, because it is faster
		std::unique_ptr<QuadTreeNode> &child = node.children[i];
		if (child == nullptr) child = std::move(std::make_unique<QuadTreeNode>(NewQuadTreeNode(child_bbox)));
		if (layers < 1) continue;
		geometry::model::multi_linestring<line_type> child_faults;
		geometry::intersection(fault, child_bbox, child_faults);
		for (auto it = child_faults.begin(); it != child_faults.end(); it++) QuadTreeAddFaultSegment(*child, *it, layers-1);
	}
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
#pragma omp parallel for
		for (int xind = min_x; xind <= max_x; xind++)
		{
			vector<box> &row_bboxes = node_bboxes[xind];
			vector<unique_ptr<QuadTreeNode>> &row_nodes = nodes[xind];
			for (int yind = min_y; yind <= max_y; yind++)
			{
				const box &node_bbox = row_bboxes[yind];
				std::unique_ptr<QuadTreeNode> &node = row_nodes[yind];
				//if there is no overlap, check the next candidate node
				//check the simpler envelope first, then the full geometry of the fault
				if (geometry::disjoint(fault_bbox, node_bbox) || geometry::disjoint(fault, node_bbox)) continue;
				//if the node doesn't already exist, make it
				if (node == nullptr) node = std::move(std::make_unique<QuadTreeNode>(NewQuadTreeNode(node_bbox)));
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

void STATS::DoBoxCount(vector<line_type> &faults)
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
	const int max_depth = 10;
	vector<std::tuple<double, long long int>> counts;
	point_type offset(0,0);
	BoxCountQuadTree(faults, counts, longest, max_depth, offset);
	
	auto t_end = std::chrono::high_resolution_clock::now();
	cout << "Boxcounting (QuadTree Method) 2: " << (clock() - startcputime) / (double)CLOCKS_PER_SEC << " seconds [CPU Clock] \n"
		 <<  std::chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;
	ofstream txtF ("BoxCounting.csv");
	txtF <<"frequency \t size" << endl;
	for (auto it = counts.begin(); it != counts.end(); it++)
		txtF << std::get<1>(*it) << "\t" << std::get<0>(*it) << endl;
}

void STATS::BoxCount(vector<line_type> faults, double &D)
{
	cout << "Box Counting" << endl;
	ofstream txtF ("BoxCounting.csv");
	txtF <<"frequency \t size" << endl;
	
	std::clock_t startcputime = std::clock();
	auto t_start = std::chrono::high_resolution_clock::now();

	GEOMETRIE geom;
	GEO georef;
	int I = 0;
	int i, j, d = 0;
	point_type Ul, Ur, Lr, Ll;
	vector<double> boxes;
	double L;
	seg_tree SG = georef.Seg_Tree(faults);
	vector <BUFFER> Boxes;
	vector <BUFFER> empty;
	vector <pair<int, double>> results;
	point_type box_centre;
	int count;
	int curBox =0;

	vec LENGTH(faults.size(),fill::zeros);
	BOOST_FOREACH(line_type l, faults)
	{
		LENGTH(I) = geometry::length(l);
		I++;
	}

	Ul.set<0>(GeoTransform[0]);
	Ul.set<1>(GeoTransform[3]);
	Ur.set<0>(GeoTransform[0] + GeoTransform[6] * GeoTransform[1]);
	Ur.set<1>(GeoTransform[3]);
	Ll.set<0>(GeoTransform[0]);
	Ll.set<1>(GeoTransform[3] + GeoTransform[7] * GeoTransform[5]);
	
	double lenX = geometry::distance(Ul, Ur) ;
	double lenY = geometry::distance(Ul, Ll) ;

	if (lenX > lenY )
		L = lenX;
	else
		L = lenY;

	do
	{
		L /= 2;
		boxes.push_back(L);
		results.push_back(make_pair(0,L));
	}
	while (L > 150);

	for(auto const& size : boxes) 
	{
		i = (int)ceil(lenX / size);
		j = (int)ceil(lenY / size);
		
		double X = GeoTransform[0] + size/2;
		double Y = GeoTransform[3] - size/2;

		for (int ii = 0; ii < i; ii++)
		{
			point_type xBox;
			box_centre.set<0>(X + ii*size);
			for (int jj = 0; jj < j; jj++)
			{
				box_centre.set<1>(Y - jj*size);
				Boxes.push_back(geom.DefineSquareBuffer(box_centre, size/2));
			}
		}
		
		count = 0;
		#pragma omp parallel for
		for (unsigned int i =0; i < Boxes.size(); i++)
		{
			BUFFER g = Boxes.at(i);
			if( SG.qbegin(!geometry::index::disjoint(g)) != SG.qend())
				count++;
		}
		results.at(curBox).first = count;
		curBox++;
		Boxes.clear();
	}
	
	auto t_end = std::chrono::high_resolution_clock::now();
	cout << "Boxcounting2: " << (clock() - startcputime) / (double)CLOCKS_PER_SEC << " seconds [CPU Clock] \n"
		 <<  std::chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << endl;
		
	for(auto it : results)
		txtF << it.first << "\t" << it.second << endl; 
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

template<typename T1, typename T2>
std::tuple<double, double> log_likelihood_ratio(const ModelHolder<T1> &model1, const T1 &params1, const ModelHolder<T2> &model2, const T2 params2, const std::vector<double> &values, const std::vector<double>::const_iterator &xmin1, const std::vector<double>::const_iterator &xmin2)
{
	std::vector<double>::const_iterator xmin = std::max(xmin1, xmin2); //use the range of data values where both models are valid
	const int N = values.end() - xmin;
	double l1 = 0, l2 = 0; //sum of l1 and l2
	std::vector<double> l1s(N), l2s(N); //l_i(x) = log-likelihood of x = ln(p_i(x)) //init them with enough values
	l1s.reserve(N);
	l2s.reserve(N);
	#pragma omp parallel for reduction(+:l1,l2)
	for (auto it = xmin; it < values.end(); it++)
	{
		const int index = it - xmin;
		double l1v = log(model1.pdf_func(*it, *xmin1, params1));
		l1 += l1v;
		l1s[index] = l1v;
		double l2v = log(model2.pdf_func(*it, *xmin2, params2));
		l2 += l2v;
		l2s[index] = l2v;
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

StatsModelData STATS::CompareStatsModels(std::vector<double> &values)
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
	cout << "Using " << NTHREADS << " threads" << endl;
	std::vector<random::mt19937> rngs(NTHREADS);
	for (auto it = rngs.begin(); it < rngs.end(); it++) *it = std::move(random::mt19937{rd}); //make one random number generator per omp thread
	ModelHolder<PowerLawParams> power_law_model = {PowerLawCalculateParams, PowerLawCDF, PowerLawPDF, PowerLawGenerateValue,
		[](const vector<double>&r, const double xmin, const PowerLawParams& params) -> vector<double> {return StandardGenerateMultiValue<PowerLawParams>(r, (const double)xmin, params, PowerLawGenerateValue);}, true, 1};
	ModelHolder<ExponentialParams> exponential_model = {ExponentialCalculateParams, ExponentialCDF, ExponentialPDF, ExponentialGenerateValue,
		[](const vector<double>&r, const double xmin, const ExponentialParams& params) -> vector<double> {return StandardGenerateMultiValue<ExponentialParams>(r, (const double)xmin, params, ExponentialGenerateValue);}, true, 1};
	ModelHolder<ExponentialParams> exponential_all_model = {ExponentialCalculateParams, ExponentialCDF, ExponentialPDF, ExponentialGenerateValue,
		[](const vector<double>&r, const double xmin, const ExponentialParams& params) -> vector<double> {return StandardGenerateMultiValue<ExponentialParams>(r, (const double)xmin, params, ExponentialGenerateValue);}, false, 1};
	ModelHolder<LogNormParams> lognorm_model = {LogNormCalculateParams, LogNormCDF, LogNormPDF, LogNormGenerateValue, LogNormGenerateMultiValue, true, 2};
	ModelHolder<LogNormParams> lognorm_all_model = {LogNormCalculateParams, LogNormCDF, LogNormPDF, LogNormGenerateValue, LogNormGenerateMultiValue, false, 2};

	vector<DistStatsType> &models = results.models; //convenience name
	
	models.push_back(DistStats<PowerLawParams>{power_law_model, std::string{"Power Law"}});
	
	models.push_back(DistStats<ExponentialParams>{exponential_model, "Exponential"});
	models.push_back(DistStats<ExponentialParams>{exponential_all_model, "Exponential All"});
	
	models.push_back(DistStats<LogNormParams>{lognorm_model, "Log Normal"});
	models.push_back(DistStats<LogNormParams>{lognorm_all_model, "Log Normal All"});

	GenerateStats generate_visitor{values, rngs};
	
	for (auto it = models.begin(); it != models.end(); it++)
	{
		 boost::apply_visitor(generate_visitor, *it);
	}
	cout << endl;
	
	const double ACCEPT_THRESHOLD = 0.1; //models are accepted if their pvalue is *ABOVE* this threshold
	
	const double DIFFERENCE_THRESHOLD = 0.1; //comparisons are considered meaningful if the pvalue is *BELOW* this threshold
	DecideBestModel(results, values, ACCEPT_THRESHOLD, DIFFERENCE_THRESHOLD); //compare the models against each other
	
	const int best = results.best_match;
	if (best >= 0)
		cout << "The best model is " << boost::apply_visitor([](auto& x)->std::string{return x.name;}, results.models[best]) << endl << endl;
	return results;
}

double STATS::ScanLineDensity(vector<line_type> faults)
{
/***********************************************************************
* Berkowitz, B. (1995). Analysis of fracture network connectivity 
* using percolation theory. 
* Mathematical Geology, 27(4), 467-483.
***********************************************************************/
	random::mt19937 gen;
	point_type p1, p2;
	line_type ScanLine1, ScanLine2;
	int line, intersec1, intersec2;
	double rX1, rX2, rY1, rY2, length;
	vector<pair<int, double>> TotalIntersec;
	
	double X1 = GeoTransform[0];
	double Y1 = GeoTransform[3];
	double X2 = GeoTransform[0] + GeoTransform[6] * GeoTransform[1];
	double Y2 = GeoTransform[3] + GeoTransform[7] * GeoTransform[5];

	line = 0;
	do{
		intersec1 = 0, intersec2 = 0;
		random::uniform_int_distribution<> distX(0, GeoTransform[6]);
		random::uniform_int_distribution<> distY(0, GeoTransform[7]);
		
		rX1 = X1 + (distX(gen) * GeoTransform[1]);  
		rX2 = X1 + (distX(gen) * GeoTransform[1]); 
		rY1 = Y1 + (distY(gen) * GeoTransform[5]);
		rY2 = Y1 + (distY(gen) * GeoTransform[5]);
		
		p1.set<0>(rX1);
		p1.set<1>(Y1);
		p2.set<0>(rX2);
		p2.set<1>(Y2);
		geometry::append(ScanLine1,p1);
		geometry::append(ScanLine1,p2);
		
		p1.set<0>(X1);
		p1.set<1>(rY1);
		p2.set<0>(X2);
		p2.set<1>(rY2);
		geometry::append(ScanLine2,p1);
		geometry::append(ScanLine2,p2);
		
		
		BOOST_FOREACH(line_type F, faults)
		{
			if (geometry::intersects(F, ScanLine1))
				intersec1++;
			
			if (geometry::intersects(ScanLine2, ScanLine2))
				intersec2++;
		}
		TotalIntersec.push_back(make_pair(intersec1, geometry::length(ScanLine1)));
		TotalIntersec.push_back(make_pair(intersec2, geometry::length(ScanLine2)));
		line++;
		geometry::clear(ScanLine1);
		geometry::clear(ScanLine2);
	}while (line < 10);
		
	int I = 0;
	vec SCAN (TotalIntersec.size(),fill::zeros);
	for ( typename vector<pair<int, double>>::const_iterator it = TotalIntersec.begin(); it != TotalIntersec.end(); it++)
		SCAN(I) = (it->first/(it->second / 1000));

	return( mean(SCAN) );

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
 
double STATS::LineExtractor(line_type L, double radius, double Transform[8], double** raster)
{
	GEOMETRIE geom;
	box AOI;
	BUFFER envelop;
	point_type point;
	int maxX, minX, maxY, minY;
	double M = 0;
	vector<double> D;

	envelop = geom.DefineLineBuffer(L, radius);
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
 
void STATS::AdditionalData(vector<line_type> faults, Graph& G)
{
	double old_GeoTransform[8];
	std::copy(std::begin(GeoTransform), std::end(GeoTransform), std::begin(old_GeoTransform));

	GEO georef;
	GEOMETRIE geom;
	BUFFER envelop;
	box AOI;
	point_type point;
	vector<BUFFER> envelopes;
	
	vertex_type U, u;
	long double Xmax, Xmin, Ymax, Ymin;
	polygon_type pl;
	box bx;
	int fault_no = 1, X, Y;
	bool prompt, accept, NEW;
	double variance = 1e38;
	double radius;
	char choice;
	string name, rasFile, Result, ResultsF, ResultsG, exten;
	double** ADD_RASTER;
	double strike, length, sinus;
	
	struct READ
	{
		string name;
		double** values;
		double transform[8];
		READ(string n, double** v, double t[8])
		{
			name	= n;
			values	= v;
			memcpy(transform, t, 8*sizeof(double));
		}
	};
	
	vector<READ>  input;
	vector<double> data;

	cout<<"Additional raster files? y/n ";
	cin>>choice;

	if (choice == 'y')
	{
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
	}

	ResultsF =  "FaultData_";
	
	for(std::vector<READ>::size_type i = 0; i != input.size(); i++) 
	{ 
		
		Result = ResultsF + input[i].name + ".txt";
		ofstream txtF (Result);
		
		ofstream Points("Points.txt");
		
		//create new box of raster
		const double raster_crop_size = 1000;
		Xmax = input[i].transform[0] + input[i].transform[1] * input[i].transform[6] - raster_crop_size;       // west
		Xmin = input[i].transform[0] + raster_crop_size; 									 				  // east
		Ymax = input[i].transform[3] - raster_crop_size; 													 // north
		Ymin = input[i].transform[3] + input[i].transform[5] * input[i].transform[7] + raster_crop_size;	// south
		bx.min_corner().set<0>( Xmin );
		bx.min_corner().set<1>( Ymin );
		bx.max_corner().set<0>( Xmax );
		bx.max_corner().set<1>( Ymax );
		geometry::convert(bx, pl);
			
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
		txtF << "Range: "	<< arma::range(rA_data) << endl;;


		//Extraction of all the data around the fault traces--------------------
		BOOST_FOREACH(line_type F, faults)
		{
			if(geometry::within(F.front(), pl) && geometry::within(F.back(),pl))
			{
				accept = false;
				radius = 2 * input[i].transform[1];

				strike = (float)(atan2(F.front().x() - F.back().x(), F.front().y() - F.back().y())) 
								* (180 / math::constants::pi<double>());
				if (strike  < 0) 
					strike  += 180;
				length = (long double) geometry::length(F) / 1000;
				sinus  =  geometry::length(F) / geometry::distance(F.front(), F.back() );

				do{
					data.clear();
					envelop = geom.DefineLineBuffer(F, radius);
					geometry::envelope(envelop, AOI);
					int maxX = (geometry::get<geometry::max_corner, 0>(AOI) - input[i].transform[0]) / input[i].transform[1];
					int minX = (geometry::get<geometry::min_corner, 0>(AOI) - input[i].transform[0]) / input[i].transform[1]; 
					int maxY = abs((geometry::get<geometry::max_corner, 1>(AOI) - input[i].transform[3]) / input[i].transform[5]);
					int minY = abs((geometry::get<geometry::min_corner, 1>(AOI) - input[i].transform[3]) / input[i].transform[5]);

					minX = std::max(minX, 0);
					maxX = std::min(maxX, (int)input[i].transform[6]);
					minY = std::max(minY, 0);
					maxY = std::min(maxY, (int)input[i].transform[7]);
					
					for (int x = minX; x < maxX; x++)
					{
						point.set<0>(input[i].transform[0] + x * input[i].transform[1]);
					
						for (int y = maxY; y < minY; y++)
						{
							point.set<1>(input[i].transform[3] + y * input[i].transform[5]);

							if (geometry::within(point, envelop))
							{
								NEW = true;
								BOOST_FOREACH(BUFFER B, envelopes)
								{
									if (geometry::within(point, B))
									{
										NEW = false;
										break;
									}
								}
								
								if (NEW) 
								{
									Points << input[i].transform[0]+ x * input[i].transform[1] << "\t" << input[i].transform[3]+ y * input[i].transform[5] << endl;;
									data.push_back(input[i].values[x][y]);
								}
							}
							geometry::clear(point);
						}
					}
					
					if (data.size() > 1)
					{
						vec DATA(data.size(),fill::zeros);
						for(int i =0; i < (int) data.size(); i++)
							DATA(i) = data.at(i);

						if (arma::var(DATA) < variance)
						{
							variance = arma::var(DATA);
							radius -= 0.1 * input[i].transform[1];
							if (radius <= 0)
								accept = true;
						}
						else
							accept = true;
					}
					else 
						accept = true;
				} while (!accept);
				
				for(int i = 0; i < (int) data.size(); i++)
					txtF << data.at(i) << endl;
					
				fault_no++;
				envelopes.push_back(envelop);
			}
		}
	}
	std::copy(std::begin(old_GeoTransform), std::end(old_GeoTransform), std::begin(GeoTransform));
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

//create statisics from input file--------------------------------------
void STATS::CreateStats(ofstream& txtF, vector<line_type> faults)
{
	size_t index = 0;
	int FaultNumber = faults.size(), FaultGenerations;
	float mean_angle, med_angle, stddev_angle, var_angle, range_angle, strike, sinuosity;
	float mean_length, med_length, stddev_length, var_length, range_length, length;
	float mean_sin, med_sin, stddev_sin, var_sin, range_sin;
	point_type centre;
	double bin_size;
	//						index  wkt   strike  length  sinu   cdf_stri cdf_len cdf_sin nearf  distance to nearest fault
	std::vector<std::tuple<int, string, double, double, double, double, double, double, int, double, int, double>> faultInfo;
	vector<p_index> points;
	vector<pl_index> points2;
	double total_dist;
	
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

		string coord = "LINESTRING(";
		BOOST_FOREACH(point_type p, F)
		{
			geometry::append(p,multiL);
			coord += to_string(p.x()) + " ";
			coord += to_string(p.y()) + ", ";
		}
		coord += ")";

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

	vec sort_strike = sort(ANGLE);
	
	//histogram of stike----------------------------------------------------
	uvec strike_hist = hist(sort_strike, 10);
	bin_size = range_angle /10;
	vector<double> maxima;
	unsigned int max_nb = strike_hist.size()-1;

	for (unsigned int i = 0; i < strike_hist.size(); i++)
	{
		if (i == 0)
			if (strike_hist(i) > strike_hist(max_nb) && 
				strike_hist(i) > strike_hist(i+1))
					maxima.push_back(i*bin_size);
		
		
		if (i != 0 && i != max_nb)
			if (strike_hist(i) > strike_hist(i-1) &&
				strike_hist(i) > strike_hist(i+1))
					maxima.push_back(i*bin_size);
					
						
		if (i == max_nb)
			if (strike_hist(i) > strike_hist(i-1) && 
				strike_hist(i) > strike_hist(0))
					maxima.push_back(i*bin_size);
	}

	FaultGenerations = maxima.size();
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
	
	double scaL  = ScanLineDensity(faults);
	double D;
// 	DoBoxCount(faults);
	//BoxCount(faults, D);
	//LSFitter(len_dir);
	
	CompareStatsModels(fault_lengths);
	
	txtF.open ("Fault_Statistics.csv", ios::out | ios::app); 
	txtF << "Fault sets:" << "\t" << FaultGenerations << endl;
	txtF << "Average Scanline density: " << "\t" << scaL  << endl;
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

	//K-MEANS Clustering----------------------------------------------------
	txtF <<"k-means length and strike (" << FaultGenerations << ") sets: " << endl;
	result = KMCluster(false, FaultGenerations, len_dir);
	for ( typename FSTATS::const_iterator it = result.begin(); it != result.end(); it++)
		txtF << it->first << "\t" << it->second << endl;
	result.clear();
	txtF << endl; 
	
	txtF <<"k-means sinuosity and strike (" << FaultGenerations << ") sets" << endl;
	result = KMCluster(false, FaultGenerations, sin_dir);
	for ( typename FSTATS::const_iterator it = result.begin(); it != result.end(); it++)
		txtF << it->first << "\t" << it->second << endl;
	result.clear();
	txtF << endl; 
	
	txtF <<"k-means sinuosity and length (" << FaultGenerations << ") sets" << endl;;
	result = KMCluster(false, FaultGenerations, len_sin);
	for ( typename FSTATS::const_iterator it = result.begin(); it != result.end(); it++)
		txtF << it->first << "\t" << it->second << endl;
	result.clear();
	txtF.close();
 
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
				G[ve].offset = arma::mean(DATA);
			}
			Eg++;
		}
		txtF2.close();
	}
	cout << " analysed edges " << endl;
	//----------------------------------------------------------------------
	cout << "DEM Analysis: " << (clock() - startcputime) / (double)CLOCKS_PER_SEC << " seconds [CPU Clock] \n" << endl;
}

void STATS::LSFitter(FSTATS Length_Strike)
{
	boost::mt19937 rng; 
	/**
	* LENGTH:
	* lognormnal: McMahon, 1971; Baecher et al., 1977; Baecher and Lanney, 1978; Long et al., 1982; Dershowitz and Einstein, 1988
	* Exponential: Robertson, 1970; Call et al., 1976; Baecher et al., 1977; Long et al., 1982; Dershowitz and Einstein, 1988
	* Power Law: Heffer and Bevan, 1990; Bonnet et al., 2001; Park et al., 2001; De Dreuzy et al., 2002
	* 
	* STRIKE:
	* Normal: Long et al., 1982
	* Fisher: Fisher, 1953; Park et al., 2001
	*/
	vector<double> length, strike, log_n_length, exp_length, pow_length;

	for (typename FSTATS::const_iterator it = Length_Strike.begin(); it != Length_Strike.end(); it++)
	{
		length.push_back(it->first);
		strike.push_back(it->second);
	}

	vec L = conv_to<vec>::from(length);
	vec S = conv_to<vec>::from(strike);

	boost::lognormal_distribution<> nd(arma::mean(L), arma::stddev(L));
	boost::variate_generator<boost::mt19937&, 
		boost::lognormal_distribution<> > var_Ln(rng, nd);


	boost::exponential_distribution<> ex(1/(arma::mean(L)));
	boost::variate_generator<boost::mt19937&, 
		boost::exponential_distribution<> > var_ex(rng, ex);

	double exponent = -(1 + (1/arma::mean(L))) ;
	double dist_min = L.min();
	double dist_max = L.max();
	double min_pow = pow (dist_min, exponent + 1);
	double max_pow = pow (dist_max, exponent + 1);

	do {
		double ln_len = var_Ln();
		double ex_len = var_ex();
		double pl_len = (size_t) pow((max_pow - min_pow) * drand48 () + min_pow, 1 / (exponent + 1)); //why (size_t)?

		if (log_n_length.size() != length.size())
		{
			if (ln_len <= L.max() && ln_len >= L.min())
				log_n_length.push_back(ln_len);
		}
		if (exp_length.size() != length.size())
		{
			if (ex_len <= L.max() && ex_len >= L.min())
			{
				exp_length.push_back(ex_len);
			}
		}
		if (pow_length.size() != length.size())
		{
			if (pl_len <= L.max() && pl_len >= L.min())
			{
				exp_length.push_back(pl_len);
			}
		}

	}while (log_n_length.size() < length.size() && exp_length.size() < length.size());


	boost::math::lognormal_distribution<> myNormal(arma::mean(L), arma::stddev(L));
	for (typename FSTATS::const_iterator it = Length_Strike.begin(); it != Length_Strike.end(); it++)
		cout  << boost::math::pdf(myNormal, it->first) << endl;


	//ML estimator for lognormal distribution ------------------------------
	double mu_logn_1 = 0, mu_logn = 0, sigma_logn = 0;

	for (auto& X : length)
	{
		mu_logn_1  += log(X) / length.size();
		mu_logn    += mu_logn_1;
		sigma_logn += (log(X) - mu_logn_1) / length.size();
	}

	cout <<mu_logn << "("<< arma::mean(L) <<") " << sigma_logn << "(" << sqrt(arma::stddev(L)) << ")" << endl; 

	boost::lognormal_distribution<> nd2(mu_logn, sqrt(sigma_logn));
	boost::variate_generator<boost::mt19937&, 
		boost::lognormal_distribution<> > var_Ln2(rng, nd2);

		
	boost::math::lognormal_distribution<> logn_estimated(mu_logn, sqrt(sigma_logn));


	for (auto& l : length)
	cout <<var_Ln2() << endl;
	//double ex_len = var_ex();
	//	cout  << boost::math::pdf(logn_estimated, l) << endl;


}


