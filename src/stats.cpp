#include "../include/stats.h"
#include "../include/geometrie.h"
#include "../include/GeoRef.h"
#include "../include/graph.h"
#include "../include/util.h"

const std::string stats_subdir = "/statisics/";

namespace FracG
{


	// create and open filestream in folder "statistics"
	std::ofstream CreateFileStream(std::string subdir_name, std::string name)
	{
		std::string stats_dirname = subdir_name + stats_subdir;
		const char* stats_dir = stats_dirname.c_str();
			if(!opendir(stats_dir))

		mkdir(stats_dir, 0777);
		std::string csv_file = stats_dir + name;
		std::ofstream txtF; 
		txtF.open (csv_file, std::ios::out); 
		return(txtF);
	}

	//struct to hold the quad tree data
	struct QuadTreeNode
	{
		box bbox; //bounding box for the square
		std::vector<box> child_bboxes; //bboxes for this node's children
		std::vector<std::unique_ptr<QuadTreeNode>> children;
		int generation;
		//how do I have this be null until i set them? if i use a std::vector, i cant guarantee that the nodes will line up with the bounding boxes? c++ is wierd
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
			std::cout << "ERROR - node " << node << " has children length " << node.children.size() << " at layer " << layers << std::endl;
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
			counts.push_back(std::make_tuple(side_length, 0));
		}
		std::get<1>(counts[index]) += 1;//count this box
		//now check its children
		for (auto it = node.children.begin(); it != node.children.end(); it++)
		{
			if (*it != nullptr) CountQuadTreeBoxes(**it, counts, index+1);
		}
	}

	void BoxCountQuadTree(const std::vector<line_type> &faults, std::vector<std::tuple<double, long long int>> &counts, const double start_scale, const int max_depth, point_type &offset)
	{
		box bbox(faults.front().front(), faults.front().front()); //bounding box for the entire fault set. start with an arbitrary point that is known to be inside the bounding box of the faults
		std::vector<box> bboxes;
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

	// 	std::cout << "FINAL bbox = x(" << bbox.min_corner().x()<<", "<<bbox.max_corner().x()<<"), y("<<bbox.min_corner().y()<<", "<<bbox.max_corner().y()<<")"<<std::endl;
	// 	std::cout << "x = " << xStart << ", " << xEnd << ", y = " << yStart << ", " << yEnd << std::endl;

		//setup the top-most layer of the Quad Tree
		std::vector<std::vector<box>> node_bboxes; //the bounding boxes for the top-most layer of QuadTreeNode s
		std::vector<std::vector<std::unique_ptr<QuadTreeNode>>> nodes; //the (nodes that are) the top-most layer of QuadTreeNode s
		node_bboxes.reserve(nX);
		nodes.reserve(nX);
		for (int ix = 0; ix < nX; ix++)
		{
			std::vector<box> bbox_row;
			bbox_row.reserve(nY);
			std::vector<std::unique_ptr<QuadTreeNode>> node_row;
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
				std::cout << "WARNING: box counting minimum x boundary returned " << min_x << std::endl;
				min_x = 0;
			}
			if (max_x >= nX)
			{
				std::cout << "WARNING: box counting maximum x boundary returned " << max_x << " when using " << nX << " boxes, from position (" << fault_bbox.min_corner().x() << ", " << fault_bbox.max_corner().x() << "), (" << fault_bbox.min_corner().y() << ", " << fault_bbox.max_corner().y() << ")" << std::endl;
				max_x = nX-1;
			}
			if (min_y < 0)
			{
				std::cout << "WARNING: box counting minimum y boundary returned " << min_y << std::endl;
				min_y = 0;
			}
			if (max_y >= nY)
			{
				std::cout << "WARNING: box counting maximum y boundary returned " << max_y << " when using " << nY << "boxes" << std::endl;
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

	void DoBoxCount(VECTOR lines)
	{
		std::cout << "Box Counting (QuadTree Method)" << std::endl;
		std::ofstream txtF = CreateFileStream(lines.out_path / stats_subdir / ("box_count.tsv"));

		std::clock_t startcputime = std::clock();
		auto t_start = std::chrono::high_resolution_clock::now();

		//use the length of the longest fault as the starting box size
		double longest = 0;
		double shortest = std::numeric_limits<double>::infinity();
		for (auto it = lines.data.begin(); it != lines.data.end(); it++)
		{
			const double length = geometry::length(*it);
			longest = std::max(longest, length);
			shortest = std::min(shortest, length);
		}

		const int max_depth = 6;
		std::vector<std::tuple<double, long long int>> counts;
		point_type offset(0,0);
		BoxCountQuadTree(lines.data, counts, longest, max_depth, offset);

		auto t_end = std::chrono::high_resolution_clock::now();
		std::cout << " CPU  time: " << (clock() - startcputime) / (double)CLOCKS_PER_SEC << " seconds\n"
			 << " Wall time: " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << " ms" << std::endl;

	//Least-Square-Fit------------------------------------------------------
		int NB = 0;
		arma::vec X(counts.size(),  arma::fill::zeros);
		arma::vec Y(counts.size(),  arma::fill::zeros);
		arma::vec XX(counts.size(), arma::fill::zeros);
		arma::vec XY(counts.size(), arma::fill::zeros);

		txtF <<" frequency \t size" << std::endl;
		for (auto it = counts.begin(); it != counts.end(); it++)
		{
			txtF << std::get<1>(*it) << "\t" << std::get<0>(*it) << std::endl;
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

		txtF <<  "No \t x(S) \t log(y(observed)) \t log(y(fitted))" << std::endl;

		for (int ii=0; ii<N; ii++)
			txtF <<"" << ii+1 << "\t"<< X[ii] << "\t" <<Y[ii] << "\t" << y_fit[ii] << std::endl;    

		txtF << " LS-Fit: "<< "\t" << a  << "x + " << b << std::endl; 
		txtF << " R^2: "   << "\t" << R2 << std::endl; 
		std::cout << "Linear fit: "<<a<<"x + " << b << " (R^2 = " << R2 << ")"  << std::endl;

		delete[] y_fit;
		txtF.close();
		std::cout << " done \n" << std::endl;
	}

	//caluclates the Kolmogorov-Smirnov test for a given std::vector of values, a function that gives the cdf for the distribution being tested, the parameters for the cdf model, and the minimum x value to iterate from
	//the data in values *must* be std::sorted (so that tha CDF for the data values is calculated in the correct order)
	template<typename T>
	double KSTest(std::vector<double> &values, const ModelHolder<T> &model, const T &params, const std::vector<double>::iterator &xmin)
	{
	// 	//for debug purposes
	// 	std::ostringstream fname;
	// 	fname << "kstest_" << (xmin - values.begin()) << ".txt";
	// 	std::ofstream outF(fname.str());
	// 	outF << "#x data_cdf_before data_cdf_after distribution_cdf" << std::endl;

		const int n_tail = values.end() - xmin;
	// 	for (auto it = values.begin(); it != values.end(); it++) if (*it >= xmin) n_tail++; //not needed, because the data is already std::sorted
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
	// 		outF << *it << " " << data_cdf_before << " " << data_cdf_after << " " << dist_cdf_value << std::endl;
		}

		return max_diff;
	}

	//Calculate the best paramaters for the distribution, for the given dataset. Iterates over each element of the values to check if it is the best option for xmin
	template<typename T>
	void CalculateBestDistParams(std::vector<double> &values, const ModelHolder<T> &model, T &params, std::vector<double>::iterator &xmin, double &best_ks)
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
		std::vector<double>::iterator xmax = values.end();
		while ((xmax > values.begin()) && (xmax == values.end() || *xmax >= values.back())) --xmax; //don't use xmin here, it isn't inisialised. this is the function that is supposed to initialise it.
		#pragma omp parallel for
		for (auto it = values.begin(); it <= xmax; it++)
		{
			//first check to see if we have any values that are strictly larger than our proposed minimum
			if (*it >= values.back()) continue; //they're std::sorted, so if this one doesn't have anything higher than it, none of the others will either
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

	//get the std::vector of values that are below xmin
	void GetLowerValues(const std::vector<double> &values, const double xmin, std::vector<double> lower, int &n_tail)
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
	PowerLawParams PowerLawCalculateParams(const std::vector<double> &data, const std::vector<double>::iterator &xmin)
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
	std::vector<double> StandardGenerateMultiValue(const std::vector<double> &r, const double xmin, const T &params, std::function<double(const double, const double, const T &params)> convert_func)
	{
		std::vector<double> samples;
		samples.reserve(r.size());
		for (auto it = r.begin(); it != r.end(); it++) samples.push_back(convert_func(*it, xmin, params));
		return std::move(samples);
	}

	//calculate exponential distribution parameters for pdf p(x,lambda) = (lambda) * exp(- lambda * x)
	//cdf(x, lambda) = 1 - exp[- lambda * x]
	//default xmin = 0
	ExponentialParams ExponentialCalculateParams(const std::vector<double> &data, const std::vector<double>::iterator &xmin)
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
	LogNormParams LogNormCalculateParams(const std::vector<double> &data, const std::vector<double>::iterator &xmin)
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

	std::vector<double> LogNormGenerateMultiValue(const std::vector<double> &r, const double xmin, const LogNormParams &params)
	{
		//this doesn't use mu, but then again the power law one doesn't use c either.
		const int N = (int)r.size() / 2; //making this many std::pairs of values
		std::vector<double> samples;
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

	WeibullParams WeibullCalculateParams(const std::vector<double> &data, const std::vector<double>::iterator &xmin)
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
				return std::make_pair(func_val, diff_val); //I think these calculations are correct, but I'm not entirely sure
			},
			//intial guess, min, max, number of digits precision, max iterations
			2.0, 0.001, 20.0, 5, max_iters);
	// 	std::cout << "used " << max_iters << " iterations" << std::endl;
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
	std::vector<double> GetSampleSingles(const ModelHolder<T> &model, const T &params, const int n_samples, const std::vector<double>::iterator &xmin, random::mt19937 &gen, const std::vector<double> &values)
	{
		std::vector<double> samples;
		samples.reserve(n_samples);
		const int max_low = std::max(0, (int) (xmin - values.begin() - 1)); //index of the final value in the data that is below xmin. Inclusive value.
		const double p_lower = (xmin - values.begin())/(double)values.size(); //probability with which to use the below-xmin portion of the samples 
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
	std::vector<double> GetSampleMulti(const ModelHolder<T> &model, const T &params, const int n_samples, const std::vector<double>::iterator &xmin, random::mt19937 &gen, const std::vector<double> &values)
	{
		const int per_call = model.samples_per_call;
		const int n_generations = std::lrint(std::ceil(n_samples/per_call));
		std::vector<double> samples;
		samples.reserve(n_generations*per_call);
		const int max_low = std::max(0, (int) (xmin - values.begin() - 1)); //index of the final value in the data that is below xmin. Inclusive value.
		const double p_lower = (xmin - values.begin())/values.size(); //probability with which to use the below-xmin portion of the samples 
		std::uniform_real_distribution<double> uniform01(0.0, 1.0);
		std::uniform_int_distribution<int> lower_selector(0, max_low);
		std::vector<double> random_values;
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
				std::vector<double> model_values = model.convert_multi_func(random_values, *xmin, params);
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
			std::ofstream &fname;
			PrintCDFSingle(const std::vector<double> &v, std::ofstream &fn) : values(v), fname(fn) {};

			template<typename T>
			void operator()(DistStats<T> &stats)
			{
				const int NV = values.size(); //total number of data values
				const int N  = NV - stats.xmin_ind; //number of data values that are above xmin
				const double ND = N; //save having to convert this multiple times
				//std::string model_fname = fname+"_"+stats.name+".txt";
				//std::replace(model_fname.begin(), model_fname.end(), ' ', '_');
				//std::transform(model_fname.begin(), model_fname.end(), model_fname.begin(), ::tolower);

				fname << "CUMMULATIVE DISTRIBUTION FUNCTION" << std::endl;
				fname << "#" << stats.name << ": " << stats.params << ", xmin = " << stats.xmin << ", xmin_ind = " << stats.xmin_ind << " of " << NV << std::endl;
				fname << "#x_value" << "\t"<< "data_cdf_low" <<"\t"<< "data_cdf_high" <<"\t"<< "model_cdf" << std::endl;

				for (int i = stats.xmin_ind; i < NV; i++)
				{
					const double data_low = (i     - stats.xmin_ind) / ND;
					const double data_hi =  (i + 1 - stats.xmin_ind) / ND;
					const double model_cdf = stats.model.cdf_func(values.at(i), stats.xmin, stats.params);
					fname << values[i] << "\t" << data_low << "\t" << data_hi << "\t" << model_cdf << std::endl;
				}
				fname << std::endl;
			}
	};

	//visitor struct to print out the pdf of a single model
	//we have to print out the different models separately, because different xmin's mean different data PDF's
	struct PrintPDFSingle : public boost::static_visitor<void>//(DistStats<T> &stats)
	{
		public:
			const std::vector<double> &values;
			std::ofstream &fname;
			PrintPDFSingle(const std::vector<double> &v, std::ofstream &fn) : values(v), fname(fn) {};

			template<typename T>
			void operator()(DistStats<T> &stats)
			{
				fname << "PROBABILITY DENSITY FUNCTION" << std::endl;
				const int NV = values.size(); //total number of data values
				const int N  = NV - stats.xmin_ind; //number of data values that are above xmin
				const double ND = N; //save having to convert this multiple times

				fname << "#" << stats.name << ": " << stats.params << ", xmin = " << stats.xmin << ", xmin_ind = " << stats.xmin_ind << " of " << NV << std::endl;
				fname << "#x_value" << "\t" << "model_pdf" << std::endl; 

				for (int i = stats.xmin_ind; i < NV; i++)
				{
					const double data_low = (i     - stats.xmin_ind) / ND;
					const double data_hi =  (i + 1 - stats.xmin_ind) / ND;
					const double model_pdf = stats.model.pdf_func(values.at(i), stats.xmin, stats.params);
					fname << values[i] << "\t" << model_pdf << std::endl;
				}
				fname << std::endl;
			}
	};

	//print out all the model CDF's, as compared to their relative data CDF's
	void PrintBest2File(const StatsModelData dists, const std::vector<double> values, std::ofstream &fname)
	{
		PrintCDFSingle print_visitor_c(values, fname);
		PrintPDFSingle print_visitor_p(values, fname);

		StatsModelData printModel;

		if (dists.best_matches.size() == 1)
		{
			const int best = dists.best_match;
			if (best >= 0)
			{
				auto s = dists.models[best];
				boost::apply_visitor(print_visitor_c, s);
				boost::apply_visitor(print_visitor_p, s);

			}

		} else if (dists.best_matches.size() > 1) 
		{
			for (auto it = dists.best_matches.begin(); it < dists.best_matches.end(); it++)
			{
				auto s = dists.models[*it];
				boost::apply_visitor(print_visitor_c, s);
				boost::apply_visitor(print_visitor_p, s);
			}
		}
	}

	template<typename T1, typename T2>
	std::tuple<double, double> LogLikelihoodRatio(const ModelHolder<T1> &model1, const T1 &params1, const ModelHolder<T2> &model2, const T2 params2, const std::vector<double> &values, const std::vector<double>::const_iterator &xmin1, const std::vector<double>::const_iterator &xmin2)
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
	// 	std::cout << "comparison: R = " << R << ", sigma = " << sigma << ", sigma_sum = " << sigma_sum <<", N = " << N << ", erfc input = " << std::abs(R)/(std::sqrt(2*N)*sigma) << std::endl;
		const double p = (sigma_sum <= 0.0) ? 1 : erfc(std::abs(R)/(std::sqrt(2*N)*sigma)); //lower <-> more confidence that R value is meaningful
		return std::make_tuple(R, p);//std::move()
	}

	//generate a number of samples from a distribution, including the values that are below xmin
	//split into single and multi calls, so as to avoid extra overhead in the (more common) single value case
	template<typename T>
	std::vector<double> GetSampleValues(const ModelHolder<T> &model, const T &params, const int n_samples, const std::vector<double>::iterator &xmin, random::mt19937 &gen, const std::vector<double> &values)
	{
		if (model.samples_per_call == 1) return std::move(GetSampleSingles<T>(model, params, n_samples, xmin, gen, values));
		return std::move(GetSampleMulti<T>(model, params, n_samples, xmin, gen, values));
	}

	template<typename T>
	double GetPValue(const ModelHolder<T> &model, const T &params, const std::vector<double>::iterator &xmin, const std::vector<double> &values, const double data_ks, std::vector<random::mt19937> &gen)
	{
		const double err = 0.01; //accuracy to two significant figures
		const int ntrials = std::lrint(std::ceil(0.25/(err*err))); //2500
		int npass = 0;
		#pragma omp parallel for reduction(+:npass)
		for (int i = 0; i < ntrials; i++)
		{
			random::mt19937 &local_generator = gen[omp_get_thread_num()]; //each thread uses its own RNG
			std::vector<double> samples = GetSampleValues<T>(model, params, (int)values.size(), xmin, local_generator, values);
			std::vector<double>::iterator sample_xmin;
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
		//cout << "Model " << stats.name << " with parameters " << stats.params << " has ks value " << stats.ks << ", pvalue = " << stats.pvalue << ", with xmin = " << stats.xmin << ", at index " << stats.xmin_ind << " of " << values.size() << std::endl;
			}
	};

	//compare the models against each other, in order to determine which one(s) are the best fit
	void DecideBestModel(StatsModelData &model_data, const std::vector<double> &values, const double ACCEPT_THRESHOLD, const double DIFFERENCE_THRESHOLD)
	{
		model_data.best_match = -1;
		model_data.best_matches.clear();
		bool have_model = false;
		std::vector<DistStatsType> &models = model_data.models; //convenience name
		const int N = model_data.models.size(); //Number of models
		//initialise comparisons as a NxN 2D std::vector, initialised with R=0 (no difference), p = 1 (no certainty that there is a meaningful difference even if R is large)//= std::vector<
		model_data.comparisons = std::vector<std::vector<std::tuple<double,double>>>(N, std::vector<std::tuple<double, double>>(N));

		//first, check to see if any models are acceptable.
		for (auto it = models.begin(); it != models.end(); it++){
			double pvalue = boost::apply_visitor([](auto& x)->double{return x.pvalue;}, *it);
			if (pvalue >= ACCEPT_THRESHOLD)
			{
				have_model = true;
				model_data.best_match = it - models.begin();
			//	std::cout << "Model " << boost::apply_visitor([](auto &x){return x.name;}, *it) << " has not been rejected" << std::endl;
			}
		}
	// 	std::cout << std::endl;

		//if none are acceptable, stop now. there is no point in comparing them
		if (!have_model)
		{
			std::cout << "All models are below the acceptance threshold" << std::endl;
			return;
		}

		//do a std::pairwise comparison between all the acceptable models
		for (auto it1 = models.begin(); it1 != models.end(); it1++)
		{
			if (boost::apply_visitor([&](auto &x){return x.pvalue < ACCEPT_THRESHOLD;}, *it1)) continue; //if we don't accept this model, don't compare it to other models
			int ind1 = it1 - models.begin();
			std::string name1 = boost::apply_visitor([](auto &x){return x.name;}, *it1);
			//need const_iterator because values is a const std::vector now. (also need const_iterator's in log_likelihood_ratio())
			const std::vector<double>::const_iterator xmin1 = std::next(values.begin(), boost::apply_visitor([](auto &x){return x.xmin_ind;}, *it1)); //we can save this one to save a bit of time
			for (auto it2 = it1; it2 != models.end(); it2++)
			{
				if (it1 == it2) continue;
				if (boost::apply_visitor([&](auto &x){return x.pvalue < ACCEPT_THRESHOLD;}, *it2)) continue;
				int ind2 = it2 - models.begin();
				std::string name2 = boost::apply_visitor([](auto &x){return x.name;}, *it2);
				std::tuple<double, double> comparison_result = boost::apply_visitor(
					[&](auto &x1, auto &x2)->std::tuple<double, double>{return LogLikelihoodRatio(x1.model, x1.params, x2.model, x2.params, values, xmin1, std::next(values.begin(), x2.xmin_ind));}
					, *it1, *it2);
				double R, p;
				std::tie(R, p) = comparison_result;
				model_data.comparisons[ind1][ind2] = std::move(comparison_result);
				if (p > DIFFERENCE_THRESHOLD)
				{
					//cout << "Models " << name1 << " and " << name2 << " cannot be distinguished, with R = " << R << ", and p = " << p << std::endl;
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
				//cout << "Model " << best_name << " is a better fit that model " << worse_name << ", with R = " << R << ", p = " << p << std::endl;
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
			std::cout << "There is a loop in the relation of the best models" << std::endl;
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
			std::cout << "The loop is: ";
			for (auto it = model_data.best_matches.begin(); it < model_data.best_matches.end(); it++)
				std::cout << ((it == model_data.best_matches.begin()) ? "" : ", ") << boost::apply_visitor([](auto &x){return x.name;}, model_data.models[*it]);
			std::cout << std::endl;
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

	StatsModelData GetLengthDist(VECTOR lines)
	{
		StatsModelData results;
		random::random_device rd{};
		std::vector <double> values;

		std::ofstream txtF = FracG::CreateFileStream(lines.out_path / stats_subdir / ("length_distributions.tsv"));

		BOOST_FOREACH(line_type l , lines.data)
			values.push_back(geometry::length(l));

		std::sort(values.begin(), values.end());
		std::cout << "Determing length distribution for values in the range " << values.front() << ", " << values.back() << std::endl;

		const int NTHREADS = std::max(1, omp_get_max_threads());

		std::cout << "Using " << NTHREADS << " threads" << std::endl;
		std::vector<random::mt19937> rngs(NTHREADS);
		for (auto it = rngs.begin(); it < rngs.end(); it++) *it = std::move(random::mt19937{rd}); //make one random number generator per omp thread
		ModelHolder<PowerLawParams> power_law_model = {PowerLawCalculateParams, PowerLawCDF, PowerLawPDF, PowerLawGenerateValue,
			[](const std::vector<double>&r, const double xmin, const PowerLawParams& params) -> std::vector<double> {return StandardGenerateMultiValue<PowerLawParams>(r, (const double)xmin, params, PowerLawGenerateValue);}, true, 1};
		ModelHolder<PowerLawParams> power_law_all_model = {PowerLawCalculateParams, PowerLawCDF, PowerLawPDF, PowerLawGenerateValue,
			[](const std::vector<double>&r, const double xmin, const PowerLawParams& params) -> std::vector<double> {return StandardGenerateMultiValue<PowerLawParams>(r, (const double)xmin, params, PowerLawGenerateValue);}, false, 1};

		ModelHolder<ExponentialParams> exponential_model = {ExponentialCalculateParams, ExponentialCDF, ExponentialPDF, ExponentialGenerateValue,
			[](const std::vector<double>&r, const double xmin, const ExponentialParams& params) -> std::vector<double> {return StandardGenerateMultiValue<ExponentialParams>(r, (const double)xmin, params, ExponentialGenerateValue);}, true, 1};
		ModelHolder<ExponentialParams> exponential_all_model = {ExponentialCalculateParams, ExponentialCDF, ExponentialPDF, ExponentialGenerateValue,
			[](const std::vector<double>&r, const double xmin, const ExponentialParams& params) -> std::vector<double> {return StandardGenerateMultiValue<ExponentialParams>(r, (const double)xmin, params, ExponentialGenerateValue);}, false, 1};

		ModelHolder<LogNormParams> lognorm_model = {LogNormCalculateParams, LogNormCDF, LogNormPDF, LogNormGenerateValue, LogNormGenerateMultiValue, true, 2};
		ModelHolder<LogNormParams> lognorm_all_model = {LogNormCalculateParams, LogNormCDF, LogNormPDF, LogNormGenerateValue, LogNormGenerateMultiValue, false, 2};

		/*
		ModelHolder<WeibullParams> weibull_model = {WeibullCalculateParams, WeibullCDF, WeibullPDF, WeibullGenerateValue,
			[](const std::vector<double>&r, const double xmin, const WeibullParams& params) -> std::vector<double> {return StandardGenerateMultiValue<WeibullParams>(r, (const double)xmin, params, WeibullGenerateValue);}, true, 1};
		ModelHolder<WeibullParams> weibull_all_model = {WeibullCalculateParams, WeibullCDF, WeibullPDF, WeibullGenerateValue,
			[](const std::vector<double>&r, const double xmin, const WeibullParams& params) -> std::vector<double> {return StandardGenerateMultiValue<WeibullParams>(r, (const double)xmin, params, WeibullGenerateValue);}, false, 1};
		*/

		std::vector<DistStatsType> &models = results.models; //convenience name

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

		const double ACCEPT_THRESHOLD = 0.1; //models are accepted if their pvalue is *ABOVE* this threshold

		const double DIFFERENCE_THRESHOLD = 0.1; //comparisons are considered meaningful if the pvalue is *BELOW* this threshold
		DecideBestModel(results, values, ACCEPT_THRESHOLD, DIFFERENCE_THRESHOLD); //compare the models against each other



		if (results.best_matches.size() == 1)
		{
			const int best = results.best_match;
			if (best >= 0)
			{
			std::cout << "The best model is " << "\t" << boost::apply_visitor([](auto& x)->std::string{return x.name;}, results.models[best]) << std::endl << std::endl;
			txtF << "The best model is " << "\t" << boost::apply_visitor([](auto& x)->std::string{return x.name;}, results.models[best]) << std::endl;
			}
		} else if (results.best_matches.size() > 1) 
		{
			std::cout << "The " << results.best_matches.size() << " best matches are: ";
			txtF << "The " << results.best_matches.size() << " best matches are: " << "\t";
			for (auto it = results.best_matches.begin(); it < results.best_matches.end(); it++)
			{
				std::cout << (it == results.best_matches.begin() ? "" : ", ") << boost::apply_visitor([](auto& x)->std::string{return x.name;}, results.models[*it]);
				txtF << (it == results.best_matches.begin() ? "" : ", ") << boost::apply_visitor([](auto& x)->std::string{return x.name;}, results.models[*it]);
			}
			std::cout << std::endl;
			txtF << std::endl;		
		}

		PrintBest2File(results, values, txtF); //print the cdf's for testing purposes
		txtF.close();
		std::cout << std::endl;
		return results;
	}

	double PointExtractor(point_type P, double radius, double Transform[8], double** raster)
	{
		box AOI;
		BUFFER envelop;
		point_type point;
		int maxX, minX, maxY, minY;
		double M = 0;
		std::vector<double> D;

		envelop = DefinePointBuffer(P, radius);
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
			arma::vec DATA(D.size(), arma::fill::zeros);
			for(int i = 0; i < (int) D.size(); i++)
				 DATA(i) = D.at(i);
			M = arma::mean(DATA);
		}
		return(M);
	}

	void RasterStatistics(VECTOR lines, double dist, std::string raster_filename)
	{
		std::cout << "Generating raster statisics " << std::endl;
		const char * name = raster_filename.c_str();
		enum {ERROR, Byte, UInt16, Int16, UInt32, Int32, Float32, Float64 };
		static std::map<std::string, int> d_type;
		if ( d_type.empty() )
		{
			d_type["Byte"] = Byte;
			d_type["UInt16"] = UInt16;
			d_type["Int16"] = Int16;
			d_type["UInt32"] = UInt32;
			d_type["Int32"] = Int32;
			d_type["Float32"] = Float32;
			d_type["Float64"] = Float64;
		}

		GDALDataset  *poDataset;
		GDALAllRegister();
		std::cout << "RasterStatistics is opening \"" << name << "\"" << std::endl;
		poDataset = (GDALDataset *) GDALOpen( name, GA_ReadOnly );
		GDALRasterBand  *poBand;
		poBand = poDataset->GetRasterBand( 1 );
		std::string datatype (GDALGetDataTypeName(poBand->GetRasterDataType()));
		GDALClose( poDataset );

		std::cout << "Raster datatype is " ;
		int type = d_type[datatype];
		 switch (type) 
		 {
			case Byte:
			{
				std::cout	<< "byte (will not perform analysis)" << std::endl;
			} break;

			case UInt16:
			{
				std::cout << "uint16" << std::endl;
				RASTER<int> R;
				R.name = raster_filename;
				AnalyseRaster<int>(lines, dist, R);
			} break;

			case Int16:
			{
				std::cout << "int16" << std::endl;
				RASTER<int> R;
				R.name = raster_filename;
				AnalyseRaster<int>(lines, dist, R);
			} break;

			case UInt32:
			{
				std::cout << "uint32" << std::endl;
				RASTER<int> R;
				R.name = raster_filename;
				AnalyseRaster<int>(lines, dist, R);

			} break;

			case Int32:
			{
				std::cout << "int32" << std::endl;
				RASTER<int> R;
				R.name = raster_filename;
				AnalyseRaster<int>(lines, dist, R);
			} break;

			case Float32:
			{
				std::cout << "Float32" << std::endl;
				RASTER<float> R;
				R.name = raster_filename;
				AnalyseRaster<float>(lines, dist, R);
			} break;

			case Float64:
			{
				std::cout << "Float64" << std::endl;
				RASTER<double> R;
				R.name = raster_filename;
				AnalyseRaster<double>(lines, dist, R);
			} break;

			default: 
			{
				std::cout << " unknown" << std::endl;
				std::cout << "ERROR: Could not determine dataype of raster" << std::endl;
			}
			break;
		 }
		 std::cout << " done \n" << std::endl;
	}

	void MA_filter(std::vector< std::pair<double,double>> &KDE, int window)
	{
		if (KDE.size() > 5 * (unsigned int) window)
		{
			unsigned int wn = (int) (window - 1) / 2;
			std::vector<float>KDE_filtered(KDE.size());

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
			for (unsigned int i =0; i < KDE.size();i++)
			{
				KDE[i].second = KDE_filtered[i];
			}
		}
	}


	//evaluate the kernel density esimation, with an gaussian kernel, at an arbitrary point
	//angles is the array of angles in degrees, smooth is a smoothing factor (>0), and x is the point to be evaluated at
	double EvaluateAngleKde(arma::vec &angles, const double smooth, double x)
	{
		x = std::fmod(x, MAX_ANGLE);
		double sum = 0;
		const double s2 = smooth*smooth;
		for (auto it = angles.begin(); it < angles.end(); it++)
		{
			//because the kde wraps around, exavluate it on either side of the gaussian
			const double x0 = *it;
			double xd1 = x - x0;
			double xd2 = MAX_ANGLE - std::abs(xd1);
			sum += std::exp(-0.5 * xd1 * xd1 / s2);
			sum += std::exp(-0.5 * xd2 * xd2 / s2);
		}
		return sum/(sqrt(2*M_PI) * angles.size() * smooth);
	}

	//return an evenly-sampled array of the kde that represents the given values and smoothing factor
	//and scale them to sum to 1, to be a probability distribution function
	arma::vec AngleKdeFillArray(int nsamples, arma::vec &angles, const double smooth)
	{
		const double dn = (double) nsamples; //convenience name, and to avoid forgetting to convert to flaoting point later
		arma::vec sampled(nsamples);
		for (int i = 0; i < nsamples; i++){
			sampled[i] = EvaluateAngleKde(angles, smooth, MAX_ANGLE * i / dn);// * MAX_ANGLE / dn //actually, don't scale it here
		}
		return sampled;
	}


	//perform a moving average of the KDE, that wraps around
	void MovingAverageFilterWraparound(std::vector< std::pair<double,double>> &KDE, int window_size)
	{
		window_size |= 1; //make the window_size odd-sized
		const unsigned int hw = window_size / 2;
		std::vector<double> saved(window_size); //saved values of the moving sum, so we know which value to remove for the next index
		std::vector<double> copy(window_size); //a copy to hold the values that get overwritten in the first pass
		//initialise values
		double sum = 0;
		for (int i = 0; i < window_size; i++)
		{
			sum += KDE[i].second;
			saved[i] = KDE[i].second;
			copy[i] = KDE[i].second;
		}
		int saved_index = 0;
		const unsigned int last_index = KDE.size() - hw - 1;
		//do the regular part of the moving average
		for (unsigned int i = hw; i < last_index; i++){
			KDE[i].second = sum / window_size;
			const int next_index = i + hw + 1;
			sum += KDE[next_index].second - saved[saved_index];
			saved[saved_index] = KDE[next_index].second;
			saved_index = (saved_index + 1) % window_size;
		}
		//now do the wrap-around part
		for (unsigned int i = last_index; (i >= last_index) || (i < hw); i = (i + 1) % KDE.size())
		{
			KDE[i].second = sum / window_size;
			const int next_index = (i + hw + 1) % KDE.size(); //the next value will wrap around
			const double next_value = next_index < window_size ? copy[next_index] : KDE[next_index].second; //either get the enxt index as normal, or get it from the saved copy if it coems from the part that is overwritten
			sum += next_value - saved[saved_index];
			saved[saved_index] = next_value;
			saved_index = (saved_index + 1) % window_size;
		}
	}

	typedef std::pair<int, bool> crossing_type; //crossing position, whether the crossing is rising or falling
	typedef std::pair<int, int> crossing_location_type; //the falling and rising edges of a gaussian candidate location
	//find zero-crossings, from a fourier-domain representation of the data
	std::vector<crossing_type> FindZerocrossings(arma::cx_vec &fd, arma::vec &freqs)
	{
		arma::cx_vec d2(fd.size());
		for (unsigned int i = 0; i < fd.size(); i++) d2[i] = -freqs[i]*freqs[i]*fd[i];
		arma::cx_vec angle = arma::ifft(d2); //arma doesn't have a real-valued fft/ifft

		std::vector<crossing_type> crossings;
		for (unsigned int i = 0; i < fd.size(); i++)
		{
			bool sign_this = !std::signbit(std::real(angle[i])); //use not, so that sign being true <-> the value is positive
			bool sign_prev = !std::signbit(std::real(angle[(i - 1 + fd.size()) % fd.size()]));
			if (sign_this != sign_prev)
			{
				const bool rising = sign_this;
				crossings.push_back(std::make_pair(i, rising));
	// 			std::cout << " found crossing at " << i << ": " << rising << std::endl;
			}
		}
	// 	std::cout << std::endl;
		return crossings;
	}


	//a simpler matching algorithm, that finds gaussian locations by matching leading edges to the next falling edge
	std::vector<crossing_location_type> SimpleLocationDetection(std::vector<crossing_type> &crossings)
	{
		std::vector<bool> used(crossings.size(), false);
		const int nGauss = crossings.size() / 2;
		std::vector<crossing_location_type> found;

		while ((int)found.size() < nGauss)
		{
			for (int i = 0; i < (int) crossings.size(); i++){
				//find a leading edge (that we haven't used yet)
				if (used[i] || crossings[i].second) continue;
				for (int j = (i + 1) % crossings.size(); j != i; j = (j + 1) % crossings.size())
				{
					//and the next falling edge (that we haven't used yet)
					if (used[j] || !crossings[j].second) continue;
					found.push_back(std::make_pair(crossings[i].first, crossings[j].first));
					used[i] = true;
					used[j] = true;
					break;
				}
			}
		}
		return found;
	}


	//data structured used for nonlinear function fitting with the Gnu Scientific Library
	struct gauss_data{
		int nData;
		int nGauss;
		arma::vec &angle;
		arma::vec &pdf;
		bool use_horizontal;
		gauss_data(int nd, int ng, arma::vec &a, arma::vec &p, bool use_horiz) : nData(nd), nGauss(ng), angle(a), pdf(p), use_horizontal(use_horiz) {};
	};

	//evaluate the gaussian functions for a particular parameter set
	int GaussSumFunc(const gsl_vector *params, void *data_ptr, gsl_vector *diff)
	{
		struct gauss_data *data = (struct gauss_data *) data_ptr;

	// 	for (int g = 0; g < data->nGauss; g++)
	// 		std::cout << "g" << g << ": amp = " << gsl_vector_get(params, PPG*g    ) << ", sigma = " << gsl_vector_get(params, PPG*g + 1) << ", pos = " << gsl_vector_get(params, PPG*g + 2) << std::endl;
	// 	std::cout << std::endl;


		for (int i = 0; i < data->nData; i++)
		{
			double sum = 0;
			for (int g = 0; g < data->nGauss; g++){
				const double amp = std::abs(gsl_vector_get(params, PPG*g    )); //amplitude
				const double sig = gsl_vector_get(params, PPG*g + 1); //sigma
					  double pos = gsl_vector_get(params, PPG*g + 2); //position
				pos = fmod(fmod(pos, MAX_ANGLE) + MAX_ANGLE, MAX_ANGLE);
				const double d1 = data->angle[i] - pos;
				//this does: d2 = d1 + 180, or d1 - 180, whichever results in a value that is in the range (-180, 180) and has the opposite sign of d1
				const double d2 = d1 - std::copysign(180, d1); //make sure that if d1 is positive, then d2 is negative, and vice versa
				const double z1 = d1/sig;
				const double z2 = d2/sig;
				sum += amp * exp(-0.5 * z1 * z1);
				sum += amp * exp(-0.5 * z2 * z2);
			}
			if (data->use_horizontal) sum += (std::sin(gsl_vector_get(params, PPG * data->nGauss)) + 1)/2;
			gsl_vector_set(diff, i, sum - data->pdf[i]); //wrap the uniform level in a sine, to map (-inf, inf) to (0,1)
			//TODO: change scaling from [0,1] to [0,1/MAX_ANGLE]
		}
		return GSL_SUCCESS;
	}

	//evaluate the derivative of the gaussian functions for a particular parameter set
	int GaussSumDerivFunc(const gsl_vector *params, void *data_ptr, gsl_matrix *jac)
	{
		struct gauss_data *data = (struct gauss_data *) data_ptr;
		//set jac[data index][parameter index] to be d(data i)/d(parameter j)
		for (int i = 0; i < data->nData; i++)
		{
			for (int j = 0; j < data->nGauss; j++)
			{
				const bool is_neg = gsl_vector_get(params, PPG*j) < 0;//false;//
				const double amp = std::abs(gsl_vector_get(params, PPG*j    ));
				const double sig = gsl_vector_get(params, PPG*j + 1);
					  double pos = gsl_vector_get(params, PPG*j + 2);
				pos = fmod(fmod(pos, MAX_ANGLE) + MAX_ANGLE, MAX_ANGLE);
				const double d1 = data->angle[i] - pos;
				const double d2 = d1 - std::copysign(180, d1);
				const double z1 = d1/sig;
				const double z2 = d2/sig;
				const double exp1 = std::exp(-0.5 * z1 * z1);
				const double exp2 = std::exp(-0.5 * z2 * z2);
				gsl_matrix_set(jac, i, PPG*j    , (exp1 + exp2) * (is_neg ? -1 : 1)); //d/d amplitude
				gsl_matrix_set(jac, i, PPG*j + 1, (amp / (sig * sig * sig)) * (d1 * d1 * exp1 + d2 * d2 * exp2)); //d/d sigma
				gsl_matrix_set(jac, i, PPG*j + 2, (amp / (sig * sig      )) * (d1 *      exp1 + d2 *      exp2)); //d/d position
			}
			if (data->use_horizontal) gsl_matrix_set(jac, i, PPG * data->nGauss, std::cos(gsl_vector_get(params, PPG*data->nGauss))/2);
		}
		return GSL_SUCCESS;
	}

	//fit the gaussians, using a set of initial values, the data values to fit against (pdf) and the angle the data is evaluated at
	AngleDistribution FitGaussians(std::vector<crossing_location_type> &initial_positions, arma::vec &pdf, arma::vec &angle, bool use_horizontal)
	{
		const int nData = pdf.size();
		const int nGaussians = initial_positions.size();
		const int nParams = PPG * nGaussians + (use_horizontal ? 1 : 0); //each gaussian has an amplitude, sigma, and location

		struct gauss_data data(nData, nGaussians, angle, pdf, use_horizontal);

		//set initial values
		gsl_vector *initial_guess = gsl_vector_alloc(nParams);
	// 	std::cout << "initial guess:" << std::endl;
		for (int i = 0; i < nGaussians; i++)
		{
			const int p1 = initial_positions[i].first, p2 = initial_positions[i].second;
			const double a1 = angle[p1], a2 = angle[p2];
			double amp, size, pos;
			if (p1 < p2)
			{
				//position is in between points 1 and 2
				amp = pdf[(p1 + p2)/2];
				size = a2 - a1;
				pos = (a1 + a2)/2;
			} else {
				//position is in the wraparound region
				amp = pdf[ ((p1 + p2 + pdf.size()) / 2) % pdf.size() ];
				size = MAX_ANGLE - abs(a1 - a2);
				pos = fmod((a1 + a2 + MAX_ANGLE)/2, MAX_ANGLE);
			}
			gsl_vector_set(initial_guess, PPG*i  ,  amp);
			gsl_vector_set(initial_guess, PPG*i+1, size);
			gsl_vector_set(initial_guess, PPG*i+2,  pos);
	// 		std::cout << a1 << ", " << a2 << " -> amp = " << amp << ", size = " << size << ", pos = " << pos << std::endl;
		}
		//if using a horizontal line, the initial guess is the minimum value of the pdf
		if (use_horizontal) gsl_vector_set(initial_guess, nParams - 1, std::asin(2*pdf.min() - 1));

		gsl_multifit_nlinear_fdf fdf;
		fdf.f = GaussSumFunc;
		fdf.df = GaussSumDerivFunc;
		fdf.fvv = NULL;
		fdf.n = nData;
		fdf.p = nParams;
		fdf.params = &data;

		const gsl_multifit_nlinear_type *type = gsl_multifit_nlinear_trust;
		gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
		gsl_multifit_nlinear_workspace *w = gsl_multifit_nlinear_alloc(type, &fdf_params, nData, nParams);
		gsl_multifit_nlinear_init(initial_guess, &fdf, w);

		const int max_iterations = 1000;
		const double xtol = std::pow(10, -4); //tolerance by change in independent variables. should be 10^-d for d digits of precision
		const double gtol = std::pow(GSL_DBL_EPSILON, 1/3.0); //tolerance by norm of gradient
		const double ftol = 0.0; // tolerance by change in cost function
		int info = 0, status = 0;

		//do the fit
		status = gsl_multifit_nlinear_driver(max_iterations, xtol, gtol, ftol, NULL, NULL, &info, w); //the NULL's are for the callback params

		if (status != GSL_SUCCESS)
		{
			std::cerr << "Warning: gaussian fitting failed, status = " << status << " -> " << gsl_strerror(status) << ", info = " << info << std::endl;
		}

		std::vector<std::tuple<double, double, double>> out_params(nGaussians);
		for (int i = 0; i < nGaussians; i++)
		{
			const double amp  =  std::abs(gsl_vector_get(w->x, PPG*i    ));
			const double size =  std::abs(gsl_vector_get(w->x, PPG*i + 1));
			const double pos  = fmod(fmod(gsl_vector_get(w->x, PPG*i + 2), MAX_ANGLE) + MAX_ANGLE, MAX_ANGLE);
			const double p = amp * std::sqrt(2*M_PI) * size;
			out_params[i] = std::make_tuple(p, size, pos);
		}

		AngleDistribution return_value;

		if (use_horizontal)
		{
			const double uniform_prob = MAX_ANGLE * (std::sin(gsl_vector_get(w->x, nParams - 1)) + 1)/2;
			return_value =  AngleDistribution(uniform_prob, out_params);
		} else {
			return_value = AngleDistribution(out_params);
		}

		gsl_vector_free(initial_guess);
		gsl_multifit_nlinear_free(w);

		return return_value;
	}

	//evaluate the sum of gaussians model at a particular angle
	double EvaluateGaussianSum(std::vector<std::tuple<double, double, double>> &gaussians, double x)
	{
	//     return -1; //filler until I replace everywhere this function is used
		AngleDistribution ad(gaussians);
		return ad.EvaluateDistribution(x);
	}

	//fit gaussians to a KDE, they are in the format <amplitude, sigma, position/angle>
	AngleDistribution FitGaussiansWraparound(std::vector<std::pair<double, double>> &KDE, std::vector<double> &fault_angles, const double param_penalty = 2, const int max_gaussians = 10)
	{
	// 	std::cout << "starting gaussian fit" << std::endl;
		arma::vec pdf(KDE.size());
		arma::vec angle(KDE.size());
		for (unsigned int i = 0; i < KDE.size(); i++)
		{
			angle[i] = KDE[i].first;
			pdf[i] = KDE[i].second;
		}

		//cout << "fitting a kde for " << fault_angles.size() << " faults" << std::endl;

	// 	const double err_thresh = 1e-6;//1e-3;

		arma::cx_vec ft = arma::fft(pdf); //conveniently, our data does actually wrap around, so we don't need a taper function
		//the frequency order isn't specified in the arma documentation, assuming that it's 0-freq -> +ve freqs -> most negative freqs -> -1/N
		arma::vec freqs(ft.size());
		for (unsigned int i = 0; i < freqs.size(); i++)
		{
			freqs[i] = i + (i <= freqs.size() / 2 ? 0 : -freqs.size()) / (double) freqs.size();
		}

		int nGaussians = 0;

		arma::vec residual(KDE.size()); //residuals, used to calculate the error

		AngleDistribution results, results_unif, previous_results;

		std::vector<std::pair<double, decltype(results)>> ic; //information criteria - used to chosse the number of gaussians

		while (nGaussians < max_gaussians)
		{
			nGaussians += 1;

			std::vector<crossing_type> base_crossings = FindZerocrossings(ft, freqs);

			//now make scale image
			const int max_peaks = base_crossings.size() / 2;
			if (max_peaks < nGaussians) break; //we don't have enough peaks to get the required amount, so use the best result before this
			int peak_count = max_peaks;
			std::vector<std::vector<crossing_type>> scale_image = {base_crossings};

			double scale = 0;
			while (peak_count > nGaussians)
			{
				double next_scale = scale + 1;
				while (true)
				{
					arma::cx_vec scaled_ft(ft.size());
					for (unsigned int i = 0; i < scaled_ft.size(); i++) scaled_ft[i] = ft[i] * std::sqrt(M_PI * next_scale) * std::exp(-M_PI * M_PI * freqs[i] * freqs[i] * next_scale);
					std::vector<crossing_type> next_crossings = FindZerocrossings(scaled_ft, freqs);
					int this_peaks = next_crossings.size() / 2;
	// 				std::cout << "At scale " << next_scale << " there are " << this_peaks << " peaks, diff = " << peak_count - this_peaks << std::endl;
					if (peak_count - this_peaks <= 1){
						scale = next_scale;
						scale_image.push_back(next_crossings);
						peak_count = this_peaks;
						break;
					}
					next_scale = (next_scale + scale) / 2;
				}
			}

			//now follow scale image back to original scale
			std::vector<crossing_location_type> positions;
		//	positions = follow_scale_image(scale_image, KDE.size(), nGaussians);
			positions = SimpleLocationDetection(scale_image.back());//base_crossings

			//now use the initial values to fit the (sum of) gaussians properly
			//debug output
	//  		for (auto it = positions.begin(); it < positions.end(); it++)
	//  		{
	//  			std::cout << "Gaussian " << it - positions.begin() << " is at " << angle[it->first] << ", " << angle[it->second] << "; ";
	//  		}
	//  		std::cout << std::endl;

			results = FitGaussians(positions, pdf, angle, false);
			results_unif = FitGaussians(positions, pdf, angle, true);

			for (int i = 0; i < (int)residual.size(); i++) residual[i] = pdf[i] - results.EvaluateDistribution(angle[i]);
			double rss = 0;
			for (int i = 0; i < (int)residual.size(); i++) rss += residual[i]*residual[i];
	//  		double err = arma::norm(residual, 2);
			if (nGaussians >= max_peaks) ft = arma::fft(residual); //If there aren't any more peaks to find, use the residuals of the fit with the current data

			double llikelihood = 1; //calculate the log likelihood
			double llikelihood_unif = 1;
			for (auto it = fault_angles.begin(); it < fault_angles.end(); it++)
			{
				double likelihood = results.EvaluateDistribution(*it);
				llikelihood += likelihood > 0 ? std::log(likelihood) : 0;
				llikelihood_unif += std::log(results_unif.EvaluateDistribution(*it));
			}
			int nParams = PPG * nGaussians; //number of free parameters
			int nParams_unif = nParams + 1;
			//2*nParams + 2 * std::log(err);// - 2 * llikelihood;
			//(corrected) akaike information criterion
			double aicc = param_penalty*nParams /*+ (2*nParams*(nParams + 1)) / (double)(PPG*fault_angles.size() - nParams - 1)*/ - 2*llikelihood;
			double aicc_unif = param_penalty * nParams_unif - 2 * llikelihood_unif;
			bool has_negative = false;
			for (auto it = results.gaussians.begin(); it < results.gaussians.end(); it++) if (std::get<0>(*it) < 0) has_negative = true; //reject fittings with negative gaussians, all the groups should be positive
			std::cout << "At ng " << nGaussians << " the log likelihood is " << llikelihood << " with AICC = " << aicc << ", unif ll is " << llikelihood_unif << ", unif aicc = " << aicc_unif << ", has neg: " << has_negative << std::endl;

	// 		 //error err
	//         has_negative = false; //debug
			if (!std::isnan(aicc) && !has_negative) ic.push_back(std::make_pair(aicc, results));
			if (!std::isnan(aicc_unif)) ic.push_back(std::make_pair(aicc_unif, results_unif));

		}
	// 	return ic[5-1].second; //debug
		decltype(ic)::iterator result_it = std::min_element(ic.begin(), ic.end(),
			[](const decltype(ic)::value_type &a, const decltype(ic)::value_type &b){ return a.first < b.first;}); //get the set of Gaussians with the minimum aicc value
		return result_it->second; //previous_results
	}


	AngleDistribution KdeEstimationStrikes(VECTOR &lines, const double param_penalty)
	{
		std::vector<line_type> &lineaments = lines.data;
		std::string out_name = FracG::AddPrefixSuffixSubdirs(lines.out_path, {stats_subdir}, "angle_distribution_KDE", ".tsv");
		std::ofstream txtF = FracG::CreateFileStream(out_name);
		int index = 0 ;
		arma::vec ANGLE(lineaments.size(), arma::fill::zeros);
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
		std::sort(ANGLE.begin(), ANGLE.end());

		std::vector< std::pair<double, double>> GAUSS;// =  kde( ANGLE, ANGLE.size(), 10); //these are (angle, pdf) std::pairs
		arma::vec est_pdf = AngleKdeFillArray(MAX_ANGLE * 10, ANGLE, 10);
		for (unsigned int i = 0; i < est_pdf.size(); i++)
		{
			GAUSS.push_back(std::make_pair(MAX_ANGLE * i / (double) est_pdf.size(), est_pdf[i]));
		}
		std::vector< std::pair<double, double>> Maximas;

		MovingAverageFilterWraparound(GAUSS, 5);

		for (unsigned int i = 0; i < GAUSS.size(); i++)
		{
			if (i == 0)
				if (GAUSS[i].second > GAUSS[GAUSS.size()-1].second && 
					GAUSS[i].second > GAUSS[i+1].second)
						Maximas.push_back(std::make_pair(GAUSS[i].first, GAUSS[i].second));

			if (i > 0 && i+1 < GAUSS.size())
				if (GAUSS[i].second > GAUSS[i-1].second &&
					GAUSS[i].second > GAUSS[i+1].second)
						Maximas.push_back(std::make_pair(GAUSS[i].first, GAUSS[i].second));


			if (i == GAUSS.size())
				if (GAUSS[i].second > GAUSS[i-1].second && 
					GAUSS[i].second > GAUSS[0].second)
						Maximas.push_back(std::make_pair(GAUSS[i].first, GAUSS[i].second));
		}

		std::vector<double> angles_vector(ANGLE.begin(), ANGLE.end());
		AngleDistribution angle_dist = FitGaussiansWraparound(GAUSS, angles_vector, param_penalty);
		gauss_params gauss_p = angle_dist.gaussians;

		txtF << "Amplitude\tSigma\tMean" << std::endl;
		std::cout << "We fit " << gauss_p.size() << " gaussians to the angle data, with parameters:" << std::endl << "Amplitude\tSigma\tMean" << std::endl;
		for (auto it = gauss_p.begin(); it < gauss_p.end(); it++)
		{
			txtF << std::get<0>(*it) << "\t" << std::get<1>(*it) << "\t" << std::get<2>(*it) << std::endl;
			std::cout << std::get<0>(*it) << "\t" << std::get<1>(*it) << "\t" << std::get<2>(*it) << std::endl;
		}

		std::cout << "Uniform Level:\t" << angle_dist.uniform_prob << std::endl;
		txtF      << "Uniform Level:\t" << angle_dist.uniform_prob << std::endl;
		std::cout << std::endl;

		txtF << "Fault sets:" << "\t" << gauss_p.size() << std::endl;
		txtF << "Angle \t Density \t Smoothed Density" << std::endl;
		for(unsigned int i =0; i < GAUSS.size(); i++)
			txtF << GAUSS[i].first << "\t " << GAUSS[i].second << "\t" << angle_dist.EvaluateDistribution(GAUSS[i].first) <<  std::endl;   
		txtF << std::endl;
		//set namespace variable if true
		return AngleDistribution(gauss_p);
	}

	std::vector<double> BootStrapping(std::vector<double>data)
	{
		const int n = data.size();
		std::vector<double> bt_data(data);

		boost::random_device dev;
		boost::mt19937 rng(dev);
		boost::random::uniform_int_distribution<> rand_bt(0,n-1);

		for (int i = 0; i < n; i++)
			bt_data[i] = data[rand_bt(rng)];

		return(bt_data);
	}

	std::vector<double> PearsonCorrelation(std::vector<double>data1, std::vector<double>data2, size_t nb_iter)
	{
		assert(data1.size() == data2.size());

		double init_pearson;
		const size_t stride = 1;
		size_t n = data1.size();
		std::vector <double> BootsTrapping_Cor;
		BootsTrapping_Cor.reserve(nb_iter + 1);

		gsl_vector_const_view D1 = gsl_vector_const_view_array( &data1[0], data1.size() );
		gsl_vector_const_view D2 = gsl_vector_const_view_array( &data2[0], data2.size() );

		//first correlation coefficinet without bootstrapping
		BootsTrapping_Cor.push_back(gsl_stats_correlation((double*) D1.vector.data, stride,
											   (double*) D2.vector.data, stride,
												n ));

		 //now perform bootstrapping analysis for iterations = nb_iter                                     
		for (int i = 0; i < nb_iter; i++)
		{
			std::vector<double> bt_data1 = BootStrapping(data1);
			std::vector<double> bt_data2 = BootStrapping(data2);

			gsl_vector_const_view bt_D1 = gsl_vector_const_view_array( &bt_data1[0], bt_data1.size() );
			gsl_vector_const_view bt_D2 = gsl_vector_const_view_array( &bt_data2[0], bt_data2.size() );

			BootsTrapping_Cor.push_back(gsl_stats_correlation((double*) bt_D1.vector.data, stride,
											   (double*) bt_D2.vector.data, stride,
																						n ));
		}
		return(BootsTrapping_Cor);
	}

	//what does this do? //I think the chooses the id of the gaussians that 
	int CheckGaussians(AngleDistribution &angle_dist, double angle)
	{
		std::vector<std::pair<int, double>> p_classes;
		double pi = boost::math::constants::pi<double>();

		if (angle_dist.gaussians.size() > 1)
		{
			int g;
			double P = 0;
			for (int i = 0; i < angle_dist.gaussians.size(); i++)
			{
	// 			std::tuple<double, double, double> params1 = angle_dist.gaussians.at(i);
	// 			double P = 1 / (std::get<1>(params1)* 2* pi ) * exp(- (  pow(angle - std::get<2>(params1), 2)) / (2 * std::get<1>(params1))  );
				double P = AngleDistribution::EvaluateGaussian(angle_dist.gaussians.at(i), angle);
				p_classes.push_back(std::make_pair(i, P));
			}
			std::pair<int, double> G = *max_element(p_classes.begin(), p_classes.end(), [](const auto& p1, const auto& p2) 
			{ 
				return p1.second < p2.second; 
			});
			if (angle_dist.with_uniform && angle_dist.uniform_prob > G.second) return -1; //the the most likely class is the uniform distrubution
			return(G.first);
		}
		else
			return(0);
	}

	line_type RandomLine(box AOI, polygon_type t_AOI, double angle)
	{
		line_type basis, r_line;
		boost::random_device dev;
		boost::mt19937 rng(dev);

		double min_x = geometry::get<geometry::min_corner, 0>(AOI);
		double min_y = geometry::get<geometry::min_corner, 1>(AOI);
		double max_x = geometry::get<geometry::max_corner, 0>(AOI);
		double max_y = geometry::get<geometry::max_corner, 1>(AOI);

		point_type ll(min_x, min_y);
		point_type lr(max_x, min_y);

		geometry::append(basis, ll);
		geometry::append(basis, lr);

		double basis_strike = (float)(atan2(basis.front().x() - basis.back().x(), basis.front().y() - basis.back().y())) 
									* (180 / math::constants::pi<double>());
		if (basis_strike  < 0) 
			basis_strike  += 180;

		boost::random::uniform_real_distribution<double> rand_x(min_x, max_x);
		boost::random::uniform_real_distribution<double> rand_y(min_y, max_y);

		bool search = true;
		double divide = 2;
		int iter = 0;
		do{
			boost::random::uniform_real_distribution<double> rand_l(geometry::length(basis)/divide, geometry::length(basis)-1);
			r_line.clear();
			double x1  = rand_x(rng);
			double y1  = rand_y(rng);

			if(geometry::within(point_type(x1, y1), t_AOI))
			{
				double len = rand_l(rng);

				double x2 = x1 + len * cos(angle + 90 * math::constants::pi<double>() / 180.0);
				double y2 = y1 + len * sin(angle + 90 * math::constants::pi<double>() / 180.0);

				if(!geometry::disjoint(point_type(x2, y2), t_AOI))
				{
					geometry::append(r_line, point_type(x1, y1));
					geometry::append(r_line, point_type(x2, y2));
					search = false;
				}

				if(geometry::disjoint(point_type(x2, y2), t_AOI));
					iter++;

				if (iter > 100)
				{
					divide += 1;
					iter = 0;
				}
			}
		} while(search);
		return(r_line);
	}

	void ScanLine(VECTOR &lines, int nb_scanlines, AngleDistribution &angle_dist)
	{
		if (angle_dist.gaussians.size() <= 0)
		{
			std::cout << "Cannot set scan line orientations (no data on prinipal orientations)" << std::endl;
			return;
		}
		std::vector <double> density;
		std::vector<std::pair<int, double>> direc_nb;
		std::ofstream txtF = FracG::CreateFileStream(lines.out_path / stats_subdir / ("scaline_analysis.tsv"));
		txtF << lines.name << std::endl;

		for (auto it = angle_dist.gaussians.begin(); it < angle_dist.gaussians.end(); it++)
		{
			float frac = std::ceil((std::get<0>(*it) - 0.05) * 10.0) / 10.0; 
			if (frac > 0)
				direc_nb.push_back(std::make_pair(frac * nb_scanlines, std::get<2>(*it)));
		}
		std::cout << "Scanline analysis for " << direc_nb.size() << " orientations " << std::endl;
		txtF << "Orientation \t Scanline begin \t Scaline end \t Length \t Intersection number \t Intesity \t Spacing"<< std::endl;

		polygon_type t_AOI = ReturnTightAOI(lines.data);
		box AOI = ReturnAOI(lines.data);

		int g = 0;
		for (auto it = direc_nb.begin(); it < direc_nb.end(); it++)
		{
			for (int i = 0; i < it->first;)
			{
				line_type scanline = RandomLine(AOI, t_AOI, it->second);
				std::vector<point_type> output;
				for (auto l = lines.data.begin(); l < lines.data.end(); l++)
				{
					line_type line = *l;
					geometry::intersection(*l, scanline, output);
				}
				if (output.size() > 0)
				{
					density.push_back(output.size() / geometry::length(scanline));
					txtF 	<< g << "\t" << std::setprecision(25)	<< scanline.front().x() << ", " 
							<< scanline.front().y()		<< "\t" << scanline.back().x() 
							<< ", " << scanline.back().y() << "\t" << std::setprecision(10)
							<< geometry::length(scanline) << "\t" << output.size() << "\t"
							<< output.size()/geometry::length(scanline) <<"\t";

					if (output.size() > 1)
					{
						std::vector<std::pair<point_type, double>> distances;
						for (auto cross = output.begin(); cross <output.end(); cross++)
							distances.push_back(std::make_pair(*cross, geometry::distance(scanline.front(), *cross)));

						std::sort(distances.begin(), distances.end(), [](const std::pair<point_type,double> &left, const std::pair<point_type, double> &right) {
							return left.second < right.second;});


						double spaceing = 0;
						for (auto d = distances.begin(); d < distances.end()-1;d++)
						{
							auto nx = std::next(d, 1);
							point_type p = d->first;
							point_type pp = nx->first;
							spaceing += geometry::distance(p, pp);
						}
						txtF << spaceing/(distances.size()-1) << std::endl;
					}
					else 	
						txtF << " " << std::endl;
				i++;
				}
			}
			g++;
		}

		std::cout << "Analysed " << nb_scanlines << " scalines \n" << std::endl;
	}

	//create statisics from input file--------------------------------------
	void CreateStats(VECTOR &lines, AngleDistribution &angle_dist)
	{
		int FaultNumber = lines.data.size();
		std::vector<line_type> faults = lines.data;
		point_type centre;

		std::vector<double> Angle, Length, Sinuosity;
		Angle.reserve(faults.size());
		Length.reserve(faults.size());
		Sinuosity.reserve(faults.size());

		std::vector<std::string> Coodinates;
		std::vector<p_index> points;
		std::vector<pl_index> points2;

		unsigned int index = 0;
		BOOST_FOREACH(line_type F, faults)
		{ 
			double strike = (float)(atan2(F.front().x() - F.back().x(), F.front().y() - F.back().y())) 
								* (180 / math::constants::pi<double>());
			if (strike  < 0) 
				strike  += 180;
			double length = (double) geometry::length(F) / 1000;
			double sinuosity =   geometry::length(F) / geometry::distance(F.front(), F.back());

			Angle.push_back(strike);
			Length.push_back(length);
			Sinuosity.push_back(sinuosity);

			geometry::unique(F);
			std::string coord = "";
			int p_nb = 0;
			BOOST_FOREACH(point_type p, F)
			{
				if (p_nb < F.size()-1)
				{
					coord += to_string(p.x()) + " ";
					coord += to_string(p.y()) + ", ";
				}
				else 
				{
					coord += to_string(p.x()) + " ";
					coord += to_string(p.y()) + "";
				}
			p_nb++;
			}
			Coodinates.push_back(coord);
			coord.clear();
			geometry::centroid(F, centre);
			points.push_back(std::make_pair(centre, index));
			points2.push_back(std::make_tuple(centre, index, geometry::length(F)));
			index++;
		}

	/***********************************************************************
	* Bour, O., & Davy, P. (1999). 
	* Clustering and size distributions of fault patterns: Theory and measurements. 
	* Geophysical Research Letters, 26(13), 2001-2004.
	***********************************************************************/
	//fault centre to  centre distance--------------------------------------
		std::vector<p_index> closest;
		PointTree(points, closest);
		arma::vec distance(points.size(), arma::fill::zeros);

	//fault centre to centre distance to larger larger fault----------------------
		std::vector<pl_index> closest2;
		PointTree2(points2, closest2, floor(arma::vec(Length).max()*1000));
		arma::vec distance2(points.size(), arma::fill::zeros);

	//write data to file----------------------------------------------------
		std::ofstream txtF = FracG::CreateFileStream(lines.out_path / stats_subdir / ("vector_properties.tsv"));//statistics
		txtF << lines.name <<std::endl;
		if (angle_dist.gaussians.size() > 0)
		{
			int i = 0;
			txtF << "No \t Amplitude \t Stddev \t Mean" << std::endl;
			for (auto it = angle_dist.gaussians.begin(); it < angle_dist.gaussians.end(); it++)
			{
				txtF << i << "\t" << std::get<0>(*it) << "\t" << std::get<1>(*it) << "\t" << std::get<2>(*it) << std::endl;
				i++;
			}
		}
		txtF << "index \t strike \t length \t sinuosity \t Generation \t nearest line" 
			 <<"\t distance to nearest line \t nearest larger line \t distance to nearest larger line \t coodinates" << std::endl;

		size_t curPos = 0;
		for (int i = 0; i< lines.data.size(); i++)
		{
			txtF << i << "\t" << Angle.at(i) << "\t" <<  Length.at(i) << "\t" << Sinuosity.at(i) << "\t" << CheckGaussians(angle_dist, Angle.at(i)) << "\t";
			if (lines.data.size() != 1)
			{
				 txtF << closest[curPos].second  << "\t" << geometry::distance(points[curPos].first, closest[curPos].first) << "\t"
					  << get<1>(closest2[curPos]) <<"\t" << geometry::distance(points[curPos].first, get<0>(closest2[curPos])) << "\t" << Coodinates.at(i) << std::endl;

			distance[i]  = geometry::distance(points[curPos].first, closest[curPos].first);
			distance2[i] = geometry::distance(points[curPos].first, get<0>(closest2[curPos]));
			}
			else
				txtF << "\t"  << "\t" << "\t" << "\t" << std::endl;
			curPos++;
		}
		txtF << "Mean:     " << "\t" << arma::mean(arma::vec(Angle)) 	<< "\t" << arma::mean(arma::vec(Length))	<< "\t" << arma::mean(arma::vec (Sinuosity)) 	<< "\t" << "\t" << "\t" << arma::mean(distance) 	<< "\t" << "\t" << arma::mean(distance2)   		<< "\n"
			 << "Median:   " << "\t" << arma::median(arma::vec(Angle)) 	<< "\t" << arma::median(arma::vec(Length))	<< "\t" << arma::median(arma::vec (Sinuosity)) 	<< "\t" << "\t" << "\t" << arma::median(distance)	<< "\t" << "\t" << arma::median(distance2) 		<< "\n"
			 << "Stddev:   " << "\t" << arma::stddev(arma::vec(Angle)) 	<< "\t" << arma::stddev(arma::vec(Length))	<< "\t" << arma::stddev(arma::vec (Sinuosity)) 	<< "\t" << "\t" << "\t" << arma::stddev(distance) 	<< "\t" << "\t" << arma::stddev(distance2) 		<< "\n"
			 << "Variance: " << "\t" << arma::var(arma::vec(Angle))  	<< "\t" << arma::var(arma::vec(Length))		<< "\t" << arma::var(arma::vec (Sinuosity))  	<< "\t" << "\t" << "\t" << arma::var(distance)		<< "\t" << "\t" << arma::var(distance2)			<< "\n"
			 << "Range:    " << "\t" << arma::range(arma::vec(Angle))	<< "\t" << arma::range(arma::vec(Length))	<< "\t" << arma::range(arma::vec (Sinuosity))	<< "\t" << "\t" << "\t" << arma::range(distance)	<< "\t" << "\t" << arma::range(distance2)	<< "\n"
			 << std::endl;
		txtF.close(); 

		txtF = FracG::CreateFileStream(lines.out_path / stats_subdir / ("length_orientations_correlations.tsv"));

	txtF << "\n Linear correlations "
		 << "\n Parameters " << "\t" << " Correlation Coefficient "<< "\t" << "StdDev" << std::endl; 

	//testing for linear correlations---------------------------------------
		std::vector<double> p_dist  = arma::conv_to< std::vector<double> >::from(distance);
		std::vector<double> p_dist2 = arma::conv_to< std::vector<double> >::from(distance2);

		std::vector<double> len_sin_cor = PearsonCorrelation(Length, Sinuosity, 50);
		txtF << "length sinuosity: \t" << arma::mean(arma::vec((len_sin_cor))) << "\t" << arma::stddev(arma::vec (len_sin_cor)) << std::endl;

		std::vector<double> len_orient_cor = PearsonCorrelation(Length, Angle, 50);
		txtF << "length orientation: \t" << arma::mean(arma::vec((len_orient_cor))) << "\t" << arma::stddev(arma::vec (len_orient_cor)) << std::endl;

		std::vector<double> orient_sin_cor = PearsonCorrelation(Length, Angle, 50);
		txtF << "orientation sinuosity: \t" << arma::mean(arma::vec((orient_sin_cor))) << "\t" << arma::stddev(arma::vec (orient_sin_cor)) << std::endl;

		std::vector<double> len_p_dist_cor = PearsonCorrelation(Length, p_dist, 50);
		txtF << "length distance: \t" << arma::mean(arma::vec((len_p_dist_cor))) << "\t" << arma::stddev(arma::vec (len_p_dist_cor)) << std::endl;

		std::vector<double> len_p_dist_cor2 = PearsonCorrelation(Length, p_dist2, 50);
		txtF << "length distanse (larger): \t" << arma::mean(arma::vec((len_p_dist_cor2))) << "\t" << arma::stddev(arma::vec (len_p_dist_cor2)) << std::endl;

		std::cout << "Wrote geometric statisics to file" << std::endl; 
	}

	int GetCluster(std::pair<double, double> point, std::vector<std::pair<double,double>> clusters)
	{
		//checking the euclidian distance between the data point and every centroid
		int  n = 0, c;
		double d = std::numeric_limits<double>::max();

		for (auto i : clusters)
		{
			double D = sqrt( pow((point.first - i.first) ,2) + pow((point.second - i.second) ,2));
			if ( D < d)
			{
				d = D; 
				c = n;
			}
			n++;
		}
		return(c);
	}


	void KMCluster(bool output, VECTOR &lines, AngleDistribution &angle_dist)
	{
		int No = angle_dist.gaussians.size();
		if (No > 1)
		{
			std::cout <<"k-means clustering of geometric properties" << std::endl;
			std::ofstream txtF;
			std::vector<double> Angle, Length, Sinuosity;
			BOOST_FOREACH(line_type F, lines.data)
			{ 
				double strike = (float)(atan2(F.front().x() - F.back().x(), F.front().y() - F.back().y())) 
									* (180 / math::constants::pi<double>());
				if (strike  < 0) 
					strike  += 180;

				double length = (double) geometry::length(F) / 1000;
				double sinuosity =   geometry::length(F) / geometry::distance(F.front(), F.back());

				Angle.push_back(strike);
				Length.push_back(length);
				Sinuosity.push_back(sinuosity);
			}

			bool status;
			arma::mat means;
			arma::mat a = arma::conv_to<arma::rowvec>::from(Angle);
			arma::mat l = arma::conv_to<arma::rowvec>::from(Length);
			arma::mat s = arma::conv_to<arma::rowvec>::from(Sinuosity);

	//Angle and Length------------------------------------------------------
			auto ang_len = arma::join_cols(a,l);
			status = kmeans(means, ang_len, No, arma::random_subset, 20, output);
			if(status == false) 
				std::cout << "clustering angle and lenght failed" << std::endl;
			else
			{
				txtF = FracG::CreateFileStream(lines.out_path / stats_subdir / ("kmeans_clustering_angle_length.tsv"));
				txtF << lines.name << std::endl;
				txtF << "km-means clustering of orientation and length" << std::endl;
				txtF << "Amplitude\tSigma\tMean" << std::endl;
				for (auto it = angle_dist.gaussians.begin(); it < angle_dist.gaussians.end(); it++)
					txtF << std::get<0>(*it) << "\t" << std::get<1>(*it) << "\t" << std::get<2>(*it) << std::endl;
				txtF << "Cluster centroids" << std::endl;
				std::vector<std::pair<double, double>> clusters;
				int j = 0;
				for (int i = 0; i < means.size()-1;i++) 
				{
					j = i + 1;
					txtF << j-1 << "\t" << means[i] << "\t" << means[j] << std::endl;
					clusters.push_back(std::make_pair(means[i], means[i++]));
				}

				txtF << "No \t Angle \t Length" << std::endl;
				for (int i = 0; i < lines.data.size(); i++)
					txtF << i << "\t" << Angle[i] << "\t" << Length[i] << "\t" << GetCluster(std::make_pair(Angle[i], Length[i]), clusters) << std::endl;
				means.clear();
			}

	//Angle and Sinuosity----------------------------------------------------
			auto ang_sin = arma::join_cols(a,s);
			status = arma::kmeans(means, ang_sin, No, arma::random_subset, 20, output);
			if(status == false) 
				std::cout << "clustering angle and lenght failed" << std::endl;
			else
			{
				txtF = FracG::CreateFileStream(lines.out_path / stats_subdir / ("kmeans_clustering_angle_sinuosity.tsv"));
				txtF << lines.name << std::endl;
				txtF << "km-means clustering of orientation and sinuosity" << std::endl;
				txtF << "Amplitude\tSigma\tMean" << std::endl;
				for (auto it = angle_dist.gaussians.begin(); it < angle_dist.gaussians.end(); it++)
					txtF << std::get<0>(*it) << "\t" << std::get<1>(*it) << "\t" << std::get<2>(*it) << std::endl;
				txtF << "Cluster centroids" << std::endl;

				std::vector<std::pair<double, double>> clusters;
				int j = 0;
				for (int i = 0; i < means.size()-1;i++) 
				{
					j = i + 1;
					txtF << j-1 << "\t" << means[i] << "\t" << means[j] << std::endl;
					clusters.push_back(std::make_pair(means[i], means[i++]));
				}

				txtF << "No \t Angle \t Sinuosity" << std::endl;
				for (int i = 0; i < lines.data.size(); i++)
					txtF << i << "\t" << Angle[i] << "\t" << Sinuosity[i] << "\t" << GetCluster(std::make_pair(Angle[i], Sinuosity[i]), clusters) << std::endl;
				means.clear();
			}

	//Length and Sinuosity--------------------------------------------------
			auto len_sin = arma::join_cols(l,s);
			status = arma::kmeans(means, len_sin, No, arma::random_subset, 20, output);
			if(status == false) 
				std::cout << "clustering angle and lenght failed" << std::endl;
			else
			{
				txtF = FracG::CreateFileStream(lines.out_path / stats_subdir / ("kmeans_clustering_length_sinuosity.tsv"));
				txtF << lines.name << std::endl;
				txtF << "km-means clustering of length and sinuosity" << std::endl;
				txtF << "Amplitude\tSigma\tMean" << std::endl;
				for (auto it = angle_dist.gaussians.begin(); it < angle_dist.gaussians.end(); it++)
					txtF << std::get<0>(*it) << "\t" << std::get<1>(*it) << "\t" << std::get<2>(*it) << std::endl;
				txtF << "Cluster centroids" << std::endl;
				std::vector<std::pair<double, double>> clusters;
				int j = 0;
				for (int i = 0; i < means.size()-1;i++) 
				{
					j = i + 1;
					txtF << j-1 << "\t" << means[i] << "\t" << means[j] << std::endl;
					clusters.push_back(std::make_pair(means[i], means[i++]));
				}
				txtF << "No \t Angle \t Sinuosity" << std::endl;
				for (int i = 0; i < lines.data.size(); i++)
					txtF << i << "\t" << Length[i] << "\t" << Sinuosity[i] << "\t" << GetCluster(std::make_pair(Length[i], Sinuosity[i]), clusters) << std::endl;
				means.clear();
			}
		}
		else
			std::cout << "Insufficinet number of clusters derived from orientation data" << std::endl;
		std::cout << " done \n" << std::endl;
	}

}
