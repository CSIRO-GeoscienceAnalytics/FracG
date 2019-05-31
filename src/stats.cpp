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
		name		= n;
		values		= v;
		memcpy(transform, t, 8*sizeof(double));
	 }
 };
 
 vector <READ> input;
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
 Xmax = input[i].transform[0] + input[i].transform[1] * input[i].transform[6] - 1000;    // west
 Xmin = input[i].transform[0] + 1000; 									  // east
 Ymax = input[i].transform[3] - 1000; 									 // north
 Ymin = input[i].transform[3] + input[i].transform[5] * input[i].transform[7] + 1000;	//south
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
		if(geometry::within(F.front(), pl) &&
		geometry::within(F.back(),pl))
		{
		accept = false;
		radius = 2 * input[i].transform[1];

		strike = (float)(atan2(F.front().x() - F.back().x(), F.front().y() - F.back().y())) 
						  * (180 / math::constants::pi<double>());
						  if (strike  < 0) 
							 strike  += 180;
		length = (long double) geometry::length(F) / 1000;
		sinus =   geometry::length(F) / geometry::distance(F.front(), F.back() );

		do{
			data.clear();
			envelop = geom.DefineLineBuffer(F, radius);
			geometry::envelope(envelop, AOI);
			int maxX = (geometry::get<geometry::max_corner, 0>(AOI) - input[i].transform[0]) / input[i].transform[1];
			int minX = (geometry::get<geometry::min_corner, 0>(AOI) - input[i].transform[0]) / input[i].transform[1]; 
			int maxY = abs((geometry::get<geometry::max_corner, 1>(AOI) - input[i].transform[3]) / input[i].transform[5]);
			int minY = abs((geometry::get<geometry::min_corner, 1>(AOI) - input[i].transform[3]) / input[i].transform[5]);

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
								NEW = false;
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
		
		for(int i =0; i < (int) data.size(); i++)
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
 
	if(status == false)
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
	} 
	while(choice == 'n');
	
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
	sinuosity =   geometry::length(F) / geometry::distance(F.front(), F.back() );

	ANGLE(index) = strike;
	LENGTH(index) = length;	
	SINUOSITY(index) = sinuosity;
	
	len_dir.push_back(make_pair(length, strike));
	len_sin.push_back(make_pair(sinuosity, length));
	sin_dir.push_back(make_pair(sinuosity, strike));

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
		if (strike_hist(i) > strike_hist(i-1) 
		&& strike_hist(i) > strike_hist(i+1))
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
 //BoxCount(faults, D);
 //LSFitter(len_dir);
 
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
						if (m_front > m_back)
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
			if (slope.size() >0 )
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

	 for ( typename FSTATS::const_iterator it = Length_Strike.begin(); it != Length_Strike.end(); it++)
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

  do
  {
	double ln_len = var_Ln();
	double ex_len = var_ex();
	double pl_len = (size_t) pow((max_pow - min_pow) * drand48 () + min_pow, 1 / (exponent + 1));

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
 for ( typename FSTATS::const_iterator it = Length_Strike.begin(); it != Length_Strike.end(); it++)
	cout  << boost::math::pdf(myNormal, it->first) << endl;


//ML estimator for lognormal distribution ------------------------------
 double mu_logn_1 = 0, mu_logn = 0, sigma_logn = 0;
 
for (auto& X : length)
 {
	mu_logn_1 += log(X) / length.size();
	mu_logn += mu_logn_1;
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







