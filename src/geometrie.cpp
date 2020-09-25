#include "../include/geometrie.h"
#include "../include/GeoRef.h"
#include "../include/util.h"

const std::string geom_subdir = "geometry";

GEOMETRIE::GEOMETRIE ()

{
}

//create a polygon around a point (square)------------------------------
BUFFER GEOMETRIE::DefineSquareBuffer(point_type POINT, const double Bdistance )
{
	BUFFER buff;
	geometry::strategy::buffer::point_square point_strategy;
	geometry::strategy::buffer::distance_symmetric<double> distance_strategy(Bdistance);
	geometry::strategy::buffer::join_round join_strategy;
	geometry::strategy::buffer::end_round end_strategy;
	geometry::strategy::buffer::side_straight side_strategy;

	geometry::buffer(POINT, buff,
	distance_strategy, side_strategy,
	join_strategy, end_strategy, point_strategy); 
	return buff; 
}

//define a polygon around a point (circle)------------------------------
BUFFER GEOMETRIE::DefinePointBuffer(point_type POINT, const double Bdistance )
{
	BUFFER buff;
	const int points_per_circle = 36;

	geometry::strategy::buffer::distance_symmetric<double> distance_strategy( Bdistance );
	geometry::strategy::buffer::join_round join_strategy(points_per_circle);
	geometry::strategy::buffer::end_round end_strategy(points_per_circle);
	geometry::strategy::buffer::point_circle circle_strategy(points_per_circle);
	geometry::strategy::buffer::side_straight side_strategy;

	geometry::buffer(POINT, buff,
	distance_strategy, side_strategy,
	join_strategy, end_strategy, circle_strategy); 
	return buff; 
}

//define a polygon around line------------------------------------------
BUFFER GEOMETRIE::DefineLineBuffer(line_type fault, const double Bdistance)
{
	BUFFER buff;
	const int points_per_circle = 36;

	geometry::strategy::buffer::distance_symmetric<double> distance_strategy( Bdistance );
	geometry::strategy::buffer::join_round join_strategy(points_per_circle);
	geometry::strategy::buffer::end_round end_strategy(points_per_circle);
	geometry::strategy::buffer::point_circle circle_strategy(points_per_circle);
	geometry::strategy::buffer::side_straight side_strategy;

	boost::geometry::buffer(fault, buff,
					distance_strategy, side_strategy,
					join_strategy, end_strategy, circle_strategy); 
	return buff; 
}

//find minimum spacing between point in  a linestring------------------ 
double GEOMETRIE::minSpacing(line_type Trace)
{
	double len = 0;
	point_type first = Trace.front();
	BOOST_FOREACH(point_type P, Trace)
	{
		if (!geometry::equals(first, P) && (len > geometry::distance(first, P)))
			len = geometry::distance(first, P);	
		first = P;
	}
	return len;
}

line_type GEOMETRIE::ShortestLine(vector <line_type> Set)
{
	double dist = std::numeric_limits<double>::max();
	line_type shortest;
	BOOST_FOREACH(line_type l, Set)
	{
		if (geometry::length(l) < dist && geometry::length(l) != 0)
		{
			dist = geometry::length(l) ;
			shortest = l;
		}
	}
	return(shortest);
}

box GEOMETRIE::ReturnAOI(vector<line_type> lines)
{
	box AOI;
	polygon_type multiL;
	BOOST_FOREACH(line_type F, lines)
		geometry::append(multiL, F);
	
	AOI = boost::geometry::return_envelope<box>(multiL);
	return(AOI);
}

polygon_type GEOMETRIE::Return_tigth_AOI(vector<line_type> lines)
{
	polygon_type hull;
	polygon_type multiL;
	BOOST_FOREACH(line_type F, lines)
		geometry::append(multiL, F);
	
	boost::geometry::convex_hull(multiL, hull);
	return(hull);
}

//return the segment of a line that is the closest to the given point
//the returned iterator points to the linsestring point that is *after* the given point
line_type::iterator get_closest_segment(line_type &Trace, point_type point)
{
	 double best_distance = std::numeric_limits<double>::infinity();
	 line_type::iterator best_index = Trace.begin();
	 point_type prev_point = Trace.front();
	 for (line_type::iterator it = Trace.begin(); it != Trace.end(); it++)
	 {
		 Segment seg(prev_point, *it);
		 double distance = geometry::distance(seg, point);
		 //~ cout << "distance between (" << prev_point.x() << ", " << prev_point.y() << ") -> (" << it->x() << ", " << it->y() << ") and point (" << point.x() << ", " << point.y() << ") is " << distance << endl;
		 if (distance < best_distance)
		 {
			 best_distance = distance;
			 best_index = it;
			 //~ cout << "best point is now (" << best_index->x() << ", " << best_index->y() << ") with distance " << best_distance << endl;
		 }
		 prev_point = *it;
	 }
	 //~ cout << "returning best index at location (" << best_index->x() << ", " <<best_index->y() << ")" << endl;
	 return best_index;
}
 
//returns the portion of the linestring that is between the two points
//assumes that the two points are on the given line_type
//the two points can be given in any order
line_type GEOMETRIE::GetSegment(line_type Trace, point_type Junction, point_type Begin)
{
	box outL;
	line_type line_seg;
	const double threshold = 1;
	point_type first, last;
	
	//determine which line segment each point is located on
	line_type::iterator jit = get_closest_segment(Trace, Junction);
	line_type::iterator bit = get_closest_segment(Trace, Begin);
	
	//sort them so that the next part works no matter which order the points are in, relative to the line
	line_type::iterator first_segment, last_segment;
	if (jit <= bit){
		first = Junction;
		last = Begin;
		first_segment = jit;
		last_segment = bit;
	} else {
		first = Begin;
		last = Junction;
		first_segment = bit;
		last_segment = jit;
	}

	//if they're in the same segment, the new line segment is just between the two points
	//we can probably take this part out
	if (jit == bit){
		geometry::append(line_seg, first);
		geometry::append(line_seg, last);
		return line_seg;
	}
	//the output line segment is the first point -> any and all entire segments that lie between the two points -> from the last whole segment to the end point
	geometry::append(line_seg, first);
	point_type previous = first;
	for (auto it = first_segment; it != last_segment; it++)
	{
		if (geometry::distance(previous, *it) > threshold) 
			geometry::append(line_seg, *it);
		previous = *it;
	}
	geometry::append(line_seg, last);
	return line_seg; 
}
 
//comparison operator that compares a particular index in a (vector of) tuples)
template<int index, template<typename> class F = std::less>
struct TupleCompare
{
	template<typename T>
	bool operator()(T const &t1, T const &t2)
	{
		return F<typename tuple_element<index, T>::type>()(std::get<index>(t1), std::get<index>(t2));
	}
};

//comparison operator that first checks indexa, then uses indexb if the two indexa values are the same
template<int indexa, int indexb, template<typename> class F = std::less, template<typename> class G = std::equal_to>
struct TupleTwoCompare
{
	template<typename T>
	bool operator()(T const &t1, T const &t2)
	{
		bool first_same = G<typename tuple_element<indexa, T>::type>()(std::get<indexa>(t1), std::get<indexa>(t2));
		if (!first_same) return F<typename tuple_element<indexa, T>::type>()(std::get<indexa>(t1), std::get<indexa>(t2));
		return F<typename tuple_element<indexb, T>::type>()(std::get<indexa>(t1), std::get<indexa>(t2));
	}
};

//need a specific comparator, not just selecting an appropriate element from the tuple
struct FaultDistanceCompare
{
	bool operator()(std::tuple<long double, point_type, AttachPoint> const &t1, std::tuple<long double, point_type, AttachPoint> const &t2)
	{
		AttachPoint const &a1 = get<2>(t1);
		AttachPoint const &a2 = get<2>(t2);
		if (a1 != a2) return a1 < a2; //if they attach to different locations, then they go in the order of those locations
		if (a1 == AttachPoint::front) return get<0>(t1) > get<0>(t2); //if they're both on the front, then the point with the *furthest* distance goes first
		return get<0>(t1) < get<0>(t2); //if they're both on the end point or (more likely) middle sectoin, add them increasing order of distance, so smaller distances first
	}
};
 
//sort a vector of points according to their distances along a line
void GEOMETRIE::SortDist(vector <std::tuple<long double, point_type, AttachPoint>>& cross)
{
	std::sort(cross.begin(), cross.end(), FaultDistanceCompare());
}

void GEOMETRIE::CentreDistanceMap(VECTOR lines, float cell_size)
{
	
	GEO georef;
	point_type centre;
	vector<p_index> result;
	geometry::index::rtree<p_index, geometry::index::rstar<16>> DistTree;

	int index = 0 ;

	BOOST_FOREACH(line_type l, lines.data)
	{
		geometry::centroid(l, centre);
		DistTree.insert(make_pair(centre, ++index));
	}
	
//now we nedd to create a georefernece system based on the bounding box
//around the features adn the size of the smallest feature--------------

	box AOI = ReturnAOI(lines.data);
	polygon_type t_AOI = Return_tigth_AOI(lines.data);
	
	double min_x = geometry::get<geometry::min_corner, 0>(AOI);
	double min_y = geometry::get<geometry::min_corner, 1>(AOI);
	double max_x = geometry::get<geometry::max_corner, 0>(AOI);
	double max_y = geometry::get<geometry::max_corner, 1>(AOI);

	point_type ll(min_x, min_y);
	point_type ul(min_x, max_y);
	point_type lr(max_x, min_y);
	
	double newGeoTransform[6] = {min_x, cell_size, 0 , max_y, 0, (cell_size * (-1))};
	
	int x_size = (long int) ceil((geometry::distance(ll, lr) / cell_size));
	int y_size = (long int) ceil((geometry::distance(ll, ul) / cell_size));

	vector<vector<double> > vec(x_size , vector<double> (y_size, 0));  
	
	double cur_y = max_y;
	double cur_x = min_x;
	
// query for distance to fault centre at every grid cell----------------
	cout << "Calculating distances to lineamnet centres for raster with size \n"
		 << vec.size()<< " x " << vec[0].size() << endl;
		 
	progress_display * show_progress =  new boost::progress_display(x_size * y_size);
	for (int i = 0; i < x_size; i++)
	{
		cur_y = max_y;
		for (int j = 0; j < y_size; j++)
		{
			point_type minBox(cur_x, (cur_y - cell_size));
			point_type maxBox((cur_x + cell_size), cur_y );
			box pixel(minBox, maxBox);
			point_type cur_pos((cur_x + cell_size/2), (cur_y - cell_size/2));
			if (!geometry::disjoint(pixel, t_AOI))
			{
				DistTree.query(geometry::index::nearest(cur_pos, 1), back_inserter(result));
				vec[i][j] = geometry::distance(cur_pos, result[0].first);
			}
			else 
				vec[i][j] = -256;
			result.clear();
			cur_y-= cell_size;
			 ++(*show_progress);
		}
		cur_x += cell_size;
	}
//write the raster file---------------------------------------------

    std::string out_name = FGraph::add_prefix_suffix_subdirs(lines.out_path, {geom_subdir}, "centre_distance_map", ".tif");
	georef.WriteRASTER(vec, lines.refWKT, newGeoTransform, lines, out_name);
	cout << " done \n" << endl;
}

void GEOMETRIE::P_Maps(VECTOR lines, float box_size)
{

	GEO georef;
	geometry::index::rtree<p_index, geometry::index::rstar<16>> DistTree;
	vector<p_index> result;

//now we nedd to create a georefernece system based on the bounding box

	box AOI = ReturnAOI(lines.data);
	polygon_type t_AOI = Return_tigth_AOI(lines.data);
	
	double min_x = geometry::get<geometry::min_corner, 0>(AOI);
	double min_y = geometry::get<geometry::min_corner, 1>(AOI);
	double max_x = geometry::get<geometry::max_corner, 0>(AOI);
	double max_y = geometry::get<geometry::max_corner, 1>(AOI);

	point_type ll(min_x, min_y);
	point_type ul(min_x, max_y);
	point_type lr(max_x, min_y);
	
	double newGeoTransform[6] = {min_x, box_size, 0 , max_y, 0, (box_size * (-1))};
	
	int x_size = (long int) ceil((geometry::distance(ll, lr) / box_size));
	int y_size = (long int) ceil((geometry::distance(ll, ul) / box_size));

	vector<vector<double> > vec_count (x_size , vector<double> (y_size, 0));
	vector<vector<double> > vec_length(x_size , vector<double> (y_size, 0));
	
	double cur_y = max_y;
	double cur_x = min_x;
	
//put the segments into an rtree, so we don't need to check each one----
	typedef std::pair<box, decltype(lines.data)::iterator> box_line; //a bounding box around the linestring, and the linestring
	geometry::index::rtree<box_line, geometry::index::rstar<16>> line_tree;
	for (auto it = lines.data.begin(); it < lines.data.end(); it++)
	{
		box fault_bounding_box = boost::geometry::return_envelope<box>(*it);
		line_tree.insert(std::make_pair(fault_bounding_box, it));
	}

// query intesity and density for every grid cell-----------------------
	cout << "Calulating P20 and P21 maps for raster with size \n"
		 << vec_count.size()<< " x " << vec_count[0].size() << endl;
	progress_display * show_progress =  new boost::progress_display(x_size * y_size);

	for (int i = 0; i < x_size ; i++)
	{
		cur_y = max_y;
		for (int j = 0; j < y_size ; j++)
		{
			point_type cur_pos(cur_x, cur_y); //cur_pos is the top-left corner of the pixel
			point_type minBox(cur_x, (cur_y - box_size));
			point_type maxBox((cur_x + box_size), cur_y );
			box pixel(minBox, maxBox);
			if (!geometry::disjoint(pixel, AOI)) //or t_AOI
			{
				double intersec = 0;
				double intersection_length = 0;

	//get the lines that have intersecting bounding boxes-------------------
				std::vector<box_line> candidates;
				line_tree.query(geometry::index::intersects(pixel), std::back_inserter(candidates));
				for (auto candidate = candidates.begin(); candidate < candidates.end(); candidate++)
				{
	//then check the full linestring to see if they intersect with this pixel
					if (!geometry::disjoint(*candidate->second, pixel))
					{
						intersec++; //intersection count for the P20 map
	//and sum the length of the intersection(s) for the P21 map-------------

						geometry::model::multi_linestring<line_type> intersecting;
						geometry::intersection(pixel, *candidate->second, intersecting);
						for (auto intersec_it = intersecting.begin(); intersec_it < intersecting.end(); intersec_it++)
						{
							intersection_length += geometry::length(*intersec_it);
						}
					}
				}
				vec_count[i][j] = intersec;
				vec_length[i][j] = intersection_length;
			}
			else
			{
				vec_count[i][j]  = -256;
				vec_length[i][j] = -256;
			}
			cur_y-= box_size;
			++(*show_progress);
		}
		cur_x += box_size;
	}
//write the raster files---------------------------------------------
    std::string p20_name = FGraph::add_prefix_suffix_subdirs(lines.out_path, {geom_subdir}, "P20_map", ".tif");
	georef.WriteRASTER(vec_count,  lines.refWKT, newGeoTransform, lines, p20_name);
    
    std::string p21_name = FGraph::add_prefix_suffix_subdirs(lines.out_path, {geom_subdir}, "P21_map", ".tif");
	georef.WriteRASTER(vec_length, lines.refWKT, newGeoTransform, lines, p21_name);
	cout << "done \n" << endl;
}
