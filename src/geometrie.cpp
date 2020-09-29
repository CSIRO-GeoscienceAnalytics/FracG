#include "../include/geometrie.h"
#include "../include/GeoRef.h"
#include "../include/util.h"

const std::string geom_subdir = "geometry";

namespace FracG
{
	namespace bgm = boost::geometry;
	//create a polygon around a point (square)------------------------------
	BUFFER DefineSquareBuffer(point_type POINT, const double Bdistance )
	{
		BUFFER buff;
		bgm::strategy::buffer::point_square point_strategy;
		bgm::strategy::buffer::distance_symmetric<double> distance_strategy(Bdistance);
		bgm::strategy::buffer::join_round join_strategy;
		bgm::strategy::buffer::end_round end_strategy;
		bgm::strategy::buffer::side_straight side_strategy;

		bgm::buffer(POINT, buff,
		distance_strategy, side_strategy,
		join_strategy, end_strategy, point_strategy); 
		return buff; 
	}

	//define a polygon around a point (circle)------------------------------
	BUFFER DefinePointBuffer(point_type POINT, const double Bdistance )
	{
		BUFFER buff;
		const int points_per_circle = 36;

		bgm::strategy::buffer::distance_symmetric<double> distance_strategy( Bdistance );
		bgm::strategy::buffer::join_round join_strategy(points_per_circle);
		bgm::strategy::buffer::end_round end_strategy(points_per_circle);
		bgm::strategy::buffer::point_circle circle_strategy(points_per_circle);
		bgm::strategy::buffer::side_straight side_strategy;

		bgm::buffer(POINT, buff,
		distance_strategy, side_strategy,
		join_strategy, end_strategy, circle_strategy); 
		return buff; 
	}

	//define a polygon around line------------------------------------------
	BUFFER DefineLineBuffer(line_type fault, const double Bdistance)
	{
		BUFFER buff;
		const int points_per_circle = 36;

		bgm::strategy::buffer::distance_symmetric<double> distance_strategy( Bdistance );
		bgm::strategy::buffer::join_round join_strategy(points_per_circle);
		bgm::strategy::buffer::end_round end_strategy(points_per_circle);
		bgm::strategy::buffer::point_circle circle_strategy(points_per_circle);
		bgm::strategy::buffer::side_straight side_strategy;

		boost::geometry::buffer(fault, buff,
						distance_strategy, side_strategy,
						join_strategy, end_strategy, circle_strategy); 
		return buff; 
	}

	//find minimum spacing between point in  a linestring------------------ 
	double MinSpacing(line_type Trace)
	{
		double len = 0;
		point_type first = Trace.front();
		BOOST_FOREACH(point_type P, Trace)
		{
			if (!bgm::equals(first, P) && (len > bgm::distance(first, P)))
				len = bgm::distance(first, P);	
			first = P;
		}
		return len;
	}

	line_type ShortestLine(std::vector <line_type> &Set)
	{
		double dist = std::numeric_limits<double>::max();
		line_type shortest;
		BOOST_FOREACH(line_type l, Set)
		{
			if (bgm::length(l) < dist && bgm::length(l) != 0)
			{
				dist = bgm::length(l) ;
				shortest = l;
			}
		}
		return(shortest);
	}

	box ReturnAOI(std::vector<line_type> &lines)
	{
		box AOI;
		polygon_type multiL;
		BOOST_FOREACH(line_type F, lines)
			bgm::append(multiL, F);

		AOI = boost::geometry::return_envelope<box>(multiL);
		return(AOI);
	}

	polygon_type ReturnTightAOI(std::vector<line_type> &lines)
	{
		polygon_type hull;
		polygon_type multiL;
		BOOST_FOREACH(line_type F, lines)
			bgm::append(multiL, F);

		boost::geometry::convex_hull(multiL, hull);
		return(hull);
	}

	//return the segment of a line that is the closest to the given point
	//the returned iterator points to the linsestring point that is *after* the given point
	line_type::iterator GetClosestSegment(line_type &Trace, point_type point)
	{
		 double best_distance = std::numeric_limits<double>::infinity();
		 line_type::iterator best_index = Trace.begin();
		 point_type prev_point = Trace.front();
		 for (line_type::iterator it = Trace.begin(); it != Trace.end(); it++)
		 {
			 Segment seg(prev_point, *it);
			 double distance = bgm::distance(seg, point);
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
	line_type GetSegment(line_type Trace, point_type Junction, point_type Begin)
	{
		box outL;
		line_type line_seg;
		const double threshold = 1;
		point_type first, last;

		//determine which line segment each point is located on
		line_type::iterator jit = GetClosestSegment(Trace, Junction);
		line_type::iterator bit = GetClosestSegment(Trace, Begin);

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
			bgm::append(line_seg, first);
			bgm::append(line_seg, last);
			return line_seg;
		}
		//the output line segment is the first point -> any and all entire segments that lie between the two points -> from the last whole segment to the end point
		bgm::append(line_seg, first);
		point_type previous = first;
		for (auto it = first_segment; it != last_segment; it++)
		{
			if (bgm::distance(previous, *it) > threshold) 
				bgm::append(line_seg, *it);
			previous = *it;
		}
		bgm::append(line_seg, last);
		return line_seg; 
	}

	//comparison operator that compares a particular index in a (vector of) tuples)
	template<int index, template<typename> class F = std::less>
	struct TupleCompare
	{
		template<typename T>
		bool operator()(T const &t1, T const &t2)
		{
			return F<typename std::tuple_element<index, T>::type>()(std::get<index>(t1), std::get<index>(t2));
		}
	};

	//comparison operator that first checks indexa, then uses indexb if the two indexa values are the same
	template<int indexa, int indexb, template<typename> class F = std::less, template<typename> class G = std::equal_to>
	struct TupleTwoCompare
	{
		template<typename T>
		bool operator()(T const &t1, T const &t2)
		{
			bool first_same = G<typename std::tuple_element<indexa, T>::type>()(std::get<indexa>(t1), std::get<indexa>(t2));
			if (!first_same) return F<typename std::tuple_element<indexa, T>::type>()(std::get<indexa>(t1), std::get<indexa>(t2));
			return F<typename std::tuple_element<indexb, T>::type>()(std::get<indexa>(t1), std::get<indexa>(t2));
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
	void SortDist(std::vector <std::tuple<long double, point_type, AttachPoint>>& cross)
	{
		std::sort(cross.begin(), cross.end(), FaultDistanceCompare());
	}

	void CentreDistanceMap(VECTOR lines, float cell_size)
	{
		point_type centre;
		std::vector<p_index> result;
		bgm::index::rtree<p_index, bgm::index::rstar<16>> DistTree;

		int index = 0 ;

		BOOST_FOREACH(line_type l, lines.data)
		{
			bgm::centroid(l, centre);
			DistTree.insert(std::make_pair(centre, ++index));
		}

	//now we nedd to create a georefernece system based on the bounding box
	//around the features adn the size of the smallest feature--------------

		box AOI = ReturnAOI(lines.data);
		polygon_type t_AOI = ReturnTightAOI(lines.data);

		double min_x = bgm::get<bgm::min_corner, 0>(AOI);
		double min_y = bgm::get<bgm::min_corner, 1>(AOI);
		double max_x = bgm::get<bgm::max_corner, 0>(AOI);
		double max_y = bgm::get<bgm::max_corner, 1>(AOI);

		point_type ll(min_x, min_y);
		point_type ul(min_x, max_y);
		point_type lr(max_x, min_y);

		double newGeoTransform[6] = {min_x, cell_size, 0 , max_y, 0, (cell_size * (-1))};

		int x_size = (long int) ceil((bgm::distance(ll, lr) / cell_size));
		int y_size = (long int) ceil((bgm::distance(ll, ul) / cell_size));

		std::vector<std::vector<double> > vec(x_size , std::vector<double> (y_size, 0));  

		double cur_y = max_y;
		double cur_x = min_x;

	// query for distance to fault centre at every grid cell----------------
		std::cout << "Calculating distances to lineamnet centres for raster with size \n"
			 << vec.size()<< " x " << vec[0].size() << std::endl;

		boost::progress_display * show_progress =  new boost::progress_display(x_size * y_size);
		for (int i = 0; i < x_size; i++)
		{
			cur_y = max_y;
			for (int j = 0; j < y_size; j++)
			{
				point_type cur_pos((cur_x + cell_size/2), (cur_y - cell_size/2));
				if (bgm::within(cur_pos, t_AOI))
				{
					DistTree.query(bgm::index::nearest(cur_pos, 1), std::back_inserter(result));
					vec[i][j] = bgm::distance(cur_pos, result[0].first);
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


		std::string out_name = FracG::AddPrefixSuffixSubdirs(lines.out_path, {geom_subdir}, "centre_distance_map", ".tif");
		WriteRASTER(vec, lines.refWKT, newGeoTransform, lines, out_name);
		std::cout << " done \n" << std::endl;
	}

	void PMaps(VECTOR lines, float box_size)
	{
		bgm::index::rtree<p_index, bgm::index::rstar<16>> DistTree;
		std::vector<p_index> result;

	//now we nedd to create a georefernece system based on the bounding box

		box AOI = ReturnAOI(lines.data);
		polygon_type t_AOI = ReturnTightAOI(lines.data);

		double min_x = bgm::get<bgm::min_corner, 0>(AOI);
		double min_y = bgm::get<bgm::min_corner, 1>(AOI);
		double max_x = bgm::get<bgm::max_corner, 0>(AOI);
		double max_y = bgm::get<bgm::max_corner, 1>(AOI);

		point_type ll(min_x, min_y);
		point_type ul(min_x, max_y);
		point_type lr(max_x, min_y);

		double newGeoTransform[6] = {min_x, box_size, 0 , max_y, 0, (box_size * (-1))};

		int x_size = (long int) ceil((bgm::distance(ll, lr) / box_size));
		int y_size = (long int) ceil((bgm::distance(ll, ul) / box_size));

		std::vector<std::vector<double> > vec_count (x_size , std::vector<double> (y_size, 0));
		std::vector<std::vector<double> > vec_length(x_size , std::vector<double> (y_size, 0));

		double cur_y = max_y;
		double cur_x = min_x;

	//put the segments into an rtree, so we don't need to check each one----
		typedef std::pair<box, decltype(lines.data)::iterator> box_line; //a bounding box around the linestring, and the linestring
		bgm::index::rtree<box_line, bgm::index::rstar<16>> line_tree;
		for (auto it = lines.data.begin(); it < lines.data.end(); it++)
		{
			box fault_bounding_box = boost::geometry::return_envelope<box>(*it);
			line_tree.insert(std::make_pair(fault_bounding_box, it));
		}

	// query intesity and density for every grid cell-----------------------
		std::cout << "Calulating P20 and P21 maps for raster with size \n"
			 << vec_count.size()<< " x " << vec_count[0].size() << std::endl;
		boost::progress_display * show_progress =  new boost::progress_display(x_size * y_size);

		for (int i = 0; i < x_size; i++)
		{
			cur_y = max_y;
			for (int j = 0; j < y_size; j++)
			{
				point_type cur_pos(cur_x, cur_y); //cur_pos is the bottom-left corner of the pixel
				if (bgm::within(cur_pos,t_AOI))
				{
					point_type minBox(cur_x, (cur_y - box_size));
					point_type maxBox((cur_x + box_size), cur_y );

					box pixel(minBox, maxBox);
					double intersec = 0;
					double intersection_length = 0;

		//get the lines that have intersecting bounding boxes-------------------
					std::vector<box_line> candidates;
					line_tree.query(bgm::index::intersects(pixel), std::back_inserter(candidates));
					for (auto candidate = candidates.begin(); candidate < candidates.end(); candidate++)
					{
		//then check the full linestring to see if they intersect with this pixel
						if (!bgm::disjoint(*candidate->second, pixel))
						{
							intersec++; //intersection count for the P20 map
		//and sum the length of the intersection(s) for the P21 map-------------

							bgm::model::multi_linestring<line_type> intersecting;
							bgm::intersection(pixel, *candidate->second, intersecting);
							for (auto intersec_it = intersecting.begin(); intersec_it < intersecting.end(); intersec_it++)
							{
								intersection_length += bgm::length(*intersec_it);
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
		std::string p20_name = FracG::AddPrefixSuffixSubdirs(lines.out_path, {geom_subdir}, "P20_map", ".tif");
		WriteRASTER(vec_count,  lines.refWKT, newGeoTransform, lines, p20_name);

		std::string p21_name = FracG::AddPrefixSuffixSubdirs(lines.out_path, {geom_subdir}, "P21_map", ".tif");
		WriteRASTER(vec_length, lines.refWKT, newGeoTransform, lines, p21_name);
		std::cout << "done \n" << std::endl;
	}
}
