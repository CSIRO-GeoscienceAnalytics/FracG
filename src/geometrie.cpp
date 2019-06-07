#include "geometrie.h"
#include "GeoRef.h"

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
	for (auto it = first_segment; it != last_segment; it++)
	{
		if (geometry::distance(*std::next(it, -1), *it) > threshold) 
			geometry::append(line_seg, *it);
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
 
//sort a vector of points according to their distances along a line
void GEOMETRIE::SortDist(vector <std::tuple<long double, point_type, edge_iter>>& cross)
{
	std::sort(cross.begin(), cross.end(), TupleCompare<0>());
}


