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
 if (!geometry::equals(first, P))
	if (len > geometry::distance(first, P))
	    len = geometry::distance(first, P);	
 first = P;
 }
 return len;	 
 }
 
//split fault trace into segments between intersection------------------ 
 line_type GEOMETRIE::GetSegment(line_type Trace, point_type Begin, point_type Junction)
 {
 bool front = false; 
 box outL;
 line_type line_seg;
 Segment seg(Begin, Junction);
 geometry::envelope(seg, outL);
 
 if (geometry::equals(Trace.back(), Begin))
	front = true;
 
 if (front)
 {
	geometry::append(line_seg, Begin);
		BOOST_FOREACH(point_type P, Trace)
			if (geometry::within(P, outL))
				geometry::append(line_seg, P);
	geometry::append(line_seg, Junction);
 }
 
 else if (!front)
 {
	geometry::append(line_seg, Begin);
		BOOST_REVERSE_FOREACH(point_type P, Trace)
			if (geometry::within(P, outL))
				geometry::append(line_seg, P);
	geometry::append(line_seg, Junction);
 }
 
 
 return line_seg; 
 }
 

//sort distances points of graph edge accoding to distrances to a point-
 void GEOMETRIE::SortDist(vector <std::tuple<long double, point_type, edge_iter>>& cross)
 {
 vector<long double> Distances;
 vector <std::tuple<long double, point_type, edge_iter>> NewCross;
 
 for (typename vector <std::tuple<long double, point_type, edge_iter>>::const_iterator I = cross.begin(); I != cross.end(); I++)
	Distances.push_back(get<0>(*I));
 sort(Distances.begin(), Distances.end()); 

 for (vector<long double>::iterator it = Distances.begin() ; it != Distances.end(); ++it)
	for (typename vector <std::tuple<long double, point_type, edge_iter>>::const_iterator I = cross.begin(); I != cross.end(); I++)
	{
		if (*it == get<0>(*I))
			 NewCross.push_back(*I);
	}
 cross = NewCross;
 }
