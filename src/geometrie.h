#ifndef _GEOMETRIE_h
#define _GEOMETRIE_h
#include "main.h"

using namespace boost;
using namespace FGraph;

enum class AttachPoint : signed char {front=-1, middle=0, back=+1};

class GEOMETRIE
{
	public:
	
	GEOMETRIE();                           
	~GEOMETRIE()
	{}
	;  
		
//construct perpendicular line through midpoint of segment--------------
	template <typename Segment>
	struct Perpencicular
	{
		line_type cross;
		point_type centre, p1, p2;
		int seg_count;
		vector<line_type> all;

		Perpencicular() : cross(), centre(), p1(), p2(), all()
		{seg_count = 0;}

		inline void operator()(Segment const& s)
		{
			double rx, ry, l;
			
			geometry::convert(s, cross);
			p1 = cross.front();
			p2 = cross.back();
			geometry::centroid(cross, centre);
				
			rx = p2.x() - p1.x();
			ry = p2.y() - p1.y();
			l = sqrt(rx*rx + ry*ry);
			geometry::clear(p1);
			geometry::clear(p2);

			p1.set<0>((centre.x() + ( ry/l) * 60 ));
			p1.set<1>((centre.y() + (-rx/l) * 60 ));
			p2.set<0>((centre.x() + (-ry/l) * 60 ));
			p2.set<1>((centre.y() + ( rx/l) * 60 ));
			geometry::clear(cross);
			
			geometry::append(cross, p1);
			geometry::append(cross, centre);
			geometry::append(cross, p2);
			all.push_back(cross);
		}
	};
	
	template <typename Segment>
	struct SegPoints
	{
		line_type cross;
		point_type centre, p1, p2;
		int seg_count;
		vector<std::tuple<point_type, point_type, point_type>> Points;
		
		SegPoints() : cross(), centre(), p1(), p2(), Points()
		{seg_count = 0;}

		inline void operator()(Segment const& s)
		{
			geometry::convert(s, cross);
			geometry::centroid(cross, centre);
			Points.push_back(make_tuple(cross.front(), centre, cross.back()));
		}
	};
	

	//find the centre of a segment------------------------------------------ 
	template <typename Segment>
	struct MidPoint
	{
		line_type cross;
		point_type centre;
		vector<point_type> all;

		MidPoint() : cross(), centre(), all()
		{}

		inline void operator()(Segment const& s)
		{
			geometry::convert(s, cross);
			geometry::centroid(cross, centre);
			geometry::clear(cross);
			all.push_back(centre);
		}
	};

	BUFFER DefineSquareBuffer(point_type POINT, const double Bdistance );
	BUFFER DefinePointBuffer(point_type POINT, const double Bdistance );
	BUFFER DefineLineBuffer(line_type fault, const double Bdistance);
	double minSpacing(line_type Trace);
	line_type GetSegment( line_type Trace, point_type Junction, point_type Begin);
	void SortDist(vector <std::tuple<long double, point_type, AttachPoint>>& cross);
	
	void CentreDistanceMap (VECTOR lines, float cell_size);
	void P_Maps(VECTOR lines, float cell_size);
	
	
	box ReturnAOI(vector<line_type> lines);
	line_type ShortestLine(vector <line_type> Set);

};
 #endif
