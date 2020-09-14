#ifndef _GEOMETRIE_h
#define _GEOMETRIE_h
#include "../include/fracg.h"

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
		vector<line_type> all;
		double d;

		Perpencicular() : cross(), centre(), p1(), p2(), all(), d()
		{}

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

			p1.set<0>((centre.x() + ( ry/l) * d ));
			p1.set<1>((centre.y() + (-rx/l) * d ));
			p2.set<0>((centre.x() + (-ry/l) * d ));
			p2.set<1>((centre.y() + ( rx/l) * d ));
			geometry::clear(cross);
			
			geometry::append(cross, p1);
			geometry::append(cross, centre);
			geometry::append(cross, p2);
			all.push_back(cross);
		}
	};
	
//gets the segmnet intersection with lne-----------------------------
	template <typename Segment>
	struct GetMidSegment
	{
		line_type cross;
		line_type midSeg;
		
		GetMidSegment() : cross(), midSeg()
		{}
		inline void operator()(Segment const& s)
		{
	
		if (geometry::intersects(s, cross))
			geometry::convert(s, midSeg);
		 
		}
	};
	
	

	BUFFER DefineSquareBuffer(point_type POINT, const double Bdistance );
	BUFFER DefinePointBuffer(point_type POINT, const double Bdistance );
	BUFFER DefineLineBuffer(line_type fault, const double Bdistance);
	double MinSpacing(line_type Trace);
	line_type GetSegment( line_type Trace, point_type Junction, point_type Begin);
	void SortDist(vector <std::tuple<long double, point_type, AttachPoint>>& cross);
	
	void CentreDistanceMap (VECTOR lines, float cell_size);
	void PMaps(VECTOR lines, float cell_size);
	
	
	box ReturnAOI(vector<line_type> lines);
	polygon_type ReturnTightAOI(vector<line_type> lines);
	line_type ShortestLine(vector <line_type> Set);

};
 #endif
