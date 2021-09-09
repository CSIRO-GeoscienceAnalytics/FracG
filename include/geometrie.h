/****************************************************************/
/*				DO NOT MODIFY THIS HEADER							*/
/*					FRACG - FRACture Graph							*/
/*				Network analysis and meshing software					*/
/*																		*/
/*						(c) 2021 CSIRO									*/
/*			GNU General Public Licence version 3 (GPLv3)				*/
/*																		*/
/*						Prepared by CSIRO								*/
/*																		*/
/*					See license for full restrictions 						*/
/****************************************************************/
#ifndef _GEOMETRIE_h
#define _GEOMETRIE_h
#include "../include/fracg.h"

namespace FracG
{

    enum class AttachPoint : signed char {front=-1, middle=0, back=+1};

    //construct perpendicular line through midpoint of segment--------------
    template <typename Segment>
    struct Perpencicular
    {
            line_type cross;
            point_type centre, p1, p2;
            std::vector<line_type> all;
            double d;

            Perpencicular() : cross(), centre(), p1(), p2(), all(), d()
            {}

            inline void operator()(Segment const& s)
            {
                    double rx, ry, l;

                    boost::geometry::convert(s, cross);
                    p1 = cross.front();
                    p2 = cross.back();
                    
                    boost::geometry::centroid(cross, centre);

                    rx = p2.x() - p1.x();
                    ry = p2.y() - p1.y();
                    l = sqrt(rx*rx + ry*ry);
                    boost::geometry::clear(p1);
                    boost::geometry::clear(p2);

                    p1.set<0>((centre.x() + ( ry/l) * d ));
                    p1.set<1>((centre.y() + (-rx/l) * d ));
                    p2.set<0>((centre.x() + (-ry/l) * d ));
                    p2.set<1>((centre.y() + ( rx/l) * d ));
                    boost::geometry::clear(cross);

                    boost::geometry::append(cross, p1);
                    boost::geometry::append(cross, centre);
                    boost::geometry::append(cross, p2);
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

            if (boost::geometry::intersects(s, cross))
                    boost::geometry::convert(s, midSeg);

            }
    };

    BUFFER DefineSquareBuffer(point_type POINT, const double Bdistance );
    BUFFER DefinePointBuffer(point_type POINT, const double Bdistance );
    BUFFER DefineLineBuffer(line_type fault, const double Bdistance);
    double MinSpacing(line_type Trace);
    line_type GetSegment( line_type Trace, point_type Junction, point_type Begin);
    void SortDist(std::vector<std::tuple<long double, point_type, AttachPoint>>& cross);

    void D_Maps (VECTOR lines, float cell_size, bool resample);
    void P_Maps(VECTOR lines, float cell_size, bool resample);

	point_type GetCentre(line_type l);

    box ReturnAOI(std::vector<line_type> &lines);
    polygon_type ReturnTightAOI(std::vector<line_type> &lines);
    line_type ShortestLine(std::vector <line_type> &Set);
}
#endif
