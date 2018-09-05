#ifndef __SMF_POLYGON_2D___
#define __SMF_POLYGON_2D___

#include "../SMF_Config.h"
#include "../math/all.h"
#include "all.h"

namespace SMF {
namespace MATH{
class CVec2D;

}

using namespace MATH;
using namespace Util;

namespace GEO{

class CAABBox2D;


/**
 * \class CPolygon2D
 *
 * \ingroup SMF_Geometric
 * \image html pics\polyedron.png
 *
 * \if pt_br
 * \brief   Representa um Polígono em espaço cartesiano 2D
 * \elseif us_en
 * \brief 	2D polygon region class
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CPolygon2D {

private:
    std::vector< CPoint2D > M_vertices;

public:
    /**
      \brief create empty polygon
    */
    CPolygon2D();

    /**
      \brief create polygon with points
      \param v array of points
    */
    CPolygon2D( const std::vector< CPoint2D > & v );

   /** \brief  Get the ith point*/
   SMF_INLINE CPoint2D & operator[](int i) { return M_vertices[i]; }


  /** \brief  Get the ith point*/
  SMF_INLINE CPoint2D const& operator[](int i) const { return M_vertices[i]; }



    /**
      \brief clear all data.
    */
    void clear();

    /**
      \brief set polygon with points
      \param v array of points
      \return const reference to itself
    */
    const CPolygon2D & set( const std::vector< CPoint2D > & v );

    /**
      \brief append point to polygon
      \param p new point
    */
    void addVertex( const CPoint2D & p )
      {
          M_vertices.push_back( p );
      }

    /**
      \brief get list of point of this polygon
      \return const reference to point list
    */
    const std::vector< CPoint2D > & vertices() const
      {
          return M_vertices;
      }
    /**
      \brief get list of point of this polygon
      \return const reference to point list
    */
    const int & numVertices() const
      {
          return M_vertices.size();
      }

    /**
      \brief check point is in this polygon or not. the point on segment lines is allowed.
      \param p point for checking
      \return true if point is in this polygon
    */
    virtual
    bool contains( const CPoint2D & p ) const
      {
          return contains( p, true );
      }

    /**
      \brief check point is in this polygon or not
      \param p point for checking
      \param allow_on_segment when point is on outline,
      if this parameter is set to true, returns true
      \return true if point is in this polygon
    */
    bool contains( const CPoint2D & p,
                   const bool allow_on_segment ) const;

    /**
      \brief get bounding box of this polygon
      \return bounding box of this polygon
    */
    CAABBox2D getBoundingBox() const;

    /**
      \brief get centor of bounding box of this polygon
      \return centor of bounding box of this polygon
    */
    CPoint2D xyCenter() const;

	CPoint2D centroid()const;
    /**
      \brief get minimum distance between this polygon and point
      \param p point
      \param closed if this parameter is set to true, handle this
      polygon as a closed polygon,
      otherwise handle this polygon as a polyline or opened polygon.
      when point is inside of this polygon, distance between plane polygon
      and point is 0,
      distance between polyline polygon and point is minimum distance
      between each segments of this polygon.
      \return minimum distance between this polygon and point
    */
    float distance( const CPoint2D & p, const bool closed = true ) const;

    /**
      \brief get area of this polygon
      \return value of area with sign.
    */
    virtual
    float area() const;

    /**
      \brief calculate signed area value
      \return value of signed area.
      If vertices are placed counterclockwise order, returns positive number.
      If vertices are placed clockwise order, returns negative number.
      Otherwise, returns 0.
    */
    float signedArea() const;

    /**
      \brief check vertexes of this polygon is placed counterclockwise ot not
      \return true if counterclockwise
    */
    bool isCounterclockwise() const;

    /**
      \brief check vertexes of this polygon is placed clockwise ot not
      \return true if clockwise
    */
    bool isClockwise() const;
	/**
      \brief get a polygon clipped by a rectangle
      \param r rectangle for clipping
      \return a polygon. if polygon is separated by edges of rectangle,
      each separated polygon is connected to one polygon.
    */
    CPolygon2D getScissoredConnectedPolygon( const CAABBox2D & r ) const;

	/** \brief  Return the point on the given polygon closest to the given point
//  If the third argument is "false", the edge from last to first point of
//  each polygon sheet is not considered part of the polygon.
// \see  CPoint2D
// \see  CPolygon2D*/
static CPoint2D closestPoint(CPolygon2D const& poly, CPoint2D const& point);

CPoint2D closestPoint(CPoint2D const& point)const;
/**
\brief return true if the polygon intersect with the segment
	**/
bool intersects( const CLineSegment2D &Segment)const;
bool intersects( const CLine2D &line)const;
bool intersects( const CRay2D &ray)const;
bool intersects( const CTriangle2D &ray)const;
bool intersects( const CCircle2D &circle)const;

	/** \brief  Pretty print*/
  ostream& toString(ostream&) const;

};
namespace _2D{

SMF_INLINE ostream& operator<< (ostream& os, CPolygon2D const& p) { return p.toString(os); }


}
} //end GEO
} //end SMF
#endif
