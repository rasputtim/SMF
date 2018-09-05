/*
  SMF -  Super Math Fabric  
  Copyright (C) 2014 Rasputtim <Rasputtim@hotmail.com>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the freeMem Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the freeMem Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

  */




#include "geometry/SMF_2Dpolygon.h"
#include "math/all.h"
#include "geometry/all.h"

//#include <cmath>
//#include <cfloat>
//#include <cstdlib>
//#include <iostream>

namespace SMF {
using namespace MATH;
using namespace Util;
namespace GEO{

/*-------------------------------------------------------------------*/
/*!

 */
CPolygon2D::CPolygon2D()
    : M_vertices()
{

}

/*-------------------------------------------------------------------*/
/*!

 */
CPolygon2D::CPolygon2D( const std::vector< CPoint2D > & v )
    : M_vertices( v )
{

}

/*-------------------------------------------------------------------*/
/*!

 */
void CPolygon2D::clear()
{
    M_vertices.clear();
}

/*-------------------------------------------------------------------*/
/*!

 */
const CPolygon2D & CPolygon2D::set( const std::vector< CPoint2D > & v )
{
    M_vertices = v;
    return *this;
}

/*-------------------------------------------------------------------*/
/*!

 */
CAABBox2D CPolygon2D::getBoundingBox() const
{
    if ( M_vertices.empty() )
    {
        return CAABBox2D();
    }

    float x_min = +IEEE_FLT_MAX;
    float x_max = -IEEE_FLT_MAX;
    float y_min = +IEEE_FLT_MAX;
    float y_max = -IEEE_FLT_MAX;

    const std::vector< CPoint2D >::const_iterator end = M_vertices.end();
    for ( std::vector< CPoint2D >::const_iterator p = M_vertices.begin();
          p != end;
          ++p )
    {
        if ( p->x() > x_max )
        {
            x_max = p->x();
        }

        if ( p->x() < x_min )
        {
            x_min = p->x();
        }

        if ( p->y() > y_max )
        {
            y_max = p->y();
        }

        if ( p->y() < y_min )
        {
            y_min = p->y();
        }
    }

    return CAABBox2D( x_min, y_min , x_max , y_max );
}

/*-------------------------------------------------------------------*/
/*!

 */
CPoint2D CPolygon2D::xyCenter() const
{
    return this -> getBoundingBox().centroid();
}
CPoint2D CPolygon2D::centroid()const
{
	const CPolygon2D & Polygon=*this;
	CPoint2D point;
	point.toInfinite();
  int i,j,n;
  float asum ,term;

  float x = 0;
  float y = 0;
  n= numVertices();
  if (n < 3) return point;

  asum = CMath::ZERO;
  j    = n - 1;

  for (i = 0; i < n;i++) {

    term = ((Polygon[j].x_ * Polygon[i].y_) - (Polygon[j].y_ * Polygon[i].x_));
    asum = asum + term;
    x = x + (Polygon[j].x_ + Polygon[i].x_) * term;
    y = y + (Polygon[j].y_ + Polygon[i].y_) * term;
    j = i;
  }

  if (!CMath::equals(asum, CMath::ZERO)) {
    x = x / (3.0 * asum);
    y = y / (3.0 * asum);
  }
  point.set(x,y);
  return point;
}

/*-------------------------------------------------------------------*/
/*!

 */
bool
CPolygon2D::contains( const CPoint2D & p, const bool allow_on_segment ) const
{
    if ( M_vertices.empty() )
    {
        return false;
    }
    else if ( M_vertices.size() == 1 )
    {
        return allow_on_segment
            //&& ( M_vertices[0].equals( p ) );
            && ( M_vertices[0] == p );

    }


    CAABBox2D r = this -> getBoundingBox();

    if ( ! r.contains( p ) )
    {
        return false;
    }


    //
    // make virtual half line
    //
    CLineSegment2D line( p, CVec2D( p.x() + ((r.maxX() - r.minX() + r.maxY() - r.minY())
                                        + (M_vertices[0] - p).getLenght()) * 3.0,
                                 p.y() ) );

    //
    // check intersection with all segments
    //
    bool inside = false;
    float min_line_x = r.maxX() + 1.0;

    for ( size_t i = 0; i < M_vertices.size(); ++i )
    {
        size_t p1_index = i + 1;

        if ( p1_index >= M_vertices.size() )
        {
            p1_index = 0;
        }

        const CPoint2D p0 = M_vertices[i];
        const CPoint2D p1 = M_vertices[p1_index];

        if ( ! allow_on_segment )
        {
            if ( CLineSegment2D( p0, p1 ).isPointOnSegment( p ) )
            {
                return false;
            }
        }

        if ( allow_on_segment
            //&& p.equalsStrictly( p0 ) )
            && p == p0 )

        {
            return true;
        }

        if ( line.intersects( CLineSegment2D( p0, p1 ) ) )
        {
            if ( p0.y() == p.y()
                 || p1.y() == p.y() )
            {
                if ( p0.y() == p.y() )
                {
                    if ( p0.x() < min_line_x )
                    {
                        min_line_x = p0.x();
                    }
                }

                if ( p1.y() == p.y() )
                {
                    if ( p1.x() < min_line_x )
                    {
                        min_line_x = p1.x();
                    }
                }


                if ( p0.y() == p1.y() )
                {
                    continue;
                }
                else if ( p0.y() < p.y()
                          || p1.y() < p.y() )
                {
                    continue;
                }
                else // bottom point on the line
                {
                    // no operation, don't skip
                }
            }

            inside = (! inside);
        }
    }

    return inside;
}


/*-------------------------------------------------------------------*/
/*!

 */
float
CPolygon2D::distance( const CPoint2D & p, bool closed ) const
{
    const size_t size = vertices().size();

    if ( size == 1 )
    {
        return ( M_vertices[0] - p ).getLenght();
    }

    if ( closed && contains( p ) )
    {
        return 0.0;
    }

    float min_dist = +IEEE_FLT_MAX;

    for ( size_t i = 0; i + 1 < size; ++i )
    {
        CLineSegment2D seg( M_vertices[i], M_vertices[i + 1] );

        float d = seg.distance( p );

        if ( d < min_dist )
        {
            min_dist = d;
        }
    }

    if ( size >= 3 )
    {
        CLineSegment2D seg( M_vertices.back(), M_vertices.front() );

        float d = seg.distance( p );

        if ( d < min_dist )
        {
            min_dist = d;
        }
    }

    // if this -> vertex().size() == 0, returns huge value

    return min_dist;
}

/*-------------------------------------------------------------------*/
/*!

 */
float CPolygon2D::area() const
{
    return CMath::fabs( signedArea() * 0.5 );
}

/*-------------------------------------------------------------------*/

float CPolygon2D::signedArea() const
{
    const size_t size = M_vertices.size();
    if ( size < 3 )
    {
        return 0.0;
    }

    float sum = 0.0;
    for ( size_t i = 0; i < size; ++i )
    {
        size_t n = i + 1;
        if ( n == size )
        {
            n = 0;
        }

        sum += (M_vertices[i].x() * M_vertices[n].y() - M_vertices[n].x() * M_vertices[i].y());
    }

    return sum;
}

/*-------------------------------------------------------------------*/
/*!

 */
bool
CPolygon2D::isCounterclockwise() const
{
    return signedArea() > 0.0;
}

/*-------------------------------------------------------------------*/
/*!

 */
bool
CPolygon2D::isClockwise() const
{
    return signedArea() < 0.0;
}


/*-------------------------------------------------------------------*/
/*!
  \brief scissorring implementation
*/
template< class Predicate >
void scissorWithLine( const Predicate & in_region,
                 const std::vector< CPoint2D > & points,
                 std::vector< CPoint2D > * new_points,
                 const CLine2D & line )
{
    new_points -> clear();

    std::vector< bool > in_rectangle( points.size() );

    for ( size_t i = 0; i < points.size(); ++i )
    {
        in_rectangle[i] = in_region( points[i] );
    }

    for ( size_t i = 0; i < points.size(); ++i )
    {
        size_t index_0 = i;
        size_t index_1 = i + 1;

        if ( index_1 >= points.size() )
        {
            index_1 = 0;
        }

        const CPoint2D & p0 = points[index_0];
        const CPoint2D & p1 = points[index_1];

        if ( in_rectangle[index_0] )
        {
            if ( in_rectangle[index_1] )
            {
                new_points -> push_back( p1 );
            }
            else
            {
                CPoint2D inter;
				const CLine2D l3(p0,p1);
				CLine2D::intersects_m2(line, l3, inter);
				CVec2D c=inter.toVec2DPos();

				SMF_ASSERT(c.isFinite())

                new_points -> push_back( c );
            }
        }
        else
        {
            if ( in_rectangle[index_1] )
            {
                CPoint2D inter;
				const CLine2D l2( p0, p1 );
				CLine2D::intersects_m2(line, l2 , inter);
				CVec2D c=inter.toVec2DPos();

				SMF_ASSERT(c.isFinite())

                new_points -> push_back( c );
                new_points -> push_back( p1 );
            }
            else
            {
                // noting to do
            }
        }
    }
}

class XLessEqual
{
private:
    float threshold;

public:
    XLessEqual( float threshold )
        : threshold( threshold )
      {
      }

    bool operator()( const CPoint2D & p ) const
      {
          return p.x() <= threshold;
      }
};

class XMoreEqual
{
private:
    float threshold;

public:
    XMoreEqual( float threshold )
        : threshold( threshold )
      {
      }

    bool operator()( const CPoint2D & p ) const
      {
          return p.x() >= threshold;
      }
};

class YLessEqual
{
private:
    float threshold;

public:
    YLessEqual( float threshold )
        : threshold( threshold )
      {
      }

    bool operator()( const CPoint2D & p ) const
      {
          return p.y() <= threshold;
      }
};

class YMoreEqual
{
private:
    float threshold;

public:
    YMoreEqual( float threshold )
        : threshold( threshold )
      {
      }

    bool operator()( const CPoint2D & p ) const
      {
          return p.y() >= threshold;
      }
};

CPolygon2D CPolygon2D::getScissoredConnectedPolygon( const CAABBox2D & r ) const
{
    if ( M_vertices.empty() )
    {
        return CPolygon2D();
    }

    std::vector< CPoint2D > p = M_vertices;
    std::vector< CPoint2D > clipped_p_1;
    std::vector< CPoint2D > clipped_p_2;
    std::vector< CPoint2D > clipped_p_3;
    std::vector< CPoint2D > clipped_p_4;

    scissorWithLine< XLessEqual >( XLessEqual( r.maxX() ),
                                   p, &clipped_p_1,
                                   CLine2D( CPoint2D( r.maxX(), 0.0 ), CVec2D(0,1) ) );

    scissorWithLine< YLessEqual >( YLessEqual( r.maxY() ),
                                   clipped_p_1, &clipped_p_2,
                                   CLine2D( CPoint2D( 0.0, r.maxY() ), CVec2D(1,0)  ) );

    scissorWithLine< XMoreEqual >( XMoreEqual( r.minX() ),
                                   clipped_p_2, &clipped_p_3,
                                   CLine2D( CPoint2D( r.minX(), 0.0 ), CVec2D(0,1)  ) );

    scissorWithLine< YMoreEqual >( YMoreEqual( r.minY() ),
                                   clipped_p_3, &clipped_p_4,
                                   CLine2D( CPoint2D( 0.0, r.minY() ), CVec2D(1,0) ) );

    return CPolygon2D( clipped_p_4 );
}



CPoint2D CPolygon2D::closestPoint(CPoint2D const& point)const{
	return closestPoint(*this,point);
}

CPoint2D CPolygon2D::closestPoint(CPolygon2D const& poly, CPoint2D const& point)
{
  bool closed=true;
  float x=point.x();
  float y=point.y();
  float dd = CPoint2D::distanceTolinesegmentSqr(poly[0].x(),poly[0].y(), poly[1].x(),poly[1].y(), x,y);
  int di = 0;
    unsigned int n = (unsigned int)(poly.numVertices());
    SMF_ASSERT( n > 1 );
    for (unsigned i=0; i+1<n; ++i)
    {
      double nd = CPoint2D::distanceTolinesegmentSqr(poly[i].x(),poly[i].y(), poly[i+1].x(),poly[i+1].y(), x,y);
      if (nd<dd) { dd=nd; di=i; }
    }
    if (closed)
    {  //distancia entre o último vértice e o primeiro
      double nd = CPoint2D::distanceTolinesegmentSqr(poly[0].x(),poly[0].y(), poly[n-1].x(),poly[n-1].y(), x,y);
      if (nd<dd) { dd=nd; di=-1;  }
    }


  float ret_x, ret_y;

  if (di == -1) {// último segmento tem a menor distância
	//CLineSegment2D seg1(CPoint2D(poly[0].x(),poly[0].y()), CPoint2D(poly[n-1].x(),poly[n-1].y()));
	//return   point.closestPoint(seg1);
    CPoint2D::closestPointToLinesegment(ret_x,ret_y, poly[0].x(),poly[0].y(), poly[n-1].x(),poly[n-1].y(), x,y);
}else{
    //CLineSegment2D seg1(CPoint2D(poly[di].x(),poly[di].y()), CPoint2D(poly[di+1].x(),poly[di+1].y()));
	//return   point.closestPoint(seg1);
	CPoint2D::closestPointToLinesegment(ret_x,ret_y, poly[di].x(),poly[di].y(), poly[di+1].x(),poly[di+1].y(), x,y);

  }
  return CPoint2D((ret_x), (ret_y));
}


bool CPolygon2D::intersects( const CLineSegment2D &Segment)const{
const CPolygon2D &Polygon=*this;
  int i,j;
  bool Result = false;
  int n = numVertices();
  if (n < 3 ) return false;
  j = n - 1;
  for (i = 0; i < n;i++ )
  {

    if (CLineSegment2D::simpleIntersects(Segment[0],Segment[1],Polygon[i],Polygon[j])){
      Result = true;
      break;
	}
    j = i;
  }
}
bool CPolygon2D::intersects( const CLine2D &line)const{
const CPolygon2D &Polygon=*this;
  int i,j;
  bool Result = false;
  int n = numVertices();
  if (n < 3 ) return false;
  j = n - 1;
  for (i = 0; i < n;i++ )
  {
    CLineSegment2D segteste(Polygon[i],Polygon[j]);
    if (line.intersects(segteste)){
      Result = true;
      break;
	}
    j = i;
  }
}

bool CPolygon2D::intersects( const CRay2D &ray)const{
const CPolygon2D &Polygon=*this;
  int i,j;
  bool Result = false;
  int n = numVertices();
  if (n < 3 ) return false;
  j = n - 1;
  for (i = 0; i < n;i++ )
  {
    CLineSegment2D segteste(Polygon[i],Polygon[j]);
    if (ray.intersects(segteste)){
      Result = true;
      break;
	}
    j = i;
  }
}
bool CPolygon2D::intersects( const CTriangle2D &triangle)const{
const CPolygon2D &Polygon=*this;
  int i,j;
  bool Result = false;
  int n = numVertices();
  if (n < 3 ) return false;
  j = n - 1;
  for (i = 0; i < n;i++ )
  {
    CLineSegment2D segteste(Polygon[i],Polygon[j]);
    if (triangle.intersects(segteste)){
      Result = true;
      break;
	}
    j = i;
  }
}
bool CPolygon2D::intersects( const CCircle2D &circle)const{
const CPolygon2D &Polygon=*this;
  int i,j;
  bool Result = false;
  int n = numVertices();
  if (n < 3 ) return false;
  j = n - 1;
  for (i = 0; i < n;i++ )
  {
    CLineSegment2D segteste(Polygon[i],Polygon[j]);
    if (circle.intersects(segteste)){
      Result = true;
      break;
	}
    j = i;
  }
}

ostream& CPolygon2D::toString (ostream& os) const
{
  if (M_vertices.size() == 0)
    os << "empty polygon\n";
  else {
    os << "Polygon with " << numVertices() << " vertices:\n";
        for (unsigned int p = 0; p < M_vertices.size(); ++p) {
        os << "Vertice " << p << ' ';
		os << '(' << M_vertices[p].x() << ',' << M_vertices[p].y() << ") ";
		os << '\n';
	   }
  }
  return os;
}


} // END GEO
} //END SMF
