#ifndef __SMF_2D_LINE_SEGMENT_
#define __SMF_2D_LINE_SEGMENT_

#include "../SMF_Config.h"
#include "../math/all.h"
#include "SMF_2DPoint.h"
#include "SMF_GeoDefs.h"
namespace SMF {
using namespace MATH;
using namespace Util;
namespace GEO{


class CLine2D;
class CAABBox2D;
/**
 * \class CLineSegment2D
 *
 * \ingroup SMF_Geometric
 *
 * \image html pics\line.png
 * \if pt_br
 * \brief  Representa um segmento de reta defnido por dois pontos no espaço cartesiano
 * \elseif us_en
 * \brief 	Represents a 2D line segment using two points.
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */



class SMF_API CLineSegment2D
{

 public:
   //: One end of line segment
  CPoint2D point1_;
  //: The other end of the line segment
  CPoint2D point2_;
 //: Default constructor - does not initialise!
  SMF_INLINE CLineSegment2D() {}

  //: Copy constructor
  SMF_INLINE CLineSegment2D(CLineSegment2D const& l)
    : point1_(l.point1_), point2_(l.point2_) {}

  //: Construct from two end points
  SMF_INLINE CLineSegment2D(CPoint2D const& p1,
                             CPoint2D const& p2)
    : point1_(p1), point2_(p2) {}

  //: Destructor
  SMF_INLINE ~CLineSegment2D() {}

  //: One end-point of the line segment.
  SMF_INLINE CPoint2D point1() const { return point1_; } // return a copy

  //: The other end-point of the line segment.
  SMF_INLINE CPoint2D point2() const { return point2_; } // return a copy

  //: The equality comparison operator
  SMF_INLINE bool operator==(CLineSegment2D const& l) const {
    return (this==&l) || (point1() == l.point1() && point2() == l.point2())
                      || (point1() == l.point2() && point2() == l.point1()); }

  //: The inequality comparison operator.
  SMF_INLINE bool operator!=(CLineSegment2D const& other)const{return !operator==(other);}

  // A consistent interface with CLine2D:

  //: Parameter a of line a*x + b*y + c = 0
  // a = (y1-y2)
  float a() const;

  //: Parameter b of line a*x + b*y + c = 0
  //b = (x2 - x1)
  float b() const;

  //: Parameter c of line a*x + b*y + c = 0
  //c = (x1*x2)-(y1*y2)
  float c() const;

  CVec2D ToVec2D_point1()const { return CVec2D(point1_.x(),point1_.y());}
  CVec2D ToVec2D_point2()const { return CVec2D(point2_.x(),point2_.y());}
  	  /** \brief Get the ith point
	  \warning if i not [0,1] return Infinite Point
	  **/
   SMF_INLINE CPoint2D operator[](int i) { if (i==0) return point1_;else if(i==1) return point1_; else return CPoint2D(CMath::NAN_FLOAT,CMath::NAN_FLOAT);}

  	  /** \brief Get the ith point
	  \warning if i not [0,1] return Infinite Point
	  **/
   SMF_INLINE CPoint2D const operator[](int i) const { if (i==0) return point1_;else if(i==1) return point1_; else return CPoint2D(CMath::NAN_FLOAT,CMath::NAN_FLOAT);}

     /*!
      \brief get 1st point of segment edge
      \return const reference to the vector object
    */
    const CPoint2D & origin() const
      {
          return point1_;
      }

    /*!
      \brief get 2nd point of segment edge
      \return const reference to the vector object
    */
    const CPoint2D & terminal() const
      {
          return point2_;
      }
  //: unit vector describing line direction
  CVec2D direction() const;

  //: unit vector orthogonal to line on the positive half Space of the segment
  CVec2D normal() const;
  //: unit vector orthogonal to line on the positive half Space of the segment
  CVec2D normalPos() const;
  //: unit vector orthogonal to line on the negative half Space of the segment
  CVec2D normalNeg() const;
  /** 
  \brief return vector orthogonal to line on the positive half Space of the segment
  \note the vector returned is not normalized
  **/
  CVec2D orthogonal() const;
  /** 
  \brief return  vector orthogonal to line on the positive half Space of the segment (Clock Wised)
  \note the vector returned is not normalized
  **/
  CVec2D orthogonalCCW() const;
  /** 
  \brief return  vector orthogonal to line on the  negative half Space of the segment (Counter Clock Wised)
  \note the vector returned is not normalized
  **/
  CVec2D orthogonalCW() const;

  /**
  \return the midpoint of this segment
  \note midpoint = (point1 + point2)/2
  **/
  CPoint2D midPoint()const;
  
  /**
  \bref return a orthogonal line that passes through the given point on the segment
  **/
  CLine2D orthogonal(CPoint2D const &p)const;

  //: angle with the oriented horizontal line y=0, measured in radians.
  //  Returns values between -pi and pi.
  double slope_radians() const;

  //: angle with the oriented horizontal line y=0, measured in 360-degrees.
  //  Returns values between -180 and 180.
  double slope_degrees() const;

  //: Assignment
  SMF_INLINE void set(CPoint2D const& p1, CPoint2D const& p2) {
    point1_ = p1; point2_ = p2; }

  SMF_INLINE CLineSegment2D & operator=(const CLineSegment2D & p)
  { point1_ = p.point1_; point2_ = p.point2_; return *this; }
  SMF_INLINE CLineSegment2D & operator=(CLineSegment2D & p)
  { point1_ = p.point1_; point2_ = p.point2_; return *this; }


  
  SMF_INLINE bool isFinite() { return point1_.isFinite() || point2_.isFinite();}
  //: Return a point on the line defined by a scalar parameter \a dist.
  // \a dist=0.0 corresponds to point1 and \a dist=1.0 to point2.
  // 0<dist<1 for points on the segment between point1 and point2.
  // dist<0 for points on the (infinite) line, outside the segment, and closer to point1 than to point2.
  // dist>1 for points on the (infinite) line, outside the segment, and closer to point2 than to point1.
  SMF_INLINE CPoint2D getPoint(const sf_u16 dist) const { return point1() + dist*(point2_-point1_); }

  /*!
      \brief check if the point is within the rectangle defined by this
      segment as a diagonal line.
      \return true if rectangle contains p
     */
    bool recContains( const CPoint2D & p ) const
      {
          return ( ( p.x() - origin().x() ) * ( p.x() - terminal().x() ) <= CMath::FLT_EPSILON 
                   && ( p.y() - origin().y() ) * ( p.y() - terminal().y() ) <= CMath::FLT_EPSILON ); 
      }

  /** \brief true if the point lies on the line segment and is between the endpoints
  \see http://lucidarme.me/?p=1952
  **/
  bool intersects(CPoint2D const &p)const;
  /**
  \brief check if this Line Segment intersects the other line
  **/
  bool intersects(CLine2D const &p)const;
  /**
  \brief check if this Line Segment intersects the other line  segment
  **/
  bool intersects(CLineSegment2D const &p)const;
  /**
  \brief check if this Line Segment intersects the other ray
  **/
  bool intersects(CRay2D const &p)const;


 /**
 * \brief Check if a point is inside the current segment
 * \param point coordinates of the point to test
 * \return  NOT_INTERESECTING if the point doesn't lay with the segment
 *          INTERESECTING_EXTREMITY_P1 if the point is merged with P1
 *          INTERESECTING_EXTREMITY_P2 if the point is merged with P2
 *          INTERESECTING if the point belongs to the segment (extremity no included)
 */
  IntersectResult isPointOnSegment_m1(CPoint2D const &p)const;
 /**
 * \brief Check if a point is inside the current segment
 * \param point coordinates of the point to test
 * \return  True if the point is on segment 
 **/
  bool isPointOnSegment_m2(CPoint2D const& p)const;
  
  bool isPointOnSegment(CPoint2D const& p, int method=0)const;

/**
\brief determining the intersection point of two lines (or line segments) in 2 dimensions. 
\return if there are any kind of interection
\param other_line the line to verify the intersection
\param [out] intersection will receive the point of intersection if it exists
\param [out] intType if the segments are PARALLEL, COINCIDENT, NOT_INTERESECTING, INTERESECTING
\see intersection point of two line segments in 2 dimensions at http://paulbourke.net/geometry/pointlineplane/
\note The algorithm on the page above has some problems: 
1 - the lack of classification on coincident lines
2 - when ua=0 the point of intersection is miscalculated. 
I solved this problems in my 
\note If the IntersectResult is NOT_INTERSECTING the point returned is the intersection point if the segments were extended until they had intersection
**/
static	bool intersects_m1(const CLineSegment2D& this_segment, const CLineSegment2D& other_segment, CPoint2D& intersection, int *intType=NULL, float EPS= 1e-3f);
/**
\brief retur true if there is an intersection
\note FAST -> just check for intersection. Do not calc intersection point
**/
static bool simpleIntersects(const CLineSegment2D &Segment1,const CLineSegment2D &Segment2);
/**
\brief intersection between segments seg1(Point1,Point2) and seg2(Point3,Point4)
**/
static bool simpleIntersects(const CPoint2D &Point1,const CPoint2D &Point2,const CPoint2D &Point3,const CPoint2D &Point4);
/**
\brief determining the intersection point of two lines (or line segments) in 2 dimensions. 
\param other_line the line to verify the intersection
\param [out] intersection will receive the point of intersection if it exists
\param [out] intType  intersects Result Bitmask
\return true if there are any kind of intersection
\note if intersection point return (NAN_FLOAT,NANFLOAT), it means it is not possible to calculate one point of intersection, or the segments are parallel or coincidents.
        if they are coincidents, there are a lot of intersection points.
\note Uses the algorith from Olivier renault on Pollycolly code

\note If the IntersectResult is NOT_INTERSECTING the point returned is the intersection point if the segments were extended until they had intersection
**/
static bool  intersects_m2(const CLineSegment2D& this_segment, const CLineSegment2D& other_segment, CPoint2D& intersection, int *intType=NULL,float EPS= 1e-3f);
static bool  intersects_m2(const CLineSegment2D& this_segment, const CLine2D& other_line, CPoint2D& intersection, int *intType=NULL,float EPS= 1e-3f);
static bool  intersects_m2(const CLineSegment2D& this_segment, const CRay2D& other_ray, CPoint2D& intersection, int *intType=NULL,float EPS= 1e-3f);
/**
\return true if there are intersection
**/
bool intersects(const CLineSegment2D& other_segment, CPoint2D& intersection);

static bool intersects(const CLineSegment2D& this_segment, const CLineSegment2D& other_segment, CPoint2D& intersection, float EPS= 1e-3f);

/** Returns true if there are intersection between this segment and the CABBox
\param ICnt [out] the number of intersections of a line segment with a box, up to two are returned in p0 and p1. 
**/
bool intersection(const CAABBox2D& box,int *ICnt,  CPoint2D* p00, CPoint2D* p01)const;

/**
  \brief return the direction vector of this line
  \see http://www.netcomuk.co.uk/~jenolive/vect3.html 
  **/
  CVec2D directionVector()const;

	/// Projects this CLineSegment2D onto the given 1D axis od carthesian.
	CLineSegment2D projectToXAxis() const;
	CLineSegment2D projectToYAxis() const;

/**  
	\brief get minimum distance from a point \a p to a line segment \a l in 2D
      \param p point
      \return minimum distance between this segment and point
    **/
float distance(CPoint2D const& p) const;
float distanceSq(CPoint2D const& p) const;

static float distanceSegmentPoint(CLineSegment2D const &l1, CPoint2D const &p);

// Given this segment  and point p, computes closest point d on the segment.
// Also returns t for the parametric position of d, d(t)=a+t*(b - a)
CPoint2D closestPoint(CPoint2D const  &p, float *t=NULL)const;
/// Computes the closest point on this segment to the given object
//CPoint2D closestPoint(const CRay2D &other, float *d = 0, float *d2 = 0) const;
//CPoint2D closestPoint(const CLine2D &other, float *d = 0, float *d2 = 0) const;
CPoint2D closestPoint(const CLineSegment2D &other) const;

/**
\brief convert this segmento to a CLine2D
**/
CLine2D toLine2D()const;
/**
\brief convert this segmento to a CABBox2D
**/
CAABBox2D toCAABBox()const;
};

namespace _2D{
//: Write to stream
// \see  CLineSegment2D

ostream&  operator<<(ostream& s, const CLineSegment2D & p);

//: Read from stream
// \see  CLineSegment2D

istream&  operator>>(istream& is,  CLineSegment2D & p);

}//end namespace _2D

} //end GEO
}  //end SMF

#endif // __SMF_2D_POINT_
