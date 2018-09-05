#ifndef __SMF_2D_LINE_
#define __SMF_2D_LINE_

#include "../SMF_Config.h"
#include "../math/all.h"
#include "SMF_GeoDefs.h"
namespace SMF {
using namespace MATH;
using namespace Util;
namespace GEO{

class CPoint2D;
class CLine2D;
class CAABBox2D;
/**
 * \class CLine2D
 *
 * \ingroup SMF_Geometric
 * \image html pics\line.png
 * \if pt_br
 * \brief Modelo de Reta em espaço cartesiano
 * \note \n fórmula:
 *       \n \f$ a*x + b*y + c = 0 \f$
 *       
 * \elseif us_en
 * \brief 	Represents a Euclidean 2D line
 * An interface for the line coefficients, [a,b,c], is provided in terms of the
 * standard implicit line equation: a*x + b*y + c = 0
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */

class SMF_API CLine2D
{
  // the data associated with this point
  float a_;
  float b_;
  float c_;

 public:
  //: Default constructor (Line 1.y==0, the X axis)
  SMF_INLINE CLine2D() : a_(0), b_(1), c_(0) {}

  //: Construct a CLine2D from its equation, three Types.
  //  The values of a and b should not be both zero.
  SMF_INLINE CLine2D(float ta, float tb, float tc) :a_(ta),b_(tb),c_(tc) { SMF_ASSERT(ta||tb); }

  //: Construct from its equation, a 3-vector.
  //  The values v[0] and v[1] should not be both zero.
  SMF_INLINE CLine2D(const float v[3]):a_(v[0]),b_(v[1]),c_(v[2]) { SMF_ASSERT(a_||b_); }

  //: Construct from homogeneous description of line
  //  The line l should not be the line at infinity.
  //CLine2D (CPoint2DHom const& l);

  //: Construct from two distinct points (join)
  //  The two points must be distinct!
  CLine2D (CPoint2D const& p1, CPoint2D const& p2);

  /** Construct from one point and one direction vector
  \param v the direction vector for this line.
  \param p one point on this line
  **/
  CLine2D (CPoint2D const& p, CVec2D const& v);

  /**
  \brief return the direction vector of this line
  \see http://www.netcomuk.co.uk/~jenolive/vect3.html 
  **/
  CVec2D directionVector()const;

#if 0 // use compiler defaults for these
  // Default destructor
  SMF_INLINE ~CLine2D () {}

  // Default assignment operator
  SMF_INLINE CLine2D& operator=(const CLine2D& l)
  { set(l.a(),l.b(),l.c()); return *this; }
#endif

  //: the comparison operator
  SMF_INLINE bool operator==(CLine2D const& l) const
  {
    return (this==&l) ||
           (a()*l.c()==c()*l.a() && b()*l.c()==c()*l.b() && b()*l.a()==a()*l.b());
  }

  SMF_INLINE bool operator!=(CLine2D const& other) const { return !operator==(other); }

  //: angle with the horizontal line y=0, measured in radians.
  // \brief coeficiente angular
  //  Returns values between -pi and pi, i.e., the lines x-y=0 and y-x=0
  //  return different values (pi/4 and -3pi/4 respectively) although these
  //  lines are identical.
  double slope_radians() const;

  //: angle with the horizontal line y=0, measured in 360-degrees.
  //  Returns values between -180 and 180, i.e., the lines x-y=0 and y-x=0
  //  return different values (45 and -135 respectively) although these
  //  lines are identical.
  double slope_degrees() const;

  // return the angular coeficient in radians
  //caso não haja coeficiente angular, retorna Nan
  float angular_coef() const;
  //reurn if ths line is parallel to the other one
  bool isParallel(CLine2D const &l, float maxDistance = 1e-6f)const;
//reurn if ths line is concurrent to the other one
  bool isConcurrent(CLine2D const &l)const {return !this->isParallel(l);}
  //return if this line is perpendicular to the other one
  bool isPerpendicular(CLine2D const &l,  float maxDistance = 1e-6f)const;

  SMF_INLINE CLine2D orthogonal(CPoint2D const &p)const{ return CLine2D( b(), -a(), (a() * p.y_ - b() * p.x_) ); };
  SMF_INLINE CLine2D perpendicuar(CPoint2D const &p)const{ return orthogonal(p); };

  /**
  \brief return true if the line is vertical on the carthesian plane
  **/
  bool isVertical()const { return !MATH::isFinite(angular_coef());}
  /**
  \brief return true if the line is horizontal on the carthesian plane
  **/
  bool isHorizontal()const { return CMath::nearZero(angular_coef());}

  // Data Access-------------------------------------------------------------

  //: Parameter a of line a*x + b*y + c = 0
  SMF_INLINE float a() const {return a_;}
  //: Parameter b of line a*x + b*y + c = 0
  SMF_INLINE float b() const {return b_;}
  //: Parameter c of line a*x + b*y + c = 0
  SMF_INLINE float c() const {return c_;}

  //: unit vector describing line direction
  SMF_INLINE CVec2D direction() const
  { CVec2D vec(b_,-a_); vec.toNormal(); return vec; }

  //: unit vector orthogonal to line
  SMF_INLINE CVec2D normal() const
  { CVec2D vec(a_,b_); vec.toNormal(); return vec; }

  //: normalize the line coefficients s.t. a^2 + b^2 = 1
  bool normalize();

  //: set a b c.
  //  The values of a and b should not be both zero.
  //  Note that it does not make sense to set a, b or c separately
  SMF_INLINE void set(float ta, float tb, float tc) { SMF_ASSERT(ta||tb); a_=ta; b_=tb; c_=tc; }

  //: Return true iff this line is the line at infinity
  //  This always returns "false"
  SMF_INLINE bool ideal(float = (float)0) const { return false; }

  
	/// Converts this CLine2D to a CLineSegment2D.
	/** \param d Specifies the position of the other endpoint along this CLine. This parameter may be negative.
		\return A CLineSegment2D with point a at pos, and point b at pos + d * dir.
		\see pos, dir, CLine2D::CLine2d, class CLineSegment2D. */
	CLineSegment2D toLineSegment2D(float d) const;
	
	/** Converts this CLine2D to a CLineSegment2D with arbitrary start ena end points.
	 	\return A CLineSegment2D .
		\see pos, dir, CLine2D::CLine2d, class CLineSegment2D. */
	CLineSegment2D toLineSegment2D() const;


  //: Get two points on the line; normally the intersection with X and Y axes.
  // When the line is parallel to one of these,
  // the point with \a y=1 or \a x=1, resp. are taken.  When the line goes
  // through the origin, the second point is (b, -a).
  void get_two_points(CPoint2D &p1, CPoint2D &p2) const;

  /**
  \brief return the point on the line with corresponding horixontal axis = x
  \param x the horizontal axis param for the point to be returned
  \note case the line is vertical return y=NAN_FLOAT, so you can ut any y you want. 
  So, don't forget to test y before use this point
  **/
  CPoint2D getPoint(float x)const;
  /**
  \brief calcula o valor de y, dado x;
  \return y caso seja possivel calculá-lo. caso seja uma reta perpendicular, retorn NAN_FLOAT
  **/
  float getY(float const &x)const;
	/**
      \brief get X-coordinate correspond to 'y'
      \param y considered Y
      \return X coordinate
    **/
    float getX( const double & y ) const;

  //: return the perpendicular distance from a point to a line in 2D
  double distance(CPoint2D const& p)const;

  // a(xo)+b(y0)+c=0 => the point is on the CLine
  bool intersects(CPoint2D const &p)const;
  // return true if the lines as concurrent
  bool intersects(CLine2D const &l)const { return isConcurrent(l);}
  // return true if the lines intersects the ray
  bool intersects(CRay2D const &l)const;
  // return true if the lines intersects the linesegment
  bool intersects(CLineSegment2D const &l)const;
  bool intersects(const CAABBox2D& box)const;

  /**
  /if pt_br
  \brief Retorna Verdadeiro se a linha intersepta o AABBox2D passado. Se sim, computa os pontos de intersecção
  \elseif us_en
  \brief Return true if line intersects box. If so, compute intersection points.
  \endif
  **/
  bool intersection(const CAABBox2D& box, int *ICnt=NULL,CPoint2D* p0=NULL, CPoint2D* p1=NULL)const;
 
  /**
  \see intersection point of two line segments in 2 dimensions at http://paulbourke.net/geometry/pointlineplane/
  \return bitmask IntersecResult
  PARALLEL, COINCIDENT, NOT_INTERESECTING, INTERESECTING, INTERESECTING_EXTREMITY_P1, INTERESECTING_EXTREMITY_P2
  **/
  static int intersects_m1(CLine2D const &l1,CLine2D const &l2, CPoint2D intersection);

  /**
  \brief check intersection between one line and other line,segment or ray
  \see CLineSegment2D::intersects_m2
  **/
  static bool  intersects_m2(const CLine2D& this_line, const CLine2D& other_line, CPoint2D& intersection, int *intType=NULL,float EPS= 1e-3f);
  static bool  intersects_m2(const CLine2D& this_line, const CLineSegment2D& other_segment, CPoint2D& intersection, int *intType=NULL,float EPS= 1e-3f);
  static bool  intersects_m2(const CLine2D& this_line, const CRay2D& other_ray, CPoint2D& intersection, int *intType=NULL,float EPS= 1e-3f);
  
  static bool intersects( const CLine2D &line0, const CLine2D &line1, CPoint2D  &intersection_point,float maxDistance = 1e-6f );
 
  /// Computes the closest point on this segment to the given object
  CPoint2D closestPoint(CLineSegment2D const &l1)const;
  CPoint2D closestPoint(CPoint2D const  &p)const;
  CPoint2D closestPoint(CRay2D const &r1)const;
  // if the lines are parallel return an arbitrary point
  CPoint2D closestPoint(CLine2D const&other) const;



  static void test_line_2d();

};

namespace _2D{
#define l CLine2D

//: Return true iff line is the line at infinity
// \see  CLine2D
SMF_INLINE
bool is_ideal(l const&, float = (float)0) { return false; }

//: Are three lines concurrent, i.e., do they pass through a common point?
// \see  CLine2D
SMF_INLINE bool concurrent(l const& l1, l const& l2, l const& l3)
{
  return l1.a()*(l2.b()*l3.c()-l3.b()*l2.c())
        +l2.a()*(l3.b()*l1.c()-l1.b()*l3.c())
        +l3.a()*(l1.b()*l2.c()-l2.b()*l1.c())==0;
}

//: Write line description to stream: "<CLine2D ax+by+c>"
// \see  CLine2D

ostream&  operator<<(ostream& s, l const& line);

//: Read in three line parameters from stream
//  Either just reads three blank-separated numbers,
//  or reads three comma-separated numbers,
//  or reads three numbers in parenthesized form "(123, 321, -456)"
//  or reads a formatted line equation "123x+321y-456=0"
// \see  CLine2D

istream&  operator>>(istream& s, l& line);

#undef l
} //end 2D
} //end GEO
}  //end SMF

#endif // __SMF_2D_POINT_
