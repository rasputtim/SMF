#ifndef __SMF_2D_POINT_
#define __SMF_2D_POINT_

#include "../SMF_Config.h"
#include "../math/all.h"

namespace SMF {
using namespace MATH;
using namespace Util;
namespace GEO{

class CPoint2D;
class CLine2D;
class CLineSegment2D;
class CRay2D;
class CPolygon2D;
/**
 * \class CPoint2D
 *
 * \ingroup SMF_Geometric
 *
 * \if pt_br
 * \brief   Representa um ponto no sistema de coordenadas cartesiana 2D
 * \elseif us_en
 * \brief 	Represents a cartesian 2D point
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CPoint2D
{

 public:
 // the data associated with this point
  float x_;
  float y_;

  // Constructors/Initializers/Destructor------------------------------------

  /// Default constructor
  SMF_INLINE CPoint2D () {};

  /// Construct from two Types.
  SMF_INLINE CPoint2D (float px, float py) : x_(px), y_(py) {};

  /// Construct from 2-array.
  SMF_INLINE CPoint2D (float const v[2]) : x_(v[0]), y_(v[1]) {};

  /// Copy constructor from CVec2D
  SMF_INLINE CPoint2D(CVec2D const& p) : x_(p.x), y_(p.y) {}

  /// Copy constructor
  SMF_INLINE CPoint2D(CPoint2D const& p) : x_(p.x()), y_(p.y()) {}
  /// Destructor
  SMF_INLINE ~CPoint2D () {}
  /// Assignment
  SMF_INLINE CPoint2D & operator=(const CPoint2D & p)
  { x_ = p.x(); y_ = p.y(); return *this; }

  /// Test for equality
  SMF_INLINE bool operator==(CPoint2D const& p) const
  { return this==&p || (x_==p.x() && y_==p.y()); };

  SMF_INLINE bool operator!=(CPoint2D const& p)const
  { return !operator==(p); };

	/** \brief The difference of two points is the vector from second to first point
	 \see  CPoint2D */
	SMF_INLINE CVec2D operator-(CPoint2D const & p2)  		{ return CVec2D((x_-p2.x()), (y_-p2.y())); }
	SMF_INLINE CVec2D operator-(CPoint2D const & p2)const	{ return CVec2D(x_-p2.x(),y_-p2.y()); };

	/** \brief Adding a vector to a point gives a new point at the end of that vector
	// Note that vector + point is not defined!  It's always point + vector.
	// \see  CPoint2D */
	SMF_INLINE CVec2D operator+(CPoint2D const &p)			{ return CVec2D (x_+p.x(), y_+p.y()); }
	SMF_INLINE CVec2D operator+(CPoint2D const &p)const   { return CVec2D   (x_+p.x(), y_+p.y()); };



	SMF_INLINE CVec2D operator/(float const& p)const   { return CVec2D((x_/p),(y_/p)); };

	SMF_INLINE CVec2D operator*(float const& p)const   { return CVec2D((x_*p),(y_*p)); };

	SMF_INLINE CPoint2D midPoint(CPoint2D const &p)const		{ return CPoint2D ((x_+p.x_)*0.5f, (y_+p.y_)*0.5f); }

  SMF_INLINE bool compare(CPoint2D const &p)const { return *this==p;}
  // Data Access-------------------------------------------------------------

  SMF_INLINE float &x() {return x_;};
  SMF_INLINE float &y() {return y_;};

  SMF_INLINE float x() const {return x_;};
  SMF_INLINE float y() const {return y_;};

  /** \brief set \a x and \a y
  //  Note that \a x and \a y can also be set individually. */
  SMF_INLINE void set (float px, float py){ x_ = px; y_ = py; };
  /** \brief set \a x and \a y
  //  Note that \a x and \a y can also be set individually.*/
  SMF_INLINE void set (float const p[2]) { x_ = p[0]; y_ = p[1]; };
  //  Note that \a x  can also be set individually.
  SMF_INLINE void setX (float px){ x_ = px;  };
  //  Note that \a x  can also be set individually.
  SMF_INLINE void setY (float py){ y_ = py;  };

  /** \brief Return true if the point is at infinity (an ideal point).
  //  Always returns false.*/
  SMF_INLINE bool ideal(float = (float)0) const { return false; };

/**
\brief return if the point is a finite one
**/
SMF_INLINE_FORCED bool	CPoint2D::isFinite() const
{
	return MATH::isFinite(x_) && MATH::isFinite(y_);
}
/**
\brief turn this point to a infinite one
**/
SMF_INLINE_FORCED void toInfinite() {
	x_=CMath::NAN_FLOAT; y_=CMath::NAN_FLOAT;
}
/** \brief distane from this point tho the line segment*/
float distance(CLineSegment2D const& l);
float distance(CLine2D const& l);
float distance(CPoint2D const &p)const;
float distance(CPolygon2D const &poly, bool closed=true)const;
float distanceSqr(CPoint2D const &p)const;

/** \brief distance between point \a P(x,y) and closest point on polygon \a (px[i],py[i])*/
static float distanceToPolygon(CPolygon2D const&poly, float x, float y, bool closed=true);

/**
\param x1,y1 Coordenadas do ponto inicial do segmento
\param x2,y2 Coordenada do ponto final do segmento
\param x,y Coordenadas do ponto que se deseja calcuar a distância
\return a menor distância ao quadrado (evita Cálculo da raiz quadrada) entre o ponto e o segmento passados
**/
static float distanceTolinesegmentSqr(float x1, float y1, float x2, float y2, float x, float y);
/**
\if pt_br
\brief retorna o vetor posicional (ponto inicial na origem) que orresponde a este ponto
\elseif us_en
\brief return the Position Vector corresponding to this point
\endif
**/
SMF_INLINE CVec2D toVec2DPos() { return CVec2D(*this);}
SMF_INLINE CVec2D toVec2DPos()const { return CVec2D(*this);}

/**
\brief check if this point is aligned to the other two
\param p1 point to check alignment
\param p2 point to check alignment
\param tolerance
\return rue if the poinst are aligned
**/
bool IsCollinear(CPoint2D const &p1,CPoint2D const &p2, const float tolerance=CMath::EPSILON_High)const;

static void test_point_2d();
  /** \brief Read from stream, possibly with formatting
  //  Either just reads two blank-separated numbers,
  //  or reads two comma-separated numbers,
  //  or reads two numbers in parenthesized form "(123, 321)"*/
   istream& read(istream& is);
/** \brief Closest point to \a this point on the line segment \a l1
\return CPoint2D the closest point
\param l1the line
**/
CPoint2D closestPoint(CLineSegment2D const &l1)const ;
CPoint2D closestPoint(CLine2D const  &li)const;
CPoint2D closestPoint(CRay2D const &r1)const;
/**
\param ret_x, ret_y  Coordenadas do ponto mais próximo ao segmento
\param x1,y1 coordenadas do ponto inicial do segmento
\param x2,y2 coordenadas do ponto final do segmento
\param x0,y0 coordenadas do ponto que se quer achar o ponto mais proximo a ele
**/
static void closestPointToLinesegment(float& ret_x, float& ret_y,
                                      float x1, float y1,
                                      float x2, float y2,
                                      float x0, float y0);
/**
\brief transforma point (x,y) into a point (x',y') by theequations:
x' = ax+ cy + m
y' = bx + cy + n
\note this is used by conic transformations where a,b,c,d are:
\f$  T = \begin{bmatrix}a & b & 0\\c & d & 0\\m & n & 1\end{bmatrix}  \f$

**/
CPoint2D transform(CMat3D const &transMat)const;

/**
Get the orientation of the point Point regrdin the line segment formed by the other two points
\return LeftHandSide if the point is at left side from the segment
\return RightHandSide if the point is at right side from the segment, CollinearOrientation if the points are collineares
\return
**/
static int CPoint2D::orientation(const CPoint2D &px1,const CPoint2D &px2,const  CPoint2D &Py);

	/** \brief compare with epsilon
	**/
SMF_INLINE_FORCED bool compare( const CPoint2D &a, const float epsilon ) const {
	if ( CMath::fabs( x_ - a.x_ ) > epsilon ) {
		return false;
	}

	if ( CMath::fabs( y_ - a.y_ ) > epsilon ) {
		return false;
	}

	return true;
}

};

namespace _2D {

#define v MATH::CVec2D

//  +-+-+ CVec2D simple I/O +-+-+

/** \brief Write "<CVec2D x,y> " to stream
// \see  CVec2D*/
ostream& operator<<(ostream& s, v const& p);

/** \brief Read from stream, possibly with formatting
//  Either just reads two blank-separated numbers,
//  or reads two comma-separated numbers,
//  or reads two numbers in parenthesized form "(123, 321)"
// \see  CVec2D*/
istream& operator>>(istream& s, v& p);
#undef v

	//  +-+-+ point_2d simple I/O +-+-+

/** \brief Write "<CPoint2D x,y>" to stream
// \see  CPoint2D*/
ostream&  operator<<(ostream& s, CPoint2D const& p);

/** \brief Read from stream, possibly with formatting
//  Either just reads two blank-separated numbers,
//  or reads two comma-separated numbers,
//  or reads two numbers in parenthesized form "(123, 321)"
// \see  CPoint2D*/
istream&  operator>>(istream& s, CPoint2D & p);

/** \brief Adding a vector to a point gives a new point at the end of that vector
// Note that vector + point is not defined!  It's always point + vector.*/
SMF_INLINE CPoint2D operator+(CPoint2D const& p,
                             CVec2D const& v)
{ return CPoint2D ((p.x()+v.x), (p.y()+v.y)); }

/** \brief The difference of two points is the vector from second to first point
// \see  CPoint2D*/
//CVec2D operator-(const CPoint2D & p1,   const CPoint2D & p2)
//{ return CVec2D(p1.x()-p2.x(), p1.y()-p2.y()); }

/** \brief Adding a vector to a point gives the point at the end of that vector
// \see  CPoint2D*/
SMF_INLINE
CPoint2D & operator+=(CPoint2D & p,
                               CVec2D const& v)
{ p.set(p.x()+v.x, p.y()+v.y); return p; }

/** \brief Subtracting a vector from a point is the same as adding the inverse vector
// \see  CPoint2D*/
SMF_INLINE
CPoint2D operator-(CPoint2D const& p,
                             CVec2D const& v)
{ return p + (-v); }

/** \brief Subtracting a vector from a point is the same as adding the inverse vector
// \see  CPoint2D*/
SMF_INLINE
CPoint2D & operator-=(CPoint2D & p,
                               CVec2D const& v)
{ return p += (-v); }



//  +-+-+ point_2d geometry +-+-+

/** \brief cross ratio of four collinear points
// This number is projectively invariant, and it is the coordinate of p4
// in the reference frame where p2 is the origin (coordinate 0), p3 is
// the unity (coordinate 1) and p1 is the point at infinity.
// This cross ratio is often denoted as ((p1, p2; p3, p4)) (which also
// equals ((p3, p4; p1, p2)) or ((p2, p1; p4, p3)) or ((p4, p3; p2, p1)) )
// and is calculated as
//  \verbatim
//                      p1 - p3   p2 - p3      (p1-p3)(p2-p4)
//                      ------- : --------  =  --------------
//                      p1 - p4   p2 - p4      (p1-p4)(p2-p3)
// \endverbatim
// If three of the given points coincide, the cross ratio is not defined.
//
// In this implementation, a least-squares result is calculated when the
// points are not exactly collinear.
//
// \see  CPoint2D*/
float cross_ratio(CPoint2D const& p1, CPoint2D const& p2,
                   CPoint2D const& p3, CPoint2D const& p4);

/** \brief Are three points collinear, i.e., do they lie on a common line?
// \see  CPoint2D*/
bool collinear(CPoint2D const& p1,
               CPoint2D const& p2,
               CPoint2D const& p3);

/** \brief Return the relative distance to p1 wrt p1-p2 of p3.
//  The three points should be collinear and p2 should not equal p1.
//  This is the coordinate of p3 in the affine 1D reference frame (p1,p2).
//  If p3=p1, the ratio is 0; if p1=p3, the ratio is 1.
//  The mid point of p1 and p2 has ratio 0.5.
//  Note that the return type is double, not float, since the ratio of e.g.
//  two CVec2D<int> need not be an int.
// \see  CPoint2D*/
//SMF_INLINE double ratio(CPoint2D const& p1, CPoint2D const& p2, CPoint2D const& p3);


/** \brief Return the point at a given ratio wrt two other points.
//  By default, the mid point (ratio=0.5) is returned.
//  Note that the third argument is float, not double, so the midpoint of e.g.
//  two CPoint2D<int> is not a valid concept.  But the reflection point
//  of p2 wrt p1 is: in that case f=-1.
// \see  CPoint2D*/
CPoint2D midpoint(CPoint2D const& p1,
                            CPoint2D const& p2,
                            float f = (float)0.5);


/** \brief Return the point at the centre of gravity of two given points.
// Identical to midpoint(p1,p2).
// \see  CPoint2D*/
CPoint2D centre(CPoint2D const& p1,
                          CPoint2D const& p2);

/** \brief Return the point at the centre of gravity of three given points.
// \see  CPoint2D*/
CPoint2D centre(CPoint2D const& p1,
                          CPoint2D const& p2,
                          CPoint2D const& p3);

/** \brief Return the point at the centre of gravity of four given points.
// \see  CPoint2D*/
CPoint2D centre(CPoint2D const& p1,
                          CPoint2D const& p2,
                          CPoint2D const& p3,
                          CPoint2D const& p4);
/** \brief Return the point at the centre of gravity of a set of given points.
// Beware of possible rounding errors when float is e.g. int.
// \see  CPoint2D*/
CPoint2D centre(std::vector<CPoint2D> const& v);


} //end 2D
} //end GEO
}  //end SMF

#endif // __SMF_2D_POINT_
