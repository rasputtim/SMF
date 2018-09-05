#ifndef __SMF_2D_POINT_HOM_
#define __SMF_2D_POINT_HOM_

#include "SMF_Config.h"
#include "math/all.h"

namespace SMF {
using namespace MATH;
using namespace Util;
namespace GEO{

class CPoint2D;
class CLine2D;
class CLine2DHom;

/**
 * \class CPoint2DHom
 *
 * \ingroup SGF_Geometry
 *
 * \if pt_br
 * \brief
 * \elseif us_en
 * \brief 	point in projective 2D space

 * \endif
 * \author (last to touch it) $Autor: Salvatore Giannotta Filho $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: a_materasu@hotmail.com
 *
 */
//: Represents a cartesian 2D point

//: Represents a homogeneous 2D point
class CPoint2DHom
{
  // the data associated with this point
  sf_s16 x_;
  sf_s16 y_;
  sf_s16 w_;

 public:

  // Constructors/Initializers/Destructor------------------------------------

  //: Default constructor with (0,0,1)
  inline CPoint2DHom() : x_(0), y_(0), w_((sf_s16)1) {}

  //: Construct from two (nonhomogeneous) or three (homogeneous) Types.
  inline CPoint2DHom(sf_s16 px, sf_s16 py, sf_s16 pw = (sf_s16)1)
    : x_(px), y_(py), w_(pw) {}

  //: Construct from homogeneous 3-array.
  inline CPoint2DHom(const sf_s16 v[3]) : x_(v[0]), y_(v[1]), w_(v[2]) {}

  //: Construct point at infinity from direction vector.
  inline CPoint2DHom(CVec2D const& v) : x_(v.x), y_(v.y), w_(0) {}

  //: Construct from (non-homogeneous) CPoint2D
  inline explicit CPoint2DHom(CPoint2D const& p)
    : x_(p.x(), y_(p.y()), w_((sf_s16)1) {}

  //: Construct from 2 lines (intersection).
  CPoint2DHom(CLine2DHom const& l1,
                    CLine2DHom const& l2);

#if 0
  // Default copy constructor
  inline CPoint2DHom(const CPoint2DHom<sf_s16>& p)
    : x_(p.x()), y_(p.y()), w_(p.w()) {}

  // Destructor
  inline ~CPoint2DHom() {}

  // Default assignment operator
  inline CPoint2DHom<sf_s16>& operator=(CPoint2DHomconst& p)
  {
    set(p.x(),p.y(),p.w());
    return *this;
  }
#endif

  //: the comparison operator
  inline bool operator==(CPoint2DHom const& p) const
  {
    return (this==&p) ||
           (x()*p.w()==w()*p.x() && y()*p.w()==w()*p.y() && y()*p.x()==x()*p.y());
  }

  inline bool operator!=(CPoint2DHom const& other)const{return !operator==(other);}

  // Data Access-------------------------------------------------------------

  inline sf_s16 x() const { return x_; }
  inline sf_s16 y() const { return y_; }
  inline sf_s16 w() const { return w_; }

  //: Set \a x,y,w
  // Note that it does not make sense to set \a x, \a y or \a w individually.
  inline void set(sf_s16 px, sf_s16 py, sf_s16 pw = (sf_s16)1)
  { x_ = px, y_ = py, w_ = pw; }

  inline void set(sf_s16 const p[3]) { x_ = p[0]; y_ = p[1]; w_ = p[2]; }

  //: Return true iff the point is at infinity (an ideal point).
  // The method checks whether |w| <= tol * max(|x|,|y|)
  inline bool ideal(sf_s16 tol = (sf_s16)0) const
  {
#define _Abs(x) (x<0?-x:x) // avoid #include of cmath.h
    return _Abs(w()) <= tol * _Abs(x()) ||
           _Abs(w()) <= tol * _Abs(y());
#undef _Abs
  }
};

#if 0
//  +-+-+ point_2d simple I/O +-+-+

//: Write "<CPoint2DHom (x,y,w) >" to stream
// \relatesalso CPoint2DHom
template <class sf_s16>
vcl_ostream& operator<<(vcl_ostream& s, CPoint2DHom const& p);

//: Read x y w from stream
// \relatesalso CPoint2DHom
template <class sf_s16>
vcl_istream& operator>>(vcl_istream& s, CPoint2DHom<sf_s16>& p);

//  +-+-+ homg_point_2d arithmetic +-+-+

//: Return true iff the point is at infinity (an ideal point).
// The method checks whether |w| <= tol * max(|x|,|y|)
// \relatesalso CPoint2DHom
template <class sf_s16> inline
bool is_ideal(CPoint2DHom const& p, sf_s16 tol=(sf_s16)0){return p.ideal(tol);}

//: The difference of two points is the vector from second to first point
// This function is only valid if the points are not at infinity.
// \relatesalso CPoint2DHom
template <class sf_s16> inline
CVec2D<sf_s16> operator-(CPoint2DHom const& p1,
                           CPoint2DHom const& p2)
{
  SMF_ASSERT(p1.w() && p2.w());
  return CVec2D<sf_s16>(p1.x()/p1.w()-p2.x()/p2.w(),
                          p1.y()/p1.w()-p2.y()/p2.w());
}

//: Adding a vector to a point gives a new point at the end of that vector
// If the point is at infinity, nothing happens.
// Note that vector + point is not defined!  It's always point + vector.
// \relatesalso CPoint2DHom
template <class sf_s16> inline
CPoint2DHom operator+(CPoint2DHom const& p,
                               CVec2D<sf_s16> const& v)
{ return CPoint2DHom<sf_s16>(p.x()+v.x()*p.w(), p.y()+v.y()*p.w(), p.w()); }

//: Adding a vector to a point gives the point at the end of that vector
// If the point is at infinity, nothing happens.
// \relatesalso CPoint2DHom
template <class sf_s16> inline
CPoint2DHom<sf_s16>& operator+=(CPoint2DHom<sf_s16>& p,
                                 CVec2D<sf_s16> const& v)
{ p.set(p.x()+v.x()*p.w(), p.y()+v.y()*p.w(), p.w()); return p; }

//: Subtracting a vector from a point is the same as adding the inverse vector
// \relatesalso CPoint2DHom
template <class sf_s16> inline
CPoint2DHom operator-(CPoint2DHom const& p,
                               CVec2D<sf_s16> const& v)
{ return p + (-v); }

//: Subtracting a vector from a point is the same as adding the inverse vector
// \relatesalso CPoint2DHom
template <class sf_s16> inline
CPoint2DHom<sf_s16>& operator-=(CPoint2DHom<sf_s16>& p,
                                 CVec2D<sf_s16> const& v)
{ return p += (-v); }

//  +-+-+ homg_point_2d geometry +-+-+

//: cross ratio of four collinear points
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
// \relatesalso CPoint2DHom
//
template <class sf_s16>
double cross_ratio(CPoint2DHomconst& p1, CPoint2DHomconst& p2,
                   CPoint2DHomconst& p3, CPoint2DHomconst& p4);

//: Are three points collinear, i.e., do they lie on a common line?
// \relatesalso CPoint2DHom
template <class sf_s16> inline
bool collinear(CPoint2DHom const& p1,
               CPoint2DHom const& p2,
               CPoint2DHom const& p3)
{
  return (p1.x()*p2.y()-p1.y()*p2.x())*p3.w()
        +(p3.x()*p1.y()-p3.y()*p1.x())*p2.w()
        +(p2.x()*p3.y()-p2.y()*p3.x())*p1.w()==0;
}

//: Return the relative distance to p1 wrt p1-p2 of p3.
//  The three points should be collinear and p2 should not equal p1.
//  This is the coordinate of p3 in the affine 1D reference frame (p1,p2).
//  If p3=p1, the ratio is 0; if p1=p3, the ratio is 1.
//  The mid point of p1 and p2 has ratio 0.5.
//  Note that the return type is double, not sf_s16, since the ratio of e.g.
//  two CVec2D<int> need not be an int.
// \relatesalso CPoint2DHom
template <class sf_s16> inline
double ratio(CPoint2DHom const& p1,
             CPoint2DHom const& p2,
             CPoint2DHom const& p3)
{ return (p3-p1)/(p2-p1); }

//: Return the point at a given ratio wrt two other points.
//  By default, the mid point (ratio=0.5) is returned.
//  Note that the third argument is sf_s16, not double, so the midpoint of e.g.
//  two CPoint2DHom<int> is not a valid concept.  But the reflection point
//  of p2 wrt p1 is: in that case f=-1.
template <class sf_s16> inline
CPoint2DHom midpoint(CPoint2DHom const& p1,
                              CPoint2DHom const& p2,
                              sf_s16 f = (sf_s16)0.5)
{ return p1 + f*(p2-p1); }


//: Return the point at the centre of gravity of two given points.
// Identical to midpoint(p1,p2).
// Invalid when both points are at infinity.
// If only one point is at infinity, that point is returned.
// \relatesalso CPoint2DHom
CPoint2DHom centre(CPoint2DHom const& p1,
                            CPoint2DHom const& p2)
{
  return CPoint2DHom<sf_s16>(p1.x()*p2.w() + p2.x()*p1.w(),
                              p1.y()*p2.w() + p2.y()*p1.w(),
                              p1.w()*p2.w()*2 );
}

//: Return the point at the centre of gravity of a set of given points.
// There are no rounding errors when sf_s16 is e.g. int, if all w() are 1.
// \relatesalso CPoint2DHom
template <class sf_s16> inline
CPoint2DHom centre(vcl_vector<CPoint2DHom > const& v)
{
  int n=v.size();
  SMF_ASSERT(n>0); // it is *not* correct to return the point (0,0) when n==0.
  sf_s16 x = 0, y = 0;
  for (int i=0; i<n; ++i) x+=v[i].x()/v[i].w(), y+=v[i].y()/v[i].w();
  return CPoint2DHom<sf_s16>(x,y,(sf_s16)n);
}
#endif
} //end GEO
}  //end SMF

#endif // __SMF_2D_POINT_
