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

#include "SMF_Config.h"
#include "geometry/SMF_2DPoint.h"
#include "geometry/SMF_2DLine.h"
#include "geometry/SMF_2DAABBox.h"
#include "geometry/SMF_2DLineSegment.h"
#include "math/SMF_Vector.h"
#include "math/all.h"
#include "geometry/all.h"
#include "util/SMF_TestLib.h"
#include <vector>
#include <string.h>
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream
#include <iostream>     // std::cout, std::right, std::endl
#include <iomanip>      // std::setw

namespace SMF {
using namespace MATH;
using namespace Util;

namespace GEO{



void CPoint2D::test_point_2d()
{
	using GEO::_2D::operator<<;
	using GEO::_2D::operator>>;
  
  float d[] = {5,5};
  CPoint2D p1(3,7), p2(d), p3(-1,-8);
   Debug::debug(Debug::math,__FUNCTION__) << p3 << endl;

  TEST("constructor", p1.x() == 3 && p1.y()==7, true);

  TEST("inequality", (p1 != p3), true);

  p3.set(3,7);
  TEST("set", p3.x() == 3 && p3.y()==7, true);

  TEST("equality", (p1 == p3), true);

  CVec2D d1 = (p1 - p2);
  
  TEST("sum; difference", CPoint2D(p2+d1), p1);
  using SMF::GEO::_2D::operator+=;
  TEST("+=", (p2+=d1), p1);
  TEST("+=", p2, p1);

  p2.set(4,5);
  p3.set(7,-1);
  bool b = _2D::collinear(p1,p2,p3);
  TEST("collinear", b, true);
  //double r = _2D::ratio(p1,p2,p3);
  //TEST("ratio", r, 4.0);
  CPoint2D m = _2D::midpoint(p1,p2,4);
  TEST("midpoint", m, p3);

  CPoint2D c = _2D::centre(p1,p3);
  CPoint2D cc(5,3);
  TEST("centre", c, cc);
  c = _2D::centre(p1,p2,c); // assignment
  TEST("centre", c, p2);
  c = _2D::centre(p2,p3,cc,p2);
  TEST("centre", c, cc);
  double r;
  r = _2D::cross_ratio(p1,p2,c,p3);
  TEST("cross_ratio", r, 1.5);

  //CLine2D l1(3,4,5), l2(3,2,1);
  //CPoint2D pi(l1,l2); // intersection
  //CPoint2D pp(1,-2);
  //TEST("intersection", pi, pp);

  //CLine2D l3(1,2,3), l4(3,2,1);
  //CPoint2D pj(l3,l4); // intersection
  //CPoint2D pq(1,-2);
  //TEST("intersection", pj, pq);

  //TEST("distance_origin", distance_origin(l1), 1);
 
  {
   std::stringstream is; 
	is << "4.4 -5 7e1";
    CPoint2D p(0.0,0.0); is >> p;
	Debug::debug(Debug::math,__FUNCTION__) << p << endl;
    TEST("istream CPoint2D (blank-separated)", p, CPoint2D(4.4,-5));
  }

  {
    stringstream is; is << "7e1, 11 , blabla";
    CPoint2D p(0.0,0.0); is >> p;
    Debug::debug(Debug::math,__FUNCTION__) << p << endl;
    TEST("istream CPoint2D (comma-separated)", p, CPoint2D(70,11));
  }
  {
    stringstream is; is << " (12,13 ) !";
    CPoint2D p(0.0,0.0); is >> p;
    Debug::debug(Debug::math,__FUNCTION__) << p << endl;
    TEST("istream CPoint2D (parenthesized)", p, CPoint2D(12,13));
  }

}

//: true if the point lies on the line segment and is between the endpoints
// \see  CPoint2D

bool testPoint2D(CPoint2D const& p,CLineSegment2D const& lseg, double epsilon=1e-14f)
{
  CPoint2D p1 = lseg.point1(), p2 = lseg.point2();
  float x1 = p1.x(), y1 = p1.y(),
    x2 = p2.x(), y2 = p2.y(),
    xp = p.x(),  yp = p.y();
  // compute squared distances
  float d1p = (xp-x1)*(xp-x1) + (yp-y1)*(yp-y1);
  float d2p = (xp-x2)*(xp-x2) + (yp-y2)*(yp-y2);
  float d12 = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1);
  double diff = CMath::sqrt(d1p) + CMath::sqrt(d2p) - CMath::sqrt(d12);
  // diff is always >= 0 (triangle inequality)
  return diff <= epsilon;
}


//: Compute discriminant function
// Returns determinant of
// \verbatim
// [ x1 x2 x3 ]
// [ y1 y2 y3 ]
// [ 1  1  1  ]
// \endverbatim

float testTriandle2Ddiscriminant(float x1, float y1,
                                 float x2, float y2,
                                 float x3, float y3)
{
  return x1*(y2-y3) - x2*(y1-y3) + x3*(y1-y2);
}
// Returns true iif (x3,y3) lies inbetween (x1,y1) and (x2,y2);
// all three points are assumed to be collinear
static 
bool inbetween(float x1, float y1, float x2, float y2, float x3, float y3)
{
  return (x1-x3)*(x2-x3)<=0 && (y1-y3)*(y2-y3)<=0;
}

bool testLineSegment2D(float x1, float y1, float x2, float y2, float x3, float y3, float x4, float y4)
{
  // reduce precision of inputs so that signs of testTriangle2Ddiscriminants
  // `a', `b', `c', and `d' are more stable when they are very close to zero.

  double px1 = x1;
  double py1 = y1; 
  double px2 = x2;
  double py2 = y2;
  double px3 = x3;
  double py3 = y3;
  double px4 = x4;
  double py4 = y4;

  px1 = (px1 + px1*1e4) - px1*1e4;
  py1 = (py1 + py1*1e4) - py1*1e4;
  px2 = (px2 + px2*1e4) - px2*1e4;
  py2 = (py2 + py2*1e4) - py2*1e4;
  px3 = (px3 + px3*1e4) - px3*1e4;
  py3 = (py3 + py3*1e4) - py3*1e4;
  px4 = (px4 + px4*1e4) - px4*1e4;
  py4 = (py4 + py4*1e4) - py4*1e4;

  // two lines intersect if p1 and p2 are on two opposite sides of the line p3-p4
  // and p3 and p4 are on two opposite sides of line p1-p2
  // degenerate cases (collinear points) are tricky to handle

  double a = testTriandle2Ddiscriminant(px1, py1, px2, py2, px3, py3);
  double b = testTriandle2Ddiscriminant(px1, py1, px2, py2, px4, py4);
  double c = testTriandle2Ddiscriminant(px3, py3, px4, py4, px1, py1);
  double d = testTriandle2Ddiscriminant(px3, py3, px4, py4, px2, py2);

  // force to be zero when they're close to zero
  a = (CMath::fabs(a) < 1e-12) ? 0 : a;
  b = (CMath::fabs(b) < 1e-12) ? 0 : b;
  c = (CMath::fabs(c) < 1e-12) ? 0 : c;
  d = (CMath::fabs(d) < 1e-12) ? 0 : d;

  return
    ( ( (a<=0 && b>0) || (a>=0 && b<0) || (a<0 && b>=0) || (a>0 && b<=0) ) &&
      ( (c<=0 && d>0) || (c>=0 && d<0) || (c<0 && d>=0) || (c>0 && d<=0) ) )
    ||
    ( // the above two conditions are only sufficient for noncollinear line segments! - PVr
      a == 0 && b == 0 && c == 0 && d == 0 &&
      ( inbetween(px1, py1, px2, py2, px3, py3) ||
        inbetween(px1, py1, px2, py2, px4, py4) ||
        inbetween(px3, py3, px4, py4, px1, py1) ||
        inbetween(px3, py3, px4, py4, px2, py2) )
    );
}

bool testLineSegment2D(CLineSegment2D const& l1,
                                     CLineSegment2D const& l2)
{
  return testLineSegment2D(l1.point1().x(),l1.point1().y(),
                                  l1.point2().x(),l1.point2().y(),
                                  l2.point1().x(),l2.point1().y(),
                                  l2.point2().x(),l2.point2().y());
}


bool TestClipLineToAABB(float a, float b, float c, // line coefficients.
                          float x1, float y1,    // bounding
                          float x2, float y2,    // box.
                          float &bx, float &by,  // start and
                          float &ex, float &ey)  // end points.
{
  if (x1>x2) Swap(x1,x2);
  if (y1>y2) Swap(y1,y2);
  // now x1 <= x2 and y1 <= y2

  if (a == 0 && b == 0) return false; // then ax+by+c=0 is the line at infinity

  bool b_set = false, // has the point (bx,by) been set to a valid point?
       e_set = false; // has the point (ex,ey) been set to a valid point?

  if (a != 0) // line is not horizontal
  {
    // intersection point with the line y=y1:
    by = y1; bx = -(b*y1+c)/a;
    // intersection point with the line y=y2:
    ey = y2; ex = -(b*y2+c)/a;

    b_set =  bx >= x1 && bx <= x2; // does this intersection point
    e_set =  ex >= x1 && ex <= x2; // lie on the bounding box?
  }

  if (b_set && e_set) return true;
  if (b_set) { Swap(bx,ex); Swap(by,ey); Swap(b_set,e_set); }
  // now b_set is false

  if (b != 0) // line is not vertical
  {
    // intersection point with the line x=x1:
    bx = x1; by = -(a*x1+c)/b;
    b_set =  by >= y1 && by <= y2;
    if (b_set && e_set) return true;
    if (b_set) { Swap(bx,ex); Swap(by,ey); e_set=true; }

    // intersection point with the line x=x2:
    bx = x2; by = -(a*x2+c)/b;
    b_set =  by >= y1 && ey <= y2;
  }

  return b_set && e_set;
}


//: clip given line to given box, and return resulting line segment
// \see  CLine2D
// \see  AABBox2D

CLineSegment2D TestClipLineToAABB(CLine2D const& l,
                                            CAABBox2D const& b)
{
  float sx, sy, ex, ey;
  bool r = TestClipLineToAABB(l.a(), l.b(), l.c(),
                                b.minX(), b.minY(), b.maxX(), b.maxY(),
                                sx, sy, ex, ey);
  return r ? CLineSegment2D(CPoint2D(sx, sy), CPoint2D(ex, ey))
           : CLineSegment2D(); // uninitialised when no intersection
}

void CLine2D::test_line_2d()
{
	using GEO::_2D::operator<<;
	using GEO::_2D::operator>>;

  float d[] = {5,5,-1};
  CLine2D l1(3,7,0), l2(d), l3(0,-1,-8);
   Debug::debug(Debug::math,__FUNCTION__) << l3 << endl;

  TEST("inequality", (l1 != l3), true);

  l3.set(3,7,0);
  TEST("equality", (l1 == l3), true);

  l2.set(3,4,0);
  l3.set(7,-1,0);
  bool b = _2D::concurrent(l1,l2,l3); // because they share the point (0,0)
  TEST("concurrent", b, true);

  CPoint2D p00(0,0), p01(0,1), p02(0,2), p03(0,3),
                       p10(1,0), p11(1,1), p12(1,2), p13(1,3),
                       p20(2,0), p21(2,1), p30(3,0), p31(3,1);
  CLine2D li(p10,p01); // line through these two points
   Debug::debug(Debug::math,__FUNCTION__) << li << endl;
  CLine2D ll(1,1,-1);
  TEST("join", li, ll);

  l3.set(0,-1,-8);
  TEST("distance(line,point)",l3.distance(p03), 11);
  TEST_NEAR("distance(line,point)", l2.distance(p03), 2.4, 1e-9);
  TEST_NEAR("distance(point,line)", p03.distance(l2), 2.4, 1e-9);

  CLineSegment2D ls(p10,p01); // line segment through these two points
  TEST("line segment first end point", ls.point1(), p10);
  TEST("line segment second end point", ls.point2(), p01);

  CLineSegment2D ls2(p01,p10); // inverse line segment through these points
  TEST("line segment first end point", ls2.point1(), p01);
  TEST("line segment second end point", ls2.point2(), p10);

  // line segment intersections:
  // 1. with points:
  CLineSegment2D ls1(p12,p30); // line segment containing p21; its support line contains p03.
  TEST("end point is coincident", testPoint2D(p12,ls1), true);
  TEST("mid point is coincident", testPoint2D(p21,ls1), true);
  TEST("point on support line is not coincident", testPoint2D(p03,ls1), false);
  TEST("arbitrary point is not coincident", testPoint2D(p11,ls1), false);
  TEST("very nearby point is not coincident", testPoint2D(CPoint2D(2.01,1.0),ls1), false);
  // 2. with each other:
  TEST("identical line segments: intersect", testLineSegment2D(ls,ls2), true);
  TEST("identical line segments: intersect",
       testLineSegment2D(CLineSegment2D(p01, p11), CLineSegment2D(p11, p01)), true);
  TEST("disjoint horizontal collinear line segments: do not intersect",
       testLineSegment2D(CLineSegment2D(p01, p11), CLineSegment2D(p21, p31)), false);
  TEST("touching horizontal collinear line segments: intersect",
       testLineSegment2D(CLineSegment2D(p01, p11), CLineSegment2D(p11, p21)), true);
  TEST("touching horizontal collinear line segments: intersect",
       testLineSegment2D(CLineSegment2D(p01, p11), CLineSegment2D(p21, p11)), true);
  TEST("overlapping horizontal collinear line segments: intersect",
       testLineSegment2D(CLineSegment2D(p01, p21), CLineSegment2D(p11, p31)), true);
  TEST("internally touching horizontal collinear line segments: intersect",
       testLineSegment2D(CLineSegment2D(p01, p21), CLineSegment2D(p21, p11)), true);
  TEST("internally overlapping horizontal collinear line segments: intersect",
       testLineSegment2D(CLineSegment2D(p01, p31), CLineSegment2D(p21, p11)), true);
  TEST("parallel horizontal line segments: do not intersect",
       testLineSegment2D(CLineSegment2D(p01, p11), CLineSegment2D(p12, p02)), false);
  TEST("disjoint vertical collinear line segments: do not intersect",
       testLineSegment2D(CLineSegment2D(p10, p11), CLineSegment2D(p12, p13)), false);
  TEST("touching vertical collinear line segments: intersect",
       testLineSegment2D(CLineSegment2D(p10, p11), CLineSegment2D(p11, p12)), true);
  TEST("touching vertical collinear line segments: intersect",
       testLineSegment2D(CLineSegment2D(p10, p11), CLineSegment2D(p12, p11)), true);
  TEST("overlapping vertical collinear line segments: intersect",
       testLineSegment2D(CLineSegment2D(p10, p12), CLineSegment2D(p11, p13)), true);
  TEST("internally touching vertical collinear line segments: intersect",
       testLineSegment2D(CLineSegment2D(p10, p12), CLineSegment2D(p12, p11)), true);
  TEST("internally overlapping vertical collinear line segments: intersect",
       testLineSegment2D(CLineSegment2D(p10, p13), CLineSegment2D(p12, p11)), true);
  TEST("parallel vertical line segments: do not intersect",
       testLineSegment2D(CLineSegment2D(p10, p11), CLineSegment2D(p21, p20)), false);
  TEST("disjoint oblique collinear line segments: do not intersect",
       testLineSegment2D(CLineSegment2D(p30, p21), CLineSegment2D(p12, p03)), false);
  TEST("touching oblique collinear line segments: intersect",
       testLineSegment2D(CLineSegment2D(p20, p11), CLineSegment2D(p11, p02)), true);
  TEST("touching oblique collinear line segments: intersect",
       testLineSegment2D(CLineSegment2D(p20, p11), CLineSegment2D(p02, p11)), true);
  TEST("overlapping oblique collinear line segments: intersect",
       testLineSegment2D(CLineSegment2D(p30, p12), CLineSegment2D(p21, p03)), true);
  TEST("internally touching oblique collinear line segments: intersect",
       testLineSegment2D(CLineSegment2D(p30, p12), CLineSegment2D(p12, p21)), true);
  TEST("internally overlapping oblique collinear line segments: intersect",
       testLineSegment2D(CLineSegment2D(p30, p03), CLineSegment2D(p12, p21)), true);
  TEST("parallel oblique line segments: do not intersect",
       testLineSegment2D(CLineSegment2D(p20, p11), CLineSegment2D(p21, p30)), false);
  TEST("disjoint noncollinear line segments: do not intersect",
       testLineSegment2D(CLineSegment2D(p10, p21), CLineSegment2D(p00, p01)), false);
  TEST("touching noncollinear line segments: intersect",
       testLineSegment2D(CLineSegment2D(p10, p21), CLineSegment2D(p10, p01)), true);
  TEST("intersecting noncollinear line segments: intersect",
       testLineSegment2D(CLineSegment2D(p00, p21), CLineSegment2D(p01, p20)), true);
  // Degenerate cases: line segments of length 0
  TEST("point on horizontal line segment: intersects",
       testLineSegment2D(CLineSegment2D(p01, p21), CLineSegment2D(p11, p11)), true);
  TEST("point collinear with horizontal line segment: does not intersect",
       testLineSegment2D(CLineSegment2D(p01, p21), CLineSegment2D(p31, p31)), false);
  TEST("point not collinear with horizontal line segment: does not intersect",
       testLineSegment2D(CLineSegment2D(p01, p21), CLineSegment2D(p12, p12)), false);
  TEST("point on vertical line segment: intersects",
       testLineSegment2D(CLineSegment2D(p10, p12), CLineSegment2D(p11, p11)), true);
  TEST("point collinear with vertical line segment: does not intersect",
       testLineSegment2D(CLineSegment2D(p10, p12), CLineSegment2D(p13, p13)), false);
  TEST("point not collinear with vertical line segment: does not intersect",
       testLineSegment2D(CLineSegment2D(p10, p12), CLineSegment2D(p21, p21)), false);
  TEST("point on oblique line segment: intersects",
       testLineSegment2D(CLineSegment2D(p03, p21), CLineSegment2D(p12, p12)), true);
  TEST("point collinear with oblique line segment: does not intersect",
       testLineSegment2D(CLineSegment2D(p03, p21), CLineSegment2D(p30, p30)), false);
  TEST("point not collinear with oblique line segment: does not intersect",
       testLineSegment2D(CLineSegment2D(p03, p21), CLineSegment2D(p11, p11)), false);
  TEST("two identical points: intersect",
       testLineSegment2D(CLineSegment2D(p01, p01), CLineSegment2D(p01, p01)), true);
  TEST("two different points: does not intersect",
       testLineSegment2D(CLineSegment2D(p01, p01), CLineSegment2D(p12, p12)), false);

  TEST("four almost collinear points: do not intersect",
       testLineSegment2D(328.99996948242187,127.99999237060547,271.98755895963365,178.20118501723579,
                                        181.43250686733217,257.93771209135798,102.99999237060547,326.99996948242187), false);

  CAABBox2D bx(0,2,0,3);
  CLineSegment2D ls3 = TestClipLineToAABB(li,bx);
  Debug::debug(Debug::math,__FUNCTION__) << ls3 << '\n';
  TEST("line segment equality", ls3, ls);
  TEST("line segment equality", ls3, ls2);
  //b must be !=zero;
  CLine2D lp1(2,-3,5);
  CLine2D lp2(9,6,1);
  CLine2D lp3(3,0,-7);
  CLine2D lp4(7,0,1);

  TEST("normalize", l2.normalize(), true);
  TEST_NEAR("normalize: a()", l2.a(), 0.6, 1e-12);
  TEST_NEAR("normalize: b()", l2.b(), 0.8, 1e-12);
  TEST("perpendicular",lp1.isPerpendicular(lp2),true);
  TEST("parallel",lp1.isParallel(lp2),false);
  TEST("concurrents",lp1.isConcurrent(lp2),true);
  TEST("perpendicular",lp3.isPerpendicular(lp4),false);
  TEST("parallel",lp3.isParallel(lp4),true);
  TEST("concurrents",lp3.isConcurrent(lp4),false);




  {
    stringstream is; is << "4.5 -5 7e1  9x+7y-8=0";
    CLine2D l; is >> l;
    TEST("istream CLine2D", l, CLine2D(4.5f,-5,70));
    is >> l;
    TEST("istream CLine2D formatted", l, CLine2D(9,7,-8));
  }

  stringstream is; is << "\n4 6 7 9";
  CLineSegment2D l_s; is >> l_s;
  TEST("istream line_segment_2d", l_s, CLineSegment2D(CPoint2D(4,6), CPoint2D(7,9)));


}

void CAABBox2D::test_box_2d()
{
  	using GEO::_2D::operator<<;
	using GEO::_2D::operator>>;
// Create empty AABBox2D
  CAABBox2D b;
  Debug::debug(Debug::math,__FUNCTION__) << b << endl;
  TEST("AABBox2D is empty", b.isEmpty(), true);
  TEST("AABBox2D has no volume", b.volume(), 0.0);

  CPoint2D p0(0,0), p1(1,0), p2(0,1), p12(1,1);
  TEST("!contains(p0)", b.contains(p0), false);
  TEST("!contains(p1)", b.contains(p1), false);
  TEST("!contains(p2)", b.contains(p2), false);
  TEST("!contains(p12)", b.contains(p12), false);

  b.add(p0); 
  Debug::debug(Debug::math,__FUNCTION__) << b << endl;
  TEST("contains(p0)", b.contains(p0), true);
  TEST("!contains(p1)", b.contains(p1), false);
  TEST("!contains(p2)", b.contains(p2), false);
  TEST("!contains(p12)", b.contains(p12), false);
  TEST("AABBox2D is not empty", b.isEmpty(), false);
  TEST("AABBox2D has no volume", b.volume(), 0.0);
  TEST("centroid", b.centroid(), p0);
  TEST("minPoint", b.minPoint(), p0);
  TEST("maxPoint", b.maxPoint(), p0);

  b.add(p1); 
  Debug::debug(Debug::math,__FUNCTION__) << b << endl;
  TEST("contains(p0)", b.contains(p0), true);
  TEST("contains(p1)", b.contains(p1), true);
  TEST("!contains(p2)", b.contains(p2), false);
  TEST("!contains(p12)", b.contains(p12), false);
  TEST("AABBox2D is not empty", b.isEmpty(), false);
  TEST("AABBox2D has no volume", b.volume(), 0.0);
  TEST("centroid", b.centroid(), p0 + (p1-p0)*0.5f);
  TEST("minPoint", b.minPoint(), p0);
  TEST("maxPoint", b.maxPoint(), p1);

  b.add(p2); 
  Debug::debug(Debug::math,__FUNCTION__) << b << endl;
  TEST("contains(p0)", b.contains(p0), true);
  TEST("contains(p1)", b.contains(p1), true);
  TEST("contains(p2)", b.contains(p2), true);
  TEST("contains(p12)", b.contains(p12), true);
  TEST("AABBox2D is not empty", b.isEmpty(), false);
  TEST("AABBox2D has volume 1", b.volume(), 1.0);
  TEST("centroid", b.centroid(), p1 + (p2-p1)*0.5f);
  TEST("minPoint", b.minPoint(), p0);
  TEST("maxPoint", b.maxPoint(), p12);

  CPoint2D p4(0.5,0.5), p5(2,2);
  CAABBox2D b2(p4,p5);
  TEST("AABBox2D has volume 2.25", b2.volume(), 2.25);
  TEST("!contains(b2)", b.contains(b2), false);
  CAABBox2D b3(p5,p4);
  TEST("boxes are equal", b3, b2);
  b.add(b2); 
  Debug::debug(Debug::math,__FUNCTION__) << b << endl;
  TEST("union AABBox2D has volume 4", b.volume(), 4.0);
  TEST("contains(b2)", b.contains(b2), true);
  TEST("centroid", b.centroid(), p12);
  TEST("minPoint", b.minPoint(), p0);
  TEST("maxPoint", b.maxPoint(), p5);

  b2=b; b2.setCentroid(p0); 
  Debug::debug(Debug::math,__FUNCTION__) << b2 << endl;
  TEST("setCentroid", b2.centroid(), p0);
  TEST("volume did not change", b2.volume(), 4.0);
  TEST("minPoint", b2.minPoint(), CPoint2D(-1,-1));
  TEST("maxPoint", b2.maxPoint(), p12);

  b2.setWidth(1.0); 
  Debug::debug(Debug::math,__FUNCTION__) << b2 << endl;
  TEST("setWidth", b2.centroid(), p0);
  TEST("volume is now 2", b2.volume(), 2.0);
  b2.setHeight(1.0); 
  Debug::debug(Debug::math,__FUNCTION__) << b2 << endl;
  TEST("setHeight", b2.centroid(), p0);
  TEST("volume is now 1", b2.volume(), 1.0);
  b2.scaleAboutCentroid(2.0);
  TEST("scaleAboutCentroid", b2.centroid(), p0);
  TEST("volume is now 4", b2.volume(), 4.0);
  b2.scaleAboutCentroid(0.5);
  TEST("scaleAboutCentroid", b2.centroid(), p0);
  TEST("volume is now 1", b2.volume(), 1.0);
  b2.expandAboutCentroid(1.0);
  TEST("expandAboutCentroid", b2.centroid(), p0);
  TEST("volume is now 4", b2.volume(), 4.0);
  b2.expandAboutCentroid(-1.0);
  TEST("expandAboutCentroid", b2.centroid(), p0);
  TEST("volume is now 1", b2.volume(), 1.0);
  b2.setCentroid(p12);
  b2.setMinPoint(p0);
  TEST("setMinPoint", b2.volume(), 2.25);
  b2.setMaxPoint(b.maxPoint());
  TEST("setMaxPoint", b2, b);

  float d0[2] = {0.0,0.0}, d1[2] = {1.0,1.0}, d2[2] = {2.0,2.0};
  b = CAABBox2D(d0,d2); 
  Debug::debug(Debug::math,__FUNCTION__) << b << endl;
  TEST("construct from two points", b, b2);
  b = CAABBox2D(d1,2,2,CAABBox2D::centre); 
  Debug::debug(Debug::math,__FUNCTION__) << b << endl;
  TEST("construct from centroid", b, b2);
  b = CAABBox2D(p12,2,2,CAABBox2D::centre);
  Debug::debug(Debug::math,__FUNCTION__) << b << endl;
  TEST("construct from centroid", b, b2);
  p12 = b2.minPoint();
  b = CAABBox2D(p12,2,2,CAABBox2D::min_pos);
  Debug::debug(Debug::math,__FUNCTION__) << b << endl;
  TEST("construct from min_pos", b, b2);
  p12 = b2.maxPoint();
  b = CAABBox2D(p12,2,2,CAABBox2D::max_pos);
  Debug::debug(Debug::math,__FUNCTION__) << b << endl;
  TEST("construct from max_pos", b, b2);

  b.empty(); 
  Debug::debug(Debug::math,__FUNCTION__) << b << endl;
  TEST("empty", b2==b, false);
  TEST("!contains(p0)", b.contains(p0), false);
  TEST("!contains(p1)", b.contains(p1), false);
  TEST("!contains(p2)", b.contains(p2), false);
  TEST("AABBox2D is empty", b.isEmpty(), true);
  TEST("AABBox2D has no volume", b.volume(), 0.0);

  CAABBox2D ib(10, 11, 10, 11);
  Debug::debug(Debug::math,__FUNCTION__) << ib << endl;
  TEST("Integer AABBox2D centroid", ib.getCentroidX() == 10.5 && ib.getCentroidY()==10.5, true);
  ib.setWidth(1); ib.setHeight(1);
  Debug::debug(Debug::math,__FUNCTION__) << ib << endl;
  TEST("Integer AABBox2D centroid drift", ib.getCentroidX() == 10.5 && ib.getCentroidY()==10.5, true);

  ib = CAABBox2D(10, 11, 10, 11);
  ib.setWidth(4); ib.setHeight(4);
  ib.setWidth(3); ib.setHeight(3);
  ib.setWidth(4); ib.setHeight(4);
  ib.setWidth(3); ib.setHeight(3);
  ib.setWidth(4); ib.setHeight(4);
  Debug::debug(Debug::math,__FUNCTION__) << ib << ib.centroid() << endl;
  TEST("Integer AABBox2D centroid drift", (int)ib.getCentroidX() == 10 && (int)ib.getCentroidY()==10, true);

  ib = CAABBox2D(9, 11, 9, 11);
  Debug::debug(Debug::math,__FUNCTION__) << ib << ib.centroid() << endl;
  ib.setWidth(3); ib.setHeight(3);
  Debug::debug(Debug::math,__FUNCTION__) << ib << ib.centroid() << endl;
  TEST("Integer AABBox2D centroid", ib.getCentroidX() == 10 && ib.getCentroidY()==10, true);

  ib = CAABBox2D(-11, -10, -11, -10);
  Debug::debug(Debug::math,__FUNCTION__) << ib << ib.centroid() << endl;
  ib.setWidth(3); ib.setHeight(3);
  ib.setWidth(4); ib.setHeight(4);
  ib.setWidth(3); ib.setHeight(3);
  ib.setWidth(4); ib.setHeight(4);
  ib.setWidth(3); ib.setHeight(3);
  Debug::debug(Debug::math,__FUNCTION__) << ib << ib.centroid() << endl;
  TEST("Integer AABBox2D negative centroid drift", (int)ib.getCentroidX() == -10 && (int)ib.getCentroidY()==-10, true);

  ib = CAABBox2D(-11, -9, -11, -9);
  Debug::debug(Debug::math,__FUNCTION__) << ib << ib.centroid() << endl;
  ib.setWidth(3); ib.setHeight(3);
  ib.setWidth(4); ib.setHeight(4);
  ib.setWidth(3); ib.setHeight(3);
  ib.setWidth(4); ib.setHeight(4);
  Debug::debug(Debug::math,__FUNCTION__) << ib << ib.centroid() << endl;
  TEST("Integer AABBox2D negative centroid drift", ib.getCentroidX() == -10 && ib.getCentroidY()==-10, true);

  CPoint2D min1(10,10), max1(20,20),
                       min2(40,40), max2(50,50),
                       min3(45,45), max3(55,55);

  CAABBox2D box1(min1, max1);
  CAABBox2D box2(min2, max2);
  CAABBox2D box3(min3, max3);

  //no intersection case
  CAABBox2D i1 = box1.intersection( box2);
  TEST("box1.intersection( box2) = false", true, i1.isEmpty());

  //intersection case
  CAABBox2D i2 = box2.intersection( box3);
  TEST("box2.intersection( box3) = true", false, i2.isEmpty());
  TEST("box2.intersection( box3) volume", 25.0, i2.volume());

  stringstream is; is << "4.4 -5 7e1 5e-1";
  CAABBox2D l; is >> l;
  TEST("istream CAABBox2D", l, CAABBox2D(4.4,70,-5,0.5)); // note different order!!
}




} //end GEO

namespace MATH {

void CVec2D::test_vector_2d()
{
	using GEO::_2D::operator<<;
	using GEO::_2D::operator>>;
  // standard constructor
  CVec2D  v(1.5f, 0.625f);
  Debug::debug(Debug::math,__FUNCTION__) << v << endl;

  // default constructor
  CVec2D v0; // == (0,0)

  // comparison
  TEST("inequality", (v != v0), true);
  // length
  TEST("length", v0.getLenght(), 0.0);

  v0.set(1.5, 0.625);
  TEST("equality", (v == v0), true);
  v0 = v; // assignment
  TEST("equality", (v == v0), true);

  TEST("length", v.getLenght(), (float)1.62500012);
  // should be "exact" ! (all these numbers are exactly representable in base 2)

  v0.set(1.5, 0);
  TEST("inequality", (v != v0), true);
  CVec2D v1 (3,0.625);
  TEST("sum", (v+v0), v1);
  TEST("difference", (v1-v0), v);
  TEST("scale", (-v1+2*v0).x, 0.0f);
  TEST("scale", (-v1+2*v0).y, -0.625f);

  TEST("dot_product", (v*v1), 4.890625);
  TEST("cross_product", v1.cross(v), 0.625*1.5);
  float angle=v1.getAngle(v);
  TEST_NEAR("angle", angle, 0.18939573, 1e-8); // 10^51'06"

  TEST("parallel", (v.isParallel(v)), true);
  TEST("parallel", (v.isParallel(CVec2D())), true); // parallel to (0,0)
  TEST("parallel", (v.isParallel(v1,0.1)), false); // not parallel, even with tol=0.1
  
//  TEST("ratio", v/v, 1);
//  TEST("ratio", (-v*3.5)/v, -3.5);

//  TEST_NEAR("ratio", v1/v, 1.852071, 1e-6);
  
  //TEST_NEAR("normalized", length(length(v1)*normalized(v1) -v1), 0.0, 1e-6);
  //CVec2D test1,test2,test3;
  //test1=v1.normalized();
  //float lenght_v1=v1.getLenght();
  //test2=lenght_v1*test1;
  //test3=test2-v1;
  //float lenght_test3=test3.getLenght();

  TEST_NEAR("normalized", ((v1.getLenght()*v1.normalized()) -v1).getLenght(), 0.0, 1e-4);
  v0=v1;
  v1.toNormal();
  //TEST_NEAR("normalize", length(length(v0)*v1 - v0), 0.0, 1e-6);
  TEST_NEAR("normalize", (v0.getLenght()*v1 - v0).getLenght(), 0.0, 1e-6);

  TEST("orthogonal", CVec2D::areOrthogonal(v,CVec2D()), true); // orthogonal to (0,0)
  TEST("!orthogonal", CVec2D::areOrthogonal(v,v1,0.1), false); // even not with tolorance
  TEST("orthogonal", CVec2D::areOrthogonal(v,CVec2D(0.625,-1.5)), true);

  {
    stringstream is; is << "4.4 -5 7e1";
    CVec2D p(0.0,0.0); is >> p;
    Debug::debug(Debug::math,__FUNCTION__) << p << endl;
    TEST("istream CVec2D (blank-separated)", p, CVec2D(4.4,-5));
  }
  {
    stringstream is; is << "7e1, 11 , blabla";
    CVec2D p(0.0,0.0); is >> p;
    Debug::debug(Debug::math,__FUNCTION__) << p << endl;
    TEST("istream CVec2D (comma-separated)", p, CVec2D(70,11));
  }
  {
    stringstream is; is << " (12,13 ) !";
    CVec2D p(0.0,0.0); is >> p;
    Debug::debug(Debug::math,__FUNCTION__) << p << endl;
    TEST("istream CVec2D (parenthesized)", p, CVec2D(12,13));
  }

}



}
}  //end SMF
