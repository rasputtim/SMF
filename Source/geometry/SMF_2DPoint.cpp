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
#include "geometry/SMF_2DLineSegment.h"
#include "geometry/SMF_GeoDefs.h"
#include "math/all.h"
#include "geometry/all.h"
#include <vector>

namespace SMF {
using namespace MATH;
using namespace Util;

namespace GEO{

float CPoint2D::distanceTolinesegmentSqr(float x1, float y1,
                                    float x2, float y2,
                                    float x, float y)
{
  // squared distance between endpoints :
  float ddh = CMath::square(x2-x1) + CMath::square(y2-y1);

  // squared distance to endpoints :
  float dd1 = CMath::square(x-x1) + CMath::square(y-y1);
  float dd2 = CMath::square(x-x2) + CMath::square(y-y2);

  // if closest to the start point :
  if (dd2 >= ddh + dd1)
    return dd1;

  // if closest to the end point :
  if (dd1 >= ddh + dd2)
    return dd2;

  // squared perpendicular distance to line :
  float a = y1-y2;
  float b = x2-x1;
  float c = x1*y2-x2*y1;
  return CMath::square(a*x + b*y + c)/float(a*a + b*b);
}
float CPoint2D::distance(CLineSegment2D const& l)
{
	CPoint2D const& p=*this;
	return CMath::sqrt(CLineSegment2D::distanceSegmentPoint(l, p));
}

float CPoint2D::distance( CLine2D const& l)
{ CPoint2D const& p=*this;
  float num = l.a()*p.x() + l.b()*p.y() + l.c();
  if (num == 0) return 0.0; // no call to sqrt if not necessary
  else return CMath::fabs(num) / CMath::sqrt(l.a()*l.a() + l.b()*l.b());
}


float CPoint2D::distance(CPoint2D const &p)const{
	float dist= CMath::sqrt(CMath::square(x()-p.x())+CMath::square(y()-p.y()));
	return dist;
}

float CPoint2D::distance(CPolygon2D const &poly, bool closed)const{
return distanceToPolygon(poly,x(),y(),closed);
}

float CPoint2D::distanceSqr(CPoint2D const &p)const{
	float dist= CMath::square(x()-p.x())+CMath::square(y()-p.y());
	return dist;
}

float CPoint2D::distanceToPolygon(CPolygon2D const&poly, float x, float y, bool closed)
{
	float dd;
	int n=poly.numVertices();
	if(closed){
	 dd = distanceTolinesegmentSqr(poly[n-1].x(), poly[n-1].y(),
											  poly[0].x(), poly[0].y(),
											  x, y);
	  for (unsigned i=0; i+1<n; ++i) {
		float nd = distanceTolinesegmentSqr(poly[i].x(), poly[i].y(),
												poly[i+1].x(), poly[i+1].y(),
												x, y);
		if (nd<dd)
		  dd = nd;
	  }

	}else {
		  dd = -1;
		  for (unsigned i=0; i+1<n; ++i) {
			float nd = distanceTolinesegmentSqr(poly[i].x(), poly[i].y(),
													poly[i+1].x(), poly[i+1].y(),
													x, y);
			if (dd<0 || nd<dd)
			  dd = nd;
		  }


	}
  return CMath::sqrt(dd);
}
//: Read from stream, possibly with formatting
//  Either just reads two blank-separated numbers,
//  or reads two comma-separated numbers,
//  or reads two numbers in parenthesized form "(123, 321)"
istream& CPoint2D::read(istream& is)
{
  if (! is.good()) return is; // (TODO: should throw an exception)
  bool paren = false;
  float tx, ty;
  is >> std::ws; // jump over any leading whitespace
  if (is.eof()) return is; // nothing to be set because of EOF (TODO: should throw an exception)
  if (is.peek() == '(') { is.ignore(); paren=true; }
  is >> std::ws >> tx >> std::ws;
  if (is.eof()) return is;
  if (is.peek() == ',') is.ignore();
  is >> std::ws >> ty >> std::ws;
  if (paren) {
    if (is.eof()) return is;
    if (is.peek() == ')') is.ignore();
    else                  return is; // closing parenthesis is missing (TODO: throw an exception)
  }
  set(tx,ty);
  return is;
}

void CPoint2D::closestPointToLinesegment(float& ret_x, float& ret_y,
                                      float x1, float y1,
                                      float x2, float y2,
                                      float x0, float y0)
{
  // squared distance between endpoints:
  float ddh = CMath::square(x2-x1) + CMath::square(y2-y1);

  // squared distance to endpoints:
  float dd1 = CMath::square(x0-x1) + CMath::square(y0-y1);
  float dd2 = CMath::square(x0-x2) + CMath::square(y0-y2);

  // if closest to the start point:
  if (dd2 > ddh + dd1) { ret_x=x1; ret_y=y1; return; }

  // if closest to the end point :
  if (dd1 > ddh + dd2) { ret_x=x2; ret_y=y2; return; }

  // line through (x0,y0) and perpendicular to the given line is
  // the line with equation (x-x0)(x2-x1)+(y-y0)(y2-y1)=0.
  // Then it just remains to intersect these two lines:
  float dx = x2-x1;
  float dy = y2-y1;
  float c = dx*dx+dy*dy;
  ret_x = float((dx*dx*x0+dy*dy*x1-dx*dy*(y1-y0))/c); // possible rounding error!
  ret_y = float((dx*dx*y1+dy*dy*y0-dx*dy*(x1-x0))/c);
}
CPoint2D CPoint2D::closestPoint( CLineSegment2D const &l1)const 
{   
	float x0 = x();
	float y0 = y();
	float x1=l1.point1().x(); float y1=l1.point1().y();

    float x2=l1.point2().x();float y2=l1.point2().y();
  // squared distance between endpoints:
  float ddh = CMath::square(x2-x1) + CMath::square(y2-y1);

  // squared distance to endpoints:
  float dd1 = CMath::square(x0-x1) + CMath::square(y0-y1);
  float dd2 = CMath::square(x0-x2) + CMath::square(y0-y2);

  // if closest to the start point:
  if (dd2 > ddh + dd1) {  return CPoint2D(x1,y1); }

  // if closest to the end point :
  if (dd1 > ddh + dd2) {  return CPoint2D(x2,y2); }

  // line through (x0,y0) and perpendicular to the given line is
  // the line with equation (x-x0)(x2-x1)+(y-y0)(y2-y1)=0.
  // Then it just remains to intersect these two lines:
  float dx = x2-x1;
  float dy = y2-y1;
  float c = dx*dx+dy*dy;
  float ret_x = float((dx*dx*x0+dy*dy*x1-dx*dy*(y1-y0))/c); // possible rounding error!
  float ret_y = float((dx*dx*y1+dy*dy*y0-dx*dy*(x1-x0))/c);
  return CPoint2D(ret_x,ret_y);
}

CPoint2D CPoint2D::closestPoint( CRay2D const &l1)const
{   
CLine2D line = l1.toLine2D();
CPoint2D result = closestPoint(line);
if (l1.contains(result)) return result;
else return l1.pos;
}

CPoint2D CPoint2D::closestPoint( CLine2D const &l)const
{   
  CPoint2D const& p = *this;
  float d = l.a()*l.a()+l.b()*l.b();
  SMF_ASSERT(d!=0); // line should not be the line at infinity
  return CPoint2D((l.b()*l.b()*p.x()-l.a()*l.b()*p.y()-l.a()*l.c())/d,
                         (l.a()*l.a()*p.y()-l.a()*l.b()*p.x()-l.b()*l.c())/d);

}
static bool IsEqual(const float Val1,const float Val2,float Epsilon){

 float  Diff;
 bool Result;
  Diff = Val1 - Val2;
  SMF_ASSERT(((-Epsilon <= Diff) && (Diff <= Epsilon)) == (CMath::fabs(Diff) <= Epsilon));
  Result = ((-Epsilon <= Diff) && (Diff <= Epsilon));
  return Result;
}
 
bool CPoint2D::IsCollinear(CPoint2D const &p2,CPoint2D const &p3, float tolerance)const{
	CPoint2D p1=*this;
	bool Result = IsEqual((x_*p2.y_) + (y_ *p3.x_)+(p2.x_*p3.y_) ,(p2.y_*p3.x_) +(p3.y_*x_)+ (p2.x_*y_),tolerance);
	return Result;
	//return (p1-p2).isParallel(p1-p3); //this way of solving  caused wrong results. do not know why
}
CPoint2D CPoint2D::transform(CMat3D const &transMat)const{
   CPoint2D TempPoint;
   float a,b,c,d,m,n;
   a = transMat[0][0];
   b = transMat[0][1];
   c = transMat[1][0];
   d = transMat[1][1];
   m = transMat[2][0];
   n = transMat[2][1];

   TempPoint.x_ = x_*a + y_*c + m;
   TempPoint.y_ = x_*b + y_*d + n;
   return TempPoint;

}

int CPoint2D::orientation(const CPoint2D &px1,const CPoint2D &px2,const CPoint2D &  point){
  
  float Orin;
  int Result;
  /* Determinant of the 3 points */
  Orin = (px2.x_ - px1.x_) * (point.y_ - px1.y_) - (point.x_ - px1.x_) * (px2.y_ - px1.y_);

  if (Orin > CMath::ZERO){
    Result = LeftHandSide;          /* Orientaion is to the left-hand side  */
  }else if (Orin < CMath::ZERO){
    Result = RightHandSide;         /* Orientaion is to the right-hand side */
  }else{
	Result = CollinearOrientation; /* Orientaion is neutral aka collinear  */
	}
	return Result;
}

namespace _2D{

//  +-+-+ point_2d geometry +-+-+

float cross_ratio(CPoint2D const& p1, CPoint2D const& p2,
                   CPoint2D const& p3, CPoint2D const& p4){
  // least squares solution: (Num_x-CR*Den_x)^2 + (Num_y-CR*Den_y)^2 minimal.
  float Num_x = (p1.x()-p3.x())*(p2.x()-p4.x());
  float Num_y = (p1.y()-p3.y())*(p2.y()-p4.y());
  float Den_x = (p1.x()-p4.x())*(p2.x()-p3.x());
  float Den_y = (p1.y()-p4.y())*(p2.y()-p3.y());
  if (Den_x == Den_y) return 0.5*(Num_x+Num_y)/Den_x;
  else return (Den_x*Num_x+Den_y*Num_y)/(Den_x*Den_x+Den_y*Den_y);
}

//: Are three points collinear, i.e., do they lie on a common line?
// \see  CPoint2D
bool collinear(const CPoint2D & p1,
               const CPoint2D & p2,
               const CPoint2D & p3)
{ 	return (p1-p2).isParallel(p1-p3); }

#if 0
//: Return the relative distance to p1 wrt p1-p2 of p3.
//  The three points should be collinear and p2 should not equal p1.
//  This is the coordinate of p3 in the affine 1D reference frame (p1,p2).
//  If p3=p1, the ratio is 0; if p1=p3, the ratio is 1.
//  The mid point of p1 and p2 has ratio 0.5.
//  Note that the return type is float, not float, since the ratio of e.g.
//  two CVec2D<int> need not be an int.
// \see  CPoint2D
float ratio(CPoint2D const& p1,
             CPoint2D const& p2,
             CPoint2D const& p3)
{
	return (p3-p1)/(p2-p1); }
#endif
//: Return the point at a given ratio wrt two other points.
//  By default, the mid point (ratio=0.5) is returned.
//  Note that the third argument is float, not float, so the midpoint of e.g.
//  two CPoint2D<int> is not a valid concept.  But the reflection point
//  of p2 wrt p1 is: in that case f=-1.
// \see  CPoint2D
CPoint2D midpoint(CPoint2D const& p1,
                            CPoint2D const& p2,
                            float f)
{
  return CPoint2D ((float)((1-f)*p1.x() + f*p2.x()),
                            (float)((1-f)*p1.y() + f*p2.y()));
}


//: Return the point at the centre of gravity of two given points.
// Identical to midpoint(p1,p2).
// \see  CPoint2D
CPoint2D centre(CPoint2D const& p1,
                          CPoint2D const& p2)
{
  return CPoint2D ((p1.x() + p2.x())/2 ,
                            (p1.y() + p2.y())/2 );
}

//: Return the point at the centre of gravity of three given points.
// \see  CPoint2D
CPoint2D centre(CPoint2D const& p1,
                          CPoint2D const& p2,
                          CPoint2D const& p3)
{
  return CPoint2D ((p1.x() + p2.x() + p3.x())/3 ,
                            (p1.y() + p2.y() + p3.y())/3 );
}

//: Return the point at the centre of gravity of four given points.
// \see  CPoint2D
CPoint2D centre(CPoint2D const& p1,
                          CPoint2D const& p2,
                          CPoint2D const& p3,
                          CPoint2D const& p4)
{
  return CPoint2D ((p1.x() + p2.x() + p3.x() + p4.x())/4 ,
                            (p1.y() + p2.y() + p3.y() + p4.y())/4 );
}

//: Return the point at the centre of gravity of a set of given points.
// Beware of possible rounding errors when float is e.g. int.
// \see  CPoint2D
CPoint2D centre(std::vector<CPoint2D> const& v)
{
  int n=v.size();
  SMF_ASSERT(n>0); // it is *not* correct to return the point (0,0) when n==0.
  float x = 0, y = 0;
  for (int i=0; i<n; ++i) x+=v[i].x(), y+=v[i].y();
  return CPoint2D (x/n,y/n);
};


//: Write "<CPoint2D x,y> " to stream
ostream&  operator<<(ostream& s, CPoint2D const& p)
{
  return s << "<CPoint2D "<< p.x() << ',' << p.y() << "> ";
}

//: Read x y from stream
istream&  operator>>(istream& is,  CPoint2D& p)
{
  return p.read(is);
}

//: Write "<CVec2D x,y> " to stream
ostream&  operator<<(ostream& s, CVec2D const& p)
{
  return s << "<CVec2D "<< p.x << ',' << p.y <<  "> ";
}


//: Read x y from stream
istream&  operator>>(istream& is, CVec2D& p)
{
  return p.read(is);
}

} //end 2D
} //end GEO
}  //end SMF
