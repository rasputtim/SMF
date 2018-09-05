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
#include "math/all.h"
#include "geometry/all.h"
#include <vector>

namespace SMF {
using namespace MATH;
using namespace Util;

namespace GEO{

	//: line through two given points
CLine2D::CLine2D (CPoint2D const& p1, CPoint2D const& p2)
: a_ ( p1.y() - p2.y() )
, b_ ( p2.x() - p1.x() )
, c_ ( p1.x() * p2.y() - p1.y() * p2.x() )
{
  SMF_ASSERT(a_||b_); // two points were distinct
}

//: line defined by one point and one vector

CLine2D::CLine2D (CPoint2D const& p, CVec2D const& v)
: a_ ( -v.y )
, b_ ( v.x )
, c_ ( -a_*p.x() - b_*p.y() )  //p.x_*(p.y_+v.y)-p.y_*(p.x_+v.x)
{
}

#if 0
CLine2D::CLine2D (CPoint2DHom const& l)
 : a_(l.a()) , b_(l.b()) , c_(l.c())
{
  //JLM I see no reason to prohibit lines through the origin
  //  SMF_ASSERT(c_);
}
#endif
//: Get two points on the line.
// These two points are normally the intersections
// with the Y axis and X axis, respectively.  When the line is parallel to one
// of these, the point with \a y=1 or \a x=1, resp. are taken.
// When the line goes through the origin, the second point is \a (b,-a).

void CLine2D::get_two_points(CPoint2D &p1, CPoint2D &p2) const
{
  if (b() == 0)       p1.set(-c()/a(),(float) 1);
  else                p1.set((float)0, -c()/b());
  if (a() == 0)       p2.set((float)1, -c()/b());
  else if ( c() == 0) p2.set(b(), -a());
  else                p2.set(-c()/a(),(float) 0);
}


CPoint2D CLine2D::getPoint(float x)const{
//	-c-ax /b
 if(b() ==0) return CPoint2D(x,(CMath::NAN_FLOAT));
 else return CPoint2D(x,((-c()-(a()*x))/b()));
}
float CLine2D::getY(float const &x)const{
	float y;
	if (b() !=0) {
	y = (-c()-(a()*x))/b();
	return y;
	}else return CMath::NAN_FLOAT;
}
float CLine2D::getX( const double & y ) const
{
     if (a() !=0) {
		 return -( b() * y + c() ) / a();
	 }else return CMath::NAN_FLOAT;
}
double CLine2D::slope_degrees() const
{
	static const double deg_per_rad = 45.0/CMath::atan64(1.0,1.0);
  // do special cases separately, to avoid rounding errors:
  if (a() == 0) return b()<0 ? 0.0 : 180.0;
  if (b() == 0) return a()<0 ? -90.0 : 90.0;
  if (a() == b()) return a()<0 ? -45.0 : 135.0;
  if (a()+b() == 0) return a()<0 ? -135.0 : 45.0;
  // general case:
  return deg_per_rad * CMath::atan64(double(a()),-double(b()));
}

CVec2D CLine2D::directionVector()const{
	CPoint2D p1,p2;
	this->get_two_points(p1,p2);
	CVec2D posVec1(p1);
	CVec2D posVec2(p2);
	return (posVec2-posVec1).normalized();
}
CPoint2D CLine2D::closestPoint(CPoint2D const  &p)const{
  return p.closestPoint(*this);
}

double CLine2D::slope_radians() const
{
	return CMath::atan64(double(a()),-double(b()));
}

float CLine2D::angular_coef() const
{

	if (b() ==0) {
		return CMath::NAN_FLOAT;
	}else return -a()/b();
}

bool CLine2D::isParallel(CLine2D const &l,  float maxDistance)const{
	float coef1=this->angular_coef();
	float coef2=l.angular_coef();
		if (MATH::isFinite(coef1)&& MATH::isFinite(coef1)){
			return this->angular_coef() - l.angular_coef() <= maxDistance;
		}else if (MATH::isFinite(coef1)&& !MATH::isFinite(coef1)){
			return false;
		}else if (!MATH::isFinite(coef1)&& MATH::isFinite(coef1)){
			return false;
		}else if (!MATH::isFinite(coef1)&& !MATH::isFinite(coef1)){
			return true;
		}
}

bool CLine2D::isPerpendicular(CLine2D const &l,  float maxDistance)const{
	float coef1=this->angular_coef();
	float coef2=l.angular_coef();
	bool result;
	if (MATH::isFinite(coef1)&& MATH::isFinite(coef1)){
				float teste = coef1*coef2 +1;
				result = teste <= maxDistance;
				return result;

		}else if (MATH::isFinite(coef1)&& !MATH::isFinite(coef1)){
			if (coef2==0) return true;
			else return false;
		}else if (!MATH::isFinite(coef1)&& MATH::isFinite(coef1)){
			if (coef1==0) return true;
			else return false;
		}else if (!MATH::isFinite(coef1)&& !MATH::isFinite(coef1)){
			return false;
		}

}

bool CLine2D::normalize()
{

  double sqrLength, invLength;
    sqrLength = a_*a_ + b_*b_;
	if (sqrLength==1.0) return true;
    if (sqrLength==0.0) return false;
    // inverse square root with 32 bits precision, returns huge number when a == 0.0
	invLength = CMath::invSqrt( sqrLength );
	a_ *= invLength;
	b_ *= invLength;
	c_ *= invLength;
	sqrLength = a_*a_ + b_*b_;
   // return false when normalisation did not succeed, e.g. when float == int:
   return sqrLength>0.99 && sqrLength<1.01;
}

double CLine2D::distance( CPoint2D const& p) const
{ CLine2D const& l=*this;
  float num = l.a()*p.x() + l.b()*p.y() + l.c();
  if (num == 0) return 0.0; // no call to sqrt if not necessary
  else return CMath::fabs(num) / CMath::sqrt(l.a()*l.a() + l.b()*l.b());
}
bool CLine2D::intersects(CPoint2D const &p)const{
	float result = (p.x()*a())+ (p.y()*b())+c();
	return result == 0;
}

CLineSegment2D CLine2D::toLineSegment2D() const{
	CPoint2D p1;
	CPoint2D p2;
	get_two_points(p1, p2);
	CLineSegment2D lx(p1,p2);
	return lx;
}

bool CLine2D::intersects(CRay2D const &l)const{
	CPoint2D p1;
	CPoint2D p2;
	CPoint2D p3;
	CPoint2D p4;
	CPoint2D intersection;
	CLine2D ly = l.toLine2D();

	bool result= intersects(*this,ly,intersection);
	if (result == true){
		bool pertence = l.intersects(intersection);
		if (pertence == true) return true;
		else return false;
	}else return false;
	
}

bool CLine2D::intersects(CLineSegment2D const &l)const{
	
	CPoint2D intersection;
	CLine2D lteste = l.toLine2D();
	bool result= intersects(*this,lteste,intersection);
	//verifica se intersection pertence ao segmento
	if(result) {
	bool result2 = l.intersects(intersection);
	return result2;
	}else return false;

}


bool CLine2D::intersects( const CLine2D &line0,
                       const CLine2D &line1,
                       CPoint2D      &intersection_point , float maxDistance)
{
  float a0, b0, c0,  a1, b1, c1;
  a0 = line0.a(); b0 = line0.b(); c0 = line0.c();
  a1 = line1.a(); b1 = line1.b(); c1 = line1.c();

  float delta, delta_x, delta_y, x, y;
  delta = a0*b1 - a1*b0;
  if ( CMath::fabs(delta) <= maxDistance ) // Lines are parallel
    return false;
  delta_x = -c0*b1 + b0*c1; delta_y = -a0*c1 + a1*c0;
  x = delta_x / delta; y = delta_y / delta;

  //   intersection_point.set( (float)x, (float)y );
  intersection_point.set( x, y );
  return true;
}

bool CLine2D::intersects(const CAABBox2D& box)const
{
	return intersection(box);
}
//: compute the intersection of an infinite line with *this box.
//  p0 and p1 are the intersection points
// In the normal case (no degeneracies) there are six possible intersection combinations:
// \verbatim
//
//                C01 /    CY     \ C11
//                   /     |       \           .
//       ymax  -----/------|--------\-----
//            |    /       |         \    |
//            |   /        |          \   |
//            |  /         |           \  | \  .
//            | /          |            \ |  \_ Bounding Box
//            |/           |             \|
//            /            |              \    .
//           /|            |              |\   .
//           ---------------------------------- CX
//          \ |            |              /
//           \|            |             /|
//            \            |            / |
//            |\           |           /  |
//            | \          |          /   |
//            |  \         |         /    |
//       xmin  ---\--------|--------/-----   xmax
//       ymin      \       |       /
//              C00 \             / C10
// \endverbatim

bool CLine2D::intersection(const CAABBox2D& box, int *ICnt,
                      CPoint2D* p00,
                      CPoint2D* p01)const
{
#define USE_IT(x,val) if(x!=NULL) { *x=val;}
	if (p00!=NULL) p00->toInfinite();
	if (p01!=NULL) p01->toInfinite();
	if (ICnt!=NULL) *ICnt=0;

	const CLine2D& line=*this;
                      
  double a = line.a(), b = line.b(), c = line.c();
  double xmin=box.minX(), xmax=box.maxX();
  double ymin=box.minY(), ymax=box.maxY();
  CPoint2D p0,p1;
  int count=0;
  // Run through the cases
  //
  if (CMath::nearZero(a))// The line is y = -c/b
  {
    float y0 = static_cast<float>(-c/b);
    // The box edge is collinear with line?
    if (CMath::equals(ymin,y0))
    {
      p0.set(static_cast<float>(xmin), static_cast<float>(ymin));
      p1.set(static_cast<float>(xmax), static_cast<float>(ymin));
	  count=2;
	  USE_IT(p00,p0)
	  USE_IT(p01,p1)
	  USE_IT(ICnt,count)
      return true;
    }
    if (CMath::equals(ymax,y0))
    {
      p0.set(static_cast<float>(xmin), static_cast<float>(ymax));
      p1.set(static_cast<float>(xmax), static_cast<float>(ymax));
      count=2;
	  USE_IT(p00,p0)
	  USE_IT(p01,p1)
	  USE_IT(ICnt,count)
	  return true;
    }

    if ((ymin > y0) || (y0 > ymax)){ // The line does not intersect the box
        count=0;
	    USE_IT(ICnt,count)
		return false;
	}else // The line does intersect
    {
      p0.set(static_cast<float>(xmin), static_cast<float>(y0));
      p1.set(static_cast<float>(xmax), static_cast<float>(y0));
	  count=2;
	  USE_IT(p00,p0)
	  USE_IT(p01,p1)
	  USE_IT(ICnt,count)
	  return true;
    }
  }

  if (CMath::nearZero(b))// The line is x = -c/a
  {
    float x0 = static_cast<float>(-c/a);
    // The box edge is collinear with l?
    if (CMath::equals(xmin,x0))
    {
      p0.set(static_cast<float>(xmin), static_cast<float>(ymin));
      p1.set(static_cast<float>(xmin), static_cast<float>(ymax));
      count=2;
	  USE_IT(p00,p0)
	  USE_IT(p01,p1)
	  USE_IT(ICnt,count)
	  return true;
    }
    if (CMath::equals(xmax,x0))
    {
      p0.set(static_cast<float>(xmax), static_cast<float>(ymin));
      p1.set(static_cast<float>(xmax), static_cast<float>(ymax));
      count=2;
	  USE_IT(p00,p0)
	  USE_IT(p01,p1)
	  USE_IT(ICnt,count)
	  return true;
    }

    if (xmin <= x0 && x0 <= xmax) // The line intersects the box
    {
      p0.set(static_cast<float>(x0), static_cast<float>(ymin));
      p1.set(static_cast<float>(x0), static_cast<float>(ymax));
      count=2;
	  USE_IT(p00,p0)
	  USE_IT(p01,p1)
	  USE_IT(ICnt,count)
	  return true;
    }
    else{
      count=0;
	  USE_IT(ICnt,count)
	  return false;

	}
  }

  // The normal case with no degeneracies
  //
  // intersection with x = xmin
  float y_xmin_int = static_cast<float>(-(c + a*xmin)/b);
  bool inside_xmin = (y_xmin_int >= ymin) && (y_xmin_int <= ymax);

  // intersection with x = xmax
  float y_xmax_int = static_cast<float>(-(c + a*xmax)/b);
  bool inside_xmax = (y_xmax_int >= ymin) && (y_xmax_int <= ymax);

  // intersection with y = ymin
  float x_ymin_int = static_cast<float>(-(c + b*ymin)/a);
  bool inside_ymin = (x_ymin_int >= xmin) && (x_ymin_int <= xmax);

  // intersection with y = ymax
  float x_ymax_int = static_cast<float>(-(c + b*ymax)/a);
  bool inside_ymax = (x_ymax_int >= xmin) && (x_ymax_int <= xmax);

  // Case CX
  if (inside_xmin && inside_xmax &&
      !(CMath::equals(y_xmin_int,ymin) && CMath::equals(y_xmax_int,ymax)))
  {
    p0.set(static_cast<float>(xmin), static_cast<float>(y_xmin_int));
    p1.set(static_cast<float>(xmax), static_cast<float>(y_xmax_int));
    count=2;
	USE_IT(p00,p0)
	USE_IT(p01,p1)
	USE_IT(ICnt,count)
	return true;
  }

  // Case CY
  if (inside_ymin && inside_ymax &&
      !(CMath::equals(x_ymin_int,xmin) && CMath::equals(x_ymax_int,xmax)))
  {
    p0.set(static_cast<float>(x_ymin_int), static_cast<float>(ymin));
    p1.set(static_cast<float>(x_ymax_int), static_cast<float>(ymax));
    count=2;
	USE_IT(p00,p0)
	USE_IT(p01,p1)
	USE_IT(ICnt,count)
	return true;
  }

  // Case C00
  if (inside_xmin && inside_ymin &&
      !(inside_xmax && inside_ymax))
  {
    p0.set(static_cast<float>(xmin), static_cast<float>(y_xmin_int));
    p1.set(static_cast<float>(x_ymin_int), static_cast<float>(ymin));
    count=2;
	USE_IT(p00,p0)
	USE_IT(p01,p1)
	USE_IT(ICnt,count)
	return true;
  }

  // Case C01
  if (inside_xmin && inside_ymax &&
      !(inside_xmax && inside_ymin))
  {
    p0.set(static_cast<float>(xmin), static_cast<float>(y_xmin_int));
    p1.set(static_cast<float>(x_ymax_int), static_cast<float>(ymax));
    count=2;
	USE_IT(p00,p0)
	USE_IT(p01,p1)
	USE_IT(ICnt,count)
	return true;
  }

  // Case C10
  if (inside_ymin && inside_xmax &&
      !(inside_xmin && inside_ymax))
  {
    p0.set(static_cast<float>(x_ymin_int), static_cast<float>(ymin));
    p1.set(static_cast<float>(xmax), static_cast<float>(y_xmax_int));
    count=2;
	  USE_IT(p00,p0)
	  USE_IT(p01,p1)
	  USE_IT(ICnt,count)
	return true;
  }

  // Case C11
  if (inside_ymax && inside_xmax &&
      !(inside_xmin && inside_ymin))
  {
    p0.set(static_cast<float>(x_ymax_int), static_cast<float>(ymax));
    p1.set(static_cast<float>(xmax), static_cast<float>(y_xmax_int));
    count=2;
	  USE_IT(p00,p0)
	  USE_IT(p01,p1)
	  USE_IT(ICnt,count)
	return true;
  }
  // Exactly passing through diagonal of BB
  if (inside_xmin && inside_xmax && inside_ymin && inside_ymax)
  {
    if (a>0) // 45 degrees
    {
      p0.set(static_cast<float>(xmin), static_cast<float>(ymin));
      p1.set(static_cast<float>(xmax), static_cast<float>(ymax));
      count=2;
	  USE_IT(p00,p0)
	  USE_IT(p01,p1)
	  USE_IT(ICnt,count)
	  return true;
    }
    else // 135 degrees
    {
      p0.set(static_cast<float>(xmin), static_cast<float>(ymax));
      p1.set(static_cast<float>(xmax), static_cast<float>(ymin));
      count=2;
	  USE_IT(p00,p0)
	  USE_IT(p01,p1)
	  USE_IT(ICnt,count)
	  return true;
    }
  }
  count=0;
  USE_IT(ICnt,count)
  return false;
#undef USE_IT(x,val)
}


int CLine2D::intersects_m1(CLine2D const &l1,CLine2D const &l2, CPoint2D intersection){
	int result=0;
	intersects_m2(l1,l2,intersection,&result);
	return result;
}

bool  CLine2D::intersects_m2(const CLine2D& t_line, const CLine2D& o_line, CPoint2D& intersection, int *intType,float EPS){
#define USE_IT(x) if(intType!=NULL) {*intType=x;}
	 int Result;
	 CLineSegment2D this_line = t_line.toLineSegment2D();
	 CLineSegment2D other_line = o_line.toLineSegment2D();
	 CVec2D A=this_line.point1_.toVec2DPos();
	 CVec2D B=this_line.point2_.toVec2DPos();
	 CVec2D V=other_line.point1_.toVec2DPos();
	 CVec2D X=other_line.point2_.toVec2DPos();
	 CVec2D N=other_line.orthogonal();
	float num = (V-A)*N;
	float denom = ((B-A)*N);
	if (denom == 0){
		//se forem rays ou segmentos, determinar se há apenas um ponto de intersecção
		// se forem linhas, não é possivel calcular um ponto de intersecção
		intersection.x_=CMath::NAN_FLOAT;
		intersection.y_=CMath::NAN_FLOAT;
		Result= PARALLEL | COINCIDENT ; //coincidentes
		USE_IT(Result)
		return false;
	}
	float t = num / denom;
	CVec2D Q= A + ((B-A)*t);
	intersection=Q;
	// like this is a line to line intersection, if they ar concurrent they intersect.
		
		Result= CONCURRENT | INTERESECTING;
		USE_IT(Result)
		return true;   

#undef USE_IT(x)
}

bool  CLine2D::intersects_m2(const CLine2D& this_line, const CLineSegment2D& other_segment, CPoint2D& intersection, int *intType,float EPS){
	return CLineSegment2D::intersects_m2(other_segment,this_line,intersection,intType,EPS);
}

bool  CLine2D::intersects_m2(const CLine2D& t_line, const CRay2D& other_ray, CPoint2D& intersection, int *intType,float EPS){
	#define USE_IT(x) if(intType!=NULL) {*intType=x;}
	int Result;
	CLineSegment2D this_line= t_line.toLineSegment2D();
	CLineSegment2D other_line= other_ray.toLineSegment2D(100);
	 CVec2D A=this_line.point1_.toVec2DPos();
	 CVec2D B=this_line.point2_.toVec2DPos();
	 CVec2D V=other_line.point1_.toVec2DPos();
	 CVec2D X=other_line.point2_.toVec2DPos();
	 CVec2D N=other_line.orthogonal();
	float num = (V-A)*N;
	float denom = ((B-A)*N);
	if (denom == 0){
		intersection.x_=CMath::NAN_FLOAT;
		intersection.y_=CMath::NAN_FLOAT;
		// verificar se há apenas um ponto de interseção entre segmentos coincidentes
		// the segmnt could be parallel or coincident. and has many intersection points
		//se forem rays ou segmentos, determinar se há apenas um ponto de intersecção
		// se forem linhas, não é possivel calcular um ponto de intersecção
		// se são coincidentes há um t que satisfaça Vt=A
		//calcular t
		// verificar se são colineares
		bool col1 = this_line.point1_.IsCollinear(other_line.point1_,other_line.point2_);
		bool col2 = this_line.point2_.IsCollinear(other_line.point1_,other_line.point2_);
		if (col1 && col2)  {//são coincidentes. //verificar então se há apenas um ponto de intersecção
			  //naõ é possivel deteminar a intersecção
				Result= COINCIDENT ;
				USE_IT(Result)
				return false;
				
			
		}else{ //são paralelos
			Result= PARALLEL | NOT_INTERESECTING;
			USE_IT(Result)
			return false;
		}
		
	}
	float t = num / denom;
	CVec2D Q= A + ((B-A)*t);
	intersection=Q;
	// like this is a line to ray intersection, chech if the point belongs to the ray.
	// if the poit belong to both then there are intersection 
	bool l1 = other_ray.isPointOnRay(Q);

	if (l1) {
		Result= CONCURRENT | INTERESECTING;
		USE_IT(Result)
		return true;
	}else {
		Result= CONCURRENT | NOT_INTERESECTING;
		USE_IT(Result)
		return false;
	}
#undef USE_IT(x)


}


CPoint2D CLine2D::closestPoint( CLine2D const&other) const
{
	
	CPoint2D intersecPoint;
	int result;
	intersects_m2(*this,other,intersecPoint,&result);
	if (SMF_CONTAINS(result,CONCURRENT) ) {
		return CVec2D(intersecPoint);
	}else {  // parallel coincidents
		CPoint2D p2;
		get_two_points(intersecPoint,p2);
		return CVec2D(intersecPoint);
	}
	 

}

CPoint2D CLine2D::closestPoint( CLineSegment2D const&other) const
{
	
	CPoint2D intersecPoint;
	int result=0;
	CLineSegment2D seg1= toLineSegment2D();
	
	 CLineSegment2D::intersects_m1(seg1,other,intersecPoint,&result);
	if ( SMF_CONTAINS(result, INTERESECTING)) {
		return CVec2D(intersecPoint);
	}else if(SMF_CONTAINS(result, NOT_INTERESECTING)) {
		// verificar se o ponto retornado pertence não pertence ao segmento
		// se não pertencer calcular a distância aos vértices do segmento
		bool intersecSegment = other.intersects(intersecPoint);

		if(!intersecSegment) {
		float pa= distance(other.point1());
		float pb= distance(other.point2());
		if ( pa <= pb) return closestPoint(other.point1());
		if ( pa > pb) return closestPoint(other.point2());
		}else return intersecPoint;
	}else {  // parallel coincidents
			CPoint2D p2;
		get_two_points(intersecPoint,p2);
		return intersecPoint;
	}
	 

}

CPoint2D CLine2D::closestPoint( CRay2D const&other) const
{
	
	CPoint2D intersecPoint;
	int result;
	CLineSegment2D seg1= toLineSegment2D();
	CLineSegment2D::intersects_m2(seg1,other,intersecPoint,&result);
	if ( SMF_CONTAINS(result, INTERESECTING)) {
		return CVec2D(intersecPoint);
	}else if(SMF_CONTAINS(result, NOT_INTERESECTING)) {
		// verificar se o ponto retornado pertence não pertence ao segmento
		// se não pertencer calcular a distância aos vértices do segmento
		bool intersecSegment = other.intersects(intersecPoint);

		if(!intersecSegment) {
		return closestPoint(other.pos);
		
		}else return CVec2D(intersecPoint);
	}else {  // parallel coincidents
			CPoint2D p2;
		get_two_points(intersecPoint,p2);
		return CVec2D(intersecPoint);
	}
	 

}


namespace _2D{
	#define vp(os,v,s) { os<<' '; if ((v)>0) os<<'+'; if ((v)&&!s[0]) os<<(v); else { \
                     if ((v)==-1) os<<'-';\
                     else if ((v)!=0&&(v)!=1) os<<(v);\
                     if ((v)!=0) os<<' '<<s; } }

//: Write line description to stream: "<CLine2D ax+by+c=0>"

ostream&  operator<<(ostream& os, CLine2D const& l)
{
  os << "<CLine2D"; vp(os,l.a(),"x"); vp(os,l.b(),"y"); vp(os,l.c(),"");
  return os << " = 0 >";
}

#undef vp

//: Read in three line parameters from stream
//  Either just reads three blank-separated numbers,
//  or reads three comma-separated numbers,
//  or reads three numbers in parenthesized form "(123, 321, 567)"
//  or reads the formatted form "123x+321y+567=0"

istream&  operator>>(istream& is, CLine2D& line)
{
  if (! is.good()) return is; // (TODO: should throw an exception)
  bool paren = false;
  bool formatted = false;
  float a, b, c;
  is >> std::ws; // jump over any leading whitespace
  if (is.eof()) return is; // nothing to be set because of EOF (TODO: should throw an exception)
  if (is.peek() == '(') { is.ignore(); paren=true; }
  is >> std::ws >> a >> std::ws;
  if (is.eof()) return is;
  if (is.peek() == ',') is.ignore();
  else if (is.peek() == 'x') { is.ignore(); formatted=true; }
  is >> std::ws >> b >> std::ws;
  if (is.eof()) return is;
  if (formatted) {
    if (is.eof()) return is;
    if (is.peek() == 'y') is.ignore();
    else                  return is; // formatted input incorrect (TODO: throw an exception)
  }
  else if (is.peek() == ',') is.ignore();
  is >> std::ws >> c >> std::ws;
  if (paren) {
    if (is.eof()) return is;
    if (is.peek() == ')') is.ignore();
    else                  return is; // closing parenthesis is missing (TODO: throw an exception)
  }
  if (formatted) {
    if (is.eof()) return is;
    if (is.peek() == '=') is.ignore();
    else                  return is; // closing parenthesis is missing (TODO: throw an exception)
    is >> std::ws;
    if (is.peek() == '0') is.ignore();
    else                  return is; // closing parenthesis is missing (TODO: throw an exception)
  }
  line.set(a,b,c);
  return is;
}

} //end 2D
} //end GEO
}  //end SMF
