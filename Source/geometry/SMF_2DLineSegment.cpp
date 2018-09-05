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
#include "geometry/SMF_2DLineSegment.h"
#include "math/all.h"
#include "geometry/all.h"
#include <vector>



namespace SMF {

using namespace MATH;
using namespace Util;

namespace GEO{

/**
how to calculate the intersection of a segment and a line?
 
to calculate the intersection of a segment and a line, we refer to the equation of a segment. P = A + (B – A).t.
 
the point of intersection Q will be a point on that line, so Q = A + (B – A).t
 
and the point of intersection Q is such that (Q – V).N = 0, where V is a vertex of the edge we are trying to intersect, and N is the normal of that edge.
 
so, we have a system of two equations
 
(1)     Q = A + (B – A).t
(2)     (Q – V).N = 0
 
now, we can calculate t, the unkown, by substituting Q from equation 1 into equation (2). the result will be .
 
t = (V – A).N / (B-A).N
 
and you find Q, the point of intersection, by plugging t back into the first equation
 
Q = A + (B – A).t

**/
bool  CLineSegment2D::intersects_m2(const CLineSegment2D& this_line, const CLineSegment2D& other_line, CPoint2D& intersection, int *intType,float EPS){
#define USE_IT(x) if(intType!=NULL) {*intType=x;}
	int Result;
	CVec2D A=this_line.point1_.toVec2DPos();
	 CVec2D B=this_line.point2_.toVec2DPos();
	 CVec2D V=other_line.point1_.toVec2DPos();
	 CVec2D X=other_line.point2_.toVec2DPos();
	 CVec2D N=(V-X).orthogonalizePos();
	float num = (V-A)*N;
	float denom = ((B-A)*N);
	if (denom == 0){
		intersection.x_=CMath::NAN_FLOAT;
		intersection.y_=CMath::NAN_FLOAT;
		// verificar se há apenas um ponto de interseção entre segmentos coincidentes
		// the segmnt could be parallel or coincident. and has many intersection points
		//se forem rays ou segmentos, determinar se há apenas um ponto de intersecção
		// se forem linhas, não é possivel calcular um ponto de intersecção
		// verificar se são colineares
		bool col1 = this_line.point1_.IsCollinear(other_line.point1_,other_line.point2_);
		bool col2 = this_line.point2_.IsCollinear(other_line.point1_,other_line.point2_);
		if (col1 && col2) {//são coincidentes. //verificar então se há apenas um ponto de intersecção
		   	       IntersectResult intersec = this_line.isPointOnSegment_m1(other_line.point1_);
		   if(intersec == INTERESECTING_EXTREMITY_P1 || intersec == INTERESECTING_EXTREMITY_P2){
			   intersection=other_line.point1_;
			   Result = COINCIDENT | INTERESECTING | INTERESECTING_EXTREMITY_P1;
			   USE_IT(Result)
			   return true;
		   }else if (intersec == INTERESECTING ){ // não é possivel calcular ponto de intersecção
			   Result =  COINCIDENT | INTERESECTING;
			   USE_IT(Result)
			   return true;
		   }else if (intersec ==NOT_INTERESECTING){
		       intersec = this_line.isPointOnSegment_m1(other_line.point2_);
			   if(intersec == INTERESECTING_EXTREMITY_P1 || intersec == INTERESECTING_EXTREMITY_P2){
					intersection=other_line.point2_;
					Result =  COINCIDENT | INTERESECTING | INTERESECTING_EXTREMITY_P2;
					USE_IT(Result)
					return true;
			   }else if (intersec == INTERESECTING ){ // não é possivel calcular ponto de intersecção
					Result =  COINCIDENT | INTERESECTING;
					USE_IT(Result)
					return true;
			   }else if (intersec ==NOT_INTERESECTING){ // não é possivel calcular ponto de intersecção
					Result =  COINCIDENT | NOT_INTERESECTING;
					// todo calcular o ponto que seria de interseção, se houvsse uma
					USE_IT(Result)
					return false;
			   }
		   }
		
		}else{ //são paralelos
			Result = PARALLEL;
			USE_IT(Result)
			return false;
		}
		
		
	}
	float t = num / denom;
	CVec2D Q= A + ((B-A)*t);
	intersection=Q;
	// like this is a segment to segment intersection, chech if the point belongs to any of the segments.
	// if the poit belong to both then there are intersection 
	bool l1 =  (t >=0 && t <=1);
	bool l2 = other_line.isPointOnSegment(Q);
	if ( l1 && l2) {
		// todo: verificar se são também coincidentes
		
		Result= CONCURRENT | INTERESECTING;
		USE_IT(Result)
		return true;   
	}else{
		
		Result =CONCURRENT | NOT_INTERESECTING;
		USE_IT(Result)
		return false;
	}
#undef USE_IT(x)
}
bool  CLineSegment2D::intersects_m2(const CLineSegment2D& this_line, const CLine2D& o_line, CPoint2D& intersection, int *intType,float EPS){
#define USE_IT(x) if(intType!=NULL) {*intType=x;}
	 int Result;
	 CLineSegment2D other_line = o_line.toLineSegment2D();
	 CVec2D A=this_line.point1_.toVec2DPos();
	 CVec2D B=this_line.point2_.toVec2DPos();
	 CVec2D V=other_line.point1_.toVec2DPos();
	 CVec2D X=other_line.point2_.toVec2DPos();
	 CVec2D N=(V-X).orthogonalizePos();
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
	// like this is a segment to line intersection, chech if the point belongs to  the segments.
	bool l1 = this_line.isPointOnSegment(Q);
	if ( l1) {
		// todo: verificar se são também coincidentes
		
		Result= CONCURRENT | INTERESECTING;
		USE_IT(Result)
		return true;   
	}else{
		
		Result =CONCURRENT | NOT_INTERESECTING;
		USE_IT(Result)
		return false;
	}
#undef USE_IT(x)
}
bool  CLineSegment2D::intersects_m2(const CLineSegment2D& this_line, const CRay2D& other_ray, CPoint2D& intersection, int *intType,float EPS){
#define USE_IT(x) if(intType!=NULL) {*intType=x;}
	int Result; 
	CLineSegment2D other_line= other_ray.toLineSegment2D(100);
	 CVec2D A=this_line.point1_.toVec2DPos();
	 CVec2D B=this_line.point2_.toVec2DPos();
	 CVec2D V=other_line.point1_.toVec2DPos();
	 CVec2D X=other_line.point2_.toVec2DPos();
	 CVec2D N=(V-X).orthogonalizePos();
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
		if (col1 && col2) {//são coincidentes. 
			//verificar se há apenas um ponto de intersecção
			//calcular o ponto
			IntersectResult intersec = this_line.isPointOnSegment_m1(other_ray.pos);
		   if(intersec == INTERESECTING_EXTREMITY_P1 || intersec == INTERESECTING_EXTREMITY_P2){
			   intersection=other_ray.pos;
			   Result = COINCIDENT | INTERESECTING | INTERESECTING_EXTREMITY_P1;
			   USE_IT(Result)
			   return true;
		   }else if (intersec == INTERESECTING ){ // não é possivel calcular ponto de intersecção
			   Result =  COINCIDENT | INTERESECTING;
			   USE_IT(Result)
			   return true;
		   }else if (intersec ==NOT_INTERESECTING){
		        // não há possivel calcular ponto de intersecção
					Result =  COINCIDENT;
					USE_IT(Result)
					return false;
			  
		   }
			
		}else{ //são paralelos
			Result= PARALLEL | NOT_INTERESECTING;
			USE_IT(Result)
			return false;
		}
		
	}
	float t = num / denom;
	CVec2D Q= A + ((B-A)*t);
	intersection=Q;
	// like this is a segment to ray intersection, chech if the point belongs to both of the segments.
	// The point is calc based on segment values, so it belongs to the segment, if t >=0 and <=t.
	// if the poit belong to both then there are intersection 
	bool l1 = (t >=0 && t <=1);
	bool l2 = other_ray.isPointOnRay(Q);

	if ( l1 && l2) {
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

bool CLineSegment2D::intersects_m1(const CLineSegment2D& this_line, const CLineSegment2D& other_line, CPoint2D& intersection, int *intType,float EPS){
#define USE_IT(x) if(intType!=NULL) {*intType=x;}
	int Result;	    
	float x1,x2,x3,x4,y1,y2,y3,y4;
	x1=this_line.point1_.x();
	x2=this_line.point2_.x();
	y1=this_line.point1_.y();
	y2=this_line.point2_.y();
	x3=other_line.point1_.x();
	x4=other_line.point2_.x();
	y3=other_line.point1_.y();
	y4=other_line.point2_.y();

	
	float denom = ((y4-y3)*(x2-x1))-((x4-x3)*(y2-y1));
    
	float nume_b = ((x2-x1)*(y1-y3))-((y2-y1)*(x2-x3));

    float nume_a = ((x4-x3)*(y1-y3))-((y4-y3)*(x1-x3));

   /* Are the line coincident? */
   if ((CMath::abs(nume_a) < EPS || CMath::abs(nume_b) < EPS) && ABS(denom) < EPS) {
	   	intersection.x_=CMath::NAN_FLOAT;
		intersection.y_=CMath::NAN_FLOAT;

	   //TODO: CLassificar os coincidentes, pois podem ser:
				// só com 1 ponto de encontro   ----.-----
				// com vários pontos de encontro   ----====-----
				//sem nenhum ponto de encontro   ---- -------
				// aqui deve verificar se a coincidencia é apenas no ponto 
				// final ou inicial do segmento e avisar.
	       IntersectResult intersec = this_line.isPointOnSegment_m1(other_line.point1_);
		   if(intersec == INTERESECTING_EXTREMITY_P1 || intersec == INTERESECTING_EXTREMITY_P2){
			   intersection=other_line.point1_;
				Result = COINCIDENT | INTERESECTING | INTERESECTING_EXTREMITY_P1;
				USE_IT(Result)
				return true;

		   }else if (intersec == INTERESECTING ){ // não é possivel calcular ponto de intersecção
			   Result= COINCIDENT | INTERESECTING;
			   USE_IT(Result)
				return true;

		   }else if (intersec ==NOT_INTERESECTING){
		       intersec = this_line.isPointOnSegment_m1(other_line.point2_);
			   if(intersec == INTERESECTING_EXTREMITY_P1 || intersec == INTERESECTING_EXTREMITY_P2){
					intersection=other_line.point2_;
					Result = COINCIDENT | INTERESECTING | INTERESECTING_EXTREMITY_P2;
					USE_IT(Result)
					return true;

			   }else if (intersec == INTERESECTING ){ // não é possivel calcular ponto de intersecção
					Result = COINCIDENT | INTERESECTING;
					USE_IT(Result)
					return true;

			   }else if (intersec ==NOT_INTERESECTING){ // não é possivel calcular ponto de intersecção
					Result= COINCIDENT;
					USE_IT(Result)
					return false;

			   }
		   }
	   
	   
		Result= COINCIDENT;  
		USE_IT(Result)
		return false;

   }


   /* Are the line parallel */
   if (CMath::fabs(denom) < EPS) {
		intersection.set(CMath::NAN_FLOAT,CMath::NAN_FLOAT);
		Result= PARALLEL;
		USE_IT(Result)
		return false;
   
   }

		/* Is the intersection along the the segments */
   
        float ua = nume_a / denom;
        float ub = nume_b / denom;

        if(ua >= 0.0f && ua <= 1.0f && ub >= 0.0f && ub <= 1.0f)
        {
			float xx,yy;
            // Get the intersection point.
			if (ua==0.0f){

			xx= other_line.point1_.x();
			yy= other_line.point1_.y();
			}else{
			xx= x1 + ua * (x2 - x1);
			yy= y1 + ua * (y2 - y1);
			}

			//converte os segmentos para retas e utiliza as formulas das retas
			
			intersection.set( xx , yy);

            Result= CONCURRENT | INTERESECTING;
			USE_IT(Result)
			return true;

        }

		if(ua>= 0.0f && ua	<=1.0f){ //the intersection is should be on line at this_line
			/* Get the intersection point.
			x = x1 + ua (x2 - x1)
			y = y1 + ua (y2 - y1)  */ 
			float xx= x1 + ua * (x2 - x1);
			float yy= y1 + ua * (y2 - y1);
			
			//float xxx=other_line.point1_.x() + ua*(other_line.point2_.x() - other_line.point1_.x());
			//float yyy=other_line.point1_.y() + ua*(other_line.point2_.y() - other_line.point1_.y());
			
			
			intersection.set( xx , yy);

		}
		if(ub>= 0.0f && ub	<=1.0f){ //the intersection is should be on line at other_line
			/* Get the intersection point.
			x = x3 + ub (x3 - x4)
			y = y3 + ub (y3 - y4)  */
			float xx=x3 + ub*(x3 - x4);
			float yy=y3 + ub*(y3 - y4);
			intersection.set( xx, yy );
            


		}
        Result= CONCURRENT | NOT_INTERESECTING;
		USE_IT(Result)
		return false;

#undef USE_IT(x)
}



bool CLineSegment2D::intersects(const CLineSegment2D& this_line, const CLineSegment2D& other_line, CPoint2D& intersection, float EPS){

CLine2D l1= this_line.toLine2D();
CLine2D l2= other_line.toLine2D();
bool result = CLine2D::intersects(l1,l2,intersection,1e-3f);
//verifica se o ponto pertence ao segmento.
bool result2 = this_line.intersects(intersection);  //problema na verificação
return result2;
}

bool CLineSegment2D::intersects(const CLineSegment2D& other_line, CPoint2D& intersection)
{
		
		return intersects_m1(*this,other_line,intersection);
}

bool CLineSegment2D::intersects(CLineSegment2D const &l)const{
#if 0

	CLine2D l1 = this->toLine2D();
	CLine2D l2 = l.toLine2D();
	CPoint2D intersects;
	bool result1 =CLine2D::intersects(l1,l2,intersects);
	if(result1==true){
		//check if the point belongs to this linesegment;
		bool result2 = this->intersects(intersects);
		return result2;
	}else return false;
#endif
	return simpleIntersects(*this,l);
}


// helper function to calc intersection
bool CLineSegment2D::simpleIntersects(const CPoint2D &Point1,const CPoint2D &Point2,const CPoint2D &Point3,const CPoint2D &Point4)
{
	bool Result = (
             ((CPoint2D::orientation(Point1,Point2,Point3) * CPoint2D::orientation(Point1,Point2,Point4)) <= 0) &&
             ((CPoint2D::orientation(Point3,Point4,Point1) * CPoint2D::orientation(Point3,Point4,Point2)) <= 0)
            );
	return Result;
}

bool CLineSegment2D::simpleIntersects(const CLineSegment2D &Segment1,const CLineSegment2D &Segment2)
{
  bool Result = simpleIntersects(Segment1.point1_,Segment1.point2_,Segment2.point1_,Segment2.point2_);
  return Result;
}
bool CLineSegment2D::intersects(CPoint2D const &p)const{


int res = isPointOnSegment_m1(p);
if (SMF_CONTAINS(res,INTERESECTING)  ||SMF_CONTAINS(res,INTERESECTING_EXTREMITY_P1) || SMF_CONTAINS(res,INTERESECTING_EXTREMITY_P2)) return true;
else return false;


}

bool CLineSegment2D::intersects(CLine2D const &l)const{
	return l.intersects(*this);
}


bool CLineSegment2D::intersects(CRay2D const &l)const{
		CLine2D l1 = this->toLine2D();
	CLine2D l2 = l.toLine2D();
	CPoint2D intersects;
	bool result1 =CLine2D::intersects(l1,l2,intersects);
	if(result1==true){
		//check if the point belongs to this linesegment;
		bool result2 = this->intersects(intersects);
		return result2;
	}else return false;

}

bool CLineSegment2D::intersection(const CAABBox2D& box,int *ICnt,  CPoint2D* p00, CPoint2D* p01)const
{
#define USE_IT(x,val) if(x!=NULL) {*x=val;}
	if (p00!=NULL) p00->toInfinite();
	if (p01!=NULL) p01->toInfinite();
	if (ICnt!=NULL) *ICnt=0;
   int nint = 0;
  CLine2D line(a(), b(), c());
  CPoint2D pi0, pi1;
  // if no intersection just return
	if (!box.intersection(line, NULL,&pi0, &pi1)){
	  USE_IT(ICnt,nint)
      return false;
	}
  
  // check if intersection points are interior to the line segment
  if (intersects(pi0)) {
    USE_IT(p00,pi0)
	nint++;
  }
  if (intersects(pi1)) {
    USE_IT(p01,pi1)
    nint++;
  }
  USE_IT(ICnt,nint)
  return nint > 0;

#undef USE_IT(x,val)
}

bool CLineSegment2D::isPointOnSegment(CPoint2D const& p, int method)const{
	if (method<1) {
		IntersectResult ret = isPointOnSegment_m1(p);
		if(SMF_CONTAINS(ret,INTERESECTING) || SMF_CONTAINS(ret,INTERESECTING_EXTREMITY_P1) || SMF_CONTAINS(ret,INTERESECTING_EXTREMITY_P2)) return true;
		else return false;
	}
	else if (method>=1) return isPointOnSegment_m2(p);
}



bool CLineSegment2D::isPointOnSegment_m2(CPoint2D const& p)const
{
  CLineSegment2D const& lseg=*this;
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
  return CMath::nearZero (diff);// <= tolerance;
}

IntersectResult CLineSegment2D::isPointOnSegment_m1(CPoint2D const &point)const{
	
    // A and B are the extremities of the current segment
    // C is the point to check
 
    // Create the vector AB
    CVec2D AB(point1_,point2_);
    // Create the vector AC
    CVec2D AC(point1_,point);
 
    // Compute the cross product of VA and PAP
    // Check if the three points are aligned (cross product is 0)
	//use a big tolerance to avoid rounding errors
    bool colinear= point.IsCollinear(point1_,point2_,0.1);
	if (!colinear) return NOT_INTERESECTING;
	//if (!( AB.cross(AC).isNull())) return NOT_INTERESECTING;
 
 
    // Compute the dot product of vectors
    double KAC = AB*AC;
    if (KAC<0) return NOT_INTERESECTING;
    if (KAC==0) return INTERESECTING_EXTREMITY_P1;
 
    // Compute the square of the segment length
    double KAB=AB*AB;
    if (KAC>KAB) return NOT_INTERESECTING;
    if (KAC==KAB) return INTERESECTING_EXTREMITY_P2;
 
    // The point is on the segment
    return INTERESECTING;
}

CLine2D CLineSegment2D::toLine2D()const{
	return CLine2D(a(),b(),c());
}

// stream operators

float CLineSegment2D::a() const
{
  return point1_.y()-point2_.y();
}


float CLineSegment2D::b() const
{
  return point2_.x()-point1_.x();
}


float CLineSegment2D::c() const
{
  return (point1_.x()*point2_.y())-(point2_.x()*point1_.y());
}

CVec2D CLineSegment2D::directionVector()const{
	CVec2D posVec1(point1_);
	CVec2D posVec2(point2_);
	return (posVec2-posVec1).normalized();
}

	/// Projects this CLineSegment2D onto the given 1D axis od carthesian.
CLineSegment2D CLineSegment2D::projectToXAxis() const{
	CPoint2D p1,p2;
	p1.set(point2_.x()-point1_.x(),0);
	p2.set(point2_.y()-point1_.y(),0);
	return CLineSegment2D(p1,p2);

}
CLineSegment2D CLineSegment2D::projectToYAxis() const{
	CPoint2D p1,p2;
	p1.set(0,point2_.x()-point1_.x());
	p2.set(0,point2_.y()-point1_.y());
	return CLineSegment2D(p1,p2);

}

CVec2D  CLineSegment2D::direction() const
{
  CVec2D v(point2_.x()-point1_.x(),point2_.y()-point1_.y());
  v.toNormal();
  return v;
}


CVec2D  CLineSegment2D::normal() const
{
  return normalPos();
}

CVec2D  CLineSegment2D::normalPos() const
{
  CVec2D v(point1_.y()-point2_.y(),point2_.x()-point1_.x());
  v.toNormal();
  return v;
}

CVec2D  CLineSegment2D::normalNeg() const
{
  CVec2D v(point2_.y()-point1_.y(),point1_.x()-point2_.x());
  v.toNormal();
  return v;
}

CVec2D CLineSegment2D::orthogonal() const{
  return orthogonalCCW();
}

CVec2D CLineSegment2D::orthogonalCCW() const{
  CVec2D v(point1_.y()-point2_.y(),point2_.x()-point1_.x());
  return v;
}

CVec2D CLineSegment2D::orthogonalCW() const{
  CVec2D v(point2_.y()-point1_.y(),point1_.x()-point2_.x());
  return v;
}

CLine2D CLineSegment2D::orthogonal(CPoint2D const &p)const{
	return CLine2D( b(), -a(), (a() * p.y_ - b() * p.x_) );

}

CPoint2D CLineSegment2D::midPoint()const{
	CPoint2D mid = point1_.midPoint(point2_);
	return mid;
}

double CLineSegment2D::slope_degrees() const
{
	static const double deg_per_rad = 45.0/CMath::atan64(1.0,1.0);
  double dy = point2_.y()-point1_.y();
  double dx = point2_.x()-point1_.x();
  // do special cases separately, to avoid rounding errors:
  if (dx == 0) return dy<0 ? -90.0 : 90.0;
  if (dy == 0) return dx<0 ? 180.0 : 0.0;
  if (dy == dx) return dy<0 ? -135.0 : 45.0;
  if (dy+dx == 0) return dy<0 ? -45.0 : 135.0;
  // general case:
  return deg_per_rad * CMath::atan64(dy,dx);
}


double CLineSegment2D::slope_radians() const
{
  double dy = point2_.y()-point1_.y();
  double dx = point2_.x()-point1_.x();
  return CMath::atan64(dy,dx);
}

float CLineSegment2D::distanceSegmentPoint(CLineSegment2D const &l1,
                                    CPoint2D const &p)
{
	float x1 = l1.point1_.x();
	float y1 = l1.point1_.y();
	float x2 = l1.point2_.x();
	float y2 = l1.point2_.y();
	float x = p.x_;
	float y = p.y_;
	// squared distance between endpoints :
	float ddh = CMath::square64(x2-x1) + CMath::square64(y2-y1);

  // squared distance to endpoints :
  float dd1 = CMath::square64(x-x1) + CMath::square64(y-y1);
  float dd2 = CMath::square64(x-x2) + CMath::square64(y-y2);

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
  return CMath::square64(a*x + b*y + c)/float(a*a + b*b);
}

float CLineSegment2D::distance(CPoint2D const& p)const
{
	CLineSegment2D const& l=*this;
	return CMath::sqrt(distanceSegmentPoint(l,p));
}

float CLineSegment2D::distanceSq(CPoint2D const& p)const
{
	CLineSegment2D const& l=*this;
	return distanceSegmentPoint(l,p);
}


CPoint2D CLineSegment2D::closestPoint(CPoint2D const  &p, float *tparam)const{
	//return p.closestPoint(*this);
	CPoint2D d;
	CVec2D b = point2_.toVec2DPos();
	CVec2D a = point1_.toVec2DPos();
	CVec2D ab= b - a;
	// project c onto ab, computing parameterized position d(t)=a+t*(b – a)
	float t = ((p - a)* ab) ;
	// If outside segment, clamp t (and therefore d) to the closest endpoint
if (t <= 0.0f) {
	// c projects outside the [a,b] interval, on the a side; clamp to a
	t = 0.0f;
	d=a;
} else {
	float denom = (ab)*(ab);
	// Always nonnegative since denom = ||ab||∧2
	if (t >= denom) {
		// c projects outside the [a,b] interval, on the b side; clamp to b
		t = 1.0f;
		d=b;
	} else {
		// c projects inside the [a,b] interval; must do deferred divide now
		t=t/denom;
		d=a+t*ab;
	}
	if(tparam) *tparam=t;
	return d;
}

}
CPoint2D CLineSegment2D::closestPoint(const CLineSegment2D &other)const{

	CPoint2D intersecPoint;
	int result=0;
	
	 CLineSegment2D::intersects_m1(*this,other,intersecPoint,&result);
	if ( SMF_CONTAINS(result, INTERESECTING)) {
		return intersecPoint;
	}else if(SMF_CONTAINS(result, NOT_INTERESECTING)) {
		// verificar se o ponto retornado pertence não pertence ao segmento
		// se não pertencer calcular a distância aos vértices do segmento
		bool intersecSegment = intersects(intersecPoint);

		if(!intersecSegment) {
		float pa= distance(other.point1());
		float pb= distance(other.point2());
		if ( pa <= pb) return closestPoint(other.point1());
		if ( pa > pb)  return closestPoint(other.point2());
		}else return intersecPoint;
	}else {  // parallel coincidents
		
		return point1();
	}

}


CAABBox2D CLineSegment2D::toCAABBox()const{
	
	float minX,minY,maxX,maxY;
	float x1,y1,x2,y2;
	x1=point1_.x();
	y1=point1_.y();
	x2=point2_.x();
	y2=point2_.y();


	if (x1==x2 && y1==y2) { //point
	minX = x1;
	minY = y1;
	maxX = x1;
	maxY = y1;

	}
	if (x1==x2) { //vertical segment
	minX = x1;
	minY = MIN(y1,y2);
	maxX = x1;
	maxY = MAX(y1,y2);
	}
	if (y1==y2){ //horizontal line
	minX = MIN(x1,x2);
	minY = y1;
	maxX = MAX(x1,x2);
	maxY = y1;
	}
	minX = MIN(x1,x2);
	minY = MIN(y1,y2);
	maxX = MAX(x1,x2);
	maxY = MAX(y1,y2);
  
  return CAABBox2D(minX,maxX,minY,maxY);
}   



namespace _2D{
ostream& operator<<(ostream& s, CLineSegment2D const & p)
{
  return s << "<CLineSegment2D " << p.point1() << " to " << p.point2() << " >";
}

istream& operator>>(istream& s, CLineSegment2D& p)
{
  CPoint2D p1, p2;
  s >> p1 >> p2;
  p.set(p1, p2);
  return s;
}

} //end _2D
} //end GEO
}  //end SMF
