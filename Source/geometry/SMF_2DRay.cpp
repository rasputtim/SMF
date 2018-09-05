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
#include "geometry/SMF_2DRay.h"
#include "math/all.h"
#include "geometry/all.h"


namespace SMF {
using namespace MATH;
using namespace Util;
namespace GEO{



CRay2D::CRay2D(const CPoint2D &pos_, const CVec2D &dir_)
:pos(pos_), dir(dir_)
{
	if(!dir.isNormalized()) dir.toNormal();
}


CRay2D::CRay2D(const CLineSegment2D &lineSegment, int pointToUse)
	:pos(pointToUse >=1? lineSegment.point1() :lineSegment.point2()), dir(lineSegment.directionVector())
{
}

CRay2D::CRay2D(const CPoint2D &pos_, const CPoint2D &p2):pos(pos_), dir(p2-pos_){
	if(!dir.isNormalized()) dir.toNormal();
}

bool CRay2D::isFinite() const
{
	return CVec2D(pos).isFinite() && dir.isFinite();
}

CVec2D CRay2D::getPoint(float d) const
{
	SMF_ASSERT(dir.isNormalized());
	return pos + d * dir;
}

void CRay2D::translate(const CVec2D &offset)
{
	using _2D::operator+=;
	pos += offset;
}
#if 0
void CRay2D::transform(const CMat3D &transform)
{
	pos = transform.transform(pos);
	dir = transform.transform(dir);
}

void CRay2D::transform(const CMatJoint3x4 &transform)
{
	pos = transform.transformPos(pos);
	dir = transform.transformDir(dir);
}

void CRay2D::transform(const CMat4D &transform)
{
	pos = transform.transformPos(pos);
	dir = transform.transformDir(dir);
}

void CRay2D::transform(const CQuaternion &transform)
{
	pos = transform.transform(pos);
	dir = transform.transform(dir);
}
#endif
bool CRay2D::contains(const CPoint2D &point, float distanceThreshold) const
{
	return closestPoint(CVec2D(point)).distanceSq(CVec2D(point)) <= distanceThreshold;
}

bool CRay2D::contains(const CLineSegment2D &lineSegment, float distanceThreshold) const
{
	return contains(lineSegment.point1(), distanceThreshold) && contains(lineSegment.point2(), distanceThreshold);
}

bool CRay2D::compare(const CRay2D &rhs, float epsilon) const
{
	return pos.compare(rhs.pos) && dir.compare(rhs.dir, epsilon);
}

float CRay2D::distance(const CVec2D &point, float *d) const
{
	return closestPoint(point, d).distance(point);
}

float CRay2D::distance(const CVec2D &point) const
{
	return distance(point, 0);
}
#if 0
float CRay2D::distance(const CRay2D &other, float *d, float *d2) const
{
	float u2;
	CVec2D c = closestPoint(other, d, &u2);
	if (d2) *d2 = u2;
	return c.distance(other.getPoint(u2));
}

float CRay2D::distance(const CRay2D &ray) const
{
	return distance(ray, 0, 0);
}
#endif
#if 0
float CRay2D::distance(const CLine2D &other, float *d, float *d2) const
{
	float u2;
	CVec2D c = closestPoint(other, d, &u2);
	if (d2) *d2 = u2;
	return c.distance(other.getPoint(u2));
}

float CRay2D::distance(const CLine2D &line) const
{
	return distance(line, 0, 0);
}
#endif
#if 0
float CRay2D::distance(const CLineSegment2D &other, float *d, float *d2) const
{
	float u2;
	CVec2D c = closestPoint(other, d, &u2);
	if (d2) *d2 = u2;
	return c.distance(other.getPoint(u2));
}

float CRay2D::distance(const CLineSegment2D &lineSegment) const
{
	return distance(lineSegment, 0, 0);
}
#endif
float CRay2D::distance(const CCircle2D &circle) const
{
	return MAX(0.f, distance(circle.getOrigin().toVec2DPos()) - circle.getRadius());
}
#if 0
float CRay2D::distance(const Capsule &capsule) const
{
	return MAX(0.f, distance(capsule.l) - capsule.r);
}
#endif

//Todo Testar os dois métodos
CPoint2D CRay2D::closestPoint(CPoint2D const  &p)const{
	return p.closestPoint(*this);
}

CVec2D CRay2D::closestPoint(const CVec2D &targetPoint, float *d) const
{
	float u = MAX(0.f, ((targetPoint - CVec2D(pos))* dir));
	if (d)
		*d = u;
	return getPoint(u);
}

bool CRay2D::intersects(CPoint2D const &p)const{
bool res = isPointOnRay(p);
return res;


}

bool CRay2D::intersects(CLine2D const &l)const{
	return l.intersects(*this);
}

bool CRay2D::intersects(CLineSegment2D const &l)const{
	CLine2D l1 = this->toLine2D();
	CLine2D l2 = l.toLine2D();
	CPoint2D intersects;
	bool result1 =CLine2D::intersects(l1,l2,intersects);
	if(result1==true){
		//check if the point belongs to this ray;
		bool result2 = this->intersects(intersects);
		return result2;
	}else return false;

}

bool CRay2D::intersects(CRay2D const &l)const{
	CLine2D l1 = this->toLine2D();
	CLine2D l2 = l.toLine2D();
	CPoint2D intersects;
	bool result1 =CLine2D::intersects(l1,l2,intersects);
	if(result1==true){
		//check if the point belongs to this ray;
		bool result2 = this->intersects(intersects);
		return result2;
	}else return false;

}


bool CRay2D::isPointOnRay(CPoint2D const &point, int *intersecResult)const{
#define DELETE_IT if(deleteresult) {delete intersecResult; intersecResult=NULL;}
	bool deleteresult=false;
	if (intersecResult==NULL) {intersecResult=new int; deleteresult=true;}
    // A and B are the extremities of the current segment
    // C is the point to check
 // o segundo ponto é bem distante...
	CPoint2D p2= toLineSegment2D(100).point2();
    // Create the vector AB
	
    CVec2D AB(pos,p2);
    // Create the vector AC
    CVec2D AC(pos,point);
 
    // Compute the cross product of VA and PAP
    // Check if the three points are aligned (cross product is 0)
    bool colinear= point.IsCollinear(pos,p2);
	if (!colinear) {*intersecResult = NOT_INTERESECTING; DELETE_IT;return false; }
	//if (!( AB.cross(AC).isNull())) return NOT_INTERESECTING;
 
 //if (SMF_CONTAINS(res,INTERESECTING) || SMF_CONTAINS(res,INTERESECTING_EXTREMITY_P1 )|| SMF_CONTAINS(res,INTERESECTING_EXTREMITY_P2)) return true;

    // Compute the dot product of vectors
    double KAC = AB*AC;
    if (KAC<0) {*intersecResult = NOT_INTERESECTING; DELETE_IT;return false; }
	if (KAC==0) {*intersecResult = INTERESECTING_EXTREMITY_P1; DELETE_IT;return true;};
 
    // Compute the square of the segment length
    double KAB=AB*AB;
	if (KAC>KAB) {*intersecResult =  INTERESECTING; DELETE_IT;return true;}
	if (KAC==KAB) {*intersecResult =  INTERESECTING; DELETE_IT;return true;}
 
    // The point is on the segment
   *intersecResult =  INTERESECTING;
   DELETE_IT;
   return true;
}



#if 0
CVec2D CRay2D::closestPoint(const CRay2D &other, float *d, float *d2) const
{
	float u, u2;
	CVec2D closestPoint = CLine2D::closestPointLineLine(pos, pos + dir, other.pos, other.pos + other.dir, &u, &u2);
	if (u < 0.f && u2 < 0.f)
	{
		closestPoint = closestPoint(other.pos, &u);

		CVec2D closestPoint2 = other.closestPoint(pos, &u2);
		if (closestPoint.distanceSq(other.pos) <= closestPoint2.distanceSq(pos))
		{
			if (d)
				*d = u;
			if (d2)
				*d2 = 0.f;
			return closestPoint;
		}
		else
		{
			if (d)
				*d = 0.f;
			if (d2)
				*d2 = u2;
			return pos;
		}
	}
	else if (u < 0.f)
	{
		if (d)
			*d = 0.f;
		if (d2)
		{
			other.closestPoint(pos, &u2);
			*d2 = MAX(0.f, u2);
		}
		return pos;
	}
	else if (u2 < 0.f)
	{
		CVec2D pt = closestPoint(other.pos, &u);
		u = MAX(0.f, u);
		if (d)
			*d = u;
		if (d2)
			*d2 = 0.f;
		return pt;
	}
	else
	{
		if (d)
			*d = u;
		if (d2)
			*d2 = u2;
		return closestPoint;
	}
}

CVec2D CRay2D::closestPoint(const CLine2D &other, float *d, float *d2) const
{
	float t;
	CVec2D closestPoint = CLine2D::closestPointLineLine(pos, pos + dir, other.pos, other.pos + other.dir, &t, d2);
	if (t <= 0.f)
	{
		if (d)
			*d = 0.f;
		if (d2)
			other.closestPoint(pos, d2);
		return pos;
	}
	else
	{
		if (d)
			*d = t;
		return closestPoint;
	}
}

CVec2D CRay2D::closestPoint(const CLineSegment2D &other, float *d, float *d2) const
{
	float u, u2;
	CPoint2D intersecPoint;
	IntersectResult result;
	CLineSegment2D seg1= ToSegment2D();
	IntersectResult result =CLineSegment2D::intersects(seg1,other,intersecPoint);
	if (result == NOT_INTERESECTING || result == INTERESECTING) {
		return CVec2D(intersecPoint);
	}else {  // parallel coincidents
		intersecPoint = other.point1();
	}
	 

}
#endif

bool  CRay2D::intersects_m2(const CRay2D& this_ray, const CLine2D& other_line, CPoint2D& intersection, int *intType,float EPS){
	return CLine2D::intersects_m2(other_line,this_ray,intersection,intType,EPS);

}
bool  CRay2D::intersects_m2(const CRay2D& this_ray, const CLineSegment2D& other_segment, CPoint2D& intersection, int *intType,float EPS){
	return CLineSegment2D::intersects_m2(other_segment,this_ray,intersection,intType,EPS);

}
bool  CRay2D::intersects_m2(const CRay2D& this_ray, const CRay2D& other_ray, CPoint2D& intersection, int *intType,float EPS){
	#define USE_IT(x) if(intType!=NULL) {*intType=x;}
	

	int Result; 
	if(this_ray.compare(other_ray)){
    intersection.x_=CMath::NAN_FLOAT;
	intersection.y_=CMath::NAN_FLOAT;
	Result = COINCIDENT;
	USE_IT(Result)
	return true;
	}
	//eles não são o mesmo
	CLineSegment2D other_line= other_ray.toLineSegment2D(100);
	CLineSegment2D this_line= this_ray.toLineSegment2D(100);
	
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
		if (col1 && col2)  {//são coincidentes. //verificar então se há apenas um ponto de intersecção
			//verificar se há apenas um ponto de intersecção
				//calcular o ponto
				if (other_ray.pos == this_ray.pos) {
					
						intersection=this_ray.pos;
						Result = COINCIDENT | INTERESECTING | INTERESECTING_EXTREMITY_P1;
						USE_IT(Result)
						return true;
					
				}else {  //naõ é possivel deteminar a intersecção
				Result= COINCIDENT;
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
	// like this is a ray to ray intersection, chech if the point belongs to any of the rays.
	// if the poit belong to both then there are intersection 
	bool l1 = this_ray.isPointOnRay(Q);
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



bool CRay2D::intersection( const CCircle2D &circle,int *ICnt,CPoint2D * I1, CPoint2D * I2 ) const{

	return circle.intersection(*this,ICnt,I1,I2);
}

bool CRay2D::intersection( const CTriangle2D &triangle,int *ICnt,CPoint2D * I1, CPoint2D * I2 ) const{

	return triangle.intersection(*this,ICnt,I1,I2);
}
bool CRay2D::intersection(const CLine2D &line,  int *ICnt, CPoint2D *I1) const{
	CPoint2D pint;
	intersects_m2(*this,line,pint);
	if(!pint.isFinite()){
		if(ICnt!=NULL) {
			*ICnt=0;
		}
		if(I1!=NULL) {
			I1->toInfinite();
		}
		return false;
	}else{
		if(ICnt!=NULL) {
			*ICnt=1;
		}
		if(I1!=NULL) {
			*I1=pint;
		}
		return true;
	}
}
bool CRay2D::intersection(const CLineSegment2D &segment,  int *ICnt, CPoint2D *I1) const{
	CPoint2D pint;
	intersects_m2(*this,segment,pint);
	if(!pint.isFinite()){
		if(ICnt!=NULL) {
			*ICnt=0;
		}
		if(I1!=NULL) {
			I1->toInfinite();
		}
		return false;
	}else{
		if(ICnt!=NULL) {
			*ICnt=1;
		}
		if(I1!=NULL) {
			*I1=pint;
		}
		return true;
	}}
bool CRay2D::intersection(const CRay2D &ray,  int *ICnt, CPoint2D *I1) const{
		CPoint2D pint;
	intersects_m2(*this,ray,pint);
	if(!pint.isFinite()){
		if(ICnt!=NULL) {
			*ICnt=0;
		}
		if(I1!=NULL) {
			I1->toInfinite();
		}
		return false;
	}else{
		if(ICnt!=NULL) {
			*ICnt=1;
		}
		if(I1!=NULL) {
			*I1=pint;
		}
		return true;
	}
}




bool CRay2D::intersects(const CCircle2D &circle) const
{

	CPoint2D point(circle.center.x_,circle.center.y_);
	float dist;
	dist= distance(point.toVec2DPos());
	return  dist <= circle.r;

}
#if 0
bool CRay2D::intersects(const CTriangle2D &triangle) const
{
	float u, v;
	float t = CTriangle2D::intersectLineTri(pos, dir, triangle.a, triangle.b, triangle.c, u, v);
	if (t < 0.f || t == CMath::INFINITY_FLOAT)
		return false;
	return true;
}


bool CRay2D::intersects(const CAABBox2D &aabb) const
{
	return aabb.intersects(*this);
}

bool CRay2D::intersects(const CAABBox2D &aabb, float &dNear, float &dFar) const
{
	return aabb.intersects(*this, dNear, dFar);
}
#endif
#if 0
bool CRay2D::intersects(const COBBox2D &obb, float &dNear, float &dFar) const
{
	return obb.intersects(*this, dNear, dFar);
}

bool CRay2D::intersects(const COBBox2D &obb) const
{
	return obb.intersects(*this);
}
#endif
#if 0
bool CRay2D::intersects(const Capsule &capsule) const
{
	return capsule.intersects(*this);
}

bool CRay2D::intersects(const CPolygon &polygon) const
{
	return polygon.intersects(*this);
}

bool CRay2D::intersects(const CPolyhedron &polyhedron) const
{
	return polyhedron.intersects(*this);
}
#endif
#if 0
bool CRay2D::intersectsDisc(const CCircle2D &disc) const
{
	return disc.intersectsDisc(*this);
}
#endif
CLine2D CRay2D::toLine2D() const
{
	return CLine2D(pos, dir);
}

CLineSegment2D CRay2D::toLineSegment2D(float d) const
{
	return CLineSegment2D(pos, getPoint(d));
}

CLineSegment2D CRay2D::toLineSegment2D(float dStart, float dEnd) const
{
	return CLineSegment2D(getPoint(dStart), getPoint(dEnd));
}

void CRay2D::projectToAxis(const CVec2D &direction, float &outMin, float &outMax) const
{

	outMin = outMax = (direction* CVec2D(pos));
	float d = (direction* dir);

	// Most of the time, the projection interval will be a half-infinite range, extending to either -inf or +inf.
	if (d > 1e-4f)
		outMax = CMath::INFINITY_FLOAT;
	else if (d < -1e4f)
		outMin = -CMath::INFINITY_FLOAT;
}
std::string CRay2D::toString() const
{
	char str[256];
	std::sprintf(str, "CRay2D(Pos:(%.2f, %.2f) getDir:(%.2f, %.2f))", pos.x(), pos.y(), dir.x, dir.y);
	return str;
}
namespace _2D{
std::ostream &operator <<(std::ostream &o, const CRay2D &ray)
{
	o << ray.toString();
	return o;
}

} //end _2D
} //end GEO
}  //end SMF


