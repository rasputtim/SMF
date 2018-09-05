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
#include "geometry/SMF_2DTriangle.h"
#include "math/all.h"
#include "geometry/all.h"

namespace SMF {
using namespace MATH;
using namespace Util;
namespace GEO{
// compute the area of a triangle using Heron's formula
static float triarea(float a, float b, float c)
{
    float s = (a + b + c)/2.0;
    return sqrt(s*(s-a)*(s-b)*(s-c));
}


CTriangle2D::CTriangle2D(const CPoint2D &a_, const CPoint2D &b_, const CPoint2D &c_)
:a(a_), b(b_), c(c_)
{
}

void CTriangle2D::translate(const CVec2D &offset)
{
	using namespace _2D;
	a += offset;
	b += offset;
	c += offset;
}







CVec3D CTriangle2D::barycentricUVW(const CPoint2D &point) const
{
	 // compute the area of the big triangle
	CVec2D A=a.toVec2DPos();
	CVec2D B=b.toVec2DPos();
	CVec2D C=c.toVec2DPos();
    float a = (B-A).getLenght();
    float b = (C-B).getLenght();
    float c = (A-C).getLenght();
    float totalarea = area();

    // compute the distances from the outer vertices to the inner vertex
    float length0 = (point-A).getLenght();
    float length1 = (point-B).getLenght();
    float length2 = (point-C).getLenght();

    // divide the area of each small triangle by the area of the big triangle
    float u = triarea(b, length1, length2)/totalarea;
    float v = triarea(c, length0, length2)/totalarea;
    float w = triarea(a, length0, length1)/totalarea;
	return CVec3D(u,v,w);
#if 0 // TODO: This version should be more SIMD-friendly, but for some reason, it doesn't return good values for all points inside the triangle.
	CVec2D v0 = b - a;
	CVec2D v1 = c - a;
	CVec2D v2 = point - a;
	float d00 = dot(v0, v0);
	float d01 = dot(v0, v1);
	float d02 = dot(v0, v2);
	float d11 = dot(v1, v1);
	float d12 = dot(v1, v2);
	float denom = 1.f / (d00 * d11 - d01 * d01);
	float v = (d11 * d02 - d01 * d12) * denom;
	float w = (d00 * d12 - d01 * d02) * denom;
	float u = 1.0f - v - w;
	return CVec3D(u, v, w);
#endif
}

CVec2D CTriangle2D::barycentricUV(const CPoint2D &point) const
{
	CVec3D uvw = barycentricUVW(point);
	return CVec2D(uvw.y, uvw.z);
}

bool CTriangle2D::barycentricInsideTriangleboundingAABB(const CVec3D &barycentric)
{
	return barycentric.x >= 0.f && barycentric.y >= 0.f && barycentric.z >= 0.f &&
		CMath::equalsAbs(barycentric.x + barycentric.y + barycentric.z, 1.f);
}

CVec2D CTriangle2D::Point(float u, float v) const
{
	return a + (b-a) * u + (c-a) * v;
}





CVec2D CTriangle2D::Point(const CVec2D &b) const
{
	return Point(b.x, b.y);
}

CVec2D CTriangle2D::centroid() const
{
	return (a + (b + c)) / 3.f;
}

CLineSegment2D CTriangle2D::getMedian(int side)const{
	SMF_ASSERT(0 <= side);
	SMF_ASSERT(side <= 2);
	if (side == 0)
		return CLineSegment2D(edge(side).midPoint(),c);
	else if (side == 1)
		return CLineSegment2D(edge(side).midPoint(),a);
	else if (side == 2)
		return CLineSegment2D(edge(side).midPoint(),b);
	else
		return CLineSegment2D(CVec2D::nan, CVec2D::nan);
}

float CTriangle2D::area() const
{
	return 0.5f * (b-a).cross(c-a);
}

float CTriangle2D::perimeter() const
{
	return a.distance(b) + b.distance(c) + c.distance(a);
}

CLineSegment2D CTriangle2D::edge(int i) const
{
	SMF_ASSERT(0 <= i);
	SMF_ASSERT(i <= 2);
	if (i == 0)
		return CLineSegment2D(a, b);
	else if (i == 1)
		return CLineSegment2D(b, c);
	else if (i == 2)
		return CLineSegment2D(c, a);
	else
		return CLineSegment2D(CVec2D::nan, CVec2D::nan);
}
CPoint2D CTriangle2D::midPointEdge(int i) const{
	SMF_ASSERT(0 <= i);
	SMF_ASSERT(i <= 2);
	CLineSegment2D side=edge(i);
	if (! side.isFinite()) return CVec2D::nan;
	else return side.midPoint();
}
CPoint2D CTriangle2D::vertex(int i) const
{
	SMF_ASSERT(0 <= i);
	SMF_ASSERT(i <= 2);
	if (i == 0)
		return a;
	else if (i == 1)
		return b;
	else if (i == 2)
		return c;
	else
		return CVec2D::nan;
}

CCircle2D CTriangle2D::circunCircle()const{
	CPoint2D center= getCircumcenter();
    float radius= center.distance(a);
	return CCircle2D(center,radius);
}

CCircle2D CTriangle2D::inCircle()const{
	CPoint2D center= getIncenter();
	float radius= 2*area()/perimeter();
	return CCircle2D(center,radius);

}


CPoint2D CTriangle2D::incenter(const CPoint2D &p1, const CPoint2D &p2,const CPoint2D &p3)
{
  float x1,y1,x2,y2,x3,y3;
  x1=p1.x();
  y1=p1.y();
  x2=p2.x();
  y2=p2.y();
  x3=p3.x();
  y3=p3.y();

  float Perim;
  float Side12;
  float Side23;
  float Side31;
  Side12 = p1.distance(p2);//distance(x1,y1,x2,y2);
  Side23 = p2.distance(p3);//distance(x2,y2,x3,y3);
  Side31 = p3.distance(p1);//distance(x3,y3,x1,y1);

  /* Using Heron's S=UR */
  Perim  = 1.0 / (Side12 + Side23 + Side31);
  float Px     = (Side23 * x1 + Side31 * x2 + Side12 * x3) * Perim;
  float Py     = (Side23 * y1 + Side31 * y2 + Side12 * y3) * Perim;
  return CPoint2D(Px,Py);
}

CPoint2D CTriangle2D::getIncenter()const{
	return incenter(a,b,c);
}

CPoint2D CTriangle2D::circumcenter(const CPoint2D &p1, const CPoint2D &p2, const CPoint2D &p3)
{
	CLineSegment2D AB,BC,CA;
	CPoint2D midAB,midBC,midCA, intersection1,intersection2;
	CVec2D normAB, normBC, normCA;

	AB = CLineSegment2D(p1, p2);
	BC = CLineSegment2D(p2, p3);
	//CA = CLineSegment2D(p3, p1);
	midAB = AB.midPoint();
	midBC = BC.midPoint();
	//midCA = CA.midPoint();
	normAB = AB.orthogonal();
	normBC = BC.orthogonal();
	//normCA = CA.orthogonal();
	CLine2D linAB(midAB, normAB);
	CLine2D linBC(midBC, normBC);
	//CLine2D linCA(midCA, normCA);

	CLine2D::intersects_m2(linAB, linBC, intersection1);
	//CLine2D::intersects_m2(linBC, linCA, intersection2);  //just to check
	return intersection1;
}


CPoint2D CTriangle2D::getCircumcenter()const{
	return circumcenter(a,b,c);
}

CPoint2D CTriangle2D::orthoCenter(const CPoint2D &p1, const CPoint2D &p2,const CPoint2D &p3){
	CTriangle2D t(p1,p2,p3);

	return t.getOrthocenter();
}

CPoint2D CTriangle2D::getOrthocenter()const{
	CLine2D l1= edge(side1).orthogonal(c);
	CLine2D l2= edge(side2).orthogonal(a);
	CPoint2D ortho;
	CLine2D::intersects_m2(l1,l2,ortho);
	return ortho;

}
CPoint2D CTriangle2D::extremePoint(const CVec2D &direction) const
{
	CVec2D mostExtreme = CVec2D::nan;
	float mostExtremeDist = -FLT_MAX;
	for(int i = 0; i < 3; ++i)
	{
		CPoint2D pt = vertex(i);
		float d = (direction* pt.toVec2DPos());
		if (d > mostExtremeDist)
		{
			mostExtremeDist = d;
			mostExtreme = pt;
		}
	}
	return mostExtreme;
}


CPolygon2D CTriangle2D::toPolygon() const
{
	CPolygon2D p;
	p.addVertex(a);
	p.addVertex(b);
	p.addVertex(c);
	return p;
}
#if 0
CPolyhedron CTriangle2D::toPolyhedron() const
{
	return toPolygon().toPolyhedron();
}
#endif
CAABBox2D CTriangle2D::boundingAABB() const
{
	CAABBox2D aabb;
	aabb.toNegativeInfinity();
	aabb.enclose(a);
	aabb.enclose(b);
	aabb.enclose(c);
	return aabb;
}

float CTriangle2D::area(const CVec2D &p1, const CVec2D &p2, const CVec2D &p3)
{
	return (p1.x - p2.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p2.y);
}

float CTriangle2D::signedArea(const CVec2D &point, const CVec2D &aa, const CVec2D &bb, const CVec2D &cc)
{
		 // compute the area of the big triangle
    float a = (bb-aa).getLenght();
    float b = (cc-bb).getLenght();
    float c = (aa-cc).getLenght();
    float totalarea = 0.5f * (bb-aa).cross(cc-aa);

    // compute the distances from the outer vertices to the inner vertex
    float length0 = (point-aa).getLenght();
    float length1 = (point-bb).getLenght();
    float length2 = (point-cc).getLenght();

    // divide the area of each small triangle by the area of the big triangle
    float u = triarea(b, length1, length2)/totalarea;
    return u;

}
bool CTriangle2D::isFinite() const
{
	return a.isFinite() && b.isFinite() && c.isFinite();
}

bool CTriangle2D::isDegenerate(float epsilon) const
{
	return isDegenerate(a, b, c, epsilon);
}

bool CTriangle2D::isDegenerate(const CPoint2D &a, const CPoint2D &b, const CPoint2D &c, float epsilon)
{
	return a.compare(b, epsilon) || a.compare(c, epsilon) || b.compare(c, epsilon);
}

bool CTriangle2D::contains(const CPoint2D &point) const
{
	CVec3D br = barycentricUVW(point);
	return br.x >= -1e-3f && br.y >= -1e-3f && br.z >= -1e-3f; // Allow for a small epsilon to properly account for points very near the edges of the triangle.
}

bool CTriangle2D::contains(const CLineSegment2D &lineSegment) const
{
	return contains(lineSegment.point1_) && contains(lineSegment.point1_);
}

bool CTriangle2D::contains(const CTriangle2D &triangle) const
{
	return contains(triangle.a) && contains(triangle.b)
	  && contains(triangle.c);
}

/*
bool CTriangle2D::contains(const CPolygon2D &polygon, float triangleThickness) const
{
	if (polygon.points.size() == 0)
		return false;
	for(int i = 0; i < polygon.points.size(); ++i)
		if (!contains(polygon.points[i], triangleThickness))
			return false;
	return true;
}
*/

float CTriangle2D::distance(const CPoint2D &point) const
{
	return closestPoint(point).distance(point);
}

float CTriangle2D::distance(const CCircle2D &circle) const
{
	return MAX(0.f, distance(circle.getOrigin()) - circle.getRadius());
}


bool CTriangle2D::intersects(const CLineSegment2D &Segment)const
{
  const CTriangle2D Triangle=*this;
  bool Result = Segment.intersects(Triangle.edge(Triangle.side1)) ||
            Segment.intersects(Triangle.edge(Triangle.side2)) ||
            Segment.intersects(Triangle.edge(Triangle.side3)) ||
            contains(Segment[0])        ||
            contains(Segment[1]);
  return Result;
}

bool CTriangle2D::intersects( const CTriangle2D &tri, const CLineSegment2D &line)
{
  int Or1;
  int Or2;
  if (tri.contains(line[0]) || tri.contains(line[1])) return true;

  bool  Result = true;

  Or1 = CPoint2D::orientation(line.point1(), line.point2(), tri.a);
  if (Or1 == 0) return Result;
  Or2 = CPoint2D::orientation(line.point1(), line.point2(), tri.b);

  if (Or2 != Or1) return Result;
  Or2 = CPoint2D::orientation(line.point1(), line.point2(), tri.c);

  return Result = (Or2 != Or1);
}

bool CTriangle2D::intersects(const CPolygon2D &polygon) const{
return polygon.intersects(*this);
}
bool CTriangle2D::intersects(const CLine2D &line) const{
 return line.intersects(edge(side1)) || line.intersects(edge(side2)) || line.intersects(edge(side3));

}
/// [groupSyntax]
bool CTriangle2D::intersection(const CLineSegment2D &Segment, int *count, CPoint2D *I1, CPoint2D *I2) const
{
#define USE_IT(x,value) if(x!=NULL) { *x=value;}
//reset the points
	if(I1!=NULL) {
		if(I1->isFinite()) I1->toInfinite();
	}
	if(I2!=NULL) {
		if(I2->isFinite()) I2->toInfinite();
	}


	int ICnt = 0;
	bool B1,B2,B3;
	CPoint2D p1,p2,p3;
	B1=CLineSegment2D::intersects_m2(edge(side1),Segment,p1);
	B2=CLineSegment2D::intersects_m2(edge(side2),Segment,p2);
	B3=CLineSegment2D::intersects_m2(edge(side3),Segment,p3);
	//check if there are equal points, this mean the ray goes through the vertice of the triangle
	//in this case turn one intersection off
	if (p1==p2) {B2=false;}
	if (p1==p3) {B3= false;}
	if (p2==p3) {B3=false;}
// if intersection is false, then ignore the point
	if(!B1) p1.toInfinite();
	if(!B2) p2.toInfinite();
	if(!B3) p3.toInfinite();

  if (B1) {
    USE_IT(I1,p1)
    ICnt++;
  }
  if (B2) {
    if (ICnt == 1) {
      USE_IT(I2,p2)
      ICnt++;
      }else {
    USE_IT(I1,p2);
    ICnt++;
	}
  }
  if (B3) {
    if (ICnt == 1) {
      USE_IT(I2,p3);
      ICnt++;
     }else {
		USE_IT(I1,p3);
		ICnt++;
	}
  }
   USE_IT(count,ICnt);
  if (ICnt > 0 ) return true;
  else return false;
#undef USE_IT(x,value)
}

bool CTriangle2D::intersection(const CLine2D &l, int *ICnt, CPoint2D *I1, CPoint2D *I2) const
{
//reset the points
	if(I1!=NULL) {
		if(I1->isFinite()) I1->toInfinite();
	}
	if(I2!=NULL) {
		if(I2->isFinite()) I2->toInfinite();
	}

	//create a sgment on the line that is on the same region of the triangle
	float xmin = MIN(MIN(a.x_,b.x_),c.x_);
	float xmax = MAX(MAX(a.x_,b.x_),c.x_);
	float ymin = MIN(MIN(a.y_,b.y_),c.y_);
	float ymax = MAX(MAX(a.y_,b.y_),c.y_);
	//verificar se a linha é vertical
	CPoint2D minPoint,maxPoint;
	if(l.isVertical()) {
		float x=-(l.c()/l.a());
		minPoint.set(x,ymin);
		maxPoint.set(x,ymax);
	}else{
	minPoint = l.getPoint(xmin);
	if (!MATH::isFinite(minPoint.y_)) minPoint.y_= ymin;
	maxPoint = l.getPoint(xmax);
	if (!MATH::isFinite(maxPoint.y_)) maxPoint.y_= ymax;
	}
	const CLineSegment2D seg(minPoint,maxPoint);
	//calcule intersection of this segment with the triangle
	bool t = intersection(seg, ICnt, I1,I2);

	return t;
}

bool CTriangle2D::intersection(const CRay2D &ray, int *ICnt, CPoint2D *I1, CPoint2D *I2) const{
#define USE_IT(x,value) if(x!=NULL) { *x=value;}
//reset the points
	if(I1!=NULL) {
		if(I1->isFinite()) I1->toInfinite();
	}
	if(I2!=NULL) {
		if(I2->isFinite()) I2->toInfinite();
	}
	int Count=0;
	bool B1,B2,B3;
	int Result1=0;
	int Result2=0;
	int Result3=0;
	CPoint2D p1,p2,p3;
	B1=CLineSegment2D::intersects_m2(edge(side1),ray,p1,&Result1);
	B2=CLineSegment2D::intersects_m2(edge(side2),ray,p2,&Result2);
	B3=CLineSegment2D::intersects_m2(edge(side3),ray,p3,&Result3);
	//check if there are equal points, this mean the ray goes through the vertice of the triangle
	//in this case turn one intersection off
	if (p1==p2) {B2=false;}
	if (p1==p3) {B3= false;}
	if (p2==p3) {B3=false;}
	// if intersection is false, then ignore the point
	if(!B1) p1.toInfinite();
	if(!B2) p2.toInfinite();
	if(!B3) p3.toInfinite();

	if(B1){
		if(p1.isFinite()) { USE_IT(I1,p1); Count++;}

	}
	if (B2){
		if(!p1.isFinite()) {
			*I1=p2;
			Count++;
		}else{
			if(p2.isFinite()) {USE_IT(I2,p2); Count++;}
		}
	}
	if (B3){
		if(!p1.isFinite()) {
			USE_IT(I1,p3);
			Count++;
		}else{
			if(p2.isFinite()) {USE_IT(I2,p3); Count++;}

		}


	}
USE_IT(ICnt,Count)
return Count>0;
#undef USE_IT(x,value)
}

bool CTriangle2D::intersects(const CCircle2D &circle, CPoint2D *closestPointOnTriangle) const
{
	CPoint2D pt = closestPoint(circle.getOrigin());

	if (closestPointOnTriangle)
		*closestPointOnTriangle = pt;

	return pt.distanceSqr(circle.getOrigin()) <= circle.getRadius() * circle.getRadius();
}

bool CTriangle2D::intersects(const CCircle2D &circle) const
{
	return intersects(circle, 0);
}



bool CTriangle2D::intersects(const CTriangle2D &Triangle2)const
{
	const CTriangle2D &Triangle1=*this;
	int i;
	bool Result = false;
  for (i = 0; i < 3;i++ )
  {

    if (CMath::equals(Triangle2.distance(Triangle1[i]),CMath::ZERO) ||
        CMath::equals(Triangle1.distance(Triangle2[i]),CMath::ZERO))
	{
      Result = true;
      return Result;
	}
  }
  return Result;
}

static void FindIntersectingLineSegments(const CTriangle2D &t, float da, float db, float dc, CLineSegment2D &l1, CLineSegment2D &l2)
{
	if (da*db > 0.f)
	{
		l1 = CLineSegment2D(t.a, t.c);
		l2 = CLineSegment2D(t.b, t.c);
	}
	else if (db*dc > 0.f)
	{
		l1 = CLineSegment2D(t.a, t.b);
		l2 = CLineSegment2D(t.a, t.c);
	}
	else
	{
		l1 = CLineSegment2D(t.a, t.b);
		l2 = CLineSegment2D(t.b, t.c);
	}
}

bool CTriangle2D::intersects(const CAABBox2D &aabb) const
{

  bool Result = aabb.intersects(edge(side1)) ||
                aabb.intersects(edge(side2)) ||
                aabb.intersects(edge(side3));

  return Result;
}


#if 0
/// [groupSyntax]
bool intersection(const CAABBox2D &aabb,  int *ICnt=NULL, CPoint2D *I1=NULL, CPoint2D *I2=NULL) const{


}
bool CTriangle2D::intersects(const COBBox &obb) const
{
	return obb.intersects(*this);
}

bool CTriangle2D::intersects(const CPolyhedron &polyhedron) const
{
	return polyhedron.intersects(*this);
}
#endif
void CTriangle2D::projectToAxis(const CVec2D &axis, float &dMin, float &dMax) const
{
	dMin = dMax = (axis* a.toVec2DPos());
	float t = (axis* b.toVec2DPos());
	dMin = MIN(t, dMin);
	dMax = MAX(t, dMax);
	t = (axis* c.toVec2DPos());
	dMin = MIN(t, dMin);
	dMax = MAX(t, dMax);
}

/// [groupSyntax]
CPoint2D CTriangle2D::closestPoint(const CPoint2D &p) const
{
	/** The code for CTriangle2D-CVec2D test is from Christer Ericson's Real-Time Collision Detection, pp. 141-142. */

	// Check if P is in vertex region outside A.
	CVec2D ab = b - a;
	CVec2D ac = c - a;
	CVec2D ap = p - a;
	float d1 = (ab* ap);
	float d2 = (ac* ap);
	if (d1 <= 0.f && d2 <= 0.f)
		return a; // Barycentric coordinates are (1,0,0).

	// Check if P is in vertex region outside B.
	CVec2D bp = p - b;
	float d3 = (ab* bp);
	float d4 = (ac* bp);
	if (d3 >= 0.f && d4 <= d3)
		return b; // Barycentric coordinates are (0,1,0).

	// Check if P is in edge region of AB, and if so, return the projection of P onto AB.
	float vc = d1*d4 - d3*d2;
	if (vc <= 0.f && d1 >= 0.f && d3 <= 0.f)
	{
		float v = d1 / (d1 - d3);
		return a + v * ab; // The barycentric coordinates are (1-v, v, 0).
	}

	// Check if P is in vertex region outside C.
	CVec2D cp = p - c;
	float d5 = (ab* cp);
	float d6 = (ac* cp);
	if (d6 >= 0.f && d5 <= d6)
		return c; // The barycentric coordinates are (0,0,1).

	// Check if P is in edge region of AC, and if so, return the projection of P onto AC.
	float vb = d5*d2 - d1*d6;
	if (vb <= 0.f && d2 >= 0.f && d6 <= 0.f)
	{
		float w = d2 / (d2 - d6);
		return a + w * ac; // The barycentric coordinates are (1-w, 0, w).
	}

	// Check if P is in edge region of BC, and if so, return the projection of P onto BC.
	float va = d3*d6 - d5*d4;
	if (va <= 0.f && d4 - d3 >= 0.f && d5 - d6 >= 0.f)
	{
		float w = (d4 - d3) / (d4 - d3 + d5 - d6);
		return b + w * (c - b); // The barycentric coordinates are (0, 1-w, w).
	}

	// P must be inside the face region. Compute the closest point through its barycentric coordinates (u,v,w).
	float denom = 1.f / (va + vb + vc);
	float v = vb * denom;
	float w = vc * denom;
	return a + ab * v + ac * w; //=u*a+v*b+w*c,u=va*denom = 1.0 f-v-w
}

CPoint2D CTriangle2D::closestPoint_m2(const CPoint2D &point)const
{
	const CPoint2D &p1=a;
	const CPoint2D &p2=b;
	const CPoint2D &p3=c;
  if (CPoint2D::orientation(p1,p2,point) != CPoint2D::orientation(p1,p2,p3)){

    return edge(0).closestPoint(point);

  }

  if (CPoint2D::orientation(p2,p3,point) != CPoint2D::orientation(p2,p3,p1) ){

    return edge(1).closestPoint(point);

  }

  if (CPoint2D::orientation(p3,p1,point) != CPoint2D::orientation(p3,p1,p2)){

    return edge(2).closestPoint(point);

  }

return point;
}

CPoint2D CTriangle2D::randomPointInside(CRandomLCG &rng) const
{
	float epsilon =CMath::EPSILON_SuperLow;
	///\todo rng.getFloat() returns [0,1[, but to be completely uniform, we'd need [0,1] here.
	float s = rng.getFloat(epsilon, 1.f - epsilon);//1e-2f, 1.f - 1e-2f);
	float t = rng.getFloat(epsilon, 1.f - epsilon);//1e-2f, 1.f - 1e-2f
	if (s + t >= 1.f)
	{
		s = 1.f - s;
		t = 1.f - t;
	}
#ifdef MATH_ASSERT_CORRECTNESS
	CVec2D pt = Point(s, t);
	CVec2D uv = barycentricUV(pt);
	SMF_ASSERT(uv.x >= 0.f);
	SMF_ASSERT(uv.y >= 0.f);
	SMF_ASSERT(uv.x + uv.y <= 1.f);
	CVec2D uvw = barycentricUVW(pt);
	SMF_ASSERT(uvw.x >= 0.f);
	SMF_ASSERT(uvw.y >= 0.f);
	SMF_ASSERT(uvw.z >= 0.f);
	SMF_ASSERT(equalsAbs(uvw.x + uvw.y + uvw.z, 1.f));
#endif
	return Point(s, t);
}

CPoint2D CTriangle2D::randomVertex(CRandomLCG &rng) const
{
	return vertex(rng.getInt(0, 2));
}

CPoint2D CTriangle2D::randomPointOnEdge(CRandomLCG &rng) const
{
	SMF_ASSERT(!isDegenerate());
	float ab = a.distance(b);
	float bc = b.distance(c);
	float ca = c.distance(a);
	float r = rng.getFloat(0, ab + bc + ca);
	if (r < ab)
		return a + (b-a) * r / ab;
	r -= ab;
	if (r < bc)
		return b + (c-b) * r / bc;
	r -= bc;
	return c + (a-c) * r / ca;
}


std::string CTriangle2D::toString() const
{
	char str[256];
	std::sprintf(str, "CTriangle2D(a:(%.2f, %.2f) b:(%.2f, %.2f) c:(%.2f, %.2f))",
		a.x_, a.y_, b.x_, b.y_, c.x_, c.y_);
	return str;
}

std::ostream &operator <<(std::ostream &o, const CTriangle2D &triangle)
{
	o << triangle.toString();
	return o;
}





} //end GEO
}  //end SMF

