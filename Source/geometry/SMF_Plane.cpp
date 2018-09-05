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
#include "geometry/SMF_Plane.h"
#include "geometry/all.h"
#include "util/SMF_StringUtils.h"


//#pragma hdrstop
namespace SMF {
namespace GEO{
using namespace MATH;
CPlane plane_origin( 0.0f, 0.0f, 0.0f, 0.0f );

/*
================
CPlane::Type
================
*/
int CPlane::getType() const {
	if ( getNormal()[0] == 0.0f ) {
		if ( getNormal()[1] == 0.0f ) {
			return getNormal()[2] > 0.0f ? PLANETYPE_Z : PLANETYPE_NEGZ;
		}
		else if ( getNormal()[2] == 0.0f ) {
			return getNormal()[1] > 0.0f ? PLANETYPE_Y : PLANETYPE_NEGY;
		}
		else {
			return PLANETYPE_ZEROX;
		}
	}
	else if ( getNormal()[1] == 0.0f ) {
		if ( getNormal()[2] == 0.0f ) {
			return getNormal()[0] > 0.0f ? PLANETYPE_X : PLANETYPE_NEGX;
		}
		else {
			return PLANETYPE_ZEROY;
		}
	}
	else if ( getNormal()[2] == 0.0f ) {
		return PLANETYPE_ZEROZ;
	}
	else {
		return PLANETYPE_NONAXIAL;
	}
}

/*
================
CPlane::heightFit
================
*/
bool CPlane::heightFit( const CVec3D *points, const int numPoints ) {
	int i;
	float sumXX = 0.0f, sumXY = 0.0f, sumXZ = 0.0f;
	float sumYY = 0.0f, sumYZ = 0.0f;
	CVec3D sum, average, dir;

	if ( numPoints == 1 ) {
		a = 0.0f;
		b = 0.0f;
		c = 1.0f;
		d = -points[0].z;
		return true;
	}
	if ( numPoints == 2 ) {
		dir = points[1] - points[0];
		getNormal() = dir.cross( CVec3D( 0, 0, 1 ) ).cross( dir );
		normalize();
		d = -( getNormal() * points[0] );
		return true;
	}

	sum.toZero();
	for ( i = 0; i < numPoints; i++) {
		sum += points[i];
	}
	average = sum / numPoints;

	for ( i = 0; i < numPoints; i++ ) {
		dir = points[i] - average;
		sumXX += dir.x * dir.x;
		sumXY += dir.x * dir.y;
		sumXZ += dir.x * dir.z;
		sumYY += dir.y * dir.y;
		sumYZ += dir.y * dir.z;
	}

	CMat2D m( sumXX, sumXY, sumXY, sumYY );
	if ( !m.inverseSelf() ) {
		return false;
	}

	a = - sumXZ * m[0][0] - sumYZ * m[0][1];
	b = - sumXZ * m[1][0] - sumYZ * m[1][1];
	c = 1.0f;
	normalize();
	d = -( a * average.x + b * average.y + c * average.z );
	return true;
}

/*
================
CPlane::planeIntersection
================
*/
bool CPlane::planeIntersection( const CPlane &plane, CVec3D &start, CVec3D &dir ) const {
	double n00, n01, n11, det, invDet, f0, f1;

	n00 = getNormal().getLengthSqr();
	n01 = getNormal() * plane.getNormal();
	n11 = plane.getNormal().getLengthSqr();
	det = n00 * n11 - n01 * n01;

	if ( CMath::fabs(det) < 1e-6f ) {
		return false;
	}

	invDet = 1.0f / det;
	f0 = ( n01 * plane.d - n11 * d ) * invDet;
	f1 = ( n01 * d - n00 * plane.d ) * invDet;

	dir = getNormal().cross( plane.getNormal() );
	start = f0 * getNormal() + f1 * plane.getNormal();
	return true;
}

/*
=============
CPlane::toString
=============
*/
const char *CPlane::toString( int precision ) const {
	return CMyString::floatArraytoString( toFloatPtr(), getDimension(), precision );
}



float CPlane::distance(const CVec3D &point) const
{
	return CMath::fabs(signedDistance(point));
}

float CPlane::distance(const CLineSegment &lineSegment) const
{
	return lineSegment.distance(*this);
}

float CPlane::distance(const CSphere &sphere) const
{
	return MAX(0.f, distance(sphere.getOrigin()) - sphere.getRadius());
}

float CPlane::signedDistance(const CVec3D &point) const
{
	return getNormal()*(point) - d;
}

template<typename T>
float  plane_SignedDistance(const CPlane &plane, const T &object)
{
	float pMin, pMax;
	 SMF_ASSERT(plane.getNormal().isNormalized());
	object.projectToAxis(plane.getNormal(), pMin, pMax);
	pMin -= plane.d;
	pMax -= plane.d;
	if (pMin * pMax <= 0.f)
		return 0.f;
	return CMath::fabs(pMin) < CMath::fabs(pMax) ? pMin : pMax;
}

float CPlane::signedDistance(const CAABBox &aabb) const { return  plane_SignedDistance(*this, aabb); }
float CPlane::signedDistance(const COBBox &obb) const { return  plane_SignedDistance(*this, obb); }
float CPlane::signedDistance(const CLine &line) const { return  plane_SignedDistance(*this, line); }
float CPlane::signedDistance(const CLineSegment &lineSegment) const { return  plane_SignedDistance(*this, lineSegment); }
float CPlane::signedDistance(const CRay &ray) const { return  plane_SignedDistance(*this, ray); }
//S float CPlane::signedDistance(const CPolygon &polygon) const { return  plane_SignedDistance(*this, polygon); }
//S float CPlane::signedDistance(const CPolyhedron &polyhedron) const { return  plane_SignedDistance(*this, polyhedron); }
float CPlane::signedDistance(const CSphere &sphere) const { return  plane_SignedDistance(*this, sphere); }
float CPlane::signedDistance(const CTriangle &triangle) const { return  plane_SignedDistance(*this, triangle); }

CVec3D CPlane::project(const CVec3D &point) const
{
	CVec3D projected = point - ((getNormal()* point) - d) * getNormal();
	//s SMF_ASSERT(projected.compare(OrthoProjection().MulPos(point)));
	return projected;
}

CLineSegment CPlane::project(const CLineSegment &lineSegment) const
{
	return CLineSegment(project(lineSegment.begin), project(lineSegment.end));
}

CLine CPlane::project(const CLine &line, bool *nonDegenerate) const
{
	CLine l;
	l.pos = project(line.pos);
	l.dir = l.dir - l.dir.projectToNorm(getNormal());
	float len = l.dir.getLenght();
	l.dir.toNormal();
	if (nonDegenerate)
		*nonDegenerate = (len > 0.f);
	return l;
}

CRay CPlane::project(const CRay &ray, bool *nonDegenerate) const
{
	CRay r;
	r.pos = project(ray.pos);
	r.dir = r.dir - r.dir.projectToNorm(getNormal());
	float len = r.dir.getLenght();
	r.dir.toNormal();
	if (nonDegenerate)
		*nonDegenerate = (len > 0.f);
	return r;
}

CTriangle CPlane::project(const CTriangle &triangle) const
{
	CTriangle t;
	t.a = project(triangle.a);
	t.b = project(triangle.b);
	t.c = project(triangle.c);
	return t;
}

CPolygon CPlane::project(const CPolygon &polygon) const
{
	CPolygon p;
	for(size_t i = 0; i < polygon.p.size(); ++i)
		p.p.push_back(project(polygon.p[i]));

	return p;
}

bool CPlane::intersects(const CPlane &plane, CLine *outLine) const
{
	CVec3D perp = getNormal().perpendicular(plane.getNormal());//CVec3D::perpendicular cross(getNormal(), plane.getNormal());

	CMat3D m;
	m.setRow(0, getNormal());
	m.setRow(1, plane.getNormal());
	m.setRow(2, perp); // This is arbitrarily chosen, to produce m invertible.
	CVec3D intersectionPos;
	bool success = m.solveAxb(CVec3D(d, plane.d, 0.f),intersectionPos);
	if (!success) // Inverse failed, so the planes must be parallel.
	{
		float normalDir = (getNormal()*plane.getNormal());
		if ((normalDir > 0.f && CMath::equalsAbs(d, plane.d)) || (normalDir < 0.f && CMath::equalsAbs(d, -plane.d)))
		{
			if (outLine)
				*outLine = CLine(getNormal()*d, plane.getNormal().perpendicular());
			return true;
		}
		else
			return false;
	}
	if (outLine)
		*outLine = CLine(intersectionPos, perp.normalized());
	return true;
}

bool CPlane::intersects(const CPlane &plane, const CPlane &plane2, CLine *outLine, CVec3D *outPoint) const
{
	CLine dummy;
	if (!outLine)
		outLine = &dummy;

	// First check all planes for parallel pairs.
	if (this->isParallel(plane) || this->isParallel(plane2))
	{
		if (CMath::equalsAbs(d, plane.d) || CMath::equalsAbs(d, plane2.d))
		{
			bool intersect = plane.intersects(plane2, outLine);
			if (intersect && outPoint)
				*outPoint = outLine->getPoint(0);
			return intersect;
		}
		else
			return false;
	}
	if (plane.isParallel(plane2))
	{
		if (CMath::equalsAbs(plane.d, plane2.d))
		{
			bool intersect = this->intersects(plane, outLine);
			if (intersect && outPoint)
				*outPoint = outLine->getPoint(0);
			return intersect;
		}
		else
			return false;
	}

	// All planes point to different directions.
	CMat3D m;
	m.setRow(0, getNormal());
	m.setRow(1, plane.getNormal());
	m.setRow(2, plane2.getNormal());
	CVec3D intersectionPos;
	bool success = m.solveAxb(CVec3D(d, plane.d, plane2.d), intersectionPos);
	if (!success)
		return false;
	if (outPoint)
		*outPoint = intersectionPos;
	return true;
}

bool CPlane::intersects(const CPolygon &polygon) const
{
	return polygon.intersects(*this);
}

#if 0
bool Plane::intersectLinePlane(const CVec3D &p, const CVec3D &n, const CVec3D &a, const CVec3D &d, float &t)
{
	/* The set of points x lying on a plane is defined by the equation

		(x - p)*n == 0, where p is a point on the plane, and n is the plane getNormal().

	The set of points x on a line is constructed explicitly by a single parameter t by

		x = a + t*d, where a is a point on the line, and d is the direction vector of the line.

	To solve the intersection of these two objects, substitute the second equation to the first above,
	and we get

		  (a + t*d - p)*n == 0, or
		t*(d*n) + (a-p)*n == 0, or
	                    t == (p-a)*n / (d*n), assuming that d*n != 0.

	If d*n == 0, then the line is parallel to the plane, and either no intersection occurs, or the whole line
	is embedded on the plane, and infinitely many intersections occur. */

	float denom = dot(d, n);
	if (equalsAbs(denom, 0.f))
	{
		t = 0.f;
		float f = dot(a-p, n);
		bool b = equalsAbs(dot(a-p, n), 0.f);
		return equalsAbs(dot(a-p, n), 0.f); // If (a-p)*n == 0, then then above equation holds for all t, and return true.
	}
	else
	{
		// Compute the distance from the line starting point to the point of intersection.
		t = dot(p - a, n) / denom;
		return true;
	}
}
#endif

bool CPlane::intersectLinePlane(const CVec3D &planeNormal, float planeD, const CVec3D &linePos, const CVec3D &lineDir, float &t)
{
	/* The set of points x lying on a plane is defined by the equation

		<planeNormal, x> = planeD.

	The set of points x on a line is constructed explicitly by a single parameter t by

		x = linePos + t*lineDir.

	To solve the intersection of these two objects, substitute the second equation to the first above,
	and we get

	                     <planeNormal, linePos + t*lineDir> == planeD, or
	    <planeNormal, linePos> + t * <planeNormal, lineDir> == planeD, or
	                                                      t == (planeD - <planeNormal, linePos>) / <planeNormal, lineDir>,
	
	                                                           assuming that <planeNormal, lineDir> != 0.

	If <planeNormal, lineDir> == 0, then the line is parallel to the plane, and either no intersection occurs, or the whole line
	is embedded on the plane, and infinitely many intersections occur. */

	float denom = planeNormal * lineDir;
	if (CMath::equalsAbs(denom, 0.f))
	{
		t = 0.f;
		return CMath::equalsAbs((planeNormal* linePos), planeD, 1e-2f);
	}
	else
	{
		// Compute the distance from the line starting point to the point of intersection.
		t = (planeD - (planeNormal* linePos)) / denom;
		return true;
	}
}


bool CPlane::intersects(const CRay &ray, float *d) const
{
	float t;
	bool success = intersectLinePlane(getNormal(), this->d, ray.pos, ray.dir, t);
	if (d)
		*d = t;
	return success && t >= 0.f;
}

bool CPlane::intersects(const CLine &line, float *d) const
{
	float t;
	bool intersects = intersectLinePlane(getNormal(), this->d, line.pos, line.dir, t);
	if (d)
		*d = t;
	return intersects;
}

bool CPlane::intersects(const CLineSegment &lineSegment, float *d) const
{
	float t;
	bool success = intersectLinePlane(getNormal(), this->d, lineSegment.begin, lineSegment.getDir(), t);
	const float lineSegmentLength = lineSegment.getLenght();
	if (d)
		*d = t / lineSegmentLength;
	return success && t >= 0.f && t <= lineSegmentLength;
}

bool CPlane::intersects(const CSphere &sphere) const
{
	return distance(sphere.getOrigin()) <= sphere.getRadius();
}


/// The Plane-CAABBox intersection is implemented according to Christer Ericson's Real-Time Collision Detection, p.164. [groupSyntax]
bool CPlane::intersects(const CAABBox &aabb) const
{
	CVec3D c = aabb.getCenter();
	CVec3D e = aabb.getHalfDiagonal();

	// Compute the projection interval radius of the CAABBox onto L(t) = aabb.center + t * plane.getNormal();
	float r = e[0]*CMath::fabs(getNormal()[0]) + e[1]*CMath::fabs(getNormal()[1]) + e[2]*CMath::fabs(getNormal()[2]);
	// Compute the distance of the box center from plane.
	float s =getNormal()* c - d;
	return CMath::fabs(s) <= r;
}

bool CPlane::intersects(const COBBox &obb) const
{
	return obb.intersects(*this);
}

bool CPlane::intersects(const CTriangle &triangle) const
{
	float a = signedDistance(triangle.a);
	float b = signedDistance(triangle.b);
	float c = signedDistance(triangle.c);
	return (a*b <= 0.f || a*c <= 0.f);
}



bool CPlane::intersects(const CPolyhedron &polyhedron) const
{
	if (polyhedron.numVertices() == 0)
		return false;
	bool sign = isOnPositiveSide(polyhedron.vertex(0));
	for(int i = 1; i < polyhedron.numVertices(); ++i)
		if (sign != isOnPositiveSide(polyhedron.vertex(i)))
			return true;
	return false;
}

int CPlane::intersects(const CCircle &circle, CVec3D *pt1, CVec3D *pt2) const
{
	CLine line;
	bool planeIntersects = intersects(circle.containingPlane(), &line);
	if (!planeIntersects)
		return false;

	// Offset both line and circle position so the circle origin is at center.
	line.pos -= circle.pos;

	float a = 1.f;
	float b = 2.f * (line.pos* line.dir);
	float c = line.pos.getLengthSqr() - circle.r * circle.r;
	float r1, r2;
	int numRoots = CPolynomial::solveQuadratic(a, b, c, r1, r2);
	if (numRoots >= 1 && pt1)
		*pt1 = circle.pos + line.getPoint(r1);
	if (numRoots >= 2 && pt2)
		*pt2 = circle.pos + line.getPoint(r2);
	return numRoots;
}

int CPlane::intersects(const CCircle &circle) const
{
	return intersects(circle, 0, 0);
}

bool CPlane::isParallel(const CPlane &plane, float epsilon) const
{
	return getNormal().compare(plane.getNormal(), epsilon);
}

bool CPlane::isInPositiveDirection(const CVec3D &directionVector) const
{
	return getNormal()*(directionVector) >= 0.f;
}

bool CPlane::isOnPositiveSide(const CVec3D &point) const
{
	return signedDistance(point) >= 0.f;
}

int CPlane::examineSide(const CTriangle &triangle) const
{
	float a = signedDistance(triangle.a);
	float b = signedDistance(triangle.b);
	float c = signedDistance(triangle.c);
	const float epsilon =CMath::EPSILON_SuperLow; // Allow a small epsilon amount for tests for floating point inaccuracies.
	if (a >= -epsilon && b >= -epsilon && c >= -epsilon)
		return 1;
	if (a <= epsilon && b <= epsilon && c <= epsilon)
		return -1;
	return 0;
}

bool CPlane::areOnSameSide(const CVec3D &p1, const CVec3D &p2) const
{
	return signedDistance(p1) * signedDistance(p2) >= 0.f;
}

CPlane::CPlane(const CRay &ray, const CVec3D &normal)
{
	CVec3D perpNormal = normal - normal.projectToNorm(ray.dir);
	set(ray.pos, perpNormal.normalized());
}

CPlane::CPlane(const CLine &line, const CVec3D &normal)
{
	CVec3D perpNormal = normal - normal.projectToNorm(line.dir);
	set(line.pos, perpNormal.normalized());
}

CPlane::CPlane(const CLineSegment &lineSegment, const CVec3D &normal)
{
	CVec3D perpNormal = normal - normal.projectTo(lineSegment.end - lineSegment.begin);
	set(lineSegment.begin, perpNormal.normalized());
}
} // END MATH
} //end SMF
