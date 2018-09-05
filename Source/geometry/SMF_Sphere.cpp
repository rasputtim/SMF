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
#include "geometry/SMF_Sphere.h"
#include "geometry/all.h"
namespace SMF{
using namespace MATH;
namespace GEO {


const CSphere CSphere::zero( CVec3D::origin, 0.0f );


/*
================
CSphere::planeDistance
================
*/
float CSphere::planeDistance( const CPlane &plane ) const {
	float d;

	d = plane.getDistance( origin );
	if ( d > radius ) {
		return d - radius;
	}
	if ( d < -radius ) {
		return d + radius;
	}
	return 0.0f;
}

/*
================
CSphere::planeSide
================
*/
int CSphere::planeSide( const CPlane &plane, const float epsilon ) const {
	float d;

	d = plane.getDistance( origin );
	if ( d > radius + epsilon ) {
		return PLANESIDE_FRONT;
	}
	if ( d < -radius - epsilon ) {
		return PLANESIDE_BACK;
	}
	return PLANESIDE_CROSS;
}

/*
============
CSphere::lineIntersection

  Returns true if the line intersects the sphere between the start and end point.
============
*/
bool CSphere::lineIntersection( const CVec3D &start, const CVec3D &end ) const {
	CVec3D r, s, e;
	float a;

	s = start - origin;
	e = end - origin;
	r = e - s;
	a = -s * r;
	if ( a <= 0 ) {
		return ( s * s < radius * radius );
	}
	else if ( a >= r * r ) {
		return ( e * e < radius * radius );
	}
	else {
		r = s + ( a / ( r * r ) ) * r;
		return ( r * r < radius * radius );
	}
}

/*
============
CSphere::rayIntersection

  Returns true if the ray intersects the sphere.
  The ray can intersect the sphere in both directions from the start point.
  If start is inside the sphere then scale1 < 0 and scale2 > 0.
============
*/
bool CSphere::rayIntersection( const CVec3D &start, const CVec3D &dir, float &scale1, float &scale2 ) const {
	double a, b, c, d, sqrtd;
	CVec3D p;

	p = start - origin;
	a = dir * dir;
	b = dir * p;
	c = p * p - radius * radius;
	d = b * b - c * a;

	if ( d < 0.0f ) {
		return false;
	}

	sqrtd = CMath::sqrt( d );
	a = 1.0f / a;

	scale1 = ( -b + sqrtd ) * a;
	scale2 = ( -b - sqrtd ) * a;

	return true;
}

/*
============
CSphere::fromPoints

  Tight sphere for a point set.
============
*/
void CSphere::fromPoints( const CVec3D *points, const int numPoints ) {
	int i;
	float radiusSqr, dist;
	CVec3D mins, maxs;

	SIMDProcessor->minMax( mins, maxs, points, numPoints );

	origin = ( mins + maxs ) * 0.5f;

	radiusSqr = 0.0f;
	for ( i = 0; i < numPoints; i++ ) {
		dist = ( points[i] - origin ).getLengthSqr();
		if ( dist > radiusSqr ) {
			radiusSqr = dist;
		}
	}
	radius = CMath::sqrt( radiusSqr );
}
bool CSphere::isFinite() const
{
	return origin.isFinite() && MATH::isFinite(radius);
}

bool CSphere::isDegenerate() const
{
	return radius < 0.f;
}

float CSphere::distance(const CVec3D &point) const
{
	return MAX(0.f, origin.distance(point) - radius);
}

float CSphere::distance(const CSphere &sphere) const
{
	return MAX(0.f, origin.distance(sphere.origin) - radius - sphere.radius);
}

float CSphere::distance(const CAABBox &aabb) const
{
	return aabb.distance(*this);
}

float CSphere::distance(const COBBox &obb) const
{
	return obb.distance(*this);
}

float CSphere::distance(const CPlane &plane) const
{
	return plane.distance(*this);
}

float CSphere::distance(const CTriangle &triangle) const
{
	return triangle.distance(*this);
}

float CSphere::distance(const CRay &ray) const
{
	return ray.distance(*this);
}

float CSphere::distance(const CLineSegment &lineSegment) const
{
	return lineSegment.distance(*this);
}

float CSphere::distance(const CLine &line) const
{
	return line.distance(*this);
}



bool CSphere::intersects(const CSphere &sphere) const
{
	return (origin- sphere.origin).getLengthSqr() <= (radius + sphere.getRadius()) * (radius + sphere.getRadius());
}



int CSphere::intersectLine(const CVec3D &linePos, const CVec3D &lineDir, const CVec3D &sphereCenter,
                          float sphereRadius, float &t1, float &t2)
{
	SMF_ASSERT(lineDir.isNormalized());
	SMF_ASSERT(sphereRadius >= 0.f);

	/* A line is represented explicitly by the set { linePos + t * lineDir }, where t is an arbitrary float.
	  A sphere is represented implictly by the set of vectors that satisfy ||v - sphereCenter|| == sphereRadius.
	  To solve which points on the line are also points on the sphere, substitute v <- linePos + t * lineDir
	  to obtain:

	    || linePos + t * lineDir - sphereCenter || == sphereRadius, and squaring both sides we get
	    || linePos + t * lineDir - sphereCenter ||^2 == sphereRadius^2, or rearranging:
	    || (linePos - sphereCenter) + t * lineDir ||^2 == sphereRadius^2. */

	// This equation represents the set of points which lie both on the line and the sphere. There is only one
	// unknown variable, t, for which we solve to get the actual points of intersection.

	// Compute variables from the above equation:
	const CVec3D a = linePos - sphereCenter;
	const float radSq = sphereRadius * sphereRadius;

	/* so now the equation looks like

	    || a + t * lineDir ||^2 == radSq.

	  Since ||x||^2 == <x,x> (i.e. the square of a vector norm equals the dot product with itself), we get
	
	    <a + t * lineDir, a + t * lineDir> == radSq,
	
	  and using the identity <a+b, a+b> == <a,a> + 2*<a,b> + <b,b> (which holds for dot product when a and b are reals),
	  we have

	    <a,a> + 2 * <a, t * lineDir> + <t * lineDir, t * lineDir> == radSq, or		
	    <a,a> - radSq + 2 * <a, lineDir> * t + <lineDir, lineDir> * t^2 == 0, or

	    C + Bt +at^2 == 0, where

	    C = <a,a> - radSq,
	    B = 2 * <a, lineDir>, and
	    A = <lineDir, lineDir> == 1, since we assumed lineDir is normalized. */

	// Warning! If dot(a,a) is large (distance between line pos and sphere center) and sphere radius very small,
	// catastrophic cancellation can occur here!
	const float C = (a*a) - radSq;
	const float B = 2.f * (a* lineDir);

	/* The equation A + Bt + Ct^2 == 0 is a second degree equation on t, which is easily solvable using the
	  known formula, and we obtain

	    t = [-B +/- sqrt(B^2 - 4AC)] / 2A. */

	float D = B*B - 4.f * C; // D = B^2 - 4AC.
	if (D < 0.f) // There is no solution to the square root, so the ray doesn't intersect the sphere.
	{
		// Output a degenerate enter-exit range so that batch processing code may use min of t1's and max of t2's to
		// compute the nearest enter and farthest exit without requiring branching on the return value of this function.
		t1 = CMath::INFINITY_FLOAT;
		t2 = -CMath::INFINITY_FLOAT;
		return 0;
	}

	if (D < 1e-4f) // The expression inside sqrt is ~ 0. The line is tangent to the sphere, and we have one solution.
	{
		t1 = t2 = -B * 0.5f;
		return 1;
	}

	// The sqrt expression is strictly positive, so we get two different solutions for t.
	D = CMath::sqrt(D);
	t1 = (-B - D) * 0.5f;
	t2 = (-B + D) * 0.5f;
	return 2;
}

int CSphere::intersects(const CRay &ray, CVec3D *intersectionPoint, CVec3D *intersectionNormal, float *d, float *d2) const
{
	float t1, t2;
	int numIntersections = intersectLine(ray.pos, ray.dir, origin, radius, t1, t2);

	// If the line of this ray intersected in two places, but the first intersection was "behind" this ray,
	// handle the second point of intersection instead. This case occurs when the origin of the ray is inside
	// the Sphere.
	if (t1 < 0.f && numIntersections == 2)
		t1 = t2;

	if (t1 < 0.f)
		return 0; // The intersection position is on the negative direction of the ray.

	CVec3D hitPoint = ray.pos + t1 * ray.dir;
	if (intersectionPoint)
		*intersectionPoint = hitPoint;
	if (intersectionNormal)
		*intersectionNormal = (hitPoint - origin).normalized();
	if (d)
		*d = t1;
	if (d2)
		*d2 = t2;

	return numIntersections;
}

int CSphere::intersects(const CLine &line, CVec3D *intersectionPoint, CVec3D *intersectionNormal, float *d, float *d2) const
{
	float t1, t2;
	int numIntersections = intersectLine(line.pos, line.dir, origin, radius, t1, t2);
	if (numIntersections == 0)
		return 0;

	CVec3D hitPoint = line.pos + t1 * line.dir;
	if (intersectionPoint)
		*intersectionPoint = hitPoint;
	if (intersectionNormal)
		*intersectionNormal = (hitPoint - origin).normalized();
	if (d)
		*d = t1;
	if (d2)
		*d2 = t2;

	return numIntersections;
}

int CSphere::intersects(const CLineSegment &l, CVec3D *intersectionPoint, CVec3D *intersectionNormal, float *d, float *d2) const
{
	float t1, t2;
	int numIntersections = intersectLine(l.begin, l.getDir(), origin, radius, t1, t2);

	if (numIntersections == 0)
		return 0;

	float lineLength = l.getLenght();
	if (t2 < 0.f || t1 > lineLength)
		return 0;
	CVec3D hitPoint = l.getPoint(t1 / lineLength);
	if (intersectionPoint)
		*intersectionPoint = hitPoint;
	if (intersectionNormal)
		*intersectionNormal = (hitPoint - origin).normalized();
	if (d)
		*d = t1 / lineLength;
	if (d2)
		*d2 = t2 / lineLength;

	return true;
}

bool CSphere::intersects(const CPlane &plane) const
{
	return plane.intersects(*this);
}

bool CSphere::intersects(const CAABBox &aabb, CVec3D *closestPointOnAABB) const
{
	return aabb.intersects(*this, closestPointOnAABB);
}

bool CSphere::intersects(const COBBox &obb, CVec3D *closestPointOnOBB) const
{
	return obb.intersects(*this, closestPointOnOBB);
}

bool CSphere::intersects(const CTriangle &triangle, CVec3D *closestPointOnTriangle) const
{
	return triangle.intersects(*this, closestPointOnTriangle);
}

bool CSphere::intersects(const CPolygon &polygon) const
{
	return polygon.intersects(*this);
}



bool CSphere::intersects(const CPolyhedron &polyhedron) const
{
	return polyhedron.intersects(*this);
}

CAABBox CSphere::minimalEnclosingAABB() const
{
	CAABBox aabb;
	aabb.setFrom(*this);
	return aabb;
}

CAABBox CSphere::maximalContainedAABB() const
{
	CAABBox aabb;
	static const float recipSqrt3 = CMath::rSqrt(3);
	float halfSideLength = radius * recipSqrt3;
	aabb.setFromCenterAndSize(origin, CVec3D(halfSideLength,halfSideLength,halfSideLength));
	return aabb;
}

void CSphere::toNegativeInfinity()
{
	origin = CVec3D(0,0,0);
	radius = CMath::NEG_INFINITY_FLOAT;
}
}  //end GEO
}  //end SMF