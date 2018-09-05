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
#include "geometry/SMF_AABBox.h"
#include "geometry/all.h"
namespace SMF {

namespace GEO {
using namespace MATH;

CAABBox bounds_zero( CVec3D::origin, CVec3D::origin );
CAABBox bounds_zeroOneCube( CVec3D( 0.0f), CVec3D( 1.0f ) );
CAABBox bounds_unitCube( CVec3D( -1.0f ), CVec3D( 1.0f ) );

/*
============
CAABBox::getRadius
============
*/
float CAABBox::getRadius() const {
	int		i;
	float	total, b0, b1;

	total = 0.0f;
	for ( i = 0; i < 3; i++ ) {
		b0 = (float)CMath::fabs( b[0][i] );
		b1 = (float)CMath::fabs( b[1][i] );
		if ( b0 > b1 ) {
			total += b0 * b0;
		} else {
			total += b1 * b1;
		}
	}
	return CMath::sqrt( total );
}

/*
============
CAABBox::getRadius
============
*/
float CAABBox::getRadius( const CVec3D &center ) const {
	int		i;
	float	total, b0, b1;

	total = 0.0f;
	for ( i = 0; i < 3; i++ ) {
		b0 = (float)CMath::fabs( center[i] - b[0][i] );
		b1 = (float)CMath::fabs( b[1][i] - center[i] );
		if ( b0 > b1 ) {
			total += b0 * b0;
		} else {
			total += b1 * b1;
		}
	}
	return CMath::sqrt( total );
}

/*
================
CAABBox::planeDistance
================
*/
float CAABBox::planeDistance( const CPlane &plane ) const {
	CVec3D center;
	float d1, d2;

	center = ( b[0] + b[1] ) * 0.5f;

	d1 = plane.getDistance( center );
	d2 = CMath::fabs( ( b[1][0] - center[0] ) * plane.getNormal()[0] ) +
			CMath::fabs( ( b[1][1] - center[1] ) * plane.getNormal()[1] ) +
				CMath::fabs( ( b[1][2] - center[2] ) * plane.getNormal()[2] );

	if ( d1 - d2 > 0.0f ) {
		return d1 - d2;
	}
	if ( d1 + d2 < 0.0f ) {
		return d1 + d2;
	}
	return 0.0f;
}

/*
================
CAABBox::planeSide
================
*/
int CAABBox::planeSide( const CPlane &plane, const float epsilon ) const {
	CVec3D center;
	float d1, d2;

	center = ( b[0] + b[1] ) * 0.5f;

	d1 = plane.getDistance( center );
	d2 = CMath::fabs( ( b[1][0] - center[0] ) * plane.getNormal()[0] ) +
			CMath::fabs( ( b[1][1] - center[1] ) * plane.getNormal()[1] ) +
				CMath::fabs( ( b[1][2] - center[2] ) * plane.getNormal()[2] );

	if ( d1 - d2 > epsilon ) {
		return PLANESIDE_FRONT;
	}
	if ( d1 + d2 < -epsilon ) {
		return PLANESIDE_BACK;
	}
	return PLANESIDE_CROSS;
}

/*
============
CAABBox::lineIntersection

  Returns true if the line intersects the bounds between the start and end point.
============
*/
bool CAABBox::lineIntersection( const CVec3D &start, const CVec3D &end ) const {
	const CVec3D center = ( b[0] + b[1] ) * 0.5f;
	const CVec3D extents = b[1] - center;
	const CVec3D lineDir = 0.5f * ( end - start );
	const CVec3D lineCenter = start + lineDir;
	const CVec3D dir = lineCenter - center;

	const float ld0 = CMath::fabs( lineDir[0] );
	if ( CMath::fabs( dir[0] ) > extents[0] + ld0 ) {
		return false;
	}

	const float ld1 = CMath::fabs( lineDir[1] );
	if ( CMath::fabs( dir[1] ) > extents[1] + ld1 ) {
		return false;
	}

	const float ld2 = CMath::fabs( lineDir[2] );
	if ( CMath::fabs( dir[2] ) > extents[2] + ld2 ) {
		return false;
	}

	const CVec3D cross = lineDir.cross( dir );

	if ( CMath::fabs( cross[0] ) > extents[1] * ld2 + extents[2] * ld1 ) {
		return false;
	}

	if ( CMath::fabs( cross[1] ) > extents[0] * ld2 + extents[2] * ld0 ) {
		return false;
	}

	if ( CMath::fabs( cross[2] ) > extents[0] * ld1 + extents[1] * ld0 ) {
		return false;
	}

	return true;
}

/*
============
CAABBox::rayIntersection

  Returns true if the ray intersects the bounds.
  The ray can intersect the bounds in both directions from the start point.
  If start is inside the bounds it is considered an intersection with scale = 0
============
*/
bool CAABBox::rayIntersection( const CVec3D &start, const CVec3D &dir, float &scale ) const {
	int i, ax0, ax1, ax2, side, inside;
	float f;
	CVec3D hit;

	ax0 = -1;
	inside = 0;
	for ( i = 0; i < 3; i++ ) {
		if ( start[i] < b[0][i] ) {
			side = 0;
		}
		else if ( start[i] > b[1][i] ) {
			side = 1;
		}
		else {
			inside++;
			continue;
		}
		if ( dir[i] == 0.0f ) {
			continue;
		}
		f = ( start[i] - b[side][i] );
		if ( ax0 < 0 || CMath::fabs( f ) > CMath::fabs( scale * dir[i] ) ) {
			scale = - ( f / dir[i] );
			ax0 = i;
		}
	}

	if ( ax0 < 0 ) {
		scale = 0.0f;
		// return true if the start point is inside the bounds
		return ( inside == 3 );
	}

	ax1 = (ax0+1)%3;
	ax2 = (ax0+2)%3;
	hit[ax1] = start[ax1] + scale * dir[ax1];
	hit[ax2] = start[ax2] + scale * dir[ax2];

	return ( hit[ax1] >= b[0][ax1] && hit[ax1] <= b[1][ax1] &&
				hit[ax2] >= b[0][ax2] && hit[ax2] <= b[1][ax2] );
}

/*
============
CAABBox::fromTransformedBounds
============
*/
void CAABBox::fromTransformedBounds( const CAABBox &bounds, const CVec3D &origin, const CMat3D &axis ) {
	int i;
	CVec3D center, extents, rotatedExtents;

	center = (bounds[0] + bounds[1]) * 0.5f;
	extents = bounds[1] - center;

	for ( i = 0; i < 3; i++ ) {
		// Remember that CMat 3d is Column-major.that is why it is trnasposed here
		rotatedExtents[i] = CMath::fabs( extents[0] * axis[i][0] ) +
							CMath::fabs( extents[1] * axis[i][1] ) +
							CMath::fabs( extents[2] * axis[i][2] );
	}

	center = origin + center * axis;
	b[0] = center - rotatedExtents;
	b[1] = center + rotatedExtents;
}

/*
============
CAABBox::fromPoints

  Most tight bounds for a point set.
============
*/
void CAABBox::fromPoints( const CVec3D *points, const int numPoints ) {
	SIMDProcessor->minMax( b[0], b[1], points, numPoints );
}

/*
============
CAABBox::fromPointTranslation

  Most tight bounds for the translational movement of the given point.
============
*/
void CAABBox::fromPointTranslation( const CVec3D &point, const CVec3D &translation ) {
	int i;

	for ( i = 0; i < 3; i++ ) {
		if ( translation[i] < 0.0f ) {
			b[0][i] = point[i] + translation[i];
			b[1][i] = point[i];
		}
		else {
			b[0][i] = point[i];
			b[1][i] = point[i] + translation[i];
		}
	}
}

/*
============
CAABBox::FromBoundsTranslation

  Most tight bounds for the translational movement of the given bounds.
============
*/
void CAABBox::FromBoundsTranslation( const CAABBox &bounds, const CVec3D &origin, const CMat3D &axis, const CVec3D &translation ) {
	int i;

	if ( axis.isRotated() ) {
		fromTransformedBounds( bounds, origin, axis );
	}
	else {
		b[0] = bounds[0] + origin;
		b[1] = bounds[1] + origin;
	}
	for ( i = 0; i < 3; i++ ) {
		if ( translation[i] < 0.0f ) {
			b[0][i] += translation[i];
		}
		else {
			b[1][i] += translation[i];
		}
	}
}
#if 0
/*
================
boundsForPointRotation

  only for rotations < 180 degrees
================
*/
CAABBox boundsForPointRotation( const CVec3D &start, const CRotation &rotation ) {
	int i;
	float radiusSqr;
	CVec3D v1, v2;
	CVec3D origin, axis, end;
	CAABBox bounds;

	end = start * rotation;
	axis = rotation.getVec();
	origin = rotation.getOrigin() + axis * ( axis * ( start - rotation.getOrigin() ) );
	radiusSqr = ( start - origin ).getLengthSqr();
	v1 = ( start - origin ).cross( axis );
	v2 = ( end - origin ).cross( axis );

	for ( i = 0; i < 3; i++ ) {
		// if the derivative changes sign along this axis during the rotation from start to end
		if ( ( v1[i] > 0.0f && v2[i] < 0.0f ) || ( v1[i] < 0.0f && v2[i] > 0.0f ) ) {
			if ( ( 0.5f * (start[i] + end[i]) - origin[i] ) > 0.0f ) {
				bounds[0][i] = min( start[i], end[i] );
				bounds[1][i] = origin[i] + CMath::sqrt( radiusSqr * ( 1.0f - axis[i] * axis[i] ) );
			}
			else {
				bounds[0][i] = origin[i] - CMath::sqrt( radiusSqr * ( 1.0f - axis[i] * axis[i] ) );
				bounds[1][i] = MAX( start[i], end[i] );
			}
		}
		else if ( start[i] > end[i] ) {
			bounds[0][i] = end[i];
			bounds[1][i] = start[i];
		}
		else {
			bounds[0][i] = start[i];
			bounds[1][i] = end[i];
		}
	}

	return bounds;
}

/*
============
CAABBox::fromPointRotation

  Most tight bounds for the rotational movement of the given point.
============
*/
void CAABBox::fromPointRotation( const CVec3D &point, const CRotation &rotation ) {
	float radius;

	if ( CMath::fabs( rotation.getAngle() ) < 180.0f ) {
		(*this) = boundsForPointRotation( point, rotation );
	}
	else {

		radius = ( point - rotation.getOrigin() ).getLenght();

		// FIXME: these bounds are usually way larger
		b[0].set( -radius, -radius, -radius );
		b[1].set( radius, radius, radius );
	}
}

/*
============
CAABBox::fromBoundsRotation

  Most tight bounds for the rotational movement of the given bounds.
============
*/
void CAABBox::fromBoundsRotation( const CAABBox &bounds, const CVec3D &origin, const CMat3D &axis, const CRotation &rotation ) {
	int i;
	float radius;
	CVec3D point;
	CAABBox rBounds;

	if ( CMath::fabs( rotation.getAngle() ) < 180.0f ) {

		(*this) = boundsForPointRotation( bounds[0] * axis + origin, rotation );
		for ( i = 1; i < 8; i++ ) {
			point[0] = bounds[(i^(i>>1))&1][0];
			point[1] = bounds[(i>>1)&1][1];
			point[2] = bounds[(i>>2)&1][2];
			(*this) += boundsForPointRotation( point * axis + origin, rotation );
		}
	}
	else {

		point = (bounds[1] - bounds[0]) * 0.5f;
		radius = (bounds[1] - point).getLenght() + (point - rotation.getOrigin()).getLenght();

		// FIXME: these bounds are usually way larger
		b[0].set( -radius, -radius, -radius );
		b[1].set( radius, radius, radius );
	}
}
#endif
/*
============
CAABBox::toPoints
============
*/
void CAABBox::toPoints( CVec3D points[8] ) const {
	for ( int i = 0; i < 8; i++ ) {
		points[i][0] = b[(i^(i>>1))&1][0];
		points[i][1] = b[(i>>1)&1][1];
		points[i][2] = b[(i>>2)&1][2];
	}
}

CVec3D CAABBox::closestPoint(const CVec3D &targetPoint) const
{
	return targetPoint.clamp(b[0], b[1]); 
}

float CAABBox::distance(const CVec3D &point) const
{
	///\todo This function could be slightly optimized. See Christer Ericson's
	/// Real-Time Collision Detection, p.131.
	return closestPoint(point).distance(point);
}

float CAABBox::distance(const CSphere &sphere) const
{
	return MAX(0.0f, distance(sphere.getOrigin()) - sphere.getRadius());
}

bool CAABBox::intersectLineAABB(const CVec3D &linePos, const CVec3D &lineDir, float &tNear, float &tFar) const
{
	// Never call the SSE version here. The SSE version does not output tNear and tFar, because
	// the memory stores have been profiled to make it slower than the CPP version. Therefore the SSE
	// version does not output tNear and tFar (profile shows it to be about 10x faster than the CPP version).
	return CSIMD::getGenProcessor()->intersectLineAABB(*this,CVec4D(linePos, 1.f), CVec4D(lineDir, 0.f), tNear, tFar);
}

bool CAABBox::intersects(const CLine &line) const
{
	float tNear = -CMath::INFINITY_FLOAT;
	float tFar = CMath::INFINITY_FLOAT;

if(CMath::MATH_AUTOMATIC_SSE){
	return CSIMD::getProcessor()->intersectLineAABB(*this,CVec4D(line.pos, 1.f), CVec4D(line.dir, 0.f), tNear, tFar);
}else{
	return CSIMD::getGenProcessor()->intersectLineAABB(*this,CVec4D(line.pos, 1.f), CVec4D(line.dir, 0.f), tNear, tFar);
}
}

bool CAABBox::intersects(const CRay &ray) const
{
	float tNear = 0;
	float tFar = CMath::INFINITY_FLOAT;

if(CMath::MATH_AUTOMATIC_SSE){
	return CSIMD::getProcessor()->intersectLineAABB(*this,CVec4D(ray.pos, 1.f), CVec4D(ray.dir, 0.f), tNear, tFar);
}else{
	return CSIMD::getGenProcessor()->intersectLineAABB(*this,CVec4D(ray.pos, 1.f), CVec4D(ray.dir, 0.f), tNear, tFar);
}
}

bool CAABBox::intersects(const CLineSegment &lineSegment) const
{
	CVec3D dir = lineSegment.end - lineSegment.begin;
	float len = dir.getLenght();
	if (len <= 1e-4f) // Degenerate line segment? Fall back to point-in-CAABBox test.
		return contains(lineSegment.begin);

	float invLen = 1.f / len;
	dir *= invLen;
	float tNear = 0.f, tFar = len;
if(CMath::MATH_AUTOMATIC_SSE){
	return CSIMD::getProcessor()->intersectLineAABB(*this,CVec4D(lineSegment.begin, 1.f), CVec4D(dir, 0.f), tNear, tFar);
}else{
	return CSIMD::getGenProcessor()->intersectLineAABB(*this,CVec4D(lineSegment.begin, 1.f),CVec4D(dir, 0.f), tNear, tFar);
}
}


bool CAABBox::intersects(const CAABBox &aabb) const
{
	// If any of the cardinal X,Y,Z axes is a separating axis, then
	// there is no intersection.
	return b[0].x < aabb.b[1].x &&
	       b[0].y < aabb.b[1].y &&
	       b[0].z < aabb.b[1].z &&
	       aabb.b[0].x < b[1].x &&
	       aabb.b[0].y < b[1].y &&
	       aabb.b[0].z < b[1].z;
}


bool CAABBox::intersects(const CRay &ray, float &dNear, float &dFar) const
{
	dNear = 0.f;
	dFar = CMath::INFINITY_FLOAT;
	return intersectLineAABB(ray.pos, ray.dir, dNear, dFar);
}

bool CAABBox::intersects(const CLine &line, float &dNear, float &dFar) const
{
	dNear = -CMath::INFINITY_FLOAT;
	dFar = CMath::INFINITY_FLOAT;
	return intersectLineAABB(line.pos, line.dir, dNear, dFar);
}

bool CAABBox::intersects(const CLineSegment &lineSegment, float &dNear, float &dFar) const
{
	CVec3D dir = lineSegment.end - lineSegment.begin;
	float len = dir.getLenght();
	if (len <= 1e-4f) // Degenerate line segment? Fall back to point-in-CAABBox test.
	{
		dNear = 0.f;
		dFar = 1.f;
		return contains(lineSegment.begin);
	}
	float invLen = 1.f / len;
	dir *= invLen;
	dNear = 0.f;
	dFar = len;
	bool hit = intersectLineAABB(lineSegment.begin, dir, dNear, dFar);
	dNear *= invLen;
	dFar *= invLen;
	return hit;
}

bool CAABBox::intersects(const CPlane &plane) const
{
	return plane.intersects(*this);
}


bool CAABBox::intersects(const COBBox &obb) const
{
	return obb.intersects(*this);
}

bool CAABBox::intersects(const CSphere &sphere, CVec3D *closestPointOnAABB) const
{
	// find the point on this CAABBox closest to the sphere center.
	CVec3D pt = closestPoint(sphere.getOrigin());

	// If that point is inside sphere, the CAABBox and sphere intersect.
	if (closestPointOnAABB)
		*closestPointOnAABB = pt;

	return pt.distanceSq(sphere.getOrigin()) <= sphere.getRadius() * sphere.getRadius();
}


bool CAABBox::intersects(const CTriangle &triangle) const
{
	return triangle.intersects(*this);
}

bool CAABBox::intersects(const CPolygon &polygon) const
{
	return toPolyhedron().intersects(polygon);
}


bool CAABBox::intersects(const CPolyhedron &polyhedron) const
{
	return polyhedron.intersects(*this);
}

CPolyhedron CAABBox::toPolyhedron() const
{
	// Note to maintainer: This function is an exact copy of COBBox:toPolyhedron() and Frustum::toPolyhedron().
	CPolyhedron p;
	// Populate the corners of this CAABBox.
	// The will be in the order 0: ---, 1: --+, 2: -+-, 3: -++, 4: +--, 5: +-+, 6: ++-, 7: +++.
	for(int i = 0; i < 8; ++i)
		p.v.push_back(cornerPoint(i));

	// Generate the 6 faces of this CAABBox.
	const int faces[6][4] =
	{
		{ 0, 1, 3, 2 }, // X-
		{ 4, 6, 7, 5 }, // X+
		{ 0, 4, 5, 1 }, // Y-
		{ 7, 6, 2, 3 }, // Y+
		{ 0, 2, 6, 4 }, // Z-
		{ 1, 5, 7, 3 }, // Z+
	};

	for(int f = 0; f < 6; ++f)
	{
		CPolyhedron::Face face;
		for(int v = 0; v < 4; ++v)
			face.v.push_back(faces[f][v]);
		p.f.push_back(face);
	}

	return p;
}


bool CAABBox::contains(const CVec3D &point) const
{
	return b[0].x <= point.x && point.x <= b[1].x &&
		   b[0].y <= point.y && point.y <= b[1].y &&
		   b[0].z <= point.z && point.z <= b[1].z;
}

bool CAABBox::contains(const CLineSegment &lineSegment) const
{
	return contains(lineSegment.begin) && contains(lineSegment.end);
}

bool CAABBox::contains(const CAABBox &aabb) const
{
	return contains(aabb.b[0]) && contains(aabb.b[1]);
}

bool CAABBox::contains(const COBBox &obb) const
{
	return contains(obb.minimalEnclosingAABB());
}

bool CAABBox::contains(const CSphere &sphere) const
{
	///\todo Optimize.
	return contains(sphere.minimalEnclosingAABB());
}

bool CAABBox::contains(const CTriangle &triangle) const
{
	return contains(triangle.a) && contains(triangle.b) && contains(triangle.c);
}

bool CAABBox::contains(const CPolygon &polygon) const
{
	return contains(polygon.minimalEnclosingAABB());
}


bool CAABBox::contains(const CPolyhedron &polyhedron) const
{
	return contains(polyhedron.minimalEnclosingAABB());
}


CLineSegment CAABBox::edge(int edgeIndex) const
{
	SMF_ASSERT(0 <= edgeIndex && edgeIndex <= 11);
	switch(edgeIndex)
	{
		default: // For release builds where SMF_ASSERT() is disabled, return always the first option if out-of-bounds.
		/* For documentation, here's the segments that are returned:
		case 0: return CLineSegment(cornerPoint(0), cornerPoint(1));
		case 1: return CLineSegment(cornerPoint(0), cornerPoint(2));
		case 2: return CLineSegment(cornerPoint(0), cornerPoint(4));
		case 3: return CLineSegment(cornerPoint(1), cornerPoint(3));
		case 4: return CLineSegment(cornerPoint(1), cornerPoint(5));
		case 5: return CLineSegment(cornerPoint(2), cornerPoint(3));
		case 6: return CLineSegment(cornerPoint(2), cornerPoint(6));
		case 7: return CLineSegment(cornerPoint(3), cornerPoint(7));
		case 8: return CLineSegment(cornerPoint(4), cornerPoint(5));
		case 9: return CLineSegment(cornerPoint(4), cornerPoint(6));
		case 10: return CLineSegment(cornerPoint(5), cornerPoint(7));
		case 11: return CLineSegment(cornerPoint(6), cornerPoint(7));
		*/
		// Force-optimize to avoid calling to cornerPoint for another switch-case statement.
		case 0: return CLineSegment(CVec3D(b[0].x, b[0].y, b[0].z), CVec3D(b[0].x, b[0].y, b[1].z));
		case 1: return CLineSegment(CVec3D(b[0].x, b[0].y, b[0].z), CVec3D(b[0].x, b[1].y, b[0].z));
		case 2: return CLineSegment(CVec3D(b[0].x, b[0].y, b[0].z), CVec3D(b[1].x, b[0].y, b[0].z));
		case 3: return CLineSegment(CVec3D(b[0].x, b[0].y, b[1].z), CVec3D(b[0].x, b[1].y, b[1].z));
		case 4: return CLineSegment(CVec3D(b[0].x, b[0].y, b[1].z), CVec3D(b[1].x, b[0].y, b[1].z));
		case 5: return CLineSegment(CVec3D(b[0].x, b[1].y, b[0].z), CVec3D(b[0].x, b[1].y, b[1].z));
		case 6: return CLineSegment(CVec3D(b[0].x, b[1].y, b[0].z), CVec3D(b[1].x, b[1].y, b[0].z));
		case 7: return CLineSegment(CVec3D(b[0].x, b[1].y, b[1].z), CVec3D(b[1].x, b[1].y, b[1].z));
		case 8: return CLineSegment(CVec3D(b[1].x, b[0].y, b[0].z), CVec3D(b[1].x, b[0].y, b[1].z));
		case 9: return CLineSegment(CVec3D(b[1].x, b[0].y, b[0].z), CVec3D(b[1].x, b[1].y, b[0].z));
		case 10: return CLineSegment(CVec3D(b[1].x, b[0].y, b[1].z), CVec3D(b[1].x, b[1].y, b[1].z));
		case 11: return CLineSegment(CVec3D(b[1].x, b[1].y, b[0].z), CVec3D(b[1].x, b[1].y, b[1].z));
	}
}

CVec3D CAABBox::cornerPoint(int cornerIndex) const
{
	SMF_ASSERT(0 <= cornerIndex && cornerIndex <= 7);
	switch(cornerIndex)
	{
		default: // For release builds where SMF_ASSERT() is disabled, return always the first option if out-of-bounds.
		case 0: return CVec3D(b[0].x, b[0].y, b[0].z);
		case 1: return CVec3D(b[0].x, b[0].y, b[1].z);
		case 2: return CVec3D(b[0].x, b[1].y, b[0].z);
		case 3: return CVec3D(b[0].x, b[1].y, b[1].z);
		case 4: return CVec3D(b[1].x, b[0].y, b[0].z);
		case 5: return CVec3D(b[1].x, b[0].y, b[1].z);
		case 6: return CVec3D(b[1].x, b[1].y, b[0].z);
		case 7: return CVec3D(b[1].x, b[1].y, b[1].z);
	}
}

CVec3D CAABBox::extremePoint(const CVec3D &direction) const
{
	CVec3D pt;
	pt.x = (direction.x >= 0.f ? b[1].x : b[0].x);
	pt.y = (direction.y >= 0.f ? b[1].y : b[0].y);
	pt.z = (direction.z >= 0.f ? b[1].z : b[0].z);
	return pt;
}

CVec3D CAABBox::pointOnEdge(int edgeIndex, float u) const
{
	SMF_ASSERT(0 <= edgeIndex && edgeIndex <= 11);
	SMF_ASSERT(0 <= u && u <= 1.f);

	CVec3D d = b[1] - b[0];
	switch(edgeIndex)
	{
	default: // For release builds where SMF_ASSERT() is disabled, return always the first option if out-of-bounds.
	case 0: return CVec3D(b[0].x, b[0].y, b[0].z + u * d.z);
	case 1: return CVec3D(b[0].x, b[1].y, b[0].z + u * d.z);
	case 2: return CVec3D(b[1].x, b[0].y, b[0].z + u * d.z);
	case 3: return CVec3D(b[1].x, b[1].y, b[0].z + u * d.z);

	case 4: return CVec3D(b[0].x, b[0].y + u * d.y, b[0].z);
	case 5: return CVec3D(b[1].x, b[0].y + u * d.y, b[0].z);
	case 6: return CVec3D(b[0].x, b[0].y + u * d.y, b[1].z);
	case 7: return CVec3D(b[1].x, b[0].y + u * d.y, b[1].z);

	case 8: return CVec3D(b[0].x + u * d.x, b[0].y, b[0].z);
	case 9: return CVec3D(b[0].x + u * d.x, b[0].y, b[1].z);
	case 10: return CVec3D(b[0].x + u * d.x, b[1].y, b[0].z);
	case 11: return CVec3D(b[0].x + u * d.x, b[1].y, b[1].z);
	}
}

CVec3D CAABBox::faceCenterPoint(int faceIndex) const
{
	SMF_ASSERT(0 <= faceIndex && faceIndex <= 5);

	CVec3D center = (b[0] + b[1]) / 2.f;
	switch(faceIndex)
	{
	default: // For release builds where SMF_ASSERT() is disabled, return always the first option if out-of-bounds.
	case 0: return CVec3D(b[0].x, center.y, center.z);
	case 1: return CVec3D(b[1].x, center.y, center.z);
	case 2: return CVec3D(center.x, b[0].y, center.z);
	case 3: return CVec3D(center.x, b[1].y, center.z);
	case 4: return CVec3D(center.x, center.y, b[0].z);
	case 5: return CVec3D(center.x, center.y, b[1].z);
	}
}

CVec3D CAABBox::facePoint(int faceIndex, float u, float v) const
{
	SMF_ASSERT(0 <= faceIndex && faceIndex <= 5);
	SMF_ASSERT(0 <= u && u <= 1.f);
	SMF_ASSERT(0 <= v && v <= 1.f);

	CVec3D d = b[1] - b[0];
	switch(faceIndex)
	{
	default: // For release builds where SMF_ASSERT() is disabled, return always the first option if out-of-bounds.
	case 0: return CVec3D(b[0].x, b[0].y + u * d.y, b[0].z + v * d.z);
	case 1: return CVec3D(b[1].x, b[0].y + u * d.y, b[0].z + v * d.z);
	case 2: return CVec3D(b[0].x + u * d.x, b[0].y, b[0].z + v * d.z);
	case 3: return CVec3D(b[0].x + u * d.x, b[1].y, b[0].z + v * d.z);
	case 4: return CVec3D(b[0].x + u * d.x, b[0].y + v * d.y, b[0].z);
	case 5: return CVec3D(b[0].x + u * d.x, b[0].y + v * d.y, b[1].z);
	}
}
void CAABBox::setFromCenterAndSize(const CVec3D &center, const CVec3D &size)
{
	CVec3D halfSize = 0.5f * size;
	b[0] = center - halfSize;
	b[1] = center + halfSize;
}
void CAABBox::setFrom(const COBBox &obb)
{
	CVec3D halfSize = (obb.axis.Row(0)*obb.extents[0]).abs() +(obb.axis.Row(1)*obb.extents[1]).abs() + (obb.axis.Row(2)*obb.extents[2]).abs();
	setFromCenterAndSize(obb.center, 2.0f*halfSize);
}

void CAABBox::setFrom(const CSphere &s)
{
	b[0] = s.getOrigin() - CVec3D(s.getRadius(), s.getRadius(),s.getRadius());
	b[1] = s.getOrigin()+ CVec3D(s.getRadius(), s.getRadius(), s.getRadius());
}

void CAABBox::enclose(const CVec3D &point)
{
	b[0] = b[0].Min( point);
	b[1] = b[1].Max( point);
}

void CAABBox::enclose(const CLineSegment &lineSegment)
{
	enclose(lineSegment.begin);
	enclose(lineSegment.end);
}

void CAABBox::enclose(const CAABBox &aabb)
{
	b[0] = b[0].Min(aabb.b[0]);
	b[1] = b[1].Max( aabb.b[1]);
}

void CAABBox::enclose(const COBBox &obb)
{
	for(int i = 0; i < 8; ++i)
		enclose(obb.cornerPoint(i));
}

void CAABBox::enclose(const CSphere &sphere)
{
	enclose(sphere.getOrigin() - CVec3D(sphere.getRadius(),sphere.getRadius(),sphere.getRadius()));
	enclose(sphere.getOrigin() + CVec3D(sphere.getRadius(),sphere.getRadius(),sphere.getRadius()));
}

void CAABBox::enclose(const CTriangle &triangle)
{
	enclose(triangle.a);
	enclose(triangle.b);
	enclose(triangle.c);
}

void CAABBox::enclose(const CPolygon &polygon)
{
	for(int i = 0; i < polygon.numVertices(); ++i)
		enclose(polygon.vertex(i));
}

void CAABBox::enclose(const CPolyhedron &polyhedron)
{
	for(int i = 0; i < polyhedron.numVertices(); ++i)
		enclose(polyhedron.vertex(i));
}

void CAABBox::enclose(const CVec3D *pointArray, int numPoints)
{
	SMF_ASSERT(pointArray || numPoints == 0);
	if (!pointArray)
		return;
	for(int i = 0; i < numPoints; ++i)
		enclose(pointArray[i]);
}

void CAABBox::toNegativeInfinity()
{
	b[0].setFromScalar(CMath::INFINITY_FLOAT);
	b[1].setFromScalar(CMath::NEG_INFINITY_FLOAT);
}

void CAABBox::scale(const CVec3D &centerPoint, float scaleFactor)
{
	return scale(centerPoint, CVec3D(scaleFactor, scaleFactor, scaleFactor));
}

void CAABBox::scale(const CVec3D &centerPoint, const CVec3D &scaleFactor)
{
	CMatJoint3x4 transform = CMatJoint3x4::scale(scaleFactor, centerPoint);
	b[0] = transform.MulPos(b[0]);
	b[1] = transform.MulPos(b[1]);
}

} // end GEO
} //end SMF
