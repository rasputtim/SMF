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
#include "math/all.h"
#include "geometry/all.h"

namespace SMF {
using namespace MATH;
using namespace Util;
namespace GEO{


/// A helper function to compute the line-line closest point.
/** This code is adapted from http://paulbourke.net/geometry/lineline3d/ .
	dmnop = (xm - xn)(xo - xp) + (ym - yn)(yo - yp) + (zm - zn)(zo - zp).
	\param v An array of four floats: [0]: line 0 start. [1]: line 0 end. [2]: line 1 start. [3]: line 1 end.
	\param m An index in the range [0, 3].
	\param n An index in the range [0, 3].
	\param o An index in the range [0, 3].
	\param p An index in the range [0, 3]. */
float Dmnop(const CVec3D *v, int m, int n, int o, int p)
{
	return (v[m].x - v[n].x) * (v[o].x - v[p].x) + (v[m].y - v[n].y) * (v[o].y - v[p].y) + (v[m].z - v[n].z) * (v[o].z - v[p].z);
}

/// Computes the closest point pair on two lines.
/** The first line is specified by two points start0 and end0. The second line is specified by
	two points start1 and end1.
	The implementation of this function follows http://paulbourke.net/geometry/lineline3d/ .
	\param start0 The starting point of the first line.
	\param end0 The ending point of the first line.
	\param start1 The starting point of the second line.
	\param end1 The ending point of the second line.
	\param d [out] If specified, receives the normalized distance of the closest point along the first line.
		This pointer may be left null.
	\param d2 [out] If specified, receives the normalized distance of the closest point along the second line.
		This pointer may be left null.
	\return Returns the closest point on line start0<->end0 to the second line.
	\note This is a low-level utility function. You probably want to use closestPoint() or distance() instead.
	\see closestPoint(), distance(). */
CVec3D CLine::closestPointLineLine(CVec3D start0, CVec3D end0, CVec3D start1, CVec3D end1, float *d, float *d2)
{
	const CVec3D v[4] = { start0, end0, start1, end1 };

	float d0232 = Dmnop(v,0,2,3,2);
	float d3210 = Dmnop(v,3,2,1,0);
	float d3232 = Dmnop(v,3,2,3,2);
	float mu = (d0232 * d3210 - Dmnop(v,0,2,1,0)*d3232) / (Dmnop(v,1,0,1,0)*Dmnop(v,3,2,3,2) - Dmnop(v,3,2,1,0)*Dmnop(v,3,2,1,0));
	if (d)
		*d = mu;

	if (d2)
		*d2 = (d0232 + mu * d3210) / d3232;

	return start0 + mu * (end0 - start0);
}

CLine::CLine(const CVec3D &pos_, const CVec3D &dir_)
:pos(pos_), dir(dir_)
{
	SMF_ASSERT(dir.isNormalized());
}

CLine::CLine(const CRay &ray)
:pos(ray.pos), dir(ray.dir)
{
	SMF_ASSERT(dir.isNormalized());
}

CLine::CLine(const CLineSegment &lineSegment)
:pos(lineSegment.begin), dir(lineSegment.getDir())
{
}

bool CLine::isFinite() const
{
	return pos.isFinite() && dir.isFinite();
}

CVec3D CLine::getPoint(float d) const
{
	SMF_ASSERT(dir.isNormalized());
	return pos + d * dir;
}

void CLine::translate(const CVec3D &offset)
{
	pos += offset;
}

void CLine::transform(const CMat3D &transform)
{
	pos = transform.transform(pos);
	dir = transform.transform(dir);
}

void CLine::transform(const CMatJoint3x4 &transform)
{
	pos = transform.transformPos(pos);
	dir = transform.transformDir(dir);
}

void CLine::transform(const CMat4D &transform)
{
	pos = transform.transformPos(pos);
	dir = transform.transformDir(dir);
}

void CLine::transform(const CQuaternion &transform)
{
	pos = transform.transform(pos);
	dir = transform.transform(dir);
}

bool CLine::contains(const CVec3D &point, float distanceThreshold) const
{
	return closestPoint(point).distanceSq(point) <= distanceThreshold;
}

bool CLine::contains(const CRay &ray, float epsilon) const
{
	return contains(ray.pos, epsilon) && dir.compare(ray.dir, epsilon);
}

bool CLine::contains(const CLineSegment &lineSegment, float epsilon) const
{
	return contains(lineSegment.begin, epsilon) && contains(lineSegment.end, epsilon);
}

bool CLine::compare(const CLine &line, float epsilon) const
{
	SMF_ASSERT(dir.isNormalized());
	SMF_ASSERT(line.dir.isNormalized());
	// If the point of the other line is on this line, and the two lines point to the same, or exactly reverse directions,
	// they must be equal.
	return contains(line.pos, epsilon) && CMath::equalsAbs(CMath::fabs(dir*(line.dir)), 1.f, epsilon);
}

float CLine::distance(const CVec3D &point, float *d) const
{
	return closestPoint(point, d).distance(point);
}

float CLine::distance(const CRay &other, float *d, float *d2) const
{
	float u2;
	CVec3D c = closestPoint(other, d, &u2);
	if (d2) *d2 = u2;
	return c.distance(other.getPoint(u2));
}

float CLine::distance(const CRay &other) const
{
	return distance(other, 0, 0);
}

float CLine::distance(const CLine &other, float *d, float *d2) const
{
	float u2;
	CVec3D c = closestPoint(other, d, &u2);
	if (d2) *d2 = u2;
	return c.distance(other.getPoint(u2));
}

float CLine::distance(const CLine &other) const
{
	return distance(other, 0, 0);
}

float CLine::distance(const CLineSegment &other, float *d, float *d2) const
{
	float u2;
	CVec3D c = closestPoint(other, d, &u2);
	if (d2) *d2 = u2;
	SMF_ASSERT(u2 >= 0.f);
	SMF_ASSERT(u2 <= 1.f);
	return c.distance(other.getPoint(u2));
}

float CLine::distance(const CLineSegment &other) const
{
	return distance(other, 0, 0);
}

float CLine::distance(const CSphere &other) const
{
	return MAX(0.f, distance(other.getOrigin()) - other.getRadius());
}
#if 0
float CLine::distance(const Capsule &other) const
{
	return MAX(0.f, distance(other.l) - other.r);
}
#endif
bool CLine::intersects(const CTriangle &triangle, float *d, CVec3D *intersectionPoint) const
{
	return triangle.intersects(*this, d, intersectionPoint);
}

bool CLine::intersects(const CPlane &plane, float *d) const
{
	return plane.intersects(*this, d);
}

bool CLine::intersects(const CSphere &s, CVec3D *intersectionPoint, CVec3D *intersectionNormal, float *d) const
{
	return s.intersects(*this, intersectionPoint, intersectionNormal, d) > 0;
}

bool CLine::intersects(const CAABBox &aabb) const
{
	return aabb.intersects(*this);
}

bool CLine::intersects(const CAABBox &aabb, float &dNear, float &dFar) const
{
	return aabb.intersects(*this, dNear, dFar);
}

bool CLine::intersects(const COBBox &obb) const
{
	return obb.intersects(*this);
}

bool CLine::intersects(const COBBox &obb, float &dNear, float &dFar) const
{
	return obb.intersects(*this, dNear, dFar);
}
#if 0
bool CLine::intersects(const Capsule &capsule) const
{
	return capsule.intersects(*this);
}
#endif
bool CLine::intersects(const CPolygon &polygon) const
{
	return polygon.intersects(*this);
}
#if 0
bool CLine::intersects(const Frustum &frustum) const
{
	return frustum.intersects(*this);
}
#endif
bool CLine::intersects(const CPolyhedron &polyhedron) const
{
	return polyhedron.intersects(*this);
}

bool CLine::intersectsDisc(const CCircle &disc) const
{
	return disc.intersectsDisc(*this);
}

CVec3D CLine::closestPoint(const CVec3D &targetPoint, float *d) const
{
	float u = ((targetPoint - pos)* dir);
	if (d)
		*d = u;
	return getPoint(u);
}

CVec3D CLine::closestPoint(const CRay &other, float *d, float *d2) const
{
	float t2;
	CVec3D closestPoint_ = closestPointLineLine(pos, pos + dir, other.pos, other.pos + other.dir, d, &t2);
	if (t2 <= 0.f)
	{
		if (d2)
			*d2 = 0.f;
		return closestPoint(other.pos, d);
	}
	else
	{
		if (d2)
			*d2 = t2;
		return closestPoint_;
	}
}

CVec3D CLine::closestPoint(const CLine &other, float *d, float *d2) const
{
	return closestPointLineLine(pos, pos + dir, other.pos, other.pos + other.dir, d, d2);
}

CVec3D CLine::closestPoint(const CLineSegment &other, float *d, float *d2) const
{
	float t2;
	CVec3D closestPoint_ = closestPointLineLine(pos, pos + dir, other.begin, other.end, d, &t2);
	if (t2 <= 0.f)
	{
		if (d2)
			*d2 = 0.f;
		return closestPoint(other.begin, d);
	}
	else if (t2 >= 1.f)
	{
		if (d2)
			*d2 = 1.f;
		return closestPoint(other.end, d);
	}
	else
	{
		if (d2)
			*d2 = t2;
		return closestPoint_;
	}
}

CVec3D CLine::closestPoint(const CTriangle &triangle, float *outU, float *outV, float *outD) const
{
	///\todo Optimize this function!
	CVec3D closestPointTriangle = triangle.closestPoint(*this);
	if (outU || outV)
	{
		CVec2D uv = triangle.barycentricUV(closestPointTriangle);
		if (outU)
			*outU = uv.x;
		if (outV)
			*outV = uv.y;
	}
	return closestPoint(closestPointTriangle, outD);
}

bool CLine::areCollinear(const CVec3D &p1, const CVec3D &p2, const CVec3D &p3, float epsilon)
{
	return CVec3D::areCollinear(p1, p2, p3, epsilon);
}

CRay CLine::toRay() const
{
	return CRay(pos, dir);
}

CLineSegment CLine::toLineSegment(float d) const
{
	return CLineSegment(pos, getPoint(d));
}

void CLine::projectToAxis(const CVec3D &direction, float &outMin, float &outMax) const
{
	// Most of the time, the projection of a line spans the whole 1D axis.
	// As a special case, if the line is perpendicular to the direction vector in question,
	// then the projection interval of this line is a single point.
	if (dir.isPerpendicular(direction))
		outMin = outMax = (direction* pos);
	else
	{
		outMin = -CMath::INFINITY_FLOAT;
		outMax = CMath::INFINITY_FLOAT;
	}
}

CLineSegment CLine::toLineSegment(float dStart, float dEnd) const
{
	return CLineSegment(getPoint(dStart), getPoint(dEnd));
}

CLine operator *(const CMat3D &transform, const CLine &l)
{
	return CLine(transform * l.pos, transform * l.dir);
}

CLine operator *(const CMatJoint3x4 &transform, const CLine &l)
{
	return CLine(transform.MulPos(l.pos), transform.MulDir(l.dir));
}

CLine operator *(const CMat4D &transform, const CLine &l)
{
	return CLine(transform.MulPos(l.pos), transform.MulDir(l.dir));
}

CLine operator *(const CQuaternion &transform, const CLine &l)
{
	return CLine(transform * l.pos, transform * l.dir);
}


std::string CLine::toString() const
{
	char str[256];
	std::sprintf(str, "CLine(pos:(%.2f, %.2f, %.2f) dir:(%.2f, %.2f, %.2f))",
		pos.x, pos.y, pos.z, dir.x, dir.y, dir.z);
	return str;
}

std::ostream &operator <<(std::ostream &o, const CLine &line)
{
	o << line.toString();
	return o;
}




} //end GEO
}  //end SMF

