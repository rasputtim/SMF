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
#include "geometry/SMF_Circle.h"
#include "math/all.h"
#include "geometry/all.h"

namespace SMF {
using namespace MATH;
using namespace Util;
namespace GEO{



CRay::CRay(const CVec3D &pos_, const CVec3D &dir_)
:pos(pos_), dir(dir_)
{
	SMF_ASSERT(dir.isNormalized());
}

CRay::CRay(const CLine &line)
:pos(line.pos), dir(line.dir)
{
	SMF_ASSERT(dir.isNormalized());
}

CRay::CRay(const CLineSegment &lineSegment)
:pos(lineSegment.begin), dir(lineSegment.getDir())
{
}

bool CRay::isFinite() const
{
	return pos.isFinite() && dir.isFinite();
}

CVec3D CRay::getPoint(float d) const
{
	SMF_ASSERT(dir.isNormalized());
	return pos + d * dir;
}

void CRay::translate(const CVec3D &offset)
{
	pos += offset;
}

void CRay::transform(const CMat3D &transform)
{
	pos = transform.transform(pos);
	dir = transform.transform(dir);
}

void CRay::transform(const CMatJoint3x4 &transform)
{
	pos = transform.transformPos(pos);
	dir = transform.transformDir(dir);
}

void CRay::transform(const CMat4D &transform)
{
	pos = transform.transformPos(pos);
	dir = transform.transformDir(dir);
}

void CRay::transform(const CQuaternion &transform)
{
	pos = transform.transform(pos);
	dir = transform.transform(dir);
}

bool CRay::contains(const CVec3D &point, float distanceThreshold) const
{
	return closestPoint(point).distanceSq(point) <= distanceThreshold;
}

bool CRay::contains(const CLineSegment &lineSegment, float distanceThreshold) const
{
	return contains(lineSegment.begin, distanceThreshold) && contains(lineSegment.end, distanceThreshold);
}

bool CRay::compare(const CRay &rhs, float epsilon) const
{
	return pos.compare(rhs.pos, epsilon) && dir.compare(rhs.dir, epsilon);
}

float CRay::distance(const CVec3D &point, float *d) const
{
	return closestPoint(point, d).distance(point);
}

float CRay::distance(const CVec3D &point) const
{
	return distance(point, 0);
}

float CRay::distance(const CRay &other, float *d, float *d2) const
{
	float u2;
	CVec3D c = closestPoint(other, d, &u2);
	if (d2) *d2 = u2;
	return c.distance(other.getPoint(u2));
}

float CRay::distance(const CRay &ray) const
{
	return distance(ray, 0, 0);
}

float CRay::distance(const CLine &other, float *d, float *d2) const
{
	float u2;
	CVec3D c = closestPoint(other, d, &u2);
	if (d2) *d2 = u2;
	return c.distance(other.getPoint(u2));
}

float CRay::distance(const CLine &line) const
{
	return distance(line, 0, 0);
}

float CRay::distance(const CLineSegment &other, float *d, float *d2) const
{
	float u2;
	CVec3D c = closestPoint(other, d, &u2);
	if (d2) *d2 = u2;
	return c.distance(other.getPoint(u2));
}

float CRay::distance(const CLineSegment &lineSegment) const
{
	return distance(lineSegment, 0, 0);
}

float CRay::distance(const CSphere &sphere) const
{
	return MAX(0.f, distance(sphere.getOrigin()) - sphere.getRadius());
}
#if 0
float CRay::distance(const Capsule &capsule) const
{
	return MAX(0.f, distance(capsule.l) - capsule.r);
}
#endif
CVec3D CRay::closestPoint(const CVec3D &targetPoint, float *d) const
{
	float u = MAX(0.f, ((targetPoint - pos)* dir));
	if (d)
		*d = u;
	return getPoint(u);
}

CVec3D CRay::closestPoint(const CRay &other, float *d, float *d2) const
{
	float u, u2;
	CVec3D closestP = CLine::closestPointLineLine(pos, pos + dir, other.pos, other.pos + other.dir, &u, &u2);
	if (u < 0.f && u2 < 0.f)
	{
		closestP = closestPoint(other.pos, &u);

		CVec3D closestPoint2 = other.closestPoint(pos, &u2);
		if (closestP.distanceSq(other.pos) <= closestPoint2.distanceSq(pos))
		{
			if (d)
				*d = u;
			if (d2)
				*d2 = 0.f;
			return closestP;
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
		CVec3D pt = closestPoint(other.pos, &u);
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
		return closestP;
	}
}

CVec3D CRay::closestPoint(const CLine &other, float *d, float *d2) const
{
	float t;
	CVec3D closestPoint = CLine::closestPointLineLine(pos, pos + dir, other.pos, other.pos + other.dir, &t, d2);
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

CVec3D CRay::closestPoint(const CLineSegment &other, float *d, float *d2) const
{
	float u, u2;
	CVec3D closestP = CLine::closestPointLineLine(pos, pos + dir, other.begin, other.end, &u, &u2);
	if (u < 0.f)
	{
		if (u2 >= 0.f && u2 <= 1.f)
		{
			if (d)
				*d = 0.f;
			if (d2)
				other.closestPoint(pos, d2);
			return pos;
		}

		CVec3D p;
		float t2;

		if (u2 < 0.f)
		{
			p = other.begin;
			t2 = 0.f;
		}
		else // u2 > 1.f
		{
			p = other.end;
			t2 = 1.f;
		}

		closestP = closestPoint(p, &u);
		CVec3D closestPoint2 = other.closestPoint(pos, &u2);
		if (closestP.distanceSq(p) <= closestPoint2.distanceSq(pos))
		{
			if (d)
				*d = u;
			if (d2)
				*d2 = t2;
			return closestP;
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
	else if (u2 < 0.f)
	{
		closestP = closestPoint(other.begin, d);
		if (d2)
			*d2 = 0.f;
		return closestP;
	}
	else if (u2 > 1.f)
	{
		closestP = closestPoint(other.end, d);
		if (d2)
			*d2 = 1.f;
		return closestP;
	}
	else
	{
		if (d)
			*d = u;
		if (d2)
			*d2 = u2;
		return closestP;
	}
}

bool CRay::intersects(const CTriangle &triangle, float *d, CVec3D *intersectionPoint) const
{
	return triangle.intersects(*this, d, intersectionPoint);
}

bool CRay::intersects(const CTriangle &triangle) const
{
	float u, v;
	float t = CTriangle::intersectLineTri(pos, dir, triangle.a, triangle.b, triangle.c, u, v);
	if (t < 0.f || t == CMath::INFINITY_FLOAT)
		return false;
	return true;
}

bool CRay::intersects(const CPlane &plane, float *d) const
{
	return plane.intersects(*this, d);
}

bool CRay::intersects(const CPlane &plane) const
{
	return plane.intersects(*this, 0);
}

bool CRay::intersects(const CSphere &sphere, CVec3D *intersectionPoint, CVec3D *intersectionNormal, float *d) const
{
	return sphere.intersects(*this, intersectionPoint, intersectionNormal, d) > 0;
}

bool CRay::intersects(const CSphere &sphere) const
{
	return sphere.intersects(*this, 0, 0, 0) > 0;
}

bool CRay::intersects(const CAABBox &aabb) const
{
	return aabb.intersects(*this);
}

bool CRay::intersects(const CAABBox &aabb, float &dNear, float &dFar) const
{
	return aabb.intersects(*this, dNear, dFar);
}

bool CRay::intersects(const COBBox &obb, float &dNear, float &dFar) const
{
	return obb.intersects(*this, dNear, dFar);
}

bool CRay::intersects(const COBBox &obb) const
{
	return obb.intersects(*this);
}
#if 0
bool CRay::intersects(const Capsule &capsule) const
{
	return capsule.intersects(*this);
}
#endif
bool CRay::intersects(const CPolygon &polygon) const
{
	return polygon.intersects(*this);
}
#if 0
bool CRay::intersects(const Frustum &frustum) const
{
	return frustum.intersects(*this);
}
#endif
bool CRay::intersects(const CPolyhedron &polyhedron) const
{
	return polyhedron.intersects(*this);
}

bool CRay::intersectsDisc(const CCircle &disc) const
{
	return disc.intersectsDisc(*this);
}

CLine CRay::toLine() const
{
	return CLine(pos, dir);
}

CLineSegment CRay::toLineSegment(float d) const
{
	return CLineSegment(pos, getPoint(d));
}

CLineSegment CRay::toLineSegment(float dStart, float dEnd) const
{
	return CLineSegment(getPoint(dStart), getPoint(dEnd));
}

void CRay::projectToAxis(const CVec3D &direction, float &outMin, float &outMax) const
{
	outMin = outMax = (direction* pos);
	float d = (direction* dir);

	// Most of the time, the projection interval will be a half-infinite range, extending to either -inf or +inf.
	if (d > 1e-4f)
		outMax = CMath::INFINITY_FLOAT;
	else if (d < -1e4f)
		outMin = -CMath::INFINITY_FLOAT;
}


std::string CRay::toString() const
{
	char str[256];
	std::sprintf(str, "CRay(Pos:(%.2f, %.2f, %.2f) getDir:(%.2f, %.2f, %.2f))", pos.x, pos.y, pos.z, dir.x, dir.y, dir.z);
	return str;
}

std::ostream &operator <<(std::ostream &o, const CRay &ray)
{
	o << ray.toString();
	return o;
}



CRay operator *(const CMat3D &transform, const CRay &ray)
{
	SMF_ASSERT(transform.isInvertible());
	return CRay(transform * ray.pos, (transform * ray.dir).normalized());
}

CRay operator *(const CMatJoint3x4 &transform, const CRay &ray)
{
	SMF_ASSERT(transform.isInvertible());
	return CRay(transform.MulPos(ray.pos), transform.MulDir(ray.dir).normalized());
}

CRay operator *(const CMat4D &transform, const CRay &ray)
{
	SMF_ASSERT(transform.isInvertible());
	return CRay(transform.MulPos(ray.pos), transform.MulDir(ray.dir).normalized());
}

CRay operator *(const CQuaternion &transform, const CRay &ray)
{
	return CRay(transform * ray.pos, transform * ray.dir);
}


} //end GEO
}  //end SMF


