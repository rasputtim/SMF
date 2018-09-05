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



#ifdef MATH_GRAPHICSENGINE_INTEROP
#include "VertexBuffer.h"
#endif



CLineSegment::CLineSegment(const CVec3D &a_, const CVec3D &b_)
:begin(a_), end(b_)
{
}

CLineSegment::CLineSegment(const CRay &ray, float d)
:begin(ray.pos), end(ray.getPoint(d))
{
}

CLineSegment::CLineSegment(const CLine &line, float d)
:begin(line.pos), end(line.getPoint(d))
{
}

CVec3D CLineSegment::getPoint(float d) const
{
	return (1.f - d) * begin + d * end;
}

CVec3D CLineSegment::centerPoint() const
{
	return (begin + end) * 0.5f;
}

void CLineSegment::reverse()
{
	Swap(begin, end);
}

CVec3D CLineSegment::getDir() const
{
	return (end - begin).normalized();
}

CVec3D CLineSegment::extremePoint(const CVec3D &direction) const
{
	return (direction*(end-begin)) >= 0.f ? end : begin;
}

void CLineSegment::translate(const CVec3D &offset)
{
	begin += offset;
	end += offset;
}

void CLineSegment::transform(const CMat3D &transform)
{
	begin = transform * begin;
	end = transform * end;
}

void CLineSegment::transform(const CMatJoint3x4 &transform)
{
	begin = transform.MulPos(begin);
	end = transform.MulPos(end);
}

void CLineSegment::transform(const CMat4D &transform)
{
	begin = transform.MulPos(begin);
	end = transform.MulPos(end);
}

void CLineSegment::transform(const CQuaternion &transform)
{
	begin = transform * begin;
	end = transform * end;
}

float CLineSegment::getLenght() const
{
	return begin.distance(end);
}

float CLineSegment::getLengthSqr() const
{
	return begin.distanceSq(end);
}

bool CLineSegment::isFinite() const
{
	return begin.isFinite() && end.isFinite();
}

bool CLineSegment::contains(const CVec3D &point, float distanceThreshold) const
{
	return closestPoint(point).distanceSq(point) <= distanceThreshold;
}

bool CLineSegment::contains(const CLineSegment &rhs, float distanceThreshold) const
{
	return contains(rhs.begin, distanceThreshold) && contains(rhs.end, distanceThreshold);
}

bool CLineSegment::compare(const CLineSegment &rhs, float e) const
{
	return (begin.compare(rhs.begin, e) && end.compare(rhs.end, e)) || (begin.compare(rhs.end, e) && end.compare(rhs.begin, e));
}

CVec3D CLineSegment::closestPoint(const CVec3D &point, float *d) const
{
	CVec3D dir = end - begin;
	float u = clamp01(((point - begin)* dir) / dir.getLengthSqr());
	if (d)
		*d = u;
	return begin + u * dir;
}

CVec3D CLineSegment::closestPoint(const CRay &other, float *d, float *d2) const
{
	float u, u2;
	other.closestPoint(*this, &u, &u2);
	if (d)
		*d = u2;
	if (d2)
		*d2 = u;
	return getPoint(u2);
}

CVec3D CLineSegment::closestPoint(const CLine &other, float *d, float *d2) const
{
	float t;
	CLine::closestPointLineLine(other.pos, other.pos + other.dir, begin, end, d2, &t);
	if (t <= 0.f)
	{
		if (d)
			*d = 0.f;
		if (d2)
			other.closestPoint(begin, d2);
		return begin;
	}
	else if (t >= 1.f)
	{
		if (d)
			*d = 1.f;
		if (d2)
			other.closestPoint(end, d2);
		return end;
	}
	else
	{
		if (d)
			*d = t;
		return getPoint(t);
	}
}

CVec3D CLineSegment::closestPoint(const CLineSegment &other, float *d, float *d2) const
{
	float u, u2;
	CVec3D closestP = CLine::closestPointLineLine(begin, end, other.begin, other.end, &u, &u2);
	if (u >= 0.f && u <= 1.f && u2 >= 0.f && u2 <= 1.f)
	{
		if (d)
			*d = u;
		if (d2)
			*d2 = u2;
		return closestP;
	}
	else if (u >= 0.f && u <= 1.f) // Only u2 is out of bounds.
	{
		CVec3D p;
		if (u2 < 0.f)
		{
			p = other.begin;
			if (d2)
				*d2 = 0.f;
		}
		else
		{
			p = other.end;
			if (d2)
				*d2 = 1.f;
		}

		return closestPoint(p, d);
	}
	else if (u2 >= 0.f && u2 <= 1.f) // Only u is out of bounds.
	{
		CVec3D p;
		if (u < 0.f)
		{
			p = begin;
			if (d)
				*d = 0.f;
		}
		else
		{
			p = end;
			if (d)
				*d = 1.f;
		}

		if (d2)
			other.closestPoint(p, d2);
		return p;
	}
	else // Both u and u2 are out of bounds.
	{
		CVec3D p;
		float t;
		if (u < 0.f)
		{
			p = begin;
			t = 0.f;
		}
		else
		{
			p = end;
			t = 1.f;
		}

		CVec3D p2;
		float t2;
		if (u2 < 0.f)
		{
			p2 = other.begin;
			t2 = 0.f;
		}
		else
		{
			p2 = other.end;
			t2 = 1.f;
		}

		float T, T2;
		closestP = closestPoint(p2, &T);
		CVec3D closestPoint2 = other.closestPoint(p, &T2);

		if (closestP.distanceSq(p2) <= closestPoint2.distanceSq(p))
		{
			if (d)
				*d = T;
			if (d2)
				*d2 = t2;
			return closestP;
		}
		else
		{
			if (d)
				*d = t;
			if (d2)
				*d2 = T2;
			return p;
		}
	}
}

float CLineSegment::distance(const CVec3D &point, float *d) const
{
	///\todo This function could be slightly optimized.
	/// See Christer Ericson's Real-Time Collision Detection, p.130.
	CVec3D closestP = closestPoint(point, d);
	return closestP.distance(point);
}

float CLineSegment::distance(const CRay &other, float *d, float *d2) const
{
	float u, u2;
	closestPoint(other, &u, &u2);
	if (d)
		*d = u;
	if (d2)
		*d2 = u2;
	return getPoint(u).distance(other.getPoint(u2));
}

float CLineSegment::distance(const CLine &other, float *d, float *d2) const
{
	float u, u2;
	CVec3D closestPoint2 = other.closestPoint(*this, &u, &u2);
	if (d)
		*d = u2;
	if (d2)
		*d2 = u;
	CVec3D closestPoint = getPoint(u2);
	return closestPoint.distance(closestPoint2);
}

float CLineSegment::distance(const CLineSegment &other, float *d, float *d2) const
{
	float u, u2;
	closestPoint(other, &u, &u2);
	if (d)
		*d = u;
	if (d2)
		*d2 = u2;
	return getPoint(u).distance(other.getPoint(u2));
}

float CLineSegment::distance(const CPlane &other) const
{
	float aDist = other.signedDistance(begin);
	float bDist = other.signedDistance(end);
	if (aDist * bDist < 0.f)
		return 0.f; // There was an intersection, so the distance is zero.
	return MIN(CMath::fabs(aDist), CMath::fabs(bDist));
}

float CLineSegment::distance(const CSphere &other) const
{
	return MAX(0.f, distance(other.getOrigin()) - other.getRadius());
}
#if 0
float CLineSegment::distance(const Capsule &other) const
{
	return MAX(0.f, distance(other.l) - other.r);
}
#endif
bool CLineSegment::intersects(const CPlane &plane) const
{
	float d = plane.signedDistance(begin);
	float d2 = plane.signedDistance(end);
	return d * d2 <= 0.f;
}
#if 0
bool CLineSegment::intersects(const Capsule &capsule) const
{
	return capsule.intersects(*this);
}
#endif
bool CLineSegment::intersects(const CPlane &plane, float *d) const
{
	return plane.intersects(*this, d);
}

bool CLineSegment::intersects(const CTriangle &triangle, float *d, CVec3D *intersectionPoint) const
{
	return triangle.intersects(*this, d, intersectionPoint);
}

bool CLineSegment::intersects(const CSphere &s, CVec3D *intersectionPoint, CVec3D *intersectionNormal, float *d) const
{
	return s.intersects(*this, intersectionPoint, intersectionNormal, d) > 0;
}

bool CLineSegment::intersects(const CAABBox &aabb) const
{
	return aabb.intersects(*this);
}

bool CLineSegment::intersects(const CAABBox &aabb, float &dNear, float &dFar) const
{
	return aabb.intersects(*this, dNear, dFar);
}

bool CLineSegment::intersects(const COBBox &obb) const
{
	return obb.intersects(*this);
}

bool CLineSegment::intersects(const COBBox &obb, float &dNear, float &dFar) const
{
	return obb.intersects(*this, dNear, dFar);
}

bool CLineSegment::intersects(const CLineSegment &lineSegment, float epsilon) const
{
	return distance(lineSegment) <= epsilon;
}

bool CLineSegment::intersects(const CPolygon &polygon) const
{
	return polygon.intersects(*this);
}
#if 0
bool CLineSegment::intersects(const Frustum &frustum) const
{
	return frustum.intersects(*this);
}
#endif
bool CLineSegment::intersects(const CPolyhedron &polyhedron) const
{
	return polyhedron.intersects(*this);
}

bool CLineSegment::intersectsDisc(const CCircle &disc) const
{
	return disc.intersectsDisc(*this);
}

CRay CLineSegment::toRay() const
{
	return CRay(begin, getDir());
}

CLine CLineSegment::toLine() const
{
	return CLine(begin, getDir());
}

void CLineSegment::projectToAxis(const CVec3D &direction, float &outMin, float &outMax) const
{
	outMin = (direction* begin);
	outMax = (direction* end);
	if (outMax < outMin)
		Swap(outMin, outMax);
}

CLineSegment operator *(const CMat3D &transform, const CLineSegment &l)
{
	return CLineSegment(transform * l.begin, transform * l.end);
}

CLineSegment operator *(const CMatJoint3x4 &transform, const CLineSegment &l)
{
	return CLineSegment(transform.MulPos(l.begin), transform.MulPos(l.end));
}

CLineSegment operator *(const CMat4D &transform, const CLineSegment &l)
{
	return CLineSegment(transform.MulPos(l.begin), transform.MulPos(l.end));
}

CLineSegment operator *(const CQuaternion &transform, const CLineSegment &l)
{
	return CLineSegment(transform * l.begin, transform * l.end);
}


std::string CLineSegment::toString() const
{
	char str[256];
	std::sprintf(str, "CLineSegment(begin:(%.2f, %.2f, %.2f) end:(%.2f, %.2f, %.2f))",
		begin.x, begin.y, begin.z, end.x, end.y, end.z);
	return str;
}

std::ostream &operator <<(std::ostream &o, const CLineSegment &lineSegment)
{
	o << lineSegment.toString();
	return o;
}



#ifdef MATH_GRAPHICSENGINE_INTEROP

void CLineSegment::toLineList(VertexBuffer &vb) const
{
	int startIndex = vb.AppendVertices(2);
	vb.set(startIndex, VDPosition, CVec4D(begin, 1.f));
	vb.set(startIndex+1, VDPosition, CVec4D(end, 1.f));
}

#endif


} //end GEO
}  //end SMF


