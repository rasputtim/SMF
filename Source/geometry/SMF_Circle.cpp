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
namespace MATH {
	class CQuaternion;
}
namespace GEO{


CCircle::CCircle(const CVec3D &center, const CVec3D &n, float radius)
:pos(center),
normal(n),
r(radius)
{
}

CVec3D CCircle::basisU() const
{
	return normal.perpendicular();
}

CVec3D CCircle::basisV() const
{
	return normal.anotherPerpendicular();
}

CVec3D CCircle::getPoint(float angleRadians) const
{
	return pos + r * (CMath::cos(angleRadians) * basisU() + CMath::sin(angleRadians) * basisV());
}

CVec3D CCircle::getPoint(float angleRadians, float d) const
{
	return pos + r * d * (CMath::cos(angleRadians) * basisU() + CMath::sin(angleRadians) * basisV());
}

CVec3D CCircle::extremePoint(const CVec3D &direction) const
{
	CVec3D d = direction - direction.projectToNorm(normal);
	if (d.isZero())
		return pos;
	else
		return pos + d.scaledToLength(r);
}

CPlane CCircle::containingPlane() const
{
	return CPlane( normal, pos);
}

void CCircle::translate(const CVec3D &offset)
{
	pos += offset;
}

void CCircle::transform(const CMat3D &transform)
{
//s	SMF_ASSERT(transform.hasUniformScale());
//s	SMF_ASSERT(transform.isColOrthogonal());
	pos = transform.mul(pos);
	normal = transform.mul(normal).normalized();
	r *= transform.Col(0).getLenght(); // scale the radius of the circle.
}

void CCircle::transform(const CMatJoint3x4 &transform)
{
//s	SMF_ASSERT(transform.hasUniformScale());
//s	SMF_ASSERT(transform.isColOrthogonal());
	pos = transform.MulPos(pos);
	normal = transform.MulDir(normal).normalized();
	r *= transform.Col(0).getLenght(); // scale the radius of the circle.
}

void CCircle::transform(const CMat4D &transform)
{
	//s SMF_ASSERT(transform.hasUniformScale());
	//s SMF_ASSERT(transform.isColOrthogonal3());
	pos = transform.MulPos(pos);
	normal = transform.MulDir(normal).normalized();
	r *= transform.Col3(0).getLenght(); // scale the radius of the circle.
}

void CCircle::transform(const CQuaternion &transform)
{
	pos = transform.mul(pos);
	normal = transform.mul(normal);
}

bool CCircle::edgeContains(const CVec3D &point, float maxDistance) const
{
	return distanceToEdge(point) <= maxDistance;
}
/*
bool CCircle::DiscContains(const CVec3D &point, float maxDistance) const
{
	return distanceToDisc(point) <= maxDistance;
}

*/
float CCircle::distanceToEdge(const CVec3D &point) const
{
	return closestPointToEdge(point).distance(point);
}
/*
float CCircle::distanceToEdge(const CRay &ray, float *d, CVec3D *closestPoint) const
{
	float t;
	CVec3D cp = closestPointToEdge(ray, &t);
	if (closestPoint)
		*closestPoint = cp;
	if (d)
		*d = t;
	return cp.distance(ray.getPoint(t));
}

float CCircle::distanceToEdge(const CLineSegment &lineSegment, float *d, CVec3D *closestPoint) const
{
	float t;
	CVec3D cp = closestPointToEdge(lineSegment, &t);
	if (closestPoint)
		*closestPoint = cp;
	if (d)
		*d = t;
	return cp.distance(lineSegment.getPoint(t));
}

float CCircle::distanceToEdge(const CLine &line, float *d, CVec3D *closestPoint) const
{
	float t;
	CVec3D cp = closestPointToEdge(line, &t);
	if (closestPoint)
		*closestPoint = cp;
	if (d)
		*d = t;
	return cp.distance(line.getPoint(t));
}
*/
float CCircle::distanceToDisc(const CVec3D &point) const
{
	return closestPointToDisc(point).distance(point);
}

CVec3D CCircle::closestPointToEdge(const CVec3D &point) const
{
	CVec3D pointOnPlane = containingPlane().project(point);
	CVec3D diff = pointOnPlane - pos;
	if (diff.isZero())
		return getPoint(0); // The point is in the center of the circle, all points are equally close.
	return pos + diff.scaledToLength(r);
}

CVec3D CCircle::closestPointToDisc(const CVec3D &point) const
{
	CVec3D pointOnPlane = containingPlane().project(point);
	CVec3D diff = pointOnPlane - pos;
	float dist = diff.getLengthSqr();
	if (dist > r*r)
		diff = diff * (r / CMath::sqrt(dist));

	return pos + diff;
}

int CCircle::intersects(const CPlane &plane, CVec3D *pt1, CVec3D *pt2) const
{
	return plane.intersects(*this, pt1, pt2);
}

int CCircle::intersects(const CPlane &plane) const
{
	return plane.intersects(*this);
}

bool CCircle::intersectsDisc(const CLine &line) const
{
	float d;
	bool intersectsPlane = line.intersects(containingPlane(), &d);
	if (intersectsPlane)
		return false;
	return line.getPoint(d).distanceSq(pos) <= r*r;
}

bool CCircle::intersectsDisc(const CLineSegment &lineSegment) const
{
	float d;
	bool intersectsPlane = lineSegment.intersects(containingPlane(), &d);
	if (intersectsPlane)
		return false;
	return lineSegment.getPoint(d).distanceSq(pos) <= r*r;
}

bool CCircle::intersectsDisc(const CRay &ray) const
{
	float d;
	bool intersectsPlane = ray.intersects(containingPlane(), &d);
	if (intersectsPlane)
		return false;
	return ray.getPoint(d).distanceSq(pos) <= r*r;
}


std::vector<CVec3D> CCircle::intersectsFaces(const CAABBox &aabb) const
{
    return intersectsFaces(aabb.toOBBox());
}

std::vector<CVec3D> CCircle::intersectsFaces(const COBBox &obb) const
{
	std::vector<CVec3D> intersectionPoints;
	for(int i = 0; i < 6; ++i)
	{		
		CPlane p = obb.facePlane(i);
		CVec3D pt1, pt2;
		int numIntersections = intersects(p, &pt1, &pt2);
		if (numIntersections >= 1 && obb.contains(pt1))
			intersectionPoints.push_back(pt1);
		if (numIntersections >= 2 && obb.contains(pt2))
			intersectionPoints.push_back(pt2);
	}
	return intersectionPoints;
}
bool CCircle::isFinite() const
{
	return pos.isFinite() && MATH::isFinite(r);
}

bool CCircle::isDegenerate() const
{
	return r < 0.f;
}
std::string CCircle::toString() const
{
	char str[256];
	std::sprintf(str, "CCircle(pos:(%.2f, %.2f, %.2f) normal:(%.2f, %.2f, %.2f), r:%.2f)",
		pos.x, pos.y, pos.z, normal.x, normal.y, normal.z, r);
	return str;
}

std::ostream &operator <<(std::ostream &o, const CCircle &circle)
{
	o << circle.toString();
	return o;
}



CCircle operator *(const CMat3D &transform, const CCircle &circle)
{
	CCircle c(circle);
	c.transform(transform);
	return c;
}

CCircle operator *(const CMatJoint3x4 &transform, const CCircle &circle)
{
	CCircle c(circle);
	c.transform(transform);
	return c;
}

CCircle operator *(const CMat4D &transform, const CCircle &circle)
{
	CCircle c(circle);
	c.transform(transform);
	return c;
}

CCircle operator *(const CQuaternion &transform, const CCircle &circle)
{
	CCircle c(circle);
	c.transform(transform);
	return c;
}



} //end GEO
}  //end SMF
