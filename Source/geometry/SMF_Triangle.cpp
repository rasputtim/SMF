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



CTriangle::CTriangle(const CVec3D &a_, const CVec3D &b_, const CVec3D &c_)
:a(a_), b(b_), c(c_)
{
}

void CTriangle::translate(const CVec3D &offset)
{
	a += offset;
	b += offset;
	c += offset;
}

void CTriangle::transform(const CMat3D &transform)
{
	transform.batchTransform(&a, 3);
}

void CTriangle::transform(const CMatJoint3x4 &transform)
{
	transform.batchTransformPos(&a, 3);
}

void CTriangle::transform(const CMat4D &transform)
{
	a = transform.MulPos(a);
	b = transform.MulPos(b);
	c = transform.MulPos(c);
}

void CTriangle::transform(const CQuaternion &transform)
{
	a = transform * a;
	b = transform * b;
	c = transform * c;
}

/// Implementation from Christer Ericson's Real-Time Collision Detection, pp. 51-52.
inline float TriArea2D(float x1, float y1, float x2, float y2, float x3, float y3)
{
	return (x1-x2)*(y2-y3) - (x2-x3)*(y1-y2);
}

CVec3D CTriangle::barycentricUVW(const CVec3D &point) const
{
	// Implementation from Christer Ericson's Real-Time Collision Detection, pp. 51-52.

	// Unnormalized triangle normal.
	CVec3D m = (b-a).cross(c-a);

	// Nominators and one-over-denominator for u and v ratios.
	float nu, nv, ood;

	// Absolute components for determining projection plane.
	float x = CMath::fabs(m.x);
	float y = CMath::fabs(m.y);
	float z = CMath::fabs(m.z);

	if (x >= y && x >= z)
	{
		// project to the yz plane.
		nu = TriArea2D(point.y, point.z, b.y, b.z, c.y, c.z); // area of PBC in yz-plane.
		nv = TriArea2D(point.y, point.z, c.y, c.z, a.y, a.z); // area OF PCA in yz-plane.
		ood = 1.f / m.x; // 1 / (2*area of ABC in yz plane)
	}
	else if (y >= z) // Note: The book has a redundant 'if (y >= x)' comparison
	{
		// y is largest, project to the xz-plane.
		nu = TriArea2D(point.x, point.z, b.x, b.z, c.x, c.z);
		nv = TriArea2D(point.x, point.z, c.x, c.z, a.x, a.z);
		ood = 1.f / -m.y;
	}
	else // z is largest, project to the xy-plane.
	{
		nu = TriArea2D(point.x, point.y, b.x, b.y, c.x, c.y);
		nv = TriArea2D(point.x, point.y, c.x, c.y, a.x, a.y);
		ood = 1.f / m.z;
	}
	float u = nu * ood;
	float v = nv * ood;
	float w = 1.f - u - v;
	return CVec3D(u,v,w);
#if 0 // TODO: This version should be more SIMD-friendly, but for some reason, it doesn't return good values for all points inside the triangle.
	CVec3D v0 = b - a;
	CVec3D v1 = c - a;
	CVec3D v2 = point - a;
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

CVec2D CTriangle::barycentricUV(const CVec3D &point) const
{
	CVec3D uvw = barycentricUVW(point);
	return CVec2D(uvw.y, uvw.z);
}

bool CTriangle::barycentricInsideTriangleboundingAABB(const CVec3D &barycentric)
{
	return barycentric.x >= 0.f && barycentric.y >= 0.f && barycentric.z >= 0.f &&
		CMath::equalsAbs(barycentric.x + barycentric.y + barycentric.z, 1.f);
}

CVec3D CTriangle::Point(float u, float v) const
{
	return a + (b-a) * u + (c-a) * v;
}

CVec3D CTriangle::Point(float u, float v, float w) const
{
	return u * a + v * b + w * c;
}

CVec3D CTriangle::Point(const CVec3D &b) const
{
	return Point(b.x, b.y, b.z);
}

CVec3D CTriangle::Point(const CVec2D &b) const
{
	return Point(b.x, b.y);
}

CVec3D CTriangle::centroid() const
{
	return (a + b + c) / 3.f;
}

float CTriangle::area() const
{
	return 0.5f * (b-a).cross(c-a).getLenght();
}

float CTriangle::perimeter() const
{
	return a.distance(b) + b.distance(c) + c.distance(a);
}

CLineSegment CTriangle::edge(int i) const
{
	SMF_ASSERT(0 <= i);
	SMF_ASSERT(i <= 2);
	if (i == 0)
		return CLineSegment(a, b);
	else if (i == 1)
		return CLineSegment(b, c);
	else if (i == 2)
		return CLineSegment(c, a);
	else
		return CLineSegment(CVec3D::nan, CVec3D::nan);
}

CVec3D CTriangle::vertex(int i) const
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
		return CVec3D::nan;
}

CPlane CTriangle::planeCCW() const
{
	return CPlane(a, b, c);
}

CPlane CTriangle::planeCW() const
{
	return CPlane(a, c, b);
}

CVec3D CTriangle::normalCCW() const
{
	return unNormalizedNormalCCW().normalized();
}

CVec3D CTriangle::normalCW() const
{
	return unNormalizedNormalCW().normalized();
}

CVec3D CTriangle::unNormalizedNormalCCW() const
{
	return (b-a).cross( c-a);
}

CVec3D CTriangle::unNormalizedNormalCW() const
{
	return (c-a).cross( b-a);
}

CVec3D CTriangle::extremePoint(const CVec3D &direction) const
{
	CVec3D mostExtreme = CVec3D::nan;
	float mostExtremeDist = -FLT_MAX;
	for(int i = 0; i < 3; ++i)
	{
		CVec3D pt = vertex(i);
		float d = (direction* pt);
		if (d > mostExtremeDist)
		{
			mostExtremeDist = d;
			mostExtreme = pt;
		}
	}
	return mostExtreme;
}

CPolygon CTriangle::toPolygon() const
{
	CPolygon p;
	p.p.push_back(a);
	p.p.push_back(b);
	p.p.push_back(c);
	return p;
}

CPolyhedron CTriangle::toPolyhedron() const
{
	return toPolygon().toPolyhedron();
}

CAABBox CTriangle::boundingAABB() const
{
	CAABBox aabb;
	aabb.toNegativeInfinity();
	aabb.enclose(a);
	aabb.enclose(b);
	aabb.enclose(c);
	return aabb;
}

float CTriangle::area2D(const CVec2D &p1, const CVec2D &p2, const CVec2D &p3)
{
	return (p1.x - p2.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p2.y);
}

float CTriangle::signedArea(const CVec3D &pt, const CVec3D &a, const CVec3D &b, const CVec3D &c)
{
	return ((b-pt).cross(c-pt))* ((b-a).cross(c-a).normalized());
}

bool CTriangle::isFinite() const
{
	return a.isFinite() && b.isFinite() && c.isFinite();
}

bool CTriangle::isDegenerate(float epsilon) const
{
	return isDegenerate(a, b, c, epsilon);
}

bool CTriangle::isDegenerate(const CVec3D &a, const CVec3D &b, const CVec3D &c, float epsilon)
{
	return a.compare(b, epsilon) || a.compare(c, epsilon) || b.compare(c, epsilon);
}

bool CTriangle::contains(const CVec3D &point, float triangleThickness) const
{
	if (planeCCW().distance(point) > triangleThickness) // The winding order of the triangle plane does not matter.
		return false; ///\todo The plane-point distance test is omitted in Real-Time Collision Detection. p. 25. A bug in the book?

	CVec3D br = barycentricUVW(point);
	return br.x >= -1e-3f && br.y >= -1e-3f && br.z >= -1e-3f; // Allow for a small epsilon to properly account for points very near the edges of the triangle.
}

bool CTriangle::contains(const CLineSegment &lineSegment, float triangleThickness) const
{
	return contains(lineSegment.begin, triangleThickness) && contains(lineSegment.end, triangleThickness);
}

bool CTriangle::contains(const CTriangle &triangle, float triangleThickness) const
{
	return contains(triangle.a, triangleThickness) && contains(triangle.b, triangleThickness)
	  && contains(triangle.c, triangleThickness);
}

/*
bool CTriangle::contains(const CPolygon &polygon, float triangleThickness) const
{
	if (polygon.points.size() == 0)
		return false;
	for(int i = 0; i < polygon.points.size(); ++i)
		if (!contains(polygon.points[i], triangleThickness))
			return false;
	return true;
}
*/
float CTriangle::distance(const CVec3D &point) const
{
	return closestPoint(point).distance(point);
}

float CTriangle::distance(const CSphere &sphere) const
{
	return MAX(0.f, distance(sphere.getOrigin()) - sphere.getRadius());
}



/** Calculates the intersection between a line and a triangle. The facing is not accounted for, so
	rays are reported to intersect triangles that are both front and backfacing.
	According to "T. M&ouml;ller, B. Trumbore. Fast, Minimum Storage CRay/CTriangle intersection. 2005."
	http://www.graphics.cornell.edu/pubs/1997/MT97.html
	\param linePos The starting point of the line.
	\param lineDir The direction vector of the line. This does not need to be normalized.
	\param v0 vertex 0 of the triangle.
	\param v1 vertex 1 of the triangle.
	\param v2 vertex 2 of the triangle.
	\param u [out] The barycentric u coordinate is returned here if an intersection occurred.
	\param v [out] The barycentric v coordinate is returned here if an intersection occurred.
	\return The distance along the ray to the point of intersection, or +inf if no intersection occurred.
		If no intersection, then u and v and t will contain undefined values. If lineDir was not normalized, then to get the
		real world-space distance, one must scale the returned value with lineDir.getLenght(). If the returned value is negative,
		then the intersection occurs 'behind' the line starting position, with respect to the direction vector lineDir. */
float CTriangle::intersectLineTri(const CVec3D &linePos, const CVec3D &lineDir,
		const CVec3D &v0, const CVec3D &v1, const CVec3D &v2,
		float &u, float &v)
{
	CVec3D vE1, vE2;
	CVec3D vT, vP, vQ;

	const float epsilon =CMath::EPSILON_SuperLow;

	// edge vectors
	vE1 = v1 - v0;
	vE2 = v2 - v0;

	// begin calculating determinant - also used to calculate U parameter
	vP = lineDir.cross(vE2);

	// If det < 0, intersecting backfacing tri, > 0, intersecting frontfacing tri, 0, parallel to plane.
	const float det = vE1*vP;

	// If determinant is near zero, ray lies in plane of triangle.
	if (fabs(det) <= epsilon)
		return CMath::INFINITY_FLOAT;
	const float recipDet = 1.f / det;

	// Calculate distance from v0 to ray origin
	vT = linePos - v0;

	// Output barycentric u
	u = (vT*vP) * recipDet;
	if (u < -epsilon || u > 1.f + epsilon)
		return CMath::INFINITY_FLOAT; // Barycentric U is outside the triangle - early out.

	// Prepare to test V parameter
	vQ = vT.cross(vE1);

	// Output barycentric v
	v = (lineDir*vQ) * recipDet;
	if (v < -epsilon || u + v > 1.f + epsilon) // Barycentric V or the combination of U and V are outside the triangle - no intersection.
		return CMath::INFINITY_FLOAT;

	// Barycentric u and v are in limits, the ray intersects the triangle.
	
	// Output signed distance from ray to triangle.
	return (vE2*vQ) * recipDet;
//	return (det < 0.f) ? IntersectBackface : IntersectFrontface;
}

/// [groupSyntax]
bool CTriangle::intersects(const CLineSegment &l, float *d, CVec3D *intersectionPoint) const
{
	/** The CTriangle-CLine/CLineSegment/CRay intersection tests are based on M&ouml;ller-Trumbore method:
		"T. M&ouml;ller, B. Trumbore. Fast, Minimum Storage CRay/CTriangle intersection. 2005."
		http://jgt.akpeters.com/papers/MollerTrumbore97/. */
	float u, v;
	float t = intersectLineTri(l.begin, l.getDir(), a, b, c, u, v);
	bool success = (t >= 0 && t != CMath::INFINITY_FLOAT);
	if (!success)
		return false;
	float length = l.getLengthSqr();
	if (t < 0.f || t*t >= length)
		return false;
	length = CMath::sqrt(length);
	if (d)
	{
		float len = t / length;
		*d = len;
		if (intersectionPoint)
			*intersectionPoint = l.getPoint(len);
	}
	else if (intersectionPoint)
		*intersectionPoint = l.getPoint(t / length);
	return true;
}

bool CTriangle::intersects(const CLine &l, float *d, CVec3D *intersectionPoint) const
{
	float u, v;
	float t = intersectLineTri(l.pos, l.dir, a, b, c, u, v);
	bool success = (t != CMath::INFINITY_FLOAT);
	if (!success)
		return false;
	if (d)
		*d = t;
	if (intersectionPoint)
		*intersectionPoint = l.getPoint(t);
	return success;
}

bool CTriangle::intersects(const CRay &r, float *d, CVec3D *intersectionPoint) const
{
	float u, v;
	float t = intersectLineTri(r.pos, r.dir, a, b, c, u, v);
	bool success = (t >= 0 && t != CMath::INFINITY_FLOAT);
	if (!success)
		return false;
	if (d)
		*d = t;
	if (intersectionPoint)
		*intersectionPoint = r.getPoint(t);
	return success;
}

bool CTriangle::intersects(const CPlane &plane) const
{
	return plane.intersects(*this);
}

/// [groupSyntax]
/** For CTriangle-CSphere intersection code, see Christer Ericson's Real-Time Collision Detection, p.167. */
bool CTriangle::intersects(const CSphere &sphere, CVec3D *closestPointOnTriangle) const
{
	CVec3D pt = closestPoint(sphere.getOrigin());

	if (closestPointOnTriangle)
		*closestPointOnTriangle = pt;

	return pt.distanceSq(sphere.getOrigin()) <= sphere.getRadius() * sphere.getRadius();
}

bool CTriangle::intersects(const CSphere &sphere) const
{
	return intersects(sphere, 0);
}

static void FindIntersectingLineSegments(const CTriangle &t, float da, float db, float dc, CLineSegment &l1, CLineSegment &l2)
{
	if (da*db > 0.f)
	{
		l1 = CLineSegment(t.a, t.c);
		l2 = CLineSegment(t.b, t.c);
	}
	else if (db*dc > 0.f)
	{
		l1 = CLineSegment(t.a, t.b);
		l2 = CLineSegment(t.a, t.c);
	}
	else
	{
		l1 = CLineSegment(t.a, t.b);
		l2 = CLineSegment(t.b, t.c);
	}
}

/// [groupSyntax]
/** The CTriangle-CTriangle test implementation is based on pseudo-code from Tomas M&ouml;ller's
	"A Fast CTriangle-CTriangle intersection Test": http://jgt.akpeters.com/papers/Moller97/.
	See also Christer Ericson's Real-Time Collision Detection, p. 172. */
bool CTriangle::intersects(const CTriangle &t2, CLineSegment *outLine) const
{
	// Is the triangle t2 completely on one side of the plane of this triangle?
	CPlane p1 = this->planeCCW();
	float t2da = p1.signedDistance(t2.a);
	float t2db = p1.signedDistance(t2.b);
	float t2dc = p1.signedDistance(t2.c);
	if (t2da*t2db > 0.f && t2da*t2dc > 0.f)
		return false;
	// Is this triangle completely on one side of the plane of the triangle t2?
	CPlane p2 = t2.planeCCW();
	float t1da = p2.signedDistance(this->a);
	float t1db = p2.signedDistance(this->b);
	float t1dc = p2.signedDistance(this->c);
	if (t1da*t1db > 0.f && t1da*t1dc > 0.f)
		return false;

	// find the intersection line of the two planes.
	CLine l;
	bool success = p1.intersects(p2, &l);
	SMF_ASSERT(success); // We already determined the two triangles have intersecting planes, so this should always succeed.
	if (!success)
		return false;

	// find the two line segments of both triangles which straddle the intersection line.
	CLineSegment l1a, l1b;
	CLineSegment l2a, l2b;
	FindIntersectingLineSegments(*this, t1da, t1db, t1dc, l1a, l1b);
	FindIntersectingLineSegments(t2, t2da, t2db, t2dc, l2a, l2b);

	// find the projection intervals on the intersection line.
	float d1a, d1b, d2a, d2b;
	l.distance(l1a, &d1a);
	l.distance(l1b, &d1b);
	l.distance(l2a, &d2a);
	l.distance(l2b, &d2b);
	if (d1a > d1b)
		Swap(d1a, d1b);
	if (d2a > d2b)
		Swap(d2a, d2b);
	float rStart = MAX(d1a, d2a);
	float rEnd = MIN(d1b, d2b);
	if (rStart <= rEnd)
	{
		if (outLine)
			*outLine = CLineSegment(l.getPoint(rStart), l.getPoint(rEnd));
		return true;
	}
	return false;
}

bool RangesOverlap(float start1, float end1, float start2, float end2)
{
	return end1 >= start2 && end2 >= start1;
}

/// [groupSyntax]
bool CTriangle::intersects(const CAABBox &aabb) const
{
/** The CAABBox-CTriangle test implementation is based on the pseudo-code in
	Christer Ericson's Real-Time Collision Detection, pp. 169-172. */
	///\todo The CTriangle-CAABBox intersection test can be greatly optimized by manually unrolling loops, trivial math and by avoiding
	/// unnecessary copying.
	float t1, t2, a1, a2;
	const CVec3D e[3] = { CVec3D(1,0,0), CVec3D(0,1,0), CVec3D(0,0,1) };

	for(int i = 0; i < 3; ++i)
	{
		projectToAxis(e[i], t1, t2);
		aabb.projectToAxis(e[i], a1, a2);
		if (!RangesOverlap(t1, t2, a1, a2))
			return false;
	}

	CVec3D n = unNormalizedNormalCCW();
	projectToAxis(n, t1, t2);
	aabb.projectToAxis(n, a1, a2);
	if (!RangesOverlap(t1, t2, a1, a2))
		return false;

	const CVec3D t[3] = { b-a, c-a, c-b };

	for(int i = 0; i < 3; ++i)
		for(int j = 0; j < 3; ++j)
		{
			CVec3D axis = e[i].cross( t[j]);
			float len = axis.getLengthSqr();
			if (len <= 1e-4f)
				continue; // Ignore tests on degenerate axes.

			projectToAxis(axis, t1, t2);
			aabb.projectToAxis(axis, a1, a2);
			if (!RangesOverlap(t1, t2, a1, a2))
				return false;
		}

	// No separating axis exists, the CAABBox and triangle intersect.
	return true;
}

bool CTriangle::intersects(const COBBox &obb) const
{
	return obb.intersects(*this);
}

bool CTriangle::intersects(const CPolygon &polygon) const
{
	return polygon.intersects(*this);
}
#if 0
bool CTriangle::intersects(const Frustum &frustum) const
{
	return frustum.intersects(*this);
}
#endif
bool CTriangle::intersects(const CPolyhedron &polyhedron) const
{
	return polyhedron.intersects(*this);
}
#if 0
bool CTriangle::intersects(const Capsule &capsule) const
{
	return capsule.intersects(*this);
}
#endif
void CTriangle::projectToAxis(const CVec3D &axis, float &dMin, float &dMax) const
{
	dMin = dMax = (axis* a);
	float t = (axis* b);
	dMin = MIN(t, dMin);
	dMax = MAX(t, dMax);
	t = (axis* c);
	dMin = MIN(t, dMin);
	dMax = MAX(t, dMax);
}

/// [groupSyntax]
CVec3D CTriangle::closestPoint(const CVec3D &p) const
{
	/** The code for CTriangle-CVec3D test is from Christer Ericson's Real-Time Collision Detection, pp. 141-142. */

	// Check if P is in vertex region outside A.
	CVec3D ab = b - a;
	CVec3D ac = c - a;
	CVec3D ap = p - a;
	float d1 = (ab* ap);
	float d2 = (ac* ap);
	if (d1 <= 0.f && d2 <= 0.f)
		return a; // Barycentric coordinates are (1,0,0).

	// Check if P is in vertex region outside B.
	CVec3D bp = p - b;
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
	CVec3D cp = p - c;
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
	return a + ab * v + ac * w;
}

CVec3D CTriangle::closestPoint(const CLineSegment &lineSegment, CVec3D *otherPt) const
{
	///\todo Optimize.
	CVec3D intersectionPoint;
	if (intersects(lineSegment, 0, &intersectionPoint))
	{
		if (otherPt)
			*otherPt = intersectionPoint;
		return intersectionPoint;
	}

	float u1,v1,d1;
	CVec3D pt1 = closestPointToTriangleEdgesignedArea(lineSegment, &u1, &v1, &d1);

	CVec3D pt2 = closestPoint(lineSegment.begin);
	CVec3D pt3 = closestPoint(lineSegment.end);
	
	float D1 = pt1.distanceSq(lineSegment.getPoint(d1));
	float D2 = pt2.distanceSq(lineSegment.begin);
	float D3 = pt3.distanceSq(lineSegment.end);

	if (D1 <= D2 && D1 <= D3)
	{
		if (otherPt)
			*otherPt = lineSegment.getPoint(d1);
		return pt1;
	}
	else if (D2 <= D3)
	{
		if (otherPt)
			*otherPt = lineSegment.begin;
		return pt2;
	}
	else
	{
		if (otherPt)
			*otherPt = lineSegment.end;
		return pt3;
	}
}

#if 0
///\todo Enable this codepath. This if rom Geometric Tools for Computer Graphics,
/// but the algorithm in the book is broken and does not take into account the
/// direction of the gradient to determine the proper region of intersection.
/// Instead using a slower code path above.
/// [groupSyntax]
CVec3D CTriangle::closestPoint(const CLineSegment &lineSegment, CVec3D *otherPt) const
{
	CVec3D e0 = b - a;
	CVec3D e1 = c - a;
	CVec3D v_p = a - lineSegment.begin;
	CVec3D d = lineSegment.end - lineSegment.begin;

	// Q(u,v) = a + u*e0 + v*e1
	// L(t)   = ls.a + t*d
	// Minimize the distance |Q(u,v) - L(t)|^2 under u >= 0, v >= 0, u+v <= 1, t >= 0, t <= 1.

	float v_p_dot_e0 = dot(v_p, e0);
	float v_p_dot_e1 = dot(v_p, e1);
	float v_p_dot_d = dot(v_p, d);

	CMat3D m;
	m[0][0] = dot(e0, e0); m[0][1] = dot(e0, e1); m[0][2] = -dot(e0, d);
	m[1][0] =     m[0][1]; m[1][1] = dot(e1, e1); m[1][2] = -dot(e1, d);
	m[2][0] =     m[0][2]; m[2][1] =     m[1][2]; m[2][2] =  dot(d, d);

	CVec3D B(-v_p_dot_e0, -v_p_dot_e1, v_p_dot_d);

	CVec3D uvt;
	bool success = m.solveAxb(B, uvt);
	if (!success)
	{
		float t1, t2, t3;
		float s1, s2, s3;
		CLineSegment e1 = edge(0);
		CLineSegment e2 = edge(1);
		CLineSegment e3 = edge(2);
		float d1 = e1.distance(lineSegment, &t1, &s1);
		float d2 = e2.distance(lineSegment, &t2, &s2);
		float d3 = e3.distance(lineSegment, &t3, &s3);
		if (d1 < d2 && d1 < d3)
		{
			if (otherPt)
				*otherPt = lineSegment.getPoint(s1);
			return e1.getPoint(t1);
		}
		else if (d2 < d3)
		{
			if (otherPt)
				*otherPt = lineSegment.getPoint(s2);
			return e2.getPoint(t2);
		}
		else
		{
			if (otherPt)
				*otherPt = lineSegment.getPoint(s3);
			return e3.getPoint(t3);
		}
	}

	if (uvt.x < 0.f)
	{
		// clamp to u == 0 and solve again.
		float m_00 = m[2][2];
		float m_01 = -m[1][2];
		float m_10 = -m[2][1];
		float m_11 = m[1][1];
		float det = m_00 * m_11 - m_01 * m_10;
		float v = m_00 * B[1] + m_01 * B[2];
		float t = m_10 * B[1] + m_11 * B[2];
		v /= det;
		t /= det;
		if (v < 0.f)
		{
			// clamp to v == 0 and solve for t.
			t = B[2] / m[2][2];
			t = clamp01(t); // The solution for t must also be in the range [0,1].
			// The solution is (u,v,t)=(0,0,t).
			if (otherPt)
				*otherPt = lineSegment.getPoint(t);
			return a;
		}
		else if (v > 1.f)
		{
			// clamp to v == 1 and solve for t.
			t = (B[2] - m[2][1]) / m[2][2];
			t = clamp01(t);
			// The solution is (u,v,t)=(0,1,t).
			if (otherPt)
				*otherPt = lineSegment.getPoint(t);
			return c; // == a + v*e1
		}
		else if (t < 0.f)
		{
			// clamp to t == 0 and solve for v.
			v = B[1] / m[1][1];
//			SMF_ASSERT(equalsAbs(v, clamp01(v)));
			v = clamp01(v); // The solution for v must also be in the range [0,1]. TODO: Is this guaranteed by the above?
			// The solution is (u,v,t)=(0,v,0).
			if (otherPt)
				*otherPt = lineSegment.begin;
			return a + v * e1;
		}
		else if (t > 1.f)
		{
			// clamp to t == 1 and solve for v.
			v = (B[1] - m[1][2]) / m[1][1];
//			SMF_ASSERT(equalsAbs(v, clamp01(v)));
			v = clamp01(v); // The solution for v must also be in the range [0,1]. TODO: Is this guaranteed by the above?
			// The solution is (u,v,t)=(0,v,1).
			if (otherPt)
				*otherPt = lineSegment.end;
			return a + v * e1;
		}
		else
		{
			// The solution is (u,v,t)=(0,v,t).
			if (otherPt)
				*otherPt = lineSegment.getPoint(t);
			return a + v * e1;
		}
	}
	else if (uvt.y < 0.f)
	{
		// clamp to v == 0 and solve again.
		float m_00 = m[2][2];
		float m_01 = -m[0][2];
		float m_10 = -m[2][0];
		float m_11 = m[0][0];
		float det = m_00 * m_11 - m_01 * m_10;
		float u = m_00 * B[0] + m_01 * B[2];
		float t = m_10 * B[0] + m_11 * B[2];
		u /= det;
		t /= det;

		if (u < 0.f)
		{
			// clamp to u == 0 and solve for t.
			t = B[2] / m[2][2];
			t = clamp01(t); // The solution for t must also be in the range [0,1].
			// The solution is (u,v,t)=(0,0,t).
			if (otherPt)
				*otherPt = lineSegment.getPoint(t);
			return a;
		}
		else if (u > 1.f)
		{
			// clamp to u == 1 and solve for t.
			t = (B[2] - m[2][0]) / m[2][2];
			t = clamp01(t); // The solution for t must also be in the range [0,1].
			// The solution is (u,v,t)=(1,0,t).
			if (otherPt)
				*otherPt = lineSegment.getPoint(t);
			return b;
		}
		else if (t < 0.f)
		{
			// clamp to t == 0 and solve for u.
			u = B[0] / m[0][0];
//			SMF_ASSERT(equalsAbs(u, clamp01(u)));
			u = clamp01(u); // The solution for u must also be in the range [0,1].
			if (otherPt)
				*otherPt = lineSegment.begin;
			return a + u * e0;
		}
		else if (t > 1.f)
		{
			// clamp to t == 1 and solve for u.
			u = (B[0] - m[0][2]) / m[0][0];
//			SMF_ASSERT(equalsAbs(u, clamp01(u)));
			u = clamp01(u); // The solution for u must also be in the range [0,1].
			if (otherPt)
				*otherPt = lineSegment.end;
			return a + u * e0;
		}
		else
		{
			// The solution is (u, 0, t).
			if (otherPt)
				*otherPt = lineSegment.getPoint(t);
			return a + u * e0;
		}
	}
	else if (uvt.z < 0.f)
	{
		if (otherPt)
			*otherPt = lineSegment.begin;
		// clamp to t == 0 and solve again.
		float m_00 = m[1][1];
		float m_01 = -m[0][1];
		float m_10 = -m[1][0];
		float m_11 = m[0][0];
		float det = m_00 * m_11 - m_01 * m_10;
		float u = m_00 * B[0] + m_01 * B[1];
		float v = m_10 * B[0] + m_11 * B[1];
		u /= det;
		v /= det;
		if (u < 0.f)
		{
			// clamp to u == 0 and solve for v.
			v = B[1] / m[1][1];
			v = clamp01(v);
			return a + v*e1;
		}
		else if (v < 0.f)
		{
			// clamp to v == 0 and solve for u.
			u = B[0] / m[0][0];
			u = clamp01(u);
			return a + u*e0;
		}
		else if (u+v > 1.f)
		{
			// set v = 1-u and solve again.
//			u = (B[0] - m[0][0]) / (m[0][0] - m[0][1]);
//			SMF_ASSERT(equalsAbs(u, clamp01(u)));
//			u = clamp01(u); // The solution for u must also be in the range [0,1].
//			return a + u*e0;

			// clamp to v = 1-u and solve again.
			float m_00 = m[2][2];
			float m_01 = m[1][2] - m[0][2];
			float m_10 = m_01;
			float m_11 = m[0][0] + m[1][1] - 2.f * m[0][1];
			float det = m_00 * m_11 - m_01 * m_10;
			float b0 = m[1][1] - m[0][1] + v_p_dot_e1 - v_p_dot_e0;
			float b1 = -m[1][2] + v_p_dot_d;
			float u = m_00 * b0 + m_01 * b1;
			u /= det;
			u = clamp01(u);

			float t = m_10 * b0 + m_11 * b1;
			t /= det;
			t = clamp01(t);
			if (otherPt)
				*otherPt = lineSegment.getPoint(t);
			return a + u*e0 + (1.f-u)*e1;
		}
		else
		{
			// The solution is (u, v, 0)
			return a + u * e0 + v * e1;
		}
	}
	else if (uvt.z > 1.f)
	{
		if (otherPt)
			*otherPt = lineSegment.end;
		// clamp to t == 1 and solve again.
		float m_00 = m[1][1];
		float m_01 = -m[0][1];
		float m_10 = -m[1][0];
		float m_11 = m[0][0];
		float det = m_00 * m_11 - m_01 * m_10;
		float u = m_00 * (B[0]-m[0][2]) + m_01 * (B[1]-m[1][2]);
		float v = m_10 * (B[0]-m[0][2]) + m_11 * (B[1]-m[1][2]);
		u /= det;
		v /= det;
		if (u < 0.f)
		{
			// clamp to u == 0 and solve again.
			v = (B[1] - m[1][2]) / m[1][1];
			v = clamp01(v);
			return a + v*e1;
		}
		else if (u > 1.f)
		{
			// clamp to u == 1 and solve again.
			v = (B[1] - m[1][0] - m[1][2]) / m[1][1];
			v = clamp01(v); // The solution for v must also be in the range [0,1]. TODO: Is this guaranteed by the above?
			// The solution is (u,v,t)=(1,v,1).
			return a + e0 + v*e1;
		}
		else if (u+v > 1.f)
		{
			// set v = 1-u and solve again.

			// Q(u,1-u) = a + u*e0 + e1 - u*e1 = a+e1 + u*(e0-e1)
			// L(1)   = ls.a + t*d = ls.b
			// Minimize the distance |Q(u,1-u) - L(1)| = |a+e1+ls.b + u*(e0-e1)|

			// |K + u*(e0-e1)|^2 = (K,K) + 2*u(K,e0-e1) + u^2 * (e0-e1,e0-e1)

			// grad = 2*(K,e0-e1) + 2*u*(e0-e1,e0-e1) == 0
			//                                      u == (K,e1-e0) / (e0-e1,e0-e1)

			u = (B[0] - m[0][1] - m[0][2]) / (m[0][0] - m[0][1]);
//			u = dot(a + e1 + lineSegment.end, e1 - e0) / dot(e0-e1, e0-e1);

//			SMF_ASSERT(equalsAbs(u, clamp01(u)));
			u = clamp01(u);
			return a + u*e0 + (1-u)*e1;
		}
		else
		{
			// The solution is (u, v, 1)
			return a + u*e0 + v*e1;
		}
	}
	else if (uvt.x + uvt.y > 1.f)
	{
		// clamp to v = 1-u and solve again.
		float m_00 = m[2][2];
		float m_01 = m[1][2] - m[0][2];
		float m_10 = m_01;
		float m_11 = m[0][0] + m[1][1] - 2.f * m[0][1];
		float det = m_00 * m_11 - m_01 * m_10;
		float b0 = m[1][1] - m[0][1] + v_p_dot_e1 - v_p_dot_e0;
		float b1 = -m[1][2] + v_p_dot_d;
		float u = m_00 * b0 + m_01 * b1;
		float t = m_10 * b0 + m_11 * b1;
		u /= det;
		t /= det;

		t = clamp01(t);
		if (otherPt)
			*otherPt = lineSegment.getPoint(t);

		if (u < 0.f)
		{
			// The solution is (u,v,t)=(0,1,t)
			return c;
		}
		if (u > 1.f)
		{
			// The solution is (u,v,t)=(1,0,t)
			return b;
		}
		SMF_ASSERT(t >= 0.f);
		SMF_ASSERT(t <= 1.f);
		return a + u*e0 + (1.f-u)*e1;
	}
	else // All parameters are within range, so the triangle and the line segment intersect, and the intersection point is the closest point.
	{
		if (otherPt)
			*otherPt = lineSegment.getPoint(uvt.z);
		return a + uvt.x * e0 + uvt.y * e1;
	}
}
#endif
CVec3D CTriangle::closestPointToTriangleEdgesignedArea(const CLine &other, float *outU, float *outV, float *outD) const
{
	///\todo Optimize!
	// The line is parallel to the triangle.
	float d1, d2, d3;
	CVec3D pt1 = edge(0).closestPoint(other, 0, &d1);
	CVec3D pt2 = edge(1).closestPoint(other, 0, &d2);
	CVec3D pt3 = edge(2).closestPoint(other, 0, &d3);
	float dist1 = pt1.distanceSq(other.getPoint(d1));
	float dist2 = pt2.distanceSq(other.getPoint(d2));
	float dist3 = pt3.distanceSq(other.getPoint(d3));
	if (dist1 <= dist2 && dist1 <= dist3)
	{
		if (outU) *outU = barycentricUV(pt1).x;
		if (outV) *outV = barycentricUV(pt1).y;
		if (outD) *outD = d1;
		return pt1;
	}
	else if (dist2 <= dist3)
	{
		if (outU) *outU = barycentricUV(pt2).x;
		if (outV) *outV = barycentricUV(pt2).y;
		if (outD) *outD = d2;
		return pt2;
	}
	else
	{
		if (outU) *outU = barycentricUV(pt3).x;
		if (outV) *outV = barycentricUV(pt3).y;
		if (outD) *outD = d3;
		return pt3;
	}
}

CVec3D CTriangle::closestPointToTriangleEdgesignedArea(const CLineSegment &lineSegment, float *outU, float *outV, float *outD) const
{
	///\todo Optimize!
	// The line is parallel to the triangle.
	float d1, d2, d3;
	CVec3D pt1 = edge(0).closestPoint(lineSegment, 0, &d1);
	CVec3D pt2 = edge(1).closestPoint(lineSegment, 0, &d2);
	CVec3D pt3 = edge(2).closestPoint(lineSegment, 0, &d3);
	float dist1 = pt1.distanceSq(lineSegment.getPoint(d1));
	float dist2 = pt2.distanceSq(lineSegment.getPoint(d2));
	float dist3 = pt3.distanceSq(lineSegment.getPoint(d3));
	if (dist1 <= dist2 && dist1 <= dist3)
	{
		if (outU) *outU = barycentricUV(pt1).x;
		if (outV) *outV = barycentricUV(pt1).y;
		if (outD) *outD = d1;
		return pt1;
	}
	else if (dist2 <= dist3)
	{
		if (outU) *outU = barycentricUV(pt2).x;
		if (outV) *outV = barycentricUV(pt2).y;
		if (outD) *outD = d2;
		return pt2;
	}
	else
	{
		if (outU) *outU = barycentricUV(pt3).x;
		if (outV) *outV = barycentricUV(pt3).y;
		if (outD) *outD = d3;
		return pt3;
	}
}

CVec3D CTriangle::closestPoint(const CLine &line, CVec3D *otherPt) const
{
	///\todo Optimize this function.
	CVec3D intersectionPoint;
	if (intersects(line, 0, &intersectionPoint))
	{
		if (otherPt)
			*otherPt = intersectionPoint;
		return intersectionPoint;
	}

	float u1,v1,d1;
	CVec3D pt1 = closestPointToTriangleEdgesignedArea(line, &u1, &v1, &d1);
	if (otherPt)
		*otherPt = line.getPoint(d1);
	return pt1;
}

#if 0
///\todo Enable this codepath. This if rom Geometric Tools for Computer Graphics,
/// but the algorithm in the book is broken and does not take into account the
/// direction of the gradient to determine the proper region of intersection.
/// Instead using a slower code path above.
CVec3D CTriangle::closestPoint(const CLine &line, CVec3D *otherPt) const
{
	CVec3D e0 = b - a;
	CVec3D e1 = c - a;
	CVec3D v_p = a - line.pos;
	CVec3D d = line.dir;

	float v_p_dot_e0 = dot(v_p, e0);
	float v_p_dot_e1 = dot(v_p, e1);
	float v_p_dot_d = dot(v_p, d);

	CMat3D m;
	m[0][0] = dot(e0, e0); m[0][1] = dot(e0, e1); m[0][2] = -dot(e0, d);
	m[1][0] =     m[0][1]; m[1][1] = dot(e1, e1); m[1][2] = -dot(e1, d);
	m[2][0] =     m[0][2]; m[2][1] =     m[1][2]; m[2][2] =  dot(d, d);

	CVec3D B(-v_p_dot_e0, -v_p_dot_e1, v_p_dot_d);

	CVec3D uvt;
	bool success = m.solveAxb(B, uvt);
	if (!success)
	{
		float t1, t2, t3;
		float s1, s2, s3;
		CLineSegment e1 = edge(0);
		CLineSegment e2 = edge(1);
		CLineSegment e3 = edge(2);
		float d1 = e1.distance(line, &t1, &s1);
		float d2 = e2.distance(line, &t2, &s2);
		float d3 = e3.distance(line, &t3, &s3);
		if (d1 < d2 && d1 < d3)
		{
			if (otherPt)
				*otherPt = line.getPoint(s1);
			return e1.getPoint(t1);
		}
		else if (d2 < d3)
		{
			if (otherPt)
				*otherPt = line.getPoint(s2);
			return e2.getPoint(t2);
		}
		else
		{
			if (otherPt)
				*otherPt = line.getPoint(s3);
			return e3.getPoint(t3);
		}
	}

	if (uvt.x < 0.f)
	{
		// clamp to u == 0 and solve again.
		float m_00 = m[2][2];
		float m_01 = -m[1][2];
		float m_10 = -m[2][1];
		float m_11 = m[1][1];
		float det = m_00 * m_11 - m_01 * m_10;
		float v = m_00 * B[1] + m_01 * B[2];
		float t = m_10 * B[1] + m_11 * B[2];
		v /= det;
		t /= det;
		if (v < 0.f)
		{
			// clamp to v == 0 and solve for t.
			t = B[2] / m[2][2];
			// The solution is (u,v,t)=(0,0,t).
			if (otherPt)
				*otherPt = line.getPoint(t);
			return a;
		}
		else if (v > 1.f)
		{
			// clamp to v == 1 and solve for t.
			t = (B[2] - m[2][1]) / m[2][2];
			// The solution is (u,v,t)=(0,1,t).
			if (otherPt)
				*otherPt = line.getPoint(t);
			return c; // == a + v*e1
		}
		else
		{
			// The solution is (u,v,t)=(0,v,t).
			if (otherPt)
				*otherPt = line.getPoint(t);
			return a + v * e1;
		}
	}
	else if (uvt.y < 0.f)
	{
		// clamp to v == 0 and solve again.
		float m_00 = m[2][2];
		float m_01 = -m[0][2];
		float m_10 = -m[2][0];
		float m_11 = m[0][0];
		float det = m_00 * m_11 - m_01 * m_10;
		float u = m_00 * B[0] + m_01 * B[2];
		float t = m_10 * B[0] + m_11 * B[2];
		u /= det;
		t /= det;

		if (u < 0.f)
		{
			// clamp to u == 0 and solve for t.
			t = B[2] / m[2][2];
			// The solution is (u,v,t)=(0,0,t).
			if (otherPt)
				*otherPt = line.getPoint(t);
			return a;
		}
		else if (u > 1.f)
		{
			// clamp to u == 1 and solve for t.
			t = (B[2] - m[2][0]) / m[2][2];
			// The solution is (u,v,t)=(1,0,t).
			if (otherPt)
				*otherPt = line.getPoint(t);
			return b;
		}
		else
		{
			// The solution is (u, 0, t).
			if (otherPt)
				*otherPt = line.getPoint(t);
			return a + u * e0;
		}
	}
	else if (uvt.x + uvt.y > 1.f)
	{
		// clamp to v = 1-u and solve again.
		float m_00 = m[2][2];
		float m_01 = m[1][2] - m[0][2];
		float m_10 = m_01;
		float m_11 = m[0][0] + m[1][1] - 2.f * m[0][1];
		float det = m_00 * m_11 - m_01 * m_10;
		float b0 = m[1][1] - m[0][1] + v_p_dot_e1 - v_p_dot_e0;
		float b1 = -m[1][2] + v_p_dot_d;
		float u = m_00 * b0 + m_01 * b1;
		float t = m_10 * b0 + m_11 * b1;
		u /= det;
		t /= det;

		if (otherPt)
			*otherPt = line.getPoint(t);

		if (u < 0.f)
		{
			// The solution is (u,v,t)=(0,1,t)
			return c;
		}
		if (u > 1.f)
		{
			// The solution is (u,v,t)=(1,0,t)
			return b;
		}
		return a + u*e0 + (1.f-u)*e1;
	}
	else // All parameters are within range, so the triangle and the line segment intersect, and the intersection point is the closest point.
	{
		if (otherPt)
			*otherPt = line.getPoint(uvt.z);
		return a + uvt.x * e0 + uvt.y * e1;
	}
}
#endif

#if 0
/// [groupSyntax]
CVec3D CTriangle::closestPoint(const CLine &other, float *outU, float *outV, float *outD) const
{
	/** The implementation of the CTriangle-CLine test is based on the pseudo-code in
		Schneider, Eberly. Geometric Tools for Computer Graphics pp. 433 - 441. */
	///\todo The CTriangle-CLine code is currently untested. Run tests to ensure the following code works properly.

	// Point on triangle: T(u,v) = a + u*b + v*c;
	// Point on line:  L(t) = p + t*d;
	// Minimize the function Q(u,v,t) = ||T(u,v) - L(t)||.

	CVec3D e0 = b-a;
	CVec3D e1 = c-a;
	CVec3D d = other.dir;

	const float d_e0e0 = dot(e0, e0);
	const float d_e0e1 = dot(e0, e1);
	const float d_e0d = dot(e0, d);
	const float d_e1e1 = dot(e1, e1);
	const float d_e1d = dot(e1, d);
	const float d_dd = dot(d, d);

	CMat3D m;
	m[0][0] = d_e0e0;  m[0][1] = d_e0e1;  m[0][2] = -d_e0d;
	m[1][0] = d_e0e1;  m[1][1] = d_e1e1;  m[1][2] = -d_e1d;
	m[2][0] = -d_e0d;  m[2][1] = -d_e1d;  m[2][2] =   d_dd;

	///\todo add optimized CMat3D::InverseSymmetric().
	bool inv = m.Inverse();
	if (!inv)
		return closestPointToTriangleEdgesignedArea(other, outU, outV, outD);

	CVec3D v_m_p = a - other.pos;
	float v_m_p_e0 = v_m_p.dot(e0);
	float v_m_p_e1 = v_m_p.dot(e1);
	float v_m_p_d = v_m_p.dot(d);
	CVec3D b = CVec3D(-v_m_p_e0, -v_m_p_e1, v_m_p_d);
	CVec3D uvt = m * b;
	// We cannot simply clamp the solution to (uv) inside the constraints, since the expression we
	// are minimizing is quadratic.
	// So, examine case-by-case which part of the space the solution lies in. Because the function is convex,
	// we can clamp the search space to the boundary planes.
	float u = uvt.x;
	float v = uvt.y;
	float t = uvt.z;
	if (u <= 0)
	{
		if (outU) *outU = 0;

		// solve 2x2 matrix for the (v,t) solution when u == 0.
		v = m[1][1]*b[1] + m[1][2]*b[2];
		t = m[2][1]*b[1] + m[2][2]*b[2];

		// Check if the solution is still out of bounds.
		if (v <= 0)
		{
			if (outV) *outV = 0;
			if (outD) *outD = v_m_p_d / d_dd;
			return Point(0, 0);
		}
		else if (v >= 1)
		{
			if (outV) *outV = 1;
			if (outD) *outD = (v_m_p_d - d_e1d) / d_dd;
			return Point(0, 1);
		}
		else // (0 <= v <= 1).
		{
			if (outV) *outV = v;
			if (outD) *outD = t;
			return Point(0, v);
		}
	}
	else if (v <= 0)
	{
		if (outV) *outV = 0;

		// solve 2x2 matrix for the (u,t) solution when v == 0.
		u = m[0][0]*b[0] + m[0][2]*b[2];
		t = m[2][0]*b[0] + m[2][2]*b[2];

		// Check if the solution is still out of bounds.
		if (u <= 0)
		{
			if (outU) *outU = 0;
			if (outD) *outD = v_m_p_d / d_dd;
			return Point(0, 0);
		}
		else if (u >= 1)
		{
			if (outU) *outU = 1;
			if (outD) *outD = (v_m_p_d - d_e0d) / d_dd;
			return Point(1, 0);
		}
		else // (0 <= u <= 1).
		{
			if (outU) *outU = u;
			if (outD) *outD = t;
			return Point(u, 0);
		}
	}
	else if (u + v >= 1.f)
	{
		// set v = 1-u.
#if 0
		float m00 = d_e0e0 + d_e1e1 - 2.f * d_e0e1;
		float m01 = -d_e0d + d_e1d;
		float m10 = -d_e0d + d_e1d;
		float m11 = d_dd;
//		float det = 1.f / (m00*m11 - m01*m10);

		float b0 = d_e1e1 - d_e0e1 + v_m_p_e0 - v_m_p_e1;
		float b1 = d_e1d + v_m_p_d;
		/*
		// Inverse 2x2 matrix.
		Swap(m00, m11);
		Swap(m01, m10);
		m01 = -m01;
		m10 = -m10;
		*/
		// 2x2 * 2 matrix*vec mul.
		u = (m00*b0 + m01*b1);// * det;
		t = (m10*b0 + m11*b1);// * det;
#endif
//		u = m[0][0]*b[0] +

		// Check if the solution is still out of bounds.
		if (u <= 0)
		{
			if (outU) *outU = 0;
			if (outV) *outV = 1;
			if (outD) *outD = (d_e1d + v_m_p_d) / d_dd;
			return Point(0, 1);
		}
		else if (u >= 1)
		{
			if (outU) *outU = 1;
			if (outV) *outV = 0;
			if (outD) *outD = (v_m_p_d + d_e0d) / d_dd;
			return Point(1, 0);
		}
		else // (0 <= u <= 1).
		{
			if (outU) *outU = u;
			if (outV) *outV = 1.f - u;
			if (outD) *outD = t;
			return Point(u, 1.f - u);
		}
	}
	else // Each u, v and t are in appropriate range.
	{
		if (outU) *outU = u;
		if (outV) *outV = v;
		if (outD) *outD = t;
		return Point(u, v);
	}
}
#endif

/// [groupSyntax]
CVec3D CTriangle::closestPoint(const CTriangle &other, CVec3D *otherPt) const
{
	/** The code for computing the closest point pair on two Triangles is based
		on pseudo-code from Christer Ericson's Real-Time Collision Detection, pp. 155-156. */

	// First detect if the two triangles are intersecting.
	CLineSegment l;
	bool success = this->intersects(other, &l);
	if (success)
	{
		CVec3D cp = l.centerPoint();
		if (otherPt)
			*otherPt = cp;
		return cp;
	}

	CVec3D closestThis = this->closestPoint(other.a);
	CVec3D closestOther = other.a;
	float closestDSq = closestThis.distanceSq(closestOther);

	CVec3D pt = this->closestPoint(other.b);
	float dSq = pt.distanceSq(other.b);
	if (dSq < closestDSq) closestThis = pt, closestOther = other.b, closestDSq = dSq;

	pt = this->closestPoint(other.c);
	dSq = pt.distanceSq(other.c);
	if (dSq < closestDSq) closestThis = pt, closestOther = other.c, closestDSq = dSq;

	CLineSegment l1[3] = { CLineSegment(a,b), CLineSegment(a,c), CLineSegment(b,c) };
	CLineSegment l2[3] = { CLineSegment(other.a,other.b), CLineSegment(other.a,other.c), CLineSegment(other.b,other.c) };
	float d, d2;
	for(int i = 0; i < 3; ++i)
		for(int j = 0; j < 3; ++j)
		{
			float dist = l1[i].distance(l2[j], &d, &d2);
			if (dist*dist < closestDSq)
			{
				closestThis = l1[i].getPoint(d);
				closestOther = l2[j].getPoint(d2);
				closestDSq = dist*dist;
			}
		}

	if (otherPt)
		*otherPt = closestOther;
	return closestThis;
}

CVec3D CTriangle::randomPointInside(CRandomLCG &rng) const
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
#ifdef TEST_FOR_CORRECTNESS
	CVec3D pt = Point(s, t);
	CVec2D uv = barycentricUV(pt);
	SMF_ASSERT(uv.x >= 0.f);
	SMF_ASSERT(uv.y >= 0.f);
	SMF_ASSERT(uv.x + uv.y <= 1.f);
	CVec3D uvw = barycentricUVW(pt);
	SMF_ASSERT(uvw.x >= 0.f);
	SMF_ASSERT(uvw.y >= 0.f);
	SMF_ASSERT(uvw.z >= 0.f);
	SMF_ASSERT(CMath::equalsAbs(uvw.x + uvw.y + uvw.z, 1.f));
#endif
	return Point(s, t);
}

CVec3D CTriangle::randomVertex(CRandomLCG &rng) const
{
	return vertex(rng.getInt(0, 2));
}

CVec3D CTriangle::randomPointOnEdge(CRandomLCG &rng) const
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

CTriangle operator *(const CMat3D &transform, const CTriangle &triangle)
{
	CTriangle t(triangle);
	t.transform(transform);
	return t;
}

CTriangle operator *(const CMatJoint3x4 &transform, const CTriangle &triangle)
{
	CTriangle t(triangle);
	t.transform(transform);
	return t;
}

CTriangle operator *(const CMat4D &transform, const CTriangle &triangle)
{
	CTriangle t(triangle);
	t.transform(transform);
	return t;
}

CTriangle operator *(const CQuaternion &transform, const CTriangle &triangle)
{
	CTriangle t(triangle);
	t.transform(transform);
	return t;
}

std::string CTriangle::toString() const
{
	char str[256];
	std::sprintf(str, "CTriangle(a:(%.2f, %.2f, %.2f) b:(%.2f, %.2f, %.2f) c:(%.2f, %.2f, %.2f))",
		a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z);
	return str;
}

std::ostream &operator <<(std::ostream &o, const CTriangle &triangle)
{
	o << triangle.toString();
	return o;
}





} //end GEO
}  //end SMF

