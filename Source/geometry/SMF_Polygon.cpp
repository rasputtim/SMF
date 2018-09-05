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
#include <float.h>
namespace SMF {
using namespace MATH;
using namespace Util;
namespace GEO{


int CPolygon::numVertices() const
{
	return (int)p.size();
}

int CPolygon::numEdges() const
{
	return (int)p.size();
}

CVec3D CPolygon::vertex(int vertexIndex) const
{
	SMF_ASSERT(vertexIndex >= 0);
	SMF_ASSERT(vertexIndex < (int)p.size());
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (vertexIndex < 0 || vertexIndex >= (int)p.size())
		return CVec3D::nan;
#endif
	return p[vertexIndex];
}

CLineSegment CPolygon::edge(int i) const
{
	if (p.empty())
		return CLineSegment(CVec3D::nan, CVec3D::nan);
	if (p.size() == 1)
		return CLineSegment(p[0], p[0]);
	return CLineSegment(p[i], p[(i+1)%p.size()]);
}

CLineSegment CPolygon::edge2D(int i) const
{
	if (p.empty())
		return CLineSegment(CVec3D::nan, CVec3D::nan);
	if (p.size() == 1)
		return CLineSegment(CVec3D::zero, CVec3D::zero);
	CVec2D temp1=mapTo2D(i);
	CVec2D temp2=mapTo2D((i+1)%p.size());
	return CLineSegment(CVec3D(temp1.x,temp1.y, 0), CVec3D(temp2.x,temp2.y, 0));
}

bool CPolygon::diagonalExists(int i, int j) const
{
	SMF_ASSERT(p.size() >= 3);
	SMF_ASSERT(i >= 0);
	SMF_ASSERT(j >= 0);
	SMF_ASSERT(i < (int)p.size());
	SMF_ASSERT(j < (int)p.size());
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (p.size() < 3 || i < 0 || j < 0 || i >= (int)p.size() || j >= (int)p.size())
		return false;
#endif
	SMF_ASSERT(isPlanar());
	SMF_ASSERT(i != j);
	if (i == j) // Degenerate if i == j.
		return false;
	if (i > j)
		Swap(i, j);
	SMF_ASSERT(i+1 != j);
	if (i+1 == j) // Is this CLineSegment an edge of this polygon?
		return false;

	CPlane polygonPlane = planeCCW();
	CLineSegment diagonal = polygonPlane.project(CLineSegment(p[i], p[j]));

	// First check that this diagonal line is not intersected by an edge of this polygon.
	for(int k = 0; k < (int)p.size(); ++k)
		if (!(k == i || k+1 == i || k == j))
			if (polygonPlane.project(CLineSegment(p[k], p[k+1])).intersects(diagonal))
				return false;

	return isConvex();
}

CVec3D CPolygon::basisU() const
{
	if (p.size() < 2)
		return CVec3D::unitX;
	CVec3D u = p[1] - p[0];
	u.toNormal(); // Always succeeds, even if u was zero (generates (1,0,0)).
	return u;
}

CVec3D CPolygon::basisV() const
{
	if (p.size() < 2)
		return CVec3D::unitY;
	return (planeCCW().getNormal()).cross( basisU()).normalized();
}

CLineSegment CPolygon::diagonal(int i, int j) const
{
	SMF_ASSERT(i >= 0);
	SMF_ASSERT(j >= 0);
	SMF_ASSERT(i < (int)p.size());
	SMF_ASSERT(j < (int)p.size());
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (i < 0 || j < 0 || i >= (int)p.size() || j >= (int)p.size())
		return CLineSegment(CVec3D::nan, CVec3D::nan);
#endif
	return CLineSegment(p[i], p[j]);
}

bool CPolygon::isConvex() const
{
	SMF_ASSERT(isPlanar());
	if (p.empty())
		return false;
	if (p.size() <= 3)
		return true;
	int i = (int)p.size()-2;
	int j = (int)p.size()-1;
	int k = 0;

	while(k < (int)p.size())
	{
		CVec2D a = mapTo2D(i);
		CVec2D b = mapTo2D(j);
		CVec2D c = mapTo2D(k);
		if (!CVec2D::orientedCCW(a, b, c))
			return false;
		i = j;
		j = k;
		++k;
	}
	return true;
}

CVec2D CPolygon::mapTo2D(int i) const
{
	SMF_ASSERT(i >= 0);
	SMF_ASSERT(i < (int)p.size());
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (i < 0 || i >= (int)p.size())
		return CVec2D::nan;
#endif
	return mapTo2D(p[i]);
}

CVec2D CPolygon::mapTo2D(const CVec3D &point) const
{
	SMF_ASSERT(!p.empty());
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (p.empty())
		return CVec2D::nan;
#endif
	CVec3D BasisU = basisU();
	CVec3D BasisV = basisV();
	CVec3D pt = point - p[0];
	return CVec2D((pt* BasisU), (pt* BasisV));
}

CVec3D CPolygon::mapFrom2D(const CVec2D &point) const
{
	SMF_ASSERT(!p.empty());
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (p.empty())
		return CVec3D::nan;
#endif
	return p[0] + point.x * basisU() + point.y * basisV();
}

bool CPolygon::isPlanar(float epsilon) const
{
	if (p.empty())
		return false;
	if (p.size() <= 3)
		return true;
	CPlane plane(p[0], p[1], p[2]);
	for(size_t i = 3; i < p.size(); ++i)
		if (plane.distance(p[i]) > epsilon)
			return false;
	return true;
}

bool CPolygon::isSimple() const
{
	SMF_ASSERT(isPlanar());
	CPlane plane = planeCCW();
	for(int i = 0; i < (int)p.size(); ++i)
	{
		CLineSegment si = plane.project(edge(i));
		for(int j = i+2; j < (int)p.size(); ++j)
		{
			if (i == 0 && j == (int)p.size() - 1)
				continue; // These two edges are consecutive and share a vertex. Don't check that pair.
			CLineSegment sj = plane.project(edge(j));
			if (si.intersects(sj))
				return false;
		}
	}
	return true;
}

bool CPolygon::isNull() const
{
	return p.empty();
}

bool CPolygon::isFinite() const
{
	for(size_t i = 0; i < p.size(); ++i)
		if (!p[i].isFinite())
			return false;

	return true;
}

bool CPolygon::isDegenerate(float epsilon) const
{
	return p.size() < 3 || area() <= epsilon;
}

CVec3D CPolygon::normalCCW() const
{
	///\todo Optimize temporaries.
	return planeCCW().getNormal();
}

CVec3D CPolygon::normalCW() const
{
	///\todo Optimize temporaries.
	return planeCW().getNormal();
}

CPlane CPolygon::planeCCW() const
{
	if (p.size() >= 3)
		return CPlane(p[0], p[1], p[2]);
	if (p.size() == 2)
		return CPlane(CLine(p[0], p[1]), (p[0]-p[1]).perpendicular());
	if (p.size() == 1)
		return CPlane(p[0], CVec3D(0,1,0));
	return CPlane();
}

CPlane CPolygon::planeCW() const
{
	if (p.size() >= 3)
		return CPlane(p[0], p[2], p[1]);
	if (p.size() == 2)
		return CPlane(CLine(p[0], p[1]), (p[0]-p[1]).perpendicular());
	if (p.size() == 1)
		return CPlane(p[0], CVec3D(0,1,0));
	return CPlane();
}

void CPolygon::translate(const CVec3D &offset)
{
	for(size_t i = 0; i < p.size(); ++i)
		p[i] += offset;
}

void CPolygon::transform(const CMat3D &transform)
{
	if (!p.empty())
		transform.batchTransform(&p[0], (int)p.size());
}

void CPolygon::transform(const CMatJoint3x4 &transform)
{
	if (!p.empty())
		transform.batchTransformPos(&p[0], (int)p.size());
}

void CPolygon::transform(const CMat4D &transform)
{
	for(size_t i = 0; i < p.size(); ++i)
		p[i] = transform.MulPos(p[i]);
}

void CPolygon::transform(const CQuaternion &transform)
{
	for(size_t i = 0; i < p.size(); ++i)
		p[i] = transform * p[i];
}

bool CPolygon::contains(const CPolygon &worldSpacePolygon, float polygonThickness) const
{
	for(int i = 0; i < worldSpacePolygon.numVertices(); ++i)
		if (!contains(worldSpacePolygon.vertex(i), polygonThickness))
			return false;
	return true;
}

bool CPolygon::contains(const CVec3D &worldSpacePoint, float polygonThickness) const
{
	// Implementation based on the description from http://erich.realtimerendering.com/ptinpoly/

	if (p.size() < 3)
		return false;

	if (planeCCW().distance(worldSpacePoint) > polygonThickness)
		return false;

	int numIntersections = 0;

	CVec3D BasisU = basisU();
	CVec3D BasisV = basisV();
	SMF_ASSERT(BasisU.isNormalized());
	SMF_ASSERT(BasisV.isNormalized());
	SMF_ASSERT(BasisU.isPerpendicular(basisV()));
	SMF_ASSERT(BasisU.isPerpendicular(planeCCW().getNormal()));
	SMF_ASSERT(BasisV.isPerpendicular(planeCCW().getNormal()));

	CVec2D localSpacePoint = CVec2D((worldSpacePoint* BasisU), (worldSpacePoint* BasisV));

	const float epsilon =CMath::EPSILON_SuperLow;

	CVec2D p0 = CVec2D((p[p.size()-1]* BasisU), (p[p.size()-1]* BasisV)) - localSpacePoint;
	if (CMath::fabs(p0.y) < epsilon)
		p0.y = -epsilon; // Robustness check - if the ray (0,0) -> (+inf, 0) would pass through a vertex, move the vertex slightly.
	for(int i = 0; i < (int)p.size(); ++i)
	{
		CVec2D p1 = CVec2D((p[i]*BasisU), (p[i]* BasisV)) - localSpacePoint;
		if (CMath::fabs(p1.y) < epsilon)
			p0.y = -epsilon; // Robustness check - if the ray (0,0) -> (+inf, 0) would pass through a vertex, move the vertex slightly.

		if (p0.y * p1.y < 0.f)
		{
			if (p0.x > 1e-3f && p1.x > 1e-3f)
				++numIntersections;
			else
			{
				// P = p0 + t*(p1-p0) == (x,0)
				//     p0.x + t*(p1.x-p0.x) == x
				//     p0.y + t*(p1.y-p0.y) == 0
				//                 t == -p0.y / (p1.y - p0.y)

				// Test whether the lines (0,0) -> (+inf,0) and p0 -> p1 intersect at a positive X-coordinate.
				CVec2D d = p1 - p0;
				if (CMath::fabs(d.y) > 1e-5f)
				{
					float t = -p0.y / d.y;
					float x = p0.x + t * d.x;
					if (t >= 0.f && t <= 1.f && x > 1e-3f)
						++numIntersections;
				}
			}
		}
		p0 = p1;
	}

	return numIntersections % 2 == 1;
}

bool CPolygon::contains(const CLineSegment &worldSpaceLineSegment, float polygonThickness) const
{
	if (p.size() < 3)
		return false;

	CPlane plane = planeCCW();
	if (plane.distance(worldSpaceLineSegment.begin) > polygonThickness ||
		plane.distance(worldSpaceLineSegment.end) > polygonThickness)
		return false;

	// For robustness, project onto the polygon plane.
	CLineSegment l = plane.project(worldSpaceLineSegment);

	if (!contains(l.begin) || !contains(l.end))
		return false;

	for(int i = 0; i < (int)p.size(); ++i)
		if (plane.project(edge(i)).intersects(l))
			return false;

	return true;
}

bool CPolygon::contains(const CTriangle &worldSpaceTriangle, float polygonThickness) const
{
	return contains(worldSpaceTriangle.edge(0), polygonThickness) &&
		contains(worldSpaceTriangle.edge(1), polygonThickness) &&
		contains(worldSpaceTriangle.edge(2), polygonThickness);
}

bool CPolygon::contains2D(const CLineSegment &localSpaceLineSegment) const
{
	if (p.size() < 3)
		return false;

///\todo Reimplement this!
//	if (!contains2D(localSpaceLineSegment.a.xy()) || !contains2D(localSpaceLineSegment.b.xy()))
//		return false;

	for(int i = 0; i < (int)p.size(); ++i)
		if (edge2D(i).intersects(localSpaceLineSegment))
			return false;

	return true;
}

bool CPolygon::intersects(const CLine &line) const
{
	float d;
	if (!planeCCW().intersects(line, &d))
		return false;
	return contains(line.getPoint(d));
}

bool CPolygon::intersects(const CRay &ray) const
{
	float d;
	if (!planeCCW().intersects(ray, &d))
		return false;
	return contains(ray.getPoint(d));
}

bool CPolygon::intersects(const CLineSegment &lineSegment) const
{
	CPlane plane = planeCCW();
	float t;
	bool intersects = CPlane::intersectLinePlane(plane.getNormal(), plane.d, lineSegment.begin, lineSegment.end - lineSegment.begin, t);
	if (!intersects || t < 0.f || t > 1.f)
		return false;

	return contains(lineSegment.getPoint(t));
}

bool CPolygon::intersects(const CPlane &plane) const
{
	// project the points of this polygon onto the 1D axis of the plane normal.
	// If there are points on both sides of the plane, then the polygon intersects the plane.
	float minD = CMath::inf;
	float maxD = -CMath::inf;
	for(size_t i = 0; i < p.size(); ++i)
	{
		float d = plane.signedDistance(p[i]);
		minD = MIN(minD, d);
		maxD = MAX(maxD, d);
	}
	// Allow a very small epsilon tolerance.
	return minD <= 1e-4f && maxD >= -1e-4f;
}

bool CPolygon::intersects(const CAABBox &aabb) const
{
	return aabb.intersects(*this);
}

bool CPolygon::intersects(const COBBox &obb) const
{
	return obb.intersects(*this);
}

bool CPolygon::intersects(const CTriangle &triangle) const
{
	return toPolyhedron().intersects(triangle);
}

bool CPolygon::intersects(const CPolygon &polygon) const
{
	return toPolyhedron().intersects(polygon);
}
#if 0
bool CPolygon::intersects(const Frustum &frustum) const
{
	return frustum.intersects(*this);
}
#endif
bool CPolygon::intersects(const CPolyhedron &polyhedron) const
{
	return polyhedron.intersects(*this);
}

bool CPolygon::intersects(const CSphere &sphere) const
{
	///\todo Optimize.
	std::vector<CTriangle> tris = triangulate();
	for(size_t i = 0; i < tris.size(); ++i)
		if (tris[i].intersects(sphere))
			return true;

	return false;
}
#if 0
bool CPolygon::intersects(const Capsule &capsule) const
{
	///\todo Optimize.
	std::vector<CTriangle> tris = triangulate();
	for(size_t i = 0; i < tris.size(); ++i)
		if (tris[i].intersects(capsule))
			return true;

	return false;
}
#endif
CVec3D CPolygon::closestPoint(const CVec3D &point) const
{
	SMF_ASSERT(isPlanar());

	std::vector<CTriangle> tris = triangulate();
	CVec3D closestPt = CVec3D::nan;
	float closestDist = FLT_MAX;
	for(size_t i = 0; i < tris.size(); ++i)
	{
		CVec3D pt = tris[i].closestPoint(point);
		float d = pt.distanceSq(point);
		if (d < closestDist)
		{
			closestPt = pt;
			closestDist = d;
		}
	}
	return closestPt;
}

CVec3D CPolygon::closestPoint(const CLineSegment &lineSegment) const
{
	return closestPoint(lineSegment, 0);
}

CVec3D CPolygon::closestPoint(const CLineSegment &lineSegment, CVec3D *lineSegmentPt) const
{
	std::vector<CTriangle> tris = triangulate();
	CVec3D closestPt = CVec3D::nan;
	CVec3D closestLineSegmentPt = CVec3D::nan;
	float closestDist = FLT_MAX;
	for(size_t i = 0; i < tris.size(); ++i)
	{
		CVec3D lineSegPt;
		CVec3D pt = tris[i].closestPoint(lineSegment, &lineSegPt);
		float d = pt.distanceSq(lineSegPt);
		if (d < closestDist)
		{
			closestPt = pt;
			closestLineSegmentPt = lineSegPt;
			closestDist = d;
		}
	}
	if (lineSegmentPt)
		*lineSegmentPt = closestLineSegmentPt;
	return closestPt;
}

float CPolygon::distance(const CVec3D &point) const
{
	CVec3D pt = closestPoint(point);
	return pt.distance(point);
}

CVec3D CPolygon::edgenormal(int edgeIndex) const
{
	return (edge(edgeIndex).getDir().cross( planeCCW().getNormal()).normalized());
}

CPlane CPolygon::edgePlane(int edgeIndex) const
{
	return CPlane(edge(edgeIndex).begin, edgenormal(edgeIndex));
}

CVec3D CPolygon::extremePoint(const CVec3D &direction) const
{
	CVec3D mostExtreme = CVec3D::nan;
	float mostExtremeDist = -FLT_MAX;
	for(int i = 0; i < numVertices(); ++i)
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

void CPolygon::projectToAxis(const CVec3D &direction, float &outMin, float &outMax) const
{
	///\todo Optimize!
	CVec3D minPt = extremePoint(-direction);
	CVec3D maxPt = extremePoint(direction);
	outMin = (minPt* direction);
	outMax = (maxPt* direction);
}

/*
/// Returns true if the edges of this polygon self-intersect.
bool IsSelfIntersecting() const;

/// Projects all vertices of this polygon to the given plane.
void ProjectToPlane(const CPlane &plane);

/// Returns true if the edges of this polygon self-intersect when viewed from the given direction.
bool IsSelfIntersecting(const CVec3D &viewDirection) const;

bool contains(const CVec3D &point, const CVec3D &viewDirection) const;

*/

/** Implementation based on Graphics Gems 2, p. 170: "IV.1. area of Planar Polygons and volume of Polyhedra." */
float CPolygon::area() const
{
	SMF_ASSERT(isPlanar());
	CVec3D area = CVec3D::zero;
	if (p.size() <= 2)
		return 0.f;

	int i = numEdges()-1;
	for(int j = 0; j < numEdges(); ++j)
	{
		area += vertex(i).cross(vertex(j));
		i = j;
	}
	return 0.5f * CMath::fabs(normalCCW()*(area));
}

float CPolygon::perimeter() const
{
	float perimeter = 0.f;
	for(int i = 0; i < numEdges(); ++i)
		perimeter += edge(i).getLenght();
	return perimeter;
}

///\bug This function does not properly compute the centroid.
CVec3D CPolygon::centroid() const
{
	if (numVertices() == 0)
		return CVec3D::nan;
	CVec3D centroid = CVec3D::zero;
	for(int i = 0; i < numVertices(); ++i)
		centroid += vertex(i);
	return centroid / (float)numVertices();
}

CVec3D CPolygon::pointOnEdge(float normalizedDistance) const
{
	if (p.empty())
		return CVec3D::nan;
	if (p.size() < 2)
		return p[0];
	normalizedDistance = CMath::frac(normalizedDistance); // Take modulo 1 so we have the range [0,1[.
	float Perimeter = perimeter();
	float d = normalizedDistance * Perimeter;
	for(int i = 0; i < numVertices(); ++i)
	{
		CLineSegment edge_ = edge(i);
		float len = edge_.getLenght();
		SMF_ASSERT(len != 0.f && "Degenerate CPolygon detected!");
		if (d <= len)
			return edge_.getPoint(d / len);
		d -= len;
	}
	SMF_ASSERT(false && "CPolygon::pointOnEdge reached end of loop which shouldn't!");
	return p[0];
}

CVec3D CPolygon::randomPointOnEdge(CRandomLCG &rng) const
{
	return pointOnEdge(rng.getFloat());
}

CVec3D CPolygon::fastRandomPointInside(CRandomLCG &rng) const
{
	std::vector<CTriangle> tris = triangulate();
	if (tris.empty())
		return CVec3D::nan;
	return tris[rng.getInt(0, (int)tris.size()-1)].randomPointInside(rng);
}

CPolyhedron CPolygon::toPolyhedron() const
{
	CPolyhedron poly;
	poly.v = p;
	poly.f.push_back(CPolyhedron::Face());
	poly.f.push_back(CPolyhedron::Face());
	for(int i = 0; i < numVertices(); ++i)
	{
		poly.f[0].v.push_back(i);
		poly.f[1].v.push_back(numVertices()-1-i);
	}
	return poly;
}

// A(u) = a1 + u * (a2-a1).
// B(v) = b1 + v * (b2-b1).
// Returns (u,v).
bool IntersectLineLine2D(const CVec2D &a1, const CVec2D &a2, const CVec2D &b1, const CVec2D &b2, CVec2D &out)
{
	float u = (b2.x - b1.x)*(a1.y - b1.y) - (b2.y - b1.y)*(a1.x - b1.x);
	float v = (a2.x - a1.x)*(a1.y - b1.y) - (a2.y - a1.y)*(a1.x - b1.x);

	float det = (b2.y - b1.y)*(a2.x - a1.x) - (b2.x - b1.x)*(a2.y - a1.y);
	if (CMath::fabs(det) < 1e-4f)
		return false;
	det = 1.f / det;
	out.x = u * det;
	out.y = v * det;

	return true;
}

bool IntersectLineSegmentLineSegment2D(const CVec2D &a1, const CVec2D &a2, const CVec2D &b1, const CVec2D &b2, CVec2D &out)
{
	bool ret = IntersectLineLine2D(a1, a2, b1, b2, out);
	return ret && out.x >= 0.f && out.x <= 1.f && out.y >= 0.f && out.y <= 1.f;
}

/// Returns true if poly[i+1] is an ear.
/// Precondition: i+2 == j (mod poly.size()).
bool IsAnEar(const std::vector<CVec2D> &poly, int i, int j)
{
	CVec2D dummy;
	int x = (int)poly.size()-1;
	for(int y = 0; y < i; ++y)
	{
		if (IntersectLineSegmentLineSegment2D(poly[i], poly[j], poly[x], poly[y], dummy))
			return false;
		x = y;
	}
	x = j+1;
	for(int y = x+1; y < (int)poly.size(); ++y)
	{
		if (IntersectLineSegmentLineSegment2D(poly[i], poly[j], poly[x], poly[y], dummy))
			return false;
		x = y;
	}
	return true;
}

/** The implementation of this function is based on the paper
	"Kong, Everett, Toussant. The Graham Scan Triangulates Simple Polygons."
	See also p. 772-775 of Geometric Tools for Computer Graphics.
	The running time of this function is O(n^2). */
std::vector<CTriangle> CPolygon::triangulate() const
{
	SMF_ASSERT(isPlanar());

	std::vector<CTriangle> t;
	// Handle degenerate cases.
	if (numVertices() < 3)
		return t;
	if (numVertices() == 3)
	{
		t.push_back(CTriangle(vertex(0), vertex(1), vertex(2)));
		return t;
	}
	std::vector<CVec2D> p2d;
	std::vector<int> polyIndices;
	for(int v = 0; v < numVertices(); ++v)
	{
		p2d.push_back(mapTo2D(v));
		polyIndices.push_back(v);
	}

	// Clip ears of the polygon until it has been reduced to a triangle.
	int i = 0;
	int j = 1;
	int k = 2;
	size_t numTries = 0; // Avoid creating an infinite loop.
	while(p2d.size() > 3 && numTries < p2d.size())
	{
		if (CVec2D::orientedCCW(p2d[i], p2d[j], p2d[k]) && IsAnEar(p2d, i, k))
		{
			// The vertex j is an ear. Clip it off.
			t.push_back(CTriangle(p[polyIndices[i]], p[polyIndices[j]], p[polyIndices[k]]));
			p2d.erase(p2d.begin() + j);
			polyIndices.erase(polyIndices.begin() + j);

			// The previous index might now have become an ear. Move back one index to see if so.
			if (i > 0)
			{
				i = (i + (int)p2d.size() - 1) % p2d.size();
				j = (j + (int)p2d.size() - 1) % p2d.size();
				k = (k + (int)p2d.size() - 1) % p2d.size();
			}
			numTries = 0;
		}
		else
		{
			// The vertex at j is not an ear. Move to test next vertex.
			i = j;
			j = k;
			k = (k+1) % p2d.size();
			++numTries;
		}
	}

	SMF_ASSERT(p2d.size() == 3);
	if (p2d.size() > 3) // If this occurs, then the polygon is NOT counter-clockwise oriented.
		return t;
/*
	{
		// For conveniency, create a copy that has the winding order fixed, and triangulate that instead.
		// (Causes a large performance hit!)
		CPolygon p2 = *this;
		for(size_t i = 0; i < p2.p.size()/2; ++i)
			std::swap(p2.p[i], p2.p[p2.p.size()-1-i]);
		return p2.triangulate();
	}
*/
	// add the last poly.
	t.push_back(CTriangle(p[polyIndices[0]], p[polyIndices[1]], p[polyIndices[2]]));

	return t;
}

CAABBox CPolygon::minimalEnclosingAABB() const
{
	CAABBox aabb;
	aabb.toNegativeInfinity();
	for(int i = 0; i < numVertices(); ++i)
		aabb.enclose(vertex(i));
	return aabb;
}

/*
/// Returns true if the given vertex is a concave vertex. Otherwise the vertex is a convex vertex.
bool IsConcaveVertex(int i) const;

/// Computes the conves hull of this polygon.
CPolygon convexHull() const;

bool IsSupportingPoint(int i) const;

bool IsSupportingPoint(const CVec3D &point) const;

/// Returns true if the quadrilateral defined by the four points is convex (and not concave or bowtie).
static bool IsConvexQuad(const CVec3D &pointA, const CVec3D &pointB, const CVec3D &pointC, const CVec3D &pointD);
*/

CPolygon operator *(const CMat3D &transform, const CPolygon &polygon)
{
	CPolygon p(polygon);
	p.transform(transform);
	return p;
}

CPolygon operator *(const CMatJoint3x4 &transform, const CPolygon &polygon)
{
	CPolygon p(polygon);
	p.transform(transform);
	return p;
}

CPolygon operator *(const CMat4D &transform, const CPolygon &polygon)
{
	CPolygon p(polygon);
	p.transform(transform);
	return p;
}

CPolygon operator *(const CQuaternion &transform, const CPolygon &polygon)
{
	CPolygon p(polygon);
	p.transform(transform);
	return p;
}


} //end GEO
}  //end SMF


