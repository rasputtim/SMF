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
#include "util/SMF_Debug.h"
#include  <sstream>
#include <set>

namespace SMF {
using namespace MATH;
using namespace Util;
namespace GEO{


#ifdef MATH_GRAPHICSENGINE_INTEROP
#include "VertexBuffer.h"
#endif


void CPolyhedron::Face::flipWindingOrder()
{
	for(size_t i = 0; i < v.size()/2; ++i)
		Swap(v[i], v[v.size()-1-i]);
}

std::string CPolyhedron::Face::toString() const
{
	std::stringstream ss;
	for(size_t i = 0; i < v.size(); ++i)
		ss << v[i] << ((i!=v.size()-1) ? ", " : "");
	return ss.str();
}

int CPolyhedron::numEdges() const
{
	return (int)edgeIndices().size();
}

CVec3D CPolyhedron::vertex(int vertexIndex) const
{
	SMF_ASSERT(vertexIndex >= 0);
	SMF_ASSERT(vertexIndex < (int)v.size());
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (vertexIndex < 0 || vertexIndex >= (int)v.size())
		return CVec3D::nan;
#endif

	return v[vertexIndex];
}

CLineSegment CPolyhedron::edge(int edgeIndex) const
{
	SMF_ASSERT(edgeIndex >= 0);
	std::vector<CLineSegment> edges_ = edges();
	SMF_ASSERT(edgeIndex < (int)edges_.size());
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (edgeIndex < 0 || edgeIndex >= (int)edges_.size())
		return CLineSegment(CVec3D::nan, CVec3D::nan);
#endif
	return edges_[edgeIndex];
}

std::vector<CLineSegment> CPolyhedron::edges() const
{
	std::vector<std::pair<int, int> > edges = edgeIndices();
	std::vector<CLineSegment> edgeLines;
	edgeLines.reserve(edges.size());
	for(size_t i = 0; i < edges.size(); ++i)
		edgeLines.push_back(CLineSegment(vertex(edges[i].first), vertex(edges[i].second)));

	return edgeLines;
}

std::vector<std::pair<int, int> > CPolyhedron::edgeIndices() const
{
	std::set<std::pair<int, int> > uniqueEdges;
	for(int i = 0; i < numFaces(); ++i)
	{
		SMF_ASSERT(f[i].v.size() >= 3);
		if (f[i].v.size() < 3)
			continue; // Degenerate face with less than three vertices, skip!
		int x = f[i].v.back();
		for(size_t j = 0; j < f[i].v.size(); ++j)
		{
			int y = f[i].v[j];
			uniqueEdges.insert(std::make_pair(MIN(x, y), MAX(x, y)));
			x = y;
		}
	}

	std::vector<std::pair<int, int> >edges;
	edges.insert(edges.end(), uniqueEdges.begin(), uniqueEdges.end());
	return edges;
}

std::vector<CPolygon> CPolyhedron::faces() const
{
	std::vector<CPolygon> faces;
	faces.reserve(numFaces());
	for(int i = 0; i < numFaces(); ++i)
		faces.push_back(FacePolygon(i));
	return faces;
}

CPolygon CPolyhedron::FacePolygon(int faceIndex) const
{
	CPolygon p;
	SMF_ASSERT(faceIndex >= 0);
	SMF_ASSERT(faceIndex < (int)f.size());
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (faceIndex < 0 || faceIndex >= (int)f.size())
		return CPolygon();
#endif

	p.p.reserve(f[faceIndex].v.size());
	for(size_t v = 0; v < f[faceIndex].v.size(); ++v)
		p.p.push_back(vertex(f[faceIndex].v[v]));
	return p;
}

CPlane CPolyhedron::facePlane(int faceIndex) const
{
	const Face &face = f[faceIndex];
	if (face.v.size() >= 3)
		return CPlane(v[face.v[0]], v[face.v[1]], v[face.v[2]]);
	else if (face.v.size() == 2)
		return CPlane(CLine(v[face.v[0]], v[face.v[1]]), (v[face.v[0]]-v[face.v[1]]).perpendicular());
	else if (face.v.size() == 1)
		return CPlane(v[face.v[0]], CVec3D(0,1,0));
	else
		return CPlane();
}

CVec3D CPolyhedron::faceNormal(int faceIndex) const
{
	const Face &face = f[faceIndex];
	if (face.v.size() >= 3)
		return (v[face.v[1]]-v[face.v[0]]).cross(v[face.v[2]]-v[face.v[0]]).normalized();
	else if (face.v.size() == 2)
		return (v[face.v[1]]-v[face.v[0]]).cross((v[face.v[0]]-v[face.v[1]]).perpendicular()-v[face.v[0]]).normalized();
	else if (face.v.size() == 1)
		return CVec3D(0,1,0);
	else
		return CVec3D::nan;
}

int CPolyhedron::extremeVertex(const CVec3D &direction) const
{
	int mostExtreme = -1;
	float mostExtremeDist = -FLT_MAX;
	for(int i = 0; i < numVertices(); ++i)
	{
		float d = (direction* vertex(i));
		if (d > mostExtremeDist)
		{
			mostExtremeDist = d;
			mostExtreme = i;
		}
	}
	return mostExtreme;
}

CVec3D CPolyhedron::extremePoint(const CVec3D &direction) const
{
	return vertex(extremeVertex(direction));
}

void CPolyhedron::projectToAxis(const CVec3D &direction, float &outMin, float &outMax) const
{
	///\todo Optimize!
	CVec3D minPt = extremePoint(-direction);
	CVec3D maxPt = extremePoint(direction);
	outMin = (minPt* direction);
	outMax = (maxPt* direction);
}

CVec3D CPolyhedron::centroid() const
{
	CVec3D centroid = CVec3D::zero;
	for(int i = 0; i < numVertices(); ++i)
		centroid += vertex(i);
	return centroid / (float)numVertices();
}

float CPolyhedron::surfaceArea() const
{
	float area = 0.f;
	for(int i = 0; i < numFaces(); ++i)
		area += FacePolygon(i).area(); ///\todo Optimize temporary copies.
	return area;
}

/** The implementation of this function is based on Graphics Gems 2, p. 170: "IV.1. area of Planar Polygons and volume of Polyhedra." */
float CPolyhedron::volume() const
{
	float volume = 0.f;
	for(int i = 0; i < numFaces(); ++i)
	{
		CPolygon face = FacePolygon(i); ///\todo Optimize temporary copies.
		volume += (face.vertex(0)*(face.normalCCW())) * face.area();
	}
	return CMath::fabs(volume) / 3.f;
}

CAABBox CPolyhedron::minimalEnclosingAABB() const
{
	CAABBox aabb;
	aabb.toNegativeInfinity();
	for(int i = 0; i < numVertices(); ++i)
		aabb.enclose(vertex(i));
	return aabb;
}

#ifdef MATH_CONTAINERLIB_SUPPORT
COBBox CPolyhedron::minimalEnclosingOBB() const
{
	return COBBox::OptimalEnclosingOBB(&v[0], (int)v.size());
}
#endif

bool CPolyhedron::faceIndicesValid() const
{
	// Test condition 1: Face indices in proper range.
	for(int i = 0; i < numFaces(); ++i)
		for(int j = 0; j < (int)f[i].v.size(); ++j)
			if (f[i].v[j] < 0 || f[i].v[j] >= (int)v.size())
				return false;

	// Test condition 2: Each face has at least three vertices.
	for(int i = 0; i < numFaces(); ++i)
		if (f[i].v.size() < 3)
			return false;

	// Test condition 3: Each face may refer to a vertex at most once. (Complexity O(n^2)).
	for(int i = 0; i < numFaces(); ++i)
		for(int j = 0; j < (int)f[i].v.size(); ++j)
			for(size_t k = j+1; k < f[i].v.size(); ++k)
				if (f[i].v[j] == f[i].v[k])
					return false;

	return true;
}

void CPolyhedron::flipWindingOrder()
{
	for(size_t i = 0; i < f.size(); ++i)
		f[i].flipWindingOrder();
}

bool CPolyhedron::IsClosed() const
{
	std::set<std::pair<int, int> > uniqueEdges;
	for(int i = 0; i < numFaces(); ++i)
	{
		SMF_ASSERT(FacePolygon(i).isPlanar());
		SMF_ASSERT(FacePolygon(i).isSimple());
		int x = f[i].v.back();
		for(size_t j = 0; j < f[i].v.size(); ++j)
		{
			int y = f[i].v[j];
			if (uniqueEdges.find(std::make_pair(x, y)) != uniqueEdges.end())
			{
			  Debug::debug(Debug::math,__FUNCTION__)<< "CPolyhedron is not simple and closed! The edge ("<< x << " , "<< y <<" ) is used twice!"<<endl;
			  return false; // This edge is being used twice! Cannot be simple and closed.
			}
			uniqueEdges.insert(std::make_pair(x, y));
			x = y;
		}
	}

	for(std::set<std::pair<int, int> >::iterator iter = uniqueEdges.begin();
		iter != uniqueEdges.end(); ++iter)
	{
		std::pair<int, int> reverse = std::make_pair(iter->second, iter->first);
		if (uniqueEdges.find(reverse) == uniqueEdges.end())
		{
			Debug::debug(Debug::math,__FUNCTION__) <<"The edge ( "<< iter->second<<" , " << iter->first<<") does not exist. CPolyhedron is not closed!"<<endl;
			return false;
		}
	}

	return true;
}

bool CPolyhedron::isConvex() const
{
	// This function is O(n^2).
	/** \todo Real-Time Collision Detection, p. 64:
		A faster O(n) approach is to compute for each face F of P the centroid C of F,
		and for all neighboring faces G of F test if C lies behind the supporting plane of
		G. If some C fails to lie behind the supporting plane of one or more neighboring
		faces, P is concave, and is otherwise assumed convex. However, note that just as the
		corresponding polygonal convexity test may fail for a pentagram this test may fail for,
		for example, a pentagram extruded out of its plane and capped at the ends. */

	for(int f = 0; f < numFaces(); ++f)
	{
		CPlane p = facePlane(f);
		for(int i = 0; i < numVertices(); ++i)
		{
			float d = p.signedDistance(vertex(i));
			if (d > 1e-3f) // Tolerate a small epsilon error.
			{
				Debug::debug(Debug::math,__FUNCTION__) <<"distance of vertex "<< i <<" from plane "<< f<<" : "<< d << endl;
				return false;
			}
		}
	}
	return true;
}

bool CPolyhedron::eulerFormulaHolds() const
{
	return numVertices() + numFaces() - numEdges() == 2;
}

bool CPolyhedron::facesAreNondegeneratePlanar(float epsilon) const
{
	for(int i = 0; i < (int)f.size(); ++i)
	{
		const Face &face = f[i];
		if (face.v.size() < 3)
			return false;
		if (face.v.size() >= 4)
		{
			CPlane facePlane_ = facePlane(i);
			for(int j = 0; j < (int)face.v.size(); ++j)
				if (facePlane_.distance(v[face.v[j]]) > epsilon)
					return false;
		}
	}

	return true;
}

bool CPolyhedron::faceContains(int faceIndex, const CVec3D &worldSpacePoint, float polygonThickness) const
{
	// N.B. This implementation is a duplicate of CPolygon::contains, but adapted to avoid dynamic memory allocation
	// related to converting the face of a CPolyhedron to a CPolygon object.

	// Implementation based on the description from http://erich.realtimerendering.com/ptinpoly/

	const Face &face = f[faceIndex];
	const std::vector<int> &vertices = face.v;

	if (vertices.size() < 3)
		return false;

	CPlane p = facePlane(faceIndex);
	if (facePlane(faceIndex).distance(worldSpacePoint) > polygonThickness)
		return false;

	int numIntersections = 0;

	CVec3D basisU = v[vertices[1]] - v[vertices[0]];
	basisU.toNormal();
	CVec3D basisV = (p.getNormal()).cross(basisU).normalized();
	SMF_ASSERT(basisU.isNormalized());
	SMF_ASSERT(basisV.isNormalized());
	SMF_ASSERT(basisU.isPerpendicular(basisV));
	SMF_ASSERT(basisU.isPerpendicular(p.getNormal()));
	SMF_ASSERT(basisV.isPerpendicular(p.getNormal()));

	CVec2D localSpacePoint = CVec2D((worldSpacePoint* basisU), (worldSpacePoint* basisV));

	const float epsilon =CMath::EPSILON_SuperLow;

	CVec2D p0 = CVec2D((v[vertices.back()]* basisU), (v[vertices.back()]* basisV)) - localSpacePoint;
	if (CMath::fabs(p0.y) < epsilon)
		p0.y = -epsilon; // Robustness check - if the ray (0,0) -> (+inf, 0) would pass through a vertex, move the vertex slightly.
	for(size_t i = 0; i < vertices.size(); ++i)
	{
		CVec2D p1 = CVec2D((v[vertices[i]]* basisU), (v[vertices[i]]* basisV)) - localSpacePoint;
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
					if (t >= 0.f && t <= 1.f && x > 1e-6f)
						++numIntersections;
				}
			}
		}
		p0 = p1;
	}

	return numIntersections % 2 == 1;
}

bool CPolyhedron::contains(const CVec3D &point) const
{
	int numIntersections = 0;
	for(int i = 0; i < (int)f.size(); ++i)
	{
		CPlane p(v[f[i].v[0]] - point, v[f[i].v[1]] - point, v[f[i].v[2]] - point);

		// find the intersection of the plane and the ray (0,0,0) -> (t,0,0), t >= 0.
		// <normal, point_on_ray> == d
		// n.x * t == d
		//       t == d / n.x
		if (CMath::fabs(p.getNormal().x) > 1e-5f)
		{
			float t = p.d / p.getNormal().x;
			// If t >= 0, the plane and the ray intersect, and the ray potentially also intersects the polygon.
			// Finish the test by checking whether the point of intersection is contained in the polygon, in
			// which case the ray-polygon intersection occurs.
			if (t >= 0.f && faceContains(i, point + CVec3D(t,0,0)))
				++numIntersections;
		}
	}

	return numIntersections % 2 == 1;
}

bool CPolyhedron::contains(const CLineSegment &lineSegment) const
{
	return contains(lineSegment.begin) && contains(lineSegment.end);
}

bool CPolyhedron::contains(const CTriangle &triangle) const
{
	return contains(triangle.a) && contains(triangle.b) && contains(triangle.c);
}

bool CPolyhedron::contains(const CPolygon &polygon) const
{
	for(int i = 0; i < polygon.numVertices(); ++i)
		if (!contains(polygon.vertex(i)))
			return false;
	return true;
}

bool CPolyhedron::contains(const CAABBox &aabb) const
{
	for(int i = 0; i < 8; ++i)
		if (!contains(aabb.cornerPoint(i)))
			return false;

	return true;
}

bool CPolyhedron::contains(const COBBox &obb) const
{
	for(int i = 0; i < 8; ++i)
		if (!contains(obb.cornerPoint(i)))
			return false;

	return true;
}


bool CPolyhedron::contains(const CPolyhedron &polyhedron) const
{
	SMF_ASSERT(polyhedron.IsClosed());
	for(int i = 0; i < polyhedron.numVertices(); ++i)
		if (!contains(polyhedron.vertex(i)))
			return false;

	return true;
}

bool CPolyhedron::containsConvex(const CVec3D &point) const
{
	SMF_ASSERT(isConvex());
	for(int i = 0; i < numFaces(); ++i)
		if (facePlane(i).signedDistance(point) > 0.f)
			return false;

	return true;
}

bool CPolyhedron::containsConvex(const CLineSegment &lineSegment) const
{
	return containsConvex(lineSegment.begin) && containsConvex(lineSegment.end);
}

bool CPolyhedron::containsConvex(const CTriangle &triangle) const
{
	return containsConvex(triangle.a) && containsConvex(triangle.b) && containsConvex(triangle.c);
}

CVec3D CPolyhedron::closestPointConvex(const CVec3D &point) const
{
	SMF_ASSERT(isConvex());
	if (containsConvex(point))
		return point;
	CVec3D closestPoint = CVec3D::nan;
	float closestDistance = FLT_MAX;
	for(int i = 0; i < numFaces(); ++i)
	{
		CVec3D closestOnPoly = FacePolygon(i).closestPoint(point);
		float d = closestOnPoly.distanceSq(point);
		if (d < closestDistance)
		{
			closestPoint = closestOnPoly;
			closestDistance = d;
		}
	}
	return closestPoint;
}

CVec3D CPolyhedron::closestPoint(const CVec3D &point) const
{
	if (contains(point))
		return point;
	CVec3D closestPoint = CVec3D::nan;
	float closestDistance = FLT_MAX;
	for(int i = 0; i < numFaces(); ++i)
	{
		CVec3D closestOnPoly = FacePolygon(i).closestPoint(point);
		float d = closestOnPoly.distanceSq(point);
		if (d < closestDistance)
		{
			closestPoint = closestOnPoly;
			closestDistance = d;
		}
	}
	return closestPoint;
}

CVec3D CPolyhedron::closestPoint(const CLineSegment &lineSegment) const
{
	return closestPoint(lineSegment, 0);
}

CVec3D CPolyhedron::closestPoint(const CLineSegment &lineSegment, CVec3D *lineSegmentPt) const
{
	if (contains(lineSegment.begin))
	{
		if (lineSegmentPt)
			*lineSegmentPt = lineSegment.begin;
		return lineSegment.begin;
	}
	if (contains(lineSegment.end))
	{
		if (lineSegmentPt)
			*lineSegmentPt = lineSegment.end;
		return lineSegment.end;
	}
	CVec3D closestPt = CVec3D::nan;
	float closestDistance = FLT_MAX;
	CVec3D closestLineSegmentPt = CVec3D::nan;
	for(int i = 0; i < numFaces(); ++i)
	{
		CVec3D lineSegPt;
		CVec3D pt = FacePolygon(i).closestPoint(lineSegment, &lineSegPt);
		float d = pt.distanceSq(lineSegPt);
		if (d < closestDistance)
		{
			closestDistance = d;
			closestPt = pt;
			closestLineSegmentPt = lineSegPt;
		}
	}
	if (lineSegmentPt)
		*lineSegmentPt = closestLineSegmentPt;
	return closestPt;
}

float CPolyhedron::distance(const CVec3D &point) const
{
	CVec3D pt = closestPoint(point);
	return pt.distance(point);
}

bool CPolyhedron::clipLineSegmentToConvexPolyhedron(const CVec3D &ptA, const CVec3D &dir,
                                                   float &tFirst, float &tLast) const
{
	SMF_ASSERT(isConvex());

	// intersects line segment against each plane.
	for(int i = 0; i < numFaces(); ++i)
	{
		/* Denoting the dot product of vectors a and b with <a,b>, we have:

		   The points P on the plane p satisfy the equation <P, p.normal> == p.d.
		   The points P on the line have the parametric equation P = ptA + dir * t.
		   Solving for the distance along the line for intersection gives

		   t = (p.d - <p.normal, ptA>) / <p.normal, dir>.
		*/

		CPlane p = facePlane(i);
		float denom = (p.getNormal()* dir);
		float dist = p.d - (p.getNormal()* ptA);

		// Avoid division by zero. In this case the line segment runs parallel to the plane.
		if (CMath::fabs(denom) < 1e-5f)
		{
			// If <P, p.normal> < p.d, then the point lies in the negative halfspace of the plane, which is inside the polyhedron.
			// If <P, p.normal> > p.d, then the point lies in the positive halfspace of the plane, which is outside the polyhedron.
			// Therefore, if p.d - <ptA, p.normal> == dist < 0, then the whole line is outside the polyhedron.
			if (dist < 0.f)
				return false;
		}
		else
		{
			float t = dist / denom;
			if (denom < 0.f) // When entering halfspace, update tFirst if t is larger.
				tFirst = MAX(t, tFirst);
			else // When exiting halfspace, updeate tLast if t is smaller.
				tLast = MIN(t, tLast);

			if (tFirst > tLast)
				return false; // We clipped the whole line segment.
		}
	}
	return true;
}

bool CPolyhedron::intersects(const CLineSegment &lineSegment) const
{
	if (contains(lineSegment))
		return true;
	for(int i = 0; i < numFaces(); ++i)
	{
		float t;
		CPlane plane = facePlane(i);
		bool intersects = CPlane::intersectLinePlane(plane.getNormal(), plane.d, lineSegment.begin, lineSegment.end - lineSegment.begin, t);
		if (intersects && t >= 0.f && t <= 1.f)
			if (faceContains(i, lineSegment.getPoint(t)))
				return true;
	}

	return false;
}

bool CPolyhedron::intersects(const CLine &line) const
{
	for(int i = 0; i < numFaces(); ++i)
		if (FacePolygon(i).intersects(line))
			return true;

	return false;
}

bool CPolyhedron::intersects(const CRay &ray) const
{
	for(int i = 0; i < numFaces(); ++i)
		if (FacePolygon(i).intersects(ray))
			return true;

	return false;
}

bool CPolyhedron::intersects(const CPlane &plane) const
{
	return plane.intersects(*this);
}

/** The algorithm for CPolyhedron-CPolyhedron intersection is from Christer Ericson's Real-Time Collision Detection, p. 384.
	As noted by the author, the algorithm is very naive (and here unoptimized), and better methods exist. [groupSyntax] */
bool CPolyhedron::intersects(const CPolyhedron &polyhedron) const
{
	if (polyhedron.contains(this->centroid()))
		return true;
	if (this->contains(polyhedron.centroid()))
		return true;

	// This test assumes that both this and the other polyhedron are closed.
	// This means that for each edge running through vertices i and j, there's a face
	// that contains the line segment (i,j) and another neighboring face that contains
	// the line segment (j,i). These represent the same line segment (but in opposite direction)
	// so we only have to test one of them for intersection. Take i < j as the canonical choice
	// and skip the other winding order.

	// Test for each edge of this polyhedron whether the other polyhedron intersects it.
	for(size_t i = 0; i < f.size(); ++i)
	{
		SMF_ASSERT(!f[i].v.empty()); // Cannot have degenerate faces here, and for performance reasons, don't start checking for this condition in release mode!
		int v0 = f[i].v.back();
		CVec3D l0 = v[v0];
		for(size_t j = 0; j < f[i].v.size(); ++j)
		{
			int v1 = f[i].v[j];
			CVec3D l1 = v[v1];
			if (v0 < v1 && polyhedron.intersects(CLineSegment(l0, l1))) // If v0 < v1, then this line segment is the canonical one.
				return true;
			l0 = l1;
			v0 = v1;
		}
	}

	// Test for each edge of the other polyhedron whether this polyhedron intersects it.
	for(size_t i = 0; i < polyhedron.f.size(); ++i)
	{
		SMF_ASSERT(!polyhedron.f[i].v.empty()); // Cannot have degenerate faces here, and for performance reasons, don't start checking for this condition in release mode!
		int v0 = polyhedron.f[i].v.back();
		CVec3D l0 = polyhedron.v[v0];
		for(size_t j = 0; j < polyhedron.f[i].v.size(); ++j)
		{
			int v1 = polyhedron.f[i].v[j];
			CVec3D l1 = polyhedron.v[v1];
			if (v0 < v1 && intersects(CLineSegment(l0, l1))) // If v0 < v1, then this line segment is the canonical one.
				return true;
			l0 = l1;
			v0 = v1;
		}
	}

	return false;
}

template<typename T>
bool PolyhedronIntersectsAABB_OBB(const CPolyhedron &p, const T &obj)
{
	if (p.contains(obj.getCenter()))
		return true;
	if (obj.contains(p.centroid()))
		return true;

	// Test for each edge of the CAABBox/COBBox whether this polyhedron intersects it.
	for(int i = 0; i < obj.numEdges(); ++i)
		if (p.intersects(obj.edge(i)))
			return true;

	// Test for each edge of this polyhedron whether the CAABBox/COBBox intersects it.
	for(size_t i = 0; i < p.f.size(); ++i)
	{
		SMF_ASSERT(!p.f[i].v.empty()); // Cannot have degenerate faces here, and for performance reasons, don't start checking for this condition in release mode!
		int v0 = p.f[i].v.back();
		CVec3D l0 = p.v[v0];
		for(size_t j = 0; j < p.f[i].v.size(); ++j)
		{
			int v1 = p.f[i].v[j];
			CVec3D l1 = p.v[v1];
			if (v0 < v1 && obj.intersects(CLineSegment(l0, l1))) // If v0 < v1, then this line segment is the canonical one.
				return true;
			l0 = l1;
			v0 = v1;
		}
	}

	return false;
}

bool CPolyhedron::intersects(const CAABBox &aabb) const
{
	return PolyhedronIntersectsAABB_OBB(*this, aabb);
}

bool CPolyhedron::intersects(const COBBox &obb) const
{
	return PolyhedronIntersectsAABB_OBB(*this, obb);
}

bool CPolyhedron::intersects(const CTriangle &triangle) const
{
	return PolyhedronIntersectsAABB_OBB(*this, triangle);
}

bool CPolyhedron::intersects(const CPolygon &polygon) const
{
	return intersects(polygon.toPolyhedron());
}



bool CPolyhedron::intersects(const CSphere &sphere) const
{
	CVec3D closestPt = closestPoint(sphere.getOrigin());
	return closestPt.distanceSq(sphere.getOrigin()) <= sphere.getRadius() * sphere.getRadius();
}


bool CPolyhedron::intersectsConvex(const CLine &line) const
{
	float tFirst = -FLT_MAX;
	float tLast = FLT_MAX;
	return clipLineSegmentToConvexPolyhedron(line.pos, line.dir, tFirst, tLast);
}

bool CPolyhedron::intersectsConvex(const CRay &ray) const
{
	float tFirst = 0.f;
	float tLast = FLT_MAX;
	return clipLineSegmentToConvexPolyhedron(ray.pos, ray.dir, tFirst, tLast);
}

bool CPolyhedron::intersectsConvex(const CLineSegment &lineSegment) const
{
	float tFirst = 0.f;
	float tLast = 1.f;
	return clipLineSegmentToConvexPolyhedron(lineSegment.begin, lineSegment.end - lineSegment.begin, tFirst, tLast);
}

void CPolyhedron::mergeConvex(const CVec3D &point)
{
//	LOGI("mergeconvex.");
	std::set<std::pair<int, int> > deletedEdges;
	std::map<std::pair<int, int>, int> remainingEdges;

	for(size_t i = 0; i < v.size(); ++i)
		if (point.distanceSq(v[i]) < 1e-3f)
			return;

//	bool hadDisconnectedHorizon = false;

	for(int i = 0; i < (int)f.size(); ++i)
	{
		// Delete all faces that don't contain the given point. (they have point in their positive side)
		CPlane p = facePlane(i);
		Face &face = f[i];
		if (p.signedDistance(point) > 1e-5f)
		{
			bool isConnected = (deletedEdges.empty());

			int v0 = face.v.back();
			for(size_t j = 0; j < face.v.size() && !isConnected; ++j)
			{
				int v1 = face.v[j];
				if (deletedEdges.find(std::make_pair(v1, v0)) != deletedEdges.end())
				{
					isConnected = true;
					break;
				}
				v0 = v1;
			}

			if (isConnected)
			{
				v0 = face.v.back();
				for(size_t j = 0; j < face.v.size(); ++j)
				{
					int v1 = face.v[j];
					deletedEdges.insert(std::make_pair(v0, v1));
			//		LOGI("edge %d,%d is to be deleted.", v0, v1);
					v0 = v1;
				}
		//		LOGI("Deleting face %d: %s. distance to vertex %f", i, face.toString().c_str(), p.signedDistance(point));
				std::swap(f[i], f.back());
				f.pop_back();
				--i;
				continue;
			}
//			else
//				hadDisconnectedHorizon = true;
		}

		int v0 = face.v.back();
		for(size_t j = 0; j < face.v.size(); ++j)
		{
			int v1 = face.v[j];
			remainingEdges[std::make_pair(v0, v1)] = i;
	//		LOGI("edge %d,%d is to be deleted.", v0, v1);
			v0 = v1;
		}

	}

	// The polyhedron contained our point, nothing to merge.
	if (deletedEdges.empty())
		return;

	// add the new point to this polyhedron.
//	if (!v.back().compare(point))
		v.push_back(point);

/*
	// Create a look-up index of all remaining uncapped edges of the polyhedron.
	std::map<std::pair<int,int>, int> edgesToFaces;
	for(size_t i = 0; i < f.size(); ++i)
	{
		Face &face = f[i];
		int v0 = face.v.back();
		for(size_t j = 0; j < face.v.size(); ++j)
		{
			int v1 = face.v[j];
			edgesToFaces[std::make_pair(v1, v0)] = i;
			v0 = v1;
		}
	}
*/
	// Now fix all edges by adding new triangular faces for the point.
//	for(size_t i = 0; i < deletedEdges.size(); ++i)
	for(std::set<std::pair<int, int> >::iterator iter = deletedEdges.begin(); iter != deletedEdges.end(); ++iter)
	{
		std::pair<int, int> opposite = std::make_pair(iter->second, iter->first);
		if (deletedEdges.find(opposite) != deletedEdges.end())
			continue;

//		std::map<std::pair<int,int>, int>::iterator iter = edgesToFaces.find(deletedEdges[i]);
//		std::map<std::pair<int,int>, int>::iterator iter = edgesToFaces.find(deletedEdges[i]);
//		if (iter != edgesToFaces.end())
		{
			// If the adjoining face is planar to the triangle we'd like to add, instead extend the face to enclose
			// this vertex.
			//CVec3D newTriangleNormal = (v[v.size()-1]-v[iter->second]).cross(v[iter->first]-v[iter->second]).normalized();

			std::map<std::pair<int, int>, int>::iterator existing = remainingEdges.find(opposite);
			SMF_ASSERT(existing != remainingEdges.end());
			MARK_UNUSED(existing);

#if 0
			int adjoiningFace = existing->second;

			if (faceNormal(adjoiningFace).dot(newTriangleNormal) >= 0.99999f) ///\todo CVec3D::IsCollinear
			{
				bool added = false;
				Face &adjoining = f[adjoiningFace];
				for(size_t i = 0; i < adjoining.v.size(); ++i)
					if (adjoining.v[i] == iter->second)
					{
						adjoining.v.insert(adjoining.v.begin() + i + 1, v.size()-1);
						added = true;
						/*
						int prev2 = (i + adjoining.v.size() - 1) % adjoining.v.size();
						int prev = i;
						int cur = i + 1;
						int next = (i + 2) % adjoining.v.size();
						int next2 = (i + 3) % adjoining.v.size();

						if (CVec3D::areCollinear(v[prev2], v[prev], v[cur]))
							adjoining.v.erase(adjoining.v.begin() + prev);
						else if (CVec3D::areCollinear(v[prev], v[cur], v[next]))
							adjoining.v.erase(adjoining.v.begin() + cur);
						else if (CVec3D::areCollinear(v[cur], v[next], v[next2]))
							adjoining.v.erase(adjoining.v.begin() + next2);
							*/

						break;
					}
				SMF_ASSERT(added);
				SMF_ASSERT(added);
			}
			else
#endif
//			if (!v[deletedEdges[i].first].compare(point) && !v[deletedEdges[i].second].compare(point))
			{
				Face tri;
				tri.v.push_back(iter->second);
				tri.v.push_back((int)v.size()-1);
				tri.v.push_back(iter->first);
				f.push_back(tri);
	//			LOGI("Added face %d: %s.", (int)f.size()-1, tri.toString().c_str());
			}
		}
	}

#define SMF_ASSERTeq(lhs, op, rhs) do { if (!((lhs) op (rhs))) { LOGE("Condition %s %s %s (%g %s %g) failed!", #lhs, #op, #rhs, (double)(lhs), #op, (double)(rhs)); SMF_ASSERT(false); } } while(0)

//	SMF_ASSERTeq(numVertices() + numFaces(), ==, 2 + numEdges());
	SMF_ASSERT(faceIndicesValid());
//	SMF_ASSERT(eulerFormulaHolds());
//	SMF_ASSERT(IsClosed());
//	SMF_ASSERT(facesAreNondegeneratePlanar());
//	SMF_ASSERT(isConvex());

//	if (hadDisconnectedHorizon)
//		mergeConvex(point);
}

void CPolyhedron::translate(const CVec3D &offset)
{
	for(size_t i = 0; i < v.size(); ++i)
		v[i] += offset;
}

void CPolyhedron::transform(const CMat3D &transform)
{
	if (!v.empty())
		transform.batchTransform(&v[0], (int)v.size());
}

void CPolyhedron::transform(const CMatJoint3x4 &transform)
{
	if (!v.empty())
		transform.batchTransformPos(&v[0], (int)v.size());
}

void CPolyhedron::transform(const CMat4D &transform)
{
	for(size_t i = 0; i < v.size(); ++i)
		v[i] = transform.MulPos(v[i]); ///\todo add CMat4D::batchTransformPos.
}

void CPolyhedron::transform(const CQuaternion &transform)
{
	for(size_t i = 0; i < v.size(); ++i)
		v[i] = transform * v[i];
}

void CPolyhedron::orientNormalsOutsideConvex()
{
	CVec3D center = v[0];
	for(size_t i = 1; i < v.size(); ++i)
		center += v[i];

	center /= (float)v.size();
	for(int i = 0; i < (int)f.size(); ++i)
		if (facePlane(i).signedDistance(center) > 0.f)
			f[i].flipWindingOrder();
}

/// edge from v1->v2.
struct AdjEdge
{
//	int v1;
//	int v2;
	int f1; // The face that has v1->v2.
	int f2; // The face that has v2->v1.
};

#include <list>

struct CHullHelp
{
	std::map<std::pair<int,int>, AdjEdge> edges;
	std::list<int> livePlanes;
};
#ifndef ARRAY_LENGTH
#define ARRAY_LENGTH(x) (sizeof((x))/sizeof((x)[0]))
#endif

CPolyhedron CPolyhedron::convexHull(const CVec3D *pointArray, int numPoints)
{
	///\todo Check input ptr and size!
	std::set<int> extremes;

	const CVec3D dirs[] =
	{
		CVec3D(1,0,0), CVec3D(0,1,0), CVec3D(0,0,1),
		CVec3D(1,1,0), CVec3D(1,0,1), CVec3D(0,1,1),
		CVec3D(1,1,1)
	};

	for(size_t i = 0; i < ARRAY_LENGTH(dirs); ++i)
	{
		int idx1, idx2;
		COBBox::extremePointsAlongDirection(dirs[i], pointArray, numPoints, idx1, idx2);
		extremes.insert(idx1);
		extremes.insert(idx2);
	}

	CPolyhedron p;
	//s SMF_ASSERT(extremes.size() >= 4); ///\todo Fix this case!
	int i = 0;
	std::set<int>::iterator iter = extremes.begin();
	for(; iter != extremes.end() && i < 4; ++iter, ++i)
		p.v.push_back(pointArray[*iter]);

	Face f;
	f.v.resize(3);
	f.v[0] = 0; f.v[1] = 1; f.v[2] = 2; p.f.push_back(f);
	f.v[0] = 0; f.v[1] = 1; f.v[2] = 3; p.f.push_back(f);
	f.v[0] = 0; f.v[1] = 2; f.v[2] = 3; p.f.push_back(f);
	f.v[0] = 1; f.v[1] = 2; f.v[2] = 3; p.f.push_back(f);
	p.orientNormalsOutsideConvex(); // Ensure that the winding order of the generated tetrahedron is correct for each face.

//	SMF_ASSERT(p.IsClosed());
	//SMF_ASSERT(p.isConvex());
	SMF_ASSERT(p.faceIndicesValid());
	SMF_ASSERT(p.eulerFormulaHolds());
//	SMF_ASSERT(p.facesAreNondegeneratePlanar());

	CHullHelp hull;
	for(int j = 0; j < (int)p.f.size(); ++j)
		hull.livePlanes.push_back(j);

	// For better performance, merge the remaining extreme points first.
	for(; iter != extremes.end(); ++iter)
	{
		p.mergeConvex(pointArray[*iter]);

//s		SMF_ASSERT(p.faceIndicesValid());
//		SMF_ASSERT(p.IsClosed());
//		SMF_ASSERT(p.facesAreNondegeneratePlanar());
//		SMF_ASSERT(p.isConvex());
	}

	// Merge all the rest of the points.
	for(int j = 0; j < numPoints; ++j)
	{
		if (p.f.size() > 5000 && (j & 255) == 0)
			Debug::debug(Debug::math,__FUNCTION__)<< "Mergeconvex: "<< j << "/"<< numPoints <<", #vertices: "<< (int)p.v.size()<<", #faces: "<<(int)p.f.size()<<endl;

		p.mergeConvex(pointArray[i]);

//s		SMF_ASSERT(p.faceIndicesValid());
//		SMF_ASSERT(p.IsClosed());
//		SMF_ASSERT(p.facesAreNondegeneratePlanar());
		//SMF_ASSERT(p.isConvex());

//		if (p.f.size() > 5000)
//			break;
	}

	return p;
}

/// See http://paulbourke.net/geometry/platonic/
CPolyhedron CPolyhedron::tetraHedron(const CVec3D &centerPos, float scale, bool ccwIsFrontFacing)
{
	const CVec3D vertices[4] = { CVec3D(1,1,1),
	                             CVec3D(-1,1,-1),
	                             CVec3D(1,-1,-1),
	                             CVec3D(-1,-1,1) };
	const int faces[4][3] = { { 0, 1, 2 },
	                          { 1, 3, 2 },
	                          { 0, 2, 3 },
	                          { 0, 3, 1 } };

	scale /= 2.f;
	CPolyhedron p;

	for(int i = 0; i < 4; ++i)
		p.v.push_back(vertices[i]*scale + centerPos);

	for(int i = 0; i < 4; ++i)
	{
		Face f;
		for(int j = 0; j < 3; ++j)
			f.v.push_back(faces[i][j]);
		p.f.push_back(f);
	}

	if (!ccwIsFrontFacing)
		p.flipWindingOrder();

	return p;
}

/// See http://paulbourke.net/geometry/platonic/
CPolyhedron CPolyhedron::octaHedron(const CVec3D &centerPos, float scale, bool ccwIsFrontFacing)
{
	float a = 1.f / (2.f * CMath::sqrt(2.f));
	float b = 0.5f;

	const CVec3D vertices[6] = { CVec3D(-a, 0, a),
	                             CVec3D(-a, 0,-a),
	                             CVec3D( 0, b, 0),
	                             CVec3D( a, 0,-a),
	                             CVec3D( 0,-b, 0),
	                             CVec3D( a, 0, a) };
	const int faces[8][3] = { { 0, 1, 2 },
	                          { 1, 3, 2 },
	                          { 3, 5, 2 },
	                          { 5, 0, 2 },
	                          { 3, 1, 4 },
	                          { 1, 0, 4 },
	                          { 5, 3, 4 },
	                          { 0, 5, 4 } };

	scale /= 2.f;
	CPolyhedron p;

	for(int i = 0; i < 6; ++i)
		p.v.push_back(vertices[i]*scale + centerPos);

	for(int i = 0; i < 8; ++i)
	{
		Face f;
		for(int j = 0; j < 3; ++j)
			f.v.push_back(faces[i][j]);
		p.f.push_back(f);
	}

	if (!ccwIsFrontFacing)
		p.flipWindingOrder();

	return p;
}

/// See http://paulbourke.net/geometry/platonic/
CPolyhedron CPolyhedron::hexaHedron(const CVec3D &centerPos, float scale, bool ccwIsFrontFacing)
{
	CAABBox aabb(CVec3D(-1,-1,-1), CVec3D(1,1,1));
	aabb.scale(CVec3D::zero, scale * 0.5f);
	aabb.translate(centerPos);
	CPolyhedron p = aabb.toPolyhedron();
	if (ccwIsFrontFacing)
		p.flipWindingOrder();
	return p;
}

/// See http://paulbourke.net/geometry/platonic/
CPolyhedron CPolyhedron::icosaHedron(const CVec3D &centerPos, float scale, bool ccwIsFrontFacing)
{
	float a = 0.5f;
	float phi = (1.f + CMath::sqrt(5.f)) / 2.f;
	float b = 1.f / (2.f * phi);

	const CVec3D vertices[12] = { CVec3D( 0,  b, -a),
	                              CVec3D( b,  a,  0),
	                              CVec3D(-b,  a,  0),
	                              CVec3D( 0,  b,  a),
	                              CVec3D( 0, -b,  a),
	                              CVec3D(-a,  0,  b),
	                              CVec3D( a,  0,  b),
	                              CVec3D( 0, -b, -a),
	                              CVec3D(-a,  0, -b),
	                              CVec3D(-b, -a,  0),
	                              CVec3D( b, -a,  0),
	                              CVec3D( a,  0, -b) };
	const int faces[20][3] = { { 0,  1,  2 },
	                           { 3,  2,  1 },
	                           { 3,  4,  5 },
	                           { 3,  6,  4 },
	                           { 0,  7, 11 },
	                           { 0,  8,  7 },
	                           { 4, 10,  9 },
	                           { 7,  9, 10 },
	                           { 2,  5,  8 },
	                           { 9,  8,  5 },
	                           { 1, 11,  6 },
	                           { 10, 6, 11 },
	                           { 3,  5,  2 },
	                           { 3,  1,  6 },
	                           { 0,  2,  8 },
	                           { 0, 11,  1 },
	                           { 7,  8,  9 },
	                           { 7, 10, 11 },
	                           { 4,  9,  5 },
	                           { 4,  6, 10 } };

	CPolyhedron p;

	for(int i = 0; i < 12; ++i)
		p.v.push_back(vertices[i]*scale + centerPos);

	for(int i = 0; i < 20; ++i)
	{
		Face f;
		for(int j = 0; j < 3; ++j)
			f.v.push_back(faces[i][j]);
		p.f.push_back(f);
	}

	if (!ccwIsFrontFacing)
		p.flipWindingOrder();

	return p;
}

/// See http://paulbourke.net/geometry/platonic/
CPolyhedron CPolyhedron::dodecaHedron(const CVec3D &centerPos, float scale, bool ccwIsFrontFacing)
{
	float phi = (1.f + CMath::sqrt(5.f)) / 2.f;
	float b = 1.f / phi;
	float c = 2.f - phi;

	const CVec3D vertices[20] = { CVec3D( c,  0,  1),
	                              CVec3D(-c,  0,  1),
	                              CVec3D(-b,  b,  b),
	                              CVec3D( 0,  1,  c),
	                              CVec3D( b,  b,  b),
	                              CVec3D( b, -b,  b),
	                              CVec3D( 0, -1,  c),
	                              CVec3D(-b, -b,  b),
	                              CVec3D( 0, -1, -c),
	                              CVec3D( b, -b, -b),
	                              CVec3D(-c,  0, -1),
	                              CVec3D( c,  0, -1),
	                              CVec3D(-b, -b, -b),
	                              CVec3D( b,  b, -b),
	                              CVec3D( 0,  1, -c),
	                              CVec3D(-b,  b, -b),
	                              CVec3D( 1,  c,  0),
	                              CVec3D(-1,  c,  0),
	                              CVec3D(-1, -c,  0),
	                              CVec3D( 1, -c,  0) };

	const int faces[12][5] = { {  0,  1,  2,  3,  4 },
	                           {  1,  0,  5,  6,  7 },
	                           { 11, 10, 12,  8,  9 },
	                           { 10, 11, 13, 14, 15 },
	                           { 13, 16,  4,  3, 14 }, // Note: The winding order of this face was flipped from PBourke's original representation.
	                           {  2, 17, 15, 14,  3 }, //       Winding order flipped.
	                           { 12, 18,  7,  6,  8 }, //       Winding order flipped.
	                           {  5, 19,  9,  8,  6 }, //       Winding order flipped.
	                           { 16, 19,  5,  0,  4 },
	                           { 19, 16, 13, 11,  9 },
	                           { 17, 18, 12, 10, 15 },
	                           { 18, 17,  2,  1,  7 } };

	scale /= 2.f;
	CPolyhedron p;

	for(int i = 0; i < 20; ++i)
		p.v.push_back(vertices[i]*scale + centerPos);

	for(int i = 0; i < 12; ++i)
	{
		Face f;
		for(int j = 0; j < 5; ++j)
			f.v.push_back(faces[i][j]);
		p.f.push_back(f);
	}

	if (!ccwIsFrontFacing)
		p.flipWindingOrder();

	return p;
}

int IntTriCmp(int a, int b)
{
	return a - b;
}

/** Does a binary search on the array list that is sorted in ascending order.
	\param list [in] A pointer to the array to search.
	\param numItems The number of elements in the array 'list'.
	\param value The element to search for.
	\param cmp The comparison operator to use. The comparison function is of form
		int CmpFunc(const T &a, const T &b), and it tests the mutual order of a and b.
		It should return -1 if a < b, +1 if a > b, and 0 if a == b.
	\return The index where a matching element lies, or -1 if not found. Note that if there are more than
	        one matching element, the first that is found is returned. */
template<typename T, typename CmpFunc>
int ArrayBinarySearch(const T *list, int numItems, const T &value, CmpFunc &cmp)
{
	int left = 0;
	int right = numItems-1;
	int order = cmp(list[left], value);
	if (order > 0) return -1;
	if (order == 0) return left;

	order = cmp(list[right], value);
	if (order < 0) return -1;
	if (order == 0) return right;

	int round = 0; // Counter to alternatingly round up or down.
	do
	{
		int middle = (left + right + round) / 2;
		round = (round+1)&1;
		order = cmp(list[middle], value);
		if (order == 0)
			return middle;
		if (order < 0)
			left = middle;
		else right = middle;
	} while(left < right);
	return -1;
}

void CPolyhedron::removeRedundantVertices()
{
	std::set<int> usedVertices;

	// Gather all used vertices.
	for(size_t i = 0; i < f.size(); ++i)
		for(size_t j = 0; j < f[i].v.size(); ++j)
			usedVertices.insert(f[i].v[j]);

	// Turn the used vertices set into a vector for random access.
	std::vector<int> usedVerticesArray;
	usedVerticesArray.reserve(usedVertices.size());
	for(std::set<int>::iterator iter = usedVertices.begin(); iter != usedVertices.end(); ++iter)
		usedVerticesArray.push_back(*iter);

	// Shift all face indices to point to the new vertex array.
	for(size_t i = 0; i < f.size(); ++i)
		for(size_t j = 0; j < f[i].v.size(); ++j)
		{
			int oldIndex = f[i].v[j];
			int newIndex = ArrayBinarySearch(&usedVerticesArray[0], (int)usedVerticesArray.size(), oldIndex, IntTriCmp);
			SMF_ASSERT(newIndex != -1);
			f[i].v[j] = newIndex;
		}

	// Delete all unused vertices from the vertex array.
	for(size_t i = 0; i < usedVerticesArray.size(); ++i)
		v[i] = v[usedVerticesArray[i]];
	v.resize(usedVerticesArray.size());

	SMF_ASSERT(faceIndicesValid());
}

void CPolyhedron::mergeAdjacentPlanarFaces()
{

}

std::vector<CTriangle> CPolyhedron::triangulate() const
{
	std::vector<CTriangle> outTriangleList;
	for(int i = 0; i < numFaces(); ++i)
	{
		CPolygon p = FacePolygon(i);
		std::vector<CTriangle> tris = p.triangulate();
		outTriangleList.insert(outTriangleList.end(), tris.begin(), tris.end());
	}
	return outTriangleList;
}

#ifdef MATH_GRAPHICSENGINE_INTEROP
void CPolyhedron::triangulate(VertexBuffer &vb, bool ccwIsFrontFacing) const
{
	for(int i = 0; i < numFaces(); ++i)
	{
		CPolygon p = FacePolygon(i);
		std::vector<CTriangle> tris = p.triangulate();
		int idx = vb.AppendVertices(3*(int)tris.size());
		for(size_t j = 0; j < tris.size(); ++j)
		{
			vb.set(idx, VDPosition, CVec4D(tris[j].a, 1.f));
			if (ccwIsFrontFacing)
			{
				vb.set(idx+1, VDPosition, CVec4D(tris[j].c, 1.f));
				vb.set(idx+2, VDPosition, CVec4D(tris[j].b, 1.f));
			}
			else
			{
				vb.set(idx+1, VDPosition, CVec4D(tris[j].b, 1.f));
				vb.set(idx+2, VDPosition, CVec4D(tris[j].c, 1.f));
			}
			// Generate flat normals if VB has space for normals.
			if (vb.Declaration()->TypeOffset(VDNormal) >= 0)
			{
				CVec3D normal = ccwIsFrontFacing ? tris[j].normalCCW() : tris[j].normalCW();
				vb.set(idx, VDNormal, CVec4D(normal, 0.f));
				vb.set(idx+1, VDNormal, CVec4D(normal, 0.f));
				vb.set(idx+2, VDNormal, CVec4D(normal, 0.f));
			}
			idx += 3;
		}
	}
}

void CPolyhedron::toLineList(VertexBuffer &vb) const
{
	std::vector<CLineSegment> edges = edges();

	int startIndex = vb.AppendVertices((int)edges.size()*2);
	for(int i = 0; i < (int)edges.size(); ++i)
	{
		vb.set(startIndex+2*i, VDPosition, CVec4D(edges[i].a, 1.f));
		vb.set(startIndex+2*i+1, VDPosition, CVec4D(edges[i].b, 1.f));
	}
}
#endif

CPolyhedron operator *(const CMat3D &transform, const CPolyhedron &polyhedron)
{
	CPolyhedron p(polyhedron);
	p.transform(transform);
	return p;
}

CPolyhedron operator *(const CMatJoint3x4 &transform, const CPolyhedron &polyhedron)
{
	CPolyhedron p(polyhedron);
	p.transform(transform);
	return p;
}

CPolyhedron operator *(const CMat4D &transform, const CPolyhedron &polyhedron)
{
	CPolyhedron p(polyhedron);
	p.transform(transform);
	return p;
}

CPolyhedron operator *(const CQuaternion &transform, const CPolyhedron &polyhedron)
{
	CPolyhedron p(polyhedron);
	p.transform(transform);
	return p;
}


} //end GEO
}  //end SMF


