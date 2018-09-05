#ifndef __SMF_POLYGON_
#define __SMF_POLYGON_

#include "../SMF_Config.h"
#include "../math/all.h"
#include "../util/SMF_RandomLCG.h"
#include "all.h"

namespace SMF {
using namespace MATH;
using namespace Util;
namespace GEO{



/**
 * \class CPolygon
 *
 * \ingroup SMF_Geometric
 * \image html pics\polygon.png
 *
 * \if pt_br
 * \brief   Representa um Polígono em espaço 3D
 * \elseif us_en
 * \brief  Represents a two-dimensional closed surface in 3D space.
 * A polygon is defined by N endpoints, or corner vertices. To be a valid polygon, there must be
   at least 3 vertices (a triangle).
   Well-formed polygons are always planar, i.e. all the vertices lie on the same plane. It is possible
   to store non-planar Polygons in this structure, but their representation is ambiguous, and for all practical
   purposes, should be avoided.
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CPolygon
{
public:
	/** 
	\brief  The default constructor creates a null polygon.
	 A null polygon has 0 vertices.
		\see isNull(). */
	CPolygon() {}

	/** 
	\brief  Stores the vertices of this polygon.
	*/
	std::vector<CVec3D> p;

	SMF_INLINE_FORCED static int numFaces() { return 1; }

	/** 
	\brief  Returns the number of edges in this polygon.
	 Since the polygon is always closed and connected, the number of edges is equal to the number of vertices.
		\see edge(), numVertices(). */
	int numEdges() const;

	/** 
	\brief  Returns the number of vertices in this polygon.
	 Since the polygon is always closed and connected, the number of edges is equal to the number of vertices.
		\see p, vertex(), numVertices(). */
	int numVertices() const;

	/** 
	\brief  Returns a pointer to an array of vertices of this polygon. The array contains numVertices() elements.
	\brief  \note Do NOT hold on to this pointer, since it is an alias to the underlying std::vector owned by this polygon. Calling any non-const CPolygon member function may invalidate the pointer!
	CVec3D *vertexArrayPtr() { return !p.empty() ? &p[0] : 0; }
	const CVec3D *vertexArrayPtr() const { return !p.empty() ? &p[0] : 0; }

	/** 
	\brief  Returns a vertex of this polygon.
	 \param vertexIndex The index of the vertex to get, in the range [0, numVertices()-1].
		\see p, numVertices(), edge(). */
	CVec3D vertex(int vertexIndex) const;

	/** 
	\brief  Returns a line segment between two adjacent vertices of this polygon.
	 \param edgeIndex The index of the edge line segment to construct, in the range [0, numEdges()-1].
		\return CLineSegment(vertex(edgeIndex), vertex((edgeIndex+1)%numVertices()).
		\see numEdges(), vertex(), edge2D(), edgenormal(), edgePlane(). */
	CLineSegment edge(int edgeIndex) const;

	/** 
	\brief  Returns a line segment between two adjacent vertices of this polygon, in the local space of this polygon.
	 In the local space of the polygon, the z-coordinate is always zero, and the polygon lies in the XY-plane, with
		the first vertex of the polygon being in the origin, and the x-axis running in the direction given by basisU() and
		the y-axis running in the direction given by basisV().
		\param edgeIndex The index of the edge line segment to construct, in the range [0, numEdges()-1].
		\see numEdges(), vertex(), edge2D(), edgenormal(), edgePlane(), basisU(), basisV(). */
	CLineSegment edge2D(int edgeIndex) const;

	/** 
	\brief  Returns the normal vector of the given edge.
	 The normal vector is perpendicular to the normal of the plane the polygon lies in, and the direction the given edge
		is pointing towards. The vector points outwards from the polygon.
		\param edgeIndex The index of the edge line segment to construct, in the range [0, numEdges()-1].
		\return A normalized direction vector perpendicular to the normal of the polygon, and the given edge.
		\see numEdges(), edge(), edge2D(), edgePlane(). */
	CVec3D edgenormal(int edgeIndex) const;

	/** 
	\brief  Returns the normal plane of the given edge.
	 The normal vector of the returned plane points towards the direction specified by edgenormal(), and the given edge
		lies inside the returned plane.
		\param edgeIndex The index of the edge line segment to construct, in the range [0, numEdges()-1]. */
	CPlane edgePlane(int edgeIndex) const;

	/** 
	\brief  Computes an extreme point of this CPolygon in the given direction.
	 An extreme point is a farthest point of this CPolygon in the given direction. Given a direction,
		this point is not necessarily unique.
		\param direction The direction vector of the direction to find the extreme point. This vector may
			be unnormalized, but may not be null.
		\return An extreme point of this CPolygon in the given direction. The returned point is always a
			vertex of this CPolygon.
		\see vertex(). */
	CVec3D extremePoint(const CVec3D &direction) const;

	/** 
	\brief  Projects this CPolygon onto the given 1D axis direction vector.
	 This function collapses this CPolygon onto an 1D axis for the purposes of e.g. separate axis test computations.
		The function returns a 1D range [outMin, outMax] denoting the interval of the projection.
		\param direction The 1D axis to project to. This vector may be unnormalized, in which case the output
			of this function gets scaled by the length of this vector.
		\param outMin [out] Returns the minimum extent of this object along the projection axis.
		\param outMax [out] Returns the maximum extent of this object along the projection axis. */
	void projectToAxis(const CVec3D &direction, float &outMin, float &outMax) const;

	/** 
	\brief  Tests if the given diagonal exists.
	 This function tests whether the diagonal that joins the two given vertices lies inside this polygon and is not intersected
		by the edges of this polygon.
		This function may only be called if this CPolygon is planar.
		\param i Index of the first endpoint of the diagonal CLineSegment, in the range [0, numVertices-1].
		\param j Index of the second endpoint of the diagonal CLineSegment, in the range [0, numVertices-1]. Do not pass in i==j
			or |i-j| == 1.
		\note If |i-j| <= 1, then SMF_ASSERT() failure is generated and false is returned.
		\return True if diagonal(vertexIndex1, vertexIndex2) exists and does not intersect the edges of this polygon.
		\see diagonal(), vertex(), edge(). */
	bool diagonalExists(int i, int j) const;

	/** 
	\brief  Returns the diagonal that joins the two given vertices.
	 If |i-j| == 1, then this returns an edge of this CPolygon.
		If i==j, then a degenerate line segment of zero length is returned.
		Otherwise, the line segment that joins the two given vertices is returned. Note that if the polygon is not planar or convex,
		this line segment might not lie inside the polygon. Use the diagonalExists() function to test whether the returned
		CLineSegment actually is a diagonal of this polygon.
		\param i Index of the first endpoint of the diagonal CLineSegment, in the range [0, numVertices-1].
		\param j Index of the second endpoint of the diagonal CLineSegment, in the range [0, numVertices-1].
		\note Whereas it is invalid to call diagonalExists() with values |i-j|<=1, it is acceptable for this function. This is to
			simplify generation of code that iterates over diagonal vertex pairs.	
		\return CLineSegment(vertex(i), vertex(j)) without checking if this actually is a valid diagonal of this polygon. If
			indices outside the valid range are passed, CLineSegment(nan, nan) is returned.
		\see vertex(), numVertices(), diagonalExists(). */
	CLineSegment diagonal(int i, int j) const;

	/** 
	\brief  Tests if this polygon is convex.
	 A polygon is convex, if for each pair of points inside the polygon, also the line segment joining those points is
		inside the polygon.
		\note This data structure can be used with both convex and non-convex polygons. In general, using convex polygons
			allows more efficient algorithms to be used with some operations. These more efficient variants are of form
			xxxConvex() in this class. Do not call those functions unless you know that the polygon is convex.
		\see isPlanar(), isSimple(), isNull(), isFinite(), isDegenerate(). */
	bool isConvex() const;

	/** 
	\brief  Tests if this polygon is planar.
	 A polygon is planar if all its vertices lie on the same plane.
		\note Almost all functions in this class require that the polygon is planar. While you can store vertices of
			non-planar polygons in this class, they are better avoided. Read the member function documentation carefully
			to avoid calling for non-planar polygons any functions which SMF_ASSERT planarity.
		\see isConvex(), isSimple(), isNull(), isFinite(), isDegenerate(). */
	bool isPlanar(float epsilon =CMath::EPSILON_SuperLow) const;

	/** Tests if this polygon is simple.
		A polygon is simple if no two nonconsecutive edges have a point in common.
		In other words, a planar polygon is simple if its edges do not self-intersect, and if each vertex is joined by
		exactly two edges.
		\note This function assumes that the polygon is planar.
		\see isConvex(), isPlanar(), isNull(), isFinite(), isDegenerate(). */
	bool isSimple() const;

	/** 
	\brief  Tests if this polygon is null.
	 A polygon is null if it has zero vertices.
		\note The null polygon is degenerate and finite.
		\see p, isConvex(), isPlanar(), isSimple(), isFinite(), isDegenerate(). */
	bool isNull() const;

	/** 
	\brief  Tests if this polygon is finite.
	 A polygon is finite if each of its vertices have finite floating point coordinates (no nans or infs).
		\note The null polygon is finite.
		\see p, isConvex(), isPlanar(), isSimple(), isNull(), isDegenerate(), ::isFinite(), isInf(), IsNan(), isFinite(), inf, negInf, nan, CVec3D::nan, CVec3D::inf. */
	bool isFinite() const;

	/** 
	\brief  Tests if this polygon is degenerate.
	 A polygon is degenerate if it has two or less vertices, or if its surface area is less or equal than the given epsilon.
		\note The null polygon is degenerate and finite.
		\see p, isConvex(), isPlanar(), isSimple(), isNull(), isDegenerate(), area(). */
	bool isDegenerate(float epsilon =CMath::EPSILON_SuperLow) const;

	/** 
	\brief  Generates the U vector of the local space of this polygon.
	 This vector specifies in global (world) space the direction of the local X axis of this polygon.
		\note This function assumes that the first two points (p[0] and p[1]) of this polygon are finite and inequal. */
	CVec3D basisU() const;
	/** 
	\brief  Generates the V vector of the local space of this polygon. [similarOverload: basisU]
	 This vector specifies in global (world) space the direction of the local Y axis of this polygon.
		\note This function assumes that the first two points (p[0] and p[1]) of this polygon are finite and inequal.
		\see mapTo2D(), mapFrom2D(), edge2D(), basisU(), basisV(). */
	CVec3D basisV() const;

	/** 
	\brief  Returns the given vertex of this polygon mapped to a local 2D space on this polygon.
	 In the local space of the polygon, the z-coordinate is always zero, and the polygon lies in the XY-plane, with
		the first vertex of the polygon being in the origin, and the x-axis running in the direction given by basisU() and
		the y-axis running in the direction given by basisV().
		\param i The index of the vertices of this polygon to generate, in the range [0, numVertices()-1].
		\see numVertices(), mapFrom2D(), edge2D(), basisU(), basisV(). */
	CVec2D mapTo2D(int i) const;

	/** 
	\brief  Maps the given global (world) space point to the local 2D space of this polygon.
	\brief  \todo Return a CVec3D to be able to read the distance of the point from the plane of the polygon? (or add an overload for that)
	\brief  \todo add mapTo2D(CLine/CLineSegment/CRay/CTriangle/CPolygon).*/
	CVec2D mapTo2D(const CVec3D &point) const;

	/** 
	\brief  Given a 2D point in the local space, returns the corresponding 3D point in the global (world) space.
	 \see mapTo2D(), basisU(), basisV(). */
	CVec3D mapFrom2D(const CVec2D &point) const;

	/** 
	\brief  Computes the normal of this polygon.
	 \return The normal of this polygon. This vector is normalized and points to the direction from which observed the
		vertices of this polygon wind in counter-clockwise order.
		\note Only call this function if this CPolygon is planar. */
	CVec3D normalCCW() const;
	/** 
	\brief  Computes the normal of this polygon in clockwise direction. [similarOverload: normalCCW]
	 \return The normal of this polygon in clockwise direction. This vector is normalized and points to the direction
		from which observed the vertices of this polygon wind in clockwise order.
		\note Only call this function if this CPolygon is planar.
		\note These two functions follow the relation normalCCW() == -normalCW().
		\see planeCW(), planeCCW(). */
	CVec3D normalCW() const;

	/** 
	\brief  Computes the plane this polygon is contained in.
	 \note Only call this function if this CPolygon is planar.
		\return The plane equation of this polygon. This normal vector of the plane points to the direction from which observed the
		vertices of this polygon wind in counter-clockwise order. */
	CPlane planeCCW() const;
	/** 
	\brief  Computes the plane this polygon is contained in, with a normal vector that points in the clockwise direction. [similarOverload: planeCCW]
	 \note Only call this function if this CPolygon is planar.
		\note The functions planeCCW() and planeCW() return the same plane, except the normals of the planes point in opposite directions.
		\see normalCCW(), normalCW(). */
	CPlane planeCW() const;

	/** 
	\brief  Translates this CPolygon in world space.
	 \param offset The amount of displacement to apply to this CPolygon, in world space coordinates.
		\see transform(). */
	void translate(const CVec3D &offset);

	/** 
	\brief  Applies a transformation to this CPolygon.
	 This function operates in-place.
		\see translate(), classes CMat3D, CMatJoint3x4, CMat4D, CQuaternion. */
	void transform(const CMat3D &transform);
	void transform(const CMatJoint3x4 &transform);
	void transform(const CMat4D &transform);
	void transform(const CQuaternion &transform);

	// Returns true if the edges of this polygon self-intersect.
//	bool IsSelfIntersecting() const;

	// Projects all vertices of this polygon to the given plane.
//	void ProjectToPlane(const CPlane &plane);

	// Returns true if the edges of this polygon self-intersect when viewed from the given direction.
//	bool IsSelfIntersecting(const CVec3D &viewDirection) const;

	// Returns true if there exists edges (p_{i-1}, p_i) and (p_i, p_{i+1}) which are collinear.
//	bool HasCollinearEdges() const;

	/** 
	\brief  Tests if the given object, expressed in global (world) space, is fully contained inside this polygon.
	 Only call this function if the polygon is planar.
		This test is performed in global space of this polygon, i.e. by specifying the other object in global (world)
		space coordinates.
		\param point The point to test for containment.
		\param polygonThickness Since a polygon is a 2D object in a 3D space, a threshold value is used to
			allow floating-point inaccuracies. This parameter defines how much "thickness" to give to the polygon
			for the purposes of the test.
		\return True if the given object is fully contained inside this polygon (and the plane of this polygon).
		\todo add containsConvex(CVec3D/etc.). See RTCD p. 202.
		\todo add contains(CCircle/Disc). */
	bool contains(const CVec3D &point, float polygonThickness = 1e-3f) const;
	bool contains(const CLineSegment &lineSegment, float polygonThickness = 1e-3f) const;
	bool contains(const CTriangle &triangle, float polygonThickness = 1e-3f) const;
	bool contains(const CPolygon &polygon, float polygonThickness = 1e-3f) const;
	//todo add RTCD, p. 202.
	//bool containsConvex(const CVec3D &worldSpacePoint, float polygonThickness = 1e-3f) const;

	/** 
	\brief  Tests if the given object, expressed in the local space of this polygon, is fully contained inside this polyhedron.
	 This test is exactly like in contains(), except it is performed in 2D in the local space of this polygon.
		\see contains(), mapTo2D().
		\todo add contains2D(CCircle/Disc/CTriangle/CPolygon). */
	bool contains2D(const CLineSegment &localSpaceLineSegment) const;

	/** 
	\brief  Tests whether this polyhedron and the given object intersect.
	 Both objects are treated as "solid", meaning that if one of the objects is fully contained inside
		another, this function still returns true.
		This test is performed in the global (world) space of this polygon.
		\return True if an intersection occurs or one of the objects is contained inside the other, false otherwise.
		\see contains(), closestPoint(), distance().
		\todo add intersects(CCircle/Disc). */
	bool intersects(const CLine &line) const;
	bool intersects(const CRay &ray) const;
	bool intersects(const CLineSegment &lineSegment) const;
	bool intersects(const CPlane &plane) const;
	bool intersects(const CAABBox &aabb) const;
	bool intersects(const COBBox &obb) const;
	bool intersects(const CTriangle &triangle) const;
	bool intersects(const CPolygon &polygon) const;
	bool intersects(const CPolyhedron &polyhedron) const;
	bool intersects(const CSphere &sphere) const;

	/** 
	\brief  Computes the closest point on this polygon to the given object.
	 If the other object intersects this polygon, this function will return an arbitrary point inside
		the region of intersection.
		\param lineSegment The line segment to find the closest point to.
		\param lineSegmentPt [out] If specified, receives the closest point on the line segment to this polygon. This
			pointer may be null.
		\see contains(), distance(), intersects().
		\todo add closestPoint(CLine/CRay/CPlane/CTriangle/CPolygon/CCircle/Disc/CAABBox/COBBox/CSphere/Capsule/Frustum/CPolyhedron). */
	CVec3D closestPoint(const CLineSegment &lineSegment, CVec3D *lineSegmentPt) const;
	CVec3D closestPoint(const CLineSegment &lineSegment) const;
	CVec3D closestPoint(const CVec3D &point) const;

	/** 
	\brief  Returns the distance between this polygon and the given point.
	 \see contains(), closestPoint(), intersects(). */
	float distance(const CVec3D &point) const;

	/** 
	\brief  Returns the surface area of this polygon.
	/** \see perimeter(), centroid(). */
	float area() const;

	/** 
	\brief  Returns the total edge length of this polygon.
	 \see area(), centroid(). */
	float perimeter() const;

	/** 
	\brief  Returns the center of mass of this polygon.
	 \see area(), perimeter(). */
	CVec3D centroid() const;
	/** 
	\brief  Identical to centerPoint(), but provided to enable common signature with CTriangle, CAABBox and COBBox to allow them to be used
	 in template classes. **/
	CVec3D centerPoint() const { return centroid(); }

	/** 
	\brief  Computes a point on the perimeter of this polygon.
	 \param normalizedDistance A value in the range [0,1[ specifying the distance along the polygon edge to travel.
		The polygon perimeter forms a closed loop, so pointOnEdge(0.f) == pointOnEdge(1.f) and is equal to the point p[0] of this
		polygon. As another example, pointOnEdge(0.5f) returns the point half-way around the polygon edge (but not necessarily the farthest
		point from p[0]).
		\see p, randomPointOnEdge(). */
	CVec3D pointOnEdge(float normalizedDistance) const;

	/** 
	\brief  Computes a random point on the perimeter of this polygon.
	 This function generates points with uniform distribution.
		\see pointOnEdge(). */
	CVec3D randomPointOnEdge(CRandomLCG &rng) const;

	CVec3D fastRandomPointInside(CRandomLCG &rng) const;

	/** 
	\brief  Converts this CPolygon to a CPolyhedron representation.
	 This function will create a CPolyhedron with two faces, one for the front face of this CPolygon,
		and one for the back face.
		\todo add toPolyhedron(float polygonThickness)
		\see triangulate(), minimalEnclosingAABB(). */
	CPolyhedron toPolyhedron() const;

	// These faces will be extruded along the CPolygon normal so that they lie polygonThickness units apart from each other.
//	CPolyhedron toPolyhedron(float polygonThickness = 0.f) const; 
	/** \todo add support for this form.*/

	/** 
	\brief Triangulates this CPolygon using the ear-clipping method.
	\see toPolyhedron(), minimalEnclosingAABB(). 
	*/
	std::vector<CTriangle> triangulate() const;

	/** 
	\brief Returns the smallest CAABBox that encloses this polygon.
	\todo add minimalEnclosingSphere() and minimalEnclosingOBB().
	\see toPolyhedron(), triangulate(). 
	*/
	CAABBox minimalEnclosingAABB() const;

	// Returns true if the given vertex is a concave vertex. Otherwise the vertex is a convex vertex.
//	bool IsConcaveVertex(int i) const;

	// Computes the conves hull of this polygon.
//	CPolygon convexHull() const;

//	bool IsSupportingPoint(int i) const;

//	bool IsSupportingPoint(const CVec3D &point) const;

	// Returns true if the quadrilateral defined by the four points is convex (and not concave or bowtie).
//	static bool IsConvexQuad(const CVec3D &pointA, const CVec3D &pointB, const CVec3D &pointC, const CVec3D &pointD);
};

CPolygon operator *(const CMat3D &transform, const CPolygon &polygon);
CPolygon operator *(const CMatJoint3x4 &transform, const CPolygon &polygon);
CPolygon operator *(const CMat4D &transform, const CPolygon &polygon);
CPolygon operator *(const CQuaternion &transform, const CPolygon &polygon);

#ifdef MATH_QT_INTEROP
Q_DECLARE_METATYPE(CPolygon)
Q_DECLARE_METATYPE(CPolygon*)
#endif

// \todo add this.
//std::ostream &operator <<(std::ostream &o, const CPolygon &polygon);


} //end GEO
}  //end SMF

#endif // __SMF_POLIGON_
