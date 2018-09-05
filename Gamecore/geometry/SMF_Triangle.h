#ifndef __SMF_TRIANGLE_
#define __SMF_TRIANGLE_

#include "../SMF_Config.h"
#include "../math/all.h"
#include "all.h"

namespace SMF {
using namespace MATH;
using namespace Util;
namespace GEO{


/**
 * \class CTriangle
 *
 * \ingroup SMF_Geometric
 * \image html pics/triangle.png
 * \if pt_br
 * \brief   Um triângulo em espaço 3D
 * \note Armazena os três pontos vértices a,b,c que definem um triângulo 
 * \elseif us_en
 * \brief 	Specifies a triangle through three points in 3D space.
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CTriangle
{
public:
	/** 
	\brief  The first triangle endpoint.*/
	CVec3D a;
	/** 
	\brief  The second triangle endpoint.
	\note [similarOverload: a] */
	CVec3D b;
	/** 
	\brief  The third triangle endpoint.
	\note [similarOverload: a]
		\note The order in which the vertices are stored in this data structure is important. The triangles
			(a,b,c) and (a,c,b) are otherwise equivalent, but their plane normals point to the opposite directions.
		\see planeCCW(), planeCW(). */
	CVec3D c;

	/** 
	\brief  The default constructor does not initialize any members of this class.
	\note This means that the values of the members a, b and c are undefined after creating a new CTriangle using this
		default constructor. Remember to assign to them before use.
		\see a, b, c. */
	CTriangle() {}

	/** 
	\brief  Constructs a triangle from three given endpoints.
	\note The normal of the plane will be constructed to point towards the halfspace where
		the vertices a, b and c wind in counter-clockwise order.
		\see a, b, c. */
	CTriangle(const CVec3D &a, const CVec3D &b, const CVec3D &c);

	SMF_INLINE_FORCED static int numFaces() { return 1; }
	SMF_INLINE_FORCED static int numEdges() { return 3; }
	SMF_INLINE_FORCED static int numVertices() { return 3; }

	/** 
	\brief  Translates this CTriangle in world space.
	\param offset The amount of displacement to apply to this CTriangle, in world space coordinates.
	\see transform(). */
	void translate(const CVec3D &offset);

	/** 
	\brief Applies a transformation to this CTriangle, in-place.
	See translate(), classes CMat3D, CMatJoint3x4, CMat4D, CQuaternion. */
	void transform(const CMat3D &transform);
	void transform(const CMatJoint3x4 &transform);
	void transform(const CMat4D &transform);
	void transform(const CQuaternion &transform);

	/** 
	\brief Expresses the given point in barycentric (u,v,w) coordinates.
	\note There are two different conventions for representing barycentric coordinates. One uses
			a (u,v,w) triplet with the equation pt == u*a + v*b + w*c, and the other uses a (u,v) pair
			with the equation pt == a + u*(b-a) + v*(c-a). These two are equivalent. Use the mappings
			(u,v) -> (1-u-v, u, v) and (u,v,w)->(v,w) to convert between these two representations.
		\param point The point of the vector space to express in barycentric coordinates. This point should
			lie in the plane formed by this triangle.
		\return The factors (u,v,w) that satisfy the weighted sum equation point == u*a + v*b + w*c.
		\see barycentricUV(), barycentricInsideTriangleboundingAABB(), Point(), http://mathworld.wolfram.com/BarycentricCoordinates.html */
	CVec3D barycentricUVW(const CVec3D &point) const;

	/** 
	\brief  Expresses the given point in barycentric (u,v) coordinates.
	\note There are two different conventions for representing barycentric coordinates. One uses
			a (u,v,w) triplet with the equation pt == u*a + v*b + w*c, and the other uses a (u,v) pair
			with the equation pt == a + u*(b-a) + v*(c-a). These two are equivalent. Use the mappings
			(u,v) -> (1-u-v, u, v) and (u,v,w)->(v,w) to convert between these two representations.
		\param point The point to express in barycentric coordinates. This point should lie in the plane
			formed by this triangle.
		\return The factors (u,v) that satisfy the weighted sum equation point == a + u*(b-a) + v*(c-a).
		\see barycentricUVW(), barycentricInsideTriangleboundingAABB(), Point(). */
	CVec2D barycentricUV(const CVec3D &point) const;

	/** 
	\brief  Tests if the given barycentric UVW coordinates lie inside a triangle.
	\note A barycentric UVW coordinate represents a point inside a triangle if
		  a) 0 <= u,v,w <= 1 and
		  b) u+v+w == 1.
	\param uvw The barycentric vector containing the barycentric (u,v,w) coordinates.
	\return True if the given coordinates lie inside a triangle.
	\see barycentricUV(), barycentricUVW(), Point(). */
	static bool barycentricInsideTriangleboundingAABB(const CVec3D &uvw);

	/** 
	\brief  Returns the point at the given barycentric coordinates.
	\note This function computes the vector space point at the given barycentric coordinates.
	\param uvw The barycentric UVW coordinate triplet. The condition u+v+w == 1 should hold for the input coordinate.
			If 0 <= u,v,w <= 1, the returned point lies inside this triangle.
	\return u*a + v*b + w*c. */
	CVec3D Point(const CVec3D &uvw) const;
	CVec3D Point(float u, float v, float w) const;
	/** These functions are an alternate form of Point(u,v,w) for the case when the barycentric coordinates are
		represented as a (u,v) pair and not as a (u,v,w) triplet. This function is provided for convenience
		and effectively just computes Point(1-u-v, u, v).
		\param uv The barycentric UV coordinates. If 0 <= u,v <= 1 and u+v <= 1, then the returned point lies inside
			this triangle.
		\return a + (b-a)*u + (c-a)*v.
		\see barycentricUV(), barycentricUVW(), barycentricInsideTriangleboundingAABB(). */
 	CVec3D Point(const CVec2D &uv) const;
	CVec3D Point(float u, float v) const;

	/** 
	\brief  Computes the center of mass of this triangle.
	\return The point (a+b+c)/3. */
	CVec3D centroid() const;
	/// Identical to centerPoint(), but provided to enable common signature with CTriangle, CAABBox and COBBox to allow them to be used
	CVec3D getCenter() const { return centroid(); }

	/** 
	\brief  Computes the surface area of this triangle.
	\return The surface area of this triangle.
	\see perimeter(), area2D(), signedArea(). */
	float area() const;

	/** 
	\brief  Computes the total edge length of this triangle.
	\return |a-b| + |b-c| + |c-a|.
	\see area(), edge(). */
	float perimeter() const;

	/** 
	\brief  Returns a pointer to the vertices of this triangle. The array contains three elements.
	*/
	CVec3D *vertexArrayPtr() { return &a; }
	const CVec3D *vertexArrayPtr() const { return &a; }

	/** 
	\brief  Returns a vertex of this triangle.
	\param i The vertex of this triangle to get: 0, 1 or 2.
	\return vertex(0) returns the point a, vertex(1) returns the point b, and vertex(2) returns the point c.
	\note If an index outside [0, 2] is passed, an SMF_ASSERT() failure occurs and CVec3D(NaN) is returned.
	\see edge(). */
	CVec3D vertex(int i) const;

	
	/** 
	\brief Returns an edge of this triangle.
	\param i The index of the edge to generate: 0, 1 or 2.
		\return A CLineSegment representing the given edge of this triangle. edge(0) returns CLineSegment(a,b), edge(1)
			returns CLineSegment(b,c) and edge(2) returns CLineSegment(c,a).
		\note If an index outside [0, 2] is passed, an SMF_ASSERT() failure occurs and CLineSegment(NaN, NaN) is returned.
		\see vertex(). */
	CLineSegment edge(int i) const;

	
	/** 
	\brief Returns the counterclockwise-oriented plane this triangle lies on.
	The normal of the returned plane points towards the halfspace in which the vertices of this triangle are winded
		in <b>counter-clockwise</b> direction. */
	CPlane planeCCW() const;

	/** 
	\brief Returns the clockwise-oriented plane this triangle lies on.
	 The normal of the returned plane points towards the halfspace in which the vertices of this triangle are winded
		in <b>clockwise</b> direction. [similarOverload: planeCCW]
		\see normalCCW(), normalCW(), unNormalizedNormalCCW(), unNormalizedNormalCW(). */
	CPlane planeCW() const;

	/** 
	\brief  Returns the normalized triangle normal pointing to the counter-clockwise direction of this triangle.
	 This function computes the normalized triangle normal vector that points towards the halfspace,
		from which observed, the vertices of this triangle wind in <b>counter-clockwise</b> (CCW) order.
		\see planeCCW(), planeCW(), unNormalizedNormalCCW(), unNormalizedNormalCW(). */
	CVec3D normalCCW() const;

	/** 
	\brief  Returns the normalized triangle normal pointing to the clockwise direction of this triangle.
	 This function computes the normalized triangle normal vector that points towards the halfspace,
		from which observed, the vertices of this triangle wind in <b>clockwise</b> (CW) order. [similarOverload: normalCCW]
		\see planeCCW(), planeCW(), unNormalizedNormalCCW(), unNormalizedNormalCW(). */
	CVec3D normalCW() const;

	/// Computes an unnormalized counter-clockwise oriented triangle normal vector.
	CVec3D unNormalizedNormalCCW() const;
	/** 
	\brief  Computes an unnormalized clockwise-oriented triangle normal vector.
	 These functions are equivalent to normalCCW() and normalCW(), except these functions do not produce
		a unit-length triangle normal vectors. Use these functions instead of normalCCW/CW() to obtain a small
		speed benefit in cases where the normalization step is not required. [similarOverload: unNormalizedNormalCCW]
		\see planeCCW(), planeCW(), normalCCW(), normalCW(). */
	CVec3D unNormalizedNormalCW() const;

	/** 
	\brief  Computes an extreme point of this CTriangle in the given direction.
	 An extreme point is a farthest point of this CTriangle in the given direction. Given a direction,
		this point is not necessarily unique.
		\param direction The direction vector of the direction to find the extreme point. This vector may
			be unnormalized, but may not be null.
		\return An extreme point of this CTriangle in the given direction. The returned point is always a
			vertex of this CTriangle.
		\see vertex(). */
	CVec3D extremePoint(const CVec3D &direction) const;

	/** 
	\brief  Returns a CPolygon representation of this CTriangle.
	 The returned polygon is identical to this CTriangle. It has three vertices a, b and c which wind in the same
		direction than in this triangle.
	\see class CPolygon, toPolyhedron(). */
	CPolygon toPolygon() const;

	/** 
	\brief  Returns a CPolyhedron representation of this CTriangle.
	 The generated polyhedron will be closed and has two triangular faces, and three vertices (a, b and c).
		The two faces share the same vertices, but in opposite winding orders. This creates a polyhedron with zero
		volume and the surface area twice of this CTriangle.
		\see class CPolyhedron, toPolygon(). */
	CPolyhedron toPolyhedron() const;

	/// Returns the tight CAABBox that encloses this CTriangle.
	CAABBox boundingAABB() const;

	/** 
	\brief  Computes the surface area of the given 2D triangle.
	 This math library does not have a separate class for 2D triangles. To compute the area of a 2D triangle,
		use this CTriangle class and set z=0 for each coordinate, or use this helper function.
		\see area(), signedArea(). */
	static float area2D(const CVec2D &p1, const CVec2D &p2, const CVec2D &p3);

	/** 
	\brief  Relates the barycentric coordinate of the given point to the surface area of triangle abc.
	 This function computes the ratio of the signed area of the triangle (point, b, c) to the signed area of
		the triangle (a, b, c). This is the same as computing the barycentric u-coordinate of the given point
		on the triangle (a, b, c).
		\see area(), area2D(), barycentricUVW(). */
	static float signedArea(const CVec3D &point, const CVec3D &a, const CVec3D &b, const CVec3D &c);

	/** 
	\brief  Tests if this CTriangle is finite.
	 A triangle is <b><i>finite</i></b> if its vertices a, b and c do not contain floating-point NaNs or +/-infs
		in them.
		\return True if each coordinate of each vertex of this triangle has a finite floating-point value.
		\see a, b, c, isDegenerate(), ::isFinite(), isInf(), IsNan(), isFinite(), inf, negInf, nan, CVec3D::nan, CVec3D::inf. */
	bool isFinite() const;

	/** 
	\brief  Returns true if this triangle is degenerate.
	 A triangle is <b><i>degenerate</i></b> if it is not finite, or if the surface area of the triangle is
		close to zero.
		\param epsilon The threshold to test against. If two vertices of this triangle are closer than this, the
		triangle is considered degenerate.
		\see a, b, c, isFinite(). */
	bool isDegenerate(float epsilon =CMath::EPSILON_SuperLow) const;

	/// Returns true if the triangle defined by the three given points is degenerate.
	static bool isDegenerate(const CVec3D &p1, const CVec3D &p2, const CVec3D &p3, float epsilon =CMath::EPSILON_SuperLow);

	/** 
	\brief  Tests if the given object is fully contained inside this triangle.
	\param triangleThickness An epsilon threshold value to use for this test.
			This specifies the maximum distance the given object can lie from the plane defined by this triangle.
	\see distance(), intersects(), closestPoint().
	\todo add CTriangle::contains(CCircle) and CTriangle::contains(Disc). */
	bool contains(const CVec3D &point, float triangleThickness = 1e-3f) const;
	bool contains(const CLineSegment &lineSegment, float triangleThickness = 1e-3f) const;
	bool contains(const CTriangle &triangle, float triangleThickness = 1e-3f) const;

	/** 
	\brief  Computes the distance between this triangle and the given object.
	 This function finds the nearest pair of points on this and the given object, and computes their distance.
		If the two objects intersect, or one object is contained inside the other, the returned distance is zero.
	\todo add CTriangle::distance(CLine/CRay/CLineSegment/CPlane/CTriangle/CPolygon/CCircle/Disc/CAABBox/COBBox/Capsule/Frustum/CPolyhedron).
	\see contains(), intersects(), closestPoint(). */
	float distance(const CVec3D &point) const;
	float distance(const CSphere &sphere) const;
//	float distance(const Capsule &capsule) const;

	/** 
	\brief  Tests whether this triangle and the given object intersect.	
	 Both objects are treated as "solid", meaning that if one of the objects is fully contained inside
		another, this function still returns true. (e.g. in case a line segment is contained inside this triangle,
		or this triangle is contained inside a sphere, etc.)
	\param d [out] If specified, this parameter will receive the parametric distance of
			the intersection point along the line object. Use the getPoint(d) function of the line class
			to get the actual point of intersection. This pointer may be null.
	\param intersectionPoint [out] If specified, receives the actual point of intersection. This pointer
			may be null.
	\return True if an intersection occurs or one of the objects is contained inside the other, false otherwise.
	\see contains(), distance(), closestPoint(), CLineSegment::getPoint(). */
	bool intersects(const CLineSegment &lineSegment, float *d = 0, CVec3D *intersectionPoint = 0) const;
	bool intersects(const CLine &line, float *d = 0, CVec3D *intersectionPoint = 0) const;
	bool intersects(const CRay &ray, float *d = 0, CVec3D *intersectionPoint = 0) const;
	bool intersects(const CPlane &plane) const;
	/** \param closestPointOnTriangle [out] If specified, receives the point of intersection between the CSphere
			and this CTriangle. Even if no intersection occurred, this parameter will receive the closest point on
			the CTriangle to the CSphere. This pointer may be null. */
	bool intersects(const CSphere &sphere, CVec3D *closestPointOnTriangle) const;
	bool intersects(const CSphere &sphere) const;
	/** \param outLine [out] If specified, receives the line segment of the common points shared by the two
			intersecting triangles. If the two triangles do not intersect, this pointer is not written to.
			This pointer may be null. */
	bool intersects(const CTriangle &triangle, CLineSegment *outLine = 0) const;
	bool intersects(const CAABBox &aabb) const;
	bool intersects(const COBBox &obb) const;
	bool intersects(const CPolygon &polygon) const;
	bool intersects(const CPolyhedron &polyhedron) const;

	/// A helper function used in line-triangle tests.
	static float intersectLineTri(const CVec3D &linePos, const CVec3D &lineDir,
		const CVec3D &v0, const CVec3D &v1, const CVec3D &v2,
		float &u, float &v);

	/** 
	\brief  Projects this CTriangle onto the given axis.
	 This function is used in SAT tests (separate axis theorem) to check the interval this triangle
		lies in on an 1D line.
		\param axis The axis to project onto. This vector can be unnormalized.
		\param dMin [out] Returns the minimum extent of this triangle on the given axis.
		\param dMax [out] Returns the maximum extent of this triangle on the given axis. */
	void projectToAxis(const CVec3D &axis, float &dMin, float &dMax) const;

	/** 
	\brief  Computes the closest point on this triangle to the given object.
	 If the other object intersects this triangle, this function will return an arbitrary point inside
		the region of intersection.
		\see contains(), distance(), intersects(), closestPointToTriangleEdgesignedArea(). */
	CVec3D closestPoint(const CVec3D &point) const;
	/** \param otherPt [out] If specified, receives the closest point on the other object to this triangle.
		This pointer may be null. */
	CVec3D closestPoint(const CLineSegment &lineSegment, CVec3D *otherPt = 0) const;
	/** 
	\param outU [out] If specified, receives the barycentric U coordinate of the returned point (in the UV convention).
			This pointer may be null. TODO add this parameter back.
	\param outV [out] If specified, receives the barycentric V coordinate of the returned point (in the UV convention).
			This pointer may be null. TODO add this parameter back.
	\param outD [out] If specified, receives the distance along the line of the closest point on the line to this triangle. TODO add this parameter back.
	\return The closest point on this triangle to the given object.
	\todo add closestPoint(CRay/CPlane/CPolygon/CCircle/Disk/CAABBox/COBBox/CSphere/Capsule/Frustum/CPolyhedron).
	\see distance(), contains(), intersects(), closestPointToTriangleEdgesignedArea(), CLine::getPoint. */
	CVec3D closestPoint(const CLine &line, CVec3D *otherPt = 0) const;
	CVec3D closestPoint(const CTriangle &triangle, CVec3D *otherPt = 0) const;

	/** 
	\brief  Computes the closest point on the edge of this triangle to the given object.
	\param outU [out] If specified, receives the barycentric U coordinate of the returned point (in the UV convention).
			This pointer may be null.
	\param outV [out] If specified, receives the barycentric V coordinate of the returned point (in the UV convention).
			This pointer may be null.
	\param outD [out] If specified, receives the distance along the line of the closest point on the line to the edge of this triangle.
	\return The closest point on the edge of this triangle to the given object.
	\todo add closestPointToTriangleEdgesignedArea(Point/CRay/CTriangle/CPlane/CPolygon/CCircle/Disk/CAABBox/COBBox/CSphere/Capsule/Frustum/CPolyhedron).
	\see distance(), contains(), intersects(), closestPointToTriangleEdgesignedArea(), CLine::getPoint. */
	CVec3D closestPointToTriangleEdgesignedArea(const CLine &line, float *outU, float *outV, float *outD) const;
	CVec3D closestPointToTriangleEdgesignedArea(const CLineSegment &lineSegment, float *outU, float *outV, float *outD) const;

	/** 
	\brief  Generates a random point inside this CTriangle.
	 The points are distributed uniformly.
		The implementation of this function is based on Graphics Gems 1, p. 25:
		"1.5 Generating random points in triangles. Method 2." The Method 1 presented in the book
		uses a sqrt() instead of the if().
		\param rng A pre-seeded random number generator object that is to be used by this function to generate random values.
		\see class CRandomLCG, randomPointOnEdge(), randomVertex(), Point(). */
	CVec3D randomPointInside(CRandomLCG &rng) const;

	/** 
	\brief  Chooses a corner vertex of this CTriangle at random.
	This function returns one of the vertices {a, b, c} at uniformly random.
	\param rng A pre-seeded random number generator object that is to be used by this function to generate random values.
	\see class CRandomLCG, randomPointInside(), randomPointOnEdge(), vertex(). */
	CVec3D randomVertex(CRandomLCG &rng) const;

	/** 
	\brief  Generates a random point on the edge of this CTriangle.
	 The points are distributed uniformly.
		This function requires that this triangle is not degenerate. If it is, an SMF_ASSERT() error will be printed,
		and the return value will be undefined.
		\param rng A pre-seeded random number generator object that is to be used by this function to generate random values.
		\see class CRandomLCG, randomPointInside(), randomVertex(), edge(), class CLineSegment, isDegenerate(). */
	CVec3D randomPointOnEdge(CRandomLCG &rng) const;

	/// Returns a human-readable representation of this CLine. Most useful for debugging purposes.
	std::string toString() const;

#ifdef MATH_QT_INTEROP
	operator QString() const { return toString(); }
	QString toString() const { return QString::fromStdString(toString()); }
#endif
};

CTriangle operator *(const CMat3D &transform, const CTriangle &t);
CTriangle operator *(const CMatJoint3x4 &transform, const CTriangle &t);
CTriangle operator *(const CMat4D &transform, const CTriangle &t);
CTriangle operator *(const CQuaternion &transform, const CTriangle &t);

#ifdef MATH_QT_INTEROP
Q_DECLARE_METATYPE(CTriangle)
Q_DECLARE_METATYPE(CTriangle*)
#endif


std::ostream &operator <<(std::ostream &o, const CTriangle &triangle);



} //end GEO
}  //end SMF

#endif // __SMF_TRIANGLE_
