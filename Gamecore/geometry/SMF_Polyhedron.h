#ifndef __SMF_POLYHEDRON_
#define __SMF_POLYHEDRON_

#include "../SMF_Config.h"
#include "../math/all.h"
#include "all.h"

namespace SMF {
using namespace MATH;
using namespace Util;
namespace GEO{



 
/**
 * \class CPolyhedron
 *
 * \ingroup SMF_Geometric
 * \image html pics\polyhedron.png
 * \if pt_br
 * \brief Representa um sólido geométrico definido por faces poligonais planas num espaço 3D
 *       
 * \elseif us_en
 * \brief 	Represents a three-dimensional closed geometric solid defined by flat polygonal faces.
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CPolyhedron
{
public:
	/// Stores a list of indices of a single face of a CPolyhedron.
	struct Face
	{
		/// Specifies the indices of the corner vertices of this polyhedron.
		/// Indices point to the polyhedron vertex array.
		/// The face vertices should all lie on the same plane.
		/// The positive direction of the plane (the direction the face outwards normal points to)
		/// is the one where the vertices are wound in counter-clockwise order.
		std::vector<int> v;

		/// Reverses the winding order of this face. This has the effect of reversing the direction
		/// the normal of this face points to.
		void flipWindingOrder();

		std::string toString() const;
	};

	/// Specifies the vertices of this polyhedron.
	std::vector<CVec3D> v;

	/// Specifies the individual faces of this polyhedron.  [similarOverload: v]
	/** Each face is described by a list of indices to the vertex array. The indices define a
		simple polygon in counter-clockwise winding order. */
	std::vector<Face> f;

	/// The default constructor creates a null polyhedron.
	/** A null polyhedron has 0 vertices and 0 faces.
		\see isNull(). */
	CPolyhedron() {}

	/// Returns the number of vertices in this polyhedron.
	/** The running time of this function is O(1).
		\see numFaces(), numEdges(), eulerFormulaHolds(). */
	int numVertices() const { return (int)v.size(); }

	/// Returns the number of faces in this polyhedron.
	/** The running time of this function is O(1).
		\see numVertices(), numEdges(), eulerFormulaHolds(), FacePolygon(), facePlane(). */
	int numFaces() const { return (int)f.size(); }

	/// Returns the number of (unique) edges in this polyhedron.
	/** This function will enumerate through all faces of this polyhedron to compute the number of unique edges.
		The running time is linear to the number of faces and vertices in this CPolyhedron.
		\see numVertices(), numFaces(), eulerFormulaHolds(), edge(), edges(), edgeIndices(). */
	int numEdges() const;

	/// Returns a pointer to an array of vertices of this polyhedron. The array contains numVertices() elements.
	/// \note Do NOT hold on to this pointer, since it is an alias to the underlying std::vector owned by this polyhedron. Calling any non-const CPolyhedron member function may invalidate the pointer!
	CVec3D *vertexArrayPtr() { return !v.empty() ? &v[0] : 0; }
	const CVec3D *vertexArrayPtr() const { return !v.empty() ? &v[0] : 0; }

	/// Returns the <i>i</i>th vertex of this polyhedron.
	/** \param vertexIndex The vertex to get, in the range [0, numVertices()-1].
		\see numVertices(). */
	CVec3D vertex(int vertexIndex) const;

	/// Returns the <i>i</i>th edge of this polyhedron.
	/** Performance warning: Use this function only if you are interested in a single edge of this CPolyhedron.
		This function calls the edges() member function to receive a list of all the edges, so has
		a complexity of O(|V|log|V|), where |V| is the number of vertices in the polyhedron.
		\param edgeIndex The index of the edge to get, in the range [0, numEdges()-1].
		\see numEdges(), edges(), edgeIndices(). */
	CLineSegment edge(int edgeIndex) const;

	/// Returns all the (unique) edges of this polyhedron.
	/** Has complexity of O(|V|log|V|), where |V| is the number of vertices in the polyhedron.
		\todo Support this in linear time.
		\see numEdges(), edge(), edgeIndices(). */
	std::vector<CLineSegment> edges() const;

	std::vector<CPolygon> faces() const;

	/// Returns all the (unique) edges of this polyhedron, as indices to the polyhedron vertex array.
	/** Has complexity of O(|V|log|V|), where |V| is the number of vertices in the polyhedron.
		\todo Support this in linear time.
		\see numEdges(), edge(), edges(). */
	std::vector<std::pair<int, int> > edgeIndices() const;

	/// Returns a polygon representing the given face.
	/** The winding order of the polygon will be the same as in the input. The normal of the polygon
		points outwards from this polyhedron, i.e. towards the space that is not part of the polyhedron.
		This function constructs a new CPolygon object, so it has a time and space complexity of
		O(|V|), where |V| is the number of vertices in this polyhedron.
		\param faceIndex The index of the face to get, in the range [0, numFaces()-1].
		\see numFaces(), facePlane(). */
	CPolygon FacePolygon(int faceIndex) const;

	/// Returns the plane of the given polyhedron face.
	/** The normal of the plane points outwards from this polyhedron, i.e. towards the space that
		is not part of the polyhedron.
		This function assumes that the given face of the polyhedron is planar, as should be for all
		well-formed polyhedron.
		\param faceIndex The index of the face to get, in the range [0, numFaces()-1].
		\see numFaces(), FacePolygon(). */
	CPlane facePlane(int faceIndex) const;

	CVec3D faceNormal(int faceIndex) const;

	/// Returns the index of the vertex of this polyhedron that reaches farthest in the given direction.
	/** \param direction The direction vector to query for. This vector can be unnormalized.
		\return The supporting point of this polyhedron that reaches farthest in the given direction.
			The supporting point for a given direction is not necessarily unique, but this function
			will always return one of the vertices of this polyhedron.
		\see v, numVertices(), vertex(). */
	int extremeVertex(const CVec3D &direction) const;
	// CVec3D SupportingPoint(const CVec3D &dir) const;
	// bool IsSupportingPlane(const CPlane &plane) const;

	/// Computes an extreme point of this CPolyhedron in the given direction.
	/** An extreme point is a farthest point of this CPolyhedron in the given direction. Given a direction,
		this point is not necessarily unique.
		\param direction The direction vector of the direction to find the extreme point. This vector may
			be unnormalized, but may not be null.
		\return An extreme point of this CPolyhedron in the given direction. The returned point is always a
			corner point of this CPolyhedron.
		\see cornerPoint(). */
	CVec3D extremePoint(const CVec3D &direction) const;

	/// Projects this CPolyhedron onto the given 1D axis direction vector.
	/** This function collapses this CPolyhedron onto an 1D axis for the purposes of e.g. separate axis test computations.
		The function returns a 1D range [outMin, outMax] denoting the interval of the projection.
		\param direction The 1D axis to project to. This vector may be unnormalized, in which case the output
			of this function gets scaled by the length of this vector.
		\param outMin [out] Returns the minimum extent of this object along the projection axis.
		\param outMax [out] Returns the maximum extent of this object along the projection axis. */
	void projectToAxis(const CVec3D &direction, float &outMin, float &outMax) const;

	/// Returns the arithmetic mean of all the corner vertices.
	/** \bug This is not the proper centroid of the polyhedron! */
	/** \see surfaceArea(), volume(). */
	CVec3D centroid() const;

	/// Computes the total surface area of the faces of this polyhedron.
	/** \see centroid(), volume(). */
	float surfaceArea() const;

	/// Computes the internal volume of this polyhedron.
	/** \see centroid(), surfaceArea(). */
	float volume() const;

	/// Returns the smallest CAABBox that encloses this polyhedron.
	/// \todo add minimalEnclosingSphere() and minimalEnclosingOBB().
	CAABBox minimalEnclosingAABB() const;

#ifdef MATH_CONTAINERLIB_SUPPORT
	COBBox minimalEnclosingOBB() const;
#endif

	void mergeAdjacentPlanarFaces();

	/// Tests if the faces in this polyhedron refer to valid existing vertices.
	/** This function performs sanity checks on the face indices array.
		1) Each vertex index for each face must be in the range [0, numVertices()-1], i.e. refer to a vertex
		   that exists in the array.
		2) Each face must contain at least three vertices. If a face contains two or one vertex, it is
		   degenerate.
		3) Each face may refer to a vertex at most once.
		\return True if the abovementioned conditions hold. Note that this does not guarantee that the
			polyhedron is completely well-formed, but if false is returned, the polyhedron is definitely
			ill-formed.
		\see isNull(), IsClosed(), isConvex(). */
	bool faceIndicesValid() const;

	/// Flips the winding order of all faces in this polyhedron.
	void flipWindingOrder();

	/// Assuming that this polyhedron is convex, reorients all faces of this polyhedron
	/// so that each face plane has its normal pointing outwards. That is, the "negative" side of the
	/// polyhedron lies inside the polyhedron, and the positive side of the polyhedron is outside the convex
	/// shape.
	void orientNormalsOutsideConvex();

	/// Removes from the vertex array all vertices that are not referred to by any of the faces of this polyhedron.
	void removeRedundantVertices();

	/// Returns true if this polyhedron has 0 vertices and 0 faces.
	/** \see faceIndicesValid(), IsClosed(), isConvex(). */
	bool isNull() const { return v.empty() && f.empty(); }

	/// Returns true if this polyhedron is closed and does not have any gaps.
	/** \note This function performs a quick check, which might not be complete.
		\see faceIndicesValid(), IsClosed(), isConvex(). */
	bool IsClosed() const;

	// Returns true if this polyhedron forms a single connected solid volume.
//	bool IsConnected() const;

	/// Returns true if this polyhedron is convex.
	/** The running time is O(F*V) ~ O(V^2).
		\see faceIndicesValid(), IsClosed(), isConvex().*/
	bool isConvex() const;

	/// Returns true if the Euler formula (V + F - E == 2) holds for this CPolyhedron.
	/** \see numVertices(), numEdges(), numFaces(). */
	bool eulerFormulaHolds() const;

	/// Tests whether all the faces of this polyhedron are non-degenerate (have at least 3 vertices)
	/// and in case they have more than 3 vertices, tests that the faces are planar.
	bool facesAreNondegeneratePlanar(float epsilon = 1e-2f) const;

	/** 
	\brief Clips the line/ray/line segment specified by L(t) = ptA + t * dir, tFirst <= t <= tLast,
	 inside this <b>convex</b> polyhedron.
	 The implementation of this function is adapted from Christer Ericson's Real-time Collision Detection, p. 199.
	\param ptA The first endpoint of the line segment.
	\param dir The direction of the line segment. This member can be unnormalized.
	\param tFirst [in, out] As input, takes the parametric distance along the line segment which
			specifies the starting point of the line segment. As output, the starting point of the line segment
			after the clipping operation is returned here.
	\param tLast [in, out] As input, takes the parametric distance along the line segment which
			specifies the ending point of the line segment. As output, the endingpoint of the line segment
			after the clipping operation is returned here.
	\note To clip a line, pass in tFirst=-CMath::INFINITY_FLOAT, tLast=CMath::INFINITY_FLOAT. To clip a ray, pass in tFirst=0 and tLast = CMath::INFINITY_FLOAT.
			To clip a line segment, pass in tFirst=0, tLast=1, and an unnormalized dir = lineSegment.b-lineSegment.a.
	\return True if the outputted range [tFirst, tLast] did not become degenerate, and these two variables contain
			valid data. If false, the whole line segment was clipped away (it was completely outside this polyhedron).
	\see CMath::INFINITY_FLOAT. */
	bool clipLineSegmentToConvexPolyhedron(const CVec3D &ptA, const CVec3D &dir, float &tFirst, float &tLast) const;

	/// Tests if the given object is fully contained inside this polyhedron.
	/** This function treats this polyhedron as a non-convex object. If you know this polyhedron
		to be convex, you can use the faster containsConvex() function.
		\see containsConvex(), closestPoint(), closestPointConvex(), distance(), intersects(), intersectsConvex().
		\todo add contains(CCircle/Disc/CSphere/Capsule). */
	bool contains(const CVec3D &point) const;
	bool contains(const CLineSegment &lineSegment) const;
	bool contains(const CTriangle &triangle) const;
	bool contains(const CPolygon &polygon) const;
	bool contains(const CAABBox &aabb) const;
	bool contains(const COBBox &obb) const;
	bool contains(const CPolyhedron &polyhedron) const;

	/// Tests if the given face of this CPolyhedron contains the given point.
	bool faceContains(int faceIndex, const CVec3D &worldSpacePoint, float polygonThickness = 1e-3f) const;

	/** 
	\brief Tests if the given object is fully contained inside this <b>convex</b> polyhedron.
	 This function behaves exactly like contains(), except this version of the containment test
	assumes this polyhedron is convex, and uses a faster method of testing containment.
	\see contains(), closestPoint(), closestPointConvex(), distance(), intersects(), intersectsConvex().
	\todo add containsConvex(CPolygon/CAABBox/COBBox/Frustum/CPolyhedron/CCircle/Disc/CSphere/Capsule). */
	bool containsConvex(const CVec3D &point) const;
	bool containsConvex(const CLineSegment &lineSegment) const;
	bool containsConvex(const CTriangle &triangle) const;

	/** 
	\brief Computes the closest point on this polyhedron to the given object.
	 If the other object intersects this polyhedron, this function will return an arbitrary point inside
		the region of intersection.
	\param lineSegment The line segment to find the closest point to.
	\param lineSegmentPt [out] If specified, returns the closest point on the line segment to this
	polyhedron. This pointer may be null.
	\todo Make lineSegmentPt an out-reference instead of an out-pointer.
	\see contains(), containsConvex(), closestPointConvex(), distance(), intersects(), intersectsConvex().
	\todo add closestPoint(CLine/CRay/CPlane/CTriangle/CPolygon/CCircle/Disc/CAABBox/COBBox/CSphere/Capsule/Frustum/CPolyhedron). */
	CVec3D closestPoint(const CLineSegment &lineSegment, CVec3D *lineSegmentPt) const;
	CVec3D closestPoint(const CLineSegment &lineSegment) const;
	/** \param point The point to find the closest point to. */
	CVec3D closestPoint(const CVec3D &point) const;

	/// Returns the closest point on this <b>convex</b> polyhedron to the given point.
	/** This function behaves exactly like closestPoint(), except this version of the test assumes
		this polyhedron is convex, and uses a faster method of finding the closest point.
		\see contains(), containsConvex(), closestPoint(), distance(), intersects(), intersectsConvex().
		\todo add closestPointConvex(CLine/CLineSegment/CRay/CPlane/CTriangle/CPolygon/CCircle/Disc/CAABBox/COBBox/CSphere/Capsule/Frustum/CPolyhedron). */
	CVec3D closestPointConvex(const CVec3D &point) const;

	/// Returns the distance between this polyhedron and the given object.
	/** This function finds the nearest pair of points on this and the given object, and computes their distance.
		If the two objects intersect, or one object is contained inside the other, the returned distance is zero.
		\see contains(), containsConvex(), closestPoint(), closestPointConvex(), intersects(), intersectsConvex().
		\todo add distance(CLine/CLineSegment/CRay/CPlane/CTriangle/CPolygon/CCircle/Disc/CAABBox/COBBox/CSphere/Capsule/Frustum/CPolyhedron). */
	float distance(const CVec3D &point) const;

	/// Tests whether this polyhedron and the given object intersect.
	/** Both objects are treated as "solid", meaning that if one of the objects is fully contained inside
		another, this function still returns true. (e.g. in case a line segment is contained inside this polyhedron,
		or this polyhedron is contained inside a sphere, etc.)
		\return True if an intersection occurs or one of the objects is contained inside the other, false otherwise.
		\see contains(), containsConvex(), closestPoint(), closestPointConvex(), distance(), intersectsConvex().
		\todo add intersects(CCircle/Disc). */
	bool intersects(const CLineSegment &lineSegment) const;
	bool intersects(const CLine &line) const;
	bool intersects(const CRay &ray) const;
	bool intersects(const CPlane &plane) const;
	bool intersects(const CPolyhedron &polyhedron) const;
	bool intersects(const CAABBox &aabb) const;
	bool intersects(const COBBox &obb) const;
	bool intersects(const CTriangle &triangle) const;
	bool intersects(const CPolygon &polygon) const;
	bool intersects(const CSphere &sphere) const;

	/// Tests whether this <b>convex</b> polyhedron and the given object intersect.
	/** This function is exactly like intersects(), but this version assumes that this polyhedron is convex,
		and uses a faster method of testing the intersection.
		\return True if an intersection occurs or one of the objects is contained inside the other, false otherwise.
		\see contains(), containsConvex(), closestPoint(), closestPointConvex(), distance(), intersects().
		\todo add intersects(CCircle/Disc). */
	bool intersectsConvex(const CLine &line) const;
	bool intersectsConvex(const CRay &ray) const;
	bool intersectsConvex(const CLineSegment &lineSegment) const;

	void mergeConvex(const CVec3D &point);

	/**
	\brief Translates this CPolyhedron in world space.
	\param offset The amount of displacement to apply to this CPolyhedron, in world space coordinates.
	\see transform(). */
	void translate(const CVec3D &offset);

	/**
	\brief Applies a transformation to this CPolyhedron.
	 This function operates in-place.
	\see translate(), classes CMat3D, CMatJoint3x4, CMat4D, CQuaternion. */
	void transform(const CMat3D &transform);
	void transform(const CMatJoint3x4 &transform);
	void transform(const CMat4D &transform);
	void transform(const CQuaternion &transform);

	/// Creates a CPolyhedron object that represents the convex hull of the given point array.
	/// \todo This function is strongly WIP!
	static CPolyhedron convexHull(const CVec3D *pointArray, int numPoints);

	static CPolyhedron tetraHedron(const CVec3D &centerPos = CVec3D(0,0,0), float scale = 1.f, bool ccwIsFrontFacing = true);
	static CPolyhedron octaHedron(const CVec3D &centerPos = CVec3D(0,0,0), float scale = 1.f, bool ccwIsFrontFacing = true);
	static CPolyhedron hexaHedron(const CVec3D &centerPos = CVec3D(0,0,0), float scale = 1.f, bool ccwIsFrontFacing = true);
	static CPolyhedron icosaHedron(const CVec3D &centerPos = CVec3D(0,0,0), float scale = 1.f, bool ccwIsFrontFacing = true);
	static CPolyhedron dodecaHedron(const CVec3D &centerPos = CVec3D(0,0,0), float scale = 1.f, bool ccwIsFrontFacing = true);

	std::vector<CTriangle> triangulate() const;

#ifdef MATH_GRAPHICSENGINE_INTEROP
	void triangulate(VertexBuffer &vb, bool ccwIsFrontFacing) const;
	void toLineList(VertexBuffer &vb) const;
#endif
};

CPolyhedron operator *(const CMat3D &transform, const CPolyhedron &polyhedron);
CPolyhedron operator *(const CMatJoint3x4 &transform, const CPolyhedron &polyhedron);
CPolyhedron operator *(const CMat4D &transform, const CPolyhedron &polyhedron);
CPolyhedron operator *(const CQuaternion &transform, const CPolyhedron &polyhedron);

#ifdef MATH_QT_INTEROP
Q_DECLARE_METATYPE(CPolyhedron)
Q_DECLARE_METATYPE(CPolyhedron*)
#endif
/*
CPolyhedron operator *(const CMat3D &m, const CPolyhedron &s);
CPolyhedron operator *(const CMatJoint3x4 &m, const CPolyhedron &s);
CPolyhedron operator *(const CMat4D &m, const CPolyhedron &s);
CPolyhedron operator *(const CQuaternion &q, const CPolyhedron &s);
*/
// \todo add this
//std::ostream &operator <<(std::ostream &o, const CPolyhedron &polyhedron);


} //end GEO
}  //end SMF

#endif // __SMF_POLYHEDRON_
