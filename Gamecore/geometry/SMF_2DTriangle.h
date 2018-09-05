#ifndef __SMF_TRIANGLE_2D_
#define __SMF_TRIANGLE_2D_

#include "../SMF_Config.h"
#include "../math/all.h"
#include "all.h"

namespace SMF {
using namespace MATH;
using namespace Util;
namespace GEO{


/**
 * \class CTriangle2D
 *
 * \ingroup SMF_Geometric
 * \image html pics/tringle.png
 * \if pt_br
 * \brief   Representa um triângulo em espaço cartesiano 2D
 * \note Armazena os três pontos vértices a,b,c que definem um triângulo 
 * \elseif us_en
 * \brief 	Specifies a triangle through three points in 3D space.
 * This class stores three member vertices a, b and c to specify the triangle.
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CTriangle2D
{
public:
	enum triangeSide{
		side1 =0,
		side2 =1,
		side3 =2
	};

	/// The first triangle endpoint.
	CPoint2D a;
	/// The second triangle endpoint.
	/** [similarOverload: a] */
	CPoint2D b;
	/// The third triangle endpoint.
	/** [similarOverload: a]
		\note The order in which the vertices are stored in this data structure is important. The triangles
			(a,b,c) and (a,c,b) are otherwise equivalent, but their plane normals point to the opposite directions.
		\see planeCCW(), planeCW(). */
	CPoint2D c;

	/// The default constructor does not initialize any members of this class.
	/** This means that the values of the members a, b and c are undefined after creating a new CTriangle2D using this
		default constructor. Remember to assign to them before use.
		\see a, b, c. */
	CTriangle2D() {}

	/// Constructs a triangle from three given endpoints.
	/** The normal of the plane will be constructed to point towards the halfspace where
		the vertices a, b and c wind in counter-clockwise order.
		\see a, b, c. */
	CTriangle2D(const CPoint2D &a, const CPoint2D &b, const CPoint2D &c);

	SMF_INLINE_FORCED static int numFaces() { return 1; }
	SMF_INLINE_FORCED static int numEdges() { return 3; }
	SMF_INLINE_FORCED static int numVertices() { return 3; }


	  	  /** \brief Get the ith point
	  \warning if i not [0,2] return Infinite Point
	  **/
  SMF_INLINE CPoint2D operator[](int i) { if (i==0) return a;else if(i==1) return b;else if(i==2) return c; else return CPoint2D(CMath::NAN_FLOAT,CMath::NAN_FLOAT);}

  	  /** \brief Get the ith point
	  \warning if i not [0,2] return Infinite Point
	  **/
   SMF_INLINE CPoint2D const operator[](int i) const { if (i==0) return a;else if(i==1) return b;else if(i==2) return c; else return CPoint2D(CMath::NAN_FLOAT,CMath::NAN_FLOAT);}


	/// Translates this CTriangle2D in world space.
	/** \param offset The amount of displacement to apply to this CTriangle2D, in world space coordinates.
		\see transform(). */
	void translate(const CVec2D &offset);


	/// Expresses the given point in barycentric (u,v,w) coordinates.
	/** \note There are two different conventions for representing barycentric coordinates. One uses
			a (u,v,w) triplet with the equation pt == u*a + v*b + w*c, and the other uses a (u,v) pair
			with the equation pt == a + u*(b-a) + v*(c-a). These two are equivalent. Use the mappings
			(u,v) -> (1-u-v, u, v) and (u,v,w)->(v,w) to convert between these two representations.
		\param point The point of the vector space to express in barycentric coordinates. This point should
			lie in the plane formed by this triangle.
		\return The factors (u,v,w) that satisfy the weighted sum equation point == u*a + v*b + w*c.
		\see barycentricUV(), barycentricInsideTriangleboundingAABB(), Point(), http://mathworld.wolfram.com/BarycentricCoordinates.html
		\see http://en.wikipedia.org/wiki/Barycentric_coordinate_system **/
	CVec3D barycentricUVW(const CPoint2D &point) const;

	/// Expresses the given point in barycentric (u,v) coordinates.
	/** \note There are two different conventions for representing barycentric coordinates. One uses
			a (u,v,w) triplet with the equation pt == u*a + v*b + w*c, and the other uses a (u,v) pair
			with the equation pt == a + u*(b-a) + v*(c-a). These two are equivalent. Use the mappings
			(u,v) -> (1-u-v, u, v) and (u,v,w)->(v,w) to convert between these two representations.
		\param point The point to express in barycentric coordinates. This point should lie in the plane
			formed by this triangle.
		\return The factors (u,v) that satisfy the weighted sum equation point == a + u*(b-a) + v*(c-a).
		\see barycentricUVW(), barycentricInsideTriangleboundingAABB(), Point(). */
	CVec2D barycentricUV(const CPoint2D &point) const;

	/// Tests if the given barycentric UVW coordinates lie inside a triangle.
	/** A barycentric UVW coordinate represents a point inside a triangle if
		  a) 0 <= u,v,w <= 1 and
		  b) u+v+w == 1.
		\param uvw The barycentric vector containing the barycentric (u,v,w) coordinates.
		\return True if the given coordinates lie inside a triangle.
		\see barycentricUV(), barycentricUVW(), Point(). */
	static bool barycentricInsideTriangleboundingAABB(const CVec3D &uvw);

	/** These functions are an alternate form of Point(u,v,w) for the case when the barycentric coordinates are
		represented as a (u,v) pair and not as a (u,v,w) triplet. This function is provided for convenience
		and effectively just computes Point(1-u-v, u, v).
		\param uv The barycentric UV coordinates. If 0 <= u,v <= 1 and u+v <= 1, then the returned point lies inside
			this triangle.
		\return a + (b-a)*u + (c-a)*v.
		\see barycentricUV(), barycentricUVW(), barycentricInsideTriangleboundingAABB(). */
 	CVec2D Point(const CVec2D &uv) const;
	CVec2D Point(float u, float v) const;


	/// Computes the surface area of this triangle.
	/** \return The surface area of this triangle.
		\see perimeter(), area2D(), signedArea(). */
	float area() const;

	/// Computes the total edge length of this triangle.
	/** \return |a-b| + |b-c| + |c-a|.
		\see area(), edge(). */
	float perimeter() const;

	/// Returns a pointer to the vertices of this triangle. The array contains three elements.
	CPoint2D *vertexArrayPtr() { return &a; }
	const CPoint2D *vertexArrayPtr() const { return &a; }

	/// Returns a vertex of this triangle.
	/** \param i The vertex of this triangle to get: 0, 1 or 2.
		\return vertex(0) returns the point a, vertex(1) returns the point b, and vertex(2) returns the point c.
		\note If an index outside [0, 2] is passed, an SMF_ASSERT() failure occurs and CVec2D(NaN) is returned.
		\see edge(). */
	CPoint2D vertex(int i) const;

	/// Returns an edge of this triangle.
	/** \param i The index of the edge to generate: 0, 1 or 2.
		\return A CLineSegment2D representing the given edge of this triangle. edge(0) returns CLineSegment2D(a,b), edge(1)
			returns CLineSegment2D(b,c) and edge(2) returns CLineSegment2D(c,a).
		\note If an index outside [0, 2] is passed, an SMF_ASSERT() failure occurs and CLineSegment2D(NaN, NaN) is returned.
		\see vertex(). */
	CLineSegment2D edge(int i) const;
	/// Returns the midPoint for an edge of this triangle.
	/** \param i The index of the edge to generate the midPoint: 0, 1 or 2.
		\return A CPoint2D representing the MidPoint for the given edge of this triangle. midPointEdge(0) returns CPoint2D(a,b)/2, midPointEdge(1)
			returns CPoint2D(b,c)/2 and midPointEdge(2) returns CPoint2D(c,a)/2.
		\note If an index outside [0, 2] is passed, an SMF_ASSERT() failure occurs and CPoint2D(NaN, NaN) is returned.
		\see vertex(). */
	CPoint2D midPointEdge(int i) const;


	/**
	\brief return the median segment for the side of the triangle
	\note There are some fascinating properties of the medians of a triangle:

    The fact that the three medians always meet at a single point is interesting in its own right
    Each median divides the triangle into two smaller triangles which have the same area
    The centroid (point where they meet) is the center of gravity of the triangle
    The three medians divide the triangle into 6 smaller triangles that all have the same area, even though they may have different shapes.
	**/
	CLineSegment2D getMedian(int side)const;
	/// Computes the center of mass of this triangle.
	/** \return The point (a+b+c)/3. */
	CVec2D centroid() const;
	/// Identical to centerPoint(), but provided to enable common signature with CTriangle2D, CAABBox2D and COBBox to allow them to be used
	/// in template classes.
	CVec2D getCenter() const { return centroid(); }
	/**
	\brief return the incenter point of the triangle represented by the three points
	\note The incenter is the center of the triangle's incircle, the largest circle that will fit inside the triangle and touch all three sides.
	\note Given the coordinates of the three vertices of a triangle ABC,
	the coordinates of the incenter O are
	\f$  Ox=((a*Ax)+(b*Bx)+(c*Cx)) / p   Oy= ((a*Ay)+(b*By)+(c*Cy)) / p  \f$
	where:
	Ax and Ay 	are the x and y coordinates of the point A etc..
	a, b and c 	are the side lengths opposite vertex A, B and C
	p 	is perimeter of the triangle (a+b+c)
	**/
	static CPoint2D incenter(const CPoint2D &p1, const CPoint2D &p2,const CPoint2D &p3);
	/**
	\return the incenter of this triangle
	**/
	CPoint2D getIncenter()const;
	/**
	\note One of several centers the triangle can have, the circumcenter is the point where the perpendicular
	bisectors of a triangle intersect. The circumcenter is also the center of the triangle's circumcircle - the
	circle that passes through all three of the triangle's vertices.

	**/
	static CPoint2D circumcenter(const CPoint2D &p1, const CPoint2D &p2, const CPoint2D &p3);
	CPoint2D getCircumcenter()const;

	static CPoint2D orthoCenter(const CPoint2D &p1, const CPoint2D &p2,const CPoint2D &p3);
	CPoint2D getOrthocenter()const;
	/// Computes an extreme point of this CTriangle2D in the given direction.
	/** An extreme point is a farthest point of this CTriangle2D in the given direction. Given a direction,
		this point is not necessarily unique.
		\param direction The direction vector of the direction to find the extreme point. This vector may
			be unnormalized, but may not be null.
		\return An extreme point of this CTriangle2D in the given direction. The returned point is always a
			vertex of this CTriangle2D.
		\see vertex(). */
	CPoint2D extremePoint(const CVec2D &direction) const;

	/// Returns a CPolygon2D representation of this CTriangle2D.
	/** The returned polygon is identical to this CTriangle2D. It has three vertices a, b and c which wind in the same
		direction than in this triangle.
	\see class CPolygon2D, toPolyhedron(). */
	CPolygon2D toPolygon() const;

	/// Returns a CPolyhedron representation of this CTriangle2D.
	/** The generated polyhedron will be closed and has two triangular faces, and three vertices (a, b and c).
		The two faces share the same vertices, but in opposite winding orders. This creates a polyhedron with zero
		volume and the surface area twice of this CTriangle2D.
		\see class CPolyhedron, toPolygon(). */
	//CPolyhedron toPolyhedron() const;

	/// Returns the tight CAABBox2D that encloses this CTriangle2D.
	CAABBox2D boundingAABB() const;

	/** Computes the surface area of the given 2D triangle.
	 \see area(), signedArea(). */
	static float area(const CVec2D &p1, const CVec2D &p2, const CVec2D &p3);

	/** \brief Relates the barycentric coordinate of the given point to the surface area of triangle abc.
	 This function computes the ratio of the signed area of the triangle (point, b, c) to the signed area of
		the triangle (a, b, c). This is the same as computing the barycentric u-coordinate of the given point
		on the triangle (a, b, c).
		\see area(), area2D(), barycentricUVW(). */
	static float signedArea(const CVec2D &point, const CVec2D &a, const CVec2D &b, const CVec2D &c);
	/// Tests if this CTriangle2D is finite.
	/** A triangle is <b><i>finite</i></b> if its vertices a, b and c do not contain floating-point NaNs or +/-infs
		in them.
		\return True if each coordinate of each vertex of this triangle has a finite floating-point value.
		\see a, b, c, isDegenerate(), ::isFinite(), isInf(), IsNan(), isFinite(), inf, negInf, nan, CVec2D::nan, CVec2D::inf. */
	bool isFinite() const;

	/// Returns true if this triangle is degenerate.
	/** A triangle is <b><i>degenerate</i></b> if it is not finite, or if the surface area of the triangle is
		close to zero.
		\param epsilon The threshold to test against. If two vertices of this triangle are closer than this, the
		triangle is considered degenerate.
		\see a, b, c, isFinite(). */
	bool isDegenerate(float epsilon =CMath::EPSILON_SuperLow) const;

	/// Returns true if the triangle defined by the three given points is degenerate.
	static bool isDegenerate(const CPoint2D &p1, const CPoint2D &p2, const CPoint2D &p3, float epsilon =CMath::EPSILON_SuperLow);

	/// Tests if the given object is fully contained inside this triangle.
	/** \see distance(), intersects(), closestPoint().
		\todo add CTriangle2D::contains(CCircle) and CTriangle2D::contains(Disc). */
	bool contains(const CPoint2D &point) const;
	bool contains(const CLineSegment2D &lineSegment) const;
	bool contains(const CTriangle2D &triangle) const;

	/// Computes the distance between this triangle and the given object.
	/** This function finds the nearest pair of points on this and the given object, and computes their distance.
		If the two objects intersect, or one object is contained inside the other, the returned distance is zero.
		\todo add CTriangle2D::distance(CLine2D/CRay2D/CLineSegment2D/CPlane/CTriangle2D/CPolygon2D/CCircle/Disc/CAABBox2D/COBBox/Capsule/Frustum/CPolyhedron).
		\see contains(), intersects(), closestPoint(). */
	float distance(const CPoint2D &point) const;
	float distance(const CCircle2D &sphere) const;

	/// Tests whether this triangle and the given object intersect.
	/** Both objects are treated as "solid", meaning that if one of the objects is fully contained inside
		another, this function still returns true. (e.g. in case a line segment is contained inside this triangle,
		or this triangle is contained inside a sphere, etc.)
		\param intersectionPoint [out] If specified, receives the actual point of intersection. This pointer
			may be null.
		\return True if an intersection occurs or one of the objects is contained inside the other, false otherwise.
		\see contains(), distance(), closestPoint(), CLineSegment2D::getPoint(). */
	bool intersection(const CLineSegment2D &Segment, int *ICnt=NULL, CPoint2D *I1=NULL, CPoint2D *I2=NULL) const;
	bool intersection(const CLine2D &line,  int *ICnt=NULL, CPoint2D *I1=NULL, CPoint2D *I2=NULL) const;
	bool intersection(const CRay2D &ray,  int *ICnt=NULL, CPoint2D *I1=NULL, CPoint2D *I2=NULL) const;
//	bool intersection(const CAABBox2D &aabb,  int *ICnt=NULL, CPoint2D *I1=NULL, CPoint2D *I2=NULL) const;

	/** \param closestPointOnTriangle [out] If specified, receives the point of intersection between the CCircle2D
			and this CTriangle2D. Even if no intersection occurred, this parameter will receive the closest point on
			the CTriangle2D to the CCircle2D. This pointer may be null. */
	bool intersects(const CCircle2D &sphere, CPoint2D *closestPointOnTriangle) const;
	/**
	\brief fast method just to know if there are intersection
	**/
    static bool intersects( const CTriangle2D &tri, const CLineSegment2D &line);

	bool intersects(const CCircle2D &circle) const;
	bool intersects(const CLineSegment2D &segment) const;
	bool intersects(const CAABBox2D &aabb) const;
	bool intersects(const CLine2D &line) const;
	bool intersects(const CPolygon2D &polygon) const;
	bool intersects(const CTriangle2D &Triangle2)const;
	//bool intersects(const COBBox &obb) const;
	//bool intersects(const CPolyhedron &polyhedron) const;


	/// Projects this CTriangle2D onto the given axis.
	/** This function is used in SAT tests (separate axis theorem) to check the interval this triangle
		lies in on an 1D line.
		\param axis The axis to project onto. This vector can be unnormalized.
		\param dMin [out] Returns the minimum extent of this triangle on the given axis.
		\param dMax [out] Returns the maximum extent of this triangle on the given axis. */
	void projectToAxis(const CVec2D &axis, float &dMin, float &dMax) const;

	/// Computes the closest point on this triangle to the given object.
	/** If the other object intersects this triangle, this function will return an arbitrary point inside
		the region of intersection.
		\note algorith on pg 142 - Real time colision detection from Christer Ericson
		\see contains(), distance(), intersects(), closestPointToTriangleEdgesignedArea(). */
	CPoint2D closestPoint(const CPoint2D &point) const;
	CPoint2D closestPoint_m2(const CPoint2D &point)const;

	/// Generates a random point inside this CTriangle2D.
	/** The points are distributed uniformly.
		The implementation of this function is based on Graphics Gems 1, p. 25:
		"1.5 Generating random points in triangles. Method 2." The Method 1 presented in the book
		uses a sqrt() instead of the if().
		\param rng A pre-seeded random number generator object that is to be used by this function to generate random values.
		\see class CRandomLCG, randomPointOnEdge(), randomVertex(), Point(). */
	CPoint2D randomPointInside(CRandomLCG &rng) const;

	/// Chooses a corner vertex of this CTriangle2D at random.
	/** This function returns one of the vertices {a, b, c} at uniformly random.
		\param rng A pre-seeded random number generator object that is to be used by this function to generate random values.
		\see class CRandomLCG, randomPointInside(), randomPointOnEdge(), vertex(). */
	CPoint2D randomVertex(CRandomLCG &rng) const;

	/// Generates a random point on the edge of this CTriangle2D.
	/** The points are distributed uniformly.
		This function requires that this triangle is not degenerate. If it is, an SMF_ASSERT() error will be printed,
		and the return value will be undefined.
		\param rng A pre-seeded random number generator object that is to be used by this function to generate random values.
		\see class CRandomLCG, randomPointInside(), randomVertex(), edge(), class CLineSegment2D, isDegenerate(). */
	CPoint2D randomPointOnEdge(CRandomLCG &rng) const;

	/**
	\brief return the circle circunscrit on this triangle
	**/
	CCircle2D circunCircle()const;
	/**
	\brief return the circle inscrit on this triangle
	\note The radius of the incircle. The radius is given by the formula:

	\f$         2a / p                      \f$
	where:
	a is the area of the triangle.
	p is the perimeter of the triangle, the sum of its sides.
	**/
	CCircle2D inCircle()const;
	/// Returns a human-readable representation of this CLine2D. Most useful for debugging purposes.
	std::string toString() const;

#ifdef MATH_QT_INTEROP
	operator QString() const { return toString(); }
	QString toString() const { return QString::fromStdString(toString()); }
#endif
};


#ifdef MATH_QT_INTEROP
Q_DECLARE_METATYPE(CTriangle2D)
Q_DECLARE_METATYPE(CTriangle2D*)
#endif


std::ostream &operator <<(std::ostream &o, const CTriangle2D &triangle);



} //end GEO
}  //end SMF

#endif // __SMF_TRIANGLE_2D_
