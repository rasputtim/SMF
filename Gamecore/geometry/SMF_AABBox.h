


#ifndef __SMF_BOUNDS_H__
#define __SMF_BOUNDS_H__

#include "../SMF_Config.h"
#include "../math/SMF_Vector.h"
#include "SMF_Plane.h"

namespace SMF{
namespace MATH{
class CMat3D;
}
namespace GEO{
using namespace MATH;
class CPolyhedron;


/**
 * \class CAABBox
 *
 * \ingroup SMF_Geometric
 * \image html pics/aabb.png
 * \if pt_br
 * \brief (CAABBox) Axis ALigned Bounding Box (Caixa delimitadora alinhada por eixos)
	\n Esta estrutura de dados pode ser usada para representar limites grosseiras de objetos, em situações em que as computações detalhadas em nível de triângulos podem ser evitadas.
	\n Em sistemas de física, caixas delimitadoras são usadas ​​como um teste pré utilização eficiente para consultas de interseção geométrica.
	\n A parte 'axis-aligned' (alinhada por eixos) no nome significa que os eixos locais desta caixa delimitadora
	são restritas para se alinhar com os eixos do sistema de coordenadas espaciais.
	\n Isto faz com que os cálculos com CAABBox's sejam muito rápidos, uma vez que CAABBox's de não pode ser arbitrariamente orientada no espaço em relação umas as outras.

	\n \bSe você precisa de representar uma caixa no espaço 3D com orientação arbitrária, consulte a classe \ref COBBox \ref. \b
 * \elseif us_en
 * \brief 	(CAABBox) Axis Aligned Bounding Box

	This data structure can be used to represent coarse bounds of objects, in situations where detailed triangle-level computations can be avoided.
	In physics systems, bounding boxes are used as an efficient early-out test for geometry intersection queries.

	The 'axis-aligned' part in the name means that the local axes of this bounding box are restricted to align with the axes of the world space coordinate system.
	This makes computations involving CAABBox's very fast, since CAABBox's cannot be arbitrarily oriented in the space with respect to each other.

	If you need to represent a box in 3D space with arbitrary orientation, see the class \ref COBBox \ref.

 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CAABBox {
public:

	SMF_INLINE_FORCED static int numFaces() { return 6; }
	SMF_INLINE_FORCED static int numEdges() { return 12; }
	SMF_INLINE_FORCED static int numVertices() { return 8; }

	/**
	\brief Usando este cosntrutor os pontos min e max ficam indefinidos.
	       Você deve se lembrar de definí-los antes de usar a classe.
	\brief This means that the values of the members minPoint and maxPoint are undefined after creating a new CAABBox using this
        default constructor. Remember to assign to them before use.
    \see b[2]
	*/
	CAABBox();
	/// Constructs this CAABBox by specifying the minimum and maximum extending corners of the box.
    /** \see b[2]
	\brief Cria um objeto CAABBox especificando os vértices extensão mínimo e máximo.
	\param min Ponto que representa o vértice mínimo
	\param max Ponto que representa o vértice máximo
	*/
	explicit CAABBox( const CVec3D &mins, const CVec3D &maxs );
	/** \see b[2]
	\brief Cria um objeto CAABBox especificando os vértices extensão mínimo e máximo através de um único ponto.
	\param point Ponto que representa o vértice mínimo e o máximo
	*/
	explicit CAABBox( const CVec3D &point );
		/// Sets this CAABBox by specifying its center and size.
	/** \param center The center point of this CAABBox.
		\param size A vector that specifies the size of this CAABBox in x, y and z directions.
		\see setFrom(), FromCenterAndSize(). */
	void setFromCenterAndSize(const CVec3D &center, const CVec3D &size);

		/// Sets this CAABBox to enclose the given COBBox.
	/** This function computes the minimal axis-aligned bounding box for the given oriented bounding box. If the orientation
		of the COBBox is not aligned with the world axes, this conversion is not exact and loosens the volume of the bounding box.
		\param obb The oriented bounding box to convert into this CAABBox.
		\todo Implement setFrom(Polyhedron).
		\see SetCenter(), class COBBox. */
	void setFrom(const COBBox &obb);

	// Computes the minimal enclosing CAABBox of the given polyhedron.
	/* This function computes the smallest CAABBox (in terms of volume) that contains the given polyhedron, and stores
		the result in this structure.
		\note An CAABBox cannot generally exactly represent a polyhedron. Converting a polyhedron to an CAABBox loses some
		features of the polyhedron.
		\return If the given polyhedron is closed, this function succeeds and returns true. If the polyhedron is uncapped
			(has infinite volume), this function does not modify this data structure, but returns false. */
    //	bool setFrom(const Polyhedron &polyhedron);

	/// Sets this CAABBox to enclose the given sphere.
	/** This function computes the smallest possible CAABBox (in terms of volume) that contains the given sphere, and stores the result in this structure. */
	void setFrom(const CSphere &s);


	const CVec3D &	operator[]( const int index ) const;
	CVec3D &		operator[]( const int index );
	CAABBox			operator+( const CVec3D &t ) const;				// returns translated bounds
	CAABBox &		operator+=( const CVec3D &t );					// translate the bounds
	CAABBox			operator*( const CMat3D &r ) const;				// returns rotated bounds
	CAABBox &		operator*=( const CMat3D &r );					// rotate the bounds
	CAABBox			operator+( const CAABBox &a ) const;
	CAABBox &		operator+=( const CAABBox &a );
	CAABBox			operator-( const CAABBox &a ) const;
	CAABBox &		operator-=( const CAABBox &a );

	bool			compare( const CAABBox &a ) const;							// exact compare, no epsilon
	bool			compare( const CAABBox &a, const float epsilon ) const;	// compare with epsilon
	bool			operator==(	const CAABBox &a ) const;						// exact compare, no epsilon
	bool			operator!=(	const CAABBox &a ) const;						// exact compare, no epsilon

	void			clear();								// inside out bounds
	void			toZero();									// single point at origin

	CVec3D			getCenter() const;						// returns center of bounds
	CVec3D			centerPoint() const{return getCenter();};						// returns center of bounds

	float			getRadius() const;						// returns the radius relative to the bounds origin
	float			getRadius( const CVec3D &center ) const;		// returns the radius relative to the given center
	float			getVolume() const;						// returns the volume of the bounds

	sf_m128 &minPoint_SSE() { return *(sf_m128*)b[0].toFloatPtr(); }
	sf_m128 &maxPoint_SSE() { return *(sf_m128*)b[1].toFloatPtr(); }
	const sf_m128 &minPoint_SSE() const { return *(sf_m128*)b[0].toFloatPtr(); }
	const sf_m128 &maxPoint_SSE() const { return *(sf_m128*)b[1].toFloatPtr(); }

		/// Returns an edge of this CAABBox.
	/** \param edgeIndex The index of the edge line segment to get, in the range [0, 11].
		\todo Specify which index generates which edge.
		\see PointInside(), cornerPoint(), pointOnEdge(), faceCenterPoint(), facePoint(). */
	CLineSegment edge(int edgeIndex) const;

	/// Returns a corner point of this CAABBox.
	/** This function generates one of the eight corner points of this CAABBox.
		\param cornerIndex The index of the corner point to generate, in the range [0, 7].
			The points are returned in the order 0: ---, 1: --+, 2: -+-, 3: -++, 4: +--, 5: +-+, 6: ++-, 7: +++. (corresponding the XYZ axis directions).
		\todo Draw which index generates which corner point.
		\see PointInside(), edge(), pointOnEdge(), faceCenterPoint(), facePoint(), GetCornerPoints(). */
	CVec3D cornerPoint(int cornerIndex) const;

	/// Computes an extreme point of this CAABBox in the given direction.
	/** An extreme point is a farthest point of this CAABBox in the given direction. Given a direction,
		this point is not necessarily unique.
		\param direction The direction vector of the direction to find the extreme point. This vector may
			be unnormalized, but may not be null.
		\return An extreme point of this CAABBox in the given direction. The returned point is always a
			corner point of this CAABBox.
		\see cornerPoint(). */
	CVec3D extremePoint(const CVec3D &direction) const;

	/// Returns a point on an edge of this CAABBox.
	/** \param edgeIndex The index of the edge to generate a point to, in the range [0, 11]. \todo Document which index generates which one.
		\param u A normalized value between [0,1]. This specifies the relative distance of the point along the edge.
		\see PointInside(), cornerPoint(), cornerPoint(), faceCenterPoint(), facePoint(). */
	CVec3D pointOnEdge(int edgeIndex, float u) const;

	/// Returns the point at the center of the given face of this CAABBox.
	/** \param faceIndex The index of the CAABBox face to generate the point at. The valid range is [0, 5].
			This index corresponds to the planes in the order (-X, +X, -Y, +Y, -Z, +Z).
		\see PointInside(), cornerPoint(), pointOnEdge(), pointOnEdge(), facePoint(). */
	CVec3D faceCenterPoint(int faceIndex) const;

	/// Generates a point at the surface of the given face of this CAABBox.
	/** \param faceIndex The index of the CAABBox face to generate the point at. The valid range is [0, 5].
			This index corresponds to the planes in the order (-X, +X, -Y, +Y, -Z, +Z).
		\param u A normalized value between [0, 1].
		\param v A normalized value between [0, 1].
		\see PointInside(), cornerPoint(), pointOnEdge(), pointOnEdge(), faceCenterPoint(). */
	CVec3D facePoint(int faceIndex, float u, float v) const;


	/**
	\brief Returns the side lengths of this CAABBox in x, y and z directions.
	 The returned vector is equal to the diagonal vector of this CAABBox, i.e. it spans from the
		minimum corner of the CAABBox to the maximum corner of the CAABBox.
		\see halfSize(), diagonal(). */
	CVec3D getSize() const;

	/// [similarOverload: size]
	/** Returns size()/2.
		\see size(), halfDiagonal(). */
	CVec3D getHalfSize() const;

	CVec3D getHalfDiagonal() const{ return getHalfSize(); };
	// returns true if bounds are inside out
	bool			isCleared() const;
	// add the point, returns true if the bounds expanded
	bool			addPoint( const CVec3D &v );
	// add the bounds, returns true if the bounds expanded
	bool			addBounds( const CAABBox &a );
	// return bounds expanded in all directions with the given value
	CAABBox		expand( const float d ) const;
	// expand bounds in all directions with the given value
	CAABBox &	expandSelf( const float d );
	// return translated bounds
	CAABBox		translate( const CVec3D &translation ) const;
	// translate this bounds
	CAABBox &	translateSelf( const CVec3D &translation );
	// return rotated bounds
	CAABBox		rotate( const CMat3D &rotation ) const;
	// rotate this bounds
	CAABBox &	rotateSelf( const CMat3D &rotation );

	float			planeDistance( const CPlane &plane ) const;
	int				planeSide( const CPlane &plane, const float epsilon = ON_EPSILON ) const;

	bool			containsPoint( const CVec3D &p ) const;			// includes touching
		/// Tests if the given object is fully contained inside this CAABBox.
	/** This function returns true if the given object lies inside this CAABBox, and false otherwise.
		\note The comparison is performed using less-or-equal, so the faces of this CAABBox count as being inside, but
			due to float inaccuracies, this cannot generally be relied upon.
		\todo add contains(Circle/Disc/Sphere/Capsule).
		\see distance(), intersects(), closestPoint(). */
	bool contains(const CVec3D &point) const;
	bool contains(const CLineSegment &lineSegment) const;
	bool contains(const CAABBox &aabb) const;
	bool contains(const COBBox &obb) const;
	bool contains(const CSphere &sphere) const;
	bool contains(const CTriangle &triangle) const;
	bool contains(const CPolygon &polygon) const;
	bool contains(const CPolyhedron &polyhedron) const;


	bool			intersectsBounds( const CAABBox &a ) const;	// includes touching
	bool			lineIntersection( const CVec3D &start, const CVec3D &end ) const;
					// intersection point is start + dir * scale
	bool			rayIntersection( const CVec3D &start, const CVec3D &dir, float &scale ) const;

					// most tight bounds for the given transformed bounds
	void			fromTransformedBounds( const CAABBox &bounds, const CVec3D &origin, const CMat3D &axis );
					// most tight bounds for a point set
	void			fromPoints( const CVec3D *points, const int numPoints );
					// most tight bounds for a translation
	void			fromPointTranslation( const CVec3D &point, const CVec3D &translation );
	void			FromBoundsTranslation( const CAABBox &bounds, const CVec3D &origin, const CMat3D &axis, const CVec3D &translation );
					// most tight bounds for a rotation
	 void			fromPointRotation( const CVec3D &point, const CRotation &rotation );
	 void			fromBoundsRotation( const CAABBox &bounds, const CVec3D &origin, const CMat3D &axis, const CRotation &rotation );

	void			toPoints( CVec3D points[8] ) const;
	//CSphere		toSphere() const;
		/// Converts this AABB to an OBB.
	/** This function returns an OBB representation of this AABB. This conversion is exact, meaning that the returned
		OBB represents the same set of points than this AABB.
		\see class OBB. */
	COBBox toOBBox() const;
		/// Projects this CAABBox onto the given axis.
	/** \param axis The axis to project onto. This vector can be unnormalized.
		\param dMin [out] Returns the minimum extent of this CAABBox on the given axis.
		\param dMax [out] Returns the maximum extent of this CAABBox on the given axis. */
	void			projectToAxis( const CVec3D &dir, float &min, float &max ) const;
	void			projectToAxis( const CVec3D &origin, const CMat3D &axis, const CVec3D &dir, float &min, float &max ) const;

	int				getDimension() const;

	const float *	toFloatPtr() const;
	float *			toFloatPtr();
		/// Converts this CAABBox to a polyhedron.
	/** This function returns a polyhedron representation of this CAABBox. This conversion is exact, meaning that the returned
		polyhedron represents the same set of points than this CAABBox.
		\see class Polyhedron. */
	CPolyhedron toPolyhedron() const;

	//====distance===============
		/// Computes the closest point inside this CAABBox to the given point.
	/** If the target point lies inside this CAABBox, then that point is returned.
		\see distance(), contains(), intersects().
		\todo add closestPoint(Line/Ray/LineSegment/Plane/Triangle/Polygon/Circle/Disc/CAABBox/COBBox/Sphere/Capsule/Frustum/Polyhedron). */
	CVec3D closestPoint(const CVec3D &targetPoint) const;

		/// Computes the distance between this CAABBox and the given object.
	/** This function finds the nearest pair of points on this and the given object, and computes their distance.
		If the two objects intersect, or one object is contained inside the other, the returned distance is zero.
		\todo add CAABBox::distance(Line/Ray/LineSegment/Plane/Triangle/Polygon/Circle/Disc/CAABBox/COBBox/Capsule/Frustum/Polyhedron).
		\see contains(), intersects(), closestPoint(). */
	float distance(const CVec3D &point) const;
	float distance(const CSphere &sphere) const;


		/// Tests whether this CAABBox and the given object intersect.
	/** Both objects are treated as "solid", meaning that if one of the objects is fully contained inside
		another, this function still returns true. (e.g. in case a line segment is contained inside this CAABBox,
		or this CAABBox is contained inside a Sphere, etc.)
		\param ray The first parameter of this function specifies the other object to test against.
		\param dNear [out] If specified, receives the parametric distance along the line denoting where the
			line entered this CAABBox.
		\param dFar [out] If specified, receives the parametric distance along the line denoting where the
			line exited this CAABBox.
		\see contains(), distance(), closestPoint().
		\note If you do not need the intersection intervals, you should call the functions without these
			parameters in the function signature for optimal performance.
		\todo add intersects(Circle/Disc). */
	bool intersects(const CRay &ray, float &dNear, float &dFar) const;
	bool intersects(const CRay &ray) const;
	bool intersects(const CLine &line, float &dNear, float &dFar) const;
	bool intersects(const CLine &line) const;
	bool intersects(const CLineSegment &lineSegment, float &dNear, float &dFar) const;
	bool intersects(const CLineSegment &lineSegment) const;
	bool intersects(const CPlane &plane) const;
	bool intersects(const CAABBox &aabb) const;
	bool intersects(const COBBox &obb) const;
		// return intersection of this bounds with the given bounds
	CAABBox		Intersects( const CAABBox &a ) const;
	// intersect this bounds with the given bounds
	CAABBox &	intersectSelf( const CAABBox &a );

	/** For reference documentation on the Sphere-CAABBox intersection test, see Christer Ericson's Real-Time Collision Detection, p. 165. [groupSyntax]
		\param sphere The first parameter of this function specifies the other object to test against.
		\param closestPointOnAABB [out] Returns the closest point on this CAABBox to the given sphere. This pointer
			may be null. */
	bool intersects(const CSphere &sphere, CVec3D *closestPointOnAABB = 0) const;
	bool intersects(const CTriangle &triangle) const;
	bool intersects(const CPolygon &polygon) const;
	bool intersects(const CPolyhedron &polyhedron) const;
		/// Computes the intersection of a line, ray or line segment and an CAABBox.
	/** Based on "T. Kay, J. Kajiya. Ray Tracing Complex Scenes. SIGGRAPH 1986 vol 20, number 4. pp. 269-"
		http://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter3.htm
		\param linePos The starting position of the line.
		\param lineDir The direction of the line. This direction vector must be normalized!
		\param tNear [in, out] For the test, the input line is treated as a line segment. Pass in the signed distance
			from the line origin to the start of the line. For a Line-CAABBox test, -CMath::INFINITY_FLOAT is typically passed here.
			For a Ray-CAABBox test, 0.0f should be inputted. If intersection occurs, the signed distance from line origin
			to the line entry point in the CAABBox is returned here.
		\param tFar [in, out] Pass in the signed distance from the line origin to the end of the line. For Line-CAABBox and
			Ray-CAABBox tests, pass in CMath::INFINITY_FLOAT. For a LineSegment-CAABBox test, pass in the length of the line segment here.
			If intersection occurs, the signed distance from line origin to the line exit point in the CAABBox
			is returned here.
		\return True if an intersection occurs, false otherwise.
		\note This is a low level utility function. It may be more convenient to use one of the CAABBox::intersects()
			functions instead.
		\see intersects(). */
	bool intersectLineAABB(const CVec3D &linePos, const CVec3D &lineDir, float &tNear, float &tFar) const;

	/// Expands this AABB to enclose the given object.
	/** This function computes an AABB that encloses both this AABB and the specified object, and stores the resulting
		AABB into this.
		\note The generated AABB is not necessarily the optimal enclosing AABB for this AABB and the given object. */
	void enclose(const CVec3D &point);
	void enclose(const CLineSegment &lineSegment);
	void enclose(const CAABBox &aabb);
	void enclose(const COBBox &obb);
	void enclose(const CSphere &sphere);
	void enclose(const CTriangle &triangle);
	void enclose(const CPolygon &polygon);
	void enclose(const CPolyhedron &polyhedron);
	void enclose(const CVec3D *pointArray, int numPoints);
		/// Sets this structure to a degenerate AABB that does not have any volume.
	/** This function is useful for initializing the AABB to "null" before a loop of calls to enclose(),
		which incrementally expands the bounds of this AABB to enclose the given objects.
		\see enclose(). */
	void toNegativeInfinity();

	/// Applies a uniform scale to this AABB.
	/** This function scales this AABB structure in-place, using the given center point as the origin
		for the scaling operation.
		\param centerPoint Specifies the center of the scaling operation, in world space.
		\param scaleFactor The uniform scale factor to apply to each world space axis.
		\see translate(), transform(). */
	void scale(const CVec3D &centerPoint, float scaleFactor);

	/// Applies a non-uniform scale to this AABB.
	/** This function scales this AABB structure in-place, using the given center point as the origin
		for the scaling operation.
		\param centerPoint Specifies the center of the scaling operation, in world space.
		\param scaleFactor The non-uniform scale factors to apply to each world space axis.
		\see translate(), transform(). */
	void scale(const CVec3D &centerPoint, const CVec3D &scaleFactor);
	SMF_INLINE CVec3D minPoint()const {return b[0];}
	SMF_INLINE CVec3D maxPoint()const {return b[1];}


private:
	/// Specifies the minimum and maximum (b[0]=min,b[1]max) extent of this CAABBox in the world space x, y and z axes.
	///Especifica as extensõs mínima e máxima da CAABBox (b[0]=min,b[1]max) nos eixos espaço do mundo x, ye z.
	CVec3D			b[2];
};

extern CAABBox	bounds_zero;
extern CAABBox bounds_zeroOneCube;
extern CAABBox bounds_unitCube;

SMF_INLINE CAABBox::CAABBox() {
}

SMF_INLINE CAABBox::CAABBox( const CVec3D &mins, const CVec3D &maxs ) {
	b[0] = mins;
	b[1] = maxs;
}

SMF_INLINE CAABBox::CAABBox( const CVec3D &point ) {
	b[0] = point;
	b[1] = point;
}

SMF_INLINE const CVec3D &CAABBox::operator[]( const int index ) const {
	return b[index];
}

SMF_INLINE CVec3D &CAABBox::operator[]( const int index ) {
	return b[index];
}

SMF_INLINE CAABBox CAABBox::operator+( const CVec3D &t ) const {
	return CAABBox( b[0] + t, b[1] + t );
}

SMF_INLINE CAABBox &CAABBox::operator+=( const CVec3D &t ) {
	b[0] += t;
	b[1] += t;
	return *this;
}

SMF_INLINE CAABBox CAABBox::operator*( const CMat3D &r ) const {
	CAABBox bounds;
	bounds.fromTransformedBounds( *this, CVec3D::origin, r );
	return bounds;
}

SMF_INLINE CAABBox &CAABBox::operator*=( const CMat3D &r ) {
	this->fromTransformedBounds( *this, CVec3D::origin, r );
	return *this;
}

SMF_INLINE CAABBox CAABBox::operator+( const CAABBox &a ) const {
	CAABBox newBounds;
	newBounds = *this;
	newBounds.addBounds( a );
	return newBounds;
}

SMF_INLINE CAABBox &CAABBox::operator+=( const CAABBox &a ) {
	CAABBox::addBounds( a );
	return *this;
}

SMF_INLINE CAABBox CAABBox::operator-( const CAABBox &a ) const {
	SMF_ASSERT( b[1][0] - b[0][0] > a.b[1][0] - a.b[0][0] &&
				b[1][1] - b[0][1] > a.b[1][1] - a.b[0][1] &&
					b[1][2] - b[0][2] > a.b[1][2] - a.b[0][2] );
	return CAABBox( CVec3D( b[0][0] + a.b[1][0], b[0][1] + a.b[1][1], b[0][2] + a.b[1][2] ),
					CVec3D( b[1][0] + a.b[0][0], b[1][1] + a.b[0][1], b[1][2] + a.b[0][2] ) );
}

SMF_INLINE CAABBox &CAABBox::operator-=( const CAABBox &a ) {
	SMF_ASSERT( b[1][0] - b[0][0] > a.b[1][0] - a.b[0][0] &&
				b[1][1] - b[0][1] > a.b[1][1] - a.b[0][1] &&
					b[1][2] - b[0][2] > a.b[1][2] - a.b[0][2] );
	b[0] += a.b[1];
	b[1] += a.b[0];
	return *this;
}

SMF_INLINE bool CAABBox::compare( const CAABBox &a ) const {
	return ( b[0].compare( a.b[0] ) && b[1].compare( a.b[1] ) );
}

SMF_INLINE bool CAABBox::compare( const CAABBox &a, const float epsilon ) const {
	return ( b[0].compare( a.b[0], epsilon ) && b[1].compare( a.b[1], epsilon ) );
}

SMF_INLINE bool CAABBox::operator==( const CAABBox &a ) const {
	return compare( a );
}

SMF_INLINE bool CAABBox::operator!=( const CAABBox &a ) const {
	return !compare( a );
}

SMF_INLINE void CAABBox::clear() {
	b[0][0] = b[0][1] = b[0][2] = CMath::INFINITY_FLOAT;
	b[1][0] = b[1][1] = b[1][2] = -CMath::INFINITY_FLOAT;
}

SMF_INLINE void CAABBox::toZero() {
	b[0][0] = b[0][1] = b[0][2] =
	b[1][0] = b[1][1] = b[1][2] = 0;
}

SMF_INLINE CVec3D CAABBox::getCenter() const {
	return CVec3D( ( b[1][0] + b[0][0] ) * 0.5f, ( b[1][1] + b[0][1] ) * 0.5f, ( b[1][2] + b[0][2] ) * 0.5f );
}

SMF_INLINE float CAABBox::getVolume() const {
	if ( b[0][0] >= b[1][0] || b[0][1] >= b[1][1] || b[0][2] >= b[1][2] ) {
		return 0.0f;
	}
	return ( ( b[1][0] - b[0][0] ) * ( b[1][1] - b[0][1] ) * ( b[1][2] - b[0][2] ) );
}
SMF_INLINE CVec3D CAABBox::getSize() const
{
	return b[1] - b[0];
}

SMF_INLINE CVec3D CAABBox::getHalfSize() const
{
	return getSize() / 2.f;
}
SMF_INLINE bool CAABBox::isCleared() const {
	return b[0][0] > b[1][0];
}

SMF_INLINE bool CAABBox::addPoint( const CVec3D &v ) {
	bool expanded = false;
	if ( v[0] < b[0][0]) {
		b[0][0] = v[0];
		expanded = true;
	}
	if ( v[0] > b[1][0]) {
		b[1][0] = v[0];
		expanded = true;
	}
	if ( v[1] < b[0][1] ) {
		b[0][1] = v[1];
		expanded = true;
	}
	if ( v[1] > b[1][1]) {
		b[1][1] = v[1];
		expanded = true;
	}
	if ( v[2] < b[0][2] ) {
		b[0][2] = v[2];
		expanded = true;
	}
	if ( v[2] > b[1][2]) {
		b[1][2] = v[2];
		expanded = true;
	}
	return expanded;
}

SMF_INLINE bool CAABBox::addBounds( const CAABBox &a ) {
	bool expanded = false;
	if ( a.b[0][0] < b[0][0] ) {
		b[0][0] = a.b[0][0];
		expanded = true;
	}
	if ( a.b[0][1] < b[0][1] ) {
		b[0][1] = a.b[0][1];
		expanded = true;
	}
	if ( a.b[0][2] < b[0][2] ) {
		b[0][2] = a.b[0][2];
		expanded = true;
	}
	if ( a.b[1][0] > b[1][0] ) {
		b[1][0] = a.b[1][0];
		expanded = true;
	}
	if ( a.b[1][1] > b[1][1] ) {
		b[1][1] = a.b[1][1];
		expanded = true;
	}
	if ( a.b[1][2] > b[1][2] ) {
		b[1][2] = a.b[1][2];
		expanded = true;
	}
	return expanded;
}

SMF_INLINE CAABBox CAABBox::Intersects( const CAABBox &a ) const {
	CAABBox n;
	n.b[0][0] = ( a.b[0][0] > b[0][0] ) ? a.b[0][0] : b[0][0];
	n.b[0][1] = ( a.b[0][1] > b[0][1] ) ? a.b[0][1] : b[0][1];
	n.b[0][2] = ( a.b[0][2] > b[0][2] ) ? a.b[0][2] : b[0][2];
	n.b[1][0] = ( a.b[1][0] < b[1][0] ) ? a.b[1][0] : b[1][0];
	n.b[1][1] = ( a.b[1][1] < b[1][1] ) ? a.b[1][1] : b[1][1];
	n.b[1][2] = ( a.b[1][2] < b[1][2] ) ? a.b[1][2] : b[1][2];
	return n;
}

SMF_INLINE CAABBox &CAABBox::intersectSelf( const CAABBox &a ) {
	if ( a.b[0][0] > b[0][0] ) {
		b[0][0] = a.b[0][0];
	}
	if ( a.b[0][1] > b[0][1] ) {
		b[0][1] = a.b[0][1];
	}
	if ( a.b[0][2] > b[0][2] ) {
		b[0][2] = a.b[0][2];
	}
	if ( a.b[1][0] < b[1][0] ) {
		b[1][0] = a.b[1][0];
	}
	if ( a.b[1][1] < b[1][1] ) {
		b[1][1] = a.b[1][1];
	}
	if ( a.b[1][2] < b[1][2] ) {
		b[1][2] = a.b[1][2];
	}
	return *this;
}

SMF_INLINE CAABBox CAABBox::expand( const float d ) const {
	return CAABBox( CVec3D( b[0][0] - d, b[0][1] - d, b[0][2] - d ),
						CVec3D( b[1][0] + d, b[1][1] + d, b[1][2] + d ) );
}

SMF_INLINE CAABBox &CAABBox::expandSelf( const float d ) {
	b[0][0] -= d;
	b[0][1] -= d;
	b[0][2] -= d;
	b[1][0] += d;
	b[1][1] += d;
	b[1][2] += d;
	return *this;
}

SMF_INLINE CAABBox CAABBox::translate( const CVec3D &translation ) const {
	return CAABBox( b[0] + translation, b[1] + translation );
}

SMF_INLINE CAABBox &CAABBox::translateSelf( const CVec3D &translation ) {
	b[0] += translation;
	b[1] += translation;
	return *this;
}

SMF_INLINE CAABBox CAABBox::rotate( const CMat3D &rotation ) const {
	CAABBox bounds;
	bounds.fromTransformedBounds( *this, CVec3D::origin, rotation );
	return bounds;
}

SMF_INLINE CAABBox &CAABBox::rotateSelf( const CMat3D &rotation ) {
	fromTransformedBounds( *this, CVec3D::origin, rotation );
	return *this;
}

SMF_INLINE bool CAABBox::containsPoint( const CVec3D &p ) const {
	if ( p[0] < b[0][0] || p[1] < b[0][1] || p[2] < b[0][2]
		|| p[0] > b[1][0] || p[1] > b[1][1] || p[2] > b[1][2] ) {
		return false;
	}
	return true;
}

SMF_INLINE bool CAABBox::intersectsBounds( const CAABBox &a ) const {
	if ( a.b[1][0] < b[0][0] || a.b[1][1] < b[0][1] || a.b[1][2] < b[0][2]
		|| a.b[0][0] > b[1][0] || a.b[0][1] > b[1][1] || a.b[0][2] > b[1][2] ) {
		return false;
	}
	return true;
}
#if 0
SMF_INLINE CSphere CAABBox::toSphere() const {
	CSphere sphere;
	sphere.setOrigin( ( b[0] + b[1] ) * 0.5f );
	sphere.SetRadius( ( b[1] - sphere.getOrigin() ).getLenght() );
	return sphere;
}
#endif
SMF_INLINE void CAABBox::projectToAxis( const CVec3D &dir, float &min, float &max ) const {
	float d1, d2;
	CVec3D center, extents;

	center = ( b[0] + b[1] ) * 0.5f;
	extents = b[1] - center;

	d1 = dir * center;
	d2 = CMath::fabs( extents[0] * dir[0] ) +
			CMath::fabs( extents[1] * dir[1] ) +
				CMath::fabs( extents[2] * dir[2] );

	min = d1 - d2;
	max = d1 + d2;
}

SMF_INLINE void CAABBox::projectToAxis( const CVec3D &origin, const CMat3D &axis, const CVec3D &dir, float &min, float &max ) const {
	float d1, d2;
	CVec3D center, extents;

	center = ( b[0] + b[1] ) * 0.5f;
	extents = b[1] - center;
	center = origin + center * axis;

	d1 = dir * center;
	d2 = CMath::fabs( extents[0] * ( dir * axis.Row(0) ) ) +
			CMath::fabs( extents[1] * ( dir * axis.Row(1) ) ) +
				CMath::fabs( extents[2] * ( dir * axis.Row(2) ) );

	min = d1 - d2;
	max = d1 + d2;
}

SMF_INLINE int CAABBox::getDimension() const {
	return 6;
}

SMF_INLINE const float *CAABBox::toFloatPtr() const {
	return &b[0].x;
}

SMF_INLINE float *CAABBox::toFloatPtr() {
	return &b[0].x;
}

} // end GEO
} //end SMF

#endif /* !__SMF_BOUNDS_H__ */
