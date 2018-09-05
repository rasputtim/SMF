

#ifndef __SMF_BOX_H__
#define __SMF_BOX_H__

#include "../SMF_Config.h"
#include "SMF_AABBox.h"
#include "SMF_Sphere.h"
#include "../math/SMF_Vector.h"
namespace SMF{
using namespace MATH;
namespace MATH {

}
namespace GEO {
class CSphere;




/**
 * \class COBBox
 * \image html pics\obb.png
 * \ingroup SMF_Geometric
 *
 * \if pt_br
 * \brief (COBBox)  Oriented Bounding Box (Caixa delimitadora Orientada)
    \n   Uma COBBox é uma Caixa delimitadora retangular numa oritentação arbitrária num espaco 3D.
	\see http://www.cs.unc.edu/~walk/papers/gottscha/sig96.pdf ,
	\see http://www.geometrictools.com/Documentation/DynamicCollisionDetection.pdf
	\see http://gamedev.stackexchange.com/questions/25397/obb-vs-obb-collision-detection
	\see http://www.flipcode.com/archives/2D_OBB_Intersection.shtml

 * \elseif us_en
 * \brief Oriented Bounding Box
    \n An COBBox is a rectangular bounding box at an arbitrary orientation in 3D-space
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API COBBox {
public:
	SMF_INLINE_FORCED static int numFaces() { return 6; }
	SMF_INLINE_FORCED static int numEdges() { return 12; }
	SMF_INLINE_FORCED static int numVertices() { return 8; }

					COBBox();   
					explicit COBBox( const CVec3D &center, const CVec3D &extents, const CMat3D &axis );
					explicit COBBox( const CVec3D &point );
					explicit COBBox( const CAABBox &bounds );
					explicit COBBox( const CAABBox &bounds, const CVec3D &origin, const CMat3D &axis );
	/// returns translated box
	COBBox			operator+( const CVec3D &t ) const;				
	/// translate the box
	COBBox &			operator+=( const CVec3D &t );					
	/// returns rotated box
	COBBox			operator*( const CMat3D &r ) const;				
	/// rotate the box
	COBBox &			operator*=( const CMat3D &r );					
	COBBox			operator+( const COBBox &a ) const;
	COBBox &			operator+=( const COBBox &a );
	COBBox			operator-( const COBBox &a ) const;
	COBBox &			operator-=( const COBBox &a );
	/// exact compare, no epsilon
	bool			compare( const COBBox &a ) const;						
	/// compare with epsilon
	bool			compare( const COBBox &a, const float epsilon ) const;	
	/// exact compare, no epsilon
	bool			operator==(	const COBBox &a ) const;						
	/// exact compare, no epsilon
	bool			operator!=(	const COBBox &a ) const;						
	/// inside out box
	void			clear();									
	/// single point at origin
	void			toZero();									
	/** 
	\brief Returns the side lengths of this COBBox in its local x, y and z directions.
	\return 2*extents. */
	CVec3D size() const;

	/** 
	\brief Returns the half-side lengths of this COBBox in its local x, y and z directions. [similarOverload: size]
	\return extents.
	\see extents, size(), halfSize(). */
	CVec3D halfSize() const;

	
	/** 
	\brief Returns a diagonal vector of this COBBox.
	This vector runs from one corner of the COBBox from the opposite corner.
	\note A box has four diagonals. This function returns the direction vector from the -X-Y-Z corner to
			the +X+Y+Z corner of the COBBox, in the global space of this COBBox. */
	CVec3D diagonal() const;
	
	/**
	\brief Returns diagonal()/2. [similarOverload: diagonal].
	\return A direction vector from the center of this COBBox to the +X+Y+Z corner of this COBBox, in global space.
	\see size(), halfSize(). */
	CVec3D halfDiagonal() const;

	 
	/** 
	\brief Computes the transformation matrix that maps from the global (world) space of this COBBox to the local space of this COBBox.
		In local space, the center of the COBBox lies at (extents.x,extents.y,extents.z), and the COBBox is aligned along the cardinal axes, i.e. is an CAABBox.
		The local +X vector runs in the direction specified by axis[0], the +Y direction is specified by axis[1], and +Z by axis[2].
		The size of the COBBox is 2*extents.
		In global (world) space, the center of the COBBox lies at the point given by the pos member variable.
	\return This global (world) to local space transform can be represented using a CMat3D matrix. This function computes
			and returns the matrix that maps from the world space of this COBBox to the local space of this COBBox.
	\see pos, extents, axis. */
	CMatJoint3x4 worldToLocal() const;

	/**
	\brief Computes the transformation matrix that maps from the local space of this COBBox to the global (world) space of this COBBox. [similarOverload: worldToLocal]
	 This mapping is the inverse of the transform computed by worldToLocal().
	\return A matrix that transforms from the local space of this COBBox to the global (world) space of this COBBox.
	\see pos, extents, axis. */
	CMatJoint3x4 localToWorld() const;
	// returns center of the box
	const CVec3D &	getCenter() const;						
	// returns extents of the box
	const CVec3D &	getExtents() const;						
	// returns the axis of the box
	const CMat3D &	getAxis() const;							
	/// returns the volume of the box
	float			getVolume() const;		
	/// returns true if box are inside out
	bool			isCleared() const;						
	/// add the point, returns true if the box expanded
	bool			addPoint( const CVec3D &v );					
	/// add the box, returns true if the box expanded
	bool			AddBox( const COBBox &a );						
	/// return box expanded in all directions with the given value
	COBBox			expand( const float d ) const;					
	/// expand box in all directions with the given value
	COBBox &			expandSelf( const float d );					
	/// return translated box
	COBBox			translate( const CVec3D &translation ) const;	
	/// translate this box
	COBBox &			translateSelf( const CVec3D &translation );		
	/// return rotated box
	COBBox			rotate( const CMat3D &rotation ) const;			
	/// rotate this box
	COBBox &			rotateSelf( const CMat3D &rotation );			

	float			planeDistance( const CPlane &plane ) const;
	int				planeSide( const CPlane &plane, const float epsilon = ON_EPSILON ) const;
	/**
	\brief Returns the plane of the given face of this OBB.
	 The normal of the plane points outwards from this OBB, i.e. towards the space that
		is not part of the OBB.
	\param faceIndex The index of the face to get, in the range [0, 5].
	\see PointInside(), edge(), cornerPoint(), pointOnEdge(), faceCenterPoint(), facePoint(). */
	CPlane facePlane(int faceIndex) const;
		/**
		\brief Returns the point at the center of the given face of this OBB.
		\param faceIndex The index of the OBB face to generate the point at. The valid range is [0, 5].
		\todo Document which index generates which face.
		\see PointInside(), edge(), cornerPoint(), pointOnEdge(), facePoint(). */
	CVec3D faceCenterPoint(int faceIndex) const;
	/// includes touching
	bool			containsPoint( const CVec3D &p ) const;			
	/// includes touching
	bool			intersectsBox( const COBBox &a ) const;			
	
	bool			lineIntersection( const CVec3D &start, const CVec3D &end ) const;
	// intersection points are (start + dir * scale1) and (start + dir * scale2)
	bool			rayIntersection( const CVec3D &start, const CVec3D &dir, float &scale1, float &scale2 ) const;


		/**
		\brief  Tests whether this COBBox and the given object intersect.
		 Both objects are treated as "solid", meaning that if one of the objects is fully contained inside
		another, this function still returns true. (e.g. in case a line segment is contained inside this COBBox,
		or this COBBox is contained inside a Sphere, etc.)
		The first parameter of this function specifies the other object to test against.
		The COBBox-COBBox intersection test is from Christer Ericson's book Real-Time Collision Detection, p. 101-106.
		See http://realtimecollisiondetection.net/ [groupSyntax]
		\param obb The other oriented bounding box to test intersection with.
		\param epsilon The COBBox-COBBox test utilizes a SAT test to detect the intersection. A robust implementation requires
			an epsilon threshold to test that the used axes are not degenerate.
		\see contains(), distance(), closestPoint().
		\todo add intersects(Circle/Disc). */
	bool intersects(const COBBox &obb, float epsilon =CMath::EPSILON_SuperLow) const;
	bool intersects(const CAABBox &aabb) const;
	bool intersects(const CPlane &plane) const;
	/** \param dNear [out] If specified, receives the parametric distance along the line denoting where the
			line entered this COBBox.
		\param dFar [out] If specified, receives the parametric distance along the line denoting where the
			line exited this COBBox. */
	bool intersects(const CRay &ray, float &dNear, float &dFar) const;
	bool intersects(const CRay &ray) const;
	bool intersects(const CLine &line, float &dNear, float &dFar) const;
	bool intersects(const CLine &line) const;
	bool intersects(const CLineSegment &lineSegment, float &dNear, float &dFar) const;
	bool intersects(const CLineSegment &lineSegment) const;
	/** \param closestPointOnOBB [out] If specified, receives the closest point on this COBBox To the given sphere. This
			pointer may be null. */
	bool intersects(const CSphere &sphere, CVec3D *closestPointOnOBB = 0) const;
	bool intersects(const CTriangle &triangle) const;
	bool intersects(const CPolygon &polygon) const;
	bool intersects(const CPolyhedron &polyhedron) const;
	/**
		\brief  Tests if the given object is fully contained inside this OBB.
	 This function returns true if the given object lies inside this OBB, and false otherwise.
		\note The comparison is performed using less-or-equal, so the faces of this OBB count as being inside, but
			due to float inaccuracies, this cannot generally be relied upon.
		\todo add contains(Circle/Disc/Sphere/Capsule).
		\see distance(), intersects(), closestPoint(). */
	bool contains(const CVec3D &point) const;
	bool contains(const CLineSegment &lineSegment) const;
	bool contains(const CAABBox &aabb) const;
	bool contains(const COBBox &obb) const;
	bool contains(const CTriangle &triangle) const;
	bool contains(const CPolygon &polygon) const;
	bool contains(const CPolyhedron &polyhedron) const;


					// tight box for a collection of points
	void			fromPoints( const CVec3D *points, const int numPoints );
					// most tight box for a translation
	void			fromPointTranslation( const CVec3D &point, const CVec3D &translation );
	void			fromBoxTranslation( const COBBox &box, const CVec3D &translation );
					// most tight box for a rotation
//s	void			fromPointRotation( const CVec3D &point, const CRotation &rotation );
//s	void			FromBoxRotation( const COBBox &box, const CRotation &rotation );

	void			toPoints( CVec3D points[8] ) const;
	CSphere			toSphere() const;
	/**
	\brief Converts this OBB to a polyhedron.
	 This function returns a polyhedron representation of this OBB. This conversion is exact, meaning that the returned
		polyhedron represents the same set of points than this OBB. */
	CPolyhedron		toPolyhedron() const;


	/// calculates the silhouette of the box
	int				getProjectionSilhouetteVerts( const CVec3D &projectionOrigin, CVec3D silVerts[6] ) const;
	int				getParallelProjectionSilhouetteVerts( const CVec3D &projectionDir, CVec3D silVerts[6] ) const;
	/**
	\brief  Computes the closest point inside this COBBox to the given point.
	 If the target point lies inside this COBBox, then that point is returned.
	\see distance(), contains(), intersects().
	\todo add closestPoint(Line/Ray/LineSegment/Plane/Triangle/Polygon/Circle/Disc/CAABBox/COBBox/Sphere/Capsule/Frustum/Polyhedron). */
	CVec3D closestPoint(const CVec3D &point) const;

	/**
	\brief  Computes the distance between this COBBox and the given object.
	 This function finds the nearest pair of points on this and the given object, and computes their distance.
	If the two objects intersect, or one object is contained inside the other, the returned distance is zero.
	\todo add COBBox::distance(Line/Ray/LineSegment/Plane/Triangle/Polygon/Circle/Disc/CAABBox/COBBox/Capsule/Frustum/Polyhedron).
	\see contains(), intersects(), closestPoint(). */
	float distance(const CVec3D &point) const;
	float distance(const CSphere &sphere) const;
	/**
	\brief  calculates the projection, or Projects this COBBox onto the given 1D axis direction vector.
	 This function collapses this COBBox onto an 1D axis for the purposes of e.g. separate axis test computations.
		The function returns a 1D range [outMin, outMax] denoting the interval of the projection.
	\param direction The 1D axis to project to. This vector may be unnormalized, in which case the output
			of this function gets scaled by the length of this vector.
	\param outMin [out] Returns the minimum extent of this object along the projection axis.
	\param outMax [out] Returns the maximum extent of this object along the projection axis. */
	void projectToAxis(const CVec3D &direction, float &outMin, float &outMax) const;
	void projectToAxis( const CMat3D &ax, CAABBox &bounds ) const;

	/**
	\brief  Returns the tightest CAABBox that contains this COBBox.
	 This function computes the optimal minimum volume CAABBox that encloses this COBBox.
	\note Since an CAABBox cannot generally represent an COBBox, this conversion is not exact, but the returned CAABBox
			specifies a larger volume.			
	\see setFrom(), maximalContainedAABB(), minimalEnclosingSphere(), maximalContainedSphere(). */
	CAABBox minimalEnclosingAABB() const;

#if 0
	/// Returns the largest CAABBox that can fit inside this COBBox.
	/** This function computes the largest CAABBox that can fit inside this COBBox. This CAABBox is unique up to the center point of the
		CAABBox. The returned CAABBox will be centered to the center point of this COBBox.
		\see minimalEnclosingAABB(), minimalEnclosingSphere(), maximalContainedSphere(). */
	CAABBox maximalContainedAABB() const;
#endif

	/// Returns the smallest sphere that contains this COBBox.
	/** This function computes the optimal minimum volume sphere that encloses this COBBox.
		\see minimalEnclosingAABB(), maximalContainedAABB(), maximalContainedSphere(). */
	CSphere minimalEnclosingSphere() const;

	/// Returns the largest sphere that can fit inside this COBBox. [similarOverload: minimalEnclosingSphere]
	/** This function computes the largest sphere that can fit inside this COBBox. This sphere is unique up to the center point
		of the sphere. The returned sphere will be positioned to the same center point as this COBBox.
		\see minimalEnclosingSphere(), maximalContainedAABB(), maximalContainedSphere(). */
	CSphere maximalContainedSphere() const;

		/// Returns an edge of this OBB.
	/** \param edgeIndex The index of the edge line segment to get, in the range [0, 11].
		\todo Draw a diagram that shows which index generates which edge.
		\see PointInside(), cornerPoint(), pointOnEdge(), faceCenterPoint(), facePoint(). */
	CLineSegment edge(int edgeIndex) const;

	/// Returns a corner point of this OBB.
	/** This function generates one of the eight corner points of this OBB.
		\param cornerIndex The index of the corner point to generate, in the range [0, 7].
		 The points are returned in the order 0: ---, 1: --+, 2: -+-, 3: -++, 4: +--, 5: +-+, 6: ++-, 7: +++. (corresponding the XYZ axis directions).
		\todo Draw a diagram that shows which index generates which edge.
		\see PointInside(), edge(), pointOnEdge(), faceCenterPoint(), facePoint(). */
	CVec3D cornerPoint(int cornerIndex) const;

	/// Computes an extreme point of this OBB in the given direction.
	/** An extreme point is a farthest point of this OBB in the given direction. Given a direction,
		this point is not necessarily unique.
		\param direction The direction vector of the direction to find the extreme point. This vector may
			be unnormalized, but may not be null.
		\return An extreme point of this OBB in the given direction. The returned point is always a
			corner point of this OBB.
		\see cornerPoint(). */
	CVec3D extremePoint(const CVec3D &direction) const;

	/// Finds the two extreme points along the given direction vector from the given point array.
	/** \param dir The direction vector to project the point array to. This vector does not need to be normalized.
		\param pointArray [in] The list of points to process.
		\param numPoints The number of elements in pointArray.
		\param idxSmallest [out] The index of the smallest point along the given direction will be received here.
			This pointer may be left null, if this information is of no interest.
		\param idxLargest [out] The index of the largest point along the given direction will be received here.
			This pointer may be left null, if this information is of no interest. */
	static void extremePointsAlongDirection(const CVec3D &dir, const CVec3D *pointArray, int numPoints, int &idxSmallest, int &idxLargest);

public:
	/**
	\if pt_br
	\brief o ponto central desta COBBox.
	\elseif us_en
	\brief The center position of this OBBox.
	/** In the local space of the COBBox, the center of this COBBox is 
	at (extents.x,extents.y,extents.z), and the COBBox is an CAABBox with size 2*extents.
	\endif
	**/
	CVec3D			center;
	/** \if pr_br 
	 \brief armazena o vetor raio.
	 \warning Estes membros devem ser sempre positivos para que esta COBBox não seja degenerada.
	 \elseif us_en
	 \brief Stores the radius vector, or half-sizes to x, y and z directions in the local space of this COBBox. 
	 \warning These members should be positive to represent a non-degenerate COBBox.
	
	 \endif
	 **/
	CVec3D			extents;
	/**
	\if pt_br
	\brief Especifica o vetor direção normalizado dos eixos locais.
	\n  a linha 0 da matriz (axis[0]) representa a direção +X.
	\n  a linha 1 da matriz (axis[1]) representa a direçao +Y.
	\n  a linha 2 da matriz (axis[2]) representa a direção +Z.

	\elseif us_en
	\brief Specifies normalized direction vectors for the local axes. [noscript] [similarOverload: pos]
	 axis[0] specifies the +X direction in the local space of this COBBox, 
	 axis[1] the +Y direction and axis[2] the +Z direction.
		The scale of these vectors is always normalized. 
		The half-length of the COBBox along its local axes is
		specified by the vector extents.
		The axis vectors must always be orthonormal. Be sure to guarantee that condition holds if you
		directly set to this member variable. 
	\endif
	**/
	
	CMat3D			axis;
};

extern COBBox	box_zero;

SMF_INLINE COBBox::COBBox() {
}

SMF_INLINE COBBox::COBBox( const CVec3D &center, const CVec3D &extents, const CMat3D &axis ) {
	this->center = center;
	this->extents = extents;
	this->axis = axis;
}

SMF_INLINE COBBox::COBBox( const CVec3D &point ) {
	this->center = point;
	this->extents.toZero();
	this->axis.toIdentity();
}

SMF_INLINE COBBox::COBBox( const CAABBox &bounds ) {
	this->center = ( bounds[0] + bounds[1] ) * 0.5f;
	this->extents = bounds[1] - this->center; 
	this->axis.toIdentity();
}

SMF_INLINE COBBox::COBBox( const CAABBox &bounds, const CVec3D &origin, const CMat3D &axis ) {
	this->center = ( bounds[0] + bounds[1] ) * 0.5f;
	//We can also refer to the “radius vector” extents of the box, which is the vector from the center to pmax:
	// = pmax - center ou =size / 2
	this->extents = bounds[1] - this->center;  
	this->center = origin + this->center * axis;
	this->axis = axis;
}

SMF_INLINE COBBox COBBox::operator+( const CVec3D &t ) const {
	return COBBox( center + t, extents, axis );
}

SMF_INLINE COBBox &COBBox::operator+=( const CVec3D &t ) {
	center += t;
	return *this;
}

SMF_INLINE COBBox COBBox::operator*( const CMat3D &r ) const {
	return COBBox( center * r, extents, axis * r );
}

SMF_INLINE COBBox &COBBox::operator*=( const CMat3D &r ) {
	center *= r;
	axis *= r;
	return *this;
}

SMF_INLINE COBBox COBBox::operator+( const COBBox &a ) const {
	COBBox newBox;
	newBox = *this;
	newBox.AddBox( a );
	return newBox;
}

SMF_INLINE COBBox &COBBox::operator+=( const COBBox &a ) {
	COBBox::AddBox( a );
	return *this;
}

SMF_INLINE COBBox COBBox::operator-( const COBBox &a ) const {
	return COBBox( center, extents - a.extents, axis );
}

SMF_INLINE COBBox &COBBox::operator-=( const COBBox &a ) {
	extents -= a.extents;
	return *this;
}

SMF_INLINE bool COBBox::compare( const COBBox &a ) const {
	return ( center.compare( a.center ) && extents.compare( a.extents ) && axis.compare( a.axis ) );
}

SMF_INLINE bool COBBox::compare( const COBBox &a, const float epsilon ) const {
	return ( center.compare( a.center, epsilon ) && extents.compare( a.extents, epsilon ) && axis.compare( a.axis, epsilon ) );
}

SMF_INLINE bool COBBox::operator==( const COBBox &a ) const {
	return compare( a );
}

SMF_INLINE bool COBBox::operator!=( const COBBox &a ) const {
	return !compare( a );
}

SMF_INLINE void COBBox::clear() {
	center.toZero();
	extents[0] = extents[1] = extents[2] = -CMath::INFINITY_FLOAT;
	axis.toIdentity();
}

SMF_INLINE void COBBox::toZero() {
	center.toZero();
	extents.toZero();
	axis.toIdentity();
}

SMF_INLINE const CVec3D &COBBox::getCenter() const {
	return center;
}

SMF_INLINE const CVec3D &COBBox::getExtents() const {
	return extents;
}

SMF_INLINE const CMat3D &COBBox::getAxis() const {
	return axis;
}

SMF_INLINE float COBBox::getVolume() const {
	return ( extents * 2.0f ).getLengthSqr();
}

SMF_INLINE bool COBBox::isCleared() const {
	return extents[0] < 0.0f;
}

SMF_INLINE COBBox COBBox::expand( const float d ) const {
	return COBBox( center, extents + CVec3D( d, d, d ), axis );
}

SMF_INLINE COBBox &COBBox::expandSelf( const float d ) {
	extents[0] += d;
	extents[1] += d;
	extents[2] += d;
	return *this;
}

SMF_INLINE COBBox COBBox::translate( const CVec3D &translation ) const {
	return COBBox( center + translation, extents, axis );
}

SMF_INLINE COBBox &COBBox::translateSelf( const CVec3D &translation ) {
	center += translation;
	return *this;
}

SMF_INLINE COBBox COBBox::rotate( const CMat3D &rotation ) const {
	return COBBox( center * rotation, extents, axis * rotation );
}

SMF_INLINE COBBox &COBBox::rotateSelf( const CMat3D &rotation ) {
	center *= rotation;
	axis *= rotation;
	return *this;
}

SMF_INLINE bool COBBox::containsPoint( const CVec3D &p ) const {
	CVec3D lp = p - center;
	if ( CMath::fabs( lp * axis.Row(0) ) > extents[0] ||
			CMath::fabs( lp * axis.Row(1) ) > extents[1] ||
				CMath::fabs( lp * axis.Row(2) ) > extents[2] ) {
		return false;
	}
	return true;
}

SMF_INLINE CSphere COBBox::toSphere() const {
	return CSphere( center, extents.getLenght() );
}

SMF_INLINE void COBBox::projectToAxis( const CVec3D &dir, float &min, float &max ) const {
	float d1 = dir * center;
	float d2 = CMath::fabs( extents[0] * ( dir * axis.Row(0) ) ) +
				CMath::fabs( extents[1] * ( dir * axis.Row(1) ) ) +
				CMath::fabs( extents[2] * ( dir * axis.Row(2) ) );
	min = d1 - d2;
	max = d1 + d2;
}
#if 0
void COBBox::projectToAxis(const CVec3D &direction, float &outMin, float &outMax) const
{
	float x = CMath::fabs((direction* axis[0]) * extents.x);
	float y = CMath::fabs((direction* axis[1]) * extents.y);
	float z = CMath::fabs((direction* axis[2]) * extents.z);
	float pt = direction* center;
	outMin = pt - x - y - z;
	outMax = pt + x + y + z;
}
#endif
SMF_INLINE void COBBox::projectToAxis( const CMat3D &ax, CAABBox &bounds ) const {
	for ( int i = 0; i < 3; i++ ) {
		float d1 = ax.Row(i) * center;
		float d2 = CMath::fabs( extents[0] * ( ax.Row(i) * axis.Row(0) ) ) +
					CMath::fabs( extents[1] * ( ax.Row(i) * axis.Row(1) ) ) +
					CMath::fabs( extents[2] * ( ax.Row(i) * axis.Row(2) ) );
		bounds[0][i] = d1 - d2;
		bounds[1][i] = d1 + d2;
	}
}
}  // end GEO
} //end SMF
#endif /* !__SMF_BOX_H__ */
