#ifndef __SMF_RAY_
#define __SMF_RAY_

#include "../SMF_Config.h"
#include "SMF_Ray.h"
#include "../math/all.h"
#include "all.h"

namespace SMF {
using namespace MATH;
using namespace Util;
namespace MATH{
class CQuaternion;
}
namespace GEO{


#ifdef MATH_OGRE_INTEROP
#include <OgreRay.h>
#endif

class CLine;
class CLineSegment;
class CPolygon;
class CPolyhedron;
class CCircle;
 


/**
 * \class CRay
 *
 * \ingroup SMF_Geometric
 * \image html pics\line.png
 * \if pt_br
 * \brief Um Raio em espaço 3D que inicia num ponto e se extende até o infinito em um direçao
 *       
 * \elseif us_en
 * \brief 	A ray in 3D space is a line that starts from an origin point and extends to infinity in one direction.
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CRay
{
public:
	/**\brief  Specifies the origin of this ray.  **/
	CVec3D pos;

	/// The normalized direction vector of this ray. [similarOverload: pos]
	/** \note For proper functionality, this direction vector needs to always be normalized. If you set to this
		member manually, remember to make sure you only assign normalized direction vectors. */
	CVec3D dir;

	/// The default constructor does not initialize any members of this class.
	/** This means that the values of the members pos and dir are undefined after creating a new CRay using this
		default constructor. Remember to assign to them before use.
		\see pos, dir. */
	CRay() {}

	/// Constructs a new ray by explicitly specifying the member variables.
	/** \param pos The origin position of the ray.
		\param dir The direction of the ray. This vector must be normalized, this function will not normalize
			the vector for you (for performance reasons).
		\see pos, dir. */
	CRay(const CVec3D &pos, const CVec3D &dir);

	/// Converts a CLine to a CRay.
	/** This conversion simply copies the members pos and dir over from the given CLine to this CRay.
		This means that the new CRay starts at the same position, but only extends to one direction in space,
		instead of two.
		\see class CLine, toLine(). */
	explicit CRay(const CLine &line);

	/// Converts a CLineSegment to a CRay.
	/** This constructor sets pos = lineSegment.a, and dir = (lineSegment.b - lineSegment.a).normalized().
		\see class CLineSegment, toLineSegment(). */
	explicit CRay(const CLineSegment &lineSegment);

	bool isFinite() const;

	/// Gets a point along the ray at the given distance.
	/** Use this function to convert a 1D parametric point along the CRay to a 3D point in the linear space.
		\param distance The point to compute. getPoint(0) will return pos. getPoint(t) will return a point
			at distance |t| from pos. Passing in negative values is allowed, but in that case, the
			returned point does not actually lie on this CRay.
		\return pos + distance * dir.
		\see pos, dir. */
	CVec3D getPoint(float distance) const;

	/// Translates this CRay in world space.
	/** \param offset The amount of displacement to apply to this CRay, in world space coordinates.
		\see transform(). */
	void translate(const CVec3D &offset);

	/// Applies a transformation to this CRay, in-place.
	/** See translate(), classes CMat3D, CMatJoint3x4, CMat4D, CQuaternion. */
	void transform(const CMat3D &transform);
	void transform(const CMatJoint3x4 &transform);
	void transform(const CMat4D &transform);
	void transform(const CQuaternion &transform);

	/// Tests if the given object is fully contained on this ray.
	/** \param distanceThreshold The magnitude of the epsilon test threshold to use. Since a CRay
		is a 1D object in a 3D space, an epsilon threshold is used to allow errors caused by floating-point
		inaccuracies.
		\return True if this ray contains the given object, up to the given distance threshold.
		\see class CLineSegment, distance(), closestPoint(), intersects(). */
	bool contains(const CVec3D &point, float distanceThreshold = 1e-3f) const;
	bool contains(const CLineSegment &lineSegment, float distanceThreshold = 1e-3f) const;

	/// Tests if two rays are equal.
	/** \return True if this and the given CRay represent the same set of points, up to the given epsilon. */
	bool compare(const CRay &otherRay, float epsilon =CMath::EPSILON_SuperLow) const;

	/// Computes the distance between this ray and the given object.
	/** This function finds the nearest pair of points on this and the given object, and computes their distance.
		If the two objects intersect, or one object is contained inside the other, the returned distance is zero.
		\param d [out] If specified, receives the parametric distance along this ray that
			specifies the closest point on this ray to the given object. The value returned here can be negative.
			This pointer may be null.
		\see contains(), intersects(), closestPoint(), getPoint(). */
	float distance(const CVec3D &point, float *d) const;
	float distance(const CVec3D &point) const;

	/** \param d2 [out] If specified, receives the parametric distance along the other line that specifies the
		closest point on that line to this ray. The value returned here can be negative. This pointer may
		be null. */
	float distance(const CRay &other, float *d, float *d2 = 0) const;
	float distance(const CRay &other) const;
	float distance(const CLine &other, float *d, float *d2 = 0) const;
	float distance(const CLine &other) const;
	float distance(const CLineSegment &other, float *d, float *d2 = 0) const;
	float distance(const CLineSegment &other) const;
	float distance(const CSphere &sphere) const;
//	float distance(const Capsule &capsule) const;

	/// Computes the closest point on this ray to the given object.
	/** If the other object intersects this ray, this function will return an arbitrary point inside
		the region of intersection.
		\param d [out] If specified, receives the parametric distance along this ray that
			specifies the closest point on this ray to the given object. The value returned here can be negative.
			This pointer may be null.
		\see contains(), distance(), intersects(), getPoint(). */
	CVec3D closestPoint(const CVec3D &targetPoint, float *d = 0) const;
	/** \param d2 [out] If specified, receives the parametric distance along the other line that specifies the
		closest point on that line to this ray. The value returned here can be negative. This pointer may
		be null. */
	CVec3D closestPoint(const CRay &other, float *d = 0, float *d2 = 0) const;
	CVec3D closestPoint(const CLine &other, float *d = 0, float *d2 = 0) const;
	CVec3D closestPoint(const CLineSegment &other, float *d = 0, float *d2 = 0) const;

	/// Tests whether this ray and the given object intersect.	
	/** Both objects are treated as "solid", meaning that if one of the objects is fully contained inside
		another, this function still returns true.
		\param d [out] If specified, this parameter will receive the parametric distance of
			the intersection point along this object. Use the getPoint(d) function
			to get the actual point of intersection. This pointer may be null.
		\param intersectionPoint [out] If specified, receives the actual point of intersection. This pointer
			may be null.
		\return True if an intersection occurs or one of the objects is contained inside the other, false otherwise.
		\see contains(), distance(), closestPoint(), getPoint(). */
	bool intersects(const CTriangle &triangle, float *d, CVec3D *intersectionPoint) const;
	bool intersects(const CTriangle &triangle) const;
	bool intersects(const CPlane &plane, float *d) const;
	bool intersects(const CPlane &plane) const;
	/** \param intersectionNormal [out] If specified, receives the surface normal of the other object at
		the point of intersection. This pointer may be null. */
	bool intersects(const CSphere &s, CVec3D *intersectionPoint, CVec3D *intersectionNormal, float *d) const;
	bool intersects(const CSphere &s) const;
	/** \param dNear [out] If specified, receives the distance along this ray to where the ray enters
		the bounding box.
		\param dFar [out] If specified, receives the distance along this ray to where the ray exits
		the bounding box. */
	bool intersects(const CAABBox &aabb, float &dNear, float &dFar) const;
	bool intersects(const CAABBox &aabb) const;
	bool intersects(const COBBox &obb, float &dNear, float &dFar) const;
	bool intersects(const COBBox &obb) const;
//	bool intersects(const Capsule &capsule) const;
	bool intersects(const CPolygon &polygon) const;
//	bool intersects(const Frustum &frustum) const;
	bool intersects(const CPolyhedron &polyhedron) const;

	bool intersects(const CPoint2D &p) const;

	/// Tests if this ray intersects the given disc.
	/// \todo This signature will be moved to bool intersects(const Disc &disc) const;
	bool intersectsDisc(const CCircle &disc) const;

	/// Converts this CRay to a CLine.
	/** The pos and dir members of the returned CLine will be equal to this CRay. The only difference is
		that a CLine extends to infinity in two directions, whereas the CRay spans only in the positive
		direction.
		\see dir, CRay::CRay, class CLine, toLineSegment(). */
	CLine toLine() const;
	/// Converts this CRay to a CLineSegment.
	/** \param d Specifies the position of the other endpoint along this CRay. This parameter may be negative,
		in which case the returned CLineSegment does not lie inside this CRay.
		\return A CLineSegment with point a at pos, and point b at pos + d * dir.
		\see pos, dir, CRay::CRay, class CLineSegment, toLine(). */
	CLineSegment toLineSegment(float d) const;

	/// Projects this CRay onto the given 1D axis direction vector.
	/** This function collapses this CRay onto an 1D axis for the purposes of e.g. separate axis test computations.
		The function returns a 1D range [outMin, outMax] denoting the interval of the projection.
		\param direction The 1D axis to project to. This vector may be unnormalized, in which case the output
			of this function gets scaled by the length of this vector.
		\param outMin [out] Returns the minimum extent of this object along the projection axis.
		\param outMax [out] Returns the maximum extent of this object along the projection axis. */
	void projectToAxis(const CVec3D &direction, float &outMin, float &outMax) const;

	/// Converts this CRay to a CLineSegment.
	/** \param dStart Specifies the position of the first endpoint along this CRay. This parameter may be negative,
		in which case the starting point lies outside this CRay to the opposite direction of the CRay.
		\param dEnd Specifies the position of the second endpoint along this CRay. This parameter may also be negative.
		\return A CLineSegment with point a at pos + dStart * dir, and point b at pos + dEnd * dir.
		\see pos, dir, CRay::CRay, class CLineSegment, toLine(). */
	CLineSegment toLineSegment(float dStart, float dEnd) const;

	/// Returns a human-readable representation of this CRay.
	/** The returned string specifies the position and direction of this CRay. */
	std::string toString() const;

#ifdef MATH_QT_INTEROP
	operator QString() const { return toString(); }
	QString toString() const { return QString::fromStdString(toString()); }
#endif

#ifdef MATH_OGRE_INTEROP
	CRay(const Ogre::CRay &other) { pos = other.getOrigin(); dir = other.getDirection(); }
	operator Ogre::CRay() const { return Ogre::CRay(pos, dir); }
#endif

};

/// \note Assumes that transform may contain scaling, and re-normalizes the ray direction
///		after the transform.
CRay operator *(const CMat3D &transform, const CRay &ray);
/// \note Assumes that transform may contain scaling, and re-normalizes the ray direction
///		after the transform.
CRay operator *(const CMatJoint3x4 &transform, const CRay &ray);
/// \note Assumes that transform may contain scaling, and re-normalizes the ray direction
///		after the transform.
CRay operator *(const CMat4D &transform, const CRay &ray);
CRay operator *(const CQuaternion &transform, const CRay &ray);

#ifdef MATH_QT_INTEROP
Q_DECLARE_METATYPE(CRay)
Q_DECLARE_METATYPE(CRay*)
#endif


std::ostream &operator <<(std::ostream &o, const CRay &ray);



} //end GEO
}  //end SMF

#endif // __SMF_RAY_
