#ifndef __SMF_LINESEGMENT_
#define __SMF_LINESEGMENT_

#include "../SMF_Config.h"
#include "../math/all.h"
#include "all.h"

namespace SMF {
using namespace MATH;
using namespace Util;
namespace GEO{


/**
 * \class CLineSegment
 *
 * \ingroup SMF_Geometric
 * \image html pics\line.png
 * \if pt_br 
 *
 * \brief Um Segmento de reta em espaço 3D, definido por dois pontos.
 *       
 * \elseif us_en
 * \brief 	A line segment in 3D space is a finite line with a start and end point.
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CLineSegment
{
public:
	/// The starting point of this line segment.
	CVec3D begin;
	/// The end point of this line segment. [similarOverload: a]
	CVec3D end;

	/// The default constructor does not initialize any members of this class.
	/** This means that the values of the members begin and end are undefined after creating a new CLineSegment using this
		default constructor. Remember to assign to them before use.
		\see begin, end. */
	CLineSegment() {}

	/// Constructs a line segment through the given end points.
	/** \see begin, end. */
	CLineSegment(const CVec3D &a, const CVec3D &b);

	/// Constructs a line segment from a ray or a line.
	/** This constructor takes the ray/line origin position as the starting point of this line segment, and defines the end point
		of the line segment using the given distance parameter.
		\param d The distance along the ray/line for the end point of this line segment. This will set b = ray.pos + d * ray.dir
			as the end point. When converting a ray to a line segment, it is possible to pass in a d value < 0, but in that case
			the resulting line segment will not lie on the ray.
		\see begin, end, classes CRay, CLine, CLine::getPoint(), CRay::getPoint(). */
	explicit CLineSegment(const CRay &ray, float d);
	explicit CLineSegment(const CLine &line, float d);

	/// Returns a point on the line.
	/** \param d The normalized distance along the line segment to compute. If a value in the range [0, 1] is passed, then the
			returned point lies along this line segment. If some other value is specified, the returned point lies on the
			line defined by this line segment, but not inside the interval from begin to  end.
		\note The meaning of d here differs from CLine::getPoint and CRay::getPoint. For the class CLineSegment,
			getPoint(0) returns begin, and getPoint(1) returns  end. This means that getPoint(1) will not generally be exactly one unit
			away from the starting point of this line segment, as is the case with CLine and CRay.
		\return (1-d)*begin + d* end.
		\see begin, end, CLine::getPoint(), CRay::getPoint(). */
	CVec3D getPoint(float d) const;

	/// Returns the center point of this line segment.
	/** This function is the same as calling getPoint(0.5f), but provided here as conveniency.
		\see getPoint(). */
	CVec3D centerPoint() const;

	/// Reverses the direction of this line segment.
	/** This function swaps the start and end points of this line segment so that it runs from  end to a.
		This does not have an effect on the set of points represented by this line segment, but it reverses
		the direction of the vector returned by getDir().
		\note This function operates in-place.
		\see begin, end, getDir(). */
	void reverse();

	/// Returns the normalized direction vector that points in the direction a-> end.
	/** \note The returned vector is normalized, meaning that its length is 1, not | end-a|.
		\see begin, end. */
	CVec3D getDir() const;

	/// Computes an extreme point of this CLineSegment in the given direction.
	/** An extreme point is a farthest point along this CLineSegment in the given direction. Given a direction,
		this point is not necessarily unique.
		\param direction The direction vector of the direction to find the extreme point. This vector may
			be unnormalized, but may not be null.
		\return An extreme point of this CLineSegment in the given direction. The returned point is always
			either a or  end.
		\see begin, end.*/
	CVec3D extremePoint(const CVec3D &direction) const;

	/// Translates this CLineSegment in world space.
	/** \param offset The amount of displacement to apply to this CLineSegment, in world space coordinates.
		\see transform(). */
	void translate(const CVec3D &offset);

	/// Applies a transformation to this line.
	/** This function operates in-place.
		\see translate(), classes CMat3D, CMatJoint3x4, CMat4D, CQuaternion, transform(). */
	void transform(const CMat3D &transform);
	void transform(const CMatJoint3x4 &transform);
	void transform(const CMat4D &transform);
	void transform(const CQuaternion &transform);

	/// Computes the length of this line segment.
	/** \return | end-begin|.
		\see begin,  end. */
	float getLenght() const;
	/// Computes the squared length of this line segment.
	/** Calling this function is faster than calling getLenght(), since this function avoids computing a square root.
		If you only need to compare lengths to each other and are not interested in the actual length values,
		you can compare by using getLengthSqr(), instead of getLenght(), since sqrt() is an order-preserving
		(monotonous and non-decreasing) function. [similarOverload: length] */
	float getLengthSqr() const;

	/// Tests if this line segment is finite.
	/** A line segment is <b><i>finite</i></b> if its endpoints begin and  end do not contain floating-point NaNs or +/-infs
		in them.
		\return True if both IntersectResult and  end have finite floating-point values. */
	bool isFinite() const;

	/// Tests if this line segment represents the same set of points than the given line segment.
	/** \param distanceThreshold Specifies how much distance threshold to allow in the comparison.
		\return True if begin == rhs.begin &&  end == rhs. end, or, begin == rhs. end &&  end = rhs.begin, within the given epsilon. */
	bool compare(const CLineSegment &rhs, float distanceThreshold = 1e-3f) const;

	/// Tests if the given point or line segment is contained on this line segment.
	/** \param distanceThreshold Because a line segment is an one-dimensional object in 3D space, an epsilon value
			is used as a threshold for this test. This effectively transforms this line segment to a capsule with
			the radius indicated by this value.
		\return True if this line segment contains the given point or line segment.
		\see intersects, closestPoint(), distance(). */
	bool contains(const CVec3D &point, float distanceThreshold = 1e-3f) const;
	bool contains(const CLineSegment &lineSegment, float distanceThreshold = 1e-3f) const;

	/// Computes the closest point on this line segment to the given object.
	/** \param d [out] If specified, this parameter receives the normalized distance along
			this line segment which specifies the closest point on this line segment to
			the specified point. This pointer may be null.
		\return The closest point on this line segment to the given object.
		\see contains(), distance(), intersects(). */
	CVec3D closestPoint(const CVec3D &point, float *d = 0) const;
	/** \param d2 [out] If specified, this parameter receives the (normalized, in case of line segment)
			distance along the other line object which specifies the closest point on that line to
			this line segment. This pointer may be null. */
	CVec3D closestPoint(const CRay &other, float *d = 0, float *d2 = 0) const;
	CVec3D closestPoint(const CLine &other, float *d = 0, float *d2 = 0) const;
	CVec3D closestPoint(const CLineSegment &other, float *d = 0, float *d2 = 0) const;

	/// Computes the distance between this line segment and the given object.
	/** \param d [out] If specified, this parameter receives the normalized distance along
			this line segment which specifies the closest point on this line segment to
			the specified point. This pointer may be null.
		\return The distance between this line segment and the given object.
		\see Constains(), closestPoint(), intersects(). */
	float distance(const CVec3D &point, float *d = 0) const;
	/** \param d2 [out] If specified, this parameter receives the (normalized, in case of line segment)
			distance along the other line object which specifies the closest point on that line to
			this line segment. This pointer may be null. */
	float distance(const CRay &other, float *d = 0, float *d2 = 0) const;
	float distance(const CLine &other, float *d = 0, float *d2 = 0) const;
	float distance(const CLineSegment &other, float *d = 0, float *d2 = 0) const;
	float distance(const CPlane &other) const;
	float distance(const CSphere &other) const;
//	float distance(const Capsule &other) const;

	/// Tests whether this line segment and the given object intersect.	
	/** Both objects are treated as "solid", meaning that if one of the objects is fully contained inside
		another, this function still returns true. (for example, if this line segment is contained inside a sphere)
		\todo Output intersection point. */
	bool intersects(const CPlane &plane) const;
	/** \param d [out] If specified, this parameter receives the normalized distance along
			this line segment which specifies the closest point on this line segment to
			the specified point. This pointer may be null. */
	bool intersects(const CPlane &plane, float *d) const;
	/** \param intersectionPoint [out] If specified, receives the point of intersection. This pointer may be null. */
	bool intersects(const CTriangle &triangle, float *d, CVec3D *intersectionPoint) const;
	/** \param intersectionNormal [out] If specified, receives the normal vector of the other object at the point of intersection.
			This pointer may be null. */
	bool intersects(const CSphere &s, CVec3D *intersectionPoint = 0, CVec3D *intersectionNormal = 0, float *d = 0) const;
	/** \param dNear [out] If specified, receives the parametric distance along this line segment denoting where the line entered the
			bounding box object.
		\param dFar [out] If specified, receives the parametric distance along this line segment denoting where the line exited the
			bounding box object. */
	bool intersects(const CAABBox &aabb, float &dNear, float &dFar) const;
	bool intersects(const CAABBox &aabb) const;
	bool intersects(const COBBox &obb, float &dNear, float &dFar) const;
	bool intersects(const COBBox &obb) const;
//	bool intersects(const Capsule &capsule) const;
	bool intersects(const CPolygon &polygon) const;
//	bool intersects(const Frustum &frustum) const;
	bool intersects(const CPolyhedron &polyhedron) const;
	/** \param epsilon If testing intersection between two line segments, a distance threshold value is used to account
			for floating-point inaccuracies. */
	bool intersects(const CLineSegment &lineSegment, float epsilon =CMath::EPSILON_SuperLow) const;
	/// Tests if this line segment intersects the given disc.
	/// \todo This signature will be moved to bool intersects(const Disc &disc) const;
	bool intersectsDisc(const CCircle &disc) const;

	/// Converts this CLineSegment to a CRay.
	/** The pos member of the returned CRay will be equal to a, and the dir member equal to getDir().
		\see class CRay, toLine(). */
	CRay toRay() const;
	/// Converts this CLineSegment to a CLine.
	/** The pos member of the returned CLine will be equal to a, and the dir member equal to getDir().
		\see class CLine, toRay(). */
	CLine toLine() const;

	/**
	\brief Projects this CLineSegment onto the given 1D axis direction vector.
	 This function collapses this CLineSegment onto an 1D axis for the purposes of e.g. separate axis test computations.
	The function returns a 1D range [outMin, outMax] denoting the interval of the projection.
	\param direction The 1D axis to project to. This vector may be unnormalized, in which case the output
			of this function gets scaled by the length of this vector.
	\param outMin [out] Returns the minimum extent of this object along the projection axis.
	\param outMax [out] Returns the maximum extent of this object along the projection axis. */
	void projectToAxis(const CVec3D &direction, float &outMin, float &outMax) const;
    

	/// Returns a human-readable representation of this CLineSegment. Most useful for debugging purposes.
	std::string toString() const;

#ifdef MATH_QT_INTEROP
	operator QString() const { return toString(); }
	QString toString() const { return QString::fromStdString(toString()); }
#endif

#ifdef MATH_GRAPHICSENGINE_INTEROP
	void toLineList(VertexBuffer &vb) const;
#endif
};

CLineSegment operator *(const CMat3D &transform, const CLineSegment &line);
CLineSegment operator *(const CMatJoint3x4 &transform, const CLineSegment &line);
CLineSegment operator *(const CMat4D &transform, const CLineSegment &line);
CLineSegment operator *(const CQuaternion &transform, const CLineSegment &line);

#ifdef MATH_QT_INTEROP
Q_DECLARE_METATYPE(CLineSegment)
Q_DECLARE_METATYPE(CLineSegment*)
#endif

std::ostream &operator <<(std::ostream &o, const CLineSegment &lineSegment);



} //end GEO
}  //end SMF

#endif // __SMF_LINESEGMENT_
