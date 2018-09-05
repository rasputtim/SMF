#ifndef __SMF_LINE_
#define __SMF_LINE_

#include "../SMF_Config.h"
#include "../math/all.h"
#include "all.h"

namespace SMF {
using namespace MATH;
using namespace Util;
namespace GEO{

/**
 * \class CLine
 *
 * \ingroup SMF_Geometric
 *
 * \image html pics\line.png
 * \if pt_br
 * \brief   Uma Linha em espaço 3D
 * \note Uma linha em espaço 3D é definida por um ponto (origem) e uma direção, e se extene ao infinito em ambas as direções 
 * \elseif us_en
 * \brief 	A Line in 3D space.
   \note  A line in 3D space is defined by an origin point and a direction, and extends to infinity in two directions.
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CLine
{
public:
	/// Specifies the origin of this line.
	CVec3D pos;

	/// The normalized direction vector of this ray. [similarOverload: pos]
	/** \note For proper functionality, this direction vector needs to always be normalized. If you set to this
		member manually, remember to make sure you only assign normalized direction vectors. */
	CVec3D dir;

	/// The default constructor does not initialize any members of this class.
	/** This means that the values of the members pos and dir are undefined after creating a new CLine using this
		default constructor. Remember to assign to them before use.
		\see pos, dir. */
	CLine() {}

	/// Constructs a new line by explicitly specifying the member variables.
	/** \param pos The origin position of the line.
		\param dir The direction of the line. This vector must be normalized, this function will not normalize
			the vector for you (for performance reasons).
		\see pos, dir. */
	CLine(const CVec3D &pos, const CVec3D &dir);

	/// Converts a CRay to a CLine.
	/** This conversion simply copies the members pos and dir over from the given CRay to this CLine.
		This means that the new CLine starts at the same position, but extends to two directions in space,
		instead of one.
		\see class CRay, toRay(). */
	explicit CLine(const CRay &ray);

	/// Converts a CLineSegment to a CLine.
	/** This constructor sets pos = lineSegment.a, and dir = (lineSegment.b - lineSegment.a).normalized().
		\see class CLineSegment, toLineSegment(). */
	explicit CLine(const CLineSegment &lineSegment);

	bool isFinite() const;

	/// Gets a point along the line at the given distance.
	/** Use this function to convert a 1D parametric point along the CLine to a 3D point in the linear space.
		\param distance The point to compute. getPoint(0) will return pos. getPoint(t) will return a point
			at distance |t| from pos, towards the direction specified by dir. If a negative value is specified,
			a point towards the direction -dir is returned.
		\return pos + distance * dir.
		\see pos, dir. */
	CVec3D getPoint(float distance) const;

	/// Translates this CLine in world space.
	/** \param offset The amount of displacement to apply to this CLine, in world space coordinates.
		\see transform(). */
	void translate(const CVec3D &offset);

	/// Applies a transformation to this line, in-place.
	/** \see translate(), classes CMat3D, CMatJoint3x4, CMat4D, CQuaternion. */
	void transform(const CMat3D &transform);
	void transform(const CMatJoint3x4 &transform);
	void transform(const CMat4D &transform);
	void transform(const CQuaternion &transform);

	/// Tests if the given object is fully contained on this line.
	/** \param distanceThreshold The magnitude of the epsilon test threshold to use. Since a CLine
		is a 1D object in a 3D space, an epsilon threshold is used to allow errors caused by floating-point
		inaccuracies.
		\return True if this line contains the given object, up to the given distance threshold.
		\see class CLineSegment, class CRay, distance(), closestPoint(), intersects(). */
	bool contains(const CVec3D &point, float distanceThreshold = 1e-3f) const;
	bool contains(const CRay &ray, float distanceThreshold = 1e-3f) const;
	bool contains(const CLineSegment &lineSegment, float distanceThreshold = 1e-3f) const;

	/// Tests if two lines are equal.
	/** This function tests for set equality (not just member value equality). This means that the pos and dir parameters
		of either line can be completely different, as long as the set of points on the both lines are equal.
		\return True if this and the given CLine represent the same set of points, up to the given epsilon. */
	bool compare(const CLine &line, float epsilon =CMath::EPSILON_SuperLow) const;

	/// Computes the distance between this line and the given object.
	/** This function finds the nearest pair of points on this and the given object, and computes their distance.
		If the two objects intersect, or one object is contained inside the other, the returned distance is zero.
		\param d [out] If specified, receives the parametric distance along this line that
			specifies the closest point on this line to the given object. The value returned here can be negative.
			This pointer may be null.
		\see contains(), intersects(), closestPoint(), getPoint(). */
	float distance(const CVec3D &point, float *d = 0) const;
	/** \param d2 [out] If specified, receives the parametric distance along the other line that specifies the
		closest point on that line to this line. The value returned here can be negative. This pointer may
		be null. */
	float distance(const CRay &other, float *d, float *d2 = 0) const;
	float distance(const CRay &other) const;
	float distance(const CLine &other, float *d, float *d2 = 0) const;
	float distance(const CLine &other) const;
	float distance(const CLineSegment &other, float *d, float *d2 = 0) const;
	float distance(const CLineSegment &other) const;
	float distance(const CSphere &other) const;
//	float distance(const Capsule &other) const;

	/// Computes the closest point on this line to the given object.
	/** If the other object intersects this line, this function will return an arbitrary point inside
		the region of intersection.
		\param d [out] If specified, receives the parametric distance along this line that
			specifies the closest point on this line to the given object. The value returned here can be negative.
			This pointer may be null.
		\see contains(), distance(), intersects(), getPoint(). */
	CVec3D closestPoint(const CVec3D &targetPoint, float *d = 0) const;
	/** \param d2 [out] If specified, receives the parametric distance along the other line that specifies the
		closest point on that line to this line. The value returned here can be negative. This pointer may
		be null. */
	CVec3D closestPoint(const CRay &other, float *d = 0, float *d2 = 0) const;
	CVec3D closestPoint(const CLine &other, float *d = 0, float *d2 = 0) const;
	CVec3D closestPoint(const CLineSegment &other, float *d = 0, float *d2 = 0) const;
	/** \param outU [out] If specified, receives the barycentric U-coordinate (in two-coordinate barycentric UV convention)
			representing the closest point on the triangle to this line. This pointer may be null.
		\param outV [out] If specified, receives the barycentric V-coordinate (in two-coordinate barycentric UV convention)
			representing the closest point on the triangle to this line. This pointer may be null.
		\see contains(), distance(), intersects(), getPoint(), CTriangle::Point(float u, float v). */
	CVec3D closestPoint(const CTriangle &triangle, float *outU = 0, float *outV = 0, float *d = 0) const;

	/// Tests whether this line and the given object intersect.	
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
	bool intersects(const CPlane &plane, float *d) const;
	/** \param intersectionNormal [out] If specified, receives the surface normal of the other object at
		the point of intersection. This pointer may be null. */
	bool intersects(const CSphere &s, CVec3D *intersectionPoint = 0, CVec3D *intersectionNormal = 0, float *d = 0) const;
	/** \param dNear [out] If specified, receives the distance along this line to where the line enters
		the bounding box.
		\param dFar [out] If specified, receives the distance along this line to where the line exits
		the bounding box. */
	bool intersects(const CAABBox &aabb, float &dNear, float &dFar) const;
	bool intersects(const CAABBox &aabb) const;
	bool intersects(const COBBox &obb, float &dNear, float &dFar) const;
	bool intersects(const COBBox &obb) const;
	//bool intersects(const Capsule &capsule) const;
	bool intersects(const CPolygon &polygon) const;
//	bool intersects(const Frustum &frustum) const;
	bool intersects(const CPolyhedron &polyhedron) const;
	/// Tests if this ray intersects the given disc.
	/// \todo This signature will be moved to bool intersects(const Disc &disc) const;
	bool intersectsDisc(const CCircle &disc) const;

	/// Converts this CLine to a CRay.
	/** The pos and dir members of the returned CRay will be equal to this CLine. The only difference is
		that a CLine extends to infinity in two directions, whereas the returned CRay spans only in
		the positive direction.
		\see dir, CLine::CLine, class CRay, toLineSegment(). */
	CRay toRay() const;

	/// Converts this CLine to a CLineSegment.
	/** \param d Specifies the position of the other endpoint along this CLine. This parameter may be negative.
		\return A CLineSegment with point a at pos, and point b at pos + d * dir.
		\see pos, dir, CLine::CLine, class CLineSegment, toRay(). */
	CLineSegment toLineSegment(float d) const;

	/// Converts this CLine to a CLineSegment.
	/** \param dStart Specifies the position of the first endpoint along this CLine. This parameter may be negative,
		in which case the starting point lies to the opposite direction of the CLine.
		\param dEnd Specifies the position of the second endpoint along this CLine. This parameter may also be negative.
		\return A CLineSegment with point a at pos + dStart * dir, and point b at pos + dEnd * dir.
		\see pos, dir, CLine::CLine, class CLineSegment, toLine(). */
	CLineSegment toLineSegment(float dStart, float dEnd) const;

	/// Projects this CLine onto the given 1D axis direction vector.
	/** This function collapses this CLine onto an 1D axis for the purposes of e.g. separate axis test computations.
		The function returns a 1D range [outMin, outMax] denoting the interval of the projection.
		\param direction The 1D axis to project to. This vector may be unnormalized, in which case the output
			of this function gets scaled by the length of this vector.
		\param outMin [out] Returns the minimum extent of this object along the projection axis.
		\param outMax [out] Returns the maximum extent of this object along the projection axis. */
	void projectToAxis(const CVec3D &direction, float &outMin, float &outMax) const;

	/// Tests if the given three points are collinear.
	/** This function tests whether the given three functions all lie on the same line.
		\param epsilon The comparison threshold to use to account for floating-point inaccuracies. */
	static bool areCollinear(const CVec3D &p1, const CVec3D &p2, const CVec3D &p3, float epsilon =CMath::EPSILON_SuperLow);

	/// Computes the closest point pair on two lines.
	/** The first line is specified by two points start0 and end0. The second line is specified by
	two points start1 and end1.
	The implementation of this function follows http://paulbourke.net/geometry/lineline3d/ .
	\param start0 The starting point of the first line.
	\param end0 The ending point of the first line.
	\param start1 The starting point of the second line.
	\param end1 The ending point of the second line.
	\param d [out] If specified, receives the normalized distance of the closest point along the first line.
		This pointer may be left null.
	\param d2 [out] If specified, receives the normalized distance of the closest point along the second line.
		This pointer may be left null.
	\return Returns the closest point on line start0<->end0 to the second line.
	\note This is a low-level utility function. You probably want to use closestPoint() or distance() instead.
	\see closestPoint(), distance(). */

	static CVec3D closestPointLineLine(CVec3D start0, CVec3D end0, CVec3D start1, CVec3D end1, float *d, float *d2);

	/**
	\brief Returns a human-readable representation of this CLine.
	 The returned string specifies the position and direction of this CLine. */
	std::string toString() const;

#ifdef MATH_QT_INTEROP
	operator QString() const { return toString(); }
	QString toString() const { return QString::fromStdString(toString()); }
#endif
};

CLine operator *(const CMat3D &transform, const CLine &line);
CLine operator *(const CMatJoint3x4 &transform, const CLine &line);
CLine operator *(const CMat4D &transform, const CLine &line);
CLine operator *(const CQuaternion &transform, const CLine &line);

#ifdef MATH_QT_INTEROP
Q_DECLARE_METATYPE(CLine)
Q_DECLARE_METATYPE(CLine*)
#endif


std::ostream &operator <<(std::ostream &o, const CLine &line);


} //end GEO
}  //end SMF

#endif // __SMF_LINE_