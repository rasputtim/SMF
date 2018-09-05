#ifndef __SMF_ELLIPSE_2D_
#define __SMF_ELLIPSE_2D_

#include "../SMF_Config.h"
#include "../math/all.h"
#include "all.h"
#include "SMF_2DConics.h"

namespace SMF {
using namespace MATH;
using namespace Util;
namespace GEO{



/**
 * \class CEllipse2D
 *
 * \ingroup SMF_Geometric
 * \image html pics\conic.png
 * \if pt_br
 * \brief Um Elipse bidimensional em espaço 2d (cartesiano)
  This class represents both a hollow ellipse (only edge) and a solid ellipse .
 An ellipse is defined by a long axis and a short axis,
 called the semi-major (a) and semi-minor axes (b), respectively.
 Usually people use the variable a to represent the length of the semi-major axis,
 and b to represent the length of the semi-minor axis
 Equations:
  Ellipse Center (h,k)
  semi-major (a)
  semi-minor axes (b)
 The general ellipse equation:
 /f$ \frac{(x-h)^2}{a^2} + \frac{(y-k)^2}{b^2} = 1 /f$

 The parametric equation of ellipse:

\f$ X(t)=h + a\,\cos t\,\cos \varphi - b\,\sin t\,\sin\varphi \f$
\f$ Y(t)=k + a\,\cos t\,\sin \varphi + b\,\sin t\,\cos\varphi \f$ 
 *
 * \elseif us_en
 *  \brief A two-dimensional ellipse in 2D space.
 This class represents both a hollow ellipse (only edge) and a solid ellipse .
 An ellipse is defined by a long axis and a short axis,
 called the semi-major (a) and semi-minor axes (b), respectively.
 Usually people use the variable a to represent the length of the semi-major axis,
 and b to represent the length of the semi-minor axis
 Equations:
  Ellipse Center (h,k)
  semi-major (a)
  semi-minor axes (b)
 The general ellipse equation:
 /f$ \frac{(x-h)^2}{a^2} + \frac{(y-k)^2}{b^2} = 1 /f$

 The parametric equation of ellipse:

\f$ X(t)=h + a\,\cos t\,\cos \varphi - b\,\sin t\,\sin\varphi \f$
\f$ Y(t)=k + a\,\cos t\,\sin \varphi + b\,\sin t\,\cos\varphi \f$ *
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CEllipse2D: public CConic2D
{

private:
	/// The radius of the ellipse on horizontal aligment. [similarOverload: center]
	///A half-axis, from the center out to the ellipse, is called a "semi-major" or a "semi-minor" axis, depending on which axis you're taking half of.
	/** This parameter must be strictly positive to specify a non-degenerate ellipse. If zero is specified, this ellipse
		is considered to be degenerate.
		\see CEllipse2D::CEllipse2D(). */
	float rx;
	/// The radius of the ellipse on vertical alignment. [similarOverload: center]
	///A half-axis, from the center out to the ellipse, is called a "semi-major" or a "semi-minor" axis, depending on which axis you're taking half of.
	/** This parameter must be strictly positive to specify a non-degenerate ellipse. If zero is specified, this ellipse
		is considered to be degenerate.
		\see CEllipse2D::CEllipse2D(). */
	float ry;
public:

	/// The center position of this ellipse.
	CPoint2D center;

	/**
	\brief ratio used ot speedup calcs on draw operations
	**/
	float xy_ratio;

	/** The default constructor does not initialize any members of this class.
	 This means that the values of the members center, phi and r are all undefined after creating a new ellipse using
		this default constructor. Remember to assign to them before use.
		\see center, Phi, rx, ry. */
	CEllipse2D();

	/** Constructs a new ellipse by explicitly specifying the member variables.
	    \param center The center point of the ellipse.
		\param rx The horizontal half-axis of the ellipse.
		\param ry The vertical half-axis of the ellipse.
		\param angle angle of rotation from the horizontal axis
		\see center, Phi, rx, ry. */
	CEllipse2D(CPoint2D center,  float rx,  float ry, float angle=0.0f);

	//\brief return (a) - The half major axis  of the ellipse
	// onde C(h,k) é o centro do elipse
	float a() { 	float a= (rx > ry ? rx: ry); return a;};
	float a()const { 	float a= (rx > ry ? rx: ry);return a;};

	//\brief return (b) - The half minor axis  of the ellipse
	float b(){  float b= (rx > ry ? ry: rx); return b;};
	float b()const {  float b= (rx > ry ? ry: rx); return b;};
	//\brief return (c) - The half distance between the focci
	float c(){  float c=CMath::sqrt(CMath::square(a())-CMath::square(b())); return c;};
	float c()const {  float c=CMath::sqrt(CMath::square(a())-CMath::square(b()));return c;};
	/**
	\brief return (g) - The flattering fattor
	\note \f$ g = 1 - b/a \f$
	**/
	float g(){  float g= 1 - (b()/a()); return g;};
	float g()const {  float g= 1 - (b()/a()); return g;};
	/**
	\brief return area of the elipse
	**/
	float getArea()const{float area= CMath::PI*a()*b(); return area;};
	/**
	\brief return directrix of the elipse from Focus 1
	\note Directrix of a conic section is a line such that ratio of the distance of the points on the conic
	section from focus to its distance from directrix is constant.
	\this assume the angle Phi is 0
	\todo consider the angle of the ellipse
	\note direc= \f$ (x) = a^2/c where c= half distance between focci and a=major half-axis**/
	CLine2D getDirectrixF1()const{return CLine2D(c(),0,-CMath::square(a()));};
	/**
	\brief return directrix of the elipse from Focus 2
	\note Directrix of a conic section is a line such that ratio of the distance of the points on the conic
	section from focus to its distance from directrix is constant.
	\this assume the angle Phi is 0
	\todo consider the angle of the ellipse
	\note direc= \f$ (x) = -a^2/c where c= half distance between focci and a=major half-axis**/
	CLine2D getDirectrixF2()const{return CLine2D(c(),0,CMath::square(a()));};

	/// Returns the center point of this ellipse.
	/** This point is also the center of mass for this ellipse. The functions centerPoint() and centroid() are equivalent.
		\see center. */
	CPoint2D centerPoint() const { return center; }
	CPoint2D centroid() const { return center; } ///< [similarOverload: centerPoint]
	/**
	\brief return ellipse excentricity
	\note For an ellipse, eccentricity is: c/a onde c= distância do centro ao foco e a= half-major
	e = Sqroot(1 - (b^2-a^2))
	\note that r/d=e
	where r is the distance from the focus to any point M(x,y) of an ellipse.
	d the distance from M(x,y) to the directrix,
	and e is the eccentricity.
	**/
	float getEccentricity()const;

	//return the focci length (the distance from focus1 t focus2)
	// focci = 2c. and c = Sqroot(a^2-b^2)
    float getFocci() const;


	void set(CPoint2D center, float rx, float ry,float anglePhi=0.0f);


	/** Returns a point at the edge of this ellipse.
	 \param angleRadians The direction of the point to get. A full ellipse is generated by the range [0, 2*pi],
			but it is ok to pass in values outside this range.
		\return A point in world space at the edge of this ellipse.
		\warning assume the elipse is Phi=0. touse this method you mus transform the ellipse to angle 0
		\Todo Consider the ellipse angle or transform the ellipse
	*/
	CPoint2D getPoint(float angleRadians) const;
	/** Returns a point at the edge of this ellipse.
	 \param direction Position Vector pointing to the direction . A full ellipse is generated by the range [0, 2*pi],
			but it is ok to pass in values outside this range.
		\return A point in world space at the edge of this ellipse.
		\warning assume the elipse is Phi=0. touse this method you mus transform the ellipse to angle 0
		\Todo Consider the ellipse angle or transform the ellipse
	*/
	CPoint2D getPoint(CVec2D direction) const;

	/** Returns a point inside this ellipse.
	 \param angleRadians The direction of the point to get. A full ellipse is generated by the range [0, 2*pi],
			but it is ok to pass in values outside this range.
	 \param dist A value in the range [0,1] that specifies the normalzied distance of the point from the center of the ellipse.
			A value of 0 returns the center point of this ellipse. A value of 1 returns a point at the edge of this ellipse.
			The range of d is not restricted, so it is ok to pass in values larger than 1 to generate a point lying completely
			outside this ellipse.
	\warning assume the elipse is Phi=0. touse this method you mus transform the ellipse to angle 0
	\Todo Consider the ellipse angle or transform the ellipse
	*/
	CPoint2D getPoint(float angleRadians, float dist) const;
	/** Returns a point inside this ellipse.
	 \param direction a 2D Position Vector that point to the direction. A full ellipse is generated by the range [0, 2*pi],
			but it is ok to pass in values outside this range.
	 \param dist A value in the range [0,1] that specifies the normalzied distance of the point from the center of the ellipse.
			A value of 0 returns the center point of this ellipse. A value of 1 returns a point at the edge of this ellipse.
			The range of d is not restricted, so it is ok to pass in values larger than 1 to generate a point lying completely
			outside this ellipse.
	\warning assume the elipse is Phi=0. touse this method you mus transform the ellipse to angle 0
	\Todo Consider the ellipse angle or transform the ellipse
	*/
	CPoint2D getPoint(CVec2D direction, float dist) const;

	/** \brief Get distance between origin and farthest edge in specified direction
	// r=ab / sqrt(a2sin2θ+b2cos2θ)  θ is the angle from major axis
	// \param angle in radians
	\warning assume the elipse is Phi=0. touse this method you mus transform the ellipse to angle 0
	\Todo Consider the ellipse angle or transform the ellipse
	**/
	double getRadius(const  float angle) const;
	/** \brief Get distance between origin and farthest edge in specified direction
	\warning assume the elipse is Phi=0. touse this method you mus transform the ellipse to angle 0
	\Todo Consider the ellipse angle or transform the ellipse
	**/
	double getRadius(const CVec2D dir) const;
	/** Return one of the focus of this elipse
	\warning assume the elipse is Phi=0. touse this method you mus transform the ellipse to angle 0
	\Todo Consider the ellipse angle or transform the ellipse
	**/
	CPoint2D getFocus1()const;
	/** Return one of the focus of this elipse
	\warning assume the elipse is Phi=0. touse this method you mus transform the ellipse to angle 0
	\Todo Consider the ellipse angle or transform the ellipse
	**/
	CPoint2D getFocus2()const;
	/**
	\brief retorna o valor de Y1 dado x, utilizando a formula da elipse;
	\warning assume the elipse is Phi=0. touse this method you mus transform the ellipse to angle 0
	\Todo Consider the ellipse angle or transform the ellipse
	**/

	float getY1(float x)const;
	/**
	\brief retorna o valor de Y2 dado x, utilizando a formula da elipse;
	\warning assume the elipse is Phi=0. touse this method you mus transform the ellipse to angle 0
	\Todo Consider the ellipse angle or transform the ellipse
	**/
	float getY2(float x)const;

	/// Translates this CEllipse2D in world space.
	/** \param offset The amount of displacement to apply to this CEllipse2D, in world space coordinates.
		\see transform(). */
	void translate(const CVec2D &offset);

	/// Applies a transformation to this CEllipse2D.
	/** \param transform The transformation to apply to this CEllipse2D. This transformation must be
		affine, and must contain an orthogonal set of column vectors (may not contain shear or projection).
		The transformation can only contain uniform scale, and may not contain mirroring.
		\see translate(), scale(), classes CMat3D, CMatJoint3x4, CMat4D, CQuaternion. */
	void transform(const CMat2D &transform);

	/// Tests if the given point is contained at the edge of this ellipse.
	/** \param point The target point to test.
		\param maxDistance The epsilon threshold to test the distance against.
		\see distanceToEdge(), distanceToDisc(), closestPointToEdge(), closestPointToDisc().
		\todo Implement DiscContains(CVec2D/CLineSegment/CTriangle). */
	bool edgeContains(const CPoint2D &point, float maxDistance = 1e-6f) const;


	/// Computes the distance of the given object to the edge of this ellipse.
	/** \todo Implement distanceToEdge(CRay/CLineSegment/CLine).
		\return The distance of the given point to the edge of this ellipse. If the point is contained on this ellipse,
			the value 0 is returned.
		\see distanceToDisc(), closestPointToEdge(), closestPointToDisc().
	\warning since use getradius, this assume the elipse is Phi=0. touse this method you mus transform the ellipse to angle 0
	\Todo Consider the ellipse angle or transform the ellipse

		*/
	float distanceToEdge(const CPoint2D &point) const;
//	float distanceToEdge(const CRay2D &ray, float *d, CVec2D *closestPoint) const;
//	float distanceToEdge(const CLineSegment2D &lineSegment, float *d, CVec2D *closestPoint) const;
//	float distanceToEdge(const CLine2D &line, float *d, CVec2D *closestPoint) const;


	/// Tests this ellipse for an intersection against the given object.
	/** \see intersects().
		\warning assume the elipse is Phi=0. touse this method you mus transform the ellipse to angle 0
	\Todo Consider the ellipse angle or transform the ellipse
    */
	bool intersectsDisc(const CLine2D &line) const;
	bool intersectsDisc(const CLineSegment2D &lineSegment) const;
//	bool intersectsDisc(const CRay2D &ray) const;

	/**
	\brief calculate the intersection points (if they exists) from a line and a ellipse;
	\return NOT_INTERESECTING if there are no intersection between the line nd the ellipse
	\return INTERESECTING if there are two intsection points
	\return INTERESECTING_EXTREMITY_P1 if there is only one intersecting point
	\note
	\see http://csharphelper.com/blog/2012/09/calculate-where-a-line-segment-and-an-ellipse-intersect-in-c/
	\note i evolved this method to use it on the ellipse with center not in the origin
	\warning assume the elipse is Phi=0. touse this method you mus transform the ellipse to angle 0
	\Todo Consider the ellipse angle or transform the ellipse

	**/
static IntersectResult	intersection(const CEllipse2D & elipse,const CLine2D &line,CPoint2D &int_point1,CPoint2D &int_point2);
		/** This function returns true if the given object lies inside this Circle, and false otherwise.
		\note The comparison is performed using less-or-equal, so the faces of this CAABBox count as being inside, but
			due to float inaccuracies, this cannot generally be relied upon.
		\see distance(), intersects(), closestPoint(). */
	bool contains(const CVec2D &point) const;
	bool contains(const CPoint2D &point) const;

    /// Returns a human-readable representation of this ellipse. Most useful for debugging purposes.
	/** The returned string specifies the center position, normal direction and the radius of this ellipse. */
	std::string toString() const;



#ifdef MATH_QT_INTEROP
	operator QString() const { return toString(); }
	QString toString() const { return QString::fromStdString(toString()); }
#endif
};

//CEllipse2D operator *(const CMat2D &transform, const CEllipse2D &ellipse);

#ifdef MATH_QT_INTEROP
Q_DECLARE_METATYPE(CEllipse2D)
Q_DECLARE_METATYPE(CEllipse2D*)
#endif

namespace _2D{
std::ostream &operator <<(std::ostream &o, const CEllipse2D &ellipse);
} //end _2D


} //end GEO
}  //end SMF

#endif // __SMF_ELLIPSE_2D_
