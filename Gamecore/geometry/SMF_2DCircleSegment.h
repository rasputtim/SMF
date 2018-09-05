#ifndef __SMF_CIRCLE_SEGMENT_2D_
#define __SMF_CIRCLE_SEGMENT_2D_

#include "../SMF_Config.h"
#include "../math/all.h"
#include "SMF_2DPoint.h"
#include "SMF_2DCircleSector.h"
#include "SMF_GeoDefs.h"

namespace SMF {
using namespace MATH;
using namespace Util;
namespace GEO{




/**
 * \class CCircleSegment2D
 *
 * \ingroup SMF_Geometric
 * \image html pics/circle.png
 * \if pt_br
 * \brief Um Setor Circular em duas dimensões num espaço cartesiano
 * \note Esta classe representa tanto um setor sólido, quanto somente a borda do setor
 * \elseif us_en
 * \brief A two-dimensional segment in 2D space.
    \note This class represents both a hollow segment (only edge) and a solid segment.
 *
 * \endif
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CCircleSegment2D: public CCircle2D
{
public:
	/// The begin of this segment from left to right, or CW.
	CPoint2D startPoint;
	/// The end of this segment from left to right, or CW.
	CPoint2D endPoint;
	/**
	\brief sagitta at the central point of the arc
	\n \image html pics/sagitta.png
	\n Calculating the sagitta
	\n \f$ s = r - \sqrt[2]{r^2 - l^2} \f$
	\n s = sagitta
	\n r = radius
	\n l = half  length of the chord (span) connecting the two ends of the segment
	**/
	float _sag;
	/**
	\brief  The default constructor does not initialize any members of this class.
	 This means that the values of the members startPoint, endPoint and supCircle  are all undefined after creating a new Arc using
		this default constructor. Remember to assign to them before use.
	\see supCircle, endPoint, startPoint. */
	CCircleSegment2D();

	/**
	\brief  Constructs a new segment by explicitly specifying the member variables.
	\param start The start point of the segment.
	\param end The endPoint of the segment
	\param circ The supporting segment of the segment.
	\warning the start and endPoint must be on the circle edge
	*/
	CCircleSegment2D(CPoint2D &start, CPoint2D &end, CCircle2D &circ);
/**
	\brief  Constructs a new segment by explicitly specifying the member variables.
	\param start The start point of the segment.
	\param end The endPoint of the segment
	\param radius of the circle.
	*/
	CCircleSegment2D(CPoint2D &start, CPoint2D &end, float radius, int orientation = RightHandSide);

	/**
  \brief  set the circleSegment
  // This method must be called by all constructors (except the default
  // constructor) and all methods that change the coefficients.
  */
	void set(CPoint2D &start, CPoint2D &end, CCircle2D &circ);

	/**
	\brief get the sagitta at specified distance from origin (ponto méio entre o start point e o end point do segmento)
	 \note s = r - (dist+square(r^2-l^2))
	 \note  distance between [0,_sag]. if dist > _sag will return _sag. if dist < 0 will return 0;
	 where:  r=  radius of the suporting circle
	         l= half lenght of the segment (startPoint - endPoint)
	**/
	float sagitta(float dist) const;




/**
\brief return the area of the segment
 The formula to find the area of the segment is given below. It can also be found by calculating the area of the whole pie-shaped sector and subtracting the area of the isosceles triangle △ACB.
 A=R^2/2 *(((π /180)-C)-sin C)

	where:
C  is the central angle in DEGREES
R  is the radius of the circle of which the segment is a part.
π  is Pi, approximately 3.142
sin  is the trigonometry Sine function.
See Trigonometry Overview

Area of a Circular Segment given its height
Definition: The number of square units it takes to fill a segment of a circle
Try this Drag one of the orange dot that defines an endpoint of the segment. Adjust the segment height. Note the number of square units it takes to fill it and the calculation.

If you know the radius of the circle and the height of the segment, you can find the segment area from the formula below.
The result will vary from zero when the height is zero, to the full area of the circle when the height is equal to the diameter.
	A= r* Arccos((r-h)/r)-(r-h)*sqrt(2rh-h^2)
	r  is the radius of the circle of which the segment is a part.
	h  is the height of the segment (sagitta).
	**/
	float area()const;
	/**
	\brief Calculating Height of an Arc at Any Point
	\note h = s + sqrt(r^2-x^2)-r
	where:
	h = the height of the arc;
	s = the sagitta of the arc;
	r = the radius of the arc;
	x = the horizontal offset from the center to the point where you want the height;
	**/
	float height(CPoint2D &point)const;

	/**
	\brief  Returns the start point at the edge of this segment.
	*/
	SMF_INLINE CPoint2D getStartPoint() const{return startPoint;};
	/**
	\brief  Returns the end point at the edge of this segment.
	*/
	SMF_INLINE CPoint2D getEndPoint() const{return endPoint;};
	/**
	\brief  Returns the supporting segment of this segment.
	*/
	CCircle2D getSupCircle() const;

	/**
	\return the central angle of the arc in radians
	*/
	float centralAngle()const;
	/**
	\return the lenght of the arc
	\param angle in radians
	\note  	ArcLen= RC
	where:
	C  is the central angle of the arc in radians.
	R  is the radius of the arc
	*/
	float getLenght()const;


    /**
	\brief  Returns a human-readable representation of this segment. Most useful for debugging purposes.
	 The returned string specifies the center position, normal direction and the radius of this segment. */
	std::string toString() const;



#ifdef MATH_QT_INTEROP
	operator QString() const { return toString(); }
	QString toString() const { return QString::fromStdString(toString()); }
#endif
};

//CCircleSegment2D operator *(const CMat2D &transform, const CCircleSegment2D &segment);

#ifdef MATH_QT_INTEROP
Q_DECLARE_METATYPE(CCircleSegment2D)
Q_DECLARE_METATYPE(CCircleSegment2D*)
#endif

std::ostream &operator <<(std::ostream &o, const CCircleSegment2D &segment);


} //end GEO
}  //end SMF

#endif // __SMF_CIRCLE_2D_
