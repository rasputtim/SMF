#ifndef __SMF_CIRCLE_ARC_2D_
#define __SMF_CIRCLE_ARC_2D_

#include "../SMF_Config.h"
#include "../math/all.h"
#include "SMF_2DPoint.h"
#include "SMF_2DCircle.h"

namespace SMF {
using namespace MATH;
using namespace Util;
namespace GEO{


 /**
 * \class CCircSector2D
 *
 * \ingroup SMF_Geometric
 * \image html pics/circle.png
 * \if pt_br
 * \brief Um Setor Circular em duas dimensões num espaço cartesiano
 * \note Esta classe representa tanto um setor sólido, quanto somente a borda do setor
 * \elseif us_en
 * \brief A two-dimensional arc in 2D space.
    \note This class represents both a hollow arc (only edge) and a solid arc (sector).
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
class SMF_API CCircSector2D : public CCircle2D
{
public:
	/// The begin of this arc from left to right, or CW.
	CPoint2D startPoint;
	/// The end of this arc from left to right, or CW.
	CPoint2D endPoint;


	/**
	\brief  The default constructor does not initialize any members of this class.
	 This means that the values of the members startPoint, endPoint and supCircle  are all undefined after creating a new Arc using
		this default constructor. Remember to assign to them before use.
	\see supCircle, endPoint, startPoint. */
	CCircSector2D();

	/**
	\brief  Constructs a new arc by explicitly specifying the member variables.
	\param start The start point of the arc.
	\param end The endPoint of the arc
	\param circ The supporting arc of the arc.
	*/
	CCircSector2D(CPoint2D &start, CPoint2D &end, CCircle2D &circ);

	/**
	\brief  Constructs a new segment by explicitly specifying the member variables.
	\param start The start point of the segment.
	\param end The endPoint of the segment
	\param radius of the circle.
	*/
	CCircSector2D(CPoint2D &start, CPoint2D &end, float radius, int orientation = RightHandSide);

	/**
  \brief  set the circleSegment
  // This method must be called by all constructors (except the default
  // constructor) and all methods that change the coefficients.
  */
	void set(CPoint2D &start, CPoint2D &end, CCircle2D &circ);
	/** Get curve sagitta at specified distance from origin.
         \param dist distance from supCircle origin */
    float sagitta(float dist) const;
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
	\brief  Returns the start point at the edge of this arc.
	*/
	SMF_INLINE CPoint2D getStartPoint() const{return startPoint;};
	/**
	\brief  Returns the end point at the edge of this arc.
	*/
	SMF_INLINE CPoint2D getEndPoint() const{return endPoint;};
	/**
	\brief  Returns the supporting arc of this arc.
	*/
	SMF_INLINE CCircle2D getSupCircle() const;


	/** This function returns true if the given object lies inside this arc, and false otherwise.
	\note The comparison is performed using less-or-equal, so the faces of this CAABBox count as being inside, but
		due to float inaccuracies, this cannot generally be relied upon.
	\see distance(), intersects(), closestPoint(). */
	bool contains(const CVec2D &point) const;
	bool contains(const CPoint2D &point) const;
	bool isInsideSector(const CPoint2D &point, const CPoint2D &center,  const CPoint2D &sectorStart, const CPoint2D &sectorEnd)const;

    /**
	\brief  Returns a human-readable representation of this arc. Most useful for debugging purposes.
	 The returned string specifies the center position, normal direction and the radius of this arc. */
	std::string toString() const;



#ifdef MATH_QT_INTEROP
	operator QString() const { return toString(); }
	QString toString() const { return QString::fromStdString(toString()); }
#endif
};

//CCircSector2D operator *(const CMat2D &transform, const CCircSector2D &arc);

#ifdef MATH_QT_INTEROP
Q_DECLARE_METATYPE(CCircSector2D)
Q_DECLARE_METATYPE(CCircSector2D*)
#endif

std::ostream &operator <<(std::ostream &o, const CCircSector2D &arc);


} //end GEO
}  //end SMF

#endif // __SMF_CIRCLE_2D_
