#ifndef __SMF_PARABOLA_2D_
#define __SMF_PARABOLA_2D_

#include "../SMF_Config.h"
#include "../math/all.h"
#include "all.h"
#include "SMF_2DConics.h"

namespace SMF {
using namespace MATH;
using namespace Util;
namespace GEO{

class CRay;

/**
 * \class CParabola
 *
 * \ingroup SMF_Geometric
 *
 * \if pt_br
 * \brief  Modelo de curva parabólica
 * \elseif us_en
 * \brief Parabola curve model
       
       This class provides an efficient parabola curve implementation.
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CParabola : public CConicBase
    {
    public:
      /** Creates a parabola curve with given radius of curvature */
      CParabola(float roc);

      bool intersect(CVec3D &point, const CRay &ray) const;

      float sagitta(float r) const;
      float derivative(float r) const;

	};

class CParabola2D{
	public:
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
Q_DECLARE_METATYPE(CParabola2D)
Q_DECLARE_METATYPE(CParabola2D*)
#endif

namespace _2D{
std::ostream &operator <<(std::ostream &o, const CParabola2D &parabola);
} //end _2D


} //end GEO
}  //end SMF

#endif // __SMF_ELLIPSE_2D_
