#ifndef __SMF_CONIC_2D_
#define __SMF_CONIC_2D_

#include "../SMF_Config.h"
#include "SMF_GeoDefs.h"
#include "../math/all.h"
#include "SMF_2DPoint.h"

namespace SMF {
using namespace MATH;
using namespace Util;
namespace GEO{

class CRay;


/**
 * \class CCurveBase
 *
 * \ingroup SMF_Geometric
 * \image html pics\curve.png
 * \if pt_br
 * \brief Classe Base para modelos de curvas
 * \note Fornece suporte para Curvas em espaço 2D e 3D
 * \elseif us_en
 * \brief CCurveBase class for surface curvature models
       
       This class defines an interface for surface curvature
       implementations. Curvature is defined as a surface curve in
       three dimensional space. 

       It provides access to sagitta (z) and gradient data on any
       curved surface point (x, y). Ability to find point of
       intersection between a given 3d ray and the curve is also
       provided.
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
    class SMF_API CCurveBase 
    {
    public:
      virtual inline ~CCurveBase();

      /** 
	  \brief Get curve sagitta at specified point 
	  \image html pics\curve.png
	  \note the sagitta of a circular arc is the distance from the center of the arc to the center of its base
	  \see http://en.wikipedia.org/wiki/Sagitta_%28geometry%29 
	  \note 
	  The famous astronomer Karl Schwarzschild found a constant K in order to make accurate squashed or elongated ellipsoid, paraboloid and even hyperboloid Cassegrain telescope mirrors, whose primary focus point (f) is half of the r or r' radius. Schwarzschild's constant is 0 for a sphere and –1 for a parabola. It is higher than 0 (positive) for squashed (oblate) ellipses and lower than –1 for hyperbolas.

		Any value between 0 and –1 for K  = –beta ^ 2  =  g ^ 2 – 1  produces elongated (prolate) ellipses. Using Schwarzschild's formula, please note that r stands for the smaller circle of reference radius for parabolas and prolate ellipses. It rather stands for the larger circle of reference radius for oblate ellipses.

		However, it turns out that Schwarzschild's constant is interconnected with the r radius. It is redundant, but one may prefer to use it anyway in order to avoid the conversion to R while working on telescope mirrors. Otherwise, the simplified sagitta formula below for the ellipse rather uses a radius R = x / 2 which is that of the ellipse's larger radius (or that of its circumscribed circle).

		The sagitta according to Schwarzschild is given by:

		sag = h ^ 2 / (r * (1 + sqr(1 – (h ^ 2 / r ^ 2) * (K + 1))))

		For the circle:  sag = R – sqr(R ^ 2 – h ^ 2) )

		For the parabola:  sag = h ^ 2 / 2 r

		Simplified formula for the ellipse:  sag = R – sqr(R 2 – (h / g) 2 ) = .2

		It should be emphasized that there is no such thing as prolate or oblate ellipse. Those words are rather related to a specific portion of the ellipse, so that there is no true sagitta for the so-called oblate ellipse. The equivalent orthogonal sagitta is more exactly the prolongation of the h (height) distance and it is given by a calculus based on the ellipse's true sagitta (h' = R – sag). However, as compared to that of a circle (see formula above) it is contracted according to g:

		Pseudo-sagitta for the oblate ellipse:  sag' = g * (R – sqr(R ^ 2 – (R – sag) ^ 2) ) = .2
	  
	  
	  */
      virtual float sagitta(const CVec2D & xy) const = 0;

      /** \brief  Get curve x and y derivative (gradient) at specified point */
      virtual void derivative(const CVec2D & xy, CVec2D & dxdy) const;

      /** \brief  Get intersection point between curve and 3d ray. Return
          false if no intersection occurred */
      virtual bool intersect(CVec3D &point, const CRay &ray) const;

      /** \brief  Get normal to curve surface at specified point */
      virtual void normal(CVec3D &normal, const CVec3D &point) const;
    };

/**
 * \class CRotationalCurv
 *
 * \ingroup SMF_Geometric
 * \image html pics\curve.png
 * \if pt_br
 * \brief Extenção da Classe CCurveBase para modelos de curvas rotacionais simétricas
 *        Define a interface para curvas rotacionais simétricas e também fornece
 *		  implementação padrão de métodos como curva genérica
 * \elseif us_en
 * \brief CCurveBase class for rotationally symmetric curves.
 *      This class defines rotationally symmetric curve interface and
 *      provide default implementation as generic non symmetric curve.
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
class SMF_API CRotationalCurv : public CCurveBase
    {
    public:
      CRotationalCurv();

      virtual void normal(CVec3D &normal, const CVec3D &point) const;

      /** Get curve sagitta at specified distance from origin.
          \param dist distance from curve origin (0, 0)
      */
      virtual float sagitta(float dist) const = 0;

      /** Get curve derivative at specified distance from origin.
          \param dist distance from curve origin (0, 0)
      */
      virtual float derivative(float dist) const;
	  /**Get curve sagitta at specified point
	  \warning this method supose the curve is at the origin (0,0) and xy is a positional vector
	   **/
      inline float sagitta(const CVec2D & xy) const;
      void derivative(const CVec2D & xy, CVec2D & dxdy) const;
      /// copy constructor
	  CRotationalCurv(CRotationalCurv const& c)
		: deriv_function(c.deriv_function){}
	  /// assignment operator
	  CRotationalCurv& operator=(CRotationalCurv const& c) {
		deriv_function=c.deriv_function; 
		return *this;
	  }
    private:
      static float deriv_function_sagitta(float x, void *params);
      CNumDEqu deriv_function;
    };

   
/**
 * \class CCurveRoc
 *
 * \ingroup SMF_Geometric
 * \image html pics\curve.png
 * \if pt_br
 * \brief Extenção da Classe CCurveBase para curvas com raio de curvatura
 * \elseif us_en
 * \brief CCurveBase class for curves with a radius of curvature 
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
class SMF_API CCurveRoc
    {
    public:
      /** Set the radius of curvature */
      inline void set_roc(float roc);

      /** Get the radius of curvature */
      inline float get_roc() const;

    protected:
      inline CCurveRoc(float roc);
	  //Remember to set theparameters of the conic
	  inline CCurveRoc(){};
	  /// copy constructor
	  CCurveRoc(CCurveRoc const& c)
		: _roc(c._roc){}
	  /// assignment operator
	  CCurveRoc& operator=(CCurveRoc const& c) {
		_roc=c._roc; 
		return *this;
	  }
	  /** curvature radius on x=0*/
      float _roc;

    };

   /**
   \brief CAxisRot class for curves with a angle of rotation 
   */
    class SMF_API CAxisRot
    {
    public:
      /** Set the angle of rotation in radians
	  \note range [0,2PI]
	  \note angle in CCW direction 
	  */
		inline void set_ang(float ang){ _phi = ang;};

      /** Get the angle of inclination */
      inline float get_ang() const { return _phi;};
	  inline float Phi() const { return _phi;}

    protected:
      inline CAxisRot(float ang);
	  //Remember to set theparameters of the conic
	  inline CAxisRot(){};
	   /// copy constructor
	  CAxisRot(CAxisRot const& c)
		: _phi(c._phi){}
	  /// assignment operator
	  CAxisRot& operator=(CAxisRot const& c) {
		_phi=c._phi; 
		return *this;
	  }
	  /** \brief rx half-axix rotation angle - the angle between rx and x-axis in radians
	  angle in radians of the rx axis to the x axis on the carthesian space
	 rx could be the major or minor axis
	 \note angle θ such that \f$ cot(2θ) = (A - C )/ B \f$
	\warning Phi range [0,2*PI]
	**/
	float _phi;

    };

	
/**
 * \class CConicBase
 *
 * \ingroup SMF_Geometric
 * \image html pics\conic.png
 * \if pt_br
 * \brief Extenção da Classe CCurveBase para curvas cônicas
          Define propriedades comuns de curvas rotacionais simétricas e cônicas
		  Estas curvas são todas definidas porum raio de curvatura e corficientes de deformação.
		  Ajuste pode ser utilizado para encontrar a cônica que melhor se ajusta de uma outra cuva 
		  simétrica rotacional com parÂmetro de deformação fixo ou livre.
	\see http://pt.wikipedia.org/wiki/Ajuste_de_curvas
 * \elseif us_en
 * \brief CCurveBase class for conic family of curves
       
       This base class defines common properties of rotationally
       symmetric conic curves. These curves are all defined by a
       radius of curvature and deformation coefficient.

       Fitting can be used to find best fit conic of an other
       rotationally symmetric curve either with fixed or free
       deformation parameter.
 * \see http://en.wikipedia.org/wiki/Curve_fitting
 *
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
 class SMF_API CConicBase : public CRotationalCurv, public CCurveRoc, public CAxisRot
    {
    public:
    
	enum conicType_t {
    no_type=0,
    real_ellipse,
    real_circle,
    imaginary_ellipse,
    imaginary_circle,
    hyperbola,
    parabola,
    real_intersecting_lines,
    complex_intersecting_lines,
    real_parallel_lines,
    complex_parallel_lines,
    coincident_lines,
    num_conic_types // is here to enable iterating through this list
  };
		
		/** get Schwarzschild constant
	  \note  conic constant K=-e^2
	  */
      float get_schwarzschild() const;

      /** get eccentricity 
	  e= -Sqr(K)*/
      float get_eccentricity() const;

      /** Adjust radius of curvature to best fit given
          curve. Deformation Won't be changed by this function. See
          Conic::fit() to adjust deformation too.

          \param curve Curve to fit
          \param radius Maximum radius used to get sample points
          \param count Number of sample points to use
      */
      float fit_roc(const CRotationalCurv &curve, float radius, unsigned int count);

	  SMF_INLINE conicType_t type() const { return type_; }

  /** 
  \brief Returns the type of the conic as a string.
  // Possible returned strings are:
  // "real ellipse", "real circle", "imaginary ellipse", "imaginary circle",
  // "hyperbola", "parabola",
  // "real intersecting lines", "complex intersecting lines",
  // "real parallel lines", "complex parallel lines", "coincident lines".
  // The default constructor sets the type to "invalid conic".
  */
  CMyString real_type() const;

  /** 
  \brief  Returns the internal enum value corresponding to the string argument.
  // Useful for comparison purposes, or for use in "case" statements.
  */
  static conicType_t type_by_name(CMyString const& name);

  /**  
  \brief  Converts the conic type from enum (internal representation) to string.
  */
  static CMyString type_by_number(conicType_t type);


      virtual bool intersect(CVec3D &point, const CRay &ray) const = 0;
      virtual float sagitta(float r) const = 0;
      virtual float derivative(float r) const = 0;
  /// copy constructor
  CConicBase(CConicBase const& c)
    : CCurveRoc(c),CAxisRot(c),CRotationalCurv(c),type_(c.type()), _sh(c._sh){}
  /// assignment operator
  CConicBase& operator=(CConicBase const& c) {
	  CCurveRoc::operator= (c);
	  CAxisRot::operator= (c);
	  CRotationalCurv::operator= (c);
    type_=c.type(); _sh=c._sh; 
    return *this;
  }

    protected:
    conicType_t type_;
  
		inline CConicBase(float roc, float sc);
	  //Remember to set theparameters of the conic
	  inline CConicBase() : type_(no_type) {};
	  float _sh;       // Schwarzschild constant + 1
    };



/**
 * \class CConic
 *
 * \ingroup SMF_Geometric
 * \image html pics\conic.png
 * \if pt_br
 * \brief Modelo de curva cônica de propósito geral
        This class models a rotationally symmetric conic curves with
       given radius of curvature and deformation coefficient. The
       later can be provided either as Schwarzschild constant or
       Eccentricity value.

       Fitting can be used to find best fit conic of an other
       rotationally symmetric curve either with fixed or free
       deformation parameter.

       \ref Parabola offer optimized implementations
       for common special cases.
    \note ellipse and hyperbola equation standard forms:

      ellipse:    \f$ (x^2)/(a^2) + (y^2)/(b^2) = 1 \f$ 
      hyperbola:  \f$ (x^2)/(a^2) - (y^2)/(b^2) = 1 \f$ 

      with bend point at (0,0):

      ellipse:    \f$ ((x-a)^2)/(a^2) + (y^2)/(b^2) = 1 \f$ 
      hyperbola:  \f$ ((x-a)^2)/(a^2) - (y^2)/(b^2) = 1 \f$ 

      and eccentricity:

      ellipse:   \f$  e = CMath::sqrt(1 - (b^2) / (a^2)) \f$ 
      hyperbola: \f$  e = CMath::sqrt(1 + (b^2) / (a^2)) \f$ 

      both can be rewritten as:

      \f$  y^2 = (e^2 - 1) * x^2 - 2 * a * (e^2 - 1) * x \f$ 

      Best fit conic through (0,0)
      ============================

      \f$  y^2 / x = (e^2 - 1) * x - 2 * a * (e^2 - 1) \f$ 

      is a line model of this form:

      \f$  Y = C0 + C1 * X \f$ 

      with

      \f$  Y = y^2 / x \f$ 
      \f$  X = x \f$ 
      \f$  C0 = - 2 * a * (e^2 - 1) \f$ 
      \f$  C1 = (e^2 - 1) \f$ 

      C0 and C1 are found by least squares fit of points (x,y)

      eccentricity can then be computed for all conic sections with:

       \f$ e = \sqrt[2]{C1 + 1} \f$ 
       \f$ sc = -C1 - 1 \f$ 

      radius of curvature can be computed:

       \f$ a = C0 / (-2 * C1) \f$ 
       \f$ roc = a * (1 + sc) \f$ 

      once simplified, works for all conic sections:

       \f$ roc = C0 / 2  \f$ 
 * \elseif us_en
 * \brief General purpose conic curve model
       
       This class models a rotationally symmetric conic curves with
       given radius of curvature and deformation coefficient. The
       later can be provided either as Schwarzschild constant or
       Eccentricity value.

       Fitting can be used to find best fit conic of an other
       rotationally symmetric curve either with fixed or free
       deformation parameter.

       \ref Parabola offer optimized implementations
       for common special cases.
    \note ellipse and hyperbola equation standard forms:

      ellipse:    \f$ (x^2)/(a^2) + (y^2)/(b^2) = 1 \f$ 
      hyperbola:  \f$ (x^2)/(a^2) - (y^2)/(b^2) = 1 \f$ 

      with bend point at (0,0):

      ellipse:    \f$ ((x-a)^2)/(a^2) + (y^2)/(b^2) = 1 \f$ 
      hyperbola:  \f$ ((x-a)^2)/(a^2) - (y^2)/(b^2) = 1 \f$ 

      and eccentricity:

      ellipse:   \f$  e = \sqrt[2]{1 - (b^2) / (a^2)} \f$ 
      hyperbola: \f$  e = \sqrt[2]{1 + (b^2) / (a^2)} \f$ 

      both can be rewritten as:

      \f$  y^2 = (e^2 - 1) * x^2 - 2 * a * (e^2 - 1) * x \f$ 

      Best fit conic through (0,0)
      ============================

      \f$  y^2 / x = (e^2 - 1) * x - 2 * a * (e^2 - 1) \f$ 

      is a line model of this form:

      \f$  Y = C0 + C1 * X \f$ 

      with

      \f$  Y = y^2 / x \f$ 
      \f$  X = x \f$ 
      \f$  C0 = - 2 * a * (e^2 - 1) \f$ 
      \f$  C1 = (e^2 - 1) \f$ 

      C0 and C1 are found by least squares fit of points (x,y)

      eccentricity can then be computed for all conic sections with:

       \f$ e = \sqrt[2]{C1 + 1} \f$ 
       \f$ sc = -C1 - 1 \f$ 

      radius of curvature can be computed:

       \f$ a = C0 / (-2 * C1) \f$ 
       \f$ roc = a * (1 + sc) \f$ 

      once simplified, works for all conic sections:

       \f$ roc = C0 / 2  \f$ 
 *
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
	class SMF_API CConic : public CConicBase
    {
    public:
		//Remember to set the parameters of the conic
		CConic(){};
		/** Creates a conic curve with given radius of curvature and
          Schwarzschild constant */
		CConic(float roc, float sc);
		/// copy constructor
		CConic(CConic const& c)
		: CConicBase(c) {}
		/// assignment operator
		CConic& operator=(CConic const& c) {
		//base class assignment
		CConicBase::operator=(c);
		return *this;
		}
	  
      /** Set Conic From Schwarzschild constant */
      inline void set_schwarzschild(float sc);

      /** Set Conic From eccentricity */
      inline void set_eccentricity(float e);

      /** Adjust radius of curvature _and_ deformation to best fit given curve

          \param curve Curve to fit
          \param radius Maximum radius used to get sample points
          \param count Number of sample points to use
		*/
		float fit(const CRotationalCurv &curve, float radius, unsigned int count);

		bool intersect(CVec3D &point, const CRay &ray) const;
		/** \brief Get curve sagitta at specified distance from origin.
		\param dist distance from curve origin (0, 0)
		**/
		float sagitta(float dist) const;
		/**
		\brief Get curve derivative at specified distance from origin.
		\param dist distance from curve origin (0, 0)
		*/
       float derivative(float dist) const;

    };


class CCircle2D;
class CPoint2D;


/**
 * \class CConic2D
 *
 * \ingroup SMF_Geometric
 * \image html pics\conic.png
 * \if pt_br
 * \brief Um Curva plana quadrática em espaço 2d(cartesiano)
Conics are degree two curves because their most general form is the following degree two implicit polynomial:

	\f$ Ax^2 + 2Bxy + Cy^2 + 2Dx + 2Ey + F = 0 \f$

	In the above polynomial, the coefficients of xy, x and y are 2B, 2D and 2E, respectively. This polynomial has six coefficients; 
	however, dividing it with a non-zero coefficient would reduce six to five. Thus, in general, five conditions can uniquely determine a conic. 
	In linear algebra, you perhaps have learned the way of reducing the above polynomial to a normal form using eigenvalues and eigenvectors.

	Frequently, we only want to know the curve type of a general second degree polynomial. 
	In this case, as long as the second degree equation represents a conic rather than two intersecting or parallel lines, it can easily be done as follows:

    If B^2 < A*C, the general equation represents an ellipse.
    IF B^2 = A*C, the general equation represents a parabola.
    If B^2 > A*C, the general equation represents a hyperbola. 

Expression B^2-A*C is called the discriminant of the general second degree polynomial. Based on the above, if the value of the discriminant is less than, equal to or greater than zero, the conic is an ellipse, a parabola, or a hyperbola.
Conics in Matrix Form
One nice thing of conics is that its general form can be rewritten compactly using matrices. 
First, each point x = (x, y) is considered as a column vector whose third component is 1 and hence the transpose is a row vector, written as xT = [ x, y, 1 ]. 
Next, the six coefficients of the general second degree polynomial are used to construct a three-by-three symmetric matrix as follows:


It is not difficult to verify that the general second degree polynomial becomes

Now what you have learned from linear algebra can be applied to this matrix form. 

//  A conic is either an ellipse (or circle), a hyperbola, or a parabola.
//  It is represented by a quadratic equation in two nonhomogeneous
//  or three homogeneous coordinates.  Conversely, every quadratic
//  equation represents a conic, be it that it can be degenerate:
//  either in two (intersecting or parallel) lines, or in two
//  coincident lines.  Also, it can have no "visible", real points,
//  when it is an imaginary ellipse, or consist of two complementary
//  imaginary lines in which case it only has one real point, which could
//  still be at infinity.
//
//  These 11 cases are the possible values of CConic2D::real_type().
//  The default constructor sets the type to "invalid conic";
//  otherwise the correct type is automatically set when the equation
//  of the conic is given to the constructor that takes 6 numeric values
//  (a,b,c,d,e,f): the cartesian equation is then
//  \f$ Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0  \f$ 
//
//  When the conic is degenerate and consists of two lines, the method
//  components() returns a list of two (possibly identical) lines.
//  Otherwise, this method returns an empty list.
// \see http://en.wikipedia.org/wiki/Conic_section#Discriminant_classification
// \see http://en.wikipedia.org/wiki/Matrix_representation_of_conic_sections
 *
 * \elseif us_en
 * \ brief A quadratic plane curve

//  This example tells you the type of the given conic equation,
//  and prints the equation in readable form:
// \code
//   CConic2D<float> c(1, 0, 2, 0, 0, -3);
//   cout << c.real_type() << '\n'; // prints "real ellipse"
//   cout << c << '\n'; // prints the equation: X^2 + 2 Y^2 - 3 = 0
// \endcode

Conics are degree two curves because their most general form is the following degree two implicit polynomial:

	\f$ Ax^2 + 2Bxy + Cy^2 + 2Dx + 2Ey + F = 0 \f$

	In the above polynomial, the coefficients of xy, x and y are 2B, 2D and 2E, respectively. This polynomial has six coefficients; 
	however, dividing it with a non-zero coefficient would reduce six to five. Thus, in general, five conditions can uniquely determine a conic. 
	In linear algebra, you perhaps have learned the way of reducing the above polynomial to a normal form using eigenvalues and eigenvectors.

	Frequently, we only want to know the curve type of a general second degree polynomial. 
	In this case, as long as the second degree equation represents a conic rather than two intersecting or parallel lines, it can easily be done as follows:

    If B^2 < A*C, the general equation represents an ellipse.
    IF B^2 = A*C, the general equation represents a parabola.
    If B^2 > A*C, the general equation represents a hyperbola. 

Expression B^2-A*C is called the discriminant of the general second degree polynomial. Based on the above, if the value of the discriminant is less than, equal to or greater than zero, the conic is an ellipse, a parabola, or a hyperbola.
Conics in Matrix Form
One nice thing of conics is that its general form can be rewritten compactly using matrices. 
First, each point x = (x, y) is considered as a column vector whose third component is 1 and hence the transpose is a row vector, written as xT = [ x, y, 1 ]. 
Next, the six coefficients of the general second degree polynomial are used to construct a three-by-three symmetric matrix as follows:


It is not difficult to verify that the general second degree polynomial becomes

Now what you have learned from linear algebra can be applied to this matrix form. 

//  A conic is either an ellipse (or circle), a hyperbola, or a parabola.
//  It is represented by a quadratic equation in two nonhomogeneous
//  or three homogeneous coordinates.  Conversely, every quadratic
//  equation represents a conic, be it that it can be degenerate:
//  either in two (intersecting or parallel) lines, or in two
//  coincident lines.  Also, it can have no "visible", real points,
//  when it is an imaginary ellipse, or consist of two complementary
//  imaginary lines in which case it only has one real point, which could
//  still be at infinity.
//
//  These 11 cases are the possible values of CConic2D::real_type().
//  The default constructor sets the type to "invalid conic";
//  otherwise the correct type is automatically set when the equation
//  of the conic is given to the constructor that takes 6 numeric values
//  (a,b,c,d,e,f): the cartesian equation is then
//  \f$ Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0 \f$
//
//  When the conic is degenerate and consists of two lines, the method
//  components() returns a list of two (possibly identical) lines.
//  Otherwise, this method returns an empty list.
// \see http://en.wikipedia.org/wiki/Conic_section#Discriminant_classification
// \see http://en.wikipedia.org/wiki/Matrix_representation_of_conic_sections
 *
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CConic2D: public CConic
{

private:
  // DATA MEMBERS

  float a_; //!< coefficient of \a x^2
  float b_; //!< coefficient of \a xy
  float c_; //!< coefficient of \a y^2
  float d_; //!< coefficient of \a x
  float e_; //!< coefficient of \a y
  float f_; //!< free term

  private:
  /** 
  \brief  Returns the coefficient of \f$X^2\f$ 
  \note A = a()  .'. a() = A
  */
  SMF_INLINE float a() const { return  a_; }

  /** 
  \brief Returns the coefficient of \f$XY\f$
  \note B = b()/2  .'. b() = 2*B
  **/
  SMF_INLINE float b() const { return  b_; }

  /** 
  \brief   Returns the coefficient of \f$Y^2\f$
  \note C= c()   .". C = c()
  */
  SMF_INLINE float c() const { return  c_; }

  /** 
  \brief Returns the coefficient of \f$X\f$
  \note 2*D= d()   .". D = d() / 2
  **/
  SMF_INLINE float d() const { return  d_; }

  /**
  \brief Returns the coefficient of \f$Y\f$
  \note E= e() / 2 .". e() = 2*E
  **/
  SMF_INLINE float e() const { return  e_; }

  /**
  \brief  Returns the free term
  */
  SMF_INLINE float f() const { return  f_; }

public:
  /** 
  \brief  Returns the parameter A from general formula 
  */
  SMF_INLINE float A() const { return  a_; }

  /** 
  \brief Returns the parameter B from general formula
  \note B = b()/2  .'. b() = 2*B
  **/
  SMF_INLINE float B() const { return  b_* 0.5f; }

  /** 
  \brief  Returns the parameter C from general formula
  */
  SMF_INLINE float C() const { return  c_; }

  /** Returns the parameter D from general formula
  \note 2*D= d()   .". D = d() / 2
  **/
  SMF_INLINE float D() const { return  d_*0.5f; }

  /** Returns the parameter E from general formula 
  \note E= e() / 2 .". e() = 2*E
  **/
  SMF_INLINE float E() const { return  e_*0.5f; }

  /** 
  \brief  Returns the parameter F from general formula
  */
  SMF_INLINE float F() const { return  f_; }

  // CONSTRUCTORS AND RELATED STUFF

  /// default constructor
  CConic2D(){}

  /// copy constructor
  CConic2D(CConic2D const& c)
    : CConic(c), a_(c.a()), b_(c.b()), c_(c.c()), d_(c.d()), e_(c.e()), f_(c.f()) {}
  /// assignment operator
  CConic2D& operator=(CConic2D const& c) {
	//base class assignment
	  CConic::operator=(c);
	a_=c.a(); b_=c.b(); c_=c.c(); d_=c.d(); e_=c.e(); f_=c.f();
    return *this;
  }
  /// destructor
  ~CConic2D() {}


  /** 
  \brief  constructor using polynomial coefficients.
    The order of the coefficients is: $X^2$, $XY$, $Y^2$, $X$, $Y$, free,
  */
  CConic2D(float a, float b, float c, float d, float e, float f);

  /** 
  \brief  constructor using polynomial coefficients, given as a C array.
    The order of the coefficients is: $X^2$, $XY$, $Y^2$, $X$, $Y$, free,
  */
  CConic2D(float const coeff[]);

  /** 
  \brief  constructor using CCircle2D .
  */
  CConic2D(CCircle2D const &circ);

  /** 
  \brief clacule the Discriminator for this conic
  Expression B2-4A*C is called the discriminant of the general second degree polynomial. 
  Based on the above, if the value of the discriminant is less than, equal to or greater than zero, 
  the conic is an ellipse, a parabola, or a hyperbola. 
  */
  virtual float getDiscriminator() const;
  /** 
  \brief  Returns true if this conic is a degenerate case, i.e., if it consists of 2 lines.
  */
  virtual bool isDegenerate( float epsilon = 1e-6f) const;
  /** 
  \brief  Returns true if a central conic, i.e., an ellipse, circle, or hyperbola.
    Also the degenerate versions of these return true.
    Returns false if a parabola or two parallel or coinciding lines.
    (This is an affine property, not a projective one.)
	Equivalent to saying that the line at infinity does not touch the conic.
	*/
	virtual bool isCentral() const;
	/**
	\brief get point given the x axis value
	\return number of points returned
	**/
	virtual int getPoint(float x, CPoint2D &p1, CPoint2D &p2)const;
	
  /** 
  \brief  set or reset the conic using polynomial coefficients.
    The order of the coefficients is: $X^2$, $XY$, $Y^2$, $X$, $Y$, free,
	*/
  void setFromPol(float a, float b, float c, float d, float e, float f);

  /** set or reset the conic using General conic coefficients.
  //  The order of the coefficients is: A x2 + 2B xy + C y2 + 2D x + 2E y + F = 0
  // \param A - the A coeficient from the general form 
  // \param B - the B coeficient from the general form 
  // \param C - the C coeficient from the general form 
  // \param D - the D coeficient from the general form 
  // \param E - the E coeficient from the general form 
  // \param F - the F coeficient from the general form 
   **/
  void setFromGen(float A, float B, float C, float D, float E, float F);


  /** 
  \brief  comparison operator.
   Comparison is on the conic, not the equation coefficients.  Hence two
   conics are identical if their coefficient vectors are multiples of
   each other.
  */
  bool operator==(CConic2D const& c) const;

  /**
  \brief return the Characteristic Matriz of this conic:
  \note: 
  \f$ S = \begin{bmatrix}A & B & D\\B & C & E\\D & E & F\end{bmatrix}  \f$
  **/
  CMat3D characMat()const;
  /**
  \brief return the  Matriz of this conic:
  \note: 
  \f$ S = \begin{bmatrix}A & B/2 & D/2\\B/2 & C & E/2\\D/2 & E/2 & F\end{bmatrix}  \f$
  **/
  CMat3D eccentMat()const;
  CMat3D toMat3()const{ return characMat();};
  // UTILITY FUNCTIONS

  // Functions related to dual space ( http://en.wikipedia.org/wiki/Dual_space )---------------------------------

  /**
  \brief Returns the dual or tangential representation of this conic.
    The homogeneous coordinates of the points belonging to the dual conic
    are the coefficients of the equations of all tangents to the original
    conic.
	\note Function related to dual space ( http://en.wikipedia.org/wiki/Dual_space )---------------------------------
   */
    CConic2D dual_conic() const;

  /**
  \brief  Returns the dual or tangential representation of this conic.
  \note Function related to dual space ( http://en.wikipedia.org/wiki/Dual_space )---------------------------------
   */
	   CConic2D tangential_form() const { return dual_conic(); }

  /**
  \brief  Modify this conic by translating it over distance \a x in the \a X direction and distance \a y in the \a Y direction.
  */
  void translate_by(float x, float y);
  /**
  \brief  Returns the curvature of the conic at point p, assuming p is on the conic.
  */
  float curvature_at(CPoint2D const& p) const;
  /**
	\brief return ellipse excentricity
	\note For an ellipse, eccentricity is: c/a onde c= distância do centro ao foco e a= half-major
	\note gives the eccentricity e if the conic section is not a parabola (which has eccentricity equal to 1), not a degenerate hyperbola or degenerate ellipse, and not an imaginary ellipse
	**/
  virtual	float getEccentricity()const;
  /**
  \brief  Converts the coefficients to a geometric description of an ellipse. 
    Returns false if the conic is not an ellipse. Double is appropriate
    since integer coefficients can produce non-integer ellipse parameters.
   The centre of the ellipse is (xc, yc)
  */
  bool ellipse_geometry(float& xc, float& yc, float& major_axis_length,
                        float& minor_axis_length, float& angle_in_radians);
  static void testEllipseIntersection();
  /**
  \brief  Returns true if the point pt belongs (is inside) to the conic.
  */
  virtual bool contains(CPoint2D const& pt) const;
  /**
  \brief  Returns true if the point pt belongs  to the edge conic.
    I.e., if it *exactly* satisfies the conic equation.
	*/
  virtual bool edgeContains(CPoint2D const& pt, float epsilon = 1e-6f) const;

  /**
  \brief check points of intersection of two conic sections
  \see http://csharphelper.com/blog/2014/11/draw-a-conic-section-from-its-polynomial-equation-in-c/
  \see http://csharphelper.com/blog/2014/11/see-where-two-ellipses-intersect-in-c-part-1/
  \see http://csharphelper.com/blog/2014/11/see-where-two-ellipses-intersect-in-c-part-2/
  \see http://csharphelper.com/blog/2014/11/see-where-two-ellipses-intersect-in-c-part-3/
  \see http://csharphelper.com/blog/2014/11/see-where-two-ellipses-intersect-in-c-part-4/
  \see http://csharphelper.com/blog/2014/11/use-newtons-method-to-find-the-roots-of-equations-in-c/
  \return true if there are any intersection, 
  \param conic1
  \param conic2
  \param [out] numIntersPoints - Numberof intersection points [0 to 4]
  \param [out] p1
  \param [out] p2
  \param [out] p3
  \param [out] p4
  \param prec set the floating point precision on points calculation. 0 means to not use the precision (full precision).
  \warning points not used return NAN_FLOATS for x,y
  \warning NOT WORKING YET
  **/
  static bool intersection(const CConic2D &conic1, const CConic2D &conic2, float numIntersPoints,  CPoint2D &p1,  CPoint2D &p2, CPoint2D &p3, CPoint2D &p4,int prec=0);

 private:
  //--------------------------------------------------------------------------
  /**
  \brief  set conic type from polynomial coefficients and store in member type_
  // This method must be called by all constructors (except the default
  // constructor) and all methods that change the coefficients.
  */
  void set_type_from_equation();


/**
  \brief  Write "<CConic2D aX^2+bXY+cY^2+dX+eY+f>" to stream
// \see  CConic2D
*/

ostream&  operator<<(ostream& s);

/**
  \brief  Read a b c d e f from stream
// \see  CConic2D
*/

istream&  operator>>(istream& s);

};




} //end GEO
}  //end SMF

#endif // __SMF_CIRCLE_2D_
