/*
  SMF -  Super Math Fabric  
  Copyright (C) 2014 Rasputtim <Rasputtim@hotmail.com>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the freeMem Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the freeMem Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

  */

#ifndef _SMF__MATH_ODE_H__
#define _SMF__MATH_ODE_H__
#include "../SMF_Config.h"
namespace SMF{
	namespace MATH{





typedef void (*deriveFunction_t)( const float t, const void *userData, const float *state, float *derivatives );
typedef float (*numdifderiveFunction_t) (float x, void * params);

#define SMF_FN_EVAL(F,x) (*((F)->function))(x,(F)->params)

/**
 * \class CNumDEqu 
 *
 * \ingroup SMF_Math
 *
 * \if pt_br
 * \brief	Solucionador de derivadas numéricos por diferenciação.
 * 
 * \elseif us_en
 * \brief 	compute numerical derivatives by finite differencing.
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
class SMF_API SMF_API CNumDEqu {

public:
	virtual				~CNumDEqu(){}

	/**
	\brief default construtor
	\warning remember to set up function and params before use the methods of the class
	**/
	CNumDEqu():params(NULL) {};

	/**
	\brief  construtor
	\param func function to devivate
	**/
	CNumDEqu(numdifderiveFunction_t func):params(NULL){function=func;};

	/**
	\brief  construtor
	\param func function to devivate
	\param par user data params
	**/
	CNumDEqu(numdifderiveFunction_t func, void * par){function=func; params=par;};
	

	/// assignment operator
	CNumDEqu& operator=(CNumDEqu const& c) {
		function=c.function;
		params=c.params;
		return *this;
	  }

	/**
    This function computes the numerical derivative of the function at the point x using an adaptive central difference algorithm with a step-size of h. 
	The derivative is returned in result and an estimate of its absolute error is returned in abserr.
    The initial value of h is used to estimate an optimal step-size, based on the scaling of the truncation error and round-off error in the derivative calculation. 
	The derivative is computed using a 5-point rule for equally spaced abscissae at x-h, x-h/2, x, x+h/2, x+h, with an error estimate taken from the difference between the 5-point rule and the corresponding 3-point rule x-h, x, x+h. Note that the value of the function at x does not contribute to the derivative calculation, so only 4-points are actually used. 
	*/
	int deriv_central (float x, float h,   float *result, float *abserr);
	/**
    This function computes the numerical derivative of the function at the point x using an adaptive backward difference algorithm with a step-size of h. 
	The function is evaluated only at points less than x, and never at x itself. The derivative is returned in result and an estimate of its absolute error is returned in abserr. 
	This function should be used if f(x) has a discontinuity at x, or is undefined for values greater than x.
    This function is equivalent to calling deriv_forward with a negative step-size. 
	*/
	int deriv_backward (float x, float h, float *result, float *abserr);
	/**
    This function computes the numerical derivative of the function at the point x using an adaptive forward difference algorithm with a step-size of h. 
	The function is evaluated only at points greater than x, and never at x itself. 
	The derivative is returned in result and an estimate of its absolute error is returned in abserr. 
	This function should be used if f(x) has a discontinuity at x, or is undefined for values less than x.
    The initial value of h is used to estimate an optimal step-size, based on the scaling of the truncation error and round-off error in the derivative calculation. The derivative at x is computed using an “open” 4-point rule for equally spaced abscissae at x+h/4, x+h/2, x+3h/4, x+h, with an error estimate taken from the difference between the 4-point rule and the corresponding 2-point rule x+h/2, x+h. 
	**/
	int deriv_forward ( float x, float h, float *result, float *abserr);


/**
The following code estimates the derivative of the function f(x) = x^{3/2} at x=2 and at x=0. 
The function f(x) is undefined for x<0 so the derivative at x=0 is computed using deriv_forward.
Here is the output of the program,

$ ./a.out

f(x) = x^(3/2)
x = 2.0
f'(x) = 2.1213203120 +/- 0.0000004064
exact = 2.1213203436

x = 0.0
f'(x) = 0.0000000160 +/- 0.0000000339
exact = 0.0000000000
\see
The algorithms used by these functions are described in the following sources:

    Abramowitz and Stegun, Handbook of Mathematical Functions, Section 25.3.4, and Table 25.5 (Coefficients for Differentiation).
    S.D. Conte and Carl de Boor, Elementary Numerical Analysis: An Algorithmic Approach, McGraw-Hill, 1972. 
**/
static int testDerivative(void);

public:

	numdifderiveFunction_t function;		// derive function
    void * params;                          // client data

};


/**
 * \class CODEqu 
 *
 * \ingroup SMF_Math
 *
 * \if pt_br
 * \brief	Solucionador numéricos para equações diferenciais ordinárias.
 * 
 * \elseif us_en
 * \brief 	Numerical solvers for ordinary differential equations (ODE).
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
class SMF_API CODEqu {

public:
	virtual				~CODEqu() {}

	virtual float		evaluate( const float *state, float *newState, float t0, float t1 ) = 0;

protected:
	int					dimension;		// dimension in floats allocated for
	deriveFunction_t	derive;			// derive function
	const void *		userData;		// client data
};


/**
 * \class CODEqu_Euler 
 *
 * \ingroup SMF_Math
 *
 * \if pt_br
 * \brief	Solucionador numéricos para equações diferenciais de Euler.
 * 
 * \elseif us_en
 * \brief 	Numerical solvers for ordinary Euler differential equations.
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
class SMF_API CODEqu_Euler : public CODEqu {

public:
						CODEqu_Euler( const int dim, const deriveFunction_t dr, const void *ud );
	virtual				~CODEqu_Euler();

	virtual float		evaluate( const float *state, float *newState, float t0, float t1 );

protected:
	float *				derivatives;	// space to store derivatives
};


/**
 * \class CODEqu_Midpoint 
 *
 * \ingroup SMF_Math
 *
 * \if pt_br
 * \brief	Solucionador numéricos para equações diferenciais de Ponto médio.
 * 
 * \elseif us_en
 * \brief 	Numerical solvers for ordinary MidPoint differential equations.
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

class SMF_API CODEqu_Midpoint : public CODEqu {

public:
						CODEqu_Midpoint( const int dim, const deriveFunction_t dr, const void *ud );
	virtual				~CODEqu_Midpoint();

	virtual float		evaluate( const float *state, float *newState, float t0, float t1 );

protected:
	float *				tmpState;
	float *				derivatives;	// space to store derivatives
};

/**
 * \class CODEqu_RK4
 *
 * \ingroup SMF_Math
 *
 * \if pt_br
 * \brief	Solucionador numéricos para equações diferenciais RK4.
 * 
 * \elseif us_en
 * \brief 	Numerical solvers for ordinary RK4 differential equations.
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

class SMF_API CODEqu_RK4 : public CODEqu {

public:
						CODEqu_RK4( const int dim, const deriveFunction_t dr, const void *ud );
	virtual				~CODEqu_RK4();

	virtual float		evaluate( const float *state, float *newState, float t0, float t1 );

protected:
	float *				tmpState;
	float *				d1;				// derivatives
	float *				d2;
	float *				d3;
	float *				d4;
};

/**
 * \class 	CODEqu_RK4Adaptive
 
 *
 * \ingroup SMF_Math
 *
 * \if pt_br
 * \brief	Solucionador numéricos para equações diferenciais RK4 adaptativas.
 * 
 * \elseif us_en
 * \brief 	Numerical solvers for ordinary RK4 differential equations.
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
class SMF_API CODEqu_RK4Adaptive : public CODEqu {

public:
						CODEqu_RK4Adaptive( const int dim, const deriveFunction_t dr, const void *ud );
	virtual				~CODEqu_RK4Adaptive();

	virtual float		evaluate( const float *state, float *newState, float t0, float t1 );
	void				setMaxError( const float err );

protected:
	float				maxError;		// maximum allowed error
	float *				tmpState;
	float *				d1;				// derivatives
	float *				d1half;
	float *				d2;
	float *				d3;
	float *				d4;
};
} //end MATH
} //end SMF
#endif /* !__MATH_ODE_H__ */
