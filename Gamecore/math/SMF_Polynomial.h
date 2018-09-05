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

#ifndef _SMF__MATH_POLYNOMIAL_H__
#define _SMF__MATH_POLYNOMIAL_H__
#include "../SMF_Config.h"
#include "SMF_Complex.h"
#include "SMF_Math.h"
#include "../util/SMF_Heap.h"
#include "../math/SMF_Vector.h"
namespace SMF {

namespace MATH {


/**
 * \class CPolyConSystem 
 *
 * \ingroup SMF_Math
 *
 * \if pt_br
 * \brief Resolve o system de duas equações polinomiais de conicas:
	\n Equação 1 ->  \f$ Ax^2, Bxy, Cy^2, Dx, Ey, F = 0 \f$ 
	\n Equação 2 ->  \f$ Gx^2, Hxy, Iy^2, Jx, Ky, L = 0 \f$ 
 * \elseif us_en
 * \brief solve a system of two equations of conic types:
	\n Equation 1 ->  \f$ Ax^2, Bxy, Cy^2, Dx, Ey, F = 0 \f$ 
	\n Equation 2 ->  \f$ Gx^2, Hxy, Iy^2, Jx, Ky, L = 0 \f$ 
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
class SMF_API CPolyConSystem{
public:
	/**
	\brief constructor
	\note this constructor does not set the system values.
	\warning remember to set up the system values
	**/
	CPolyConSystem();
	/**
	\brief constructor
	\param a,b,c,d,e,f coeficients from Ax^2, Bxy, Cy^2, Dx, Ey, F = 0
	\param g,h,i,j,k,l coeficients from Gx^2, Hxy, Iy^2, Jx, Ky, L = 0
	**/
	CPolyConSystem(float a,float b, float c,float d, float e, float f,
		float g,float h, float i,float j, float k, float l);
	void setEqu1(float a,float b, float c,float d, float e, float f);
	void setEqu2(float g,float h, float i,float j, float k, float l);
	SMF_INLINE int getNumRealSol()const{return numRealSol;}
	/**
	\brief get real solutions for the system
	\param [out] sol1,sol2,sol3,sol4 real solutions
	\return number of real solutions
	**/
	int getSol(CVec2D &sol1, CVec2D &sol2, CVec2D &sol3, CVec2D &sol4,int prec=0);
	
	/**
	\brief return the first pair (x,y) solution for the system 
	\warning if there is no solution it will return a CVec2D::nan;
	**/
	SMF_INLINE CVec2D getSol1()const{ return sol1;}
	/**
	\brief return the first pair (x,y) solution for the system 
	\warning if there is no solution it will return a CVec2D::nan;
	**/
	SMF_INLINE CVec2D getSol2()const{ return sol2;}
	/**
	\brief return the first pair (x,y) solution for the system 
	\warning if there is no solution it will return a CVec2D::nan;
	**/
	SMF_INLINE CVec2D getSol3()const { return sol3;};
	/**
	\brief return the first pair (x,y) solution for the system 
	\warning if there is no solution it will return a CVec2D::nan;
	**/
	SMF_INLINE CVec2D getSol4()const { return sol4;};
	/**
	\param precision precisionto round the result. 0 means no use of precision
	**/
	void solveSystem(int precision=0);
	CVec2D sol1,sol2,sol3,sol4;
private:
	float A,B,C,D,E,F,G,H,I,J,K,L;
	float alpha, beta, gamma,delta, epsilon;
	float Y1,Y2,Y3,Y4;
	float X1A,X2A,X3A,X4A;
	float X1B,X2B,X3B,X4B;
	///number of realsolutions of the system
	int numRealSol;
	//indica que o systema já foi resolvido
	bool solved;
	/**
	\brief if one solution is found, put it on the right solution variable
	**/
	void setSol(float x,float y);
	/**
	\encontra os valores de y
	\return numero de raizes reais, ou soluções reais
	**/
	int solveY();
	void calcCoeff();
	/**
	\brief calcula o valor de x para o dado y quando x1==x2
	**/
	float calcX_a0(float y);
	/**
	\brief calcula o valor de x1 para o dado y
	**/
	float calcX1(float y);
	/**
	\brief calcula o valor de x2 para o dado y
	**/
	float calcX2(float y);
	/*
	\brief reseta os valores dos coeficientes e resultados
	*/
	void resetCoeff();
	/*
	\brief testa o valor nas duas equações para verificar se o sistema foi solucionado
	*/
	bool testValues( float x1, float y);
	
};



/**
 * \class CPolynomial
 *
 * \ingroup SMF_Math
 *
 * \if pt_br
 * \brief polinônio de grau arbitrário com coeficientes reais
   \note Um dos tipos mais simples de funções que se constrói mediante a aplicação repetida das operações elementares, adição e multiplicação ,
são as funções racionais ou polinômios.

    a_k x^k + ... + a_1 x^1 + a_0, 

Aplicando-se estas operações a uma variável independente x e a um conjunto de números reais 
ou complexos obtem-se os polinômios:

Onde n é um número natural e x também chamado de variável independente , pode assumir valores reais ou complexos .

Portanto:

onde ak, ... a1, a0 são coeficientes
para generalizar essa expressão em todos os casos, é necessário admitir que 0 seja um coeficiente possível.

Para o maior i com ai ≠ 0 (se houver), ai é chamado de coeficiente líder do polinômio.

Por exemplo, o coeficiente líder do polinômio seguinte:

     4x^5 + x^3 + 2x^2

é 4.

 * \elseif us_en
 * \brief Polynomial of arbitrary degree with real coefficients.
 * \endif
 *
 * \see http://www.mathsisfun.com/algebra/polynomials.html ,http://www.cursinhovirtual.com.br/Polinom/Mat03.htm
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CPolynomial {
public:
					CPolynomial();
					/**
					\brief cria uma função polinomial com grau d
					\param d grau do polinomio
					**/
					explicit CPolynomial( int d );
					/**
					\brief cria uma função polinomial com dois coeficientes
					\param a coeficiente do grau 1 (x^1)
					\param b coeficiente do grau 0 (x^0)
					**/
					explicit CPolynomial( float a, float b );
					/**
					\brief cria uma função polinomial com três coeficientes
					\param a coeficiente do grau 2 (x^2)
					\param b coeficiente do grau 1 (x^1)
					\param c coeficiente do grau 0 (x^0)
					**/
					explicit CPolynomial( float a, float b, float c );
					/**
					\brief cria uma função polinomial com quatro coeficientes
					\param a coeficiente do grau 3 (x^3)
					\param b coeficiente do grau 2 (x^2)
					\param c coeficiente do grau 1 (x^1)
					\param d coeficiente do grau 0 (x^0)
					**/
					explicit CPolynomial( float a, float b, float c, float d );
					/**
					\brief cria uma função polinomial com cinco coeficientes
					\param a coeficiente do grau 4 (x^4)
					\param b coeficiente do grau 3 (x^3)
					\param c coeficiente do grau 2 (x^2)
					\param d coeficiente do grau 1 (x^1)
					\param e coeficiente do grau 0 (x^0)
					**/
					explicit CPolynomial( float a, float b, float c, float d, float e );

	float			operator[]( int index ) const;
	float &			operator[]( int index );

	CPolynomial	operator-() const;
	CPolynomial &	operator=( const CPolynomial &p );

	CPolynomial	operator+( const CPolynomial &p ) const;
	CPolynomial	operator-( const CPolynomial &p ) const;
	CPolynomial	operator*( const float s ) const;
	CPolynomial	operator/( const float s ) const;

	CPolynomial &	operator+=( const CPolynomial &p );
	CPolynomial &	operator-=( const CPolynomial &p );
	CPolynomial &	operator*=( const float s );
	CPolynomial &	operator/=( const float s );

	bool			compare( const CPolynomial &p ) const;						// exact compare, no epsilon
	bool			compare( const CPolynomial &p, const float epsilon ) const;// compare with epsilon
	bool			operator==(	const CPolynomial &p ) const;					// exact compare, no epsilon
	bool			operator!=(	const CPolynomial &p ) const;					// exact compare, no epsilon

	void			toZero();
	void			zero( int d );
	/**
	\brief retorna o grau do polinômio
	\return o grau do polinômio
	**/
	int				getDimension() const;									// get the degree of the polynomial
	/**
	\brief retorna o grau do polinômio
	\return o grau do polinômio
	**/
	int				getDegree() const;									// get the degree of the polynomial
	/**
	\brief calcula o valor do polinônio, dada o valor da variável independente (x)
	\param x valor da variável independente (Real)
	\return o valor calculado do polinômio (real)
	**/
	float			getValue( const float x ) const;							// evaluate the polynomial with the given real value
	/**
	\brief calcula o valor do polinônio, dada o valor da variável independente (x)
	\param x valor da variável independente (Complexo)
	\return o valor calculado do polinômio (Complexo)
	**/
	CComplex		getValue( const CComplex &x ) const;						// evaluate the polynomial with the given complex value
	/**
	\brief calcula a primeira derivada do polinômio
	\return retorna um polinômio que corresponde à primeira derivada do polinômio original
	**/
	CPolynomial	getDerivative() const;								// get the first derivative of the polynomial
	/**
	\brief calcula a integral, ou anti-derivada do polinômio
	\return retorna um polinômio que corresponde à integral, ou anti-derivada do polinômio original
	**/
	CPolynomial	getAntiDerivative() const;							// get the anti derivative of the polynomial

	/**
	\brief calcula todas as raizes complexas do polinômio
	\param roots ponteiro para o vetor de complexos que será preenchido com as raizes complexas
	\return retorna o numero de raizes
	**/
	int				getRoots( CComplex *roots ) const;							// get all roots
	/**
	\brief calcula todas as raizes reais do polinômio
	\param roots ponteiro para o vetor de floats que será preenchido com as raizes reais
	\return retorna o numero de raizes
	**/
	int				getRoots( float *roots ) const;								// get the real roots
	

	/**
	\brief calcula todas as raizes reais do polinômio de grau 1
	\param roots ponteiro para o vetor de floats que será preenchido com as raizes reais
	\return retorna o numero de raizes
	**/
	static int		getRoots1( float a, float b, float *roots );
	/**
	\brief calcula todas as raizes reais do polinômio de grau 2
	\param roots ponteiro para o vetor de floats que será preenchido com as raizes reais
	\return retorna o numero de raizes
	**/
	static int		getRoots2( float a, float b, float c, float *roots );
	/**
	\brief calcula todas as raizes reais do polinômio de grau 3
	\param roots ponteiro para o vetor de floats que será preenchido com as raizes reais
	\return retorna o numero de raizes
	**/
	static int		getRoots3( float a, float b, float c, float d, float *roots );
	/**
	\brief calcula todas as raizes reais do polinômio de grau 4
	\param roots ponteiro para o vetor de floats que será preenchido com as raizes reais
	\return retorna o numero de raizes
	**/
	static int		getRoots4( float a, float b, float c, float d, float e, float *roots );

	const float *	toFloatPtr() const;
	float *			toFloatPtr();
	const char *	toString( int precision = 2 ) const;

	static void		Test();

		/// solves a quadratic equation ax^2 + bx + c = 0 for x. Returns the number of roots found.
	static int solveQuadratic(float a, float b, float c, float &root1, float &root2);

	/** solves a cubic quadratic  c2x^2 + c1x + c0 = 0 for x. Returns the number of roots found.
	 \params     c[0] + c[1]*x + c[2]*x^2 + c[3]*x^3 + c[4]*x^4 = 0
	 **/
	static int solveQuadratic(float c[ 3 ], float s[ 2 ]);
	/** solves a cubic equation c3x^3 + c2x^2 + c1x + c0 = 0 for x. Returns the number of roots found.
	 **/
	 static int solveCubic(float c3, float c2, float c1, float c0, float &root1, float &root2, float &root3);
	/** solves a cubic equation c3x^3 + c2x^2 + c1x + c0 = 0 for x. Returns the number of roots found.
	 \params     c[0] + c[1]*x + c[2]*x^2 + c[3]*x^3 + c[4]*x^4 = 0
	 **/
	 static int solveCubic(float c[ 4 ], float s[ 3 ]);
	
	/// solves a quartic equation c4x^4 +c3x^3 + c2x^2 + c1x + c0 = 0 for x. Returns the number of roots found.
	static int solveQuartic(float c4,float c3, float c2, float c1, float c0, float &root1, float &root2, float &root3, float &root4);
	/**
	\brief method to solve quartic equation 
	\note created by my son Arthur
	**/
	static void solveQuarticArthur(const float& ALPHA, const float& BETA, const float& GAMMA, const float& DELTA, const float& EPSILON, CComplex &RR1, CComplex &RR2, CComplex &RR3, CComplex &RR4);
	/**
	\brief method to solve cubic equation 
	\note created by my son Arthur
	**/
	static void solveCubicArthur(const float& ALPHA, const float& BETA, const float& GAMMA, const float& DELTA, CComplex &RR1, CComplex &RR2, CComplex &RR3);
	/**
	\brief method to solve quadratic equation 
	\note created by my son Arthur
	**/
	static void solveQuadraticArthur(const float& ALPHA, const float& BETA, const float& GAMMA, CComplex &R1, CComplex &R2);
	/**
	\brief method to solve simple equation 
	\note created by my son Arthur
	**/
	static void solveSimple(const float& ALPHA, const float& BETA, CComplex &R1);
	/**
	\brief calcula todas as raizes reais do polinômio
	\param c0,c1,c2,c3,c4,c5 coeficientes do polinomio c1 x^4+ c2 x^3 + c3 x^2 + c4 x^1 + c5 x^0
	\param [out] roots ponteiro para o vetor de CComplex que será preenchido com as raizes reais
	\return retorna o numero de raizes
	\note se você deseja só as raizes reais, basta retirar a parte real do resultado complexo
	**/
	static int		solveQuarticLaguer(float c1,float c2,float c3, float c4, float c5, CComplex *roots );
	/**
	\brief numeric method to find roots
	\see http://en.wikipedia.org/wiki/Laguerre%27s_method
	**/
	static	int				Laguer( const CComplex *coef, const int degree, CComplex &r );

private:
	int				degree;
	int				allocated;
	float *			coefficient;

	void			resize( int d, bool keep );

};

SMF_INLINE_FORCED CPolynomial::CPolynomial() {
	degree = -1;
	allocated = 0;
	coefficient = NULL;
}

SMF_INLINE_FORCED CPolynomial::CPolynomial( int d ) {
	degree = -1;
	allocated = 0;
	coefficient = NULL;
	resize( d, false );
}

SMF_INLINE_FORCED CPolynomial::CPolynomial( float a, float b ) {
	degree = -1;
	allocated = 0;
	coefficient = NULL;
	resize( 1, false );
	coefficient[0] = b;
	coefficient[1] = a;
}

SMF_INLINE_FORCED CPolynomial::CPolynomial( float a, float b, float c ) {
	degree = -1;
	allocated = 0;
	coefficient = NULL;
	resize( 2, false );
	coefficient[0] = c;
	coefficient[1] = b;
	coefficient[2] = a;
}

SMF_INLINE_FORCED CPolynomial::CPolynomial( float a, float b, float c, float d ) {
	degree = -1;
	allocated = 0;
	coefficient = NULL;
	resize( 3, false );
	coefficient[0] = d;
	coefficient[1] = c;
	coefficient[2] = b;
	coefficient[3] = a;
}

SMF_INLINE_FORCED CPolynomial::CPolynomial( float a, float b, float c, float d, float e ) {
	degree = -1;
	allocated = 0;
	coefficient = NULL;
	resize( 4, false );
	coefficient[0] = e;
	coefficient[1] = d;
	coefficient[2] = c;
	coefficient[3] = b;
	coefficient[4] = a;
}

SMF_INLINE_FORCED float CPolynomial::operator[]( int index ) const {
	SMF_ASSERT( index >= 0 && index <= degree );
	return coefficient[ index ];
}

SMF_INLINE_FORCED float& CPolynomial::operator[]( int index ) {
	SMF_ASSERT( index >= 0 && index <= degree );
	return coefficient[ index ];
}

SMF_INLINE_FORCED CPolynomial CPolynomial::operator-() const {
	int i;
	CPolynomial n;

	n = *this;
	for ( i = 0; i <= degree; i++ ) {
		n[i] = -n[i];
	}
	return n;
}

SMF_INLINE_FORCED CPolynomial &CPolynomial::operator=( const CPolynomial &p ) { 
	resize( p.degree, false );
	for ( int i = 0; i <= degree; i++ ) {
		coefficient[i] = p.coefficient[i];
	}
	return *this;
}

SMF_INLINE_FORCED CPolynomial CPolynomial::operator+( const CPolynomial &p ) const {
	int i;
	CPolynomial n;

	if ( degree > p.degree ) {
		n.resize( degree, false );
		for ( i = 0; i <= p.degree; i++ ) {
			n.coefficient[i] = coefficient[i] + p.coefficient[i];
		}
		for ( ; i <= degree; i++ ) {
			n.coefficient[i] = coefficient[i];
		}
		n.degree = degree;
	} else if ( p.degree > degree ) {
		n.resize( p.degree, false );
		for ( i = 0; i <= degree; i++ ) {
			n.coefficient[i] = coefficient[i] + p.coefficient[i];
		}
		for ( ; i <= p.degree; i++ ) {
			n.coefficient[i] = p.coefficient[i];
		}
		n.degree = p.degree;
	} else {
		n.resize( degree, false );
		n.degree = 0;
		for ( i = 0; i <= degree; i++ ) {
			n.coefficient[i] = coefficient[i] + p.coefficient[i];
			if ( n.coefficient[i] != 0.0f ) {
				n.degree = i;
			}
		}
	}
	return n;
}

SMF_INLINE_FORCED CPolynomial CPolynomial::operator-( const CPolynomial &p ) const {
	int i;
	CPolynomial n;

	if ( degree > p.degree ) {
		n.resize( degree, false );
		for ( i = 0; i <= p.degree; i++ ) {
			n.coefficient[i] = coefficient[i] - p.coefficient[i];
		}
		for ( ; i <= degree; i++ ) {
			n.coefficient[i] = coefficient[i];
		}
		n.degree = degree;
	} else if ( p.degree >= degree ) {
		n.resize( p.degree, false );
		for ( i = 0; i <= degree; i++ ) {
			n.coefficient[i] = coefficient[i] - p.coefficient[i];
		}
		for ( ; i <= p.degree; i++ ) {
			n.coefficient[i] = - p.coefficient[i];
		}
		n.degree = p.degree;
	} else {
		n.resize( degree, false );
		n.degree = 0;
		for ( i = 0; i <= degree; i++ ) {
			n.coefficient[i] = coefficient[i] - p.coefficient[i];
			if ( n.coefficient[i] != 0.0f ) {
				n.degree = i;
			}
		}
	}
	return n;
}

SMF_INLINE_FORCED CPolynomial CPolynomial::operator*( const float s ) const {
	CPolynomial n;

	if ( s == 0.0f ) {
		n.degree = 0;
	} else {
		n.resize( degree, false );
		for ( int i = 0; i <= degree; i++ ) {
			n.coefficient[i] = coefficient[i] * s;
		}
	}
	return n;
}

SMF_INLINE_FORCED CPolynomial CPolynomial::operator/( const float s ) const {
	float invs;
	CPolynomial n;

	SMF_ASSERT( s != 0.0f );
	n.resize( degree, false );
	invs = 1.0f / s;
	for ( int i = 0; i <= degree; i++ ) {
		n.coefficient[i] = coefficient[i] * invs;
	}
	return n;
}

SMF_INLINE_FORCED CPolynomial &CPolynomial::operator+=( const CPolynomial &p ) {
	int i;

	if ( degree > p.degree ) {
		for ( i = 0; i <= p.degree; i++ ) {
			coefficient[i] += p.coefficient[i];
		}
	} else if ( p.degree > degree ) {
		resize( p.degree, true );
		for ( i = 0; i <= degree; i++ ) {
			coefficient[i] += p.coefficient[i];
		}
		for ( ; i <= p.degree; i++ ) {
			coefficient[i] = p.coefficient[i];
		}
	} else {
		for ( i = 0; i <= degree; i++ ) {
			coefficient[i] += p.coefficient[i];
			if ( coefficient[i] != 0.0f ) {
				degree = i;
			}
		}
	}
	return *this;
}

SMF_INLINE_FORCED CPolynomial &CPolynomial::operator-=( const CPolynomial &p ) {
	int i;

	if ( degree > p.degree ) {
		for ( i = 0; i <= p.degree; i++ ) {
			coefficient[i] -= p.coefficient[i];
		}
	} else if ( p.degree > degree ) {
		resize( p.degree, true );
		for ( i = 0; i <= degree; i++ ) {
			coefficient[i] -= p.coefficient[i];
		}
		for ( ; i <= p.degree; i++ ) {
			coefficient[i] = - p.coefficient[i];
		}
	} else {
		for ( i = 0; i <= degree; i++ ) {
			coefficient[i] -= p.coefficient[i];
			if ( coefficient[i] != 0.0f ) {
				degree = i;
			}
		}
	}
	return *this;
}

SMF_INLINE_FORCED CPolynomial &CPolynomial::operator*=( const float s ) {
	if ( s == 0.0f ) {
		degree = 0;
	} else {
		for ( int i = 0; i <= degree; i++ ) {
			coefficient[i] *= s;
		}
	}
	return *this;
}

SMF_INLINE_FORCED CPolynomial &CPolynomial::operator/=( const float s ) {
	float invs;

	SMF_ASSERT( s != 0.0f );
	invs = 1.0f / s;
	for ( int i = 0; i <= degree; i++ ) {
		coefficient[i] = invs;
	}
	return *this;;
}

SMF_INLINE_FORCED bool CPolynomial::compare( const CPolynomial &p ) const {
	if ( degree != p.degree ) {
		return false;
	}
	for ( int i = 0; i <= degree; i++ ) {
		if ( coefficient[i] != p.coefficient[i] ) {
			return false;
		}
	}
	return true;
}

SMF_INLINE_FORCED bool CPolynomial::compare( const CPolynomial &p, const float epsilon ) const {
	if ( degree != p.degree ) {
		return false;
	}
	for ( int i = 0; i <= degree; i++ ) {
		if ( CMath::fabs( coefficient[i] - p.coefficient[i] ) > epsilon ) {
			return false;
		}
	}
	return true;
}

SMF_INLINE_FORCED bool CPolynomial::operator==( const CPolynomial &p ) const {
	return compare( p );
}

SMF_INLINE_FORCED bool CPolynomial::operator!=( const CPolynomial &p ) const {
	return !compare( p );
}

SMF_INLINE_FORCED void CPolynomial::toZero() {
	degree = 0;
}

SMF_INLINE_FORCED void CPolynomial::zero( int d ) {
	resize( d, false );
	for ( int i = 0; i <= degree; i++ ) {
		coefficient[i] = 0.0f;
	}
}

SMF_INLINE_FORCED int CPolynomial::getDimension() const {
	return degree;
}

SMF_INLINE_FORCED int CPolynomial::getDegree() const {
	return degree;
}

SMF_INLINE_FORCED float CPolynomial::getValue( const float x ) const {
	float y, z;
	y = coefficient[0];
	z = x;
	for ( int i = 1; i <= degree; i++ ) {
		y += coefficient[i] * z;
		z *= x;
	}
	return y;
}

SMF_INLINE_FORCED CComplex CPolynomial::getValue( const CComplex &x ) const {
	CComplex y, z;
	y.set( coefficient[0], 0.0f );
	z = x;
	for ( int i = 1; i <= degree; i++ ) {
		y += coefficient[i] * z;
		z *= x;
	}
	return y;
}

SMF_INLINE_FORCED CPolynomial CPolynomial::getDerivative() const {
	CPolynomial n;

	if ( degree == 0 ) {
		return n;
	}
	n.resize( degree - 1, false );
	for ( int i = 1; i <= degree; i++ ) {
		n.coefficient[i-1] = i * coefficient[i];
	}
	return n;
}

SMF_INLINE_FORCED CPolynomial CPolynomial::getAntiDerivative() const {
	CPolynomial n;

	if ( degree == 0 ) {
		return n;
	}
	n.resize( degree + 1, false );
	n.coefficient[0] = 0.0f;
	for ( int i = 0; i <= degree; i++ ) {
		n.coefficient[i+1] = coefficient[i] / ( i + 1 );
	}
	return n;
}

SMF_INLINE_FORCED int CPolynomial::getRoots1( float a, float b, float *roots ) {
	SMF_ASSERT( a != 0.0f );
	roots[0] = - b / a;
	return 1;
}

SMF_INLINE_FORCED int CPolynomial::getRoots2( float a, float b, float c, float *roots ) {
	float inva, ds;  //ds=delta

	if ( a != 1.0f ) {
		SMF_ASSERT( a != 0.0f );
		inva = 1.0f / a;
		c *= inva;
		b *= inva;
	}
	ds = b * b - 4.0f * c;
	if ( ds < 0.0f ) {
		return 0;
	} else if ( ds > 0.0f ) {
		ds = CMath::sqrt( ds );
		roots[0] = 0.5f * ( -b - ds );
		roots[1] = 0.5f * ( -b + ds );
		return 2;
	} else {
		roots[0] = 0.5f * -b;
		return 1;
	}
}

SMF_INLINE_FORCED int CPolynomial::getRoots3( float a, float b, float c, float d, float *roots ) {
	float inva, f, g, halfg, ofs, ds, dist, angle, cs, ss, t;

	if ( a != 1.0f ) {
		SMF_ASSERT( a != 0.0f );
		inva = 1.0f / a;
		d *= inva;
		c *= inva;
		b *= inva;
	}

	f = ( 1.0f / 3.0f ) * ( 3.0f * c - b * b );
	g = ( 1.0f / 27.0f ) * ( 2.0f * b * b * b - 9.0f * c * b + 27.0f * d );
	halfg = 0.5f * g;
	ofs = ( 1.0f / 3.0f ) * b;
	ds = 0.25f * g * g + ( 1.0f / 27.0f ) * f * f * f;

	if ( ds < 0.0f ) {
		dist = CMath::sqrt( ( -1.0f / 3.0f ) * f );
		angle = ( 1.0f / 3.0f ) * CMath::atan( CMath::sqrt( -ds ), -halfg );
		cs = CMath::cos( angle );
		ss = CMath::sin( angle );
		roots[0] = 2.0f * dist * cs - ofs;
		roots[1] = -dist * ( cs + CMath::SQRT_THREE * ss ) - ofs;
		roots[2] = -dist * ( cs - CMath::SQRT_THREE * ss ) - ofs;
		return 3;
	} else if ( ds > 0.0f )  {
		ds = CMath::sqrt( ds );
		t = -halfg + ds;
		if ( t >= 0.0f ) {
			roots[0] = CMath::pow( t, ( 1.0f / 3.0f ) );
		} else {
			roots[0] = -CMath::pow( -t, ( 1.0f / 3.0f ) );
		}
		t = -halfg - ds;
		if ( t >= 0.0f ) {
			roots[0] += CMath::pow( t, ( 1.0f / 3.0f ) );
		} else {
			roots[0] -= CMath::pow( -t, ( 1.0f / 3.0f ) );
		}
		roots[0] -= ofs;
		return 1;
	} else {
		if ( halfg >= 0.0f ) {
			t = -CMath::pow( halfg, ( 1.0f / 3.0f ) );
		} else {
			t = CMath::pow( -halfg, ( 1.0f / 3.0f ) );
		}
		roots[0] = 2.0f * t - ofs;
		roots[1] = -t - ofs;
		roots[2] = roots[1];
		return 3;
	}
}

SMF_INLINE_FORCED int CPolynomial::getRoots4( float a, float b, float c, float d, float e, float *roots ) {
	int count;
	float inva, y, ds, r, s1, s2, t1, t2, tp, tm;
	float roots3[3];

	if ( a != 1.0f ) {
		SMF_ASSERT( a != 0.0f );
		inva = 1.0f / a;
		e *= inva;
		d *= inva;
		c *= inva;
		b *= inva;
	}

	count = 0;

	getRoots3( 1.0f, -c, b * d - 4.0f * e, -b * b * e + 4.0f * c * e - d * d, roots3 );
	y = roots3[0];
	ds = 0.25f * b * b - c + y;

	if ( ds < 0.0f ) {
		return 0;
	} else if ( ds > 0.0f ) {
		r = CMath::sqrt( ds );
		t1 = 0.75f * b * b - r * r - 2.0f * c;
		t2 = ( 4.0f * b * c - 8.0f * d - b * b * b ) / ( 4.0f * r );
		tp = t1 + t2;
		tm = t1 - t2;

		if ( tp >= 0.0f ) {
			s1 = CMath::sqrt( tp );
			roots[count++] = -0.25f * b + 0.5f * ( r + s1 );
			roots[count++] = -0.25f * b + 0.5f * ( r - s1 );
		}
		if ( tm >= 0.0f ) {
			s2 = CMath::sqrt( tm );
			roots[count++] = -0.25f * b + 0.5f * ( s2 - r );
			roots[count++] = -0.25f * b - 0.5f * ( s2 + r );
		}
		return count;
	} else {
		t2 = y * y - 4.0f * e;
		if ( t2 >= 0.0f ) {
			t2 = 2.0f * CMath::sqrt( t2 );
			t1 = 0.75f * b * b - 2.0f * c;
			if ( t1 + t2 >= 0.0f ) {
				s1 = CMath::sqrt( t1 + t2 );
				roots[count++] = -0.25f * b + 0.5f * s1;
				roots[count++] = -0.25f * b - 0.5f * s1;
			}
			if ( t1 - t2 >= 0.0f ) {
				s2 = CMath::sqrt( t1 - t2 );
				roots[count++] = -0.25f * b + 0.5f * s2;
				roots[count++] = -0.25f * b - 0.5f * s2;
			}
		}
		return count;
	}
}

SMF_INLINE_FORCED const float *CPolynomial::toFloatPtr() const {
	return coefficient;
}

SMF_INLINE_FORCED float *CPolynomial::toFloatPtr() {
	return coefficient;
}

SMF_INLINE_FORCED void CPolynomial::resize( int d, bool keep ) {
	int alloc = ( d + 1 + 3 ) & ~3;
	if ( alloc > allocated ) {
		float *ptr = (float *) mem_Alloc16( alloc * sizeof( float ) );
		if ( coefficient != NULL ) {
			if ( keep ) {
				for ( int i = 0; i <= degree; i++ ) {
					ptr[i] = coefficient[i];
				}
			}
			mem_Free16( coefficient );
		}
		allocated = alloc;
		coefficient = ptr;
	}
	degree = d;
}



} //end MATH
} //end SMF
#endif /* !__MATH_POLYNOMIAL_H__ */
