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

#ifndef _SMF__MATH_COMPLEX_H__
#define _SMF__MATH_COMPLEX_H__
#include "../SMF_Config.h"
#include "SMF_Math.h"
#include "exceptions/all.h"
namespace SMF {
namespace MATH{




/**
 * \class CComplex
 *
 * \ingroup SMF_Math
 *
 * \if pt_br
 * \brief Numeros Coomplexos
 * \elseif us_en
 * \brief Complex Number
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CComplex {
public:
    /// real part
	float				r;		
	/// imaginary part
	float				i;		

						CComplex();
						CComplex( const float r, const float i );
						// \brief imaginary part is 0.f
						CComplex( const float r);
	void 				set( const float r, const float i );
	void				toZero();
	float				real()const { return r;}
	float				operator[]( int index ) const;
	float &				operator[]( int index );

	CComplex			operator-() const;
	CComplex &			operator=( const CComplex &a );

	CComplex			operator*( const CComplex &a ) const;
	CComplex			operator/( const CComplex &a ) const;
	CComplex			operator+( const CComplex &a ) const;
	CComplex			operator-( const CComplex &a ) const;

	CComplex &			operator*=( const CComplex &a );
	CComplex &			operator/=( const CComplex &a );
	CComplex &			operator+=( const CComplex &a );
	CComplex &			operator-=( const CComplex &a );

	CComplex			operator*( const float a ) const;
	CComplex			operator/( const float a ) const;
	CComplex			operator+( const float a ) const;
	CComplex			operator-( const float a ) const;

	CComplex &			operator*=( const float a );
	CComplex &			operator/=( const float a );
	CComplex &			operator+=( const float a );
	CComplex &			operator-=( const float a );

	friend CComplex	operator*( const float a, const CComplex &b );
	friend CComplex	operator/( const float a, const CComplex &b );
	friend CComplex	operator+( const float a, const CComplex &b );
	friend CComplex	operator-( const float a, const CComplex &b );
    /// exact compare, no epsilon
	bool				compare( const CComplex &a ) const;						
	/// compare with epsilon
	bool				compare( const CComplex &a, const float epsilon ) const;	
	/// exact compare, no epsilon
	bool				operator==(	const CComplex &a ) const;						
	// exact compare, no epsilon
	bool				operator!=(	const CComplex &a ) const;						
	/* z=1/a */
	CComplex			Reciprocal() const;
	CComplex			inverse()const { return Reciprocal();}
	
	CComplex			sqrt() const;
	float				abs() const;
	/**
	\see http://www.milefoot.com/math/complex/functionsofi.htm
	**/
	CComplex			cos() const;
	/**
	\see http://www.milefoot.com/math/complex/functionsofi.htm
	**/
	CComplex			sin() const;
	/**
	\note \f$ \cos^{-1} i = -i\ln\left[i+i\sqrt{1-i^2}\right] = -i\ln\left[(1+\sqrt{2})i\right]  \f$
	**/
	CComplex			acos() const;
	/**
	\note \f$[ \begin{equation*} 
	\sin^{-1} x = -i\ln\left[ix\pm \sqrt{1-x^2}\right]
	\end{equation*}
	\f$]
	**/
	CComplex			asin() const;



	int					getDimension() const;

	const float *		toFloatPtr() const;
	float *				toFloatPtr();
	const char *		toString( int precision = 2 ) const;


	/* r=sqrt(z) */
	static CComplex sqrt (CComplex z);  
	/* r=sqrt(x) (x<0 ok) */
	static CComplex sqrt_real (float x);  
	/* r=a^b */
	static CComplex pow (CComplex a, CComplex b); 
	/* r=a^b */
	static CComplex pow_real (CComplex a, float b);  
	/* r=exp(a) */
	static CComplex exp (CComplex a); 

	/* r=log(a) (base e) */
	static CComplex log (CComplex a); 
	/* r=log10(a) (base 10) */
	static CComplex log10 (CComplex a);  
	 /* r=log_b(a) (base=b) */
	static CComplex log_b (CComplex a, CComplex b);  

	/* Complex Trigonometric Functions */
	/* r=sin(a) 
	\return the complex sine of the complex number z, \sin(z) = (\exp(iz) - \exp(-iz))/(2i). 
	*/
	static CComplex sin (CComplex a);
	/* r=cos(a) 
	\return the complex cosine of the complex number z, \cos(z) = (\exp(iz) + \exp(-iz))/2. 
	*/
	static CComplex cos (CComplex a);
	/* r=sec(a) */
	static CComplex sec (CComplex a);
	/* r=csc(a) */
	static CComplex csc (CComplex a);
	/* r=tan(a) */
	static CComplex tan (CComplex a);
	/* r=cot(a) */
	static CComplex cot (CComplex a);  

	/* Inverse Complex Trigonometric Functions */
	/* r=arcsin(a) */
	static CComplex arcsin (CComplex a);
	/* r=arcsin(a) */
	static CComplex arcsin_real (float a);
	/* r=arccos(a) */
	static CComplex arccos (CComplex a);
	/* r=arccos(a) */
	static CComplex arccos_real (float a); 
	/* r=arcsec(a) */
	static CComplex arcsec (CComplex a); 
	/* r=arcsec(a) */
	static CComplex arcsec_real (float a); 
	/* r=arccsc(a) */
	static CComplex arccsc (CComplex a);  
	/* r=arccsc(a) */
	static CComplex arccsc_real (float a); 
	/* r=arctan(a) */
	static CComplex arctan (CComplex a);
	/* r=arccot(a) */
	static CComplex arccot (CComplex a);  

	/* Complex Hyperbolic Functions */
	/* r=sinh(a) */
	static CComplex sinh (CComplex a);
	/* r=coshh(a) */
	static CComplex cosh (CComplex a);
	/* r=sech(a) */
	static CComplex sech (CComplex a); 
	/* r=csch(a) */
	static CComplex csch (CComplex a);  
	/* r=tanh(a) */
	static CComplex tanh (CComplex a);
	/* r=coth(a) */
	static CComplex coth (CComplex a);  

	/* Inverse Complex Hyperbolic Functions */
	/* r=arcsinh(a) */
	static CComplex arcsinh (CComplex a); 
	/* r=arccosh(a) */
	static CComplex arccosh (CComplex a);
	/* r=arccosh(a) */
	static CComplex arccosh_real (float a);
	/* r=arcsech(a) */
	static CComplex arcsech (CComplex a);
	/* r=arccsch(a) */
	static CComplex arccsch (CComplex a);
	/* r=arctanh(a) */
	static CComplex arctanh (CComplex a); 
	/* r=arctanh(a) */
	static CComplex arctanh_real (float a);
	/* r=arccoth(a) */
	static CComplex arccoth (CComplex a);  
	/* a Nan CComplex representation*/
	static const CComplex nan;

};

extern CComplex complex_origin;
#define complex_zero complex_origin

SMF_INLINE_FORCED CComplex::CComplex() {
}

SMF_INLINE_FORCED CComplex::CComplex( const float r) {
	this->r = r;
	this->i = 0.0f;
}
SMF_INLINE_FORCED CComplex::CComplex( const float r, const float i ) {
	this->r = r;
	this->i = i;
}

SMF_INLINE_FORCED void CComplex::set( const float r, const float i ) {
	this->r = r;
	this->i = i;
}

SMF_INLINE_FORCED void CComplex::toZero() {
	r = i = 0.0f;
}

SMF_INLINE_FORCED float CComplex::operator[]( int index ) const {
	SMF_ASSERT( index >= 0 && index < 2 );
	return ( &r )[ index ];
}

SMF_INLINE_FORCED float& CComplex::operator[]( int index ) {
	SMF_ASSERT( index >= 0 && index < 2 );
	return ( &r )[ index ];
}

SMF_INLINE_FORCED CComplex CComplex::operator-() const {
	return CComplex( -r, -i );
}

SMF_INLINE_FORCED CComplex &CComplex::operator=( const CComplex &a ) {
	r = a.r;
	i = a.i;
	return *this;
}

SMF_INLINE_FORCED CComplex CComplex::operator*( const CComplex &a ) const {
	return CComplex( r * a.r - i * a.i, i * a.r + r * a.i );
}

SMF_INLINE_FORCED CComplex CComplex::operator/( const CComplex &a ) const {
	float s, t;
	if ( CMath::fabs( a.r ) >= CMath::fabs( a.i ) ) {
		s = a.i / a.r;
		t = 1.0f / ( a.r + s * a.i );
		return CComplex( ( r + s * i ) * t, ( i - s * r ) * t );
	} else {
		s = a.r / a.i;
		t = 1.0f / ( s * a.r + a.i );
		return CComplex( ( r * s + i ) * t, ( i * s - r ) * t );
	}
}

SMF_INLINE_FORCED CComplex CComplex::operator+( const CComplex &a ) const {
	return CComplex( r + a.r, i + a.i );
}

SMF_INLINE_FORCED CComplex CComplex::operator-( const CComplex &a ) const {
	return CComplex( r - a.r, i - a.i );
}

SMF_INLINE_FORCED CComplex &CComplex::operator*=( const CComplex &a ) {
	*this = CComplex( r * a.r - i * a.i, i * a.r + r * a.i );
	return *this;
}

SMF_INLINE_FORCED CComplex &CComplex::operator/=( const CComplex &a ) {
	float s, t;
	if ( CMath::fabs( a.r ) >= CMath::fabs( a.i ) ) {
		s = a.i / a.r;
		t = 1.0f / ( a.r + s * a.i );
		*this = CComplex( ( r + s * i ) * t, ( i - s * r ) * t );
	} else {
		s = a.r / a.i;
		t = 1.0f / ( s * a.r + a.i );
		*this = CComplex( ( r * s + i ) * t, ( i * s - r ) * t );
	}
	return *this;
}

SMF_INLINE_FORCED CComplex &CComplex::operator+=( const CComplex &a ) {
	r += a.r;
	i += a.i;
	return *this;
}

SMF_INLINE_FORCED CComplex &CComplex::operator-=( const CComplex &a ) {
	r -= a.r;
	i -= a.i;
	return *this;
}

SMF_INLINE_FORCED CComplex CComplex::operator*( const float a ) const {
	return CComplex( r * a, i * a );
}

SMF_INLINE_FORCED CComplex CComplex::operator/( const float a ) const {
	float s = 1.0f / a;
	return CComplex( r * s, i * s );
}

SMF_INLINE_FORCED CComplex CComplex::operator+( const float a ) const {
	return CComplex( r + a, i );
}

SMF_INLINE_FORCED CComplex CComplex::operator-( const float a ) const {
	return CComplex( r - a, i );
}

SMF_INLINE_FORCED CComplex &CComplex::operator*=( const float a ) {
	r *= a;
	i *= a;
	return *this;
}

SMF_INLINE_FORCED CComplex &CComplex::operator/=( const float a ) {
	float s = 1.0f / a;
	r *= s;
	i *= s;
	return *this;
}

SMF_INLINE_FORCED CComplex &CComplex::operator+=( const float a ) {
	r += a;
	return *this;
}

SMF_INLINE_FORCED CComplex &CComplex::operator-=( const float a ) {
	r -= a;
	return *this;
}

SMF_INLINE_FORCED CComplex operator*( const float a, const CComplex &b ) {
	return CComplex( a * b.r, a * b.i );
}

SMF_INLINE_FORCED CComplex operator/( const float a, const CComplex &b ) {
	float s, t;
	if ( CMath::fabs( b.r ) >= CMath::fabs( b.i ) ) {
		s = b.i / b.r;
		t = a / ( b.r + s * b.i );
		return CComplex( t, - s * t );
	} else {
		s = b.r / b.i;
		t = a / ( s * b.r + b.i );
		return CComplex( s * t, - t );
	}
}

SMF_INLINE_FORCED CComplex operator+( const float a, const CComplex &b ) {
	return CComplex( a + b.r, b.i );
}

SMF_INLINE_FORCED CComplex operator-( const float a, const CComplex &b ) {
	return CComplex( a - b.r, -b.i );
}

SMF_INLINE_FORCED CComplex CComplex::Reciprocal() const {
	float s, t;
	if ( CMath::fabs( r ) >= CMath::fabs( i ) ) {
		s = i / r;
		t = 1.0f / ( r + s * i );
		return CComplex( t, - s * t );
	} else {
		s = r / i;
		t = 1.0f / ( s * r + i );
		return CComplex( s * t, - t );
	}
}

SMF_INLINE_FORCED CComplex CComplex::sqrt() const {
	float x, y, w;

	if ( r == 0.0f && i == 0.0f ) {
		return CComplex( 0.0f, 0.0f );
	}
	x = CMath::fabs( r );
	y = CMath::fabs( i );
	if ( x >= y ) {
		w = y / x;
		w = CMath::sqrt( x ) * CMath::sqrt( 0.5f * ( 1.0f + CMath::sqrt( 1.0f + w * w ) ) );
	} else {
		w = x / y;
		w = CMath::sqrt( y ) * CMath::sqrt( 0.5f * ( w + CMath::sqrt( 1.0f + w * w ) ) );
	}
	if ( w == 0.0f ) {
		return CComplex( 0.0f, 0.0f );
	}
	if ( r >= 0.0f ) {
		return CComplex( w, 0.5f * i / w );
	} else {
		return CComplex( 0.5f * y / w, ( i >= 0.0f ) ? w : -w );
	}
}



SMF_INLINE_FORCED CComplex CComplex::cos() const{
	//cos(a+bi)=cosa cosh b − i sina sinh b
	return this->cos (*this);
}
SMF_INLINE_FORCED CComplex CComplex::sin() const{
	//sin(a+bi)=sina coshb + i cosa sinhb
	return this->sin (*this);
}

SMF_INLINE_FORCED CComplex CComplex::acos() const{
	return this->arcsin(*this);
}

SMF_INLINE_FORCED CComplex CComplex::asin() const{
	return this->arccos(*this);
}

/* return |z| */
SMF_INLINE_FORCED float CComplex::abs() const {
	float x, y, t;
	x = CMath::fabs( r );
	y = CMath::fabs( i );
	if ( x == 0.0f ) {
		return y;
	} else if ( y == 0.0f ) {
		return x;
	} else if ( x > y ) {
		t = y / x;
		return x * CMath::sqrt( 1.0f + t * t );
	} else {
		t = x / y;
		return y * CMath::sqrt( 1.0f + t * t );
	}
}

SMF_INLINE_FORCED bool CComplex::compare( const CComplex &a ) const {
	return ( ( r == a.r ) && ( i == a.i ) );
}

SMF_INLINE_FORCED bool CComplex::compare( const CComplex &a, const float epsilon ) const {
	if ( CMath::fabs( r - a.r ) > epsilon ) {
		return false;
	}
	if ( CMath::fabs( i - a.i ) > epsilon ) {
		return false;
	}
	return true;
}

SMF_INLINE_FORCED bool CComplex::operator==( const CComplex &a ) const {
	return compare( a );
}

SMF_INLINE_FORCED bool CComplex::operator!=( const CComplex &a ) const {
	return !compare( a );
}

SMF_INLINE_FORCED int CComplex::getDimension() const {
	return 2;
}

SMF_INLINE_FORCED const float *CComplex::toFloatPtr() const {
	return &r;
}

SMF_INLINE_FORCED float *CComplex::toFloatPtr() {
	return &r;
}
} //end MATH
} //end SMF
#endif /* !__MATH_COMPLEX_H__ */
