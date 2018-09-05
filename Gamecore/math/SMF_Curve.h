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

#ifndef _SMF__MATH_CURVE_H__
#define _SMF__MATH_CURVE_H__
#include "../SMF_Config.h"
#include  "../structures/SMF_List.h"
#include "SMF_Matriz.h"

namespace SMF{
namespace MATH{
/*
===============================================================================

	Curve base template.

===============================================================================
*/
/**
 * \class  CCurve
 *
 * \ingroup SMF_Math
 *
 * \brief Template para Classe Curva. Cria uma lista de valores e tempos
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 
**/
template< class type >
class  SMF_API CCurve {
public:
	
	/**
	\brief construtor: cria um objeto CCurve
	**/
	CCurve();
	virtual			~CCurve();

	virtual int			addValue( const float time, const type &value );
	virtual void		removeIndex( const int index ) { values.removeIndex(index); times.removeIndex(index); changed = true; }
	virtual void		clear() { values.clear(); times.clear(); currentIndex = -1; changed = true; }
	/**
	\brief retorna o valor da curva no tempo solicitado
	\param time tempo no qual se deseja o valor da curva
	**/
	virtual type		getCurrentValue( const float time ) const;
	/**
	\brief retorna o valor da primeira derivada no tempo solicitado
	\param time tempo no qual se deseja o valor da derivada
	**/
	virtual type		getCurrentFirstDerivative( const float time ) const;
	/**
	\brief retorna o valor da segunda derivada no tempo solicitado
	\param time tempo no qual se deseja o valor da derivada
	**/
	virtual type		getCurrentSecondDerivative( const float time ) const;
	/**
	\brief verifica se o tempo passado está no intervalo de tempo da curva
	\return verdadeiro caso o tempo passado esteja dentro do iintervalo de tempo desta curve, falso, caso contrário
	**/
	virtual bool		isDone( const float time ) const;

	int					getNumValues() const { return values.Num(); }
	void				setValue( const int index, const type &value ) { values[index] = value; changed = true; }
	type				getValue( const int index ) const { return values[index]; }
	type *				getValueAddress( const int index ) { return &values[index]; }
	float				getTime( const int index ) const { return times[index]; }

	float				getLengthForTime( const float time ) const;
	float				getTimeForLength( const float length, const float epsilon = 0.1f ) const;
	float				getLengthBetweenKnots( const int i0, const int i1 ) const;

	void				makeUniform( const float totalTime );
	void				setConstantSpeed( const float totalTime );
	void				shiftTime( const float deltaTime );
	void				translate( const type &translation );

protected:
	/// knots
	CList<float>		times;
	/// knot values
	CList<type>			values;			
	/// cached index for fast lookup
	mutable int			currentIndex;	
	/// set whenever the curve changes
	mutable bool		changed;		

	int					indexForTime( const float time ) const;
	float				timeForIndex( const int index ) const;
	type				valueForIndex( const int index ) const;

	float				getSpeed( const float time ) const;
	/**
	\brief resolve a integral pelo método de Romberg

	\note http://pt.wikipedia.org/wiki/Método_de_Romberg
	\note http://www.sawp.com.br/blog/?p=1665
	**/
	float				rombergIntegral( const float t0, const float t1, const int order ) const;
};

/*
====================
CCurve::CCurve
====================
*/
template< class type >
SMF_INLINE_FORCED CCurve<type>::CCurve() {
	currentIndex = -1;
	changed = false;
}

/*
====================
CCurve::~CCurve
====================
*/
template< class type >
SMF_INLINE_FORCED CCurve<type>::~CCurve() {
}

/*
====================
CCurve::addValue

  add a timed/value pair to the spline
  returns the index to the inserted pair
====================
*/
template< class type >
SMF_INLINE_FORCED int CCurve<type>::addValue( const float time, const type &value ) {
	int i;

	i = indexForTime( time );
	times.insert( time, i );
	values.insert( value, i );
	changed = true;
	return i;
}

/*
====================
CCurve::GetCurrentValue

  get the value for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve<type>::getCurrentValue( const float time ) const {
	int i;

	i = indexForTime( time );
	if ( i >= values.Num() ) {
		return values[values.Num() - 1];
	} else {
		return values[i];
	}
}

/*
====================
CCurve::GetCurrentFirstDerivative

  get the first derivative for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve<type>::getCurrentFirstDerivative( const float time ) const {
	return ( values[0] - values[0] );
}

/*
====================
CCurve::GetCurrentSecondDerivative

  get the second derivative for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve<type>::getCurrentSecondDerivative( const float time ) const {
	return ( values[0] - values[0] );
}

/*
====================
CCurve::IsDone
====================
*/
template< class type >
SMF_INLINE_FORCED bool CCurve<type>::isDone( const float time ) const {
	return ( time >= times[ times.getNum() - 1 ] );
}

/*
====================
CCurve::getSpeed
====================
*/
template< class type >
SMF_INLINE_FORCED float CCurve<type>::getSpeed( const float time ) const {
	int i;
	float speed;
	type value;

	value = getCurrentFirstDerivative( time );
	for ( speed = 0.0f, i = 0; i < value.getDimension(); i++ ) {
		speed += value[i] * value[i];
	}
	return CMath::sqrt( speed );
}

/*
====================
CCurve::rombergIntegral
====================
*/
template< class type >
SMF_INLINE_FORCED float CCurve<type>::rombergIntegral( const float t0, const float t1, const int order ) const {
	int i, j, k, m, n;
	float sum, delta;
	float *temp[2];

	temp[0] = (float *) _alloca16( order * sizeof( float ) );
	temp[1] = (float *) _alloca16( order * sizeof( float ) );

	delta = t1 - t0;
	temp[0][0] = 0.5f * delta * ( getSpeed( t0 ) + getSpeed( t1 ) );

	for ( i = 2, m = 1; i <= order; i++, m *= 2, delta *= 0.5f ) {

		// approximate using the trapezoid rule
		sum = 0.0f;
		for ( j = 1; j <= m; j++ ) {
			sum += getSpeed( t0 + delta * ( j - 0.5f ) );
		}

		// Richardson extrapolation
		temp[1][0] = 0.5f * ( temp[0][0] + delta * sum );
		for ( k = 1, n = 4; k < i; k++, n *= 4 ) {
			temp[1][k] = ( n * temp[1][k-1] - temp[0][k-1] ) / ( n - 1 );
		}

		for ( j = 0; j < i; j++ ) {
			temp[0][j] = temp[1][j];
		}
	}
	return temp[0][order-1];
}

/*
====================
CCurve::getLengthBetweenKnots
====================
*/
template< class type >
SMF_INLINE_FORCED float CCurve<type>::getLengthBetweenKnots( const int i0, const int i1 ) const {
	float length = 0.0f;
	for ( int i = i0; i < i1; i++ ) {
		length += rombergIntegral( times[i], times[i+1], 5 );
	}
	return length;
}

/*
====================
CCurve::getLengthForTime
====================
*/
template< class type >
SMF_INLINE_FORCED float CCurve<type>::getLengthForTime( const float time ) const {
	float length = 0.0f;
	int index = indexForTime( time );
	for ( int i = 0; i < index; i++ ) {
		length += rombergIntegral( times[i], times[i+1], 5 );
	}
	length += rombergIntegral( times[index], time, 5 );
	return length;
}

/*
====================
CCurve::getTimeForLength
====================
*/
template< class type >
SMF_INLINE_FORCED float CCurve<type>::getTimeForLength( const float length, const float epsilon ) const {
	int i, index;
	float *accumLength, totalLength, len0, len1, t, diff;

	if ( length <= 0.0f ) {
		return times[0];
	}

	accumLength = (float *) _alloca16( values.Num() * sizeof( float ) );
	totalLength = 0.0f;
	for ( index = 0; index < values.Num() - 1; index++ ) {
		totalLength += getLengthBetweenKnots( index, index + 1 );
		accumLength[index] = totalLength;
		if ( length < accumLength[index] ) {
			break;
		}
	}

	if ( index >= values.Num() - 1 ) {
		return times[times.getNum() - 1];
	}

	if ( index == 0 ) {
		len0 = length;
		len1 = accumLength[0];
	} else {
		len0 = length - accumLength[index-1];
		len1 = accumLength[index] - accumLength[index-1];
	}

	// invert the arc length integral using Newton's method
	t = ( times[index+1] - times[index] ) * len0 / len1;
	for ( i = 0; i < 32; i++ ) {
		diff = rombergIntegral( times[index], times[index] + t, 5 ) - len0;
		if ( CMath::fabs( diff ) <= epsilon ) {
			return times[index] + t;
		}
		t -= diff / getSpeed( times[index] + t );
	}
	return times[index] + t;
}

/*
====================
CCurve::makeUniform
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve<type>::makeUniform( const float totalTime ) {
	int i, n;

	n = times.getNum() - 1;
	for ( i = 0; i <= n; i++ ) {
		times[i] = i * totalTime / n;
	}
	changed = true;
}

/*
====================
CCurve::setConstantSpeed
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve<type>::setConstantSpeed( const float totalTime ) {
	int i;
	float *length, totalLength, scale, t;

	length = (float *) _alloca16( values.Num() * sizeof( float ) );
	totalLength = 0.0f;
	for ( i = 0; i < values.Num() - 1; i++ ) {
		length[i] = getLengthBetweenKnots( i, i + 1 );
		totalLength += length[i];
	}
	scale = totalTime / totalLength;
	for ( t = 0.0f, i = 0; i < times.getNum() - 1; i++ ) {
		times[i] = t;
		t += scale * length[i];
	}
	times[times.getNum() - 1] = totalTime;
	changed = true;
}

/*
====================
CCurve::shiftTime
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve<type>::shiftTime( const float deltaTime ) {
	for ( int i = 0; i < times.getNum(); i++ ) {
		times[i] += deltaTime;
	}
	changed = true;
}

/*
====================
CCurve::translate
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve<type>::translate( const type &translation ) {
	for ( int i = 0; i < values.Num(); i++ ) {
		values[i] += translation;
	}
	changed = true;
}

/*
====================
CCurve::indexForTime

  find the index for the first time greater than or equal to the given time
====================
*/
template< class type >
SMF_INLINE_FORCED int CCurve<type>::indexForTime( const float time ) const {
	int len, mid, offset, res;

	if ( currentIndex >= 0 && currentIndex <= times.getNum() ) {
		// use the cached index if it is still valid
		if ( currentIndex == 0 ) {
			if ( time <= times[currentIndex] ) {
				return currentIndex;
			}
		} else if ( currentIndex == times.getNum() ) {
			if ( time > times[currentIndex-1] ) {
				return currentIndex;
			}
		} else if ( time > times[currentIndex-1] && time <= times[currentIndex] ) {
			return currentIndex;
		} else if ( time > times[currentIndex] && ( currentIndex+1 == times.getNum() || time <= times[currentIndex+1] ) ) {
			// use the next index
			currentIndex++;
			return currentIndex;
		}
	}

	// use binary search to find the index for the given time
	len = times.getNum();
	mid = len;
	offset = 0;
	res = 0;
	while( mid > 0 ) {
		mid = len >> 1;
		if ( time == times[offset+mid] ) {
			return offset+mid;
		} else if ( time > times[offset+mid] ) {
			offset += mid;
			len -= mid;
			res = 1;
		} else {
			len -= mid;
			res = 0;
		}
	}
	currentIndex = offset+res;
	return currentIndex;
}

/*
====================
CCurve::valueForIndex

  get the value for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve<type>::valueForIndex( const int index ) const {
	int n = values.Num()-1;

	if ( index < 0 ) {
		return values[0] + index * ( values[1] - values[0] );
	} else if ( index > n ) {
		return values[n] + ( index - n ) * ( values[n] - values[n-1] );
	}
	return values[index];
}

/*
====================
CCurve::timeForIndex

  get the value for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED float CCurve<type>::timeForIndex( const int index ) const {
	int n = times.getNum()-1;

	if ( index < 0 ) {
		return times[0] + index * ( times[1] - times[0] );
	} else if ( index > n ) {
		return times[n] + ( index - n ) * ( times[n] - times[n-1] );
	}
	return times[index];
}



/**
 * \class CCurve_Bezier
 *
 * \ingroup SMF_Math
 * 
 * \if pt_br
 * \brief Template de Curva de Bzier
 * \note O grau do polinômio é igual ao numero de nós menos um.
 * \elseif us_en
 * \brief 	Bezier Curve template.
 * \note	The degree of the polynomial equals the number of knots minus one.
 *
 * \endif
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \see http://pt.wikipedia.org/wiki/Curva_de_B%C3%A9zier
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
template< class type >
class CCurve_Bezier : public CCurve<type> {
public:
						CCurve_Bezier();

	virtual type		getCurrentValue( const float time ) const;
	virtual type		getCurrentFirstDerivative( const float time ) const;
	virtual type		getCurrentSecondDerivative( const float time ) const;

protected:
	void				Basis( const int order, const float t, float *bvals ) const;
	void				BasisFirstDerivative( const int order, const float t, float *bvals ) const;
	void				BasisSecondDerivative( const int order, const float t, float *bvals ) const;
};

/*
====================
CCurve_Bezier::CCurve_Bezier
====================
*/
template< class type >
SMF_INLINE_FORCED CCurve_Bezier<type>::CCurve_Bezier() {
}

/*
====================
CCurve_Bezier::GetCurrentValue

  get the value for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_Bezier<type>::getCurrentValue( const float time ) const {
	int i;
	float *bvals;
	type v;

	bvals = (float *) _alloca16( this->values.Num() * sizeof( float ) );

	Basis( this->values.Num(), time, bvals );
	v = bvals[0] * this->values[0];
	for ( i = 1; i < this->values.Num(); i++ ) {
		v += bvals[i] * this->values[i];
	}
	return v;
}

/*
====================
CCurve_Bezier::GetCurrentFirstDerivative

  get the first derivative for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_Bezier<type>::getCurrentFirstDerivative( const float time ) const {
	int i;
	float *bvals, d;
	type v;

	bvals = (float *) _alloca16( this->values.Num() * sizeof( float ) );

	BasisFirstDerivative( this->values.Num(), time, bvals );
	v = bvals[0] * this->values[0];
	for ( i = 1; i < this->values.Num(); i++ ) {
		v += bvals[i] * this->values[i];
	}
	d = ( this->times[this->times.Num()-1] - this->times[0] );
	return ( (float) (this->values.Num()-1) / d ) * v;
}

/*
====================
CCurve_Bezier::GetCurrentSecondDerivative

  get the second derivative for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_Bezier<type>::getCurrentSecondDerivative( const float time ) const {
	int i;
	float *bvals, d;
	type v;

	bvals = (float *) _alloca16( this->values.Num() * sizeof( float ) );

	BasisSecondDerivative( this->values.Num(), time, bvals );
	v = bvals[0] * this->values[0];
	for ( i = 1; i < this->values.Num(); i++ ) {
		v += bvals[i] * this->values[i];
	}
	d = ( this->times[this->times.Num()-1] - this->times[0] );
	return ( (float) (this->values.Num()-2) * (this->values.Num()-1) / ( d * d ) ) * v;
}

/*
====================
CCurve_Bezier::Basis

  bezier basis functions
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_Bezier<type>::Basis( const int order, const float t, float *bvals ) const {
	int i, j, d;
	float *c, c1, c2, s, o, ps, po;

	bvals[0] = 1.0f;
	d = order - 1;
	if ( d <= 0 ) {
		return;
	}

	c = (float *) _alloca16( (d+1) * sizeof( float ) );
	s = (float) ( t - this->times[0] ) / ( this->times[this->times.Num()-1] - this->times[0] );
    o = 1.0f - s;
	ps = s;
	po = o;

	for ( i = 1; i < d; i++ ) {
		c[i] = 1.0f;
	}
	for ( i = 1; i < d; i++ ) {
		c[i-1] = 0.0f;
		c1 = c[i];
		c[i] = 1.0f;
		for ( j = i+1; j <= d; j++ ) {
			c2 = c[j];
			c[j] = c1 + c[j-1];
			c1 = c2;
		}
		bvals[i] = c[d] * ps;
		ps *= s;
	}
	for ( i = d-1; i >= 0; i-- ) {
		bvals[i] *= po;
		po *= o;
	}
	bvals[d] = ps;
}

/*
====================
CCurve_Bezier::BasisFirstDerivative

  first derivative of bezier basis functions
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_Bezier<type>::BasisFirstDerivative( const int order, const float t, float *bvals ) const {
	int i;

	Basis( order-1, t, bvals+1 );
	bvals[0] = 0.0f;
	for ( i = 0; i < order-1; i++ ) {
		bvals[i] -= bvals[i+1];
	}
}

/*
====================
CCurve_Bezier::BasisSecondDerivative

  second derivative of bezier basis functions
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_Bezier<type>::BasisSecondDerivative( const int order, const float t, float *bvals ) const {
	int i;

	BasisFirstDerivative( order-1, t, bvals+1 );
	bvals[0] = 0.0f;
	for ( i = 0; i < order-1; i++ ) {
		bvals[i] -= bvals[i+1];
	}
}


/**
 * \class CCurve_QuadraticBezier
 *
 * \ingroup SMF_Math
 * 
 * \if pt_br
 * \brief Implementa Curve de Bezier Quadrática
 * \elseif us_en
 * \brief 	Quadratic Bezier Curve template.
 *	\note   Should always have exactly three knots.
 *
 * \endif
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \see http://en.wikipedia.org/wiki/B%C3%A9zier_curve#Quadratic_curves
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
template< class type >
class CCurve_QuadraticBezier : public CCurve<type> {

public:
						CCurve_QuadraticBezier();

	virtual type		getCurrentValue( const float time ) const;
	virtual type		getCurrentFirstDerivative( const float time ) const;
	virtual type		getCurrentSecondDerivative( const float time ) const;

protected:
	void				Basis( const float t, float *bvals ) const;
	void				BasisFirstDerivative( const float t, float *bvals ) const;
	void				BasisSecondDerivative( const float t, float *bvals ) const;
};

/*
====================
CCurve_QuadraticBezier::CCurve_QuadraticBezier
====================
*/
template< class type >
SMF_INLINE_FORCED CCurve_QuadraticBezier<type>::CCurve_QuadraticBezier() {
}


/*
====================
CCurve_QuadraticBezier::GetCurrentValue

  get the value for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_QuadraticBezier<type>::getCurrentValue( const float time ) const {
	float bvals[3];
	SMF_ASSERT( this->values.Num() == 3 );
	Basis( time, bvals );
	return ( bvals[0] * this->values[0] + bvals[1] * this->values[1] + bvals[2] * this->values[2] );
}

/*
====================
CCurve_QuadraticBezier::GetCurrentFirstDerivative

  get the first derivative for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_QuadraticBezier<type>::getCurrentFirstDerivative( const float time ) const {
	float bvals[3], d;
	SMF_ASSERT( this->values.Num() == 3 );
	BasisFirstDerivative( time, bvals );
	d = ( this->times[2] - this->times[0] );
	return ( bvals[0] * this->values[0] + bvals[1] * this->values[1] + bvals[2] * this->values[2] ) / d;
}

/*
====================
CCurve_QuadraticBezier::GetCurrentSecondDerivative

  get the second derivative for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_QuadraticBezier<type>::getCurrentSecondDerivative( const float time ) const {
	float bvals[3], d;
	SMF_ASSERT( this->values.Num() == 3 );
	BasisSecondDerivative( time, bvals );
	d = ( this->times[2] - this->times[0] );
	return ( bvals[0] * this->values[0] + bvals[1] * this->values[1] + bvals[2] * this->values[2] ) / ( d * d );
}

/*
====================
CCurve_QuadraticBezier::Basis

  quadratic bezier basis functions
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_QuadraticBezier<type>::Basis( const float t, float *bvals ) const {
	float s1 = (float) ( t - this->times[0] ) / ( this->times[2] - this->times[0] );
	float s2 = s1 * s1;
	bvals[0] = s2 - 2.0f * s1 + 1.0f;
	bvals[1] = -2.0f * s2 + 2.0f * s1;
	bvals[2] = s2;
}

/*
====================
CCurve_QuadraticBezier::BasisFirstDerivative

  first derivative of quadratic bezier basis functions
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_QuadraticBezier<type>::BasisFirstDerivative( const float t, float *bvals ) const {
	float s1 = (float) ( t - this->times[0] ) / ( this->times[2] - this->times[0] );
	bvals[0] = 2.0f * s1 - 2.0f;
	bvals[1] = -4.0f * s1 + 2.0f;
	bvals[2] = 2.0f * s1;
}

/*
====================
CCurve_QuadraticBezier::BasisSecondDerivative

  second derivative of quadratic bezier basis functions
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_QuadraticBezier<type>::BasisSecondDerivative( const float t, float *bvals ) const {
//	float s1 = (float) ( t - this->times[0] ) / ( this->times[2] - this->times[0] );
	bvals[0] = 2.0f;
	bvals[1] = -4.0f;
	bvals[2] = 2.0f;
}



/**
 * \class CCurve_CubicBezier
 *
 * \ingroup SMF_Math
 * 
 * \if pt_br
 * \brief Template de Curva de Bezier Cúbica
 * \note O grau do polinômio é igual ao numero de nós menos um.
 * \elseif us_en
 * \brief 	Cubic Bezier Curve template.
 *	\note   Should always have exactly four knots.
 *
 * \endif
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \see http://en.wikipedia.org/wiki/B%C3%A9zier_curve#Cubic_B.C3.A9zier_curves
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
template< class type >
class CCurve_CubicBezier : public CCurve<type> {

public:
						CCurve_CubicBezier();

	virtual type		getCurrentValue( const float time ) const;
	virtual type		getCurrentFirstDerivative( const float time ) const;
	virtual type		getCurrentSecondDerivative( const float time ) const;

protected:
	void				Basis( const float t, float *bvals ) const;
	void				BasisFirstDerivative( const float t, float *bvals ) const;
	void				BasisSecondDerivative( const float t, float *bvals ) const;
};

/*
====================
CCurve_CubicBezier::CCurve_CubicBezier
====================
*/
template< class type >
SMF_INLINE_FORCED CCurve_CubicBezier<type>::CCurve_CubicBezier() {
}


/*
====================
CCurve_CubicBezier::GetCurrentValue

  get the value for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_CubicBezier<type>::getCurrentValue( const float time ) const {
	float bvals[4];
	SMF_ASSERT( this->values.Num() == 4 );
	Basis( time, bvals );
	return ( bvals[0] * this->values[0] + bvals[1] * this->values[1] + bvals[2] * this->values[2] + bvals[3] * this->values[3] );
}

/*
====================
CCurve_CubicBezier::GetCurrentFirstDerivative

  get the first derivative for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_CubicBezier<type>::getCurrentFirstDerivative( const float time ) const {
	float bvals[4], d;
	SMF_ASSERT( this->values.Num() == 4 );
	BasisFirstDerivative( time, bvals );
	d = ( this->times[3] - this->times[0] );
	return ( bvals[0] * this->values[0] + bvals[1] * this->values[1] + bvals[2] * this->values[2] + bvals[3] * this->values[3] ) / d;
}

/*
====================
CCurve_CubicBezier::GetCurrentSecondDerivative

  get the second derivative for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_CubicBezier<type>::getCurrentSecondDerivative( const float time ) const {
	float bvals[4], d;
	SMF_ASSERT( this->values.Num() == 4 );
	BasisSecondDerivative( time, bvals );
	d = ( this->times[3] - this->times[0] );
	return ( bvals[0] * this->values[0] + bvals[1] * this->values[1] + bvals[2] * this->values[2] + bvals[3] * this->values[3] ) / ( d * d );
}

/*
====================
CCurve_CubicBezier::Basis

  cubic bezier basis functions
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_CubicBezier<type>::Basis( const float t, float *bvals ) const {
	float s1 = (float) ( t - this->times[0] ) / ( this->times[3] - this->times[0] );
	float s2 = s1 * s1;
	float s3 = s2 * s1;
	bvals[0] = -s3 + 3.0f * s2 - 3.0f * s1 + 1.0f;
	bvals[1] = 3.0f * s3 - 6.0f * s2 + 3.0f * s1;
	bvals[2] = -3.0f * s3 + 3.0f * s2;
	bvals[3] = s3;
}

/*
====================
CCurve_CubicBezier::BasisFirstDerivative

  first derivative of cubic bezier basis functions
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_CubicBezier<type>::BasisFirstDerivative( const float t, float *bvals ) const {
	float s1 = (float) ( t - this->times[0] ) / ( this->times[3] - this->times[0] );
	float s2 = s1 * s1;
	bvals[0] = -3.0f * s2 + 6.0f * s1 - 3.0f;
	bvals[1] = 9.0f * s2 - 12.0f * s1 + 3.0f;
	bvals[2] = -9.0f * s2 + 6.0f * s1;
	bvals[3] = 3.0f * s2;
}

/*
====================
CCurve_CubicBezier::BasisSecondDerivative

  second derivative of cubic bezier basis functions
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_CubicBezier<type>::BasisSecondDerivative( const float t, float *bvals ) const {
	float s1 = (float) ( t - this->times[0] ) / ( this->times[3] - this->times[0] );
	bvals[0] = -6.0f * s1 + 6.0f;
	bvals[1] = 18.0f * s1 - 12.0f;
	bvals[2] = -18.0f * s1 + 6.0f;
	bvals[3] = 6.0f * s1;
}


/*
===============================================================================

	

===============================================================================
*/


/**
 * \class CCurve_Spline
 *
 * \ingroup SMF_Math
 * 
 * \if pt_br
 * \brief	Template base de curva Spline
 * \elseif us_en
 * \brief 	Spline base template.
 * \endif
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \see http://en.wikipedia.org/wiki/Spline_%28mathematics%29, http://pt.wikipedia.org/wiki/Spline
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
template< class type >
class CCurve_Spline : public CCurve<type> {

public:
	enum				boundary_t { BT_FREE, BT_CLAMPED, BT_CLOSED };

						CCurve_Spline();

	virtual bool		isDone( const float time ) const;

	virtual void		SetBoundaryType( const boundary_t bt ) { boundaryType = bt; this->changed = true; }
	virtual boundary_t	GetBoundaryType() const { return boundaryType; }

	virtual void		SetCloseTime( const float t ) { closeTime = t; this->changed = true; }
	virtual float		GetCloseTime() { return boundaryType == BT_CLOSED ? closeTime : 0.0f; }

protected:
	boundary_t			boundaryType;
	float				closeTime;

	type				valueForIndex( const int index ) const;
	float				timeForIndex( const int index ) const;
	float				ClampedTime( const float t ) const;
};

/*
====================
CCurve_Spline::CCurve_Spline
====================
*/
template< class type >
SMF_INLINE_FORCED CCurve_Spline<type>::CCurve_Spline() {
	boundaryType = BT_FREE;
	closeTime = 0.0f;
}

/*
====================
CCurve_Spline::valueForIndex

  get the value for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_Spline<type>::valueForIndex( const int index ) const {
	int n = this->values.Num()-1;

	if ( index < 0 ) {
		if ( boundaryType == BT_CLOSED ) {
			return this->values[ this->values.Num() + index % this->values.Num() ];
		}
		else {
			return this->values[0] + index * ( this->values[1] - this->values[0] );
		}
	}
	else if ( index > n ) {
		if ( boundaryType == BT_CLOSED ) {
			return this->values[ index % this->values.Num() ];
		}
		else {
			return this->values[n] + ( index - n ) * ( this->values[n] - this->values[n-1] );
		}
	}
	return this->values[index];
}

/*
====================
CCurve_Spline::timeForIndex

  get the value for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED float CCurve_Spline<type>::timeForIndex( const int index ) const {
	int n = this->times.Num()-1;

	if ( index < 0 ) {
		if ( boundaryType == BT_CLOSED ) {
			return ( index / this->times.Num() ) * ( this->times[n] + closeTime ) - ( this->times[n] + closeTime - this->times[this->times.Num() + index % this->times.Num()] );
		}
		else {
			return this->times[0] + index * ( this->times[1] - this->times[0] );
		}
	}
	else if ( index > n ) {
		if ( boundaryType == BT_CLOSED ) {
			return ( index / this->times.Num() ) * ( this->times[n] + closeTime ) + this->times[index % this->times.Num()];
		}
		else {
			return this->times[n] + ( index - n ) * ( this->times[n] - this->times[n-1] );
		}
	}
	return this->times[index];
}

/*
====================
CCurve_Spline::ClampedTime

  return the clamped time based on the boundary type
====================
*/
template< class type >
SMF_INLINE_FORCED float CCurve_Spline<type>::ClampedTime( const float t ) const {
	if ( boundaryType == BT_CLAMPED ) {
		if ( t < this->times[0] ) {
			return this->times[0];
		}
		else if ( t >= this->times[this->times.Num()-1] ) {
			return this->times[this->times.Num()-1];
		}
	}
	return t;
}

/*
====================
CCurve_Spline::IsDone
====================
*/
template< class type >
SMF_INLINE_FORCED bool CCurve_Spline<type>::isDone( const float time ) const {
	return ( boundaryType != BT_CLOSED && time >= this->times[ this->times.Num() - 1 ] );
}


/*
===============================================================================


===============================================================================
*/


/**
 * \class CCurve_NaturalCubicSpline
 *
 * \ingroup SMF_Math
 * 
 * \if pt_br
 * \brief	Template curva Spline Cubica Interpolar
 * \elseif us_en
 * \brief 	Cubic Interpolating Spline template.
 * \note	The curve goes through all the knots.
 * \endif
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \see http://en.wikipedia.org/wiki/Spline_%28mathematics%29, http://pt.wikipedia.org/wiki/Spline
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
template< class type >
class CCurve_NaturalCubicSpline : public CCurve_Spline<type> {
public:
						CCurve_NaturalCubicSpline();

	virtual void		clear() { CCurve_Spline<type>::clear(); this->values.clear(); b.clear(); c.clear(); d.clear(); }

	virtual type		getCurrentValue( const float time ) const;
	virtual type		getCurrentFirstDerivative( const float time ) const;
	virtual type		getCurrentSecondDerivative( const float time ) const;

protected:
	mutable CList<type>b;
	mutable CList<type>c;
	mutable CList<type>d;

	void				Setup() const;
	void				SetupFree() const;
	void				SetupClamped() const;
	void				SetupClosed() const;
};

/*
====================
CCurve_NaturalCubicSpline::CCurve_NaturalCubicSpline
====================
*/
template< class type >
SMF_INLINE_FORCED CCurve_NaturalCubicSpline<type>::CCurve_NaturalCubicSpline() {
}

/*
====================
CCurve_NaturalCubicSpline::GetCurrentValue

  get the value for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_NaturalCubicSpline<type>::getCurrentValue( const float time ) const {
	float clampedTime = this->ClampedTime( time );
	int i = this->indexForTime( clampedTime );
	float s = time - this->timeForIndex( i );
	Setup();
	return ( this->values[i] + s * ( b[i] + s * ( c[i] + s * d[i] ) ) );
}

/*
====================
CCurve_NaturalCubicSpline::GetCurrentFirstDerivative

  get the first derivative for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_NaturalCubicSpline<type>::getCurrentFirstDerivative( const float time ) const {
	float clampedTime = this->ClampedTime( time );
	int i = this->indexForTime( clampedTime );
	float s = time - this->timeForIndex( i );
	Setup();
	return ( b[i] + s * ( 2.0f * c[i] + 3.0f * s * d[i] ) );
}


/*
====================
CCurve_NaturalCubicSpline::GetCurrentSecondDerivative

  get the second derivative for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_NaturalCubicSpline<type>::getCurrentSecondDerivative( const float time ) const {
	float clampedTime = this->ClampedTime( time );
	int i = this->indexForTime( clampedTime );
	float s = time - this->timeForIndex( i );
	Setup();
	return ( 2.0f * c[i] + 6.0f * s * d[i] );
}

/*
====================
CCurve_NaturalCubicSpline::Setup
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_NaturalCubicSpline<type>::Setup() const {
	if ( this->changed ) {
		switch( this->boundaryType ) {
			case CCurve_Spline<type>::BT_FREE:		SetupFree(); break;
			case CCurve_Spline<type>::BT_CLAMPED:	SetupClamped(); break;
			case CCurve_Spline<type>::BT_CLOSED:		SetupClosed(); break;
		}
		this->changed = false;
	}
}

/*
====================
CCurve_NaturalCubicSpline::SetupFree
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_NaturalCubicSpline<type>::SetupFree() const {
	int i;
	float inv;
	float *d0, *d1, *beta, *gamma;
	type *alpha, *delta;

	d0 = (float *) _alloca16( ( this->values.Num() - 1 ) * sizeof( float ) );
	d1 = (float *) _alloca16( ( this->values.Num() - 1 ) * sizeof( float ) );
	alpha = (type *) _alloca16( ( this->values.Num() - 1 ) * sizeof( type ) );
	beta = (float *) _alloca16( this->values.Num() * sizeof( float ) );
	gamma = (float *) _alloca16( ( this->values.Num() - 1 ) * sizeof( float ) );
	delta = (type *) _alloca16( this->values.Num() * sizeof( type ) );

	for ( i = 0; i < this->values.Num() - 1; i++ ) {
		d0[i] = this->times[i+1] - this->times[i];
	}

	for ( i = 1; i < this->values.Num() - 1; i++ ) {
		d1[i] = this->times[i+1] - this->times[i-1];
	}

	for ( i = 1; i < this->values.Num() - 1; i++ ) {
		type sum = 3.0f * ( d0[i-1] * this->values[i+1] - d1[i] * this->values[i] + d0[i] * this->values[i-1] );
		inv = 1.0f / ( d0[i-1] * d0[i] );
		alpha[i] = inv * sum;
	}

	beta[0] = 1.0f;
	gamma[0] = 0.0f;
	delta[0] = this->values[0] - this->values[0];

	for ( i = 1; i < this->values.Num() - 1; i++ ) {
		beta[i] = 2.0f * d1[i] - d0[i-1] * gamma[i-1];
		inv = 1.0f / beta[i];
		gamma[i] = inv * d0[i];
		delta[i] = inv * ( alpha[i] - d0[i-1] * delta[i-1] );
	}
	beta[this->values.Num() - 1] = 1.0f;
	delta[this->values.Num() - 1] = this->values[0] - this->values[0];

	b.AssureSize( this->values.Num() );
	c.AssureSize( this->values.Num() );
	d.AssureSize( this->values.Num() );

	c[this->values.Num() - 1] = this->values[0] - this->values[0];

	for ( i = this->values.Num() - 2; i >= 0; i-- ) {
		c[i] = delta[i] - gamma[i] * c[i+1];
		inv = 1.0f / d0[i];
		b[i] = inv * ( this->values[i+1] - this->values[i] ) - ( 1.0f / 3.0f ) * d0[i] * ( c[i+1] + 2.0f * c[i] );
		d[i] = ( 1.0f / 3.0f ) * inv * ( c[i+1] - c[i] );
	}
}

/*
====================
CCurve_NaturalCubicSpline::SetupClamped
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_NaturalCubicSpline<type>::SetupClamped() const {
	int i;
	float inv;
	float *d0, *d1, *beta, *gamma;
	type *alpha, *delta;

	d0 = (float *) _alloca16( ( this->values.Num() - 1 ) * sizeof( float ) );
	d1 = (float *) _alloca16( ( this->values.Num() - 1 ) * sizeof( float ) );
	alpha = (type *) _alloca16( ( this->values.Num() - 1 ) * sizeof( type ) );
	beta = (float *) _alloca16( this->values.Num() * sizeof( float ) );
	gamma = (float *) _alloca16( ( this->values.Num() - 1 ) * sizeof( float ) );
	delta = (type *) _alloca16( this->values.Num() * sizeof( type ) );

	for ( i = 0; i < this->values.Num() - 1; i++ ) {
		d0[i] = this->times[i+1] - this->times[i];
	}

	for ( i = 1; i < this->values.Num() - 1; i++ ) {
		d1[i] = this->times[i+1] - this->times[i-1];
	}

	inv = 1.0f / d0[0];
	alpha[0] = 3.0f * ( inv - 1.0f ) * ( this->values[1] - this->values[0] );
	inv = 1.0f / d0[this->values.Num() - 2];
	alpha[this->values.Num() - 1] = 3.0f * ( 1.0f - inv ) * ( this->values[this->values.Num() - 1] - this->values[this->values.Num() - 2] );

	for ( i = 1; i < this->values.Num() - 1; i++ ) {
		type sum = 3.0f * ( d0[i-1] * this->values[i+1] - d1[i] * this->values[i] + d0[i] * this->values[i-1] );
		inv = 1.0f / ( d0[i-1] * d0[i] );
		alpha[i] = inv * sum;
	}

	beta[0] = 2.0f * d0[0];
	gamma[0] = 0.5f;
	inv = 1.0f / beta[0];
	delta[0] = inv * alpha[0];

	for ( i = 1; i < this->values.Num() - 1; i++ ) {
		beta[i] = 2.0f * d1[i] - d0[i-1] * gamma[i-1];
		inv = 1.0f / beta[i];
		gamma[i] = inv * d0[i];
		delta[i] = inv * ( alpha[i] - d0[i-1] * delta[i-1] );
	}

	beta[this->values.Num() - 1] = d0[this->values.Num() - 2] * ( 2.0f - gamma[this->values.Num() - 2] );
	inv = 1.0f / beta[this->values.Num() - 1];
	delta[this->values.Num() - 1] = inv * ( alpha[this->values.Num() - 1] - d0[this->values.Num() - 2] * delta[this->values.Num() - 2] );

	b.AssureSize( this->values.Num() );
	c.AssureSize( this->values.Num() );
	d.AssureSize( this->values.Num() );

	c[this->values.Num() - 1] = delta[this->values.Num() - 1];

	for ( i = this->values.Num() - 2; i >= 0; i-- ) {
		c[i] = delta[i] - gamma[i] * c[i+1];
		inv = 1.0f / d0[i];
		b[i] = inv * ( this->values[i+1] - this->values[i] ) - ( 1.0f / 3.0f ) * d0[i]* ( c[i+1] + 2.0f * c[i] );
		d[i] = ( 1.0f / 3.0f ) * inv * ( c[i+1] - c[i] );
	}
}

/*
====================
CCurve_NaturalCubicSpline::SetupClosed
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_NaturalCubicSpline<type>::SetupClosed() const {
	int i, j;
	float c0, c1;
	float *d0;
	CMatXD mat;
	CVecXD x;

	d0 = (float *) _alloca16( ( this->values.Num() - 1 ) * sizeof( float ) );
	x.setData( this->values.Num(), VECX_ALLOCA( this->values.Num() ) );
	mat.setData( this->values.Num(), this->values.Num(), MATX_ALLOCA( this->values.Num() * this->values.Num() ) );

	b.AssureSize( this->values.Num() );
	c.AssureSize( this->values.Num() );
	d.AssureSize( this->values.Num() );

	for ( i = 0; i < this->values.Num() - 1; i++ ) {
		d0[i] = this->times[i+1] - this->times[i];
	}

	// matrix of system
	mat[0][0] = 1.0f;
	mat[0][this->values.Num() - 1] = -1.0f;
	for ( i = 1; i <= this->values.Num() - 2; i++ ) {
		mat[i][i-1] = d0[i-1];
		mat[i][i  ] = 2.0f * ( d0[i-1] + d0[i] );
		mat[i][i+1] = d0[i];
	}
	mat[this->values.Num() - 1][this->values.Num() - 2] = d0[this->values.Num() - 2];
	mat[this->values.Num() - 1][0] = 2.0f * ( d0[this->values.Num() - 2] + d0[0] );
	mat[this->values.Num() - 1][1] = d0[0];

	// right-hand side
	c[0].toZero();
	for ( i = 1; i <= this->values.Num() - 2; i++ ) {
		c0 = 1.0f / d0[i];
		c1 = 1.0f / d0[i-1];
		c[i] = 3.0f * ( c0 * ( this->values[i + 1] - this->values[i] ) - c1 * ( this->values[i] - this->values[i - 1] ) );
	}
	c0 = 1.0f / d0[0];
	c1 = 1.0f / d0[this->values.Num() - 2];
	c[this->values.Num() - 1] = 3.0f * ( c0 * ( this->values[1] - this->values[0] ) - c1 * ( this->values[0] - this->values[this->values.Num() - 2] ) );

	// solve system for each dimension
	mat.lu_Factor( NULL );
	for ( i = 0; i < this->values[0].getDimension(); i++ ) {
		for ( j = 0; j < this->values.Num(); j++ ) {
			x[j] = c[j][i];
		}
		mat.lu_solve( x, x, NULL );
		for ( j = 0; j < this->values.Num(); j++ ) {
			c[j][i] = x[j];
		}
	}

	for ( i = 0; i < this->values.Num() - 1; i++ ) {
		c0 = 1.0f / d0[i];
		b[i] = c0 * ( this->values[i + 1] - this->values[i] ) - ( 1.0f / 3.0f ) * ( c[i+1] + 2.0f * c[i] ) * d0[i];
		d[i] = ( 1.0f / 3.0f ) * c0 * ( c[i + 1] - c[i] );
	}
}



/**
 * \class CCurve_CatmullRomSpline
 *
 * \ingroup SMF_Math
 * 
 * \if pt_br
 * \brief	Template de curva Spline Cubica Uniforme
 * \elseif us_en
 * \brief 	Uniform Cubic Interpolating Spline template.
 * \note	The curve goes through all the knots.
 * \endif
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \see http://en.wikipedia.org/wiki/Spline_%28mathematics%29, http://pt.wikipedia.org/wiki/Spline
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
template< class type >
class CCurve_CatmullRomSpline : public CCurve_Spline<type> {

public:
						CCurve_CatmullRomSpline();

	virtual type		getCurrentValue( const float time ) const;
	virtual type		getCurrentFirstDerivative( const float time ) const;
	virtual type		getCurrentSecondDerivative( const float time ) const;

protected:
	void				Basis( const int index, const float t, float *bvals ) const;
	void				BasisFirstDerivative( const int index, const float t, float *bvals ) const;
	void				BasisSecondDerivative( const int index, const float t, float *bvals ) const;
};

/*
====================
CCurve_CatmullRomSpline::CCurve_CatmullRomSpline
====================
*/
template< class type >
SMF_INLINE_FORCED CCurve_CatmullRomSpline<type>::CCurve_CatmullRomSpline() {
}

/*
====================
CCurve_CatmullRomSpline::GetCurrentValue

  get the value for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_CatmullRomSpline<type>::getCurrentValue( const float time ) const {
	int i, j, k;
	float bvals[4], clampedTime;
	type v;

	if ( this->times.Num() == 1 ) {
		return this->values[0];
	}

	clampedTime = this->ClampedTime( time );
	i = this->indexForTime( clampedTime );
	Basis( i-1, clampedTime, bvals );
	v = this->values[0] - this->values[0];
	for ( j = 0; j < 4; j++ ) {
		k = i + j - 2;
		v += bvals[j] * this->valueForIndex( k );
	}
	return v;
}

/*
====================
CCurve_CatmullRomSpline::GetCurrentFirstDerivative

  get the first derivative for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_CatmullRomSpline<type>::getCurrentFirstDerivative( const float time ) const {
	int i, j, k;
	float bvals[4], d, clampedTime;
	type v;

	if ( this->times.Num() == 1 ) {
		return ( this->values[0] - this->values[0] );
	}

	clampedTime = this->ClampedTime( time );
	i = this->indexForTime( clampedTime );
	BasisFirstDerivative( i-1, clampedTime, bvals );
	v = this->values[0] - this->values[0];
	for ( j = 0; j < 4; j++ ) {
		k = i + j - 2;
		v += bvals[j] * this->valueForIndex( k );
	}
	d = ( this->timeForIndex( i ) - this->timeForIndex( i-1 ) );
	return v / d;
}

/*
====================
CCurve_CatmullRomSpline::GetCurrentSecondDerivative

  get the second derivative for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_CatmullRomSpline<type>::getCurrentSecondDerivative( const float time ) const {
	int i, j, k;
	float bvals[4], d, clampedTime;
	type v;

	if ( this->times.Num() == 1 ) {
		return ( this->values[0] - this->values[0] );
	}

	clampedTime = this->ClampedTime( time );
	i = this->indexForTime( clampedTime );
	BasisSecondDerivative( i-1, clampedTime, bvals );
	v = this->values[0] - this->values[0];
	for ( j = 0; j < 4; j++ ) {
		k = i + j - 2;
		v += bvals[j] * this->valueForIndex( k );
	}
	d = ( this->timeForIndex( i ) - this->timeForIndex( i-1 ) );
	return v / ( d * d );
}

/*
====================
CCurve_CatmullRomSpline::Basis

  spline basis functions
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_CatmullRomSpline<type>::Basis( const int index, const float t, float *bvals ) const {
	float s = (float) ( t - this->timeForIndex( index ) ) / ( this->timeForIndex( index+1 ) - this->timeForIndex( index ) );
	bvals[0] = ( ( -s + 2.0f ) * s - 1.0f ) * s * 0.5f;				// -0.5f s * s * s + s * s - 0.5f * s
	bvals[1] = ( ( ( 3.0f * s - 5.0f ) * s ) * s + 2.0f ) * 0.5f;	// 1.5f * s * s * s - 2.5f * s * s + 1.0f
	bvals[2] = ( ( -3.0f * s + 4.0f ) * s + 1.0f ) * s * 0.5f;		// -1.5f * s * s * s - 2.0f * s * s + 0.5f s
	bvals[3] = ( ( s - 1.0f ) * s * s ) * 0.5f;						// 0.5f * s * s * s - 0.5f * s * s
}

/*
====================
CCurve_CatmullRomSpline::BasisFirstDerivative

  first derivative of spline basis functions
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_CatmullRomSpline<type>::BasisFirstDerivative( const int index, const float t, float *bvals ) const {
	float s = (float) ( t - this->timeForIndex( index ) ) / ( this->timeForIndex( index+1 ) - this->timeForIndex( index ) );
	bvals[0] = ( -1.5f * s + 2.0f ) * s - 0.5f;						// -1.5f * s * s + 2.0f * s - 0.5f
	bvals[1] = ( 4.5f * s - 5.0f ) * s;								// 4.5f * s * s - 5.0f * s
	bvals[2] = ( -4.5 * s + 4.0f ) * s + 0.5f;						// -4.5 * s * s + 4.0f * s + 0.5f
	bvals[3] = 1.5f * s * s - s;									// 1.5f * s * s - s
}

/*
====================
CCurve_CatmullRomSpline::BasisSecondDerivative

  second derivative of spline basis functions
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_CatmullRomSpline<type>::BasisSecondDerivative( const int index, const float t, float *bvals ) const {
	float s = (float) ( t - this->timeForIndex( index ) ) / ( this->timeForIndex( index+1 ) - this->timeForIndex( index ) );
	bvals[0] = -3.0f * s + 2.0f;
	bvals[1] = 9.0f * s - 5.0f;
	bvals[2] = -9.0f * s + 4.0f;
	bvals[3] = 3.0f * s - 1.0f;
}



/**
 * \class CCurve_KochanekBartelsSpline
 *
 * \ingroup SMF_Math
 * 
 * \if pt_br
 * \brief	Template de curva Spline Cubica Interpolar
 * \elseif us_en
 * \brief 	Cubic Interpolating Spline template.
 * \note	The curve goes through all the knots.
 * \note	The curve becomes the Catmull-Rom spline if the tension,
 * \note	continuity and bias are all set to zero.
 * \endif
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \see http://en.wikipedia.org/wiki/Spline_%28mathematics%29, http://pt.wikipedia.org/wiki/Spline
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
template< class type >
class CCurve_KochanekBartelsSpline : public CCurve_Spline<type> {

public:
						CCurve_KochanekBartelsSpline();

	virtual int			addValue( const float time, const type &value );
	virtual int			addValue( const float time, const type &value, const float tension, const float continuity, const float bias );
	virtual void		removeIndex( const int index ) { this->values.removeIndex(index); this->times.removeIndex(index); tension.removeIndex(index); continuity.removeIndex(index); bias.removeIndex(index); }
	virtual void		clear() { this->values.clear(); this->times.clear(); tension.clear(); continuity.clear(); bias.clear(); this->currentIndex = -1; }

	virtual type		getCurrentValue( const float time ) const;
	virtual type		getCurrentFirstDerivative( const float time ) const;
	virtual type		getCurrentSecondDerivative( const float time ) const;

protected:
	CList<float>		tension;
	CList<float>		continuity;
	CList<float>		bias;

	void				TangentsForIndex( const int index, type &t0, type &t1 ) const;

	void				Basis( const int index, const float t, float *bvals ) const;
	void				BasisFirstDerivative( const int index, const float t, float *bvals ) const;
	void				BasisSecondDerivative( const int index, const float t, float *bvals ) const;
};

/*
====================
CCurve_KochanekBartelsSpline::CCurve_KochanekBartelsSpline
====================
*/
template< class type >
SMF_INLINE_FORCED CCurve_KochanekBartelsSpline<type>::CCurve_KochanekBartelsSpline() {
}

/*
====================
CCurve_KochanekBartelsSpline::addValue

  add a timed/value pair to the spline
  returns the index to the inserted pair
====================
*/
template< class type >
SMF_INLINE_FORCED int CCurve_KochanekBartelsSpline<type>::addValue( const float time, const type &value ) {
	int i;

	i = this->indexForTime( time );
	this->times.insert( time, i );
	this->values.insert( value, i );
	tension.insert( 0.0f, i );
	continuity.insert( 0.0f, i );
	bias.insert( 0.0f, i );
	return i;
}

/*
====================
CCurve_KochanekBartelsSpline::addValue

  add a timed/value pair to the spline
  returns the index to the inserted pair
====================
*/
template< class type >
SMF_INLINE_FORCED int CCurve_KochanekBartelsSpline<type>::addValue( const float time, const type &value, const float tension, const float continuity, const float bias ) {
	int i;

	i = this->indexForTime( time );
	this->times.insert( time, i );
	this->values.insert( value, i );
	this->tension.insert( tension, i );
	this->continuity.insert( continuity, i );
	this->bias.insert( bias, i );
	return i;
}

/*
====================
CCurve_KochanekBartelsSpline::GetCurrentValue

  get the value for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_KochanekBartelsSpline<type>::getCurrentValue( const float time ) const {
	int i;
	float bvals[4], clampedTime;
	type v, t0, t1;

	if ( this->times.Num() == 1 ) {
		return this->values[0];
	}

	clampedTime = this->ClampedTime( time );
	i = this->indexForTime( clampedTime );
	TangentsForIndex( i - 1, t0, t1 );
	Basis( i - 1, clampedTime, bvals );
	v = bvals[0] * this->valueForIndex( i - 1 );
	v += bvals[1] * this->valueForIndex( i );
	v += bvals[2] * t0;
	v += bvals[3] * t1;
	return v;
}

/*
====================
CCurve_KochanekBartelsSpline::GetCurrentFirstDerivative

  get the first derivative for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_KochanekBartelsSpline<type>::getCurrentFirstDerivative( const float time ) const {
	int i;
	float bvals[4], d, clampedTime;
	type v, t0, t1;

	if ( this->times.Num() == 1 ) {
		return ( this->values[0] - this->values[0] );
	}

	clampedTime = this->ClampedTime( time );
	i = this->indexForTime( clampedTime );
	TangentsForIndex( i - 1, t0, t1 );
	BasisFirstDerivative( i - 1, clampedTime, bvals );
	v = bvals[0] * this->valueForIndex( i - 1 );
	v += bvals[1] * this->valueForIndex( i );
	v += bvals[2] * t0;
	v += bvals[3] * t1;
	d = ( this->timeForIndex( i ) - this->timeForIndex( i-1 ) );
	return v / d;
}

/*
====================
CCurve_KochanekBartelsSpline::GetCurrentSecondDerivative

  get the second derivative for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_KochanekBartelsSpline<type>::getCurrentSecondDerivative( const float time ) const {
	int i;
	float bvals[4], d, clampedTime;
	type v, t0, t1;

	if ( this->times.Num() == 1 ) {
		return ( this->values[0] - this->values[0] );
	}

	clampedTime = this->ClampedTime( time );
	i = this->indexForTime( clampedTime );
	TangentsForIndex( i - 1, t0, t1 );
	BasisSecondDerivative( i - 1, clampedTime, bvals );
	v = bvals[0] * this->valueForIndex( i - 1 );
	v += bvals[1] * this->valueForIndex( i );
	v += bvals[2] * t0;
	v += bvals[3] * t1;
	d = ( this->timeForIndex( i ) - this->timeForIndex( i-1 ) );
	return v / ( d * d );
}

/*
====================
CCurve_KochanekBartelsSpline::TangentsForIndex
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_KochanekBartelsSpline<type>::TangentsForIndex( const int index, type &t0, type &t1 ) const {
	float dt, omt, omc, opc, omb, opb, adj, s0, s1;
	type delta;

	delta = this->valueForIndex( index + 1 ) - this->valueForIndex( index );
	dt = this->timeForIndex( index + 1 ) - this->timeForIndex( index );

	omt = 1.0f - tension[index];
	omc = 1.0f - continuity[index];
	opc = 1.0f + continuity[index];
	omb = 1.0f - bias[index];
	opb = 1.0f + bias[index];
	adj = 2.0f * dt / ( this->timeForIndex( index + 1 ) - this->timeForIndex( index - 1 ) );
	s0 = 0.5f * adj * omt * opc * opb;
	s1 = 0.5f * adj * omt * omc * omb;

	// outgoing tangent at first point
	t0 = s1 * delta + s0 * ( this->valueForIndex( index ) - this->valueForIndex( index - 1 ) );

	omt = 1.0f - tension[index + 1];
	omc = 1.0f - continuity[index + 1];
	opc = 1.0f + continuity[index + 1];
	omb = 1.0f - bias[index + 1];
	opb = 1.0f + bias[index + 1];
	adj = 2.0f * dt / ( this->timeForIndex( index + 2 ) - this->timeForIndex( index ) );
	s0 = 0.5f * adj * omt * omc * opb;
	s1 = 0.5f * adj * omt * opc * omb;

	// incoming tangent at second point
	t1 = s1 * ( this->valueForIndex( index + 2 ) - this->valueForIndex( index + 1 ) ) + s0 * delta;
}

/*
====================
CCurve_KochanekBartelsSpline::Basis

  spline basis functions
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_KochanekBartelsSpline<type>::Basis( const int index, const float t, float *bvals ) const {
	float s = (float) ( t - this->timeForIndex( index ) ) / ( this->timeForIndex( index+1 ) - this->timeForIndex( index ) );
	bvals[0] = ( ( 2.0f * s - 3.0f ) * s ) * s + 1.0f;				// 2.0f * s * s * s - 3.0f * s * s + 1.0f
	bvals[1] = ( ( -2.0f * s + 3.0f ) * s ) * s;					// -2.0f * s * s * s + 3.0f * s * s
	bvals[2] = ( ( s - 2.0f ) * s ) * s + s;						// s * s * s - 2.0f * s * s + s
	bvals[3] = ( ( s - 1.0f ) * s ) * s;							// s * s * s - s * s
}

/*
====================
CCurve_KochanekBartelsSpline::BasisFirstDerivative

  first derivative of spline basis functions
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_KochanekBartelsSpline<type>::BasisFirstDerivative( const int index, const float t, float *bvals ) const {
	float s = (float) ( t - this->timeForIndex( index ) ) / ( this->timeForIndex( index+1 ) - this->timeForIndex( index ) );
	bvals[0] = ( 6.0f * s - 6.0f ) * s;								// 6.0f * s * s - 6.0f * s
	bvals[1] = ( -6.0f * s + 6.0f ) * s;							// -6.0f * s * s + 6.0f * s
	bvals[2] = ( 3.0f * s - 4.0f ) * s + 1.0f;						// 3.0f * s * s - 4.0f * s + 1.0f
	bvals[3] = ( 3.0f * s - 2.0f ) * s;								// 3.0f * s * s - 2.0f * s
}

/*
====================
CCurve_KochanekBartelsSpline::BasisSecondDerivative

  second derivative of spline basis functions
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_KochanekBartelsSpline<type>::BasisSecondDerivative( const int index, const float t, float *bvals ) const {
	float s = (float) ( t - this->timeForIndex( index ) ) / ( this->timeForIndex( index+1 ) - this->timeForIndex( index ) );
	bvals[0] = 12.0f * s - 6.0f;
	bvals[1] = -12.0f * s + 6.0f;
	bvals[2] = 6.0f * s - 4.0f;
	bvals[3] = 6.0f * s - 2.0f;
}


/*
===============================================================================


===============================================================================
*/

/**
 * \class CCurve_BSpline
 *
 * \ingroup SMF_Math
 * 
 * \if pt_br
 * \brief	Template de curva B-Spline
 * \note    é muito lento e usa definições recursivas. Use CCurve_UniformCubicBSpline ou CCurve_NonUniformBSpline.
 * \elseif us_en
 * \brief 	B-Spline base template. Uses recursive definition and is slow.
 * \note	Use CCurve_UniformCubicBSpline or CCurve_NonUniformBSpline instead.
 *
 * \endif
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \see http://en.wikipedia.org/wiki/Spline_%28mathematics%29, http://pt.wikipedia.org/wiki/Spline
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
template< class type >
class CCurve_BSpline : public CCurve_Spline<type> {

public:
						CCurve_BSpline();

	virtual int			GetOrder() const { return order; }
	virtual void		SetOrder( const int i ) { SMF_ASSERT( i > 0 && i < 10 ); order = i; }

	virtual type		getCurrentValue( const float time ) const;
	virtual type		getCurrentFirstDerivative( const float time ) const;
	virtual type		getCurrentSecondDerivative( const float time ) const;

protected:
	int					order;

	float				Basis( const int index, const int order, const float t ) const;
	float				BasisFirstDerivative( const int index, const int order, const float t ) const;
	float				BasisSecondDerivative( const int index, const int order, const float t ) const;
};

/*
====================
CCurve_BSpline::CCurve_NaturalCubicSpline
====================
*/
template< class type >
SMF_INLINE_FORCED CCurve_BSpline<type>::CCurve_BSpline() {
	order = 4;	// default to cubic
}

/*
====================
CCurve_BSpline::GetCurrentValue

  get the value for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_BSpline<type>::getCurrentValue( const float time ) const {
	int i, j, k;
	float clampedTime;
	type v;

	if ( this->times.Num() == 1 ) {
		return this->values[0];
	}

	clampedTime = this->ClampedTime( time );
	i = this->indexForTime( clampedTime );
	v = this->values[0] - this->values[0];
	for ( j = 0; j < order; j++ ) {
		k = i + j - ( order >> 1 );
		v += Basis( k-2, order, clampedTime ) * this->valueForIndex( k );
	}
	return v;
}

/*
====================
CCurve_BSpline::GetCurrentFirstDerivative

  get the first derivative for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_BSpline<type>::getCurrentFirstDerivative( const float time ) const {
	int i, j, k;
	float clampedTime;
	type v;

	if ( this->times.Num() == 1 ) {
		return this->values[0];
	}

	clampedTime = this->ClampedTime( time );
	i = this->indexForTime( clampedTime );
	v = this->values[0] - this->values[0];
	for ( j = 0; j < order; j++ ) {
		k = i + j - ( order >> 1 );
		v += BasisFirstDerivative( k-2, order, clampedTime ) * this->valueForIndex( k );
	}
	return v;
}

/*
====================
CCurve_BSpline::GetCurrentSecondDerivative

  get the second derivative for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_BSpline<type>::getCurrentSecondDerivative( const float time ) const {
	int i, j, k;
	float clampedTime;
	type v;

	if ( this->times.Num() == 1 ) {
		return this->values[0];
	}

	clampedTime = this->ClampedTime( time );
	i = this->indexForTime( clampedTime );
	v = this->values[0] - this->values[0];
	for ( j = 0; j < order; j++ ) {
		k = i + j - ( order >> 1 );
		v += BasisSecondDerivative( k-2, order, clampedTime ) * this->valueForIndex( k );
	}
	return v;
}

/*
====================
CCurve_BSpline::Basis

  spline basis function
====================
*/
template< class type >
SMF_INLINE_FORCED float CCurve_BSpline<type>::Basis( const int index, const int order, const float t ) const {
	if ( order <= 1 ) {
		if ( this->timeForIndex( index ) < t && t <= this->timeForIndex( index + 1 ) ) {
			return 1.0f;
		} else {
			return 0.0f;
		}
	} else {
		float sum = 0.0f;
		float d1 = this->timeForIndex( index+order-1 ) - this->timeForIndex( index );
		if ( d1 != 0.0f ) {
			sum += (float) ( t - this->timeForIndex( index ) ) * Basis( index, order-1, t ) / d1;
		}

		float d2 = this->timeForIndex( index+order ) - this->timeForIndex( index+1 );
		if ( d2 != 0.0f ) {
			sum += (float) ( this->timeForIndex( index+order ) - t ) * Basis( index+1, order-1, t ) / d2;
		}
		return sum;
	}
}

/*
====================
CCurve_BSpline::BasisFirstDerivative

  first derivative of spline basis function
====================
*/
template< class type >
SMF_INLINE_FORCED float CCurve_BSpline<type>::BasisFirstDerivative( const int index, const int order, const float t ) const {
	return ( Basis( index, order-1, t ) - Basis( index+1, order-1, t ) ) *
			(float) ( order - 1 ) / ( this->timeForIndex( index + ( order - 1 ) - 2 ) - this->timeForIndex( index - 2 ) );
}

/*
====================
CCurve_BSpline::BasisSecondDerivative

  second derivative of spline basis function
====================
*/
template< class type >
SMF_INLINE_FORCED float CCurve_BSpline<type>::BasisSecondDerivative( const int index, const int order, const float t ) const {
	return ( BasisFirstDerivative( index, order-1, t ) - BasisFirstDerivative( index+1, order-1, t ) ) *
			(float) ( order - 1 ) / ( this->timeForIndex( index + ( order - 1 ) - 2 ) - this->timeForIndex( index - 2 ) );
}


/*
===============================================================================

	

===============================================================================
*/
/**
 * \class CCurve_UniformCubicBSpline
 *
 * \ingroup SMF_Math
 * 
 * \if pt_br
 * \brief	Template de curva B-Spline Cúbica Uniforme Não-racional
 * \note    é muito lento e usa definições recursivas. Use CCurve_UniformCubicBSpline ou CCurve_NonUniformBSpline.
 * \elseif us_en
 * \brief 	Uniform Non-Rational Cubic B-Spline template.
 * \endif
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \see http://en.wikipedia.org/wiki/Spline_%28mathematics%29, http://pt.wikipedia.org/wiki/Spline
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
template< class type >
class CCurve_UniformCubicBSpline : public CCurve_BSpline<type> {
	
public:
						CCurve_UniformCubicBSpline();

	virtual type		getCurrentValue( const float time ) const;
	virtual type		getCurrentFirstDerivative( const float time ) const;
	virtual type		getCurrentSecondDerivative( const float time ) const;

protected:
	void				Basis( const int index, const float t, float *bvals ) const;
	void				BasisFirstDerivative( const int index, const float t, float *bvals ) const;
	void				BasisSecondDerivative( const int index, const float t, float *bvals ) const;
};

/*
====================
CCurve_UniformCubicBSpline::CCurve_UniformCubicBSpline
====================
*/
template< class type >
SMF_INLINE_FORCED CCurve_UniformCubicBSpline<type>::CCurve_UniformCubicBSpline() {
	this->order = 4;	// always cubic
}

/*
====================
CCurve_UniformCubicBSpline::GetCurrentValue

  get the value for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_UniformCubicBSpline<type>::getCurrentValue( const float time ) const {
	int i, j, k;
	float bvals[4], clampedTime;
	type v;

	if ( this->times.Num() == 1 ) {
		return this->values[0];
	}

	clampedTime = this->ClampedTime( time );
	i = this->indexForTime( clampedTime );
	Basis( i-1, clampedTime, bvals );
	v = this->values[0] - this->values[0];
	for ( j = 0; j < 4; j++ ) {
		k = i + j - 2;
		v += bvals[j] * this->valueForIndex( k );
	}
	return v;
}

/*
====================
CCurve_UniformCubicBSpline::GetCurrentFirstDerivative

  get the first derivative for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_UniformCubicBSpline<type>::getCurrentFirstDerivative( const float time ) const {
	int i, j, k;
	float bvals[4], d, clampedTime;
	type v;

	if ( this->times.Num() == 1 ) {
		return ( this->values[0] - this->values[0] );
	}

	clampedTime = this->ClampedTime( time );
	i = this->indexForTime( clampedTime );
	BasisFirstDerivative( i-1, clampedTime, bvals );
	v = this->values[0] - this->values[0];
	for ( j = 0; j < 4; j++ ) {
		k = i + j - 2;
		v += bvals[j] * this->valueForIndex( k );
	}
	d = ( this->timeForIndex( i ) - this->timeForIndex( i-1 ) );
	return v / d;
}

/*
====================
CCurve_UniformCubicBSpline::GetCurrentSecondDerivative

  get the second derivative for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_UniformCubicBSpline<type>::getCurrentSecondDerivative( const float time ) const {
	int i, j, k;
	float bvals[4], d, clampedTime;
	type v;

	if ( this->times.Num() == 1 ) {
		return ( this->values[0] - this->values[0] );
	}

	clampedTime = this->ClampedTime( time );
	i = this->indexForTime( clampedTime );
	BasisSecondDerivative( i-1, clampedTime, bvals );
	v = this->values[0] - this->values[0];
	for ( j = 0; j < 4; j++ ) {
		k = i + j - 2;
		v += bvals[j] * this->valueForIndex( k );
	}
	d = ( this->timeForIndex( i ) - this->timeForIndex( i-1 ) );
	return v / ( d * d );
}

/*
====================
CCurve_UniformCubicBSpline::Basis

  spline basis functions
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_UniformCubicBSpline<type>::Basis( const int index, const float t, float *bvals ) const {
	float s = (float) ( t - this->timeForIndex( index ) ) / ( this->timeForIndex( index+1 ) - this->timeForIndex( index ) );
	bvals[0] = ( ( ( -s + 3.0f ) * s - 3.0f ) * s + 1.0f ) * ( 1.0f / 6.0f );
	bvals[1] = ( ( ( 3.0f * s - 6.0f ) * s ) * s + 4.0f ) * ( 1.0f / 6.0f );
	bvals[2] = ( ( ( -3.0f * s + 3.0f ) * s + 3.0f ) * s + 1.0f ) * ( 1.0f / 6.0f );
	bvals[3] = ( s * s * s ) * ( 1.0f / 6.0f );
}

/*
====================
CCurve_UniformCubicBSpline::BasisFirstDerivative

  first derivative of spline basis functions
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_UniformCubicBSpline<type>::BasisFirstDerivative( const int index, const float t, float *bvals ) const {
	float s = (float) ( t - this->timeForIndex( index ) ) / ( this->timeForIndex( index+1 ) - this->timeForIndex( index ) );
	bvals[0] = -0.5f * s * s + s - 0.5f;
	bvals[1] = 1.5f * s * s - 2.0f * s;
	bvals[2] = -1.5f * s * s + s + 0.5f;
	bvals[3] = 0.5f * s * s;
}

/*
====================
CCurve_UniformCubicBSpline::BasisSecondDerivative

  second derivative of spline basis functions
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_UniformCubicBSpline<type>::BasisSecondDerivative( const int index, const float t, float *bvals ) const {
	float s = (float) ( t - this->timeForIndex( index ) ) / ( this->timeForIndex( index+1 ) - this->timeForIndex( index ) );
	bvals[0] = -s + 1.0f;
	bvals[1] = 3.0f * s - 2.0f;
	bvals[2] = -3.0f * s + 1.0f;
	bvals[3] = s;
}



/**
 * \class CCurve_NonUniformBSpline
 *
 * \ingroup SMF_Math
 * 
 * \if pt_br
 * \brief	Template de curva B-Spline Cúbica Não-Uniforme Não-racional
 * \note    é muito lento e usa definições recursivas. Use CCurve_UniformCubicBSpline ou CCurve_NonUniformBSpline.
 * \elseif us_en
 * \brief 	Non-Uniform Non-Rational B-Spline (NUBS) template.
 * \endif
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \see http://en.wikipedia.org/wiki/Spline_%28mathematics%29, http://pt.wikipedia.org/wiki/Spline
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
template< class type >
class CCurve_NonUniformBSpline : public CCurve_BSpline<type> {
	
public:
						CCurve_NonUniformBSpline();

	virtual type		getCurrentValue( const float time ) const;
	virtual type		getCurrentFirstDerivative( const float time ) const;
	virtual type		getCurrentSecondDerivative( const float time ) const;

protected:
	void				Basis( const int index, const int order, const float t, float *bvals ) const;
	void				BasisFirstDerivative( const int index, const int order, const float t, float *bvals ) const;
	void				BasisSecondDerivative( const int index, const int order, const float t, float *bvals ) const;
};

/*
====================
CCurve_NonUniformBSpline::CCurve_NonUniformBSpline
====================
*/
template< class type >
SMF_INLINE_FORCED CCurve_NonUniformBSpline<type>::CCurve_NonUniformBSpline() {
}

/*
====================
CCurve_NonUniformBSpline::GetCurrentValue

  get the value for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_NonUniformBSpline<type>::getCurrentValue( const float time ) const {
	int i, j, k;
	float clampedTime;
	type v;
	float *bvals = (float *) _alloca16( this->order * sizeof(float) );

	if ( this->times.Num() == 1 ) {
		return this->values[0];
	}

	clampedTime = this->ClampedTime( time );
	i = this->indexForTime( clampedTime );
	Basis( i-1, this->order, clampedTime, bvals );
	v = this->values[0] - this->values[0];
	for ( j = 0; j < this->order; j++ ) {
		k = i + j - ( this->order >> 1 );
		v += bvals[j] * this->valueForIndex( k );
	}
	return v;
}

/*
====================
CCurve_NonUniformBSpline::GetCurrentFirstDerivative

  get the first derivative for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_NonUniformBSpline<type>::getCurrentFirstDerivative( const float time ) const {
	int i, j, k;
	float clampedTime;
	type v;
	float *bvals = (float *) _alloca16( this->order * sizeof(float) );

	if ( this->times.Num() == 1 ) {
		return ( this->values[0] - this->values[0] );
	}

	clampedTime = this->ClampedTime( time );
	i = this->indexForTime( clampedTime );
	BasisFirstDerivative( i-1, this->order, clampedTime, bvals );
	v = this->values[0] - this->values[0];
	for ( j = 0; j < this->order; j++ ) {
		k = i + j - ( this->order >> 1 );
		v += bvals[j] * this->valueForIndex( k );
	}
	return v;
}

/*
====================
CCurve_NonUniformBSpline::GetCurrentSecondDerivative

  get the second derivative for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_NonUniformBSpline<type>::getCurrentSecondDerivative( const float time ) const {
	int i, j, k;
	float clampedTime;
	type v;
	float *bvals = (float *) _alloca16( this->order * sizeof(float) );

	if ( this->times.Num() == 1 ) {
		return ( this->values[0] - this->values[0] );
	}

	clampedTime = this->ClampedTime( time );
	i = this->indexForTime( clampedTime );
	BasisSecondDerivative( i-1, this->order, clampedTime, bvals );
	v = this->values[0] - this->values[0];
	for ( j = 0; j < this->order; j++ ) {
		k = i + j - ( this->order >> 1 );
		v += bvals[j] * this->valueForIndex( k );
	}
	return v;
}

/*
====================
CCurve_NonUniformBSpline::Basis

  spline basis functions
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_NonUniformBSpline<type>::Basis( const int index, const int order, const float t, float *bvals ) const {
    int r, s, i;
    float omega;

    bvals[order-1] = 1.0f;
    for ( r = 2; r <= order; r++ ) {
		i = index - r + 1;
		bvals[order - r] = 0.0f;
		for ( s = order - r + 1; s < order; s++ ) {
			i++;
			omega = (float) ( t - this->timeForIndex( i ) ) / ( this->timeForIndex( i + r - 1 ) - this->timeForIndex( i ) );
			bvals[s - 1] += ( 1.0f - omega ) * bvals[s];
			bvals[s] *= omega;
		}
    }
}

/*
====================
CCurve_NonUniformBSpline::BasisFirstDerivative

  first derivative of spline basis functions
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_NonUniformBSpline<type>::BasisFirstDerivative( const int index, const int order, const float t, float *bvals ) const {
	int i;

	Basis( index, order-1, t, bvals+1 );
	bvals[0] = 0.0f;
	for ( i = 0; i < order-1; i++ ) {
		bvals[i] -= bvals[i+1];
		bvals[i] *= (float) ( order - 1) / ( this->timeForIndex( index + i + (order-1) - 2 ) - this->timeForIndex( index + i - 2 ) );
	}
	bvals[i] *= (float) ( order - 1) / ( this->timeForIndex( index + i + (order-1) - 2 ) - this->timeForIndex( index + i - 2 ) );
}

/*
====================
CCurve_NonUniformBSpline::BasisSecondDerivative

  second derivative of spline basis functions
====================
*/
template< class type >
SMF_INLINE_FORCED void CCurve_NonUniformBSpline<type>::BasisSecondDerivative( const int index, const int order, const float t, float *bvals ) const {
	int i;

	BasisFirstDerivative( index, order-1, t, bvals+1 );
	bvals[0] = 0.0f;
	for ( i = 0; i < order-1; i++ ) {
		bvals[i] -= bvals[i+1];
		bvals[i] *= (float) ( order - 1) / ( this->timeForIndex( index + i + (order-1) - 2 ) - this->timeForIndex( index + i - 2 ) );
	}
	bvals[i] *= (float) ( order - 1) / ( this->timeForIndex( index + i + (order-1) - 2 ) - this->timeForIndex( index + i - 2 ) );
}



/**
 * \class CCurve_NURBS
 *
 * \ingroup SMF_Math
 * 
 * \if pt_br
 * \brief	Template de curva B-Spline  Não-Uniforme racional (NURBS)
 * \note    é muito lento e usa definições recursivas. Use CCurve_UniformCubicBSpline ou CCurve_NonUniformBSpline.
 * \elseif us_en
 * \brief 	Non-Uniform Rational B-Spline (NURBS) template.
 * \endif
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \see http://en.wikipedia.org/wiki/Spline_%28mathematics%29, http://pt.wikipedia.org/wiki/Spline
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
template< class type >
class CCurve_NURBS : public CCurve_NonUniformBSpline<type> {
	
public:
						CCurve_NURBS();

	virtual int			addValue( const float time, const type &value );
	virtual int			addValue( const float time, const type &value, const float weight );
	virtual void		removeIndex( const int index ) { this->values.removeIndex(index); this->times.removeIndex(index); weights.removeIndex(index); }
	virtual void		clear() { this->values.clear(); this->times.clear(); weights.clear(); this->currentIndex = -1; }

	virtual type		getCurrentValue( const float time ) const;
	virtual type		getCurrentFirstDerivative( const float time ) const;
	virtual type		getCurrentSecondDerivative( const float time ) const;

protected:
	CList<float>		weights;

	float				WeightForIndex( const int index ) const;
};

/*
====================
CCurve_NURBS::CCurve_NURBS
====================
*/
template< class type >
SMF_INLINE_FORCED CCurve_NURBS<type>::CCurve_NURBS() {
}

/*
====================
CCurve_NURBS::addValue

  add a timed/value pair to the spline
  returns the index to the inserted pair
====================
*/
template< class type >
SMF_INLINE_FORCED int CCurve_NURBS<type>::addValue( const float time, const type &value ) {
	int i;

	i = this->indexForTime( time );
	this->times.insert( time, i );
	this->values.insert( value, i );
	weights.insert( 1.0f, i );
	return i;
}

/*
====================
CCurve_NURBS::addValue

  add a timed/value pair to the spline
  returns the index to the inserted pair
====================
*/
template< class type >
SMF_INLINE_FORCED int CCurve_NURBS<type>::addValue( const float time, const type &value, const float weight ) {
	int i;

	i = this->indexForTime( time );
	this->times.insert( time, i );
	this->values.insert( value, i );
	weights.insert( weight, i );
	return i;
}

/*
====================
CCurve_NURBS::GetCurrentValue

  get the value for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_NURBS<type>::getCurrentValue( const float time ) const {
	int i, j, k;
	float w, b, *bvals, clampedTime;
	type v;

	if ( this->times.Num() == 1 ) {
		return this->values[0];
	}

	bvals = (float *) _alloca16( this->order * sizeof(float) );

	clampedTime = this->ClampedTime( time );
	i = this->indexForTime( clampedTime );
	this->Basis( i-1, this->order, clampedTime, bvals );
	v = this->values[0] - this->values[0];
	w = 0.0f;
	for ( j = 0; j < this->order; j++ ) {
		k = i + j - ( this->order >> 1 );
		b = bvals[j] * WeightForIndex( k );
		w += b;
		v += b * this->valueForIndex( k );
	}
	return v / w;
}

/*
====================
CCurve_NURBS::GetCurrentFirstDerivative

  get the first derivative for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_NURBS<type>::getCurrentFirstDerivative( const float time ) const {
	int i, j, k;
	float w, wb, wd1, b, d1, *bvals, *d1vals, clampedTime;
	type v, vb, vd1;

	if ( this->times.Num() == 1 ) {
		return this->values[0];
	}

	bvals = (float *) _alloca16( this->order * sizeof(float) );
	d1vals = (float *) _alloca16( this->order * sizeof(float) );

	clampedTime = this->ClampedTime( time );
	i = this->indexForTime( clampedTime );
	this->Basis( i-1, this->order, clampedTime, bvals );
	this->BasisFirstDerivative( i-1, this->order, clampedTime, d1vals );
	vb = vd1 = this->values[0] - this->values[0];
	wb = wd1 = 0.0f;
	for ( j = 0; j < this->order; j++ ) {
		k = i + j - ( this->order >> 1 );
		w = WeightForIndex( k );
		b = bvals[j] * w;
		d1 = d1vals[j] * w;
		wb += b;
		wd1 += d1;
		v = this->valueForIndex( k );
		vb += b * v;
		vd1 += d1 * v;
	}
	return ( wb * vd1 - vb * wd1 ) / ( wb * wb );
}

/*
====================
CCurve_NURBS::GetCurrentSecondDerivative

  get the second derivative for the given time
====================
*/
template< class type >
SMF_INLINE_FORCED type CCurve_NURBS<type>::getCurrentSecondDerivative( const float time ) const {
	int i, j, k;
	float w, wb, wd1, wd2, b, d1, d2, *bvals, *d1vals, *d2vals, clampedTime;
	type v, vb, vd1, vd2;

	if ( this->times.Num() == 1 ) {
		return this->values[0];
	}

	bvals = (float *) _alloca16( this->order * sizeof(float) );
	d1vals = (float *) _alloca16( this->order * sizeof(float) );
	d2vals = (float *) _alloca16( this->order * sizeof(float) );

	clampedTime = this->ClampedTime( time );
	i = this->indexForTime( clampedTime );
	this->Basis( i-1, this->order, clampedTime, bvals );
	this->BasisFirstDerivative( i-1, this->order, clampedTime, d1vals );
	this->BasisSecondDerivative( i-1, this->order, clampedTime, d2vals );
	vb = vd1 = vd2 = this->values[0] - this->values[0];
	wb = wd1 = wd2 = 0.0f;
	for ( j = 0; j < this->order; j++ ) {
		k = i + j - ( this->order >> 1 );
		w = WeightForIndex( k );
		b = bvals[j] * w;
		d1 = d1vals[j] * w;
		d2 = d2vals[j] * w;
		wb += b;
		wd1 += d1;
		wd2 += d2;
		v = this->valueForIndex( k );
		vb += b * v;
		vd1 += d1 * v;
		vd2 += d2 * v;
	}
	return ( ( wb * wb ) * ( wb * vd2 - vb * wd2 ) - ( wb * vd1 - vb * wd1 ) * 2.0f * wb * wd1 ) / ( wb * wb * wb * wb );
}

/*
====================
CCurve_NURBS::WeightForIndex

  get the weight for the given index
====================
*/
template< class type >
SMF_INLINE_FORCED float CCurve_NURBS<type>::WeightForIndex( const int index ) const {
	int n = weights.getNum()-1;

	if ( index < 0 ) {
		if ( this->boundaryType == CCurve_Spline<type>::BT_CLOSED ) {
			return weights[ weights.getNum() + index % weights.getNum() ];
		} else {
			return weights[0] + index * ( weights[1] - weights[0] );
		}
	} else if ( index > n ) {
		if ( this->boundaryType == CCurve_Spline<type>::BT_CLOSED ) {
			return weights[ index % weights.getNum() ];
		} else {
			return weights[n] + ( index - n ) * ( weights[n] - weights[n-1] );
		}
	}
	return weights[index];
}
} // end MATH
} //end SMF
#endif /* !__MATH_CURVE_H__ */
