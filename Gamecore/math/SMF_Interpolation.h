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

#ifndef _SMF__MATH_INTERPOLATE_H__
#define _SMF__MATH_INTERPOLATE_H__
#include "../SMF_Config.h"
namespace SMF {
namespace MATH{


/**
 * \class CInterpolation
 *
 * \ingroup SMF_Math
 *
 * \brief Implementa interpolação linear, com fórmula de Lagrange
 * 	Linear interpolation.
	Interpolação de Lagrange
	\see http://pt.wikipedia.org/wiki/Interpolação

 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *

 \note http://pt.wikipedia.org/wiki/Interpolação
 \note ya = y1 + ( ((xa-x1)/(x2-x1))*(y2-y1) ) sendo:
       ya = valor que se deseja conhecer (interpolar)
	   xa = parâmetro de tempo fornecido
	   x1 = tempo inicial
	   x2 = tempo final
	   y1 = valor inicial da interpolação
	   y2 = valor final da interpolação
 */
template< class type >
class  SMF_API CInterpolation {
public:
	/**
	\brief construtor: cria um objeto interpolação linear
	\note o eixo x é o eixo do tempo
	**/
	CInterpolation();
	/**
	\brief inicializa os valores conhecidos da interpolação
	\param startTime valor x1 da interpolação
	\param duration valor x1+ tempo de duração
	\param startValue valor y1
	\param endValue valor y2
	\note currentValue será proporcional (Interpolação de Lagrange)
	**/
	void				init( const float startTime, const float duration, const type &startValue, const type &endValue );
	void				setStartTime( float time ) { this->startTime = time; }
	void				setDuration( float duration ) { this->duration = duration; }
	void				setStartValue( const type &startValue ) { this->startValue = startValue; }
	void				setEndValue( const type &endValue ) { this->endValue = endValue; }
	/**
	\brief retorna o valor da interpolação, no parÂmetro passado
	\param time tempo no qual se deseja saber o valor da interpolação
	\return retorna o valor da interpolação no tempo dado
	**/
	type				getCurrentValue( float time ) const;
	/**
	\brief verifica se o tempo passado está dentro do tempo da interpolação
	\return verdadeiro se o tempo passado está dentro do intervalo da interpolaçãp, falso em caso contrário
	**/
	bool				isDone( float time ) const { return ( time >= startTime + duration ); }

	float				getStartTime() const { return startTime; }
	float				getEndTime() const { return startTime + duration; }
	float				getDuration() const { return duration; }
	const type &		getStartValue() const { return startValue; }
	const type &		getEndValue() const { return endValue; }

private:
	float				startTime;
	float				duration;
	type				startValue;
	type				endValue;
	mutable float		currentTime;
	mutable type		currentValue;
};

/*
====================
CInterpolation::CInterpolation
====================
*/
template< class type >
SMF_INLINE_FORCED CInterpolation<type>::CInterpolation() {
	currentTime = startTime = duration = 0;
	memset( &currentValue, 0, sizeof( currentValue ) );
	startValue = endValue = currentValue;
}

/*
====================
CInterpolation::init
====================
*/
template< class type >
SMF_INLINE_FORCED void CInterpolation<type>::init( const float startTime, const float duration, const type &startValue, const type &endValue ) {
	this->startTime = startTime;
	this->duration = duration;
	this->startValue = startValue;
	this->endValue = endValue;
	this->currentTime = startTime - 1;
	this->currentValue = startValue;
}

/*
====================
CInterpolation::GetCurrentValue
====================
*/
template< class type >
SMF_INLINE_FORCED type CInterpolation<type>::getCurrentValue( float time ) const {
	float deltaTime;

	deltaTime = time - startTime;
	if ( time != currentTime ) {
		currentTime = time;
		if ( deltaTime <= 0 ) {
			currentValue = startValue;
		} else if ( deltaTime >= duration ) {
			currentValue = endValue;
		} else {
			currentValue = startValue + ( endValue - startValue ) * ( (float) deltaTime / duration );
		}
	}
	return currentValue;
}


/*
==============================================================================================

	Continuous interpolation with linear acceleration and deceleration phase.
	The velocity is continuous but the acceleration is not.

==============================================================================================
*/

template< class type >
class CInterpolateAccelDecelLinear  {
public:
						CInterpolateAccelDecelLinear();

	void				init( const float startTime, const float accelTime, const float decelTime, const float duration, const type &startValue, const type &endValue );
	void				SetStartTime( float time ) { startTime = time; Invalidate(); }
	void				SetStartValue( const type &startValue ) { this->startValue = startValue; Invalidate(); }
	void				SetEndValue( const type &endValue ) { this->endValue = endValue; Invalidate(); }

	type				GetCurrentValue( float time ) const;
	type				GetCurrentSpeed( float time ) const;
	bool				IsDone( float time ) const { return ( time >= startTime + accelTime + linearTime + decelTime ); }

	float				GetStartTime() const { return startTime; }
	float				GetEndTime() const { return startTime + accelTime + linearTime + decelTime; }
	float				GetDuration() const { return accelTime + linearTime + decelTime; }
	float				GetAcceleration() const { return accelTime; }
	float				GetDeceleration() const { return decelTime; }
	const type &		GetStartValue() const { return startValue; }
	const type &		GetEndValue() const { return endValue; }

private:
	float				startTime;
	float				accelTime;
	float				linearTime;
	float				decelTime;
	type				startValue;
	type				endValue;
	mutable CExtrapolation<type> extrapolate;

	void				Invalidate();
	void				SetPhase( float time ) const;
};

/*
====================
CInterpolateAccelDecelLinear::CInterpolateAccelDecelLinear
====================
*/
template< class type >
SMF_INLINE_FORCED CInterpolateAccelDecelLinear<type>::CInterpolateAccelDecelLinear() {
	startTime = accelTime = linearTime = decelTime = 0;
	memset( &startValue, 0, sizeof( startValue ) );
	endValue = startValue;
}

/*
====================
CInterpolateAccelDecelLinear::init
====================
*/
template< class type >
SMF_INLINE_FORCED void CInterpolateAccelDecelLinear<type>::init( const float startTime, const float accelTime, const float decelTime, const float duration, const type &startValue, const type &endValue ) {
	type speed;

	this->startTime = startTime;
	this->accelTime = accelTime;
	this->decelTime = decelTime;
	this->startValue = startValue;
	this->endValue = endValue;

	if ( duration <= 0.0f ) {
		return;
	}

	if ( this->accelTime + this->decelTime > duration ) {
		this->accelTime = this->accelTime * duration / ( this->accelTime + this->decelTime );
		this->decelTime = duration - this->accelTime;
	}
	this->linearTime = duration - this->accelTime - this->decelTime;
	speed = ( endValue - startValue ) * ( 1000.0f / ( (float) this->linearTime + ( this->accelTime + this->decelTime ) * 0.5f ) );

	if ( this->accelTime ) {
		extrapolate.init( startTime, this->accelTime, startValue, ( startValue - startValue ), speed, EXTRAP_ACCELLINEAR );
	} else if ( this->linearTime ) {
		extrapolate.init( startTime, this->linearTime, startValue, ( startValue - startValue ), speed, EXTRAP_LINEAR );
	} else {
		extrapolate.init( startTime, this->decelTime, startValue, ( startValue - startValue ), speed, EXTRAP_DECELLINEAR );
	}
}

/*
====================
CInterpolateAccelDecelLinear::Invalidate
====================
*/
template< class type >
SMF_INLINE_FORCED void CInterpolateAccelDecelLinear<type>::Invalidate() {
	extrapolate.init( 0, 0, extrapolate.GetStartValue(), extrapolate.GetBaseSpeed(), extrapolate.getSpeed(), EXTRAP_NONE );
}

/*
====================
CInterpolateAccelDecelLinear::SetPhase
====================
*/
template< class type >
SMF_INLINE_FORCED void CInterpolateAccelDecelLinear<type>::SetPhase( float time ) const {
	float deltaTime;

	deltaTime = time - startTime;
	if ( deltaTime < accelTime ) {
		if ( extrapolate.GetExtrapolationType() != EXTRAP_ACCELLINEAR ) {
			extrapolate.init( startTime, accelTime, startValue, extrapolate.GetBaseSpeed(), extrapolate.getSpeed(), EXTRAP_ACCELLINEAR );
		}
	} else if ( deltaTime < accelTime + linearTime ) {
		if ( extrapolate.GetExtrapolationType() != EXTRAP_LINEAR ) {
			extrapolate.init( startTime + accelTime, linearTime, startValue + extrapolate.getSpeed() * ( accelTime * 0.001f * 0.5f ), extrapolate.GetBaseSpeed(), extrapolate.getSpeed(), EXTRAP_LINEAR );
		}
	} else {
		if ( extrapolate.GetExtrapolationType() != EXTRAP_DECELLINEAR ) {
			extrapolate.init( startTime + accelTime + linearTime, decelTime, endValue - ( extrapolate.getSpeed() * ( decelTime * 0.001f * 0.5f ) ), extrapolate.GetBaseSpeed(), extrapolate.getSpeed(), EXTRAP_DECELLINEAR );
		}
	}
}

/*
====================
CInterpolateAccelDecelLinear::GetCurrentValue
====================
*/
template< class type >
SMF_INLINE_FORCED type CInterpolateAccelDecelLinear<type>::GetCurrentValue( float time ) const {
	SetPhase( time );
	return extrapolate.GetCurrentValue( time );
}

/*
====================
CInterpolateAccelDecelLinear::GetCurrentSpeed
====================
*/
template< class type >
SMF_INLINE_FORCED type CInterpolateAccelDecelLinear<type>::GetCurrentSpeed( float time ) const {
	SetPhase( time );
	return extrapolate.GetCurrentSpeed( time );
}


/*
==============================================================================================

	Continuous interpolation with sinusoidal acceleration and deceleration phase.
	Both the velocity and acceleration are continuous.

==============================================================================================
*/

template< class type >
class CInterpolateAccelDecelSine  {
public:
						CInterpolateAccelDecelSine();

	void				init( const float startTime, const float accelTime, const float decelTime, const float duration, const type &startValue, const type &endValue );
	void				SetStartTime( float time ) { startTime = time; Invalidate(); }
	void				SetStartValue( const type &startValue ) { this->startValue = startValue; Invalidate(); }
	void				SetEndValue( const type &endValue ) { this->endValue = endValue; Invalidate(); }

	type				GetCurrentValue( float time ) const;
	type				GetCurrentSpeed( float time ) const;
	bool				IsDone( float time ) const { return ( time >= startTime + accelTime + linearTime + decelTime ); }

	float				GetStartTime() const { return startTime; }
	float				GetEndTime() const { return startTime + accelTime + linearTime + decelTime; }
	float				GetDuration() const { return accelTime + linearTime + decelTime; }
	float				GetAcceleration() const { return accelTime; }
	float				GetDeceleration() const { return decelTime; }
	const type &		GetStartValue() const { return startValue; }
	const type &		GetEndValue() const { return endValue; }

private:
	float				startTime;
	float				accelTime;
	float				linearTime;
	float				decelTime;
	type				startValue;
	type				endValue;
	mutable CExtrapolation<type> extrapolate;

	void				Invalidate();
	void				SetPhase( float time ) const;
};

/*
====================
CInterpolateAccelDecelSine::CInterpolateAccelDecelSine
====================
*/
template< class type >
SMF_INLINE_FORCED CInterpolateAccelDecelSine<type>::CInterpolateAccelDecelSine() {
	startTime = accelTime = linearTime = decelTime = 0;
	memset( &startValue, 0, sizeof( startValue ) );
	endValue = startValue;
}

/*
====================
CInterpolateAccelDecelSine::init
====================
*/
template< class type >
SMF_INLINE_FORCED void CInterpolateAccelDecelSine<type>::init( const float startTime, const float accelTime, const float decelTime, const float duration, const type &startValue, const type &endValue ) {
	type speed;

	this->startTime = startTime;
	this->accelTime = accelTime;
	this->decelTime = decelTime;
	this->startValue = startValue;
	this->endValue = endValue;

	if ( duration <= 0.0f ) {
		return;
	}

	if ( this->accelTime + this->decelTime > duration ) {
		this->accelTime = this->accelTime * duration / ( this->accelTime + this->decelTime );
		this->decelTime = duration - this->accelTime;
	}
	this->linearTime = duration - this->accelTime - this->decelTime;
	speed = ( endValue - startValue ) * ( 1000.0f / ( (float) this->linearTime + ( this->accelTime + this->decelTime ) * CMath::SQRT_1OVER2 ) );

	if ( this->accelTime ) {
		extrapolate.init( startTime, this->accelTime, startValue, ( startValue - startValue ), speed, EXTRAP_ACCELSINE );
	} else if ( this->linearTime ) {
		extrapolate.init( startTime, this->linearTime, startValue, ( startValue - startValue ), speed, EXTRAP_LINEAR );
	} else {
		extrapolate.init( startTime, this->decelTime, startValue, ( startValue - startValue ), speed, EXTRAP_DECELSINE );
	}
}

/*
====================
CInterpolateAccelDecelSine::Invalidate
====================
*/
template< class type >
SMF_INLINE_FORCED void CInterpolateAccelDecelSine<type>::Invalidate() {
	extrapolate.init( 0, 0, extrapolate.GetStartValue(), extrapolate.GetBaseSpeed(), extrapolate.getSpeed(), EXTRAP_NONE );
}

/*
====================
CInterpolateAccelDecelSine::SetPhase
====================
*/
template< class type >
SMF_INLINE_FORCED void CInterpolateAccelDecelSine<type>::SetPhase( float time ) const {
	float deltaTime;

	deltaTime = time - startTime;
	if ( deltaTime < accelTime ) {
		if ( extrapolate.GetExtrapolationType() != EXTRAP_ACCELSINE ) {
			extrapolate.init( startTime, accelTime, startValue, extrapolate.GetBaseSpeed(), extrapolate.getSpeed(), EXTRAP_ACCELSINE );
		}
	} else if ( deltaTime < accelTime + linearTime ) {
		if ( extrapolate.GetExtrapolationType() != EXTRAP_LINEAR ) {
			extrapolate.init( startTime + accelTime, linearTime, startValue + extrapolate.getSpeed() * ( accelTime * 0.001f * CMath::SQRT_1OVER2 ), extrapolate.GetBaseSpeed(), extrapolate.getSpeed(), EXTRAP_LINEAR );
		}
	} else {
		if ( extrapolate.GetExtrapolationType() != EXTRAP_DECELSINE ) {
			extrapolate.init( startTime + accelTime + linearTime, decelTime, endValue - ( extrapolate.getSpeed() * ( decelTime * 0.001f * CMath::SQRT_1OVER2 ) ), extrapolate.GetBaseSpeed(), extrapolate.getSpeed(), EXTRAP_DECELSINE );
		}
	}
}

/*
====================
CInterpolateAccelDecelSine::GetCurrentValue
====================
*/
template< class type >
SMF_INLINE_FORCED type CInterpolateAccelDecelSine<type>::GetCurrentValue( float time ) const {
	SetPhase( time );
	return extrapolate.GetCurrentValue( time );
}

/*
====================
CInterpolateAccelDecelSine::GetCurrentSpeed
====================
*/
template< class type >
SMF_INLINE_FORCED type CInterpolateAccelDecelSine<type>::GetCurrentSpeed( float time ) const {
	SetPhase( time );
	return extrapolate.GetCurrentSpeed( time );
}
} //end MATH
} //end SMF
#endif /* !__MATH_INTERPOLATE_H__ */
