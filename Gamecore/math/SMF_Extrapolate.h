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

#ifndef _SMF__MATH_EXTRAPOLATE_H__
#define _SMF__MATH_EXTRAPOLATE_H__
#include "../SMF_Config.h"
namespace SMF {
namespace MATH{


/// type of extrapolation
typedef enum {
	/// no extrapolation, covered distance = duration * 0.001 * ( baseSpeed )
	EXTRAP_NONE			= 0x01,
	/// linear extrapolation, covered distance = duration * 0.001 * ( baseSpeed + speed )
	EXTRAP_LINEAR		= 0x02,
	/// linear acceleration, covered distance = duration * 0.001 * ( baseSpeed + 0.5 * speed )
	EXTRAP_ACCELLINEAR	= 0x04,
	/// linear deceleration, covered distance = duration * 0.001 * ( baseSpeed + 0.5 * speed )
	EXTRAP_DECELLINEAR	= 0x08,
	/// sinusoidal acceleration, covered distance = duration * 0.001 * ( baseSpeed + sqrt( 0.5 ) * speed )
	EXTRAP_ACCELSINE		= 0x10,
	/// sinusoidal deceleration, covered distance = duration * 0.001 * ( baseSpeed + sqrt( 0.5 ) * speed )
	EXTRAP_DECELSINE		= 0x20,
	/// do not stop at startTime + duration
	EXTRAP_NOSTOP		= 0x40
} extrapolation_t;
/**
 * \class CExtrapolation
 *
 * \ingroup SMF_Math
 *
 * \brief Implementa extrapolação linear
 * \see http://en.wikipedia.org/wiki/Extrapolation
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *

 \note http://en.wikipedia.org/wiki/Extrapolation
 **/
template< class type >
class  SMF_API CExtrapolation {
public:

	/**
	\brief construtor: cria um objeto extrapolation
	**/
	CExtrapolation();


	void				init( const float startTime, const float duration, const type &startValue, const type &baseSpeed, const type &speed, const extrapolation_t extrapolationType );
	type				getCurrentValue( float time ) const;
	type				getCurrentSpeed( float time ) const;
	bool				isDone( float time ) const { return ( !( extrapolationType & EXTRAP_NOSTOP ) && time >= startTime + duration ); }
	void				setStartTime( float time ) { startTime = time; currentTime = -1; }
	float				getStartTime() const { return startTime; }
	float				getEndTime() const { return ( !( extrapolationType & EXTRAP_NOSTOP ) && duration > 0 ) ? startTime + duration : 0; }
	float				getDuration() const { return duration; }
	void				setStartValue( const type &value ) { startValue = value; currentTime = -1; }
	const type &		getStartValue() const { return startValue; }
	const type &		getBaseSpeed() const { return baseSpeed; }
	const type &		getSpeed() const { return speed; }
	extrapolation_t		getExtrapolationType() const { return extrapolationType; }

private:
	extrapolation_t		extrapolationType;
	float				startTime;
	float				duration;
	type				startValue;
	type				baseSpeed;
	type				speed;
	mutable float		currentTime;
	mutable type		currentValue;
};

/*
====================
CExtrapolation::CExtrapolation
====================
*/
template< class type >
SMF_INLINE_FORCED CExtrapolation<type>::CExtrapolation() {
	extrapolationType = EXTRAP_NONE;
	startTime = duration = 0.0f;
	memset( &startValue, 0, sizeof( startValue ) );
	memset( &baseSpeed, 0, sizeof( baseSpeed ) );
	memset( &speed, 0, sizeof( speed ) );
	currentTime = -1;
	currentValue = startValue;
}

/*
====================
CExtrapolation::init
====================
*/
template< class type >
SMF_INLINE_FORCED void CExtrapolation<type>::init( const float startTime, const float duration, const type &startValue, const type &baseSpeed, const type &speed, const extrapolation_t extrapolationType ) {
	this->extrapolationType = extrapolationType;
	this->startTime = startTime;
	this->duration = duration;
	this->startValue = startValue;
	this->baseSpeed = baseSpeed;
	this->speed = speed;
	currentTime = -1;
	currentValue = startValue;
}

/*
====================
CExtrapolation::GetCurrentValue
====================
*/
template< class type >
SMF_INLINE_FORCED type CExtrapolation<type>::getCurrentValue( float time ) const {
	float deltaTime, s;

	if ( time == currentTime ) {
		return currentValue;
	}

	currentTime = time;

	if ( time < startTime ) {
		return startValue;
	}

	if ( !( extrapolationType &	EXTRAP_NOSTOP ) && ( time > startTime + duration ) ) {
		time = startTime + duration;
	}

	switch( extrapolationType & ~EXTRAP_NOSTOP ) {
		case EXTRAP_NONE: {
			deltaTime = ( time - startTime ) * 0.001f;
			currentValue = startValue + deltaTime * baseSpeed;
			break;
		}
		case EXTRAP_LINEAR: {
			deltaTime = ( time - startTime ) * 0.001f;
			currentValue = startValue + deltaTime * ( baseSpeed + speed );
			break;
		}
		case EXTRAP_ACCELLINEAR: {
			if ( !duration ) {
				currentValue = startValue;
			} else {
				deltaTime = ( time - startTime ) / duration;
				s = ( 0.5f * deltaTime * deltaTime ) * ( duration * 0.001f );
				currentValue = startValue + deltaTime * baseSpeed + s * speed;
			}
			break;
		}
		case EXTRAP_DECELLINEAR: {
			if ( !duration ) {
				currentValue = startValue;
			} else {
				deltaTime = ( time - startTime ) / duration;
				s = ( deltaTime - ( 0.5f * deltaTime * deltaTime ) ) * ( duration * 0.001f );
				currentValue = startValue + deltaTime * baseSpeed + s * speed;
			}
			break;
		}
		case EXTRAP_ACCELSINE: {
			if ( !duration ) {
				currentValue = startValue;
			} else {
				deltaTime = ( time - startTime ) / duration;
				s = ( 1.0f - CMath::cos( deltaTime * CMath::HALF_PI ) ) * duration * 0.001f * CMath::SQRT_1OVER2;
				currentValue = startValue + deltaTime * baseSpeed + s * speed;
			}
			break;
		}
		case EXTRAP_DECELSINE: {
			if ( !duration ) {
				currentValue = startValue;
			} else {
				deltaTime = ( time - startTime ) / duration;
				s = CMath::sin( deltaTime * CMath::HALF_PI ) * duration * 0.001f * CMath::SQRT_1OVER2;
				currentValue = startValue + deltaTime * baseSpeed + s * speed;
			}
			break;
		}
	}
	return currentValue;
}

/*
====================
CExtrapolation::GetCurrentSpeed
====================
*/
template< class type >
SMF_INLINE_FORCED type CExtrapolation<type>::getCurrentSpeed( float time ) const {
	float deltaTime, s;

	if ( time < startTime || !duration ) {
		return ( startValue - startValue );
	}

	if ( !( extrapolationType &	EXTRAP_NOSTOP ) && ( time > startTime + duration ) ) {
		return ( startValue - startValue );
	}

	switch( extrapolationType & ~EXTRAP_NOSTOP ) {
		case EXTRAP_NONE: {
			return baseSpeed;
		}
		case EXTRAP_LINEAR: {
			return baseSpeed + speed;
		}
		case EXTRAP_ACCELLINEAR: {
			deltaTime = ( time - startTime ) / duration;
			s = deltaTime;
			return baseSpeed + s * speed;
		}
		case EXTRAP_DECELLINEAR: {
			deltaTime = ( time - startTime ) / duration;
			s = 1.0f - deltaTime;
			return baseSpeed + s * speed;
		}
		case EXTRAP_ACCELSINE: {
			deltaTime = ( time - startTime ) / duration;
			s = CMath::sin( deltaTime * CMath::HALF_PI );
			return baseSpeed + s * speed;
		}
		case EXTRAP_DECELSINE: {
			deltaTime = ( time - startTime ) / duration;
			s = CMath::cos( deltaTime * CMath::HALF_PI );
			return baseSpeed + s * speed;
		}
		default: {
			return baseSpeed;
		}
	}
}

} //end MATH
} //end SMF
#endif /* !__MATH_EXTRAPOLATE_H__ */
