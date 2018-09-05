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

#ifndef _SMF__MATH_ANGLES_H__
#define _SMF__MATH_ANGLES_H__
#include "../SMF_Config.h"
#include "SMF_Vector.h"


namespace SMF {
namespace MATH{
// angle indexes
#define	PITCH				0		// up / down
#define	YAW					1		// left / right
#define	ROLL				2		// fall over


class CVec3D;
class CQuaternion;
class CRotation;
class CMat3D;
class CMat4D;

/**
 * \class CEulerAngles
 *
 * \ingroup SMF_Math
 *
 * \brief Implementa Ângulos de Euler
 *
 * \note http://pt.wikipedia.org/wiki/%C3%82ngulos_de_Euler
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
  */
class SMF_API  CEulerAngles {
public:
	float			pitch;
	float			yaw;
	float			roll;

					CEulerAngles();
					/**
					\brief cria um objeto através das variáveis passadas
					\param pitch ângulo de pith em graus
					\param yaw ângulo de yaw em graus
					\param roll ângulo de roll em graus
					**/
					CEulerAngles( float pitch, float yaw, float roll );
					/**
					\brief cria um objeto através do vetor passado
					\param v vetor contendo os ângulos yaw pitch e roll
					**/
					explicit CEulerAngles( const CVec3D &v );
	/**
	\brief seta os ângulos do objeto atual
	\param pitch ângulo de pith em graus
	\param yaw ângulo de yaw em graus
	\param roll ângulo de roll em graus
	**/
	void 			setAngles( float pitch, float yaw, float roll );
	/**
	\brief zera os angulos do objeto
	**/
	CEulerAngles &	toZero();
	float				operator[]( int index ) const;
	float &				operator[]( int index );
	CEulerAngles		operator-() const;			// negate angles, in general not the inverse rotation
	CEulerAngles &		operator=( const CEulerAngles &a );
	CEulerAngles		operator+( const CEulerAngles &a ) const;
	CEulerAngles &		operator+=( const CEulerAngles &a );
	CEulerAngles		operator-( const CEulerAngles &a ) const;
	CEulerAngles &		operator-=( const CEulerAngles &a );
	CEulerAngles		operator*( const float a ) const;
	CEulerAngles &		operator*=( const float a );
	CEulerAngles		operator/( const float a ) const;
	CEulerAngles &		operator/=( const float a );

	friend CEulerAngles	operator*( const float a, const CEulerAngles &b );
	/**
	\brief comparação exata, sem episilon
	\return retorna verdadeiro se os dois ângulos são iguais e falso caso contrário
	**/
	bool			compare( const CEulerAngles &a ) const;							// exact compare, no epsilon
	/**
	\brief comparação com episilon
	\return retorna verdadeiro se os dois ângulos são iguais e falso caso contrário
	**/
	bool			compare( const CEulerAngles &a, const float epsilon ) const;	// compare with epsilon
	/**
	\brief comparação exata, sem episilon
	\return retorna verdadeiro se os dois ângulos são iguais e falso caso contrário
	**/
	bool			operator==(	const CEulerAngles &a ) const;						// exact compare, no epsilon
	/**
	\brief comparação exata, sem episilon
	\return retorna verdadeiro se os dois Ângulos são diferentes e falso caso contrário
	**/
	bool			operator!=(	const CEulerAngles &a ) const;						// exact compare, no epsilon
	/**
	\brief  retorna ângulos normalizados ao intervalo [0 <= angle < 360]
	**/
	CEulerAngles &		normalize360();	// normalizes 'this'
	/**
	\brief retorna ângulos normalizados ao intervalo [-180 < angle <= 180]
	**/
	CEulerAngles &		normalize180();	// normalizes 'this'

	void			clamp( const CEulerAngles &min, const CEulerAngles &max );
	/**
	\brief retorna a dimensão do ângulo de Euler
	**/
	int				getDimension() const;
	/**
	\brief calcula os vetores correspondentes que descrevem o Angulo
	\param forward vetor que será preenchido com os dados do forward
	\param right vetor que será preenchido com os dados do right
	\param up vetor que será preenchido com os dados do up
	**/
	void			toVectors( CVec3D *forward, CVec3D *right = NULL, CVec3D *up = NULL ) const;
	/**
	\brief calcula o vetor forward
	\return retorna um vetor que descreve o forward
	**/
	CVec3D		toForward() const;
	/**
	\brief converte para Quaternon
	**/
	CQuaternion		toQuat() const;
	/**
	\brief converte para Rotation
	**/
	CRotation		toRotation() const;
	/**
	\brief converte para matriz
	**/
	CMat3D		toMat3() const;
	/**
	\brief converte para matriz
	**/
	CMat4D		toMat4() const;
	/**
	\brief converte para vetor
	**/
	CVec3D		toAngularVelocity() const;
	/**
	\brief retorna um ponteiro float para o primeiro ângulo (pith)
	**/
	const float *	toFloatPtr() const;
	/**
	\brief retorna um ponteiro float para o primeiro ângulo (pith)
	**/
	float *			toFloatPtr();
	const char *	toString( int precision = 2 ) const;
};

extern CEulerAngles ang_zero;

SMF_INLINE_FORCED CEulerAngles::CEulerAngles() {
}

SMF_INLINE_FORCED CEulerAngles::CEulerAngles( float pitch, float yaw, float roll ) {
	this->pitch = pitch;
	this->yaw	= yaw;
	this->roll	= roll;
}

SMF_INLINE_FORCED CEulerAngles::CEulerAngles( const CVec3D &v ) {
	this->pitch = v[0];
	this->yaw	= v[1];
	this->roll	= v[2];
}

SMF_INLINE_FORCED void CEulerAngles::setAngles( float pitch, float yaw, float roll ) {
	this->pitch = pitch;
	this->yaw	= yaw;
	this->roll	= roll;
}

SMF_INLINE_FORCED CEulerAngles &CEulerAngles::toZero() {
	pitch = yaw = roll = 0.0f;
	return *this;
}

SMF_INLINE_FORCED float CEulerAngles::operator[]( int index ) const {
	SMF_ASSERT( ( index >= 0 ) && ( index < 3 ) );
	return ( &pitch )[ index ];
}

SMF_INLINE_FORCED float &CEulerAngles::operator[]( int index ) {
	SMF_ASSERT( ( index >= 0 ) && ( index < 3 ) );
	return ( &pitch )[ index ];
}

SMF_INLINE_FORCED CEulerAngles CEulerAngles::operator-() const {
	return CEulerAngles( -pitch, -yaw, -roll );
}

SMF_INLINE_FORCED CEulerAngles &CEulerAngles::operator=( const CEulerAngles &a ) {
	pitch	= a.pitch;
	yaw		= a.yaw;
	roll	= a.roll;
	return *this;
}

SMF_INLINE_FORCED CEulerAngles CEulerAngles::operator+( const CEulerAngles &a ) const {
	return CEulerAngles( pitch + a.pitch, yaw + a.yaw, roll + a.roll );
}

SMF_INLINE_FORCED CEulerAngles& CEulerAngles::operator+=( const CEulerAngles &a ) {
	pitch	+= a.pitch;
	yaw		+= a.yaw;
	roll	+= a.roll;

	return *this;
}

SMF_INLINE_FORCED CEulerAngles CEulerAngles::operator-( const CEulerAngles &a ) const {
	return CEulerAngles( pitch - a.pitch, yaw - a.yaw, roll - a.roll );
}

SMF_INLINE_FORCED CEulerAngles& CEulerAngles::operator-=( const CEulerAngles &a ) {
	pitch	-= a.pitch;
	yaw		-= a.yaw;
	roll	-= a.roll;

	return *this;
}

SMF_INLINE_FORCED CEulerAngles CEulerAngles::operator*( const float a ) const {
	return CEulerAngles( pitch * a, yaw * a, roll * a );
}

SMF_INLINE_FORCED CEulerAngles& CEulerAngles::operator*=( float a ) {
	pitch	*= a;
	yaw		*= a;
	roll	*= a;
	return *this;
}

SMF_INLINE_FORCED CEulerAngles CEulerAngles::operator/( const float a ) const {
	float inva = 1.0f / a;
	return CEulerAngles( pitch * inva, yaw * inva, roll * inva );
}

SMF_INLINE_FORCED CEulerAngles& CEulerAngles::operator/=( float a ) {
	float inva = 1.0f / a;
	pitch	*= inva;
	yaw		*= inva;
	roll	*= inva;
	return *this;
}

SMF_INLINE_FORCED CEulerAngles operator*( const float a, const CEulerAngles &b ) {
	return CEulerAngles( a * b.pitch, a * b.yaw, a * b.roll );
}

SMF_INLINE_FORCED bool CEulerAngles::compare( const CEulerAngles &a ) const {
	return ( ( a.pitch == pitch ) && ( a.yaw == yaw ) && ( a.roll == roll ) );
}

SMF_INLINE_FORCED bool CEulerAngles::compare( const CEulerAngles &a, const float epsilon ) const {
	if ( CMath::fabs( pitch - a.pitch ) > epsilon ) {
		return false;
	}
			
	if ( CMath::fabs( yaw - a.yaw ) > epsilon ) {
		return false;
	}

	if ( CMath::fabs( roll - a.roll ) > epsilon ) {
		return false;
	}

	return true;
}

SMF_INLINE_FORCED bool CEulerAngles::operator==( const CEulerAngles &a ) const {
	return compare( a );
}

SMF_INLINE_FORCED bool CEulerAngles::operator!=( const CEulerAngles &a ) const {
	return !compare( a );
}

SMF_INLINE_FORCED void CEulerAngles::clamp( const CEulerAngles &min, const CEulerAngles &max ) {
	if ( pitch < min.pitch ) {
		pitch = min.pitch;
	} else if ( pitch > max.pitch ) {
		pitch = max.pitch;
	}
	if ( yaw < min.yaw ) {
		yaw = min.yaw;
	} else if ( yaw > max.yaw ) {
		yaw = max.yaw;
	}
	if ( roll < min.roll ) {
		roll = min.roll;
	} else if ( roll > max.roll ) {
		roll = max.roll;
	}
}

SMF_INLINE_FORCED int CEulerAngles::getDimension() const {
	return 3;
}

SMF_INLINE_FORCED const float *CEulerAngles::toFloatPtr() const {
	return &pitch;
}

SMF_INLINE_FORCED float *CEulerAngles::toFloatPtr() {
	return &pitch;
}

} //end MATH
} //end SMF
#endif /* !__MATH_ANGLES_H__ */
