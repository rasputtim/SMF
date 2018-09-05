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

#ifndef _SMF__MATH_QUAT_H__
#define _SMF__MATH_QUAT_H__
#include "../SMF_Config.h"
#include "SMF_Vector.h"
#include "exceptions\all.h"

namespace SMF{
namespace MATH{



class CVec3D;
class CEulerAngles;
class CRotation;
class CMat3D;
class CMat4D;
class CCompQuaternion;
/**
 * \class CQuaternion
 *
 * \ingroup SMF_Math
 *
 * \if pt_br
 * \brief Quaternion
   \n
   Funções para acelerar operações com quaternions, comulmente utilizadas em
   calculos com números imaginários (complexos) \ref CComplex \ref e rotações 3D em animações de
   caracteres já que quaternions não podem são atingidas pelo travamento da rótula de 
   gimbal (gimbal lock) como nos ângulos de euler. 
   
   Quaternions representam rotação ao longo de um eixo arbitrário com um certo ângulo.
   Especificamente, podem ser definidos como:
	\f$ Qx = Ax * sin(r/2) \f$  ;
	\f$ Qy = Ay * sin(r/2) \f$  ;
	\f$ Qz = Az * sin(r/2) \f$  ;
	\f$ Qw = cos(r/2) \f$  ;
	Onde A representa um eixo de rotação 3D, r representa o ângulo em radianos e Q representa o quaternion
   \note Quaternions são armazenados no formato xyzw, e NÃO wxyz
 * \see  http://www.anticz.com/eularqua.htm 
 * \elseif us_en
 * \brief Quaternion
	The Quaternion module of SMF gives functions to accelerate quaternion operations commonly used 
	in calculations of imaginary numbers and 3D rotations in character animation that cannot suffer
	from Gimbal Lock. Quaternions are stored in xyzw format, NOT wxyz. 
	Quaternions represent rotation along an arbitrary axis with a certain angle. 
	Specifically, a quaternion can be defined as:
	 \f$  Qx = Ax * sin(r/2) \f$  ;
	 \f$  Qy = Ay * sin(r/2) \f$  ;
	 \f$  Qz = Az * sin(r/2) \f$  ;
	 \f$  Qw = cos(r/2) \f$  ;
	Where A represents the 3D axis of rotation, r represents the angle in radians, and Q represents the quaternion.

* \see http://www.youtube.com/watch?v=zc8b2Jo7mno , http://3dgep.com/understanding-quaternions/

 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CQuaternion {
public:
	union
	{
		struct
		{
			
			float x; ///< The factor of i.
			float y; ///< The factor of j. [similarOverload: x]
			float z; ///< The factor of k. [similarOverload: x]
			float w; ///< The scalar part. Sometimes also referred to as 'r'. [similarOverload: x]

		};
		sf_m128 q;
	};
	//float			x;
	//float			y;
	//float			z;
	//float			w;

					CQuaternion();
					CQuaternion( float x, float y, float z, float w );
	/** 
	\brief Constructs this quaternion by specifying a rotation axis and the amount of rotation to be performed
	about that axis.
	\param rotationAxis The normalized rotation axis to rotate about.
	\param rotationAngleRadians The angle to rotate by, in radians. For example, Pi/4.f equals to 45 degrees, Pi/2.f is 90 degrees, and Pi is 180 degrees.
	**/
	CQuaternion(const CVec3D &rotationAxis, float rotationAngleRadians);

	void 			set( float x, float y, float z, float w );

	float			operator[]( int index ) const;
	float &			operator[]( int index );
	CQuaternion			operator-() const;
	CQuaternion &		operator=( const CQuaternion &a );
	CQuaternion			operator+( const CQuaternion &a ) const;
	CQuaternion &		operator+=( const CQuaternion &a );
	CQuaternion			operator-( const CQuaternion &a ) const;
	CQuaternion &		operator-=( const CQuaternion &a );
	CQuaternion			operator*( const CQuaternion &a ) const;
	CVec3D			operator*( const CVec3D &a ) const;
	CQuaternion			operator*( float a ) const;
	CQuaternion &		operator*=( const CQuaternion &a );
	CQuaternion &		operator*=( float a );
	/// Divides a quaternion by another. Division "a / b" results in a quaternion that rotates the orientation b to coincide with the orientation a.
	CQuaternion operator /(const CQuaternion &rhs) const;

	friend CQuaternion	operator*( const float a, const CQuaternion &b );
	friend CVec3D	operator*( const CVec3D &a, const CQuaternion &b );

	bool			compare( const CQuaternion &a ) const;						// exact compare, no epsilon
	bool			compare( const CQuaternion &a, const float epsilon ) const;	// compare with epsilon
	bool			operator==(	const CQuaternion &a ) const;					// exact compare, no epsilon
	bool			operator!=(	const CQuaternion &a ) const;					// exact compare, no epsilon

	CQuaternion			Inverse() const;
	float			getLenght() const;
	CQuaternion &		toNormal();

	float			calcW() const;
	int				getDimension() const;

	CEulerAngles	toAngles() const;
	CRotation		toRotation() const;
	CMat3D			toMat3() const;
	CMat4D			toMat4() const;
	CMatJoint3x4	toMat3x4()const;
	CCompQuaternion	toQuat() const;
	CVec3D			toAngularVelocity() const;

	const float *	toFloatPtr() const;
	float *			toFloatPtr();
	const char *	toString( int precision = 2 ) const;

	CQuaternion &		slerp( const CQuaternion &from, const CQuaternion &to, float t );
	/// Rotates the given vector by this quaternion.
	MUST_USE_RESULT CVec3D transform(float x, float y, float z) const;
	MUST_USE_RESULT CVec3D transform(const CVec3D &vec) const;

	/// Rotates the given vector by this quaternion. The w component of the vector is assumed to be zero or one.
	MUST_USE_RESULT CVec4D transform(const CVec4D &vec) const;

		/// Transforms the given vector by this Quaternion.
	/// \note Technically, this function does not perform a simple multiplication of 'q * v',
	/// but instead performs a conjugation operation 'q*v*q^-1'. This corresponds to transforming
	/// the given vector by this Quaternion.
	CVec3D MUST_USE_RESULT mul(const CVec3D &vector) const{ return this->transform(vector); };
	CVec4D MUST_USE_RESULT mul(const CVec4D &vector) const{ return this->transform(vector); };
	/** Returns the rotation axis and angle of this quaternion.
	\param rotationAxis [out] Received the normalized axis of the rotation.
	\param rotationAngleRadians [out] Receives the angle of rotation around the given axis. This parameter is returned in the range [0, 2pi].
	**/
	void toAxisAngle(CVec3D &rotationAxis, float &rotationAngleRadians) const;
	/**
	\brief Sets this quaternion by specifying the axis about which the rotation is performed, and the angle of rotation.
	\param rotationAxis The axis of rotation. This vector must be normalized to call this function.
	\param rotationAngleRadians The angle of rotation in radians.
	**/
	void setFromAxisAngle(const CVec3D &rotationAxis, float rotationAngleRadians);

		/// Returns true if the length of this quaternion is one.
	bool isNormalized(float epsilon =CMath::EPSILON_SuperLow) const;

	bool isInvertible(float epsilon =CMath::EPSILON_SuperLow) const;
	float getLengthSqr() const;

	/// Returns true if the entries of this quaternion are all finite.
	bool isFinite() const;

};

SMF_INLINE_FORCED CQuaternion::CQuaternion() {
}

SMF_INLINE_FORCED CQuaternion::CQuaternion( float x, float y, float z, float w ) {
	this->x = x;
	this->y = y;
	this->z = z;
	this->w = w;
}
SMF_INLINE_FORCED CQuaternion::CQuaternion(const CVec3D &rotationAxis, float rotationAngle)
{
	setFromAxisAngle(rotationAxis, rotationAngle);
}
SMF_INLINE_FORCED float CQuaternion::operator[]( int index ) const {
	SMF_ASSERT( ( index >= 0 ) && ( index < 4 ) );
	return ( &x )[ index ];
}

SMF_INLINE_FORCED float& CQuaternion::operator[]( int index ) {
	SMF_ASSERT( ( index >= 0 ) && ( index < 4 ) );
	return ( &x )[ index ];
}

SMF_INLINE_FORCED CQuaternion CQuaternion::operator-() const {
	return CQuaternion( -x, -y, -z, -w );
}

SMF_INLINE_FORCED CQuaternion &CQuaternion::operator=( const CQuaternion &a ) {
	x = a.x;
	y = a.y;
	z = a.z;
	w = a.w;

	return *this;
}

SMF_INLINE_FORCED CQuaternion CQuaternion::operator+( const CQuaternion &a ) const {
	return CQuaternion( x + a.x, y + a.y, z + a.z, w + a.w );
}

SMF_INLINE_FORCED CQuaternion& CQuaternion::operator+=( const CQuaternion &a ) {
	x += a.x;
	y += a.y;
	z += a.z;
	w += a.w;

	return *this;
}

SMF_INLINE_FORCED CQuaternion CQuaternion::operator-( const CQuaternion &a ) const {
	return CQuaternion( x - a.x, y - a.y, z - a.z, w - a.w );
}

SMF_INLINE_FORCED CQuaternion& CQuaternion::operator-=( const CQuaternion &a ) {
	x -= a.x;
	y -= a.y;
	z -= a.z;
	w -= a.w;

	return *this;
}


SMF_INLINE_FORCED CVec3D CQuaternion::operator*( const CVec3D &a ) const {
#if 0
	// it's faster to do the conversion to a 3x3 matrix and multiply the vector by this 3x3 matrix
	return ( toMat3() * a );
#else
	// result = this->Inverse() * CQuaternion( a.x, a.y, a.z, 0.0f ) * (*this)
	float xxzz = x*x - z*z;
	float wwyy = w*w - y*y;

	float xw2 = x*w*2.0f;
	float xy2 = x*y*2.0f;
	float xz2 = x*z*2.0f;
	float yw2 = y*w*2.0f;
	float yz2 = y*z*2.0f;
	float zw2 = z*w*2.0f;

	return CVec3D(
		(xxzz + wwyy)*a.x		+ (xy2 + zw2)*a.y		+ (xz2 - yw2)*a.z,
		(xy2 - zw2)*a.x			+ (y*y+w*w-x*x-z*z)*a.y	+ (yz2 + xw2)*a.z,
		(xz2 + yw2)*a.x			+ (yz2 - xw2)*a.y		+ (wwyy - xxzz)*a.z
	);
#endif
}

SMF_INLINE_FORCED CQuaternion CQuaternion::operator*( float a ) const {
	return CQuaternion( x * a, y * a, z * a, w * a );
}

SMF_INLINE_FORCED CQuaternion operator*( const float a, const CQuaternion &b ) {
	return b * a;
}

SMF_INLINE_FORCED CVec3D operator*( const CVec3D &a, const CQuaternion &b ) {
	return b * a;
}

SMF_INLINE_FORCED CQuaternion& CQuaternion::operator*=( const CQuaternion &a ) {
	*this = *this * a;

	return *this;
}

SMF_INLINE_FORCED CQuaternion& CQuaternion::operator*=( float a ) {
	x *= a;
	y *= a;
	z *= a;
	w *= a;

	return *this;
}

SMF_INLINE_FORCED bool CQuaternion::compare( const CQuaternion &a ) const {
	return ( ( x == a.x ) && ( y == a.y ) && ( z == a.z ) && ( w == a.w ) );
}

SMF_INLINE_FORCED bool CQuaternion::compare( const CQuaternion &a, const float epsilon ) const {
	if ( CMath::fabs( x - a.x ) > epsilon ) {
		return false;
	}
	if ( CMath::fabs( y - a.y ) > epsilon ) {
		return false;
	}
	if ( CMath::fabs( z - a.z ) > epsilon ) {
		return false;
	}
	if ( CMath::fabs( w - a.w ) > epsilon ) {
		return false;
	}
	return true;
}

SMF_INLINE_FORCED bool CQuaternion::operator==( const CQuaternion &a ) const {
	return compare( a );
}

SMF_INLINE_FORCED bool CQuaternion::operator!=( const CQuaternion &a ) const {
	return !compare( a );
}

SMF_INLINE_FORCED void CQuaternion::set( float x, float y, float z, float w ) {
	this->x = x;
	this->y = y;
	this->z = z;
	this->w = w;
}

SMF_INLINE_FORCED CQuaternion CQuaternion::Inverse() const {
	return CQuaternion( -x, -y, -z, w );
}

SMF_INLINE_FORCED float CQuaternion::getLenght() const {
	float len;

	len = x * x + y * y + z * z + w * w;
	return CMath::sqrt( len );
}

SMF_INLINE_FORCED CQuaternion& CQuaternion::toNormal() {
	float len;
	float ilength;

	len = this->getLenght();
	if ( len ) {
		ilength = 1 / len;
		x *= ilength;
		y *= ilength;
		z *= ilength;
		w *= ilength;
	}
	return *this;
}

SMF_INLINE_FORCED float CQuaternion::calcW() const {
	// take the absolute value because floating point rounding may cause the dot of x,y,z to be larger than 1
	return sqrt( CMath::fabs( 1.0f - ( x * x + y * y + z * z ) ) );
}

SMF_INLINE_FORCED int CQuaternion::getDimension() const {
	return 4;
}

SMF_INLINE_FORCED const float *CQuaternion::toFloatPtr() const {
	return &x;
}

SMF_INLINE_FORCED float *CQuaternion::toFloatPtr() {
	return &x;
}



/**
 * \class CCompQuaternion
 *
 * \ingroup SMF_Math
 *
 * \if pt_br
 * \brief Quaternion Comprimido
   \n
   Funções para acelerar operações com quaternions, comulmente utilizadas em
   calculos com números imaginários (complexos) \ref CComplex \ref e rotações 3D em animações de
   caracteres que não podem são atingidas pelo travamento da rótula de gimbal (gimbal lock). 
   Quaternions representam rotação ao longo de um eixo arbitrário com um certo ângulo.
   Especificamente, podem ser definidos como:
	Qx = Ax * sin(r/2);
	Qy = Ay * sin(r/2);
	Qz = Az * sin(r/2);
	Qw = cos(r/2);
	Onde A representa um eixo de rotação 3D, r representa o ângulo em radianos e Q representa o quaternion
   \note Quaternions são armazenados no formato xyzw, e NÃO wxyz
 * \see  http://www.anticz.com/eularqua.htm 
 * \elseif us_en
 * \brief Compressed Quaternion 
	The Quaternion module of SMF gives functions to accelerate quaternion operations commonly used 
	in calculations of imaginary numbers and 3D rotations in character animation that cannot suffer
	from Gimbal Lock. Quaternions are stored in xyzw format, NOT wxyz. 
	Quaternions represent rotation along an arbitrary axis with a certain angle. 
	Specifically, a quaternion can be defined as:
	Qx = Ax * sin(r/2);
	Qy = Ay * sin(r/2);
	Qz = Az * sin(r/2);
	Qw = cos(r/2);
	Where A represents the 3D axis of rotation, r represents the angle in radians, and Q represents the quaternion.

* \see http://3dgep.com/understanding-quaternions/

 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */

class SMF_API CCompQuaternion {
public:
	float			x;
	float			y;
	float			z;

					CCompQuaternion();
					CCompQuaternion( float x, float y, float z );

	void 			set( float x, float y, float z );

	float			operator[]( int index ) const;
	float &			operator[]( int index );

	bool			compare( const CCompQuaternion &a ) const;						// exact compare, no epsilon
	bool			compare( const CCompQuaternion &a, const float epsilon ) const;	// compare with epsilon
	bool			operator==(	const CCompQuaternion &a ) const;					// exact compare, no epsilon
	bool			operator!=(	const CCompQuaternion &a ) const;					// exact compare, no epsilon

	int				getDimension() const;

	CEulerAngles	toAngles() const;
	CRotation		toRotation() const;
	CMat3D			toMat3() const;
	CMat4D			toMat4() const;
	CQuaternion			toQuat() const;
	const float *	toFloatPtr() const;
	float *			toFloatPtr();
	const char *	toString( int precision = 2 ) const;
};

SMF_INLINE_FORCED CCompQuaternion::CCompQuaternion() {
}

SMF_INLINE_FORCED CCompQuaternion::CCompQuaternion( float x, float y, float z ) {
	this->x = x;
	this->y = y;
	this->z = z;
}

SMF_INLINE_FORCED void CCompQuaternion::set( float x, float y, float z ) {
	this->x = x;
	this->y = y;
	this->z = z;
}

SMF_INLINE_FORCED float CCompQuaternion::operator[]( int index ) const {
	SMF_ASSERT( ( index >= 0 ) && ( index < 3 ) );
	return ( &x )[ index ];
}

SMF_INLINE_FORCED float& CCompQuaternion::operator[]( int index ) {
	SMF_ASSERT( ( index >= 0 ) && ( index < 3 ) );
	return ( &x )[ index ];
}

SMF_INLINE_FORCED bool CCompQuaternion::compare( const CCompQuaternion &a ) const {
	return ( ( x == a.x ) && ( y == a.y ) && ( z == a.z ) );
}

SMF_INLINE_FORCED bool CCompQuaternion::compare( const CCompQuaternion &a, const float epsilon ) const {
	if ( CMath::fabs( x - a.x ) > epsilon ) {
		return false;
	}
	if ( CMath::fabs( y - a.y ) > epsilon ) {
		return false;
	}
	if ( CMath::fabs( z - a.z ) > epsilon ) {
		return false;
	}
	return true;
}

SMF_INLINE_FORCED bool CCompQuaternion::operator==( const CCompQuaternion &a ) const {
	return compare( a );
}

SMF_INLINE_FORCED bool CCompQuaternion::operator!=( const CCompQuaternion &a ) const {
	return !compare( a );
}

SMF_INLINE_FORCED int CCompQuaternion::getDimension() const {
	return 3;
}

SMF_INLINE_FORCED CQuaternion CCompQuaternion::toQuat() const {
	// take the absolute value because floating point rounding may cause the dot of x,y,z to be larger than 1
	return CQuaternion( x, y, z, sqrt( fabs( 1.0f - ( x * x + y * y + z * z ) ) ) );
}

SMF_INLINE_FORCED const float *CCompQuaternion::toFloatPtr() const {
	return &x;
}

SMF_INLINE_FORCED float *CCompQuaternion::toFloatPtr() {
	return &x;
}

} //end MATH
} //end SMF
#endif /* !__MATH_QUAT_H__ */
