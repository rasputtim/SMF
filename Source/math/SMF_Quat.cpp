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


#include "math/SMF_Quaternion.h"
#include "math/SMF_Vector.h"
#include "math/SMF_Matriz.h"
#include "math/SMF_EulerAngles.h"
#include "util/SMF_StringUtils.h"
#include "math/SMF_Rotation.h"
#include "math/SMF_JointTransform.h"


namespace SMF {
namespace MATH{

CQuaternion CQuaternion::operator /(const CQuaternion &r) const
{
if (CMath::MATH_AUTOMATIC_SSE){
    sf_m128 ret=CSIMD::getProcessor()->quat_div_quat(q, r.q);
	return *reinterpret_cast<CQuaternion *>(&ret);
}else{
    sf_m128 ret=CSIMD::getGenProcessor()->quat_div_quat(q, r.q);
	return *reinterpret_cast<CQuaternion *>(&ret);
}
}

CQuaternion CQuaternion::operator*( const CQuaternion &a ) const {
    sf_m128 ret;
	if (CMath::MATH_AUTOMATIC_SSE) {
	    ret=CSIMD::getProcessor()->quat_mul_quat(q, a.q);
	return *reinterpret_cast<CQuaternion *>(&ret);
}else{
    ret=CSIMD::getGenProcessor()->quat_mul_quat(q, a.q);
	return *reinterpret_cast<CQuaternion *>(&ret);
}

}


CEulerAngles CQuaternion::toAngles( void ) const {
	return toMat3().toAngles();
}

/*
=====================
CQuaternion::toRotation
=====================
*/
CRotation CQuaternion::toRotation( void ) const {
	CVec3D vec;
	float angle;

	vec.x = x;
	vec.y = y;
	vec.z = z;
	angle = CMath::acos( w );
	if ( angle == 0.0f ) {
		vec.set( 0.0f, 0.0f, 1.0f );
	} else {
		//vec *= (1.0f / sin( angle ));
		vec.toNormal();
		vec.fixDegenerateNormal();
		angle *= 2.0f * CMath::M_RAD2DEG;
	}
	return CRotation( CVec3D::origin, vec, angle );
}

/*
=====================
CQuaternion::toMat3
=====================
*/
CMat3D CQuaternion::toMat3( void ) const {
	CMat3D	mat;
	float	wx, wy, wz;
	float	xx, yy, yz;
	float	xy, xz, zz;
	float	x2, y2, z2;

	x2 = x + x;
	y2 = y + y;
	z2 = z + z;

	xx = x * x2;
	xy = x * y2;
	xz = x * z2;

	yy = y * y2;
	yz = y * z2;
	zz = z * z2;

	wx = w * x2;
	wy = w * y2;
	wz = w * z2;
	//col 0
	mat[ 0 ][ 0 ] = 1.0f - ( yy + zz );
	mat[ 0 ][ 1 ] = xy - wz;
	mat[ 0 ][ 2 ] = xz + wy;
	//col 1
	mat[ 1 ][ 0 ] = xy + wz;
	mat[ 1 ][ 1 ] = 1.0f - ( xx + zz );
	mat[ 1 ][ 2 ] = yz - wx;
	//col 2
	mat[ 2 ][ 0 ] = xz - wy;
	mat[ 2 ][ 1 ] = yz + wx;
	mat[ 2 ][ 2 ] = 1.0f - ( xx + yy );

	return mat;
}

void CQuaternion::toAxisAngle(CVec3D &axis, float &angle) const
{
	angle = CMath::acos(w) * 2.f;
	float sinz = CMath::sin(angle/2.f);
	if (fabs(sinz) > 1e-4f)
	{
		sinz = 1.f / sinz;
		axis = CVec3D(x * sinz, y * sinz, z * sinz);
	}
	else
	{
		// The quaternion does not produce any rotation. Still, explicitly
		// set the axis so that the user gets a valid normalized vector back.
		angle = 0.f;
		axis = CVec3D(1.f, 0.f, 0.f);
	}
}

void CQuaternion::setFromAxisAngle(const CVec3D &axis, float angle)
{
	//s assume(axis.isNormalized());
	//s assume(MATH_NS::isFinite(angle));
	float cosz = CMath::cos(angle/2.f);
	float sinz = CMath::sin(angle/2.f);
	x = axis.x * sinz;
	y = axis.y * sinz;
	z = axis.z * sinz;
	w = cosz;
}



CMatJoint3x4	CQuaternion::toMat3x4()const{
CMatJoint3x4 m;
if(CMath::MATH_AUTOMATIC_SSE){
	CSIMD::getProcessor()->quat_to_mat3x4(q, MATH::SET_PS(1,0,0,0), m.Row(0).toM128Ptr());
}else{
	CSIMD::getGenProcessor()->quat_to_mat3x4(q, MATH::SET_PS(1,0,0,0), m.Row(0).toM128Ptr());
}
return m;

}
CMat4D CQuaternion::toMat4( void ) const {
	return toMat3().toMat4();
}


CCompQuaternion CQuaternion::toQuat( void ) const {
	if ( w < 0.0f ) {
		return CCompQuaternion( -x, -y, -z );
	}
	return CCompQuaternion( x, y, z );
}

/*
============
CQuaternion::toAngularVelocity
============
*/
CVec3D CQuaternion::toAngularVelocity( void ) const {
	CVec3D vec;

	vec.x = x;
	vec.y = y;
	vec.z = z;
	vec.toNormal();
	return vec * CMath::acos( w );
}

/*
=============
CQuaternion::toString
=============
*/
const char *CQuaternion::toString( int precision ) const {
	return CMyString::floatArraytoString( toFloatPtr(), getDimension(), precision );
}

/*
=====================
CQuaternion::slerp

Spherical linear interpolation between two quaternions.
=====================
*/
CQuaternion &CQuaternion::slerp( const CQuaternion &from, const CQuaternion &to, float t ) {
	CQuaternion	temp;
	float	omega, cosom, sinom, scale0, scale1;

	if ( t <= 0.0f ) {
		*this = from;
		return *this;
	}

	if ( t >= 1.0f ) {
		*this = to;
		return *this;
	}

	if ( from == to ) {
		*this = to;
		return *this;
	}

	cosom = from.x * to.x + from.y * to.y + from.z * to.z + from.w * to.w;
	if ( cosom < 0.0f ) {
		temp = -to;
		cosom = -cosom;
	} else {
		temp = to;
	}

	if ( ( 1.0f - cosom ) > 1e-6f ) {
#if 0
		omega = CMath::acos( cosom );
		sinom = 1.0f / sin( omega );
		scale0 = sin( ( 1.0f - t ) * omega ) * sinom;
		scale1 = sin( t * omega ) * sinom;
#else
		scale0 = 1.0f - cosom * cosom;
		sinom = CMath::invSqrt( scale0 );
		omega = CMath::atan16( scale0 * sinom, cosom );
		scale0 = CMath::sin16( ( 1.0f - t ) * omega ) * sinom;
		scale1 = CMath::sin16( t * omega ) * sinom;
#endif
	} else {
		scale0 = 1.0f - t;
		scale1 = t;
	}

	*this = ( scale0 * from ) + ( scale1 * temp );
	return *this;
}

/*
=============
CCompQuaternion::toAngles
=============
*/
CEulerAngles CCompQuaternion::toAngles( void ) const {
	return toQuat().toAngles();
}

/*
=============
CCompQuaternion::toRotation
=============
*/
CRotation CCompQuaternion::toRotation( void ) const {
	return toQuat().toRotation();
}

/*
=============
CCompQuaternion::toMat3
=============
*/
CMat3D CCompQuaternion::toMat3( void ) const {
	return toQuat().toMat3();
}

/*
=============
CCompQuaternion::toMat4
=============
*/
CMat4D CCompQuaternion::toMat4( void ) const {
	return toQuat().toMat4();
}

CVec3D MUST_USE_RESULT CQuaternion::transform(const CVec3D &vec) const
{
	SMF_ASSERT(this->isNormalized());
if (CMath::MATH_AUTOMATIC_SSE) {
	///\todo Check the generation of temporaries here!
	return CVec4D(CSIMD::getProcessor()->quat_transform_vec4(q, CVec4D(vec,0.f).v)).xyz();
}else{
	///\todo Optimize/benchmark the scalar path not to generate a matrix!
	//CMat3D mat = this->toMat3();
	//return mat * vec;
	return CVec4D(CSIMD::getGenProcessor()->quat_transform_vec4(q, CVec4D(vec,0.f).v)).xyz();
}
}

CVec3D MUST_USE_RESULT CQuaternion::transform(float x, float y, float z) const
{
if (CMath::MATH_AUTOMATIC_SSE) {
	///\todo Check the generation of temporaries here!
	return CVec4D(CSIMD::getProcessor()->quat_transform_vec4(q, CVec4D(x,y,z,0.f).v)).xyz();
}else{
	return transform(CVec3D(x, y, z));
}
}

CVec4D MUST_USE_RESULT CQuaternion::transform(const CVec4D &vec) const
{
	SMF_ASSERT(vec.isWZeroOrOne());

if (CMath::MATH_AUTOMATIC_SSE) {
	return CSIMD::getProcessor()->quat_transform_vec4(q, vec.v);
}else{

	return CVec4D(transform(vec.x, vec.y, vec.z), vec.w);
}
}



/*
=============
CCompQuaternion::toString
=============
*/
const char *CCompQuaternion::toString( int precision ) const {
	return CMyString::floatArraytoString( toFloatPtr(), getDimension(), precision );
}
float CQuaternion::getLengthSqr() const
{
	return x*x + y*y + z*z + w*w;

}

bool CQuaternion::isFinite() const
{
	return MATH::isFinite(x) && MATH::isFinite(y) && MATH::isFinite(z) && MATH::isFinite(w);
}

bool CQuaternion::isNormalized(float epsilon) const
{
	return CMath::equalsAbs(getLengthSqr(), 1.f, epsilon);
}

bool CQuaternion::isInvertible(float epsilon) const
{
	return getLengthSqr() > epsilon && isFinite();
}

} //END MATH
} //end SMF
