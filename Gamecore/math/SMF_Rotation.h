/*
  SMF -  Super Math Fabric  (https:/** \briefsourceforge.net/promects/sMfabric/)
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

#ifndef _SMF__MATH_ROTATION_H__

#define _SMF__MATH_ROTATION_H__

#include "../SMF_Config.h"
#include "SMF_Vector.h"
#include "SMF_Matriz.h"
#include "SMF_EulerAngles.h"
#include "../util/SMF_StringUtils.h"



namespace SMF {
	namespace MATH{
class CEulerAngles;
class CQuaternion;
class CMat3D;

/**
 * \class CRotation
 *
 * \ingroup SMF_Math
 * 
 * \if pt_br
 * \brief    Descreve uma rotação completa em graus em torno de um eixo arbitrário em 3D.
 *		  \n Uma matriz de rotação local é armazenada para rotação rápida de multiplos pontos
 * \elseif us_en
 * \brief 	Describes a complete 3D rotation in degrees about an abritray axis.
	        \n A local rotation matrix is stored for fast rotation of multiple points.
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */

class SMF_API CRotation {

	friend class CEulerAngles;
	friend class CQuaternion;
	friend class CMat3D;

public:
						CRotation();
						CRotation( const CVec3D &rotationOrigin, const CVec3D &rotationVec, const float rotationAngle );

	void				set( const CVec3D &rotationOrigin, const CVec3D &rotationVec, const float rotationAngle );
	void				setOrigin( const CVec3D &rotationOrigin );
    /** \brief setup rotationvector
	\note has to be normalized
	*/
	void				setVec( const CVec3D &rotationVec );					
	/** \brief setup rotationvector
	\note has to be normalized
	*/
	void				setVec( const float x, const float y, const float z );
	void				setAngle( const float rotationAngle );
	void				scale( const float s );
	void				reCalculateMatrix();
	const CVec3D &		getOrigin() const;
	const CVec3D &		getVec() const;
	float				getAngle() const;
	/** \brief flips rotation
	*/
	CRotation			operator-() const;
	/** \brief scale rotation */
	CRotation			operator*( const float s ) const;						
	/** \brief scale rotation*/
	CRotation			operator/( const float s ) const;						
	/** \brief scale rotation*/
	CRotation &			operator*=( const float s );							
	/** \brief scale rotation*/
	CRotation &			operator/=( const float s );							
	/** \brief rotate vector*/
	CVec3D				operator*( const CVec3D &v ) const;						
    /** \brief scale rotation*/
	friend CRotation	operator*( const float s, const CRotation &r );			
	/** \brief rotate vector*/
	friend CVec3D		operator*( const CVec3D &v, const CRotation &r );		
    /** \brief rotate vector*/
	friend CVec3D &		operator*=( CVec3D &v, const CRotation &r );			

	CEulerAngles		toAngles() const;
	CQuaternion			toQuat() const;
	const CMat3D &		toMat3() const;
	CMat4D				toMat4() const;
	CVec3D				toAngularVelocity() const;
	/**
	\brief rotate the point
	\note rotation formula: \f$ axix * (point - origin) + origin  \f$
	\n where:
	\n axis as the matrix: 
	\n  \f[   \left( \begin{array}{ccc}
			1-(yy+zz) & xy-wz & xz+wy \\
			zy+wz & 1-(xx+zz) & yz-wx \\
			xz-wy & yz+wx & 1-(xx+yy) \end{array} \right) 
	\f]
	\n point to be rotated
	\n origin origin of rotation
	\n \f$ xx = [vec.x * sin(angle)]*[2*(vec.x * sin(angle) )] \f$ 
	\n \f$  xy = [vec.x * sin(angle)]*[2*(vec.y * sin(angle) )] \f$ 
	\n \f$  xz = [vec.x * sin(angle)]*[2*(vec.z * sin(angle) )] \f$ 
	\n \f$  yy = [vec.y * sin(angle)]*[2*(vec.x * sin(angle) )] \f$ 
	\n \f$  yz = [vec.y * sin(angle)]*[2*(vec.z * sin(angle) )] \f$ 
	\n \f$  zz = [vec.z * sin(angle)]*[2*(vec.z * sin(angle) )] \f$ 
	\n \f$  wx = cos(angle)*[2*(vec.x * sin(angle) )] \f$ 
	\n \f$  wy = cos(angle)*[2*(vec.y * sin(angle) )] \f$ 
	\n \f$  wz = cos(angle)*[2*(vec.z * sin(angle) )] \f$ 
	\n 
	\note factor out the -origin part, it becomes:

	\n \f$ rotation * point + ( - rotation * origin + origin ); \f$
	\n where the second part is point-invariant and can be precomputed.
	*/
	void				rotatePoint( CVec3D &point ) const;

	void				normalize180();
	void				normalize360();

private:
    /** \brief origin of rotation*/
	CVec3D			origin;		
    /** \brief normalized vector to rotate around*/
	CVec3D			vec;	
    /** \brief angle of rotation in degrees*/
	float			angle;		
	/** \brief rotation axis*/
	mutable CMat3D	axis;			
	/** \brief true if rotation axis is valid*/
	mutable bool	axisValid;		
};


SMF_INLINE_FORCED CRotation::CRotation() {
}

SMF_INLINE_FORCED CRotation::CRotation( const CVec3D &rotationOrigin, const CVec3D &rotationVec, const float rotationAngle ) {
	origin = rotationOrigin;
	vec = rotationVec;
	angle = rotationAngle;
	axisValid = false;
}

SMF_INLINE_FORCED void CRotation::set( const CVec3D &rotationOrigin, const CVec3D &rotationVec, const float rotationAngle ) {
	origin = rotationOrigin;
	vec = rotationVec;
	angle = rotationAngle;
	axisValid = false;
}

SMF_INLINE_FORCED void CRotation::setOrigin( const CVec3D &rotationOrigin ) {
	origin = rotationOrigin;
}

SMF_INLINE_FORCED void CRotation::setVec( const CVec3D &rotationVec ) {
	vec = rotationVec;
	axisValid = false;
}

SMF_INLINE_FORCED void CRotation::setVec( float x, float y, float z ) {
	vec[0] = x;
	vec[1] = y;
	vec[2] = z;
	axisValid = false;
}

SMF_INLINE_FORCED void CRotation::setAngle( const float rotationAngle ) {
	angle = rotationAngle;
	axisValid = false;
}

SMF_INLINE_FORCED void CRotation::scale( const float s ) {
	angle *= s;
	axisValid = false;
}

SMF_INLINE_FORCED void CRotation::reCalculateMatrix() {
	axisValid = false;
	toMat3();
}

SMF_INLINE_FORCED const CVec3D &CRotation::getOrigin() const {
	return origin;
}

SMF_INLINE_FORCED const CVec3D &CRotation::getVec() const  {
	return vec;
}

SMF_INLINE_FORCED float CRotation::getAngle() const  {
	return angle;
}

SMF_INLINE_FORCED CRotation CRotation::operator-() const {
	return CRotation( origin, vec, -angle );
}

SMF_INLINE_FORCED CRotation CRotation::operator*( const float s ) const {
	return CRotation( origin, vec, angle * s );
}

SMF_INLINE_FORCED CRotation CRotation::operator/( const float s ) const {
	SMF_ASSERT( s != 0.0f );
	return CRotation( origin, vec, angle / s );
}

SMF_INLINE_FORCED CRotation &CRotation::operator*=( const float s ) {
	angle *= s;
	axisValid = false;
	return *this;
}

SMF_INLINE_FORCED CRotation &CRotation::operator/=( const float s ) {
	SMF_ASSERT( s != 0.0f );
	angle /= s;
	axisValid = false;
	return *this;
}

SMF_INLINE_FORCED CVec3D CRotation::operator*( const CVec3D &v ) const {
	if ( !axisValid ) {
		toMat3();
	}
	return ((v - origin) * axis + origin);
}

SMF_INLINE_FORCED CRotation operator*( const float s, const CRotation &r ) {
	return r * s;
}

SMF_INLINE_FORCED CVec3D operator*( const CVec3D &v, const CRotation &r ) {
	return r * v;
}

SMF_INLINE_FORCED CVec3D &operator*=( CVec3D &v, const CRotation &r ) {
	v = r * v;
	return v;
}

SMF_INLINE_FORCED void CRotation::rotatePoint( CVec3D &point ) const {
	if ( !axisValid ) {
		toMat3();
	}
	point = ((point - origin) * axis + origin);
}

} /**  end MATH */
} /** end SMF */
#endif /* _SMF__MATH_ROTATION_H__ */
