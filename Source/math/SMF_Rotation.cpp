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

#include "math/SMF_Rotation.h"
#include "math/SMF_Vector.h"
#include "math/SMF_Matriz.h"
#include "math/SMF_EulerAngles.h"
#include "math/SMF_Quaternion.h"
#include "util/SMF_StringUtils.h"
#include "math/all.h"


namespace SMF {
namespace MATH{


CEulerAngles CRotation::toAngles() const {
	return toMat3().toAngles();
}


CQuaternion CRotation::toQuat() const {
	float a, s, c;

	a = angle * ( CMath::M_DEG2RAD * 0.5f );
	CMath::sincos( a, s, c );
	return CQuaternion( vec.x * s, vec.y * s, vec.z * s, c );
}


const CMat3D &CRotation::toMat3() const {
	float wx, wy, wz;
	float xx, yy, yz;
	float xy, xz, zz;
	float x2, y2, z2;
	float a, cos_a, sin_a, x, y, z;

	if ( axisValid ) {
		return axis;
	}

	a = angle * ( CMath::M_DEG2RAD * 0.5f );  //convert to radians and divide by 2
	CMath::sincos( a, sin_a, cos_a );

	   
	 /** vec is a normalized vector to rotate around*/
	x = vec[0] * sin_a;
	y = vec[1] * sin_a;
	z = vec[2] * sin_a;

	x2 = x + x;  // 2*(vec.x * sin(angle) )
	y2 = y + y;  // 2*(vec.y * sin(angle) )
	z2 = z + z;  // 2*(vec.z * sin(angle) )

	xx = x * x2; // (vec.x * sin(angle))*[2*(vec.x * sin(angle) )]
	xy = x * y2; // (vec.x * sin(angle))*[2*(vec.y * sin(angle) )]
	xz = x * z2; // (vec.x * sin(angle))*[2*(vec.z * sin(angle) )]

	yy = y * y2; // (vec.y * sin(angle))*[2*(vec.x * sin(angle) )]
	yz = y * z2; // (vec.y * sin(angle))*[2*(vec.z * sin(angle) )]
	zz = z * z2; // (vec.z * sin(angle))*[2*(vec.z * sin(angle) )]

	wx = cos_a * x2; // cos(angle)*[2*(vec.x * sin(angle) )]
	wy = cos_a * y2; // cos(angle)*[2*(vec.y * sin(angle) )]
	wz = cos_a * z2; // cos(angle)*[2*(vec.z * sin(angle) )]
	// axis [column][row]. remember CMat3D is column major
	axis[ 0 ][ 0 ] = 1.0f - ( yy + zz );
	axis[ 0 ][ 1 ] = xy - wz;
	axis[ 0 ][ 2 ] = xz + wy;

	axis[ 1 ][ 0 ] = xy + wz;
	axis[ 1 ][ 1 ] = 1.0f - ( xx + zz );
	axis[ 1 ][ 2 ] = yz - wx;

	axis[ 2 ][ 0 ] = xz - wy;
	axis[ 2 ][ 1 ] = yz + wx;
	axis[ 2 ][ 2 ] = 1.0f - ( xx + yy );

	axisValid = true;

	return axis;
}

/*
============
CRotation::toMat4
============
*/
CMat4D CRotation::toMat4() const {
	return toMat3().toMat4();
}

/*
============
CRotation::toAngularVelocity
============
*/
CVec3D CRotation::toAngularVelocity() const {
	return vec * DEG2RAD( angle );
}

/*
============
CRotation::normalize180
============
*/
void CRotation::normalize180() {
	angle -= CMath::floor( angle / 360.0f ) * 360.0f;
	if ( angle > 180.0f ) {
		angle -= 360.0f;
	}
	else if ( angle < -180.0f ) {
		angle += 360.0f;
	}
}

/*
============
CRotation::normalize360
============
*/
void CRotation::normalize360() {
	angle -= CMath::floor( angle / 360.0f ) * 360.0f;
	if ( angle > 360.0f ) {
		angle -= 360.0f;
	}
	else if ( angle < 0.0f ) {
		angle += 360.0f;
	}
}
} //END MATH
} //end SMF
