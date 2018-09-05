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

#include "math/SMF_EulerAngles.h"
#include "math/SMF_Quaternion.h"
#include "math/SMF_Rotation.h"

//#pragma hdrstop

#include <float.h>
namespace SMF {
namespace MATH{
CEulerAngles ang_zero( 0.0f, 0.0f, 0.0f );


/*
=================
CEulerAngles::normalize360

returns angles normalized to the range [0 <= angle < 360]
=================
*/
CEulerAngles& CEulerAngles::normalize360() {
	int i;

	for ( i = 0; i < 3; i++ ) {
		if ( ( (*this)[i] >= 360.0f ) || ( (*this)[i] < 0.0f ) ) {
			(*this)[i] -= CMath::floor( (*this)[i] / 360.0f ) * 360.0f;

			if ( (*this)[i] >= 360.0f ) {
				(*this)[i] -= 360.0f;
			}
			if ( (*this)[i] < 0.0f ) {
				(*this)[i] += 360.0f;
			}
		}
	}

	return *this;
}

/*
=================
CEulerAngles::normalize180

returns angles normalized to the range [-180 < angle <= 180]
=================
*/
CEulerAngles& CEulerAngles::normalize180() {
	normalize360();

	if ( pitch > 180.0f ) {
		pitch -= 360.0f;
	}
	
	if ( yaw > 180.0f ) {
		yaw -= 360.0f;
	}

	if ( roll > 180.0f ) {
		roll -= 360.0f;
	}
	return *this;
}

/*
=================
CEulerAngles::ToVectors
=================
*/
void CEulerAngles::toVectors( CVec3D *forward, CVec3D *right, CVec3D *up ) const {
	float sr, sp, sy, cr, cp, cy;
	
	CMath::sincos( DEG2RAD( yaw ), sy, cy );
	CMath::sincos( DEG2RAD( pitch ), sp, cp );
	CMath::sincos( DEG2RAD( roll ), sr, cr );

	if ( forward ) {
		forward->set( cp * cy, cp * sy, -sp );
	}

	if ( right ) {
		right->set( -sr * sp * cy + cr * sy, -sr * sp * sy + -cr * cy, -sr * cp );
	}

	if ( up ) {
		up->set( cr * sp * cy + -sr * -sy, cr * sp * sy + -sr * cy, cr * cp );
	}
}

/*
=================
CEulerAngles::ToForward
=================
*/
CVec3D CEulerAngles::toForward() const {
	float sp, sy, cp, cy;
	
	CMath::sincos( DEG2RAD( yaw ), sy, cy );
	CMath::sincos( DEG2RAD( pitch ), sp, cp );

	return CVec3D( cp * cy, cp * sy, -sp );
}

/*
=================
CEulerAngles::toQuat
=================
*/
CQuaternion CEulerAngles::toQuat() const {
	float sx, cx, sy, cy, sz, cz;
	float sxcy, cxcy, sxsy, cxsy;

	CMath::sincos( DEG2RAD( yaw ) * 0.5f, sz, cz );
	CMath::sincos( DEG2RAD( pitch ) * 0.5f, sy, cy );
	CMath::sincos( DEG2RAD( roll ) * 0.5f, sx, cx );

	sxcy = sx * cy;
	cxcy = cx * cy;
	sxsy = sx * sy;
	cxsy = cx * sy;

	return CQuaternion( cxsy*sz - sxcy*cz, -cxsy*cz - sxcy*sz, sxsy*cz - cxcy*sz, cxcy*cz + sxsy*sz );
}

/*
=================
CEulerAngles::toRotation
=================
*/
CRotation CEulerAngles::toRotation() const {
	CVec3D vec;
	float angle, w;
	float sx, cx, sy, cy, sz, cz;
	float sxcy, cxcy, sxsy, cxsy;

	if ( pitch == 0.0f ) {
		if ( yaw == 0.0f ) {
			return CRotation( CVec3D::origin, CVec3D( -1.0f, 0.0f, 0.0f ), roll );
		}
		if ( roll == 0.0f ) {
			return CRotation( CVec3D::origin, CVec3D( 0.0f, 0.0f, -1.0f ), yaw );
		}
	} else if ( yaw == 0.0f && roll == 0.0f ) {
		return CRotation( CVec3D::origin, CVec3D( 0.0f, -1.0f, 0.0f ), pitch );
	}

	CMath::sincos( DEG2RAD( yaw ) * 0.5f, sz, cz );
	CMath::sincos( DEG2RAD( pitch ) * 0.5f, sy, cy );
	CMath::sincos( DEG2RAD( roll ) * 0.5f, sx, cx );

	sxcy = sx * cy;
	cxcy = cx * cy;
	sxsy = sx * sy;
	cxsy = cx * sy;

	vec.x =  cxsy * sz - sxcy * cz;
	vec.y = -cxsy * cz - sxcy * sz;
	vec.z =  sxsy * cz - cxcy * sz;
	w =		 cxcy * cz + sxsy * sz;
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
=================
CEulerAngles::toMat3
=================
*/
CMat3D CEulerAngles::toMat3() const {
	CMat3D mat;
	float sr, sp, sy, cr, cp, cy;

	CMath::sincos( DEG2RAD( yaw ), sy, cy );
	CMath::sincos( DEG2RAD( pitch ), sp, cp );
	CMath::sincos( DEG2RAD( roll ), sr, cr );

	mat[ 0 ].set( cp * cy, cp * sy, -sp );
	mat[ 1 ].set( sr * sp * cy + cr * -sy, sr * sp * sy + cr * cy, sr * cp );
	mat[ 2 ].set( cr * sp * cy + -sr * -sy, cr * sp * sy + -sr * cy, cr * cp );

	return mat;
}

/*
=================
CEulerAngles::toMat4
=================
*/
CMat4D CEulerAngles::toMat4() const {
	return toMat3().toMat4();
}

/*
=================
CEulerAngles::toAngularVelocity
=================
*/
CVec3D CEulerAngles::toAngularVelocity() const {
	CRotation rotation = CEulerAngles::toRotation();
	return rotation.getVec() * DEG2RAD( rotation.getAngle() );
}

/*
=============
CEulerAngles::toString
=============
*/
const char *CEulerAngles::toString( int precision ) const {
	return CMyString::floatArraytoString( toFloatPtr(), getDimension(), precision );
}
} //end MATH
} //end SMF
