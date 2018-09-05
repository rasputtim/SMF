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

#ifndef _SMF__MATH_PLUECKER_H__
#define _SMF__MATH_PLUECKER_H__
#include "../SMF_Config.h"
#include "SMF_Vector.h"
#include "SMF_Matriz.h"
#include "SMF_EulerAngles.h"
#include "../util/SMF_StringUtils.h"

namespace SMF{
namespace MATH{

/**
 * \class CPluecker
 *
 * \ingroup SMF_Math
 *
 * \brief Implementação do sistema de coordenadas 3d Plucker
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 
 * \see http://en.wikipedia.org/wiki/Plücker_coordinates
 **/
class SMF_API CPluecker {
public:	
					CPluecker();
					explicit CPluecker( const float *a );
					explicit CPluecker( const CVec3D &start, const CVec3D &end );
					explicit CPluecker( const float a1, const float a2, const float a3, const float a4, const float a5, const float a6 );

	float			operator[]( const int index ) const;
	float &			operator[]( const int index );
	CPluecker		operator-() const;											// flips the direction
	CPluecker		operator*( const float a ) const;
	CPluecker		operator/( const float a ) const;
	float			operator*( const CPluecker &a ) const;						// permuted inner product
	CPluecker		operator-( const CPluecker &a ) const;
	CPluecker		operator+( const CPluecker &a ) const;
	CPluecker &	operator*=( const float a );
	CPluecker &	operator/=( const float a );
	CPluecker &	operator+=( const CPluecker &a );
	CPluecker &	operator-=( const CPluecker &a );

	bool			compare( const CPluecker &a ) const;						// exact compare, no epsilon
	bool			compare( const CPluecker &a, const float epsilon ) const;	// compare with epsilon
	bool			operator==(	const CPluecker &a ) const;					// exact compare, no epsilon
	bool			operator!=(	const CPluecker &a ) const;					// exact compare, no epsilon

	void 			set( const float a1, const float a2, const float a3, const float a4, const float a5, const float a6 );
	void			toZero();

	void			fromLine( const CVec3D &start, const CVec3D &end );			// pluecker from line
	void			fromRay( const CVec3D &start, const CVec3D &dir );			// pluecker from ray
	bool			fromPlanes( const CPlane &p1, const CPlane &p2 );			// pluecker from intersection of planes
	bool			toLine( CVec3D &start, CVec3D &end ) const;					// pluecker to line
	bool			toRay( CVec3D &start, CVec3D &dir ) const;					// pluecker to ray
	void			toDir( CVec3D &dir ) const;									// pluecker to direction
	float			permutedInnerProduct( const CPluecker &a ) const;			// pluecker permuted inner product
	float			distance3DSqr( const CPluecker &a ) const;					// pluecker line distance

	float			getLenght() const;										// pluecker length
	float			getLengthSqr() const;									// pluecker squared length
	CPluecker		toNormal() const;									// pluecker normalize
	float			normalizeSelf();										// pluecker normalize

	int				getDimension() const;

	const float *	toFloatPtr() const;
	float *			toFloatPtr();
	const char *	toString( int precision = 2 ) const;

private:
	float			p[6];
};

extern CPluecker pluecker_origin;
#define pluecker_zero pluecker_origin

SMF_INLINE_FORCED CPluecker::CPluecker() {
}

SMF_INLINE_FORCED CPluecker::CPluecker( const float *a ) {
	memcpy( p, a, 6 * sizeof( float ) );
}

SMF_INLINE_FORCED CPluecker::CPluecker( const CVec3D &start, const CVec3D &end ) {
	fromLine( start, end );
}

SMF_INLINE_FORCED CPluecker::CPluecker( const float a1, const float a2, const float a3, const float a4, const float a5, const float a6 ) {
	p[0] = a1;
	p[1] = a2;
	p[2] = a3;
	p[3] = a4;
	p[4] = a5;
	p[5] = a6;
}

SMF_INLINE_FORCED CPluecker CPluecker::operator-() const {
	return CPluecker( -p[0], -p[1], -p[2], -p[3], -p[4], -p[5] );
}

SMF_INLINE_FORCED float CPluecker::operator[]( const int index ) const {
	return p[index];
}

SMF_INLINE_FORCED float &CPluecker::operator[]( const int index ) {
	return p[index];
}

SMF_INLINE_FORCED CPluecker CPluecker::operator*( const float a ) const {
	return CPluecker( p[0]*a, p[1]*a, p[2]*a, p[3]*a, p[4]*a, p[5]*a );
}

SMF_INLINE_FORCED float CPluecker::operator*( const CPluecker &a ) const {
	return p[0] * a.p[4] + p[1] * a.p[5] + p[2] * a.p[3] + p[4] * a.p[0] + p[5] * a.p[1] + p[3] * a.p[2];
}

SMF_INLINE_FORCED CPluecker CPluecker::operator/( const float a ) const {
	float inva;

	SMF_ASSERT( a != 0.0f );
	inva = 1.0f / a;
	return CPluecker( p[0]*inva, p[1]*inva, p[2]*inva, p[3]*inva, p[4]*inva, p[5]*inva );
}

SMF_INLINE_FORCED CPluecker CPluecker::operator+( const CPluecker &a ) const {
	return CPluecker( p[0] + a[0], p[1] + a[1], p[2] + a[2], p[3] + a[3], p[4] + a[4], p[5] + a[5] );
}

SMF_INLINE_FORCED CPluecker CPluecker::operator-( const CPluecker &a ) const {
	return CPluecker( p[0] - a[0], p[1] - a[1], p[2] - a[2], p[3] - a[3], p[4] - a[4], p[5] - a[5] );
}

SMF_INLINE_FORCED CPluecker &CPluecker::operator*=( const float a ) {
	p[0] *= a;
	p[1] *= a;
	p[2] *= a;
	p[3] *= a;
	p[4] *= a;
	p[5] *= a;
	return *this;
}

SMF_INLINE_FORCED CPluecker &CPluecker::operator/=( const float a ) {
	float inva;

	SMF_ASSERT( a != 0.0f );
	inva = 1.0f / a;
	p[0] *= inva;
	p[1] *= inva;
	p[2] *= inva;
	p[3] *= inva;
	p[4] *= inva;
	p[5] *= inva;
	return *this;
}

SMF_INLINE_FORCED CPluecker &CPluecker::operator+=( const CPluecker &a ) {
	p[0] += a[0];
	p[1] += a[1];
	p[2] += a[2];
	p[3] += a[3];
	p[4] += a[4];
	p[5] += a[5];
	return *this;
}

SMF_INLINE_FORCED CPluecker &CPluecker::operator-=( const CPluecker &a ) {
	p[0] -= a[0];
	p[1] -= a[1];
	p[2] -= a[2];
	p[3] -= a[3];
	p[4] -= a[4];
	p[5] -= a[5];
	return *this;
}

SMF_INLINE_FORCED bool CPluecker::compare( const CPluecker &a ) const {
	return ( ( p[0] == a[0] ) && ( p[1] == a[1] ) && ( p[2] == a[2] ) &&
			( p[3] == a[3] ) && ( p[4] == a[4] ) && ( p[5] == a[5] ) );
}

SMF_INLINE_FORCED bool CPluecker::compare( const CPluecker &a, const float epsilon ) const {
	if ( CMath::fabs( p[0] - a[0] ) > epsilon ) {
		return false;
	}
			
	if ( CMath::fabs( p[1] - a[1] ) > epsilon ) {
		return false;
	}

	if ( CMath::fabs( p[2] - a[2] ) > epsilon ) {
		return false;
	}

	if ( CMath::fabs( p[3] - a[3] ) > epsilon ) {
		return false;
	}

	if ( CMath::fabs( p[4] - a[4] ) > epsilon ) {
		return false;
	}

	if ( CMath::fabs( p[5] - a[5] ) > epsilon ) {
		return false;
	}

	return true;
}

SMF_INLINE_FORCED bool CPluecker::operator==( const CPluecker &a ) const {
	return compare( a );
}

SMF_INLINE_FORCED bool CPluecker::operator!=( const CPluecker &a ) const {
	return !compare( a );
}

SMF_INLINE_FORCED void CPluecker::set( const float a1, const float a2, const float a3, const float a4, const float a5, const float a6 ) {
	p[0] = a1;
	p[1] = a2;
	p[2] = a3;
	p[3] = a4;
	p[4] = a5;
	p[5] = a6;
}

SMF_INLINE_FORCED void CPluecker::toZero() {
	p[0] = p[1] = p[2] = p[3] = p[4] = p[5] = 0.0f;
}

SMF_INLINE_FORCED void CPluecker::fromLine( const CVec3D &start, const CVec3D &end ) {
	p[0] = start[0] * end[1] - end[0] * start[1];
	p[1] = start[0] * end[2] - end[0] * start[2];
	p[2] = start[0] - end[0];
	p[3] = start[1] * end[2] - end[1] * start[2];
	p[4] = start[2] - end[2];
	p[5] = end[1] - start[1];
}

SMF_INLINE_FORCED void CPluecker::fromRay( const CVec3D &start, const CVec3D &dir ) {
	p[0] = start[0] * dir[1] - dir[0] * start[1];
	p[1] = start[0] * dir[2] - dir[0] * start[2];
	p[2] = -dir[0];
	p[3] = start[1] * dir[2] - dir[1] * start[2];
	p[4] = -dir[2];
	p[5] = dir[1];
}

SMF_INLINE_FORCED bool CPluecker::toLine( CVec3D &start, CVec3D &end ) const {
	CVec3D dir1, dir2;
	float d;

	dir1[0] = p[3];
	dir1[1] = -p[1];
	dir1[2] = p[0];

	dir2[0] = -p[2];
	dir2[1] = p[5];
	dir2[2] = -p[4];

	d = dir2 * dir2;
	if ( d == 0.0f ) {
		return false; // pluecker coordinate does not represent a line
	}

	start = dir2.cross(dir1) * (1.0f / d);
	end = start + dir2;
	return true;
}

SMF_INLINE_FORCED bool CPluecker::toRay( CVec3D &start, CVec3D &dir ) const {
	CVec3D dir1;
	float d;

	dir1[0] = p[3];
	dir1[1] = -p[1];
	dir1[2] = p[0];

	dir[0] = -p[2];
	dir[1] = p[5];
	dir[2] = -p[4];

	d = dir * dir;
	if ( d == 0.0f ) {
		return false; // pluecker coordinate does not represent a line
	}

	start = dir.cross(dir1) * (1.0f / d);
	return true;
}

SMF_INLINE_FORCED void CPluecker::toDir( CVec3D &dir ) const {
	dir[0] = -p[2];
	dir[1] = p[5];
	dir[2] = -p[4];
}

SMF_INLINE_FORCED float CPluecker::permutedInnerProduct( const CPluecker &a ) const {
	return p[0] * a.p[4] + p[1] * a.p[5] + p[2] * a.p[3] + p[4] * a.p[0] + p[5] * a.p[1] + p[3] * a.p[2];
}

SMF_INLINE_FORCED float CPluecker::getLenght() const {
	return ( float )CMath::sqrt( p[5] * p[5] + p[4] * p[4] + p[2] * p[2] );
}

SMF_INLINE_FORCED float CPluecker::getLengthSqr() const {
	return ( p[5] * p[5] + p[4] * p[4] + p[2] * p[2] );
}

SMF_INLINE_FORCED float CPluecker::normalizeSelf() {
	float l, d;

	l = getLengthSqr();
	if ( l == 0.0f ) {
		return l; // pluecker coordinate does not represent a line
	}
	d = CMath::invSqrt( l );
	p[0] *= d;
	p[1] *= d;
	p[2] *= d;
	p[3] *= d;
	p[4] *= d;
	p[5] *= d;
	return d * l;
}

SMF_INLINE_FORCED CPluecker CPluecker::toNormal() const {
	float d;

	d = getLengthSqr();
	if ( d == 0.0f ) {
		return *this; // pluecker coordinate does not represent a line
	}
	d = CMath::invSqrt( d );
	return CPluecker( p[0]*d, p[1]*d, p[2]*d, p[3]*d, p[4]*d, p[5]*d );
}

SMF_INLINE_FORCED int CPluecker::getDimension() const {
	return 6;
}

SMF_INLINE_FORCED const float *CPluecker::toFloatPtr() const {
	return p;
}

SMF_INLINE_FORCED float *CPluecker::toFloatPtr() {
	return p;
}


} //end MATH
} //end SMF
#endif /* !__MATH_PLUECKER_H__ */
