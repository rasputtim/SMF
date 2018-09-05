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

//#pragma hdrstop

#include "math/SMF_SimdGeneric.h"
#include "math/SMF_Vector.h"
#include "math/SMF_JointTransform.h"

#include "geometry/SMF_Plane.h"
#include "math/SMF_Math.h"
#include "geometry/SMF_DrawVert.h"
#include "math/SMF_Math.h"
#include "geometry/SMF_AABBox.h"

namespace SMF{
namespace MATH{
//===============================================================
//
//	Generic implementation of CSIMDProcessor
//
//===============================================================

#define UNROLL1(Y) { int _IX; for (_IX=0;_IX<count;_IX++) {Y(_IX);} }
#define UNROLL2(Y) { int _IX, _NM = count&0xfffffffe; for (_IX=0;_IX<_NM;_IX+=2){Y(_IX+0);Y(_IX+1);} if (_IX < count) {Y(_IX);}}
#define UNROLL4(Y) { int _IX, _NM = count&0xfffffffc; for (_IX=0;_IX<_NM;_IX+=4){Y(_IX+0);Y(_IX+1);Y(_IX+2);Y(_IX+3);}for(;_IX<count;_IX++){Y(_IX);}}
#define UNROLL8(Y) { int _IX, _NM = count&0xfffffff8; for (_IX=0;_IX<_NM;_IX+=8){Y(_IX+0);Y(_IX+1);Y(_IX+2);Y(_IX+3);Y(_IX+4);Y(_IX+5);Y(_IX+6);Y(_IX+7);} _NM = count&0xfffffffe; for(;_IX<_NM;_IX+=2){Y(_IX); Y(_IX+1);} if (_IX < count) {Y(_IX);} }

#ifdef _DEBUG
#define NODEFAULT	default: SMF_ASSERT( 0 )
#elif _WIN32
#define NODEFAULT	default: __assume( 0 )
#else
#define NODEFAULT
#endif


/*
============
CSIMD_Generic::getName
============
*/
const char * CSIMD_Generic::getName() const {
	return "generic code";
}

/*
============
CSIMD_Generic::add

  dst[i] = constant + src[i];
============
*/
void VPCALL CSIMD_Generic::add( float *dst, const float constant, const float *src, const int count ) {
#define OPER(X) dst[(X)] = src[(X)] + constant;
	UNROLL4(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::add

  dst[i] = src0[i] + src1[i];
============
*/
void VPCALL CSIMD_Generic::add( float *dst, const float *src0, const float *src1, const int count ) {
#define OPER(X) dst[(X)] = src0[(X)] + src1[(X)];
	UNROLL4(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::sub

  dst[i] = constant - src[i];
============
*/
void VPCALL CSIMD_Generic::sub( float *dst, const float constant, const float *src, const int count ) {
	double c = constant;
#define OPER(X) dst[(X)] = c - src[(X)];
	UNROLL4(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::sub

  dst[i] = src0[i] - src1[i];
============
*/
void VPCALL CSIMD_Generic::sub( float *dst, const float *src0, const float *src1, const int count ) {
#define OPER(X) dst[(X)] = src0[(X)] - src1[(X)];
	UNROLL4(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::mul

  dst[i] = constant * src[i];
============
*/
void VPCALL CSIMD_Generic::mul( float *dst, const float constant, const float *src0, const int count) {
	double c = constant;
#define OPER(X) (dst[(X)] = (c * src0[(X)]))
	UNROLL4(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::mul

  dst[i] = src0[i] * src1[i];
============
*/
void VPCALL CSIMD_Generic::mul( float *dst, const float *src0, const float *src1, const int count ) {
#define OPER(X) (dst[(X)] = src0[(X)] * src1[(X)])
	UNROLL4(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::div

  dst[i] = constant / divisor[i];
============
*/
void VPCALL CSIMD_Generic::div( float *dst, const float constant, const float *divisor, const int count ) {
	double c = constant;
#define OPER(X) (dst[(X)] = (c / divisor[(X)]))
	UNROLL4(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::div

  dst[i] = src0[i] / src1[i];
============
*/
void VPCALL CSIMD_Generic::div( float *dst, const float *src0, const float *src1, const int count ) {
#define OPER(X) (dst[(X)] = src0[(X)] / src1[(X)])
	UNROLL4(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::mulAdd

  dst[i] += constant * src[i];
============
*/
void VPCALL CSIMD_Generic::mulAdd( float *dst, const float constant, const float *src, const int count ) {
	double c = constant;
#define OPER(X) (dst[(X)] += c * src[(X)])
	UNROLL4(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::mulAdd

  dst[i] += src0[i] * src1[i];
============
*/
void VPCALL CSIMD_Generic::mulAdd( float *dst, const float *src0, const float *src1, const int count ) {
#define OPER(X) (dst[(X)] += src0[(X)] * src1[(X)])
	UNROLL4(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::mulSub

  dst[i] -= constant * src[i];
============
*/
void VPCALL CSIMD_Generic::mulSub( float *dst, const float constant, const float *src, const int count ) {
	double c = constant;
#define OPER(X) (dst[(X)] -= c * src[(X)])
	UNROLL4(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::mulSub

  dst[i] -= src0[i] * src1[i];
============
*/
void VPCALL CSIMD_Generic::mulSub( float *dst, const float *src0, const float *src1, const int count ) {
#define OPER(X) (dst[(X)] -= src0[(X)] * src1[(X)])
	UNROLL4(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::dot

  dst[i] = constant * src[i];
============
*/
void VPCALL CSIMD_Generic::dot( float *dst, const CVec3D &constant, const CVec3D *src, const int count ) {
#define OPER(X) dst[(X)] = constant * src[(X)];
	UNROLL1(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::dot

  dst[i] = constant * src[i].Normal() + src[i][3];
============
*/
void VPCALL CSIMD_Generic::dot( float *dst, const CVec3D &constant, const CPlane *src, const int count ) {
#define OPER(X) dst[(X)] = constant * src[(X)].getNormal() + src[(X)][3];
	UNROLL1(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::dot

  dst[i] = constant * src[i].xyz;
============
*/
void VPCALL CSIMD_Generic::dot( float *dst, const CVec3D &constant, const CVertex *src, const int count ) {
#define OPER(X) dst[(X)] = constant * src[(X)].xyz;
	UNROLL1(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::dot

  dst[i] = constant.Normal() * src[i] + constant[3];
============
*/
void VPCALL CSIMD_Generic::dot( float *dst, const CPlane &constant, const CVec3D *src, const int count ) {
#define OPER(X) dst[(X)] = constant.getNormal() * src[(X)] + constant[3];
	UNROLL1(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::dot

  dst[i] = constant.Normal() * src[i].Normal() + constant[3] * src[i][3];
============
*/
void VPCALL CSIMD_Generic::dot( float *dst, const CPlane &constant, const CPlane *src, const int count ) {
#define OPER(X) dst[(X)] = constant.getNormal() * src[(X)].getNormal() + constant[3] * src[(X)][3];
	UNROLL1(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::dot

  dst[i] = constant.getNormal() * src[i].xyz + constant[3];
============
*/
void VPCALL CSIMD_Generic::dot( float *dst, const CPlane &constant, const CVertex *src, const int count ) {
#define OPER(X) dst[(X)] = constant.getNormal() * src[(X)].xyz + constant[3];
	UNROLL1(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::dot

  dst[i] = src0[i] * src1[i];
============
*/
void VPCALL CSIMD_Generic::dot( float *dst, const CVec3D *src0, const CVec3D *src1, const int count ) {
#define OPER(X) dst[(X)] = src0[(X)] * src1[(X)];
	UNROLL1(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::dot

  dot = src1[0] * src2[0] + src1[1] * src2[1] + src1[2] * src2[2] + ...
============
*/
void VPCALL CSIMD_Generic::dot( float &dot, const float *src1, const float *src2, const int count ) {
#if 1

	switch( count ) {
		case 0: {
			dot = 0.0f;
			return;
		}
		case 1: {
			dot = src1[0] * src2[0];
			return;
		}
		case 2: {
			dot = src1[0] * src2[0] + src1[1] * src2[1];
			return;
		}
		case 3: {
			dot = src1[0] * src2[0] + src1[1] * src2[1] + src1[2] * src2[2];
			return;
		}
		default: {
			int i;
			double s0, s1, s2, s3;
			s0 = src1[0] * src2[0];
			s1 = src1[1] * src2[1];
			s2 = src1[2] * src2[2];
			s3 = src1[3] * src2[3];
			for ( i = 4; i < count-7; i += 8 ) {
				s0 += src1[i+0] * src2[i+0];
				s1 += src1[i+1] * src2[i+1];
				s2 += src1[i+2] * src2[i+2];
				s3 += src1[i+3] * src2[i+3];
				s0 += src1[i+4] * src2[i+4];
				s1 += src1[i+5] * src2[i+5];
				s2 += src1[i+6] * src2[i+6];
				s3 += src1[i+7] * src2[i+7];
			}
			switch( count - i ) {
				NODEFAULT;
				case 7: s0 += src1[i+6] * src2[i+6];
				case 6: s1 += src1[i+5] * src2[i+5];
				case 5: s2 += src1[i+4] * src2[i+4];
				case 4: s3 += src1[i+3] * src2[i+3];
				case 3: s0 += src1[i+2] * src2[i+2];
				case 2: s1 += src1[i+1] * src2[i+1];
				case 1: s2 += src1[i+0] * src2[i+0];
				case 0: break;
			}
			double sum;
			sum = s3;
			sum += s2;
			sum += s1;
			sum += s0;
			dot = sum;
		}
	}

#else

	dot = 0.0f;
	for ( i = 0; i < count; i++ ) {
		dot += src1[i] * src2[i];
	}

#endif
}

/*
============
CSIMD_Generic::cmpGT

  dst[i] = src0[i] > constant;
============
*/
void VPCALL CSIMD_Generic::cmpGT( sf_u8 *dst, const float *src0, const float constant, const int count ) {
#define OPER(X) dst[(X)] = src0[(X)] > constant;
	UNROLL4(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::cmpGT

  dst[i] |= ( src0[i] > constant ) << bitNum;
============
*/
void VPCALL CSIMD_Generic::cmpGT( sf_u8 *dst, const sf_u8 bitNum, const float *src0, const float constant, const int count ) {
#define OPER(X) dst[(X)] |= ( src0[(X)] > constant ) << bitNum;
	UNROLL4(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::cmpGE

  dst[i] = src0[i] >= constant;
============
*/
void VPCALL CSIMD_Generic::cmpGE( sf_u8 *dst, const float *src0, const float constant, const int count ) {
#define OPER(X) dst[(X)] = src0[(X)] >= constant;
	UNROLL4(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::cmpGE

  dst[i] |= ( src0[i] >= constant ) << bitNum;
============
*/
void VPCALL CSIMD_Generic::cmpGE( sf_u8 *dst, const sf_u8 bitNum, const float *src0, const float constant, const int count ) {
#define OPER(X) dst[(X)] |= ( src0[(X)] >= constant ) << bitNum;
	UNROLL4(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::cmpLT

  dst[i] = src0[i] < constant;
============
*/
void VPCALL CSIMD_Generic::cmpLT( sf_u8 *dst, const float *src0, const float constant, const int count ) {
#define OPER(X) dst[(X)] = src0[(X)] < constant;
	UNROLL4(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::cmpLT

  dst[i] |= ( src0[i] < constant ) << bitNum;
============
*/
void VPCALL CSIMD_Generic::cmpLT( sf_u8 *dst, const sf_u8 bitNum, const float *src0, const float constant, const int count ) {
#define OPER(X) dst[(X)] |= ( src0[(X)] < constant ) << bitNum;
	UNROLL4(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::cmpLE

  dst[i] = src0[i] <= constant;
============
*/
void VPCALL CSIMD_Generic::cmpLE( sf_u8 *dst, const float *src0, const float constant, const int count ) {
#define OPER(X) dst[(X)] = src0[(X)] <= constant;
	UNROLL4(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::cmpLE

  dst[i] |= ( src0[i] <= constant ) << bitNum;
============
*/
void VPCALL CSIMD_Generic::cmpLE( sf_u8 *dst, const sf_u8 bitNum, const float *src0, const float constant, const int count ) {
#define OPER(X) dst[(X)] |= ( src0[(X)] <= constant ) << bitNum;
	UNROLL4(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::minMax
============
*/
void VPCALL CSIMD_Generic::minMax( float &min, float &max, const float *src, const int count ) {
	min = CMath::INFINITY_FLOAT; max = -CMath::INFINITY_FLOAT;
#define OPER(X) if ( src[(X)] < min ) {min = src[(X)];} if ( src[(X)] > max ) {max = src[(X)];}
	UNROLL1(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::minMax
============
*/
void VPCALL CSIMD_Generic::minMax( CVec2D &min, CVec2D &max, const CVec2D *src, const int count ) {
	min[0] = min[1] = CMath::INFINITY_FLOAT; max[0] = max[1] = -CMath::INFINITY_FLOAT;
#define OPER(X) const CVec2D &v = src[(X)]; if ( v[0] < min[0] ) { min[0] = v[0]; } if ( v[0] > max[0] ) { max[0] = v[0]; } if ( v[1] < min[1] ) { min[1] = v[1]; } if ( v[1] > max[1] ) { max[1] = v[1]; }
	UNROLL1(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::minMax
============
*/
void VPCALL CSIMD_Generic::minMax( CVec3D &min, CVec3D &max, const CVec3D *src, const int count ) {
	min[0] = min[1] = min[2] = CMath::INFINITY_FLOAT; max[0] = max[1] = max[2] = -CMath::INFINITY_FLOAT;
#define OPER(X) const CVec3D &v = src[(X)]; if ( v[0] < min[0] ) { min[0] = v[0]; } if ( v[0] > max[0] ) { max[0] = v[0]; } if ( v[1] < min[1] ) { min[1] = v[1]; } if ( v[1] > max[1] ) { max[1] = v[1]; } if ( v[2] < min[2] ) { min[2] = v[2]; } if ( v[2] > max[2] ) { max[2] = v[2]; }
	UNROLL1(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::minMax
============
*/
void VPCALL CSIMD_Generic::minMax( CVec3D &min, CVec3D &max, const CVertex *src, const int count ) {
	min[0] = min[1] = min[2] = CMath::INFINITY_FLOAT; max[0] = max[1] = max[2] = -CMath::INFINITY_FLOAT;
#define OPER(X) const CVec3D &v = src[(X)].xyz; if ( v[0] < min[0] ) { min[0] = v[0]; } if ( v[0] > max[0] ) { max[0] = v[0]; } if ( v[1] < min[1] ) { min[1] = v[1]; } if ( v[1] > max[1] ) { max[1] = v[1]; } if ( v[2] < min[2] ) { min[2] = v[2]; } if ( v[2] > max[2] ) { max[2] = v[2]; }
	UNROLL1(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::minMax
============
*/
void VPCALL CSIMD_Generic::minMax( CVec3D &min, CVec3D &max, const CVertex *src, const int *indexes, const int count ) {
	min[0] = min[1] = min[2] = CMath::INFINITY_FLOAT; max[0] = max[1] = max[2] = -CMath::INFINITY_FLOAT;
#define OPER(X) const CVec3D &v = src[indexes[(X)]].xyz; if ( v[0] < min[0] ) { min[0] = v[0]; } if ( v[0] > max[0] ) { max[0] = v[0]; } if ( v[1] < min[1] ) { min[1] = v[1]; } if ( v[1] > max[1] ) { max[1] = v[1]; } if ( v[2] < min[2] ) { min[2] = v[2]; } if ( v[2] > max[2] ) { max[2] = v[2]; }
	UNROLL1(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::clamp
============
*/
void VPCALL CSIMD_Generic::clamp( float *dst, const float *src, const float min, const float max, const int count ) {
#define OPER(X) dst[(X)] = src[(X)] < min ? min : src[(X)] > max ? max : src[(X)];
	UNROLL1(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::clampMin
============
*/
void VPCALL CSIMD_Generic::clampMin( float *dst, const float *src, const float min, const int count ) {
#define OPER(X) dst[(X)] = src[(X)] < min ? min : src[(X)];
	UNROLL1(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::clampMax
============
*/
void VPCALL CSIMD_Generic::clampMax( float *dst, const float *src, const float max, const int count ) {
#define OPER(X) dst[(X)] = src[(X)] > max ? max : src[(X)];
	UNROLL1(OPER)
#undef OPER
}

/*
================
CSIMD_Generic::memCopy
================
*/
void VPCALL CSIMD_Generic::memCopy( void *dst, const void *src, const int count ) {
	memcpy( dst, src, count );
}

/*
================
CSIMD_Generic::memSet
================
*/
void VPCALL CSIMD_Generic::memSet( void *dst, const int val, const int count ) {
	memset( dst, val, count );
}

/*
============
CSIMD_Generic::zero16
============
*/
void VPCALL CSIMD_Generic::zero16( float *dst, const int count ) {
	memset( dst, 0, count * sizeof( float ) );
}

/*
============
CSIMD_Generic::negate16
============
*/
void VPCALL CSIMD_Generic::negate16( float *dst, const int count ) {
	unsigned int *ptr = reinterpret_cast<unsigned int *>(dst);
#define OPER(X) ptr[(X)] ^= ( 1 << 31 )		// IEEE 32 bits float sign bit
	UNROLL1(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::copy16
============
*/
void VPCALL CSIMD_Generic::copy16( float *dst, const float *src, const int count ) {
#define OPER(X) dst[(X)] = src[(X)]
	UNROLL1(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::add16
============
*/
void VPCALL CSIMD_Generic::add16( float *dst, const float *src1, const float *src2, const int count ) {
#define OPER(X) dst[(X)] = src1[(X)] + src2[(X)]
	UNROLL1(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::sub16
============
*/
void VPCALL CSIMD_Generic::sub16( float *dst, const float *src1, const float *src2, const int count ) {
#define OPER(X) dst[(X)] = src1[(X)] - src2[(X)]
	UNROLL1(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::mul16
============
*/
void VPCALL CSIMD_Generic::mul16( float *dst, const float *src1, const float constant, const int count ) {
#define OPER(X) dst[(X)] = src1[(X)] * constant
	UNROLL1(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::addAssign16
============
*/
void VPCALL CSIMD_Generic::addAssign16( float *dst, const float *src, const int count ) {
#define OPER(X) dst[(X)] += src[(X)]
	UNROLL1(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::subAssign16
============
*/
void VPCALL CSIMD_Generic::subAssign16( float *dst, const float *src, const int count ) {
#define OPER(X) dst[(X)] -= src[(X)]
	UNROLL1(OPER)
#undef OPER
}

/*
============
CSIMD_Generic::mulAssign16
============
*/
void VPCALL CSIMD_Generic::mulAssign16( float *dst, const float constant, const int count ) {
#define OPER(X) dst[(X)] *= constant
	UNROLL1(OPER)
#undef OPER
}


bool VPCALL CSIMD_Generic::matX_inverse_4x4(float* mat)
{

#if 0
float d, di;
di = mat[0];
mat[0] = d = 1.0f / di;  //mat 11
mat[4] *= -d;  //mat 21
mat[8] *= -d;  //mat 31
mat[12] *= -d;  //mat 31
mat[1] *= d;    //mat 12
mat[2] *= d;    //mat 13
mat[3] *= d;    //mat 14
mat[5] += mat[4] * mat[1] * di;
mat[6] += mat[4] * mat[2] * di;
mat[7] += mat[4] * mat[3] * di;
mat[9] += mat[8] * mat[1] * di;
mat[10] += mat[8] * mat[2] * di;
mat[11] += mat[8] * mat[3] * di;
mat[13] += mat[12] * mat[1] * di;
mat[14] += mat[12] * mat[2] * di;
mat[15] += mat[12] * mat[3] * di;
di = mat[5];
mat[5] = d = 1.0f / di;
mat[1] *= -d;
mat[9] *= -d;
mat[13] *= -d;
mat[4] *= d;
mat[6] *= d;
mat[7] *= d;
mat[0] += mat[1] * mat[4] * di;
mat[2] += mat[1] * mat[6] * di;
mat[3] += mat[1] * mat[7] * di;
mat[8] += mat[9] * mat[4] * di;
mat[10] += mat[9] * mat[6] * di;
mat[11] += mat[9] * mat[7] * di;
mat[12] += mat[13] * mat[4] * di;
mat[14] += mat[13] * mat[6] * di;
mat[15] += mat[13] * mat[7] * di;
di = mat[10];
mat[10] = d = 1.0f / di;
mat[2] *= -d;
mat[6] *= -d;
mat[14] *= -d;
mat[8] *= d;
mat[9] *= d;
mat[11] *= d;
mat[0] += mat[2] * mat[8] * di;
mat[1] += mat[2] * mat[9] * di;
mat[3] += mat[2] * mat[11] * di;
mat[4] += mat[6] * mat[8] * di;
mat[5] += mat[6] * mat[9] * di;
mat[7] += mat[6] * mat[11] * di;
mat[12] += mat[14] * mat[8] * di;
mat[13] += mat[14] * mat[9] * di;
mat[15] += mat[14] * mat[11] * di;
di = mat[15];
mat[15] = d = 1.0f / di;
mat[3] *= -d;
mat[7] *= -d;
mat[11] *= -d;
mat[12] *= d;
mat[13] *= d;
mat[14] *= d;
mat[0] += mat[3] * mat[12] * di;
mat[1] += mat[3] * mat[13] * di;
mat[2] += mat[3] * mat[14] * di;
mat[4] += mat[7] * mat[12] * di;
mat[5] += mat[7] * mat[13] * di;
mat[6] += mat[7] * mat[14] * di;
mat[8] += mat[11] * mat[12] * di;
mat[9] += mat[11] * mat[13] * di;
mat[10] += mat[11] * mat[14] * di;
#else
//	6*8+2*6 = 60 multiplications
	//		2*1 =  2 divisions
	CMat2D r0, r1, r2, r3;
	float a, det, invDet;
	//float *mat = reinterpret_cast<float *>(this);

	// r0 = m0.Inverse();
	det = mat[0*4+0] * mat[1*4+1] - mat[0*4+1] * mat[1*4+0];

	if ( CMath::fabs( det ) < MATRIX_INVERSE_EPSILON ) {
		return false;
	}

	invDet = 1.0f / det;

	r0[0][0] =   mat[1*4+1] * invDet;
	r0[0][1] = - mat[0*4+1] * invDet;
	r0[1][0] = - mat[1*4+0] * invDet;
	r0[1][1] =   mat[0*4+0] * invDet;

	// r1 = r0 * m1;
	r1[0][0] = r0[0][0] * mat[0*4+2] + r0[0][1] * mat[1*4+2];
	r1[0][1] = r0[0][0] * mat[0*4+3] + r0[0][1] * mat[1*4+3];
	r1[1][0] = r0[1][0] * mat[0*4+2] + r0[1][1] * mat[1*4+2];
	r1[1][1] = r0[1][0] * mat[0*4+3] + r0[1][1] * mat[1*4+3];

	// r2 = m2 * r1;
	r2[0][0] = mat[2*4+0] * r1[0][0] + mat[2*4+1] * r1[1][0];
	r2[0][1] = mat[2*4+0] * r1[0][1] + mat[2*4+1] * r1[1][1];
	r2[1][0] = mat[3*4+0] * r1[0][0] + mat[3*4+1] * r1[1][0];
	r2[1][1] = mat[3*4+0] * r1[0][1] + mat[3*4+1] * r1[1][1];

	// r3 = r2 - m3;
	r3[0][0] = r2[0][0] - mat[2*4+2];
	r3[0][1] = r2[0][1] - mat[2*4+3];
	r3[1][0] = r2[1][0] - mat[3*4+2];
	r3[1][1] = r2[1][1] - mat[3*4+3];

	// r3.InverseSelf();
	det = r3[0][0] * r3[1][1] - r3[0][1] * r3[1][0];

	if ( CMath::fabs( det ) < MATRIX_INVERSE_EPSILON ) {
		return false;
	}

	invDet = 1.0f / det;

	a = r3[0][0];
	r3[0][0] =   r3[1][1] * invDet;
	r3[0][1] = - r3[0][1] * invDet;
	r3[1][0] = - r3[1][0] * invDet;
	r3[1][1] =   a * invDet;

	// r2 = m2 * r0;
	r2[0][0] = mat[2*4+0] * r0[0][0] + mat[2*4+1] * r0[1][0];
	r2[0][1] = mat[2*4+0] * r0[0][1] + mat[2*4+1] * r0[1][1];
	r2[1][0] = mat[3*4+0] * r0[0][0] + mat[3*4+1] * r0[1][0];
	r2[1][1] = mat[3*4+0] * r0[0][1] + mat[3*4+1] * r0[1][1];

	// m2 = r3 * r2;
	mat[2*4+0] = r3[0][0] * r2[0][0] + r3[0][1] * r2[1][0];
	mat[2*4+1] = r3[0][0] * r2[0][1] + r3[0][1] * r2[1][1];
	mat[3*4+0] = r3[1][0] * r2[0][0] + r3[1][1] * r2[1][0];
	mat[3*4+1] = r3[1][0] * r2[0][1] + r3[1][1] * r2[1][1];

	// m0 = r0 - r1 * m2;
	mat[0*4+0] = r0[0][0] - r1[0][0] * mat[2*4+0] - r1[0][1] * mat[3*4+0];
	mat[0*4+1] = r0[0][1] - r1[0][0] * mat[2*4+1] - r1[0][1] * mat[3*4+1];
	mat[1*4+0] = r0[1][0] - r1[1][0] * mat[2*4+0] - r1[1][1] * mat[3*4+0];
	mat[1*4+1] = r0[1][1] - r1[1][0] * mat[2*4+1] - r1[1][1] * mat[3*4+1];

	// m1 = r1 * r3;
	mat[0*4+2] = r1[0][0] * r3[0][0] + r1[0][1] * r3[1][0];
	mat[0*4+3] = r1[0][0] * r3[0][1] + r1[0][1] * r3[1][1];
	mat[1*4+2] = r1[1][0] * r3[0][0] + r1[1][1] * r3[1][0];
	mat[1*4+3] = r1[1][0] * r3[0][1] + r1[1][1] * r3[1][1];

	// m3 = -r3;
	mat[2*4+2] = -r3[0][0];
	mat[2*4+3] = -r3[0][1];
	mat[3*4+2] = -r3[1][0];
	mat[3*4+3] = -r3[1][1];

	return true;

#endif

}



/*
============
CSIMD_Generic::matX_MultiplyVecX
============
*/
void VPCALL CSIMD_Generic::matX_MultiplyVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec ) {
	int i, j, numRows;
	const float *mPtr, *vPtr;
	float *dstPtr;

	SMF_ASSERT( vec.getSize() >= mat.getNumColumns() );
	SMF_ASSERT( dst.getSize() >= mat.getNumRows() );

	mPtr = mat.toFloatPtr();
	vPtr = vec.toFloatPtr();
	dstPtr = dst.toFloatPtr();
	numRows = mat.getNumRows();
	switch( mat.getNumColumns() ) {
		case 1:
			for ( i = 0; i < numRows; i++ ) {
				dstPtr[i] = mPtr[0] * vPtr[0];
				mPtr++;
			}
			break;
		case 2:
			for ( i = 0; i < numRows; i++ ) {
				dstPtr[i] = mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1];
				mPtr += 2;
			}
			break;
		case 3:
			for ( i = 0; i < numRows; i++ ) {
				dstPtr[i] = mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1] + mPtr[2] * vPtr[2];
				mPtr += 3;
			}
			break;
		case 4:
			for ( i = 0; i < numRows; i++ ) {
				dstPtr[i] = mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1] + mPtr[2] * vPtr[2] +
							mPtr[3] * vPtr[3];
				mPtr += 4;
			}
			break;
		case 5:
			for ( i = 0; i < numRows; i++ ) {
				dstPtr[i] = mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1] + mPtr[2] * vPtr[2] +
							mPtr[3] * vPtr[3] + mPtr[4] * vPtr[4];
				mPtr += 5;
			}
			break;
		case 6:
			for ( i = 0; i < numRows; i++ ) {
				dstPtr[i] = mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1] + mPtr[2] * vPtr[2] +
							mPtr[3] * vPtr[3] + mPtr[4] * vPtr[4] + mPtr[5] * vPtr[5];
				mPtr += 6;
			}
			break;
		default:
			int numColumns = mat.getNumColumns();
			for ( i = 0; i < numRows; i++ ) {
				float sum = mPtr[0] * vPtr[0];
				for ( j = 1; j < numColumns; j++ ) {
					sum += mPtr[j] * vPtr[j];
				}
				dstPtr[i] = sum;
				mPtr += numColumns;
			}
			break;
	}
}

/*
============
CSIMD_Generic::matX_MultiplyAddVecX
============
*/
void VPCALL CSIMD_Generic::matX_MultiplyAddVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec ) {
	int i, j, numRows;
	const float *mPtr, *vPtr;
	float *dstPtr;

	SMF_ASSERT( vec.getSize() >= mat.getNumColumns() );
	SMF_ASSERT( dst.getSize() >= mat.getNumRows() );

	mPtr = mat.toFloatPtr();
	vPtr = vec.toFloatPtr();
	dstPtr = dst.toFloatPtr();
	numRows = mat.getNumRows();
	switch( mat.getNumColumns() ) {
		case 1:
			for ( i = 0; i < numRows; i++ ) {
				dstPtr[i] += mPtr[0] * vPtr[0];
				mPtr++;
			}
			break;
		case 2:
			for ( i = 0; i < numRows; i++ ) {
				dstPtr[i] += mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1];
				mPtr += 2;
			}
			break;
		case 3:
			for ( i = 0; i < numRows; i++ ) {
				dstPtr[i] += mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1] + mPtr[2] * vPtr[2];
				mPtr += 3;
			}
			break;
		case 4:
			for ( i = 0; i < numRows; i++ ) {
				dstPtr[i] += mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1] + mPtr[2] * vPtr[2] +
							mPtr[3] * vPtr[3];
				mPtr += 4;
			}
			break;
		case 5:
			for ( i = 0; i < numRows; i++ ) {
				dstPtr[i] += mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1] + mPtr[2] * vPtr[2] +
							mPtr[3] * vPtr[3] + mPtr[4] * vPtr[4];
				mPtr += 5;
			}
			break;
		case 6:
			for ( i = 0; i < numRows; i++ ) {
				dstPtr[i] += mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1] + mPtr[2] * vPtr[2] +
							mPtr[3] * vPtr[3] + mPtr[4] * vPtr[4] + mPtr[5] * vPtr[5];
				mPtr += 6;
			}
			break;
		default:
			int numColumns = mat.getNumColumns();
			for ( i = 0; i < numRows; i++ ) {
				float sum = mPtr[0] * vPtr[0];
				for ( j = 1; j < numColumns; j++ ) {
					sum += mPtr[j] * vPtr[j];
				}
				dstPtr[i] += sum;
				mPtr += numColumns;
			}
			break;
	}
}

/*
============
CSIMD_Generic::matX_MultiplySubVecX
============
*/
void VPCALL CSIMD_Generic::matX_MultiplySubVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec ) {
	int i, j, numRows;
	const float *mPtr, *vPtr;
	float *dstPtr;

	SMF_ASSERT( vec.getSize() >= mat.getNumColumns() );
	SMF_ASSERT( dst.getSize() >= mat.getNumRows() );

	mPtr = mat.toFloatPtr();
	vPtr = vec.toFloatPtr();
	dstPtr = dst.toFloatPtr();
	numRows = mat.getNumRows();
	switch( mat.getNumColumns() ) {
		case 1:
			for ( i = 0; i < numRows; i++ ) {
				dstPtr[i] -= mPtr[0] * vPtr[0];
				mPtr++;
			}
			break;
		case 2:
			for ( i = 0; i < numRows; i++ ) {
				dstPtr[i] -= mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1];
				mPtr += 2;
			}
			break;
		case 3:
			for ( i = 0; i < numRows; i++ ) {
				dstPtr[i] -= mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1] + mPtr[2] * vPtr[2];
				mPtr += 3;
			}
			break;
		case 4:
			for ( i = 0; i < numRows; i++ ) {
				dstPtr[i] -= mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1] + mPtr[2] * vPtr[2] +
							mPtr[3] * vPtr[3];
				mPtr += 4;
			}
			break;
		case 5:
			for ( i = 0; i < numRows; i++ ) {
				dstPtr[i] -= mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1] + mPtr[2] * vPtr[2] +
							mPtr[3] * vPtr[3] + mPtr[4] * vPtr[4];
				mPtr += 5;
			}
			break;
		case 6:
			for ( i = 0; i < numRows; i++ ) {
				dstPtr[i] -= mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1] + mPtr[2] * vPtr[2] +
							mPtr[3] * vPtr[3] + mPtr[4] * vPtr[4] + mPtr[5] * vPtr[5];
				mPtr += 6;
			}
			break;
		default:
			int numColumns = mat.getNumColumns();
			for ( i = 0; i < numRows; i++ ) {
				float sum = mPtr[0] * vPtr[0];
				for ( j = 1; j < numColumns; j++ ) {
					sum += mPtr[j] * vPtr[j];
				}
				dstPtr[i] -= sum;
				mPtr += numColumns;
			}
			break;
	}
}

/*
============
CSIMD_Generic::matX_TransposeMultiplyVecX
============
*/
void VPCALL CSIMD_Generic::matX_TransposeMultiplyVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec ) {
	int i, j, numColumns;
	const float *mPtr, *vPtr;
	float *dstPtr;

	SMF_ASSERT( vec.getSize() >= mat.getNumRows() );
	SMF_ASSERT( dst.getSize() >= mat.getNumColumns() );

	mPtr = mat.toFloatPtr();
	vPtr = vec.toFloatPtr();
	dstPtr = dst.toFloatPtr();
	numColumns = mat.getNumColumns();
	switch( mat.getNumRows() ) {
		case 1:
			for ( i = 0; i < numColumns; i++ ) {
				dstPtr[i] = *(mPtr) * vPtr[0];
				mPtr++;
			}
			break;
		case 2:
			for ( i = 0; i < numColumns; i++ ) {
				dstPtr[i] = *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1];
				mPtr++;
			}
			break;
		case 3:
			for ( i = 0; i < numColumns; i++ ) {
				dstPtr[i] = *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1] + *(mPtr+2*numColumns) * vPtr[2];
				mPtr++;
			}
			break;
		case 4:
			for ( i = 0; i < numColumns; i++ ) {
				dstPtr[i] = *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1] + *(mPtr+2*numColumns) * vPtr[2] +
						*(mPtr+3*numColumns) * vPtr[3];
				mPtr++;
			}
			break;
		case 5:
			for ( i = 0; i < numColumns; i++ ) {
				dstPtr[i] = *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1] + *(mPtr+2*numColumns) * vPtr[2] +
						*(mPtr+3*numColumns) * vPtr[3] + *(mPtr+4*numColumns) * vPtr[4];
				mPtr++;
			}
			break;
		case 6:
			for ( i = 0; i < numColumns; i++ ) {
				dstPtr[i] = *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1] + *(mPtr+2*numColumns) * vPtr[2] +
						*(mPtr+3*numColumns) * vPtr[3] + *(mPtr+4*numColumns) * vPtr[4] + *(mPtr+5*numColumns) * vPtr[5];
				mPtr++;
			}
			break;
		default:
			int numRows = mat.getNumRows();
			for ( i = 0; i < numColumns; i++ ) {
				mPtr = mat.toFloatPtr() + i;
				float sum = mPtr[0] * vPtr[0];
				for ( j = 1; j < numRows; j++ ) {
					mPtr += numColumns;
					sum += mPtr[0] * vPtr[j];
				}
				dstPtr[i] = sum;
			}
			break;
	}
}

/*
============
CSIMD_Generic::matX_TransposeMultiplyAddVecX
============
*/
void VPCALL CSIMD_Generic::matX_TransposeMultiplyAddVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec ) {
	int i, j, numColumns;
	const float *mPtr, *vPtr;
	float *dstPtr;

	SMF_ASSERT( vec.getSize() >= mat.getNumRows() );
	SMF_ASSERT( dst.getSize() >= mat.getNumColumns() );

	mPtr = mat.toFloatPtr();
	vPtr = vec.toFloatPtr();
	dstPtr = dst.toFloatPtr();
	numColumns = mat.getNumColumns();
	switch( mat.getNumRows() ) {
		case 1:
			for ( i = 0; i < numColumns; i++ ) {
				dstPtr[i] += *(mPtr) * vPtr[0];
				mPtr++;
			}
			break;
		case 2:
			for ( i = 0; i < numColumns; i++ ) {
				dstPtr[i] += *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1];
				mPtr++;
			}
			break;
		case 3:
			for ( i = 0; i < numColumns; i++ ) {
				dstPtr[i] += *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1] + *(mPtr+2*numColumns) * vPtr[2];
				mPtr++;
			}
			break;
		case 4:
			for ( i = 0; i < numColumns; i++ ) {
				dstPtr[i] += *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1] + *(mPtr+2*numColumns) * vPtr[2] +
						*(mPtr+3*numColumns) * vPtr[3];
				mPtr++;
			}
			break;
		case 5:
			for ( i = 0; i < numColumns; i++ ) {
				dstPtr[i] += *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1] + *(mPtr+2*numColumns) * vPtr[2] +
						*(mPtr+3*numColumns) * vPtr[3] + *(mPtr+4*numColumns) * vPtr[4];
				mPtr++;
			}
			break;
		case 6:
			for ( i = 0; i < numColumns; i++ ) {
				dstPtr[i] += *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1] + *(mPtr+2*numColumns) * vPtr[2] +
						*(mPtr+3*numColumns) * vPtr[3] + *(mPtr+4*numColumns) * vPtr[4] + *(mPtr+5*numColumns) * vPtr[5];
				mPtr++;
			}
			break;
		default:
			int numRows = mat.getNumRows();
			for ( i = 0; i < numColumns; i++ ) {
				mPtr = mat.toFloatPtr() + i;
				float sum = mPtr[0] * vPtr[0];
				for ( j = 1; j < numRows; j++ ) {
					mPtr += numColumns;
					sum += mPtr[0] * vPtr[j];
				}
				dstPtr[i] += sum;
			}
			break;
	}
}

/*
============
CSIMD_Generic::matX_TransposeMultiplySubVecX
============
*/
void VPCALL CSIMD_Generic::matX_TransposeMultiplySubVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec ) {
	int i, numColumns;
	const float *mPtr, *vPtr;
	float *dstPtr;

	SMF_ASSERT( vec.getSize() >= mat.getNumRows() );
	SMF_ASSERT( dst.getSize() >= mat.getNumColumns() );

	mPtr = mat.toFloatPtr();
	vPtr = vec.toFloatPtr();
	dstPtr = dst.toFloatPtr();
	numColumns = mat.getNumColumns();
	switch( mat.getNumRows() ) {
		case 1:
			for ( i = 0; i < numColumns; i++ ) {
				dstPtr[i] -= *(mPtr) * vPtr[0];
				mPtr++;
			}
			break;
		case 2:
			for ( i = 0; i < numColumns; i++ ) {
				dstPtr[i] -= *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1];
				mPtr++;
			}
			break;
		case 3:
			for ( i = 0; i < numColumns; i++ ) {
				dstPtr[i] -= *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1] + *(mPtr+2*numColumns) * vPtr[2];
				mPtr++;
			}
			break;
		case 4:
			for ( i = 0; i < numColumns; i++ ) {
				dstPtr[i] -= *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1] + *(mPtr+2*numColumns) * vPtr[2] +
						*(mPtr+3*numColumns) * vPtr[3];
				mPtr++;
			}
			break;
		case 5:
			for ( i = 0; i < numColumns; i++ ) {
				dstPtr[i] -= *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1] + *(mPtr+2*numColumns) * vPtr[2] +
						*(mPtr+3*numColumns) * vPtr[3] + *(mPtr+4*numColumns) * vPtr[4];
				mPtr++;
			}
			break;
		case 6:
			for ( i = 0; i < numColumns; i++ ) {
				dstPtr[i] -= *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1] + *(mPtr+2*numColumns) * vPtr[2] +
						*(mPtr+3*numColumns) * vPtr[3] + *(mPtr+4*numColumns) * vPtr[4] + *(mPtr+5*numColumns) * vPtr[5];
				mPtr++;
			}
			break;
		default:
			int numRows = mat.getNumRows();
			for ( i = 0; i < numColumns; i++ ) {
				mPtr = mat.toFloatPtr() + i;
				float sum = mPtr[0] * vPtr[0];
				for ( int j = 1; j < numRows; j++ ) {
					mPtr += numColumns;
					sum += mPtr[0] * vPtr[j];
				}
				dstPtr[i] -= sum;
			}
			break;
	}
}

/*
============
CSIMD_Generic::matX_MultiplyMatX

	optimizes the following matrix multiplications:

	NxN * Nx6
	6xN * Nx6
	Nx6 * 6xN
	6x6 * 6xN

	with N in the range [1-6].
============
*/
void VPCALL CSIMD_Generic::matX_MultiplyMatX( CMatXD &dst, const CMatXD &m1, const CMatXD &m2 ) {
	int i, j, k, l, n;
	float *dstPtr;
	const float *m1Ptr, *m2Ptr;
	double sum;

	SMF_ASSERT( m1.getNumColumns() == m2.getNumRows() );

	dstPtr = dst.toFloatPtr();
	m1Ptr = m1.toFloatPtr();
	m2Ptr = m2.toFloatPtr();
	k = m1.getNumRows();
	l = m2.getNumColumns();

	switch( m1.getNumColumns() ) {
		case 1: {
			if ( l == 6 ) {
				for ( i = 0; i < k; i++ ) {		// Nx1 * 1x6
					*dstPtr++ = m1Ptr[i] * m2Ptr[0];
					*dstPtr++ = m1Ptr[i] * m2Ptr[1];
					*dstPtr++ = m1Ptr[i] * m2Ptr[2];
					*dstPtr++ = m1Ptr[i] * m2Ptr[3];
					*dstPtr++ = m1Ptr[i] * m2Ptr[4];
					*dstPtr++ = m1Ptr[i] * m2Ptr[5];
				}
				return;
			}
			for ( i = 0; i < k; i++ ) {
				m2Ptr = m2.toFloatPtr();
				for ( j = 0; j < l; j++ ) {
					*dstPtr++ = m1Ptr[0] * m2Ptr[0];
					m2Ptr++;
				}
				m1Ptr++;
			}
			break;
		}
		case 2: {
			if ( l == 6 ) {
				for ( i = 0; i < k; i++ ) {		// Nx2 * 2x6
					*dstPtr++ = m1Ptr[0] * m2Ptr[0] + m1Ptr[1] * m2Ptr[6];
					*dstPtr++ = m1Ptr[0] * m2Ptr[1] + m1Ptr[1] * m2Ptr[7];
					*dstPtr++ = m1Ptr[0] * m2Ptr[2] + m1Ptr[1] * m2Ptr[8];
					*dstPtr++ = m1Ptr[0] * m2Ptr[3] + m1Ptr[1] * m2Ptr[9];
					*dstPtr++ = m1Ptr[0] * m2Ptr[4] + m1Ptr[1] * m2Ptr[10];
					*dstPtr++ = m1Ptr[0] * m2Ptr[5] + m1Ptr[1] * m2Ptr[11];
					m1Ptr += 2;
				}
				return;
			}
			for ( i = 0; i < k; i++ ) {
				m2Ptr = m2.toFloatPtr();
				for ( j = 0; j < l; j++ ) {
					*dstPtr++ = m1Ptr[0] * m2Ptr[0] + m1Ptr[1] * m2Ptr[l];
					m2Ptr++;
				}
				m1Ptr += 2;
			}
			break;
		}
		case 3: {
			if ( l == 6 ) {
				for ( i = 0; i < k; i++ ) {		// Nx3 * 3x6
					*dstPtr++ = m1Ptr[0] * m2Ptr[0] + m1Ptr[1] * m2Ptr[6] + m1Ptr[2] * m2Ptr[12];
					*dstPtr++ = m1Ptr[0] * m2Ptr[1] + m1Ptr[1] * m2Ptr[7] + m1Ptr[2] * m2Ptr[13];
					*dstPtr++ = m1Ptr[0] * m2Ptr[2] + m1Ptr[1] * m2Ptr[8] + m1Ptr[2] * m2Ptr[14];
					*dstPtr++ = m1Ptr[0] * m2Ptr[3] + m1Ptr[1] * m2Ptr[9] + m1Ptr[2] * m2Ptr[15];
					*dstPtr++ = m1Ptr[0] * m2Ptr[4] + m1Ptr[1] * m2Ptr[10] + m1Ptr[2] * m2Ptr[16];
					*dstPtr++ = m1Ptr[0] * m2Ptr[5] + m1Ptr[1] * m2Ptr[11] + m1Ptr[2] * m2Ptr[17];
					m1Ptr += 3;
				}
				return;
			}
			for ( i = 0; i < k; i++ ) {
				m2Ptr = m2.toFloatPtr();
				for ( j = 0; j < l; j++ ) {
					*dstPtr++ = m1Ptr[0] * m2Ptr[0] + m1Ptr[1] * m2Ptr[l] + m1Ptr[2] * m2Ptr[2*l];
					m2Ptr++;
				}
				m1Ptr += 3;
			}
			break;
		}
		case 4: {
			if ( l == 6 ) {
				for ( i = 0; i < k; i++ ) {		// Nx4 * 4x6
					*dstPtr++ = m1Ptr[0] * m2Ptr[0] + m1Ptr[1] * m2Ptr[6] + m1Ptr[2] * m2Ptr[12] + m1Ptr[3] * m2Ptr[18];
					*dstPtr++ = m1Ptr[0] * m2Ptr[1] + m1Ptr[1] * m2Ptr[7] + m1Ptr[2] * m2Ptr[13] + m1Ptr[3] * m2Ptr[19];
					*dstPtr++ = m1Ptr[0] * m2Ptr[2] + m1Ptr[1] * m2Ptr[8] + m1Ptr[2] * m2Ptr[14] + m1Ptr[3] * m2Ptr[20];
					*dstPtr++ = m1Ptr[0] * m2Ptr[3] + m1Ptr[1] * m2Ptr[9] + m1Ptr[2] * m2Ptr[15] + m1Ptr[3] * m2Ptr[21];
					*dstPtr++ = m1Ptr[0] * m2Ptr[4] + m1Ptr[1] * m2Ptr[10] + m1Ptr[2] * m2Ptr[16] + m1Ptr[3] * m2Ptr[22];
					*dstPtr++ = m1Ptr[0] * m2Ptr[5] + m1Ptr[1] * m2Ptr[11] + m1Ptr[2] * m2Ptr[17] + m1Ptr[3] * m2Ptr[23];
					m1Ptr += 4;
				}
				return;
			}
			for ( i = 0; i < k; i++ ) {
				m2Ptr = m2.toFloatPtr();
				for ( j = 0; j < l; j++ ) {
					*dstPtr++ = m1Ptr[0] * m2Ptr[0] + m1Ptr[1] * m2Ptr[l] + m1Ptr[2] * m2Ptr[2*l] +
									 m1Ptr[3] * m2Ptr[3*l];
					m2Ptr++;
				}
				m1Ptr += 4;
			}
			break;
		}
		case 5: {
			if ( l == 6 ) {
				for ( i = 0; i < k; i++ ) {		// Nx5 * 5x6
					*dstPtr++ = m1Ptr[0] * m2Ptr[0] + m1Ptr[1] * m2Ptr[6] + m1Ptr[2] * m2Ptr[12] + m1Ptr[3] * m2Ptr[18] + m1Ptr[4] * m2Ptr[24];
					*dstPtr++ = m1Ptr[0] * m2Ptr[1] + m1Ptr[1] * m2Ptr[7] + m1Ptr[2] * m2Ptr[13] + m1Ptr[3] * m2Ptr[19] + m1Ptr[4] * m2Ptr[25];
					*dstPtr++ = m1Ptr[0] * m2Ptr[2] + m1Ptr[1] * m2Ptr[8] + m1Ptr[2] * m2Ptr[14] + m1Ptr[3] * m2Ptr[20] + m1Ptr[4] * m2Ptr[26];
					*dstPtr++ = m1Ptr[0] * m2Ptr[3] + m1Ptr[1] * m2Ptr[9] + m1Ptr[2] * m2Ptr[15] + m1Ptr[3] * m2Ptr[21] + m1Ptr[4] * m2Ptr[27];
					*dstPtr++ = m1Ptr[0] * m2Ptr[4] + m1Ptr[1] * m2Ptr[10] + m1Ptr[2] * m2Ptr[16] + m1Ptr[3] * m2Ptr[22] + m1Ptr[4] * m2Ptr[28];
					*dstPtr++ = m1Ptr[0] * m2Ptr[5] + m1Ptr[1] * m2Ptr[11] + m1Ptr[2] * m2Ptr[17] + m1Ptr[3] * m2Ptr[23] + m1Ptr[4] * m2Ptr[29];
					m1Ptr += 5;
				}
				return;
			}
			for ( i = 0; i < k; i++ ) {
				m2Ptr = m2.toFloatPtr();
				for ( j = 0; j < l; j++ ) {
					*dstPtr++ = m1Ptr[0] * m2Ptr[0] + m1Ptr[1] * m2Ptr[l] + m1Ptr[2] * m2Ptr[2*l] +
									 m1Ptr[3] * m2Ptr[3*l] + m1Ptr[4] * m2Ptr[4*l];
					m2Ptr++;
				}
				m1Ptr += 5;
			}
			break;
		}
		case 6: {
			switch( k ) {
				case 1: {
					if ( l == 1 ) {		// 1x6 * 6x1
						dstPtr[0] = m1Ptr[0] * m2Ptr[0] + m1Ptr[1] * m2Ptr[1] + m1Ptr[2] * m2Ptr[2] +
									 m1Ptr[3] * m2Ptr[3] + m1Ptr[4] * m2Ptr[4] + m1Ptr[5] * m2Ptr[5];
						return;
					}
					break;
				}
				case 2: {
					if ( l == 2 ) {		// 2x6 * 6x2
						for ( i = 0; i < 2; i++ ) {
							for ( j = 0; j < 2; j++ ) {
								*dstPtr = m1Ptr[0] * m2Ptr[ 0 * 2 + j ]
										+ m1Ptr[1] * m2Ptr[ 1 * 2 + j ]
										+ m1Ptr[2] * m2Ptr[ 2 * 2 + j ]
										+ m1Ptr[3] * m2Ptr[ 3 * 2 + j ]
										+ m1Ptr[4] * m2Ptr[ 4 * 2 + j ]
										+ m1Ptr[5] * m2Ptr[ 5 * 2 + j ];
								dstPtr++;
							}
							m1Ptr += 6;
						}
						return;
					}
					break;
				}
				case 3: {
					if ( l == 3 ) {		// 3x6 * 6x3
						for ( i = 0; i < 3; i++ ) {
							for ( j = 0; j < 3; j++ ) {
								*dstPtr = m1Ptr[0] * m2Ptr[ 0 * 3 + j ]
										+ m1Ptr[1] * m2Ptr[ 1 * 3 + j ]
										+ m1Ptr[2] * m2Ptr[ 2 * 3 + j ]
										+ m1Ptr[3] * m2Ptr[ 3 * 3 + j ]
										+ m1Ptr[4] * m2Ptr[ 4 * 3 + j ]
										+ m1Ptr[5] * m2Ptr[ 5 * 3 + j ];
								dstPtr++;
							}
							m1Ptr += 6;
						}
						return;
					}
					break;
				}
				case 4: {
					if ( l == 4 ) {		// 4x6 * 6x4
						for ( i = 0; i < 4; i++ ) {
							for ( j = 0; j < 4; j++ ) {
								*dstPtr = m1Ptr[0] * m2Ptr[ 0 * 4 + j ]
										+ m1Ptr[1] * m2Ptr[ 1 * 4 + j ]
										+ m1Ptr[2] * m2Ptr[ 2 * 4 + j ]
										+ m1Ptr[3] * m2Ptr[ 3 * 4 + j ]
										+ m1Ptr[4] * m2Ptr[ 4 * 4 + j ]
										+ m1Ptr[5] * m2Ptr[ 5 * 4 + j ];
								dstPtr++;
							}
							m1Ptr += 6;
						}
						return;
					}
				}
				case 5: {
					if ( l == 5 ) {		// 5x6 * 6x5
						for ( i = 0; i < 5; i++ ) {
							for ( j = 0; j < 5; j++ ) {
								*dstPtr = m1Ptr[0] * m2Ptr[ 0 * 5 + j ]
										+ m1Ptr[1] * m2Ptr[ 1 * 5 + j ]
										+ m1Ptr[2] * m2Ptr[ 2 * 5 + j ]
										+ m1Ptr[3] * m2Ptr[ 3 * 5 + j ]
										+ m1Ptr[4] * m2Ptr[ 4 * 5 + j ]
										+ m1Ptr[5] * m2Ptr[ 5 * 5 + j ];
								dstPtr++;
							}
							m1Ptr += 6;
						}
						return;
					}
				}
				case 6: {
					switch( l ) {
						case 1: {		// 6x6 * 6x1
							for ( i = 0; i < 6; i++ ) {
								*dstPtr = m1Ptr[0] * m2Ptr[ 0 * 1 ]
										+ m1Ptr[1] * m2Ptr[ 1 * 1 ]
										+ m1Ptr[2] * m2Ptr[ 2 * 1 ]
										+ m1Ptr[3] * m2Ptr[ 3 * 1 ]
										+ m1Ptr[4] * m2Ptr[ 4 * 1 ]
										+ m1Ptr[5] * m2Ptr[ 5 * 1 ];
								dstPtr++;
								m1Ptr += 6;
							}
							return;
						}
						case 2: {		// 6x6 * 6x2
							for ( i = 0; i < 6; i++ ) {
								for ( j = 0; j < 2; j++ ) {
									*dstPtr = m1Ptr[0] * m2Ptr[ 0 * 2 + j ]
											+ m1Ptr[1] * m2Ptr[ 1 * 2 + j ]
											+ m1Ptr[2] * m2Ptr[ 2 * 2 + j ]
											+ m1Ptr[3] * m2Ptr[ 3 * 2 + j ]
											+ m1Ptr[4] * m2Ptr[ 4 * 2 + j ]
											+ m1Ptr[5] * m2Ptr[ 5 * 2 + j ];
									dstPtr++;
								}
								m1Ptr += 6;
							}
							return;
						}
						case 3: {		// 6x6 * 6x3
							for ( i = 0; i < 6; i++ ) {
								for ( j = 0; j < 3; j++ ) {
									*dstPtr = m1Ptr[0] * m2Ptr[ 0 * 3 + j ]
											+ m1Ptr[1] * m2Ptr[ 1 * 3 + j ]
											+ m1Ptr[2] * m2Ptr[ 2 * 3 + j ]
											+ m1Ptr[3] * m2Ptr[ 3 * 3 + j ]
											+ m1Ptr[4] * m2Ptr[ 4 * 3 + j ]
											+ m1Ptr[5] * m2Ptr[ 5 * 3 + j ];
									dstPtr++;
								}
								m1Ptr += 6;
							}
							return;
						}
						case 4: {		// 6x6 * 6x4
							for ( i = 0; i < 6; i++ ) {
								for ( j = 0; j < 4; j++ ) {
									*dstPtr = m1Ptr[0] * m2Ptr[ 0 * 4 + j ]
											+ m1Ptr[1] * m2Ptr[ 1 * 4 + j ]
											+ m1Ptr[2] * m2Ptr[ 2 * 4 + j ]
											+ m1Ptr[3] * m2Ptr[ 3 * 4 + j ]
											+ m1Ptr[4] * m2Ptr[ 4 * 4 + j ]
											+ m1Ptr[5] * m2Ptr[ 5 * 4 + j ];
									dstPtr++;
								}
								m1Ptr += 6;
							}
							return;
						}
						case 5: {		// 6x6 * 6x5
							for ( i = 0; i < 6; i++ ) {
								for ( j = 0; j < 5; j++ ) {
									*dstPtr = m1Ptr[0] * m2Ptr[ 0 * 5 + j ]
											+ m1Ptr[1] * m2Ptr[ 1 * 5 + j ]
											+ m1Ptr[2] * m2Ptr[ 2 * 5 + j ]
											+ m1Ptr[3] * m2Ptr[ 3 * 5 + j ]
											+ m1Ptr[4] * m2Ptr[ 4 * 5 + j ]
											+ m1Ptr[5] * m2Ptr[ 5 * 5 + j ];
									dstPtr++;
								}
								m1Ptr += 6;
							}
							return;
						}
						case 6: {		// 6x6 * 6x6
							for ( i = 0; i < 6; i++ ) {
								for ( j = 0; j < 6; j++ ) {
									*dstPtr = m1Ptr[0] * m2Ptr[ 0 * 6 + j ]
											+ m1Ptr[1] * m2Ptr[ 1 * 6 + j ]
											+ m1Ptr[2] * m2Ptr[ 2 * 6 + j ]
											+ m1Ptr[3] * m2Ptr[ 3 * 6 + j ]
											+ m1Ptr[4] * m2Ptr[ 4 * 6 + j ]
											+ m1Ptr[5] * m2Ptr[ 5 * 6 + j ];
									dstPtr++;
								}
								m1Ptr += 6;
							}
							return;
						}
					}
				}
			}
			for ( i = 0; i < k; i++ ) {
				m2Ptr = m2.toFloatPtr();
				for ( j = 0; j < l; j++ ) {
					*dstPtr++ = m1Ptr[0] * m2Ptr[0] + m1Ptr[1] * m2Ptr[l] + m1Ptr[2] * m2Ptr[2*l] +
									 m1Ptr[3] * m2Ptr[3*l] + m1Ptr[4] * m2Ptr[4*l] + m1Ptr[5] * m2Ptr[5*l];
					m2Ptr++;
				}
				m1Ptr += 6;
			}
			break;
		}
		default: {
			for ( i = 0; i < k; i++ ) {
				for ( j = 0; j < l; j++ ) {
					m2Ptr = m2.toFloatPtr() + j;
					sum = m1Ptr[0] * m2Ptr[0];
					for ( n = 1; n < m1.getNumColumns(); n++ ) {
						m2Ptr += l;
						sum += m1Ptr[n] * m2Ptr[0];
					}
					*dstPtr++ = sum;
				}
				m1Ptr += m1.getNumColumns();
			}
			break;
		}
	}
}

/*
============
CSIMD_Generic::matX_TransposeMultiplyMatX

	optimizes the following tranpose matrix multiplications:

	Nx6 * NxN
	6xN * 6x6

	with N in the range [1-6].
============
*/
void VPCALL CSIMD_Generic::matX_TransposeMultiplyMatX( CMatXD &dst, const CMatXD &m1, const CMatXD &m2 ) {
	int i, j, k, l, n;
	float *dstPtr;
	const float *m1Ptr, *m2Ptr;
	double sum;

	SMF_ASSERT( m1.getNumRows() == m2.getNumRows() );

	m1Ptr = m1.toFloatPtr();
	m2Ptr = m2.toFloatPtr();
	dstPtr = dst.toFloatPtr();
	k = m1.getNumColumns();
	l = m2.getNumColumns();

	switch( m1.getNumRows() ) {
		case 1:
			if ( k == 6 && l == 1 ) {			// 1x6 * 1x1
				for ( i = 0; i < 6; i++ ) {
					*dstPtr++ = m1Ptr[0] * m2Ptr[0];
					m1Ptr++;
				}
				return;
			}
			for ( i = 0; i < k; i++ ) {
				m2Ptr = m2.toFloatPtr();
				for ( j = 0; j < l; j++ ) {
					*dstPtr++ = m1Ptr[0] * m2Ptr[0];
					m2Ptr++;
				}
				m1Ptr++;
			}
			break;
		case 2:
			if ( k == 6 && l == 2 ) {			// 2x6 * 2x2
				for ( i = 0; i < 6; i++ ) {
					*dstPtr++ = m1Ptr[0*6] * m2Ptr[0*2+0] + m1Ptr[1*6] * m2Ptr[1*2+0];
					*dstPtr++ = m1Ptr[0*6] * m2Ptr[0*2+1] + m1Ptr[1*6] * m2Ptr[1*2+1];
					m1Ptr++;
				}
				return;
			}
			for ( i = 0; i < k; i++ ) {
				m2Ptr = m2.toFloatPtr();
				for ( j = 0; j < l; j++ ) {
					*dstPtr++ = m1Ptr[0] * m2Ptr[0] + m1Ptr[k] * m2Ptr[l];
					m2Ptr++;
				}
				m1Ptr++;
			}
			break;
		case 3:
			if ( k == 6 && l == 3 ) {			// 3x6 * 3x3
				for ( i = 0; i < 6; i++ ) {
					*dstPtr++ = m1Ptr[0*6] * m2Ptr[0*3+0] + m1Ptr[1*6] * m2Ptr[1*3+0] + m1Ptr[2*6] * m2Ptr[2*3+0];
					*dstPtr++ = m1Ptr[0*6] * m2Ptr[0*3+1] + m1Ptr[1*6] * m2Ptr[1*3+1] + m1Ptr[2*6] * m2Ptr[2*3+1];
					*dstPtr++ = m1Ptr[0*6] * m2Ptr[0*3+2] + m1Ptr[1*6] * m2Ptr[1*3+2] + m1Ptr[2*6] * m2Ptr[2*3+2];
					m1Ptr++;
				}
				return;
			}
			for ( i = 0; i < k; i++ ) {
				m2Ptr = m2.toFloatPtr();
				for ( j = 0; j < l; j++ ) {
					*dstPtr++ = m1Ptr[0] * m2Ptr[0] + m1Ptr[k] * m2Ptr[l] + m1Ptr[2*k] * m2Ptr[2*l];
					m2Ptr++;
				}
				m1Ptr++;
			}
			break;
		case 4:
			if ( k == 6 && l == 4 ) {			// 4x6 * 4x4
				for ( i = 0; i < 6; i++ ) {
					*dstPtr++ = m1Ptr[0*6] * m2Ptr[0*4+0] + m1Ptr[1*6] * m2Ptr[1*4+0] + m1Ptr[2*6] * m2Ptr[2*4+0] + m1Ptr[3*6] * m2Ptr[3*4+0];
					*dstPtr++ = m1Ptr[0*6] * m2Ptr[0*4+1] + m1Ptr[1*6] * m2Ptr[1*4+1] + m1Ptr[2*6] * m2Ptr[2*4+1] + m1Ptr[3*6] * m2Ptr[3*4+1];
					*dstPtr++ = m1Ptr[0*6] * m2Ptr[0*4+2] + m1Ptr[1*6] * m2Ptr[1*4+2] + m1Ptr[2*6] * m2Ptr[2*4+2] + m1Ptr[3*6] * m2Ptr[3*4+2];
					*dstPtr++ = m1Ptr[0*6] * m2Ptr[0*4+3] + m1Ptr[1*6] * m2Ptr[1*4+3] + m1Ptr[2*6] * m2Ptr[2*4+3] + m1Ptr[3*6] * m2Ptr[3*4+3];
					m1Ptr++;
				}
				return;
			}
			for ( i = 0; i < k; i++ ) {
				m2Ptr = m2.toFloatPtr();
				for ( j = 0; j < l; j++ ) {
					*dstPtr++ = m1Ptr[0] * m2Ptr[0] + m1Ptr[k] * m2Ptr[l] + m1Ptr[2*k] * m2Ptr[2*l] +
									m1Ptr[3*k] * m2Ptr[3*l];
					m2Ptr++;
				}
				m1Ptr++;
			}
			break;
		case 5:
			if ( k == 6 && l == 5 ) {			// 5x6 * 5x5
				for ( i = 0; i < 6; i++ ) {
					*dstPtr++ = m1Ptr[0*6] * m2Ptr[0*5+0] + m1Ptr[1*6] * m2Ptr[1*5+0] + m1Ptr[2*6] * m2Ptr[2*5+0] + m1Ptr[3*6] * m2Ptr[3*5+0] + m1Ptr[4*6] * m2Ptr[4*5+0];
					*dstPtr++ = m1Ptr[0*6] * m2Ptr[0*5+1] + m1Ptr[1*6] * m2Ptr[1*5+1] + m1Ptr[2*6] * m2Ptr[2*5+1] + m1Ptr[3*6] * m2Ptr[3*5+1] + m1Ptr[4*6] * m2Ptr[4*5+1];
					*dstPtr++ = m1Ptr[0*6] * m2Ptr[0*5+2] + m1Ptr[1*6] * m2Ptr[1*5+2] + m1Ptr[2*6] * m2Ptr[2*5+2] + m1Ptr[3*6] * m2Ptr[3*5+2] + m1Ptr[4*6] * m2Ptr[4*5+2];
					*dstPtr++ = m1Ptr[0*6] * m2Ptr[0*5+3] + m1Ptr[1*6] * m2Ptr[1*5+3] + m1Ptr[2*6] * m2Ptr[2*5+3] + m1Ptr[3*6] * m2Ptr[3*5+3] + m1Ptr[4*6] * m2Ptr[4*5+3];
					*dstPtr++ = m1Ptr[0*6] * m2Ptr[0*5+4] + m1Ptr[1*6] * m2Ptr[1*5+4] + m1Ptr[2*6] * m2Ptr[2*5+4] + m1Ptr[3*6] * m2Ptr[3*5+4] + m1Ptr[4*6] * m2Ptr[4*5+4];
					m1Ptr++;
				}
				return;
			}
			for ( i = 0; i < k; i++ ) {
				m2Ptr = m2.toFloatPtr();
				for ( j = 0; j < l; j++ ) {
					*dstPtr++ = m1Ptr[0] * m2Ptr[0] + m1Ptr[k] * m2Ptr[l] + m1Ptr[2*k] * m2Ptr[2*l] +
									m1Ptr[3*k] * m2Ptr[3*l] + m1Ptr[4*k] * m2Ptr[4*l];
					m2Ptr++;
				}
				m1Ptr++;
			}
			break;
		case 6:
			if ( l == 6 ) {
				switch( k ) {
					case 1:						// 6x1 * 6x6
						m2Ptr = m2.toFloatPtr();
						for ( j = 0; j < 6; j++ ) {
							*dstPtr++ = m1Ptr[0*1] * m2Ptr[0*6] +
										m1Ptr[1*1] * m2Ptr[1*6] +
										m1Ptr[2*1] * m2Ptr[2*6] +
										m1Ptr[3*1] * m2Ptr[3*6] +
										m1Ptr[4*1] * m2Ptr[4*6] +
										m1Ptr[5*1] * m2Ptr[5*6];
							m2Ptr++;
						}
						return;
					case 2:						// 6x2 * 6x6
						for ( i = 0; i < 2; i++ ) {
							m2Ptr = m2.toFloatPtr();
							for ( j = 0; j < 6; j++ ) {
								*dstPtr++ = m1Ptr[0*2] * m2Ptr[0*6] +
											m1Ptr[1*2] * m2Ptr[1*6] +
											m1Ptr[2*2] * m2Ptr[2*6] +
											m1Ptr[3*2] * m2Ptr[3*6] +
											m1Ptr[4*2] * m2Ptr[4*6] +
											m1Ptr[5*2] * m2Ptr[5*6];
								m2Ptr++;
							}
							m1Ptr++;
						}
						return;
					case 3:						// 6x3 * 6x6
						for ( i = 0; i < 3; i++ ) {
							m2Ptr = m2.toFloatPtr();
							for ( j = 0; j < 6; j++ ) {
								*dstPtr++ = m1Ptr[0*3] * m2Ptr[0*6] +
											m1Ptr[1*3] * m2Ptr[1*6] +
											m1Ptr[2*3] * m2Ptr[2*6] +
											m1Ptr[3*3] * m2Ptr[3*6] +
											m1Ptr[4*3] * m2Ptr[4*6] +
											m1Ptr[5*3] * m2Ptr[5*6];
								m2Ptr++;
							}
							m1Ptr++;
						}
						return;
					case 4:						// 6x4 * 6x6
						for ( i = 0; i < 4; i++ ) {
							m2Ptr = m2.toFloatPtr();
							for ( j = 0; j < 6; j++ ) {
								*dstPtr++ = m1Ptr[0*4] * m2Ptr[0*6] +
											m1Ptr[1*4] * m2Ptr[1*6] +
											m1Ptr[2*4] * m2Ptr[2*6] +
											m1Ptr[3*4] * m2Ptr[3*6] +
											m1Ptr[4*4] * m2Ptr[4*6] +
											m1Ptr[5*4] * m2Ptr[5*6];
								m2Ptr++;
							}
							m1Ptr++;
						}
						return;
					case 5:						// 6x5 * 6x6
						for ( i = 0; i < 5; i++ ) {
							m2Ptr = m2.toFloatPtr();
							for ( j = 0; j < 6; j++ ) {
								*dstPtr++ = m1Ptr[0*5] * m2Ptr[0*6] +
											m1Ptr[1*5] * m2Ptr[1*6] +
											m1Ptr[2*5] * m2Ptr[2*6] +
											m1Ptr[3*5] * m2Ptr[3*6] +
											m1Ptr[4*5] * m2Ptr[4*6] +
											m1Ptr[5*5] * m2Ptr[5*6];
								m2Ptr++;
							}
							m1Ptr++;
						}
						return;
					case 6:						// 6x6 * 6x6
						for ( i = 0; i < 6; i++ ) {
							m2Ptr = m2.toFloatPtr();
							for ( j = 0; j < 6; j++ ) {
								*dstPtr++ = m1Ptr[0*6] * m2Ptr[0*6] +
											m1Ptr[1*6] * m2Ptr[1*6] +
											m1Ptr[2*6] * m2Ptr[2*6] +
											m1Ptr[3*6] * m2Ptr[3*6] +
											m1Ptr[4*6] * m2Ptr[4*6] +
											m1Ptr[5*6] * m2Ptr[5*6];
								m2Ptr++;
							}
							m1Ptr++;
						}
						return;
				}
			}
			for ( i = 0; i < k; i++ ) {
				m2Ptr = m2.toFloatPtr();
				for ( j = 0; j < l; j++ ) {
					*dstPtr++ = m1Ptr[0] * m2Ptr[0] + m1Ptr[k] * m2Ptr[l] + m1Ptr[2*k] * m2Ptr[2*l] +
									m1Ptr[3*k] * m2Ptr[3*l] + m1Ptr[4*k] * m2Ptr[4*l] + m1Ptr[5*k] * m2Ptr[5*l];
					m2Ptr++;
				}
				m1Ptr++;
			}
			break;
		default:
			for ( i = 0; i < k; i++ ) {
				for ( j = 0; j < l; j++ ) {
					m1Ptr = m1.toFloatPtr() + i;
					m2Ptr = m2.toFloatPtr() + j;
					sum = m1Ptr[0] * m2Ptr[0];
					for ( n = 1; n < m1.getNumRows(); n++ ) {
						m1Ptr += k;
						m2Ptr += l;
						sum += m1Ptr[0] * m2Ptr[0];
					}
					*dstPtr++ = sum;
				}
			}
		break;
	}
}

/*
============
CSIMD_Generic::matX_LowerTriangularsolve

  solves x in Lx = b for the n * n sub-matrix of L
  if skip > 0 the first skip elements of x are assumed to be valid already
  L has to be a lower triangular matrix with (implicit) ones on the diagonal
  x == b is allowed
============
*/
void VPCALL CSIMD_Generic::matX_LowerTriangularsolve( const CMatXD &L, float *x, const float *b, const int n, int skip ) {
#if 1

	int nc;
	const float *lptr;

	if ( skip >= n ) {
		return;
	}

	lptr = L.toFloatPtr();
	nc = L.getNumColumns();

	// unrolled cases for n < 8
	if ( n < 8 ) {
		#define NSKIP( n, s )	((n<<3)|(s&7))
		switch( NSKIP( n, skip ) ) {
			case NSKIP( 1, 0 ): x[0] = b[0];
				return;
			case NSKIP( 2, 0 ): x[0] = b[0];
			case NSKIP( 2, 1 ): x[1] = b[1] - lptr[1*nc+0] * x[0];
				return;
			case NSKIP( 3, 0 ): x[0] = b[0];
			case NSKIP( 3, 1 ): x[1] = b[1] - lptr[1*nc+0] * x[0];
			case NSKIP( 3, 2 ): x[2] = b[2] - lptr[2*nc+0] * x[0] - lptr[2*nc+1] * x[1];
				return;
			case NSKIP( 4, 0 ): x[0] = b[0];
			case NSKIP( 4, 1 ): x[1] = b[1] - lptr[1*nc+0] * x[0];
			case NSKIP( 4, 2 ): x[2] = b[2] - lptr[2*nc+0] * x[0] - lptr[2*nc+1] * x[1];
			case NSKIP( 4, 3 ): x[3] = b[3] - lptr[3*nc+0] * x[0] - lptr[3*nc+1] * x[1] - lptr[3*nc+2] * x[2];
				return;
			case NSKIP( 5, 0 ): x[0] = b[0];
			case NSKIP( 5, 1 ): x[1] = b[1] - lptr[1*nc+0] * x[0];
			case NSKIP( 5, 2 ): x[2] = b[2] - lptr[2*nc+0] * x[0] - lptr[2*nc+1] * x[1];
			case NSKIP( 5, 3 ): x[3] = b[3] - lptr[3*nc+0] * x[0] - lptr[3*nc+1] * x[1] - lptr[3*nc+2] * x[2];
			case NSKIP( 5, 4 ): x[4] = b[4] - lptr[4*nc+0] * x[0] - lptr[4*nc+1] * x[1] - lptr[4*nc+2] * x[2] - lptr[4*nc+3] * x[3];
				return;
			case NSKIP( 6, 0 ): x[0] = b[0];
			case NSKIP( 6, 1 ): x[1] = b[1] - lptr[1*nc+0] * x[0];
			case NSKIP( 6, 2 ): x[2] = b[2] - lptr[2*nc+0] * x[0] - lptr[2*nc+1] * x[1];
			case NSKIP( 6, 3 ): x[3] = b[3] - lptr[3*nc+0] * x[0] - lptr[3*nc+1] * x[1] - lptr[3*nc+2] * x[2];
			case NSKIP( 6, 4 ): x[4] = b[4] - lptr[4*nc+0] * x[0] - lptr[4*nc+1] * x[1] - lptr[4*nc+2] * x[2] - lptr[4*nc+3] * x[3];
			case NSKIP( 6, 5 ): x[5] = b[5] - lptr[5*nc+0] * x[0] - lptr[5*nc+1] * x[1] - lptr[5*nc+2] * x[2] - lptr[5*nc+3] * x[3] - lptr[5*nc+4] * x[4];
				return;
			case NSKIP( 7, 0 ): x[0] = b[0];
			case NSKIP( 7, 1 ): x[1] = b[1] - lptr[1*nc+0] * x[0];
			case NSKIP( 7, 2 ): x[2] = b[2] - lptr[2*nc+0] * x[0] - lptr[2*nc+1] * x[1];
			case NSKIP( 7, 3 ): x[3] = b[3] - lptr[3*nc+0] * x[0] - lptr[3*nc+1] * x[1] - lptr[3*nc+2] * x[2];
			case NSKIP( 7, 4 ): x[4] = b[4] - lptr[4*nc+0] * x[0] - lptr[4*nc+1] * x[1] - lptr[4*nc+2] * x[2] - lptr[4*nc+3] * x[3];
			case NSKIP( 7, 5 ): x[5] = b[5] - lptr[5*nc+0] * x[0] - lptr[5*nc+1] * x[1] - lptr[5*nc+2] * x[2] - lptr[5*nc+3] * x[3] - lptr[5*nc+4] * x[4];
			case NSKIP( 7, 6 ): x[6] = b[6] - lptr[6*nc+0] * x[0] - lptr[6*nc+1] * x[1] - lptr[6*nc+2] * x[2] - lptr[6*nc+3] * x[3] - lptr[6*nc+4] * x[4] - lptr[6*nc+5] * x[5];
				return;
		}
		return;
	}

	// process first 4 rows
	switch( skip ) {
		case 0: x[0] = b[0];
		case 1: x[1] = b[1] - lptr[1*nc+0] * x[0];
		case 2: x[2] = b[2] - lptr[2*nc+0] * x[0] - lptr[2*nc+1] * x[1];
		case 3: x[3] = b[3] - lptr[3*nc+0] * x[0] - lptr[3*nc+1] * x[1] - lptr[3*nc+2] * x[2];
				skip = 4;
	}

	lptr = L[skip];

	int i, j;
	register double s0, s1, s2, s3;

	for ( i = skip; i < n; i++ ) {
		s0 = lptr[0] * x[0];
		s1 = lptr[1] * x[1];
		s2 = lptr[2] * x[2];
		s3 = lptr[3] * x[3];
		for ( j = 4; j < i-7; j += 8 ) {
			s0 += lptr[j+0] * x[j+0];
			s1 += lptr[j+1] * x[j+1];
			s2 += lptr[j+2] * x[j+2];
			s3 += lptr[j+3] * x[j+3];
			s0 += lptr[j+4] * x[j+4];
			s1 += lptr[j+5] * x[j+5];
			s2 += lptr[j+6] * x[j+6];
			s3 += lptr[j+7] * x[j+7];
		}
		switch( i - j ) {
			NODEFAULT;
			case 7: s0 += lptr[j+6] * x[j+6];
			case 6: s1 += lptr[j+5] * x[j+5];
			case 5: s2 += lptr[j+4] * x[j+4];
			case 4: s3 += lptr[j+3] * x[j+3];
			case 3: s0 += lptr[j+2] * x[j+2];
			case 2: s1 += lptr[j+1] * x[j+1];
			case 1: s2 += lptr[j+0] * x[j+0];
			case 0: break;
		}
		double sum;
		sum = s3;
		sum += s2;
		sum += s1;
		sum += s0;
		sum -= b[i];
		x[i] = -sum;
		lptr += nc;
	}

#else

	int i, j;
	const float *lptr;
	double sum;

	for ( i = skip; i < n; i++ ) {
		sum = b[i];
		lptr = L[i];
		for ( j = 0; j < i; j++ ) {
			sum -= lptr[j] * x[j];
		}
		x[i] = sum;
	}

#endif
}

/*
============
CSIMD_Generic::matX_LowerTriangularsolveTranspose

  solves x in L'x = b for the n * n sub-matrix of L
  L has to be a lower triangular matrix with (implicit) ones on the diagonal
  x == b is allowed
============
*/
void VPCALL CSIMD_Generic::matX_LowerTriangularsolveTranspose( const CMatXD &L, float *x, const float *b, const int n ) {
#if 1

	int nc;
	const float *lptr;

	lptr = L.toFloatPtr();
	nc = L.getNumColumns();

	// unrolled cases for n < 8
	if ( n < 8 ) {
		switch( n ) {
			case 0:
				return;
			case 1:
				x[0] = b[0];
				return;
			case 2:
				x[1] = b[1];
				x[0] = b[0] - lptr[1*nc+0] * x[1];
				return;
			case 3:
				x[2] = b[2];
				x[1] = b[1] - lptr[2*nc+1] * x[2];
				x[0] = b[0] - lptr[2*nc+0] * x[2] - lptr[1*nc+0] * x[1];
				return;
			case 4:
				x[3] = b[3];
				x[2] = b[2] - lptr[3*nc+2] * x[3];
				x[1] = b[1] - lptr[3*nc+1] * x[3] - lptr[2*nc+1] * x[2];
				x[0] = b[0] - lptr[3*nc+0] * x[3] - lptr[2*nc+0] * x[2] - lptr[1*nc+0] * x[1];
				return;
			case 5:
				x[4] = b[4];
				x[3] = b[3] - lptr[4*nc+3] * x[4];
				x[2] = b[2] - lptr[4*nc+2] * x[4] - lptr[3*nc+2] * x[3];
				x[1] = b[1] - lptr[4*nc+1] * x[4] - lptr[3*nc+1] * x[3] - lptr[2*nc+1] * x[2];
				x[0] = b[0] - lptr[4*nc+0] * x[4] - lptr[3*nc+0] * x[3] - lptr[2*nc+0] * x[2] - lptr[1*nc+0] * x[1];
				return;
			case 6:
				x[5] = b[5];
				x[4] = b[4] - lptr[5*nc+4] * x[5];
				x[3] = b[3] - lptr[5*nc+3] * x[5] - lptr[4*nc+3] * x[4];
				x[2] = b[2] - lptr[5*nc+2] * x[5] - lptr[4*nc+2] * x[4] - lptr[3*nc+2] * x[3];
				x[1] = b[1] - lptr[5*nc+1] * x[5] - lptr[4*nc+1] * x[4] - lptr[3*nc+1] * x[3] - lptr[2*nc+1] * x[2];
				x[0] = b[0] - lptr[5*nc+0] * x[5] - lptr[4*nc+0] * x[4] - lptr[3*nc+0] * x[3] - lptr[2*nc+0] * x[2] - lptr[1*nc+0] * x[1];
				return;
			case 7:
				x[6] = b[6];
				x[5] = b[5] - lptr[6*nc+5] * x[6];
				x[4] = b[4] - lptr[6*nc+4] * x[6] - lptr[5*nc+4] * x[5];
				x[3] = b[3] - lptr[6*nc+3] * x[6] - lptr[5*nc+3] * x[5] - lptr[4*nc+3] * x[4];
				x[2] = b[2] - lptr[6*nc+2] * x[6] - lptr[5*nc+2] * x[5] - lptr[4*nc+2] * x[4] - lptr[3*nc+2] * x[3];
				x[1] = b[1] - lptr[6*nc+1] * x[6] - lptr[5*nc+1] * x[5] - lptr[4*nc+1] * x[4] - lptr[3*nc+1] * x[3] - lptr[2*nc+1] * x[2];
				x[0] = b[0] - lptr[6*nc+0] * x[6] - lptr[5*nc+0] * x[5] - lptr[4*nc+0] * x[4] - lptr[3*nc+0] * x[3] - lptr[2*nc+0] * x[2] - lptr[1*nc+0] * x[1];
				return;
		}
		return;
	}

	int i, j;
	register double s0, s1, s2, s3;
	float *xptr;

	lptr = L.toFloatPtr() + n * nc + n - 4;
	xptr = x + n;

	// process 4 rows at a time
	for ( i = n; i >= 4; i -= 4 ) {
		s0 = b[i-4];
		s1 = b[i-3];
		s2 = b[i-2];
		s3 = b[i-1];
		// process 4x4 blocks
		for ( j = 0; j < n-i; j += 4 ) {
			s0 -= lptr[(j+0)*nc+0] * xptr[j+0];
			s1 -= lptr[(j+0)*nc+1] * xptr[j+0];
			s2 -= lptr[(j+0)*nc+2] * xptr[j+0];
			s3 -= lptr[(j+0)*nc+3] * xptr[j+0];
			s0 -= lptr[(j+1)*nc+0] * xptr[j+1];
			s1 -= lptr[(j+1)*nc+1] * xptr[j+1];
			s2 -= lptr[(j+1)*nc+2] * xptr[j+1];
			s3 -= lptr[(j+1)*nc+3] * xptr[j+1];
			s0 -= lptr[(j+2)*nc+0] * xptr[j+2];
			s1 -= lptr[(j+2)*nc+1] * xptr[j+2];
			s2 -= lptr[(j+2)*nc+2] * xptr[j+2];
			s3 -= lptr[(j+2)*nc+3] * xptr[j+2];
			s0 -= lptr[(j+3)*nc+0] * xptr[j+3];
			s1 -= lptr[(j+3)*nc+1] * xptr[j+3];
			s2 -= lptr[(j+3)*nc+2] * xptr[j+3];
			s3 -= lptr[(j+3)*nc+3] * xptr[j+3];
		}
		// process left over of the 4 rows
		s0 -= lptr[0-1*nc] * s3;
		s1 -= lptr[1-1*nc] * s3;
		s2 -= lptr[2-1*nc] * s3;
		s0 -= lptr[0-2*nc] * s2;
		s1 -= lptr[1-2*nc] * s2;
		s0 -= lptr[0-3*nc] * s1;
		// store result
		xptr[-4] = s0;
		xptr[-3] = s1;
		xptr[-2] = s2;
		xptr[-1] = s3;
		// update pointers for next four rows
		lptr -= 4 + 4 * nc;
		xptr -= 4;
	}
	// process left over rows
	for ( i--; i >= 0; i-- ) {
		s0 = b[i];
		lptr = L[0] + i;
		for ( j = i + 1; j < n; j++ ) {
			s0 -= lptr[j*nc] * x[j];
		}
		x[i] = s0;
	}

#else

	int i, j, nc;
	const float *ptr;
	double sum;

	nc = L.getNumColumns();
	for ( i = n - 1; i >= 0; i-- ) {
		sum = b[i];
		ptr = L[0] + i;
		for ( j = i + 1; j < n; j++ ) {
			sum -= ptr[j*nc] * x[j];
		}
		x[i] = sum;
	}

#endif
}

/*
============
CSIMD_Generic::matX_LDLTFactor

  in-place factorization LDL' of the n * n sub-matrix of mat
  the reciprocal of the diagonal elements are stored in invDiag
============
*/
bool VPCALL CSIMD_Generic::matX_LDLTFactor( CMatXD &mat, CVecXD &invDiag, const int n ) {
#if 1

	int i, j, k, nc;
	float *v, *diag, *mptr;
	double s0, s1, s2, s3, sum, d;

	v = (float *) _allocafloat16( n * sizeof( float ) );
	diag = (float *) _allocafloat16( n * sizeof( float ) );

	nc = mat.getNumColumns();

	if ( n <= 0 ) {
		return true;
	}

	mptr = mat[0];

	sum = mptr[0];

	if ( sum == 0.0f ) {
		return false;
	}

	diag[0] = sum;
	invDiag[0] = d = 1.0f / sum;

	if ( n <= 1 ) {
		return true;
	}

	mptr = mat[0];
	for ( j = 1; j < n; j++ ) {
		mptr[j*nc+0] = ( mptr[j*nc+0] ) * d;
	}

	mptr = mat[1];

	v[0] = diag[0] * mptr[0]; s0 = v[0] * mptr[0];
	sum = mptr[1] - s0;

	if ( sum == 0.0f ) {
		return false;
	}

	mat[1][1] = sum;
	diag[1] = sum;
	invDiag[1] = d = 1.0f / sum;

	if ( n <= 2 ) {
		return true;
	}

	mptr = mat[0];
	for ( j = 2; j < n; j++ ) {
		mptr[j*nc+1] = ( mptr[j*nc+1] - v[0] * mptr[j*nc+0] ) * d;
	}

	mptr = mat[2];

	v[0] = diag[0] * mptr[0]; s0 = v[0] * mptr[0];
	v[1] = diag[1] * mptr[1]; s1 = v[1] * mptr[1];
	sum = mptr[2] - s0 - s1;

	if ( sum == 0.0f ) {
		return false;
	}

	mat[2][2] = sum;
	diag[2] = sum;
	invDiag[2] = d = 1.0f / sum;

	if ( n <= 3 ) {
		return true;
	}

	mptr = mat[0];
	for ( j = 3; j < n; j++ ) {
		mptr[j*nc+2] = ( mptr[j*nc+2] - v[0] * mptr[j*nc+0] - v[1] * mptr[j*nc+1] ) * d;
	}

	mptr = mat[3];

	v[0] = diag[0] * mptr[0]; s0 = v[0] * mptr[0];
	v[1] = diag[1] * mptr[1]; s1 = v[1] * mptr[1];
	v[2] = diag[2] * mptr[2]; s2 = v[2] * mptr[2];
	sum = mptr[3] - s0 - s1 - s2;

	if ( sum == 0.0f ) {
		return false;
	}

	mat[3][3] = sum;
	diag[3] = sum;
	invDiag[3] = d = 1.0f / sum;

	if ( n <= 4 ) {
		return true;
	}

	mptr = mat[0];
	for ( j = 4; j < n; j++ ) {
		mptr[j*nc+3] = ( mptr[j*nc+3] - v[0] * mptr[j*nc+0] - v[1] * mptr[j*nc+1] - v[2] * mptr[j*nc+2] ) * d;
	}

	for ( i = 4; i < n; i++ ) {

		mptr = mat[i];

		v[0] = diag[0] * mptr[0]; s0 = v[0] * mptr[0];
		v[1] = diag[1] * mptr[1]; s1 = v[1] * mptr[1];
		v[2] = diag[2] * mptr[2]; s2 = v[2] * mptr[2];
		v[3] = diag[3] * mptr[3]; s3 = v[3] * mptr[3];
		for ( k = 4; k < i-3; k += 4 ) {
			v[k+0] = diag[k+0] * mptr[k+0]; s0 += v[k+0] * mptr[k+0];
			v[k+1] = diag[k+1] * mptr[k+1]; s1 += v[k+1] * mptr[k+1];
			v[k+2] = diag[k+2] * mptr[k+2]; s2 += v[k+2] * mptr[k+2];
			v[k+3] = diag[k+3] * mptr[k+3]; s3 += v[k+3] * mptr[k+3];
		}
		switch( i - k ) {
			NODEFAULT;
			case 3: v[k+2] = diag[k+2] * mptr[k+2]; s0 += v[k+2] * mptr[k+2];
			case 2: v[k+1] = diag[k+1] * mptr[k+1]; s1 += v[k+1] * mptr[k+1];
			case 1: v[k+0] = diag[k+0] * mptr[k+0]; s2 += v[k+0] * mptr[k+0];
			case 0: break;
		}
		sum = s3;
		sum += s2;
		sum += s1;
		sum += s0;
		sum = mptr[i] - sum;

		if ( sum == 0.0f ) {
			return false;
		}

		mat[i][i] = sum;
		diag[i] = sum;
		invDiag[i] = d = 1.0f / sum;

		if ( i + 1 >= n ) {
			return true;
		}

		mptr = mat[i+1];
		for ( j = i+1; j < n; j++ ) {
			s0 = mptr[0] * v[0];
			s1 = mptr[1] * v[1];
			s2 = mptr[2] * v[2];
			s3 = mptr[3] * v[3];
			for ( k = 4; k < i-7; k += 8 ) {
				s0 += mptr[k+0] * v[k+0];
				s1 += mptr[k+1] * v[k+1];
				s2 += mptr[k+2] * v[k+2];
				s3 += mptr[k+3] * v[k+3];
				s0 += mptr[k+4] * v[k+4];
				s1 += mptr[k+5] * v[k+5];
				s2 += mptr[k+6] * v[k+6];
				s3 += mptr[k+7] * v[k+7];
			}
			switch( i - k ) {
				NODEFAULT;
				case 7: s0 += mptr[k+6] * v[k+6];
				case 6: s1 += mptr[k+5] * v[k+5];
				case 5: s2 += mptr[k+4] * v[k+4];
				case 4: s3 += mptr[k+3] * v[k+3];
				case 3: s0 += mptr[k+2] * v[k+2];
				case 2: s1 += mptr[k+1] * v[k+1];
				case 1: s2 += mptr[k+0] * v[k+0];
				case 0: break;
			}
			sum = s3;
			sum += s2;
			sum += s1;
			sum += s0;
			mptr[i] = ( mptr[i] - sum ) * d;
			mptr += nc;
		}
	}

	return true;

#else

	int i, j, k, nc;
	float *v, *ptr, *diagPtr;
	double d, sum;

	v = (float *) _alloca16( n * sizeof( float ) );
	nc = mat.getNumColumns();

	for ( i = 0; i < n; i++ ) {

		ptr = mat[i];
		diagPtr = mat[0];
		sum = ptr[i];
		for ( j = 0; j < i; j++ ) {
			d = ptr[j];
		    v[j] = diagPtr[0] * d;
		    sum -= v[j] * d;
			diagPtr += nc + 1;
		}

		if ( sum == 0.0f ) {
			return false;
		}

		diagPtr[0] = sum;
		invDiag[i] = d = 1.0f / sum;

		if ( i + 1 >= n ) {
			continue;
		}

		ptr = mat[i+1];
		for ( j = i + 1; j < n; j++ ) {
			sum = ptr[i];
			for ( k = 0; k < i; k++ ) {
				sum -= ptr[k] * v[k];
			}
			ptr[i] = sum * d;
			ptr += nc;
		}
	}

	return true;

#endif
}

/*
============
CSIMD_Generic::blendJoints
============
*/
void VPCALL CSIMD_Generic::blendJoints( CJointQuaternion *joints, const CJointQuaternion *blendJoints, const float lerp, const int *index, const int numJoints ) {
	int i;

	for ( i = 0; i < numJoints; i++ ) {
		int j = index[i];
		joints[j].q.slerp( joints[j].q, blendJoints[j].q, lerp );
		joints[j].t.lerp( joints[j].t, blendJoints[j].t, lerp );
	}
}

/*
============
CSIMD_Generic::convertJointQuatsToJointMats
============
*/
void VPCALL CSIMD_Generic::convertJointQuatsToJointMats( CMatJoint3x4 *jointMats, const CJointQuaternion *jointQuats, const int numJoints ) {
	int i;

	for ( i = 0; i < numJoints; i++ ) {
		jointMats[i].setRotation( jointQuats[i].q.toMat3() );
		jointMats[i].setTranslation( jointQuats[i].t );
	}
}

/*
============
CSIMD_Generic::convertJointMatsToJointQuats
============
*/
void VPCALL CSIMD_Generic::convertJointMatsToJointQuats( CJointQuaternion *jointQuats, const CMatJoint3x4 *jointMats, const int numJoints ) {
	int i;

	for ( i = 0; i < numJoints; i++ ) {
		jointQuats[i] = jointMats[i].ToJointQuat();
	}
}

/*
============
CSIMD_Generic::transformJoints
============
*/
void VPCALL CSIMD_Generic::transformJoints( CMatJoint3x4 *jointMats, const int *parents, const int firstJoint, const int lastJoint ) {
	int i;

	for( i = firstJoint; i <= lastJoint; i++ ) {
		SMF_ASSERT( parents[i] < i );
		jointMats[i] *= jointMats[parents[i]];
	}
}

/*
============
CSIMD_Generic::untransformJoints
============
*/
void VPCALL CSIMD_Generic::untransformJoints( CMatJoint3x4 *jointMats, const int *parents, const int firstJoint, const int lastJoint ) {
	int i;

	for( i = lastJoint; i >= firstJoint; i-- ) {
		SMF_ASSERT( parents[i] < i );
		jointMats[i] /= jointMats[parents[i]];
	}
}

/*
============
CSIMD_Generic::transformVerts
============
*/
void VPCALL CSIMD_Generic::transformVerts( CVertex *verts, const int numVerts, const CMatJoint3x4 *joints, const CVec4D *weights, const int *index, int numWeights ) {
	int i, j;
	const sf_u8 *jointsPtr = (sf_u8 *)joints;

	for( j = i = 0; i < numVerts; i++ ) {
		CVec3D v;

		v = ( *(CMatJoint3x4 *) ( jointsPtr + index[j*2+0] ) ) * weights[j];
		while( index[j*2+1] == 0 ) {
			j++;
			v += ( *(CMatJoint3x4 *) ( jointsPtr + index[j*2+0] ) ) * weights[j];
		}
		j++;

		verts[i].xyz = v;
	}
}

/*
============
CSIMD_Generic::tracePointCull
============
*/
void VPCALL CSIMD_Generic::tracePointCull( sf_u8 *cullBits, sf_u8 &totalOr, const float radius, const CPlane *planes, const CVertex *verts, const int numVerts ) {
	int i;
	sf_u8 tOr;

	tOr = 0;

	for ( i = 0; i < numVerts; i++ ) {
		sf_u8 bits;
		float d0, d1, d2, d3, t;
		const CVec3D &v = verts[i].xyz;

		d0 = planes[0].getDistance( v );
		d1 = planes[1].getDistance( v );
		d2 = planes[2].getDistance( v );
		d3 = planes[3].getDistance( v );

		t = d0 + radius;
		bits  = FLOATSIGNBITSET( t ) << 0;
		t = d1 + radius;
		bits |= FLOATSIGNBITSET( t ) << 1;
		t = d2 + radius;
		bits |= FLOATSIGNBITSET( t ) << 2;
		t = d3 + radius;
		bits |= FLOATSIGNBITSET( t ) << 3;

		t = d0 - radius;
		bits |= FLOATSIGNBITSET( t ) << 4;
		t = d1 - radius;
		bits |= FLOATSIGNBITSET( t ) << 5;
		t = d2 - radius;
		bits |= FLOATSIGNBITSET( t ) << 6;
		t = d3 - radius;
		bits |= FLOATSIGNBITSET( t ) << 7;

		bits ^= 0x0F;		// flip lower four bits

		tOr |= bits;
		cullBits[i] = bits;
	}

	totalOr = tOr;
}

/*
============
CSIMD_Generic::decalPointCull
============
*/
void VPCALL CSIMD_Generic::decalPointCull( sf_u8 *cullBits, const CPlane *planes, const CVertex *verts, const int numVerts ) {
	int i;

	for ( i = 0; i < numVerts; i++ ) {
		sf_u8 bits;
		float d0, d1, d2, d3, d4, d5;
		const CVec3D &v = verts[i].xyz;

		d0 = planes[0].getDistance( v );
		d1 = planes[1].getDistance( v );
		d2 = planes[2].getDistance( v );
		d3 = planes[3].getDistance( v );
		d4 = planes[4].getDistance( v );
		d5 = planes[5].getDistance( v );

		bits  = FLOATSIGNBITSET( d0 ) << 0;
		bits |= FLOATSIGNBITSET( d1 ) << 1;
		bits |= FLOATSIGNBITSET( d2 ) << 2;
		bits |= FLOATSIGNBITSET( d3 ) << 3;
		bits |= FLOATSIGNBITSET( d4 ) << 4;
		bits |= FLOATSIGNBITSET( d5 ) << 5;

		cullBits[i] = bits ^ 0x3F;		// flip lower 6 bits
	}
}

/*
============
CSIMD_Generic::overlayPointCull
============
*/
void VPCALL CSIMD_Generic::overlayPointCull( sf_u8 *cullBits, CVec2D *texCoords, const CPlane *planes, const CVertex *verts, const int numVerts ) {
	int i;

	for ( i = 0; i < numVerts; i++ ) {
		sf_u8 bits;
		float d0, d1;
		const CVec3D &v = verts[i].xyz;

		texCoords[i][0] = d0 = planes[0].getDistance( v );
		texCoords[i][1] = d1 = planes[1].getDistance( v );

		bits  = FLOATSIGNBITSET( d0 ) << 0;
		d0 = 1.0f - d0;
		bits |= FLOATSIGNBITSET( d1 ) << 1;
		d1 = 1.0f - d1;
		bits |= FLOATSIGNBITSET( d0 ) << 2;
		bits |= FLOATSIGNBITSET( d1 ) << 3;

		cullBits[i] = bits;
	}
}

/*
============
CSIMD_Generic::deriveTriPlanes

	Derives a plane equation for each triangle.
============
*/
void VPCALL CSIMD_Generic::deriveTriPlanes( CPlane *planes, const CVertex *verts, const int numVerts, const int *indexes, const int numIndexes ) {
	int i;

	for ( i = 0; i < numIndexes; i += 3 ) {
		const CVertex *a, *b, *c;
		float d0[3], d1[3], f;
		CVec3D n;

		a = verts + indexes[i + 0];
		b = verts + indexes[i + 1];
		c = verts + indexes[i + 2];

		d0[0] = b->xyz[0] - a->xyz[0];
		d0[1] = b->xyz[1] - a->xyz[1];
		d0[2] = b->xyz[2] - a->xyz[2];

		d1[0] = c->xyz[0] - a->xyz[0];
		d1[1] = c->xyz[1] - a->xyz[1];
		d1[2] = c->xyz[2] - a->xyz[2];

		n[0] = d1[1] * d0[2] - d1[2] * d0[1];
		n[1] = d1[2] * d0[0] - d1[0] * d0[2];
		n[2] = d1[0] * d0[1] - d1[1] * d0[0];

		f = CMath::rSqrt( n.x * n.x + n.y * n.y + n.z * n.z );

		n.x *= f;
		n.y *= f;
		n.z *= f;

		planes->setNormal( n );
		planes->fitThroughPoint( a->xyz );
		planes++;
	}
}

/*
============
CSIMD_Generic:: deriveTangents

	Derives the normal and orthogonal tangent vectors for the triangle vertices.
	For each vertex the normal and tangent vectors are derived from all triangles
	using the vertex which results in smooth tangents across the mesh.
	In the process the triangle planes are calculated as well.
============
*/
void VPCALL CSIMD_Generic:: deriveTangents( CPlane *planes, CVertex *verts, const int numVerts, const int *indexes, const int numIndexes ) {
	int i;

	bool *used = (bool *)_allocafloat16( numVerts * sizeof( used[0] ) );
	memset( used, 0, numVerts * sizeof( used[0] ) );

	CPlane *planesPtr = planes;
	for ( i = 0; i < numIndexes; i += 3 ) {
		CVertex *a, *b, *c;
		unsigned long signBit;
		float d0[5], d1[5], f, area;
		CVec3D n, t0, t1;

		int v0 = indexes[i + 0];
		int v1 = indexes[i + 1];
		int v2 = indexes[i + 2];

		a = verts + v0;
		b = verts + v1;
		c = verts + v2;

		d0[0] = b->xyz[0] - a->xyz[0];
		d0[1] = b->xyz[1] - a->xyz[1];
		d0[2] = b->xyz[2] - a->xyz[2];
		d0[3] = b->st[0] - a->st[0];
		d0[4] = b->st[1] - a->st[1];

		d1[0] = c->xyz[0] - a->xyz[0];
		d1[1] = c->xyz[1] - a->xyz[1];
		d1[2] = c->xyz[2] - a->xyz[2];
		d1[3] = c->st[0] - a->st[0];
		d1[4] = c->st[1] - a->st[1];

		// normal
		n[0] = d1[1] * d0[2] - d1[2] * d0[1];
		n[1] = d1[2] * d0[0] - d1[0] * d0[2];
		n[2] = d1[0] * d0[1] - d1[1] * d0[0];

		f = CMath::rSqrt( n.x * n.x + n.y * n.y + n.z * n.z );

		n.x *= f;
		n.y *= f;
		n.z *= f;

		planesPtr->setNormal( n );
		planesPtr->fitThroughPoint( a->xyz );
		planesPtr++;

		// area sign bit
		area = d0[3] * d1[4] - d0[4] * d1[3];
		signBit = ( *(unsigned long *)&area ) & ( 1 << 31 );

		// first tangent
		t0[0] = d0[0] * d1[4] - d0[4] * d1[0];
		t0[1] = d0[1] * d1[4] - d0[4] * d1[1];
		t0[2] = d0[2] * d1[4] - d0[4] * d1[2];

		f = CMath::rSqrt( t0.x * t0.x + t0.y * t0.y + t0.z * t0.z );
		*(unsigned long *)&f ^= signBit;

		t0.x *= f;
		t0.y *= f;
		t0.z *= f;

		// second tangent
		t1[0] = d0[3] * d1[0] - d0[0] * d1[3];
		t1[1] = d0[3] * d1[1] - d0[1] * d1[3];
		t1[2] = d0[3] * d1[2] - d0[2] * d1[3];

		f = CMath::rSqrt( t1.x * t1.x + t1.y * t1.y + t1.z * t1.z );
		*(unsigned long *)&f ^= signBit;

		t1.x *= f;
		t1.y *= f;
		t1.z *= f;

		if ( used[v0] ) {
			a->normal += n;
			a->tangents[0] += t0;
			a->tangents[1] += t1;
		} else {
			a->normal = n;
			a->tangents[0] = t0;
			a->tangents[1] = t1;
			used[v0] = true;
		}

		if ( used[v1] ) {
			b->normal += n;
			b->tangents[0] += t0;
			b->tangents[1] += t1;
		} else {
			b->normal = n;
			b->tangents[0] = t0;
			b->tangents[1] = t1;
			used[v1] = true;
		}

		if ( used[v2] ) {
			c->normal += n;
			c->tangents[0] += t0;
			c->tangents[1] += t1;
		} else {
			c->normal = n;
			c->tangents[0] = t0;
			c->tangents[1] = t1;
			used[v2] = true;
		}
	}
}

/*
============
CSIMD_Generic::deriveUnsmoothedTangents

	Derives the normal and orthogonal tangent vectors for the triangle vertices.
	For each vertex the normal and tangent vectors are derived from a single dominant triangle.
============
*/
#define DERIVE_UNSMOOTHED_BITANGENT
/* s
void VPCALL CSIMD_Generic::deriveUnsmoothedTangents( CVertex *verts, const dominantTri_s *dominantTris, const int numVerts ) {
	int i;

	for ( i = 0; i < numVerts; i++ ) {
		CVertex *a, *b, *c;
		float d0, d1, d2, d3, d4;
		float d5, d6, d7, d8, d9;
		float s0, s1, s2;
		float n0, n1, n2;
		float t0, t1, t2;
		float t3, t4, t5;

		const dominantTri_s &dt = dominantTris[i];

		a = verts + i;
		b = verts + dt.v2;
		c = verts + dt.v3;

		d0 = b->xyz[0] - a->xyz[0];
		d1 = b->xyz[1] - a->xyz[1];
		d2 = b->xyz[2] - a->xyz[2];
		d3 = b->st[0] - a->st[0];
		d4 = b->st[1] - a->st[1];

		d5 = c->xyz[0] - a->xyz[0];
		d6 = c->xyz[1] - a->xyz[1];
		d7 = c->xyz[2] - a->xyz[2];
		d8 = c->st[0] - a->st[0];
		d9 = c->st[1] - a->st[1];

		s0 = dt.normalizationScale[0];
		s1 = dt.normalizationScale[1];
		s2 = dt.normalizationScale[2];

		n0 = s2 * ( d6 * d2 - d7 * d1 );
		n1 = s2 * ( d7 * d0 - d5 * d2 );
		n2 = s2 * ( d5 * d1 - d6 * d0 );

		t0 = s0 * ( d0 * d9 - d4 * d5 );
		t1 = s0 * ( d1 * d9 - d4 * d6 );
		t2 = s0 * ( d2 * d9 - d4 * d7 );

#ifndef DERIVE_UNSMOOTHED_BITANGENT
		t3 = s1 * ( d3 * d5 - d0 * d8 );
		t4 = s1 * ( d3 * d6 - d1 * d8 );
		t5 = s1 * ( d3 * d7 - d2 * d8 );
#else
		t3 = s1 * ( n2 * t1 - n1 * t2 );
		t4 = s1 * ( n0 * t2 - n2 * t0 );
		t5 = s1 * ( n1 * t0 - n0 * t1 );
#endif

		a->normal[0] = n0;
		a->normal[1] = n1;
		a->normal[2] = n2;

		a->tangents[0][0] = t0;
		a->tangents[0][1] = t1;
		a->tangents[0][2] = t2;

		a->tangents[1][0] = t3;
		a->tangents[1][1] = t4;
		a->tangents[1][2] = t5;
	}
}
*/
/*
============
CSIMD_Generic:: normalizeTangents

	Normalizes each vertex normal and projects and normalizes the
	tangent vectors onto the plane orthogonal to the vertex normal.
============
*/

void VPCALL CSIMD_Generic:: normalizeTangents( CVertex *verts, const int numVerts ) {

	for ( int i = 0; i < numVerts; i++ ) {
		CVec3D &v = verts[i].normal;
		float f;

		f = CMath::rSqrt( v.x * v.x + v.y * v.y + v.z * v.z );
		v.x *= f; v.y *= f; v.z *= f;

		for ( int j = 0; j < 2; j++ ) {
			CVec3D &t = verts[i].tangents[j];

			t -= ( t * v ) * v;
			f = CMath::rSqrt( t.x * t.x + t.y * t.y + t.z * t.z );
			t.x *= f; t.y *= f; t.z *= f;
		}
	}
}

/*
============
CSIMD_Generic:: createTextureSpaceLightVectors

	Calculates light vectors in texture space for the given triangle vertices.
	For each vertex the direction towards the light origin is projected onto texture space.
	The light vectors are only calculated for the vertices referenced by the indexes.
============
*/
void VPCALL CSIMD_Generic:: createTextureSpaceLightVectors( CVec3D *lightVectors, const CVec3D &lightOrigin, const CVertex *verts, const int numVerts, const int *indexes, const int numIndexes ) {

	bool *used = (bool *)_allocafloat16( numVerts * sizeof( used[0] ) );
	memset( used, 0, numVerts * sizeof( used[0] ) );

	for ( int i = numIndexes - 1; i >= 0; i-- ) {
		used[indexes[i]] = true;
	}

	for ( int i = 0; i < numVerts; i++ ) {
		if ( !used[i] ) {
			continue;
		}

		const CVertex *v = &verts[i];

		CVec3D lightDir = lightOrigin - v->xyz;

		lightVectors[i][0] = lightDir * v->tangents[0];
		lightVectors[i][1] = lightDir * v->tangents[1];
		lightVectors[i][2] = lightDir * v->normal;
	}
}

/*
============
CSIMD_Generic:: createSpecularTextureCoords

	Calculates specular texture coordinates for the given triangle vertices.
	For each vertex the normalized direction towards the light origin is added to the
	normalized direction towards the view origin and the result is projected onto texture space.
	The texture coordinates are only calculated for the vertices referenced by the indexes.
============
*/
void VPCALL CSIMD_Generic:: createSpecularTextureCoords( CVec4D *texCoords, const CVec3D &lightOrigin, const CVec3D &viewOrigin, const CVertex *verts, const int numVerts, const int *indexes, const int numIndexes ) {

	bool *used = (bool *)_allocafloat16( numVerts * sizeof( used[0] ) );
	memset( used, 0, numVerts * sizeof( used[0] ) );

	for ( int i = numIndexes - 1; i >= 0; i-- ) {
		used[indexes[i]] = true;
	}

	for ( int i = 0; i < numVerts; i++ ) {
		if ( !used[i] ) {
			continue;
		}

		const CVertex *v = &verts[i];

		CVec3D lightDir = lightOrigin - v->xyz;
		CVec3D viewDir = viewOrigin - v->xyz;

		float ilength;

		ilength = CMath::rSqrt( lightDir * lightDir );
		lightDir[0] *= ilength;
		lightDir[1] *= ilength;
		lightDir[2] *= ilength;

		ilength = CMath::rSqrt( viewDir * viewDir );
		viewDir[0] *= ilength;
		viewDir[1] *= ilength;
		viewDir[2] *= ilength;

		lightDir += viewDir;

		texCoords[i][0] = lightDir * v->tangents[0];
		texCoords[i][1] = lightDir * v->tangents[1];
		texCoords[i][2] = lightDir * v->normal;
		texCoords[i][3] = 1.0f;
	}
}

/*
============
CSIMD_Generic:: createShadowCache
============
*/
int VPCALL CSIMD_Generic:: createShadowCache( CVec4D *vertexCache, int *vertRemap, const CVec3D &lightOrigin, const CVertex *verts, const int numVerts ) {
	int outVerts = 0;

	for ( int i = 0; i < numVerts; i++ ) {
		if ( vertRemap[i] ) {
			continue;
		}
		const float *v = verts[i].xyz.toFloatPtr();
		vertexCache[outVerts+0][0] = v[0];
		vertexCache[outVerts+0][1] = v[1];
		vertexCache[outVerts+0][2] = v[2];
		vertexCache[outVerts+0][3] = 1.0f;

		// R_SetupProjection() builds the projection matrix with a slight crunch
		// for depth, which keeps this w=0 division from rasterizing right at the
		// wrap around point and causing depth fighting with the rear caps
		vertexCache[outVerts+1][0] = v[0] - lightOrigin[0];
		vertexCache[outVerts+1][1] = v[1] - lightOrigin[1];
		vertexCache[outVerts+1][2] = v[2] - lightOrigin[2];
		vertexCache[outVerts+1][3] = 0.0f;
		vertRemap[i] = outVerts;
		outVerts += 2;
	}
	return outVerts;
}

/*
============
CSIMD_Generic:: createVertexProgramShadowCache
============
*/
int VPCALL CSIMD_Generic:: createVertexProgramShadowCache( CVec4D *vertexCache, const CVertex *verts, const int numVerts ) {
	for ( int i = 0; i < numVerts; i++ ) {
		const float *v = verts[i].xyz.toFloatPtr();
		vertexCache[i*2+0][0] = v[0];
		vertexCache[i*2+1][0] = v[0];
		vertexCache[i*2+0][1] = v[1];
		vertexCache[i*2+1][1] = v[1];
		vertexCache[i*2+0][2] = v[2];
		vertexCache[i*2+1][2] = v[2];
		vertexCache[i*2+0][3] = 1.0f;
		vertexCache[i*2+1][3] = 0.0f;
	}
	return numVerts * 2;
}

/*
============
CSIMD_Generic:: upSamplePCMTo44kHz

  Duplicate samples for 44kHz output.
============
*/
void CSIMD_Generic:: upSamplePCMTo44kHz( float *dest, const short *src, const int numSamples, const int kHz, const int numChannels ) {
	if ( kHz == 11025 ) {
		if ( numChannels == 1 ) {
			for ( int i = 0; i < numSamples; i++ ) {
				dest[i*4+0] = dest[i*4+1] = dest[i*4+2] = dest[i*4+3] = (float) src[i+0];
			}
		} else {
			for ( int i = 0; i < numSamples; i += 2 ) {
				dest[i*4+0] = dest[i*4+2] = dest[i*4+4] = dest[i*4+6] = (float) src[i+0];
				dest[i*4+1] = dest[i*4+3] = dest[i*4+5] = dest[i*4+7] = (float) src[i+1];
			}
		}
	} else if ( kHz == 22050 ) {
		if ( numChannels == 1 ) {
			for ( int i = 0; i < numSamples; i++ ) {
				dest[i*2+0] = dest[i*2+1] = (float) src[i+0];
			}
		} else {
			for ( int i = 0; i < numSamples; i += 2 ) {
				dest[i*2+0] = dest[i*2+2] = (float) src[i+0];
				dest[i*2+1] = dest[i*2+3] = (float) src[i+1];
			}
		}
	} else if ( kHz == 44100 ) {
		for ( int i = 0; i < numSamples; i++ ) {
			dest[i] = (float) src[i];
		}
	} else {
		SMF_ASSERT( 0 );
	}
}

/*
============
CSIMD_Generic:: upSampleOGGTo44kHz

  Duplicate samples for 44kHz output.
============
*/
void CSIMD_Generic:: upSampleOGGTo44kHz( float *dest, const float * const *ogg, const int numSamples, const int kHz, const int numChannels ) {
	if ( kHz == 11025 ) {
		if ( numChannels == 1 ) {
			for ( int i = 0; i < numSamples; i++ ) {
				dest[i*4+0] = dest[i*4+1] = dest[i*4+2] = dest[i*4+3] = ogg[0][i] * 32768.0f;
			}
		} else {
			for ( int i = 0; i < numSamples >> 1; i++ ) {
				dest[i*8+0] = dest[i*8+2] = dest[i*8+4] = dest[i*8+6] = ogg[0][i] * 32768.0f;
				dest[i*8+1] = dest[i*8+3] = dest[i*8+5] = dest[i*8+7] = ogg[1][i] * 32768.0f;
			}
		}
	} else if ( kHz == 22050 ) {
		if ( numChannels == 1 ) {
			for ( int i = 0; i < numSamples; i++ ) {
				dest[i*2+0] = dest[i*2+1] = ogg[0][i] * 32768.0f;
			}
		} else {
			for ( int i = 0; i < numSamples >> 1; i++ ) {
				dest[i*4+0] = dest[i*4+2] = ogg[0][i] * 32768.0f;
				dest[i*4+1] = dest[i*4+3] = ogg[1][i] * 32768.0f;
			}
		}
	} else if ( kHz == 44100 ) {
		if ( numChannels == 1 ) {
			for ( int i = 0; i < numSamples; i++ ) {
				dest[i*1+0] = ogg[0][i] * 32768.0f;
			}
		} else {
			for ( int i = 0; i < numSamples >> 1; i++ ) {
				dest[i*2+0] = ogg[0][i] * 32768.0f;
				dest[i*2+1] = ogg[1][i] * 32768.0f;
			}
		}
	} else {
		SMF_ASSERT( 0 );
	}
}

/*
============
CSIMD_Generic:: mixSoundTwoSpeakerMono
============
*/
void VPCALL CSIMD_Generic:: mixSoundTwoSpeakerMono( float *mixBuffer, const float *samples, const int numSamples, const float lastV[2], const float currentV[2] ) {
	float sL = lastV[0];
	float sR = lastV[1];
	float incL = ( currentV[0] - lastV[0] ) / MIXBUFFER_SAMPLES;
	float incR = ( currentV[1] - lastV[1] ) / MIXBUFFER_SAMPLES;

	SMF_ASSERT( numSamples == MIXBUFFER_SAMPLES );

	for( int j = 0; j < MIXBUFFER_SAMPLES; j++ ) {
		mixBuffer[j*2+0] += samples[j] * sL;
		mixBuffer[j*2+1] += samples[j] * sR;
		sL += incL;
		sR += incR;
	}
}

/*
============
CSIMD_Generic:: mixSoundTwoSpeakerStereo
============
*/
void VPCALL CSIMD_Generic:: mixSoundTwoSpeakerStereo( float *mixBuffer, const float *samples, const int numSamples, const float lastV[2], const float currentV[2] ) {
	float sL = lastV[0];
	float sR = lastV[1];
	float incL = ( currentV[0] - lastV[0] ) / MIXBUFFER_SAMPLES;
	float incR = ( currentV[1] - lastV[1] ) / MIXBUFFER_SAMPLES;

	SMF_ASSERT( numSamples == MIXBUFFER_SAMPLES );

	for( int j = 0; j < MIXBUFFER_SAMPLES; j++ ) {
		mixBuffer[j*2+0] += samples[j*2+0] * sL;
		mixBuffer[j*2+1] += samples[j*2+1] * sR;
		sL += incL;
		sR += incR;
	}
}

/*
============
CSIMD_Generic:: mixSoundSixSpeakerMono
============
*/
void VPCALL CSIMD_Generic:: mixSoundSixSpeakerMono( float *mixBuffer, const float *samples, const int numSamples, const float lastV[6], const float currentV[6] ) {
	float sL0 = lastV[0];
	float sL1 = lastV[1];
	float sL2 = lastV[2];
	float sL3 = lastV[3];
	float sL4 = lastV[4];
	float sL5 = lastV[5];

	float incL0 = ( currentV[0] - lastV[0] ) / MIXBUFFER_SAMPLES;
	float incL1 = ( currentV[1] - lastV[1] ) / MIXBUFFER_SAMPLES;
	float incL2 = ( currentV[2] - lastV[2] ) / MIXBUFFER_SAMPLES;
	float incL3 = ( currentV[3] - lastV[3] ) / MIXBUFFER_SAMPLES;
	float incL4 = ( currentV[4] - lastV[4] ) / MIXBUFFER_SAMPLES;
	float incL5 = ( currentV[5] - lastV[5] ) / MIXBUFFER_SAMPLES;

	SMF_ASSERT( numSamples == MIXBUFFER_SAMPLES );

	for( int i = 0; i < MIXBUFFER_SAMPLES; i++ ) {
		mixBuffer[i*6+0] += samples[i] * sL0;
		mixBuffer[i*6+1] += samples[i] * sL1;
		mixBuffer[i*6+2] += samples[i] * sL2;
		mixBuffer[i*6+3] += samples[i] * sL3;
		mixBuffer[i*6+4] += samples[i] * sL4;
		mixBuffer[i*6+5] += samples[i] * sL5;
		sL0 += incL0;
		sL1 += incL1;
		sL2 += incL2;
		sL3 += incL3;
		sL4 += incL4;
		sL5 += incL5;
	}
}

/*
============
CSIMD_Generic:: mixSoundSixSpeakerStereo
============
*/
void VPCALL CSIMD_Generic:: mixSoundSixSpeakerStereo( float *mixBuffer, const float *samples, const int numSamples, const float lastV[6], const float currentV[6] ) {
	float sL0 = lastV[0];
	float sL1 = lastV[1];
	float sL2 = lastV[2];
	float sL3 = lastV[3];
	float sL4 = lastV[4];
	float sL5 = lastV[5];

	float incL0 = ( currentV[0] - lastV[0] ) / MIXBUFFER_SAMPLES;
	float incL1 = ( currentV[1] - lastV[1] ) / MIXBUFFER_SAMPLES;
	float incL2 = ( currentV[2] - lastV[2] ) / MIXBUFFER_SAMPLES;
	float incL3 = ( currentV[3] - lastV[3] ) / MIXBUFFER_SAMPLES;
	float incL4 = ( currentV[4] - lastV[4] ) / MIXBUFFER_SAMPLES;
	float incL5 = ( currentV[5] - lastV[5] ) / MIXBUFFER_SAMPLES;

	SMF_ASSERT( numSamples == MIXBUFFER_SAMPLES );

	for( int i = 0; i < MIXBUFFER_SAMPLES; i++ ) {
		mixBuffer[i*6+0] += samples[i*2+0] * sL0;
		mixBuffer[i*6+1] += samples[i*2+1] * sL1;
		mixBuffer[i*6+2] += samples[i*2+0] * sL2;
		mixBuffer[i*6+3] += samples[i*2+0] * sL3;
		mixBuffer[i*6+4] += samples[i*2+0] * sL4;
		mixBuffer[i*6+5] += samples[i*2+1] * sL5;
		sL0 += incL0;
		sL1 += incL1;
		sL2 += incL2;
		sL3 += incL3;
		sL4 += incL4;
		sL5 += incL5;
	}
}

/*
============
CSIMD_Generic:: mixedSoundToSamples
============
*/
void VPCALL CSIMD_Generic:: mixedSoundToSamples( short *samples, const float *mixBuffer, const int numSamples ) {

	for ( int i = 0; i < numSamples; i++ ) {
		if ( mixBuffer[i] <= -32768.0f ) {
			samples[i] = -32768;
		} else if ( mixBuffer[i] >= 32767.0f ) {
			samples[i] = 32767;
		} else {
			samples[i] = (short) mixBuffer[i];
		}
	}
}

//=================GNU C CODE==========================================



/*
	TODO:
	Optimize distance() -- It currently uses already implemented functions
*/


void VPCALL CSIMD_Generic::vector3D_Sum(CVec3D* pOut, const CVec3D* pIn)
{

	pOut->x += pIn->x;
	pOut->y += pIn->y;
	pOut->z += pIn->z;

}

void VPCALL CSIMD_Generic::vector3D_SumOf(CVec3D* pOut, const CVec3D* pIn1, const CVec3D* pIn2)
{
	pOut->x = pIn1->x + pIn2->x;
	pOut->y = pIn1->y + pIn2->y;
	pOut->z = pIn1->z + pIn2->z;

}


void VPCALL CSIMD_Generic::vector3D_Diff(CVec3D* pLeft, CVec3D* pRight)
{

	pLeft->x -= pRight->x;
	pLeft->y -= pRight->y;
	pLeft->z -= pRight->z;

}

void VPCALL CSIMD_Generic::vector3D_DiffOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight)
{

	pOut->x = pLeft->x - pRight->x;
	pOut->y = pLeft->y - pRight->y;
	pOut->z = pLeft->z - pRight->z;

}

void VPCALL CSIMD_Generic::vector3D_Scale(CVec3D* pOut, float scalar)
{

	pOut->x *= scalar;
	pOut->y *= scalar;
	pOut->z *= scalar;


}

void VPCALL CSIMD_Generic::vector3D_ScaleOf(CVec3D* pOut, const CVec3D* pIn, float scalar)
{

	pOut->x = pIn->x * scalar;
	pOut->y = pIn->y * scalar;
	pOut->z = pIn->z * scalar;


}

float  VPCALL CSIMD_Generic::vector3D_Dot(const CVec3D* pSrc1, const CVec3D* pSrc2)
{

	return (pSrc1->x * pSrc2->x) + (pSrc1->y * pSrc2->y) + (pSrc1->z * pSrc2->z);


}


float  VPCALL CSIMD_Generic::vector3D_LengthSq(const CVec3D* pVec)
{

	return (pVec->x*pVec->x) + (pVec->y*pVec->y) + (pVec->z*pVec->z);

}

float  VPCALL CSIMD_Generic::vector3D_Length(const CVec3D* pVec)
{

	return sqrtf((pVec->x*pVec->x) + (pVec->y*pVec->y) + (pVec->z*pVec->z));


}

void VPCALL CSIMD_Generic::vector3D_Cross(CVec3D* pLeft, const CVec3D* pRight)
{

	CVec3D tmp;

	tmp.x = pLeft->x;
	tmp.y = pLeft->y;
	tmp.z = pLeft->z;

	pLeft->x = (tmp.y * pRight->z) - (tmp.z * pRight->y);
	pLeft->y = (tmp.z * pRight->x) - (tmp.x * pRight->z);
	pLeft->z = (tmp.x * pRight->y) - (tmp.y * pRight->x);


}

void VPCALL CSIMD_Generic::vector3D_CrossOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight)
{

	pOut->x = (pLeft->y * pRight->z) - (pLeft->z * pRight->y);
	pOut->y = (pLeft->z * pRight->x) - (pLeft->x * pRight->z);
	pOut->z = (pLeft->x * pRight->y) - (pLeft->y * pRight->x);


}
void VPCALL CSIMD_Generic::vector3D_Normalize(CVec3D* pVec)
{

	float inv_mag = 1.0f / (sqrtf((pVec->x*pVec->x) + (pVec->y*pVec->y) + (pVec->z*pVec->z))+CMath::FLT_EPSILON);

	pVec->x *= inv_mag;
	pVec->y *= inv_mag;
	pVec->z *= inv_mag;

}
void VPCALL CSIMD_Generic::vector3D_NormalizeOf(CVec3D* pOut, const CVec3D* pVec)
{
	float inv_mag = 1.0f / (sqrtf((pVec->x*pVec->x) + (pVec->y*pVec->y) + (pVec->z*pVec->z))+CMath::FLT_EPSILON);

	pOut->x = pVec->x * inv_mag;
	pOut->y = pVec->y * inv_mag;
	pOut->z = pVec->z * inv_mag;

}

float  VPCALL CSIMD_Generic::vector3D_Distance(const CVec3D* pVec1, const CVec3D* pVec2)
{
	/* TODO: Optimize me completely */
	CVec3D diff;
	vector3D_DiffOf(&diff, pVec1, pVec2);
	return vector3D_Length(&diff);
}

/*
	TODO:
	Optimize distance() -- It currently uses already implemented functions
	Optimize 3DNow! cross product -- Too many memory references
*/

void VPCALL CSIMD_Generic::vector3D_AlignedSum(CVec3D* pOut, const CVec3D* pIn)
{

	pOut->x += pIn->x;
	pOut->y += pIn->y;
	pOut->z += pIn->z;

}

void VPCALL CSIMD_Generic::vector3D_AlignedSumOf(CVec3D* pOut, const CVec3D* pIn1, const CVec3D* pIn2)
{
	pOut->x = pIn1->x + pIn2->x;
	pOut->y = pIn1->y + pIn2->y;
	pOut->z = pIn1->z + pIn2->z;

}


void VPCALL CSIMD_Generic::vector3D_AlignedDiff(CVec3D* pLeft, CVec3D* pRight)
{

	pLeft->x -= pRight->x;
	pLeft->y -= pRight->y;
	pLeft->z -= pRight->z;

}

void VPCALL CSIMD_Generic::vector3D_AlignedDiffOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight)
{

	pOut->x = pLeft->x - pRight->x;
	pOut->y = pLeft->y - pRight->y;
	pOut->z = pLeft->z - pRight->z;

}

void VPCALL CSIMD_Generic::vector3D_AlignedScale(CVec3D* pOut, float scalar)
{

	pOut->x *= scalar;
	pOut->y *= scalar;
	pOut->z *= scalar;

}

void VPCALL CSIMD_Generic::vector3D_AlignedScaleOf(CVec3D* pOut, const CVec3D* pIn, float scalar)
{

	pOut->x = pIn->x * scalar;
	pOut->y = pIn->y * scalar;
	pOut->z = pIn->z * scalar;


}

float  VPCALL CSIMD_Generic::vector3D_AlignedDot(const CVec3D* pSrc1, const CVec3D* pSrc2)
{

	return (pSrc1->x * pSrc2->x) + (pSrc1->y * pSrc2->y) + (pSrc1->z * pSrc2->z);


}



float  VPCALL CSIMD_Generic::vector3D_AlignedLengthSq(const CVec3D* pVec)
{

	return (pVec->x*pVec->x) + (pVec->y*pVec->y) + (pVec->z*pVec->z);

}

float  VPCALL CSIMD_Generic::vector3D_AlignedLength(const CVec3D* pVec)
{


	return sqrtf((pVec->x*pVec->x) + (pVec->y*pVec->y) + (pVec->z*pVec->z));


}

void VPCALL CSIMD_Generic::vector3D_AlignedCross(CVec3D* pLeft, const CVec3D* pRight)
{

	CVec3D tmp;

	tmp.x = pLeft->x;
	tmp.y = pLeft->y;
	tmp.z = pLeft->z;

	pLeft->x = (tmp.y * pRight->z) - (tmp.z * pRight->y);
	pLeft->y = (tmp.z * pRight->x) - (tmp.x * pRight->z);
	pLeft->z = (tmp.x * pRight->y) - (tmp.y * pRight->x);


}

void VPCALL CSIMD_Generic::vector3D_AlignedCrossOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight)
{

	pOut->x = (pLeft->y * pRight->z) - (pLeft->z * pRight->y);
	pOut->y = (pLeft->z * pRight->x) - (pLeft->x * pRight->z);
	pOut->z = (pLeft->x * pRight->y) - (pLeft->y * pRight->x);


}
void VPCALL CSIMD_Generic::vector3D_AlignedNormalize(CVec3D* pVec)
{

	float inv_mag = 1.0f / (sqrtf((pVec->x*pVec->x) + (pVec->y*pVec->y) + (pVec->z*pVec->z))+CMath::FLT_EPSILON);

	pVec->x *= inv_mag;
	pVec->y *= inv_mag;
	pVec->z *= inv_mag;

}
void VPCALL CSIMD_Generic::vector3D_AlignedNormalizeOf(CVec3D* pOut, const CVec3D* pVec)
{
	float inv_mag = 1.0f / (sqrtf((pVec->x*pVec->x) + (pVec->y*pVec->y) + (pVec->z*pVec->z))+CMath::FLT_EPSILON);

	pOut->x = pVec->x * inv_mag;
	pOut->y = pVec->y * inv_mag;
	pOut->z = pVec->z * inv_mag;

}

float  VPCALL CSIMD_Generic::vector3D_AlignedDistance(const CVec3D* pVec1, const CVec3D* pVec2)
{


	/* TODO: Optimize me completely */
	CVec3D diff;
	vector3D_DiffOf(&diff, pVec1, pVec2);
	return vector3D_Length(&diff);

}

//======================  CVec4D============================

void VPCALL CSIMD_Generic::vector4D_Sum(CVec4D* pOut, const CVec4D* pIn)
{
	pOut->x += pIn->x;
	pOut->y += pIn->y;
	pOut->z += pIn->z;
	pOut->w += pIn->w;
}

void VPCALL CSIMD_Generic::vector4D_SumOf(CVec4D* pOut, const CVec4D* pIn1, const CVec4D* pIn2)
{
	pOut->x = pIn1->x + pIn2->x;
	pOut->y = pIn1->y + pIn2->y;
	pOut->z = pIn1->z + pIn2->z;
	pOut->w = pIn1->w + pIn2->w;
}


//Subtracts pIn from pOut and stores the result in pOut. POut = pIn-pOut
void VPCALL CSIMD_Generic::vector4D_Diff(CVec4D* pOut, CVec4D* pIn)
{
	pOut->x -= pIn->x;
	pOut->y -= pIn->y;
	pOut->z -= pIn->z;
	pOut->w -= pIn->w;
}


//Subtracts pIn2 from pIn1 and stores the result in pOut.  pout = pIn1-Pin2
void VPCALL CSIMD_Generic::vector4D_DiffOf(CVec4D* pOut, const CVec4D* pIn1, const CVec4D* pIn2)
{
	pOut->x = pIn1->x - pIn2->x;
	pOut->y = pIn1->y - pIn2->y;
	pOut->z = pIn1->z - pIn2->z;
	pOut->w = pIn1->w - pIn2->w;
}
//Scales the components of pOut by scalar and stores the result in pOut.
void VPCALL CSIMD_Generic::vector4D_Scale(CVec4D* pOut, float scalar)
{
	pOut->x *= scalar;
	pOut->y *= scalar;
	pOut->z *= scalar;
	pOut->w *= scalar;

}
//Scales the components of pIn by scalar and stores the result in pOut.
void VPCALL CSIMD_Generic::vector4D_ScaleOf(CVec4D* pOut, const CVec4D* pIn, float scalar)
{
	pOut->x = pIn->x * scalar;
	pOut->y = pIn->y * scalar;
	pOut->z = pIn->z * scalar;
	pOut->w = pIn->w * scalar;

}
//Calculates the dot product (known as the inner product) of pSrc1 and pSrc2, that represents the cosine of the angle between them, if both are normalized.

float  VPCALL CSIMD_Generic::vector4D_Dot(const CVec4D* pSrc1, const CVec4D* pSrc2)
{
	return (pSrc1->x * pSrc2->x) + (pSrc1->y * pSrc2->y) + (pSrc1->z * pSrc2->z)+ (pSrc1->w * pSrc2->w);

}

float  VPCALL CSIMD_Generic::vector4D_LengthSq(const CVec4D* pVec)
{
		return (pVec->x*pVec->x) + (pVec->y*pVec->y) + (pVec->z*pVec->z)+ (pVec->w*pVec->w);
}

float  VPCALL CSIMD_Generic::vector4D_Length(const CVec4D* pVec)
{
	return sqrtf((pVec->x*pVec->x) + (pVec->y*pVec->y) + (pVec->z*pVec->z)+ (pVec->w*pVec->w));

}




void VPCALL CSIMD_Generic::vector4D_Normalize(CVec4D* pVec)
{

	float inv_mag = 1.0f / (sqrtf((pVec->x*pVec->x) + (pVec->y*pVec->y) + (pVec->z*pVec->z)+ (pVec->w*pVec->w))+CMath::FLT_EPSILON);

	pVec->x *= inv_mag;
	pVec->y *= inv_mag;
	pVec->z *= inv_mag;
	pVec->w *= inv_mag;

}

void VPCALL CSIMD_Generic::vector4D_NormalizeOf(CVec4D* pOut, const CVec4D* pVec)
{
	float inv_mag = 1.0f / (sqrtf((pVec->x*pVec->x) + (pVec->y*pVec->y) + (pVec->z*pVec->z)+ (pVec->w*pVec->w))+CMath::FLT_EPSILON);

	pOut->x = pVec->x * inv_mag;
	pOut->y = pVec->y * inv_mag;
	pOut->z = pVec->z * inv_mag;
	pOut->w = pVec->w * inv_mag;


}

float  VPCALL CSIMD_Generic::vector4D_Distance(const CVec4D* pVec1, const CVec4D* pVec2)
{
	/* TODO: Optimize me completely */
	CVec4D diff;
	vector4D_DiffOf(&diff, pVec1, pVec2);
	return vector4D_Length(&diff);
}


void VPCALL CSIMD_Generic::vector4D_AlignedSum(CVec4D* pOut, const CVec4D* pIn)
{
	pOut->x += pIn->x;
	pOut->y += pIn->y;
	pOut->z += pIn->z;
	pOut->w += pIn->w;
}

void VPCALL CSIMD_Generic::vector4D_AlignedSumOf(CVec4D* pOut, const CVec4D* pIn1, const CVec4D* pIn2)
{
	pOut->x = pIn1->x + pIn2->x;
	pOut->y = pIn1->y + pIn2->y;
	pOut->z = pIn1->z + pIn2->z;
	pOut->w = pIn1->w + pIn2->w;
}


//Subtracts pIn from pOut and stores the result in pOut. POut = pIn-pOut
void VPCALL CSIMD_Generic::vector4D_AlignedDiff(CVec4D* pOut, CVec4D* pIn)
{
	pOut->x -= pIn->x;
	pOut->y -= pIn->y;
	pOut->z -= pIn->z;
	pOut->w -= pIn->w;
}


//Subtracts pIn2 from pIn1 and stores the result in pOut.  pout = pIn1-Pin2
void VPCALL CSIMD_Generic::vector4D_AlignedDiffOf(CVec4D* pOut, const CVec4D* pIn1, const CVec4D* pIn2)
{
	pOut->x = pIn1->x - pIn2->x;
	pOut->y = pIn1->y - pIn2->y;
	pOut->z = pIn1->z - pIn2->z;
	pOut->w = pIn1->w - pIn2->w;
}
//Scales the components of pOut by scalar and stores the result in pOut.
void VPCALL CSIMD_Generic::vector4D_AlignedScale(CVec4D* pOut, float scalar)
{
	pOut->x *= scalar;
	pOut->y *= scalar;
	pOut->z *= scalar;
	pOut->w *= scalar;

}
//Scales the components of pIn by scalar and stores the result in pOut.
void VPCALL CSIMD_Generic::vector4D_AlignedScaleOf(CVec4D* pOut, const CVec4D* pIn, float scalar)
{
	pOut->x = pIn->x * scalar;
	pOut->y = pIn->y * scalar;
	pOut->z = pIn->z * scalar;
	pOut->w = pIn->w * scalar;

}
//Calculates the dot product (known as the inner product) of pSrc1 and pSrc2, that represents the cosine of the angle between them, if both are normalized.

float  VPCALL CSIMD_Generic::vector4D_AlignedDot(const CVec4D* pSrc1, const CVec4D* pSrc2)
{
	return (pSrc1->x * pSrc2->x) + (pSrc1->y * pSrc2->y) + (pSrc1->z * pSrc2->z)+ (pSrc1->w * pSrc2->w);

}

float  VPCALL CSIMD_Generic::vector4D_AlignedLengthSq(const CVec4D* pVec)
{
		return (pVec->x*pVec->x) + (pVec->y*pVec->y) + (pVec->z*pVec->z)+ (pVec->w*pVec->w);
}

float  VPCALL CSIMD_Generic::vector4D_AlignedLength(const CVec4D* pVec)
{
	return sqrtf((pVec->x*pVec->x) + (pVec->y*pVec->y) + (pVec->z*pVec->z)+ (pVec->w*pVec->w));

}




void VPCALL CSIMD_Generic::vector4D_AlignedNormalize(CVec4D* pVec)
{

	float inv_mag = 1.0f / (sqrtf((pVec->x*pVec->x) + (pVec->y*pVec->y) + (pVec->z*pVec->z)+ (pVec->w*pVec->w))+CMath::FLT_EPSILON);

	pVec->x *= inv_mag;
	pVec->y *= inv_mag;
	pVec->z *= inv_mag;
	pVec->w *= inv_mag;

}

void VPCALL CSIMD_Generic::vector4D_AlignedNormalizeOf(CVec4D* pOut, const CVec4D* pVec)
{
	float inv_mag = 1.0f / (sqrtf((pVec->x*pVec->x) + (pVec->y*pVec->y) + (pVec->z*pVec->z)+ (pVec->w*pVec->w))+CMath::FLT_EPSILON);

	pOut->x = pVec->x * inv_mag;
	pOut->y = pVec->y * inv_mag;
	pOut->z = pVec->z * inv_mag;
	pOut->w = pVec->w * inv_mag;


}

float  VPCALL CSIMD_Generic::vector4D_AlignedDistance(const CVec4D* pVec1, const CVec4D* pVec2)
{
	/* TODO: Optimize me completely */
	CVec4D diff;
	vector4D_DiffOf(&diff, pVec1, pVec2);
	return vector4D_Length(&diff);
}

//===========trigonometry=====================================
	float VPCALL  CSIMD_Generic::invSqrt( float x ){
	return CMath::invSqrt16(x);
	};
	void  VPCALL  CSIMD_Generic::InvSqrt4( float x[4] ){
        for (int a=0;a<4;a++){
        x[a]=CMath::invSqrt16(x[a]);
       }
	};
	float VPCALL  CSIMD_Generic::sinZeroHalfPI( float x ){
	return CMath::sin16(x);
	};
	void  VPCALL  CSIMD_Generic::sin4ZeroHalfPI( float entrada[4], float saida[4] ){
        for (int a=0;a<4;a++){
        saida[a]=CMath::sin16(entrada[a]);
       }
	};
	float VPCALL  CSIMD_Generic::sin( float a ){
	return CMath::sin16(a);
	};
	void  VPCALL  CSIMD_Generic::sin4( float entrada[4], float saida[4] ){
       for (int a=0;a<4;a++){
        saida[a]=CMath::sin16(entrada[a]);
       }

	};
	float VPCALL  CSIMD_Generic::cosZeroHalfPI( float x ){
	return CMath::cos16(x);
	};
	void  VPCALL  CSIMD_Generic::cos4ZeroHalfPI( float entrada[4], float saida[4] ){
        for (int a=0;a<4;a++){
        saida[a]=CMath::cos16(entrada[a]);
       }
	};
	float VPCALL  CSIMD_Generic::cos( float a ){
	return CMath::cos16(a);
	};
	void  VPCALL  CSIMD_Generic::cos4( float entrada[4], float saida[4] ){
    for (int a=0;a<4;a++){
        saida[a]=CMath::cos16(entrada[a]);
       }
	};
	void  VPCALL  CSIMD_Generic::sincos( float a, float &sin, float &cos ){
        cos = CMath::cos16(a);
        sin = CMath::sin16(a);
	};
	void  VPCALL  CSIMD_Generic::sincos4( float entrada[4], float sin[4], float cos[4] ){
        for (int a=0;a<4;a++){
        cos[a]=CMath::cos16(entrada[a]);
        sin[a]=CMath::sin16(entrada[a]);
       }
	};
	float VPCALL  CSIMD_Generic::aTanPositive( float y, float x ){
	return CMath::atan16(y,x);
	};
	void  VPCALL  CSIMD_Generic::aTan4Positive( float y[4], float x[4], float resultado[4] ){
        for (int a=0;a<4;a++){
        resultado[a]=CMath::atan16(y[a],x[a]);
       }
	};
	float VPCALL  CSIMD_Generic::atan( float y, float x ){
	return CMath::atan16(y,x);
	};
	void  VPCALL  CSIMD_Generic::aTan4( float y[4], float x[4], float resultado[4] ){
        for (int a=0;a<4;a++){
        resultado[a]=CMath::atan16(y[a],x[a]);
       }

	};

//========CMat4D==============================

void VPCALL CSIMD_Generic::mat4D_Sum(CMat4D* Out, const CMat4D* In)
{


	Out->mat[0].x += In->mat[0].x; Out->mat[0].y += In->mat[0].y; Out->mat[0].z += In->mat[0].z; Out->mat[0].w += In->mat[0].w;
	Out->mat[1].x += In->mat[1].x; Out->mat[1].y += In->mat[1].y; Out->mat[1].z += In->mat[1].z; Out->mat[1].w += In->mat[1].w;
	Out->mat[2].x += In->mat[2].x; Out->mat[2].y += In->mat[2].y; Out->mat[2].z += In->mat[2].z; Out->mat[2].w += In->mat[2].w;
	Out->mat[3].x += In->mat[3].x; Out->mat[3].y += In->mat[3].y; Out->mat[3].z += In->mat[3].z; Out->mat[3].w += In->mat[3].w;

}


void VPCALL CSIMD_Generic::mat4D_SumOf(CMat4D* Out, const CMat4D* In1, const CMat4D* In2)
{

		Out->mat[0].x=In1->mat[0].x + In2->mat[0].x;
		Out->mat[0].y=In1->mat[0].y + In2->mat[0].y;
		Out->mat[0].z=In1->mat[0].z + In2->mat[0].z;
		Out->mat[0].w=In1->mat[0].w + In2->mat[0].w,
		Out->mat[1].x=In1->mat[1].x + In2->mat[1].x;
		Out->mat[1].y=In1->mat[1].y + In2->mat[1].y;
		Out->mat[1].z=In1->mat[1].z + In2->mat[1].z;
		Out->mat[1].w=In1->mat[1].w + In2->mat[1].w,
		Out->mat[2].x=In1->mat[2].x + In2->mat[2].x;
		Out->mat[2].y=In1->mat[2].y + In2->mat[2].y;
		Out->mat[2].z=In1->mat[2].z + In2->mat[2].z;
		Out->mat[2].w=In1->mat[2].w + In2->mat[2].w,
		Out->mat[3].x=In1->mat[3].x + In2->mat[3].x;
		Out->mat[3].y=In1->mat[3].y + In2->mat[3].y;
		Out->mat[3].z=In1->mat[3].z + In2->mat[3].z;
		Out->mat[3].w=In1->mat[3].w + In2->mat[3].w ;

}

void VPCALL CSIMD_Generic::mat4D_Diff(CMat4D* Out, const CMat4D* In)
{


	Out->mat[0].x -= In->mat[0].x; Out->mat[0].y -= In->mat[0].y; Out->mat[0].z -= In->mat[0].z; Out->mat[0].w -= In->mat[0].w;
	Out->mat[1].x -= In->mat[1].x; Out->mat[1].y -= In->mat[1].y; Out->mat[1].z -= In->mat[1].z; Out->mat[1].w -= In->mat[1].w;
	Out->mat[2].x -= In->mat[2].x; Out->mat[2].y -= In->mat[2].y; Out->mat[2].z -= In->mat[2].z; Out->mat[2].w -= In->mat[2].w;
	Out->mat[3].x -= In->mat[3].x; Out->mat[3].y -= In->mat[3].y; Out->mat[3].z -= In->mat[3].z; Out->mat[3].w -= In->mat[3].w;

}

void VPCALL CSIMD_Generic::mat4D_DiffOf(CMat4D* Out, const CMat4D* In1, const CMat4D* In2)
{

		Out->mat[0].x=In1->mat[0].x - In2->mat[0].x;
		Out->mat[0].y=In1->mat[0].y - In2->mat[0].y;
		Out->mat[0].z=In1->mat[0].z - In2->mat[0].z;
		Out->mat[0].w=In1->mat[0].w - In2->mat[0].w,
		Out->mat[1].x=In1->mat[1].x - In2->mat[1].x;
		Out->mat[1].y=In1->mat[1].y - In2->mat[1].y;
		Out->mat[1].z=In1->mat[1].z - In2->mat[1].z;
		Out->mat[1].w=In1->mat[1].w - In2->mat[1].w,
		Out->mat[2].x=In1->mat[2].x - In2->mat[2].x;
		Out->mat[2].y=In1->mat[2].y - In2->mat[2].y;
		Out->mat[2].z=In1->mat[2].z - In2->mat[2].z;
		Out->mat[2].w=In1->mat[2].w - In2->mat[2].w,
		Out->mat[3].x=In1->mat[3].x - In2->mat[3].x;
		Out->mat[3].y=In1->mat[3].y - In2->mat[3].y;
		Out->mat[3].z=In1->mat[3].z - In2->mat[3].z;
		Out->mat[3].w=In1->mat[3].w - In2->mat[3].w ;


}

void VPCALL CSIMD_Generic::mat4D_Scale(CMat4D* mtx, float scalar)
{


	/* Don't waste CPU time on this */
	if(scalar == 1.0f)
		return;
	mtx->mat[0].x *= scalar; mtx->mat[0].y *= scalar; mtx->mat[0].z *= scalar; mtx->mat[0].w *= scalar;
	mtx->mat[1].x *= scalar; mtx->mat[1].y *= scalar; mtx->mat[1].z *= scalar; mtx->mat[1].w *= scalar;
	mtx->mat[2].x *= scalar; mtx->mat[2].y *= scalar; mtx->mat[2].z *= scalar; mtx->mat[2].w *= scalar;
	mtx->mat[3].x *= scalar; mtx->mat[3].y *= scalar; mtx->mat[3].z *= scalar; mtx->mat[3].w *= scalar;


}

void VPCALL CSIMD_Generic::mat4D_ScaleOf(CMat4D* pOut, const CMat4D* pIn, float scalar)
{

	/* Don't waste CPU time */
	if(scalar == 1.0f)
	{
		memcpy((void*)pOut, (const void*)pIn, sizeof(CMat4D));
		return;
	}

	pOut->mat[0].x=  pIn->mat[0].x * scalar;
	pOut->mat[0].y=  pIn->mat[0].y * scalar;
	pOut->mat[0].z=  pIn->mat[0].z * scalar;
	pOut->mat[0].w=  pIn->mat[0].w * scalar;
	pOut->mat[1].x=  pIn->mat[1].x * scalar;
	pOut->mat[1].y=  pIn->mat[1].y * scalar;
	pOut->mat[1].z=  pIn->mat[1].z * scalar;
	pOut->mat[1].w=  pIn->mat[1].w * scalar;
	pOut->mat[2].x=  pIn->mat[2].x * scalar;
	pOut->mat[2].y=  pIn->mat[2].y * scalar;
	pOut->mat[2].z=  pIn->mat[2].z * scalar;
	pOut->mat[2].w=  pIn->mat[2].w * scalar;
	pOut->mat[3].x=  pIn->mat[3].x * scalar;
	pOut->mat[3].y=  pIn->mat[3].y * scalar;
	pOut->mat[3].z=  pIn->mat[3].z * scalar;
	pOut->mat[3].w=  pIn->mat[3].w * scalar;



}

void VPCALL CSIMD_Generic::mat4D_Multiply(CMat4D* pLeft, const CMat4D* pRight)
{


	int i, j;
	const float *m1Ptr, *m2Ptr;
	float *dstPtr;
	CMat4D dst;

	m1Ptr = pLeft->toFloatPtr();
	m2Ptr = pRight->toFloatPtr();
	dstPtr = reinterpret_cast<float *>(&dst);

	for ( i = 0; i < 4; i++ ) {
		for ( j = 0; j < 4; j++ ) {
			*dstPtr = m1Ptr[0] * m2Ptr[ 0 * 4 + j ]
					+ m1Ptr[1] * m2Ptr[ 1 * 4 + j ]
					+ m1Ptr[2] * m2Ptr[ 2 * 4 + j ]
					+ m1Ptr[3] * m2Ptr[ 3 * 4 + j ];
			dstPtr++;
		}
		m1Ptr += 4;
	}
	*pLeft=dst;

}

void VPCALL CSIMD_Generic::mat4D_MultiplyOf(CMat4D* pOut, const CMat4D* pLeft, const CMat4D* pRight)
{

	int i, j;
	const float *m1Ptr, *m2Ptr;
	float *dstPtr;

	m1Ptr = pLeft->toFloatPtr();
	m2Ptr = pRight->toFloatPtr();
	dstPtr = pOut->toFloatPtr();

	for ( i = 0; i < 4; i++ ) {
		for ( j = 0; j < 4; j++ ) {
			*dstPtr = m1Ptr[0] * m2Ptr[ 0 * 4 + j ]
					+ m1Ptr[1] * m2Ptr[ 1 * 4 + j ]
					+ m1Ptr[2] * m2Ptr[ 2 * 4 + j ]
					+ m1Ptr[3] * m2Ptr[ 3 * 4 + j ];
			dstPtr++;
		}
		m1Ptr += 4;
	}

}



void VPCALL CSIMD_Generic::mat4D_Transpose(CMat4D* pIn)
{

	CMat4D	transpose;
	int		i, j;

	for( i = 0; i < 4; i++ ) {
		for( j = 0; j < 4; j++ ) {
			transpose[ i ][ j ] = pIn->mat[ j ][ i ];
        }
	}
	*pIn= transpose;

}



void VPCALL CSIMD_Generic::mat4D_TransposeOf(CMat4D* pOut, const CMat4D* pIn)
{
	*pOut=pIn->transpose();
	CMat4D	transpose;
	int		i, j;

	for( i = 0; i < 4; i++ ) {
		for( j = 0; j < 4; j++ ) {
			transpose[ i ][ j ] = pIn->mat[ j ][ i ];
        }
	}
	*pOut= transpose;


}

void VPCALL CSIMD_Generic::mat4D_VectorMultiply(CVec3D* pVec, const CMat4D* pMat)
{
	/*
		Does a normal vector x matrix, but since w = 1.0f, the final column of the matrix is
		merely added.
	*/


	float In[3] = { pVec->x, pVec->y, pVec->z };

	pVec->x = (In[0]*pMat->mat[0].x) + (In[1]*pMat->mat[1].x) + (In[2]*pMat->mat[2].x) + pMat->mat[3].x;
	pVec->y = (In[0]*pMat->mat[0].y) + (In[1]*pMat->mat[1].y) + (In[2]*pMat->mat[2].y) + pMat->mat[3].y;
	pVec->z = (In[0]*pMat->mat[0].z) + (In[1]*pMat->mat[1].z) + (In[2]*pMat->mat[2].z) + pMat->mat[3].z;


}

void VPCALL CSIMD_Generic::mat4D_VectorMultiplyOf(CVec3D* pOut, const CVec3D* pIn, const CMat4D* pMat)
{
	/*
		Does a normal vector x matrix, but since w = 1.0f, the final column of the matrix is
		merely added.
	*/


	pOut->x = (pIn->x*pMat->mat[0].x) + (pIn->y*pMat->mat[1].x) + (pIn->z*pMat->mat[2].x) + pMat->mat[3].x;
	pOut->y = (pIn->x*pMat->mat[0].y) + (pIn->y*pMat->mat[1].y) + (pIn->z*pMat->mat[2].y) + pMat->mat[3].y;
	pOut->z = (pIn->x*pMat->mat[0].z) + (pIn->y*pMat->mat[1].z) + (pIn->z*pMat->mat[2].z) + pMat->mat[3].z;


}

void VPCALL CSIMD_Generic::mat4D_VectorMultiply(CVec4D* pOut4D, const CMat4D* pMat)
{

	float In[4] = { pOut4D->x, pOut4D->y, pOut4D->z, pOut4D->w };

	pOut4D->x = (In[0]*pMat->mat[0].x) + (In[1]*pMat->mat[1].x) + (In[2]*pMat->mat[2].x) + (In[3]*pMat->mat[3].x);
	pOut4D->y = (In[0]*pMat->mat[0].y) + (In[1]*pMat->mat[1].y) + (In[2]*pMat->mat[2].y) + (In[3]*pMat->mat[3].y);
	pOut4D->z = (In[0]*pMat->mat[0].z) + (In[1]*pMat->mat[1].z) + (In[2]*pMat->mat[2].z) + (In[3]*pMat->mat[3].z);
	pOut4D->w = (In[0]*pMat->mat[0].w) + (In[1]*pMat->mat[1].w) + (In[2]*pMat->mat[2].w) + (In[3]*pMat->mat[3].w);



}


void VPCALL CSIMD_Generic::mat4D_VectorMultiplyOf(CVec4D* pOut4D, const CVec4D* pIn4D, const CMat4D* pMat)
{


	pOut4D->x = (pIn4D->x*pMat->mat[0].x) + (pIn4D->y*pMat->mat[1].x) + (pIn4D->z*pMat->mat[2].x) + (pIn4D->w*pMat->mat[3].x);
	pOut4D->y = (pIn4D->x*pMat->mat[0].y) + (pIn4D->y*pMat->mat[1].y) + (pIn4D->z*pMat->mat[2].y) + (pIn4D->w*pMat->mat[3].y);
	pOut4D->z = (pIn4D->x*pMat->mat[0].z) + (pIn4D->y*pMat->mat[1].z) + (pIn4D->z*pMat->mat[2].z) + (pIn4D->w*pMat->mat[3].z);
	pOut4D->w = (pIn4D->x*pMat->mat[0].w) + (pIn4D->y*pMat->mat[1].w) + (pIn4D->z*pMat->mat[2].w) + (pIn4D->w*pMat->mat[3].w);

}


void VPCALL CSIMD_Generic::mat4D_ToRotate(CMat4D* pMat, float yaw, float pitch, float roll)
{

	float cos_yaw = cos(yaw);
	float sin_yaw = sin(yaw);
	float cos_pitch = cos(pitch);
	float sin_pitch = sin(pitch);
	float cos_roll = cos(roll);
	float sin_roll = sin(roll);

	pMat->mat[0].x = (cos_yaw * cos_pitch);
	pMat->mat[0].y = (cos_yaw * sin_pitch * sin_roll) - (sin_yaw * cos_roll);
	pMat->mat[0].z = (cos_yaw * sin_pitch * cos_roll) + (sin_yaw * sin_roll);
	pMat->mat[0].w = 0.0;

	pMat->mat[1].x = (sin_yaw * cos_pitch);
	pMat->mat[1].y = (cos_yaw * cos_roll) + (sin_yaw * sin_pitch * sin_roll);
	pMat->mat[1].z = (sin_yaw * sin_pitch * cos_roll) - (cos_yaw * sin_roll);
	pMat->mat[1].w = 0.0;

	pMat->mat[2].x = -sin_pitch;
	pMat->mat[2].y = (cos_pitch * sin_roll);
	pMat->mat[2].z = (cos_pitch * cos_roll);
	pMat->mat[2].w = 0.0;

	pMat->mat[3].x = 0.0;
	pMat->mat[3].y = 0.0;
	pMat->mat[3].z = 0.0;
	pMat->mat[3].w = 1.0;
}

void VPCALL CSIMD_Generic::mat4D_ToRotateOf(CMat4D* pMat, const CVec3D* pYawPitchRoll)
{
	float cos_yaw = cos(pYawPitchRoll->x);
	float sin_yaw = sin(pYawPitchRoll->x);
	float cos_pitch = cos(pYawPitchRoll->y);
	float sin_pitch = sin(pYawPitchRoll->y);
	float cos_roll = cos(pYawPitchRoll->z);
	float sin_roll = sin(pYawPitchRoll->z);

	pMat->mat[0].x = (cos_yaw * cos_pitch);
	pMat->mat[0].y = (cos_yaw * sin_pitch * sin_roll) - (sin_yaw * cos_roll);
	pMat->mat[0].z = (cos_yaw * sin_pitch * cos_roll) + (sin_yaw * sin_roll);
	pMat->mat[0].w = 0.0;

	pMat->mat[1].x = (sin_yaw * cos_pitch);
	pMat->mat[1].y = (cos_yaw * cos_roll) + (sin_yaw * sin_pitch * sin_roll);
	pMat->mat[1].z = (sin_yaw * sin_pitch * cos_roll) - (cos_yaw * sin_roll);
	pMat->mat[1].w = 0.0;

	pMat->mat[2].x = -sin_pitch;
	pMat->mat[2].y = (cos_pitch * sin_roll);
	pMat->mat[2].z = (cos_pitch * cos_roll);
	pMat->mat[2].w = 0.0;

	pMat->mat[3].x = 0.0;
	pMat->mat[3].y = 0.0;
	pMat->mat[3].z = 0.0;
	pMat->mat[3].w = 1.0;
}

void  VPCALL CSIMD_Generic::mat4D_AlignedSum(CMat4D* pMat, const CMat4D* pIn){
	mat4D_Sum(pMat, pIn);
}
void  VPCALL CSIMD_Generic::mat4D_AlignedSumOf(CMat4D* pMat, const CMat4D* pIn1, const CMat4D* pIn2){
	 mat4D_SumOf(pMat, pIn1, pIn2);
}
void  VPCALL CSIMD_Generic::mat4D_AlignedDiff(CMat4D* pMat, const CMat4D* pIn){
	mat4D_Diff( pMat, pIn);
}
void  VPCALL CSIMD_Generic::mat4D_AlignedDiffOf(CMat4D* pMat, const CMat4D* pLeft, const CMat4D* pRight){
	mat4D_DiffOf(pMat, pLeft, pRight);
}
void  VPCALL CSIMD_Generic::mat4D_AlignedScale(CMat4D* pMat, float scalar){
	mat4D_Scale(pMat, scalar);
}
void  VPCALL CSIMD_Generic::mat4D_AlignedScaleOf(CMat4D* pOut, const CMat4D* pIn, float scalar){
	mat4D_ScaleOf(pOut, pIn, scalar);
}
void  VPCALL CSIMD_Generic::mat4D_AlignedMultiply(CMat4D* pLeft, const CMat4D* pRight){
	mat4D_Multiply(pLeft, pRight);
}
void  VPCALL CSIMD_Generic::mat4D_AlignedMultiplyOf(CMat4D* pOut, const CMat4D* pLeft, const CMat4D* pRight){
	mat4D_MultiplyOf(pOut, pLeft, pRight);
}
void  VPCALL CSIMD_Generic::mat4D_AlignedTranspose(CMat4D* pIn){
	mat4D_Transpose(pIn);
}
void  VPCALL CSIMD_Generic::mat4D_AlignedTransposeOf(CMat4D* pOut, const CMat4D* pIn){
	mat4D_TransposeOf(pOut, pIn);
}
void  VPCALL CSIMD_Generic::mat4D_AlignedVectorMultiply(CVec3D* pOut, const CMat4D* pIn){
	mat4D_VectorMultiply(pOut, pIn);
}
void  VPCALL CSIMD_Generic::mat4D_AlignedVectorMultiplyOf(CVec3D* pOut, const CVec3D* pIn, const CMat4D* pMat){
	mat4D_VectorMultiplyOf(pOut, pIn, pMat);
}
void  VPCALL CSIMD_Generic::mat4D_AlignedVectorMultiply(CVec4D* pOut4D, const CMat4D* pIn){
	mat4D_VectorMultiply(pOut4D, pIn);
}
void  VPCALL CSIMD_Generic::mat4D_AlignedVectorMultiplyOf(CVec4D* pOut4D, const CVec4D* pIn4D, const CMat4D* pMat){
	mat4D_VectorMultiplyOf(pOut4D, pIn4D, pMat);
}
void  VPCALL CSIMD_Generic::mat4D_AlignedToRotate(CMat4D* pMat, float yaw, float pitch, float roll){
	mat4D_ToRotate(pMat, yaw, pitch, roll);
}
void  VPCALL CSIMD_Generic::mat4D_AlignedToRotateOf(CMat4D* pMat, const CVec3D* pYawPitchRoll){
	mat4D_ToRotateOf(pMat, pYawPitchRoll);
}

void  VPCALL CSIMD_Generic::quat_to_mat4x4(sf_m128 q, sf_m128 t, sf_m128 *mat){

	CVec3D vec= *reinterpret_cast<CVec3D *>(&t);
	CQuaternion orientation= *reinterpret_cast<CQuaternion *>(&q);
	CMat4D m1 = *reinterpret_cast<CMat4D *>(mat);
	m1.setRotatePart(orientation);
	m1.setRow(3, 0, 0, 0, 1);
	m1.setTranslatePart(vec);
}

/// Compute the product M*v, where M is a 4x4 matrix denoted by an array of 4 sf_m128's, and v is a 4x1 vector.
// If we have SSE 4.1, we can use the dpps (dot product) instruction, _mm_dp_ps intrinsic.
/// Compute the product M*v, where M is a 4x4 matrix denoted by an array of 4 sf_m128's, and v is a 4x1 vector.
// If we have SSE3, we can repeatedly use haddps to accumulate the result.
sf_m128  VPCALL CSIMD_Generic::mat4x4_mul_sse(const sf_m128 *matrix, sf_m128 vector){
	CMat4D  mat = *reinterpret_cast<CMat4D *>(&matrix);
	return CVec4D(mat.Row(0)* vector,
	              mat.Row(1)* vector,
	              mat.Row(2)* vector,
	              mat.Row(3)* vector).v;

}
 void   VPCALL CSIMD_Generic::mat4x4_mul_sse(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2){
	int i, j;
	const float *m1Ptr, *m2Ptr;
	float *dstPtr;


	m1Ptr = reinterpret_cast<const float *>(m1);
	m2Ptr = reinterpret_cast<const float *>(m2);
	dstPtr = reinterpret_cast<float *>(out);

	for ( i = 0; i < 4; i++ ) {
		for ( j = 0; j < 4; j++ ) {
			*dstPtr = m1Ptr[0] * m2Ptr[ 0 * 4 + j ]
					+ m1Ptr[1] * m2Ptr[ 1 * 4 + j ]
					+ m1Ptr[2] * m2Ptr[ 2 * 4 + j ]
					+ m1Ptr[3] * m2Ptr[ 3 * 4 + j ];
			dstPtr++;
		}
		m1Ptr += 4;
	}


 }

sf_m128  VPCALL CSIMD_Generic::colmajor_mat4x4_mul_sse(const sf_m128 *matrix, sf_m128 vector){
//todo:
	sf_m128 teste;
	return teste;
}


/**
\brief Compute the product M*v, where M is a 3x4 matrix denoted by an array of 4 sf_m128's, and v is a 4x1 vector.
\return sf_m128 representing Vector 4D
*/
 sf_m128  VPCALL CSIMD_Generic::mat3x4_mul_sse(const sf_m128 *matrix, sf_m128 vec){
	CMatJoint3x4 m1 = *reinterpret_cast<const CMatJoint3x4 *>(matrix);
	CVec4D vector = *reinterpret_cast<CVec4D *>(&vec);
	CVec4D lin0(m1.mat[0*4+0],m1.mat[0*4+1],m1.mat[0*4+2],m1.mat[0*4+3]);
	CVec4D lin1(m1.mat[1*4+0],m1.mat[1*4+1],m1.mat[1*4+2],m1.mat[1*4+3]);
	CVec4D lin2(m1.mat[2*4+0],m1.mat[2*4+1],m1.mat[2*4+2],m1.mat[2*4+3]);


	CVec4D ret(lin0 * vector,
				  lin1 * vector,
				  lin2 * vector,
				  vector.w);
	return *reinterpret_cast<sf_m128 *>(&ret);
 }

 void  VPCALL CSIMD_Generic::mat3x4_mul_sse(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2){
	CMatJoint3x4 r= *reinterpret_cast<CMatJoint3x4 *>(out);
	CMatJoint3x4 mulmat= *reinterpret_cast<const CMatJoint3x4 *>(m1);
	CMatJoint3x4 rhs= *reinterpret_cast<const CMatJoint3x4 *>(m2);


	CVec3D lin0_3(mulmat.mat[0 * 4 + 0],mulmat.mat[0 * 4 + 1],mulmat.mat[0 * 4 + 2]);
	CVec3D lin1_3(mulmat.mat[1 * 4 + 0],mulmat.mat[1 * 4 + 1],mulmat.mat[1 * 4 + 2]);
	CVec3D lin2_3(mulmat.mat[2 * 4 + 0],mulmat.mat[2 * 4 + 1],mulmat.mat[2 * 4 + 2]);

	CVec3D rhsCol0_3=rhs.Col(0);
	CVec3D rhsCol1_3=rhs.Col(1);
	CVec3D rhsCol2_3=rhs.Col(2);
	CVec3D rhsCol3_3=rhs.Col(3);

	r.mat[0 * 4 + 0] = lin0_3*rhsCol0_3;
	r.mat[0 * 4 + 1] = lin0_3*rhsCol1_3;
	r.mat[0 * 4 + 2] = lin0_3*rhsCol2_3;
	r.mat[0 * 4 + 3] = (lin0_3*rhsCol3_3) + mulmat.mat[0 * 4 + 3];

	r.mat[1 * 4 + 0] = lin1_3*rhsCol0_3;
	r.mat[1 * 4 + 1] = lin1_3*rhsCol1_3;
	r.mat[1 * 4 + 2] = lin1_3*rhsCol2_3;
	r.mat[1 * 4 + 3] = (lin1_3*rhsCol3_3) + mulmat.mat[1 * 4 + 3];

	r.mat[1 * 4 + 0] = lin2_3*rhsCol0_3;
	r.mat[1 * 4 + 1] = lin2_3*rhsCol1_3;
	r.mat[1 * 4 + 2] = lin2_3*rhsCol2_3;
	r.mat[1 * 4 + 3] = (lin2_3*rhsCol3_3) + mulmat.mat[2 * 4 + 3];
}

CVec3D  VPCALL CSIMD_Generic::mat3x4_mul_vec(const sf_m128 *matrix, sf_m128 vector){
	CMat4D  mat = *reinterpret_cast<CMat4D *>(&matrix);
	CVec3D  vec = *reinterpret_cast<CVec3D *>(&vector);
return CVec3D(mat.Row(0).toVec3()*vec,
					  mat.Row(0).toVec3()* vec,
					  mat.Row(0).toVec3()* vec);
}

/**
\brief multiplica duas matrizes CMat4D  m1*m2
\param [out] out resultado da multiplicação
\param m1 matriz a ser mutiplicada
\param m2 matria a ser multiplicada
**/
void  VPCALL CSIMD_Generic::mat4x4_mul_dpps(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2){
	mat4x4_mul_sse(out,m1,m2);
}
void  VPCALL CSIMD_Generic::mat4x4_mul_dpps_2(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2){
	mat4x4_mul_sse(out,m1,m2);

}
void  VPCALL CSIMD_Generic::mat4x4_mul_dpps_3(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2){
	mat4x4_mul_sse(out,m1,m2);

}



float  VPCALL CSIMD_Generic::mat4x4_inverse(const CMat4D *matrix, CMat4D *out){
		// 84+4+16 = 104 multiplications
	//			   1 division
	double det, invDet;

	// 2x2 sub-determinants required to calculate 4x4 determinant
	float det2_01_01 = matrix->mat[0][0] * matrix->mat[1][1] - matrix->mat[0][1] * matrix->mat[1][0];
	float det2_01_02 = matrix->mat[0][0] * matrix->mat[1][2] - matrix->mat[0][2] * matrix->mat[1][0];
	float det2_01_03 = matrix->mat[0][0] * matrix->mat[1][3] - matrix->mat[0][3] * matrix->mat[1][0];
	float det2_01_12 = matrix->mat[0][1] * matrix->mat[1][2] - matrix->mat[0][2] * matrix->mat[1][1];
	float det2_01_13 = matrix->mat[0][1] * matrix->mat[1][3] - matrix->mat[0][3] * matrix->mat[1][1];
	float det2_01_23 = matrix->mat[0][2] * matrix->mat[1][3] - matrix->mat[0][3] * matrix->mat[1][2];

	// 3x3 sub-determinants required to calculate 4x4 determinant
	float det3_201_012 = matrix->mat[2][0] * det2_01_12 - matrix->mat[2][1] * det2_01_02 + matrix->mat[2][2] * det2_01_01;
	float det3_201_013 = matrix->mat[2][0] * det2_01_13 - matrix->mat[2][1] * det2_01_03 + matrix->mat[2][3] * det2_01_01;
	float det3_201_023 = matrix->mat[2][0] * det2_01_23 - matrix->mat[2][2] * det2_01_03 + matrix->mat[2][3] * det2_01_02;
	float det3_201_123 = matrix->mat[2][1] * det2_01_23 - matrix->mat[2][2] * det2_01_13 + matrix->mat[2][3] * det2_01_12;

	det = ( - det3_201_123 * matrix->mat[3][0] + det3_201_023 * matrix->mat[3][1] - det3_201_013 * matrix->mat[3][2] + det3_201_012 * matrix->mat[3][3] );

	if ( CMath::fabs( det ) < MATRIX_INVERSE_EPSILON ) {
		return false;
	}

	invDet = 1.0f / det;

	// remaining 2x2 sub-determinants
	float det2_03_01 = matrix->mat[0][0] * matrix->mat[3][1] - matrix->mat[0][1] * matrix->mat[3][0];
	float det2_03_02 = matrix->mat[0][0] * matrix->mat[3][2] - matrix->mat[0][2] * matrix->mat[3][0];
	float det2_03_03 = matrix->mat[0][0] * matrix->mat[3][3] - matrix->mat[0][3] * matrix->mat[3][0];
	float det2_03_12 = matrix->mat[0][1] * matrix->mat[3][2] - matrix->mat[0][2] * matrix->mat[3][1];
	float det2_03_13 = matrix->mat[0][1] * matrix->mat[3][3] - matrix->mat[0][3] * matrix->mat[3][1];
	float det2_03_23 = matrix->mat[0][2] * matrix->mat[3][3] - matrix->mat[0][3] * matrix->mat[3][2];

	float det2_13_01 = matrix->mat[1][0] * matrix->mat[3][1] - matrix->mat[1][1] * matrix->mat[3][0];
	float det2_13_02 = matrix->mat[1][0] * matrix->mat[3][2] - matrix->mat[1][2] * matrix->mat[3][0];
	float det2_13_03 = matrix->mat[1][0] * matrix->mat[3][3] - matrix->mat[1][3] * matrix->mat[3][0];
	float det2_13_12 = matrix->mat[1][1] * matrix->mat[3][2] - matrix->mat[1][2] * matrix->mat[3][1];
	float det2_13_13 = matrix->mat[1][1] * matrix->mat[3][3] - matrix->mat[1][3] * matrix->mat[3][1];
	float det2_13_23 = matrix->mat[1][2] * matrix->mat[3][3] - matrix->mat[1][3] * matrix->mat[3][2];

	// remaining 3x3 sub-determinants
	float det3_203_012 = matrix->mat[2][0] * det2_03_12 - matrix->mat[2][1] * det2_03_02 + matrix->mat[2][2] * det2_03_01;
	float det3_203_013 = matrix->mat[2][0] * det2_03_13 - matrix->mat[2][1] * det2_03_03 + matrix->mat[2][3] * det2_03_01;
	float det3_203_023 = matrix->mat[2][0] * det2_03_23 - matrix->mat[2][2] * det2_03_03 + matrix->mat[2][3] * det2_03_02;
	float det3_203_123 = matrix->mat[2][1] * det2_03_23 - matrix->mat[2][2] * det2_03_13 + matrix->mat[2][3] * det2_03_12;

	float det3_213_012 = matrix->mat[2][0] * det2_13_12 - matrix->mat[2][1] * det2_13_02 + matrix->mat[2][2] * det2_13_01;
	float det3_213_013 = matrix->mat[2][0] * det2_13_13 - matrix->mat[2][1] * det2_13_03 + matrix->mat[2][3] * det2_13_01;
	float det3_213_023 = matrix->mat[2][0] * det2_13_23 - matrix->mat[2][2] * det2_13_03 + matrix->mat[2][3] * det2_13_02;
	float det3_213_123 = matrix->mat[2][1] * det2_13_23 - matrix->mat[2][2] * det2_13_13 + matrix->mat[2][3] * det2_13_12;

	float det3_301_012 = matrix->mat[3][0] * det2_01_12 - matrix->mat[3][1] * det2_01_02 + matrix->mat[3][2] * det2_01_01;
	float det3_301_013 = matrix->mat[3][0] * det2_01_13 - matrix->mat[3][1] * det2_01_03 + matrix->mat[3][3] * det2_01_01;
	float det3_301_023 = matrix->mat[3][0] * det2_01_23 - matrix->mat[3][2] * det2_01_03 + matrix->mat[3][3] * det2_01_02;
	float det3_301_123 = matrix->mat[3][1] * det2_01_23 - matrix->mat[3][2] * det2_01_13 + matrix->mat[3][3] * det2_01_12;

	out->mat[0][0] =	- det3_213_123 * invDet;
	out->mat[1][0] = + det3_213_023 * invDet;
	out->mat[2][0] = - det3_213_013 * invDet;
	out->mat[3][0] = + det3_213_012 * invDet;

	out->mat[0][1] = + det3_203_123 * invDet;
	out->mat[1][1] = - det3_203_023 * invDet;
	out->mat[2][1] = + det3_203_013 * invDet;
	out->mat[3][1] = - det3_203_012 * invDet;

	out->mat[0][2] = + det3_301_123 * invDet;
	out->mat[1][2] = - det3_301_023 * invDet;
	out->mat[2][2] = + det3_301_013 * invDet;
	out->mat[3][2] = - det3_301_012 * invDet;

	out->mat[0][3] = - det3_201_123 * invDet;
	out->mat[1][3] = + det3_201_023 * invDet;
	out->mat[2][3] = - det3_201_013 * invDet;
	out->mat[3][3] = + det3_201_012 * invDet;

	return det;

}

/// Inverts a 3x4 affine transformation matrix (in row-major format) that only consists of rotation (+possibly mirroring) and translation.
void  VPCALL CSIMD_Generic::mat3x4_inverse_orthonormal(sf_m128 *matrix, sf_m128 *out){

		/* In this function, we seek to optimize the matrix inverse in the case this
		   matrix is orthonormal, i.e. it can be written in the following form:

					  [ R | T ]
				  M = [---+---]
					  [ 0 | 1 ]

		   where R is a 3x3 orthonormal (orthogonal vectors, normalized columns) rotation
		   matrix, and T is a 3x1 vector representing the translation performed by
		   this matrix.

		   In this form, the inverse of this matrix is simple to compute and will not
		   require the calculation of determinants or expensive Gaussian elimination. The
		   inverse is of form

						 [ R^t | R^t(-T) ]
				  M^-1 = [-----+---------]
						 [  0  |    1    ]

		   which can be seen by multiplying out M * M^(-1) in block form. Especially the top-
		   right cell turns out to (remember that R^(-1) == R^t since R is orthonormal)

				R * R^t(-T) + T * 1 == (R * R^t)(-T) + T == -T + T == 0, as expected.

		   Therefore the inversion requires only two steps: */

		// a) transpose the top-left 3x3 part in-place to produce R^t.
		CMatJoint3x4  mat = *reinterpret_cast<CMatJoint3x4 *>(&matrix);
		Swap(mat.Row(0).y , mat.Row(1).x);
		Swap(mat.Row(0).z , mat.Row(2).x );
		Swap(mat.Row(1).z , mat.Row(2).y );

		// b) replace the top-right 3x1 part by computing R^t(-T).
		mat.setTranslatePart(mat.transformDir(-mat.Row(0).w , -mat.Row(1).w, -mat.Row(2).w));
		out= reinterpret_cast<sf_m128 *>(&mat);

}

sf_m128  VPCALL CSIMD_Generic:: newtonRhapsonRecipStep(sf_m128 recip, sf_m128 estimate){
	//todo
	sf_m128 ret;
	return ret;
}

sf_m128  VPCALL CSIMD_Generic:: newtonRhapsonRecip(sf_m128 recip){
	//todo
	sf_m128 ret;
	return ret;
}

/// Computes the determinant of a 4x4 matrix.
float  VPCALL CSIMD_Generic::mat4x4_determinant(const CMat4D *matrix){
	float a = matrix->Row(0).x;
	float b = matrix->Row(0).y;
	float c = matrix->Row(0).z;
	float d = matrix->Row(0).w;

	float e = matrix->Row(1).x;
	float f = matrix->Row(1).y;
	float g = matrix->Row(1).z;
	float h = matrix->Row(1).w;

	// 2x2 sub-determinants
	float det2_01_01 = a * f - b * e;
	float det2_01_02 = a * g - c * e;
	float det2_01_03 = a * h - d * e;
	float det2_01_12 = b * g - c * f;
	float det2_01_13 = b * h - d * f;
	float det2_01_23 = c * h - d * g;

	float i = matrix->Row(2).x;
	float j = matrix->Row(2).y;
	float k = matrix->Row(2).z;
	float l = matrix->Row(2).w;


	// 3x3 sub-determinants
	float det3_201_012 = i * det2_01_12 - j * det2_01_02 + k * det2_01_01;
	float det3_201_013 = i * det2_01_13 - j * det2_01_03 + l * det2_01_01;
	float det3_201_023 = i * det2_01_23 - k * det2_01_03 + l * det2_01_02;
	float det3_201_123 = j * det2_01_23 - k * det2_01_13 + l * det2_01_12;

	float m = matrix->Row(3).x;
	float n = matrix->Row(3).y;
	float o = matrix->Row(3).z;
	float p = matrix->Row(3).w;

	return ( - det3_201_123 * m + det3_201_023 * n - det3_201_013 * o + det3_201_012 * p );


}

/// Computes the determinant of a 3x4 matrix stored in row-major format. (Treated as a square matrix with last row [0,0,0,1])
float  VPCALL CSIMD_Generic::mat3x4_determinant(const sf_m128 *row){

	   CMatJoint3x4 matj;
     	const float a = matj.Row(0).x;
		const float b = matj.Row(0).y;
		const float c = matj.Row(0).z;
		const float d = matj.Row(1).x;
		const float e = matj.Row(1).y;
		const float f = matj.Row(1).z;
		const float g = matj.Row(2).x;
		const float h = matj.Row(2).y;
		const float i = matj.Row(2).z;

		return a*e*i + b*f*g + c*d*h - a*f*h - b*d*i - c*e*g;

}

void  VPCALL CSIMD_Generic::mat3x4_transpose(const sf_m128 *src, sf_m128 *dst){
	//\todo

}




//=============Quaternion====================


void  VPCALL CSIMD_Generic:: quaternion_Normalize(CQuaternion* pQuat)
{


	float inv_len = 1.0f / (sqrtf((pQuat->x*pQuat->x) + (pQuat->y*pQuat->y) + (pQuat->z*pQuat->z) + (pQuat->w*pQuat->w)) + MATH::CMath::FLT_EPSILON);

	pQuat->x *= inv_len;
	pQuat->y *= inv_len;
	pQuat->z *= inv_len;
	pQuat->w *= inv_len;

}

void  VPCALL CSIMD_Generic:: quaternion_NormalizeOf(CQuaternion* pOut, const CQuaternion* pQuat){

	float inv_len = 1.0f / (sqrtf((pQuat->x*pQuat->x) + (pQuat->y*pQuat->y) + (pQuat->z*pQuat->z) + (pQuat->w*pQuat->w)) +  CMath::FLT_EPSILON);

	pOut->x = pQuat->x * inv_len;
	pOut->y = pQuat->y * inv_len;
	pOut->z = pQuat->z * inv_len;
	pOut->w = pQuat->w * inv_len;



}

void  VPCALL CSIMD_Generic:: quaternion_Multiply(CQuaternion* pLeft, const CQuaternion* pRight)
{
	float quat[4];

	quat[0] = pLeft->x;
	quat[1] = pLeft->y;
	quat[2] = pLeft->z;
	quat[3] = pLeft->w;

	pLeft->x = (quat[3] * pRight->x) + (quat[0] * pRight->w) + (quat[1] * pRight->z) - (quat[2] * pRight->y);
	pLeft->y = (quat[3] * pRight->y) - (quat[0] * pRight->z) + (quat[1] * pRight->w) + (quat[2] * pRight->x);
	pLeft->z = (quat[3] * pRight->z) + (quat[0] * pRight->y) - (quat[1] * pRight->x) + (quat[2] * pRight->w);
	pLeft->w = (quat[3] * pRight->w) - (quat[0] * pRight->x) - (quat[1] * pRight->y) - (quat[2] * pRight->z);


}

void  VPCALL CSIMD_Generic:: quaternion_MultiplyOf(CQuaternion* pOut, const CQuaternion* pLeft, const CQuaternion* pRight)
{

	pOut->x = (pLeft->w * pRight->x) + (pLeft->x * pRight->w) + (pLeft->y * pRight->z) - (pLeft->z * pRight->y);
	pOut->y = (pLeft->w * pRight->y) - (pLeft->x * pRight->z) + (pLeft->y * pRight->w) + (pLeft->z * pRight->x);
	pOut->z = (pLeft->w * pRight->z) + (pLeft->x * pRight->y) - (pLeft->y * pRight->x) + (pLeft->z * pRight->w);
	pOut->w = (pLeft->w * pRight->w) - (pLeft->x * pRight->x) - (pLeft->y * pRight->y) - (pLeft->z * pRight->z);


}

void VPCALL CSIMD_Generic::quat_to_mat3x4(sf_m128 q, sf_m128 t, sf_m128 *m){
	CVec3D vec= *reinterpret_cast<CVec3D *>(&t);
	CQuaternion quat= *reinterpret_cast<CQuaternion *>(&q);
	CMatJoint3x4 mat(quat,vec);
	m = reinterpret_cast<sf_m128 *>(&mat) ;
}

sf_m128 VPCALL CSIMD_Generic::quat_transform_vec4(sf_m128 quat, sf_m128 vec){
	CMat3D  mat = *reinterpret_cast<CMat3D *>(&quat);
	CVec3D vec3 = *reinterpret_cast<CVec3D *>(&vec) ;
	CVec3D vec4 = mat * vec3;
	return *reinterpret_cast<sf_m128 *>(&vec4) ;

}
sf_m128 VPCALL CSIMD_Generic::quat_mul_quat(sf_m128 q1, sf_m128 q2){

CQuaternion one=	*reinterpret_cast<CQuaternion *>(&q1) ;
CQuaternion a=	*reinterpret_cast<CQuaternion *>(&q2) ;

CQuaternion quat(	one.w*a.x + one.x*a.w + one.y*a.z - one.z*a.y,
					one.w*a.y + one.y*a.w + one.z*a.x - one.x*a.z,
					one.w*a.z + one.z*a.w + one.x*a.y - one.y*a.x,
					one.w*a.w - one.x*a.x - one.y*a.y - one.z*a.z );
return *reinterpret_cast<sf_m128 *>(&quat) ;
}
sf_m128 VPCALL CSIMD_Generic::quat_div_quat(sf_m128 q1, sf_m128 q2){
	CQuaternion one=	*reinterpret_cast<CQuaternion *>(&q1) ;
CQuaternion r=	*reinterpret_cast<CQuaternion *>(&q2) ;

	CQuaternion quat(one.x*r.w - one.y*r.z + one.z*r.y - one.w*r.x,
	            one.x*r.z + one.y*r.w - one.z*r.x - one.w*r.y,
	           -one.x*r.y + one.y*r.x + one.z*r.w - one.w*r.z,
	            one.x*r.x + one.y*r.y + one.z*r.z + one.w*r.w);
	return *reinterpret_cast<sf_m128 *>(&quat) ;
}






//==================CPlane=================================

void CSIMD_Generic:: plane_FromPoints(CPlane* pOut, const CVec3D* pA, const CVec3D* pB,const CVec3D* pC)
{

	CVec3D v0, v1, cross;

	vector3D_DiffOf(&v0, pA, pB);
	vector3D_DiffOf(&v1, pC, pB);


	vector3D_CrossOf(&cross, &v0, &v1);
	/*
			Calculate d...
			1) ax + by + cz + d = 0
			2) ax + by + cz = -d
			3) -(ax+by+cz) = d
			-dot(p,v) = d
	*/

	pOut->a = cross.x;
	pOut->b = cross.y;
	pOut->c = cross.z;

	/* solve for d with scalar product */
	pOut->d = -(vector3D_Dot(&cross, pA));

}

float CSIMD_Generic:: plane_DistToPoint(const CPlane* pPlane, const CVec3D* pPoint)
{

	return fabsf((pPlane->a*pPoint->x) + (pPlane->b*pPoint->y) + (pPlane->c*pPoint->z) + pPlane->d);

}

float CSIMD_Generic:: plane_Dot(const CPlane* pPlane, const CVec3D* pVec)
{


	return (pPlane->a * pVec->x) + (pPlane->b * pVec->y) + (pPlane->c * pVec->z) + pPlane->d;

}

float CSIMD_Generic:: plane_Dot4(const CPlane* pPlane, const CVec4D* pVec4)
{

	return ((pPlane->a *(*pVec4)[0] ) + (pPlane->b *(*pVec4)[1]) + (pPlane->c *(*pVec4)[2]) + (pPlane->d *(*pVec4)[3]));

}

float CSIMD_Generic:: plane_DotNormal(const CPlane* pPlane, const CVec3D* pVec)
{

	return (pPlane->a * pVec->x) + (pPlane->b * pVec->y) + (pPlane->c * pVec->z);


}

float CSIMD_Generic:: plane_DotPlane(const CPlane* pPlane1, const CPlane* pPlane2)
{


	/*
		equation: cos(theta) = dot(P1, P2) / (mag(P1)*mag(p2))

		Thus, the only way the below code works if mag(P1) == mag(P2) == 1.0f
	*/

	return (pPlane1->a*pPlane2->a) + (pPlane1->b*pPlane2->b) + (pPlane1->c*pPlane2->c);

}

void CSIMD_Generic:: plane_Normalize(CPlane* pOut)
{
	float invmag = 1.0f / sqrtf((pOut->a * pOut->a) + (pOut->b * pOut->b) + (pOut->c * pOut->c));

	/* Note that you have to also divide the distance from origin (d) by magnitude of vector as well */
	pOut->a *= invmag;
	pOut->b *= invmag;
	pOut->c *= invmag;
	pOut->d *= invmag;

}

void CSIMD_Generic::plane_NormalizeOf(CPlane* pOut, CPlane* pIn)
{
	float invmag = 1.0f / sqrtf((pIn->a * pIn->a) + (pIn->b * pIn->b) + (pIn->c * pIn->c));

	/* Note that you have to also divide the distance from origin (d) by magnitude of vector as well */
	pOut->a = pIn->a * invmag;
	pOut->b = pIn->b * invmag;
	pOut->c = pIn->c * invmag;
	pOut->d = pIn->d * invmag;

}

//=====================AABBox===================================

bool  VPCALL  CSIMD_Generic::intersectLineAABB(const CAABBox &box, const CVec4D &rayPos, const CVec4D &rayDir, float tNear, float tFar)
{
	CVec3D linePos=rayPos.toVec3();
	CVec3D lineDir=rayDir.toVec3();
	//s SMF_ASSERT(lineDir.isNormalized());
	//s SMF_ASSERT(tNear <= tFar && "CAABBox::intersectLineAABB: User gave a degenerate line as input for the intersection test!");
	// The user should have inputted values for tNear and tFar to specify the desired subrange [tNear, tFar] of the line
	// for this intersection test.
	// For a Line-CAABBox test, pass in
	//    tNear = -CMath::INFINITY_FLOAT;
	//    tFar = CMath::INFINITY_FLOAT;
	// For a Ray-CAABBox test, pass in
	//    tNear = 0.f;
	//    tFar = CMath::INFINITY_FLOAT;
	// For a LineSegment-CAABBox test, pass in
	//    tNear = 0.f;
	//    tFar = LineSegment.getLenght();

	// Test each cardinal plane (X, Y and Z) in turn.
	if (!CMath::equalsAbs(lineDir.x, 0.f))
	{
		float recipDir = CMath::recipFast(lineDir.x);
		float t1 = (box.minPoint().x - linePos.x) * recipDir;
		float t2 = (box.maxPoint().x - linePos.x) * recipDir;

		// tNear tracks distance to intersect (enter) the CAABBox.
		// tFar tracks the distance to exit the CAABBox.
		if (t1 < t2)
			tNear = MAX(t1, tNear), tFar = MIN(t2, tFar);
		else // Swap t1 and t2.
			tNear = MAX(t2, tNear), tFar = MIN(t1, tFar);

		if (tNear > tFar)
			return false; // Box is missed since we "exit" before entering it.
	}
	else if (linePos.x < box.minPoint().x || linePos.x > box.maxPoint().x)
		return false; // The ray can't possibly enter the box, abort.

	if (!CMath::equalsAbs(lineDir.y, 0.f))
	{
		float recipDir = CMath::recipFast(lineDir.y);
		float t1 = (box.minPoint().y - linePos.y) * recipDir;
		float t2 = (box.maxPoint().y - linePos.y) * recipDir;

		if (t1 < t2)
			tNear = MAX(t1, tNear), tFar = MIN(t2, tFar);
		else // Swap t1 and t2.
			tNear = MAX(t2, tNear), tFar = MIN(t1, tFar);

		if (tNear > tFar)
			return false; // Box is missed since we "exit" before entering it.
	}
	else if (linePos.y < box.minPoint().y || linePos.y > box.maxPoint().y)
		return false; // The ray can't possibly enter the box, abort.

	if (!CMath::equalsAbs(lineDir.z, 0.f)) // ray is parallel to plane in question
	{
		float recipDir = CMath::recipFast(lineDir.z);
		float t1 = (box.minPoint().z - linePos.z) * recipDir;
		float t2 = (box.maxPoint().z - linePos.z) * recipDir;

		if (t1 < t2)
			tNear =MAX(t1, tNear), tFar = MIN(t2, tFar);
		else // Swap t1 and t2.
			tNear = MAX(t2, tNear), tFar = MIN(t1, tFar);
	}
	else if (linePos.z < box.minPoint().z || linePos.z > box.maxPoint().z)
		return false; // The ray can't possibly enter the box, abort.

	return tNear <= tFar;
}


} //end MATH

} //end SMF
