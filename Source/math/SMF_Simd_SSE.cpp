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
#include "math/SMF_Math.h"
#include "math/SMF_Matriz.h"
#include "math/SMF_Vector.h"
#include "math/SMF_SimdGeneric.h"
#include "math/SMF_SimdMMX.h"
#include "math/SMF_SimdSSE.h"
#include "geometry/SMF_Plane.h"
#include "geometry/SMF_DrawVert.h"
#include "math/SMF_JointTransform.h"
#include "util/SMF_Debug.h"

namespace SMF {
namespace MATH{

ALIGN4_INIT1( unsigned long SIMD_DW_mat2quatShuffle0, (3<<0)|(2<<8)|(1<<16)|(0<<24) );
ALIGN4_INIT1( unsigned long SIMD_DW_mat2quatShuffle1, (0<<0)|(1<<8)|(2<<16)|(3<<24) );
ALIGN4_INIT1( unsigned long SIMD_DW_mat2quatShuffle2, (1<<0)|(0<<8)|(3<<16)|(2<<24) );
ALIGN4_INIT1( unsigned long SIMD_DW_mat2quatShuffle3, (2<<0)|(3<<8)|(0<<16)|(1<<24) );

ALIGN4_INIT4( unsigned int  _SIMDx86_float_POSPOSPOSNEG, 0x00000000, 0x00000000, 0x00000000, 0x80000000 );
ALIGN4_INIT4( unsigned int  _SIMDx86_float_SSE_NO_W_MASK, 0xFFFFFFFF,  0xFFFFFFFF, 0xFFFFFFFF, 0x00000000 );
ALIGN4_INIT4( unsigned int  _SIMDx86_float_SSE_NO_XYZ_MASK, 0x00000000, 0x00000000, 0x00000000, 0xFFFFFFFF );

ALIGN4_INIT4( unsigned int POSNEGPOSNEG, 0x00000000, 0x80000000, 0x00000000, 0x80000000 );
ALIGN4_INIT4( unsigned int POSPOSNEGNEG, 0x00000000, 0x00000000, 0x80000000, 0x80000000 );
ALIGN4_INIT4( unsigned int NEGPOSPOSNEG, 0x80000000, 0x00000000, 0x00000000, 0x80000000 );

//(1 << 31 ) == 0x80000000  ==
//(1 << 23 ) == 0x00800000  ==
//(1 << 15 ) == 0x00008000  ==
ALIGN4_INIT4( unsigned int  _SIMDx86_float_NEGPOSPOSPOS, 0x80000000, 0x00000000, 0x00000000, 0x00000000 );

ALIGN4_INIT4( unsigned long SIMD_SP_singleSignBitMask, (unsigned long) ( 1 << 31 ), 0, 0, 0 );
ALIGN4_INIT1( unsigned long SIMD_SP_signBitMask, (unsigned long) ( 1 << 31 ) );
ALIGN4_INIT1( unsigned long SIMD_SP_absMask, (unsigned long) ~( 1 << 31 ) );
ALIGN4_INIT1( unsigned long SIMD_SP_infinityMask, (unsigned long) ~( 1 << 23 ) );
ALIGN4_INIT1( unsigned long SIMD_SP_not, 0xFFFFFFFF );

ALIGN4_INIT1( float SIMD_SP_zero, 0.0f );
ALIGN4_INIT1( float SIMD_SP_half, 0.5f );
ALIGN4_INIT1( float SIMD_SP_one, 1.0f );
ALIGN4_INIT1( float SIMD_SP_two, 2.0f );
ALIGN4_INIT1( float SIMD_SP_three, 3.0f );
ALIGN4_INIT1( float SIMD_SP_four, 4.0f );
ALIGN4_INIT1( float SIMD_SP_maxShort, (1<<15) );
ALIGN4_INIT1( float SIMD_SP_tiny, 1e-10f );
ALIGN4_INIT1( float SIMD_SP_PI, CMath::PI );
ALIGN4_INIT1( float SIMD_SP_halfPI, CMath::HALF_PI );
ALIGN4_INIT1( float SIMD_SP_twoPI, CMath::TWO_PI );
ALIGN4_INIT1( float SIMD_SP_oneOverTwoPI, 1.0f / CMath::TWO_PI );
ALIGN4_INIT1( float SIMD_SP_infinity, CMath::INFINITY_FLOAT );
ALIGN4_INIT4( float SIMD_SP_lastOne, 0.0f, 0.0f, 0.0f, 1.0f );

ALIGN4_INIT1( float SIMD_SP_rsqrt_c0,  3.0f );
ALIGN4_INIT1( float SIMD_SP_rsqrt_c1, -0.5f );
ALIGN4_INIT1( float SIMD_SP_mat2quat_rsqrt_c1, -0.5f*0.5f );

ALIGN4_INIT1( float SIMD_SP_sin_c0, -2.39e-08f );
ALIGN4_INIT1( float SIMD_SP_sin_c1,  2.7526e-06f );
ALIGN4_INIT1( float SIMD_SP_sin_c2, -1.98409e-04f );
ALIGN4_INIT1( float SIMD_SP_sin_c3,  8.3333315e-03f );
ALIGN4_INIT1( float SIMD_SP_sin_c4, -1.666666664e-01f );

ALIGN4_INIT1( float SIMD_SP_cos_c0, -2.605e-07f );
ALIGN4_INIT1( float SIMD_SP_cos_c1,  2.47609e-05f );
ALIGN4_INIT1( float SIMD_SP_cos_c2, -1.3888397e-03f );
ALIGN4_INIT1( float SIMD_SP_cos_c3,  4.16666418e-02f );
ALIGN4_INIT1( float SIMD_SP_cos_c4, -4.999999963e-01f );

ALIGN4_INIT1( float SIMD_SP_atan_c0,  0.0028662257f );
ALIGN4_INIT1( float SIMD_SP_atan_c1, -0.0161657367f );
ALIGN4_INIT1( float SIMD_SP_atan_c2,  0.0429096138f );
ALIGN4_INIT1( float SIMD_SP_atan_c3, -0.0752896400f );
ALIGN4_INIT1( float SIMD_SP_atan_c4,  0.1065626393f );
ALIGN4_INIT1( float SIMD_SP_atan_c5, -0.1420889944f );
ALIGN4_INIT1( float SIMD_SP_atan_c6,  0.1999355085f );
ALIGN4_INIT1( float SIMD_SP_atan_c7, -0.3333314528f );

ALIGN4_INIT1(double _SIMDx86_double_one, 1.0);
ALIGN4_INIT1(float  _SIMDx86_float_one, 1.0f);
ALIGN4_INIT4(float  _SIMDx86_float_one_w, 0.0f, 0.0f, 0.0f, 1.0f );
ALIGN4_INIT4(float  _SIMDx86_float_one_y, 0.0f, 1.0f, 0.0f, 0.0f );

ALIGN2_INIT2(unsigned int  _SIMDx86_float_POSNEG, 0x00000000, 0x80000000 );	/* XOR Changes the signs to + - */
ALIGN2_INIT2(unsigned int  _SIMDx86_float_NEGPOS, 0x80000000, 0x00000000 );	/* XOR Changes the signs to - + */
ALIGN2_INIT2(unsigned int  _SIMDx86_float_NEGNEG, 0x80000000, 0x80000000 );	/* XOR Changes the signs to - - */
ALIGN2_INIT2(unsigned int  _SIMDx86_float_3DNOW_NO_W_MASK, 0xFFFFFFFF, 0x00000000 );
ALIGN4_INIT4(unsigned int  _SIMDx86_float_ABS, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF,0x7FFFFFFF );




//===============================================================
//                                                        M
//  SSE implementation of CSIMDProcessor                MrE
//                                                        E
//===============================================================


#if defined(_WIN32) && defined(_MSC_VER)

/*
============
SSE_InvSqrt
============
*/
float CSIMD_SSE::invSqrt( float x ) {
	float y;

	__asm {
		movss		xmm0, x
		rsqrtss		xmm1, xmm0
		mulss		xmm0, xmm1
		mulss		xmm0, xmm1
		subss		xmm0, SIMD_SP_rsqrt_c0
		mulss		xmm1, SIMD_SP_rsqrt_c1
		mulss		xmm0, xmm1
		movss		y, xmm0
	}
	return y;
}

/*
============
SSE_InvSqrt4
============
*/
void CSIMD_SSE::InvSqrt4( float x[4] ) {
	__asm {
		mov			edi, x
		movaps		xmm0, [edi]
		rsqrtps		xmm1, xmm0
		mulps		xmm0, xmm1
		mulps		xmm0, xmm1
		subps		xmm0, SIMD_SP_rsqrt_c0
		mulps		xmm1, SIMD_SP_rsqrt_c1
		mulps		xmm0, xmm1
		movaps		[edi], xmm0
	}
}

/*
============
SSE_SinZeroHalfPI

  The angle must be between zero and half PI. (0 a 90 graus)
  \return o seno do ângulo passado
  \param x ângulo em radianos que se deseja calcular o seno
============
*/
float CSIMD_SSE::sinZeroHalfPI( float a ) {
#if 1

	float t;

	SMF_ASSERT( a >= 0.0f && a <= CMath::HALF_PI );

	__asm {
		movss		xmm0, a
		movss		xmm1, xmm0
		mulss		xmm1, xmm1
		movss		xmm2, SIMD_SP_sin_c0
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_sin_c1
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_sin_c2
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_sin_c3
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_sin_c4
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_one
		mulss		xmm2, xmm0
		movss		t, xmm2
	}

	return t;

#else

	float s, t;

	SMF_ASSERT( a >= 0.0f && a <= CMath::HALF_PI );

	s = a * a;
	t = -2.39e-08f;
	t *= s;
	t += 2.7526e-06f;
	t *= s;
	t += -1.98409e-04f;
	t *= s;
	t += 8.3333315e-03f;
	t *= s;
	t += -1.666666664e-01f;
	t *= s;
	t += 1.0f;
	t *= a;

	return t;

#endif
}

/*
============
SSE_Sin4ZeroHalfPI

  The angle must be between zero and half PI.

  \param a entrada
  \param s saida
============
*/
void CSIMD_SSE::sin4ZeroHalfPI( float a[4], float s[4] ) {
	__asm {
		mov			edi, a
		mov			esi, s
		movaps		xmm0, [edi]
		movaps		xmm1, xmm0
		mulps		xmm1, xmm1
		movaps		xmm2, SIMD_SP_sin_c0
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_sin_c1
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_sin_c2
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_sin_c3
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_sin_c4
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_one
		mulps		xmm2, xmm0
		movaps		[esi], xmm2
	}
}

/*
============
SSE_Sin
\param a angulo em radianos
============
*/
float CSIMD_SSE::sin( float a ) {
#if 1

	float t;

	__asm {
		movss		xmm1, a
		movss		xmm2, xmm1
		movss		xmm3, xmm1
		mulss		xmm2, SIMD_SP_oneOverTwoPI
		cvttss2si	ecx, xmm2
		cmpltss		xmm3, SIMD_SP_zero
		andps		xmm3, SIMD_SP_one
		cvtsi2ss	xmm2, ecx
		subss		xmm2, xmm3
		mulss		xmm2, SIMD_SP_twoPI
		subss		xmm1, xmm2

		movss		xmm0, SIMD_SP_PI			// xmm0 = PI
		subss		xmm0, xmm1					// xmm0 = PI - a
		movss		xmm1, xmm0					// xmm1 = PI - a
		andps		xmm1, SIMD_SP_signBitMask	// xmm1 = signbit( PI - a )
		movss		xmm2, xmm0					// xmm2 = PI - a
		xorps		xmm2, xmm1					// xmm2 = fabs( PI - a )
		cmpnltss	xmm2, SIMD_SP_halfPI		// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? 0xFFFFFFFF : 0x00000000
		movss		xmm3, SIMD_SP_PI			// xmm3 = PI
		xorps		xmm3, xmm1					// xmm3 = PI ^ signbit( PI - a )
		andps		xmm3, xmm2					// xmm3 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? ( PI ^ signbit( PI - a ) ) : 0.0f
		andps		xmm2, SIMD_SP_signBitMask	// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? SIMD_SP_signBitMask : 0.0f
		xorps		xmm0, xmm2
		addps		xmm0, xmm3

		movss		xmm1, xmm0
		mulss		xmm1, xmm1
		movss		xmm2, SIMD_SP_sin_c0
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_sin_c1
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_sin_c2
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_sin_c3
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_sin_c4
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_one
		mulss		xmm2, xmm0
		movss		t, xmm2
	}

	return t;

#else

	float s, t;

	if ( ( a < 0.0f ) || ( a >= CMath::TWO_PI ) ) {
		a -= floorf( a / CMath::TWO_PI ) * CMath::TWO_PI;
	}

	a = CMath::PI - a;
	if ( fabs( a ) >= CMath::HALF_PI ) {
		a = ( ( a < 0.0f ) ? -CMath::PI : CMath::PI ) - a;
	}

	s = a * a;
	t = -2.39e-08f;
	t *= s;
	t += 2.7526e-06f;
	t *= s;
	t += -1.98409e-04f;
	t *= s;
	t += 8.3333315e-03f;
	t *= s;
	t += -1.666666664e-01f;
	t *= s;
	t += 1.0f;
	t *= a;

	return t;

#endif
}

/*
============
SSE_Sin4
============
*/
void CSIMD_SSE::sin4( float a[4], float s[4] ) {
	__asm {
		mov			edi, a
		mov			esi, s
		movaps		xmm1, [edi]
		movaps		xmm2, xmm1
		mulps		xmm2, SIMD_SP_oneOverTwoPI
		movhlps		xmm3, xmm2
		cvttss2si	ecx, xmm2
		cvtsi2ss	xmm2, ecx
		cvttss2si	edx, xmm3
		cvtsi2ss	xmm3, edx
		shufps		xmm2, xmm2, R_SHUFFLEPS( 1, 0, 0, 0 )
		shufps		xmm3, xmm3, R_SHUFFLEPS( 1, 0, 0, 0 )
		cvttss2si	ecx, xmm2
		cvtsi2ss	xmm2, ecx
		cvttss2si	edx, xmm3
		cvtsi2ss	xmm3, edx
		shufps		xmm2, xmm3, R_SHUFFLEPS( 1, 0, 1, 0 )
		movaps		xmm3, xmm1
		cmpltps		xmm3, SIMD_SP_zero
		andps		xmm3, SIMD_SP_one
		subps		xmm2, xmm3
		mulps		xmm2, SIMD_SP_twoPI
		subps		xmm1, xmm2

		movaps		xmm0, SIMD_SP_PI			// xmm0 = PI
		subps		xmm0, xmm1					// xmm0 = PI - a
		movaps		xmm1, xmm0					// xmm1 = PI - a
		andps		xmm1, SIMD_SP_signBitMask	// xmm1 = signbit( PI - a )
		movaps		xmm2, xmm0					// xmm2 = PI - a
		xorps		xmm2, xmm1					// xmm2 = fabs( PI - a )
		cmpnltps	xmm2, SIMD_SP_halfPI		// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? 0xFFFFFFFF : 0x00000000
		movaps		xmm3, SIMD_SP_PI			// xmm3 = PI
		xorps		xmm3, xmm1					// xmm3 = PI ^ signbit( PI - a )
		andps		xmm3, xmm2					// xmm3 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? ( PI ^ signbit( PI - a ) ) : 0.0f
		andps		xmm2, SIMD_SP_signBitMask	// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? SIMD_SP_signBitMask : 0.0f
		xorps		xmm0, xmm2
		addps		xmm0, xmm3

		movaps		xmm1, xmm0
		mulps		xmm1, xmm1
		movaps		xmm2, SIMD_SP_sin_c0
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_sin_c1
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_sin_c2
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_sin_c3
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_sin_c4
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_one
		mulps		xmm2, xmm0
		movaps		[esi], xmm2
	}
}

/*
============
SSE_CosZeroHalfPI

  The angle must be between zero and half PI.
============
*/
float CSIMD_SSE::cosZeroHalfPI( float a ) {
#if 1

	float t;

	SMF_ASSERT( a >= 0.0f && a <= CMath::HALF_PI );

	__asm {
		movss		xmm0, a
		mulss		xmm0, xmm0
		movss		xmm1, SIMD_SP_cos_c0
		mulss		xmm1, xmm0
		addss		xmm1, SIMD_SP_cos_c1
		mulss		xmm1, xmm0
		addss		xmm1, SIMD_SP_cos_c2
		mulss		xmm1, xmm0
		addss		xmm1, SIMD_SP_cos_c3
		mulss		xmm1, xmm0
		addss		xmm1, SIMD_SP_cos_c4
		mulss		xmm1, xmm0
		addss		xmm1, SIMD_SP_one
		movss		t, xmm1
	}

	return t;

#else

	float s, t;

	SMF_ASSERT( a >= 0.0f && a <= CMath::HALF_PI );

	s = a * a;
	t = -2.605e-07f;
	t *= s;
	t += 2.47609e-05f;
	t *= s;
	t += -1.3888397e-03f;
	t *= s;
	t += 4.16666418e-02f;
	t *= s;
	t += -4.999999963e-01f;
	t *= s;
	t += 1.0f;

	return t;

#endif
}

/*
============
SSE_Cos4ZeroHalfPI

  The angle must be between zero and half PI.
============
*/
void CSIMD_SSE::cos4ZeroHalfPI( float a[4], float c[4] ) {
	__asm {
		mov			edi, a
		mov			esi, c
		movaps		xmm0, [edi]
		mulps		xmm0, xmm0
		movaps		xmm1, SIMD_SP_cos_c0
		mulps		xmm1, xmm0
		addps		xmm1, SIMD_SP_cos_c1
		mulps		xmm1, xmm0
		addps		xmm1, SIMD_SP_cos_c2
		mulps		xmm1, xmm0
		addps		xmm1, SIMD_SP_cos_c3
		mulps		xmm1, xmm0
		addps		xmm1, SIMD_SP_cos_c4
		mulps		xmm1, xmm0
		addps		xmm1, SIMD_SP_one
		movaps		[esi], xmm2
	}
}

/*
============
SSE_Cos
============
*/
float CSIMD_SSE::cos( float a ) {
#if 1

	float t;

	__asm {
		movss		xmm1, a
		movss		xmm2, xmm1
		movss		xmm3, xmm1
		mulss		xmm2, SIMD_SP_oneOverTwoPI
		cvttss2si	ecx, xmm2
		cmpltss		xmm3, SIMD_SP_zero
		andps		xmm3, SIMD_SP_one
		cvtsi2ss	xmm2, ecx
		subss		xmm2, xmm3
		mulss		xmm2, SIMD_SP_twoPI
		subss		xmm1, xmm2

		movss		xmm0, SIMD_SP_PI			// xmm0 = PI
		subss		xmm0, xmm1					// xmm0 = PI - a
		movss		xmm1, xmm0					// xmm1 = PI - a
		andps		xmm1, SIMD_SP_signBitMask	// xmm1 = signbit( PI - a )
		movss		xmm2, xmm0					// xmm2 = PI - a
		xorps		xmm2, xmm1					// xmm2 = fabs( PI - a )
		cmpnltss	xmm2, SIMD_SP_halfPI		// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? 0xFFFFFFFF : 0x00000000
		movss		xmm3, SIMD_SP_PI			// xmm3 = PI
		xorps		xmm3, xmm1					// xmm3 = PI ^ signbit( PI - a )
		andps		xmm3, xmm2					// xmm3 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? ( PI ^ signbit( PI - a ) ) : 0.0f
		andps		xmm2, SIMD_SP_signBitMask	// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? SIMD_SP_signBitMask : 0.0f
		xorps		xmm0, xmm2
		addps		xmm0, xmm3

		mulss		xmm0, xmm0
		movss		xmm1, SIMD_SP_cos_c0
		mulss		xmm1, xmm0
		addss		xmm1, SIMD_SP_cos_c1
		mulss		xmm1, xmm0
		addss		xmm1, SIMD_SP_cos_c2
		mulss		xmm1, xmm0
		addss		xmm1, SIMD_SP_cos_c3
		mulss		xmm1, xmm0
		addss		xmm1, SIMD_SP_cos_c4
		mulss		xmm1, xmm0
		addss		xmm1, SIMD_SP_one
		xorps		xmm2, SIMD_SP_signBitMask
		xorps		xmm1, xmm2
		movss		t, xmm1
	}

	return t;

#else

	float s, t;

	if ( ( a < 0.0f ) || ( a >= CMath::TWO_PI ) ) {
		a -= floorf( a / CMath::TWO_PI ) * CMath::TWO_PI;
	}

	a = CMath::PI - a;
	if ( fabs( a ) >= CMath::HALF_PI ) {
		a = ( ( a < 0.0f ) ? -CMath::PI : CMath::PI ) - a;
		d = 1.0f;
	} else {
		d = -1.0f;
	}

	s = a * a;
	t = -2.605e-07f;
	t *= s;
	t += 2.47609e-05f;
	t *= s;
	t += -1.3888397e-03f;
	t *= s;
	t += 4.16666418e-02f;
	t *= s;
	t += -4.999999963e-01f;
	t *= s;
	t += 1.0f;
	t *= d;

	return t;

#endif
}

/*
============
SSE_Cos4
============
*/
void CSIMD_SSE::cos4( float a[4], float c[4] ) {
	__asm {
		mov			edi, a
		mov			esi, c
		movaps		xmm1, [edi]
		movaps		xmm2, xmm1
		mulps		xmm2, SIMD_SP_oneOverTwoPI
		movhlps		xmm3, xmm2
		cvttss2si	ecx, xmm2
		cvtsi2ss	xmm2, ecx
		cvttss2si	edx, xmm3
		cvtsi2ss	xmm3, edx
		shufps		xmm2, xmm2, R_SHUFFLEPS( 1, 0, 0, 0 )
		shufps		xmm3, xmm3, R_SHUFFLEPS( 1, 0, 0, 0 )
		cvttss2si	ecx, xmm2
		cvtsi2ss	xmm2, ecx
		cvttss2si	edx, xmm3
		cvtsi2ss	xmm3, edx
		shufps		xmm2, xmm3, R_SHUFFLEPS( 1, 0, 1, 0 )
		movaps		xmm3, xmm1
		cmpltps		xmm3, SIMD_SP_zero
		andps		xmm3, SIMD_SP_one
		subps		xmm2, xmm3
		mulps		xmm2, SIMD_SP_twoPI
		subps		xmm1, xmm2

		movaps		xmm0, SIMD_SP_PI			// xmm0 = PI
		subps		xmm0, xmm1					// xmm0 = PI - a
		movaps		xmm1, xmm0					// xmm1 = PI - a
		andps		xmm1, SIMD_SP_signBitMask	// xmm1 = signbit( PI - a )
		movaps		xmm2, xmm0					// xmm2 = PI - a
		xorps		xmm2, xmm1					// xmm2 = fabs( PI - a )
		cmpnltps	xmm2, SIMD_SP_halfPI		// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? 0xFFFFFFFF : 0x00000000
		movaps		xmm3, SIMD_SP_PI			// xmm3 = PI
		xorps		xmm3, xmm1					// xmm3 = PI ^ signbit( PI - a )
		andps		xmm3, xmm2					// xmm3 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? ( PI ^ signbit( PI - a ) ) : 0.0f
		andps		xmm2, SIMD_SP_signBitMask	// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? SIMD_SP_signBitMask : 0.0f
		xorps		xmm0, xmm2
		addps		xmm0, xmm3

		mulps		xmm0, xmm0
		movaps		xmm1, SIMD_SP_cos_c0
		mulps		xmm1, xmm0
		addps		xmm1, SIMD_SP_cos_c1
		mulps		xmm1, xmm0
		addps		xmm1, SIMD_SP_cos_c2
		mulps		xmm1, xmm0
		addps		xmm1, SIMD_SP_cos_c3
		mulps		xmm1, xmm0
		addps		xmm1, SIMD_SP_cos_c4
		mulps		xmm1, xmm0
		addps		xmm1, SIMD_SP_one
		xorps		xmm2, SIMD_SP_signBitMask
		xorps		xmm1, xmm2
		movaps		[esi], xmm1
	}
}

/*
============
SSE_SinCos
============
*/
void CSIMD_SSE::sincos( float a, float &s, float &c ) {
	__asm {
		mov			edi, s
		mov			esi, c
		movss		xmm1, a
		movss		xmm2, xmm1
		movss		xmm3, xmm1
		mulss		xmm2, SIMD_SP_oneOverTwoPI
		cvttss2si	ecx, xmm2
		cmpltss		xmm3, SIMD_SP_zero
		andps		xmm3, SIMD_SP_one
		cvtsi2ss	xmm2, ecx
		subss		xmm2, xmm3
		mulss		xmm2, SIMD_SP_twoPI
		subss		xmm1, xmm2

		movss		xmm0, SIMD_SP_PI			// xmm0 = PI
		subss		xmm0, xmm1					// xmm0 = PI - a
		movss		xmm1, xmm0					// xmm1 = PI - a
		andps		xmm1, SIMD_SP_signBitMask	// xmm1 = signbit( PI - a )
		movss		xmm2, xmm0					// xmm2 = PI - a
		xorps		xmm2, xmm1					// xmm2 = fabs( PI - a )
		cmpnltss	xmm2, SIMD_SP_halfPI		// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? 0xFFFFFFFF : 0x00000000
		movss		xmm3, SIMD_SP_PI			// xmm3 = PI
		xorps		xmm3, xmm1					// xmm3 = PI ^ signbit( PI - a )
		andps		xmm3, xmm2					// xmm3 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? ( PI ^ signbit( PI - a ) ) : 0.0f
		andps		xmm2, SIMD_SP_signBitMask	// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? SIMD_SP_signBitMask : 0.0f
		xorps		xmm0, xmm2
		addps		xmm0, xmm3

		movss		xmm1, xmm0
		mulss		xmm1, xmm1
		movss		xmm3, SIMD_SP_sin_c0
		movss		xmm4, SIMD_SP_cos_c0
		mulss		xmm3, xmm1
		mulss		xmm4, xmm1
		addss		xmm3, SIMD_SP_sin_c1
		addss		xmm4, SIMD_SP_cos_c1
		mulss		xmm3, xmm1
		mulss		xmm4, xmm1
		addss		xmm3, SIMD_SP_sin_c2
		addss		xmm4, SIMD_SP_cos_c2
		mulss		xmm3, xmm1
		mulss		xmm4, xmm1
		addss		xmm3, SIMD_SP_sin_c3
		addss		xmm4, SIMD_SP_cos_c3
		mulss		xmm3, xmm1
		mulss		xmm4, xmm1
		addss		xmm3, SIMD_SP_sin_c4
		addss		xmm4, SIMD_SP_cos_c4
		mulss		xmm3, xmm1
		mulss		xmm4, xmm1
		addss		xmm3, SIMD_SP_one
		addss		xmm4, SIMD_SP_one
		mulss		xmm3, xmm0
		xorps		xmm2, SIMD_SP_signBitMask
		xorps		xmm4, xmm2
		movss		[edi], xmm2
		movss		[esi], xmm3
	}
}

/*
============
SSE_SinCos4
============
*/
void CSIMD_SSE::sincos4( float a[4], float s[4], float c[4] ) {
	__asm {
		mov			eax, a
		mov			edi, s
		mov			esi, c
		movaps		xmm1, [eax]
		movaps		xmm2, xmm1
		mulps		xmm2, SIMD_SP_oneOverTwoPI
		movhlps		xmm3, xmm2
		cvttss2si	ecx, xmm2
		cvtsi2ss	xmm2, ecx
		cvttss2si	edx, xmm3
		cvtsi2ss	xmm3, edx
		shufps		xmm2, xmm2, R_SHUFFLEPS( 1, 0, 0, 0 )
		shufps		xmm3, xmm3, R_SHUFFLEPS( 1, 0, 0, 0 )
		cvttss2si	ecx, xmm2
		cvtsi2ss	xmm2, ecx
		cvttss2si	edx, xmm3
		cvtsi2ss	xmm3, edx
		shufps		xmm2, xmm3, R_SHUFFLEPS( 1, 0, 1, 0 )
		movaps		xmm3, xmm1
		cmpltps		xmm3, SIMD_SP_zero
		andps		xmm3, SIMD_SP_one
		subps		xmm2, xmm3
		mulps		xmm2, SIMD_SP_twoPI
		subps		xmm1, xmm2

		movaps		xmm0, SIMD_SP_PI			// xmm0 = PI
		subps		xmm0, xmm1					// xmm0 = PI - a
		movaps		xmm1, xmm0					// xmm1 = PI - a
		andps		xmm1, SIMD_SP_signBitMask	// xmm1 = signbit( PI - a )
		movaps		xmm2, xmm0					// xmm2 = PI - a
		xorps		xmm2, xmm1					// xmm2 = fabs( PI - a )
		cmpnltps	xmm2, SIMD_SP_halfPI		// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? 0xFFFFFFFF : 0x00000000
		movaps		xmm3, SIMD_SP_PI			// xmm3 = PI
		xorps		xmm3, xmm1					// xmm3 = PI ^ signbit( PI - a )
		andps		xmm3, xmm2					// xmm3 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? ( PI ^ signbit( PI - a ) ) : 0.0f
		andps		xmm2, SIMD_SP_signBitMask	// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? SIMD_SP_signBitMask : 0.0f
		xorps		xmm0, xmm2
		addps		xmm0, xmm3

		movaps		xmm0, [eax]
		movaps		xmm1, xmm0
		mulps		xmm1, xmm1
		movaps		xmm3, SIMD_SP_sin_c0
		movaps		xmm4, SIMD_SP_cos_c0
		mulps		xmm3, xmm1
		mulps		xmm4, xmm1
		addps		xmm3, SIMD_SP_sin_c1
		addps		xmm4, SIMD_SP_cos_c1
		mulps		xmm3, xmm1
		mulps		xmm4, xmm1
		addps		xmm3, SIMD_SP_sin_c2
		addps		xmm4, SIMD_SP_cos_c2
		mulps		xmm3, xmm1
		mulps		xmm4, xmm1
		addps		xmm3, SIMD_SP_sin_c3
		addps		xmm4, SIMD_SP_cos_c3
		mulps		xmm3, xmm1
		mulps		xmm4, xmm1
		addps		xmm3, SIMD_SP_sin_c4
		addps		xmm4, SIMD_SP_cos_c4
		mulps		xmm3, xmm1
		mulps		xmm4, xmm1
		addps		xmm3, SIMD_SP_one
		addps		xmm4, SIMD_SP_one
		mulps		xmm3, xmm0
		xorps		xmm2, SIMD_SP_signBitMask
		xorps		xmm4, xmm2
		movaps		[edi], xmm3
		movaps		[esi], xmm4
	}
}

/*
============
SSE_ATanPositive

  Both 'x' and 'y' must be positive.
============
*/
float CSIMD_SSE::aTanPositive( float y, float x ) {
#if 1

	float t;

	SMF_ASSERT( y >= 0.0f && x >= 0.0f );

	__asm {
		movss		xmm0, x
		movss		xmm3, xmm0
		movss		xmm1, y
		minss		xmm0, xmm1
		maxss		xmm1, xmm3
		cmpeqss		xmm3, xmm0
		rcpss		xmm2, xmm1
		mulss		xmm1, xmm2
		mulss		xmm1, xmm2
		addss		xmm2, xmm2
		subss		xmm2, xmm1				// xmm2 = 1 / y or 1 / x
		mulss		xmm0, xmm2				// xmm0 = x / y or y / x
		movss		xmm1, xmm3
		andps		xmm1, SIMD_SP_signBitMask
		xorps		xmm0, xmm1				// xmm0 = -x / y or y / x
		andps		xmm3, SIMD_SP_halfPI	// xmm3 = HALF_PI or 0.0f
		movss		xmm1, xmm0
		mulss		xmm1, xmm1				// xmm1 = s
		movss		xmm2, SIMD_SP_atan_c0
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_atan_c1
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_atan_c2
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_atan_c3
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_atan_c4
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_atan_c5
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_atan_c6
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_atan_c7
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_one
		mulss		xmm2, xmm0
		addss		xmm2, xmm3
		movss		t, xmm2
	}

	return t;

#else

	float a, d, s, t;

	SMF_ASSERT( y >= 0.0f && x >= 0.0f );

	if ( y > x ) {
		a = -x / y;
		d = CMath::HALF_PI;
	} else {
		a = y / x;
		d = 0.0f;
	}
	s = a * a;
	t = 0.0028662257f;
	t *= s;
	t += -0.0161657367f;
	t *= s;
	t += 0.0429096138f;
	t *= s;
	t += -0.0752896400f;
	t *= s;
	t += 0.1065626393f;
	t *= s;
	t += -0.1420889944f;
	t *= s;
	t += 0.1999355085f;
	t *= s;
	t += -0.3333314528f;
	t *= s;
	t += 1.0f;
	t *= a;
	t += d;

	return t;

#endif
}

/*
============
SSE_ATan4Positive

  Both 'x' and 'y' must be positive.
============
*/
void CSIMD_SSE::aTan4Positive( float y[4], float x[4], float at[4] ) {
	__asm {
		mov			esi, x
		mov			edi, y
		mov			edx, at
		movaps		xmm0, [esi]
		movaps		xmm3, xmm0
		movaps		xmm1, [edi]
		minps		xmm0, xmm1
		maxps		xmm1, xmm3
		cmpeqps		xmm3, xmm0
		rcpps		xmm2, xmm1
		mulps		xmm1, xmm2
		mulps		xmm1, xmm2
		addps		xmm2, xmm2
		subps		xmm2, xmm1				// xmm2 = 1 / y or 1 / x
		mulps		xmm0, xmm2				// xmm0 = x / y or y / x
		movaps		xmm1, xmm3
		andps		xmm1, SIMD_SP_signBitMask
		xorps		xmm0, xmm1				// xmm0 = -x / y or y / x
		andps		xmm3, SIMD_SP_halfPI	// xmm3 = HALF_PI or 0.0f
		movaps		xmm1, xmm0
		mulps		xmm1, xmm1				// xmm1 = s
		movaps		xmm2, SIMD_SP_atan_c0
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_atan_c1
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_atan_c2
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_atan_c3
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_atan_c4
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_atan_c5
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_atan_c6
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_atan_c7
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_one
		mulps		xmm2, xmm0
		addps		xmm2, xmm3
		movaps		[edx], xmm2
	}
}

/*
============
SSE_ATan
============
*/
float CSIMD_SSE::atan( float y, float x ) {
#if 1

	float t;

	__asm {
		movss		xmm0, x
		movss		xmm3, xmm0
		movss		xmm4, xmm0
		andps		xmm0, SIMD_SP_absMask
		movss		xmm1, y
		xorps		xmm4, xmm1
		andps		xmm1, SIMD_SP_absMask
		andps		xmm4, SIMD_SP_signBitMask
		minss		xmm0, xmm1
		maxss		xmm1, xmm3
		cmpeqss		xmm3, xmm0
		rcpss		xmm2, xmm1
		mulss		xmm1, xmm2
		mulss		xmm1, xmm2
		addss		xmm2, xmm2
		subss		xmm2, xmm1				// xmm2 = 1 / y or 1 / x
		mulss		xmm0, xmm2				// xmm0 = x / y or y / x
		xorps		xmm0, xmm4
		movss		xmm1, xmm3
		andps		xmm1, SIMD_SP_signBitMask
		xorps		xmm0, xmm1				// xmm0 = -x / y or y / x
		orps		xmm4, SIMD_SP_halfPI	// xmm4 = +/- HALF_PI
		andps		xmm3, xmm4				// xmm3 = +/- HALF_PI or 0.0f
		movss		xmm1, xmm0
		mulss		xmm1, xmm1				// xmm1 = s
		movss		xmm2, SIMD_SP_atan_c0
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_atan_c1
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_atan_c2
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_atan_c3
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_atan_c4
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_atan_c5
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_atan_c6
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_atan_c7
		mulss		xmm2, xmm1
		addss		xmm2, SIMD_SP_one
		mulss		xmm2, xmm0
		addss		xmm2, xmm3
		movss		t, xmm2
	}

	return t;

#else

	float a, d, s, t;

	if ( fabs( y ) > fabs( x ) ) {
		a = -x / y;
		d = CMath::HALF_PI;
		*((unsigned long *)&d) ^= ( *((unsigned long *)&x) ^ *((unsigned long *)&y) ) & (1<<31);
	} else {
		a = y / x;
		d = 0.0f;
	}

	s = a * a;
	t = 0.0028662257f;
	t *= s;
	t += -0.0161657367f;
	t *= s;
	t += 0.0429096138f;
	t *= s;
	t += -0.0752896400f;
	t *= s;
	t += 0.1065626393f;
	t *= s;
	t += -0.1420889944f;
	t *= s;
	t += 0.1999355085f;
	t *= s;
	t += -0.3333314528f;
	t *= s;
	t += 1.0f;
	t *= a;
	t += d;

	return t;

#endif
}

/*
============
SSE_ATan4
============
*/
void CSIMD_SSE::aTan4( float y[4], float x[4], float at[4] ) {
	__asm {
		mov			esi, x
		mov			edi, y
		mov			edx, at
		movaps		xmm0, [esi]
		movaps		xmm3, xmm0
		movaps		xmm4, xmm0
		andps		xmm0, SIMD_SP_absMask
		movaps		xmm1, [edi]
		xorps		xmm4, xmm1
		andps		xmm1, SIMD_SP_absMask
		andps		xmm4, SIMD_SP_signBitMask
		minps		xmm0, xmm1
		maxps		xmm1, xmm3
		cmpeqps		xmm3, xmm0
		rcpps		xmm2, xmm1
		mulps		xmm1, xmm2
		mulps		xmm1, xmm2
		addps		xmm2, xmm2
		subps		xmm2, xmm1				// xmm2 = 1 / y or 1 / x
		mulps		xmm0, xmm2				// xmm0 = x / y or y / x
		xorps		xmm0, xmm4
		movaps		xmm1, xmm3
		andps		xmm1, SIMD_SP_signBitMask
		xorps		xmm0, xmm1				// xmm0 = -x / y or y / x
		orps		xmm4, SIMD_SP_halfPI	// xmm4 = +/- HALF_PI
		andps		xmm3, xmm4				// xmm3 = +/- HALF_PI or 0.0f
		movaps		xmm1, xmm0
		mulps		xmm1, xmm1				// xmm1 = s
		movaps		xmm2, SIMD_SP_atan_c0
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_atan_c1
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_atan_c2
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_atan_c3
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_atan_c4
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_atan_c5
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_atan_c6
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_atan_c7
		mulps		xmm2, xmm1
		addps		xmm2, SIMD_SP_one
		mulps		xmm2, xmm0
		addps		xmm2, xmm3
		movaps		[edx], xmm2
	}
}


/*
============
CSIMD_SSE::getName
============
*/
const char * CSIMD_SSE::getName() const {
	return "MMX & SSE";
}

/*
============
CSIMD_SSE::add

  dst[i] = constant + src[i];
============
*/
void VPCALL CSIMD_SSE::add( float *dst, const float constant, const float *src, const int count ) {
	KFLOAT_CA( add, dst, src, constant, count )
}

/*
============
CSIMD_SSE::add

  dst[i] = src0[i] + src1[i];
============
*/
void VPCALL CSIMD_SSE::add( float *dst, const float *src0, const float *src1, const int count ) {
	KFLOAT_AA( add, dst, src0, src1, count )
}

/*
============
CSIMD_SSE::sub

  dst[i] = constant - src[i];
============
*/
void VPCALL CSIMD_SSE::sub( float *dst, const float constant, const float *src, const int count ) {
	KFLOAT_CA( sub, dst, src, constant, count )
}

/*
============
CSIMD_SSE::sub

  dst[i] = src0[i] - src1[i];
============
*/
void VPCALL CSIMD_SSE::sub( float *dst, const float *src0, const float *src1, const int count ) {
	KFLOAT_AA( sub, dst, src0, src1, count )
}

/*
============
CSIMD_SSE::mul

  dst[i] = constant * src[i];
============
*/
void VPCALL CSIMD_SSE::mul( float *dst, const float constant, const float *src, const int count ) {
	KFLOAT_CA( mul, dst, src, constant, count )
}

/*
============
CSIMD_SSE::mul

  dst[i] = src0[i] * src1[i];
============
*/
void VPCALL CSIMD_SSE::mul( float *dst, const float *src0, const float *src1, const int count ) {
	KFLOAT_AA( mul, dst, src0, src1, count )
}

/*
============
CSIMD_SSE::div

  dst[i] = constant / src[i];
============
*/
void VPCALL CSIMD_SSE::div( float *dst, const float constant, const float *src, const int count ) {
	int pre, post;

	//	1 / x = 2 * rcpps(x) - (x * rcpps(x) * rcpps(x));
	__asm
	{
		movss	xmm1,constant
		shufps	xmm1,xmm1,0

		KFLOATINITDS( dst, src, count, pre, post )
		and		eax,15
		jne		lpNA
		jmp		lpA
		align	16
lpA:
		movaps	xmm2,[edx+ebx]
		movaps	xmm3,[edx+ebx+16]
		rcpps	xmm4,xmm2
		rcpps	xmm5,xmm3
		prefetchnta	[edx+ebx+64]
		mulps	xmm2,xmm4
		mulps	xmm2,xmm4
		mulps	xmm3,xmm5
		mulps	xmm3,xmm5
		addps	xmm4,xmm4
		addps	xmm5,xmm5
		subps	xmm4,xmm2
		subps	xmm5,xmm3
		mulps	xmm4,xmm1
		mulps	xmm5,xmm1
		movaps	[edi+ebx],xmm4
		movaps	[edi+ebx+16],xmm5
		add		ebx,16*2
		jl		lpA
		jmp		done
		align	16
lpNA:
		movups	xmm2,[edx+ebx]
		movups	xmm3,[edx+ebx+16]
		rcpps	xmm4,xmm2
		rcpps	xmm5,xmm3
		prefetchnta	[edx+ebx+64]
		mulps	xmm2,xmm4
		mulps	xmm2,xmm4
		mulps	xmm3,xmm5
		mulps	xmm3,xmm5
		addps	xmm4,xmm4
		addps	xmm5,xmm5
		subps	xmm4,xmm2
		subps	xmm5,xmm3
		mulps	xmm4,xmm1
		mulps	xmm5,xmm1
		movaps	[edi+ebx],xmm4
		movaps	[edi+ebx+16],xmm5
		add		ebx,16*2
		jl		lpNA
done:
		mov		edx,src
		mov		edi,dst
		KFLOATOPER( KDIVDSS1( [edi+ebx],xmm1,[edx+ebx] ),
					KDIVDSS4( [edi+ebx],xmm1,[edx+ebx] ), count )
	}
}

/*
============
CSIMD_SSE::div

  dst[i] = src0[i] / src1[i];
============
*/
void VPCALL CSIMD_SSE::div( float *dst, const float *src0, const float *src1, const int count ) {
	int		pre,post;

	//	1 / x = 2 * rcpps(x) - (x * rcpps(x) * rcpps(x));
	__asm
	{
		KFLOATINITDSS( dst, src0, src1, count, pre, post )
		and		eax,15
		jne		lpNA
		jmp		lpA
		align	16
lpA:
		movaps	xmm2,[esi+ebx]
		movaps	xmm3,[esi+ebx+16]
		rcpps	xmm4,xmm2
		rcpps	xmm5,xmm3
		prefetchnta	[esi+ebx+64]
		mulps	xmm2,xmm4
		mulps	xmm2,xmm4
		mulps	xmm3,xmm5
		mulps	xmm3,xmm5
		addps	xmm4,xmm4
		addps	xmm5,xmm5
		subps	xmm4,xmm2
		subps	xmm5,xmm3
		mulps	xmm4,[edx+ebx]
		mulps	xmm5,[edx+ebx+16]
		movaps	[edi+ebx],xmm4
		movaps	[edi+ebx+16],xmm5
		add		ebx,16*2
		jl		lpA
		jmp		done
		align	16
lpNA:
		movups	xmm2,[esi+ebx]
		movups	xmm3,[esi+ebx+16]
		rcpps	xmm4,xmm2
		rcpps	xmm5,xmm3
		prefetchnta	[esi+ebx+64]
		mulps	xmm2,xmm4
		mulps	xmm2,xmm4
		mulps	xmm3,xmm5
		mulps	xmm3,xmm5
		addps	xmm4,xmm4
		addps	xmm5,xmm5
		subps	xmm4,xmm2
		subps	xmm5,xmm3
		movups	xmm2,[edx+ebx]
		movups	xmm3,[edx+ebx+16]
		mulps	xmm4,xmm2
		mulps	xmm5,xmm3
		movaps	[edi+ebx],xmm4
		movaps	[edi+ebx+16],xmm5
		add		ebx,16*2
		jl		lpNA
done:
		mov		edx,src0
		mov		esi,src1
		mov		edi,dst
		KFLOATOPER( KDIVDSS1( [edi+ebx],[edx+ebx],[esi+ebx] ),
					KDIVDSS4( [edi+ebx],[edx+ebx],[esi+ebx] ), count )
	}
}
/*
============
Simd_MulAdd

 assumes count >= 7
============
*/
static void Simd_MulAdd( float *dst, const float constant, const float *src, const int count ) {
	__asm	mov			esi, dst
	__asm	mov			edi, src
	__asm	mov			eax, count
	__asm	shl			eax, 2
	__asm	mov			ecx, esi
	__asm	mov			edx, eax
	__asm	or			ecx, edi
	__asm	fld			constant
	__asm	and			ecx, 15
	__asm	jz			SimdMulAdd16
	__asm	and			ecx, 3
	__asm	jnz			SimdMulAdd8
	__asm	mov			ecx, esi
	__asm	xor			ecx, edi
	__asm	and			ecx, 15
	__asm	jnz			MulAdd8
	__asm	mov			ecx, esi
	__asm	and			ecx, 15
	__asm	neg			ecx
	__asm	add			ecx, 16
	__asm	sub			eax, ecx
	__asm	add			edi, ecx
	__asm	add			esi, ecx
	__asm	neg			ecx
	__asm	mov			edx, eax
	__asm loopPreMulAdd16:
	__asm	fld			st
	__asm	fmul		dword ptr [edi+ecx]
	__asm	fadd		dword ptr [esi+ecx]
	__asm	fstp		dword ptr [esi+ecx]
	__asm	add			ecx, 4
	__asm	jl			loopPreMulAdd16
	__asm SimdMulAdd16:
	__asm	and			eax, ~15
	__asm	movss		xmm1, constant
	__asm	shufps		xmm1, xmm1, 0x00
	__asm	add			esi, eax
	__asm	add			edi, eax
	__asm	neg			eax
	__asm	align		16
	__asm loopMulAdd16:
	__asm	movaps		xmm0, [edi+eax]
	__asm	mulps		xmm0, xmm1
	__asm	addps		xmm0, [esi+eax]
	__asm	movaps		[esi+eax], xmm0
	__asm	add			eax, 16
	__asm	jl			loopMulAdd16
	__asm	jmp			postMulAdd
	__asm MulAdd8:
	__asm	mov			ecx, esi
	__asm	and			ecx, 7
	__asm	jz			SimdMulAdd8
	__asm	sub			eax, ecx
	__asm	add			esi, ecx
	__asm	add			edi, ecx
	__asm	neg			ecx
	__asm	mov			edx, eax
	__asm loopPreMulAdd8:
	__asm	fld			st
	__asm	fmul		dword ptr [edi+ecx]
	__asm	fadd		dword ptr [esi+ecx]
	__asm	fstp		dword ptr [esi+ecx]
	__asm	add			ecx, 4
	__asm	jl			loopPreMulAdd8
	__asm SimdMulAdd8:
	__asm	and			eax, ~15
	__asm	movss		xmm1, constant
	__asm	shufps		xmm1, xmm1, 0x00
	__asm	add			esi, eax
	__asm	add			edi, eax
	__asm	neg			eax
	__asm	align		16
	__asm loopMulAdd8:
	__asm	movlps		xmm0, [edi+eax]
	__asm	movhps		xmm0, [edi+eax+8]
	__asm	mulps		xmm0, xmm1
	__asm	movlps		xmm2, [esi+eax]
	__asm	movhps		xmm2, [esi+eax+8]
	__asm	addps		xmm0, xmm2
	__asm	movlps		[esi+eax], xmm0
	__asm	movhps		[esi+eax+8], xmm0
	__asm	add			eax, 16
	__asm	jl			loopMulAdd8
	__asm	jmp			postMulAdd
	__asm postMulAdd:
	__asm	and			edx, 15
	__asm	jz			MulAddDone
	__asm	add			esi, edx
	__asm	add			edi, edx
	__asm	neg			edx
	__asm loopPostMulAdd:
	__asm	fld			st
	__asm	fmul		dword ptr [edi+edx]
	__asm	fadd		dword ptr [esi+edx]
	__asm	fstp		dword ptr [esi+edx]
	__asm	add			edx, 4
	__asm	jl			loopPostMulAdd
	__asm MulAddDone:
	__asm	fstp		st
}

#define MULADD_FEW( OPER )																				\
switch( count ) {																						\
	case 0:																								\
		return;																							\
	case 1:																								\
		dst[0] OPER c * src[0];																			\
		return;																							\
	case 2:																								\
		dst[0] OPER c * src[0]; dst[1] OPER c * src[1];													\
		return;																							\
	case 3:																								\
		dst[0] OPER c * src[0]; dst[1] OPER c * src[1]; dst[2] OPER c * src[2];							\
		return;																							\
	case 4:																								\
		dst[0] OPER c * src[0]; dst[1] OPER c * src[1]; dst[2] OPER c * src[2]; dst[3] OPER c * src[3];	\
		return;																							\
	case 5:																								\
		dst[0] OPER c * src[0]; dst[1] OPER c * src[1]; dst[2] OPER c * src[2]; dst[3] OPER c * src[3];	\
		dst[4] OPER c * src[4];																			\
		return;																							\
	case 6:																								\
		dst[0] OPER c * src[0]; dst[1] OPER c * src[1]; dst[2] OPER c * src[2]; dst[3] OPER c * src[3];	\
		dst[4] OPER c * src[4]; dst[5] OPER c * src[5];													\
		return;																							\
	case 7:																								\
		dst[0] OPER c * src[0]; dst[1] OPER c * src[1]; dst[2] OPER c * src[2]; dst[3] OPER c * src[3];	\
		dst[4] OPER c * src[4]; dst[5] OPER c * src[5]; dst[6] OPER c * src[6];							\
		return;																							\
	case 8:																								\
		dst[0] OPER c * src[0]; dst[1] OPER c * src[1]; dst[2] OPER c * src[2]; dst[3] OPER c * src[3];	\
		dst[4] OPER c * src[4]; dst[5] OPER c * src[5]; dst[6] OPER c * src[6]; dst[7] OPER c * src[7];	\
		return;																							\
	case 9:																								\
		dst[0] OPER c * src[0]; dst[1] OPER c * src[1]; dst[2] OPER c * src[2]; dst[3] OPER c * src[3];	\
		dst[4] OPER c * src[4]; dst[5] OPER c * src[5]; dst[6] OPER c * src[6]; dst[7] OPER c * src[7];	\
		dst[8] OPER c * src[8];																			\
		return;																							\
	case 10:																							\
		dst[0] OPER c * src[0]; dst[1] OPER c * src[1]; dst[2] OPER c * src[2]; dst[3] OPER c * src[3];	\
		dst[4] OPER c * src[4]; dst[5] OPER c * src[5]; dst[6] OPER c * src[6]; dst[7] OPER c * src[7];	\
		dst[8] OPER c * src[8]; dst[9] OPER c * src[9];													\
		return;																							\
	case 11:																							\
		dst[0] OPER c * src[0]; dst[1] OPER c * src[1]; dst[2] OPER c * src[2]; dst[3] OPER c * src[3];	\
		dst[4] OPER c * src[4]; dst[5] OPER c * src[5]; dst[6] OPER c * src[6]; dst[7] OPER c * src[7];	\
		dst[8] OPER c * src[8]; dst[9] OPER c * src[9]; dst[10] OPER c * src[10];						\
		return;																							\
}

/*
============
CSIMD_SSE::mulAdd

  dst[i] += constant * src[i];
============
*/
void VPCALL CSIMD_SSE::mulAdd( float *dst, const float constant, const float *src, const int count ) {
	float c = constant;
	MULADD_FEW( += )
	Simd_MulAdd( dst, constant, src, count );
}

/*
============
CSIMD_SSE::mulAdd

  dst[i] += src0[i] * src1[i];
============
*/
void VPCALL CSIMD_SSE::mulAdd( float *dst, const float *src0, const float *src1, const int count ) {
	for ( int i = 0; i < count; i++ ) {
		dst[i] += src0[i] + src1[i];
	}
}

/*
============
CSIMD_SSE::mulSub

  dst[i] -= constant * src[i];
============
*/
void VPCALL CSIMD_SSE::mulSub( float *dst, const float constant, const float *src, const int count ) {
	float c = constant;
	MULADD_FEW( -= )
	Simd_MulAdd( dst, -constant, src, count );
}

/*
============
CSIMD_SSE::mulSub

  dst[i] -= src0[i] * src1[i];
============
*/
void VPCALL CSIMD_SSE::mulSub( float *dst, const float *src0, const float *src1, const int count ) {
	for ( int i = 0; i < count; i++ ) {
		dst[i] -= src0[i] + src1[i];
	}
}

/*
============
CSIMD_SSE::dot

  dst[i] = constant * src[i];
============
*/
void VPCALL CSIMD_SSE::dot( float *dst, const CVec3D &constant, const CVec3D *src, const int count ) {
	__asm
	{
		mov			eax, count
		mov			edi, constant
		mov			edx, eax
		mov			esi, src
		mov			ecx, dst
		and			eax, ~3

		movss		xmm4, [edi+0]                          //carrega matriz constante
		shufps		xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )
		movss		xmm5, [edi+4]
		shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )
		movss		xmm6, [edi+8]
		shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )

		jz			done4
		imul		eax, 12
		add			esi, eax
		neg			eax

	loop4:
		movlps		xmm1, [esi+eax+ 0]
		movlps		xmm2, [esi+eax+ 8]
		movlps		xmm3, [esi+eax+16]
		movhps		xmm1, [esi+eax+24]
		movhps		xmm2, [esi+eax+32]
		movhps		xmm3, [esi+eax+40]
		movaps		xmm0, xmm1
		shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 2, 1, 3 )
		shufps		xmm1, xmm3, R_SHUFFLEPS( 1, 3, 0, 2 )
		shufps		xmm2, xmm3, R_SHUFFLEPS( 0, 2, 1, 3 )
		add			ecx, 16
		add			eax, 4*12
		mulps		xmm0, xmm4
		mulps		xmm1, xmm5
		mulps		xmm2, xmm6
		addps		xmm0, xmm1
		addps		xmm0, xmm2
		shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 2, 1, 3 )
		movlps		[ecx-16+0], xmm0
		movhps		[ecx-16+8], xmm0
		jl			loop4

	done4:
		and			edx, 3
		jz			done1

	loop1:
		movss		xmm0, [esi+eax+0]
		movss		xmm1, [esi+eax+4]
		movss		xmm2, [esi+eax+8]
		mulss		xmm0, xmm4
		mulss		xmm1, xmm5
		mulss		xmm2, xmm6
		add			ecx, 4
		addss		xmm0, xmm1
		add			eax, 12
		addss		xmm0, xmm2
		dec			edx
		movss		[ecx-4], xmm0
		jnz			loop1

	done1:
	}
}

/*
============
CSIMD_SSE::dot

  dst[i] = constant * src[i].Normal() + src[i][3];
============
*/
void VPCALL CSIMD_SSE::dot( float *dst, const CVec3D &constant, const CPlane *src, const int count ) {
	__asm {
		mov			eax, count
		mov			edi, constant
		mov			edx, eax
		mov			esi, src
		mov			ecx, dst
		and			eax, ~3

		movss		xmm5, [edi+0]
		shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )
		movss		xmm6, [edi+4]
		shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
		movss		xmm7, [edi+8]
		shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )

		jz			startVert1
		imul		eax, 16
		add			esi, eax
		neg			eax

	loopVert4:

		movlps		xmm1, [esi+eax+ 0]
		movlps		xmm3, [esi+eax+ 8]
		movhps		xmm1, [esi+eax+16]
		movhps		xmm3, [esi+eax+24]
		movlps		xmm2, [esi+eax+32]
		movlps		xmm4, [esi+eax+40]
		movhps		xmm2, [esi+eax+48]
		movhps		xmm4, [esi+eax+56]
		movaps		xmm0, xmm1
		shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 2, 0, 2 )
		shufps		xmm1, xmm2, R_SHUFFLEPS( 1, 3, 1, 3 )
		movaps		xmm2, xmm3
		shufps		xmm2, xmm4, R_SHUFFLEPS( 0, 2, 0, 2 )
		shufps		xmm3, xmm4, R_SHUFFLEPS( 1, 3, 1, 3 )

		add			ecx, 16
		add			eax, 4*16

		mulps		xmm0, xmm5
		mulps		xmm1, xmm6
		mulps		xmm2, xmm7
		addps		xmm0, xmm3
		addps		xmm0, xmm1
		addps		xmm0, xmm2

		movlps		[ecx-16+0], xmm0
		movhps		[ecx-16+8], xmm0
		jl			loopVert4

	startVert1:
		and			edx, 3
		jz			done

	loopVert1:
		movss		xmm0, [esi+eax+0]
		movss		xmm1, [esi+eax+4]
		movss		xmm2, [esi+eax+8]
		mulss		xmm0, xmm5
		mulss		xmm1, xmm6
		mulss		xmm2, xmm7
		addss		xmm0, [esi+eax+12]
		add			ecx, 4
		addss		xmm0, xmm1
		add			eax, 16
		addss		xmm0, xmm2
		dec			edx
		movss		[ecx-4], xmm0
		jnz			loopVert1

	done:
	}
}

/*
============
CSIMD_SSE::dot

  dst[i] = constant * src[i].xyz;
============
*/
void VPCALL CSIMD_SSE::dot( float *dst, const CVec3D &constant, const CVertex *src, const int count ) {

	SMF_ASSERT( sizeof( CVertex ) == DRAWVERT_SIZE );
	SMF_ASSERT( (int)&((CVertex *)0)->xyz == DRAWVERT_XYZ_OFFSET );

	// 0,  1,  2
	// 3,  4,  5
	// 6,  7,  8
	// 9, 10, 11

	__asm {
		mov			eax, count
		mov			edi, constant
		mov			edx, eax
		mov			esi, src
		mov			ecx, dst
		and			eax, ~3

		movss		xmm4, [edi+0]
		shufps		xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )
		movss		xmm5, [edi+4]
		shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )
		movss		xmm6, [edi+8]
		shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )

		jz			startVert1
		imul		eax, DRAWVERT_SIZE
		add			esi, eax
		neg			eax

	loopVert4:
		movss		xmm0, [esi+eax+1*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]	//  3,  X,  X,  X
		movss		xmm2, [esi+eax+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+8]	//  2,  X,  X,  X
		movhps		xmm0, [esi+eax+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]	//  3,  X,  0,  1
		movaps		xmm1, xmm0												//  3,  X,  0,  1

		movlps		xmm1, [esi+eax+1*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+4]	//  4,  5,  0,  1
		shufps		xmm2, xmm1, R_SHUFFLEPS( 0, 1, 0, 1 )					//  2,  X,  4,  5

		movss		xmm3, [esi+eax+3*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]	//  9,  X,  X,  X
		movhps		xmm3, [esi+eax+2*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]	//  9,  X,  6,  7
		shufps		xmm0, xmm3, R_SHUFFLEPS( 2, 0, 2, 0 )					//  0,  3,  6,  9

		movlps		xmm3, [esi+eax+3*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+4]	// 10, 11,  6,  7
		shufps		xmm1, xmm3, R_SHUFFLEPS( 3, 0, 3, 0 )					//  1,  4,  7, 10

		movhps		xmm3, [esi+eax+2*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+8]	// 10, 11,  8,  X
		shufps		xmm2, xmm3, R_SHUFFLEPS( 0, 3, 2, 1 )					//  2,  5,  8, 11

		add			ecx, 16
		add			eax, 4*DRAWVERT_SIZE

		mulps		xmm0, xmm4
		mulps		xmm1, xmm5
		mulps		xmm2, xmm6
		addps		xmm0, xmm1
		addps		xmm0, xmm2

		movlps		[ecx-16+0], xmm0
		movhps		[ecx-16+8], xmm0
		jl			loopVert4

	startVert1:
		and			edx, 3
		jz			done

	loopVert1:
		movss		xmm0, [esi+eax+DRAWVERT_XYZ_OFFSET+0]
		movss		xmm1, [esi+eax+DRAWVERT_XYZ_OFFSET+4]
		movss		xmm2, [esi+eax+DRAWVERT_XYZ_OFFSET+8]
		mulss		xmm0, xmm4
		mulss		xmm1, xmm5
		mulss		xmm2, xmm6
		add			ecx, 4
		addss		xmm0, xmm1
		add			eax, DRAWVERT_SIZE
		addss		xmm0, xmm2
		dec			edx
		movss		[ecx-4], xmm0
		jnz			loopVert1

	done:
	}
}

/*
============
CSIMD_SSE::dot

  dst[i] = constant.Normal() * src[i] + constant[3];
============
*/
void VPCALL CSIMD_SSE::dot( float *dst, const CPlane &constant, const CVec3D *src, const int count ) {
	__asm
	{
		mov			eax, count
		mov			edi, constant
		mov			edx, eax
		mov			esi, src
		mov			ecx, dst
		and			eax, ~3

		movss		xmm4, [edi+0]
		shufps		xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )
		movss		xmm5, [edi+4]
		shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )
		movss		xmm6, [edi+8]
		shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
		movss		xmm7, [edi+12]
		shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )

		jz			done4
		imul		eax, 12
		add			esi, eax
		neg			eax

	loop4:
		movlps		xmm1, [esi+eax+ 0]
		movlps		xmm2, [esi+eax+ 8]
		movlps		xmm3, [esi+eax+16]
		movhps		xmm1, [esi+eax+24]
		movhps		xmm2, [esi+eax+32]
		movhps		xmm3, [esi+eax+40]
		movaps		xmm0, xmm1
		shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 2, 1, 3 )
		shufps		xmm1, xmm3, R_SHUFFLEPS( 1, 3, 0, 2 )
		shufps		xmm2, xmm3, R_SHUFFLEPS( 0, 2, 1, 3 )

		add			ecx, 16
		add			eax, 4*12

		mulps		xmm0, xmm4
		mulps		xmm1, xmm5
		mulps		xmm2, xmm6
		addps		xmm0, xmm7
		addps		xmm0, xmm1
		addps		xmm0, xmm2
		shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 2, 1, 3 )

		movlps		[ecx-16+0], xmm0
		movhps		[ecx-16+8], xmm0
		jl			loop4

	done4:
		and			edx, 3
		jz			done1

	loop1:
		movss		xmm0, [esi+eax+0]
		movss		xmm1, [esi+eax+4]
		movss		xmm2, [esi+eax+8]
		mulss		xmm0, xmm4
		mulss		xmm1, xmm5
		mulss		xmm2, xmm6
		addss		xmm0, xmm7
		add			ecx, 4
		addss		xmm0, xmm1
		add			eax, 12
		addss		xmm0, xmm2
		dec			edx
		movss		[ecx-4], xmm0
		jnz			loop1

	done1:
	}
}

/*
============
CSIMD_SSE::dot

  dst[i] = constant.Normal() * src[i].Normal() + constant[3] * src[i][3];
============
*/
void VPCALL CSIMD_SSE::dot( float *dst, const CPlane &constant, const CPlane *src, const int count ) {

#define SINGLE_OP(SRC, DEST)							\
	__asm	movlps		xmm0,[SRC]						\
	__asm	movlps		xmm1,[SRC+8]					\
	__asm	mulps		xmm0,xmm4						\
	__asm	mulps		xmm1,xmm5						\
	__asm	addps		xmm0,xmm1						\
	__asm	movaps		xmm1,xmm0						\
	__asm	shufps		xmm1,xmm1,SHUFFLEPS(1,1,1,1)	\
	__asm	addss		xmm0,xmm1						\
	__asm	movss		[DEST],xmm0						\
	__asm	add			SRC,16							\
	__asm	add			DEST,4

#define DUAL_OP(SRC, DEST)								\
	__asm	movlps		xmm0,[SRC]						\
	__asm	movlps		xmm1,[SRC+8]					\
	__asm	movhps		xmm0,[SRC+16]					\
	__asm	movhps		xmm1,[SRC+24]					\
	__asm	mulps		xmm0,xmm4						\
	__asm	mulps		xmm1,xmm5						\
	__asm	addps		xmm0,xmm1						\
	__asm	shufps		xmm1,xmm0,SHUFFLEPS(2,0,1,0)	\
	__asm	shufps		xmm0,xmm0,SHUFFLEPS(3,1,2,0)	\
	__asm	addps		xmm0,xmm1						\
	__asm	movhps		[DEST],xmm0						\
	__asm	add			SRC,32							\
	__asm	add			DEST,8

	__asm {
		mov			edx, dst
		mov			eax, src
		mov			ebx, constant
		mov			ecx, count

		movlps		xmm4, [ebx]
		shufps		xmm4, xmm4, SHUFFLEPS(1,0,1,0)
		movlps		xmm5, [ebx+8]
		shufps		xmm5, xmm5, SHUFFLEPS(1,0,1,0)

		xorps		xmm0, xmm0
		xorps		xmm1, xmm1

	_lpAlignDest:
		test		edx, 0x0f
		jz			_destAligned
		SINGLE_OP(eax,edx)
		dec			ecx
		jnz			_lpAlignDest
		jmp			_vpExit

	_destAligned:
		push		ecx

		cmp			ecx, 4
		jl			_post

		and			ecx, ~3
		shl			ecx, 2
		lea			eax, [eax+ecx*4]
		add			edx, ecx
		neg			ecx

		movlps		xmm0, [eax+ecx*4]
		movhps		xmm0, [eax+ecx*4+16]
		movlps		xmm2, [eax+ecx*4+32]
		movhps		xmm2, [eax+ecx*4+48]
		jmp			_lpStart

		align	16
	_lp:
		prefetchnta	[eax+ecx*4+128]
		addps		xmm1, xmm0
		movlps		xmm0, [eax+ecx*4]
		movhps		xmm0, [eax+ecx*4+16]
		movlps		xmm2, [eax+ecx*4+32]
		movhps		xmm2, [eax+ecx*4+48]
		movaps		[edx+ecx-16],xmm1
	_lpStart:
		movlps		xmm1, [eax+ecx*4+8]
		movhps		xmm1, [eax+ecx*4+24]
		movlps		xmm3, [eax+ecx*4+40]
		movhps		xmm3, [eax+ecx*4+56]
		add			ecx, 16
		mulps		xmm1, xmm5
		mulps		xmm2, xmm4
		mulps		xmm3, xmm5
		addps		xmm2, xmm3						// y3+w3 x3+z3 y2+w2 x2+z2
		mulps		xmm0, xmm4
		addps		xmm0, xmm1						// y1+w1 x1+z1 y0+w0 x0+z0
		movaps		xmm1, xmm0
		shufps		xmm0, xmm2, SHUFFLEPS(2,0,2,0)	// x3+z3 x2+z2 x1+z1 x0+z0
		shufps		xmm1, xmm2, SHUFFLEPS(3,1,3,1)	// y3+w3 y2+w2 y1+w1 y0+w0
		js			_lp
		addps		xmm1, xmm0
		movaps		[edx+ecx-16], xmm1
	_post:
		pop			ecx
		and			ecx, 0x3
		cmp			ecx, 2
		jl			_post1
		DUAL_OP(eax,edx)
		sub			ecx, 2
	_post1:
		cmp			ecx, 1
		jne			_vpExit
		SINGLE_OP(eax,edx)
	_vpExit:
	}

#undef DUAL_OP
#undef SINGLE_OP

}

/*
============
CSIMD_SSE::dot

  dst[i] = constant.Normal() * src[i].xyz + constant[3];
============
*/
void VPCALL CSIMD_SSE::dot( float *dst, const CPlane &constant, const CVertex *src, const int count ) {

	SMF_ASSERT( sizeof( CVertex ) == DRAWVERT_SIZE );
	SMF_ASSERT( (int)&((CVertex *)0)->xyz == DRAWVERT_XYZ_OFFSET );

	// 0,  1,  2
	// 3,  4,  5
	// 6,  7,  8
	// 9, 10, 11

	__asm {
		mov			eax, count
		mov			edi, constant
		mov			edx, eax
		mov			esi, src
		mov			ecx, dst
		and			eax, ~3

		movss		xmm4, [edi+0]
		shufps		xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )
		movss		xmm5, [edi+4]
		shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )
		movss		xmm6, [edi+8]
		shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
		movss		xmm7, [edi+12]
		shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )

		jz			startVert1
		imul		eax, DRAWVERT_SIZE
		add			esi, eax
		neg			eax

	loopVert4:
		movss		xmm0, [esi+eax+1*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]	//  3,  X,  X,  X
		movss		xmm2, [esi+eax+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+8]	//  2,  X,  X,  X
		movhps		xmm0, [esi+eax+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]	//  3,  X,  0,  1
		movaps		xmm1, xmm0												//  3,  X,  0,  1

		movlps		xmm1, [esi+eax+1*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+4]	//  4,  5,  0,  1
		shufps		xmm2, xmm1, R_SHUFFLEPS( 0, 1, 0, 1 )					//  2,  X,  4,  5

		movss		xmm3, [esi+eax+3*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]	//  9,  X,  X,  X
		movhps		xmm3, [esi+eax+2*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]	//  9,  X,  6,  7
		shufps		xmm0, xmm3, R_SHUFFLEPS( 2, 0, 2, 0 )					//  0,  3,  6,  9

		movlps		xmm3, [esi+eax+3*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+4]	// 10, 11,  6,  7
		shufps		xmm1, xmm3, R_SHUFFLEPS( 3, 0, 3, 0 )					//  1,  4,  7, 10

		movhps		xmm3, [esi+eax+2*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+8]	// 10, 11,  8,  X
		shufps		xmm2, xmm3, R_SHUFFLEPS( 0, 3, 2, 1 )					//  2,  5,  8, 11

		add			ecx, 16
		add			eax, 4*DRAWVERT_SIZE

		mulps		xmm0, xmm4
		mulps		xmm1, xmm5
		mulps		xmm2, xmm6
		addps		xmm0, xmm7
		addps		xmm0, xmm1
		addps		xmm0, xmm2

		movlps		[ecx-16+0], xmm0
		movhps		[ecx-16+8], xmm0
		jl			loopVert4

	startVert1:
		and			edx, 3
		jz			done

	loopVert1:
		movss		xmm0, [esi+eax+DRAWVERT_XYZ_OFFSET+0]
		movss		xmm1, [esi+eax+DRAWVERT_XYZ_OFFSET+4]
		movss		xmm2, [esi+eax+DRAWVERT_XYZ_OFFSET+8]
		mulss		xmm0, xmm4
		mulss		xmm1, xmm5
		mulss		xmm2, xmm6
		addss		xmm0, xmm7
		add			ecx, 4
		addss		xmm0, xmm1
		add			eax, DRAWVERT_SIZE
		addss		xmm0, xmm2
		dec			edx
		movss		[ecx-4], xmm0
		jnz			loopVert1

	done:
	}
}

/*
============
CSIMD_SSE::dot

  dst[i] = src0[i] * src1[i];
============
*/
void VPCALL CSIMD_SSE::dot( float *dst, const CVec3D *src0, const CVec3D *src1, const int count ) {
	__asm
	{
		mov			eax, count
		mov			edi, src0
		mov			edx, eax
		mov			esi, src1
		mov			ecx, dst
		and			eax, ~3

		jz			done4
		imul		eax, 12
		add			edi, eax
		add			esi, eax
		neg			eax

	loop4:
		movlps		xmm0, [esi+eax]						// 0, 1, X, X
		movlps		xmm3, [edi+eax]						// 0, 1, X, X
		movlps		xmm1, [esi+eax+8]					// 2, 3, X, X
		movlps		xmm4, [edi+eax+8]					// 2, 3, X, X
		movhps		xmm0, [esi+eax+24]					// 0, 1, 6, 7
		movhps		xmm3, [edi+eax+24]					// 0, 1, 6, 7
		movhps		xmm1, [esi+eax+32]					// 2, 3, 8, 9
		movhps		xmm4, [edi+eax+32]					// 2, 3, 8, 9
		movlps		xmm2, [esi+eax+16]					// 4, 5, X, X
		movlps		xmm5, [edi+eax+16]					// 4, 5, X, X
		movhps		xmm2, [esi+eax+40]					// 4, 5, 10, 11
		movhps		xmm5, [edi+eax+40]					// 4, 5, 10, 11

		add			ecx, 16
		add			eax, 48

		mulps		xmm0, xmm3
		mulps		xmm1, xmm4
		mulps		xmm2, xmm5
		movaps		xmm7, xmm0
		shufps		xmm7, xmm1, R_SHUFFLEPS( 0, 2, 1, 3 )	// 0, 6, 3, 9
		shufps		xmm0, xmm2, R_SHUFFLEPS( 1, 3, 0, 2 )	// 1, 7, 4, 10
		shufps		xmm1, xmm2, R_SHUFFLEPS( 0, 2, 1, 3 )	// 2, 8, 5, 11
		addps		xmm7, xmm0
		addps		xmm7, xmm1
		shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 2, 1, 3 )

		movlps		[ecx-16+0], xmm7
		movhps		[ecx-16+8], xmm7
		jl			loop4

	done4:
		and			edx, 3
		jz			done1

	loop1:
		movss		xmm0, [esi+eax+0]
		movss		xmm3, [edi+eax+0]
		movss		xmm1, [esi+eax+4]
		movss		xmm4, [edi+eax+4]
		movss		xmm2, [esi+eax+8]
		movss		xmm5, [edi+eax+8]
		mulss		xmm0, xmm3
		mulss		xmm1, xmm4
		mulss		xmm2, xmm5
		add			ecx, 4
		addss		xmm0, xmm1
		add			eax, 12
		addss		xmm0, xmm2
		dec			edx
		movss		[ecx-4], xmm0
		jnz			loop1

	done1:
	}
}

/*
============
CSIMD_SSE::dot

  dot = src1[0] * src2[0] + src1[1] * src2[1] + src1[2] * src2[2] + ...
============
*/
void VPCALL CSIMD_SSE::dot( float &dot, const float *src1, const float *src2, const int count ) {
	switch( count ) {
		case 0:
			dot = 0.0f;
			return;
		case 1:
			dot = src1[0] * src2[0];
			return;
		case 2:
			dot = src1[0] * src2[0] + src1[1] * src2[1];
			return;
		case 3:
			dot = src1[0] * src2[0] + src1[1] * src2[1] + src1[2] * src2[2];
			return;
		default:
			__asm {
				mov			ecx, src1
				mov			edx, src2
				mov			eax, ecx
				or			eax, edx
				and			eax, 15
				jz			alignedDot
				// unaligned
				mov			eax, count
				shr			eax, 2
				shl			eax, 4
				add			ecx, eax
				add			edx, eax
				neg			eax
				movups		xmm0, [ecx+eax]
				movups		xmm1, [edx+eax]
				mulps		xmm0, xmm1
				add			eax, 16
				jz			doneDot
			loopUnalignedDot:
				movups		xmm1, [ecx+eax]
				movups		xmm2, [edx+eax]
				mulps		xmm1, xmm2
				addps		xmm0, xmm1
				add			eax, 16
				jl			loopUnalignedDot
				jmp			doneDot
				// aligned
			alignedDot:
				mov			eax, count
				shr			eax, 2
				shl			eax, 4
				add			ecx, eax
				add			edx, eax
				neg			eax
				movaps		xmm0, [ecx+eax]
				movaps		xmm1, [edx+eax]
				mulps		xmm0, xmm1
				add			eax, 16
				jz			doneDot
			loopAlignedDot:
				movaps		xmm1, [ecx+eax]
				movaps		xmm2, [edx+eax]
				mulps		xmm1, xmm2
				addps		xmm0, xmm1
				add			eax, 16
				jl			loopAlignedDot
			doneDot:
			}
			switch( count & 3 ) {
				case 1:
					__asm {
						movss	xmm1, [ecx]
						movss	xmm2, [edx]
						mulss	xmm1, xmm2
						addss	xmm0, xmm1
					}
					break;
				case 2:
					__asm {
						xorps	xmm2, xmm2
						movlps	xmm1, [ecx]
						movlps	xmm2, [edx]
						mulps	xmm1, xmm2
						addps	xmm0, xmm1
					}
					break;
				case 3:
					__asm {
						movss	xmm1, [ecx]
						movhps	xmm1, [ecx+4]
						movss	xmm2, [edx]
						movhps	xmm2, [edx+4]
						mulps	xmm1, xmm2
						addps	xmm0, xmm1
					}
					break;
			}
			__asm {
				movhlps		xmm1, xmm0
				addps		xmm0, xmm1
				movaps		xmm1, xmm0
				shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 0, 0, 0 )
				addss		xmm0, xmm1
				mov			eax, dot
				movss		[eax], xmm0
			}
			return;
	}
}

//
//	cmpeqps		==		Equal
//	cmpneqps	!=		Not Equal
//	cmpltps		<		Less Than
//  cmpnltps	>=		Not Less Than
//	cmpnleps	>		Not Less Or Equal
//
#define FLIP	not al
#define NOFLIP

#define COMPARECONSTANT( DST, SRC0, CONSTANT, COUNT, CMP, CMPSIMD, DOFLIP )				\
	int i, cnt, pre, post;																\
	float *aligned;																		\
																						\
	/* if the float array is not aligned on a 4 sf_u8 boundary */						\
	if ( ((int) SRC0) & 3 ) {															\
		/* unaligned memory access */													\
		pre = 0;																		\
		cnt = COUNT >> 2;																\
		post = COUNT - (cnt<<2);														\
		__asm	mov			edx, cnt													\
		__asm	test		edx, edx													\
		__asm	je			doneCmp														\
		__asm	push		ebx															\
		__asm	neg			edx															\
		__asm	mov			esi, SRC0													\
		__asm	prefetchnta	[esi+64]													\
		__asm	movss		xmm1, CONSTANT												\
		__asm	shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )						\
		__asm	mov			edi, DST													\
		__asm	mov			ecx, 0x01010101												\
		__asm loopNA:																	\
		__asm	movups		xmm0, [esi]													\
		__asm	prefetchnta	[esi+128]													\
		__asm	CMPSIMD		xmm0, xmm1													\
		__asm	movmskps	eax, xmm0													\
		__asm	DOFLIP																	\
		__asm	mov			ah, al														\
		__asm	shr			ah, 1														\
		__asm	mov			bx, ax														\
		__asm	shl			ebx, 14														\
		__asm	mov			bx, ax														\
		__asm	and			ebx, ecx													\
		__asm	mov			dword ptr [edi], ebx										\
		__asm	add			esi, 16														\
		__asm	add			edi, 4														\
		__asm	inc			edx															\
		__asm	jl			loopNA														\
		__asm	pop			ebx															\
	}																					\
	else {																				\
		/* aligned memory access */														\
		aligned = (float *) ((((int) SRC0) + 15) & ~15);								\
		if ( (int)aligned > ((int)src0) + COUNT ) {										\
			pre = COUNT;																\
			post = 0;																	\
		}																				\
		else {																			\
			pre = aligned - SRC0;														\
			cnt = (COUNT - pre) >> 2;													\
			post = COUNT - pre - (cnt<<2);												\
			__asm	mov			edx, cnt												\
			__asm	test		edx, edx												\
			__asm	je			doneCmp													\
			__asm	push		ebx														\
			__asm	neg			edx														\
			__asm	mov			esi, aligned											\
			__asm	prefetchnta	[esi+64]												\
			__asm	movss		xmm1, CONSTANT											\
			__asm	shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )					\
			__asm	mov			edi, DST												\
			__asm	add			edi, pre												\
			__asm	mov			ecx, 0x01010101											\
			__asm loopA:																\
			__asm	movaps		xmm0, [esi]												\
			__asm	prefetchnta	[esi+128]												\
			__asm	CMPSIMD		xmm0, xmm1												\
			__asm	movmskps	eax, xmm0												\
			__asm	DOFLIP																\
			__asm	mov			ah, al													\
			__asm	shr			ah, 1													\
			__asm	mov			bx, ax													\
			__asm	shl			ebx, 14													\
			__asm	mov			bx, ax													\
			__asm	and			ebx, ecx												\
			__asm	mov			dword ptr [edi], ebx									\
			__asm	add			esi, 16													\
			__asm	add			edi, 4													\
			__asm	inc			edx														\
			__asm	jl			loopA													\
			__asm	pop			ebx														\
		}																				\
	}																					\
	doneCmp:																			\
	double c = constant;																\
	for ( i = 0; i < pre; i++ ) {														\
		dst[i] = src0[i] CMP c;															\
	}																					\
 	for ( i = count - post; i < count; i++ ) {											\
		dst[i] = src0[i] CMP c;															\
	}

#define COMPAREBITCONSTANT( DST, BITNUM, SRC0, CONSTANT, COUNT, CMP, CMPSIMD, DOFLIP )	\
	int i, cnt, pre, post;																\
	float *aligned;																		\
																						\
	/* if the float array is not aligned on a 4 sf_u8 boundary */						\
	if ( ((int) SRC0) & 3 ) {															\
		/* unaligned memory access */													\
		pre = 0;																		\
		cnt = COUNT >> 2;																\
		post = COUNT - (cnt<<2);														\
		__asm	mov			edx, cnt													\
		__asm	test		edx, edx													\
		__asm	je			doneCmp														\
		__asm	push		ebx															\
		__asm	neg			edx															\
		__asm	mov			esi, SRC0													\
		__asm	prefetchnta	[esi+64]													\
		__asm	movss		xmm1, CONSTANT												\
		__asm	shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )						\
		__asm	mov			edi, DST													\
		__asm	mov			cl, bitNum													\
		__asm loopNA:																	\
		__asm	movups		xmm0, [esi]													\
		__asm	prefetchnta	[esi+128]													\
		__asm	CMPSIMD		xmm0, xmm1													\
		__asm	movmskps	eax, xmm0													\
		__asm	DOFLIP																	\
		__asm	mov			ah, al														\
		__asm	shr			ah, 1														\
		__asm	mov			bx, ax														\
		__asm	shl			ebx, 14														\
		__asm	mov			bx, ax														\
		__asm	and			ebx, 0x01010101												\
		__asm	shl			ebx, cl														\
		__asm	or			ebx, dword ptr [edi]										\
		__asm	mov			dword ptr [edi], ebx										\
		__asm	add			esi, 16														\
		__asm	add			edi, 4														\
		__asm	inc			edx															\
		__asm	jl			loopNA														\
		__asm	pop			ebx															\
	}																					\
	else {																				\
		/* aligned memory access */														\
		aligned = (float *) ((((int) SRC0) + 15) & ~15);								\
		if ( (int)aligned > ((int)src0) + COUNT ) {										\
			pre = COUNT;																\
			post = 0;																	\
		}																				\
		else {																			\
			pre = aligned - SRC0;														\
			cnt = (COUNT - pre) >> 2;													\
			post = COUNT - pre - (cnt<<2);												\
			__asm	mov			edx, cnt												\
			__asm	test		edx, edx												\
			__asm	je			doneCmp													\
			__asm	push		ebx														\
			__asm	neg			edx														\
			__asm	mov			esi, aligned											\
			__asm	prefetchnta	[esi+64]												\
			__asm	movss		xmm1, CONSTANT											\
			__asm	shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )					\
			__asm	mov			edi, DST												\
			__asm	add			edi, pre												\
			__asm	mov			cl, bitNum												\
			__asm loopA:																\
			__asm	movaps		xmm0, [esi]												\
			__asm	prefetchnta	[esi+128]												\
			__asm	CMPSIMD		xmm0, xmm1												\
			__asm	movmskps	eax, xmm0												\
			__asm	DOFLIP																\
			__asm	mov			ah, al													\
			__asm	shr			ah, 1													\
			__asm	mov			bx, ax													\
			__asm	shl			ebx, 14													\
			__asm	mov			bx, ax													\
			__asm	and			ebx, 0x01010101											\
			__asm	shl			ebx, cl													\
			__asm	or			ebx, dword ptr [edi]									\
			__asm	mov			dword ptr [edi], ebx									\
			__asm	add			esi, 16													\
			__asm	add			edi, 4													\
			__asm	inc			edx														\
			__asm	jl			loopA													\
			__asm	pop			ebx														\
		}																				\
	}																					\
	doneCmp:																			\
	float c = constant;																	\
	for ( i = 0; i < pre; i++ ) {														\
		dst[i] |= ( src0[i] CMP c ) << BITNUM;											\
	}																					\
 	for ( i = count - post; i < count; i++ ) {											\
		dst[i] |= ( src0[i] CMP c ) << BITNUM;											\
	}

/*
============
CSIMD_SSE::cmpGT

  dst[i] = src0[i] > constant;
============
*/
void VPCALL CSIMD_SSE::cmpGT( sf_u8 *dst, const float *src0, const float constant, const int count ) {
	COMPARECONSTANT( dst, src0, constant, count, >, cmpnleps, NOFLIP )
}

/*
============
CSIMD_SSE::cmpGT

  dst[i] |= ( src0[i] > constant ) << bitNum;
============
*/
void VPCALL CSIMD_SSE::cmpGT( sf_u8 *dst, const sf_u8 bitNum, const float *src0, const float constant, const int count ) {
	COMPAREBITCONSTANT( dst, bitNum, src0, constant, count, >, cmpnleps, NOFLIP )
}

/*
============
CSIMD_SSE::cmpGE

  dst[i] = src0[i] >= constant;
============
*/
void VPCALL CSIMD_SSE::cmpGE( sf_u8 *dst, const float *src0, const float constant, const int count ) {
	COMPARECONSTANT( dst, src0, constant, count, >=, cmpnltps, NOFLIP )
}

/*
============
CSIMD_SSE::cmpGE

  dst[i] |= ( src0[i] >= constant ) << bitNum;
============
*/
void VPCALL CSIMD_SSE::cmpGE( sf_u8 *dst, const sf_u8 bitNum, const float *src0, const float constant, const int count ) {
	COMPAREBITCONSTANT( dst, bitNum, src0, constant, count, >=, cmpnltps, NOFLIP )
}

/*
============
CSIMD_SSE::cmpLT

  dst[i] = src0[i] < constant;
============
*/
void VPCALL CSIMD_SSE::cmpLT( sf_u8 *dst, const float *src0, const float constant, const int count ) {
	COMPARECONSTANT( dst, src0, constant, count, <, cmpltps, NOFLIP )
}

/*
============
CSIMD_SSE::cmpLT

  dst[i] |= ( src0[i] < constant ) << bitNum;
============
*/
void VPCALL CSIMD_SSE::cmpLT( sf_u8 *dst, const sf_u8 bitNum, const float *src0, const float constant, const int count ) {
	COMPAREBITCONSTANT( dst, bitNum, src0, constant, count, <, cmpltps, NOFLIP )
}

/*
============
CSIMD_SSE::cmpLE

  dst[i] = src0[i] <= constant;
============
*/
void VPCALL CSIMD_SSE::cmpLE( sf_u8 *dst, const float *src0, const float constant, const int count ) {
	COMPARECONSTANT( dst, src0, constant, count, <=, cmpnleps, FLIP )
}

/*
============
CSIMD_SSE::cmpLE

  dst[i] |= ( src0[i] <= constant ) << bitNum;
============
*/
void VPCALL CSIMD_SSE::cmpLE( sf_u8 *dst, const sf_u8 bitNum, const float *src0, const float constant, const int count ) {
	COMPAREBITCONSTANT( dst, bitNum, src0, constant, count, <=, cmpnleps, FLIP )
}

/*
============
CSIMD_SSE::minMax
============
*/
void VPCALL CSIMD_SSE::minMax( float &min, float &max, const float *src, const int count ) {
	int i, pre, post;

	min = CMath::INFINITY_FLOAT; max = -CMath::INFINITY_FLOAT;

	__asm
	{
		push		ebx
		mov			eax, min
		mov			ebx, max
		movss		xmm0, [eax]
		movss		xmm1, [ebx]
		shufps		xmm0, xmm0, 0
		shufps		xmm1, xmm1, 0

		KFLOATINITS( src, count, pre, post )
		and			eax, 15
		jz			lpA
		jmp			lpNA
		align		16
lpNA:
		movups		xmm2, [edx+ebx]
		movups		xmm3, [edx+ebx+16]
		minps		xmm0, xmm2
		maxps		xmm1, xmm2
		prefetchnta	[edx+ebx+64]
		minps		xmm0, xmm3
		maxps		xmm1, xmm3
		add			ebx, 16*2
		jl			lpNA
		jmp			done2
lpA:
		movaps		xmm2, [edx+ebx]
		movaps		xmm3, [edx+ebx+16]
		minps		xmm0, xmm2
		maxps		xmm1, xmm2
		prefetchnta	[edx+ebx+64]
		minps		xmm0, xmm3
		maxps		xmm1, xmm3
		add			ebx, 16*2
		jl			lpA
		jmp			done2
		align		16
done2:
		movaps		xmm2, xmm0
		movaps		xmm3, xmm1
		shufps		xmm2, xmm2, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm3, xmm3, R_SHUFFLEPS( 1, 2, 3, 0 )
		minss		xmm0, xmm2
		maxss		xmm1, xmm3
		shufps		xmm2, xmm2, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm3, xmm3, R_SHUFFLEPS( 1, 2, 3, 0 )
		minss		xmm0, xmm2
		maxss		xmm1, xmm3
		shufps		xmm2, xmm2, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm3, xmm3, R_SHUFFLEPS( 1, 2, 3, 0 )
		minss		xmm0, xmm2
		maxss		xmm1, xmm3
		mov			eax, min
		mov			ebx, max
		movss		[eax], xmm0
		movss		[ebx], xmm1
done:
		pop			ebx
	}

	for ( i = 0; i < pre; i++ ) {
		float tmp = src[i];
		if ( tmp > max ) {
			max = tmp;
		}
		if ( tmp < min ) {
			min = tmp;
		}
	}
 	for ( i = count - post; i < count; i++ ) {
		float tmp = src[i];
		if ( tmp > max ) {
			max = tmp;
		}
		if ( tmp < min ) {
			min = tmp;
		}
	}
}

/*
============
CSIMD_SSE::minMax
============
*/
void VPCALL CSIMD_SSE::minMax( CVec2D &min, CVec2D &max, const CVec2D *src, const int count ) {
	__asm {
		mov			eax, count
		test		eax, eax
		movss		xmm0, CMath::INFINITY_FLOAT
		xorps		xmm1, xmm1
		shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
		subps		xmm1, xmm0
		jz			done
		mov			ecx, eax
		and			ecx, 1
		mov			esi, src
		jz			startLoop
		movlps		xmm2, [esi]
		shufps		xmm2, xmm2, R_SHUFFLEPS( 0, 1, 0, 1 )
		dec			eax
		add			esi, 2*4
		minps		xmm0, xmm2
		maxps		xmm1, xmm2
	startLoop:
		imul		eax, 2*4
		add			esi, eax
		neg			eax
	loopVert:
		movlps		xmm2, [esi+eax]
		movhps		xmm2, [esi+eax+8]
		add			eax, 4*4
		minps		xmm0, xmm2
		maxps		xmm1, xmm2
		jl			loopVert
	done:
		movaps		xmm2, xmm0
		shufps		xmm2, xmm2, R_SHUFFLEPS( 2, 3, 0, 1 )
		minps		xmm0, xmm2
		mov			esi, min
		movlps		[esi], xmm0
		movaps		xmm3, xmm1
		shufps		xmm3, xmm3, R_SHUFFLEPS( 2, 3, 0, 1 )
		maxps		xmm1, xmm3
		mov			edi, max
		movlps		[edi], xmm1
	}
}

/*
============
CSIMD_SSE::minMax
============
*/
void VPCALL CSIMD_SSE::minMax( CVec3D &min, CVec3D &max, const CVec3D *src, const int count ) {
	__asm {

		movss		xmm0, CMath::INFINITY_FLOAT
		xorps		xmm1, xmm1
		shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
		subps		xmm1, xmm0
		movaps		xmm2, xmm0
		movaps		xmm3, xmm1

		mov			esi, src
		mov			eax, count
		and			eax, ~3
		jz			done4
		imul		eax, 12
		add			esi, eax
		neg			eax

	loop4:
//		prefetchnta	[esi+4*12]

		movss		xmm4, [esi+eax+0*12+8]
		movhps		xmm4, [esi+eax+0*12+0]
		minps		xmm0, xmm4
		maxps		xmm1, xmm4

		movss		xmm5, [esi+eax+1*12+0]
		movhps		xmm5, [esi+eax+1*12+4]
		minps		xmm2, xmm5
		maxps		xmm3, xmm5

		movss		xmm6, [esi+eax+2*12+8]
		movhps		xmm6, [esi+eax+2*12+0]
		minps		xmm0, xmm6
		maxps		xmm1, xmm6

		movss		xmm7, [esi+eax+3*12+0]
		movhps		xmm7, [esi+eax+3*12+4]
		minps		xmm2, xmm7
		maxps		xmm3, xmm7

		add			eax, 4*12
		jl			loop4

	done4:
		mov			eax, count
		and			eax, 3
		jz			done1
		imul		eax, 12
		add			esi, eax
		neg			eax

	loop1:
		movss		xmm4, [esi+eax+0*12+8]
		movhps		xmm4, [esi+eax+0*12+0]
		minps		xmm0, xmm4
		maxps		xmm1, xmm4

		add			eax, 12
		jl			loop1

	done1:
		shufps		xmm2, xmm2, R_SHUFFLEPS( 3, 1, 0, 2 )
		shufps		xmm3, xmm3, R_SHUFFLEPS( 3, 1, 0, 2 )
		minps		xmm0, xmm2
		maxps		xmm1, xmm3
		mov			esi, min
		movhps		[esi], xmm0
		movss		[esi+8], xmm0
		mov			edi, max
		movhps		[edi], xmm1
		movss		[edi+8], xmm1
	}
}

/*
============
CSIMD_SSE::minMax
============
*/
void VPCALL CSIMD_SSE::minMax( CVec3D &min, CVec3D &max, const CVertex *src, const int count ) {

	SMF_ASSERT( sizeof( CVertex ) == DRAWVERT_SIZE );
	SMF_ASSERT( (int)&((CVertex *)0)->xyz == DRAWVERT_XYZ_OFFSET );

	__asm {

		movss		xmm0, CMath::INFINITY_FLOAT
		xorps		xmm1, xmm1
		shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
		subps		xmm1, xmm0
		movaps		xmm2, xmm0
		movaps		xmm3, xmm1

		mov			esi, src
		mov			eax, count
		and			eax, ~3
		jz			done4
		imul		eax, DRAWVERT_SIZE
		add			esi, eax
		neg			eax

	loop4:
//		prefetchnta	[esi+4*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET]

		movss		xmm4, [esi+eax+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+8]
		movhps		xmm4, [esi+eax+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]
		minps		xmm0, xmm4
		maxps		xmm1, xmm4

		movss		xmm5, [esi+eax+1*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]
		movhps		xmm5, [esi+eax+1*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+4]
		minps		xmm2, xmm5
		maxps		xmm3, xmm5

		movss		xmm6, [esi+eax+2*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+8]
		movhps		xmm6, [esi+eax+2*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]
		minps		xmm0, xmm6
		maxps		xmm1, xmm6

		movss		xmm7, [esi+eax+3*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]
		movhps		xmm7, [esi+eax+3*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+4]
		minps		xmm2, xmm7
		maxps		xmm3, xmm7

		add			eax, 4*DRAWVERT_SIZE
		jl			loop4

	done4:
		mov			eax, count
		and			eax, 3
		jz			done1
		imul		eax, DRAWVERT_SIZE
		add			esi, eax
		neg			eax

	loop1:
		movss		xmm4, [esi+eax+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+8]
		movhps		xmm4, [esi+eax+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]
		minps		xmm0, xmm4
		maxps		xmm1, xmm4

		add			eax, DRAWVERT_SIZE
		jl			loop1

	done1:
		shufps		xmm2, xmm2, R_SHUFFLEPS( 3, 1, 0, 2 )
		shufps		xmm3, xmm3, R_SHUFFLEPS( 3, 1, 0, 2 )
		minps		xmm0, xmm2
		maxps		xmm1, xmm3
		mov			esi, min
		movhps		[esi], xmm0
		movss		[esi+8], xmm0
		mov			edi, max
		movhps		[edi], xmm1
		movss		[edi+8], xmm1
	}
}

/*
============
CSIMD_SSE::minMax
============
*/
void VPCALL CSIMD_SSE::minMax( CVec3D &min, CVec3D &max, const CVertex *src, const int *indexes, const int count ) {

	SMF_ASSERT( sizeof( CVertex ) == DRAWVERT_SIZE );
	SMF_ASSERT( (int)&((CVertex *)0)->xyz == DRAWVERT_XYZ_OFFSET );

	__asm {

		movss		xmm0, CMath::INFINITY_FLOAT
		xorps		xmm1, xmm1
		shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
		subps		xmm1, xmm0
		movaps		xmm2, xmm0
		movaps		xmm3, xmm1

		mov			edi, indexes
		mov			esi, src
		mov			eax, count
		and			eax, ~3
		jz			done4
		shl			eax, 2
		add			edi, eax
		neg			eax

	loop4:
//		prefetchnta	[edi+128]
//		prefetchnta	[esi+4*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET]

		mov			edx, [edi+eax+0]
		imul		edx, DRAWVERT_SIZE
		movss		xmm4, [esi+edx+DRAWVERT_XYZ_OFFSET+8]
		movhps		xmm4, [esi+edx+DRAWVERT_XYZ_OFFSET+0]
		minps		xmm0, xmm4
		maxps		xmm1, xmm4

		mov			edx, [edi+eax+4]
		imul		edx, DRAWVERT_SIZE
		movss		xmm5, [esi+edx+DRAWVERT_XYZ_OFFSET+0]
		movhps		xmm5, [esi+edx+DRAWVERT_XYZ_OFFSET+4]
		minps		xmm2, xmm5
		maxps		xmm3, xmm5

		mov			edx, [edi+eax+8]
		imul		edx, DRAWVERT_SIZE
		movss		xmm6, [esi+edx+DRAWVERT_XYZ_OFFSET+8]
		movhps		xmm6, [esi+edx+DRAWVERT_XYZ_OFFSET+0]
		minps		xmm0, xmm6
		maxps		xmm1, xmm6

		mov			edx, [edi+eax+12]
		imul		edx, DRAWVERT_SIZE
		movss		xmm7, [esi+edx+DRAWVERT_XYZ_OFFSET+0]
		movhps		xmm7, [esi+edx+DRAWVERT_XYZ_OFFSET+4]
		minps		xmm2, xmm7
		maxps		xmm3, xmm7

		add			eax, 4*4
		jl			loop4

	done4:
		mov			eax, count
		and			eax, 3
		jz			done1
		shl			eax, 2
		add			edi, eax
		neg			eax

	loop1:
		mov			edx, [edi+eax+0]
		imul		edx, DRAWVERT_SIZE;
		movss		xmm4, [esi+edx+DRAWVERT_XYZ_OFFSET+8]
		movhps		xmm4, [esi+edx+DRAWVERT_XYZ_OFFSET+0]
		minps		xmm0, xmm4
		maxps		xmm1, xmm4

		add			eax, 4
		jl			loop1

	done1:
		shufps		xmm2, xmm2, R_SHUFFLEPS( 3, 1, 0, 2 )
		shufps		xmm3, xmm3, R_SHUFFLEPS( 3, 1, 0, 2 )
		minps		xmm0, xmm2
		maxps		xmm1, xmm3
		mov			esi, min
		movhps		[esi], xmm0
		movss		[esi+8], xmm0
		mov			edi, max
		movhps		[edi], xmm1
		movss		[edi+8], xmm1
	}
}

/*
============
CSIMD_SSE::clamp
============
*/
void VPCALL CSIMD_SSE::clamp( float *dst, const float *src, const float min, const float max, const int count ) {
	int	i, pre, post;

	__asm
	{
		movss	xmm0,min
		movss	xmm1,max
		shufps	xmm0,xmm0,0
		shufps	xmm1,xmm1,0

		KFLOATINITDS( dst, src, count, pre, post )
		and		eax,15
		jne		lpNA
		jmp		lpA
		align	16
lpA:
		movaps	xmm2,[edx+ebx]
		movaps	xmm3,[edx+ebx+16]
		maxps	xmm2,xmm0
		maxps	xmm3,xmm0
		prefetchnta	[edx+ebx+64]
		minps	xmm2,xmm1
		minps	xmm3,xmm1
		movaps	[edi+ebx],xmm2
		movaps	[edi+ebx+16],xmm3
		add		ebx,16*2
		jl		lpA
		jmp		done

		align	16
lpNA:
		movups	xmm2,[edx+ebx]
		movups	xmm3,[edx+ebx+16]
		maxps	xmm2,xmm0
		maxps	xmm3,xmm0
		prefetchnta	[edx+ebx+64]
		minps	xmm2,xmm1
		minps	xmm3,xmm1
		movaps	[edi+ebx],xmm2
		movaps	[edi+ebx+16],xmm3
		add		ebx,16*2
		jl		lpNA
done:
	}

	for ( i = 0; i < pre; i++ ) {
		if ( src[i] < min )
			dst[i] = min;
		else if ( src[i] > max )
			dst[i] = max;
		else
			dst[i] = src[i];
	}

	for( i = count - post; i < count; i++ ) {
		if ( src[i] < min )
			dst[i] = min;
		else if ( src[i] > max )
			dst[i] = max;
		else
			dst[i] = src[i];
	}
}

/*
============
CSIMD_SSE::clampMin
============
*/
void VPCALL CSIMD_SSE::clampMin( float *dst, const float *src, const float min, const int count ) {
	int	i, pre, post;

	__asm
	{
		movss	xmm0,min
		shufps	xmm0,xmm0,0

		KFLOATINITDS( dst, src, count, pre, post )
		and		eax,15
		jne		lpNA
		jmp		lpA
		align	16
lpA:
		movaps	xmm2,[edx+ebx]
		movaps	xmm3,[edx+ebx+16]
		maxps	xmm2,xmm0
		prefetchnta	[edx+ebx+64]
		maxps	xmm3,xmm0
		movaps	[edi+ebx],xmm2
		movaps	[edi+ebx+16],xmm3
		add		ebx,16*2
		jl		lpA
		jmp		done

		align	16
lpNA:
		movups	xmm2,[edx+ebx]
		movups	xmm3,[edx+ebx+16]
		maxps	xmm2,xmm0
		prefetchnta	[edx+ebx+64]
		maxps	xmm3,xmm0
		movaps	[edi+ebx],xmm2
		movaps	[edi+ebx+16],xmm3
		add		ebx,16*2
		jl		lpNA
done:
	}

	for( i = 0; i < pre; i++ ) {
		if ( src[i] < min )
			dst[i] = min;
		else
			dst[i] = src[i];
	}
	for( i = count - post; i < count; i++ ) {
		if ( src[i] < min )
			dst[i] = min;
		else
			dst[i] = src[i];
	}
}

/*
============
CSIMD_SSE::clampMax
============
*/
void VPCALL CSIMD_SSE::clampMax( float *dst, const float *src, const float max, const int count ) {
	int	i, pre, post;

	__asm
	{
		movss	xmm1,max
		shufps	xmm1,xmm1,0

		KFLOATINITDS( dst, src, count, pre, post )
		and		eax,15
		jne		lpNA
		jmp		lpA
		align	16
lpA:
		movaps	xmm2,[edx+ebx]
		movaps	xmm3,[edx+ebx+16]
		minps	xmm2,xmm1
		prefetchnta	[edx+ebx+64]
		minps	xmm3,xmm1
		movaps	[edi+ebx],xmm2
		movaps	[edi+ebx+16],xmm3
		add		ebx,16*2
		jl		lpA
		jmp		done

		align	16
lpNA:
		movups	xmm2,[edx+ebx]
		movups	xmm3,[edx+ebx+16]
		minps	xmm2,xmm1
		prefetchnta	[edx+ebx+64]
		minps	xmm3,xmm1
		movaps	[edi+ebx],xmm2
		movaps	[edi+ebx+16],xmm3
		add		ebx,16*2
		jl		lpNA
done:
	}

	for( i = 0; i < pre; i++ ) {
		if ( src[i] > max )
			dst[i] = max;
		else
			dst[i] = src[i];
	}

	for( i = count - post; i < count; i++ ) {
		if ( src[i] > max )
			dst[i] = max;
		else
			dst[i] = src[i];
	}
}

/*
============
CSIMD_SSE::zero16
============
*/
void VPCALL CSIMD_SSE::zero16( float *dst, const int count ) {
	__asm {
		mov		edx, dst
		mov		eax, count
		add		eax, 3
		shr		eax, 2
		jz		doneZero16
		shl		eax, 4
		add		edx, eax
		neg		eax
		xorps	xmm0, xmm0
	loopZero16:
		movaps	[edx+eax], xmm0
		add		eax, 16
		jl		loopZero16
	doneZero16:
	}
}

/*
============
CSIMD_SSE::negate16
============
*/
void VPCALL CSIMD_SSE::negate16( float *dst, const int count ) {
	__asm {
		mov		edx, dst
		mov		eax, count
		add		eax, 3
		shr		eax, 2
		jz		doneNegate16
		shl		eax, 4
		add		edx, eax
		neg		eax
		movss	xmm0, SIMD_SP_signBitMask
		shufps	xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
	loopNegate16:
		movaps	xmm1, [edx+eax]
		xorps	xmm1, xmm0
		movaps	[edx+eax], xmm1
		add		eax, 16
		jl		loopNegate16
	doneNegate16:
	}
}

/*
============
CSIMD_SSE::copy16
============
*/



void VPCALL CSIMD_SSE::copy16( float *dst, const float *src, const int count ) {
//	Debug::debug(Debug::math, __FUNCTION__) << "size of float: " << sizeof(float) << endl;

	__asm {
		mov		ecx, src
		mov		edx, dst
		mov		eax, count
		add		eax, 3
		shr		eax, 2
		jz		doneCopy16
		shl		eax, 4
		add		ecx, eax
		add		edx, eax
		neg		eax
	loopCopy16:

		movaps	xmm0, [ecx+eax]
		movaps	[edx+eax], xmm0
		add		eax, 16
		jl		loopCopy16
	doneCopy16:
	}
}

/*
============
CSIMD_SSE::add16
============
*/
void VPCALL CSIMD_SSE::add16( float *dst, const float *src1, const float *src2, const int count ) {
	__asm {
		mov		ecx, src1
		mov		edx, src2
		mov		esi, dst
		mov		eax, count
		add		eax, 3
		shr		eax, 2
		jz		doneAdd16
		shl		eax, 4
		add		esi, eax
		add		ecx, eax
		add		edx, eax
		neg		eax
	loopAdd16:
		movaps	xmm0, [ecx+eax]
		addps	xmm0, [edx+eax]
		movaps	[esi+eax], xmm0
		add		eax, 16
		jl		loopAdd16
	doneAdd16:
	}
}

/*
============
CSIMD_SSE::sub16
============
*/
void VPCALL CSIMD_SSE::sub16( float *dst, const float *src1, const float *src2, const int count ) {
	__asm {
		mov		ecx, src1
		mov		edx, src2
		mov		esi, dst
		mov		eax, count
		add		eax, 3
		shr		eax, 2
		jz		doneSub16
		shl		eax, 4
		add		esi, eax
		add		ecx, eax
		add		edx, eax
		neg		eax
	loopSub16:
		movaps	xmm0, [ecx+eax]
		subps	xmm0, [edx+eax]
		movaps	[esi+eax], xmm0
		add		eax, 16
		jl		loopSub16
	doneSub16:
	}
}

/*
============
CSIMD_SSE::mul16
============
*/
void VPCALL CSIMD_SSE::mul16( float *dst, const float *src1, const float constant, const int count ) {
	__asm {
		mov		ecx, dst
		mov		edx, src1
		mov		eax, count
		add		eax, 3
		shr		eax, 2
		jz		doneMulScalar16
		movss	xmm1, constant
		shl		eax, 4
		add		ecx, eax
		add		edx, eax
		neg		eax
		shufps	xmm1, xmm1, 0x00
	loopMulScalar16:
		movaps	xmm0, [edx+eax]
		mulps	xmm0, xmm1
		movaps	[ecx+eax], xmm0
		add		eax, 16
		jl		loopMulScalar16
	doneMulScalar16:
	}
}

/*
============
CSIMD_SSE::addAssign16
============
*/
void VPCALL CSIMD_SSE::addAssign16( float *dst, const float *src, const int count ) {
	__asm {
		mov		ecx, dst
		mov		edx, src
		mov		eax, count
		add		eax, 3
		shr		eax, 2
		jz		doneAddAssign16
		shl		eax, 4
		add		ecx, eax
		add		edx, eax
		neg		eax
	loopAddAssign16:
		movaps	xmm0, [ecx+eax]
		addps	xmm0, [edx+eax]
		movaps	[ecx+eax], xmm0
		add		eax, 16
		jl		loopAddAssign16
	doneAddAssign16:
	}
}

/*
============
CSIMD_SSE::subAssign16
============
*/
void VPCALL CSIMD_SSE::subAssign16( float *dst, const float *src, const int count ) {
	__asm {
		mov		ecx, dst
		mov		edx, src
		mov		eax, count
		add		eax, 3
		shr		eax, 2
		jz		doneSubAssign16
		shl		eax, 4
		add		ecx, eax
		add		edx, eax
		neg		eax
	loopSubAssign16:
		movaps	xmm0, [ecx+eax]
		subps	xmm0, [edx+eax]
		movaps	[ecx+eax], xmm0
		add		eax, 16
		jl		loopSubAssign16
	doneSubAssign16:
	}
}

/*
============
CSIMD_SSE::mulAssign16
============
*/
void VPCALL CSIMD_SSE::mulAssign16( float *dst, const float constant, const int count ) {
	__asm {
		mov		ecx, dst
		mov		eax, count
		add		eax, 3
		shr		eax, 2
		jz		doneMulAssign16
		movss	xmm1, constant
		shl		eax, 4
		add		ecx, eax
		neg		eax
		shufps	xmm1, xmm1, 0x00
	loopMulAssign16:
		movaps	xmm0, [ecx+eax]
		mulps	xmm0, xmm1
		movaps	[ecx+eax], xmm0
		add		eax, 16
		jl		loopMulAssign16
	doneMulAssign16:
	}
}

bool VPCALL CSIMD_SSE::matX_inverse_4x4(float* src)
{
sf_m128 minor0, minor1, minor2, minor3;
sf_m128 row0, row1, row2, row3;
sf_m128 det, tmp1;
tmp1 = _mm_loadh_pi(_mm_loadl_pi(tmp1, (__m64*)(src) ), (__m64*)(src+ 4));
row1 = _mm_loadh_pi(_mm_loadl_pi(row1, (__m64*)(src+8)), (__m64*)(src+12));
row0 = _mm_shuffle_ps(tmp1, row1, 0x88);
row1 = _mm_shuffle_ps(row1, tmp1, 0xDD);
tmp1 = _mm_loadh_pi(_mm_loadl_pi(tmp1, (__m64*)(src+ 2)), (__m64*)(src+ 6));
row3 = _mm_loadh_pi(_mm_loadl_pi(row3, (__m64*)(src+10)), (__m64*)(src+14));
row2 = _mm_shuffle_ps(tmp1, row3, 0x88);
row3 = _mm_shuffle_ps(row3, tmp1, 0xDD);

// -----------------------------------------------

tmp1 = _mm_mul_ps(row2, row3);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
minor0 = _mm_mul_ps(row1, tmp1);
minor1 = _mm_mul_ps(row0, tmp1);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor0 = _mm_sub_ps(_mm_mul_ps(row1, tmp1), minor0);
minor1 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor1);
minor1 = _mm_shuffle_ps(minor1, minor1, 0x4E);

// -----------------------------------------------

tmp1 = _mm_mul_ps(row1, row2);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
minor0 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor0);
minor3 = _mm_mul_ps(row0, tmp1);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor0 = _mm_sub_ps(minor0, _mm_mul_ps(row3, tmp1));
minor3 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor3);
minor3 = _mm_shuffle_ps(minor3, minor3, 0x4E);

// -----------------------------------------------

tmp1 = _mm_mul_ps(_mm_shuffle_ps(row1, row1, 0x4E), row3);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
row2 = _mm_shuffle_ps(row2, row2, 0x4E);
minor0 = _mm_add_ps(_mm_mul_ps(row2, tmp1), minor0);
minor2 = _mm_mul_ps(row0, tmp1);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor0 = _mm_sub_ps(minor0, _mm_mul_ps(row2, tmp1));
minor2 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor2);
minor2 = _mm_shuffle_ps(minor2, minor2, 0x4E);

// -----------------------------------------------

tmp1 = _mm_mul_ps(row0, row1);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
minor2 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor2);
minor3 = _mm_sub_ps(_mm_mul_ps(row2, tmp1), minor3);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor2 = _mm_sub_ps(_mm_mul_ps(row3, tmp1), minor2);
minor3 = _mm_sub_ps(minor3, _mm_mul_ps(row2, tmp1));

// -----------------------------------------------

tmp1 = _mm_mul_ps(row0, row3);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
minor1 = _mm_sub_ps(minor1, _mm_mul_ps(row2, tmp1));
minor2 = _mm_add_ps(_mm_mul_ps(row1, tmp1), minor2);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor1 = _mm_add_ps(_mm_mul_ps(row2, tmp1), minor1);
minor2 = _mm_sub_ps(minor2, _mm_mul_ps(row1, tmp1));

// -----------------------------------------------

tmp1 = _mm_mul_ps(row0, row2);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
minor1 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor1);
minor3 = _mm_sub_ps(minor3, _mm_mul_ps(row1, tmp1));
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor1 = _mm_sub_ps(minor1, _mm_mul_ps(row3, tmp1));
minor3 = _mm_add_ps(_mm_mul_ps(row1, tmp1), minor3);

// -----------------------------------------------

det = _mm_mul_ps(row0, minor0);
det = _mm_add_ps(_mm_shuffle_ps(det, det, 0x4E), det);
det = _mm_add_ss(_mm_shuffle_ps(det, det, 0xB1), det);
tmp1 = _mm_rcp_ss(det);
det = _mm_sub_ss(_mm_add_ss(tmp1, tmp1), _mm_mul_ss(det, _mm_mul_ss(tmp1, tmp1)));
det = _mm_shuffle_ps(det, det, 0x00);
minor0 = _mm_mul_ps(det, minor0);
_mm_storel_pi((__m64*)(src), minor0);
_mm_storeh_pi((__m64*)(src+2), minor0);
minor1 = _mm_mul_ps(det, minor1);
_mm_storel_pi((__m64*)(src+4), minor1);
_mm_storeh_pi((__m64*)(src+6), minor1);
minor2 = _mm_mul_ps(det, minor2);
_mm_storel_pi((__m64*)(src+ 8), minor2);
_mm_storeh_pi((__m64*)(src+10), minor2);
minor3 = _mm_mul_ps(det, minor3);
_mm_storel_pi((__m64*)(src+12), minor3);
_mm_storeh_pi((__m64*)(src+14), minor3);
return true; //por enquanto
}





/*
============
CSIMD_SSE::matX_MultiplyVecX

	optimizes the following matrix multiplications:

	NxN * Nx1
	Nx6 * 6x1
	6xN * Nx1

	with N in the range [1-6]
============
*/
void VPCALL CSIMD_SSE::matX_MultiplyVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec ) {
#define STORE1( offset, reg1, reg2 )		\
	__asm movss		[eax+offset], reg1
#define STORE2LO( offset, reg1, reg2 )		\
	__asm movlps	[eax+offset], reg1
#define STORE2HI( offset, reg1, reg2 )		\
	__asm movhps	[eax+offset], reg1
#define STORE4( offset, reg1, reg2 )		\
	__asm movlps	[eax+offset], reg1		\
	__asm movhps	[eax+offset+8], reg1
#define STOREC		=

	int numRows;
	const float *mPtr, *vPtr;
	float *dstPtr;

	SMF_ASSERT( vec.getSize() >= mat.getNumColumns() );
	SMF_ASSERT( dst.getSize() >= mat.getNumRows() );

	mPtr = mat.toFloatPtr();
	vPtr = vec.toFloatPtr();
	dstPtr = dst.toFloatPtr();
	numRows = mat.getNumRows();
	switch( mat.getNumColumns() ) {
		case 1: {
			switch( numRows ) {
				case 1: {		// 1x1 * 1x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movss		xmm0, [esi]
						mulss		xmm0, [edi]
						STORE1( 0, xmm0, xmm1 )
					}
					return;
				}
				case 6: {		// 6x1 * 1x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movss		xmm0, [esi]
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
						movaps		xmm1, xmm0
						mulps		xmm0, [edi]
						mulps		xmm1, [edi+16]
						STORE4( 0, xmm0, xmm2 )
						STORE2LO( 16, xmm1, xmm2 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numRows; i++ ) {
						dstPtr[i] STOREC mPtr[0] * vPtr[0];
						mPtr++;
					}
					return;
				}
			}
			break;
		}
		case 2: {
			switch( numRows ) {
				case 2: {		// 2x2 * 2x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movss		xmm0, [esi]
						movss		xmm1, [esi+4]
						movss		xmm2, [edi]
						mulss		xmm2, xmm0
						movss		xmm3, [edi+4]
						mulss		xmm3, xmm1
						addss		xmm2, xmm3
						STORE1( 0, xmm2, xmm4 )
						mulss		xmm0, [edi+8]
						mulss		xmm1, [edi+8+4]
						addss		xmm0, xmm1
						STORE1( 4, xmm0, xmm4 )
					}
					return;
				}
				case 6: {		// 6x2 * 2x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm7, [esi]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 1, 0, 1 )
						movaps		xmm0, [edi]
						mulps		xmm0, xmm7
						movaps		xmm1, [edi+16]
						mulps		xmm1, xmm7
						movaps		xmm2, xmm0
						shufps		xmm0, xmm1, R_SHUFFLEPS( 0, 2, 0, 2 )
						shufps		xmm2, xmm1, R_SHUFFLEPS( 1, 3, 1, 3 )
						movaps		xmm3, [edi+32]
						addps		xmm0, xmm2
						mulps		xmm3, xmm7
						STORE4( 0, xmm0, xmm4 )
						shufps		xmm3, xmm3, R_SHUFFLEPS( 0, 2, 1, 3 )
						movhlps		xmm1, xmm3
						addps		xmm3, xmm1
						STORE2LO( 16, xmm3, xmm4 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numRows; i++ ) {
						dstPtr[i] STOREC mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1];
						mPtr += 2;
					}
					return;
				}
			}
			break;
		}
		case 3: {
			switch( numRows ) {
				case 3: {		// 3x3 * 3x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movss		xmm0, [esi]
						movss		xmm4, [edi]
						mulss		xmm4, xmm0
						movss		xmm1, [esi+4]
						movss		xmm5, [edi+4]
						mulss		xmm5, xmm1
						addss		xmm4, xmm5
						movss		xmm2, [esi+8]
						movss		xmm6, [edi+8]
						mulss		xmm6, xmm2
						addss		xmm4, xmm6
						movss		xmm3, [edi+12]
						mulss		xmm3, xmm0
						STORE1( 0, xmm4, xmm7 );
						movss		xmm5, [edi+12+4]
						mulss		xmm5, xmm1
						addss		xmm3, xmm5
						movss		xmm6, [edi+12+8]
						mulss		xmm6, xmm2
						addss		xmm3, xmm6
						mulss		xmm0, [edi+24]
						mulss		xmm1, [edi+24+4]
						STORE1( 4, xmm3, xmm7 );
						addss		xmm0, xmm1
						mulss		xmm2, [edi+24+8]
						addss		xmm0, xmm2
						STORE1( 8, xmm0, xmm7 );
					}
					return;
				}
				case 6: {		// 6x3 * 3x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movss		xmm5, [esi]
						shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )
						movss		xmm6, [esi+4]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						movss		xmm7, [esi+8]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
						movaps		xmm0, [edi]								// xmm0 = 0, 1, 2, 3
						movlps		xmm1, [edi+4*4]
						shufps		xmm1, xmm0, R_SHUFFLEPS( 0, 1, 1, 2 )	// xmm1 = 4, 5, 1, 2
						movlps		xmm2, [edi+6*4]
						movhps		xmm2, [edi+8*4]							// xmm2 = 6, 7, 8, 9
						shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 3, 0, 3 )	// xmm0 = 0, 3, 6, 9
						mulps		xmm0, xmm5
						movlps		xmm3, [edi+10*4]
						shufps		xmm2, xmm3, R_SHUFFLEPS( 1, 2, 0, 1 )	// xmm2 = 7, 8, 10, 11
						movaps		xmm3, xmm1
						shufps		xmm1, xmm2, R_SHUFFLEPS( 2, 0, 0, 2 )	// xmm1 = 1, 4, 7, 10
						mulps		xmm1, xmm6
						shufps		xmm3, xmm2, R_SHUFFLEPS( 3, 1, 1, 3 )	// xmm3 = 2, 5, 8, 11
						mulps		xmm3, xmm7
						addps		xmm0, xmm1
						addps		xmm0, xmm3
						STORE4( 0, xmm0, xmm4 )
						movss		xmm1, [edi+12*4]
						mulss		xmm1, xmm5
						movss		xmm2, [edi+13*4]
						mulss		xmm2, xmm6
						movss		xmm3, [edi+14*4]
						mulss		xmm3, xmm7
						addss		xmm1, xmm2
						addss		xmm1, xmm3
						STORE1( 16, xmm1, xmm4 )
						mulss		xmm5, [edi+15*4]
						mulss		xmm6, [edi+16*4]
						mulss		xmm7, [edi+17*4]
						addss		xmm5, xmm6
						addss		xmm5, xmm7
						STORE1( 20, xmm5, xmm4 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numRows; i++ ) {
						dstPtr[i] STOREC mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1] + mPtr[2] * vPtr[2];
						mPtr += 3;
					}
					return;
				}
			}
			break;
		}
		case 4: {
			switch( numRows ) {
				case 4: {		// 4x4 * 4x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm6, qword ptr [esi ]
						movlps		xmm0, qword ptr [edi ]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 1, 0, 1 )
						movhps		xmm0, qword ptr [edi+16]
						mulps		xmm0, xmm6
						movlps		xmm7, qword ptr [esi+ 8]
						movlps		xmm2, qword ptr [edi+ 8]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 1, 0, 1 )
						movhps		xmm2, qword ptr [edi+24]
						mulps		xmm2, xmm7
						movlps		xmm1, qword ptr [edi+32]
						movhps		xmm1, qword ptr [edi+48]
						mulps		xmm1, xmm6
						movlps		xmm3, qword ptr [edi+40]
						addps		xmm0, xmm2
						movhps		xmm3, qword ptr [edi+56]
						mulps		xmm3, xmm7
						movaps		xmm4, xmm0
						addps		xmm1, xmm3
						shufps		xmm4, xmm1, R_SHUFFLEPS( 0, 2, 0, 2 )
						shufps		xmm0, xmm1, R_SHUFFLEPS( 1, 3, 1, 3 )
						addps		xmm0, xmm4
						STORE4( 0, xmm0, xmm2 )
					}
					return;
				}
				case 6: {		// 6x4 * 4x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm6, qword ptr [esi+ 0]
						movlps		xmm0, qword ptr [edi+ 0]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 1, 0, 1 )
						movhps		xmm0, qword ptr [edi+16]
						mulps		xmm0, xmm6
						movlps		xmm7, qword ptr [esi+ 8]
						movlps		xmm2, qword ptr [edi+ 8]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 1, 0, 1 )
						movhps		xmm2, qword ptr [edi+24]
						mulps		xmm2, xmm7
						movlps		xmm1, qword ptr [edi+32]
						movhps		xmm1, qword ptr [edi+48]
						mulps		xmm1, xmm6
						movlps		xmm3, qword ptr [edi+40]
						addps		xmm0, xmm2
						movhps		xmm3, qword ptr [edi+56]
						mulps		xmm3, xmm7
						movaps		xmm4, xmm0
						addps		xmm1, xmm3
						shufps		xmm4, xmm1, R_SHUFFLEPS( 0, 2, 0, 2 )
						shufps		xmm0, xmm1, R_SHUFFLEPS( 1, 3, 1, 3 )
						addps		xmm0, xmm4
						movlps		xmm1, qword ptr [edi+64]
						movhps		xmm1, qword ptr [edi+80]
						STORE4( 0, xmm0, xmm4 )
						mulps		xmm1, xmm6
						movlps		xmm2, qword ptr [edi+72]
						movhps		xmm2, qword ptr [edi+88]
						mulps		xmm2, xmm7
						addps		xmm1, xmm2
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 2, 1, 3 )
						movhlps		xmm3, xmm1
						addps		xmm1, xmm3
						STORE2LO( 16, xmm1, xmm4 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numRows; i++ ) {
						dstPtr[i] STOREC mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1] + mPtr[2] * vPtr[2] + mPtr[3] * vPtr[3];
						mPtr += 4;
					}
					return;
				}
			}
			break;
		}
		case 5: {
			switch( numRows ) {
				case 5: {		// 5x5 * 5x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movss		xmm0, [edi+5*4]							// xmm0 =  5,  X,  X,  X
						movhps		xmm0, [edi+0*4]							// xmm0 =  5,  X,  0,  1
						movss		xmm5, [edi+15*4]						// xmm4 = 15,  X,  X,  X
						movhps		xmm5, [edi+10*4]						// xmm5 = 15,  X, 10, 11
						movaps		xmm1, xmm0								// xmm1 =  5,  X,  0,  1
						shufps		xmm0, xmm5, R_SHUFFLEPS( 2, 0, 2, 0 )	// xmm0 =  0,  5, 10, 15
						movlps		xmm1, [edi+6*4]							// xmm1 =  6,  7,  0,  1
						movlps		xmm5, [edi+16*4]						// xmm5 = 16, 17, 10, 11
						movaps		xmm2, xmm1								// xmm2 =  6,  7,  0,  1
						shufps		xmm1, xmm5, R_SHUFFLEPS( 3, 0, 3, 0 )	// xmm1 =  1,  6, 11, 16
						movhps		xmm2, [edi+2*4]							// xmm2 =  6,  7,  2,  3
						movhps		xmm5, [edi+12*4]						// xmm5 = 16, 17, 12, 13
						movaps		xmm3, xmm2								// xmm3 =  6,  7,  2,  3
						shufps		xmm2, xmm5, R_SHUFFLEPS( 2, 1, 2, 1 )	// xmm2 =  2,  7, 12, 17
						movlps		xmm3, [edi+8*4]							// xmm3 =  8,  9,  2,  3
						movlps		xmm5, [edi+18*4]						// xmm5 = 18, 19, 12, 13
						movss		xmm4, [edi+4*4]							// xmm4 =  4,  X,  X,  X
						movlhps		xmm4, xmm3								// xmm4 =  4,  X,  8,  9
						shufps		xmm3, xmm5, R_SHUFFLEPS( 3, 0, 3, 0 )	// xmm3 =  3,  8, 13, 18
						movhps		xmm5, [edi+14*4]						// xmm6 = 18, 19, 14, 15
						shufps		xmm4, xmm5, R_SHUFFLEPS( 0, 3, 2, 1 )	// xmm4 =  4,  9, 14, 19
						movss		xmm7, [esi+0*4]
						shufps		xmm7, xmm7, 0
						mulps		xmm0, xmm7
						movss		xmm5, [esi+1*4]
						shufps		xmm5, xmm5, 0
						mulps		xmm1, xmm5
						addps		xmm0, xmm1
						movss		xmm6, [esi+2*4]
						shufps		xmm6, xmm6, 0
						mulps		xmm2, xmm6
						addps		xmm0, xmm2
						movss		xmm1, [esi+3*4]
						shufps		xmm1, xmm1, 0
						mulps		xmm3, xmm1
						addps		xmm0, xmm3
						movss		xmm2, [esi+4*4]
						shufps		xmm2, xmm2, 0
						mulps		xmm4, xmm2
						addps		xmm0, xmm4
						mulss		xmm7, [edi+20*4]
						mulss		xmm5, [edi+21*4]
						addps		xmm7, xmm5
						mulss		xmm6, [edi+22*4]
						addps		xmm7, xmm6
						mulss		xmm1, [edi+23*4]
						addps		xmm7, xmm1
						mulss		xmm2, [edi+24*4]
						addps		xmm7, xmm2
						STORE4( 0, xmm0, xmm3 )
						STORE1( 16, xmm7, xmm4 )
					}
					return;
				}
				case 6: {		// 6x5 * 5x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm6, [esi]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 1, 0, 1 )
						movlps		xmm7, [esi+8]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 1, 0, 1 )
						movlps		xmm0, [edi]
						movhps		xmm3, [edi+8]
						movaps		xmm1, [edi+16]
						movlps		xmm2, [edi+32]
						shufps		xmm0, xmm1, R_SHUFFLEPS( 0, 1, 1, 2 )	// xmm0 = 0, 1, 5, 6
						shufps		xmm1, xmm2, R_SHUFFLEPS( 0, 3, 0, 1 )	// xmm1 = 4, 7, 8, 9
						shufps		xmm3, xmm1, R_SHUFFLEPS( 2, 3, 1, 2 )	// xmm3 = 2, 3, 7, 8
						mulps		xmm0, xmm6
						mulps		xmm3, xmm7
						movlps		xmm2, [edi+40]
						addps		xmm0, xmm3								// xmm0 + xmm1
						movhps		xmm5, [edi+40+8]
						movlps		xmm3, [edi+40+16]
						movhps		xmm3, [edi+40+24]
						movlps		xmm4, [edi+40+32]
						shufps		xmm2, xmm3, R_SHUFFLEPS( 0, 1, 1, 2 )	// xmm2 = 10, 11, 15, 16
						shufps		xmm3, xmm4, R_SHUFFLEPS( 0, 3, 0, 1 )	// xmm3 = 14, 17, 18, 19
						shufps		xmm5, xmm3, R_SHUFFLEPS( 2, 3, 1, 2 )	// xmm5 = 12, 13, 17, 18
						mulps		xmm2, xmm6
						mulps		xmm5, xmm7
						addps		xmm2, xmm5								// xmm2 + xmm3
						movss		xmm5, [esi+16]
						shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )
						movaps		xmm4, xmm0
						shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 2, 0, 2 )
						shufps		xmm4, xmm2, R_SHUFFLEPS( 1, 3, 1, 3 )
						shufps		xmm1, xmm3, R_SHUFFLEPS( 0, 3, 0, 3 )
						addps		xmm0, xmm4
						mulps		xmm1, xmm5
						addps		xmm0, xmm1
						STORE4( 0, xmm0, xmm2 )
						movlps		xmm4, [edi+80]
						movhps		xmm3, [edi+80+8]
						movaps		xmm1, [edi+80+16]
						movlps		xmm2, [edi+80+32]
						shufps		xmm4, xmm1, R_SHUFFLEPS( 0, 1, 1, 2 )	// xmm4 = 20, 21, 25, 26
						shufps		xmm1, xmm2, R_SHUFFLEPS( 0, 3, 0, 1 )	// xmm1 = 24, 27, 28, 29
						shufps		xmm3, xmm1, R_SHUFFLEPS( 2, 3, 1, 2 )	// xmm3 = 22, 23, 27, 28
						mulps		xmm4, xmm6
						mulps		xmm3, xmm7
						mulps		xmm1, xmm5
						addps		xmm4, xmm3								// xmm4 + xmm1
						shufps		xmm1, xmm4, R_SHUFFLEPS( 0, 3, 0, 2 )
						shufps		xmm4, xmm4, R_SHUFFLEPS( 1, 3, 0, 0 )
						addps		xmm4, xmm1
						shufps		xmm1, xmm1, R_SHUFFLEPS( 2, 3, 0, 1 )
						addps		xmm4, xmm1
						STORE2LO( 16, xmm4, xmm2 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numRows; i++ ) {
						dstPtr[i] STOREC mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1] + mPtr[2] * vPtr[2] + mPtr[3] * vPtr[3] + mPtr[4] * vPtr[4];
						mPtr += 5;
					}
					return;
				}
			}
			break;
		}
		case 6: {
			switch( numRows ) {
				case 1: {		// 1x6 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
                        mov			eax, dstPtr
						movss		xmm0, [esi]
						mulss		xmm0, [edi]
						movss		xmm1, [esi+4]
						mulss		xmm1, [edi+4]
						movss		xmm2, [esi+8]
						addss		xmm0, xmm1
						mulss		xmm2, [edi+8]
						movss		xmm3, [esi+12]
						addss		xmm0, xmm2
						mulss		xmm3, [edi+12]
						movss		xmm4, [esi+16]
						addss		xmm0, xmm3
						mulss		xmm4, [edi+16]
						movss		xmm5, [esi+20]
						addss		xmm0, xmm4
						mulss		xmm5, [edi+20]
						movss		xmm6, [esi+24]
						addss		xmm0, xmm5
						mulss		xmm6, [edi+24]
						addss		xmm0, xmm6
						STORE1( 0, xmm0, xmm7 )
					}
					return;
				}
				case 2: {		// 2x6 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						// load idVecX
						movlps		xmm4, [esi]
						movhps		xmm4, [esi+8]
						movlps		xmm5, [esi+16]
						movlhps		xmm5, xmm4
						movhlps		xmm6, xmm4
						movlhps		xmm6, xmm5
						// row 0 and 1
						movaps		xmm0, [edi]
						movaps		xmm1, [edi+16]
						movaps		xmm2, [edi+32]
						mulps		xmm0, xmm4
						mulps		xmm1, xmm5
						mulps		xmm2, xmm6
						movhlps		xmm3, xmm0
						movlhps		xmm3, xmm2
						addps		xmm1, xmm3
						shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 1, 2, 3 )
						addps		xmm1, xmm0
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 2, 1, 3 )
						movhlps		xmm0, xmm1
						addps		xmm0, xmm1
						STORE2LO( 0, xmm0, xmm3 )
					}
					return;
				}
				case 3: {		// 3x6 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						// load idVecX
						movlps		xmm4, [esi]
						movhps		xmm4, [esi+8]
						movlps		xmm5, [esi+16]
						movlhps		xmm5, xmm4
						movhlps		xmm6, xmm4
						movlhps		xmm6, xmm5
						// row 0 and 1
						movaps		xmm0, [edi]
						movaps		xmm1, [edi+16]
						movaps		xmm2, [edi+32]
						mulps		xmm0, xmm4
						mulps		xmm1, xmm5
						mulps		xmm2, xmm6
						movhlps		xmm3, xmm0
						movlhps		xmm3, xmm2
						addps		xmm1, xmm3
						shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 1, 2, 3 )
						addps		xmm1, xmm0
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 2, 1, 3 )
						movhlps		xmm0, xmm1
						addps		xmm0, xmm1
						STORE2LO( 0, xmm0, xmm3 )
						// row 2
						movaps		xmm0, [edi+48]
						movaps		xmm1, [edi+48+16]
						mulps		xmm0, xmm4
						mulps		xmm1, xmm5
						addps		xmm0, xmm1
						movhlps		xmm1, xmm0
						addps		xmm0, xmm1
						movaps		xmm1, xmm0
						shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 0, 0, 0 )
						addss		xmm0, xmm1
						STORE1( 8, xmm0, xmm3 )
					}
					return;
				}
				case 4: {		// 4x6 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						// load idVecX
						movlps		xmm4, [esi]
						movhps		xmm4, [esi+8]
						movlps		xmm5, [esi+16]
						movlhps		xmm5, xmm4
						movhlps		xmm6, xmm4
						movlhps		xmm6, xmm5
						// row 0 and 1
						movaps		xmm0, [edi]
						movaps		xmm1, [edi+16]
						movaps		xmm2, [edi+32]
						mulps		xmm0, xmm4
						mulps		xmm1, xmm5
						mulps		xmm2, xmm6
						movhlps		xmm7, xmm0
						movlhps		xmm7, xmm2
						addps		xmm7, xmm1
						shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 1, 2, 3 )
						addps		xmm7, xmm0
						// row 2 and 3
						movaps		xmm0, [edi+48]
						movaps		xmm1, [edi+48+16]
						movaps		xmm2, [edi+48+32]
						mulps		xmm0, xmm4
						mulps		xmm1, xmm5
						mulps		xmm2, xmm6
						movhlps		xmm3, xmm0
						movlhps		xmm3, xmm2
						addps		xmm1, xmm3
						shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 1, 2, 3 )
						addps		xmm1, xmm0
						// last 4 additions for the first 4 rows and store result
						movaps		xmm0, xmm7
						shufps		xmm7, xmm1, R_SHUFFLEPS( 0, 2, 0, 2 )
						shufps		xmm0, xmm1, R_SHUFFLEPS( 1, 3, 1, 3 )
						addps		xmm0, xmm7
						STORE4( 0, xmm0, xmm4 )
					}
					return;
				}
				case 5: {		// 5x6 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						// load idVecX
						movlps		xmm4, [esi]
						movhps		xmm4, [esi+8]
						movlps		xmm5, [esi+16]
						movlhps		xmm5, xmm4
						movhlps		xmm6, xmm4
						movlhps		xmm6, xmm5
						// row 0 and 1
						movaps		xmm0, [edi]
						movaps		xmm1, [edi+16]
						movaps		xmm2, [edi+32]
						mulps		xmm0, xmm4
						mulps		xmm1, xmm5
						mulps		xmm2, xmm6
						movhlps		xmm7, xmm0
						movlhps		xmm7, xmm2
						addps		xmm7, xmm1
						shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 1, 2, 3 )
						addps		xmm7, xmm0
						// row 2 and 3
						movaps		xmm0, [edi+48]
						movaps		xmm1, [edi+48+16]
						movaps		xmm2, [edi+48+32]
						mulps		xmm0, xmm4
						mulps		xmm1, xmm5
						mulps		xmm2, xmm6
						movhlps		xmm3, xmm0
						movlhps		xmm3, xmm2
						addps		xmm1, xmm3
						shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 1, 2, 3 )
						addps		xmm1, xmm0
						// last 4 additions for the first 4 rows and store result
						movaps		xmm0, xmm7
						shufps		xmm7, xmm1, R_SHUFFLEPS( 0, 2, 0, 2 )
						shufps		xmm0, xmm1, R_SHUFFLEPS( 1, 3, 1, 3 )
						addps		xmm0, xmm7
						STORE4( 0, xmm0, xmm3 )
						// row 5
						movaps		xmm0, [edi+96]
						movaps		xmm1, [edi+96+16]
						mulps		xmm0, xmm4
						mulps		xmm1, xmm5
						addps		xmm0, xmm1
						movhlps		xmm1, xmm0
						addps		xmm0, xmm1
						movaps		xmm1, xmm0
						shufps		xmm1, xmm1, 0x01
						addss		xmm0, xmm1
						STORE1( 16, xmm0, xmm3 )
					}
					return;
				}
				case 6: {		// 6x6 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm7, qword ptr [esi]
						movlps		xmm6, qword ptr [esi+8]
						shufps		xmm7, xmm7, 0x44
						shufps		xmm6, xmm6, 0x44
						movlps		xmm0, qword ptr [edi    ]
						movhps		xmm0, qword ptr [edi+ 24]
						mulps		xmm0, xmm7
						movlps		xmm3, qword ptr [edi+  8]
						movhps		xmm3, qword ptr [edi+ 32]
						mulps		xmm3, xmm6
						movlps		xmm1, qword ptr [edi+ 48]
						movhps		xmm1, qword ptr [edi+ 72]
						mulps		xmm1, xmm7
						movlps		xmm2, qword ptr [edi+ 96]
						movhps		xmm2, qword ptr [edi+120]
						mulps		xmm2, xmm7
						movlps		xmm4, qword ptr [edi+ 56]
						movhps		xmm4, qword ptr [edi+ 80]
						movlps		xmm5, qword ptr [edi+104]
						movhps		xmm5, qword ptr [edi+128]
						mulps		xmm4, xmm6
						movlps		xmm7, qword ptr [esi+16]
						addps		xmm0, xmm3
						shufps		xmm7, xmm7, 0x44
						mulps		xmm5, xmm6
						addps		xmm1, xmm4
						movlps		xmm3, qword ptr [edi+ 16]
						movhps		xmm3, qword ptr [edi+ 40]
						addps		xmm2, xmm5
						movlps		xmm4, qword ptr [edi+ 64]
						movhps		xmm4, qword ptr [edi+ 88]
						mulps		xmm3, xmm7
						movlps		xmm5, qword ptr [edi+112]
						movhps		xmm5, qword ptr [edi+136]
						addps		xmm0, xmm3
						mulps		xmm4, xmm7
						mulps		xmm5, xmm7
						addps		xmm1, xmm4
						addps		xmm2, xmm5
						movaps		xmm6, xmm0
						shufps		xmm0, xmm1, 0x88
						shufps		xmm6, xmm1, 0xDD
						movaps		xmm7, xmm2
						shufps		xmm7, xmm2, 0x88
						shufps		xmm2, xmm2, 0xDD
						addps		xmm0, xmm6
						addps		xmm2, xmm7
						STORE4( 0, xmm0, xmm3 )
						STORE2LO( 16, xmm2, xmm4 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numRows; i++ ) {
						dstPtr[i] STOREC mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1] + mPtr[2] * vPtr[2] +
									mPtr[3] * vPtr[3] + mPtr[4] * vPtr[4] + mPtr[5] * vPtr[5];
						mPtr += 6;
					}
					return;
				}
			}
			break;
		}
		default: {
			int numColumns = mat.getNumColumns();
			for ( int i = 0; i < numRows; i++ ) {
				float sum = mPtr[0] * vPtr[0];
				for ( int j = 1; j < numColumns; j++ ) {
					sum += mPtr[j] * vPtr[j];
				}
				dstPtr[i] STOREC sum;
				mPtr += numColumns;
			}
			break;
		}
	}

#undef STOREC
#undef STORE4
#undef STORE2HI
#undef STORE2LO
#undef STORE1

}

/*
============
CSIMD_SSE::matX_MultiplyAddVecX

	optimizes the following matrix multiplications:

	NxN * Nx1
	Nx6 * 6x1
	6xN * Nx1

	with N in the range [1-6]
============
*/
void VPCALL CSIMD_SSE::matX_MultiplyAddVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec ) {
#define STORE1( offset, reg1, reg2 )		\
	__asm movss		reg2, [eax+offset]		\
	__asm addss		reg2, reg1				\
	__asm movss		[eax+offset], reg2
#define STORE2LO( offset, reg1, reg2 )		\
	__asm movlps	reg2, [eax+offset]		\
	__asm addps		reg2, reg1				\
	__asm movlps	[eax+offset], reg2
#define STORE2HI( offset, reg1, reg2 )		\
	__asm movhps	reg2, [eax+offset]		\
	__asm addps		reg2, reg1				\
	__asm movhps	[eax+offset], reg2
#define STORE4( offset, reg1, reg2 )		\
	__asm movlps	reg2, [eax+offset]		\
	__asm movhps	reg2, [eax+offset+8]	\
	__asm addps		reg2, reg1				\
	__asm movlps	[eax+offset], reg2		\
	__asm movhps	[eax+offset+8], reg2
#define STOREC		+=

	int numRows;
	const float *mPtr, *vPtr;
	float *dstPtr;

	SMF_ASSERT( vec.getSize() >= mat.getNumColumns() );
	SMF_ASSERT( dst.getSize() >= mat.getNumRows() );

	mPtr = mat.toFloatPtr();
	vPtr = vec.toFloatPtr();
	dstPtr = dst.toFloatPtr();
	numRows = mat.getNumRows();
	switch( mat.getNumColumns() ) {
		case 1: {
			switch( numRows ) {
				case 1: {		// 1x1 * 1x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movss		xmm0, [esi]
						mulss		xmm0, [edi]
						STORE1( 0, xmm0, xmm1 )
					}
					return;
				}
				case 6: {		// 6x1 * 1x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movss		xmm0, [esi]
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
						movaps		xmm1, xmm0
						mulps		xmm0, [edi]
						mulps		xmm1, [edi+16]
						STORE4( 0, xmm0, xmm2 )
						STORE2LO( 16, xmm1, xmm2 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numRows; i++ ) {
						dstPtr[i] STOREC mPtr[0] * vPtr[0];
						mPtr++;
					}
					return;
				}
			}
			break;
		}
		case 2: {
			switch( numRows ) {
				case 2: {		// 2x2 * 2x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movss		xmm0, [esi]
						movss		xmm1, [esi+4]
						movss		xmm2, [edi]
						mulss		xmm2, xmm0
						movss		xmm3, [edi+4]
						mulss		xmm3, xmm1
						addss		xmm2, xmm3
						STORE1( 0, xmm2, xmm4 )
						mulss		xmm0, [edi+8]
						mulss		xmm1, [edi+8+4]
						addss		xmm0, xmm1
						STORE1( 4, xmm0, xmm4 )
					}
					return;
				}
				case 6: {		// 6x2 * 2x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm7, [esi]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 1, 0, 1 )
						movaps		xmm0, [edi]
						mulps		xmm0, xmm7
						movaps		xmm1, [edi+16]
						mulps		xmm1, xmm7
						movaps		xmm2, xmm0
						shufps		xmm0, xmm1, R_SHUFFLEPS( 0, 2, 0, 2 )
						shufps		xmm2, xmm1, R_SHUFFLEPS( 1, 3, 1, 3 )
						movaps		xmm3, [edi+32]
						addps		xmm0, xmm2
						mulps		xmm3, xmm7
						STORE4( 0, xmm0, xmm4 )
						shufps		xmm3, xmm3, R_SHUFFLEPS( 0, 2, 1, 3 )
						movhlps		xmm1, xmm3
						addps		xmm3, xmm1
						STORE2LO( 16, xmm3, xmm4 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numRows; i++ ) {
						dstPtr[i] STOREC mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1];
						mPtr += 2;
					}
					return;
				}
			}
			break;
		}
		case 3: {
			switch( numRows ) {
				case 3: {		// 3x3 * 3x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movss		xmm0, [esi]
						movss		xmm4, [edi]
						mulss		xmm4, xmm0
						movss		xmm1, [esi+4]
						movss		xmm5, [edi+4]
						mulss		xmm5, xmm1
						addss		xmm4, xmm5
						movss		xmm2, [esi+8]
						movss		xmm6, [edi+8]
						mulss		xmm6, xmm2
						addss		xmm4, xmm6
						movss		xmm3, [edi+12]
						mulss		xmm3, xmm0
						STORE1( 0, xmm4, xmm7 );
						movss		xmm5, [edi+12+4]
						mulss		xmm5, xmm1
						addss		xmm3, xmm5
						movss		xmm6, [edi+12+8]
						mulss		xmm6, xmm2
						addss		xmm3, xmm6
						mulss		xmm0, [edi+24]
						mulss		xmm1, [edi+24+4]
						STORE1( 4, xmm3, xmm7 );
						addss		xmm0, xmm1
						mulss		xmm2, [edi+24+8]
						addss		xmm0, xmm2
						STORE1( 8, xmm0, xmm7 );
					}
					return;
				}
				case 6: {		// 6x3 * 3x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movss		xmm5, [esi]
						shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )
						movss		xmm6, [esi+4]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						movss		xmm7, [esi+8]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
						movaps		xmm0, [edi]								// xmm0 = 0, 1, 2, 3
						movlps		xmm1, [edi+4*4]
						shufps		xmm1, xmm0, R_SHUFFLEPS( 0, 1, 1, 2 )	// xmm1 = 4, 5, 1, 2
						movlps		xmm2, [edi+6*4]
						movhps		xmm2, [edi+8*4]							// xmm2 = 6, 7, 8, 9
						shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 3, 0, 3 )	// xmm0 = 0, 3, 6, 9
						mulps		xmm0, xmm5
						movlps		xmm3, [edi+10*4]
						shufps		xmm2, xmm3, R_SHUFFLEPS( 1, 2, 0, 1 )	// xmm2 = 7, 8, 10, 11
						movaps		xmm3, xmm1
						shufps		xmm1, xmm2, R_SHUFFLEPS( 2, 0, 0, 2 )	// xmm1 = 1, 4, 7, 10
						mulps		xmm1, xmm6
						shufps		xmm3, xmm2, R_SHUFFLEPS( 3, 1, 1, 3 )	// xmm3 = 2, 5, 8, 11
						mulps		xmm3, xmm7
						addps		xmm0, xmm1
						addps		xmm0, xmm3
						STORE4( 0, xmm0, xmm4 )
						movss		xmm1, [edi+12*4]
						mulss		xmm1, xmm5
						movss		xmm2, [edi+13*4]
						mulss		xmm2, xmm6
						movss		xmm3, [edi+14*4]
						mulss		xmm3, xmm7
						addss		xmm1, xmm2
						addss		xmm1, xmm3
						STORE1( 16, xmm1, xmm4 )
						mulss		xmm5, [edi+15*4]
						mulss		xmm6, [edi+16*4]
						mulss		xmm7, [edi+17*4]
						addss		xmm5, xmm6
						addss		xmm5, xmm7
						STORE1( 20, xmm5, xmm4 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numRows; i++ ) {
						dstPtr[i] STOREC mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1] + mPtr[2] * vPtr[2];
						mPtr += 3;
					}
					return;
				}
			}
			break;
		}
		case 4: {
			switch( numRows ) {
				case 4: {		// 4x4 * 4x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm6, qword ptr [esi ]
						movlps		xmm0, qword ptr [edi ]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 1, 0, 1 )
						movhps		xmm0, qword ptr [edi+16]
						mulps		xmm0, xmm6
						movlps		xmm7, qword ptr [esi+ 8]
						movlps		xmm2, qword ptr [edi+ 8]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 1, 0, 1 )
						movhps		xmm2, qword ptr [edi+24]
						mulps		xmm2, xmm7
						movlps		xmm1, qword ptr [edi+32]
						movhps		xmm1, qword ptr [edi+48]
						mulps		xmm1, xmm6
						movlps		xmm3, qword ptr [edi+40]
						addps		xmm0, xmm2
						movhps		xmm3, qword ptr [edi+56]
						mulps		xmm3, xmm7
						movaps		xmm4, xmm0
						addps		xmm1, xmm3
						shufps		xmm4, xmm1, R_SHUFFLEPS( 0, 2, 0, 2 )
						shufps		xmm0, xmm1, R_SHUFFLEPS( 1, 3, 1, 3 )
						addps		xmm0, xmm4
						STORE4( 0, xmm0, xmm2 )
					}
					return;
				}
				case 6: {		// 6x4 * 4x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm6, qword ptr [esi+ 0]
						movlps		xmm0, qword ptr [edi+ 0]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 1, 0, 1 )
						movhps		xmm0, qword ptr [edi+16]
						mulps		xmm0, xmm6
						movlps		xmm7, qword ptr [esi+ 8]
						movlps		xmm2, qword ptr [edi+ 8]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 1, 0, 1 )
						movhps		xmm2, qword ptr [edi+24]
						mulps		xmm2, xmm7
						movlps		xmm1, qword ptr [edi+32]
						movhps		xmm1, qword ptr [edi+48]
						mulps		xmm1, xmm6
						movlps		xmm3, qword ptr [edi+40]
						addps		xmm0, xmm2
						movhps		xmm3, qword ptr [edi+56]
						mulps		xmm3, xmm7
						movaps		xmm4, xmm0
						addps		xmm1, xmm3
						shufps		xmm4, xmm1, R_SHUFFLEPS( 0, 2, 0, 2 )
						shufps		xmm0, xmm1, R_SHUFFLEPS( 1, 3, 1, 3 )
						addps		xmm0, xmm4
						movlps		xmm1, qword ptr [edi+64]
						movhps		xmm1, qword ptr [edi+80]
						STORE4( 0, xmm0, xmm4 )
						mulps		xmm1, xmm6
						movlps		xmm2, qword ptr [edi+72]
						movhps		xmm2, qword ptr [edi+88]
						mulps		xmm2, xmm7
						addps		xmm1, xmm2
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 2, 1, 3 )
						movhlps		xmm3, xmm1
						addps		xmm1, xmm3
						STORE2LO( 16, xmm1, xmm4 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numRows; i++ ) {
						dstPtr[i] STOREC mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1] + mPtr[2] * vPtr[2] + mPtr[3] * vPtr[3];
						mPtr += 4;
					}
					return;
				}
			}
			break;
		}
		case 5: {
			switch( numRows ) {
				case 5: {		// 5x5 * 5x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movss		xmm0, [edi+5*4]							// xmm0 =  5,  X,  X,  X
						movhps		xmm0, [edi+0*4]							// xmm0 =  5,  X,  0,  1
						movss		xmm5, [edi+15*4]						// xmm4 = 15,  X,  X,  X
						movhps		xmm5, [edi+10*4]						// xmm5 = 15,  X, 10, 11
						movaps		xmm1, xmm0								// xmm1 =  5,  X,  0,  1
						shufps		xmm0, xmm5, R_SHUFFLEPS( 2, 0, 2, 0 )	// xmm0 =  0,  5, 10, 15
						movlps		xmm1, [edi+6*4]							// xmm1 =  6,  7,  0,  1
						movlps		xmm5, [edi+16*4]						// xmm5 = 16, 17, 10, 11
						movaps		xmm2, xmm1								// xmm2 =  6,  7,  0,  1
						shufps		xmm1, xmm5, R_SHUFFLEPS( 3, 0, 3, 0 )	// xmm1 =  1,  6, 11, 16
						movhps		xmm2, [edi+2*4]							// xmm2 =  6,  7,  2,  3
						movhps		xmm5, [edi+12*4]						// xmm5 = 16, 17, 12, 13
						movaps		xmm3, xmm2								// xmm3 =  6,  7,  2,  3
						shufps		xmm2, xmm5, R_SHUFFLEPS( 2, 1, 2, 1 )	// xmm2 =  2,  7, 12, 17
						movlps		xmm3, [edi+8*4]							// xmm3 =  8,  9,  2,  3
						movlps		xmm5, [edi+18*4]						// xmm5 = 18, 19, 12, 13
						movss		xmm4, [edi+4*4]							// xmm4 =  4,  X,  X,  X
						movlhps		xmm4, xmm3								// xmm4 =  4,  X,  8,  9
						shufps		xmm3, xmm5, R_SHUFFLEPS( 3, 0, 3, 0 )	// xmm3 =  3,  8, 13, 18
						movhps		xmm5, [edi+14*4]						// xmm6 = 18, 19, 14, 15
						shufps		xmm4, xmm5, R_SHUFFLEPS( 0, 3, 2, 1 )	// xmm4 =  4,  9, 14, 19
						movss		xmm7, [esi+0*4]
						shufps		xmm7, xmm7, 0
						mulps		xmm0, xmm7
						movss		xmm5, [esi+1*4]
						shufps		xmm5, xmm5, 0
						mulps		xmm1, xmm5
						addps		xmm0, xmm1
						movss		xmm6, [esi+2*4]
						shufps		xmm6, xmm6, 0
						mulps		xmm2, xmm6
						addps		xmm0, xmm2
						movss		xmm1, [esi+3*4]
						shufps		xmm1, xmm1, 0
						mulps		xmm3, xmm1
						addps		xmm0, xmm3
						movss		xmm2, [esi+4*4]
						shufps		xmm2, xmm2, 0
						mulps		xmm4, xmm2
						addps		xmm0, xmm4
						mulss		xmm7, [edi+20*4]
						mulss		xmm5, [edi+21*4]
						addps		xmm7, xmm5
						mulss		xmm6, [edi+22*4]
						addps		xmm7, xmm6
						mulss		xmm1, [edi+23*4]
						addps		xmm7, xmm1
						mulss		xmm2, [edi+24*4]
						addps		xmm7, xmm2
						STORE4( 0, xmm0, xmm3 )
						STORE1( 16, xmm7, xmm4 )
					}
					return;
				}
				case 6: {		// 6x5 * 5x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm6, [esi]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 1, 0, 1 )
						movlps		xmm7, [esi+8]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 1, 0, 1 )
						movlps		xmm0, [edi]
						movhps		xmm3, [edi+8]
						movaps		xmm1, [edi+16]
						movlps		xmm2, [edi+32]
						shufps		xmm0, xmm1, R_SHUFFLEPS( 0, 1, 1, 2 )	// xmm0 = 0, 1, 5, 6
						shufps		xmm1, xmm2, R_SHUFFLEPS( 0, 3, 0, 1 )	// xmm1 = 4, 7, 8, 9
						shufps		xmm3, xmm1, R_SHUFFLEPS( 2, 3, 1, 2 )	// xmm3 = 2, 3, 7, 8
						mulps		xmm0, xmm6
						mulps		xmm3, xmm7
						movlps		xmm2, [edi+40]
						addps		xmm0, xmm3								// xmm0 + xmm1
						movhps		xmm5, [edi+40+8]
						movlps		xmm3, [edi+40+16]
						movhps		xmm3, [edi+40+24]
						movlps		xmm4, [edi+40+32]
						shufps		xmm2, xmm3, R_SHUFFLEPS( 0, 1, 1, 2 )	// xmm2 = 10, 11, 15, 16
						shufps		xmm3, xmm4, R_SHUFFLEPS( 0, 3, 0, 1 )	// xmm3 = 14, 17, 18, 19
						shufps		xmm5, xmm3, R_SHUFFLEPS( 2, 3, 1, 2 )	// xmm5 = 12, 13, 17, 18
						mulps		xmm2, xmm6
						mulps		xmm5, xmm7
						addps		xmm2, xmm5								// xmm2 + xmm3
						movss		xmm5, [esi+16]
						shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )
						movaps		xmm4, xmm0
						shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 2, 0, 2 )
						shufps		xmm4, xmm2, R_SHUFFLEPS( 1, 3, 1, 3 )
						shufps		xmm1, xmm3, R_SHUFFLEPS( 0, 3, 0, 3 )
						addps		xmm0, xmm4
						mulps		xmm1, xmm5
						addps		xmm0, xmm1
						STORE4( 0, xmm0, xmm2 )
						movlps		xmm4, [edi+80]
						movhps		xmm3, [edi+80+8]
						movaps		xmm1, [edi+80+16]
						movlps		xmm2, [edi+80+32]
						shufps		xmm4, xmm1, R_SHUFFLEPS( 0, 1, 1, 2 )	// xmm4 = 20, 21, 25, 26
						shufps		xmm1, xmm2, R_SHUFFLEPS( 0, 3, 0, 1 )	// xmm1 = 24, 27, 28, 29
						shufps		xmm3, xmm1, R_SHUFFLEPS( 2, 3, 1, 2 )	// xmm3 = 22, 23, 27, 28
						mulps		xmm4, xmm6
						mulps		xmm3, xmm7
						mulps		xmm1, xmm5
						addps		xmm4, xmm3								// xmm4 + xmm1
						shufps		xmm1, xmm4, R_SHUFFLEPS( 0, 3, 0, 2 )
						shufps		xmm4, xmm4, R_SHUFFLEPS( 1, 3, 0, 0 )
						addps		xmm4, xmm1
						shufps		xmm1, xmm1, R_SHUFFLEPS( 2, 3, 0, 1 )
						addps		xmm4, xmm1
						STORE2LO( 16, xmm4, xmm2 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numRows; i++ ) {
						dstPtr[i] STOREC mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1] + mPtr[2] * vPtr[2] + mPtr[3] * vPtr[3] + mPtr[4] * vPtr[4];
						mPtr += 5;
					}
					return;
				}
			}
			break;
		}
		case 6: {
			switch( numRows ) {
				case 1: {		// 1x6 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
                        mov			eax, dstPtr
						movss		xmm0, [esi]
						mulss		xmm0, [edi]
						movss		xmm1, [esi+4]
						mulss		xmm1, [edi+4]
						movss		xmm2, [esi+8]
						addss		xmm0, xmm1
						mulss		xmm2, [edi+8]
						movss		xmm3, [esi+12]
						addss		xmm0, xmm2
						mulss		xmm3, [edi+12]
						movss		xmm4, [esi+16]
						addss		xmm0, xmm3
						mulss		xmm4, [edi+16]
						movss		xmm5, [esi+20]
						addss		xmm0, xmm4
						mulss		xmm5, [edi+20]
						movss		xmm6, [esi+24]
						addss		xmm0, xmm5
						mulss		xmm6, [edi+24]
						addss		xmm0, xmm6
						STORE1( 0, xmm0, xmm7 )
					}
					return;
				}
				case 2: {		// 2x6 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						// load CVecXD
						movlps		xmm4, [esi]
						movhps		xmm4, [esi+8]
						movlps		xmm5, [esi+16]
						movlhps		xmm5, xmm4
						movhlps		xmm6, xmm4
						movlhps		xmm6, xmm5
						// row 0 and 1
						movaps		xmm0, [edi]
						movaps		xmm1, [edi+16]
						movaps		xmm2, [edi+32]
						mulps		xmm0, xmm4
						mulps		xmm1, xmm5
						mulps		xmm2, xmm6
						movhlps		xmm3, xmm0
						movlhps		xmm3, xmm2
						addps		xmm1, xmm3
						shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 1, 2, 3 )
						addps		xmm1, xmm0
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 2, 1, 3 )
						movhlps		xmm0, xmm1
						addps		xmm0, xmm1
						STORE2LO( 0, xmm0, xmm3 )
					}
					return;
				}
				case 3: {		// 3x6 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						// load CVecXD
						movlps		xmm4, [esi]
						movhps		xmm4, [esi+8]
						movlps		xmm5, [esi+16]
						movlhps		xmm5, xmm4
						movhlps		xmm6, xmm4
						movlhps		xmm6, xmm5
						// row 0 and 1
						movaps		xmm0, [edi]
						movaps		xmm1, [edi+16]
						movaps		xmm2, [edi+32]
						mulps		xmm0, xmm4
						mulps		xmm1, xmm5
						mulps		xmm2, xmm6
						movhlps		xmm3, xmm0
						movlhps		xmm3, xmm2
						addps		xmm1, xmm3
						shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 1, 2, 3 )
						addps		xmm1, xmm0
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 2, 1, 3 )
						movhlps		xmm0, xmm1
						addps		xmm0, xmm1
						STORE2LO( 0, xmm0, xmm3 )
						// row 2
						movaps		xmm0, [edi+48]
						movaps		xmm1, [edi+48+16]
						mulps		xmm0, xmm4
						mulps		xmm1, xmm5
						addps		xmm0, xmm1
						movhlps		xmm1, xmm0
						addps		xmm0, xmm1
						movaps		xmm1, xmm0
						shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 0, 0, 0 )
						addss		xmm0, xmm1
						STORE1( 8, xmm0, xmm3 )
					}
					return;
				}
				case 4: {		// 4x6 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						// load CVecXD
						movlps		xmm4, [esi]
						movhps		xmm4, [esi+8]
						movlps		xmm5, [esi+16]
						movlhps		xmm5, xmm4
						movhlps		xmm6, xmm4
						movlhps		xmm6, xmm5
						// row 0 and 1
						movaps		xmm0, [edi]
						movaps		xmm1, [edi+16]
						movaps		xmm2, [edi+32]
						mulps		xmm0, xmm4
						mulps		xmm1, xmm5
						mulps		xmm2, xmm6
						movhlps		xmm7, xmm0
						movlhps		xmm7, xmm2
						addps		xmm7, xmm1
						shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 1, 2, 3 )
						addps		xmm7, xmm0
						// row 2 and 3
						movaps		xmm0, [edi+48]
						movaps		xmm1, [edi+48+16]
						movaps		xmm2, [edi+48+32]
						mulps		xmm0, xmm4
						mulps		xmm1, xmm5
						mulps		xmm2, xmm6
						movhlps		xmm3, xmm0
						movlhps		xmm3, xmm2
						addps		xmm1, xmm3
						shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 1, 2, 3 )
						addps		xmm1, xmm0
						// last 4 additions for the first 4 rows and store result
						movaps		xmm0, xmm7
						shufps		xmm7, xmm1, R_SHUFFLEPS( 0, 2, 0, 2 )
						shufps		xmm0, xmm1, R_SHUFFLEPS( 1, 3, 1, 3 )
						addps		xmm0, xmm7
						STORE4( 0, xmm0, xmm4 )
					}
					return;
				}
				case 5: {		// 5x6 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						// load CVecXD
						movlps		xmm4, [esi]
						movhps		xmm4, [esi+8]
						movlps		xmm5, [esi+16]
						movlhps		xmm5, xmm4
						movhlps		xmm6, xmm4
						movlhps		xmm6, xmm5
						// row 0 and 1
						movaps		xmm0, [edi]
						movaps		xmm1, [edi+16]
						movaps		xmm2, [edi+32]
						mulps		xmm0, xmm4
						mulps		xmm1, xmm5
						mulps		xmm2, xmm6
						movhlps		xmm7, xmm0
						movlhps		xmm7, xmm2
						addps		xmm7, xmm1
						shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 1, 2, 3 )
						addps		xmm7, xmm0
						// row 2 and 3
						movaps		xmm0, [edi+48]
						movaps		xmm1, [edi+48+16]
						movaps		xmm2, [edi+48+32]
						mulps		xmm0, xmm4
						mulps		xmm1, xmm5
						mulps		xmm2, xmm6
						movhlps		xmm3, xmm0
						movlhps		xmm3, xmm2
						addps		xmm1, xmm3
						shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 1, 2, 3 )
						addps		xmm1, xmm0
						// last 4 additions for the first 4 rows and store result
						movaps		xmm0, xmm7
						shufps		xmm7, xmm1, R_SHUFFLEPS( 0, 2, 0, 2 )
						shufps		xmm0, xmm1, R_SHUFFLEPS( 1, 3, 1, 3 )
						addps		xmm0, xmm7
						STORE4( 0, xmm0, xmm3 )
						// row 5
						movaps		xmm0, [edi+96]
						movaps		xmm1, [edi+96+16]
						mulps		xmm0, xmm4
						mulps		xmm1, xmm5
						addps		xmm0, xmm1
						movhlps		xmm1, xmm0
						addps		xmm0, xmm1
						movaps		xmm1, xmm0
						shufps		xmm1, xmm1, 0x01
						addss		xmm0, xmm1
						STORE1( 16, xmm0, xmm3 )
					}
					return;
				}
				case 6: {		// 6x6 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm7, qword ptr [esi]
						movlps		xmm6, qword ptr [esi+8]
						shufps		xmm7, xmm7, 0x44
						shufps		xmm6, xmm6, 0x44
						movlps		xmm0, qword ptr [edi    ]
						movhps		xmm0, qword ptr [edi+ 24]
						mulps		xmm0, xmm7
						movlps		xmm3, qword ptr [edi+  8]
						movhps		xmm3, qword ptr [edi+ 32]
						mulps		xmm3, xmm6
						movlps		xmm1, qword ptr [edi+ 48]
						movhps		xmm1, qword ptr [edi+ 72]
						mulps		xmm1, xmm7
						movlps		xmm2, qword ptr [edi+ 96]
						movhps		xmm2, qword ptr [edi+120]
						mulps		xmm2, xmm7
						movlps		xmm4, qword ptr [edi+ 56]
						movhps		xmm4, qword ptr [edi+ 80]
						movlps		xmm5, qword ptr [edi+104]
						movhps		xmm5, qword ptr [edi+128]
						mulps		xmm4, xmm6
						movlps		xmm7, qword ptr [esi+16]
						addps		xmm0, xmm3
						shufps		xmm7, xmm7, 0x44
						mulps		xmm5, xmm6
						addps		xmm1, xmm4
						movlps		xmm3, qword ptr [edi+ 16]
						movhps		xmm3, qword ptr [edi+ 40]
						addps		xmm2, xmm5
						movlps		xmm4, qword ptr [edi+ 64]
						movhps		xmm4, qword ptr [edi+ 88]
						mulps		xmm3, xmm7
						movlps		xmm5, qword ptr [edi+112]
						movhps		xmm5, qword ptr [edi+136]
						addps		xmm0, xmm3
						mulps		xmm4, xmm7
						mulps		xmm5, xmm7
						addps		xmm1, xmm4
						addps		xmm2, xmm5
						movaps		xmm6, xmm0
						shufps		xmm0, xmm1, 0x88
						shufps		xmm6, xmm1, 0xDD
						movaps		xmm7, xmm2
						shufps		xmm7, xmm2, 0x88
						shufps		xmm2, xmm2, 0xDD
						addps		xmm0, xmm6
						addps		xmm2, xmm7
						STORE4( 0, xmm0, xmm3 )
						STORE2LO( 16, xmm2, xmm4 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numRows; i++ ) {
						dstPtr[i] STOREC mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1] + mPtr[2] * vPtr[2] +
									mPtr[3] * vPtr[3] + mPtr[4] * vPtr[4] + mPtr[5] * vPtr[5];
						mPtr += 6;
					}
					return;
				}
			}
			break;
		}
		default: {
			int numColumns = mat.getNumColumns();
			for ( int i = 0; i < numRows; i++ ) {
				float sum = mPtr[0] * vPtr[0];
				for ( int j = 1; j < numColumns; j++ ) {
					sum += mPtr[j] * vPtr[j];
				}
				dstPtr[i] STOREC sum;
				mPtr += numColumns;
			}
			break;
		}
	}

#undef STOREC
#undef STORE4
#undef STORE2HI
#undef STORE2LO
#undef STORE1
}

/*
============
CSIMD_SSE::matX_MultiplySubVecX

	optimizes the following matrix multiplications:

	NxN * Nx1
	Nx6 * 6x1
	6xN * Nx1

	with N in the range [1-6]
============
*/
void VPCALL CSIMD_SSE::matX_MultiplySubVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec ) {
#define STORE1( offset, reg1, reg2 )		\
	__asm movss		reg2, [eax+offset]		\
	__asm subss		reg2, reg1				\
	__asm movss		[eax+offset], reg2
#define STORE2LO( offset, reg1, reg2 )		\
	__asm movlps	reg2, [eax+offset]		\
	__asm subps		reg2, reg1				\
	__asm movlps	[eax+offset], reg2
#define STORE2HI( offset, reg1, reg2 )		\
	__asm movhps	reg2, [eax+offset]		\
	__asm subps		reg2, reg1				\
	__asm movhps	[eax+offset], reg2
#define STORE4( offset, reg1, reg2 )		\
	__asm movlps	reg2, [eax+offset]		\
	__asm movhps	reg2, [eax+offset+8]	\
	__asm subps		reg2, reg1				\
	__asm movlps	[eax+offset], reg2		\
	__asm movhps	[eax+offset+8], reg2
#define STOREC		-=

	int numRows;
	const float *mPtr, *vPtr;
	float *dstPtr;

	SMF_ASSERT( vec.getSize() >= mat.getNumColumns() );
	SMF_ASSERT( dst.getSize() >= mat.getNumRows() );

	mPtr = mat.toFloatPtr();
	vPtr = vec.toFloatPtr();
	dstPtr = dst.toFloatPtr();
	numRows = mat.getNumRows();
	switch( mat.getNumColumns() ) {
		case 1: {
			switch( numRows ) {
				case 1: {		// 1x1 * 1x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movss		xmm0, [esi]
						mulss		xmm0, [edi]
						STORE1( 0, xmm0, xmm1 )
					}
					return;
				}
				case 6: {		// 6x1 * 1x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movss		xmm0, [esi]
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
						movaps		xmm1, xmm0
						mulps		xmm0, [edi]
						mulps		xmm1, [edi+16]
						STORE4( 0, xmm0, xmm2 )
						STORE2LO( 16, xmm1, xmm2 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numRows; i++ ) {
						dstPtr[i] STOREC mPtr[0] * vPtr[0];
						mPtr++;
					}
					return;
				}
			}
			break;
		}
		case 2: {
			switch( numRows ) {
				case 2: {		// 2x2 * 2x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movss		xmm0, [esi]
						movss		xmm1, [esi+4]
						movss		xmm2, [edi]
						mulss		xmm2, xmm0
						movss		xmm3, [edi+4]
						mulss		xmm3, xmm1
						addss		xmm2, xmm3
						STORE1( 0, xmm2, xmm4 )
						mulss		xmm0, [edi+8]
						mulss		xmm1, [edi+8+4]
						addss		xmm0, xmm1
						STORE1( 4, xmm0, xmm4 )
					}
					return;
				}
				case 6: {		// 6x2 * 2x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm7, [esi]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 1, 0, 1 )
						movaps		xmm0, [edi]
						mulps		xmm0, xmm7
						movaps		xmm1, [edi+16]
						mulps		xmm1, xmm7
						movaps		xmm2, xmm0
						shufps		xmm0, xmm1, R_SHUFFLEPS( 0, 2, 0, 2 )
						shufps		xmm2, xmm1, R_SHUFFLEPS( 1, 3, 1, 3 )
						movaps		xmm3, [edi+32]
						addps		xmm0, xmm2
						mulps		xmm3, xmm7
						STORE4( 0, xmm0, xmm4 )
						shufps		xmm3, xmm3, R_SHUFFLEPS( 0, 2, 1, 3 )
						movhlps		xmm1, xmm3
						addps		xmm3, xmm1
						STORE2LO( 16, xmm3, xmm4 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numRows; i++ ) {
						dstPtr[i] STOREC mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1];
						mPtr += 2;
					}
					return;
				}
			}
			break;
		}
		case 3: {
			switch( numRows ) {
				case 3: {		// 3x3 * 3x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movss		xmm0, [esi]
						movss		xmm4, [edi]
						mulss		xmm4, xmm0
						movss		xmm1, [esi+4]
						movss		xmm5, [edi+4]
						mulss		xmm5, xmm1
						addss		xmm4, xmm5
						movss		xmm2, [esi+8]
						movss		xmm6, [edi+8]
						mulss		xmm6, xmm2
						addss		xmm4, xmm6
						movss		xmm3, [edi+12]
						mulss		xmm3, xmm0
						STORE1( 0, xmm4, xmm7 );
						movss		xmm5, [edi+12+4]
						mulss		xmm5, xmm1
						addss		xmm3, xmm5
						movss		xmm6, [edi+12+8]
						mulss		xmm6, xmm2
						addss		xmm3, xmm6
						mulss		xmm0, [edi+24]
						mulss		xmm1, [edi+24+4]
						STORE1( 4, xmm3, xmm7 );
						addss		xmm0, xmm1
						mulss		xmm2, [edi+24+8]
						addss		xmm0, xmm2
						STORE1( 8, xmm0, xmm7 );
					}
					return;
				}
				case 6: {		// 6x3 * 3x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movss		xmm5, [esi]
						shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )
						movss		xmm6, [esi+4]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						movss		xmm7, [esi+8]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
						movaps		xmm0, [edi]								// xmm0 = 0, 1, 2, 3
						movlps		xmm1, [edi+4*4]
						shufps		xmm1, xmm0, R_SHUFFLEPS( 0, 1, 1, 2 )	// xmm1 = 4, 5, 1, 2
						movlps		xmm2, [edi+6*4]
						movhps		xmm2, [edi+8*4]							// xmm2 = 6, 7, 8, 9
						shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 3, 0, 3 )	// xmm0 = 0, 3, 6, 9
						mulps		xmm0, xmm5
						movlps		xmm3, [edi+10*4]
						shufps		xmm2, xmm3, R_SHUFFLEPS( 1, 2, 0, 1 )	// xmm2 = 7, 8, 10, 11
						movaps		xmm3, xmm1
						shufps		xmm1, xmm2, R_SHUFFLEPS( 2, 0, 0, 2 )	// xmm1 = 1, 4, 7, 10
						mulps		xmm1, xmm6
						shufps		xmm3, xmm2, R_SHUFFLEPS( 3, 1, 1, 3 )	// xmm3 = 2, 5, 8, 11
						mulps		xmm3, xmm7
						addps		xmm0, xmm1
						addps		xmm0, xmm3
						STORE4( 0, xmm0, xmm4 )
						movss		xmm1, [edi+12*4]
						mulss		xmm1, xmm5
						movss		xmm2, [edi+13*4]
						mulss		xmm2, xmm6
						movss		xmm3, [edi+14*4]
						mulss		xmm3, xmm7
						addss		xmm1, xmm2
						addss		xmm1, xmm3
						STORE1( 16, xmm1, xmm4 )
						mulss		xmm5, [edi+15*4]
						mulss		xmm6, [edi+16*4]
						mulss		xmm7, [edi+17*4]
						addss		xmm5, xmm6
						addss		xmm5, xmm7
						STORE1( 20, xmm5, xmm4 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numRows; i++ ) {
						dstPtr[i] STOREC mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1] + mPtr[2] * vPtr[2];
						mPtr += 3;
					}
					return;
				}
			}
			break;
		}
		case 4: {
			switch( numRows ) {
				case 4: {		// 4x4 * 4x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm6, qword ptr [esi ]
						movlps		xmm0, qword ptr [edi ]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 1, 0, 1 )
						movhps		xmm0, qword ptr [edi+16]
						mulps		xmm0, xmm6
						movlps		xmm7, qword ptr [esi+ 8]
						movlps		xmm2, qword ptr [edi+ 8]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 1, 0, 1 )
						movhps		xmm2, qword ptr [edi+24]
						mulps		xmm2, xmm7
						movlps		xmm1, qword ptr [edi+32]
						movhps		xmm1, qword ptr [edi+48]
						mulps		xmm1, xmm6
						movlps		xmm3, qword ptr [edi+40]
						addps		xmm0, xmm2
						movhps		xmm3, qword ptr [edi+56]
						mulps		xmm3, xmm7
						movaps		xmm4, xmm0
						addps		xmm1, xmm3
						shufps		xmm4, xmm1, R_SHUFFLEPS( 0, 2, 0, 2 )
						shufps		xmm0, xmm1, R_SHUFFLEPS( 1, 3, 1, 3 )
						addps		xmm0, xmm4
						STORE4( 0, xmm0, xmm2 )
					}
					return;
				}
				case 6: {		// 6x4 * 4x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm6, qword ptr [esi+ 0]
						movlps		xmm0, qword ptr [edi+ 0]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 1, 0, 1 )
						movhps		xmm0, qword ptr [edi+16]
						mulps		xmm0, xmm6
						movlps		xmm7, qword ptr [esi+ 8]
						movlps		xmm2, qword ptr [edi+ 8]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 1, 0, 1 )
						movhps		xmm2, qword ptr [edi+24]
						mulps		xmm2, xmm7
						movlps		xmm1, qword ptr [edi+32]
						movhps		xmm1, qword ptr [edi+48]
						mulps		xmm1, xmm6
						movlps		xmm3, qword ptr [edi+40]
						addps		xmm0, xmm2
						movhps		xmm3, qword ptr [edi+56]
						mulps		xmm3, xmm7
						movaps		xmm4, xmm0
						addps		xmm1, xmm3
						shufps		xmm4, xmm1, R_SHUFFLEPS( 0, 2, 0, 2 )
						shufps		xmm0, xmm1, R_SHUFFLEPS( 1, 3, 1, 3 )
						addps		xmm0, xmm4
						movlps		xmm1, qword ptr [edi+64]
						movhps		xmm1, qword ptr [edi+80]
						STORE4( 0, xmm0, xmm4 )
						mulps		xmm1, xmm6
						movlps		xmm2, qword ptr [edi+72]
						movhps		xmm2, qword ptr [edi+88]
						mulps		xmm2, xmm7
						addps		xmm1, xmm2
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 2, 1, 3 )
						movhlps		xmm3, xmm1
						addps		xmm1, xmm3
						STORE2LO( 16, xmm1, xmm4 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numRows; i++ ) {
						dstPtr[i] STOREC mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1] + mPtr[2] * vPtr[2] + mPtr[3] * vPtr[3];
						mPtr += 4;
					}
					return;
				}
			}
			break;
		}
		case 5: {
			switch( numRows ) {
				case 5: {		// 5x5 * 5x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movss		xmm0, [edi+5*4]							// xmm0 =  5,  X,  X,  X
						movhps		xmm0, [edi+0*4]							// xmm0 =  5,  X,  0,  1
						movss		xmm5, [edi+15*4]						// xmm4 = 15,  X,  X,  X
						movhps		xmm5, [edi+10*4]						// xmm5 = 15,  X, 10, 11
						movaps		xmm1, xmm0								// xmm1 =  5,  X,  0,  1
						shufps		xmm0, xmm5, R_SHUFFLEPS( 2, 0, 2, 0 )	// xmm0 =  0,  5, 10, 15
						movlps		xmm1, [edi+6*4]							// xmm1 =  6,  7,  0,  1
						movlps		xmm5, [edi+16*4]						// xmm5 = 16, 17, 10, 11
						movaps		xmm2, xmm1								// xmm2 =  6,  7,  0,  1
						shufps		xmm1, xmm5, R_SHUFFLEPS( 3, 0, 3, 0 )	// xmm1 =  1,  6, 11, 16
						movhps		xmm2, [edi+2*4]							// xmm2 =  6,  7,  2,  3
						movhps		xmm5, [edi+12*4]						// xmm5 = 16, 17, 12, 13
						movaps		xmm3, xmm2								// xmm3 =  6,  7,  2,  3
						shufps		xmm2, xmm5, R_SHUFFLEPS( 2, 1, 2, 1 )	// xmm2 =  2,  7, 12, 17
						movlps		xmm3, [edi+8*4]							// xmm3 =  8,  9,  2,  3
						movlps		xmm5, [edi+18*4]						// xmm5 = 18, 19, 12, 13
						movss		xmm4, [edi+4*4]							// xmm4 =  4,  X,  X,  X
						movlhps		xmm4, xmm3								// xmm4 =  4,  X,  8,  9
						shufps		xmm3, xmm5, R_SHUFFLEPS( 3, 0, 3, 0 )	// xmm3 =  3,  8, 13, 18
						movhps		xmm5, [edi+14*4]						// xmm6 = 18, 19, 14, 15
						shufps		xmm4, xmm5, R_SHUFFLEPS( 0, 3, 2, 1 )	// xmm4 =  4,  9, 14, 19
						movss		xmm7, [esi+0*4]
						shufps		xmm7, xmm7, 0
						mulps		xmm0, xmm7
						movss		xmm5, [esi+1*4]
						shufps		xmm5, xmm5, 0
						mulps		xmm1, xmm5
						addps		xmm0, xmm1
						movss		xmm6, [esi+2*4]
						shufps		xmm6, xmm6, 0
						mulps		xmm2, xmm6
						addps		xmm0, xmm2
						movss		xmm1, [esi+3*4]
						shufps		xmm1, xmm1, 0
						mulps		xmm3, xmm1
						addps		xmm0, xmm3
						movss		xmm2, [esi+4*4]
						shufps		xmm2, xmm2, 0
						mulps		xmm4, xmm2
						addps		xmm0, xmm4
						mulss		xmm7, [edi+20*4]
						mulss		xmm5, [edi+21*4]
						addps		xmm7, xmm5
						mulss		xmm6, [edi+22*4]
						addps		xmm7, xmm6
						mulss		xmm1, [edi+23*4]
						addps		xmm7, xmm1
						mulss		xmm2, [edi+24*4]
						addps		xmm7, xmm2
						STORE4( 0, xmm0, xmm3 )
						STORE1( 16, xmm7, xmm4 )
					}
					return;
				}
				case 6: {		// 6x5 * 5x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm6, [esi]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 1, 0, 1 )
						movlps		xmm7, [esi+8]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 1, 0, 1 )
						movlps		xmm0, [edi]
						movhps		xmm3, [edi+8]
						movaps		xmm1, [edi+16]
						movlps		xmm2, [edi+32]
						shufps		xmm0, xmm1, R_SHUFFLEPS( 0, 1, 1, 2 )	// xmm0 = 0, 1, 5, 6
						shufps		xmm1, xmm2, R_SHUFFLEPS( 0, 3, 0, 1 )	// xmm1 = 4, 7, 8, 9
						shufps		xmm3, xmm1, R_SHUFFLEPS( 2, 3, 1, 2 )	// xmm3 = 2, 3, 7, 8
						mulps		xmm0, xmm6
						mulps		xmm3, xmm7
						movlps		xmm2, [edi+40]
						addps		xmm0, xmm3								// xmm0 + xmm1
						movhps		xmm5, [edi+40+8]
						movlps		xmm3, [edi+40+16]
						movhps		xmm3, [edi+40+24]
						movlps		xmm4, [edi+40+32]
						shufps		xmm2, xmm3, R_SHUFFLEPS( 0, 1, 1, 2 )	// xmm2 = 10, 11, 15, 16
						shufps		xmm3, xmm4, R_SHUFFLEPS( 0, 3, 0, 1 )	// xmm3 = 14, 17, 18, 19
						shufps		xmm5, xmm3, R_SHUFFLEPS( 2, 3, 1, 2 )	// xmm5 = 12, 13, 17, 18
						mulps		xmm2, xmm6
						mulps		xmm5, xmm7
						addps		xmm2, xmm5								// xmm2 + xmm3
						movss		xmm5, [esi+16]
						shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )
						movaps		xmm4, xmm0
						shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 2, 0, 2 )
						shufps		xmm4, xmm2, R_SHUFFLEPS( 1, 3, 1, 3 )
						shufps		xmm1, xmm3, R_SHUFFLEPS( 0, 3, 0, 3 )
						addps		xmm0, xmm4
						mulps		xmm1, xmm5
						addps		xmm0, xmm1
						STORE4( 0, xmm0, xmm2 )
						movlps		xmm4, [edi+80]
						movhps		xmm3, [edi+80+8]
						movaps		xmm1, [edi+80+16]
						movlps		xmm2, [edi+80+32]
						shufps		xmm4, xmm1, R_SHUFFLEPS( 0, 1, 1, 2 )	// xmm4 = 20, 21, 25, 26
						shufps		xmm1, xmm2, R_SHUFFLEPS( 0, 3, 0, 1 )	// xmm1 = 24, 27, 28, 29
						shufps		xmm3, xmm1, R_SHUFFLEPS( 2, 3, 1, 2 )	// xmm3 = 22, 23, 27, 28
						mulps		xmm4, xmm6
						mulps		xmm3, xmm7
						mulps		xmm1, xmm5
						addps		xmm4, xmm3								// xmm4 + xmm1
						shufps		xmm1, xmm4, R_SHUFFLEPS( 0, 3, 0, 2 )
						shufps		xmm4, xmm4, R_SHUFFLEPS( 1, 3, 0, 0 )
						addps		xmm4, xmm1
						shufps		xmm1, xmm1, R_SHUFFLEPS( 2, 3, 0, 1 )
						addps		xmm4, xmm1
						STORE2LO( 16, xmm4, xmm2 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numRows; i++ ) {
						dstPtr[i] STOREC mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1] + mPtr[2] * vPtr[2] + mPtr[3] * vPtr[3] + mPtr[4] * vPtr[4];
						mPtr += 5;
					}
					return;
				}
			}
			break;
		}
		case 6: {
			switch( numRows ) {
				case 1: {		// 1x6 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
                        mov			eax, dstPtr
						movss		xmm0, [esi]
						mulss		xmm0, [edi]
						movss		xmm1, [esi+4]
						mulss		xmm1, [edi+4]
						movss		xmm2, [esi+8]
						addss		xmm0, xmm1
						mulss		xmm2, [edi+8]
						movss		xmm3, [esi+12]
						addss		xmm0, xmm2
						mulss		xmm3, [edi+12]
						movss		xmm4, [esi+16]
						addss		xmm0, xmm3
						mulss		xmm4, [edi+16]
						movss		xmm5, [esi+20]
						addss		xmm0, xmm4
						mulss		xmm5, [edi+20]
						movss		xmm6, [esi+24]
						addss		xmm0, xmm5
						mulss		xmm6, [edi+24]
						addss		xmm0, xmm6
						STORE1( 0, xmm0, xmm7 )
					}
					return;
				}
				case 2: {		// 2x6 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						// load CVecXD
						movlps		xmm4, [esi]
						movhps		xmm4, [esi+8]
						movlps		xmm5, [esi+16]
						movlhps		xmm5, xmm4
						movhlps		xmm6, xmm4
						movlhps		xmm6, xmm5
						// row 0 and 1
						movaps		xmm0, [edi]
						movaps		xmm1, [edi+16]
						movaps		xmm2, [edi+32]
						mulps		xmm0, xmm4
						mulps		xmm1, xmm5
						mulps		xmm2, xmm6
						movhlps		xmm3, xmm0
						movlhps		xmm3, xmm2
						addps		xmm1, xmm3
						shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 1, 2, 3 )
						addps		xmm1, xmm0
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 2, 1, 3 )
						movhlps		xmm0, xmm1
						addps		xmm0, xmm1
						STORE2LO( 0, xmm0, xmm3 )
					}
					return;
				}
				case 3: {		// 3x6 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						// load CVecXD
						movlps		xmm4, [esi]
						movhps		xmm4, [esi+8]
						movlps		xmm5, [esi+16]
						movlhps		xmm5, xmm4
						movhlps		xmm6, xmm4
						movlhps		xmm6, xmm5
						// row 0 and 1
						movaps		xmm0, [edi]
						movaps		xmm1, [edi+16]
						movaps		xmm2, [edi+32]
						mulps		xmm0, xmm4
						mulps		xmm1, xmm5
						mulps		xmm2, xmm6
						movhlps		xmm3, xmm0
						movlhps		xmm3, xmm2
						addps		xmm1, xmm3
						shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 1, 2, 3 )
						addps		xmm1, xmm0
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 2, 1, 3 )
						movhlps		xmm0, xmm1
						addps		xmm0, xmm1
						STORE2LO( 0, xmm0, xmm3 )
						// row 2
						movaps		xmm0, [edi+48]
						movaps		xmm1, [edi+48+16]
						mulps		xmm0, xmm4
						mulps		xmm1, xmm5
						addps		xmm0, xmm1
						movhlps		xmm1, xmm0
						addps		xmm0, xmm1
						movaps		xmm1, xmm0
						shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 0, 0, 0 )
						addss		xmm0, xmm1
						STORE1( 8, xmm0, xmm3 )
					}
					return;
				}
				case 4: {		// 4x6 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						// load CVecXD
						movlps		xmm4, [esi]
						movhps		xmm4, [esi+8]
						movlps		xmm5, [esi+16]
						movlhps		xmm5, xmm4
						movhlps		xmm6, xmm4
						movlhps		xmm6, xmm5
						// row 0 and 1
						movaps		xmm0, [edi]
						movaps		xmm1, [edi+16]
						movaps		xmm2, [edi+32]
						mulps		xmm0, xmm4
						mulps		xmm1, xmm5
						mulps		xmm2, xmm6
						movhlps		xmm7, xmm0
						movlhps		xmm7, xmm2
						addps		xmm7, xmm1
						shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 1, 2, 3 )
						addps		xmm7, xmm0
						// row 2 and 3
						movaps		xmm0, [edi+48]
						movaps		xmm1, [edi+48+16]
						movaps		xmm2, [edi+48+32]
						mulps		xmm0, xmm4
						mulps		xmm1, xmm5
						mulps		xmm2, xmm6
						movhlps		xmm3, xmm0
						movlhps		xmm3, xmm2
						addps		xmm1, xmm3
						shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 1, 2, 3 )
						addps		xmm1, xmm0
						// last 4 additions for the first 4 rows and store result
						movaps		xmm0, xmm7
						shufps		xmm7, xmm1, R_SHUFFLEPS( 0, 2, 0, 2 )
						shufps		xmm0, xmm1, R_SHUFFLEPS( 1, 3, 1, 3 )
						addps		xmm0, xmm7
						STORE4( 0, xmm0, xmm4 )
					}
					return;
				}
				case 5: {		// 5x6 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						// load CVecXD
						movlps		xmm4, [esi]
						movhps		xmm4, [esi+8]
						movlps		xmm5, [esi+16]
						movlhps		xmm5, xmm4
						movhlps		xmm6, xmm4
						movlhps		xmm6, xmm5
						// row 0 and 1
						movaps		xmm0, [edi]
						movaps		xmm1, [edi+16]
						movaps		xmm2, [edi+32]
						mulps		xmm0, xmm4
						mulps		xmm1, xmm5
						mulps		xmm2, xmm6
						movhlps		xmm7, xmm0
						movlhps		xmm7, xmm2
						addps		xmm7, xmm1
						shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 1, 2, 3 )
						addps		xmm7, xmm0
						// row 2 and 3
						movaps		xmm0, [edi+48]
						movaps		xmm1, [edi+48+16]
						movaps		xmm2, [edi+48+32]
						mulps		xmm0, xmm4
						mulps		xmm1, xmm5
						mulps		xmm2, xmm6
						movhlps		xmm3, xmm0
						movlhps		xmm3, xmm2
						addps		xmm1, xmm3
						shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 1, 2, 3 )
						addps		xmm1, xmm0
						// last 4 additions for the first 4 rows and store result
						movaps		xmm0, xmm7
						shufps		xmm7, xmm1, R_SHUFFLEPS( 0, 2, 0, 2 )
						shufps		xmm0, xmm1, R_SHUFFLEPS( 1, 3, 1, 3 )
						addps		xmm0, xmm7
						STORE4( 0, xmm0, xmm3 )
						// row 5
						movaps		xmm0, [edi+96]
						movaps		xmm1, [edi+96+16]
						mulps		xmm0, xmm4
						mulps		xmm1, xmm5
						addps		xmm0, xmm1
						movhlps		xmm1, xmm0
						addps		xmm0, xmm1
						movaps		xmm1, xmm0
						shufps		xmm1, xmm1, 0x01
						addss		xmm0, xmm1
						STORE1( 16, xmm0, xmm3 )
					}
					return;
				}
				case 6: {		// 6x6 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm7, qword ptr [esi]
						movlps		xmm6, qword ptr [esi+8]
						shufps		xmm7, xmm7, 0x44
						shufps		xmm6, xmm6, 0x44
						movlps		xmm0, qword ptr [edi    ]
						movhps		xmm0, qword ptr [edi+ 24]
						mulps		xmm0, xmm7
						movlps		xmm3, qword ptr [edi+  8]
						movhps		xmm3, qword ptr [edi+ 32]
						mulps		xmm3, xmm6
						movlps		xmm1, qword ptr [edi+ 48]
						movhps		xmm1, qword ptr [edi+ 72]
						mulps		xmm1, xmm7
						movlps		xmm2, qword ptr [edi+ 96]
						movhps		xmm2, qword ptr [edi+120]
						mulps		xmm2, xmm7
						movlps		xmm4, qword ptr [edi+ 56]
						movhps		xmm4, qword ptr [edi+ 80]
						movlps		xmm5, qword ptr [edi+104]
						movhps		xmm5, qword ptr [edi+128]
						mulps		xmm4, xmm6
						movlps		xmm7, qword ptr [esi+16]
						addps		xmm0, xmm3
						shufps		xmm7, xmm7, 0x44
						mulps		xmm5, xmm6
						addps		xmm1, xmm4
						movlps		xmm3, qword ptr [edi+ 16]
						movhps		xmm3, qword ptr [edi+ 40]
						addps		xmm2, xmm5
						movlps		xmm4, qword ptr [edi+ 64]
						movhps		xmm4, qword ptr [edi+ 88]
						mulps		xmm3, xmm7
						movlps		xmm5, qword ptr [edi+112]
						movhps		xmm5, qword ptr [edi+136]
						addps		xmm0, xmm3
						mulps		xmm4, xmm7
						mulps		xmm5, xmm7
						addps		xmm1, xmm4
						addps		xmm2, xmm5
						movaps		xmm6, xmm0
						shufps		xmm0, xmm1, 0x88
						shufps		xmm6, xmm1, 0xDD
						movaps		xmm7, xmm2
						shufps		xmm7, xmm2, 0x88
						shufps		xmm2, xmm2, 0xDD
						addps		xmm0, xmm6
						addps		xmm2, xmm7
						STORE4( 0, xmm0, xmm3 )
						STORE2LO( 16, xmm2, xmm4 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numRows; i++ ) {
						dstPtr[i] STOREC mPtr[0] * vPtr[0] + mPtr[1] * vPtr[1] + mPtr[2] * vPtr[2] +
									mPtr[3] * vPtr[3] + mPtr[4] * vPtr[4] + mPtr[5] * vPtr[5];
						mPtr += 6;
					}
					return;
				}
			}
			break;
		}
		default: {
			int numColumns = mat.getNumColumns();
			for ( int i = 0; i < numRows; i++ ) {
				float sum = mPtr[0] * vPtr[0];
				for ( int j = 1; j < numColumns; j++ ) {
					sum += mPtr[j] * vPtr[j];
				}
				dstPtr[i] STOREC sum;
				mPtr += numColumns;
			}
			break;
		}
	}

#undef STOREC
#undef STORE4
#undef STORE2HI
#undef STORE2LO
#undef STORE1
}

/*
============
CSIMD_SSE::matX_TransposeMultiplyVecX

	optimizes the following matrix multiplications:

	Nx6 * Nx1
	6xN * 6x1

	with N in the range [1-6]
============
*/
void VPCALL CSIMD_SSE::matX_TransposeMultiplyVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec ) {
#define STORE1( offset, reg1, reg2 )		\
	__asm movss		[eax+offset], reg1
#define STORE2LO( offset, reg1, reg2 )		\
	__asm movlps	[eax+offset], reg1
#define STORE2HI( offset, reg1, reg2 )		\
	__asm movhps	[eax+offset], reg1
#define STORE4( offset, reg1, reg2 )		\
	__asm movlps	[eax+offset], reg1		\
	__asm movhps	[eax+offset+8], reg1
#define STOREC		=

	int numColumns;
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
			switch( numColumns ) {
				case 6: {		// 1x6 * 1x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movss		xmm0, [esi]
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
						movaps		xmm1, xmm0
						mulps		xmm0, [edi]
						mulps		xmm1, [edi+16]
						STORE4( 0, xmm0, xmm2 )
						STORE2LO( 16, xmm1, xmm3 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numColumns; i++ ) {
						dstPtr[i] STOREC *(mPtr) * vPtr[0];
						mPtr++;
					}
					return;
				}
			}
			break;
		case 2:
			switch( numColumns ) {
				case 6: {		// 2x6 * 2x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm0, [esi]
						movaps		xmm1, xmm0
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
						shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 1, 1, 1 )
						movaps		xmm2, [edi]
						mulps		xmm2, xmm0
						movlps		xmm3, [edi+24]
						movhps		xmm3, [edi+32]
						mulps		xmm3, xmm1
						addps		xmm2, xmm3
						shufps		xmm0, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
						movlps		xmm4, [edi+16]
						movhps		xmm4, [edi+40]
						mulps		xmm4, xmm0
						movhlps		xmm3, xmm4
						addps		xmm3, xmm4
						STORE4( 0, xmm2, xmm5 )
						STORE2LO( 16, xmm3, xmm6 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numColumns; i++ ) {
						dstPtr[i] STOREC *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1];
						mPtr++;
					}
					return;
				}
			}
			break;
		case 3:
			switch( numColumns ) {
				case 6: {		// 3x6 * 3x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm0, [esi+0*4]
						movss		xmm1, [esi+2*4]
						movlps		xmm3, [edi+(0*6+0)*4]
						movhps		xmm3, [edi+(0*6+2)*4]
						movaps		xmm4, xmm0
						shufps		xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm3, xmm4
						movlps		xmm5, [edi+(1*6+0)*4]
						movhps		xmm5, [edi+(1*6+2)*4]
						movaps		xmm6, xmm0
						shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 1, 1, 1 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						movlps		xmm4, [edi+(2*6+0)*4]
						movhps		xmm4, [edi+(2*6+2)*4]
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm4, xmm1
						addps		xmm3, xmm4
						STORE4( 0, xmm3, xmm7 )
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 1, 1 )
						movlps		xmm3, [edi+(0*6+4)*4]
						movhps		xmm3, [edi+(1*6+4)*4]
						mulps		xmm3, xmm0
						movhlps		xmm4, xmm3
						addps		xmm3, xmm4
						movlps		xmm5, [edi+(2*6+4)*4]
						mulps		xmm5, xmm1
						addps		xmm3, xmm5
						STORE2LO( 16, xmm3, xmm7 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numColumns; i++ ) {
						dstPtr[i] STOREC *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1] + *(mPtr+2*numColumns) * vPtr[2];
						mPtr++;
					}
					return;
				}
			}
			break;
		case 4:
			switch( numColumns ) {
				case 6: {		// 4x6 * 4x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm0, [esi+0*4]
						movlps		xmm1, [esi+2*4]
						movaps		xmm3, xmm0
						shufps		xmm3, xmm3, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm3, [edi+(0*6+0)*4]
						movlps		xmm5, [edi+(1*6+0)*4]
						movhps		xmm5, [edi+(1*6+2)*4]
						movaps		xmm6, xmm0
						shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 1, 1, 1 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						movlps		xmm4, [edi+(2*6+0)*4]
						movhps		xmm4, [edi+(2*6+2)*4]
						movaps		xmm6, xmm1
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm4, xmm6
						addps		xmm3, xmm4
						movlps		xmm5, [edi+(3*6+0)*4]
						movhps		xmm5, [edi+(3*6+2)*4]
						movaps		xmm6, xmm1
						shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 1, 1, 1 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						STORE4( 0, xmm3, xmm7 )
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 1, 1 )
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 1, 1 )
						movlps		xmm3, [edi+(0*6+4)*4]
						movhps		xmm3, [edi+(1*6+4)*4]
						mulps		xmm3, xmm0
						movlps		xmm4, [edi+(2*6+4)*4]
						movhps		xmm4, [edi+(3*6+4)*4]
						mulps		xmm4, xmm1
						addps		xmm3, xmm4
						movhlps		xmm4, xmm3
						addps		xmm3, xmm4
						STORE2LO( 16, xmm3, xmm7 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numColumns; i++ ) {
						dstPtr[i] STOREC *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1] + *(mPtr+2*numColumns) * vPtr[2] +
								*(mPtr+3*numColumns) * vPtr[3];
						mPtr++;
					}
					return;
				}
			}
			break;
		case 5:
			switch( numColumns ) {
				case 6: {		// 5x6 * 5x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm0, [esi+0*4]
						movlps		xmm1, [esi+2*4]
						movss		xmm2, [esi+4*4]
						movaps		xmm3, xmm0
						shufps		xmm3, xmm3, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm3, [edi+(0*6+0)*4]
						movlps		xmm5, [edi+(1*6+0)*4]
						movhps		xmm5, [edi+(1*6+2)*4]
						movaps		xmm6, xmm0
						shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 1, 1, 1 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						movaps		xmm6, xmm1
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm6, [edi+(2*6+0)*4]
						addps		xmm3, xmm6
						movlps		xmm5, [edi+(3*6+0)*4]
						movhps		xmm5, [edi+(3*6+2)*4]
						movaps		xmm6, xmm1
						shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 1, 1, 1 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						shufps		xmm2, xmm2, R_SHUFFLEPS( 0, 0, 0, 0 )
						movaps		xmm4, xmm2
						mulps		xmm4, [edi+(4*6+0)*4]
						addps		xmm3, xmm4
						STORE4( 0, xmm3, xmm7 )
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 1, 1 )
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 1, 1 )
						movlps		xmm3, [edi+(0*6+4)*4]
						movhps		xmm3, [edi+(1*6+4)*4]
						mulps		xmm3, xmm0
						movlps		xmm4, [edi+(2*6+4)*4]
						movhps		xmm4, [edi+(3*6+4)*4]
						mulps		xmm4, xmm1
						addps		xmm3, xmm4
						movhlps		xmm4, xmm3
						addps		xmm3, xmm4
						movlps		xmm5, [edi+(4*6+4)*4]
						mulps		xmm5, xmm2
						addps		xmm3, xmm5
						STORE2LO( 16, xmm3, xmm7 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numColumns; i++ ) {
						dstPtr[i] STOREC *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1] + *(mPtr+2*numColumns) * vPtr[2] +
								*(mPtr+3*numColumns) * vPtr[3] + *(mPtr+4*numColumns) * vPtr[4];
						mPtr++;
					}
					return;
				}
			}
			break;
		case 6:
			switch( numColumns ) {
				case 1: {		// 6x1 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm0, [esi]
						movhps		xmm0, [esi+8]
						movlps		xmm1, [esi+16]
						mulps		xmm0, [edi]
						mulps		xmm1, [edi+16]
						shufps		xmm1, xmm0, R_SHUFFLEPS( 0, 1, 3, 2 )
						addps		xmm0, xmm1
						movhlps		xmm2, xmm0
						addss		xmm2, xmm0
						shufps		xmm0, xmm0, R_SHUFFLEPS( 1, 0, 0, 0 )
						addss		xmm2, xmm0
						STORE1( 0, xmm2, xmm3 )
					}
					return;
				}
				case 2: {		// 6x2 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm0, [esi+0*4]
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 1, 1 )
						movaps		xmm6, [edi+0*4]
						mulps		xmm6, xmm0
						movlps		xmm1, [esi+2*4]
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 1, 1 )
						movaps		xmm7, [edi+4*4]
						mulps		xmm7, xmm1
						addps		xmm6, xmm7
						movlps		xmm2, [esi+4*4]
						shufps		xmm2, xmm2, R_SHUFFLEPS( 0, 0, 1, 1 )
						movaps		xmm7, [edi+8*4]
						mulps		xmm7, xmm2
						addps		xmm6, xmm7
						movhlps		xmm3, xmm6
						addps		xmm3, xmm6
						STORE2LO( 0, xmm3, xmm7 )
					}
					return;
				}
				case 3: {		// 6x3 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movss		xmm0, [edi+(0*3+2)*4]
						movhps		xmm0, [edi+(0*3+0)*4]
						shufps		xmm0, xmm0, R_SHUFFLEPS( 2, 1, 3, 0 )
						movss		xmm6, [esi+0*4]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm6, xmm0
						movss		xmm1, [edi+(1*3+0)*4]
						movhps		xmm1, [edi+(1*3+1)*4]
						movss		xmm7, [esi+1*4]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm1
						addps		xmm6, xmm7
						movss		xmm2, [edi+(2*3+2)*4]
						movhps		xmm2, [edi+(2*3+0)*4]
						shufps		xmm2, xmm2, R_SHUFFLEPS( 2, 1, 3, 0 )
						movss		xmm7, [esi+2*4]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm2
						addps		xmm6, xmm7
						movss		xmm3, [edi+(3*3+0)*4]
						movhps		xmm3, [edi+(3*3+1)*4]
						movss		xmm7, [esi+3*4]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm3
						addps		xmm6, xmm7
						movss		xmm4, [edi+(4*3+2)*4]
						movhps		xmm4, [edi+(4*3+0)*4]
						shufps		xmm4, xmm4, R_SHUFFLEPS( 2, 1, 3, 0 )
						movss		xmm7, [esi+4*4]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm4
						addps		xmm6, xmm7
						movss		xmm5, [edi+(5*3+0)*4]
						movhps		xmm5, [edi+(5*3+1)*4]
						movss		xmm7, [esi+5*4]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm5
						addps		xmm6, xmm7
						STORE1( 0, xmm6, xmm7 )
						STORE2HI( 4, xmm6, xmm7 )
					}
					return;
				}
				case 4: {		// 6x4 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm3, [edi+(0*4+0)*4]
						movhps		xmm3, [edi+(0*4+2)*4]
						movss		xmm4, [esi+0*4]
						shufps		xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm3, xmm4
						movlps		xmm5, [edi+(1*4+0)*4]
						movhps		xmm5, [edi+(1*4+2)*4]
						movss		xmm6, [esi+1*4]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						movlps		xmm4, [edi+(2*4+0)*4]
						movhps		xmm4, [edi+(2*4+2)*4]
						movss		xmm6, [esi+2*4]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm4, xmm6
						addps		xmm3, xmm4
						movlps		xmm5, [edi+(3*4+0)*4]
						movhps		xmm5, [edi+(3*4+2)*4]
						movss		xmm6, [esi+3*4]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						movlps		xmm4, [edi+(4*4+0)*4]
						movhps		xmm4, [edi+(4*4+2)*4]
						movss		xmm6, [esi+4*4]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm4, xmm6
						addps		xmm3, xmm4
						movlps		xmm5, [edi+(5*4+0)*4]
						movhps		xmm5, [edi+(5*4+2)*4]
						movss		xmm6, [esi+5*4]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						STORE4( 0, xmm3, xmm7 )
					}
					return;
				}
				case 5: {		// 6x5 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm6, [edi+(0*5+0)*4]
						movhps		xmm6, [edi+(0*5+2)*4]
						movss		xmm0, [esi+0*4]
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm6, xmm0
						movlps		xmm7, [edi+(1*5+0)*4]
						movhps		xmm7, [edi+(1*5+2)*4]
						movss		xmm1, [esi+1*4]
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm1
						addps		xmm6, xmm7
						movlps		xmm7, [edi+(2*5+0)*4]
						movhps		xmm7, [edi+(2*5+2)*4]
						movss		xmm2, [esi+2*4]
						shufps		xmm2, xmm2, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm2
						addps		xmm6, xmm7
						movlps		xmm7, [edi+(3*5+0)*4]
						movhps		xmm7, [edi+(3*5+2)*4]
						movss		xmm3, [esi+3*4]
						shufps		xmm3, xmm3, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm3
						addps		xmm6, xmm7
						movlps		xmm7, [edi+(4*5+0)*4]
						movhps		xmm7, [edi+(4*5+2)*4]
						movss		xmm4, [esi+4*4]
						shufps		xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm4
						addps		xmm6, xmm7
						movlps		xmm7, [edi+(5*5+0)*4]
						movhps		xmm7, [edi+(5*5+2)*4]
						movss		xmm5, [esi+5*4]
						shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm5
						addps		xmm6, xmm7
						STORE4( 0, xmm6, xmm7 )
						movss		xmm6, [edi+(0*5+4)*4]
						mulss		xmm6, xmm0
						movss		xmm7, [edi+(1*5+4)*4]
						mulss		xmm7, xmm1
						addss		xmm6, xmm7
						movss		xmm7, [edi+(2*5+4)*4]
						mulss		xmm7, xmm2
						addss		xmm6, xmm7
						movss		xmm7, [edi+(3*5+4)*4]
						mulss		xmm7, xmm3
						addss		xmm6, xmm7
						movss		xmm7, [edi+(4*5+4)*4]
						mulss		xmm7, xmm4
						addss		xmm6, xmm7
						movss		xmm7, [edi+(5*5+4)*4]
						mulss		xmm7, xmm5
						addss		xmm6, xmm7
						STORE1( 16, xmm6, xmm7 )
					}
					return;
				}
				case 6: {		// 6x6 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm0, [esi+0*4]
						movlps		xmm1, [esi+2*4]
						movlps		xmm2, [esi+4*4]
						movaps		xmm3, xmm0
						shufps		xmm3, xmm3, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm3, [edi+(0*6+0)*4]
						movlps		xmm5, [edi+(1*6+0)*4]
						movhps		xmm5, [edi+(1*6+2)*4]
						movaps		xmm6, xmm0
						shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 1, 1, 1 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						movaps		xmm6, xmm1
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm6, [edi+(2*6+0)*4]
						addps		xmm3, xmm6
						movaps		xmm6, xmm1
						shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 1, 1, 1 )
						movlps		xmm5, [edi+(3*6+0)*4]
						movhps		xmm5, [edi+(3*6+2)*4]
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						movaps		xmm6, xmm2
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm6, [edi+(4*6+0)*4]
						addps		xmm3, xmm6
						movaps		xmm6, xmm2
						shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 1, 1, 1 )
						movlps		xmm5, [edi+(5*6+0)*4]
						movhps		xmm5, [edi+(5*6+2)*4]
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						STORE4( 0, xmm3, xmm7 )
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 1, 1 )
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 1, 1 )
						shufps		xmm2, xmm2, R_SHUFFLEPS( 0, 0, 1, 1 )
						movlps		xmm3, [edi+(0*6+4)*4]
						movhps		xmm3, [edi+(1*6+4)*4]
						mulps		xmm3, xmm0
						movlps		xmm4, [edi+(2*6+4)*4]
						movhps		xmm4, [edi+(3*6+4)*4]
						mulps		xmm4, xmm1
						addps		xmm3, xmm4
						movlps		xmm5, [edi+(4*6+4)*4]
						movhps		xmm5, [edi+(5*6+4)*4]
						mulps		xmm5, xmm2
						addps		xmm3, xmm5
						movhlps		xmm4, xmm3
						addps		xmm3, xmm4
						STORE2LO( 16, xmm3, xmm7 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numColumns; i++ ) {
						dstPtr[i] STOREC *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1] + *(mPtr+2*numColumns) * vPtr[2] +
								*(mPtr+3*numColumns) * vPtr[3] + *(mPtr+4*numColumns) * vPtr[4] + *(mPtr+5*numColumns) * vPtr[5];
						mPtr++;
					}
					return;
				}
			}
			break;
		default:
			int numRows = mat.getNumRows();
			for ( int i = 0; i < numColumns; i++ ) {
				mPtr = mat.toFloatPtr() + i;
				float sum = mPtr[0] * vPtr[0];
				for ( int j = 1; j < numRows; j++ ) {
					mPtr += numColumns;
					sum += mPtr[0] * vPtr[j];
				}
				dstPtr[i] STOREC sum;
			}
			break;
	}

#undef STOREC
#undef STORE4
#undef STORE2HI
#undef STORE2LO
#undef STORE1
}

/*
============
CSIMD_SSE::matX_TransposeMultiplyAddVecX

	optimizes the following matrix multiplications:

	Nx6 * Nx1
	6xN * 6x1

	with N in the range [1-6]
============
*/
void VPCALL CSIMD_SSE::matX_TransposeMultiplyAddVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec ) {
#define STORE1( offset, reg1, reg2 )		\
	__asm movss		reg2, [eax+offset]		\
	__asm addss		reg2, reg1				\
	__asm movss		[eax+offset], reg2
#define STORE2LO( offset, reg1, reg2 )		\
	__asm movlps	reg2, [eax+offset]		\
	__asm addps		reg2, reg1				\
	__asm movlps	[eax+offset], reg2
#define STORE2HI( offset, reg1, reg2 )		\
	__asm movhps	reg2, [eax+offset]		\
	__asm addps		reg2, reg1				\
	__asm movhps	[eax+offset], reg2
#define STORE4( offset, reg1, reg2 )		\
	__asm movlps	reg2, [eax+offset]		\
	__asm movhps	reg2, [eax+offset+8]	\
	__asm addps		reg2, reg1				\
	__asm movlps	[eax+offset], reg2		\
	__asm movhps	[eax+offset+8], reg2
#define STOREC		+=

	int numColumns;
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
			switch( numColumns ) {
				case 6: {		// 1x6 * 1x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movss		xmm0, [esi]
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
						movaps		xmm1, xmm0
						mulps		xmm0, [edi]
						mulps		xmm1, [edi+16]
						STORE4( 0, xmm0, xmm2 )
						STORE2LO( 16, xmm1, xmm3 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numColumns; i++ ) {
						dstPtr[i] STOREC *(mPtr) * vPtr[0];
						mPtr++;
					}
					return;
				}
			}
			break;
		case 2:
			switch( numColumns ) {
				case 6: {		// 2x6 * 2x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm0, [esi]
						movaps		xmm1, xmm0
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
						shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 1, 1, 1 )
						movaps		xmm2, [edi]
						mulps		xmm2, xmm0
						movlps		xmm3, [edi+24]
						movhps		xmm3, [edi+32]
						mulps		xmm3, xmm1
						addps		xmm2, xmm3
						shufps		xmm0, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
						movlps		xmm4, [edi+16]
						movhps		xmm4, [edi+40]
						mulps		xmm4, xmm0
						movhlps		xmm3, xmm4
						addps		xmm3, xmm4
						STORE4( 0, xmm2, xmm5 )
						STORE2LO( 16, xmm3, xmm6 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numColumns; i++ ) {
						dstPtr[i] STOREC *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1];
						mPtr++;
					}
					return;
				}
			}
			break;
		case 3:
			switch( numColumns ) {
				case 6: {		// 3x6 * 3x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm0, [esi+0*4]
						movss		xmm1, [esi+2*4]
						movlps		xmm3, [edi+(0*6+0)*4]
						movhps		xmm3, [edi+(0*6+2)*4]
						movaps		xmm4, xmm0
						shufps		xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm3, xmm4
						movlps		xmm5, [edi+(1*6+0)*4]
						movhps		xmm5, [edi+(1*6+2)*4]
						movaps		xmm6, xmm0
						shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 1, 1, 1 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						movlps		xmm4, [edi+(2*6+0)*4]
						movhps		xmm4, [edi+(2*6+2)*4]
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm4, xmm1
						addps		xmm3, xmm4
						STORE4( 0, xmm3, xmm7 )
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 1, 1 )
						movlps		xmm3, [edi+(0*6+4)*4]
						movhps		xmm3, [edi+(1*6+4)*4]
						mulps		xmm3, xmm0
						movhlps		xmm4, xmm3
						addps		xmm3, xmm4
						movlps		xmm5, [edi+(2*6+4)*4]
						mulps		xmm5, xmm1
						addps		xmm3, xmm5
						STORE2LO( 16, xmm3, xmm7 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numColumns; i++ ) {
						dstPtr[i] STOREC *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1] + *(mPtr+2*numColumns) * vPtr[2];
						mPtr++;
					}
					return;
				}
			}
			break;
		case 4:
			switch( numColumns ) {
				case 6: {		// 4x6 * 4x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm0, [esi+0*4]
						movlps		xmm1, [esi+2*4]
						movaps		xmm3, xmm0
						shufps		xmm3, xmm3, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm3, [edi+(0*6+0)*4]
						movlps		xmm5, [edi+(1*6+0)*4]
						movhps		xmm5, [edi+(1*6+2)*4]
						movaps		xmm6, xmm0
						shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 1, 1, 1 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						movlps		xmm4, [edi+(2*6+0)*4]
						movhps		xmm4, [edi+(2*6+2)*4]
						movaps		xmm6, xmm1
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm4, xmm6
						addps		xmm3, xmm4
						movlps		xmm5, [edi+(3*6+0)*4]
						movhps		xmm5, [edi+(3*6+2)*4]
						movaps		xmm6, xmm1
						shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 1, 1, 1 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						STORE4( 0, xmm3, xmm7 )
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 1, 1 )
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 1, 1 )
						movlps		xmm3, [edi+(0*6+4)*4]
						movhps		xmm3, [edi+(1*6+4)*4]
						mulps		xmm3, xmm0
						movlps		xmm4, [edi+(2*6+4)*4]
						movhps		xmm4, [edi+(3*6+4)*4]
						mulps		xmm4, xmm1
						addps		xmm3, xmm4
						movhlps		xmm4, xmm3
						addps		xmm3, xmm4
						STORE2LO( 16, xmm3, xmm7 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numColumns; i++ ) {
						dstPtr[i] STOREC *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1] + *(mPtr+2*numColumns) * vPtr[2] +
								*(mPtr+3*numColumns) * vPtr[3];
						mPtr++;
					}
					return;
				}
			}
			break;
		case 5:
			switch( numColumns ) {
				case 6: {		// 5x6 * 5x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm0, [esi+0*4]
						movlps		xmm1, [esi+2*4]
						movss		xmm2, [esi+4*4]
						movaps		xmm3, xmm0
						shufps		xmm3, xmm3, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm3, [edi+(0*6+0)*4]
						movlps		xmm5, [edi+(1*6+0)*4]
						movhps		xmm5, [edi+(1*6+2)*4]
						movaps		xmm6, xmm0
						shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 1, 1, 1 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						movaps		xmm6, xmm1
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm6, [edi+(2*6+0)*4]
						addps		xmm3, xmm6
						movlps		xmm5, [edi+(3*6+0)*4]
						movhps		xmm5, [edi+(3*6+2)*4]
						movaps		xmm6, xmm1
						shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 1, 1, 1 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						shufps		xmm2, xmm2, R_SHUFFLEPS( 0, 0, 0, 0 )
						movaps		xmm4, xmm2
						mulps		xmm4, [edi+(4*6+0)*4]
						addps		xmm3, xmm4
						STORE4( 0, xmm3, xmm7 )
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 1, 1 )
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 1, 1 )
						movlps		xmm3, [edi+(0*6+4)*4]
						movhps		xmm3, [edi+(1*6+4)*4]
						mulps		xmm3, xmm0
						movlps		xmm4, [edi+(2*6+4)*4]
						movhps		xmm4, [edi+(3*6+4)*4]
						mulps		xmm4, xmm1
						addps		xmm3, xmm4
						movhlps		xmm4, xmm3
						addps		xmm3, xmm4
						movlps		xmm5, [edi+(4*6+4)*4]
						mulps		xmm5, xmm2
						addps		xmm3, xmm5
						STORE2LO( 16, xmm3, xmm7 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numColumns; i++ ) {
						dstPtr[i] STOREC *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1] + *(mPtr+2*numColumns) * vPtr[2] +
								*(mPtr+3*numColumns) * vPtr[3] + *(mPtr+4*numColumns) * vPtr[4];
						mPtr++;
					}
					return;
				}
			}
			break;
		case 6:
			switch( numColumns ) {
				case 1: {		// 6x1 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm0, [esi]
						movhps		xmm0, [esi+8]
						movlps		xmm1, [esi+16]
						mulps		xmm0, [edi]
						mulps		xmm1, [edi+16]
						shufps		xmm1, xmm0, R_SHUFFLEPS( 0, 1, 3, 2 )
						addps		xmm0, xmm1
						movhlps		xmm2, xmm0
						addss		xmm2, xmm0
						shufps		xmm0, xmm0, R_SHUFFLEPS( 1, 0, 0, 0 )
						addss		xmm2, xmm0
						STORE1( 0, xmm2, xmm3 )
					}
					return;
				}
				case 2: {		// 6x2 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm0, [esi+0*4]
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 1, 1 )
						movaps		xmm6, [edi+0*4]
						mulps		xmm6, xmm0
						movlps		xmm1, [esi+2*4]
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 1, 1 )
						movaps		xmm7, [edi+4*4]
						mulps		xmm7, xmm1
						addps		xmm6, xmm7
						movlps		xmm2, [esi+4*4]
						shufps		xmm2, xmm2, R_SHUFFLEPS( 0, 0, 1, 1 )
						movaps		xmm7, [edi+8*4]
						mulps		xmm7, xmm2
						addps		xmm6, xmm7
						movhlps		xmm3, xmm6
						addps		xmm3, xmm6
						STORE2LO( 0, xmm3, xmm7 )
					}
					return;
				}
				case 3: {		// 6x3 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movss		xmm0, [edi+(0*3+2)*4]
						movhps		xmm0, [edi+(0*3+0)*4]
						shufps		xmm0, xmm0, R_SHUFFLEPS( 2, 1, 3, 0 )
						movss		xmm6, [esi+0*4]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm6, xmm0
						movss		xmm1, [edi+(1*3+0)*4]
						movhps		xmm1, [edi+(1*3+1)*4]
						movss		xmm7, [esi+1*4]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm1
						addps		xmm6, xmm7
						movss		xmm2, [edi+(2*3+2)*4]
						movhps		xmm2, [edi+(2*3+0)*4]
						shufps		xmm2, xmm2, R_SHUFFLEPS( 2, 1, 3, 0 )
						movss		xmm7, [esi+2*4]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm2
						addps		xmm6, xmm7
						movss		xmm3, [edi+(3*3+0)*4]
						movhps		xmm3, [edi+(3*3+1)*4]
						movss		xmm7, [esi+3*4]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm3
						addps		xmm6, xmm7
						movss		xmm4, [edi+(4*3+2)*4]
						movhps		xmm4, [edi+(4*3+0)*4]
						shufps		xmm4, xmm4, R_SHUFFLEPS( 2, 1, 3, 0 )
						movss		xmm7, [esi+4*4]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm4
						addps		xmm6, xmm7
						movss		xmm5, [edi+(5*3+0)*4]
						movhps		xmm5, [edi+(5*3+1)*4]
						movss		xmm7, [esi+5*4]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm5
						addps		xmm6, xmm7
						STORE1( 0, xmm6, xmm7 )
						STORE2HI( 4, xmm6, xmm7 )
					}
					return;
				}
				case 4: {		// 6x4 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm3, [edi+(0*4+0)*4]
						movhps		xmm3, [edi+(0*4+2)*4]
						movss		xmm4, [esi+0*4]
						shufps		xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm3, xmm4
						movlps		xmm5, [edi+(1*4+0)*4]
						movhps		xmm5, [edi+(1*4+2)*4]
						movss		xmm6, [esi+1*4]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						movlps		xmm4, [edi+(2*4+0)*4]
						movhps		xmm4, [edi+(2*4+2)*4]
						movss		xmm6, [esi+2*4]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm4, xmm6
						addps		xmm3, xmm4
						movlps		xmm5, [edi+(3*4+0)*4]
						movhps		xmm5, [edi+(3*4+2)*4]
						movss		xmm6, [esi+3*4]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						movlps		xmm4, [edi+(4*4+0)*4]
						movhps		xmm4, [edi+(4*4+2)*4]
						movss		xmm6, [esi+4*4]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm4, xmm6
						addps		xmm3, xmm4
						movlps		xmm5, [edi+(5*4+0)*4]
						movhps		xmm5, [edi+(5*4+2)*4]
						movss		xmm6, [esi+5*4]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						STORE4( 0, xmm3, xmm7 )
					}
					return;
				}
				case 5: {		// 6x5 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm6, [edi+(0*5+0)*4]
						movhps		xmm6, [edi+(0*5+2)*4]
						movss		xmm0, [esi+0*4]
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm6, xmm0
						movlps		xmm7, [edi+(1*5+0)*4]
						movhps		xmm7, [edi+(1*5+2)*4]
						movss		xmm1, [esi+1*4]
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm1
						addps		xmm6, xmm7
						movlps		xmm7, [edi+(2*5+0)*4]
						movhps		xmm7, [edi+(2*5+2)*4]
						movss		xmm2, [esi+2*4]
						shufps		xmm2, xmm2, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm2
						addps		xmm6, xmm7
						movlps		xmm7, [edi+(3*5+0)*4]
						movhps		xmm7, [edi+(3*5+2)*4]
						movss		xmm3, [esi+3*4]
						shufps		xmm3, xmm3, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm3
						addps		xmm6, xmm7
						movlps		xmm7, [edi+(4*5+0)*4]
						movhps		xmm7, [edi+(4*5+2)*4]
						movss		xmm4, [esi+4*4]
						shufps		xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm4
						addps		xmm6, xmm7
						movlps		xmm7, [edi+(5*5+0)*4]
						movhps		xmm7, [edi+(5*5+2)*4]
						movss		xmm5, [esi+5*4]
						shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm5
						addps		xmm6, xmm7
						STORE4( 0, xmm6, xmm7 )
						movss		xmm6, [edi+(0*5+4)*4]
						mulss		xmm6, xmm0
						movss		xmm7, [edi+(1*5+4)*4]
						mulss		xmm7, xmm1
						addss		xmm6, xmm7
						movss		xmm7, [edi+(2*5+4)*4]
						mulss		xmm7, xmm2
						addss		xmm6, xmm7
						movss		xmm7, [edi+(3*5+4)*4]
						mulss		xmm7, xmm3
						addss		xmm6, xmm7
						movss		xmm7, [edi+(4*5+4)*4]
						mulss		xmm7, xmm4
						addss		xmm6, xmm7
						movss		xmm7, [edi+(5*5+4)*4]
						mulss		xmm7, xmm5
						addss		xmm6, xmm7
						STORE1( 16, xmm6, xmm7 )
					}
					return;
				}
				case 6: {		// 6x6 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm0, [esi+0*4]
						movlps		xmm1, [esi+2*4]
						movlps		xmm2, [esi+4*4]
						movaps		xmm3, xmm0
						shufps		xmm3, xmm3, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm3, [edi+(0*6+0)*4]
						movlps		xmm5, [edi+(1*6+0)*4]
						movhps		xmm5, [edi+(1*6+2)*4]
						movaps		xmm6, xmm0
						shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 1, 1, 1 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						movaps		xmm6, xmm1
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm6, [edi+(2*6+0)*4]
						addps		xmm3, xmm6
						movaps		xmm6, xmm1
						shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 1, 1, 1 )
						movlps		xmm5, [edi+(3*6+0)*4]
						movhps		xmm5, [edi+(3*6+2)*4]
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						movaps		xmm6, xmm2
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm6, [edi+(4*6+0)*4]
						addps		xmm3, xmm6
						movaps		xmm6, xmm2
						shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 1, 1, 1 )
						movlps		xmm5, [edi+(5*6+0)*4]
						movhps		xmm5, [edi+(5*6+2)*4]
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						STORE4( 0, xmm3, xmm7 )
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 1, 1 )
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 1, 1 )
						shufps		xmm2, xmm2, R_SHUFFLEPS( 0, 0, 1, 1 )
						movlps		xmm3, [edi+(0*6+4)*4]
						movhps		xmm3, [edi+(1*6+4)*4]
						mulps		xmm3, xmm0
						movlps		xmm4, [edi+(2*6+4)*4]
						movhps		xmm4, [edi+(3*6+4)*4]
						mulps		xmm4, xmm1
						addps		xmm3, xmm4
						movlps		xmm5, [edi+(4*6+4)*4]
						movhps		xmm5, [edi+(5*6+4)*4]
						mulps		xmm5, xmm2
						addps		xmm3, xmm5
						movhlps		xmm4, xmm3
						addps		xmm3, xmm4
						STORE2LO( 16, xmm3, xmm7 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numColumns; i++ ) {
						dstPtr[i] STOREC *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1] + *(mPtr+2*numColumns) * vPtr[2] +
								*(mPtr+3*numColumns) * vPtr[3] + *(mPtr+4*numColumns) * vPtr[4] + *(mPtr+5*numColumns) * vPtr[5];
						mPtr++;
					}
					return;
				}
			}
			break;
		default:
			int numRows = mat.getNumRows();
			for ( int i = 0; i < numColumns; i++ ) {
				mPtr = mat.toFloatPtr() + i;
				float sum = mPtr[0] * vPtr[0];
				for ( int j = 1; j < numRows; j++ ) {
					mPtr += numColumns;
					sum += mPtr[0] * vPtr[j];
				}
				dstPtr[i] STOREC sum;
			}
			break;
	}

#undef STOREC
#undef STORE4
#undef STORE2HI
#undef STORE2LO
#undef STORE1
}

/*
============
void CSIMD_SSE::matX_TransposeMultiplySubVecX

	optimizes the following matrix multiplications:

	Nx6 * Nx1
	6xN * 6x1

	with N in the range [1-6]
============
*/
void VPCALL CSIMD_SSE::matX_TransposeMultiplySubVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec ) {
#define STORE1( offset, reg1, reg2 )		\
	__asm movss		reg2, [eax+offset]		\
	__asm subss		reg2, reg1				\
	__asm movss		[eax+offset], reg2
#define STORE2LO( offset, reg1, reg2 )		\
	__asm movlps	reg2, [eax+offset]		\
	__asm subps		reg2, reg1				\
	__asm movlps	[eax+offset], reg2
#define STORE2HI( offset, reg1, reg2 )		\
	__asm movhps	reg2, [eax+offset]		\
	__asm subps		reg2, reg1				\
	__asm movhps	[eax+offset], reg2
#define STORE4( offset, reg1, reg2 )		\
	__asm movlps	reg2, [eax+offset]		\
	__asm movhps	reg2, [eax+offset+8]	\
	__asm subps		reg2, reg1				\
	__asm movlps	[eax+offset], reg2		\
	__asm movhps	[eax+offset+8], reg2
#define STOREC		-=

	int numColumns;
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
			switch( numColumns ) {
				case 6: {		// 1x6 * 1x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movss		xmm0, [esi]
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
						movaps		xmm1, xmm0
						mulps		xmm0, [edi]
						mulps		xmm1, [edi+16]
						STORE4( 0, xmm0, xmm2 )
						STORE2LO( 16, xmm1, xmm3 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numColumns; i++ ) {
						dstPtr[i] STOREC *(mPtr) * vPtr[0];
						mPtr++;
					}
					return;
				}
			}
			break;
		case 2:
			switch( numColumns ) {
				case 6: {		// 2x6 * 2x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm0, [esi]
						movaps		xmm1, xmm0
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
						shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 1, 1, 1 )
						movaps		xmm2, [edi]
						mulps		xmm2, xmm0
						movlps		xmm3, [edi+24]
						movhps		xmm3, [edi+32]
						mulps		xmm3, xmm1
						addps		xmm2, xmm3
						shufps		xmm0, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
						movlps		xmm4, [edi+16]
						movhps		xmm4, [edi+40]
						mulps		xmm4, xmm0
						movhlps		xmm3, xmm4
						addps		xmm3, xmm4
						STORE4( 0, xmm2, xmm5 )
						STORE2LO( 16, xmm3, xmm6 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numColumns; i++ ) {
						dstPtr[i] STOREC *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1];
						mPtr++;
					}
					return;
				}
			}
			break;
		case 3:
			switch( numColumns ) {
				case 6: {		// 3x6 * 3x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm0, [esi+0*4]
						movss		xmm1, [esi+2*4]
						movlps		xmm3, [edi+(0*6+0)*4]
						movhps		xmm3, [edi+(0*6+2)*4]
						movaps		xmm4, xmm0
						shufps		xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm3, xmm4
						movlps		xmm5, [edi+(1*6+0)*4]
						movhps		xmm5, [edi+(1*6+2)*4]
						movaps		xmm6, xmm0
						shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 1, 1, 1 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						movlps		xmm4, [edi+(2*6+0)*4]
						movhps		xmm4, [edi+(2*6+2)*4]
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm4, xmm1
						addps		xmm3, xmm4
						STORE4( 0, xmm3, xmm7 )
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 1, 1 )
						movlps		xmm3, [edi+(0*6+4)*4]
						movhps		xmm3, [edi+(1*6+4)*4]
						mulps		xmm3, xmm0
						movhlps		xmm4, xmm3
						addps		xmm3, xmm4
						movlps		xmm5, [edi+(2*6+4)*4]
						mulps		xmm5, xmm1
						addps		xmm3, xmm5
						STORE2LO( 16, xmm3, xmm7 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numColumns; i++ ) {
						dstPtr[i] STOREC *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1] + *(mPtr+2*numColumns) * vPtr[2];
						mPtr++;
					}
					return;
				}
			}
			break;
		case 4:
			switch( numColumns ) {
				case 6: {		// 4x6 * 4x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm0, [esi+0*4]
						movlps		xmm1, [esi+2*4]
						movaps		xmm3, xmm0
						shufps		xmm3, xmm3, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm3, [edi+(0*6+0)*4]
						movlps		xmm5, [edi+(1*6+0)*4]
						movhps		xmm5, [edi+(1*6+2)*4]
						movaps		xmm6, xmm0
						shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 1, 1, 1 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						movlps		xmm4, [edi+(2*6+0)*4]
						movhps		xmm4, [edi+(2*6+2)*4]
						movaps		xmm6, xmm1
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm4, xmm6
						addps		xmm3, xmm4
						movlps		xmm5, [edi+(3*6+0)*4]
						movhps		xmm5, [edi+(3*6+2)*4]
						movaps		xmm6, xmm1
						shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 1, 1, 1 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						STORE4( 0, xmm3, xmm7 )
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 1, 1 )
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 1, 1 )
						movlps		xmm3, [edi+(0*6+4)*4]
						movhps		xmm3, [edi+(1*6+4)*4]
						mulps		xmm3, xmm0
						movlps		xmm4, [edi+(2*6+4)*4]
						movhps		xmm4, [edi+(3*6+4)*4]
						mulps		xmm4, xmm1
						addps		xmm3, xmm4
						movhlps		xmm4, xmm3
						addps		xmm3, xmm4
						STORE2LO( 16, xmm3, xmm7 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numColumns; i++ ) {
						dstPtr[i] STOREC *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1] + *(mPtr+2*numColumns) * vPtr[2] +
								*(mPtr+3*numColumns) * vPtr[3];
						mPtr++;
					}
					return;
				}
			}
			break;
		case 5:
			switch( numColumns ) {
				case 6: {		// 5x6 * 5x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm0, [esi+0*4]
						movlps		xmm1, [esi+2*4]
						movss		xmm2, [esi+4*4]
						movaps		xmm3, xmm0
						shufps		xmm3, xmm3, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm3, [edi+(0*6+0)*4]
						movlps		xmm5, [edi+(1*6+0)*4]
						movhps		xmm5, [edi+(1*6+2)*4]
						movaps		xmm6, xmm0
						shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 1, 1, 1 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						movaps		xmm6, xmm1
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm6, [edi+(2*6+0)*4]
						addps		xmm3, xmm6
						movlps		xmm5, [edi+(3*6+0)*4]
						movhps		xmm5, [edi+(3*6+2)*4]
						movaps		xmm6, xmm1
						shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 1, 1, 1 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						shufps		xmm2, xmm2, R_SHUFFLEPS( 0, 0, 0, 0 )
						movaps		xmm4, xmm2
						mulps		xmm4, [edi+(4*6+0)*4]
						addps		xmm3, xmm4
						STORE4( 0, xmm3, xmm7 )
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 1, 1 )
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 1, 1 )
						movlps		xmm3, [edi+(0*6+4)*4]
						movhps		xmm3, [edi+(1*6+4)*4]
						mulps		xmm3, xmm0
						movlps		xmm4, [edi+(2*6+4)*4]
						movhps		xmm4, [edi+(3*6+4)*4]
						mulps		xmm4, xmm1
						addps		xmm3, xmm4
						movhlps		xmm4, xmm3
						addps		xmm3, xmm4
						movlps		xmm5, [edi+(4*6+4)*4]
						mulps		xmm5, xmm2
						addps		xmm3, xmm5
						STORE2LO( 16, xmm3, xmm7 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numColumns; i++ ) {
						dstPtr[i] STOREC *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1] + *(mPtr+2*numColumns) * vPtr[2] +
								*(mPtr+3*numColumns) * vPtr[3] + *(mPtr+4*numColumns) * vPtr[4];
						mPtr++;
					}
					return;
				}
			}
			break;
		case 6:
			switch( numColumns ) {
				case 1: {		// 6x1 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm0, [esi]
						movhps		xmm0, [esi+8]
						movlps		xmm1, [esi+16]
						mulps		xmm0, [edi]
						mulps		xmm1, [edi+16]
						shufps		xmm1, xmm0, R_SHUFFLEPS( 0, 1, 3, 2 )
						addps		xmm0, xmm1
						movhlps		xmm2, xmm0
						addss		xmm2, xmm0
						shufps		xmm0, xmm0, R_SHUFFLEPS( 1, 0, 0, 0 )
						addss		xmm2, xmm0
						STORE1( 0, xmm2, xmm3 )
					}
					return;
				}
				case 2: {		// 6x2 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm0, [esi+0*4]
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 1, 1 )
						movaps		xmm6, [edi+0*4]
						mulps		xmm6, xmm0
						movlps		xmm1, [esi+2*4]
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 1, 1 )
						movaps		xmm7, [edi+4*4]
						mulps		xmm7, xmm1
						addps		xmm6, xmm7
						movlps		xmm2, [esi+4*4]
						shufps		xmm2, xmm2, R_SHUFFLEPS( 0, 0, 1, 1 )
						movaps		xmm7, [edi+8*4]
						mulps		xmm7, xmm2
						addps		xmm6, xmm7
						movhlps		xmm3, xmm6
						addps		xmm3, xmm6
						STORE2LO( 0, xmm3, xmm7 )
					}
					return;
				}
				case 3: {		// 6x3 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movss		xmm0, [edi+(0*3+2)*4]
						movhps		xmm0, [edi+(0*3+0)*4]
						shufps		xmm0, xmm0, R_SHUFFLEPS( 2, 1, 3, 0 )
						movss		xmm6, [esi+0*4]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm6, xmm0
						movss		xmm1, [edi+(1*3+0)*4]
						movhps		xmm1, [edi+(1*3+1)*4]
						movss		xmm7, [esi+1*4]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm1
						addps		xmm6, xmm7
						movss		xmm2, [edi+(2*3+2)*4]
						movhps		xmm2, [edi+(2*3+0)*4]
						shufps		xmm2, xmm2, R_SHUFFLEPS( 2, 1, 3, 0 )
						movss		xmm7, [esi+2*4]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm2
						addps		xmm6, xmm7
						movss		xmm3, [edi+(3*3+0)*4]
						movhps		xmm3, [edi+(3*3+1)*4]
						movss		xmm7, [esi+3*4]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm3
						addps		xmm6, xmm7
						movss		xmm4, [edi+(4*3+2)*4]
						movhps		xmm4, [edi+(4*3+0)*4]
						shufps		xmm4, xmm4, R_SHUFFLEPS( 2, 1, 3, 0 )
						movss		xmm7, [esi+4*4]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm4
						addps		xmm6, xmm7
						movss		xmm5, [edi+(5*3+0)*4]
						movhps		xmm5, [edi+(5*3+1)*4]
						movss		xmm7, [esi+5*4]
						shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm5
						addps		xmm6, xmm7
						STORE1( 0, xmm6, xmm7 )
						STORE2HI( 4, xmm6, xmm7 )
					}
					return;
				}
				case 4: {		// 6x4 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm3, [edi+(0*4+0)*4]
						movhps		xmm3, [edi+(0*4+2)*4]
						movss		xmm4, [esi+0*4]
						shufps		xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm3, xmm4
						movlps		xmm5, [edi+(1*4+0)*4]
						movhps		xmm5, [edi+(1*4+2)*4]
						movss		xmm6, [esi+1*4]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						movlps		xmm4, [edi+(2*4+0)*4]
						movhps		xmm4, [edi+(2*4+2)*4]
						movss		xmm6, [esi+2*4]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm4, xmm6
						addps		xmm3, xmm4
						movlps		xmm5, [edi+(3*4+0)*4]
						movhps		xmm5, [edi+(3*4+2)*4]
						movss		xmm6, [esi+3*4]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						movlps		xmm4, [edi+(4*4+0)*4]
						movhps		xmm4, [edi+(4*4+2)*4]
						movss		xmm6, [esi+4*4]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm4, xmm6
						addps		xmm3, xmm4
						movlps		xmm5, [edi+(5*4+0)*4]
						movhps		xmm5, [edi+(5*4+2)*4]
						movss		xmm6, [esi+5*4]
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						STORE4( 0, xmm3, xmm7 )
					}
					return;
				}
				case 5: {		// 6x5 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm6, [edi+(0*5+0)*4]
						movhps		xmm6, [edi+(0*5+2)*4]
						movss		xmm0, [esi+0*4]
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm6, xmm0
						movlps		xmm7, [edi+(1*5+0)*4]
						movhps		xmm7, [edi+(1*5+2)*4]
						movss		xmm1, [esi+1*4]
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm1
						addps		xmm6, xmm7
						movlps		xmm7, [edi+(2*5+0)*4]
						movhps		xmm7, [edi+(2*5+2)*4]
						movss		xmm2, [esi+2*4]
						shufps		xmm2, xmm2, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm2
						addps		xmm6, xmm7
						movlps		xmm7, [edi+(3*5+0)*4]
						movhps		xmm7, [edi+(3*5+2)*4]
						movss		xmm3, [esi+3*4]
						shufps		xmm3, xmm3, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm3
						addps		xmm6, xmm7
						movlps		xmm7, [edi+(4*5+0)*4]
						movhps		xmm7, [edi+(4*5+2)*4]
						movss		xmm4, [esi+4*4]
						shufps		xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm4
						addps		xmm6, xmm7
						movlps		xmm7, [edi+(5*5+0)*4]
						movhps		xmm7, [edi+(5*5+2)*4]
						movss		xmm5, [esi+5*4]
						shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm7, xmm5
						addps		xmm6, xmm7
						STORE4( 0, xmm6, xmm7 )
						movss		xmm6, [edi+(0*5+4)*4]
						mulss		xmm6, xmm0
						movss		xmm7, [edi+(1*5+4)*4]
						mulss		xmm7, xmm1
						addss		xmm6, xmm7
						movss		xmm7, [edi+(2*5+4)*4]
						mulss		xmm7, xmm2
						addss		xmm6, xmm7
						movss		xmm7, [edi+(3*5+4)*4]
						mulss		xmm7, xmm3
						addss		xmm6, xmm7
						movss		xmm7, [edi+(4*5+4)*4]
						mulss		xmm7, xmm4
						addss		xmm6, xmm7
						movss		xmm7, [edi+(5*5+4)*4]
						mulss		xmm7, xmm5
						addss		xmm6, xmm7
						STORE1( 16, xmm6, xmm7 )
					}
					return;
				}
				case 6: {		// 6x6 * 6x1
					__asm {
						mov			esi, vPtr
						mov			edi, mPtr
						mov			eax, dstPtr
						movlps		xmm0, [esi+0*4]
						movlps		xmm1, [esi+2*4]
						movlps		xmm2, [esi+4*4]
						movaps		xmm3, xmm0
						shufps		xmm3, xmm3, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm3, [edi+(0*6+0)*4]
						movlps		xmm5, [edi+(1*6+0)*4]
						movhps		xmm5, [edi+(1*6+2)*4]
						movaps		xmm6, xmm0
						shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 1, 1, 1 )
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						movaps		xmm6, xmm1
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm6, [edi+(2*6+0)*4]
						addps		xmm3, xmm6
						movaps		xmm6, xmm1
						shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 1, 1, 1 )
						movlps		xmm5, [edi+(3*6+0)*4]
						movhps		xmm5, [edi+(3*6+2)*4]
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						movaps		xmm6, xmm2
						shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
						mulps		xmm6, [edi+(4*6+0)*4]
						addps		xmm3, xmm6
						movaps		xmm6, xmm2
						shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 1, 1, 1 )
						movlps		xmm5, [edi+(5*6+0)*4]
						movhps		xmm5, [edi+(5*6+2)*4]
						mulps		xmm5, xmm6
						addps		xmm3, xmm5
						STORE4( 0, xmm3, xmm7 )
						shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 1, 1 )
						shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 1, 1 )
						shufps		xmm2, xmm2, R_SHUFFLEPS( 0, 0, 1, 1 )
						movlps		xmm3, [edi+(0*6+4)*4]
						movhps		xmm3, [edi+(1*6+4)*4]
						mulps		xmm3, xmm0
						movlps		xmm4, [edi+(2*6+4)*4]
						movhps		xmm4, [edi+(3*6+4)*4]
						mulps		xmm4, xmm1
						addps		xmm3, xmm4
						movlps		xmm5, [edi+(4*6+4)*4]
						movhps		xmm5, [edi+(5*6+4)*4]
						mulps		xmm5, xmm2
						addps		xmm3, xmm5
						movhlps		xmm4, xmm3
						addps		xmm3, xmm4
						STORE2LO( 16, xmm3, xmm7 )
					}
					return;
				}
				default: {
					for ( int i = 0; i < numColumns; i++ ) {
						dstPtr[i] STOREC *(mPtr) * vPtr[0] + *(mPtr+numColumns) * vPtr[1] + *(mPtr+2*numColumns) * vPtr[2] +
								*(mPtr+3*numColumns) * vPtr[3] + *(mPtr+4*numColumns) * vPtr[4] + *(mPtr+5*numColumns) * vPtr[5];
						mPtr++;
					}
					return;
				}
			}
			break;
		default:
			int numRows = mat.getNumRows();
			for ( int i = 0; i < numColumns; i++ ) {
				mPtr = mat.toFloatPtr() + i;
				float sum = mPtr[0] * vPtr[0];
				for ( int j = 1; j < numRows; j++ ) {
					mPtr += numColumns;
					sum += mPtr[0] * vPtr[j];
				}
				dstPtr[i] STOREC sum;
			}
			break;
	}

#undef STOREC
#undef STORE4
#undef STORE2HI
#undef STORE2LO
#undef STORE1
}

/*
============
CSIMD_SSE::matX_MultiplyMatX

	optimizes the following matrix multiplications:

	NxN * Nx6
	6xN * Nx6
	Nx6 * 6xN
	6x6 * 6xN

	with N in the range [1-6].

	The hot cache clock cycle counts are generally better for the SIMD version than the
	FPU version.at times up to 40% less clock cycles on a P3. In practise however,
	the results are poor probably due to memory access.
============
*/
void VPCALL CSIMD_SSE::matX_MultiplyMatX( CMatXD &dst, const CMatXD &m1, const CMatXD &m2 ) {
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
	n = m1.getNumColumns();

	switch( n ) {
		case 1: {
			if ( !(l^6) ) {
				switch( k ) {
					case 1:	{			// 1x1 * 1x6, no precision loss compared to FPU version
						__asm {
							mov			esi, m2Ptr
							mov			edi, m1Ptr
							mov			eax, dstPtr
							movss		xmm0, [edi]
							shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
							movaps		xmm1, [esi]
							mulps		xmm1, xmm0
							movaps		[eax], xmm1
							movlps		xmm2, [esi+16]
							mulps		xmm2, xmm0
							movlps		[eax+16], xmm2
						}
						return;
					}
					case 6: {			// 6x1 * 1x6, no precision loss compared to FPU version
						__asm {
							mov			esi, m2Ptr
							mov			edi, m1Ptr
							mov			eax, dstPtr
							xorps		xmm1, xmm1
							movaps		xmm0, [edi]
							movlps		xmm1, [edi+16]
							movlhps		xmm1, xmm0
							movhlps		xmm2, xmm0
							movlhps		xmm2, xmm1
							// row 0 and 1
							movaps		xmm3, [esi]
							movaps		xmm4, xmm3
							shufps		xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )
							movaps		xmm5, xmm3
							shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 1, 1 )
							movaps		xmm6, xmm3
							shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 1, 1, 1 )
							mulps		xmm4, xmm0
							mulps		xmm5, xmm1
							mulps		xmm6, xmm2
							movaps		[eax], xmm4
							movaps		[eax+16], xmm5
							movaps		[eax+32], xmm6
							// row 2 and 3
							movaps		xmm4, xmm3
							shufps		xmm4, xmm4, R_SHUFFLEPS( 2, 2, 2, 2 )
							movaps		xmm5, xmm3
							shufps		xmm5, xmm5, R_SHUFFLEPS( 2, 2, 3, 3 )
							shufps		xmm3, xmm3, R_SHUFFLEPS( 3, 3, 3, 3 )
							mulps		xmm4, xmm0
							mulps		xmm5, xmm1
							mulps		xmm3, xmm2
							movaps		[eax+48], xmm4
							movaps		[eax+64], xmm5
							movaps		[eax+80], xmm3
							// row 4 and 5
							movlps		xmm3, [esi+16]
							movaps		xmm4, xmm3
							shufps		xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )
							movaps		xmm5, xmm3
							shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 1, 1 )
							shufps		xmm3, xmm3, R_SHUFFLEPS( 1, 1, 1, 1 )
							mulps		xmm4, xmm0
							mulps		xmm5, xmm1
							mulps		xmm3, xmm2
							movaps		[eax+96], xmm4
							movaps		[eax+112], xmm5
							movaps		[eax+128], xmm3
						}
						return;
					}
				}
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
			if ( !(l^6) ) {
				switch( k ) {
					case 2: {			// 2x2 * 2x6

						#define MUL_Nx2_2x6_INIT								\
						__asm mov		esi, m2Ptr								\
						__asm mov		edi, m1Ptr								\
						__asm mov		eax, dstPtr								\
						__asm movaps	xmm0, [esi]								\
						__asm movlps	xmm1, [esi+16]							\
						__asm movhps	xmm1, [esi+40]							\
						__asm movlps	xmm2, [esi+24]							\
						__asm movhps	xmm2, [esi+32]

						#define MUL_Nx2_2x6_ROW2( row )							\
						__asm movaps	xmm3, [edi+row*16]						\
						__asm movaps	xmm5, xmm0								\
						__asm movaps	xmm4, xmm3								\
						__asm shufps	xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm5, xmm4								\
						__asm movaps	xmm4, xmm3								\
						__asm movaps	xmm6, xmm2								\
						__asm shufps	xmm4, xmm4, R_SHUFFLEPS( 1, 1, 1, 1 )	\
						__asm mulps		xmm6, xmm4								\
						__asm addps		xmm5, xmm6								\
						__asm movaps	[eax+row*48], xmm5						\
						__asm movaps	xmm4, xmm3								\
						__asm shufps	xmm4, xmm4, R_SHUFFLEPS( 0, 0, 1, 1 )	\
						__asm movaps	xmm7, xmm1								\
						__asm mulps		xmm7, xmm4								\
						__asm movaps	xmm4, xmm3								\
						__asm movaps	xmm5, xmm0								\
						__asm shufps	xmm4, xmm4, R_SHUFFLEPS( 2, 2, 2, 2 )	\
						__asm mulps		xmm5, xmm4								\
						__asm movaps	xmm4, xmm3								\
						__asm movaps	xmm6, xmm2								\
						__asm shufps	xmm4, xmm4, R_SHUFFLEPS( 3, 3, 3, 3 )	\
						__asm mulps		xmm6, xmm4								\
						__asm addps		xmm5, xmm6								\
						__asm shufps	xmm3, xmm3, R_SHUFFLEPS( 2, 2, 3, 3 )	\
						__asm movaps	xmm6, xmm1								\
						__asm mulps		xmm6, xmm3								\
						__asm movaps	xmm4, xmm7								\
						__asm movlhps	xmm7, xmm6								\
						__asm movhlps	xmm6, xmm4								\
						__asm addps		xmm6, xmm7								\
						__asm movlps	[eax+row*48+16], xmm6					\
						__asm movlps	[eax+row*48+24], xmm5					\
						__asm movhps	[eax+row*48+32], xmm5					\
						__asm movhps	[eax+row*48+40], xmm6

						MUL_Nx2_2x6_INIT
						MUL_Nx2_2x6_ROW2( 0 )

						return;
					}
					case 6: {			// 6x2 * 2x6

						MUL_Nx2_2x6_INIT
						MUL_Nx2_2x6_ROW2( 0 )
						MUL_Nx2_2x6_ROW2( 1 )
						MUL_Nx2_2x6_ROW2( 2 )

						return;
					}
				}
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
			if ( !(l^6) ) {
				switch( k ) {
					case 3: {			// 3x3 * 3x6
						__asm {
							mov		esi, m2Ptr
							mov		edi, m1Ptr
							mov		eax, dstPtr
							movaps	xmm5, xmmword ptr [esi]
							movlps	xmm6, qword ptr [esi+24]
							movhps	xmm6, qword ptr [esi+32]
							movaps	xmm7, xmmword ptr [esi+48]
							movss	xmm0, dword ptr [edi]
							shufps	xmm0, xmm0, 0
							mulps	xmm0, xmm5
							movss	xmm1, dword ptr [edi+4]
							shufps	xmm1, xmm1, 0
							mulps	xmm1, xmm6
							movss	xmm2, dword ptr [edi+8]
							shufps	xmm2, xmm2, 0
							mulps	xmm2, xmm7
							addps	xmm0, xmm1
							addps	xmm0, xmm2
							movaps	xmmword ptr [eax], xmm0
							movss	xmm3, dword ptr [edi+12]
							shufps	xmm3, xmm3, 0
							mulps	xmm3, xmm5
							movss	xmm4, dword ptr [edi+16]
							shufps	xmm4, xmm4, 0
							mulps	xmm4, xmm6
							movss	xmm0, dword ptr [edi+20]
							shufps	xmm0, xmm0, 0
							mulps	xmm0, xmm7
							addps	xmm3, xmm4
							addps	xmm0, xmm3
							movlps	qword ptr [eax+24], xmm0
							movhps	qword ptr [eax+32], xmm0
							movss	xmm1, dword ptr [edi+24]
							shufps	xmm1, xmm1, 0
							mulps	xmm1, xmm5
							movss	xmm2, dword ptr [edi+28]
							shufps	xmm2, xmm2, 0
							mulps	xmm2, xmm6
							movss	xmm3, dword ptr [edi+32]
							shufps	xmm3, xmm3, 0
							mulps	xmm3, xmm7
							addps	xmm1, xmm2
							addps	xmm1, xmm3
							movaps	xmmword ptr [eax+48], xmm1
							movlps	xmm5, qword ptr [esi+16]
							movlps	xmm6, qword ptr [esi+40]
							movlps	xmm7, qword ptr [esi+64]
							shufps	xmm5, xmm5, 0x44
							shufps	xmm6, xmm6, 0x44
							shufps	xmm7, xmm7, 0x44
							movaps	xmm3, xmmword ptr [edi]
							movlps	xmm4, qword ptr [edi+16]
							movaps	xmm0, xmm3
							shufps	xmm0, xmm0, 0xF0
							mulps	xmm0, xmm5
							movaps	xmm1, xmm3
							shufps	xmm1, xmm4, 0x05
							mulps	xmm1, xmm6
							shufps	xmm3, xmm4, 0x5A
							mulps	xmm3, xmm7
							addps	xmm1, xmm0
							addps	xmm1, xmm3
							movlps	qword ptr [eax+16], xmm1
							movhps	qword ptr [eax+40], xmm1
							movss	xmm0, dword ptr [edi+24]
							shufps	xmm0, xmm0, 0
							mulps	xmm0, xmm5
							movss	xmm2, dword ptr [edi+28]
							shufps	xmm2, xmm2, 0
							mulps	xmm2, xmm6
							movss	xmm4, dword ptr [edi+32]
							shufps	xmm4, xmm4, 0
							mulps	xmm4, xmm7
							addps	xmm0, xmm2
							addps	xmm0, xmm4
							movlps	qword ptr [eax+64], xmm0
						}
						return;
					}
					case 6: {			// 6x3 * 3x6
						#define MUL_Nx3_3x6_FIRST4COLUMNS_INIT						\
						__asm mov			esi, m2Ptr								\
						__asm mov			edi, m1Ptr								\
						__asm mov			eax, dstPtr								\
						__asm movlps		xmm0, [esi+ 0*4]						\
						__asm movhps		xmm0, [esi+ 2*4]						\
						__asm movlps		xmm1, [esi+ 6*4]						\
						__asm movhps		xmm1, [esi+ 8*4]						\
						__asm movlps		xmm2, [esi+12*4]						\
						__asm movhps		xmm2, [esi+14*4]

						#define MUL_Nx3_3x6_FIRST4COLUMNS_ROW( row )				\
						__asm movss			xmm3, [edi+(row*3+0)*4]					\
						__asm shufps		xmm3, xmm3, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps			xmm3, xmm0								\
						__asm movss			xmm4, [edi+(row*3+1)*4]					\
						__asm shufps		xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps			xmm4, xmm1								\
						__asm addps			xmm3, xmm4								\
						__asm movss			xmm5, [edi+(row*3+2)*4]					\
						__asm shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps			xmm5, xmm2								\
						__asm addps			xmm3, xmm5								\
						__asm movlps		[eax+(row*6+0)*4], xmm3					\
						__asm movhps		[eax+(row*6+2)*4], xmm3

						#define MUL_Nx3_3x6_LAST2COLUMNS_ROW6						\
						__asm movlps		xmm0, [esi+ 4*4]						\
						__asm movlps		xmm1, [esi+10*4]						\
						__asm movlps		xmm2, [esi+16*4]						\
						__asm shufps		xmm0, xmm0, 0x44						\
						__asm shufps		xmm1, xmm1, 0x44						\
						__asm shufps		xmm2, xmm2, 0x44						\
						__asm movlps		xmm3, [edi+0*4]							\
						__asm movhps		xmm3, [edi+2*4]							\
						__asm movaps		xmm4, xmm3								\
						__asm movaps		xmm5, xmm3								\
						__asm shufps		xmm3, xmm3, 0xF0						\
						__asm mulps			xmm3, xmm0								\
						__asm movlps		xmm6, [edi+4*4]							\
						__asm movhps		xmm6, [edi+6*4]							\
						__asm shufps		xmm4, xmm6, 0x05						\
						__asm mulps			xmm4, xmm1								\
						__asm addps			xmm3, xmm4								\
						__asm shufps		xmm5, xmm6, 0x5A						\
						__asm mulps			xmm5, xmm2								\
						__asm addps			xmm3, xmm5								\
						__asm movlps		[eax+4*4], xmm3							\
						__asm movhps		[eax+10*4], xmm3						\
						__asm movaps		xmm5, xmm6								\
						__asm movlps		xmm3, [edi+8*4]							\
						__asm movhps		xmm3, [edi+10*4]						\
						__asm movaps		xmm4, xmm3								\
						__asm shufps		xmm5, xmm3, 0x5A						\
						__asm mulps			xmm5, xmm0								\
						__asm shufps		xmm6, xmm3, 0xAF						\
						__asm mulps			xmm6, xmm1								\
						__asm addps			xmm5, xmm6								\
						__asm shufps		xmm4, xmm4, 0xF0						\
						__asm mulps			xmm4, xmm2								\
						__asm addps			xmm4, xmm5								\
						__asm movlps		[eax+16*4], xmm4						\
						__asm movhps		[eax+22*4], xmm4						\
						__asm movlps		xmm6, [edi+12*4]						\
						__asm movhps		xmm6, [edi+14*4]						\
						__asm movaps		xmm5, xmm6								\
						__asm movaps		xmm4, xmm6								\
						__asm shufps		xmm6, xmm6, 0xF0						\
						__asm mulps			xmm6, xmm0								\
						__asm movlps		xmm3, [edi+16*4]						\
						__asm shufps		xmm5, xmm3, 0x05						\
						__asm mulps			xmm5, xmm1								\
						__asm addps			xmm5, xmm6								\
						__asm shufps		xmm4, xmm3, 0x5A						\
						__asm mulps			xmm4, xmm2								\
						__asm addps			xmm4, xmm5								\
						__asm movlps		[eax+28*4], xmm4						\
						__asm movhps		[eax+34*4], xmm4

						MUL_Nx3_3x6_FIRST4COLUMNS_INIT
						MUL_Nx3_3x6_FIRST4COLUMNS_ROW( 0 )
						MUL_Nx3_3x6_FIRST4COLUMNS_ROW( 1 )
						MUL_Nx3_3x6_FIRST4COLUMNS_ROW( 2 )
						MUL_Nx3_3x6_FIRST4COLUMNS_ROW( 3 )
						MUL_Nx3_3x6_FIRST4COLUMNS_ROW( 4 )
						MUL_Nx3_3x6_FIRST4COLUMNS_ROW( 5 )
						MUL_Nx3_3x6_LAST2COLUMNS_ROW6

						return;
					}
				}
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
			if ( !(l^6) ) {
				switch( k ) {
					case 4: {			// 4x4 * 4x6

						#define MUL_Nx4_4x6_FIRST4COLUMNS_INIT						\
						__asm mov			esi, m2Ptr								\
						__asm mov			edi, m1Ptr								\
						__asm mov			eax, dstPtr								\
						__asm movlps		xmm0, [esi+ 0*4]						\
						__asm movhps		xmm0, [esi+ 2*4]						\
						__asm movlps		xmm1, [esi+ 6*4]						\
						__asm movhps		xmm1, [esi+ 8*4]						\
						__asm movlps		xmm2, [esi+12*4]						\
						__asm movhps		xmm2, [esi+14*4]						\
						__asm movlps		xmm3, [esi+18*4]						\
						__asm movhps		xmm3, [esi+20*4]

						#define MUL_Nx4_4x6_FIRST4COLUMNS_ROW( row )				\
						__asm movss			xmm4, [edi+row*16+0*4]					\
						__asm shufps		xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps			xmm4, xmm0								\
						__asm movss			xmm5, [edi+row*16+1*4]					\
						__asm shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps			xmm5, xmm1								\
						__asm addps			xmm4, xmm5								\
						__asm movss			xmm6, [edi+row*16+2*4]					\
						__asm shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps			xmm6, xmm2								\
						__asm addps			xmm4, xmm6								\
						__asm movss			xmm7, [edi+row*16+3*4]					\
						__asm shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps			xmm7, xmm3								\
						__asm addps			xmm4, xmm7								\
						__asm movlps		[eax+row*24+0], xmm4					\
						__asm movhps		[eax+row*24+8], xmm4

						#define MUL_Nx4_4x6_LAST2COLUMNS_INIT						\
						__asm movlps		xmm0, [esi+ 4*4]						\
						__asm movlps		xmm1, [esi+10*4]						\
						__asm movlps		xmm2, [esi+16*4]						\
						__asm movlps		xmm3, [esi+22*4]						\
						__asm shufps		xmm0, xmm1, R_SHUFFLEPS( 0, 1, 0, 1 )	\
						__asm shufps		xmm1, xmm0, R_SHUFFLEPS( 0, 1, 0, 1 )	\
						__asm shufps		xmm2, xmm3, R_SHUFFLEPS( 0, 1, 0, 1 )	\
						__asm shufps		xmm3, xmm2, R_SHUFFLEPS( 0, 1, 0, 1 )

						#define MUL_Nx4_4x6_LAST2COLUMNS_ROW2( row )				\
						__asm movlps		xmm7, [edi+row*32+ 0*4]					\
						__asm movhps		xmm7, [edi+row*32+ 4*4]					\
						__asm movaps		xmm6, xmm7								\
						__asm shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 3, 3 )	\
						__asm mulps			xmm6, xmm0								\
						__asm shufps		xmm7, xmm7, R_SHUFFLEPS( 1, 1, 2, 2 )	\
						__asm mulps			xmm7, xmm1								\
						__asm addps			xmm6, xmm7								\
						__asm movlps		xmm4, [edi+row*32+ 2*4]					\
						__asm movhps		xmm4, [edi+row*32+ 6*4]					\
						__asm movaps		xmm5, xmm4								\
						__asm shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 3, 3 )	\
						__asm mulps			xmm5, xmm2								\
						__asm addps			xmm6, xmm5								\
						__asm shufps		xmm4, xmm4, R_SHUFFLEPS( 1, 1, 2, 2 )	\
						__asm mulps			xmm4, xmm3								\
						__asm addps			xmm6, xmm4								\
						__asm movlps		[eax+row*48+ 4*4], xmm6					\
						__asm movhps		[eax+row*48+10*4], xmm6

						MUL_Nx4_4x6_FIRST4COLUMNS_INIT
						MUL_Nx4_4x6_FIRST4COLUMNS_ROW( 0 )
						MUL_Nx4_4x6_FIRST4COLUMNS_ROW( 1 )
						MUL_Nx4_4x6_FIRST4COLUMNS_ROW( 2 )
						MUL_Nx4_4x6_FIRST4COLUMNS_ROW( 3 )
						MUL_Nx4_4x6_LAST2COLUMNS_INIT
						MUL_Nx4_4x6_LAST2COLUMNS_ROW2( 0 )
						MUL_Nx4_4x6_LAST2COLUMNS_ROW2( 1 )

						return;
					}
					case 6: {			// 6x4 * 4x6

						MUL_Nx4_4x6_FIRST4COLUMNS_INIT
						MUL_Nx4_4x6_FIRST4COLUMNS_ROW( 0 )
						MUL_Nx4_4x6_FIRST4COLUMNS_ROW( 1 )
						MUL_Nx4_4x6_FIRST4COLUMNS_ROW( 2 )
						MUL_Nx4_4x6_FIRST4COLUMNS_ROW( 3 )
						MUL_Nx4_4x6_FIRST4COLUMNS_ROW( 4 )
						MUL_Nx4_4x6_FIRST4COLUMNS_ROW( 5 )
						MUL_Nx4_4x6_LAST2COLUMNS_INIT
						MUL_Nx4_4x6_LAST2COLUMNS_ROW2( 0 )
						MUL_Nx4_4x6_LAST2COLUMNS_ROW2( 1 )
						MUL_Nx4_4x6_LAST2COLUMNS_ROW2( 2 )

						return;
					}
				}
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
			if ( !(l^6) ) {
				switch( k ) {
					case 5: {			// 5x5 * 5x6

						#define MUL_Nx5_5x6_FIRST4COLUMNS_INIT						\
						__asm mov			esi, m2Ptr								\
						__asm mov			edi, m1Ptr								\
						__asm mov			eax, dstPtr								\
						__asm movlps		xmm0, [esi+ 0*4]						\
						__asm movhps		xmm0, [esi+ 2*4]						\
						__asm movlps		xmm1, [esi+ 6*4]						\
						__asm movhps		xmm1, [esi+ 8*4]						\
						__asm movlps		xmm2, [esi+12*4]						\
						__asm movhps		xmm2, [esi+14*4]						\
						__asm movlps		xmm3, [esi+18*4]						\
						__asm movhps		xmm3, [esi+20*4]						\
						__asm movlps		xmm4, [esi+24*4]						\
						__asm movhps		xmm4, [esi+26*4]

						#define MUL_Nx5_5x6_FIRST4COLUMNS_ROW( row )				\
						__asm movss			xmm6, [edi+row*20+0*4]					\
						__asm shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps			xmm6, xmm0								\
						__asm movss			xmm5, [edi+row*20+1*4]					\
						__asm shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps			xmm5, xmm1								\
						__asm addps			xmm6, xmm5								\
						__asm movss			xmm5, [edi+row*20+2*4]					\
						__asm shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps			xmm5, xmm2								\
						__asm addps			xmm6, xmm5								\
						__asm movss			xmm5, [edi+row*20+3*4]					\
						__asm shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps			xmm5, xmm3								\
						__asm addps			xmm6, xmm5								\
						__asm movss			xmm5, [edi+row*20+4*4]					\
						__asm shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps			xmm5, xmm4								\
						__asm addps			xmm6, xmm5								\
						__asm movlps		[eax+row*24+0], xmm6					\
						__asm movhps		[eax+row*24+8], xmm6

						#define MUL_Nx5_5x6_LAST2COLUMNS_INIT						\
						__asm movlps		xmm0, [esi+ 4*4]						\
						__asm movlps		xmm1, [esi+10*4]						\
						__asm movlps		xmm2, [esi+16*4]						\
						__asm movlps		xmm3, [esi+22*4]						\
						__asm movlps		xmm4, [esi+28*4]						\
						__asm shufps		xmm0, xmm1, R_SHUFFLEPS( 0, 1, 0, 1 )	\
						__asm shufps		xmm1, xmm2, R_SHUFFLEPS( 0, 1, 0, 1 )	\
						__asm shufps		xmm2, xmm3, R_SHUFFLEPS( 0, 1, 0, 1 )	\
						__asm shufps		xmm3, xmm4, R_SHUFFLEPS( 0, 1, 0, 1 )	\
						__asm shufps		xmm4, xmm0, R_SHUFFLEPS( 0, 1, 0, 1 )

						#define MUL_Nx5_5x6_LAST2COLUMNS_ROW2( row )				\
						__asm movlps		xmm7, [edi+row*40+ 0*4]					\
						__asm movhps		xmm7, [edi+row*40+ 6*4]					\
						__asm movaps		xmm6, xmm7								\
						__asm shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 2, 2 )	\
						__asm mulps			xmm6, xmm0								\
						__asm movaps		xmm5, xmm7								\
						__asm shufps		xmm5, xmm5, R_SHUFFLEPS( 1, 1, 3, 3 )	\
						__asm mulps			xmm5, xmm1								\
						__asm addps			xmm6, xmm5								\
						__asm movlps		xmm7, [edi+row*40+ 2*4]					\
						__asm movhps		xmm7, [edi+row*40+ 8*4]					\
						__asm movaps		xmm5, xmm7								\
						__asm shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 2, 2 )	\
						__asm mulps			xmm5, xmm2								\
						__asm addps			xmm6, xmm5								\
						__asm movaps		xmm5, xmm7								\
						__asm shufps		xmm5, xmm5, R_SHUFFLEPS( 1, 1, 3, 3 )	\
						__asm mulps			xmm5, xmm3								\
						__asm addps			xmm6, xmm5								\
						__asm movlps		xmm5, [edi+row*40+ 4*4]					\
						__asm shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 1, 1 )	\
						__asm mulps			xmm5, xmm4								\
						__asm addps			xmm6, xmm5								\
						__asm movlps		[eax+row*48+ 4*4], xmm6					\
						__asm movhps		[eax+row*48+10*4], xmm6

						#define MUL_Nx5_5x6_LAST2COLUMNS_ROW( row )					\
						__asm movlps		xmm6, [edi+20*4+0*4]					\
						__asm unpcklps		xmm6, xmm6								\
						__asm mulps			xmm6, xmm0								\
						__asm movlps		xmm5, [edi+20*4+2*4]					\
						__asm unpcklps		xmm5, xmm5								\
						__asm mulps			xmm5, xmm2								\
						__asm addps			xmm6, xmm5								\
						__asm movss			xmm5, [edi+20*4+4*4]					\
						__asm unpcklps		xmm5, xmm5								\
						__asm mulps			xmm5, xmm4								\
						__asm addps			xmm6, xmm5								\
						__asm movhlps		xmm7, xmm6								\
						__asm addps			xmm6, xmm7								\
						__asm movlps		[eax+row*24+4*4], xmm6

						MUL_Nx5_5x6_FIRST4COLUMNS_INIT
						MUL_Nx5_5x6_FIRST4COLUMNS_ROW( 0 )
						MUL_Nx5_5x6_FIRST4COLUMNS_ROW( 1 )
						MUL_Nx5_5x6_FIRST4COLUMNS_ROW( 2 )
						MUL_Nx5_5x6_FIRST4COLUMNS_ROW( 3 )
						MUL_Nx5_5x6_FIRST4COLUMNS_ROW( 4 )
						MUL_Nx5_5x6_LAST2COLUMNS_INIT
						MUL_Nx5_5x6_LAST2COLUMNS_ROW2( 0 )
						MUL_Nx5_5x6_LAST2COLUMNS_ROW2( 1 )
						MUL_Nx5_5x6_LAST2COLUMNS_ROW( 4 )

						return;
					}
					case 6: {			// 6x5 * 5x6

						MUL_Nx5_5x6_FIRST4COLUMNS_INIT
						MUL_Nx5_5x6_FIRST4COLUMNS_ROW( 0 )
						MUL_Nx5_5x6_FIRST4COLUMNS_ROW( 1 )
						MUL_Nx5_5x6_FIRST4COLUMNS_ROW( 2 )
						MUL_Nx5_5x6_FIRST4COLUMNS_ROW( 3 )
						MUL_Nx5_5x6_FIRST4COLUMNS_ROW( 4 )
						MUL_Nx5_5x6_FIRST4COLUMNS_ROW( 5 )
						MUL_Nx5_5x6_LAST2COLUMNS_INIT
						MUL_Nx5_5x6_LAST2COLUMNS_ROW2( 0 )
						MUL_Nx5_5x6_LAST2COLUMNS_ROW2( 1 )
						MUL_Nx5_5x6_LAST2COLUMNS_ROW2( 2 )

						return;
					}
				}
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
					if ( !(l^1) ) {		// 1x6 * 6x1
						dstPtr[0] = m1Ptr[0] * m2Ptr[0] + m1Ptr[1] * m2Ptr[1] + m1Ptr[2] * m2Ptr[2] +
									 m1Ptr[3] * m2Ptr[3] + m1Ptr[4] * m2Ptr[4] + m1Ptr[5] * m2Ptr[5];
						return;
					}
					break;
				}
				case 2: {
					if ( !(l^2) ) {		// 2x6 * 6x2

						#define MUL_Nx6_6x2_INIT								\
						__asm mov		esi, m2Ptr								\
						__asm mov		edi, m1Ptr								\
						__asm mov		eax, dstPtr								\
						__asm movaps	xmm0, [esi]								\
						__asm movaps	xmm1, [esi+16]							\
						__asm movaps	xmm2, [esi+32]

						#define MUL_Nx6_6x2_ROW2( row )							\
						__asm movaps	xmm7, [edi+row*48+0*4]					\
						__asm movaps	xmm6, xmm7								\
						__asm shufps	xmm7, xmm7, R_SHUFFLEPS( 0, 0, 1, 1 )	\
						__asm mulps		xmm7, xmm0								\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 2, 2, 3, 3 )	\
						__asm mulps		xmm6, xmm1								\
						__asm addps		xmm7, xmm6								\
						__asm movaps	xmm6, [edi+row*48+4*4]					\
						__asm movaps	xmm5, xmm6								\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 1, 1 )	\
						__asm mulps		xmm6, xmm2								\
						__asm addps		xmm7, xmm6								\
						__asm shufps	xmm5, xmm5, R_SHUFFLEPS( 2, 2, 3, 3 )	\
						__asm mulps		xmm5, xmm0								\
						__asm movaps	xmm6, [edi+row*48+24+2*4]				\
						__asm movaps	xmm4, xmm6								\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 1, 1 )	\
						__asm mulps		xmm6, xmm1								\
						__asm addps		xmm5, xmm6								\
						__asm shufps	xmm4, xmm4, R_SHUFFLEPS( 2, 2, 3, 3 )	\
						__asm mulps		xmm4, xmm2								\
						__asm addps		xmm5, xmm4								\
						__asm movaps	xmm4, xmm5								\
						__asm movhlps	xmm5, xmm7								\
						__asm movlhps	xmm7, xmm4								\
						__asm addps		xmm7, xmm5								\
						__asm movaps	[eax+row*16], xmm7

						MUL_Nx6_6x2_INIT
						MUL_Nx6_6x2_ROW2( 0 )

						return;
					}
					break;
				}
				case 3: {
					if ( !(l^3) ) {		// 3x6 * 6x3

						#define MUL_Nx6_6x3_INIT								\
						__asm mov		esi, m2Ptr								\
						__asm mov		edi, m1Ptr								\
						__asm mov		eax, dstPtr								\
						__asm movss		xmm0, [esi+ 0*4]						\
						__asm movhps	xmm0, [esi+ 1*4]						\
						__asm movss		xmm1, [esi+ 3*4]						\
						__asm movhps	xmm1, [esi+ 4*4]						\
						__asm movss		xmm2, [esi+ 6*4]						\
						__asm movhps	xmm2, [esi+ 7*4]						\
						__asm movss		xmm3, [esi+ 9*4]						\
						__asm movhps	xmm3, [esi+10*4]						\
						__asm movss		xmm4, [esi+12*4]						\
						__asm movhps	xmm4, [esi+13*4]						\
						__asm movss		xmm5, [esi+15*4]						\
						__asm movhps	xmm5, [esi+16*4]

						#define MUL_Nx6_6x3_ROW( row )							\
						__asm movss		xmm7, [edi+row*24+0]					\
						__asm shufps	xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm7, xmm0								\
						__asm movss		xmm6, [edi+row*24+4]					\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm6, xmm1								\
						__asm addps		xmm7, xmm6								\
						__asm movss		xmm6, [edi+row*24+8]					\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm6, xmm2								\
						__asm addps		xmm7, xmm6								\
						__asm movss		xmm6, [edi+row*24+12]					\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm6, xmm3								\
						__asm addps		xmm7, xmm6								\
						__asm movss		xmm6, [edi+row*24+16]					\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm6, xmm4								\
						__asm addps		xmm7, xmm6								\
						__asm movss		xmm6, [edi+row*24+20]					\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm6, xmm5								\
						__asm addps		xmm7, xmm6								\
						__asm movss		[eax+row*12+0], xmm7					\
						__asm movhps	[eax+row*12+4], xmm7

						MUL_Nx6_6x3_INIT
						MUL_Nx6_6x3_ROW( 0 )
						MUL_Nx6_6x3_ROW( 1 )
						MUL_Nx6_6x3_ROW( 2 )

						return;
					}
					break;
				}
				case 4: {
					if ( !(l^4) ) {		// 4x6 * 6x4

						#define MUL_Nx6_6x4_INIT								\
						__asm mov		esi, m2Ptr								\
						__asm mov		edi, m1Ptr								\
						__asm mov		eax, dstPtr								\
						__asm movaps	xmm0, [esi]								\
						__asm movaps	xmm1, [esi+16]							\
						__asm movaps	xmm2, [esi+32]							\
						__asm movaps	xmm3, [esi+48]							\
						__asm movaps	xmm4, [esi+64]							\
						__asm movaps	xmm5, [esi+80]

						#define MUL_Nx6_6x4_ROW( row )							\
						__asm movss		xmm7, [edi+row*24+0]					\
						__asm shufps	xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm7, xmm0								\
						__asm movss		xmm6, [edi+row*24+4]					\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm6, xmm1								\
						__asm addps		xmm7, xmm6								\
						__asm movss		xmm6, [edi+row*24+8]					\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm6, xmm2								\
						__asm addps		xmm7, xmm6								\
						__asm movss		xmm6, [edi+row*24+12]					\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm6, xmm3								\
						__asm addps		xmm7, xmm6								\
						__asm movss		xmm6, [edi+row*24+16]					\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm6, xmm4								\
						__asm addps		xmm7, xmm6								\
						__asm movss		xmm6, [edi+row*24+20]					\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm6, xmm5								\
						__asm addps		xmm7, xmm6								\
						__asm movaps	[eax+row*16], xmm7

						MUL_Nx6_6x4_INIT
						MUL_Nx6_6x4_ROW( 0 )
						MUL_Nx6_6x4_ROW( 1 )
						MUL_Nx6_6x4_ROW( 2 )
						MUL_Nx6_6x4_ROW( 3 )

						return;
					}
					break;
				}
				case 5: {
					if ( !(l^5) ) {		// 5x6 * 6x5

						#define MUL_Nx6_6x5_INIT								\
						__asm mov		esi, m2Ptr								\
						__asm mov		edi, m1Ptr								\
						__asm mov		eax, dstPtr								\
						__asm movaps	xmm0, [esi]								\
						__asm movlps	xmm1, [esi+20]							\
						__asm movhps	xmm1, [esi+28]							\
						__asm movlps	xmm2, [esi+40]							\
						__asm movhps	xmm2, [esi+48]							\
						__asm movlps	xmm3, [esi+60]							\
						__asm movhps	xmm3, [esi+68]							\
						__asm movaps	xmm4, [esi+80]							\
						__asm movlps	xmm5, [esi+100]							\
						__asm movhps	xmm5, [esi+108]

						#define MUL_Nx6_6x5_ROW( row )							\
						__asm movss		xmm7, [edi+row*24+0]					\
						__asm shufps	xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm7, xmm0								\
						__asm fld		dword ptr [edi+(row*6+0)*4]				\
						__asm fmul		dword ptr [esi+(4+0*5)*4]				\
						__asm movss		xmm6, [edi+row*24+4]					\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm6, xmm1								\
						__asm addps		xmm7, xmm6								\
						__asm fld		dword ptr [edi+(row*6+1)*4]				\
						__asm fmul		dword ptr [esi+(4+1*5)*4]				\
						__asm faddp		st(1),st								\
						__asm movss		xmm6, [edi+row*24+8]					\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm6, xmm2								\
						__asm addps		xmm7, xmm6								\
						__asm fld		dword ptr [edi+(row*6+2)*4]				\
						__asm fmul		dword ptr [esi+(4+2*5)*4]				\
						__asm faddp		st(1),st								\
						__asm movss		xmm6, [edi+row*24+12]					\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm6, xmm3								\
						__asm addps		xmm7, xmm6								\
						__asm fld		dword ptr [edi+(row*6+3)*4]				\
						__asm fmul		dword ptr [esi+(4+3*5)*4]				\
						__asm faddp		st(1),st								\
						__asm movss		xmm6, [edi+row*24+16]					\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm6, xmm4								\
						__asm addps		xmm7, xmm6								\
						__asm fld		dword ptr [edi+(row*6+4)*4]				\
						__asm fmul		dword ptr [esi+(4+4*5)*4]				\
						__asm faddp		st(1),st								\
						__asm movss		xmm6, [edi+row*24+20]					\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm6, xmm5								\
						__asm addps		xmm7, xmm6								\
						__asm fld		dword ptr [edi+(row*6+5)*4]				\
						__asm fmul		dword ptr [esi+(4+5*5)*4]				\
						__asm faddp		st(1),st								\
						__asm fstp		dword ptr [eax+(row*5+4)*4]				\
						__asm movlps	[eax+row*20], xmm7						\
						__asm movhps	[eax+row*20+8], xmm7

						MUL_Nx6_6x5_INIT
						MUL_Nx6_6x5_ROW( 0 )
						MUL_Nx6_6x5_ROW( 1 )
						MUL_Nx6_6x5_ROW( 2 )
						MUL_Nx6_6x5_ROW( 3 )
						MUL_Nx6_6x5_ROW( 4 )

						return;
					}
					break;
				}
				case 6: {
					switch( l ) {
						case 1: {		// 6x6 * 6x1
							__asm {
								mov			esi, m2Ptr
								mov			edi, m1Ptr
								mov			eax, dstPtr
								movlps		xmm7, qword ptr [esi]
								movlps		xmm6, qword ptr [esi+8]
								shufps		xmm7, xmm7, 0x44
								shufps		xmm6, xmm6, 0x44
								movlps		xmm0, qword ptr [edi    ]
								movhps		xmm0, qword ptr [edi+ 24]
								mulps		xmm0, xmm7
								movlps		xmm3, qword ptr [edi+  8]
								movhps		xmm3, qword ptr [edi+ 32]
								mulps		xmm3, xmm6
								movlps		xmm1, qword ptr [edi+ 48]
								movhps		xmm1, qword ptr [edi+ 72]
								mulps		xmm1, xmm7
								movlps		xmm2, qword ptr [edi+ 96]
								movhps		xmm2, qword ptr [edi+120]
								mulps		xmm2, xmm7
								movlps		xmm4, qword ptr [edi+ 56]
								movhps		xmm4, qword ptr [edi+ 80]
								movlps		xmm5, qword ptr [edi+104]
								movhps		xmm5, qword ptr [edi+128]
								mulps		xmm4, xmm6
								movlps		xmm7, qword ptr [esi+16]
								addps		xmm0, xmm3
								shufps		xmm7, xmm7, 0x44
								mulps		xmm5, xmm6
								addps		xmm1, xmm4
								movlps		xmm3, qword ptr [edi+ 16]
								movhps		xmm3, qword ptr [edi+ 40]
								addps		xmm2, xmm5
								movlps		xmm4, qword ptr [edi+ 64]
								movhps		xmm4, qword ptr [edi+ 88]
								mulps		xmm3, xmm7
								movlps		xmm5, qword ptr [edi+112]
								movhps		xmm5, qword ptr [edi+136]
								addps		xmm0, xmm3
								mulps		xmm4, xmm7
								mulps		xmm5, xmm7
								addps		xmm1, xmm4
								addps		xmm2, xmm5
								movaps		xmm6, xmm0
								shufps		xmm0, xmm1, 0x88
								shufps		xmm6, xmm1, 0xDD
								movaps		xmm7, xmm2
								shufps		xmm7, xmm2, 0x88
								shufps		xmm2, xmm2, 0xDD
								addps		xmm0, xmm6
								addps		xmm2, xmm7
								movlps		[eax], xmm0
								movhps		[eax+8], xmm0
								movlps		[eax+16], xmm2
							}
							return;
						}
						case 2: {		// 6x6 * 6x2

							MUL_Nx6_6x2_INIT
							MUL_Nx6_6x2_ROW2( 0 )
							MUL_Nx6_6x2_ROW2( 1 )
							MUL_Nx6_6x2_ROW2( 2 )

							return;
						}
						case 3: {		// 6x6 * 6x3

							MUL_Nx6_6x3_INIT
							MUL_Nx6_6x3_ROW( 0 )
							MUL_Nx6_6x3_ROW( 1 )
							MUL_Nx6_6x3_ROW( 2 )
							MUL_Nx6_6x3_ROW( 3 )
							MUL_Nx6_6x3_ROW( 4 )
							MUL_Nx6_6x3_ROW( 5 )

							return;
						}
						case 4: {		// 6x6 * 6x4

							MUL_Nx6_6x4_INIT
							MUL_Nx6_6x4_ROW( 0 )
							MUL_Nx6_6x4_ROW( 1 )
							MUL_Nx6_6x4_ROW( 2 )
							MUL_Nx6_6x4_ROW( 3 )
							MUL_Nx6_6x4_ROW( 4 )
							MUL_Nx6_6x4_ROW( 5 )

							return;
						}
						case 5: {		// 6x6 * 6x5

							MUL_Nx6_6x5_INIT
							MUL_Nx6_6x5_ROW( 0 )
							MUL_Nx6_6x5_ROW( 1 )
							MUL_Nx6_6x5_ROW( 2 )
							MUL_Nx6_6x5_ROW( 3 )
							MUL_Nx6_6x5_ROW( 4 )
							MUL_Nx6_6x5_ROW( 5 )

							return;
						}
						case 6: {		// 6x6 * 6x6
							__asm {
								mov			ecx, dword ptr m2Ptr
								movlps		xmm3, qword ptr [ecx+72]
								mov			edx, dword ptr m1Ptr
								// Loading first 4 columns (upper 4 rows) of m2Ptr.
								movaps		xmm0, xmmword ptr [ecx]
								movlps		xmm1, qword ptr [ecx+24]
								movhps		xmm1, qword ptr [ecx+32]
								movaps		xmm2, xmmword ptr [ecx+48]
								movhps		xmm3, qword ptr [ecx+80]
								// Calculating first 4 elements in the first row of the destination matrix.
								movss		xmm4, dword ptr [edx]
								movss		xmm5, dword ptr [edx+4]
								mov			eax, dword ptr dstPtr
								shufps		xmm4, xmm4, 0
								movss		xmm6, dword ptr [edx+8]
								shufps		xmm5, xmm5, 0
								movss		xmm7, dword ptr [edx+12]
								mulps		xmm4, xmm0
								shufps		xmm6, xmm6, 0
								shufps		xmm7, xmm7, 0
								mulps		xmm5, xmm1
								mulps		xmm6, xmm2
								addps		xmm5, xmm4
								mulps		xmm7, xmm3
								addps		xmm6, xmm5
								addps		xmm7, xmm6
								movaps		xmmword ptr [eax], xmm7
								// Calculating first 4 elements in the second row of the destination matrix.
								movss		xmm4, dword ptr [edx+24]
								shufps		xmm4, xmm4, 0
								mulps		xmm4, xmm0
								movss		xmm5, dword ptr [edx+28]
								shufps		xmm5, xmm5, 0
								mulps		xmm5, xmm1
								movss		xmm6, dword ptr [edx+32]
								shufps		xmm6, xmm6, 0
								movss		xmm7, dword ptr [edx+36]
								shufps		xmm7, xmm7, 0
								mulps		xmm6, xmm2
								mulps		xmm7, xmm3
								addps		xmm7, xmm6
								addps		xmm5, xmm4
								addps		xmm7, xmm5
								// Calculating first 4 elements in the third row of the destination matrix.
								movss		xmm4, dword ptr [edx+48]
								movss		xmm5, dword ptr [edx+52]
								movlps		qword ptr [eax+24], xmm7 ; save 2nd
								movhps		qword ptr [eax+32], xmm7 ; row
								movss		xmm6, dword ptr [edx+56]
								movss		xmm7, dword ptr [edx+60]
								shufps		xmm4, xmm4, 0
								shufps		xmm5, xmm5, 0
								shufps		xmm6, xmm6, 0
								shufps		xmm7, xmm7, 0
								mulps		xmm4, xmm0
								mulps		xmm5, xmm1
								mulps		xmm6, xmm2
								mulps		xmm7, xmm3
								addps		xmm5, xmm4
								addps		xmm7, xmm6
								addps		xmm7, xmm5
								movaps		xmmword ptr [eax+48], xmm7
								// Calculating first 4 elements in the fourth row of the destination matrix.
								movss		xmm4, dword ptr [edx+72]
								movss		xmm5, dword ptr [edx+76]
								movss		xmm6, dword ptr [edx+80]
								movss		xmm7, dword ptr [edx+84]
								shufps		xmm4, xmm4, 0
								shufps		xmm5, xmm5, 0
								shufps		xmm6, xmm6, 0
								shufps		xmm7, xmm7, 0
								mulps		xmm4, xmm0
								mulps		xmm5, xmm1
								mulps		xmm6, xmm2
								mulps		xmm7, xmm3
								addps		xmm4, xmm5
								addps		xmm6, xmm4
								addps		xmm7, xmm6
								movlps		qword ptr [eax+72], xmm7
								movhps		qword ptr [eax+80], xmm7
								// Calculating first 4 elements in the fifth row of the destination matrix.
								movss		xmm4, dword ptr [edx+96]
								movss		xmm5, dword ptr [edx+100]
								movss		xmm6, dword ptr [edx+104]
								movss		xmm7, dword ptr [edx+108]
								shufps		xmm4, xmm4, 0
								shufps		xmm5, xmm5, 0
								shufps		xmm6, xmm6, 0
								shufps		xmm7, xmm7, 0
								mulps		xmm4, xmm0
								mulps		xmm5, xmm1
								mulps		xmm6, xmm2
								mulps		xmm7, xmm3
								addps		xmm5, xmm4
								addps		xmm7, xmm6
								addps		xmm7, xmm5
								movaps		xmmword ptr [eax+96], xmm7
								// Calculating first 4 elements in the sixth row of the destination matrix.
								movss		xmm4, dword ptr [edx+120]
								movss		xmm5, dword ptr [edx+124]
								movss		xmm6, dword ptr [edx+128]
								movss		xmm7, dword ptr [edx+132]
								shufps		xmm4, xmm4, 0
								shufps		xmm5, xmm5, 0
								shufps		xmm6, xmm6, 0
								shufps		xmm7, xmm7, 0
								mulps		xmm4, xmm0
								mulps		xmm5, xmm1
								mulps		xmm6, xmm2
								mulps		xmm7, xmm3
								addps		xmm4, xmm5
								addps		xmm6, xmm4
								addps		xmm7, xmm6
								movhps		qword ptr [eax+128], xmm7
								movlps		qword ptr [eax+120], xmm7
								// Loading first 4 columns (lower 2 rows) of m2Ptr.
								movlps		xmm0, qword ptr [ecx+96]
								movhps		xmm0, qword ptr [ecx+104]
								movlps		xmm1, qword ptr [ecx+120]
								movhps		xmm1, qword ptr [ecx+128]
								// Calculating first 4 elements in the first row of the destination matrix.
								movss		xmm2, dword ptr [edx+16]
								shufps		xmm2, xmm2, 0
								movss		xmm4, dword ptr [edx+40]
								movss		xmm3, dword ptr [edx+20]
								movss		xmm5, dword ptr [edx+44]
								movaps		xmm6, xmmword ptr [eax]
								movlps		xmm7, qword ptr [eax+24]
								shufps		xmm3, xmm3, 0
								shufps		xmm5, xmm5, 0
								movhps		xmm7, qword ptr [eax+32]
								shufps		xmm4, xmm4, 0
								mulps		xmm5, xmm1
								mulps		xmm2, xmm0
								mulps		xmm3, xmm1
								mulps		xmm4, xmm0
								addps		xmm6, xmm2
								addps		xmm7, xmm4
								addps		xmm7, xmm5
								addps		xmm6, xmm3
								movlps		qword ptr [eax+24], xmm7
								movaps		xmmword ptr [eax], xmm6
								movhps		qword ptr [eax+32], xmm7
								// Calculating first 4 elements in the third row of the destination matrix.
								movss		xmm2, dword ptr [edx+64]
								movss		xmm4, dword ptr [edx+88]
								movss		xmm5, dword ptr [edx+92]
								movss		xmm3, dword ptr [edx+68]
								movaps		xmm6, xmmword ptr [eax+48]
								movlps		xmm7, qword ptr [eax+72]
								movhps		xmm7, qword ptr [eax+80]
								shufps		xmm2, xmm2, 0
								shufps		xmm4, xmm4, 0
								shufps		xmm5, xmm5, 0
								shufps		xmm3, xmm3, 0
								mulps		xmm2, xmm0
								mulps		xmm4, xmm0
								mulps		xmm5, xmm1
								mulps		xmm3, xmm1
								addps		xmm6, xmm2
								addps		xmm6, xmm3
								addps		xmm7, xmm4
								addps		xmm7, xmm5
								movlps		qword ptr [eax+72], xmm7
								movaps		xmmword ptr [eax+48], xmm6
								movhps		qword ptr [eax+80], xmm7
								// Calculating first 4 elements in the fifth row of the destination matrix.
								movss		xmm2, dword ptr [edx+112]
								movss		xmm3, dword ptr [edx+116]
								movaps		xmm6, xmmword ptr [eax+96]
								shufps		xmm2, xmm2, 0
								shufps		xmm3, xmm3, 0
								mulps		xmm2, xmm0
								mulps		xmm3, xmm1
								addps		xmm6, xmm2
								addps		xmm6, xmm3
								movaps		xmmword ptr [eax+96], xmm6
								// Calculating first 4 elements in the sixth row of the destination matrix.
								movss		xmm4, dword ptr [edx+136]
								movss		xmm5, dword ptr [edx+140]
								movhps		xmm7, qword ptr [eax+128]
								movlps		xmm7, qword ptr [eax+120]
								shufps		xmm4, xmm4, 0
								shufps		xmm5, xmm5, 0
								mulps		xmm4, xmm0
								mulps		xmm5, xmm1
								addps		xmm7, xmm4
								addps		xmm7, xmm5
								// Calculating last 2 columns of the destination matrix.
								movlps		xmm0, qword ptr [ecx+16]
								movhps		xmm0, qword ptr [ecx+40]
								movhps		qword ptr [eax+128], xmm7
								movlps		qword ptr [eax+120], xmm7
								movlps		xmm2, qword ptr [ecx+64]
								movhps		xmm2, qword ptr [ecx+88]
								movaps		xmm3, xmm2
								shufps		xmm3, xmm3, 4Eh
								movlps		xmm4, qword ptr [ecx+112]
								movhps		xmm4, qword ptr [ecx+136]
								movaps		xmm5, xmm4
								shufps		xmm5, xmm5, 4Eh
								movlps		xmm6, qword ptr [edx]
								movhps		xmm6, qword ptr [edx+24]
								movaps		xmm7, xmm6
								shufps		xmm7, xmm7, 0F0h
								mulps		xmm7, xmm0
								shufps		xmm6, xmm6, 0A5h
								movaps		xmm1, xmm0
								shufps		xmm1, xmm1, 4Eh
								mulps		xmm1, xmm6
								addps		xmm7, xmm1
								movlps		xmm6, qword ptr [edx+8]
								movhps		xmm6, qword ptr [edx+32]
								movaps		xmm1, xmm6
								shufps		xmm1, xmm1, 0F0h
								shufps		xmm6, xmm6, 0A5h
								mulps		xmm1, xmm2
								mulps		xmm6, xmm3
								addps		xmm7, xmm1
								addps		xmm7, xmm6
								movhps		xmm6, qword ptr [edx+40]
								movlps		xmm6, qword ptr [edx+16]
								movaps		xmm1, xmm6
								shufps		xmm1, xmm1, 0F0h
								shufps		xmm6, xmm6, 0A5h
								mulps		xmm1, xmm4
								mulps		xmm6, xmm5
								addps		xmm7, xmm1
								addps		xmm7, xmm6
								movlps		qword ptr [eax+16], xmm7
								movhps		qword ptr [eax+40], xmm7
								movlps		xmm6, qword ptr [edx+48]
								movhps		xmm6, qword ptr [edx+72]
								movaps		xmm7, xmm6
								shufps		xmm7, xmm7, 0F0h
								mulps		xmm7, xmm0
								shufps		xmm6, xmm6, 0A5h
								movaps		xmm1, xmm0
								shufps		xmm1, xmm1, 4Eh
								mulps		xmm1, xmm6
								addps		xmm7, xmm1
								movhps		xmm6, qword ptr [edx+80]
								movlps		xmm6, qword ptr [edx+56]
								movaps		xmm1, xmm6
								shufps		xmm1, xmm1, 0F0h
								shufps		xmm6, xmm6, 0A5h
								mulps		xmm1, xmm2
								mulps		xmm6, xmm3
								addps		xmm7, xmm1
								addps		xmm7, xmm6
								movlps		xmm6, qword ptr [edx+64]
								movhps		xmm6, qword ptr [edx+88]
								movaps		xmm1, xmm6
								shufps		xmm1, xmm1, 0F0h
								shufps		xmm6, xmm6, 0A5h
								mulps		xmm1, xmm4
								mulps		xmm6, xmm5
								addps		xmm7, xmm1
								addps		xmm7, xmm6
								movlps		qword ptr [eax+64], xmm7
								movhps		qword ptr [eax+88], xmm7
								movlps		xmm6, qword ptr [edx+96]
								movhps		xmm6, qword ptr [edx+120]
								movaps		xmm7, xmm6
								shufps		xmm7, xmm7, 0F0h
								mulps		xmm7, xmm0
								shufps		xmm6, xmm6, 0A5h
								movaps		xmm1, xmm0
								shufps		xmm1, xmm1, 4Eh
								mulps		xmm1, xmm6
								addps		xmm7, xmm1
								movlps		xmm6, qword ptr [edx+104]
								movhps		xmm6, qword ptr [edx+128]
								movaps		xmm1, xmm6
								shufps		xmm1, xmm1, 0F0h
								shufps		xmm6, xmm6, 0A5h
								mulps		xmm1, xmm2
								mulps		xmm6, xmm3
								addps		xmm7, xmm1
								addps		xmm7, xmm6
								movlps		xmm6, qword ptr [edx+112]
								movhps		xmm6, qword ptr [edx+136]
								movaps		xmm1, xmm6
								shufps		xmm1, xmm1, 0F0h
								shufps		xmm6, xmm6, 0A5h
								mulps		xmm1, xmm4
								mulps		xmm6, xmm5
								addps		xmm7, xmm1
								addps		xmm7, xmm6
								movlps		qword ptr [eax+112], xmm7
								movhps		qword ptr [eax+136], xmm7
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
CSIMD_SSE::matX_TransposeMultiplyMatX

	optimizes the following transpose matrix multiplications:

	Nx6 * NxN
	6xN * 6x6

	with N in the range [1-6].
============
*/
void VPCALL CSIMD_SSE::matX_TransposeMultiplyMatX( CMatXD &dst, const CMatXD &m1, const CMatXD &m2 ) {
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
			if ( !((k^6)|(l^1)) ) {			// 1x6 * 1x1
				__asm {
					mov		esi, m2Ptr
					mov		edi, m1Ptr
					mov		eax, dstPtr
					movss	xmm0, [esi]
					shufps	xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
					movaps	xmm1, xmm0
					mulps	xmm0, [edi]
					mulps	xmm1, [edi+16]
					movaps	[eax], xmm0
					movlps	[eax+16], xmm1
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
			if ( !((k^6)|(l^2)) ) {			// 2x6 * 2x2
				#define MUL_2xN_2x2_INIT								\
				__asm mov		esi, m2Ptr								\
				__asm mov		edi, m1Ptr								\
				__asm mov		eax, dstPtr								\
				__asm movlps	xmm0, [esi]								\
				__asm shufps	xmm0, xmm0, R_SHUFFLEPS( 0, 1, 0, 1 )	\
				__asm movlps	xmm1, [esi+8]							\
				__asm shufps	xmm1, xmm1, R_SHUFFLEPS( 0, 1, 0, 1 )

				#define MUL_2xN_2x2_ROW2( N, row )						\
				__asm movlps	xmm6, [edi+(row+0*N)*4]					\
				__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 1, 1 )	\
				__asm movlps	xmm7, [edi+(row+1*N)*4]					\
				__asm shufps	xmm7, xmm7, R_SHUFFLEPS( 0, 0, 1, 1 )	\
				__asm mulps		xmm6, xmm0								\
				__asm mulps		xmm7, xmm1								\
				__asm addps		xmm6, xmm7								\
				__asm movaps	[eax+(row*2)*4], xmm6

				MUL_2xN_2x2_INIT
				MUL_2xN_2x2_ROW2( 6, 0 )
				MUL_2xN_2x2_ROW2( 6, 2 )
				MUL_2xN_2x2_ROW2( 6, 4 )

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
			if ( !((k^6)|(l^3)) ) {			// 3x6 * 3x3

				#define MUL_3xN_3x3_INIT								\
				__asm mov		esi, m2Ptr								\
				__asm mov		edi, m1Ptr								\
				__asm mov		eax, dstPtr								\
				__asm movss		xmm0, [esi+(0*3+0)*4]					\
				__asm movhps	xmm0, [esi+(0*3+1)*4]					\
				__asm movss		xmm1, [esi+(1*3+0)*4]					\
				__asm movhps	xmm1, [esi+(1*3+1)*4]					\
				__asm movss		xmm2, [esi+(2*3+0)*4]					\
				__asm movhps	xmm2, [esi+(2*3+1)*4]

				#define MUL_3xN_3x3_INIT_ROW4							\
				__asm shufps	xmm0, xmm0, R_SHUFFLEPS( 0, 2, 3, 0 )	\
				__asm shufps	xmm1, xmm1, R_SHUFFLEPS( 0, 2, 3, 0 )	\
				__asm shufps	xmm2, xmm2, R_SHUFFLEPS( 0, 2, 3, 0 )

				#define MUL_3xN_3x3_ROW4( N, row )						\
				__asm movlps	xmm3, [edi+(row+0*N+0)*4]				\
				__asm shufps	xmm3, xmm3, R_SHUFFLEPS( 0, 0, 0, 1 )	\
				__asm movlps	xmm4, [edi+(row+1*N+0)*4]				\
				__asm shufps	xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 1 )	\
				__asm movlps	xmm5, [edi+(row+2*N+0)*4]				\
				__asm shufps	xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 1 )	\
				__asm mulps		xmm3, xmm0								\
				__asm mulps		xmm4, xmm1								\
				__asm mulps		xmm5, xmm2								\
				__asm addps		xmm3, xmm4								\
				__asm addps		xmm3, xmm5								\
				__asm movaps	[eax+(row*3+0)*4], xmm3					\
				__asm shufps	xmm0, xmm0, R_SHUFFLEPS( 1, 2, 3, 1 )	\
				__asm shufps	xmm1, xmm1, R_SHUFFLEPS( 1, 2, 3, 1 )	\
				__asm shufps	xmm2, xmm2, R_SHUFFLEPS( 1, 2, 3, 1 )	\
				__asm movlps	xmm3, [edi+(row+0*N+1)*4]				\
				__asm shufps	xmm3, xmm3, R_SHUFFLEPS( 0, 0, 1, 1 )	\
				__asm movlps	xmm4, [edi+(row+1*N+1)*4]				\
				__asm shufps	xmm4, xmm4, R_SHUFFLEPS( 0, 0, 1, 1 )	\
				__asm movlps	xmm5, [edi+(row+2*N+1)*4]				\
				__asm shufps	xmm5, xmm5, R_SHUFFLEPS( 0, 0, 1, 1 )	\
				__asm mulps		xmm3, xmm0								\
				__asm mulps		xmm4, xmm1								\
				__asm mulps		xmm5, xmm2								\
				__asm addps		xmm3, xmm4								\
				__asm addps		xmm3, xmm5								\
				__asm movaps	[eax+(row*3+4)*4], xmm3					\
				__asm shufps	xmm0, xmm0, R_SHUFFLEPS( 1, 2, 3, 1 )	\
				__asm shufps	xmm1, xmm1, R_SHUFFLEPS( 1, 2, 3, 1 )	\
				__asm shufps	xmm2, xmm2, R_SHUFFLEPS( 1, 2, 3, 1 )	\
				__asm movlps	xmm3, [edi+(row+0*N+2)*4]				\
				__asm shufps	xmm3, xmm3, R_SHUFFLEPS( 0, 1, 1, 1 )	\
				__asm movlps	xmm4, [edi+(row+1*N+2)*4]				\
				__asm shufps	xmm4, xmm4, R_SHUFFLEPS( 0, 1, 1, 1 )	\
				__asm movlps	xmm5, [edi+(row+2*N+2)*4]				\
				__asm shufps	xmm5, xmm5, R_SHUFFLEPS( 0, 1, 1, 1 )	\
				__asm mulps		xmm3, xmm0								\
				__asm mulps		xmm4, xmm1								\
				__asm mulps		xmm5, xmm2								\
				__asm addps		xmm3, xmm4								\
				__asm addps		xmm3, xmm5								\
				__asm movaps	[eax+(row*3+8)*4], xmm3

				#define MUL_3xN_3x3_INIT_ROW4_ROW4						\
				__asm shufps	xmm0, xmm0, R_SHUFFLEPS( 1, 2, 3, 0 )	\
				__asm shufps	xmm1, xmm1, R_SHUFFLEPS( 1, 2, 3, 0 )	\
				__asm shufps	xmm2, xmm2, R_SHUFFLEPS( 1, 2, 3, 0 )

				#define MUL_3xN_3x3_INIT_ROW4_ROW						\
				__asm shufps	xmm0, xmm0, R_SHUFFLEPS( 1, 1, 2, 3 )	\
				__asm shufps	xmm1, xmm1, R_SHUFFLEPS( 1, 1, 2, 3 )	\
				__asm shufps	xmm2, xmm2, R_SHUFFLEPS( 1, 1, 2, 3 )

				#define MUL_3xN_3x3_ROW( N, row )						\
				__asm movss		xmm3, [edi+(row+0*N)*4]					\
				__asm shufps	xmm3, xmm3, R_SHUFFLEPS( 0, 0, 0, 0 )	\
				__asm movss		xmm4, [edi+(row+1*N)*4]					\
				__asm shufps	xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )	\
				__asm movss		xmm5, [edi+(row+2*N)*4]					\
				__asm shufps	xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )	\
				__asm mulps		xmm3, xmm0								\
				__asm mulps		xmm4, xmm1								\
				__asm mulps		xmm5, xmm2								\
				__asm addps		xmm3, xmm4								\
				__asm addps		xmm3, xmm5								\
				__asm movss		[eax+(row*3+0)*4], xmm3					\
				__asm movhps	[eax+(row*3+1)*4], xmm3

				MUL_3xN_3x3_INIT
				MUL_3xN_3x3_INIT_ROW4
				MUL_3xN_3x3_ROW4( 6, 0 )
				MUL_3xN_3x3_INIT_ROW4_ROW
				MUL_3xN_3x3_ROW( 6, 4 )
				MUL_3xN_3x3_ROW( 6, 5 )

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
			if ( !((k^6)|(l^4)) ) {			// 4x6 * 4x4

				#define MUL_4xN_4x4_INIT								\
				__asm mov		esi, m2Ptr								\
				__asm mov		edi, m1Ptr								\
				__asm mov		eax, dstPtr								\
				__asm movaps	xmm0, [esi]								\
				__asm movaps	xmm1, [esi+16]							\
				__asm movaps	xmm2, [esi+32]							\
				__asm movaps	xmm3, [esi+48]

				#define MUL_4xN_4x4_ROW( N, row )						\
				__asm movss		xmm7, [edi+(row+0*N)*4]					\
				__asm shufps	xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )	\
				__asm mulps		xmm7, xmm0								\
				__asm movss		xmm6, [edi+(row+1*N)*4]					\
				__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
				__asm mulps		xmm6, xmm1								\
				__asm addps		xmm7, xmm6								\
				__asm movss		xmm6, [edi+(row+2*N)*4]					\
				__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
				__asm mulps		xmm6, xmm2								\
				__asm addps		xmm7, xmm6								\
				__asm movss		xmm6, [edi+(row+3*N)*4]					\
				__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
				__asm mulps		xmm6, xmm3								\
				__asm addps		xmm7, xmm6								\
				__asm movaps	[eax+row*16], xmm7

				MUL_4xN_4x4_INIT
				MUL_4xN_4x4_ROW( 6, 0 )
				MUL_4xN_4x4_ROW( 6, 1 )
				MUL_4xN_4x4_ROW( 6, 2 )
				MUL_4xN_4x4_ROW( 6, 3 )
				MUL_4xN_4x4_ROW( 6, 4 )
				MUL_4xN_4x4_ROW( 6, 5 )

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
			if ( !((k^6)|(l^5)) ) {			// 5x6 * 5x5

				#define MUL_5xN_5x5_INIT								\
				__asm mov		esi, m2Ptr								\
				__asm mov		edi, m1Ptr								\
				__asm mov		eax, dstPtr								\
				__asm movlps	xmm0, [esi+ 0*4]						\
				__asm movhps	xmm0, [esi+ 2*4]						\
				__asm movlps	xmm1, [esi+ 5*4]						\
				__asm movhps	xmm1, [esi+ 7*4]						\
				__asm movlps	xmm2, [esi+10*4]						\
				__asm movhps	xmm2, [esi+12*4]						\
				__asm movlps	xmm3, [esi+15*4]						\
				__asm movhps	xmm3, [esi+17*4]						\
				__asm movlps	xmm4, [esi+20*4]						\
				__asm movhps	xmm4, [esi+22*4]

				#define MUL_5xN_5x5_ROW( N, row )						\
				__asm movss		xmm6, [edi+(row+0*N)*4]					\
				__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
				__asm mulps		xmm6, xmm0								\
				__asm fld		dword ptr [edi+(row+0*N)*4]				\
				__asm fmul		dword ptr [esi+ 4*4]					\
				__asm movss		xmm5, [edi+(row+1*N)*4]					\
				__asm shufps	xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )	\
				__asm mulps		xmm5, xmm1								\
				__asm addps		xmm6, xmm5								\
				__asm fld		dword ptr [edi+(row+1*N)*4]				\
				__asm fmul		dword ptr [esi+ 9*4]					\
				__asm faddp		st(1),st								\
				__asm movss		xmm5, [edi+(row+2*N)*4]					\
				__asm shufps	xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )	\
				__asm mulps		xmm5, xmm2								\
				__asm addps		xmm6, xmm5								\
				__asm fld		dword ptr [edi+(row+2*N)*4]				\
				__asm fmul		dword ptr [esi+14*4]					\
				__asm faddp		st(1),st								\
				__asm movss		xmm5, [edi+(row+3*N)*4]					\
				__asm shufps	xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )	\
				__asm mulps		xmm5, xmm3								\
				__asm addps		xmm6, xmm5								\
				__asm fld		dword ptr [edi+(row+3*N)*4]				\
				__asm fmul		dword ptr [esi+19*4]					\
				__asm faddp		st(1),st								\
				__asm movss		xmm5, [edi+(row+4*N)*4]					\
				__asm shufps	xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )	\
				__asm mulps		xmm5, xmm4								\
				__asm addps		xmm6, xmm5								\
				__asm fld		dword ptr [edi+(row+4*N)*4]				\
				__asm fmul		dword ptr [esi+24*4]					\
				__asm faddp		st(1),st								\
				__asm fstp		dword ptr [eax+(row*5+4)*4]				\
				__asm movlps	[eax+(row*5+0)*4], xmm6					\
				__asm movhps	[eax+(row*5+2)*4], xmm6

				MUL_5xN_5x5_INIT
				MUL_5xN_5x5_ROW( 6, 0 )
				MUL_5xN_5x5_ROW( 6, 1 )
				MUL_5xN_5x5_ROW( 6, 2 )
				MUL_5xN_5x5_ROW( 6, 3 )
				MUL_5xN_5x5_ROW( 6, 4 )
				MUL_5xN_5x5_ROW( 6, 5 )

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
			if ( !(l^6) ) {
				switch( k ) {
					case 1: {					// 6x1 * 6x6
						#define MUL_6xN_6x6_FIRST4COLUMNS_INIT					\
						__asm mov		esi, m2Ptr								\
						__asm mov		edi, m1Ptr								\
						__asm mov		eax, dstPtr								\
						__asm movlps	xmm0, [esi+ 0*4]						\
						__asm movhps	xmm0, [esi+ 2*4]						\
						__asm movlps	xmm1, [esi+ 6*4]						\
						__asm movhps	xmm1, [esi+ 8*4]						\
						__asm movlps	xmm2, [esi+12*4]						\
						__asm movhps	xmm2, [esi+14*4]						\
						__asm movlps	xmm3, [esi+18*4]						\
						__asm movhps	xmm3, [esi+20*4]						\
						__asm movlps	xmm4, [esi+24*4]						\
						__asm movhps	xmm4, [esi+26*4]						\
						__asm movlps	xmm5, [esi+30*4]						\
						__asm movhps	xmm5, [esi+32*4]

						#define MUL_6xN_6x6_FIRST4COLUMNS_ROW( N, row )			\
						__asm movss		xmm7, [edi+(row+0*N)*4]					\
						__asm shufps	xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm7, xmm0								\
						__asm movss		xmm6, [edi+(row+1*N)*4]					\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm6, xmm1								\
						__asm addps		xmm7, xmm6								\
						__asm movss		xmm6, [edi+(row+2*N)*4]					\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm6, xmm2								\
						__asm addps		xmm7, xmm6								\
						__asm movss		xmm6, [edi+(row+3*N)*4]					\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm6, xmm3								\
						__asm addps		xmm7, xmm6								\
						__asm movss		xmm6, [edi+(row+4*N)*4]					\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm6, xmm4								\
						__asm addps		xmm7, xmm6								\
						__asm movss		xmm6, [edi+(row+5*N)*4]					\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm6, xmm5								\
						__asm addps		xmm7, xmm6								\
						__asm movlps	[eax+(row*6+0)*4], xmm7					\
						__asm movhps	[eax+(row*6+2)*4], xmm7

						#define MUL_6xN_6x6_LAST2COLUMNS_INIT					\
						__asm movlps	xmm0, [esi+ 4*4]						\
						__asm movlps	xmm1, [esi+10*4]						\
						__asm shufps	xmm0, xmm0, R_SHUFFLEPS( 0, 1, 0, 1 )	\
						__asm shufps	xmm1, xmm1, R_SHUFFLEPS( 0, 1, 0, 1 )	\
						__asm movlps	xmm2, [esi+16*4]						\
						__asm movlps	xmm3, [esi+22*4]						\
						__asm shufps	xmm2, xmm2, R_SHUFFLEPS( 0, 1, 0, 1 )	\
						__asm shufps	xmm3, xmm3, R_SHUFFLEPS( 0, 1, 0, 1 )	\
						__asm movlps	xmm4, [esi+28*4]						\
						__asm movlps	xmm5, [esi+34*4]						\
						__asm shufps	xmm4, xmm4, R_SHUFFLEPS( 0, 1, 0, 1 )	\
						__asm shufps	xmm5, xmm5, R_SHUFFLEPS( 0, 1, 0, 1 )

						#define MUL_6xN_6x6_LAST2COLUMNS_ROW2( N, row )			\
						__asm movlps	xmm7, [edi+(row*2+0*N)*4]				\
						__asm shufps	xmm7, xmm7, R_SHUFFLEPS( 0, 0, 1, 1 )	\
						__asm mulps		xmm7, xmm0								\
						__asm movlps	xmm6, [edi+(row*2+1*N)*4]				\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 1, 1 )	\
						__asm mulps		xmm6, xmm1								\
						__asm addps		xmm7, xmm6								\
						__asm movlps	xmm6, [edi+(row*2+2*N)*4]				\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 1, 1 )	\
						__asm mulps		xmm6, xmm2								\
						__asm addps		xmm7, xmm6								\
						__asm movlps	xmm6, [edi+(row*2+3*N)*4]				\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 1, 1 )	\
						__asm mulps		xmm6, xmm3								\
						__asm addps		xmm7, xmm6								\
						__asm movlps	xmm6, [edi+(row*2+4*N)*4]				\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 1, 1 )	\
						__asm mulps		xmm6, xmm4								\
						__asm addps		xmm7, xmm6								\
						__asm movlps	xmm6, [edi+(row*2+5*N)*4]				\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 1, 1 )	\
						__asm mulps		xmm6, xmm5								\
						__asm addps		xmm7, xmm6								\
						__asm movlps	[eax+(row*12+ 4)*4], xmm7				\
						__asm movhps	[eax+(row*12+10)*4], xmm7

						#define MUL_6xN_6x6_LAST2COLUMNS_ROW( N, row )			\
						__asm movss		xmm7, [edi+(1*N-1)*4]					\
						__asm shufps	xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm7, xmm0								\
						__asm movss		xmm6, [edi+(2*N-1)*4]					\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm6, xmm1								\
						__asm addps		xmm7, xmm6								\
						__asm movss		xmm6, [edi+(3*N-1)*4]					\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm6, xmm2								\
						__asm addps		xmm7, xmm6								\
						__asm movss		xmm6, [edi+(4*N-1)*4]					\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm6, xmm3								\
						__asm addps		xmm7, xmm6								\
						__asm movss		xmm6, [edi+(5*N-1)*4]					\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm6, xmm4								\
						__asm addps		xmm7, xmm6								\
						__asm movss		xmm6, [edi+(6*N-1)*4]					\
						__asm shufps	xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )	\
						__asm mulps		xmm6, xmm5								\
						__asm addps		xmm7, xmm6								\
						__asm movlps	[eax+(row*6+4)*4], xmm7

						MUL_6xN_6x6_FIRST4COLUMNS_INIT
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 1, 0 )
						MUL_6xN_6x6_LAST2COLUMNS_INIT
						MUL_6xN_6x6_LAST2COLUMNS_ROW( 1, 0 )

						return;
					}
					case 2: {					// 6x2 * 6x6

						MUL_6xN_6x6_FIRST4COLUMNS_INIT
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 2, 0 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 2, 1 )
						MUL_6xN_6x6_LAST2COLUMNS_INIT
						MUL_6xN_6x6_LAST2COLUMNS_ROW2( 2, 0 )

						return;
					}
					case 3: {					// 6x3 * 6x6

						MUL_6xN_6x6_FIRST4COLUMNS_INIT
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 3, 0 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 3, 1 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 3, 2 )
						MUL_6xN_6x6_LAST2COLUMNS_INIT
						MUL_6xN_6x6_LAST2COLUMNS_ROW2( 3, 0 )
						MUL_6xN_6x6_LAST2COLUMNS_ROW( 3, 2 )

						return;
					}
					case 4: {					// 6x4 * 6x6

						MUL_6xN_6x6_FIRST4COLUMNS_INIT
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 4, 0 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 4, 1 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 4, 2 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 4, 3 )
						MUL_6xN_6x6_LAST2COLUMNS_INIT
						MUL_6xN_6x6_LAST2COLUMNS_ROW2( 4, 0 )
						MUL_6xN_6x6_LAST2COLUMNS_ROW2( 4, 1 )

						return;
					}
					case 5: {					// 6x5 * 6x6

						MUL_6xN_6x6_FIRST4COLUMNS_INIT
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 5, 0 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 5, 1 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 5, 2 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 5, 3 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 5, 4 )
						MUL_6xN_6x6_LAST2COLUMNS_INIT
						MUL_6xN_6x6_LAST2COLUMNS_ROW2( 5, 0 )
						MUL_6xN_6x6_LAST2COLUMNS_ROW2( 5, 1 )
						MUL_6xN_6x6_LAST2COLUMNS_ROW( 5, 4 )

						return;
					}
					case 6: {					// 6x6 * 6x6

						MUL_6xN_6x6_FIRST4COLUMNS_INIT
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 6, 0 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 6, 1 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 6, 2 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 6, 3 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 6, 4 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 6, 5 )
						MUL_6xN_6x6_LAST2COLUMNS_INIT
						MUL_6xN_6x6_LAST2COLUMNS_ROW2( 6, 0 )
						MUL_6xN_6x6_LAST2COLUMNS_ROW2( 6, 1 )
						MUL_6xN_6x6_LAST2COLUMNS_ROW2( 6, 2 )

						return;
					}
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
CSIMD_SSE::matX_LowerTriangularsolve

  solves x in Lx = b for the n * n sub-matrix of L
  if skip > 0 the first skip elements of x are assumed to be valid already
  L has to be a lower triangular matrix with (implicit) ones on the diagonal
  x == b is allowed
============
*/
void VPCALL CSIMD_SSE::matX_LowerTriangularSolve( const CMatXD &L, float *x, const float *b, const int n, int skip ) {
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

	// this code assumes n > 4
	__asm {
		push		ebx
		mov			eax, skip				// eax = i
		shl			eax, 2					// eax = i*4
		mov			edx, n					// edx = n
		shl			edx, 2					// edx = n*4
		mov			esi, x					// esi = x
		mov			edi, lptr				// edi = lptr
		add			esi, eax
		add			edi, eax
		mov			ebx, b					// ebx = b

		// check for aligned memory
		mov			ecx, nc
		shl			ecx, 2
		or			ecx, esi
		or			ecx, edi
		and			ecx, 15
		jnz			loopurow

		// aligned
	looprow:
		mov			ecx, eax
		neg			ecx
		movaps		xmm0, [esi+ecx]
		mulps		xmm0, [edi+ecx]
		add			ecx, 12*4
		jg			donedot8
	dot8:
		movaps		xmm1, [esi+ecx-(8*4)]
		mulps		xmm1, [edi+ecx-(8*4)]
		addps		xmm0, xmm1
		movaps		xmm3, [esi+ecx-(4*4)]
		mulps		xmm3, [edi+ecx-(4*4)]
		addps		xmm0, xmm3
		add			ecx, 8*4
		jle			dot8
	donedot8:
		sub			ecx, 4*4
		jg			donedot4
	//dot4:
		movaps		xmm1, [esi+ecx-(4*4)]
		mulps		xmm1, [edi+ecx-(4*4)]
		addps		xmm0, xmm1
		add			ecx, 4*4
	donedot4:
		movhlps		xmm1, xmm0
		addps		xmm0, xmm1
		movaps		xmm1, xmm0
		shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 0, 0, 0 )
		addss		xmm0, xmm1
		sub			ecx, 4*4
		jz			dot0
		add			ecx, 4
		jz			dot1
		add			ecx, 4
		jz			dot2
	//dot3:
		movss		xmm1, [esi-(3*4)]
		mulss		xmm1, [edi-(3*4)]
		addss		xmm0, xmm1
	dot2:
		movss		xmm3, [esi-(2*4)]
		mulss		xmm3, [edi-(2*4)]
		addss		xmm0, xmm3
	dot1:
		movss		xmm5, [esi-(1*4)]
		mulss		xmm5, [edi-(1*4)]
		addss		xmm0, xmm5
	dot0:
		movss		xmm1, [ebx+eax]
		subss		xmm1, xmm0
		movss		[esi], xmm1
		add			eax, 4
		cmp			eax, edx
		jge			done
		add			esi, 4
		mov			ecx, nc
		shl			ecx, 2
		add			edi, ecx
		add			edi, 4
		jmp			looprow

		// unaligned
	loopurow:
		mov			ecx, eax
		neg			ecx
		movups		xmm0, [esi+ecx]
		movups		xmm1, [edi+ecx]
		mulps		xmm0, xmm1
		add			ecx, 12*4
		jg			doneudot8
	udot8:
		movups		xmm1, [esi+ecx-(8*4)]
		movups		xmm2, [edi+ecx-(8*4)]
		mulps		xmm1, xmm2
		addps		xmm0, xmm1
		movups		xmm3, [esi+ecx-(4*4)]
		movups		xmm4, [edi+ecx-(4*4)]
		mulps		xmm3, xmm4
		addps		xmm0, xmm3
		add			ecx, 8*4
		jle			udot8
	doneudot8:
		sub			ecx, 4*4
		jg			doneudot4
	//udot4:
		movups		xmm1, [esi+ecx-(4*4)]
		movups		xmm2, [edi+ecx-(4*4)]
		mulps		xmm1, xmm2
		addps		xmm0, xmm1
		add			ecx, 4*4
	doneudot4:
		movhlps		xmm1, xmm0
		addps		xmm0, xmm1
		movaps		xmm1, xmm0
		shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 0, 0, 0 )
		addss		xmm0, xmm1
		sub			ecx, 4*4
		jz			udot0
		add			ecx, 4
		jz			udot1
		add			ecx, 4
		jz			udot2
	//udot3:
		movss		xmm1, [esi-(3*4)]
		movss		xmm2, [edi-(3*4)]
		mulss		xmm1, xmm2
		addss		xmm0, xmm1
	udot2:
		movss		xmm3, [esi-(2*4)]
		movss		xmm4, [edi-(2*4)]
		mulss		xmm3, xmm4
		addss		xmm0, xmm3
	udot1:
		movss		xmm5, [esi-(1*4)]
		movss		xmm6, [edi-(1*4)]
		mulss		xmm5, xmm6
		addss		xmm0, xmm5
	udot0:
		movss		xmm1, [ebx+eax]
		subss		xmm1, xmm0
		movss		[esi], xmm1
		add			eax, 4
		cmp			eax, edx
		jge			done
		add			esi, 4
		mov			ecx, nc
		shl			ecx, 2
		add			edi, ecx
		add			edi, 4
		jmp			loopurow
	done:
		pop			ebx
	}
}

/*
============
CSIMD_SSE::matX_LowerTriangularsolveTranspose

  solves x in L'x = b for the n * n sub-matrix of L
  L has to be a lower triangular matrix with (implicit) ones on the diagonal
  x == b is allowed
============
*/
void VPCALL CSIMD_SSE::matX_LowerTriangularSolveTranspose( const CMatXD &L, float *x, const float *b, const int n ) {
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

#if 1

	int i, j, m;
	float *xptr;
	double s0;

	// if the number of columns is not a multiple of 2 we're screwed for alignment.
	// however, if the number of columns is a multiple of 2 but the number of to be
	// processed rows is not a multiple of 2 we can still run 8 sf_u8 aligned
	m = n;
	if ( m & 1 ) {

		m--;
		x[m] = b[m];

		lptr = L.toFloatPtr() + m * nc + m - 4;
		xptr = x + m;
		__asm {
			push		ebx
			mov			eax, m					// eax = i
			mov			esi, xptr				// esi = xptr
			mov			edi, lptr				// edi = lptr
			mov			ebx, b					// ebx = b
			mov			edx, nc					// edx = nc*sizeof(float)
			shl			edx, 2
		process4rows_1:
			movlps		xmm0, [ebx+eax*4-16]	// load b[i-2], b[i-1]
			movhps		xmm0, [ebx+eax*4-8]		// load b[i-4], b[i-3]
			xor			ecx, ecx
			sub			eax, m
			neg			eax
			jz			done4x4_1
		process4x4_1:	// process 4x4 blocks
			movlps		xmm2, [edi+0]
			movhps		xmm2, [edi+8]
			add			edi, edx
			movss		xmm1, [esi+4*ecx+0]
			shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
			movlps		xmm3, [edi+0]
			movhps		xmm3, [edi+8]
			add			edi, edx
			mulps		xmm1, xmm2
			subps		xmm0, xmm1
			movss		xmm1, [esi+4*ecx+4]
			shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
			movlps		xmm4, [edi+0]
			movhps		xmm4, [edi+8]
			add			edi, edx
			mulps		xmm1, xmm3
			subps		xmm0, xmm1
			movss		xmm1, [esi+4*ecx+8]
			shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
			movlps		xmm5, [edi+0]
			movhps		xmm5, [edi+8]
			add			edi, edx
			mulps		xmm1, xmm4
			subps		xmm0, xmm1
			movss		xmm1, [esi+4*ecx+12]
			shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
			add			ecx, 4
			cmp			ecx, eax
			mulps		xmm1, xmm5
			subps		xmm0, xmm1
			jl			process4x4_1
		done4x4_1:		// process left over of the 4 rows
			movlps		xmm2, [edi+0]
			movhps		xmm2, [edi+8]
			movss		xmm1, [esi+4*ecx]
			shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
			mulps		xmm1, xmm2
			subps		xmm0, xmm1
			imul		ecx, edx
			sub			edi, ecx
			neg			eax

			add			eax, m
			sub			eax, 4
			movaps		xmm1, xmm0
			shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 1, 1, 1 )
			movaps		xmm2, xmm0
			shufps		xmm2, xmm2, R_SHUFFLEPS( 2, 2, 2, 2 )
			movaps		xmm3, xmm0
			shufps		xmm3, xmm3, R_SHUFFLEPS( 3, 3, 3, 3 )
			sub			edi, edx
			movss		[esi-4], xmm3			// xptr[-1] = s3
			movss		xmm4, xmm3
			movss		xmm5, xmm3
			mulss		xmm3, [edi+8]			// lptr[-1*nc+2] * s3
			mulss		xmm4, [edi+4]			// lptr[-1*nc+1] * s3
			mulss		xmm5, [edi+0]			// lptr[-1*nc+0] * s3
			subss		xmm2, xmm3
			movss		[esi-8], xmm2			// xptr[-2] = s2
			movss		xmm6, xmm2
			sub			edi, edx
			subss		xmm0, xmm5
			subss		xmm1, xmm4
			mulss		xmm2, [edi+4]			// lptr[-2*nc+1] * s2
			mulss		xmm6, [edi+0]			// lptr[-2*nc+0] * s2
			subss		xmm1, xmm2
			movss		[esi-12], xmm1			// xptr[-3] = s1
			subss		xmm0, xmm6
			sub			edi, edx
			cmp			eax, 4
			mulss		xmm1, [edi+0]			// lptr[-3*nc+0] * s1
			subss		xmm0, xmm1
			movss		[esi-16], xmm0			// xptr[-4] = s0
			jl			done4rows_1
			sub			edi, edx
			sub			edi, 16
			sub			esi, 16
			jmp			process4rows_1
		done4rows_1:
			pop			ebx
		}

	} else {

		lptr = L.toFloatPtr() + m * nc + m - 4;
		xptr = x + m;
		__asm {
			push		ebx
			mov			eax, m					// eax = i
			mov			esi, xptr				// esi = xptr
			mov			edi, lptr				// edi = lptr
			mov			ebx, b					// ebx = b
			mov			edx, nc					// edx = nc*sizeof(float)
			shl			edx, 2
		process4rows:
			movlps		xmm0, [ebx+eax*4-16]	// load b[i-2], b[i-1]
			movhps		xmm0, [ebx+eax*4-8]		// load b[i-4], b[i-3]
			sub			eax, m
			jz			done4x4
			neg			eax
			xor			ecx, ecx
		process4x4:		// process 4x4 blocks
			movlps		xmm2, [edi+0]
			movhps		xmm2, [edi+8]
			add			edi, edx
			movss		xmm1, [esi+4*ecx+0]
			shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
			movlps		xmm3, [edi+0]
			movhps		xmm3, [edi+8]
			add			edi, edx
			mulps		xmm1, xmm2
			subps		xmm0, xmm1
			movss		xmm1, [esi+4*ecx+4]
			shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
			movlps		xmm4, [edi+0]
			movhps		xmm4, [edi+8]
			add			edi, edx
			mulps		xmm1, xmm3
			subps		xmm0, xmm1
			movss		xmm1, [esi+4*ecx+8]
			shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
			movlps		xmm5, [edi+0]
			movhps		xmm5, [edi+8]
			add			edi, edx
			mulps		xmm1, xmm4
			subps		xmm0, xmm1
			movss		xmm1, [esi+4*ecx+12]
			shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
			add			ecx, 4
			cmp			ecx, eax
			mulps		xmm1, xmm5
			subps		xmm0, xmm1
			jl			process4x4
			imul		ecx, edx
			sub			edi, ecx
			neg			eax
		done4x4:		// process left over of the 4 rows
			add			eax, m
			sub			eax, 4
			movaps		xmm1, xmm0
			shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 1, 1, 1 )
			movaps		xmm2, xmm0
			shufps		xmm2, xmm2, R_SHUFFLEPS( 2, 2, 2, 2 )
			movaps		xmm3, xmm0
			shufps		xmm3, xmm3, R_SHUFFLEPS( 3, 3, 3, 3 )
			sub			edi, edx
			movss		[esi-4], xmm3			// xptr[-1] = s3
			movss		xmm4, xmm3
			movss		xmm5, xmm3
			mulss		xmm3, [edi+8]			// lptr[-1*nc+2] * s3
			mulss		xmm4, [edi+4]			// lptr[-1*nc+1] * s3
			mulss		xmm5, [edi+0]			// lptr[-1*nc+0] * s3
			subss		xmm2, xmm3
			movss		[esi-8], xmm2			// xptr[-2] = s2
			movss		xmm6, xmm2
			sub			edi, edx
			subss		xmm0, xmm5
			subss		xmm1, xmm4
			mulss		xmm2, [edi+4]			// lptr[-2*nc+1] * s2
			mulss		xmm6, [edi+0]			// lptr[-2*nc+0] * s2
			subss		xmm1, xmm2
			movss		[esi-12], xmm1			// xptr[-3] = s1
			subss		xmm0, xmm6
			sub			edi, edx
			cmp			eax, 4
			mulss		xmm1, [edi+0]			// lptr[-3*nc+0] * s1
			subss		xmm0, xmm1
			movss		[esi-16], xmm0			// xptr[-4] = s0
			jl			done4rows
			sub			edi, edx
			sub			edi, 16
			sub			esi, 16
			jmp			process4rows
		done4rows:
			pop			ebx
		}
	}

	// process left over rows
	for ( i = (m&3)-1; i >= 0; i-- ) {
		s0 = b[i];
		lptr = L[0] + i;
		for ( j = i + 1; j < n; j++ ) {
			s0 -= lptr[j*nc] * x[j];
		}
		x[i] = s0;
	}

#else

	int i, j, m;
	double s0, s1, s2, s3, t;
	const float *lptr2;
	float *xptr, *xptr2;

	m = n;
	if ( m & 1 ) {

		m--;
		x[m] = b[m];

		lptr = L.toFloatPtr() + m * nc + m - 4;
		xptr = x + m;
		// process 4 rows at a time
		for ( i = m; i >= 4; i -= 4 ) {
			s0 = b[i-4];
			s1 = b[i-3];
			s2 = b[i-2];
			s3 = b[i-1];
			// process 4x4 blocks
			xptr2 = xptr;	// x + i;
			lptr2 = lptr;	// ptr = L[i] + i - 4;
			for ( j = 0; j < m-i; j += 4 ) {
				t = xptr2[0];
				s0 -= lptr2[0] * t;
				s1 -= lptr2[1] * t;
				s2 -= lptr2[2] * t;
				s3 -= lptr2[3] * t;
				lptr2 += nc;
				xptr2++;
				t = xptr2[0];
				s0 -= lptr2[0] * t;
				s1 -= lptr2[1] * t;
				s2 -= lptr2[2] * t;
				s3 -= lptr2[3] * t;
				lptr2 += nc;
				xptr2++;
				t = xptr2[0];
				s0 -= lptr2[0] * t;
				s1 -= lptr2[1] * t;
				s2 -= lptr2[2] * t;
				s3 -= lptr2[3] * t;
				lptr2 += nc;
				xptr2++;
				t = xptr2[0];
				s0 -= lptr2[0] * t;
				s1 -= lptr2[1] * t;
				s2 -= lptr2[2] * t;
				s3 -= lptr2[3] * t;
				lptr2 += nc;
				xptr2++;
			}
			t = xptr2[0];
			s0 -= lptr2[0] * t;
			s1 -= lptr2[1] * t;
			s2 -= lptr2[2] * t;
			s3 -= lptr2[3] * t;
			// process left over of the 4 rows
			lptr -= nc;
			s0 -= lptr[0] * s3;
			s1 -= lptr[1] * s3;
			s2 -= lptr[2] * s3;
			lptr -= nc;
			s0 -= lptr[0] * s2;
			s1 -= lptr[1] * s2;
			lptr -= nc;
			s0 -= lptr[0] * s1;
			lptr -= nc;
			// store result
			xptr[-4] = s0;
			xptr[-3] = s1;
			xptr[-2] = s2;
			xptr[-1] = s3;
			// update pointers for next four rows
			lptr -= 4;
			xptr -= 4;
		}

	} else {

		lptr = L.toFloatPtr() + m * nc + m - 4;
		xptr = x + m;
		// process 4 rows at a time
		for ( i = m; i >= 4; i -= 4 ) {
			s0 = b[i-4];
			s1 = b[i-3];
			s2 = b[i-2];
			s3 = b[i-1];
			// process 4x4 blocks
			xptr2 = xptr;	// x + i;
			lptr2 = lptr;	// ptr = L[i] + i - 4;
			for ( j = 0; j < m-i; j += 4 ) {
				t = xptr2[0];
				s0 -= lptr2[0] * t;
				s1 -= lptr2[1] * t;
				s2 -= lptr2[2] * t;
				s3 -= lptr2[3] * t;
				lptr2 += nc;
				xptr2++;
				t = xptr2[0];
				s0 -= lptr2[0] * t;
				s1 -= lptr2[1] * t;
				s2 -= lptr2[2] * t;
				s3 -= lptr2[3] * t;
				lptr2 += nc;
				xptr2++;
				t = xptr2[0];
				s0 -= lptr2[0] * t;
				s1 -= lptr2[1] * t;
				s2 -= lptr2[2] * t;
				s3 -= lptr2[3] * t;
				lptr2 += nc;
				xptr2++;
				t = xptr2[0];
				s0 -= lptr2[0] * t;
				s1 -= lptr2[1] * t;
				s2 -= lptr2[2] * t;
				s3 -= lptr2[3] * t;
				lptr2 += nc;
				xptr2++;
			}
			// process left over of the 4 rows
			lptr -= nc;
			s0 -= lptr[0] * s3;
			s1 -= lptr[1] * s3;
			s2 -= lptr[2] * s3;
			lptr -= nc;
			s0 -= lptr[0] * s2;
			s1 -= lptr[1] * s2;
			lptr -= nc;
			s0 -= lptr[0] * s1;
			lptr -= nc;
			// store result
			xptr[-4] = s0;
			xptr[-3] = s1;
			xptr[-2] = s2;
			xptr[-1] = s3;
			// update pointers for next four rows
			lptr -= 4;
			xptr -= 4;
		}
	}
	// process left over rows
	for ( i--; i >= 0; i-- ) {
		s0 = b[i];
		lptr = L[0] + i;
		for ( j = i + 1; j < m; j++ ) {
			s0 -= lptr[j*nc] * x[j];
		}
		x[i] = s0;
	}

#endif
}

/*
============
CSIMD_SSE::matX_LDLTFactor

  in-place factorization LDL' of the n * n sub-matrix of mat
  the reciprocal of the diagonal elements are stored in invDiag
  currently assumes the number of columns of mat is a multiple of 4
============
*/
bool VPCALL CSIMD_SSE::matX_LDLTFactor( CMatXD &mat, CVecXD &invDiag, const int n ) {
#if 1

	int j, nc;
	float *v, *diag, *invDiagPtr, *mptr;
	double s0, s1, s2, sum, d;

	v = (float *) _alloca16( n * sizeof( float ) );
	diag = (float *) _alloca16( n * sizeof( float ) );
	invDiagPtr = invDiag.toFloatPtr();

	nc = mat.getNumColumns();

	SMF_ASSERT( ( nc & 3 ) == 0 );

	if ( n <= 0 ) {
		return true;
	}

	mptr = mat[0];

	sum = mptr[0];

	if ( sum == 0.0f ) {
		return false;
	}

	diag[0] = sum;
	invDiagPtr[0] = d = 1.0f / sum;

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
	invDiagPtr[1] = d = 1.0f / sum;

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
	invDiagPtr[2] = d = 1.0f / sum;

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
	invDiagPtr[3] = d = 1.0f / sum;

	if ( n <= 4 ) {
		return true;
	}

	mptr = mat[0];
	for ( j = 4; j < n; j++ ) {
		mptr[j*nc+3] = ( mptr[j*nc+3] - v[0] * mptr[j*nc+0] - v[1] * mptr[j*nc+1] - v[2] * mptr[j*nc+2] ) * d;
	}

	int ncf = nc * sizeof( float );
	mptr = mat[0];

	__asm {
		xorps		xmm2, xmm2
		xorps		xmm3, xmm3
		xorps		xmm4, xmm4

		push		ebx
		mov			ebx, 4

	loopRow:
			cmp			ebx, n
			jge			done

			mov			ecx, ebx				// esi = i
			shl			ecx, 2					// esi = i * 4
			mov			edx, diag				// edx = diag
			add			edx, ecx				// edx = &diag[i]
			mov			edi, ebx				// edi = i
			imul		edi, ncf				// edi = i * nc * sizeof( float )
			add			edi, mptr				// edi = mat[i]
			add			edi, ecx				// edi = &mat[i][i]
			mov			esi, v					// ecx = v
			add			esi, ecx				// ecx = &v[i]
			mov			eax, invDiagPtr			// eax = invDiagPtr
			add			eax, ecx				// eax = &invDiagPtr[i]
			neg			ecx

			movaps		xmm0, [edx+ecx]
			mulps		xmm0, [edi+ecx]
			movaps		[esi+ecx], xmm0
			mulps		xmm0, [edi+ecx]
			add			ecx, 12*4
			jg			doneDot8
		dot8:
			movaps		xmm1, [edx+ecx-(8*4)]
			mulps		xmm1, [edi+ecx-(8*4)]
			movaps		[esi+ecx-(8*4)], xmm1
			mulps		xmm1, [edi+ecx-(8*4)]
			addps		xmm0, xmm1
			movaps		xmm2, [edx+ecx-(4*4)]
			mulps		xmm2, [edi+ecx-(4*4)]
			movaps		[esi+ecx-(4*4)], xmm2
			mulps		xmm2, [edi+ecx-(4*4)]
			addps		xmm0, xmm2
			add			ecx, 8*4
			jle			dot8
		doneDot8:
			sub			ecx, 4*4
			jg			doneDot4
			movaps		xmm1, [edx+ecx-(4*4)]
			mulps		xmm1, [edi+ecx-(4*4)]
			movaps		[esi+ecx-(4*4)], xmm1
			mulps		xmm1, [edi+ecx-(4*4)]
			addps		xmm0, xmm1
			add			ecx, 4*4
		doneDot4:
			sub			ecx, 2*4
			jg			doneDot2
			movlps		xmm3, [edx+ecx-(2*4)]
			movlps		xmm4, [edi+ecx-(2*4)]
			mulps		xmm3, xmm4
			movlps		[esi+ecx-(2*4)], xmm3
			mulps		xmm3, xmm4
			addps		xmm0, xmm3
			add			ecx, 2*4
		doneDot2:
			sub			ecx, 1*4
			jg			doneDot1
			movss		xmm3, [edx+ecx-(1*4)]
			movss		xmm4, [edi+ecx-(1*4)]
			mulss		xmm3, xmm4
			movss		[esi+ecx-(1*4)], xmm3
			mulss		xmm3, xmm4
			addss		xmm0, xmm3
		doneDot1:
			movhlps		xmm2, xmm0
			addps		xmm0, xmm2
			movaps		xmm2, xmm0
			shufps		xmm2, xmm2, R_SHUFFLEPS( 1, 0, 0, 0 )
			addss		xmm0, xmm2
			movss		xmm1, [edi]
			subss		xmm1, xmm0
			movss		[edi], xmm1				// mptr[i] = sum;
			movss		[edx], xmm1				// diag[i] = sum;

			// if ( sum == 0.0f ) return false;
			movaps		xmm2, xmm1
			cmpeqss		xmm2, SIMD_SP_zero
			andps		xmm2, SIMD_SP_tiny
			orps		xmm1, xmm2

			rcpss		xmm7, xmm1
			mulss		xmm1, xmm7
			mulss		xmm1, xmm7
			addss		xmm7, xmm7
			subss		xmm7, xmm1
			movss		[eax], xmm7				// invDiagPtr[i] = 1.0f / sum;

			mov			edx, n					// edx = n
			sub			edx, ebx				// edx = n - i
			dec			edx						// edx = n - i - 1
			jle			doneSubRow				// if ( i + 1 >= n ) return true;

			mov			eax, ebx				// eax = i
			shl			eax, 2					// eax = i * 4
			neg			eax

		loopSubRow:
				add			edi, ncf
				mov			ecx, eax
				movaps		xmm0, [esi+ecx]
				mulps		xmm0, [edi+ecx]
				add			ecx, 12*4
				jg			doneSubDot8
			subDot8:
				movaps		xmm1, [esi+ecx-(8*4)]
				mulps		xmm1, [edi+ecx-(8*4)]
				addps		xmm0, xmm1
				movaps		xmm2, [esi+ecx-(4*4)]
				mulps		xmm2, [edi+ecx-(4*4)]
				addps		xmm0, xmm2
				add			ecx, 8*4
				jle			subDot8
			doneSubDot8:
				sub			ecx, 4*4
				jg			doneSubDot4
				movaps		xmm1, [esi+ecx-(4*4)]
				mulps		xmm1, [edi+ecx-(4*4)]
				addps		xmm0, xmm1
				add			ecx, 4*4
			doneSubDot4:
				sub			ecx, 2*4
				jg			doneSubDot2
				movlps		xmm3, [esi+ecx-(2*4)]
				movlps		xmm4, [edi+ecx-(2*4)]
				mulps		xmm3, xmm4
				addps		xmm0, xmm3
				add			ecx, 2*4
			doneSubDot2:
				sub			ecx, 1*4
				jg			doneSubDot1
				movss		xmm3, [esi+ecx-(1*4)]
				movss		xmm4, [edi+ecx-(1*4)]
				mulss		xmm3, xmm4
				addss		xmm0, xmm3
			doneSubDot1:
				movhlps		xmm2, xmm0
				addps		xmm0, xmm2
				movaps		xmm2, xmm0
				shufps		xmm2, xmm2, R_SHUFFLEPS( 1, 0, 0, 0 )
				addss		xmm0, xmm2
				movss		xmm1, [edi]
				subss		xmm1, xmm0
				mulss		xmm1, xmm7
				movss		[edi], xmm1
				dec			edx
				jg			loopSubRow
		doneSubRow:
			inc		ebx
			jmp		loopRow
	done:
		pop		ebx
	}

	return true;

#else

	int i, j, k, nc;
	float *v, *diag, *mptr;
	double s0, s1, s2, s3, sum, d;

	v = (float *) _alloca16( n * sizeof( float ) );
	diag = (float *) _alloca16( n * sizeof( float ) );

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
			case 3: v[k+2] = diag[k+2] * mptr[k+2]; s0 += v[k+2] * mptr[k+2];
			case 2: v[k+1] = diag[k+1] * mptr[k+1]; s1 += v[k+1] * mptr[k+1];
			case 1: v[k+0] = diag[k+0] * mptr[k+0]; s2 += v[k+0] * mptr[k+0];
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
				case 7: s0 += mptr[k+6] * v[k+6];
				case 6: s1 += mptr[k+5] * v[k+5];
				case 5: s2 += mptr[k+4] * v[k+4];
				case 4: s3 += mptr[k+3] * v[k+3];
				case 3: s0 += mptr[k+2] * v[k+2];
				case 2: s1 += mptr[k+1] * v[k+1];
				case 1: s2 += mptr[k+0] * v[k+0];
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

#endif
}

/*
============
CSIMD_SSE::blendJoints
============
*/
#define REFINE_BLENDJOINTS_RECIPROCAL

void VPCALL CSIMD_SSE::blendJoints( CJointQuaternion *joints, const CJointQuaternion *blendJoints, const float lerp, const int *index, const int numJoints ) {
	int i;

	if ( lerp <= 0.0f ) {
		return;
	} else if ( lerp >= 1.0f ) {
		for ( i = 0; i < numJoints; i++ ) {
			int j = index[i];
			joints[j] = blendJoints[j];
		}
		return;
	}

	for ( i = 0; i <= numJoints - 4; i += 4 ) {
		ALIGN16( float jointVert0[4] );
		ALIGN16( float jointVert1[4] );
		ALIGN16( float jointVert2[4] );
		ALIGN16( float blendVert0[4] );
		ALIGN16( float blendVert1[4] );
		ALIGN16( float blendVert2[4] );
		ALIGN16( float jointQuat0[4] );
		ALIGN16( float jointQuat1[4] );
		ALIGN16( float jointQuat2[4] );
		ALIGN16( float jointQuat3[4] );
		ALIGN16( float blendQuat0[4] );
		ALIGN16( float blendQuat1[4] );
		ALIGN16( float blendQuat2[4] );
		ALIGN16( float blendQuat3[4] );

		for ( int j = 0; j < 4; j++ ) {
			int n = index[i+j];

			jointVert0[j] = joints[n].t[0];
			jointVert1[j] = joints[n].t[1];
			jointVert2[j] = joints[n].t[2];

			blendVert0[j] = blendJoints[n].t[0];
			blendVert1[j] = blendJoints[n].t[1];
			blendVert2[j] = blendJoints[n].t[2];

			jointQuat0[j] = joints[n].q[0];
			jointQuat1[j] = joints[n].q[1];
			jointQuat2[j] = joints[n].q[2];
			jointQuat3[j] = joints[n].q[3];

			blendQuat0[j] = blendJoints[n].q[0];
			blendQuat1[j] = blendJoints[n].q[1];
			blendQuat2[j] = blendJoints[n].q[2];
			blendQuat3[j] = blendJoints[n].q[3];
		}

#if 1
		__asm {
			// lerp translation
			movss		xmm7, lerp
			shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
			movaps		xmm0, blendVert0
			subps		xmm0, jointVert0
			mulps		xmm0, xmm7
			addps		xmm0, jointVert0
			movaps		jointVert0, xmm0
			movaps		xmm1, blendVert1
			subps		xmm1, jointVert1
			mulps		xmm1, xmm7
			addps		xmm1, jointVert1
			movaps		jointVert1, xmm1
			movaps		xmm2, blendVert2
			subps		xmm2, jointVert2
			mulps		xmm2, xmm7
			addps		xmm2, jointVert2
			movaps		jointVert2, xmm2

			// lerp quaternions
			movaps		xmm0, jointQuat0
			mulps		xmm0, blendQuat0
			movaps		xmm1, jointQuat1
			mulps		xmm1, blendQuat1
			addps		xmm0, xmm1
			movaps		xmm2, jointQuat2
			mulps		xmm2, blendQuat2
			addps		xmm0, xmm2
			movaps		xmm3, jointQuat3
			mulps		xmm3, blendQuat3
			addps		xmm0, xmm3					// xmm0 = cosom

			movaps		xmm1, xmm0
			movaps		xmm2, xmm0
			andps		xmm1, SIMD_SP_signBitMask	// xmm1 = signBit
			xorps		xmm0, xmm1
			mulps		xmm2, xmm2

			xorps		xmm4, xmm4
			movaps		xmm3, SIMD_SP_one
			subps		xmm3, xmm2					// xmm3 = scale0
			cmpeqps		xmm4, xmm3
			andps		xmm4, SIMD_SP_tiny			// if values are zero replace them with a tiny number
			andps		xmm3, SIMD_SP_absMask		// make sure the values are positive
			orps		xmm3, xmm4

#ifdef REFINE_BLENDJOINTS_RECIPROCAL
			movaps		xmm2, xmm3
			rsqrtps		xmm4, xmm2
			mulps		xmm2, xmm4
			mulps		xmm2, xmm4
			subps		xmm2, SIMD_SP_rsqrt_c0
			mulps		xmm4, SIMD_SP_rsqrt_c1
			mulps		xmm2, xmm4
#else
			rsqrtps		xmm2, xmm3					// xmm2 = sinom
#endif
			mulps		xmm3, xmm2					// xmm3 = sqrt( scale0 )

			// omega0 = atan2( xmm3, xmm0 )
			movaps		xmm4, xmm0
			minps		xmm0, xmm3
			maxps		xmm3, xmm4
			cmpeqps		xmm4, xmm0

#ifdef REFINE_BLENDJOINTS_RECIPROCAL
			rcpps		xmm5, xmm3
			mulps		xmm3, xmm5
			mulps		xmm3, xmm5
			addps		xmm5, xmm5
			subps		xmm5, xmm3					// xmm5 = 1 / y or 1 / x
			mulps		xmm0, xmm5					// xmm0 = x / y or y / x
#else
			rcpps		xmm3, xmm3					// xmm3 = 1 / y or 1 / x
			mulps		xmm0, xmm3					// xmm0 = x / y or y / x
#endif
			movaps		xmm3, xmm4
			andps		xmm3, SIMD_SP_signBitMask
			xorps		xmm0, xmm3					// xmm0 = -x / y or y / x
			andps		xmm4, SIMD_SP_halfPI		// xmm4 = HALF_PI or 0.0f
			movaps		xmm3, xmm0
			mulps		xmm3, xmm3					// xmm3 = s
			movaps		xmm5, SIMD_SP_atan_c0
			mulps		xmm5, xmm3
			addps		xmm5, SIMD_SP_atan_c1
			mulps		xmm5, xmm3
			addps		xmm5, SIMD_SP_atan_c2
			mulps		xmm5, xmm3
			addps		xmm5, SIMD_SP_atan_c3
			mulps		xmm5, xmm3
			addps		xmm5, SIMD_SP_atan_c4
			mulps		xmm5, xmm3
			addps		xmm5, SIMD_SP_atan_c5
			mulps		xmm5, xmm3
			addps		xmm5, SIMD_SP_atan_c6
			mulps		xmm5, xmm3
			addps		xmm5, SIMD_SP_atan_c7
			mulps		xmm5, xmm3
			addps		xmm5, SIMD_SP_one
			mulps		xmm5, xmm0
			addps		xmm5, xmm4					// xmm5 = omega0

			movaps		xmm6, xmm7					// xmm6 = lerp
			mulps		xmm6, xmm5					// xmm6 = omega1
			subps		xmm5, xmm6					// xmm5 = omega0

			// scale0 = sin( xmm5 ) * xmm2
			// scale1 = sin( xmm6 ) * xmm2
			movaps		xmm3, xmm5
			movaps		xmm7, xmm6
			mulps		xmm3, xmm3
			mulps		xmm7, xmm7
			movaps		xmm4, SIMD_SP_sin_c0
			movaps		xmm0, SIMD_SP_sin_c0
			mulps		xmm4, xmm3
			mulps		xmm0, xmm7
			addps		xmm4, SIMD_SP_sin_c1
			addps		xmm0, SIMD_SP_sin_c1
			mulps		xmm4, xmm3
			mulps		xmm0, xmm7
			addps		xmm4, SIMD_SP_sin_c2
			addps		xmm0, SIMD_SP_sin_c2
			mulps		xmm4, xmm3
			mulps		xmm0, xmm7
			addps		xmm4, SIMD_SP_sin_c3
			addps		xmm0, SIMD_SP_sin_c3
			mulps		xmm4, xmm3
			mulps		xmm0, xmm7
			addps		xmm4, SIMD_SP_sin_c4
			addps		xmm0, SIMD_SP_sin_c4
			mulps		xmm4, xmm3
			mulps		xmm0, xmm7
			addps		xmm4, SIMD_SP_one
			addps		xmm0, SIMD_SP_one
			mulps		xmm5, xmm4
			mulps		xmm6, xmm0
			mulps		xmm5, xmm2					// xmm5 = scale0
			mulps		xmm6, xmm2					// xmm6 = scale1

			xorps		xmm6, xmm1

			movaps		xmm0, jointQuat0
			mulps		xmm0, xmm5
			movaps		xmm1, blendQuat0
			mulps		xmm1, xmm6
			addps		xmm0, xmm1
			movaps		jointQuat0, xmm0

			movaps		xmm1, jointQuat1
			mulps		xmm1, xmm5
			movaps		xmm2, blendQuat1
			mulps		xmm2, xmm6
			addps		xmm1, xmm2
			movaps		jointQuat1, xmm1

			movaps		xmm2, jointQuat2
			mulps		xmm2, xmm5
			movaps		xmm3, blendQuat2
			mulps		xmm3, xmm6
			addps		xmm2, xmm3
			movaps		jointQuat2, xmm2

			movaps		xmm3, jointQuat3
			mulps		xmm3, xmm5
			movaps		xmm4, blendQuat3
			mulps		xmm4, xmm6
			addps		xmm3, xmm4
			movaps		jointQuat3, xmm3
		}

#else

		jointVert0[0] += lerp * ( blendVert0[0] - jointVert0[0] );
		jointVert0[1] += lerp * ( blendVert0[1] - jointVert0[1] );
		jointVert0[2] += lerp * ( blendVert0[2] - jointVert0[2] );
		jointVert0[3] += lerp * ( blendVert0[3] - jointVert0[3] );

		jointVert1[0] += lerp * ( blendVert1[0] - jointVert1[0] );
		jointVert1[1] += lerp * ( blendVert1[1] - jointVert1[1] );
		jointVert1[2] += lerp * ( blendVert1[2] - jointVert1[2] );
		jointVert1[3] += lerp * ( blendVert1[3] - jointVert1[3] );

		jointVert2[0] += lerp * ( blendVert2[0] - jointVert2[0] );
		jointVert2[1] += lerp * ( blendVert2[1] - jointVert2[1] );
		jointVert2[2] += lerp * ( blendVert2[2] - jointVert2[2] );
		jointVert2[3] += lerp * ( blendVert2[3] - jointVert2[3] );

		ALIGN16( float cosom[4] );
		ALIGN16( float sinom[4] );
		ALIGN16( float omega0[4] );
		ALIGN16( float omega1[4] );
		ALIGN16( float scale0[4] );
		ALIGN16( float scale1[4] );
		ALIGN16( unsigned long signBit[4] );

		cosom[0] = jointQuat0[0] * blendQuat0[0];
		cosom[1] = jointQuat0[1] * blendQuat0[1];
		cosom[2] = jointQuat0[2] * blendQuat0[2];
		cosom[3] = jointQuat0[3] * blendQuat0[3];

		cosom[0] += jointQuat1[0] * blendQuat1[0];
		cosom[1] += jointQuat1[1] * blendQuat1[1];
		cosom[2] += jointQuat1[2] * blendQuat1[2];
		cosom[3] += jointQuat1[3] * blendQuat1[3];

		cosom[0] += jointQuat2[0] * blendQuat2[0];
		cosom[1] += jointQuat2[1] * blendQuat2[1];
		cosom[2] += jointQuat2[2] * blendQuat2[2];
		cosom[3] += jointQuat2[3] * blendQuat2[3];

		cosom[0] += jointQuat3[0] * blendQuat3[0];
		cosom[1] += jointQuat3[1] * blendQuat3[1];
		cosom[2] += jointQuat3[2] * blendQuat3[2];
		cosom[3] += jointQuat3[3] * blendQuat3[3];

		signBit[0] = (*(unsigned long *)&cosom[0]) & ( 1 << 31 );
		signBit[1] = (*(unsigned long *)&cosom[1]) & ( 1 << 31 );
		signBit[2] = (*(unsigned long *)&cosom[2]) & ( 1 << 31 );
		signBit[3] = (*(unsigned long *)&cosom[3]) & ( 1 << 31 );

		(*(unsigned long *)&cosom[0]) ^= signBit[0];
		(*(unsigned long *)&cosom[1]) ^= signBit[1];
		(*(unsigned long *)&cosom[2]) ^= signBit[2];
		(*(unsigned long *)&cosom[3]) ^= signBit[3];

		scale0[0] = 1.0f - cosom[0] * cosom[0];
		scale0[1] = 1.0f - cosom[1] * cosom[1];
		scale0[2] = 1.0f - cosom[2] * cosom[2];
		scale0[3] = 1.0f - cosom[3] * cosom[3];

		scale0[0] = ( scale0[0] <= 0.0f ) ? SIMD_SP_tiny[0] : scale0[0];
		scale0[1] = ( scale0[1] <= 0.0f ) ? SIMD_SP_tiny[1] : scale0[1];
		scale0[2] = ( scale0[2] <= 0.0f ) ? SIMD_SP_tiny[2] : scale0[2];
		scale0[3] = ( scale0[3] <= 0.0f ) ? SIMD_SP_tiny[3] : scale0[3];

		sinom[0] = CMath::rSqrt( scale0[0] );
		sinom[1] = CMath::rSqrt( scale0[1] );
		sinom[2] = CMath::rSqrt( scale0[2] );
		sinom[3] = CMath::rSqrt( scale0[3] );

		scale0[0] *= sinom[0];
		scale0[1] *= sinom[1];
		scale0[2] *= sinom[2];
		scale0[3] *= sinom[3];

		omega0[0] = SSE_ATanPositive( scale0[0], cosom[0] );
		omega0[1] = SSE_ATanPositive( scale0[1], cosom[1] );
		omega0[2] = SSE_ATanPositive( scale0[2], cosom[2] );
		omega0[3] = SSE_ATanPositive( scale0[3], cosom[3] );

		omega1[0] = lerp * omega0[0];
		omega1[1] = lerp * omega0[1];
		omega1[2] = lerp * omega0[2];
		omega1[3] = lerp * omega0[3];

		omega0[0] -= omega1[0];
		omega0[1] -= omega1[1];
		omega0[2] -= omega1[2];
		omega0[3] -= omega1[3];

		scale0[0] = SSE_SinZeroHalfPI( omega0[0] ) * sinom[0];
		scale0[1] = SSE_SinZeroHalfPI( omega0[1] ) * sinom[1];
		scale0[2] = SSE_SinZeroHalfPI( omega0[2] ) * sinom[2];
		scale0[3] = SSE_SinZeroHalfPI( omega0[3] ) * sinom[3];

		scale1[0] = SSE_SinZeroHalfPI( omega1[0] ) * sinom[0];
		scale1[1] = SSE_SinZeroHalfPI( omega1[1] ) * sinom[1];
		scale1[2] = SSE_SinZeroHalfPI( omega1[2] ) * sinom[2];
		scale1[3] = SSE_SinZeroHalfPI( omega1[3] ) * sinom[3];

		(*(unsigned long *)&scale1[0]) ^= signBit[0];
		(*(unsigned long *)&scale1[1]) ^= signBit[1];
		(*(unsigned long *)&scale1[2]) ^= signBit[2];
		(*(unsigned long *)&scale1[3]) ^= signBit[3];

		jointQuat0[0] = scale0[0] * jointQuat0[0] + scale1[0] * blendQuat0[0];
		jointQuat0[1] = scale0[1] * jointQuat0[1] + scale1[1] * blendQuat0[1];
		jointQuat0[2] = scale0[2] * jointQuat0[2] + scale1[2] * blendQuat0[2];
		jointQuat0[3] = scale0[3] * jointQuat0[3] + scale1[3] * blendQuat0[3];

		jointQuat1[0] = scale0[0] * jointQuat1[0] + scale1[0] * blendQuat1[0];
		jointQuat1[1] = scale0[1] * jointQuat1[1] + scale1[1] * blendQuat1[1];
		jointQuat1[2] = scale0[2] * jointQuat1[2] + scale1[2] * blendQuat1[2];
		jointQuat1[3] = scale0[3] * jointQuat1[3] + scale1[3] * blendQuat1[3];

		jointQuat2[0] = scale0[0] * jointQuat2[0] + scale1[0] * blendQuat2[0];
		jointQuat2[1] = scale0[1] * jointQuat2[1] + scale1[1] * blendQuat2[1];
		jointQuat2[2] = scale0[2] * jointQuat2[2] + scale1[2] * blendQuat2[2];
		jointQuat2[3] = scale0[3] * jointQuat2[3] + scale1[3] * blendQuat2[3];

		jointQuat3[0] = scale0[0] * jointQuat3[0] + scale1[0] * blendQuat3[0];
		jointQuat3[1] = scale0[1] * jointQuat3[1] + scale1[1] * blendQuat3[1];
		jointQuat3[2] = scale0[2] * jointQuat3[2] + scale1[2] * blendQuat3[2];
		jointQuat3[3] = scale0[3] * jointQuat3[3] + scale1[3] * blendQuat3[3];

#endif

		for ( int j = 0; j < 4; j++ ) {
			int n = index[i+j];

			joints[n].t[0] = jointVert0[j];
			joints[n].t[1] = jointVert1[j];
			joints[n].t[2] = jointVert2[j];

			joints[n].q[0] = jointQuat0[j];
			joints[n].q[1] = jointQuat1[j];
			joints[n].q[2] = jointQuat2[j];
			joints[n].q[3] = jointQuat3[j];
		}
	}

	for ( ; i < numJoints; i++ ) {
		int n = index[i];

		CVec3D &jointVert = joints[n].t;
		const CVec3D &blendVert = blendJoints[n].t;

		jointVert[0] += lerp * ( blendVert[0] - jointVert[0] );
		jointVert[1] += lerp * ( blendVert[1] - jointVert[1] );
		jointVert[2] += lerp * ( blendVert[2] - jointVert[2] );

		CQuaternion &jointQuat = joints[n].q;
		const CQuaternion &blendQuat = blendJoints[n].q;

		float cosom;
		float sinom;
		float omega;
		float scale0;
		float scale1;
		unsigned long signBit;

		cosom = jointQuat.x * blendQuat.x + jointQuat.y * blendQuat.y + jointQuat.z * blendQuat.z + jointQuat.w * blendQuat.w;

		signBit = (*(unsigned long *)&cosom) & ( 1 << 31 );

		(*(unsigned long *)&cosom) ^= signBit;

		scale0 = 1.0f - cosom * cosom;
		scale0 = ( scale0 <= 0.0f ) ? SIMD_SP_tiny[0] : scale0;
		sinom = CMath::invSqrt( scale0 );
		omega = CMath::atan16( scale0 * sinom, cosom );
		scale0 = CMath::sin16( ( 1.0f - lerp ) * omega ) * sinom;
		scale1 = CMath::sin16( lerp * omega ) * sinom;

		(*(unsigned long *)&scale1) ^= signBit;

		jointQuat.x = scale0 * jointQuat.x + scale1 * blendQuat.x;
		jointQuat.y = scale0 * jointQuat.y + scale1 * blendQuat.y;
		jointQuat.z = scale0 * jointQuat.z + scale1 * blendQuat.z;
		jointQuat.w = scale0 * jointQuat.w + scale1 * blendQuat.w;
	}
}

/*
============
CSIMD_SSE::convertJointQuatsToJointMats
============
*/
void VPCALL CSIMD_SSE::convertJointQuatsToJointMats( CMatJoint3x4 *jointMats, const CJointQuaternion *jointQuats, const int numJoints ) {

	//size_t jq = sizeof( CJointQuaternion );
	//int jq2 = JOINTQUAT_SIZE;
	SMF_ASSERT( sizeof( CJointQuaternion ) == JOINTQUAT_SIZE );
	SMF_ASSERT( sizeof( CMatJoint3x4 ) == JOINTMAT_SIZE );
	SMF_ASSERT( (int)(&((CJointQuaternion *)0)->t) == (int)(&((CJointQuaternion *)0)->q) + (int)sizeof( ((CJointQuaternion *)0)->q ) );

	for ( int i = 0; i < numJoints; i++ ) {

		const float *q = jointQuats[i].q.toFloatPtr();
		float *m = jointMats[i].toFloatPtr();

		m[0*4+3] = q[4];
		m[1*4+3] = q[5];
		m[2*4+3] = q[6];

		float x2 = q[0] + q[0];
		float y2 = q[1] + q[1];
		float z2 = q[2] + q[2];

		{
			float xx = q[0] * x2;
			float yy = q[1] * y2;
			float zz = q[2] * z2;

			m[0*4+0] = 1.0f - yy - zz;
			m[1*4+1] = 1.0f - xx - zz;
			m[2*4+2] = 1.0f - xx - yy;
		}

		{
			float yz = q[1] * z2;
			float wx = q[3] * x2;

			m[2*4+1] = yz - wx;
			m[1*4+2] = yz + wx;
		}

		{
			float xy = q[0] * y2;
			float wz = q[3] * z2;

			m[1*4+0] = xy - wz;
			m[0*4+1] = xy + wz;
		}

		{
			float xz = q[0] * z2;
			float wy = q[3] * y2;

			m[0*4+2] = xz - wy;
			m[2*4+0] = xz + wy;
		}
	}
}

/*
============
CSIMD_SSE::convertJointMatsToJointQuats
============
*/
void VPCALL CSIMD_SSE::convertJointMatsToJointQuats( CJointQuaternion *jointQuats, const CMatJoint3x4 *jointMats, const int numJoints ) {

	SMF_ASSERT( sizeof( CJointQuaternion ) == JOINTQUAT_SIZE );
	SMF_ASSERT( sizeof( CMatJoint3x4 ) == JOINTMAT_SIZE );
	SMF_ASSERT( (int)(&((CJointQuaternion *)0)->t) == (int)(&((CJointQuaternion *)0)->q) + (int)sizeof( ((CJointQuaternion *)0)->q ) );

#if 1

	ALIGN16( sf_u8 shuffle[16] );

	__asm {
		mov			eax, numJoints
		mov			esi, jointMats
		mov			edi, jointQuats
		and			eax, ~3
		jz			done4
		imul		eax, JOINTMAT_SIZE
		add			esi, eax
		neg			eax

	loopMat4:
		movss		xmm5, [esi+eax+3*JOINTMAT_SIZE+0*16+0*4]
		movss		xmm6, [esi+eax+3*JOINTMAT_SIZE+1*16+1*4]
		movss		xmm7, [esi+eax+3*JOINTMAT_SIZE+2*16+2*4]

		shufps		xmm5, xmm5, R_SHUFFLEPS( 3, 0, 1, 2 )
		shufps		xmm6, xmm6, R_SHUFFLEPS( 3, 0, 1, 2 )
		shufps		xmm7, xmm7, R_SHUFFLEPS( 3, 0, 1, 2 )

		movss		xmm0, [esi+eax+2*JOINTMAT_SIZE+0*16+0*4]
		movss		xmm1, [esi+eax+2*JOINTMAT_SIZE+1*16+1*4]
		movss		xmm2, [esi+eax+2*JOINTMAT_SIZE+2*16+2*4]

		movss		xmm5, xmm0
		movss		xmm6, xmm1
		movss		xmm7, xmm2

		shufps		xmm5, xmm5, R_SHUFFLEPS( 3, 0, 1, 2 )
		shufps		xmm6, xmm6, R_SHUFFLEPS( 3, 0, 1, 2 )
		shufps		xmm7, xmm7, R_SHUFFLEPS( 3, 0, 1, 2 )

		movss		xmm0, [esi+eax+1*JOINTMAT_SIZE+0*16+0*4]
		movss		xmm1, [esi+eax+1*JOINTMAT_SIZE+1*16+1*4]
		movss		xmm2, [esi+eax+1*JOINTMAT_SIZE+2*16+2*4]

		movss		xmm5, xmm0
		movss		xmm6, xmm1
		movss		xmm7, xmm2

		shufps		xmm5, xmm5, R_SHUFFLEPS( 3, 0, 1, 2 )
		shufps		xmm6, xmm6, R_SHUFFLEPS( 3, 0, 1, 2 )
		shufps		xmm7, xmm7, R_SHUFFLEPS( 3, 0, 1, 2 )

		movss		xmm0, [esi+eax+0*JOINTMAT_SIZE+0*16+0*4]
		movss		xmm1, [esi+eax+0*JOINTMAT_SIZE+1*16+1*4]
		movss		xmm2, [esi+eax+0*JOINTMAT_SIZE+2*16+2*4]

		movss		xmm5, xmm0
		movss		xmm6, xmm1
		movss		xmm7, xmm2

		// -------------------

		movaps		xmm0, xmm5
		addps		xmm0, xmm6
		addps		xmm0, xmm7
		cmpnltps	xmm0, SIMD_SP_zero						// xmm0 = m[0 * 4 + 0] + m[1 * 4 + 1] + m[2 * 4 + 2] > 0.0f

		movaps		xmm1, xmm5
		movaps		xmm2, xmm5
		cmpnltps	xmm1, xmm6
		cmpnltps	xmm2, xmm7
		andps		xmm2, xmm1								// xmm2 = m[0 * 4 + 0] > m[1 * 4 + 1] && m[0 * 4 + 0] > m[2 * 4 + 2]

		movaps		xmm4, xmm6
		cmpnltps	xmm4, xmm7								// xmm3 = m[1 * 4 + 1] > m[2 * 4 + 2]

		movaps		xmm1, xmm0
		andnps		xmm1, xmm2
		orps		xmm2, xmm0
		movaps		xmm3, xmm2
		andnps		xmm2, xmm4
		orps		xmm3, xmm2
		xorps		xmm3, SIMD_SP_not

		andps		xmm0, SIMD_DW_mat2quatShuffle0
		movaps		xmm4, xmm1
		andps		xmm4, SIMD_DW_mat2quatShuffle1
		orps		xmm0, xmm4
		movaps		xmm4, xmm2
		andps		xmm4, SIMD_DW_mat2quatShuffle2
		orps		xmm0, xmm4
		movaps		xmm4, xmm3
		andps		xmm4, SIMD_DW_mat2quatShuffle3
		orps		xmm4, xmm0

		movaps		shuffle, xmm4

		movaps		xmm0, xmm2
		orps		xmm0, xmm3								// xmm0 = xmm2 | xmm3	= s0
		orps		xmm2, xmm1								// xmm2 = xmm1 | xmm2	= s2
		orps		xmm1, xmm3								// xmm1 = xmm1 | xmm3	= s1

		andps		xmm0, SIMD_SP_signBitMask
		andps		xmm1, SIMD_SP_signBitMask
		andps		xmm2, SIMD_SP_signBitMask

		xorps		xmm5, xmm0
		xorps		xmm6, xmm1
		xorps		xmm7, xmm2
		addps		xmm5, xmm6
		addps		xmm7, SIMD_SP_one
		addps		xmm5, xmm7								// xmm5 = t

		movaps		xmm7, xmm5								// xmm7 = t
		rsqrtps		xmm6, xmm5
		mulps		xmm5, xmm6
		mulps		xmm5, xmm6
		subps		xmm5, SIMD_SP_rsqrt_c0
		mulps		xmm6, SIMD_SP_mat2quat_rsqrt_c1
		mulps		xmm6, xmm5								// xmm5 = s

		mulps		xmm7, xmm6								// xmm7 = s * t
		xorps		xmm6, SIMD_SP_signBitMask				// xmm6 = -s

		// -------------------

		add			edi, 4*JOINTQUAT_SIZE

		movzx		ecx, byte ptr shuffle[0*4+0]			// ecx = k0
		movss		[edi+ecx*4-4*JOINTQUAT_SIZE], xmm7		// q[k0] = s * t;

		movzx		edx, byte ptr shuffle[0*4+1]			// edx = k1
		movss		xmm4, [esi+eax+0*JOINTMAT_SIZE+1*16+0*4]
		xorps		xmm4, xmm2
		subss		xmm4, [esi+eax+0*JOINTMAT_SIZE+0*16+1*4]
		mulss		xmm4, xmm6
		movss		[edi+edx*4-4*JOINTQUAT_SIZE], xmm4		// q[k1] = ( m[0 * 4 + 1] - s2 * m[1 * 4 + 0] ) * s;

		movzx		ecx, byte ptr shuffle[0*4+2]			// ecx = k2
		movss		xmm3, [esi+eax+0*JOINTMAT_SIZE+0*16+2*4]
		xorps		xmm3, xmm1
		subss		xmm3, [esi+eax+0*JOINTMAT_SIZE+2*16+0*4]
		mulss		xmm3, xmm6
		movss		[edi+ecx*4-4*JOINTQUAT_SIZE], xmm3		// q[k2] = ( m[2 * 4 + 0] - s1 * m[0 * 4 + 2] ) * s;

		movzx		edx, byte ptr shuffle[0*4+3]			// edx = k3
		movss		xmm4, [esi+eax+0*JOINTMAT_SIZE+2*16+1*4]
		xorps		xmm4, xmm0
		subss		xmm4, [esi+eax+0*JOINTMAT_SIZE+1*16+2*4]
		mulss		xmm4, xmm6
		movss		[edi+edx*4-4*JOINTQUAT_SIZE], xmm4		// q[k3] = ( m[1 * 4 + 2] - s0 * m[2 * 4 + 1] ) * s;

		mov			ecx, [esi+eax+0*JOINTMAT_SIZE+0*16+3*4]
		mov			[edi-4*JOINTQUAT_SIZE+16], ecx			// q[4] = m[0 * 4 + 3];
		mov			edx, [esi+eax+0*JOINTMAT_SIZE+1*16+3*4]
		mov			[edi-4*JOINTQUAT_SIZE+20], edx			// q[5] = m[1 * 4 + 3];
		mov			ecx, [esi+eax+0*JOINTMAT_SIZE+2*16+3*4]
		mov			[edi-4*JOINTQUAT_SIZE+24], ecx			// q[6] = m[2 * 4 + 3];

		shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm7, xmm7, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm0, xmm0, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm2, xmm2, R_SHUFFLEPS( 1, 2, 3, 0 )

		movzx		ecx, byte ptr shuffle[1*4+0]			// ecx = k0
		movss		[edi+ecx*4-3*JOINTQUAT_SIZE], xmm7		// q[k0] = s * t;

		movzx		edx, byte ptr shuffle[1*4+1]			// edx = k1
		movss		xmm4, [esi+eax+1*JOINTMAT_SIZE+1*16+0*4]
		xorps		xmm4, xmm2
		subss		xmm4, [esi+eax+1*JOINTMAT_SIZE+0*16+1*4]
		mulss		xmm4, xmm6
		movss		[edi+edx*4-3*JOINTQUAT_SIZE], xmm4		// q[k1] = ( m[0 * 4 + 1] - s2 * m[1 * 4 + 0] ) * s;

		movzx		ecx, byte ptr shuffle[1*4+2]			// ecx = k2
		movss		xmm3, [esi+eax+1*JOINTMAT_SIZE+0*16+2*4]
		xorps		xmm3, xmm1
		subss		xmm3, [esi+eax+1*JOINTMAT_SIZE+2*16+0*4]
		mulss		xmm3, xmm6
		movss		[edi+ecx*4-3*JOINTQUAT_SIZE], xmm3		// q[k2] = ( m[2 * 4 + 0] - s1 * m[0 * 4 + 2] ) * s;

		movzx		edx, byte ptr shuffle[1*4+3]			// edx = k3
		movss		xmm4, [esi+eax+1*JOINTMAT_SIZE+2*16+1*4]
		xorps		xmm4, xmm0
		subss		xmm4, [esi+eax+1*JOINTMAT_SIZE+1*16+2*4]
		mulss		xmm4, xmm6
		movss		[edi+edx*4-3*JOINTQUAT_SIZE], xmm4		// q[k3] = ( m[1 * 4 + 2] - s0 * m[2 * 4 + 1] ) * s;

		mov			ecx, [esi+eax+1*JOINTMAT_SIZE+0*16+3*4]
		mov			[edi-3*JOINTQUAT_SIZE+16], ecx			// q[4] = m[0 * 4 + 3];
		mov			edx, [esi+eax+1*JOINTMAT_SIZE+1*16+3*4]
		mov			[edi-3*JOINTQUAT_SIZE+20], edx			// q[5] = m[1 * 4 + 3];
		mov			ecx, [esi+eax+1*JOINTMAT_SIZE+2*16+3*4]
		mov			[edi-3*JOINTQUAT_SIZE+24], ecx			// q[6] = m[2 * 4 + 3];

		shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm7, xmm7, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm0, xmm0, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm2, xmm2, R_SHUFFLEPS( 1, 2, 3, 0 )

		movzx		ecx, byte ptr shuffle[2*4+0]			// ecx = k0
		movss		[edi+ecx*4-2*JOINTQUAT_SIZE], xmm7		// q[k0] = s * t;

		movzx		edx, byte ptr shuffle[2*4+1]			// edx = k1
		movss		xmm4, [esi+eax+2*JOINTMAT_SIZE+1*16+0*4]
		xorps		xmm4, xmm2
		subss		xmm4, [esi+eax+2*JOINTMAT_SIZE+0*16+1*4]
		mulss		xmm4, xmm6
		movss		[edi+edx*4-2*JOINTQUAT_SIZE], xmm4		// q[k1] = ( m[0 * 4 + 1] - s2 * m[1 * 4 + 0] ) * s;

		movzx		ecx, byte ptr shuffle[2*4+2]			// ecx = k2
		movss		xmm3, [esi+eax+2*JOINTMAT_SIZE+0*16+2*4]
		xorps		xmm3, xmm1
		subss		xmm3, [esi+eax+2*JOINTMAT_SIZE+2*16+0*4]
		mulss		xmm3, xmm6
		movss		[edi+ecx*4-2*JOINTQUAT_SIZE], xmm3		// q[k2] = ( m[2 * 4 + 0] - s1 * m[0 * 4 + 2] ) * s;

		movzx		edx, byte ptr shuffle[2*4+3]			// edx = k3
		movss		xmm4, [esi+eax+2*JOINTMAT_SIZE+2*16+1*4]
		xorps		xmm4, xmm0
		subss		xmm4, [esi+eax+2*JOINTMAT_SIZE+1*16+2*4]
		mulss		xmm4, xmm6
		movss		[edi+edx*4-2*JOINTQUAT_SIZE], xmm4		// q[k3] = ( m[1 * 4 + 2] - s0 * m[2 * 4 + 1] ) * s;

		mov			ecx, [esi+eax+2*JOINTMAT_SIZE+0*16+3*4]
		mov			[edi-2*JOINTQUAT_SIZE+16], ecx			// q[4] = m[0 * 4 + 3];
		mov			edx, [esi+eax+2*JOINTMAT_SIZE+1*16+3*4]
		mov			[edi-2*JOINTQUAT_SIZE+20], edx			// q[5] = m[1 * 4 + 3];
		mov			ecx, [esi+eax+2*JOINTMAT_SIZE+2*16+3*4]
		mov			[edi-2*JOINTQUAT_SIZE+24], ecx			// q[6] = m[2 * 4 + 3];

		shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm7, xmm7, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm0, xmm0, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm2, xmm2, R_SHUFFLEPS( 1, 2, 3, 0 )

		movzx		ecx, byte ptr shuffle[3*4+0]			// ecx = k0
		movss		[edi+ecx*4-1*JOINTQUAT_SIZE], xmm7		// q[k0] = s * t;

		movzx		edx, byte ptr shuffle[3*4+1]			// edx = k1
		movss		xmm4, [esi+eax+3*JOINTMAT_SIZE+1*16+0*4]
		xorps		xmm4, xmm2
		subss		xmm4, [esi+eax+3*JOINTMAT_SIZE+0*16+1*4]
		mulss		xmm4, xmm6
		movss		[edi+edx*4-1*JOINTQUAT_SIZE], xmm4		// q[k1] = ( m[0 * 4 + 1] - s2 * m[1 * 4 + 0] ) * s;

		movzx		ecx, byte ptr shuffle[3*4+2]			// ecx = k2
		movss		xmm3, [esi+eax+3*JOINTMAT_SIZE+0*16+2*4]
		xorps		xmm3, xmm1
		subss		xmm3, [esi+eax+3*JOINTMAT_SIZE+2*16+0*4]
		mulss		xmm3, xmm6
		movss		[edi+ecx*4-1*JOINTQUAT_SIZE], xmm3		// q[k2] = ( m[2 * 4 + 0] - s1 * m[0 * 4 + 2] ) * s;

		movzx		edx, byte ptr shuffle[3*4+3]			// edx = k3
		movss		xmm4, [esi+eax+3*JOINTMAT_SIZE+2*16+1*4]
		xorps		xmm4, xmm0
		subss		xmm4, [esi+eax+3*JOINTMAT_SIZE+1*16+2*4]
		mulss		xmm4, xmm6
		movss		[edi+edx*4-1*JOINTQUAT_SIZE], xmm4		// q[k3] = ( m[1 * 4 + 2] - s0 * m[2 * 4 + 1] ) * s;

		mov			ecx, [esi+eax+3*JOINTMAT_SIZE+0*16+3*4]
		mov			[edi-1*JOINTQUAT_SIZE+16], ecx			// q[4] = m[0 * 4 + 3];
		mov			edx, [esi+eax+3*JOINTMAT_SIZE+1*16+3*4]
		mov			[edi-1*JOINTQUAT_SIZE+20], edx			// q[5] = m[1 * 4 + 3];
		mov			ecx, [esi+eax+3*JOINTMAT_SIZE+2*16+3*4]
		mov			[edi-1*JOINTQUAT_SIZE+24], ecx			// q[6] = m[2 * 4 + 3];

		add			eax, 4*JOINTMAT_SIZE
		jl			loopMat4

	done4:
		mov			eax, numJoints
		and			eax, 3
		jz			done1
		imul		eax, JOINTMAT_SIZE
		add			esi, eax
		neg			eax

	loopMat1:
		movss		xmm5, [esi+eax+0*JOINTMAT_SIZE+0*16+0*4]
		movss		xmm6, [esi+eax+0*JOINTMAT_SIZE+1*16+1*4]
		movss		xmm7, [esi+eax+0*JOINTMAT_SIZE+2*16+2*4]

		// -------------------

		movaps		xmm0, xmm5
		addss		xmm0, xmm6
		addss		xmm0, xmm7
		cmpnltss	xmm0, SIMD_SP_zero						// xmm0 = m[0 * 4 + 0] + m[1 * 4 + 1] + m[2 * 4 + 2] > 0.0f

		movaps		xmm1, xmm5
		movaps		xmm2, xmm5
		cmpnltss	xmm1, xmm6
		cmpnltss	xmm2, xmm7
		andps		xmm2, xmm1								// xmm2 = m[0 * 4 + 0] > m[1 * 4 + 1] && m[0 * 4 + 0] > m[2 * 4 + 2]

		movaps		xmm4, xmm6
		cmpnltss	xmm4, xmm7								// xmm3 = m[1 * 4 + 1] > m[2 * 4 + 2]

		movaps		xmm1, xmm0
		andnps		xmm1, xmm2
		orps		xmm2, xmm0
		movaps		xmm3, xmm2
		andnps		xmm2, xmm4
		orps		xmm3, xmm2
		xorps		xmm3, SIMD_SP_not

		andps		xmm0, SIMD_DW_mat2quatShuffle0
		movaps		xmm4, xmm1
		andps		xmm4, SIMD_DW_mat2quatShuffle1
		orps		xmm0, xmm4
		movaps		xmm4, xmm2
		andps		xmm4, SIMD_DW_mat2quatShuffle2
		orps		xmm0, xmm4
		movaps		xmm4, xmm3
		andps		xmm4, SIMD_DW_mat2quatShuffle3
		orps		xmm4, xmm0

		movss		shuffle, xmm4

		movaps		xmm0, xmm2
		orps		xmm0, xmm3								// xmm0 = xmm2 | xmm3	= s0
		orps		xmm2, xmm1								// xmm2 = xmm1 | xmm2	= s2
		orps		xmm1, xmm3								// xmm1 = xmm1 | xmm3	= s1

		andps		xmm0, SIMD_SP_signBitMask
		andps		xmm1, SIMD_SP_signBitMask
		andps		xmm2, SIMD_SP_signBitMask

		xorps		xmm5, xmm0
		xorps		xmm6, xmm1
		xorps		xmm7, xmm2
		addss		xmm5, xmm6
		addss		xmm7, SIMD_SP_one
		addss		xmm5, xmm7								// xmm5 = t

		movss		xmm7, xmm5								// xmm7 = t
		rsqrtss		xmm6, xmm5
		mulss		xmm5, xmm6
		mulss		xmm5, xmm6
		subss		xmm5, SIMD_SP_rsqrt_c0
		mulss		xmm6, SIMD_SP_mat2quat_rsqrt_c1
		mulss		xmm6, xmm5								// xmm5 = s

		mulss		xmm7, xmm6								// xmm7 = s * t
		xorps		xmm6, SIMD_SP_signBitMask				// xmm6 = -s

		// -------------------

		movzx		ecx, byte ptr shuffle[0]				// ecx = k0
		add			edi, JOINTQUAT_SIZE
		movss		[edi+ecx*4-1*JOINTQUAT_SIZE], xmm7		// q[k0] = s * t;

		movzx		edx, byte ptr shuffle[1]				// edx = k1
		movss		xmm4, [esi+eax+0*JOINTMAT_SIZE+1*16+0*4]
		xorps		xmm4, xmm2
		subss		xmm4, [esi+eax+0*JOINTMAT_SIZE+0*16+1*4]
		mulss		xmm4, xmm6
		movss		[edi+edx*4-1*JOINTQUAT_SIZE], xmm4		// q[k1] = ( m[0 * 4 + 1] - s2 * m[1 * 4 + 0] ) * s;

		movzx		ecx, byte ptr shuffle[2]				// ecx = k2
		movss		xmm3, [esi+eax+0*JOINTMAT_SIZE+0*16+2*4]
		xorps		xmm3, xmm1
		subss		xmm3, [esi+eax+0*JOINTMAT_SIZE+2*16+0*4]
		mulss		xmm3, xmm6
		movss		[edi+ecx*4-1*JOINTQUAT_SIZE], xmm3		// q[k2] = ( m[2 * 4 + 0] - s1 * m[0 * 4 + 2] ) * s;

		movzx		edx, byte ptr shuffle[3]				// edx = k3
		movss		xmm4, [esi+eax+0*JOINTMAT_SIZE+2*16+1*4]
		xorps		xmm4, xmm0
		subss		xmm4, [esi+eax+0*JOINTMAT_SIZE+1*16+2*4]
		mulss		xmm4, xmm6
		movss		[edi+edx*4-1*JOINTQUAT_SIZE], xmm4		// q[k3] = ( m[1 * 4 + 2] - s0 * m[2 * 4 + 1] ) * s;

		mov			ecx, [esi+eax+0*JOINTMAT_SIZE+0*16+3*4]
		mov			[edi-1*JOINTQUAT_SIZE+16], ecx			// q[4] = m[0 * 4 + 3];
		mov			edx, [esi+eax+0*JOINTMAT_SIZE+1*16+3*4]
		mov			[edi-1*JOINTQUAT_SIZE+20], edx			// q[5] = m[1 * 4 + 3];
		mov			ecx, [esi+eax+0*JOINTMAT_SIZE+2*16+3*4]
		mov			[edi-1*JOINTQUAT_SIZE+24], ecx			// q[6] = m[2 * 4 + 3];

		add			eax, JOINTMAT_SIZE
		jl			loopMat1

	done1:
	}

#elif 0

	for ( int i = 0; i < numJoints; i++ ) {
		float s0, s1, s2;
		int k0, k1, k2, k3;

		float *q = jointQuats[i].q.toFloatPtr();
		const float *m = jointMats[i].toFloatPtr();

		if ( m[0 * 4 + 0] + m[1 * 4 + 1] + m[2 * 4 + 2] > 0.0f ) {

			k0 = 3;
			k1 = 2;
			k2 = 1;
			k3 = 0;
			s0 = 1.0f;
			s1 = 1.0f;
			s2 = 1.0f;

		} else if ( m[0 * 4 + 0] > m[1 * 4 + 1] && m[0 * 4 + 0] > m[2 * 4 + 2] ) {

			k0 = 0;
			k1 = 1;
			k2 = 2;
			k3 = 3;
			s0 = 1.0f;
			s1 = -1.0f;
			s2 = -1.0f;

		} else if ( m[1 * 4 + 1] > m[2 * 4 + 2] ) {

			k0 = 1;
			k1 = 0;
			k2 = 3;
			k3 = 2;
			s0 = -1.0f;
			s1 = 1.0f;
			s2 = -1.0f;

		} else {

			k0 = 2;
			k1 = 3;
			k2 = 0;
			k3 = 1;
			s0 = -1.0f;
			s1 = -1.0f;
			s2 = 1.0f;

		}

		float t = s0 * m[0 * 4 + 0] + s1 * m[1 * 4 + 1] + s2 * m[2 * 4 + 2] + 1.0f;
		float s = CMath::invSqrt( t ) * 0.5f;

		q[k0] = s * t;
		q[k1] = ( m[0 * 4 + 1] - s2 * m[1 * 4 + 0] ) * s;
		q[k2] = ( m[2 * 4 + 0] - s1 * m[0 * 4 + 2] ) * s;
		q[k3] = ( m[1 * 4 + 2] - s0 * m[2 * 4 + 1] ) * s;

		q[4] = m[0 * 4 + 3];
		q[5] = m[1 * 4 + 3];
		q[6] = m[2 * 4 + 3];
	}

#elif 1

	for ( int i = 0; i < numJoints; i++ ) {

		float *q = jointQuats[i].q.toFloatPtr();
		const float *m = jointMats[i].toFloatPtr();

		if ( m[0 * 4 + 0] + m[1 * 4 + 1] + m[2 * 4 + 2] > 0.0f ) {

			float t = + m[0 * 4 + 0] + m[1 * 4 + 1] + m[2 * 4 + 2] + 1.0f;
			float s = CMath::invSqrt( t ) * 0.5f;

			q[3] = s * t;
			q[2] = ( m[0 * 4 + 1] - m[1 * 4 + 0] ) * s;
			q[1] = ( m[2 * 4 + 0] - m[0 * 4 + 2] ) * s;
			q[0] = ( m[1 * 4 + 2] - m[2 * 4 + 1] ) * s;

		} else if ( m[0 * 4 + 0] > m[1 * 4 + 1] && m[0 * 4 + 0] > m[2 * 4 + 2] ) {

			float t = + m[0 * 4 + 0] - m[1 * 4 + 1] - m[2 * 4 + 2] + 1.0f;
			float s = CMath::invSqrt( t ) * 0.5f;

			q[0] = s * t;
			q[1] = ( m[0 * 4 + 1] + m[1 * 4 + 0] ) * s;
			q[2] = ( m[2 * 4 + 0] + m[0 * 4 + 2] ) * s;
			q[3] = ( m[1 * 4 + 2] - m[2 * 4 + 1] ) * s;

		} else if ( m[1 * 4 + 1] > m[2 * 4 + 2] ) {

			float t = - m[0 * 4 + 0] + m[1 * 4 + 1] - m[2 * 4 + 2] + 1.0f;
			float s = CMath::invSqrt( t ) * 0.5f;

			q[1] = s * t;
			q[0] = ( m[0 * 4 + 1] + m[1 * 4 + 0] ) * s;
			q[3] = ( m[2 * 4 + 0] - m[0 * 4 + 2] ) * s;
			q[2] = ( m[1 * 4 + 2] + m[2 * 4 + 1] ) * s;

		} else {

			float t = - m[0 * 4 + 0] - m[1 * 4 + 1] + m[2 * 4 + 2] + 1.0f;
			float s = CMath::invSqrt( t ) * 0.5f;

			q[2] = s * t;
			q[3] = ( m[0 * 4 + 1] - m[1 * 4 + 0] ) * s;
			q[0] = ( m[2 * 4 + 0] + m[0 * 4 + 2] ) * s;
			q[1] = ( m[1 * 4 + 2] + m[2 * 4 + 1] ) * s;

		}

		q[4] = m[0 * 4 + 3];
		q[5] = m[1 * 4 + 3];
		q[6] = m[2 * 4 + 3];
	}

#endif
}

/*
============
CSIMD_SSE::transformJoints
============
*/
void VPCALL CSIMD_SSE::transformJoints( CMatJoint3x4 *jointMats, const int *parents, const int firstJoint, const int lastJoint ) {
#if 1

	SMF_ASSERT( sizeof( CMatJoint3x4 ) == JOINTMAT_SIZE );

	__asm {

		mov			ecx, firstJoint
		mov			eax, lastJoint
		sub			eax, ecx
		jl			done
		imul		ecx, 4
		mov			edi, parents
		add			edi, ecx
		imul		ecx, 12
		mov			esi, jointMats
		imul		eax, 4
		add			edi, eax
		neg			eax

	loopJoint:

		movaps		xmm0, [esi+ecx+ 0]						// xmm0 = m0, m1, m2, t0
		mov			edx, [edi+eax]
		movaps		xmm1, [esi+ecx+16]						// xmm1 = m2, m3, m4, t1
		imul		edx, JOINTMAT_SIZE
		movaps		xmm2, [esi+ecx+32]						// xmm2 = m5, m6, m7, t2

		movss		xmm4, [esi+edx+ 0]
		shufps		xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm4, xmm0

		movss		xmm5, [esi+edx+ 4]
		shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm5, xmm1
		addps		xmm4, xmm5
		movss		xmm6, [esi+edx+ 8]
		shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm6, xmm2
		addps		xmm4, xmm6

		movss		xmm5, [esi+edx+16]
		shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm5, xmm0

		movss		xmm7, [esi+edx+12]
		shufps		xmm7, xmm7, R_SHUFFLEPS( 1, 2, 3, 0 )
		addps		xmm4, xmm7

		movaps		[esi+ecx+ 0], xmm4

		movss		xmm6, [esi+edx+20]
		shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm6, xmm1
		addps		xmm5, xmm6
		movss		xmm7, [esi+edx+24]
		shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm7, xmm2
		addps		xmm5, xmm7

		movss		xmm6, [esi+edx+32]
		shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm6, xmm0

		movss		xmm3, [esi+edx+28]
		shufps		xmm3, xmm3, R_SHUFFLEPS( 1, 2, 3, 0 )
		addps		xmm5, xmm3

		movaps		[esi+ecx+16], xmm5

		movss		xmm7, [esi+edx+36]
		shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm7, xmm1
		addps		xmm6, xmm7
		movss		xmm3, [esi+edx+40]
		shufps		xmm3, xmm3, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm3, xmm2
		addps		xmm6, xmm3

		movss		xmm7, [esi+edx+44]
		shufps		xmm7, xmm7, R_SHUFFLEPS( 1, 2, 3, 0 )
		addps		xmm6, xmm7

		movaps		[esi+ecx+32], xmm6

		add			ecx, JOINTMAT_SIZE
		add			eax, 4
		jle			loopJoint
	done:
	}

#else

	int i;

	for( i = firstJoint; i <= lastJoint; i++ ) {
		SMF_ASSERT( parents[i] < i );
		jointMats[i] *= jointMats[parents[i]];
	}

#endif
}

/*
============
CSIMD_SSE::untransformJoints
============
*/
void VPCALL CSIMD_SSE::untransformJoints( CMatJoint3x4 *jointMats, const int *parents, const int firstJoint, const int lastJoint ) {
#if 1

	SMF_ASSERT( sizeof( CMatJoint3x4 ) == JOINTMAT_SIZE );

	__asm {

		mov			edx, firstJoint
		mov			eax, lastJoint
		mov			ecx, eax
		sub			eax, edx
		jl			done
		mov			esi, jointMats
		imul		ecx, JOINTMAT_SIZE
		imul		edx, 4
		mov			edi, parents
		add			edi, edx
		imul		eax, 4

	loopJoint:

		movaps		xmm0, [esi+ecx+ 0]						// xmm0 = m0, m1, m2, t0
		mov			edx, [edi+eax]
		movaps		xmm1, [esi+ecx+16]						// xmm1 = m2, m3, m4, t1
		imul		edx, JOINTMAT_SIZE
		movaps		xmm2, [esi+ecx+32]						// xmm2 = m5, m6, m7, t2

		movss		xmm6, [esi+edx+12]
		shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 2, 3, 0 )
		subps		xmm0, xmm6
		movss		xmm7, [esi+edx+28]
		shufps		xmm7, xmm7, R_SHUFFLEPS( 1, 2, 3, 0 )
		subps		xmm1, xmm7
		movss		xmm3, [esi+edx+44]
		shufps		xmm3, xmm3, R_SHUFFLEPS( 1, 2, 3, 0 )
		subps		xmm2, xmm3

		movss		xmm4, [esi+edx+ 0]
		shufps		xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm4, xmm0
		movss		xmm5, [esi+edx+16]
		shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm5, xmm1
		addps		xmm4, xmm5
		movss		xmm6, [esi+edx+32]
		shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm6, xmm2
		addps		xmm4, xmm6

		movaps		[esi+ecx+ 0], xmm4

		movss		xmm5, [esi+edx+ 4]
		shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm5, xmm0
		movss		xmm6, [esi+edx+20]
		shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm6, xmm1
		addps		xmm5, xmm6
		movss		xmm7, [esi+edx+36]
		shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm7, xmm2
		addps		xmm5, xmm7

		movaps		[esi+ecx+16], xmm5

		movss		xmm6, [esi+edx+ 8]
		shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm6, xmm0
		movss		xmm7, [esi+edx+24]
		shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm7, xmm1
		addps		xmm6, xmm7
		movss		xmm3, [esi+edx+40]
		shufps		xmm3, xmm3, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm3, xmm2
		addps		xmm6, xmm3

		movaps		[esi+ecx+32], xmm6

		sub			ecx, JOINTMAT_SIZE
		sub			eax, 4
		jge			loopJoint
	done:
	}

#else

	int i;

	for( i = lastJoint; i >= firstJoint; i-- ) {
		SMF_ASSERT( parents[i] < i );
		jointMats[i] /= jointMats[parents[i]];
	}

#endif
}

/*
============
CSIMD_SSE::transformVerts
============
*/
void VPCALL CSIMD_SSE::transformVerts( CVertex *verts, const int numVerts, const CMatJoint3x4 *joints, const CVec4D *weights, const int *index, const int numWeights ) {
#if 1

	SMF_ASSERT( sizeof( CVertex ) == DRAWVERT_SIZE );
	SMF_ASSERT( (int)&((CVertex *)0)->xyz == DRAWVERT_XYZ_OFFSET );
	SMF_ASSERT( sizeof( CVec4D ) == JOINTWEIGHT_SIZE );
	SMF_ASSERT( sizeof( CMatJoint3x4 ) == JOINTMAT_SIZE );

	__asm
	{
		mov			eax, numVerts
		test		eax, eax
		jz			done
		imul		eax, DRAWVERT_SIZE

		mov			ecx, verts
		mov			edx, index
		mov			esi, weights
		mov			edi, joints

		add			ecx, eax
		neg			eax

	loopVert:
		mov			ebx, [edx]
		movaps		xmm2, [esi]
		add			edx, 8
		movaps		xmm0, xmm2
		add			esi, JOINTWEIGHT_SIZE
		movaps		xmm1, xmm2

		mulps		xmm0, [edi+ebx+ 0]						// xmm0 = m0, m1, m2, t0
		mulps		xmm1, [edi+ebx+16]						// xmm1 = m3, m4, m5, t1
		mulps		xmm2, [edi+ebx+32]						// xmm2 = m6, m7, m8, t2

		cmp			dword ptr [edx-4], 0

		jne			doneWeight

	loopWeight:
		mov			ebx, [edx]
		movaps		xmm5, [esi]
		add			edx, 8
		movaps		xmm3, xmm5
		add			esi, JOINTWEIGHT_SIZE
		movaps		xmm4, xmm5

		mulps		xmm3, [edi+ebx+ 0]						// xmm3 = m0, m1, m2, t0
		mulps		xmm4, [edi+ebx+16]						// xmm4 = m3, m4, m5, t1
		mulps		xmm5, [edi+ebx+32]						// xmm5 = m6, m7, m8, t2

		cmp			dword ptr [edx-4], 0

		addps		xmm0, xmm3
		addps		xmm1, xmm4
		addps		xmm2, xmm5

		je			loopWeight

	doneWeight:
		add			eax, DRAWVERT_SIZE

		movaps		xmm6, xmm0								// xmm6 =    m0,    m1,          m2,          t0
		unpcklps	xmm6, xmm1								// xmm6 =    m0,    m3,          m1,          m4
		unpckhps	xmm0, xmm1								// xmm1 =    m2,    m5,          t0,          t1
		addps		xmm6, xmm0								// xmm6 = m0+m2, m3+m5,       m1+t0,       m4+t1

		movaps		xmm7, xmm2								// xmm7 =    m6,    m7,          m8,          t2
		movlhps		xmm2, xmm6								// xmm2 =    m6,    m7,       m0+m2,       m3+m5
		movhlps		xmm6, xmm7								// xmm6 =    m8,    t2,       m1+t0,       m4+t1
		addps		xmm6, xmm2								// xmm6 = m6+m8, m7+t2, m0+m1+m2+t0, m3+m4+m5+t1

		movhps		[ecx+eax-DRAWVERT_SIZE+0], xmm6

		movaps		xmm5, xmm6								// xmm5 = m6+m8, m7+t2
		shufps		xmm5, xmm5, R_SHUFFLEPS( 1, 0, 2, 3 )	// xmm5 = m7+t2, m6+m8
		addss		xmm5, xmm6								// xmm5 = m6+m8+m7+t2

		movss		[ecx+eax-DRAWVERT_SIZE+8], xmm5

		jl			loopVert
	done:
	}

#else

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

#endif
}

/*
============
CSIMD_SSE::tracePointCull
============
*/
void VPCALL CSIMD_SSE::tracePointCull( sf_u8 *cullBits, sf_u8 &totalOr, const float radius, const CPlane *planes, const CVertex *verts, const int numVerts ) {
#if 1

	SMF_ASSERT( sizeof( CVertex ) == DRAWVERT_SIZE );
	SMF_ASSERT( (int)&((CVertex *)0)->xyz == DRAWVERT_XYZ_OFFSET );

	__asm {
		push		ebx
		mov			eax, numVerts
		test		eax, eax
		jz			done

		mov			edi, planes
		movlps		xmm1, [edi]								// xmm1 =  0,  1,  X,  X
		movhps		xmm1, [edi+16]							// xmm1 =  0,  1,  4,  5
		movlps		xmm3, [edi+8]							// xmm3 =  2,  3,  X,  X
		movhps		xmm3, [edi+24]							// xmm3 =  2,  3,  6,  7
		movlps		xmm4, [edi+32]							// xmm4 =  8,  9,  X,  X
		movhps		xmm4, [edi+48]							// xmm4 =  8,  9, 12, 13
		movlps		xmm5, [edi+40]							// xmm5 = 10, 11,  X,  X
		movhps		xmm5, [edi+56]							// xmm5 = 10, 11, 14, 15
		movaps		xmm0, xmm1								// xmm0 =  0,  1,  4,  5
		shufps		xmm0, xmm4, R_SHUFFLEPS( 0, 2, 0, 2 )	// xmm0 =  0,  4,  8, 12
		shufps		xmm1, xmm4, R_SHUFFLEPS( 1, 3, 1, 3 )	// xmm1 =  1,  5,  9, 13
		movaps		xmm2, xmm3								// xmm2 =  2,  3,  6,  7
		shufps		xmm2, xmm5, R_SHUFFLEPS( 0, 2, 0, 2 )	// xmm2 =  2,  6, 10, 14
		shufps		xmm3, xmm5, R_SHUFFLEPS( 1, 3, 1, 3 )	// xmm3 =  3,  7, 11, 15
		movss		xmm7, radius
		shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )

		xor			edx, edx
		mov			esi, verts
		mov			edi, cullBits
		imul		eax, DRAWVERT_SIZE
		add			esi, eax
		neg			eax

	loopVert:
		movss		xmm4, [esi+eax+DRAWVERT_XYZ_OFFSET+0]
		shufps		xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )
		movss		xmm5, [esi+eax+DRAWVERT_XYZ_OFFSET+4]
		mulps		xmm4, xmm0
		shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )
		movss		xmm6, [esi+eax+DRAWVERT_XYZ_OFFSET+8]
		mulps		xmm5, xmm1
		shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
		addps		xmm4, xmm5
		mulps		xmm6, xmm2
		addps		xmm4, xmm3
		addps		xmm4, xmm6
		movaps		xmm5, xmm4
		xorps		xmm5, SIMD_SP_signBitMask
		cmpltps		xmm4, xmm7
		movmskps	ecx, xmm4
		cmpltps		xmm5, xmm7
		movmskps	ebx, xmm5
		shl			cx, 4
		or			cl, bl
		inc			edi
		or			dl, cl
		add			eax, DRAWVERT_SIZE
		mov			byte ptr [edi-1], cl
		jl			loopVert

	done:
		mov			esi, totalOr
        mov			byte ptr [esi], dl
		pop			ebx
	}

#else

	int i;
	sf_u8 tOr;

	tOr = 0;

	for ( i = 0; i < numVerts; i++ ) {
		sf_u8 bits;
		float d0, d1, d2, d3, t;
		const CVec3D &v = verts[i].xyz;

		d0 = planes[0][0] * v[0] + planes[0][1] * v[1] + planes[0][2] * v[2] + planes[0][3];
		d1 = planes[1][0] * v[0] + planes[1][1] * v[1] + planes[1][2] * v[2] + planes[1][3];
		d2 = planes[2][0] * v[0] + planes[2][1] * v[1] + planes[2][2] * v[2] + planes[2][3];
		d3 = planes[3][0] * v[0] + planes[3][1] * v[1] + planes[3][2] * v[2] + planes[3][3];

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

#endif
}

/*
============
CSIMD_SSE::decalPointCull
============
*/
void VPCALL CSIMD_SSE::decalPointCull( sf_u8 *cullBits, const CPlane *planes, const CVertex *verts, const int numVerts ) {
#if 1

	ALIGN16( float p0[4] );
	ALIGN16( float p1[4] );
	ALIGN16( float p2[4] );
	ALIGN16( float p3[4] );
	ALIGN16( float p4[4] );
	ALIGN16( float p5[4] );
	ALIGN16( float p6[4] );
	ALIGN16( float p7[4] );

	SMF_ASSERT( sizeof( CVertex ) == DRAWVERT_SIZE );
	SMF_ASSERT( (int)&((CVertex *)0)->xyz == DRAWVERT_XYZ_OFFSET );

	__asm {
		mov			ecx, planes
		movlps		xmm1, [ecx]								// xmm1 =  0,  1,  X,  X
		movhps		xmm1, [ecx+16]							// xmm1 =  0,  1,  4,  5
		movlps		xmm3, [ecx+8]							// xmm3 =  2,  3,  X,  X
		movhps		xmm3, [ecx+24]							// xmm3 =  2,  3,  6,  7
		movlps		xmm4, [ecx+32]							// xmm4 =  8,  9,  X,  X
		movhps		xmm4, [ecx+48]							// xmm4 =  8,  9, 12, 13
		movlps		xmm5, [ecx+40]							// xmm5 = 10, 11,  X,  X
		movhps		xmm5, [ecx+56]							// xmm5 = 10, 11, 14, 15
		movaps		xmm0, xmm1								// xmm0 =  0,  1,  4,  5
		shufps		xmm0, xmm4, R_SHUFFLEPS( 0, 2, 0, 2 )	// xmm0 =  0,  4,  8, 12
		shufps		xmm1, xmm4, R_SHUFFLEPS( 1, 3, 1, 3 )	// xmm1 =  1,  5,  9, 13
		movaps		xmm2, xmm3								// xmm2 =  2,  3,  6,  7
		shufps		xmm2, xmm5, R_SHUFFLEPS( 0, 2, 0, 2 )	// xmm2 =  2,  6, 10, 14
		shufps		xmm3, xmm5, R_SHUFFLEPS( 1, 3, 1, 3 )	// xmm3 =  3,  7, 11, 15

		movaps		p0, xmm0
		movaps		p1, xmm1
		movaps		p2, xmm2
		movaps		p3, xmm3

		movlps		xmm4, [ecx+64]							// xmm4 = p40, p41,   X,   X
		movhps		xmm4, [ecx+80]							// xmm4 = p40, p41, p50, p51
		movaps		xmm5, xmm4								// xmm5 = p40, p41, p50, p51
		shufps		xmm4, xmm4, R_SHUFFLEPS( 0, 2, 0, 2 )	// xmm4 = p40, p50, p40, p50
		shufps		xmm5, xmm5, R_SHUFFLEPS( 1, 3, 1, 3 )	// xmm5 = p41, p51, p41, p51
		movlps		xmm6, [ecx+72]							// xmm6 = p42, p43,   X,   X
		movhps		xmm6, [ecx+88]							// xmm6 = p42, p43, p52, p53
		movaps		xmm7, xmm6								// xmm7 = p42, p43, p52, p53
		shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 2, 0, 2 )	// xmm6 = p42, p52, p42, p52
		shufps		xmm7, xmm7, R_SHUFFLEPS( 1, 3, 1, 3 )	// xmm7 = p43, p53, p43, p53

		movaps		p4, xmm4
		movaps		p5, xmm5
		movaps		p6, xmm6
		movaps		p7, xmm7

		mov			esi, verts
		mov			edi, cullBits
		mov			eax, numVerts
		and			eax, ~1
		jz			done2
		imul		eax, DRAWVERT_SIZE
		add			esi, eax
		neg			eax

	loopVert2:
		movaps		xmm6, p0
		movss		xmm0, [esi+eax+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]
		shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm6, xmm0
		movaps		xmm7, p1
		movss		xmm1, [esi+eax+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+4]
		shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm7, xmm1
		addps		xmm6, xmm7
		movaps		xmm7, p2
		movss		xmm2, [esi+eax+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+8]
		shufps		xmm2, xmm2, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm7, xmm2
		addps		xmm6, xmm7
		addps		xmm6, p3

		cmpnltps	xmm6, SIMD_SP_zero
		movmskps	ecx, xmm6

		movaps		xmm6, p0
		movss		xmm3, [esi+eax+1*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]
		shufps		xmm3, xmm3, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm6, xmm3
		movaps		xmm7, p1
		movss		xmm4, [esi+eax+1*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+4]
		shufps		xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm7, xmm4
		addps		xmm6, xmm7
		movaps		xmm7, p2
		movss		xmm5, [esi+eax+1*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+8]
		shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm7, xmm5
		addps		xmm6, xmm7
		addps		xmm6, p3

		cmpnltps	xmm6, SIMD_SP_zero
		movmskps	edx, xmm6
		mov			ch, dl

		shufps		xmm0, xmm3, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm0, p4
		shufps		xmm1, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm1, p5
		addps		xmm0, xmm1
		shufps		xmm2, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm2, p6
		addps		xmm0, xmm2
		addps		xmm0, p7

		cmpnltps	xmm0, SIMD_SP_zero
		movmskps	edx, xmm0

		add			edi, 2

		mov			dh, dl
		shl			dl, 4
		shl			dh, 2
		and			edx, (3<<4)|(3<<12)
		or			ecx, edx

		add			eax, 2*DRAWVERT_SIZE
		mov			word ptr [edi-2], cx
		jl			loopVert2

	done2:

		mov			eax, numVerts
		and			eax, 1
		jz			done

		movaps		xmm6, p0
		movss		xmm0, [esi+DRAWVERT_XYZ_OFFSET+0]
		shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm6, xmm0
		movaps		xmm7, p1
		movss		xmm1, [esi+DRAWVERT_XYZ_OFFSET+4]
		shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm7, xmm1
		addps		xmm6, xmm7
		movaps		xmm7, p2
		movss		xmm2, [esi+DRAWVERT_XYZ_OFFSET+8]
		shufps		xmm2, xmm2, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm7, xmm2
		addps		xmm6, xmm7
		addps		xmm6, p3

		cmpnltps	xmm6, SIMD_SP_zero
		movmskps	ecx, xmm6

		mulps		xmm0, p4
		mulps		xmm1, p5
		addps		xmm0, xmm1
		mulps		xmm2, p6
		addps		xmm0, xmm2
		addps		xmm0, p7

		cmpnltps	xmm0, SIMD_SP_zero
		movmskps	edx, xmm0

		and			edx, 3
		shl			edx, 4
		or			ecx, edx

		mov			byte ptr [edi], cl

	done:
	}


#else

	int i;

	for ( i = 0; i < numVerts; i += 2 ) {
		unsigned short bits0, bits1;
		float d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11;
		const CVec3D &v0 = verts[i+0].xyz;
		const CVec3D &v1 = verts[i+1].xyz;

		d0  = planes[0][0] * v0[0] + planes[0][1] * v0[1] + planes[0][2] * v0[2] + planes[0][3];
		d1  = planes[1][0] * v0[0] + planes[1][1] * v0[1] + planes[1][2] * v0[2] + planes[1][3];
		d2  = planes[2][0] * v0[0] + planes[2][1] * v0[1] + planes[2][2] * v0[2] + planes[2][3];
		d3  = planes[3][0] * v0[0] + planes[3][1] * v0[1] + planes[3][2] * v0[2] + planes[3][3];

		d4  = planes[4][0] * v0[0] + planes[4][1] * v0[1] + planes[4][2] * v0[2] + planes[4][3];
		d5  = planes[5][0] * v0[0] + planes[5][1] * v0[1] + planes[5][2] * v0[2] + planes[5][3];
		d10 = planes[4][0] * v1[0] + planes[4][1] * v1[1] + planes[4][2] * v1[2] + planes[4][3];
		d11 = planes[5][0] * v1[0] + planes[5][1] * v1[1] + planes[5][2] * v1[2] + planes[5][3];

		d6  = planes[0][0] * v1[0] + planes[0][1] * v1[1] + planes[0][2] * v1[2] + planes[0][3];
		d7  = planes[1][0] * v1[0] + planes[1][1] * v1[1] + planes[1][2] * v1[2] + planes[1][3];
		d8  = planes[2][0] * v1[0] + planes[2][1] * v1[1] + planes[2][2] * v1[2] + planes[2][3];
		d9  = planes[3][0] * v1[0] + planes[3][1] * v1[1] + planes[3][2] * v1[2] + planes[3][3];

		bits0  = FLOATSIGNBITSET( d0  ) << (0+0);
		bits0 |= FLOATSIGNBITSET( d1  ) << (0+1);
		bits0 |= FLOATSIGNBITSET( d2  ) << (0+2);
		bits0 |= FLOATSIGNBITSET( d3  ) << (0+3);
		bits0 |= FLOATSIGNBITSET( d4  ) << (0+4);
		bits0 |= FLOATSIGNBITSET( d5  ) << (0+5);

		bits1  = FLOATSIGNBITSET( d6  ) << (8+0);
		bits1 |= FLOATSIGNBITSET( d7  ) << (8+1);
		bits1 |= FLOATSIGNBITSET( d8  ) << (8+2);
		bits1 |= FLOATSIGNBITSET( d9  ) << (8+3);
		bits1 |= FLOATSIGNBITSET( d10 ) << (8+4);
		bits1 |= FLOATSIGNBITSET( d11 ) << (8+5);

		*(unsigned short *)(cullBits + i) = ( bits0 | bits1 ) ^ 0x3F3F;
	}

	if ( numVerts & 1 ) {
		sf_u8 bits;
		float d0, d1, d2, d3, d4, d5;
		const CVec3D &v = verts[numVerts - 1].xyz;

		d0 = planes[0][0] * v[0] + planes[0][1] * v[1] + planes[0][2] * v[2] + planes[0][3];
		d1 = planes[1][0] * v[0] + planes[1][1] * v[1] + planes[1][2] * v[2] + planes[1][3];
		d2 = planes[2][0] * v[0] + planes[2][1] * v[1] + planes[2][2] * v[2] + planes[2][3];
		d3 = planes[3][0] * v[0] + planes[3][1] * v[1] + planes[3][2] * v[2] + planes[3][3];

		d4 = planes[4][0] * v[0] + planes[4][1] * v[1] + planes[4][2] * v[2] + planes[4][3];
		d5 = planes[5][0] * v[0] + planes[5][1] * v[1] + planes[5][2] * v[2] + planes[5][3];

		bits  = FLOATSIGNBITSET( d0 ) << 0;
		bits |= FLOATSIGNBITSET( d1 ) << 1;
		bits |= FLOATSIGNBITSET( d2 ) << 2;
		bits |= FLOATSIGNBITSET( d3 ) << 3;

		bits |= FLOATSIGNBITSET( d4 ) << 4;
		bits |= FLOATSIGNBITSET( d5 ) << 5;

		cullBits[numVerts - 1] = bits ^ 0x3F;		// flip lower 6 bits
	}

#endif
}

/*
============
CSIMD_SSE::overlayPointCull
============
*/
void VPCALL CSIMD_SSE::overlayPointCull( sf_u8 *cullBits, CVec2D *texCoords, const CPlane *planes, const CVertex *verts, const int numVerts ) {
#if 1

	SMF_ASSERT( sizeof( CVertex ) == DRAWVERT_SIZE );
	SMF_ASSERT( (int)&((CVertex *)0)->xyz == DRAWVERT_XYZ_OFFSET );

	__asm {
		mov			eax, numVerts
		mov			edx, verts
		mov			esi, texCoords
		mov			edi, cullBits

		mov			ecx, planes
		movss		xmm4, [ecx+ 0]
		movss		xmm5, [ecx+16]
		shufps		xmm4, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )
		shufps		xmm4, xmm4, R_SHUFFLEPS( 0, 2, 0, 2 )
		movss		xmm5, [ecx+ 4]
		movss		xmm6, [ecx+20]
		shufps		xmm5, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
		shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 2, 0, 2 )
		movss		xmm6, [ecx+ 8]
		movss		xmm7, [ecx+24]
		shufps		xmm6, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
		shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 2, 0, 2 )
		movss		xmm7, [ecx+12]
		movss		xmm0, [ecx+28]
		shufps		xmm7, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
		shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 2, 0, 2 )

		and			eax, ~1
		jz			done2
		add			edi, eax
		neg			eax

	loopVert2:
		movss		xmm0, [edx+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]
		movss		xmm1, [edx+1*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]
		shufps		xmm0, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm0, xmm4
		movss		xmm1, [edx+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+4]
		movss		xmm2, [edx+1*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+4]
		shufps		xmm1, xmm2, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm1, xmm5
		movss		xmm2, [edx+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+8]
		movss		xmm3, [edx+1*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+8]
		shufps		xmm2, xmm3, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm2, xmm6
		addps		xmm0, xmm1
		addps		xmm0, xmm2
		addps		xmm0, xmm7
		movaps		[esi], xmm0
		movaps		xmm1, xmm0
		movaps		xmm2, SIMD_SP_one
		subps		xmm2, xmm0
		shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 1, 0, 1 )
		shufps		xmm1, xmm2, R_SHUFFLEPS( 2, 3, 2, 3 )
		add			edx, 2*DRAWVERT_SIZE
		movmskps	ecx, xmm0
		mov			byte ptr [edi+eax+0], cl
		add			esi, 4*4
		movmskps	ecx, xmm1
		mov			byte ptr [edi+eax+1], cl
		add			eax, 2
		jl			loopVert2

	done2:
		mov			eax, numVerts
		and			eax, 1
		jz			done

		movss		xmm0, [edx+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]
		shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm0, xmm4
		movss		xmm1, [edx+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+4]
		shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm1, xmm5
		movss		xmm2, [edx+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+8]
		shufps		xmm2, xmm2, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm2, xmm6
		addps		xmm0, xmm1
		addps		xmm0, xmm2
		addps		xmm0, xmm7
		movlps		[esi], xmm0
		movaps		xmm1, xmm0
		movaps		xmm2, SIMD_SP_one
		subps		xmm2, xmm0
		shufps		xmm0, xmm2, R_SHUFFLEPS( 0, 1, 0, 1 )
		movmskps	ecx, xmm0
		mov			byte ptr [edi], cl

	done:
	}

#else

	const CPlane &p0 = planes[0];
	const CPlane &p1 = planes[1];

	for ( int i = 0; i < numVerts - 1; i += 2 ) {
		unsigned short bits;
		float d0, d1, d2, d3;

		const CVec3D &v0 = verts[i+0].xyz;
		const CVec3D &v1 = verts[i+1].xyz;

		d0 = p0[0] * v0[0] + p0[1] * v0[1] + p0[2] * v0[2] + p0[3];
		d1 = p1[0] * v0[0] + p1[1] * v0[1] + p1[2] * v0[2] + p1[3];
		d2 = p0[0] * v1[0] + p0[1] * v1[1] + p0[2] * v1[2] + p0[3];
		d3 = p1[0] * v1[0] + p1[1] * v1[1] + p1[2] * v1[2] + p1[3];

		texCoords[i+0][0] = d0;
		texCoords[i+0][1] = d1;
		texCoords[i+1][0] = d2;
		texCoords[i+1][1] = d3;

		bits  = FLOATSIGNBITSET( d0 ) << 0;
		bits |= FLOATSIGNBITSET( d1 ) << 1;
		bits |= FLOATSIGNBITSET( d2 ) << 8;
		bits |= FLOATSIGNBITSET( d3 ) << 9;

		d0 = 1.0f - d0;
		d1 = 1.0f - d1;
		d2 = 1.0f - d2;
		d3 = 1.0f - d3;

		bits |= FLOATSIGNBITSET( d0 ) << 2;
		bits |= FLOATSIGNBITSET( d1 ) << 3;
		bits |= FLOATSIGNBITSET( d2 ) << 10;
		bits |= FLOATSIGNBITSET( d3 ) << 11;

		*(unsigned short *)(cullBits + i) = bits;
	}

	if ( numVerts & 1 ) {
		sf_u8 bits;
		float d0, d1;

		const CPlane &p0 = planes[0];
		const CPlane &p1 = planes[1];
		const CVec3D &v0 = verts[numVerts - 1].xyz;

		d0 = p0[0] * v0[0] + p0[1] * v0[1] + p0[2] * v0[2] + p0[3];
		d1 = p1[0] * v0[0] + p1[1] * v0[1] + p1[2] * v0[2] + p1[3];

		texCoords[i][0] = d0;
		texCoords[i][1] = d1;

		bits  = FLOATSIGNBITSET( d0 ) << 0;
		bits |= FLOATSIGNBITSET( d1 ) << 1;

		d0 = 1.0f - d0;
		d1 = 1.0f - d1;

		bits |= FLOATSIGNBITSET( d0 ) << 2;
		bits |= FLOATSIGNBITSET( d1 ) << 3;

		cullBits[numVerts - 1] = bits;
	}

#endif
}

/*
============
CSIMD_SSE::deriveTriPlanes
============
*/
void VPCALL CSIMD_SSE::deriveTriPlanes( CPlane *planes, const CVertex *verts, const int numVerts, const int *indexes, const int numIndexes ) {
#if 1

	SMF_ASSERT( sizeof( CVertex ) == DRAWVERT_SIZE );
	SMF_ASSERT( (int)&((CVertex *)0)->xyz == DRAWVERT_XYZ_OFFSET );

	__asm {
		mov			eax, numIndexes
		shl			eax, 2
		mov			esi, verts
		mov			edi, indexes
		mov			edx, planes

		add			edi, eax
		neg			eax

		add			eax, 4*12
		jge			done4

	loopPlane4:
		mov			ebx, [edi+eax-4*12+4]
		imul		ebx, DRAWVERT_SIZE
		mov			ecx, [edi+eax-4*12+0]
		imul		ecx, DRAWVERT_SIZE

		movss		xmm0, [esi+ebx+DRAWVERT_XYZ_OFFSET+0]
		subss		xmm0, [esi+ecx+DRAWVERT_XYZ_OFFSET+0]

		movss		xmm1, [esi+ebx+DRAWVERT_XYZ_OFFSET+4]
		subss		xmm1, [esi+ecx+DRAWVERT_XYZ_OFFSET+4]

		movss		xmm2, [esi+ebx+DRAWVERT_XYZ_OFFSET+8]
		subss		xmm2, [esi+ecx+DRAWVERT_XYZ_OFFSET+8]

		mov			ebx, [edi+eax-4*12+8]
		imul		ebx, DRAWVERT_SIZE

		shufps		xmm0, xmm0, R_SHUFFLEPS( 3, 0, 1, 2 )
		shufps		xmm1, xmm1, R_SHUFFLEPS( 3, 0, 1, 2 )
		shufps		xmm2, xmm2, R_SHUFFLEPS( 3, 0, 1, 2 )

		movss		xmm3, [esi+ebx+DRAWVERT_XYZ_OFFSET+0]
		subss		xmm3, [esi+ecx+DRAWVERT_XYZ_OFFSET+0]

		movss		xmm4, [esi+ebx+DRAWVERT_XYZ_OFFSET+4]
		subss		xmm4, [esi+ecx+DRAWVERT_XYZ_OFFSET+4]

		movss		xmm5, [esi+ebx+DRAWVERT_XYZ_OFFSET+8]
		subss		xmm5, [esi+ecx+DRAWVERT_XYZ_OFFSET+8]

		mov			ebx, [edi+eax-3*12+4]
		imul		ebx, DRAWVERT_SIZE
		mov			ecx, [edi+eax-3*12+0]
		imul		ecx, DRAWVERT_SIZE

		shufps		xmm3, xmm3, R_SHUFFLEPS( 3, 0, 1, 2 )
		shufps		xmm4, xmm4, R_SHUFFLEPS( 3, 0, 1, 2 )
		shufps		xmm5, xmm5, R_SHUFFLEPS( 3, 0, 1, 2 )

		movss		xmm6, [esi+ebx+DRAWVERT_XYZ_OFFSET+0]
		subss		xmm6, [esi+ecx+DRAWVERT_XYZ_OFFSET+0]
		movss		xmm0, xmm6

		movss		xmm7, [esi+ebx+DRAWVERT_XYZ_OFFSET+4]
		subss		xmm7, [esi+ecx+DRAWVERT_XYZ_OFFSET+4]
		movss		xmm1, xmm7

		movss		xmm6, [esi+ebx+DRAWVERT_XYZ_OFFSET+8]
		subss		xmm6, [esi+ecx+DRAWVERT_XYZ_OFFSET+8]
		movss		xmm2, xmm6

		mov			ebx, [edi+eax-3*12+8]
		imul		ebx, DRAWVERT_SIZE

		shufps		xmm0, xmm0, R_SHUFFLEPS( 3, 0, 1, 2 )
		shufps		xmm1, xmm1, R_SHUFFLEPS( 3, 0, 1, 2 )
		shufps		xmm2, xmm2, R_SHUFFLEPS( 3, 0, 1, 2 )

		movss		xmm7, [esi+ebx+DRAWVERT_XYZ_OFFSET+0]
		subss		xmm7, [esi+ecx+DRAWVERT_XYZ_OFFSET+0]
		movss		xmm3, xmm7

		movss		xmm6, [esi+ebx+DRAWVERT_XYZ_OFFSET+4]
		subss		xmm6, [esi+ecx+DRAWVERT_XYZ_OFFSET+4]
		movss		xmm4, xmm6

		movss		xmm7, [esi+ebx+DRAWVERT_XYZ_OFFSET+8]
		subss		xmm7, [esi+ecx+DRAWVERT_XYZ_OFFSET+8]
		movss		xmm5, xmm7

		mov			ebx, [edi+eax-2*12+4]
		imul		ebx, DRAWVERT_SIZE
		mov			ecx, [edi+eax-2*12+0]
		imul		ecx, DRAWVERT_SIZE

		shufps		xmm3, xmm3, R_SHUFFLEPS( 3, 0, 1, 2 )
		shufps		xmm4, xmm4, R_SHUFFLEPS( 3, 0, 1, 2 )
		shufps		xmm5, xmm5, R_SHUFFLEPS( 3, 0, 1, 2 )

		movss		xmm6, [esi+ebx+DRAWVERT_XYZ_OFFSET+0]
		subss		xmm6, [esi+ecx+DRAWVERT_XYZ_OFFSET+0]
		movss		xmm0, xmm6

		movss		xmm7, [esi+ebx+DRAWVERT_XYZ_OFFSET+4]
		subss		xmm7, [esi+ecx+DRAWVERT_XYZ_OFFSET+4]
		movss		xmm1, xmm7

		movss		xmm6, [esi+ebx+DRAWVERT_XYZ_OFFSET+8]
		subss		xmm6, [esi+ecx+DRAWVERT_XYZ_OFFSET+8]
		movss		xmm2, xmm6

		mov			ebx, [edi+eax-2*12+8]
		imul		ebx, DRAWVERT_SIZE

		shufps		xmm0, xmm0, R_SHUFFLEPS( 3, 0, 1, 2 )
		shufps		xmm1, xmm1, R_SHUFFLEPS( 3, 0, 1, 2 )
		shufps		xmm2, xmm2, R_SHUFFLEPS( 3, 0, 1, 2 )

		movss		xmm7, [esi+ebx+DRAWVERT_XYZ_OFFSET+0]
		subss		xmm7, [esi+ecx+DRAWVERT_XYZ_OFFSET+0]
		movss		xmm3, xmm7

		movss		xmm6, [esi+ebx+DRAWVERT_XYZ_OFFSET+4]
		subss		xmm6, [esi+ecx+DRAWVERT_XYZ_OFFSET+4]
		movss		xmm4, xmm6

		movss		xmm7, [esi+ebx+DRAWVERT_XYZ_OFFSET+8]
		subss		xmm7, [esi+ecx+DRAWVERT_XYZ_OFFSET+8]
		movss		xmm5, xmm7

		mov			ebx, [edi+eax-1*12+4]
		imul		ebx, DRAWVERT_SIZE
		mov			ecx, [edi+eax-1*12+0]
		imul		ecx, DRAWVERT_SIZE

		shufps		xmm3, xmm3, R_SHUFFLEPS( 3, 0, 1, 2 )
		shufps		xmm4, xmm4, R_SHUFFLEPS( 3, 0, 1, 2 )
		shufps		xmm5, xmm5, R_SHUFFLEPS( 3, 0, 1, 2 )

		movss		xmm6, [esi+ebx+DRAWVERT_XYZ_OFFSET+0]
		subss		xmm6, [esi+ecx+DRAWVERT_XYZ_OFFSET+0]
		movss		xmm0, xmm6

		movss		xmm7, [esi+ebx+DRAWVERT_XYZ_OFFSET+4]
		subss		xmm7, [esi+ecx+DRAWVERT_XYZ_OFFSET+4]
		movss		xmm1, xmm7

		movss		xmm6, [esi+ebx+DRAWVERT_XYZ_OFFSET+8]
		subss		xmm6, [esi+ecx+DRAWVERT_XYZ_OFFSET+8]
		movss		xmm2, xmm6

		mov			ebx, [edi+eax-1*12+8]
		imul		ebx, DRAWVERT_SIZE

		movss		xmm7, [esi+ebx+DRAWVERT_XYZ_OFFSET+0]
		subss		xmm7, [esi+ecx+DRAWVERT_XYZ_OFFSET+0]
		movss		xmm3, xmm7

		movss		xmm6, [esi+ebx+DRAWVERT_XYZ_OFFSET+4]
		subss		xmm6, [esi+ecx+DRAWVERT_XYZ_OFFSET+4]
		movss		xmm4, xmm6

		movss		xmm7, [esi+ebx+DRAWVERT_XYZ_OFFSET+8]
		subss		xmm7, [esi+ecx+DRAWVERT_XYZ_OFFSET+8]
		movss		xmm5, xmm7

		movaps		xmm6, xmm4
		mulps		xmm6, xmm2
		movaps		xmm7, xmm5
		mulps		xmm7, xmm1
		subps		xmm6, xmm7

		mulps		xmm5, xmm0
		mulps		xmm2, xmm3
		subps		xmm5, xmm2

		mulps		xmm3, xmm1
		mulps		xmm4, xmm0
		subps		xmm3, xmm4

		movaps		xmm0, xmm6
		mulps		xmm6, xmm6
		movaps		xmm1, xmm5
		mulps		xmm5, xmm5
		movaps		xmm2, xmm3
		mulps		xmm3, xmm3

		addps		xmm3, xmm5
		addps		xmm3, xmm6
		rsqrtps		xmm3, xmm3

		add			edx, 4*16
		mov			ecx, [edi+eax-1*12+0]
		imul		ecx, DRAWVERT_SIZE

		mulps		xmm0, xmm3
		mulps		xmm1, xmm3
		mulps		xmm2, xmm3

		movss		[edx-1*16+0], xmm0
		movss		[edx-1*16+4], xmm1
		movss		[edx-1*16+8], xmm2

		mulss		xmm0, [esi+ecx+DRAWVERT_XYZ_OFFSET+0]
		mulss		xmm1, [esi+ecx+DRAWVERT_XYZ_OFFSET+4]
		mulss		xmm2, [esi+ecx+DRAWVERT_XYZ_OFFSET+8]

		xorps		xmm0, SIMD_SP_singleSignBitMask
		subss		xmm0, xmm1
		subss		xmm0, xmm2
		movss		[edx-1*16+12], xmm0

		mov			ecx, [edi+eax-2*12+0]
		imul		ecx, DRAWVERT_SIZE

		shufps		xmm0, xmm0, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm2, xmm2, R_SHUFFLEPS( 1, 2, 3, 0 )

		movss		[edx-2*16+0], xmm0
		movss		[edx-2*16+4], xmm1
		movss		[edx-2*16+8], xmm2

		mulss		xmm0, [esi+ecx+DRAWVERT_XYZ_OFFSET+0]
		mulss		xmm1, [esi+ecx+DRAWVERT_XYZ_OFFSET+4]
		mulss		xmm2, [esi+ecx+DRAWVERT_XYZ_OFFSET+8]

		xorps		xmm0, SIMD_SP_singleSignBitMask
		subss		xmm0, xmm1
		subss		xmm0, xmm2
		movss		[edx-2*16+12], xmm0

		mov			ecx, [edi+eax-3*12+0]
		imul		ecx, DRAWVERT_SIZE

		shufps		xmm0, xmm0, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm2, xmm2, R_SHUFFLEPS( 1, 2, 3, 0 )

		movss		[edx-3*16+0], xmm0
		movss		[edx-3*16+4], xmm1
		movss		[edx-3*16+8], xmm2

		mulss		xmm0, [esi+ecx+DRAWVERT_XYZ_OFFSET+0]
		mulss		xmm1, [esi+ecx+DRAWVERT_XYZ_OFFSET+4]
		mulss		xmm2, [esi+ecx+DRAWVERT_XYZ_OFFSET+8]

		xorps		xmm0, SIMD_SP_singleSignBitMask
		subss		xmm0, xmm1
		subss		xmm0, xmm2
		movss		[edx-3*16+12], xmm0

		mov			ecx, [edi+eax-4*12+0]
		imul		ecx, DRAWVERT_SIZE

		shufps		xmm0, xmm0, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm2, xmm2, R_SHUFFLEPS( 1, 2, 3, 0 )

		movss		[edx-4*16+0], xmm0
		movss		[edx-4*16+4], xmm1
		movss		[edx-4*16+8], xmm2

		mulss		xmm0, [esi+ecx+DRAWVERT_XYZ_OFFSET+0]
		mulss		xmm1, [esi+ecx+DRAWVERT_XYZ_OFFSET+4]
		mulss		xmm2, [esi+ecx+DRAWVERT_XYZ_OFFSET+8]

		xorps		xmm0, SIMD_SP_singleSignBitMask
		subss		xmm0, xmm1
		subss		xmm0, xmm2
		movss		[edx-4*16+12], xmm0

		add			eax, 4*12
		jle			loopPlane4

	done4:

		sub			eax, 4*12
		jge			done

	loopPlane1:
		mov			ebx, [edi+eax+4]
		imul		ebx, DRAWVERT_SIZE
		mov			ecx, [edi+eax+0]
		imul		ecx, DRAWVERT_SIZE

		movss		xmm0, [esi+ebx+DRAWVERT_XYZ_OFFSET+0]
		subss		xmm0, [esi+ecx+DRAWVERT_XYZ_OFFSET+0]

		movss		xmm1, [esi+ebx+DRAWVERT_XYZ_OFFSET+4]
		subss		xmm1, [esi+ecx+DRAWVERT_XYZ_OFFSET+4]

		movss		xmm2, [esi+ebx+DRAWVERT_XYZ_OFFSET+8]
		subss		xmm2, [esi+ecx+DRAWVERT_XYZ_OFFSET+8]

		mov			ebx, [edi+eax+8]
		imul		ebx, DRAWVERT_SIZE

		movss		xmm3, [esi+ebx+DRAWVERT_XYZ_OFFSET+0]
		subss		xmm3, [esi+ecx+DRAWVERT_XYZ_OFFSET+0]

		movss		xmm4, [esi+ebx+DRAWVERT_XYZ_OFFSET+4]
		subss		xmm4, [esi+ecx+DRAWVERT_XYZ_OFFSET+4]

		movss		xmm5, [esi+ebx+DRAWVERT_XYZ_OFFSET+8]
		subss		xmm5, [esi+ecx+DRAWVERT_XYZ_OFFSET+8]

		movss		xmm6, xmm4
		mulss		xmm6, xmm2
		movss		xmm7, xmm5
		mulss		xmm7, xmm1
		subss		xmm6, xmm7

		mulss		xmm5, xmm0
		mulss		xmm2, xmm3
		subss		xmm5, xmm2

		mulss		xmm3, xmm1
		mulss		xmm4, xmm0
		subss		xmm3, xmm4

		movss		xmm0, xmm6
		mulss		xmm6, xmm6
		movss		xmm1, xmm5
		mulss		xmm5, xmm5
		movss		xmm2, xmm3
		mulss		xmm3, xmm3

		addss		xmm3, xmm5
		addss		xmm3, xmm6
		rsqrtss		xmm3, xmm3

		add			edx, 1*16

		mulss		xmm0, xmm3
		mulss		xmm1, xmm3
		mulss		xmm2, xmm3

		movss		[edx-1*16+0], xmm0
		movss		[edx-1*16+4], xmm1
		movss		[edx-1*16+8], xmm2

		mulss		xmm0, [esi+ecx+DRAWVERT_XYZ_OFFSET+0]
		mulss		xmm1, [esi+ecx+DRAWVERT_XYZ_OFFSET+4]
		mulss		xmm2, [esi+ecx+DRAWVERT_XYZ_OFFSET+8]

		xorps		xmm0, SIMD_SP_singleSignBitMask
		subss		xmm0, xmm1
		subss		xmm0, xmm2
		movss		[edx-1*16+12], xmm0

		add			eax, 1*12
		jl			loopPlane1

	done:
	}

#else

	int i, j;

	for ( i = 0; i <= numIndexes - 12; i += 12 ) {
		ALIGN16( float d0[4] );
		ALIGN16( float d1[4] );
		ALIGN16( float d2[4] );
		ALIGN16( float d3[4] );
		ALIGN16( float d4[4] );
		ALIGN16( float d5[4] );
		ALIGN16( float n0[4] );
		ALIGN16( float n1[4] );
		ALIGN16( float n2[4] );

		for ( j = 0; j < 4; j++ ) {
			const CVertex *a, *b, *c;

			a = verts + indexes[i + j * 3 + 0];
			b = verts + indexes[i + j * 3 + 1];
			c = verts + indexes[i + j * 3 + 2];

			d0[j] = b->xyz[0] - a->xyz[0];
			d1[j] = b->xyz[1] - a->xyz[1];
			d2[j] = b->xyz[2] - a->xyz[2];

			d3[j] = c->xyz[0] - a->xyz[0];
			d4[j] = c->xyz[1] - a->xyz[1];
			d5[j] = c->xyz[2] - a->xyz[2];
		}

		ALIGN16( float tmp[4] );

		n0[0] = d4[0] * d2[0];
		n0[1] = d4[1] * d2[1];
		n0[2] = d4[2] * d2[2];
		n0[3] = d4[3] * d2[3];

		n0[0] -= d5[0] * d1[0];
		n0[1] -= d5[1] * d1[1];
		n0[2] -= d5[2] * d1[2];
		n0[3] -= d5[3] * d1[3];

		n1[0] = d5[0] * d0[0];
		n1[1] = d5[1] * d0[1];
		n1[2] = d5[2] * d0[2];
		n1[3] = d5[3] * d0[3];

		n1[0] -= d3[0] * d2[0];
		n1[1] -= d3[1] * d2[1];
		n1[2] -= d3[2] * d2[2];
		n1[3] -= d3[3] * d2[3];

		n2[0] = d3[0] * d1[0];
		n2[1] = d3[1] * d1[1];
		n2[2] = d3[2] * d1[2];
		n2[3] = d3[3] * d1[3];

		n2[0] -= d4[0] * d0[0];
		n2[1] -= d4[1] * d0[1];
		n2[2] -= d4[2] * d0[2];
		n2[3] -= d4[3] * d0[3];

		tmp[0] = n0[0] * n0[0];
		tmp[1] = n0[1] * n0[1];
		tmp[2] = n0[2] * n0[2];
		tmp[3] = n0[3] * n0[3];

		tmp[0] += n1[0] * n1[0];
		tmp[1] += n1[1] * n1[1];
		tmp[2] += n1[2] * n1[2];
		tmp[3] += n1[3] * n1[3];

		tmp[0] += n2[0] * n2[0];
		tmp[1] += n2[1] * n2[1];
		tmp[2] += n2[2] * n2[2];
		tmp[3] += n2[3] * n2[3];

		tmp[0] = CMath::rSqrt( tmp[0] );
		tmp[1] = CMath::rSqrt( tmp[1] );
		tmp[2] = CMath::rSqrt( tmp[2] );
		tmp[3] = CMath::rSqrt( tmp[3] );

		n0[0] *= tmp[0];
		n0[1] *= tmp[1];
		n0[2] *= tmp[2];
		n0[3] *= tmp[3];

		n1[0] *= tmp[0];
		n1[1] *= tmp[1];
		n1[2] *= tmp[2];
		n1[3] *= tmp[3];

		n2[0] *= tmp[0];
		n2[1] *= tmp[1];
		n2[2] *= tmp[2];
		n2[3] *= tmp[3];


		for ( j = 0; j < 4; j++ ) {
			const CVertex *a;

			a = verts + indexes[i + j * 3];

			planes->Normal()[0] = n0[j];
			planes->Normal()[1] = n1[j];
			planes->Normal()[2] = n2[j];
			planes->fitThroughPoint( a->xyz );
			planes++;
		}
	}

	for ( ; i < numIndexes; i += 3 ) {
		const CVertex *a, *b, *c;
		float d0, d1, d2, d3, d4, d5;
		float n0, n1, n2;

		a = verts + indexes[i + 0];
		b = verts + indexes[i + 1];
		c = verts + indexes[i + 2];

		d0 = b->xyz[0] - a->xyz[0];
		d1 = b->xyz[1] - a->xyz[1];
		d2 = b->xyz[2] - a->xyz[2];

		d3 = c->xyz[0] - a->xyz[0];
		d4 = c->xyz[1] - a->xyz[1];
		d5 = c->xyz[2] - a->xyz[2];

		float tmp;

		n0 = d4 * d2 - d5 * d1;
		n1 = d5 * d0 - d3 * d2;
		n2 = d3 * d1 - d4 * d0;

		tmp = CMath::rSqrt( n0 * n0 + n1 * n1 + n2 * n2 );

		n0 *= tmp;
		n1 *= tmp;
		n2 *= tmp;

		planes->Normal()[0] = n0;
		planes->Normal()[1] = n1;
		planes->Normal()[2] = n2;
		planes->fitThroughPoint( a->xyz );
		planes++;
	}

#endif
}

/*
============
CSIMD_SSE:: deriveTangents
============
*/
//#define REFINE_TANGENT_SQUAREROOT
#define FIX_DEGENERATE_TANGENT

void VPCALL CSIMD_SSE:: deriveTangents( CPlane *planes, CVertex *verts, const int numVerts, const int *indexes, const int numIndexes ) {
	int i;

	SMF_ASSERT( sizeof( CVertex ) == DRAWVERT_SIZE );
	SMF_ASSERT( (int)&((CVertex *)0)->normal == DRAWVERT_NORMAL_OFFSET );
	SMF_ASSERT( (int)&((CVertex *)0)->tangents[0] == DRAWVERT_TANGENT0_OFFSET );
	SMF_ASSERT( (int)&((CVertex *)0)->tangents[1] == DRAWVERT_TANGENT1_OFFSET );

	SMF_ASSERT( planes != NULL );
	SMF_ASSERT( verts != NULL );
	SMF_ASSERT( numVerts >= 0 );

#ifdef REFINE_TANGENT_SQUAREROOT
	__asm {
		movaps		xmm6, SIMD_SP_rsqrt_c0
		movaps		xmm7, SIMD_SP_rsqrt_c1
	}
#endif

	bool *used = (bool *)_alloca16( numVerts * sizeof( used[0] ) );
	memset( used, 0, numVerts * sizeof( used[0] ) );

	for ( i = 0; i <= numIndexes - 12; i += 12 ) {
		CVertex *a, *b, *c;
		ALIGN16( unsigned long signBit[4] );
		ALIGN16( float d0[4] );
		ALIGN16( float d1[4] );
		ALIGN16( float d2[4] );
		ALIGN16( float d3[4] );
		ALIGN16( float d4[4] );
		ALIGN16( float d5[4] );
		ALIGN16( float d6[4] );
		ALIGN16( float d7[4] );
		ALIGN16( float d8[4] );
		ALIGN16( float d9[4] );
		ALIGN16( float n0[4] );
		ALIGN16( float n1[4] );
		ALIGN16( float n2[4] );
		ALIGN16( float t0[4] );
		ALIGN16( float t1[4] );
		ALIGN16( float t2[4] );
		ALIGN16( float t3[4] );
		ALIGN16( float t4[4] );
		ALIGN16( float t5[4] );

		for ( int j = 0; j < 4; j++ ) {

			a = verts + indexes[i + j * 3 + 0];
			b = verts + indexes[i + j * 3 + 1];
			c = verts + indexes[i + j * 3 + 2];

			d0[j] = b->xyz[0] - a->xyz[0];
			d1[j] = b->xyz[1] - a->xyz[1];
			d2[j] = b->xyz[2] - a->xyz[2];
			d3[j] = b->st[0] - a->st[0];
			d4[j] = b->st[1] - a->st[1];

			d5[j] = c->xyz[0] - a->xyz[0];
			d6[j] = c->xyz[1] - a->xyz[1];
			d7[j] = c->xyz[2] - a->xyz[2];
			d8[j] = c->st[0] - a->st[0];
			d9[j] = c->st[1] - a->st[1];
		}

#if 1

		__asm {
			// normal
			movaps		xmm0, d6
			mulps		xmm0, d2
			movaps		xmm1, d7
			mulps		xmm1, d1
			subps		xmm0, xmm1

			movaps		xmm1, d7
			mulps		xmm1, d0
			movaps		xmm2, d5
			mulps		xmm2, d2
			subps		xmm1, xmm2

			movaps		xmm2, d5
			mulps		xmm2, d1
			movaps		xmm3, d6
			mulps		xmm3, d0
			subps		xmm2, xmm3

			movaps		xmm3, xmm0
			movaps		xmm4, xmm1
			movaps		xmm5, xmm2

			mulps		xmm3, xmm3
			mulps		xmm4, xmm4
			mulps		xmm5, xmm5

			addps		xmm3, xmm4
			addps		xmm3, xmm5

#ifdef FIX_DEGENERATE_TANGENT
			xorps		xmm4, xmm4
			cmpeqps		xmm4, xmm3
			andps		xmm4, SIMD_SP_tiny			// if values are zero replace them with a tiny number
			andps		xmm3, SIMD_SP_absMask		// make sure the values are positive
			orps		xmm3, xmm4
#endif

#ifdef REFINE_TANGENT_SQUAREROOT
			rsqrtps		xmm4, xmm3
			mulps		xmm3, xmm4
			mulps		xmm3, xmm4
			subps		xmm3, xmm6
			mulps		xmm4, xmm7
			mulps		xmm3, xmm4
#else
			rsqrtps		xmm3, xmm3
#endif
			mulps		xmm0, xmm3
			movaps		n0, xmm0
			mulps		xmm1, xmm3
			movaps		n1, xmm1
			mulps		xmm2, xmm3
			movaps		n2, xmm2

			// area sign bit
			movaps		xmm0, d3
			mulps		xmm0, d9
			movaps		xmm1, d4
			mulps		xmm1, d8
			subps		xmm0, xmm1
			andps		xmm0, SIMD_SP_signBitMask
			movaps		signBit, xmm0

			// first tangent
			movaps		xmm0, d0
			mulps		xmm0, d9
			movaps		xmm1, d4
			mulps		xmm1, d5
			subps		xmm0, xmm1

			movaps		xmm1, d1
			mulps		xmm1, d9
			movaps		xmm2, d4
			mulps		xmm2, d6
			subps		xmm1, xmm2

			movaps		xmm2, d2
			mulps		xmm2, d9
			movaps		xmm3, d4
			mulps		xmm3, d7
			subps		xmm2, xmm3

			movaps		xmm3, xmm0
			movaps		xmm4, xmm1
			movaps		xmm5, xmm2

			mulps		xmm3, xmm3
			mulps		xmm4, xmm4
			mulps		xmm5, xmm5

			addps		xmm3, xmm4
			addps		xmm3, xmm5

#ifdef FIX_DEGENERATE_TANGENT
			xorps		xmm4, xmm4
			cmpeqps		xmm4, xmm3
			andps		xmm4, SIMD_SP_tiny			// if values are zero replace them with a tiny number
			andps		xmm3, SIMD_SP_absMask		// make sure the values are positive
			orps		xmm3, xmm4
#endif

#ifdef REFINE_TANGENT_SQUAREROOT
			rsqrtps		xmm4, xmm3
			mulps		xmm3, xmm4
			mulps		xmm3, xmm4
			subps		xmm3, xmm6
			mulps		xmm4, xmm7
			mulps		xmm3, xmm4
#else
			rsqrtps		xmm3, xmm3
#endif
			xorps		xmm3, signBit

			mulps		xmm0, xmm3
			movaps		t0, xmm0
			mulps		xmm1, xmm3
			movaps		t1, xmm1
			mulps		xmm2, xmm3
			movaps		t2, xmm2

			// second tangent
			movaps		xmm0, d3
			mulps		xmm0, d5
			movaps		xmm1, d0
			mulps		xmm1, d8
			subps		xmm0, xmm1

			movaps		xmm1, d3
			mulps		xmm1, d6
			movaps		xmm2, d1
			mulps		xmm2, d8
			subps		xmm1, xmm2

			movaps		xmm2, d3
			mulps		xmm2, d7
			movaps		xmm3, d2
			mulps		xmm3, d8
			subps		xmm2, xmm3

			movaps		xmm3, xmm0
			movaps		xmm4, xmm1
			movaps		xmm5, xmm2

			mulps		xmm3, xmm3
			mulps		xmm4, xmm4
			mulps		xmm5, xmm5

			addps		xmm3, xmm4
			addps		xmm3, xmm5

#ifdef FIX_DEGENERATE_TANGENT
			xorps		xmm4, xmm4
			cmpeqps		xmm4, xmm3
			andps		xmm4, SIMD_SP_tiny			// if values are zero replace them with a tiny number
			andps		xmm3, SIMD_SP_absMask		// make sure the values are positive
			orps		xmm3, xmm4
#endif

#ifdef REFINE_TANGENT_SQUAREROOT
			rsqrtps		xmm4, xmm3
			mulps		xmm3, xmm4
			mulps		xmm3, xmm4
			subps		xmm3, xmm6
			mulps		xmm4, xmm7
			mulps		xmm3, xmm4
#else
			rsqrtps		xmm3, xmm3
#endif
			xorps		xmm3, signBit

			mulps		xmm0, xmm3
			movaps		t3, xmm0
			mulps		xmm1, xmm3
			movaps		t4, xmm1
			mulps		xmm2, xmm3
			movaps		t5, xmm2
		}

#else

		ALIGN16( float tmp[4] );

		// normal
		n0[0] = d6[0] * d2[0];
		n0[1] = d6[1] * d2[1];
		n0[2] = d6[2] * d2[2];
		n0[3] = d6[3] * d2[3];

		n0[0] -= d7[0] * d1[0];
		n0[1] -= d7[1] * d1[1];
		n0[2] -= d7[2] * d1[2];
		n0[3] -= d7[3] * d1[3];

		n1[0] = d7[0] * d0[0];
		n1[1] = d7[1] * d0[1];
		n1[2] = d7[2] * d0[2];
		n1[3] = d7[3] * d0[3];

		n1[0] -= d5[0] * d2[0];
		n1[1] -= d5[1] * d2[1];
		n1[2] -= d5[2] * d2[2];
		n1[3] -= d5[3] * d2[3];

		n2[0] = d5[0] * d1[0];
		n2[1] = d5[1] * d1[1];
		n2[2] = d5[2] * d1[2];
		n2[3] = d5[3] * d1[3];

		n2[0] -= d6[0] * d0[0];
		n2[1] -= d6[1] * d0[1];
		n2[2] -= d6[2] * d0[2];
		n2[3] -= d6[3] * d0[3];

		tmp[0] = n0[0] * n0[0];
		tmp[1] = n0[1] * n0[1];
		tmp[2] = n0[2] * n0[2];
		tmp[3] = n0[3] * n0[3];

		tmp[0] += n1[0] * n1[0];
		tmp[1] += n1[1] * n1[1];
		tmp[2] += n1[2] * n1[2];
		tmp[3] += n1[3] * n1[3];

		tmp[0] += n2[0] * n2[0];
		tmp[1] += n2[1] * n2[1];
		tmp[2] += n2[2] * n2[2];
		tmp[3] += n2[3] * n2[3];

		tmp[0] = CMath::rSqrt( tmp[0] );
		tmp[1] = CMath::rSqrt( tmp[1] );
		tmp[2] = CMath::rSqrt( tmp[2] );
		tmp[3] = CMath::rSqrt( tmp[3] );

		n0[0] *= tmp[0];
		n0[1] *= tmp[1];
		n0[2] *= tmp[2];
		n0[3] *= tmp[3];

		n1[0] *= tmp[0];
		n1[1] *= tmp[1];
		n1[2] *= tmp[2];
		n1[3] *= tmp[3];

		n2[0] *= tmp[0];
		n2[1] *= tmp[1];
		n2[2] *= tmp[2];
		n2[3] *= tmp[3];

		// area sign bit
		tmp[0] = d3[0] * d9[0];
		tmp[1] = d3[1] * d9[1];
		tmp[2] = d3[2] * d9[2];
		tmp[3] = d3[3] * d9[3];

		tmp[0] -= d4[0] * d8[0];
		tmp[1] -= d4[1] * d8[1];
		tmp[2] -= d4[2] * d8[2];
		tmp[3] -= d4[3] * d8[3];

		signBit[0] = ( *(unsigned long *)&tmp[0] ) & ( 1 << 31 );
		signBit[1] = ( *(unsigned long *)&tmp[1] ) & ( 1 << 31 );
		signBit[2] = ( *(unsigned long *)&tmp[2] ) & ( 1 << 31 );
		signBit[3] = ( *(unsigned long *)&tmp[3] ) & ( 1 << 31 );

		// first tangent
		t0[0] = d0[0] * d9[0];
		t0[1] = d0[1] * d9[1];
		t0[2] = d0[2] * d9[2];
		t0[3] = d0[3] * d9[3];

		t0[0] -= d4[0] * d5[0];
		t0[1] -= d4[1] * d5[1];
		t0[2] -= d4[2] * d5[2];
		t0[3] -= d4[3] * d5[3];

		t1[0] = d1[0] * d9[0];
		t1[1] = d1[1] * d9[1];
		t1[2] = d1[2] * d9[2];
		t1[3] = d1[3] * d9[3];

		t1[0] -= d4[0] * d6[0];
		t1[1] -= d4[1] * d6[1];
		t1[2] -= d4[2] * d6[2];
		t1[3] -= d4[3] * d6[3];

		t2[0] = d2[0] * d9[0];
		t2[1] = d2[1] * d9[1];
		t2[2] = d2[2] * d9[2];
		t2[3] = d2[3] * d9[3];

		t2[0] -= d4[0] * d7[0];
		t2[1] -= d4[1] * d7[1];
		t2[2] -= d4[2] * d7[2];
		t2[3] -= d4[3] * d7[3];

		tmp[0] = t0[0] * t0[0];
		tmp[1] = t0[1] * t0[1];
		tmp[2] = t0[2] * t0[2];
		tmp[3] = t0[3] * t0[3];

		tmp[0] += t1[0] * t1[0];
		tmp[1] += t1[1] * t1[1];
		tmp[2] += t1[2] * t1[2];
		tmp[3] += t1[3] * t1[3];

		tmp[0] += t2[0] * t2[0];
		tmp[1] += t2[1] * t2[1];
		tmp[2] += t2[2] * t2[2];
		tmp[3] += t2[3] * t2[3];

		tmp[0] = CMath::rSqrt( tmp[0] );
		tmp[1] = CMath::rSqrt( tmp[1] );
		tmp[2] = CMath::rSqrt( tmp[2] );
		tmp[3] = CMath::rSqrt( tmp[3] );

		*(unsigned long *)&tmp[0] ^= signBit[0];
		*(unsigned long *)&tmp[1] ^= signBit[1];
		*(unsigned long *)&tmp[2] ^= signBit[2];
		*(unsigned long *)&tmp[3] ^= signBit[3];

		t0[0] *= tmp[0];
		t0[1] *= tmp[1];
		t0[2] *= tmp[2];
		t0[3] *= tmp[3];

		t1[0] *= tmp[0];
		t1[1] *= tmp[1];
		t1[2] *= tmp[2];
		t1[3] *= tmp[3];

		t2[0] *= tmp[0];
		t2[1] *= tmp[1];
		t2[2] *= tmp[2];
		t2[3] *= tmp[3];

		// second tangent
		t3[0] = d3[0] * d5[0];
		t3[1] = d3[1] * d5[1];
		t3[2] = d3[2] * d5[2];
		t3[3] = d3[3] * d5[3];

		t3[0] -= d0[0] * d8[0];
		t3[1] -= d0[1] * d8[1];
		t3[2] -= d0[2] * d8[2];
		t3[3] -= d0[3] * d8[3];

		t4[0] = d3[0] * d6[0];
		t4[1] = d3[1] * d6[1];
		t4[2] = d3[2] * d6[2];
		t4[3] = d3[3] * d6[3];

		t4[0] -= d1[0] * d8[0];
		t4[1] -= d1[1] * d8[1];
		t4[2] -= d1[2] * d8[2];
		t4[3] -= d1[3] * d8[3];

		t5[0] = d3[0] * d7[0];
		t5[1] = d3[1] * d7[1];
		t5[2] = d3[2] * d7[2];
		t5[3] = d3[3] * d7[3];

		t5[0] -= d2[0] * d8[0];
		t5[1] -= d2[1] * d8[1];
		t5[2] -= d2[2] * d8[2];
		t5[3] -= d2[3] * d8[3];

		tmp[0] = t3[0] * t3[0];
		tmp[1] = t3[1] * t3[1];
		tmp[2] = t3[2] * t3[2];
		tmp[3] = t3[3] * t3[3];

		tmp[0] += t4[0] * t4[0];
		tmp[1] += t4[1] * t4[1];
		tmp[2] += t4[2] * t4[2];
		tmp[3] += t4[3] * t4[3];

		tmp[0] += t5[0] * t5[0];
		tmp[1] += t5[1] * t5[1];
		tmp[2] += t5[2] * t5[2];
		tmp[3] += t5[3] * t5[3];

		tmp[0] = CMath::rSqrt( tmp[0] );
		tmp[1] = CMath::rSqrt( tmp[1] );
		tmp[2] = CMath::rSqrt( tmp[2] );
		tmp[3] = CMath::rSqrt( tmp[3] );

		*(unsigned long *)&tmp[0] ^= signBit[0];
		*(unsigned long *)&tmp[1] ^= signBit[1];
		*(unsigned long *)&tmp[2] ^= signBit[2];
		*(unsigned long *)&tmp[3] ^= signBit[3];

		t3[0] *= tmp[0];
		t3[1] *= tmp[1];
		t3[2] *= tmp[2];
		t3[3] *= tmp[3];

		t4[0] *= tmp[0];
		t4[1] *= tmp[1];
		t4[2] *= tmp[2];
		t4[3] *= tmp[3];

		t5[0] *= tmp[0];
		t5[1] *= tmp[1];
		t5[2] *= tmp[2];
		t5[3] *= tmp[3];

#endif

		for ( int j = 0; j < 4; j++ ) {

			const int v0 = indexes[i + j * 3 + 0];
			const int v1 = indexes[i + j * 3 + 1];
			const int v2 = indexes[i + j * 3 + 2];

			a = verts + v0;
			b = verts + v1;
			c = verts + v2;

			planes->getNormal()[0] = n0[j];
			planes->getNormal()[1] = n1[j];
			planes->getNormal()[2] = n2[j];
			planes->fitThroughPoint( a->xyz );
			planes++;

			if ( used[v0] ) {
				a->normal[0] += n0[j];
				a->normal[1] += n1[j];
				a->normal[2] += n2[j];

				a->tangents[0][0] += t0[j];
				a->tangents[0][1] += t1[j];
				a->tangents[0][2] += t2[j];

				a->tangents[1][0] += t3[j];
				a->tangents[1][1] += t4[j];
				a->tangents[1][2] += t5[j];
			} else {
				a->normal[0] = n0[j];
				a->normal[1] = n1[j];
				a->normal[2] = n2[j];

				a->tangents[0][0] = t0[j];
				a->tangents[0][1] = t1[j];
				a->tangents[0][2] = t2[j];

				a->tangents[1][0] = t3[j];
				a->tangents[1][1] = t4[j];
				a->tangents[1][2] = t5[j];

				used[v0] = true;
			}

			if ( used[v1] ) {
				b->normal[0] += n0[j];
				b->normal[1] += n1[j];
				b->normal[2] += n2[j];

				b->tangents[0][0] += t0[j];
				b->tangents[0][1] += t1[j];
				b->tangents[0][2] += t2[j];

				b->tangents[1][0] += t3[j];
				b->tangents[1][1] += t4[j];
				b->tangents[1][2] += t5[j];
			} else {
				b->normal[0] = n0[j];
				b->normal[1] = n1[j];
				b->normal[2] = n2[j];

				b->tangents[0][0] = t0[j];
				b->tangents[0][1] = t1[j];
				b->tangents[0][2] = t2[j];

				b->tangents[1][0] = t3[j];
				b->tangents[1][1] = t4[j];
				b->tangents[1][2] = t5[j];

				used[v1] = true;
			}

			if ( used[v2] ) {
				c->normal[0] += n0[j];
				c->normal[1] += n1[j];
				c->normal[2] += n2[j];

				c->tangents[0][0] += t0[j];
				c->tangents[0][1] += t1[j];
				c->tangents[0][2] += t2[j];

				c->tangents[1][0] += t3[j];
				c->tangents[1][1] += t4[j];
				c->tangents[1][2] += t5[j];
			} else {
				c->normal[0] = n0[j];
				c->normal[1] = n1[j];
				c->normal[2] = n2[j];

				c->tangents[0][0] = t0[j];
				c->tangents[0][1] = t1[j];
				c->tangents[0][2] = t2[j];

				c->tangents[1][0] = t3[j];
				c->tangents[1][1] = t4[j];
				c->tangents[1][2] = t5[j];

				used[v2] = true;
			}
		}
	}

	for ( ; i < numIndexes; i += 3 ) {
		CVertex *a, *b, *c;
		ALIGN16( unsigned long signBit[4] );
		float d0, d1, d2, d3, d4;
		float d5, d6, d7, d8, d9;
		float n0, n1, n2;
		float t0, t1, t2;
		float t3, t4, t5;

		const int v0 = indexes[i + 0];
		const int v1 = indexes[i + 1];
		const int v2 = indexes[i + 2];

		a = verts + v0;
		b = verts + v1;
		c = verts + v2;

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

#if 1

		__asm {
			// normal
			movss		xmm0, d6
			mulss		xmm0, d2
			movss		xmm1, d7
			mulss		xmm1, d1
			subss		xmm0, xmm1

			movss		xmm1, d7
			mulss		xmm1, d0
			movss		xmm2, d5
			mulss		xmm2, d2
			subss		xmm1, xmm2

			movss		xmm2, d5
			mulss		xmm2, d1
			movss		xmm3, d6
			mulss		xmm3, d0
			subss		xmm2, xmm3

			movss		xmm3, xmm0
			movss		xmm4, xmm1
			movss		xmm5, xmm2

			mulss		xmm3, xmm3
			mulss		xmm4, xmm4
			mulss		xmm5, xmm5

			addss		xmm3, xmm4
			addss		xmm3, xmm5

#ifdef FIX_DEGENERATE_TANGENT
			xorps		xmm4, xmm4
			cmpeqps		xmm4, xmm3
			andps		xmm4, SIMD_SP_tiny			// if values are zero replace them with a tiny number
			andps		xmm3, SIMD_SP_absMask		// make sure the values are positive
			orps		xmm3, xmm4
#endif

#ifdef REFINE_TANGENT_SQUAREROOT
			rsqrtss		xmm4, xmm3
			mulss		xmm3, xmm4
			mulss		xmm3, xmm4
			subss		xmm3, xmm6
			mulss		xmm4, xmm7
			mulss		xmm3, xmm4
#else
			rsqrtss		xmm3, xmm3
#endif
			mulss		xmm0, xmm3
			movss		n0, xmm0
			mulss		xmm1, xmm3
			movss		n1, xmm1
			mulss		xmm2, xmm3
			movss		n2, xmm2

			// area sign bit
			movss		xmm0, d3
			mulss		xmm0, d9
			movss		xmm1, d4
			mulss		xmm1, d8
			subss		xmm0, xmm1
			andps		xmm0, SIMD_SP_signBitMask
			movaps		signBit, xmm0

			// first tangent
			movss		xmm0, d0
			mulss		xmm0, d9
			movss		xmm1, d4
			mulss		xmm1, d5
			subss		xmm0, xmm1

			movss		xmm1, d1
			mulss		xmm1, d9
			movss		xmm2, d4
			mulss		xmm2, d6
			subss		xmm1, xmm2

			movss		xmm2, d2
			mulss		xmm2, d9
			movss		xmm3, d4
			mulss		xmm3, d7
			subss		xmm2, xmm3

			movss		xmm3, xmm0
			movss		xmm4, xmm1
			movss		xmm5, xmm2

			mulss		xmm3, xmm3
			mulss		xmm4, xmm4
			mulss		xmm5, xmm5

			addss		xmm3, xmm4
			addss		xmm3, xmm5

#ifdef FIX_DEGENERATE_TANGENT
			xorps		xmm4, xmm4
			cmpeqps		xmm4, xmm3
			andps		xmm4, SIMD_SP_tiny			// if values are zero replace them with a tiny number
			andps		xmm3, SIMD_SP_absMask		// make sure the values are positive
			orps		xmm3, xmm4
#endif

#ifdef REFINE_TANGENT_SQUAREROOT
			rsqrtss		xmm4, xmm3
			mulss		xmm3, xmm4
			mulss		xmm3, xmm4
			subss		xmm3, xmm6
			mulss		xmm4, xmm7
			mulss		xmm3, xmm4
#else
			rsqrtss		xmm3, xmm3
#endif
			xorps		xmm3, signBit

			mulss		xmm0, xmm3
			movss		t0, xmm0
			mulss		xmm1, xmm3
			movss		t1, xmm1
			mulss		xmm2, xmm3
			movss		t2, xmm2

			// second tangent
			movss		xmm0, d3
			mulss		xmm0, d5
			movss		xmm1, d0
			mulss		xmm1, d8
			subss		xmm0, xmm1

			movss		xmm1, d3
			mulss		xmm1, d6
			movss		xmm2, d1
			mulss		xmm2, d8
			subss		xmm1, xmm2

			movss		xmm2, d3
			mulss		xmm2, d7
			movss		xmm3, d2
			mulss		xmm3, d8
			subss		xmm2, xmm3

			movss		xmm3, xmm0
			movss		xmm4, xmm1
			movss		xmm5, xmm2

			mulss		xmm3, xmm3
			mulss		xmm4, xmm4
			mulss		xmm5, xmm5

			addss		xmm3, xmm4
			addss		xmm3, xmm5

#ifdef FIX_DEGENERATE_TANGENT
			xorps		xmm4, xmm4
			cmpeqps		xmm4, xmm3
			andps		xmm4, SIMD_SP_tiny			// if values are zero replace them with a tiny number
			andps		xmm3, SIMD_SP_absMask		// make sure the values are positive
			orps		xmm3, xmm4
#endif

#ifdef REFINE_TANGENT_SQUAREROOT
			rsqrtss		xmm4, xmm3
			mulss		xmm3, xmm4
			mulss		xmm3, xmm4
			subss		xmm3, xmm6
			mulss		xmm4, xmm7
			mulss		xmm3, xmm4
#else
			rsqrtss		xmm3, xmm3
#endif
			xorps		xmm3, signBit

			mulss		xmm0, xmm3
			movss		t3, xmm0
			mulss		xmm1, xmm3
			movss		t4, xmm1
			mulss		xmm2, xmm3
			movss		t5, xmm2
		}

#else

		float tmp;

		// normal
		n0 = d6 * d2 - d7 * d1;
		n1 = d7 * d0 - d5 * d2;
		n2 = d5 * d1 - d6 * d0;

		tmp = CMath::rSqrt( n0 * n0 + n1 * n1 + n2 * n2 );

		n0 *= tmp;
		n1 *= tmp;
		n2 *= tmp;

		// area sign bit
		tmp = d3 * d9 - d4 * d8;
		signBit[0] = ( *(unsigned long *)&tmp ) & ( 1 << 31 );

		// first tangent
		t0 = d0 * d9 - d4 * d5;
		t1 = d1 * d9 - d4 * d6;
		t2 = d2 * d9 - d4 * d7;

		tmp = CMath::rSqrt( t0 * t0 + t1 * t1 + t2 * t2 );
		*(unsigned long *)&tmp ^= signBit[0];

		t0 *= tmp;
		t1 *= tmp;
		t2 *= tmp;

		// second tangent
		t3 = d3 * d5 - d0 * d8;
		t4 = d3 * d6 - d1 * d8;
		t5 = d3 * d7 - d2 * d8;

		tmp = CMath::rSqrt( t3 * t3 + t4 * t4 + t5 * t5 );
		*(unsigned long *)&tmp ^= signBit[0];

		t3 *= tmp;
		t4 *= tmp;
		t5 *= tmp;

#endif

		planes->getNormal()[0] = n0;
		planes->getNormal()[1] = n1;
		planes->getNormal()[2] = n2;
		planes->fitThroughPoint( a->xyz );
		planes++;

		if ( used[v0] ) {
			a->normal[0] += n0;
			a->normal[1] += n1;
			a->normal[2] += n2;

			a->tangents[0][0] += t0;
			a->tangents[0][1] += t1;
			a->tangents[0][2] += t2;

			a->tangents[1][0] += t3;
			a->tangents[1][1] += t4;
			a->tangents[1][2] += t5;
		} else {
			a->normal[0] = n0;
			a->normal[1] = n1;
			a->normal[2] = n2;

			a->tangents[0][0] = t0;
			a->tangents[0][1] = t1;
			a->tangents[0][2] = t2;

			a->tangents[1][0] = t3;
			a->tangents[1][1] = t4;
			a->tangents[1][2] = t5;

			used[v0] = true;
		}

		if ( used[v1] ) {
			b->normal[0] += n0;
			b->normal[1] += n1;
			b->normal[2] += n2;

			b->tangents[0][0] += t0;
			b->tangents[0][1] += t1;
			b->tangents[0][2] += t2;

			b->tangents[1][0] += t3;
			b->tangents[1][1] += t4;
			b->tangents[1][2] += t5;
		} else {
			b->normal[0] = n0;
			b->normal[1] = n1;
			b->normal[2] = n2;

			b->tangents[0][0] = t0;
			b->tangents[0][1] = t1;
			b->tangents[0][2] = t2;

			b->tangents[1][0] = t3;
			b->tangents[1][1] = t4;
			b->tangents[1][2] = t5;

			used[v1] = true;
		}

		if ( used[v2] ) {
			c->normal[0] += n0;
			c->normal[1] += n1;
			c->normal[2] += n2;

			c->tangents[0][0] += t0;
			c->tangents[0][1] += t1;
			c->tangents[0][2] += t2;

			c->tangents[1][0] += t3;
			c->tangents[1][1] += t4;
			c->tangents[1][2] += t5;
		} else {
			c->normal[0] = n0;
			c->normal[1] = n1;
			c->normal[2] = n2;

			c->tangents[0][0] = t0;
			c->tangents[0][1] = t1;
			c->tangents[0][2] = t2;

			c->tangents[1][0] = t3;
			c->tangents[1][1] = t4;
			c->tangents[1][2] = t5;

			used[v2] = true;
		}
	}
}

/*
============
CSIMD_SSE::DeriveUnsmoothedTangents
============
*/
#define DERIVE_UNSMOOTHED_BITANGENT
/* s
void VPCALL CSIMD_SSE::DeriveUnsmoothedTangents( CVertex *verts, const dominantTri_s *dominantTris, const int numVerts ) {
	int i, j;

	for ( i = 0; i <= numVerts - 4; i += 4 ) {
		ALIGN16( float s0[4] );
		ALIGN16( float s1[4] );
		ALIGN16( float s2[4] );
		ALIGN16( float d0[4] );
		ALIGN16( float d1[4] );
		ALIGN16( float d2[4] );
		ALIGN16( float d3[4] );
		ALIGN16( float d4[4] );
		ALIGN16( float d5[4] );
		ALIGN16( float d6[4] );
		ALIGN16( float d7[4] );
		ALIGN16( float d8[4] );
		ALIGN16( float d9[4] );
		ALIGN16( float n0[4] );
		ALIGN16( float n1[4] );
		ALIGN16( float n2[4] );
		ALIGN16( float t0[4] );
		ALIGN16( float t1[4] );
		ALIGN16( float t2[4] );
		ALIGN16( float t3[4] );
		ALIGN16( float t4[4] );
		ALIGN16( float t5[4] );

		for ( j = 0; j < 4; j++ ) {
			const CVertex *a, *b, *c;

			const dominantTri_s &dt = dominantTris[i+j];

			s0[j] = dt.normalizationScale[0];
			s1[j] = dt.normalizationScale[1];
			s2[j] = dt.normalizationScale[2];

			a = verts + i + j;
			b = verts + dt.v2;
			c = verts + dt.v3;

			d0[j] = b->xyz[0] - a->xyz[0];
			d1[j] = b->xyz[1] - a->xyz[1];
			d2[j] = b->xyz[2] - a->xyz[2];
			d3[j] = b->st[0] - a->st[0];
			d4[j] = b->st[1] - a->st[1];

			d5[j] = c->xyz[0] - a->xyz[0];
			d6[j] = c->xyz[1] - a->xyz[1];
			d7[j] = c->xyz[2] - a->xyz[2];
			d8[j] = c->st[0] - a->st[0];
			d9[j] = c->st[1] - a->st[1];
		}

#if 1

		__asm {

			movaps		xmm0, d6
			mulps		xmm0, d2
			movaps		xmm1, d7
			mulps		xmm1, d1

			movaps		xmm2, d7
			mulps		xmm2, d0
			movaps		xmm3, d5
			mulps		xmm3, d2

			movaps		xmm4, d5
			mulps		xmm4, d1
			movaps		xmm5, d6
			mulps		xmm5, d0

			subps		xmm0, xmm1
			subps		xmm2, xmm3
			movaps		xmm7, s2
			subps		xmm4, xmm5

			mulps		xmm0, xmm7
			movaps		n0, xmm0
			mulps		xmm2, xmm7
			movaps		n1, xmm2
			mulps		xmm4, xmm7
			movaps		n2, xmm4

			movaps		xmm0, d0
			mulps		xmm0, d9
			movaps		xmm1, d4
			mulps		xmm1, d5

			movaps		xmm2, d1
			mulps		xmm2, d9
			movaps		xmm3, d4
			mulps		xmm3, d6

			movaps		xmm4, d2
			mulps		xmm4, d9
			movaps		xmm5, d4
			mulps		xmm5, d7

			subps		xmm0, xmm1
			subps		xmm2, xmm3
			movaps		xmm7, s0
			subps		xmm4, xmm5

			mulps		xmm0, xmm7
			movaps		t0, xmm0
			mulps		xmm2, xmm7
			movaps		t1, xmm2
			mulps		xmm4, xmm7
			movaps		t2, xmm4

#ifndef DERIVE_UNSMOOTHED_BITANGENT
			movaps		xmm0, d3
			mulps		xmm0, d5
			movaps		xmm1, d0
			mulps		xmm1, d8

			movaps		xmm2, d3
			mulps		xmm2, d6
			movaps		xmm3, d1
			mulps		xmm3, d8

			movaps		xmm4, d3
			mulps		xmm4, d7
			movaps		xmm5, d2
			mulps		xmm5, d8
#else
			movaps		xmm0, n2
			mulps		xmm0, t1
			movaps		xmm1, n1
			mulps		xmm1, t2

			movaps		xmm2, n0
			mulps		xmm2, t2
			movaps		xmm3, n2
			mulps		xmm3, t0

			movaps		xmm4, n1
			mulps		xmm4, t0
			movaps		xmm5, n0
			mulps		xmm5, t1
#endif
			subps		xmm0, xmm1
			subps		xmm2, xmm3
			movaps		xmm7, s1
			subps		xmm4, xmm5

			mulps		xmm0, xmm7
			movaps		t3, xmm0
			mulps		xmm2, xmm7
			movaps		t4, xmm2
			mulps		xmm4, xmm7
			movaps		t5, xmm4
		}

#else

		n0[0] = d6[0] * d2[0];
		n0[1] = d6[1] * d2[1];
		n0[2] = d6[2] * d2[2];
		n0[3] = d6[3] * d2[3];

		n1[0] = d7[0] * d0[0];
		n1[1] = d7[1] * d0[1];
		n1[2] = d7[2] * d0[2];
		n1[3] = d7[3] * d0[3];

		n2[0] = d5[0] * d1[0];
		n2[1] = d5[1] * d1[1];
		n2[2] = d5[2] * d1[2];
		n2[3] = d5[3] * d1[3];

		n0[0] -= d7[0] * d1[0];
		n0[1] -= d7[1] * d1[1];
		n0[2] -= d7[2] * d1[2];
		n0[3] -= d7[3] * d1[3];

		n1[0] -= d5[0] * d2[0];
		n1[1] -= d5[1] * d2[1];
		n1[2] -= d5[2] * d2[2];
		n1[3] -= d5[3] * d2[3];

		n2[0] -= d6[0] * d0[0];
		n2[1] -= d6[1] * d0[1];
		n2[2] -= d6[2] * d0[2];
		n2[3] -= d6[3] * d0[3];

		n0[0] *= s2[0];
		n0[1] *= s2[1];
		n0[2] *= s2[2];
		n0[3] *= s2[3];

		n1[0] *= s2[0];
		n1[1] *= s2[1];
		n1[2] *= s2[2];
		n1[3] *= s2[3];

		n2[0] *= s2[0];
		n2[1] *= s2[1];
		n2[2] *= s2[2];
		n2[3] *= s2[3];

		t0[0] = d0[0] * d9[0];
		t0[1] = d0[1] * d9[1];
		t0[2] = d0[2] * d9[2];
		t0[3] = d0[3] * d9[3];

		t1[0] = d1[0] * d9[0];
		t1[1] = d1[1] * d9[1];
		t1[2] = d1[2] * d9[2];
		t1[3] = d1[3] * d9[3];

		t2[0] = d2[0] * d9[0];
		t2[1] = d2[1] * d9[1];
		t2[2] = d2[2] * d9[2];
		t2[3] = d2[3] * d9[3];

		t0[0] -= d4[0] * d5[0];
		t0[1] -= d4[1] * d5[1];
		t0[2] -= d4[2] * d5[2];
		t0[3] -= d4[3] * d5[3];

		t1[0] -= d4[0] * d6[0];
		t1[1] -= d4[1] * d6[1];
		t1[2] -= d4[2] * d6[2];
		t1[3] -= d4[3] * d6[3];

		t2[0] -= d4[0] * d7[0];
		t2[1] -= d4[1] * d7[1];
		t2[2] -= d4[2] * d7[2];
		t2[3] -= d4[3] * d7[3];

		t0[0] *= s0[0];
		t0[1] *= s0[1];
		t0[2] *= s0[2];
		t0[3] *= s0[3];

		t1[0] *= s0[0];
		t1[1] *= s0[1];
		t1[2] *= s0[2];
		t1[3] *= s0[3];

		t2[0] *= s0[0];
		t2[1] *= s0[1];
		t2[2] *= s0[2];
		t2[3] *= s0[3];

#ifndef DERIVE_UNSMOOTHED_BITANGENT
		t3[0] = d3[0] * d5[0];
		t3[1] = d3[1] * d5[1];
		t3[2] = d3[2] * d5[2];
		t3[3] = d3[3] * d5[3];

		t4[0] = d3[0] * d6[0];
		t4[1] = d3[1] * d6[1];
		t4[2] = d3[2] * d6[2];
		t4[3] = d3[3] * d6[3];

		t5[0] = d3[0] * d7[0];
		t5[1] = d3[1] * d7[1];
		t5[2] = d3[2] * d7[2];
		t5[3] = d3[3] * d7[3];

		t3[0] -= d0[0] * d8[0];
		t3[1] -= d0[1] * d8[1];
		t3[2] -= d0[2] * d8[2];
		t3[3] -= d0[3] * d8[3];

		t4[0] -= d1[0] * d8[0];
		t4[1] -= d1[1] * d8[1];
		t4[2] -= d1[2] * d8[2];
		t4[3] -= d1[3] * d8[3];

		t5[0] -= d2[0] * d8[0];
		t5[1] -= d2[1] * d8[1];
		t5[2] -= d2[2] * d8[2];
		t5[3] -= d2[3] * d8[3];
#else
		t3[0] = n2[0] * t1[0];
		t3[1] = n2[1] * t1[1];
		t3[2] = n2[2] * t1[2];
		t3[3] = n2[3] * t1[3];

		t4[0] = n0[0] * t2[0];
		t4[1] = n0[1] * t2[1];
		t4[2] = n0[2] * t2[2];
		t4[3] = n0[3] * t2[3];

		t5[0] = n1[0] * t0[0];
		t5[1] = n1[1] * t0[1];
		t5[2] = n1[2] * t0[2];
		t5[3] = n1[3] * t0[3];

		t3[0] -= n1[0] * t2[0];
		t3[1] -= n1[1] * t2[1];
		t3[2] -= n1[2] * t2[2];
		t3[3] -= n1[3] * t2[3];

		t4[0] -= n2[0] * t0[0];
		t4[1] -= n2[1] * t0[1];
		t4[2] -= n2[2] * t0[2];
		t4[3] -= n2[3] * t0[3];

		t5[0] -= n0[0] * t1[0];
		t5[1] -= n0[1] * t1[1];
		t5[2] -= n0[2] * t1[2];
		t5[3] -= n0[3] * t1[3];
#endif
		t3[0] *= s1[0];
		t3[1] *= s1[1];
		t3[2] *= s1[2];
		t3[3] *= s1[3];

		t4[0] *= s1[0];
		t4[1] *= s1[1];
		t4[2] *= s1[2];
		t4[3] *= s1[3];

		t5[0] *= s1[0];
		t5[1] *= s1[1];
		t5[2] *= s1[2];
		t5[3] *= s1[3];

#endif

		for ( j = 0; j < 4; j++ ) {
			CVertex *a;

			a = verts + i + j;

			a->normal[0] = n0[j];
			a->normal[1] = n1[j];
			a->normal[2] = n2[j];

			a->tangents[0][0] = t0[j];
			a->tangents[0][1] = t1[j];
			a->tangents[0][2] = t2[j];

			a->tangents[1][0] = t3[j];
			a->tangents[1][1] = t4[j];
			a->tangents[1][2] = t5[j];
		}
	}

	for ( ; i < numVerts; i++ ) {
		CVertex *a, *b, *c;
		float d0, d1, d2, d3, d4;
		float d5, d6, d7, d8, d9;
		float s0, s1, s2;
		float n0, n1, n2;
		float t0, t1, t2;
		float t3, t4, t5;

		const dominantTri_s &dt = dominantTris[i];

		s0 = dt.normalizationScale[0];
		s1 = dt.normalizationScale[1];
		s2 = dt.normalizationScale[2];

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

#if 1

		__asm {

			movss		xmm0, d6
			mulss		xmm0, d2
			movss		xmm1, d7
			mulss		xmm1, d1

			movss		xmm2, d7
			mulss		xmm2, d0
			movss		xmm3, d5
			mulss		xmm3, d2

			movss		xmm4, d5
			mulss		xmm4, d1
			movss		xmm5, d6
			mulss		xmm5, d0

			subss		xmm0, xmm1
			subss		xmm2, xmm3
			movss		xmm7, s2
			subss		xmm4, xmm5

			mulss		xmm0, xmm7
			movss		n0, xmm0
			mulss		xmm2, xmm7
			movss		n1, xmm2
			mulss		xmm4, xmm7
			movss		n2, xmm4

			movss		xmm0, d0
			mulss		xmm0, d9
			movss		xmm1, d4
			mulss		xmm1, d5

			movss		xmm2, d1
			mulss		xmm2, d9
			movss		xmm3, d4
			mulss		xmm3, d6

			movss		xmm4, d2
			mulss		xmm4, d9
			movss		xmm5, d4
			mulss		xmm5, d7

			subss		xmm0, xmm1
			subss		xmm2, xmm3
			movss		xmm7, s0
			subss		xmm4, xmm5

			mulss		xmm0, xmm7
			movss		t0, xmm0
			mulss		xmm2, xmm7
			movss		t1, xmm2
			mulss		xmm4, xmm7
			movss		t2, xmm4

#ifndef DERIVE_UNSMOOTHED_BITANGENT
			movss		xmm0, d3
			mulss		xmm0, d5
			movss		xmm1, d0
			mulss		xmm1, d8

			movss		xmm2, d3
			mulss		xmm2, d6
			movss		xmm3, d1
			mulss		xmm3, d8

			movss		xmm4, d3
			mulss		xmm4, d7
			movss		xmm5, d2
			mulss		xmm5, d8
#else
			movss		xmm0, n2
			mulss		xmm0, t1
			movss		xmm1, n1
			mulss		xmm1, t2

			movss		xmm2, n0
			mulss		xmm2, t2
			movss		xmm3, n2
			mulss		xmm3, t0

			movss		xmm4, n1
			mulss		xmm4, t0
			movss		xmm5, n0
			mulss		xmm5, t1
#endif
			subss		xmm0, xmm1
			subss		xmm2, xmm3
			movss		xmm7, s1
			subss		xmm4, xmm5

			mulss		xmm0, xmm7
			movss		t3, xmm0
			mulss		xmm2, xmm7
			movss		t4, xmm2
			mulss		xmm4, xmm7
			movss		t5, xmm4
		}

#else

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
CSIMD_SSE:: normalizeTangents
============
*/
void VPCALL CSIMD_SSE:: normalizeTangents( CVertex *verts, const int numVerts ) {
	ALIGN16( float normal[12] );

	SMF_ASSERT( sizeof( CVertex ) == DRAWVERT_SIZE );
	SMF_ASSERT( (int)&((CVertex *)0)->normal == DRAWVERT_NORMAL_OFFSET );
	SMF_ASSERT( (int)&((CVertex *)0)->tangents[0] == DRAWVERT_TANGENT0_OFFSET );
	SMF_ASSERT( (int)&((CVertex *)0)->tangents[1] == DRAWVERT_TANGENT1_OFFSET );

	SMF_ASSERT( verts != NULL );
	SMF_ASSERT( numVerts >= 0 );

	__asm {
		mov			eax, numVerts
		test		eax, eax
		jz			done
#ifdef REFINE_TANGENT_SQUAREROOT
		movaps		xmm6, SIMD_SP_rsqrt_c0
		movaps		xmm7, SIMD_SP_rsqrt_c1
#endif
		mov			esi, verts
		imul		eax, DRAWVERT_SIZE
		add			esi, eax
		neg			eax
		add			eax, DRAWVERT_SIZE*4
		jle			loopVert4

		sub			eax, DRAWVERT_SIZE*4
		jl			loopVert1

	loopVert4:

		sub			eax, DRAWVERT_SIZE*4

		// normalize 4 CVertex::normal

		movss		xmm0, [esi+eax+DRAWVERT_SIZE*0+DRAWVERT_NORMAL_OFFSET+0]	//  0,  X,  X,  X
		movhps		xmm0, [esi+eax+DRAWVERT_SIZE*1+DRAWVERT_NORMAL_OFFSET+0]	//  0,  X,  3,  4
		movss		xmm2, [esi+eax+DRAWVERT_SIZE*1+DRAWVERT_NORMAL_OFFSET+8]	//  5,  X,  X,  X
		movhps		xmm2, [esi+eax+DRAWVERT_SIZE*0+DRAWVERT_NORMAL_OFFSET+4]	//	5,  X,  1,  2
		movss		xmm4, [esi+eax+DRAWVERT_SIZE*2+DRAWVERT_NORMAL_OFFSET+0]	//  6,  X,  X,  X
		movhps		xmm4, [esi+eax+DRAWVERT_SIZE*3+DRAWVERT_NORMAL_OFFSET+0]	//  6,  X,  9, 10
		movss		xmm3, [esi+eax+DRAWVERT_SIZE*3+DRAWVERT_NORMAL_OFFSET+8]	// 11,  X,  X,  X
		movhps		xmm3, [esi+eax+DRAWVERT_SIZE*2+DRAWVERT_NORMAL_OFFSET+4]	// 11,  X,  7,  8

		movaps		xmm1, xmm0
		movaps		xmm5, xmm2
		shufps		xmm0, xmm4, R_SHUFFLEPS( 0, 2, 0, 2 )		//  0,  3,  6,  9
		shufps		xmm2, xmm3, R_SHUFFLEPS( 3, 0, 3, 0 )		//  2,  5,  8, 11
		shufps		xmm1, xmm5, R_SHUFFLEPS( 3, 3, 2, 2 )		//  4,  4,  1,  1
		shufps		xmm4, xmm3, R_SHUFFLEPS( 3, 3, 2, 2 )		// 10, 10,  7,  7
		shufps		xmm1, xmm4, R_SHUFFLEPS( 2, 0, 2, 0 )		//  1,  4,  7, 10

		movaps		xmm3, xmm0
		movaps		xmm4, xmm1
		movaps		xmm5, xmm2

		mulps		xmm3, xmm3
		mulps		xmm4, xmm4
		mulps		xmm5, xmm5
		addps		xmm3, xmm4
		addps		xmm3, xmm5

#ifdef REFINE_TANGENT_SQUAREROOT
		rsqrtps		xmm4, xmm3
		mulps		xmm3, xmm4
		mulps		xmm3, xmm4
		subps		xmm3, xmm6
		mulps		xmm4, xmm7
		mulps		xmm3, xmm4
#else
		rsqrtps		xmm3, xmm3
#endif

		mulps		xmm0, xmm3
		mulps		xmm1, xmm3
		mulps		xmm2, xmm3

		// save the 4 CVertex::normal to project the tangents

		movaps		[normal+ 0], xmm0
		movaps		[normal+16], xmm1
		movaps		[normal+32], xmm2

		movss		[esi+eax+DRAWVERT_SIZE*0+DRAWVERT_NORMAL_OFFSET+0], xmm0
		movss		[esi+eax+DRAWVERT_SIZE*0+DRAWVERT_NORMAL_OFFSET+4], xmm1
		movss		[esi+eax+DRAWVERT_SIZE*0+DRAWVERT_NORMAL_OFFSET+8], xmm2

		shufps		xmm0, xmm0, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm2, xmm2, R_SHUFFLEPS( 1, 2, 3, 0 )

		movss		[esi+eax+DRAWVERT_SIZE*1+DRAWVERT_NORMAL_OFFSET+0], xmm0
		movss		[esi+eax+DRAWVERT_SIZE*1+DRAWVERT_NORMAL_OFFSET+4], xmm1
		movss		[esi+eax+DRAWVERT_SIZE*1+DRAWVERT_NORMAL_OFFSET+8], xmm2

		shufps		xmm0, xmm0, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm2, xmm2, R_SHUFFLEPS( 1, 2, 3, 0 )

		movss		[esi+eax+DRAWVERT_SIZE*2+DRAWVERT_NORMAL_OFFSET+0], xmm0
		movss		[esi+eax+DRAWVERT_SIZE*2+DRAWVERT_NORMAL_OFFSET+4], xmm1
		movss		[esi+eax+DRAWVERT_SIZE*2+DRAWVERT_NORMAL_OFFSET+8], xmm2

		shufps		xmm0, xmm0, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm2, xmm2, R_SHUFFLEPS( 1, 2, 3, 0 )

		movss		[esi+eax+DRAWVERT_SIZE*3+DRAWVERT_NORMAL_OFFSET+0], xmm0
		movss		[esi+eax+DRAWVERT_SIZE*3+DRAWVERT_NORMAL_OFFSET+4], xmm1
		movss		[esi+eax+DRAWVERT_SIZE*3+DRAWVERT_NORMAL_OFFSET+8], xmm2

		// project and normalize 4 CVertex::tangent[0]

		movss		xmm0, [esi+eax+DRAWVERT_SIZE*0+DRAWVERT_TANGENT0_OFFSET+0]	//  0,  X,  X,  X
		movhps		xmm0, [esi+eax+DRAWVERT_SIZE*1+DRAWVERT_TANGENT0_OFFSET+0]	//  0,  X,  3,  4
		movss		xmm2, [esi+eax+DRAWVERT_SIZE*1+DRAWVERT_TANGENT0_OFFSET+8]	//  5,  X,  X,  X
		movhps		xmm2, [esi+eax+DRAWVERT_SIZE*0+DRAWVERT_TANGENT0_OFFSET+4]	//	5,  X,  1,  2
		movss		xmm4, [esi+eax+DRAWVERT_SIZE*2+DRAWVERT_TANGENT0_OFFSET+0]	//  6,  X,  X,  X
		movhps		xmm4, [esi+eax+DRAWVERT_SIZE*3+DRAWVERT_TANGENT0_OFFSET+0]	//  6,  X,  9, 10
		movss		xmm3, [esi+eax+DRAWVERT_SIZE*3+DRAWVERT_TANGENT0_OFFSET+8]	// 11,  X,  X,  X
		movhps		xmm3, [esi+eax+DRAWVERT_SIZE*2+DRAWVERT_TANGENT0_OFFSET+4]	// 11,  X,  7,  8

		movaps		xmm1, xmm0
		movaps		xmm5, xmm2
		shufps		xmm0, xmm4, R_SHUFFLEPS( 0, 2, 0, 2 )		//  0,  3,  6,  9
		shufps		xmm2, xmm3, R_SHUFFLEPS( 3, 0, 3, 0 )		//  2,  5,  8, 11
		shufps		xmm1, xmm5, R_SHUFFLEPS( 3, 3, 2, 2 )		//  4,  4,  1,  1
		shufps		xmm4, xmm3, R_SHUFFLEPS( 3, 3, 2, 2 )		// 10, 10,  7,  7
		shufps		xmm1, xmm4, R_SHUFFLEPS( 2, 0, 2, 0 )		//  1,  4,  7, 10

		movaps		xmm3, xmm0
		movaps		xmm4, xmm1
		movaps		xmm5, xmm2

		mulps		xmm3, [normal+ 0]
		mulps		xmm4, [normal+16]
		mulps		xmm5, [normal+32]
		addps		xmm3, xmm4
		addps		xmm3, xmm5

		movaps		xmm4, xmm3
		movaps		xmm5, xmm3
		mulps		xmm3, [normal+ 0]
		mulps		xmm4, [normal+16]
		mulps		xmm5, [normal+32]
		subps		xmm0, xmm3
		subps		xmm1, xmm4
		subps		xmm2, xmm5

		movaps		xmm3, xmm0
		movaps		xmm4, xmm1
		movaps		xmm5, xmm2

		mulps		xmm3, xmm3
		mulps		xmm4, xmm4
		mulps		xmm5, xmm5
		addps		xmm3, xmm4
		addps		xmm3, xmm5

#ifdef REFINE_TANGENT_SQUAREROOT
		rsqrtps		xmm4, xmm3
		mulps		xmm3, xmm4
		mulps		xmm3, xmm4
		subps		xmm3, xmm6
		mulps		xmm4, xmm7
		mulps		xmm3, xmm4
#else
		rsqrtps		xmm3, xmm3
#endif

		mulps		xmm0, xmm3
		mulps		xmm1, xmm3
		mulps		xmm2, xmm3

		movss		[esi+eax+DRAWVERT_SIZE*0+DRAWVERT_TANGENT0_OFFSET+0], xmm0
		movss		[esi+eax+DRAWVERT_SIZE*0+DRAWVERT_TANGENT0_OFFSET+4], xmm1
		movss		[esi+eax+DRAWVERT_SIZE*0+DRAWVERT_TANGENT0_OFFSET+8], xmm2

		shufps		xmm0, xmm0, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm2, xmm2, R_SHUFFLEPS( 1, 2, 3, 0 )

		movss		[esi+eax+DRAWVERT_SIZE*1+DRAWVERT_TANGENT0_OFFSET+0], xmm0
		movss		[esi+eax+DRAWVERT_SIZE*1+DRAWVERT_TANGENT0_OFFSET+4], xmm1
		movss		[esi+eax+DRAWVERT_SIZE*1+DRAWVERT_TANGENT0_OFFSET+8], xmm2

		shufps		xmm0, xmm0, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm2, xmm2, R_SHUFFLEPS( 1, 2, 3, 0 )

		movss		[esi+eax+DRAWVERT_SIZE*2+DRAWVERT_TANGENT0_OFFSET+0], xmm0
		movss		[esi+eax+DRAWVERT_SIZE*2+DRAWVERT_TANGENT0_OFFSET+4], xmm1
		movss		[esi+eax+DRAWVERT_SIZE*2+DRAWVERT_TANGENT0_OFFSET+8], xmm2

		shufps		xmm0, xmm0, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm2, xmm2, R_SHUFFLEPS( 1, 2, 3, 0 )

		movss		[esi+eax+DRAWVERT_SIZE*3+DRAWVERT_TANGENT0_OFFSET+0], xmm0
		movss		[esi+eax+DRAWVERT_SIZE*3+DRAWVERT_TANGENT0_OFFSET+4], xmm1
		movss		[esi+eax+DRAWVERT_SIZE*3+DRAWVERT_TANGENT0_OFFSET+8], xmm2

		// project and normalize 4 CVertex::tangent[1]

		movss		xmm0, [esi+eax+DRAWVERT_SIZE*0+DRAWVERT_TANGENT1_OFFSET+0]	//  0,  X,  X,  X
		movhps		xmm0, [esi+eax+DRAWVERT_SIZE*1+DRAWVERT_TANGENT1_OFFSET+0]	//  0,  X,  3,  4
		movss		xmm2, [esi+eax+DRAWVERT_SIZE*1+DRAWVERT_TANGENT1_OFFSET+8]	//  5,  X,  X,  X
		movhps		xmm2, [esi+eax+DRAWVERT_SIZE*0+DRAWVERT_TANGENT1_OFFSET+4]	//	5,  X,  1,  2
		movss		xmm4, [esi+eax+DRAWVERT_SIZE*2+DRAWVERT_TANGENT1_OFFSET+0]	//  6,  X,  X,  X
		movhps		xmm4, [esi+eax+DRAWVERT_SIZE*3+DRAWVERT_TANGENT1_OFFSET+0]	//  6,  X,  9, 10
		movss		xmm3, [esi+eax+DRAWVERT_SIZE*3+DRAWVERT_TANGENT1_OFFSET+8]	// 11,  X,  X,  X
		movhps		xmm3, [esi+eax+DRAWVERT_SIZE*2+DRAWVERT_TANGENT1_OFFSET+4]	// 11,  X,  7,  8

		movaps		xmm1, xmm0
		movaps		xmm5, xmm2
		shufps		xmm0, xmm4, R_SHUFFLEPS( 0, 2, 0, 2 )		//  0,  3,  6,  9
		shufps		xmm2, xmm3, R_SHUFFLEPS( 3, 0, 3, 0 )		//  2,  5,  8, 11
		shufps		xmm1, xmm5, R_SHUFFLEPS( 3, 3, 2, 2 )		//  4,  4,  1,  1
		shufps		xmm4, xmm3, R_SHUFFLEPS( 3, 3, 2, 2 )		// 10, 10,  7,  7
		shufps		xmm1, xmm4, R_SHUFFLEPS( 2, 0, 2, 0 )		//  1,  4,  7, 10

		movaps		xmm3, xmm0
		movaps		xmm4, xmm1
		movaps		xmm5, xmm2

		mulps		xmm3, [normal+ 0]
		mulps		xmm4, [normal+16]
		mulps		xmm5, [normal+32]
		addps		xmm3, xmm4
		addps		xmm3, xmm5

		movaps		xmm4, xmm3
		movaps		xmm5, xmm3
		mulps		xmm3, [normal+ 0]
		mulps		xmm4, [normal+16]
		mulps		xmm5, [normal+32]
		subps		xmm0, xmm3
		subps		xmm1, xmm4
		subps		xmm2, xmm5

		movaps		xmm3, xmm0
		movaps		xmm4, xmm1
		movaps		xmm5, xmm2

		mulps		xmm3, xmm3
		mulps		xmm4, xmm4
		mulps		xmm5, xmm5
		addps		xmm3, xmm4
		addps		xmm3, xmm5

#ifdef REFINE_TANGENT_SQUAREROOT
		rsqrtps		xmm4, xmm3
		mulps		xmm3, xmm4
		mulps		xmm3, xmm4
		subps		xmm3, xmm6
		mulps		xmm4, xmm7
		mulps		xmm3, xmm4
#else
		rsqrtps		xmm3, xmm3
#endif

		mulps		xmm0, xmm3
		mulps		xmm1, xmm3
		mulps		xmm2, xmm3

		movss		[esi+eax+DRAWVERT_SIZE*0+DRAWVERT_TANGENT1_OFFSET+0], xmm0
		movss		[esi+eax+DRAWVERT_SIZE*0+DRAWVERT_TANGENT1_OFFSET+4], xmm1
		movss		[esi+eax+DRAWVERT_SIZE*0+DRAWVERT_TANGENT1_OFFSET+8], xmm2

		shufps		xmm0, xmm0, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm2, xmm2, R_SHUFFLEPS( 1, 2, 3, 0 )

		movss		[esi+eax+DRAWVERT_SIZE*1+DRAWVERT_TANGENT1_OFFSET+0], xmm0
		movss		[esi+eax+DRAWVERT_SIZE*1+DRAWVERT_TANGENT1_OFFSET+4], xmm1
		movss		[esi+eax+DRAWVERT_SIZE*1+DRAWVERT_TANGENT1_OFFSET+8], xmm2

		shufps		xmm0, xmm0, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm2, xmm2, R_SHUFFLEPS( 1, 2, 3, 0 )

		movss		[esi+eax+DRAWVERT_SIZE*2+DRAWVERT_TANGENT1_OFFSET+0], xmm0
		movss		[esi+eax+DRAWVERT_SIZE*2+DRAWVERT_TANGENT1_OFFSET+4], xmm1
		movss		[esi+eax+DRAWVERT_SIZE*2+DRAWVERT_TANGENT1_OFFSET+8], xmm2

		shufps		xmm0, xmm0, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 2, 3, 0 )
		shufps		xmm2, xmm2, R_SHUFFLEPS( 1, 2, 3, 0 )

		movss		[esi+eax+DRAWVERT_SIZE*3+DRAWVERT_TANGENT1_OFFSET+0], xmm0
		movss		[esi+eax+DRAWVERT_SIZE*3+DRAWVERT_TANGENT1_OFFSET+4], xmm1
		movss		[esi+eax+DRAWVERT_SIZE*3+DRAWVERT_TANGENT1_OFFSET+8], xmm2

		add			eax, DRAWVERT_SIZE*8

		jle			loopVert4

		sub			eax, DRAWVERT_SIZE*4
		jge			done

	loopVert1:

		// normalize one CVertex::normal

		movss		xmm0, [esi+eax+DRAWVERT_NORMAL_OFFSET+0]
		movss		xmm1, [esi+eax+DRAWVERT_NORMAL_OFFSET+4]
		movss		xmm2, [esi+eax+DRAWVERT_NORMAL_OFFSET+8]
		movss		xmm3, xmm0
		movss		xmm4, xmm1
		movss		xmm5, xmm2

		mulss		xmm3, xmm3
		mulss		xmm4, xmm4
		mulss		xmm5, xmm5
		addss		xmm3, xmm4
		addss		xmm3, xmm5

#ifdef REFINE_TANGENT_SQUAREROOT
		rsqrtss		xmm4, xmm3
		mulss		xmm3, xmm4
		mulss		xmm3, xmm4
		subss		xmm3, xmm6
		mulss		xmm4, xmm7
		mulss		xmm3, xmm4
#else
		rsqrtss		xmm3, xmm3
#endif

		mulss		xmm0, xmm3
		mulss		xmm1, xmm3
		mulss		xmm2, xmm3

		movss		[esi+eax+DRAWVERT_NORMAL_OFFSET+0], xmm0
		movss		[esi+eax+DRAWVERT_NORMAL_OFFSET+4], xmm1
		movss		[esi+eax+DRAWVERT_NORMAL_OFFSET+8], xmm2

		// project and normalize one CVertex::tangent[0]

		movss		xmm0, [esi+eax+DRAWVERT_TANGENT0_OFFSET+0]
		movss		xmm1, [esi+eax+DRAWVERT_TANGENT0_OFFSET+4]
		movss		xmm2, [esi+eax+DRAWVERT_TANGENT0_OFFSET+8]
		movss		xmm3, xmm0
		movss		xmm4, xmm1
		movss		xmm5, xmm2

		mulss		xmm3, [esi+eax+DRAWVERT_NORMAL_OFFSET+0]
		mulss		xmm4, [esi+eax+DRAWVERT_NORMAL_OFFSET+4]
		mulss		xmm5, [esi+eax+DRAWVERT_NORMAL_OFFSET+8]
		addss		xmm3, xmm4
		addss		xmm3, xmm5

		movss		xmm4, xmm3
		movss		xmm5, xmm3
		mulss		xmm3, [esi+eax+DRAWVERT_NORMAL_OFFSET+0]
		mulss		xmm4, [esi+eax+DRAWVERT_NORMAL_OFFSET+4]
		mulss		xmm5, [esi+eax+DRAWVERT_NORMAL_OFFSET+8]
		subss		xmm0, xmm3
		subss		xmm1, xmm4
		subss		xmm2, xmm5

		movss		xmm3, xmm0
		movss		xmm4, xmm1
		movss		xmm5, xmm2

		mulss		xmm3, xmm3
		mulss		xmm4, xmm4
		mulss		xmm5, xmm5
		addss		xmm3, xmm4
		addss		xmm3, xmm5

#ifdef REFINE_TANGENT_SQUAREROOT
		rsqrtss		xmm4, xmm3
		mulss		xmm3, xmm4
		mulss		xmm3, xmm4
		subss		xmm3, xmm6
		mulss		xmm4, xmm7
		mulss		xmm3, xmm4
#else
		rsqrtss		xmm3, xmm3
#endif

		mulss		xmm0, xmm3
		mulss		xmm1, xmm3
		mulss		xmm2, xmm3

		movss		[esi+eax+DRAWVERT_TANGENT0_OFFSET+0], xmm0
		movss		[esi+eax+DRAWVERT_TANGENT0_OFFSET+4], xmm1
		movss		[esi+eax+DRAWVERT_TANGENT0_OFFSET+8], xmm2

		// project and normalize one CVertex::tangent[1]

		movss		xmm0, [esi+eax+DRAWVERT_TANGENT1_OFFSET+0]
		movss		xmm1, [esi+eax+DRAWVERT_TANGENT1_OFFSET+4]
		movss		xmm2, [esi+eax+DRAWVERT_TANGENT1_OFFSET+8]
		movss		xmm3, xmm0
		movss		xmm4, xmm1
		movss		xmm5, xmm2

		mulss		xmm3, [esi+eax+DRAWVERT_NORMAL_OFFSET+0]
		mulss		xmm4, [esi+eax+DRAWVERT_NORMAL_OFFSET+4]
		mulss		xmm5, [esi+eax+DRAWVERT_NORMAL_OFFSET+8]
		addss		xmm3, xmm4
		addss		xmm3, xmm5

		movss		xmm4, xmm3
		movss		xmm5, xmm3
		mulss		xmm3, [esi+eax+DRAWVERT_NORMAL_OFFSET+0]
		mulss		xmm4, [esi+eax+DRAWVERT_NORMAL_OFFSET+4]
		mulss		xmm5, [esi+eax+DRAWVERT_NORMAL_OFFSET+8]
		subss		xmm0, xmm3
		subss		xmm1, xmm4
		subss		xmm2, xmm5

		movss		xmm3, xmm0
		movss		xmm4, xmm1
		movss		xmm5, xmm2

		mulss		xmm3, xmm3
		mulss		xmm4, xmm4
		mulss		xmm5, xmm5
		addss		xmm3, xmm4
		addss		xmm3, xmm5

#ifdef REFINE_TANGENT_SQUAREROOT
		rsqrtss		xmm4, xmm3
		mulss		xmm3, xmm4
		mulss		xmm3, xmm4
		subss		xmm3, xmm6
		mulss		xmm4, xmm7
		mulss		xmm3, xmm4
#else
		rsqrtss		xmm3, xmm3
#endif

		mulss		xmm0, xmm3
		mulss		xmm1, xmm3
		mulss		xmm2, xmm3

		movss		[esi+eax+DRAWVERT_TANGENT1_OFFSET+0], xmm0
		movss		[esi+eax+DRAWVERT_TANGENT1_OFFSET+4], xmm1
		movss		[esi+eax+DRAWVERT_TANGENT1_OFFSET+8], xmm2

		add			eax, DRAWVERT_SIZE

		jl			loopVert1
	done:
	}
}

/*
============
CSIMD_SSE:: createTextureSpaceLightVectors
============
*/
void VPCALL CSIMD_SSE:: createTextureSpaceLightVectors( CVec3D *lightVectors, const CVec3D &lightOrigin, const CVertex *verts, const int numVerts, const int *indexes, const int numIndexes ) {

	SMF_ASSERT( sizeof( CVertex ) == DRAWVERT_SIZE );
	SMF_ASSERT( (int)&((CVertex *)0)->xyz == DRAWVERT_XYZ_OFFSET );
	SMF_ASSERT( (int)&((CVertex *)0)->normal == DRAWVERT_NORMAL_OFFSET );
	SMF_ASSERT( (int)&((CVertex *)0)->tangents[0] == DRAWVERT_TANGENT0_OFFSET );
	SMF_ASSERT( (int)&((CVertex *)0)->tangents[1] == DRAWVERT_TANGENT1_OFFSET );

	bool *used = (bool *)_alloca16( numVerts * sizeof( used[0] ) );
	memset( used, 0, numVerts * sizeof( used[0] ) );

	for ( int i = numIndexes - 1; i >= 0; i-- ) {
		used[indexes[i]] = true;
	}

#if 0

	__asm {

		mov			eax, numVerts

		mov			esi, used
		add			esi, eax

		mov			edi, verts
		sub			edi, DRAWVERT_SIZE

		neg			eax
		dec			eax

		mov			ecx, lightOrigin
		movss		xmm7, [ecx+0]
		movhps		xmm7, [ecx+4]

		mov			ecx, lightVectors
		sub			ecx, 3*4

	loopVert:
		inc			eax
		jge			done

		add			edi, DRAWVERT_SIZE
		add			ecx, 3*4

		cmp			sf_u8 ptr [esi+eax], 0
		je			loopVert

		movaps		xmm0, xmm7
		movss		xmm1, [edi+DRAWVERT_XYZ_OFFSET+0]
		movhps		xmm1, [edi+DRAWVERT_XYZ_OFFSET+4]
		subps		xmm0, xmm1

		// 0,  X,  1,  2
		// 3,  X,  4,  5
		// 6,  X,  7,  8

		movss		xmm2, [edi+DRAWVERT_TANGENT0_OFFSET+0]
		movhps		xmm2, [edi+DRAWVERT_TANGENT0_OFFSET+4]
		mulps		xmm2, xmm0

		movss		xmm3, [edi+DRAWVERT_TANGENT1_OFFSET+0]
		movhps		xmm3, [edi+DRAWVERT_TANGENT1_OFFSET+4]
		mulps		xmm3, xmm0

		movaps		xmm5, xmm2								// xmm5 = 0,  X,  1,  2
		unpcklps	xmm5, xmm3								// xmm5 = 0,  3,  X,  X
		unpckhps	xmm2, xmm3								// xmm2 = 1,  4,  2,  5

		movss		xmm4, [edi+DRAWVERT_NORMAL_OFFSET+0]
		movhps		xmm4, [edi+DRAWVERT_NORMAL_OFFSET+4]
		mulps		xmm4, xmm0

		movlhps		xmm5, xmm4								// xmm5 = 0,  3,  6,  X
		movhlps		xmm4, xmm2								// xmm4 = 2,  5,  7,  8
		shufps		xmm2, xmm4, R_SHUFFLEPS( 0, 1, 3, 2 )	// xmm2 = 2,  5,  8,  7

		addps		xmm5, xmm4
		addps		xmm5, xmm2
		movlps		[ecx+0], xmm5
		shufps		xmm5, xmm5, R_SHUFFLEPS( 2, 3, 0, 1 )
		movss		[ecx+8], xmm5

		jmp			loopVert

	done:
	}

#elif 1

	for ( int i = 0; i < numVerts; i++ ) {
		if ( !used[i] ) {
			continue;
		}

		const CVertex *v = &verts[i];
		CVec3D lightDir;

		lightDir[0] = lightOrigin[0] - v->xyz[0];
		lightDir[1] = lightOrigin[1] - v->xyz[1];
		lightDir[2] = lightOrigin[2] - v->xyz[2];

		lightVectors[i][0] = lightDir[0] * v->tangents[0][0] + lightDir[1] * v->tangents[0][1] + lightDir[2] * v->tangents[0][2];
		lightVectors[i][1] = lightDir[0] * v->tangents[1][0] + lightDir[1] * v->tangents[1][1] + lightDir[2] * v->tangents[1][2];
		lightVectors[i][2] = lightDir[0] * v->normal[0] + lightDir[1] * v->normal[1] + lightDir[2] * v->normal[2];
	}

#elif 1

	ALIGN16( int usedVertNums[4] );
	ALIGN16( float lightDir0[4] );
	ALIGN16( float lightDir1[4] );
	ALIGN16( float lightDir2[4] );
	ALIGN16( float normal0[4] );
	ALIGN16( float normal1[4] );
	ALIGN16( float normal2[4] );
	ALIGN16( float tangent0[4] );
	ALIGN16( float tangent1[4] );
	ALIGN16( float tangent2[4] );
	ALIGN16( float tangent3[4] );
	ALIGN16( float tangent4[4] );
	ALIGN16( float tangent5[4] );
	CVec3D localLightOrigin = lightOrigin;

	__asm {

		xor			ecx, ecx
		mov			eax, numVerts

		mov			esi, used
		add			esi, eax

		mov			edi, verts
		sub			edi, DRAWVERT_SIZE

		neg			eax
		dec			eax

	loopVert4:
		inc			eax
		jge			done4

		add			edi, DRAWVERT_SIZE

		cmp			byte ptr [esi+eax], 0
		je			loopVert4

		mov			usedVertNums[ecx*4], eax

		inc			ecx
		cmp			ecx, 4

		movss		xmm0, localLightOrigin[0]
		movss		xmm1, localLightOrigin[4]
		movss		xmm2, localLightOrigin[8]

		subss		xmm0, [edi+DRAWVERT_XYZ_OFFSET+0]
		subss		xmm1, [edi+DRAWVERT_XYZ_OFFSET+4]
		subss		xmm2, [edi+DRAWVERT_XYZ_OFFSET+8]

		movss		lightDir0[ecx*4-4], xmm0
		movss		lightDir1[ecx*4-4], xmm1
		movss		lightDir2[ecx*4-4], xmm2

		movss		xmm3, [edi+DRAWVERT_NORMAL_OFFSET+0]
		movss		xmm4, [edi+DRAWVERT_NORMAL_OFFSET+4]
		movss		xmm5, [edi+DRAWVERT_NORMAL_OFFSET+8]

		movss		normal0[ecx*4-4], xmm3
		movss		normal1[ecx*4-4], xmm4
		movss		normal2[ecx*4-4], xmm5

		movss		xmm0, [edi+DRAWVERT_TANGENT0_OFFSET+0]
		movss		xmm1, [edi+DRAWVERT_TANGENT0_OFFSET+4]
		movss		xmm2, [edi+DRAWVERT_TANGENT0_OFFSET+8]

		movss		tangent0[ecx*4-4], xmm0
		movss		tangent1[ecx*4-4], xmm1
		movss		tangent2[ecx*4-4], xmm2

		movss		xmm3, [edi+DRAWVERT_TANGENT1_OFFSET+0]
		movss		xmm4, [edi+DRAWVERT_TANGENT1_OFFSET+4]
		movss		xmm5, [edi+DRAWVERT_TANGENT1_OFFSET+8]

		movss		tangent3[ecx*4-4], xmm3
		movss		tangent4[ecx*4-4], xmm4
		movss		tangent5[ecx*4-4], xmm5

		jl			loopVert4

		movaps		xmm0, lightDir0
		movaps		xmm1, lightDir1
		movaps		xmm2, lightDir2

		movaps		xmm3, tangent0
		mulps		xmm3, xmm0
		movaps		xmm4, tangent1
		mulps		xmm4, xmm1
		movaps		xmm5, tangent2
		mulps		xmm5, xmm2

		addps		xmm3, xmm4
		addps		xmm5, xmm3

		movaps		xmm3, tangent3
		mulps		xmm3, xmm0
		movaps		xmm4, tangent4
		mulps		xmm4, xmm1
		movaps		xmm6, tangent5
		mulps		xmm6, xmm2

		addps		xmm3, xmm4
		addps		xmm6, xmm3

		mulps		xmm0, normal0
		mulps		xmm1, normal1
		mulps		xmm2, normal2

		addps		xmm0, xmm1
		addps		xmm0, xmm2

		mov			ecx, numVerts
		imul		ecx, 12
		mov			edx, usedVertNums[0]
		add			ecx, lightVectors
		imul		edx, 12

		movss		[ecx+edx+0], xmm5
		movss		[ecx+edx+4], xmm6
		movss		[ecx+edx+8], xmm0

		shufps		xmm5, xmm5, R_SHUFFLEPS( 1, 2, 3, 0 )
		mov			edx, usedVertNums[4]
		shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 2, 3, 0 )
		imul		edx, 12
		shufps		xmm0, xmm0, R_SHUFFLEPS( 1, 2, 3, 0 )

		movss		[ecx+edx+0], xmm5
		movss		[ecx+edx+4], xmm6
		movss		[ecx+edx+8], xmm0

		shufps		xmm5, xmm5, R_SHUFFLEPS( 1, 2, 3, 0 )
		mov			edx, usedVertNums[8]
		shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 2, 3, 0 )
		imul		edx, 12
		shufps		xmm0, xmm0, R_SHUFFLEPS( 1, 2, 3, 0 )

		movss		[ecx+edx+0], xmm5
		movss		[ecx+edx+4], xmm6
		movss		[ecx+edx+8], xmm0

		shufps		xmm5, xmm5, R_SHUFFLEPS( 1, 2, 3, 0 )
		mov			edx, usedVertNums[12]
		shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 2, 3, 0 )
		imul		edx, 12
		shufps		xmm0, xmm0, R_SHUFFLEPS( 1, 2, 3, 0 )

		movss		[ecx+edx+0], xmm5
		movss		[ecx+edx+4], xmm6
		movss		[ecx+edx+8], xmm0

		xor			ecx, ecx
		jmp			loopVert4

	done4:
		test		ecx, ecx
		jz			done
		xor			eax, eax
		mov			edi, numVerts
		imul		edi, 12
		add			edi, lightVectors

	loopVert1:
		movss		xmm0, lightDir0[eax*4]
		movss		xmm1, lightDir1[eax*4]
		movss		xmm2, lightDir2[eax*4]

		mov			edx, usedVertNums[eax*4]
		imul		edx, 12

		movss		xmm3, tangent0[eax*4]
		mulss		xmm3, xmm0
		movss		xmm4, tangent1[eax*4]
		mulss		xmm4, xmm1
		movss		xmm5, tangent2[eax*4]
		mulss		xmm5, xmm2

		addss		xmm3, xmm4
		addss		xmm5, xmm3
		movss		[edi+edx+0], xmm5

		movss		xmm3, tangent3[eax*4]
		mulss		xmm3, xmm0
		movss		xmm4, tangent4[eax*4]
		mulss		xmm4, xmm1
		movss		xmm6, tangent5[eax*4]
		mulss		xmm6, xmm2

		addss		xmm3, xmm4
		addss		xmm6, xmm3
		movss		[edi+edx+4], xmm6

		mulss		xmm0, normal0[eax*4]
		mulss		xmm1, normal1[eax*4]
		mulss		xmm2, normal2[eax*4]

		addss		xmm0, xmm1
		addss		xmm0, xmm2
		movss		[edi+edx+8], xmm0

		inc			eax
		dec			ecx
		jg			loopVert1

	done:
	}

#else

	ALIGN16( float lightVectors0[4] );
	ALIGN16( float lightVectors1[4] );
	ALIGN16( float lightVectors2[4] );
	int numUsedVerts = 0;

	for ( int i = 0; i < numVerts; i++ ) {
		if ( !used[i] ) {
			continue;
		}

		const CVertex *v = &verts[i];

		lightDir0[numUsedVerts] = lightOrigin[0] - v->xyz[0];
		lightDir1[numUsedVerts] = lightOrigin[1] - v->xyz[1];
		lightDir2[numUsedVerts] = lightOrigin[2] - v->xyz[2];

		normal0[numUsedVerts] = v->normal[0];
		normal1[numUsedVerts] = v->normal[1];
		normal2[numUsedVerts] = v->normal[2];

		tangent0[numUsedVerts] = v->tangents[0][0];
		tangent1[numUsedVerts] = v->tangents[0][1];
		tangent2[numUsedVerts] = v->tangents[0][2];

		tangent3[numUsedVerts] = v->tangents[1][0];
		tangent4[numUsedVerts] = v->tangents[1][1];
		tangent5[numUsedVerts] = v->tangents[1][2];

		usedVertNums[numUsedVerts++] = i;
		if ( numUsedVerts < 4 ) {
			continue;
		}

		lightVectors0[0] = lightDir0[0] * tangent0[0];
		lightVectors0[1] = lightDir0[1] * tangent0[1];
		lightVectors0[2] = lightDir0[2] * tangent0[2];
		lightVectors0[3] = lightDir0[3] * tangent0[3];

		lightVectors0[0] += lightDir1[0] * tangent1[0];
		lightVectors0[1] += lightDir1[1] * tangent1[1];
		lightVectors0[2] += lightDir1[2] * tangent1[2];
		lightVectors0[3] += lightDir1[3] * tangent1[3];

		lightVectors0[0] += lightDir2[0] * tangent2[0];
		lightVectors0[1] += lightDir2[1] * tangent2[1];
		lightVectors0[2] += lightDir2[2] * tangent2[2];
		lightVectors0[3] += lightDir2[3] * tangent2[3];

		lightVectors1[0] = lightDir0[0] * tangent3[0];
		lightVectors1[1] = lightDir0[1] * tangent3[1];
		lightVectors1[2] = lightDir0[2] * tangent3[2];
		lightVectors1[3] = lightDir0[3] * tangent3[3];

		lightVectors1[0] += lightDir1[0] * tangent4[0];
		lightVectors1[1] += lightDir1[1] * tangent4[1];
		lightVectors1[2] += lightDir1[2] * tangent4[2];
		lightVectors1[3] += lightDir1[3] * tangent4[3];

		lightVectors1[0] += lightDir2[0] * tangent5[0];
		lightVectors1[1] += lightDir2[1] * tangent5[1];
		lightVectors1[2] += lightDir2[2] * tangent5[2];
		lightVectors1[3] += lightDir2[3] * tangent5[3];

		lightVectors2[0] = lightDir0[0] * normal0[0];
		lightVectors2[1] = lightDir0[1] * normal0[1];
		lightVectors2[2] = lightDir0[2] * normal0[2];
		lightVectors2[3] = lightDir0[3] * normal0[3];

		lightVectors2[0] += lightDir1[0] * normal1[0];
		lightVectors2[1] += lightDir1[1] * normal1[1];
		lightVectors2[2] += lightDir1[2] * normal1[2];
		lightVectors2[3] += lightDir1[3] * normal1[3];

		lightVectors2[0] += lightDir2[0] * normal2[0];
		lightVectors2[1] += lightDir2[1] * normal2[1];
		lightVectors2[2] += lightDir2[2] * normal2[2];
		lightVectors2[3] += lightDir2[3] * normal2[3];


		for ( int j = 0; j < 4; j++ ) {
			int n = usedVertNums[j];

			lightVectors[n][0] = lightVectors0[j];
			lightVectors[n][1] = lightVectors1[j];
			lightVectors[n][2] = lightVectors2[j];
		}

		numUsedVerts = 0;
	}

	for ( int i = 0; i < numUsedVerts; i++ ) {

		lightVectors0[i] = lightDir0[i] * tangent0[i] + lightDir1[i] * tangent1[i] + lightDir2[i] * tangent2[i];
		lightVectors1[i] = lightDir0[i] * tangent3[i] + lightDir1[i] * tangent4[i] + lightDir2[i] * tangent5[i];
		lightVectors2[i] = lightDir0[i] * normal0[i] + lightDir1[i] * normal1[i] + lightDir2[i] * normal2[i];

		int n = usedVertNums[i];
		lightVectors[n][0] = lightVectors0[i];
		lightVectors[n][1] = lightVectors1[i];
		lightVectors[n][2] = lightVectors2[i];
	}

#endif
}

/*
============
CSIMD_SSE:: createSpecularTextureCoords
============
*/
void VPCALL CSIMD_SSE:: createSpecularTextureCoords( CVec4D *texCoords, const CVec3D &lightOrigin, const CVec3D &viewOrigin, const CVertex *verts, const int numVerts, const int *indexes, const int numIndexes ) {

	SMF_ASSERT( sizeof( CVertex ) == DRAWVERT_SIZE );
	SMF_ASSERT( (int)&((CVertex *)0)->xyz == DRAWVERT_XYZ_OFFSET );
	SMF_ASSERT( (int)&((CVertex *)0)->normal == DRAWVERT_NORMAL_OFFSET );
	SMF_ASSERT( (int)&((CVertex *)0)->tangents[0] == DRAWVERT_TANGENT0_OFFSET );
	SMF_ASSERT( (int)&((CVertex *)0)->tangents[1] == DRAWVERT_TANGENT1_OFFSET );

	bool *used = (bool *)_alloca16( numVerts * sizeof( used[0] ) );
	memset( used, 0, numVerts * sizeof( used[0] ) );

	for ( int i = numIndexes - 1; i >= 0; i-- ) {
		used[indexes[i]] = true;
	}

#if 0

	__asm {

		mov			eax, numVerts

		mov			esi, used
		add			esi, eax

		mov			edi, verts
		sub			edi, DRAWVERT_SIZE

		neg			eax
		dec			eax

		mov			ecx, viewOrigin
		movss		xmm6, [ecx+0]
		movhps		xmm6, [ecx+4]

		mov			ecx, lightOrigin
		movss		xmm7, [ecx+0]
		movhps		xmm7, [ecx+4]

		mov			ecx, texCoords
		sub			ecx, 4*4

	loopVert:
		inc			eax
		jge			done

		add			edi, DRAWVERT_SIZE
		add			ecx, 4*4

		cmp			byte ptr [esi+eax], 0
		je			loopVert

		movaps		xmm0, xmm7
		movaps		xmm1, xmm6
		movss		xmm2, [edi+DRAWVERT_XYZ_OFFSET+0]
		movhps		xmm2, [edi+DRAWVERT_XYZ_OFFSET+4]
		subps		xmm0, xmm2
		subps		xmm1, xmm2

		movaps		xmm3, xmm0
		movaps		xmm4, xmm1
		mulps		xmm3, xmm3
		mulps		xmm4, xmm4

		// 0,  X,  1,  2
		// 3,  X,  4,  5

		movaps		xmm5, xmm3								// xmm5 = 0,  X,  1,  2
		unpcklps	xmm5, xmm4								// xmm5 = 0,  3,  X,  X
		unpckhps	xmm3, xmm4								// xmm3 = 1,  4,  2,  5
		movhlps		xmm4, xmm3								// xmm4 = 2,  5,  4,  5

		addps		xmm5, xmm3
		addps		xmm5, xmm4
		shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 1, 0, 1 )
		rsqrtps		xmm5, xmm5

		movaps		xmm4, xmm5
		shufps		xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )
		shufps		xmm5, xmm5, R_SHUFFLEPS( 1, 1, 1, 1 )

		mulps		xmm0, xmm4
		mulps		xmm1, xmm5
		addps		xmm0, xmm1

		movss		xmm2, [edi+DRAWVERT_TANGENT0_OFFSET+0]
		movhps		xmm2, [edi+DRAWVERT_TANGENT0_OFFSET+4]
		mulps		xmm2, xmm0

		movss		xmm3, [edi+DRAWVERT_TANGENT1_OFFSET+0]
		movhps		xmm3, [edi+DRAWVERT_TANGENT1_OFFSET+4]
		mulps		xmm3, xmm0

		movss		xmm4, [edi+DRAWVERT_NORMAL_OFFSET+0]
		movhps		xmm4, [edi+DRAWVERT_NORMAL_OFFSET+4]
		mulps		xmm4, xmm0

		movaps		xmm5, xmm2								// xmm5 = 0,  X,  1,  2
		unpcklps	xmm5, xmm3								// xmm5 = 0,  3,  X,  X
		unpckhps	xmm2, xmm3								// xmm2 = 1,  4,  2,  5

		movlhps		xmm5, xmm4								// xmm5 = 0,  3,  6,  X
		movhlps		xmm4, xmm2								// xmm4 = 2,  5,  7,  8
		shufps		xmm2, xmm4, R_SHUFFLEPS( 0, 1, 3, 2 )	// xmm2 = 2,  5,  8,  7

		movaps		xmm3, SIMD_SP_one

		addps		xmm5, xmm4
		addps		xmm5, xmm2
		movaps		[ecx+0], xmm5
		movss		[ecx+12], xmm3

		jmp			loopVert

	done:
	}

#elif 0

	for ( int i = 0; i < numVerts; i++ ) {
		if ( !used[i] ) {
			continue;
		}

		const CVertex *v = &verts[i];

		CVec3D lightDir = lightOrigin - v->xyz;
		CVec3D viewDir = viewOrigin - v->xyz;

		float ilength;

		ilength = CMath::rSqrt( lightDir[0] * lightDir[0] + lightDir[1] * lightDir[1] + lightDir[2] * lightDir[2] );
		lightDir[0] *= ilength;
		lightDir[1] *= ilength;
		lightDir[2] *= ilength;

		ilength = CMath::rSqrt( viewDir[0] * viewDir[0] + viewDir[1] * viewDir[1] + viewDir[2] * viewDir[2] );
		viewDir[0] *= ilength;
		viewDir[1] *= ilength;
		viewDir[2] *= ilength;

		lightDir += viewDir;

		texCoords[i][0] = lightDir[0] * v->tangents[0][0] + lightDir[1] * v->tangents[0][1] + lightDir[2] * v->tangents[0][2];
		texCoords[i][1] = lightDir[0] * v->tangents[1][0] + lightDir[1] * v->tangents[1][1] + lightDir[2] * v->tangents[1][2];
		texCoords[i][2] = lightDir[0] * v->normal[0] + lightDir[1] * v->normal[1] + lightDir[2] * v->normal[2];
		texCoords[i][3] = 1.0f;
	}


#elif 1

	ALIGN16( int usedVertNums[4] );
	ALIGN16( float lightDir0[4] );
	ALIGN16( float lightDir1[4] );
	ALIGN16( float lightDir2[4] );
	ALIGN16( float viewDir0[4] );
	ALIGN16( float viewDir1[4] );
	ALIGN16( float viewDir2[4] );
	ALIGN16( float normal0[4] );
	ALIGN16( float normal1[4] );
	ALIGN16( float normal2[4] );
	ALIGN16( float tangent0[4] );
	ALIGN16( float tangent1[4] );
	ALIGN16( float tangent2[4] );
	ALIGN16( float tangent3[4] );
	ALIGN16( float tangent4[4] );
	ALIGN16( float tangent5[4] );
	CVec3D localLightOrigin = lightOrigin;
	CVec3D localViewOrigin = viewOrigin;

	__asm {

		xor			ecx, ecx
		mov			eax, numVerts

		mov			esi, used
		add			esi, eax

		mov			edi, verts
		sub			edi, DRAWVERT_SIZE

		neg			eax
		dec			eax

	loopVert4:
		inc			eax
		jge			done4

		add			edi, DRAWVERT_SIZE

		cmp			byte ptr [esi+eax], 0
		je			loopVert4

		mov			usedVertNums[ecx*4], eax

		inc			ecx
		cmp			ecx, 4

		movss		xmm3, localLightOrigin[0]
		movss		xmm4, localLightOrigin[4]
		movss		xmm5, localLightOrigin[8]

		subss		xmm3, [edi+DRAWVERT_XYZ_OFFSET+0]
		subss		xmm4, [edi+DRAWVERT_XYZ_OFFSET+4]
		subss		xmm5, [edi+DRAWVERT_XYZ_OFFSET+8]

		movss		lightDir0[ecx*4-4], xmm3
		movss		lightDir1[ecx*4-4], xmm4
		movss		lightDir2[ecx*4-4], xmm5

		movss		xmm0, localViewOrigin[0]
		movss		xmm1, localViewOrigin[4]
		movss		xmm2, localViewOrigin[8]

		subss		xmm0, [edi+DRAWVERT_XYZ_OFFSET+0]
		subss		xmm1, [edi+DRAWVERT_XYZ_OFFSET+4]
		subss		xmm2, [edi+DRAWVERT_XYZ_OFFSET+8]

		movss		viewDir0[ecx*4-4], xmm0
		movss		viewDir1[ecx*4-4], xmm1
		movss		viewDir2[ecx*4-4], xmm2

		movss		xmm3, [edi+DRAWVERT_NORMAL_OFFSET+0]
		movss		xmm4, [edi+DRAWVERT_NORMAL_OFFSET+4]
		movss		xmm5, [edi+DRAWVERT_NORMAL_OFFSET+8]

		movss		normal0[ecx*4-4], xmm3
		movss		normal1[ecx*4-4], xmm4
		movss		normal2[ecx*4-4], xmm5

		movss		xmm0, [edi+DRAWVERT_TANGENT0_OFFSET+0]
		movss		xmm1, [edi+DRAWVERT_TANGENT0_OFFSET+4]
		movss		xmm2, [edi+DRAWVERT_TANGENT0_OFFSET+8]

		movss		tangent0[ecx*4-4], xmm0
		movss		tangent1[ecx*4-4], xmm1
		movss		tangent2[ecx*4-4], xmm2

		movss		xmm3, [edi+DRAWVERT_TANGENT1_OFFSET+0]
		movss		xmm4, [edi+DRAWVERT_TANGENT1_OFFSET+4]
		movss		xmm5, [edi+DRAWVERT_TANGENT1_OFFSET+8]

		movss		tangent3[ecx*4-4], xmm3
		movss		tangent4[ecx*4-4], xmm4
		movss		tangent5[ecx*4-4], xmm5

		jl			loopVert4

		movaps		xmm6, lightDir0
		movaps		xmm0, xmm6
		mulps		xmm6, xmm6
		movaps		xmm7, lightDir1
		movaps		xmm1, xmm7
		mulps		xmm7, xmm7
		addps		xmm6, xmm7
		movaps		xmm5, lightDir2
		movaps		xmm2, xmm5
		mulps		xmm5, xmm5
		addps		xmm6, xmm5
		rsqrtps		xmm6, xmm6

		mulps		xmm0, xmm6
		mulps		xmm1, xmm6
		mulps		xmm2, xmm6

		movaps		xmm3, viewDir0
		movaps		xmm7, xmm3
		mulps		xmm7, xmm7
		movaps		xmm4, viewDir1
		movaps		xmm6, xmm4
		mulps		xmm6, xmm6
		addps		xmm7, xmm6
		movaps		xmm5, viewDir2
		movaps		xmm6, xmm5
		mulps		xmm6, xmm6
		addps		xmm7, xmm6
		rsqrtps		xmm7, xmm7

		mulps		xmm3, xmm7
		addps		xmm0, xmm3
		mulps		xmm4, xmm7
		addps		xmm1, xmm4
		mulps		xmm5, xmm7
		addps		xmm2, xmm5

		movaps		xmm3, tangent0
		mulps		xmm3, xmm0
		movaps		xmm4, tangent1
		mulps		xmm4, xmm1
		addps		xmm3, xmm4
		movaps		xmm5, tangent2
		mulps		xmm5, xmm2
		addps		xmm5, xmm3

		movaps		xmm3, tangent3
		mulps		xmm3, xmm0
		movaps		xmm4, tangent4
		mulps		xmm4, xmm1
		addps		xmm3, xmm4
		movaps		xmm6, tangent5
		mulps		xmm6, xmm2
		addps		xmm6, xmm3

		mulps		xmm0, normal0
		mulps		xmm1, normal1
		addps		xmm0, xmm1
		mulps		xmm2, normal2
		addps		xmm0, xmm2

		mov			ecx, numVerts
		shl			ecx, 4
		mov			edx, usedVertNums[0]
		add			ecx, texCoords
		shl			edx, 4
		movss		xmm3, SIMD_SP_one

		movss		[ecx+edx+0], xmm5
		movss		[ecx+edx+4], xmm6
		movss		[ecx+edx+8], xmm0
		movss		[ecx+edx+12], xmm3

		shufps		xmm5, xmm5, R_SHUFFLEPS( 1, 2, 3, 0 )
		mov			edx, usedVertNums[4]
		shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 2, 3, 0 )
		shl			edx, 4
		shufps		xmm0, xmm0, R_SHUFFLEPS( 1, 2, 3, 0 )

		movss		[ecx+edx+0], xmm5
		movss		[ecx+edx+4], xmm6
		movss		[ecx+edx+8], xmm0
		movss		[ecx+edx+12], xmm3

		shufps		xmm5, xmm5, R_SHUFFLEPS( 1, 2, 3, 0 )
		mov			edx, usedVertNums[8]
		shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 2, 3, 0 )
		shl			edx, 4
		shufps		xmm0, xmm0, R_SHUFFLEPS( 1, 2, 3, 0 )

		movss		[ecx+edx+0], xmm5
		movss		[ecx+edx+4], xmm6
		movss		[ecx+edx+8], xmm0
		movss		[ecx+edx+12], xmm3

		shufps		xmm5, xmm5, R_SHUFFLEPS( 1, 2, 3, 0 )
		mov			edx, usedVertNums[12]
		shufps		xmm6, xmm6, R_SHUFFLEPS( 1, 2, 3, 0 )
		shl			edx, 4
		shufps		xmm0, xmm0, R_SHUFFLEPS( 1, 2, 3, 0 )

		movss		[ecx+edx+0], xmm5
		movss		[ecx+edx+4], xmm6
		movss		[ecx+edx+8], xmm0
		movss		[ecx+edx+12], xmm3

		xor			ecx, ecx
		jmp			loopVert4

	done4:
		test		ecx, ecx
		jz			done
		xor			eax, eax
		mov			edi, numVerts
		shl			edi, 4
		add			edi, texCoords

	loopVert1:
		movss		xmm6, lightDir0[eax*4]
		movss		xmm0, xmm6
		mulss		xmm6, xmm6
		movss		xmm7, lightDir1[eax*4]
		movss		xmm1, xmm7
		mulss		xmm7, xmm7
		addss		xmm6, xmm7
		movss		xmm5, lightDir2[eax*4]
		movss		xmm2, xmm5
		mulss		xmm5, xmm5
		addss		xmm6, xmm5
		rsqrtss		xmm6, xmm6

		mulss		xmm0, xmm6
		mulss		xmm1, xmm6
		mulss		xmm2, xmm6

		movss		xmm3, viewDir0[eax*4]
		movss		xmm7, xmm3
		mulss		xmm7, xmm7
		movss		xmm4, viewDir1[eax*4]
		movss		xmm6, xmm4
		mulss		xmm6, xmm6
		addss		xmm7, xmm6
		movss		xmm5, viewDir2[eax*4]
		movss		xmm6, xmm5
		mulss		xmm6, xmm6
		addss		xmm7, xmm6
		rsqrtss		xmm7, xmm7

		mulss		xmm3, xmm7
		addss		xmm0, xmm3
		mulss		xmm4, xmm7
		addss		xmm1, xmm4
		mulss		xmm5, xmm7
		addss		xmm2, xmm5

		mov			edx, usedVertNums[eax*4]
		shl			edx, 4

		movss		xmm3, tangent0[eax*4]
		mulss		xmm3, xmm0
		movss		xmm4, tangent1[eax*4]
		mulss		xmm4, xmm1
		addss		xmm3, xmm4
		movss		xmm5, tangent2[eax*4]
		mulss		xmm5, xmm2
		addss		xmm5, xmm3
		movss		[edi+edx+0], xmm5

		movss		xmm3, tangent3[eax*4]
		mulss		xmm3, xmm0
		movss		xmm4, tangent4[eax*4]
		mulss		xmm4, xmm1
		addss		xmm3, xmm4
		movss		xmm6, tangent5[eax*4]
		mulss		xmm6, xmm2
		addss		xmm6, xmm3
		movss		[edi+edx+4], xmm6

		mulss		xmm0, normal0[eax*4]
		mulss		xmm1, normal1[eax*4]
		addss		xmm0, xmm1
		mulss		xmm2, normal2[eax*4]
		addss		xmm0, xmm2
		movss		[edi+edx+8], xmm0

		movss		xmm3, SIMD_SP_one
		movss		[edi+edx+12], xmm3

		inc			eax
		dec			ecx
		jg			loopVert1

	done:
	}

#else

	ALIGN16( int usedVertNums[4] );
	ALIGN16( float lightDir0[4] );
	ALIGN16( float lightDir1[4] );
	ALIGN16( float lightDir2[4] );
	ALIGN16( float viewDir0[4] );
	ALIGN16( float viewDir1[4] );
	ALIGN16( float viewDir2[4] );
	ALIGN16( float normal0[4] );
	ALIGN16( float normal1[4] );
	ALIGN16( float normal2[4] );
	ALIGN16( float tangent0[4] );
	ALIGN16( float tangent1[4] );
	ALIGN16( float tangent2[4] );
	ALIGN16( float tangent3[4] );
	ALIGN16( float tangent4[4] );
	ALIGN16( float tangent5[4] );
	ALIGN16( float texCoords0[4] );
	ALIGN16( float texCoords1[4] );
	ALIGN16( float texCoords2[4] );
	CVec3D localLightOrigin = lightOrigin;
	CVec3D localViewOrigin = viewOrigin;
	int numUsedVerts = 0;

	for ( int i = 0; i < numVerts; i++ ) {
		if ( !used[i] ) {
			continue;
		}

		const CVertex *v = &verts[i];

		lightDir0[numUsedVerts] = localLightOrigin[0] - v->xyz[0];
		lightDir1[numUsedVerts] = localLightOrigin[1] - v->xyz[1];
		lightDir2[numUsedVerts] = localLightOrigin[2] - v->xyz[2];

		viewDir0[numUsedVerts] = localViewOrigin[0] - v->xyz[0];
		viewDir1[numUsedVerts] = localViewOrigin[1] - v->xyz[1];
		viewDir2[numUsedVerts] = localViewOrigin[2] - v->xyz[2];

		normal0[numUsedVerts] = v->normal[0];
		normal1[numUsedVerts] = v->normal[1];
		normal2[numUsedVerts] = v->normal[2];

		tangent0[numUsedVerts] = v->tangents[0][0];
		tangent1[numUsedVerts] = v->tangents[0][1];
		tangent2[numUsedVerts] = v->tangents[0][2];

		tangent3[numUsedVerts] = v->tangents[1][0];
		tangent4[numUsedVerts] = v->tangents[1][1];
		tangent5[numUsedVerts] = v->tangents[1][2];

		usedVertNums[numUsedVerts++] = i;
		if ( numUsedVerts < 4 ) {
			continue;
		}

		ALIGN16( float temp[4] );

		temp[0] = lightDir0[0] * lightDir0[0];
		temp[1] = lightDir0[1] * lightDir0[1];
		temp[2] = lightDir0[2] * lightDir0[2];
		temp[3] = lightDir0[3] * lightDir0[3];

		temp[0] += lightDir1[0] * lightDir1[0];
		temp[1] += lightDir1[1] * lightDir1[1];
		temp[2] += lightDir1[2] * lightDir1[2];
		temp[3] += lightDir1[3] * lightDir1[3];

		temp[0] += lightDir2[0] * lightDir2[0];
		temp[1] += lightDir2[1] * lightDir2[1];
		temp[2] += lightDir2[2] * lightDir2[2];
		temp[3] += lightDir2[3] * lightDir2[3];

		temp[0] = CMath::rSqrt( temp[0] );
		temp[1] = CMath::rSqrt( temp[1] );
		temp[2] = CMath::rSqrt( temp[2] );
		temp[3] = CMath::rSqrt( temp[3] );

		lightDir0[0] *= temp[0];
		lightDir0[1] *= temp[1];
		lightDir0[2] *= temp[2];
		lightDir0[3] *= temp[3];

		lightDir1[0] *= temp[0];
		lightDir1[1] *= temp[1];
		lightDir1[2] *= temp[2];
		lightDir1[3] *= temp[3];

		lightDir2[0] *= temp[0];
		lightDir2[1] *= temp[1];
		lightDir2[2] *= temp[2];
		lightDir2[3] *= temp[3];

		temp[0] = viewDir0[0] * viewDir0[0];
		temp[1] = viewDir0[1] * viewDir0[1];
		temp[2] = viewDir0[2] * viewDir0[2];
		temp[3] = viewDir0[3] * viewDir0[3];

		temp[0] += viewDir1[0] * viewDir1[0];
		temp[1] += viewDir1[1] * viewDir1[1];
		temp[2] += viewDir1[2] * viewDir1[2];
		temp[3] += viewDir1[3] * viewDir1[3];

		temp[0] += viewDir2[0] * viewDir2[0];
		temp[1] += viewDir2[1] * viewDir2[1];
		temp[2] += viewDir2[2] * viewDir2[2];
		temp[3] += viewDir2[3] * viewDir2[3];

		temp[0] = CMath::rSqrt( temp[0] );
		temp[1] = CMath::rSqrt( temp[1] );
		temp[2] = CMath::rSqrt( temp[2] );
		temp[3] = CMath::rSqrt( temp[3] );

		viewDir0[0] *= temp[0];
		viewDir0[1] *= temp[1];
		viewDir0[2] *= temp[2];
		viewDir0[3] *= temp[3];

		viewDir1[0] *= temp[0];
		viewDir1[1] *= temp[1];
		viewDir1[2] *= temp[2];
		viewDir1[3] *= temp[3];

		viewDir2[0] *= temp[0];
		viewDir2[1] *= temp[1];
		viewDir2[2] *= temp[2];
		viewDir2[3] *= temp[3];

		lightDir0[0] += viewDir0[0];
		lightDir0[1] += viewDir0[1];
		lightDir0[2] += viewDir0[2];
		lightDir0[3] += viewDir0[3];

		lightDir1[0] += viewDir1[0];
		lightDir1[1] += viewDir1[1];
		lightDir1[2] += viewDir1[2];
		lightDir1[3] += viewDir1[3];

		lightDir2[0] += viewDir2[0];
		lightDir2[1] += viewDir2[1];
		lightDir2[2] += viewDir2[2];
		lightDir2[3] += viewDir2[3];

		texCoords0[0] = lightDir0[0] * tangent0[0];
		texCoords0[1] = lightDir0[1] * tangent0[1];
		texCoords0[2] = lightDir0[2] * tangent0[2];
		texCoords0[3] = lightDir0[3] * tangent0[3];

		texCoords0[0] += lightDir1[0] * tangent1[0];
		texCoords0[1] += lightDir1[1] * tangent1[1];
		texCoords0[2] += lightDir1[2] * tangent1[2];
		texCoords0[3] += lightDir1[3] * tangent1[3];

		texCoords0[0] += lightDir2[0] * tangent2[0];
		texCoords0[1] += lightDir2[1] * tangent2[1];
		texCoords0[2] += lightDir2[2] * tangent2[2];
		texCoords0[3] += lightDir2[3] * tangent2[3];

		texCoords1[0] = lightDir0[0] * tangent3[0];
		texCoords1[1] = lightDir0[1] * tangent3[1];
		texCoords1[2] = lightDir0[2] * tangent3[2];
		texCoords1[3] = lightDir0[3] * tangent3[3];

		texCoords1[0] += lightDir1[0] * tangent4[0];
		texCoords1[1] += lightDir1[1] * tangent4[1];
		texCoords1[2] += lightDir1[2] * tangent4[2];
		texCoords1[3] += lightDir1[3] * tangent4[3];

		texCoords1[0] += lightDir2[0] * tangent5[0];
		texCoords1[1] += lightDir2[1] * tangent5[1];
		texCoords1[2] += lightDir2[2] * tangent5[2];
		texCoords1[3] += lightDir2[3] * tangent5[3];

		texCoords2[0] = lightDir0[0] * normal0[0];
		texCoords2[1] = lightDir0[1] * normal0[1];
		texCoords2[2] = lightDir0[2] * normal0[2];
		texCoords2[3] = lightDir0[3] * normal0[3];

		texCoords2[0] += lightDir1[0] * normal1[0];
		texCoords2[1] += lightDir1[1] * normal1[1];
		texCoords2[2] += lightDir1[2] * normal1[2];
		texCoords2[3] += lightDir1[3] * normal1[3];

		texCoords2[0] += lightDir2[0] * normal2[0];
		texCoords2[1] += lightDir2[1] * normal2[1];
		texCoords2[2] += lightDir2[2] * normal2[2];
		texCoords2[3] += lightDir2[3] * normal2[3];

		for ( int j = 0; j < 4; j++ ) {
			int n = usedVertNums[j];

			texCoords[n][0] = texCoords0[j];
			texCoords[n][1] = texCoords1[j];
			texCoords[n][2] = texCoords2[j];
			texCoords[n][3] = 1.0f;
		}

		numUsedVerts = 0;
	}

	for ( int i = 0; i < numUsedVerts; i++ ) {
		float temp;

		temp = lightDir0[i] * lightDir0[i] + lightDir1[i] * lightDir1[i] + lightDir2[i] * lightDir2[i];
		temp = CMath::rSqrt( temp );

		lightDir0[i] *= temp;
		lightDir1[i] *= temp;
		lightDir2[i] *= temp;

		temp = viewDir0[i] * viewDir0[i] + viewDir1[i] * viewDir1[i] + viewDir2[i] * viewDir2[i];
		temp = CMath::rSqrt( temp );

		viewDir0[i] *= temp;
		viewDir1[i] *= temp;
		viewDir2[i] *= temp;

		lightDir0[i] += viewDir0[i];
		lightDir1[i] += viewDir1[i];
		lightDir2[i] += viewDir2[i];

		texCoords0[i] = lightDir0[i] * tangent0[i] + lightDir1[i] * tangent1[i] + lightDir2[i] * tangent2[i];
		texCoords1[i] = lightDir0[i] * tangent3[i] + lightDir1[i] * tangent4[i] + lightDir2[i] * tangent5[i];
		texCoords2[i] = lightDir0[i] * normal0[i] + lightDir1[i] * normal1[i] + lightDir2[i] * normal2[i];

		int n = usedVertNums[i];
		texCoords[n][0] = texCoords0;
		texCoords[n][1] = texCoords1;
		texCoords[n][2] = texCoords2;
		texCoords[n][3] = 1.0f;
	}

#endif
}

/*
============
CSIMD_SSE:: createShadowCache
============
*/
int VPCALL CSIMD_SSE:: createShadowCache( CVec4D *vertexCache, int *vertRemap, const CVec3D &lightOrigin, const CVertex *verts, const int numVerts ) {
#if 1
	int outVerts;

	__asm {
		push		ebx

		mov			esi, lightOrigin
		movaps		xmm5, SIMD_SP_lastOne
		movss		xmm6, [esi+0]
		movhps		xmm6, [esi+4]
		shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 2, 3, 1 )
		orps		xmm6, SIMD_SP_lastOne
		movaps		xmm7, xmm6

		xor			ebx, ebx
		xor			ecx, ecx

		mov			edx, vertRemap
		mov			esi, verts
		mov			edi, vertexCache
		mov			eax, numVerts
		and			eax, ~3
		jz			done4
		shl			eax, 2
		add			edx, eax
		neg			eax

	loop4:
		prefetchnta	[edx+128]
		prefetchnta	[esi+4*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET]

		cmp         dword ptr [edx+eax+0], ebx
		jne         skip1

		mov			dword ptr [edx+eax+0], ecx
		movss		xmm0, [esi+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+8]
		movhps		xmm0, [esi+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]
		add			ecx, 2
		shufps		xmm0, xmm0, R_SHUFFLEPS( 2, 3, 0, 1 );
		orps		xmm0, xmm5
		movaps		[edi+0*16], xmm0
		subps		xmm0, xmm6
		movaps		[edi+1*16], xmm0
		add			edi, 2*16

	skip1:
		cmp         dword ptr [edx+eax+4], ebx
		jne         skip2

		mov			dword ptr [edx+eax+4], ecx
		movss		xmm1, [esi+1*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]
		movhps		xmm1, [esi+1*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+4]
		add			ecx, 2
		shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 2, 3, 1 )
		orps		xmm1, xmm5
		movaps		[edi+0*16], xmm1
		subps		xmm1, xmm7
		movaps		[edi+1*16], xmm1
		add			edi, 2*16

	skip2:
		cmp         dword ptr [edx+eax+8], ebx
		jne         skip3

		mov			dword ptr [edx+eax+8], ecx
		movss		xmm2, [esi+2*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+8]
		movhps		xmm2, [esi+2*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]
		add			ecx, 2
		shufps		xmm2, xmm2, R_SHUFFLEPS( 2, 3, 0, 1 );
		orps		xmm2, xmm5
		movaps		[edi+0*16], xmm2
		subps		xmm2, xmm6
		movaps		[edi+1*16], xmm2
		add			edi, 2*16

	skip3:
		cmp         dword ptr [edx+eax+12], ebx
		jne         skip4

		mov			dword ptr [edx+eax+12], ecx
		movss		xmm3, [esi+3*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]
		movhps		xmm3, [esi+3*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+4]
		add			ecx, 2
		shufps		xmm3, xmm3, R_SHUFFLEPS( 0, 2, 3, 1 )
		orps		xmm3, xmm5
		movaps		[edi+0*16], xmm3
		subps		xmm3, xmm7
		movaps		[edi+1*16], xmm3
		add			edi, 2*16

	skip4:
		add			esi, 4*DRAWVERT_SIZE
		add			eax, 4*4
		jl			loop4

	done4:
		mov			eax, numVerts
		and			eax, 3
		jz			done1
		shl			eax, 2
		add			edx, eax
		neg			eax

	loop1:
		cmp         dword ptr [edx+eax+0], ebx
		jne         skip0

		mov			dword ptr [edx+eax+0], ecx
		movss		xmm0, [esi+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+8]
		movhps		xmm0, [esi+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]
		add			ecx, 2
		shufps		xmm0, xmm0, R_SHUFFLEPS( 2, 3, 0, 1 )
		orps		xmm0, xmm5
		movaps		[edi+0*16], xmm0
		subps		xmm0, xmm6
		movaps		[edi+1*16], xmm0
		add			edi, 2*16

	skip0:

		add			esi, DRAWVERT_SIZE
		add			eax, 4
		jl			loop1

	done1:
		pop			ebx
		mov			outVerts, ecx
	}
	return outVerts;

#else

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

#endif
}

/*
============
CSIMD_SSE:: createVertexProgramShadowCache
============
*/
int VPCALL CSIMD_SSE:: createVertexProgramShadowCache( CVec4D *vertexCache, const CVertex *verts, const int numVerts ) {
#if 1

	__asm {
		movaps		xmm4, SIMD_SP_lastOne
		movaps		xmm5, xmm4
		movaps		xmm6, xmm4
		movaps		xmm7, xmm4

		mov			esi, verts
		mov			edi, vertexCache
		mov			eax, numVerts
		and			eax, ~3
		jz			done4
		shl			eax, 5
		add			edi, eax
		neg			eax

	loop4:
		prefetchnta	[esi+4*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET]

		movss		xmm0, [esi+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+8]
		movhps		xmm0, [esi+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]
		shufps		xmm0, xmm0, R_SHUFFLEPS( 2, 3, 0, 1 );
		movaps		[edi+eax+1*16], xmm0
		orps		xmm0, xmm4
		movaps		[edi+eax+0*16], xmm0

		movss		xmm1, [esi+1*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]
		movhps		xmm1, [esi+1*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+4]
		shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 2, 3, 1 )
		movaps		[edi+eax+3*16], xmm1
		orps		xmm1, xmm5
		movaps		[edi+eax+2*16], xmm1

		movss		xmm2, [esi+2*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+8]
		movhps		xmm2, [esi+2*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]
		shufps		xmm2, xmm2, R_SHUFFLEPS( 2, 3, 0, 1 );
		movaps		[edi+eax+5*16], xmm2
		orps		xmm2, xmm6
		movaps		[edi+eax+4*16], xmm2

		movss		xmm3, [esi+3*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]
		movhps		xmm3, [esi+3*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+4]
		shufps		xmm3, xmm3, R_SHUFFLEPS( 0, 2, 3, 1 )
		movaps		[edi+eax+7*16], xmm3
		orps		xmm3, xmm7
		movaps		[edi+eax+6*16], xmm3

		add			esi, 4*DRAWVERT_SIZE
		add			eax, 4*8*4
		jl			loop4

	done4:
		mov			eax, numVerts
		and			eax, 3
		jz			done1
		shl			eax, 5
		add			edi, eax
		neg			eax

	loop1:
		movss		xmm0, [esi+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+8]
		movhps		xmm0, [esi+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]
		shufps		xmm0, xmm0, R_SHUFFLEPS( 2, 3, 0, 1 );
		movaps		[edi+eax+1*16], xmm0
		orps		xmm0, xmm4
		movaps		[edi+eax+0*16], xmm0

		add			esi, DRAWVERT_SIZE
		add			eax, 8*4
		jl			loop1

	done1:
	}
	return numVerts * 2;

#else

	for ( int i = 0; i < numVerts; i++ ) {
		const float *v = verts[i].xyz.toFloatPtr();
		vertexCache[i*2+0][0] = v[0];
		vertexCache[i*2+0][1] = v[1];
		vertexCache[i*2+0][2] = v[2];
		vertexCache[i*2+0][3] = 1.0f;

		vertexCache[i*2+1][0] = v[0];
		vertexCache[i*2+1][1] = v[1];
		vertexCache[i*2+1][2] = v[2];
		vertexCache[i*2+1][3] = 0.0f;
	}
	return numVerts * 2;

#endif
}

/*
============
SSE_UpSample11kHzMonoPCMTo44kHz
============
*/
static void SSE_UpSample11kHzMonoPCMTo44kHz( float *dest, const short *src, const int numSamples ) {
	__asm {
		mov			esi, src
		mov			edi, dest

		mov			eax, numSamples
		and			eax, ~1
		jz			done2
		shl			eax, 1
		add			esi, eax
		neg			eax

		align		16
	loop2:
		add			edi, 2*4*4

		movsx		ecx, word ptr [esi+eax+0]
		cvtsi2ss	xmm0, ecx
		shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
		movlps		[edi-2*4*4+0], xmm0
		movhps		[edi-2*4*4+8], xmm0

		movsx		edx, word ptr [esi+eax+2]
		cvtsi2ss	xmm1, edx
		shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
		movlps		[edi-1*4*4+0], xmm1
		movhps		[edi-1*4*4+8], xmm1

		add			eax, 2*2
		jl			loop2

	done2:
		mov			eax, numSamples
		and			eax, 1
		jz			done

		movsx		ecx, word ptr [esi]
		cvtsi2ss	xmm0, ecx
		shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
		movlps		[edi+0], xmm0
		movhps		[edi+8], xmm0

	done:
	}
}

/*
============
SSE_UpSample11kHzStereoPCMTo44kHz
============
*/
static void SSE_UpSample11kHzStereoPCMTo44kHz( float *dest, const short *src, const int numSamples ) {
	__asm {
		mov			esi, src
		mov			edi, dest

		mov			eax, numSamples
		test		eax, ~1
		jz			done2
		shl			eax, 1
		add			esi, eax
		neg			eax

		align		16
	loop2:
		add			edi, 8*4

		movsx		ecx, word ptr [esi+eax+0]
		cvtsi2ss	xmm0, ecx

		movsx		edx, word ptr [esi+eax+2]
		cvtsi2ss	xmm1, edx

		unpcklps	xmm0, xmm1

		movlps		[edi-8*4+0], xmm0
		movlps		[edi-8*4+8], xmm0
		movlps		[edi-4*4+0], xmm0
		movlps		[edi-4*4+8], xmm0

		add			eax, 2*2
		jl			loop2

	done2:
	}
}

/*
============
SSE_UpSample22kHzMonoPCMTo44kHz
============
*/
static void SSE_UpSample22kHzMonoPCMTo44kHz( float *dest, const short *src, const int numSamples ) {
	__asm {
		mov			esi, src
		mov			edi, dest

		mov			eax, numSamples
		and			eax, ~1
		jz			done2
		shl			eax, 1
		add			esi, eax
		neg			eax

		align		16
	loop2:
		add			edi, 4*4

		movsx		ecx, word ptr [esi+eax+0]
		cvtsi2ss	xmm0, ecx

		movsx		edx, word ptr [esi+eax+2]
		cvtsi2ss	xmm1, edx

		shufps		xmm0, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
		movlps		[edi-4*4+0], xmm0
		movhps		[edi-4*4+8], xmm0

		add			eax, 2*2
		jl			loop2

	done2:
		mov			eax, numSamples
		and			eax, 1
		jz			done

		movsx		ecx, word ptr [esi]
		cvtsi2ss	xmm0, ecx
		shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
		movlps		[edi], xmm0

	done:
	}
}

/*
============
SSE_UpSample22kHzStereoPCMTo44kHz
============
*/
static void SSE_UpSample22kHzStereoPCMTo44kHz( float *dest, const short *src, const int numSamples ) {
	__asm {
		mov			esi, src
		mov			edi, dest

		mov			eax, numSamples
		test		eax, ~1
		jz			done2
		shl			eax, 1
		add			esi, eax
		neg			eax

		align		16
	loop2:
		add			edi, 4*4

		movsx		ecx, word ptr [esi+eax+0]
		cvtsi2ss	xmm0, ecx
		movss		[edi-4*4], xmm0
		movss		[edi-2*4], xmm0

		movsx		edx, word ptr [esi+eax+2]
		cvtsi2ss	xmm1, edx
		movss		[edi-3*4], xmm1
		movss		[edi-1*4], xmm1

		add			eax, 2*2
		jl			loop2

	done2:
	}
}

/*
============
SSE_UpSample44kHzMonoPCMTo44kHz
============
*/
static void SSE_UpSample44kHzMonoPCMTo44kHz( float *dest, const short *src, const int numSamples ) {
	__asm {
		mov			esi, src
		mov			edi, dest

		mov			eax, numSamples
		and			eax, ~1
		jz			done2
		shl			eax, 1
		add			esi, eax
		neg			eax

		align		16
	loop2:
		add			edi, 2*4

		movsx		ecx, word ptr [esi+eax+0]
		cvtsi2ss	xmm0, ecx
		movss		[edi-2*4], xmm0

		movsx		edx, word ptr [esi+eax+2]
		cvtsi2ss	xmm1, edx
		movss		[edi-1*4], xmm1

		add			eax, 2*2
		jl			loop2

	done2:
		mov			eax, numSamples
		and			eax, 1
		jz			done

		movsx		ecx, word ptr [esi]
		cvtsi2ss	xmm0, ecx
		movss		[edi], xmm0

	done:
	}
}

/*
============
CSIMD_SSE:: upSamplePCMTo44kHz

  Duplicate samples for 44kHz output.
============
*/
void CSIMD_SSE:: upSamplePCMTo44kHz( float *dest, const short *src, const int numSamples, const int kHz, const int numChannels ) {
	if ( kHz == 11025 ) {
		if ( numChannels == 1 ) {
			SSE_UpSample11kHzMonoPCMTo44kHz( dest, src, numSamples );
		} else {
			SSE_UpSample11kHzStereoPCMTo44kHz( dest, src, numSamples );
		}
	} else if ( kHz == 22050 ) {
		if ( numChannels == 1 ) {
			SSE_UpSample22kHzMonoPCMTo44kHz( dest, src, numSamples );
		} else {
			SSE_UpSample22kHzStereoPCMTo44kHz( dest, src, numSamples );
		}
	} else if ( kHz == 44100 ) {
		SSE_UpSample44kHzMonoPCMTo44kHz( dest, src, numSamples );
	} else {
		SMF_ASSERT( 0 );
	}
}

/*
============
SSE_UpSample11kHzMonoOGGTo44kHz
============
*/
static void SSE_UpSample11kHzMonoOGGTo44kHz( float *dest, const float *src, const int numSamples ) {
	float constant = 32768.0f;
	__asm {
		mov			esi, src
		mov			edi, dest
		movss		xmm7, constant
		shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )

		mov			eax, numSamples
		and			eax, ~1
		jz			done2
		shl			eax, 2
		add			esi, eax
		neg			eax

		align		16
	loop2:
		add			edi, 2*16

		movss		xmm0, [esi+eax+0]
		mulss		xmm0, xmm7
		shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
		movlps		[edi-32], xmm0
		movlps		[edi-24], xmm0

		movss		xmm1, [esi+eax+4]
		mulss		xmm1, xmm7
		shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
		movlps		[edi-16], xmm1
		movlps		[edi- 8], xmm1

		add			eax, 2*4
		jl			loop2

	done2:
		mov			eax, numSamples
		and			eax, 1
		jz			done

		movss		xmm0, [esi]
		mulss		xmm0, xmm7
		shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
		movlps		[edi+0], xmm0
		movlps		[edi+8], xmm0

	done:
	}
}

/*
============
SSE_UpSample11kHzStereoOGGTo44kHz
============
*/
static void SSE_UpSample11kHzStereoOGGTo44kHz( float *dest, const float * const *src, const int numSamples ) {
	float constant = 32768.0f;
	__asm {
		mov			esi, src
		mov			ecx, [esi+0]
		mov			edx, [esi+4]
		mov			edi, dest
		movss		xmm7, constant
		shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )

		mov			eax, numSamples
		and			eax, ~1
		jz			done2
		shl			eax, 1
		add			ecx, eax
		add			edx, eax
		neg			eax

		align		16
	loop2:
		add			edi, 4*16

		movlps		xmm0, [ecx+eax]
		movlps		xmm1, [edx+eax]
		unpcklps	xmm0, xmm1
		mulps		xmm0, xmm7
		movlps		[edi-8*8], xmm0
		movlps		[edi-7*8], xmm0
		movlps		[edi-6*8], xmm0
		movlps		[edi-5*8], xmm0
		movhps		[edi-4*8], xmm0
		movhps		[edi-3*8], xmm0
		movhps		[edi-2*8], xmm0
		movhps		[edi-1*8], xmm0

		add			eax, 2*4
		jl			loop2

	done2:
		mov			eax, numSamples
		and			eax, 1
		jz			done

		movss		xmm0, [ecx]
		movss		xmm1, [edx]
		unpcklps	xmm0, xmm1
		mulps		xmm0, xmm7
		movlps		[edi+0*8], xmm0
		movlps		[edi+1*8], xmm0
		movlps		[edi+2*8], xmm0
		movlps		[edi+3*8], xmm0

	done:
	}
}

/*
============
SSE_UpSample22kHzMonoOGGTo44kHz
============
*/
static void SSE_UpSample22kHzMonoOGGTo44kHz( float *dest, const float *src, const int numSamples ) {
	float constant = 32768.0f;
	__asm {
		mov			esi, src
		mov			edi, dest
		movss		xmm7, constant
		shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )

		mov			eax, numSamples
		and			eax, ~1
		jz			done2
		shl			eax, 2
		add			esi, eax
		neg			eax

		align		16
	loop2:
		add			edi, 2*8

		movss		xmm0, [esi+eax+0]
		movss		xmm1, [esi+eax+4]
		shufps		xmm0, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm0, xmm7
		movlps		[edi-16], xmm0
		movhps		[edi- 8], xmm0

		add			eax, 2*4
		jl			loop2

	done2:
		mov			eax, numSamples
		and			eax, 1
		jz			done

		movss		xmm0, [esi]
		mulss		xmm0, xmm7
		shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
		movlps		[edi+0], xmm0

	done:
	}
}

/*
============
SSE_UpSample22kHzStereoOGGTo44kHz
============
*/
static void SSE_UpSample22kHzStereoOGGTo44kHz( float *dest, const float * const *src, const int numSamples ) {
	float constant = 32768.0f;
	__asm {
		mov			esi, src
		mov			ecx, [esi+0]
		mov			edx, [esi+4]
		mov			edi, dest
		movss		xmm7, constant
		shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )

		mov			eax, numSamples
		and			eax, ~1
		jz			done2
		shl			eax, 1
		add			ecx, eax
		add			edx, eax
		neg			eax

		align		16
	loop2:
		add			edi, 2*16

		movlps		xmm0, [ecx+eax]
		movlps		xmm1, [edx+eax]
		unpcklps	xmm0, xmm1
		mulps		xmm0, xmm7
		movlps		[edi-4*8], xmm0
		movlps		[edi-3*8], xmm0
		movhps		[edi-2*8], xmm0
		movhps		[edi-1*8], xmm0

		add			eax, 2*4
		jl			loop2

	done2:
		mov			eax, numSamples
		and			eax, 1
		jz			done

		movss		xmm0, [ecx]
		movss		xmm1, [edx]
		unpcklps	xmm0, xmm1
		mulps		xmm0, xmm7
		movlps		[edi+0*8], xmm0
		movlps		[edi+1*8], xmm0

	done:
	}
}

/*
============
SSE_UpSample44kHzMonoOGGTo44kHz
============
*/
static void SSE_UpSample44kHzMonoOGGTo44kHz( float *dest, const float *src, const int numSamples ) {
	float constant = 32768.0f;
	KFLOAT_CA( mul, dest, src, constant, numSamples )
}

/*
============
SSE_UpSample44kHzStereoOGGTo44kHz
============
*/
static void SSE_UpSample44kHzStereoOGGTo44kHz( float *dest, const float * const *src, const int numSamples ) {
	float constant = 32768.0f;
	__asm {
		mov			esi, src
		mov			ecx, [esi+0]
		mov			edx, [esi+4]
		mov			edi, dest
		movss		xmm7, constant
		shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )

		mov			eax, numSamples
		and			eax, ~1
		jz			done2
		shl			eax, 1
		add			ecx, eax
		add			edx, eax
		neg			eax

		align		16
	loop2:
		add			edi, 16

		movlps		xmm0, [ecx+eax]
		movlps		xmm1, [edx+eax]
		unpcklps	xmm0, xmm1
		mulps		xmm0, xmm7
		movlps		[edi-2*8], xmm0
		movhps		[edi-1*8], xmm0

		add			eax, 2*4
		jl			loop2

	done2:
		mov			eax, numSamples
		and			eax, 1
		jz			done

		movss		xmm0, [ecx]
		movss		xmm1, [edx]
		unpcklps	xmm0, xmm1
		mulps		xmm0, xmm7
		movlps		[edi+0*8], xmm0

	done:
	}
}

/*
============
CSIMD_SSE:: upSampleOGGTo44kHz

  Duplicate samples for 44kHz output.
============
*/
void CSIMD_SSE:: upSampleOGGTo44kHz( float *dest, const float * const *ogg, const int numSamples, const int kHz, const int numChannels ) {
	if ( kHz == 11025 ) {
		if ( numChannels == 1 ) {
			SSE_UpSample11kHzMonoOGGTo44kHz( dest, ogg[0], numSamples );
		} else {
			SSE_UpSample11kHzStereoOGGTo44kHz( dest, ogg, numSamples );
		}
	} else if ( kHz == 22050 ) {
		if ( numChannels == 1 ) {
			SSE_UpSample22kHzMonoOGGTo44kHz( dest, ogg[0], numSamples );
		} else {
			SSE_UpSample22kHzStereoOGGTo44kHz( dest, ogg, numSamples );
		}
	} else if ( kHz == 44100 ) {
		if ( numChannels == 1 ) {
			SSE_UpSample44kHzMonoOGGTo44kHz( dest, ogg[0], numSamples );
		} else {
			SSE_UpSample44kHzStereoOGGTo44kHz( dest, ogg, numSamples );
		}
	} else {
		SMF_ASSERT( 0 );
	}
}

/*
============
CSIMD_SSE:: mixSoundTwoSpeakerMono
============
*/
void VPCALL CSIMD_SSE:: mixSoundTwoSpeakerMono( float *mixBuffer, const float *samples, const int numSamples, const float lastV[2], const float currentV[2] ) {
#if 1

	ALIGN16( float incs[2] );

	SMF_ASSERT( numSamples == MIXBUFFER_SAMPLES );

	incs[0] = ( currentV[0] - lastV[0] ) / MIXBUFFER_SAMPLES;
	incs[1] = ( currentV[1] - lastV[1] ) / MIXBUFFER_SAMPLES;

	__asm {
		mov			eax, MIXBUFFER_SAMPLES
		mov			edi, mixBuffer
		mov			esi, samples
		shl			eax, 2
		add			esi, eax
		neg			eax

		mov			ecx, lastV
		movlps		xmm6, [ecx]
		xorps		xmm7, xmm7
		movhps		xmm7, incs
		shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 1, 0, 1 )
		addps		xmm6, xmm7
		shufps		xmm7, xmm7, R_SHUFFLEPS( 2, 3, 2, 3 )
		addps		xmm7, xmm7

	loop16:
		add			edi, 4*4*4

		movaps		xmm0, [esi+eax+0*4*4]
		movaps		xmm1, xmm0
		shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 1, 1 )
		mulps		xmm0, xmm6
		addps		xmm0, [edi-4*4*4]
		addps		xmm6, xmm7
		movaps		[edi-4*4*4], xmm0

		shufps		xmm1, xmm1, R_SHUFFLEPS( 2, 2, 3, 3 )
		mulps		xmm1, xmm6
		addps		xmm1, [edi-3*4*4]
		addps		xmm6, xmm7
		movaps		[edi-3*4*4], xmm1

		movaps		xmm2, [esi+eax+1*4*4]
		movaps		xmm3, xmm2
		shufps		xmm2, xmm2, R_SHUFFLEPS( 0, 0, 1, 1 )
		mulps		xmm2, xmm6
		addps		xmm2, [edi-2*4*4]
		addps		xmm6, xmm7
		movaps		[edi-2*4*4], xmm2

		shufps		xmm3, xmm3, R_SHUFFLEPS( 2, 2, 3, 3 )
		mulps		xmm3, xmm6
		addps		xmm3, [edi-1*4*4]
		addps		xmm6, xmm7
		movaps		[edi-1*4*4], xmm3

		add			eax, 2*4*4

		jl			loop16
	}

#else

	int i;
	float incL;
	float incR;
	float sL0, sL1;
	float sR0, sR1;

	SMF_ASSERT( numSamples == MIXBUFFER_SAMPLES );

	incL = ( currentV[0] - lastV[0] ) / MIXBUFFER_SAMPLES;
	incR = ( currentV[1] - lastV[1] ) / MIXBUFFER_SAMPLES;

	sL0 = lastV[0];
	sR0 = lastV[1];
	sL1 = lastV[0] + incL;
	sR1 = lastV[1] + incR;

	incL *= 2;
	incR *= 2;

	for( i = 0; i < MIXBUFFER_SAMPLES; i += 2 ) {
		mixBuffer[i*2+0] += samples[i+0] * sL0;
		mixBuffer[i*2+1] += samples[i+0] * sR0;
		mixBuffer[i*2+2] += samples[i+1] * sL1;
		mixBuffer[i*2+3] += samples[i+1] * sR1;
		sL0 += incL;
		sR0 += incR;
		sL1 += incL;
		sR1 += incR;
	}

#endif
}

/*
============
CSIMD_SSE:: mixSoundTwoSpeakerStereo
============
*/
void VPCALL CSIMD_SSE:: mixSoundTwoSpeakerStereo( float *mixBuffer, const float *samples, const int numSamples, const float lastV[2], const float currentV[2] ) {
#if 1

	ALIGN16( float incs[2] );

	SMF_ASSERT( numSamples == MIXBUFFER_SAMPLES );

	incs[0] = ( currentV[0] - lastV[0] ) / MIXBUFFER_SAMPLES;
	incs[1] = ( currentV[1] - lastV[1] ) / MIXBUFFER_SAMPLES;

	__asm {
		mov			eax, MIXBUFFER_SAMPLES
		mov			edi, mixBuffer
		mov			esi, samples
		shl			eax, 3
		add			esi, eax
		neg			eax

		mov			ecx, lastV
		movlps		xmm6, [ecx]
		xorps		xmm7, xmm7
		movhps		xmm7, incs
		shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 1, 0, 1 )
		addps		xmm6, xmm7
		shufps		xmm7, xmm7, R_SHUFFLEPS( 2, 3, 2, 3 )
		addps		xmm7, xmm7

	loop16:
		add			edi, 4*4*4

		movaps		xmm0, [esi+eax+0*4*4]
		mulps		xmm0, xmm6
		addps		xmm0, [edi-4*4*4]
		addps		xmm6, xmm7
		movaps		[edi-4*4*4], xmm0

		movaps		xmm2, [esi+eax+1*4*4]
		mulps		xmm2, xmm6
		addps		xmm2, [edi-3*4*4]
		addps		xmm6, xmm7
		movaps		[edi-3*4*4], xmm2

		movaps		xmm3, [esi+eax+2*4*4]
		mulps		xmm3, xmm6
		addps		xmm3, [edi-2*4*4]
		addps		xmm6, xmm7
		movaps		[edi-2*4*4], xmm3

		movaps		xmm4, [esi+eax+3*4*4]
		mulps		xmm4, xmm6
		addps		xmm4, [edi-1*4*4]
		addps		xmm6, xmm7
		movaps		[edi-1*4*4], xmm4

		add			eax, 4*4*4

		jl			loop16
	}

#else

	int i;
	float incL;
	float incR;
	float sL0, sL1;
	float sR0, sR1;

	SMF_ASSERT( numSamples == MIXBUFFER_SAMPLES );

	incL = ( currentV[0] - lastV[0] ) / MIXBUFFER_SAMPLES;
	incR = ( currentV[1] - lastV[1] ) / MIXBUFFER_SAMPLES;

	sL0 = lastV[0];
	sR0 = lastV[1];
	sL1 = lastV[0] + incL;
	sR1 = lastV[1] + incR;

	incL *= 2;
	incR *= 2;

	for( i = 0; i < MIXBUFFER_SAMPLES; i += 2 ) {
		mixBuffer[i*2+0] += samples[i*2+0] * sL0;
		mixBuffer[i*2+1] += samples[i*2+1] * sR0;
		mixBuffer[i*2+2] += samples[i*2+2] * sL1;
		mixBuffer[i*2+3] += samples[i*2+3] * sR1;
		sL0 += incL;
		sR0 += incR;
		sL1 += incL;
		sR1 += incR;
	}

#endif
}

/*
============
CSIMD_SSE:: mixSoundSixSpeakerMono
============
*/
void VPCALL CSIMD_SSE:: mixSoundSixSpeakerMono( float *mixBuffer, const float *samples, const int numSamples, const float lastV[6], const float currentV[6] ) {
#if 1

	ALIGN16( float incs[6] );

	SMF_ASSERT( numSamples == MIXBUFFER_SAMPLES );

	incs[0] = ( currentV[0] - lastV[0] ) / MIXBUFFER_SAMPLES;
	incs[1] = ( currentV[1] - lastV[1] ) / MIXBUFFER_SAMPLES;
	incs[2] = ( currentV[2] - lastV[2] ) / MIXBUFFER_SAMPLES;
	incs[3] = ( currentV[3] - lastV[3] ) / MIXBUFFER_SAMPLES;
	incs[4] = ( currentV[4] - lastV[4] ) / MIXBUFFER_SAMPLES;
	incs[5] = ( currentV[5] - lastV[5] ) / MIXBUFFER_SAMPLES;

	__asm {
		mov			eax, MIXBUFFER_SAMPLES
		mov			edi, mixBuffer
		mov			esi, samples
		shl			eax, 2
		add			esi, eax
		neg			eax

		mov			ecx, lastV
		movlps		xmm2, [ecx+ 0]
		movhps		xmm2, [ecx+ 8]
		movlps		xmm3, [ecx+16]
		movaps		xmm4, xmm2
		shufps		xmm3, xmm2, R_SHUFFLEPS( 0, 1, 0, 1 )
		shufps		xmm4, xmm3, R_SHUFFLEPS( 2, 3, 0, 1 )

		xorps		xmm5, xmm5
		movhps		xmm5, incs
		movlps		xmm7, incs+8
		movhps		xmm7, incs+16
		addps		xmm3, xmm5
		addps		xmm4, xmm7
		shufps		xmm5, xmm7, R_SHUFFLEPS( 2, 3, 0, 1 )
		movaps		xmm6, xmm7
		shufps		xmm6, xmm5, R_SHUFFLEPS( 2, 3, 0, 1 )
		addps		xmm5, xmm5
		addps		xmm6, xmm6
		addps		xmm7, xmm7

	loop24:
		add			edi, 6*16

		movaps		xmm0, [esi+eax]

		movaps		xmm1, xmm0
		shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 0, 0 )
		mulps		xmm1, xmm2
		addps		xmm1, [edi-6*16]
		addps		xmm2, xmm5
		movaps		[edi-6*16], xmm1

		movaps		xmm1, xmm0
		shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 0, 1, 1 )
		mulps		xmm1, xmm3
		addps		xmm1, [edi-5*16]
		addps		xmm3, xmm6
		movaps		[edi-5*16], xmm1

		movaps		xmm1, xmm0
		shufps		xmm1, xmm1, R_SHUFFLEPS( 1, 1, 1, 1 )
		mulps		xmm1, xmm4
		addps		xmm1, [edi-4*16]
		addps		xmm4, xmm7
		movaps		[edi-4*16], xmm1

		movaps		xmm1, xmm0
		shufps		xmm1, xmm1, R_SHUFFLEPS( 2, 2, 2, 2 )
		mulps		xmm1, xmm2
		addps		xmm1, [edi-3*16]
		addps		xmm2, xmm5
		movaps		[edi-3*16], xmm1

		movaps		xmm1, xmm0
		shufps		xmm1, xmm1, R_SHUFFLEPS( 2, 2, 3, 3 )
		mulps		xmm1, xmm3
		addps		xmm1, [edi-2*16]
		addps		xmm3, xmm6
		movaps		[edi-2*16], xmm1

		shufps		xmm0, xmm0, R_SHUFFLEPS( 3, 3, 3, 3 )
		mulps		xmm0, xmm4
		addps		xmm0, [edi-1*16]
		addps		xmm4, xmm7
		movaps		[edi-1*16], xmm0

		add			eax, 4*4

		jl			loop24
	}

#else

	int i;
	float sL0, sL1, sL2, sL3, sL4, sL5, sL6, sL7, sL8, sL9, sL10, sL11;
	float incL0, incL1, incL2, incL3, incL4, incL5;

	SMF_ASSERT( numSamples == MIXBUFFER_SAMPLES );

	incL0 = ( currentV[0] - lastV[0] ) / MIXBUFFER_SAMPLES;
	incL1 = ( currentV[1] - lastV[1] ) / MIXBUFFER_SAMPLES;
	incL2 = ( currentV[2] - lastV[2] ) / MIXBUFFER_SAMPLES;
	incL3 = ( currentV[3] - lastV[3] ) / MIXBUFFER_SAMPLES;
	incL4 = ( currentV[4] - lastV[4] ) / MIXBUFFER_SAMPLES;
	incL5 = ( currentV[5] - lastV[5] ) / MIXBUFFER_SAMPLES;

	sL0  = lastV[0];
	sL1  = lastV[1];
	sL2  = lastV[2];
	sL3  = lastV[3];
	sL4  = lastV[4];
	sL5  = lastV[5];

	sL6  = lastV[0] + incL0;
	sL7  = lastV[1] + incL1;
	sL8  = lastV[2] + incL2;
	sL9  = lastV[3] + incL3;
	sL10 = lastV[4] + incL4;
	sL11 = lastV[5] + incL5;

	incL0 *= 2;
	incL1 *= 2;
	incL2 *= 2;
	incL3 *= 2;
	incL4 *= 2;
	incL5 *= 2;

	for( i = 0; i <= MIXBUFFER_SAMPLES - 2; i += 2 ) {
		mixBuffer[i*6+ 0] += samples[i+0] * sL0;
		mixBuffer[i*6+ 1] += samples[i+0] * sL1;
		mixBuffer[i*6+ 2] += samples[i+0] * sL2;
		mixBuffer[i*6+ 3] += samples[i+0] * sL3;

		mixBuffer[i*6+ 4] += samples[i+0] * sL4;
		mixBuffer[i*6+ 5] += samples[i+0] * sL5;
		mixBuffer[i*6+ 6] += samples[i+1] * sL6;
		mixBuffer[i*6+ 7] += samples[i+1] * sL7;

		mixBuffer[i*6+ 8] += samples[i+1] * sL8;
		mixBuffer[i*6+ 9] += samples[i+1] * sL9;
		mixBuffer[i*6+10] += samples[i+1] * sL10;
		mixBuffer[i*6+11] += samples[i+1] * sL11;

		sL0  += incL0;
		sL1  += incL1;
		sL2  += incL2;
		sL3  += incL3;

		sL4  += incL4;
		sL5  += incL5;
		sL6  += incL0;
		sL7  += incL1;

		sL8  += incL2;
		sL9  += incL3;
		sL10 += incL4;
		sL11 += incL5;
	}

#endif
}

/*
============
CSIMD_SSE:: mixSoundSixSpeakerStereo
============
*/
void VPCALL CSIMD_SSE:: mixSoundSixSpeakerStereo( float *mixBuffer, const float *samples, const int numSamples, const float lastV[6], const float currentV[6] ) {
#if 1

	ALIGN16( float incs[6] );

	SMF_ASSERT( numSamples == MIXBUFFER_SAMPLES );
	SMF_ASSERT( SPEAKER_RIGHT == 1 );
	SMF_ASSERT( SPEAKER_BACKRIGHT == 5 );

	incs[0] = ( currentV[0] - lastV[0] ) / MIXBUFFER_SAMPLES;
	incs[1] = ( currentV[1] - lastV[1] ) / MIXBUFFER_SAMPLES;
	incs[2] = ( currentV[2] - lastV[2] ) / MIXBUFFER_SAMPLES;
	incs[3] = ( currentV[3] - lastV[3] ) / MIXBUFFER_SAMPLES;
	incs[4] = ( currentV[4] - lastV[4] ) / MIXBUFFER_SAMPLES;
	incs[5] = ( currentV[5] - lastV[5] ) / MIXBUFFER_SAMPLES;

	__asm {
		mov			eax, MIXBUFFER_SAMPLES
		mov			edi, mixBuffer
		mov			esi, samples
		shl			eax, 3
		add			esi, eax
		neg			eax

		mov			ecx, lastV
		movlps		xmm2, [ecx+ 0]
		movhps		xmm2, [ecx+ 8]
		movlps		xmm3, [ecx+16]
		movaps		xmm4, xmm2
		shufps		xmm3, xmm2, R_SHUFFLEPS( 0, 1, 0, 1 )
		shufps		xmm4, xmm3, R_SHUFFLEPS( 2, 3, 0, 1 )

		xorps		xmm5, xmm5
		movhps		xmm5, incs
		movlps		xmm7, incs+ 8
		movhps		xmm7, incs+16
		addps		xmm3, xmm5
		addps		xmm4, xmm7
		shufps		xmm5, xmm7, R_SHUFFLEPS( 2, 3, 0, 1 )
		movaps		xmm6, xmm7
		shufps		xmm6, xmm5, R_SHUFFLEPS( 2, 3, 0, 1 )
		addps		xmm5, xmm5
		addps		xmm6, xmm6
		addps		xmm7, xmm7

	loop12:
		add			edi, 3*16

		movaps		xmm0, [esi+eax+0]

		movaps		xmm1, xmm0
		shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 1, 0, 0 )
		mulps		xmm1, xmm2
		addps		xmm1, [edi-3*16]
		addps		xmm2, xmm5
		movaps		[edi-3*16], xmm1

		movaps		xmm1, xmm0
		shufps		xmm1, xmm1, R_SHUFFLEPS( 0, 1, 2, 3 )
		mulps		xmm1, xmm3
		addps		xmm1, [edi-2*16]
		addps		xmm3, xmm6
		movaps		[edi-2*16], xmm1

		add			eax, 4*4

		shufps		xmm0, xmm0, R_SHUFFLEPS( 2, 2, 2, 3 )
		mulps		xmm0, xmm4
		addps		xmm0, [edi-1*16]
		addps		xmm4, xmm7
		movaps		[edi-1*16], xmm0

		jl			loop12

		emms
	}

#else

	int i;
	float sL0, sL1, sL2, sL3, sL4, sL5, sL6, sL7, sL8, sL9, sL10, sL11;
	float incL0, incL1, incL2, incL3, incL4, incL5;

	SMF_ASSERT( numSamples == MIXBUFFER_SAMPLES );
	SMF_ASSERT( SPEAKER_RIGHT == 1 );
	SMF_ASSERT( SPEAKER_BACKRIGHT == 5 );

	incL0 = ( currentV[0] - lastV[0] ) / MIXBUFFER_SAMPLES;
	incL1 = ( currentV[1] - lastV[1] ) / MIXBUFFER_SAMPLES;
	incL2 = ( currentV[2] - lastV[2] ) / MIXBUFFER_SAMPLES;
	incL3 = ( currentV[3] - lastV[3] ) / MIXBUFFER_SAMPLES;
	incL4 = ( currentV[4] - lastV[4] ) / MIXBUFFER_SAMPLES;
	incL5 = ( currentV[5] - lastV[5] ) / MIXBUFFER_SAMPLES;

	sL0  = lastV[0];
	sL1  = lastV[1];
	sL2  = lastV[2];
	sL3  = lastV[3];
	sL4  = lastV[4];
	sL5  = lastV[5];

	sL6  = lastV[0] + incL0;
	sL7  = lastV[1] + incL1;
	sL8  = lastV[2] + incL2;
	sL9  = lastV[3] + incL3;
	sL10 = lastV[4] + incL4;
	sL11 = lastV[5] + incL5;

	incL0 *= 2;
	incL1 *= 2;
	incL2 *= 2;
	incL3 *= 2;
	incL4 *= 2;
	incL5 *= 2;

	for( i = 0; i <= MIXBUFFER_SAMPLES - 2; i += 2 ) {
		mixBuffer[i*6+ 0] += samples[i*2+0+0] * sL0;
		mixBuffer[i*6+ 1] += samples[i*2+0+1] * sL1;
		mixBuffer[i*6+ 2] += samples[i*2+0+0] * sL2;
		mixBuffer[i*6+ 3] += samples[i*2+0+0] * sL3;

		mixBuffer[i*6+ 4] += samples[i*2+0+0] * sL4;
		mixBuffer[i*6+ 5] += samples[i*2+0+1] * sL5;
		mixBuffer[i*6+ 6] += samples[i*2+2+0] * sL6;
		mixBuffer[i*6+ 7] += samples[i*2+2+1] * sL7;

		mixBuffer[i*6+ 8] += samples[i*2+2+0] * sL8;
		mixBuffer[i*6+ 9] += samples[i*2+2+0] * sL9;
		mixBuffer[i*6+10] += samples[i*2+2+0] * sL10;
		mixBuffer[i*6+11] += samples[i*2+2+1] * sL11;

		sL0  += incL0;
		sL1  += incL1;
		sL2  += incL2;
		sL3  += incL3;

		sL4  += incL4;
		sL5  += incL5;
		sL6  += incL0;
		sL7  += incL1;

		sL8  += incL2;
		sL9  += incL3;
		sL10 += incL4;
		sL11 += incL5;
	}

#endif
}

/*
============
CSIMD_SSE:: mixedSoundToSamples
============
*/
void VPCALL CSIMD_SSE:: mixedSoundToSamples( short *samples, const float *mixBuffer, const int numSamples ) {
#if 1

	SMF_ASSERT( ( numSamples % MIXBUFFER_SAMPLES ) == 0 );

	__asm {

		mov			eax, numSamples
		mov			edi, mixBuffer
		mov			esi, samples
		shl			eax, 2
		add			edi, eax
		neg			eax

	loop16:

		movaps		xmm0, [edi+eax+0*16]
		movaps		xmm2, [edi+eax+1*16]
		movaps		xmm4, [edi+eax+2*16]
		movaps		xmm6, [edi+eax+3*16]

		add			esi, 4*4*2

		movhlps		xmm1, xmm0
		movhlps		xmm3, xmm2
		movhlps		xmm5, xmm4
		movhlps		xmm7, xmm6

		prefetchnta	[edi+eax+64]

		cvtps2pi	mm0, xmm0
		cvtps2pi	mm2, xmm2
		cvtps2pi	mm4, xmm4
		cvtps2pi	mm6, xmm6

		prefetchnta	[edi+eax+128]

		cvtps2pi	mm1, xmm1
		cvtps2pi	mm3, xmm3
		cvtps2pi	mm5, xmm5
		cvtps2pi	mm7, xmm7

		add			eax, 4*16

		packssdw	mm0, mm1
		packssdw	mm2, mm3
		packssdw	mm4, mm5
		packssdw	mm6, mm7

		movq		[esi-4*4*2], mm0
		movq		[esi-3*4*2], mm2
		movq		[esi-2*4*2], mm4
		movq		[esi-1*4*2], mm6

		jl			loop16

		emms
	}

#else

	for ( int i = 0; i < numSamples; i++ ) {
		if ( mixBuffer[i] <= -32768.0f ) {
			samples[i] = -32768;
		} else if ( mixBuffer[i] >= 32767.0f ) {
			samples[i] = 32767;
		} else {
			samples[i] = (short) mixBuffer[i];
		}
	}

#endif
}
//======================  CVec3D============================
#define  HIPREC
void VPCALL CSIMD_SSE::vector3D_Sum(CVec3D* pOut, const CVec3D* pIn)
{
	/* SSE/SSE2/SSE3 Implementation */
	 __asm
        {       
                mov eax, pIn                       // Load pointers into CPU regs
                mov ebx, pOut
				mov ecx, pOut
                movups xmm0, [eax]                 // Move unaligned vectors to SSE regs
                movups xmm1, [ebx]
				andps xmm0, _SIMDx86_float_SSE_NO_W_MASK //zera o elemento w
                andps xmm1, _SIMDx86_float_SSE_NO_W_MASK  //remove w component
				addps xmm0, xmm1                   // add vector elements
                movlpd	[ecx+ 0], xmm0  //retira os tres bytes
				shufps xmm0,xmm0, R_SHUFFLEPS( 2, 0, 0, 0 ) //shuffle colocar nos bits 0-31 o 3o byte
				movss	[ecx+ 8], xmm0
        }
		
}

void VPCALL CSIMD_SSE::vector3D_SumOf(CVec3D* pOut, const CVec3D* pIn1, const CVec3D* pIn2)
{
	/* SSE/SSE2/SSE3 Implementation */
		 __asm
        {       
                mov eax, pIn1                       // Load pointers into CPU regs
                mov ebx, pIn2
				mov ecx, pOut
                movups xmm0, [eax]                 // Move unaligned vectors to SSE regs
                movups xmm1, [ebx]
				andps xmm0, _SIMDx86_float_SSE_NO_W_MASK //zera o elemento w
                andps xmm1, _SIMDx86_float_SSE_NO_W_MASK  //remove w component
				addps xmm0, xmm1                   // add vector elements
                movlpd	[ecx+ 0], xmm0  //retira os tres bytes
				shufps xmm0,xmm0, R_SHUFFLEPS( 2, 0, 0, 0 ) //shuffle colocar nos bits 0-31 o 3o byte
				movss	[ecx+ 8], xmm0
        }	
}


//Subtracts pIn from pOut and stores the result in pOut. POut = pIn-pOut
void VPCALL CSIMD_SSE::vector3D_Diff(CVec3D* pOut, CVec3D* pIn)
{
	/* SSE/SSE2/SSE3 Implementation */
		 __asm
        {       
                mov eax, pOut                      // Load pointers into CPU regs
                mov ebx, pIn
				mov ecx, pOut
                movups xmm0, [eax]                 // Move unaligned vectors to SSE regs
                movups xmm1, [ebx]
				andps xmm0, _SIMDx86_float_SSE_NO_W_MASK //zera o elemento w
                andps xmm1, _SIMDx86_float_SSE_NO_W_MASK  //remove w component
				subps xmm0, xmm1                   // xmm0 = xmm0 - xmm1
                movlpd	[ecx+ 0], xmm0  //retira os tres bytes
				shufps xmm0,xmm0, R_SHUFFLEPS( 2, 0, 0, 0 ) //shuffle colocar nos bits 0-31 o 3o byte
				movss	[ecx+ 8], xmm0
        }
}


//Subtracts pIn2 from pIn1 and stores the result in pOut.  pout = pIn1-Pin2
void VPCALL CSIMD_SSE::vector3D_DiffOf(CVec3D* pOut, const CVec3D* pIn1, const CVec3D* pIn2)
{

	__asm
        {       
                mov eax, pIn1                       // Load pointers into CPU regs
                mov ebx, pIn2
				mov ecx, pOut
                movups xmm0, [eax]                 // Move unaligned vectors to SSE regs
                movups xmm1, [ebx]
				andps xmm0, _SIMDx86_float_SSE_NO_W_MASK //zera o elemento w
                andps xmm1, _SIMDx86_float_SSE_NO_W_MASK  //remove w component
				subps xmm0, xmm1                   // xmm0 = xmm0 - xmm1
                movlpd	[ecx+ 0], xmm0  //retira os tres bytes
				shufps xmm0,xmm0, R_SHUFFLEPS( 2, 0, 0, 0 ) //shuffle colocar nos bits 0-31 o 3o byte
				movss	[ecx+ 8], xmm0
        }
}
//Scales the components of pOut by scalar and stores the result in pOut.
void VPCALL CSIMD_SSE::vector3D_Scale(CVec3D* pOut, float scalar)
{
	float *test=&scalar;
	__asm
        {
			mov eax, pOut                       // Load pointers into CPU regs
			mov ebx, test
			mov ecx, pOut
            movss xmm0, [ebx]
			movups xmm1, [eax]
			andps xmm1, _SIMDx86_float_SSE_NO_W_MASK  //remove w component
			shufps xmm0,xmm0,R_SHUFFLEPS( 0, 0, 0, 0 )	//repete o valor da scalar para os demais bytes de xmm0
			mulps xmm0,xmm1
		    movlpd	[ecx+ 0], xmm0  //retira os tres bytes
			shufps xmm0,xmm0, R_SHUFFLEPS( 2, 0, 0, 0 ) //shuffle colocar nos bits 0-31 o 3o byte
			movss	[ecx+ 8], xmm0
		}


}
//Scales the components of pIn by scalar and stores the result in pOut.
void VPCALL CSIMD_SSE::vector3D_ScaleOf(CVec3D* pOut, const CVec3D* pIn, float scalar)
{
	float *test=&scalar;
	__asm
        {
			mov eax, pIn                       // Load pointers into CPU regs
			mov ebx, test
			mov ecx, pOut
            movss xmm0, [ebx]
			movups xmm1, [eax]
			andps xmm1, _SIMDx86_float_SSE_NO_W_MASK  //remove w component /* xmm0 = x | y | z | 0.0f */
			shufps xmm0,xmm0,R_SHUFFLEPS( 0, 0, 0, 0 )	//repete o valor da scalar para os demais bytes de xmm0
			mulps xmm0,xmm1
		    movlpd	[ecx+ 0], xmm0  //retira os tres bytes
			shufps xmm0,xmm0, R_SHUFFLEPS( 2, 0, 0, 0 ) //shuffle colocar nos bits 0-31 o 3o byte
			movss	[ecx+ 8], xmm0
		}



}
//Calculates the dot product (known as the inner product) of pSrc1 and pSrc2, that represents the cosine of the angle between them, if both are normalized.
float  VPCALL CSIMD_SSE::vector3D_Dot(const CVec3D* pSrc1, const CVec3D* pSrc2)
{

		/* SSE/SSE2 Implementation */
		float dummy;
		float *min = &dummy;
		__asm
        {
			mov eax, pSrc1                       // Load pointers into CPU regs and multiply
			mov ebx, pSrc2
			movups xmm0, [eax]
			movups xmm1, [ebx]
			andps xmm0, _SIMDx86_float_SSE_NO_W_MASK  //remove w component /* xmm0 = x | y | z | 0.0f */
			andps xmm1, _SIMDx86_float_SSE_NO_W_MASK  //remove w component /* xmm1 = x | y | z | 0.0f */
			mulps xmm1, xmm0
			/* Shift data around in the registers (I loves me some haddps right now!!) */
			movhlps xmm2, xmm1    /* xmm2 = ?   | ?   | 0   |  z's */
		    addss   xmm2, xmm1    /* xmm2 = ?   | ?   | 0   | z's + x's */
			shufps  xmm1, xmm1, R_SHUFFLEPS( 1,1,1,1 )  /*0x55=01010101 / xmm1 = y's | y's | y's | y's */
			addss   xmm2, xmm1 /* xmm2 = ?   | ?   | ?   | x's+y's+z's */
			//coloca 0-31 de xmm2 em dummy
			mov			esi, min
			movss		[esi], xmm2
		
		}
		
		return dummy;
		

}

//Returns the squared length of pSrc1. This can be useful when doing sphere-sphere collision tests where the exact distance is not needed and a squared one can be used, for example.
//It is signifigantly faster than using the square root version SIMDx86Vector_Length().

float  VPCALL CSIMD_SSE::vector3D_LengthSq(const CVec3D* pSrc1)
{

		float dummy;

		
		/* SSE/SSE2 Implementation */
		float *min = &dummy;
		__asm
        {
			mov eax, pSrc1                       // Load pointers into CPU regs and multiply
			movups xmm0, [eax]
			andps xmm0, _SIMDx86_float_SSE_NO_W_MASK  //remove w component /* xmm0 = x | y | z | 0.0f */
			mulps xmm0, xmm0    /* xmm0 = xmm0^2 */

			/* Shift data around in the registers (I loves me some haddps right now!!) */
			movhlps xmm1, xmm0		/* xmm1 = ?   | ?   | 0   |  z's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | 0   | z's + x's */
			shufps  xmm0, xmm0, R_SHUFFLEPS( 1,1,1,1 )  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			
		    //coloca 0-31 de xmm1 em dummy
			mov			esi, min
			movss		[esi], xmm1


		}

		/* USE_SSE == 1 ||  USE_SSE == 2 */

		return dummy;

}
//Returns the length (magnitude) of pSrc1. normalized vectors (e.g. the output of SIMDx86Vector_Normalize() ) are of unit length, or 1.0)
float  VPCALL CSIMD_SSE::vector3D_Length(const CVec3D* pSrc1)
{

		float dummy;
		/* SSE/SSE2 Implementation */
		float *min = &dummy;
		__asm
        {
			mov eax, pSrc1                       // Load pointers into CPU regs and multiply
			movups xmm0, [eax]
			andps xmm0, _SIMDx86_float_SSE_NO_W_MASK  //remove w component /* xmm0 = x | y | z | 0.0f */
			mulps xmm0, xmm0    /* xmm0 = xmm0^2 */

				/* Shift data around in the registers (I loves me some haddps right now!!) */
			movhlps xmm1, xmm0		/* xmm1 = ?   | ?   | 0   |  z's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | 0   | z's + x's */
			shufps  xmm0, xmm0, R_SHUFFLEPS( 1,1,1,1 )  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			#ifdef HIPREC
			/* Full square root */
			sqrtss xmm1, xmm1
			#else
			/* rcp( rsqrt(value) ) (This may be very inaccurate) */
			rsqrtss xmm1, xmm1
			rcpss xmm1, xmm1
			#endif

		    //coloca 0-31 de xmm1 em dummy
			mov			esi, min
			movss		[esi], xmm1


		}
		
		


		return dummy;
		
	
}


void VPCALL CSIMD_SSE::vector3D_Cross(CVec3D* pLeft, const CVec3D* pRight)
{
	
	__asm
        {
			mov eax, pLeft 
			mov ebx, pRight // Load pointers into CPU regs and multiply
			mov ecx, pLeft
			movups xmm0, [eax]  /* xmm0 = pLeft */
			movups xmm1, [ebx]  /* xmm1 = pRight */
			andps xmm0, _SIMDx86_float_SSE_NO_W_MASK  //remove w component /* xmm0 = x | y | z | 0.0f */
			andps xmm1, _SIMDx86_float_SSE_NO_W_MASK  //remove w component /* xmm1 = x | y | z | 0.0f */
			movaps xmm2,xmm0  /* xmm2 = pLeft */
			movaps xmm3,xmm1  /* xmm3 = pRight */
			/*
			Shuffles
			w | x | z | y  == 1100 1001 = 0xC9
			w | y | x | z  == 1101 0010 = 0xD2
			*/
			shufps xmm0,xmm0, SHUFFLEPS(3,0,2,1)  /* 0xC9=11001001 */
			shufps xmm1,xmm1, SHUFFLEPS(3,1,0,2)  /* 0xD2=11010010 */
			shufps xmm2,xmm2, SHUFFLEPS(3,1,0,2)  /* 0xD2=11010010 */
			shufps xmm3,xmm3, SHUFFLEPS(3,0,2,1)  /* 0xC9=11001001 */
			/* multiply columns 1&2 and 3&4 */
			mulps xmm1, xmm0
			mulps xmm3, xmm2

			/* Subtract products to get the cross product! */
			subps xmm1, xmm3

			/* Store */
	        movlpd	[ecx+ 0], xmm1  //retira os tres bytes
			shufps xmm1,xmm1, R_SHUFFLEPS( 2, 0, 0, 0 ) //shuffle colocar nos bits 0-31 o 3o byte
			movss	[ecx+ 8], xmm1

		}

}

void VPCALL CSIMD_SSE::vector3D_CrossOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight)
{

	__asm
        {
			mov eax, pLeft 
			mov ebx, pRight // Load pointers into CPU regs and multiply
			mov ecx, pOut
			movups xmm0, [eax]  /* xmm0 = pLeft */
			movups xmm1, [ebx]  /* xmm1 = pRight */
			andps xmm0, _SIMDx86_float_SSE_NO_W_MASK  //remove w component /* xmm0 = x | y | z | 0.0f */
			andps xmm1, _SIMDx86_float_SSE_NO_W_MASK  //remove w component /* xmm1 = x | y | z | 0.0f */
			movaps xmm2,xmm0  /* xmm2 = pLeft */
			movaps xmm3,xmm1  /* xmm3 = pRight */
			/*
			Shuffles
			w | x | z | y  == 1100 1001 = 0xC9
			w | y | x | z  == 1101 0010 = 0xD2
			*/
			shufps xmm0,xmm0, SHUFFLEPS(3,0,2,1)  /* 0xC9=11001001 */
			shufps xmm1,xmm1, SHUFFLEPS(3,1,0,2)  /* 0xD2=11010010 */
			shufps xmm2,xmm2, SHUFFLEPS(3,1,0,2)  /* 0xD2=11010010 */
			shufps xmm3,xmm3, SHUFFLEPS(3,0,2,1)  /* 0xC9=11001001 */
			/* multiply columns 1&2 and 3&4 */
			mulps xmm1, xmm0
			mulps xmm3, xmm2

			/* Subtract products to get the cross product! */
			subps xmm1, xmm3

			/* Store */
	        movlpd	[ecx+ 0], xmm1  //retira os tres bytes
			shufps xmm1,xmm1, R_SHUFFLEPS( 2, 0, 0, 0 ) //shuffle colocar nos bits 0-31 o 3o byte
			movss	[ecx+ 8], xmm1

		}


}

void VPCALL CSIMD_SSE::vector3D_Normalize(CVec3D* pVec)
{

		
		/* SSE/SSE2 Implementation */
			__asm
        {
			mov eax, pVec 
			mov ecx, pVec
			movups xmm0, [eax]  /* xmm0 = pVec */
			andps xmm0, _SIMDx86_float_SSE_NO_W_MASK  //remove w component /* xmm0 = x | y | z | 0.0f */
			movaps xmm7, xmm0			/* Save for the division by length */
		    mulps xmm0, xmm0    /* xmm0 = xmm0^2 */

				/* Shift data around in the registers (I loves me some haddps right now!!) */
			movhlps xmm1, xmm0		/* xmm1 = ?   | ?   | 0   |  z's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | 0   | z's + x's */
			shufps  xmm0, xmm0, R_SHUFFLEPS( 1,1,1,1 )  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			
		
			#ifdef HIPREC

			/* Divide by magnitude */
			sqrtss xmm1, xmm1			/* xmm1 = ? | ? | ? | mag(pVec) */
			shufps xmm1, xmm1,	SHUFFLEPS( 0,0,0,0 ) /*$0x00=00000000,  xmm1 = mag(pVec) | mag(pVec) | mag(pVec) | mag(pVec) */
			divps xmm7, xmm1			/* xmm1 = ? | z/mag | y/mag | x/mag */

			#else

			/* multiply by reciprocal magnitude */
			rsqrtss xmm1 ,xmm1			/* xmm1 = ? | ? | ? | rcp(mag) */
			shufps xmm1, xmm1, SHUFFLEPS( 0,0,0,0 )	        /*$0x00=00000000,  xmm3 = rcp(mag) | rcp(mag) | rcp(mag) | rcp(mag) */
			mulps xmm7, xmm1			/* xmm7 = ? | z/mag | y/mag | x/mag */

			#endif

			/* Store */
	        movlpd	[ecx+ 0], xmm7  //retira os tres bytes
			shufps xmm7,xmm7, R_SHUFFLEPS( 2, 0, 0, 0 ) //shuffle colocar nos bits 0-31 o 3o byte
			movss	[ecx+ 8], xmm7

		}


}

void VPCALL CSIMD_SSE::vector3D_NormalizeOf(CVec3D* pOut, const CVec3D* pVec)
{

		/* SSE/SSE2 Implementation */
			__asm
        {
			mov eax, pVec 
			mov ecx, pOut
			movups xmm0, [eax]  /* xmm0 = pVec */
			andps xmm0, _SIMDx86_float_SSE_NO_W_MASK  //remove w component /* xmm0 = x | y | z | 0.0f */
			movaps xmm7, xmm0			/* Save for the division by length */
		    mulps xmm0, xmm0    /* xmm0 = xmm0^2 */

				/* Shift data around in the registers (I loves me some haddps right now!!) */
			movhlps xmm1, xmm0		/* xmm1 = ?   | ?   | 0   |  z's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | 0   | z's + x's */
			shufps  xmm0, xmm0, R_SHUFFLEPS( 1,1,1,1 )  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			
		
			#ifdef HIPREC

			/* Divide by magnitude */
			sqrtss xmm1, xmm1			/* xmm1 = ? | ? | ? | mag(pVec) */
			shufps xmm1, xmm1,	SHUFFLEPS( 0,0,0,0 ) /*$0x00=00000000,  xmm1 = mag(pVec) | mag(pVec) | mag(pVec) | mag(pVec) */
			divps xmm7, xmm1			/* xmm1 = ? | z/mag | y/mag | x/mag */

			#else

			/* multiply by reciprocal magnitude */
			rsqrtss xmm1 ,xmm1			/* xmm1 = ? | ? | ? | rcp(mag) */
			shufps xmm1, xmm1, SHUFFLEPS( 0,0,0,0 )	        /*$0x00=00000000,  xmm3 = rcp(mag) | rcp(mag) | rcp(mag) | rcp(mag) */
			mulps xmm7, xmm1			/* xmm7 = ? | z/mag | y/mag | x/mag */

			#endif

			/* Store */
	        movlpd	[ecx+ 0], xmm7  //retira os tres bytes
			shufps xmm7,xmm7, R_SHUFFLEPS( 2, 0, 0, 0 ) //shuffle colocar nos bits 0-31 o 3o byte
			movss	[ecx+ 8], xmm7

		}

	/* USE_SSE == 1 || USE_SSE == 2 */


}

float  VPCALL CSIMD_SSE::vector3D_Distance(const CVec3D* pVec1, const CVec3D* pVec2)
{
	/* TODO: Optimize me completely */
	CVec3D diff;
	vector3D_DiffOf(&diff, pVec1, pVec2);
	return vector3D_Length(&diff);
}

/*
	TODO:
	Optimize distance() -- It currently uses already implemented functions
*/

void VPCALL CSIMD_SSE::vector3D_AlignedSum(CVec3D* pOut, const CVec3D* pIn)
{


		/* SSE/SSE2/SSE3 Implementation */
	 __asm
        {       
                mov eax, pIn                       // Load pointers into CPU regs
                mov ebx, pOut
				mov ecx, pOut
                movaps xmm0, [eax]                 // Move unaligned vectors to SSE regs
                movaps xmm1, [ebx]
				andps xmm0, _SIMDx86_float_SSE_NO_W_MASK //zera o elemento w
                andps xmm1, _SIMDx86_float_SSE_NO_W_MASK  //remove w component
				addps xmm0, xmm1                   // add vector elements
                movlpd	[ecx+ 0], xmm0  //retira os tres bytes
				shufps xmm0,xmm0, R_SHUFFLEPS( 2, 0, 0, 0 ) //shuffle colocar nos bits 0-31 o 3o byte
				movss	[ecx+ 8], xmm0
        }

}

void VPCALL CSIMD_SSE::vector3D_AlignedSumOf(CVec3D* pOut, const CVec3D* pIn1, const CVec3D* pIn2)
{
		/* SSE/SSE2/SSE3 Implementation */
		 __asm
        {       
                mov eax, pIn1                       // Load pointers into CPU regs
                mov ebx, pIn2
				mov ecx, pOut
                movaps xmm0, [eax]                 // Move unaligned vectors to SSE regs
                movaps xmm1, [ebx]
				andps xmm0, _SIMDx86_float_SSE_NO_W_MASK //zera o elemento w
                andps xmm1, _SIMDx86_float_SSE_NO_W_MASK  //remove w component
				addps xmm0, xmm1                   // add vector elements
                movlpd	[ecx+ 0], xmm0  //retira os tres bytes
				shufps xmm0,xmm0, R_SHUFFLEPS( 2, 0, 0, 0 ) //shuffle colocar nos bits 0-31 o 3o byte
				movss	[ecx+ 8], xmm0
        }


}


void VPCALL CSIMD_SSE::vector3D_AlignedDiff(CVec3D* pOut, CVec3D* pIn)
{

		/* SSE/SSE2/SSE3 Implementation */
		 __asm
        {       
                mov eax, pOut                      // Load pointers into CPU regs
                mov ebx, pIn
				mov ecx, pOut
                movaps xmm0, [eax]                 // Move unaligned vectors to SSE regs
                movaps xmm1, [ebx]
				andps xmm0, _SIMDx86_float_SSE_NO_W_MASK //zera o elemento w
                andps xmm1, _SIMDx86_float_SSE_NO_W_MASK  //remove w component
				subps xmm0, xmm1                   // xmm0 = xmm0 - xmm1
                movlpd	[ecx+ 0], xmm0  //retira os tres bytes
				shufps xmm0,xmm0, R_SHUFFLEPS( 2, 0, 0, 0 ) //shuffle colocar nos bits 0-31 o 3o byte
				movss	[ecx+ 8], xmm0
        }

}

void VPCALL CSIMD_SSE::vector3D_AlignedDiffOf(CVec3D* pOut, const CVec3D* pIn1, const CVec3D* pIn2)
{

		__asm
        {       
                mov eax, pIn1                       // Load pointers into CPU regs
                mov ebx, pIn2
				mov ecx, pOut
                movaps xmm0, [eax]                 // Move unaligned vectors to SSE regs
                movaps xmm1, [ebx]
				andps xmm0, _SIMDx86_float_SSE_NO_W_MASK //zera o elemento w
                andps xmm1, _SIMDx86_float_SSE_NO_W_MASK  //remove w component
				subps xmm0, xmm1                   // xmm0 = xmm0 - xmm1
                movlpd	[ecx+ 0], xmm0  //retira os tres bytes
				shufps xmm0,xmm0, R_SHUFFLEPS( 2, 0, 0, 0 ) //shuffle colocar nos bits 0-31 o 3o byte
				movss	[ecx+ 8], xmm0
        }



}

void VPCALL CSIMD_SSE::vector3D_AlignedScale(CVec3D* pOut, float scalar)
{

	float *test=&scalar;
	__asm
        {
			mov eax, pOut                       // Load pointers into CPU regs
			mov ebx, test
			mov ecx, pOut
            movss xmm0, [ebx]
			movaps xmm1, [eax]
			andps xmm1, _SIMDx86_float_SSE_NO_W_MASK  //remove w component
			shufps xmm0,xmm0,R_SHUFFLEPS( 0, 0, 0, 0 )	//repete o valor da scalar para os demais bytes de xmm0
			mulps xmm0,xmm1
		    movlpd	[ecx+ 0], xmm0  //retira os tres bytes
			shufps xmm0,xmm0, R_SHUFFLEPS( 2, 0, 0, 0 ) //shuffle colocar nos bits 0-31 o 3o byte
			movss	[ecx+ 8], xmm0
		}



}

void VPCALL CSIMD_SSE::vector3D_AlignedScaleOf(CVec3D* pOut, const CVec3D* pIn, float scalar)
{

	float *test=&scalar;
	__asm
        {
			mov eax, pIn                       // Load pointers into CPU regs
			mov ebx, test
			mov ecx, pOut
            movss xmm0, [ebx]
			movaps xmm1, [eax]
			andps xmm1, _SIMDx86_float_SSE_NO_W_MASK  //remove w component /* xmm0 = x | y | z | 0.0f */
			shufps xmm0,xmm0,R_SHUFFLEPS( 0, 0, 0, 0 )	//repete o valor da scalar para os demais bytes de xmm0
			mulps xmm0,xmm1
		    movlpd	[ecx+ 0], xmm0  //retira os tres bytes
			shufps xmm0,xmm0, R_SHUFFLEPS( 2, 0, 0, 0 ) //shuffle colocar nos bits 0-31 o 3o byte
			movss	[ecx+ 8], xmm0
		}


}

float  VPCALL CSIMD_SSE::vector3D_AlignedDot(const CVec3D* pSrc1, const CVec3D* pSrc2)
{

		/* SSE/SSE2 Implementation */
		float dummy;
		float *min = &dummy;
		__asm
        {
			mov eax, pSrc1                       // Load pointers into CPU regs and multiply
			mov ebx, pSrc2
			movaps xmm0, [eax]
			movaps xmm1, [ebx]
			andps xmm0, _SIMDx86_float_SSE_NO_W_MASK  //remove w component /* xmm0 = x | y | z | 0.0f */
			andps xmm1, _SIMDx86_float_SSE_NO_W_MASK  //remove w component /* xmm1 = x | y | z | 0.0f */
			mulps xmm1, xmm0
			/* Shift data around in the registers (I loves me some haddps right now!!) */
			movhlps xmm2, xmm1    /* xmm2 = ?   | ?   | 0   |  z's */
		    addss   xmm2, xmm1    /* xmm2 = ?   | ?   | 0   | z's + x's */
			shufps  xmm1, xmm1, R_SHUFFLEPS( 1,1,1,1 )  /*0x55=01010101 / xmm1 = y's | y's | y's | y's */
			addss   xmm2, xmm1 /* xmm2 = ?   | ?   | ?   | x's+y's+z's */
			//coloca 0-31 de xmm2 em dummy
			mov			esi, min
			movss		[esi], xmm2
		
		}


}


float  VPCALL CSIMD_SSE::vector3D_AlignedLengthSq(const CVec3D* pSrc1)
{


		float dummy;

		
		/* SSE/SSE2 Implementation */
		float *min = &dummy;
		__asm
        {
			mov eax, pSrc1                       // Load pointers into CPU regs and multiply
			movaps xmm0, [eax]
			andps xmm0, _SIMDx86_float_SSE_NO_W_MASK  //remove w component /* xmm0 = x | y | z | 0.0f */
			mulps xmm0, xmm0    /* xmm0 = xmm0^2 */

			/* Shift data around in the registers (I loves me some haddps right now!!) */
			movhlps xmm1, xmm0		/* xmm1 = ?   | ?   | 0   |  z's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | 0   | z's + x's */
			shufps  xmm0, xmm0, R_SHUFFLEPS( 1,1,1,1 )  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			
		    //coloca 0-31 de xmm1 em dummy
			mov			esi, min
			movss		[esi], xmm1


		}

		/* USE_SSE == 1 ||  USE_SSE == 2 */

		return dummy;

}

float  VPCALL CSIMD_SSE::vector3D_AlignedLength(const CVec3D* pSrc1)
{


		float dummy;
		/* SSE/SSE2 Implementation */
		float *min = &dummy;
		__asm
        {
			mov eax, pSrc1                       // Load pointers into CPU regs and multiply
			movaps xmm0, [eax]
			andps xmm0, _SIMDx86_float_SSE_NO_W_MASK  //remove w component /* xmm0 = x | y | z | 0.0f */
			mulps xmm0, xmm0    /* xmm0 = xmm0^2 */

				/* Shift data around in the registers (I loves me some haddps right now!!) */
			movhlps xmm1, xmm0		/* xmm1 = ?   | ?   | 0   |  z's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | 0   | z's + x's */
			shufps  xmm0, xmm0, R_SHUFFLEPS( 1,1,1,1 )  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			#ifdef HIPREC
			/* Full square root */
			sqrtss xmm1, xmm1
			#else
			/* rcp( rsqrt(value) ) (This may be very inaccurate) */
			rsqrtss xmm1, xmm1
			rcpss xmm1, xmm1
			#endif

		    //coloca 0-31 de xmm1 em dummy
			mov			esi, min
			movss		[esi], xmm1


		}
		
		


		return dummy;


}

void VPCALL CSIMD_SSE::vector3D_AlignedCross(CVec3D* pLeft, const CVec3D* pRight)
{

	__asm
        {
			mov eax, pLeft 
			mov ebx, pRight // Load pointers into CPU regs and multiply
			mov ecx, pLeft
			movaps xmm0, [eax]  /* xmm0 = pLeft */
			movaps xmm1, [ebx]  /* xmm1 = pRight */
			andps xmm0, _SIMDx86_float_SSE_NO_W_MASK  //remove w component /* xmm0 = x | y | z | 0.0f */
			andps xmm1, _SIMDx86_float_SSE_NO_W_MASK  //remove w component /* xmm1 = x | y | z | 0.0f */
			movaps xmm2,xmm0  /* xmm2 = pLeft */
			movaps xmm3,xmm1  /* xmm3 = pRight */
			/*
			Shuffles
			w | x | z | y  == 1100 1001 = 0xC9
			w | y | x | z  == 1101 0010 = 0xD2
			*/
			shufps xmm0,xmm0, SHUFFLEPS(3,0,2,1)  /* 0xC9=11001001 */
			shufps xmm1,xmm1, SHUFFLEPS(3,1,0,2)  /* 0xD2=11010010 */
			shufps xmm2,xmm2, SHUFFLEPS(3,1,0,2)  /* 0xD2=11010010 */
			shufps xmm3,xmm3, SHUFFLEPS(3,0,2,1)  /* 0xC9=11001001 */
			/* multiply columns 1&2 and 3&4 */
			mulps xmm1, xmm0
			mulps xmm3, xmm2

			/* Subtract products to get the cross product! */
			subps xmm1, xmm3

			/* Store */
	        movlpd	[ecx+ 0], xmm1  //retira os tres bytes
			shufps xmm1,xmm1, R_SHUFFLEPS( 2, 0, 0, 0 ) //shuffle colocar nos bits 0-31 o 3o byte
			movss	[ecx+ 8], xmm1

		}


}

void VPCALL CSIMD_SSE::vector3D_AlignedCrossOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight)
{

		__asm
        {
			mov eax, pLeft 
			mov ebx, pRight // Load pointers into CPU regs and multiply
			mov ecx, pOut
			movaps xmm0, [eax]  /* xmm0 = pLeft */
			movaps xmm1, [ebx]  /* xmm1 = pRight */
			andps xmm0, _SIMDx86_float_SSE_NO_W_MASK  //remove w component /* xmm0 = x | y | z | 0.0f */
			andps xmm1, _SIMDx86_float_SSE_NO_W_MASK  //remove w component /* xmm1 = x | y | z | 0.0f */
			movaps xmm2,xmm0  /* xmm2 = pLeft */
			movaps xmm3,xmm1  /* xmm3 = pRight */
			/*
			Shuffles
			w | x | z | y  == 1100 1001 = 0xC9
			w | y | x | z  == 1101 0010 = 0xD2
			*/
			shufps xmm0,xmm0, SHUFFLEPS(3,0,2,1)  /* 0xC9=11001001 */
			shufps xmm1,xmm1, SHUFFLEPS(3,1,0,2)  /* 0xD2=11010010 */
			shufps xmm2,xmm2, SHUFFLEPS(3,1,0,2)  /* 0xD2=11010010 */
			shufps xmm3,xmm3, SHUFFLEPS(3,0,2,1)  /* 0xC9=11001001 */
			/* multiply columns 1&2 and 3&4 */
			mulps xmm1, xmm0
			mulps xmm3, xmm2

			/* Subtract products to get the cross product! */
			subps xmm1, xmm3

			/* Store */
	        movlpd	[ecx+ 0], xmm1  //retira os tres bytes
			shufps xmm1,xmm1, R_SHUFFLEPS( 2, 0, 0, 0 ) //shuffle colocar nos bits 0-31 o 3o byte
			movss	[ecx+ 8], xmm1

		}


}
void VPCALL CSIMD_SSE::vector3D_AlignedNormalize(CVec3D* pVec)
{

		/* SSE/SSE2 Implementation */
			__asm
        {
			mov eax, pVec 
			mov ecx, pVec
			movaps xmm0, [eax]  /* xmm0 = pVec */
			andps xmm0, _SIMDx86_float_SSE_NO_W_MASK  //remove w component /* xmm0 = x | y | z | 0.0f */
			movaps xmm7, xmm0			/* Save for the division by length */
		    mulps xmm0, xmm0    /* xmm0 = xmm0^2 */

				/* Shift data around in the registers (I loves me some haddps right now!!) */
			movhlps xmm1, xmm0		/* xmm1 = ?   | ?   | 0   |  z's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | 0   | z's + x's */
			shufps  xmm0, xmm0, R_SHUFFLEPS( 1,1,1,1 )  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			
		
			#ifdef HIPREC

			/* Divide by magnitude */
			sqrtss xmm1, xmm1			/* xmm1 = ? | ? | ? | mag(pVec) */
			shufps xmm1, xmm1,	SHUFFLEPS( 0,0,0,0 ) /*$0x00=00000000,  xmm1 = mag(pVec) | mag(pVec) | mag(pVec) | mag(pVec) */
			divps xmm7, xmm1			/* xmm1 = ? | z/mag | y/mag | x/mag */

			#else

			/* multiply by reciprocal magnitude */
			rsqrtss xmm1 ,xmm1			/* xmm1 = ? | ? | ? | rcp(mag) */
			shufps xmm1, xmm1, SHUFFLEPS( 0,0,0,0 )	        /*$0x00=00000000,  xmm3 = rcp(mag) | rcp(mag) | rcp(mag) | rcp(mag) */
			mulps xmm7, xmm1			/* xmm7 = ? | z/mag | y/mag | x/mag */

			#endif

			/* Store */
	        movlpd	[ecx+ 0], xmm7  //retira os tres bytes
			shufps xmm7,xmm7, R_SHUFFLEPS( 2, 0, 0, 0 ) //shuffle colocar nos bits 0-31 o 3o byte
			movss	[ecx+ 8], xmm7

		}


}
void VPCALL CSIMD_SSE::vector3D_AlignedNormalizeOf(CVec3D* pOut, const CVec3D* pVec)
{

		/* SSE/SSE2 Implementation */
			__asm
        {
			mov eax, pVec 
			mov ecx, pOut
			movups xmm0, [eax]  /* xmm0 = pVec */
			andps xmm0, _SIMDx86_float_SSE_NO_W_MASK  //remove w component /* xmm0 = x | y | z | 0.0f */
			movaps xmm7, xmm0			/* Save for the division by length */
		    mulps xmm0, xmm0    /* xmm0 = xmm0^2 */

				/* Shift data around in the registers (I loves me some haddps right now!!) */
			movhlps xmm1, xmm0		/* xmm1 = ?   | ?   | 0   |  z's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | 0   | z's + x's */
			shufps  xmm0, xmm0, R_SHUFFLEPS( 1,1,1,1 )  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			
		
			#ifdef HIPREC

			/* Divide by magnitude */
			sqrtss xmm1, xmm1			/* xmm1 = ? | ? | ? | mag(pVec) */
			shufps xmm1, xmm1,	SHUFFLEPS( 0,0,0,0 ) /*$0x00=00000000,  xmm1 = mag(pVec) | mag(pVec) | mag(pVec) | mag(pVec) */
			divps xmm7, xmm1			/* xmm1 = ? | z/mag | y/mag | x/mag */

			#else

			/* multiply by reciprocal magnitude */
			rsqrtss xmm1 ,xmm1			/* xmm1 = ? | ? | ? | rcp(mag) */
			shufps xmm1, xmm1, SHUFFLEPS( 0,0,0,0 )	        /*$0x00=00000000,  xmm3 = rcp(mag) | rcp(mag) | rcp(mag) | rcp(mag) */
			mulps xmm7, xmm1			/* xmm7 = ? | z/mag | y/mag | x/mag */

			#endif

			/* Store */
	        movlpd	[ecx+ 0], xmm7  //retira os tres bytes
			shufps xmm7,xmm7, R_SHUFFLEPS( 2, 0, 0, 0 ) //shuffle colocar nos bits 0-31 o 3o byte
			movss	[ecx+ 8], xmm7

		}



}

float  VPCALL CSIMD_SSE::vector3D_AlignedDistance(const CVec3D* pVec1, const CVec3D* pVec2)
{
#if 0 /* Not so good just yet... */
	#if defined(USE_SSE)

	#if USE_SSE == 3
		float dummy;
		__asm{
			movaps (%1), xmm0
			subps (%2), xmm0
			"andps _SIMDx86_float_SSE_NO_W_MASK, xmm0
			"haddps xmm0, xmm0
			"haddps xmm0, xmm0
			"subl $4, esp
			movss xmm0, (esp)
			"flds (esp)
			"addl $4, esp
			: "=t" (dummy)
			: "r" (pVec1), "r" (pVec2)
			);
		return dummy;
	#else
	#endif
	#endif

#endif

	/* TODO: Optimize me completely */
	CVec3D diff;
	vector3D_DiffOf(&diff, pVec1, pVec2);
	return vector3D_Length(&diff);
	
}



//======================  CVec4D============================

void VPCALL CSIMD_SSE::vector4D_Sum(CVec4D* pOut, const CVec4D* pIn)
{
	/* SSE/SSE2/SSE3 Implementation */
	 __asm
        {       
                mov eax, pIn                       // Load pointers into CPU regs
                mov ebx, pOut
				mov ecx, pOut
                movups xmm0, [eax]                 // Move unaligned vectors to SSE regs
                movups xmm1, [ebx]
				addps xmm0, xmm1                   // add vector elements
                movups	[ecx], xmm0  //retira os tres bytes
				
        }
		
}

void VPCALL CSIMD_SSE::vector4D_SumOf(CVec4D* pOut, const CVec4D* pIn1, const CVec4D* pIn2)
{
	/* SSE/SSE2/SSE3 Implementation */
		 __asm
        {       
                mov eax, pIn1                       // Load pointers into CPU regs
                mov ebx, pIn2
				mov ecx, pOut
                movups xmm0, [eax]                 // Move unaligned vectors to SSE regs
                movups xmm1, [ebx]
				addps xmm0, xmm1                   // add vector elements
                movups	[ecx], xmm0  //retira os tres bytes
				
        }	
}


//Subtracts pIn from pOut and stores the result in pOut. POut = pIn-pOut
void VPCALL CSIMD_SSE::vector4D_Diff(CVec4D* pOut, CVec4D* pIn)
{
	/* SSE/SSE2/SSE3 Implementation */
		 __asm
        {       
                mov eax, pOut                      // Load pointers into CPU regs
                mov ebx, pIn
				mov ecx, pOut
                movups xmm0, [eax]                 // Move unaligned vectors to SSE regs
                movups xmm1, [ebx]
				subps xmm0, xmm1                   // xmm0 = xmm0 - xmm1
                movups	[ecx], xmm0  //retira os tres bytes
				
        }
}


//Subtracts pIn2 from pIn1 and stores the result in pOut.  pout = pIn1-Pin2
void VPCALL CSIMD_SSE::vector4D_DiffOf(CVec4D* pOut, const CVec4D* pIn1, const CVec4D* pIn2)
{

	__asm
        {       
                mov eax, pIn1                       // Load pointers into CPU regs
                mov ebx, pIn2
				mov ecx, pOut
                movups xmm0, [eax]                 // Move unaligned vectors to SSE regs
                movups xmm1, [ebx]
				subps xmm0, xmm1                   // xmm0 = xmm0 - xmm1
                movups	[ecx+ 0], xmm0  //retira os tres bytes
				
        }
}
//Scales the components of pOut by scalar and stores the result in pOut.
void VPCALL CSIMD_SSE::vector4D_Scale(CVec4D* pOut, float scalar)
{
	float *test=&scalar;
	__asm
        {
			mov eax, pOut                       // Load pointers into CPU regs
			mov ebx, test
			mov ecx, pOut
            movss xmm0, [ebx]
			movups xmm1, [eax]
			shufps xmm0,xmm0,R_SHUFFLEPS( 0, 0, 0, 0 )	//repete o valor da scalar para os demais bytes de xmm0
			mulps xmm0,xmm1
		    movups	[ecx+ 0], xmm0  //retira os tres bytes
			
		}


}
//Scales the components of pIn by scalar and stores the result in pOut.
void VPCALL CSIMD_SSE::vector4D_ScaleOf(CVec4D* pOut, const CVec4D* pIn, float scalar)
{
	float *test=&scalar;
	__asm
        {
			mov eax, pIn                       // Load pointers into CPU regs
			mov ebx, test
			mov ecx, pOut
            movss xmm0, [ebx]
			movups xmm1, [eax]
			shufps xmm0,xmm0,R_SHUFFLEPS( 0, 0, 0, 0 )	//repete o valor da scalar para os demais bytes de xmm0
			mulps xmm0,xmm1
		    movups	[ecx+ 0], xmm0  //retira os tres bytes
			
		}



}
//Calculates the dot product (known as the inner product) of pSrc1 and pSrc2, that represents the cosine of the angle between them, if both are normalized.

float  VPCALL CSIMD_SSE::vector4D_Dot(const CVec4D* pSrc1, const CVec4D* pSrc2)
{

	/* SSE/SSE2 Implementation */	
	float dummy;
	float *min = &dummy;
		__asm
        {
			mov eax, pSrc1                       // Load pointers into CPU regs and multiply
			mov ebx, pSrc2
			movups xmm0, [eax]
			movups xmm1, [ebx]
			mulps xmm1, xmm0
			movaps xmm2,xmm1
			shufps xmm2, xmm2,  SHUFFLEPS( 0,1,2,3 ) /*0x1B= 00011011 / xmm2 = x | y | z | w */
			
			addps  xmm2, xmm1			/* xmm2 = w+x | y+z | y+z | w+x */
			movss  xmm3, xmm2			/* xmm3 = ??? | ??? | ??? | w+x */
			shufps xmm1, xmm1, SHUFFLEPS( 0,0,0,1 )  /* 0x01= 00000001 / xmm3 = ??? | ??? | ??? | w+x */
			addss  xmm2, xmm3			/* xmm2 = ??? | ??? | ??? | dot4 */
	
			
			//coloca 0-31 de xmm2 em dummy
			mov			esi, min
			movss		[esi], xmm2
		
		}

	return dummy;
	


}
//Returns the squared length of pSrc1. This can be useful when doing sphere-sphere collision tests where the exact distance is not needed and a squared one can be used, for example.
//It is signifigantly faster than using the square root version SIMDx86Vector_Length().

float  VPCALL CSIMD_SSE::vector4D_LengthSq(const CVec4D* pSrc1)
{

		float dummy;
		/* SSE/SSE2 Implementation */
		float *min = &dummy;
		__asm
        {
			mov eax, pSrc1                       // Load pointers into CPU regs and multiply
			movups xmm0, [eax]  /* xmm0 = x | y | z | w */
			mulps xmm0, xmm0    /* xmm0 = xmm0^2 */

			/* Shift data around in the registers (I loves me some haddps right now!!) */
			movaps  xmm7,xmm0
			shufps  xmm7,xmm7, R_SHUFFLEPS( 3,0,0,0 )  /* xmm7 = ?   | ?   | ?  |  w's *///colocar nos bis 0-31 de xmm7 os bits 96-127 de xmm7
			
			movhlps xmm1, xmm0		/* xmm1 = ?   | ?   | w's  |  z's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | w's+y's   | z's + x's */
			shufps  xmm0, xmm0, R_SHUFFLEPS( 1,1,1,1 )  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			addss xmm1, xmm7		/* xmm1 = ?   | ?   | ?   | x's+y's+z's+w's */
			
		    //coloca 0-31 de xmm1 em dummy
			mov			esi, min
			movss		[esi], xmm1


		}

		/* USE_SSE == 1 ||  USE_SSE == 2 */

		return dummy;

}
//Returns the length (magnitude) of pSrc1. normalized vectors (e.g. the output of SIMDx86Vector_Normalize() ) are of unit length, or 1.0)
float  VPCALL CSIMD_SSE::vector4D_Length(const CVec4D* pSrc1)
{

		float dummy;
		/* SSE/SSE2 Implementation */
		float *min = &dummy;
		__asm
        {
			mov eax, pSrc1                       // Load pointers into CPU regs and multiply
			movups xmm0, [eax]
			mulps xmm0, xmm0    /* xmm0 = xmm0^2 */

			/* Shift data around in the registers (I loves me some haddps right now!!) */
			movaps  xmm7,xmm0
			shufps  xmm7,xmm7, R_SHUFFLEPS( 3,0,0,0 )  /* xmm7 = ?   | ?   | ?  |  w's *///colocar nos bis 0-31 de xmm7 os bits 96-127 de xmm7
			
			movhlps xmm1, xmm0		/* xmm1 = ?   | ?   | w's  |  z's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | w's+y's   | z's + x's */
			shufps  xmm0, xmm0, R_SHUFFLEPS( 1,1,1,1 )  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			addss xmm1, xmm7		/* xmm1 = ?   | ?   | ?   | x's+y's+z's+w's */
			#ifdef HIPREC
			/* Full square root */
			sqrtss xmm1, xmm1
			#else
			/* rcp( rsqrt(value) ) (This may be very inaccurate) */
			rsqrtss xmm1, xmm1
			rcpss xmm1, xmm1
			#endif

		    //coloca 0-31 de xmm1 em dummy
			mov			esi, min
			movss		[esi], xmm1


		}
		
		


		return dummy;
		
	
}




void VPCALL CSIMD_SSE::vector4D_Normalize(CVec4D* pVec)
{

		
		/* SSE/SSE2 Implementation */
			__asm
        {
			mov eax, pVec 
			mov ecx, pVec
			movups xmm0, [eax]  /* xmm0 = pVec */
			movaps xmm7, xmm0			/* Save for the division by length */
		    mulps xmm0, xmm0    /* xmm0 = xmm0^2 */

			/* Shift data around in the registers (I loves me some haddps right now!!) */
			movaps  xmm2,xmm0
			shufps  xmm2,xmm2, R_SHUFFLEPS( 3,0,0,0 )  /* xmm2 = ?   | ?   | ?  |  w's *///colocar nos bis 0-31 de xmm7 os bits 96-127 de xmm7
			
			movhlps xmm1, xmm0		/* xmm1 = ?   | ?   | w's  |  z's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | w's+y's   | z's + x's */
			shufps  xmm0, xmm0, R_SHUFFLEPS( 1,1,1,1 )  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			addss xmm1, xmm2		/* xmm1 = ?   | ?   | ?   | x's+y's+z's+w's */

		
			#ifdef HIPREC

			/* Divide by magnitude */
			sqrtss xmm1, xmm1			/* xmm1 = ? | ? | ? | mag(pVec) */
			shufps xmm1, xmm1,	SHUFFLEPS( 0,0,0,0 ) /*$0x00=00000000,  xmm1 = mag(pVec) | mag(pVec) | mag(pVec) | mag(pVec) */
			divps xmm7, xmm1			/* xmm1 = w/mag | z/mag | y/mag | x/mag */

			#else

			/* multiply by reciprocal magnitude */
			rsqrtss xmm1 ,xmm1			/* xmm1 = ? | ? | ? | rcp(mag) */
			shufps xmm1, xmm1, SHUFFLEPS( 0,0,0,0 )	        /*$0x00=00000000,  xmm3 = rcp(mag) | rcp(mag) | rcp(mag) | rcp(mag) */
			mulps xmm7, xmm1			/* xmm7 = w/mag | z/mag | y/mag | x/mag */

			#endif

			/* Store */
	        movups	[ecx+ 0], xmm7  //retira os tres bytes

		}


}

void VPCALL CSIMD_SSE::vector4D_NormalizeOf(CVec4D* pOut, const CVec4D* pVec)
{

		/* SSE/SSE2 Implementation */
			__asm
        {
			mov eax, pVec 
			mov ecx, pOut
			movups xmm0, [eax]  /* xmm0 = pVec */
			movaps xmm7, xmm0			/* Save for the division by length */
		    mulps xmm0, xmm0    /* xmm0 = xmm0^2 */

			/* Shift data around in the registers (I loves me some haddps right now!!) */
			movaps  xmm2,xmm0
			shufps  xmm2,xmm2, R_SHUFFLEPS( 3,0,0,0 )  /* xmm7 = ?   | ?   | ?  |  w's *///colocar nos bis 0-31 de xmm7 os bits 96-127 de xmm7
			
			movhlps xmm1, xmm0		/* xmm1 = ?   | ?   | w's  |  z's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | w's+y's   | z's + x's */
			shufps  xmm0, xmm0, R_SHUFFLEPS( 1,1,1,1 )  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			addss xmm1, xmm2		/* xmm1 = ?   | ?   | ?   | x's+y's+z's+w's */

		
			#ifdef HIPREC

			/* Divide by magnitude */
			sqrtss xmm1, xmm1			/* xmm1 = ? | ? | ? | mag(pVec) */
			shufps xmm1, xmm1,	SHUFFLEPS( 0,0,0,0 ) /*$0x00=00000000,  xmm1 = mag(pVec) | mag(pVec) | mag(pVec) | mag(pVec) */
			divps xmm7, xmm1			/* xmm1 = w/mag | z/mag | y/mag | x/mag */

			#else

			/* multiply by reciprocal magnitude */
			rsqrtss xmm1 ,xmm1			/* xmm1 = ? | ? | ? | rcp(mag) */
			shufps xmm1, xmm1, SHUFFLEPS( 0,0,0,0 )	        /*$0x00=00000000,  xmm3 = rcp(mag) | rcp(mag) | rcp(mag) | rcp(mag) */
			mulps xmm7, xmm1			/* xmm7 = w/mag | z/mag | y/mag | x/mag */

			#endif

			/* Store */
	        movups	[ecx+ 0], xmm7  //retira os tres bytes

		}

	/* USE_SSE == 1 || USE_SSE == 2 */


}

float  VPCALL CSIMD_SSE::vector4D_Distance(const CVec4D* pVec1, const CVec4D* pVec2)
{
	/* TODO: Optimize me completely */
	CVec4D diff;
	vector4D_DiffOf(&diff, pVec1, pVec2);
	return vector4D_Length(&diff);
}








void VPCALL CSIMD_SSE::vector4D_AlignedSum(CVec4D* pOut, const CVec4D* pIn)
{
	/* SSE/SSE2/SSE3 Implementation */
	 __asm
        {       
                mov eax, pIn                       // Load pointers into CPU regs
                mov ebx, pOut
				mov ecx, pOut
                movaps xmm0, [eax]                 // Move unaligned vectors to SSE regs
                movaps xmm1, [ebx]
				addps xmm0, xmm1                   // add vector elements
                movups	[ecx], xmm0  //retira os tres bytes
				
        }
		
}

void VPCALL CSIMD_SSE::vector4D_AlignedSumOf(CVec4D* pOut, const CVec4D* pIn1, const CVec4D* pIn2)
{
	/* SSE/SSE2/SSE3 Implementation */
		 __asm
        {       
                mov eax, pIn1                       // Load pointers into CPU regs
                mov ebx, pIn2
				mov ecx, pOut
                movaps xmm0, [eax]                 // Move unaligned vectors to SSE regs
                movaps xmm1, [ebx]
				addps xmm0, xmm1                   // add vector elements
                movups	[ecx], xmm0  //retira os tres bytes
				
        }	
}


//Subtracts pIn from pOut and stores the result in pOut. POut = pIn-pOut
void VPCALL CSIMD_SSE::vector4D_AlignedDiff(CVec4D* pOut, CVec4D* pIn)
{
	/* SSE/SSE2/SSE3 Implementation */
		 __asm
        {       
                mov eax, pOut                      // Load pointers into CPU regs
                mov ebx, pIn
				mov ecx, pOut
                movaps xmm0, [eax]                 // Move unaligned vectors to SSE regs
                movaps xmm1, [ebx]
				subps xmm0, xmm1                   // xmm0 = xmm0 - xmm1
                movups	[ecx], xmm0  //retira os tres bytes
				
        }
}


//Subtracts pIn2 from pIn1 and stores the result in pOut.  pout = pIn1-Pin2
void VPCALL CSIMD_SSE::vector4D_AlignedDiffOf(CVec4D* pOut, const CVec4D* pIn1, const CVec4D* pIn2)
{

	__asm
        {       
                mov eax, pIn1                       // Load pointers into CPU regs
                mov ebx, pIn2
				mov ecx, pOut
                movaps xmm0, [eax]                 // Move unaligned vectors to SSE regs
                movaps xmm1, [ebx]
				subps xmm0, xmm1                   // xmm0 = xmm0 - xmm1
                movups	[ecx+ 0], xmm0  //retira os tres bytes
				
        }
}
//Scales the components of pOut by scalar and stores the result in pOut.
void VPCALL CSIMD_SSE::vector4D_AlignedScale(CVec4D* pOut, float scalar)
{
	float *test=&scalar;
	__asm
        {
			mov eax, pOut                       // Load pointers into CPU regs
			mov ebx, test
			mov ecx, pOut
            movss xmm0, [ebx]
			movaps xmm1, [eax]
			shufps xmm0,xmm0,R_SHUFFLEPS( 0, 0, 0, 0 )	//repete o valor da scalar para os demais bytes de xmm0
			mulps xmm0,xmm1
		    movups	[ecx+ 0], xmm0  //retira os tres bytes
			
		}


}
//Scales the components of pIn by scalar and stores the result in pOut.
void VPCALL CSIMD_SSE::vector4D_AlignedScaleOf(CVec4D* pOut, const CVec4D* pIn, float scalar)
{
	float *test=&scalar;
	__asm
        {
			mov eax, pIn                       // Load pointers into CPU regs
			mov ebx, test
			mov ecx, pOut
            movss xmm0, [ebx]
			movaps xmm1, [eax]
			shufps xmm0,xmm0,R_SHUFFLEPS( 0, 0, 0, 0 )	//repete o valor da scalar para os demais bytes de xmm0
			mulps xmm0,xmm1
		    movups	[ecx+ 0], xmm0  //retira os tres bytes
			
		}



}
//Calculates the dot product (known as the inner product) of pSrc1 and pSrc2, that represents the cosine of the angle between them, if both are normalized.

float  VPCALL CSIMD_SSE::vector4D_AlignedDot(const CVec4D* pSrc1, const CVec4D* pSrc2)
{

	/* SSE/SSE2 Implementation */	
	float dummy;
	float *min = &dummy;
		__asm
        {
			mov eax, pSrc1                       // Load pointers into CPU regs and multiply
			mov ebx, pSrc2
			movaps xmm0, [eax]
			movaps xmm1, [ebx]
			mulps xmm1, xmm0
			movaps xmm2,xmm1
			shufps xmm2, xmm2,  SHUFFLEPS( 0,1,2,3 ) /*0x1B= 00011011 / xmm2 = x | y | z | w */
			
			addps  xmm2, xmm1			/* xmm2 = w+x | y+z | y+z | w+x */
			movss  xmm3, xmm2			/* xmm3 = ??? | ??? | ??? | w+x */
			shufps xmm1, xmm1, SHUFFLEPS( 0,0,0,1 )  /* 0x01= 00000001 / xmm3 = ??? | ??? | ??? | w+x */
			addss  xmm2, xmm3			/* xmm2 = ??? | ??? | ??? | dot4 */
	
			
			//coloca 0-31 de xmm2 em dummy
			mov			esi, min
			movss		[esi], xmm2
		
		}

	return dummy;
	


}
//Returns the squared length of pSrc1. This can be useful when doing sphere-sphere collision tests where the exact distance is not needed and a squared one can be used, for example.
//It is signifigantly faster than using the square root version SIMDx86Vector_Length().

float  VPCALL CSIMD_SSE::vector4D_AlignedLengthSq(const CVec4D* pSrc1)
{

		float dummy;
		/* SSE/SSE2 Implementation */
		float *min = &dummy;
		__asm
        {
			mov eax, pSrc1                       // Load pointers into CPU regs and multiply
			movaps xmm0, [eax]  /* xmm0 = x | y | z | w */
			mulps xmm0, xmm0    /* xmm0 = xmm0^2 */

			/* Shift data around in the registers (I loves me some haddps right now!!) */
			movaps  xmm7,xmm0
			shufps  xmm7,xmm7, R_SHUFFLEPS( 3,0,0,0 )  /* xmm7 = ?   | ?   | ?  |  w's *///colocar nos bis 0-31 de xmm7 os bits 96-127 de xmm7
			
			movhlps xmm1, xmm0		/* xmm1 = ?   | ?   | w's  |  z's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | w's+y's   | z's + x's */
			shufps  xmm0, xmm0, R_SHUFFLEPS( 1,1,1,1 )  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			addss xmm1, xmm7		/* xmm1 = ?   | ?   | ?   | x's+y's+z's+w's */
			
		    //coloca 0-31 de xmm1 em dummy
			mov			esi, min
			movss		[esi], xmm1


		}

		/* USE_SSE == 1 ||  USE_SSE == 2 */

		return dummy;

}
//Returns the length (magnitude) of pSrc1. normalized vectors (e.g. the output of SIMDx86Vector_Normalize() ) are of unit length, or 1.0)
float  VPCALL CSIMD_SSE::vector4D_AlignedLength(const CVec4D* pSrc1)
{

		float dummy;
		/* SSE/SSE2 Implementation */
		float *min = &dummy;
		__asm
        {
			mov eax, pSrc1                       // Load pointers into CPU regs and multiply
			movaps xmm0, [eax]
			mulps xmm0, xmm0    /* xmm0 = xmm0^2 */

			/* Shift data around in the registers (I loves me some haddps right now!!) */
			movaps  xmm7,xmm0
			shufps  xmm7,xmm7, R_SHUFFLEPS( 3,0,0,0 )  /* xmm7 = ?   | ?   | ?  |  w's *///colocar nos bis 0-31 de xmm7 os bits 96-127 de xmm7
			
			movhlps xmm1, xmm0		/* xmm1 = ?   | ?   | w's  |  z's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | w's+y's   | z's + x's */
			shufps  xmm0, xmm0, R_SHUFFLEPS( 1,1,1,1 )  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			addss xmm1, xmm7		/* xmm1 = ?   | ?   | ?   | x's+y's+z's+w's */
			#ifdef HIPREC
			/* Full square root */
			sqrtss xmm1, xmm1
			#else
			/* rcp( rsqrt(value) ) (This may be very inaccurate) */
			rsqrtss xmm1, xmm1
			rcpss xmm1, xmm1
			#endif

		    //coloca 0-31 de xmm1 em dummy
			mov			esi, min
			movss		[esi], xmm1


		}
		
		


		return dummy;
		
	
}




void VPCALL CSIMD_SSE::vector4D_AlignedNormalize(CVec4D* pVec)
{

		
		/* SSE/SSE2 Implementation */
			__asm
        {
			mov eax, pVec 
			mov ecx, pVec
			movaps xmm0, [eax]  /* xmm0 = pVec */
			movaps xmm7, xmm0			/* Save for the division by length */
		    mulps xmm0, xmm0    /* xmm0 = xmm0^2 */

			/* Shift data around in the registers (I loves me some haddps right now!!) */
			movaps  xmm2,xmm0
			shufps  xmm2,xmm2, R_SHUFFLEPS( 3,0,0,0 )  /* xmm2 = ?   | ?   | ?  |  w's *///colocar nos bis 0-31 de xmm7 os bits 96-127 de xmm7
			
			movhlps xmm1, xmm0		/* xmm1 = ?   | ?   | w's  |  z's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | w's+y's   | z's + x's */
			shufps  xmm0, xmm0, R_SHUFFLEPS( 1,1,1,1 )  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			addss xmm1, xmm2		/* xmm1 = ?   | ?   | ?   | x's+y's+z's+w's */

		
			#ifdef HIPREC

			/* Divide by magnitude */
			sqrtss xmm1, xmm1			/* xmm1 = ? | ? | ? | mag(pVec) */
			shufps xmm1, xmm1,	SHUFFLEPS( 0,0,0,0 ) /*$0x00=00000000,  xmm1 = mag(pVec) | mag(pVec) | mag(pVec) | mag(pVec) */
			divps xmm7, xmm1			/* xmm1 = w/mag | z/mag | y/mag | x/mag */

			#else

			/* multiply by reciprocal magnitude */
			rsqrtss xmm1 ,xmm1			/* xmm1 = ? | ? | ? | rcp(mag) */
			shufps xmm1, xmm1, SHUFFLEPS( 0,0,0,0 )	        /*$0x00=00000000,  xmm3 = rcp(mag) | rcp(mag) | rcp(mag) | rcp(mag) */
			mulps xmm7, xmm1			/* xmm7 = w/mag | z/mag | y/mag | x/mag */

			#endif

			/* Store */
	        movups	[ecx+ 0], xmm7  //retira os tres bytes

		}


}

void VPCALL CSIMD_SSE::vector4D_AlignedNormalizeOf(CVec4D* pOut, const CVec4D* pVec)
{

		/* SSE/SSE2 Implementation */
			__asm
        {
			mov eax, pVec 
			mov ecx, pOut
			movaps xmm0, [eax]  /* xmm0 = pVec */
			movaps xmm7, xmm0			/* Save for the division by length */
		    mulps xmm0, xmm0    /* xmm0 = xmm0^2 */

			/* Shift data around in the registers (I loves me some haddps right now!!) */
			movaps  xmm2,xmm0
			shufps  xmm2,xmm2, R_SHUFFLEPS( 3,0,0,0 )  /* xmm7 = ?   | ?   | ?  |  w's *///colocar nos bis 0-31 de xmm7 os bits 96-127 de xmm7
			
			movhlps xmm1, xmm0		/* xmm1 = ?   | ?   | w's  |  z's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | w's+y's   | z's + x's */
			shufps  xmm0, xmm0, R_SHUFFLEPS( 1,1,1,1 )  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			addss xmm1, xmm0		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			addss xmm1, xmm2		/* xmm1 = ?   | ?   | ?   | x's+y's+z's+w's */

		
			#ifdef HIPREC

			/* Divide by magnitude */
			sqrtss xmm1, xmm1			/* xmm1 = ? | ? | ? | mag(pVec) */
			shufps xmm1, xmm1,	SHUFFLEPS( 0,0,0,0 ) /*$0x00=00000000,  xmm1 = mag(pVec) | mag(pVec) | mag(pVec) | mag(pVec) */
			divps xmm7, xmm1			/* xmm1 = w/mag | z/mag | y/mag | x/mag */

			#else

			/* multiply by reciprocal magnitude */
			rsqrtss xmm1 ,xmm1			/* xmm1 = ? | ? | ? | rcp(mag) */
			shufps xmm1, xmm1, SHUFFLEPS( 0,0,0,0 )	        /*$0x00=00000000,  xmm3 = rcp(mag) | rcp(mag) | rcp(mag) | rcp(mag) */
			mulps xmm7, xmm1			/* xmm7 = w/mag | z/mag | y/mag | x/mag */

			#endif

			/* Store */
	        movups	[ecx+ 0], xmm7  //retira os tres bytes

		}

	/* USE_SSE == 1 || USE_SSE == 2 */


}

float  VPCALL CSIMD_SSE::vector4D_AlignedDistance(const CVec4D* pVec1, const CVec4D* pVec2)
{
	/* TODO: Optimize me completely */
	CVec4D diff;
	vector4D_DiffOf(&diff, pVec1, pVec2);
	return vector4D_Length(&diff);
}

#define LOADREGS0( pOut ) \
	__asm mov ecx, pOut    

#define LOADREGS( pOut , pIn ) \
	__asm mov ecx, pOut    \
	__asm mov eax, pIn

#define LOADREGS2( pOut , pIn1, pIn2 ) \
	__asm mov ecx, pOut    \
	__asm mov eax, pIn1		\
	__asm mov ebx, pIn2		\


//========CMat4D==============================
void VPCALL CSIMD_SSE::mat4D_Sum(CMat4D* Out, const CMat4D* In)
{
		
			LOADREGS( Out , In )
		__asm
        {	
			movups  xmm0,[ecx]
			movups  xmm1,[ecx+16]
			movups  xmm2,[ecx+32]
			movups  xmm3,[ecx+48]
			movups  xmm4,[eax]
			movups  xmm5,[eax+16]	
			movups  xmm6,[eax+32]	
			movups  xmm7,[eax+48]
			addps	xmm4,xmm0
			addps	xmm5, xmm1
			addps	xmm6,xmm2
			addps	xmm7,xmm3
			movups   [ecx],xmm4  
			movups  [ecx+16],xmm5
			movups  [ecx+32],xmm6
			movups  [ecx+48],xmm7
		}

}


void VPCALL CSIMD_SSE::mat4D_SumOf(CMat4D* Out, const CMat4D* In1, const CMat4D* In2)
{

	LOADREGS2( Out , In1, In2 )
	 __asm {
	movups  xmm0,[eax]
	movups  xmm1,[eax+16]
	movups  xmm2,[eax+32]
	movups  xmm3,[eax+48]
	movups  xmm4,[ebx]
	movups  xmm5,[ebx+16]
	movups  xmm6,[ebx+32]
	movups  xmm7,[ebx+48]
	addps	xmm4,xmm0
	addps	xmm5,xmm1
	addps	xmm6,xmm2
	addps	xmm7,xmm3
	movups  [ecx],xmm4
	movups  [ecx+16],xmm5
	movups  [ecx+32],xmm6
	movups  [ecx+48],xmm7
	}

}

void VPCALL CSIMD_SSE::mat4D_Diff(CMat4D* Out, const CMat4D* In)
{
	
	LOADREGS( Out , In )
	 __asm {
	movups    xmm0,[ecx]
	movups  xmm1,[ecx+16]
	movups  xmm2,[ecx+32]
	movups  xmm3,[ecx+48]
	movups    xmm4,[eax]
	movups  xmm5,[eax+16]
	movups  xmm6,[eax+32]
	movups  xmm7,[eax+48]
	subps  xmm0,xmm4
	subps  xmm1,xmm5
	subps  xmm2,xmm6
	subps  xmm3,xmm7
	movups   [ecx], xmm0
	movups  [ecx+16],xmm1
	movups  [ecx+32],xmm2
	movups  [ecx+48],xmm3
	}
}

void VPCALL CSIMD_SSE::mat4D_DiffOf(CMat4D* Out, const CMat4D* In1, const CMat4D* In2)
{
	LOADREGS2( Out , In1, In2 ) 
	__asm {
	movups  xmm0,[eax]
	movups  xmm1,[eax+16]
	movups  xmm2,[eax+32]
	movups  xmm3,[eax+48]
	movups  xmm4, [ebx]
	movups  xmm5,[ebx+16]
	movups  xmm6,[ebx+32]
	movups	xmm7, [ebx+48]
	subps	xmm0,xmm4
	subps	xmm1,xmm5
	subps	xmm2,xmm6
	subps	xmm3,xmm7
	movups  [ecx],xmm0
	movups  [ecx+16],xmm1
	movups  [ecx+32],xmm2
	movups	[ecx+48],xmm3 
	}
}

void VPCALL CSIMD_SSE::mat4D_Scale(CMat4D* mtx, float scalar)
{
	LOADREGS( mtx, scalar )
	 __asm {
	/* Store scalar in xmm4.x */
	movss  xmm4,eax

	/* Get the matrix into registers */
	movups  xmm0,[ecx]
	movups  xmm1,[ecx+16]
	movups  xmm2,[ecx+32]
	movups  xmm3,[ecx+48]

	/* Broadcast element x to yzw to make a duplicated scalar register */
	shufps  xmm4, xmm4, 0x00

	/* scale the matrix in parallel */
	mulps	xmm0, xmm4
	mulps	xmm1, xmm4
	mulps	xmm2, xmm4
	mulps	xmm3, xmm4

	/* Store results */
	movups  [ecx],xmm0
	movups  [ecx+16],xmm1
	movups  [ecx+32],xmm2
	movups  [ecx+48],xmm3
	}
}

void VPCALL CSIMD_SSE::mat4D_ScaleOf(CMat4D* pOut, const CMat4D* pIn, float scalar)
{
	
	 __asm {
	mov ecx, pOut    
	mov eax, pIn		
	movss	xmm4, scalar
	movups  xmm0, [eax]
	movups  xmm1,[eax+16]
	movups  xmm2,[eax+32]
	movups  xmm3,[eax+48]
	shufps  xmm4, xmm4,0x00
	mulps	xmm0,xmm4
	mulps	xmm1,xmm4
	mulps	xmm2,xmm4
	mulps	xmm3,xmm4
	movups  [ecx],xmm0
	movups  [ecx+16],xmm1
	movups  [ecx+32],xmm2
	movups  [ecx+48],xmm3
	}
}

void VPCALL CSIMD_SSE::mat4D_Multiply(CMat4D* pLeft, const CMat4D* pRight)
{
	LOADREGS( pLeft , pRight ) 
	 
	__asm {
	movups    xmm0,[eax]	/* xmm0 = pRight[0..3] */
	movups  xmm1,[eax+16]	/* xmm1 = pRight[5..7] */
	movups  xmm2,[eax+32]	/* xmm2 = pRight[8..11] */
	movups  xmm3,[eax+48]	/* xmm3 = pRight[12..15] */

	/* Processes 1/2 of the matrix at a time (2x4), unrolled loop */
	movss    xmm4, [ecx]
	movss    xmm6,[ecx+4]
	movss  xmm5,[ecx+16]
	movss   xmm7,[ecx+20]
	shufps xmm4, xmm4, 0x00
	shufps  xmm5, xmm5,0x00
	shufps  xmm6, xmm6,0x00
	shufps  xmm7, xmm7,0x00
	mulps  xmm4,xmm0
	mulps  xmm5,xmm0
	mulps  xmm6,xmm1
	mulps  xmm7,xmm1
	addps  xmm5,xmm7
	addps  xmm4,xmm6


	movss   xmm6,[ecx+8]
	movss	xmm7,[ecx+24]
	shufps  xmm6, xmm6,0x00
	shufps xmm7, xmm7, 0x00
	mulps  xmm6,xmm2
	mulps xmm7,xmm2 
	addps xmm4,xmm6 
	addps xmm5,xmm7 

	movss   xmm6,[ecx+12]
	movss   xmm7,[ecx+28]
	shufps  xmm6, xmm6,0x00
	shufps  xmm7, xmm7,0x00
	mulps  xmm6,xmm3
	mulps  xmm7,xmm3
	addps  xmm4,xmm6
	addps  xmm5,xmm7

	movups  [ecx],xmm4
	movups  [ecx+16],xmm5

	/* second half of the matrix */
	movss   xmm4,[ecx+32]
	movss   xmm6,[ecx+36]
	movss   xmm5,[ecx+48]
	movss   xmm7,[ecx+52]
	shufps  xmm4, xmm4,0x00
	shufps  xmm5, xmm5,0x00
	mulps  xmm4,xmm0
	mulps  xmm5,xmm0

	shufps xmm6, xmm6, 0x00
	shufps  xmm7, xmm7,0x00
	mulps  xmm6,xmm1
	mulps  xmm7,xmm1
	addps  xmm4,xmm6
	addps  xmm5,xmm7


	movss  xmm6,[ecx+40]
	movss  xmm7,[ecx+56]
	shufps  xmm6, xmm6,0x00
	shufps  xmm7, xmm7,0x00
	mulps  xmm6,xmm2
	mulps  xmm7,xmm2
	addps  xmm4,xmm6
	addps  xmm5,xmm7

	movss   xmm6,[ecx+44]
	movss   xmm7,[ecx+60]
	shufps  xmm6, xmm6,0x00
	shufps  xmm7, xmm7,0x00
	mulps  xmm6,xmm3
	mulps  xmm7,xmm3
	addps  xmm4,xmm6
	addps  xmm5,xmm7

	movups  [ecx+32],xmm4
	movups  [ecx+48],xmm5
	}

}

void VPCALL CSIMD_SSE::mat4D_MultiplyOf(CMat4D* pOut, const CMat4D* pLeft, const CMat4D* pRight)
{
	
	//LOADREGS2( pOut , pLeft, pRight ) 
	
	
	__asm {
	mov eax, pRight
	mov ebx, pLeft
	mov ecx, pOut
	movups  xmm0, [eax]  	/* xmm0 = pRight[0..3] */
	movups  xmm1, [eax+16]	/* xmm1 = pRight[5..7] */
	movups  xmm2, [eax+32]	/* xmm2 = pRight[8..11] */
	movups  xmm3, [eax+48]	/* xmm3 = pRight[12..15] */

	/* Processes 1/2 of the matrix at a time [2x4], unrolled loop */
	movss     xmm4, [ebx]
	movss    xmm6, [ebx+4]
	movss   xmm5, [ebx+16]
	movss   xmm7, [ebx+20]
	shufps  xmm4, xmm4,0x00
	shufps  xmm5, xmm5,0x00
	shufps  xmm6, xmm6,0x00
	shufps  xmm7, xmm7,0x00
	mulps  xmm4, xmm0
	mulps  xmm5, xmm0
	mulps  xmm6, xmm1
	mulps  xmm7, xmm1
	addps  xmm5, xmm7
	addps  xmm4, xmm6


	movss   xmm6, [ebx+8]
	movss  xmm7, [ebx+24]
	shufps  xmm6, xmm6,0x00
	shufps  xmm7, xmm7,0x00
	mulps  xmm6, xmm2
	mulps  xmm7, xmm2
	addps  xmm4, xmm6
	addps  xmm5, xmm7

	movss   xmm6, [ebx+12]
	movss   xmm7, [ebx+28]
	shufps  xmm6, xmm6,0x00
	shufps  xmm7, xmm7,0x00
	mulps  xmm6, xmm3
	mulps  xmm7, xmm3
	addps  xmm4, xmm6
	addps  xmm5, xmm7
	
	movups  [ecx], xmm4
	movups [ecx+16], xmm5

	/* second half of the matrix */
	movss   xmm4, [ebx+32]
	movss   xmm6, [ebx+36]
	movss   xmm5, [ebx+48]
	movss   xmm7, [ebx+52]
	shufps  xmm4, xmm4,0x00
	shufps  xmm5, xmm5,0x00
	mulps  xmm4, xmm0
	mulps  xmm5, xmm0

	shufps  xmm6, xmm6,0x00
	shufps  xmm7, xmm7,0x00
	mulps  xmm6, xmm1
	mulps  xmm7, xmm1
	addps  xmm4, xmm6
	addps  xmm5, xmm7


	movss  xmm6, [ebx+40]
	movss  xmm7, [ebx+56]
	shufps  xmm6, xmm6,0x00
	shufps  xmm7, xmm7,0x00
	mulps  xmm6,xmm2
	mulps  xmm7,xmm2
	addps  xmm4,xmm6
	addps  xmm5,xmm7
	movss   xmm6, [ebx+44]
	movss   xmm7, [ebx+60]
	shufps  xmm6, xmm6,0x00
	shufps  xmm7, xmm7,0x00
	mulps  xmm6,xmm3
	mulps  xmm7,xmm3
	addps  xmm4,xmm6
	addps  xmm5,xmm7

	movups  32[ecx],xmm4
	movups  48[ecx],xmm5

	}
}



void VPCALL CSIMD_SSE::mat4D_Transpose(CMat4D* pIn)
{
	LOADREGS0(pIn ) 
	
	 __asm {
		movlps	 xmm1,[ecx]
		movlps	 xmm3,[ecx+8]
		movhps	 xmm1,[ecx+16]
		movhps	 xmm3,[ecx+24]
		movlps	 xmm5,[ecx+32]
		movlps	 xmm4,[ecx+40]
		movhps	 xmm5,[ecx+48]
		movhps	 xmm4,[ecx+56]
		movaps	 xmm0,xmm1
		movaps	 xmm2,xmm3
		shufps	 xmm1,xmm5,0xDD
		shufps	 xmm3,xmm4,0xDD
		shufps	 xmm0,xmm5,0x88
		shufps	 xmm2,xmm4,0x88
		movups   [ecx],xmm0
		movups	 [ecx+16],xmm1
		movups   [ecx+32],xmm2
		movups   [ecx+48],xmm3
	}


}



void VPCALL CSIMD_SSE::mat4D_TransposeOf(CMat4D* pOut, const CMat4D* pIn)
{
	LOADREGS( pOut , pIn ) 
	
	__asm {
		movlps	 xmm1,[eax]
		movlps	 xmm3,[eax+8]
		movhps	 xmm1,[eax+16]
		movhps	 xmm3,[eax+24]
		movlps	 xmm5,[eax+32]
		movlps	 xmm4,[eax+40]
		movhps	 xmm5,[eax+48]
		movhps	 xmm4,[eax+56]
		movaps	 xmm0,xmm1
		movaps	 xmm2,xmm3
		shufps	 xmm1,xmm5,0xDD
		shufps	 xmm3,xmm4,0xDD
		shufps	 xmm0,xmm5,0x88
		shufps	 xmm2,xmm4,0x88
		movups   [ecx],xmm0
		movups   [ecx+16],xmm1
		movups   [ecx+32],xmm2
		movups   [ecx+48],xmm3
}

}

void VPCALL CSIMD_SSE::mat4D_VectorMultiply(CVec3D* pVec, const CMat4D* pMat)
{
	/*
		Does a normal vector x matrix, but since w = 1.0f, the final column of the matrix is
		merely added.
	*/

	LOADREGS( pVec,pMat ) 
	
	__asm {
	movups    xmm4,[ecx]
	movups    xmm0,[eax]
	movups  xmm1,[eax+16]
	movups  xmm2,[eax+32]
	movups  xmm3,[eax+48]
	movaps  xmm5,xmm4
	movaps  xmm6,xmm4
	shufps  xmm4, xmm4,0x00    /* xmm4 = x | x | x | x */
	shufps  xmm5, xmm5,0x55   /* xmm5 = y | y | y | y */
	shufps  xmm6, xmm6,0xAA    /* xmm6 = z | z | z | z */

	/* multiply with each row */
	mulps  xmm0,xmm4
	mulps  xmm1,xmm5
	mulps  xmm2,xmm6

	/* Sum results, xmm0-2 = transformed x,y,z colums, xmm3 = w column */
	addps  xmm1,xmm0    /* xmm1 = tx + ty */
	addps  xmm3,xmm2    /* xmm3 = tz + w */
	addps  xmm1,xmm3    /* xmm1 = tx + ty + tz + w*/

	/* Store translated vector */
	movups  [ecx],xmm1
}

}

void VPCALL CSIMD_SSE::mat4D_VectorMultiplyOf(CVec3D* pOut, const CVec3D* pIn, const CMat4D* pMat)
{
	/*
		Does a normal vector x matrix, but since w = 1.0f, the final column of the matrix is
		merely added.
	*/

	LOADREGS2( pOut , pIn, pMat ) 
	
	__asm {
	movups    xmm4,[eax]
	movups    xmm0,[ebx]
	movups  xmm1,[ebx+16]
	movups  xmm2,[ebx+32]
	movups  xmm3,[ebx+48]
	movaps  xmm5,xmm4
	movaps  xmm6,xmm4
	shufps  xmm4, xmm4,0x00    /* xmm4 = x | x | x | x */
	shufps  xmm5, xmm5,0x55    /* xmm5 = y | y | y | y */
	shufps  xmm6, xmm6,0xAA    /* xmm6 = z | z | z | z */

	/* multiply with each row */
	mulps  xmm0,xmm4
	mulps  xmm1,xmm5
	mulps  xmm2,xmm6

	/* Sum results, xmm0-2 = transformed x,y,z colums, xmm3 = w column */
	addps  xmm1,xmm0    /* xmm1 = tx + ty */
	addps  xmm3,xmm2    /* xmm3 = tz + w */
	addps  xmm1,xmm3    /* xmm1 = tx + ty + tz + w*/

	/* Store translated vector */
	movups  [ecx],xmm1
};
	
}

void VPCALL CSIMD_SSE::mat4D_VectorMultiply(CVec4D* pOut4D, const CMat4D* pMat)
{
	
	LOADREGS( pOut4D , pMat ) 
	
	__asm {

	/*
	    IMPLICIT:
		movl 4(%esp), %1
		movl 8(%esp), %0

		Access %1 first, it has finished and has no dependencies, unlike loading from %0

		Surprisingly results in 9% speedup on Celeron, 6% on Athlon64
	*/

	movups  xmm4,[ecx]
	movups  xmm0,[eax]
	movups  xmm1,[eax+16]
	movups  xmm2,[eax+32]
	movups  xmm3,[eax+48]
	movaps  xmm5,xmm4
	movaps  xmm6,xmm4
	movaps  xmm7,xmm4
	shufps  xmm4, xmm4,0x00
	shufps  xmm5, xmm5,0x55
	shufps  xmm6, xmm6,0xAA
	shufps  xmm7, xmm7,0xFF

	/* multiply with each row */
	mulps  xmm0,xmm4
	mulps  xmm1,xmm5
	mulps  xmm2,xmm6
	mulps  xmm3,xmm7

	/* Sum results */
	addps  xmm1,xmm0
	addps  xmm3,xmm2
	addps  xmm1,xmm3

	/* Store translated vector */
	movups  [ecx],xmm1
}

}


void VPCALL CSIMD_SSE::mat4D_VectorMultiplyOf(CVec4D* pOut4D, const CVec4D* pIn4D, const CMat4D* pMat)
{
	LOADREGS2(pOut4D , pIn4D, pMat ) 
	
	__asm {
	/*
	    IMPLICIT:
		movl 4(%esp), %1
		movl 8(%esp), %0

		Access %1 first, it has finished and has no dependencies, unlike loading from %0

		Surprisingly results in 9% speedup on Celeron, 6% on Athlon64
	*/
	movups  xmm4,[eax]
	movups    xmm0,[ebx]
	movups  xmm1,[ebx+16]
	movups  xmm2,[ebx+32]
	movups  xmm3,[ebx+48]
	movaps  xmm5,xmm4
	movaps  xmm6,xmm4
	movaps   xmm7,xmm4
	shufps  xmm4, xmm4,0x00
	shufps  xmm5, xmm5,0x55
	shufps  xmm6, xmm6,0xAA
	shufps  xmm7, xmm7,0xFF
	mulps  xmm0,xmm4
	mulps  xmm1,xmm5
	mulps  xmm2,xmm6
	mulps  xmm3,xmm7
	addps  xmm1,xmm0
	addps  xmm3,xmm2
	addps  xmm1,xmm3
	movups  [ecx],xmm1
}

}


void VPCALL CSIMD_SSE::mat4D_ToRotate(CMat4D* pMat, float yaw, float pitch, float roll)
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

void VPCALL CSIMD_SSE::mat4D_ToRotateOf(CMat4D* pMat, const CVec3D* pYawPitchRoll)
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


void VPCALL CSIMD_SSE::mat4D_AlignedSum(CMat4D* Out, const CMat4D* In)
{
		
			LOADREGS( Out , In )
		__asm
        {	
			movaps  xmm0,[ecx]
			movaps  xmm1,[ecx+16]
			movaps  xmm2,[ecx+32]
			movaps  xmm3,[ecx+48]
			movaps  xmm4,[eax]
			movaps  xmm5,[eax+16]	
			movaps  xmm6,[eax+32]	
			movaps  xmm7,[eax+48]
			addps	xmm4,xmm0
			addps	xmm5, xmm1
			addps	xmm6,xmm2
			addps	xmm7,xmm3
			movaps   [ecx],xmm4  
			movaps  [ecx+16],xmm5
			movaps  [ecx+32],xmm6
			movaps  [ecx+48],xmm7
		}

}


void VPCALL CSIMD_SSE::mat4D_AlignedSumOf(CMat4D* Out, const CMat4D* In1, const CMat4D* In2)
{

	LOADREGS2( Out , In1, In2 )
	 __asm {
	movaps  xmm0,[eax]
	movaps  xmm1,[eax+16]
	movaps  xmm2,[eax+32]
	movaps  xmm3,[eax+48]
	movaps  xmm4,[ebx]
	movaps  xmm5,[ebx+16]
	movaps  xmm6,[ebx+32]
	movaps  xmm7,[ebx+48]
	addps	xmm4,xmm0
	addps	xmm5,xmm1
	addps	xmm6,xmm2
	addps	xmm7,xmm3
	movaps  [ecx],xmm4
	movaps  [ecx+16],xmm5
	movaps  [ecx+32],xmm6
	movaps  [ecx+48],xmm7
	}

}

void VPCALL CSIMD_SSE::mat4D_AlignedDiff(CMat4D* Out, const CMat4D* In)
{
	
	LOADREGS( Out , In )
	 __asm {
	movaps    xmm0,[ecx]
	movaps  xmm1,[ecx+16]
	movaps  xmm2,[ecx+32]
	movaps  xmm3,[ecx+48]
	movaps    xmm4,[eax]
	movaps  xmm5,[eax+16]
	movaps  xmm6,[eax+32]
	movaps  xmm7,[eax+48]
	subps  xmm0,xmm4
	subps  xmm1,xmm5
	subps  xmm2,xmm6
	subps  xmm3,xmm7
	movaps   [ecx], xmm0
	movaps  [ecx+16],xmm1
	movaps  [ecx+32],xmm2
	movaps  [ecx+48],xmm3
	}
}

void VPCALL CSIMD_SSE::mat4D_AlignedDiffOf(CMat4D* Out, const CMat4D* In1, const CMat4D* In2)
{
	LOADREGS2( Out , In1, In2 ) 
	__asm {
	movaps  xmm0,[eax]
	movaps  xmm1,[eax+16]
	movaps  xmm2,[eax+32]
	movaps  xmm3,[eax+48]
	movaps  xmm4, [ebx]
	movaps  xmm5,[ebx+16]
	movaps  xmm6,[ebx+32]
	movaps	xmm7, [ebx+48]
	subps	xmm0,xmm4
	subps	xmm1,xmm5
	subps	xmm2,xmm6
	subps	xmm3,xmm7
	movaps  [ecx],xmm0
	movaps  [ecx+16],xmm1
	movaps  [ecx+32],xmm2
	movaps	[ecx+48],xmm3 
	}
}

void VPCALL CSIMD_SSE::mat4D_AlignedScale(CMat4D* mtx, float scalar)
{
	LOADREGS( mtx, scalar )
	 __asm {
	/* Store scalar in xmm4.x */
	movss  xmm4,eax

	/* Get the matrix into registers */
	movaps  xmm0,[ecx]
	movaps  xmm1,[ecx+16]
	movaps  xmm2,[ecx+32]
	movaps  xmm3,[ecx+48]

	/* Broadcast element x to yzw to make a duplicated scalar register */
	shufps  xmm4, xmm4, 0x00

	/* scale the matrix in parallel */
	mulps	xmm0, xmm4
	mulps	xmm1, xmm4
	mulps	xmm2, xmm4
	mulps	xmm3, xmm4

	/* Store results */
	movaps  [ecx],xmm0
	movaps  [ecx+16],xmm1
	movaps  [ecx+32],xmm2
	movaps  [ecx+48],xmm3
	}
}

void VPCALL CSIMD_SSE::mat4D_AlignedScaleOf(CMat4D* pOut, const CMat4D* pIn, float scalar)
{
	
	LOADREGS2( pOut , pIn, scalar ) 
	
		 __asm {
	mov ecx, pOut    
	mov eax, pIn		
	movss	xmm4, scalar
	movaps  xmm0, [eax]
	movaps  xmm1,[eax+16]
	movaps  xmm2,[eax+32]
	movaps  xmm3,[eax+48]
	shufps  xmm4, xmm4,0x00
	mulps	xmm0,xmm4
	mulps	xmm1,xmm4
	mulps	xmm2,xmm4
	mulps	xmm3,xmm4
	movaps  [ecx],xmm0
	movaps  [ecx+16],xmm1
	movaps  [ecx+32],xmm2
	movaps  [ecx+48],xmm3
	}

}

void VPCALL CSIMD_SSE::mat4D_AlignedMultiply(CMat4D* pLeft, const CMat4D* pRight)
{
	LOADREGS( pLeft , pRight ) 
	 
	__asm {
	movaps    xmm0,[eax]	/* xmm0 = pRight[0..3] */
	movaps  xmm1,[eax+16]	/* xmm1 = pRight[5..7] */
	movaps  xmm2,[eax+32]	/* xmm2 = pRight[8..11] */
	movaps  xmm3,[eax+48]	/* xmm3 = pRight[12..15] */

	/* Processes 1/2 of the matrix at a time (2x4), unrolled loop */
	movss    xmm4, [ecx]
	movss    xmm6,[ecx+4]
	movss  xmm5,[ecx+16]
	movss   xmm7,[ecx+20]
	shufps xmm4, xmm4, 0x00
	shufps  xmm5, xmm5,0x00
	shufps  xmm6, xmm6,0x00
	shufps  xmm7, xmm7,0x00
	mulps  xmm4,xmm0
	mulps  xmm5,xmm0
	mulps  xmm6,xmm1
	mulps  xmm7,xmm1
	addps  xmm5,xmm7
	addps  xmm4,xmm6


	movss   xmm6,[ecx+8]
	movss	xmm7,[ecx+24]
	shufps  xmm6, xmm6,0x00
	shufps xmm7, xmm7, 0x00
	mulps  xmm6,xmm2
	mulps xmm7,xmm2 
	addps xmm4,xmm6 
	addps xmm5,xmm7 

	movss   xmm6,[ecx+12]
	movss   xmm7,[ecx+28]
	shufps  xmm6, xmm6,0x00
	shufps  xmm7, xmm7,0x00
	mulps  xmm6,xmm3
	mulps  xmm7,xmm3
	addps  xmm4,xmm6
	addps  xmm5,xmm7

	movaps  [ecx],xmm4
	movaps  [ecx+16],xmm5

	/* second half of the matrix */
	movss   xmm4,[ecx+32]
	movss   xmm6,[ecx+36]
	movss   xmm5,[ecx+48]
	movss   xmm7,[ecx+52]
	shufps  xmm4, xmm4,0x00
	shufps  xmm5, xmm5,0x00
	mulps  xmm4,xmm0
	mulps  xmm5,xmm0

	shufps xmm6, xmm6, 0x00
	shufps  xmm7, xmm7,0x00
	mulps  xmm6,xmm1
	mulps  xmm7,xmm1
	addps  xmm4,xmm6
	addps  xmm5,xmm7


	movss  xmm6,[ecx+40]
	movss  xmm7,[ecx+56]
	shufps  xmm6, xmm6,0x00
	shufps  xmm7, xmm7,0x00
	mulps  xmm6,xmm2
	mulps  xmm7,xmm2
	addps  xmm4,xmm6
	addps  xmm5,xmm7

	movss   xmm6,[ecx+44]
	movss   xmm7,[ecx+60]
	shufps  xmm6, xmm6,0x00
	shufps  xmm7, xmm7,0x00
	mulps  xmm6,xmm3
	mulps  xmm7,xmm3
	addps  xmm4,xmm6
	addps  xmm5,xmm7

	movaps  [ecx+32],xmm4
	movaps  [ecx+48],xmm5
	}

}

void VPCALL CSIMD_SSE::mat4D_AlignedMultiplyOf(CMat4D* pOut, const CMat4D* pLeft, const CMat4D* pRight)
{
	
	//LOADREGS2( pOut , pLeft, pRight ) 
	
	
	__asm {
	mov eax, pRight
	mov ebx, pLeft
	mov ecx, pOut
	movaps  xmm0, [eax]  	/* xmm0 = pRight[0..3] */
	movaps  xmm1, [eax+16]	/* xmm1 = pRight[5..7] */
	movaps  xmm2, [eax+32]	/* xmm2 = pRight[8..11] */
	movaps  xmm3, [eax+48]	/* xmm3 = pRight[12..15] */

	/* Processes 1/2 of the matrix at a time [2x4], unrolled loop */
	movss     xmm4, [ebx]
	movss    xmm6, [ebx+4]
	movss   xmm5, [ebx+16]
	movss   xmm7, [ebx+20]
	shufps  xmm4, xmm4,0x00
	shufps  xmm5, xmm5,0x00
	shufps  xmm6, xmm6,0x00
	shufps  xmm7, xmm7,0x00
	mulps  xmm4, xmm0
	mulps  xmm5, xmm0
	mulps  xmm6, xmm1
	mulps  xmm7, xmm1
	addps  xmm5, xmm7
	addps  xmm4, xmm6


	movss   xmm6, [ebx+8]
	movss  xmm7, [ebx+24]
	shufps  xmm6, xmm6,0x00
	shufps  xmm7, xmm7,0x00
	mulps  xmm6, xmm2
	mulps  xmm7, xmm2
	addps  xmm4, xmm6
	addps  xmm5, xmm7

	movss   xmm6, [ebx+12]
	movss   xmm7, [ebx+28]
	shufps  xmm6, xmm6,0x00
	shufps  xmm7, xmm7,0x00
	mulps  xmm6, xmm3
	mulps  xmm7, xmm3
	addps  xmm4, xmm6
	addps  xmm5, xmm7
	
	movups  [ecx], xmm4
	movups [ecx+16], xmm5

	/* second half of the matrix */
	movss   xmm4, [ebx+32]
	movss   xmm6, [ebx+36]
	movss   xmm5, [ebx+48]
	movss   xmm7, [ebx+52]
	shufps  xmm4, xmm4,0x00
	shufps  xmm5, xmm5,0x00
	mulps  xmm4, xmm0
	mulps  xmm5, xmm0

	shufps  xmm6, xmm6,0x00
	shufps  xmm7, xmm7,0x00
	mulps  xmm6, xmm1
	mulps  xmm7, xmm1
	addps  xmm4, xmm6
	addps  xmm5, xmm7


	movss  xmm6, [ebx+40]
	movss  xmm7, [ebx+56]
	shufps  xmm6, xmm6,0x00
	shufps  xmm7, xmm7,0x00
	mulps  xmm6,xmm2
	mulps  xmm7,xmm2
	addps  xmm4,xmm6
	addps  xmm5,xmm7
	movss   xmm6, [ebx+44]
	movss   xmm7, [ebx+60]
	shufps  xmm6, xmm6,0x00
	shufps  xmm7, xmm7,0x00
	mulps  xmm6,xmm3
	mulps  xmm7,xmm3
	addps  xmm4,xmm6
	addps  xmm5,xmm7

	movaps  32[ecx],xmm4
	movaps  48[ecx],xmm5

	}
}



void VPCALL CSIMD_SSE::mat4D_AlignedTranspose(CMat4D* pIn)
{
	LOADREGS0(pIn ) 
	
	 __asm {
		movlps	 xmm1,[ecx]
		movlps	 xmm3,[ecx+8]
		movhps	 xmm1,[ecx+16]
		movhps	 xmm3,[ecx+24]
		movlps	 xmm5,[ecx+32]
		movlps	 xmm4,[ecx+40]
		movhps	 xmm5,[ecx+48]
		movhps	 xmm4,[ecx+56]
		movaps	 xmm0,xmm1
		movaps	 xmm2,xmm3
		shufps	 xmm1,xmm5,0xDD
		shufps	 xmm3,xmm4,0xDD
		shufps	 xmm0,xmm5,0x88
		shufps	 xmm2,xmm4,0x88
		movaps   [ecx],xmm0
		movaps	 [ecx+16],xmm1
		movaps   [ecx+32],xmm2
		movaps   [ecx+48],xmm3
	}


}



void VPCALL CSIMD_SSE::mat4D_AlignedTransposeOf(CMat4D* pOut, const CMat4D* pIn)
{
	LOADREGS( pOut , pIn ) 
	
	__asm {
		movlps	 xmm1,[eax]
		movlps	 xmm3,[eax+8]
		movhps	 xmm1,[eax+16]
		movhps	 xmm3,[eax+24]
		movlps	 xmm5,[eax+32]
		movlps	 xmm4,[eax+40]
		movhps	 xmm5,[eax+48]
		movhps	 xmm4,[eax+56]
		movaps	 xmm0,xmm1
		movaps	 xmm2,xmm3
		shufps	 xmm1,xmm5,0xDD
		shufps	 xmm3,xmm4,0xDD
		shufps	 xmm0,xmm5,0x88
		shufps	 xmm2,xmm4,0x88
		movaps   [ecx],xmm0
		movaps   [ecx+16],xmm1
		movaps   [ecx+32],xmm2
		movaps   [ecx+48],xmm3
}

}

void VPCALL CSIMD_SSE::mat4D_AlignedVectorMultiply(CVec3D* pVec, const CMat4D* pMat)
{
	/*
		Does a normal vector x matrix, but since w = 1.0f, the final column of the matrix is
		merely added.
	*/

	LOADREGS( pVec,pMat ) 
	
	__asm {
	movaps    xmm4,[ecx]
	movaps    xmm0,[eax]
	movaps  xmm1,[eax+16]
	movaps  xmm2,[eax+32]
	movaps  xmm3,[eax+48]
	movaps  xmm5,xmm4
	movaps  xmm6,xmm4
	shufps  xmm4, xmm4,0x00    /* xmm4 = x | x | x | x */
	shufps  xmm5, xmm5,0x55   /* xmm5 = y | y | y | y */
	shufps  xmm6, xmm6,0xAA    /* xmm6 = z | z | z | z */

	/* multiply with each row */
	mulps  xmm0,xmm4
	mulps  xmm1,xmm5
	mulps  xmm2,xmm6

	/* Sum results, xmm0-2 = transformed x,y,z colums, xmm3 = w column */
	addps  xmm1,xmm0    /* xmm1 = tx + ty */
	addps  xmm3,xmm2    /* xmm3 = tz + w */
	addps  xmm1,xmm3    /* xmm1 = tx + ty + tz + w*/

	/* Store translated vector */
	movaps  [ecx],xmm1
}

}

void VPCALL CSIMD_SSE::mat4D_AlignedVectorMultiplyOf(CVec3D* pOut, const CVec3D* pIn, const CMat4D* pMat)
{
	/*
		Does a normal vector x matrix, but since w = 1.0f, the final column of the matrix is
		merely added.
	*/

	LOADREGS2( pOut , pIn, pMat ) 
	
	__asm {
	movaps    xmm4,[eax]
	movaps    xmm0,[ebx]
	movups  xmm1,[ebx+16]
	movups  xmm2,[ebx+32]
	movups  xmm3,[ebx+48]
	movaps  xmm5,xmm4
	movaps  xmm6,xmm4
	shufps  xmm4, xmm4,0x00    /* xmm4 = x | x | x | x */
	shufps  xmm5, xmm5,0x55    /* xmm5 = y | y | y | y */
	shufps  xmm6, xmm6,0xAA    /* xmm6 = z | z | z | z */

	/* multiply with each row */
	mulps  xmm0,xmm4
	mulps  xmm1,xmm5
	mulps  xmm2,xmm6

	/* Sum results, xmm0-2 = transformed x,y,z colums, xmm3 = w column */
	addps  xmm1,xmm0    /* xmm1 = tx + ty */
	addps  xmm3,xmm2    /* xmm3 = tz + w */
	addps  xmm1,xmm3    /* xmm1 = tx + ty + tz + w*/

	/* Store translated vector */
	movups  [ecx],xmm1
};
	
}

void VPCALL CSIMD_SSE::mat4D_AlignedVectorMultiply(CVec4D* pOut4D, const CMat4D* pMat)
{
	
	LOADREGS( pOut4D , pMat ) 
	
	__asm {

	/*
	    IMPLICIT:
		movl 4(%esp), %1
		movl 8(%esp), %0

		Access %1 first, it has finished and has no dependencies, unlike loading from %0

		Surprisingly results in 9% speedup on Celeron, 6% on Athlon64
	*/

	movups  xmm4,[ecx]
	movups  xmm0,[eax]
	movups  xmm1,[eax+16]
	movups  xmm2,[eax+32]
	movups  xmm3,[eax+48]
	movaps  xmm5,xmm4
	movaps  xmm6,xmm4
	movaps  xmm7,xmm4
	shufps  xmm4, xmm4,0x00
	shufps  xmm5, xmm5,0x55
	shufps  xmm6, xmm6,0xAA
	shufps  xmm7, xmm7,0xFF

	/* multiply with each row */
	mulps  xmm0,xmm4
	mulps  xmm1,xmm5
	mulps  xmm2,xmm6
	mulps  xmm3,xmm7

	/* Sum results */
	addps  xmm1,xmm0
	addps  xmm3,xmm2
	addps  xmm1,xmm3

	/* Store translated vector */
	movaps  [ecx],xmm1
}

}


void VPCALL CSIMD_SSE::mat4D_AlignedVectorMultiplyOf(CVec4D* pOut4D, const CVec4D* pIn4D, const CMat4D* pMat)
{
	LOADREGS2(pOut4D , pIn4D, pMat ) 
	
	__asm {
	/*
	    IMPLICIT:
		movl 4(%esp), %1
		movl 8(%esp), %0

		Access %1 first, it has finished and has no dependencies, unlike loading from %0

		Surprisingly results in 9% speedup on Celeron, 6% on Athlon64
	*/
	movaps  xmm4,[eax]
	movaps  xmm0,[ebx]
	movaps  xmm1,[ebx+16]
	movaps  xmm2,[ebx+32]
	movaps  xmm3,[ebx+48]
	movaps  xmm5,xmm4
	movaps  xmm6,xmm4
	movaps   xmm7,xmm4
	shufps  xmm4, xmm4,0x00
	shufps  xmm5, xmm5,0x55
	shufps  xmm6, xmm6,0xAA
	shufps  xmm7, xmm7,0xFF
	mulps  xmm0,xmm4
	mulps  xmm1,xmm5
	mulps  xmm2,xmm6
	mulps  xmm3,xmm7
	addps  xmm1,xmm0
	addps  xmm3,xmm2
	addps  xmm1,xmm3
	movaps  [ecx],xmm1
}

}


void VPCALL CSIMD_SSE::mat4D_AlignedToRotate(CMat4D* pMat, float yaw, float pitch, float roll)
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

void VPCALL CSIMD_SSE::mat4D_AlignedToRotateOf(CMat4D* pMat, const CVec3D* pYawPitchRoll)
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

//=====================Quaternion===================================

void  VPCALL CSIMD_SSE:: quaternion_Normalize(CQuaternion* pQuat)
{

		__asm{
		mov eax , pQuat
		movups  xmm0,[eax]		/* xmm0 = w | z | y | x  */
		movaps  xmm1,xmm0	/* xmm1 = w | z | y | x */
		mulps xmm0, xmm0	/* xmm0 = w*w | z*z | y*y | x*x */
		movhlps  xmm2,xmm0	/* xmm2 = ? | ? | w*w | z*z */
		addps  xmm2,xmm0	/* xmm2 = ? | ? | w*w+y*y | z*z+x*x */
		movss  xmm0,xmm2	/* xmm0 = ? | ? | ? | z*z+x*x */
		shufps  xmm2, xmm2,0x55	/* xmm2 = w*w+y*y |  w*w+y*y |  w*w+y*y | w*w+y*y */
		addss  xmm2,xmm0
		#ifdef HIPREC
		/* Full division by magnitude */
		sqrtss  xmm2,xmm2
		shufps  xmm2, xmm2,0x00
		divps  xmm1,xmm2
		movups  [eax],xmm1
		#else
		/* multiply by reciprocal root approximation */
		rsqrtss  xmm2,xmm2
		shufps  xmm2, xmm2,0x00
		mulps  xmm1,xmm2
		movups  [eax],xmm1
		#endif
		}


}

void  VPCALL CSIMD_SSE:: quaternion_NormalizeOf(CQuaternion* pOut, const CQuaternion* pQuat)
{

		__asm{
		mov eax, pQuat
		mov ecx, pOut
		movups  xmm0,[eax]		/* xmm0 = w | z | y | x  */
		movaps  xmm1,xmm0	/* xmm1 = w | z | y | x */
		mulps  xmm0,xmm0	/* xmm0 = w*w | z*z | y*y | x*x */
		movhlps  xmm2,xmm0	/* xmm2 = ? | ? | w*w | z*z */
		addps  xmm2,xmm0	/* xmm2 = ? | ? | w*w+y*y | z*z+x*x */
		movss  xmm0,xmm2	/* xmm0 = ? | ? | ? | z*z+x*x */
		shufps  xmm2, xmm2,0x55	/* xmm2 = w*w+y*y |  w*w+y*y |  w*w+y*y | w*w+y*y */
		addss  xmm2,xmm0
		#ifdef HIPREC
		/* Full divide by square root */
		sqrtss  xmm2,xmm2			/* xmm2 = ??? | ??? | ??? | mag */
		shufps  xmm2, xmm2,0x00	/* xmm2 = mag | mag | mag | mag */
		divps  xmm1,xmm2			/* xmm1 = w/mag | z/mag | y/mag | x/mag */
		movups  [ecx],xmm1
		#else
		/* multiply by reciprocal root approximation */
		rsqrtss  xmm2,xmm2
		shufps xmm2, xmm2, 0x00
		mulps  xmm1,xmm2
		movups [ecx],xmm1 
		#endif
		}




}

void  VPCALL CSIMD_SSE:: quaternion_Multiply(CQuaternion* pLeft, const CQuaternion* pRight)
{
	__asm{
	mov eax,pRight
	mov ecx,pLeft
	movups  xmm0,[eax] 	/* xmm0 = Rw | Rz | Ry | Rx */
	movups xmm4, [ecx] 	/* xmm1 = Lw | Lz | Ly | Lx */


	/* Duplicate right throughout xmm0-xmm3 */
	movaps  xmm1,xmm0	/* xmm1 = xmm0 */
	movaps  xmm2,xmm0	/* xmm2 = xmm0 */
	movaps  xmm3,xmm0	/* xmm3 = xmm0 */

	/* Duplicate left throughout xmm4-xmm7 */
	movaps  xmm5,xmm4	/* xmm5 = xmm4 */
	movaps  xmm6,xmm4	/* xmm6 = xmm4 */
	movaps  xmm7,xmm4	/* xmm7 = xmm4 */

	/*
		Broadcast elements in xmm4-xmm7 to create scalar registers
		==================
		0000 0000 = xxxx = 0x00
	    0101 0101 = yyyy = 0x55
		1010 1010 = zzzz = 0xAA
		1111 1111 = wwww = 0xFF
	*/
	shufps  xmm4, xmm4,0x00		/* xmm4 = Rx | Rx | Rx | Rx */
	shufps  xmm5, xmm5,0x55		/* xmm5 = Ry | Ry | Ry | Ry */
	shufps  xmm6, xmm6,0xAA		/* xmm6 = Rz | Rz | Rz | Rz */
	shufps  xmm7, xmm7,0xFF		/* xmm7 = Rw | Rw | Rw | Rw */

	/*
		set up columns
		==============
		C1 = w | z | y | x = 1110 0100 =
		C2 = x | y | z | w = 0001 1011 = 0x1B
		C3 = y | x | w | z = 0100 1110 = 0x4E
		C4 = z | w | x | y = 1011 0001 = 0xB1
	*/

	/* C1 is already w | z | y | x  format, no shufps needed */
	shufps  xmm1, xmm1,0x1B
	shufps  xmm2, xmm2,0x4E
	shufps  xmm3, xmm3,0xB1

	/* multiply columns */
	mulps  xmm7,xmm0		/* C1 *= Lw */
	mulps  xmm4,xmm1		/* C2 *= Lx */
	mulps  xmm5,xmm2		/* C3 *= Ly */
	mulps  xmm6,xmm3		/* C4 *= Lz */

	/* Change the signs of the columns (C1, aka: xmm4, doesnt need it)*/
	xorps  xmm4,POSNEGPOSNEG		/* C2 = { + - + - } */
	xorps  xmm5,POSPOSNEGNEG		/* C3 = { + + - - } */
	xorps  xmm6,NEGPOSPOSNEG		/* C4 = { - + + - } */

	addps  xmm5,xmm4		/* C2 += C1 */
	addps  xmm7,xmm6		/* C4 += C3 */
	addps  xmm7,xmm5		/* C4 += C2 */

	movups  [ecx],xmm7			/* xmm7 = new quaternion, write it out */
	}


}

void  VPCALL CSIMD_SSE:: quaternion_MultiplyOf(CQuaternion* pOut, const CQuaternion* pLeft, const CQuaternion* pRight)
{

	__asm{
	mov eax,pRight
	mov ebx,pLeft
	mov ecx,pOut
	movups  xmm0,[eax] 	/* xmm0 = Rw | Rz | Ry | Rx */
	movups  xmm4,[ebx] 	/* xmm1 = Lw | Lz | Ly | Lx */


	/* Duplicate right throughout xmm0-xmm3 */
	movaps  xmm1,xmm0	/* xmm1 = xmm0 */
	movaps  xmm2,xmm0	/* xmm2 = xmm0 */
	movaps  xmm3,xmm0	/* xmm3 = xmm0 */

	/* Duplicate left throughout xmm4-xmm7 */
	movaps  xmm5,xmm4	/* xmm5 = xmm4 */
	movaps  xmm6,xmm4	/* xmm6 = xmm4 */
	movaps  xmm7,xmm4	/* xmm7 = xmm4 */

	/*
		Broadcast elements in xmm4-xmm7 to create scalar registers
		==================
		0000 0000 = xxxx = 0x00
	    0101 0101 = yyyy = 0x55
		1010 1010 = zzzz = 0xAA
		1111 1111 = wwww = 0xFF
	*/
	shufps  xmm4, xmm4,0x00		/* xmm4 = Rx | Rx | Rx | Rx */
	shufps  xmm5, xmm5,0x55		/* xmm5 = Ry | Ry | Ry | Ry */
	shufps  xmm6, xmm6,0xAA		/* xmm6 = Rz | Rz | Rz | Rz */
	shufps  xmm7, xmm7,0xFF	/* xmm7 = Rw | Rw | Rw | Rw */

	/*
		set up columns
		==============
		C1 = w | z | y | x = 1110 0100 =
		C2 = x | y | z | w = 0001 1011 = 0x1B
		C3 = y | x | w | z = 0100 1110 = 0x4E
		C4 = z | w | x | y = 1011 0001 = 0xB1
	*/

	/* C1 is already w | z | y | x  format, no shufps needed */
	shufps  xmm1, xmm1,0x1B
	shufps  xmm2, xmm2,0x4E
	shufps  xmm3, xmm3,0xB1

	/* multiply columns */
	mulps  xmm7,xmm0		/* C1 *= Lw */
	mulps  xmm4,xmm1		/* C2 *= Lx */
	mulps  xmm5,xmm2		/* C3 *= Ly */
	mulps  xmm6,xmm3		/* C4 *= Lz */

	/* Change the signs of the columns (C1, aka: xmm4, doesnt need it)*/
	xorps  xmm4, POSNEGPOSNEG		/* C2 = { + - + - } */
	xorps  xmm5, POSPOSNEGNEG		/* C3 = { + + - - } */
	xorps  xmm6, NEGPOSPOSNEG		/* C4 = { - + + - } */

	addps  xmm5,xmm4		/* C2 += C1 */
	addps  xmm7,xmm6		/* C4 += C3 */
	addps  xmm7,xmm5		/* C4 += C2 */

	movups  [ecx],xmm7			/* xmm7 = new quaternion, write it out */
	}



}

//========CPlane==========================================



void CSIMD_SSE:: plane_FromPoints(CPlane* pOut, const CVec3D* pA, const CVec3D* pB,const CVec3D* pC)
{
	__asm {
	mov eax, pA
	mov ebx, pB
	mov edx, pC
	mov ecx, pOut
	movups  xmm0,[eax]	/* pA */
	movups  xmm1,[ebx]	/* pB */
	movups  xmm2,[edx]	/* pC */
	movaps  xmm7,xmm0   /* Save 'a' into xmm7 */
	subps  xmm0,xmm1
	subps  xmm2,xmm1

	/* Now, just need the cross product of xmm0 and xmm2 */
	movaps  xmm1,xmm0	/* xmm0 = xmm1 = left */
	movaps  xmm3,xmm2	/* xmm2 = xmm3 = right */
	shufps  xmm0, xmm0,0xC9	/* left.yxz */
	shufps  xmm1, xmm1,0xD2	/* left.xzy */
	shufps  xmm2, xmm2,0xD2	/* right.xzy */
	shufps  xmm3, xmm3,0xC9	/* right.yxz */

	/* multiply columns 1&2 and 3&4 */
	mulps  xmm2,xmm0
	mulps  xmm3,xmm1

	/* Got the cross product, OK */
	subps  xmm2,xmm3   /* xmm2 = cross (left x right)

	/* Begin calculation of 'd' component */


	/* AND off bits 96-127 (w component) */
	andps  xmm2,_SIMDx86_float_SSE_NO_W_MASK

	/* save xmm4 = 0 | z | y | x */
	movaps  xmm4,xmm2
	/* d = -(pa x cross)
	/* multiply with point 'a' on the polygon (saved in xmm7): xmm2 = 0 | a.z*z | a.y*y | a.x*x */
	mulps  xmm2,xmm7

	
	movhlps  xmm3,xmm2		/* xmm3 = ?   | ?   | 0   |  z^2 */
	addss  xmm3,xmm2		/* xmm3 = ?   | ?   | 0   | z^2 + x^2 */
	shufps  xmm2, xmm2,0x55	/* xmm2 = y^2 | y^2 | y^2 | y^2 */
	addss  xmm3,xmm2		/* xmm3 = ?   | ?   | ?   | x^2+y^2+z^2 */
	/* bytes 0-31 de xmm3 has d*/
	/* Change sign */
	xorps  xmm3, _SIMDx86_float_NEGPOSPOSPOS /* xmm3 = ? | ? | ? | -(x^2+y^2+z^2)*/

	/* Move to w component location, mask off xyz, and OR with saved portions */
	shufps  xmm3, xmm3,0x00				/* xmm3 = -(x^2+y^2+z^2) | -(x^2+y^2+z^2) | -(x^2+y^2+z^2) | -(x^2+y^2+z^2)*/
	andps  xmm3,_SIMDx86_float_SSE_NO_XYZ_MASK	/* xmm3 = -(x^2+y^2+z^2) | 0 | 0 | 0 */
	orps  xmm4,xmm3							/* xmm4 = -(x^2+y^2+z^2) | z | y | x */

	

	/* Save plane coefficients */
	movups  [ecx],xmm4
}



}

float CSIMD_SSE:: plane_DistToPoint(const CPlane* pPlane, const CVec3D* pPoint)
{
	float dummy;
	#if 0
	__asm{
mov     eax, pPlane
mov     ecx, pPoint
movss   xmm0, dword ptr [eax]
mulss   xmm0, dword ptr [ecx]
movss   xmm1, dword ptr [eax+4]
mulss   xmm1, dword ptr [ecx+4]
addss   xmm0, xmm1
movss   xmm1, dword ptr [eax+8]
mulss   xmm1, dword ptr [ecx+8]
addss   xmm0, xmm1
addss   xmm0, dword ptr [eax+0Ch]
andps  xmm0,_SIMDx86_float_ABS	/* xmm1 = ??? | ??? | ??? | fabsf(dot(pPlane, pPoint)) */
	
	movss  [esp-4],xmm0
	
	}
#endif	
	__asm {
	mov eax, pPlane
	mov ebx, pPoint
	movups  xmm0,[eax] /* xmm0 = pPlane->d | pPlane->c | pPlane->b | pPlane->a */
	movups  xmm1,[ebx] /* xmm1 = ????????? | pPoint->z | pPoint->y | pPoint->x */
	
	andps  xmm1,_SIMDx86_float_SSE_NO_W_MASK   /* xmm1 = 0 | ... */
	movaps  xmm7,xmm0						/* xmm7 = pPlane... */

	mulps  xmm1,xmm0                        /* xmm1 = d*0.0 | c*z | b*y | a*x */
	shufps  xmm7, xmm7,0xFF				/* xmm7 = d | d | d | d */
	
	
    movhlps  xmm2,xmm1      /* xmm2 = ???? | ???? | d*0.0 | z*c */
	addss  xmm2,xmm1        /* xmm2 = ???? | ???? | ????  | x*a + z*c*/
	shufps  xmm1, xmm1,0x55 /* xmm1 = ???? | ???? | ????  | y*b */
	andps  xmm1,_SIMDx86_float_ABS	/* xmm1 = ??? | ??? | ??? | fabsf(dot(pPlane, pPoint)) */
	addss  xmm1,xmm2        /* xmm1 = ???? | ???? | ????? |  fabsf(dot(pPlane, pPoint)) + pPlane->d */
	
	movss  [esp-4],xmm1
	//flds   [esp-4]
	}

}

float CSIMD_SSE:: plane_Dot(const CPlane* pPlane, const CVec3D* pVec)
{
    float dummy;
	__asm {
	mov eax, pPlane
	mov ebx, pVec
	movups  xmm1,[ebx] /* xmm1 = ?    | V->z | V->x | V->y */
	movups  xmm0,[eax] /* xmm0 = P->d | P->c | P->b | P->a */
	andps  xmm1,_SIMDx86_float_SSE_NO_W_MASK /* 0 | z | y | x */
	mulps  xmm1,xmm0    /* xmm0 = 0 | P->c*V->z | P->b*V->y | P->a*V->x */

	
		shufps  xmm0, xmm0,0xFF
		movhlps  xmm7,xmm1		/* xmm7 = ?   | ?   | 0   |  z's */
		addss  xmm7,xmm1		/* xmm7 = ?   | ?   | 0   | z's + x's */
		shufps  xmm1, xmm1,0x55	/* xmm1 = y's | y's | y's | y's */
		addss  xmm7,xmm0        /* xmm7  = ?  | ?   | ?   | z's + x's + d*/
		addss  xmm7,xmm1		/* xmm7 = ?   | ?   | ?   | x's + y's + z's + d */
		movss  [esp-4],xmm7
		//flds -4(esp)
		}

	
}

float CSIMD_SSE:: plane_Dot4(const CPlane* pPlane, const CVec4D* pVec4)
{
    float dummy;
	__asm {
	mov eax, pPlane
	mov ebx, pVec4
	movups  xmm0,[eax] /* xmm0 = d | c | b | a */
	movups  xmm1,[ebx] /* xmm1 = w | z | y | z */
	mulps  xmm1,xmm0	/* xmm1 = dw | cz | by | az */

	
		movhlps  xmm7,xmm1		/* xmm7 = ? | ? | dw | cz */
		addps  xmm7,xmm1        /* xmm7 = ? | ? | dw+by | cz+ax */
		movaps  xmm1,xmm7       /* xmm1 = ? | ? | dw+by | cz+ax */
		shufps  xmm1, xmm1,0x55	/* xmm1 = ? | ? | ? | dw+by */
		addss  xmm7,xmm1        /* xmm7  = ?  | ?   | ?   | dw+cz+by+ax */
		movss  [esp-4],xmm7
		//flds -4(esp)
		}

	
}

float CSIMD_SSE:: plane_DotNormal(const CPlane* pPlane, const CVec3D* pVec)
{
	float dummy;
	__asm {
	mov eax, pPlane
	mov ebx, pVec
	movups  xmm0,[eax] /* xmm0 = P->d | P->c | P->b | P->a */
	movups  xmm1,[ebx] /* xmm1 = ?    | V->z | V->b | V->a */
	andps  xmm0,_SIMDx86_float_SSE_NO_W_MASK	 /* 0 | P1.c | P1.b | P1.a */

	mulps  xmm0,xmm1    /* xmm0 = 0 | P->c*V->z | P->b*V->y | P->a*V->x */

		movhlps  xmm1,xmm0		/* xmm1 = ?   | ?   | 0   |  z's */
		addss  xmm1,xmm0		/* xmm1 = ?   | ?   | 0   | z's + x's */
		shufps  xmm0, xmm0,0x55	/* xmm0 = y's | y's | y's | y's */
		addss  xmm1,xmm0		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */

		movss  [esp-4],xmm1
		//flds	[esp-4]
		}

	
}

float CSIMD_SSE:: plane_DotPlane(const CPlane* pPlane1, const CPlane* pPlane2)
{
	float dummy;
	__asm {
	mov eax, pPlane1
	mov ebx, pPlane2
	movups  xmm0,[eax] /* xmm0 = P1.d | P1.c | P1.b | P1.a */
	movups  xmm1,[ebx] /* xmm1 = P2.d | P2.c | P2.b | P2.a */
	andps  xmm0,_SIMDx86_float_SSE_NO_W_MASK	 /* 0 | P1.c | P1.b | P1.a */
	
	mulps  xmm0,xmm1    /* xmm0 = 0 | P1.c*P2.c | P1.b*P2.b | P1.a*P2.a */
	
	movhlps  xmm1,xmm0		/* xmm1 = ?   | ?   | 0   |  z's */
	addss  xmm1,xmm0		/* xmm1 = ?   | ?   | 0   | z's + x's */
	shufps  xmm0, xmm0,0x55	/* xmm0 = y's | y's | y's | y's */
	addss  xmm1,xmm0		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */

	movss	[esp-4],xmm1
	//flds	[esp-4]
	}
	

}

void CSIMD_SSE:: plane_Normalize(CPlane* pOut)
{
	__asm {
	mov eax, pOut
	movups  xmm0,[eax]								/* xmm0 = d | c | b | a */
	movaps  xmm1,xmm0							/* xmm1 = d | c | b | a */
	andps  xmm0,_SIMDx86_float_SSE_NO_W_MASK		/* xmm0 = 0 | c | b | a */

	/* Perform dot product */
	mulps  xmm0,xmm0							/* xmm0 = 0 | c*c | b*b | a*a */
	
	movhlps  xmm3,xmm0							/* xmm3 = ? | ? | 0 | c*c */
	addss  xmm3,xmm0							/* xmm3 = ? | ? | ? | a*a + c*c */
	shufps  xmm0, xmm0,0x55					/* xmm0 = b*b | b*b | b*b | b*b */
	addss  xmm0,xmm3							/* xmm0 = ? | ? | ? | a*b + b*b + c*c */
	

	/* Divide it all by the square root */
	#if defined(HIPREC)
	sqrtss  xmm0,xmm0	/* xmm0 = sqrtf(a*a+b*b+c*c) */
	shufps  xmm0,xmm0,0x00	/* xmm0 = mag | mag | mag | mag */
	divps  xmm1,xmm0	/* xmm1 = d/mag | c/mag | b/mag | a/mag */
	#else
	rsqrtss  xmm0,xmm0	/* xmm0 = ? | ? | ? | 1.0f / sqrtf(a*a+b*b+c*c) */
	mulps  xmm1,xmm0	/* xmm1 = d*invmag | c*invmag | b*invmag | a*invmag */	
	
	#endif
	movups  [eax],xmm1
	}

	return;

	

}

void CSIMD_SSE:: plane_NormalizeOf(CPlane* pOut, CPlane* pIn)
{

	__asm {
	mov eax, pIn
	mov ecx, pOut
	movups  xmm0,[eax]								/* xmm0 = d | c | b | a */
	movaps  xmm1,xmm0							/* xmm1 = d | c | b | a */
	andps  xmm0,_SIMDx86_float_SSE_NO_W_MASK		/* xmm0 = 0 | c | b | a */

	/* Perform dot product */
	mulps xmm0, xmm0							/* xmm0 = 0 | c*c | b*b | a*a */
	movhlps  xmm3,xmm0							/* xmm3 = ? | ? | 0 | c*c */
	addss  xmm3,xmm0							/* xmm3 = ? | ? | ? | a*a + c*c */
	shufps  xmm0, xmm0,0x55					/* xmm0 = b*b | b*b | b*b | b*b */
	addss  xmm0,xmm3							/* xmm0 = ? | ? | ? | a*b + b*b + c*c */
	
	/* Divide it all by the square root */
	#if defined(HIPREC)
	sqrtss xmm0, xmm0	/* xmm0 = sqrtf(a*a+b*b+c*c) */
	shufps  xmm0,xmm0,0x00	/* xmm0 = mag | mag | mag | mag */
	divps  xmm1, xmm0	/* xmm1 = d/mag | c/mag | b/mag | a/mag */
	#else
	rsqrtss  xmm0,xmm0	/* xmm0 = ? | ? | ? | 1.0f / sqrtf(a*a+b*b+c*c) */
	mulps  xmm1,xmm0	/* xmm1 = d*invmag | c*invmag | b*invmag | a*invmag */
	#endif

	movups  [ecx],xmm1
	}
	return;

}


#endif /* WIN32 */

} //end MATH

} //end SMF
