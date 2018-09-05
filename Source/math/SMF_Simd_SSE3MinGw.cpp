/*
  SMF - Salvathor Math Fabric  (https://sourceforge.net/projects/smfabric/)
  Copyright (C) 2010-2011 Salvatore Giannotta Filho <a_materasu@hotmail.com>

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
#include "math/SMF_SimdMMX.h"
#include "math/SMF_SimdSSE.h"
#include "math/SMF_SimdSSE2.h"
#include "math/SMF_SimdSSE3.h"
#include "math/SMF_Matriz.h"
#include "math/SMF_Vector.h"
#include "geometry/SMF_DrawVert.h"
#include "math/SMF_JointTransform.h"

namespace SMF{
namespace MATH{

ALIGN4_INIT1( unsigned long SIMD_DW_mat2quatShuffle0, (3<<0)|(2<<8)|(1<<16)|(0<<24) );
ALIGN4_INIT1( unsigned long SIMD_DW_mat2quatShuffle1, (0<<0)|(1<<8)|(2<<16)|(3<<24) );
ALIGN4_INIT1( unsigned long SIMD_DW_mat2quatShuffle2, (1<<0)|(0<<8)|(3<<16)|(2<<24) );
ALIGN4_INIT1( unsigned long SIMD_DW_mat2quatShuffle3, (2<<0)|(3<<8)|(0<<16)|(1<<24) );

ALIGN4_INIT4( unsigned int  _SIMDx86_float_POSPOSPOSNEG, 0x00000000, 0x00000000, 0x00000000, 0x80000000 );
ALIGN4_INIT4( unsigned int  _SIMDx86_float_SSE_NO_W_MASK, 0xFFFFFFFF,  0xFFFFFFFF, 0xFFFFFFFF, 0x00000000 );
ALIGN4_INIT4( unsigned int  _SIMDx86_float_SSE_NO_XYZ_MASK, 0x00000000, 0x00000000, 0x00000000, 0xFFFFFFFF );
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
//
//	SSE3 implementation of CSIMDProcessor
//
//===============================================================

#if defined(MACOS_X) && defined(__i386__)

/*
============
CSIMD_SSE3::getName
============
*/
const char * CSIMD_SSE3::getName() const {
	return "MMX & SSE & SSE2 & SSE3";
}

#elif defined(WIN32)

#include <xmmintrin.h>


/*

	The first argument of an instruction macro is the destination
	and the second argument is the source operand. The destination
	operand can be _xmm0 to _xmm7 only. The source operand can be
	any one of the registers _xmm0 to _xmm7 or _eax, _ecx, _edx, _esp,
	_ebp, _ebx, _esi, or _edi that contains the effective address.

	For instance:  haddps   xmm0, xmm1
	becomes:       haddps( _xmm0, _xmm1 )
	and:           haddps   xmm0, [esi]
	becomes:       haddps( _xmm0, _esi )

	The ADDRESS_ADDC macro can be used when the effective source address
	is formed by adding a constant to a general purpose register.
	For instance:  haddps   xmm0, [esi+48]
	becomes:       haddps( _xmm0, ADDRESS_ADDC( _esi, 48 ) )

	The ADDRESS_ADDR macro can be used when the effective source address
	is formed by adding two general purpose registers.
	For instance:  haddps   xmm0, [esi+eax]
	becomes:       haddps( _xmm0, ADDRESS_ADDR( _esi, _eax ) )

	The ADDRESS_ADDRC macro can be used when the effective source address
	is formed by adding two general purpose registers and a constant.
	The constant must be in the range [-128, 127].
	For instance:  haddps   xmm0, [esi+eax+48]
	becomes:       haddps( _xmm0, ADDRESS_ADDRC( _esi, _eax, 48 ) )

	The ADDRESS_SCALEADDR macro can be used when the effective source address is formed
	by adding a scaled general purpose register to another general purpose register.
	The scale must be either 1, 2, 4 or 8.
	For instance:  haddps   xmm0, [esi+eax*4]
	becomes:       haddps( _xmm0, ADDRESS_SCALEADDR( _esi, _eax, 4 ) )

	The ADDRESS_SCALEADDRC macro can be used when the effective source address is formed
	by adding a scaled general purpose register to another general purpose register and
	also adding a constant. The scale must be either 1, 2, 4 or 8. The constant must
	be in the range [-128, 127].
	For instance:  haddps   xmm0, [esi+eax*4+64]
	becomes:       haddps( _xmm0, ADDRESS_SCALEADDRC( _esi, _eax, 4, 64 ) )

*/

#define _eax	0x00
#define _ecx	0x01
#define _edx	0x02
#define _ebx	0x03
#define _esp	0x04
#define _ebp	0x05
#define _esi	0x06
#define _edi	0x07

#define _xmm0	0xC0
#define _xmm1	0xC1
#define _xmm2	0xC2
#define _xmm3	0xC3
#define _xmm4	0xC4
#define _xmm5	0xC5
#define _xmm6	0xC6
#define _xmm7	0xC7

#define RSCALE( s )		( (s&2)<<5 ) | ( (s&4)<<5 ) | ( (s&8)<<3 ) | ( (s&8)<<4 )

#if defined(_MSC_VER)
#define ADDRESS_ADDC( reg0, constant )						0x40 | ( reg0 & 7 )	\
	_asm _emit constant

#define ADDRESS_ADDR( reg0, reg1 )							0x04				\
	_asm _emit ( ( reg1 & 7 ) << 3 ) | ( reg0 & 7 )

#define ADDRESS_ADDRC( reg0, reg1, constant )				0x44				\
	_asm _emit ( ( reg1 & 7 ) << 3 ) | ( reg0 & 7 )								\
	_asm _emit constant

#define ADDRESS_SCALEADDR( reg0, reg1, scale )				0x04				\
	_asm _emit ( ( reg1 & 7 ) << 3 ) | ( reg0 & 7 ) | RSCALE( scale )

#define ADDRESS_SCALEADDRC( reg0, reg1, scale, constant )	0x44				\
	_asm _emit ( ( reg1 & 7 ) << 3 ) | ( reg0 & 7 ) | RSCALE( scale )			\
	_asm _emit constant


// Packed Single-FP add/Subtract ( dst[0]=dst[0]+src[0], dst[1]=dst[1]-src[1], dst[2]=dst[2]+src[2], dst[3]=dst[3]-src[3] )
#define addsubps( dst, src )						\
	_asm _emit 0xF2									\
	_asm _emit 0x0F									\
	_asm _emit 0xD0									\
	_asm _emit ( ( dst & 7 ) << 3 ) | src

// Packed Double-FP add/Subtract ( dst[0]=dst[0]+src[0], dst[1]=dst[1]-src[1] )
#define addsubpd( dst, src )						\
	_asm _emit 0x66									\
	_asm _emit 0x0F									\
	_asm _emit 0xD0									\
	_asm _emit ( ( dst & 7 ) << 3 ) | src

// Packed Single-FP Horizontal add ( dst[0]=dst[0]+dst[1], dst[1]=dst[2]+dst[3], dst[2]=src[0]+src[1], dst[3]=src[2]+src[3] )
#define haddps( dst, src )							\
	_asm _emit 0xF2									\
	_asm _emit 0x0F									\
	_asm _emit 0x7C									\
	_asm _emit ( ( dst & 7 ) << 3 ) | src

// Packed Double-FP Horizontal add ( dst[0]=dst[0]+dst[1], dst[1]=src[0]+src[1] )
#define haddpd( dst, src )							\
	_asm _emit 0x66									\
	_asm _emit 0x0F									\
	_asm _emit 0x7C									\
	_asm _emit ( ( dst & 7 ) << 3 ) | src

// Packed Single-FP Horizontal Subtract ( dst[0]=dst[0]-dst[1], dst[1]=dst[2]-dst[3], dst[2]=src[0]-src[1], dst[3]=src[2]-src[3] )
#define hsubps( dst, src )							\
	_asm _emit 0xF2									\
	_asm _emit 0x0F									\
	_asm _emit 0x7D									\
	_asm _emit ( ( dst & 7 ) << 3 ) | src

// Packed Double-FP Horizontal Subtract ( dst[0]=dst[0]-dst[1], dst[1]=src[0]-src[1] )
#define hsubpd( dst, src )							\
	_asm _emit 0x66									\
	_asm _emit 0x0F									\
	_asm _emit 0x7D									\
	_asm _emit ( ( dst & 7 ) << 3 ) | src

// Move Packed Single-FP Low and Duplicate ( dst[0]=src[0], dst[1]=src[0], dst[2]=src[2], dst[3]=src[2] )
#define movsldup( dst, src )						\
	_asm _emit 0xF3									\
	_asm _emit 0x0F									\
	_asm _emit 0x12									\
	_asm _emit ( ( dst & 7 ) << 3 ) | src

// Move One Double-FP Low and Duplicate ( dst[0]=src[0], dst[1]=src[0] )
#define movdldup( dst, src )						\
	_asm _emit 0xF2									\
	_asm _emit 0x0F									\
	_asm _emit 0x12									\
	_asm _emit ( ( dst & 7 ) << 3 ) | src

// Move Packed Single-FP High and Duplicate ( dst[0]=src[1], dst[1]=src[1], dst[2]=src[3], dst[3]=src[3] )
#define movshdup( dst, src )						\
	_asm _emit 0xF3									\
	_asm _emit 0x0F									\
	_asm _emit 0x16									\
	_asm _emit ( ( dst & 7 ) << 3 ) | src

// Move One Double-FP High and Duplicate ( dst[0]=src[1], dst[1]=src[1] )
#define movdhdup( dst, src )						\
	_asm _emit 0xF2									\
	_asm _emit 0x0F									\
	_asm _emit 0x16									\
	_asm _emit ( ( dst & 7 ) << 3 ) | src

// Load Unaligned Integer 128 bits
#define lddqu( dst, src )							\
	_asm _emit 0xF2									\
	_asm _emit 0x0F									\
	_asm _emit 0xF0									\
	_asm _emit ( ( dst & 7 ) << 3 ) | src
#else   //mingw
#define STRINGIZE(X) #X
#define ADDRESS_ADDC( reg0, constant )						0x40 | ( reg0 & 7 )	\
asm volatile (".byte constant")

#define ADDRESS_ADDR( reg0, reg1 )							0x04				\
asm volatile (".byte ( ( reg1 & 7 ) << 3 ) | ( reg0 & 7 )")

#define ADDRESS_ADDRC( reg0, reg1, constant )				0x44				\
asm volatile (".byte ( ( reg1 & 7 ) << 3 ) | ( reg0 & 7 )")							\
asm volatile (".byte  constant")

#define ADDRESS_SCALEADDR( reg0, reg1, scale )				0x04				\
asm volatile (".byte ( ( reg1 & 7 ) << 3 ) | ( reg0 & 7 ) | RSCALE( scale )")

#define ADDRESS_SCALEADDRC( reg0, reg1, scale, constant )	0x44				\
asm volatile (".byte ( ( reg1 & 7 ) << 3 ) | ( reg0 & 7 ) | RSCALE( scale )	")		\
asm volatile (".byte constant")


// Packed Single-FP add/Subtract ( dst[0]=dst[0]+src[0], dst[1]=dst[1]-src[1], dst[2]=dst[2]+src[2], dst[3]=dst[3]-src[3] )
#define addsubps( dst, src )						\
".byte 0xF2\n"                    \
".byte 0x0F\n"                    \
".byte 0xD0\n"                    \
".byte  ( ( dst & 7 ) << 3 ) | src\n"

// Packed Double-FP add/Subtract ( dst[0]=dst[0]+src[0], dst[1]=dst[1]-src[1] )
#define addsubpd( dst, src )					\
".byte 0x66\n"									\
".byte 0x0F\n"									\
".byte 0xD0\n"									\
".byte ( ( dst & 7 ) << 3 ) | src\n"

// Packed Single-FP Horizontal add ( dst[0]=dst[0]+dst[1], dst[1]=dst[2]+dst[3], dst[2]=src[0]+src[1], dst[3]=src[2]+src[3] )
#define haddps( dst, src )						\
".byte 0xF2\n"									\
".byte 0x0F\n"									\
".byte 0x7C\n"									\
".byte  ( ( "STRINGIZE(dst)" & 7 ) << 3 ) | "STRINGIZE(src)"\n"

// Packed Double-FP Horizontal add ( dst[0]=dst[0]+dst[1], dst[1]=src[0]+src[1] )
#define haddpd( dst, src )						\
".byte 0x66\n"									\
".byte 0x0F\n"									\
".byte 0x7C\n"									\
".byte  ( ( "STRINGIZE(dst)" & 7 ) << 3 ) | "STRINGIZE(src)"\n"

// Packed Single-FP Horizontal Subtract ( dst[0]=dst[0]-dst[1], dst[1]=dst[2]-dst[3], dst[2]=src[0]-src[1], dst[3]=src[2]-src[3] )
#define hsubps( dst, src )							\
".byte 0xF2\n"									\
".byte 0x0F\n"									\
".byte 0x7D\n"									\
".byte   ( ( "STRINGIZE(dst)" & 7 ) << 3 ) | "STRINGIZE(src)"\n"

// Packed Double-FP Horizontal Subtract ( dst[0]=dst[0]-dst[1], dst[1]=src[0]-src[1] )
#define hsubpd( dst, src )							\
".byte 0x66\n"									\
".byte 0x0F\n"									\
".byte 0x7D\n"									\
".byte  ( ( "STRINGIZE(dst)" & 7 ) << 3 ) | "STRINGIZE(src)"\n"

// Move Packed Single-FP Low and Duplicate ( dst[0]=src[0], dst[1]=src[0], dst[2]=src[2], dst[3]=src[2] )
#define movsldup( dst, src )						\
".byte 0xF3\n"									\
".byte 0x0F\n"									\
".byte 0x12\n"								\
".byte  ( ( "STRINGIZE(dst)" & 7 ) << 3 ) | "STRINGIZE(src)"\n"

// Move One Double-FP Low and Duplicate ( dst[0]=src[0], dst[1]=src[0] )
#define movdldup( dst, src )						\
".byte 0xF2\n"									\
".byte 0x0F\n"									\
".byte 0x12\n"									\
".byte   ( ( "STRINGIZE(dst)" & 7 ) << 3 ) | "STRINGIZE(src)"\n"

// Move Packed Single-FP High and Duplicate ( dst[0]=src[1], dst[1]=src[1], dst[2]=src[3], dst[3]=src[3] )
#define movshdup( dst, src )						\
".byte 0xF3\n"									\
".byte 0x0F\n"									\
".byte 0x16\n"									\
".byte  ( ( "STRINGIZE(dst)" & 7 ) << 3 ) | "STRINGIZE(src)"\n"

// Move One Double-FP High and Duplicate ( dst[0]=src[1], dst[1]=src[1] )
#define movdhdup( dst, src )						\
".byte 0xF2\n"									\
".byte 0x0F\n"									\
".byte 0x16\n"									\
".byte  ( ( "STRINGIZE(dst)" & 7 ) << 3 ) | "STRINGIZE(src)"\n"

// Load Unaligned Integer 128 bits
#define lddqu( dst, src )							\
".byte 0xF2\n"									\
".byte 0x0F\n"									\
".byte 0xF0\n"									\
".byte   ( ( "STRINGIZE(dst)" & 7 ) << 3 ) | "STRINGIZE(src)"\n"


#endif


/*
============
SSE3_Dot
============
*/
float SSE3_Dot( const CVec4D &v1, const CVec4D &v2 ) {
	float d;
	asm volatile (
		"mov		esi, %0\n"
		"mov		edi, %1\n"
		"movups	xmm0, [esi]\n"
		"mulps	xmm0, [edi]\n"

		haddps(	_xmm0, _xmm0 )
		haddps(	_xmm0, _xmm0 )

		"movss	%2, xmm0\n"
	:
    :"r"(&v1),"r"(&v2),"m"(d)
    :"esi","edi");

	return d;
}

/*
============
CSIMD_SSE3::getName
============
*/
const char * CSIMD_SSE3::getName() const {
	return "MMX & SSE & SSE2 & SSE3";
}

/*
============
CSIMD_SSE3::transformVerts
============
*/
void VPCALL CSIMD_SSE3::transformVerts( CVertex *verts, const int numVerts, const CMatJoint3x4 *joints, const CVec4D *weights, const int *index, const int numWeights ) {
#if 1
    SMF_ASSERT( sizeof( CVertex ) == DRAWVERT_SIZE );
	SMF_ASSERT( (int)&((CVertex *)0)->xyz == DRAWVERT_XYZ_OFFSET );
	SMF_ASSERT( sizeof( CVec4D ) == JOINTWEIGHT_SIZE );
	SMF_ASSERT( sizeof( CMatJoint3x4 ) == JOINTMAT_SIZE );

	asm (
		"mov		eax, %1\n"   /*numVerts*/
		"test		eax, eax\n"
		"jz			done%=\n"
		"imul		eax, "DRAWVERT_SIZE_STR"\n"

		"mov		ecx, %0\n"  /*verts*/
		"mov		edx, %4\n"  /*index*/
		"mov		esi, %3\n"  /*weights*/
		"mov		edi, %2\n"  /*joints*/

		"add			ecx, eax\n"
		"neg			eax\n"

	"loopVert%=:\n"
		"mov			ebx, [edx]\n"
		"movups		xmm2, [esi]\n"
		"add			edx, 8\n"
		"movups		xmm0, xmm2\n"
		"add			esi, "JOINTWEIGHT_SIZE_STR"\n"
		"movups		xmm1, xmm2\n"

		"mulps		xmm0, [edi+ebx+ 0]\n"						// xmm0 = m0, m1, m2, t0
		"mulps		xmm1, [edi+ebx+16]\n"						// xmm1 = m3, m4, m5, t1
		"mulps		xmm2, [edi+ebx+32]\n"						// xmm2 = m6, m7, m8, t2

		"cmp			dword ptr [edx-4], 0\n"

		"jne			doneWeight%=\n"

	"loopWeight%=:\n"
		"mov		ebx, [edx]\n"
		"movups		xmm5, [esi]\n"
		"add		edx, 8\n"
		"movups		xmm3, xmm5\n"
		"add		esi, "JOINTWEIGHT_SIZE_STR"\n"
		"movups		xmm4, xmm5\n"

		"mulps		xmm3, [edi+ebx+ 0]\n"						// xmm3 = m0, m1, m2, t0
		"mulps		xmm4, [edi+ebx+16]\n"						// xmm4 = m3, m4, m5, t1
		"mulps		xmm5, [edi+ebx+32]\n"						// xmm5 = m6, m7, m8, t2

		"cmp			dword ptr [edx-4], 0\n"

		"addps		xmm0, xmm3\n"
		"addps		xmm1, xmm4\n"
		"addps		xmm2, xmm5\n"

		"je			loopWeight%=\n"

	"doneWeight%=:\n"
		"add			eax, "DRAWVERT_SIZE_STR"\n"

		haddps(		_xmm0, _xmm1 )
		haddps(		_xmm2, _xmm0 )

		"movhps		[ecx+eax-"DRAWVERT_SIZE_STR"+0], xmm2\n"

		haddps(		_xmm2, _xmm2 )

		"movss		[ecx+eax-"DRAWVERT_SIZE_STR"+8], xmm2\n"

		"jl			loopVert%=\n"
	"done%=:\n"
	:
    :"m"(verts), "m"(numVerts), "m"(joints), "m"(weights), "m"(index), "m"(numWeights)
    :);


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

#if !defined(MASM_INTEL)  //masm=at&t
//======================  CVec3D============================
#define  HIPREC
//Calculates the dot product (known as the inner product) of pSrc1 and pSrc2, that represents the cosine of the angle between them, if both are normalized.
float  VPCALL CSIMD_SSE3::vector3D_(const CVec3D* pSrc1, const CVec3D* pSrc2)
{
		/* SSE3 Implementation */

		asm(

		/*
			NB: Make sure that haddps does not add the 'w' component into equation.

			Do so with AND

			x & 0xFFFFFFFF -> x
			y & 0xFFFFFFFF -> y
			z & 0xFFFFFFFF -> z
			w & 0x00000000 -> 0.0f

		*/

		/* Load the vectors */
		"movups (%0), %%xmm0\n"
		"movups (%1), %%xmm1\n"

		/* Remove w component from one with AND mask */
		"andps %2, %%xmm0\n"
		"mulps %%xmm0, %%xmm1\n"

		/* Sum components */
		"haddps %%xmm1, %%xmm1\n"
		"haddps %%xmm1, %%xmm1\n"

		"movss %%xmm1, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		:
		: "r" (pSrc1), "r" (pSrc2),"m"(_SIMDx86_float_SSE_NO_W_MASK)
		);



}
float  VPCALL CSIMD_SSE3::vector4D_Dot(const CVec4D* pSrc1, const CVec4D* pSrc2)
{

	float dummy;
	/* SSE3 Implementation */
	asm(
	"movups (%1), %%xmm0\n"
	"movups (%2), %%xmm1\n"
	"mulps %%xmm0, %%xmm1\n"
	"haddps %%xmm1, %%xmm1\n"
	"haddps %%xmm1, %%xmm1\n"
	"movss %%xmm1, -4(%%esp)\n"
	"flds -4(%%esp)\n"
	: "=t" (dummy)
	: "r" (pSrc1), "r" (pSrc2)
	);
	return dummy;

}
//Returns the squared length of pSrc1. This can be useful when doing sphere-sphere collision tests where the exact distance is not needed and a squared one can be used, for example.
//It is signifigantly faster than using the square root version SIMDx86Vector_Length().

float  VPCALL CSIMD_SSE3::vector3D_LengthSq(const CVec3D* pVec)
{

		float dummy;


		/* SSE3 Implementation */
		asm(
		"movups (%1), %%xmm0\n"             /* xmm0 = x | y | z | w*/
		"andps %2, %%xmm0\n"    /* xmm0 = x | y | z | 0.0f */
		"mulps %%xmm0, %%xmm0\n"            /* xmm0 = x*x | y*y | z*z | 0.0f */
		"haddps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"           /* xmm0 = lensq | lensq | lensq | lensq */

		/* Store */
		"movss %%xmm0, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pVec),"m"(_SIMDx86_float_SSE_NO_W_MASK)
		);


		return dummy;

}
//Returns the length (magnitude) of pVec. Normalized vectors (e.g. the output of SIMDx86Vector_Normalize() ) are of unit length, or 1.0)
float  VPCALL CSIMD_SSE3::vector3D_Length(const CVec3D* pVec)
{
				float dummy;

		/* SSE3 Implementation */
		asm(
		"movups (%1), %%xmm0\n"
		"andps %2, %%xmm0\n"    /* xmm0 = x | y | z | 0.0f */
		"mulps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"

		#ifdef HIPREC
		/* Full square root */
		"sqrtss %%xmm0, %%xmm0\n"
		#else
		/* rcp( rsqrt(value) ) (This may be very inaccurate) */
		"rsqrtss %%xmm0, %%xmm0\n"
		"rcpss %%xmm0, %%xmm0\n"
		#endif


		/* Store */
		"movss %%xmm0, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pVec),"m"(_SIMDx86_float_SSE_NO_W_MASK)
		);

		return dummy;

}
void VPCALL CSIMD_SSE3::vector3D_Normalize(CVec3D* pVec)
{

		/* SSE3 Implementation */
		asm(
		"movups (%0), %%xmm0\n"
		"andps %1, %%xmm0\n"    /* xmm0 = x | y | z | 0.0f */
		"movaps %%xmm0, %%xmm2\n"
		"mulps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"

		/* High Precision SSE3 */
		#ifdef HIPREC

			/* Do full divide by the square root in parallel */
			"sqrtps %%xmm0, %%xmm0\n"				/* xmm0 = mag(pVec) */
			"divps %%xmm0, %%xmm2\n"				/* xmm2 = ? | z / mag | y / mag | x / mag */

		#else /* Low Precision SSE3 */

			/* Multiply by reciprocal square root */
			"rsqrtps %%xmm0, %%xmm1\n"				/* xmm1 = rcp(pVec) | rcp(pVec) | rcp(pVec) | rcp(pVec) */
			"mulps %%xmm1, %%xmm2\n"				/* xmm2 = norm(pVecv) */

		#endif

		/* Store*/
        "movlpd %%xmm2, (%0)\n"
        "shufps $0x02, %%xmm2,%%xmm2\n"  //shuffle colocar nos bits 0-31 o 3o byte
        "movss	%%xmm2,8(%0)\n"

		:
		: "r" (pVec),"m"(_SIMDx86_float_SSE_NO_W_MASK)
		);



}
void VPCALL CSIMD_SSE3::vector3D_NormalizeOf(CVec3D* pOut, const CVec3D* pVec)
{

		/* SSE3 Implementation */
		asm(
		"movups (%0), %%xmm0\n"
		"andps %2, %%xmm0\n"    /* xmm0 = x | y | z | 0.0f */
		"movaps %%xmm0, %%xmm2\n"
		"mulps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"

		/* High Precision SSE3 */
		#ifdef HIPREC

			/* Do full divide by the square root in parallel */
			"sqrtps %%xmm0, %%xmm0\n"				/* xmm0 = mag(pVec) */
			"divps %%xmm0, %%xmm2\n"				/* xmm2 = ? | z / mag | y / mag | x / mag */

		#else /* Low Precision SSE3 */

			/* Multiply by reciprocal square root */
			"rsqrtps %%xmm0, %%xmm1\n"				/* xmm1 = rcp(pVec) | rcp(pVec) | rcp(pVec) | rcp(pVec) */
			"mulps %%xmm1, %%xmm2\n"				/* xmm2 = norm(pVecv) */

		#endif

		/* Store*/
        "movlpd %%xmm2, (%1)\n"
        "shufps $0x02, %%xmm2,%%xmm2\n"  //shuffle colocar nos bits 0-31 o 3o byte
        "movss	%%xmm1,8(%1)\n"

		:
		: "r" (pVec), "r" (pOut),"m"(_SIMDx86_float_SSE_NO_W_MASK)
		);


}
float  VPCALL CSIMD_SSE3::vector3D_Distance(const CVec3D* pVec1, const CVec3D* pVec2)
{
	/* TODO: Optimize me completely */
	CVec3D diff;
	vector3D_DiffOf(&diff, pVec1, pVec2);
	return vector3D_Length(&diff);
}
#if 0
/*
	TODO:
	Optimize Distance() -- It currently uses already implemented functions
*/


float  VPCALL CSIMD_SSE3::vector3D_AlignedDot(const CVec3D* pSrc1, const CVec3D* pSrc2)
{

		/* SSE3 Implementation */
		#if USE_SSE >= 3

		asm(

		/*
			NB: Make sure that haddps does not add the 'w' component into equation.

			Do so with AND

			x & 0xFFFFFFFF -> x
			y & 0xFFFFFFFF -> y
			z & 0xFFFFFFFF -> z
			w & 0x00000000 -> 0.0f

		*/

		/* Load the vectors */
		"movaps (%0), %%xmm0\n"
		"mulps (%1), %%xmm0\n"

		/* Remove w component from one with AND mask */
		"andps _SIMDx86_float_SSE_NO_W_MASK, %%xmm0\n"

		/* Sum components */
		"haddps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"

		"movss %%xmm0, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		:
		: "r" (pSrc1), "r" (pSrc2)
		);


		#else
		/* SSE/SSE2 Implementation */
		float dummy;
		asm(

		/* Move vectors and multiply across */
		"movaps (%1), %%xmm0\n"
		"mulps (%2), %%xmm0\n"
		"andps _SIMDx86_float_SSE_NO_W_MASK, %%xmm0\n"

		/* Shift data around in the registers (I loves me some haddps right now!!) */
		"movhlps %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   |  z's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
		"shufps $0x55, %%xmm0, %%xmm0\n"/* xmm0 = y's | y's | y's | y's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */


		"movss %%xmm1, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pSrc1), "r" (pSrc2)
		);
		return dummy;
		#endif


}

float  VPCALL CSIMD_SSE3::vector3D_AlignedDot4(const float* pSrc4D1, const float* pSrc4D2)
{

	#if USE_SSE >= 3

	float dummy;
	/* SSE3 Implementation */
	asm(
	"movaps (%1), %%xmm0\n"
	"mulps (%2), %%xmm0\n"
	"haddps %%xmm0, %%xmm0\n"
	"haddps %%xmm0, %%xmm0\n"
	"movss %%xmm0, -4(%%esp)\n"
	"flds -4(%%esp)\n"
	: "=t" (dummy)
	: "r" (pSrc4D1), "r" (pSrc4D2)
	);
	return dummy;
	#else
	/* SSE/SSE2 Implementation */
	float dummy;
	asm(

	"movaps (%1), %%xmm0\n"
	"mulps (%2), %%xmm0\n"				/* xmm0 = w1*w2 | z1*z2 | y1*y2 | x1*x2 */
	"movaps %%xmm0, %%xmm1\n"			/* xmm1 = xmm0 */
	"shufps $0x1B, %%xmm1, %%xmm1\n"	/* xmm1 = x | y | z | w */
	"addps %%xmm0, %%xmm1\n"			/* xmm1 = w+x | y+z | y+z | w+x */
	"movss %%xmm1, %%xmm2\n"			/* xmm2 = ??? | ??? | ??? | w+x */
	"shufps $0x01, %%xmm1, %%xmm1\n"	/* xmm2 = ??? | ??? | ??? | y+z */
	"addss %%xmm1, %%xmm2\n"			/* xmm2 = ??? | ??? | ??? | dot4 */
	"movss %%xmm2, -4(%%esp)\n"
	"flds -4(%%esp)\n"
	: "=t" (dummy)
	: "r" (pSrc4D1), "r" (pSrc4D2)
	);
	return dummy;
	#endif


}

float  VPCALL CSIMD_SSE3::vector3D_AlignedLengthSq(const CVec3D* pVec)
{


		float dummy;

		#if (USE_SSE >= 3)

		/* SSE3 Implementation */
		asm(
		"movaps (%1), %%xmm0\n"             /* xmm0 = x | y | z | w*/
		"andps _SIMDx86_float_SSE_NO_W_MASK, %%xmm0\n"    /* xmm0 = x | y | z | 0.0f */
		"mulps %%xmm0, %%xmm0\n"            /* xmm0 = x*x | y*y | z*z | 0.0f */
		"haddps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"           /* xmm0 = lensq | lensq | lensq | lensq */

		/* Store */
		"movss %%xmm0, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pVec)
		);

		#else /* Using SSE/SSE2 */

		/* SSE/SSE2 Implementation */
		asm(
		"movaps (%1), %%xmm0\n"
		"andps _SIMDx86_float_SSE_NO_W_MASK, %%xmm0\n"    /* xmm0 = x | y | z | 0.0f */
		"mulps %%xmm0, %%xmm0\n"



		/* Shift data around in the registers (I loves me some haddps right now!!) */
		"movhlps %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   |  z's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
		"shufps $0x55, %%xmm0, %%xmm0\n"/* xmm0 = y's | y's | y's | y's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */


		"movss %%xmm1, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pVec)
		);


		#endif /* USE_SSE == 1 ||  USE_SSE == 2 */

		return dummy;

}

float  VPCALL CSIMD_SSE3::vector3D_AlignedLength(const CVec3D* pVec)
{


		#if USE_SSE >= 3
		float dummy;

		/* SSE3 Implementation */
		asm(
		"movaps (%1), %%xmm0\n"
		"andps _SIMDx86_float_SSE_NO_W_MASK, %%xmm0\n"    /* xmm0 = x | y | z | 0.0f */
		"mulps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"

		#ifdef HIPREC
		/* Full square root */
		"sqrtss %%xmm0, %%xmm0\n"
		#else
		/* rcp( rsqrt(value) ) (This may be very inaccurate) */
		"rsqrtss %%xmm0, %%xmm0\n"
		"rcpss %%xmm0, %%xmm0\n"
		#endif

		/* Store */
		"movss %%xmm0, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pVec)
		);

		return dummy;
		#else

		float dummy;
		/* SSE/SSE2 Implementation */
		asm(
		"movaps (%1), %%xmm0\n"
		"andps _SIMDx86_float_SSE_NO_W_MASK, %%xmm0\n"    /* xmm0 = x | y | z | 0.0f */
		"mulps %%xmm0, %%xmm0\n"

		/* Shift data around in the registers (I loves me some haddps right now!!) */
		"movhlps %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   |  z's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
		"shufps $0x55, %%xmm0, %%xmm0\n"/* xmm0 = y's | y's | y's | y's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */


		#ifdef HIPREC
		/* Full square root */
		"sqrtss %%xmm1, %%xmm1\n"
		#else
		/* rcp( rsqrt(value) ) (This may be very inaccurate) */
		"rsqrtss %%xmm1, %%xmm1\n"
		"rcpss %%xmm1, %%xmm1\n"
		#endif

		"movss %%xmm1, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pVec)
		);

		return dummy;
		#endif

}

void VPCALL CSIMD_SSE3::vector3D_AlignedCross(CVec3D* pLeft, const CVec3D* pRight)
{

	asm(
	"movaps (%0), %%xmm0\n"		/* xmm0 = pLeft */
	"movaps (%1), %%xmm1\n"		/* xmm1 = pRight */
	"movaps %%xmm0, %%xmm2\n"	/* xmm2 = pLeft */
	"movaps %%xmm1, %%xmm3\n"	/* xmm3 = pRight */

	/*
		Shuffles
		w | x | z | y  == 1100 1001 = 0xC9
		w | y | x | z  == 1101 0010 = 0xD2
	*/
	"shufps $0xC9, %%xmm0, %%xmm0\n"
	"shufps $0xD2, %%xmm1, %%xmm1\n"
	"shufps $0xD2, %%xmm2, %%xmm2\n"
	"shufps $0xC9, %%xmm3, %%xmm3\n"

	/* Multiply columns 1&2 and 3&4 */
	"mulps %%xmm0, %%xmm1\n"
	"mulps %%xmm2, %%xmm3\n"

	/* Subtract products to get the cross product! */
	"subps %%xmm3, %%xmm1\n"

	/* Store */
	"movaps %%xmm1, (%0)\n"
	:
	: "r" (pLeft), "r" (pRight)
	);


}

void VPCALL CSIMD_SSE3::vector3D_AlignedCrossOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight)
{

	asm(
	"movaps (%1), %%xmm0\n"		/* xmm0 = pLeft */
	"movaps (%2), %%xmm1\n"		/* xmm1 = pRight */
	"movaps %%xmm0, %%xmm2\n"	/* xmm2 = pLeft */
	"movaps %%xmm1, %%xmm3\n"	/* xmm3 = pRight */

	/*
		Shuffles
		w | x | z | y  == 1100 1001 = 0xC9
		w | y | x | z  == 1101 0010 = 0xD2
	*/
	"shufps $0xC9, %%xmm0, %%xmm0\n"
	"shufps $0xD2, %%xmm1, %%xmm1\n"
	"shufps $0xD2, %%xmm2, %%xmm2\n"
	"shufps $0xC9, %%xmm3, %%xmm3\n"

	/* Multiply columns 1&2 and 3&4 */
	"mulps %%xmm0, %%xmm1\n"
	"mulps %%xmm2, %%xmm3\n"

	/* Subtract products to get the cross product! */
	"subps %%xmm3, %%xmm1\n"
	/* Store */
	"movaps %%xmm1, (%0)\n"
	:
	: "r" (pOut), "r" (pLeft), "r" (pRight)
	);


}
void VPCALL CSIMD_SSE3::vector3D_AlignedNormalize(CVec3D* pVec)
{

		#if USE_SSE >= 3

		/* SSE3 Implementation */
		asm(
		"movaps (%0), %%xmm0\n"
		"andps _SIMDx86_float_SSE_NO_W_MASK, %%xmm0\n"    /* xmm0 = x | y | z | 0.0f */
		"movaps %%xmm0, %%xmm2\n"
		"mulps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"

		/* High Precision SSE3 */
		#ifdef HIPREC

			/* Do full divide by the square root in parallel */
			"sqrtps %%xmm0, %%xmm0\n"				/* xmm0 = mag(pVec) */
			"divps %%xmm0, %%xmm2\n"				/* xmm2 = ? | z / mag | y / mag | x / mag */

		#else /* Low Precision SSE3 */

			/* Multiply by reciprocal square root */
			"rsqrtps %%xmm0, %%xmm1\n"				/* xmm1 = rcp(pVec) | rcp(pVec) | rcp(pVec) | rcp(pVec) */
			"mulps %%xmm1, %%xmm2\n"				/* xmm2 = norm(pVecv) */

		#endif

		/* Store */
		"movaps %%xmm2, (%0)\n"

		:
		: "r" (pVec)
		);
		#else

		/* SSE/SSE2 Implementation */
		asm(
		"movaps (%0), %%xmm0\n"
		"andps _SIMDx86_float_SSE_NO_W_MASK, %%xmm0\n"    /* xmm0 = x | y | z | 0.0f */
		"movaps %%xmm0, %%xmm7\n"			/* Save for the division by length */
		"mulps %%xmm0, %%xmm0\n"

		/* Shift data around in the registers (I loves me some haddps right now!!) */
		"movhlps %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   |  z's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
		"shufps $0x55, %%xmm0, %%xmm0\n"/* xmm0 = y's | y's | y's | y's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */

		#ifdef HIPREC

			/* Divide by magnitude */
			"sqrtss %%xmm1, %%xmm1\n"			/* xmm3 = ? | ? | ? | mag(pVec) */
			"shufps $0x00, %%xmm1, %%xmm1\n"	/* xmm3 = mag(pVec) | mag(pVec) | mag(pVec) | mag(pVec) */
			"divps %%xmm1, %%xmm7\n"			/* xmm1 = ? | z/mag | y/mag | x/mag */

		#else

			/* Multiply by reciprocal magnitude */
			"rsqrtss %%xmm1 ,%%xmm1\n"			/* xmm3 = ? | ? | ? | rcp(mag) */
			"shufps $0x00, %%xmm1, %%xmm1\n"	/* xmm3 = rcp(mag) | rcp(mag) | rcp(mag) | rcp(mag) */
			"mulps %%xmm1, %%xmm7\n"			/* xmm1 = ? | z/mag | y/mag | x/mag */

		#endif
		/* Store */
		"movaps %%xmm7, (%0)\n"
				:
		: "r" (pVec)
		);


		#endif


}
void VPCALL CSIMD_SSE3::vector3D_AlignedNormalizeOf(CVec3D* pOut, const CVec3D* pVec)
{

		#if USE_SSE >= 3

		/* SSE3 Implementation */
		asm(
		"movaps (%0), %%xmm0\n"
		"andps _SIMDx86_float_SSE_NO_W_MASK, %%xmm0\n"    /* xmm0 = x | y | z | 0.0f */
		"movaps %%xmm0, %%xmm2\n"
		"mulps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"

		/* High Precision SSE3 */
		#ifdef HIPREC

			/* Do full divide by the square root in parallel */
			"sqrtps %%xmm0, %%xmm0\n"				/* xmm0 = mag(pVec) */
			"divps %%xmm0, %%xmm2\n"				/* xmm2 = ? | z / mag | y / mag | x / mag */

		#else /* Low Precision SSE3 */

			/* Multiply by reciprocal square root */
			"rsqrtps %%xmm0, %%xmm1\n"				/* xmm1 = rcp(pVec) | rcp(pVec) | rcp(pVec) | rcp(pVec) */
			"mulps %%xmm1, %%xmm2\n"				/* xmm2 = norm(pVecv) */

		#endif

		/* Store */
		"movaps %%xmm2, (%1)\n"

		:
		: "r" (pVec), "r" (pOut)
		);
		#else /* SSE/SSE2 Implementation */
		asm(
		"movaps (%0), %%xmm0\n"
		"andps _SIMDx86_float_SSE_NO_W_MASK, %%xmm0\n"    /* xmm0 = x | y | z | 0.0f */
		"movaps %%xmm0, %%xmm7\n"			/* Save for the division by length */
		"mulps %%xmm0, %%xmm0\n"

		/* Shift data around in the registers (I loves me some haddps right now!!) */
		"movhlps %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   |  z's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
		"shufps $0x55, %%xmm0, %%xmm0\n"/* xmm0 = y's | y's | y's | y's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */

		#ifdef HIPREC

			/* Divide by magnitude */
			"sqrtss %%xmm1, %%xmm1\n"			/* xmm3 = ? | ? | ? | mag(pVec) */
			"shufps $0x00, %%xmm1, %%xmm1\n"	/* xmm3 = mag(pVec) | mag(pVec) | mag(pVec) | mag(pVec) */
			"divps %%xmm1, %%xmm7\n"			/* xmm1 = ? | z/mag | y/mag | x/mag */

		#else

			/* Multiply by reciprocal magnitude */
			"rsqrtss %%xmm1 ,%%xmm1\n"			/* xmm3 = ? | ? | ? | rcp(mag) */
			"shufps $0x00, %%xmm1, %%xmm1\n"	/* xmm3 = rcp(mag) | rcp(mag) | rcp(mag) | rcp(mag) */
			"mulps %%xmm1, %%xmm7\n"			/* xmm1 = ? | z/mag | y/mag | x/mag */

		#endif
		/* Store */
		"movaps %%xmm7, (%1)\n"
				:
		: "r" (pVec), "r" (pOut)
		);

	#endif /* USE_SSE == 1 || USE_SSE == 2 */


}

float  VPCALL CSIMD_SSE3::vector3D_AlignedDistance(const CVec3D* pVec1, const CVec3D* pVec2)
{
#if 0 /* Not so good just yet... */
	#if defined(USE_SSE)

	#if USE_SSE == 3
		float dummy;
		asm(
			"movaps (%1), %%xmm0\n"
			"subps (%2), %%xmm0\n"
			"andps _SIMDx86_float_SSE_NO_W_MASK, %%xmm0\n"
			"haddps %%xmm0, %%xmm0\n"
			"haddps %%xmm0, %%xmm0\n"
			"subl $4, %%esp\n"
			"movss %%xmm0, (%%esp)\n"
			"flds (%%esp)\n"
			"addl $4, %%esp\n"
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

#endif //if 0

//=================== CVec4D======================


float  VPCALL CSIMD_SSE::vector4D_Dot(const CVec4D* pSrc4D1, const CVec4D* pSrc4D2)
{

	#if USE_SSE >= 3

	float dummy;
	/* SSE3 Implementation */
	asm(
	"movups (%1), %%xmm0\n"
	"movups (%2), %%xmm1\n"
	"mulps %%xmm0, %%xmm1\n"
	"haddps %%xmm1, %%xmm1\n"
	"haddps %%xmm1, %%xmm1\n"
	"movss %%xmm1, -4(%%esp)\n"
	"flds -4(%%esp)\n"
	: "=t" (dummy)
	: "r" (pSrc4D1), "r" (pSrc4D2)
	);
	return dummy;
	#else
	/* SSE/SSE2 Implementation */
	float dummy;
	asm(
	"movups (%1), %%xmm0\n"
	"movups (%2), %%xmm1\n"
	"mulps %%xmm0, %%xmm1\n"
	"movaps %%xmm1, %%xmm2\n"
	"shufps $0x1B, %%xmm2, %%xmm2\n"	/* xmm2 = x | y | z | w */
	"addps %%xmm1, %%xmm2\n"			/* xmm2 = w+x | y+z | y+z | w+x */
	"movss %%xmm2, %%xmm3\n"				/* xmm3 = ??? | ??? | ??? | w+x */
	"shufps $0x01, %%xmm2, %%xmm2\n"			/* xmm2 = ??? | ??? | ??? | y+z */
	"addss %%xmm3, %%xmm2\n"			/* xmm2 = ??? | ??? | ??? | dot4 */
	"movss %%xmm2, -4(%%esp)\n"
	"flds -4(%%esp)\n"
	: "=t" (dummy)
	: "r" (pSrc4D1), "r" (pSrc4D2)
	);
	return dummy;
	#endif


}

float  VPCALL CSIMD_SSE::vector4D_LengthSq(const CVec4D* pVec)
{

		float dummy;

		#if (USE_SSE >= 3)

		/* SSE3 Implementation */
		asm(
		"movups (%1), %%xmm0\n"             /* xmm0 = x | y | z | w*/
		"mulps %%xmm0, %%xmm0\n"            /* xmm0 = x*x | y*y | z*z | 0.0f */
		"haddps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"           /* xmm0 = lensq | lensq | lensq | lensq */

		/* Store */
		"movss %%xmm0, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pVec)
		);

		#else /* Using SSE/SSE2 */

		/* SSE/SSE2 Implementation */
		asm(
		"movups (%1), %%xmm0\n"
		"mulps %%xmm0, %%xmm0\n"

		/* Shift data around in the registers (I loves me some haddps right now!!) */
		"movhlps %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   |  z's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
		"shufps $0x55, %%xmm0, %%xmm0\n"/* xmm0 = y's | y's | y's | y's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */

		"movss %%xmm1, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pVec)
		);


		#endif /* USE_SSE == 1 ||  USE_SSE == 2 */

		return dummy;

}

float  VPCALL CSIMD_SSE::vector4D_Length(const CVec4D* pVec)
{

		#if USE_SSE >= 3
		float dummy;

		/* SSE3 Implementation */
		asm(
		"movups (%1), %%xmm0\n"
		"mulps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"

		#ifdef HIPREC
		/* Full square root */
		"sqrtss %%xmm0, %%xmm0\n"
		#else
		/* rcp( rsqrt(value) ) (This may be very inaccurate) */
		"rsqrtss %%xmm0, %%xmm0\n"
		"rcpss %%xmm0, %%xmm0\n"
		#endif


		/* Store */
		"movss %%xmm0, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pVec)
		);

		return dummy;
		#else

		float dummy;
		/* SSE/SSE2 Implementation */
		asm(
		"movups (%1), %%xmm0\n"
		"mulps %%xmm0, %%xmm0\n"

		/* Shift data around in the registers (I loves me some haddps right now!!) */
		"movhlps %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   |  z's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
		"shufps $0x55, %%xmm0, %%xmm0\n"/* xmm0 = y's | y's | y's | y's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */

		#ifdef HIPREC
		/* Full square root */
		"sqrtss %%xmm1, %%xmm1\n"
		#else
		/* rcp( rsqrt(value) ) (This may be very inaccurate) */
		"rsqrtss %%xmm1, %%xmm1\n"
		"rcpss %%xmm1, %%xmm1\n"
		#endif

		"movss %%xmm1, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pVec)
		);

		return dummy;
		#endif

}


void VPCALL CSIMD_SSE::vector4D_Normalize(CVec4D* pVec)
{

		#if USE_SSE >= 3

		/* SSE3 Implementation */
		asm(
		"movups (%0), %%xmm0\n"
		"movaps %%xmm0, %%xmm2\n"
		"mulps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"

		/* High Precision SSE3 */
		#ifdef HIPREC

			/* Do full divide by the square root in parallel */
			"sqrtps %%xmm0, %%xmm0\n"				/* xmm0 = mag(pVec) */
			"divps %%xmm0, %%xmm2\n"				/* xmm2 = ? | z / mag | y / mag | x / mag */

		#else /* Low Precision SSE3 */

			/* Multiply by reciprocal square root */
			"rsqrtps %%xmm0, %%xmm1\n"				/* xmm1 = rcp(pVec) | rcp(pVec) | rcp(pVec) | rcp(pVec) */
			"mulps %%xmm1, %%xmm2\n"				/* xmm2 = norm(pVecv) */

		#endif

		/* Store */
		"movups %%xmm2, (%0)\n"

		:
		: "r" (pVec)
		);
		#else

		/* SSE/SSE2 Implementation */
		asm(
		"movups (%0), %%xmm0\n"
		"movaps %%xmm0, %%xmm7\n"			/* Save for the division by length */
		"mulps %%xmm0, %%xmm0\n"

		/* Shift data around in the registers (I loves me some haddps right now!!) */
		"movhlps %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   |  z's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
		"shufps $0x55, %%xmm0, %%xmm0\n"/* xmm0 = y's | y's | y's | y's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */

		#ifdef HIPREC

			/* Divide by magnitude */
			"sqrtss %%xmm1, %%xmm1\n"			/* xmm3 = ? | ? | ? | mag(pVec) */
			"shufps $0x00, %%xmm1, %%xmm1\n"	/* xmm3 = mag(pVec) | mag(pVec) | mag(pVec) | mag(pVec) */
			"divps %%xmm1, %%xmm7\n"			/* xmm1 = ? | z/mag | y/mag | x/mag */

		#else

			/* Multiply by reciprocal magnitude */
			"rsqrtss %%xmm1 ,%%xmm1\n"			/* xmm3 = ? | ? | ? | rcp(mag) */
			"shufps $0x00, %%xmm1, %%xmm1\n"	/* xmm3 = rcp(mag) | rcp(mag) | rcp(mag) | rcp(mag) */
			"mulps %%xmm1, %%xmm7\n"			/* xmm1 = ? | z/mag | y/mag | x/mag */

		#endif
		/* Store */
		"movups %%xmm7, (%0)\n"
				:
		: "r" (pVec)
		);


		#endif


}
void VPCALL CSIMD_SSE::vector4D_NormalizeOf(CVec4D* pOut, const CVec4D* pVec)
{

		#if USE_SSE >= 3

		/* SSE3 Implementation */
		asm(
		"movups (%0), %%xmm0\n"
		"movaps %%xmm0, %%xmm2\n"
		"mulps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"

		/* High Precision SSE3 */
		#ifdef HIPREC

			/* Do full divide by the square root in parallel */
			"sqrtps %%xmm0, %%xmm0\n"				/* xmm0 = mag(pVec) */
			"divps %%xmm0, %%xmm2\n"				/* xmm2 = ? | z / mag | y / mag | x / mag */

		#else /* Low Precision SSE3 */

			/* Multiply by reciprocal square root */
			"rsqrtps %%xmm0, %%xmm1\n"				/* xmm1 = rcp(pVec) | rcp(pVec) | rcp(pVec) | rcp(pVec) */
			"mulps %%xmm1, %%xmm2\n"				/* xmm2 = norm(pVecv) */

		#endif

		/* Store */
		"movups %%xmm2, (%1)\n"

		:
		: "r" (pVec), "r" (pOut)
		);
		#else /* SSE/SSE2 Implementation */
		asm(
		"movups (%0), %%xmm0\n"
		"movaps %%xmm0, %%xmm7\n"			/* Save for the division by length */
		"mulps %%xmm0, %%xmm0\n"

		/* Shift data around in the registers (I loves me some haddps right now!!) */
		"movhlps %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   |  z's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
		"shufps $0x55, %%xmm0, %%xmm0\n"/* xmm0 = y's | y's | y's | y's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */

		#ifdef HIPREC

			/* Divide by magnitude */
			"sqrtss %%xmm1, %%xmm1\n"			/* xmm3 = ? | ? | ? | mag(pVec) */
			"shufps $0x00, %%xmm1, %%xmm1\n"	/* xmm3 = mag(pVec) | mag(pVec) | mag(pVec) | mag(pVec) */
			"divps %%xmm1, %%xmm7\n"			/* xmm1 = ? | z/mag | y/mag | x/mag */

		#else

			/* Multiply by reciprocal magnitude */
			"rsqrtss %%xmm1 ,%%xmm1\n"			/* xmm3 = ? | ? | ? | rcp(mag) */
			"shufps $0x00, %%xmm1, %%xmm1\n"	/* xmm3 = rcp(mag) | rcp(mag) | rcp(mag) | rcp(mag) */
			"mulps %%xmm1, %%xmm7\n"			/* xmm1 = ? | z/mag | y/mag | x/mag */

		#endif
		/* Store */
		"movups %%xmm7, (%1)\n"
				:
		: "r" (pVec), "r" (pOut)
		);

	#endif /* USE_SSE == 1 || USE_SSE == 2 */


}

float  VPCALL CSIMD_SSE::vector4D_Distance(const CVec4D* pVec1, const CVec4D* pVec2)
{
	/* TODO: Optimize me completely */
	CVec4D diff;
	vector4D_DiffOf(&diff, pVec1, pVec2);
	return vector4D_Length(&diff);
}


float  VPCALL CSIMD_SSE::vector4D_AlignedDot(const CVec4D* pSrc4D1, const CVec4D* pSrc4D2)
{

	#if USE_SSE >= 3

	float dummy;
	/* SSE3 Implementation */
	asm(
	"movaps (%1), %%xmm0\n"
	"mulps (%2), %%xmm0\n"
	"haddps %%xmm0, %%xmm0\n"
	"haddps %%xmm0, %%xmm0\n"
	"movss %%xmm0, -4(%%esp)\n"
	"flds -4(%%esp)\n"
	: "=t" (dummy)
	: "r" (pSrc4D1), "r" (pSrc4D2)
	);
	return dummy;
	#else
	/* SSE/SSE2 Implementation */
	float dummy;
	asm(

	"movaps (%1), %%xmm0\n"
	"mulps (%2), %%xmm0\n"				/* xmm0 = w1*w2 | z1*z2 | y1*y2 | x1*x2 */
	"movaps %%xmm0, %%xmm1\n"			/* xmm1 = xmm0 */
	"shufps $0x1B, %%xmm1, %%xmm1\n"	/* xmm1 = x | y | z | w */
	"addps %%xmm0, %%xmm1\n"			/* xmm1 = w+x | y+z | y+z | w+x */
	"movss %%xmm1, %%xmm2\n"			/* xmm2 = ??? | ??? | ??? | w+x */
	"shufps $0x01, %%xmm1, %%xmm1\n"	/* xmm2 = ??? | ??? | ??? | y+z */
	"addss %%xmm1, %%xmm2\n"			/* xmm2 = ??? | ??? | ??? | dot4 */
	"movss %%xmm2, -4(%%esp)\n"
	"flds -4(%%esp)\n"
	: "=t" (dummy)
	: "r" (pSrc4D1), "r" (pSrc4D2)
	);
	return dummy;
	#endif


}

float  VPCALL CSIMD_SSE::vector4D_AlignedLengthSq(const CVec4D* pVec)
{


		float dummy;

		#if (USE_SSE >= 3)

		/* SSE3 Implementation */
		asm(
		"movaps (%1), %%xmm0\n"             /* xmm0 = x | y | z | w*/
		"mulps %%xmm0, %%xmm0\n"            /* xmm0 = x*x | y*y | z*z | 0.0f */
		"haddps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"           /* xmm0 = lensq | lensq | lensq | lensq */

		/* Store */
		"movss %%xmm0, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pVec)
		);

		#else /* Using SSE/SSE2 */

		/* SSE/SSE2 Implementation */
		asm(
		"movaps (%1), %%xmm0\n"
		"mulps %%xmm0, %%xmm0\n"



		/* Shift data around in the registers (I loves me some haddps right now!!) */
		"movhlps %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   |  z's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
		"shufps $0x55, %%xmm0, %%xmm0\n"/* xmm0 = y's | y's | y's | y's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */


		"movss %%xmm1, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pVec)
		);


		#endif /* USE_SSE == 1 ||  USE_SSE == 2 */

		return dummy;

}

float  VPCALL CSIMD_SSE::vector4D_AlignedLength(const CVec4D* pVec)
{


		#if USE_SSE >= 3
		float dummy;

		/* SSE3 Implementation */
		asm(
		"movaps (%1), %%xmm0\n"
		"mulps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"

		#ifdef HIPREC
		/* Full square root */
		"sqrtss %%xmm0, %%xmm0\n"
		#else
		/* rcp( rsqrt(value) ) (This may be very inaccurate) */
		"rsqrtss %%xmm0, %%xmm0\n"
		"rcpss %%xmm0, %%xmm0\n"
		#endif

		/* Store */
		"movss %%xmm0, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pVec)
		);

		return dummy;
		#else

		float dummy;
		/* SSE/SSE2 Implementation */
		asm(
		"movaps (%1), %%xmm0\n"
		"mulps %%xmm0, %%xmm0\n"

		/* Shift data around in the registers (I loves me some haddps right now!!) */
		"movhlps %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   |  z's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
		"shufps $0x55, %%xmm0, %%xmm0\n"/* xmm0 = y's | y's | y's | y's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */


		#ifdef HIPREC
		/* Full square root */
		"sqrtss %%xmm1, %%xmm1\n"
		#else
		/* rcp( rsqrt(value) ) (This may be very inaccurate) */
		"rsqrtss %%xmm1, %%xmm1\n"
		"rcpss %%xmm1, %%xmm1\n"
		#endif

		"movss %%xmm1, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pVec)
		);

		return dummy;
		#endif

}


void VPCALL CSIMD_SSE::vector4D_AlignedNormalize(CVec4D* pVec)
{

		#if USE_SSE >= 3

		/* SSE3 Implementation */
		asm(
		"movaps (%0), %%xmm0\n"
		"movaps %%xmm0, %%xmm2\n"
		"mulps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"

		/* High Precision SSE3 */
		#ifdef HIPREC

			/* Do full divide by the square root in parallel */
			"sqrtps %%xmm0, %%xmm0\n"				/* xmm0 = mag(pVec) */
			"divps %%xmm0, %%xmm2\n"				/* xmm2 = ? | z / mag | y / mag | x / mag */

		#else /* Low Precision SSE3 */

			/* Multiply by reciprocal square root */
			"rsqrtps %%xmm0, %%xmm1\n"				/* xmm1 = rcp(pVec) | rcp(pVec) | rcp(pVec) | rcp(pVec) */
			"mulps %%xmm1, %%xmm2\n"				/* xmm2 = norm(pVecv) */

		#endif

		/* Store */
		"movaps %%xmm2, (%0)\n"

		:
		: "r" (pVec)
		);
		#else

		/* SSE/SSE2 Implementation */
		asm(
		"movaps (%0), %%xmm0\n"
		"movaps %%xmm0, %%xmm7\n"			/* Save for the division by length */
		"mulps %%xmm0, %%xmm0\n"

		/* Shift data around in the registers (I loves me some haddps right now!!) */
		"movhlps %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   |  z's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
		"shufps $0x55, %%xmm0, %%xmm0\n"/* xmm0 = y's | y's | y's | y's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */

		#ifdef HIPREC

			/* Divide by magnitude */
			"sqrtss %%xmm1, %%xmm1\n"			/* xmm3 = ? | ? | ? | mag(pVec) */
			"shufps $0x00, %%xmm1, %%xmm1\n"	/* xmm3 = mag(pVec) | mag(pVec) | mag(pVec) | mag(pVec) */
			"divps %%xmm1, %%xmm7\n"			/* xmm1 = ? | z/mag | y/mag | x/mag */

		#else

			/* Multiply by reciprocal magnitude */
			"rsqrtss %%xmm1 ,%%xmm1\n"			/* xmm3 = ? | ? | ? | rcp(mag) */
			"shufps $0x00, %%xmm1, %%xmm1\n"	/* xmm3 = rcp(mag) | rcp(mag) | rcp(mag) | rcp(mag) */
			"mulps %%xmm1, %%xmm7\n"			/* xmm1 = ? | z/mag | y/mag | x/mag */

		#endif
		/* Store */
		"movaps %%xmm7, (%0)\n"
				:
		: "r" (pVec)
		);


		#endif


}
void VPCALL CSIMD_SSE::vector4D_AlignedNormalizeOf(CVec4D* pOut, const CVec4D* pVec)
{

		#if USE_SSE >= 3

		/* SSE3 Implementation */
		asm(
		"movaps (%0), %%xmm0\n"
		"movaps %%xmm0, %%xmm2\n"
		"mulps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"

		/* High Precision SSE3 */
		#ifdef HIPREC

			/* Do full divide by the square root in parallel */
			"sqrtps %%xmm0, %%xmm0\n"				/* xmm0 = mag(pVec) */
			"divps %%xmm0, %%xmm2\n"				/* xmm2 = ? | z / mag | y / mag | x / mag */

		#else /* Low Precision SSE3 */

			/* Multiply by reciprocal square root */
			"rsqrtps %%xmm0, %%xmm1\n"				/* xmm1 = rcp(pVec) | rcp(pVec) | rcp(pVec) | rcp(pVec) */
			"mulps %%xmm1, %%xmm2\n"				/* xmm2 = norm(pVecv) */

		#endif

		/* Store */
		"movaps %%xmm2, (%1)\n"

		:
		: "r" (pVec), "r" (pOut)
		);
		#else /* SSE/SSE2 Implementation */
		asm(
		"movaps (%0), %%xmm0\n"
		"movaps %%xmm0, %%xmm7\n"			/* Save for the division by length */
		"mulps %%xmm0, %%xmm0\n"

		/* Shift data around in the registers (I loves me some haddps right now!!) */
		"movhlps %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   |  z's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
		"shufps $0x55, %%xmm0, %%xmm0\n"/* xmm0 = y's | y's | y's | y's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */

		#ifdef HIPREC

			/* Divide by magnitude */
			"sqrtss %%xmm1, %%xmm1\n"			/* xmm3 = ? | ? | ? | mag(pVec) */
			"shufps $0x00, %%xmm1, %%xmm1\n"	/* xmm3 = mag(pVec) | mag(pVec) | mag(pVec) | mag(pVec) */
			"divps %%xmm1, %%xmm7\n"			/* xmm1 = ? | z/mag | y/mag | x/mag */

		#else

			/* Multiply by reciprocal magnitude */
			"rsqrtss %%xmm1 ,%%xmm1\n"			/* xmm3 = ? | ? | ? | rcp(mag) */
			"shufps $0x00, %%xmm1, %%xmm1\n"	/* xmm3 = rcp(mag) | rcp(mag) | rcp(mag) | rcp(mag) */
			"mulps %%xmm1, %%xmm7\n"			/* xmm1 = ? | z/mag | y/mag | x/mag */

		#endif
		/* Store */
		"movaps %%xmm7, (%1)\n"
				:
		: "r" (pVec), "r" (pOut)
		);

	#endif /* USE_SSE == 1 || USE_SSE == 2 */


}

float  VPCALL CSIMD_SSE::vector4D_AlignedDistance(const CVec4D* pVec1, const CVec4D* pVec2)
{
#if 0 /* Not so good just yet... */
	#if defined(USE_SSE)

	#if USE_SSE == 3
		float dummy;
		asm(
			"movaps (%1), %%xmm0\n"
			"subps (%2), %%xmm0\n"
			"haddps %%xmm0, %%xmm0\n"
			"haddps %%xmm0, %%xmm0\n"
			"subl $4, %%esp\n"
			"movss %%xmm0, (%%esp)\n"
			"flds (%%esp)\n"
			"addl $4, %%esp\n"
			: "=t" (dummy)
			: "r" (pVec1), "r" (pVec2)
			);
		return dummy;
	#else
	#endif
	#endif

#endif

	/* TODO: Optimize me completely */
	CVec4D diff;
	vector4D_DiffOf(&diff, pVec1, pVec2);
	return vector4D_Length(&diff);

}

//==============END CVec4D============================================

//=====================Quaternion===================================

void  VPCALL CSIMD_SSE3::quaternion_Normalize(CQuaternion* pQuat)
{
		asm(
		"movups (%0), %%xmm0\n"			/* xmm0 = w | z | y | x  */
		"movaps %%xmm0, %%xmm1\n"		/* xmm1 = xmm0 */
		"mulps %%xmm0, %%xmm0\n"		/* xmm0 = w*w | z*z | y*y | x*x */
		"haddps %%xmm0, %%xmm0\n"		/* xmm0 = */
		"haddps %%xmm0, %%xmm0\n"		/* xmm0 = magsq | magsq | magsq | magsq */
		"rsqrtps %%xmm0, %%xmm0\n"		/* xmm0 = rcp(mag) | rcp(mag) | rcp(mag) | rcp(mag) */
		"mulps %%xmm0, %%xmm1\n"
		"movups %%xmm1, (%0)\n"
		:
		: "r" (pQuat)
		);

}

void  VPCALL CSIMD_SSE3::quaternion_NormalizeOf(CQuaternion* pOut, const CQuaternion* pQuat)
{

		asm(
		"movups (%0), %%xmm0\n"			/* xmm0 = w | z | y | x  */
		"movaps %%xmm0, %%xmm1\n"		/* xmm1 = xmm0 */
		"mulps %%xmm0, %%xmm0\n"		/* xmm0 = w*w | z*z | y*y | x*x */
		"haddps %%xmm0, %%xmm0\n"		/* xmm0 = */
		"haddps %%xmm0, %%xmm0\n"		/* xmm0 = magsq | magsq | magsq | magsq */
		"rsqrtps %%xmm0, %%xmm0\n"		/* xmm0 = rcp(mag) | rcp(mag) | rcp(mag) | rcp(mag) */
		"mulps %%xmm0, %%xmm1\n"
		"movups %%xmm1, (%1)\n"
		:
		: "r" (pQuat), "r" (pOut)
		);


}

#else //masm = intel
//============CVec3D===================================



//======================  CVec3D============================
#define  HIPREC
//Calculates the dot product (known as the inner product) of pSrc1 and pSrc2, that represents the cosine of the angle between them, if both are normalized.
float  VPCALL CSIMD_SSE3::vector3D_Dot(const CVec3D* pSrc1, const CVec3D* pSrc2)
{
	    float dummy;
		/* SSE3 Implementation */
		float *min = &dummy;
		asm
        (
			"mov eax, %0\n"                       // Load pointers into CPU regs and multiply
			"mov ebx, %1\n"
			"movups xmm0, [eax]\n"
			"movups xmm1, [ebx]\n"
			"andps xmm0, %2\n"  //remove w component /* xmm0 = x | y | z | 0.0f */
			"andps xmm1, %2\n"  //remove w component /* xmm1 = x | y | z | 0.0f */
			"mulps xmm1, xmm0\n"
			/* Sum components */
			"haddps xmm1, xmm1\n"
			"haddps xmm1, xmm1\n"

			//coloca 0-31 de xmm2 em dummy
			"mov		esi, %3\n"
			"movss		[esi], xmm1\n"
		:
		: "r" (pSrc1), "r" (pSrc2), "m" (_SIMDx86_float_SSE_NO_W_MASK),"m" (min)
		);


		return dummy;
}
float  VPCALL CSIMD_SSE3::vector4D_Dot(const CVec4D* pSrc1, const CVec4D* pSrc2)
{
	float dummy;

			/* SSE3 Implementation */
		float *min = &dummy;
		asm
        (
			"mov eax, %0\n"                       // Load pointers into CPU regs and multiply
			"mov ebx, %1\n"
			"movups xmm0, [eax]\n"
			"movups xmm1, [ebx]\n"
			"mulps xmm1, xmm0\n"
			/* Sum components */
			"haddps xmm1, xmm1\n"
			"haddps xmm1, xmm1\n"

			//coloca 0-31 de xmm2 em dummy
			"mov		esi, %2\n"
			"movss		[esi], xmm1\n"
		:
		: "r" (pSrc1), "r" (pSrc2), "m" (min)
		);

	return dummy;


}
//Returns the squared length of pSrc1. This can be useful when doing sphere-sphere collision tests where the exact distance is not needed and a squared one can be used, for example.
//It is signifigantly faster than using the square root version SIMDx86Vector_Length().

float  VPCALL CSIMD_SSE3::vector3D_LengthSq(const CVec3D* pSrc1)
{

		float dummy;

		float *min = &dummy;
		asm
        (
			"mov eax, %0\n"                       // Load pointers into CPU regs and multiply
			"movups xmm0, [eax]\n"
			"andps xmm0, %1\n"  //remove w component /* xmm0 = x | y | z | 0.0f */
			"mulps xmm0, xmm0\n"    /* xmm0 = xmm0^2 */
			/* Sum components */
			"haddps xmm0, xmm0\n"
			"haddps xmm0, xmm0\n"

			/* Store */
		    //coloca 0-31 de xmm1 em dummy
			"mov			esi, %2\n"
			"movss		[esi], xmm0\n"


		:
		: "r" (pSrc1), "m" (_SIMDx86_float_SSE_NO_W_MASK),"m" (min)
		);
		return dummy;

}
//Returns the length (magnitude) of pSrc1. Normalized vectors (e.g. the output of SIMDx86Vector_Normalize() ) are of unit length, or 1.0)
float  VPCALL CSIMD_SSE3::vector3D_Length(const CVec3D* pSrc1)
{
		float dummy;

float *min = &dummy;
		asm
        (
			"mov eax, %0\n"                       // Load pointers into CPU regs and multiply
			"movups xmm0, [eax]\n"
			"andps xmm0, %1\n"  //remove w component /* xmm0 = x | y | z | 0.0f */
			"mulps xmm0, xmm0\n"    /* xmm0 = xmm0^2 */
			/* Sum components */
			"haddps xmm0, xmm0\n"
			"haddps xmm0, xmm0\n"

			#ifdef HIPREC
			/* Full square root */
			"sqrtss xmm0, xmm0\n"
			#else
			/* rcp( rsqrt(value) ) (This may be very inaccurate) */
			"rsqrtss xmm0, xmm0\n"
			"rcpss xmm0, xmm0\n"
			#endif
			/* Store */
		    //coloca 0-31 de xmm1 em dummy
			"mov			esi, %2\n"
			"movss		[esi], xmm0\n"

		:
		: "r" (pSrc1), "m" (_SIMDx86_float_SSE_NO_W_MASK),"m" (min)
		);
		return dummy;

}
void VPCALL CSIMD_SSE3::vector3D_Normalize(CVec3D* pVec)
{

		/* SSE3 Implementation */
		asm
        (
			"mov eax, %0\n"
			"mov ecx, %0\n"
			"movups xmm0, [eax]\n"  /* xmm0 = pVec */
			"andps xmm0, %1\n"  //remove w component /* xmm0 = x | y | z | 0.0f */
			"movaps xmm2, xmm0\n"			/* Save for the division by length */
		    "mulps xmm0, xmm0\n"    /* xmm0 = xmm0^2 */
			/* Sum components */
			"haddps xmm0, xmm0\n"
			"haddps xmm0, xmm0\n"

		/* High Precision SSE3 */
		#ifdef HIPREC

			/* Do full divide by the square root in parallel */
			"sqrtps xmm0, xmm0\n"				/* xmm0 = mag(pVec) */
			"divps xmm2, xmm0\n"				/* xmm2 = ? | z / mag | y / mag | x / mag */

		#else /* Low Precision SSE3 */

			/* Multiply by reciprocal square root */
			"rsqrtps xmm1, xmm0				/* xmm1 = rcp(pVec) | rcp(pVec) | rcp(pVec) | rcp(pVec) */
			"mulps xmm2, xmm1				/* xmm2 = norm(pVecv) */

		#endif

		/* Store */
		"movlpd	[ecx+ 0], xmm2\n"  //retira os tres bytes
		"shufps xmm2,xmm2, 0X02\n"   /*R_SHUFFLEPS( 2, 0, 0, 0 )*/  //shuffle colocar nos bits 0-31 o 3o byte
		"movss	[ecx+ 8], xmm2\n"
			:
		: "r" (pVec), "m" (_SIMDx86_float_SSE_NO_W_MASK)
		);


}
void VPCALL CSIMD_SSE3::vector3D_NormalizeOf(CVec3D* pOut, const CVec3D* pVec)
{

		/* SSE3 Implementation */
		asm
        (
			"mov eax, %0\n"
			"mov ecx, %1\n"
			"movups xmm0, [eax]\n"  /* xmm0 = pVec */
			"andps xmm0, %2\n"  //remove w component /* xmm0 = x | y | z | 0.0f */
			"movaps xmm2, xmm0\n"			/* Save for the division by length */
		    "mulps xmm0, xmm0\n"    /* xmm0 = xmm0^2 */
			/* Sum components */
			"haddps xmm0, xmm0\n"
			"haddps xmm0, xmm0\n"

		/* High Precision SSE3 */
		#ifdef HIPREC

			/* Do full divide by the square root in parallel */
			"sqrtps xmm0, xmm0\n"				/* xmm0 = mag(pVec) */
			"divps xmm2, xmm0\n"				/* xmm2 = ? | z / mag | y / mag | x / mag */

		#else /* Low Precision SSE3 */

			/* Multiply by reciprocal square root */
			"rsqrtps xmm1, xmm0\n"				/* xmm1 = rcp(pVec) | rcp(pVec) | rcp(pVec) | rcp(pVec) */
			"mulps xmm2, xmm1\n"				/* xmm2 = norm(pVecv) */

		#endif

		/* Store */
		"movlpd	[ecx+ 0], xmm2\n"  //retira os tres bytes
		"shufps xmm2,xmm2, 0x02\n"   /*R_SHUFFLEPS( 2, 0, 0, 0 )\n" //shuffle colocar nos bits 0-31 o 3o byte*/
		"movss	[ecx+ 8], xmm2\n"
		:
		: "r" (pVec), "r" (pOut), "m" (_SIMDx86_float_SSE_NO_W_MASK)
		);

}
float  VPCALL CSIMD_SSE3::vector3D_Distance(const CVec3D* pVec1, const CVec3D* pVec2)
{
	/* TODO: Optimize me completely */
	CVec3D diff;
	vector3D_DiffOf(&diff, pVec1, pVec2);
	return vector3D_Length(&diff);
}


#endif //masm_ATT
#endif /* _WIN32 */

} //end MATH
}//END SGF
