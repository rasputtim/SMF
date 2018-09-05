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
#include "math/SMF_Math.h"
#include "math/SMF_Matriz.h"
#include "math/SMF_Vector.h"
#include "math/SMF_SimdGeneric.h"
#include "math/SMF_SimdMMX.h"
#include "math/SMF_SimdSSE.h"
#include "geometry/SMF_Plane.h"
#include "geometry/SMF_DrawVert.h"
#include "math/SMF_JointTransform.h"

static unsigned long testeVar=0xFFFFFFF;

namespace SMF {
namespace MATH{
#define STRINGIZE(X) #X
ALIGN4_INIT1( unsigned long SIMD_DW_mat2quatShuffle0, (3<<0)|(2<<8)|(1<<16)|(0<<24) );
ALIGN4_INIT1( unsigned long SIMD_DW_mat2quatShuffle1, (0<<0)|(1<<8)|(2<<16)|(3<<24) );
ALIGN4_INIT1( unsigned long SIMD_DW_mat2quatShuffle2, (1<<0)|(0<<8)|(3<<16)|(2<<24) );
ALIGN4_INIT1( unsigned long SIMD_DW_mat2quatShuffle3, (2<<0)|(3<<8)|(0<<16)|(1<<24) );
ALIGN4_INIT4( unsigned int  SIMDx86_float_POSPOSPOSNEG, 0x00000000, 0x00000000, 0x00000000, 0x80000000 );
ALIGN4_INIT4( unsigned int  _SIMDx86_float_SSE_NO_W_MASK, 0xFFFFFFFF,  0xFFFFFFFF, 0xFFFFFFFF, 0x00000000 );
ALIGN4_INIT4( unsigned int  _SIMDx86_float_SSE_NO_XYZ_MASK, 0x00000000, 0x00000000, 0x00000000, 0xFFFFFFFF );
//(1 << 31 ) == 0x80000000  ==
//(1 << 23 ) == 0x00800000  ==
//(1 << 15 ) == 0x00008000  ==
ALIGN4_INIT4( unsigned int  _SIMDx86_float_NEGPOSPOSPOS, 0x80000000, 0x00000000, 0x00000000, 0x00000000 );
ALIGN4_INIT4( unsigned int POSNEGPOSNEG, 0x00000000, 0x80000000, 0x00000000, 0x80000000 );
ALIGN4_INIT4( unsigned int POSPOSNEGNEG, 0x00000000, 0x00000000, 0x80000000, 0x80000000 );
ALIGN4_INIT4( unsigned int NEGPOSPOSNEG, 0x80000000, 0x00000000, 0x00000000, 0x80000000 );

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
//TODO IMPLEMENT GNUC

/*
============
CSIMD_SSE::getName
============
*/
const char * CSIMD_SSE::getName() const {
	return "MMX & SSE";
}
#if defined(MASM_INTEL)
#include <xmmintrin.h>


// transpose a 4x4 matrix loaded into 4 xmm registers (reg4 is temporary)
#define TRANSPOSE_4x4( reg0, reg1, reg2, reg3, reg4 )											\
	__asm__ ("	movaps		reg4, reg2\n");								/* reg4 =  8,  9, 10, 11 */		\
	__asm__ ("	unpcklps	reg2, reg3\n");									/* reg2 =  8, 12,  9, 13 */		\
	__asm__ ("	unpckhps	reg4, reg3\n");									/* reg4 = 10, 14, 11, 15 */		\
	__asm__ ("	movaps		reg3, reg0\n");									/* reg3 =  0,  1,  2,  3 */		\
	__asm__ ("	unpcklps	reg0, reg1\n");									/* reg0 =  0,  4,  1,  5 */		\
	__asm__ ("	unpckhps	reg3, reg1\n");									/* reg3 =  2,  6,  3,  7 */		\
	__asm__ ("	movaps		reg1, reg0\n");									/* reg1 =  0,  4,  1,  5 */		\
	__asm__ ("	shufps		reg0, reg2, 0x44\n");		/* reg0 =  0,  4,  8, 12 */		\
	__asm__ ("	shufps		reg1, reg2, R_SHUFFLEPS( 2, 3, 2, 3 )\n");		/* reg1 =  1,  5,  9, 13 */		\
	__asm__ ("	movaps		reg2, reg3\n");									/* reg2 =  2,  6,  3,  7 */		\
	__asm__ ("	shufps		reg2, reg4, 0x44\n");		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/		/* reg2 =  2,  6, 10, 14 */		\
	__asm__ ("	shufps		reg3, reg4, R_SHUFFLEPS( 2, 3, 2, 3 )\n");		/* reg3 =  3,  7, 11, 15 */



// transpose a 4x4 matrix from memory into 4 xmm registers (reg4 is temporary)
#define TRANPOSE_4x4_FROM_MEMORY( address, reg0, reg1, reg2, reg3, reg4 )						\
	__asm__ ("	movlps		reg1, [address+ 0]\n");						/* reg1 =  0,  1,  X,  X */		\
	__asm__ ("	movlps		reg3, [address+ 8]\n");						/* reg3 =  2,  3,  X,  X */		\
	__asm__ ("	movhps		reg1, [address+16]\n");						/* reg1 =  0,  1,  4,  5 */		\
	__asm__ ("	movhps		reg3, [address+24]\n");						/* reg3 =  2,  3,  6,  7 */		\
	__asm__ ("	movlps		reg2, [address+32]\n");						/* reg2 =  8,  9,  X,  X */		\
	__asm__ ("	movlps		reg4, [address+40]\n");						/* reg4 = 10, 11,  X,  X */		\
	__asm__ ("	movhps		reg2, [address+48]\n");						/* reg2 =  8,  9, 12, 13 */		\
	__asm__ ("	movhps		reg4, [address+56]\n");						/* reg4 = 10, 11, 14, 15 */		\
	__asm__ ("	movaps		reg0, reg1\n");								/* reg0 =  0,  1,  4,  5 */		\
	__asm__ ("	shufps		reg0, reg2, R_SHUFFLEPS( 0, 2, 0, 2 )\n");	/* reg0 =  0,  4,  8, 12 */		\
	__asm__ ("	shufps		reg1, reg2, R_SHUFFLEPS( 1, 3, 1, 3 )\n");	/* reg1 =  1,  5,  9, 13 */		\
	__asm__ ("	movaps		reg2, reg3\n");								/* reg2 =  2,  3,  6,  7 */		\
	__asm__ ("	shufps		reg2, reg4, R_SHUFFLEPS( 0, 2, 0, 2 )\n");	/* reg2 =  2,  6, 10, 14 */		\
	__asm__ ("	shufps		reg3, reg4, R_SHUFFLEPS( 1, 3, 1, 3 )\n");	/* reg3 =  3,  7, 11, 15 */

// transpose a 4x4 matrix to memory from 4 xmm registers (reg4 is temporary)
#define TRANPOSE_4x4_TO_MEMORY( address, reg0, reg1, reg2, reg3, reg4 )							\
	__asm__ ("	movaps		reg4, reg0\n");								/* reg4 =  0,  4,  8, 12 */		\
	__asm__ ("	unpcklps	reg0, reg1\n");								/* reg0 =  0,  1,  4,  5 */		\
	__asm__ ("	unpckhps	reg4, reg1\n");								/* reg4 =  8,  9, 12, 13 */		\
	__asm__ ("	movaps		reg1, reg2\n");								/* reg1 =  2,  6, 10, 14 */		\
	__asm__ ("	unpcklps	reg2, reg3\n");								/* reg2 =  2,  3,  6,  7 */		\
	__asm__ ("	unpckhps	reg1, reg3\n");								/* reg1 = 10, 11, 14, 15 */		\
	__asm__ ("	movlps		[address+ 0], reg0\n");						/* mem0 =  0,  1,  X,  X */		\
	__asm__ ("	movlps		[address+ 8], reg2\n");						/* mem0 =  0,  1,  2,  3 */		\
	__asm__ ("	movhps		[address+16], reg0\n");						/* mem1 =  4,  5,  X,  X */		\
	__asm__ ("	movhps		[address+24], reg2\n");						/* mem1 =  4,  5,  6,  7 */		\
	__asm__ ("	movlps		[address+32], reg4\n");						/* mem2 =  8,  9,  X,  X */		\
	__asm__ ("	movlps		[address+40], reg1\n");						/* mem2 =  8,  9, 10, 11 */		\
	__asm__ ("	movhps		[address+48], reg4\n");						/* mem3 = 12, 13,  X,  X */		\
	__asm__ ("	movhps		[address+56], reg1\n");						/* mem3 = 12, 13, 14, 15 */

// transpose a 4x3 matrix loaded into 3 xmm registers (reg3 is temporary)
#define TRANSPOSE_4x3( reg0, reg1, reg2, reg3 )													\
	__asm__ ("	movaps		reg3, reg2\n");									/* reg3 =  8,  9, 10, 11 */		\
	__asm__ ("	shufps		reg3, reg1, 0x4E\n" /*R_SHUFFLEPS( 2, 3, 0, 1 )=01001110*/\n");		/* reg3 = 10, 11,  4,  5 */		\
	__asm__ ("	shufps		reg2, reg0, R_SHUFFLEPS( 0, 1, 2, 3 )\n");		/* reg2 =  8,  9,  2,  3 */		\
	__asm__ ("	shufps		reg1, reg0, 0x4E\n" /*R_SHUFFLEPS( 2, 3, 0, 1 )=01001110*/\n");		/* reg1 =  6,  7,  0,  1 */		\
	__asm__ ("	movaps		reg0, reg1\n");									/* reg0 =  6,  7,  0,  1 */		\
	__asm__ ("	shufps		reg0, reg2, R_SHUFFLEPS( 2, 0, 3, 1 )\n");		/* reg0 =  0,  6,  3,  9 */		\
	__asm__ ("	shufps		reg1, reg3, R_SHUFFLEPS( 3, 1, 2, 0 )\n");		/* reg1 =  1,  7,  4, 10 */		\
	__asm__ ("	shufps		reg2, reg3, R_SHUFFLEPS( 2, 0, 3, 1 )\n");		/* reg2 =  2,  8,  5, 11 */

// transpose a 4x3 matrix from memory into 3 xmm registers (reg3 is temporary)
#define TRANSPOSE_4x3_FROM_MEMORY( address, reg0, reg1, reg2, reg3 )							\
	__asm__ ("	movlps		reg1, [address+ 0]\n");							/* reg1 =  0,  1,  X,  X */		\
	__asm__ ("	movlps		reg2, [address+ 8]\n");							/* reg2 =  2,  3,  X,  X */		\
	__asm__ ("	movlps		reg3, [address+16]\n");							/* reg3 =  4,  5,  X,  X */		\
	__asm__ ("	movhps		reg1, [address+24]\n");							/* reg1 =  0,  1,  6,  7 */		\
	__asm__ ("	movhps		reg2, [address+32]\n");							/* reg2 =  2,  3,  8,  9 */		\
	__asm__ ("	movhps		reg3, [address+40]\n");							/* reg3 =  4,  5, 10, 11 */		\
	__asm__ ("	movaps		reg0, reg1\n");									/* reg0 =  0,  1,  6,  7 */		\
	__asm__ ("	shufps		reg0, reg2, 0xd8 \n" /*R_SHUFFLEPS( 0, 2, 1, 3 )=11011000*/\n");		/* reg0 =  0,  6,  3,  9 */		\
	__asm__ ("	shufps		reg1, reg3, 0x8d\n"  /*R_SHUFFLEPS( 1, 3, 0, 2 )*/\n");		/* reg1 =  1,  7,  4, 10 */		\
	__asm__ ("	shufps		reg2, reg3, 0xd8 \n" /*R_SHUFFLEPS( 0, 2, 1, 3 )=11011000*/\n");		/* reg2 =  2,  8,  5, 11 */

// transpose a 4x3 matrix to memory from 3 xmm registers (reg3 is temporary)
#define TRANSPOSE_4x3_TO_MEMORY( address, reg0, reg1, reg2, reg3 )								\
	__asm__ ("	movhlps		reg3, reg0\n");									/* reg3 =  3,  9,  X,  X */		\
	__asm__ ("	unpcklps	reg0, reg1\n");									/* reg0 =  0,  1,  6,  7 */		\
	__asm__ ("	unpckhps	reg1, reg2\n");									/* reg1 =  4,  5, 10, 11 */		\
	__asm__ ("	unpcklps	reg2, reg3\n");									/* reg2 =  2,  3,  8,  9 */		\
	__asm__ ("	movlps		[address+ 0], reg0\n");							/* mem0 =  0,  1,  X,  X */		\
	__asm__ ("	movlps		[address+ 8], reg2\n");							/* mem0 =  0,  1,  2,  3 */		\
	__asm__ ("	movlps		[address+16], reg1\n");							/* mem1 =  4,  5,  X,  X */		\
	__asm__ ("	movhps		[address+24], reg0\n");							/* mem1 =  4,  5,  6,  7 */		\
	__asm__ ("	movhps		[address+32], reg2\n");							/* mem2 =  8,  9,  X,  X */		\
	__asm__ ("	movhps		[address+40], reg1\n");							/* mem2 =  8,  9, 10, 11 */


// with alignment
#define KFLOATINITS(   SRC0, COUNT, PRE, POST )				KFLOATINITDSS( SRC0,SRC0,SRC0,COUNT,PRE,POST )
#define KFLOATINITD(   DST, COUNT, PRE, POST )				KFLOATINITDSS( DST,DST,DST,COUNT,PRE,POST )
#define KFLOATINITDS(  DST, SRC0, COUNT, PRE, POST )		KFLOATINITDSS( DST,SRC0,SRC0,COUNT,PRE,POST )
//dica label: http://stackoverflow.com/questions/3898435/labels-in-gcc-inline-assembly
#define KFLOATINITDSS( DST, SRC0, SRC1, COUNT, PRE, POST )\
	"mov		ecx,"STRINGIZE(DST)"\n"								\
	"shr		ecx,2\n"								\
	"mov		ebx,"STRINGIZE(COUNT)"\n"							\
	"neg		ecx\n"									\
	"mov		edx,"STRINGIZE(SRC0)"\n"							\
	"and		ecx,3\n"								\
	"mov		esi,"STRINGIZE(SRC1)"\n"							\
	"sub		ebx,ecx\n"								\
	"jge		noUnderFlow%=\n"							\
	"xor		ebx,ebx\n"								\
	"mov		ecx,"STRINGIZE(COUNT)"\n"							\
	"noUnderFlow%=:\n"								\
	"mov		"STRINGIZE(PRE)",ecx\n"								\
	"mov		eax,ebx\n"								\
	"mov		edi,"STRINGIZE(DST)"\n"					\
	"and		eax,8-1\n"								\
	"mov		"STRINGIZE(POST)",eax\n"							\
	"and		ebx,0xfffffff8\n"						\
	"jle		done%=\n"								\
	"shl		ebx,2\n"								\
	"lea		ecx,[(ecx*4)+ebx]\n"						\
	"neg		ebx\n"									\
	"add		edx,ecx\n"								\
	"add		esi,ecx\n"								\
	"add		edi,ecx\n"								\
	"mov		eax,edx\n"								\
	"or		    eax,esi\n"




// without alignment (pre==0)
#define KFLOATINITS_NA(   SRC0, COUNT, PRE, POST )				KFLOATINITDSS_NA( SRC0,SRC0,SRC0,COUNT,PRE,POST )
#define KFLOATINITD_NA(   DST, COUNT, PRE, POST )				KFLOATINITDSS_NA( DST,DST,DST,COUNT,PRE,POST )
#define KFLOATINITDS_NA(  DST, SRC0, COUNT, PRE, POST )			KFLOATINITDSS_NA( DST,SRC0,SRC0,COUNT,PRE,POST )
#define KFLOATINITDSS_NA( DST, SRC0, SRC1, COUNT, PRE, POST )\
        "mov		eax,"STRINGIZE(COUNT)"\n"							\
		"mov		"STRINGIZE(PRE)",0\n"								\
		"and		eax,8-1\n"							\
		"mov		ebx,"STRINGIZE(COUNT)"\n"							\
		"mov		POST,eax\n"							\
		"and		ebx,0xfffffff8\n"						\
		"je		done%=\n"								\
		"shl		ebx,2\n"								\
		"mov		edx,"STRINGIZE(SRC0)"\n"							\
		"mov		esi,"STRINGIZE(SRC1)"\n"							\
		"mov		edi,"STRINGIZE(DST)"\n"								\
		"add		edx,ebx\n"								\
		"add		esi,ebx\n"								\
		"add		edi,ebx\n"								\
		"mov		eax,edx\n"								\
		"or		eax,esi\n"							\
		"or		eax,edi\n"								\
		"neg		ebx\n"									\


/*
	when OPER is called:
	edx = s0
	esi	= s1
	edi	= d
	ebx	= index*4

	xmm0 & xmm1	must not be trashed
*/
#define KMOVDS1( DST, SRC0 )							\
	"movss	xmm2,"STRINGIZE(SRC0)"\n"							\
	"movss	"STRINGIZE(DST)",xmm2\n"
#define KMOVDS4( DST, SRC0 )							\
	"movups	xmm2,"STRINGIZE(SRC0)"\n"							\
	"movups	"STRINGIZE(DST)",xmm2\n"
#define KMINDS1( DST, SRC0 )							\
	"movss	xmm2,"STRINGIZE(SRC0)"\n"							\
	"minss	"STRINGIZE(DST)",xmm2\n"
#define KMAXDS1( DST, SRC0 )							\
	"movss	xmm2,"STRINGIZE(SRC0)"\n"							\
	"maxss	"STRINGIZE(DST)",xmm2\n"
// general ALU operation
#define KALUDSS1( OP, DST, SRC0, SRC1 )					\
	"movss	xmm2,"STRINGIZE(SRC0)"\n"							\
	STRINGIZE(OP##ss) "  xmm2,"STRINGIZE(SRC1)"\n"						\
	"movss	"STRINGIZE(DST)",xmm2\n"
#define KALUDSS4( OP, DST, SRC0, SRC1 )					\
	"movups	xmm2,"STRINGIZE(SRC0)"\n"							\
	"movups	xmm3,"STRINGIZE(SRC1)"\n"							\
	STRINGIZE(OP##ss) "  xmm2,xmm3\n"							\
	"movups	"STRINGIZE(DST)",xmm2\n"

#define KADDDSS1( DST, SRC0, SRC1 )		KALUDSS1( add, DST,SRC0,SRC1 )
#define KADDDSS4( DST, SRC0, SRC1 )		KALUDSS4( add, DST,SRC0,SRC1 )
#define KSUBDSS1( DST, SRC0, SRC1 )		KALUDSS1( sub, DST,SRC0,SRC1 )
#define KSUBDSS4( DST, SRC0, SRC1 )		KALUDSS4( sub, DST,SRC0,SRC1 )
#define KMULDSS1( DST, SRC0, SRC1 )		KALUDSS1( mul, DST,SRC0,SRC1 )
#define KMULDSS4( DST, SRC0, SRC1 )		KALUDSS4( mul, DST,SRC0,SRC1 )

#define KDIVDSS1( DST, SRC0, SRC1 )						\
	"movss	xmm2,"STRINGIZE(SRC1)"\n"							\
	"rcpss	xmm3,xmm2\n"							\
	"mulss	xmm2,xmm3\n"							\
	"mulss	xmm2,xmm3\n"							\
	"addss	xmm3,xmm3\n"							\
	"subss	xmm3,xmm2\n"							\
	"mulss	xmm3,"STRINGIZE(SRC0)"\n"							\
	"movss	"STRINGIZE(DST)",xmm3\n"
#define KDIVDSS4( DST, SRC0, SRC1 )						\
	"movups	xmm2,"STRINGIZE(SRC1)"\n"							\
	"rcpps	xmm3,xmm2\n"							\
	"mulps	xmm2,xmm3\n"							\
	"mulps	xmm2,xmm3\n"							\
	"addps	xmm3,xmm3\n"							\
	"subps	xmm3,xmm2\n"							\
	"movups	xmm2,"STRINGIZE(SRC0)"\n"							\
	"mulps	xmm3,xmm2\n"							\
	"movups	"STRINGIZE(DST)",xmm3\n"
#define	KF2IDS1( SRC0 )									\
	"movss		xmm2,"STRINGIZE(SRC0)"\n"						\
	"cvttps2pi	mm2,xmm2\n"						\
	"movd		[edi+ebx],mm2\n"
#define	KF2IDS4( SRC0 )									\
	"movups		xmm2,"STRINGIZE(SRC0)"\n"						\
	"cvttps2pi	mm2,xmm2\n"						\
	"movq		[edi+ebx+0],mm2\n"					\
	"shufps		xmm2,xmm2,SHUFFLEPS(1,0,3,2)\n"	\
	"cvttps2pi	mm2,xmm2\n");						\
	"movq		[edi+ebx+8],mm2\n"
#define	KISQRTDS1( DST,SRC0 )							\
	"movss	xmm2,"STRINGIZE(SRC0)"\n"							\
	"rsqrtss	xmm3,xmm2\n"							\
	"mulss	xmm2,xmm3\n"							\
	"mulss	xmm2,xmm3\n"							\
	"subss	xmm2,xmm1\n"							\
	"mulss	xmm3,xmm0\n"							\
	"mulss	xmm3,xmm2\n"							\
	"movss	"STRINGIZE(DST)",xmm3\n"
#define	KISQRTDS4( DST,SRC0 )							\
	"movups	xmm2,"STRINGIZE(SRC0)"\n"							\
	"rsqrtps	xmm3,xmm2\n"							\
	"mulps	xmm2,xmm3\n"							\
	"mulps	xmm2,xmm3\n"							\
	"subps	xmm2,xmm1\n"							\
	"mulps	xmm3,xmm0\n"							\
	"mulps	xmm3,xmm2\n"							\
	"movups	"STRINGIZE(DST)",xmm3\n"

// this is used in vector4 implementation to shift constant V4
#define KANDREGDSV( DST, SRC0, VALUE )					\
	"mov		"STRINGIZE(DST)","STRINGIZE(SRC0)"\n"							\
	"and		"STRINGIZE(DST)","STRINGIZE(VALUE)"\n"

// this is used in vector4 code to operate with float arrays as sources
#define KEXPANDFLOAT( DST, SRC )						\
	"movss	"STRINGIZE(DST)","STRINGIZE(SRC)"\n"								\
	"shufps  "STRINGIZE(DST)","STRINGIZE(DST)",0\n"

#define	KADDDS1( DST,SRC )		KADDDSS1( DST,DST,SRC )
#define	KADDDS4( DST,SRC )		KADDDSS4( DST,DST,SRC )
#define	KSUBDS1( DST,SRC )		KSUBDSS1( DST,DST,SRC )
#define	KSUBDS4( DST,SRC )		KSUBDSS4( DST,DST,SRC )
#define	KMULDS1( DST,SRC )		KMULDSS1( DST,DST,SRC )
#define	KMULDS4( DST,SRC )		KMULDSS4( DST,DST,SRC )
#define	KDIVDS1( DST,SRC )		KDIVDSS1( DST,DST,SRC )
#define	KDIVDS4( DST,SRC )		KDIVDSS4( DST,DST,SRC )

// handles pre & post leftovers
#define	KFLOATOPER( OPER, OPER4, COUNT, PRE,POST )				\
	"mov		ecx,"STRINGIZE(PRE)"\n"							\
	"mov		ebx,"STRINGIZE(COUNT)"\n"						\
	"cmp		ebx,ecx\n"							\
	"cmovl	ecx,"STRINGIZE(COUNT)"\n"						\
	"test	ecx,ecx\n"							\
	"je		preDone%=\n"							\
	"xor		ebx,ebx\n"							\
	"lpPre%=:\n"										\
    OPER									\
	"add		ebx,4\n"							\
	"dec		ecx\n"								\
	"jg		lpPre%=\n"							\
	"preDone%=:\n"									\
	"mov		ecx,"STRINGIZE(POST)"\n"						\
	"mov		ebx,"STRINGIZE(COUNT)"\n"						\
	"sub		ebx,ecx	\n"						\
	"shl		ebx,2\n"							\
	"cmp		ecx,4\n"							\
	"jl		post4Done%=\n"						\
    OPER4									\
	"sub		ecx,4\n"							\
	"add		ebx,4*4\n"							\
	"post4Done%=:\n"									\
	"test	ecx,ecx\n"							\
	"je		postDone%=\n"						\
	"lpPost%=:\n"									\
    OPER									\
	"add		ebx,4\n"							\
	"dec		ecx\n"								\
	"jg		lpPost%=\n"							\
	"postDone%=:\n"

// operate on a constant and a float array
//#define KFLOAT_CA( ALUOP, DST_p, SRC_p, CONSTANT, COUNT )	\
//: "m"(constant),"m"(count),"r"(dstP),"r"(srcP),"m"(pre),"m"(post)
//KFLOATINITDS( dstP, srcP, count, pre, post )
//
#define KFLOAT_CA( ALUOP )	\
    "movss	xmm0,%0\n"					\
	"shufps	xmm0,xmm0,0\n"						\
	KFLOATINITDS( %2, %3, %1, %4, %5 )	\
	"and	eax,15\n"							\
	"jne	lpNA%=\n"							\
	"jmp	lpA%=\n"							\
	/*"align	16\n"	*/						\
	"lpA%=:\n"									\
	"prefetchnta	[edx+ebx+64]\n"				\
	"movaps	xmm1,xmm0\n"						\
	"movaps	xmm2,xmm0\n"						\
	STRINGIZE(ALUOP##ps) "	xmm1,[edx+ebx]\n"				\
	STRINGIZE(ALUOP##ps) "	xmm2,[edx+ebx+16]\n"			\
	"movaps	[edi+ebx],xmm1\n"					\
	"movaps	[edi+ebx+16],xmm2\n"				\
	"add		ebx,16*2\n"						\
	"jl		lpA%=\n"								\
	"jmp		done%=\n"							\
	/*"align	16\n"*/							\
	"lpNA%=:\n"								\
    "prefetchnta	[edx+ebx+64]\n"				\
	"movaps	xmm1,xmm0\n"						\
	"movaps	xmm2,xmm0\n"						\
	"movups	xmm3,[edx+ebx]\n"					\
	"movups	xmm4,[edx+ebx+16]\n"				\
	STRINGIZE(ALUOP##ps) "	xmm1,xmm3\n"					\
	STRINGIZE(ALUOP##ps) "	xmm2,xmm4\n"					\
	"movaps	[edi+ebx],xmm1\n"					\
	"movaps	[edi+ebx+16],xmm2\n"				\
	"add	ebx,16*2\n"						    \
	"jl		lpNA%=\n"							    \
	"done%=:\n"									\
	"mov		edx,%3\n"						\
    "mov		edi,%2\n"						\
	              KFLOATOPER( KALUDSS1( ALUOP, [edi+ebx],xmm0,[edx+ebx] ),	\
                    KALUDSS4( ALUOP, [edi+ebx],xmm0,[edx+ebx] ), %1 ,%4,%5)


// operate on two float arrays
//#define KFLOAT_AA( ALUOP, DST, SRC0, SRC1, COUNT )		\

#define KFLOAT_AA( ALUOP )		\
	KFLOATINITDSS( %1, %2, %3, %0, %4, %5 )	\
	"and		eax,15\n"							\
	"jne		lpNA%=\n"							\
	"jmp		lpA%=\n"								\
	/*"align	16\n"	*/							\
	"lpA%=:\n"										\
	"movaps	xmm1,[edx+ebx]\n"					\
	"movaps	xmm2,[edx+ebx+16]\n"				\
	STRINGIZE(ALUOP##ps) "	xmm1,[esi+ebx]\n"				\
	STRINGIZE(ALUOP##ps) "	xmm2,[esi+ebx+16]\n"			\
	"prefetchnta	[edx+ebx+64]\n"				\
	"prefetchnta	[esi+ebx+64]\n"				\
	"movaps	[edi+ebx],xmm1\n"					\
	"movaps	[edi+ebx+16],xmm2\n"				\
	"add		ebx,16*2\n"						\
	"jl		lpA%=	\n"							\
	"jmp		done%=\n"							\
	/*"align	16\n"	*/							\
	"lpNA%=:\n"										\
	"movups	xmm1,[edx+ebx]\n"					\
	"movups	xmm2,[edx+ebx+16]\n"				\
	"movups	xmm3,[esi+ebx]\n"					\
	"movups	xmm4,[esi+ebx+16]\n"				\
	"prefetchnta	[edx+ebx+64]\n"				\
	"prefetchnta	[esi+ebx+64]\n"				\
	STRINGIZE(ALUOP##ps) "	xmm1,xmm3\n"				\
	STRINGIZE(ALUOP##ps) "	xmm2,xmm4\n"					\
	"movaps	[edi+ebx],xmm1\n"					\
	"movaps	[edi+ebx+16],xmm2\n"				\
	"add	ebx,16*2\n"						\
	"jl		lpNA%=\n"							\
	"done%=:\n"										\
	"mov		edx,%2\n"						\
	"mov		esi,%3\n"						\
	"mov		edi,%1\n"					\
	KFLOATOPER( KALUDSS1( ALUOP, [edi+ebx],[edx+ebx],[esi+ebx] ),		\
				KALUDSS4( ALUOP, [edi+ebx],[edx+ebx],[esi+ebx] ), %0,%4,%5 )


#if defined(WIN32) && defined(_MSC_VER)

#define ALIGN4_INIT1( X, INIT )				ALIGN16( static X[4] ) = { INIT, INIT, INIT, INIT }
#define ALIGN4_INIT4( X, I0, I1, I2, I3 )	ALIGN16( static X[4] ) = { I0, I1, I2, I3 }
#define ALIGN8_INIT1( X, INIT )				ALIGN16( static X[8] ) = { INIT, INIT, INIT, INIT, INIT, INIT, INIT, INIT }

ALIGN8_INIT1( unsigned short SIMD_W_zero, 0 );
ALIGN8_INIT1( unsigned short SIMD_W_maxShort, 1<<15 );

ALIGN4_INIT1( unsigned long SIMD_DW_mat2quatShuffle0, (3<<0)|(2<<8)|(1<<16)|(0<<24) );
ALIGN4_INIT1( unsigned long SIMD_DW_mat2quatShuffle1, (0<<0)|(1<<8)|(2<<16)|(3<<24) );
ALIGN4_INIT1( unsigned long SIMD_DW_mat2quatShuffle2, (1<<0)|(0<<8)|(3<<16)|(2<<24) );
ALIGN4_INIT1( unsigned long SIMD_DW_mat2quatShuffle3, (2<<0)|(3<<8)|(0<<16)|(1<<24) );

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
#endif // defined

/*
============
SSE_InvSqrt
============
*/
float VPCALL  CSIMD_SSE::invSqrt( float x ) {
    ALIGNTO16 float y;
	asm (
		"movss		xmm0, %0\n"
		"rsqrtss	xmm1, xmm0\n"
		"mulss		xmm0, xmm1\n"
		"mulss		xmm0, xmm1\n"
		"subss		xmm0, %2\n"
		"mulss		xmm1, %3\n"
		"mulss		xmm0, xmm1\n"
		"movss		%1, xmm0\n"
        :
	    :"m"(x),"m"(y),"m"(SIMD_SP_rsqrt_c0),"m"(SIMD_SP_rsqrt_c1)
	    :);

	return y;
}
/*
============
SSE_InvSqrt4
============
*/
void VPCALL  CSIMD_SSE::InvSqrt4( float x[4] ) {
	asm(
		"mov		edi, %0\n"
		"movaps		xmm0, [edi]\n"
		"rsqrtps	xmm1, xmm0\n"
		"mulps		xmm0, xmm1\n"
		"mulps		xmm0, xmm1\n"
		"subps		xmm0, %1\n"
		"mulps		xmm1, %2\n"
		"mulps		xmm0, xmm1\n"
		"movaps		[edi], xmm0\n"
	:
	:"r"(x),"m"(SIMD_SP_rsqrt_c0),"m"(SIMD_SP_rsqrt_c1)
	:
	);



}

/*
============
SSE_SinZeroHalfPI

  The angle must be between zero and half PI.  (0 a 90 graus)
  \param x ângulo em radianos que se deseja calcular o seno
============
*/
float VPCALL  CSIMD_SSE::sinZeroHalfPI( float x ) {
#if 1
	ALIGNTO16 float resultado;

	SMF_ASSERT( x >= 0.0f && x <= MATH::CMath::HALF_PI );
	asm(
		//"mov        eax, %7\n"
		"movss		xmm0, %7\n"
		"movss		xmm1, xmm0\n"
		"mulss		xmm1, xmm1\n"
		"movss		xmm2, %0\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %1\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %2\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %3\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %4\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %5\n"
		"mulss		xmm2, xmm0\n"
		"movss		%6, xmm2\n"
	:
    : "m"(SIMD_SP_sin_c0),"m"(SIMD_SP_sin_c1),"m"(SIMD_SP_sin_c2),"m"(SIMD_SP_sin_c3),"m"(SIMD_SP_sin_c4),"m"(SIMD_SP_one),"m"(resultado),"m"(x)
    :
     );

	return resultado;

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
============
*/
void VPCALL  CSIMD_SSE::sin4ZeroHalfPI( float entrada[4], float saida[4] ) {
	SMF_ASSERT( entrada[0] >= 0.0f && entrada[0] <= MATH::CMath::HALF_PI );
	SMF_ASSERT( entrada[1] >= 0.0f && entrada[1] <= MATH::CMath::HALF_PI );
	SMF_ASSERT( entrada[2] >= 0.0f && entrada[2] <= MATH::CMath::HALF_PI );
	SMF_ASSERT( entrada[3] >= 0.0f && entrada[3] <= MATH::CMath::HALF_PI );
   	asm(
		"mov		edi, %7\n"
		"mov		esi, %6\n"
		"movaps		xmm0, [edi]\n"
		"movaps		xmm1, xmm0\n"
		"mulps		xmm1, xmm1\n"
		"movaps		xmm2, %0\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %1\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %2\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %3\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %4\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %5\n"
		"mulps		xmm2, xmm0\n"
		"movaps		[esi], xmm2\n"
        :
		: "m"(SIMD_SP_sin_c0),"m"(SIMD_SP_sin_c1),"m"(SIMD_SP_sin_c2),"m"(SIMD_SP_sin_c3),"m"(SIMD_SP_sin_c4),"m"(SIMD_SP_one),"r"(saida),"r"(entrada)
    );
}

/*
============sin
SSE_Sin
\param a ângulo em radianos
============
*/
float  VPCALL CSIMD_SSE::sin( float a ) {
#if 1

	float t;

	asm (
		"movss		xmm1, %13\n"
		"movss		xmm2, xmm1\n"
		"movss		xmm3, xmm1\n"
		"mulss		xmm2, %0\n"
		"cvttss2si	ecx, xmm2\n"
		"cmpltss	xmm3, %1\n"
		"andps		xmm3, %2\n"
		"cvtsi2ss	xmm2, ecx\n"
		"subss		xmm2, xmm3\n"
		"mulss		xmm2, %3\n"
		"subss		xmm1, xmm2\n"

		"movss		xmm0, %4\n"			// xmm0 = PI
		"subss		xmm0, xmm1\n"					// xmm0 = PI - a
		"movss		xmm1, xmm0\n"					// xmm1 = PI - a
		"andps		xmm1, %5\n"	// xmm1 = signbit( PI - a )
		"movss		xmm2, xmm0\n"					// xmm2 = PI - a
		"xorps		xmm2, xmm1\n"					// xmm2 = fabs( PI - a )
		"cmpnltss	xmm2, %6\n"		// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? 0xFFFFFFFF : 0x00000000
		"movss		xmm3, %4\n"			// xmm3 = PI
		"xorps		xmm3, xmm1\n"					// xmm3 = PI ^ signbit( PI - a )
		"andps		xmm3, xmm2\n"					// xmm3 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? ( PI ^ signbit( PI - a ) ) : 0.0f
		"andps		xmm2, %5\n"	// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? SIMD_SP_signBitMask : 0.0f
		"xorps		xmm0, xmm2\n"
		"addps		xmm0, xmm3\n"

		"movss		xmm1, xmm0\n"
		"mulss		xmm1, xmm1\n"
		"movss		xmm2, %7\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %8\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %9\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %10\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %11\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %12\n"
		"mulss		xmm2, xmm0\n"
		/* Store */
		"movss		%14, xmm2\n"
	    :
        : "m" (SIMD_SP_oneOverTwoPI), "m"(SIMD_SP_zero),"m"(SIMD_SP_one),"m"(SIMD_SP_twoPI),
          "m"(SIMD_SP_PI),"m"(SIMD_SP_signBitMask), "m"(SIMD_SP_halfPI),
          "m"(SIMD_SP_sin_c0),"m"(SIMD_SP_sin_c1),"m"(SIMD_SP_sin_c2),
          "m"(SIMD_SP_sin_c3),"m"(SIMD_SP_sin_c4),"m"(SIMD_SP_one),"m"(a),"m"(t)
        );
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
\param entrada  vetor com angulos em radianos
\param saida vetor com angulos em radianos
============
*/
void  VPCALL CSIMD_SSE::sin4( float entrada[4], float saida[4] ) {

    asm (
		"mov		edi, %13\n"
		"mov		esi, %14\n"
		"movaps		xmm1, [edi]\n"
		"movaps		xmm2, xmm1\n"
		"mulps		xmm2, %0\n"
		"movhlps	xmm3, xmm2\n"
		"cvttss2si	ecx, xmm2\n"
		"cvtsi2ss	xmm2, ecx\n"
		"cvttss2si	edx, xmm3\n"
		"cvtsi2ss	xmm3, edx\n"
		"shufps		xmm2, xmm2, 0x01\n"
		"shufps		xmm3, xmm3, 0x01\n"
		"cvttss2si	ecx, xmm2\n"
		"cvtsi2ss	xmm2, ecx\n"
		"cvttss2si	edx, xmm3\n"
		"cvtsi2ss	xmm3, edx\n"
		"shufps		xmm2, xmm3,0x11 \n"   /*R_SHUFFLEPS( 1, 0, 1, 0 )*/
		"movaps		xmm3, xmm1\n"
		"cmpltps	xmm3, %1\n"
		"andps		xmm3, %2\n"
		"subps		xmm2, xmm3\n"
		"mulps		xmm2, %3\n"
		"subps		xmm1, xmm2\n"

		"movaps		xmm0, %4\n"			// xmm0 = PI
		"subps		xmm0, xmm1\n"					// xmm0 = PI - a
		"movaps		xmm1, xmm0\n"					// xmm1 = PI - a
		"andps		xmm1, %5\n"	// xmm1 = signbit( PI - a )
		"movaps		xmm2, xmm0\n"					// xmm2 = PI - a
		"xorps		xmm2, xmm1\n"					// xmm2 = fabs( PI - a )
		"cmpnltps	xmm2, %6\n"	// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? 0xFFFFFFFF : 0x00000000
		"movaps		xmm3, %4\n"			// xmm3 = PI
		"xorps		xmm3, xmm1\n"					// xmm3 = PI ^ signbit( PI - a )
		"andps		xmm3, xmm2\n"					// xmm3 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? ( PI ^ signbit( PI - a ) ) : 0.0f
		"andps		xmm2, %5\n"	// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? SIMD_SP_signBitMask : 0.0f
		"xorps		xmm0, xmm2\n"
		"addps		xmm0, xmm3\n"

		"movaps		xmm1, xmm0\n"
		"mulps		xmm1, xmm1\n"
		"movaps		xmm2, %7\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %8\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %9\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %10\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %11\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %12\n"
		"mulps		xmm2, xmm0\n"
		"movaps		[esi], xmm2\n"
        :
	    : "m" (SIMD_SP_oneOverTwoPI), "m"(SIMD_SP_zero),"m"(SIMD_SP_one),"m"(SIMD_SP_twoPI),
          "m"(SIMD_SP_PI),"m"(SIMD_SP_signBitMask), "m"(SIMD_SP_halfPI),
          "m"(SIMD_SP_sin_c0),"m"(SIMD_SP_sin_c1),"m"(SIMD_SP_sin_c2),
          "m"(SIMD_SP_sin_c3),"m"(SIMD_SP_sin_c4),"m"(SIMD_SP_one),"r"(entrada),"r"(saida)
        );

}

/*
============
SSE_CosZeroHalfPI

  The angle must be between zero and half PI.
  \param x ângulo em radianos
============
*/
float  VPCALL CSIMD_SSE::cosZeroHalfPI( float x ) {

#if 1

	float resultado;

	SMF_ASSERT( x >= 0.0f && x <= CMath::HALF_PI );

	asm (
		"movss		xmm0, %7\n"
		"mulss		xmm0, xmm0\n"
		"movss		xmm1, %0\n"
		"mulss		xmm1, xmm0\n"
		"addss		xmm1, %1\n"
		"mulss		xmm1, xmm0\n"
		"addss		xmm1, %2\n"
		"mulss		xmm1, xmm0\n"
		"addss		xmm1, %3\n"
		"mulss		xmm1, xmm0\n"
		"addss		xmm1, %4\n"
		"mulss		xmm1, xmm0\n"
		"addss		xmm1, %5\n"
		"movss		%6, xmm1\n"
	:
    : "m"(SIMD_SP_cos_c0),"m"(SIMD_SP_cos_c1),"m"(SIMD_SP_cos_c2),"m"(SIMD_SP_cos_c3),"m"(SIMD_SP_cos_c4),"m"(SIMD_SP_one),"m"(resultado),"m"(x)
    :
     );

	return resultado;

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
void  VPCALL CSIMD_SSE::cos4ZeroHalfPI( float entrada[4], float saida[4] ) {
	SMF_ASSERT( entrada[0] >= 0.0f && entrada[0] <= MATH::CMath::HALF_PI );
	SMF_ASSERT( entrada[1] >= 0.0f && entrada[1] <= MATH::CMath::HALF_PI );
	SMF_ASSERT( entrada[2] >= 0.0f && entrada[2] <= MATH::CMath::HALF_PI );
	SMF_ASSERT( entrada[3] >= 0.0f && entrada[3] <= MATH::CMath::HALF_PI );
   	asm(
		"mov			edi, %7\n"
		"mov			esi, %6\n"
		"movaps		xmm0, [edi]\n"
		"mulps		xmm0, xmm0\n"
		"movaps		xmm1, %0\n"
		"mulps		xmm1, xmm0\n"
		"addps		xmm1, %1\n"
		"mulps		xmm1, xmm0\n"
		"addps		xmm1, %2\n"
		"mulps		xmm1, xmm0\n"
		"addps		xmm1, %3\n"
		"mulps		xmm1, xmm0\n"
		"addps		xmm1, %4\n"
		"mulps		xmm1, xmm0\n"
		"addps		xmm1, %5\n"
		"movaps		[esi], xmm2\n"
        :
		: "m"(SIMD_SP_cos_c0),"m"(SIMD_SP_cos_c1),"m"(SIMD_SP_cos_c2),"m"(SIMD_SP_cos_c3),"m"(SIMD_SP_cos_c4),"m"(SIMD_SP_one),"r"(saida),"r"(entrada)
    );

}
/*
============
SSE_Cos
============
*/
float  VPCALL CSIMD_SSE::cos( float a ) {
#if 1

	float t;

	asm (
		"movss		xmm1, %13\n"
		"movss		xmm2, xmm1\n"
		"movss		xmm3, xmm1\n"
		"mulss		xmm2, %0\n"
		"cvttss2si	ecx, xmm2\n"
		"cmpltss	xmm3, %1\n"
		"andps		xmm3, %2\n"
		"cvtsi2ss	xmm2, ecx\n"
		"subss		xmm2, xmm3\n"
		"mulss		xmm2, %3\n"
		"subss		xmm1, xmm2\n"

		"movss		xmm0, %4\n"			// xmm0 = PI
		"subss		xmm0, xmm1\n"					// xmm0 = PI - a
		"movss		xmm1, xmm0\n"					// xmm1 = PI - a
		"andps		xmm1, %5\n"	// xmm1 = signbit( PI - a )
		"movss		xmm2, xmm0\n"					// xmm2 = PI - a
		"xorps		xmm2, xmm1\n"					// xmm2 = fabs( PI - a )
		"cmpnltss	xmm2, %6\n"		// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? 0xFFFFFFFF : 0x00000000
		"movss		xmm3, %4\n"			// xmm3 = PI
		"xorps		xmm3, xmm1\n"					// xmm3 = PI ^ signbit( PI - a )
		"andps		xmm3, xmm2\n"					// xmm3 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? ( PI ^ signbit( PI - a ) ) : 0.0f
		"andps		xmm2, %5\n"	// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? SIMD_SP_signBitMask : 0.0f
		"xorps		xmm0, xmm2\n"
		"addps		xmm0, xmm3\n"

		"mulss		xmm0, xmm0\n"
		"movss		xmm1, %7\n"
		"mulss		xmm1, xmm0\n"
		"addss		xmm1, %8\n"
		"mulss		xmm1, xmm0\n"
		"addss		xmm1, %9\n"
		"mulss		xmm1, xmm0\n"
		"addss		xmm1, %10\n"
		"mulss		xmm1, xmm0\n"
		"addss		xmm1, %11\n"
		"mulss		xmm1, xmm0\n"
		"addss		xmm1, %12\n"
		"xorps		xmm2, %5\n"
		"xorps		xmm1, xmm2\n"
		"movss		%14, xmm1\n"
        :
        : "m" (SIMD_SP_oneOverTwoPI), "m"(SIMD_SP_zero),"m"(SIMD_SP_one),"m"(SIMD_SP_twoPI),
          "m"(SIMD_SP_PI),"m"(SIMD_SP_signBitMask), "m"(SIMD_SP_halfPI),
          "m"(SIMD_SP_cos_c0),"m"(SIMD_SP_cos_c1),"m"(SIMD_SP_cos_c2),
          "m"(SIMD_SP_cos_c3),"m"(SIMD_SP_cos_c4),"m"(SIMD_SP_one),"m"(a),"m"(t)
        );

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
void  VPCALL CSIMD_SSE::cos4( float entrada[4], float saida[4] ) {
	asm (
		"mov		edi, %13\n"
		"mov		esi, %14\n"
		"movaps		xmm1, [edi]\n"
		"movaps		xmm2, xmm1\n"
		"mulps		xmm2, %0\n"
		"movhlps	xmm3, xmm2\n"
		"cvttss2si	ecx, xmm2\n"
		"cvtsi2ss	xmm2, ecx\n"
		"cvttss2si	edx, xmm3\n"
		"cvtsi2ss	xmm3, edx\n"
		"shufps		xmm2, xmm2, 0x01\n"     /*R_SHUFFLEPS( 1, 0, 0, 0 )*/
		"shufps		xmm3, xmm3, 0x01\n"     /*R_SHUFFLEPS( 1, 0, 0, 0 )*/
		"cvttss2si	ecx, xmm2\n"
		"cvtsi2ss	xmm2, ecx\n"
		"cvttss2si	edx, xmm3\n"
		"cvtsi2ss	xmm3, edx\n"
		"shufps		xmm2, xmm3, 0x11\n"     /*R_SHUFFLEPS( 1, 0, 1, 0 )*/
		"movaps		xmm3, xmm1\n"
		"cmpltps	xmm3, %1\n"
		"andps		xmm3, %2\n"
		"subps		xmm2, xmm3\n"
		"mulps		xmm2, %3\n"
		"subps		xmm1, xmm2\n"

		"movaps		xmm0, %4\n"			// xmm0 = PI
		"subps		xmm0, xmm1\n"					// xmm0 = PI - a
		"movaps		xmm1, xmm0\n"					// xmm1 = PI - a
		"andps		xmm1, %5\n"	// xmm1 = signbit( PI - a )
		"movaps		xmm2, xmm0\n"					// xmm2 = PI - a
		"xorps		xmm2, xmm1\n"					// xmm2 = fabs( PI - a )
		"cmpnltps	xmm2, %6\n"		// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? 0xFFFFFFFF : 0x00000000
		"movaps		xmm3, %4\n"			// xmm3 = PI
		"xorps		xmm3, xmm1\n"					// xmm3 = PI ^ signbit( PI - a )
		"andps		xmm3, xmm2\n"					// xmm3 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? ( PI ^ signbit( PI - a ) ) : 0.0f
		"andps		xmm2, %5\n"	// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? SIMD_SP_signBitMask : 0.0f
		"xorps		xmm0, xmm2\n"
		"addps		xmm0, xmm3\n"

		"mulps		xmm0, xmm0\n"
		"movaps		xmm1, %7\n"
		"mulps		xmm1, xmm0\n"
		"addps		xmm1, %8\n"
		"mulps		xmm1, xmm0\n"
		"addps		xmm1, %9\n"
		"mulps		xmm1, xmm0\n"
		"addps		xmm1, %10\n"
		"mulps		xmm1, xmm0\n"
		"addps		xmm1, %11\n"
		"mulps		xmm1, xmm0\n"
		"addps		xmm1, %12\n"
		"xorps		xmm2, %5\n"
		"xorps		xmm1, xmm2\n"
		"movaps		[esi], xmm1\n"
        :
	    : "m" (SIMD_SP_oneOverTwoPI), "m"(SIMD_SP_zero),"m"(SIMD_SP_one),"m"(SIMD_SP_twoPI),
          "m"(SIMD_SP_PI),"m"(SIMD_SP_signBitMask), "m"(SIMD_SP_halfPI),
          "m"(SIMD_SP_cos_c0),"m"(SIMD_SP_cos_c1),"m"(SIMD_SP_cos_c2),
          "m"(SIMD_SP_cos_c3),"m"(SIMD_SP_cos_c4),"m"(SIMD_SP_one),"r"(entrada),"r"(saida)
        );
}

/*
============
SSE_SinCossincos
\param a ângulo em radianos
\param sin seno do angulo
\param cos cosseno do ângulo
============
*/
void  VPCALL CSIMD_SSE::sincos( float a, float &sin, float &cos ) {
	asm (
		"mov		edi, %18\n"
		"mov		esi, %19\n"
		"movss		xmm1, %17\n"
		"movss		xmm2, xmm1\n"
		"movss		xmm3, xmm1\n"
		"mulss		xmm2, %0\n"
		"cvttss2si	ecx, xmm2\n"
		"cmpltss	xmm3, %1\n"
		"andps		xmm3, %2\n"
		"cvtsi2ss	xmm2, ecx\n"
		"subss		xmm2, xmm3\n"
		"mulss		xmm2, %3\n"
		"subss		xmm1, xmm2\n"

		"movss		xmm0, %4\n"			// xmm0 = PI
		"subss		xmm0, xmm1\n"					// xmm0 = PI - a
		"movss		xmm1, xmm0\n"					// xmm1 = PI - a
		"andps		xmm1, %5\n"	// xmm1 = signbit( PI - a )
		"movss		xmm2, xmm0\n"					// xmm2 = PI - a
		"xorps		xmm2, xmm1\n"					// xmm2 = fabs( PI - a )
		"cmpnltss	xmm2, %6\n"		// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? 0xFFFFFFFF : 0x00000000
		"movss		xmm3, %4\n"			// xmm3 = PI
		"xorps		xmm3, xmm1\n"					// xmm3 = PI ^ signbit( PI - a )
		"andps		xmm3, xmm2\n"					// xmm3 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? ( PI ^ signbit( PI - a ) ) : 0.0f
		"andps		xmm2, %5\n"	// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? SIMD_SP_signBitMask : 0.0f
		"xorps		xmm0, xmm2\n"
		"addps		xmm0, xmm3\n"

		"movss		xmm1, xmm0\n"
		"mulss		xmm1, xmm1\n"
		"movss		xmm3, %12\n"
		"movss		xmm4, %7\n"
		"mulss		xmm3, xmm1\n"
		"mulss		xmm4, xmm1\n"
		"addss		xmm3, %13\n"
		"addss		xmm4, %8\n"
		"mulss		xmm3, xmm1\n"
		"mulss		xmm4, xmm1\n"
		"addss		xmm3, %14\n"
		"addss		xmm4, %9\n"
		"mulss		xmm3, xmm1\n"
		"mulss		xmm4, xmm1\n"
		"addss		xmm3, %15\n"
		"addss		xmm4, %10\n"
		"mulss		xmm3, xmm1\n"
		"mulss		xmm4, xmm1\n"
		"addss		xmm3, %16\n"
		"addss		xmm4, %11\n"
		"mulss		xmm3, xmm1\n"
		"mulss		xmm4, xmm1\n"
		"addss		xmm3, %2\n"
		"addss		xmm4, %2\n"
		"mulss		xmm3, xmm0\n"
		"xorps		xmm2, %5\n"
		"xorps		xmm4, xmm2\n"
		"movss		[edi], xmm2\n"
		"movss		[esi], xmm3\n"
	:
	    : "m" (SIMD_SP_oneOverTwoPI), "m"(SIMD_SP_zero),"m"(SIMD_SP_one),"m"(SIMD_SP_twoPI),
          "m"(SIMD_SP_PI),"m"(SIMD_SP_signBitMask), "m"(SIMD_SP_halfPI),
          "m"(SIMD_SP_cos_c0),"m"(SIMD_SP_cos_c1),"m"(SIMD_SP_cos_c2),
          "m"(SIMD_SP_cos_c3),"m"(SIMD_SP_cos_c4),
          "m"(SIMD_SP_sin_c0),"m"(SIMD_SP_sin_c1),"m"(SIMD_SP_sin_c2),
          "m"(SIMD_SP_sin_c3),"m"(SIMD_SP_sin_c4),"m"(a),"r"(sin),"r"(cos)
);
}

/*
============
SSE_SinCos4
============
*/
void  VPCALL CSIMD_SSE::sincos4( float a[4], float sin[4], float cos[4] ) {
	asm (
		"mov		eax, %17\n"
		"mov		edi, %18\n"
		"mov		esi, %19\n"
		"movaps		xmm1, [eax]\n"
		"movaps		xmm2, xmm1\n"
		"mulps		xmm2, %0\n"
		"movhlps	xmm3, xmm2\n"
		"cvttss2si	ecx, xmm2\n"
		"cvtsi2ss	xmm2, ecx\n"
		"cvttss2si	edx, xmm3\n"
		"cvtsi2ss	xmm3, edx\n"
		"shufps		xmm2, xmm2, 0x01\n" /*R_SHUFFLEPS( 1, 0, 0, 0 )*/
		"shufps		xmm3, xmm3, 0x01\n" /*R_SHUFFLEPS( 1, 0, 0, 0 )*/
		"cvttss2si	ecx, xmm2\n"
		"cvtsi2ss	xmm2, ecx\n"
		"cvttss2si	edx, xmm3\n"
		"cvtsi2ss	xmm3, edx\n"
		"shufps		xmm2, xmm3, 0x11\n"  /*R_SHUFFLEPS( 1, 0, 1, 0 )*/
		"movaps		xmm3, xmm1\n"
		"cmpltps	xmm3, %1\n"
		"andps		xmm3, %2\n"
		"subps		xmm2, xmm3\n"
		"mulps		xmm2, %3\n"
		"subps		xmm1, xmm2\n"

		"movaps		xmm0, %4\n"			// xmm0 = PI
		"subps		xmm0, xmm1\n"				// xmm0 = PI - a
		"movaps		xmm1, xmm0\n"					// xmm1 = PI - a
		"andps		xmm1, %5\n"	// xmm1 = signbit( PI - a )
		"movaps		xmm2, xmm0\n"					// xmm2 = PI - a
		"xorps		xmm2, xmm1\n"					// xmm2 = fabs( PI - a )
		"cmpnltps	xmm2, %6\n"		// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? 0xFFFFFFFF : 0x00000000
		"movaps		xmm3, %4\n"			// xmm3 = PI
		"xorps		xmm3, xmm1\n"					// xmm3 = PI ^ signbit( PI - a )
		"andps		xmm3, xmm2\n"					// xmm3 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? ( PI ^ signbit( PI - a ) ) : 0.0f
		"andps		xmm2, %5\n"	// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? SIMD_SP_signBitMask : 0.0f
		"xorps		xmm0, xmm2\n"
		"addps		xmm0, xmm3\n"

		"movaps		xmm0, [eax]\n"
		"movaps		xmm1, xmm0\n"
		"mulps		xmm1, xmm1\n"
		"movaps		xmm3, %12\n"
		"movaps		xmm4, %7\n"
		"mulps		xmm3, xmm1\n"
		"mulps		xmm4, xmm1\n"
		"addps		xmm3, %13\n"
		"addps		xmm4, %8\n"
		"mulps		xmm3, xmm1\n"
		"mulps		xmm4, xmm1\n"
		"addps		xmm3, %14\n"
		"addps		xmm4, %9\n"
		"mulps		xmm3, xmm1\n"
		"mulps		xmm4, xmm1\n"
		"addps		xmm3, %15\n"
		"addps		xmm4, %10\n"
		"mulps		xmm3, xmm1\n"
		"mulps		xmm4, xmm1\n"
		"addps		xmm3, %16\n"
		"addps		xmm4, %11\n"
		"mulps		xmm3, xmm1\n"
		"mulps		xmm4, xmm1\n"
		"addps		xmm3, %2\n"
		"addps		xmm4, %2\n"
		"mulps		xmm3, xmm0\n"
		"xorps		xmm2, %5\n"
		"xorps		xmm4, xmm2\n"
		"movaps		[edi], xmm3\n"
		"movaps		[esi], xmm4\n"
	:
	    : "m" (SIMD_SP_oneOverTwoPI), "m"(SIMD_SP_zero),"m"(SIMD_SP_one),"m"(SIMD_SP_twoPI),
          "m"(SIMD_SP_PI),"m"(SIMD_SP_signBitMask), "m"(SIMD_SP_halfPI),
          "m"(SIMD_SP_cos_c0),"m"(SIMD_SP_cos_c1),"m"(SIMD_SP_cos_c2),
          "m"(SIMD_SP_cos_c3),"m"(SIMD_SP_cos_c4),
          "m"(SIMD_SP_sin_c0),"m"(SIMD_SP_sin_c1),"m"(SIMD_SP_sin_c2),
          "m"(SIMD_SP_sin_c3),"m"(SIMD_SP_sin_c4),"r"(a),"r"(sin),"r"(cos)
);
}

/*
============
SSE_ATanPositive

  Both 'x' and 'y' must be positive.
\return aTan
============
*/
float  VPCALL CSIMD_SSE::aTanPositive( float y, float x ) {
#if 1

	float t;

	SMF_ASSERT( y >= 0.0f && x >= 0.0f );

	asm (
		"movss		xmm0, %12\n"
		"movss		xmm3, xmm0\n"
		"movss		xmm1, %13\n"
		"minss		xmm0, xmm1\n"
		"maxss		xmm1, xmm3\n"
		"cmpeqss	xmm3, xmm0\n"
		"rcpss		xmm2, xmm1\n"
		"mulss		xmm1, xmm2\n"
		"mulss		xmm1, xmm2\n"
		"addss		xmm2, xmm2\n"
		"subss		xmm2, xmm1\n"				// xmm2 = 1 / y or 1 / x
		"mulss		xmm0, xmm2\n"				// xmm0 = x / y or y / x
		"movss		xmm1, xmm3\n"
		"andps		xmm1, %1\n"
		"xorps		xmm0, xmm1\n"				// xmm0 = -x / y or y / x
		"andps		xmm3, %2\n"	// xmm3 = HALF_PI or 0.0f
		"movss		xmm1, xmm0\n"
		"mulss		xmm1, xmm1\n"				// xmm1 = s
		"movss		xmm2, %3\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %4\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %5\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %6\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %7\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %8\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %9\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %10\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %0\n"
		"mulss		xmm2, xmm0\n"
		"addss		xmm2, xmm3\n"
		"movss		%11, xmm2\n"
        :
        : "m"(SIMD_SP_one),"m"(SIMD_SP_signBitMask), "m"(SIMD_SP_halfPI),
          "m"(SIMD_SP_atan_c0),"m"(SIMD_SP_atan_c1),"m"(SIMD_SP_atan_c2),
          "m"(SIMD_SP_atan_c3),"m"(SIMD_SP_atan_c4),"m"(SIMD_SP_atan_c5),
          "m"(SIMD_SP_atan_c5), "m"(SIMD_SP_atan_c7),"m"(t),"m"(x),"m"(y)
        );

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
void  VPCALL CSIMD_SSE::aTan4Positive( float y[4], float x[4], float resultado[4] ) {
	asm (
		"mov		esi, %12\n"
		"mov		edi, %13\n"
		"mov		edx, %11\n"
		"movaps		xmm0, [esi]\n"
		"movaps		xmm3, xmm0\n"
		"movaps		xmm1, [edi]\n"
		"minps		xmm0, xmm1\n"
		"maxps		xmm1, xmm3\n"
		"cmpeqps	xmm3, xmm0\n"
		"rcpps		xmm2, xmm1\n"
		"mulps		xmm1, xmm2\n"
		"mulps		xmm1, xmm2\n"
		"addps		xmm2, xmm2\n"
		"subps		xmm2, xmm1\n"				// xmm2 = 1 / y or 1 / x
		"mulps		xmm0, xmm2\n"				// xmm0 = x / y or y / x
		"movaps		xmm1, xmm3\n"
		"andps		xmm1, %1\n"
		"xorps		xmm0, xmm1\n"				// xmm0 = -x / y or y / x
		"andps		xmm3, %2\n"	// xmm3 = HALF_PI or 0.0f
		"movaps		xmm1, xmm0\n"
		"mulps		xmm1, xmm1\n"				// xmm1 = s
		"movaps		xmm2, %3\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %4\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %5\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %6\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %7\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %8\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %9\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %10\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %0\n"
		"mulps		xmm2, xmm0\n"
		"addps		xmm2, xmm3\n"
		"movaps		[edx], xmm2\n"
        :
        : "m"(SIMD_SP_one),"m"(SIMD_SP_signBitMask), "m"(SIMD_SP_halfPI),
          "m"(SIMD_SP_atan_c0),"m"(SIMD_SP_atan_c1),"m"(SIMD_SP_atan_c2),
          "m"(SIMD_SP_atan_c3),"m"(SIMD_SP_atan_c4),"m"(SIMD_SP_atan_c5),
          "m"(SIMD_SP_atan_c5), "m"(SIMD_SP_atan_c7),"r"(resultado),"r"(x),"r"(y)
        );
}

/*
============
SSE_ATan
============
*/
float  VPCALL CSIMD_SSE::atan( float y, float x ) {
#if 1

	float resultado;

	asm (
		"movss		xmm0, %12\n"
		"movss		xmm3, xmm0\n"
		"movss		xmm4, xmm0\n"
		"andps		xmm0, %14\n"
		"movss		xmm1, %13\n"
		"xorps		xmm4, xmm1\n"
		"andps		xmm1, %14\n"
		"andps		xmm4, %1\n"
		"minss		xmm0, xmm1\n"
		"maxss		xmm1, xmm3\n"
		"cmpeqss	xmm3, xmm0\n"
		"rcpss		xmm2, xmm1\n"
		"mulss		xmm1, xmm2\n"
		"mulss		xmm1, xmm2\n"
		"addss		xmm2, xmm2\n"
		"subss		xmm2, xmm1\n"				// xmm2 = 1 / y or 1 / x
		"mulss		xmm0, xmm2\n"				// xmm0 = x / y or y / x
		"xorps		xmm0, xmm4\n"
		"movss		xmm1, xmm3\n"
		"andps		xmm1, %1\n"
		"xorps		xmm0, xmm1\n"				// xmm0 = -x / y or y / x
		"orps		xmm4, %2\n"	// xmm4 = +/- HALF_PI
		"andps		xmm3, xmm4\n"				// xmm3 = +/- HALF_PI or 0.0f
		"movss		xmm1, xmm0\n"
		"mulss		xmm1, xmm1\n"				// xmm1 = s
		"movss		xmm2, %3\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %4\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %5\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %6\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %7\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %8\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %9\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %10\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %0\n"
		"mulss		xmm2, xmm0\n"
		"addss		xmm2, xmm3\n"
		"movss		%11, xmm2\n"
        :
        : "m"(SIMD_SP_one),"m"(SIMD_SP_signBitMask), "m"(SIMD_SP_halfPI),
          "m"(SIMD_SP_atan_c0),"m"(SIMD_SP_atan_c1),"m"(SIMD_SP_atan_c2),
          "m"(SIMD_SP_atan_c3),"m"(SIMD_SP_atan_c4),"m"(SIMD_SP_atan_c5),
          "m"(SIMD_SP_atan_c5), "m"(SIMD_SP_atan_c7),"m"(resultado),"m"(x),"m"(y),
          "m"(SIMD_SP_absMask)
        );

	return resultado;

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
void  VPCALL CSIMD_SSE::aTan4( float y[4], float x[4], float resultado[4] ) {
	asm (
		"mov		esi, %12\n"
		"mov		edi, %13\n"
		"mov		edx, %11\n"
		"movaps		xmm0, [esi]\n"
		"movaps		xmm3, xmm0\n"
		"movaps		xmm4, xmm0\n"
		"andps		xmm0, %14\n"
		"movaps		xmm1, [edi]\n"
		"xorps		xmm4, xmm1\n"
		"andps		xmm1, %14\n"
		"andps		xmm4, %1\n"
		"minps		xmm0, xmm1\n"
		"maxps		xmm1, xmm3\n"
		"cmpeqps	xmm3, xmm0\n"
		"rcpps		xmm2, xmm1\n"
		"mulps		xmm1, xmm2\n"
		"mulps		xmm1, xmm2\n"
		"addps		xmm2, xmm2\n"
		"subps		xmm2, xmm1\n"				// xmm2 = 1 / y or 1 / x
		"mulps		xmm0, xmm2\n"				// xmm0 = x / y or y / x
		"xorps		xmm0, xmm4\n"
		"movaps		xmm1, xmm3\n"
		"andps		xmm1, %1\n"
		"xorps		xmm0, xmm1\n"				// xmm0 = -x / y or y / x
		"orps		xmm4, %2\n"	// xmm4 = +/- HALF_PI
		"andps		xmm3, xmm4\n"				// xmm3 = +/- HALF_PI or 0.0f
		"movaps		xmm1, xmm0\n"
		"mulps		xmm1, xmm1\n"				// xmm1 = s
		"movaps		xmm2, %3\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %4\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %5\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %6\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %7\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %8\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %9\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %10\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %0\n"
		"mulps		xmm2, xmm0\n"
		"addps		xmm2, xmm3\n"
		"movaps		[edx], xmm2\n"
        :
        : "m"(SIMD_SP_one),"m"(SIMD_SP_signBitMask), "m"(SIMD_SP_halfPI),
          "m"(SIMD_SP_atan_c0),"m"(SIMD_SP_atan_c1),"m"(SIMD_SP_atan_c2),
          "m"(SIMD_SP_atan_c3),"m"(SIMD_SP_atan_c4),"m"(SIMD_SP_atan_c5),
          "m"(SIMD_SP_atan_c5), "m"(SIMD_SP_atan_c7),"r"(resultado),"r"(x),"r"(y),
          "m"(SIMD_SP_absMask)
        );
}


/*
============
CSIMD_SSE::add

  dst[i] = constant + src[i];
============
*/
void VPCALL CSIMD_SSE::add( float *dstP, const float constant, const float *srcP, const int count ) {

	//KFLOAT_CA( add, dstP, srcP, constant, count )
	int	pre,post;
    //printf("KFLOAT_CA");
    asm(
        KFLOAT_CA( add )
    :
    : "m"(constant),"m"(count),"m"(dstP),"m"(srcP),"m"(pre),"m"(post)
    : "eax","ecx","ebx","edx");
}

/*
============
CSIMD_SSE::add

  dst[i] = src0[i] + src1[i];
============
*/
void VPCALL CSIMD_SSE::add( float *dst, const float *src0, const float *src1, const int count ) {

	//KFLOAT_AA( add, dst, src0, src1, count )

	int	pre,post;

    asm(
        KFLOAT_AA( add )
    :
    : "m"(count),"m"(dst),"m"(src0),"m"(src1),"m"(pre),"m"(post)
    : "eax","ecx","ebx","edx" );

}


/*
============
CSIMD_SSE::sub

  dst[i] = constant - src[i];
============
*/
void VPCALL CSIMD_SSE::sub( float *dst, const float constant, const float *src, const int count ) {


	//KFLOAT_CA( sub, dst, src, constant, count )
	int	pre,post;
    asm(
        KFLOAT_CA( sub )
    :
    : "m"(constant),"m"(count),"m"(dst),"m"(src),"m"(pre),"m"(post)
    : "eax","ecx","ebx","edx" );

}

/*
============
CSIMD_SSE::sub

  dst[i] = src0[i] - src1[i];
============
*/
void VPCALL CSIMD_SSE::sub( float *dst, const float *src0, const float *src1, const int count ) {
	//KFLOAT_AA( sub, dst, src0, src1, count )
	int	pre,post;
    asm(
        KFLOAT_AA( sub )
    :
    : "m"(count),"m"(dst),"m"(src0),"m"(src1),"m"(pre),"m"(post)
    : "eax","ecx","ebx","edx" );
}

/*
============
CSIMD_SSE::mul

  dst[i] = constant * src[i];
============
*/
void VPCALL CSIMD_SSE::mul( float *dst, const float constant, const float *src, const int count ) {
	//KFLOAT_CA( mul, dst, src, constant, count )

	int	pre,post;
    asm(
        KFLOAT_CA( mul )
    :
    : "m"(constant),"m"(count),"m"(dst),"m"(src),"m"(pre),"m"(post)
    : "eax","ecx","ebx","edx" );
}

/*
============
CSIMD_SSE::mul

  dst[i] = src0[i] * src1[i];
============
*/
void VPCALL CSIMD_SSE::mul( float *dst, const float *src0, const float *src1, const int count ) {
	//KFLOAT_AA( mul, dst, src0, src1, count )
	int	pre,post;
    asm(
        KFLOAT_AA( mul )
    :
    : "m"(count),"m"(dst),"m"(src0),"m"(src1),"m"(pre),"m"(post)
    : "eax","ecx","ebx","edx" );
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
	asm(
		"movss	xmm1,%0\n"
		"shufps	xmm1,xmm1,0\n"

		KFLOATINITDS( %2, %3, %1, %4, %5 )
		"and		eax,15\n"
		"jne		lpNA%=\n"
		"jmp		lpA%=\n"
		/*align	16*/
"lpA%=:\n"
		"movaps	xmm2,[edx+ebx]\n"
		"movaps	xmm3,[edx+ebx+16]\n"
		"rcpps	xmm4,xmm2\n"
		"rcpps	xmm5,xmm3\n"
		"prefetchnta	[edx+ebx+64]\n"
		"mulps	xmm2,xmm4\n"
		"mulps	xmm2,xmm4\n"
		"mulps	xmm3,xmm5\n"
		"mulps	xmm3,xmm5\n"
		"addps	xmm4,xmm4\n"
		"addps	xmm5,xmm5\n"
		"subps	xmm4,xmm2\n"
		"subps	xmm5,xmm3\n"
		"mulps	xmm4,xmm1\n"
		"mulps	xmm5,xmm1\n"
		"movaps	[edi+ebx],xmm4\n"
		"movaps	[edi+ebx+16],xmm5\n"
		"add		ebx,16*2\n"
		"jl		lpA%=\n"
		"jmp		done%=\n"
		/*align	16*/
"lpNA%=:\n"
		"movups	xmm2,[edx+ebx]\n"
		"movups	xmm3,[edx+ebx+16]\n"
		"rcpps	xmm4,xmm2\n"
		"rcpps	xmm5,xmm3\n"
		"prefetchnta	[edx+ebx+64]\n"
		"mulps	xmm2,xmm4\n"
		"mulps	xmm2,xmm4\n"
		"mulps	xmm3,xmm5\n"
		"mulps	xmm3,xmm5\n"
		"addps	xmm4,xmm4\n"
		"addps	xmm5,xmm5\n"
		"subps	xmm4,xmm2\n"
		"subps	xmm5,xmm3\n"
		"mulps	xmm4,xmm1\n"
		"mulps	xmm5,xmm1\n"
		"movaps	[edi+ebx],xmm4\n"
		"movaps	[edi+ebx+16],xmm5\n"
		"add	ebx,16*2\n"
		"jl		lpNA%=\n"
"done%=:\n"
		"mov		edx,%3\n"
		"mov		edi,%2\n"
		KFLOATOPER( KDIVDSS1( [edi+ebx],xmm1,[edx+ebx] ),
					KDIVDSS4( [edi+ebx],xmm1,[edx+ebx] ), %1,%4,%5 )
	:
    :"m"(constant),"m"(count),"m"(dst),"m"(src),"m"(pre),"m"(post)
    : "eax","ecx","ebx","edx","esi","edi");
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
	asm(
		KFLOATINITDSS( %1, %2, %3, %0,%4,%5 )
		"and		eax,15\n"
		"jne		lpNA%=\n"
		"jmp		lpA%=\n"
		/*align	16*/
"lpA%=:\n"
		"movaps	xmm2,[esi+ebx]\n"
		"movaps	xmm3,[esi+ebx+16]\n"
		"rcpps	xmm4,xmm2\n"
		"rcpps	xmm5,xmm3\n"
		"prefetchnta	[esi+ebx+64]\n"
		"mulps	xmm2,xmm4\n"
		"mulps	xmm2,xmm4\n"
		"mulps	xmm3,xmm5\n"
		"mulps	xmm3,xmm5\n"
		"addps	xmm4,xmm4\n"
		"addps	xmm5,xmm5\n"
		"subps	xmm4,xmm2\n"
		"subps	xmm5,xmm3\n"
		"mulps	xmm4,[edx+ebx]\n"
		"mulps	xmm5,[edx+ebx+16]\n"
		"movaps	[edi+ebx],xmm4\n"
		"movaps	[edi+ebx+16],xmm5\n"
		"add	ebx,16*2\n"
		"jl		lpA%=\n"
		"jmp		done%=\n"
		/*align	16*/
"lpNA%=:\n"
		"movups	xmm2,[esi+ebx]\n"
		"movups	xmm3,[esi+ebx+16]\n"
		"rcpps	xmm4,xmm2\n"
		"rcpps	xmm5,xmm3\n"
		"prefetchnta	[esi+ebx+64]\n"
		"mulps	xmm2,xmm4\n"
		"mulps	xmm2,xmm4\n"
		"mulps	xmm3,xmm5\n"
		"mulps	xmm3,xmm5\n"
		"addps	xmm4,xmm4\n"
		"addps	xmm5,xmm5\n"
		"subps	xmm4,xmm2\n"
		"subps	xmm5,xmm3\n"
		"movups	xmm2,[edx+ebx]\n"
		"movups	xmm3,[edx+ebx+16]\n"
		"mulps	xmm4,xmm2\n"
		"mulps	xmm5,xmm3\n"
		"movaps	[edi+ebx],xmm4\n"
		"movaps	[edi+ebx+16],xmm5\n"
		"add	ebx,16*2\n"
		"jl		lpNA%=\n"
"done%=:\n"
		"mov		edx,%2\n"
		"mov		esi,%3\n"
		"mov		edi,%1\n"
		KFLOATOPER( KDIVDSS1( [edi+ebx],[edx+ebx],[esi+ebx] ),
					KDIVDSS4( [edi+ebx],[edx+ebx],[esi+ebx] ), %0,%4,%5 )
        :
	    :"m"(count),"m"(dst), "m"(src0), "m"(src1),"m"(pre),"m"(post)
	        : "eax","ecx","ebx","edx");
}

/*
============
Simd_MulAdd

 assumes count >= 7
============
*/
static void Simd_MulAdd( float *dst, const float constant, const float *src, const int count ) {
	asm(
    "mov			esi, %3\n"
	"mov			edi, %2\n"
	"mov			eax, %0\n"
	"shl			eax, 2\n"
	"mov			ecx, esi\n"
	"mov			edx, eax\n"
	"or			ecx, edi\n"
    "fld			%1\n"
	"and			ecx, 15\n"
	"jz			SimdMuladd16\n"
	"and			ecx, 3\n"
	"jnz			SimdMulAdd8\n"
	"mov			ecx, esi\n"
	"xor			ecx, edi\n"
	"and			ecx, 15\n"
	"jnz			MulAdd8\n"
	"mov			ecx, esi\n"
	"and			ecx, 15\n"
	"neg			ecx\n"
	"add			ecx, 16\n"
	"sub			eax, ecx\n"
	"add			edi, ecx\n"
	"add			esi, ecx\n"
	"neg			ecx\n"
	"mov			edx, eax\n"
	"loopPreMuladd16:\n"
	"fld			st\n"
	"fmul		dword ptr [edi+ecx]\n"
	"fadd		dword ptr [esi+ecx]\n"
	"fstp		dword ptr [esi+ecx]\n"
	"add			ecx, 4\n"
	"jl			loopPreMuladd16\n"
	"SimdMuladd16:\n"
	"and			eax, ~15\n"
	"movss		xmm1, %1\n"
	"shufps		xmm1, xmm1, 0x00\n"
	"add			esi, eax\n"
	"add			edi, eax\n"
	"neg			eax\n"
	/*align		16*/
	"loopMuladd16:\n"
	"movaps		xmm0, [edi+eax]\n"
	"mulps		xmm0, xmm1\n"
	"addps		xmm0, [esi+eax]\n"
	"movaps		[esi+eax], xmm0\n"
	"add			eax, 16\n"
	"jl			loopMuladd16\n"
	"jmp			postMulAdd\n"
	"MulAdd8:\n"
	"mov			ecx, esi\n"
	"and			ecx, 7\n"
	"jz			SimdMulAdd8\n"
	"sub			eax, ecx\n"
	"add			esi, ecx\n"
	"add			edi, ecx\n"
	"neg			ecx\n"
	"mov			edx, eax\n"
	"loopPreMulAdd8:\n"
	"fld			st\n"
	"fmul		dword ptr [edi+ecx]\n"
	"fadd		dword ptr [esi+ecx]\n"
	"fstp		dword ptr [esi+ecx]\n"
	"add		ecx, 4\n"
	"jl			loopPreMulAdd8\n"
	"SimdMulAdd8:\n"
	"and		eax, ~15\n"
	"movss		xmm1, %1\n"
	"shufps		xmm1, xmm1, 0x00\n"
	"add			esi, eax\n"
	"add			edi, eax\n"
	"neg			eax\n"
	/*align		16*/
	"loopMulAdd8:\n"
	"movlps		xmm0, [edi+eax]\n"
	"movhps		xmm0, [edi+eax+8]\n"
	"mulps		xmm0, xmm1\n"
	"movlps		xmm2, [esi+eax]\n"
	"movhps		xmm2, [esi+eax+8]\n"
	"addps		xmm0, xmm2\n"
	"movlps		[esi+eax], xmm0\n"
	"movhps		[esi+eax+8], xmm0\n"
	"add			eax, 16\n"
	"jl			loopMulAdd8\n"
	"jmp			postMulAdd\n"
	"postMulAdd:\n"
	"and			edx, 15\n"
	"jz			MulAddDone\n"
	"add			esi, edx\n"
	"add			edi, edx\n"
	"neg			edx\n"
	"loopPostMulAdd:\n"
	"fld			st\n"
	"fmul		dword ptr [edi+edx]\n"
	"fadd		dword ptr [esi+edx]\n"
	"fstp		dword ptr [esi+edx]\n"
	"add		edx, 4\n"
	"jl			loopPostMulAdd\n"
	"MulAddDone:\n"
	"fstp		st\n"
    :
    :"m"(count),"m"(constant),"r"(src),"r"(dst)
    :);
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

    float * teste = constant.toFloatPtr();
	asm(
		"mov			eax, %0\n"  //count
		"mov			edi, %1\n"  //constant_Ptr
		"mov			edx, eax\n"
		"mov			esi, %2\n"  //src
		"mov			ecx, %3\n"   //dst
		"and			eax, ~3\n"

		"movss		xmm4, [edi+0]\n"
		"shufps		xmm4, xmm4, 0x00\n"
		"movss		xmm5, [edi+4]\n"
		"shufps		xmm5, xmm5, 0x00\n"
		"movss		xmm6, [edi+8]\n"
		"shufps		xmm6, xmm6, 0x00\n"

		"jz			done4%=\n"
		"imul		eax, 12\n"
		"add		esi, eax\n"
		"neg		eax\n"

	"loop4%=:\n"
		"movlps		xmm1, [esi+eax+ 0]\n"
		"movlps		xmm2, [esi+eax+ 8]\n"
		"movlps		xmm3, [esi+eax+16]\n"
		"movhps		xmm1, [esi+eax+24]\n"
		"movhps		xmm2, [esi+eax+32]\n"
		"movhps		xmm3, [esi+eax+40]\n"
		"movaps		xmm0, xmm1\n"
		"shufps		xmm0, xmm2, 0xd8\n" /*R_SHUFFLEPS( 0, 2, 1, 3 )*/
		"shufps		xmm1, xmm3, 0x8d\n" /*R_SHUFFLEPS( 1, 3, 0, 2 )*/
		"shufps		xmm2, xmm3, 0xd8\n"  /*R_SHUFFLEPS( 0, 2, 1, 3 )*/
		"add		ecx, 16\n"
		"add		eax, 4*12\n"
		"mulps		xmm0, xmm4\n"
		"mulps		xmm1, xmm5\n"
		"mulps		xmm2, xmm6\n"
		"addps		xmm0, xmm1\n"
		"addps		xmm0, xmm2\n"
		"shufps		xmm0, xmm0, 0xd8\n" /*R_SHUFFLEPS( 0, 2, 1, 3 )*/
		"movlps		[ecx-16+0], xmm0\n"
		"movhps		[ecx-16+8], xmm0\n"
		"jl			loop4%=\n"

	"done4%=:\n"
		"and		edx, 3\n"
		"jz			done1%=\n"

	"loop1%=:\n"
		"movss		xmm0, [esi+eax+0]\n"
		"movss		xmm1, [esi+eax+4]\n"
		"movss		xmm2, [esi+eax+8]\n"
		"mulss		xmm0, xmm4\n"
		"mulss		xmm1, xmm5\n"
		"mulss		xmm2, xmm6\n"
		"add		ecx, 4\n"
		"addss		xmm0, xmm1\n"
		"add		eax, 12\n"
		"addss		xmm0, xmm2\n"
		"dec		edx\n"
		"movss		[ecx-4], xmm0\n"
		"jnz			loop1%=\n"

	"done1%=:\n"
	:
    :"m"(count),"m"(teste),"m"(src),"m"(dst)
    :"edi","eax","esi");

}

/*
============
CSIMD_SSE::dot

  dst[i] = constant * src[i].Normal() + src[i][3];
============
*/
void VPCALL CSIMD_SSE::dot( float *dst, const CVec3D &constant, const CPlane *src, const int count ) {
    float * teste = constant.toFloatPtr();

	asm (
		"mov			eax, %3\n"
		"mov			edi, %1\n"
		"mov			edx, eax\n"
		"mov			esi, %2\n"
		"mov			ecx, %0\n"
		"and			eax, ~3\n"

		"movss		xmm5, [edi+0]\n"
		"shufps		xmm5, xmm5, 0x00\n"
		"movss		xmm6, [edi+4]\n"
		"shufps		xmm6, xmm6, 0x00\n"
		"movss		xmm7, [edi+8]\n"
		"shufps		xmm7, xmm7, 0x00\n"

		"jz			startVert1\n"
		"imul		eax, 16\n"
		"add		esi, eax\n"
		"neg		eax\n"

	"loopVert4:\n"

		"movlps		xmm1, [esi+eax+ 0]\n"
		"movlps		xmm3, [esi+eax+ 8]\n"
		"movhps		xmm1, [esi+eax+16]\n"
		"movhps		xmm3, [esi+eax+24]\n"
		"movlps		xmm2, [esi+eax+32]\n"
		"movlps		xmm4, [esi+eax+40]\n"
		"movhps		xmm2, [esi+eax+48]\n"
		"movhps		xmm4, [esi+eax+56]\n"
		"movaps		xmm0, xmm1\n"
		"shufps		xmm0, xmm2, 0x88 \n" /*R_SHUFFLEPS( 0, 2, 0, 2 )*/
		"shufps		xmm1, xmm2, 0xdd\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )*/
		"movaps		xmm2, xmm3\n"
		"shufps		xmm2, xmm4, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )*/
		"shufps		xmm3, xmm4, 0xdd\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )=11011101*/

		"add		ecx, 16\n"
		"add		eax, 4*16\n"

		"mulps		xmm0, xmm5\n"
		"mulps		xmm1, xmm6\n"
		"mulps		xmm2, xmm7\n"
		"addps		xmm0, xmm3\n"
		"addps		xmm0, xmm1\n"
		"addps		xmm0, xmm2\n"

		"movlps		[ecx-16+0], xmm0\n"
		"movhps		[ecx-16+8], xmm0\n"
		"jl			loopVert4\n"

	"startVert1:\n"
		"and		edx, 3\n"
		"jz			done\n"

	"loopVert1:\n"
		"movss		xmm0, [esi+eax+0]\n"
		"movss		xmm1, [esi+eax+4]\n"
		"movss		xmm2, [esi+eax+8]\n"
		"mulss		xmm0, xmm5\n"
		"mulss		xmm1, xmm6\n"
		"mulss		xmm2, xmm7\n"
		"addss		xmm0, [esi+eax+12]\n"
		"add		ecx, 4\n"
		"addss		xmm0, xmm1\n"
		"add		eax, 16\n"
		"addss		xmm0, xmm2\n"
		"dec		edx\n"
		"movss		[ecx-4], xmm0\n"
		"jnz		loopVert1\n"

	"done:\n"
	:
    :"m"(dst), "m"(teste), "m"(src), "m"(count)
    :"edi","eax","esi");

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
    float * constPtr = constant.toFloatPtr();

    // 0,  1,  2
	// 3,  4,  5
	// 6,  7,  8
	// 9, 10, 11

	asm (
		"mov			eax, %3\n"
		"mov			edi, %1\n"
		"mov			edx, eax\n"
		"mov			esi, %2\n"
		"mov			ecx, %0\n"
		"and			eax, ~3\n"

		"movss		xmm4, [edi+0]\n"
		"shufps		xmm4, xmm4, 0x00\n"
		"movss		xmm5, [edi+4]\n"
		"shufps		xmm5, xmm5, 0x00\n"
		"movss		xmm6, [edi+8]\n"
		"shufps		xmm6, xmm6, 0x00\n"

		"jz			startVert1\n"
		"imul		eax, "DRAWVERT_SIZE_STR"\n"
		"add		esi, eax\n"
		"neg		eax\n"

	"loopVert4%=:\n"
		"movss		xmm0, [esi+eax+1*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"	//  3,  X,  X,  X
		"movss		xmm2, [esi+eax+0*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"	//  2,  X,  X,  X
		"movhps		xmm0, [esi+eax+0*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"	//  3,  X,  0,  1
		"movaps		xmm1, xmm0\n"												//  3,  X,  0,  1

		"movlps		xmm1, [esi+eax+1*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"	//  4,  5,  0,  1
		"shufps		xmm2, xmm1, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/			//  2,  X,  4,  5

		"movss		xmm3, [esi+eax+3*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"	//  9,  X,  X,  X
		"movhps		xmm3, [esi+eax+2*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"	//  9,  X,  6,  7
		"shufps		xmm0, xmm3, 0x22\n" /*R_SHUFFLEPS( 2, 0, 2, 0 )=00100010*/					//  0,  3,  6,  9

		"movlps		xmm3, [esi+eax+3*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"	// 10, 11,  6,  7
		"shufps		xmm1, xmm3, 0x33\n" /*R_SHUFFLEPS( 3, 0, 3, 0 )=00110011*/					//  1,  4,  7, 10

		"movhps		xmm3, [esi+eax+2*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"	// 10, 11,  8,  X
		"shufps		xmm2, xmm3, 0x6c\n" /*R_SHUFFLEPS( 0, 3, 2, 1 )=01101100*/					//  2,  5,  8, 11

		"add		ecx, 16\n"
		"add		eax, 4 *" DRAWVERT_SIZE_STR "\n"

		"mulps		xmm0, xmm4\n"
		"mulps		xmm1, xmm5\n"
		"mulps		xmm2, xmm6\n"
		"addps		xmm0, xmm1\n"
		"addps		xmm0, xmm2\n"

		"movlps		[ecx-16+0], xmm0\n"
		"movhps		[ecx-16+8], xmm0\n"
		"jl			loopVert4%=\n"

	"startVert1%=:\n"
		"and		edx, 3\n"
		"jz			done%=\n"

	"loopVert1%=:\n"
		"movss		xmm0, [esi+eax+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"movss		xmm1, [esi+eax+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"movss		xmm2, [esi+eax+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"mulss		xmm0, xmm4\n"
		"mulss		xmm1, xmm5\n"
		"mulss		xmm2, xmm6\n"
		"add		ecx, 4\n"
		"addss		xmm0, xmm1\n"
		"add		eax, "DRAWVERT_SIZE_STR"\n"
		"addss		xmm0, xmm2\n"
		"dec		edx\n"
		"movss		[ecx-4], xmm0\n"
		"jnz		loopVert1%=\n"

	"done%=:\n"
	:
    :"m"(dst), "m"(constPtr), "m"(src), "m"(count)
    :"esi","edi","eax");


}

/*
============
CSIMD_SSE::dot

  dst[i] = constant.Normal() * src[i] + constant[3];
============
*/
void VPCALL CSIMD_SSE::dot( float *dst, const CPlane &constant, const CVec3D *src, const int count ) {

    float * constPtr = constant.toFloatPtr();

	asm(
		"mov			eax, %3\n"
		"mov			edi, %1\n"
		"mov			edx, eax\n"
		"mov			esi, %2\n"
		"mov			ecx, %0\n"
		"and			eax, ~3\n"

		"movss		xmm4, [edi+0]\n"
		"shufps		xmm4, xmm4, 0x00\n"
		"movss		xmm5, [edi+4]\n"
		"shufps		xmm5, xmm5, 0x00\n"
		"movss		xmm6, [edi+8]\n"
		"shufps		xmm6, xmm6, 0x00\n"
		"movss		xmm7, [edi+12]\n"
		"shufps		xmm7, xmm7, 0x00\n"

		"jz			done4\n"
		"imul		eax, 12\n"
		"add		esi, eax\n"
		"neg		eax\n"

	"loop4%=:\n"
		"movlps		xmm1, [esi+eax+ 0]\n"
		"movlps		xmm2, [esi+eax+ 8]\n"
		"movlps		xmm3, [esi+eax+16]\n"
		"movhps		xmm1, [esi+eax+24]\n"
		"movhps		xmm2, [esi+eax+32]\n"
		"movhps		xmm3, [esi+eax+40]\n"
		"movaps		xmm0, xmm1\n"
		"shufps		xmm0, xmm2, 0xd8\n" /*R_SHUFFLEPS( 0, 2, 1, 3 )=11011000*/
		"shufps		xmm1, xmm3, 0x8d\n"  /*R_SHUFFLEPS( 1, 3, 0, 2 )*/
		"shufps		xmm2, xmm3, 0xd8\n" /*R_SHUFFLEPS( 0, 2, 1, 3 )=11011000*/

		"add			ecx, 16\n"
		"add			eax, 4*12\n"

		"mulps		xmm0, xmm4\n"
		"mulps		xmm1, xmm5\n"
		"mulps		xmm2, xmm6\n"
		"addps		xmm0, xmm7\n"
		"addps		xmm0, xmm1\n"
		"addps		xmm0, xmm2\n"
		"shufps		xmm0, xmm0, 0xd8 \n" /*R_SHUFFLEPS( 0, 2, 1, 3 )=11011000*/

		"movlps		[ecx-16+0], xmm0\n"
		"movhps		[ecx-16+8], xmm0\n"
		"jl			loop4%=\n"

	"done4:%=:\n"
		"and		edx, 3\n"
		"jz			done1%=\n"

	"loop1%=:\n"
		"movss		xmm0, [esi+eax+0]\n"
		"movss		xmm1, [esi+eax+4]\n"
		"movss		xmm2, [esi+eax+8]\n"
		"mulss		xmm0, xmm4\n"
		"mulss		xmm1, xmm5\n"
		"mulss		xmm2, xmm6\n"
		"addss		xmm0, xmm7\n"
		"add		ecx, 4\n"
		"addss		xmm0, xmm1\n"
		"add		eax, 12\n"
		"addss		xmm0, xmm2\n"
		"dec		edx\n"
		"movss		[ecx-4], xmm0\n"
		"jnz			loop1%=\n"
	"done1%=:\n"
	:
    :"m"(dst), "m"(constPtr), "m"(src), "m"(count)
    :);

}

/*
============
CSIMD_SSE::dot

  dst[i] = constant.Normal() * src[i].Normal() + constant[3] * src[i][3];
============
*/
void VPCALL CSIMD_SSE::dot( float *dst, const CPlane &constant, const CPlane *src, const int count ) {
    float * constPtr = constant.toFloatPtr();

#define SINGLE_OP(SRC, DEST)							\
	"movlps		xmm0,["STRINGIZE(SRC)"]\n"						\
	"movlps		xmm1,["STRINGIZE(SRC)"+8]\n"					\
	"mulps		xmm0,xmm4\n"						\
	"mulps		xmm1,xmm5\n"						\
	"addps		xmm0,xmm1\n"						\
	"movaps		xmm1,xmm0\n"						\
	"shufps		xmm1,xmm1, 0x55 \n" /*SHUFFLEPS(1,1,1,1) */	\
	"addss		xmm0,xmm1\n"						\
	"movss		["STRINGIZE(DEST)"],xmm0\n"						\
	"add		"STRINGIZE(SRC)", 16 \n"							\
	"add		"STRINGIZE(DEST)", 4 \n"

#define DUAL_OP(SRC, DEST)								\
	"movlps		xmm0,["STRINGIZE(SRC)"]\n"						\
	"movlps		xmm1,["STRINGIZE(SRC)"+8]\n"					\
	"movhps		xmm0,["STRINGIZE(SRC)"+16]\n"					\
	"movhps		xmm1,["STRINGIZE(SRC)"+24]\n"					\
	"mulps		xmm0,xmm4\n"						\
	"mulps		xmm1,xmm5\n"						\
	"addps		xmm0,xmm1\n"						\
	"shufps		xmm1,xmm0, 0x12\n" 	\
	"shufps		xmm0,xmm0, 0x27\n"  	\
	"addps		xmm0,xmm1\n"						\
	"movhps		["STRINGIZE(DEST)"],xmm0\n"						\
	"add			"STRINGIZE(SRC)",32\n"							\
	"add			"STRINGIZE(DEST)",8\n"

	asm (
		"mov			edx, %0\n"   //dst
		"mov			eax, %2\n"  //SRC
		"mov			ebx, %1\n"  //constant
		"mov			ecx, %3\n"  //count

		"movlps		xmm4, [ebx]\n"
		"shufps		xmm4, xmm4, 0x44 \n"  /*SHUFFLEPS(1,0,1,0)=01000100 */
		"movlps		xmm5, [ebx+8]\n"
		"shufps		xmm5, xmm5, 0x44 \n" /*SHUFFLEPS(1,0,1,0)=01000100 */

		"xorps		xmm0, xmm0\n"
		"xorps		xmm1, xmm1\n"

	"_lpAlignDest%=:\n"
		"test		edx, 0x0f\n"
		"jz			_destAligned%=\n"
		SINGLE_OP(eax,edx)
 		"dec			ecx\n"
		"jnz			_lpAlignDest%=\n"
		"jmp			_vpExit%=\n"

	"_destAligned%=:\n"
		"push		ecx\n"

		"cmp			ecx, 4\n"
		"jl			_post%=\n"

		"and			ecx, ~3\n"
		"shl			ecx, 2\n"
		"lea			eax, [eax+ecx*4]\n"
		"add			edx, ecx\n"
		"neg			ecx\n"

		"movlps		xmm0, [eax+ecx*4]\n"
		"movhps		xmm0, [eax+ecx*4+16]\n"
		"movlps		xmm2, [eax+ecx*4+32]\n"
		"movhps		xmm2, [eax+ecx*4+48]\n"
		"jmp		_lpStart%=\n"

		/*align	16*/
	"_lp%=:\n"
		"prefetchnta	[eax+ecx*4+128]\n"
		"addps		xmm1, xmm0\n"
		"movlps		xmm0, [eax+ecx*4]\n"
		"movhps		xmm0, [eax+ecx*4+16]\n"
		"movlps		xmm2, [eax+ecx*4+32]\n"
		"movhps		xmm2, [eax+ecx*4+48]\n"
		"movaps		[edx+ecx-16],xmm1\n"
	"_lpStart%=:\n"
		"movlps		xmm1, [eax+ecx*4+8]\n"
		"movhps		xmm1, [eax+ecx*4+24]\n"
		"movlps		xmm3, [eax+ecx*4+40]\n"
		"movhps		xmm3, [eax+ecx*4+56]\n"
		"add		ecx, 16\n"
		"mulps		xmm1, xmm5\n"
		"mulps		xmm2, xmm4\n"
		"mulps		xmm3, xmm5\n"
		"addps		xmm2, xmm3\n"						// y3+w3 x3+z3 y2+w2 x2+z2
		"mulps		xmm0, xmm4\n"
		"addps		xmm0, xmm1\n"						// y1+w1 x1+z1 y0+w0 x0+z0
		"movaps		xmm1, xmm0\n"
		"shufps		xmm0, xmm2, 0x88\n" /*SHUFFLEPS(2,0,2,0)=10001000 */	// x3+z3 x2+z2 x1+z1 x0+z0
		"shufps		xmm1, xmm2, 0xdd\n"  /*SHUFFLEPS(3,1,3,1)= 11011101*/	// y3+w3 y2+w2 y1+w1 y0+w0
		"js			_lp%=\n"
		"addps		xmm1, xmm0\n"
		"movaps		[edx+ecx-16], xmm1\n"
	"_post%=:\n"
		"pop			ecx\n"
		"and			ecx, 0x3\n"
		"cmp			ecx, 2\n"
		"jl			_post1%=\n"
        DUAL_OP(eax,edx)
        "sub			ecx, 2\n"
	"_post1%=:\n"
		"cmp			ecx, 1\n"
		"jne			_vpExit%=\n"
		SINGLE_OP(eax,edx)
	"_vpExit%=:\n"
	:
    :"m"(dst), "m"(constPtr), "m"(src), "m"(count)
    :"esi","edi","eax","ebx","ecx");

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
   float * constPtr = constant.toFloatPtr();

	// 0,  1,  2
	// 3,  4,  5
	// 6,  7,  8
	// 9, 10, 11

	asm (
		"mov			eax, %3\n"
		"mov			edi, %1\n"
		"mov			edx, eax\n"
		"mov			esi, %2\n"
		"mov			ecx, %0\n"
		"and			eax, ~3\n"

		"movss		xmm4, [edi+0]\n"
		"shufps		xmm4, xmm4, 0x00\n"
		"movss		xmm5, [edi+4]\n"
		"shufps		xmm5, xmm5, 0x00\n"
		"movss		xmm6, [edi+8]\n"
		"shufps		xmm6, xmm6, 0x00\n"
		"movss		xmm7, [edi+12]\n"
		"shufps		xmm7, xmm7, 0x00\n"

		"jz			startVert1%=\n"
		"imul		eax, "DRAWVERT_SIZE_STR"\n"
		"add			esi, eax\n"
		"neg			eax\n"

	"loopVert4%=:\n"
		"movss		xmm0, [esi+eax+1*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"	//  3,  X,  X,  X
		"movss		xmm2, [esi+eax+0*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"	//  2,  X,  X,  X
		"movhps		xmm0, [esi+eax+0*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"	//  3,  X,  0,  1
		"movaps		xmm1, xmm0\n"												//  3,  X,  0,  1

		"movlps		xmm1, [esi+eax+1*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"	//  4,  5,  0,  1
		"shufps		xmm2, xmm1, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/					//  2,  X,  4,  5

		"movss		xmm3, [esi+eax+3*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"	//  9,  X,  X,  X
		"movhps		xmm3, [esi+eax+2*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"	//  9,  X,  6,  7
		"shufps		xmm0, xmm3, 0x22\n" /*R_SHUFFLEPS( 2, 0, 2, 0 )=00100010*/					//  0,  3,  6,  9

		"movlps		xmm3, [esi+eax+3*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"	// 10, 11,  6,  7
		"shufps		xmm1, xmm3, 0x33\n" /*R_SHUFFLEPS( 3, 0, 3, 0 )=00110011*/					//  1,  4,  7, 10

		"movhps		xmm3, [esi+eax+2*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"	// 10, 11,  8,  X
		"shufps		xmm2, xmm3, 0x6c\n" /*R_SHUFFLEPS( 0, 3, 2, 1 )=01101100*/					//  2,  5,  8, 11

		"add			ecx, 16\n"
		"add			eax, 4*"DRAWVERT_SIZE_STR"\n"

		"mulps		xmm0, xmm4\n"
		"mulps		xmm1, xmm5\n"
		"mulps		xmm2, xmm6\n"
		"addps		xmm0, xmm7\n"
		"addps		xmm0, xmm1\n"
		"addps		xmm0, xmm2\n"

		"movlps		[ecx-16+0], xmm0\n"
		"movhps		[ecx-16+8], xmm0\n"
		"jl			loopVert4%=\n"

	"startVert1%=:\n"
		"and		edx, 3\n"
		"jz			done%=\n"

	"loopVert1%=:\n"
		"movss		xmm0, [esi+eax+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"movss		xmm1, [esi+eax+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"movss		xmm2, [esi+eax+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"mulss		xmm0, xmm4\n"
		"mulss		xmm1, xmm5\n"
		"mulss		xmm2, xmm6\n"
		"addss		xmm0, xmm7\n"
		"add		ecx, 4\n"
		"addss		xmm0, xmm1\n"
		"add		eax, "DRAWVERT_SIZE_STR"\n"
		"addss		xmm0, xmm2\n"
		"dec		edx\n"
		"movss		[ecx-4], xmm0\n"
		"jnz		loopVert1%=\n"

	"done%=:\n"
	::"m"(dst), "m"(constPtr), "m"(src),"m"(count)
    :"esi","edi","eax","edx","ecx");

}

/*
============
CSIMD_SSE::dot

  dst[i] = src0[i] * src1[i];
============
*/
void VPCALL CSIMD_SSE::dot( float *dst, const CVec3D *src0, const CVec3D *src1, const int count ) {

	asm(
		"mov			eax, %3\n"      //count
		"mov			edi, %1\n"     //src
		"mov			edx, eax\n"
		"mov			esi, %2\n"      //src
		"mov			ecx, %0\n"    //dest
		"and			eax, ~3\n"

		"jz			done4%=\n"
		"imul		eax, 12\n"
		"add			edi, eax\n"
		"add			esi, eax\n"
		"neg			eax\n"

	"loop4%=:\n"
		"movlps		xmm0, [esi+eax]\n"					// 0, 1, X, X
		"movlps		xmm3, [edi+eax]\n"						// 0, 1, X, X
		"movlps		xmm1, [esi+eax+8]\n"					// 2, 3, X, X
		"movlps		xmm4, [edi+eax+8]\n"					// 2, 3, X, X
		"movhps		xmm0, [esi+eax+24]\n"					// 0, 1, 6, 7
		"movhps		xmm3, [edi+eax+24]\n"					// 0, 1, 6, 7
		"movhps		xmm1, [esi+eax+32]\n"					// 2, 3, 8, 9
		"movhps		xmm4, [edi+eax+32]\n"					// 2, 3, 8, 9
		"movlps		xmm2, [esi+eax+16]\n"					// 4, 5, X, X
		"movlps		xmm5, [edi+eax+16]\n"					// 4, 5, X, X
		"movhps		xmm2, [esi+eax+40]\n"					// 4, 5, 10, 11
		"movhps		xmm5, [edi+eax+40]\n"					// 4, 5, 10, 11

		"add			ecx, 16\n"
		"add			eax, 48\n"

		"mulps		xmm0, xmm3\n"
		"mulps		xmm1, xmm4\n"
		"mulps		xmm2, xmm5\n"
		"movaps		xmm7, xmm0\n"
		"shufps		xmm7, xmm1, 0xd8 \n" /*R_SHUFFLEPS( 0, 2, 1, 3 )=11011000*/	// 0, 6, 3, 9
		"shufps		xmm0, xmm2, 0x8d\n"  /*R_SHUFFLEPS( 1, 3, 0, 2 )*/	// 1, 7, 4, 10
		"shufps		xmm1, xmm2, 0xd8 \n" /*R_SHUFFLEPS( 0, 2, 1, 3 )=11011000*/	// 2, 8, 5, 11
		"addps		xmm7, xmm0\n"
		"addps		xmm7, xmm1\n"
		"shufps		xmm7, xmm7, 0xd8 \n" /*R_SHUFFLEPS( 0, 2, 1, 3 )=11011000*/

		"movlps		[ecx-16+0], xmm7\n"
		"movhps		[ecx-16+8], xmm7\n"
		"jl			loop4%=\n"

	"done4%=:\n"
		"and			edx, 3\n"
		"jz			done1%=\n"

	"loop1%=:\n"
		"movss		xmm0, [esi+eax+0]\n"
		"movss		xmm3, [edi+eax+0]\n"
		"movss		xmm1, [esi+eax+4]\n"
		"movss		xmm4, [edi+eax+4]\n"
		"movss		xmm2, [esi+eax+8]\n"
		"movss		xmm5, [edi+eax+8]\n"
		"mulss		xmm0, xmm3\n"
		"mulss		xmm1, xmm4\n"
		"mulss		xmm2, xmm5\n"
		"add		ecx, 4\n"
		"addss		xmm0, xmm1\n"
		"add		eax, 12\n"
		"addss		xmm0, xmm2\n"
		"dec		edx\n"
		"movss		[ecx-4], xmm0\n"
		"jnz		loop1%=\n"

	"done1%=:\n"
	::"m"(dst), "m"(src0), "m"(src1), "m"(count)
    :);
}

/*
============
CSIMD_SSE::dot

  dot = src1[0] * src2[0] + src1[1] * src2[1] + src1[2] * src2[2] + ...
============
*/
void VPCALL CSIMD_SSE::dot( float &dot, const float *src1, const float *src2, const int count ) {

float * dotPtr = &dot;

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

			asm (
				"mov		ecx, %1\n"
				"mov		edx, %2\n"
				"mov		eax, ecx\n"
				"or			eax, edx\n"
				"and		eax, 15\n"
				"jz			alignedDot%=\n"
				// unaligned
				"mov			eax, %3\n"
				"shr			eax, 2\n"
				"shl			eax, 4\n"
				"add			ecx, eax\n"
				"add			edx, eax\n"
				"neg			eax\n"
				"movups		xmm0, [ecx+eax]\n"
				"movups		xmm1, [edx+eax]\n"
				"mulps		xmm0, xmm1\n"
				"add			eax, 16\n"
				"jz			doneDot%=\n"
			"loopUnalignedDot%=:\n"
				"movups		xmm1, [ecx+eax]\n"
				"movups		xmm2, [edx+eax]\n"
				"mulps		xmm1, xmm2\n"
				"addps		xmm0, xmm1\n"
				"add		eax, 16\n"
				"jl			loopUnalignedDot%=\n"
				"jmp		doneDot%=\n"
				// aligned
			"alignedDot%=:\n"
				"mov		eax, %3\n"
				"shr		eax, 2\n"
				"shl		eax, 4\n"
				"add		ecx, eax\n"
				"add		edx, eax\n"
				"neg		eax\n"
				"movaps		xmm0, [ecx+eax]\n"
				"movaps		xmm1, [edx+eax]\n"
				"mulps		xmm0, xmm1\n"
				"add		eax, 16\n"
				"jz			doneDot%=\n"
			"loopAlignedDot%=:\n"
				"movaps		xmm1, [ecx+eax]\n"
				"movaps		xmm2, [edx+eax]\n"
				"mulps		xmm1, xmm2\n"
				"addps		xmm0, xmm1\n"
				"add		eax, 16\n"
				"jl			loopAlignedDot%=\n"
			"doneDot%=:\n"
			:
            :"m"(dot), "m"(src1), "m"(src2), "m"(count)
            :"ecx","edx","eax");
			}
			switch( count & 3 ) {
				case 1:
					asm volatile(
						"movss	xmm1, [ecx]\n"
						"movss	xmm2, [edx]\n"
						"mulss	xmm1, xmm2\n"
						"addss	xmm0, xmm1\n"
					);
					break;
				case 2:
					asm volatile(
						"xorps	xmm2, xmm2\n"
						"movlps	xmm1, [ecx]\n"
						"movlps	xmm2, [edx]\n"
						"mulps	xmm1, xmm2\n"
						"addps	xmm0, xmm1\n"
					);
					break;
				case 3:
					asm volatile(
						"movss	xmm1, [ecx]\n"
						"movhps	xmm1, [ecx+4]\n"
						"movss	xmm2, [edx]\n"
						"movhps	xmm2, [edx+4]\n"
						"mulps	xmm1, xmm2\n"
						"addps	xmm0, xmm1\n"
					);
					break;
			}
			asm volatile(
				"movhlps	xmm1, xmm0\n"
				"addps		xmm0, xmm1\n"
				"movaps		xmm1, xmm0\n"
				"shufps		xmm1, xmm1, 0x01 \n"  /*R_SHUFFLEPS( 1, 0, 0, 0 )*/
				"addss		xmm0, xmm1\n"
				"mov		eax, %0\n"
				"movss		[eax], xmm0\n"
			::"m"(dotPtr)
			    :"eax");
			return;
}

//
//	"cmpeqps		==		Equal
//	"cmpneqps	!=		Not Equal
//	"cmpltps		<		Less Than
//  "cmpnltps	>=		Not Less Than
//	"cmpnleps	>		Not Less Or Equal
//
#define FLIP	"not al\n"
#define NOFLIP

#define COMPARECONSTANT( DST, SRC0, CONSTANT, COUNT, CMP, CMPSIMD, DOFLIP )				\
	int i, cnt, pre, post;																\
	float *aligned;																		\
	bool exit=0;                                                            \
	/* if the float array is not aligned on a 4 sf_u8 boundary */					\
	if ( ((int) SRC0) & 3 ) {															\
		/* unaligned memory access */													\
		pre = 0;																		\
		cnt = COUNT >> 2;																\
		post = COUNT - (cnt<<2);														\
		asm(                                                                           \
        "mov		edx, %3\n"													\
		"test		edx, edx\n"													\
		"je			ExitLabel%=\n"														\
		"push		ebx\n"															\
		"neg			edx\n"															\
		"mov			esi, %1\n"													\
		"prefetchnta	[esi+64]\n"													\
		"movss		xmm1, %2\n"												\
		"shufps		xmm1, xmm1, 0x00\n"						\
		"mov		edi, %0\n"													\
		"mov		ecx, 0x01010101\n"												\
    "loopNA%=:\n"																	\
		"movups		xmm0, [esi]	\n"												\
		"prefetchnta	[esi+128]\n"													\
		CMPSIMD"	xmm0, xmm1\n"													\
		"movmskps	eax, xmm0\n"													\
		DOFLIP																	\
		"mov			ah, al\n"														\
		"shr			ah, 1\n"														\
		"mov			bx, ax\n"														\
		"shl			ebx, 14\n"														\
		"mov			bx, ax\n"														\
		"and			ebx, ecx\n"													\
		"mov			dword ptr [edi], ebx\n"										\
		"add			esi, 16\n"														\
		"add			edi, 4\n"														\
		"inc			edx\n"															\
		"jl			loopNA%=\n"														\
		"pop			ebx\n"                                              \
	"ExitLabel%=:"	     \
	     "mov      %4, 0x01\n"                                          \
	    ::"m"(DST), "m"(SRC0),"m" (CONSTANT),"m"(cnt),"m"(exit)                                                                         \
		    :);                                     \
		if(exit==1)  goto doneCmp;    															\
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
			asm(                                            \
            "mov		edx, %1\n"												\
			"test		edx, edx\n"												\
			"je			ExitLabel%=\n"													\
			"push		ebx\n"														\
			"neg			edx\n"														\
			"mov			esi, %3\n"											\
			"prefetchnta	[esi+64]\n"												\
			"movss		xmm1, %0\n"											\
			"shufps		xmm1, xmm1, 0x00\n"					\
			"mov			edi, %4\n"												\
			"add			edi, %2\n"												\
			"mov			ecx, 0x01010101\n"											\
			"loopA%=:\n"																\
			"movaps		xmm0, [esi]\n"												\
			"prefetchnta	[esi+128]\n"												\
			CMPSIMD"		xmm0, xmm1\n"												\
			"movmskps	eax, xmm0\n"												\
			DOFLIP																\
			"mov			ah, al\n"													\
			"shr			ah, 1\n"													\
			"mov			bx, ax\n"													\
			"shl			ebx, 14\n"													\
			"mov			bx, ax\n"													\
			"and			ebx, ecx\n"												\
			"mov			dword ptr [edi], ebx\n"									\
			"add			esi, 16\n"													\
			"add			edi, 4\n"													\
			"inc			edx\n"														\
			"jl			loopA%=\n"													\
			"pop			ebx\n"	\
		"ExitLabel%=:"	     \
            "mov      %5, 0x01\n"                                          \
	    	::"m" (CONSTANT),"m"(cnt),"m"(pre),"m"(aligned),"m"(DST),"m"(exit)                                                                         \
		    :"esi","ebx","edi","edx");                                     \
		if(exit==1)  goto doneCmp;                                              \
                                                                                        \
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
	bool exit;																					\
	/* if the float array is not aligned on a 4 sf_u8 boundary */						\
	if ( ((int) SRC0) & 3 ) {															\
		/* unaligned memory access */													\
		pre = 0;																		\
		cnt = COUNT >> 2;																\
		post = COUNT - (cnt<<2);														\
		asm("mov	edx, %1\n"													\
		"test		edx, edx\n"													\
		"je			ExitLabel%=\n"													\
		"push		ebx	\n"														\
		"neg			edx\n"															\
		"mov			esi, %3\n"													\
		"prefetchnta	[esi+64]\n"													\
		"movss		xmm1, %0\n"												\
		"shufps		xmm1, xmm1, 0x00\n"						\
		"mov			edi, %2\n"													\
		"mov			cl, %4\n"													\
		"loopNA%=:\n"																	\
		"movups		xmm0, [esi]\n"													\
		"prefetchnta	[esi+128]\n"													\
		CMPSIMD"		xmm0, xmm1\n"													\
		"movmskps	eax, xmm0\n"													\
		DOFLIP																	\
		"mov			ah, al\n"														\
		"shr			ah, 1\n"														\
		"mov			bx, ax\n"														\
		"shl			ebx, 14\n"														\
		"mov			bx, ax\n"														\
		"and			ebx, 0x01010101\n"												\
		"shl			ebx, cl\n"														\
		"or			ebx, dword ptr [edi]\n"										\
		"mov			dword ptr [edi], ebx\n"										\
		"add			esi, 16\n"														\
		"add			edi, 4\n"														\
		"inc			edx	\n"														\
		"jl			loopNA%=	\n"													\
		"pop			ebx\n"                                  \
	"ExitLabel%=:"	     \
        "mov      %5, 0x01\n"                                          \
        ::"m"(CONSTANT),"m"(cnt),"m"(DST),"m"(SRC0), "m"(BITNUM), "m"(exit)          \
        :"esi","ebx","edi","edx");	                                                                    \
        if(exit==1)  goto doneCmp;                                              \
                               														\
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
			asm("mov			edx, %1\n"												\
			"test		edx, edx\n"												\
			"je			ExitLabel%=\n"													\
			"push		ebx\n"														\
			"neg			edx\n"														\
			"mov			esi, %3\n"											\
			"prefetchnta	[esi+64]\n"												\
			"movss		xmm1, %0\n"											\
			"shufps		xmm1, xmm1, 0x00\n"					\
			"mov			edi, %4\n"												\
			"add			edi, %2\n"												\
			"mov			cl, %5\n"												\
			"loopA%=:\n"																\
			"movaps		xmm0, [esi]\n"												\
			"prefetchnta	[esi+128]\n"												\
			CMPSIMD"		xmm0, xmm1\n"												\
			"movmskps	eax, xmm0\n"												\
			DOFLIP																\
			"mov			ah, al\n"													\
			"shr			ah, 1\n"													\
			"mov			bx, ax\n"													\
			"shl			ebx, 14\n"													\
			"mov			bx, ax\n"													\
			"and			ebx, 0x01010101\n"											\
			"shl			ebx, cl\n"													\
			"or			ebx, dword ptr [edi]\n"									\
			"mov			dword ptr [edi], ebx\n"									\
			"add			esi, 16\n"													\
			"add			edi, 4\n"													\
			"inc			edx\n"														\
			"jl			loopA%=\n"													\
			"pop			ebx\n"                     \
	"ExitLabel%=:"	     \
            "mov      %6, 0x01\n"                                          \
			 ::"m"(CONSTANT),"m"(cnt),"m"(pre),"m"(aligned),"m"(DST), "m"(BITNUM),"m"(exit)                                                                                                                                        \
        :"esi","ebx","edi","edx");															\
		if(exit==1)  goto doneCmp;                                              \
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
CSIMD_SSE::cmp

  dst[i] = src0[i] > constant;
============
*/
void VPCALL CSIMD_SSE::cmpGT( sf_u8 *dst, const float *src0, const float constant, const int count ) {
	 COMPARECONSTANT( dst, src0, constant, count, >, "cmpnleps", NOFLIP )
}

/*
============
CSIMD_SSE::cmp

  dst[i] |= ( src0[i] > constant ) << bitNum;
============
*/
void VPCALL CSIMD_SSE::cmpGT( sf_u8 *dst, const sf_u8 bitNum, const float *src0, const float constant, const int count ) {
	 COMPAREBITCONSTANT( dst, bitNum, src0, constant, count, >, "cmpnleps", NOFLIP )
}

/*
============
CSIMD_SSE::"cmpGE

  dst[i] = src0[i] >= constant;
============
*/
void VPCALL CSIMD_SSE::cmpGE( sf_u8 *dst, const float *src0, const float constant, const int count ) {
	COMPARECONSTANT( dst, src0, constant, count, >=, "cmpnltps", NOFLIP )
}

/*
============
CSIMD_SSE::cmpGE

  dst[i] |= ( src0[i] >= constant ) << bitNum;
============
*/
void VPCALL CSIMD_SSE::cmpGE( sf_u8 *dst, const sf_u8 bitNum, const float *src0, const float constant, const int count ) {
	COMPAREBITCONSTANT( dst, bitNum, src0, constant, count, >=, "cmpnltps", NOFLIP )
}

/*
============
CSIMD_SSE::cmpLT

  dst[i] = src0[i] < constant;
============
*/
void VPCALL CSIMD_SSE::cmpLT( sf_u8 *dst, const float *src0, const float constant, const int count ) {
	COMPARECONSTANT( dst, src0, constant, count, <, "cmpltps", NOFLIP )
}

/*
============
CSIMD_SSE::cmpLT

  dst[i] |= ( src0[i] < constant ) << bitNum;
============
*/
void VPCALL CSIMD_SSE::cmpLT( sf_u8 *dst, const sf_u8 bitNum, const float *src0, const float constant, const int count ) {
	COMPAREBITCONSTANT( dst, bitNum, src0, constant, count, <, "cmpltps", NOFLIP )
}

/*
============
CSIMD_SSE::cmpLE

  dst[i] = src0[i] <= constant;
============
*/
void VPCALL CSIMD_SSE::cmpLE( sf_u8 *dst, const float *src0, const float constant, const int count ) {
	COMPARECONSTANT( dst, src0, constant, count, <=, "cmpnleps", FLIP )
}

/*
============
CSIMD_SSE::cmpLE

  dst[i] |= ( src0[i] <= constant ) << bitNum;
============
*/
void VPCALL CSIMD_SSE::cmpLE( sf_u8 *dst, const sf_u8 bitNum, const float *src0, const float constant, const int count ) {
	COMPAREBITCONSTANT( dst, bitNum, src0, constant, count, <=, "cmpnleps", FLIP )
}

/*
============
CSIMD_SSE::minMax
============
*/
void VPCALL CSIMD_SSE::minMax( float &min, float &max, const float *src, const int count ) {

	int i, pre, post;

	min = CMath::INFINITY_FLOAT; max = -CMath::INFINITY_FLOAT;

asm(
		"push		ebx\n"
		"mov		eax, %0\n"
		"mov		ebx, %1\n"
		"movss		xmm0, [eax]\n"
		"movss		xmm1, [ebx]\n"
		"shufps		xmm0, xmm0, 0\n"
		"shufps		xmm1, xmm1, 0\n"

		KFLOATINITS( %2,%3,%4,%5 )
		"and			eax, 15\n"
		"jz			lpA%=\n"
		"jmp			lpNA%=\n"
		/*align		16*/
"lpNA%=:\n"
		"movups		xmm2, [edx+ebx]\n"
		"movups		xmm3, [edx+ebx+16]\n"
		"minps		xmm0, xmm2\n"
		"maxps		xmm1, xmm2\n"
		"prefetchnta	[edx+ebx+64]\n"
		"minps		xmm0, xmm3\n"
		"maxps		xmm1, xmm3\n"
		"add			ebx, 16*2\n"
		"jl			lpNA%=\n"
		"jmp			done2%=\n"
"lpA%=:\n"
		"movaps		xmm2, [edx+ebx]\n"
		"movaps		xmm3, [edx+ebx+16]\n"
		"minps		xmm0, xmm2\n"
		"maxps		xmm1, xmm2\n"
		"prefetchnta	[edx+ebx+64]\n"
		"minps		xmm0, xmm3\n"
		"maxps		xmm1, xmm3\n"
		"add		ebx, 16*2\n"
		"jl			lpA%=\n"
		"jmp		done2%=\n"
		/*align		16*/
"done2%=:\n"
		"movaps		xmm2, xmm0\n"
		"movaps		xmm3, xmm1\n"
		"shufps		xmm2, xmm2, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm3, xmm3, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"minss		xmm0, xmm2\n"
		"maxss		xmm1, xmm3\n"
		"shufps		xmm2, xmm2, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm3, xmm3, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"minss		xmm0, xmm2\n"
		"maxss		xmm1, xmm3\n"
		"shufps		xmm2, xmm2, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm3, xmm3, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"minss		xmm0, xmm2\n"
		"maxss		xmm1, xmm3\n"
		"mov		eax, %0\n"
		"mov		ebx, %1\n"
		"movss		[eax], xmm0\n"
		"movss		[ebx], xmm1\n"
"done%=:\n"
		"pop			ebx\n"
	:
    :"m"(&min), "m"(&max),"m"(src),"m"(count), "m"(pre), "m"(post)
	:);

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

	asm (
		"mov		eax, %3\n"
		"test		eax, eax\n"
		"movss		xmm0, %4\n"
		"xorps		xmm1, xmm1\n"
		"shufps		xmm0, xmm0, 0x00\n"
		"subps		xmm1, xmm0\n"
		"jz			done%=\n"
		"mov		ecx, eax\n"
		"and		ecx, 1\n"
		"mov		esi, %2\n"
		"jz			startLoop\n"
		"movlps		xmm2, [esi]\n"
		"shufps		xmm2, xmm2, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
		"dec		eax\n"
		"add		esi, 2*4\n"
		"minps		xmm0, xmm2\n"
		"maxps		xmm1, xmm2\n"
	"startLoop:\n"
		"imul		eax, 2*4\n"
		"add		esi, eax\n"
		"neg		eax\n"
	"loopVert:\n"
		"movlps		xmm2, [esi+eax]\n"
		"movhps		xmm2, [esi+eax+8]\n"
		"add		eax, 4*4\n"
		"minps		xmm0, xmm2\n"
		"maxps		xmm1, xmm2\n"
		"jl			loopVert\n"
	"done%=:\n"
		"movaps		xmm2, xmm0\n"
		"shufps		xmm2, xmm2, 0x4E\n"   /*R_SHUFFLEPS( 2, 3, 0, 1 )=01001110*/
		"minps		xmm0, xmm2\n"
		"mov		esi, %0\n"
		"movlps		[esi], xmm0\n"
		"movaps		xmm3, xmm1\n"
		"shufps		xmm3, xmm3, 0x4E\n"   /*R_SHUFFLEPS( 2, 3, 0, 1 )=01001110*/
		"maxps		xmm1, xmm3\n"
		"mov		edi, %1\n"
		"movlps		[edi], xmm1\n"
	:
    :"m"(&min),"m"(&max),"m"(src), "m"(count), "m"(MATH::CMath::INFINITY_FLOAT)
    :);

}

/*
============
CSIMD_SSE::minMax
============
*/
void VPCALL CSIMD_SSE::minMax( CVec3D &min, CVec3D &max, const CVec3D *src, const int count ) {

	asm (

		"movss		xmm0, %4\n"
		"xorps		xmm1, xmm1\n"
		"shufps		xmm0, xmm0, 0x00\n"
		"subps		xmm1, xmm0\n"
		"movaps		xmm2, xmm0\n"
		"movaps		xmm3, xmm1\n"

		"mov		esi, %2\n"
		"mov		eax, %3\n"
		"and		eax, ~3\n"
		"jz			done4%=\n"
		"imul		eax, 12\n"
		"add		esi, eax\n"
		"neg		eax\n"

	"loop4:\n"
//		"prefetchnta	[esi+4*12]

		"movss		xmm4, [esi+eax+0*12+8]\n"
		"movhps		xmm4, [esi+eax+0*12+0]\n"
		"minps		xmm0, xmm4\n"
		"maxps		xmm1, xmm4\n"

		"movss		xmm5, [esi+eax+1*12+0]\n"
		"movhps		xmm5, [esi+eax+1*12+4]\n"
		"minps		xmm2, xmm5\n"
		"maxps		xmm3, xmm5\n"

		"movss		xmm6, [esi+eax+2*12+8]\n"
		"movhps		xmm6, [esi+eax+2*12+0]\n"
		"minps		xmm0, xmm6\n"
		"maxps		xmm1, xmm6\n"

		"movss		xmm7, [esi+eax+3*12+0]\n"
		"movhps		xmm7, [esi+eax+3*12+4]\n"
		"minps		xmm2, xmm7\n"
		"maxps		xmm3, xmm7\n"

		"add		eax, 4*12\n"
		"jl			loop4\n"

	"done4%=:\n"
		"mov		eax, %3\n"
		"and		eax, 3\n"
		"jz			done1%=\n"
		"imul		eax, 12\n"
		"add		esi, eax\n"
		"neg		eax\n"

	"loop1%=:\n"
		"movss		xmm4, [esi+eax+0*12+8]\n"
		"movhps		xmm4, [esi+eax+0*12+0]\n"
		"minps		xmm0, xmm4\n"
		"maxps		xmm1, xmm4\n"

		"add		eax, 12\n"
		"jl			loop1%=\n"

	"done1%=:\n"
		"shufps		xmm2, xmm2, 0x87\n"  /*R_SHUFFLEPS( 3, 1, 0, 2 )=10000111*/
		"shufps		xmm3, xmm3, 0x87\n"  /*R_SHUFFLEPS( 3, 1, 0, 2 )=10000111*/
		"minps		xmm0, xmm2\n"
		"maxps		xmm1, xmm3\n"
		"mov		esi, %0\n"
		"movhps		[esi], xmm0\n"
		"movss		[esi+8], xmm0\n"
		"mov		edi, %1\n"
		"movhps		[edi], xmm1\n"
		"movss		[edi+8], xmm1\n"
    :
    :"m"(&min),"m"(&max),"m"(src), "m"(count), "m"(MATH::CMath::INFINITY_FLOAT)
    :);
}

/*
============
CSIMD_SSE::minMax
============
*/
void VPCALL CSIMD_SSE::minMax( CVec3D &min, CVec3D &max, const CVertex *src, const int count ) {
	SMF_ASSERT( sizeof( CVertex ) == DRAWVERT_SIZE );
	SMF_ASSERT( (int)&((CVertex *)0)->xyz == DRAWVERT_XYZ_OFFSET );

	asm (

		"movss		xmm0, %4\n"
		"xorps		xmm1, xmm1\n"
		"shufps		xmm0, xmm0, 0x00\n"
		"subps		xmm1, xmm0\n"
		"movaps		xmm2, xmm0\n"
		"movaps		xmm3, xmm1\n"

		"mov		esi, %2\n"
		"mov		eax, %3\n"
		"and		eax, ~3\n"
		"jz			done4%=\n"
		"imul		eax, "DRAWVERT_SIZE_STR"\n"
		"add		esi, eax\n"
		"neg		eax\n"

	"loop4%=:\n"
//		"prefetchnta	[esi+4*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"]

		"movss		xmm4, [esi+eax+0*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"movhps		xmm4, [esi+eax+0*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"minps		xmm0, xmm4\n"
		"maxps		xmm1, xmm4\n"

		"movss		xmm5, [esi+eax+1*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"movhps		xmm5, [esi+eax+1*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"minps		xmm2, xmm5\n"
		"maxps		xmm3, xmm5\n"

		"movss		xmm6, [esi+eax+2*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"movhps		xmm6, [esi+eax+2*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"minps		xmm0, xmm6\n"
		"maxps		xmm1, xmm6\n"

		"movss		xmm7, [esi+eax+3*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"movhps		xmm7, [esi+eax+3*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"minps		xmm2, xmm7\n"
		"maxps		xmm3, xmm7\n"

		"add		eax, 4*"DRAWVERT_SIZE_STR"\n"
		"jl			loop4%=\n"

	"done4%=:\n"
		"mov		eax, %3\n"
		"and		eax, 3\n"
		"jz			done1%=\n"
		"imul		eax, "DRAWVERT_SIZE_STR"\n"
		"add		esi, eax\n"
		"neg		eax\n"

	"loop1%=:\n"
		"movss		xmm4, [esi+eax+0*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"movhps		xmm4, [esi+eax+0*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"minps		xmm0, xmm4\n"
		"maxps		xmm1, xmm4\n"

		"add		eax, "DRAWVERT_SIZE_STR"\n"
		"jl			loop1%=\n"

	"done1%=:\n"
		"shufps		xmm2, xmm2, 0x87\n"  /*R_SHUFFLEPS( 3, 1, 0, 2 )=10000111*/
		"shufps		xmm3, xmm3, 0x87\n"  /*R_SHUFFLEPS( 3, 1, 0, 2 )=10000111*/
		"minps		xmm0, xmm2\n"
		"maxps		xmm1, xmm3\n"
		"mov		esi, %0\n"
		"movhps		[esi], xmm0\n"
		"movss		[esi+8], xmm0\n"
		"mov		edi, %1\n"
		"movhps		[edi], xmm1\n"
		"movss		[edi+8], xmm1\n"
    :
    :"m"(&min),"m"(&max),"m"(src), "m"(count), "m"(MATH::CMath::INFINITY_FLOAT)
    :);
}

/*
============
CSIMD_SSE::minMax
============
*/
void VPCALL CSIMD_SSE::minMax( CVec3D &min, CVec3D &max, const CVertex *src, const int *indexes, const int count ) {

	SMF_ASSERT( sizeof( CVertex ) == DRAWVERT_SIZE );
	SMF_ASSERT( (int)&((CVertex *)0)->xyz == DRAWVERT_XYZ_OFFSET );

	asm (

		"movss		xmm0, %4\n"
		"xorps		xmm1, xmm1\n"
		"shufps		xmm0, xmm0, 0x00\n"
		"subps		xmm1, xmm0\n"
		"movaps		xmm2, xmm0\n"
		"movaps		xmm3, xmm1\n"

		"mov			edi, %5\n"
		"mov			esi, %2\n"
		"mov			eax, %3\n"
		"and			eax, ~3\n"
		"jz			done4%=\n"
		"shl			eax, 2\n"
		"add			edi, eax\n"
		"neg			eax\n"

	"loop4%=:\n"
//		"prefetchnta	[edi+128]
//		"prefetchnta	[esi+4*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"]

		"mov		edx, [edi+eax+0]\n"
		"imul		edx, "DRAWVERT_SIZE_STR"\n"
		"movss		xmm4, [esi+edx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"movhps		xmm4, [esi+edx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"minps		xmm0, xmm4\n"
		"maxps		xmm1, xmm4\n"

		"mov		edx, [edi+eax+4]\n"
		"imul		edx, "DRAWVERT_SIZE_STR"\n"
		"movss		xmm5, [esi+edx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"movhps		xmm5, [esi+edx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"minps		xmm2, xmm5\n"
		"maxps		xmm3, xmm5\n"

		"mov		edx, [edi+eax+8]\n"
		"imul		edx, "DRAWVERT_SIZE_STR"\n"
		"movss		xmm6, [esi+edx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"movhps		xmm6, [esi+edx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"minps		xmm0, xmm6\n"
		"maxps		xmm1, xmm6\n"

		"mov		edx, [edi+eax+12]\n"
		"imul		edx, "DRAWVERT_SIZE_STR"\n"
		"movss		xmm7, [esi+edx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"movhps		xmm7, [esi+edx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"minps		xmm2, xmm7\n"
		"maxps		xmm3, xmm7\n"

		"add		eax, 4*4\n"
		"jl			loop4%=\n"

	"done4%=:\n"
		"mov		eax, %3\n"
		"and		eax, 3\n"
		"jz			done1%=\n"
		"shl		eax, 2\n"
		"add		edi, eax\n"
		"neg		eax\n"

	"loop1%=:\n"
		"mov		edx, [edi+eax+0]\n"
		"imul		edx, "DRAWVERT_SIZE_STR"\n"
		"movss		xmm4, [esi+edx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"movhps		xmm4, [esi+edx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"minps		xmm0, xmm4\n"
		"maxps		xmm1, xmm4\n"

		"add		eax, 4\n"
		"jl			loop1%=\n"

	"done1%=:\n"
		"shufps		xmm2, xmm2, 0x87\n"  /*R_SHUFFLEPS( 3, 1, 0, 2 )=10000111*/
		"shufps		xmm3, xmm3, 0x87\n"  /*R_SHUFFLEPS( 3, 1, 0, 2 )=10000111*/
		"minps		xmm0, xmm2\n"
		"maxps		xmm1, xmm3\n"
		"mov		esi, %0\n"
		"movhps		[esi], xmm0\n"
		"movss		[esi+8], xmm0\n"
		"mov		edi, %1\n"
		"movhps		[edi], xmm1\n"
		"movss		[edi+8], xmm1\n"
	:
    :"m"(&min),"m"(&max),"m"(src), "m"(count), "m"(MATH::CMath::INFINITY_FLOAT),"m"(indexes)
    :);
}

/*
============
CSIMD_SSE::clamp
============
*/
void VPCALL CSIMD_SSE::clamp( float *dst, const float *src, const float min, const float max, const int count ) {

	int	i, pre, post;

	asm(
		"movss	xmm0,%2\n"
		"movss	xmm1,%3\n"
		"shufps	xmm0,xmm0,0\n"
		"shufps	xmm1,xmm1,0\n"

		KFLOATINITDS( %0, %1, %4, %5, %6 )
		"and		eax,15\n"
		"jne		lpNA%=\n"
		"jmp		lpA%=\n"
		/*align	16*/
"lpA%=:\n"
		"movaps	xmm2,[edx+ebx]\n"
		"movaps	xmm3,[edx+ebx+16]\n"
		"maxps	xmm2,xmm0\n"
		"maxps	xmm3,xmm0\n"
		"prefetchnta	[edx+ebx+64]\n"
		"minps	xmm2,xmm1\n"
		"minps	xmm3,xmm1\n"
		"movaps	[edi+ebx],xmm2\n"
		"movaps	[edi+ebx+16],xmm3\n"
		"add	ebx,16*2\n"
		"jl		lpA%=\n"
		"jmp	done%=\n"

		/*align	16*/
"lpNA%=:\n"
		"movups	xmm2,[edx+ebx]\n"
		"movups	xmm3,[edx+ebx+16]\n"
		"maxps	xmm2,xmm0\n"
		"maxps	xmm3,xmm0\n"
		"prefetchnta	[edx+ebx+64]\n"
		"minps	xmm2,xmm1\n"
		"minps	xmm3,xmm1\n"
		"movaps	[edi+ebx],xmm2\n"
		"movaps	[edi+ebx+16],xmm3\n"
		"add	ebx,16*2\n"
		"jl		lpNA%=\n"
"done%=:\n"
    :
    :"m"(dst), "m"(src), "m"(min), "m"(max), "m"(count),"m"(pre), "m"(post)
    :);


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

	asm(
		"movss	xmm0,%2\n"
		"shufps	xmm0,xmm0,0\n"

		KFLOATINITDS( %0, %1, %3, %4, %5 )
		"and		eax,15\n"
		"jne		lpNA%=\n"
		"jmp		lpA%=\n"
		/*align	16*/
"lpA%=:\n"
		"movaps	xmm2,[edx+ebx]\n"
		"movaps	xmm3,[edx+ebx+16]\n"
		"maxps	xmm2,xmm0\n"
		"prefetchnta	[edx+ebx+64]\n"
		"maxps	xmm3,xmm0\n"
		"movaps	[edi+ebx],xmm2\n"
		"movaps	[edi+ebx+16],xmm3\n"
		"add		ebx,16*2\n"
		"jl		lpA%=\n"
		"jmp		done%=\n"

		/*align	16*/
"lpNA%=:\n"
		"movups	xmm2,[edx+ebx]\n"
		"movups	xmm3,[edx+ebx+16]\n"
		"maxps	xmm2,xmm0\n"
		"prefetchnta	[edx+ebx+64]\n"
		"maxps	xmm3,xmm0\n"
		"movaps	[edi+ebx],xmm2\n"
		"movaps	[edi+ebx+16],xmm3\n"
		"add		ebx,16*2\n"
		"jl		lpNA%=\n"
"done%=:\n"
    :
    :"m"(dst), "m"(src), "m"(min),  "m"(count),"m"(pre), "m"(post)
    :);

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

	asm (
		"movss	xmm1,%5\n"
		"shufps	xmm1,xmm1,0\n"
        KFLOATINITDS( %0,%1,%2,%3,%4 )
		"and		eax,15\n"
		"jne		lpNA%=\n"
		"jmp		lpA%=\n"
		/*align	16*/
"lpA%=:\n"
		"movaps	xmm2,[edx+ebx]\n"
		"movaps	xmm3,[edx+ebx+16]\n"
		"minps	xmm2,xmm1\n"
		"prefetchnta	[edx+ebx+64]\n"
		"minps	xmm3,xmm1\n"
		"movaps	[edi+ebx],xmm2\n"
		"movaps	[edi+ebx+16],xmm3\n"
		"add		ebx,16*2\n"
		"jl		lpA%=\n"
		"jmp		done%=\n"

		/*align	16*/
"lpNA%=:\n"
		"movups	xmm2,[edx+ebx]\n"
		"movups	xmm3,[edx+ebx+16]\n"
		"minps	xmm2,xmm1\n"
		"prefetchnta	[edx+ebx+64]\n"
		"minps	xmm3,xmm1\n"
		"movaps	[edi+ebx],xmm2\n"
		"movaps	[edi+ebx+16],xmm3\n"
		"add		ebx,16*2\n"
		"jl		lpNA%=\n"
"done%=:\n"
	:
	:"m"(dst), "m"(src),"m"(count), "m"(pre), "m"(post),"m"(max)
	:);

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

	asm (
		"mov		edx, %0\n"
		"mov		eax, %1\n"
		"add		eax, 3\n"
		"shr		eax, 2\n"
		"jz		doneZero16\n"
		"shl		eax, 4\n"
		"add		edx, eax\n"
		"neg		eax\n"
		"xorps	xmm0, xmm0\n"
	"loopZero16:\n"
		"movaps	[edx+eax], xmm0\n"
		"add		eax, 16\n"
		"jl		loopZero16\n"
	"doneZero16:\n"
	::"m"(dst),"m"(count)
	    :);

}

/*
============
CSIMD_SSE::negate16
============
*/
void VPCALL CSIMD_SSE::negate16( float *dst, const int count ) {

	asm (
		"mov		edx, %2\n"
		"mov		eax, %1\n"
		"add		eax, 3\n"
		"shr		eax, 2\n"
		"jz		donenegate16\n"
		"shl		eax, 4\n"
		"add		edx, eax\n"
		"neg		eax\n"
		"movss	xmm0, %0\n"
		"shufps	xmm0, xmm0, 0x00\n"
	"loopnegate16:\n"
		"movaps	xmm1, [edx+eax]\n"
		"xorps	xmm1, xmm0\n"
		"movaps	[edx+eax], xmm1\n"
		"add		eax, 16\n"
		"jl		loopnegate16\n"
	"donenegate16:\n"
	::"m"(SIMD_SP_signBitMask),"m"(count),"r"(dst)
	    :);

}

/*
============
CSIMD_SSE::copy16
============
*/
void VPCALL CSIMD_SSE::copy16( float *dst, const float *src, const int count ) {

	asm (
		"mov		ecx, %1\n"
		"mov		edx, %0\n"
		"mov		eax, %2\n"
		"add		eax, 3\n"
		"shr		eax, 2\n"
		"jz		donecopy16\n"
		"shl		eax, 4\n"
		"add		ecx, eax\n"
		"add		edx, eax\n"
		"neg		eax\n"
	"loopcopy16:\n"
		"movaps	xmm0, [ecx+eax]\n"
		"movaps	[edx+eax], xmm0\n"
		"add		eax, 16\n"
		"jl		loopcopy16\n"
	"donecopy16:\n"
	::"r"(dst), "r"(src), "m"(count)
    :);

/*
board(,%ecx,1)
section:disp(base,index,scale)
where base "and index are the optional 32-bit base "and index registers, disp is the optional displacement, "and scale
, taking the values 1, 2, 4, "and 8, multiplies index to calculate the address of the oper"and. If no scale is specified, scale is taken to be 1
+------------------------------+------------------------------------+
|       Intel Code             |      AT&T Code                     |
+------------------------------+------------------------------------+
| mov     eax,1                |  movl    $1,%eax                   |
| mov     ebx,0ffh             |  movl    $0xff,%ebx                |
| int     80h                  |  int     $0x80                     |
| mov     ebx, eax             |  movl    %eax, %ebx                |
| mov     eax,[ecx]            |  movl    (%ecx),%eax               |
| mov     eax,[ebx+3]          |  movl    3(%ebx),%eax              |
| mov     eax,[ebx+20h]        |  movl    0x20(%ebx),%eax           |
| add     eax,[ebx+ecx*2h]     |  addl    (%ebx,%ecx,0x2),%eax      |
| "lea     eax,[ebx+ecx]        |  "leal    (%ebx,%ecx),%eax          |
| "sub     eax,[ebx+ecx*4h-20h] |  "subl    -0x20(%ebx,%ecx,0x4),%eax |
+------------------------------+------------------------------------+
*/
#if 0
   __asm__("\
	            xorl %eax, %eax\n\
	            xorl %ebx, %ebx\n\
	            xorl %ecx, %ecx\n\
	            xorl %edx, %edx\n\
	        Init:\n\
	            movl src, %ecx\n\
	            movl dst, %edx\n\
	            movl count, %eax\n\
	            add $3, %eax\n\
	            "shr	$2, %eax\n\
                "jz	donecopy16\n\
                "shl	$4, %eax\n\
                add	%eax, %ecx\n\
                add	%eax, %edx\n\
                "neg	%eax\n\
            loopcopy16:\n\
                "movaps	 (%ecx,%eax,1), %xmm0\n\
                "movaps	 %xmm0,(%edx,%eax,1)\n\
                add		 $16, %eax\n\
                "jl		loopcopy16\n\
	        donecopy16:\n\
            ");
            #endif
}

/*
============
CSIMD_SSE::add16
============
*/


void VPCALL CSIMD_SSE::add16( float *dst, const float *src1, const float *src2, const int count ) {
	asm (
		"mov	ecx, %1\n"
		"mov	edx, %2\n"
		"mov	esi, %0\n"
		"mov	eax, %3\n"
		"add	eax, 3\n"
		"shr	eax, 2\n"
		"jz		doneadd16\n"
		"shl	eax, 4\n"
		"add	esi, eax\n"
		"add	ecx, eax\n"
		"add	edx, eax\n"
		"neg	eax\n"
	"loopadd16:\n"
		"movaps	xmm0, [ecx+eax]\n"
		"addps	xmm0, [edx+eax]\n"
		"movaps	[esi+eax], xmm0\n"
		"add	eax, 16\n"
		"jl		loopadd16\n"
	"doneadd16:\n"
	:: "r"(dst), "r"(src1), "r"(src2), "m"(count)
    :);
}

/*
============
CSIMD_SSE::"sub16
============
*/

void VPCALL CSIMD_SSE::sub16( float *dst, const float *src1, const float *src2, const int count ) {
 asm (
		"mov		ecx, %1\n"
		"mov		edx, %2\n"
		"mov		esi, %0\n"
		"mov		eax, %3\n"
		"add		eax, 3\n"
		"shr		eax, 2\n"
		"jz		donesub16\n"
		"shl		eax, 4\n"
		"add		esi, eax\n"
		"add		ecx, eax\n"
		"add		edx, eax\n"
		"neg		eax\n"
	"loopsub16:\n"
		"movaps	xmm0, [ecx+eax]\n"
		"subps	xmm0, [edx+eax]\n"
		"movaps	[esi+eax], xmm0\n"
		"add		eax, 16\n"
		"jl		loopsub16\n"
	"donesub16:\n"
	:: "r"(dst), "r"(src1), "r"(src2), "m"(count)
    :);

}





/*
============
CSIMD_SSE::mul16
============
*/
void VPCALL CSIMD_SSE::mul16( float *dst, const float *src1, const float constant, const int count ) {
	asm (
		"mov		ecx, %0\n"
		"mov		edx, %1\n"
		"mov		eax, %3\n"
		"add		eax, 3\n"
		"shr		eax, 2\n"
		"jz		doneMulScalar16\n"
		"movss	xmm1, %2\n"
		"shl		eax, 4\n"
		"add		ecx, eax\n"
		"add		edx, eax\n"
		"neg		eax\n"
		"shufps	xmm1, xmm1, 0x00\n"
	"loopMulScalar16:\n"
		"movaps	xmm0, [edx+eax]\n"
		"mulps	xmm0, xmm1\n"
		"movaps	[ecx+eax], xmm0\n"
		"add		eax, 16\n"
		"jl		loopMulScalar16\n"
	"doneMulScalar16:\n"
	::"r"(dst), "r"(src1), "m"(constant), "m"(count)
    :);
}

/*
============
CSIMD_SSE::addAssign16
============
*/
void VPCALL CSIMD_SSE::addAssign16( float *dst, const float *src, const int count ) {
	asm (
		"mov		ecx, %0\n"
		"mov		edx, %1\n"
		"mov		eax, %2\n"
		"add		eax, 3\n"
		"shr		eax, 2\n"
		"jz		doneAddAssign16\n"
		"shl		eax, 4\n"
		"add		ecx, eax\n"
		"add		edx, eax\n"
		"neg		eax\n"
	"loopAddAssign16:\n"
		"movaps	xmm0, [ecx+eax]\n"
		"addps	xmm0, [edx+eax]\n"
		"movaps	[ecx+eax], xmm0\n"
		"add	eax, 16\n"
		"jl		loopAddAssign16\n"
	"doneAddAssign16:\n"
	::"r"(dst), "r"(src), "m"(count)
    :"eax","ecx","edx"
	);
}

/*
============
CSIMD_SSE::"subAssign16
============
*/
void VPCALL CSIMD_SSE::subAssign16( float *dst, const float *src, const int count ) {
	asm (
		"mov		ecx, %0\n"
		"mov		edx, %1\n"
		"mov		eax, %2\n"
		"add		eax, 3\n"
		"shr		eax, 2\n"
		"jz		donesubAssign16\n"
		"shl		eax, 4\n"
		"add		ecx, eax\n"
		"add		edx, eax\n"
		"neg		eax\n"
	"loopsubAssign16:\n"
		"movaps	xmm0, [ecx+eax]\n"
		"subps	xmm0, [edx+eax]\n"
		"movaps	[ecx+eax], xmm0\n"
		"add		eax, 16\n"
		"jl		loopsubAssign16\n"
	"donesubAssign16:\n"
	::"r"(dst), "r"(src), "m"(count)
    :"eax","ecx","edx"
	);
}

/*
============
CSIMD_SSE::mulAssign16
============
*/
void VPCALL CSIMD_SSE::mulAssign16( float *dst, const float constant, const int count ) {
	asm (
		"mov		ecx, %0\n"
		"mov		eax, %2\n"
		"add		eax, 3\n"
		"shr		eax, 2\n"
		"jz		donemulAssign16\n"
		"movss	xmm1, %1\n"
		"shl		eax, 4\n"
		"add		ecx, eax\n"
		"neg		eax\n"
		"shufps	xmm1, xmm1, 0x00\n"
	"loopmulAssign16:\n"
		"movaps	xmm0, [ecx+eax]\n"
		"mulps	xmm0, xmm1\n"
		"movaps	[ecx+eax], xmm0\n"
		"add	eax, 16\n"
		"jl		loopmulAssign16\n"
	"donemulAssign16:\n"
	::"r"(dst), "m"(constant), "m"(count)
    :"eax","ecx","edx"
	);
}

bool VPCALL CSIMD_SSE::matX_inverse_4x4(float* src)
{
__m128 minor0, minor1, minor2, minor3;
__m128 row0, row1, row2, row3;
__m128 det, tmp1;
#if 0
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
#endif
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
	"movss		[eax+offset], reg1\n"

#define STORE2LO( offset, reg1, reg2 )		\
	"movlps	[eax+" STRINGIZE(offset)"], "STRINGIZE(reg1)"\n"

#define STORE2HI( offset, reg1, reg2 )		\
	"movhps	[eax+" STRINGIZE(offset)"], "STRINGIZE(reg1)"\n"

#define STORE4( offset, reg1, reg2 )		\
	"movlps	[eax+" STRINGIZE(offset)"], "STRINGIZE(reg1)"\n"	\
	"movhps	[eax+" STRINGIZE(offset)"+8], "STRINGIZE(reg1)"\n"

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
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						"movss		xmm0, [esi]\n"
						"mulss		xmm0, [edi]\n"

						//STORE1( 0, xmm0, xmm1 )
                        "movss		[eax+0], xmm0\n"
                        ::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);

					return;
				}
				case 6: {		// 6x1 * 1x1
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						"movss		xmm0, [esi]\n"
						"shufps		xmm0, xmm0, 0x00\n"
						"movaps		xmm1, xmm0\n"
						"mulps		xmm0, [edi]\n"
						"mulps		xmm1, [edi+16]\n"
						STORE4( 0, xmm0, xmm2 )
						STORE2LO( 16, xmm1, xmm2 )
                        ::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);

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
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						"movss		xmm0, [esi]\n"
						"movss		xmm1, [esi+4]\n"
						"movss		xmm2, [edi]\n"
						"mulss		xmm2, xmm0\n"
						"movss		xmm3, [edi+4]\n"
						"mulss		xmm3, xmm1\n"
						"addss		xmm2, xmm3\n"
						//STORE1( 0, xmm2, xmm4 )
                        "movss		[eax+0], xmm2\n"
						"mulss		xmm0, [edi+8]\n"
						"mulss		xmm1, [edi+8+4]\n"
						"addss		xmm0, xmm1\n"
						//STORE1( 4, xmm0, xmm4 )
                        "movss		[eax+4], xmm0\n"
                        ::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);

					return;
				}
				case 6: {		// 6x2 * 2x1
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						"movlps		xmm7, [esi]\n"
						"shufps		xmm7, xmm7, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
						"movaps		xmm0, [edi]\n"
						"mulps		xmm0, xmm7\n"
						"movaps		xmm1, [edi+16]\n"
						"mulps		xmm1, xmm7\n"
						"movaps		xmm2, xmm0\n"
						"shufps		xmm0, xmm1, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/
						"shufps		xmm2, xmm1, 0xDD\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )=11011101*/
						"movaps		xmm3, [edi+32]\n"
						"addps		xmm0, xmm2\n"
						"mulps		xmm3, xmm7\n"
						STORE4( 0, xmm0, xmm4 )
						"shufps		xmm3, xmm3, 0xd8 \n" /*R_SHUFFLEPS( 0, 2, 1, 3 )=11011000*/
						"movhlps		xmm1, xmm3\n"
						"addps		xmm3, xmm1\n"
						STORE2LO( 16, xmm3, xmm4 )
                        ::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);

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
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						"movss		xmm0, [esi]\n"
						"movss		xmm4, [edi]\n"
						"mulss		xmm4, xmm0\n"
						"movss		xmm1, [esi+4]\n"
						"movss		xmm5, [edi+4]\n"
						"mulss		xmm5, xmm1\n"
						"addss		xmm4, xmm5\n"
						"movss		xmm2, [esi+8]\n"
						"movss		xmm6, [edi+8]\n"
						"mulss		xmm6, xmm2\n"
						"addss		xmm4, xmm6\n"
						"movss		xmm3, [edi+12]\n"
						"mulss		xmm3, xmm0\n"
						//STORE1( 0, xmm4, xmm7 )
                        "movss		[eax+0], xmm4\n"
						"movss		xmm5, [edi+12+4]\n"
						"mulss		xmm5, xmm1\n"
						"addss		xmm3, xmm5\n"
						"movss		xmm6, [edi+12+8]\n"
						"mulss		xmm6, xmm2\n"
						"addss		xmm3, xmm6\n"
						"mulss		xmm0, [edi+24]\n"
						"mulss		xmm1, [edi+24+4]\n"
						//STORE1( 4, xmm3, xmm7 )
						"movss		[eax+4], xmm3\n"
						"addss		xmm0, xmm1\n"
						"mulss		xmm2, [edi+24+8]\n"
						"addss		xmm0, xmm2\n"
						//STORE1( 8, xmm0, xmm7 )
                        "movss		[eax+8], xmm0\n"
						::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);

					return;
				}
				case 6: {		// 6x3 * 3x1
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						"movss		xmm5, [esi]\n"
						"shufps		xmm5, xmm5, 0x00\n"
						"movss		xmm6, [esi+4]\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"movss		xmm7, [esi+8]\n"
						"shufps		xmm7, xmm7, 0x00\n"
						"movaps		xmm0, [edi]	\n"							// xmm0 = 0, 1, 2, 3
						"movlps		xmm1, [edi+4*4]\n"
						"shufps		xmm1, xmm0, 0x94\n"  /*R_SHUFFLEPS( 0, 1, 1, 2 )=10010100 */	// xmm1 = 4, 5, 1, 2
						"movlps		xmm2, [edi+6*4]\n"
						"movhps		xmm2, [edi+8*4]\n"							// xmm2 = 6, 7, 8, 9
						"shufps		xmm0, xmm2, 0xCC \n" /*R_SHUFFLEPS( 0, 3, 0, 3 )=11001100 */	// xmm0 = 0, 3, 6, 9
						"mulps		xmm0, xmm5\n"
						"movlps		xmm3, [edi+10*4]\n"
						"shufps		xmm2, xmm3, 0x49\n"  /* R_SHUFFLEPS( 1, 2, 0, 1 )=01001001*/	// xmm2 = 7, 8, 10, 11
						"movaps		xmm3, xmm1\n"
						"shufps		xmm1, xmm2, 0x82\n" /*R_SHUFFLEPS( 2, 0, 0, 2 )=10000010 */	// xmm1 = 1, 4, 7, 10
						"mulps		xmm1, xmm6\n"
						"shufps		xmm3, xmm2, 0xD7\n" /*R_SHUFFLEPS( 3, 1, 1, 3 )=11010111  */	// xmm3 = 2, 5, 8, 11
						"mulps		xmm3, xmm7\n"
						"addps		xmm0, xmm1\n"
						"addps		xmm0, xmm3\n"
						STORE4( 0, xmm0, xmm4 )
						"movss		xmm1, [edi+12*4]\n"
						"mulss		xmm1, xmm5\n"
						"movss		xmm2, [edi+13*4]\n"
						"mulss		xmm2, xmm6\n"
						"movss		xmm3, [edi+14*4]\n"
						"mulss		xmm3, xmm7\n"
						"addss		xmm1, xmm2\n"
						"addss		xmm1, xmm3\n"
						//STORE1( 16, xmm1, xmm4 )
                        "movss		[eax+16], xmm1\n"
						"mulss		xmm5, [edi+15*4]\n"
						"mulss		xmm6, [edi+16*4]\n"
						"mulss		xmm7, [edi+17*4]\n"
						"addss		xmm5, xmm6\n"
						"addss		xmm5, xmm7\n"
						//STORE1( 20, xmm5, xmm4 )
                        "movss		[eax+20], xmm5\n"
                        ::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);

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
					asm (
						"mov			esi, %0\n"
						"mov			edi,%1\n"
						"mov			eax, %2\n"
						"movlps		xmm6, qword ptr [esi ]\n"
						"movlps		xmm0, qword ptr [edi ]\n"
						"shufps		xmm6, xmm6, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
						"movhps		xmm0, qword ptr [edi+16]\n"
						"mulps		xmm0, xmm6\n"
						"movlps		xmm7, qword ptr [esi+ 8]\n"
						"movlps		xmm2, qword ptr [edi+ 8]\n"
						"shufps		xmm7, xmm7, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
						"movhps		xmm2, qword ptr [edi+24]\n"
						"mulps		xmm2, xmm7\n"
						"movlps		xmm1, qword ptr [edi+32]\n"
						"movhps		xmm1, qword ptr [edi+48]\n"
						"mulps		xmm1, xmm6\n"
						"movlps		xmm3, qword ptr [edi+40]\n"
						"addps		xmm0, xmm2\n"
						"movhps		xmm3, qword ptr [edi+56]\n"
						"mulps		xmm3, xmm7\n"
						"movaps		xmm4, xmm0\n"
						"addps		xmm1, xmm3\n"
						"shufps		xmm4, xmm1, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/
						"shufps		xmm0, xmm1, 0xDD\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )=11011101*/
						"addps		xmm0, xmm4\n"
						STORE4( 0, xmm0, xmm2 )
                        ::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);

					return;
				}
				case 6: {		// 6x4 * 4x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm6, qword ptr [esi+ 0]\n"
						"movlps		xmm0, qword ptr [edi+ 0]\n"
						"shufps		xmm6, xmm6, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
						"movhps		xmm0, qword ptr [edi+16]\n"
						"mulps		xmm0, xmm6\n"
						"movlps		xmm7, qword ptr [esi+ 8]\n"
						"movlps		xmm2, qword ptr [edi+ 8]\n"
						"shufps		xmm7, xmm7, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
						"movhps		xmm2, qword ptr [edi+24]\n"
						"mulps		xmm2, xmm7\n"
						"movlps		xmm1, qword ptr [edi+32]\n"
						"movhps		xmm1, qword ptr [edi+48]\n"
						"mulps		xmm1, xmm6\n"
						"movlps		xmm3, qword ptr [edi+40]\n"
						"addps		xmm0, xmm2\n"
						"movhps		xmm3, qword ptr [edi+56]\n"
						"mulps		xmm3, xmm7\n"
						"movaps		xmm4, xmm0\n"
						"addps		xmm1, xmm3\n"
						"shufps		xmm4, xmm1, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/
						"shufps		xmm0, xmm1, 0xDD\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )=11011101*/
						"addps		xmm0, xmm4\n"
						"movlps		xmm1, qword ptr [edi+64]\n"
						"movhps		xmm1, qword ptr [edi+80]\n"
						STORE4( 0, xmm0, xmm4 )
						"mulps		xmm1, xmm6\n"
						"movlps		xmm2, qword ptr [edi+72]\n"
						"movhps		xmm2, qword ptr [edi+88]\n"
						"mulps		xmm2, xmm7\n"
						"addps		xmm1, xmm2\n"
						"shufps		xmm1, xmm1, 0xd8 \n" /*R_SHUFFLEPS( 0, 2, 1, 3 )=11011000*/
						"movhlps		xmm3, xmm1\n"
						"addps		xmm1, xmm3\n"
						STORE2LO( 16, xmm1, xmm4 )
                        ::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);

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
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movss		xmm0, [edi+5*4]\n"							// xmm0 =  5,  X,  X,  X
						"movhps		xmm0, [edi+0*4]\n"							// xmm0 =  5,  X,  0,  1
						"movss		xmm5, [edi+15*4]\n"						// xmm4 = 15,  X,  X,  X
						"movhps		xmm5, [edi+10*4]\n"						// xmm5 = 15,  X, 10, 11
						"movaps		xmm1, xmm0\n"								// xmm1 =  5,  X,  0,  1
						"shufps		xmm0, xmm5, 0x22\n" /*R_SHUFFLEPS( 2, 0, 2, 0 )=00100010*/	// xmm0 =  0,  5, 10, 15
						"movlps		xmm1, [edi+6*4]\n"							// xmm1 =  6,  7,  0,  1
						"movlps		xmm5, [edi+16*4]\n"						// xmm5 = 16, 17, 10, 11
						"movaps		xmm2, xmm1\n"								// xmm2 =  6,  7,  0,  1
						"shufps		xmm1, xmm5, 0x33\n" /*R_SHUFFLEPS( 3, 0, 3, 0 )=00110011*/	// xmm1 =  1,  6, 11, 16
						"movhps		xmm2, [edi+2*4]\n"							// xmm2 =  6,  7,  2,  3
						"movhps		xmm5, [edi+12*4]\n"						// xmm5 = 16, 17, 12, 13
						"movaps		xmm3, xmm2\n"								// xmm3 =  6,  7,  2,  3
						"shufps		xmm2, xmm5, 0x66 \n" /*R_SHUFFLEPS( 2, 1, 2, 1 )=01100110*/	// xmm2 =  2,  7, 12, 17
						"movlps		xmm3, [edi+8*4]\n"							// xmm3 =  8,  9,  2,  3
						"movlps		xmm5, [edi+18*4]\n"						// xmm5 = 18, 19, 12, 13
						"movss		xmm4, [edi+4*4]\n"							// xmm4 =  4,  X,  X,  X
						"movlhps		xmm4, xmm3\n"								// xmm4 =  4,  X,  8,  9
						"shufps		xmm3, xmm5, 0x33\n" /*R_SHUFFLEPS( 3, 0, 3, 0 )=00110011*/	// xmm3 =  3,  8, 13, 18
						"movhps		xmm5, [edi+14*4]\n"						// xmm6 = 18, 19, 14, 15
						"shufps		xmm4, xmm5, 0x6c\n" /*R_SHUFFLEPS( 0, 3, 2, 1 )=01101100*/	// xmm4 =  4,  9, 14, 19
						"movss		xmm7, [esi+0*4]\n"
						"shufps		xmm7, xmm7, 0\n"
						"mulps		xmm0, xmm7\n"
						"movss		xmm5, [esi+1*4]\n"
						"shufps		xmm5, xmm5, 0\n"
						"mulps		xmm1, xmm5\n"
						"addps		xmm0, xmm1\n"
						"movss		xmm6, [esi+2*4]\n"
						"shufps		xmm6, xmm6, 0\n"
						"mulps		xmm2, xmm6\n"
						"addps		xmm0, xmm2\n"
						"movss		xmm1, [esi+3*4]\n"
						"shufps		xmm1, xmm1, 0\n"
						"mulps		xmm3, xmm1\n"
						"addps		xmm0, xmm3\n"
						"movss		xmm2, [esi+4*4]\n"
						"shufps		xmm2, xmm2, 0\n"
						"mulps		xmm4, xmm2\n"
						"addps		xmm0, xmm4\n"
						"mulss		xmm7, [edi+20*4]\n"
						"mulss		xmm5, [edi+21*4]\n"
						"addps		xmm7, xmm5\n"
						"mulss		xmm6, [edi+22*4]\n"
						"addps		xmm7, xmm6\n"
						"mulss		xmm1, [edi+23*4]\n"
						"addps		xmm7, xmm1\n"
						"mulss		xmm2, [edi+24*4]\n"
						"addps		xmm7, xmm2\n"
						STORE4( 0, xmm0, xmm3 )
						//STORE1( 16, xmm7, xmm4 )
                        "movss		[eax+16], xmm7\n"
                        ::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);

					return;
				}
				case 6: {		// 6x5 * 5x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm6, [esi]\n"
						"shufps		xmm6, xmm6, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
						"movlps		xmm7, [esi+8]\n"
						"shufps		xmm7, xmm7, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
						"movlps		xmm0, [edi]\n"
						"movhps		xmm3, [edi+8]\n"
						"movaps		xmm1, [edi+16]\n"
						"movlps		xmm2, [edi+32]\n"
						"shufps		xmm0, xmm1, 0x94\n"  /*R_SHUFFLEPS( 0, 1, 1, 2 )=10010100 */	// xmm0 = 0, 1, 5, 6
						"shufps		xmm1, xmm2, 0x4C \n" /*R_SHUFFLEPS( 0, 3, 0, 1 )=01001100*/	// xmm1 = 4, 7, 8, 9
						"shufps		xmm3, xmm1, 0x9E \n" /*R_SHUFFLEPS( 2, 3, 1, 2 )=10011110 */	// xmm3 = 2, 3, 7, 8
						"mulps		xmm0, xmm6\n"
						"mulps		xmm3, xmm7\n"
						"movlps		xmm2, [edi+40]\n"
						"addps		xmm0, xmm3\n"								// xmm0 + xmm1
						"movhps		xmm5, [edi+40+8]\n"
						"movlps		xmm3, [edi+40+16]\n"
						"movhps		xmm3, [edi+40+24]\n"
						"movlps		xmm4, [edi+40+32]\n"
						"shufps		xmm2, xmm3, 0x94\n"  /*R_SHUFFLEPS( 0, 1, 1, 2 )=10010100 */	// xmm2 = 10, 11, 15, 16
						"shufps		xmm3, xmm4, 0x4C \n" /*R_SHUFFLEPS( 0, 3, 0, 1 )=01001100*/	// xmm3 = 14, 17, 18, 19
						"shufps		xmm5, xmm3, 0x9E \n" /*R_SHUFFLEPS( 2, 3, 1, 2 )=10011110 */	// xmm5 = 12, 13, 17, 18
						"mulps		xmm2, xmm6\n"
						"mulps		xmm5, xmm7\n"
						"addps		xmm2, xmm5\n"								// xmm2 + xmm3
						"movss		xmm5, [esi+16]\n"
						"shufps		xmm5, xmm5, 0x00\n"
						"movaps		xmm4, xmm0\n"
						"shufps		xmm0, xmm2, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/
						"shufps		xmm4, xmm2, 0xDD\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )=11011101*/
						"shufps		xmm1, xmm3, 0xCC \n" /*R_SHUFFLEPS( 0, 3, 0, 3 )=11001100 */
						"addps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"addps		xmm0, xmm1\n"
						STORE4( 0, xmm0, xmm2 )
						"movlps		xmm4, [edi+80]\n"
						"movhps		xmm3, [edi+80+8]\n"
						"movaps		xmm1, [edi+80+16]\n"
						"movlps		xmm2, [edi+80+32]\n"
						"shufps		xmm4, xmm1, 0x94\n"  /*R_SHUFFLEPS( 0, 1, 1, 2 )=10010100 */	// xmm4 = 20, 21, 25, 26
						"shufps		xmm1, xmm2, 0x4C \n" /*R_SHUFFLEPS( 0, 3, 0, 1 )=01001100*/	// xmm1 = 24, 27, 28, 29
						"shufps		xmm3, xmm1, 0x9E \n" /*R_SHUFFLEPS( 2, 3, 1, 2 )=10011110 */	// xmm3 = 22, 23, 27, 28
						"mulps		xmm4, xmm6\n"
						"mulps		xmm3, xmm7\n"
						"mulps		xmm1, xmm5\n"
						"addps		xmm4, xmm3\n"								// xmm4 + xmm1
						"shufps		xmm1, xmm4, 0x8C \n"  /* R_SHUFFLEPS( 0, 3, 0, 2 )=10001100 */
						"shufps		xmm4, xmm4, 0x0D \n" /* R_SHUFFLEPS( 1, 3, 0, 0 )=00001101 */
						"addps		xmm4, xmm1\n"
						"shufps		xmm1, xmm1, 0x4E\n"   /*R_SHUFFLEPS( 2, 3, 0, 1 )=01001110*/
						"addps		xmm4, xmm1\n"
						STORE2LO( 16, xmm4, xmm2 )
                        ::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);

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
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
                        "mov		eax, %2\n"
						"movss		xmm0, [esi]\n"
						"mulss		xmm0, [edi]\n"
						"movss		xmm1, [esi+4]\n"
						"mulss		xmm1, [edi+4]\n"
						"movss		xmm2, [esi+8]\n"
						"addss		xmm0, xmm1\n"
						"mulss		xmm2, [edi+8]\n"
						"movss		xmm3, [esi+12]\n"
						"addss		xmm0, xmm2\n"
						"mulss		xmm3, [edi+12]\n"
						"movss		xmm4, [esi+16]\n"
						"addss		xmm0, xmm3\n"
						"mulss		xmm4, [edi+16]\n"
						"movss		xmm5, [esi+20]\n"
						"addss		xmm0, xmm4\n"
						"mulss		xmm5, [edi+20]\n"
						"movss		xmm6, [esi+24]\n"
						"addss		xmm0, xmm5\n"
						"mulss		xmm6, [edi+24]\n"
						"addss		xmm0, xmm6\n"
						//STORE1( 0, xmm0, xmm7 )
                        "movss		[eax+0], xmm0\n"
                        ::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);

					return;
				}
				case 2: {		// 2x6 * 6x1
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						// load CVecXD
						"movlps		xmm4, [esi]\n"
						"movhps		xmm4, [esi+8]\n"
						"movlps		xmm5, [esi+16]\n"
						"movlhps		xmm5, xmm4\n"
						"movhlps		xmm6, xmm4\n"
						"movlhps		xmm6, xmm5\n"
						// row 0 "and 1
						"movaps		xmm0, [edi]\n"
						"movaps		xmm1, [edi+16]\n"
						"movaps		xmm2, [edi+32]\n"
						"mulps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"mulps		xmm2, xmm6\n"
						"movhlps	xmm3, xmm0\n"
						"movlhps	xmm3, xmm2\n"
						"addps		xmm1, xmm3\n"
						"shufps		xmm0, xmm2, 0xE4 \n"  /* R_SHUFFLEPS( 0, 1, 2, 3 )=11100100*/
						"addps		xmm1, xmm0\n"
						"shufps		xmm1, xmm1, 0xd8 \n"  /* R_SHUFFLEPS( 0, 2, 1, 3 )=11011000*/
						"movhlps	xmm0, xmm1\n"
						"addps		xmm0, xmm1\n"
						STORE2LO( 0, xmm0, xmm3 )
                        ::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);

					return;
				}
				case 3: {		// 3x6 * 6x1
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						// load CVecXD
						"movlps		xmm4, [esi]\n"
						"movhps		xmm4, [esi+8]\n"
						"movlps		xmm5, [esi+16]\n"
						"movlhps		xmm5, xmm4\n"
						"movhlps		xmm6, xmm4\n"
						"movlhps		xmm6, xmm5\n"
						// row 0 "and 1
						"movaps		xmm0, [edi]\n"
						"movaps		xmm1, [edi+16]\n"
						"movaps		xmm2, [edi+32]\n"
						"mulps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"mulps		xmm2, xmm6\n"
						"movhlps		xmm3, xmm0\n"
						"movlhps		xmm3, xmm2\n"
						"addps		xmm1, xmm3\n"
						"shufps		xmm0, xmm2, 0xE4 \n"  /* R_SHUFFLEPS( 0, 1, 2, 3 )=11100100*/
						"addps		xmm1, xmm0\n"
						"shufps		xmm1, xmm1, 0xd8 \n"  /* R_SHUFFLEPS( 0, 2, 1, 3 )=11011000*/
						"movhlps		xmm0, xmm1\n"
						"addps		xmm0, xmm1\n"
						STORE2LO( 0, xmm0, xmm3 )
						// row 2
						"movaps		xmm0, [edi+48]\n"
						"movaps		xmm1, [edi+48+16]\n"
						"mulps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"addps		xmm0, xmm1\n"
						"movhlps	xmm1, xmm0\n"
						"addps		xmm0, xmm1\n"
						"movaps		xmm1, xmm0\n"
						"shufps		xmm1, xmm1, 0x01 \n"  /* R_SHUFFLEPS( 1, 0, 0, 0 )=00000001 */
						"addss		xmm0, xmm1\n"
						//STORE1( 8, xmm0, xmm3 )
                        "movss		[eax+8], xmm0\n"
                        ::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);

					return;
				}
				case 4: {		// 4x6 * 6x1
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						// load CVecXD
						"movlps		xmm4, [esi]\n"
						"movhps		xmm4, [esi+8]\n"
						"movlps		xmm5, [esi+16]\n"
						"movlhps		xmm5, xmm4\n"
						"movhlps		xmm6, xmm4\n"
						"movlhps		xmm6, xmm5\n"
						// row 0 "and 1
						"movaps		xmm0, [edi]\n"
						"movaps		xmm1, [edi+16]\n"
						"movaps		xmm2, [edi+32]\n"
						"mulps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"mulps		xmm2, xmm6\n"
						"movhlps		xmm7, xmm0\n"
						"movlhps		xmm7, xmm2\n"
						"addps		xmm7, xmm1\n"
						"shufps		xmm0, xmm2, 0xE4 \n"  /* R_SHUFFLEPS( 0, 1, 2, 3 )=11100100*/
						"addps		xmm7, xmm0\n"
						// row 2 "and 3
						"movaps		xmm0, [edi+48]\n"
						"movaps		xmm1, [edi+48+16]\n"
						"movaps		xmm2, [edi+48+32]\n"
						"mulps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"mulps		xmm2, xmm6\n"
						"movhlps		xmm3, xmm0\n"
						"movlhps		xmm3, xmm2\n"
						"addps		xmm1, xmm3\n"
						"shufps		xmm0, xmm2, 0xE4 \n"  /* R_SHUFFLEPS( 0, 1, 2, 3 )=11100100*/
						"addps		xmm1, xmm0\n"
						// last 4 additions for the first 4 rows "and store result
						"movaps		xmm0, xmm7\n"
						"shufps		xmm7, xmm1, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/
						"shufps		xmm0, xmm1, 0xDD\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )=11011101*/
						"addps		xmm0, xmm7\n"
						STORE4( 0, xmm0, xmm4 )
                        ::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);

					return;
				}
				case 5: {		// 5x6 * 6x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						// load CVecXD
						"movlps		xmm4, [esi]\n"
						"movhps		xmm4, [esi+8]\n"
						"movlps		xmm5, [esi+16]\n"
						"movlhps		xmm5, xmm4\n"
						"movhlps		xmm6, xmm4\n"
						"movlhps		xmm6, xmm5\n"
						// row 0 "and 1
						"movaps		xmm0, [edi]\n"
						"movaps		xmm1, [edi+16]\n"
						"movaps		xmm2, [edi+32]\n"
						"mulps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"mulps		xmm2, xmm6\n"
						"movhlps	xmm7, xmm0\n"
						"movlhps	xmm7, xmm2\n"
						"addps		xmm7, xmm1\n"
						"shufps		xmm0, xmm2, 0xE4 \n"  /* R_SHUFFLEPS( 0, 1, 2, 3 )=11100100*/
						"addps		xmm7, xmm0\n"
						// row 2 "and 3
						"movaps		xmm0, [edi+48]\n"
						"movaps		xmm1, [edi+48+16]\n"
						"movaps		xmm2, [edi+48+32]\n"
						"mulps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"mulps		xmm2, xmm6\n"
						"movhlps	xmm3, xmm0\n"
						"movlhps	xmm3, xmm2\n"
						"addps		xmm1, xmm3\n"
						"shufps		xmm0, xmm2, 0xE4 \n"  /* R_SHUFFLEPS( 0, 1, 2, 3 )=11100100*/
						"addps		xmm1, xmm0\n"
						// last 4 additions for the first 4 rows "and store result
						"movaps		xmm0, xmm7\n"
						"shufps		xmm7, xmm1, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/
						"shufps		xmm0, xmm1, 0xDD\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )=11011101*/
						"addps		xmm0, xmm7\n"
						STORE4( 0, xmm0, xmm3 )
						// row 5
						"movaps		xmm0, [edi+96]\n"
						"movaps		xmm1, [edi+96+16]\n"
						"mulps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"addps		xmm0, xmm1\n"
						"movhlps	xmm1, xmm0\n"
						"addps		xmm0, xmm1\n"
						"movaps		xmm1, xmm0\n"
						"shufps		xmm1, xmm1, 0x01\n"
						"addss		xmm0, xmm1\n"
						//STORE1( 16, xmm0, xmm3 )
                        "movss		[eax+16], xmm0\n"
                        ::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);

					return;
				}
				case 6: {		// 6x6 * 6x1
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						"movlps		xmm7, qword ptr [esi]\n"
						"movlps		xmm6, qword ptr [esi+8]\n"
						"shufps		xmm7, xmm7, 0x44\n"
						"shufps		xmm6, xmm6, 0x44\n"
						"movlps		xmm0, qword ptr [edi    ]\n"
						"movhps		xmm0, qword ptr [edi+ 24]\n"
						"mulps		xmm0, xmm7\n"
						"movlps		xmm3, qword ptr [edi+  8]\n"
						"movhps		xmm3, qword ptr [edi+ 32]\n"
						"mulps		xmm3, xmm6\n"
						"movlps		xmm1, qword ptr [edi+ 48]\n"
						"movhps		xmm1, qword ptr [edi+ 72]\n"
						"mulps		xmm1, xmm7\n"
						"movlps		xmm2, qword ptr [edi+ 96]\n"
						"movhps		xmm2, qword ptr [edi+120]\n"
						"mulps		xmm2, xmm7\n"
						"movlps		xmm4, qword ptr [edi+ 56]\n"
						"movhps		xmm4, qword ptr [edi+ 80]\n"
						"movlps		xmm5, qword ptr [edi+104]\n"
						"movhps		xmm5, qword ptr [edi+128]\n"
						"mulps		xmm4, xmm6\n"
						"movlps		xmm7, qword ptr [esi+16]\n"
						"addps		xmm0, xmm3\n"
						"shufps		xmm7, xmm7, 0x44\n"
						"mulps		xmm5, xmm6\n"
						"addps		xmm1, xmm4\n"
						"movlps		xmm3, qword ptr [edi+ 16]\n"
						"movhps		xmm3, qword ptr [edi+ 40]\n"
						"addps		xmm2, xmm5\n"
						"movlps		xmm4, qword ptr [edi+ 64]\n"
						"movhps		xmm4, qword ptr [edi+ 88]\n"
						"mulps		xmm3, xmm7\n"
						"movlps		xmm5, qword ptr [edi+112]\n"
						"movhps		xmm5, qword ptr [edi+136]\n"
						"addps		xmm0, xmm3\n"
						"mulps		xmm4, xmm7\n"
						"mulps		xmm5, xmm7\n"
						"addps		xmm1, xmm4\n"
						"addps		xmm2, xmm5\n"
						"movaps		xmm6, xmm0\n"
						"shufps		xmm0, xmm1, 0x88\n"
						"shufps		xmm6, xmm1, 0xDD\n"
						"movaps		xmm7, xmm2\n"
						"shufps		xmm7, xmm2, 0x88\n"
						"shufps		xmm2, xmm2, 0xDD\n"
						"addps		xmm0, xmm6\n"
						"addps		xmm2, xmm7\n"
						STORE4( 0, xmm0, xmm3 )
						STORE2LO( 16, xmm2, xmm4 )
                        ::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);

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
	"movss		"STRINGIZE(reg2)", [eax+"STRINGIZE(offset)"]\n"		\
	"addss		"STRINGIZE(reg2)", "STRINGIZE(reg1)"\n"				\
	"movss		[eax+"STRINGIZE(offset)"], "STRINGIZE(reg2)"\n"
#define STORE2LO( offset, reg1, reg2 )		\
	"movlps	"STRINGIZE(reg2)", [eax+"STRINGIZE(offset)"]\n"		\
	"addps		"STRINGIZE(reg2)", "STRINGIZE(reg1)"\n"				\
	"movlps	[eax+"STRINGIZE(offset)"], "STRINGIZE(reg2)"\n"
#define STORE2HI( offset, reg1, reg2 )		\
	"movhps	"STRINGIZE(reg2)", [eax+"STRINGIZE(offset)"]\n"		\
	"addps		"STRINGIZE(reg2)", "STRINGIZE(reg1)"\n"				\
	"movhps	[eax+"STRINGIZE(offset)"], "STRINGIZE(reg2)"\n"
#define STORE4( offset, reg1, reg2 )		\
	"movlps	"STRINGIZE(reg2)", [eax+"STRINGIZE(offset)"]\n"		\
	"movhps	"STRINGIZE(reg2)", [eax+"STRINGIZE(offset)"+8]\n"	\
	"addps		"STRINGIZE(reg2)", "STRINGIZE(reg1)"\n"				\
	"movlps	[eax+"STRINGIZE(offset)"], "STRINGIZE(reg2)"\n"		\
	"movhps	[eax+"STRINGIZE(offset)"+8], "STRINGIZE(reg2)"\n"
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
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						"movss		xmm0, [esi]\n"
						"mulss		xmm0, [edi]\n"
						STORE1( 0, xmm0, xmm1 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 6: {		// 6x1 * 1x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movss		xmm0, [esi]\n"
						"shufps		xmm0, xmm0, 0x00\n"
						"movaps		xmm1, xmm0\n"
						"mulps		xmm0, [edi]\n"
						"mulps		xmm1, [edi+16]\n"
						STORE4( 0, xmm0, xmm2)
						STORE2LO( 16, xmm1, xmm2 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						"movss		xmm0, [esi]\n"
						"movss		xmm1, [esi+4]\n"
						"movss		xmm2, [edi]\n"
						"mulss		xmm2, xmm0\n"
						"movss		xmm3, [edi+4]\n"
						"mulss		xmm3, xmm1\n"
						"addss		xmm2, xmm3\n"
						STORE1( 0, xmm2, xmm4 )
						"mulss		xmm0, [edi+8]\n"
						"mulss		xmm1, [edi+8+4]\n"
						"addss		xmm0, xmm1\n"
						STORE1( 4, xmm0, xmm4 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 6: {		// 6x2 * 2x1
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						"movlps		xmm7, [esi]\n"
						"shufps		xmm7, xmm7, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
						"movaps		xmm0, [edi]\n"
						"mulps		xmm0, xmm7\n"
						"movaps		xmm1, [edi+16]\n"
						"mulps		xmm1, xmm7\n"
						"movaps		xmm2, xmm0\n"
						"shufps		xmm0, xmm1, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/
						"shufps		xmm2, xmm1, 0xDD\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )=11011101*/
						"movaps		xmm3, [edi+32]\n"
						"addps		xmm0, xmm2\n"
						"mulps		xmm3, xmm7\n"
						STORE4( 0, xmm0, xmm4 )
						"shufps		xmm3, xmm3, 0xd8 \n" /*R_SHUFFLEPS( 0, 2, 1, 3 )=11011000*/
						"movhlps		xmm1, xmm3\n"
						"addps		xmm3, xmm1\n"
						STORE2LO( 16, xmm3, xmm4 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						"movss		xmm0, [esi]\n"
						"movss		xmm4, [edi]\n"
						"mulss		xmm4, xmm0\n"
						"movss		xmm1, [esi+4]\n"
						"movss		xmm5, [edi+4]\n"
						"mulss		xmm5, xmm1\n"
						"addss		xmm4, xmm5\n"
						"movss		xmm2, [esi+8]\n"
						"movss		xmm6, [edi+8]\n"
						"mulss		xmm6, xmm2\n"
						"addss		xmm4, xmm6\n"
						"movss		xmm3, [edi+12]\n"
						"mulss		xmm3, xmm0\n"
						STORE1( 0, xmm4, xmm7 )
						"movss		xmm5, [edi+12+4]\n"
						"mulss		xmm5, xmm1\n"
						"addss		xmm3, xmm5\n"
						"movss		xmm6, [edi+12+8]\n"
						"mulss		xmm6, xmm2\n"
						"addss		xmm3, xmm6\n"
						"mulss		xmm0, [edi+24]\n"
						"mulss		xmm1, [edi+24+4]\n"
						STORE1( 4, xmm3, xmm7 )
						"addss		xmm0, xmm1\n"
						"mulss		xmm2, [edi+24+8]\n"
						"addss		xmm0, xmm2\n"
						STORE1( 8, xmm0, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 6: {		// 6x3 * 3x1
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						"movss		xmm5, [esi]\n"
						"shufps		xmm5, xmm5, 0x00\n"
						"movss		xmm6, [esi+4]\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"movss		xmm7, [esi+8]\n"
						"shufps		xmm7, xmm7, 0x00\n"
						"movaps		xmm0, [edi]\n"								// xmm0 = 0, 1, 2, 3
						"movlps		xmm1, [edi+4*4]\n"
						"shufps		xmm1, xmm0, 0x94\n"  /*R_SHUFFLEPS( 0, 1, 1, 2 )=10010100 */	// xmm1 = 4, 5, 1, 2
						"movlps		xmm2, [edi+6*4]\n"
						"movhps		xmm2, [edi+8*4]\n"							// xmm2 = 6, 7, 8, 9
						"shufps		xmm0, xmm2, 0xCC \n" /*R_SHUFFLEPS( 0, 3, 0, 3 )=11001100 */	// xmm0 = 0, 3, 6, 9
						"mulps		xmm0, xmm5\n"
						"movlps		xmm3, [edi+10*4]\n"
						"shufps		xmm2, xmm3, 0x49\n"  /* R_SHUFFLEPS( 1, 2, 0, 1 )=01001001*/	// xmm2 = 7, 8, 10, 11
						"movaps		xmm3, xmm1\n"
						"shufps		xmm1, xmm2, 0x82 \n"  /* R_SHUFFLEPS( 2, 0, 0, 2 )=10000010	*/ // xmm1 = 1, 4, 7, 10
						"mulps		xmm1, xmm6\n"
						"shufps		xmm3, xmm2, 0xD7 \n" /*R_SHUFFLEPS( 3, 1, 1, 3 )=11010111  */	// xmm3 = 2, 5, 8, 11
						"mulps		xmm3, xmm7\n"
						"addps		xmm0, xmm1\n"
						"addps		xmm0, xmm3\n"
						STORE4( 0, xmm0, xmm4 )
						"movss		xmm1, [edi+12*4]\n"
						"mulss		xmm1, xmm5\n"
						"movss		xmm2, [edi+13*4]\n"
						"mulss		xmm2, xmm6\n"
						"movss		xmm3, [edi+14*4]\n"
						"mulss		xmm3, xmm7\n"
						"addss		xmm1, xmm2\n"
						"addss		xmm1, xmm3\n"
						STORE1( 16, xmm1, xmm4 )
						"mulss		xmm5, [edi+15*4]\n"
						"mulss		xmm6, [edi+16*4]\n"
						"mulss		xmm7, [edi+17*4]\n"
						"addss		xmm5, xmm6\n"
						"addss		xmm5, xmm7\n"
						STORE1( 20, xmm5, xmm4 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						"movlps		xmm6, qword ptr [esi ]\n"
						"movlps		xmm0, qword ptr [edi ]\n"
						"shufps		xmm6, xmm6, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
						"movhps		xmm0, qword ptr [edi+16]\n"
						"mulps		xmm0, xmm6\n"
						"movlps		xmm7, qword ptr [esi+ 8]\n"
						"movlps		xmm2, qword ptr [edi+ 8]\n"
						"shufps		xmm7, xmm7, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
						"movhps		xmm2, qword ptr [edi+24]\n"
						"mulps		xmm2, xmm7\n"
						"movlps		xmm1, qword ptr [edi+32]\n"
						"movhps		xmm1, qword ptr [edi+48]\n"
						"mulps		xmm1, xmm6\n"
						"movlps		xmm3, qword ptr [edi+40]\n"
						"addps		xmm0, xmm2\n"
						"movhps		xmm3, qword ptr [edi+56]\n"
						"mulps		xmm3, xmm7\n"
						"movaps		xmm4, xmm0\n"
						"addps		xmm1, xmm3\n"
						"shufps		xmm4, xmm1, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/
						"shufps		xmm0, xmm1, 0xDD\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )=11011101*/
						"addps		xmm0, xmm4\n"
						STORE4( 0, xmm0, xmm2 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 6: {		// 6x4 * 4x1
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						"movlps		xmm6, qword ptr [esi+ 0]\n"
						"movlps		xmm0, qword ptr [edi+ 0]\n"
						"shufps		xmm6, xmm6, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
						"movhps		xmm0, qword ptr [edi+16]\n"
						"mulps		xmm0, xmm6\n"
						"movlps		xmm7, qword ptr [esi+ 8]\n"
						"movlps		xmm2, qword ptr [edi+ 8]\n"
						"shufps		xmm7, xmm7, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
						"movhps		xmm2, qword ptr [edi+24]\n"
						"mulps		xmm2, xmm7\n"
						"movlps		xmm1, qword ptr [edi+32]\n"
						"movhps		xmm1, qword ptr [edi+48]\n"
						"mulps		xmm1, xmm6\n"
						"movlps		xmm3, qword ptr [edi+40]\n"
						"addps		xmm0, xmm2\n"
						"movhps		xmm3, qword ptr [edi+56]\n"
						"mulps		xmm3, xmm7\n"
						"movaps		xmm4, xmm0\n"
						"addps		xmm1, xmm3\n"
						"shufps		xmm4, xmm1, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/
						"shufps		xmm0, xmm1, 0xDD\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )=11011101*/
						"addps		xmm0, xmm4\n"
						"movlps		xmm1, qword ptr [edi+64]\n"
						"movhps		xmm1, qword ptr [edi+80]\n"
						STORE4( 0, xmm0, xmm4 )
						"mulps		xmm1, xmm6\n"
						"movlps		xmm2, qword ptr [edi+72]\n"
						"movhps		xmm2, qword ptr [edi+88]\n"
						"mulps		xmm2, xmm7\n"
						"addps		xmm1, xmm2\n"
						"shufps		xmm1, xmm1, 0xd8 \n" /*R_SHUFFLEPS( 0, 2, 1, 3 )=11011000*/
						"movhlps		xmm3, xmm1\n"
						"addps		xmm1, xmm3\n"
						STORE2LO( 16, xmm1, xmm4 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						"movss		xmm0, [edi+5*4]\n"							// xmm0 =  5,  X,  X,  X
						"movhps		xmm0, [edi+0*4]\n"							// xmm0 =  5,  X,  0,  1
						"movss		xmm5, [edi+15*4]\n"						// xmm4 = 15,  X,  X,  X
						"movhps		xmm5, [edi+10*4]\n"						// xmm5 = 15,  X, 10, 11
						"movaps		xmm1, xmm0\n"								// xmm1 =  5,  X,  0,  1
						"shufps		xmm0, xmm5, 0x22\n" /*R_SHUFFLEPS( 2, 0, 2, 0 )=00100010*/	// xmm0 =  0,  5, 10, 15
						"movlps		xmm1, [edi+6*4]\n"							// xmm1 =  6,  7,  0,  1
						"movlps		xmm5, [edi+16*4]\n"						// xmm5 = 16, 17, 10, 11
						"movaps		xmm2, xmm1\n"								// xmm2 =  6,  7,  0,  1
						"shufps		xmm1, xmm5, 0x33\n" /*R_SHUFFLEPS( 3, 0, 3, 0 )=00110011*/	// xmm1 =  1,  6, 11, 16
						"movhps		xmm2, [edi+2*4]\n"							// xmm2 =  6,  7,  2,  3
						"movhps		xmm5, [edi+12*4]\n"						// xmm5 = 16, 17, 12, 13
						"movaps		xmm3, xmm2\n"								// xmm3 =  6,  7,  2,  3
						"shufps		xmm2, xmm5, 0x66 \n" /*R_SHUFFLEPS( 2, 1, 2, 1 )=01100110*/	// xmm2 =  2,  7, 12, 17
						"movlps		xmm3, [edi+8*4]\n"							// xmm3 =  8,  9,  2,  3
						"movlps		xmm5, [edi+18*4]\n"						// xmm5 = 18, 19, 12, 13
						"movss		xmm4, [edi+4*4]\n"							// xmm4 =  4,  X,  X,  X
						"movlhps		xmm4, xmm3\n"								// xmm4 =  4,  X,  8,  9
						"shufps		xmm3, xmm5, 0x33\n" /*R_SHUFFLEPS( 3, 0, 3, 0 )=00110011*/	// xmm3 =  3,  8, 13, 18
						"movhps		xmm5, [edi+14*4]\n"						// xmm6 = 18, 19, 14, 15
						"shufps		xmm4, xmm5, 0x6c\n" /*R_SHUFFLEPS( 0, 3, 2, 1 )=01101100*/	// xmm4 =  4,  9, 14, 19
						"movss		xmm7, [esi+0*4]\n"
						"shufps		xmm7, xmm7, 0\n"
						"mulps		xmm0, xmm7\n"
						"movss		xmm5, [esi+1*4]\n"
						"shufps		xmm5, xmm5, 0\n"
						"mulps		xmm1, xmm5\n"
						"addps		xmm0, xmm1\n"
						"movss		xmm6, [esi+2*4]\n"
						"shufps		xmm6, xmm6, 0\n"
						"mulps		xmm2, xmm6\n"
						"addps		xmm0, xmm2\n"
						"movss		xmm1, [esi+3*4]\n"
						"shufps		xmm1, xmm1, 0x00\n"
						"mulps		xmm3, xmm1\n"
						"addps		xmm0, xmm3\n"
						"movss		xmm2, [esi+4*4]\n"
						"shufps		xmm2, xmm2, 0\n"
						"mulps		xmm4, xmm2\n"
						"addps		xmm0, xmm4\n"
						"mulss		xmm7, [edi+20*4]\n"
						"mulss		xmm5, [edi+21*4]\n"
						"addps		xmm7, xmm5\n"
						"mulss		xmm6, [edi+22*4]\n"
						"addps		xmm7, xmm6\n"
						"mulss		xmm1, [edi+23*4]\n"
						"addps		xmm7, xmm1\n"
						"mulss		xmm2, [edi+24*4]\n"
						"addps		xmm7, xmm2\n"
						STORE4( 0, xmm0, xmm3 )
						STORE1( 16, xmm7, xmm4 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 6: {		// 6x5 * 5x1
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						"movlps		xmm6, [esi]\n"
						"shufps		xmm6, xmm6, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
						"movlps		xmm7, [esi+8]\n"
						"shufps		xmm7, xmm7, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
						"movlps		xmm0, [edi]\n"
						"movhps		xmm3, [edi+8]\n"
						"movaps		xmm1, [edi+16]\n"
						"movlps		xmm2, [edi+32]\n"
						"shufps		xmm0, xmm1, 0x94\n"  /*R_SHUFFLEPS( 0, 1, 1, 2 )=10010100 */	// xmm0 = 0, 1, 5, 6
						"shufps		xmm1, xmm2, 0x4C \n" /*R_SHUFFLEPS( 0, 3, 0, 1 )=01001100*/	// xmm1 = 4, 7, 8, 9
						"shufps		xmm3, xmm1, 0x9E \n" /*R_SHUFFLEPS( 2, 3, 1, 2 )=10011110 */	// xmm3 = 2, 3, 7, 8
						"mulps		xmm0, xmm6\n"
						"mulps		xmm3, xmm7\n"
						"movlps		xmm2, [edi+40]\n"
						"addps		xmm0, xmm3\n"								// xmm0 + xmm1
						"movhps		xmm5, [edi+40+8]\n"
						"movlps		xmm3, [edi+40+16]\n"
						"movhps		xmm3, [edi+40+24]\n"
						"movlps		xmm4, [edi+40+32]\n"
						"shufps		xmm2, xmm3, 0x94\n"  /*R_SHUFFLEPS( 0, 1, 1, 2 )=10010100 */	// xmm2 = 10, 11, 15, 16
						"shufps		xmm3, xmm4, 0x4C \n" /*R_SHUFFLEPS( 0, 3, 0, 1 )=01001100*/	// xmm3 = 14, 17, 18, 19
						"shufps		xmm5, xmm3, 0x9E \n" /*R_SHUFFLEPS( 2, 3, 1, 2 )=10011110 */	// xmm5 = 12, 13, 17, 18
						"mulps		xmm2, xmm6\n"
						"mulps		xmm5, xmm7\n"
						"addps		xmm2, xmm5\n"								// xmm2 + xmm3
						"movss		xmm5, [esi+16]\n"
						"shufps		xmm5, xmm5, 0x00\n"
						"movaps		xmm4, xmm0\n"
						"shufps		xmm0, xmm2, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/
						"shufps		xmm4, xmm2, 0xDD\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )=11011101*/
						"shufps		xmm1, xmm3, 0xCC \n" /*R_SHUFFLEPS( 0, 3, 0, 3 )=11001100 */
						"addps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"addps		xmm0, xmm1\n"
						STORE4( 0, xmm0, xmm2 )
						"movlps		xmm4, [edi+80]\n"
						"movhps		xmm3, [edi+80+8]\n"
						"movaps		xmm1, [edi+80+16]\n"
						"movlps		xmm2, [edi+80+32]\n"
						"shufps		xmm4, xmm1, 0x94\n"  /*R_SHUFFLEPS( 0, 1, 1, 2 )=10010100 */	// xmm4 = 20, 21, 25, 26
						"shufps		xmm1, xmm2, 0x4C \n" /*R_SHUFFLEPS( 0, 3, 0, 1 )=01001100*/	// xmm1 = 24, 27, 28, 29
						"shufps		xmm3, xmm1, 0x9E \n" /*R_SHUFFLEPS( 2, 3, 1, 2 )=10011110 */	// xmm3 = 22, 23, 27, 28
						"mulps		xmm4, xmm6\n"
						"mulps		xmm3, xmm7\n"
						"mulps		xmm1, xmm5\n"
						"addps		xmm4, xmm3\n"								// xmm4 + xmm1
						"shufps		xmm1, xmm4, 0x8C \n"  /* R_SHUFFLEPS( 0, 3, 0, 2 )=10001100 */
						"shufps		xmm4, xmm4, 0x0D \n" /* R_SHUFFLEPS( 1, 3, 0, 0 )=00001101 */
						"addps		xmm4, xmm1\n"
						"shufps		xmm1, xmm1, 0x4E\n"   /*R_SHUFFLEPS( 2, 3, 0, 1 )=01001110*/
						"addps		xmm4, xmm1\n"
						STORE2LO( 16, xmm4, xmm2 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
                        "mov		eax, %2\n"
						"movss		xmm0, [esi]\n"
						"mulss		xmm0, [edi]\n"
						"movss		xmm1, [esi+4]\n"
						"mulss		xmm1, [edi+4]\n"
						"movss		xmm2, [esi+8]\n"
						"addss		xmm0, xmm1\n"
						"mulss		xmm2, [edi+8]\n"
						"movss		xmm3, [esi+12]\n"
						"addss		xmm0, xmm2\n"
						"mulss		xmm3, [edi+12]\n"
						"movss		xmm4, [esi+16]\n"
						"addss		xmm0, xmm3\n"
						"mulss		xmm4, [edi+16]\n"
						"movss		xmm5, [esi+20]\n"
						"addss		xmm0, xmm4\n"
						"mulss		xmm5, [edi+20]\n"
						"movss		xmm6, [esi+24]\n"
						"addss		xmm0, xmm5\n"
						"mulss		xmm6, [edi+24]\n"
						"addss		xmm0, xmm6\n"
						STORE1( 0, xmm0, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 2: {		// 2x6 * 6x1
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						// load CVecXD
						"movlps		xmm4, [esi]\n"
						"movhps		xmm4, [esi+8]\n"
						"movlps		xmm5, [esi+16]\n"
						"movlhps		xmm5, xmm4\n"
						"movhlps		xmm6, xmm4\n"
						"movlhps		xmm6, xmm5\n"
						// row 0 "and 1
						"movaps		xmm0, [edi]\n"
						"movaps		xmm1, [edi+16]\n"
						"movaps		xmm2, [edi+32]\n"
						"mulps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"mulps		xmm2, xmm6\n"
						"movhlps		xmm3, xmm0\n"
						"movlhps		xmm3, xmm2\n"
						"addps		xmm1, xmm3\n"
						"shufps		xmm0, xmm2, 0xE4 \n"  /* R_SHUFFLEPS( 0, 1, 2, 3 )=11100100*/
						"addps		xmm1, xmm0\n"
						"shufps		xmm1, xmm1, 0xd8 \n" /*R_SHUFFLEPS( 0, 2, 1, 3 )=11011000*/
						"movhlps		xmm0, xmm1\n"
						"addps		xmm0, xmm1\n"
						STORE2LO( 0, xmm0, xmm3 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 3: {		// 3x6 * 6x1
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						// load CVecXD
						"movlps		xmm4, [esi]\n"
						"movhps		xmm4, [esi+8]\n"
						"movlps		xmm5, [esi+16]\n"
						"movlhps		xmm5, xmm4\n"
						"movhlps		xmm6, xmm4\n"
						"movlhps		xmm6, xmm5\n"
						// row 0 "and 1
						"movaps		xmm0, [edi]\n"
						"movaps		xmm1, [edi+16]\n"
						"movaps		xmm2, [edi+32]\n"
						"mulps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"mulps		xmm2, xmm6\n"
						"movhlps		xmm3, xmm0\n"
						"movlhps		xmm3, xmm2\n"
						"addps		xmm1, xmm3\n"
						"shufps		xmm0, xmm2, 0xE4 \n"  /* R_SHUFFLEPS( 0, 1, 2, 3 )=11100100*/
						"addps		xmm1, xmm0\n"
						"shufps		xmm1, xmm1, 0xd8 \n" /*R_SHUFFLEPS( 0, 2, 1, 3 )=11011000*/
						"movhlps		xmm0, xmm1\n"
						"addps		xmm0, xmm1\n"
						STORE2LO( 0, xmm0, xmm3 )
						// row 2
						"movaps		xmm0, [edi+48]\n"
						"movaps		xmm1, [edi+48+16]\n"
						"mulps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"addps		xmm0, xmm1\n"
						"movhlps		xmm1, xmm0\n"
						"addps		xmm0, xmm1\n"
						"movaps		xmm1, xmm0\n"
						"shufps		xmm1, xmm1, 0x01 \n"  /* R_SHUFFLEPS( 1, 0, 0, 0 )=00000001 */
						"addss		xmm0, xmm1\n"
						STORE1( 8, xmm0, xmm3 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 4: {		// 4x6 * 6x1
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						// load CVecXD
						"movlps		xmm4, [esi]\n"
						"movhps		xmm4, [esi+8]\n"
						"movlps		xmm5, [esi+16]\n"
						"movlhps		xmm5, xmm4\n"
						"movhlps		xmm6, xmm4\n"
						"movlhps		xmm6, xmm5\n"
						// row 0 "and 1
						"movaps		xmm0, [edi]\n"
						"movaps		xmm1, [edi+16]\n"
						"movaps		xmm2, [edi+32]\n"
						"mulps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"mulps		xmm2, xmm6\n"
						"movhlps		xmm7, xmm0\n"
						"movlhps		xmm7, xmm2\n"
						"addps		xmm7, xmm1\n"
						"shufps		xmm0, xmm2, 0xE4 \n"  /* R_SHUFFLEPS( 0, 1, 2, 3 )=11100100*/
						"addps		xmm7, xmm0\n"
						// row 2 "and 3
						"movaps		xmm0, [edi+48]\n"
						"movaps		xmm1, [edi+48+16]\n"
						"movaps		xmm2, [edi+48+32]\n"
						"mulps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"mulps		xmm2, xmm6\n"
						"movhlps		xmm3, xmm0\n"
						"movlhps		xmm3, xmm2\n"
						"addps		xmm1, xmm3\n"
						"shufps		xmm0, xmm2, 0xE4 \n"  /* R_SHUFFLEPS( 0, 1, 2, 3 )=11100100*/
						"addps		xmm1, xmm0\n"
						// last 4 additions for the first 4 rows "and store result
						"movaps		xmm0, xmm7\n"
						"shufps		xmm7, xmm1, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/
						"shufps		xmm0, xmm1, 0xDD\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )=11011101*/
						"addps		xmm0, xmm7\n"
						STORE4( 0, xmm0, xmm4 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 5: {		// 5x6 * 6x1
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						// load CVecXD
						"movlps		xmm4, [esi]\n"
						"movhps		xmm4, [esi+8]\n"
						"movlps		xmm5, [esi+16]\n"
						"movlhps		xmm5, xmm4\n"
						"movhlps		xmm6, xmm4\n"
						"movlhps		xmm6, xmm5\n"
						// row 0 "and 1
						"movaps		xmm0, [edi]\n"
						"movaps		xmm1, [edi+16]\n"
						"movaps		xmm2, [edi+32]\n"
						"mulps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"mulps		xmm2, xmm6\n"
						"movhlps		xmm7, xmm0\n"
						"movlhps		xmm7, xmm2\n"
						"addps		xmm7, xmm1\n"
						"shufps		xmm0, xmm2, 0xE4 \n"  /* R_SHUFFLEPS( 0, 1, 2, 3 )=11100100*/
						"addps		xmm7, xmm0\n"
						// row 2 "and 3
						"movaps		xmm0, [edi+48]\n"
						"movaps		xmm1, [edi+48+16]\n"
						"movaps		xmm2, [edi+48+32]\n"
						"mulps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"mulps		xmm2, xmm6\n"
						"movhlps		xmm3, xmm0\n"
						"movlhps		xmm3, xmm2\n"
						"addps		xmm1, xmm3\n"
						"shufps		xmm0, xmm2, 0xE4 \n"  /* R_SHUFFLEPS( 0, 1, 2, 3 )=11100100*/
						"addps		xmm1, xmm0\n"
						// last 4 additions for the first 4 rows "and store result
						"movaps		xmm0, xmm7\n"
						"shufps		xmm7, xmm1, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/
						"shufps		xmm0, xmm1, 0xDD\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )=11011101*/
						"addps		xmm0, xmm7\n"
						STORE4( 0, xmm0, xmm3 )
						// row 5
						"movaps		xmm0, [edi+96]\n"
						"movaps		xmm1, [edi+96+16]\n"
						"mulps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"addps		xmm0, xmm1\n"
						"movhlps		xmm1, xmm0\n"
						"addps		xmm0, xmm1\n"
						"movaps		xmm1, xmm0\n"
						"shufps		xmm1, xmm1, 0x01\n"
						"addss		xmm0, xmm1\n"
						STORE1( 16, xmm0, xmm3 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 6: {		// 6x6 * 6x1
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						"movlps		xmm7, qword ptr [esi]\n"
						"movlps		xmm6, qword ptr [esi+8]\n"
						"shufps		xmm7, xmm7, 0x44\n"
						"shufps		xmm6, xmm6, 0x44\n"
						"movlps		xmm0, qword ptr [edi    ]\n"
						"movhps		xmm0, qword ptr [edi+ 24]\n"
						"mulps		xmm0, xmm7\n"
						"movlps		xmm3, qword ptr [edi+  8]\n"
						"movhps		xmm3, qword ptr [edi+ 32]\n"
						"mulps		xmm3, xmm6\n"
						"movlps		xmm1, qword ptr [edi+ 48]\n"
						"movhps		xmm1, qword ptr [edi+ 72]\n"
						"mulps		xmm1, xmm7\n"
						"movlps		xmm2, qword ptr [edi+ 96]\n"
						"movhps		xmm2, qword ptr [edi+120]\n"
						"mulps		xmm2, xmm7\n"
						"movlps		xmm4, qword ptr [edi+ 56]\n"
						"movhps		xmm4, qword ptr [edi+ 80]\n"
						"movlps		xmm5, qword ptr [edi+104]\n"
						"movhps		xmm5, qword ptr [edi+128]\n"
						"mulps		xmm4, xmm6\n"
						"movlps		xmm7, qword ptr [esi+16]\n"
						"addps		xmm0, xmm3\n"
						"shufps		xmm7, xmm7, 0x44\n"
						"mulps		xmm5, xmm6\n"
						"addps		xmm1, xmm4\n"
						"movlps		xmm3, qword ptr [edi+ 16]\n"
						"movhps		xmm3, qword ptr [edi+ 40]\n"
						"addps		xmm2, xmm5\n"
						"movlps		xmm4, qword ptr [edi+ 64]\n"
						"movhps		xmm4, qword ptr [edi+ 88]\n"
						"mulps		xmm3, xmm7\n"
						"movlps		xmm5, qword ptr [edi+112]\n"
						"movhps		xmm5, qword ptr [edi+136]\n"
						"addps		xmm0, xmm3\n"
						"mulps		xmm4, xmm7\n"
						"mulps		xmm5, xmm7\n"
						"addps		xmm1, xmm4\n"
						"addps		xmm2, xmm5\n"
						"movaps		xmm6, xmm0\n"
						"shufps		xmm0, xmm1, 0x88\n"
						"shufps		xmm6, xmm1, 0xDD\n"
						"movaps		xmm7, xmm2\n"
						"shufps		xmm7, xmm2, 0x88\n"
						"shufps		xmm2, xmm2, 0xDD\n"
						"addps		xmm0, xmm6\n"
						"addps		xmm2, xmm7\n"
						STORE4( 0, xmm0, xmm3 )
						STORE2LO( 16, xmm2, xmm4 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
	"movss		"STRINGIZE(reg2)", [eax+"STRINGIZE(offset)"]\n"		\
	"subss		"STRINGIZE(reg2)", "STRINGIZE(reg1)"\n"				\
	"movss		[eax+"STRINGIZE(offset)"], "STRINGIZE(reg2)"\n"
#define STORE2LO( offset, reg1, reg2 )		\
	"movlps	"STRINGIZE(reg2)", [eax+"STRINGIZE(offset)"]\n"		\
	"subps		"STRINGIZE(reg2)", "STRINGIZE(reg1)"\n"				\
	"movlps	[eax+"STRINGIZE(offset)"], "STRINGIZE(reg2)"\n"
#define STORE2HI( offset, reg1, reg2 )		\
	"movhps	"STRINGIZE(reg2)", [eax+"STRINGIZE(offset)"]\n"		\
	"subps		"STRINGIZE(reg2)", "STRINGIZE(reg1)"\n"				\
	"movhps	[eax+"STRINGIZE(offset)"], "STRINGIZE(reg2)"\n"
#define STORE4( offset, reg1, reg2 )		\
	"movlps	"STRINGIZE(reg2)", [eax+"STRINGIZE(offset)"]\n"		\
	"movhps	"STRINGIZE(reg2)", [eax+"STRINGIZE(offset)"+8]\n"	\
	"subps		"STRINGIZE(reg2)", "STRINGIZE(reg1)"\n"				\
	"movlps	[eax+"STRINGIZE(offset)"], "STRINGIZE(reg2)"\n"		\
	"movhps	[eax+"STRINGIZE(offset)"+8], "STRINGIZE(reg2)"\n"
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
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						"movss		xmm0, [esi]\n"
						"mulss		xmm0, [edi]\n"
						STORE1( 0, xmm0, xmm1 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 6: {		// 6x1 * 1x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movss		xmm0, [esi]\n"
						"shufps		xmm0, xmm0, 0x00\n"
						"movaps		xmm1, xmm0\n"
						"mulps		xmm0, [edi]\n"
						"mulps		xmm1, [edi+16]\n"
						STORE4( 0, xmm0, xmm2 )
						STORE2LO( 16, xmm1, xmm2 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movss		xmm0, [esi]\n"
						"movss		xmm1, [esi+4]\n"
						"movss		xmm2, [edi]\n"
						"mulss		xmm2, xmm0\n"
						"movss		xmm3, [edi+4]\n"
						"mulss		xmm3, xmm1\n"
						"addss		xmm2, xmm3\n"
						STORE1( 0, xmm2, xmm4 )
						"mulss		xmm0, [edi+8]\n"
						"mulss		xmm1, [edi+8+4]\n"
						"addss		xmm0, xmm1\n"
						STORE1( 4, xmm0, xmm4 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 6: {		// 6x2 * 2x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm7, [esi]\n"
						"shufps		xmm7, xmm7, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
						"movaps		xmm0, [edi]\n"
						"mulps		xmm0, xmm7\n"
						"movaps		xmm1, [edi+16]\n"
						"mulps		xmm1, xmm7\n"
						"movaps		xmm2, xmm0\n"
						"shufps		xmm0, xmm1, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/
						"shufps		xmm2, xmm1, 0xDD\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )=11011101*/
						"movaps		xmm3, [edi+32]\n"
						"addps		xmm0, xmm2\n"
						"mulps		xmm3, xmm7\n"
						STORE4( 0, xmm0, xmm4 )
						"shufps		xmm3, xmm3, 0xd8 \n" /*R_SHUFFLEPS( 0, 2, 1, 3 )=11011000*/
						"movhlps		xmm1, xmm3\n"
						"addps		xmm3, xmm1\n"
						STORE2LO( 16, xmm3, xmm4 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movss		xmm0, [esi]\n"
						"movss		xmm4, [edi]\n"
						"mulss		xmm4, xmm0\n"
						"movss		xmm1, [esi+4]\n"
						"movss		xmm5, [edi+4]\n"
						"mulss		xmm5, xmm1\n"
						"addss		xmm4, xmm5\n"
						"movss		xmm2, [esi+8]\n"
						"movss		xmm6, [edi+8]\n"
						"mulss		xmm6, xmm2\n"
						"addss		xmm4, xmm6\n"
						"movss		xmm3, [edi+12]\n"
						"mulss		xmm3, xmm0\n"
						STORE1( 0, xmm4, xmm7 )
						"movss		xmm5, [edi+12+4]\n"
						"mulss		xmm5, xmm1\n"
						"addss		xmm3, xmm5\n"
						"movss		xmm6, [edi+12+8]\n"
						"mulss		xmm6, xmm2\n"
						"addss		xmm3, xmm6\n"
						"mulss		xmm0, [edi+24]\n"
						"mulss		xmm1, [edi+24+4]\n"
						STORE1( 4, xmm3, xmm7 )
						"addss		xmm0, xmm1\n"
						"mulss		xmm2, [edi+24+8]\n"
						"addss		xmm0, xmm2\n"
						STORE1( 8, xmm0, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 6: {		// 6x3 * 3x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movss		xmm5, [esi]\n"
						"shufps		xmm5, xmm5, 0x00\n"
						"movss		xmm6, [esi+4]\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"movss		xmm7, [esi+8]\n"
						"shufps		xmm7, xmm7, 0x00\n"
						"movaps		xmm0, [edi]\n"								// xmm0 = 0, 1, 2, 3
						"movlps		xmm1, [edi+4*4]\n"
						"shufps		xmm1, xmm0, 0x94\n"  /*R_SHUFFLEPS( 0, 1, 1, 2 )=10010100 */	// xmm1 = 4, 5, 1, 2
						"movlps		xmm2, [edi+6*4]\n"
						"movhps		xmm2, [edi+8*4]\n"							// xmm2 = 6, 7, 8, 9
						"shufps		xmm0, xmm2, 0xCC \n" /*R_SHUFFLEPS( 0, 3, 0, 3 )=11001100 */	// xmm0 = 0, 3, 6, 9
						"mulps		xmm0, xmm5\n"
						"movlps		xmm3, [edi+10*4]\n"
						"shufps		xmm2, xmm3, 0x49\n"  /* R_SHUFFLEPS( 1, 2, 0, 1 )=01001001*/	// xmm2 = 7, 8, 10, 11
						"movaps		xmm3, xmm1\n"
						"shufps		xmm1, xmm2, 0x82 \n"  /* R_SHUFFLEPS( 2, 0, 0, 2 )=10000010 */	// xmm1 = 1, 4, 7, 10
						"mulps		xmm1, xmm6\n"
						"shufps		xmm3, xmm2, 0xD7 \n" /*R_SHUFFLEPS( 3, 1, 1, 3 )=11010111  */	// xmm3 = 2, 5, 8, 11
						"mulps		xmm3, xmm7\n"
						"addps		xmm0, xmm1\n"
						"addps		xmm0, xmm3\n"
						STORE4( 0, xmm0, xmm4 )
						"movss		xmm1, [edi+12*4]\n"
						"mulss		xmm1, xmm5\n"
						"movss		xmm2, [edi+13*4]\n"
						"mulss		xmm2, xmm6\n"
						"movss		xmm3, [edi+14*4]\n"
						"mulss		xmm3, xmm7\n"
						"addss		xmm1, xmm2\n"
						"addss		xmm1, xmm3\n"
						STORE1( 16, xmm1, xmm4 )
						"mulss		xmm5, [edi+15*4]\n"
						"mulss		xmm6, [edi+16*4]\n"
						"mulss		xmm7, [edi+17*4]\n"
						"addss		xmm5, xmm6\n"
						"addss		xmm5, xmm7\n"
						STORE1( 20, xmm5, xmm4 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm6, qword ptr [esi ]\n"
						"movlps		xmm0, qword ptr [edi ]\n"
						"shufps		xmm6, xmm6, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
						"movhps		xmm0, qword ptr [edi+16]\n"
						"mulps		xmm0, xmm6\n"
						"movlps		xmm7, qword ptr [esi+ 8]\n"
						"movlps		xmm2, qword ptr [edi+ 8]\n"
						"shufps		xmm7, xmm7, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
						"movhps		xmm2, qword ptr [edi+24]\n"
						"mulps		xmm2, xmm7\n"
						"movlps		xmm1, qword ptr [edi+32]\n"
						"movhps		xmm1, qword ptr [edi+48]\n"
						"mulps		xmm1, xmm6\n"
						"movlps		xmm3, qword ptr [edi+40]\n"
						"addps		xmm0, xmm2\n"
						"movhps		xmm3, qword ptr [edi+56]\n"
						"mulps		xmm3, xmm7\n"
						"movaps		xmm4, xmm0\n"
						"addps		xmm1, xmm3\n"
						"shufps		xmm4, xmm1, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/
						"shufps		xmm0, xmm1, 0xDD\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )=11011101*/
						"addps		xmm0, xmm4\n"
						STORE4( 0, xmm0, xmm2)
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 6: {		// 6x4 * 4x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm6, qword ptr [esi+ 0]\n"
						"movlps		xmm0, qword ptr [edi+ 0]\n"
						"shufps		xmm6, xmm6, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
						"movhps		xmm0, qword ptr [edi+16]\n"
						"mulps		xmm0, xmm6\n"
						"movlps		xmm7, qword ptr [esi+ 8]\n"
						"movlps		xmm2, qword ptr [edi+ 8]\n"
						"shufps		xmm7, xmm7, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
						"movhps		xmm2, qword ptr [edi+24]\n"
						"mulps		xmm2, xmm7\n"
						"movlps		xmm1, qword ptr [edi+32]\n"
						"movhps		xmm1, qword ptr [edi+48]\n"
						"mulps		xmm1, xmm6\n"
						"movlps		xmm3, qword ptr [edi+40]\n"
						"addps		xmm0, xmm2\n"
						"movhps		xmm3, qword ptr [edi+56]\n"
						"mulps		xmm3, xmm7\n"
						"movaps		xmm4, xmm0\n"
						"addps		xmm1, xmm3\n"
						"shufps		xmm4, xmm1, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/
						"shufps		xmm0, xmm1, 0xDD\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )=11011101*/
						"addps		xmm0, xmm4\n"
						"movlps		xmm1, qword ptr [edi+64]\n"
						"movhps		xmm1, qword ptr [edi+80]\n"
						STORE4( 0, xmm0, xmm4 )
						"mulps		xmm1, xmm6\n"
						"movlps		xmm2, qword ptr [edi+72]\n"
						"movhps		xmm2, qword ptr [edi+88]\n"
						"mulps		xmm2, xmm7\n"
						"addps		xmm1, xmm2\n"
						"shufps		xmm1, xmm1, 0xd8 \n" /*R_SHUFFLEPS( 0, 2, 1, 3 )=11011000*/
						"movhlps		xmm3, xmm1\n"
						"addps		xmm1, xmm3\n"
						STORE2LO( 16, xmm1, xmm4 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movss		xmm0, [edi+5*4]	\n"						// xmm0 =  5,  X,  X,  X
						"movhps		xmm0, [edi+0*4]	\n"						// xmm0 =  5,  X,  0,  1
						"movss		xmm5, [edi+15*4]\n"						// xmm4 = 15,  X,  X,  X
						"movhps		xmm5, [edi+10*4]\n"						// xmm5 = 15,  X, 10, 11
						"movaps		xmm1, xmm0\n"								// xmm1 =  5,  X,  0,  1
						"shufps		xmm0, xmm5, 0x22\n" /*R_SHUFFLEPS( 2, 0, 2, 0 )=00100010*/	// xmm0 =  0,  5, 10, 15
						"movlps		xmm1, [edi+6*4]	\n"						// xmm1 =  6,  7,  0,  1
						"movlps		xmm5, [edi+16*4]\n"						// xmm5 = 16, 17, 10, 11
						"movaps		xmm2, xmm1\n"								// xmm2 =  6,  7,  0,  1
						"shufps		xmm1, xmm5, 0x33\n" /*R_SHUFFLEPS( 3, 0, 3, 0 )=00110011*/	// xmm1 =  1,  6, 11, 16
						"movhps		xmm2, [edi+2*4]	\n"						// xmm2 =  6,  7,  2,  3
						"movhps		xmm5, [edi+12*4]\n"						// xmm5 = 16, 17, 12, 13
						"movaps		xmm3, xmm2\n"								// xmm3 =  6,  7,  2,  3
						"shufps		xmm2, xmm5, 0x66 \n" /*R_SHUFFLEPS( 2, 1, 2, 1 )=01100110*/	// xmm2 =  2,  7, 12, 17
						"movlps		xmm3, [edi+8*4]	\n"						// xmm3 =  8,  9,  2,  3
						"movlps		xmm5, [edi+18*4]\n"						// xmm5 = 18, 19, 12, 13
						"movss		xmm4, [edi+4*4]	\n"						// xmm4 =  4,  X,  X,  X
						"movlhps		xmm4, xmm3\n"								// xmm4 =  4,  X,  8,  9
						"shufps		xmm3, xmm5, 0x33\n" /*R_SHUFFLEPS( 3, 0, 3, 0 )=00110011*/	// xmm3 =  3,  8, 13, 18
						"movhps		xmm5, [edi+14*4]\n"						// xmm6 = 18, 19, 14, 15
						"shufps		xmm4, xmm5, 0x6c\n" /*R_SHUFFLEPS( 0, 3, 2, 1 )=01101100*/	// xmm4 =  4,  9, 14, 19
						"movss		xmm7, [esi+0*4]\n"
						"shufps		xmm7, xmm7, 0\n"
						"mulps		xmm0, xmm7\n"
						"movss		xmm5, [esi+1*4]\n"
						"shufps		xmm5, xmm5, 0\n"
						"mulps		xmm1, xmm5\n"
						"addps		xmm0, xmm1\n"
						"movss		xmm6, [esi+2*4]\n"
						"shufps		xmm6, xmm6, 0\n"
						"mulps		xmm2, xmm6\n"
						"addps		xmm0, xmm2\n"
						"movss		xmm1, [esi+3*4]\n"
						"shufps		xmm1, xmm1, 0\n"
						"mulps		xmm3, xmm1\n"
						"addps		xmm0, xmm3\n"
						"movss		xmm2, [esi+4*4]\n"
						"shufps		xmm2, xmm2, 0\n"
						"mulps		xmm4, xmm2\n"
						"addps		xmm0, xmm4\n"
						"mulss		xmm7, [edi+20*4]\n"
						"mulss		xmm5, [edi+21*4]\n"
						"addps		xmm7, xmm5\n"
						"mulss		xmm6, [edi+22*4]\n"
						"addps		xmm7, xmm6\n"
						"mulss		xmm1, [edi+23*4]\n"
						"addps		xmm7, xmm1\n"
						"mulss		xmm2, [edi+24*4]\n"
						"addps		xmm7, xmm2\n"
						STORE4( 0, xmm0, xmm3 )
						STORE1( 16, xmm7, xmm4 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 6: {		// 6x5 * 5x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm6, [esi]\n"
						"shufps		xmm6, xmm6, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
						"movlps		xmm7, [esi+8]\n"
						"shufps		xmm7, xmm7, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
						"movlps		xmm0, [edi]\n"
						"movhps		xmm3, [edi+8]\n"
						"movaps		xmm1, [edi+16]\n"
						"movlps		xmm2, [edi+32]\n"
						"shufps		xmm0, xmm1, 0x94\n"  /*R_SHUFFLEPS( 0, 1, 1, 2 )=10010100 */	// xmm0 = 0, 1, 5, 6
						"shufps		xmm1, xmm2, 0x4C \n" /*R_SHUFFLEPS( 0, 3, 0, 1 )=01001100*/	// xmm1 = 4, 7, 8, 9
						"shufps		xmm3, xmm1, 0x9E \n" /*R_SHUFFLEPS( 2, 3, 1, 2 )=10011110 */	// xmm3 = 2, 3, 7, 8
						"mulps		xmm0, xmm6\n"
						"mulps		xmm3, xmm7\n"
						"movlps		xmm2, [edi+40]\n"
						"addps		xmm0, xmm3\n"								// xmm0 + xmm1
						"movhps		xmm5, [edi+40+8]\n"
						"movlps		xmm3, [edi+40+16]\n"
						"movhps		xmm3, [edi+40+24]\n"
						"movlps		xmm4, [edi+40+32]\n"
						"shufps		xmm2, xmm3, 0x94\n"  /*R_SHUFFLEPS( 0, 1, 1, 2 )=10010100 */	// xmm2 = 10, 11, 15, 16
						"shufps		xmm3, xmm4, 0x4C \n" /*R_SHUFFLEPS( 0, 3, 0, 1 )=01001100 */	// xmm3 = 14, 17, 18, 19
						"shufps		xmm5, xmm3, 0x9E \n" /*R_SHUFFLEPS( 2, 3, 1, 2 )=10011110 */	// xmm5 = 12, 13, 17, 18
						"mulps		xmm2, xmm6\n"
						"mulps		xmm5, xmm7\n"
						"addps		xmm2, xmm5\n"								// xmm2 + xmm3
						"movss		xmm5, [esi+16]\n"
						"shufps		xmm5, xmm5, 0x00\n"
						"movaps		xmm4, xmm0\n"
						"shufps		xmm0, xmm2, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/
						"shufps		xmm4, xmm2, 0xDD\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )=11011101*/
						"shufps		xmm1, xmm3, 0xCC \n" /*R_SHUFFLEPS( 0, 3, 0, 3 )=11001100 */
						"addps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"addps		xmm0, xmm1\n"
						STORE4( 0, xmm0, xmm2)
						"movlps		xmm4, [edi+80]\n"
						"movhps		xmm3, [edi+80+8]\n"
						"movaps		xmm1, [edi+80+16]\n"
						"movlps		xmm2, [edi+80+32]\n"
						"shufps		xmm4, xmm1, 0x94\n"  /*R_SHUFFLEPS( 0, 1, 1, 2 )=10010100 */	// xmm4 = 20, 21, 25, 26
						"shufps		xmm1, xmm2, 0x4C \n" /*R_SHUFFLEPS( 0, 3, 0, 1 )=01001100*/	// xmm1 = 24, 27, 28, 29
						"shufps		xmm3, xmm1, 0x9E \n" /*R_SHUFFLEPS( 2, 3, 1, 2 )=10011110 */	// xmm3 = 22, 23, 27, 28
						"mulps		xmm4, xmm6\n"
						"mulps		xmm3, xmm7\n"
						"mulps		xmm1, xmm5\n"
						"addps		xmm4, xmm3\n"								// xmm4 + xmm1
						"shufps		xmm1, xmm4, 0x8C \n"  /* R_SHUFFLEPS( 0, 3, 0, 2 )=10001100 */
						"shufps		xmm4, xmm4, 0x0D \n" /* R_SHUFFLEPS( 1, 3, 0, 0 )=00001101 */
						"addps		xmm4, xmm1\n"
						"shufps		xmm1, xmm1, 0x4E\n"   /*R_SHUFFLEPS( 2, 3, 0, 1 )=01001110*/
						"addps		xmm4, xmm1\n"
						STORE2LO( 16, xmm4, xmm2 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
                        "mov			eax, %2\n"
						"movss		xmm0, [esi]\n"
						"mulss		xmm0, [edi]\n"
						"movss		xmm1, [esi+4]\n"
						"mulss		xmm1, [edi+4]\n"
						"movss		xmm2, [esi+8]\n"
						"addss		xmm0, xmm1\n"
						"mulss		xmm2, [edi+8]\n"
						"movss		xmm3, [esi+12]\n"
						"addss		xmm0, xmm2\n"
						"mulss		xmm3, [edi+12]\n"
						"movss		xmm4, [esi+16]\n"
						"addss		xmm0, xmm3\n"
						"mulss		xmm4, [edi+16]\n"
						"movss		xmm5, [esi+20]\n"
						"addss		xmm0, xmm4\n"
						"mulss		xmm5, [edi+20]\n"
						"movss		xmm6, [esi+24]\n"
						"addss		xmm0, xmm5\n"
						"mulss		xmm6, [edi+24]\n"
						"addss		xmm0, xmm6\n"
						STORE1( 0, xmm0, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 2: {		// 2x6 * 6x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						// load CVecXD
						"movlps		xmm4, [esi]\n"
						"movhps		xmm4, [esi+8]\n"
						"movlps		xmm5, [esi+16]\n"
						"movlhps		xmm5, xmm4\n"
						"movhlps		xmm6, xmm4\n"
						"movlhps		xmm6, xmm5\n"
						// row 0 "and 1
						"movaps		xmm0, [edi]\n"
						"movaps		xmm1, [edi+16]\n"
						"movaps		xmm2, [edi+32]\n"
						"mulps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"mulps		xmm2, xmm6\n"
						"movhlps		xmm3, xmm0\n"
						"movlhps		xmm3, xmm2\n"
						"addps		xmm1, xmm3\n"
						"shufps		xmm0, xmm2, 0xE4 \n"  /* R_SHUFFLEPS( 0, 1, 2, 3 )=11100100*/
						"addps		xmm1, xmm0\n"
						"shufps		xmm1, xmm1, 0xd8 \n" /*R_SHUFFLEPS( 0, 2, 1, 3 )=11011000*/
						"movhlps		xmm0, xmm1\n"
						"addps		xmm0, xmm1\n"
						STORE2LO( 0, xmm0, xmm3 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 3: {		// 3x6 * 6x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						// load CVecXD
						"movlps		xmm4, [esi]\n"
						"movhps		xmm4, [esi+8]\n"
						"movlps		xmm5, [esi+16]\n"
						"movlhps		xmm5, xmm4\n"
						"movhlps		xmm6, xmm4\n"
						"movlhps		xmm6, xmm5\n"
						// row 0 "and 1
						"movaps		xmm0, [edi]\n"
						"movaps		xmm1, [edi+16]\n"
						"movaps		xmm2, [edi+32]\n"
						"mulps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"mulps		xmm2, xmm6\n"
						"movhlps		xmm3, xmm0\n"
						"movlhps		xmm3, xmm2\n"
						"addps		xmm1, xmm3\n"
						"shufps		xmm0, xmm2, 0xE4 \n"  /* R_SHUFFLEPS( 0, 1, 2, 3 )=11100100*/
						"addps		xmm1, xmm0\n"
						"shufps		xmm1, xmm1, 0xd8 \n" /*R_SHUFFLEPS( 0, 2, 1, 3 )=11011000*/
						"movhlps		xmm0, xmm1\n"
						"addps		xmm0, xmm1\n"
						STORE2LO( 0, xmm0, xmm3 )
						// row 2
						"movaps		xmm0, [edi+48]\n"
						"movaps		xmm1, [edi+48+16]\n"
						"mulps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"addps		xmm0, xmm1\n"
						"movhlps		xmm1, xmm0\n"
						"addps		xmm0, xmm1\n"
						"movaps		xmm1, xmm0\n"
						"shufps		xmm1, xmm1, 0x01 \n"  /* R_SHUFFLEPS( 1, 0, 0, 0 )=00000001 */
						"addss		xmm0, xmm1\n"
						STORE1( 8, xmm0, xmm3 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 4: {		// 4x6 * 6x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						// load CVecXD
						"movlps		xmm4, [esi]\n"
						"movhps		xmm4, [esi+8]\n"
						"movlps		xmm5, [esi+16]\n"
						"movlhps		xmm5, xmm4\n"
						"movhlps		xmm6, xmm4\n"
						"movlhps		xmm6, xmm5\n"
						// row 0 "and 1
						"movaps		xmm0, [edi]\n"
						"movaps		xmm1, [edi+16]\n"
						"movaps		xmm2, [edi+32]\n"
						"mulps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"mulps		xmm2, xmm6\n"
						"movhlps		xmm7, xmm0\n"
						"movlhps		xmm7, xmm2\n"
						"addps		xmm7, xmm1\n"
						"shufps		xmm0, xmm2, 0xE4 \n"  /* R_SHUFFLEPS( 0, 1, 2, 3 )=11100100*/
						"addps		xmm7, xmm0\n"
						// row 2 "and 3
						"movaps		xmm0, [edi+48]\n"
						"movaps		xmm1, [edi+48+16]\n"
						"movaps		xmm2, [edi+48+32]\n"
						"mulps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"mulps		xmm2, xmm6\n"
						"movhlps		xmm3, xmm0\n"
						"movlhps		xmm3, xmm2\n"
						"addps		xmm1, xmm3\n"
						"shufps		xmm0, xmm2, 0xE4 \n"  /* R_SHUFFLEPS( 0, 1, 2, 3 )=11100100*/
						"addps		xmm1, xmm0\n"
						// last 4 additions for the first 4 rows "and store result
						"movaps		xmm0, xmm7\n"
						"shufps		xmm7, xmm1, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/
						"shufps		xmm0, xmm1, 0xDD\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )=11011101*/
						"addps		xmm0, xmm7\n"
						STORE4( 0, xmm0, xmm4 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 5: {		// 5x6 * 6x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						// load CVecXD
						"movlps		xmm4, [esi]\n"
						"movhps		xmm4, [esi+8]\n"
						"movlps		xmm5, [esi+16]\n"
						"movlhps		xmm5, xmm4\n"
						"movhlps		xmm6, xmm4\n"
						"movlhps		xmm6, xmm5\n"
						// row 0 "and 1
						"movaps		xmm0, [edi]\n"
						"movaps		xmm1, [edi+16]\n"
						"movaps		xmm2, [edi+32]\n"
						"mulps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"mulps		xmm2, xmm6\n"
						"movhlps		xmm7, xmm0\n"
						"movlhps		xmm7, xmm2\n"
						"addps		xmm7, xmm1\n"
						"shufps		xmm0, xmm2, 0xE4 \n"  /* R_SHUFFLEPS( 0, 1, 2, 3 )=11100100*/
						"addps		xmm7, xmm0\n"
						// row 2 "and 3
						"movaps		xmm0, [edi+48]\n"
						"movaps		xmm1, [edi+48+16]\n"
						"movaps		xmm2, [edi+48+32]\n"
						"mulps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"mulps		xmm2, xmm6\n"
						"movhlps		xmm3, xmm0\n"
						"movlhps		xmm3, xmm2\n"
						"addps		xmm1, xmm3\n"
						"shufps		xmm0, xmm2, 0xE4 \n"  /* R_SHUFFLEPS( 0, 1, 2, 3 )=11100100*/
						"addps		xmm1, xmm0\n"
						// last 4 additions for the first 4 rows "and store result
						"movaps		xmm0, xmm7\n"
						"shufps		xmm7, xmm1, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/
						"shufps		xmm0, xmm1, 0xDD\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )=11011101*/
						"addps		xmm0, xmm7\n"
						STORE4( 0, xmm0, xmm3 )
						// row 5
						"movaps		xmm0, [edi+96]\n"
						"movaps		xmm1, [edi+96+16]\n"
						"mulps		xmm0, xmm4\n"
						"mulps		xmm1, xmm5\n"
						"addps		xmm0, xmm1\n"
						"movhlps		xmm1, xmm0\n"
						"addps		xmm0, xmm1\n"
						"movaps		xmm1, xmm0\n"
						"shufps		xmm1, xmm1, 0x01\n"
						"addss		xmm0, xmm1\n"
						STORE1( 16, xmm0, xmm3 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 6: {		// 6x6 * 6x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm7, qword ptr [esi]\n"
						"movlps		xmm6, qword ptr [esi+8]\n"
						"shufps		xmm7, xmm7, 0x44\n"
						"shufps		xmm6, xmm6, 0x44\n"
						"movlps		xmm0, qword ptr [edi    ]\n"
						"movhps		xmm0, qword ptr [edi+ 24]\n"
						"mulps		xmm0, xmm7\n"
						"movlps		xmm3, qword ptr [edi+  8]\n"
						"movhps		xmm3, qword ptr [edi+ 32]\n"
						"mulps		xmm3, xmm6\n"
						"movlps		xmm1, qword ptr [edi+ 48]\n"
						"movhps		xmm1, qword ptr [edi+ 72]\n"
						"mulps		xmm1, xmm7\n"
						"movlps		xmm2, qword ptr [edi+ 96]\n"
						"movhps		xmm2, qword ptr [edi+120]\n"
						"mulps		xmm2, xmm7\n"
						"movlps		xmm4, qword ptr [edi+ 56]\n"
						"movhps		xmm4, qword ptr [edi+ 80]\n"
						"movlps		xmm5, qword ptr [edi+104]\n"
						"movhps		xmm5, qword ptr [edi+128]\n"
						"mulps		xmm4, xmm6\n"
						"movlps		xmm7, qword ptr [esi+16]\n"
						"addps		xmm0, xmm3\n"
						"shufps		xmm7, xmm7, 0x44\n"
						"mulps		xmm5, xmm6\n"
						"addps		xmm1, xmm4\n"
						"movlps		xmm3, qword ptr [edi+ 16]\n"
						"movhps		xmm3, qword ptr [edi+ 40]\n"
						"addps		xmm2, xmm5\n"
						"movlps		xmm4, qword ptr [edi+ 64]\n"
						"movhps		xmm4, qword ptr [edi+ 88]\n"
						"mulps		xmm3, xmm7\n"
						"movlps		xmm5, qword ptr [edi+112]\n"
						"movhps		xmm5, qword ptr [edi+136]\n"
						"addps		xmm0, xmm3\n"
						"mulps		xmm4, xmm7\n"
						"mulps		xmm5, xmm7\n"
						"addps		xmm1, xmm4\n"
						"addps		xmm2, xmm5\n"
						"movaps		xmm6, xmm0\n"
						"shufps		xmm0, xmm1, 0x88\n"
						"shufps		xmm6, xmm1, 0xDD\n"
						"movaps		xmm7, xmm2\n"
						"shufps		xmm7, xmm2, 0x88\n"
						"shufps		xmm2, xmm2, 0xDD\n"
						"addps		xmm0, xmm6\n"
						"addps		xmm2, xmm7\n"
						STORE4( 0, xmm0, xmm3 )
						STORE2LO( 16, xmm2, xmm4 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
	"movss		[eax+"STRINGIZE(offset)"], "STRINGIZE(reg1)"\n"
#define STORE2LO( offset, reg1, reg2 )		\
	"movlps	[eax+"STRINGIZE(offset)"], "STRINGIZE(reg1)"\n"
#define STORE2HI( offset, reg1, reg2 )		\
	"movhps	[eax+"STRINGIZE(offset)"], "STRINGIZE(reg1)"\n"
#define STORE4( offset, reg1, reg2 )		\
	"movlps	[eax+"STRINGIZE(offset)"], "STRINGIZE(reg1)"\n"		\
	"movhps	[eax+"STRINGIZE(offset)"+8], "STRINGIZE(reg1)"\n"
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
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movss		xmm0, [esi]\n"
						"shufps		xmm0, xmm0, 0x00\n"
						"movaps		xmm1, xmm0\n"
						"mulps		xmm0, [edi]\n"
						"mulps		xmm1, [edi+16]\n"
						STORE4( 0, xmm0, xmm2)
						STORE2LO( 16, xmm1, xmm3 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm0, [esi]\n"
						"movaps		xmm1, xmm0\n"
						"shufps		xmm0, xmm0, 0x00\n"
						"shufps		xmm1, xmm1, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"movaps		xmm2, [edi]\n"
						"mulps		xmm2, xmm0\n"
						"movlps		xmm3, [edi+24]\n"
						"movhps		xmm3, [edi+32]\n"
						"mulps		xmm3, xmm1\n"
						"addps		xmm2, xmm3\n"
						"shufps		xmm0, xmm1, 0x00\n"
						"movlps		xmm4, [edi+16]\n"
						"movhps		xmm4, [edi+40]\n"
						"mulps		xmm4, xmm0\n"
						"movhlps	xmm3, xmm4\n"
						"addps		xmm3, xmm4\n"
						STORE4( 0, xmm2, xmm5 )
						STORE2LO( 16, xmm3, xmm6 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm0, [esi+0*4]\n"
						"movss		xmm1, [esi+2*4]\n"
						"movlps		xmm3, [edi+(0*6+0)*4]\n"
						"movhps		xmm3, [edi+(0*6+2)*4]\n"
						"movaps		xmm4, xmm0\n"
						"shufps		xmm4, xmm4, 0x00\n"
						"mulps		xmm3, xmm4\n"
						"movlps		xmm5, [edi+(1*6+0)*4]\n"
						"movhps		xmm5, [edi+(1*6+2)*4]\n"
						"movaps		xmm6, xmm0\n"
						"shufps		xmm6, xmm6, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						"movlps		xmm4, [edi+(2*6+0)*4]\n"
						"movhps		xmm4, [edi+(2*6+2)*4]\n"
						"shufps		xmm1, xmm1, 0x00\n"
						"mulps		xmm4, xmm1\n"
						"addps		xmm3, xmm4\n"
						STORE4( 0, xmm3, xmm7 )
						"shufps		xmm0, xmm0, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"movlps		xmm3, [edi+(0*6+4)*4]\n"
						"movhps		xmm3, [edi+(1*6+4)*4]\n"
						"mulps		xmm3, xmm0\n"
						"movhlps		xmm4, xmm3\n"
						"addps		xmm3, xmm4\n"
						"movlps		xmm5, [edi+(2*6+4)*4]\n"
						"mulps		xmm5, xmm1\n"
						"addps		xmm3, xmm5\n"
						STORE2LO( 16, xmm3, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm0, [esi+0*4]\n"
						"movlps		xmm1, [esi+2*4]\n"
						"movaps		xmm3, xmm0\n"
						"shufps		xmm3, xmm3, 0x00\n"
						"mulps		xmm3, [edi+(0*6+0)*4]\n"
						"movlps		xmm5, [edi+(1*6+0)*4]\n"
						"movhps		xmm5, [edi+(1*6+2)*4]\n"
						"movaps		xmm6, xmm0\n"
						"shufps		xmm6, xmm6, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						"movlps		xmm4, [edi+(2*6+0)*4]\n"
						"movhps		xmm4, [edi+(2*6+2)*4]\n"
						"movaps		xmm6, xmm1\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm4, xmm6\n"
						"addps		xmm3, xmm4\n"
						"movlps		xmm5, [edi+(3*6+0)*4]\n"
						"movhps		xmm5, [edi+(3*6+2)*4]\n"
						"movaps		xmm6, xmm1\n"
						"shufps		xmm6, xmm6, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						STORE4( 0, xmm3, xmm7 )
						"shufps		xmm0, xmm0, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"shufps		xmm1, xmm1, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"movlps		xmm3, [edi+(0*6+4)*4]\n"
						"movhps		xmm3, [edi+(1*6+4)*4]\n"
						"mulps		xmm3, xmm0\n"
						"movlps		xmm4, [edi+(2*6+4)*4]\n"
						"movhps		xmm4, [edi+(3*6+4)*4]\n"
						"mulps		xmm4, xmm1\n"
						"addps		xmm3, xmm4\n"
						"movhlps		xmm4, xmm3\n"
						"addps		xmm3, xmm4\n"
						STORE2LO( 16, xmm3, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm0, [esi+0*4]\n"
						"movlps		xmm1, [esi+2*4]\n"
						"movss		xmm2, [esi+4*4]\n"
						"movaps		xmm3, xmm0\n"
						"shufps		xmm3, xmm3, 0x00\n"
						"mulps		xmm3, [edi+(0*6+0)*4]\n"
						"movlps		xmm5, [edi+(1*6+0)*4]\n"
						"movhps		xmm5, [edi+(1*6+2)*4]\n"
						"movaps		xmm6, xmm0\n"
						"shufps		xmm6, xmm6, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						"movaps		xmm6, xmm1\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm6, [edi+(2*6+0)*4]\n"
						"addps		xmm3, xmm6\n"
						"movlps		xmm5, [edi+(3*6+0)*4]\n"
						"movhps		xmm5, [edi+(3*6+2)*4]\n"
						"movaps		xmm6, xmm1\n"
						"shufps		xmm6, xmm6, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						"shufps		xmm2, xmm2, 0x00\n"
						"movaps		xmm4, xmm2\n"
						"mulps		xmm4, [edi+(4*6+0)*4]\n"
						"addps		xmm3, xmm4\n"
						STORE4( 0, xmm3, xmm7 )
						"shufps		xmm0, xmm0, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"shufps		xmm1, xmm1, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"movlps		xmm3, [edi+(0*6+4)*4]\n"
						"movhps		xmm3, [edi+(1*6+4)*4]\n"
						"mulps		xmm3, xmm0\n"
						"movlps		xmm4, [edi+(2*6+4)*4]\n"
						"movhps		xmm4, [edi+(3*6+4)*4]\n"
						"mulps		xmm4, xmm1\n"
						"addps		xmm3, xmm4\n"
						"movhlps		xmm4, xmm3\n"
						"addps		xmm3, xmm4\n"
						"movlps		xmm5, [edi+(4*6+4)*4]\n"
						"mulps		xmm5, xmm2\n"
						"addps		xmm3, xmm5\n"
						STORE2LO( 16, xmm3, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm0, [esi]\n"
						"movhps		xmm0, [esi+8]\n"
						"movlps		xmm1, [esi+16]\n"
						"mulps		xmm0, [edi]\n"
						"mulps		xmm1, [edi+16]\n"
						"shufps		xmm1, xmm0, 0xB4 \n" /*R_SHUFFLEPS( 0, 1, 3, 2 )=10110100 */
						"addps		xmm0, xmm1\n"
						"movhlps		xmm2, xmm0\n"
						"addss		xmm2, xmm0\n"
						"shufps		xmm0, xmm0, 0x01 \n"  /* R_SHUFFLEPS( 1, 0, 0, 0 )=00000001 */
						"addss		xmm2, xmm0\n"
						STORE1( 0, xmm2, xmm3 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 2: {		// 6x2 * 6x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm0, [esi+0*4]\n"
						"shufps		xmm0, xmm0, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"movaps		xmm6, [edi+0*4]\n"
						"mulps		xmm6, xmm0\n"
						"movlps		xmm1, [esi+2*4]\n"
						"shufps		xmm1, xmm1, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"movaps		xmm7, [edi+4*4]\n"
						"mulps		xmm7, xmm1\n"
						"addps		xmm6, xmm7\n"
						"movlps		xmm2, [esi+4*4]\n"
						"shufps		xmm2, xmm2, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"movaps		xmm7, [edi+8*4]\n"
						"mulps		xmm7, xmm2\n"
						"addps		xmm6, xmm7\n"
						"movhlps		xmm3, xmm6\n"
						"addps		xmm3, xmm6\n"
						STORE2LO( 0, xmm3, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 3: {		// 6x3 * 6x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movss		xmm0, [edi+(0*3+2)*4]\n"
						"movhps		xmm0, [edi+(0*3+0)*4]\n"
						"shufps		xmm0, xmm0, 0x36\n"  /*R_SHUFFLEPS( 2, 1, 3, 0 )=00110110*/
						"movss		xmm6, [esi+0*4]\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm6, xmm0\n"
						"movss		xmm1, [edi+(1*3+0)*4]\n"
						"movhps		xmm1, [edi+(1*3+1)*4]\n"
						"movss		xmm7, [esi+1*4]\n"
						"shufps		xmm7, xmm7, 0x00\n"
						"mulps		xmm7, xmm1\n"
						"addps		xmm6, xmm7\n"
						"movss		xmm2, [edi+(2*3+2)*4]\n"
						"movhps		xmm2, [edi+(2*3+0)*4]\n"
						"shufps		xmm2, xmm2, 0x36\n"  /*R_SHUFFLEPS( 2, 1, 3, 0 )=00110110*/
						"movss		xmm7, [esi+2*4]\n"
						"shufps		xmm7, xmm7, 0x00\n"
						"mulps		xmm7, xmm2\n"
						"addps		xmm6, xmm7\n"
						"movss		xmm3, [edi+(3*3+0)*4]\n"
						"movhps		xmm3, [edi+(3*3+1)*4]\n"
						"movss		xmm7, [esi+3*4]\n"
						"shufps		xmm7, xmm7, 0x00\n"
						"mulps		xmm7, xmm3\n"
						"addps		xmm6, xmm7\n"
						"movss		xmm4, [edi+(4*3+2)*4]\n"
						"movhps		xmm4, [edi+(4*3+0)*4]\n"
						"shufps		xmm4, xmm4, 0x36\n"  /*R_SHUFFLEPS( 2, 1, 3, 0 )=00110110*/
						"movss		xmm7, [esi+4*4]\n"
						"shufps		xmm7, xmm7, 0x00\n"
						"mulps		xmm7, xmm4\n"
						"addps		xmm6, xmm7\n"
						"movss		xmm5, [edi+(5*3+0)*4]\n"
						"movhps		xmm5, [edi+(5*3+1)*4]\n"
						"movss		xmm7, [esi+5*4]\n"
						"shufps		xmm7, xmm7, 0x00\n"
						"mulps		xmm7, xmm5\n"
						"addps		xmm6, xmm7\n"
						STORE1( 0, xmm6, xmm7 )
						STORE2HI( 4, xmm6, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 4: {		// 6x4 * 6x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm3, [edi+(0*4+0)*4]\n"
						"movhps		xmm3, [edi+(0*4+2)*4]\n"
						"movss		xmm4, [esi+0*4]\n"
						"shufps		xmm4, xmm4, 0x00\n"
						"mulps		xmm3, xmm4\n"
						"movlps		xmm5, [edi+(1*4+0)*4]\n"
						"movhps		xmm5, [edi+(1*4+2)*4]\n"
						"movss		xmm6, [esi+1*4]\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						"movlps		xmm4, [edi+(2*4+0)*4]\n"
						"movhps		xmm4, [edi+(2*4+2)*4]\n"
						"movss		xmm6, [esi+2*4]\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm4, xmm6\n"
						"addps		xmm3, xmm4\n"
						"movlps		xmm5, [edi+(3*4+0)*4]\n"
						"movhps		xmm5, [edi+(3*4+2)*4]\n"
						"movss		xmm6, [esi+3*4]\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						"movlps		xmm4, [edi+(4*4+0)*4]\n"
						"movhps		xmm4, [edi+(4*4+2)*4]\n"
						"movss		xmm6, [esi+4*4]\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm4, xmm6\n"
						"addps		xmm3, xmm4\n"
						"movlps		xmm5, [edi+(5*4+0)*4]\n"
						"movhps		xmm5, [edi+(5*4+2)*4]\n"
						"movss		xmm6, [esi+5*4]\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						STORE4( 0, xmm3, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 5: {		// 6x5 * 6x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm6, [edi+(0*5+0)*4]\n"
						"movhps		xmm6, [edi+(0*5+2)*4]\n"
						"movss		xmm0, [esi+0*4]\n"
						"shufps		xmm0, xmm0, 0x00\n"
						"mulps		xmm6, xmm0\n"
						"movlps		xmm7, [edi+(1*5+0)*4]\n"
						"movhps		xmm7, [edi+(1*5+2)*4]\n"
						"movss		xmm1, [esi+1*4]\n"
						"shufps		xmm1, xmm1, 0x00\n"
						"mulps		xmm7, xmm1\n"
						"addps		xmm6, xmm7\n"
						"movlps		xmm7, [edi+(2*5+0)*4]\n"
						"movhps		xmm7, [edi+(2*5+2)*4]\n"
						"movss		xmm2, [esi+2*4]\n"
						"shufps		xmm2, xmm2, 0x00\n"
						"mulps		xmm7, xmm2\n"
						"addps		xmm6, xmm7\n"
						"movlps		xmm7, [edi+(3*5+0)*4]\n"
						"movhps		xmm7, [edi+(3*5+2)*4]\n"
						"movss		xmm3, [esi+3*4]\n"
						"shufps		xmm3, xmm3, 0x00\n"
						"mulps		xmm7, xmm3\n"
						"addps		xmm6, xmm7\n"
						"movlps		xmm7, [edi+(4*5+0)*4]\n"
						"movhps		xmm7, [edi+(4*5+2)*4]\n"
						"movss		xmm4, [esi+4*4]\n"
						"shufps		xmm4, xmm4, 0x00\n"
						"mulps		xmm7, xmm4\n"
						"addps		xmm6, xmm7\n"
						"movlps		xmm7, [edi+(5*5+0)*4]\n"
						"movhps		xmm7, [edi+(5*5+2)*4]\n"
						"movss		xmm5, [esi+5*4]\n"
						"shufps		xmm5, xmm5, 0x00\n"
						"mulps		xmm7, xmm5\n"
						"addps		xmm6, xmm7\n"
						STORE4( 0, xmm6, xmm7 )
						"movss		xmm6, [edi+(0*5+4)*4]\n"
						"mulss		xmm6, xmm0\n"
						"movss		xmm7, [edi+(1*5+4)*4]\n"
						"mulss		xmm7, xmm1\n"
						"addss		xmm6, xmm7\n"
						"movss		xmm7, [edi+(2*5+4)*4]\n"
						"mulss		xmm7, xmm2\n"
						"addss		xmm6, xmm7\n"
						"movss		xmm7, [edi+(3*5+4)*4]\n"
						"mulss		xmm7, xmm3\n"
						"addss		xmm6, xmm7\n"
						"movss		xmm7, [edi+(4*5+4)*4]\n"
						"mulss		xmm7, xmm4\n"
						"addss		xmm6, xmm7\n"
						"movss		xmm7, [edi+(5*5+4)*4]\n"
						"mulss		xmm7, xmm5\n"
						"addss		xmm6, xmm7\n"
						STORE1( 16, xmm6, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 6: {		// 6x6 * 6x1
					asm (
						"mov		esi, %0\n"
						"mov		edi, %1\n"
						"mov		eax, %2\n"
						"movlps		xmm0, [esi+0*4]\n"
						"movlps		xmm1, [esi+2*4]\n"
						"movlps		xmm2, [esi+4*4]\n"
						"movaps		xmm3, xmm0\n"
						"shufps		xmm3, xmm3, 0x00\n"
						"mulps		xmm3, [edi+(0*6+0)*4]\n"
						"movlps		xmm5, [edi+(1*6+0)*4]\n"
						"movhps		xmm5, [edi+(1*6+2)*4]\n"
						"movaps		xmm6, xmm0\n"
						"shufps		xmm6, xmm6, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						"movaps		xmm6, xmm1\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm6, [edi+(2*6+0)*4]\n"
						"addps		xmm3, xmm6\n"
						"movaps		xmm6, xmm1\n"
						"shufps		xmm6, xmm6, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"movlps		xmm5, [edi+(3*6+0)*4]\n"
						"movhps		xmm5, [edi+(3*6+2)*4]\n"
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						"movaps		xmm6, xmm2\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm6, [edi+(4*6+0)*4]\n"
						"addps		xmm3, xmm6\n"
						"movaps		xmm6, xmm2\n"
						"shufps		xmm6, xmm6, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"movlps		xmm5, [edi+(5*6+0)*4]\n"
						"movhps		xmm5, [edi+(5*6+2)*4]\n"
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						STORE4( 0, xmm3, xmm7 )
						"shufps		xmm0, xmm0, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"shufps		xmm1, xmm1, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"shufps		xmm2, xmm2, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"movlps		xmm3, [edi+(0*6+4)*4]\n"
						"movhps		xmm3, [edi+(1*6+4)*4]\n"
						"mulps		xmm3, xmm0\n"
						"movlps		xmm4, [edi+(2*6+4)*4]\n"
						"movhps		xmm4, [edi+(3*6+4)*4]\n"
						"mulps		xmm4, xmm1\n"
						"addps		xmm3, xmm4\n"
						"movlps		xmm5, [edi+(4*6+4)*4]\n"
						"movhps		xmm5, [edi+(5*6+4)*4]\n"
						"mulps		xmm5, xmm2\n"
						"addps		xmm3, xmm5\n"
						"movhlps	xmm4, xmm3\n"
						"addps		xmm3, xmm4\n"
						STORE2LO( 16, xmm3, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
	"movss		"STRINGIZE(reg2)", [eax+"STRINGIZE(offset)"]\n"		\
	"addss		"STRINGIZE(reg2)", "STRINGIZE(reg1)"\n"				\
	"movss		[eax+"STRINGIZE(offset)"], "STRINGIZE(reg2)"\n"
#define STORE2LO( offset, reg1, reg2 )		\
	"movlps	"STRINGIZE(reg2)", [eax+"STRINGIZE(offset)"]\n"		\
	"addps		"STRINGIZE(reg2)", "STRINGIZE(reg1)"\n"				\
	"movlps	[eax+"STRINGIZE(offset)"], "STRINGIZE(reg2)"\n"
#define STORE2HI( offset, reg1, reg2 )		\
	"movhps	"STRINGIZE(reg2)", [eax+"STRINGIZE(offset)"]\n"		\
	"addps		"STRINGIZE(reg2)", "STRINGIZE(reg1)"\n"				\
	"movhps	[eax+"STRINGIZE(offset)"], "STRINGIZE(reg2)"\n"

#define STORE4( offset, reg1, reg2 )		\
	"movlps	"STRINGIZE(reg2)", [eax+"STRINGIZE(offset)"]\n"		\
	"movhps	"STRINGIZE(reg2)", [eax+"STRINGIZE(offset)"+8]\n"	\
	"addps		"STRINGIZE(reg2)", "STRINGIZE(reg1)"\n"			\
	"movlps	[eax+"STRINGIZE(offset)"], "STRINGIZE(reg2)"\n"		\
	"movhps	[eax+"STRINGIZE(offset)"+8], "STRINGIZE(reg2)"\n"

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
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movss		xmm0, [esi]\n"
						"shufps		xmm0, xmm0, 0x00\n"
						"movaps		xmm1, xmm0\n"
						"mulps		xmm0, [edi]\n"
						"mulps		xmm1, [edi+16]\n"
						STORE4( 0, xmm0, xmm2 )
						STORE2LO( 16, xmm1, xmm3 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm0, [esi]\n"
						"movaps		xmm1, xmm0\n"
						"shufps		xmm0, xmm0, 0x00\n"
						"shufps		xmm1, xmm1, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"movaps		xmm2, [edi]\n"
						"mulps		xmm2, xmm0\n"
						"movlps		xmm3, [edi+24]\n"
						"movhps		xmm3, [edi+32]\n"
						"mulps		xmm3, xmm1\n"
						"addps		xmm2, xmm3\n"
						"shufps		xmm0, xmm1, 0x00\n"
						"movlps		xmm4, [edi+16]\n"
						"movhps		xmm4, [edi+40]\n"
						"mulps		xmm4, xmm0\n"
						"movhlps	xmm3, xmm4\n"
						"addps		xmm3, xmm4\n"
						STORE4( 0, xmm2, xmm5)
						STORE2LO( 16, xmm3, xmm6 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm0, [esi+0*4]\n"
						"movss		xmm1, [esi+2*4]\n"
						"movlps		xmm3, [edi+(0*6+0)*4]\n"
						"movhps		xmm3, [edi+(0*6+2)*4]\n"
						"movaps		xmm4, xmm0\n"
						"shufps		xmm4, xmm4, 0x00\n"
						"mulps		xmm3, xmm4\n"
						"movlps		xmm5, [edi+(1*6+0)*4]\n"
						"movhps		xmm5, [edi+(1*6+2)*4]\n"
						"movaps		xmm6, xmm0\n"
						"shufps		xmm6, xmm6, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						"movlps		xmm4, [edi+(2*6+0)*4]\n"
						"movhps		xmm4, [edi+(2*6+2)*4]\n"
						"shufps		xmm1, xmm1, 0x00\n"
						"mulps		xmm4, xmm1\n"
						"addps		xmm3, xmm4\n"
						STORE4( 0, xmm3, xmm7 )
						"shufps		xmm0, xmm0, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"movlps		xmm3, [edi+(0*6+4)*4]\n"
						"movhps		xmm3, [edi+(1*6+4)*4]\n"
						"mulps		xmm3, xmm0\n"
						"movhlps		xmm4, xmm3\n"
						"addps		xmm3, xmm4\n"
						"movlps		xmm5, [edi+(2*6+4)*4]\n"
						"mulps		xmm5, xmm1\n"
						"addps		xmm3, xmm5\n"
						STORE2LO( 16, xmm3, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm0, [esi+0*4]\n"
						"movlps		xmm1, [esi+2*4]\n"
						"movaps		xmm3, xmm0\n"
						"shufps		xmm3, xmm3, 0x00\n"
						"mulps		xmm3, [edi+(0*6+0)*4]\n"
						"movlps		xmm5, [edi+(1*6+0)*4]\n"
						"movhps		xmm5, [edi+(1*6+2)*4]\n"
						"movaps		xmm6, xmm0\n"
						"shufps		xmm6, xmm6, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						"movlps		xmm4, [edi+(2*6+0)*4]\n"
						"movhps		xmm4, [edi+(2*6+2)*4]\n"
						"movaps		xmm6, xmm1\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm4, xmm6\n"
						"addps		xmm3, xmm4\n"
						"movlps		xmm5, [edi+(3*6+0)*4]\n"
						"movhps		xmm5, [edi+(3*6+2)*4]\n"
						"movaps		xmm6, xmm1\n"
						"shufps		xmm6, xmm6, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						STORE4( 0, xmm3, xmm7 )
						"shufps		xmm0, xmm0, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"shufps		xmm1, xmm1, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"movlps		xmm3, [edi+(0*6+4)*4]\n"
						"movhps		xmm3, [edi+(1*6+4)*4]\n"
						"mulps		xmm3, xmm0\n"
						"movlps		xmm4, [edi+(2*6+4)*4]\n"
						"movhps		xmm4, [edi+(3*6+4)*4]\n"
						"mulps		xmm4, xmm1\n"
						"addps		xmm3, xmm4\n"
						"movhlps		xmm4, xmm3\n"
						"addps		xmm3, xmm4\n"
						STORE2LO( 16, xmm3, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm0, [esi+0*4]\n"
						"movlps		xmm1, [esi+2*4]\n"
						"movss		xmm2, [esi+4*4]\n"
						"movaps		xmm3, xmm0\n"
						"shufps		xmm3, xmm3, 0x00\n"
						"mulps		xmm3, [edi+(0*6+0)*4]\n"
						"movlps		xmm5, [edi+(1*6+0)*4]\n"
						"movhps		xmm5, [edi+(1*6+2)*4]\n"
						"movaps		xmm6, xmm0\n"
						"shufps		xmm6, xmm6, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						"movaps		xmm6, xmm1\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm6, [edi+(2*6+0)*4]\n"
						"addps		xmm3, xmm6\n"
						"movlps		xmm5, [edi+(3*6+0)*4]\n"
						"movhps		xmm5, [edi+(3*6+2)*4]\n"
						"movaps		xmm6, xmm1\n"
						"shufps		xmm6, xmm6, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						"shufps		xmm2, xmm2, 0x00\n"
						"movaps		xmm4, xmm2\n"
						"mulps		xmm4, [edi+(4*6+0)*4]\n"
						"addps		xmm3, xmm4\n"
						STORE4( 0, xmm3, xmm7 )
						"shufps		xmm0, xmm0, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"shufps		xmm1, xmm1, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"movlps		xmm3, [edi+(0*6+4)*4]\n"
						"movhps		xmm3, [edi+(1*6+4)*4]\n"
						"mulps		xmm3, xmm0\n"
						"movlps		xmm4, [edi+(2*6+4)*4]\n"
						"movhps		xmm4, [edi+(3*6+4)*4]\n"
						"mulps		xmm4, xmm1\n"
						"addps		xmm3, xmm4\n"
						"movhlps		xmm4, xmm3\n"
						"addps		xmm3, xmm4\n"
						"movlps		xmm5, [edi+(4*6+4)*4]\n"
						"mulps		xmm5, xmm2\n"
						"addps		xmm3, xmm5\n"
						STORE2LO( 16, xmm3, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm0, [esi]\n"
						"movhps		xmm0, [esi+8]\n"
						"movlps		xmm1, [esi+16]\n"
						"mulps		xmm0, [edi]\n"
						"mulps		xmm1, [edi+16]\n"
						"shufps		xmm1, xmm0, 0xB4 \n" /*R_SHUFFLEPS( 0, 1, 3, 2 )=10110100 */
						"addps		xmm0, xmm1\n"
						"movhlps		xmm2, xmm0\n"
						"addss		xmm2, xmm0\n"
						"shufps		xmm0, xmm0, 0x01 \n"  /* R_SHUFFLEPS( 1, 0, 0, 0 )=00000001 */
						"addss		xmm2, xmm0\n"
						STORE1( 0, xmm2, xmm3 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 2: {		// 6x2 * 6x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm0, [esi+0*4]\n"
						"shufps		xmm0, xmm0, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"movaps		xmm6, [edi+0*4]\n"
						"mulps		xmm6, xmm0\n"
						"movlps		xmm1, [esi+2*4]\n"
						"shufps		xmm1, xmm1, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"movaps		xmm7, [edi+4*4]\n"
						"mulps		xmm7, xmm1\n"
						"addps		xmm6, xmm7\n"
						"movlps		xmm2, [esi+4*4]\n"
						"shufps		xmm2, xmm2, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"movaps		xmm7, [edi+8*4]\n"
						"mulps		xmm7, xmm2\n"
						"addps		xmm6, xmm7\n"
						"movhlps		xmm3, xmm6\n"
						"addps		xmm3, xmm6\n"
						STORE2LO( 0, xmm3, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 3: {		// 6x3 * 6x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movss		xmm0, [edi+(0*3+2)*4]\n"
						"movhps		xmm0, [edi+(0*3+0)*4]\n"
						"shufps		xmm0, xmm0, 0x36\n"  /*R_SHUFFLEPS( 2, 1, 3, 0 )=00110110*/
						"movss		xmm6, [esi+0*4]\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm6, xmm0\n"
						"movss		xmm1, [edi+(1*3+0)*4]\n"
						"movhps		xmm1, [edi+(1*3+1)*4]\n"
						"movss		xmm7, [esi+1*4]\n"
						"shufps		xmm7, xmm7, 0x00\n"
						"mulps		xmm7, xmm1\n"
						"addps		xmm6, xmm7\n"
						"movss		xmm2, [edi+(2*3+2)*4]\n"
						"movhps		xmm2, [edi+(2*3+0)*4]\n"
						"shufps		xmm2, xmm2, 0x36\n"  /*R_SHUFFLEPS( 2, 1, 3, 0 )=00110110*/
						"movss		xmm7, [esi+2*4]\n"
						"shufps		xmm7, xmm7, 0x00\n"
						"mulps		xmm7, xmm2\n"
						"addps		xmm6, xmm7\n"
						"movss		xmm3, [edi+(3*3+0)*4]\n"
						"movhps		xmm3, [edi+(3*3+1)*4]\n"
						"movss		xmm7, [esi+3*4]\n"
						"shufps		xmm7, xmm7, 0x00\n"
						"mulps		xmm7, xmm3\n"
						"addps		xmm6, xmm7\n"
						"movss		xmm4, [edi+(4*3+2)*4]\n"
						"movhps		xmm4, [edi+(4*3+0)*4]\n"
						"shufps		xmm4, xmm4, 0x36\n"  /*R_SHUFFLEPS( 2, 1, 3, 0 )=00110110*/
						"movss		xmm7, [esi+4*4]\n"
						"shufps		xmm7, xmm7, 0x00\n"
						"mulps		xmm7, xmm4\n"
						"addps		xmm6, xmm7\n"
						"movss		xmm5, [edi+(5*3+0)*4]\n"
						"movhps		xmm5, [edi+(5*3+1)*4]\n"
						"movss		xmm7, [esi+5*4]\n"
						"shufps		xmm7, xmm7, 0x00\n"
						"mulps		xmm7, xmm5\n"
						"addps		xmm6, xmm7\n"
						STORE1( 0, xmm6, xmm7 )
						STORE2HI( 4, xmm6, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 4: {		// 6x4 * 6x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm3, [edi+(0*4+0)*4]\n"
						"movhps		xmm3, [edi+(0*4+2)*4]\n"
						"movss		xmm4, [esi+0*4]\n"
						"shufps		xmm4, xmm4, 0x00\n"
						"mulps		xmm3, xmm4\n"
						"movlps		xmm5, [edi+(1*4+0)*4]\n"
						"movhps		xmm5, [edi+(1*4+2)*4]\n"
						"movss		xmm6, [esi+1*4]\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						"movlps		xmm4, [edi+(2*4+0)*4]\n"
						"movhps		xmm4, [edi+(2*4+2)*4]\n"
						"movss		xmm6, [esi+2*4]\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm4, xmm6\n"
						"addps		xmm3, xmm4\n"
						"movlps		xmm5, [edi+(3*4+0)*4]\n"
						"movhps		xmm5, [edi+(3*4+2)*4]\n"
						"movss		xmm6, [esi+3*4]\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						"movlps		xmm4, [edi+(4*4+0)*4]\n"
						"movhps		xmm4, [edi+(4*4+2)*4]\n"
						"movss		xmm6, [esi+4*4]\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm4, xmm6\n"
						"addps		xmm3, xmm4\n"
						"movlps		xmm5, [edi+(5*4+0)*4]\n"
						"movhps		xmm5, [edi+(5*4+2)*4]\n"
						"movss		xmm6, [esi+5*4]\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						STORE4( 0, xmm3, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 5: {		// 6x5 * 6x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm6, [edi+(0*5+0)*4]\n"
						"movhps		xmm6, [edi+(0*5+2)*4]\n"
						"movss		xmm0, [esi+0*4]\n"
						"shufps		xmm0, xmm0, 0x00\n"
						"mulps		xmm6, xmm0\n"
						"movlps		xmm7, [edi+(1*5+0)*4]\n"
						"movhps		xmm7, [edi+(1*5+2)*4]\n"
						"movss		xmm1, [esi+1*4]\n"
						"shufps		xmm1, xmm1, 0x00\n"
						"mulps		xmm7, xmm1\n"
						"addps		xmm6, xmm7\n"
						"movlps		xmm7, [edi+(2*5+0)*4]\n"
						"movhps		xmm7, [edi+(2*5+2)*4]\n"
						"movss		xmm2, [esi+2*4]\n"
						"shufps		xmm2, xmm2, 0x00\n"
						"mulps		xmm7, xmm2\n"
						"addps		xmm6, xmm7\n"
						"movlps		xmm7, [edi+(3*5+0)*4]\n"
						"movhps		xmm7, [edi+(3*5+2)*4]\n"
						"movss		xmm3, [esi+3*4]\n"
						"shufps		xmm3, xmm3, 0x00\n"
						"mulps		xmm7, xmm3\n"
						"addps		xmm6, xmm7\n"
						"movlps		xmm7, [edi+(4*5+0)*4]\n"
						"movhps		xmm7, [edi+(4*5+2)*4]\n"
						"movss		xmm4, [esi+4*4]\n"
						"shufps		xmm4, xmm4, 0x00\n"
						"mulps		xmm7, xmm4\n"
						"addps		xmm6, xmm7\n"
						"movlps		xmm7, [edi+(5*5+0)*4]\n"
						"movhps		xmm7, [edi+(5*5+2)*4]\n"
						"movss		xmm5, [esi+5*4]\n"
						"shufps		xmm5, xmm5, 0x00\n"
						"mulps		xmm7, xmm5\n"
						"addps		xmm6, xmm7\n"
						STORE4( 0, xmm6, xmm7 )
						"movss		xmm6, [edi+(0*5+4)*4]\n"
						"mulss		xmm6, xmm0\n"
						"movss		xmm7, [edi+(1*5+4)*4]\n"
						"mulss		xmm7, xmm1\n"
						"addss		xmm6, xmm7\n"
						"movss		xmm7, [edi+(2*5+4)*4]\n"
						"mulss		xmm7, xmm2\n"
						"addss		xmm6, xmm7\n"
						"movss		xmm7, [edi+(3*5+4)*4]\n"
						"mulss		xmm7, xmm3\n"
						"addss		xmm6, xmm7\n"
						"movss		xmm7, [edi+(4*5+4)*4]\n"
						"mulss		xmm7, xmm4\n"
						"addss		xmm6, xmm7\n"
						"movss		xmm7, [edi+(5*5+4)*4]\n"
						"mulss		xmm7, xmm5\n"
						"addss		xmm6, xmm7\n"
						STORE1( 16, xmm6, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 6: {		// 6x6 * 6x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm0, [esi+0*4]\n"
						"movlps		xmm1, [esi+2*4]\n"
						"movlps		xmm2, [esi+4*4]\n"
						"movaps		xmm3, xmm0\n"
						"shufps		xmm3, xmm3, 0x00\n"
						"mulps		xmm3, [edi+(0*6+0)*4]\n"
						"movlps		xmm5, [edi+(1*6+0)*4]\n"
						"movhps		xmm5, [edi+(1*6+2)*4]\n"
						"movaps		xmm6, xmm0\n"
						"shufps		xmm6, xmm6, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						"movaps		xmm6, xmm1\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm6, [edi+(2*6+0)*4]\n"
						"addps		xmm3, xmm6\n"
						"movaps		xmm6, xmm1\n"
						"shufps		xmm6, xmm6, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"movlps		xmm5, [edi+(3*6+0)*4]\n"
						"movhps		xmm5, [edi+(3*6+2)*4]\n"
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						"movaps		xmm6, xmm2\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm6, [edi+(4*6+0)*4]\n"
						"addps		xmm3, xmm6\n"
						"movaps		xmm6, xmm2\n"
						"shufps		xmm6, xmm6, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"movlps		xmm5, [edi+(5*6+0)*4]\n"
						"movhps		xmm5, [edi+(5*6+2)*4]\n"
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						STORE4( 0, xmm3, xmm7 )
						"shufps		xmm0, xmm0, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"shufps		xmm1, xmm1, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"shufps		xmm2, xmm2, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"movlps		xmm3, [edi+(0*6+4)*4]\n"
						"movhps		xmm3, [edi+(1*6+4)*4]\n"
						"mulps		xmm3, xmm0\n"
						"movlps		xmm4, [edi+(2*6+4)*4]\n"
						"movhps		xmm4, [edi+(3*6+4)*4]\n"
						"mulps		xmm4, xmm1\n"
						"addps		xmm3, xmm4\n"
						"movlps		xmm5, [edi+(4*6+4)*4]\n"
						"movhps		xmm5, [edi+(5*6+4)*4]\n"
						"mulps		xmm5, xmm2\n"
						"addps		xmm3, xmm5\n"
						"movhlps	xmm4, xmm3\n"
						"addps		xmm3, xmm4\n"
						STORE2LO( 16, xmm3, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
void CSIMD_SSE::matX_TransposeMultiply"subVecX

	optimizes the following matrix multiplications:

	Nx6 * Nx1
	6xN * 6x1

	with N in the range [1-6]
============
*/
void VPCALL CSIMD_SSE::matX_TransposeMultiplySubVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec ) {
#define STORE1( offset, reg1, reg2 )		\
	"movss		"STRINGIZE(reg2)", [eax+"STRINGIZE(offset)"]\n"		\
	"subss		"STRINGIZE(reg2)", "STRINGIZE(reg1)"\n"			\
	"movss		[eax+"STRINGIZE(offset)"], "STRINGIZE(reg2)"\n"
#define STORE2LO( offset, reg1, reg2 )		\
	"movlps	"STRINGIZE(reg2)", [eax+"STRINGIZE(offset)"]\n"		\
	"subps		"STRINGIZE(reg2)", "STRINGIZE(reg1)"\n"				\
	"movlps	[eax+"STRINGIZE(offset)"], "STRINGIZE(reg2)"\n"
#define STORE2HI( offset, reg1, reg2 )		\
	"movhps	"STRINGIZE(reg2)", [eax+"STRINGIZE(offset)"]\n"		\
	"subps		"STRINGIZE(reg2)", "STRINGIZE(reg1)"\n"				\
	"movhps	[eax+"STRINGIZE(offset)"], "STRINGIZE(reg2)"\n"
#define STORE4( offset, reg1, reg2 )		\
	"movlps	"STRINGIZE(reg2)", [eax+"STRINGIZE(offset)"]\n"		\
	"movhps	"STRINGIZE(reg2)", [eax+"STRINGIZE(offset)"+8]\n"	\
	"subps		"STRINGIZE(reg2)", "STRINGIZE(reg1)"\n"				\
	"movlps	[eax+"STRINGIZE(offset)"], "STRINGIZE(reg2)"\n"		\
	"movhps	[eax+"STRINGIZE(offset)"+8], "STRINGIZE(reg2)"\n"
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
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movss		xmm0, [esi]\n"
						"shufps		xmm0, xmm0, 0x00\n"
						"movaps		xmm1, xmm0\n"
						"mulps		xmm0, [edi]\n"
						"mulps		xmm1, [edi+16]\n"
						STORE4( 0, xmm0, xmm2)
						STORE2LO( 16, xmm1, xmm3 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm0, [esi]\n"
						"movaps		xmm1, xmm0\n"
						"shufps		xmm0, xmm0, 0x00\n"
						"shufps		xmm1, xmm1, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"movaps		xmm2, [edi]\n"
						"mulps		xmm2, xmm0\n"
						"movlps		xmm3, [edi+24]\n"
						"movhps		xmm3, [edi+32]\n"
						"mulps		xmm3, xmm1\n"
						"addps		xmm2, xmm3\n"
						"shufps		xmm0, xmm1, 0x00\n"
						"movlps		xmm4, [edi+16]\n"
						"movhps		xmm4, [edi+40]\n"
						"mulps		xmm4, xmm0\n"
						"movhlps		xmm3, xmm4\n"
						"addps		xmm3, xmm4\n"
						STORE4( 0, xmm2, xmm5 )
						STORE2LO( 16, xmm3, xmm6 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm0, [esi+0*4]\n"
						"movss		xmm1, [esi+2*4]\n"
						"movlps		xmm3, [edi+(0*6+0)*4]\n"
						"movhps		xmm3, [edi+(0*6+2)*4]\n"
						"movaps		xmm4, xmm0\n"
						"shufps		xmm4, xmm4, 0x00\n"
						"mulps		xmm3, xmm4\n"
						"movlps		xmm5, [edi+(1*6+0)*4]\n"
						"movhps		xmm5, [edi+(1*6+2)*4]\n"
						"movaps		xmm6, xmm0\n"
						"shufps		xmm6, xmm6, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						"movlps		xmm4, [edi+(2*6+0)*4]\n"
						"movhps		xmm4, [edi+(2*6+2)*4]\n"
						"shufps		xmm1, xmm1, 0x00\n"
						"mulps		xmm4, xmm1\n"
						"addps		xmm3, xmm4\n"
						STORE4( 0, xmm3, xmm7 )
						"shufps		xmm0, xmm0, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"movlps		xmm3, [edi+(0*6+4)*4]\n"
						"movhps		xmm3, [edi+(1*6+4)*4]\n"
						"mulps		xmm3, xmm0\n"
						"movhlps		xmm4, xmm3\n"
						"addps		xmm3, xmm4\n"
						"movlps		xmm5, [edi+(2*6+4)*4]\n"
						"mulps		xmm5, xmm1\n"
						"addps		xmm3, xmm5\n"
						STORE2LO( 16, xmm3, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm0, [esi+0*4]\n"
						"movlps		xmm1, [esi+2*4]\n"
						"movaps		xmm3, xmm0\n"
						"shufps		xmm3, xmm3, 0x00\n"
						"mulps		xmm3, [edi+(0*6+0)*4]\n"
						"movlps		xmm5, [edi+(1*6+0)*4]\n"
						"movhps		xmm5, [edi+(1*6+2)*4]\n"
						"movaps		xmm6, xmm0\n"
						"shufps		xmm6, xmm6, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						"movlps		xmm4, [edi+(2*6+0)*4]\n"
						"movhps		xmm4, [edi+(2*6+2)*4]\n"
						"movaps		xmm6, xmm1\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm4, xmm6\n"
						"addps		xmm3, xmm4\n"
						"movlps		xmm5, [edi+(3*6+0)*4]\n"
						"movhps		xmm5, [edi+(3*6+2)*4]\n"
						"movaps		xmm6, xmm1\n"
						"shufps		xmm6, xmm6, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						STORE4( 0, xmm3, xmm7 )
						"shufps		xmm0, xmm0, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"shufps		xmm1, xmm1, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"movlps		xmm3, [edi+(0*6+4)*4]\n"
						"movhps		xmm3, [edi+(1*6+4)*4]\n"
						"mulps		xmm3, xmm0\n"
						"movlps		xmm4, [edi+(2*6+4)*4]\n"
						"movhps		xmm4, [edi+(3*6+4)*4]\n"
						"mulps		xmm4, xmm1\n"
						"addps		xmm3, xmm4\n"
						"movhlps		xmm4, xmm3\n"
						"addps		xmm3, xmm4\n"
						STORE2LO( 16, xmm3, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm0, [esi+0*4]\n"
						"movlps		xmm1, [esi+2*4]\n"
						"movss		xmm2, [esi+4*4]\n"
						"movaps		xmm3, xmm0\n"
						"shufps		xmm3, xmm3, 0x00\n"
						"mulps		xmm3, [edi+(0*6+0)*4]\n"
						"movlps		xmm5, [edi+(1*6+0)*4]\n"
						"movhps		xmm5, [edi+(1*6+2)*4]\n"
						"movaps		xmm6, xmm0\n"
						"shufps		xmm6, xmm6, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						"movaps		xmm6, xmm1\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm6, [edi+(2*6+0)*4]\n"
						"addps		xmm3, xmm6\n"
						"movlps		xmm5, [edi+(3*6+0)*4]\n"
						"movhps		xmm5, [edi+(3*6+2)*4]\n"
						"movaps		xmm6, xmm1\n"
						"shufps		xmm6, xmm6, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						"shufps		xmm2, xmm2, 0x00\n"
						"movaps		xmm4, xmm2\n"
						"mulps		xmm4, [edi+(4*6+0)*4]\n"
						"addps		xmm3, xmm4\n"
						STORE4( 0, xmm3, xmm7 )
						"shufps		xmm0, xmm0, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"shufps		xmm1, xmm1, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"movlps		xmm3, [edi+(0*6+4)*4]\n"
						"movhps		xmm3, [edi+(1*6+4)*4]\n"
						"mulps		xmm3, xmm0\n"
						"movlps		xmm4, [edi+(2*6+4)*4]\n"
						"movhps		xmm4, [edi+(3*6+4)*4]\n"
						"mulps		xmm4, xmm1\n"
						"addps		xmm3, xmm4\n"
						"movhlps		xmm4, xmm3\n"
						"addps		xmm3, xmm4\n"
						"movlps		xmm5, [edi+(4*6+4)*4]\n"
						"mulps		xmm5, xmm2\n"
						"addps		xmm3, xmm5\n"
						STORE2LO( 16, xmm3, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm0, [esi]\n"
						"movhps		xmm0, [esi+8]\n"
						"movlps		xmm1, [esi+16]\n"
						"mulps		xmm0, [edi]\n"
						"mulps		xmm1, [edi+16]\n"
						"shufps		xmm1, xmm0, 0xB4 \n" /*R_SHUFFLEPS( 0, 1, 3, 2 )=10110100 */
						"addps		xmm0, xmm1\n"
						"movhlps		xmm2, xmm0\n"
						"addss		xmm2, xmm0\n"
						"shufps		xmm0, xmm0, 0x01 \n"  /* R_SHUFFLEPS( 1, 0, 0, 0 )=00000001 */
						"addss		xmm2, xmm0\n"
						STORE1( 0, xmm2, xmm3 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 2: {		// 6x2 * 6x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm0, [esi+0*4]\n"
						"shufps		xmm0, xmm0, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"movaps		xmm6, [edi+0*4]\n"
						"mulps		xmm6, xmm0\n"
						"movlps		xmm1, [esi+2*4]\n"
						"shufps		xmm1, xmm1, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"movaps		xmm7, [edi+4*4]\n"
						"mulps		xmm7, xmm1\n"
						"addps		xmm6, xmm7\n"
						"movlps		xmm2, [esi+4*4]\n"
						"shufps		xmm2, xmm2, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"movaps		xmm7, [edi+8*4]\n"
						"mulps		xmm7, xmm2\n"
						"addps		xmm6, xmm7\n"
						"movhlps		xmm3, xmm6\n"
						"addps		xmm3, xmm6\n"
						STORE2LO( 0, xmm3, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 3: {		// 6x3 * 6x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movss		xmm0, [edi+(0*3+2)*4]\n"
						"movhps		xmm0, [edi+(0*3+0)*4]\n"
						"shufps		xmm0, xmm0, 0x36\n"  /*R_SHUFFLEPS( 2, 1, 3, 0 )=00110110*/
						"movss		xmm6, [esi+0*4]\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm6, xmm0\n"
						"movss		xmm1, [edi+(1*3+0)*4]\n"
						"movhps		xmm1, [edi+(1*3+1)*4]\n"
						"movss		xmm7, [esi+1*4]\n"
						"shufps		xmm7, xmm7, 0x00\n"
						"mulps		xmm7, xmm1\n"
						"addps		xmm6, xmm7\n"
						"movss		xmm2, [edi+(2*3+2)*4]\n"
						"movhps		xmm2, [edi+(2*3+0)*4]\n"
						"shufps		xmm2, xmm2, 0x36\n"  /*R_SHUFFLEPS( 2, 1, 3, 0 )=00110110*/
						"movss		xmm7, [esi+2*4]\n"
						"shufps		xmm7, xmm7, 0x00\n"
						"mulps		xmm7, xmm2\n"
						"addps		xmm6, xmm7\n"
						"movss		xmm3, [edi+(3*3+0)*4]\n"
						"movhps		xmm3, [edi+(3*3+1)*4]\n"
						"movss		xmm7, [esi+3*4]\n"
						"shufps		xmm7, xmm7, 0x00\n"
						"mulps		xmm7, xmm3\n"
						"addps		xmm6, xmm7\n"
						"movss		xmm4, [edi+(4*3+2)*4]\n"
						"movhps		xmm4, [edi+(4*3+0)*4]\n"
						"shufps		xmm4, xmm4, 0x36\n"  /*R_SHUFFLEPS( 2, 1, 3, 0 )=00110110*/
						"movss		xmm7, [esi+4*4]\n"
						"shufps		xmm7, xmm7, 0x00\n"
						"mulps		xmm7, xmm4\n"
						"addps		xmm6, xmm7\n"
						"movss		xmm5, [edi+(5*3+0)*4]\n"
						"movhps		xmm5, [edi+(5*3+1)*4]\n"
						"movss		xmm7, [esi+5*4]\n"
						"shufps		xmm7, xmm7, 0x00\n"
						"mulps		xmm7, xmm5\n"
						"addps		xmm6, xmm7\n"
						STORE1( 0, xmm6, xmm7 )
						STORE2HI( 4, xmm6, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 4: {		// 6x4 * 6x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm3, [edi+(0*4+0)*4]\n"
						"movhps		xmm3, [edi+(0*4+2)*4]\n"
						"movss		xmm4, [esi+0*4]\n"
						"shufps		xmm4, xmm4, 0x00\n"
						"mulps		xmm3, xmm4\n"
						"movlps		xmm5, [edi+(1*4+0)*4]\n"
						"movhps		xmm5, [edi+(1*4+2)*4]\n"
						"movss		xmm6, [esi+1*4]\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						"movlps		xmm4, [edi+(2*4+0)*4]\n"
						"movhps		xmm4, [edi+(2*4+2)*4]\n"
						"movss		xmm6, [esi+2*4]\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm4, xmm6\n"
						"addps		xmm3, xmm4\n"
						"movlps		xmm5, [edi+(3*4+0)*4]\n"
						"movhps		xmm5, [edi+(3*4+2)*4]\n"
						"movss		xmm6, [esi+3*4]\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						"movlps		xmm4, [edi+(4*4+0)*4]\n"
						"movhps		xmm4, [edi+(4*4+2)*4]\n"
						"movss		xmm6, [esi+4*4]\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm4, xmm6\n"
						"addps		xmm3, xmm4\n"
						"movlps		xmm5, [edi+(5*4+0)*4]\n"
						"movhps		xmm5, [edi+(5*4+2)*4]\n"
						"movss		xmm6, [esi+5*4]\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						STORE4( 0, xmm3, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 5: {		// 6x5 * 6x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm6, [edi+(0*5+0)*4]\n"
						"movhps		xmm6, [edi+(0*5+2)*4]\n"
						"movss		xmm0, [esi+0*4]\n"
						"shufps		xmm0, xmm0, 0x00\n"
						"mulps		xmm6, xmm0\n"
						"movlps		xmm7, [edi+(1*5+0)*4]\n"
						"movhps		xmm7, [edi+(1*5+2)*4]\n"
						"movss		xmm1, [esi+1*4]\n"
						"shufps		xmm1, xmm1, 0x00\n"
						"mulps		xmm7, xmm1\n"
						"addps		xmm6, xmm7\n"
						"movlps		xmm7, [edi+(2*5+0)*4]\n"
						"movhps		xmm7, [edi+(2*5+2)*4]\n"
						"movss		xmm2, [esi+2*4]\n"
						"shufps		xmm2, xmm2, 0x00\n"
						"mulps		xmm7, xmm2\n"
						"addps		xmm6, xmm7\n"
						"movlps		xmm7, [edi+(3*5+0)*4]\n"
						"movhps		xmm7, [edi+(3*5+2)*4]\n"
						"movss		xmm3, [esi+3*4]\n"
						"shufps		xmm3, xmm3, 0x00\n"
						"mulps		xmm7, xmm3\n"
						"addps		xmm6, xmm7\n"
						"movlps		xmm7, [edi+(4*5+0)*4]\n"
						"movhps		xmm7, [edi+(4*5+2)*4]\n"
						"movss		xmm4, [esi+4*4]\n"
						"shufps		xmm4, xmm4, 0x00\n"
						"mulps		xmm7, xmm4\n"
						"addps		xmm6, xmm7\n"
						"movlps		xmm7, [edi+(5*5+0)*4]\n"
						"movhps		xmm7, [edi+(5*5+2)*4]\n"
						"movss		xmm5, [esi+5*4]\n"
						"shufps		xmm5, xmm5, 0x00\n"
						"mulps		xmm7, xmm5\n"
						"addps		xmm6, xmm7\n"
						STORE4( 0, xmm6, xmm7 )
						"movss		xmm6, [edi+(0*5+4)*4]\n"
						"mulss		xmm6, xmm0\n"
						"movss		xmm7, [edi+(1*5+4)*4]\n"
						"mulss		xmm7, xmm1\n"
						"addss		xmm6, xmm7\n"
						"movss		xmm7, [edi+(2*5+4)*4]\n"
						"mulss		xmm7, xmm2\n"
						"addss		xmm6, xmm7\n"
						"movss		xmm7, [edi+(3*5+4)*4]\n"
						"mulss		xmm7, xmm3\n"
						"addss		xmm6, xmm7\n"
						"movss		xmm7, [edi+(4*5+4)*4]\n"
						"mulss		xmm7, xmm4\n"
						"addss		xmm6, xmm7\n"
						"movss		xmm7, [edi+(5*5+4)*4]\n"
						"mulss		xmm7, xmm5\n"
						"addss		xmm6, xmm7\n"
						STORE1( 16, xmm6, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
					return;
				}
				case 6: {		// 6x6 * 6x1
					asm (
						"mov			esi, %0\n"
						"mov			edi, %1\n"
						"mov			eax, %2\n"
						"movlps		xmm0, [esi+0*4]\n"
						"movlps		xmm1, [esi+2*4]\n"
						"movlps		xmm2, [esi+4*4]\n"
						"movaps		xmm3, xmm0\n"
						"shufps		xmm3, xmm3, 0x00\n"
						"mulps		xmm3, [edi+(0*6+0)*4]\n"
						"movlps		xmm5, [edi+(1*6+0)*4]\n"
						"movhps		xmm5, [edi+(1*6+2)*4]\n"
						"movaps		xmm6, xmm0\n"
						"shufps		xmm6, xmm6, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						"movaps		xmm6, xmm1\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm6, [edi+(2*6+0)*4]\n"
						"addps		xmm3, xmm6\n"
						"movaps		xmm6, xmm1\n"
						"shufps		xmm6, xmm6, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"movlps		xmm5, [edi+(3*6+0)*4]\n"
						"movhps		xmm5, [edi+(3*6+2)*4]\n"
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						"movaps		xmm6, xmm2\n"
						"shufps		xmm6, xmm6, 0x00\n"
						"mulps		xmm6, [edi+(4*6+0)*4]\n"
						"addps		xmm3, xmm6\n"
						"movaps		xmm6, xmm2\n"
						"shufps		xmm6, xmm6, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
						"movlps		xmm5, [edi+(5*6+0)*4]\n"
						"movhps		xmm5, [edi+(5*6+2)*4]\n"
						"mulps		xmm5, xmm6\n"
						"addps		xmm3, xmm5\n"
						STORE4( 0, xmm3, xmm7 )
						"shufps		xmm0, xmm0, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"shufps		xmm1, xmm1, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"shufps		xmm2, xmm2, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
						"movlps		xmm3, [edi+(0*6+4)*4]\n"
						"movhps		xmm3, [edi+(1*6+4)*4]\n"
						"mulps		xmm3, xmm0\n"
						"movlps		xmm4, [edi+(2*6+4)*4]\n"
						"movhps		xmm4, [edi+(3*6+4)*4]\n"
						"mulps		xmm4, xmm1\n"
						"addps		xmm3, xmm4\n"
						"movlps		xmm5, [edi+(4*6+4)*4]\n"
						"movhps		xmm5, [edi+(5*6+4)*4]\n"
						"mulps		xmm5, xmm2\n"
						"addps		xmm3, xmm5\n"
						"movhlps		xmm4, xmm3\n"
						"addps		xmm3, xmm4\n"
						STORE2LO( 16, xmm3, xmm7 )
					::"r"(vPtr),"r"(mPtr),"r"(dstPtr)
                        :);
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
	FPU version. At times up to 40% less clock cycles on a P3. In practise however,
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
						asm (
							"mov			esi, %0\n"
							"mov			edi, %1\n"
							"mov			eax, %2\n"
							"movss		xmm0, [edi]\n"
							"shufps		xmm0, xmm0, 0x00\n"
							"movaps		xmm1, [esi]\n"
							"mulps		xmm1, xmm0\n"
							"movaps		[eax], xmm1\n"
							"movlps		xmm2, [esi+16]\n"
							"mulps		xmm2, xmm0\n"
							"movlps		[eax+16], xmm2\n"
						::"r"(m2Ptr),
							"r"(m1Ptr),
							"r"(dstPtr)
                        :"esi","edi","eax"
                            );
						return;
					}
					case 6: {			// 6x1 * 1x6, no precision loss compared to FPU version
						asm (
							"mov			esi, %0\n"
							"mov			edi, %1\n"
							"mov			eax, %2\n"
							"xorps		xmm1, xmm1\n"
							"movaps		xmm0, [edi]\n"
							"movlps		xmm1, [edi+16]\n"
							"movlhps		xmm1, xmm0\n"
							"movhlps		xmm2, xmm0\n"
							"movlhps		xmm2, xmm1\n"
							// row 0 "and 1
							"movaps		xmm3, [esi]\n"
							"movaps		xmm4, xmm3\n"
							"shufps		xmm4, xmm4, 0x00\n"
							"movaps		xmm5, xmm3\n"
							"shufps		xmm5, xmm5, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
							"movaps		xmm6, xmm3\n"
							"shufps		xmm6, xmm6, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
							"mulps		xmm4, xmm0\n"
							"mulps		xmm5, xmm1\n"
							"mulps		xmm6, xmm2\n"
							"movaps		[eax], xmm4\n"
							"movaps		[eax+16], xmm5\n"
							"movaps		[eax+32], xmm6\n"
							// row 2 "and 3
							"movaps		xmm4, xmm3\n"
							"shufps		xmm4, xmm4, 0xAA \n" /*R_SHUFFLEPS( 2, 2, 2, 2 )=10101010*/
							"movaps		xmm5, xmm3\n"
							"shufps		xmm5, xmm5, 0xFA\n" /*R_SHUFFLEPS( 2, 2, 3, 3 )=11111010*/
							"shufps		xmm3, xmm3, 0xFF\n" /*R_SHUFFLEPS( 3, 3, 3, 3 )=11111111*/
							"mulps		xmm4, xmm0\n"
							"mulps		xmm5, xmm1\n"
							"mulps		xmm3, xmm2\n"
							"movaps		[eax+48], xmm4\n"
							"movaps		[eax+64], xmm5\n"
							"movaps		[eax+80], xmm3\n"
							// row 4 "and 5
							"movlps		xmm3, [esi+16]\n"
							"movaps		xmm4, xmm3\n"
							"shufps		xmm4, xmm4, 0x00\n"
							"movaps		xmm5, xmm3\n"
							"shufps		xmm5, xmm5, 0x50\n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
							"shufps		xmm3, xmm3, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
							"mulps		xmm4, xmm0\n"
							"mulps		xmm5, xmm1\n"
							"mulps		xmm3, xmm2\n"
							"movaps		[eax+96], xmm4\n"
							"movaps		[eax+112], xmm5\n"
							"movaps		[eax+128], xmm3\n"
						::"r"(m2Ptr),
							"r"(m1Ptr),
							"r"(dstPtr)

                        :"esi","edi","eax");
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

						#define MUL_Nx2_2x6_INIT( A, B , C )								\
						"mov		esi, "STRINGIZE(A)"\n"								\
						"mov		edi, "STRINGIZE(B)"\n"								\
						"mov		eax, "STRINGIZE(C)"	\n"							\
						"movaps	xmm0, [esi]\n"								\
						"movlps	xmm1, [esi+16]\n"							\
						"movhps	xmm1, [esi+40]\n"							\
						"movlps	xmm2, [esi+24]\n"							\
						"movhps	xmm2, [esi+32]\n"

						#define MUL_Nx2_2x6_ROW2( row )							\
						"movaps	xmm3, [edi+"STRINGIZE(row)"*16]\n"						\
						"movaps	xmm5, xmm0\n"								\
						"movaps	xmm4, xmm3\n"								\
						"shufps	xmm4, xmm4, 0x00\n"	\
						"mulps		xmm5, xmm4\n"								\
						"movaps	xmm4, xmm3\n"								\
						"movaps	xmm6, xmm2\n"								\
						"shufps	xmm4, xmm4, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */	\
						"mulps		xmm6, xmm4\n"								\
						"addps		xmm5, xmm6\n"								\
						"movaps	[eax+"STRINGIZE(row)"*48], xmm5\n"						\
						"movaps	xmm4, xmm3\n"								\
						"shufps	xmm4, xmm4, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/	\
						"movaps	xmm7, xmm1\n"								\
						"mulps		xmm7, xmm4\n"								\
						"movaps	xmm4, xmm3\n"								\
						"movaps	xmm5, xmm0\n"								\
						"shufps	xmm4, xmm4, 0xAA \n" /*R_SHUFFLEPS( 2, 2, 2, 2 )=10101010*/	\
						"mulps		xmm5, xmm4\n"								\
						"movaps	xmm4, xmm3\n"								\
						"movaps	xmm6, xmm2\n"								\
						"shufps	xmm4, xmm4, 0xFF \n" /*R_SHUFFLEPS( 3, 3, 3, 3 )=11111111*/	\
						"mulps		xmm6, xmm4\n"								\
						"addps		xmm5, xmm6\n"								\
						"shufps	xmm3, xmm3, 0xFA\n" /*R_SHUFFLEPS( 2, 2, 3, 3 )=11111010*/	\
						"movaps	xmm6, xmm1\n"								\
						"mulps		xmm6, xmm3\n"								\
						"movaps	xmm4, xmm7\n"								\
						"movlhps	xmm7, xmm6\n"								\
						"movhlps	xmm6, xmm4\n"								\
						"addps		xmm6, xmm7\n"								\
						"movlps	[eax+"STRINGIZE(row)"*48+16], xmm6\n"					\
						"movlps	[eax+"STRINGIZE(row)"*48+24], xmm5\n"					\
						"movhps	[eax+"STRINGIZE(row)"*48+32], xmm5\n"					\
						"movhps	[eax+"STRINGIZE(row)"*48+40], xmm6\n"
                         asm(
						MUL_Nx2_2x6_INIT(%0,%1,%2)
						MUL_Nx2_2x6_ROW2( 0 )
                        ::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:"esi","edi","eax");
						return;
					}
					case 6: {			// 6x2 * 2x6
                        asm(
						MUL_Nx2_2x6_INIT(%0,%1,%2)
						MUL_Nx2_2x6_ROW2( 0 )
						MUL_Nx2_2x6_ROW2( 1 )
						MUL_Nx2_2x6_ROW2( 2 )
                        ::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:"esi","edi","eax");
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
						asm (
							"mov		esi, %0\n"
							"mov		edi, %1\n"
							"mov		eax, %2\n"
							"movaps	xmm5, xmmword ptr [esi]\n"
							"movlps	xmm6, qword ptr [esi+24]\n"
							"movhps	xmm6, qword ptr [esi+32]\n"
							"movaps	xmm7, xmmword ptr [esi+48]\n"
							"movss	xmm0, dword ptr [edi]\n"
							"shufps	xmm0, xmm0, 0\n"
							"mulps	xmm0, xmm5\n"
							"movss	xmm1, dword ptr [edi+4]\n"
							"shufps	xmm1, xmm1, 0\n"
							"mulps	xmm1, xmm6\n"
							"movss	xmm2, dword ptr [edi+8]\n"
							"shufps	xmm2, xmm2, 0\n"
							"mulps	xmm2, xmm7\n"
							"addps	xmm0, xmm1\n"
							"addps	xmm0, xmm2\n"
							"movaps	xmmword ptr [eax], xmm0\n"
							"movss	xmm3, dword ptr [edi+12]\n"
							"shufps	xmm3, xmm3, 0\n"
							"mulps	xmm3, xmm5\n"
							"movss	xmm4, dword ptr [edi+16]\n"
							"shufps	xmm4, xmm4, 0\n"
							"mulps	xmm4, xmm6\n"
							"movss	xmm0, dword ptr [edi+20]\n"
							"shufps	xmm0, xmm0, 0\n"
							"mulps	xmm0, xmm7\n"
							"addps	xmm3, xmm4\n"
							"addps	xmm0, xmm3\n"
							"movlps	qword ptr [eax+24], xmm0\n"
							"movhps	qword ptr [eax+32], xmm0\n"
							"movss	xmm1, dword ptr [edi+24]\n"
							"shufps	xmm1, xmm1, 0\n"
							"mulps	xmm1, xmm5\n"
							"movss	xmm2, dword ptr [edi+28]\n"
							"shufps	xmm2, xmm2, 0\n"
							"mulps	xmm2, xmm6\n"
							"movss	xmm3, dword ptr [edi+32]\n"
							"shufps	xmm3, xmm3, 0\n"
							"mulps	xmm3, xmm7\n"
							"addps	xmm1, xmm2\n"
							"addps	xmm1, xmm3\n"
							"movaps	xmmword ptr [eax+48], xmm1\n"
							"movlps	xmm5, qword ptr [esi+16]\n"
							"movlps	xmm6, qword ptr [esi+40]\n"
							"movlps	xmm7, qword ptr [esi+64]\n"
							"shufps	xmm5, xmm5, 0x44\n"
							"shufps	xmm6, xmm6, 0x44\n"
							"shufps	xmm7, xmm7, 0x44\n"
							"movaps	xmm3, xmmword ptr [edi]\n"
							"movlps	xmm4, qword ptr [edi+16]\n"
							"movaps	xmm0, xmm3\n"
							"shufps	xmm0, xmm0, 0xF0\n"
							"mulps	xmm0, xmm5\n"
							"movaps	xmm1, xmm3\n"
							"shufps	xmm1, xmm4, 0x05\n"
							"mulps	xmm1, xmm6\n"
							"shufps	xmm3, xmm4, 0x5A\n"
							"mulps	xmm3, xmm7\n"
							"addps	xmm1, xmm0\n"
							"addps	xmm1, xmm3\n"
							"movlps	qword ptr [eax+16], xmm1\n"
							"movhps	qword ptr [eax+40], xmm1\n"
							"movss	xmm0, dword ptr [edi+24]\n"
							"shufps	xmm0, xmm0, 0\n"
							"mulps	xmm0, xmm5\n"
							"movss	xmm2, dword ptr [edi+28]\n"
							"shufps	xmm2, xmm2, 0\n"
							"mulps	xmm2, xmm6\n"
							"movss	xmm4, dword ptr [edi+32]\n"
							"shufps	xmm4, xmm4, 0\n"
							"mulps	xmm4, xmm7\n"
							"addps	xmm0, xmm2\n"
							"addps	xmm0, xmm4\n"
							"movlps	qword ptr [eax+64], xmm0\n"
						::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:"esi","edi","eax");
						return;
					}
					case 6: {			// 6x3 * 3x6
						#define MUL_Nx3_3x6_FIRST4COLUMNS_INIT( A , B , C )						\
						"mov			esi, "STRINGIZE(A)"\n"								\
						"mov			edi, "STRINGIZE(B)"\n"								\
						"mov			eax, "STRINGIZE(C)"\n"							\
						"movlps		xmm0, [esi+ 0*4]\n"						\
						"movhps		xmm0, [esi+ 2*4]\n"						\
						"movlps		xmm1, [esi+ 6*4]\n"						\
						"movhps		xmm1, [esi+ 8*4]\n"						\
						"movlps		xmm2, [esi+12*4]\n"						\
						"movhps		xmm2, [esi+14*4]\n"

						#define MUL_Nx3_3x6_FIRST4COLUMNS_ROW( row )				\
						"movss			xmm3, [edi+("STRINGIZE(row)"*3+0)*4]\n"					\
						"shufps		xmm3, xmm3, 0x00\n"	\
						"mulps			xmm3, xmm0\n"								\
						"movss			xmm4, [edi+("STRINGIZE(row)"*3+1)*4]\n"					\
						"shufps		xmm4, xmm4, 0x00\n"	\
						"mulps			xmm4, xmm1\n"								\
						"addps			xmm3, xmm4\n"								\
						"movss			xmm5, [edi+("STRINGIZE(row)"*3+2)*4]\n"					\
						"shufps		xmm5, xmm5, 0x00\n"	\
						"mulps			xmm5, xmm2\n"								\
						"addps			xmm3, xmm5\n"								\
						"movlps		[eax+("STRINGIZE(row)"*6+0)*4], xmm3\n"					\
						"movhps		[eax+("STRINGIZE(row)"*6+2)*4], xmm3\n"

						#define MUL_Nx3_3x6_LAST2COLUMNS_ROW6						\
						"movlps		xmm0, [esi+ 4*4]\n"						\
						"movlps		xmm1, [esi+10*4]\n"						\
						"movlps		xmm2, [esi+16*4]\n"						\
						"shufps		xmm0, xmm0, 0x44\n"						\
						"shufps		xmm1, xmm1, 0x44\n"						\
						"shufps		xmm2, xmm2, 0x44\n"						\
						"movlps		xmm3, [edi+0*4]\n"							\
						"movhps		xmm3, [edi+2*4]\n"							\
						"movaps		xmm4, xmm3\n"								\
						"movaps		xmm5, xmm3\n"								\
						"shufps		xmm3, xmm3, 0xF0\n"						\
						"mulps			xmm3, xmm0	\n"							\
						"movlps		xmm6, [edi+4*4]\n"							\
						"movhps		xmm6, [edi+6*4]\n"							\
						"shufps		xmm4, xmm6, 0x05\n"						\
						"mulps			xmm4, xmm1\n"								\
						"addps			xmm3, xmm4\n"								\
						"shufps		xmm5, xmm6, 0x5A\n"						\
						"mulps			xmm5, xmm2\n"								\
						"addps			xmm3, xmm5\n"								\
						"movlps		[eax+4*4], xmm3\n"							\
						"movhps		[eax+10*4], xmm3\n"						\
						"movaps		xmm5, xmm6\n"								\
						"movlps		xmm3, [edi+8*4]\n"							\
						"movhps		xmm3, [edi+10*4]\n"						\
						"movaps		xmm4, xmm3\n"								\
						"shufps		xmm5, xmm3, 0x5A\n"						\
						"mulps			xmm5, xmm0\n"								\
						"shufps		xmm6, xmm3, 0xAF\n"						\
						"mulps			xmm6, xmm1\n"								\
						"addps			xmm5, xmm6\n"								\
						"shufps		xmm4, xmm4, 0xF0\n"						\
						"mulps			xmm4, xmm2\n"								\
						"addps			xmm4, xmm5\n"								\
						"movlps		[eax+16*4], xmm4\n"						\
						"movhps		[eax+22*4], xmm4\n"						\
						"movlps		xmm6, [edi+12*4]\n"						\
						"movhps		xmm6, [edi+14*4]\n"						\
						"movaps		xmm5, xmm6\n"								\
						"movaps		xmm4, xmm6\n"								\
						"shufps		xmm6, xmm6, 0xF0\n"						\
						"mulps			xmm6, xmm0\n"								\
						"movlps		xmm3, [edi+16*4]\n"						\
						"shufps		xmm5, xmm3, 0x05\n"						\
						"mulps			xmm5, xmm1\n"								\
						"addps			xmm5, xmm6\n"								\
						"shufps		xmm4, xmm3, 0x5A\n"						\
						"mulps			xmm4, xmm2\n"								\
						"addps			xmm4, xmm5\n"								\
						"movlps		[eax+28*4], xmm4\n"						\
						"movhps		[eax+34*4], xmm4\n"
                        asm(
						MUL_Nx3_3x6_FIRST4COLUMNS_INIT(%0,%1,%2)
						MUL_Nx3_3x6_FIRST4COLUMNS_ROW( 0 )
						MUL_Nx3_3x6_FIRST4COLUMNS_ROW( 1 )
						MUL_Nx3_3x6_FIRST4COLUMNS_ROW( 2 )
						MUL_Nx3_3x6_FIRST4COLUMNS_ROW( 3 )
						MUL_Nx3_3x6_FIRST4COLUMNS_ROW( 4 )
						MUL_Nx3_3x6_FIRST4COLUMNS_ROW( 5 )
						MUL_Nx3_3x6_LAST2COLUMNS_ROW6
                        ::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:"esi","edi","eax");
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

						#define MUL_Nx4_4x6_FIRST4COLUMNS_INIT( A , B , C )						\
						"mov			esi, "STRINGIZE(A)"\n"								\
						"mov			edi, "STRINGIZE(B)"\n"								\
						"mov			eax, "STRINGIZE(C)"	\n"							\
						"movlps		xmm0, [esi+ 0*4]\n"						\
						"movhps		xmm0, [esi+ 2*4]\n"						\
						"movlps		xmm1, [esi+ 6*4]\n"						\
						"movhps		xmm1, [esi+ 8*4]\n"						\
						"movlps		xmm2, [esi+12*4]\n"						\
						"movhps		xmm2, [esi+14*4]\n"						\
						"movlps		xmm3, [esi+18*4]\n"						\
						"movhps		xmm3, [esi+20*4]\n"

						#define MUL_Nx4_4x6_FIRST4COLUMNS_ROW( row )				\
						"movss			xmm4, [edi+"STRINGIZE(row)"*16+0*4]\n"					\
						"shufps		xmm4, xmm4, 0x00\n"	\
						"mulps			xmm4, xmm0	\n"							\
						"movss			xmm5, [edi+"STRINGIZE(row)"*16+1*4]\n"					\
						"shufps		xmm5, xmm5, 0x00\n"	\
						"mulps			xmm5, xmm1\n"								\
						"addps			xmm4, xmm5\n"								\
						"movss			xmm6, [edi+"STRINGIZE(row)"*16+2*4]\n"					\
						"shufps		xmm6, xmm6, 0x00\n"	\
						"mulps			xmm6, xmm2\n"								\
						"addps			xmm4, xmm6\n"								\
						"movss			xmm7, [edi+"STRINGIZE(row)"*16+3*4]\n"					\
						"shufps		xmm7, xmm7, 0x00\n"	\
						"mulps			xmm7, xmm3\n"								\
						"addps			xmm4, xmm7\n"								\
						"movlps		[eax+"STRINGIZE(row)"*24+0], xmm4\n"					\
						"movhps		[eax+"STRINGIZE(row)"*24+8], xmm4\n"

						#define MUL_Nx4_4x6_LAST2COLUMNS_INIT						\
						"movlps		xmm0, [esi+ 4*4]\n"						\
						"movlps		xmm1, [esi+10*4]\n"						\
						"movlps		xmm2, [esi+16*4]\n"						\
						"movlps		xmm3, [esi+22*4]\n"						\
						"shufps		xmm0, xmm1, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/	\
						"shufps		xmm1, xmm0, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/	\
						"shufps		xmm2, xmm3, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/	\
						"shufps		xmm3, xmm2, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/

						#define MUL_Nx4_4x6_LAST2COLUMNS_ROW2( row )				\
						"movlps		xmm7, [edi+"STRINGIZE(row)"*32+ 0*4]\n"					\
						"movhps		xmm7, [edi+"STRINGIZE(row)"*32+ 4*4]\n"					\
						"movaps		xmm6, xmm7\n"								\
						"shufps		xmm6, xmm6, 0xA0 \n"  /*R_SHUFFLEPS( 0, 0, 3, 3 )=10100000*/	\
						"mulps			xmm6, xmm0\n"								\
						"shufps		xmm7, xmm7, 0xA5 \n"  /*R_SHUFFLEPS( 1, 1, 2, 2 )=10100101	*/\
						"mulps			xmm7, xmm1\n"								\
						"addps			xmm6, xmm7\n"								\
						"movlps		xmm4, [edi+"STRINGIZE(row)"*32+ 2*4]\n"				\
						"movhps		xmm4, [edi+"STRINGIZE(row)"*32+ 6*4]\n"					\
						"movaps		xmm5, xmm4\n"								\
						"shufps		xmm5, xmm5, 0xA0 \n"  /*R_SHUFFLEPS( 0, 0, 3, 3 )=10100000*/	\
						"mulps			xmm5, xmm2\n"								\
						"addps			xmm6, xmm5\n"								\
						"shufps		xmm4, xmm4, 0xA5 \n"  /*R_SHUFFLEPS( 1, 1, 2, 2 )=10100101	*/	\
						"mulps			xmm4, xmm3\n"								\
						"addps			xmm6, xmm4\n"								\
						"movlps		[eax+"STRINGIZE(row)"*48+ 4*4], xmm6\n"					\
						"movhps		[eax+"STRINGIZE(row)"*48+10*4], xmm6\n"
                        asm(
						MUL_Nx4_4x6_FIRST4COLUMNS_INIT(%0 , %1, %2)
						MUL_Nx4_4x6_FIRST4COLUMNS_ROW( 0 )
						MUL_Nx4_4x6_FIRST4COLUMNS_ROW( 1 )
						MUL_Nx4_4x6_FIRST4COLUMNS_ROW( 2 )
						MUL_Nx4_4x6_FIRST4COLUMNS_ROW( 3 )
						MUL_Nx4_4x6_LAST2COLUMNS_INIT
						MUL_Nx4_4x6_LAST2COLUMNS_ROW2( 0 )
						MUL_Nx4_4x6_LAST2COLUMNS_ROW2( 1 )
                        ::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:"esi","edi","eax");
						return;
					}
					case 6: {			// 6x4 * 4x6
                        asm(
						MUL_Nx4_4x6_FIRST4COLUMNS_INIT(%0 , %1, %2)
						MUL_Nx4_4x6_FIRST4COLUMNS_ROW( 0 )
						MUL_Nx4_4x6_FIRST4COLUMNS_ROW( 1 )
						MUL_Nx4_4x6_FIRST4COLUMNS_ROW( 2 )
						MUL_Nx4_4x6_FIRST4COLUMNS_ROW( 3 )
						MUL_Nx4_4x6_FIRST4COLUMNS_ROW( 4 )
						MUL_Nx4_4x6_FIRST4COLUMNS_ROW( 5 )
						::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:);
						asm(
						MUL_Nx4_4x6_LAST2COLUMNS_INIT
						MUL_Nx4_4x6_LAST2COLUMNS_ROW2( 0 )
						MUL_Nx4_4x6_LAST2COLUMNS_ROW2( 1 )
						MUL_Nx4_4x6_LAST2COLUMNS_ROW2( 2 )
                        ::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:"esi","edi","eax");
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

						#define MUL_Nx5_5x6_FIRST4COLUMNS_INIT( A  , B  , C )						\
						"mov		esi, "STRINGIZE(A)"\n"								\
						"mov		edi, "STRINGIZE(B)"\n"								\
						"mov		eax, "STRINGIZE(C)"	\n"							\
						"movlps		xmm0, [esi+ 0*4]\n"						\
						"movhps		xmm0, [esi+ 2*4]\n"						\
						"movlps		xmm1, [esi+ 6*4]\n"						\
						"movhps		xmm1, [esi+ 8*4]\n"						\
						"movlps		xmm2, [esi+12*4]\n"						\
						"movhps		xmm2, [esi+14*4]\n"						\
						"movlps		xmm3, [esi+18*4]\n"						\
						"movhps		xmm3, [esi+20*4]\n"						\
						"movlps		xmm4, [esi+24*4]\n"						\
						"movhps		xmm4, [esi+26*4]\n"

						#define MUL_Nx5_5x6_FIRST4COLUMNS_ROW( row )				\
						"movss			xmm6, [edi+"STRINGIZE(row)"*20+0*4]\n"					\
						"shufps		xmm6, xmm6, 0x00\n"	\
						"mulps			xmm6, xmm0\n"								\
						"movss			xmm5, [edi+"STRINGIZE(row)"*20+1*4]\n"					\
						"shufps		xmm5, xmm5, 0x00\n"	\
						"mulps			xmm5, xmm1\n"								\
						"addps			xmm6, xmm5\n"								\
						"movss			xmm5, [edi+"STRINGIZE(row)"*20+2*4]\n"					\
						"shufps		xmm5, xmm5, 0x00\n"	\
						"mulps			xmm5, xmm2\n"								\
						"addps			xmm6, xmm5\n"								\
						"movss			xmm5, [edi+"STRINGIZE(row)"*20+3*4]\n"					\
						"shufps		xmm5, xmm5, 0x00\n"	\
						"mulps			xmm5, xmm3\n"								\
						"addps			xmm6, xmm5\n"								\
						"movss			xmm5, [edi+"STRINGIZE(row)"*20+4*4]\n"					\
						"shufps		xmm5, xmm5, 0x00\n"	\
						"mulps			xmm5, xmm4\n"								\
						"addps			xmm6, xmm5\n"								\
						"movlps		[eax+"STRINGIZE(row)"*24+0], xmm6\n"					\
						"movhps		[eax+"STRINGIZE(row)"*24+8], xmm6\n"

						#define MUL_Nx5_5x6_LAST2COLUMNS_INIT						\
						"movlps		xmm0, [esi+ 4*4]\n"						\
						"movlps		xmm1, [esi+10*4]\n"						\
						"movlps		xmm2, [esi+16*4]\n"						\
						"movlps		xmm3, [esi+22*4]\n"						\
						"movlps		xmm4, [esi+28*4]\n"						\
						"shufps		xmm0, xmm1, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/	\
						"shufps		xmm1, xmm2, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/	\
						"shufps		xmm2, xmm3, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/	\
						"shufps		xmm3, xmm4, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/	\
						"shufps		xmm4, xmm0, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/

						#define MUL_Nx5_5x6_LAST2COLUMNS_ROW2( row )				\
						"movlps		xmm7, [edi+"STRINGIZE(row)"*40+ 0*4]\n"					\
						"movhps		xmm7, [edi+"STRINGIZE(row)"*40+ 6*4]\n"					\
						"movaps		xmm6, xmm7\n"								\
						"shufps		xmm6, xmm6, 0xA0 \n"  /*R_SHUFFLEPS( 0, 0, 2, 2 )=10100000*/	\
						"mulps		xmm6, xmm0\n"								\
						"movaps		xmm5, xmm7\n"								\
						"shufps		xmm5, xmm5, 0xF5 \n" /*R_SHUFFLEPS( 1, 1, 3, 3 )=11110101*/	\
						"mulps		xmm5, xmm1\n"								\
						"addps		xmm6, xmm5\n"								\
						"movlps		xmm7, [edi+"STRINGIZE(row)"*40+ 2*4]\n"					\
						"movhps		xmm7, [edi+"STRINGIZE(row)"*40+ 8*4]\n"					\
						"movaps		xmm5, xmm7\n"								\
						"shufps		xmm5, xmm5, 0xA0 \n"  /*R_SHUFFLEPS( 0, 0, 2, 2 )=10100000*/	\
						"mulps		xmm5, xmm2\n"								\
						"addps		xmm6, xmm5\n"								\
						"movaps		xmm5, xmm7\n"								\
						"shufps		xmm5, xmm5, 0xF5 \n" /*R_SHUFFLEPS( 1, 1, 3, 3 )=11110101*/	\
						"mulps		xmm5, xmm3\n"								\
						"addps		xmm6, xmm5\n"								\
						"movlps		xmm5, [edi+"STRINGIZE(row)"*40+ 4*4]\n"					\
						"shufps		xmm5, xmm5, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/	\
						"mulps		xmm5, xmm4\n"								\
						"addps		xmm6, xmm5\n"								\
						"movlps		[eax+"STRINGIZE(row)"*48+ 4*4], xmm6\n"					\
						"movhps		[eax+"STRINGIZE(row)"*48+10*4], xmm6\n"

						#define MUL_Nx5_5x6_LAST2COLUMNS_ROW( row )					\
						"movlps		xmm6, [edi+20*4+0*4]\n"					\
						"unpcklps		xmm6, xmm6\n"								\
						"mulps			xmm6, xmm0\n"								\
						"movlps		xmm5, [edi+20*4+2*4]\n"					\
						"unpcklps		xmm5, xmm5\n"								\
						"mulps			xmm5, xmm2\n"								\
						"addps			xmm6, xmm5\n"								\
						"movss			xmm5, [edi+20*4+4*4]\n"					\
						"unpcklps		xmm5, xmm5\n"								\
						"mulps			xmm5, xmm4\n"								\
						"addps			xmm6, xmm5\n"								\
						"movhlps		xmm7, xmm6\n"								\
						"addps			xmm6, xmm7\n"								\
						"movlps		[eax+"STRINGIZE(row)"*24+4*4], xmm6\n"
                        asm volatile(
						MUL_Nx5_5x6_FIRST4COLUMNS_INIT(%0 , %1, %2)
						MUL_Nx5_5x6_FIRST4COLUMNS_ROW( 0 )
						MUL_Nx5_5x6_FIRST4COLUMNS_ROW( 1 )
						MUL_Nx5_5x6_FIRST4COLUMNS_ROW( 2 )
						MUL_Nx5_5x6_FIRST4COLUMNS_ROW( 3 )
						MUL_Nx5_5x6_FIRST4COLUMNS_ROW( 4 )

						MUL_Nx5_5x6_LAST2COLUMNS_INIT
						MUL_Nx5_5x6_LAST2COLUMNS_ROW2( 0 )
						MUL_Nx5_5x6_LAST2COLUMNS_ROW2( 1 )
						MUL_Nx5_5x6_LAST2COLUMNS_ROW( 4 )
                        ::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:"esi","edi","eax");
						return;
					}
					case 6: {			// 6x5 * 5x6
                        asm volatile(
						MUL_Nx5_5x6_FIRST4COLUMNS_INIT(%0 , %1, %2)
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
                        ::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:"esi","edi","eax");

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

						#define MUL_Nx6_6x2_INIT( A , B , C )								\
						"mov		esi, "STRINGIZE(A)"\n"								\
						"mov		edi, "STRINGIZE(B)"\n"								\
						"mov		eax, "STRINGIZE(C)"\n"								\
						"movaps	xmm0, [esi]\n"								\
						"movaps	xmm1, [esi+16]\n"							\
						"movaps	xmm2, [esi+32]\n"

						#define MUL_Nx6_6x2_ROW2( row )							\
						"movaps	xmm7, [edi+"STRINGIZE(row)"*48+0*4]\n"					\
						"movaps	xmm6, xmm7\n"								\
						"shufps	xmm7, xmm7, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/	\
						"mulps		xmm7, xmm0\n"								\
						"shufps	xmm6, xmm6, 0xFA \n"  /*R_SHUFFLEPS( 2, 2, 3, 3 )=11111010 */	\
						"mulps		xmm6, xmm1\n"								\
						"addps		xmm7, xmm6\n"								\
						"movaps	xmm6, [edi+"STRINGIZE(row)"*48+4*4]\n"					\
						"movaps	xmm5, xmm6\n"								\
						"shufps	xmm6, xmm6, 0x50\n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/	\
						"mulps		xmm6, xmm2\n"								\
						"addps		xmm7, xmm6\n"								\
						"shufps	xmm5, xmm5, 0xFA\n" /*R_SHUFFLEPS( 2, 2, 3, 3 )=11111010*/	\
						"mulps		xmm5, xmm0\n"								\
						"movaps	xmm6, [edi+"STRINGIZE(row)"*48+24+2*4]\n"				\
						"movaps	xmm4, xmm6\n"								\
						"shufps	xmm6, xmm6, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/	\
						"mulps		xmm6, xmm1\n"								\
						"addps		xmm5, xmm6\n"								\
						"shufps	xmm4, xmm4, 0xFA\n" /*R_SHUFFLEPS( 2, 2, 3, 3 )=11111010*/	\
						"mulps		xmm4, xmm2\n"								\
						"addps		xmm5, xmm4\n"								\
						"movaps	xmm4, xmm5\n"								\
						"movhlps	xmm5, xmm7\n"								\
						"movlhps	xmm7, xmm4\n"								\
						"addps		xmm7, xmm5\n"								\
						"movaps	[eax+"STRINGIZE(row)"*16], xmm7\n"
                        asm(
						MUL_Nx6_6x2_INIT(%0 , %1, %2)
						MUL_Nx6_6x2_ROW2( 0 )
                        ::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:"esi","edi","eax");
						return;
					}
					break;
				}
				case 3: {
					if ( !(l^3) ) {		// 3x6 * 6x3

						#define MUL_Nx6_6x3_INIT( A , B , C )								\
						"mov		esi, "STRINGIZE(A)"\n"								\
						"mov		edi, "STRINGIZE(B)"\n"								\
						"mov		eax, "STRINGIZE(C)"\n"								\
						"movss		xmm0, [esi+ 0*4]\n"						\
						"movhps	xmm0, [esi+ 1*4]\n"						\
						"movss		xmm1, [esi+ 3*4]\n"						\
						"movhps	xmm1, [esi+ 4*4]\n"						\
						"movss		xmm2, [esi+ 6*4]\n"						\
						"movhps	xmm2, [esi+ 7*4]\n"						\
						"movss		xmm3, [esi+ 9*4]\n"						\
						"movhps	xmm3, [esi+10*4]\n"						\
						"movss		xmm4, [esi+12*4]\n"						\
						"movhps	xmm4, [esi+13*4]\n"						\
						"movss		xmm5, [esi+15*4]\n"						\
						"movhps	xmm5, [esi+16*4]\n"

						#define MUL_Nx6_6x3_ROW( row )							\
						"movss		xmm7, [edi+"STRINGIZE(row)"*24+0]\n"					\
						"shufps	xmm7, xmm7, 0x00\n"	\
						"mulps		xmm7, xmm0\n"								\
						"movss		xmm6, [edi+"STRINGIZE(row)"*24+4]\n"					\
						"shufps	xmm6, xmm6, 0x00\n"	\
						"mulps		xmm6, xmm1\n"								\
						"addps		xmm7, xmm6\n"								\
						"movss		xmm6, [edi+"STRINGIZE(row)"*24+8]\n"					\
						"shufps	xmm6, xmm6, 0x00\n"	\
						"mulps		xmm6, xmm2\n"								\
						"addps		xmm7, xmm6\n"								\
						"movss		xmm6, [edi+"STRINGIZE(row)"*24+12]\n"					\
						"shufps	xmm6, xmm6, 0x00\n"	\
						"mulps		xmm6, xmm3\n"								\
						"addps		xmm7, xmm6\n"								\
						"movss		xmm6, [edi+"STRINGIZE(row)"*24+16]\n"					\
						"shufps	xmm6, xmm6, 0x00\n"	\
						"mulps		xmm6, xmm4\n"								\
						"addps		xmm7, xmm6\n"								\
						"movss		xmm6, [edi+"STRINGIZE(row)"*24+20]\n"					\
						"shufps	xmm6, xmm6, 0x00\n"	\
						"mulps		xmm6, xmm5\n"								\
						"addps		xmm7, xmm6\n"								\
						"movss		[eax+"STRINGIZE(row)"*12+0], xmm7\n"					\
						"movhps	[eax+"STRINGIZE(row)"*12+4], xmm7\n"
                        asm(
						MUL_Nx6_6x3_INIT(%0 , %1, %2)
						MUL_Nx6_6x3_ROW( 0 )
						MUL_Nx6_6x3_ROW( 1 )
						MUL_Nx6_6x3_ROW( 2 )
                        ::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:"esi","edi","eax");
						return;
					}
					break;
				}
				case 4: {
					if ( !(l^4) ) {		// 4x6 * 6x4

						#define MUL_Nx6_6x4_INIT( A , B , C )								\
						"mov		esi, "STRINGIZE(A)"\n"								\
						"mov		edi, "STRINGIZE(B)"\n"								\
						"mov		eax, "STRINGIZE(C)"\n"								\
						"movaps	xmm0, [esi]\n"								\
						"movaps	xmm1, [esi+16]\n"							\
						"movaps	xmm2, [esi+32]\n"							\
						"movaps	xmm3, [esi+48]\n"							\
						"movaps	xmm4, [esi+64]\n"							\
						"movaps	xmm5, [esi+80]\n"

						#define MUL_Nx6_6x4_ROW( row )							\
						"movss		xmm7, [edi+"STRINGIZE(row)"*24+0]\n"					\
						"shufps	xmm7, xmm7, 0x00\n"	\
						"mulps		xmm7, xmm0\n"								\
						"movss		xmm6, [edi+"STRINGIZE(row)"*24+4]\n"					\
						"shufps	xmm6, xmm6, 0x00\n"	\
						"mulps		xmm6, xmm1\n"								\
						"addps		xmm7, xmm6\n"								\
						"movss		xmm6, [edi+"STRINGIZE(row)"*24+8]\n"					\
						"shufps	xmm6, xmm6, 0x00\n"	\
						"mulps		xmm6, xmm2\n"								\
						"addps		xmm7, xmm6\n"								\
						"movss		xmm6, [edi+"STRINGIZE(row)"*24+12]\n"					\
						"shufps	xmm6, xmm6, 0x00\n"	\
						"mulps		xmm6, xmm3\n"								\
						"addps		xmm7, xmm6\n"								\
						"movss		xmm6, [edi+"STRINGIZE(row)"*24+16]\n"					\
						"shufps	xmm6, xmm6, 0x00\n"	\
						"mulps		xmm6, xmm4\n"								\
						"addps		xmm7, xmm6\n"								\
						"movss		xmm6, [edi+"STRINGIZE(row)"*24+20]\n"					\
						"shufps	xmm6, xmm6, 0x00\n"	\
						"mulps		xmm6, xmm5\n"								\
						"addps		xmm7, xmm6\n"								\
						"movaps	[eax+"STRINGIZE(row)"*16], xmm7\n"
                        asm(
						MUL_Nx6_6x4_INIT(%0 , %1, %2)
						MUL_Nx6_6x4_ROW( 0 )
						MUL_Nx6_6x4_ROW( 1 )
						MUL_Nx6_6x4_ROW( 2 )
						MUL_Nx6_6x4_ROW( 3 )
                        ::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:"esi","edi","eax");
						return;
					}
					break;
				}
				case 5: {
					if ( !(l^5) ) {		// 5x6 * 6x5

						#define MUL_Nx6_6x5_INIT( A , B , C )								\
						"mov		esi, "STRINGIZE(A)"\n"								\
						"mov		edi, "STRINGIZE(B)"\n"								\
						"mov		eax, "STRINGIZE(C)"\n"								\
						"movaps	xmm0, [esi]\n"								\
						"movlps	xmm1, [esi+20]\n"							\
						"movhps	xmm1, [esi+28]\n"							\
						"movlps	xmm2, [esi+40]\n"							\
						"movhps	xmm2, [esi+48]\n"							\
						"movlps	xmm3, [esi+60]\n"							\
						"movhps	xmm3, [esi+68]\n"							\
						"movaps	xmm4, [esi+80]\n"							\
						"movlps	xmm5, [esi+100]\n"							\
						"movhps	xmm5, [esi+108]\n"

						#define MUL_Nx6_6x5_ROW( row )							\
						"movss		xmm7, [edi+"STRINGIZE(row)"*24+0]\n"					\
						"shufps	xmm7, xmm7, 0x00\n"	\
						"mulps		xmm7, xmm0\n"								\
						"fld		dword ptr [edi+("STRINGIZE(row)"*6+0)*4]\n"				\
						"fmul		dword ptr [esi+(4+0*5)*4]\n"				\
						"movss		xmm6, [edi+"STRINGIZE(row)"*24+4]\n"					\
						"shufps	xmm6, xmm6, 0x00\n"	\
						"mulps		xmm6, xmm1\n"								\
						"addps		xmm7, xmm6\n"								\
						"fld		dword ptr [edi+("STRINGIZE(row)"*6+1)*4]\n"				\
						"fmul		dword ptr [esi+(4+1*5)*4]\n"				\
						"faddp		st(1),st\n"								\
						"movss		xmm6, [edi+"STRINGIZE(row)"*24+8]\n"					\
						"shufps	xmm6, xmm6, 0x00\n"	\
						"mulps		xmm6, xmm2\n"								\
						"addps		xmm7, xmm6\n"								\
						"fld		dword ptr [edi+("STRINGIZE(row)"*6+2)*4]\n"				\
						"fmul		dword ptr [esi+(4+2*5)*4]\n"				\
						"faddp		st(1),st\n"								\
						"movss		xmm6, [edi+"STRINGIZE(row)"*24+12]\n"					\
						"shufps	xmm6, xmm6, 0x00\n"	\
						"mulps		xmm6, xmm3\n"								\
						"addps		xmm7, xmm6\n"								\
						"fld		dword ptr [edi+("STRINGIZE(row)"*6+3)*4]\n"				\
						"fmul		dword ptr [esi+(4+3*5)*4]\n"				\
						"faddp		st(1),st\n"								\
						"movss		xmm6, [edi+"STRINGIZE(row)"*24+16]\n"					\
						"shufps	xmm6, xmm6, 0x00\n"	\
						"mulps		xmm6, xmm4\n"								\
						"addps		xmm7, xmm6\n"								\
						"fld		dword ptr [edi+("STRINGIZE(row)"*6+4)*4]\n"				\
						"fmul		dword ptr [esi+(4+4*5)*4]\n"				\
						"faddp		st(1),st\n"								\
						"movss		xmm6, [edi+"STRINGIZE(row)"*24+20]\n"					\
						"shufps	xmm6, xmm6, 0x00\n"	\
						"mulps		xmm6, xmm5\n"								\
						"addps		xmm7, xmm6\n"								\
						"fld		dword ptr [edi+("STRINGIZE(row)"*6+5)*4]\n"				\
						"fmul		dword ptr [esi+(4+5*5)*4]\n"				\
						"faddp		st(1),st\n"								\
						"fstp		dword ptr [eax+("STRINGIZE(row)"*5+4)*4]\n"				\
						"movlps	[eax+"STRINGIZE(row)"*20], xmm7\n"						\
						"movhps	[eax+"STRINGIZE(row)"*20+8], xmm7\n"

                        asm(
						MUL_Nx6_6x5_INIT(%0,%1,%2)
						MUL_Nx6_6x5_ROW( 0 )
						MUL_Nx6_6x5_ROW( 1 )
						MUL_Nx6_6x5_ROW( 2 )
						MUL_Nx6_6x5_ROW( 3 )
						MUL_Nx6_6x5_ROW( 4 )
                        ::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:"esi","edi","eax");
						return;
					}
					break;
				}
				case 6: {
					switch( l ) {
						case 1: {		// 6x6 * 6x1
							asm (
								"mov			esi, %0\n"
								"mov			edi, %1\n"
								"mov			eax, %2\n"
								"movlps		xmm7, qword ptr [esi]\n"
								"movlps		xmm6, qword ptr [esi+8]\n"
								"shufps		xmm7, xmm7, 0x44\n"
								"shufps		xmm6, xmm6, 0x44\n"
								"movlps		xmm0, qword ptr [edi    ]\n"
								"movhps		xmm0, qword ptr [edi+ 24]\n"
								"mulps		xmm0, xmm7\n"
								"movlps		xmm3, qword ptr [edi+  8]\n"
								"movhps		xmm3, qword ptr [edi+ 32]\n"
								"mulps		xmm3, xmm6\n"
								"movlps		xmm1, qword ptr [edi+ 48]\n"
								"movhps		xmm1, qword ptr [edi+ 72]\n"
								"mulps		xmm1, xmm7\n"
								"movlps		xmm2, qword ptr [edi+ 96]\n"
								"movhps		xmm2, qword ptr [edi+120]\n"
								"mulps		xmm2, xmm7\n"
								"movlps		xmm4, qword ptr [edi+ 56]\n"
								"movhps		xmm4, qword ptr [edi+ 80]\n"
								"movlps		xmm5, qword ptr [edi+104]\n"
								"movhps		xmm5, qword ptr [edi+128]\n"
								"mulps		xmm4, xmm6\n"
								"movlps		xmm7, qword ptr [esi+16]\n"
								"addps		xmm0, xmm3\n"
								"shufps		xmm7, xmm7, 0x44\n"
								"mulps		xmm5, xmm6\n"
								"addps		xmm1, xmm4\n"
								"movlps		xmm3, qword ptr [edi+ 16]\n"
								"movhps		xmm3, qword ptr [edi+ 40]\n"
								"addps		xmm2, xmm5\n"
								"movlps		xmm4, qword ptr [edi+ 64]\n"
								"movhps		xmm4, qword ptr [edi+ 88]\n"
								"mulps		xmm3, xmm7\n"
								"movlps		xmm5, qword ptr [edi+112]\n"
								"movhps		xmm5, qword ptr [edi+136]\n"
								"addps		xmm0, xmm3\n"
								"mulps		xmm4, xmm7\n"
								"mulps		xmm5, xmm7\n"
								"addps		xmm1, xmm4\n"
								"addps		xmm2, xmm5\n"
								"movaps		xmm6, xmm0\n"
								"shufps		xmm0, xmm1, 0x88\n"
								"shufps		xmm6, xmm1, 0xDD\n"
								"movaps		xmm7, xmm2\n"
								"shufps		xmm7, xmm2, 0x88\n"
								"shufps		xmm2, xmm2, 0xDD\n"
								"addps		xmm0, xmm6\n"
								"addps		xmm2, xmm7\n"
								"movlps		[eax], xmm0\n"
								"movhps		[eax+8], xmm0\n"
								"movlps		[eax+16], xmm2\n"
							::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:"esi","edi","eax");
							return;
						}
						case 2: {		// 6x6 * 6x2
                            asm(
							MUL_Nx6_6x2_INIT( %0 , %1 , %2 )
							MUL_Nx6_6x2_ROW2( 0 )
							MUL_Nx6_6x2_ROW2( 1 )
							MUL_Nx6_6x2_ROW2( 2 )
                            ::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:"esi","edi","eax");
							return;
						}
						case 3: {		// 6x6 * 6x3
                            asm(
							MUL_Nx6_6x3_INIT(%0,%1,%2)
							MUL_Nx6_6x3_ROW( 0 )
							MUL_Nx6_6x3_ROW( 1 )
							MUL_Nx6_6x3_ROW( 2 )
							MUL_Nx6_6x3_ROW( 3 )
							MUL_Nx6_6x3_ROW( 4 )
							MUL_Nx6_6x3_ROW( 5 )
                            ::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:"esi","edi","eax");
							return;
						}
						case 4: {		// 6x6 * 6x4
                            asm(
							MUL_Nx6_6x4_INIT(%0,%1,%2)
							MUL_Nx6_6x4_ROW( 0 )
							MUL_Nx6_6x4_ROW( 1 )
							MUL_Nx6_6x4_ROW( 2 )
							MUL_Nx6_6x4_ROW( 3 )
							MUL_Nx6_6x4_ROW( 4 )
							MUL_Nx6_6x4_ROW( 5 )
                            ::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:"esi","edi","eax");
							return;
						}
						case 5: {		// 6x6 * 6x5
                            asm(
							MUL_Nx6_6x5_INIT(%0,%1,%2)
							MUL_Nx6_6x5_ROW( 0 )
							MUL_Nx6_6x5_ROW( 1 )
							MUL_Nx6_6x5_ROW( 2 )
							MUL_Nx6_6x5_ROW( 3 )
							MUL_Nx6_6x5_ROW( 4 )
							MUL_Nx6_6x5_ROW( 5 )
                            ::"r"(m2Ptr),
                            "r"(m1Ptr),
                            "r"(dstPtr)
                            :"esi","edi","eax");
							return;
						}
						case 6: {		// 6x6 * 6x6
							asm (
								"mov		ecx,  %0\n"
								"movlps		xmm3, qword ptr [ecx+72]\n"
								"mov		edx,%1\n"
								// Loading first 4 columns (upper 4 rows) of m2Ptr.
								"movaps		xmm0, xmmword ptr [ecx]\n"
								"movlps		xmm1, qword ptr [ecx+24]\n"
								"movhps		xmm1, qword ptr [ecx+32]\n"
								"movaps		xmm2, xmmword ptr [ecx+48]\n"
								"movhps		xmm3, qword ptr [ecx+80]\n"
								// Calculating first 4 elements in the first row of the destination matrix.
								"movss		xmm4, dword ptr [edx]\n"
								"movss		xmm5, dword ptr [edx+4]\n"
								"mov		eax, %2\n"
								"shufps		xmm4, xmm4, 0\n"
								"movss		xmm6, dword ptr [edx+8]\n"
								"shufps		xmm5, xmm5, 0\n"
								"movss		xmm7, dword ptr [edx+12]\n"
								"mulps		xmm4, xmm0\n"
								"shufps		xmm6, xmm6, 0\n"
								"shufps		xmm7, xmm7, 0\n"
								"mulps		xmm5, xmm1\n"
								"mulps		xmm6, xmm2\n"
								"addps		xmm5, xmm4\n"
								"mulps		xmm7, xmm3\n"
								"addps		xmm6, xmm5\n"
								"addps		xmm7, xmm6\n"
								"movaps		xmmword ptr [eax], xmm7\n"
								// Calculating first 4 elements in the second row of the destination matrix.
								"movss		xmm4, dword ptr [edx+24]\n"
								"shufps		xmm4, xmm4, 0\n"
								"mulps		xmm4, xmm0\n"
								"movss		xmm5, dword ptr [edx+28]\n"
								"shufps		xmm5, xmm5, 0\n"
								"mulps		xmm5, xmm1\n"
								"movss		xmm6, dword ptr [edx+32]\n"
								"shufps		xmm6, xmm6, 0\n"
								"movss		xmm7, dword ptr [edx+36]\n"
								"shufps		xmm7, xmm7, 0\n"
								"mulps		xmm6, xmm2\n"
								"mulps		xmm7, xmm3\n"
								"addps		xmm7, xmm6\n"
								"addps		xmm5, xmm4\n"
								"addps		xmm7, xmm5\n"
								// Calculating first 4 elements in the third row of the destination matrix.
								"movss		xmm4, dword ptr [edx+48]\n"
								"movss		xmm5, dword ptr [edx+52]\n"
								"movlps		qword ptr [eax+24], xmm7\n"// ; save 2nd
								"movhps		qword ptr [eax+32], xmm7\n" //; row
								"movss		xmm6, dword ptr [edx+56]\n"
								"movss		xmm7, dword ptr [edx+60]\n"
								"shufps		xmm4, xmm4, 0\n"
								"shufps		xmm5, xmm5, 0\n"
								"shufps		xmm6, xmm6, 0\n"
								"shufps		xmm7, xmm7, 0\n"
								"mulps		xmm4, xmm0\n"
								"mulps		xmm5, xmm1\n"
								"mulps		xmm6, xmm2\n"
								"mulps		xmm7, xmm3\n"
								"addps		xmm5, xmm4\n"
								"addps		xmm7, xmm6\n"
								"addps		xmm7, xmm5\n"
								"movaps		xmmword ptr [eax+48], xmm7\n"
								// Calculating first 4 elements in the fourth row of the destination matrix.
								"movss		xmm4, dword ptr [edx+72]\n"
								"movss		xmm5, dword ptr [edx+76]\n"
								"movss		xmm6, dword ptr [edx+80]\n"
								"movss		xmm7, dword ptr [edx+84]\n"
								"shufps		xmm4, xmm4, 0\n"
								"shufps		xmm5, xmm5, 0\n"
								"shufps		xmm6, xmm6, 0\n"
								"shufps		xmm7, xmm7, 0\n"
								"mulps		xmm4, xmm0\n"
								"mulps		xmm5, xmm1\n"
								"mulps		xmm6, xmm2\n"
								"mulps		xmm7, xmm3\n"
								"addps		xmm4, xmm5\n"
								"addps		xmm6, xmm4\n"
								"addps		xmm7, xmm6\n"
								"movlps		qword ptr [eax+72], xmm7\n"
								"movhps		qword ptr [eax+80], xmm7\n"
								// Calculating first 4 elements in the fifth row of the destination matrix.
								"movss		xmm4, dword ptr [edx+96]\n"
								"movss		xmm5, dword ptr [edx+100]\n"
								"movss		xmm6, dword ptr [edx+104]\n"
								"movss		xmm7, dword ptr [edx+108]\n"
								"shufps		xmm4, xmm4, 0\n"
								"shufps		xmm5, xmm5, 0\n"
								"shufps		xmm6, xmm6, 0\n"
								"shufps		xmm7, xmm7, 0\n"
								"mulps		xmm4, xmm0\n"
								"mulps		xmm5, xmm1\n"
								"mulps		xmm6, xmm2\n"
								"mulps		xmm7, xmm3\n"
								"addps		xmm5, xmm4\n"
								"addps		xmm7, xmm6\n"
								"addps		xmm7, xmm5\n"
								"movaps		xmmword ptr [eax+96], xmm7\n"
								// Calculating first 4 elements in the sixth row of the destination matrix.
								"movss		xmm4, dword ptr [edx+120]\n"
								"movss		xmm5, dword ptr [edx+124]\n"
								"movss		xmm6, dword ptr [edx+128]\n"
								"movss		xmm7, dword ptr [edx+132]\n"
								"shufps		xmm4, xmm4, 0\n"
								"shufps		xmm5, xmm5, 0\n"
								"shufps		xmm6, xmm6, 0\n"
								"shufps		xmm7, xmm7, 0\n"
								"mulps		xmm4, xmm0\n"
								"mulps		xmm5, xmm1\n"
								"mulps		xmm6, xmm2\n"
								"mulps		xmm7, xmm3\n"
								"addps		xmm4, xmm5\n"
								"addps		xmm6, xmm4\n"
								"addps		xmm7, xmm6\n"
								"movhps		qword ptr [eax+128], xmm7\n"
								"movlps		qword ptr [eax+120], xmm7\n"
								// Loading first 4 columns (lower 2 rows) of m2Ptr.
								"movlps		xmm0, qword ptr [ecx+96]\n"
								"movhps		xmm0, qword ptr [ecx+104]\n"
								"movlps		xmm1, qword ptr [ecx+120]\n"
								"movhps		xmm1, qword ptr [ecx+128]\n"
								// Calculating first 4 elements in the first row of the destination matrix.
								"movss		xmm2, dword ptr [edx+16]\n"
								"shufps		xmm2, xmm2, 0\n"
								"movss		xmm4, dword ptr [edx+40]\n"
								"movss		xmm3, dword ptr [edx+20]\n"
								"movss		xmm5, dword ptr [edx+44]\n"
								"movaps		xmm6, xmmword ptr [eax]\n"
								"movlps		xmm7, qword ptr [eax+24]\n"
								"shufps		xmm3, xmm3, 0\n"
								"shufps		xmm5, xmm5, 0\n"
								"movhps		xmm7, qword ptr [eax+32]\n"
								"shufps		xmm4, xmm4, 0\n"
								"mulps		xmm5, xmm1\n"
								"mulps		xmm2, xmm0\n"
								"mulps		xmm3, xmm1\n"
								"mulps		xmm4, xmm0\n"
								"addps		xmm6, xmm2\n"
								"addps		xmm7, xmm4\n"
								"addps		xmm7, xmm5\n"
								"addps		xmm6, xmm3\n"
								"movlps		qword ptr [eax+24], xmm7\n"
								"movaps		xmmword ptr [eax], xmm6\n"
								"movhps		qword ptr [eax+32], xmm7\n"
								// Calculating first 4 elements in the third row of the destination matrix.
								"movss		xmm2, dword ptr [edx+64]\n"
								"movss		xmm4, dword ptr [edx+88]\n"
								"movss		xmm5, dword ptr [edx+92]\n"
								"movss		xmm3, dword ptr [edx+68]\n"
								"movaps		xmm6, xmmword ptr [eax+48]\n"
								"movlps		xmm7, qword ptr [eax+72]\n"
								"movhps		xmm7, qword ptr [eax+80]\n"
								"shufps		xmm2, xmm2, 0\n"
								"shufps		xmm4, xmm4, 0\n"
								"shufps		xmm5, xmm5, 0\n"
								"shufps		xmm3, xmm3, 0\n"
								"mulps		xmm2, xmm0\n"
								"mulps		xmm4, xmm0\n"
								"mulps		xmm5, xmm1\n"
								"mulps		xmm3, xmm1\n"
								"addps		xmm6, xmm2\n"
								"addps		xmm6, xmm3\n"
								"addps		xmm7, xmm4\n"
								"addps		xmm7, xmm5\n"
								"movlps		qword ptr [eax+72], xmm7\n"
								"movaps		xmmword ptr [eax+48], xmm6\n"
								"movhps		qword ptr [eax+80], xmm7\n"
								// Calculating first 4 elements in the fifth row of the destination matrix.
								"movss		xmm2, dword ptr [edx+112]\n"
								"movss		xmm3, dword ptr [edx+116]\n"
								"movaps		xmm6, xmmword ptr [eax+96]\n"
								"shufps		xmm2, xmm2, 0\n"
								"shufps		xmm3, xmm3, 0\n"
								"mulps		xmm2, xmm0\n"
								"mulps		xmm3, xmm1\n"
								"addps		xmm6, xmm2\n"
								"addps		xmm6, xmm3\n"
								"movaps		xmmword ptr [eax+96], xmm6\n"
								// Calculating first 4 elements in the sixth row of the destination matrix.
								"movss		xmm4, dword ptr [edx+136]\n"
								"movss		xmm5, dword ptr [edx+140]\n"
								"movhps		xmm7, qword ptr [eax+128]\n"
								"movlps		xmm7, qword ptr [eax+120]\n"
								"shufps		xmm4, xmm4, 0\n"
								"shufps		xmm5, xmm5, 0\n"
								"mulps		xmm4, xmm0\n"
								"mulps		xmm5, xmm1\n"
								"addps		xmm7, xmm4\n"
								"addps		xmm7, xmm5\n"
								// Calculating last 2 columns of the destination matrix.
								"movlps		xmm0, qword ptr [ecx+16]\n"
								"movhps		xmm0, qword ptr [ecx+40]\n"
								"movhps		qword ptr [eax+128], xmm7\n"
								"movlps		qword ptr [eax+120], xmm7\n"
								"movlps		xmm2, qword ptr [ecx+64]\n"
								"movhps		xmm2, qword ptr [ecx+88]\n"
								"movaps		xmm3, xmm2\n"
								"shufps		xmm3, xmm3, 0x4E\n"
								"movlps		xmm4, qword ptr [ecx+112]\n"
								"movhps		xmm4, qword ptr [ecx+136]\n"
								"movaps		xmm5, xmm4\n"
								"shufps		xmm5, xmm5, 0x4E\n"
								"movlps		xmm6, qword ptr [edx]\n"
								"movhps		xmm6, qword ptr [edx+24]\n"
								"movaps		xmm7, xmm6\n"
								"shufps		xmm7, xmm7, 0x0F0\n"
								"mulps		xmm7, xmm0\n"
								"shufps		xmm6, xmm6, 0x0A5\n"
								"movaps		xmm1, xmm0\n"
								"shufps		xmm1, xmm1, 0x4E\n"
								"mulps		xmm1, xmm6\n"
								"addps		xmm7, xmm1\n"
								"movlps		xmm6, qword ptr [edx+8]\n"
								"movhps		xmm6, qword ptr [edx+32]\n"
								"movaps		xmm1, xmm6\n"
								"shufps		xmm1, xmm1, 0x0F0\n"
								"shufps		xmm6, xmm6, 0x0A5\n"
								"mulps		xmm1, xmm2\n"
								"mulps		xmm6, xmm3\n"
								"addps		xmm7, xmm1\n"
								"addps		xmm7, xmm6\n"
								"movhps		xmm6, qword ptr [edx+40]\n"
								"movlps		xmm6, qword ptr [edx+16]\n"
								"movaps		xmm1, xmm6\n"
								"shufps		xmm1, xmm1, 0x0F0\n"
								"shufps		xmm6, xmm6, 0x0A5\n"
								"mulps		xmm1, xmm4\n"
								"mulps		xmm6, xmm5\n"
								"addps		xmm7, xmm1\n"
								"addps		xmm7, xmm6\n"
								"movlps		qword ptr [eax+16], xmm7\n"
								"movhps		qword ptr [eax+40], xmm7\n"
								"movlps		xmm6, qword ptr [edx+48]\n"
								"movhps		xmm6, qword ptr [edx+72]\n"
								"movaps		xmm7, xmm6\n"
								"shufps		xmm7, xmm7, 0x0F0\n"
								"mulps		xmm7, xmm0\n"
								"shufps		xmm6, xmm6, 0x0A5\n"
								"movaps		xmm1, xmm0\n"
								"shufps		xmm1, xmm1, 0x4E\n"
								"mulps		xmm1, xmm6\n"
								"addps		xmm7, xmm1\n"
								"movhps		xmm6, qword ptr [edx+80]\n"
								"movlps		xmm6, qword ptr [edx+56]\n"
								"movaps		xmm1, xmm6\n"
								"shufps		xmm1, xmm1, 0x0F0\n"
								"shufps		xmm6, xmm6, 0x0A5\n"
								"mulps		xmm1, xmm2\n"
								"mulps		xmm6, xmm3\n"
								"addps		xmm7, xmm1\n"
								"addps		xmm7, xmm6\n"
								"movlps		xmm6, qword ptr [edx+64]\n"
								"movhps		xmm6, qword ptr [edx+88]\n"
								"movaps		xmm1, xmm6\n"
								"shufps		xmm1, xmm1, 0x0F0\n"
								"shufps		xmm6, xmm6, 0x0A5\n"
								"mulps		xmm1, xmm4\n"
								"mulps		xmm6, xmm5\n"
								"addps		xmm7, xmm1\n"
								"addps		xmm7, xmm6\n"
								"movlps		qword ptr [eax+64], xmm7\n"
								"movhps		qword ptr [eax+88], xmm7\n"
								"movlps		xmm6, qword ptr [edx+96]\n"
								"movhps		xmm6, qword ptr [edx+120]\n"
								"movaps		xmm7, xmm6\n"
								"shufps		xmm7, xmm7, 0x0F0\n"
								"mulps		xmm7, xmm0\n"
								"shufps		xmm6, xmm6, 0x0A5\n"
								"movaps		xmm1, xmm0\n"
								"shufps		xmm1, xmm1, 0x4E\n"
								"mulps		xmm1, xmm6\n"
								"addps		xmm7, xmm1\n"
								"movlps		xmm6, qword ptr [edx+104]\n"
								"movhps		xmm6, qword ptr [edx+128]\n"
								"movaps		xmm1, xmm6\n"
								"shufps		xmm1, xmm1, 0x0F0\n"
								"shufps		xmm6, xmm6, 0x0A5\n"
								"mulps		xmm1, xmm2\n"
								"mulps		xmm6, xmm3\n"
								"addps		xmm7, xmm1\n"
								"addps		xmm7, xmm6\n"
								"movlps		xmm6, qword ptr [edx+112]\n"
								"movhps		xmm6, qword ptr [edx+136]\n"
								"movaps		xmm1, xmm6\n"
								"shufps		xmm1, xmm1, 0x0F0\n"
								"shufps		xmm6, xmm6, 0x0A5\n"
								"mulps		xmm1, xmm4\n"
								"mulps		xmm6, xmm5\n"
								"addps		xmm7, xmm1\n"
								"addps		xmm7, xmm6\n"
								"movlps		qword ptr [eax+112], xmm7\n"
								"movhps		qword ptr [eax+136], xmm7\n"
							::"r"(m2Ptr),
                            "r"(m1Ptr),
                            "r"(dstPtr)
                            :"esi","edi","eax");
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
				asm (
					"mov		esi, %0\n"
					"mov		edi, %1\n"
					"mov		eax, %2\n"
					"movss	xmm0, [esi]\n"
					"shufps	xmm0, xmm0, 0x00\n"
					"movaps	xmm1, xmm0\n"
					"mulps	xmm0, [edi]\n"
					"mulps	xmm1, [edi+16]\n"
					"movaps	[eax], xmm0\n"
					"movlps	[eax+16], xmm1\n"
					::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:);

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
				#define MUL_2xN_2x2_INIT( A , B , C )								\
				"mov		esi, "STRINGIZE(A)"\n"								\
				"mov		edi, "STRINGIZE(B)"\n"								\
				"mov		eax, "STRINGIZE(C)"\n"								\
				"movlps	xmm0, [esi]\n"								\
				"shufps	xmm0, xmm0, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/	\
				"movlps	xmm1, [esi+8]\n"							\
				"shufps	xmm1, xmm1, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/

				#define MUL_2xN_2x2_ROW2( N, row )						\
				"movlps	xmm6, [edi+("STRINGIZE(row)"+0*"STRINGIZE(N)")*4]\n"					\
				"shufps	xmm6, xmm6, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/	\
				"movlps	xmm7, [edi+("STRINGIZE(row)"+1*"STRINGIZE(N)")*4]\n"					\
				"shufps	xmm7, xmm7, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/	\
				"mulps		xmm6, xmm0\n"								\
				"mulps		xmm7, xmm1\n"								\
				"addps		xmm6, xmm7\n"								\
				"movaps	[eax+("STRINGIZE(row)"*2)*4], xmm6\n"

                asm(
				MUL_2xN_2x2_INIT(%0 , %1, %2)
				MUL_2xN_2x2_ROW2( 6, 0 )
				MUL_2xN_2x2_ROW2( 6, 2 )
				MUL_2xN_2x2_ROW2( 6, 4 )
                ::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:);
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

				#define MUL_3xN_3x3_INIT( A , B , C )								\
				"mov		esi, "STRINGIZE(A)"\n"								\
				"mov		edi, "STRINGIZE(B)"\n"								\
				"mov		eax, "STRINGIZE(C)"\n"								\
				"movss		xmm0, [esi+(0*3+0)*4]\n"					\
				"movhps	xmm0, [esi+(0*3+1)*4]\n"					\
				"movss		xmm1, [esi+(1*3+0)*4]\n"					\
				"movhps	xmm1, [esi+(1*3+1)*4]\n"					\
				"movss		xmm2, [esi+(2*3+0)*4]\n"					\
				"movhps	xmm2, [esi+(2*3+1)*4]\n"

				#define MUL_3xN_3x3_INIT_ROW4							\
				"shufps	xmm0, xmm0, 0x38\n"  /*R_SHUFFLEPS( 0, 2, 3, 0 )=00111000*/	\
				"shufps	xmm1, xmm1, 0x38\n"  /*R_SHUFFLEPS( 0, 2, 3, 0 )=00111000*/	\
				"shufps	xmm2, xmm2, 0x38\n"  /*R_SHUFFLEPS( 0, 2, 3, 0 )=00111000*/

				#define MUL_3xN_3x3_ROW4( N, row )						\
				"movlps	xmm3, [edi+("STRINGIZE(row)"+0*"STRINGIZE(N)"+0)*4]\n"				\
				"shufps	xmm3, xmm3, 0x40\n"  /*R_SHUFFLEPS( 0, 0, 0, 1 )=01000000 */	\
				"movlps	xmm4, [edi+("STRINGIZE(row)"+1*"STRINGIZE(N)"+0)*4]\n"				\
				"shufps	xmm4, xmm4, 0x40\n"  /*R_SHUFFLEPS( 0, 0, 0, 1 )=01000000 */	\
				"movlps	xmm5, [edi+("STRINGIZE(row)"+2*"STRINGIZE(N)"+0)*4]\n"				\
				"shufps	xmm5, xmm5, 0x40\n"  /*R_SHUFFLEPS( 0, 0, 0, 1 )=01000000 */	\
				"mulps		xmm3, xmm0\n"								\
				"mulps		xmm4, xmm1\n"								\
				"mulps		xmm5, xmm2\n"								\
				"addps		xmm3, xmm4\n"								\
				"addps		xmm3, xmm5\n"								\
				"movaps	[eax+("STRINGIZE(row)"*3+0)*4], xmm3\n"					\
				"shufps	xmm0, xmm0, 0x79\n" /*R_SHUFFLEPS( 1, 2, 3, 1 )=01111001 */	\
				"shufps	xmm1, xmm1, 0x79\n" /*R_SHUFFLEPS( 1, 2, 3, 1 )=01111001 */	\
				"shufps	xmm2, xmm2, 0x79\n" /*R_SHUFFLEPS( 1, 2, 3, 1 )=01111001 */	\
				"movlps	xmm3, [edi+("STRINGIZE(row)"+0*"STRINGIZE(N)"+1)*4]\n"				\
				"shufps	xmm3, xmm3, 0x50\n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/	\
				"movlps	xmm4, [edi+("STRINGIZE(row)"+1*"STRINGIZE(N)"+1)*4]\n"				\
				"shufps	xmm4, xmm4, 0x50\n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/	\
				"movlps	xmm5, [edi+("STRINGIZE(row)"+2*"STRINGIZE(N)"+1)*4]\n"				\
				"shufps	xmm5, xmm5, 0x50\n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/	\
				"mulps		xmm3, xmm0\n"								\
				"mulps		xmm4, xmm1\n"								\
				"mulps		xmm5, xmm2\n"								\
				"addps		xmm3, xmm4\n"								\
				"addps		xmm3, xmm5\n"								\
				"movaps	[eax+("STRINGIZE(row)"*3+4)*4], xmm3\n"					\
				"shufps	xmm0, xmm0, 0x79\n" /*R_SHUFFLEPS( 1, 2, 3, 1 )=01111001 */	\
				"shufps	xmm1, xmm1, 0x79\n" /*R_SHUFFLEPS( 1, 2, 3, 1 )=01111001 */	\
				"shufps	xmm2, xmm2, 0x79\n" /*R_SHUFFLEPS( 1, 2, 3, 1 )=01111001 */	\
				"movlps	xmm3, [edi+("STRINGIZE(row)"+0*"STRINGIZE(N)"+2)*4]\n"				\
				"shufps	xmm3, xmm3, 0x54\n"  /*R_SHUFFLEPS( 0, 1, 1, 1 )=01010100*/	\
				"movlps	xmm4, [edi+("STRINGIZE(row)"+1*"STRINGIZE(N)"+2)*4]\n"				\
				"shufps	xmm4, xmm4, 0x54\n"  /*R_SHUFFLEPS( 0, 1, 1, 1 )=01010100*/		\
				"movlps	xmm5, [edi+("STRINGIZE(row)"+2*"STRINGIZE(N)"+2)*4]\n"				\
				"shufps	xmm5, xmm5, 0x54\n"  /*R_SHUFFLEPS( 0, 1, 1, 1 )=01010100*/		\
				"mulps		xmm3, xmm0\n"								\
				"mulps		xmm4, xmm1\n"								\
				"mulps		xmm5, xmm2\n"								\
				"addps		xmm3, xmm4\n"								\
				"addps		xmm3, xmm5\n"								\
				"movaps	[eax+("STRINGIZE(row)"*3+8)*4], xmm3\n"

				#define MUL_3xN_3x3_INIT_ROW4_ROW4						\
				"shufps	xmm0, xmm0, 0x39\n"   /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/	\
				"shufps	xmm1, xmm1, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/	\
				"shufps	xmm2, xmm2, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/

				#define MUL_3xN_3x3_INIT_ROW4_ROW						\
				"shufps	xmm0, xmm0, 0xe5\n"   /*R_SHUFFLEPS( 1, 1, 2, 3 )=11100101*/	\
				"shufps	xmm1, xmm1, 0xe5\n"  /*R_SHUFFLEPS( 1, 1, 2, 3 )=11100101*/	\
				"shufps	xmm2, xmm2, 0xe5\n"  /*R_SHUFFLEPS( 1, 1, 2, 3 )=11100101*/

				#define MUL_3xN_3x3_ROW( N, row )						\
				"movss		xmm3, [edi+("STRINGIZE(row)"+0*"STRINGIZE(N)")*4]\n"					\
				"shufps	xmm3, xmm3, 0x00\n"	\
				"movss		xmm4, [edi+("STRINGIZE(row)"+1*"STRINGIZE(N)")*4]\n"					\
				"shufps	xmm4, xmm4, 0x00\n"	\
				"movss		xmm5, [edi+("STRINGIZE(row)"+2*"STRINGIZE(N)")*4]\n"					\
				"shufps	xmm5, xmm5, 0x00\n"	\
				"mulps		xmm3, xmm0\n"								\
				"mulps		xmm4, xmm1\n"								\
				"mulps		xmm5, xmm2\n"								\
				"addps		xmm3, xmm4\n"								\
				"addps		xmm3, xmm5\n"								\
				"movss		[eax+("STRINGIZE(row)"*3+0)*4], xmm3\n"					\
				"movhps	[eax+("STRINGIZE(row)"*3+1)*4], xmm3\n"

                asm(
				MUL_3xN_3x3_INIT(%0 , %1, %2)
				MUL_3xN_3x3_INIT_ROW4
				MUL_3xN_3x3_ROW4( 6, 0 )
                ::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:);
                asm(
				MUL_3xN_3x3_INIT_ROW4_ROW
				MUL_3xN_3x3_ROW( 6, 4 )
				MUL_3xN_3x3_ROW( 6, 5 )
                );
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

				#define MUL_4xN_4x4_INIT( A , B , C )								\
				"mov		esi, "STRINGIZE(A)"\n"								\
				"mov		edi, "STRINGIZE(B)"\n"								\
				"mov		eax, "STRINGIZE(C)"\n"								\
				"movaps	xmm0, [esi]	\n"							\
				"movaps	xmm1, [esi+16]\n"							\
				"movaps	xmm2, [esi+32]\n"							\
				"movaps	xmm3, [esi+48]\n"

				#define MUL_4xN_4x4_ROW( N, row )						\
				"movss		xmm7, [edi+("STRINGIZE(row)"+0*"STRINGIZE(N)")*4]\n"					\
				"shufps	xmm7, xmm7, 0x00\n"	\
				"mulps		xmm7, xmm0\n"								\
				"movss		xmm6, [edi+("STRINGIZE(row)"+1*"STRINGIZE(N)")*4]\n"					\
				"shufps	xmm6, xmm6, 0x00\n"	\
				"mulps		xmm6, xmm1\n"								\
				"addps		xmm7, xmm6\n"								\
				"movss		xmm6, [edi+("STRINGIZE(row)"+2*"STRINGIZE(N)")*4]\n"					\
				"shufps	xmm6, xmm6, 0x00\n"	\
				"mulps		xmm6, xmm2\n"								\
				"addps		xmm7, xmm6\n"								\
				"movss		xmm6, [edi+("STRINGIZE(row)"+3*"STRINGIZE(N)")*4]\n"					\
				"shufps	xmm6, xmm6, 0x00\n"	\
				"mulps		xmm6, xmm3\n"								\
				"addps		xmm7, xmm6\n"								\
				"movaps	[eax+"STRINGIZE(row)"*16], xmm7\n"
                asm(
				MUL_4xN_4x4_INIT(%0 , %1, %2)
				MUL_4xN_4x4_ROW( 6, 0 )
				MUL_4xN_4x4_ROW( 6, 1 )
				MUL_4xN_4x4_ROW( 6, 2 )
				MUL_4xN_4x4_ROW( 6, 3 )
				MUL_4xN_4x4_ROW( 6, 4 )
				MUL_4xN_4x4_ROW( 6, 5 )
                ::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:);
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

				#define MUL_5xN_5x5_INIT( A , B , C )								\
				"mov		esi, "STRINGIZE(A)"\n"								\
				"mov		edi, "STRINGIZE(B)"\n"								\
				"mov		eax, "STRINGIZE(C)"\n"								\
				"movlps	xmm0, [esi+ 0*4]\n"						\
				"movhps	xmm0, [esi+ 2*4]\n"						\
				"movlps	xmm1, [esi+ 5*4]\n"						\
				"movhps	xmm1, [esi+ 7*4]\n"						\
				"movlps	xmm2, [esi+10*4]\n"						\
				"movhps	xmm2, [esi+12*4]\n"						\
				"movlps	xmm3, [esi+15*4]\n"						\
				"movhps	xmm3, [esi+17*4]\n"						\
				"movlps	xmm4, [esi+20*4]\n"						\
				"movhps	xmm4, [esi+22*4]\n"

				#define MUL_5xN_5x5_ROW( N, row )						\
				"movss		xmm6, [edi+("STRINGIZE(row)"+0*"STRINGIZE(N)")*4]\n"					\
				"shufps	xmm6, xmm6, 0x00\n"	\
				"mulps		xmm6, xmm0\n"								\
				"fld		dword ptr [edi+("STRINGIZE(row)"+0*"STRINGIZE(N)")*4]\n"				\
				"fmul		dword ptr [esi+ 4*4]\n"					\
				"movss		xmm5, [edi+("STRINGIZE(row)"+1*"STRINGIZE(N)")*4]\n"					\
				"shufps	xmm5, xmm5, 0x00\n"	\
				"mulps		xmm5, xmm1\n"								\
				"addps		xmm6, xmm5\n"								\
				"fld		dword ptr [edi+("STRINGIZE(row)"+1*"STRINGIZE(N)")*4]\n"				\
				"fmul		dword ptr [esi+ 9*4]\n"					\
				"faddp		st(1),st\n"								\
				"movss		xmm5, [edi+("STRINGIZE(row)"+2*"STRINGIZE(N)")*4]\n"					\
				"shufps	xmm5, xmm5, 0x00\n"	\
				"mulps		xmm5, xmm2\n"								\
				"addps		xmm6, xmm5\n"								\
				"fld		dword ptr [edi+("STRINGIZE(row)"+2*"STRINGIZE(N)")*4]\n"				\
				"fmul		dword ptr [esi+14*4]\n"					\
				"faddp		st(1),st\n"								\
				"movss		xmm5, [edi+("STRINGIZE(row)"+3*"STRINGIZE(N)")*4]\n"					\
				"shufps	xmm5, xmm5, 0x00\n"	\
				"mulps		xmm5, xmm3\n"								\
				"addps		xmm6, xmm5\n"								\
				"fld		dword ptr [edi+("STRINGIZE(row)"+3*"STRINGIZE(N)")*4]\n"				\
				"fmul		dword ptr [esi+19*4]\n"					\
				"faddp		st(1),st\n"								\
				"movss		xmm5, [edi+("STRINGIZE(row)"+4*"STRINGIZE(N)")*4]\n"					\
				"shufps	xmm5, xmm5, 0x00\n"	\
				"mulps		xmm5, xmm4\n"								\
				"addps		xmm6, xmm5\n"								\
				"fld		dword ptr [edi+("STRINGIZE(row)"+4*"STRINGIZE(N)")*4]\n"				\
				"fmul		dword ptr [esi+24*4]\n"					\
				"faddp		st(1),st\n"								\
				"fstp		dword ptr [eax+("STRINGIZE(row)"*5+4)*4]\n"				\
				"movlps	[eax+("STRINGIZE(row)"*5+0)*4], xmm6\n"					\
				"movhps	[eax+("STRINGIZE(row)"*5+2)*4], xmm6\n"
                asm(
				MUL_5xN_5x5_INIT(%0 , %1, %2)
				MUL_5xN_5x5_ROW( 6, 0 )
				MUL_5xN_5x5_ROW( 6, 1 )
				MUL_5xN_5x5_ROW( 6, 2 )
				MUL_5xN_5x5_ROW( 6, 3 )
				MUL_5xN_5x5_ROW( 6, 4 )
				MUL_5xN_5x5_ROW( 6, 5 )
                ::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:);
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
						#define MUL_6xN_6x6_FIRST4COLUMNS_INIT( A , B , C )					\
						"mov		esi, "STRINGIZE(A)"\n"								\
						"mov		edi, "STRINGIZE(B)"\n"								\
						"mov		eax, "STRINGIZE(C)"\n"								\
						"movlps	xmm0, [esi+ 0*4]\n"						\
						"movhps	xmm0, [esi+ 2*4]\n"						\
						"movlps	xmm1, [esi+ 6*4]\n"						\
						"movhps	xmm1, [esi+ 8*4]\n"						\
						"movlps	xmm2, [esi+12*4]\n"						\
						"movhps	xmm2, [esi+14*4]\n"						\
						"movlps	xmm3, [esi+18*4]\n"						\
						"movhps	xmm3, [esi+20*4]\n"						\
						"movlps	xmm4, [esi+24*4]\n"						\
						"movhps	xmm4, [esi+26*4]\n"						\
						"movlps	xmm5, [esi+30*4]\n"						\
						"movhps	xmm5, [esi+32*4]\n"

						#define MUL_6xN_6x6_FIRST4COLUMNS_ROW( N, row )			\
						"movss		xmm7, [edi+("STRINGIZE(row)"+0*"STRINGIZE(N)")*4]\n"					\
						"shufps	xmm7, xmm7, 0x00\n"	\
						"mulps		xmm7, xmm0\n"								\
						"movss		xmm6, [edi+("STRINGIZE(row)"+1*"STRINGIZE(N)")*4]\n"					\
						"shufps	xmm6, xmm6, 0x00\n"	\
						"mulps		xmm6, xmm1\n"								\
						"addps		xmm7, xmm6\n"								\
						"movss		xmm6, [edi+("STRINGIZE(row)"+2*"STRINGIZE(N)")*4]\n"					\
						"shufps	xmm6, xmm6, 0x00\n"	\
						"mulps		xmm6, xmm2\n"								\
						"addps		xmm7, xmm6\n"								\
						"movss		xmm6, [edi+("STRINGIZE(row)"+3*"STRINGIZE(N)")*4]\n"					\
						"shufps	xmm6, xmm6, 0x00\n"	\
						"mulps		xmm6, xmm3\n"								\
						"addps		xmm7, xmm6\n"								\
						"movss		xmm6, [edi+("STRINGIZE(row)"+4*"STRINGIZE(N)")*4]\n"					\
						"shufps	xmm6, xmm6, 0x00\n"	\
						"mulps		xmm6, xmm4\n"								\
						"addps		xmm7, xmm6\n"								\
						"movss		xmm6, [edi+("STRINGIZE(row)"+5*"STRINGIZE(N)")*4]\n"					\
						"shufps	xmm6, xmm6, 0x00\n"	\
						"mulps		xmm6, xmm5\n"								\
						"addps		xmm7, xmm6\n"								\
						"movlps	[eax+("STRINGIZE(row)"*6+0)*4], xmm7\n"					\
						"movhps	[eax+("STRINGIZE(row)"*6+2)*4], xmm7\n"

						#define MUL_6xN_6x6_LAST2COLUMNS_INIT					\
						"movlps	xmm0, [esi+ 4*4]\n"						\
						"movlps	xmm1, [esi+10*4]\n"						\
						"shufps	xmm0, xmm0, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/	\
						"shufps	xmm1, xmm1, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/	\
						"movlps	xmm2, [esi+16*4]\n"						\
						"movlps	xmm3, [esi+22*4]\n"						\
						"shufps	xmm2, xmm2, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/	\
						"shufps	xmm3, xmm3, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/	\
						"movlps	xmm4, [esi+28*4]\n"						\
						"movlps	xmm5, [esi+34*4]\n"						\
						"shufps	xmm4, xmm4, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/	\
						"shufps	xmm5, xmm5, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/

						#define MUL_6xN_6x6_LAST2COLUMNS_ROW2( N , row )			\
						"movlps	xmm7, [edi+("STRINGIZE(row)"*2+0*"STRINGIZE(N)")*4]\n"				\
						"shufps	xmm7, xmm7, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/	\
						"mulps		xmm7, xmm0\n"								\
						"movlps	xmm6, [edi+("STRINGIZE(row)"*2+1*"STRINGIZE(N)")*4]\n"				\
						"shufps	xmm6, xmm6, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/	\
						"mulps		xmm6, xmm1\n"								\
						"addps		xmm7, xmm6\n"								\
						"movlps	xmm6, [edi+("STRINGIZE(row)"*2+2*"STRINGIZE(N)")*4]\n"				\
						"shufps	xmm6, xmm6, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/	\
						"mulps		xmm6, xmm2\n"								\
						"addps		xmm7, xmm6\n"								\
						"movlps	xmm6, [edi+("STRINGIZE(row)"*2+3*"STRINGIZE(N)")*4]\n"				\
						"shufps	xmm6, xmm6, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/	\
						"mulps		xmm6, xmm3\n"								\
						"addps		xmm7, xmm6\n"								\
						"movlps	xmm6, [edi+("STRINGIZE(row)"*2+4*"STRINGIZE(N)")*4]\n"				\
						"shufps	xmm6, xmm6, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/	\
						"mulps		xmm6, xmm4\n"								\
						"addps		xmm7, xmm6\n"								\
						"movlps	xmm6, [edi+("STRINGIZE(row)"*2+5*"STRINGIZE(N)")*4]\n"				\
						"shufps	xmm6, xmm6, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/	\
						"mulps		xmm6, xmm5\n"								\
						"addps		xmm7, xmm6\n"								\
						"movlps	[eax+("STRINGIZE(row)"*12+ 4)*4], xmm7\n"				\
						"movhps	[eax+("STRINGIZE(row)"*12+10)*4], xmm7\n"

						#define MUL_6xN_6x6_LAST2COLUMNS_ROW( N, row )			\
						"movss		xmm7, [edi+(1*"STRINGIZE(N)"-1)*4]\n"					\
						"shufps	xmm7, xmm7, 0x00\n"	\
						"mulps		xmm7, xmm0\n"								\
						"movss		xmm6, [edi+(2*"STRINGIZE(N)"-1)*4]\n"					\
						"shufps	xmm6, xmm6, 0x00\n"	\
						"mulps		xmm6, xmm1\n"								\
						"addps		xmm7, xmm6\n"								\
						"movss		xmm6, [edi+(3*"STRINGIZE(N)"-1)*4]\n"					\
						"shufps	xmm6, xmm6, 0x00\n"	\
						"mulps		xmm6, xmm2\n"								\
						"addps		xmm7, xmm6\n"								\
						"movss		xmm6, [edi+(4*"STRINGIZE(N)"-1)*4]\n"					\
						"shufps	xmm6, xmm6, 0x00\n"	\
						"mulps		xmm6, xmm3\n"								\
						"addps		xmm7, xmm6\n"								\
						"movss		xmm6, [edi+(5*"STRINGIZE(N)"-1)*4]\n"					\
						"shufps	xmm6, xmm6, 0x00\n"	\
						"mulps		xmm6, xmm4\n"								\
						"addps		xmm7, xmm6\n"								\
						"movss		xmm6, [edi+(6*"STRINGIZE(N)"-1)*4]\n"					\
						"shufps	xmm6, xmm6, 0x00\n"	\
						"mulps		xmm6, xmm5\n"								\
						"addps		xmm7, xmm6\n"								\
						"movlps	[eax+("STRINGIZE(row)"*6+4)*4], xmm7\n"

                        asm(
						MUL_6xN_6x6_FIRST4COLUMNS_INIT(%0 , %1, %2)
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 1, 0 )
						MUL_6xN_6x6_LAST2COLUMNS_INIT
						MUL_6xN_6x6_LAST2COLUMNS_ROW( 1, 0 )
                        ::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:);
						return;
					}
					case 2: {					// 6x2 * 6x6
                        asm(
						MUL_6xN_6x6_FIRST4COLUMNS_INIT(%0 , %1, %2)
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 2, 0 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 2, 1 )
						MUL_6xN_6x6_LAST2COLUMNS_INIT
						MUL_6xN_6x6_LAST2COLUMNS_ROW2( 2, 0 )
                        ::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:);
						return;
					}
					case 3: {					// 6x3 * 6x6
                        asm(
						MUL_6xN_6x6_FIRST4COLUMNS_INIT(%0 , %1, %2)
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 3, 0 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 3, 1 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 3, 2 )
						MUL_6xN_6x6_LAST2COLUMNS_INIT
						MUL_6xN_6x6_LAST2COLUMNS_ROW2( 3, 0 )
						MUL_6xN_6x6_LAST2COLUMNS_ROW( 3, 2 )
                        ::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:);
						return;
					}
					case 4: {					// 6x4 * 6x6
                        asm(
						MUL_6xN_6x6_FIRST4COLUMNS_INIT(%0 , %1, %2)
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 4, 0 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 4, 1 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 4, 2 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 4, 3 )
						::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:);
						asm(
						MUL_6xN_6x6_LAST2COLUMNS_INIT
						MUL_6xN_6x6_LAST2COLUMNS_ROW2( 4, 0 )
						MUL_6xN_6x6_LAST2COLUMNS_ROW2( 4, 1 )
                        );
						return;
					}
					case 5: {					// 6x5 * 6x6
                        asm(
						MUL_6xN_6x6_FIRST4COLUMNS_INIT(%0 , %1, %2)
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 5, 0 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 5, 1 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 5, 2 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 5, 3 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 5, 4 )
						::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:);
						asm(
						MUL_6xN_6x6_LAST2COLUMNS_INIT
						MUL_6xN_6x6_LAST2COLUMNS_ROW2( 5, 0 )
						MUL_6xN_6x6_LAST2COLUMNS_ROW2( 5, 1 )
						MUL_6xN_6x6_LAST2COLUMNS_ROW( 5, 4 )
                        );
						return;
					}
					case 6: {					// 6x6 * 6x6
                        asm(
						MUL_6xN_6x6_FIRST4COLUMNS_INIT(%0 , %1, %2)
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 6, 0 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 6, 1 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 6, 2 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 6, 3 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 6, 4 )
						MUL_6xN_6x6_FIRST4COLUMNS_ROW( 6, 5 )
						::"r"(m2Ptr),
						"r"(m1Ptr),
						"r"(dstPtr)
						:);
						asm(
						MUL_6xN_6x6_LAST2COLUMNS_INIT
						MUL_6xN_6x6_LAST2COLUMNS_ROW2( 6, 0 )
						MUL_6xN_6x6_LAST2COLUMNS_ROW2( 6, 1 )
						MUL_6xN_6x6_LAST2COLUMNS_ROW2( 6, 2 )
                        );
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
CSIMD_SSE::matX_LowerTriangularSolve

  solves x in Lx = b for the n * n "sub-matrix of L
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
	asm (
		"push		ebx\n"
		"mov		eax, %4\n"				// eax = i
		"shl		eax, 2\n"					// eax = i*4
		"mov		edx, %3\n"					// edx = n
		"shl		edx, 2\n"					// edx = n*4
		"mov		esi, %1\n"					// esi = x
		"mov		edi, %0\n"				// edi = lptr
		"add		esi, eax\n"
		"add		edi, eax\n"
		"mov		ebx, %2\n"					// ebx = b

		// check for aligned memory
		"mov		ecx, %5\n"
		"shl		ecx, 2\n"
		"or			ecx, esi\n"
		"or			ecx, edi\n"
		"and		ecx, 15\n"
		"jnz		loopurow\n"

		// aligned
	"looprow:\n"
		"mov		ecx, eax\n"
		"neg		ecx\n"
		"movaps		xmm0, [esi+ecx]\n"
		"mulps		xmm0, [edi+ecx]\n"
		"add		ecx, 12*4\n"
		"jg			donedot8\n"
	"dot8:\n"
		"movaps		xmm1, [esi+ecx-(8*4)]\n"
		"mulps		xmm1, [edi+ecx-(8*4)]\n"
		"addps		xmm0, xmm1\n"
		"movaps		xmm3, [esi+ecx-(4*4)]\n"
		"mulps		xmm3, [edi+ecx-(4*4)]\n"
		"addps		xmm0, xmm3\n"
		"add		ecx, 8*4\n"
		"jle		dot8\n"
	"donedot8:\n"
		"sub		ecx, 4*4\n"
		"jg			donedot4\n"
	//dot4:
		"movaps		xmm1, [esi+ecx-(4*4)]\n"
		"mulps		xmm1, [edi+ecx-(4*4)]\n"
		"addps		xmm0, xmm1\n"
		"add		ecx, 4*4\n"
	"donedot4:\n"
		"movhlps	xmm1, xmm0\n"
		"addps		xmm0, xmm1\n"
		"movaps		xmm1, xmm0\n"
		"shufps		xmm1, xmm1, 0x01 \n"  /* R_SHUFFLEPS( 1, 0, 0, 0 )=00000001 */
		"addss		xmm0, xmm1\n"
		"sub		ecx, 4*4\n"
		"jz			dot0\n"
		"add		ecx, 4\n"
		"jz			dot1\n"
		"add		ecx, 4\n"
		"jz			dot2\n"
	//dot3:
		"movss		xmm1, [esi-(3*4)]\n"
		"mulss		xmm1, [edi-(3*4)]\n"
		"addss		xmm0, xmm1\n"
	"dot2:\n"
		"movss		xmm3, [esi-(2*4)]\n"
		"mulss		xmm3, [edi-(2*4)]\n"
		"addss		xmm0, xmm3\n"
	"dot1:\n"
		"movss		xmm5, [esi-(1*4)]\n"
		"mulss		xmm5, [edi-(1*4)]\n"
		"addss		xmm0, xmm5\n"
	"dot0:\n"
		"movss		xmm1, [ebx+eax]\n"
		"subss		xmm1, xmm0\n"
		"movss		[esi], xmm1\n"
		"add		eax, 4\n"
		"cmp		eax, edx\n"
		"jge		doneLocal\n"
		"add		esi, 4\n"
		"mov		ecx, %5\n"
		"shl		ecx, 2\n"
		"add		edi, ecx\n"
		"add		edi, 4\n"
		"jmp		looprow\n"

		// unaligned
	"loopurow:\n"
		"mov		ecx, eax\n"
		"neg		ecx\n"
		"movups		xmm0, [esi+ecx]\n"
		"movups		xmm1, [edi+ecx]\n"
		"mulps		xmm0, xmm1\n"
		"add		ecx, 12*4\n"
		"jg			doneudot8\n"
	"udot8:\n"
		"movups		xmm1, [esi+ecx-(8*4)]\n"
		"movups		xmm2, [edi+ecx-(8*4)]\n"
		"mulps		xmm1, xmm2\n"
		"addps		xmm0, xmm1\n"
		"movups		xmm3, [esi+ecx-(4*4)]\n"
		"movups		xmm4, [edi+ecx-(4*4)]\n"
		"mulps		xmm3, xmm4\n"
		"addps		xmm0, xmm3\n"
		"add		ecx, 8*4\n"
		"jle		udot8\n"
	"doneudot8:\n"
		"sub		ecx, 4*4\n"
		"jg			doneudot4\n"
	//udot4:
		"movups		xmm1, [esi+ecx-(4*4)]\n"
		"movups		xmm2, [edi+ecx-(4*4)]\n"
		"mulps		xmm1, xmm2\n"
		"addps		xmm0, xmm1\n"
		"add		ecx, 4*4\n"
	"doneudot4:\n"
		"movhlps	xmm1, xmm0\n"
		"addps		xmm0, xmm1\n"
		"movaps		xmm1, xmm0\n"
		"shufps		xmm1, xmm1, 0x01 \n"  /* R_SHUFFLEPS( 1, 0, 0, 0 )=00000001 */
		"addss		xmm0, xmm1\n"
		"sub		ecx, 4*4\n"
		"jz			udot0\n"
		"add		ecx, 4\n"
		"jz			udot1\n"
		"add		ecx, 4\n"
		"jz			udot2\n"
	//udot3:
		"movss		xmm1, [esi-(3*4)]\n"
		"movss		xmm2, [edi-(3*4)]\n"
		"mulss		xmm1, xmm2\n"
		"addss		xmm0, xmm1\n"
	"udot2:\n"
		"movss		xmm3, [esi-(2*4)]\n"
		"movss		xmm4, [edi-(2*4)]\n"
		"mulss		xmm3, xmm4\n"
		"addss		xmm0, xmm3\n"
	"udot1:\n"
		"movss		xmm5, [esi-(1*4)]\n"
		"movss		xmm6, [edi-(1*4)]\n"
		"mulss		xmm5, xmm6\n"
		"addss		xmm0, xmm5\n"
	"udot0:\n"
		"movss		xmm1, [ebx+eax]\n"
		"subss		xmm1, xmm0\n"
		"movss		[esi], xmm1\n"
		"add		eax, 4\n"
		"cmp		eax, edx\n"
		"jge		doneLocal\n"
		"add		esi, 4\n"
		"mov		ecx, %5\n"
		"shl		ecx, 2\n"
		"add		edi, ecx\n"
		"add		edi, 4\n"
		"jmp		loopurow\n"
	"doneLocal:\n"
		"pop			ebx\n"
	::"m"(lptr),"m"(x), "m"(b), "m"(n), "m"(skip),"m"(nc)
    :"eax","ebx","ecx","edx","esi","edi");
}

/*
============
CSIMD_SSE::matX_LowerTriangularSolveTranspose

  solves x in L'x = b for the n * n "sub-matrix of L
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
		asm (
			"push		ebx\n"
			"mov			eax, %0\n"					// eax = i
			"mov			esi, %1\n"				// esi = xptr
			"mov			edi, %3\n"				// edi = lptr
			"mov			ebx, %4\n"				// ebx = b
			"mov			edx, %2\n"					// edx = nc*sizeof(float)
			"shl			edx, 2\n"
		"process4rows_1:\n"
			"movlps		xmm0, [ebx+eax*4-16]\n"	// load b[i-2], b[i-1]
			"movhps		xmm0, [ebx+eax*4-8]\n"		// load b[i-4], b[i-3]
			"xor		ecx, ecx\n"
			"sub		eax, %0\n"
			"neg		eax\n"
			"jz			done4x4_1\n"
		"process4x4_1:\n"	// process 4x4 blocks
			"movlps		xmm2, [edi+0]\n"
			"movhps		xmm2, [edi+8]\n"
			"add		edi, edx\n"
			"movss		xmm1, [esi+4*ecx+0]\n"
			"shufps		xmm1, xmm1, 0x00\n"
			"movlps		xmm3, [edi+0]\n"
			"movhps		xmm3, [edi+8]\n"
			"add		edi, edx\n"
			"mulps		xmm1, xmm2\n"
			"subps		xmm0, xmm1\n"
			"movss		xmm1, [esi+4*ecx+4]\n"
			"shufps		xmm1, xmm1, 0x00\n"
			"movlps		xmm4, [edi+0]\n"
			"movhps		xmm4, [edi+8]\n"
			"add		edi, edx\n"
			"mulps		xmm1, xmm3\n"
			"subps		xmm0, xmm1\n"
			"movss		xmm1, [esi+4*ecx+8]\n"
			"shufps		xmm1, xmm1, 0x00\n"
			"movlps		xmm5, [edi+0]\n"
			"movhps		xmm5, [edi+8]\n"
			"add		edi, edx\n"
			"mulps		xmm1, xmm4\n"
			"subps		xmm0, xmm1\n"
			"movss		xmm1, [esi+4*ecx+12]\n"
			"shufps		xmm1, xmm1, 0x00\n"
			"add		ecx, 4\n"
			"cmp		ecx, eax\n"
			"mulps		xmm1, xmm5\n"
			"subps		xmm0, xmm1\n"
			"jl			process4x4_1\n"
		"done4x4_1:\n"		// process left over of the 4 rows
			"movlps		xmm2, [edi+0]\n"
			"movhps		xmm2, [edi+8]\n"
			"movss		xmm1, [esi+4*ecx]\n"
			"shufps		xmm1, xmm1, 0x00\n"
			"mulps		xmm1, xmm2\n"
			"subps		xmm0, xmm1\n"
			"imul		ecx, edx\n"
			"sub		edi, ecx\n"
			"neg		eax\n"

			"add		eax, %0\n"
			"sub		eax, 4\n"
			"movaps		xmm1, xmm0\n"
			"shufps		xmm1, xmm1, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
			"movaps		xmm2, xmm0\n"
			"shufps		xmm2, xmm2, 0xAA \n" /*R_SHUFFLEPS( 2, 2, 2, 2 )=10101010*/
			"movaps		xmm3, xmm0\n"
			"shufps		xmm3, xmm3, 0xFF \n" /*R_SHUFFLEPS( 3, 3, 3, 3 )=11111111*/
			"sub		edi, edx\n"
			"movss		[esi-4], xmm3\n"			// xptr[-1] = s3
			"movss		xmm4, xmm3\n"
			"movss		xmm5, xmm3\n"
			"mulss		xmm3, [edi+8]\n"			// lptr[-1*nc+2] * s3
			"mulss		xmm4, [edi+4]\n"			// lptr[-1*nc+1] * s3
			"mulss		xmm5, [edi+0]\n"			// lptr[-1*nc+0] * s3
			"subss		xmm2, xmm3\n"
			"movss		[esi-8], xmm2\n"			// xptr[-2] = s2
			"movss		xmm6, xmm2\n"
			"sub		edi, edx\n"
			"subss		xmm0, xmm5\n"
			"subss		xmm1, xmm4\n"
			"mulss		xmm2, [edi+4]\n"			// lptr[-2*nc+1] * s2
			"mulss		xmm6, [edi+0]\n"			// lptr[-2*nc+0] * s2
			"subss		xmm1, xmm2\n"
			"movss		[esi-12], xmm1\n"			// xptr[-3] = s1
			"subss		xmm0, xmm6\n"
			"sub		edi, edx\n"
			"cmp		eax, 4\n"
			"mulss		xmm1, [edi+0]\n"			// lptr[-3*nc+0] * s1
			"subss		xmm0, xmm1\n"
			"movss		[esi-16], xmm0\n"			// xptr[-4] = s0
			"jl			done4rows_1\n"
			"sub		edi, edx\n"
			"sub		edi, 16\n"
			"sub		esi, 16\n"
			"jmp		process4rows_1\n"
		"done4rows_1:\n"
			"pop			ebx\n"
		::	"m"(m),"m"(xptr),"m"(nc),"m"(lptr),"m"(b)
        :"eax","ebx","ecx","edx","esi","edi");

	} else {

		lptr = L.toFloatPtr() + m * nc + m - 4;
		xptr = x + m;
		asm (
			"push		ebx\n"
			"mov			eax, %0\n"					// eax = i
			"mov			esi, %1\n"				// esi = xptr
			"mov			edi, %3\n"				// edi = lptr
			"mov			ebx, %4\n"					// ebx = b
			"mov			edx, %2\n"					// edx = nc*sizeof(float)
			"shl			edx, 2\n"
		"process4rows:\n"
			"movlps		xmm0, [ebx+eax*4-16]\n"	// load b[i-2], b[i-1]
			"movhps		xmm0, [ebx+eax*4-8]\n"		// load b[i-4], b[i-3]
			"sub			eax, %0\n"
			"jz			done4x4\n"
			"neg			eax\n"
			"xor			ecx, ecx\n"
		"process4x4:\n"		// process 4x4 blocks
			"movlps		xmm2, [edi+0]\n"
			"movhps		xmm2, [edi+8]\n"
			"add			edi, edx\n"
			"movss		xmm1, [esi+4*ecx+0]\n"
			"shufps		xmm1, xmm1, 0x00\n"
			"movlps		xmm3, [edi+0]\n"
			"movhps		xmm3, [edi+8]\n"
			"add			edi, edx\n"
			"mulps		xmm1, xmm2\n"
			"subps		xmm0, xmm1\n"
			"movss		xmm1, [esi+4*ecx+4]\n"
			"shufps		xmm1, xmm1, 0x00\n"
			"movlps		xmm4, [edi+0]\n"
			"movhps		xmm4, [edi+8]\n"
			"add			edi, edx\n"
			"mulps		xmm1, xmm3\n"
			"subps		xmm0, xmm1\n"
			"movss		xmm1, [esi+4*ecx+8]\n"
			"shufps		xmm1, xmm1, 0x00\n"
			"movlps		xmm5, [edi+0]\n"
			"movhps		xmm5, [edi+8]\n"
			"add			edi, edx\n"
			"mulps		xmm1, xmm4\n"
			"subps		xmm0, xmm1\n"
			"movss		xmm1, [esi+4*ecx+12]\n"
			"shufps		xmm1, xmm1, 0x00\n"
			"add			ecx, 4\n"
			"cmp			ecx, eax\n"
			"mulps		xmm1, xmm5\n"
			"subps		xmm0, xmm1\n"
			"jl			process4x4\n"
			"imul		ecx, edx\n"
			"sub			edi, ecx\n"
			"neg			eax\n"
		"done4x4:\n"		// process left over of the 4 rows
			"add			eax, %0\n"
			"sub			eax, 4\n"
			"movaps		xmm1, xmm0\n"
			"shufps		xmm1, xmm1, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
			"movaps		xmm2, xmm0\n"
			"shufps		xmm2, xmm2, 0xAA \n" /*R_SHUFFLEPS( 2, 2, 2, 2 )=10101010*/
			"movaps		xmm3, xmm0\n"
			"shufps		xmm3, xmm3, 0xFF \n" /*R_SHUFFLEPS( 3, 3, 3, 3 )=11111111*/
			"sub			edi, edx\n"
			"movss		[esi-4], xmm3\n"			// xptr[-1] = s3
			"movss		xmm4, xmm3\n"
			"movss		xmm5, xmm3\n"
			"mulss		xmm3, [edi+8]\n"			// lptr[-1*nc+2] * s3
			"mulss		xmm4, [edi+4]\n"			// lptr[-1*nc+1] * s3
			"mulss		xmm5, [edi+0]\n"			// lptr[-1*nc+0] * s3
			"subss		xmm2, xmm3\n"
			"movss		[esi-8], xmm2\n"			// xptr[-2] = s2
			"movss		xmm6, xmm2\n"
			"sub			edi, edx\n"
			"subss		xmm0, xmm5\n"
			"subss		xmm1, xmm4\n"
			"mulss		xmm2, [edi+4]\n"			// lptr[-2*nc+1] * s2
			"mulss		xmm6, [edi+0]\n"			// lptr[-2*nc+0] * s2
			"subss		xmm1, xmm2\n"
			"movss		[esi-12], xmm1\n"			// xptr[-3] = s1
			"subss		xmm0, xmm6\n"
			"sub			edi, edx\n"
			"cmp			eax, 4\n"
			"mulss		xmm1, [edi+0]\n"			// lptr[-3*nc+0] * s1
			"subss		xmm0, xmm1\n"
			"movss		[esi-16], xmm0\n"			// xptr[-4] = s0
			"jl			done4rows\n"
			"sub			edi, edx\n"
			"sub			edi, 16\n"
			"sub			esi, 16\n"
			"jmp			process4rows\n"
		"done4rows:\n"
			"pop			ebx\n"
		::	"m"(m),"m"(xptr),"m"(nc),"m"(lptr),"m"(b)
        :"eax","ebx","ecx","edx","esi","edi");
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

  in-place factorization LDL' of the n * n "sub-matrix of mat
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

	asm (
		"xorps		xmm2, xmm2\n"
		"xorps		xmm3, xmm3\n"
		"xorps		xmm4, xmm4\n"

		"push		ebx\n"
		"mov		ebx, 4\n"

	"loopRow:\n"
			"cmp			ebx, %7\n"
			"jge			doneLDLTFactor\n"

			"mov			ecx, ebx\n"				// esi = i
			"shl			ecx, 2\n"					// esi = i * 4
			"mov			edx, %4\n"				// edx = diag
			"add			edx, ecx\n"				// edx = &diag[i]
			"mov			edi, ebx\n"				// edi = i
			"imul		edi, %2\n"				// edi = i * nc * sizeof( float )
			"add		edi, %5\n"				// edi = mat[i]
			"add		edi, ecx\n"				// edi = &mat[i][i]
			"mov		esi, %6\n"					// ecx = v
			"add		esi, ecx\n"				// ecx = &v[i]
			"mov		eax, %3\n"			// eax = invDiagPtr
			"add		eax, ecx\n"				// eax = &invDiagPtr[i]
			"neg		ecx\n"

			"movaps		xmm0, [edx+ecx]\n"
			"mulps		xmm0, [edi+ecx]\n"
			"movaps		[esi+ecx], xmm0\n"
			"mulps		xmm0, [edi+ecx]\n"
			"add		ecx, 12*4\n"
			"jg			doneDot8\n"
		"dot8%=:\n"
			"movaps		xmm1, [edx+ecx-(8*4)]\n"
			"mulps		xmm1, [edi+ecx-(8*4)]\n"
			"movaps		[esi+ecx-(8*4)], xmm1\n"
			"mulps		xmm1, [edi+ecx-(8*4)]\n"
			"addps		xmm0, xmm1\n"
			"movaps		xmm2, [edx+ecx-(4*4)]\n"
			"mulps		xmm2, [edi+ecx-(4*4)]\n"
			"movaps		[esi+ecx-(4*4)], xmm2\n"
			"mulps		xmm2, [edi+ecx-(4*4)]\n"
			"addps		xmm0, xmm2\n"
			"add		ecx, 8*4\n"
			"jle		dot8%=\n"
		"doneDot8:\n"
			"sub		ecx, 4*4\n"
			"jg			doneDot4\n"
			"movaps		xmm1, [edx+ecx-(4*4)]\n"
			"mulps		xmm1, [edi+ecx-(4*4)]\n"
			"movaps		[esi+ecx-(4*4)], xmm1\n"
			"mulps		xmm1, [edi+ecx-(4*4)]\n"
			"addps		xmm0, xmm1\n"
			"add		ecx, 4*4\n"
		"doneDot4:\n"
			"sub		ecx, 2*4\n"
			"jg			doneDot2\n"
			"movlps		xmm3, [edx+ecx-(2*4)]\n"
			"movlps		xmm4, [edi+ecx-(2*4)]\n"
			"mulps		xmm3, xmm4\n"
			"movlps		[esi+ecx-(2*4)], xmm3\n"
			"mulps		xmm3, xmm4\n"
			"addps		xmm0, xmm3\n"
			"add		ecx, 2*4\n"
		"doneDot2:\n"
			"sub		ecx, 1*4\n"
			"jg			doneDot1\n"
			"movss		xmm3, [edx+ecx-(1*4)]\n"
			"movss		xmm4, [edi+ecx-(1*4)]\n"
			"mulss		xmm3, xmm4\n"
			"movss		[esi+ecx-(1*4)], xmm3\n"
			"mulss		xmm3, xmm4\n"
			"addss		xmm0, xmm3\n"
		"doneDot1:\n"
			"movhlps	xmm2, xmm0\n"
			"addps		xmm0, xmm2\n"
			"movaps		xmm2, xmm0\n"
			"shufps		xmm2, xmm2, 0x01 \n"  /* R_SHUFFLEPS( 1, 0, 0, 0 )=00000001 */
			"addss		xmm0, xmm2\n"
			"movss		xmm1, [edi]\n"
			"subss		xmm1, xmm0\n"
			"movss		[edi], xmm1\n"				// mptr[i] = sum;
			"movss		[edx], xmm1\n"				// diag[i] = sum;

			// if ( sum == 0.0f ) return false;
			"movaps		xmm2, xmm1\n"
			"cmpeqss	xmm2, %0\n"
			"andps		xmm2, %1\n"
			"orps		xmm1, xmm2\n"

			"rcpss		xmm7, xmm1\n"
			"mulss		xmm1, xmm7\n"
			"mulss		xmm1, xmm7\n"
			"addss		xmm7, xmm7\n"
			"subss		xmm7, xmm1\n"
			"movss		[eax], xmm7\n"				// invDiagPtr[i] = 1.0f / sum;

			"mov			edx, %7\n"					// edx = n
			"sub			edx, ebx\n"				// edx = n - i
			"dec			edx\n"						// edx = n - i - 1
			"jle			donesubRow\n"				// if ( i + 1 >= n ) return true;

			"mov			eax, ebx\n"				// eax = i
			"shl			eax, 2\n"			// eax = i * 4
			"neg			eax\n"

		"loopsubRow:\n"
				"add		edi, %2\n"
				"mov		ecx, eax\n"
				"movaps		xmm0, [esi+ecx]\n"
				"mulps		xmm0, [edi+ecx]\n"
				"add		ecx, 12*4\n"
				"jg			donesubDot8\n"
			"subDot8:\n"
				"movaps		xmm1, [esi+ecx-(8*4)]\n"
				"mulps		xmm1, [edi+ecx-(8*4)]\n"
				"addps		xmm0, xmm1\n"
				"movaps		xmm2, [esi+ecx-(4*4)]\n"
				"mulps		xmm2, [edi+ecx-(4*4)]\n"
				"addps		xmm0, xmm2\n"
				"add		ecx, 8*4\n"
				"jle		subDot8\n"
			"donesubDot8:\n"
				"sub		ecx, 4*4\n"
				"jg			donesubDot4\n"
				"movaps		xmm1, [esi+ecx-(4*4)]\n"
				"mulps		xmm1, [edi+ecx-(4*4)]\n"
				"addps		xmm0, xmm1\n"
				"add		ecx, 4*4\n"
			"donesubDot4:\n"
				"sub		ecx, 2*4\n"
				"jg			donesubDot2\n"
				"movlps		xmm3, [esi+ecx-(2*4)]\n"
				"movlps		xmm4, [edi+ecx-(2*4)]\n"
				"mulps		xmm3, xmm4\n"
				"addps		xmm0, xmm3\n"
				"add		ecx, 2*4\n"
			"donesubDot2:\n"
				"sub		ecx, 1*4\n"
				"jg			donesubDot1\n"
				"movss		xmm3, [esi+ecx-(1*4)]\n"
				"movss		xmm4, [edi+ecx-(1*4)]\n"
				"mulss		xmm3, xmm4\n"
				"addss		xmm0, xmm3\n"
			"donesubDot1:\n"
				"movhlps	xmm2, xmm0\n"
				"addps		xmm0, xmm2\n"
				"movaps		xmm2, xmm0\n"
				"shufps		xmm2, xmm2, 0x01 \n"  /* R_SHUFFLEPS( 1, 0, 0, 0 )=00000001 */
				"addss		xmm0, xmm2\n"
				"movss		xmm1, [edi]\n"
				"subss		xmm1, xmm0\n"
				"mulss		xmm1, xmm7\n"
				"movss		[edi], xmm1\n"
				"dec		edx\n"
				"jg			loopsubRow\n"
		"donesubRow:\n"
			"inc		ebx\n"
			"jmp		loopRow\n"
	"doneLDLTFactor:\n"
		"pop		ebx\n"
	::"m"(SIMD_SP_zero),"m"(SIMD_SP_tiny),"m"(ncf),"m"(invDiagPtr),"m"(diag),"m"(mptr),"m"(v),"m"(n)
    :"eax","ebx","ecx","edx","esi","edi");

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

        ALIGN16( void * jointVert0Ptr =&jointVert0[4] );
		ALIGN16( void * jointVert1Ptr =&jointVert1[4] );
		ALIGN16( void * jointVert2Ptr =&jointVert2[4] );
		ALIGN16( void * blendVert0Ptr =&blendVert0[4] );
		ALIGN16( void * blendVert1Ptr =&blendVert1[4] );
		ALIGN16( void * blendVert2Ptr =&blendVert2[4] );
		ALIGN16( void * jointQuat0Ptr =&jointQuat0[4] );
		ALIGN16( void * jointQuat1Ptr =&jointQuat1[4] );
		ALIGN16( void * jointQuat2Ptr =&jointQuat2[4] );
		ALIGN16( void * jointQuat3Ptr =&jointQuat3[4] );
		ALIGN16( void * blendQuat0Ptr =&blendQuat0[4] );
		ALIGN16( void * blendQuat1Ptr =&blendQuat1[4] );
		ALIGN16( void * blendQuat2Ptr =&blendQuat2[4] );
		ALIGN16( void * blendQuat3Ptr =&blendQuat3[4] );

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

#if 0 //todo: resolver problema dos 30 parâmetros do asm
		asm (
			// lerp translation
			"movss		xmm7, %33\n"
			"shufps		xmm7, xmm7, 0x00\n"
			"movaps		xmm0, %3\n"     /*blendVert0*/
			"subps		xmm0, %0\n"     /*jointVert0*/
			"mulps		xmm0, xmm7\n"
			"addps		xmm0, %0\n"     /*jointVert0*/
			"movaps		%0, xmm0\n"     /*jointVert0*/
			"movaps		xmm1, %4\n"     /*blendVert1*/
			"subps		xmm1, %1\n"     /*jointVert1*/
			"mulps		xmm1, xmm7\n"
			"addps		xmm1, %1\n"     /*jointVert1*/
			"movaps		%1, xmm1\n"     /*jointVert1*/
			"movaps		xmm2, %5\n"     /*blendVert2*/
			"subps		xmm2, %2\n"     /*jointVert2*/
			"mulps		xmm2, xmm7\n"
			"addps		xmm2, %2\n"     /*jointVert2*/
			"movaps		%2, xmm2\n"     /*jointVert2*/

			// lerp quaternions
			"movaps		xmm0, %6\n"      /*jointQuat0*/
			"mulps		xmm0, %10\n"     /*blendQuat0*/
			"movaps		xmm1, %7\n"       /*jointQuat1*/
			"mulps		xmm1, %11\n"     /*blendQuat1*/
			"addps		xmm0, xmm1\n"
			"movaps		xmm2, %8\n"      /*jointQuat2\*/
			"mulps		xmm2, %12\n"     /*blendQuat2*/
			"addps		xmm0, xmm2\n"
			"movaps		xmm3, %9\n"
			"mulps		xmm3, %13\n"
			"addps		xmm0, xmm3\n"					// xmm0 = cosom

			"movaps		xmm1, xmm0\n"
			"movaps		xmm2, xmm0\n"
			"andps		xmm1, %15\n"        /*SIMD_SP_signBitMask*/	// xmm1 = signBit
			"xorps		xmm0, xmm1\n"
			"mulps		xmm2, xmm2\n"

			"xorps		xmm4, xmm4\n"
			"movaps		xmm3, %14\n"    /*SIMD_SP_one*/
			"subps		xmm3, xmm2\n"					// xmm3 = scale0
			"cmpeqps	xmm4, xmm3\n"
			"andps		xmm4, %30\n"	/*SIMD_SP_tiny*/		// if values are zero replace them with a tiny number
			"andps		xmm3, %31\n"	/*SIMD_SP_absMask*/	// make sure the values are positive
			"orps		xmm3, xmm4\n"

#ifdef REFINE_BLENDJOINTS_RECIPROCAL
			"movaps		xmm2, xmm3\n"
			"rsqrtps	xmm4, xmm2\n"
			"mulps		xmm2, xmm4\n"
			"mulps		xmm2, xmm4\n"
			"subps		xmm2, %32\n"   /*SIMD_SP_rsqrt_c0*/
			"mulps		xmm4, %33\n"    /*SIMD_SP_rsqrt_c1*/
			"mulps		xmm2, xmm4\n"
#else
			"rsqrtps		xmm2, xmm3\n"					// xmm2 = sinom
#endif
			"mulps		xmm3, xmm2\n"					// xmm3 = sqrt( scale0 )

			// omega0 = atan2( xmm3, xmm0 )
			"movaps		xmm4, xmm0\n"
			"minps		xmm0, xmm3\n"
			"maxps		xmm3, xmm4\n"
			"cmpeqps	xmm4, xmm0\n"

#ifdef REFINE_BLENDJOINTS_RECIPROCAL
			"rcpps		xmm5, xmm3\n"
			"mulps		xmm3, xmm5\n"
			"mulps		xmm3, xmm5\n"
			"addps		xmm5, xmm5\n"
			"subps		xmm5, xmm3\n"					// xmm5 = 1 / y "or 1 / x
			"mulps		xmm0, xmm5\n"					// xmm0 = x / y "or y / x
#else
			"rcpps		xmm3, xmm3\n"					// xmm3 = 1 / y "or 1 / x
			"mulps		xmm0, xmm3\n"					// xmm0 = x / y "or y / x
#endif
			"movaps		xmm3, xmm4\n"
			"andps		xmm3, %15\n"        /*SIMD_SP_signBitMask*/
			"xorps		xmm0, xmm3\n"					// xmm0 = -x / y "or y / x
			"andps		xmm4, %24\n"    /*SIMD_SP_halfPI*/		// xmm4 = HALF_PI "or 0.0f
			"movaps		xmm3, xmm0\n"
			"mulps		xmm3, xmm3\n"					// xmm3 = s
			"movaps		xmm5, %16\n"        /*SIMD_SP_atan_c0*/
			"mulps		xmm5, xmm3\n"
			"addps		xmm5, %17\n"        /*SIMD_SP_atan_c1*/
			"mulps		xmm5, xmm3\n"
			"addps		xmm5, %18\n"        /*SIMD_SP_atan_c2*/
			"mulps		xmm5, xmm3\n"
			"addps		xmm5, %19\n"        /*SIMD_SP_atan_c3*/
			"mulps		xmm5, xmm3\n"
			"addps		xmm5, %20\n"        /*SIMD_SP_atan_c4*/
			"mulps		xmm5, xmm3\n"
			"addps		xmm5, %21\n"        /*SIMD_SP_atan_c5*/
			"mulps		xmm5, xmm3\n"
			"addps		xmm5, %22\n"        /*SIMD_SP_atan_c6*/
			"mulps		xmm5, xmm3\n"
			"addps		xmm5, %23\n"        /*SIMD_SP_atan_c7*/
			"mulps		xmm5, xmm3\n"
			"addps		xmm5, %14\n"    /*SIMD_SP_one*/
			"mulps		xmm5, xmm0\n"
			"addps		xmm5, xmm4\n"					// xmm5 = omega0

			"movaps		xmm6, xmm7\n"					// xmm6 = lerp
			"mulps		xmm6, xmm5\n"					// xmm6 = omega1
			"subps		xmm5, xmm6\n"					// xmm5 = omega0

			// scale0 = sin( xmm5 ) * xmm2
			// scale1 = sin( xmm6 ) * xmm2
			"movaps		xmm3, xmm5\n"
			"movaps		xmm7, xmm6\n"
			"mulps		xmm3, xmm3\n"
			"mulps		xmm7, xmm7\n"
			"movaps		xmm4, %25\n"     /*SIMD_SP_sin_c0*/
			"movaps		xmm0, %25\n"     /*SIMD_SP_sin_c0*/
			"mulps		xmm4, xmm3\n"
			"mulps		xmm0, xmm7\n"
			"addps		xmm4, %26\n"     /*SIMD_SP_sin_c1*/
			"addps		xmm0, %26\n"     /*SIMD_SP_sin_c1*/
			"mulps		xmm4, xmm3\n"
			"mulps		xmm0, xmm7\n"
			"addps		xmm4, %27\n"     /*SIMD_SP_sin_c2*/
			"addps		xmm0, %27\n"     /*SIMD_SP_sin_c2*/
			"mulps		xmm4, xmm3\n"
			"mulps		xmm0, xmm7\n"
			"addps		xmm4, %28\n"     /*SIMD_SP_sin_c3*/
			"addps		xmm0, %28\n"     /*SIMD_SP_sin_c3*/
			"mulps		xmm4, xmm3\n"
			"mulps		xmm0, xmm7\n"
			"addps		xmm4, %29\n"     /*SIMD_SP_sin_c4*/
			"addps		xmm0, %29\n"     /*SIMD_SP_sin_c4*/
			"mulps		xmm4, xmm3\n"
			"mulps		xmm0, xmm7\n"
			"addps		xmm4, %14\n"    /*SIMD_SP_one*/
			"addps		xmm0, %14\n"    /*SIMD_SP_one*/
			"mulps		xmm5, xmm4\n"
			"mulps		xmm6, xmm0\n"
			"mulps		xmm5, xmm2\n"					// xmm5 = scale0
			"mulps		xmm6, xmm2\n"					// xmm6 = scale1

			"xorps		xmm6, xmm1\n"

			"movaps		xmm0, %6\n"     /*jointQuat0*/
			"mulps		xmm0, xmm5\n"
			"movaps		xmm1, %10\n"     /*blendQuat0*/
			"mulps		xmm1, xmm6\n"
			"addps		xmm0, xmm1\n"
			"movaps		%6, xmm0\n"

			"movaps		xmm1, %7\n"     /*jointQuat1*/
			"mulps		xmm1, xmm5\n"
			"movaps		xmm2, %11\n"     /*blendQuat1*/
			"mulps		xmm2, xmm6\n"
			"addps		xmm1, xmm2\n"
			"movaps		%7, xmm1\n"     /*jointQuat1*/

			"movaps		xmm2, %8\n"     /*jointQuat2*/
			"mulps		xmm2, xmm5\n"
			"movaps		xmm3, %12\n"     /*blendQuat2*/
			"mulps		xmm3, xmm6\n"
			"addps		xmm2, xmm3\n"
			"movaps		%8, xmm2\n"     /*jointQuat2*/

			"movaps		xmm3, %9\n"     /*jointQuat3*/
			"mulps		xmm3, xmm5\n"
			"movaps		xmm4, %13\n"    /*blendQuat3*/
			"mulps		xmm4, xmm6\n"
			"addps		xmm3, xmm4\n"
			"movaps		%9, xmm3\n"
		::"m"(jointVert0),
		"m"(jointVert1),
		"m"(jointVert2),
		"m"(blendVert0),
		"m"(blendVert1),
		"m"(blendVert2),
		"m"(jointQuat0),
		"m"(jointQuat1),
		"m"(jointQuat2),
		"m"(jointQuat3),
		"m"(blendQuat0),
		"m"(blendQuat1),
		"m"(blendQuat2),
		"m"(blendQuat3),
        "m"(SIMD_SP_one),
        "m"(SIMD_SP_signBitMask),
        "m"(SIMD_SP_atan_c0),
        "m"(SIMD_SP_atan_c1),
        "m"(SIMD_SP_atan_c2),
        "m"(SIMD_SP_atan_c3),
        "m"(SIMD_SP_atan_c4),
        "m"(SIMD_SP_atan_c5),
        "m"(SIMD_SP_atan_c6),
        "m"(SIMD_SP_atan_c7),
        "m"(SIMD_SP_halfPI),   //%24
        "m"(SIMD_SP_sin_c0), //%25
        "m"(SIMD_SP_sin_c1),
        "m"(SIMD_SP_sin_c2),
        "m"(SIMD_SP_sin_c3),
        "m"(SIMD_SP_sin_c4),
        "m"(SIMD_SP_tiny),  //%30
        "m"(SIMD_SP_absMask),
        "m"(SIMD_SP_rsqrt_c0),
        "m"(lerp)
       :);

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

		omega0[0] = this->aTanPositive( scale0[0], cosom[0] );
		omega0[1] = this->aTanPositive( scale0[1], cosom[1] );
		omega0[2] = this->aTanPositive( scale0[2], cosom[2] );
		omega0[3] = this->aTanPositive( scale0[3], cosom[3] );

		omega1[0] = lerp * omega0[0];
		omega1[1] = lerp * omega0[1];
		omega1[2] = lerp * omega0[2];
		omega1[3] = lerp * omega0[3];

		omega0[0] -= omega1[0];
		omega0[1] -= omega1[1];
		omega0[2] -= omega1[2];
		omega0[3] -= omega1[3];

		scale0[0] = this->sinZeroHalfPI( omega0[0] ) * sinom[0];
		scale0[1] = this->sinZeroHalfPI( omega0[1] ) * sinom[1];
		scale0[2] = this->sinZeroHalfPI( omega0[2] ) * sinom[2];
		scale0[3] = this->sinZeroHalfPI( omega0[3] ) * sinom[3];

		scale1[0] = this->sinZeroHalfPI( omega1[0] ) * sinom[0];
		scale1[1] = this->sinZeroHalfPI( omega1[1] ) * sinom[1];
		scale1[2] = this->sinZeroHalfPI( omega1[2] ) * sinom[2];
		scale1[3] = this->sinZeroHalfPI( omega1[3] ) * sinom[3];

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

	asm (
		"mov			eax, %2\n"
		"mov			esi, %1\n"
		"mov			edi, %0\n"
		"and			eax, ~3\n"
		"jz			done4local\n"
		"imul		eax, "JOINTMAT_SIZE_STR"\n"
		"add			esi, eax\n"
		"neg			eax\n"

	"loopMat4:\n"
		"movss		xmm5, [esi+eax+3*"JOINTMAT_SIZE_STR"+0*16+0*4]\n"
		"movss		xmm6, [esi+eax+3*"JOINTMAT_SIZE_STR"+1*16+1*4]\n"
		"movss		xmm7, [esi+eax+3*"JOINTMAT_SIZE_STR"+2*16+2*4]\n"

		"shufps		xmm5, xmm5, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */
		"shufps		xmm6, xmm6, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */
		"shufps		xmm7, xmm7, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */

		"movss		xmm0, [esi+eax+2*"JOINTMAT_SIZE_STR"+0*16+0*4]\n"
		"movss		xmm1, [esi+eax+2*"JOINTMAT_SIZE_STR"+1*16+1*4]\n"
		"movss		xmm2, [esi+eax+2*"JOINTMAT_SIZE_STR"+2*16+2*4]\n"

		"movss		xmm5, xmm0\n"
		"movss		xmm6, xmm1\n"
		"movss		xmm7, xmm2\n"

		"shufps		xmm5, xmm5, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */
		"shufps		xmm6, xmm6, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */
		"shufps		xmm7, xmm7, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */

		"movss		xmm0, [esi+eax+1*"JOINTMAT_SIZE_STR"+0*16+0*4]\n"
		"movss		xmm1, [esi+eax+1*"JOINTMAT_SIZE_STR"+1*16+1*4]\n"
		"movss		xmm2, [esi+eax+1*"JOINTMAT_SIZE_STR"+2*16+2*4]\n"

		"movss		xmm5, xmm0\n"
		"movss		xmm6, xmm1\n"
		"movss		xmm7, xmm2\n"

		"shufps		xmm5, xmm5, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */
		"shufps		xmm6, xmm6, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */
		"shufps		xmm7, xmm7, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */

		"movss		xmm0, [esi+eax+0*"JOINTMAT_SIZE_STR"+0*16+0*4]\n"
		"movss		xmm1, [esi+eax+0*"JOINTMAT_SIZE_STR"+1*16+1*4]\n"
		"movss		xmm2, [esi+eax+0*"JOINTMAT_SIZE_STR"+2*16+2*4]\n"

		"movss		xmm5, xmm0\n"
		"movss		xmm6, xmm1\n"
		"movss		xmm7, xmm2\n"

		// -------------------

		"movaps		xmm0, xmm5\n"
		"addps		xmm0, xmm6\n"
		"addps		xmm0, xmm7\n"
		"cmpnltps	xmm0, %8\n"	  /*SIMD_SP_zero*/					// xmm0 = m[0 * 4 + 0] + m[1 * 4 + 1] + m[2 * 4 + 2] > 0.0f

		"movaps		xmm1, xmm5\n"
		"movaps		xmm2, xmm5\n"
		"cmpnltps	xmm1, xmm6\n"
		"cmpnltps	xmm2, xmm7\n"
		"andps		xmm2, xmm1\n"								// xmm2 = m[0 * 4 + 0] > m[1 * 4 + 1] && m[0 * 4 + 0] > m[2 * 4 + 2]

		"movaps		xmm4, xmm6\n"
		"cmpnltps	xmm4, xmm7\n"								// xmm3 = m[1 * 4 + 1] > m[2 * 4 + 2]

		"movaps		xmm1, xmm0\n"
		"andnps		xmm1, xmm2\n"
		"orps		xmm2, xmm0\n"
		"movaps		xmm3, xmm2\n"
		"andnps		xmm2, xmm4\n"
		"orps		xmm3, xmm2\n"
		"xorps		xmm3, %10\n"    /*SIMD_SP_not*/

		"andps		xmm0, %4\n"   /*SIMD_DW_mat2quatShuffle0*/
		"movaps		xmm4, xmm1\n"
		"andps		xmm4, %5\n"   /*SIMD_DW_mat2quatShuffle1*/
		"orps		xmm0, xmm4\n"
		"movaps		xmm4, xmm2\n"
		"andps		xmm4, %6\n"   /*SIMD_DW_mat2quatShuffle2*/
		"orps		xmm0, xmm4\n"
		"movaps		xmm4, xmm3\n"
		"andps		xmm4, %7\n"   /*SIMD_DW_mat2quatShuffle3*/
		"orps		xmm4, xmm0\n"

		"movaps		%13, xmm4\n"

		"movaps		xmm0, xmm2\n"
		"orps		xmm0, xmm3\n"								// xmm0 = xmm2 | xmm3	= s0
		"orps		xmm2, xmm1\n"								// xmm2 = xmm1 | xmm2	= s2
		"orps		xmm1, xmm3\n"								// xmm1 = xmm1 | xmm3	= s1

		"andps		xmm0, %3\n"
		"andps		xmm1, %3\n"
		"andps		xmm2, %3\n"

		"xorps		xmm5, xmm0\n"
		"xorps		xmm6, xmm1\n"
		"xorps		xmm7, xmm2\n"
		"addps		xmm5, xmm6\n"
		"addps		xmm7, %9\n"   /*SIMD_SP_one*/
		"addps		xmm5, xmm7\n"								// xmm5 = t

		"movaps		xmm7, xmm5\n"								// xmm7 = t
		"rsqrtps	xmm6, xmm5\n"
		"mulps		xmm5, xmm6\n"
		"mulps		xmm5, xmm6\n"
		"subps		xmm5, %11\n"    /*SIMD_SP_rsqrt_c0*/
		"mulps		xmm6, %12\n"    /*SIMD_SP_mat2quat_rsqrt_c1*/
		"mulps		xmm6, xmm5\n"								// xmm5 = s

		"mulps		xmm7, xmm6\n"								// xmm7 = s * t
		"xorps		xmm6, %3\n"			// xmm6 = -s

		// -------------------

		"add			edi, 4*"JOINTQUAT_SIZE_STR"\n"

		"movzx		ecx, byte ptr %13[0*4+0]\n"			// ecx = k0
		"movss		[edi+ecx*4-4*"JOINTQUAT_SIZE_STR"], xmm7\n"		// q[k0] = s * t;

		"movzx		edx, byte ptr %13[0*4+1]\n"			// edx = k1
		"movss		xmm4, [esi+eax+0*"JOINTMAT_SIZE_STR"+1*16+0*4]\n"
		"xorps		xmm4, xmm2\n"
		"subss		xmm4, [esi+eax+0*"JOINTMAT_SIZE_STR"+0*16+1*4]\n"
		"mulss		xmm4, xmm6\n"
		"movss		[edi+edx*4-4*"JOINTQUAT_SIZE_STR"], xmm4\n"		// q[k1] = ( m[0 * 4 + 1] - s2 * m[1 * 4 + 0] ) * s;

		"movzx		ecx, byte ptr %13[0*4+2]\n"			// ecx = k2
		"movss		xmm3, [esi+eax+0*"JOINTMAT_SIZE_STR"+0*16+2*4]\n"
		"xorps		xmm3, xmm1\n"
		"subss		xmm3, [esi+eax+0*"JOINTMAT_SIZE_STR"+2*16+0*4]\n"
		"mulss		xmm3, xmm6\n"
		"movss		[edi+ecx*4-4*"JOINTQUAT_SIZE_STR"], xmm3\n"		// q[k2] = ( m[2 * 4 + 0] - s1 * m[0 * 4 + 2] ) * s;

		"movzx		edx, byte ptr %13[0*4+3]\n"			// edx = k3
		"movss		xmm4, [esi+eax+0*"JOINTMAT_SIZE_STR"+2*16+1*4]\n"
		"xorps		xmm4, xmm0\n"
		"subss		xmm4, [esi+eax+0*"JOINTMAT_SIZE_STR"+1*16+2*4]\n"
		"mulss		xmm4, xmm6\n"
		"movss		[edi+edx*4-4*"JOINTQUAT_SIZE_STR"], xmm4\n"		// q[k3] = ( m[1 * 4 + 2] - s0 * m[2 * 4 + 1] ) * s;

		"mov			ecx, [esi+eax+0*"JOINTMAT_SIZE_STR"+0*16+3*4]\n"
		"mov			[edi-4*"JOINTQUAT_SIZE_STR"+16], ecx\n"			// q[4] = m[0 * 4 + 3];
		"mov			edx, [esi+eax+0*"JOINTMAT_SIZE_STR"+1*16+3*4]\n"
		"mov			[edi-4*"JOINTQUAT_SIZE_STR"+20], edx\n"			// q[5] = m[1 * 4 + 3];
		"mov			ecx, [esi+eax+0*"JOINTMAT_SIZE_STR"+2*16+3*4]\n"
		"mov			[edi-4*"JOINTQUAT_SIZE_STR"+24], ecx\n"			// q[6] = m[2 * 4 + 3];

		"shufps		xmm6, xmm6, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm7, xmm7, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm0, xmm0, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm1, xmm1, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm2, xmm2, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/

		"movzx		ecx, byte ptr %13[1*4+0]\n"			// ecx = k0
		"movss		[edi+ecx*4-3*"JOINTQUAT_SIZE_STR"], xmm7\n"		// q[k0] = s * t;

		"movzx		edx, byte ptr %13[1*4+1]\n"			// edx = k1
		"movss		xmm4, [esi+eax+1*"JOINTMAT_SIZE_STR"+1*16+0*4]\n"
		"xorps		xmm4, xmm2\n"
		"subss		xmm4, [esi+eax+1*"JOINTMAT_SIZE_STR"+0*16+1*4]\n"
		"mulss		xmm4, xmm6\n"
		"movss		[edi+edx*4-3*"JOINTQUAT_SIZE_STR"], xmm4\n"		// q[k1] = ( m[0 * 4 + 1] - s2 * m[1 * 4 + 0] ) * s;

		"movzx		ecx, byte ptr %13[1*4+2]\n"			// ecx = k2
		"movss		xmm3, [esi+eax+1*"JOINTMAT_SIZE_STR"+0*16+2*4]\n"
		"xorps		xmm3, xmm1\n"
		"subss		xmm3, [esi+eax+1*"JOINTMAT_SIZE_STR"+2*16+0*4]\n"
		"mulss		xmm3, xmm6\n"
		"movss		[edi+ecx*4-3*"JOINTQUAT_SIZE_STR"], xmm3\n"		// q[k2] = ( m[2 * 4 + 0] - s1 * m[0 * 4 + 2] ) * s;

		"movzx		edx, byte ptr %13[1*4+3]\n"			// edx = k3
		"movss		xmm4, [esi+eax+1*"JOINTMAT_SIZE_STR"+2*16+1*4]\n"
		"xorps		xmm4, xmm0\n"
		"subss		xmm4, [esi+eax+1*"JOINTMAT_SIZE_STR"+1*16+2*4]\n"
		"mulss		xmm4, xmm6\n"
		"movss		[edi+edx*4-3*"JOINTQUAT_SIZE_STR"], xmm4\n"		// q[k3] = ( m[1 * 4 + 2] - s0 * m[2 * 4 + 1] ) * s;

		"mov			ecx, [esi+eax+1*"JOINTMAT_SIZE_STR"+0*16+3*4]\n"
		"mov			[edi-3*"JOINTQUAT_SIZE_STR"+16], ecx\n"			// q[4] = m[0 * 4 + 3];
		"mov			edx, [esi+eax+1*"JOINTMAT_SIZE_STR"+1*16+3*4]\n"
		"mov			[edi-3*"JOINTQUAT_SIZE_STR"+20], edx\n"			// q[5] = m[1 * 4 + 3];
		"mov			ecx, [esi+eax+1*"JOINTMAT_SIZE_STR"+2*16+3*4]\n"
		"mov			[edi-3*"JOINTQUAT_SIZE_STR"+24], ecx\n"			// q[6] = m[2 * 4 + 3];

		"shufps		xmm6, xmm6, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm7, xmm7, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm0, xmm0, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm1, xmm1, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm2, xmm2, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/

		"movzx		ecx, byte ptr %13[2*4+0]\n"			// ecx = k0
		"movss		[edi+ecx*4-2*"JOINTQUAT_SIZE_STR"], xmm7\n"		// q[k0] = s * t;

		"movzx		edx, byte ptr %13[2*4+1]\n"			// edx = k1
		"movss		xmm4, [esi+eax+2*"JOINTMAT_SIZE_STR"+1*16+0*4]\n"
		"xorps		xmm4, xmm2\n"
		"subss		xmm4, [esi+eax+2*"JOINTMAT_SIZE_STR"+0*16+1*4]\n"
		"mulss		xmm4, xmm6\n"
		"movss		[edi+edx*4-2*"JOINTQUAT_SIZE_STR"], xmm4\n"		// q[k1] = ( m[0 * 4 + 1] - s2 * m[1 * 4 + 0] ) * s;

		"movzx		ecx, byte ptr %13[2*4+2]\n"			// ecx = k2
		"movss		xmm3, [esi+eax+2*"JOINTMAT_SIZE_STR"+0*16+2*4]\n"
		"xorps		xmm3, xmm1\n"
		"subss		xmm3, [esi+eax+2*"JOINTMAT_SIZE_STR"+2*16+0*4]\n"
		"mulss		xmm3, xmm6\n"
		"movss		[edi+ecx*4-2*"JOINTQUAT_SIZE_STR"], xmm3\n"		// q[k2] = ( m[2 * 4 + 0] - s1 * m[0 * 4 + 2] ) * s;

		"movzx		edx, byte ptr %13[2*4+3]\n"			// edx = k3
		"movss		xmm4, [esi+eax+2*"JOINTMAT_SIZE_STR"+2*16+1*4]\n"
		"xorps		xmm4, xmm0\n"
		"subss		xmm4, [esi+eax+2*"JOINTMAT_SIZE_STR"+1*16+2*4]\n"
		"mulss		xmm4, xmm6\n"
		"movss		[edi+edx*4-2*"JOINTQUAT_SIZE_STR"], xmm4\n"		// q[k3] = ( m[1 * 4 + 2] - s0 * m[2 * 4 + 1] ) * s;

		"mov			ecx, [esi+eax+2*"JOINTMAT_SIZE_STR"+0*16+3*4]\n"
		"mov			[edi-2*"JOINTQUAT_SIZE_STR"+16], ecx\n"			// q[4] = m[0 * 4 + 3];
		"mov			edx, [esi+eax+2*"JOINTMAT_SIZE_STR"+1*16+3*4]\n"
		"mov			[edi-2*"JOINTQUAT_SIZE_STR"+20], edx\n"			// q[5] = m[1 * 4 + 3];
		"mov			ecx, [esi+eax+2*"JOINTMAT_SIZE_STR"+2*16+3*4]\n"
		"mov			[edi-2*"JOINTQUAT_SIZE_STR"+24], ecx\n"			// q[6] = m[2 * 4 + 3];

		"shufps		xmm6, xmm6, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm7, xmm7, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm0, xmm0, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm1, xmm1, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm2, xmm2, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/

		"movzx		ecx, byte ptr %13[3*4+0]\n"			// ecx = k0
		"movss		[edi+ecx*4-1*"JOINTQUAT_SIZE_STR"], xmm7\n"		// q[k0] = s * t;

		"movzx		edx, byte ptr %13[3*4+1]\n"			// edx = k1
		"movss		xmm4, [esi+eax+3*"JOINTMAT_SIZE_STR"+1*16+0*4]\n"
		"xorps		xmm4, xmm2\n"
		"subss		xmm4, [esi+eax+3*"JOINTMAT_SIZE_STR"+0*16+1*4]\n"
		"mulss		xmm4, xmm6\n"
		"movss		[edi+edx*4-1*"JOINTQUAT_SIZE_STR"], xmm4\n"		// q[k1] = ( m[0 * 4 + 1] - s2 * m[1 * 4 + 0] ) * s;

		"movzx		ecx, byte ptr %13[3*4+2]\n"			// ecx = k2
		"movss		xmm3, [esi+eax+3*"JOINTMAT_SIZE_STR"+0*16+2*4]\n"
		"xorps		xmm3, xmm1\n"
		"subss		xmm3, [esi+eax+3*"JOINTMAT_SIZE_STR"+2*16+0*4]\n"
		"mulss		xmm3, xmm6\n"
		"movss		[edi+ecx*4-1*"JOINTQUAT_SIZE_STR"], xmm3\n"		// q[k2] = ( m[2 * 4 + 0] - s1 * m[0 * 4 + 2] ) * s;

		"movzx		edx, byte ptr %13[3*4+3]\n"			// edx = k3
		"movss		xmm4, [esi+eax+3*"JOINTMAT_SIZE_STR"+2*16+1*4]\n"
		"xorps		xmm4, xmm0\n"
		"subss		xmm4, [esi+eax+3*"JOINTMAT_SIZE_STR"+1*16+2*4]\n"
		"mulss		xmm4, xmm6\n"
		"movss		[edi+edx*4-1*"JOINTQUAT_SIZE_STR"], xmm4\n"		// q[k3] = ( m[1 * 4 + 2] - s0 * m[2 * 4 + 1] ) * s;

		"mov			ecx, [esi+eax+3*"JOINTMAT_SIZE_STR"+0*16+3*4]\n"
		"mov			[edi-1*"JOINTQUAT_SIZE_STR"+16], ecx\n"			// q[4] = m[0 * 4 + 3];
		"mov			edx, [esi+eax+3*"JOINTMAT_SIZE_STR"+1*16+3*4]\n"
		"mov			[edi-1*"JOINTQUAT_SIZE_STR"+20], edx\n"			// q[5] = m[1 * 4 + 3];
		"mov			ecx, [esi+eax+3*"JOINTMAT_SIZE_STR"+2*16+3*4]\n"
		"mov			[edi-1*"JOINTQUAT_SIZE_STR"+24], ecx\n"			// q[6] = m[2 * 4 + 3];

		"add			eax, 4*"JOINTMAT_SIZE_STR"\n"
		"jl			loopMat4\n"

	"done4local:\n"
		"mov			eax, %2\n"
		"and			eax, 3\n"
		"jz			done1\n"
		"imul		eax, "JOINTMAT_SIZE_STR"\n"
		"add			esi, eax\n"
		"neg			eax\n"

	"loopMat1:\n"
		"movss		xmm5, [esi+eax+0*"JOINTMAT_SIZE_STR"+0*16+0*4]\n"
		"movss		xmm6, [esi+eax+0*"JOINTMAT_SIZE_STR"+1*16+1*4]\n"
		"movss		xmm7, [esi+eax+0*"JOINTMAT_SIZE_STR"+2*16+2*4]\n"

		// -------------------

		"movaps		xmm0, xmm5\n"
		"addss		xmm0, xmm6\n"
		"addss		xmm0, xmm7\n"
		"cmpnltss	xmm0, %8	\n"					// xmm0 = m[0 * 4 + 0] + m[1 * 4 + 1] + m[2 * 4 + 2] > 0.0f

		"movaps		xmm1, xmm5\n"
		"movaps		xmm2, xmm5\n"
		"cmpnltss	xmm1, xmm6\n"
		"cmpnltss	xmm2, xmm7\n"
		"andps		xmm2, xmm1\n"								// xmm2 = m[0 * 4 + 0] > m[1 * 4 + 1] && m[0 * 4 + 0] > m[2 * 4 + 2]

		"movaps		xmm4, xmm6\n"
		"cmpnltss	xmm4, xmm7\n"								// xmm3 = m[1 * 4 + 1] > m[2 * 4 + 2]

		"movaps		xmm1, xmm0\n"
		"andnps		xmm1, xmm2\n"
		"orps		xmm2, xmm0\n"
		"movaps		xmm3, xmm2\n"
		"andnps		xmm2, xmm4\n"
		"orps		xmm3, xmm2\n"
		"xorps		xmm3, %10\n"  /*SIMD_SP_not*/

		"andps		xmm0, %4\n"  /*SIMD_DW_mat2quatShuffle0*/
		"movaps		xmm4, xmm1\n"
		"andps		xmm4, %5\n"
		"orps		xmm0, xmm4\n"
		"movaps		xmm4, xmm2\n"
		"andps		xmm4, %6\n"
		"orps		xmm0, xmm4\n"
		"movaps		xmm4, xmm3\n"
		"andps		xmm4, %7\n"
		"orps		xmm4, xmm0\n"

		"movss		%13, xmm4\n"

		"movaps		xmm0, xmm2\n"
		"orps		xmm0, xmm3\n"								// xmm0 = xmm2 | xmm3	= s0
		"orps		xmm2, xmm1\n"								// xmm2 = xmm1 | xmm2	= s2
		"orps		xmm1, xmm3\n"								// xmm1 = xmm1 | xmm3	= s1

		"andps		xmm0, %3\n"  /*SIMD_SP_signBitMask*/
		"andps		xmm1, %3\n"
		"andps		xmm2, %3\n"

		"xorps		xmm5, xmm0\n"
		"xorps		xmm6, xmm1\n"
		"xorps		xmm7, xmm2\n"
		"addss		xmm5, xmm6\n"
		"addss		xmm7, %9\n" /*SIMD_SP_one*/
		"addss		xmm5, xmm7\n"								// xmm5 = t

		"movss		xmm7, xmm5\n"								// xmm7 = t
		"rsqrtss		xmm6, xmm5\n"
		"mulss		xmm5, xmm6\n"
		"mulss		xmm5, xmm6\n"
		"subss		xmm5, %11\n"    /*SIMD_SP_rsqrt_c0*/
		"mulss		xmm6, %12\n"    /*SIMD_SP_mat2quat_rsqrt_c1*/
		"mulss		xmm6, xmm5\n"								// xmm5 = s

		"mulss		xmm7, xmm6\n"								// xmm7 = s * t
		"xorps		xmm6, %3\n"  /*SIMD_SP_signBitMask*/				// xmm6 = -s

		// -------------------

		"movzx		ecx, byte ptr %13[0]\n"				// ecx = k0
		"add			edi, "JOINTQUAT_SIZE_STR"\n"
		"movss		[edi+ecx*4-1*"JOINTQUAT_SIZE_STR"], xmm7\n"		// q[k0] = s * t;

		"movzx		edx, byte ptr %13[1]\n"			// edx = k1
		"movss		xmm4, [esi+eax+0*"JOINTMAT_SIZE_STR"+1*16+0*4]\n"
		"xorps		xmm4, xmm2\n"
		"subss		xmm4, [esi+eax+0*"JOINTMAT_SIZE_STR"+0*16+1*4]\n"
		"mulss		xmm4, xmm6\n"
		"movss		[edi+edx*4-1*"JOINTQUAT_SIZE_STR"], xmm4\n"		// q[k1] = ( m[0 * 4 + 1] - s2 * m[1 * 4 + 0] ) * s;

		"movzx		ecx, byte ptr %13[2]\n"				// ecx = k2
		"movss		xmm3, [esi+eax+0*"JOINTMAT_SIZE_STR"+0*16+2*4]\n"
		"xorps		xmm3, xmm1\n"
		"subss		xmm3, [esi+eax+0*"JOINTMAT_SIZE_STR"+2*16+0*4]\n"
		"mulss		xmm3, xmm6\n"
		"movss		[edi+ecx*4-1*"JOINTQUAT_SIZE_STR"], xmm3\n"		// q[k2] = ( m[2 * 4 + 0] - s1 * m[0 * 4 + 2] ) * s;

		"movzx		edx, byte ptr %13[3]\n"				// edx = k3
		"movss		xmm4, [esi+eax+0*"JOINTMAT_SIZE_STR"+2*16+1*4]\n"
		"xorps		xmm4, xmm0\n"
		"subss		xmm4, [esi+eax+0*"JOINTMAT_SIZE_STR"+1*16+2*4]\n"
		"mulss		xmm4, xmm6\n"
		"movss		[edi+edx*4-1*"JOINTQUAT_SIZE_STR"], xmm4\n"		// q[k3] = ( m[1 * 4 + 2] - s0 * m[2 * 4 + 1] ) * s;

		"mov			ecx, [esi+eax+0*"JOINTMAT_SIZE_STR"+0*16+3*4]\n"
		"mov			[edi-1*"JOINTQUAT_SIZE_STR"+16], ecx\n"			// q[4] = m[0 * 4 + 3];
		"mov			edx, [esi+eax+0*"JOINTMAT_SIZE_STR"+1*16+3*4]\n"
		"mov			[edi-1*"JOINTQUAT_SIZE_STR"+20], edx\n"			// q[5] = m[1 * 4 + 3];
		"mov			ecx, [esi+eax+0*"JOINTMAT_SIZE_STR"+2*16+3*4]\n"
		"mov			[edi-1*"JOINTQUAT_SIZE_STR"+24], ecx\n"			// q[6] = m[2 * 4 + 3];

		"add			eax, "JOINTMAT_SIZE_STR"\n"
		"jl			loopMat1\n"

	"done1:\n"
	:: "m"(jointQuats), "m"(jointMats), "m"(numJoints),"m"(SIMD_SP_signBitMask),"m"(SIMD_DW_mat2quatShuffle0),
	"m"(SIMD_DW_mat2quatShuffle1),"m"(SIMD_DW_mat2quatShuffle2),"m"(SIMD_DW_mat2quatShuffle3),
	"m"(SIMD_SP_zero),"m"(SIMD_SP_one), "m"(SIMD_SP_not),"m"(SIMD_SP_rsqrt_c0),"m"(SIMD_SP_mat2quat_rsqrt_c1),"m"(shuffle)
    :"eax","ebx","ecx","edx","esi","edi");

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

	asm (

		"mov		ecx, %2\n"
		"mov		eax, %3\n"
		"sub		eax, ecx\n"
		"jl			done%=\n"
		"imul		ecx, 4\n"
		"mov		edi, %1\n"
		"add		edi, ecx\n"
		"imul		ecx, 12\n"
		"mov		esi, %0\n"
		"imul		eax, 4\n"
		"add		edi, eax\n"
		"neg		eax\n"

	"loopJoint%=:\n"

		"movaps		xmm0, [esi+ecx+ 0]\n"						// xmm0 = m0, m1, m2, t0
		"mov		edx, [edi+eax]\n"
		"movaps		xmm1, [esi+ecx+16]	\n"					// xmm1 = m2, m3, m4, t1
		"imul		edx, "JOINTMAT_SIZE_STR"\n"
		"movaps		xmm2, [esi+ecx+32]	\n"					// xmm2 = m5, m6, m7, t2

		"movss		xmm4, [esi+edx+ 0]\n"
		"shufps		xmm4, xmm4, 0x00\n"
		"mulps		xmm4, xmm0\n"

		"movss		xmm5, [esi+edx+ 4]\n"
		"shufps		xmm5, xmm5, 0x00\n"
		"mulps		xmm5, xmm1\n"
		"addps		xmm4, xmm5\n"
		"movss		xmm6, [esi+edx+ 8]\n"
		"shufps		xmm6, xmm6, 0x00\n"
		"mulps		xmm6, xmm2\n"
		"addps		xmm4, xmm6\n"

		"movss		xmm5, [esi+edx+16]\n"
		"shufps		xmm5, xmm5, 0x00\n"
		"mulps		xmm5, xmm0\n"

		"movss		xmm7, [esi+edx+12]\n"
		"shufps		xmm7, xmm7, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"addps		xmm4, xmm7\n"

		"movaps		[esi+ecx+ 0], xmm4\n"

		"movss		xmm6, [esi+edx+20]\n"
		"shufps		xmm6, xmm6, 0x00\n"
		"mulps		xmm6, xmm1\n"
		"addps		xmm5, xmm6\n"
		"movss		xmm7, [esi+edx+24]\n"
		"shufps		xmm7, xmm7, 0x00\n"
		"mulps		xmm7, xmm2\n"
		"addps		xmm5, xmm7\n"

		"movss		xmm6, [esi+edx+32]\n"
		"shufps		xmm6, xmm6, 0x00\n"
		"mulps		xmm6, xmm0\n"

		"movss		xmm3, [esi+edx+28]\n"
		"shufps		xmm3, xmm3, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"addps		xmm5, xmm3\n"

		"movaps		[esi+ecx+16], xmm5\n"

		"movss		xmm7, [esi+edx+36]\n"
		"shufps		xmm7, xmm7, 0x00\n"
		"mulps		xmm7, xmm1\n"
		"addps		xmm6, xmm7\n"
		"movss		xmm3, [esi+edx+40]\n"
		"shufps		xmm3, xmm3, 0x00\n"
		"mulps		xmm3, xmm2\n"
		"addps		xmm6, xmm3\n"

		"movss		xmm7, [esi+edx+44]\n"
		"shufps		xmm7, xmm7, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"addps		xmm6, xmm7\n"

		"movaps		[esi+ecx+32], xmm6\n"

		"add		ecx, "JOINTMAT_SIZE_STR"\n"
		"add		eax, 4\n"
		"jle		loopJoint%=\n"
	"done%=:\n"
	::"m"(jointMats), "m"(parents), "m"(firstJoint), "m"(lastJoint)
    :"eax","ebx","ecx","edx","esi","edi");

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

	asm (

		"mov			edx, %2\n"
		"mov			eax, %3\n"
		"mov			ecx, eax\n"
		"sub			eax, edx\n"
		"jl			done%=\n"
		"mov			esi, %0\n"
		"imul		ecx, "JOINTMAT_SIZE_STR"\n"
		"imul		edx, 4\n"
		"mov			edi, %1\n"
		"add			edi, edx\n"
		"imul		eax, 4\n"

	"loopJoint%=:\n"

		"movaps		xmm0, [esi+ecx+ 0]\n"						// xmm0 = m0, m1, m2, t0
		"mov		edx, [edi+eax]\n"
		"movaps		xmm1, [esi+ecx+16]\n"						// xmm1 = m2, m3, m4, t1
		"imul		edx, "JOINTMAT_SIZE_STR"\n"
		"movaps		xmm2, [esi+ecx+32]\n"						// xmm2 = m5, m6, m7, t2

		"movss		xmm6, [esi+edx+12]\n"
		"shufps		xmm6, xmm6, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"subps		xmm0, xmm6\n"
		"movss		xmm7, [esi+edx+28]\n"
		"shufps		xmm7, xmm7, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"subps		xmm1, xmm7\n"
		"movss		xmm3, [esi+edx+44]\n"
		"shufps		xmm3, xmm3, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"subps		xmm2, xmm3\n"

		"movss		xmm4, [esi+edx+ 0]\n"
		"shufps		xmm4, xmm4, 0x00\n"
		"mulps		xmm4, xmm0\n"
		"movss		xmm5, [esi+edx+16]\n"
		"shufps		xmm5, xmm5, 0x00\n"
		"mulps		xmm5, xmm1\n"
		"addps		xmm4, xmm5\n"
		"movss		xmm6, [esi+edx+32]\n"
		"shufps		xmm6, xmm6, 0x00\n"
		"mulps		xmm6, xmm2\n"
		"addps		xmm4, xmm6\n"

		"movaps		[esi+ecx+ 0], xmm4\n"

		"movss		xmm5, [esi+edx+ 4]\n"
		"shufps		xmm5, xmm5, 0x00\n"
		"mulps		xmm5, xmm0\n"
		"movss		xmm6, [esi+edx+20]\n"
		"shufps		xmm6, xmm6, 0x00\n"
		"mulps		xmm6, xmm1\n"
		"addps		xmm5, xmm6\n"
		"movss		xmm7, [esi+edx+36]\n"
		"shufps		xmm7, xmm7, 0x00\n"
		"mulps		xmm7, xmm2\n"
		"addps		xmm5, xmm7\n"

		"movaps		[esi+ecx+16], xmm5\n"

		"movss		xmm6, [esi+edx+ 8]\n"
		"shufps		xmm6, xmm6, 0x00\n"
		"mulps		xmm6, xmm0\n"
		"movss		xmm7, [esi+edx+24]\n"
		"shufps		xmm7, xmm7, 0x00\n"
		"mulps		xmm7, xmm1\n"
		"addps		xmm6, xmm7\n"
		"movss		xmm3, [esi+edx+40]\n"
		"shufps		xmm3, xmm3, 0x00\n"
		"mulps		xmm3, xmm2\n"
		"addps		xmm6, xmm3\n"

		"movaps		[esi+ecx+32], xmm6\n"

		"sub			ecx, "JOINTMAT_SIZE_STR"\n"
		"sub			eax, 4\n"
		"jge			loopJoint%=\n"
	"done%=:\n"
	::"m"(jointMats), "m"(parents), "m"(firstJoint), "m"(lastJoint)
    :"eax","ebx","ecx","edx","esi","edi");

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

	asm(
		"mov		eax, %1\n"
		"test		eax, eax\n"
		"jz			done%=\n"
		"imul		eax, "DRAWVERT_SIZE_STR"\n"

		"mov			ecx, %0\n"
		"mov			edx, %4\n"
		"mov			esi, %3\n"
		"mov			edi, %2\n"

		"add			ecx, eax\n"
		"neg			eax\n"

	"loopVert%=:\n"
		"mov			ebx, [edx]\n"
		"movaps		xmm2, [esi]\n"
		"add			edx, 8\n"
		"movaps		xmm0, xmm2\n"
		"add			esi, "JOINTWEIGHT_SIZE_STR"\n"
		"movaps		xmm1, xmm2\n"

		"mulps		xmm0, [edi+ebx+ 0]\n"						// xmm0 = m0, m1, m2, t0
		"mulps		xmm1, [edi+ebx+16]\n"						// xmm1 = m3, m4, m5, t1
		"mulps		xmm2, [edi+ebx+32]\n"						// xmm2 = m6, m7, m8, t2

		"cmp			dword ptr [edx-4], 0\n"

		"jne			doneWeight%=\n"

	"loopWeight%=:\n"
		"mov			ebx, [edx]\n"
		"movaps		xmm5, [esi]\n"
		"add			edx, 8\n"
		"movaps		xmm3, xmm5\n"
		"add			esi, "JOINTWEIGHT_SIZE_STR"\n"
		"movaps		xmm4, xmm5\n"

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

		"movaps		xmm6, xmm0\n"								// xmm6 =    m0,    m1,          m2,          t0
		"unpcklps	xmm6, xmm1\n"								// xmm6 =    m0,    m3,          m1,          m4
		"unpckhps	xmm0, xmm1\n"								// xmm1 =    m2,    m5,          t0,          t1
		"addps		xmm6, xmm0\n"								// xmm6 = m0+m2, m3+m5,       m1+t0,       m4+t1

		"movaps		xmm7, xmm2\n"								// xmm7 =    m6,    m7,          m8,          t2
		"movlhps		xmm2, xmm6\n"								// xmm2 =    m6,    m7,       m0+m2,       m3+m5
		"movhlps		xmm6, xmm7\n"								// xmm6 =    m8,    t2,       m1+t0,       m4+t1
		"addps		xmm6, xmm2\n"								// xmm6 = m6+m8, m7+t2, m0+m1+m2+t0, m3+m4+m5+t1

		"movhps		[ecx+eax-"DRAWVERT_SIZE_STR"+0], xmm6\n"

		"movaps		xmm5, xmm6\n"								// xmm5 = m6+m8, m7+t2
		"shufps		xmm5, xmm5, 0xE1\n"  /*R_SHUFFLEPS( 1, 0, 2, 3 )=11100001 */	// xmm5 = m7+t2, m6+m8
		"addss		xmm5, xmm6\n"								// xmm5 = m6+m8+m7+t2

		"movss		[ecx+eax-"DRAWVERT_SIZE_STR"+8], xmm5\n"

		"jl			loopVert%=\n"
	"done%=:\n"
	:
    :"m"(verts), "m"(numVerts), "m"(joints), "m"(weights), "m"(index), "m"(numWeights)
    :"eax","ebx","ecx","edx","esi","edi");

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
void VPCALL CSIMD_SSE::tracePointCull( sf_u8 *cullBits, sf_u8 &totalOr2, const float radius, const CPlane *planes, const CVertex *verts, const int numVerts ) {
#if 1
    void *totalOr = &totalOr2;
	SMF_ASSERT( sizeof( CVertex ) == DRAWVERT_SIZE );
	SMF_ASSERT( (int)&((CVertex *)0)->xyz == DRAWVERT_XYZ_OFFSET );

	asm (
		"push		ebx\n"
		"mov		eax, %5\n"
		"test		eax, eax\n"
		"jz			done%=\n"

		"mov		edi, %3\n"
		"movlps		xmm1, [edi]\n"							// xmm1 =  0,  1,  X,  X
		"movhps		xmm1, [edi+16]\n"							// xmm1 =  0,  1,  4,  5
		"movlps		xmm3, [edi+8]\n"							// xmm3 =  2,  3,  X,  X
		"movhps		xmm3, [edi+24]\n"							// xmm3 =  2,  3,  6,  7
		"movlps		xmm4, [edi+32]\n"							// xmm4 =  8,  9,  X,  X
		"movhps		xmm4, [edi+48]\n"							// xmm4 =  8,  9, 12, 13
		"movlps		xmm5, [edi+40]\n"							// xmm5 = 10, 11,  X,  X
		"movhps		xmm5, [edi+56]\n"							// xmm5 = 10, 11, 14, 15
		"movaps		xmm0, xmm1\n"								// xmm0 =  0,  1,  4,  5
		"shufps		xmm0, xmm4, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/	// xmm0 =  0,  4,  8, 12
		"shufps		xmm1, xmm4, 0xDD\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )=11011101*/	// xmm1 =  1,  5,  9, 13
		"movaps		xmm2, xmm3\n"								// xmm2 =  2,  3,  6,  7
		"shufps		xmm2, xmm5, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/	// xmm2 =  2,  6, 10, 14
		"shufps		xmm3, xmm5, 0xDD\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )=11011101*/	// xmm3 =  3,  7, 11, 15
		"movss		xmm7, %2\n"
		"shufps		xmm7, xmm7, 0x00\n"

		"xor		edx, edx\n"
		"mov		esi, %4\n"
		"mov		edi, %0\n"
		"imul		eax, "DRAWVERT_SIZE_STR"\n"
		"add		esi, eax\n"
		"neg		eax\n"

	"loopVert%=:\n"
		"movss		xmm4, [esi+eax+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"shufps		xmm4, xmm4, 0x00\n"
		"movss		xmm5, [esi+eax+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"mulps		xmm4, xmm0\n"
		"shufps		xmm5, xmm5, 0x00\n"
		"movss		xmm6, [esi+eax+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"mulps		xmm5, xmm1\n"
		"shufps		xmm6, xmm6, 0x00\n"
		"addps		xmm4, xmm5\n"
		"mulps		xmm6, xmm2\n"
		"addps		xmm4, xmm3\n"
		"addps		xmm4, xmm6\n"
		"movaps		xmm5, xmm4\n"
		"xorps		xmm5, %6\n"
		"cmpltps	xmm4, xmm7\n"
		"movmskps	ecx, xmm4\n"
		"cmpltps	xmm5, xmm7\n"
		"movmskps	ebx, xmm5\n"
		"shl		cx, 4\n"
		"or			cl, bl\n"
		"inc		edi\n"
		"or			dl, cl\n"
		"add		eax, "DRAWVERT_SIZE_STR"\n"
		"mov		byte ptr [edi-1], cl\n"
		"jl			loopVert%=\n"

	"done%=:\n"
		"mov		esi, %1\n"
        "mov		byte ptr [esi], dl\n"
		"pop		ebx\n"
	:
    : "m"(cullBits), "m"(totalOr), "m"(radius), "m"(planes), "m"(verts), "m"(numVerts),"m"(SIMD_SP_signBitMask)
    :"eax","ebx","ecx","edx","esi","edi");

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
    asm (
		"mov		ecx, %1\n"   //planes
		"movlps		xmm1, [ecx]\n"								// xmm1 =  0,  1,  X,  X
		"movhps		xmm1, [ecx+16]\n"							// xmm1 =  0,  1,  4,  5
		"movlps		xmm3, [ecx+8]\n"							// xmm3 =  2,  3,  X,  X
		"movhps		xmm3, [ecx+24]\n"							// xmm3 =  2,  3,  6,  7
		"movlps		xmm4, [ecx+32]\n"							// xmm4 =  8,  9,  X,  X
		"movhps		xmm4, [ecx+48]\n"							// xmm4 =  8,  9, 12, 13
		"movlps		xmm5, [ecx+40]\n"							// xmm5 = 10, 11,  X,  X
		"movhps		xmm5, [ecx+56]\n"							// xmm5 = 10, 11, 14, 15
		"movaps		xmm0, xmm1\n"								// xmm0 =  0,  1,  4,  5
		"shufps		xmm0, xmm4, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/	// xmm0 =  0,  4,  8, 12
		"shufps		xmm1, xmm4, 0xDD\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )=11011101*/	// xmm1 =  1,  5,  9, 13
		"movaps		xmm2, xmm3\n"								// xmm2 =  2,  3,  6,  7
		"shufps		xmm2, xmm5, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/	// xmm2 =  2,  6, 10, 14
		"shufps		xmm3, xmm5, 0xDD\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )=11011101*/	// xmm3 =  3,  7, 11, 15

		"movaps		%5, xmm0\n"
		"movaps		%6, xmm1\n"
		"movaps		%7, xmm2\n"
		"movaps		%8, xmm3\n"

		"movlps		xmm4, [ecx+64]\n"							// xmm4 = p40, p41,   X,   X
		"movhps		xmm4, [ecx+80]\n"							// xmm4 = p40, p41, p50, p51
		"movaps		xmm5, xmm4\n"								// xmm5 = p40, p41, p50, p51
		"shufps		xmm4, xmm4, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/	// xmm4 = p40, p50, p40, p50
		"shufps		xmm5, xmm5, 0xDD\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )=11011101*/	// xmm5 = p41, p51, p41, p51
		"movlps		xmm6, [ecx+72]\n"							// xmm6 = p42, p43,   X,   X
		"movhps		xmm6, [ecx+88]\n"							// xmm6 = p42, p43, p52, p53
		"movaps		xmm7, xmm6\n"								// xmm7 = p42, p43, p52, p53
		"shufps		xmm6, xmm6, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/	// xmm6 = p42, p52, p42, p52
		"shufps		xmm7, xmm7, 0xDD\n" /*R_SHUFFLEPS( 1, 3, 1, 3 )=11011101*/	// xmm7 = p43, p53, p43, p53

		"movaps		%9, xmm4\n"
		"movaps		%10, xmm5\n"
		"movaps		%11, xmm6\n"
		"movaps		%12, xmm7\n"

		"mov		esi, %2\n"   //verts
		"mov		edi, %0\n"
		"mov		eax, %3\n"  //numverts
		"and		eax, ~1\n"
		"jz			done2%=\n"
		"imul		eax, "DRAWVERT_SIZE_STR"\n"
		"add		esi, eax\n"
		"neg		eax\n"

	"loopVert2%=:\n"
		"movaps		xmm6, %5\n"
		"movss		xmm0, [esi+eax+0*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"shufps		xmm0, xmm0, 0x00\n"
		"mulps		xmm6, xmm0\n"
		"movaps		xmm7, %6\n"
		"movss		xmm1, [esi+eax+0*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"shufps		xmm1, xmm1, 0x00\n"
		"mulps		xmm7, xmm1\n"
		"addps		xmm6, xmm7\n"
		"movaps		xmm7, %7\n"
		"movss		xmm2, [esi+eax+0*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"shufps		xmm2, xmm2, 0x00\n"
		"mulps		xmm7, xmm2\n"
		"addps		xmm6, xmm7\n"
		"addps		xmm6, %8\n"

		"cmpnltps	xmm6, %4\n"
		"movmskps	ecx, xmm6\n"

		"movaps		xmm6, %5\n"
		"movss		xmm3, [esi+eax+1*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"shufps		xmm3, xmm3, 0x00\n"
		"mulps		xmm6, xmm3\n"
		"movaps		xmm7, %6\n"
		"movss		xmm4, [esi+eax+1*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"shufps		xmm4, xmm4, 0x00\n"
		"mulps		xmm7, xmm4\n"
		"addps		xmm6, xmm7\n"
		"movaps		xmm7, %7\n"
		"movss		xmm5, [esi+eax+1*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"shufps		xmm5, xmm5, 0x00\n"
		"mulps		xmm7, xmm5\n"
		"addps		xmm6, xmm7\n"
		"addps		xmm6, %8\n"

		"cmpnltps	xmm6, %4\n"
		"movmskps	edx, xmm6\n"
		"mov		ch, dl\n"

		"shufps		xmm0, xmm3, 0x00\n"
		"mulps		xmm0, %9\n"
		"shufps		xmm1, xmm4, 0x00\n"
		"mulps		xmm1, %10\n"
		"addps		xmm0, xmm1\n"
		"shufps		xmm2, xmm5, 0x00\n"
		"mulps		xmm2, %11\n"
		"addps		xmm0, xmm2\n"
		"addps		xmm0, %12\n"

		"cmpnltps	xmm0, %4\n"
		"movmskps	edx, xmm0\n"

		"add			edi, 2\n"

		"mov			dh, dl\n"
		"shl			dl, 4\n"
		"shl			dh, 2\n"
		"and			edx, (3<<4)|(3<<12)\n"
		"or			ecx, edx\n"

		"add			eax, 2*"DRAWVERT_SIZE_STR"\n"
		"mov			word ptr [edi-2], cx\n"
		"jl			loopVert2%=\n"

	"done2%=:\n"

		"mov			eax, %3\n"   //numverts
		"and			eax, 1\n"
		"jz			done%=\n"

		"movaps		xmm6, %5\n"
		"movss		xmm0, [esi+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"shufps		xmm0, xmm0, 0x00\n"
		"mulps		xmm6, xmm0\n"
		"movaps		xmm7, %6\n"
		"movss		xmm1, [esi+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"shufps		xmm1, xmm1, 0x00\n"
		"mulps		xmm7, xmm1\n"
		"addps		xmm6, xmm7\n"
		"movaps		xmm7, %7\n"
		"movss		xmm2, [esi+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"shufps		xmm2, xmm2, 0x00\n"
		"mulps		xmm7, xmm2\n"
		"addps		xmm6, xmm7\n"
		"addps		xmm6, %8\n"

		"cmpnltps	xmm6, %4\n"
		"movmskps	ecx, xmm6\n"

		"mulps		xmm0, %9\n"
		"mulps		xmm1, %10\n"
		"addps		xmm0, xmm1\n"
		"mulps		xmm2, %11\n"
		"addps		xmm0, xmm2\n"
		"addps		xmm0, %12\n"

		"cmpnltps	xmm0, %4\n"
		"movmskps	edx, xmm0\n"

		"and		edx, 3\n"
		"shl		edx, 4\n"
		"or			ecx, edx\n"

		"mov		byte ptr [edi], cl\n"

	"done%=:\n":
    :"m"(cullBits), "m"(planes), "m"(verts), "m"(numVerts),"m"(SIMD_SP_zero),
    "m"(p0), //%5
	"m"(p1),
	"m"(p2),
	"m"(p3),
	"m"(p4),
	"m"(p5),
	"m"(p6),
	"m"(p7)

	:"eax","ecx","edx","esi","edi"
	);


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

	asm (
		"mov		eax, %4\n"
		"mov		edx, %3\n"
		"mov		esi, %1\n"
		"mov		edi, %0\n"

		"mov		ecx, %2\n"
		"movss		xmm4, [ecx+ 0]\n"
		"movss		xmm5, [ecx+16]\n"
		"shufps		xmm4, xmm5, 0x00\n"
		"shufps		xmm4, xmm4, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/
		"movss		xmm5, [ecx+ 4]\n"
		"movss		xmm6, [ecx+20]\n"
		"shufps		xmm5, xmm6, 0x00\n"
		"shufps		xmm5, xmm5, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/
		"movss		xmm6, [ecx+ 8]\n"
		"movss		xmm7, [ecx+24]\n"
		"shufps		xmm6, xmm7, 0x00\n"
		"shufps		xmm6, xmm6, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/
		"movss		xmm7, [ecx+12]\n"
		"movss		xmm0, [ecx+28]\n"
		"shufps		xmm7, xmm0, 0x00\n"
		"shufps		xmm7, xmm7, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/

		"and		eax, ~1\n"
		"jz			done2%=\n"
		"add		edi, eax\n"
		"neg		eax\n"

	"loopVert2%=:\n"
		"movss		xmm0, [edx+0*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"movss		xmm1, [edx+1*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"shufps		xmm0, xmm1, 0x00\n"
		"mulps		xmm0, xmm4\n"
		"movss		xmm1, [edx+0*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"movss		xmm2, [edx+1*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"shufps		xmm1, xmm2, 0x00\n"
		"mulps		xmm1, xmm5\n"
		"movss		xmm2, [edx+0*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"movss		xmm3, [edx+1*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"shufps		xmm2, xmm3, 0x00\n"
		"mulps		xmm2, xmm6\n"
		"addps		xmm0, xmm1\n"
		"addps		xmm0, xmm2\n"
		"addps		xmm0, xmm7\n"
		"movaps		[esi], xmm0\n"
		"movaps		xmm1, xmm0\n"
		"movaps		xmm2, %5\n"
		"subps		xmm2, xmm0\n"
		"shufps		xmm0, xmm2, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
		"shufps		xmm1, xmm2, 0xEE\n"  /*R_SHUFFLEPS( 2, 3, 2, 3 )=11101110*/
		"add		edx, 2*"DRAWVERT_SIZE_STR"\n"
		"movmskps	ecx, xmm0\n"
		"mov		byte ptr [edi+eax+0], cl\n"
		"add		esi, 4*4\n"
		"movmskps	ecx, xmm1\n"
		"mov		byte ptr [edi+eax+1], cl\n"
		"add		eax, 2\n"
		"jl			loopVert2%=\n"

	"done2%=:\n"
		"mov		eax, %4\n"
		"and		eax, 1\n"
		"jz			done%=\n"

		"movss		xmm0, [edx+0*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"shufps		xmm0, xmm0, 0x00\n"
		"mulps		xmm0, xmm4\n"
		"movss		xmm1, [edx+0*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"shufps		xmm1, xmm1, 0x00\n"
		"mulps		xmm1, xmm5\n"
		"movss		xmm2, [edx+0*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"shufps		xmm2, xmm2, 0x00\n"
		"mulps		xmm2, xmm6\n"
		"addps		xmm0, xmm1\n"
		"addps		xmm0, xmm2\n"
		"addps		xmm0, xmm7\n"
		"movlps		[esi], xmm0\n"
		"movaps		xmm1, xmm0\n"
		"movaps		xmm2, %5\n"
		"subps		xmm2, xmm0\n"
		"shufps		xmm0, xmm2, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
		"movmskps	ecx, xmm0\n"
		"mov		byte ptr [edi], cl\n"

	"done%=:\n"
	::"m"(cullBits), "m"(texCoords), "m"(planes), "m"(verts), "m"(numVerts),"m"(SIMD_SP_one)
    :"eax","ebx","ecx","edx","esi","edi");

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

	asm (
		"mov			eax, %4\n"
		"shl			eax, 2\n"
		"mov			esi, %1\n"
		"mov			edi, %3\n"
		"mov			edx, %0\n"

		"add			edi, eax\n"
		"neg			eax\n"

		"add			eax, 4*12\n"
		"jge			done4%=\n"

	"loopPlane4%=:\n"
		"mov			ebx, [edi+eax-4*12+4]\n"
		"imul		ebx, "DRAWVERT_SIZE_STR"\n"
		"mov			ecx, [edi+eax-4*12+0]\n"
		"imul		ecx, "DRAWVERT_SIZE_STR"\n"

		"movss		xmm0, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"subss		xmm0, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"

		"movss		xmm1, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"subss		xmm1, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"

		"movss		xmm2, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"subss		xmm2, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"

		"mov			ebx, [edi+eax-4*12+8]\n"
		"imul		ebx, "DRAWVERT_SIZE_STR"\n"

		"shufps		xmm0, xmm0, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */
		"shufps		xmm1, xmm1, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */
		"shufps		xmm2, xmm2, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */

		"movss		xmm3, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"subss		xmm3, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"

		"movss		xmm4, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"subss		xmm4, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"

		"movss		xmm5, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"subss		xmm5, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"

		"mov			ebx, [edi+eax-3*12+4]\n"
		"imul		ebx, "DRAWVERT_SIZE_STR"\n"
		"mov			ecx, [edi+eax-3*12+0]\n"
		"imul		ecx, "DRAWVERT_SIZE_STR"\n"

		"shufps		xmm3, xmm3, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */
		"shufps		xmm4, xmm4, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */
		"shufps		xmm5, xmm5, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */

		"movss		xmm6, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"subss		xmm6, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"movss		xmm0, xmm6\n"

		"movss		xmm7, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"subss		xmm7, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"movss		xmm1, xmm7\n"

		"movss		xmm6, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"subss		xmm6, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"movss		xmm2, xmm6\n"

		"mov			ebx, [edi+eax-3*12+8]\n"
		"imul		ebx, "DRAWVERT_SIZE_STR"\n"

		"shufps		xmm0, xmm0, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */
		"shufps		xmm1, xmm1, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */
		"shufps		xmm2, xmm2, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */

		"movss		xmm7, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"subss		xmm7, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"movss		xmm3, xmm7\n"

		"movss		xmm6, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"subss		xmm6, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"movss		xmm4, xmm6\n"

		"movss		xmm7, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"subss		xmm7, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"movss		xmm5, xmm7\n"

		"mov			ebx, [edi+eax-2*12+4]\n"
		"imul		ebx, "DRAWVERT_SIZE_STR"\n"
		"mov			ecx, [edi+eax-2*12+0]\n"
		"imul		ecx, "DRAWVERT_SIZE_STR"\n"

		"shufps		xmm3, xmm3, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */
		"shufps		xmm4, xmm4, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */
		"shufps		xmm5, xmm5, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */

		"movss		xmm6, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"subss		xmm6, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"movss		xmm0, xmm6\n"

		"movss		xmm7, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"subss		xmm7, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"movss		xmm1, xmm7\n"

		"movss		xmm6, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"subss		xmm6, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"movss		xmm2, xmm6\n"

		"mov			ebx, [edi+eax-2*12+8]\n"
		"imul		ebx, "DRAWVERT_SIZE_STR"\n"

		"shufps		xmm0, xmm0, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */
		"shufps		xmm1, xmm1, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */
		"shufps		xmm2, xmm2, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */

		"movss		xmm7, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"subss		xmm7, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"movss		xmm3, xmm7\n"

		"movss		xmm6, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"subss		xmm6, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"movss		xmm4, xmm6\n"

		"movss		xmm7, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"subss		xmm7, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"movss		xmm5, xmm7\n"

		"mov			ebx, [edi+eax-1*12+4]\n"
		"imul		ebx, "DRAWVERT_SIZE_STR"\n"
		"mov			ecx, [edi+eax-1*12+0]\n"
		"imul		ecx, "DRAWVERT_SIZE_STR"\n"

		"shufps		xmm3, xmm3, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */
		"shufps		xmm4, xmm4, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */
		"shufps		xmm5, xmm5, 0x93\n"  /*R_SHUFFLEPS( 3, 0, 1, 2 )=10010011 */

		"movss		xmm6, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"subss		xmm6, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"movss		xmm0, xmm6\n"

		"movss		xmm7, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"subss		xmm7, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"movss		xmm1, xmm7\n"

		"movss		xmm6, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"subss		xmm6, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"movss		xmm2, xmm6\n"

		"mov			ebx, [edi+eax-1*12+8]\n"
		"imul		ebx, "DRAWVERT_SIZE_STR"\n"

		"movss		xmm7, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"subss		xmm7, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"movss		xmm3, xmm7\n"

		"movss		xmm6, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"subss		xmm6, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"movss		xmm4, xmm6\n"

		"movss		xmm7, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"subss		xmm7, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"movss		xmm5, xmm7\n"

		"movaps		xmm6, xmm4\n"
		"mulps		xmm6, xmm2\n"
		"movaps		xmm7, xmm5\n"
		"mulps		xmm7, xmm1\n"
		"subps		xmm6, xmm7\n"

		"mulps		xmm5, xmm0\n"
		"mulps		xmm2, xmm3\n"
		"subps		xmm5, xmm2\n"

		"mulps		xmm3, xmm1\n"
		"mulps		xmm4, xmm0\n"
		"subps		xmm3, xmm4\n"

		"movaps		xmm0, xmm6\n"
		"mulps		xmm6, xmm6\n"
		"movaps		xmm1, xmm5\n"
		"mulps		xmm5, xmm5\n"
		"movaps		xmm2, xmm3\n"
		"mulps		xmm3, xmm3\n"

		"addps		xmm3, xmm5\n"
		"addps		xmm3, xmm6\n"
		"rsqrtps		xmm3, xmm3\n"

		"add			edx, 4*16\n"
		"mov			ecx, [edi+eax-1*12+0]\n"
		"imul		ecx, "DRAWVERT_SIZE_STR"\n"

		"mulps		xmm0, xmm3\n"
		"mulps		xmm1, xmm3\n"
		"mulps		xmm2, xmm3\n"

		"movss		[edx-1*16+0], xmm0\n"
		"movss		[edx-1*16+4], xmm1\n"
		"movss		[edx-1*16+8], xmm2\n"

		"mulss		xmm0, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"mulss		xmm1, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"mulss		xmm2, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"

		"xorps		xmm0, %5\n"
		"subss		xmm0, xmm1\n"
		"subss		xmm0, xmm2\n"
		"movss		[edx-1*16+12], xmm0\n"

		"mov			ecx, [edi+eax-2*12+0]\n"
		"imul		ecx, "DRAWVERT_SIZE_STR"\n"

		"shufps		xmm0, xmm0, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm1, xmm1, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm2, xmm2, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/

		"movss		[edx-2*16+0], xmm0\n"
		"movss		[edx-2*16+4], xmm1\n"
		"movss		[edx-2*16+8], xmm2\n"

		"mulss		xmm0, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"mulss		xmm1, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"mulss		xmm2, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"

		"xorps		xmm0, %5\n"
		"subss		xmm0, xmm1\n"
		"subss		xmm0, xmm2\n"
		"movss		[edx-2*16+12], xmm0\n"

		"mov			ecx, [edi+eax-3*12+0]\n"
		"imul		ecx, "DRAWVERT_SIZE_STR"\n"

		"shufps		xmm0, xmm0, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm1, xmm1, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm2, xmm2, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/

		"movss		[edx-3*16+0], xmm0\n"
		"movss		[edx-3*16+4], xmm1\n"
		"movss		[edx-3*16+8], xmm2\n"

		"mulss		xmm0, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"mulss		xmm1, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"mulss		xmm2, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"

		"xorps		xmm0, %5\n"
		"subss		xmm0, xmm1\n"
		"subss		xmm0, xmm2\n"
		"movss		[edx-3*16+12], xmm0\n"

		"mov			ecx, [edi+eax-4*12+0]\n"
		"imul		ecx, "DRAWVERT_SIZE_STR"\n"

		"shufps		xmm0, xmm0, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm1, xmm1, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm2, xmm2, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/

		"movss		[edx-4*16+0], xmm0\n"
		"movss		[edx-4*16+4], xmm1\n"
		"movss		[edx-4*16+8], xmm2\n"

		"mulss		xmm0, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"mulss		xmm1, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"mulss		xmm2, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"

		"xorps		xmm0, %5\n"
		"subss		xmm0, xmm1\n"
		"subss		xmm0, xmm2\n"
		"movss		[edx-4*16+12], xmm0\n"

		"add			eax, 4*12\n"
		"jle			loopPlane4%=\n"

	"done4%=:\n"

		"sub			eax, 4*12\n"
		"jge			done%=\n"

	"loopPlane1%=:\n"
		"mov			ebx, [edi+eax+4]\n"
		"imul		ebx, "DRAWVERT_SIZE_STR"\n"
		"mov			ecx, [edi+eax+0]\n"
		"imul		ecx, "DRAWVERT_SIZE_STR"\n"

		"movss		xmm0, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"subss		xmm0, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"

		"movss		xmm1, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"subss		xmm1, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"

		"movss		xmm2, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"subss		xmm2, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"

		"mov			ebx, [edi+eax+8]\n"
		"imul		ebx, "DRAWVERT_SIZE_STR"\n"

		"movss		xmm3, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"subss		xmm3, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"

		"movss		xmm4, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"subss		xmm4, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"

		"movss		xmm5, [esi+ebx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"subss		xmm5, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"

		"movss		xmm6, xmm4\n"
		"mulss		xmm6, xmm2\n"
		"movss		xmm7, xmm5\n"
		"mulss		xmm7, xmm1\n"
		"subss		xmm6, xmm7\n"

		"mulss		xmm5, xmm0\n"
		"mulss		xmm2, xmm3\n"
		"subss		xmm5, xmm2\n"

		"mulss		xmm3, xmm1\n"
		"mulss		xmm4, xmm0\n"
		"subss		xmm3, xmm4\n"

		"movss		xmm0, xmm6\n"
		"mulss		xmm6, xmm6\n"
		"movss		xmm1, xmm5\n"
		"mulss		xmm5, xmm5\n"
		"movss		xmm2, xmm3\n"
		"mulss		xmm3, xmm3\n"

		"addss		xmm3, xmm5\n"
		"addss		xmm3, xmm6\n"
		"rsqrtss	xmm3, xmm3\n"

		"add		edx, 1*16\n"

		"mulss		xmm0, xmm3\n"
		"mulss		xmm1, xmm3\n"
		"mulss		xmm2, xmm3\n"

		"movss		[edx-1*16+0], xmm0\n"
		"movss		[edx-1*16+4], xmm1\n"
		"movss		[edx-1*16+8], xmm2\n"

		"mulss		xmm0, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"mulss		xmm1, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"mulss		xmm2, [esi+ecx+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"

		"xorps		xmm0, %5\n"
		"subss		xmm0, xmm1\n"
		"subss		xmm0, xmm2\n"
		"movss		[edx-1*16+12], xmm0\n"

		"add		eax, 1*12\n"
		"jl			loopPlane1%=\n"

	"done%=:\n"
	::"m"(planes), "m"(verts), "m"(numVerts), "m"(indexes), "m"(numIndexes),"m"(SIMD_SP_singleSignBitMask)
    :"eax","ebx","ecx","edx","esi","edi");

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
			planes->FitThroughPoint( a->xyz );
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
		planes->FitThroughPoint( a->xyz );
		planes++;
	}

#endif
}

/*
============
CSIMD_SSE::deriveTangents
============
*/
//#define REFINE_TANGENT_SQUAREROOT
#define FIX_DEGENERATE_TANGENT

void VPCALL CSIMD_SSE::deriveTangents( CPlane *planes, CVertex *verts, const int numVerts, const int *indexes, const int numIndexes ) {
	int i;

	SMF_ASSERT( sizeof( CVertex ) == DRAWVERT_SIZE );
	SMF_ASSERT( (int)&((CVertex *)0)->normal == DRAWVERT_NORMAL_OFFSET );
	SMF_ASSERT( (int)&((CVertex *)0)->tangents[0] == DRAWVERT_TANGENT0_OFFSET );
	SMF_ASSERT( (int)&((CVertex *)0)->tangents[1] == DRAWVERT_TANGENT1_OFFSET );

	SMF_ASSERT( planes != NULL );
	SMF_ASSERT( verts != NULL );
	SMF_ASSERT( numVerts >= 0 );

#ifdef REFINE_TANGENT_SQUAREROOT
	asm volatile(
		"movaps		xmm6, SIMD_SP_rsqrt_c0
		"movaps		xmm7, SIMD_SP_rsqrt_c1
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

		asm (
			// normal
			"movaps		xmm0, %6\n"
			"mulps		xmm0, %2\n"
			"movaps		xmm1, %7\n"
			"mulps		xmm1, %1\n"
			"subps		xmm0, xmm1\n"

			"movaps		xmm1, %7\n"
			"mulps		xmm1, %0\n"
			"movaps		xmm2, %5\n"
			"mulps		xmm2, %2\n"
			"subps		xmm1, xmm2\n"

			"movaps		xmm2, %5\n"
			"mulps		xmm2, %1\n"
			"movaps		xmm3, %6\n"
			"mulps		xmm3, %0\n"
			"subps		xmm2, xmm3\n"

			"movaps		xmm3, xmm0\n"
			"movaps		xmm4, xmm1\n"
			"movaps		xmm5, xmm2\n"

			"mulps		xmm3, xmm3\n"
			"mulps		xmm4, xmm4\n"
			"mulps		xmm5, xmm5\n"

			"addps		xmm3, xmm4\n"
			"addps		xmm3, xmm5\n"

#ifdef FIX_DEGENERATE_TANGENT
			"xorps		xmm4, xmm4\n"
			"cmpeqps	xmm4, xmm3\n"
			"andps		xmm4, %21\n"    /*SIMD_SP_tiny*/			// if values are zero replace them with a tiny number
			"andps		xmm3, %20\n"		// make sure the values are positive
			"orps		xmm3, xmm4\n"
#endif

#ifdef REFINE_TANGENT_SQUAREROOT
			"rsqrtps		xmm4, xmm3\n"
			"mulps		xmm3, xmm4\n"
			"mulps		xmm3, xmm4\n"
			"subps		xmm3, xmm6\n"
			"mulps		xmm4, xmm7\n"
			"mulps		xmm3, xmm4\n"
#else
			"rsqrtps		xmm3, xmm3\n"
#endif
			"mulps		xmm0, xmm3\n"
			"movaps		%10, xmm0\n"
			"mulps		xmm1, xmm3\n"
			"movaps		%11, xmm1\n"
			"mulps		xmm2, xmm3\n"
			"movaps		%12, xmm2\n"

			// area sign bit
			"movaps		xmm0, %3\n"
			"mulps		xmm0, %9\n"
			"movaps		xmm1, %4\n"
			"mulps		xmm1, %8\n"
			"subps		xmm0, xmm1\n"
			"andps		xmm0, %22\n"
			"movaps		%19, xmm0\n"

			// first tangent
			"movaps		xmm0, %0\n"
			"mulps		xmm0, %9\n"
			"movaps		xmm1, %4\n"
			"mulps		xmm1, %5\n"
			"subps		xmm0, xmm1\n"

			"movaps		xmm1, %1\n"
			"mulps		xmm1, %9\n"
			"movaps		xmm2, %4\n"
			"mulps		xmm2, %6\n"
			"subps		xmm1, xmm2\n"

			"movaps		xmm2, %2\n"
			"mulps		xmm2, %9\n"
			"movaps		xmm3, %4\n"
			"mulps		xmm3, %7\n"
			"subps		xmm2, xmm3\n"

			"movaps		xmm3, xmm0\n"
			"movaps		xmm4, xmm1\n"
			"movaps		xmm5, xmm2\n"

			"mulps		xmm3, xmm3\n"
			"mulps		xmm4, xmm4\n"
			"mulps		xmm5, xmm5\n"

			"addps		xmm3, xmm4\n"
			"addps		xmm3, xmm5\n"

#ifdef FIX_DEGENERATE_TANGENT
			"xorps		xmm4, xmm4\n"
			"cmpeqps	xmm4, xmm3\n"
			"andps		xmm4, %21\n"			// if values are zero replace them with a tiny number
			"andps		xmm3, %20\n"		// make sure the values are positive
			"orps		xmm3, xmm4\n"
#endif

#ifdef REFINE_TANGENT_SQUAREROOT
			"rsqrtps	xmm4, xmm3\n"
			"mulps		xmm3, xmm4\n"
			"mulps		xmm3, xmm4\n"
			"subps		xmm3, xmm6\n"
			"mulps		xmm4, xmm7\n"
			"mulps		xmm3, xmm4\n"
#else
			"rsqrtps		xmm3, xmm3\n"
#endif
			"xorps		xmm3, %19\n"

			"mulps		xmm0, xmm3\n"
			"movaps		%13, xmm0\n"
			"mulps		xmm1, xmm3\n"
			"movaps		%14, xmm1\n"
			"mulps		xmm2, xmm3\n"
			"movaps		%15, xmm2\n"

			// second tangent
			"movaps		xmm0, %3\n"
			"mulps		xmm0, %5\n"
			"movaps		xmm1, %0\n"
			"mulps		xmm1, %8\n"
			"subps		xmm0, xmm1\n"

			"movaps		xmm1, %3\n"
			"mulps		xmm1, %6\n"
			"movaps		xmm2, %1\n"
			"mulps		xmm2, %8\n"
			"subps		xmm1, xmm2\n"

			"movaps		xmm2, %3\n"
			"mulps		xmm2, %7\n"
			"movaps		xmm3, %2\n"
			"mulps		xmm3, %8\n"
			"subps		xmm2, xmm3\n"

			"movaps		xmm3, xmm0\n"
			"movaps		xmm4, xmm1\n"
			"movaps		xmm5, xmm2\n"

			"mulps		xmm3, xmm3\n"
			"mulps		xmm4, xmm4\n"
			"mulps		xmm5, xmm5\n"

			"addps		xmm3, xmm4\n"
			"addps		xmm3, xmm5\n"

#ifdef FIX_DEGENERATE_TANGENT
			"xorps		xmm4, xmm4\n"
			"cmpeqps		xmm4, xmm3\n"
			"andps		xmm4, %21\n"	/*SIMD_SP_tiny*/		// if values are zero replace them with a tiny number
			"andps		xmm3, %20\n"    /*SIMD_SP_absMask*/		// make sure the values are positive
			"orps		xmm3, xmm4\n"
#endif

#ifdef REFINE_TANGENT_SQUAREROOT
			"rsqrtps		xmm4, xmm3\n"
			"mulps		xmm3, xmm4\n"
			"mulps		xmm3, xmm4\n"
			"subps		xmm3, xmm6\n"
			"mulps		xmm4, xmm7\n"
			"mulps		xmm3, xmm4\n"
#else
			"rsqrtps		xmm3, xmm3\n"
#endif
			"xorps		xmm3, %19\n"

			"mulps		xmm0, xmm3\n"
			"movaps		%16, xmm0\n"
			"mulps		xmm1, xmm3\n"
			"movaps		%17, xmm1\n"
			"mulps		xmm2, xmm3\n"
			"movaps		%18, xmm2\n"
		::		"m"(d0),
		"m"(d1),
		"m"(d2),
		"m"(d3),
		"m"(d4),
		"m"(d5),
		"m"(d6),
		"m"(d7),
		"m"(d8),
		"m"(d9),
		"m"(n0),   //%10
		"m"(n1),
		"m"(n2),
		"m"(t0),   //%13
		"m"(t1),
		"m"(t2),
		"m"(t3),
		"m"(t4),
		"m"(t5),
        "m"(signBit),   //%19
        "m"(SIMD_SP_absMask),"m"(SIMD_SP_tiny),"m"(SIMD_SP_signBitMask)
        :"eax","ebx","ecx","edx","esi","edi");

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

		asm (
			// normal
			"movss		xmm0, %6\n"
			"mulss		xmm0, %2\n"
			"movss		xmm1, %7\n"
			"mulss		xmm1, %1\n"
			"subss		xmm0, xmm1\n"

			"movss		xmm1, %7\n"
			"mulss		xmm1, %0\n"
			"movss		xmm2, %5\n"
			"mulss		xmm2, %2\n"
			"subss		xmm1, xmm2\n"

			"movss		xmm2, %5\n"
			"mulss		xmm2, %1\n"
			"movss		xmm3, %6\n"
			"mulss		xmm3, %0\n"
			"subss		xmm2, xmm3\n"

			"movss		xmm3, xmm0\n"
			"movss		xmm4, xmm1\n"
			"movss		xmm5, xmm2\n"

			"mulss		xmm3, xmm3\n"
			"mulss		xmm4, xmm4\n"
			"mulss		xmm5, xmm5\n"

			"addss		xmm3, xmm4\n"
			"addss		xmm3, xmm5\n"

#ifdef FIX_DEGENERATE_TANGENT
			"xorps		xmm4, xmm4\n"
			"cmpeqps	xmm4, xmm3\n"
			"andps		xmm4, %21\n"			// if values are zero replace them with a tiny number
			"andps		xmm3, %20\n"		// make sure the values are positive
			"orps		xmm3, xmm4\n"
#endif

#ifdef REFINE_TANGENT_SQUAREROOT
			"rsqrtss	xmm4, xmm3\n"
			"mulss		xmm3, xmm4\n"
			"mulss		xmm3, xmm4\n"
			"subss		xmm3, xmm6\n"
			"mulss		xmm4, xmm7\n"
			"mulss		xmm3, xmm4\n"
#else
			"rsqrtss	xmm3, xmm3\n"
#endif
			"mulss		xmm0, xmm3\n"
			"movss		%10, xmm0\n"
			"mulss		xmm1, xmm3\n"
			"movss		%11, xmm1\n"
			"mulss		xmm2, xmm3\n"
			"movss		%12, xmm2\n"

			// area sign bit
			"movss		xmm0, %3\n"
			"mulss		xmm0, %9\n"
			"movss		xmm1, %4\n"
			"mulss		xmm1, %8\n"
			"subss		xmm0, xmm1\n"
			"andps		xmm0, %22\n"
			"movaps		%19, xmm0\n"

			// first tangent
			"movss		xmm0, %0\n"
			"mulss		xmm0, %9\n"
			"movss		xmm1, %4\n"
			"mulss		xmm1, %5\n"
			"subss		xmm0, xmm1\n"

			"movss		xmm1, %1\n"
			"mulss		xmm1, %9\n"
			"movss		xmm2, %4\n"
			"mulss		xmm2, %6\n"
			"subss		xmm1, xmm2\n"

			"movss		xmm2, %2\n"
			"mulss		xmm2, %9\n"
			"movss		xmm3, %4\n"
			"mulss		xmm3, %7\n"
			"subss		xmm2, xmm3\n"

			"movss		xmm3, xmm0\n"
			"movss		xmm4, xmm1\n"
			"movss		xmm5, xmm2\n"

			"mulss		xmm3, xmm3\n"
			"mulss		xmm4, xmm4\n"
			"mulss		xmm5, xmm5\n"

			"addss		xmm3, xmm4\n"
			"addss		xmm3, xmm5\n"

#ifdef FIX_DEGENERATE_TANGENT
			"xorps		xmm4, xmm4\n"
			"cmpeqps	xmm4, xmm3\n"
			"andps		xmm4, %21\n"			// if values are zero replace them with a tiny number
			"andps		xmm3, %20\n"		// make sure the values are positive
			"orps		xmm3, xmm4\n"
#endif

#ifdef REFINE_TANGENT_SQUAREROOT
			"rsqrtss	xmm4, xmm3\n"
			"mulss		xmm3, xmm4\n"
			"mulss		xmm3, xmm4\n"
			"subss		xmm3, xmm6\n"
			"mulss		xmm4, xmm7\n"
			"mulss		xmm3, xmm4\n"
#else
			"rsqrtss	xmm3, xmm3\n"
#endif
			"xorps		xmm3, %19\n"

			"mulss		xmm0, xmm3\n"
			"movss		%13, xmm0\n"
			"mulss		xmm1, xmm3\n"
			"movss		%14, xmm1\n"
			"mulss		xmm2, xmm3\n"
			"movss		%15, xmm2\n"

			// second tangent
			"movss		xmm0, %3\n"
			"mulss		xmm0, %5\n"
			"movss		xmm1, %0\n"
			"mulss		xmm1, %8\n"
			"subss		xmm0, xmm1\n"

			"movss		xmm1, %3\n"
			"mulss		xmm1, %6\n"
			"movss		xmm2, %1\n"
			"mulss		xmm2, %8\n"
			"subss		xmm1, xmm2\n"

			"movss		xmm2, %3\n"
			"mulss		xmm2, %7\n"
			"movss		xmm3, %2\n"
			"mulss		xmm3, %8\n"
			"subss		xmm2, xmm3\n"

			"movss		xmm3, xmm0\n"
			"movss		xmm4, xmm1\n"
			"movss		xmm5, xmm2\n"

			"mulss		xmm3, xmm3\n"
			"mulss		xmm4, xmm4\n"
			"mulss		xmm5, xmm5\n"

			"addss		xmm3, xmm4\n"
			"addss		xmm3, xmm5\n"

#ifdef FIX_DEGENERATE_TANGENT
			"xorps		xmm4, xmm4\n"
			"cmpeqps		xmm4, xmm3\n"
			"andps		xmm4, %21\n"			// if values are zero replace them with a tiny number
			"andps		xmm3, %20\n"		// make sure the values are positive
			"orps		xmm3, xmm4\n"
#endif

#ifdef REFINE_TANGENT_SQUAREROOT
			"rsqrtss	xmm4, xmm3\n"
			"mulss		xmm3, xmm4\n"
			"mulss		xmm3, xmm4\n"
			"subss		xmm3, xmm6\n"
			"mulss		xmm4, xmm7\n"
			"mulss		xmm3, xmm4\n"
#else
			"rsqrtss	xmm3, xmm3\n"
#endif
			"xorps		xmm3, %19\n"

			"mulss		xmm0, xmm3\n"
			"movss		%16, xmm0\n"
			"mulss		xmm1, xmm3\n"
			"movss		%17, xmm1\n"
			"mulss		xmm2, xmm3\n"
			"movss		%18, xmm2\n"
		::		"m"(d0),
		"m"(d1),
		"m"(d2),
		"m"(d3),
		"m"(d4),
		"m"(d5),
		"m"(d6),
		"m"(d7),
		"m"(d8),
		"m"(d9),
		"m"(n0),   //%10
		"m"(n1),
		"m"(n2),
		"m"(t0),   //%13
		"m"(t1),
		"m"(t2),
		"m"(t3),
		"m"(t4),
		"m"(t5),
        "m"(signBit),   //%19
        "m"(SIMD_SP_absMask),"m"(SIMD_SP_tiny),"m"(SIMD_SP_signBitMask)
        :"eax","ebx","ecx","edx","esi","edi");

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
CSIMD_SSE::deriveUnsmoothedTangents
============
*/
#define DERIVE_UNSMOOTHED_BITANGENT
/* s
void VPCALL CSIMD_SSE::deriveUnsmoothedTangents( CVertex *verts, const dominantTri_s *dominantTris, const int numVerts ) {
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

		asm (

			"movaps		xmm0, d6
			"mulps		xmm0, d2
			"movaps		xmm1, d7
			"mulps		xmm1, d1

			"movaps		xmm2, d7
			"mulps		xmm2, d0
			"movaps		xmm3, d5
			"mulps		xmm3, d2

			"movaps		xmm4, d5
			"mulps		xmm4, d1
			"movaps		xmm5, d6
			"mulps		xmm5, d0

			"subps		xmm0, xmm1\n"
			"subps		xmm2, xmm3\n"
			"movaps		xmm7, s2
			"subps		xmm4, xmm5\n"

			"mulps		xmm0, xmm7\n"
			"movaps		n0, xmm0
			"mulps		xmm2, xmm7\n"
			"movaps		n1, xmm2\n"
			"mulps		xmm4, xmm7\n"
			"movaps		n2, xmm4\n"

			"movaps		xmm0, d0
			"mulps		xmm0, d9
			"movaps		xmm1, d4
			"mulps		xmm1, d5

			"movaps		xmm2, d1
			"mulps		xmm2, d9
			"movaps		xmm3, d4
			"mulps		xmm3, d6

			"movaps		xmm4, d2
			"mulps		xmm4, d9
			"movaps		xmm5, d4
			"mulps		xmm5, d7

			"subps		xmm0, xmm1\n"
			"subps		xmm2, xmm3\n"
			"movaps		xmm7, s0
			"subps		xmm4, xmm5\n"

			"mulps		xmm0, xmm7\n"
			"movaps		t0, xmm0
			"mulps		xmm2, xmm7\n"
			"movaps		t1, xmm2\n"
			"mulps		xmm4, xmm7\n"
			"movaps		t2, xmm4\n"

#ifndef DERIVE_UNSMOOTHED_BITANGENT
			"movaps		xmm0, d3
			"mulps		xmm0, d5
			"movaps		xmm1, d0
			"mulps		xmm1, d8

			"movaps		xmm2, d3
			"mulps		xmm2, d6
			"movaps		xmm3, d1
			"mulps		xmm3, d8

			"movaps		xmm4, d3
			"mulps		xmm4, d7
			"movaps		xmm5, d2
			"mulps		xmm5, d8
#else
			"movaps		xmm0, n2
			"mulps		xmm0, t1
			"movaps		xmm1, n1
			"mulps		xmm1, t2

			"movaps		xmm2, n0
			"mulps		xmm2, t2
			"movaps		xmm3, n2
			"mulps		xmm3, t0

			"movaps		xmm4, n1
			"mulps		xmm4, t0
			"movaps		xmm5, n0
			"mulps		xmm5, t1
#endif
			"subps		xmm0, xmm1\n"
			"subps		xmm2, xmm3\n"
			"movaps		xmm7, s1
			"subps		xmm4, xmm5\n"

			"mulps		xmm0, xmm7\n"
			"movaps		t3, xmm0
			"mulps		xmm2, xmm7\n"
			"movaps		t4, xmm2\n"
			"mulps		xmm4, xmm7\n"
			"movaps		t5, xmm4\n"
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

		asm (

			"movss		xmm0, d6
			"mulss		xmm0, d2
			"movss		xmm1, d7
			"mulss		xmm1, d1

			"movss		xmm2, d7
			"mulss		xmm2, d0
			"movss		xmm3, d5
			"mulss		xmm3, d2

			"movss		xmm4, d5
			"mulss		xmm4, d1
			"movss		xmm5, d6
			"mulss		xmm5, d0

			"subss		xmm0, xmm1\n"
			"subss		xmm2, xmm3\n"
			"movss		xmm7, s2
			"subss		xmm4, xmm5\n"

			"mulss		xmm0, xmm7\n"
			"movss		n0, xmm0
			"mulss		xmm2, xmm7\n"
			"movss		n1, xmm2\n"
			"mulss		xmm4, xmm7\n"
			"movss		n2, xmm4\n"

			"movss		xmm0, d0
			"mulss		xmm0, d9
			"movss		xmm1, d4
			"mulss		xmm1, d5

			"movss		xmm2, d1
			"mulss		xmm2, d9
			"movss		xmm3, d4
			"mulss		xmm3, d6

			"movss		xmm4, d2
			"mulss		xmm4, d9
			"movss		xmm5, d4
			"mulss		xmm5, d7

			"subss		xmm0, xmm1\n"
			"subss		xmm2, xmm3\n"
			"movss		xmm7, s0
			"subss		xmm4, xmm5\n"

			"mulss		xmm0, xmm7\n"
			"movss		t0, xmm0
			"mulss		xmm2, xmm7\n"
			"movss		t1, xmm2\n"
			"mulss		xmm4, xmm7\n"
			"movss		t2, xmm4\n"

#ifndef DERIVE_UNSMOOTHED_BITANGENT
			"movss		xmm0, d3
			"mulss		xmm0, d5
			"movss		xmm1, d0
			"mulss		xmm1, d8

			"movss		xmm2, d3
			"mulss		xmm2, d6
			"movss		xmm3, d1
			"mulss		xmm3, d8

			"movss		xmm4, d3
			"mulss		xmm4, d7
			"movss		xmm5, d2
			"mulss		xmm5, d8
#else
			"movss		xmm0, n2
			"mulss		xmm0, t1
			"movss		xmm1, n1
			"mulss		xmm1, t2

			"movss		xmm2, n0
			"mulss		xmm2, t2
			"movss		xmm3, n2
			"mulss		xmm3, t0

			"movss		xmm4, n1
			"mulss		xmm4, t0
			"movss		xmm5, n0
			"mulss		xmm5, t1
#endif
			"subss		xmm0, xmm1\n"
			"subss		xmm2, xmm3\n"
			"movss		xmm7, s1
			"subss		xmm4, xmm5\n"

			"mulss		xmm0, xmm7\n"
			"movss		t3, xmm0
			"mulss		xmm2, xmm7\n"
			"movss		t4, xmm2\n"
			"mulss		xmm4, xmm7\n"
			"movss		t5, xmm4\n"
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
CSIMD_SSE::normalizeTangents
============
*/
void VPCALL CSIMD_SSE::normalizeTangents( CVertex *verts, const int numVerts ) {
	ALIGN16( float normal[12] );

	SMF_ASSERT( sizeof( CVertex ) == DRAWVERT_SIZE );
	SMF_ASSERT( (int)&((CVertex *)0)->normal == DRAWVERT_NORMAL_OFFSET );
	SMF_ASSERT( (int)&((CVertex *)0)->tangents[0] == DRAWVERT_TANGENT0_OFFSET );
	SMF_ASSERT( (int)&((CVertex *)0)->tangents[1] == DRAWVERT_TANGENT1_OFFSET );

	SMF_ASSERT( verts != NULL );
	SMF_ASSERT( numVerts >= 0 );

	asm (
		"mov		eax, %1\n"
		"test		eax, eax\n"
		"jz			done%=\n"
#ifdef REFINE_TANGENT_SQUAREROOT
		"movaps		xmm6, SIMD_SP_rsqrt_c0\n"
		"movaps		xmm7, SIMD_SP_rsqrt_c1\n"
#endif
		"mov		esi, %0\n"
		"imul		eax, "DRAWVERT_SIZE_STR"\n"
		"add			esi, eax\n"
		"neg			eax\n"
		"add			eax, "DRAWVERT_SIZE_STR"*4\n"
		"jle			loopVert4%=\n"

		"sub			eax, "DRAWVERT_SIZE_STR"*4\n"
		"jl			loopVert1%=\n"

	"loopVert4%=:\n"

		"sub			eax, "DRAWVERT_SIZE_STR"*4\n"

		// normalize 4 CVertex::normal

		"movss		xmm0, [esi+eax+"DRAWVERT_SIZE_STR"*0+"DRAWVERT_NORMAL_OFFSET_STR"+0]\n"	//  0,  X,  X,  X
		"movhps		xmm0, [esi+eax+"DRAWVERT_SIZE_STR"*1+"DRAWVERT_NORMAL_OFFSET_STR"+0]\n"	//  0,  X,  3,  4
		"movss		xmm2, [esi+eax+"DRAWVERT_SIZE_STR"*1+"DRAWVERT_NORMAL_OFFSET_STR"+8]\n"	//  5,  X,  X,  X
		"movhps		xmm2, [esi+eax+"DRAWVERT_SIZE_STR"*0+"DRAWVERT_NORMAL_OFFSET_STR"+4]\n"	//	5,  X,  1,  2
		"movss		xmm4, [esi+eax+"DRAWVERT_SIZE_STR"*2+"DRAWVERT_NORMAL_OFFSET_STR"+0]\n"	//  6,  X,  X,  X
		"movhps		xmm4, [esi+eax+"DRAWVERT_SIZE_STR"*3+"DRAWVERT_NORMAL_OFFSET_STR"+0]\n"	//  6,  X,  9, 10
		"movss		xmm3, [esi+eax+"DRAWVERT_SIZE_STR"*3+"DRAWVERT_NORMAL_OFFSET_STR"+8]\n"	// 11,  X,  X,  X
		"movhps		xmm3, [esi+eax+"DRAWVERT_SIZE_STR"*2+"DRAWVERT_NORMAL_OFFSET_STR"+4]\n"	// 11,  X,  7,  8

		"movaps		xmm1, xmm0\n"
		"movaps		xmm5, xmm2\n"
		"shufps		xmm0, xmm4, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/		//  0,  3,  6,  9
		"shufps		xmm2, xmm3, 0x33\n" /*R_SHUFFLEPS( 3, 0, 3, 0 )=00110011*/		//  2,  5,  8, 11
		"shufps		xmm1, xmm5, 0XAF\n"  /*R_SHUFFLEPS( 3, 3, 2, 2 )=10101111	*/	//  4,  4,  1,  1
		"shufps		xmm4, xmm3, 0XAF\n"  /*R_SHUFFLEPS( 3, 3, 2, 2 )=10101111	*/		// 10, 10,  7,  7
		"shufps		xmm1, xmm4, 0x22\n" /*R_SHUFFLEPS( 2, 0, 2, 0 )=00100010*/		//  1,  4,  7, 10

		"movaps		xmm3, xmm0\n"
		"movaps		xmm4, xmm1\n"
		"movaps		xmm5, xmm2\n"

		"mulps		xmm3, xmm3\n"
		"mulps		xmm4, xmm4\n"
		"mulps		xmm5, xmm5\n"
		"addps		xmm3, xmm4\n"
		"addps		xmm3, xmm5\n"

#ifdef REFINE_TANGENT_SQUAREROOT
		"rsqrtps	xmm4, xmm3\n"
		"mulps		xmm3, xmm4\n"
		"mulps		xmm3, xmm4\n"
		"subps		xmm3, xmm6\n"
		"mulps		xmm4, xmm7\n"
		"mulps		xmm3, xmm4\n"
#else
		"rsqrtps	xmm3, xmm3\n"
#endif

		"mulps		xmm0, xmm3\n"
		"mulps		xmm1, xmm3\n"
		"mulps		xmm2, xmm3\n"

		// save the 4 CVertex::normal to project the tangents

		"movaps		[%2+ 0], xmm0\n"
		"movaps		[%2+16], xmm1\n"
		"movaps		[%2+32], xmm2\n"

		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*0+"DRAWVERT_NORMAL_OFFSET_STR"+0], xmm0\n"
		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*0+"DRAWVERT_NORMAL_OFFSET_STR"+4], xmm1\n"
		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*0+"DRAWVERT_NORMAL_OFFSET_STR"+8], xmm2\n"

		"shufps		xmm0, xmm0, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm1, xmm1, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm2, xmm2, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/

		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*1+"DRAWVERT_NORMAL_OFFSET_STR"+0], xmm0\n"
		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*1+"DRAWVERT_NORMAL_OFFSET_STR"+4], xmm1\n"
		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*1+"DRAWVERT_NORMAL_OFFSET_STR"+8], xmm2\n"

		"shufps		xmm0, xmm0, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm1, xmm1, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm2, xmm2, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/

		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*2+"DRAWVERT_NORMAL_OFFSET_STR"+0], xmm0\n"
		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*2+"DRAWVERT_NORMAL_OFFSET_STR"+4], xmm1\n"
		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*2+"DRAWVERT_NORMAL_OFFSET_STR"+8], xmm2\n"

		"shufps		xmm0, xmm0, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm1, xmm1, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm2, xmm2, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/

		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*3+"DRAWVERT_NORMAL_OFFSET_STR"+0], xmm0\n"
		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*3+"DRAWVERT_NORMAL_OFFSET_STR"+4], xmm1\n"
		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*3+"DRAWVERT_NORMAL_OFFSET_STR"+8], xmm2\n"

		// project "and normalize 4 CVertex::tangent[0]

		"movss		xmm0, [esi+eax+"DRAWVERT_SIZE_STR"*0+"DRAWVERT_TANGENT0_OFFSET_STR"+0]\n"	//  0,  X,  X,  X
		"movhps		xmm0, [esi+eax+"DRAWVERT_SIZE_STR"*1+"DRAWVERT_TANGENT0_OFFSET_STR"+0]\n"	//  0,  X,  3,  4
		"movss		xmm2, [esi+eax+"DRAWVERT_SIZE_STR"*1+"DRAWVERT_TANGENT0_OFFSET_STR"+8]\n"	//  5,  X,  X,  X
		"movhps		xmm2, [esi+eax+"DRAWVERT_SIZE_STR"*0+"DRAWVERT_TANGENT0_OFFSET_STR"+4]\n"	//	5,  X,  1,  2
		"movss		xmm4, [esi+eax+"DRAWVERT_SIZE_STR"*2+"DRAWVERT_TANGENT0_OFFSET_STR"+0]\n"	//  6,  X,  X,  X
		"movhps		xmm4, [esi+eax+"DRAWVERT_SIZE_STR"*3+"DRAWVERT_TANGENT0_OFFSET_STR"+0]\n"	//  6,  X,  9, 10
		"movss		xmm3, [esi+eax+"DRAWVERT_SIZE_STR"*3+"DRAWVERT_TANGENT0_OFFSET_STR"+8]\n"	// 11,  X,  X,  X
		"movhps		xmm3, [esi+eax+"DRAWVERT_SIZE_STR"*2+"DRAWVERT_TANGENT0_OFFSET_STR"+4]\n"	// 11,  X,  7,  8

		"movaps		xmm1, xmm0\n"
		"movaps		xmm5, xmm2\n"
		"shufps		xmm0, xmm4, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/		//  0,  3,  6,  9
		"shufps		xmm2, xmm3, 0x33\n" /*R_SHUFFLEPS( 3, 0, 3, 0 )=00110011*/		//  2,  5,  8, 11
		"shufps		xmm1, xmm5, 0XAF \n"  /*R_SHUFFLEPS( 3, 3, 2, 2 )=10101111	*/		//  4,  4,  1,  1
		"shufps		xmm4, xmm3, 0XAF \n"  /*R_SHUFFLEPS( 3, 3, 2, 2 )=10101111	*/		// 10, 10,  7,  7
		"shufps		xmm1, xmm4, 0x22\n" /*R_SHUFFLEPS( 2, 0, 2, 0 )=00100010*/		//  1,  4,  7, 10

		"movaps		xmm3, xmm0\n"
		"movaps		xmm4, xmm1\n"
		"movaps		xmm5, xmm2\n"

		"mulps		xmm3, [%2+ 0]\n"
		"mulps		xmm4, [%2+16]\n"
		"mulps		xmm5, [%2+32]\n"
		"addps		xmm3, xmm4\n"
		"addps		xmm3, xmm5\n"

		"movaps		xmm4, xmm3\n"
		"movaps		xmm5, xmm3\n"
		"mulps		xmm3, [%2+ 0]\n"
		"mulps		xmm4, [%2+16]\n"
		"mulps		xmm5, [%2+32]\n"
		"subps		xmm0, xmm3\n"
		"subps		xmm1, xmm4\n"
		"subps		xmm2, xmm5\n"

		"movaps		xmm3, xmm0\n"
		"movaps		xmm4, xmm1\n"
		"movaps		xmm5, xmm2\n"

		"mulps		xmm3, xmm3\n"
		"mulps		xmm4, xmm4\n"
		"mulps		xmm5, xmm5\n"
		"addps		xmm3, xmm4\n"
		"addps		xmm3, xmm5\n"

#ifdef REFINE_TANGENT_SQUAREROOT
		"rsqrtps		xmm4, xmm3\n"
		"mulps		xmm3, xmm4\n"
		"mulps		xmm3, xmm4\n"
		"subps		xmm3, xmm6\n"
		"mulps		xmm4, xmm7\n"
		"mulps		xmm3, xmm4\n"
#else
		"rsqrtps		xmm3, xmm3\n"
#endif

		"mulps		xmm0, xmm3\n"
		"mulps		xmm1, xmm3\n"
		"mulps		xmm2, xmm3\n"

		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*0+"DRAWVERT_TANGENT0_OFFSET_STR"+0], xmm0\n"
		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*0+"DRAWVERT_TANGENT0_OFFSET_STR"+4], xmm1\n"
		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*0+"DRAWVERT_TANGENT0_OFFSET_STR"+8], xmm2\n"

		"shufps		xmm0, xmm0, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm1, xmm1, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm2, xmm2, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/

		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*1+"DRAWVERT_TANGENT0_OFFSET_STR"+0], xmm0\n"
		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*1+"DRAWVERT_TANGENT0_OFFSET_STR"+4], xmm1\n"
		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*1+"DRAWVERT_TANGENT0_OFFSET_STR"+8], xmm2\n"

		"shufps		xmm0, xmm0, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm1, xmm1, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm2, xmm2, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/

		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*2+"DRAWVERT_TANGENT0_OFFSET_STR"+0], xmm0\n"
		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*2+"DRAWVERT_TANGENT0_OFFSET_STR"+4], xmm1\n"
		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*2+"DRAWVERT_TANGENT0_OFFSET_STR"+8], xmm2\n"

		"shufps		xmm0, xmm0, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm1, xmm1, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm2, xmm2, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/

		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*3+"DRAWVERT_TANGENT0_OFFSET_STR"+0], xmm0\n"
		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*3+"DRAWVERT_TANGENT0_OFFSET_STR"+4], xmm1\n"
		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*3+"DRAWVERT_TANGENT0_OFFSET_STR"+8], xmm2\n"

		// project "and normalize 4 CVertex::tangent[1]

		"movss		xmm0, [esi+eax+"DRAWVERT_SIZE_STR"*0+"DRAWVERT_TANGENT1_OFFSET_STR"+0]\n"	//  0,  X,  X,  X
		"movhps		xmm0, [esi+eax+"DRAWVERT_SIZE_STR"*1+"DRAWVERT_TANGENT1_OFFSET_STR"+0]\n"	//  0,  X,  3,  4
		"movss		xmm2, [esi+eax+"DRAWVERT_SIZE_STR"*1+"DRAWVERT_TANGENT1_OFFSET_STR"+8]\n"	//  5,  X,  X,  X
		"movhps		xmm2, [esi+eax+"DRAWVERT_SIZE_STR"*0+"DRAWVERT_TANGENT1_OFFSET_STR"+4]\n"	//	5,  X,  1,  2
		"movss		xmm4, [esi+eax+"DRAWVERT_SIZE_STR"*2+"DRAWVERT_TANGENT1_OFFSET_STR"+0]\n"	//  6,  X,  X,  X
		"movhps		xmm4, [esi+eax+"DRAWVERT_SIZE_STR"*3+"DRAWVERT_TANGENT1_OFFSET_STR"+0]\n"	//  6,  X,  9, 10
		"movss		xmm3, [esi+eax+"DRAWVERT_SIZE_STR"*3+"DRAWVERT_TANGENT1_OFFSET_STR"+8]\n"	// 11,  X,  X,  X
		"movhps		xmm3, [esi+eax+"DRAWVERT_SIZE_STR"*2+"DRAWVERT_TANGENT1_OFFSET_STR"+4]\n"	// 11,  X,  7,  8

		"movaps		xmm1, xmm0\n"
		"movaps		xmm5, xmm2\n"
		"shufps		xmm0, xmm4, 0x88\n" /*R_SHUFFLEPS( 0, 2, 0, 2 )=10001000*/		//  0,  3,  6,  9
		"shufps		xmm2, xmm3, 0x33\n" /*R_SHUFFLEPS( 3, 0, 3, 0 )=00110011*/		//  2,  5,  8, 11
		"shufps		xmm1, xmm5, 0XAF \n"  /*R_SHUFFLEPS( 3, 3, 2, 2 )=10101111	*/		//  4,  4,  1,  1
		"shufps		xmm4, xmm3, 0XAF \n"  /*R_SHUFFLEPS( 3, 3, 2, 2 )=10101111	*/		// 10, 10,  7,  7
		"shufps		xmm1, xmm4, 0x22\n" /*R_SHUFFLEPS( 2, 0, 2, 0 )=00100010*/		//  1,  4,  7, 10

		"movaps		xmm3, xmm0\n"
		"movaps		xmm4, xmm1\n"
		"movaps		xmm5, xmm2\n"

		"mulps		xmm3, [%2+ 0]\n"
		"mulps		xmm4, [%2+16]\n"
		"mulps		xmm5, [%2+32]\n"
		"addps		xmm3, xmm4\n"
		"addps		xmm3, xmm5\n"

		"movaps		xmm4, xmm3\n"
		"movaps		xmm5, xmm3\n"
		"mulps		xmm3, [%2+ 0]\n"
		"mulps		xmm4, [%2+16]\n"
		"mulps		xmm5, [%2+32]\n"
		"subps		xmm0, xmm3\n"
		"subps		xmm1, xmm4\n"
		"subps		xmm2, xmm5\n"

		"movaps		xmm3, xmm0\n"
		"movaps		xmm4, xmm1\n"
		"movaps		xmm5, xmm2\n"

		"mulps		xmm3, xmm3\n"
		"mulps		xmm4, xmm4\n"
		"mulps		xmm5, xmm5\n"
		"addps		xmm3, xmm4\n"
		"addps		xmm3, xmm5\n"

#ifdef REFINE_TANGENT_SQUAREROOT
		"rsqrtps		xmm4, xmm3\n"
		"mulps		xmm3, xmm4\n"
		"mulps		xmm3, xmm4\n"
		"subps		xmm3, xmm6\n"
		"mulps		xmm4, xmm7\n"
		"mulps		xmm3, xmm4\n"
#else
		"rsqrtps		xmm3, xmm3\n"
#endif

		"mulps		xmm0, xmm3\n"
		"mulps		xmm1, xmm3\n"
		"mulps		xmm2, xmm3\n"

		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*0+"DRAWVERT_TANGENT1_OFFSET_STR"+0], xmm0\n"
		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*0+"DRAWVERT_TANGENT1_OFFSET_STR"+4], xmm1\n"
		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*0+"DRAWVERT_TANGENT1_OFFSET_STR"+8], xmm2\n"

		"shufps		xmm0, xmm0, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm1, xmm1, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm2, xmm2, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/

		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*1+"DRAWVERT_TANGENT1_OFFSET_STR"+0], xmm0\n"
		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*1+"DRAWVERT_TANGENT1_OFFSET_STR"+4], xmm1\n"
		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*1+"DRAWVERT_TANGENT1_OFFSET_STR"+8], xmm2\n"

		"shufps		xmm0, xmm0, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm1, xmm1, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm2, xmm2, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/

		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*2+"DRAWVERT_TANGENT1_OFFSET_STR"+0], xmm0\n"
		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*2+"DRAWVERT_TANGENT1_OFFSET_STR"+4], xmm1\n"
		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*2+"DRAWVERT_TANGENT1_OFFSET_STR"+8], xmm2\n"

		"shufps		xmm0, xmm0, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm1, xmm1, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shufps		xmm2, xmm2, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/

		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*3+"DRAWVERT_TANGENT1_OFFSET_STR"+0], xmm0\n"
		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*3+"DRAWVERT_TANGENT1_OFFSET_STR"+4], xmm1\n"
		"movss		[esi+eax+"DRAWVERT_SIZE_STR"*3+"DRAWVERT_TANGENT1_OFFSET_STR"+8], xmm2\n"

		"add			eax, "DRAWVERT_SIZE_STR"*8\n"

		"jle			loopVert4%=\n"

		"sub			eax, "DRAWVERT_SIZE_STR"*4\n"
		"jge			done%=\n"

	"loopVert1%=:\n"

		// normalize one CVertex::normal

		"movss		xmm0, [esi+eax+"DRAWVERT_NORMAL_OFFSET_STR"+0]\n"
		"movss		xmm1, [esi+eax+"DRAWVERT_NORMAL_OFFSET_STR"+4]\n"
		"movss		xmm2, [esi+eax+"DRAWVERT_NORMAL_OFFSET_STR"+8]\n"
		"movss		xmm3, xmm0\n"
		"movss		xmm4, xmm1\n"
		"movss		xmm5, xmm2\n"

		"mulss		xmm3, xmm3\n"
		"mulss		xmm4, xmm4\n"
		"mulss		xmm5, xmm5\n"
		"addss		xmm3, xmm4\n"
		"addss		xmm3, xmm5\n"

#ifdef REFINE_TANGENT_SQUAREROOT
		"rsqrtss		xmm4, xmm3\n"
		"mulss		xmm3, xmm4\n"
		"mulss		xmm3, xmm4\n"
		"subss		xmm3, xmm6\n"
		"mulss		xmm4, xmm7\n"
		"mulss		xmm3, xmm4\n"
#else
		"rsqrtss		xmm3, xmm3\n"
#endif

		"mulss		xmm0, xmm3\n"
		"mulss		xmm1, xmm3\n"
		"mulss		xmm2, xmm3\n"

		"movss		[esi+eax+"DRAWVERT_NORMAL_OFFSET_STR"+0], xmm0\n"
		"movss		[esi+eax+"DRAWVERT_NORMAL_OFFSET_STR"+4], xmm1\n"
		"movss		[esi+eax+"DRAWVERT_NORMAL_OFFSET_STR"+8], xmm2\n"

		// project "and normalize one CVertex::tangent[0]

		"movss		xmm0, [esi+eax+"DRAWVERT_TANGENT0_OFFSET_STR"+0]\n"
		"movss		xmm1, [esi+eax+"DRAWVERT_TANGENT0_OFFSET_STR"+4]\n"
		"movss		xmm2, [esi+eax+"DRAWVERT_TANGENT0_OFFSET_STR"+8]\n"
		"movss		xmm3, xmm0\n"
		"movss		xmm4, xmm1\n"
		"movss		xmm5, xmm2\n"

		"mulss		xmm3, [esi+eax+"DRAWVERT_NORMAL_OFFSET_STR"+0]\n"
		"mulss		xmm4, [esi+eax+"DRAWVERT_NORMAL_OFFSET_STR"+4]\n"
		"mulss		xmm5, [esi+eax+"DRAWVERT_NORMAL_OFFSET_STR"+8]\n"
		"addss		xmm3, xmm4\n"
		"addss		xmm3, xmm5\n"

		"movss		xmm4, xmm3\n"
		"movss		xmm5, xmm3\n"
		"mulss		xmm3, [esi+eax+"DRAWVERT_NORMAL_OFFSET_STR"+0]\n"
		"mulss		xmm4, [esi+eax+"DRAWVERT_NORMAL_OFFSET_STR"+4]\n"
		"mulss		xmm5, [esi+eax+"DRAWVERT_NORMAL_OFFSET_STR"+8]\n"
		"subss		xmm0, xmm3\n"
		"subss		xmm1, xmm4\n"
		"subss		xmm2, xmm5\n"

		"movss		xmm3, xmm0\n"
		"movss		xmm4, xmm1\n"
		"movss		xmm5, xmm2\n"

		"mulss		xmm3, xmm3\n"
		"mulss		xmm4, xmm4\n"
		"mulss		xmm5, xmm5\n"
		"addss		xmm3, xmm4\n"
		"addss		xmm3, xmm5\n"

#ifdef REFINE_TANGENT_SQUAREROOT
		"rsqrtss		xmm4, xmm3\n"
		"mulss		xmm3, xmm4\n"
		"mulss		xmm3, xmm4\n"
		"subss		xmm3, xmm6\n"
		"mulss		xmm4, xmm7\n"
		"mulss		xmm3, xmm4\n"
#else
		"rsqrtss		xmm3, xmm3\n"
#endif

		"mulss		xmm0, xmm3\n"
		"mulss		xmm1, xmm3\n"
		"mulss		xmm2, xmm3\n"

		"movss		[esi+eax+"DRAWVERT_TANGENT0_OFFSET_STR"+0], xmm0\n"
		"movss		[esi+eax+"DRAWVERT_TANGENT0_OFFSET_STR"+4], xmm1\n"
		"movss		[esi+eax+"DRAWVERT_TANGENT0_OFFSET_STR"+8], xmm2\n"

		// project "and normalize one CVertex::tangent[1]

		"movss		xmm0, [esi+eax+"DRAWVERT_TANGENT1_OFFSET_STR"+0]\n"
		"movss		xmm1, [esi+eax+"DRAWVERT_TANGENT1_OFFSET_STR"+4]\n"
		"movss		xmm2, [esi+eax+"DRAWVERT_TANGENT1_OFFSET_STR"+8]\n"
		"movss		xmm3, xmm0\n"
		"movss		xmm4, xmm1\n"
		"movss		xmm5, xmm2\n"

		"mulss		xmm3, [esi+eax+"DRAWVERT_NORMAL_OFFSET_STR"+0]\n"
		"mulss		xmm4, [esi+eax+"DRAWVERT_NORMAL_OFFSET_STR"+4]\n"
		"mulss		xmm5, [esi+eax+"DRAWVERT_NORMAL_OFFSET_STR"+8]\n"
		"addss		xmm3, xmm4\n"
		"addss		xmm3, xmm5\n"

		"movss		xmm4, xmm3\n"
		"movss		xmm5, xmm3\n"
		"mulss		xmm3, [esi+eax+"DRAWVERT_NORMAL_OFFSET_STR"+0]\n"
		"mulss		xmm4, [esi+eax+"DRAWVERT_NORMAL_OFFSET_STR"+4]\n"
		"mulss		xmm5, [esi+eax+"DRAWVERT_NORMAL_OFFSET_STR"+8]\n"
		"subss		xmm0, xmm3\n"
		"subss		xmm1, xmm4\n"
		"subss		xmm2, xmm5\n"

		"movss		xmm3, xmm0\n"
		"movss		xmm4, xmm1\n"
		"movss		xmm5, xmm2\n"

		"mulss		xmm3, xmm3\n"
		"mulss		xmm4, xmm4\n"
		"mulss		xmm5, xmm5\n"
		"addss		xmm3, xmm4\n"
		"addss		xmm3, xmm5\n"

#ifdef REFINE_TANGENT_SQUAREROOT
		"rsqrtss		xmm4, xmm3\n"
		"mulss		xmm3, xmm4\n"
		"mulss		xmm3, xmm4\n"
		"subss		xmm3, xmm6\n"
		"mulss		xmm4, xmm7\n"
		"mulss		xmm3, xmm4\n"
#else
		"rsqrtss		xmm3, xmm3\n"
#endif

		"mulss		xmm0, xmm3\n"
		"mulss		xmm1, xmm3\n"
		"mulss		xmm2, xmm3\n"

		"movss		[esi+eax+"DRAWVERT_TANGENT1_OFFSET_STR"+0], xmm0\n"
		"movss		[esi+eax+"DRAWVERT_TANGENT1_OFFSET_STR"+4], xmm1\n"
		"movss		[esi+eax+"DRAWVERT_TANGENT1_OFFSET_STR"+8], xmm2\n"

		"add			eax, "DRAWVERT_SIZE_STR"\n"

		"jl			loopVert1%=\n"
	"done%=:\n"
	:: "m"(verts), "m"(numVerts),"m"(normal)
    :"eax","ebx","ecx","edx","esi","edi");
}

/*
============
CSIMD_SSE::createTextureSpaceLightVectors
============
*/
void VPCALL CSIMD_SSE::createTextureSpaceLightVectors( CVec3D *lightVectors, const CVec3D &lightOrigin, const CVertex *verts, const int numVerts, const int *indexes, const int numIndexes ) {

	SMF_ASSERT( sizeof( CVertex ) == DRAWVERT_SIZE );
	SMF_ASSERT( (int)&((CVertex *)0)->xyz == DRAWVERT_XYZ_OFFSET );
	SMF_ASSERT( (int)&((CVertex *)0)->normal == DRAWVERT_NORMAL_OFFSET);
	SMF_ASSERT( (int)&((CVertex *)0)->tangents[0] == DRAWVERT_TANGENT0_OFFSET );
	SMF_ASSERT( (int)&((CVertex *)0)->tangents[1] == DRAWVERT_TANGENT1_OFFSET );

	bool *used = (bool *)_alloca16( numVerts * sizeof( used[0] ) );
	memset( used, 0, numVerts * sizeof( used[0] ) );

	for ( int i = numIndexes - 1; i >= 0; i-- ) {
		used[indexes[i]] = true;
	}

#if 0

	asm (

		"mov			eax, numVerts

		"mov			esi, used
		"add			esi, eax

		"mov			edi, verts
		"sub			edi, "DRAWVERT_SIZE_STR"

		"neg			eax
		"dec			eax

		"mov			ecx, lightOrigin
		"movss		xmm7, [ecx+0]
		"movhps		xmm7, [ecx+4]

		"mov			ecx, lightVectors
		"sub			ecx, 3*4

	loopVert:
		"inc			eax
		jge			done

		"add			edi, "DRAWVERT_SIZE_STR"
		"add			ecx, 3*4

		"cmp			sf_u8 ptr [esi+eax], 0
		"je			loopVert

		"movaps		xmm0, xmm7\n"
		"movss		xmm1, [edi+"DRAWVERT_XYZ_OFFSET_STR"+0]
		"movhps		xmm1, [edi+"DRAWVERT_XYZ_OFFSET_STR"+4]
		"subps		xmm0, xmm1\n"

		// 0,  X,  1,  2
		// 3,  X,  4,  5
		// 6,  X,  7,  8

		"movss		xmm2, [edi+"DRAWVERT_TANGENT0_OFFSET_STR"+0]
		"movhps		xmm2, [edi+"DRAWVERT_TANGENT0_OFFSET_STR"+4]
		"mulps		xmm2, xmm0

		"movss		xmm3, [edi+"DRAWVERT_TANGENT1_OFFSET_STR"+0]
		"movhps		xmm3, [edi+"DRAWVERT_TANGENT1_OFFSET_STR"+4]
		"mulps		xmm3, xmm0

		"movaps		xmm5, xmm2\n"								// xmm5 = 0,  X,  1,  2
		unpcklps	xmm5, xmm3\n"								// xmm5 = 0,  3,  X,  X
		unpckhps	xmm2, xmm3\n"								// xmm2 = 1,  4,  2,  5

		"movss		xmm4, [edi+"DRAWVERT_NORMAL_OFFSET_STR"+0]
		"movhps		xmm4, [edi+"DRAWVERT_NORMAL_OFFSET_STR"+4]
		"mulps		xmm4, xmm0

		"movlhps		xmm5, xmm4\n"								// xmm5 = 0,  3,  6,  X
		"movhlps		xmm4, xmm2\n"								// xmm4 = 2,  5,  7,  8
		"shufps		xmm2, xmm4, 0xB4 \n" /*R_SHUFFLEPS( 0, 1, 3, 2 )=10110100 */	// xmm2 = 2,  5,  8,  7

		"addps		xmm5, xmm4\n"
		"addps		xmm5, xmm2\n"
		"movlps		[ecx+0], xmm5
		"shufps		xmm5, xmm5, 0x4E\n"   /*R_SHUFFLEPS( 2, 3, 0, 1 )=01001110*/
		"movss		[ecx+8], xmm5

		"jmp			loopVert

	"done%=:\n"
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

	asm (

		xor			ecx, ecx\n"
		"mov			eax, numVerts

		"mov			esi, used
		"add			esi, eax

		"mov			edi, verts
		"sub			edi, "DRAWVERT_SIZE_STR"

		"neg			eax
		"dec			eax

	loopVert4:
		"inc			eax
		jge			done4

		"add			edi, "DRAWVERT_SIZE_STR"

		"cmp			byte ptr [esi+eax], 0
		"je			loopVert4

		mov			usedVertNums[ecx*4], eax

		"inc			ecx
		"cmp			ecx, 4

		"movss		xmm0, localLightOrigin[0]
		"movss		xmm1, localLightOrigin[4]
		"movss		xmm2, localLightOrigin[8]

		"subss		xmm0, [edi+"DRAWVERT_XYZ_OFFSET_STR"+0]
		"subss		xmm1, [edi+"DRAWVERT_XYZ_OFFSET_STR"+4]
		"subss		xmm2, [edi+"DRAWVERT_XYZ_OFFSET_STR"+8]

		"movss		lightDir0[ecx*4-4], xmm0\n"
		"movss		lightDir1[ecx*4-4], xmm1
		"movss		lightDir2[ecx*4-4], xmm2

		"movss		xmm3, [edi+"DRAWVERT_NORMAL_OFFSET_STR"+0]
		"movss		xmm4, [edi+"DRAWVERT_NORMAL_OFFSET_STR"+4]
		"movss		xmm5, [edi+"DRAWVERT_NORMAL_OFFSET_STR"+8]

		"movss		normal0[ecx*4-4], xmm3
		"movss		normal1[ecx*4-4], xmm4
		"movss		normal2[ecx*4-4], xmm5

		"movss		xmm0, [edi+"DRAWVERT_TANGENT0_OFFSET_STR"+0]
		"movss		xmm1, [edi+"DRAWVERT_TANGENT0_OFFSET_STR"+4]
		"movss		xmm2, [edi+"DRAWVERT_TANGENT0_OFFSET_STR"+8]

		"movss		tangent0[ecx*4-4], xmm0\n"
		"movss		tangent1[ecx*4-4], xmm1
		"movss		tangent2[ecx*4-4], xmm2

		"movss		xmm3, [edi+"DRAWVERT_TANGENT1_OFFSET_STR"+0]
		"movss		xmm4, [edi+"DRAWVERT_TANGENT1_OFFSET_STR"+4]
		"movss		xmm5, [edi+"DRAWVERT_TANGENT1_OFFSET_STR"+8]

		"movss		tangent3[ecx*4-4], xmm3
		"movss		tangent4[ecx*4-4], xmm4
		"movss		tangent5[ecx*4-4], xmm5

		"jl			loopVert4

		"movaps		xmm0, lightDir0
		"movaps		xmm1, lightDir1
		"movaps		xmm2, lightDir2

		"movaps		xmm3, tangent0
		"mulps		xmm3, xmm0
		"movaps		xmm4, tangent1
		"mulps		xmm4, xmm1\n"
		"movaps		xmm5, tangent2
		"mulps		xmm5, xmm2\n"

		"addps		xmm3, xmm4\n"
		"addps		xmm5, xmm3\n"

		"movaps		xmm3, tangent3
		"mulps		xmm3, xmm0
		"movaps		xmm4, tangent4
		"mulps		xmm4, xmm1\n"
		"movaps		xmm6, tangent5
		"mulps		xmm6, xmm2\n"

		"addps		xmm3, xmm4\n"
		"addps		xmm6, xmm3\n"

		"mulps		xmm0, normal0
		"mulps		xmm1, normal1
		"mulps		xmm2, normal2

		"addps		xmm0, xmm1\n"
		"addps		xmm0, xmm2\n"

		"mov			ecx, numVerts
		"imul		ecx, 12
		"mov			edx, usedVertNums[0]
		"add			ecx, lightVectors
		"imul		edx, 12

		"movss		[ecx+edx+0], xmm5
		"movss		[ecx+edx+4], xmm6
		"movss		[ecx+edx+8], xmm0\n"

		"shufps		xmm5, xmm5, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"mov			edx, usedVertNums[4]
		"shufps		xmm6, xmm6, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"imul		edx, 12
		"shufps		xmm0, xmm0, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/

		"movss		[ecx+edx+0], xmm5
		"movss		[ecx+edx+4], xmm6
		"movss		[ecx+edx+8], xmm0\n"

		"shufps		xmm5, xmm5, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"mov			edx, usedVertNums[8]
		"shufps		xmm6, xmm6, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"imul		edx, 12
		"shufps		xmm0, xmm0, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/

		"movss		[ecx+edx+0], xmm5
		"movss		[ecx+edx+4], xmm6
		"movss		[ecx+edx+8], xmm0\n"

		"shufps		xmm5, xmm5, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"mov			edx, usedVertNums[12]
		"shufps		xmm6, xmm6, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"imul		edx, 12
		"shufps		xmm0, xmm0, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/

		"movss		[ecx+edx+0], xmm5
		"movss		[ecx+edx+4], xmm6
		"movss		[ecx+edx+8], xmm0\n"

		xor			ecx, ecx\n"
		"jmp			loopVert4

	"done4:%=:\n"
		"test		ecx, ecx\n"
		"jz			done
		xor			eax, eax
		"mov			edi, numVerts
		"imul		edi, 12
		"add			edi, lightVectors

	"loopVert1%=:\n"
		"movss		xmm0, lightDir0[eax*4]
		"movss		xmm1, lightDir1[eax*4]
		"movss		xmm2, lightDir2[eax*4]

		"mov			edx, usedVertNums[eax*4]
		"imul		edx, 12

		"movss		xmm3, tangent0[eax*4]
		"mulss		xmm3, xmm0
		"movss		xmm4, tangent1[eax*4]
		"mulss		xmm4, xmm1\n"
		"movss		xmm5, tangent2[eax*4]
		"mulss		xmm5, xmm2\n"

		"addss		xmm3, xmm4\n"
		"addss		xmm5, xmm3\n"
		"movss		[edi+edx+0], xmm5

		"movss		xmm3, tangent3[eax*4]
		"mulss		xmm3, xmm0
		"movss		xmm4, tangent4[eax*4]
		"mulss		xmm4, xmm1\n"
		"movss		xmm6, tangent5[eax*4]
		"mulss		xmm6, xmm2\n"

		"addss		xmm3, xmm4\n"
		"addss		xmm6, xmm3\n"
		"movss		[edi+edx+4], xmm6

		"mulss		xmm0, normal0[eax*4]
		"mulss		xmm1, normal1[eax*4]
		"mulss		xmm2, normal2[eax*4]

		"addss		xmm0, xmm1\n"
		"addss		xmm0, xmm2\n"
		"movss		[edi+edx+8], xmm0\n"

		"inc			eax
		"dec			ecx
		jg			loopVert1

	"done%=:\n"
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
CSIMD_SSE::createSpecularTextureCoords
============
*/
void VPCALL CSIMD_SSE::createSpecularTextureCoords( CVec4D *texCoords, const CVec3D &lightOrigin, const CVec3D &viewOrigin, const CVertex *verts, const int numVerts, const int *indexes, const int numIndexes ) {

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

	asm (

		"mov			eax, numVerts

		"mov			esi, used
		"add			esi, eax

		"mov			edi, verts
		"sub			edi, "DRAWVERT_SIZE_STR"

		"neg			eax
		"dec			eax

		"mov			ecx, viewOrigin
		"movss		xmm6, [ecx+0]
		"movhps		xmm6, [ecx+4]

		"mov			ecx, lightOrigin
		"movss		xmm7, [ecx+0]
		"movhps		xmm7, [ecx+4]

		"mov			ecx, texCoords
		"sub			ecx, 4*4

	loopVert:
		"inc			eax
		jge			done

		"add			edi, "DRAWVERT_SIZE_STR"
		"add			ecx, 4*4

		"cmp			byte ptr [esi+eax], 0
		"je			loopVert

		"movaps		xmm0, xmm7\n"
		"movaps		xmm1, xmm6\n"
		"movss		xmm2, [edi+"DRAWVERT_XYZ_OFFSET_STR"+0]
		"movhps		xmm2, [edi+"DRAWVERT_XYZ_OFFSET_STR"+4]
		"subps		xmm0, xmm2\n"
		"subps		xmm1, xmm2\n"

		"movaps		xmm3, xmm0
		"movaps		xmm4, xmm1\n"
		"mulps		xmm3, xmm3\n"
		"mulps		xmm4, xmm4\n"

		// 0,  X,  1,  2
		// 3,  X,  4,  5

		"movaps		xmm5, xmm3\n"								// xmm5 = 0,  X,  1,  2
		unpcklps	xmm5, xmm4\n"			0x00					// xmm5 = 0,  3,  X,  X
		unpckhps	xmm3, xmm4\n"								// xmm3 = 1,  4,  2,  5
		"movhlps		xmm4, xmm3\n"								// xmm4 = 2,  5,  4,  5

		"addps		xmm5, xmm3\n"
		"addps		xmm5, xmm4\n"
		"shufps		xmm5, xmm5, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
		rsqrtps		xmm5, xmm5\n"

		"movaps		xmm4, xmm5\n"
		"shufps		xmm4, xmm4, 0x00\n"
		"shufps		xmm5, xmm5, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */

		"mulps		xmm0, xmm4\n"
		"mulps		xmm1, xmm5\n"
		"addps		xmm0, xmm1\n"

		"movss		xmm2, [edi+"DRAWVERT_TANGENT0_OFFSET_STR"+0]
		"movhps		xmm2, [edi+"DRAWVERT_TANGENT0_OFFSET_STR"+4]
		"mulps		xmm2, xmm0

		"movss		xmm3, [edi+"DRAWVERT_TANGENT1_OFFSET_STR"+0]
		"movhps		xmm3, [edi+"DRAWVERT_TANGENT1_OFFSET_STR"+4]
		"mulps		xmm3, xmm0

		"movss		xmm4, [edi+"DRAWVERT_NORMAL_OFFSET_STR"+0]
		"movhps		xmm4, [edi+"DRAWVERT_NORMAL_OFFSET_STR"+4]
		"mulps		xmm4, xmm0

		"movaps		xmm5, xmm2\n"								// xmm5 = 0,  X,  1,  2
		unpcklps	xmm5, xmm3\n"								// xmm5 = 0,  3,  X,  X
		unpckhps	xmm2, xmm3\n"								// xmm2 = 1,  4,  2,  5

		"movlhps		xmm5, xmm4\n"								// xmm5 = 0,  3,  6,  X
		"movhlps		xmm4, xmm2\n"								// xmm4 = 2,  5,  7,  8
		"shufps		xmm2, xmm4, 0xB4 \n" /*R_SHUFFLEPS( 0, 1, 3, 2 )=10110100 */	// xmm2 = 2,  5,  8,  7

		"movaps		xmm3, SIMD_SP_one

		"addps		xmm5, xmm4\n"
		"addps		xmm5, xmm2\n"
		"movaps		[ecx+0], xmm5
		"movss		[ecx+12], xmm3

		"jmp			loopVert

	"done%=:\n"
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

	asm (

		"xor		ecx, ecx\n"
		"mov		eax, %2\n"  /*numVerts*/

		"mov		esi, %23\n"
		"add		esi, eax\n"

		"mov		edi, %1\n"  /*verts*/
		"sub		edi, "DRAWVERT_SIZE_STR"\n"

		"neg		eax\n"
		"dec		eax\n"

	"loopVert4%=:\n"
		"inc		eax\n"
		"jge		done4%=\n"

		"add		edi, "DRAWVERT_SIZE_STR"\n"

		"cmp		byte ptr [esi+eax], 0\n"
		"je			loopVert4%=\n"

		"mov		%6[ecx*4], eax\n"  /*usedVertNums*/

		"inc		ecx\n"
		"cmp		ecx, 4\n"

		"movss		xmm3, %3[0]\n"  /*localLightOrigin*/
		"movss		xmm4, %3[4]\n"  /*localLightOrigin*/
		"movss		xmm5, %3[8]\n"  /*localLightOrigin*/

		"subss		xmm3, [edi+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"subss		xmm4, [edi+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"subss		xmm5, [edi+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"

		"movss		%7[ecx*4-4], xmm3\n" /*lightDir0*/
		"movss		%8[ecx*4-4], xmm4\n"  /*lightDir1*/
		"movss		%9[ecx*4-4], xmm5\n"  /*lightDir2*/

		"movss		xmm0, %4[0]\n"  /*localViewOrigin*/
		"movss		xmm1, %4[4]\n"  /*localViewOrigin*/
		"movss		xmm2, %4[8]\n"  /*localViewOrigin*/

		"subss		xmm0, [edi+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"subss		xmm1, [edi+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"subss		xmm2, [edi+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"

		"movss		%10[ecx*4-4], xmm0\n"  /*viewDir0*/
		"movss		%11[ecx*4-4], xmm1\n"  /*viewDir1*/
		"movss		%12[ecx*4-4], xmm2\n"  /*viewDir2*/

		"movss		xmm3, [edi+"DRAWVERT_NORMAL_OFFSET_STR"+0]\n"
		"movss		xmm4, [edi+"DRAWVERT_NORMAL_OFFSET_STR"+4]\n"
		"movss		xmm5, [edi+"DRAWVERT_NORMAL_OFFSET_STR"+8]\n"

		"movss		%13[ecx*4-4], xmm3\n"    /*normal0*/
		"movss		%14[ecx*4-4], xmm4\n"    /*normal1*/
		"movss		%15[ecx*4-4], xmm5\n"    /*normal2*/

		"movss		xmm0, [edi+"DRAWVERT_TANGENT0_OFFSET_STR"+0]\n"
		"movss		xmm1, [edi+"DRAWVERT_TANGENT0_OFFSET_STR"+4]\n"
		"movss		xmm2, [edi+"DRAWVERT_TANGENT0_OFFSET_STR"+8]\n"

		"movss		%16[ecx*4-4], xmm0\n"      /*tangent0*/
		"movss		%17[ecx*4-4], xmm1\n"      /*tangent1*/
		"movss		%18[ecx*4-4], xmm2\n"      /*tangent2*/

		"movss		xmm3, [edi+"DRAWVERT_TANGENT1_OFFSET_STR"+0]\n"
		"movss		xmm4, [edi+"DRAWVERT_TANGENT1_OFFSET_STR"+4]\n"
		"movss		xmm5, [edi+"DRAWVERT_TANGENT1_OFFSET_STR"+8]\n"

		"movss		%19[ecx*4-4], xmm3\n"      /*tangent3*/
		"movss		%20[ecx*4-4], xmm4\n"      /*tangent4*/
		"movss		%21[ecx*4-4], xmm5\n"      /*tangent5*/

		"jl			loopVert4%=\n"

		"movaps		xmm6, %7\n"  /*lightDir0*/
		"movaps		xmm0, xmm6\n"
		"mulps		xmm6, xmm6\n"
		"movaps		xmm7, %8\n"  /*lightDir1*/
		"movaps		xmm1, xmm7\n"
		"mulps		xmm7, xmm7\n"
		"addps		xmm6, xmm7\n"
		"movaps		xmm5, %9\n"  /*lightDir2*/
		"movaps		xmm2, xmm5\n"
		"mulps		xmm5, xmm5\n"
		"addps		xmm6, xmm5\n"
		"rsqrtps		xmm6, xmm6\n"

		"mulps		xmm0, xmm6\n"
		"mulps		xmm1, xmm6\n"
		"mulps		xmm2, xmm6\n"

		"movaps		xmm3, %10\n"  /*viewDir0*/
		"movaps		xmm7, xmm3\n"
		"mulps		xmm7, xmm7\n"
		"movaps		xmm4, %11\n"  /*viewDir0*/
		"movaps		xmm6, xmm4\n"
		"mulps		xmm6, xmm6\n"
		"addps		xmm7, xmm6\n"
		"movaps		xmm5, %12\n"  /*viewDir0*/
		"movaps		xmm6, xmm5\n"
		"mulps		xmm6, xmm6\n"
		"addps		xmm7, xmm6\n"
		"rsqrtps	xmm7, xmm7\n"

		"mulps		xmm3, xmm7\n"
		"addps		xmm0, xmm3\n"
		"mulps		xmm4, xmm7\n"
		"addps		xmm1, xmm4\n"
		"mulps		xmm5, xmm7\n"
		"addps		xmm2, xmm5\n"

		"movaps		xmm3, %16\n"      /*tangent0*/
		"mulps		xmm3, xmm0\n"
		"movaps		xmm4, %17\n"      /*tangent1*/
		"mulps		xmm4, xmm1\n"
		"addps		xmm3, xmm4\n"
		"movaps		xmm5, %18\n"      /*tangent2*/
		"mulps		xmm5, xmm2\n"
		"addps		xmm5, xmm3\n"

		"movaps		xmm3, %19\n"      /*tangent3*/
		"mulps		xmm3, xmm0\n"
		"movaps		xmm4, %20\n"      /*tangent4*/
		"mulps		xmm4, xmm1\n"
		"addps		xmm3, xmm4\n"
		"movaps		xmm6, %21\n"      /*tangent5*/
		"mulps		xmm6, xmm2\n"
		"addps		xmm6, xmm3\n"

		"mulps		xmm0, %13\n"    /*normal0*/
		"mulps		xmm1, %14\n"    /*normal1*/
		"addps		xmm0, xmm1\n"
		"mulps		xmm2, %15\n"    /*normal2*/
		"addps		xmm0, xmm2\n"

		"mov		ecx, %2\n"   /*numVerts*/
		"shl		ecx, 4\n"
		"mov		edx, %6[0]\n"  /*usedVertNums*/
		"add		ecx, %0\n" /*texCoords*/
		"shl		edx, 4\n"
		"movss		xmm3, %22\n"

		"movss		[ecx+edx+0], xmm5\n"
		"movss		[ecx+edx+4], xmm6\n"
		"movss		[ecx+edx+8], xmm0\n"
		"movss		[ecx+edx+12], xmm3\n"

		"shufps		xmm5, xmm5, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"mov		edx, %6[4]\n"  /*usedVertNums*/
		"shufps		xmm6, xmm6, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shl		edx, 4\n"
		"shufps		xmm0, xmm0, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/

		"movss		[ecx+edx+0], xmm5\n"
		"movss		[ecx+edx+4], xmm6\n"
		"movss		[ecx+edx+8], xmm0\n"
		"movss		[ecx+edx+12], xmm3\n"

		"shufps		xmm5, xmm5, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"mov		edx, %6[8]\n"  /*usedVertNums*/
		"shufps		xmm6, xmm6, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shl		edx, 4\n"
		"shufps		xmm0, xmm0, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/

		"movss		[ecx+edx+0], xmm5\n"
		"movss		[ecx+edx+4], xmm6\n"
		"movss		[ecx+edx+8], xmm0\n"
		"movss		[ecx+edx+12], xmm3\n"

		"shufps		xmm5, xmm5, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"mov		edx, %6[12]\n"  /*usedVertNums*/
		"shufps		xmm6, xmm6, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/
		"shl		edx, 4\n"
		"shufps		xmm0, xmm0, 0x39\n" /* R_SHUFFLEPS( 1, 2, 3, 0 )=00111001*/

		"movss		[ecx+edx+0], xmm5\n"
		"movss		[ecx+edx+4], xmm6\n"
		"movss		[ecx+edx+8], xmm0\n"
		"movss		[ecx+edx+12], xmm3\n"

		"xor		ecx, ecx\n"
		"jmp		loopVert4%=\n"

	"done4%=:\n"
		"test		ecx, ecx\n"
		"jz			done%=\n"
		"xor		eax, eax\n"
		"mov		edi, %2\n"   /*numVerts*/
		"shl		edi, 4\n"
		"add		edi, %0\n"  /*texCoords*/

	"loopVert1%=:\n"
		"movss		xmm6, %7[eax*4]\n"  /*lightDir0*/
		"movss		xmm0, xmm6\n"
		"mulss		xmm6, xmm6\n"
		"movss		xmm7, %8[eax*4]\n"  /*lightDir1*/
		"movss		xmm1, xmm7\n"
		"mulss		xmm7, xmm7\n"
		"addss		xmm6, xmm7\n"
		"movss		xmm5, %9[eax*4]\n"  /*lightDir2*/
		"movss		xmm2, xmm5\n"
		"mulss		xmm5, xmm5\n"
		"addss		xmm6, xmm5\n"
		"rsqrtss	xmm6, xmm6\n"

		"mulss		xmm0, xmm6\n"
		"mulss		xmm1, xmm6\n"
		"mulss		xmm2, xmm6\n"

		"movss		xmm3, %10[eax*4]\n"  /*viewDir0*/
		"movss		xmm7, xmm3\n"
		"mulss		xmm7, xmm7\n"
		"movss		xmm4, %11[eax*4]\n"  /*viewDir1*/
		"movss		xmm6, xmm4\n"
		"mulss		xmm6, xmm6\n"
		"addss		xmm7, xmm6\n"
		"movss		xmm5, %12[eax*4]\n"  /*viewDir2*/
		"movss		xmm6, xmm5\n"
		"mulss		xmm6, xmm6\n"
		"addss		xmm7, xmm6\n"
		"rsqrtss	xmm7, xmm7\n"

		"mulss		xmm3, xmm7\n"
		"addss		xmm0, xmm3\n"
		"mulss		xmm4, xmm7\n"
		"addss		xmm1, xmm4\n"
		"mulss		xmm5, xmm7\n"
		"addss		xmm2, xmm5\n"

		"mov		edx, %6[eax*4]\n"  /**/
		"shl		edx, 4\n"

		"movss		xmm3, %16[eax*4]\n"      /*tangent0*/
		"mulss		xmm3, xmm0\n"
		"movss		xmm4, %17[eax*4]\n"      /*tangent1*/
		"mulss		xmm4, xmm1\n"
		"addss		xmm3, xmm4\n"
		"movss		xmm5, %18[eax*4]\n"      /*tangent2*/
		"mulss		xmm5, xmm2\n"
		"addss		xmm5, xmm3\n"
		"movss		[edi+edx+0], xmm5\n"

		"movss		xmm3, %19[eax*4]\n"      /*tangent3*/
		"mulss		xmm3, xmm0\n"
		"movss		xmm4, %20[eax*4]\n"      /*tangent4*/
		"mulss		xmm4, xmm1\n"
		"addss		xmm3, xmm4\n"
		"movss		xmm6, %21[eax*4]\n"      /*tangent5*/
		"mulss		xmm6, xmm2\n"
		"addss		xmm6, xmm3\n"
		"movss		[edi+edx+4], xmm6\n"

		"mulss		xmm0, %13[eax*4]\n"     /*normal0*/
		"mulss		xmm1, %14[eax*4]\n"    /*normal1*/
		"addss		xmm0, xmm1\n"
		"mulss		xmm2, %15[eax*4]\n"    /*normal2*/
		"addss		xmm0, xmm2\n"
		"movss		[edi+edx+8], xmm0\n"

		"movss		xmm3, %5\n"
		"movss		[edi+edx+12], xmm3\n"

		"inc		eax\n"
		"dec		ecx\n"
		"jg			loopVert1%=\n"

	"done%=:\n"
	::"m"(texCoords), "m"(verts), "m"(numVerts),
    "m"(localLightOrigin),"m"(localViewOrigin),"m"(SIMD_SP_one),
    "m"(usedVertNums),
	"m"(lightDir0),
	"m"(lightDir1),
	"m"(lightDir2),
	"m"(viewDir0), //%10
	"m"(viewDir1),
	"m"(viewDir2),
	"m"(normal0),  //%13
	"m"(normal1),
	"m"(normal2),
	"m"(tangent0),
	"m"(tangent1),
	"m"(tangent2),
	"m"(tangent3),
	"m"(tangent4),
	"m"(tangent5),
    "m"(SIMD_SP_one), //%22
    "m"(used)
    :"eax","ebx","ecx","edx","esi","edi");

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
CSIMD_SSE::createShadowCache
============
*/
int VPCALL CSIMD_SSE::createShadowCache( CVec4D *vertexCache, int *vertRemap, const CVec3D &lightOrigin, const CVertex *verts, const int numVerts ) {
#if 1
	int outVerts;
    float *lighOriginPtr = lightOrigin.toFloatPtr();
	asm (
		"push		ebx\n"

		"mov		esi, %3\n"    /*lightOrigin*/
		"movaps		xmm5, %0\n"  /*SIMD_SP_lastOne*/
		"movss		xmm6, [esi+0]\n"
		"movhps		xmm6, [esi+4]\n"
		"shufps		xmm6, xmm6, 0x78\n"   /*R_SHUFFLEPS( 0, 2, 3, 1 )=01111000*/
		"orps		xmm6, %0\n"
		"movaps		xmm7, xmm6\n"

		"xor		ebx, ebx\n"
		"xor		ecx, ecx\n"

		"mov		edx, %2\n"     /*vertRemap*/
		"mov		esi, %4\n"        /*verts*/
		"mov		edi, %1\n"      /*vertexCache*/
		"mov		eax, %5\n"        /*numVerts*/
		"and		eax, ~3\n"
		"jz			done4%=\n"
		"shl		eax, 2\n"
		"add		edx, eax\n"
		"neg		eax\n"

	"loop4%=:\n"
		"prefetchnta	[edx+128]\n"
		"prefetchnta	[esi+4*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"]\n"

		"cmp        dword ptr [edx+eax+0], ebx\n"
		"jne        skip1%=\n"

		"mov		dword ptr [edx+eax+0], ecx\n"
		"movss		xmm0, [esi+0*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"movhps		xmm0, [esi+0*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"add		ecx, 2\n"
		"shufps		xmm0, xmm0, 0x4E\n" /*R_SHUFFLEPS( 2, 3, 0, 1 )=01001110*/
		"orps		xmm0, xmm5\n"
		"movaps		[edi+0*16], xmm0\n"
		"subps		xmm0, xmm6\n"
		"movaps		[edi+1*16], xmm0\n"
		"add		edi, 2*16\n"

	"skip1%=:\n"
		"cmp        dword ptr [edx+eax+4], ebx\n"
		"jne        skip2%=\n"

		"mov		dword ptr [edx+eax+4], ecx\n"
		"movss		xmm1, [esi+1*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"movhps		xmm1, [esi+1*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"add		ecx, 2\n"
		"shufps		xmm1, xmm1, 0x78\n"   /*R_SHUFFLEPS( 0, 2, 3, 1 )=01111000*/
		"orps		xmm1, xmm5\n"
		"movaps		[edi+0*16], xmm1\n"
		"subps		xmm1, xmm7\n"
		"movaps		[edi+1*16], xmm1\n"
		"add		edi, 2*16\n"

	"skip2%=:\n"
		"cmp        dword ptr [edx+eax+8], ebx\n"
		"jne        skip3%=\n"

		"mov		dword ptr [edx+eax+8], ecx\n"
		"movss		xmm2, [esi+2*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"movhps		xmm2, [esi+2*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"add		ecx, 2\n"
		"shufps		xmm2, xmm2, 0x4E\n" /*R_SHUFFLEPS( 2, 3, 0, 1 )=01001110*/
		"orps		xmm2, xmm5\n"
		"movaps		[edi+0*16], xmm2\n"
		"subps		xmm2, xmm6\n"
		"movaps		[edi+1*16], xmm2\n"
		"add		edi, 2*16\n"

	"skip3%=:\n"
		"cmp        dword ptr [edx+eax+12], ebx\n"
		"jne        skip4%=\n"

		"mov		dword ptr [edx+eax+12], ecx\n"
		"movss		xmm3, [esi+3*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"movhps		xmm3, [esi+3*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"add		ecx, 2\n"
		"shufps		xmm3, xmm3, 0x78\n"   /*R_SHUFFLEPS( 0, 2, 3, 1 )=01111000*/
		"orps		xmm3, xmm5\n"
		"movaps		[edi+0*16], xmm3\n"
		"subps		xmm3, xmm7\n"
		"movaps		[edi+1*16], xmm3\n"
		"add		edi, 2*16\n"

	"skip4%=:\n"
		"add		esi, 4*"DRAWVERT_SIZE_STR"\n"
		"add		eax, 4*4\n"
		"jl			loop4%=\n"

	"done4%=:\n"
		"mov		eax, %5\n"        /*numVerts*/
		"and		eax, 3\n"
		"jz			done1%=\n"
		"shl		eax, 2\n"
		"add		edx, eax\n"
		"neg		eax\n"

	"loop1%=:\n"
		"cmp         dword ptr [edx+eax+0], ebx\n"
		"jne         skip0%=\n"

		"mov			dword ptr [edx+eax+0], ecx\n"
		"movss		xmm0, [esi+0*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"movhps		xmm0, [esi+0*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"add		ecx, 2\n"
		"shufps		xmm0, xmm0, 0x4E\n"   /*R_SHUFFLEPS( 2, 3, 0, 1 )=01001110*/
		"orps		xmm0, xmm5\n"
		"movaps		[edi+0*16], xmm0\n"
		"subps		xmm0, xmm6\n"
		"movaps		[edi+1*16], xmm0\n"
		"add		edi, 2*16\n"

	"skip0%=:\n"

		"add		esi, "DRAWVERT_SIZE_STR"\n"
		"add		eax, 4\n"
		"jl			loop1%=\n"

	"done1%=:\n"
		"pop		ebx\n"
		"mov		%6, ecx\n"
	:
    :"m"(SIMD_SP_lastOne), "m"(vertexCache), "m"(vertRemap), "m"(lighOriginPtr), "m"(verts), "m"(numVerts) ,"m"(outVerts)
    :);
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
		// wrap around point "and causing depth fighting with the rear caps
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
CSIMD_SSE::createVertexProgramShadowCache
============
*/
int VPCALL CSIMD_SSE::createVertexProgramShadowCache( CVec4D *vertexCache, const CVertex *verts, const int numVerts ) {
#if 1

	asm (
		"movaps		xmm4, %0\n"
		"movaps		xmm5, xmm4\n"
		"movaps		xmm6, xmm4\n"
		"movaps		xmm7, xmm4\n"

		"mov		esi, %2\n"
		"mov		edi, %1\n"
		"mov		eax, %3\n"
		"and		eax, ~3\n"
		"jz			done4%=\n"
		"shl		eax, 5\n"
		"add		edi, eax\n"
		"neg		eax\n"

	"loop4%=:\n"
		"prefetchnta	[esi+4*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"]\n"

		"movss		xmm0, [esi+0*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"movhps		xmm0, [esi+0*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"shufps		xmm0, xmm0, 0x4E\n" /*R_SHUFFLEPS( 2, 3, 0, 1 )=01001110*/
		"movaps		[edi+eax+1*16], xmm0\n"
		"orps		xmm0, xmm4\n"
		"movaps		[edi+eax+0*16], xmm0\n"

		"movss		xmm1, [esi+1*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"movhps		xmm1, [esi+1*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"shufps		xmm1, xmm1, 0x78\n"   /*R_SHUFFLEPS( 0, 2, 3, 1 )=01111000*/
		"movaps		[edi+eax+3*16], xmm1\n"
		"orps		xmm1, xmm5\n"
		"movaps		[edi+eax+2*16], xmm1\n"

		"movss		xmm2, [esi+2*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"movhps		xmm2, [esi+2*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"shufps		xmm2, xmm2, 0x4E\n" /*R_SHUFFLEPS( 2, 3, 0, 1 )=01001110*/
		"movaps		[edi+eax+5*16], xmm2\n"
		"orps		xmm2, xmm6\n"
		"movaps		[edi+eax+4*16], xmm2\n"

		"movss		xmm3, [esi+3*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"movhps		xmm3, [esi+3*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+4]\n"
		"shufps		xmm3, xmm3, 0x78\n"   /*R_SHUFFLEPS( 0, 2, 3, 1 )=01111000*/
		"movaps		[edi+eax+7*16], xmm3\n"
		"orps		xmm3, xmm7\n"
		"movaps		[edi+eax+6*16], xmm3\n"

		"add		esi, 4*"DRAWVERT_SIZE_STR"\n"
		"add		eax, 4*8*4\n"
		"jl			loop4%=\n"

	"done4%=:\n"
		"mov		eax, %3\n"
		"and		eax, 3\n"
		"jz			done1%=\n"
		"shl		eax, 5\n"
		"add		edi, eax\n"
		"neg		eax\n"

	"loop1%=:\n"
		"movss		xmm0, [esi+0*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+8]\n"
		"movhps		xmm0, [esi+0*"DRAWVERT_SIZE_STR"+"DRAWVERT_XYZ_OFFSET_STR"+0]\n"
		"shufps		xmm0, xmm0, 0x4E\n" /*R_SHUFFLEPS( 2, 3, 0, 1 )=01001110*/
		"movaps		[edi+eax+1*16], xmm0\n"
		"orps		xmm0, xmm4\n"
		"movaps		[edi+eax+0*16], xmm0\n"

		"add		esi, "DRAWVERT_SIZE_STR"\n"
		"add		eax, 8*4\n"
		"jl			loop1%=\n"

	"done1%=:\n"
	::"m"(SIMD_SP_lastOne), "r"(vertexCache), "r"(verts), "m"(numVerts)
    :);
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
	asm (
		"mov			esi, %1\n"
		"mov			edi, %0\n"

		"mov			eax, %2\n"
		"and			eax, ~1\n"
		"jz			done2%=\n"
		"shl			eax, 1\n"
		"add			esi, eax\n"
		"neg			eax\n"

		/*align		16  */
	"loop2%=:\n"
		"add			edi, 2*4*4\n"

		"movsx		ecx, word ptr [esi+eax+0]\n"
		"cvtsi2ss	xmm0, ecx\n"
		"shufps		xmm0, xmm0, 0x00\n"
		"movlps		[edi-2*4*4+0], xmm0\n"
		"movhps		[edi-2*4*4+8], xmm0\n"

		"movsx		edx, word ptr [esi+eax+2]\n"
		"cvtsi2ss	xmm1, edx\n"
		"shufps		xmm1, xmm1, 0x00\n"
		"movlps		[edi-1*4*4+0], xmm1\n"
		"movhps		[edi-1*4*4+8], xmm1\n"

		"add			eax, 2*2\n"
		"jl			loop2%=\n"

	"done2%=:\n"
		"mov			eax, %2\n"
		"and			eax, 1\n"
		"jz			done%=\n"

		"movsx		ecx, word ptr [esi]\n"
		"cvtsi2ss	xmm0, ecx\n"
		"shufps		xmm0, xmm0, 0x00\n"
		"movlps		[edi+0], xmm0\n"
		"movhps		[edi+8], xmm0\n"

	"done%=:\n"
	::"r"(dest), "r"(src), "m"(numSamples)
    :);
}

/*
============
SSE_UpSample11kHzStereoPCMTo44kHz
============
*/
static void SSE_UpSample11kHzStereoPCMTo44kHz( float *dest, const short *src, const int numSamples ) {
	asm (
		"mov			esi, %1\n"
		"mov			edi, %0\n"

		"mov			eax, %2\n"
		"test		eax, ~1\n"
		"jz			done2%=\n"
		"shl			eax, 1\n"
		"add			esi, eax\n"
		"neg			eax\n"

		/*align		16*/
	"loop2%=:\n"
		"add			edi, 8*4\n"

		"movsx		ecx, word ptr [esi+eax+0]\n"
		"cvtsi2ss	xmm0, ecx\n"

		"movsx		edx, word ptr [esi+eax+2]\n"
		"cvtsi2ss	xmm1, edx\n"

		"unpcklps	xmm0, xmm1\n"

		"movlps		[edi-8*4+0], xmm0\n"
		"movlps		[edi-8*4+8], xmm0\n"
		"movlps		[edi-4*4+0], xmm0\n"
		"movlps		[edi-4*4+8], xmm0\n"

		"add			eax, 2*2\n"
		"jl			loop2%=\n"

	"done2%=:\n"
	::"r"(dest), "r"(src), "m"(numSamples)
    :);
}

/*
============
SSE_UpSample22kHzMonoPCMTo44kHz
============
*/
static void SSE_UpSample22kHzMonoPCMTo44kHz( float *dest, const short *src, const int numSamples ) {
	asm (
		"mov			esi, %1\n"
		"mov			edi, %0\n"

		"mov			eax, %2\n"
		"and			eax, ~1\n"
		"jz			done2%=\n"
		"shl			eax, 1\n"
		"add			esi, eax\n"
		"neg			eax\n"

		/*align		16*/
	"loop2%=:\n"
		"add			edi, 4*4\n"

		"movsx		ecx, word ptr [esi+eax+0]\n"
		"cvtsi2ss	xmm0, ecx\n"

		"movsx		edx, word ptr [esi+eax+2]\n"
		"cvtsi2ss	xmm1, edx\n"

		"shufps		xmm0, xmm1, 0x00\n"
		"movlps		[edi-4*4+0], xmm0\n"
		"movhps		[edi-4*4+8], xmm0\n"

		"add			eax, 2*2\n"
		"jl			loop2%=\n"

	"done2%=:\n"
		"mov			eax, %2\n"
		"and			eax, 1\n"
		"jz			done%=\n"

		"movsx		ecx, word ptr [esi]\n"
		"cvtsi2ss	xmm0, ecx\n"
		"shufps		xmm0, xmm0, 0x00\n"
		"movlps		[edi], xmm0\n"

	"done%=:\n"
	::"r"(dest), "r"(src), "m"(numSamples)
    :);
}

/*
============
SSE_UpSample22kHzStereoPCMTo44kHz
============
*/
static void SSE_UpSample22kHzStereoPCMTo44kHz( float *dest, const short *src, const int numSamples ) {
	asm (
		"mov			esi, %1\n"
		"mov			edi, %0\n"

		"mov			eax, %2\n"
		"test		eax, ~1\n"
		"jz			done2%=\n"
		"shl			eax, 1\n"
		"add			esi, eax\n"
		"neg			eax\n"

		/*align		16*/
	"loop2%=:\n"
		"add		edi, 4*4\n"

		"movsx		ecx, word ptr [esi+eax+0]\n"
		"cvtsi2ss	xmm0, ecx\n"
		"movss		[edi-4*4], xmm0\n"
		"movss		[edi-2*4], xmm0\n"

		"movsx		edx, word ptr [esi+eax+2]\n"
		"cvtsi2ss	xmm1, edx\n"
		"movss		[edi-3*4], xmm1\n"
		"movss		[edi-1*4], xmm1\n"

		"add		eax, 2*2\n"
		"jl			loop2%=\n"

	"done2%=:\n"
	::"r"(dest), "r"(src), "m"(numSamples)
    :);
}

/*
============
SSE_UpSample44kHzMonoPCMTo44kHz
============
*/
static void SSE_UpSample44kHzMonoPCMTo44kHz( float *dest, const short *src, const int numSamples ) {
	asm (
		"mov			esi, %1\n"
		"mov			edi, %0\n"

		"mov			eax, %2\n"
		"and			eax, ~1\n"
		"jz			done2%=\n"
		"shl			eax, 1\n"
		"add			esi, eax\n"
		"neg			eax\n"

		/*align		16*/
	"loop2%=:\n"
		"add			edi, 2*4\n"

		"movsx		ecx, word ptr [esi+eax+0]\n"
		"cvtsi2ss	xmm0, ecx\n"
		"movss		[edi-2*4], xmm0\n"

		"movsx		edx, word ptr [esi+eax+2]\n"
		"cvtsi2ss	xmm1, edx\n"
		"movss		[edi-1*4], xmm1\n"

		"add		eax, 2*2\n"
		"jl			loop2%=\n"

	"done2%=:\n"
		"mov		eax, %2\n"
		"and		eax, 1\n"
		"jz			done%=\n"

		"movsx		ecx, word ptr [esi]\n"
		"cvtsi2ss	xmm0, ecx\n"
		"movss		[edi], xmm0\n"

	"done%=:\n"
	::"r"(dest), "r"(src), "m"(numSamples)
    :);
}

/*
============
CSIMD_SSE::upSamplePCMTo44kHz

  Duplicate samples for 44kHz output.
============
*/
void CSIMD_SSE::upSamplePCMTo44kHz( float *dest, const short *src, const int numSamples, const int kHz, const int numChannels ) {
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
	asm (
		"mov			esi,%1\n"
		"mov			edi, %0\n"
		"movss		xmm7, %3\n"
		"shufps		xmm7, xmm7, 0x00\n"

		"mov			eax, %2\n"
		"and			eax, ~1\n"
		"jz			done2%=\n"
		"shl			eax, 2\n"
		"add			esi, eax\n"
		"neg			eax\n"

		/*align		16*/
	"loop2%=:\n"
		"add			edi, 2*16\n"

		"movss		xmm0, [esi+eax+0]\n"
		"mulss		xmm0, xmm7\n"
		"shufps		xmm0, xmm0, 0x00\n"
		"movlps		[edi-32], xmm0\n"
		"movlps		[edi-24], xmm0\n"

		"movss		xmm1, [esi+eax+4]\n"
		"mulss		xmm1, xmm7\n"
		"shufps		xmm1, xmm1, 0x00\n"
		"movlps		[edi-16], xmm1\n"
		"movlps		[edi- 8], xmm1\n"

		"add			eax, 2*4\n"
		"jl			loop2%=\n"

	"done2%=:\n"
		"mov			eax, %2\n"
		"and			eax, 1\n"
		"jz			done%=\n"

		"movss		xmm0, [esi]\n"
		"mulss		xmm0, xmm7\n"
		"shufps		xmm0, xmm0, 0x00\n"
		"movlps		[edi+0], xmm0\n"
		"movlps		[edi+8], xmm0\n"

	"done%=:\n"
	::"r"(dest), "r"(src), "m"(numSamples),"m"(constant)
    :);
}

/*
============
SSE_UpSample11kHzStereoOGGTo44kHz
============
*/
static void SSE_UpSample11kHzStereoOGGTo44kHz( float *dest, const float * const *src, const int numSamples ) {
	float constant = 32768.0f;
	asm (
		"mov			esi,%1\n"
		"mov			ecx, [esi+0]\n"
		"mov			edx, [esi+4]\n"
		"mov			edi, %0\n"
		"movss		xmm7, %3\n"
		"shufps		xmm7, xmm7, 0x00\n"

		"mov			eax, %2\n"
		"and			eax, ~1\n"
		"jz			done2%=\n"
		"shl			eax, 1\n"
		"add			ecx, eax\n"
		"add			edx, eax\n"
		"neg			eax\n"

		/*align		16*/
	"loop2%=:\n"
		"add			edi, 4*16\n"

		"movlps		xmm0, [ecx+eax]\n"
		"movlps		xmm1, [edx+eax]\n"
		"unpcklps	xmm0, xmm1\n"
		"mulps		xmm0, xmm7\n"
		"movlps		[edi-8*8], xmm0\n"
		"movlps		[edi-7*8], xmm0\n"
		"movlps		[edi-6*8], xmm0\n"
		"movlps		[edi-5*8], xmm0\n"
		"movhps		[edi-4*8], xmm0\n"
		"movhps		[edi-3*8], xmm0\n"
		"movhps		[edi-2*8], xmm0\n"
		"movhps		[edi-1*8], xmm0\n"

		"add			eax, 2*4\n"
		"jl			loop2%=\n"

	"done2%=:\n"
		"mov			eax, %2\n"
		"and			eax, 1\n"
		"jz			done%=\n"

		"movss		xmm0, [ecx]\n"
		"movss		xmm1, [edx]\n"
		"unpcklps	xmm0, xmm1\n"
		"mulps		xmm0, xmm7\n"
		"movlps		[edi+0*8], xmm0\n"
		"movlps		[edi+1*8], xmm0\n"
		"movlps		[edi+2*8], xmm0\n"
		"movlps		[edi+3*8], xmm0\n"

	"done%=:\n"
	::"r"(dest), "r"(src), "m"(numSamples),"m"(constant)
    :);
}

/*
============
SSE_UpSample22kHzMonoOGGTo44kHz
============
*/
static void SSE_UpSample22kHzMonoOGGTo44kHz( float *dest, const float *src, const int numSamples ) {
	float constant = 32768.0f;
	asm (
		"mov			esi,%1\n"
		"mov			edi, %0\n"
		"movss		xmm7, %3\n"
		"shufps		xmm7, xmm7, 0x00\n"

		"mov			eax, %2\n"
		"and			eax, ~1\n"
		"jz			done2%=\n"
		"shl			eax, 2\n"
		"add			esi, eax\n"
		"neg			eax\n"

		/*align		16*/
	"loop2%=:\n"
		"add			edi, 2*8\n"

		"movss		xmm0, [esi+eax+0]\n"
		"movss		xmm1, [esi+eax+4]\n"
		"shufps		xmm0, xmm1, 0x00\n"
		"mulps		xmm0, xmm7\n"
		"movlps		[edi-16], xmm0\n"
		"movhps		[edi- 8], xmm0\n"

		"add			eax, 2*4\n"
		"jl			loop2%=\n"

	"done2%=:\n"
		"mov			eax, %2\n"
		"and			eax, 1\n"
		"jz			done%=\n"

		"movss		xmm0, [esi]\n"
		"mulss		xmm0, xmm7\n"
		"shufps		xmm0, xmm0, 0x00\n"
		"movlps		[edi+0], xmm0\n"

	"done%=:\n"
	::"r"(dest), "r"(src), "m"(numSamples),"m"(constant)
    :);
}

/*
============
SSE_UpSample22kHzStereoOGGTo44kHz
============
*/
static void SSE_UpSample22kHzStereoOGGTo44kHz( float *dest, const float * const *src, const int numSamples ) {
	float constant = 32768.0f;
	asm (
		"mov			esi,%1\n"
		"mov			ecx, [esi+0]\n"
		"mov			edx, [esi+4]\n"
		"mov			edi, %0\n"
		"movss		xmm7, %3\n"
		"shufps		xmm7, xmm7, 0x00\n"

		"mov			eax, %2\n"
		"and			eax, ~1\n"
		"jz			done2%=\n"
		"shl			eax, 1\n"
		"add			ecx, eax\n"
		"add			edx, eax\n"
		"neg			eax\n"

		/*align		16*/
	"loop2%=:\n"
		"add			edi, 2*16\n"

		"movlps		xmm0, [ecx+eax]\n"
		"movlps		xmm1, [edx+eax]\n"
		"unpcklps	xmm0, xmm1\n"
		"mulps		xmm0, xmm7\n"
		"movlps		[edi-4*8], xmm0\n"
		"movlps		[edi-3*8], xmm0\n"
		"movhps		[edi-2*8], xmm0\n"
		"movhps		[edi-1*8], xmm0\n"

		"add			eax, 2*4\n"
		"jl			loop2%=\n"

	"done2%=:\n"
		"mov			eax, %2\n"
		"and			eax, 1\n"
		"jz			done%=\n"

		"movss		xmm0, [ecx]\n"
		"movss		xmm1, [edx]\n"
		"unpcklps	xmm0, xmm1\n"
		"mulps		xmm0, xmm7\n"
		"movlps		[edi+0*8], xmm0\n"
		"movlps		[edi+1*8], xmm0\n"

	"done%=:\n"
	::"r"(dest), "r"(src), "m"(numSamples),"m"(constant)
    :);
}

/*
============
SSE_UpSample44kHzMonoOGGTo44kHz
============
*/
static void SSE_UpSample44kHzMonoOGGTo44kHz( float *dst, const float *src, const int numSamples ) {
	float constant = 32768.0f;
	//KFLOAT_CA( mul, dst, src, constant, numSamples )
	int	pre,post;
    asm(
        KFLOAT_CA( mul )
    :
    : "m"(constant),"m"(numSamples),"m"(dst),"m"(src),"m"(pre),"m"(post)
    : );
}

/*
============
SSE_UpSample44kHzStereoOGGTo44kHz
============
*/
static void SSE_UpSample44kHzStereoOGGTo44kHz( float *dest, const float * const *src, const int numSamples ) {
	float constant = 32768.0f;
	asm (
		"mov			esi,%1\n"
		"mov			ecx, [esi+0]\n"
		"mov			edx, [esi+4]\n"
		"mov			edi, %0\n"
		"movss		xmm7, %3\n"
		"shufps		xmm7, xmm7, 0x00\n"

		"mov			eax, %2\n"
		"and			eax, ~1\n"
		"jz			done2%=\n"
		"shl			eax, 1\n"
		"add			ecx, eax\n"
		"add			edx, eax\n"
		"neg			eax\n"

		/*align		16*/
	"loop2%=:\n"
		"add			edi, 16\n"

		"movlps		xmm0, [ecx+eax]\n"
		"movlps		xmm1, [edx+eax]\n"
		"unpcklps	xmm0, xmm1\n"
		"mulps		xmm0, xmm7\n"
		"movlps		[edi-2*8], xmm0\n"
		"movhps		[edi-1*8], xmm0\n"

		"add			eax, 2*4\n"
		"jl			loop2%=\n"

	"done2%=:\n"
		"mov			eax, %2\n"
		"and			eax, 1\n"
		"jz			done%=\n"

		"movss		xmm0, [ecx]\n"
		"movss		xmm1, [edx]\n"
		"unpcklps	xmm0, xmm1\n"
		"mulps		xmm0, xmm7\n"
		"movlps		[edi+0*8], xmm0\n"

	"done%=:\n"
	::"m"(dest), "m"(src), "m"(numSamples),"m"(constant)
    :);
}

/*
============
CSIMD_SSE::upSampleOGGTo44kHz

  Duplicate samples for 44kHz output.
============
*/
void CSIMD_SSE::upSampleOGGTo44kHz( float *dest, const float * const *ogg, const int numSamples, const int kHz, const int numChannels ) {
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
CSIMD_SSE::mixSoundTwoSpeakerMono
============
*/
void VPCALL CSIMD_SSE::mixSoundTwoSpeakerMono( float *mixBuffer, const float *samples, const int numSamples, const float lastV[2], const float currentV[2] ) {
#if 1

	ALIGN16( float incs[2] );

	SMF_ASSERT( numSamples == MIXBUFFER_SAMPLES );

	incs[0] = ( currentV[0] - lastV[0] ) / MIXBUFFER_SAMPLES;
	incs[1] = ( currentV[1] - lastV[1] ) / MIXBUFFER_SAMPLES;

	asm (
		"mov			eax, %3\n"
		"mov			edi, %0\n"
		"mov			esi, %1\n"
		"shl			eax, 2\n"
		"add			esi, eax\n"
		"neg			eax\n"

		"mov			ecx, %2\n"
		"movlps		xmm6, [ecx]\n"
		"xorps		xmm7, xmm7\n"
		"movhps		xmm7, %4\n"
		"shufps		xmm6, xmm6, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
		"addps		xmm6, xmm7\n"
		"shufps		xmm7, xmm7, 0xEE\n"     /*R_SHUFFLEPS( 2, 3, 2, 3 )=11101110*/
		"addps		xmm7, xmm7\n"

	"loop16%=:\n"
		"add		edi, 4*4*4\n"

		"movaps		xmm0, [esi+eax+0*4*4]\n"
		"movaps		xmm1, xmm0\n"
		"shufps		xmm0, xmm0, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
		"mulps		xmm0, xmm6\n"
		"addps		xmm0, [edi-4*4*4]\n"
		"addps		xmm6, xmm7\n"
		"movaps		[edi-4*4*4], xmm0\n"

		"shufps		xmm1, xmm1, 0xFA\n" /*R_SHUFFLEPS( 2, 2, 3, 3 )=11111010*/
		"mulps		xmm1, xmm6\n"
		"addps		xmm1, [edi-3*4*4]\n"
		"addps		xmm6, xmm7\n"
		"movaps		[edi-3*4*4], xmm1\n"

		"movaps		xmm2, [esi+eax+1*4*4]\n"
		"movaps		xmm3, xmm2\n"
		"shufps		xmm2, xmm2, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
		"mulps		xmm2, xmm6\n"
		"addps		xmm2, [edi-2*4*4]\n"
		"addps		xmm6, xmm7\n"
		"movaps		[edi-2*4*4], xmm2\n"

		"shufps		xmm3, xmm3, 0xFA\n" /*R_SHUFFLEPS( 2, 2, 3, 3 )=11111010*/
		"mulps		xmm3, xmm6\n"
		"addps		xmm3, [edi-1*4*4]\n"
		"addps		xmm6, xmm7\n"
		"movaps		[edi-1*4*4], xmm3\n"

		"add		eax, 2*4*4\n"

		"jl			loop16%=\n"
	:: "m"(mixBuffer), "m"(samples),  "m"(lastV), "m"(MIXBUFFER_SAMPLES),"m"(incs)
    :);

#else

	int i;
	float "incL;
	float "incR;
	float sL0, sL1;
	float sR0, sR1;

	SMF_ASSERT( numSamples == MIXBUFFER_SAMPLES );

	"incL = ( currentV[0] - lastV[0] ) / MIXBUFFER_SAMPLES;
	"incR = ( currentV[1] - lastV[1] ) / MIXBUFFER_SAMPLES;

	sL0 = lastV[0];
	sR0 = lastV[1];
	sL1 = lastV[0] + "incL;
	sR1 = lastV[1] + "incR;

	"incL *= 2;
	"incR *= 2;

	for( i = 0; i < MIXBUFFER_SAMPLES; i += 2 ) {
		mixBuffer[i*2+0] += samples[i+0] * sL0;
		mixBuffer[i*2+1] += samples[i+0] * sR0;
		mixBuffer[i*2+2] += samples[i+1] * sL1;
		mixBuffer[i*2+3] += samples[i+1] * sR1;
		sL0 += "incL;
		sR0 += "incR;
		sL1 += "incL;
		sR1 += "incR;
	}

#endif
}

/*
============
CSIMD_SSE::mixSoundTwoSpeakerStereo
============
*/
void VPCALL CSIMD_SSE::mixSoundTwoSpeakerStereo( float *mixBuffer, const float *samples, const int numSamples, const float lastV[2], const float currentV[2] ) {
#if 1

	ALIGN16( float incs[2] );

	SMF_ASSERT( numSamples == MIXBUFFER_SAMPLES );

	incs[0] = ( currentV[0] - lastV[0] ) / MIXBUFFER_SAMPLES;
	incs[1] = ( currentV[1] - lastV[1] ) / MIXBUFFER_SAMPLES;

	asm (
		"mov			eax, %3\n"
		"mov			edi, %0\n"
		"mov			esi, %1\n"
		"shl			eax, 3\n"
		"add			esi, eax\n"
		"neg			eax\n"

		"mov			ecx, %2\n"
		"movlps		xmm6, [ecx]\n"
		"xorps		xmm7, xmm7\n"
		"movhps		xmm7, %4\n"
		"shufps		xmm6, xmm6, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
		"addps		xmm6, xmm7\n"
		"shufps		xmm7, xmm7, 0XEE\n"  /*R_SHUFFLEPS( 2, 3, 2, 3 )=11101110*/
		"addps		xmm7, xmm7\n"

	"loop16%=:\n"
		"add		edi, 4*4*4\n"

		"movaps		xmm0, [esi+eax+0*4*4]\n"
		"mulps		xmm0, xmm6\n"
		"addps		xmm0, [edi-4*4*4]\n"
		"addps		xmm6, xmm7\n"
		"movaps		[edi-4*4*4], xmm0\n"

		"movaps		xmm2, [esi+eax+1*4*4]\n"
		"mulps		xmm2, xmm6\n"
		"addps		xmm2, [edi-3*4*4]\n"
		"addps		xmm6, xmm7\n"
		"movaps		[edi-3*4*4], xmm2\n"

		"movaps		xmm3, [esi+eax+2*4*4]\n"
		"mulps		xmm3, xmm6\n"
		"addps		xmm3, [edi-2*4*4]\n"
		"addps		xmm6, xmm7\n"
		"movaps		[edi-2*4*4], xmm3\n"

		"movaps		xmm4, [esi+eax+3*4*4]\n"
		"mulps		xmm4, xmm6\n"
		"addps		xmm4, [edi-1*4*4]\n"
		"addps		xmm6, xmm7\n"
		"movaps		[edi-1*4*4], xmm4\n"

		"add			eax, 4*4*4\n"

		"jl			loop16%=\n"
	:: "m"(mixBuffer), "m"(samples),  "m"(lastV), "m"(MIXBUFFER_SAMPLES),"m"(incs)
    :);

#else

	int i;
	float "incL;
	float "incR;
	float sL0, sL1;
	float sR0, sR1;

	SMF_ASSERT( numSamples == MIXBUFFER_SAMPLES );

	"incL = ( currentV[0] - lastV[0] ) / MIXBUFFER_SAMPLES;
	"incR = ( currentV[1] - lastV[1] ) / MIXBUFFER_SAMPLES;

	sL0 = lastV[0];
	sR0 = lastV[1];
	sL1 = lastV[0] + "incL;
	sR1 = lastV[1] + "incR;

	"incL *= 2;
	"incR *= 2;

	for( i = 0; i < MIXBUFFER_SAMPLES; i += 2 ) {
		mixBuffer[i*2+0] += samples[i*2+0] * sL0;
		mixBuffer[i*2+1] += samples[i*2+1] * sR0;
		mixBuffer[i*2+2] += samples[i*2+2] * sL1;
		mixBuffer[i*2+3] += samples[i*2+3] * sR1;
		sL0 += "incL;
		sR0 += "incR;
		sL1 += "incL;
		sR1 += "incR;
	}

#endif
}

/*
============
CSIMD_SSE::mixSoundSixSpeakerMono
============
*/
void VPCALL CSIMD_SSE::mixSoundSixSpeakerMono( float *mixBuffer, const float *samples, const int numSamples, const float lastV[6], const float currentV[6] ) {
#if 1

	ALIGN16( float incs[6] );

	SMF_ASSERT( numSamples == MIXBUFFER_SAMPLES );

	incs[0] = ( currentV[0] - lastV[0] ) / MIXBUFFER_SAMPLES;
	incs[1] = ( currentV[1] - lastV[1] ) / MIXBUFFER_SAMPLES;
	incs[2] = ( currentV[2] - lastV[2] ) / MIXBUFFER_SAMPLES;
	incs[3] = ( currentV[3] - lastV[3] ) / MIXBUFFER_SAMPLES;
	incs[4] = ( currentV[4] - lastV[4] ) / MIXBUFFER_SAMPLES;
	incs[5] = ( currentV[5] - lastV[5] ) / MIXBUFFER_SAMPLES;

	asm (
		"mov		eax, %4\n"
		"mov		edi, %0\n"
		"mov		esi, %1\n"
		"shl		eax, 2\n"
		"add		esi, eax\n"
		"neg		eax\n"

		"mov		ecx, %2\n"
		"movlps		xmm2, [ecx+ 0]\n"
		"movhps		xmm2, [ecx+ 8]\n"
		"movlps		xmm3, [ecx+16]\n"
		"movaps		xmm4, xmm2\n"
		"shufps		xmm3, xmm2, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
		"shufps		xmm4, xmm3, 0x4E\n"   /*R_SHUFFLEPS( 2, 3, 0, 1 )=01001110*/

		"xorps		xmm5, xmm5\n"
		"movhps		xmm5, %3\n"
		"movlps		xmm7, %3+8\n"
		"movhps		xmm7, %3+16\n"
		"addps		xmm3, xmm5\n"
		"addps		xmm4, xmm7\n"
		"shufps		xmm5, xmm7, 0x4E\n"   /*R_SHUFFLEPS( 2, 3, 0, 1 )=01001110*/
		"movaps		xmm6, xmm7\n"
		"shufps		xmm6, xmm5, 0x4E\n"   /*R_SHUFFLEPS( 2, 3, 0, 1 )=01001110*/
		"addps		xmm5, xmm5\n"
		"addps		xmm6, xmm6\n"
		"addps		xmm7, xmm7\n"

	"loop24%=:\n"
		"add			edi, 6*16\n"

		"movaps		xmm0, [esi+eax]\n"

		"movaps		xmm1, xmm0\n"
		"shufps		xmm1, xmm1, 0x00\n"
		"mulps		xmm1, xmm2\n"
		"addps		xmm1, [edi-6*16]\n"
		"addps		xmm2, xmm5\n"
		"movaps		[edi-6*16], xmm1\n"

		"movaps		xmm1, xmm0\n"
		"shufps		xmm1, xmm1, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
		"mulps		xmm1, xmm3\n"
		"addps		xmm1, [edi-5*16]\n"
		"addps		xmm3, xmm6\n"
		"movaps		[edi-5*16], xmm1\n"

		"movaps		xmm1, xmm0\n"
		"shufps		xmm1, xmm1, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
		"mulps		xmm1, xmm4\n"
		"addps		xmm1, [edi-4*16]\n"
		"addps		xmm4, xmm7\n"
		"movaps		[edi-4*16], xmm1\n"

		"movaps		xmm1, xmm0\n"
		"shufps		xmm1, xmm1, 0xAA \n" /*R_SHUFFLEPS( 2, 2, 2, 2 )=10101010*/
		"mulps		xmm1, xmm2\n"
		"addps		xmm1, [edi-3*16]\n"
		"addps		xmm2, xmm5\n"
		"movaps		[edi-3*16], xmm1\n"

		"movaps		xmm1, xmm0\n"
		"shufps		xmm1, xmm1, 0xFA\n" /*R_SHUFFLEPS( 2, 2, 3, 3 )=11111010*/
		"mulps		xmm1, xmm3\n"
		"addps		xmm1, [edi-2*16]\n"
		"addps		xmm3, xmm6\n"
		"movaps		[edi-2*16], xmm1\n"

		"shufps		xmm0, xmm0, 0xFF \n" /*R_SHUFFLEPS( 3, 3, 3, 3 )=11111111*/
		"mulps		xmm0, xmm4\n"
		"addps		xmm0, [edi-1*16]\n"
		"addps		xmm4, xmm7\n"
		"movaps		[edi-1*16], xmm0\n"

		"add			eax, 4*4\n"

		"jl			loop24%=\n"
	:
    : "m"(mixBuffer), "m"(samples), "m"(lastV),"m"(incs),"m"(MIXBUFFER_SAMPLES)
    :);

#else

	int i;
	float sL0, sL1, sL2, sL3, sL4, sL5, sL6, sL7, sL8, sL9, sL10, sL11;
	float "incL0, "incL1, "incL2, "incL3, "incL4, "incL5;

	SMF_ASSERT( numSamples == MIXBUFFER_SAMPLES );

	"incL0 = ( currentV[0] - lastV[0] ) / MIXBUFFER_SAMPLES;
	"incL1 = ( currentV[1] - lastV[1] ) / MIXBUFFER_SAMPLES;
	"incL2 = ( currentV[2] - lastV[2] ) / MIXBUFFER_SAMPLES;
	"incL3 = ( currentV[3] - lastV[3] ) / MIXBUFFER_SAMPLES;
	"incL4 = ( currentV[4] - lastV[4] ) / MIXBUFFER_SAMPLES;
	"incL5 = ( currentV[5] - lastV[5] ) / MIXBUFFER_SAMPLES;

	sL0  = lastV[0];
	sL1  = lastV[1];
	sL2  = lastV[2];
	sL3  = lastV[3];
	sL4  = lastV[4];
	sL5  = lastV[5];

	sL6  = lastV[0] + "incL0;
	sL7  = lastV[1] + "incL1;
	sL8  = lastV[2] + "incL2;
	sL9  = lastV[3] + "incL3;
	sL10 = lastV[4] + "incL4;
	sL11 = lastV[5] + "incL5;

	"incL0 *= 2;
	"incL1 *= 2;
	"incL2 *= 2;
	"incL3 *= 2;
	"incL4 *= 2;
	"incL5 *= 2;

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

		sL0  += "incL0;
		sL1  += "incL1;
		sL2  += "incL2;
		sL3  += "incL3;

		sL4  += "incL4;
		sL5  += "incL5;
		sL6  += "incL0;
		sL7  += "incL1;

		sL8  += "incL2;
		sL9  += "incL3;
		sL10 += "incL4;
		sL11 += "incL5;
	}

#endif
}

/*
============
CSIMD_SSE::mixSoundSixSpeakerStereo
============
*/
void VPCALL CSIMD_SSE::mixSoundSixSpeakerStereo( float *mixBuffer, const float *samples, const int numSamples, const float lastV[6], const float currentV[6] ) {
#if 1

	ALIGN16( float incs[6] );

	SMF_ASSERT( numSamples == MIXBUFFER_SAMPLES );
	SMF_ASSERT( SPEAKER_RIGHT == 1 );
	SMF_ASSERT( SPEAKER_BACKRIGHT == 5 );
    float * lastVPtr = (float*) &lastV;
    float * incsPtr = (float *) &incs;

    incs[0] = ( currentV[0] - lastV[0] ) / MIXBUFFER_SAMPLES;
	incs[1] = ( currentV[1] - lastV[1] ) / MIXBUFFER_SAMPLES;
	incs[2] = ( currentV[2] - lastV[2] ) / MIXBUFFER_SAMPLES;
	incs[3] = ( currentV[3] - lastV[3] ) / MIXBUFFER_SAMPLES;
	incs[4] = ( currentV[4] - lastV[4] ) / MIXBUFFER_SAMPLES;
	incs[5] = ( currentV[5] - lastV[5] ) / MIXBUFFER_SAMPLES;

	asm (
		"mov		eax, %4\n"  //MIXBUFFER_SAMPLES
		"mov		edi, %0\n" //mixbuffer
		"mov		esi, %1\n" //samples
		"shl		eax, 2\n"
		"add		esi, eax\n"
		"neg		eax\n"

		"mov		ecx, %2\n"   //lastVptr
		"movlps		xmm2, [ecx+ 0]\n"
		"movhps		xmm2, [ecx+ 8]\n"
		"movlps		xmm3, [ecx+16]\n"
		"movaps		xmm4, xmm2\n"
		"shufps		xmm3, xmm2, 0x44\n"		/*R_SHUFFLEPS( 0, 1, 0, 1 )=1000100*/
		"shufps		xmm4, xmm3, 0x4E\n"   /*R_SHUFFLEPS( 2, 3, 0, 1 )=01001110*/

		"xorps		xmm5, xmm5\n"
		"mov        ebx,  %3\n"
		"movhps		xmm5, [ebx]\n"
		"movlps		xmm7, [ebx+8]\n"
		"movhps		xmm7, [ebx+16]\n"
		"addps		xmm3, xmm5\n"
		"addps		xmm4, xmm7\n"
		"shufps		xmm5, xmm7, 0x4E\n"   /*R_SHUFFLEPS( 2, 3, 0, 1 )=01001110*/
		"movaps		xmm6, xmm7\n"
		"shufps		xmm6, xmm5, 0x4E\n"   /*R_SHUFFLEPS( 2, 3, 0, 1 )=01001110*/
		"addps		xmm5, xmm5\n"
		"addps		xmm6, xmm6\n"
		"addps		xmm7, xmm7\n"

	"loop24%=:\n"
		"add			edi, 6*16\n"

		"movaps		xmm0, [esi+eax]\n"

		"movaps		xmm1, xmm0\n"
		"shufps		xmm1, xmm1, 0x00\n"
		"mulps		xmm1, xmm2\n"
		"addps		xmm1, [edi-6*16]\n"
		"addps		xmm2, xmm5\n"
		"movaps		[edi-6*16], xmm1\n"

		"movaps		xmm1, xmm0\n"
		"shufps		xmm1, xmm1, 0x50 \n" /*R_SHUFFLEPS( 0, 0, 1, 1 )=01010000*/
		"mulps		xmm1, xmm3\n"
		"addps		xmm1, [edi-5*16]\n"
		"addps		xmm3, xmm6\n"
		"movaps		[edi-5*16], xmm1\n"

		"movaps		xmm1, xmm0\n"
		"shufps		xmm1, xmm1, 0x55\n"  /* R_SHUFFLEPS( 1, 1, 1, 1 )=01010101 */
		"mulps		xmm1, xmm4\n"
		"addps		xmm1, [edi-4*16]\n"
		"addps		xmm4, xmm7\n"
		"movaps		[edi-4*16], xmm1\n"

		"movaps		xmm1, xmm0\n"
		"shufps		xmm1, xmm1, 0xAA \n" /*R_SHUFFLEPS( 2, 2, 2, 2 )=10101010*/
		"mulps		xmm1, xmm2\n"
		"addps		xmm1, [edi-3*16]\n"
		"addps		xmm2, xmm5\n"
		"movaps		[edi-3*16], xmm1\n"

		"movaps		xmm1, xmm0\n"
		"shufps		xmm1, xmm1, 0xFA\n" /*R_SHUFFLEPS( 2, 2, 3, 3 )=11111010*/
		"mulps		xmm1, xmm3\n"
		"addps		xmm1, [edi-2*16]\n"
		"addps		xmm3, xmm6\n"
		"movaps		[edi-2*16], xmm1\n"

		"shufps		xmm0, xmm0, 0xFF \n" /*R_SHUFFLEPS( 3, 3, 3, 3 )=11111111*/
		"mulps		xmm0, xmm4\n"
		"addps		xmm0, [edi-1*16]\n"
		"addps		xmm4, xmm7\n"
		"movaps		[edi-1*16], xmm0\n"

		"add			eax, 4*4\n"

		"jl			loop24%=\n"
	:
    : "m"(mixBuffer), "m"(samples), "m"(lastVPtr),"m"(incsPtr),"m"(MIXBUFFER_SAMPLES)
    :"eax","ebx","ecx","edx","edi","esi");

#else

	int i;
	float sL0, sL1, sL2, sL3, sL4, sL5, sL6, sL7, sL8, sL9, sL10, sL11;
	float "incL0, "incL1, "incL2, "incL3, "incL4, "incL5;

	SMF_ASSERT( numSamples == MIXBUFFER_SAMPLES );
	SMF_ASSERT( SPEAKER_RIGHT == 1 );
	SMF_ASSERT( SPEAKER_BACKRIGHT == 5 );

	"incL0 = ( currentV[0] - lastV[0] ) / MIXBUFFER_SAMPLES;
	"incL1 = ( currentV[1] - lastV[1] ) / MIXBUFFER_SAMPLES;
	"incL2 = ( currentV[2] - lastV[2] ) / MIXBUFFER_SAMPLES;
	"incL3 = ( currentV[3] - lastV[3] ) / MIXBUFFER_SAMPLES;
	"incL4 = ( currentV[4] - lastV[4] ) / MIXBUFFER_SAMPLES;
	"incL5 = ( currentV[5] - lastV[5] ) / MIXBUFFER_SAMPLES;

	sL0  = lastV[0];
	sL1  = lastV[1];
	sL2  = lastV[2];
	sL3  = lastV[3];
	sL4  = lastV[4];
	sL5  = lastV[5];

	sL6  = lastV[0] + "incL0;
	sL7  = lastV[1] + "incL1;
	sL8  = lastV[2] + "incL2;
	sL9  = lastV[3] + "incL3;
	sL10 = lastV[4] + "incL4;
	sL11 = lastV[5] + "incL5;

	"incL0 *= 2;
	"incL1 *= 2;
	"incL2 *= 2;
	"incL3 *= 2;
	"incL4 *= 2;
	"incL5 *= 2;

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

		sL0  += "incL0;
		sL1  += "incL1;
		sL2  += "incL2;
		sL3  += "incL3;

		sL4  += "incL4;
		sL5  += "incL5;
		sL6  += "incL0;
		sL7  += "incL1;

		sL8  += "incL2;
		sL9  += "incL3;
		sL10 += "incL4;
		sL11 += "incL5;
	}

#endif
}

/*
============
CSIMD_SSE::mixedSoundToSamples
============
*/
void VPCALL CSIMD_SSE::mixedSoundToSamples( short *samples, const float *mixBuffer, const int numSamples ) {
#if 1

	SMF_ASSERT( ( numSamples % MIXBUFFER_SAMPLES ) == 0 );

	asm (

		"mov			eax, %2\n"
		"mov			edi, %0\n"
		"mov			esi, %1\n"
		"shl			eax, 2\n"
		"add			edi, eax\n"
		"neg			eax\n"

	"loop16%=:\n"

		"movaps		xmm0, [edi+eax+0*16]\n"
		"movaps		xmm2, [edi+eax+1*16]\n"
		"movaps		xmm4, [edi+eax+2*16]\n"
		"movaps		xmm6, [edi+eax+3*16]\n"

		"add			esi, 4*4*2\n"

		"movhlps		xmm1, xmm0\n"
		"movhlps		xmm3, xmm2\n"
		"movhlps		xmm5, xmm4\n"
		"movhlps		xmm7, xmm6\n"

		"prefetchnta	[edi+eax+64]\n"

		"cvtps2pi	mm0, xmm0\n"
		"cvtps2pi	mm2, xmm2\n"
		"cvtps2pi	mm4, xmm4\n"
		"cvtps2pi	mm6, xmm6\n"

		"prefetchnta	[edi+eax+128]\n"

		"cvtps2pi	mm1, xmm1\n"
		"cvtps2pi	mm3, xmm3\n"
		"cvtps2pi	mm5, xmm5\n"
		"cvtps2pi	mm7, xmm7\n"

		"add		eax, 4*16\n"

		"packssdw	mm0, mm1\n"
		"packssdw	mm2, mm3\n"
		"packssdw	mm4, mm5\n"
		"packssdw	mm6, mm7\n"

		"movq		[esi-4*4*2], mm0\n"
		"movq		[esi-3*4*2], mm2\n"
		"movq		[esi-2*4*2], mm4\n"
		"movq		[esi-1*4*2], mm6\n"

		"jl			loop16%=\n"

		"emms\n"
	:
    : "m"(mixBuffer), "m"(samples) ,"m"(numSamples)
    :"eax","esi","edi");

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

#define  HIPREC

void VPCALL CSIMD_SSE::vector3D_Sum(CVec3D* pOut, const CVec3D* pIn)
{
	/* SSE/SSE2/SSE3 Implementation */
	 asm(
                "mov ecx, %0\n"
                "movups xmm0, [%1]\n"                 // Move unaligned vectors to SSE regs
                "movups xmm1, [%0]\n"
				"andps xmm0, %2\n" //zera o elemento w
                "andps xmm1, %2\n"  //remove w component
				"addps xmm0, xmm1\n"                   // add vector elements
                "movlpd	[ecx+ 0], xmm0\n"  //retira os tres bytes
				"shufps xmm0,xmm0, 0x02\n" //shuffle colocar nos bits 0-31 o 3o byte
				"movss	[ecx+ 8], xmm0\n"
    :
	: "r" (pOut), "r" (pIn), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);

}

void VPCALL CSIMD_SSE::vector3D_SumOf(CVec3D* pOut, const CVec3D* pIn1, const CVec3D* pIn2)
{
	/* SSE/SSE2/SSE3 Implementation */
asm(
                "mov eax, %0\n"                       // Load pointers into CPU regs
                "mov ebx, %1\n"
				"mov ecx, %2\n"
                "movups xmm0, [eax]\n"                 // Move unaligned vectors to SSE regs
                "movups xmm1, [ebx]\n"
				"andps xmm0, %3\n" //zera o elemento w
                "andps xmm1, %3\n"  //re"move w component
				"addps xmm0, xmm1\n"                   // add vector elements
                "movlpd	[ecx+ 0], xmm0\n"  //retira os tres bytes
				"shufps xmm0,xmm0, 0x02\n" //shuffle colocar nos bits 0-31 o 3o byte
				"movss	[ecx+ 8], xmm0\n"

	:
	: "r" (pIn1), "r" (pIn2), "r" (pOut), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);
}

//Subtracts pIn from pOut and stores the result in pOut. POut = pIn-pOut
void VPCALL CSIMD_SSE::vector3D_Diff(CVec3D* pOut, CVec3D* pIn)
{
	/* SSE/SSE2/SSE3 Implementation */
		 asm(
                "mov eax, %0\n"                      // Load pointers into CPU regs
                "mov ebx, %1\n"
				"mov ecx, %0\n"
                "movups xmm0, [eax]\n"                 // Move unaligned vectors to SSE regs
                "movups xmm1, [ebx]\n"
				"andps xmm0, %2\n" //zera o elemento w
                "andps xmm1, %2\n"  //re"move w component
				"subps xmm0, xmm1\n"                   // xmm0 = xmm0 - xmm1
                "movlpd	[ecx+ 0], xmm0\n"  //retira os tres bytes
				"shufps xmm0,xmm0, 0x02\n" //shuffle colocar nos bits 0-31 o 3o byte
				"movss	[ecx+ 8], xmm0\n"
        :
	: "r" (pOut), "r" (pIn), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);
}


//Subtracts pIn2 from pIn1 and stores the result in pOut.  pout = pIn1-Pin2
void VPCALL CSIMD_SSE::vector3D_DiffOf(CVec3D* pOut, const CVec3D* pIn1, const CVec3D* pIn2)
{

	asm(
                "mov eax, %0\n"                       // Load pointers into CPU regs
                "mov ebx, %1\n"
				"mov ecx, %2\n"
                "movups xmm0, [eax]\n"                 // Move unaligned vectors to SSE regs
                "movups xmm1, [ebx]\n"
				"andps xmm0, %3\n" //zera o elemento w
                "andps xmm1, %3\n"  //re"move w component
				"subps xmm0, xmm1\n"                   // xmm0 = xmm0 - xmm1
                "movlpd	[ecx+ 0], xmm0\n"  //retira os tres bytes
				"shufps xmm0,xmm0, 0x02\n" //shuffle colocar nos bits 0-31 o 3o byte
				"movss	[ecx+ 8], xmm0\n"
	:
	: "r" (pIn1), "r" (pIn2), "r" (pOut), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);
}
//Scales the components of pOut by scalar and stores the result in pOut.
void VPCALL CSIMD_SSE::vector3D_Scale(CVec3D* pOut, float scalar)
{
	float *test=&scalar;
	asm(
			"mov eax, %0\n"                       // Load pointers into CPU regs
			"mov ebx, %1\n"
			"mov ecx, %0\n"
            "movss xmm0, [ebx]\n"
			"movups xmm1, [eax]\n"
			"andps xmm1, %2\n"  //re"move w component
			"shufps xmm0,xmm0,0x00\n"	//repete o valor da scalar para os demais bytes de xmm0
			"mulps xmm0,xmm1\n"
		    "movlpd	[ecx+ 0], xmm0\n"  //retira os tres bytes
			"shufps xmm0,xmm0, 0x02\n" //shuffle colocar nos bits 0-31 o 3o byte
			"movss	[ecx+ 8], xmm0\n"
	:
	: "r" (pOut), "m" (test), "m"(_SIMDx86_float_SSE_NO_W_MASK)
	);

}
//Scales the components of pIn by scalar and stores the result in pOut.
void VPCALL CSIMD_SSE::vector3D_ScaleOf(CVec3D* pOut, const CVec3D* pIn, float scalar)
{
	float *test=&scalar;
	asm(
			"mov ecx, %0\n"
            "movss xmm0, [%2]\n"
			"movups xmm1, [%1]\n"
			"andps xmm1, %3\n"  //re"move w component /* xmm0 = x | y | z | 0.0f */
			"shufps xmm0,xmm0, 0x00\n"	//repete o valor da scalar para os demais bytes de xmm0
			"mulps xmm0,xmm1\n"
		    "movlpd	[ecx+ 0], xmm0\n"  //retira os tres bytes
			"shufps xmm0,xmm0, 0x02\n" //shuffle colocar nos bits 0-31 o 3o byte
			"movss	[ecx+ 8], xmm0\n"
	:
	: "r" (pOut), "r" (pIn),  "m" (scalar), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);



}
//Calculates the dot product (known as the inner product) of pSrc1 and pSrc2, that represents the cosine of the angle between them, if both are normalized.
float  VPCALL CSIMD_SSE::vector3D_Dot(const CVec3D* pSrc1, const CVec3D* pSrc2)
{

		/* SSE/SSE2 Implementation */
		float dummy;
		float *min = &dummy;
		asm(
			"mov eax, %0\n"                       // Load pointers into CPU regs and multiply
			"mov ebx, %3\n"
			"movups xmm0, [eax]\n"
			"movups xmm1, [ebx]\n"
			"andps xmm0, %1\n"  //remove w component /* xmm0 = x | y | z | 0.0f */
			"andps xmm1, %1\n"  //remove w component /* xmm1 = x | y | z | 0.0f */
			"mulps xmm1, xmm0\n"
			/* Shift data around in the registers (I loves me some haddps right now!!) */
			"movhlps xmm2, xmm1\n"    /* xmm2 = ?   | ?   | 0   |  z's */
		    "addss   xmm2, xmm1\n"    /* xmm2 = ?   | ?   | 0   | z's + x's */
			"shufps  xmm1, xmm1,0x55 \n"  /*R_SHUFFLEPS( 1,1,1,1 )*/ /*0x55=01010101 / xmm1 = y's | y's | y's | y's */
			"addss   xmm2, xmm1\n" /* xmm2 = ?   | ?   | ?   | x's+y's+z's */
			//coloca 0-31 de xmm2 em dummy
			"mov			esi, %2\n"
			"movss		[esi], xmm2\n"
            :
            :"r"(pSrc1), "m" (_SIMDx86_float_SSE_NO_W_MASK),"m" (min),"r"(pSrc2)
            :
		);

		return dummy;


}

//Returns the squared length of pSrc1. This can be useful when doing sphere-sphere collision tests where the exact distance is not needed and a squared one can be used, for example.
//It is signifigantly faster than using the square root version SIMDx86Vector_Length().

float  VPCALL CSIMD_SSE::vector3D_LengthSq(const CVec3D* pSrc1)
{

		float dummy;


		/* SSE/SSE2 Implementation */
		float *min = &dummy;
		asm(
			"mov eax, %0\n"                       // Load pointers into CPU regs and multiply
			"movups xmm0, [eax]\n"
			"andps xmm0, %1\n"  //remove w component /* xmm0 = x | y | z | 0.0f */
			"mulps xmm0, xmm0\n"    /* xmm0 = xmm0^2 */

			/* Shift data around in the registers (I loves me some haddps right now!!) */
			"movhlps xmm1, xmm0\n"		/* xmm1 = ?   | ?   | 0   |  z's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
			"shufps  xmm0, xmm0, 0x55\n"  /*R_SHUFFLEPS( 1,1,1,1 )*/  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */

		    //coloca 0-31 de xmm1 em dummy
			"mov			esi, %2\n"
			"movss		[esi], xmm1\n"
            :
            :"r"(pSrc1), "m" (_SIMDx86_float_SSE_NO_W_MASK),"m" (min)
            :
		);

		/* USE_SSE == 1 ||  USE_SSE == 2 */

		return dummy;

}

//Returns the length (magnitude) of pSrc1. Normalized vectors (e.g. the output of SIMDx86Vector_Normalize() ) are of unit length, or 1.0)
float  VPCALL CSIMD_SSE::vector3D_Length(const CVec3D* pSrc1)
{

		float dummy;
		/* SSE/SSE2 Implementation */
		float *min = &dummy;
		asm(
			"mov eax, %0\n"                       // Load pointers into CPU regs and multiply
			"movups xmm0, [eax]\n"
			"andps xmm0, %1\n"  //re"move w component /* xmm0 = x | y | z | 0.0f */
			"mulps xmm0, xmm0\n"    /* xmm0 = xmm0^2 */

				/* Shift data around in the registers (I loves me some haddps right now!!) */
			"movhlps xmm1, xmm0\n"		/* xmm1 = ?   | ?   | 0   |  z's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
			"shufps  xmm0, xmm0, 0x55\n"  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			#ifdef HIPREC
			/* Full square root */
			"sqrtss xmm1, xmm1\n"
			#else
			/* rcp( rSqrt(value) ) (This may be very inaccurate) */
			"rsqrtss xmm1, xmm1\n"
			"rcpss xmm1, xmm1\n"
			#endif

		    //coloca 0-31 de xmm1 em dummy
			"mov	esi, %2\n"
			"movss	[esi], xmm1\n"


		:
		: "r" (pSrc1), "m" (_SIMDx86_float_SSE_NO_W_MASK),"m" (min)
		);




		return dummy;


}

void VPCALL CSIMD_SSE::vector3D_Cross(CVec3D* pLeft, const CVec3D* pRight)
{

	asm(
			"mov eax, %0\n"
			"mov ebx, %1\n" // Load pointers into CPU regs and multiply
			"mov ecx, %0\n"
			"movups xmm0, [eax]\n"  /* xmm0 = pLeft */
			"movups xmm1, [ebx]\n"  /* xmm1 = pRight */
			"andps xmm0, %2\n"  //re"move w component /* xmm0 = x | y | z | 0.0f */
			"andps xmm1, %2\n"  //re"move w component /* xmm1 = x | y | z | 0.0f */
			"movaps xmm2,xmm0\n"  /* xmm0=xmm2 = pLeft */
			"movaps xmm3,xmm1\n"  /* xmm1=xmm3 = pRight */
			/*
			Shuffles
			w | x | z | y  == 1100 1001 = 0xC9
			w | y | x | z  == 1101 0010 = 0xD2
			*/
			"shufps xmm0,xmm0, 0xC9\n"  /* 0xC9=11001001 */
			"shufps xmm1,xmm1, 0xD2\n"  /* 0xD2=11010010 */
			"shufps xmm2,xmm2, 0xD2\n"  /* 0xD2=11010010 */
			"shufps xmm3,xmm3, 0xC9\n"  /* 0xC9=11001001 */
			/* Multiply columns 1&2 and 3&4 */
			"mulps xmm1, xmm0\n"
			"mulps xmm3, xmm2\n"

			/* Subtract products to get the cross product! */
			"subps xmm1, xmm3\n"

			/* Store */
	        "movlpd	[ecx+ 0], xmm1\n"  //retira os tres bytes
			"shufps xmm1,xmm1, 0x02\n" //shuffle colocar nos bits 0-31 o 3o byte
			"movss	[ecx+ 8], xmm1\n"

	:
	: "r" (pLeft), "r" (pRight), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);

}

void VPCALL CSIMD_SSE::vector3D_CrossOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight)
{

	asm(
			"mov eax, %0\n"
			"mov ebx, %1\n" // Load pointers into CPU regs and multiply
			"mov ecx, %3\n"
			"movups xmm0, [eax]\n"  /* xmm0 = pLeft */
			"movups xmm1, [ebx]\n"  /* xmm1 = pRight */
			"andps xmm0, %2\n"  //re"move w component /* xmm0 = x | y | z | 0.0f */
			"andps xmm1, %2\n"  //re"move w component /* xmm1 = x | y | z | 0.0f */
			"movaps xmm2,xmm0\n"  /* xmm2 = pLeft */
			"movaps xmm3,xmm1\n"  /* xmm3 = pRight */
			/*
			Shuffles
			w | x | z | y  == 1100 1001 = 0xC9
			w | y | x | z  == 1101 0010 = 0xD2
			*/
			"shufps xmm0,xmm0, 0xC9\n"  /* 0xC9=11001001 */
			"shufps xmm1,xmm1, 0xD2\n"  /* 0xD2=11010010 */
			"shufps xmm2,xmm2, 0xD2\n"  /* 0xD2=11010010 */
			"shufps xmm3,xmm3, 0xC9\n"  /* 0xC9=11001001 */
			/* Multiply columns 1&2 and 3&4 */
			"mulps xmm1, xmm0\n"
			"mulps xmm3, xmm2\n"

			/* Subtract products to get the cross product! */
			"subps xmm1, xmm3\n"

			/* Store */
	        "movlpd	[ecx+ 0], xmm1\n"  //retira os tres bytes
			"shufps xmm1,xmm1, 0x02\n" //shuffle colocar nos bits 0-31 o 3o byte
			"movss	[ecx+ 8], xmm1\n"

	:
	: "r" (pLeft), "r" (pRight), "m" (_SIMDx86_float_SSE_NO_W_MASK),"r"(pOut)
	);


}

void VPCALL CSIMD_SSE::vector3D_Normalize(CVec3D* pVec)
{


		/* SSE/SSE2 Implementation */
		asm(
			"mov eax, %0\n"
			"mov ecx, %0\n"
			"movups xmm0, [eax]\n"  /* xmm0 = pVec */
			"andps xmm0, %1\n"  //re"move w component /* xmm0 = x | y | z | 0.0f */
			"movaps xmm7, xmm0\n"			/* Save for the division by length */
		    "mulps xmm0, xmm0\n"    /* xmm0 = xmm0^2 */

				/* Shift data around in the registers (I loves me some haddps right now!!) */
			"movhlps xmm1, xmm0\n"		/* xmm1 = ?   | ?   | 0   |  z's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
			"shufps  xmm0, xmm0, 0x55\n"  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */


			#ifdef HIPREC

			/* Divide by magnitude */
			"sqrtss xmm1, xmm1\n"			/* xmm1 = ? | ? | ? | mag(pVec) */
			"shufps xmm1, xmm1,	0x00\n" /*$0x00=00000000,  xmm1 = mag(pVec) | mag(pVec) | mag(pVec) | mag(pVec) */
			"divps xmm7, xmm1\n"			/* xmm1 = ? | z/mag | y/mag | x/mag */

			#else

			/* Multiply by reciprocal magnitude */
			"rsqrtss xmm1 ,xmm1\n"			/* xmm1 = ? | ? | ? | rcp(mag) */
			"shufps xmm1, xmm1, 0x00\n"	        /*$0x00=00000000,  xmm3 = rcp(mag) | rcp(mag) | rcp(mag) | rcp(mag) */
			"mulps xmm7, xmm1\n"			/* xmm7 = ? | z/mag | y/mag | x/mag */

			#endif

			/* Store */
	        "movlpd	[ecx+ 0], xmm7\n"  //retira os tres bytes
			"shufps xmm7,xmm7, 0x02\n" //shuffle colocar nos bits 0-31 o 3o byte
			"movss	[ecx+ 8], xmm7\n"

			:
		: "r" (pVec), "m" (_SIMDx86_float_SSE_NO_W_MASK)
		);


}

void VPCALL CSIMD_SSE::vector3D_NormalizeOf(CVec3D* pOut, const CVec3D* pVec)
{

		/* SSE/SSE2 Implementation */
		asm(
			"mov eax, %0\n"
			"mov ecx, %1\n"
			"movups xmm0, [eax]\n"  /* xmm0 = pVec */
			"andps xmm0, %2\n"  //re"move w component /* xmm0 = x | y | z | 0.0f */
			"movaps xmm7, xmm0\n"			/* Save for the division by length */
		    "mulps xmm0, xmm0\n"   /* xmm0 = xmm0^2 */

				/* Shift data around in the registers (I loves me some haddps right now!!) */
			"movhlps xmm1, xmm0\n"		/* xmm1 = ?   | ?   | 0   |  z's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
			"shufps  xmm0, xmm0, 0x55\n"  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */


			#ifdef HIPREC

			/* Divide by magnitude */
			"sqrtss xmm1, xmm1\n"			/* xmm1 = ? | ? | ? | mag(pVec) */
			"shufps xmm1, xmm1,	0x00\n" /*$0x00=00000000,  xmm1 = mag(pVec) | mag(pVec) | mag(pVec) | mag(pVec) */
			"divps xmm7, xmm1\n"			/* xmm1 = ? | z/mag | y/mag | x/mag */

			#else

			/* Multiply by reciprocal magnitude */
			"rsqrtss xmm1 ,xmm1\n"			/* xmm1 = ? | ? | ? | rcp(mag) */
			"shufps xmm1, xmm1, $0x00\n"	        /*$0x00=00000000,  xmm3 = rcp(mag) | rcp(mag) | rcp(mag) | rcp(mag) */
			"mulps xmm7, xmm1\n"			/* xmm7 = ? | z/mag | y/mag | x/mag */

			#endif

			/* Store */
	        "movlpd	[ecx+ 0], xmm7\n"  //retira os tres bytes
			"shufps xmm7,xmm7, 0x02\n" //shuffle colocar nos bits 0-31 o 3o byte
			"movss	[ecx+ 8], xmm7\n"

		:
		: "r" (pVec), "r" (pOut),"m" (_SIMDx86_float_SSE_NO_W_MASK)
		);

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
	Optimize Distance() -- It currently uses already implemented functions
*/

void VPCALL CSIMD_SSE::vector3D_AlignedSum(CVec3D* pOut, const CVec3D* pIn)
{


	/* SSE/SSE2/SSE3 Implementation */
	 asm(
                "mov ecx, %0\n"
                "movaps xmm0, [%1]\n"                 // Move unaligned vectors to SSE regs
                "movaps xmm1, [%0]\n"
				"andps xmm0, %2\n" //zera o elemento w
                "andps xmm1, %2\n"  //remove w component
				"addps xmm0, xmm1\n"                   // add vector elements
                "movlpd	[ecx+ 0], xmm0\n"  //retira os tres bytes
				"shufps xmm0,xmm0, 0x02\n" //shuffle colocar nos bits 0-31 o 3o byte
				"movss	[ecx+ 8], xmm0\n"
    :
	: "r" (pOut), "r" (pIn), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);



}

void VPCALL CSIMD_SSE::vector3D_AlignedSumOf(CVec3D* pOut, const CVec3D* pIn1, const CVec3D* pIn2)
{
		/* SSE/SSE2/SSE3 Implementation */
		asm(
                //"mov eax, %0\n"                       // Load pointers into CPU regs
                //"mov ebx, %1\n"
				//"mov ecx, %2\n"
                "movaps xmm0, [eax]\n"                 // Move unaligned vectors to SSE regs
                "movaps xmm1, [ebx]\n"
				"andps xmm0, %3\n" //zera o elemento w
                "andps xmm1, %3\n"  //re"move w component
				"addps xmm0, xmm1\n"                   // add vector elements
                "movlpd	[ecx+ 0], xmm0\n"  //retira os tres bytes
				"shufps xmm0,xmm0, 0x02\n" //shuffle colocar nos bits 0-31 o 3o byte
				"movss	[ecx+ 8], xmm0\n"
	:
	: "a" (pIn1), "b" (pIn2), "c" (pOut), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);


}


void VPCALL CSIMD_SSE::vector3D_AlignedDiff(CVec3D* pOut, CVec3D* pIn)
{

		/* SSE/SSE2/SSE3 Implementation */
		asm(
                "mov eax, %0\n"                      // Load pointers into CPU regs
                "mov ebx, %1\n"
				"mov ecx, %0\n"
                "movaps xmm0, [eax]\n"                 // Move unaligned vectors to SSE regs
                "movaps xmm1, [ebx]\n"
				"andps xmm0, %2\n" //zera o elemento w
                "andps xmm1, %2\n"  //re"move w component
				"subps xmm0, xmm1\n"                   // xmm0 = xmm0 - xmm1
                "movlpd	[ecx+ 0], xmm0\n"  //retira os tres bytes
				"shufps xmm0,xmm0, 0x02\n" //shuffle colocar nos bits 0-31 o 3o byte
				"movss	[ecx+ 8], xmm0\n"
  	:
	: "r" (pOut), "r" (pIn), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);

}

void VPCALL CSIMD_SSE::vector3D_AlignedDiffOf(CVec3D* pOut, const CVec3D* pIn1, const CVec3D* pIn2)
{

		asm(
                "mov eax, %0\n"                       // Load pointers into CPU regs
                "mov ebx, %1\n"
				"mov ecx, %2\n"
                "movaps xmm0, [eax]\n"                 // Move unaligned vectors to SSE regs
                "movaps xmm1, [ebx]\n"
				"andps xmm0, %3\n" //zera o elemento w
                "andps xmm1, %3\n"  //re"move w component
				"subps xmm0, xmm1\n"                   // xmm0 = xmm0 - xmm1
                "movlpd	[ecx+ 0], xmm0\n"  //retira os tres bytes
				"shufps xmm0,xmm0, 0x02\n" //shuffle colocar nos bits 0-31 o 3o byte
				"movss	[ecx+ 8], xmm0\n"
	:
	: "r" (pIn1), "r" (pIn2), "r" (pOut), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);



}

void VPCALL CSIMD_SSE::vector3D_AlignedScale(CVec3D* pOut, float scalar)
{

	float *test=&scalar;
	asm(
			"mov eax, %0\n"                       // Load pointers into CPU regs
			"mov ebx, %1\n"
			"mov ecx, %0\n"
            "movss xmm0, [ebx]\n"
			"movaps xmm1, [eax]\n"
			"andps xmm1, %2\n"  //re"move w component
			"shufps xmm0,xmm0,0x00\n"	//repete o valor da scalar para os demais bytes de xmm0
			"mulps xmm0,xmm1\n"
		    "movlpd	[ecx+ 0], xmm0\n"  //retira os tres bytes
			"shufps xmm0,xmm0, 0x02\n" //shuffle colocar nos bits 0-31 o 3o byte
			"movss	[ecx+ 8], xmm0\n"
	:
	: "r" (pOut), "m" (test), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);



}

void VPCALL CSIMD_SSE::vector3D_AlignedScaleOf(CVec3D* pOut, const CVec3D* pIn, float scalar)
{

	float *test=&scalar;
	asm(
			"mov eax, %1\n"                       // Load pointers into CPU regs
			"mov ebx, %2\n"
			"mov ecx, %0\n"
            "movss xmm0, [ebx]\n"
			"movaps xmm1, [eax]\n"
			"andps xmm1, %3\n"  //re"move w component /* xmm0 = x | y | z | 0.0f */
			"shufps xmm0,xmm0, 0x00\n"	//repete o valor da scalar para os demais bytes de xmm0
			"mulps xmm0,xmm1\n"
		    "movlpd	[ecx+ 0], xmm0\n"  //retira os tres bytes
			"shufps xmm0,xmm0, 0x02\n" //shuffle colocar nos bits 0-31 o 3o byte
			"movss	[ecx+ 8], xmm0\n"
	:
	: "r" (pOut), "r" (pIn),  "m" (test), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);


}

float  VPCALL CSIMD_SSE::vector3D_AlignedDot(const CVec3D* pSrc1, const CVec3D* pSrc2)
{

		/* SSE/SSE2 Implementation */
		float dummy;
		float *min = &dummy;
		asm(
			"mov eax, %0\n"                       // Load pointers into CPU regs and multiply
			"mov ebx, %1\n"
			"movaps xmm0, [eax]\n"
			"movaps xmm1, [ebx]\n"
			"andps xmm0, %2\n"  //re"move w component /* xmm0 = x | y | z | 0.0f */
			"andps xmm1, %2\n"  //re"move w component /* xmm1 = x | y | z | 0.0f */
			"mulps xmm1, xmm0\n"
			/* Shift data around in the registers (I loves me some haddps right now!!) */
			"movhlps xmm2, xmm1\n"    /* xmm2 = ?   | ?   | 0   |  z's */
		    "addss   xmm2, xmm1\n"    /* xmm2 = ?   | ?   | 0   | z's + x's */
			"shufps  xmm1, xmm1, 0x55\n"  /*0x55=01010101 / xmm1 = y's | y's | y's | y's */
			"addss   xmm2, xmm1\n" /* xmm2 = ?   | ?   | ?   | x's+y's+z's */
			//coloca 0-31 de xmm2 em dummy
			"mov		esi, %3\n"
			"movss		[esi], xmm2\n"

		:
		: "r" (pSrc1), "r" (pSrc2), "m" (_SIMDx86_float_SSE_NO_W_MASK),"m" (min)
		);


}


float  VPCALL CSIMD_SSE::vector3D_AlignedLengthSq(const CVec3D* pSrc1)
{


		float dummy;


		/* SSE/SSE2 Implementation */
		float *min = &dummy;
		asm(
			"mov eax, %0\n"                       // Load pointers into CPU regs and multiply
			"movaps xmm0, [eax]\n"
			"andps xmm0, %1\n"  //re"move w component /* xmm0 = x | y | z | 0.0f */
			"mulps xmm0, xmm0\n"    /* xmm0 = xmm0^2 */

			/* Shift data around in the registers (I loves me some haddps right now!!) */
			"movhlps xmm1, xmm0\n"		/* xmm1 = ?   | ?   | 0   |  z's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
			"shufps  xmm0, xmm0, 0x55\n"  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */

		    //coloca 0-31 de xmm1 em dummy
			"mov		esi, %2\n"
			"movss		[esi], xmm1\n"


		:
		: "r" (pSrc1), "m" (_SIMDx86_float_SSE_NO_W_MASK),"m" (min)
		);
		/* USE_SSE == 1 ||  USE_SSE == 2 */

		return dummy;

}

float  VPCALL CSIMD_SSE::vector3D_AlignedLength(const CVec3D* pSrc1)
{


		float dummy;
		/* SSE/SSE2 Implementation */
		float *min = &dummy;
		asm(
			"mov eax, %0\n"                       // Load pointers into CPU regs and multiply
			"movaps xmm0, [eax]\n"
			"andps xmm0, %1\n" //re"move w component /* xmm0 = x | y | z | 0.0f */
			"mulps xmm0, xmm0\n"    /* xmm0 = xmm0^2 */

				/* Shift data around in the registers (I loves me some haddps right now!!) */
			"movhlps xmm1, xmm0\n"		/* xmm1 = ?   | ?   | 0   |  z's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
			"shufps  xmm0, xmm0, 0x55\n"  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			#ifdef HIPREC
			/* Full square root */
			"sqrtss xmm1, xmm1\n"
			#else
			/* rcp( rSqrt(value) ) (This may be very inaccurate) */
			"rsqrtss xmm1, xmm1\n"
			"rcpss xmm1, xmm1\n"
			#endif

		    //coloca 0-31 de xmm1 em dummy
			"mov		esi, %2\n"
			"movss		[esi], xmm1\n"


		:
		: "r" (pSrc1), "m" (_SIMDx86_float_SSE_NO_W_MASK),"m" (min)
		);




		return dummy;


}

void VPCALL CSIMD_SSE::vector3D_AlignedCross(CVec3D* pLeft, const CVec3D* pRight)
{

	asm(
			"mov eax, %0\n"
			"mov ebx, %1\n" // Load pointers into CPU regs and multiply
			"mov ecx, %0\n"
			"movaps xmm0, [eax]\n"  /* xmm0 = pLeft */
			"movaps xmm1, [ebx]\n"  /* xmm1 = pRight */
			"andps xmm0, %2\n"  //re"move w component /* xmm0 = x | y | z | 0.0f */
			"andps xmm1, %2\n"  //re"move w component /* xmm1 = x | y | z | 0.0f */
			"movaps xmm2,xmm0\n"  /* xmm2 = pLeft */
			"movaps xmm3,xmm1\n"  /* xmm3 = pRight */
			/*
			Shuffles
			w | x | z | y  == 1100 1001 = 0xC9
			w | y | x | z  == 1101 0010 = 0xD2
			*/
			"shufps xmm0,xmm0, 0xC9\n"  /* 0xC9=11001001 */
			"shufps xmm1,xmm1, 0xD2\n"  /* 0xD2=11010010 */
			"shufps xmm2,xmm2, 0xD2\n"  /* 0xD2=11010010 */
			"shufps xmm3,xmm3, 0xC9\n"  /* 0xC9=11001001 */
			/* Multiply columns 1&2 and 3&4 */
			"mulps xmm1, xmm0\n"
			"mulps xmm3, xmm2\n"

			/* Subtract products to get the cross product! */
			"subps xmm1, xmm3\n"

			/* Store */
	        "movlpd	[ecx+ 0], xmm1\n"  //retira os tres bytes
			"shufps xmm1,xmm1, 0x02\n" //shuffle colocar nos bits 0-31 o 3o byte
			"movss	[ecx+ 8], xmm1\n"

	:
	: "r" (pLeft), "r" (pRight), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);


}

void VPCALL CSIMD_SSE::vector3D_AlignedCrossOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight)
{

		asm(
			"mov eax, %1\n"
			"mov ebx, %2\n" // Load pointers into CPU regs and multiply
			"mov ecx, %0\n"
			"movaps xmm0, [eax]\n"  /* xmm0 = pLeft */
			"movaps xmm1, [ebx]\n"  /* xmm1 = pRight */
			"andps xmm0, %3\n"  //re"move w component /* xmm0 = x | y | z | 0.0f */
			"andps xmm1, %3\n"  //re"move w component /* xmm1 = x | y | z | 0.0f */
			"movaps xmm2,xmm0\n"  /* xmm2 = pLeft */
			"movaps xmm3,xmm1\n"  /* xmm3 = pRight */
			/*
			Shuffles
			w | x | z | y  == 1100 1001 = 0xC9
			w | y | x | z  == 1101 0010 = 0xD2
			*/
			"shufps xmm0,xmm0, 0xC9\n"  /* 0xC9=11001001 */
			"shufps xmm1,xmm1, 0xD2\n"  /* 0xD2=11010010 */
			"shufps xmm2,xmm2, 0xD2\n"  /* 0xD2=11010010 */
			"shufps xmm3,xmm3, 0xC9\n"  /* 0xC9=11001001 */
			/* Multiply columns 1&2 and 3&4 */
			"mulps xmm1, xmm0\n"
			"mulps xmm3, xmm2\n"

			/* Subtract products to get the cross product! */
			"subps xmm1, xmm3\n"

			/* Store */
	        "movlpd	[ecx+ 0], xmm1\n"  //retira os tres bytes
			"shufps xmm1,xmm1, 0x02\n" //shuffle colocar nos bits 0-31 o 3o byte
			"movss	[ecx+ 8], xmm1\n"

	:
	: "r" (pOut), "r" (pLeft), "r" (pRight), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);


}
void VPCALL CSIMD_SSE::vector3D_AlignedNormalize(CVec3D* pVec)
{

		/* SSE/SSE2 Implementation */
		asm(
			"mov eax, %0\n"
			"mov ecx, %0\n"
			"movaps xmm0, [eax]\n"  /* xmm0 = pVec */
			"andps xmm0, %1\n"  //re"move w component /* xmm0 = x | y | z | 0.0f */
			"movaps xmm7, xmm0\n"			/* Save for the division by length */
		    "mulps xmm0, xmm0\n"    /* xmm0 = xmm0^2 */

				/* Shift data around in the registers (I loves me some haddps right now!!) */
			"movhlps xmm1, xmm0\n"		/* xmm1 = ?   | ?   | 0   |  z's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
			"shufps  xmm0, xmm0, 0x55\n"  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */


			#ifdef HIPREC

			/* Divide by magnitude */
			"sqrtss xmm1, xmm1\n"			/* xmm1 = ? | ? | ? | mag(pVec) */
			"shufps xmm1, xmm1,	0x00\n" /*$0x00=00000000,  xmm1 = mag(pVec) | mag(pVec) | mag(pVec) | mag(pVec) */
			"divps xmm7, xmm1\n"			/* xmm1 = ? | z/mag | y/mag | x/mag */

			#else

			/* Multiply by reciprocal magnitude */
			"rsqrtss xmm1 ,xmm1\n"			/* xmm1 = ? | ? | ? | rcp(mag) */
			"shufps xmm1, xmm1, $0x00\n"	        /*$0x00=00000000,  xmm3 = rcp(mag) | rcp(mag) | rcp(mag) | rcp(mag) */
			"mulps xmm7, xmm1\n"			/* xmm7 = ? | z/mag | y/mag | x/mag */

			#endif

			/* Store */
	        "movlpd	[ecx+ 0], xmm7\n"  //retira os tres bytes
			"shufps xmm7,xmm7, 0x02\n" //shuffle colocar nos bits 0-31 o 3o byte
			"movss	[ecx+ 8], xmm7\n"
			:
		: "r" (pVec), "m" (_SIMDx86_float_SSE_NO_W_MASK)
		);


}
void VPCALL CSIMD_SSE::vector3D_AlignedNormalizeOf(CVec3D* pOut, const CVec3D* pVec)
{

		/* SSE/SSE2 Implementation */
		asm(
			"mov eax, %0\n"
			"mov ecx, %1\n"
			"movups xmm0, [eax]\n"  /* xmm0 = pVec */
			"andps xmm0, %2\n"  //re"move w component /* xmm0 = x | y | z | 0.0f */
			"movaps xmm7, xmm0\n"			/* Save for the division by length */
		    "mulps xmm0, xmm0\n"    /* xmm0 = xmm0^2 */

				/* Shift data around in the registers (I loves me some haddps right now!!) */
			"movhlps xmm1, xmm0\n"		/* xmm1 = ?   | ?   | 0   |  z's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
			"shufps  xmm0, xmm0, 0x55\n"  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */


			#ifdef HIPREC

			/* Divide by magnitude */
			"sqrtss xmm1, xmm1\n"			/* xmm1 = ? | ? | ? | mag(pVec) */
			"shufps xmm1, xmm1,	0x00\n" /*$0x00=00000000,  xmm1 = mag(pVec) | mag(pVec) | mag(pVec) | mag(pVec) */
			"divps xmm7, xmm1\n"			/* xmm1 = ? | z/mag | y/mag | x/mag */

			#else

			/* Multiply by reciprocal magnitude */
			"rsqrtss xmm1 ,xmm1\n"			/* xmm1 = ? | ? | ? | rcp(mag) */
			"shufps xmm1, xmm1, $0x00\n"	        /*$0x00=00000000,  xmm3 = rcp(mag) | rcp(mag) | rcp(mag) | rcp(mag) */
			"mulps xmm7, xmm1\n"			/* xmm7 = ? | z/mag | y/mag | x/mag */

			#endif

			/* Store */
	        "movlpd	[ecx+ 0], xmm7\n"  //retira os tres bytes
			"shufps xmm7,xmm7, 0x02\n" //shuffle colocar nos bits 0-31 o 3o byte
			"movss	[ecx+ 8], xmm7\n"

			:
		: "r" (pVec), "r" (pOut), "m" (_SIMDx86_float_SSE_NO_W_MASK)
		);



}

float  VPCALL CSIMD_SSE::vector3D_AlignedDistance(const CVec3D* pVec1, const CVec3D* pVec2)
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



//======================  CVec4D============================

void VPCALL CSIMD_SSE::vector4D_Sum(CVec4D* pOut, const CVec4D* pIn)
{
	/* SSE/SSE2/SSE3 Implementation */
	 asm(
                "mov eax, %1\n"                      // Load pointers into CPU regs
                "mov ebx, %0\n"
				"mov ecx, %0\n"
                "movups xmm0, [eax]\n"                // Move unaligned vectors to SSE regs
                "movups xmm1, [ebx]\n"
				"addps xmm0, xmm1\n"                   // add vector elements
                "movups	[%0], xmm0\n"  //retira os tres bytes

	:
	: "r" (pOut), "r" (pIn)
	:"eax","ebx","ecx");

}

void VPCALL CSIMD_SSE::vector4D_SumOf(CVec4D* pOut, const CVec4D* pIn1, const CVec4D* pIn2)
{
	/* SSE/SSE2/SSE3 Implementation */
		 asm(
                "mov eax, %0\n"                       // Load pointers into CPU regs
                "mov ebx, %1\n"
				"mov ecx, %2\n"
                "movups xmm0, [eax]\n"                 // Move unaligned vectors to SSE regs
                "movups xmm1, [ebx]\n"
				"addps xmm0, xmm1\n"                   // add vector elements
                "movups	[ecx], xmm0\n"  //retira os tres bytes

	:
	: "r" (pIn1), "r" (pIn2), "r" (pOut)
	);
}


//Subtracts pIn from pOut and stores the result in pOut. POut = pIn-pOut
void VPCALL CSIMD_SSE::vector4D_Diff(CVec4D* pOut, CVec4D* pIn)
{
	/* SSE/SSE2/SSE3 Implementation */
		 asm(
                "mov eax, %0\n"                      // Load pointers into CPU regs
                "mov ebx, %1\n"
				"mov ecx, %0\n"
                "movups xmm0, [eax]\n"                 // Move unaligned vectors to SSE regs
                "movups xmm1, [ebx]\n"
				"subps xmm0, xmm1\n"                   // xmm0 = xmm0 - xmm1
                "movups	[ecx], xmm0\n"  //retira os tres bytes

	:
	: "r" (pOut), "r" (pIn)
	);
}


//Subtracts pIn2 from pIn1 and stores the result in pOut.  pout = pIn1-Pin2
void VPCALL CSIMD_SSE::vector4D_DiffOf(CVec4D* pOut, const CVec4D* pIn1, const CVec4D* pIn2)
{

	asm(
                "mov eax, %1\n"                       // Load pointers into CPU regs
                "mov ebx, %2\n"
				"mov ecx, %0\n"
                "movups xmm0, [eax]\n"                 // Move unaligned vectors to SSE regs
                "movups xmm1, [ebx]\n"
				"subps xmm0, xmm1\n"                   // xmm0 = xmm0 - xmm1
                "movups	[%0], xmm0\n"  //retira os tres bytes
	:
	: "r" (pOut), "r" (pIn1), "r" (pIn2)
	:"eax","ebx","ecx");
}
//Scales the components of pOut by scalar and stores the result in pOut.
void VPCALL CSIMD_SSE::vector4D_Scale(CVec4D* pOut, float scalar)
{
	float *test=&scalar;
	asm(
			"mov eax, %0\n"                       // Load pointers into CPU regs
			"mov ebx, %1\n"
			"mov ecx, %0\n"
            "movss xmm0, [ebx]\n"
			"movups xmm1, [eax]\n"
			"shufps xmm0,xmm0,0x00\n"	//repete o valor da scalar para os demais bytes de xmm0
			"mulps xmm0,xmm1\n"
		    "movups	[ecx+ 0], xmm0\n"  //retira os tres bytes\n"

	:
	: "r" (pOut), "m" (test)
	);


}
//Scales the components of pIn by scalar and stores the result in pOut.
void VPCALL CSIMD_SSE::vector4D_ScaleOf(CVec4D* pOut, const CVec4D* pIn, float scalar)
{
	float *test=&scalar;
	asm(
			"mov eax, %0\n"                       // Load pointers into CPU regs
			"mov ebx, %1\n"
			"mov ecx, %2\n"
            "movss xmm0, [ebx]\n"
			"movups xmm1, [eax]\n"
			"shufps xmm0,xmm0,0x00\n"	//repete o valor da scalar para os demais bytes de xmm0
			"mulps xmm0,xmm1\n"
		    "movups	[ecx+ 0], xmm0\n"  //retira os tres bytes\n"

	:
	: "r" (pIn), "m" (test),"r"(pOut)
	);


}
//Calculates the dot product (known as the inner product) of pSrc1 and pSrc2, that represents the cosine of the angle between them, if both are normalized.

float  VPCALL CSIMD_SSE::vector4D_Dot(const CVec4D* pSrc1, const CVec4D* pSrc2)
{

	/* SSE/SSE2 Implementation */
	float dummy;
	float *min = &dummy;
		asm(
			"mov eax, %0\n"                       // Load pointers into CPU regs and multiply
			"mov ebx, %1\n"
			"movups xmm0, [eax]\n"
			"movups xmm1, [ebx]\n"
			"mulps xmm1, xmm0\n"
			"movaps xmm2,xmm1\n"
			"shufps xmm2, xmm2,  0x1B\n" /*0x1B= 00011011 / xmm2 = x | y | z | w */

			"addps  xmm2, xmm1\n"			/* xmm2 = w+x | y+z | y+z | w+x */
			"movss  xmm3, xmm2\n"			/* xmm3 = ??? | ??? | ??? | w+x */
			"shufps xmm1, xmm1, 0x01\n"  /* 0x01= 00000001 / xmm3 = ??? | ??? | ??? | w+x */
			"addss  xmm2, xmm3\n"			/* xmm2 = ??? | ??? | ??? | dot4 */


			//coloca 0-31 de xmm2 em dummy
			"mov			esi, %2\n"
			"movss		[esi], xmm2\n"

	:
	: "r" (pSrc1), "r" (pSrc2),"m" (min)
	);

	return dummy;



}
//Returns the squared length of pSrc1. This can be useful when doing sphere-sphere collision tests where the exact distance is not needed and a squared one can be used, for example.
//It is signifigantly faster than using the square root version SIMDx86Vector_Length().

float  VPCALL CSIMD_SSE::vector4D_LengthSq(const CVec4D* pSrc1)
{

		float dummy;
		/* SSE/SSE2 Implementation */
		float *min = &dummy;
		asm(
			"mov eax, %0\n"                       // Load pointers into CPU regs and multiply
			"movups xmm0, [eax]\n"  /* xmm0 = x | y | z | w */
			"mulps xmm0, xmm0\n"    /* xmm0 = xmm0^2 */

			/* Shift data around in the registers (I loves me some haddps right now!!) */
			"movaps  xmm7,xmm0\n"
			"shufps  xmm7,xmm7, 0x03\n"  /* xmm7 = ?   | ?   | ?  |  w's *///colocar nos bis 0-31 de xmm7 os bits 96-127 de xmm7

			"movhlps xmm1, xmm0\n"		/* xmm1 = ?   | ?   | w's  |  z's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | w's+y's   | z's + x's */
			"shufps  xmm0, xmm0, 0x55\n"  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			"addss xmm1, xmm7\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's+w's */

		    //coloca 0-31 de xmm1 em dummy
			"mov		esi, %1\n"
			"movss		[esi], xmm1\n"


		:
		: "r" (pSrc1),"m" (min)
		);

		/* USE_SSE == 1 ||  USE_SSE == 2 */

		return dummy;

}
//Returns the length (magnitude) of pSrc1. Normalized vectors (e.g. the output of SIMDx86Vector_Normalize() ) are of unit length, or 1.0)
float  VPCALL CSIMD_SSE::vector4D_Length(const CVec4D* pSrc1)
{

		float dummy;
		/* SSE/SSE2 Implementation */
		float *min = &dummy;
		asm(
			"mov eax, %0\n"                       // Load pointers into CPU regs and multiply
			"movups xmm0, [eax]\n"
			"mulps xmm0, xmm0\n"    /* xmm0 = xmm0^2 */

			/* Shift data around in the registers (I loves me some haddps right now!!) */
			"movaps  xmm7,xmm0\n"
			"shufps  xmm7,xmm7, 0x03\n"  /* xmm7 = ?   | ?   | ?  |  w's *///colocar nos bis 0-31 de xmm7 os bits 96-127 de xmm7

			"movhlps xmm1, xmm0\n"		/* xmm1 = ?   | ?   | w's  |  z's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | w's+y's   | z's + x's */
			"shufps  xmm0, xmm0, 0x55\n"  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			"addss xmm1, xmm7\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's+w's */
			#ifdef HIPREC
			/* Full square root */
			"sqrtss xmm1, xmm1\n"
			#else
			/* rcp( rSqrt(value) ) (This may be very inaccurate) */
			"rsqrtss xmm1, xmm1\n"
			"rcpss xmm1, xmm1\n"
			#endif

		    //coloca 0-31 de xmm1 em dummy
			"mov			esi, %1\n"
			"movss		[esi], xmm1\n"


		:
		: "r" (pSrc1), "m" (min)
		);




		return dummy;


}




void VPCALL CSIMD_SSE::vector4D_Normalize(CVec4D* pVec)
{


		/* SSE/SSE2 Implementation */
		asm(
			"mov eax, %0\n"
			"mov ecx, %0\n"
			"movups xmm0, [eax]\n"  /* xmm0 = pVec */
			"movaps xmm7, xmm0\n"			/* Save for the division by length */
		    "mulps xmm0, xmm0\n"    /* xmm0 = xmm0^2 */

			/* Shift data around in the registers (I loves me some haddps right now!!) */
			"movaps  xmm2,xmm0\n"
			"shufps  xmm2,xmm2, 0x03\n"  /* xmm2 = ?   | ?   | ?  |  w's *///colocar nos bis 0-31 de xmm7 os bits 96-127 de xmm7

			"movhlps xmm1, xmm0\n"		/* xmm1 = ?   | ?   | w's  |  z's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | w's+y's   | z's + x's */
			"shufps  xmm0, xmm0, 0x55\n"  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			"addss xmm1, xmm2\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's+w's */


			#ifdef HIPREC

			/* Divide by magnitude */
			"sqrtss xmm1, xmm1\n"			/* xmm1 = ? | ? | ? | mag(pVec) */
			"shufps xmm1, xmm1,	0x00\n" /*$0x00=00000000,  xmm1 = mag(pVec) | mag(pVec) | mag(pVec) | mag(pVec) */
			"divps xmm7, xmm1\n"			/* xmm1 = w/mag | z/mag | y/mag | x/mag */

			#else

			/* Multiply by reciprocal magnitude */
			"rsqrtss xmm1 ,xmm1\n"			/* xmm1 = ? | ? | ? | rcp(mag) */
			"shufps xmm1, xmm1, 0x00\n"	        /*$0x00=00000000,  xmm3 = rcp(mag) | rcp(mag) | rcp(mag) | rcp(mag) */
			"mulps xmm7, xmm1\n"			/* xmm7 = w/mag | z/mag | y/mag | x/mag */

			#endif

			/* Store */
	        "movups	[ecx+ 0], xmm7\n"  //retira os tres bytes

		:
		: "r" (pVec)
		);


}

void VPCALL CSIMD_SSE::vector4D_NormalizeOf(CVec4D* pOut, const CVec4D* pVec)
{

		/* SSE/SSE2 Implementation */
		asm(
			"mov eax, %0\n"
			"mov ecx, %1\n"
			"movups xmm0, [eax]\n"  /* xmm0 = pVec */
			"movaps xmm7, xmm0\n"			/* Save for the division by length */
		    "mulps xmm0, xmm0\n"    /* xmm0 = xmm0^2 */

			/* Shift data around in the registers (I loves me some haddps right now!!) */
			"movaps  xmm2,xmm0\n"
			"shufps  xmm2,xmm2, 0x03\n"  /* xmm7 = ?   | ?   | ?  |  w's *///colocar nos bis 0-31 de xmm7 os bits 96-127 de xmm7

			"movhlps xmm1, xmm0\n"		/* xmm1 = ?   | ?   | w's  |  z's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | w's+y's   | z's + x's */
			"shufps  xmm0, xmm0, 0x55\n"  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			"addss xmm1, xmm2\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's+w's */


			#ifdef HIPREC

			/* Divide by magnitude */
			"sqrtss xmm1, xmm1\n"			/* xmm1 = ? | ? | ? | mag(pVec) */
			"shufps xmm1, xmm1,	0x00\n" /*$0x00=00000000,  xmm1 = mag(pVec) | mag(pVec) | mag(pVec) | mag(pVec) */
			"divps xmm7, xmm1\n"			/* xmm1 = w/mag | z/mag | y/mag | x/mag */

			#else

			/* Multiply by reciprocal magnitude */
			"rsqrtss xmm1 ,xmm1\n"			/* xmm1 = ? | ? | ? | rcp(mag) */
			"shufps xmm1, xmm1, 0x00\n"	        /*$0x00=00000000,  xmm3 = rcp(mag) | rcp(mag) | rcp(mag) | rcp(mag) */
			"mulps xmm7, xmm1\n"			/* xmm7 = w/mag | z/mag | y/mag | x/mag */

			#endif

			/* Store */
	        "movups	[ecx+ 0], xmm7\n"  //retira os tres bytes

		:
		: "r" (pVec), "r" (pOut)
		);

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
	 asm(
                "mov eax, %1\n"                      // Load pointers into CPU regs
                "mov ebx, %0\n"
				"mov ecx, %0\n"
                "movaps xmm0, [eax]\n"                // Move unaligned vectors to SSE regs
                "movaps xmm1, [ebx]\n"
				"addps xmm0, xmm1\n"                   // add vector elements
                "movaps	[%0], xmm0\n"  //retira os tres bytes

	:
	: "r" (pOut), "r" (pIn)
	:"eax","ebx","ecx");

}

void VPCALL CSIMD_SSE::vector4D_AlignedSumOf(CVec4D* pOut, const CVec4D* pIn1, const CVec4D* pIn2)
{
	/* SSE/SSE2/SSE3 Implementation */
		 asm(
                "mov eax, %0\n"                       // Load pointers into CPU regs
                "mov ebx, %1\n"
				"mov ecx, %2\n"
                "movaps xmm0, [eax]\n"                 // Move unaligned vectors to SSE regs
                "movaps xmm1, [ebx]\n"
				"addps xmm0, xmm1\n"                   // add vector elements
                "movaps	[ecx], xmm0\n"  //retira os tres bytes

	:
	: "r" (pIn1), "r" (pIn2), "r" (pOut)
	);
}


//Subtracts pIn from pOut and stores the result in pOut. POut = pIn-pOut
void VPCALL CSIMD_SSE::vector4D_AlignedDiff(CVec4D* pOut, CVec4D* pIn)
{
	/* SSE/SSE2/SSE3 Implementation */
		 asm(
                "mov eax, %0\n"                      // Load pointers into CPU regs
                "mov ebx, %1\n"
				"mov ecx, %0\n"
                "movaps xmm0, [eax]\n"                 // Move unaligned vectors to SSE regs
                "movaps xmm1, [ebx]\n"
				"subps xmm0, xmm1\n"                   // xmm0 = xmm0 - xmm1
                "movaps	[ecx], xmm0\n"  //retira os tres bytes

	:
	: "r" (pOut), "r" (pIn)
	);
}


//Subtracts pIn2 from pIn1 and stores the result in pOut.  pout = pIn1-Pin2
void VPCALL CSIMD_SSE::vector4D_AlignedDiffOf(CVec4D* pOut, const CVec4D* pIn1, const CVec4D* pIn2)
{

	asm(
                "mov eax, %1\n"                       // Load pointers into CPU regs
                "mov ebx, %2\n"
				"mov ecx, %0\n"
                "movaps xmm0, [eax]\n"                 // Move unaligned vectors to SSE regs
                "movaps xmm1, [ebx]\n"
				"subps xmm0, xmm1\n"                   // xmm0 = xmm0 - xmm1
                "movaps	[%0], xmm0\n"  //retira os tres bytes
	:
	: "r" (pOut), "r" (pIn1), "r" (pIn2)
	:"eax","ebx","ecx");
}
//Scales the components of pOut by scalar and stores the result in pOut.
void VPCALL CSIMD_SSE::vector4D_AlignedScale(CVec4D* pOut, float scalar)
{
	float *test=&scalar;
	asm(
			"mov eax, %0\n"                       // Load pointers into CPU regs
			"mov ebx, %1\n"
			"mov ecx, %0\n"
            "movss xmm0, [ebx]\n"
			"movaps xmm1, [eax]\n"
			"shufps xmm0,xmm0,0x00\n"	//repete o valor da scalar para os demais bytes de xmm0
			"mulps xmm0,xmm1\n"
		    "movaps	[ecx+ 0], xmm0\n"  //retira os tres bytes\n"

	:
	: "r" (pOut), "m" (test)
	);


}
//Scales the components of pIn by scalar and stores the result in pOut.
void VPCALL CSIMD_SSE::vector4D_AlignedScaleOf(CVec4D* pOut, const CVec4D* pIn, float scalar)
{
	float *test=&scalar;
	asm(
			"mov eax, %0\n"                       // Load pointers into CPU regs
			"mov ebx, %1\n"
			"mov ecx, %2\n"
            "movss xmm0, [ebx]\n"
			"movaps xmm1, [eax]\n"
			"shufps xmm0,xmm0,0x00\n"	//repete o valor da scalar para os demais bytes de xmm0
			"mulps xmm0,xmm1\n"
		    "movaps	[ecx+ 0], xmm0\n"  //retira os tres bytes\n"

	:
	: "r" (pIn), "m" (test),"r"(pOut)
	);


}
//Calculates the dot product (known as the inner product) of pSrc1 and pSrc2, that represents the cosine of the angle between them, if both are normalized.

float  VPCALL CSIMD_SSE::vector4D_AlignedDot(const CVec4D* pSrc1, const CVec4D* pSrc2)
{

	/* SSE/SSE2 Implementation */
	float dummy;
	float *min = &dummy;
		asm(
			"mov eax, %0\n"                       // Load pointers into CPU regs and multiply
			"mov ebx, %1\n"
			"movaps xmm0, [eax]\n"
			"movaps xmm1, [ebx]\n"
			"mulps xmm1, xmm0\n"
			"movaps xmm2,xmm1\n"
			"shufps xmm2, xmm2,  0x1B\n" /*0x1B= 00011011 / xmm2 = x | y | z | w */

			"addps  xmm2, xmm1\n"			/* xmm2 = w+x | y+z | y+z | w+x */
			"movss  xmm3, xmm2\n"			/* xmm3 = ??? | ??? | ??? | w+x */
			"shufps xmm1, xmm1, 0x01\n"  /* 0x01= 00000001 / xmm3 = ??? | ??? | ??? | w+x */
			"addss  xmm2, xmm3\n"			/* xmm2 = ??? | ??? | ??? | dot4 */


			//coloca 0-31 de xmm2 em dummy
			"mov			esi, %2\n"
			"movss		[esi], xmm2\n"

	:
	: "r" (pSrc1), "r" (pSrc2),"m" (min)
	);

	return dummy;



}
//Returns the squared length of pSrc1. This can be useful when doing sphere-sphere collision tests where the exact distance is not needed and a squared one can be used, for example.
//It is signifigantly faster than using the square root version SIMDx86Vector_Length().

float  VPCALL CSIMD_SSE::vector4D_AlignedLengthSq(const CVec4D* pSrc1)
{

		float dummy;
		/* SSE/SSE2 Implementation */
		float *min = &dummy;
		asm(
			"mov eax, %0\n"                       // Load pointers into CPU regs and multiply
			"movaps xmm0, [eax]\n"  /* xmm0 = x | y | z | w */
			"mulps xmm0, xmm0\n"    /* xmm0 = xmm0^2 */

			/* Shift data around in the registers (I loves me some haddps right now!!) */
			"movaps  xmm7,xmm0\n"
			"shufps  xmm7,xmm7, 0x03\n"  /* xmm7 = ?   | ?   | ?  |  w's *///colocar nos bis 0-31 de xmm7 os bits 96-127 de xmm7

			"movhlps xmm1, xmm0\n"		/* xmm1 = ?   | ?   | w's  |  z's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | w's+y's   | z's + x's */
			"shufps  xmm0, xmm0, 0x55\n"  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			"addss xmm1, xmm7\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's+w's */

		    //coloca 0-31 de xmm1 em dummy
			"mov		esi, %1\n"
			"movss		[esi], xmm1\n"


		:
		: "r" (pSrc1),"m" (min)
		);

		/* USE_SSE == 1 ||  USE_SSE == 2 */

		return dummy;

}
//Returns the length (magnitude) of pSrc1. Normalized vectors (e.g. the output of SIMDx86Vector_Normalize() ) are of unit length, or 1.0)
float  VPCALL CSIMD_SSE::vector4D_AlignedLength(const CVec4D* pSrc1)
{

		float dummy;
		/* SSE/SSE2 Implementation */
		float *min = &dummy;
		asm(
			"mov eax, %0\n"                       // Load pointers into CPU regs and multiply
			"movaps xmm0, [eax]\n"
			"mulps xmm0, xmm0\n"    /* xmm0 = xmm0^2 */

			/* Shift data around in the registers (I loves me some haddps right now!!) */
			"movaps  xmm7,xmm0\n"
			"shufps  xmm7,xmm7, 0x03\n"  /* xmm7 = ?   | ?   | ?  |  w's *///colocar nos bis 0-31 de xmm7 os bits 96-127 de xmm7

			"movhlps xmm1, xmm0\n"		/* xmm1 = ?   | ?   | w's  |  z's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | w's+y's   | z's + x's */
			"shufps  xmm0, xmm0, 0x55\n"  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			"addss xmm1, xmm7\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's+w's */
			#ifdef HIPREC
			/* Full square root */
			"sqrtss xmm1, xmm1\n"
			#else
			/* rcp( rSqrt(value) ) (This may be very inaccurate) */
			"rsqrtss xmm1, xmm1\n"
			"rcpss xmm1, xmm1\n"
			#endif

		    //coloca 0-31 de xmm1 em dummy
			"mov			esi, %1\n"
			"movss		[esi], xmm1\n"


		:
		: "r" (pSrc1), "m" (min)
		);




		return dummy;


}




void VPCALL CSIMD_SSE::vector4D_AlignedNormalize(CVec4D* pVec)
{


		/* SSE/SSE2 Implementation */
		asm(
			"mov eax, %0\n"
			"mov ecx, %0\n"
			"movaps xmm0, [eax]\n"  /* xmm0 = pVec */
			"movaps xmm7, xmm0\n"			/* Save for the division by length */
		    "mulps xmm0, xmm0\n"    /* xmm0 = xmm0^2 */

			/* Shift data around in the registers (I loves me some haddps right now!!) */
			"movaps  xmm2,xmm0\n"
			"shufps  xmm2,xmm2, 0x03\n"  /* xmm2 = ?   | ?   | ?  |  w's *///colocar nos bis 0-31 de xmm7 os bits 96-127 de xmm7

			"movhlps xmm1, xmm0\n"		/* xmm1 = ?   | ?   | w's  |  z's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | w's+y's   | z's + x's */
			"shufps  xmm0, xmm0, 0x55\n"  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			"addss xmm1, xmm2\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's+w's */


			#ifdef HIPREC

			/* Divide by magnitude */
			"sqrtss xmm1, xmm1\n"			/* xmm1 = ? | ? | ? | mag(pVec) */
			"shufps xmm1, xmm1,	0x00\n" /*$0x00=00000000,  xmm1 = mag(pVec) | mag(pVec) | mag(pVec) | mag(pVec) */
			"divps xmm7, xmm1\n"			/* xmm1 = w/mag | z/mag | y/mag | x/mag */

			#else

			/* Multiply by reciprocal magnitude */
			"rsqrtss xmm1 ,xmm1\n"			/* xmm1 = ? | ? | ? | rcp(mag) */
			"shufps xmm1, xmm1, 0x00\n"	        /*$0x00=00000000,  xmm3 = rcp(mag) | rcp(mag) | rcp(mag) | rcp(mag) */
			"mulps xmm7, xmm1\n"			/* xmm7 = w/mag | z/mag | y/mag | x/mag */

			#endif

			/* Store */
	        "movaps	[ecx+ 0], xmm7\n"  //retira os tres bytes

		:
		: "r" (pVec)
		);


}

void VPCALL CSIMD_SSE::vector4D_AlignedNormalizeOf(CVec4D* pOut, const CVec4D* pVec)
{

		/* SSE/SSE2 Implementation */
		asm(
			"mov eax, %0\n"
			"mov ecx, %1\n"
			"movaps xmm0, [eax]\n"  /* xmm0 = pVec */
			"movaps xmm7, xmm0\n"			/* Save for the division by length */
		    "mulps xmm0, xmm0\n"    /* xmm0 = xmm0^2 */

			/* Shift data around in the registers (I loves me some haddps right now!!) */
			"movaps  xmm2,xmm0\n"
			"shufps  xmm2,xmm2, 0x03\n"  /* xmm7 = ?   | ?   | ?  |  w's *///colocar nos bis 0-31 de xmm7 os bits 96-127 de xmm7

			"movhlps xmm1, xmm0\n"		/* xmm1 = ?   | ?   | w's  |  z's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | w's+y's   | z's + x's */
			"shufps  xmm0, xmm0, 0x55\n"  /*0x55=01010101 / xmm0 = y's | y's | y's | y's */
			"addss xmm1, xmm0\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
			"addss xmm1, xmm2\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's+w's */


			#ifdef HIPREC

			/* Divide by magnitude */
			"sqrtss xmm1, xmm1\n"			/* xmm1 = ? | ? | ? | mag(pVec) */
			"shufps xmm1, xmm1,	0x00\n" /*$0x00=00000000,  xmm1 = mag(pVec) | mag(pVec) | mag(pVec) | mag(pVec) */
			"divps xmm7, xmm1\n"			/* xmm1 = w/mag | z/mag | y/mag | x/mag */

			#else

			/* Multiply by reciprocal magnitude */
			"rsqrtss xmm1 ,xmm1\n"			/* xmm1 = ? | ? | ? | rcp(mag) */
			"shufps xmm1, xmm1, 0x00\n"	        /*$0x00=00000000,  xmm3 = rcp(mag) | rcp(mag) | rcp(mag) | rcp(mag) */
			"mulps xmm7, xmm1\n"			/* xmm7 = w/mag | z/mag | y/mag | x/mag */

			#endif

			/* Store */
	        "movaps	[ecx+ 0], xmm7\n"  //retira os tres bytes

		:
		: "r" (pVec), "r" (pOut)
		);

	/* USE_SSE == 1 || USE_SSE == 2 */


}

float  VPCALL CSIMD_SSE::vector4D_AlignedDistance(const CVec4D* pVec1, const CVec4D* pVec2)
{
	/* TODO: Optimize me completely */
	CVec4D diff;
	vector4D_DiffOf(&diff, pVec1, pVec2);
	return vector4D_Length(&diff);
}

//========CMat4D==============================
#define LOADREGS0( pOut ) \
	 "mov ecx, "STRINGIZE(pOut)"\n"

#define LOADREGS( pOut , pIn ) \
	"mov ecx, "STRINGIZE(pOut)"\n"    \
	"mov eax, "STRINGIZE(pIn)"\n"

#define LOADREGS2( pOut , pIn1, pIn2 ) \
	"mov ecx, "STRINGIZE(pOut)"\n"    \
	"mov eax, "STRINGIZE(pIn1)"\n"		\
	"mov ebx, "STRINGIZE(pIn2)"\n"

void VPCALL CSIMD_SSE::mat4D_Sum(CMat4D* Out, const CMat4D* In)
{

		asm volatile (
            LOADREGS( %0, %1 )
            "movups  xmm0,[ecx]\n"
			"movups  xmm1,[ecx+16]\n"
			"movups  xmm2,[ecx+32]\n"
			"movups  xmm3,[ecx+48]\n"
			"movups  xmm4,[eax]\n"
			"movups  xmm5,[eax+16]\n"
			"movups  xmm6,[eax+32]\n"
			"movups  xmm7,[eax+48]\n"
			"addps	xmm4,xmm0\n"
			"addps	xmm5, xmm1\n"
			"addps	xmm6,xmm2\n"
			"addps	xmm7,xmm3\n"
			"movups   [ecx],xmm4\n"
			"movups  [ecx+16],xmm5\n"
			"movups  [ecx+32],xmm6\n"
			"movups  [ecx+48],xmm7\n"
	:
	: "r" (Out), "r" (In)
	:"eax","ebx","ecx"
	);

}


void VPCALL CSIMD_SSE::mat4D_SumOf(CMat4D* Out, const CMat4D* In1, const CMat4D* In2)
{

	asm volatile (
	LOADREGS2( %2 , %0, %1 )
	"movups  xmm0,[eax]\n"
	"movups  xmm1,[eax+16]\n"
	"movups  xmm2,[eax+32]\n"
	"movups  xmm3,[eax+48]\n"
	"movups  xmm4,[ebx]\n"
	"movups  xmm5,[ebx+16]\n"
	"movups  xmm6,[ebx+32]\n"
	"movups  xmm7,[ebx+48]\n"
	"addps	xmm4,xmm0\n"
	"addps	xmm5,xmm1\n"
	"addps	xmm6,xmm2\n"
	"addps	xmm7,xmm3\n"
	"movups  [ecx],xmm4\n"
	"movups  [ecx+16],xmm5\n"
	"movups  [ecx+32],xmm6\n"
	"movups  [ecx+48],xmm7\n"
	:
	: "r" (In1), "r" (In2), "r" (Out)
	:"eax","ebx","ecx"
	);

}

void VPCALL CSIMD_SSE::mat4D_Diff(CMat4D* Out, const CMat4D* In)
{

	asm volatile (
	LOADREGS( %0 , %1 )
	"movups    xmm0,[ecx]\n"
	"movups  xmm1,[ecx+16]\n"
	"movups  xmm2,[ecx+32]\n"
	"movups  xmm3,[ecx+48]\n"
	"movups    xmm4,[eax]\n"
	"movups  xmm5,[eax+16]\n"
	"movups  xmm6,[eax+32]\n"
	"movups  xmm7,[eax+48]\n"
	"subps  xmm0,xmm4\n"
	"subps  xmm1,xmm5\n"
	"subps  xmm2,xmm6\n"
	"subps  xmm3,xmm7\n"
	"movups   [ecx], xmm0\n"
	"movups  [ecx+16],xmm1\n"
	"movups  [ecx+32],xmm2\n"
	"movups  [ecx+48],xmm3\n"
	:
	: "r" (Out), "r" (In)
	:"eax","ebx","ecx"
	);
}

void VPCALL CSIMD_SSE::mat4D_DiffOf(CMat4D* Out, const CMat4D* In1, const CMat4D* In2)
{
	asm volatile(
	LOADREGS2( %2, %0, %1 )
	"movups  xmm0,[eax]\n"
	"movups  xmm1,[eax+16]\n"
	"movups  xmm2,[eax+32]\n"
	"movups  xmm3,[eax+48]\n"
	"movups  xmm4, [ebx]\n"
	"movups  xmm5,[ebx+16]\n"
	"movups  xmm6,[ebx+32]\n"
	"movups	xmm7, [ebx+48]\n"
	"subps	xmm0,xmm4\n"
	"subps	xmm1,xmm5\n"
	"subps	xmm2,xmm6\n"
	"subps	xmm3,xmm7\n"
	"movups  [ecx],xmm0\n"
	"movups  [ecx+16],xmm1\n"
	"movups  [ecx+32],xmm2\n"
	"movups	[ecx+48],xmm3\n"
	:
	: "r" (In1), "r" (In2), "r" (Out)
	:"eax","ebx","ecx"
	);
}

void VPCALL CSIMD_SSE::mat4D_Scale(CMat4D* mtx, float scalar)
{
	asm volatile (
	LOADREGS( %0, %1 )
	/* Store scalar in xmm4.x */
	"movss  xmm4, [eax]\n"

	/* Get the matrix into registers */
	"movups  xmm0,[ecx]\n"
	"movups  xmm1,[ecx+16]\n"
	"movups  xmm2,[ecx+32]\n"
	"movups  xmm3,[ecx+48]\n"

	/* Broadcast element x to yzw to make a duplicated scalar register */
	"shufps  xmm4, xmm4, 0x00\n"

	/* Scale the matrix in parallel */
	"mulps	xmm0, xmm4\n"
	"mulps	xmm1, xmm4\n"
	"mulps	xmm2, xmm4\n"
	"mulps	xmm3, xmm4\n"

	/* Store results */
	"movups  [ecx],xmm0\n"
	"movups  [ecx+16],xmm1\n"
	"movups  [ecx+32],xmm2\n"
	"movups  [ecx+48],xmm3\n"
	:
	: "r" (mtx), "r" (&scalar)
	:"eax","ebx","ecx"
	);
}

void VPCALL CSIMD_SSE::mat4D_ScaleOf(CMat4D* pOut, const CMat4D* pIn, float scalar)
{

	//LOADREGS2( pOut , pIn, scalar )

	asm volatile (
	"mov ecx, %2\n"
	"mov eax, %0\n"
	"movss	xmm4, %1\n"
	"movups  xmm0, [eax]\n"
	"movups  xmm1,[eax+16]\n"
	"movups  xmm2,[eax+32]\n"
	"movups  xmm3,[eax+48]\n"
	"shufps  xmm4, xmm4,0x00\n"
	"mulps	xmm0,xmm4\n"
	"mulps	xmm1,xmm4\n"
	"mulps	xmm2,xmm4\n"
	"mulps	xmm3,xmm4\n"
	"movups  [ecx],xmm0\n"
	"movups  [ecx+16],xmm1\n"
	"movups  [ecx+32],xmm2\n"
	"movups  [ecx+48],xmm3\n"
	:
	: "r" (pIn), "m" (scalar), "r" (pOut)
	:"eax","ebx","ecx"
	);
}

void VPCALL CSIMD_SSE::mat4D_Multiply(CMat4D* pLeft, const CMat4D* pRight)
{

	asm volatile(
    "mov ecx, %1\n"
    "mov eax, %0\n"
    "movups    xmm0,[eax]\n"	/* xmm0 = pRight[0..3] */
	"movups  xmm1,[eax+16]\n"	/* xmm1 = pRight[5..7] */
	"movups  xmm2,[eax+32]\n"	/* xmm2 = pRight[8..11] */
	"movups  xmm3,[eax+48]\n"	/* xmm3 = pRight[12..15] */

	/* Processes 1/2 of the matrix at a time (2x4), unrolled loop */
	"movss    xmm4, [ecx]\n"
	"movss    xmm6,[ecx+4]\n"
	"movss  xmm5,[ecx+16]\n"
	"movss   xmm7,[ecx+20]\n"
	"shufps xmm4, xmm4, 0x00\n"
	"shufps  xmm5, xmm5,0x00\n"
	"shufps  xmm6, xmm6,0x00\n"
	"shufps  xmm7, xmm7,0x00\n"
	"mulps  xmm4,xmm0\n"
	"mulps  xmm5,xmm0\n"
	"mulps  xmm6,xmm1\n"
	"mulps  xmm7,xmm1\n"
	"addps  xmm5,xmm7\n"
	"addps  xmm4,xmm6\n"


	"movss   xmm6,[ecx+8]\n"
	"movss	xmm7,[ecx+24]\n"
	"shufps  xmm6, xmm6,0x00\n"
	"shufps xmm7, xmm7, 0x00\n"
	"mulps  xmm6,xmm2\n"
	"mulps xmm7,xmm2\n"
	"addps xmm4,xmm6\n"
	"addps xmm5,xmm7\n"

	"movss   xmm6,[ecx+12]\n"
	"movss   xmm7,[ecx+28]\n"
	"shufps  xmm6, xmm6,0x00\n"
	"shufps  xmm7, xmm7,0x00\n"
	"mulps  xmm6,xmm3\n"
	"mulps  xmm7,xmm3\n"
	"addps  xmm4,xmm6\n"
	"addps  xmm5,xmm7\n"

	"movups  [ecx],xmm4\n"
	"movups  [ecx+16],xmm5\n"

	/* second half of the matrix */
	"movss   xmm4,[ecx+32]\n"
	"movss   xmm6,[ecx+36]\n"
	"movss   xmm5,[ecx+48]\n"
	"movss   xmm7,[ecx+52]\n"
	"shufps  xmm4, xmm4,0x00\n"
	"shufps  xmm5, xmm5,0x00\n"
	"mulps  xmm4,xmm0\n"
	"mulps  xmm5,xmm0\n"

	"shufps xmm6, xmm6, 0x00\n"
	"shufps  xmm7, xmm7,0x00\n"
	"mulps  xmm6,xmm1\n"
	"mulps  xmm7,xmm1\n"
	"addps  xmm4,xmm6\n"
	"addps  xmm5,xmm7\n"


	"movss  xmm6,[ecx+40]\n"
	"movss  xmm7,[ecx+56]\n"
	"shufps  xmm6, xmm6,0x00\n"
	"shufps  xmm7, xmm7,0x00\n"
	"mulps  xmm6,xmm2\n"
	"mulps  xmm7,xmm2\n"
	"addps  xmm4,xmm6\n"
	"addps  xmm5,xmm7\n"

	"movss   xmm6,[ecx+44]\n"
	"movss   xmm7,[ecx+60]\n"
	"shufps  xmm6, xmm6,0x00\n"
	"shufps  xmm7, xmm7,0x00\n"
	"mulps  xmm6,xmm3\n"
	"mulps  xmm7,xmm3\n"
	"addps  xmm4,xmm6\n"
	"addps  xmm5,xmm7\n"

	"movups  [ecx+32],xmm4\n"
	"movups  [ecx+48],xmm5\n"

	:
	: "r" (pRight), "r" (pLeft)
	:"eax","ebx","ecx"
	);

}

void VPCALL CSIMD_SSE::mat4D_MultiplyOf(CMat4D* pOut, const CMat4D* pLeft, const CMat4D* pRight)
{

	asm volatile(
    "mov ecx, %1\n"
    "mov eax, %0\n"
    "mov ebx, %2\n"
    "movups    xmm0,[eax]\n"	/* xmm0 = pRight[0..3] */
	"movups  xmm1,[eax+16]\n"	/* xmm1 = pRight[5..7] */
	"movups  xmm2,[eax+32]\n"	/* xmm2 = pRight[8..11] */
	"movups  xmm3,[eax+48]\n"	/* xmm3 = pRight[12..15] */

	/* Processes 1/2 of the matrix at a time (2x4), unrolled loop */
	"movss    xmm4, [ecx]\n"
	"movss    xmm6,[ecx+4]\n"
	"movss  xmm5,[ecx+16]\n"
	"movss   xmm7,[ecx+20]\n"
	"shufps xmm4, xmm4, 0x00\n"
	"shufps  xmm5, xmm5,0x00\n"
	"shufps  xmm6, xmm6,0x00\n"
	"shufps  xmm7, xmm7,0x00\n"
	"mulps  xmm4,xmm0\n"
	"mulps  xmm5,xmm0\n"
	"mulps  xmm6,xmm1\n"
	"mulps  xmm7,xmm1\n"
	"addps  xmm5,xmm7\n"
	"addps  xmm4,xmm6\n"


	"movss   xmm6,[ecx+8]\n"
	"movss	xmm7,[ecx+24]\n"
	"shufps  xmm6, xmm6,0x00\n"
	"shufps xmm7, xmm7, 0x00\n"
	"mulps  xmm6,xmm2\n"
	"mulps xmm7,xmm2\n"
	"addps xmm4,xmm6\n"
	"addps xmm5,xmm7\n"

	"movss   xmm6,[ecx+12]\n"
	"movss   xmm7,[ecx+28]\n"
	"shufps  xmm6, xmm6,0x00\n"
	"shufps  xmm7, xmm7,0x00\n"
	"mulps  xmm6,xmm3\n"
	"mulps  xmm7,xmm3\n"
	"addps  xmm4,xmm6\n"
	"addps  xmm5,xmm7\n"

	"movups  [ebx],xmm4\n"
	"movups  [ebx+16],xmm5\n"

	/* second half of the matrix */
	"movss   xmm4,[ecx+32]\n"
	"movss   xmm6,[ecx+36]\n"
	"movss   xmm5,[ecx+48]\n"
	"movss   xmm7,[ecx+52]\n"
	"shufps  xmm4, xmm4,0x00\n"
	"shufps  xmm5, xmm5,0x00\n"
	"mulps  xmm4,xmm0\n"
	"mulps  xmm5,xmm0\n"

	"shufps xmm6, xmm6, 0x00\n"
	"shufps  xmm7, xmm7,0x00\n"
	"mulps  xmm6,xmm1\n"
	"mulps  xmm7,xmm1\n"
	"addps  xmm4,xmm6\n"
	"addps  xmm5,xmm7\n"


	"movss  xmm6,[ecx+40]\n"
	"movss  xmm7,[ecx+56]\n"
	"shufps  xmm6, xmm6,0x00\n"
	"shufps  xmm7, xmm7,0x00\n"
	"mulps  xmm6,xmm2\n"
	"mulps  xmm7,xmm2\n"
	"addps  xmm4,xmm6\n"
	"addps  xmm5,xmm7\n"

	"movss   xmm6,[ecx+44]\n"
	"movss   xmm7,[ecx+60]\n"
	"shufps  xmm6, xmm6,0x00\n"
	"shufps  xmm7, xmm7,0x00\n"
	"mulps  xmm6,xmm3\n"
	"mulps  xmm7,xmm3\n"
	"addps  xmm4,xmm6\n"
	"addps  xmm5,xmm7\n"

	"movups  [ebx+32],xmm4\n"
	"movups  [ebx+48],xmm5\n"

	:
	: "r" (pRight), "r" (pLeft), "r"(pOut)
	:"eax","ebx","ecx"
	);


}



void VPCALL CSIMD_SSE::mat4D_Transpose(CMat4D* pIn)
{

	asm volatile (
	LOADREGS0(%0 )
	"movlps	 xmm1,[ecx]\n"
		"movlps	 xmm3,[ecx+8]\n"
		"movhps	 xmm1,[ecx+16]\n"
		"movhps	 xmm3,[ecx+24]\n"
		"movlps	 xmm5,[ecx+32]\n"
		"movlps	 xmm4,[ecx+40]\n"
		"movhps	 xmm5,[ecx+48]\n"
		"movhps	 xmm4,[ecx+56]\n"
		"movaps	 xmm0,xmm1\n"
		"movaps	 xmm2,xmm3\n"
		"shufps	 xmm1,xmm5,0xDD\n"
		"shufps	 xmm3,xmm4,0xDD\n"
		"shufps	 xmm0,xmm5,0x88\n"
		"shufps	 xmm2,xmm4,0x88\n"
		"movups   [ecx],xmm0\n"
		"movups	 [ecx+16],xmm1\n"
		"movups   [ecx+32],xmm2\n"
		"movups   [ecx+48],xmm3\n"
		:
		: "r" (pIn)
		:"eax","ebx","ecx"
	);


}



void VPCALL CSIMD_SSE::mat4D_TransposeOf(CMat4D* pOut, const CMat4D* pIn)
{

	asm volatile(
	LOADREGS( %1 , %0 )
	"movlps	 xmm1,[eax]\n"
		"movlps	 xmm3,[eax+8]\n"
		"movhps	 xmm1,[eax+16]\n"
		"movhps	 xmm3,[eax+24]\n"
		"movlps	 xmm5,[eax+32]\n"
		"movlps	 xmm4,[eax+40]\n"
		"movhps	 xmm5,[eax+48]\n"
		"movhps	 xmm4,[eax+56]\n"
		"movaps	 xmm0,xmm1\n"
		"movaps	 xmm2,xmm3\n"
		"shufps	 xmm1,xmm5,0xDD\n"
		"shufps	 xmm3,xmm4,0xDD\n"
		"shufps	 xmm0,xmm5,0x88\n"
		"shufps	 xmm2,xmm4,0x88\n"
		"movups   [ecx],xmm0\n"
		"movups   [ecx+16],xmm1\n"
		"movups   [ecx+32],xmm2\n"
		"movups   [ecx+48],xmm3\n"
		:
		: "r" (pIn), "r" (pOut)
		:"eax","ebx","ecx"
	);

}

void VPCALL CSIMD_SSE::mat4D_VectorMultiply(CVec3D* pVec, const CMat4D* pMat)
{
	/*
		Does a normal vector x matrix, but since w = 1.0f, the final column of the matrix is
		merely added.
	*/

	asm volatile(
	LOADREGS( %1,%0 )
	"movups    xmm4,[ecx]\n"
	"movups    xmm0,[eax]\n"
	"movups  xmm1,[eax+16]\n"
	"movups  xmm2,[eax+32]\n"
	"movups  xmm3,[eax+48]\n"
	"movaps  xmm5,xmm4\n"
	"movaps  xmm6,xmm4\n"
	"shufps  xmm4, xmm4,0x00\n"    /* xmm4 = x | x | x | x */
	"shufps  xmm5, xmm5,0x55\n"   /* xmm5 = y | y | y | y */
	"shufps  xmm6, xmm6,0xAA\n"    /* xmm6 = z | z | z | z */

	/* Multiply with each row */
	"mulps  xmm0,xmm4\n"
	"mulps  xmm1,xmm5\n"
	"mulps  xmm2,xmm6\n"

	/* Sum results, xmm0-2 = transformed x,y,z colums, xmm3 = w column */
	"addps  xmm1,xmm0\n"    /* xmm1 = tx + ty */
	"addps  xmm3,xmm2\n"    /* xmm3 = tz + w */
	"addps  xmm1,xmm3\n"    /* xmm1 = tx + ty + tz + w*/

	/* Store translated vector */
	"movups  [ecx],xmm1\n"
	:
	: "r" (pMat), "r" (pVec)
	:"eax","ebx","ecx"
	);

}

void VPCALL CSIMD_SSE::mat4D_VectorMultiplyOf(CVec3D* pOut, const CVec3D* pIn, const CMat4D* pMat)
{
	/*
		Does a normal vector x matrix, but since w = 1.0f, the final column of the matrix is
		merely added.
	*/


	asm volatile(
    LOADREGS2( %2, %1, %0 )
    "movups    xmm4,[eax]\n"
	"movups    xmm0,[ebx]\n"
	"movups  xmm1,[ebx+16]\n"
	"movups  xmm2,[ebx+32]\n"
	"movups  xmm3,[ebx+48]\n"
	"movaps  xmm5,xmm4\n"
	"movaps  xmm6,xmm4\n"
	"shufps  xmm4, xmm4,0x00\n"    /* xmm4 = x | x | x | x */
	"shufps  xmm5, xmm5,0x55\n"    /* xmm5 = y | y | y | y */
	"shufps  xmm6, xmm6,0xAA\n"    /* xmm6 = z | z | z | z */

	/* Multiply with each row */
	"mulps  xmm0,xmm4\n"
	"mulps  xmm1,xmm5\n"
	"mulps  xmm2,xmm6\n"

	/* Sum results, xmm0-2 = transformed x,y,z colums, xmm3 = w column */
	"addps  xmm1,xmm0\n"    /* xmm1 = tx + ty */
	"addps  xmm3,xmm2\n"    /* xmm3 = tz + w */
	"addps  xmm1,xmm3\n"    /* xmm1 = tx + ty + tz + w*/

	/* Store translated vector */
	"movups  [ecx],xmm1\n"
	:
	: "r" (pMat), "r" (pIn), "r" (pOut)
	:"eax","ebx","ecx"
	);

}

void VPCALL CSIMD_SSE::mat4D_VectorMultiply(CVec4D* pOut4D, const CMat4D* pMat)
{


	asm volatile(
    LOADREGS( %1 , %0 )

	/*
	    IMPLICIT:
		movl 4(%esp), %1
		movl 8(%esp), %0

		Access %1 first, it has finished and has no dependencies, unlike loading from %0

		Surprisingly results in 9% speedup on Celeron, 6% on Athlon64
	*/

	"movups  xmm4,[ecx]\n"
	"movups  xmm0,[eax]\n"
	"movups  xmm1,[eax+16]\n"
	"movups  xmm2,[eax+32]\n"
	"movups  xmm3,[eax+48]\n"
	"movaps  xmm5,xmm4\n"
	"movaps  xmm6,xmm4\n"
	"movaps  xmm7,xmm4\n"
	"shufps  xmm4, xmm4,0x00\n"
	"shufps  xmm5, xmm5,0x55\n"
	"shufps  xmm6, xmm6,0xAA\n"
	"shufps  xmm7, xmm7,0xFF\n"

	/* Multiply with each row */
	"mulps  xmm0,xmm4\n"
	"mulps  xmm1,xmm5\n"
	"mulps  xmm2,xmm6\n"
	"mulps  xmm3,xmm7\n"

	/* Sum results */
	"addps  xmm1,xmm0\n"
	"addps  xmm3,xmm2\n"
	"addps  xmm1,xmm3\n"

	/* Store translated vector */
	"movups  [ecx],xmm1\n"
	:
	: "r" (pMat), "r" (pOut4D)
	:"eax","ebx","ecx"
	);

}


void VPCALL CSIMD_SSE::mat4D_VectorMultiplyOf(CVec4D* pOut4D, const CVec4D* pIn4D, const CMat4D* pMat)
{

	asm volatile(
	LOADREGS2(%2 , %1, %0 )
    /*
	    IMPLICIT:
		movl 4(%esp), %1
		movl 8(%esp), %0

		Access %1 first, it has finished and has no dependencies, unlike loading from %0

		Surprisingly results in 9% speedup on Celeron, 6% on Athlon64
	*/
	"movups  xmm4,[eax]\n"
	"movups    xmm0,[ebx]\n"
	"movups  xmm1,[ebx+16]\n"
	"movups  xmm2,[ebx+32]\n"
	"movups  xmm3,[ebx+48]\n"
	"movaps  xmm5,xmm4\n"
	"movaps  xmm6,xmm4\n"
	"movaps   xmm7,xmm4\n"
	"shufps  xmm4, xmm4,0x00\n"
	"shufps  xmm5, xmm5,0x55\n"
	"shufps  xmm6, xmm6,0xAA\n"
	"shufps  xmm7, xmm7,0xFF\n"
	"mulps  xmm0,xmm4\n"
	"mulps  xmm1,xmm5\n"
	"mulps  xmm2,xmm6\n"
	"mulps  xmm3,xmm7\n"
	"addps  xmm1,xmm0\n"
	"addps  xmm3,xmm2\n"
	"addps  xmm1,xmm3\n"
	"movups  [ecx],xmm1\n"
	:
	: "r" (pMat), "r" (pIn4D),"r" (pOut4D)
	:"eax","ebx","ecx"
	);

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

		asm volatile (
            LOADREGS( %0, %1 )
            "movups  xmm0,[ecx]\n"
			"movups  xmm1,[ecx+16]\n"
			"movups  xmm2,[ecx+32]\n"
			"movups  xmm3,[ecx+48]\n"
			"movups  xmm4,[eax]\n"
			"movups  xmm5,[eax+16]\n"
			"movups  xmm6,[eax+32]\n"
			"movups  xmm7,[eax+48]\n"
			"addps	xmm4,xmm0\n"
			"addps	xmm5, xmm1\n"
			"addps	xmm6,xmm2\n"
			"addps	xmm7,xmm3\n"
			"movups   [ecx],xmm4\n"
			"movups  [ecx+16],xmm5\n"
			"movups  [ecx+32],xmm6\n"
			"movups  [ecx+48],xmm7\n"
	:
	: "r" (Out), "r" (In)
	:"eax","ebx","ecx"
	);

}


void VPCALL CSIMD_SSE::mat4D_AlignedSumOf(CMat4D* Out, const CMat4D* In1, const CMat4D* In2)
{

	asm volatile (
	LOADREGS2( %2 , %0, %1 )
	"movaps  xmm0,[eax]\n"
	"movaps  xmm1,[eax+16]\n"
	"movaps  xmm2,[eax+32]\n"
	"movaps  xmm3,[eax+48]\n"
	"movaps  xmm4,[ebx]\n"
	"movaps  xmm5,[ebx+16]\n"
	"movaps  xmm6,[ebx+32]\n"
	"movaps  xmm7,[ebx+48]\n"
	"addps	xmm4,xmm0\n"
	"addps	xmm5,xmm1\n"
	"addps	xmm6,xmm2\n"
	"addps	xmm7,xmm3\n"
	"movaps  [ecx],xmm4\n"
	"movaps  [ecx+16],xmm5\n"
	"movaps  [ecx+32],xmm6\n"
	"movaps  [ecx+48],xmm7\n"
	:
	: "r" (In1), "r" (In2), "r" (Out)
	:"eax","ebx","ecx"
	);

}

void VPCALL CSIMD_SSE::mat4D_AlignedDiff(CMat4D* Out, const CMat4D* In)
{

	asm volatile (
	LOADREGS( %0 , %1 )
	"movaps    xmm0,[ecx]\n"
	"movaps  xmm1,[ecx+16]\n"
	"movaps  xmm2,[ecx+32]\n"
	"movaps  xmm3,[ecx+48]\n"
	"movaps    xmm4,[eax]\n"
	"movaps  xmm5,[eax+16]\n"
	"movaps  xmm6,[eax+32]\n"
	"movaps  xmm7,[eax+48]\n"
	"subps  xmm0,xmm4\n"
	"subps  xmm1,xmm5\n"
	"subps  xmm2,xmm6\n"
	"subps  xmm3,xmm7\n"
	"movaps   [ecx], xmm0\n"
	"movaps  [ecx+16],xmm1\n"
	"movaps  [ecx+32],xmm2\n"
	"movaps  [ecx+48],xmm3\n"
	:
	: "r" (Out), "r" (In)
	:"eax","ebx","ecx"
	);
}

void VPCALL CSIMD_SSE::mat4D_AlignedDiffOf(CMat4D* Out, const CMat4D* In1, const CMat4D* In2)
{
	asm volatile(
	LOADREGS2( %2, %0, %1 )
	"movaps  xmm0,[eax]\n"
	"movaps  xmm1,[eax+16]\n"
	"movaps  xmm2,[eax+32]\n"
	"movaps  xmm3,[eax+48]\n"
	"movaps  xmm4, [ebx]\n"
	"movaps  xmm5,[ebx+16]\n"
	"movaps  xmm6,[ebx+32]\n"
	"movaps	xmm7, [ebx+48]\n"
	"subps	xmm0,xmm4\n"
	"subps	xmm1,xmm5\n"
	"subps	xmm2,xmm6\n"
	"subps	xmm3,xmm7\n"
	"movaps  [ecx],xmm0\n"
	"movaps  [ecx+16],xmm1\n"
	"movaps  [ecx+32],xmm2\n"
	"movaps	[ecx+48],xmm3\n"
	:
	: "r" (In1), "r" (In2), "r" (Out)
	:"eax","ebx","ecx"
	);
}

void VPCALL CSIMD_SSE::mat4D_AlignedScale(CMat4D* mtx, float scalar)
{
	asm volatile (
	LOADREGS( %0, %1 )
	/* Store scalar in xmm4.x */
	"movss  xmm4, [eax]\n"

	/* Get the matrix into registers */
	"movaps  xmm0,[ecx]\n"
	"movaps  xmm1,[ecx+16]\n"
	"movaps  xmm2,[ecx+32]\n"
	"movaps  xmm3,[ecx+48]\n"

	/* Broadcast element x to yzw to make a duplicated scalar register */
	"shufps  xmm4, xmm4, 0x00\n"

	/* Scale the matrix in parallel */
	"mulps	xmm0, xmm4\n"
	"mulps	xmm1, xmm4\n"
	"mulps	xmm2, xmm4\n"
	"mulps	xmm3, xmm4\n"

	/* Store results */
	"movaps  [ecx],xmm0\n"
	"movaps  [ecx+16],xmm1\n"
	"movaps  [ecx+32],xmm2\n"
	"movaps  [ecx+48],xmm3\n"
	:
	: "r" (mtx), "r" (&scalar)
	:"eax","ebx","ecx"
	);
}

void VPCALL CSIMD_SSE::mat4D_AlignedScaleOf(CMat4D* pOut, const CMat4D* pIn, float scalar)
{

	//LOADREGS2( pOut , pIn, scalar )

	asm volatile (
	"mov ecx, %2\n"
	"mov eax, %0\n"
	"movss	xmm4, %1\n"
	"movaps  xmm0, [eax]\n"
	"movaps  xmm1,[eax+16]\n"
	"movaps  xmm2,[eax+32]\n"
	"movaps  xmm3,[eax+48]\n"
	"shufps  xmm4, xmm4,0x00\n"
	"mulps	xmm0,xmm4\n"
	"mulps	xmm1,xmm4\n"
	"mulps	xmm2,xmm4\n"
	"mulps	xmm3,xmm4\n"
	"movaps  [ecx],xmm0\n"
	"movaps  [ecx+16],xmm1\n"
	"movaps  [ecx+32],xmm2\n"
	"movaps  [ecx+48],xmm3\n"
	:
	: "r" (pIn), "m" (scalar), "r" (pOut)
	:"eax","ebx","ecx"
	);
}

void VPCALL CSIMD_SSE::mat4D_AlignedMultiply(CMat4D* pLeft, const CMat4D* pRight)
{

	asm volatile(
    "mov ecx, %1\n"
    "mov eax, %0\n"
    "movaps    xmm0,[eax]\n"	/* xmm0 = pRight[0..3] */
	"movaps  xmm1,[eax+16]\n"	/* xmm1 = pRight[5..7] */
	"movaps  xmm2,[eax+32]\n"	/* xmm2 = pRight[8..11] */
	"movaps  xmm3,[eax+48]\n"	/* xmm3 = pRight[12..15] */

	/* Processes 1/2 of the matrix at a time (2x4), unrolled loop */
	"movss    xmm4, [ecx]\n"
	"movss    xmm6,[ecx+4]\n"
	"movss  xmm5,[ecx+16]\n"
	"movss   xmm7,[ecx+20]\n"
	"shufps xmm4, xmm4, 0x00\n"
	"shufps  xmm5, xmm5,0x00\n"
	"shufps  xmm6, xmm6,0x00\n"
	"shufps  xmm7, xmm7,0x00\n"
	"mulps  xmm4,xmm0\n"
	"mulps  xmm5,xmm0\n"
	"mulps  xmm6,xmm1\n"
	"mulps  xmm7,xmm1\n"
	"addps  xmm5,xmm7\n"
	"addps  xmm4,xmm6\n"


	"movss   xmm6,[ecx+8]\n"
	"movss	xmm7,[ecx+24]\n"
	"shufps  xmm6, xmm6,0x00\n"
	"shufps xmm7, xmm7, 0x00\n"
	"mulps  xmm6,xmm2\n"
	"mulps xmm7,xmm2\n"
	"addps xmm4,xmm6\n"
	"addps xmm5,xmm7\n"

	"movss   xmm6,[ecx+12]\n"
	"movss   xmm7,[ecx+28]\n"
	"shufps  xmm6, xmm6,0x00\n"
	"shufps  xmm7, xmm7,0x00\n"
	"mulps  xmm6,xmm3\n"
	"mulps  xmm7,xmm3\n"
	"addps  xmm4,xmm6\n"
	"addps  xmm5,xmm7\n"

	"movaps  [ecx],xmm4\n"
	"movaps  [ecx+16],xmm5\n"

	/* second half of the matrix */
	"movss   xmm4,[ecx+32]\n"
	"movss   xmm6,[ecx+36]\n"
	"movss   xmm5,[ecx+48]\n"
	"movss   xmm7,[ecx+52]\n"
	"shufps  xmm4, xmm4,0x00\n"
	"shufps  xmm5, xmm5,0x00\n"
	"mulps  xmm4,xmm0\n"
	"mulps  xmm5,xmm0\n"

	"shufps xmm6, xmm6, 0x00\n"
	"shufps  xmm7, xmm7,0x00\n"
	"mulps  xmm6,xmm1\n"
	"mulps  xmm7,xmm1\n"
	"addps  xmm4,xmm6\n"
	"addps  xmm5,xmm7\n"


	"movss  xmm6,[ecx+40]\n"
	"movss  xmm7,[ecx+56]\n"
	"shufps  xmm6, xmm6,0x00\n"
	"shufps  xmm7, xmm7,0x00\n"
	"mulps  xmm6,xmm2\n"
	"mulps  xmm7,xmm2\n"
	"addps  xmm4,xmm6\n"
	"addps  xmm5,xmm7\n"

	"movss   xmm6,[ecx+44]\n"
	"movss   xmm7,[ecx+60]\n"
	"shufps  xmm6, xmm6,0x00\n"
	"shufps  xmm7, xmm7,0x00\n"
	"mulps  xmm6,xmm3\n"
	"mulps  xmm7,xmm3\n"
	"addps  xmm4,xmm6\n"
	"addps  xmm5,xmm7\n"

	"movaps  [ecx+32],xmm4\n"
	"movaps  [ecx+48],xmm5\n"

	:
	: "r" (pRight), "r" (pLeft)
	:"eax","ebx","ecx"
	);

}

void VPCALL CSIMD_SSE::mat4D_AlignedMultiplyOf(CMat4D* pOut, const CMat4D* pLeft, const CMat4D* pRight)
{


	asm volatile(
    "mov ecx, %1\n"
    "mov eax, %0\n"
    "mov ebx, %2\n"
    "movaps    xmm0,[eax]\n"	/* xmm0 = pRight[0..3] */
	"movaps  xmm1,[eax+16]\n"	/* xmm1 = pRight[5..7] */
	"movaps  xmm2,[eax+32]\n"	/* xmm2 = pRight[8..11] */
	"movaps  xmm3,[eax+48]\n"	/* xmm3 = pRight[12..15] */

	/* Processes 1/2 of the matrix at a time (2x4), unrolled loop */
	"movss    xmm4, [ecx]\n"
	"movss    xmm6,[ecx+4]\n"
	"movss  xmm5,[ecx+16]\n"
	"movss   xmm7,[ecx+20]\n"
	"shufps xmm4, xmm4, 0x00\n"
	"shufps  xmm5, xmm5,0x00\n"
	"shufps  xmm6, xmm6,0x00\n"
	"shufps  xmm7, xmm7,0x00\n"
	"mulps  xmm4,xmm0\n"
	"mulps  xmm5,xmm0\n"
	"mulps  xmm6,xmm1\n"
	"mulps  xmm7,xmm1\n"
	"addps  xmm5,xmm7\n"
	"addps  xmm4,xmm6\n"


	"movss   xmm6,[ecx+8]\n"
	"movss	xmm7,[ecx+24]\n"
	"shufps  xmm6, xmm6,0x00\n"
	"shufps xmm7, xmm7, 0x00\n"
	"mulps  xmm6,xmm2\n"
	"mulps xmm7,xmm2\n"
	"addps xmm4,xmm6\n"
	"addps xmm5,xmm7\n"

	"movss   xmm6,[ecx+12]\n"
	"movss   xmm7,[ecx+28]\n"
	"shufps  xmm6, xmm6,0x00\n"
	"shufps  xmm7, xmm7,0x00\n"
	"mulps  xmm6,xmm3\n"
	"mulps  xmm7,xmm3\n"
	"addps  xmm4,xmm6\n"
	"addps  xmm5,xmm7\n"

	"movaps  [ebx],xmm4\n"
	"movaps  [ebx+16],xmm5\n"

	/* second half of the matrix */
	"movss   xmm4,[ecx+32]\n"
	"movss   xmm6,[ecx+36]\n"
	"movss   xmm5,[ecx+48]\n"
	"movss   xmm7,[ecx+52]\n"
	"shufps  xmm4, xmm4,0x00\n"
	"shufps  xmm5, xmm5,0x00\n"
	"mulps  xmm4,xmm0\n"
	"mulps  xmm5,xmm0\n"

	"shufps xmm6, xmm6, 0x00\n"
	"shufps  xmm7, xmm7,0x00\n"
	"mulps  xmm6,xmm1\n"
	"mulps  xmm7,xmm1\n"
	"addps  xmm4,xmm6\n"
	"addps  xmm5,xmm7\n"


	"movss  xmm6,[ecx+40]\n"
	"movss  xmm7,[ecx+56]\n"
	"shufps  xmm6, xmm6,0x00\n"
	"shufps  xmm7, xmm7,0x00\n"
	"mulps  xmm6,xmm2\n"
	"mulps  xmm7,xmm2\n"
	"addps  xmm4,xmm6\n"
	"addps  xmm5,xmm7\n"

	"movss   xmm6,[ecx+44]\n"
	"movss   xmm7,[ecx+60]\n"
	"shufps  xmm6, xmm6,0x00\n"
	"shufps  xmm7, xmm7,0x00\n"
	"mulps  xmm6,xmm3\n"
	"mulps  xmm7,xmm3\n"
	"addps  xmm4,xmm6\n"
	"addps  xmm5,xmm7\n"

	"movaps  [ebx+32],xmm4\n"
	"movaps  [ebx+48],xmm5\n"

	:
	: "r" (pRight), "r" (pLeft), "r"(pOut)
	:"eax","ebx","ecx"
	);


}



void VPCALL CSIMD_SSE::mat4D_AlignedTranspose(CMat4D* pIn)
{

	asm volatile (
	LOADREGS0(%0 )
	"movlps	 xmm1,[ecx]\n"
		"movlps	 xmm3,[ecx+8]\n"
		"movhps	 xmm1,[ecx+16]\n"
		"movhps	 xmm3,[ecx+24]\n"
		"movlps	 xmm5,[ecx+32]\n"
		"movlps	 xmm4,[ecx+40]\n"
		"movhps	 xmm5,[ecx+48]\n"
		"movhps	 xmm4,[ecx+56]\n"
		"movaps	 xmm0,xmm1\n"
		"movaps	 xmm2,xmm3\n"
		"shufps	 xmm1,xmm5,0xDD\n"
		"shufps	 xmm3,xmm4,0xDD\n"
		"shufps	 xmm0,xmm5,0x88\n"
		"shufps	 xmm2,xmm4,0x88\n"
		"movaps   [ecx],xmm0\n"
		"movaps	 [ecx+16],xmm1\n"
		"movaps   [ecx+32],xmm2\n"
		"movaps   [ecx+48],xmm3\n"
		:
		: "r" (pIn)
		:"eax","ebx","ecx"
	);


}



void VPCALL CSIMD_SSE::mat4D_AlignedTransposeOf(CMat4D* pOut, const CMat4D* pIn)
{

	asm volatile(
	LOADREGS( %1 , %0 )
	"movlps	 xmm1,[eax]\n"
		"movlps	 xmm3,[eax+8]\n"
		"movhps	 xmm1,[eax+16]\n"
		"movhps	 xmm3,[eax+24]\n"
		"movlps	 xmm5,[eax+32]\n"
		"movlps	 xmm4,[eax+40]\n"
		"movhps	 xmm5,[eax+48]\n"
		"movhps	 xmm4,[eax+56]\n"
		"movaps	 xmm0,xmm1\n"
		"movaps	 xmm2,xmm3\n"
		"shufps	 xmm1,xmm5,0xDD\n"
		"shufps	 xmm3,xmm4,0xDD\n"
		"shufps	 xmm0,xmm5,0x88\n"
		"shufps	 xmm2,xmm4,0x88\n"
		"movaps   [ecx],xmm0\n"
		"movaps   [ecx+16],xmm1\n"
		"movaps   [ecx+32],xmm2\n"
		"movaps   [ecx+48],xmm3\n"
		:
		: "r" (pIn), "r" (pOut)
		:"eax","ebx","ecx"
	);

}

void VPCALL CSIMD_SSE::mat4D_AlignedVectorMultiply(CVec3D* pVec, const CMat4D* pMat)
{
	/*
		Does a normal vector x matrix, but since w = 1.0f, the final column of the matrix is
		merely added.
	*/

	asm volatile(
	LOADREGS( %1,%0 )
	"movaps    xmm4,[ecx]\n"
	"movaps    xmm0,[eax]\n"
	"movaps  xmm1,[eax+16]\n"
	"movaps  xmm2,[eax+32]\n"
	"movaps  xmm3,[eax+48]\n"
	"movaps  xmm5,xmm4\n"
	"movaps  xmm6,xmm4\n"
	"shufps  xmm4, xmm4,0x00\n"    /* xmm4 = x | x | x | x */
	"shufps  xmm5, xmm5,0x55\n"   /* xmm5 = y | y | y | y */
	"shufps  xmm6, xmm6,0xAA\n"    /* xmm6 = z | z | z | z */

	/* Multiply with each row */
	"mulps  xmm0,xmm4\n"
	"mulps  xmm1,xmm5\n"
	"mulps  xmm2,xmm6\n"

	/* Sum results, xmm0-2 = transformed x,y,z colums, xmm3 = w column */
	"addps  xmm1,xmm0\n"    /* xmm1 = tx + ty */
	"addps  xmm3,xmm2\n"    /* xmm3 = tz + w */
	"addps  xmm1,xmm3\n"    /* xmm1 = tx + ty + tz + w*/

	/* Store translated vector */
	"movaps  [ecx],xmm1\n"
	:
	: "r" (pMat), "r" (pVec)
	:"eax","ebx","ecx"
	);

}

void VPCALL CSIMD_SSE::mat4D_AlignedVectorMultiplyOf(CVec3D* pOut, const CVec3D* pIn, const CMat4D* pMat)
{
	/*
		Does a normal vector x matrix, but since w = 1.0f, the final column of the matrix is
		merely added.
	*/


	asm volatile(
    LOADREGS2( %2, %1, %0 )
    "movaps    xmm4,[eax]\n"
	"movaps    xmm0,[ebx]\n"
	"movaps  xmm1,[ebx+16]\n"
	"movaps  xmm2,[ebx+32]\n"
	"movaps  xmm3,[ebx+48]\n"
	"movaps  xmm5,xmm4\n"
	"movaps  xmm6,xmm4\n"
	"shufps  xmm4, xmm4,0x00\n"    /* xmm4 = x | x | x | x */
	"shufps  xmm5, xmm5,0x55\n"    /* xmm5 = y | y | y | y */
	"shufps  xmm6, xmm6,0xAA\n"    /* xmm6 = z | z | z | z */

	/* Multiply with each row */
	"mulps  xmm0,xmm4\n"
	"mulps  xmm1,xmm5\n"
	"mulps  xmm2,xmm6\n"

	/* Sum results, xmm0-2 = transformed x,y,z colums, xmm3 = w column */
	"addps  xmm1,xmm0\n"    /* xmm1 = tx + ty */
	"addps  xmm3,xmm2\n"    /* xmm3 = tz + w */
	"addps  xmm1,xmm3\n"    /* xmm1 = tx + ty + tz + w*/

	/* Store translated vector */
	"movaps  [ecx],xmm1\n"
	:
	: "r" (pMat), "r" (pIn), "r" (pOut)
	:"eax","ebx","ecx"
	);

}

void VPCALL CSIMD_SSE::mat4D_AlignedVectorMultiply(CVec4D* pOut4D, const CMat4D* pMat)
{


	asm volatile(
    LOADREGS( %1 , %0 )

	/*
	    IMPLICIT:
		movl 4(%esp), %1
		movl 8(%esp), %0

		Access %1 first, it has finished and has no dependencies, unlike loading from %0

		Surprisingly results in 9% speedup on Celeron, 6% on Athlon64
	*/

	"movaps  xmm4,[ecx]\n"
	"movaps  xmm0,[eax]\n"
	"movaps  xmm1,[eax+16]\n"
	"movaps  xmm2,[eax+32]\n"
	"movaps  xmm3,[eax+48]\n"
	"movaps  xmm5,xmm4\n"
	"movaps  xmm6,xmm4\n"
	"movaps  xmm7,xmm4\n"
	"shufps  xmm4, xmm4,0x00\n"
	"shufps  xmm5, xmm5,0x55\n"
	"shufps  xmm6, xmm6,0xAA\n"
	"shufps  xmm7, xmm7,0xFF\n"

	/* Multiply with each row */
	"mulps  xmm0,xmm4\n"
	"mulps  xmm1,xmm5\n"
	"mulps  xmm2,xmm6\n"
	"mulps  xmm3,xmm7\n"

	/* Sum results */
	"addps  xmm1,xmm0\n"
	"addps  xmm3,xmm2\n"
	"addps  xmm1,xmm3\n"

	/* Store translated vector */
	"movaps  [ecx],xmm1\n"
	:
	: "r" (pMat), "r" (pOut4D)
	:"eax","ebx","ecx"
	);

}


void VPCALL CSIMD_SSE::mat4D_AlignedVectorMultiplyOf(CVec4D* pOut4D, const CVec4D* pIn4D, const CMat4D* pMat)
{

	asm volatile(
	LOADREGS2(%2 , %1, %0 )
    /*
	    IMPLICIT:
		movl 4(%esp), %1
		movl 8(%esp), %0

		Access %1 first, it has finished and has no dependencies, unlike loading from %0

		Surprisingly results in 9% speedup on Celeron, 6% on Athlon64
	*/
	"movaps  xmm4,[eax]\n"
	"movaps    xmm0,[ebx]\n"
	"movaps  xmm1,[ebx+16]\n"
	"movaps  xmm2,[ebx+32]\n"
	"movaps  xmm3,[ebx+48]\n"
	"movaps  xmm5,xmm4\n"
	"movaps  xmm6,xmm4\n"
	"movaps   xmm7,xmm4\n"
	"shufps  xmm4, xmm4,0x00\n"
	"shufps  xmm5, xmm5,0x55\n"
	"shufps  xmm6, xmm6,0xAA\n"
	"shufps  xmm7, xmm7,0xFF\n"
	"mulps  xmm0,xmm4\n"
	"mulps  xmm1,xmm5\n"
	"mulps  xmm2,xmm6\n"
	"mulps  xmm3,xmm7\n"
	"addps  xmm1,xmm0\n"
	"addps  xmm3,xmm2\n"
	"addps  xmm1,xmm3\n"
	"movaps  [ecx],xmm1\n"
	:
	: "r" (pMat), "r" (pIn4D),"r" (pOut4D)
	:"eax","ebx","ecx"
	);

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

void  VPCALL CSIMD_SSE::quaternion_Normalize(CQuaternion* pQuat)
{

		asm(
		"movups  xmm0,[%0]\n"		/* xmm0 = w | z | y | x  */
		"movaps  xmm1,xmm0\n"	/* xmm1 = w | z | y | x */
		"mulps xmm0, xmm0\n"	/* xmm0 = w*w | z*z | y*y | x*x */
		"movhlps  xmm2,xmm0\n"	/* xmm2 = ? | ? | w*w | z*z */
		"addps  xmm2,xmm0\n"	/* xmm2 = ? | ? | w*w+y*y | z*z+x*x */
		"movss  xmm0,xmm2\n"	/* xmm0 = ? | ? | ? | z*z+x*x */
		"shufps  xmm2, xmm2,0x55\n"	/* xmm2 = w*w+y*y |  w*w+y*y |  w*w+y*y | w*w+y*y */
		"addss  xmm2,xmm0\n"
		#ifdef HIPREC
		/* Full division by magnitude */
		"sqrtss  xmm2,xmm2\n"
		"shufps  xmm2, xmm2,0x00\n"
		"divps  xmm1,xmm2\n"
		"movups  [%0],xmm1\n"
		#else
		/* Multiply by reciprocal root approximation */
		"rsqrtss  xmm2,xmm2\n"
		"shufps  xmm2, xmm2,0x00\n"
		"mulps  xmm1,xmm2\n"
		"movups  [%0],xmm1\n"
		#endif

		:
		: "r" (pQuat)
		);



}

void  VPCALL CSIMD_SSE::quaternion_NormalizeOf(CQuaternion* pOut, const CQuaternion* pQuat)
{

		asm(
		"movups  xmm0,[%0]\n"		/* xmm0 = w | z | y | x  */
		"movaps  xmm1,xmm0\n"	/* xmm1 = w | z | y | x */
		"mulps  xmm0,xmm0\n"	/* xmm0 = w*w | z*z | y*y | x*x */
		"movhlps  xmm2,xmm0\n"	/* xmm2 = ? | ? | w*w | z*z */
		"addps  xmm2,xmm0\n"	/* xmm2 = ? | ? | w*w+y*y | z*z+x*x */
		"movss  xmm0,xmm2\n"	/* xmm0 = ? | ? | ? | z*z+x*x */
		"shufps  xmm2, xmm2,0x55\n"	/* xmm2 = w*w+y*y |  w*w+y*y |  w*w+y*y | w*w+y*y */
		"addss  xmm2,xmm0\n"
		#ifdef HIPREC
		/* Full divide by square root */
		"sqrtss  xmm2,xmm2\n"			/* xmm2 = ??? | ??? | ??? | mag */
		"shufps  xmm2, xmm2,0x00\n"	/* xmm2 = mag | mag | mag | mag */
		"divps  xmm1,xmm2\n"			/* xmm1 = w/mag | z/mag | y/mag | x/mag */
		"movups  [%1],xmm1\n"
		#else
		/* Multiply by reciprocal root approximation */
		"rsqrtss  xmm2,xmm2\n"
		"shufps xmm2, xmm2, 0x00\n"
		"mulps  xmm1,xmm2\n"
		"movups [%1],xmm1\n"
		#endif

		:
		: "r" (pQuat), "r" (pOut)
		);




}

void  VPCALL CSIMD_SSE::quaternion_Multiply(CQuaternion* pLeft, const CQuaternion* pRight)
{
	asm(
	"movups  xmm0,[%0]\n" 	/* xmm0 = Rw | Rz | Ry | Rx */
	"movups xmm4, [%1]\n" 	/* xmm1 = Lw | Lz | Ly | Lx */


	/* Duplicate Right throughout xmm0-xmm3 */
	"movaps  xmm1,xmm0\n"	/* xmm1 = xmm0 */
	"movaps  xmm2,xmm0\n"	/* xmm2 = xmm0 */
	"movaps  xmm3,xmm0\n"	/* xmm3 = xmm0 */

	/* Duplicate Left throughout xmm4-xmm7 */
	"movaps  xmm5,xmm4\n"	/* xmm5 = xmm4 */
	"movaps  xmm6,xmm4\n"	/* xmm6 = xmm4 */
	"movaps  xmm7,xmm4\n"	/* xmm7 = xmm4 */

	/*
		Broadcast elements in xmm4-xmm7 to create scalar registers
		==================
		0000 0000 = xxxx = 0x00
	    0101 0101 = yyyy = 0x55
		1010 1010 = zzzz = 0xAA
		1111 1111 = wwww = 0xFF
	*/
	"shufps  xmm4, xmm4,0x00\n"		/* xmm4 = Rx | Rx | Rx | Rx */
	"shufps  xmm5, xmm5,0x55\n"		/* xmm5 = Ry | Ry | Ry | Ry */
	"shufps  xmm6, xmm6,0xAA\n"		/* xmm6 = Rz | Rz | Rz | Rz */
	"shufps  xmm7, xmm7,0xFF\n"		/* xmm7 = Rw | Rw | Rw | Rw */

	/*
		Set up columns
		==============
		C1 = w | z | y | x = 1110 0100 =
		C2 = x | y | z | w = 0001 1011 = 0x1B
		C3 = y | x | w | z = 0100 1110 = 0x4E
		C4 = z | w | x | y = 1011 0001 = 0xB1
	*/

	/* C1 is already w | z | y | x  format, no shufps needed */
	"shufps  xmm1, xmm1,0x1B\n"
	"shufps  xmm2, xmm2,0x4E\n"
	"shufps  xmm3, xmm3,0xB1\n"

	/* Multiply columns */
	"mulps  xmm7,xmm0\n"		/* C1 *= Lw */
	"mulps  xmm4,xmm1\n"		/* C2 *= Lx */
	"mulps  xmm5,xmm2\n"		/* C3 *= Ly */
	"mulps  xmm6,xmm3\n"		/* C4 *= Lz */

	/* Change the signs of the columns (C1, aka: xmm4, doesnt need it)*/
	"xorps  xmm4,%4\n"		/* C2 = { + - + - } */
	"xorps  xmm5,%3\n"		/* C3 = { + + - - } */
	"xorps  xmm6,%2\n"		/* C4 = { - + + - } */

	"addps  xmm5,xmm4\n"		/* C2 += C1 */
	"addps  xmm7,xmm6\n"		/* C4 += C3 */
	"addps  xmm7,xmm5\n"		/* C4 += C2 */

	"movups  [%1],xmm7\n"			/* xmm7 = new quaternion, write it out */

	:
	: "r" (pRight), "r" (pLeft),"m"(NEGPOSPOSNEG),"m"(POSPOSNEGNEG),"m"(POSNEGPOSNEG)
	);


}

void  VPCALL CSIMD_SSE::quaternion_MultiplyOf(CQuaternion* pOut, const CQuaternion* pLeft, const CQuaternion* pRight)
{

	asm(
	"movups  xmm0,[%0]\n" 	/* xmm0 = Rw | Rz | Ry | Rx */
	"movups  xmm4,[%1]\n" 	/* xmm1 = Lw | Lz | Ly | Lx */


	/* Duplicate Right throughout xmm0-xmm3 */
	"movaps  xmm1,xmm0\n"	/* xmm1 = xmm0 */
	"movaps  xmm2,xmm0\n"	/* xmm2 = xmm0 */
	"movaps  xmm3,xmm0\n"	/* xmm3 = xmm0 */

	/* Duplicate Left throughout xmm4-xmm7 */
	"movaps  xmm5,xmm4\n"	/* xmm5 = xmm4 */
	"movaps  xmm6,xmm4\n"	/* xmm6 = xmm4 */
	"movaps  xmm7,xmm4\n"	/* xmm7 = xmm4 */

	/*
		Broadcast elements in xmm4-xmm7 to create scalar registers
		==================
		0000 0000 = xxxx = 0x00
	    0101 0101 = yyyy = 0x55
		1010 1010 = zzzz = 0xAA
		1111 1111 = wwww = 0xFF
	*/
	"shufps  xmm4, xmm4,0x00\n"		/* xmm4 = Rx | Rx | Rx | Rx */
	"shufps  xmm5, xmm5,0x55\n"		/* xmm5 = Ry | Ry | Ry | Ry */
	"shufps  xmm6, xmm6,0xAA\n"		/* xmm6 = Rz | Rz | Rz | Rz */
	"shufps  xmm7, xmm7,0xFF\n"	/* xmm7 = Rw | Rw | Rw | Rw */

	/*
		Set up columns
		==============
		C1 = w | z | y | x = 1110 0100 =
		C2 = x | y | z | w = 0001 1011 = 0x1B
		C3 = y | x | w | z = 0100 1110 = 0x4E
		C4 = z | w | x | y = 1011 0001 = 0xB1
	*/

	/* C1 is already w | z | y | x  format, no shufps needed */
	"shufps  xmm1, xmm1,0x1B\n"
	"shufps  xmm2, xmm2,0x4E\n"
	"shufps  xmm3, xmm3,0xB1\n"

	/* Multiply columns */
	"mulps  xmm7,xmm0\n"		/* C1 *= Lw */
	"mulps  xmm4,xmm1\n"		/* C2 *= Lx */
	"mulps  xmm5,xmm2\n"		/* C3 *= Ly */
	"mulps  xmm6,xmm3\n"		/* C4 *= Lz */

	/* Change the signs of the columns (C1, aka: xmm4, doesnt need it)*/
	"xorps  xmm4, %5\n"		/* C2 = { + - + - } */
	"xorps  xmm5, %4\n"		/* C3 = { + + - - } */
	"xorps  xmm6, %3\n"		/* C4 = { - + + - } */

	"addps  xmm5,xmm4\n"		/* C2 += C1 */
	"addps  xmm7,xmm6\n"		/* C4 += C3 */
	"addps  xmm7,xmm5\n"		/* C4 += C2 */

	"movups  [%2],xmm7\n"			/* xmm7 = new quaternion, write it out */
	:
	: "r" (pRight), "r" (pLeft), "r" (pOut),"m"(NEGPOSPOSNEG),"m"(POSPOSNEGNEG),"m"(POSNEGPOSNEG)
	);



}
//================end Quaternion==========================
//==================PLANE============================

//========CPlane==========================================



void  VPCALL CSIMD_SSE::plane_FromPoints(CPlane* pOut, const CVec3D* pA, const CVec3D* pB,const CVec3D* pC)
{
	asm(
	"movups  xmm0,[%0]\n"	/* pA */
	"movups  xmm1,[%1]\n"	/* pB */
	"movups  xmm2,[%2]\n"	/* pC */
	"movaps  xmm7,xmm0\n"  /* Save 'a' into xmm7 */
	"subps  xmm0,xmm1\n"
	"subps  xmm2,xmm1\n"

	/* Now, just need the cross product of xmm0 and xmm2 */
	"movaps  xmm1,xmm0\n"	/* xmm0 = xmm1 = Left */
	"movaps  xmm3,xmm2\n"	/* xmm2 = xmm3 = Right */
	"shufps  xmm0, xmm0,0xC9\n"	/* Left.yxz */
	"shufps  xmm1, xmm1,0xD2\n"	/* Left.xzy */
	"shufps  xmm2, xmm2,0xD2\n"	/* Right.xzy */
	"shufps  xmm3, xmm3,0xC9\n"	/* Right.yxz */

	/* Multiply columns 1&2 and 3&4 */
	"mulps  xmm2,xmm0\n"
	"mulps  xmm3,xmm1\n"

	/* Got the cross product, OK */
	"subps  xmm2,xmm3\n"

	/* Begin calculation of 'd' component */


	/* AND off bits 96-127 (w component) */
	"andps  xmm2,%4\n"

	/* save xmm4 = 0 | z | y | x */
	"movaps  xmm4,xmm2\n"

	/* Multiply with point 'a' on the polygon (saved in xmm7): xmm2 = 0 | a.z*z | a.y*y | a.x*x */
	"mulps  xmm2,xmm7\n"


	"movhlps  xmm3,xmm2\n"		/* xmm3 = ?   | ?   | 0   |  z^2 */
	"addss  xmm3,xmm2\n"		/* xmm3 = ?   | ?   | 0   | z^2 + x^2 */
	"shufps  xmm2,xmm2,0x55\n"	/* xmm2 = y^2 | y^2 | y^2 | y^2 */
	"addss  xmm3,xmm2\n"		/* xmm3 = ?   | ?   | ?   | x^2+y^2+z^2 */

	/* Change sign */
	"xorps  xmm3, %5\n" /* xmm3 = ? | ? | ? | -(x^2+y^2+z^2)*/

	/* Move to w component location, mask off xyz, and OR with saved portions */
	"shufps  xmm3, xmm3,0x00\n"				/* xmm3 = -(x^2+y^2+z^2) | -(x^2+y^2+z^2) | -(x^2+y^2+z^2) | -(x^2+y^2+z^2)*/
	"andps  xmm3,%6\n"	/* xmm3 = -(x^2+y^2+z^2) | 0 | 0 | 0 */
	"orps  xmm4,xmm3\n"							/* xmm4 = -(x^2+y^2+z^2) | z | y | x */



	/* Save plane coefficients */
	"movups  [%3],xmm4\n"
	:
	: "r" (pA), "r" (pB), "r" (pC), "r" (pOut),"m"(_SIMDx86_float_SSE_NO_W_MASK),"m"(_SIMDx86_float_NEGPOSPOSPOS),"m"(_SIMDx86_float_SSE_NO_XYZ_MASK)
	);



}

float VPCALL CSIMD_SSE::plane_DistToPoint(const CPlane* pPlane, const CVec3D* pPoint)
{
	float dummy;
	asm(
	"movups  xmm0,[%1]\n" /* xmm0 = pPlane->d | pPlane->c | pPlane->b | pPlane->a */
	"movups  xmm1,[%2]\n" /* xmm1 = ????????? | pPoint->z | pPoint->y | pPoint->x */

	"andps  xmm1,%3\n"   /* xmm1 = 0 | ... */
	"movaps  xmm7,xmm0\n"						/* xmm7 = pPlane... */

	"mulps  xmm1,xmm0\n"                        /* xmm1 = d*0.0 | c*z | b*y | a*x */
	"shufps  xmm7, xmm7,0xFF\n"				/* xmm7 = d | d | d | d */


    "movhlps  xmm2,xmm1\n"      /* xmm2 = ???? | ???? | d*0.0 | z*c */
	"addss  xmm2,xmm1\n"        /* xmm2 = ???? | ???? | ????  | x*a + z*c*/
	"shufps  xmm1, xmm1,0x55\n" /* xmm1 = ???? | ???? | ????  | y*b */
	"andps  xmm1,%4\n"	/* xmm1 = ??? | ??? | ??? | fabsf(dot(pPlane, pPoint)) */
	"addss  xmm1,xmm2\n"        /* xmm1 = ???? | ???? | ????? |  fabsf(dot(pPlane, pPoint)) + pPlane->d */

	"movss  [esp-4],xmm1\n"
	//"flds (esp)\n"

	: "=t" (dummy)
	: "r" (pPlane), "r" (pPoint),"m"(_SIMDx86_float_SSE_NO_W_MASK),"m"(_SIMDx86_float_ABS)
	);

}

float  VPCALL CSIMD_SSE::plane_Dot(const CPlane* pPlane, const CVec3D* pVec)
{
    float dummy;
	asm(
	"movups  xmm1,[%2]\n" /* xmm1 = ?    | V->z | V->x | V->y */
	"movups  xmm0,[%1]\n" /* xmm0 = P->d | P->c | P->b | P->a */
	"andps  xmm1,%3\n" /* 0 | z | y | x */
	"mulps  xmm1,xmm0\n"    /* xmm0 = 0 | P->c*V->z | P->b*V->y | P->a*V->x */


		"shufps  xmm0, xmm0,0xFF\n"
		"movhlps  xmm7,xmm1\n"		/* xmm7 = ?   | ?   | 0   |  z's */
		"addss  xmm7,xmm1\n"		/* xmm7 = ?   | ?   | 0   | z's + x's */
		"shufps  xmm1, xmm1,0x55\n"	/* xmm1 = y's | y's | y's | y's */
		"addss  xmm7,xmm0\n"        /* xmm7  = ?  | ?   | ?   | z's + x's + d*/
		"addss  xmm7,xmm1\n"		/* xmm7 = ?   | ?   | ?   | x's + y's + z's + d */
		"movss  [esp-4],xmm7\n"
		//"flds -4(esp)\n"


	: "=t" (dummy)
	: "r" (pPlane), "r" (pVec),"m"(_SIMDx86_float_SSE_NO_W_MASK)
	);


}

float  VPCALL CSIMD_SSE::plane_Dot4(const CPlane* pPlane, const CVec4D* pVec4)
{
    float dummy;
	asm(
	"movups  xmm0,[%1]\n" /* xmm0 = d | c | b | a */
	"movups  xmm1,[%2]\n" /* xmm1 = w | z | y | z */
	"mulps  xmm1,xmm0\n"	/* xmm1 = dw | cz | by | az */


		"movhlps  xmm7,xmm1\n"		/* xmm7 = ? | ? | dw | cz */
		"addps  xmm7,xmm1\n"        /* xmm7 = ? | ? | dw+by | cz+ax */
		"movaps  xmm1,xmm7\n"       /* xmm1 = ? | ? | dw+by | cz+ax */
		"shufps  xmm1, xmm1,0x55\n"	/* xmm1 = ? | ? | ? | dw+by */
		"addss  xmm7,xmm1\n"        /* xmm7  = ?  | ?   | ?   | dw+cz+by+ax */
		"movss  [esp-4],xmm7\n"
		//"flds -4(esp)\n"


	: "=t" (dummy)
	: "r" (pPlane), "r" (pVec4)
	);


}

float  VPCALL CSIMD_SSE::plane_DotNormal(const CPlane* pPlane, const CVec3D* pVec)
{
	float dummy;
	asm(
	"movups  xmm0,[%1]\n" /* xmm0 = P->d | P->c | P->b | P->a */
	"movups  xmm1,[%2]\n" /* xmm1 = ?    | V->z | V->b | V->a */
	"andps  xmm0,%3\n"	 /* 0 | P1.c | P1.b | P1.a */

	"mulps  xmm0,xmm1\n"    /* xmm0 = 0 | P->c*V->z | P->b*V->y | P->a*V->x */

		"movhlps  xmm1,xmm0\n"		/* xmm1 = ?   | ?   | 0   |  z's */
		"addss  xmm1,xmm0\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
		"shufps  xmm0, xmm0,0x55\n"	/* xmm0 = y's | y's | y's | y's */
		"addss  xmm1,xmm0\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */

		"movss  [esp-4],xmm1\n"
		//"flds -4(esp)\n"

	: "=t" (dummy)
	: "r" (pPlane), "r" (pVec),"m"(_SIMDx86_float_SSE_NO_W_MASK)
	);


}

float  VPCALL CSIMD_SSE::plane_DotPlane(const CPlane* pPlane1, const CPlane* pPlane2)
{
	float dummy;
	asm(
	"movups  xmm0,[%1]\n" /* xmm0 = P1.d | P1.c | P1.b | P1.a */
	"movups  xmm1,[%2]\n" /* xmm1 = P2.d | P2.c | P2.b | P2.a */
	"andps  xmm0,%3\n"	 /* 0 | P1.c | P1.b | P1.a */

	"mulps  xmm0,xmm1\n"    /* xmm0 = 0 | P1.c*P2.c | P1.b*P2.b | P1.a*P2.a */

	"movhlps  xmm1,xmm0\n"		/* xmm1 = ?   | ?   | 0   |  z's */
	"addss  xmm1,xmm0\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
	"shufps  xmm0, xmm0,0x55\n"	/* xmm0 = y's | y's | y's | y's */
	"addss  xmm1,xmm0\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */

	"movss  [esp-4],xmm1\n"
	//"flds -4(esp)\n"
	: "=t" (dummy)
	: "r" (pPlane1), "r" (pPlane2),"m"(_SIMDx86_float_SSE_NO_W_MASK)
	);


}

void  VPCALL CSIMD_SSE::plane_Normalize(CPlane* pOut)
{

	asm(
	"movups  xmm0,[%0]\n"								/* xmm0 = d | c | b | a */
	"movaps  xmm1,xmm0\n"							/* xmm1 = d | c | b | a */
	"andps  xmm0,%1\n"		/* xmm0 = 0 | c | b | a */

	/* Perform dot product */
	"mulps  xmm0,xmm0\n"							/* xmm0 = 0 | c*c | b*b | a*a */

	"movhlps  xmm3,xmm0\n"							/* xmm3 = ? | ? | 0 | c*c */
	"addss  xmm3,xmm0\n"							/* xmm3 = ? | ? | ? | a*a + c*c */
	"shufps  xmm0, xmm0,0x55\n"					/* xmm0 = b*b | b*b | b*b | b*b */
	"addss  xmm0,xmm3\n"							/* xmm0 = ? | ? | ? | a*b + b*b + c*c */


	/* Divide it all by the square root */
	#if defined(HIPREC)
	"sqrtss  xmm0,xmm0\n"	/* xmm0 = sqrtf(a*a+b*b+c*c) */
	"shufps  xmm0,xmm0,0x00\n"	/* xmm0 = mag | mag | mag | mag */
	"divps  xmm1, xmm0\n"	/* xmm1 = d/mag | c/mag | b/mag | a/mag */
	#else
	"rsqrtss  xmm0,xmm0\n"	/* xmm0 = ? | ? | ? | 1.0f / sqrtf(a*a+b*b+c*c) */
	"mulps    xmm1, xmm0\n"	/* xmm1 = d*invmag | c*invmag | b*invmag | a*invmag */
	#endif

	"movups  [%0],xmm1\n"
	:
	: "r" (pOut),"m"(_SIMDx86_float_SSE_NO_W_MASK)
	);
	return;



}

void  VPCALL CSIMD_SSE::plane_NormalizeOf(CPlane* pOut, CPlane* pIn)
{

	asm(
	"movups  xmm0,[%0]\n"								/* xmm0 = d | c | b | a */
	"movaps  xmm1,xmm0\n"							/* xmm1 = d | c | b | a */
	"andps  xmm0,%2\n"		/* xmm0 = 0 | c | b | a */

	/* Perform dot product */
	"mulps  xmm0,xmm0\n"							/* xmm0 = 0 | c*c | b*b | a*a */

	"movhlps  xmm3,xmm0\n"							/* xmm3 = ? | ? | 0 | c*c */
	"addss  xmm3,xmm0\n"							/* xmm3 = ? | ? | ? | a*a + c*c */
	"shufps  xmm0, xmm0,0x55\n"					/* xmm0 = b*b | b*b | b*b | b*b */
	"addss  xmm0,xmm3\n"							/* xmm0 = ? | ? | ? | a*b + b*b + c*c */


	/* Divide it all by the square root */
	#if defined(HIPREC)
	"sqrtss  xmm0,xmm0\n"	/* xmm0 = sqrtf(a*a+b*b+c*c) */
	"shufps  xmm0,xmm0,0x00\n"	/* xmm0 = mag | mag | mag | mag */
	"divps  xmm1, xmm0\n"	/* xmm1 = d/mag | c/mag | b/mag | a/mag */
	#else
	"rsqrtss  xmm0,xmm0\n"	/* xmm0 = ? | ? | ? | 1.0f / sqrtf(a*a+b*b+c*c) */
	"mulps    xmm1, xmm0\n"	/* xmm1 = d*invmag | c*invmag | b*invmag | a*invmag */
	#endif

	"movups  [%1],xmm1\n"
	:
	: "r" (pIn), "r" (pOut),"m"(_SIMDx86_float_SSE_NO_W_MASK)
	);
	return;

}


//=============END PLANE=========================

#else     //masm_ATT



//=================GNU C CODE==========================================

//#define  HIPREC

/*
	TODO:
	Optimize Distance() -- It currently uses already implemented functions
*/


void VPCALL CSIMD_SSE::vector3D_Sum(CVec3D* pOut, const CVec3D* pIn)
{

//testeVar
	/* SSE/SSE2/SSE3 Implementation */

	asm(
	"movups (%0), %%xmm0\n"
	"movups (%1), %%xmm1\n"
	//"movaps  (%2), %%xmm2\n"
	/* Remove w component from one with AND mask */
	"andps %2, %%xmm1\n"
	"addps %%xmm0, %%xmm1\n"
	/* Store*/
	"movlpd %%xmm1, (%0)\n"
	"shufps $0x02, %%xmm1,%%xmm1\n"  //shuffle colocar nos bits 0-31 o 3o byte
	"movss	%%xmm1,8(%0)\n"
	:
	: "r" (pOut), "r" (pIn), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);



}

void VPCALL CSIMD_SSE::vector3D_SumOf(CVec3D* pOut, const CVec3D* pIn1, const CVec3D* pIn2)
{
	/* SSE/SSE2/SSE3 Implementation */
   asm(
	"movups (%0), %%xmm0\n"
	"movups (%1), %%xmm1\n"
	"movaps  %3, %%xmm2\n"
	/* Remove w component from one with AND mask */
	"andps %%xmm2, %%xmm1\n"
	"addps %%xmm0, %%xmm1\n"
	/* Store*/
	"movlpd %%xmm1, (%2)\n"
	"shufps $0x02, %%xmm1,%%xmm1\n"  //shuffle colocar nos bits 0-31 o 3o byte
	"movss	%%xmm1,8(%2)\n"

	:
	: "r" (pIn1), "r" (pIn2), "r" (pOut), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);


}


void VPCALL CSIMD_SSE::vector3D_Diff(CVec3D* pLeft, CVec3D* pRight)
{

	/* SSE/SSE2/SSE3 Implementation */
   asm(
	"movups (%0), %%xmm0\n"		/* xmm0 = pLeft */
	"movups (%1), %%xmm1\n"		/* xmm1 = pRight */
	"movaps  (%2), %%xmm2\n"
    /* Remove w component from one with AND mask */
	"andps %%xmm2, %%xmm1\n"
	"subps %%xmm1, %%xmm0\n"	/* xmm0 = pLeft - pRight */
	/* Store*/
	"movlpd %%xmm0, (%0)\n"
	"shufps $0x02, %%xmm0,%%xmm0\n"  //shuffle colocar nos bits 0-31 o 3o byte
	"movss	%%xmm0,8(%0)\n"
	:
	: "r" (pLeft), "r" (pRight), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);


}

void VPCALL CSIMD_SSE::vector3D_DiffOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight)
{
   asm(
	"movups (%0), %%xmm0\n"		/* xmm0 = pLeft */
	"movups (%1), %%xmm1\n"		/* xmm1 = pRight */
	"movaps  (%3), %%xmm2\n"
    /* Remove w component from one with AND mask */
	"andps %%xmm2, %%xmm1\n"
	"subps %%xmm1, %%xmm0\n"	/* xmm0 = pLeft - pRight */
	/* Store*/
	"movlpd %%xmm0, (%2)\n"
	"shufps $0x02, %%xmm0,%%xmm0\n"  //shuffle colocar nos bits 0-31 o 3o byte
	"movss	%%xmm0,8(%2)\n"
	:
	: "r" (pLeft), "r" (pRight), "r" (pOut), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);



}

void VPCALL CSIMD_SSE::vector3D_Scale(CVec3D* pOut, float scalar)
{
   asm(
	"movss %1, %%xmm0\n"		/* using movss, faster might be movlps since it doesnt clear top 96 bits */
	"movups (%0), %%xmm1\n"
	"movaps  %2, %%xmm2\n"
	/* Remove w component from one with AND mask */
	"andps %%xmm2, %%xmm1\n"
	"shufps $0x00, %%xmm0, %%xmm0\n"
	"mulps %%xmm1, %%xmm0\n"
	/* Store*/
	"movlpd %%xmm0, (%0)\n"
	"shufps $0x02, %%xmm0,%%xmm0\n"  //shuffle colocar nos bits 0-31 o 3o byte
	"movss	%%xmm0,8(%0)\n"
	:
	: "r" (pOut), "m" (scalar), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);


}

void VPCALL CSIMD_SSE::vector3D_ScaleOf(CVec3D* pOut, const CVec3D* pIn, float scalar)
{
   asm(
	"movss %2, %%xmm0\n"		/* using movss, faster might be movlps since it doesnt clear top 96 bits */
	"movups (%1), %%xmm1\n"
	"movaps  %3, %%xmm2\n"
	/* Remove w component from one with AND mask */
	"andps %%xmm2, %%xmm1\n"
	"shufps $0x00, %%xmm0, %%xmm0\n"
	"mulps %%xmm1, %%xmm0\n"
	/* Store*/
	"movlpd %%xmm0, (%0)\n"
	"shufps $0x02, %%xmm0,%%xmm0\n"  //shuffle colocar nos bits 0-31 o 3o byte
	"movss	%%xmm0,8(%0)\n"
	:
	: "r" (pOut), "r" (pIn),  "m" (scalar), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);


}

float  VPCALL CSIMD_SSE::vector3D_(const CVec3D* pSrc1, const CVec3D* pSrc2)
{
   	/* SSE3 Implementation */
   float dummy;
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
		"movups (%0), %%xmm0\n"
		"movups (%1), %%xmm1\n"
        "movaps  %2, %%xmm2\n"

		/* Remove w component from one with AND mask */
		"andps %%xmm2, %%xmm0\n"
        "mulps %%xmm0, %%xmm1\n"

		/* Sum components */
		"haddps %%xmm1, %%xmm1\n"
		"haddps %%xmm1, %%xmm1\n"

		"movss %%xmm1, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		:
		: "r" (pSrc1), "r" (pSrc2), "m" (_SIMDx86_float_SSE_NO_W_MASK)
		);


		#else
		/* SSE/SSE2 Implementation */
		asm(

		/* Move vectors and multiply across */
		"movups (%1), %%xmm0\n"
		"movups (%2), %%xmm1\n"
		"movaps  %3, %%xmm2\n"
        "andps %%xmm2, %%xmm0\n"
        "mulps %%xmm0, %%xmm1\n"

		/* Shift data around in the registers (I loves me some haddps right now!!) */
		"movhlps %%xmm1, %%xmm2\n"		/* xmm2 = ?   | ?   | 0   |  z's */
		"addss %%xmm1, %%xmm2\n"		/* xmm2 = ?   | ?   | 0   | z's + x's */
		"shufps $0x55, %%xmm1, %%xmm1\n"/* xmm1 = y's | y's | y's | y's */
		"addss %%xmm1, %%xmm2\n"		/* xmm2 = ?   | ?   | ?   | x's+y's+z's */

		"movss %%xmm2, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pSrc1), "r" (pSrc2), "m" (_SIMDx86_float_SSE_NO_W_MASK)
		);
		#endif
return dummy;

}


float  VPCALL CSIMD_SSE::vector3D_LengthSq(const CVec3D* pVec)
{
		float dummy;

		#if (USE_SSE >= 3)

		/* SSE3 Implementation */
		asm(
		"movups (%1), %%xmm0\n"             /* xmm0 = x | y | z | w*/
		"movaps  %2, %%xmm2\n"
        "andps %%xmm2, %%xmm0\n"
        "mulps %%xmm0, %%xmm0\n"            /* xmm0 = x*x | y*y | z*z | 0.0f */
		"haddps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"           /* xmm0 = lensq | lensq | lensq | lensq */

		/* Store */
		//"movss %%xmm0, -4(%%esp)\n"
		"movss %%xmm0, (%0)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pVec), "m" (_SIMDx86_float_SSE_NO_W_MASK)
		);

		#else /* Using SSE/SSE2 */

		/* SSE/SSE2 Implementation */
		asm(
		"movups (%1), %%xmm0\n"
		"movaps  (%2), %%xmm2\n"
        "andps %%xmm2, %%xmm0\n"
        "mulps %%xmm0, %%xmm0\n"

		/* Shift data around in the registers (I loves me some haddps right now!!) */
		"movhlps %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   |  z's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
		"shufps $0x55, %%xmm0, %%xmm0\n"/* xmm0 = y's | y's | y's | y's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */

		"movss %%xmm1, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pVec), "m" (_SIMDx86_float_SSE_NO_W_MASK)
		);


		#endif /* USE_SSE == 1 ||  USE_SSE == 2 */

		return dummy;

}
#define HIPREC
#include <xmmintrin.h>
//-----------------------------------------------------------------------------
// SSE implementations of optimized routines:
//-----------------------------------------------------------------------------
float _SSE_Sqrt(float x)
{
	float	root = 0.f;

	__asm__ __volatile__(
		"movss %1,%%xmm2\n"
		"sqrtss %%xmm2,%%xmm1\n"
		//"movss %%xmm1,%0"
       	"movss %%xmm2, -4(%%esp)\n"
		"flds -4(%%esp)\n": "=t" (root)
		: "m" (x)
	);

	return root;
}

// Single iteration NewtonRaphson reciprocal square root:
// 0.5 * rsqrtps * (3 - x * rsqrtps(x) * rsqrtps(x))
// Very low error, and fine to use in place of 1.f / sqrtf(x).

// Intel / Kipps SSE rSqrt.  Significantly faster than above.
float _SSE_RSqrtAccurate(float a)
{
	float half = 0.5f;
	float three = 3.f;
    float x;

#if defined (MSCVER)
	__asm
	{
		movss   xmm3, a;
		movss   xmm1, half;
		movss   xmm2, three;
		rsqrtss xmm0, xmm3;

		mulss   xmm3, xmm0;
		mulss   xmm1, xmm0;
		mulss   xmm3, xmm0;
		subss   xmm2, xmm3;
		mulss   xmm1, xmm2;

		movss   x,    xmm1;
	}
#endif // _WIN32
#if !defined(MASM_INTEL)
	__asm__ __volatile__(
		"movss   %1, %%xmm3 \n"
        "movss   %2, %%xmm1 \n"
        "movss   %3, %%xmm2 \n"
        "rsqrtss %%xmm3, %%xmm0 \n"
        "mulss   %%xmm0, %%xmm3 \n"
        "mulss   %%xmm0, %%xmm1 \n"
        "mulss   %%xmm0, %%xmm3 \n"
        "subss   %%xmm3, %%xmm2 \n"
        "mulss   %%xmm2, %%xmm1 \n"
        "movss   %%xmm1, %0 \n"
		: "=m" (x)
		: "m" (a), "m" (half), "m" (three)
);
#else  //masm=intel
	#error "Not Implemented"
#endif

	return x;
}


// Simple SSE rSqrt.  Usually accurate to around 6 (relative) decimal places
// or so, so ok for closed transforms.  (ie, computing lighting normals)
float _SSE_RSqrtFast(float x)
{

	float rroot;
#if defined(MSCVER)
	_asm
	{
		rsqrtss	xmm0, x
		movss	rroot, xmm0
	}
#endif // _WIN32
#if !defined(MASM_INTEL)
	 asm(
		"rsqrtss %1, %%xmm0 \n"
		//"rsqrtss %%xmm1, %%xmm1\n"
		//"movss %%xmm0, %0 \n\t"
		"movss %%xmm0, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (x)
		: "m" (rroot)
		:
	);
#else  //masm_intel

#endif

	return rroot;
}
inline float
sqrtf ( float x )
{
#ifdef PURE_VECTOR
    __m128  xx = _mm_load_ss( & x );
    __m128  xr = _mm_rsqrt_ss( xx );
    __m128  xt;

    xt = _mm_mul_ss( xr, xr );
    xt = _mm_mul_ss( xt, xx );
    xt = _mm_sub_ss( f3, xt );
    xt = _mm_mul_ss( xt, f05 );
    xr = _mm_mul_ss( xr, xt );

    _mm_store_ss( & x, xr );

    return x;
#else
    float   r;

    _mm_store_ss( & r, _mm_sqrt_ss( _mm_load_ss( & x ) ) );

    r *= ((3.0f - r * r * x) * 0.5f);

    return r;
#endif
}
float SqRt_SSE( const float* x )
    {
        ALIGNTO16  __m128 r0;
        ALIGNTO16  float V[4]={14.0f,0.0f,0.0f,0.0f};
        ALIGNTO16  __m128 v_i;
        //v_i = _mm_loadu_ps(x);
        ALIGNTO16 float teste=14.f;
        r0 = _mm_load_ss( &teste );
        r0 = _mm_sqrt_ss( r0 );

        float y;

        _mm_store_ss( &y, r0 );

        return y;
    }
float  VPCALL CSIMD_SSE::vector3D_Length(const CVec3D* pVec)
{




		ALIGNTO16 float dummy;
		/* SSE/SSE2 Implementation */
		asm(
		"movups (%1), %%xmm0\n"
		"andps %2, %%xmm0\n"    /* xmm0 = x | y | z | 0.0f */
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
		/* rcp( rSqrt(value) ) (This may be very inaccurate) */
		"rsqrtss %%xmm1, %%xmm1\n"
		"rcpss %%xmm1, %%xmm1\n"
		#endif

		"movss %%xmm1, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pVec),"m"(_SIMDx86_float_SSE_NO_W_MASK)
		);
        //float test = __builtin_sqrt (dummy);
		#ifdef HIPREC
//		 __builtin_ia32_sqrtss(raiz);
		//return _SSE_Sqrt(dummy);
        return dummy;
        #else
        return _SSE_RSqrtAccurate(dummy);
        //return dummy;
        #endif
}

void VPCALL CSIMD_SSE::vector3D_Cross(CVec3D* pLeft, const CVec3D* pRight)
{
	asm("movups (%0), %%xmm0\n"		/* xmm0 = pLeft */
	"movups (%1), %%xmm1\n"		/* xmm1 = pRight */
	"movaps (%2), %%xmm4\n"
	/* Remove w component from one with AND mask */
	"andps %%xmm4, %%xmm0\n"
	"andps %%xmm4, %%xmm1\n"
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

	/* Store*/
	"movlpd %%xmm1, (%0)\n"
	"shufps $0x02, %%xmm1,%%xmm1\n"  //shuffle colocar nos bits 0-31 o 3o byte
	"movss	%%xmm1,8(%0)\n"
	:
	: "r" (pLeft), "r" (pRight), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	:);


}

void VPCALL CSIMD_SSE::vector3D_CrossOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight)
{
    asm(
	"movups (%1), %%xmm0\n"		/* xmm0 = pLeft */
	"movups (%2), %%xmm1\n"		/* xmm1 = pRight */
	"movaps %3, %%xmm4\n"
	/* Remove w component from one with AND mask */
	"andps %%xmm4, %%xmm0\n"
	"andps %%xmm4, %%xmm1\n"
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
	/* Store*/
	"movlpd %%xmm1, (%0)\n"
	"shufps $0x02, %%xmm1,%%xmm1\n"  //shuffle colocar nos bits 0-31 o 3o byte
	"movss	%%xmm1,8(%0)\n"
	:
	: "r" (pOut), "r" (pLeft), "r" (pRight), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);


}
void VPCALL CSIMD_SSE::vector3D_Normalize(CVec3D* pVec)
{
		#if USE_SSE >= 3

		/* SSE3 Implementation */
		asm(
		"movups (%0), %%xmm0\n"
		"movaps %1, %%xmm4\n"
        "andps %%xmm4, %%xmm0\n"
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
		: "r" (pVec), "m" (_SIMDx86_float_SSE_NO_W_MASK)
		);
		#else

		/* SSE/SSE2 Implementation */
		asm(
		"movups (%0), %%xmm0\n"
		"movaps %1, %%xmm4\n"
        "andps %%xmm4, %%xmm0\n"
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
		/* Store*/
        "movlpd %%xmm7, (%0)\n"
        "shufps $0x02, %%xmm7,%%xmm7\n"  //shuffle colocar nos bits 0-31 o 3o byte
        "movss	%%xmm7,8(%0)\n"
			:
		: "r" (pVec), "m" (_SIMDx86_float_SSE_NO_W_MASK)
		);


		#endif


}
void VPCALL CSIMD_SSE::vector3D_NormalizeOf(CVec3D* pOut, const CVec3D* pVec)
{

		#if USE_SSE >= 3

		/* SSE3 Implementation */
		asm(
		"movups (%0), %%xmm0\n"
		"movaps %2, %%xmm4\n"
        "andps %%xmm4, %%xmm0\n"
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
		: "r" (pVec), "r" (pOut),"m" (_SIMDx86_float_SSE_NO_W_MASK)
		);
		#else /* SSE/SSE2 Implementation */
		asm(
		"movups (%0), %%xmm0\n"
		"movaps %2, %%xmm4\n"
        "andps %%xmm4, %%xmm0\n"    /* xmm0 = x | y | z | 0.0f */
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

		/* Store*/
        "movlpd %%xmm7, (%1)\n"
        "shufps $0x02, %%xmm7,%%xmm7\n"  //shuffle colocar nos bits 0-31 o 3o byte
        "movss	%%xmm7,8(%1)\n"
			:
		: "r" (pVec), "r" (pOut), "m" (_SIMDx86_float_SSE_NO_W_MASK)
		);

	#endif /* USE_SSE == 1 || USE_SSE == 2 */


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
	Optimize Distance() -- It currently uses already implemented functions
*/

void VPCALL CSIMD_SSE::vector3D_AlignedSum(CVec3D* pOut, const CVec3D* pIn)
{
	asm(
	"movaps (%0), %%xmm0\n"
	"movaps (%1),%%xmm1\n"
	"movaps %2,%%xmm4\n"
	/* Remove w component from one with AND mask */
	"andps %%xmm4, %%xmm1\n"
	"addps %%xmm1, %%xmm0\n"
	/* Store*/
	"movlpd %%xmm0, (%0)\n"
	"shufps $0x02, %%xmm0,%%xmm0\n"  //shuffle colocar nos bits 0-31 o 3o byte
	"movss	%%xmm0,8(%0)\n"
	:
	: "r" (pOut), "r" (pIn), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);

}

void VPCALL CSIMD_SSE::vector3D_AlignedSumOf(CVec3D* pOut, const CVec3D* pIn1, const CVec3D* pIn2)
{
	/* SSE/SSE2/SSE3 Implementation */

	asm(
	"movaps (%0), %%xmm0\n"
	"movaps (%1),%%xmm1\n"
	"movaps %3,%%xmm4\n"
	/* Remove w component from one with AND mask */
	"andps %%xmm4, %%xmm1\n"
	"addps %%xmm1, %%xmm0\n"
	/* Store*/
	"movlpd %%xmm0, (%2)\n"
	"shufps $0x02, %%xmm0,%%xmm0\n"  //shuffle colocar nos bits 0-31 o 3o byte
	"movss	%%xmm0,8(%2)\n"
	:
	: "r" (pIn1), "r" (pIn2), "r" (pOut), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);


}


void VPCALL CSIMD_SSE::vector3D_AlignedDiff(CVec3D* pLeft, CVec3D* pRight)
{
	/* SSE/SSE2/SSE3 Implementation */

	asm(
	"movaps (%0), %%xmm0\n"		/* xmm0 = pLeft */
	"movaps (%1),%%xmm1\n"
	"movaps %2,%%xmm4\n"
	/* Remove w component from one with AND mask */
	"andps %%xmm4, %%xmm1\n"
	"subps %%xmm1, %%xmm0\n"	/* xmm0 = pLeft - pRight */
	/* Store*/
	"movlpd %%xmm0, (%0)\n"
	"shufps $0x02, %%xmm0,%%xmm0\n"  //shuffle colocar nos bits 0-31 o 3o byte
	"movss	%%xmm0,8(%0)\n"
	:
	: "r" (pLeft), "r" (pRight), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);


}

void VPCALL CSIMD_SSE::vector3D_AlignedDiffOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight)
{

	asm(
	"movaps (%0), %%xmm0\n"		/* xmm0 = pLeft */
	"movaps (%1),%%xmm1\n"
	"movaps %3,%%xmm4\n"
	/* Remove w component from one with AND mask */
	"andps %%xmm4, %%xmm1\n"
	"subps %%xmm1, %%xmm0\n"	/* xmm0 = pLeft - pRight */
	/* Store*/
	"movlpd %%xmm0, (%2)\n"
	"shufps $0x02, %%xmm0,%%xmm0\n"  //shuffle colocar nos bits 0-31 o 3o byte
	"movss	%%xmm0,8(%2)\n"
	:
	: "r" (pLeft), "r" (pRight), "r" (pOut), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);



}

void VPCALL CSIMD_SSE::vector3D_AlignedScale(CVec3D* pOut, float scalar)
{

	asm(
	"movss %1, %%xmm0\n"		/* using movss, faster might be movlps since it doesnt clear top 96 bits */
	"shufps $0x00, %%xmm0, %%xmm0\n"
	"movaps (%0),%%xmm1\n"
	"movaps %2,%%xmm4\n"
	/* Remove w component from one with AND mask */
	"andps %%xmm4, %%xmm1\n"
	"mulps %%xmm1, %%xmm0\n"
	/* Store*/
	"movlpd %%xmm0, (%0)\n"
	"shufps $0x02, %%xmm0,%%xmm0\n"  //shuffle colocar nos bits 0-31 o 3o byte
	"movss	%%xmm0,8(%0)\n"
	:
	: "r" (pOut), "m" (scalar), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);


}

void VPCALL CSIMD_SSE::vector3D_AlignedScaleOf(CVec3D* pOut, const CVec3D* pIn, float scalar)
{

	asm(
	"movss %2, %%xmm0\n"		/* using movss, faster might be movlps since it doesnt clear top 96 bits */
	"shufps $0x00, %%xmm0, %%xmm0\n"
	"movaps (%1),%%xmm1\n"
	"movaps %3,%%xmm4\n"
	/* Remove w component from one with AND mask */
	"andps %%xmm4, %%xmm1\n"
	"mulps %%xmm1, %%xmm0\n"
	/* Store*/
	"movlpd %%xmm0, (%0)\n"
	"shufps $0x02, %%xmm0,%%xmm0\n"  //shuffle colocar nos bits 0-31 o 3o byte
	"movss	%%xmm0,8(%0)\n"
	:
	: "r" (pOut), "r" (pIn),  "m" (scalar), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);


}

float  VPCALL CSIMD_SSE::vector3D_AlignedDot(const CVec3D* pSrc1, const CVec3D* pSrc2)
{
float dummy;
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
		"movaps %2,%%xmm4\n"
        "mulps (%1), %%xmm0\n"

		/* Remove w component from one with AND mask */
		"andps %%xmm4, %%xmm0\n"

		/* Sum components */
		"haddps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"

		"movss %%xmm0, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		:
		: "r" (pSrc1), "r" (pSrc2), "m" (_SIMDx86_float_SSE_NO_W_MASK)
		);


		#else
		/* SSE/SSE2 Implementation */
		asm(

		/* Move vectors and multiply across */
		"movaps (%1), %%xmm0\n"
		"movaps %3,%%xmm4\n"
        "mulps (%2), %%xmm0\n"
		"andps %%xmm4, %%xmm0\n"

		/* Shift data around in the registers (I loves me some haddps right now!!) */
		"movhlps %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   |  z's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
		"shufps $0x55, %%xmm0, %%xmm0\n"/* xmm0 = y's | y's | y's | y's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */


		"movss %%xmm1, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pSrc1), "r" (pSrc2), "m" (_SIMDx86_float_SSE_NO_W_MASK)
		);
		return dummy;
		#endif

return dummy;
}


float  VPCALL CSIMD_SSE::vector3D_AlignedLengthSq(const CVec3D* pVec)
{


		float dummy;

		#if (USE_SSE >= 3)

		/* SSE3 Implementation */
		asm(
		"movaps (%1), %%xmm0\n"             /* xmm0 = x | y | z | w*/
		"movaps %2,%%xmm4\n"
        "andps %%xmm4, %%xmm0\n"    /* xmm0 = x | y | z | 0.0f */
		"mulps %%xmm0, %%xmm0\n"            /* xmm0 = x*x | y*y | z*z | 0.0f */
		"haddps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"           /* xmm0 = lensq | lensq | lensq | lensq */

		/* Store */
		"movss %%xmm0, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pVec), "m" (_SIMDx86_float_SSE_NO_W_MASK)
		);

		#else /* Using SSE/SSE2 */

		/* SSE/SSE2 Implementation */
		asm(
		"movaps (%1), %%xmm0\n"
		"movaps %2,%%xmm4\n"
        "andps %%xmm4, %%xmm0\n"    /* xmm0 = x | y | z | 0.0f */
		"mulps %%xmm0, %%xmm0\n"



		/* Shift data around in the registers (I loves me some haddps right now!!) */
		"movhlps %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   |  z's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
		"shufps $0x55, %%xmm0, %%xmm0\n"/* xmm0 = y's | y's | y's | y's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */


		"movss %%xmm1, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pVec), "m" (_SIMDx86_float_SSE_NO_W_MASK)
		);


		#endif /* USE_SSE == 1 ||  USE_SSE == 2 */

		return dummy;

}

float  VPCALL CSIMD_SSE::vector3D_AlignedLength(const CVec3D* pVec)
{


		#if USE_SSE >= 3
		float dummy;

		/* SSE3 Implementation */
		asm(
		"movaps (%1), %%xmm0\n"
		"movaps %2,%%xmm4\n"
        "andps %%xmm4, %%xmm0\n"    /* xmm0 = x | y | z | 0.0f */
		"mulps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"
		"haddps %%xmm0, %%xmm0\n"

		#ifdef HIPREC
		/* Full square root */
		"sqrtss %%xmm0, %%xmm0\n"
		#else
		/* rcp( rSqrt(value) ) (This may be very inaccurate) */
		"rsqrtss %%xmm0, %%xmm0\n"
		"rcpss %%xmm0, %%xmm0\n"
		#endif

		/* Store */
		"movss %%xmm0, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pVec), "m" (_SIMDx86_float_SSE_NO_W_MASK)
		);

		return dummy;
		#else

		float dummy;
		/* SSE/SSE2 Implementation */
		asm(
		"movaps (%1), %%xmm0\n"
		"movaps %2,%%xmm4\n"
        "andps %%xmm4, %%xmm0\n"    /* xmm0 = x | y | z | 0.0f */
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
		/* rcp( rSqrt(value) ) (This may be very inaccurate) */
		"rsqrtss %%xmm1, %%xmm1\n"
		"rcpss %%xmm1, %%xmm1\n"
		#endif

		"movss %%xmm1, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pVec), "m" (_SIMDx86_float_SSE_NO_W_MASK)
		);

		return dummy;
		#endif

}

void VPCALL CSIMD_SSE::vector3D_AlignedCross(CVec3D* pLeft, const CVec3D* pRight)
{

	asm(
	"movaps (%0), %%xmm0\n"		/* xmm0 = pLeft */
	"movaps (%1), %%xmm1\n"		/* xmm1 = pRight */
	"movaps %2,%%xmm4\n"
	"andps %%xmm4, %%xmm0\n"
	"andps %%xmm4, %%xmm1\n"
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

	/* Store*/
	"movlpd %%xmm1, (%0)\n"
	"shufps $0x02, %%xmm1,%%xmm1\n"  //shuffle colocar nos bits 0-31 o 3o byte
	"movss	%%xmm1,8(%0)\n"
	:
	: "r" (pLeft), "r" (pRight), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);


}

void VPCALL CSIMD_SSE::vector3D_AlignedCrossOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight)
{

	asm(
	"movaps (%1), %%xmm0\n"		/* xmm0 = pLeft */
	"movaps (%2), %%xmm1\n"		/* xmm1 = pRight */
	"movaps %3,%%xmm4\n"
	"andps %%xmm4, %%xmm0\n"
	"andps %%xmm4, %%xmm1\n"
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
	/* Store*/
	"movlpd %%xmm1, (%0)\n"
	"shufps $0x02, %%xmm1,%%xmm1\n"  //shuffle colocar nos bits 0-31 o 3o byte
	"movss	%%xmm1,8(%0)\n"
	:
	: "r" (pOut), "r" (pLeft), "r" (pRight), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);


}
void VPCALL CSIMD_SSE::vector3D_AlignedNormalize(CVec3D* pVec)
{
		#if USE_SSE >= 3

		/* SSE3 Implementation */
		asm(
		"movaps (%0), %%xmm0\n"
		"movaps %1,%%xmm4\n"
        "andps %%xmm4, %%xmm0\n"    /* xmm0 = x | y | z | 0.0f */
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
		: "r" (pVec), "m" (_SIMDx86_float_SSE_NO_W_MASK)
		);
		#else

		/* SSE/SSE2 Implementation */
		asm(
		"movaps (%0), %%xmm0\n"
		"movaps %1,%%xmm4\n"
        "andps %%xmm4, %%xmm0\n"    /* xmm0 = x | y | z | 0.0f */
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
		/* Store*/
        "movlpd %%xmm7, (%0)\n"
        "shufps $0x02, %%xmm7,%%xmm7\n"  //shuffle colocar nos bits 0-31 o 3o byte
        "movss	%%xmm7,8(%0)\n"
			:
		: "r" (pVec), "m" (_SIMDx86_float_SSE_NO_W_MASK)
		);


		#endif


}
void VPCALL CSIMD_SSE::vector3D_AlignedNormalizeOf(CVec3D* pOut, const CVec3D* pVec)
{

		#if USE_SSE >= 3

		/* SSE3 Implementation */
		asm(
		"movaps (%0), %%xmm0\n"
		"movaps %2,%%xmm4\n"
        "andps %%xmm4, %%xmm0\n"    /* xmm0 = x | y | z | 0.0f */
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
        "movss	%%xmm2,8(%1)\n"

		:
		: "r" (pVec), "r" (pOut), "m" (_SIMDx86_float_SSE_NO_W_MASK)
		);
		#else /* SSE/SSE2 Implementation */
		asm(
		"movaps (%0), %%xmm0\n"
		"movaps %2,%%xmm4\n"
        "andps %%xmm4, %%xmm0\n"    /* xmm0 = x | y | z | 0.0f */
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
		/* Store*/
        "movlpd %%xmm7, (%1)\n"
        "shufps $0x02, %%xmm7,%%xmm7\n"  //shuffle colocar nos bits 0-31 o 3o byte
        "movss	%%xmm7,8(%1)\n"
			:
		: "r" (pVec), "r" (pOut), "m" (_SIMDx86_float_SSE_NO_W_MASK)
		);

	#endif /* USE_SSE == 1 || USE_SSE == 2 */


}

float  VPCALL CSIMD_SSE::vector3D_AlignedDistance(const CVec3D* pVec1, const CVec3D* pVec2)
{
#if 0 /* Not so good just yet... */

	#if defined(USE_SSE)

	#if USE_SSE == 3
		float dummy;
		asm(
			"movaps (%1), %%xmm0\n"
			"movaps %3,%%xmm4\n"
            "subps (%2), %%xmm0\n"
			"andps %%xmm4, %%xmm0\n"
			"haddps %%xmm0, %%xmm0\n"
			"haddps %%xmm0, %%xmm0\n"
			"subl $4, %%esp\n"
			"movss %%xmm0, (%%esp)\n"
			"flds (%%esp)\n"
			"addl $4, %%esp\n"
			: "=t" (dummy)
			: "r" (pVec1), "r" (pVec2), "m" (_SIMDx86_float_SSE_NO_W_MASK)
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
//======================CVec4D ===============================


/*
	TODO:
	Optimize Distance() -- It currently uses already implemented functions
*/


void VPCALL CSIMD_SSE::vector4D_Sum(CVec4D* pOut, const CVec4D* pIn)
{

	/* SSE/SSE2/SSE3 Implementation */

	asm(
	"movups (%0), %%xmm0\n"
	"movups (%1), %%xmm1\n"
	"addps %%xmm0, %%xmm1\n"
	"movups %%xmm1, (%0)\n"
	:
	: "r" (pOut), "r" (pIn)
	);

}

void VPCALL CSIMD_SSE::vector4D_SumOf(CVec4D* pOut, const CVec4D* pIn1, const CVec4D* pIn2)
{
	/* SSE/SSE2/SSE3 Implementation */

	asm(
	"movups (%0), %%xmm0\n"
	"movups (%1), %%xmm1\n"
	"addps %%xmm0, %%xmm1\n"
	"movups %%xmm1, (%2)\n"
	:
	: "r" (pIn1), "r" (pIn2), "r" (pOut)
	);


}


void VPCALL CSIMD_SSE::vector4D_Diff(CVec4D* pLeft, CVec4D* pRight)
{

	/* SSE/SSE2/SSE3 Implementation */

	asm(
	"movups (%0), %%xmm0\n"		/* xmm0 = pLeft */
	"movups (%1), %%xmm1\n"		/* xmm1 = pRight */
	"subps %%xmm1, %%xmm0\n"	/* xmm0 = pLeft - pRight */
	"movups %%xmm0, (%0)\n"		/* store */
	:
	: "r" (pLeft), "r" (pRight)
	);


}

void VPCALL CSIMD_SSE::vector4D_DiffOf(CVec4D* pOut, const CVec4D* pLeft, const CVec4D* pRight)
{

	asm(
	"movups (%0), %%xmm0\n"		/* xmm0 = pLeft */
	"movups (%1), %%xmm1\n"		/* xmm1 = pRight */
	"subps %%xmm1, %%xmm0\n"	/* xmm0 = pLeft - pRight */
	"movups %%xmm0, (%2)\n"		/* store */
	:
	: "r" (pLeft), "r" (pRight), "r" (pOut)
	);



}

void VPCALL CSIMD_SSE::vector4D_Scale(CVec4D* pOut, float scalar)
{

	asm(
	"movss %1, %%xmm0\n"		/* using movss, faster might be movlps since it doesnt clear top 96 bits */
	"movups (%0), %%xmm1\n"
	"shufps $0x00, %%xmm0, %%xmm0\n"
	"mulps %%xmm1, %%xmm0\n"
	"movups %%xmm0, (%0)\n"
	:
	: "r" (pOut), "m" (scalar)
	);


}

void VPCALL CSIMD_SSE::vector4D_ScaleOf(CVec4D* pOut, const CVec4D* pIn, float scalar)
{

	asm(
	"movss %2, %%xmm0\n"		/* using movss, faster might be movlps since it doesnt clear top 96 bits */
	"movups (%1), %%xmm1\n"
	"shufps $0x00, %%xmm0, %%xmm0\n"
	"mulps %%xmm1, %%xmm0\n"
	"movups %%xmm0, (%0)\n"
	:
	: "r" (pOut), "r" (pIn),  "m" (scalar)
	);


}


float  VPCALL CSIMD_SSE::vector4D_Dot(const CVec4D* pSrc4D1, const CVec4D* pSrc4D2)
{


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


}

float  VPCALL CSIMD_SSE::vector4D_LengthSq(const CVec4D* pVec)
{

		float dummy;
        /* Using SSE/SSE2 */

		/* SSE/SSE2 Implementation */
		asm(
		"movl   %1,%%eax\n"
		"movups (%%eax), %%xmm0\n"
		"mulps %%xmm0, %%xmm0\n"
        /* Shift data around in the registers (I loves me some haddps right now!!) */
		"movaps  %%xmm0, %%xmm7\n"
		"shufps  $0x03,%%xmm7,%%xmm7\n"  /* xmm7 = ?   | ?   | ?  |  w's *///colocar nos bis 0-31 de xmm7 os bits 96-127 de xmm7

		/* Shift data around in the registers (I loves me some haddps right now!!) */
		"movhlps %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   |  z's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
		"shufps $0x55, %%xmm0, %%xmm0\n"/* xmm0 = y's | y's | y's | y's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
        "addss %%xmm7, %%xmm1\n"        /* xmm1 = ?   | ?   | ?   | x's+y's+z's+w's */

		"movss %%xmm1, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pVec)
		);


		 /* USE_SSE == 1 ||  USE_SSE == 2 */

		return dummy;

}

float  VPCALL CSIMD_SSE::vector4D_Length(const CVec4D* pVec)
{


		float dummy;
		/* SSE/SSE2 Implementation */
		asm(
		"movups (%1), %%xmm0\n"
		"mulps %%xmm0, %%xmm0\n"

        /* Shift data around in the registers (I loves me some haddps right now!!) */
		"movaps  %%xmm0, %%xmm7\n"
		"shufps  $0x03,%%xmm7,%%xmm7\n"  /* xmm7 = ?   | ?   | ?  |  w's *///colocar nos bis 0-31 de xmm7 os bits 96-127 de xmm7

		/* Shift data around in the registers (I loves me some haddps right now!!) */
		"movhlps %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   |  z's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
		"shufps $0x55, %%xmm0, %%xmm0\n"/* xmm0 = y's | y's | y's | y's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
        "addss %%xmm7, %%xmm1\n"        /* xmm1 = ?   | ?   | ?   | x's+y's+z's+w's */

		#ifdef HIPREC
		/* Full square root */
		"sqrtss %%xmm1, %%xmm1\n"
		#else
		/* rcp( rSqrt(value) ) (This may be very inaccurate) */
		"rsqrtss %%xmm1, %%xmm1\n"
		"rcpss %%xmm1, %%xmm1\n"
		#endif

		"movss %%xmm1, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pVec)
		);

		return dummy;

}


void VPCALL CSIMD_SSE::vector4D_Normalize(CVec4D* pVec)
{


		/* SSE/SSE2 Implementation */
		asm(
		"movups (%0), %%xmm0\n"
		"movaps %%xmm0, %%xmm7\n"			/* Save for the division by length */
		"mulps %%xmm0, %%xmm0\n"

        /* Shift data around in the registers (I loves me some haddps right now!!) */
		"movaps  %%xmm0, %%xmm2\n"
		"shufps  $0x03,%%xmm2,%%xmm2\n"  /* xmm7 = ?   | ?   | ?  |  w's *///colocar nos bis 0-31 de xmm7 os bits 96-127 de xmm7

		/* Shift data around in the registers (I loves me some haddps right now!!) */
		"movhlps %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   |  z's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
		"shufps $0x55, %%xmm0, %%xmm0\n"/* xmm0 = y's | y's | y's | y's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
        "addss %%xmm2, %%xmm1\n"        /* xmm1 = ?   | ?   | ?   | x's+y's+z's+w's */

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



}
void VPCALL CSIMD_SSE::vector4D_NormalizeOf(CVec4D* pOut, const CVec4D* pVec)
{

 /* SSE/SSE2 Implementation */
		asm(
        "movl   %0,%%eax\n"
        "movl   %1,%%ecx\n"
		"movups (%%eax), %%xmm0\n"
		"movaps %%xmm0, %%xmm7\n"			/* Save for the division by length */
		"mulps %%xmm0, %%xmm0\n"

        /* Shift data around in the registers (I loves me some haddps right now!!) */
		"movaps  %%xmm0, %%xmm2\n"
		"shufps  $0x03,%%xmm2,%%xmm2\n"  /* xmm7 = ?   | ?   | ?  |  w's *///colocar nos bis 0-31 de xmm7 os bits 96-127 de xmm7

		/* Shift data around in the registers (I loves me some haddps right now!!) */
		"movhlps %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   |  z's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
		"shufps $0x55, %%xmm0, %%xmm0\n"/* xmm0 = y's | y's | y's | y's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */
        "addss %%xmm2, %%xmm1\n"        /* xmm1 = ?   | ?   | ?   | x's+y's+z's+w's */

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
		"movups %%xmm7, (%%ecx)\n"
				:
		: "r" (pVec), "r" (pOut)
		);

 /* USE_SSE == 1 || USE_SSE == 2 */


}

float  VPCALL CSIMD_SSE::vector4D_Distance(const CVec4D* pVec1, const CVec4D* pVec2)
{
	/* TODO: Optimize me completely */
	CVec4D diff;
	vector4D_DiffOf(&diff, pVec1, pVec2);
	return vector4D_Length(&diff);
}

/*
	TODO:
	Optimize Distance() -- It currently uses already implemented functions
*/

void VPCALL CSIMD_SSE::vector4D_AlignedSum(CVec4D* pOut, const CVec4D* pIn)
{


	asm(
	"movaps (%0), %%xmm0\n"
	"addps (%1), %%xmm0\n"
	"movaps %%xmm0, (%0)\n"
	:
	: "r" (pOut), "r" (pIn)
	);

}

void VPCALL CSIMD_SSE::vector4D_AlignedSumOf(CVec4D* pOut, const CVec4D* pIn1, const CVec4D* pIn2)
{
	/* SSE/SSE2/SSE3 Implementation */

	asm(
	"movaps (%0), %%xmm0\n"
	"addps (%1), %%xmm0\n"
	"movaps %%xmm0, (%2)\n"
	:
	: "r" (pIn1), "r" (pIn2), "r" (pOut)
	);


}


void VPCALL CSIMD_SSE::vector4D_AlignedDiff(CVec4D* pLeft, CVec4D* pRight)
{

	/* SSE/SSE2/SSE3 Implementation */

	asm(
	"movaps (%0), %%xmm0\n"		/* xmm0 = pLeft */
	"subps (%1), %%xmm0\n"	/* xmm0 = pLeft - pRight */
	"movaps %%xmm0, (%0)\n"		/* store */
	:
	: "r" (pLeft), "r" (pRight)
	);


}

void VPCALL CSIMD_SSE::vector4D_AlignedDiffOf(CVec4D* pOut, const CVec4D* pLeft, const CVec4D* pRight)
{

	asm(
	"movaps (%0), %%xmm0\n"		/* xmm0 = pLeft */
	"subps (%1), %%xmm0\n"	/* xmm0 = pLeft - pRight */
	"movaps %%xmm0, (%2)\n"		/* store */
	:
	: "r" (pLeft), "r" (pRight), "r" (pOut)
	);



}

void VPCALL CSIMD_SSE::vector4D_AlignedScale(CVec4D* pOut, float scalar)
{

	asm(
	"movss %1, %%xmm0\n"		/* using movss, faster might be movlps since it doesnt clear top 96 bits */
	"shufps $0x00, %%xmm0, %%xmm0\n"
	"mulps (%0), %%xmm0\n"
	"movaps %%xmm0, (%0)\n"
	:
	: "r" (pOut), "m" (scalar)
	);


}

void VPCALL CSIMD_SSE::vector4D_AlignedScaleOf(CVec4D* pOut, const CVec4D* pIn, float scalar)
{

	asm(
	"movss %2, %%xmm0\n"		/* using movss, faster might be movlps since it doesnt clear top 96 bits */
	"shufps $0x00, %%xmm0, %%xmm0\n"
	"mulps (%1), %%xmm0\n"
	"movaps %%xmm0, (%0)\n"
	:
	: "r" (pOut), "r" (pIn),  "m" (scalar)
	);


}


float  VPCALL CSIMD_SSE::vector4D_AlignedDot(const CVec4D* pSrc4D1, const CVec4D* pSrc4D2)
{


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



}

float  VPCALL CSIMD_SSE::vector4D_AlignedLengthSq(const CVec4D* pVec)
{


		float dummy;
        /* Using SSE/SSE2 */

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



		return dummy;

}

float  VPCALL CSIMD_SSE::vector4D_AlignedLength(const CVec4D* pVec)
{




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
		/* rcp( rSqrt(value) ) (This may be very inaccurate) */
		"rsqrtss %%xmm1, %%xmm1\n"
		"rcpss %%xmm1, %%xmm1\n"
		#endif

		"movss %%xmm1, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pVec)
		);

		return dummy;


}


void VPCALL CSIMD_SSE::vector4D_AlignedNormalize(CVec4D* pVec)
{


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




}
void VPCALL CSIMD_SSE::vector4D_AlignedNormalizeOf(CVec4D* pOut, const CVec4D* pVec)
{

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

//========CMat==============================

void VPCALL CSIMD_SSE::mat4D_Sum(CMat4D* Out, const CMat4D* In)
{
	asm(
	"movups   (%0), %%xmm0\n"
	"movups 16(%0), %%xmm1\n"
	"movups 32(%0), %%xmm2\n"
	"movups 48(%0), %%xmm3\n"
	"movups   (%1), %%xmm4\n"
	"movups 16(%1), %%xmm5\n"
	"movups 32(%1), %%xmm6\n"
	"movups 48(%1), %%xmm7\n"
	"addps %%xmm0, %%xmm4\n"
	"addps %%xmm1, %%xmm5\n"
	"addps %%xmm2, %%xmm6\n"
	"addps %%xmm3, %%xmm7\n"
	"movups %%xmm4,   (%0)\n"
	"movups %%xmm5, 16(%0)\n"
	"movups %%xmm6, 32(%0)\n"
	"movups %%xmm7, 48(%0)\n"
	:
	: "r" (Out), "r" (In)
	);

}


void VPCALL CSIMD_SSE::mat4D_SumOf(CMat4D* Out, const CMat4D* In1, const CMat4D* In2)
{


	asm(
	"movups   (%0), %%xmm0\n"
	"movups 16(%0), %%xmm1\n"
	"movups 32(%0), %%xmm2\n"
	"movups 48(%0), %%xmm3\n"
	"movups   (%1), %%xmm4\n"
	"movups 16(%1), %%xmm5\n"
	"movups 32(%1), %%xmm6\n"
	"movups 48(%1), %%xmm7\n"
	"addps %%xmm0, %%xmm4\n"
	"addps %%xmm1, %%xmm5\n"
	"addps %%xmm2, %%xmm6\n"
	"addps %%xmm3, %%xmm7\n"
	"movups %%xmm4,   (%2)\n"
	"movups %%xmm5, 16(%2)\n"
	"movups %%xmm6, 32(%2)\n"
	"movups %%xmm7, 48(%2)\n"
	:
	: "r" (In1), "r" (In2), "r" (Out)
	);

}

void VPCALL CSIMD_SSE::mat4D_Diff(CMat4D* Out, const CMat4D* In)
{

	asm(
	"movups   (%0), %%xmm0\n"
	"movups 16(%0), %%xmm1\n"
	"movups 32(%0), %%xmm2\n"
	"movups 48(%0), %%xmm3\n"
	"movups   (%1), %%xmm4\n"
	"movups 16(%1), %%xmm5\n"
	"movups 32(%1), %%xmm6\n"
	"movups 48(%1), %%xmm7\n"
	"subps %%xmm4, %%xmm0\n"
	"subps %%xmm5, %%xmm1\n"
	"subps %%xmm6, %%xmm2\n"
	"subps %%xmm7, %%xmm3\n"
	"movups %%xmm0,   (%0)\n"
	"movups %%xmm1, 16(%0)\n"
	"movups %%xmm2, 32(%0)\n"
	"movups %%xmm3, 48(%0)\n"
	:
	: "r" (Out), "r" (In)
	);
}

void VPCALL CSIMD_SSE::mat4D_DiffOf(CMat4D* Out, const CMat4D* In1, const CMat4D* In2)
{
	asm(
	"movups   (%0), %%xmm0\n"
	"movups 16(%0), %%xmm1\n"
	"movups 32(%0), %%xmm2\n"
	"movups 48(%0), %%xmm3\n"
	"movups   (%1), %%xmm4\n"
	"movups 16(%1), %%xmm5\n"
	"movups 32(%1), %%xmm6\n"
	"movups 48(%1), %%xmm7\n"
	"subps %%xmm4, %%xmm0\n"
	"subps %%xmm5, %%xmm1\n"
	"subps %%xmm6, %%xmm2\n"
	"subps %%xmm7, %%xmm3\n"
	"movups %%xmm0,   (%2)\n"
	"movups %%xmm1, 16(%2)\n"
	"movups %%xmm2, 32(%2)\n"
	"movups %%xmm3, 48(%2)\n"
	:
	: "r" (In1), "r" (In2), "r" (Out)
	);
}

void VPCALL CSIMD_SSE::mat4D_Scale(CMat4D* mtx, float scalar)
{

	asm(
	/* Store scalar in xmm4.x */
	"movss %1, %%xmm4\n"

	/* Get the matrix into registers */
	"movups   (%0), %%xmm0\n"
	"movups 16(%0), %%xmm1\n"
	"movups 32(%0), %%xmm2\n"
	"movups 48(%0), %%xmm3\n"

	/* Broadcast element x to yzw to make a duplicated scalar register */
	"shufps $0x00, %%xmm4, %%xmm4\n"

	/* Scale the matrix in parallel */
	"mulps %%xmm4, %%xmm0\n"
	"mulps %%xmm4, %%xmm1\n"
	"mulps %%xmm4, %%xmm2\n"
	"mulps %%xmm4, %%xmm3\n"

	/* Store results */
	"movups %%xmm0,   (%0)\n"
	"movups %%xmm1, 16(%0)\n"
	"movups %%xmm2, 32(%0)\n"
	"movups %%xmm3, 48(%0)\n"
	:
	: "r" (mtx), "m" (scalar)
	);
}

void VPCALL CSIMD_SSE::mat4D_ScaleOf(CMat4D* pOut, const CMat4D* pIn, float scalar)
{

	asm(
	"movss %1, %%xmm4\n"
	"movups   (%0), %%xmm0\n"
	"movups 16(%0), %%xmm1\n"
	"movups 32(%0), %%xmm2\n"
	"movups 48(%0), %%xmm3\n"
	"shufps $0x00, %%xmm4, %%xmm4\n"
	"mulps %%xmm4, %%xmm0\n"
	"mulps %%xmm4, %%xmm1\n"
	"mulps %%xmm4, %%xmm2\n"
	"mulps %%xmm4, %%xmm3\n"
	"movups %%xmm0,   (%2)\n"
	"movups %%xmm1, 16(%2)\n"
	"movups %%xmm2, 32(%2)\n"
	"movups %%xmm3, 48(%2)\n"
	:
	: "r" (pIn), "m" (scalar), "r" (pOut)
	);
}

void VPCALL CSIMD_SSE::mat4D_Multiply(CMat4D* pLeft, const CMat4D* pRight)
{
	asm(
	"movups   (%0), %%xmm0\n"	/* xmm0 = pRight[0..3] */
	"movups 16(%0), %%xmm1\n"	/* xmm1 = pRight[5..7] */
	"movups 32(%0), %%xmm2\n"	/* xmm2 = pRight[8..11] */
	"movups 48(%0), %%xmm3\n"	/* xmm3 = pRight[12..15] */

	/* Processes 1/2 of the matrix at a time (2x4), unrolled loop */
	"movss    (%1), %%xmm4\n"
	"movss   4(%1), %%xmm6\n"
	"movss  16(%1), %%xmm5\n"
	"movss  20(%1), %%xmm7\n"
	"shufps $0x00, %%xmm4, %%xmm4\n"
	"shufps $0x00, %%xmm5, %%xmm5\n"
	"shufps $0x00, %%xmm6, %%xmm6\n"
	"shufps $0x00, %%xmm7, %%xmm7\n"
	"mulps %%xmm0, %%xmm4\n"
	"mulps %%xmm0, %%xmm5\n"
	"mulps %%xmm1, %%xmm6\n"
	"mulps %%xmm1, %%xmm7\n"
	"addps %%xmm7, %%xmm5\n"
	"addps %%xmm6, %%xmm4\n"


	"movss  8(%1), %%xmm6\n"
	"movss 24(%1), %%xmm7\n"
	"shufps $0x00, %%xmm6, %%xmm6\n"
	"shufps $0x00, %%xmm7, %%xmm7\n"
	"mulps %%xmm2, %%xmm6\n"
	"mulps %%xmm2, %%xmm7\n"
	"addps %%xmm6, %%xmm4\n"
	"addps %%xmm7, %%xmm5\n"

	"movss  12(%1), %%xmm6\n"
	"movss  28(%1), %%xmm7\n"
	"shufps $0x00, %%xmm6, %%xmm6\n"
	"shufps $0x00, %%xmm7, %%xmm7\n"
	"mulps %%xmm3, %%xmm6\n"
	"mulps %%xmm3, %%xmm7\n"
	"addps %%xmm6, %%xmm4\n"
	"addps %%xmm7, %%xmm5\n"

	"movups %%xmm4, (%1)\n"
	"movups %%xmm5, 16(%1)\n"

	/* second half of the matrix */
	"movss  32(%1), %%xmm4\n"
	"movss  36(%1), %%xmm6\n"
	"movss  48(%1), %%xmm5\n"
	"movss  52(%1), %%xmm7\n"
	"shufps $0x00, %%xmm4, %%xmm4\n"
	"shufps $0x00, %%xmm5, %%xmm5\n"
	"mulps %%xmm0, %%xmm4\n"
	"mulps %%xmm0, %%xmm5\n"

	"shufps $0x00, %%xmm6, %%xmm6\n"
	"shufps $0x00, %%xmm7, %%xmm7\n"
	"mulps %%xmm1, %%xmm6\n"
	"mulps %%xmm1, %%xmm7\n"
	"addps %%xmm6, %%xmm4\n"
	"addps %%xmm7, %%xmm5\n"


	"movss 40(%1), %%xmm6\n"
	"movss 56(%1), %%xmm7\n"
	"shufps $0x00, %%xmm6, %%xmm6\n"
	"shufps $0x00, %%xmm7, %%xmm7\n"
	"mulps %%xmm2, %%xmm6\n"
	"mulps %%xmm2, %%xmm7\n"
	"addps %%xmm6, %%xmm4\n"
	"addps %%xmm7, %%xmm5\n"

	"movss  44(%1), %%xmm6\n"
	"movss  60(%1), %%xmm7\n"
	"shufps $0x00, %%xmm6, %%xmm6\n"
	"shufps $0x00, %%xmm7, %%xmm7\n"
	"mulps %%xmm3, %%xmm6\n"
	"mulps %%xmm3, %%xmm7\n"
	"addps %%xmm6, %%xmm4\n"
	"addps %%xmm7, %%xmm5\n"

	"movups %%xmm4, 32(%1)\n"
	"movups %%xmm5, 48(%1)\n"

	:
	: "r" (pRight), "r" (pLeft)
	);

}

void VPCALL CSIMD_SSE::mat4D_MultiplyOf(CMat4D* pOut, const CMat4D* pLeft, const CMat4D* pRight)
{

	asm(
	"movups   (%0), %%xmm0\n"	/* xmm0 = pRight[0..3] */
	"movups 16(%0), %%xmm1\n"	/* xmm1 = pRight[5..7] */
	"movups 32(%0), %%xmm2\n"	/* xmm2 = pRight[8..11] */
	"movups 48(%0), %%xmm3\n"	/* xmm3 = pRight[12..15] */

	/* Processes 1/2 of the matrix at a time (2x4), unrolled loop */
	"movss    (%1), %%xmm4\n"
	"movss   4(%1), %%xmm6\n"
	"movss  16(%1), %%xmm5\n"
	"movss  20(%1), %%xmm7\n"
	"shufps $0x00, %%xmm4, %%xmm4\n"
	"shufps $0x00, %%xmm5, %%xmm5\n"
	"shufps $0x00, %%xmm6, %%xmm6\n"
	"shufps $0x00, %%xmm7, %%xmm7\n"
	"mulps %%xmm0, %%xmm4\n"
	"mulps %%xmm0, %%xmm5\n"
	"mulps %%xmm1, %%xmm6\n"
	"mulps %%xmm1, %%xmm7\n"
	"addps %%xmm7, %%xmm5\n"
	"addps %%xmm6, %%xmm4\n"


	"movss  8(%1), %%xmm6\n"
	"movss 24(%1), %%xmm7\n"
	"shufps $0x00, %%xmm6, %%xmm6\n"
	"shufps $0x00, %%xmm7, %%xmm7\n"
	"mulps %%xmm2, %%xmm6\n"
	"mulps %%xmm2, %%xmm7\n"
	"addps %%xmm6, %%xmm4\n"
	"addps %%xmm7, %%xmm5\n"

	"movss  12(%1), %%xmm6\n"
	"movss  28(%1), %%xmm7\n"
	"shufps $0x00, %%xmm6, %%xmm6\n"
	"shufps $0x00, %%xmm7, %%xmm7\n"
	"mulps %%xmm3, %%xmm6\n"
	"mulps %%xmm3, %%xmm7\n"
	"addps %%xmm6, %%xmm4\n"
	"addps %%xmm7, %%xmm5\n"

	"movups %%xmm4, (%2)\n"
	"movups %%xmm5, 16(%2)\n"

	/* second half of the matrix */
	"movss  32(%1), %%xmm4\n"
	"movss  36(%1), %%xmm6\n"
	"movss  48(%1), %%xmm5\n"
	"movss  52(%1), %%xmm7\n"
	"shufps $0x00, %%xmm4, %%xmm4\n"
	"shufps $0x00, %%xmm5, %%xmm5\n"
	"mulps %%xmm0, %%xmm4\n"
	"mulps %%xmm0, %%xmm5\n"

	"shufps $0x00, %%xmm6, %%xmm6\n"
	"shufps $0x00, %%xmm7, %%xmm7\n"
	"mulps %%xmm1, %%xmm6\n"
	"mulps %%xmm1, %%xmm7\n"
	"addps %%xmm6, %%xmm4\n"
	"addps %%xmm7, %%xmm5\n"


	"movss 40(%1), %%xmm6\n"
	"movss 56(%1), %%xmm7\n"
	"shufps $0x00, %%xmm6, %%xmm6\n"
	"shufps $0x00, %%xmm7, %%xmm7\n"
	"mulps %%xmm2, %%xmm6\n"
	"mulps %%xmm2, %%xmm7\n"
	"addps %%xmm6, %%xmm4\n"
	"addps %%xmm7, %%xmm5\n"

	"movss  44(%1), %%xmm6\n"
	"movss  60(%1), %%xmm7\n"
	"shufps $0x00, %%xmm6, %%xmm6\n"
	"shufps $0x00, %%xmm7, %%xmm7\n"
	"mulps %%xmm3, %%xmm6\n"
	"mulps %%xmm3, %%xmm7\n"
	"addps %%xmm6, %%xmm4\n"
	"addps %%xmm7, %%xmm5\n"

	"movups %%xmm4, 32(%2)\n"
	"movups %%xmm5, 48(%2)\n"

	:
	: "r" (pRight), "r" (pLeft), "r" (pOut)
	);

}



void VPCALL CSIMD_SSE::mat4D_Transpose(CMat4D* pIn)
{

	asm(
		"movlps	  (%0), %%xmm1\n"
		"movlps	 8(%0), %%xmm3\n"
		"movhps	16(%0), %%xmm1\n"
		"movhps	24(%0), %%xmm3\n"
		"movlps	32(%0), %%xmm5\n"
		"movlps	40(%0), %%xmm4\n"
		"movhps	48(%0), %%xmm5\n"
		"movhps	56(%0), %%xmm4\n"
		"movaps	%%xmm1, %%xmm0\n"
		"movaps	%%xmm3, %%xmm2\n"
		"shufps	$0xDD, %%xmm5,%%xmm1\n"
		"shufps	$0xDD, %%xmm4,%%xmm3\n"
		"shufps	$0x88, %%xmm5,%%xmm0\n"
		"shufps	$0x88, %%xmm4,%%xmm2\n"
		"movups %%xmm0,   (%0)\n"
		"movups %%xmm1, 16(%0)\n"
		"movups %%xmm2, 32(%0)\n"
		"movups %%xmm3, 48(%0)\n"
		:
		: "r" (pIn)
		);


}



void VPCALL CSIMD_SSE::mat4D_TransposeOf(CMat4D* pOut, const CMat4D* pIn)
{
	asm(
		"movlps	  (%0), %%xmm1\n"
		"movlps	 8(%0), %%xmm3\n"
		"movhps	16(%0), %%xmm1\n"
		"movhps	24(%0), %%xmm3\n"
		"movlps	32(%0), %%xmm5\n"
		"movlps	40(%0), %%xmm4\n"
		"movhps	48(%0), %%xmm5\n"
		"movhps	56(%0), %%xmm4\n"
		"movaps	%%xmm1, %%xmm0\n"
		"movaps	%%xmm3, %%xmm2\n"
		"shufps	$0xDD, %%xmm5,%%xmm1\n"
		"shufps	$0xDD, %%xmm4,%%xmm3\n"
		"shufps	$0x88, %%xmm5,%%xmm0\n"
		"shufps	$0x88, %%xmm4,%%xmm2\n"
		"movups %%xmm0,   (%1)\n"
		"movups %%xmm1, 16(%1)\n"
		"movups %%xmm2, 32(%1)\n"
		"movups %%xmm3, 48(%1)\n"
		:
		: "r" (pIn), "r" (pOut)
		);

}

void VPCALL CSIMD_SSE::mat4D_VectorMultiply(CVec3D* pVec, const CMat4D* pMat)
{
	/*
		Does a normal vector x matrix, but since w = 1.0f, the final column of the matrix is
		merely added.
	*/


	asm(
	"movups   (%1), %%xmm4\n"
	"movups   (%0), %%xmm0\n"
	"movups 16(%0), %%xmm1\n"
	"movups 32(%0), %%xmm2\n"
	"movups 48(%0), %%xmm3\n"
	"movaps %%xmm4, %%xmm5\n"
	"movaps %%xmm4, %%xmm6\n"
	"shufps $0x00, %%xmm4, %%xmm4\n"    /* xmm4 = x | x | x | x */
	"shufps $0x55, %%xmm5, %%xmm5\n"    /* xmm5 = y | y | y | y */
	"shufps $0xAA, %%xmm6, %%xmm6\n"    /* xmm6 = z | z | z | z */

	/* Multiply with each row */
	"mulps %%xmm4, %%xmm0\n"
	"mulps %%xmm5, %%xmm1\n"
	"mulps %%xmm6, %%xmm2\n"

	/* Sum results, xmm0-2 = transformed x,y,z colums, xmm3 = w column */
	"addps %%xmm0, %%xmm1\n"    /* xmm1 = tx + ty */
	"addps %%xmm2, %%xmm3\n"    /* xmm3 = tz + w */
	"addps %%xmm3, %%xmm1\n"    /* xmm1 = tx + ty + tz + w*/

	/* Store translated vector */
	"movups %%xmm1, (%1)\n"
	:
	: "r" (pMat), "r" (pVec)
	);

}

void VPCALL CSIMD_SSE::mat4D_VectorMultiplyOf(CVec3D* pOut, const CVec3D* pIn, const CMat4D* pMat)
{
	/*
		Does a normal vector x matrix, but since w = 1.0f, the final column of the matrix is
		merely added.
	*/


	asm(
	"movups   (%1), %%xmm4\n"
	"movups   (%0), %%xmm0\n"
	"movups 16(%0), %%xmm1\n"
	"movups 32(%0), %%xmm2\n"
	"movups 48(%0), %%xmm3\n"
	"movaps %%xmm4, %%xmm5\n"
	"movaps %%xmm4, %%xmm6\n"
	"shufps $0x00, %%xmm4, %%xmm4\n"    /* xmm4 = x | x | x | x */
	"shufps $0x55, %%xmm5, %%xmm5\n"    /* xmm5 = y | y | y | y */
	"shufps $0xAA, %%xmm6, %%xmm6\n"    /* xmm6 = z | z | z | z */

	/* Multiply with each row */
	"mulps %%xmm4, %%xmm0\n"
	"mulps %%xmm5, %%xmm1\n"
	"mulps %%xmm6, %%xmm2\n"

	/* Sum results, xmm0-2 = transformed x,y,z colums, xmm3 = w column */
	"addps %%xmm0, %%xmm1\n"    /* xmm1 = tx + ty */
	"addps %%xmm2, %%xmm3\n"    /* xmm3 = tz + w */
	"addps %%xmm3, %%xmm1\n"    /* xmm1 = tx + ty + tz + w*/

	/* Store translated vector */
	"movups %%xmm1, (%2)\n"
	:
	: "r" (pMat), "r" (pIn), "r" (pOut)
	);

}

void VPCALL CSIMD_SSE::mat4D_VectorMultiply(CVec4D* pOut4D, const CMat4D* pMat)
{

	asm(

	/*
	    IMPLICIT:
		movl 4(%esp), %1
		movl 8(%esp), %0

		Access %1 first, it has finished and has no dependencies, unlike loading from %0

		Surprisingly results in 9% speedup on Celeron, 6% on Athlon64
	*/

	"movups (%1), %%xmm4\n"
	"movups   (%0), %%xmm0\n"
	"movups 16(%0), %%xmm1\n"
	"movups 32(%0), %%xmm2\n"
	"movups 48(%0), %%xmm3\n"
	"movaps %%xmm4, %%xmm5\n"
	"movaps %%xmm4, %%xmm6\n"
	"movaps  %%xmm4, %%xmm7\n"
	"shufps $0x00, %%xmm4, %%xmm4\n"
	"shufps $0x55, %%xmm5, %%xmm5\n"
	"shufps $0xAA, %%xmm6, %%xmm6\n"
	"shufps $0xFF, %%xmm7, %%xmm7\n"

	/* Multiply with each row */
	"mulps %%xmm4, %%xmm0\n"
	"mulps %%xmm5, %%xmm1\n"
	"mulps %%xmm6, %%xmm2\n"
	"mulps %%xmm7, %%xmm3\n"

	/* Sum results */
	"addps %%xmm0, %%xmm1\n"
	"addps %%xmm2, %%xmm3\n"
	"addps %%xmm3, %%xmm1\n"

	/* Store translated vector */
	"movups %%xmm1, (%1)\n"
	:
	: "r" (pMat), "r" (pOut4D)
	);

}


void VPCALL CSIMD_SSE::mat4D_VectorMultiplyOf(CVec4D* pOut4D, const CVec4D* pIn4D, const CMat4D* pMat)
{
	asm(
	/*
	    IMPLICIT:
		movl 4(%esp), %1
		movl 8(%esp), %0

		Access %1 first, it has finished and has no dependencies, unlike loading from %0

		Surprisingly results in 9% speedup on Celeron, 6% on Athlon64
	*/
	"movups (%1), %%xmm4\n"
	"movups   (%0), %%xmm0\n"
	"movups 16(%0), %%xmm1\n"
	"movups 32(%0), %%xmm2\n"
	"movups 48(%0), %%xmm3\n"
	"movaps %%xmm4, %%xmm5\n"
	"movaps %%xmm4, %%xmm6\n"
	"movaps  %%xmm4, %%xmm7\n"
	"shufps $0x00, %%xmm4, %%xmm4\n"
	"shufps $0x55, %%xmm5, %%xmm5\n"
	"shufps $0xAA, %%xmm6, %%xmm6\n"
	"shufps $0xFF, %%xmm7, %%xmm7\n"
	"mulps %%xmm4, %%xmm0\n"
	"mulps %%xmm5, %%xmm1\n"
	"mulps %%xmm6, %%xmm2\n"
	"mulps %%xmm7, %%xmm3\n"
	"addps %%xmm0, %%xmm1\n"
	"addps %%xmm2, %%xmm3\n"
	"addps %%xmm3, %%xmm1\n"
	"movups %%xmm1, (%2)\n"
	:
	: "r" (pMat), "r" (pIn4D),"r" (pOut4D)
	);

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

//=====================Quaternion===================================

void  VPCALL CSIMD_SSE::quaternion_Normalize(CQuaternion* pQuat)
{

		asm(
		"movups (%0), %%xmm0\n"		/* xmm0 = w | z | y | x  */
		"movaps %%xmm0, %%xmm1\n"	/* xmm1 = w | z | y | x */
		"mulps %%xmm0, %%xmm0\n"	/* xmm0 = w*w | z*z | y*y | x*x */
		"movhlps %%xmm0, %%xmm2\n"	/* xmm2 = ? | ? | w*w | z*z */
		"addps %%xmm0, %%xmm2\n"	/* xmm2 = ? | ? | w*w+y*y | z*z+x*x */
		"movss %%xmm2, %%xmm0\n"	/* xmm0 = ? | ? | ? | z*z+x*x */
		"shufps $0x55, %%xmm2, %%xmm2\n"	/* xmm2 = w*w+y*y |  w*w+y*y |  w*w+y*y | w*w+y*y */
		"addss %%xmm0, %%xmm2\n"
		#ifdef HIPREC
		/* Full division by magnitude */
		"sqrtss %%xmm2, %%xmm2\n"
		"shufps $0x00, %%xmm2, %%xmm2\n"
		"divps %%xmm2, %%xmm1\n"
		"movups %%xmm1, (%0)\n"
		#else
		/* Multiply by reciprocal root approximation */
		"rsqrtss %%xmm2, %%xmm2\n"
		"shufps $0x00, %%xmm2, %%xmm2\n"
		"mulps %%xmm2, %%xmm1\n"
		"movups %%xmm1, (%0)\n"
		#endif

		:
		: "r" (pQuat)
		);



}

void  VPCALL CSIMD_SSE::quaternion_NormalizeOf(CQuaternion* pOut, const CQuaternion* pQuat)
{

		asm(
		"movups (%0), %%xmm0\n"		/* xmm0 = w | z | y | x  */
		"movaps %%xmm0, %%xmm1\n"	/* xmm1 = w | z | y | x */
		"mulps %%xmm0, %%xmm0\n"	/* xmm0 = w*w | z*z | y*y | x*x */
		"movhlps %%xmm0, %%xmm2\n"	/* xmm2 = ? | ? | w*w | z*z */
		"addps %%xmm0, %%xmm2\n"	/* xmm2 = ? | ? | w*w+y*y | z*z+x*x */
		"movss %%xmm2, %%xmm0\n"	/* xmm0 = ? | ? | ? | z*z+x*x */
		"shufps $0x55, %%xmm2, %%xmm2\n"	/* xmm2 = w*w+y*y |  w*w+y*y |  w*w+y*y | w*w+y*y */
		"addss %%xmm0, %%xmm2\n"
		#ifdef HIPREC
		/* Full divide by square root */
		"sqrtss %%xmm2, %%xmm2\n"			/* xmm2 = ??? | ??? | ??? | mag */
		"shufps $0x00, %%xmm2, %%xmm2\n"	/* xmm2 = mag | mag | mag | mag */
		"divps %%xmm2, %%xmm1\n"			/* xmm1 = w/mag | z/mag | y/mag | x/mag */
		"movups %%xmm1, (%1)\n"
		#else
		/* Multiply by reciprocal root approximation */
		"rsqrtss %%xmm2, %%xmm2\n"
		"shufps $0x00, %%xmm2, %%xmm2\n"
		"mulps %%xmm2, %%xmm1\n"
		"movups %%xmm1, (%1)\n"
		#endif

		:
		: "r" (pQuat), "r" (pOut)
		);




}

void  VPCALL CSIMD_SSE::quaternion_Multiply(CQuaternion* pLeft, const CQuaternion* pRight)
{
	asm(
	"movups (%0), %%xmm0\n" 	/* xmm0 = Rw | Rz | Ry | Rx */
	"movups (%1), %%xmm4\n" 	/* xmm1 = Lw | Lz | Ly | Lx */


	/* Duplicate Right throughout xmm0-xmm3 */
	"movaps %%xmm0, %%xmm1\n"	/* xmm1 = xmm0 */
	"movaps %%xmm0, %%xmm2\n"	/* xmm2 = xmm0 */
	"movaps %%xmm0, %%xmm3\n"	/* xmm3 = xmm0 */

	/* Duplicate Left throughout xmm4-xmm7 */
	"movaps %%xmm4, %%xmm5\n"	/* xmm5 = xmm4 */
	"movaps %%xmm4, %%xmm6\n"	/* xmm6 = xmm4 */
	"movaps %%xmm4, %%xmm7\n"	/* xmm7 = xmm4 */

	/*
		Broadcast elements in xmm4-xmm7 to create scalar registers
		==================
		0000 0000 = xxxx = 0x00
	    0101 0101 = yyyy = 0x55
		1010 1010 = zzzz = 0xAA
		1111 1111 = wwww = 0xFF
	*/
	"shufps $0x00, %%xmm4, %%xmm4\n"		/* xmm4 = Rx | Rx | Rx | Rx */
	"shufps $0x55, %%xmm5, %%xmm5\n"		/* xmm5 = Ry | Ry | Ry | Ry */
	"shufps $0xAA, %%xmm6, %%xmm6\n"		/* xmm6 = Rz | Rz | Rz | Rz */
	"shufps $0xFF, %%xmm7, %%xmm7\n"		/* xmm7 = Rw | Rw | Rw | Rw */

	/*
		Set up columns
		==============
		C1 = w | z | y | x = 1110 0100 =
		C2 = x | y | z | w = 0001 1011 = 0x1B
		C3 = y | x | w | z = 0100 1110 = 0x4E
		C4 = z | w | x | y = 1011 0001 = 0xB1
	*/

	/* C1 is already w | z | y | x  format, no shufps needed */
	"shufps $0x1B, %%xmm1, %%xmm1\n"
	"shufps $0x4E, %%xmm2, %%xmm2\n"
	"shufps $0xB1, %%xmm3, %%xmm3\n"

	/* Multiply columns */
	"mulps %%xmm0, %%xmm7\n"		/* C1 *= Lw */
	"mulps %%xmm1, %%xmm4\n"		/* C2 *= Lx */
	"mulps %%xmm2, %%xmm5\n"		/* C3 *= Ly */
	"mulps %%xmm3, %%xmm6\n"		/* C4 *= Lz */

	/* Change the signs of the columns (C1, aka: xmm4, doesnt need it)*/
	"xorps _POSNEGPOSNEG, %%xmm4\n"		/* C2 = { + - + - } */
	"xorps _POSPOSNEGNEG, %%xmm5\n"		/* C3 = { + + - - } */
	"xorps _NEGPOSPOSNEG, %%xmm6\n"		/* C4 = { - + + - } */

	"addps %%xmm4, %%xmm5\n"		/* C2 += C1 */
	"addps %%xmm6, %%xmm7\n"		/* C4 += C3 */
	"addps %%xmm5, %%xmm7\n"		/* C4 += C2 */

	"movups %%xmm7, (%1)\n"			/* xmm7 = new quaternion, write it out */

	:
	: "r" (pRight), "r" (pLeft)
	);


}

void  VPCALL CSIMD_SSE::quaternion_MultiplyOf(CQuaternion* pOut, const CQuaternion* pLeft, const CQuaternion* pRight)
{

	asm(
	"movups (%0), %%xmm0\n" 	/* xmm0 = Rw | Rz | Ry | Rx */
	"movups (%1), %%xmm4\n" 	/* xmm1 = Lw | Lz | Ly | Lx */


	/* Duplicate Right throughout xmm0-xmm3 */
	"movaps %%xmm0, %%xmm1\n"	/* xmm1 = xmm0 */
	"movaps %%xmm0, %%xmm2\n"	/* xmm2 = xmm0 */
	"movaps %%xmm0, %%xmm3\n"	/* xmm3 = xmm0 */

	/* Duplicate Left throughout xmm4-xmm7 */
	"movaps %%xmm4, %%xmm5\n"	/* xmm5 = xmm4 */
	"movaps %%xmm4, %%xmm6\n"	/* xmm6 = xmm4 */
	"movaps %%xmm4, %%xmm7\n"	/* xmm7 = xmm4 */

	/*
		Broadcast elements in xmm4-xmm7 to create scalar registers
		==================
		0000 0000 = xxxx = 0x00
	    0101 0101 = yyyy = 0x55
		1010 1010 = zzzz = 0xAA
		1111 1111 = wwww = 0xFF
	*/
	"shufps $0x00, %%xmm4, %%xmm4\n"		/* xmm4 = Rx | Rx | Rx | Rx */
	"shufps $0x55, %%xmm5, %%xmm5\n"		/* xmm5 = Ry | Ry | Ry | Ry */
	"shufps $0xAA, %%xmm6, %%xmm6\n"		/* xmm6 = Rz | Rz | Rz | Rz */
	"shufps $0xFF, %%xmm7, %%xmm7\n"		/* xmm7 = Rw | Rw | Rw | Rw */

	/*
		Set up columns
		==============
		C1 = w | z | y | x = 1110 0100 =
		C2 = x | y | z | w = 0001 1011 = 0x1B
		C3 = y | x | w | z = 0100 1110 = 0x4E
		C4 = z | w | x | y = 1011 0001 = 0xB1
	*/

	/* C1 is already w | z | y | x  format, no shufps needed */
	"shufps $0x1B, %%xmm1, %%xmm1\n"
	"shufps $0x4E, %%xmm2, %%xmm2\n"
	"shufps $0xB1, %%xmm3, %%xmm3\n"

	/* Multiply columns */
	"mulps %%xmm0, %%xmm7\n"		/* C1 *= Lw */
	"mulps %%xmm1, %%xmm4\n"		/* C2 *= Lx */
	"mulps %%xmm2, %%xmm5\n"		/* C3 *= Ly */
	"mulps %%xmm3, %%xmm6\n"		/* C4 *= Lz */

	/* Change the signs of the columns (C1, aka: xmm4, doesnt need it)*/
	"xorps _POSNEGPOSNEG, %%xmm4\n"		/* C2 = { + - + - } */
	"xorps _POSPOSNEGNEG, %%xmm5\n"		/* C3 = { + + - - } */
	"xorps _NEGPOSPOSNEG, %%xmm6\n"		/* C4 = { - + + - } */

	"addps %%xmm4, %%xmm5\n"		/* C2 += C1 */
	"addps %%xmm6, %%xmm7\n"		/* C4 += C3 */
	"addps %%xmm5, %%xmm7\n"		/* C4 += C2 */

	"movups %%xmm7, (%2)\n"			/* xmm7 = new quaternion, write it out */
	:
	: "r" (pRight), "r" (pLeft), "r" (pOut)
	);



}

//=============PLANE=======================

void plane_FromPoints(CPlane* pOut, const CVec3D* pA, const CVec3D* pB,const CVec3D* pC)
{
	asm(
	"movups (%0), %%xmm0\n"	/* pA */
	"movups (%1), %%xmm1\n"	/* pB */
	"movups (%2), %%xmm2\n"	/* pC */
	"movaps %%xmm0, %%xmm7\n"   /* Save 'a' into xmm7 */
	"subps %%xmm1, %%xmm0\n"
	"subps %%xmm1, %%xmm2\n"

	/* Now, just need the cross product of xmm0 and xmm2 */
	"movaps %%xmm0, %%xmm1\n"	/* xmm0 = xmm1 = Left */
	"movaps %%xmm2, %%xmm3\n"	/* xmm2 = xmm3 = Right */
	"shufps $0xC9, %%xmm0, %%xmm0\n"	/* Left.yxz */
	"shufps $0xD2, %%xmm1, %%xmm1\n"	/* Left.xzy */
	"shufps $0xD2, %%xmm2, %%xmm2\n"	/* Right.xzy */
	"shufps $0xC9, %%xmm3, %%xmm3\n"	/* Right.yxz */

	/* Multiply columns 1&2 and 3&4 */
	"mulps %%xmm0, %%xmm2\n"
	"mulps %%xmm1, %%xmm3\n"

	/* Got the cross product, OK */
	"subps %%xmm3, %%xmm2\n"

	/* Begin calculation of 'd' component */


	/* AND off bits 96-127 (w component) */
	"andps _SIMDx86_float_SSE_NO_W_MASK, %%xmm2\n"

	/* save xmm4 = 0 | z | y | x */
	"movaps %%xmm2, %%xmm4\n"

	/* Multiply with point 'a' on the polygon (saved in xmm7): xmm2 = 0 | a.z*z | a.y*y | a.x*x */
	"mulps %%xmm7, %%xmm2\n"


	"movhlps %%xmm2, %%xmm3\n"		/* xmm3 = ?   | ?   | 0   |  z^2 */
	"addss %%xmm2, %%xmm3\n"		/* xmm3 = ?   | ?   | 0   | z^2 + x^2 */
	"shufps $0x55, %%xmm2, %%xmm2\n"/* xmm2 = y^2 | y^2 | y^2 | y^2 */
	"addss %%xmm2, %%xmm3\n"		/* xmm3 = ?   | ?   | ?   | x^2+y^2+z^2 */

	/* Change sign */
	"xorps _SIMDx86_float_NEGPOSPOSPOS, %%xmm3\n" /* xmm3 = ? | ? | ? | -(x^2+y^2+z^2)*/

	/* Move to w component location, mask off xyz, and OR with saved portions */
	"shufps $0x00, %%xmm3, %%xmm3\n"				/* xmm3 = -(x^2+y^2+z^2) | -(x^2+y^2+z^2) | -(x^2+y^2+z^2) | -(x^2+y^2+z^2)*/
	"andps _SIMDx86_float_SSE_NO_XYZ_MASK, %%xmm3\n"/* xmm3 = -(x^2+y^2+z^2) | 0 | 0 | 0 */
	"orps %%xmm3, %%xmm4\n"							/* xmm4 = -(x^2+y^2+z^2) | z | y | x */



	/* Save plane coefficients */
	"movups %%xmm4, (%3)\n"
	:
	: "r" (pA), "r" (pB), "r" (pC), "r" (pOut)
	);



}

float plane_DistToPoint(const CPlane* pPlane, const CVec3D* pPoint)
{
	float dummy;
	asm(
	"movups (%1), %%xmm0\n" /* xmm0 = pPlane->d | pPlane->c | pPlane->b | pPlane->a */
	"movups (%2), %%xmm1\n" /* xmm1 = ????????? | pPoint->z | pPoint->y | pPoint->x */

	"andps SIMDx86_float_SSE_NO_W_MASK, %%xmm1\n"   /* xmm1 = 0 | ... */
	"movaps %%xmm0, %%xmm7\n"						/* xmm7 = pPlane... */

	"mulps %%xmm0, %%xmm1\n"                        /* xmm1 = d*0.0 | c*z | b*y | a*x */
	"shufps $0xFF, %%xmm7, %%xmm7\n"				/* xmm7 = d | d | d | d */


    "movhlps %%xmm1, %%xmm2\n"      /* xmm2 = ???? | ???? | d*0.0 | z*c */
	"addss %%xmm1, %%xmm2\n"        /* xmm2 = ???? | ???? | ????  | x*a + z*c*/
	"shufps $x55, %%xmm1, %%xmm1\n" /* xmm1 = ???? | ???? | ????  | y*b */
	"andps _SIMDx86_float_ABS, %%xmm1\n"	/* xmm1 = ??? | ??? | ??? | fabsf(dot(pPlane, pPoint)) */
	"addss %%xmm2, %%xmm1\n"        /* xmm1 = ???? | ???? | ????? |  fabsf(dot(pPlane, pPoint)) + pPlane->d */

	"movss %%xmm1, -4(%%esp)\n"
	"flds -4(%%esp)\n"

	: "=t" (dummy)
	: "r" (pPlane), "r" (pPoint)
	);

}

float plane_Dot(const CPlane* pPlane, const CVec3D* pVec)
{
    float dummy;
	asm(
	"movups (%2), %%xmm1\n" /* xmm1 = ?    | V->z | V->x | V->y */
	"movups (%1), %%xmm0\n" /* xmm0 = P->d | P->c | P->b | P->a */
	"andps _SIMDx86_float_SSE_NO_W_MASK, %%xmm1\n" /* 0 | z | y | x */
	"mulps %%xmm0, %%xmm1\n"    /* xmm0 = 0 | P->c*V->z | P->b*V->y | P->a*V->x */


		"shufps $0xFF, %%xmm0, %%xmm0"
		"movhlps %%xmm1, %%xmm7\n"		/* xmm7 = ?   | ?   | 0   |  z's */
		"addss %%xmm1, %%xmm7\n"		/* xmm7 = ?   | ?   | 0   | z's + x's */
		"shufps $0x55, %%xmm1, %%xmm1\n"/* xmm1 = y's | y's | y's | y's */
		"addss %%xmm0, %%xmm7\n"        /* xmm7  = ?  | ?   | ?   | z's + x's + d*/
		"addss %%xmm1, %%xmm7\n"		/* xmm7 = ?   | ?   | ?   | x's + y's + z's + d */
		"movss %%xmm7, -4(%%esp)\n"
		"flds -4(%%esp)\n"


	: "=t" (dummy)
	: "r" (pPlane), "r" (pVec)
	);


}

float plane_Dot4(const CPlane* pPlane, const CVec4D* pVec4)
{
    float dummy;
	asm(
	"movups (%1), %%xmm0\n" /* xmm0 = d | c | b | a */
	"movups (%2), %%xmm1\n" /* xmm1 = w | z | y | z */
	"mulps %%xmm0, %%xmm1\n"/* xmm1 = dw | cz | by | az */


		"movhlps %%xmm1, %%xmm7\n"		/* xmm7 = ? | ? | dw | cz */
		"addps %%xmm1, %%xmm7\n"        /* xmm7 = ? | ? | dw+by | cz+ax */
		"movaps %%xmm7, %%xmm1\n"       /* xmm1 = ? | ? | dw+by | cz+ax */
		"shufps $0x55, %%xmm1, %%xmm1\n"/* xmm1 = ? | ? | ? | dw+by */
		"addss %%xmm1, %%xmm7\n"        /* xmm7  = ?  | ?   | ?   | dw+cz+by+ax */
		"movss %%xmm7, -4(%%esp)\n"
		"flds -4(%%esp)\n"


	: "=t" (dummy)
	: "r" (pPlane), "r" (pVec4)
	);


}

float plane_DotNormal(const CPlane* pPlane, const CVec3D* pVec)
{
	float dummy;
	asm(
	"movups (%1), %%xmm0\n" /* xmm0 = P->d | P->c | P->b | P->a */
	"movups (%2), %%xmm1\n" /* xmm1 = ?    | V->z | V->b | V->a */
	"andps _SIMDx86_float_SSE_NO_W_MASK, %%xmm0\n" /* 0 | P1.c | P1.b | P1.a */

	"mulps %%xmm1, %%xmm0\n"    /* xmm0 = 0 | P->c*V->z | P->b*V->y | P->a*V->x */

		"movhlps %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   |  z's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
		"shufps $0x55, %%xmm0, %%xmm0\n"/* xmm0 = y's | y's | y's | y's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */

		"movss %%xmm1, -4(%%esp)\n"
		"flds -4(%%esp)\n"

	: "=t" (dummy)
	: "r" (pPlane), "r" (pVec)
	);


}

float plane_DotPlane(const CPlane* pPlane1, const CPlane* pPlane2)
{
	float dummy;
	asm(
	"movups (%1), %%xmm0\n" /* xmm0 = P1.d | P1.c | P1.b | P1.a */
	"movups (%2), %%xmm1\n" /* xmm1 = P2.d | P2.c | P2.b | P2.a */
	"andps _SIMDx86_float_SSE_NO_W_MASK, %%xmm0\n" /* 0 | P1.c | P1.b | P1.a */

	"mulps %%xmm1, %%xmm0\n"    /* xmm0 = 0 | P1.c*P2.c | P1.b*P2.b | P1.a*P2.a */

	"movhlps %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   |  z's */
	"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
	"shufps $0x55, %%xmm0, %%xmm0\n"/* xmm0 = y's | y's | y's | y's */
	"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */

	"movss %%xmm1, -4(%%esp)\n"
	"flds -4(%%esp)\n"
	: "=t" (dummy)
	: "r" (pPlane1), "r" (pPlane2)
	);


}

void plane_Normalize(CPlane* pOut)
{

	asm(
	"movups (%0), %%xmm0\n"								/* xmm0 = d | c | b | a */
	"movaps %%xmm0, %%xmm1\n"							/* xmm1 = d | c | b | a */
	"andps _SIMDx86_float_SSE_NO_W_MASK, %%xmm0\n"		/* xmm0 = 0 | c | b | a */

	/* Perform dot product */
	"mulps %%xmm0, %%xmm0\n"							/* xmm0 = 0 | c*c | b*b | a*a */

	"movhlps %%xmm0, %%xmm3\n"							/* xmm3 = ? | ? | 0 | c*c */
	"addss %%xmm0, %%xmm3\n"							/* xmm3 = ? | ? | ? | a*a + c*c */
	"shufps $0x55, %%xmm0, %%xmm0\n"					/* xmm0 = b*b | b*b | b*b | b*b */
	"addss %%xmm3, %%xmm0\n"							/* xmm0 = ? | ? | ? | a*b + b*b + c*c */


	/* Divide it all by the square root */
	#if defined(HIPREC)
	"sqrtss %%xmm0, %%xmm0\n"	/* xmm0 = sqrtf(a*a+b*b+c*c) */
	"shufps %%xmm0, %%xmm0\n"	/* xmm0 = mag | mag | mag | mag */
	"divps %%xmm0, %%xmm1\n"	/* xmm1 = d/mag | c/mag | b/mag | a/mag */
	#else
	"rsqrtss %%xmm0, %%xmm0\n"	/* xmm0 = ? | ? | ? | 1.0f / sqrtf(a*a+b*b+c*c) */
	"mulps %%xmm0, %%xmm1\n"	/* xmm1 = d*invmag | c*invmag | b*invmag | a*invmag */
	#endif

	"movups %%xmm1, (%0)\n"
	:
	: "r" (pOut)
	);
	return;



}

void plane_NormalizeOf(CPlane* pOut, CPlane* pIn)
{

	asm(
	"movups (%0), %%xmm0\n"								/* xmm0 = d | c | b | a */
	"movaps %%xmm0, %%xmm1\n"							/* xmm1 = d | c | b | a */
	"andps _SIMDx86_float_SSE_NO_W_MASK, %%xmm0\n"		/* xmm0 = 0 | c | b | a */

	/* Perform dot product */
	"mulps %%xmm0, %%xmm0\n"							/* xmm0 = 0 | c*c | b*b | a*a */
	"movhlps %%xmm0, %%xmm3\n"							/* xmm3 = ? | ? | 0 | c*c */
	"addss %%xmm0, %%xmm3\n"							/* xmm3 = ? | ? | ? | a*a + c*c */
	"shufps $0x55, %%xmm0, %%xmm0\n"					/* xmm0 = b*b | b*b | b*b | b*b */
	"addss %%xmm3, %%xmm0\n"							/* xmm0 = ? | ? | ? | a*b + b*b + c*c */

	/* Divide it all by the square root */
	#if defined(HIPREC)
	"sqrtss %%xmm0, %%xmm0\n"	/* xmm0 = sqrtf(a*a+b*b+c*c) */
	"shufps %%xmm0, %%xmm0\n"	/* xmm0 = mag | mag | mag | mag */
	"divps  %%xmm0,%%xmm1\n"	/* xmm1 = d/mag | c/mag | b/mag | a/mag */
	#else
	"rsqrtss %%xmm0, %%xmm0\n"	/* xmm0 = ? | ? | ? | 1.0f / sqrtf(a*a+b*b+c*c) */
	"mulps  %%xmm0,%%xmm1\n"	/* xmm1 = d*invmag | c*invmag | b*invmag | a*invmag */
	#endif

	"movups %%xmm1, (%1)\n"
	:
	: "r" (pOut), "r" (pIn)
	);
	return;

}



#endif // -masm=AT&T





} //end MATH
} //end SGF
