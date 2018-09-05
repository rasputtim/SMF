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

#ifndef _SMF__MATH_SIMD_SSE_H__
#define _SMF__MATH_SIMD_SSE_H__
#include "../SMF_Config.h"
#include "SMF_Math.h"
#include "SMF_SimdMMX.h"
#include "SMF_Quaternion.h"
#include "../geometry/SMF_DrawVert.h"

namespace SMF {
namespace MATH{
#if defined(MACOS_X) && defined(__i686__)

#include <xmmintrin.h>


//#define SHUFFLEPS( x, y, z, w )		(( (x) & 3 ) << 6 | ( (y) & 3 ) << 4 | ( (z) & 3 ) << 2 | ( (w) & 3 ))
//#define R_SHUFFLEPS( x, y, z, w )	(( (w) & 3 ) << 6 | ( (z) & 3 ) << 4 | ( (y) & 3 ) << 2 | ( (x) & 3 ))
#define SHUFFLEPS( x, y, z, w )		(( (x) & 3 ) << 6 | ( (y) & 3 ) << 4 | ( (z) & 3 ) << 2 | ( (w) & 3 ))
//#define R_SHUFFLEPS( x, y, z, w )	(( (w) & 3 ) << 6 | ( (z) & 3 ) << 4 | ( (y) & 3 ) << 2 | ( (x) & 3 ))
#define R_SHUFFLEPS( x, y, z, w )	 w & 3  << 6 |  z & 3  << 4 |  y & 3  << 2 |  x & 3

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
CSIMD_SSE::dot

  dst[i] = constant.Normal() * src[i].xyz + constant[3];
============
*/
void VPCALL CSIMD_SSE::dot( float *dst, const CPlane &constant, const CVertex *src, const int count ) {
	// 0,  1,  2
	// 3,  4,  5
	// 6,  7,  8
	// 9, 10, 11

	/*
		mov			eax, count
		mov			edi, constant
		mov			edx, eax
		mov			esi, src
		mov			ecx, dst
	*/
	sf_m128 xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7; 	// Declare 8 xmm registers.
	int count_l4 = count;                                   // count_l4 = eax
	int count_l1 = count;                                   // count_l1 = edx
	char *constant_p = (char *)&constant;                   // constant_p = edi
	char *src_p = (char *) src;                             // src_p = esi
	char *dst_p = (char *) dst;                             // dst_p = ecx

	SMF_ASSERT( sizeof( CVertex ) == DRAWVERT_SIZE );
	SMF_ASSERT( (int)&((CVertex *)0)->xyz == DRAWVERT_XYZ_OFFSET );

	/*
		and			eax, ~3
		movss		xmm4, [edi+0]
		shufps		xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 )
		movss		xmm5, [edi+4]
		shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )
		movss		xmm6, [edi+8]
		shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
		movss		xmm7, [edi+12]
		shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
	*/
	count_l4 = count_l4 & ~3;
	xmm4 = _mm_load_ss((float *) (constant_p));
	xmm4 = _mm_shuffle_ps(xmm4, xmm4, R_SHUFFLEPS( 0, 0, 0, 0 ));
	xmm5 = _mm_load_ss((float *) (constant_p + 4));
	xmm5 = _mm_shuffle_ps(xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 ));
	xmm6 = _mm_load_ss((float *) (constant_p + 8));
	xmm6 = _mm_shuffle_ps(xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 ));
	xmm7 = _mm_load_ss((float *) (constant_p + 12));
 	xmm7 = _mm_shuffle_ps(xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 ));

	/*
		jz			startVert1
	*/
	if(count_l4 != 0) {
	/*
		imul		eax, DRAWVERT_SIZE
		add			esi, eax
		neg			eax
	*/
		count_l4 = count_l4 * DRAWVERT_SIZE;
		src_p = src_p + count_l4;
		count_l4 = -count_l4;
	/*
	loopVert4:
	*/
		do {
	/*
		movss		xmm0, [esi+eax+1*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]	//  3,  X,  X,  X
		movss		xmm2, [esi+eax+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+8]	//  2,  X,  X,  X
		movhps		xmm0, [esi+eax+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]	//  3,  X,  0,  1
		movaps		xmm1, xmm0												//  3,  X,  0,  1
	*/
			xmm0 = _mm_load_ss((float *) (src_p+count_l4+1*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0));        // 3,  X,  X,  X
			xmm2 = _mm_load_ss((float *) (src_p+count_l4+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+8));        // 2,  X,  X,  X
			xmm0 = _mm_loadh_pi(xmm0, (__m64 *) (src_p+count_l4+0*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0)); // 3,  X,  0,  1
			xmm1 = xmm0;							                                                    // 3,  X,  0,  1

	/*
		movlps		xmm1, [esi+eax+1*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+4]	//  4,  5,  0,  1
		shufps		xmm2, xmm1, R_SHUFFLEPS( 0, 1, 0, 1 )					//  2,  X,  4,  5
	*/
			xmm1 = _mm_loadl_pi(xmm1, (__m64 *) (src_p+count_l4+1*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+4)); // 4,  5,  0,  1
			xmm2 = _mm_shuffle_ps(xmm2, xmm1, R_SHUFFLEPS( 0, 1, 0, 1 ));                               // 2,  X,  4,  5

	/*
		movss		xmm3, [esi+eax+3*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]	//  9,  X,  X,  X
		movhps		xmm3, [esi+eax+2*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0]	//  9,  X,  6,  7
		shufps		xmm0, xmm3, R_SHUFFLEPS( 2, 0, 2, 0 )					//  0,  3,  6,  9
	*/
			xmm3 = _mm_load_ss((float *) (src_p+count_l4+3*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0));        // 9,  X,  X,  X
			xmm3 = _mm_loadh_pi(xmm3, (__m64 *) (src_p+count_l4+2*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+0)); // 9,  X,  6,  7
			xmm0 = _mm_shuffle_ps(xmm0, xmm3, R_SHUFFLEPS( 2, 0, 2, 0 ));                               // 0,  3,  6,  9
	/*
		movlps		xmm3, [esi+eax+3*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+4]	// 10, 11,  6,  7
		shufps		xmm1, xmm3, R_SHUFFLEPS( 3, 0, 3, 0 )					//  1,  4,  7, 10
	*/
			xmm3 = _mm_loadl_pi(xmm3, (__m64 *)(src_p+count_l4+3*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+4));  // 10, 11, 6,  7
			xmm1 = _mm_shuffle_ps(xmm1, xmm3, R_SHUFFLEPS( 3, 0, 3, 0 ));                               // 1,  4,  7,  10
	/*
		movhps		xmm3, [esi+eax+2*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+8]	// 10, 11,  8,  X
		shufps		xmm2, xmm3, R_SHUFFLEPS( 0, 3, 2, 1 )					//  2,  5,  8, 11
	*/
			xmm3 = _mm_loadh_pi(xmm3, (__m64 *)(src_p+count_l4+2*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET+8));  // 10, 11, 8,  X
			xmm2 = _mm_shuffle_ps(xmm2, xmm3, R_SHUFFLEPS( 0, 3, 2, 1 ));                               // 2,  5,  8,  11

	/*
		add			ecx, 16
		add			eax, 4*DRAWVERT_SIZE
	*/
			dst_p = dst_p + 16;
			count_l4 = count_l4 + 4*DRAWVERT_SIZE;

	/*
		mulps		xmm0, xmm4
		mulps		xmm1, xmm5
		mulps		xmm2, xmm6
		addps		xmm0, xmm7
		addps		xmm0, xmm1
		addps		xmm0, xmm2
	*/
			xmm0 = _mm_mul_ps(xmm0, xmm4);
			xmm1 = _mm_mul_ps(xmm1, xmm5);
			xmm2 = _mm_mul_ps(xmm2, xmm6);
			xmm0 = _mm_add_ps(xmm0, xmm7);
			xmm0 = _mm_add_ps(xmm0, xmm1);
			xmm0 = _mm_add_ps(xmm0, xmm2);

	/*
		movlps		[ecx-16+0], xmm0
		movhps		[ecx-16+8], xmm0
		jl			loopVert4
	*/
			_mm_storel_pi((__m64 *) (dst_p-16+0), xmm0);
			_mm_storeh_pi((__m64 *) (dst_p-16+8), xmm0);
		} while(count_l4 < 0);
	}

	/*
	startVert1:
		and			edx, 3
		jz			done
	*/
	count_l1 = count_l1 & 3;
	if(count_l1 != 0) {
	/*
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
	*/
		do {
			xmm0 = _mm_load_ss((float *) (src_p+count_l4+DRAWVERT_XYZ_OFFSET+0));
			xmm1 = _mm_load_ss((float *) (src_p+count_l4+DRAWVERT_XYZ_OFFSET+4));
			xmm2 = _mm_load_ss((float *) (src_p+count_l4+DRAWVERT_XYZ_OFFSET+8));
			xmm0 = _mm_mul_ss(xmm0, xmm4);
			xmm1 = _mm_mul_ss(xmm1, xmm5);
			xmm2 = _mm_mul_ss(xmm2, xmm6);
			xmm0 = _mm_add_ss(xmm0, xmm7);
			dst_p = dst_p + 4;
			xmm0 = _mm_add_ss(xmm0, xmm1);
			count_l4 = count_l4 + DRAWVERT_SIZE;
			xmm0 = _mm_add_ss(xmm0, xmm2);
			count_l1 = count_l1 - 1;
			_mm_store_ss((float *) (dst_p-4), xmm0);
		} while( count_l1 != 0);
	}
	/*
		done:
	*/
}

/*
============
CSIMD_SSE::minMax
============
*/
void VPCALL CSIMD_SSE::minMax( CVec3D &min, CVec3D &max, const CVertex *src, const int *indexes, const int count ) {

	SMF_ASSERT( sizeof( CVertex ) == DRAWVERT_SIZE );
	SMF_ASSERT( (int)&((CVertex *)0)->xyz == DRAWVERT_XYZ_OFFSET );

	sf_m128 xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7;
	char *indexes_p;
	char *src_p;
	int count_l;
	int edx;
	char *min_p;
	char *max_p;

	/*
		movss		xmm0, CMath::INFINITY_FLOAT
		xorps		xmm1, xmm1
		shufps		xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 )
		subps		xmm1, xmm0
		movaps		xmm2, xmm0
		movaps		xmm3, xmm1
	*/
		xmm0 = _mm_load_ss(&CMath::INFINITY_FLOAT);
		// To satisfy the compiler use xmm0 instead.
		xmm1 = _mm_xor_ps(xmm0, xmm0);
		xmm0 = _mm_shuffle_ps(xmm0, xmm0, R_SHUFFLEPS( 0, 0, 0, 0 ));
		xmm1 = _mm_sub_ps(xmm1, xmm0);
		xmm2 = xmm0;
		xmm3 = xmm1;

	/*
		mov			edi, indexes
		mov			esi, src
		mov			eax, count
		and			eax, ~3
		jz			done4
	*/
		indexes_p = (char *) indexes;
		src_p = (char *) src;
		count_l = count;
		count_l = count_l & ~3;
		if(count_l != 0) {
	/*
		shl			eax, 2
		add			edi, eax
		neg			eax
	*/
			count_l = count_l << 2;
			indexes_p = indexes_p + count_l;
			count_l = -count_l;
	/*
	loop4:
//		prefetchnta	[edi+128]
//		prefetchnta	[esi+4*DRAWVERT_SIZE+DRAWVERT_XYZ_OFFSET]
	*/
		do {
	/*
		mov			edx, [edi+eax+0]
		imul		edx, DRAWVERT_SIZE
		movss		xmm4, [esi+edx+DRAWVERT_XYZ_OFFSET+8]
		movhps		xmm4, [esi+edx+DRAWVERT_XYZ_OFFSET+0]
		minps		xmm0, xmm4
		maxps		xmm1, xmm4
	*/
			edx = *((int*)(indexes_p+count_l+0));
			edx = edx * DRAWVERT_SIZE;
			xmm4 = _mm_load_ss((float *) (src_p+edx+DRAWVERT_XYZ_OFFSET+8));
			xmm4 = _mm_loadh_pi(xmm4, (__m64 *) (src_p+edx+DRAWVERT_XYZ_OFFSET+0) );
			xmm0 = _mm_min_ps(xmm0, xmm4);
			xmm1 = _mm_max_ps(xmm1, xmm4);

	/*
		mov			edx, [edi+eax+4]
		imul		edx, DRAWVERT_SIZE
		movss		xmm5, [esi+edx+DRAWVERT_XYZ_OFFSET+0]
		movhps		xmm5, [esi+edx+DRAWVERT_XYZ_OFFSET+4]
		minps		xmm2, xmm5
		maxps		xmm3, xmm5
	*/
			edx = *((int*)(indexes_p+count_l+4));
			edx = edx * DRAWVERT_SIZE;
			xmm5 = _mm_load_ss((float *) (src_p+edx+DRAWVERT_XYZ_OFFSET+0));
			xmm5 = _mm_loadh_pi(xmm5, (__m64 *) (src_p+edx+DRAWVERT_XYZ_OFFSET+4) );
			xmm2 = _mm_min_ps(xmm2, xmm5);
			xmm3 = _mm_max_ps(xmm3, xmm5);

	/*
		mov			edx, [edi+eax+8]
		imul		edx, DRAWVERT_SIZE
		movss		xmm6, [esi+edx+DRAWVERT_XYZ_OFFSET+8]
		movhps		xmm6, [esi+edx+DRAWVERT_XYZ_OFFSET+0]
		minps		xmm0, xmm6
		maxps		xmm1, xmm6
	*/
			edx = *((int*)(indexes_p+count_l+8));
			edx = edx * DRAWVERT_SIZE;
			xmm6 = _mm_load_ss((float *) (src_p+edx+DRAWVERT_XYZ_OFFSET+8));
			xmm6 = _mm_loadh_pi(xmm6, (__m64 *) (src_p+edx+DRAWVERT_XYZ_OFFSET+0) );
			xmm0 = _mm_min_ps(xmm0, xmm6);
			xmm1 = _mm_max_ps(xmm1, xmm6);

	/*
		mov			edx, [edi+eax+12]
		imul		edx, DRAWVERT_SIZE
		movss		xmm7, [esi+edx+DRAWVERT_XYZ_OFFSET+0]
		movhps		xmm7, [esi+edx+DRAWVERT_XYZ_OFFSET+4]
		minps		xmm2, xmm7
		maxps		xmm3, xmm7
	*/
			edx = *((int*)(indexes_p+count_l+12));
			edx = edx * DRAWVERT_SIZE;
			xmm7 = _mm_load_ss((float *) (src_p+edx+DRAWVERT_XYZ_OFFSET+0));
			xmm7 = _mm_loadh_pi(xmm7, (__m64 *) (src_p+edx+DRAWVERT_XYZ_OFFSET+4) );
			xmm2 = _mm_min_ps(xmm2, xmm7);
			xmm3 = _mm_max_ps(xmm3, xmm7);

	/*
		add			eax, 4*4
		jl			loop4
	*/
			count_l = count_l + 4*4;
		} while (count_l < 0);
	}
	/*
	done4:
		mov			eax, count
		and			eax, 3
		jz			done1
	*/
	count_l = count;
	count_l = count_l & 3;
	if(count_l != 0) {
	/*
		shl			eax, 2
		add			edi, eax
		neg			eax
	*/
		count_l = count_l << 2;
		indexes_p = indexes_p + count_l;
		count_l = -count_l;
	/*
	loop1:
	*/
		do{
	/*
		mov			edx, [edi+eax+0]
		imul		edx, DRAWVERT_SIZE;
		movss		xmm4, [esi+edx+DRAWVERT_XYZ_OFFSET+8]
		movhps		xmm4, [esi+edx+DRAWVERT_XYZ_OFFSET+0]
		minps		xmm0, xmm4
		maxps		xmm1, xmm4
	*/
			edx = *((int*)(indexes_p+count_l+0));
			edx = edx * DRAWVERT_SIZE;
			xmm4 = _mm_load_ss((float *) (src_p+edx+DRAWVERT_XYZ_OFFSET+8));
			xmm4 = _mm_loadh_pi(xmm4, (__m64 *) (src_p+edx+DRAWVERT_XYZ_OFFSET+0) );
			xmm0 = _mm_min_ps(xmm0, xmm4);
			xmm1 = _mm_max_ps(xmm1, xmm4);

	/*
		add			eax, 4
		jl			loop1
	*/
			count_l = count_l + 4;
		} while (count_l < 0);

	}

	/*
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
	*/
	xmm2 = _mm_shuffle_ps(xmm2, xmm2, R_SHUFFLEPS( 3, 1, 0, 2 ));
	xmm3 = _mm_shuffle_ps(xmm3, xmm3, R_SHUFFLEPS( 3, 1, 0, 2 ));
	xmm0 = _mm_min_ps(xmm0, xmm2);
	xmm1 = _mm_max_ps(xmm1, xmm3);
	min_p = (char *) &min;
	_mm_storeh_pi((__m64 *)(min_p), xmm0);
	_mm_store_ss((float *)(min_p+8), xmm0);
	max_p = (char *) &max;
	_mm_storeh_pi((__m64 *)(max_p), xmm1);
	_mm_store_ss((float *)(max_p+8), xmm1);
}

/*
============
CSIMD_SSE::dot

  dst[i] = constant * src[i].Normal() + src[i][3];
============
*/
void VPCALL CSIMD_SSE::dot( float *dst, const CVec3D &constant, const CPlane *src, const int count ) {
	int count_l4;
	int count_l1;
	char *constant_p;
	char *src_p;
	char *dst_p;
	sf_m128 xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7;

	/*
		mov			eax, count
		mov			edi, constant
		mov			edx, eax
		mov			esi, src
		mov			ecx, dst
		and			eax, ~3
	*/
	count_l4 = count;
	constant_p = (char *) &constant;
	count_l1 = count_l4;
	src_p = (char *) src;
	dst_p = (char *) dst;
	count_l4 = count_l4 & ~3;

	/*
		movss		xmm5, [edi+0]
		shufps		xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 )
		movss		xmm6, [edi+4]
		shufps		xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 )
		movss		xmm7, [edi+8]
		shufps		xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 )
	*/
	xmm5 = _mm_load_ss((float *) (constant_p+0));
	xmm5 = _mm_shuffle_ps(xmm5, xmm5, R_SHUFFLEPS( 0, 0, 0, 0 ));
	xmm6 = _mm_load_ss((float *) (constant_p+4));
	xmm6 = _mm_shuffle_ps(xmm6, xmm6, R_SHUFFLEPS( 0, 0, 0, 0 ));
	xmm7 = _mm_load_ss((float *) (constant_p+8));
	xmm7 = _mm_shuffle_ps(xmm7, xmm7, R_SHUFFLEPS( 0, 0, 0, 0 ));

	/*
		jz			startVert1
	*/
	if (count != 0) {
	/*
		imul		eax, 16
		add			esi, eax
		neg			eax
	*/
		count_l4 = count_l4 * 16;
		src_p = src_p + count_l4;
		count_l4 = -count_l4;
	/*
	loopVert4:
	*/
		do {
	/*
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
	*/
			xmm1 = _mm_loadl_pi(xmm1, (__m64 *)(src_p+count_l4+ 0));
			xmm3 = _mm_loadl_pi(xmm3, (__m64 *)(src_p+count_l4+ 8));
			xmm1 = _mm_loadh_pi(xmm1, (__m64 *)(src_p+count_l4+16));
			xmm3 = _mm_loadh_pi(xmm3, (__m64 *)(src_p+count_l4+24));
			xmm2 = _mm_loadl_pi(xmm2, (__m64 *)(src_p+count_l4+32));
			xmm4 = _mm_loadl_pi(xmm4, (__m64 *)(src_p+count_l4+40));
			xmm2 = _mm_loadh_pi(xmm2, (__m64 *)(src_p+count_l4+48));
			xmm4 = _mm_loadh_pi(xmm4, (__m64 *)(src_p+count_l4+56));

			xmm0 = xmm1;
			xmm0 = _mm_shuffle_ps(xmm0, xmm2, R_SHUFFLEPS( 0, 2, 0, 2 ));
			xmm1 = _mm_shuffle_ps(xmm1, xmm2, R_SHUFFLEPS( 1, 3, 1, 3 ));
			xmm2 = xmm3;
			xmm2 = _mm_shuffle_ps(xmm2, xmm4, R_SHUFFLEPS( 0, 2, 0, 2 ));
			xmm3 = _mm_shuffle_ps(xmm3, xmm4, R_SHUFFLEPS( 1, 3, 1, 3 ));

	/*
		add			ecx, 16
		add			eax, 4*16
	*/
			dst_p = dst_p + 16;
			count_l4 = count_l4 + 4*16;

	/*
		mulps		xmm0, xmm5
		mulps		xmm1, xmm6
		mulps		xmm2, xmm7
		addps		xmm0, xmm3
		addps		xmm0, xmm1
		addps		xmm0, xmm2
	*/
			xmm0 = _mm_mul_ps(xmm0, xmm5);
			xmm1 = _mm_mul_ps(xmm1, xmm6);
			xmm2 = _mm_mul_ps(xmm2, xmm7);
			xmm0 = _mm_add_ps(xmm0, xmm3);
			xmm0 = _mm_add_ps(xmm0, xmm1);
			xmm0 = _mm_add_ps(xmm0, xmm2);

	/*
		movlps		[ecx-16+0], xmm0
		movhps		[ecx-16+8], xmm0
		jl			loopVert4
	*/
			_mm_storel_pi((__m64 *) (dst_p-16+0), xmm0);
			_mm_storeh_pi((__m64 *) (dst_p-16+8), xmm0);
		} while (count_l4 < 0);
	}

	/*
	startVert1:
		and			edx, 3
		jz			done
	*/
	count_l1 = count_l1 & 3;

	if(count_l1 != 0) {
	/*
	loopVert1:
	*/
		do {
	/*
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
	*/
			xmm0 = _mm_load_ss((float *) (src_p+count_l4+ 0));
			xmm1 = _mm_load_ss((float *) (src_p+count_l4+ 4));
			xmm2 = _mm_load_ss((float *) (src_p+count_l4+ 8));
			xmm3 = _mm_load_ss((float *) (src_p+count_l4+12));

			xmm0 = _mm_mul_ss(xmm0, xmm5);
			xmm1 = _mm_mul_ss(xmm1, xmm6);
			xmm2 = _mm_mul_ss(xmm2, xmm7);

			xmm0 = _mm_add_ss(xmm0, xmm3);
			dst_p = dst_p + 4;
			xmm0 = _mm_add_ss(xmm0, xmm1);
			count_l4 = count_l4 + 16;
			xmm0 = _mm_add_ss(xmm0, xmm2);
			count_l1 = count_l1 - 1;
			_mm_store_ss((float *) (dst_p-4), xmm0);
		} while (count_l1 != 0);
	}
	/*
	done:
	*/
}

#elif defined(_MSC_VER)

#include <xmmintrin.h>
//http://www.jaist.ac.jp/iscenter-new/mpc/altix/altixdata/opt/intel/vtune/doc/users_guide/mergedProjects/analyzer_ec/mergedProjects/reference_olh/mergedProjects/instructions/instruct32_hh/vc293.htm
/**
shufps oper1, oper2, control
xyzw são valores que variam de 0 a 3 (dec) que determinal qual vai ser o byte de controle da instrução

x determina qual parte de oper2 será colocada na quarta doubleword (32 bits) do oper1

CASE (x) OF
0: oper1[127-96] oper2[31-0];
1: oper1[127-96] oper2[63-32];
2: oper1[127-96] oper2[95-64];
3: oper1[127-96] oper2[127-96];
ESAC

y determina qual parte de oper2 será colocada na terceira doubleword (32 bits) do oper1
CASE (y) OF
0: oper1[95-64] oper2[31-0];
1: oper1[95-64] oper2[63-32];
2: oper1[95-64] oper2[95-64];
3: oper1[95-64] oper2[127-96];
ESAC;
z determina qual parte de oper1 será colocada na segunda doubleword (32 bits) do oper1
CASE (z) OF
0: oper1[63-32] oper1[31-0];
1: oper1[63-32] oper1[63-32];
2: oper1[63-32] oper1[95-64];
3: oper1[63-32] oper1[127-96];
ESAC;
w determina qual parte de oper1 será colocada na primeira doubleword (32 bits) do oper1
CASE (w) OF
0: oper1[31-0] oper1[31-0];
1: oper1[31-0] oper1[63-32];
2: oper1[31-0] oper1[95-64];
3: oper1[31-0] oper1[127-96];
ESAC;
**/
#define SHUFFLEPS( x, y, z, w )		(( (x) & 3 ) << 6 | ( (y) & 3 ) << 4 | ( (z) & 3 ) << 2 | ( (w) & 3 ))
/**
shufps oper1, oper2, control
w determina qual parte de oper2 será colocada na quarta doubleword (32 bits) do oper1
CASE (w) OF
0: oper1[127-96] oper2[31-0];
1: oper1[127-96] oper2[63-32];
2: oper1[127-96] oper2[95-64];
3: oper1[127-96] oper2[127-96];
ESAC
z determina qual parte de oper2 será colocada na terceira doubleword (32 bits) do oper1
CASE (z) OF
0: oper1[95-64] oper2[31-0];
1: oper1[95-64] oper2[63-32];
2: oper1[95-64] oper2[95-64];
3: oper1[95-64] oper2[127-96];
ESAC;
y determina qual parte de oper1 será colocada na segunda doubleword (32 bits) do oper1
CASE (y) OF
0: oper1[63-32] oper1[31-0];
1: oper1[63-32] oper1[63-32];
2: oper1[63-32] oper1[95-64];
3: oper1[63-32] oper1[127-96];
ESAC;
x determina qual parte de oper1 será colocada na primeira doubleword (32 bits) do oper1
CASE (x) OF
0: oper1[31-0] oper1[31-0];
1: oper1[31-0] oper1[63-32];
2: oper1[31-0] oper1[95-64];
3: oper1[31-0] oper1[127-96];
ESAC;
**/

// transpose a 4x4 matrix loaded into 4 xmm registers (reg4 is temporary)
#define TRANSPOSE_4x4( reg0, reg1, reg2, reg3, reg4 )											\
	__asm	movaps		reg4, reg2								/* reg4 =  8,  9, 10, 11 */		\
	__asm	unpcklps	reg2, reg3								/* reg2 =  8, 12,  9, 13 */		\
	__asm	unpckhps	reg4, reg3								/* reg4 = 10, 14, 11, 15 */		\
	__asm	movaps		reg3, reg0								/* reg3 =  0,  1,  2,  3 */		\
	__asm	unpcklps	reg0, reg1								/* reg0 =  0,  4,  1,  5 */		\
	__asm	unpckhps	reg3, reg1								/* reg3 =  2,  6,  3,  7 */		\
	__asm	movaps		reg1, reg0								/* reg1 =  0,  4,  1,  5 */		\
	__asm	shufps		reg0, reg2, R_SHUFFLEPS( 0, 1, 0, 1 )	/* reg0 =  0,  4,  8, 12 */		\
	__asm	shufps		reg1, reg2, R_SHUFFLEPS( 2, 3, 2, 3 )	/* reg1 =  1,  5,  9, 13 */		\
	__asm	movaps		reg2, reg3								/* reg2 =  2,  6,  3,  7 */		\
	__asm	shufps		reg2, reg4, R_SHUFFLEPS( 0, 1, 0, 1 )	/* reg2 =  2,  6, 10, 14 */		\
	__asm	shufps		reg3, reg4, R_SHUFFLEPS( 2, 3, 2, 3 )	/* reg3 =  3,  7, 11, 15 */

// transpose a 4x4 matrix from memory into 4 xmm registers (reg4 is temporary)
#define TRANPOSE_4x4_FROM_MEMORY( address, reg0, reg1, reg2, reg3, reg4 )						\
	__asm	movlps		reg1, [address+ 0]						/* reg1 =  0,  1,  X,  X */		\
	__asm	movlps		reg3, [address+ 8]						/* reg3 =  2,  3,  X,  X */		\
	__asm	movhps		reg1, [address+16]						/* reg1 =  0,  1,  4,  5 */		\
	__asm	movhps		reg3, [address+24]						/* reg3 =  2,  3,  6,  7 */		\
	__asm	movlps		reg2, [address+32]						/* reg2 =  8,  9,  X,  X */		\
	__asm	movlps		reg4, [address+40]						/* reg4 = 10, 11,  X,  X */		\
	__asm	movhps		reg2, [address+48]						/* reg2 =  8,  9, 12, 13 */		\
	__asm	movhps		reg4, [address+56]						/* reg4 = 10, 11, 14, 15 */		\
	__asm	movaps		reg0, reg1								/* reg0 =  0,  1,  4,  5 */		\
	__asm	shufps		reg0, reg2, R_SHUFFLEPS( 0, 2, 0, 2 )	/* reg0 =  0,  4,  8, 12 */		\
	__asm	shufps		reg1, reg2, R_SHUFFLEPS( 1, 3, 1, 3 )	/* reg1 =  1,  5,  9, 13 */		\
	__asm	movaps		reg2, reg3								/* reg2 =  2,  3,  6,  7 */		\
	__asm	shufps		reg2, reg4, R_SHUFFLEPS( 0, 2, 0, 2 )	/* reg2 =  2,  6, 10, 14 */		\
	__asm	shufps		reg3, reg4, R_SHUFFLEPS( 1, 3, 1, 3 )	/* reg3 =  3,  7, 11, 15 */

// transpose a 4x4 matrix to memory from 4 xmm registers (reg4 is temporary)
#define TRANPOSE_4x4_TO_MEMORY( address, reg0, reg1, reg2, reg3, reg4 )							\
	__asm	movaps		reg4, reg0								/* reg4 =  0,  4,  8, 12 */		\
	__asm	unpcklps	reg0, reg1								/* reg0 =  0,  1,  4,  5 */		\
	__asm	unpckhps	reg4, reg1								/* reg4 =  8,  9, 12, 13 */		\
	__asm	movaps		reg1, reg2								/* reg1 =  2,  6, 10, 14 */		\
	__asm	unpcklps	reg2, reg3								/* reg2 =  2,  3,  6,  7 */		\
	__asm	unpckhps	reg1, reg3								/* reg1 = 10, 11, 14, 15 */		\
	__asm	movlps		[address+ 0], reg0						/* mem0 =  0,  1,  X,  X */		\
	__asm	movlps		[address+ 8], reg2						/* mem0 =  0,  1,  2,  3 */		\
	__asm	movhps		[address+16], reg0						/* mem1 =  4,  5,  X,  X */		\
	__asm	movhps		[address+24], reg2						/* mem1 =  4,  5,  6,  7 */		\
	__asm	movlps		[address+32], reg4						/* mem2 =  8,  9,  X,  X */		\
	__asm	movlps		[address+40], reg1						/* mem2 =  8,  9, 10, 11 */		\
	__asm	movhps		[address+48], reg4						/* mem3 = 12, 13,  X,  X */		\
	__asm	movhps		[address+56], reg1						/* mem3 = 12, 13, 14, 15 */

// transpose a 4x3 matrix loaded into 3 xmm registers (reg3 is temporary)
#define TRANSPOSE_4x3( reg0, reg1, reg2, reg3 )													\
	__asm	movaps		reg3, reg2								/* reg3 =  8,  9, 10, 11 */		\
	__asm	shufps		reg3, reg1, R_SHUFFLEPS( 2, 3, 0, 1 )	/* reg3 = 10, 11,  4,  5 */		\
	__asm	shufps		reg2, reg0, R_SHUFFLEPS( 0, 1, 2, 3 )	/* reg2 =  8,  9,  2,  3 */		\
	__asm	shufps		reg1, reg0, R_SHUFFLEPS( 2, 3, 0, 1 )	/* reg1 =  6,  7,  0,  1 */		\
	__asm	movaps		reg0, reg1								/* reg0 =  6,  7,  0,  1 */		\
	__asm	shufps		reg0, reg2, R_SHUFFLEPS( 2, 0, 3, 1 )	/* reg0 =  0,  6,  3,  9 */		\
	__asm	shufps		reg1, reg3, R_SHUFFLEPS( 3, 1, 2, 0 )	/* reg1 =  1,  7,  4, 10 */		\
	__asm	shufps		reg2, reg3, R_SHUFFLEPS( 2, 0, 3, 1 )	/* reg2 =  2,  8,  5, 11 */

// transpose a 4x3 matrix from memory into 3 xmm registers (reg3 is temporary)
#define TRANSPOSE_4x3_FROM_MEMORY( address, reg0, reg1, reg2, reg3 )							\
	__asm	movlps		reg1, [address+ 0]						/* reg1 =  0,  1,  X,  X */		\
	__asm	movlps		reg2, [address+ 8]						/* reg2 =  2,  3,  X,  X */		\
	__asm	movlps		reg3, [address+16]						/* reg3 =  4,  5,  X,  X */		\
	__asm	movhps		reg1, [address+24]						/* reg1 =  0,  1,  6,  7 */		\
	__asm	movhps		reg2, [address+32]						/* reg2 =  2,  3,  8,  9 */		\
	__asm	movhps		reg3, [address+40]						/* reg3 =  4,  5, 10, 11 */		\
	__asm	movaps		reg0, reg1								/* reg0 =  0,  1,  6,  7 */		\
	__asm	shufps		reg0, reg2, R_SHUFFLEPS( 0, 2, 1, 3 )	/* reg0 =  0,  6,  3,  9 */		\
	__asm	shufps		reg1, reg3, R_SHUFFLEPS( 1, 3, 0, 2 )	/* reg1 =  1,  7,  4, 10 */		\
	__asm	shufps		reg2, reg3, R_SHUFFLEPS( 0, 2, 1, 3 )	/* reg2 =  2,  8,  5, 11 */

// transpose a 4x3 matrix to memory from 3 xmm registers (reg3 is temporary)
#define TRANSPOSE_4x3_TO_MEMORY( address, reg0, reg1, reg2, reg3 )								\
	__asm	movhlps		reg3, reg0 								/* reg3 =  3,  9,  X,  X */		\
	__asm	unpcklps	reg0, reg1								/* reg0 =  0,  1,  6,  7 */		\
	__asm	unpckhps	reg1, reg2								/* reg1 =  4,  5, 10, 11 */		\
	__asm	unpcklps	reg2, reg3								/* reg2 =  2,  3,  8,  9 */		\
	__asm	movlps		[address+ 0], reg0						/* mem0 =  0,  1,  X,  X */		\
	__asm	movlps		[address+ 8], reg2						/* mem0 =  0,  1,  2,  3 */		\
	__asm	movlps		[address+16], reg1						/* mem1 =  4,  5,  X,  X */		\
	__asm	movhps		[address+24], reg0						/* mem1 =  4,  5,  6,  7 */		\
	__asm	movhps		[address+32], reg2						/* mem2 =  8,  9,  X,  X */		\
	__asm	movhps		[address+40], reg1						/* mem2 =  8,  9, 10, 11 */


// with alignment
#define KFLOATINITS(   SRC0, COUNT, PRE, POST )				KFLOATINITDSS( SRC0,SRC0,SRC0,COUNT,PRE,POST )
#define KFLOATINITD(   DST, COUNT, PRE, POST )				KFLOATINITDSS( DST,DST,DST,COUNT,PRE,POST )
#define KFLOATINITDS(  DST, SRC0, COUNT, PRE, POST )		KFLOATINITDSS( DST,SRC0,SRC0,COUNT,PRE,POST )

#define KFLOATINITDSS( DST, SRC0, SRC1, COUNT, PRE, POST )\
	__asm	mov		ecx,DST								\
	__asm	shr		ecx,2								\
	__asm	mov		ebx,COUNT							\
	__asm	neg		ecx									\
	__asm	mov		edx,SRC0							\
	__asm	and		ecx,3								\
	__asm	mov		esi,SRC1							\
	__asm	sub		ebx,ecx								\
	__asm	jge		noUnderFlow							\
	__asm	xor		ebx,ebx								\
	__asm	mov		ecx,COUNT							\
	__asm	noUnderFlow:								\
	__asm	mov		PRE,ecx								\
	__asm	mov		eax,ebx								\
	__asm	mov		edi,DST								\
	__asm	and		eax,8-1								\
	__asm	mov		POST,eax							\
	__asm	and		ebx,0xfffffff8						\
	__asm	jle		done								\
	__asm	shl		ebx,2								\
	__asm	lea		ecx,[ecx*4+ebx]						\
	__asm	neg		ebx									\
	__asm	add		edx,ecx								\
	__asm	add		esi,ecx								\
	__asm	add		edi,ecx								\
	__asm	mov		eax,edx								\
	__asm	or		eax,esi

// without alignment (pre==0)
#define KFLOATINITS_NA(   SRC0, COUNT, PRE, POST )				KFLOATINITDSS_NA( SRC0,SRC0,SRC0,COUNT,PRE,POST )
#define KFLOATINITD_NA(   DST, COUNT, PRE, POST )				KFLOATINITDSS_NA( DST,DST,DST,COUNT,PRE,POST )
#define KFLOATINITDS_NA(  DST, SRC0, COUNT, PRE, POST )			KFLOATINITDSS_NA( DST,SRC0,SRC0,COUNT,PRE,POST )
#define KFLOATINITDSS_NA( DST, SRC0, SRC1, COUNT, PRE, POST )\
	__asm	mov		eax,COUNT							\
	__asm	mov		PRE,0								\
	__asm	and		eax,8-1								\
	__asm	mov		ebx,COUNT							\
	__asm	mov		POST,eax							\
	__asm	and		ebx,0xfffffff8						\
	__asm	je		done								\
	__asm	shl		ebx,2								\
	__asm	mov		edx,SRC0							\
	__asm	mov		esi,SRC1							\
	__asm	mov		edi,DST								\
	__asm	add		edx,ebx								\
	__asm	add		esi,ebx								\
	__asm	add		edi,ebx								\
	__asm	mov		eax,edx								\
	__asm	or		eax,esi								\
	__asm	or		eax,edi								\
	__asm	neg		ebx									\

/*
	when OPER is called:
	edx = s0
	esi	= s1
	edi	= d
	ebx	= index*4

	xmm0 & xmm1	must not be trashed
*/
#define KMOVDS1( DST, SRC0 )							\
	__asm	movss	xmm2,SRC0							\
	__asm	movss	DST,xmm2
#define KMOVDS4( DST, SRC0 )							\
	__asm	movups	xmm2,SRC0							\
	__asm	movups	DST,xmm2
#define KMINDS1( DST, SRC0 )							\
	__asm	movss	xmm2,SRC0							\
	__asm	minss	DST,xmm2
#define KMAXDS1( DST, SRC0 )							\
	__asm	movss	xmm2,SRC0							\
	__asm	maxss	DST,xmm2

// general ALU operation
#define KALUDSS1( OP, DST, SRC0, SRC1 )					\
	__asm	movss	xmm2,SRC0							\
	__asm	OP##ss	xmm2,SRC1							\
	__asm	movss	DST,xmm2
#define KALUDSS4( OP, DST, SRC0, SRC1 )					\
	__asm	movups	xmm2,SRC0							\
	__asm	movups	xmm3,SRC1							\
	__asm	OP##ps	xmm2,xmm3							\
	__asm	movups	DST,xmm2

#define KADDDSS1( DST, SRC0, SRC1 )		KALUDSS1( add, DST,SRC0,SRC1 )
#define KADDDSS4( DST, SRC0, SRC1 )		KALUDSS4( add, DST,SRC0,SRC1 )
#define KSUBDSS1( DST, SRC0, SRC1 )		KALUDSS1( sub, DST,SRC0,SRC1 )
#define KSUBDSS4( DST, SRC0, SRC1 )		KALUDSS4( sub, DST,SRC0,SRC1 )
#define KMULDSS1( DST, SRC0, SRC1 )		KALUDSS1( mul, DST,SRC0,SRC1 )
#define KMULDSS4( DST, SRC0, SRC1 )		KALUDSS4( mul, DST,SRC0,SRC1 )

#define KDIVDSS1( DST, SRC0, SRC1 )						\
	__asm	movss	xmm2,SRC1							\
	__asm	rcpss	xmm3,xmm2							\
	__asm	mulss	xmm2,xmm3							\
	__asm	mulss	xmm2,xmm3							\
	__asm	addss	xmm3,xmm3							\
	__asm	subss	xmm3,xmm2							\
	__asm	mulss	xmm3,SRC0							\
	__asm	movss	DST,xmm3
#define KDIVDSS4( DST, SRC0, SRC1 )						\
	__asm	movups	xmm2,SRC1							\
	__asm	rcpps	xmm3,xmm2							\
	__asm	mulps	xmm2,xmm3							\
	__asm	mulps	xmm2,xmm3							\
	__asm	addps	xmm3,xmm3							\
	__asm	subps	xmm3,xmm2							\
	__asm	movups	xmm2,SRC0							\
	__asm	mulps	xmm3,xmm2							\
	__asm	movups	DST,xmm3
#define	KF2IDS1( SRC0 )									\
	__asm	movss		xmm2,SRC0						\
	__asm	cvttps2pi	mm2,xmm2						\
	__asm	movd		[edi+ebx],mm2
#define	KF2IDS4( SRC0 )									\
	__asm	movups		xmm2,SRC0						\
	__asm	cvttps2pi	mm2,xmm2						\
	__asm	movq		[edi+ebx+0],mm2					\
	__asm	shufps		xmm2,xmm2,SHUFFLEPS(1,0,3,2)	\
	__asm	cvttps2pi	mm2,xmm2						\
	__asm	movq		[edi+ebx+8],mm2
#define	KISQRTDS1( DST,SRC0 )							\
	__asm	movss	xmm2,SRC0							\
	__asm	rsqrtss	xmm3,xmm2							\
	__asm	mulss	xmm2,xmm3							\
	__asm	mulss	xmm2,xmm3							\
	__asm	subss	xmm2,xmm1							\
	__asm	mulss	xmm3,xmm0							\
	__asm	mulss	xmm3,xmm2							\
	__asm	movss	DST,xmm3
#define	KISQRTDS4( DST,SRC0 )							\
	__asm	movups	xmm2,SRC0							\
	__asm	rsqrtps	xmm3,xmm2							\
	__asm	mulps	xmm2,xmm3							\
	__asm	mulps	xmm2,xmm3							\
	__asm	subps	xmm2,xmm1							\
	__asm	mulps	xmm3,xmm0							\
	__asm	mulps	xmm3,xmm2							\
	__asm	movups	DST,xmm3

// this is used in vector4 implementation to shift constant V4
#define KANDREGDSV( DST, SRC0, VALUE )					\
	__asm	mov		DST,SRC0							\
	__asm	and		DST,VALUE

// this is used in vector4 code to operate with float arrays as sources
#define KEXPANDFLOAT( DST, SRC )						\
	__asm	movss	DST,SRC								\
	__asm	shufps  DST,DST,0

#define	KADDDS1( DST,SRC )		KADDDSS1( DST,DST,SRC )
#define	KADDDS4( DST,SRC )		KADDDSS4( DST,DST,SRC )
#define	KSUBDS1( DST,SRC )		KSUBDSS1( DST,DST,SRC )
#define	KSUBDS4( DST,SRC )		KSUBDSS4( DST,DST,SRC )
#define	KMULDS1( DST,SRC )		KMULDSS1( DST,DST,SRC )
#define	KMULDS4( DST,SRC )		KMULDSS4( DST,DST,SRC )
#define	KDIVDS1( DST,SRC )		KDIVDSS1( DST,DST,SRC )
#define	KDIVDS4( DST,SRC )		KDIVDSS4( DST,DST,SRC )

// handles pre & post leftovers
#define	KFLOATOPER( OPER, OPER4, COUNT )				\
	__asm		mov		ecx,pre							\
	__asm		mov		ebx,COUNT						\
	__asm		cmp		ebx,ecx							\
	__asm		cmovl	ecx,COUNT						\
	__asm		test	ecx,ecx							\
	__asm		je		preDone							\
	__asm		xor		ebx,ebx							\
	__asm	lpPre:										\
				OPER									\
	__asm		add		ebx,4							\
	__asm		dec		ecx								\
	__asm		jg		lpPre							\
	__asm	preDone:									\
	__asm		mov		ecx,post						\
	__asm		mov		ebx,COUNT						\
	__asm		sub		ebx,ecx							\
	__asm		shl		ebx,2							\
	__asm		cmp		ecx,4							\
	__asm		jl		post4Done						\
				OPER4									\
	__asm		sub		ecx,4							\
	__asm		add		ebx,4*4							\
	__asm	post4Done:									\
	__asm		test	ecx,ecx							\
	__asm		je		postDone						\
	__asm	lpPost:										\
				OPER									\
	__asm		add		ebx,4							\
	__asm		dec		ecx								\
	__asm		jg		lpPost							\
	__asm	postDone:

// operate on a constant and a float array
//http://stackoverflow.com/questions/11277652/what-is-the-meaning-of-align-an-the-start-of-a-section

#define KFLOAT_CA( ALUOP, DST, SRC, CONSTANT, COUNT )	\
	int	pre,post;										\
	__asm		movss	xmm0,CONSTANT					\
	__asm		shufps	xmm0,xmm0,0						\
			KFLOATINITDS( DST, SRC, COUNT, pre, post )	\
	__asm		and		eax,15							\
	__asm		jne		lpNA							\
	__asm		jmp		lpA								\
	__asm		align	16								\
	__asm	lpA:										\
	__asm		prefetchnta	[edx+ebx+64]				\
	__asm		movaps	xmm1,xmm0						\
	__asm		movaps	xmm2,xmm0						\
	__asm		ALUOP##ps	xmm1,[edx+ebx]				\
	__asm		ALUOP##ps	xmm2,[edx+ebx+16]			\
	__asm		movaps	[edi+ebx],xmm1					\
	__asm		movaps	[edi+ebx+16],xmm2				\
	__asm		add		ebx,16*2						\
	__asm		jl		lpA								\
	__asm		jmp		done							\
	__asm		align	16								\
	__asm	lpNA:										\
	__asm		prefetchnta	[edx+ebx+64]				\
	__asm		movaps	xmm1,xmm0						\
	__asm		movaps	xmm2,xmm0						\
	__asm		movups	xmm3,[edx+ebx]					\
	__asm		movups	xmm4,[edx+ebx+16]				\
	__asm		ALUOP##ps	xmm1,xmm3					\
	__asm		ALUOP##ps	xmm2,xmm4					\
	__asm		movaps	[edi+ebx],xmm1					\
	__asm		movaps	[edi+ebx+16],xmm2				\
	__asm		add		ebx,16*2						\
	__asm		jl		lpNA							\
	__asm	done:										\
	__asm		mov		edx,SRC							\
	__asm		mov		edi,DST							\
	__asm		KFLOATOPER( KALUDSS1( ALUOP, [edi+ebx],xmm0,[edx+ebx] ),	\
	__asm					KALUDSS4( ALUOP, [edi+ebx],xmm0,[edx+ebx] ), COUNT )

// operate on two float arrays
#define KFLOAT_AA( ALUOP, DST, SRC0, SRC1, COUNT )		\
	int	pre,post;										\
	KFLOATINITDSS( DST, SRC0, SRC1, COUNT, pre, post )	\
	__asm		and		eax,15							\
	__asm		jne		lpNA							\
	__asm		jmp		lpA								\
	__asm		align	16								\
	__asm	lpA:										\
	__asm		movaps	xmm1,[edx+ebx]					\
	__asm		movaps	xmm2,[edx+ebx+16]				\
	__asm		ALUOP##ps	xmm1,[esi+ebx]				\
	__asm		ALUOP##ps	xmm2,[esi+ebx+16]			\
	__asm		prefetchnta	[edx+ebx+64]				\
	__asm		prefetchnta	[esi+ebx+64]				\
	__asm		movaps	[edi+ebx],xmm1					\
	__asm		movaps	[edi+ebx+16],xmm2				\
	__asm		add		ebx,16*2						\
	__asm		jl		lpA								\
	__asm		jmp		done							\
	__asm		align	16								\
	__asm	lpNA:										\
	__asm		movups	xmm1,[edx+ebx]					\
	__asm		movups	xmm2,[edx+ebx+16]				\
	__asm		movups	xmm3,[esi+ebx]					\
	__asm		movups	xmm4,[esi+ebx+16]				\
	__asm		prefetchnta	[edx+ebx+64]				\
	__asm		prefetchnta	[esi+ebx+64]				\
	__asm		ALUOP##ps	xmm1,xmm3					\
	__asm		ALUOP##ps	xmm2,xmm4					\
	__asm		movaps	[edi+ebx],xmm1					\
	__asm		movaps	[edi+ebx+16],xmm2				\
	__asm		add		ebx,16*2						\
	__asm		jl		lpNA							\
	__asm	done:										\
	__asm		mov		edx,SRC0						\
	__asm		mov		esi,SRC1						\
	__asm		mov		edi,DST							\
	KFLOATOPER( KALUDSS1( ALUOP, [edi+ebx],[edx+ebx],[esi+ebx] ),		\
				KALUDSS4( ALUOP, [edi+ebx],[edx+ebx],[esi+ebx] ), COUNT )


#endif
#if  defined(__GNUC__) || defined(_MSC_VER)

#define JOINTQUAT_SIZE				(8*4)
#define JOINTQUAT_SIZE_STR			"(8*4)"
#define JOINTMAT_SIZE				(4*3*4)
#define JOINTMAT_SIZE_STR			"(4*3*4)"

#define JOINTWEIGHT_SIZE			(4*4)
#define JOINTWEIGHT_SIZE_STR		"(4*4)"

#define ALIGN2_INIT1( X, INIT )				ALIGN16( static X[2] ) = { INIT, INIT }
#define ALIGN2_INIT2( X, I0, I1 )			ALIGN16( static X[2] ) = { I0, I1 }
#define ALIGN4_INIT1( X, INIT )				ALIGN16( static X[4] ) = { INIT, INIT, INIT, INIT }
#define ALIGN4_INIT4( X, I0, I1, I2, I3 )	ALIGN16( static X[4] ) = { I0, I1, I2, I3 }
#define ALIGN8_INIT1( X, INIT )				ALIGN16( static X[8] ) = { INIT, INIT, INIT, INIT, INIT, INIT, INIT, INIT }

#define SHUFFLEPS( x, y, z, w )		(( (x) & 3 ) << 6 | ( (y) & 3 ) << 4 | ( (z) & 3 ) << 2 | ( (w) & 3 ))
#define R_SHUFFLEPS( x, y, z, w )	(( (w) & 3 ) << 6 | ( (z) & 3 ) << 4 | ( (y) & 3 ) << 2 | ( (x) & 3 ))
#define R_SHUFFLEPS2( x )	1
#define SHUFFLEPS( x, y, z, w )		(( (x) & 3 ) << 6 | ( (y) & 3 ) << 4 | ( (z) & 3 ) << 2 | ( (w) & 3 ))
#define R_SHUFFLEPS( x, y, z, w )	(( (w) & 3 ) << 6 | ( (z) & 3 ) << 4 | ( (y) & 3 ) << 2 | ( (x) & 3 ))
#define SHUFFLEPD( x, y )			(( (x) & 1 ) << 1 | ( (y) & 1 ))
#define R_SHUFFLEPD( x, y )			(( (y) & 1 ) << 1 | ( (x) & 1 ))


ALIGN8_INIT1( unsigned short SIMD_W_zero, 0 );
ALIGN8_INIT1( unsigned short SIMD_W_maxShort, 1<<15 );
//(1 << 31 ) == 0x80000000  ==
//(1 << 23 ) == 0x00800000  ==
//(1 << 15 ) == 0x00008000  ==
// for (MATH_SSE2) && (MATH_AVX)  We can use the pshufd instruction, which was introduced in SSE2 32-bit integer ops.
/// Swizzles/permutes a single SSE register into another SSE register. Requires SSE2.
#define shuffle2_ps(reg, shuffle) _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128((reg)), (shuffle)))
 // We only have SSE 1, so must use the slightly worse shufps instruction, which always destroys the input operand - or we have AVX where we can use this operation without destroying input
#define shuffle1_ps(reg, shuffle) _mm_shuffle_ps((reg), (reg), (shuffle))


#endif

/*
===============================================================================

	SSE implementation of CSIMDProcessor

===============================================================================
*/
/**
 * \class CSIMD_SSE
 *
 * \ingroup SMF_Math
 *
 * \brief Implementação SSE do processador SIMD
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *

 * \note http://pt.wikipedia.org/wiki/SSE
 * \note http://neilkemp.us/src/sse_tutorial/sse_tutorial.html
 * \note http://arstechnica.com/features/2000/03/simd/
 * \note http://www.popoloski.com/posts/sse_move_instructions/
 **/
class SMF_API CSIMD_SSE : public CSIMD_MMX {
public:
	CSIMD_SSE():CSIMD_MMX(CPUID_SSE){};
	CSIMD_SSE(cpuid_t id):CSIMD_MMX(id){};
	~CSIMD_SSE(){};


#if defined(MACOS_X) && defined(__i386__)
	virtual const char * VPCALL getName() const;
	virtual void VPCALL dot( float *dst,			const CPlane &constant,const CVertex *src,	const int count );
	virtual	void VPCALL minMax( CVec3D &min,		CVec3D &max,			const CVertex *src,	const int *indexes,		const int count );
	virtual void VPCALL dot( float *dst,			const CVec3D &constant,	const CPlane *src,		const int count );
#else
	virtual const char * VPCALL getName() const;
	//Trigonometry=====================

	virtual float VPCALL  invSqrt( float x );
	virtual void  VPCALL  InvSqrt4( float x[4] );
	virtual float VPCALL  sinZeroHalfPI( float x );
	virtual void  VPCALL  sin4ZeroHalfPI( float entrada[4], float saida[4] );
	virtual float VPCALL  sin( float a );
	virtual void  VPCALL  sin4( float entrada[4], float saida[4] );
	virtual float VPCALL  cosZeroHalfPI( float x );
	virtual void  VPCALL  cos4ZeroHalfPI( float entrada[4], float saida[4] );
	virtual float VPCALL  cos( float a );
	virtual void  VPCALL  cos4( float entrada[4], float saida[4] );
	virtual void  VPCALL  sincos( float a, float &sin, float &cos );
	virtual void  VPCALL  sincos4( float a[4], float sin[4], float cos[4] );
	virtual float VPCALL  aTanPositive( float y, float x );
	virtual void  VPCALL  aTan4Positive( float y[4], float x[4], float resultado[4] );
	virtual float VPCALL  atan( float y, float x );
	virtual void  VPCALL  aTan4( float y[4], float x[4], float resultado[4] );

	//===================
	virtual void VPCALL add( float *dst,			const float constant,	const float *src,		const int count );
	virtual void VPCALL add( float *dst,			const float *src0,		const float *src1,		const int count );
	virtual void VPCALL sub( float *dst,			const float constant,	const float *src,		const int count );
	virtual void VPCALL sub( float *dst,			const float *src0,		const float *src1,		const int count );
	virtual void VPCALL mul( float *dst,			const float constant,	const float *src,		const int count );
	virtual void VPCALL mul( float *dst,			const float *src0,		const float *src1,		const int count );
	virtual void VPCALL div( float *dst,			const float constant,	const float *src,		const int count );
	virtual void VPCALL div( float *dst,			const float *src0,		const float *src1,		const int count );
	virtual void VPCALL mulAdd( float *dst,			const float constant,	const float *src,		const int count );
	virtual void VPCALL mulAdd( float *dst,			const float *src0,		const float *src1,		const int count );
	virtual void VPCALL mulSub( float *dst,			const float constant,	const float *src,		const int count );
	virtual void VPCALL mulSub( float *dst,			const float *src0,		const float *src1,		const int count );

	virtual void VPCALL dot( float *dst,			const CVec3D &constant,	const CVec3D *src,		const int count );
	virtual void VPCALL dot( float *dst,			const CVec3D &constant,	const CPlane *src,		const int count );
	virtual void VPCALL dot( float *dst,			const CVec3D &constant,	const CVertex *src,	const int count );
	virtual void VPCALL dot( float *dst,			const CPlane &constant,const CVec3D *src,		const int count );
	virtual void VPCALL dot( float *dst,			const CPlane &constant,const CPlane *src,		const int count );
	virtual void VPCALL dot( float *dst,			const CPlane &constant,const CVertex *src,	const int count );
	virtual void VPCALL dot( float *dst,			const CVec3D *src0,		const CVec3D *src1,		const int count );
	virtual void VPCALL dot( float &dot,			const float *src1,		const float *src2,		const int count );

	virtual void VPCALL cmpGT( sf_u8 *dst,			const float *src0,		const float constant,	const int count );
	virtual void VPCALL cmpGT( sf_u8 *dst,			const sf_u8 bitNum,		const float *src0,		const float constant,	const int count );
	virtual void VPCALL cmpGE( sf_u8 *dst,			const float *src0,		const float constant,	const int count );
	virtual void VPCALL cmpGE( sf_u8 *dst,			const sf_u8 bitNum,		const float *src0,		const float constant,	const int count );
	virtual void VPCALL cmpLT( sf_u8 *dst,			const float *src0,		const float constant,	const int count );
	virtual void VPCALL cmpLT( sf_u8 *dst,			const sf_u8 bitNum,		const float *src0,		const float constant,	const int count );
	virtual void VPCALL cmpLE( sf_u8 *dst,			const float *src0,		const float constant,	const int count );
	virtual void VPCALL cmpLE( sf_u8 *dst,			const sf_u8 bitNum,		const float *src0,		const float constant,	const int count );

	virtual void VPCALL minMax( float &min,			float &max,				const float *src,		const int count );
	virtual	void VPCALL minMax( CVec2D &min,		CVec2D &max,			const CVec2D *src,		const int count );
	virtual void VPCALL minMax( CVec3D &min,		CVec3D &max,			const CVec3D *src,		const int count );
	virtual	void VPCALL minMax( CVec3D &min,		CVec3D &max,			const CVertex *src,	const int count );
	virtual	void VPCALL minMax( CVec3D &min,		CVec3D &max,			const CVertex *src,	const int *indexes,		const int count );

	virtual void VPCALL clamp( float *dst,			const float *src,		const float min,		const float max,		const int count );
	virtual void VPCALL clampMin( float *dst,		const float *src,		const float min,		const int count );
	virtual void VPCALL clampMax( float *dst,		const float *src,		const float max,		const int count );

	virtual void VPCALL zero16( float *dst,			const int count );
	virtual void VPCALL negate16( float *dst,		const int count );
	virtual void VPCALL copy16( float *dst,			const float *src,		const int count );
	virtual void VPCALL add16( float *dst,			const float *src1,		const float *src2,		const int count );
	virtual void VPCALL sub16( float *dst,			const float *src1,		const float *src2,		const int count );
	virtual void VPCALL mul16( float *dst,			const float *src1,		const float constant,	const int count );
	virtual void VPCALL addAssign16( float *dst,	const float *src,		const int count );
	virtual void VPCALL subAssign16( float *dst,	const float *src,		const int count );
	virtual void VPCALL mulAssign16( float *dst,	const float constant,	const int count );

	virtual bool VPCALL matX_inverse_4x4(float* src);
	virtual void VPCALL matX_MultiplyVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec );
	virtual void VPCALL matX_MultiplyAddVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec );
	virtual void VPCALL matX_MultiplySubVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec );
	virtual void VPCALL matX_TransposeMultiplyVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec );
	virtual void VPCALL matX_TransposeMultiplyAddVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec );
	virtual void VPCALL matX_TransposeMultiplySubVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec );
	virtual void VPCALL matX_MultiplyMatX( CMatXD &dst, const CMatXD &m1, const CMatXD &m2 );
	virtual void VPCALL matX_TransposeMultiplyMatX( CMatXD &dst, const CMatXD &m1, const CMatXD &m2 );
	virtual void VPCALL matX_LowerTriangularSolve( const CMatXD &L, float *x, const float *b, const int n, int skip = 0 );
	virtual void VPCALL matX_LowerTriangularSolveTranspose( const CMatXD &L, float *x, const float *b, const int n );
	virtual bool VPCALL matX_LDLTFactor( CMatXD &mat, CVecXD &invDiag, const int n );

	virtual void VPCALL blendJoints( CJointQuaternion *joints, const CJointQuaternion *blendJoints, const float lerp, const int *index, const int numJoints );
	virtual void VPCALL convertJointQuatsToJointMats( CMatJoint3x4 *jointMats, const CJointQuaternion *jointQuats, const int numJoints );
	virtual void VPCALL convertJointMatsToJointQuats( CJointQuaternion *jointQuats, const CMatJoint3x4 *jointMats, const int numJoints );
	virtual void VPCALL transformJoints( CMatJoint3x4 *jointMats, const int *parents, const int firstJoint, const int lastJoint );
	virtual void VPCALL untransformJoints( CMatJoint3x4 *jointMats, const int *parents, const int firstJoint, const int lastJoint );
	virtual void VPCALL transformVerts( CVertex *verts, const int numVerts, const CMatJoint3x4 *joints, const CVec4D *weights, const int *index, const int numWeights );
	virtual void VPCALL tracePointCull( sf_u8 *cullBits, sf_u8 &totalOr, const float radius, const CPlane *planes, const CVertex *verts, const int numVerts );
	virtual void VPCALL decalPointCull( sf_u8 *cullBits, const CPlane *planes, const CVertex *verts, const int numVerts );
	virtual void VPCALL overlayPointCull( sf_u8 *cullBits, CVec2D *texCoords, const CPlane *planes, const CVertex *verts, const int numVerts );
	virtual void VPCALL deriveTriPlanes( CPlane *planes, const CVertex *verts, const int numVerts, const int *indexes, const int numIndexes );
	virtual void VPCALL  deriveTangents( CPlane *planes, CVertex *verts, const int numVerts, const int *indexes, const int numIndexes );
	//virtual void VPCALL deriveUnsmoothedTangents( CVertex *verts, const dominantTri_s *dominantTris, const int numVerts );
	virtual void VPCALL  normalizeTangents( CVertex *verts, const int numVerts );
	virtual void VPCALL  createTextureSpaceLightVectors( CVec3D *lightVectors, const CVec3D &lightOrigin, const CVertex *verts, const int numVerts, const int *indexes, const int numIndexes );
	virtual void VPCALL  createSpecularTextureCoords( CVec4D *texCoords, const CVec3D &lightOrigin, const CVec3D &viewOrigin, const CVertex *verts, const int numVerts, const int *indexes, const int numIndexes );
	virtual int  VPCALL  createShadowCache( CVec4D *vertexCache, int *vertRemap, const CVec3D &lightOrigin, const CVertex *verts, const int numVerts );
	virtual int  VPCALL  createVertexProgramShadowCache( CVec4D *vertexCache, const CVertex *verts, const int numVerts );

	virtual void VPCALL  upSamplePCMTo44kHz( float *dest, const short *pcm, const int numSamples, const int kHz, const int numChannels );
	virtual void VPCALL  upSampleOGGTo44kHz( float *dest, const float * const *ogg, const int numSamples, const int kHz, const int numChannels );
	virtual void VPCALL  mixSoundTwoSpeakerMono( float *mixBuffer, const float *samples, const int numSamples, const float lastV[2], const float currentV[2] );
	virtual void VPCALL  mixSoundTwoSpeakerStereo( float *mixBuffer, const float *samples, const int numSamples, const float lastV[2], const float currentV[2] );
	virtual void VPCALL  mixSoundSixSpeakerMono( float *mixBuffer, const float *samples, const int numSamples, const float lastV[6], const float currentV[6] );
	virtual void VPCALL  mixSoundSixSpeakerStereo( float *mixBuffer, const float *samples, const int numSamples, const float lastV[6], const float currentV[6] );
	virtual void VPCALL  mixedSoundToSamples( short *samples, const float *mixBuffer, const int numSamples );

//=================================


virtual void VPCALL vector3D_Sum(CVec3D* pOut, const CVec3D* pIn);
virtual void VPCALL vector3D_SumOf(CVec3D* pOut, const CVec3D* pIn1, const CVec3D* pIn2);
virtual void VPCALL vector3D_Diff(CVec3D* pLeft, CVec3D* pRight);
virtual void VPCALL vector3D_DiffOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight);
virtual void VPCALL vector3D_Scale(CVec3D* pOut, float scalar);
virtual void VPCALL vector3D_ScaleOf(CVec3D* pOut, const CVec3D* pIn, float scalar);
virtual float VPCALL vector3D_Dot(const CVec3D* pSrc1, const CVec3D* pSrc2);
virtual float VPCALL vector3D_LengthSq(const CVec3D* pVec);
virtual float VPCALL vector3D_Length(const CVec3D* pVec);
virtual void VPCALL vector3D_Cross(CVec3D* pLeft, const CVec3D* pRight);
virtual void VPCALL vector3D_CrossOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight);
virtual void VPCALL vector3D_Normalize(CVec3D* pVec);
virtual void VPCALL vector3D_NormalizeOf(CVec3D* pOut, const CVec3D* pVec);
virtual float VPCALL vector3D_Distance(const CVec3D* pVec1, const CVec3D* pVec2);

virtual float VPCALL vector4D_Dot(const CVec4D* pSrc1, const CVec4D* pSrc2);



virtual void VPCALL vector4D_Sum(CVec4D* pOut, const CVec4D* pIn);
virtual void VPCALL vector4D_SumOf(CVec4D* pOut, const CVec4D* pIn1, const CVec4D* pIn2);
virtual void VPCALL vector4D_Diff(CVec4D* pLeft, CVec4D* pRight);
virtual void VPCALL vector4D_DiffOf(CVec4D* pOut, const CVec4D* pLeft, const CVec4D* pRight);
virtual void VPCALL vector4D_Scale(CVec4D* pOut, float scalar);
virtual void VPCALL vector4D_ScaleOf(CVec4D* pOut, const CVec4D* pIn, float scalar);
virtual float VPCALL vector4D_LengthSq(const CVec4D* pVec);
virtual float VPCALL vector4D_Length(const CVec4D* pVec);
virtual void VPCALL vector4D_Normalize(CVec4D* pVec);
virtual void VPCALL vector4D_NormalizeOf(CVec4D* pOut, const CVec4D* pVec);
virtual float VPCALL vector4D_Distance(const CVec4D* pVec1, const CVec4D* pVec2);

virtual void VPCALL vector3D_AlignedSum(CVec3D* pOut, const CVec3D* pIn);
virtual void VPCALL vector3D_AlignedSumOf(CVec3D* pOut, const CVec3D* pIn1, const CVec3D* pIn2);
virtual void VPCALL vector3D_AlignedDiff(CVec3D* pLeft, CVec3D* pRight);
virtual void VPCALL vector3D_AlignedDiffOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight);
virtual void VPCALL vector3D_AlignedScale(CVec3D* pOut, float scalar);
virtual void VPCALL vector3D_AlignedScaleOf(CVec3D* pOut, const CVec3D* pIn, float scalar);
virtual float VPCALL vector3D_AlignedDot(const CVec3D* pSrc1, const CVec3D* pSrc2);
virtual float VPCALL vector3D_AlignedLengthSq(const CVec3D* pVec);
virtual float VPCALL vector3D_AlignedLength(const CVec3D* pVec);
virtual void VPCALL vector3D_AlignedCross(CVec3D* pLeft, const CVec3D* pRight);
virtual void VPCALL vector3D_AlignedCrossOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight);
virtual void VPCALL vector3D_AlignedNormalize(CVec3D* pVec);
virtual void VPCALL vector3D_AlignedNormalizeOf(CVec3D* pOut, const CVec3D* pVec);
virtual float VPCALL vector3D_AlignedDistance(const CVec3D* pVec1, const CVec3D* pVec2);



//====== CVec4D========================================================

virtual void VPCALL vector4D_AlignedSum(CVec4D* pOut, const CVec4D* pIn);
virtual void VPCALL vector4D_AlignedSumOf(CVec4D* pOut, const CVec4D* pIn1, const CVec4D* pIn2);
virtual void VPCALL vector4D_AlignedDiff(CVec4D* pLeft, CVec4D* pRight);
virtual void VPCALL vector4D_AlignedDiffOf(CVec4D* pOut, const CVec4D* pLeft, const CVec4D* pRight);
virtual void VPCALL vector4D_AlignedScale(CVec4D* pOut, float scalar);
virtual void VPCALL vector4D_AlignedScaleOf(CVec4D* pOut, const CVec4D* pIn, float scalar);
virtual float VPCALL vector4D_AlignedLengthSq(const CVec4D* pVec);
virtual float VPCALL vector4D_AlignedDot(const CVec4D* pSrc4D1, const CVec4D* pSrc4D2);
virtual float VPCALL vector4D_AlignedLength(const CVec4D* pVec);
virtual void VPCALL vector4D_AlignedNormalize(CVec4D* pVec);
virtual void VPCALL vector4D_AlignedNormalizeOf(CVec4D* pOut, const CVec4D* pVec);
virtual float VPCALL vector4D_AlignedDistance(const CVec4D* pVec1, const CVec4D* pVec2);

//============CMat4D===================================================

/**
\brief soma pMat e pIn e guarda o resultado em pMat
\param pMat Matriz que será somada a pIn e guardará o resultado da soma
\param pIn Matriz que será somada a pMat
\note Sums pMat and pIn and stores the result in pMat.
**/
virtual void  VPCALL mat4D_Sum(CMat4D* pMat, const CMat4D* pIn);
/**
\brief soma pIn1 e pIn2 e armazena o resultado em pMat
\brief Sums pIn1 and pIn2 and stores the result in pMat.
\param pMat Matriz que guardará o resultado da soma
\param pIn1 Matriz que será somada a pIn2
\param pIn2 Matriz que será somada a pIn1
**/
virtual void  VPCALL mat4D_SumOf(CMat4D* pMat, const CMat4D* pIn1, const CMat4D* pIn2);
/**
\brief subtrai pIn de pMat e armazena o resultado em pMat. (pMat = pMat - pIn)
\brief Subtracts pIn from pMat and stores the result in pMat.
\param pMat matriz que será utilizada para subtrair pIn e armazenará o resultado da subtração
\param pIn matriz que será subtrída de pMat
**/
virtual void  VPCALL mat4D_Diff(CMat4D* pMat, const CMat4D* pIn);
/**
\brief subtrai pRight de pLeft e armazena o resultado em pMat. (pMat = pLeft - pRight)
\param pMat matriz que armazenará o resultado da subtração.
\param pLeft matriz que será utilizada na subtração (pMat = pLeft - pRight)
\param pRight matriz que será utilizada na subtração (pMat = pLeft - pRight)
\brief Subtracts pRight from pLeft and stores the result in pMat.
**/
virtual void  VPCALL mat4D_DiffOf(CMat4D* pMat, const CMat4D* pLeft, const CMat4D* pRight);
/**
\brief multiplica pMat por scalare armazena o resultado em pMat. (pMat = pMat * scalar)
\param pMat matriz que será multiplicada e armazenará o resultado
\param scalar númerom que será multiplicador de pMat.
\note Scales the components of pMat by scalar and stores the result in pMat.
**/
virtual void  VPCALL mat4D_Scale(CMat4D* pMat, float scalar);
/**
\brief multiplica pMat por scalare armazena o resultado em pMat. (pMat = pMat * scalar)
\param pOut matriz que  armazenará o resultado
\param pIn matriz que será multiplicada por scalar
\param scalar número que será multiplicador de pMat.
\note Scales the components of pIn by scalar and stores the result in pMat.
**/
virtual void  VPCALL mat4D_ScaleOf(CMat4D* pOut, const CMat4D* pIn, float scalar);
/**
\brief Multiplica pLeft por pRight e armazena o resultado em pLeft. (pLeft = pRight * pLeft)
\param pLeft matriz
\param pRight matriz
\note Multiplies pLeft by pRight and stores the result in pLeft.
**/
virtual void  VPCALL mat4D_Multiply(CMat4D* pLeft, const CMat4D* pRight);
/**
\brief Multiplica pLeft por pRighte armazena o resultado em pOut. (pOut = pRight * pLeft)
\param pLeft matriz
\param pRight matriz
\param pOut Matriz que armazena o resultado da multiplicação
\brief Multiplies pLeft by pRight and stores the result in pLeft.
**/
virtual void  VPCALL mat4D_MultiplyOf(CMat4D* pOut, const CMat4D* pLeft, const CMat4D* pRight);
/**
\brief calcula a matriz transposta de pIn e armazena o resultado em pIn. (pIn = pInt)
\param pIn matriz que se vai calcular a transposta e armazenará o resultado
\note Transposes the matrix pIn stores the result in pIn.
**/
virtual void  VPCALL mat4D_Transpose(CMat4D* pIn);
/**
\brief calcula a matriz transposta de pIn e armazena o resultado em pOut. (pOut = pInt)
\param pIn matriz que se vai calcular a transposta
\param pOut matriz que armazenará o resultado
\note Transposes the matrix pIn stores the result in pOut.
**/
virtual void  VPCALL mat4D_TransposeOf(CMat4D* pOut, const CMat4D* pIn);
/**
\brief Multiplica o vetor pela matriz (pOut = pOut * pIn)
\param pOut Vetor multiplicador
\param pIn Matriz multiplicadora
\note Transforms the 3D vector (as 4D with w = 1.0) pVec by the matrix pMat and stores the result in pVec.
The w component is discarded. If you need to save the w component (for example, a software renderer uses it for perspective correction)
then use 4D vectors with mat4D_VectorMultiply(CVec4D* pOut4D, const CMat4D* pIn).

**/
virtual void  VPCALL mat4D_VectorMultiply(CVec3D* pOut, const CMat4D* pIn);
/**
\brief Multiplica o vetor pela matriz (pOut = pOut * pIn)
\param pOut Vetor que armazenará o resultado
\param pIn Vetor multiplicador
\param pMat Matriz multiplicadora
\note Transforms the 3D vector (as 4D with w = 1.0) pIn by the matrix pMat and stores the result in pOut.
The w component is discarded. If you need to save the w component (for example, a software renderer uses it for perspective correction)
then use 4D vectors with mat4D_VectorMultiply(CVec4D* pOut4D, const CVec4D* pIn4D, const CMat4D* pMat).

**/
virtual void  VPCALL mat4D_VectorMultiplyOf(CVec3D* pOut, const CVec3D* pIn, const CMat4D* pMat);
/**
\brief Multiplica o vetor pela matriz (pOut4D = pOut4D * pIn)
\param pOut4D Vetor multiplicador
\param pIn Matriz multiplicadora
\note Transforms the 4D vector pVec4D by the matrix pMat and stores the result in pVec4D.
**/
virtual void  VPCALL mat4D_VectorMultiply(CVec4D* pOut4D, const CMat4D* pIn);
/**
\brief Multiplica o vetor pela matriz (pOut4D = pIn4D * pMat)
\param pOut4D Vetor que armazenará o resultado
\param pIn4D Vetor multiplicador
\param pMat Matriz multiplicadora
\note Transforms the 4D vector pIn4D by the matrix pMat and stores the result in pOut4D.
**/
virtual void  VPCALL mat4D_VectorMultiplyOf(CVec4D* pOut4D, const CVec4D* pIn4D, const CMat4D* pMat);
/**
\brief rotaciona pMat
\param yaw  parâmetro de rotação no sentido yaw
\param pitch parâmetro de rotação no sentido pitch
\param roll parâmetro de rotação no sentido row
**/
virtual void  VPCALL mat4D_ToRotate(CMat4D* pMat, float yaw, float pitch, float roll);
/**
\brief rotaciona pMat
\param pYawPitchRoll  Vetor que contém os parâmetros de rotação no sentido yaw, pitch, roll de rotação

**/
virtual void  VPCALL mat4D_ToRotateOf(CMat4D* pMat, const CVec3D* pYawPitchRoll);

/**
\brief soma pMat e pIn e guarda o resultado em pMat
\note pMat,pIn devem estar alinhados em 16bytes
\param pMat Matriz que será somada a pIn e guardará o resultado da soma
\param pIn Matria que será somada a pMat
\note Sums pMat and pIn and stores the result in pMat.
\note pMat, pIn must be 16 bytes aligned
**/
virtual void  VPCALL mat4D_AlignedSum(CMat4D* pMat, const CMat4D* pIn);
/**
\brief subtrai pIn de pMat e armazena o resultado em pMat. (pMat = pMat - pIn)
\brief Sums pIn1 and pIn2 and stores the result in pMat.
\note  pMat,pIn1,pIn2 devem estar alinhados em 16bytes.
\param pMat matriz que será utilizada para subtrair pIn e armazenará o resultado da subtração
\param pIn matriz que será subtrída de pMat
\note  pMat,pIn1,pIn2 must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedSumOf(CMat4D* pMat, const CMat4D* pIn1, const CMat4D* pIn2);
/**
\brief subtrai pIn de pMat e armazena o resultado em pMat. (pMat = pMat - pIn)
\brief Subtracts pIn from pMat and stores the result in pMat.
\param pMat matriz que será utilizada para subtrair pIn e armazenará o resultado da subtração
\param pIn matriz que será subtrída de pMat
\note  pMat,pIn devem estar alinhados em 16bytes.
\note pMat,pIn must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedDiff(CMat4D* pMat, const CMat4D* pIn);
/**
\brief subtrai pRight de pLeft e armazena o resultado em pMat. (pMat = pLeft - pRight)
\param pMat matriz que armazenará o resultado da subtração.
\param pLeft matriz que será utilizada na subtração (pMat = pLeft - pRight)
\param pRight matriz que será utilizada na subtração (pMat = pLeft - pRight)
\note Subtracts pRight from pLeft and stores the result in pMat.
\note pMat,pLeft,pRight devem estar alinhados em 16bytes.
\note pMat,pLeft,pRight must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedDiffOf(CMat4D* pMat, const CMat4D* pLeft, const CMat4D* pRight);
/**
\brief multiplica pMat por scalare armazena o resultado em pMat. (pMat = pMat * scalar)
\param pMat matriz que será multiplicada e armazenará o resultado
\param scalar númerom que será multiplicador de pMat.
\note Scales the components of pMat by scalar and stores the result in pMat.
\note pMat deve estar alinhada em 16bytes.
\note pMat must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedScale(CMat4D* pMat, float scalar);
/**
\brief multiplica pMat por scalare armazena o resultado em pMat. (pMat = pMat * scalar)
\param pOut matriz que  armazenará o resultado
\param pIn matriz que será multiplicada por scalar
\param scalar número que será multiplicador de pMat.
\note Scales the components of pIn by scalar and stores the result in pMat.
\note pOut, pIn devem estar alinhados em 16bytes.
\note pOut, pIn must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedScaleOf(CMat4D* pOut, const CMat4D* pIn, float scalar);
/**
\brief Multiplica pLeft por pRight e armazena o resultado em pLeft. (pLeft = pRight * pLeft)
\param pLeft matriz
\param pRight matriz
\note Multiplies pLeft by pRight and stores the result in pLeft.
\note pLeft,pRight devem estar alinhados em 16bytes.
\note pLeft,pRight must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedMultiply(CMat4D* pLeft, const CMat4D* pRight);
/**
\brief Multiplica pLeft por pRighte armazena o resultado em pOut. (pOut = pRight * pLeft)
\param pLeft matriz
\param pRight matriz
\param pOut Matriz que armazena o resultado da multiplicação
\note Multiplies pLeft by pRight and stores the result in pLeft.
\note pOut,pLeft,pRight devem estar alinhados em 16bytes.
\note pOut,pLeft,pRight must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedMultiplyOf(CMat4D* pOut, const CMat4D* pLeft, const CMat4D* pRight);
/**
\brief calcula a matriz transposta de pIn e armazena o resultado em pIn. (pIn = pInt)
\param pIn matriz que se vai calcular a transposta e armazenará o resultado
\brief Transposes the matrix pIn stores the result in pIn.
\note pIn deve estar alinhada em 16bytes.
\note pIn must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedTranspose(CMat4D* pIn);
/**
\brief calcula a matriz transposta de pIn e armazena o resultado em pOut. (pOut = pInt)
\param pIn matriz que se vai calcular a transposta
\param pOut matriz que armazenará o resultado
\note Transposes the matrix pIn stores the result in pOut.
\note pOut, pIn devem estar alinhados em 16bytes.
\note pOut, pIn must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedTransposeOf(CMat4D* pOut, const CMat4D* pIn);
/**
\brief Multiplica o vetor pela matriz (pOut = pOut4D * pIn)
\brief Transforms the 4D vector pVec4D by the matrix pMat and stores the result in pVec4D.
\param pOut Vetor multiplicador e que armazenará o resultado
\param pIn Matriz multiplicadora
\brief Transforms the 3D vector (as 4D with w = 1.0) pVec by the matrix pMat and stores the result in pVec. The w component is discarded. If you need to save the w component (for example, a software renderer uses it for perspective correction) then use 4D vectors with Matrix_Vector4Multiply().
\note pOut, pIn devem estar alinhados em 16bytes.
\note pOut, pIn must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedVectorMultiply(CVec3D* pOut, const CMat4D* pIn);
/**
\brief Multiplica o vetor pela matriz (pOut4D = pIn4D * pMat)
\param pOut Vetor que armazenará o resultado
\param pIn Vetor multiplicador
\param pMat Matriz multiplicadora
\brief Transforms the 3D vector (as 4D with w = 1.0) pIn by the matrix pMat and stores the result in pOut. The w component is discarded. If you need to save the w component (for example, a software renderer uses it for perspective correction) then use 4D vectors with Matrix_Vector4MultiplyOf().
\note pOut, pIn, pMat devem estar alinhados em 16bytes.
\note pOut, pIn, pMat must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedVectorMultiplyOf(CVec3D* pOut, const CVec3D* pIn, const CMat4D* pMat);
/**
\brief Multiplica o vetor pela matriz (pOut4D = pOut4D * pMat)
\brief Transforms the 4D vector pVec4D by the matrix pMat and stores the result in pVec4D.
\param pOut4D Vetor multiplicador e que armazenará o resultado
\param pMat Matriz multiplicadora
\note pOut4D, pMat devem estar alinhados em 16bytes.
\note pOut4D, pMat must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedVectorMultiply(CVec4D* pOut4D, const CMat4D* pMat);
/**
\brief Multiplica o vetor pela matriz (pOut4D = pIn4D * pMat)
\param pOut4D Vetor que armazenará o resultado
\param pIn4D Vetor multiplicador
\param pMat Matriz multiplicadora
\note Transforms the 4D vector pIn4D by the matrix pMat and stores the result in pOut4D.
\note pOut4D, pIn4D, pMat must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedVectorMultiplyOf(CVec4D* pOut4D, const CVec4D* pIn4D, const CMat4D* pMat);
/**
\brief rotaciona pMat
\param yaw  parâmetro de rotação no sentido yaw
\param pitch parâmetro de rotação no sentido pitch
\param roll parâmetro de rotação no sentido row
\note pMat deve estar alinhada em 16bytes.
\note pMat must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedToRotate(CMat4D* pMat, float yaw, float pitch, float roll);
/**
\brief rotaciona pMat
\param pYawPitchRoll  Vetor que contém os parâmetros de rotação no sentido yaw, pitch, roll de rotação
\note pMat deve estar alinhada em 16bytes.
\note pMat must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedToRotateOf(CMat4D* pMat, const CVec3D* pYawPitchRoll);

virtual void  VPCALL quat_to_mat4x4(sf_m128 q, sf_m128 t, sf_m128 *m);

/// Compute the product M*v, where M is a 4x4 matrix denoted by an array of 4 sf_m128's, and v is a 4x1 vector.
// If we have SSE 4.1, we can use the dpps (dot product) instruction, _mm_dp_ps intrinsic.
/// Compute the product M*v, where M is a 4x4 matrix denoted by an array of 4 sf_m128's, and v is a 4x1 vector.
// If we have SSE3, we can repeatedly use haddps to accumulate the result.
virtual sf_m128  VPCALL mat4x4_mul_sse(const sf_m128 *matrix, sf_m128 vector);
virtual  void   VPCALL mat4x4_mul_sse(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2);


/// Compute the product M*v, where M is a 3x4 matrix denoted by an array of 4 sf_m128's, and v is a 4x1 vector.
virtual sf_m128  VPCALL mat3x4_mul_sse(const sf_m128 *matrix, sf_m128 vector);

virtual CVec3D  VPCALL mat3x4_mul_vec(const sf_m128 *matrix, sf_m128 vector);

/**
\brief multiplica duas matrizes CMat4D  m1*m2
\param [out] out resultado da multiplicação
\param m1 matriz a ser mutiplicada
\param m2 matria a ser multiplicada
**/
virtual void  VPCALL mat4x4_mul_dpps(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2);
virtual void  VPCALL mat4x4_mul_dpps_2(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2);
virtual void  VPCALL mat4x4_mul_dpps_3(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2);


virtual void  VPCALL mat3x4_mul_sse(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2);
virtual float  VPCALL mat4x4_inverse(const CMat4D *mat, CMat4D *out);

/// Inverts a 3x4 affine transformation matrix (in row-major format) that only consists of rotation (+possibly mirroring) and translation.
virtual void  VPCALL mat3x4_inverse_orthonormal(sf_m128 *mat, sf_m128 *out);
virtual sf_m128  VPCALL  newtonRhapsonRecipStep(sf_m128 recip, sf_m128 estimate);

virtual sf_m128  VPCALL  newtonRhapsonRecip(sf_m128 recip);

/// Computes the determinant of a 4x4 matrix.
virtual float  VPCALL mat4x4_determinant(const CMat4D *row);

/// Computes the determinant of a 3x4 matrix stored in row-major format. (Treated as a square matrix with last row [0,0,0,1])
virtual float  VPCALL mat3x4_determinant(const sf_m128 *row);

virtual void  VPCALL mat3x4_transpose(const sf_m128 *src, sf_m128 *dst);

virtual sf_m128  VPCALL colmajor_mat4x4_mul_sse(const sf_m128 *matrix, sf_m128 vector);

//=========================Quaternion========================================
/**
\brief Normaiza o Quaternion e armazena o resultado em pQuat
\param pQuat Quaternion a ser normalizado
\brief Normalizes pQuat and stores it in pQuat.
**/
virtual void VPCALL  quaternion_Normalize(CQuaternion* pQuat);
/**
\brief Normaiza o Quaternion e armazena o resultado em pOut
\param pQuat Quaternion a ser normalizado
\param pOut Quaternion que recebera o resultado
\brief Normalizes pQuat and stores it in pOut.
**/
virtual void VPCALL  quaternion_NormalizeOf(CQuaternion* pOut, const CQuaternion* pQuat);
/**
\brief Multiplica os Quaternions e armazena o resultado em pLeft (pLeft = pLeft * pRight)
\param pLeft Quaternion que será multiplicado e receberá o resultado
\param pRight Quaternion que será multiplicado
\brief Multiplies pLeft by pRight and stores the result in pLeft. The result is not explicitly normalized,
       but will be normal if two normal quaternions are used as inputs.
**/
virtual void VPCALL  quaternion_Multiply(CQuaternion* pLeft, const CQuaternion* pRight);
/**
\if pt_br
\brief Multiplica os Quaternions e armazena o resultado em pLeft (pLeft = pLeft * pRight)
\elseif us_en
\brief Multiplies pLeft by pRight and stores the result in pOut. The result is not explicitly normalized,
       but will be normal if two normal quaternions are used as inputs.
\endif
\param pLeft Quaternion que será multiplicado
\param pRight Quaternion que será multiplicado
\param pOut Quaternion que receberá o resultado
**/
virtual void VPCALL  quaternion_MultiplyOf(CQuaternion* pOut, const CQuaternion* pLeft, const CQuaternion* pRight);
/**
\brief create object CMatJoint3x4   from a quaternion
\param q quaternion
\param t
\param [out] matriz CMatJoint3x4
**/
virtual void VPCALL quat_to_mat3x4(sf_m128 q, sf_m128 t, sf_m128 *m);
virtual sf_m128 VPCALL quat_transform_vec4(sf_m128 quat, sf_m128 vec);
virtual sf_m128 VPCALL quat_mul_quat(sf_m128 q1, sf_m128 q2);
virtual sf_m128 VPCALL quat_div_quat(sf_m128 q1, sf_m128 q2);

//================CPlane===========================================
/**
\brief constrói um plano através de três pontos. O Plano não é normalizado
\param pOut O Plano que será criado
\param pA primeiro ponto
\param pB segundo ponto
\param pC terceiro ponto
\brief This function constructs a plane from three points. The plane is not normalized.
**/
virtual void   VPCALL  plane_FromPoints(CPlane* pOut, const CVec3D* pA, const CVec3D* pB,const CVec3D* pC);
/**
\brief retorna a distância sem sinal(módulo), entre um ponto e um plano
\param pPlane Plano a seu utilizado no cálculo
\param pPoint ponto a ser utilizado no cálculo
\brief This function returns the unsigned distance between a point and a plane.
**/
virtual float  VPCALL  plane_DistToPoint(const CPlane* pPlane, const CVec3D* pPoint);
/**
\brief calcula o produto escalar entre o plano e um ponto com w=1.0f. Util para classificar
       um ponto em relação ao plano: um valor menor que zero significa que o ponto está atraz do plano
	   e um valor maoir que zero significa que o ponto está na frente do plano
\param pPlane Plano que será utilizado no cálculo
\param pVec Ponto que será utilizado no cálculo
\brief This function returns the dot product between a plane and a point with w = 1.0.
This is useful for classifying a point in relation to a plane: a value of less than zero
implies the point it behind the plane, a value of zero implies the point is on the plane,
and value of greater than zero implies the point is in front of the plane.

**/
virtual float  VPCALL  plane_Dot(const CPlane* pPlane, const CVec3D* pVec);
/**
\brief retorna o produto escalar entre o plano e um ponto 4D
\param pPlane plano que será utilizado no cálculo do produto
\param pVec4 ponto de 4 dimensões que será utilizado no calculo do produto
\brief This function returns the dot product between a plane and a 4D point.
This is useful for classifying a point in relation to a plane: a value of less than zero implies the point it behind the plane, a value of zero implies the point is on the plane, and value of greater than zero implies the point is in front of the plane.

**/
virtual float  VPCALL  plane_Dot4(const CPlane* pPlane, const CVec4D* pVec4);
/**
\brief retorna o produoto escalar entre o plano e outra normal. representa o cosseno do angulo entre os dois, se ambos estiverem normalizados
\param pPlane Plano que será utilizado para calcular o produto
\param pVec normal do segundo plano
\return O resultado do produto
\brief This functions returns the dot product between the plane's normal and another normal.
This represents the cosine of the angle between the two, if both are normalized.

**/
virtual float  VPCALL  plane_DotNormal(const CPlane* pPlane, const CVec3D* pVec);
/**
\brief calcula o produto escalar entre dois planos normalizados.
\param pPlane1 Plano que será multiplicado
\param pPlane2 Plano que será multiplicado
\return o resultado do produto
\brief This functions returns the dot product between the two planes' normals.
This represents the cosine of the dihedral angle between the two, if both are normalized.
**/
virtual float  VPCALL  plane_DotPlane(const CPlane* pPlane1, const CPlane* pPlane2);
/**
\brief Normaliza o plano pOut e armazena o resultado nele mesmo
\param pOut plano que será normalizado e armazenará o resultado
\brief This function normalizes the plane such that the magnitude of its normal is 1 and modifies 'd' component appropriately.

**/
virtual void   VPCALL  plane_Normalize(CPlane* pOut);
/**
\brief Normaliza o plano e armazena o resultado em pOut
\param pOut plano que armazenará o resultado
\param pIn Plano que será normalizado
\brief This function stores the normalized plane pIn into pOut, while pIn is unaffected.
**/
virtual void   VPCALL  plane_NormalizeOf(CPlane* pOut, CPlane* pIn);

virtual	bool   VPCALL  intersectLineAABB(const CAABBox &box, const CVec4D &rayPos, const CVec4D &rayDir, float tNear, float tFar);



#endif
};

} //end MATH
} // end SMF
#endif /* !__MATH_SIMD_SSE_H__ */
