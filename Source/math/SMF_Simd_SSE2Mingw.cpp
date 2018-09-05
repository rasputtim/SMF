/*
  SMF -  Salvathor Math Fabric  (https://sourceforge.net/promects/sMfabric/)
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
#include "math/SMF_SimdGeneric.h"
#include "math/SMF_SimdMMX.h"
#include "math/SMF_SimdSSE.h"
#include "math/SMF_SimdSSE2.h"
#include "math/SMF_Matriz.h"
namespace SMF {

namespace MATH{
//===============================================================
//
//	SSE2 implementation of CSIMDProcessor
//
//===============================================================

#if defined(WIN32)&& defined(MASM_INTEL) //defined(__i386__)


/*
============
CSIMD_SSE2::GetName
============
*/
const char * CSIMD_SSE2::GetName() const {
	return "MMX & SSE & SSE2";
}

		// the SSE2 code is ungodly slow

/*
============
CSIMD_SSE2::MatX_LowerTriangularSolve

  solves x in Lx = b for the n * n sub-matrix of L
  if skip > 0 the first skip elements of x are assumed to be valid already
  L has to be a lower triangular matrix with (implicit) ones on the diagonal
  x == b is allowed
============
*/
void VPCALL CSIMD_SSE2::MatX_LowerTriangularSolve( const CMatrizXD &L, float *x, const float *b, const int n, int skip ) {
	int nc;
	const float *lptr;

	if ( skip >= n ) {
		return;
	}

	lptr = L[skip];
	nc = L.GetNumColumns();

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

	asm(
		"push		ebx\n"
		"mov			eax, %4\n"				// eax = i
		"shl			eax, 2\n"					// eax = i*4
		"mov			edx, %3\n"					// edx = n
		"shl			edx, 2\n"					// edx = n*4
		"mov			esi, %1\n"					// esi = x
		"mov			edi, %0\n"				// edi = lptr
		"add			esi, eax\n"
		"add			edi, eax\n"
		"mov			ebx, %2\n"					// ebx = b
		// aligned
	"looprow%=:\n"
		"mov			ecx, eax\n"
		"neg			ecx\n"
		"cvtps2pd	xmm0, [esi+ecx]\n"
		"cvtps2pd	xmm2, [edi+ecx]\n"
		"mulpd		xmm0, xmm2\n"
		"cvtps2pd	xmm1, [esi+ecx+8]\n"
		"cvtps2pd	xmm3, [edi+ecx+8]\n"
		"mulpd		xmm1, xmm3\n"
		"add			ecx, 20*4\n"
		"jg			donedot16%=\n"
	"dot16%=:\n"
		"cvtps2pd	xmm2, [esi+ecx-(16*4)]\n"
		"cvtps2pd	xmm3, [edi+ecx-(16*4)]\n"
		"cvtps2pd	xmm4, [esi+ecx-(14*4)]\n"
		"mulpd		xmm2, xmm3\n"
		"cvtps2pd	xmm5, [edi+ecx-(14*4)]\n"
		"addpd		xmm0, xmm2\n"
		"cvtps2pd	xmm2, [esi+ecx-(12*4)]\n"
		"mulpd		xmm4, xmm5\n"
		"cvtps2pd	xmm3, [edi+ecx-(12*4)]\n"
		"addpd		xmm1, xmm4\n"
		"cvtps2pd	xmm4, [esi+ecx-(10*4)]\n"
		"mulpd		xmm2, xmm3\n"
		"cvtps2pd	xmm5, [edi+ecx-(10*4)]\n"
		"addpd		xmm0, xmm2\n"
		"cvtps2pd	xmm2, [esi+ecx-(8*4)]\n"
		"mulpd		xmm4, xmm5\n"
		"cvtps2pd	xmm3, [edi+ecx-(8*4)]\n"
		"addpd		xmm1, xmm4\n"
		"cvtps2pd	xmm4, [esi+ecx-(6*4)]\n"
		"mulpd		xmm2, xmm3\n"
		"cvtps2pd	xmm5, [edi+ecx-(6*4)]\n"
		"addpd		xmm0, xmm2\n"
		"cvtps2pd	xmm2, [esi+ecx-(4*4)]\n"
		"mulpd		xmm4, xmm5\n"
		"cvtps2pd	xmm3, [edi+ecx-(4*4)]\n"
		"addpd		xmm1, xmm4\n"
		"cvtps2pd	xmm4, [esi+ecx-(2*4)]\n"
		"mulpd		xmm2, xmm3\n"
		"cvtps2pd	xmm5, [edi+ecx-(2*4)]\n"
		"addpd		xmm0, xmm2\n"
		"add		ecx, 16*4\n"
		"mulpd		xmm4, xmm5\n"
		"addpd		xmm1, xmm4\n"
		"jle			dot16%=\n"
	"donedot16%=:\n"
		"sub			ecx, 8*4\n"
		"jg			donedot8%=\n"
	"dot8%=:\n"
		"cvtps2pd	xmm2, [esi+ecx-(8*4)]\n"
		"cvtps2pd	xmm3, [edi+ecx-(8*4)]\n"
		"cvtps2pd	xmm7, [esi+ecx-(6*4)]\n"
		"mulpd		xmm2, xmm3\n"
		"cvtps2pd	xmm5, [edi+ecx-(6*4)]\n"
		"addpd		xmm0, xmm2\n"
		"cvtps2pd	xmm6, [esi+ecx-(4*4)]\n"
		"mulpd		xmm7, xmm5\n"
		"cvtps2pd	xmm3, [edi+ecx-(4*4)]\n"
		"addpd		xmm1, xmm7\n"
		"cvtps2pd	xmm4, [esi+ecx-(2*4)]\n"
		"mulpd		xmm6, xmm3\n"
		"cvtps2pd	xmm7, [edi+ecx-(2*4)]\n"
		"addpd		xmm0, xmm6\n"
		"add			ecx, 8*4\n"
		"mulpd		xmm4, xmm7\n"
		"addpd		xmm1, xmm4\n"
	"donedot8%=:\n"
		"sub			ecx, 4*4\n"
		"jg			donedot4%=\n"
	"dot4%=:\n"
		"cvtps2pd	xmm2, [esi+ecx-(4*4)]\n"
		"cvtps2pd	xmm3, [edi+ecx-(4*4)]\n"
		"cvtps2pd	xmm4, [esi+ecx-(2*4)]\n"
		"mulpd		xmm2, xmm3\n"
		"cvtps2pd	xmm5, [edi+ecx-(2*4)]\n"
		"addpd		xmm0, xmm2\n"
		"add			ecx, 4*4\n"
		"mulpd		xmm4, xmm5\n"
		"addpd		xmm1, xmm4\n"
	"donedot4%=:\n"
		"addpd		xmm0, xmm1\n"
		"movups		xmm1, xmm0\n"
		"shufpd		xmm1, xmm1, 0x01\n" /*R_SHUFFLEPD( 1, 0 )=01*/
		"addsd		xmm0, xmm1\n"
		"sub		ecx, 4*4\n"
		"jz			dot0%=\n"
		"add			ecx, 4\n"
		"jz			dot1%=\n"
		"add			ecx, 4\n"
		"jz			dot2%=\n"
	//dot3%=:\n"
		"cvtss2sd	xmm1, [esi-(3*4)]\n"
		"cvtss2sd	xmm2, [edi-(3*4)]\n"
		"mulsd		xmm1, xmm2\n"
		"addsd		xmm0, xmm1\n"
	"dot2%=:\n"
		"cvtss2sd	xmm3, [esi-(2*4)]\n"
		"cvtss2sd	xmm4, [edi-(2*4)]\n"
		"mulsd		xmm3, xmm4\n"
		"addsd		xmm0, xmm3\n"
	"dot1%=:\n"
		"cvtss2sd	xmm5, [esi-(1*4)]\n"
		"cvtss2sd	xmm6, [edi-(1*4)]\n"
		"mulsd		xmm5, xmm6\n"
		"addsd		xmm0, xmm5\n"
	"dot0%=:\n"
		"cvtss2sd	xmm1, [ebx+eax]\n"
		"subsd		xmm1, xmm0\n"
		"cvtsd2ss	xmm0, xmm1\n"
		"movss		[esi], xmm0\n"
		"add			eax, 4\n"
		"cmp			eax, edx\n"
		"jge			done%=\n"
		"add			esi, 4\n"
		"mov			ecx, %5\n"
		"shl			ecx, 2\n"
		"add			edi, ecx\n"
		"add			edi, 4\n"
		"jmp			looprow%=\n"
		// done
	"done%=:\n"
		"pop			ebx\n"
	::"m"(lptr), "m"(x), "m"(b), "m"(n), "m"(skip),"m"(nc)
    :"eax","esi","edi","ebx","ecx");
}

/*
============
CSIMD_SSE2::MatX_LowerTriangularSolveTranspose

  solves x in L'x = b for the n * n sub-matrix of L
  L has to be a lower triangular matrix with (implicit) ones on the diagonal
  x == b is allowed
============
*/
void VPCALL CSIMD_SSE2::MatX_LowerTriangularSolveTranspose( const CMatrizXD &L, float *x, const float *b, const int n ) {
	int nc;
	const float *lptr;

	lptr = L.ToFloatPtr();
	nc = L.GetNumColumns();

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

	int i, j, m;
	float *xptr;
	double s0;

	// if the number of columns is not a multiple of 2 we're screwed for alignment.
	// however, if the number of columns is a multiple of 2 but the number of to be
	// processed rows is not a multiple of 2 we can still run 8 SGF_Byte aligned
	m = n;
	if ( m & 1 ) {
		m--;
		x[m] = b[m];

		lptr = L[m] + m - 4;
		xptr = x + m;
		asm(
			"push		ebx\n"
			"mov			eax, %3\n"					// eax = i
			"mov			esi, %1\n"				// esi = xptr
			"mov			edi, %0\n"				// edi = lptr
			"mov			ebx, %2\n"					// ebx = b
			"mov			edx, %4	\n"				// edx = nc*sizeof(float)
			"shl			edx, 2\n"
		"process4rows_1%=:\n"
			"cvtps2pd	xmm0, [ebx+eax*4-16]\n"	// load b[i-2], b[i-1]
			"cvtps2pd	xmm2, [ebx+eax*4-8]\n"		// load b[i-4], b[i-3]
			"xor			ecx, ecx\n"
			"sub			eax, %3\n"
			"neg			eax\n"
			"jz			done4x4_1%=\n"
		"process4x4_1%=:\n"	// process 4x4 blocks
			"cvtps2pd	xmm3, [edi]\n"
			"cvtps2pd	xmm4, [edi+8]\n"
			"add		edi, edx\n"
			"cvtss2sd	xmm5, [esi+4*ecx+0]\n"
			"shufpd		xmm5, xmm5, 0x00\n"  /*R_SHUFFLEPD( 0, 0 )*/
			"mulpd		xmm3, xmm5\n"
			"cvtps2pd	xmm1, [edi]\n"
			"mulpd		xmm4, xmm5\n"
			"cvtps2pd	xmm6, [edi+8]\n"
			"subpd		xmm0, xmm3\n"
			"subpd		xmm2, xmm4\n"
			"add		edi, edx\n"
			"cvtss2sd	xmm7, [esi+4*ecx+4]\n"
			"shufpd		xmm7, xmm7, 0x00\n"  /*R_SHUFFLEPD( 0, 0 )*/
			"mulpd		xmm1, xmm7\n"
			"cvtps2pd	xmm3, [edi]\n"
			"mulpd		xmm6, xmm7\n"
			"cvtps2pd	xmm4, [edi+8]\n"
			"subpd		xmm0, xmm1\n"
			"subpd		xmm2, xmm6\n"
			"add		edi, edx\n"
			"cvtss2sd	xmm5, [esi+4*ecx+8]\n"
			"shufpd		xmm5, xmm5, 0x00\n"  /*R_SHUFFLEPD( 0, 0 )*/
			"mulpd		xmm3, xmm5\n"
			"cvtps2pd	xmm1, [edi]\n"
			"mulpd		xmm4, xmm5\n"
			"cvtps2pd	xmm6, [edi+8]\n"
			"subpd		xmm0, xmm3\n"
			"subpd		xmm2, xmm4\n"
			"add			edi, edx\n"
			"cvtss2sd	xmm7, [esi+4*ecx+12]\n"
			"shufpd		xmm7, xmm7, 0x00\n"  /*R_SHUFFLEPD( 0, 0 )*/
			"mulpd		xmm1, xmm7\n"
			"add			ecx, 4\n"
			"mulpd		xmm6, xmm7\n"
			"cmp			ecx, eax\n"
			"subpd		xmm0, xmm1\n"
			"subpd		xmm2, xmm6\n"
			"jl			process4x4_1%=\n"
		"done4x4_1%=:\n"		// process left over of the 4 rows
			"cvtps2pd	xmm3, [edi]\n"
			"cvtps2pd	xmm4, [edi+8]\n"
			"cvtss2sd	xmm5, [esi+4*ecx]\n"
			"shufpd		xmm5, xmm5, 0x00\n"  /*R_SHUFFLEPD( 0, 0 )*/
			"mulpd		xmm3, xmm5\n"
			"mulpd		xmm4, xmm5\n"
			"subpd		xmm0, xmm3\n"
			"subpd		xmm2, xmm4\n"
			"imul		ecx, edx\n"
			"sub			edi, ecx\n"
			"neg			eax\n"

			"add			eax, %3\n"
			"sub			eax, 4\n"
			"movapd		xmm1, xmm0\n"
			"shufpd		xmm1, xmm1, 0X03\n"  /*R_SHUFFLEPD( 1, 1 )*/
			"movapd		xmm3, xmm2\n"
			"shufpd		xmm3, xmm3, 0X03\n"  /*R_SHUFFLEPD( 1, 1 )*/
			"sub			edi, edx\n"
			"cvtsd2ss	xmm7, xmm3\n"
			"movss		[esi-4], xmm7\n"			// xptr[-1] = s3
			"movsd		xmm4, xmm3\n"
			"movsd		xmm5, xmm3\n"
			"cvtss2sd	xmm7, [edi+8]\n"
			"mulsd		xmm3, xmm7\n"				// lptr[-1*nc+2] * s3
			"cvtss2sd	xmm7, [edi+4]\n"
			"mulsd		xmm4, xmm7\n"				// lptr[-1*nc+1] * s3
			"cvtss2sd	xmm7, [edi]\n"
			"mulsd		xmm5, xmm7\n"				// lptr[-1*nc+0] * s3
			"subsd		xmm2, xmm3\n"
			"cvtsd2ss	xmm7, xmm2\n"
			"movss		[esi-8], xmm7\n"			// xptr[-2] = s2
			"movsd		xmm6, xmm2\n"
			"sub			edi, edx\n"
			"subsd		xmm0, xmm5\n"
			"subsd		xmm1, xmm4\n"
			"cvtss2sd	xmm7, [edi+4]\n"
			"mulsd		xmm2, xmm7\n"				// lptr[-2*nc+1] * s2
			"cvtss2sd	xmm7, [edi]\n"
			"mulsd		xmm6, xmm7\n"				// lptr[-2*nc+0] * s2
			"subsd		xmm1, xmm2\n"
			"cvtsd2ss	xmm7, xmm1\n"
			"movss		[esi-12], xmm7\n"			// xptr[-3] = s1
			"subsd		xmm0, xmm6\n"
			"sub			edi, edx\n"
			"cmp			eax, 4\n"
			"cvtss2sd	xmm7, [edi]\n"
			"mulsd		xmm1, xmm7\n"				// lptr[-3*nc+0] * s1
			"subsd		xmm0, xmm1\n"
			"cvtsd2ss	xmm7, xmm0\n"
			"movss		[esi-16], xmm7\n"			// xptr[-4] = s0
			"jl			done4rows_1%=\n"
			"sub			edi, edx\n"
			"sub			edi, 16\n"
			"sub			esi, 16\n"
			"jmp			process4rows_1%=\n"
		"done4rows_1%=:\n"
			"pop			ebx\n"
		::"m"(lptr), "m"(xptr), "m"(b), "m"(m),"m"(nc)
        :);
	}
	else {
		lptr = L.ToFloatPtr() + m * L.GetNumColumns() + m - 4;
		xptr = x + m;
		asm(
			"push		ebx\n"
			"mov			eax, %3\n"					// eax = i
			"mov			esi, %1\n"				// esi = xptr
			"mov			edi, %0\n"				// edi = lptr
			"mov			ebx, %2\n"					// ebx = b
			"mov			edx, %4\n"					// edx = nc*sizeof(float)
			"shl			edx, 2\n"
		"process4rows%=:\n"
			"cvtps2pd	xmm0, [ebx+eax*4-16]\n"	// load b[i-2], b[i-1]
			"cvtps2pd	xmm2, [ebx+eax*4-8]\n"		// load b[i-4], b[i-3]
			"sub			eax, %3\n"
			"jz			done4x4%=\n"
			"neg			eax\n"
			"xor			ecx, ecx\n"
		"process4x4%=:\n"		// process 4x4 blocks
			"cvtps2pd	xmm3, [edi]\n"
			"cvtps2pd	xmm4, [edi+8]\n"
			"add			edi, edx\n"
			"cvtss2sd	xmm5, [esi+4*ecx+0]\n"
			"shufpd		xmm5, xmm5, 0x00\n"  /*R_SHUFFLEPD( 0, 0 )*/
			"mulpd		xmm3, xmm5\n"
			"cvtps2pd	xmm1, [edi]\n"
			"mulpd		xmm4, xmm5\n"
			"cvtps2pd	xmm6, [edi+8]\n"
			"subpd		xmm0, xmm3\n"
			"subpd		xmm2, xmm4\n"
			"add			edi, edx\n"
			"cvtss2sd	xmm7, [esi+4*ecx+4]\n"
			"shufpd		xmm7, xmm7, 0x00\n"  /*R_SHUFFLEPD( 0, 0 )*/
			"mulpd		xmm1, xmm7\n"
			"cvtps2pd	xmm3, [edi]\n"
			"mulpd		xmm6, xmm7\n"
			"cvtps2pd	xmm4, [edi+8]\n"
			"subpd		xmm0, xmm1\n"
			"subpd		xmm2, xmm6\n"
			"add			edi, edx\n"
			"cvtss2sd	xmm5, [esi+4*ecx+8]\n"
			"shufpd		xmm5, xmm5, 0x00\n"  /*R_SHUFFLEPD( 0, 0 )*/
			"mulpd		xmm3, xmm5\n"
			"cvtps2pd	xmm1, [edi]\n"
			"mulpd		xmm4, xmm5\n"
			"cvtps2pd	xmm6, [edi+8]\n"
			"subpd		xmm0, xmm3\n"
			"subpd		xmm2, xmm4\n"
			"add			edi, edx\n"
			"cvtss2sd	xmm7, [esi+4*ecx+12]\n"
			"shufpd		xmm7, xmm7, 0x00\n"  /*R_SHUFFLEPD( 0, 0 )*/
			"mulpd		xmm1, xmm7\n"
			"add			ecx, 4\n"
			"mulpd		xmm6, xmm7\n"
			"cmp			ecx, eax\n"
			"subpd		xmm0, xmm1\n"
			"subpd		xmm2, xmm6\n"
			"jl			process4x4%=\n"
			"imul		ecx, edx\n"
			"sub			edi, ecx\n"
			"neg			eax\n"
		"done4x4%=:\n"		// process left over of the 4 rows
			"add			eax, %3\n"
			"sub			eax, 4\n"
			"movapd		xmm1, xmm0\n"
			"shufpd		xmm1, xmm1, 0X03\n"  /*R_SHUFFLEPD( 1, 1 )*/
			"movapd		xmm3, xmm2\n"
			"shufpd		xmm3, xmm3, 0X03\n"  /*R_SHUFFLEPD( 1, 1 )*/
			"sub			edi, edx\n"
			"cvtsd2ss	xmm7, xmm3\n"
			"movss		[esi-4], xmm7\n"			// xptr[-1] = s3
			"movsd		xmm4, xmm3\n"
			"movsd		xmm5, xmm3\n"
			"cvtss2sd	xmm7, [edi+8]\n"
			"mulsd		xmm3, xmm7\n"				// lptr[-1*nc+2] * s3
			"cvtss2sd	xmm7, [edi+4]\n"
			"mulsd		xmm4, xmm7\n"				// lptr[-1*nc+1] * s3
			"cvtss2sd	xmm7, [edi]\n"
			"mulsd		xmm5, xmm7\n"				// lptr[-1*nc+0] * s3
			"subsd		xmm2, xmm3\n"
			"cvtsd2ss	xmm7, xmm2\n"
			"movss		[esi-8], xmm7\n"			// xptr[-2] = s2
			"movsd		xmm6, xmm2\n"
			"sub			edi, edx\n"
			"subsd		xmm0, xmm5\n"
			"subsd		xmm1, xmm4\n"
			"cvtss2sd	xmm7, [edi+4]\n"
			"mulsd		xmm2, xmm7\n"				// lptr[-2*nc+1] * s2
			"cvtss2sd	xmm7, [edi]\n"
			"mulsd		xmm6, xmm7\n"				// lptr[-2*nc+0] * s2
			"subsd		xmm1, xmm2\n"
			"cvtsd2ss	xmm7, xmm1\n"
			"movss		[esi-12], xmm7\n"			// xptr[-3] = s1
			"subsd		xmm0, xmm6\n"
			"sub			edi, edx\n"
			"cmp			eax, 4\n"
			"cvtss2sd	xmm7, [edi]\n"
			"mulsd		xmm1, xmm7	\n"			// lptr[-3*nc+0] * s1
			"subsd		xmm0, xmm1\n"
			"cvtsd2ss	xmm7, xmm0\n"
			"movss		[esi-16], xmm7\n"			// xptr[-4] = s0
			"jl			done4rows%=\n"
			"sub			edi, edx\n"
			"sub			edi, 16\n"
			"sub			esi, 16\n"
			"jmp			process4rows%=\n"
		"done4rows%=:\n"
			"pop			ebx\n"
		::"m"(lptr), "m"(xptr), "m"(b), "m"(m),"m"(nc)
        :);
}

	// process left over rows
	for ( i = (m&3)-1; i >= 0; i-- ) {
		s0 = b[i];
		lptr = L[i+1] + i;
		for ( j = i + 1; j < m; j++ ) {
			s0 -= lptr[0] * x[j];
			lptr += nc;
		}
		x[i] = s0;
	}
}


/*
============
CSIMD_SSE2::MixedSoundToSamples
============
*/
void VPCALL CSIMD_SSE2::MixedSoundToSamples( short *samples, const float *mixBuffer, const int numSamples ) {


	assert( ( numSamples % MIXBUFFER_SAMPLES ) == 0 );

	asm(

		"mov			eax, %2\n" //numSamples
		"mov			edi, %0\n" //mixBuffer
		"mov			esi, %1\n" //samples
		"shl			eax, 2\n"
		"add			edi, eax\n"
		"neg			eax\n"

	"loop16%=:\n"

		"movups		xmm0, [edi+eax+0*16]\n"
		"movups		xmm1, [edi+eax+1*16]\n"
		"movups		xmm2, [edi+eax+2*16]\n"
		"movups		xmm3, [edi+eax+3*16]\n"

		"add			esi, 4*4*2\n"

		"cvtps2dq	xmm4, xmm0\n"
		"cvtps2dq	xmm5, xmm1\n"
		"cvtps2dq	xmm6, xmm2\n"
		"cvtps2dq	xmm7, xmm3\n"

		"prefetchnta	[edi+eax+128]\n"

		"packssdw	xmm4, xmm5\n"
		"packssdw	xmm6, xmm7\n"

		"add			eax, 4*16\n"

		"movlps		[esi-4*4*2], xmm4\n"		// FIXME: should not use movlps/movhps to move integer data
		"movhps		[esi-3*4*2], xmm4\n"
		"movlps		[esi-2*4*2], xmm6\n"
		"movhps		[esi-1*4*2], xmm6\n"

		"jl			loop16%=\n"
		:
    : "m"(mixBuffer), "m"(samples) ,"m"(numSamples)
    :"eax","esi","edi");


}

#endif /* _WIN32 */

} //end MATH
}//end SMF
