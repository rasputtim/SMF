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

//#include "../precompiled.h"
//#pragma hdrstop

#include "math/SMF_SimdGeneric.h"
#include "math/SMF_SimdMMX.h"
#include "math/SMF_Simd3DNow.h"


//===============================================================
//
//	3DNow! implementation of CSIMDProcessor
//
//===============================================================
namespace SMF {
namespace MATH{
#if !defined (__GNUC__)
#if defined(_MSC_VER)

/*
============
CSIMD_3DNow::getName
============
*/
const char * CSIMD_3DNow::getName() const {
	return "MMX & 3DNow!";
}

// Very optimized memcpy() routine for all AMD Athlon and Duron family.
// This code uses any of FOUR different basic copy methods, depending
// on the transfer size.
// NOTE:  Since this code uses MOVNTQ (also known as "Non-Temporal MOV" or
// "Streaming Store"), and also uses the software prefetchnta instructions,
// be sure you're running on Athlon/Duron or other recent CPU before calling!

#define TINY_BLOCK_COPY 64       // upper limit for movsd type copy
// The smallest copy uses the X86 "movsd" instruction, in an optimized
// form which is an "unrolled loop".

#define IN_CACHE_COPY 64 * 1024  // upper limit for movq/movq copy w/SW prefetch
// Next is a copy that uses the MMX registers to copy 8 bytes at a time,
// also using the "unrolled loop" optimization.   This code uses
// the software prefetch instruction to get the data into the cache.

#define UNCACHED_COPY 197 * 1024 // upper limit for movq/movntq w/SW prefetch
// For larger blocks, which will spill beyond the cache, it's faster to
// use the Streaming Store instruction MOVNTQ.   This write instruction
// bypasses the cache and writes straight to main memory.  This code also
// uses the software prefetch instruction to pre-read the data.
// USE 64 * 1024 FOR THIS VALUE IF YOU'RE ALWAYS FILLING A "CLEAN CACHE"

#define BLOCK_PREFETCH_COPY  infinity // no limit for movq/movntq w/block prefetch
#define CACHEBLOCK 80h // number of 64-sf_u8 blocks (cache lines) for block prefetch
// For the largest size blocks, a special technique called Block Prefetch
// can be used to accelerate the read operations.   Block Prefetch reads
// one address per cache line, for a series of cache lines, in a short loop.
// This is faster than using software prefetch.  The technique is great for
// getting maximum read bandwidth, especially in DDR memory systems.

/*
================
CSIMD_3DNow::memCopy

  optimized memory copy routine that handles all alignment cases and block sizes efficiently
================
*/
void VPCALL CSIMD_3DNow::memCopy( void *dest, const void *src, const int n ) {



  __asm {

	mov		ecx, [n]					// number of bytes to copy
	mov		edi, [dest]					// destination
	mov		esi, [src]					// source
	mov		ebx, ecx					// keep a copy of count

	cld
	cmp		ecx, TINY_BLOCK_COPY
	jb		$memcpy_ic_3				// tiny? skip mmx copy

	cmp		ecx, 32*1024				// don't align between 32k-64k because
	jbe		$memcpy_do_align			//  it appears to be slower
	cmp		ecx, 64*1024
	jbe		$memcpy_align_done
$memcpy_do_align:
	mov		ecx, 8						// a trick that's faster than rep movsb...
	sub		ecx, edi					// align destination to qword
	and		ecx, 111b					// get the low bits
	sub		ebx, ecx					// update copy count
	neg		ecx							// set up to jump into the array
	add		ecx, offset $memcpy_align_done
	jmp		ecx							// jump to array of movsb's

align 4
	movsb
	movsb
	movsb
	movsb
	movsb
	movsb
	movsb
	movsb

$memcpy_align_done:						// destination is dword aligned
	mov		ecx, ebx					// number of bytes left to copy
	shr		ecx, 6						// get 64-sf_u8 block count
	jz		$memcpy_ic_2				// finish the last few bytes

	cmp		ecx, IN_CACHE_COPY/64		// too big 4 cache? use uncached copy
	jae		$memcpy_uc_test

// This is small block copy that uses the MMX registers to copy 8 bytes
// at a time.  It uses the "unrolled loop" optimization, and also uses
// the software prefetch instruction to get the data into the cache.
align 16
$memcpy_ic_1:							// 64-sf_u8 block copies, in-cache copy

	prefetchnta [esi + (200*64/34+192)]	// start reading ahead

	movq	mm0, [esi+0]				// read 64 bits
	movq	mm1, [esi+8]
	movq	[edi+0], mm0				// write 64 bits
	movq	[edi+8], mm1				//    note:  the normal movq writes the
	movq	mm2, [esi+16]				//    data to cache; a cache line will be
	movq	mm3, [esi+24]				//    allocated as needed, to store the data
	movq	[edi+16], mm2
	movq	[edi+24], mm3
	movq	mm0, [esi+32]
	movq	mm1, [esi+40]
	movq	[edi+32], mm0
	movq	[edi+40], mm1
	movq	mm2, [esi+48]
	movq	mm3, [esi+56]
	movq	[edi+48], mm2
	movq	[edi+56], mm3

	add		esi, 64						// update source pointer
	add		edi, 64						// update destination pointer
	dec		ecx							// count down
	jnz		$memcpy_ic_1				// last 64-sf_u8 block?

$memcpy_ic_2:
	mov		ecx, ebx					// has valid low 6 bits of the sf_u8 count
$memcpy_ic_3:
	shr		ecx, 2						// dword count
	and		ecx, 1111b					// only look at the "remainder" bits
	neg		ecx							// set up to jump into the array
	add		ecx, offset $memcpy_last_few
	jmp		ecx							// jump to array of movsd's

$memcpy_uc_test:
	cmp		ecx, UNCACHED_COPY/64		// big enough? use block prefetch copy
	jae		$memcpy_bp_1

$memcpy_64_test:
	or		ecx, ecx					// tail end of block prefetch will jump here
	jz		$memcpy_ic_2				// no more 64-sf_u8 blocks left

// For larger blocks, which will spill beyond the cache, it's faster to
// use the Streaming Store instruction MOVNTQ.   This write instruction
// bypasses the cache and writes straight to main memory.  This code also
// uses the software prefetch instruction to pre-read the data.
align 16
$memcpy_uc_1:							// 64-sf_u8 blocks, uncached copy

	prefetchnta [esi + (200*64/34+192)]	// start reading ahead

	movq	mm0,[esi+0]					// read 64 bits
	add		edi,64						// update destination pointer
	movq	mm1,[esi+8]
	add		esi,64						// update source pointer
	movq	mm2,[esi-48]
	movntq	[edi-64], mm0				// write 64 bits, bypassing the cache
	movq	mm0,[esi-40]				//    note: movntq also prevents the CPU
	movntq	[edi-56], mm1				//    from READING the destination address
	movq	mm1,[esi-32]				//    into the cache, only to be over-written
	movntq	[edi-48], mm2				//    so that also helps performance
	movq	mm2,[esi-24]
	movntq	[edi-40], mm0
	movq	mm0,[esi-16]
	movntq	[edi-32], mm1
	movq	mm1,[esi-8]
	movntq	[edi-24], mm2
	movntq	[edi-16], mm0
	dec		ecx
	movntq	[edi-8], mm1
	jnz		$memcpy_uc_1				// last 64-sf_u8 block?

	jmp		$memcpy_ic_2				// almost done

// For the largest size blocks, a special technique called Block Prefetch
// can be used to accelerate the read operations.   Block Prefetch reads
// one address per cache line, for a series of cache lines, in a short loop.
// This is faster than using software prefetch, in this case.
// The technique is great for getting maximum read bandwidth,
// especially in DDR memory systems.
$memcpy_bp_1:							// large blocks, block prefetch copy

	cmp		ecx, CACHEBLOCK				// big enough to run another prefetch loop?
	jl		$memcpy_64_test				// no, back to regular uncached copy

	mov		eax, CACHEBLOCK / 2			// block prefetch loop, unrolled 2X
	add		esi, CACHEBLOCK * 64		// move to the top of the block
align 16
$memcpy_bp_2:
	mov		edx, [esi-64]				// grab one address per cache line
	mov		edx, [esi-128]				// grab one address per cache line
	sub		esi, 128					// go reverse order
	dec		eax							// count down the cache lines
	jnz		$memcpy_bp_2				// keep grabbing more lines into cache

	mov		eax, CACHEBLOCK				// now that it's in cache, do the copy
align 16
$memcpy_bp_3:
	movq	mm0, [esi   ]				// read 64 bits
	movq	mm1, [esi+ 8]
	movq	mm2, [esi+16]
	movq	mm3, [esi+24]
	movq	mm4, [esi+32]
	movq	mm5, [esi+40]
	movq	mm6, [esi+48]
	movq	mm7, [esi+56]
	add		esi, 64						// update source pointer
	movntq	[edi   ], mm0				// write 64 bits, bypassing cache
	movntq	[edi+ 8], mm1				//    note: movntq also prevents the CPU
	movntq	[edi+16], mm2				//    from READING the destination address
	movntq	[edi+24], mm3				//    into the cache, only to be over-written,
	movntq	[edi+32], mm4				//    so that also helps performance
	movntq	[edi+40], mm5
	movntq	[edi+48], mm6
	movntq	[edi+56], mm7
	add		edi, 64						// update dest pointer

	dec		eax							// count down

	jnz		$memcpy_bp_3				// keep copying
	sub		ecx, CACHEBLOCK				// update the 64-sf_u8 block count
	jmp		$memcpy_bp_1				// keep processing chunks

// The smallest copy uses the X86 "movsd" instruction, in an optimized
// form which is an "unrolled loop".   Then it handles the last few bytes.
align 4
	movsd
	movsd								// perform last 1-15 dword copies
	movsd
	movsd
	movsd
	movsd
	movsd
	movsd
	movsd
	movsd								// perform last 1-7 dword copies
	movsd
	movsd
	movsd
	movsd
	movsd
	movsd

$memcpy_last_few:						// dword aligned from before movsd's
	mov		ecx, ebx					// has valid low 2 bits of the sf_u8 count
	and		ecx, 11b					// the last few cows must come home
	jz		$memcpy_final				// no more, let's leave
	rep		movsb						// the last 1, 2, or 3 bytes

$memcpy_final:
	emms								// clean up the MMX state
	sfence								// flush the write buffer
	mov		eax, [dest]					// ret value = destination pointer

    }

}

#elif defined (__GNUC__)/* _GNUC_VER */

//=================GNU C CODE==========================================



/*
	TODO:
	Optimize distance() -- It currently uses already implemented functions
	Optimize 3DNow! cross product -- Too many memory references
*/

#if 0 //TODO: Suport for 3DNow
void VPCALL CSIMD_3DNow::vector3D_Sum(CVec3D* pOut, const CVec3D* pIn)
{


	/* 3DNow! Implementation */
	asm(
	"movq  (%0), %%mm0\n"
	"movq 8(%0), %%mm1\n"
	"pfadd  (%1), %%mm0\n"
	"pfadd 8(%1), %%mm1\n"
	"movq %%mm0,  (%0)\n"
	"movq %%mm1, 8(%0)\n"
	:
	: "r" (pOut), "r" (pIn)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif




}

void VPCALL CSIMD_3DNow::vector3D_SumOf(CVec3D* pOut, const CVec3D* pIn1, const CVec3D* pIn2)
{

	/* 3DNow! Implementation */
	asm(
	"movq  (%0), %%mm0\n"
	"movq 8(%0), %%mm1\n"
	"pfadd  (%1), %%mm0\n"
	"pfadd 8(%1), %%mm1\n"
	"movq %%mm0,  (%2)\n"
	"movq %%mm1, 8(%2)\n"
	:
	: "r" (pIn1), "r" (pIn2), "r" (pOut)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif



}


void VPCALL CSIMD_3DNow::vector3D_Diff(CVec3D* pLeft, CVec3D* pRight)
{

	/* 3DNow! Implementation */
	asm(
	"movq  (%0), %%mm0\n"
	"movq 8(%0), %%mm1\n"
	"pfsub  (%1), %%mm0\n"
	"pfsub 8(%1), %%mm1\n"
	"movq %%mm0,  (%0)\n"
	"movq %%mm1, 8(%0)\n"
	:
	: "r" (pLeft), "r" (pRight)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif

}

void VPCALL CSIMD_3DNow::vector3D_DiffOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight)
{
	/* 3DNow! Implementation */
	asm(
	"movq  (%0), %%mm0\n"
	"movq 8(%0), %%mm1\n"
	"pfsub  (%1), %%mm0\n"
	"pfsub 8(%1), %%mm1\n"
	"movq %%mm0,  (%2)\n"
	"movq %%mm1, 8(%2)\n"
	:
	: "r" (pLeft), "r" (pRight), "r" (pOut)

	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif

}

void VPCALL CSIMD_3DNow::vector3D_Scale(CVec3D* pOut, float scalar)
{
	asm(
	"movd %1, %%mm0\n"
	"punpckldq %%mm0, %%mm0\n"
	"movq %%mm0, %%mm1\n"
	"pfmul  (%0), %%mm0\n"
	"pfmul 8(%0), %%mm1\n"
	"movq %%mm0,  (%0)\n"
	"movq %%mm1, 8(%0)\n"
	:
	: "r" (pOut), "m" (scalar)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif


}

void VPCALL CSIMD_3DNow::vector3D_ScaleOf(CVec3D* pOut, const CVec3D* pIn, float scalar)
{

	asm(
	"movd %2, %%mm0\n"
	"punpckldq %%mm0, %%mm0\n"
	"movq %%mm0, %%mm1\n"
	"pfmul  (%1), %%mm0\n"
	"pfmul 8(%1), %%mm1\n"
	"movq %%mm0,  (%0)\n"
	"movq %%mm1, 8(%0)\n"
	:
	: "r" (pOut), "r" (pIn), "m" (scalar)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif

}

float  VPCALL CSIMD_3DNow::vector3D_Dot(const CVec3D* pSrc1, const CVec3D* pSrc2)
{
	float dummy;
	asm(

	"movq  (%2), %%mm0\n"   /* mm0 = Y1 | X1*/
	"movq 8(%2), %%mm1\n"   /* mm1 = ?? | Z1 */

	"pfmul  (%1), %%mm0\n"  /* mm0 = Y1*Y2 | X1*X2 */
	"pand _SIMDx86_float_3DNOW_NO_W_MASK, %%mm1\n" /* mm1 = ?? & 0 | Z1 & 0xFFFFFFFF */
	"pfmul 8(%1), %%mm1\n"  /* mm1 =   0.0 | Z1*Z2 */

	/* Sum the pieces */
	"pfadd %%mm0, %%mm1\n"  /* mm1 = Y1*Y2 | X1*X2 + Z1*Z2 */
	"pfacc %%mm1, %%mm1\n"  /* mm0 = ? | X1*X2 + Y1*Y2 + Z1*Z2 */

	/* Store result */
	"movd %%mm1, -4(%%esp)\n"

	/* Yikes,'femms' must be done, sorry! Even if NO_EMMS is defined, the x87 register holds the result... */
	"femms\n"
	"flds -4(%%esp)\n"
	: "=t" (dummy)
	: "r" (pSrc1), "r" (pSrc2)
	);
	return dummy;

}

float  VPCALL CSIMD_3DNow::vector3D_Dot4(const float* pSrc4D1, const float* pSrc4D2)
{
	/* 3DNow!/3DNow+ Implementation */

	float dummy;
	asm(
	"movq  (%1), %%mm0\n"	/* mm0 = y1 | x1 */
	"movq 8(%1), %%mm1\n"	/* mm1 = w1 | z1 */
	"pfmul  (%2), %%mm0\n"	/* mm0 = y1*y2 | x1*x2 */
	"pfmul 8(%2), %%mm1\n"	/* mm1 = w1*w2 | z1*z2 */
	"pfadd %%mm0, %%mm1\n"	/* mm1 = y1*y2+w1*w2 | x1*x2+z1*z2 */
	"pfacc %%mm1, %%mm1\n"			/* mm1 = dot4 | dot4 */
	"movd %%mm1, -4(%%esp)\n"
	"femms\n"
	"flds -4(%%esp)\n"
	: "=t" (dummy)
	: "r" (pSrc4D1), "r" (pSrc4D2)
	);
	return dummy;


}

float  VPCALL CSIMD_3DNow::vector3D_LengthSq(const CVec3D* pVec)
{
	float dummy;
	asm(
	"movq  (%1), %%mm0\n"
	"movd 8(%1), %%mm1\n"

	/* square each component */
	"pfmul %%mm0, %%mm0\n"
	"pfmul %%mm1, %%mm1\n"

	/* Sum components */
	"pfacc %%mm0, %%mm0\n"
	"pfadd %%mm1, %%mm0\n"

	/* Move from MMX register to stack, then 'femms' (sorry!) and load x87 register */
	"movd %%mm0, -4(%%esp)\n"
	"femms\n"
	"flds -4(%%esp)\n"
	: "=t" (dummy)
	: "r" (pVec)
	);
	return dummy;

}

float  VPCALL CSIMD_3DNow::vector3D_Length(const CVec3D* pVec)
{

	float dummy;
	asm(
	"movq  (%1), %%mm0\n"
	"movd 8(%1), %%mm1\n"

	/* square each component */
	"pfmul %%mm0, %%mm0\n"
	"pfmul %%mm1, %%mm1\n"

	/* Sum components */
	"pfacc %%mm0, %%mm0\n"
	"pfadd %%mm1, %%mm0\n"

	#if defined(HIPREC)
	"movq %%mm0, %%mm3\n"
	"pfrsqrt %%mm0, %%mm1\n"        /* mm1 = x0 = rsqrt(mag*mag) */
	"movq %%mm1, %%mm2\n"           /* mm2 = x0 */
	"pfmul %%mm1, %%mm1\n"          /* mm1 = x1 = x0*x0 */
	"pfrsqit1 %%mm1, %%mm0\n"       /* mm0 = x2 = pfrsqit1(val, x1)*/
	"pfrcpit2 %%mm2, %%mm0\n"       /* mm2 = 1/sqrt(val) = pfrcpit2(x2, x0)*/
	"pfmul %%mm0, %%mm3\n"

	"movd %%mm3, -4(%%esp)\n"
	"femms\n"
	"flds -4(%%esp)\n"
	#else
	"pfrsqrt %%mm0, %%mm1\n"
	"pfmul %%mm0, %%mm1\n"
	"movd %%mm1, -4(%%esp)\n"
	"femms\n"
	"flds -4(%%esp)\n"

	#endif

	: "=t" (dummy)
	: "r" (pVec)
	);

	return dummy;


}

void VPCALL CSIMD_3DNow::vector3D_Cross(CVec3D* pLeft, const CVec3D* pRight)
{
	/* UGH! There has to be a much better way of accessing the memory then these fragmented ones */
	asm(
	"movq 4(%1), %%mm0\n"		/* mm0 = Rz | Ry */
	"movd  (%1), %%mm1\n"		/* mm1 =  ? | Rz */
	"movd 4(%1), %%mm5\n"		/* mm5 =  ? | Ry */
	"movd 8(%1), %%mm4\n"		/* mm4 =  ? | Rz */
	"punpckldq (%1), %%mm4\n"	/* mm4 = Rx | Rz */
	"movd 4(%0), %%mm3\n"		/* mm3 =  ? | Ly */
	"movd 8(%0), %%mm2\n"		/* mm2 =  ? | Lz */
	"punpckldq (%0), %%mm2\n"	/* mm2 = Lx | Lz */
	"movq 4(%0), %%mm6\n"		/* mm6 = Lz | Ly */
	"movd  (%0), %%mm7\n"		/* mm7 =  ? | Lx */


/*
	pLeft->x = (tmp.y * pRight->z) - (tmp.z * pRight->y);
	pLeft->y = (tmp.z * pRight->x) - (tmp.x * pRight->z);
	pLeft->z = (tmp.x * pRight->y) - (tmp.y * pRight->x);
*/

	/* SO UGLY, 10 memory accesses compared to SSE's 2 */

	/* multiply Column 1 & 2 */
	"pfmul %%mm4, %%mm6\n"
	"pfmul %%mm5, %%mm7\n"

	/* multiply Column 3 & 4 */
	"pfmul %%mm0, %%mm2\n"
	"pfmul %%mm1, %%mm3\n"

	/* Sutract... (1*2) - (3*4) */
	"pfsub %%mm2, %%mm6\n"
	"pfsub %%mm3, %%mm7\n"

	/* And store */
	"movq %%mm6,  (%0)\n"
	"movd %%mm7, 8(%0)\n"
	:
	: "r" (pLeft), "r" (pRight)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif


}

void VPCALL CSIMD_3DNow::vector3D_CrossOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight)
{
	/* UGH! There has to be a much better way of accessing the memory then these fragmented ones */
	asm(
	"movq 4(%2), %%mm0\n"		/* mm0 = Rz | Ry */
	"movd  (%2), %%mm1\n"		/* mm1 =  ? | Rz */
	"movd 4(%2), %%mm5\n"		/* mm5 =  ? | Ry */
	"movd 8(%2), %%mm4\n"		/* mm4 =  ? | Rz */
	"punpckldq (%2), %%mm4\n"	/* mm4 = Rx | Rz */
	"movd 4(%1), %%mm3\n"		/* mm3 =  ? | Ly */
	"movd 8(%1), %%mm2\n"		/* mm2 =  ? | Lz */
	"punpckldq (%1), %%mm2\n"	/* mm2 = Lx | Lz */
	"movq 4(%1), %%mm6\n"		/* mm6 = Lz | Ly */
	"movd  (%1), %%mm7\n"		/* mm7 =  ? | Lx */


/*
		pOut->x = (pLeft->y * pRight->z) - (pLeft->z * pRight->y);
		pOut->y = (pLeft->z * pRight->x) - (pLeft->x * pRight->z);
		pOut->z = (pLeft->x * pRight->y) - (pLeft->y * pRight->x);
*/
	/* SO UGLY, 10 memory accesses compared to SSE's 2 */

	/* multiply Column 1 & 2 */
	"pfmul %%mm4, %%mm6\n"	/* mm4 = Lz*Rx | Ly*Rz */
	"pfmul %%mm5, %%mm7\n"	/* mm5 = ????? | Lx*Ry */

	/* multiply Column 3 & 4 */
	"pfmul %%mm0, %%mm2\n"
	"pfmul %%mm1, %%mm3\n"

	/* Sutract... (1*2) - (3*4) */
	"pfsub %%mm2, %%mm6\n"
	"pfsub %%mm3, %%mm7\n"

	/* And store */
	"movq %%mm6,  (%0)\n"
	"movd %%mm7, 8(%0)\n"
	:
	: "r" (pOut), "r" (pLeft), "r" (pRight)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif

}
void VPCALL CSIMD_3DNow::vector3D_Normalize(CVec3D* pVec)
{
	asm (

	/* mm0-mm1 = pVec (duplicated in mm6-mm7) */
	"movq  (%0), %%mm0\n"
	"movq 8(%0), %%mm1\n"
	"pand _SIMDx86_float_3DNOW_NO_W_MASK, %%mm1\n"    /* mm1 = 0.0 | z */

	"movq %%mm0, %%mm6\n"
	"movq %%mm1, %%mm7\n"

	/* square components */
	"pfmul %%mm0, %%mm0\n"	/* y*y | x*x */
	"pfmul %%mm1, %%mm1\n"	/* 0   | z*z */

	/* Sum */
	"pfadd %%mm0, %%mm1\n"	/* 0+y*y | z*z+x*x */
	"pfacc %%mm1, %%mm1\n"	/* mag*mag | mag*mag */

	#ifdef HIPREC

	/* Get reciprocal, use Newton-Raphson HW to improve accuracy, then multiply by the high precision reciprocal. */
	"movq %%mm1, %%mm3\n"
	"pfrsqrt %%mm1, %%mm0\n"        /* mm1 = x0 = rsqrt(mag*mag) */
	"movq %%mm0, %%mm2\n"           /* mm2 = x0 */
	"pfmul %%mm0, %%mm0\n"          /* mm1 = x1 = x0*x0 */
	"pfrsqit1 %%mm0, %%mm1\n"       /* mm0 = x2 = pfrsqit1(val, x1)*/
	"pfrcpit2 %%mm2, %%mm1\n"       /* mm2 = 1/sqrt(val) = pfrcpit2(x2, x0)*/
	"pfmul %%mm1, %%mm6\n"
	"pfmul %%mm1, %%mm7\n"

	#else

	/* Get fast reciprocal and return */
	"pfrsqrt %%mm1, %%mm1\n"	/* mm1 = rcp(mag) | rcp(mag) (approx) */
	"pfmul %%mm1, %%mm6\n"		/* mm2 = y/mag | x/mag */
	"pfmul %%mm1, %%mm7\n"		/* mm3 = 0 | z/mag */
	#endif

	/* Store results */
	"movq %%mm6,  (%0)\n"
	"movq %%mm7, 8(%0)\n"
	:
	: "r" (pVec)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif

}
void VPCALL CSIMD_3DNow::vector3D_NormalizeOf(CVec3D* pOut, const CVec3D* pVec)
{
	asm (

	/* mm0-mm1 = pVec (duplicated in mm2-mm3) */
	"movq  (%0), %%mm0\n"
	"movq 8(%0), %%mm1\n"
	"pand _SIMDx86_float_3DNOW_NO_W_MASK, %%mm1\n"    /* mm1 = 0.0 | z */
	"movq %%mm0, %%mm6\n"
	"movq %%mm1, %%mm7\n"

	/* square components */
	"pfmul %%mm0, %%mm0\n"	/* y*y | x*x */
	"pfmul %%mm1, %%mm1\n"	/* 0   | z*z */

	/* Sum */
	"pfadd %%mm0, %%mm1\n"	/* 0+y*y | z*z+x*x */
	"pfacc %%mm1, %%mm1\n"	/* mag*mag | mag*mag */

	#ifdef HIPREC

	/* Get reciprocal, use Newton-Raphson HW to improve accuracy, then multiply by the high precision reciprocal. */
	"movq %%mm1, %%mm3\n"
	"pfrsqrt %%mm1, %%mm0\n"        /* mm1 = x0 = rsqrt(mag*mag) */
	"movq %%mm0, %%mm2\n"           /* mm2 = x0 */
	"pfmul %%mm0, %%mm0\n"          /* mm1 = x1 = x0*x0 */
	"pfrsqit1 %%mm0, %%mm1\n"       /* mm0 = x2 = pfrsqit1(val, x1)*/
	"pfrcpit2 %%mm2, %%mm1\n"       /* mm2 = 1/sqrt(val) = pfrcpit2(x2, x0)*/
	"pfmul %%mm1, %%mm6\n"
	"pfmul %%mm1, %%mm7\n"

	#else

	/* Get fast reciprocal and return */
	"pfrsqrt %%mm1, %%mm1\n"	/* mm1 = rcp(mag) | rcp(mag) (approx) */
	"pfmul %%mm1, %%mm6\n"		/* mm2 = y/mag | x/mag */
	"pfmul %%mm1, %%mm7\n"		/* mm3 = 0 | z/mag */
	#endif

	/* Store results */
	"movq %%mm6,  (%1)\n"
	"movq %%mm7, 8(%1)\n"
	:
	: "r" (pVec), "r" (pOut)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif


}

float  VPCALL CSIMD_3DNow::vector3D_Distance(const CVec3D* pVec1, const CVec3D* pVec2)
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

void VPCALL CSIMD_3DNow::vector3D_AlignedSum(CVec3D* pOut, const CVec3D* pIn)
{

	/* 3DNow! Implementation */
	asm(
	"movq  (%0), %%mm0\n"
	"movq 8(%0), %%mm1\n"
	"pfadd  (%1), %%mm0\n"
	"pfadd 8(%1), %%mm1\n"
	"movq %%mm0,  (%0)\n"
	"movq %%mm1, 8(%0)\n"
	:
	: "r" (pOut), "r" (pIn)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif

}

void VPCALL CSIMD_3DNow::vector3D_AlignedSumOf(CVec3D* pOut, const CVec3D* pIn1, const CVec3D* pIn2)
{
	/* 3DNow! Implementation */
	asm(
	"movq  (%0), %%mm0\n"
	"movq 8(%0), %%mm1\n"
	"pfadd  (%1), %%mm0\n"
	"pfadd 8(%1), %%mm1\n"
	"movq %%mm0,  (%2)\n"
	"movq %%mm1, 8(%2)\n"
	:
	: "r" (pIn1), "r" (pIn2), "r" (pOut)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif

}


void VPCALL CSIMD_3DNow::vector3D_AlignedDiff(CVec3D* pLeft, CVec3D* pRight)
{

	/* 3DNow! Implementation */
	asm(
	"movq  (%0), %%mm0\n"
	"movq 8(%0), %%mm1\n"
	"pfsub  (%1), %%mm0\n"
	"pfsub 8(%1), %%mm1\n"
	"movq %%mm0,  (%0)\n"
	"movq %%mm1, 8(%0)\n"
	:
	: "r" (pLeft), "r" (pRight)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif

}

void VPCALL CSIMD_3DNow::vector3D_AlignedDiffOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight)
{
	/* 3DNow! Implementation */
	asm(
	"movq  (%0), %%mm0\n"
	"movq 8(%0), %%mm1\n"
	"pfsub  (%1), %%mm0\n"
	"pfsub 8(%1), %%mm1\n"
	"movq %%mm0,  (%2)\n"
	"movq %%mm1, 8(%2)\n"
	:
	: "r" (pLeft), "r" (pRight), "r" (pOut)

	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif

}

void VPCALL CSIMD_3DNow::vector3D_AlignedScale(CVec3D* pOut, float scalar)
{
	asm(
	"movd %1, %%mm0\n"
	"punpckldq %%mm0, %%mm0\n"
	"movq %%mm0, %%mm1\n"
	"pfmul  (%0), %%mm0\n"
	"pfmul 8(%0), %%mm1\n"
	"movq %%mm0,  (%0)\n"
	"movq %%mm1, 8(%0)\n"
	:
	: "r" (pOut), "m" (scalar)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif
}

void VPCALL CSIMD_3DNow::vector3D_AlignedScaleOf(CVec3D* pOut, const CVec3D* pIn, float scalar)
{
	asm(
	"movd %2, %%mm0\n"
	"punpckldq %%mm0, %%mm0\n"
	"movq %%mm0, %%mm1\n"
	"pfmul  (%1), %%mm0\n"
	"pfmul 8(%1), %%mm1\n"
	"movq %%mm0,  (%0)\n"
	"movq %%mm1, 8(%0)\n"
	:
	: "r" (pOut), "r" (pIn), "m" (scalar)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif

}

float  VPCALL CSIMD_3DNow::vector3D_AlignedDot(const CVec3D* pSrc1, const CVec3D* pSrc2)
{

	float dummy;
	asm(

	"movq  (%2), %%mm0\n"   /* mm0 = Y1 | X1*/
	"movq 8(%2), %%mm1\n"   /* mm1 = ?? | Z1 */

	"pfmul  (%1), %%mm0\n"  /* mm0 = Y1*Y2 | X1*X2 */
	"pand _SIMDx86_float_3DNOW_NO_W_MASK, %%mm1\n" /* mm1 = ?? & 0 | Z1 & 0xFFFFFFFF */
	"pfmul 8(%1), %%mm1\n"  /* mm1 =   0.0 | Z1*Z2 */

	/* Sum the pieces */
	"pfadd %%mm0, %%mm1\n"  /* mm1 = Y1*Y2 | X1*X2 + Z1*Z2 */
	"pfacc %%mm1, %%mm1\n"  /* mm0 = ? | X1*X2 + Y1*Y2 + Z1*Z2 */

	/* Store result */
	"movd %%mm1, -4(%%esp)\n"

	/* Yikes,'femms' must be done, sorry! Even if NO_EMMS is defined, the x87 register holds the result... */
	"femms\n"
	"flds -4(%%esp)\n"
	: "=t" (dummy)
	: "r" (pSrc1), "r" (pSrc2)
	);
	return dummy;

}

float  VPCALL CSIMD_3DNow::vector3D_AlignedDot4(const float* pSrc4D1, const float* pSrc4D2)
{
	/* 3DNow!/3DNow+ Implementation */

	float dummy;
	asm(
	"movq  (%1), %%mm0\n"	/* mm0 = y1 | x1 */
	"movq 8(%1), %%mm1\n"	/* mm1 = w1 | z1 */
	"pfmul  (%2), %%mm0\n"	/* mm0 = y1*y2 | x1*x2 */
	"pfmul 8(%2), %%mm1\n"	/* mm1 = w1*w2 | z1*z2 */
	"pfadd %%mm0, %%mm1\n"	/* mm1 = y1*y2+w1*w2 | x1*x2+z1*z2 */
	"pfacc %%mm1, %%mm1\n"			/* mm1 = dot4 | dot4 */
	"movd %%mm1, -4(%%esp)\n"
	"femms\n"
	"flds -4(%%esp)\n"
	: "=t" (dummy)
	: "r" (pSrc4D1), "r" (pSrc4D2)
	);
	return dummy;


}

float  VPCALL CSIMD_3DNow::vector3D_AlignedLengthSq(const CVec3D* pVec)
{
	float dummy;
	asm(
	"movq  (%1), %%mm0\n"
	"movd 8(%1), %%mm1\n"

	/* square each component */
	"pfmul %%mm0, %%mm0\n"
	"pfmul %%mm1, %%mm1\n"

	/* Sum components */
	"pfacc %%mm0, %%mm0\n"
	"pfadd %%mm1, %%mm0\n"

	/* Move from MMX register to stack, then 'femms' (sorry!) and load x87 register */
	"movd %%mm0, -4(%%esp)\n"
	"femms\n"
	"flds -4(%%esp)\n"
	: "=t" (dummy)
	: "r" (pVec)
	);
	return dummy;

}

float  VPCALL CSIMD_3DNow::vector3D_AlignedLength(const CVec3D* pVec)
{
	float dummy;
	asm(
	"movq  (%1), %%mm0\n"
	"movd 8(%1), %%mm1\n"

	/* square each component */
	"pfmul %%mm0, %%mm0\n"
	"pfmul %%mm1, %%mm1\n"

	/* Sum components */
	"pfacc %%mm0, %%mm0\n"
	"pfadd %%mm1, %%mm0\n"

	#if defined(HIPREC)
	"movq %%mm0, %%mm3\n"
	"pfrsqrt %%mm0, %%mm1\n"        /* mm1 = x0 = rsqrt(mag*mag) */
	"movq %%mm1, %%mm2\n"           /* mm2 = x0 */
	"pfmul %%mm1, %%mm1\n"          /* mm1 = x1 = x0*x0 */
	"pfrsqit1 %%mm1, %%mm0\n"       /* mm0 = x2 = pfrsqit1(val, x1)*/
	"pfrcpit2 %%mm2, %%mm0\n"       /* mm2 = 1/sqrt(val) = pfrcpit2(x2, x0)*/
	"pfmul %%mm0, %%mm3\n"
	"movd %%mm3, -4(%%esp)\n"
	"femms\n"
	"flds -4(%%esp)\n"
	#else
	"pfrsqrt %%mm0, %%mm1\n"
	"pfmul %%mm0, %%mm1\n"
	"movd %%mm1, -4(%%esp)\n"
	"femms\n"
	"flds -4(%%esp)\n"

	#endif

	: "=t" (dummy)
	: "r" (pVec)
	);

	return dummy;


}

void VPCALL CSIMD_3DNow::vector3D_AlignedCross(CVec3D* pLeft, const CVec3D* pRight)
{
/* UGH! There has to be a much better way of accessing the memory then these fragmented ones */
	asm(
	"movq 4(%1), %%mm0\n"		/* mm0 = Rz | Ry */
	"movd  (%1), %%mm1\n"		/* mm1 =  ? | Rz */
	"movd 4(%1), %%mm5\n"		/* mm5 =  ? | Ry */
	"movd 8(%1), %%mm4\n"		/* mm4 =  ? | Rz */
	"punpckldq (%1), %%mm4\n"	/* mm4 = Rx | Rz */
	"movd 4(%0), %%mm3\n"		/* mm3 =  ? | Ly */
	"movd 8(%0), %%mm2\n"		/* mm2 =  ? | Lz */
	"punpckldq (%0), %%mm2\n"	/* mm2 = Lx | Lz */
	"movq 4(%0), %%mm6\n"		/* mm6 = Lz | Ly */
	"movd  (%0), %%mm7\n"		/* mm7 =  ? | Lx */


/*
	pLeft->x = (tmp.y * pRight->z) - (tmp.z * pRight->y);
	pLeft->y = (tmp.z * pRight->x) - (tmp.x * pRight->z);
	pLeft->z = (tmp.x * pRight->y) - (tmp.y * pRight->x);
*/

	/* SO UGLY, 10 memory accesses compared to SSE's 2 */

	/* multiply Column 1 & 2 */
	"pfmul %%mm4, %%mm6\n"
	"pfmul %%mm5, %%mm7\n"

	/* multiply Column 3 & 4 */
	"pfmul %%mm0, %%mm2\n"
	"pfmul %%mm1, %%mm3\n"

	/* Sutract... (1*2) - (3*4) */
	"pfsub %%mm2, %%mm6\n"
	"pfsub %%mm3, %%mm7\n"

	/* And store */
	"movq %%mm6,  (%0)\n"
	"movd %%mm7, 8(%0)\n"
	:
	: "r" (pLeft), "r" (pRight)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif
}

void VPCALL CSIMD_3DNow::vector3D_AlignedCrossOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight)
{
	/* UGH! There has to be a much better way of accessing the memory then these fragmented ones */
	asm(
	"movq 4(%2), %%mm0\n"		/* mm0 = Rz | Ry */
	"movd  (%2), %%mm1\n"		/* mm1 =  ? | Rz */
	"movd 4(%2), %%mm5\n"		/* mm5 =  ? | Ry */
	"movd 8(%2), %%mm4\n"		/* mm4 =  ? | Rz */
	"punpckldq (%2), %%mm4\n"	/* mm4 = Rx | Rz */
	"movd 4(%1), %%mm3\n"		/* mm3 =  ? | Ly */
	"movd 8(%1), %%mm2\n"		/* mm2 =  ? | Lz */
	"punpckldq (%1), %%mm2\n"	/* mm2 = Lx | Lz */
	"movq 4(%1), %%mm6\n"		/* mm6 = Lz | Ly */
	"movd  (%1), %%mm7\n"		/* mm7 =  ? | Lx */


/*
		pOut->x = (pLeft->y * pRight->z) - (pLeft->z * pRight->y);
		pOut->y = (pLeft->z * pRight->x) - (pLeft->x * pRight->z);
		pOut->z = (pLeft->x * pRight->y) - (pLeft->y * pRight->x);
*/
	/* SO UGLY, 10 memory accesses compared to SSE's 2 */

	/* multiply Column 1 & 2 */
	"pfmul %%mm4, %%mm6\n"	/* mm4 = Lz*Rx | Ly*Rz */
	"pfmul %%mm5, %%mm7\n"	/* mm5 = ????? | Lx*Ry */

	/* multiply Column 3 & 4 */
	"pfmul %%mm0, %%mm2\n"
	"pfmul %%mm1, %%mm3\n"

	/* Sutract... (1*2) - (3*4) */
	"pfsub %%mm2, %%mm6\n"
	"pfsub %%mm3, %%mm7\n"

	/* And store */
	"movq %%mm6,  (%0)\n"
	"movd %%mm7, 8(%0)\n"
	:
	: "r" (pOut), "r" (pLeft), "r" (pRight)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif


}
void VPCALL CSIMD_3DNow::vector3D_AlignedNormalize(CVec3D* pVec)
{

	asm (

	/* mm0-mm1 = pVec (duplicated in mm6-mm7) */
	"movq  (%0), %%mm0\n"
	"movq 8(%0), %%mm1\n"
	"pand _SIMDx86_float_3DNOW_NO_W_MASK, %%mm1\n"    /* mm1 = 0.0 | z */

	"movq %%mm0, %%mm6\n"
	"movq %%mm1, %%mm7\n"

	/* square components */
	"pfmul %%mm0, %%mm0\n"	/* y*y | x*x */
	"pfmul %%mm1, %%mm1\n"	/* 0   | z*z */

	/* Sum */
	"pfadd %%mm0, %%mm1\n"	/* 0+y*y | z*z+x*x */
	"pfacc %%mm1, %%mm1\n"	/* mag*mag | mag*mag */

	#ifdef HIPREC

	/* Get reciprocal, use Newton-Raphson HW to improve accuracy, then multiply by the high precision reciprocal. */
	"movq %%mm1, %%mm3\n"
	"pfrsqrt %%mm1, %%mm0\n"        /* mm1 = x0 = rsqrt(mag*mag) */
	"movq %%mm0, %%mm2\n"           /* mm2 = x0 */
	"pfmul %%mm0, %%mm0\n"          /* mm1 = x1 = x0*x0 */
	"pfrsqit1 %%mm0, %%mm1\n"       /* mm0 = x2 = pfrsqit1(val, x1)*/
	"pfrcpit2 %%mm2, %%mm1\n"       /* mm2 = 1/sqrt(val) = pfrcpit2(x2, x0)*/
	"pfmul %%mm1, %%mm6\n"
	"pfmul %%mm1, %%mm7\n"

	#else

	/* Get fast reciprocal and return */
	"pfrsqrt %%mm1, %%mm1\n"	/* mm1 = rcp(mag) | rcp(mag) (approx) */
	"pfmul %%mm1, %%mm6\n"		/* mm2 = y/mag | x/mag */
	"pfmul %%mm1, %%mm7\n"		/* mm3 = 0 | z/mag */
	#endif

	/* Store results */
	"movq %%mm6,  (%0)\n"
	"movq %%mm7, 8(%0)\n"
	:
	: "r" (pVec)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif

}
void VPCALL CSIMD_3DNow::vector3D_AlignedNormalizeOf(CVec3D* pOut, const CVec3D* pVec)
{

	asm (

	/* mm0-mm1 = pVec (duplicated in mm2-mm3) */
	"movq  (%0), %%mm0\n"
	"movq 8(%0), %%mm1\n"
	"pand _SIMDx86_float_3DNOW_NO_W_MASK, %%mm1\n"    /* mm1 = 0.0 | z */
	"movq %%mm0, %%mm6\n"
	"movq %%mm1, %%mm7\n"

	/* square components */
	"pfmul %%mm0, %%mm0\n"	/* y*y | x*x */
	"pfmul %%mm1, %%mm1\n"	/* 0   | z*z */

	/* Sum */
	"pfadd %%mm0, %%mm1\n"	/* 0+y*y | z*z+x*x */
	"pfacc %%mm1, %%mm1\n"	/* mag*mag | mag*mag */

	#ifdef HIPREC

	/* Get reciprocal, use Newton-Raphson HW to improve accuracy, then multiply by the high precision reciprocal. */
	"movq %%mm1, %%mm3\n"
	"pfrsqrt %%mm1, %%mm0\n"        /* mm1 = x0 = rsqrt(mag*mag) */
	"movq %%mm0, %%mm2\n"           /* mm2 = x0 */
	"pfmul %%mm0, %%mm0\n"          /* mm1 = x1 = x0*x0 */
	"pfrsqit1 %%mm0, %%mm1\n"       /* mm0 = x2 = pfrsqit1(val, x1)*/
	"pfrcpit2 %%mm2, %%mm1\n"       /* mm2 = 1/sqrt(val) = pfrcpit2(x2, x0)*/
	"pfmul %%mm1, %%mm6\n"
	"pfmul %%mm1, %%mm7\n"

	#else

	/* Get fast reciprocal and return */
	"pfrsqrt %%mm1, %%mm1\n"	/* mm1 = rcp(mag) | rcp(mag) (approx) */
	"pfmul %%mm1, %%mm6\n"		/* mm2 = y/mag | x/mag */
	"pfmul %%mm1, %%mm7\n"		/* mm3 = 0 | z/mag */
	#endif

	/* Store results */
	"movq %%mm6,  (%1)\n"
	"movq %%mm7, 8(%1)\n"
	:
	: "r" (pVec), "r" (pOut)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif

}

float  VPCALL CSIMD_3DNow::vector3D_AlignedDistance(const CVec3D* pVec1, const CVec3D* pVec2)
{

	/* TODO: Optimize me completely */
	CVec3D diff;
	vector3D_DiffOf(&diff, pVec1, pVec2);
	return vector3D_Length(&diff);

}

//========CMat4D==============================

void VPCALL CSIMD_3DNow::mat4D_Sum(CMat4D* Out, const CMat4D* In)
{

	asm(
	"movq   (%0), %%mm0\n"
	"movq  8(%0), %%mm1\n"
	"movq 16(%0), %%mm2\n"
	"movq 24(%0), %%mm3\n"
	"movq 32(%0), %%mm4\n"
	"movq 40(%0), %%mm5\n"
	"movq 48(%0), %%mm6\n"
	"movq 56(%0), %%mm7\n"
	"pfadd   (%1), %%mm0\n"
	"pfadd  8(%1), %%mm1\n"
	"pfadd 16(%1), %%mm2\n"
	"pfadd 24(%1), %%mm3\n"
	"pfadd 32(%1), %%mm4\n"
	"pfadd 40(%1), %%mm5\n"
	"pfadd 48(%1), %%mm6\n"
	"pfadd 56(%1), %%mm7\n"
	"movq %%mm0,   (%0)\n"
	"movq %%mm1,  8(%0)\n"
	"movq %%mm2, 16(%0)\n"
	"movq %%mm3, 24(%0)\n"
	"movq %%mm4, 32(%0)\n"
	"movq %%mm5, 40(%0)\n"
	"movq %%mm6, 48(%0)\n"
	"movq %%mm7, 56(%0)\n"
	:
	: "r" (Out), "r" (In)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif



}


void VPCALL CSIMD_3DNow::mat4D_SumOf(CMat4D* Out, const CMat4D* In1, const CMat4D* In2)
{

	

	asm(

	/* Get matrix */
	"movq   (%0), %%mm0\n"
	"movq  8(%0), %%mm1\n"
	"movq 16(%0), %%mm2\n"
	"movq 24(%0), %%mm3\n"
	"movq 32(%0), %%mm4\n"
	"movq 40(%0), %%mm5\n"
	"movq 48(%0), %%mm6\n"
	"movq 56(%0), %%mm7\n"

	/* add with other */
	"pfadd   (%1), %%mm0\n"
	"pfadd  8(%1), %%mm1\n"
	"pfadd 16(%1), %%mm2\n"
	"pfadd 24(%1), %%mm3\n"
	"pfadd 32(%1), %%mm4\n"
	"pfadd 40(%1), %%mm5\n"
	"pfadd 48(%1), %%mm6\n"
	"pfadd 56(%1), %%mm7\n"

	/* Write out results */
	"movq %%mm0,   (%2)\n"
	"movq %%mm1,  8(%2)\n"
	"movq %%mm2, 16(%2)\n"
	"movq %%mm3, 24(%2)\n"
	"movq %%mm4, 32(%2)\n"
	"movq %%mm5, 40(%2)\n"
	"movq %%mm6, 48(%2)\n"
	"movq %%mm7, 56(%2)\n"
	:
	: "r" (In1), "r" (In2), "r" (Out)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif


}

void VPCALL CSIMD_3DNow::mat4D_Diff(CMat4D* Out, const CMat4D* In)
{
	

	asm(
	"movq   (%0), %%mm0\n"
	"movq  8(%0), %%mm1\n"
	"movq 16(%0), %%mm2\n"
	"movq 24(%0), %%mm3\n"
	"movq 32(%0), %%mm4\n"
	"movq 40(%0), %%mm5\n"
	"movq 48(%0), %%mm6\n"
	"movq 56(%0), %%mm7\n"
	"pfsub   (%1), %%mm0\n"
	"pfsub  8(%1), %%mm1\n"
	"pfsub 16(%1), %%mm2\n"
	"pfsub 24(%1), %%mm3\n"
	"pfsub 32(%1), %%mm4\n"
	"pfsub 40(%1), %%mm5\n"
	"pfsub 48(%1), %%mm6\n"
	"pfsub 56(%1), %%mm7\n"
	"movq %%mm0,   (%0)\n"
	"movq %%mm1,  8(%0)\n"
	"movq %%mm2, 16(%0)\n"
	"movq %%mm3, 24(%0)\n"
	"movq %%mm4, 32(%0)\n"
	"movq %%mm5, 40(%0)\n"
	"movq %%mm6, 48(%0)\n"
	"movq %%mm7, 56(%0)\n"
	:
	: "r" (Out), "r" (In)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif

}

void VPCALL CSIMD_3DNow::mat4D_DiffOf(CMat4D* Out, const CMat4D* In1, const CMat4D* In2)
{

	asm(

	/* Get 'In1'*/
	"movq   (%0), %%mm0\n"
	"movq  8(%0), %%mm1\n"
	"movq 16(%0), %%mm2\n"
	"movq 24(%0), %%mm3\n"
	"movq 32(%0), %%mm4\n"
	"movq 40(%0), %%mm5\n"
	"movq 48(%0), %%mm6\n"
	"movq 56(%0), %%mm7\n"

	/* Subtract 'In2' */
	"pfsub   (%1), %%mm0\n"
	"pfsub  8(%1), %%mm1\n"
	"pfsub 16(%1), %%mm2\n"
	"pfsub 24(%1), %%mm3\n"
	"pfsub 32(%1), %%mm4\n"
	"pfsub 40(%1), %%mm5\n"
	"pfsub 48(%1), %%mm6\n"
	"pfsub 56(%1), %%mm7\n"

	/* Store in 'Out' */
	"movq %%mm0,   (%2)\n"
	"movq %%mm1,  8(%2)\n"
	"movq %%mm2, 16(%2)\n"
	"movq %%mm3, 24(%2)\n"
	"movq %%mm4, 32(%2)\n"
	"movq %%mm5, 40(%2)\n"
	"movq %%mm6, 48(%2)\n"
	"movq %%mm7, 56(%2)\n"

	:
	: "r" (In1), "r" (In2), "r" (Out)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif


}

void VPCALL CSIMD_3DNow::mat4D_Scale(CMat4D* mtx, float scalar)
{
	
	asm(

	/* set mm0 = scalar | scalar */
	"movd %0, %%mm0\n"
	"punpckldq %%mm0, %%mm0\n"

	/* Duplicate dual-scalar register */
	"movq %%mm0, %%mm1\n"
	"movq %%mm0, %%mm2\n"
	"movq %%mm0, %%mm3\n"
	"movq %%mm0, %%mm4\n"
	"movq %%mm0, %%mm5\n"
	"movq %%mm0, %%mm6\n"
	"movq %%mm0, %%mm7\n"

	/* multiply matrix by scalar */
	"pfmul   (%1), %%mm0\n"
	"pfmul  8(%1), %%mm1\n"
	"pfmul 16(%1), %%mm2\n"
	"pfmul 24(%1), %%mm3\n"
	"pfmul 32(%1), %%mm4\n"
	"pfmul 40(%1), %%mm5\n"
	"pfmul 48(%1), %%mm6\n"
	"pfmul 56(%1), %%mm7\n"

	/* Store results */
	"movq %%mm0,   (%1)\n"
	"movq %%mm1,  8(%1)\n"
	"movq %%mm2, 16(%1)\n"
	"movq %%mm3, 24(%1)\n"
	"movq %%mm4, 32(%1)\n"
	"movq %%mm5, 40(%1)\n"
	"movq %%mm6, 48(%1)\n"
	"movq %%mm7, 56(%1)\n"
	:
	: "m" (scalar), "r" (mtx)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif

}

void VPCALL CSIMD_3DNow::mat4D_ScaleOf(CMat4D* pOut, const CMat4D* pIn, float scalar)
{
	
	asm(
	"movd %0, %%mm0\n"
	"punpckldq %%mm0, %%mm0\n"
	"movq %%mm0, %%mm1\n"
	"movq %%mm0, %%mm2\n"
	"movq %%mm0, %%mm3\n"
	"movq %%mm0, %%mm4\n"
	"movq %%mm0, %%mm5\n"
	"movq %%mm0, %%mm6\n"
	"movq %%mm0, %%mm7\n"
	"pfmul   (%1), %%mm0\n"
	"pfmul  8(%1), %%mm1\n"
	"pfmul 16(%1), %%mm2\n"
	"pfmul 24(%1), %%mm3\n"
	"pfmul 32(%1), %%mm4\n"
	"pfmul 40(%1), %%mm5\n"
	"pfmul 48(%1), %%mm6\n"
	"pfmul 56(%1), %%mm7\n"
	"movq %%mm0,   (%2)\n"
	"movq %%mm1,  8(%2)\n"
	"movq %%mm2, 16(%2)\n"
	"movq %%mm3, 24(%2)\n"
	"movq %%mm4, 32(%2)\n"
	"movq %%mm5, 40(%2)\n"
	"movq %%mm6, 48(%2)\n"
	"movq %%mm7, 56(%2)\n"
	:
	: "m" (scalar), "r" (pIn), "r" (pOut)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif

}

void VPCALL CSIMD_3DNow::mat4D_Multiply(CMat4D* pLeft, const CMat4D* pRight)
{
	

	asm(
	"movl $4, %%ecx\n"
	"MatrixMultiply_3DNow_Loop%=:\n"
	"movd   (%0), %%mm0\n"			/* mm0 = ? | x */
	"movd  4(%0), %%mm2\n"			/* mm2 = ? | y */
	"movd  8(%0), %%mm4\n"			/* mm4 = ? | z */
	"movd 12(%0), %%mm6\n"			/* mm6 = ? | w */
	"prefetchw 16(%0)\n"			/* prefetch_for_writing(m[4]...m[7]);  */
	"punpckldq %%mm0, %%mm0\n"		/* mm0 = x | x */
	"punpckldq %%mm2, %%mm2\n"		/* mm2 = y | y */
	"punpckldq %%mm4, %%mm4\n"		/* mm4 = z | z */
	"punpckldq %%mm6, %%mm6\n"		/* mm6 = w | w */
	"movq %%mm0, %%mm1\n"			/* mm1 = x | x */
	"movq %%mm2, %%mm3\n"			/* mm3 = y | y */
	"movq %%mm4, %%mm5\n"			/* mm5 = z | z */
	"movq %%mm6, %%mm7\n"			/* mm7 = w | w */
	"pfmul   (%1), %%mm0\n"			/* mm0 = x*m[1]  | x*m[0] */
	"pfmul  8(%1), %%mm1\n"			/* mm1 = x*m[3]  | x*m[2] */
	"pfmul 16(%1), %%mm2\n" 		/* mm2 = y*m[5]  | y*m[4] */
	"pfmul 24(%1), %%mm3\n" 		/* mm3 = y*m[7]  | y*m[6] */
	"pfmul 32(%1), %%mm4\n" 		/* mm4 = z*m[9]  | z*m[8] */
	"pfmul 40(%1), %%mm5\n" 		/* mm5 = z*m[11] | z*m[10] */
	"pfmul 48(%1), %%mm6\n" 		/* mm6 = w*m[13] | w*m[12] */
	"pfmul 56(%1), %%mm7\n" 		/* mm7 = w*m[15] | w*m[14] */
	"pfadd %%mm0, %%mm2\n"
	"pfadd %%mm1, %%mm3\n"
	"pfadd %%mm4, %%mm6\n"
	"pfadd %%mm5, %%mm7\n"
	"pfadd %%mm2, %%mm6\n"
	"pfadd %%mm3, %%mm7\n"
	"movq %%mm6,  (%0)\n"
	"movq %%mm7, 8(%0)\n"
	"addl $16, %0\n"
	"loop MatrixMultiply_3DNow_Loop%=\n"
	:
	: "r" (pLeft), "r" (pRight)
	: "%ecx"
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif


}

void VPCALL CSIMD_3DNow::mat4D_MultiplyOf(CMat4D* pOut, const CMat4D* pLeft, const CMat4D* pRight)
{
	
	asm(
	"movl $4, %%ecx\n"
	"MatrixMultiplyOf_3DNow_Loop%=:\n"
	"movd   (%0), %%mm0\n"			/* mm0 = ? | x */
	"movd  4(%0), %%mm2\n"			/* mm2 = ? | y */
	"movd  8(%0), %%mm4\n"			/* mm4 = ? | z */
	"movd 12(%0), %%mm6\n"			/* mm6 = ? | w */
	"prefetch 16(%0)\n"				/* prefetch_for_reading(in+4);  */
	"prefetchw 16(%2)\n"			/* prefetch_for_writing(out+4); */
	"punpckldq %%mm0, %%mm0\n"		/* mm0 = x | x */
	"punpckldq %%mm2, %%mm2\n"		/* mm2 = y | y */
	"punpckldq %%mm4, %%mm4\n"		/* mm4 = z | z */
	"punpckldq %%mm6, %%mm6\n"		/* mm6 = w | w */
	"movq %%mm0, %%mm1\n"			/* mm1 = x | x */
	"movq %%mm2, %%mm3\n"			/* mm3 = y | y */
	"movq %%mm4, %%mm5\n"			/* mm5 = z | z */
	"movq %%mm6, %%mm7\n"			/* mm7 = w | w */
	"pfmul   (%1), %%mm0\n"			/* mm0 = x*m[1]  | x*m[0] */
	"pfmul  8(%1), %%mm1\n"			/* mm1 = x*m[3]  | x*m[2] */
	"pfmul 16(%1), %%mm2\n" 		/* mm2 = y*m[5]  | y*m[4] */
	"pfmul 24(%1), %%mm3\n" 		/* mm3 = y*m[7]  | y*m[6] */
	"pfmul 32(%1), %%mm4\n" 		/* mm4 = z*m[9]  | z*m[8] */
	"pfmul 40(%1), %%mm5\n" 		/* mm5 = z*m[11] | z*m[10] */
	"pfmul 48(%1), %%mm6\n" 		/* mm6 = w*m[13] | w*m[12] */
	"pfmul 56(%1), %%mm7\n" 		/* mm7 = w*m[15] | w*m[14] */
	"pfadd %%mm0, %%mm2\n"
	"pfadd %%mm1, %%mm3\n"
	"pfadd %%mm4, %%mm6\n"
	"pfadd %%mm5, %%mm7\n"
	"pfadd %%mm2, %%mm6\n"
	"pfadd %%mm3, %%mm7\n"
	"movq %%mm6,  (%2)\n"
	"movq %%mm7, 8(%2)\n"
	"addl $16, %0\n"
	"addl $16, %2\n"
	"loop MatrixMultiplyOf_3DNow_Loop%=\n"
	:
	: "r" (pLeft), "r" (pRight), "r" (pOut)
	: "%ecx"
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif


}


void VPCALL CSIMD_3DNow::mat4D_VectorMultiply(CVec3D* pVec, const CMat4D* pMat)
{
	/*
		Does a normal vector x matrix, but since w = 1.0f, the final column of the matrix is
		merely added.
	*/

	

	asm(
	"movd   (%0), %%mm0\n"			/* mm0 = ? | x */
	"movd  4(%0), %%mm2\n"			/* mm2 = ? | y */
	"movd  8(%0), %%mm4\n"			/* mm4 = ? | z */
	"punpckldq %%mm0, %%mm0\n"		/* mm0 = x | x */
	"punpckldq %%mm2, %%mm2\n"		/* mm2 = y | y */
	"punpckldq %%mm4, %%mm4\n"		/* mm4 = z | z */
	"movq %%mm0, %%mm1\n"			/* mm1 = x | x */
	"movq %%mm2, %%mm3\n"			/* mm3 = y | y */
	"movq %%mm4, %%mm5\n"			/* mm5 = z | z */
	"pfmul   (%1), %%mm0\n"			/* mm0 = x*m[1]  | x*m[0] */
	"pfmul  8(%1), %%mm1\n"			/* mm1 = x*m[3]  | x*m[2] */
	"pfmul 16(%1), %%mm2\n" 		/* mm2 = y*m[5]  | y*m[4] */
	"pfmul 24(%1), %%mm3\n" 		/* mm3 = y*m[7]  | y*m[6] */
	"pfmul 32(%1), %%mm4\n" 		/* mm4 = z*m[9]  | z*m[8] */
	"pfmul 40(%1), %%mm5\n" 		/* mm5 = z*m[11] | z*m[10] */
	"pfadd 48(%1), %%mm4\n"         /* mm4 = z*m[9] + m[13]  | z*m[8] + m[12] */
	"pfadd 56(%1), %%mm5\n"			/* mm5 = z*m[11] + m[15] | z*m[10] + m[14] */

	/* Sum it... (not displaying register contents.... */

	/* Combine X and Y column into mm2-3 */
	"pfadd %%mm0, %%mm2\n"
	"pfadd %%mm1, %%mm3\n"

	/* Combine with ZW column into mm2-3 */
	"pfadd %%mm4, %%mm2\n"
	"pfadd %%mm5, %%mm3\n"

	"movq %%mm2,  (%0)\n"
	"movq %%mm3, 8(%0)\n"
	:
	: "r" (pVec), "r" (pMat)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif


}

void VPCALL CSIMD_3DNow::mat4D_VectorMultiplyOf(CVec3D* pOut, const CVec3D* pIn, const CMat4D* pMat)
{
	/*
		Does a normal vector x matrix, but since w = 1.0f, the final column of the matrix is
		merely added.
	*/

	
	asm(
	"movd   (%0), %%mm0\n"			/* mm0 = ? | x */
	"movd  4(%0), %%mm2\n"			/* mm2 = ? | y */
	"movd  8(%0), %%mm4\n"			/* mm4 = ? | z */
	"punpckldq %%mm0, %%mm0\n"		/* mm0 = x | x */
	"punpckldq %%mm2, %%mm2\n"		/* mm2 = y | y */
	"punpckldq %%mm4, %%mm4\n"		/* mm4 = z | z */
	"movq %%mm0, %%mm1\n"			/* mm1 = x | x */
	"movq %%mm2, %%mm3\n"			/* mm3 = y | y */
	"movq %%mm4, %%mm5\n"			/* mm5 = z | z */
	"pfmul   (%1), %%mm0\n"			/* mm0 = x*m[1]  | x*m[0] */
	"pfmul  8(%1), %%mm1\n"			/* mm1 = x*m[3]  | x*m[2] */
	"pfmul 16(%1), %%mm2\n" 		/* mm2 = y*m[5]  | y*m[4] */
	"pfmul 24(%1), %%mm3\n" 		/* mm3 = y*m[7]  | y*m[6] */
	"pfmul 32(%1), %%mm4\n" 		/* mm4 = z*m[9]  | z*m[8] */
	"pfmul 40(%1), %%mm5\n" 		/* mm5 = z*m[11] | z*m[10] */
	"pfadd 48(%1), %%mm4\n"         /* mm4 = z*m[9] + m[13]  | z*m[8] + m[12] */
	"pfadd 56(%1), %%mm5\n"			/* mm5 = z*m[11] + m[15] | z*m[10] + m[14] */

	/* Sum it... (not displaying register contents.... */

	/* Combine X and Y column into mm2-3 */
	"pfadd %%mm0, %%mm2\n"
	"pfadd %%mm1, %%mm3\n"

	/* Combine with ZW column into mm2-3 */
	"pfadd %%mm4, %%mm2\n"
	"pfadd %%mm5, %%mm3\n"

	"movq %%mm2,  (%2)\n"
	"movq %%mm3, 8(%2)\n"
	:
	: "r" (pIn), "r" (pMat), "r" (pOut)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif


}

void VPCALL CSIMD_3DNow::mat4D_VectorMultiply(CVec4D* pOut4D, const CMat4D* pMat)
{
	
	asm(
	"movd   (%0), %%mm0\n"			/* mm0 = ? | x */
	"movd  4(%0), %%mm2\n"			/* mm2 = ? | y */
	"movd  8(%0), %%mm4\n"			/* mm4 = ? | z */
	"movd 12(%0), %%mm6\n"			/* mm6 = ? | w */
	"punpckldq %%mm0, %%mm0\n"		/* mm0 = x | x */
	"punpckldq %%mm2, %%mm2\n"		/* mm2 = y | y */
	"punpckldq %%mm4, %%mm4\n"		/* mm4 = z | z */
	"punpckldq %%mm6, %%mm6\n"		/* mm6 = w | w */
	"movq %%mm0, %%mm1\n"			/* mm1 = x | x */
	"movq %%mm2, %%mm3\n"			/* mm3 = y | y */
	"movq %%mm4, %%mm5\n"			/* mm5 = z | z */
	"movq %%mm6, %%mm7\n"			/* mm7 = w | w */
	"pfmul   (%1), %%mm0\n"			/* mm0 = x*m[1]  | x*m[0] */
	"pfmul  8(%1), %%mm1\n"			/* mm1 = x*m[3]  | x*m[2] */
	"pfmul 16(%1), %%mm2\n" 		/* mm2 = y*m[5]  | y*m[4] */
	"pfmul 24(%1), %%mm3\n" 		/* mm3 = y*m[7]  | y*m[6] */
	"pfmul 32(%1), %%mm4\n" 		/* mm4 = z*m[9]  | z*m[8] */
	"pfmul 40(%1), %%mm5\n" 		/* mm5 = z*m[11] | z*m[10] */
	"pfmul 48(%1), %%mm6\n" 		/* mm6 = w*m[13] | w*m[12] */
	"pfmul 56(%1), %%mm7\n" 		/* mm7 = w*m[15] | w*m[14] */

	/* Sum it... (not displaying register contents.... */
	"pfadd %%mm0, %%mm2\n"
	"pfadd %%mm1, %%mm3\n"
	"pfadd %%mm4, %%mm6\n"
	"pfadd %%mm5, %%mm7\n"
	"pfadd %%mm2, %%mm6\n"
	"pfadd %%mm3, %%mm7\n"
	"movq %%mm6,  (%0)\n"
	"movq %%mm7, 8(%0)\n"
	:
	: "r" (pOut4D), "r" (pMat)
	);

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif



}


void VPCALL CSIMD_3DNow::mat4D_VectorMultiplyOf(CVec4D* pOut4D, const CVec4D* pIn4D, const CMat4D* pMat)
{
#if 0
	asm(
	"movd   (%0), %%mm0\n"			/* mm0 = ? | x */
	"movd  4(%0), %%mm2\n"			/* mm2 = ? | y */
	"movd  8(%0), %%mm4\n"			/* mm4 = ? | z */
	"movd 12(%0), %%mm6\n"			/* mm6 = ? | w */
	"punpckldq %%mm0, %%mm0\n"		/* mm0 = x | x */
	"punpckldq %%mm2, %%mm2\n"		/* mm2 = y | y */
	"punpckldq %%mm4, %%mm4\n"		/* mm4 = z | z */
	"punpckldq %%mm6, %%mm6\n"		/* mm6 = w | w */
	"movq %%mm0, %%mm1\n"			/* mm1 = x | x */
	"movq %%mm2, %%mm3\n"			/* mm3 = y | y */
	"movq %%mm4, %%mm5\n"			/* mm5 = z | z */
	"movq %%mm6, %%mm7\n"			/* mm7 = w | w */
	"pfmul   (%1), %%mm0\n"			/* mm0 = x*m[1]  | x*m[0] */
	"pfmul  8(%1), %%mm1\n"			/* mm1 = x*m[3]  | x*m[2] */
	"pfmul 16(%1), %%mm2\n" 		/* mm2 = y*m[5]  | y*m[4] */
	"pfmul 24(%1), %%mm3\n" 		/* mm3 = y*m[7]  | y*m[6] */
	"pfmul 32(%1), %%mm4\n" 		/* mm4 = z*m[9]  | z*m[8] */
	"pfmul 40(%1), %%mm5\n" 		/* mm5 = z*m[11] | z*m[10] */
	"pfmul 48(%1), %%mm6\n" 		/* mm6 = w*m[13] | w*m[12] */
	"pfmul 56(%1), %%mm7\n" 		/* mm7 = w*m[15] | w*m[14] */
	"pfadd %%mm0, %%mm2\n"
	"pfadd %%mm1, %%mm3\n"
	"pfadd %%mm4, %%mm6\n"
	"pfadd %%mm5, %%mm7\n"
	"pfadd %%mm2, %%mm6\n"
	"pfadd %%mm3, %%mm7\n"
	"movq %%mm6,  (%2)\n"
	"movq %%mm7, 8(%2)\n"
	:
	: "r" (pIn4D), "r" (pMat), "r" (pOut4D)
	   );
#else
	asm(
	"addl $16, %2\n" /* res++ */
	"movq  (%1), %%mm0\n"		/* mm0 = y | x */
	"movq 8(%1), %%mm1\n"		/* mm1 = w | z */
	"movq %%mm0, %%mm2\n"		/* mm2 = y | x */
	"movq (%0), %%mm3\n"		/* mm3 = m[1] | m[0] */
	"punpckldq %%mm0, %%mm0\n"	/* mm0 = x | x */
	"movq 16(%0), %%mm4\n"		/* mm4 = m[5] | m[4] */
	"pfmul %%mm0, %%mm3\n"		/* mm3 = x*m[1] | x*m[0] */
	"punpckhdq %%mm2, %%mm2\n"	/* mm2 = y | y */
	"pfmul %%mm2, %%mm4\n"		/* mm4 = y*m[5] | y*m[6] */
	"movq  8(%0), %%mm5\n"		/* mm5 = m[3] | m[2] */
	"movq 24(%0), %%mm7\n"		/* mm7 = m[7] | m[6] */
	"movq %%mm1, %%mm6\n"		/* mm6 = w | z */
	"pfmul %%mm0, %%mm5\n"		/* mm5 = m[3]*x | m[2]*x */

	/* Finished multiplying x's */

	"movq 32(%0), %%mm0\n"		/* mm0 = m[9] | m[8] */
	"punpckldq %%mm1, %%mm1\n"	/* mm1 = z | z */
	"pfmul %%mm2, %%mm7\n"		/* mm7 = m[7]*y | m[6]*y */

	/* Finished multiplying y's */

	"movq 40(%0), %%mm2\n"		/* mm2 = m[10] | m[9] */
	"pfmul %%mm1, %%mm0\n"		/* mm0 = m[9]*z | m[8]*z */
	"pfadd %%mm4, %%mm3\n"		/* mm3 = x*m[1] + y*m[5] | x*m[0] + y*m[4] */
	"movq 48(%0), %%mm4\n"		/* mm4 = m[13] | m[12] */
	"pfmul %%mm1, %%mm2\n"		/* mm2 = m[11]*z | m[10]*z */

	/* Finished multiplying z's */

	"pfadd %%mm7, %%mm5\n"		/* mm5 = m[3]*x + m[7]*y | m[2]*x + m[6]*y */

	"movq 56(%0), %%mm1\n"		/* mm1 = m[15] | m[14] */
	"punpckhdq %%mm6, %%mm6\n"	/* mm6 = w | w */
	"pfadd %%mm0, %%mm3\n"		/* mm3 = x*m[1] + y*m[5] + m[9]*z | x*m[0] + y*m[4] + m[8]*z */

	"pfmul %%mm6, %%mm4\n"		/* mm4 = w*m[13] | w*m[12] */
	"pfadd %%mm2, %%mm5\n"		/* mm5 = m[3]*x + m[7]*y + m[11]*z | m[2]*x + m[6]*y + m[10]*z */
	"pfmul %%mm6, %%mm1\n"		/* mm1 = w*m[15] | w*m[14] */
	"pfadd %%mm4, %%mm3\n"		/* mm3 = x*m[1] + y*m[5] + m[9]*z + w*m[13] | x*m[0] + y*m[4] + m[8]*z + w*m[12] */

	"movq %%mm3, -16(%2)\n"
	"pfadd %%mm1, %%mm5\n"		/* mm5 = m[3]*x + m[7]*y + m[11]*z + m[15]*w | m[2]*x + m[6]*y + m[10]*z + m[14]*w */
	"movq %%mm5, -8(%2)\n"
	:
	: "r" (pMat), "r" (pIn4D), "r" (pOut4D)
	   );
#endif

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif

}

//=====================Quaternion===================================

void  VPCALL CSIMD_3DNow:: quaternion_Normalize(CQuaternion* pQuat)
{

	asm(
	"movq  (%0), %%mm0\n"		/* mm0 = y | x */
	"movq 8(%0), %%mm1\n"		/* mm1 = w | z */
	"movq %%mm0, %%mm2\n"		/* mm2 = y | x */
	"movq %%mm1, %%mm3\n"		/* mm3 = w | z */
	"pfmul %%mm0, %%mm0\n"		/* y*y | x*x */
	"pfmul %%mm1, %%mm1\n"		/* w*w | z*z */
	"pfadd %%mm1, %%mm0\n"		/* y*y + w*w | z*z + x*x */
	"pfacc %%mm0, %%mm0\n"		/* magsq | magsq */
	"pfrsqrt %%mm0, %%mm1\n"	/* mm1 = rcp(mag) | rcp(mag) (approx) */

	#ifdef HIPREC /* High precision */
	"movq %%mm1, %%mm4\n"		/* mm4 = mm1 = x0 = rcp(mag) */
	"pfmul %%mm1, %%mm1\n"		/* mm1 = x1 = x0*x0 */
	"pfrsqit1 %%mm1, %%mm0\n"	/* mm1 = x2 = pfrsqit1(val, x1) */
	"pfrcpit2 %%mm4, %%mm0\n"	/* mm2 = 1/sqrt(mag*mag) = pfrcpit2(x2, x0)*/
	"pfmul %%mm0, %%mm2\n"
	"pfmul %%mm0, %%mm3\n"

	#else /* Low precision */
	"pfmul %%mm1, %%mm2\n"
	"pfmul %%mm1, %%mm3\n"
	#endif
	"movq %%mm2,  (%0)\n"
	"movq %%mm3, 8(%0)\n"
	:
	: "r" (pQuat)
	);

	/* Execute FEMMS if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif


}

void  VPCALL CSIMD_3DNow:: quaternion_NormalizeOf(CQuaternion* pOut, const CQuaternion* pQuat)
{

	asm(
	"movq  (%0), %%mm0\n"		/* mm0 = y | x */
	"movq 8(%0), %%mm1\n"		/* mm1 = w | z */
	"movq %%mm0, %%mm2\n"		/* mm2 = y | x */
	"movq %%mm1, %%mm3\n"		/* mm3 = w | z */
	"pfmul %%mm0, %%mm0\n"		/* y*y | x*x */
	"pfmul %%mm1, %%mm1\n"		/* w*w | z*z */
	"pfadd %%mm1, %%mm0\n"		/* y*y + w*w | z*z + x*x */
	"pfacc %%mm0, %%mm0\n"		/* magsq | magsq */
	"pfrsqrt %%mm0, %%mm1\n"	/* mm1 = rcp(mag) | rcp(mag) (approx) */

	#ifdef HIPREC /* High precision */
	"movq %%mm1, %%mm4\n"		/* mm4 = mm1 = x0 = rcp(mag) */
	"pfmul %%mm1, %%mm1\n"		/* mm1 = x1 = x0*x0 */
	"pfrsqit1 %%mm1, %%mm0\n"	/* mm1 = x2 = pfrsqit1(val, x1) */
	"pfrcpit2 %%mm4, %%mm0\n"	/* mm2 = 1/sqrt(mag*mag) = pfrcpit2(x2, x0)*/
	"pfmul %%mm0, %%mm2\n"
	"pfmul %%mm0, %%mm3\n"

	#else /* Low precision */
	"pfmul %%mm1, %%mm2\n"
	"pfmul %%mm1, %%mm3\n"
	#endif
	"movq %%mm2,  (%1)\n"
	"movq %%mm3, 8(%1)\n"
	:
	: "r" (pQuat), "r" (pOut)
	);

	/* Execute FEMMS if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif



}

void  VPCALL CSIMD_3DNow:: quaternion_Multiply(CQuaternion* pLeft, const CQuaternion* pRight)
{
	/* 3DNow+ Implementation (uses pswapd) */
//	#if defined(USE_3DNOW) && (USE_3DNOW >= 2)
	asm(
	/* set up mm0 and mm1 as second column RIGHT */
	"pswapd  (%0), %%mm1\n"		/* mm1 = Rx | Ry */
	"pswapd 8(%0), %%mm0\n"		/* mm0 = Rz | Rw */

	/* Get the scalars */
	"movq  (%1), %%mm2\n"		/* mm2 = Ly | Lx */
	"movq 8(%1), %%mm7\n"		/* mm7 = Lw | Lz */

	/* set up mm4 and mm5 as third column RIGHT */
	"pswapd %%mm0, %%mm4\n"		/* mm4 = Rw | Rz */
	"pswapd %%mm1, %%mm5\n"		/* mm5 = Ry | Rx */

	/* Copy Rz|Rw value from mm0 into mm6 for later */
	"movq %%mm0, %%mm6\n"

	/* set up mm2 and mm3 as scalars */
	"movq %%mm2, %%mm3\n"		/* mm3 = mm2 = Ly | Lx */
	"punpckldq %%mm2, %%mm2\n"	/* mm2 = Lx | Lx */
	"punpckhdq %%mm3, %%mm3\n"	/* mm3 = Ly | Ly */

	/* multiply Column 2 (X's, frees register mm2) */
	"pfmul %%mm2, %%mm0\n"
	"pfmul %%mm2, %%mm1\n"


	/* Change signs appropriately (+, -, +, -) */
	"pxor _SIMDx86_float_POSNEG, %%mm0\n"
	"pxor _SIMDx86_float_POSNEG, %%mm1\n"

	/* Move Rx|Ry value from mm5 to mm2 */
	"pswapd %%mm5, %%mm2\n"


	/* multiply Colum 3 (Y's, frees register mm3) */
	"pfmul %%mm3, %%mm4\n"
	"pfmul %%mm3, %%mm5\n"

	/* Combine columns 2 and 3 (frees registers mm4 and mm5) */
	"pfadd %%mm4, %%mm0\n"
	"pfsub %%mm5, %%mm1\n"

	/*
		Current register status:
		========================
		mm0 = C2+C3
		mm1 = C2+C3
		mm2 = Rx | Ry
		mm3 = free
		mm4 = free
		mm5 = free
		mm6 = Rz | Rw
		mm7 = Lw | Lz

		OK! This is current:
		-------------------
		Column 4 = mm2 | mm6


		This is needed:
		--------------
		Column 1 = mm3{swap(mm2)} | mm4{swap(mm6)}
		mm5 = Lz | Lz
		mm7 = Lw | Lw

		After that is done, just multiply, invert signs as needed, then add and store

	*/

	/* set up registeres */
	"pswapd %%mm2, %%mm3\n"		/* mm3 = Ry | Rx */
	"movq %%mm7, %%mm5\n"		/* mm5 = Lw | Lz */
	"pswapd %%mm6, %%mm4\n"		/* mm4 = Rw | Rz */
	"punpckhdq %%mm7, %%mm7\n"	/* mm7 = Lw | Lw */
	"punpckldq %%mm5, %%mm5\n"	/* mm5 = Lz | Lz */


	/* All registers set up! Lets bust out some fp-ops */

	/* multiply w by first column */
	"pfmul %%mm7, %%mm3\n"
	"pfmul %%mm7, %%mm4\n"

	/* multiply z by fourth column (remember 4th = { mm2 | mm6 } )*/
	"pfmul %%mm5, %%mm2\n"
	"pfmul %%mm5, %%mm6\n"

	/* add results of C1 to C2+C3 */
	"pfadd %%mm3, %%mm0\n"
	"pfadd %%mm4, %%mm1\n"

	/* Change signs of fourth column to -, +, +, - */
	"pxor _SIMDx86_float_NEGPOS, %%mm2\n"
	"pxor _SIMDx86_float_POSNEG, %%mm6\n"

	/* Finally, get C1+C2+C3+C4, or the multiplication of the quaternions */
	"pfadd %%mm2, %%mm0\n"
	"pfadd %%mm6, %%mm1\n"

	/* Store */
	"movq %%mm0,  (%1)\n"
	"movq %%mm1, 8(%1)\n"


	:
	: "r" (pRight), "r" (pLeft)
	);

	/* 36 3DNow! Instructions (including movements) vs 70 x87 Instructions (-O3 -fomit-frame-pointer) */

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif

}

void  VPCALL CSIMD_3DNow:: quaternion_MultiplyOf(CQuaternion* pOut, const CQuaternion* pLeft, const CQuaternion* pRight)
{
	//#if defined(USE_3DNOW) && (USE_3DNOW >= 2) /* 3DNow+ Implementation (uses pswapd) */
	asm(
	/* set up mm0 and mm1 as second column RIGHT */
	"pswapd  (%0), %%mm1\n"		/* mm1 = Rx | Ry */
	"pswapd 8(%0), %%mm0\n"		/* mm0 = Rz | Rw */

	/* Get the scalars */
	"movq  (%1), %%mm2\n"		/* mm2 = Ly | Lx */
	"movq 8(%1), %%mm7\n"		/* mm7 = Lw | Lz */

	/* set up mm4 and mm5 as third column RIGHT */
	"pswapd %%mm0, %%mm4\n"		/* mm4 = Rw | Rz */
	"pswapd %%mm1, %%mm5\n"		/* mm5 = Ry | Rx */

	/* Copy Rz|Rw value from mm0 into mm6 for later */
	"movq %%mm0, %%mm6\n"

	/* set up mm2 and mm3 as scalars */
	"movq %%mm2, %%mm3\n"		/* mm3 = mm2 = Ly | Lx */
	"punpckldq %%mm2, %%mm2\n"	/* mm2 = Lx | Lx */
	"punpckhdq %%mm3, %%mm3\n"	/* mm3 = Ly | Ly */

	/* multiply Column 2 (X's, frees register mm2) */
	"pfmul %%mm2, %%mm0\n"
	"pfmul %%mm2, %%mm1\n"


	/* Change signs appropriately (+, -, +, -) */
	"pxor _SIMDx86_float_POSNEG, %%mm0\n"
	"pxor _SIMDx86_float_POSNEG, %%mm1\n"

	/* Move Rx|Ry value from mm5 to mm2 */
	"pswapd %%mm5, %%mm2\n"


	/* multiply Colum 3 (Y's, frees register mm3) */
	"pfmul %%mm3, %%mm4\n"
	"pfmul %%mm3, %%mm5\n"

	/* Combine columns 2 and 3 (frees registers mm4 and mm5) */
	"pfadd %%mm4, %%mm0\n"
	"pfsub %%mm5, %%mm1\n"

	/*
		Current register status:
		========================
		mm0 = C2+C3
		mm1 = C2+C3
		mm2 = Rx | Ry
		mm3 = free
		mm4 = free
		mm5 = free
		mm6 = Rz | Rw
		mm7 = Lw | Lz

		OK! This is current:
		-------------------
		Column 4 = mm2 | mm6


		This is needed:
		--------------
		Column 1 = mm3{swap(mm2)} | mm4{swap(mm6)}
		mm5 = Lz | Lz
		mm7 = Lw | Lw

		After that is done, just multiply, invert signs as needed, then add and store

	*/

	/* set up registeres */
	"pswapd %%mm2, %%mm3\n"		/* mm3 = Ry | Rx */
	"movq %%mm7, %%mm5\n"		/* mm5 = Lw | Lz */
	"pswapd %%mm6, %%mm4\n"		/* mm4 = Rw | Rz */
	"punpckhdq %%mm7, %%mm7\n"	/* mm7 = Lw | Lw */
	"punpckldq %%mm5, %%mm5\n"	/* mm5 = Lz | Lz */


	/* All registers set up! Lets bust out some fp-ops */

	/* multiply w by first column */
	"pfmul %%mm7, %%mm3\n"
	"pfmul %%mm7, %%mm4\n"

	/* multiply z by fourth column (remember 4th = { mm2 | mm6 } )*/
	"pfmul %%mm5, %%mm2\n"
	"pfmul %%mm5, %%mm6\n"

	/* add results of C1 to C2+C3 */
	"pfadd %%mm3, %%mm0\n"
	"pfadd %%mm4, %%mm1\n"

	/* Change signs of fourth column to -, +, +, - */
	"pxor _SIMDx86_float_NEGPOS, %%mm2\n"
	"pxor _SIMDx86_float_POSNEG, %%mm6\n"

	/* add fourth column. Result! */
	"pfadd %%mm2, %%mm0\n"
	"pfadd %%mm6, %%mm1\n"

	/* Store */
	"movq %%mm0,  (%2)\n"
	"movq %%mm1, 8(%2)\n"

	:
	: "r" (pRight), "r" (pLeft), "r" (pOut)
	);


	/* 36 3DNow! Instructions (including movements) vs 70 x87 Instructions (-O3 -fomit-frame-pointer) */

	/* Execute 'femms' if desired */
	#ifndef NO_EMMS
	asm("femms\n");
	#endif



}


//========CPlane==========================================


float  plane_DotNormal(const CPlane* pPlane, const CVec3D* pVec)
{
	float dummy;

	asm(
	"movq  (%1), %%mm0\n"
	"movd 8(%1), %%mm1\n"
	"pfmul  (%2), %%mm0\n"
	"pfmul 8(%2), %%mm1\n"
	"pfadd %%mm0, %%mm1\n"
	"pfacc %%mm1, %%mm1\n"

	"movd %%mm1, -4(%%esp)\n"
	"femms\n"
	"flds -4(%%esp)\n"

	: "=t" (dummy)
	: "r" (pPlane), "r" (pVec)
	);

	return dummy;
	
}

float  plane_DotPlane(const CPlane* pPlane1, const CPlane* pPlane2)
{
	float dummy;
	
	asm(
	"movq  (%1), %%mm0\n"       /* mm0 = P1.y | P1.x */
	"movd 8(%1), %%mm1\n"       /* mm1 = 0    | P1.z */
	"pfmul  (%2), %%mm0\n"      /* mm0 = P1.y*P2.y | P1.x*P2.x */
	"pfmul 8(%2), %%mm1\n"      /* mm1 = 0         | P1.z*P2.z */
	"pfadd %%mm0, %%mm1\n"      /* mm1 = P1.y*P2.y | P1.x*P2.x+P1.z*P2.z */
	"pfacc %%mm1, %%mm1\n"      /* mm1 = ?         | P1.x*P2.x+P1.y*P2.y+P1.z*P2.z */
	
	"movd %%mm1, -4(%%esp)\n"
	"femms\n"
	"flds -4(%%esp)\n"
	
	: "=t" (dummy)
	: "r" (pPlane1), "r" (pPlane2)
	);
	
	return dummy;
}

void  plane_Normalize(CPlane* pOut)
{

	/*
		K6-2 and Athlon can interleave pfmul/pfadd every cycle "4 flops/cycle" (not considering latency...), so
		it is best to make sure both get issued in parallel rather than serialize them.
	*/

	asm(
	"movq  (%0), %%mm0\n"   /* mm0 = b | a */
	"movq 8(%0), %%mm1\n"   /* mm1 = d | c */
	"movq %%mm0, %%mm3\n"	/* mm3 = b | a, saved for later */
	"pfmul %%mm0, %%mm0\n"	/* mm0 = b*b | a*a */
	"movq %%mm1, %%mm4\n"	/* mm4 = d | c, saved for later */
	"pfacc %%mm0, %%mm0\n"	/* mm0 = a*a+b*b | a*a+b*b */
	"pfmul %%mm1, %%mm1\n"	/* mm1 = d*d | c*c (get ready for some dependency chain stalls...) */
	"pand _SIMDx86_float_3DNOW_NO_W_MASK, %%mm1\n" /* mm1 = 0 | c */
	"pfadd %%mm1, %%mm0\n"      /* mm0 = ? | a*a+b*b+c*c */
	"pfrsqrt %%mm0, %%mm1\n"	/* mm1 = 1.0f/sqrtf(a*a+b*b+c*c) | 1/sqrtf(a*a+b*b+c*c) */
	"pfmul %%mm1, %%mm3\n"      /* mm3 = b*invmag | a*invmag */
	"pfmul %%mm1, %%mm4\n"      /* mm4 = d*invmag | c*invmag */
	"movq %%mm3,  (%0)\n"		/* write b and a */
	"movq %%mm4, 8(%0)\n"       /* write d and c */
	
	#if !defined(NO_EMMS)
	"femms"
	#endif
	:
	: "r" (pOut)
	);
	

}

void  plane_NormalizeOf(CPlane* pOut, CPlane* pIn)
{

	/*
		K6-2 and Athlon can interleave pfmul/pfadd every cycle "4 flops/cycle" (not considering latency...), so
		it is best to make sure both get issued in parallel rather than serialize them.
	*/

	asm(
	"movq  (%0), %%mm0\n"   /* mm0 = b | a */
	"movq 8(%0), %%mm1\n"   /* mm1 = d | c */
	"movq %%mm0, %%mm3\n"	/* mm3 = b | a, saved for later */
	"pfmul %%mm0, %%mm0\n"	/* mm0 = b*b | a*a */
	"movq %%mm1, %%mm4\n"	/* mm4 = d | c, saved for later */
	"pfacc %%mm0, %%mm0\n"	/* mm0 = a*a+b*b | a*a+b*b */
	"pfmul %%mm1, %%mm1\n"	/* mm1 = d*d | c*c (get ready for some dependency chain stalls...) */
	"pand _SIMDx86_float_3DNOW_NO_W_MASK, %%mm1\n" /* mm1 = 0 | c */
	"pfadd %%mm1, %%mm0\n"      /* mm0 = ? | a*a+b*b+c*c */
	"pfrsqrt %%mm0, %%mm1\n"	/* mm1 = 1.0f/sqrtf(a*a+b*b+c*c) | 1/sqrtf(a*a+b*b+c*c) */
	"pfmul %%mm1, %%mm3\n"      /* mm3 = b*invmag | a*invmag */
	"pfmul %%mm1, %%mm4\n"      /* mm4 = d*invmag | c*invmag */
	"movq %%mm3,  (%1)\n"		/* write b and a */
	"movq %%mm4, 8(%1)\n"       /* write d and c */
	
	/* Execute FEMMS if desired */
	#if !defined(NO_EMMS)
	"femms"
	#endif
	:
	: "r" (pOut), "r" (pIn)
	);


}

#endif //Endif 0
#endif /*__GNUC__ */
#endif // not _GCC__
} //end MATH
} //end SMF
