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

#include "math/SMF_SimdGeneric.h"
#include "math/SMF_SimdMMX.h"

namespace SMF{
namespace MATH{
//===============================================================
//
//	MMX implementation of CSIMDProcessor
//
//===============================================================

#if defined(MACOS_X) && defined(__i386__)
/*
============
CSIMD_MMX::getName
============
*/
const char * CSIMD_MMX::getName() const {
	return "MMX";
}

#elif defined(_WIN32) && defined(_MSC_VER)
#define EMMS_INSTRUCTION		__asm emms

/*
============
CSIMD_MMX::getName
============
*/
const char * CSIMD_MMX::getName() const {
	return "MMX";
}

/*
================
MMX_Memcpy8B
================
*/
void MMX_Memcpy8B( void *dest, const void *src, const int count ) {
	_asm {
        mov		esi, src
        mov		edi, dest
        mov		ecx, count
        shr		ecx, 3			// 8 bytes per iteration

loop1:
        movq	mm1,  0[ESI]	// Read in source data
        movntq	0[EDI], mm1		// Non-temporal stores

        add		esi, 8
        add		edi, 8
        dec		ecx
        jnz		loop1

	}
	EMMS_INSTRUCTION
}

/*
================
MMX_Memcpy64B

  165MB/sec
================
*/
void MMX_Memcpy64B( void *dest, const void *src, const int count ) {
	_asm {
        mov		esi, src
        mov		edi, dest
        mov		ecx, count
        shr		ecx, 6		// 64 bytes per iteration

loop1:
        prefetchnta 64[ESI]	// Prefetch next loop, non-temporal
        prefetchnta 96[ESI]

        movq mm1,  0[ESI]	// Read in source data
        movq mm2,  8[ESI]
        movq mm3, 16[ESI]
        movq mm4, 24[ESI]
        movq mm5, 32[ESI]
        movq mm6, 40[ESI]
        movq mm7, 48[ESI]
        movq mm0, 56[ESI]

        movntq  0[EDI], mm1	// Non-temporal stores
        movntq  8[EDI], mm2
        movntq 16[EDI], mm3
        movntq 24[EDI], mm4
        movntq 32[EDI], mm5
        movntq 40[EDI], mm6
        movntq 48[EDI], mm7
        movntq 56[EDI], mm0

        add		esi, 64
        add		edi, 64
        dec		ecx
        jnz		loop1
	}
	EMMS_INSTRUCTION
}

/*
================
MMX_Memcpy2kB

  240MB/sec
================
*/
void MMX_Memcpy2kB( void *dest, const void *src, const int count ) {
	sf_u8 *tbuf = (sf_u8 *)_alloca16(2048);
	__asm {
		push	ebx
        mov		esi, src
        mov		ebx, count
        shr		ebx, 11		// 2048 bytes at a time
        mov		edi, dest

loop2k:
        push	edi			// copy 2k into temporary buffer
        mov		edi, tbuf
        mov		ecx, 32

loopMemToL1:
        prefetchnta 64[ESI] // Prefetch next loop, non-temporal
        prefetchnta 96[ESI]

        movq mm1,  0[ESI]	// Read in source data
        movq mm2,  8[ESI]
        movq mm3, 16[ESI]
        movq mm4, 24[ESI]
        movq mm5, 32[ESI]
        movq mm6, 40[ESI]
        movq mm7, 48[ESI]
        movq mm0, 56[ESI]

        movq  0[EDI], mm1	// Store into L1
        movq  8[EDI], mm2
        movq 16[EDI], mm3
        movq 24[EDI], mm4
        movq 32[EDI], mm5
        movq 40[EDI], mm6
        movq 48[EDI], mm7
        movq 56[EDI], mm0
        add		esi, 64
        add		edi, 64
        dec		ecx
        jnz		loopMemToL1

        pop		edi			// Now copy from L1 to system memory
        push	esi
        mov		esi, tbuf
        mov		ecx, 32

loopL1ToMem:
        movq mm1, 0[ESI]	// Read in source data from L1
        movq mm2, 8[ESI]
        movq mm3, 16[ESI]
        movq mm4, 24[ESI]
        movq mm5, 32[ESI]
        movq mm6, 40[ESI]
        movq mm7, 48[ESI]
        movq mm0, 56[ESI]

        movntq 0[EDI], mm1	// Non-temporal stores
        movntq 8[EDI], mm2
        movntq 16[EDI], mm3
        movntq 24[EDI], mm4
        movntq 32[EDI], mm5
        movntq 40[EDI], mm6
        movntq 48[EDI], mm7
        movntq 56[EDI], mm0

        add		esi, 64
        add		edi, 64
        dec		ecx
        jnz		loopL1ToMem

        pop		esi			// Do next 2k block
        dec		ebx
        jnz		loop2k
		pop		ebx
	}
	EMMS_INSTRUCTION
}


/*
================
CSIMD_MMX::Memcpy

  optimized memory copy routine that handles all alignment cases and block sizes efficiently
================
*/
void VPCALL CSIMD_MMX::Memcpy( void *dest0, const void *src0, const int count0 ) {
	// if copying more than 16 bytes and we can copy 8 sf_u8 aligned
	if ( count0 > 16 && !( ( (int)dest0 ^ (int)src0 ) & 7 ) ) {
		sf_u8 *dest = (sf_u8 *)dest0;
		sf_u8 *src = (sf_u8 *)src0;

		// copy up to the first 8 sf_u8 aligned boundary
		int count = ((int)dest) & 7;
		memcpy( dest, src, count );
		dest += count;
		src += count;
		count = count0 - count;

		// if there are multiple blocks of 2kB
		if ( count & ~4095 ) {
			MMX_Memcpy2kB( dest, src, count );
			src += (count & ~2047);
			dest += (count & ~2047);
			count &= 2047;
		}

		// if there are blocks of 64 bytes
		if ( count & ~63 ) {
			MMX_Memcpy64B( dest, src, count );
			src += (count & ~63);
			dest += (count & ~63);
			count &= 63;
		}

		// if there are blocks of 8 bytes
		if ( count & ~7 ) {
			MMX_Memcpy8B( dest, src, count );
			src += (count & ~7);
			dest += (count & ~7);
			count &= 7;
		}

		// copy any remaining bytes
		memcpy( dest, src, count );
	} else {
		// use the regular one if we cannot copy 8 sf_u8 aligned
		memcpy( dest0, src0, count0 );
	}

	// the MMX_Memcpy* functions use MOVNTQ, issue a fence operation
	__asm {
		sfence
	}
}

/*
================
CSIMD_MMX::Memset
================
*/
void VPCALL CSIMD_MMX::Memset( void* dest0, const int val, const int count0 ) {
	union {
		sf_u8	bytes[8];
		sf_u16	words[4];
		sf_u32	dwords[2];
	} dat;

	sf_u8 *dest = (sf_u8 *)dest0;
	int count = count0;

	while ( count > 0 && (((int)dest) & 7) ) {
		*dest = val;
		dest++;
		count--;
	}
	if ( !count ) {
		return;
	}

	dat.bytes[0] = val;
	dat.bytes[1] = val;
	dat.words[1] = dat.words[0];
	dat.dwords[1] = dat.dwords[0];

	if ( count >= 64 ) {
		__asm {
			mov edi, dest
			mov ecx, count
			shr ecx, 6				// 64 bytes per iteration
			movq mm1, dat			// Read in source data
			movq mm2, mm1
			movq mm3, mm1
			movq mm4, mm1
			movq mm5, mm1
			movq mm6, mm1
			movq mm7, mm1
			movq mm0, mm1
loop1:
			movntq  0[EDI], mm1		// Non-temporal stores
			movntq  8[EDI], mm2
			movntq 16[EDI], mm3
			movntq 24[EDI], mm4
			movntq 32[EDI], mm5
			movntq 40[EDI], mm6
			movntq 48[EDI], mm7
			movntq 56[EDI], mm0

			add edi, 64
			dec ecx
			jnz loop1
		}
		dest += ( count & ~63 );
		count &= 63;
	}

	if ( count >= 8 ) {
		__asm {
			mov edi, dest
			mov ecx, count
			shr ecx, 3				// 8 bytes per iteration
			movq mm1, dat			// Read in source data
loop2:
			movntq  0[EDI], mm1		// Non-temporal stores

			add edi, 8
			dec ecx
			jnz loop2
		}
		dest += (count & ~7);
		count &= 7;
	}

	while ( count > 0 ) {
		*dest = val;
		dest++;
		count--;
	}

	EMMS_INSTRUCTION

	// the MMX_Memcpy* functions use MOVNTQ, issue a fence operation
	__asm {
		sfence
	}
}

#else if defined(__GNUC__) && defined (MASM_INTEL)   //gcc asm intel version

#define EMMS_INSTRUCTION		asm("emms");

#define FEMMS_INSTRUCTION  asm volatile ("femms");

/*
============
CSIMD_MMX::getName
============
*/
const char * CSIMD_MMX::getName() const {
	return "MMX";
}

/*
================
MMX_Memcpy8B
================
*/
void MMX_Memcpy8B( void *dest, const void *src, const int count ) {

	asm(
        "mov		esi, %0\n"
        "mov		edi, %1\n"
        "mov		ecx, %2\n"
        "shr		ecx, 3\n"			// 8 bytes per iteration

"loop1%=:\n"
        "movq	mm1,  0[ESI]\n"	// Read in source data
        "movntq	0[EDI], mm1\n"		// Non-temporal stores

        "add		esi, 8\n"
        "add		edi, 8\n"
        "dec		ecx\n"
        "jnz		loop1%=\n"

	::"r"(src),"r"(dest),"m"(count)
    :"esi","edi","ecx");
	EMMS_INSTRUCTION
}

/*
================
MMX_Memcpy64B

  165MB/sec
================
*/
void MMX_Memcpy64B( void *dest, const void *src, const int count ) {

	asm(
        "mov		esi, %1\n"
        "mov		edi, %0\n"
        "mov		ecx, %2\n"
        "shr		ecx, 6\n"		// 64 bytes per iteration

"loop1%=:\n"
        "prefetchnta 64[ESI]\n"	// Prefetch next loop, non-temporal
        "prefetchnta 96[ESI]\n"

        "movq mm1,  0[ESI]\n"	// Read in source data
        "movq mm2,  8[ESI]\n"
        "movq mm3, 16[ESI]\n"
        "movq mm4, 24[ESI]\n"
        "movq mm5, 32[ESI]\n"
        "movq mm6, 40[ESI]\n"
        "movq mm7, 48[ESI]\n"
        "movq mm0, 56[ESI]\n"

        "movntq  0[EDI], mm1\n"	// Non-temporal stores
        "movntq  8[EDI], mm2\n"
        "movntq 16[EDI], mm3\n"
        "movntq 24[EDI], mm4\n"
        "movntq 32[EDI], mm5\n"
        "movntq 40[EDI], mm6\n"
        "movntq 48[EDI], mm7\n"
        "movntq 56[EDI], mm0\n"

        "add		esi, 64\n"
        "add		edi, 64\n"
        "dec		ecx\n"
        "jnz		loop1%=\n"
	::"r"(dest), "r"(src), "m"(count)
    :"esi","edi","ecx");
	EMMS_INSTRUCTION
}

/*
================
MMX_Memcpy2kB

  240MB/sec
================
*/
void MMX_Memcpy2kB( void *dest, const void *src, const int count ) {
	sf_u8 *tbuf = (sf_u8 *)_alloca16(2048);

	asm(
		"push	ebx\n"
        "mov	esi, %1\n"      //src
        "mov	ebx, %2\n"      //count
        "shr	ebx, 11\n"		//2048 bytes at a time
        "mov	edi, %0\n"      //dst

"loop2k%=:\n"
        "push	edi\n"			// copy 2k into temporary buffer
        "mov	edi, %3\n"      //tbuf
        "mov	ecx, 32\n"

"loopMemToL1%=:\n"
        "prefetchnta 64[ESI]\n" // Prefetch next loop, non-temporal
        "prefetchnta 96[ESI]\n"

        "movq mm1,  0[ESI]\n"	// Read in source data
        "movq mm2,  8[ESI]\n"
        "movq mm3, 16[ESI]\n"
        "movq mm4, 24[ESI]\n"
        "movq mm5, 32[ESI]\n"
        "movq mm6, 40[ESI]\n"
        "movq mm7, 48[ESI]\n"
        "movq mm0, 56[ESI]\n"

        "movq  0[EDI], mm1\n"	// Store into L1
        "movq  8[EDI], mm2\n"
        "movq 16[EDI], mm3\n"
        "movq 24[EDI], mm4\n"
        "movq 32[EDI], mm5\n"
        "movq 40[EDI], mm6\n"
        "movq 48[EDI], mm7\n"
        "movq 56[EDI], mm0\n"
        "add	esi, 64\n"
        "add	edi, 64\n"
        "dec	ecx\n"
        "jnz	loopMemToL1%=\n"

        "pop	edi\n"			// Now copy from L1 to system memory
        "push	esi\n"
        "mov	esi, %3\n"
        "mov	ecx, 32\n"

"loopL1ToMem%=:\n"
        "movq mm1, 0[ESI]\n"	// Read in source data from L1
        "movq mm2, 8[ESI]\n"
        "movq mm3, 16[ESI]\n"
        "movq mm4, 24[ESI]\n"
        "movq mm5, 32[ESI]\n"
        "movq mm6, 40[ESI]\n"
        "movq mm7, 48[ESI]\n"
        "movq mm0, 56[ESI]\n"

        "movntq 0[EDI], mm1\n"	// Non-temporal stores
        "movntq 8[EDI], mm2\n"
        "movntq 16[EDI], mm3\n"
        "movntq 24[EDI], mm4\n"
        "movntq 32[EDI], mm5\n"
        "movntq 40[EDI], mm6\n"
        "movntq 48[EDI], mm7\n"
        "movntq 56[EDI], mm0\n"

        "add	esi, 64\n"
        "add	edi, 64\n"
        "dec	ecx\n"
        "jnz	loopL1ToMem%=\n"

        "pop	esi\n"			// Do next 2k block
        "dec	ebx\n"
        "jnz	loop2k%=\n"
		"pop	ebx\n"
	::"m"(dest), "m"(src), "m"(count),"m"(tbuf)
    :"esi","edi","ecx");

	EMMS_INSTRUCTION
}


/*
================
CSIMD_MMX::Memcpy

  optimized memory copy routine that handles all alignment cases and block sizes efficiently
================
*/
void VPCALL CSIMD_MMX::memCopy( void *dest0, const void *src0, const int count0 ) {
	// if copying more than 16 bytes and we can copy 8 sf_u8 aligned
	if ( count0 > 16 && !( ( (int)dest0 ^ (int)src0 ) & 7 ) ) {
		sf_u8 *dest = (sf_u8 *)dest0;
		sf_u8 *src = (sf_u8 *)src0;

		// copy up to the first 8 sf_u8 aligned boundary
		int count = ((int)dest) & 7;
		memcpy( dest, src, count );
		dest += count;
		src += count;
		count = count0 - count;

		// if there are multiple blocks of 2kB
		if ( count & ~4095 ) {
			MMX_Memcpy2kB( dest, src, count );
			src += (count & ~2047);
			dest += (count & ~2047);
			count &= 2047;
		}

		// if there are blocks of 64 bytes
		if ( count & ~63 ) {
			MMX_Memcpy64B( dest, src, count );
			src += (count & ~63);
			dest += (count & ~63);
			count &= 63;
		}

		// if there are blocks of 8 bytes
		if ( count & ~7 ) {
			MMX_Memcpy8B( dest, src, count );
			src += (count & ~7);
			dest += (count & ~7);
			count &= 7;
		}

		// copy any remaining bytes
		memcpy( dest, src, count );
	} else {
		// use the regular one if we cannot copy 8 sf_u8 aligned
		memcpy( dest0, src0, count0 );
	}

	// the MMX_Memcpy* functions use MOVNTQ, issue a fence operation
	asm volatile("sfence");

}

/*
================
CSIMD_MMX::Memset
================
*/
void VPCALL CSIMD_MMX::memSet( void* dest0, const int val, const int count0 ) {
	union {
		sf_u8	bytes[8];
		sf_u16	words[4];
		sf_u32	dwords[2];
	} dat;

	sf_u8 *dest = (sf_u8 *)dest0;
	int count = count0;

	while ( count > 0 && (((int)dest) & 7) ) {
		*dest = val;
		dest++;
		count--;
	}
	if ( !count ) {
		return;
	}

	dat.bytes[0] = val;
	dat.bytes[1] = val;
	dat.words[1] = dat.words[0];
	dat.dwords[1] = dat.dwords[0];

	if ( count >= 64 ) {

		asm(
			"mov edi, %0\n"
			"mov ecx, %2\n"
			"shr ecx, 6\n"				// 64 bytes per iteration
			"movq mm1, %1\n"			// Read in source data
			"movq mm2, mm1\n"
			"movq mm3, mm1\n"
			"movq mm4, mm1\n"
			"movq mm5, mm1\n"
			"movq mm6, mm1\n"
			"movq mm7, mm1\n"
			"movq mm0, mm1\n"
"loop1%=:\n"
			"movntq  0[EDI], mm1\n"		// Non-temporal stores
			"movntq  8[EDI], mm2\n"
			"movntq 16[EDI], mm3\n"
			"movntq 24[EDI], mm4\n"
			"movntq 32[EDI], mm5\n"
			"movntq 40[EDI], mm6\n"
			"movntq 48[EDI], mm7\n"
			"movntq 56[EDI], mm0\n"

			"add edi, 64\n"
			"dec ecx\n"
			"jnz loop1%=\n"
		::"r"(dest), "m"(dat), "m"(count)
        :"esi","edi","ecx");

		dest += ( count & ~63 );
		count &= 63;
	}

	if ( count >= 8 ) {

		asm(
			"mov edi, %0\n"
			"mov ecx, %2\n"
			"shr ecx, 3\n"				// 8 bytes per iteration
			"movq mm1, %1\n"			// Read in source data
"loop2%=:\n"
			"movntq  0[EDI], mm1\n"		// Non-temporal stores

			"add edi, 8\n"
			"dec ecx\n"
			"jnz loop2%=\n"
		::"r"(dest), "m"(dat), "m"(count)
        :"esi","edi","ecx");

		dest += (count & ~7);
		count &= 7;
	}

	while ( count > 0 ) {
		*dest = val;
		dest++;
		count--;
	}

	EMMS_INSTRUCTION

	// the MMX_Memcpy* functions use MOVNTQ, issue a fence operation
	asm volatile("sfence" : : : "memory");
}


#endif // defined

} //end MATH
} //end SMF
