/*
  SMF -  Super Math Fabric  
  Copyright (C) 2014 Rasputtim <Rasputtim@hotmail.com>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

  */

#ifndef _SMF__UTIL_STRUCTS_H_
#define _SMF__UTIL_STRUCTS_H_

#ifdef _WIN32
#include <io.h>
#else
#include <sys/io.h>
#endif
#include <string>
#include <vector>
#include <list>
#include <iostream>
#include "../SMF_Config.h"

#ifndef BYTE
typedef unsigned char BYTE;
#endif


using namespace std;
namespace SMF {



typedef struct {
	int		num;
	int		minSize;
	int		maxSize;
	int		totalSize;
} memoryStats_t;
#if 0
typedef struct
CMathState{
	//mem heap
	CHeap *			mem_heap;
	memoryStats_t	mem_total_allocs;
	memoryStats_t	mem_frame_allocs;
	memoryStats_t	mem_frame_frees;
	//math init

}CMathState_s;
#endif	
typedef enum {
	CPUID_NONE							= 0x00000,
	CPUID_UNSUPPORTED					= 0x00001,	// unsupported (386/486)
	CPUID_GENERIC						= 0x00002,	// unrecognized processor
	CPUID_INTEL							= 0x00004,	// Intel
	CPUID_AMD							= 0x00008,	// AMD
	CPUID_MMX							= 0x00010,	// Multi Media Extensions
	CPUID_3DNOW							= 0x00020,	// 3DNow!
	CPUID_SSE							= 0x00040,	// Streaming SIMD Extensions
	CPUID_SSE2							= 0x00080,	// Streaming SIMD Extensions 2
	CPUID_SSE3							= 0x00100,	// Streaming SIMD Extentions 3 aka Prescott's New Instructions
	CPUID_ALTIVEC						= 0x00200,	// AltiVec
	CPUID_HTT							= 0x01000,	// Hyper-Threading Technology
	CPUID_CMOV							= 0x02000,	// Conditional Move (CMOV) and fast floating point comparison (FCOMI) instructions
	CPUID_FTZ							= 0x04000,	// Flush-To-zero mode (denormal results are flushed to zero)
	CPUID_DAZ							= 0x08000,	// Denormals-Are-zero mode (denormal source operands are set to zero)
	CPUID_SSE4							= 0x10000,	// Streaming SIMD Extentions 4
	CPUID_SSE41							= 0x20000,	// Streaming SIMD Extentions 4.1
	CPUID_SSE42							= 0x40000,	// Streaming SIMD Extentions 4.2
	CPUID_AVX							= 0x80000,	// Streaming AVX
} cpuid_t;
} //end namespace SMF
#endif
