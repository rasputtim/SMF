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

#ifndef _SMF__TYPES_H_
#define _SMF__TYPES_H_
#include <stdio.h>
#include <stdarg.h>
#include "sys/SMF_HalfFloat.hpp"
#if defined (_MSC_VER)
#include <string>

#else
#include <stdint.h>

#endif

#include <assert.h>
//include this one for MinGW
#include <stdint.h>
///////////////////////////////////////////////////////////
// define alguns tipos
///////////////////////////////////////////////////////////
#ifndef Float64
#define Float64 double
#endif
#ifndef Float32
#define Float32 float
#endif

#define SMF_OK 1
#define SMF_ERROR -1

#if (__STDC_VERSION__ >= 199901L) || (_MSC_VER >= 1600) || (__GNUC__)
typedef int8_t					sf_s8;  ///< a single byte: -128 - 127.
typedef uint8_t					sf_u8; ///< a single byte: 0-255.
typedef int16_t					sf_s16;  ///< 2 bytes: -32768 - 32767.
typedef uint16_t				sf_u16; ///< 2 bytes: 0 - 65535.
typedef int32_t					sf_s32;  ///< 4 bytes signed: max 2,147,483,647 ~ 2000 million or 2e9.
typedef uint32_t				sf_u32; ///< 4 bytes: 0 - 4,294,967,295 ~ 4000 million or 4e9.
typedef int64_t					sf_s64;  ///< 8 bytes signed. 9,223,372,036,854,775,807 ~ 9e18.
typedef uint64_t				sf_u64; ///< 8 bytes: 18,446,744,073,709,551,615 ~1.8e19.
typedef sf_u16                  halfFloat_t;
typedef half_float::half		sf_hf16;

typedef unsigned int			sf_Uint;
typedef unsigned long			sf_Ulong;
typedef __m128					sf_m128;
typedef double					sf_double;

#else // No boost or unknown if we have C99. Have to guess the following are correct.

#include <limits.h>

//#pragma warning "Not using boost and C99 not defined. Guessing the built-ins for fixed-width types!"

typedef unsigned char sf_u8; ///< a single byte: 0-255.
typedef unsigned short sf_u16; ///< 2 bytes: 0 - 65535.
typedef unsigned long long sf_u64; ///< 8 bytes: 18,446,744,073,709,551,615 ~1.8e19.

typedef signed char sf_s8; ///< a single byte: -128 - 127.
typedef signed short sf_s16; ///< 2 bytes: -32768 - 32767.
typedef __m128					sf_m128;

#if ULONG_MAX == 0xffffffff
typedef unsigned long sf_s32; ///< 4 bytes: 0 - 4,294,967,295 ~ 4000 million or 4e9.
typedef long sf_s32; ///< 4 bytes signed: max 2,147,483,647 ~ 2000 million or 2e9.
#elif UINT_MAX == 0xffffffff
typedef unsigned int sf_s32; ///< 4 bytes: 0 - 4,294,967,295 ~ 4000 million or 4e9.
typedef int sf_s32; ///< 4 bytes signed: max 2,147,483,647 ~ 2000 million or 2e9.
#endif

typedef signed long long sf_s64; ///< 8 bytes signed. 9,223,372,036,854,775,807 ~ 9e18.
typedef double					sf_double;
typedef sf_u16                  halfFloat_t;

#endif
/// As per C99, union-reinterpret should now be safe: http://stackoverflow.com/questions/8511676/portable-data-reinterpretation
union FloatIntReinterpret
{
	float f;
	sf_s32 i;
};

union DoubleU64Reinterpret
{
	double d;
	sf_u64 i;
};


#endif
