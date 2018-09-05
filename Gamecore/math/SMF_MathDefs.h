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
#include "../SMF_Config.h"

#ifndef __SMF_MATH_DEFS__
#define __SMF_MATH_DEFS__

#include <math.h>


#ifdef MACOS_X
// for square root estimate instruction
#ifdef PPC_INTRINSICS
#include <ppc_intrinsics.h>
#endif
// for FLT_MIN
#include <float.h>
#endif


#ifdef INFINITY_FLOAT
#undef INFINITY_FLOAT
#endif

#ifdef FLT_EPSILON
#undef FLT_EPSILON
#endif

//converte graus para radianos e vice versa
#define DEG2RAD(a)				( (a) * CMath::M_DEG2RAD )
#define RAD2DEG(a)				( (a) * CMath::M_RAD2DEG )

#define SEC2MS(t)				( CMath::ftoiFast( (t) * CMath::M_SEC2MS ) )
#define MS2SEC(t)				( (t) * CMath::M_MS2SEC )

#define	ANGLE2SHORT(x)			( CMath::ftoiFast( (x) * 65536.0f / 360.0f ) & 65535 )
#define	SHORT2ANGLE(x)			( (x) * ( 360.0f / 65536.0f ) )

#define	ANGLE2BYTE(x)			( CMath::ftoiFast( (x) * 256.0f / 360.0f ) & 255 )
#define	BYTE2ANGLE(x)			( (x) * ( 360.0f / 256.0f ) )

#define C_FLOAT_TO_INT( x )		(int)(x)




//  Matriz Definitions


#define MATX_MAX_TEMP		1024
#define MATX_QUAD( x )		( ( ( ( x ) + 3 ) & ~3 ) * sizeof( float ) )
#define MATX_CLEAREND()		int s = numRows * numColumns; while( s < ( ( s + 3 ) & ~3 ) ) { mat[s++] = 0.0f; }
#define MATX_ALLOCA( n )	( (float *) _alloca16( MATX_QUAD( n ) ) )
#define MATX_ALLOCAFLOAT( n )	( (float *) _allocafloat16( MATX_QUAD( n ) ) )
#define MATX_SIMD

// Vector Definitions

#define VECX_MAX_TEMP		1024
#define VECX_QUAD( x )		( ( ( ( x ) + 3 ) & ~3 ) * sizeof( float ) )
#define VECX_CLEAREND()		int s = size; while( s < ( ( s + 3) & ~3 ) ) { p[s++] = 0.0f; }
#define VECX_ALLOCA( n )	( (float *) _alloca16( VECX_QUAD( n ) ) )
#define VECX_ALLOCAFLOAT( n )	( (float *) _allocafloat16( VECX_QUAD( n ) ) )

/// Computes the dot product of a CVec3D and another vector given by three floats.
/// \param v1 A vector of type CVec3D, or a C array of three elements.
/// \param x The x component of a second vector.
/// \param y The y component of a second vector.
/// \param z The z component of a second vector.
/// \see DOT3(), ABSDOT3(), DOT3STRIDED().
#define DOT3_xyz(v1, x, y, z) ((v1)[0] * (x) + (v1)[1] * (y) + (v1)[2] * (z))
/*===============================================================================?8


/*===============================================================================?8





/*
================================================================================================

	floating point bit layouts according to the IEEE 754-1985 and 754-2008 standard

================================================================================================
*/

#define IEEE_FLT16_MANTISSA_BITS	10
#define IEEE_FLT16_EXPONENT_BITS	5
#define IEEE_FLT16_EXPONENT_BIAS	15
#define IEEE_FLT16_SIGN_BIT			15
#define IEEE_FLT16_SIGN_MASK		( 1U << IEEE_FLT16_SIGN_BIT )

#define IEEE_FLT_MANTISSA_BITS		23
#define IEEE_FLT_EXPONENT_BITS		8
#define IEEE_FLT_EXPONENT_BIAS		127
#define IEEE_FLT_SIGN_BIT			31
#define IEEE_FLT_SIGN_MASK			( 1UL << IEEE_FLT_SIGN_BIT )

#define IEEE_DBL_MANTISSA_BITS		52
#define IEEE_DBL_EXPONENT_BITS		11
#define IEEE_DBL_EXPONENT_BIAS		1023
#define IEEE_DBL_SIGN_BIT			63
#define IEEE_DBL_SIGN_MASK			( 1ULL << IEEE_DBL_SIGN_BIT )

#define IEEE_DBLE_MANTISSA_BITS		63
#define IEEE_DBLE_EXPONENT_BITS		15
#define IEEE_DBLE_EXPONENT_BIAS		0
#define IEEE_DBLE_SIGN_BIT			79


/*
================================================================================================

	floating point sign bit tests

================================================================================================
*/

#define IEEE_FLT_SIGNBITSET( a )	(reinterpret_cast<const unsigned int &>(a) >> IEEE_FLT_SIGN_BIT)
#define IEEE_FLT_SIGNBITNOTSET( a )	((~reinterpret_cast<const unsigned int &>(a)) >> IEEE_FLT_SIGN_BIT)
#define IEEE_FLT_ISNOTZERO( a )		(reinterpret_cast<const unsigned int &>(a) & ~(1u<<IEEE_FLT_SIGN_BIT))

//===IEEE==================

#define DBL_DIG         15                      /* # of decimal digits of precision */
#define DBL_EPSILON     2.2204460492503131e-016 /* smallest such that 1.0+DBL_EPSILON != 1.0 */
#define DBL_MANT_DIG    53                      /* # of bits in mantissa */
#define IEEE_DBL_MAX         1.7976931348623158e+308 /* max value */
#define DBL_MAX_10_EXP  308                     /* max decimal exponent */
#define DBL_MAX_EXP     1024                    /* max binary exponent */
#define IEEE_DBL_MIN         2.2250738585072014e-308 /* min positive value */
#define DBL_MIN_10_EXP  (-307)                  /* min decimal exponent */
#define DBL_MIN_EXP     (-1021)                 /* min binary exponent */
#define _DBL_RADIX      2                       /* exponent radix */
#define _DBL_ROUNDS     1                       /* addition rounding: near */
#define DBL_SQRT_EPSILON   1.4901161193847656e-08

#define FLT_DIG         6                       /* # of decimal digits of precision */
//#define FLT_EPSILON     1.192092896e-07F        /* smallest such that 1.0+FLT_EPSILON != 1.0 */
#define FLT_GUARD       0
#define FLT_MANT_DIG    24                      /* # of bits in mantissa */
#define IEEE_FLT_MAX         3.402823466e+38F        /* max value */
#define FLT_MAX_10_EXP  38                      /* max decimal exponent */
#define FLT_MAX_EXP     128                     /* max binary exponent */
#define IEEE_FLT_MIN         1.175494351e-38F        /* min positive value */
#define FLT_MIN_10_EXP  (-37)                   /* min decimal exponent */
#define FLT_MIN_EXP     (-125)                  /* min binary exponent */
#define FLT_NORMALIZE   0
#define FLT_RADIX       2                       /* exponent radix */
#define FLT_ROUNDS      1                       /* addition rounding: near */
#define FLT_SQRT_EPSILON   3.4526698300124393e-04




#define FLOATSIGNBITSET(f)		((*(const unsigned long *)&(f)) >> 31)
#define FLOATSIGNBITNOTSET(f)	((~(*(const unsigned long *)&(f))) >> 31)
#define FLOATNOTZERO(f)			((*(const unsigned long *)&(f)) & ~(1<<31) )
#define INTSIGNBITSET(i)		(((const unsigned long)(i)) >> 31)
#define INTSIGNBITNOTSET(i)		((~((const unsigned long)(i))) >> 31)

#define	FLOAT_IS_NAN(x)			(((*(const unsigned long *)&x) & 0x7f800000) == 0x7f800000)
#define FLOAT_IS_INF(x)			(((*(const unsigned long *)&x) & 0x7fffffff) == 0x7f800000)
#define FLOAT_IS_IND(x)			((*(const unsigned long *)&x) == 0xffc00000)
#define	FLOAT_IS_DENORMAL(x)	(((*(const unsigned long *)&x) & 0x7f800000) == 0x00000000 && \
								 ((*(const unsigned long *)&x) & 0x007fffff) != 0x00000000 )

#ifndef M_PI_2
#define M_PI_2     1.57079632679489661923132169164      /* pi/2 */
#endif
#ifndef M_LN2
#define M_LN2      0.69314718055994530941723212146      /* ln(2) */
#endif




enum {
  SMF_IEEE_TYPE_NAN = 1,
  SMF_IEEE_TYPE_INF = 2,
  SMF_IEEE_TYPE_NORMAL = 3,
  SMF_IEEE_TYPE_DENORMAL = 4,
  SMF_IEEE_TYPE_ZERO = 5
} ;

typedef struct  {
  int sign ;
  char mantissa[24] ; /* Actual bits are 0..22, element 23 is \0 */
  int exponent ;
  int type ;
} ieee_float_rep_s ;

typedef struct  {
  int sign ;
  char mantissa[53] ; /* Actual bits are 0..51, element 52 is \0 */
  int exponent ;
  int type ;
} ieee_double_rep_s ;


enum {
  SMF_IEEE_SINGLE_PRECISION = 1,
  SMF_IEEE_DOUBLE_PRECISION = 2,
  SMF_IEEE_EXTENDED_PRECISION = 3
} ;

enum {
  SMF_IEEE_ROUND_TO_NEAREST = 1,
  SMF_IEEE_ROUND_DOWN = 2,
  SMF_IEEE_ROUND_UP = 3,
  SMF_IEEE_ROUND_TO_ZERO = 4
} ;

enum {
  SMF_IEEE_MASK_INVALID = 1,
  SMF_IEEE_MASK_DENORMALIZED = 2,
  SMF_IEEE_MASK_DIVISION_BY_ZERO = 4,
  SMF_IEEE_MASK_OVERFLOW = 8,
  SMF_IEEE_MASK_UNDERFLOW = 16,
  SMF_IEEE_MASK_ALL = 31,
  SMF_IEEE_TRAP_INEXACT = 32
} ;



// Functions annotated with MUST_USE_RESULT require that the user stores the return value, or otherwise
// a warning is printed.
#if _MSC_VER >= 1700
// http://msdn.microsoft.com/en-us/library/jj159529.aspx
#define MUST_USE_RESULT _Check_return_
#elif defined(__clang__) || (defined(__GNUC__) && ((__GNUC__*10000+__GNUC_MINOR*100) >= 30400))
// http://gcc.gnu.org/onlinedocs/gcc-3.4.0/gcc/Function-Attributes.html
#define MUST_USE_RESULT __attribute__((warn_unused_result))
#else
#define MUST_USE_RESULT
#endif




#endif
