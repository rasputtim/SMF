#include "SMF_Config.h"

#ifndef __SMF_MATH_DEFS__
#define __SMF_MATH_DEFS__
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
/// @param v1 A vector of type CVec3D, or a C array of three elements.
/// @param x The x component of a second vector.
/// @param y The y component of a second vector.
/// @param z The z component of a second vector.
/// @see DOT3(), ABSDOT3(), DOT3STRIDED().
#define DOT3_xyz(v1, x, y, z) ((v1)[0] * (x) + (v1)[1] * (y) + (v1)[2] * (z))

/*
================================================================================================

	floating point sign bit tests

================================================================================================
*/

#define IEEE_FLT_SIGNBITSET( a )	(reinterpret_cast<const unsigned int &>(a) >> IEEE_FLT_SIGN_BIT)
#define IEEE_FLT_SIGNBITNOTSET( a )	((~reinterpret_cast<const unsigned int &>(a)) >> IEEE_FLT_SIGN_BIT)
#define IEEE_FLT_ISNOTZERO( a )		(reinterpret_cast<const unsigned int &>(a) & ~(1u<<IEEE_FLT_SIGN_BIT))

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
