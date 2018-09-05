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

#ifndef _SMF__MATH_MATH_H___
#define _SMF__MATH_MATH_H___

#include "../SMF_Config.h"
#include "SMF_MathDefs.h"
#include "../exceptions/all.h"



template<class T> SMF_INLINE_FORCED int	MaxIndex( T x, T y ) { return  ( x > y ) ? 0 : 1; }
template<class T> SMF_INLINE_FORCED int	MinIndex( T x, T y ) { return ( x < y ) ? 0 : 1; }

template<class T> SMF_INLINE_FORCED T	Max3( T x, T y, T z ) { return ( x > y ) ? ( ( x > z ) ? x : z ) : ( ( y > z ) ? y : z ); }
template<class T> SMF_INLINE_FORCED T	Min3( T x, T y, T z ) { return ( x < y ) ? ( ( x < z ) ? x : z ) : ( ( y < z ) ? y : z ); }
template<class T> SMF_INLINE_FORCED int	Max3Index( T x, T y, T z ) { return ( x > y ) ? ( ( x > z ) ? 0 : 2 ) : ( ( y > z ) ? 1 : 2 ); }
template<class T> SMF_INLINE_FORCED int	Min3Index( T x, T y, T z ) { return ( x < y ) ? ( ( x < z ) ? 0 : 2 ) : ( ( y < z ) ? 1 : 2 ); }

template<class T> SMF_INLINE_FORCED T	sign( T f ) { return ( f > 0 ) ? 1 : ( ( f < 0 ) ? -1 : 0 ); }
template<class T> SMF_INLINE_FORCED T	square( T x ) { return x * x; }
template<class T> SMF_INLINE_FORCED T	Cube( T x ) { return x * x * x; }

#define MIN(x,y)     (((x) < (y)) ? (x) : (y))
#define MAX(x,y)     (((x) > (y)) ? (x) : (y))
#define MID(x,y,z)   MAX((x), MIN((y), (z)))
/// Computes the smallest of three values.
/** \see clamp(), clamp01(), MAX(). */
#define CLAMP(x,y,z) MIN(MIN(x,y),z)
/** \return The absolute value of a. */
#define  ABS(a) ((a) >= 0 ? (a) : -(a))

/// Clamps the given input value to the range [min, max].
/** \see clamp01(), min(), max(). */
template<typename T>
inline T clamp(const T &val, const T &floor, const T &ceil)
{
	//SMF_ASSERT(floor <= ceil);
	return val <= ceil ? (val >= floor ? val : floor) : ceil;
}

/// Clamps the given input value to the range [0, 1].
/** \see clamp(), min(), max(). */
template<typename T>
inline T clamp01(const T &val) { return clamp(val, T(0), T(1)); }

#ifdef __GNUC__
#define NOT_NECESSARILY_USED __attribute__ ((unused))
#else
#define NOT_NECESSARILY_USED
#endif



namespace SMF {
namespace MATH{


/*
================================================
halfFloat_t
================================================
*/


// GPU half-float bit patterns
#define HF_MANTISSA(x)	(x&1023)
#define HF_EXP(x)		((x&32767)>>10)
#define HF_SIGN(x)		((x&32768)?-1:1)

/*
========================
F16toF32
========================
*/
SMF_INLINE float F16toF32( halfFloat_t x ) {
	int e = HF_EXP( x );
	int m = HF_MANTISSA( x );
	int s = HF_SIGN( x );

	if ( 0 < e && e < 31 ) {
		return s * powf( 2.0f, ( e - 15.0f ) ) * ( 1 + m / 1024.0f );
	} else if ( m == 0 ) {
        return s * 0.0f;
	}
    return s * powf( 2.0f, -14.0f ) * ( m / 1024.0f );
}

/*
========================
F32toF16
========================
*/
SMF_INLINE halfFloat_t F32toF16( float a ) {
	unsigned int f = *(unsigned *)( &a );
	unsigned int signbit  = ( f & 0x80000000 ) >> 16;
	int exponent = ( ( f & 0x7F800000 ) >> 23 ) - 112;
	unsigned int mantissa = ( f & 0x007FFFFF );

	if ( exponent <= 0 ) {
		return 0;
	}
	if ( exponent > 30 ) {
		return (halfFloat_t)( signbit | 0x7BFF );
	}

	return (halfFloat_t)( signbit | ( exponent << 10 ) | ( mantissa >> 13 ) );
}


/*
================================================================================================

	floating point special value tests

================================================================================================
*/

/// Returns the bit pattern of the given float as a unsigned int 32.
sf_s32 ReinterpretAsU32(float f);
/// Returns the bit pattern of the given float as a unsigned int 64.
sf_u64 ReinterpretAsU64(double d);

/// Converts the bit pattern specified by the given integer to a floating point (this is a binary conversion, not numeral!).
float ReinterpretAsFloat(sf_s32 i);
double ReinterpretAsDouble(sf_u64 i);

/// Returns true if the given value is not an inf or a nan.
template<typename T> SMF_INLINE_FORCED bool isFinite(T /*value*/) { return true; }

template<> SMF_INLINE_FORCED bool isFinite<float>(float f) { return (ReinterpretAsU32(f) << 1) < 0xFF000000u; }
template<> SMF_INLINE_FORCED bool isFinite<double>(double d) { return (ReinterpretAsU64(d) << 1) < 0xFFE0000000000000ULL; }

/// Returns true if the given value is a not-a-number.
SMF_INLINE_FORCED bool isNan(float f) { return (ReinterpretAsU32(f) << 1) > 0xFF000000u; }
SMF_INLINE_FORCED bool isNan(double d) { return (ReinterpretAsU64(d) << 1) > 0xFFE0000000000000ULL; }

/// Returns true if the given value is +inf or -inf.
SMF_INLINE_FORCED bool isInf(float f) { return (ReinterpretAsU32(f) << 1) == 0xFF000000u; }
SMF_INLINE_FORCED bool isInf(double d) { return (ReinterpretAsU64(d) << 1) == 0xFFE0000000000000ULL; }

/*
\brief
*/
SMF_INLINE_EXTERN bool IEEE_FLT_IS_NAN( float x ) {
	return isNan(x);
}

/*
\brief
*/
SMF_INLINE_EXTERN bool IEEE_FLT_IS_INF( float x ) {
	return isInf(x);
}

/*
\brief
*/
SMF_INLINE_EXTERN bool IEEE_FLT_IS_INF_NAN( float x ) {
	return isInf(x) || isNan(x);
}

/*
\brief
*/
SMF_INLINE_EXTERN bool IEEE_FLT_IS_IND( float x ) {
	return	(reinterpret_cast<const unsigned int &>(x) == 0xffc00000);
}

/*
\brief
*/
SMF_INLINE_EXTERN bool IEEE_FLT_IS_DENORMAL( float x ) {
	return ((reinterpret_cast<const unsigned int &>(x) & 0x7f800000) == 0x00000000 &&
			(reinterpret_cast<const unsigned int &>(x) & 0x007fffff) != 0x00000000 );
}


/*
\brief
*/
template<class type>
SMF_INLINE_EXTERN bool isNAN( const type &v ) {
	for ( int i = 0; i < v.getDimension(); i++ ) {
		const float f = v.toFloatPtr()[i];
		if ( IEEE_FLT_IS_NAN( f ) || IEEE_FLT_IS_INF( f ) || IEEE_FLT_IS_IND( f ) ) {
			return true;
		}
	}
	return false;
}

/*
\brief
*/
template<class type>
SMF_INLINE_EXTERN bool isValid( const type &v ) {
	for ( int i = 0; i < v.getDimension(); i++ ) {
		const float f = v.toFloatPtr()[i];
		if ( IEEE_FLT_IS_NAN( f ) || IEEE_FLT_IS_INF( f ) || IEEE_FLT_IS_IND( f ) || IEEE_FLT_IS_DENORMAL( f ) ) {
			return false;
		}
	}
	return true;
}

/*
\brief
*/
template<>
SMF_INLINE_EXTERN bool isValid( const float & f ) {	// these parameter must be a reference for the function to be considered a specialization
bool t1 = IEEE_FLT_IS_NAN( f );
bool t2 = IEEE_FLT_IS_INF( f );
bool t3 = IEEE_FLT_IS_IND( f );
bool t4 = IEEE_FLT_IS_DENORMAL( f );
bool ret = !(t1 || t2 || t3 || t4);
return ret;
//return !( IEEE_FLT_IS_NAN( f ) || IEEE_FLT_IS_INF( f ) || IEEE_FLT_IS_IND( f ) || IEEE_FLT_IS_DENORMAL( f ) );
}

/*
\brief
*/
template<>
SMF_INLINE_EXTERN bool isNAN( const float & f ) {	// these parameter must be a reference for the function to be considered a specialization
	if ( IEEE_FLT_IS_NAN( f ) || IEEE_FLT_IS_INF( f ) || IEEE_FLT_IS_IND( f ) ) {
		return true;
	}
	return false;
}

/*
\brief
\return Returns true if any scalar is greater than the range or less than the negative range.
*/
template<class type>
SMF_INLINE_EXTERN bool isInRange( const type &v, const float range ) {
	for ( int i = 0; i < v.getDimension(); i++ ) {
		const float f = v.toFloatPtr()[i];
		if ( f > range || f < -range ) {
			return false;
		}
	}
	return true;
}




/// Swaps the two values.
template<typename T>
void Swap(T &a, T &b)
{
	T temp = a;
	a = b;
	b = temp;
}
/* For the IEEE float format the bits are found from the following
   masks,

   sign      = 0x80000000
   exponent  = 0x7f800000
   mantisssa = 0x007fffff

   For the IEEE double format the masks are,

   sign      = 0x8000000000000000
   exponent  = 0x7ff0000000000000
   mantissa  = 0x000fffffffffffff

   */
void ieee_float_to_rep (const float * x, ieee_float_rep_s * r) ;
void ieee_double_to_rep (const double * x, ieee_double_rep_s * r) ;

void ieee_printf_float (const float * x) ;
void ieee_printf_double (const double * x) ;

void ieee_fprintf_float (FILE * stream, const float * x) ;
void ieee_fprintf_double (FILE * stream, const double * x) ;
/**
\brief get the iee modo of work frok the string passed
\param description the string to be availated
\param [out] precison precision configuration
\param [out] rounding rounding configuration
\param [out] exception_mask exception mask configuration
**/
int ieee_read_mode_string (const char * description,
                           int * precision,
                           int * rounding,
                           int * exception_mask);

/**
\brief setup the fpu according the environment variable SMF_IEEE_MODE
**/
void ieee_env_setup (void);
/**
 * \class CMath
 *
 * \ingroup SMF_Math
 *
 * \if pt_br
 * \brief Funções matemáticas diversas
 * \elseif us_en
 * \brief Mathematical functions
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CMath {
public:


	static void					init();



	/**
	\brief verifica se a classe CMath foi inicializada
	\note existem métodos que requerem a inicialização
	\return verdadeiro caso esteja inicializada, falso, caso contrário
	**/
	static bool					isInitialized();
	/// reciprocal square root, returns huge number when x == 0.0
	static float				rSqrt( float x );
	/// inverse square root with 32 bits precision, returns huge number when x == 0.0
	static float				invSqrt( float x );
	/// inverse square root with 16 bits precision, returns huge number when x == 0.0
	static float				invSqrt16( float x );
	/// inverse square root with 64 bits precision, returns huge number when x == 0.0
	static double				invSqrt64( double x );
	/// Cubic root with 32 bits precision
	static float				cbrt( float x );
	/// Cubic root with 16 bits precision
	static float				cbrt16( float x );
	/// Cubic root with 64 bits precision
	static double				cbrt64( double x );


	/// square root with 32 bits precision
	static float				sqrt( float x );
	/// square root with 16 bits precision
	static float				sqrt16( float x );
	/// square root with 64 bits precision
	static double				sqrt64( double x );


	/// sine with 32 bits precision
	static float				sin( float a );
	/// sine with 16 bits precision, maximum absolute error is 2.3082e-09
	static float				sin16( float a );
	/// sine with 64 bits precision
	static double				sin64( double a );
	/// cosine with 32 bits precision
	static float				cos( float a );
	/// cosine with 16 bits precision, maximum absolute error is 2.3082e-09
	static float				cos16( float a );
	/// cosine with 64 bits precision
	static double				cos64( double a );

	/// sine and cosine with 32 bits precision
	static void					sincos( float a, float &s, float &c );
	/// sine and cosine with 16 bits precision
	static void					sincos16( float a, float &s, float &c );
	/// sine and cosine with 64 bits precision
	static void					sinCos64( double a, double &s, double &c );


	/// tangent with 32 bits precision
	static float				tan( float a );
	/// tangent with 16 bits precision, maximum absolute error is 1.8897e-08
	static float				tan16( float a );
	/// tangent with 64 bits precision
	static double				tan64( double a );

	/// arc sine with 32 bits precision, input is clamped to [-1, 1] to avoid a silent NaN
	static float				asin( float a );
	/// arc sine with 16 bits precision, maximum absolute error is 6.7626e-05
	static float				asin16( float a );
	/// arc sine with 64 bits precision
	static double				asin64( double a );

	/// arc cosine with 32 bits precision, input is clamped to [-1, 1] to avoid a silent NaN
	static float				acos( float a );
	/// arc cosine with 16 bits precision, maximum absolute error is 6.7626e-05
	static float				acos16( float a );
	/// arc cosine with 64 bits precision
	static double				acos64( double a );

	/// arc tangent with 32 bits precision
	static float				atan( float a );
	/// arc tangent with 16 bits precision, maximum absolute error is 1.3593e-08
	static float				atan16( float a );
	/// arc tangent with 64 bits precision
	static double				atan64( double a );

	/// arc tangent with 32 bits precision
	static float				atan( float y, float x );
	/// arc tangent with 16 bits precision, maximum absolute error is 1.3593e-08
	static float				atan16( float y, float x );
	/// arc tangent with 64 bits precision
	static double				atan64( double y, double x );
	/**
	\brief hyperbolic sin
	**/
	static float				asinh (const float x);
	/**
	\brief hyperbolic cos
	**/
	static float				acosh (const float x);
	/**
	\brief hyperbolic tan
	**/
	static float				atanh (const float x);
	/// x raised to the power y with 32 bits precision
	static float				pow( float x, float y );
	/// x raised to the power y with 16 bits precision
	static float				pow16( float x, float y );
	/// x raised to the power y with 64 bits precision
	static double				pow64( double x, double y );

	/// x raised to the power 2 with 32 bits precision
	static float				square( float x );
	/// x raised to the power 2 with 16 bits precision
	static float				square16( float x );
	/// x raised to the power 2 with 64 bits precision
	static double				square64( double x );

	/// e raised to the power f with 32 bits precision
	static float				exp( float f );
	/// e raised to the power f with 16 bits precision
	static float				exp16( float f );
	/// e raised to the power f with 64 bits precision
	static double				exp64( double f );
	/**
	\brief computes the value of \log(1+x) in a way that is accurate for small x. It provides an alternative to the BSD math function log1p(x).
	**/
	static float				log1p (const float x);
	/// natural logarithm with 32 bits precision
	static float				log( float f );
	/// natural logarithm with 16 bits precision
	static float				log16( float f );
	/// natural logarithm with 64 bits precision
	static double				log64( double f );

	/// integral x raised to the power y
	static int					ipow( int x, int y );
	/// integral base-2 logarithm of the floating point value
	static int					ilog2( float f );
	/// integral base-2 logarithm of the integer value
	static int					ilog2( int i );

	/// minumum number of bits required to represent ceil( f )
	static int					bitsForFloat( float f );
	/// minumum number of bits required to represent i
	static int					bitsForInteger( int i );
	/// returns 0x00000000 if x >= 0.0f and returns 0xFFFFFFFF if x <= -0.0f
	static int					maskForFloatSign( float f );
	/// returns 0x00000000 if x >= 0 and returns 0xFFFFFFFF if x < 0
	static int					maskForIntegerSign( int i );
	/// round x down to the nearest power of 2
	static int					floorPowerOfTwo( int x );
	/// round x up to the nearest power of 2
	static int					ceilPowerOfTwo( int x );
	/// returns true if x is a power of 2
	static bool					isPowerOfTwo( int x );
	/// returns the number of 1 bits in x
	static int					bitCount( int x );
	/// returns the bit reverse of x
	static int					bitReverse( int x );

	/// returns the absolute value of the integer value (for reference only)
	static int					abs( int x );
	/// returns the absolute value of the floating point value
	static float				fabs( float f );
	/// returns the largest integer that is less than or equal to the given value
	static float				floor( float f );
	/// returns the smallest integer that is greater than or equal to the given value
	static float				ceil( float f );
	/// returns the nearest integer
	static float				rint( float f );

	/// float to int conversion
	static int					ftoi( float f );
	/// fast float to int conversion but uses current FPU round mode (default round nearest)
	static int					ftoiFast( float f );
	/// float to char conversion
	static char					ftoi8( float f );
	/// float to short conversion
	static short				ftoi16( float f );
	/// float to unsigned short conversion
	static unsigned short		ftoui16( float f );
	/// float to byte conversion, the result is clamped to the range [0-255]
	static sf_u8				ftob( float f );
	/// float to long conversion
	static unsigned long		ftol( float f );
	/// fast float to long conversion but uses current FPU round mode (default round nearest)
	static unsigned long		ftolFast( float );


	/// Returns f rounded up to the next integer, as integer.
	/** \see ceil(), floor(), floorInt(), round(), roundInt(). */
	static int ceilInt(float f);
	/// Returns f rounded down to the previous integer, as integer.
	/** \see ceil(), ceilInt(), floor(), round(), roundInt(). */
	static int floorInt(float f);
	/// Returns f rounded to the nearest integer, as float.
	/** \see ceil(), ceilInt(), floor(), floorInt(), roundInt(). */
	static float round(float f);
	/**
	\brief return f rounded to prec precision
	\param prec precision to round
	\param f number to be rounded
	**/
	static float round(float f,float prec);

	/// Returns f rounded to the nearest integer, as integer.
	/** \see ceil(), ceilInt(), floor(), floorInt(), round(). */
	static int roundInt(float f);

	/// Returns -1 or 1 depending on the sign of f.
	/** \see signOrZero(). */
	static float sign(float f);
	/// Returns 0 if f is zero up to the given epsilon. Otherwise returns -1 or 1 depending on the sign of f.
	/** \see sign(). */
	static float signOrZero(float f, float epsilon = 1e-8f);

	/// Returns the fractional part of x.
	static float				frac(float x);

	/** Compares the two values for equality, allowing the given amount of absolute error. */
	static bool equalsAbs(float a, float b, float epsilon =CMath::EPSILON_SuperLow);
	/** Compares the two numbers for equality, allowing  a given ammount of error (epsilon)
	**/
	static bool equals(float a, float b);
	/// Return true if the number a  is near zero
	static bool nearZero(float a);

	/// Returns 1/x, the reciprocal of x, using a fast approximation (SSE rcp instruction).

	SMF_INLINE_FORCED static float recipFast(float x)
{

	return 1.f / x;

}
	static signed char			clampChar( int i );
	static signed short			clampShort( int i );
	static int					clampInt( int min, int max, int value );
	static float				clampFloat( float min, float max, float value );

	static float				angleNormalize360( float angle );
	static float				angleNormalize180( float angle );
	static float				angleDelta( float angle1, float angle2 );

	static int					floatToBits( float f, int exponentBits, int mantissaBits );
	static float				bitsToFloat( int i, int exponentBits, int mantissaBits );

	static int					floatHash( const float *array, const int numFloats );

	/**
	\brief SIMD Operation for math calculation
	false means not to use SSE system for math calculation
	true means use SSE system for math calculations
	\warning this shoulb be always true. use false only for benchmarking
	**/
	static bool					MATH_AUTOMATIC_SSE;
	/// pi
	static const float			PI;
	/// pi * 2
	static const float			TWO_PI;
	/// pi / 2
	static const float			HALF_PI;
	/// pi / 4
	static const float			ONEFOURTH_PI;
	/// e
	static const float			E;
	/// sqrt( 2 )
	static const float			SQRT_TWO;
	/// sqrt( 3 )
	static const float			SQRT_THREE;
	/// sqrt( 1 / 2 )
	static const float			SQRT_1OVER2;
	/// sqrt( 1 / 3 )
	static const float			SQRT_1OVER3;
	/// degrees to radians multiplier
	static const float			M_DEG2RAD;
	/// radians to degrees multiplier
	static const float			M_RAD2DEG;
	/// seconds to milliseconds multiplier
	static const float			M_SEC2MS;
	/// milliseconds to seconds multiplier
	static const float			M_MS2SEC;
	/// huge number which should be larger than any valid number used
	static const float			INFINITY_FLOAT;
	/// negative INFINITE_FLOAT
	static const float			NEG_INFINITY_FLOAT;
	/// Represents a floating-point not-a-number. \note Never compare a float against nan, use isFinite() instead!
	static const float			NAN_FLOAT;
	/// smallest positive number such that 1.0+FLT_EPSILON != 1.0
	static const float		    EPSILON_High;
	static const float          EPSILON_Medium;
	static const float          EPSILON_Low;
	static const float          EPSILON_SuperLow;
	static const float			FLT_EPSILON;
	/// smallest non-denormal 32-bit floating point value
	static const float			FLT_SMALLEST_NON_DENORMAL;
	/// 0.0f
	static const float			ZERO;
	/// The floating point representation for +\f$\inf\f$.
	static const float NOT_NECESSARILY_USED inf;
	/// The floating point representation for -\f$\inf\f$.
	static const float NOT_NECESSARILY_USED negInf;
	/// Represents a floating-point not-a-number. \note Never compare a float against nan, use isFinite() or isValid() instead!
	static const float NOT_NECESSARILY_USED nan;

private:
	enum {
		LOOKUP_BITS				= 8,
		EXP_POS					= 23,
		EXP_BIAS				= 127,
		LOOKUP_POS				= (EXP_POS-LOOKUP_BITS),
		SEED_POS				= (EXP_POS-8),
		SQRT_TABLE_SIZE			= (2<<LOOKUP_BITS),
		LOOKUP_MASK				= (SQRT_TABLE_SIZE-1)
	};

	union _flint {
		sf_u32				i;
		float				f;
	};

	static sf_u32			iSqrt[SQRT_TABLE_SIZE];
	static bool					initialized;
};


template< typename T >
SMF_INLINE T lerp( const T from, const T to, float f ) {
	return from + ( ( to - from ) * f );
}

template<>
SMF_INLINE int lerp( const int from, const int to, float f ) {
	return CMath::ftoi( (float) from + ( ( (float) to - (float) from ) * f ) );
}

SMF_INLINE_FORCED float CMath::rSqrt( float x ) {
	if (x == 0.0f )return 0.0f;
	if (! x > FLT_SMALLEST_NON_DENORMAL ) return INFINITY_FLOAT;
	long i;
	float y, r;

	y = x * 0.5f;
	i = *reinterpret_cast<long *>( &x );
	i = 0x5f3759df - ( i >> 1 );
	r = *reinterpret_cast<float *>( &i );
	r = r * ( 1.5f - r * r * y );
	return r;
}

SMF_INLINE_FORCED float CMath::invSqrt16( float x ) {
	if (x == 0.0f )return 0.0f;
	if (! x > FLT_SMALLEST_NON_DENORMAL ) return INFINITY_FLOAT;
	sf_u32 a = ((union _flint*)(&x))->i;
	union _flint seed;

	if (!initialized) init();

	double y = x * 0.5f;
	seed.i = (( ( (3*EXP_BIAS-1) - ( (a >> EXP_POS) & 0xFF) ) >> 1)<<EXP_POS) | iSqrt[(a >> (EXP_POS-LOOKUP_BITS)) & LOOKUP_MASK];
	double r = seed.f;
	r = r * ( 1.5f - r * r * y );
	return (float) r;
}

SMF_INLINE_FORCED float CMath::invSqrt( float x ) {
	if (x == 0.0f )return 0.0f;
	if (! x > FLT_SMALLEST_NON_DENORMAL ) return INFINITY_FLOAT;

	sf_u32 a = ((union _flint*)(&x))->i;
	union _flint seed;

	if (!initialized) init();

	double y = x * 0.5f;
	seed.i = (( ( (3*EXP_BIAS-1) - ( (a >> EXP_POS) & 0xFF) ) >> 1)<<EXP_POS) | iSqrt[(a >> (EXP_POS-LOOKUP_BITS)) & LOOKUP_MASK];
	double r = seed.f;
	r = r * ( 1.5f - r * r * y );
	r = r * ( 1.5f - r * r * y );
	return (float) r;
}

SMF_INLINE_FORCED double CMath::invSqrt64( double x ) {

	sf_u32 a = ((union _flint*)(&x))->i;
	union _flint seed;

	if (!initialized) init();

	double y = x * 0.5f;
	seed.i = (( ( (3*EXP_BIAS-1) - ( (a >> EXP_POS) & 0xFF) ) >> 1)<<EXP_POS) | iSqrt[(a >> (EXP_POS-LOOKUP_BITS)) & LOOKUP_MASK];
	double r = seed.f;
	r = r * ( 1.5f - r * r * y );
	r = r * ( 1.5f - r * r * y );
	r = r * ( 1.5f - r * r * y );
	return r;
}

SMF_INLINE_FORCED float CMath::cbrt( float x ) {
return ((x) > 0.0 ? std::pow((float)(x), 1.0/3.0) : \
                          ((x) < 0.0 ? -std::pow((float)-(x), 1.0/3.0) : 0.0));
}
//! todo: convert code to 16 byte precision
SMF_INLINE_FORCED float CMath::cbrt16( float x ) {
return ((x) > 0.0 ? std::pow((float)(x), 1.0/3.0) : \
                          ((x) < 0.0 ? -std::pow((float)-(x), 1.0/3.0) : 0.0));
}

SMF_INLINE_FORCED double CMath::cbrt64( double x ) {
return ((x) > 0.0 ? std::pow((double)(x), 1.0/3.0) : \
                          ((x) < 0.0 ? -std::pow((double)-(x), 1.0/3.0) : 0.0));
}

SMF_INLINE_FORCED float CMath::sqrt16( float x ) {
	if (x<0.0f) return NAN_FLOAT;
	return x * invSqrt16( x );
}

SMF_INLINE_FORCED float CMath::sqrt( float x ) {
	if (x<0.0f) return NAN_FLOAT;
	return x * invSqrt( x );
}

SMF_INLINE_FORCED double CMath::sqrt64( double x ) {
	if (x<0.0f) return 0.0f;
	return x * invSqrt64( x );
}

SMF_INLINE_FORCED float CMath::sin( float a ) {

    return sinf(a);

}

SMF_INLINE_FORCED float CMath::sin16( float a ) {
	float s;

	if ( ( a < 0.0f ) || ( a >= TWO_PI ) ) {
		a -= floorf( a / TWO_PI ) * TWO_PI;

	}
#if 1
	if ( a < PI ) {
		if ( a > HALF_PI ) {
			a = PI - a;
		}
	} else {
		if ( a > PI + HALF_PI ) {
			a = a - TWO_PI;
		} else {
			a = PI - a;
		}
	}
#else
	a = PI - a;
	if ( std::fabs( a ) >= HALF_PI ) {
		a = ( ( a < 0.0f ) ? -PI : PI ) - a;
	}
#endif
	s = a * a;
	return a * ( ( ( ( ( -2.39e-08f * s + 2.7526e-06f ) * s - 1.98409e-04f ) * s + 8.3333315e-03f ) * s - 1.666666664e-01f ) * s + 1.0f );
}

SMF_INLINE_FORCED double CMath::sin64( double a ) {
	return std::sin( a );
}

SMF_INLINE_FORCED float CMath::cos( float a ) {

    return cosf( a );

}

SMF_INLINE_FORCED float CMath::cos16( float a ) {
	float s, d;

	if ( ( a < 0.0f ) || ( a >= TWO_PI ) ) {
		a -= floorf( a / TWO_PI ) * TWO_PI;
	}
#if 1
	if ( a < PI ) {
		if ( a > HALF_PI ) {
			a = PI - a;
			d = -1.0f;
		} else {
			d = 1.0f;
		}
	} else {
		if ( a > PI + HALF_PI ) {
			a = a - TWO_PI;
			d = 1.0f;
		} else {
			a = PI - a;
			d = -1.0f;
		}
	}
#else
	a = PI - a;
	if ( std::fabs( a ) >= HALF_PI ) {
		a = ( ( a < 0.0f ) ? -PI : PI ) - a;
		d = 1.0f;
	} else {
		d = -1.0f;
	}
#endif
	s = a * a;
	return d * ( ( ( ( ( -2.605e-07f * s + 2.47609e-05f ) * s - 1.3888397e-03f ) * s + 4.16666418e-02f ) * s - 4.999999963e-01f ) * s + 1.0f );
}

SMF_INLINE_FORCED double CMath::cos64( double a ) {
	return std::cos( a );
}

SMF_INLINE_FORCED void CMath::sincos( float a, float &s, float &c ) {
#if defined(_MSC_VER)
	_asm {
		fld		a
		fsincos
		mov		ecx, c
		mov		edx, s
		fstp	dword ptr [ecx]
		fstp	dword ptr [edx]
	}
#else
	s = sinf( a );
	c = cosf( a );
#endif
}

SMF_INLINE_FORCED void CMath::sincos16( float a, float &s, float &c ) {
	float t, d;

	if ( ( a < 0.0f ) || ( a >= CMath::TWO_PI ) ) {
		a -= floorf( a / CMath::TWO_PI ) * CMath::TWO_PI;
	}
#if 1
	if ( a < PI ) {
		if ( a > HALF_PI ) {
			a = PI - a;
			d = -1.0f;
		} else {
			d = 1.0f;
		}
	} else {
		if ( a > PI + HALF_PI ) {
			a = a - TWO_PI;
			d = 1.0f;
		} else {
			a = PI - a;
			d = -1.0f;
		}
	}
#else
	a = PI - a;
	if ( std::fabs( a ) >= HALF_PI ) {
		a = ( ( a < 0.0f ) ? -PI : PI ) - a;
		d = 1.0f;
	} else {
		d = -1.0f;
	}
#endif
	t = a * a;
	s = a * ( ( ( ( ( -2.39e-08f * t + 2.7526e-06f ) * t - 1.98409e-04f ) * t + 8.3333315e-03f ) * t - 1.666666664e-01f ) * t + 1.0f );
	c = d * ( ( ( ( ( -2.605e-07f * t + 2.47609e-05f ) * t - 1.3888397e-03f ) * t + 4.16666418e-02f ) * t - 4.999999963e-01f ) * t + 1.0f );
}

SMF_INLINE_FORCED void CMath::sinCos64( double a, double &s, double &c ) {
#if defined(_MSC_VER)
	_asm {
		fld		a
		fsincos
		mov		ecx, c
		mov		edx, s
		fstp	qword ptr [ecx]
		fstp	qword ptr [edx]
	}
#else
	s = std::sin( a );
	c = std::cos( a );
#endif
}

SMF_INLINE_FORCED float CMath::tan( float a ) {
	return tanf( a );
}

SMF_INLINE_FORCED float CMath::tan16( float a ) {
	float s;
	bool reciprocal;

	if ( ( a < 0.0f ) || ( a >= PI ) ) {
		a -= floorf( a / PI ) * PI;
	}
#if 1
	if ( a < HALF_PI ) {
		if ( a > ONEFOURTH_PI ) {
			a = HALF_PI - a;
			reciprocal = true;
		} else {
			reciprocal = false;
		}
	} else {
		if ( a > HALF_PI + ONEFOURTH_PI ) {
			a = a - PI;
			reciprocal = false;
		} else {
			a = HALF_PI - a;
			reciprocal = true;
		}
	}
#else
	a = HALF_PI - a;
	if ( st::fabs( a ) >= ONEFOURTH_PI ) {
		a = ( ( a < 0.0f ) ? -HALF_PI : HALF_PI ) - a;
		reciprocal = false;
	} else {
		reciprocal = true;
	}
#endif
	s = a * a;
	s = a * ( ( ( ( ( ( 9.5168091e-03f * s + 2.900525e-03f ) * s + 2.45650893e-02f ) * s + 5.33740603e-02f ) * s + 1.333923995e-01f ) * s + 3.333314036e-01f ) * s + 1.0f );
	if ( reciprocal ) {
		return 1.0f / s;
	} else {
		return s;
	}
}

SMF_INLINE_FORCED double CMath::tan64( double a ) {
	return std::tan( a );
}

SMF_INLINE_FORCED float CMath::asin( float a ) {
	if ( a <= -1.0f ) {
		return -HALF_PI;
	}
	if ( a >= 1.0f ) {
		return HALF_PI;
	}
	return asinf( a );
}

SMF_INLINE_FORCED float CMath::asin16( float a ) {
	if ( FLOATSIGNBITSET( a ) ) {
		if ( a <= -1.0f ) {
			return -HALF_PI;
		}
		a = std::fabs( a );
		return ( ( ( -0.0187293f * a + 0.0742610f ) * a - 0.2121144f ) * a + 1.5707288f ) * sqrt( 1.0f - a ) - HALF_PI;
	} else {
		if ( a >= 1.0f ) {
			return HALF_PI;
		}
		return HALF_PI - ( ( ( -0.0187293f * a + 0.0742610f ) * a - 0.2121144f ) * a + 1.5707288f ) * sqrt( 1.0f - a );
	}
}

SMF_INLINE_FORCED double CMath::asin64( double a ) {
	if ( a <= -1.0f ) {
		return -HALF_PI;
	}
	if ( a >= 1.0f ) {
		return HALF_PI;
	}
	return std::asin( a );
}

SMF_INLINE_FORCED float CMath::acos( float a ) {
	if ( a <= -1.0f ) {
		return PI;
	}
	if ( a >= 1.0f ) {
		return 0.0f;
	}
	return acosf( a );
}

SMF_INLINE_FORCED float CMath::acos16( float a ) {
	if ( FLOATSIGNBITSET( a ) ) {
		if ( a <= -1.0f ) {
			return PI;
		}
		a = std::fabs( a );
		return PI - ( ( ( -0.0187293f * a + 0.0742610f ) * a - 0.2121144f ) * a + 1.5707288f ) * sqrt( 1.0f - a );
	} else {
		if ( a >= 1.0f ) {
			return 0.0f;
		}
		return ( ( ( -0.0187293f * a + 0.0742610f ) * a - 0.2121144f ) * a + 1.5707288f ) * sqrt( 1.0f - a );
	}
}

SMF_INLINE_FORCED double CMath::acos64( double a ) {
	if ( a <= -1.0f ) {
		return PI;
	}
	if ( a >= 1.0f ) {
		return 0.0f;
	}
	return std::acos( a );
}

SMF_INLINE_FORCED float CMath::atan( float a ) {
	return atanf( a );
}

SMF_INLINE_FORCED float CMath::atan16( float a ) {
	float s;

	if ( std::fabs( a ) > 1.0f ) {
		a = 1.0f / a;
		s = a * a;
		s = - ( ( ( ( ( ( ( ( ( 0.0028662257f * s - 0.0161657367f ) * s + 0.0429096138f ) * s - 0.0752896400f )
				* s + 0.1065626393f ) * s - 0.1420889944f ) * s + 0.1999355085f ) * s - 0.3333314528f ) * s ) + 1.0f ) * a;
		if ( FLOATSIGNBITSET( a ) ) {
			return s - HALF_PI;
		} else {
			return s + HALF_PI;
		}
	} else {
		s = a * a;
		return ( ( ( ( ( ( ( ( ( 0.0028662257f * s - 0.0161657367f ) * s + 0.0429096138f ) * s - 0.0752896400f )
			* s + 0.1065626393f ) * s - 0.1420889944f ) * s + 0.1999355085f ) * s - 0.3333314528f ) * s ) + 1.0f ) * a;
	}
}

SMF_INLINE_FORCED double CMath::atan64( double a ) {
	return std::atan( a );
}

SMF_INLINE_FORCED float CMath::atan( float y, float x ) {
	SMF_ASSERT( fabs( y ) > CMath::FLT_SMALLEST_NON_DENORMAL || fabs( x ) > CMath::FLT_SMALLEST_NON_DENORMAL );

	return atan2f( y, x );
}

SMF_INLINE_FORCED float CMath::atan16( float y, float x ) {
	SMF_ASSERT( fabs( y ) > CMath::FLT_SMALLEST_NON_DENORMAL || fabs( x ) > CMath::FLT_SMALLEST_NON_DENORMAL );

	float a, s;

	if ( std::fabs( y ) > fabs( x ) ) {
		a = x / y;
		s = a * a;
		s = - ( ( ( ( ( ( ( ( ( 0.0028662257f * s - 0.0161657367f ) * s + 0.0429096138f ) * s - 0.0752896400f )
				* s + 0.1065626393f ) * s - 0.1420889944f ) * s + 0.1999355085f ) * s - 0.3333314528f ) * s ) + 1.0f ) * a;
		if ( FLOATSIGNBITSET( a ) ) {
			return s - HALF_PI;
		} else {
			return s + HALF_PI;
		}
	} else {
		a = y / x;
		s = a * a;
		return ( ( ( ( ( ( ( ( ( 0.0028662257f * s - 0.0161657367f ) * s + 0.0429096138f ) * s - 0.0752896400f )
			* s + 0.1065626393f ) * s - 0.1420889944f ) * s + 0.1999355085f ) * s - 0.3333314528f ) * s ) + 1.0f ) * a;
	}
}

SMF_INLINE_FORCED double CMath::atan64( double y, double x ) {
	return std::atan2( y, x );
}

//hyperbolic===============

SMF_INLINE_FORCED float CMath::log1p (const float x)
{
  volatile float y;
  y = 1 + x;
  return CMath::log(y) - ((y-1)-x)/y ;  /* cancels errors with IEEE arithmetic */
}

SMF_INLINE_FORCED float CMath::acosh (const float x)
{
  if (x > 1.0 / FLT_SQRT_EPSILON)
    {
      return CMath::log (x) + M_LN2;
    }
  else if (x > 2)
    {
      return CMath::log (2 * x - 1 / (CMath::sqrt (x * x - 1) + x));
    }
  else if (x > 1)
    {
      float t = x - 1;
      return CMath::log1p (t + CMath::sqrt (2 * t + t * t));
    }
  else if (x == 1)
    {
      return 0;
    }
  else
    {
      return CMath::NAN_FLOAT;
    }
}

SMF_INLINE_FORCED float CMath::asinh (const float x)
{
  float a = CMath::fabs(x);
  float s = (x < 0) ? -1 : 1;

  if (a > 1 / FLT_SQRT_EPSILON)
    {
      return s * (CMath::log (a) + M_LN2);
    }
  else if (a > 2)
    {
      return s * CMath::log (2 * a + 1 / (a + CMath::sqrt (a * a + 1)));
    }
  else if (a > FLT_SQRT_EPSILON)
    {
      float a2 = a * a;
      return s * CMath::log1p (a + a2 / (1 + CMath::sqrt (1 + a2)));
    }
  else
    {
      return x;
    }
}

SMF_INLINE_FORCED float CMath::atanh (const float x)
{
  float a = CMath::fabs (x);
  float s = (x < 0) ? -1 : 1;

  if (a > 1)
    {
      return CMath::NAN_FLOAT;
    }
  else if (a == 1)
    {
      return (x < 0) ? CMath::NEG_INFINITY_FLOAT : CMath::INFINITY_FLOAT;
    }
  else if (a >= 0.5)
    {
      return s * 0.5 * CMath::log1p (2 * a / (1 - a));
    }
  else if (a > CMath::FLT_EPSILON)
    {
      return s * 0.5 * CMath::log1p (2 * a + 2 * a * a / (1 - a));
    }
  else
    {
      return x;
    }
}

SMF_INLINE_FORCED float CMath::pow( float x, float y ) {
	return powf( x, y );
}

SMF_INLINE_FORCED float CMath::pow16( float x, float y ) {
	return exp16( y * log16( x ) );
}

SMF_INLINE_FORCED double CMath::pow64( double x, double y ) {
	return std::pow( x, y );
}

SMF_INLINE_FORCED float CMath::square( float x ) {
	return powf( x, 2 );
}

SMF_INLINE_FORCED float CMath::square16( float x ) {
	return exp16( 2 * log16( x ) );
}

SMF_INLINE_FORCED double CMath::square64( double x ) {
	return std::pow( x, 2 );
}
SMF_INLINE_FORCED float CMath::exp( float f ) {
	return expf( f );
}

SMF_INLINE_FORCED float CMath::exp16( float f ) {
	int i, s, e, m, exponent;
	float x, x2, y, p, q;

	x = f * 1.44269504088896340f;		// multiply with ( 1 / log( 2 ) )
#if 1
	i = *reinterpret_cast<int *>(&x);
	s = ( i >> IEEE_FLT_SIGN_BIT );
	e = ( ( i >> IEEE_FLT_MANTISSA_BITS ) & ( ( 1 << IEEE_FLT_EXPONENT_BITS ) - 1 ) ) - IEEE_FLT_EXPONENT_BIAS;
	m = ( i & ( ( 1 << IEEE_FLT_MANTISSA_BITS ) - 1 ) ) | ( 1 << IEEE_FLT_MANTISSA_BITS );
	i = ( ( m >> ( IEEE_FLT_MANTISSA_BITS - e ) ) & ~( e >> 31 ) ) ^ s;
#else
	i = (int) x;
	if ( x < 0.0f ) {
		i--;
	}
#endif
	exponent = ( i + IEEE_FLT_EXPONENT_BIAS ) << IEEE_FLT_MANTISSA_BITS;
	y = *reinterpret_cast<float *>(&exponent);
	x -= (float) i;
	if ( x >= 0.5f ) {
		x -= 0.5f;
		y *= 1.4142135623730950488f;	// multiply with sqrt( 2 )
	}
	x2 = x * x;
	p = x * ( 7.2152891511493f + x2 * 0.0576900723731f );
	q = 20.8189237930062f + x2;
	x = y * ( q + p ) / ( q - p );
	return x;
}

SMF_INLINE_FORCED double CMath::exp64( double f ) {
	return std::exp( f );
}

SMF_INLINE_FORCED float CMath::log( float f ) {
	return logf( f );
}

SMF_INLINE_FORCED float CMath::log16( float f ) {
	int i, exponent;
	float y, y2;

	i = *reinterpret_cast<int *>(&f);
	exponent = ( ( i >> IEEE_FLT_MANTISSA_BITS ) & ( ( 1 << IEEE_FLT_EXPONENT_BITS ) - 1 ) ) - IEEE_FLT_EXPONENT_BIAS;
	i -= ( exponent + 1 ) << IEEE_FLT_MANTISSA_BITS;	// get value in the range [.5, 1>
	y = *reinterpret_cast<float *>(&i);
	y *= 1.4142135623730950488f;						// multiply with sqrt( 2 )
	y = ( y - 1.0f ) / ( y + 1.0f );
	y2 = y * y;
	y = y * ( 2.000000000046727f + y2 * ( 0.666666635059382f + y2 * ( 0.4000059794795f + y2 * ( 0.28525381498f + y2 * 0.2376245609f ) ) ) );
	y += M_LN2 * ( (float)exponent + 0.5f );
	return y;
}

SMF_INLINE_FORCED double CMath::log64( double f ) {
	return std::log( f );
}

SMF_INLINE_FORCED int CMath::ipow( int x, int y ) {
	int r; for( r = x; y > 1; y-- ) { r *= x; } return r;
}

SMF_INLINE_FORCED int CMath::ilog2( float f ) {
	return ( ( (*reinterpret_cast<int *>(&f)) >> IEEE_FLT_MANTISSA_BITS ) & ( ( 1 << IEEE_FLT_EXPONENT_BITS ) - 1 ) ) - IEEE_FLT_EXPONENT_BIAS;
}

SMF_INLINE_FORCED int CMath::ilog2( int i ) {
	return ilog2( (float)i );
}

SMF_INLINE_FORCED int CMath::bitsForFloat( float f ) {
	return ilog2( f ) + 1;
}

SMF_INLINE_FORCED int CMath::bitsForInteger( int i ) {
	return ilog2( (float)i ) + 1;
}

SMF_INLINE_FORCED int CMath::maskForFloatSign( float f ) {
	return ( (*reinterpret_cast<int *>(&f)) >> 31 );
}

SMF_INLINE_FORCED int CMath::maskForIntegerSign( int i ) {
	return ( i >> 31 );
}

SMF_INLINE_FORCED int CMath::floorPowerOfTwo( int x ) {
	return ceilPowerOfTwo( x ) >> 1;
}

SMF_INLINE_FORCED int CMath::ceilPowerOfTwo( int x ) {
	x--;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	x++;
	return x;
}

SMF_INLINE_FORCED bool CMath::isPowerOfTwo( int x ) {
	return ( x & ( x - 1 ) ) == 0 && x > 0;
}

SMF_INLINE_FORCED int CMath::bitCount( int x ) {
	x -= ( ( x >> 1 ) & 0x55555555 );
	x = ( ( ( x >> 2 ) & 0x33333333 ) + ( x & 0x33333333 ) );
	x = ( ( ( x >> 4 ) + x ) & 0x0f0f0f0f );
	x += ( x >> 8 );
	return ( ( x + ( x >> 16 ) ) & 0x0000003f );
}

SMF_INLINE_FORCED int CMath::bitReverse( int x ) {
	x = ( ( ( x >> 1 ) & 0x55555555 ) | ( ( x & 0x55555555 ) << 1 ) );
	x = ( ( ( x >> 2 ) & 0x33333333 ) | ( ( x & 0x33333333 ) << 2 ) );
	x = ( ( ( x >> 4 ) & 0x0f0f0f0f ) | ( ( x & 0x0f0f0f0f ) << 4 ) );
	x = ( ( ( x >> 8 ) & 0x00ff00ff ) | ( ( x & 0x00ff00ff ) << 8 ) );
	return ( ( x >> 16 ) | ( x << 16 ) );
}

SMF_INLINE_FORCED int CMath::abs( int x ) {
   int y = x >> 31;
   return ( ( x ^ y ) - y );
}

SMF_INLINE_FORCED float CMath::fabs( float f ) {
	int tmp = *reinterpret_cast<int *>( &f );
	tmp &= 0x7FFFFFFF;
	return *reinterpret_cast<float *>( &tmp );
}

SMF_INLINE_FORCED float CMath::floor( float f ) {
	return floorf( f );
}

SMF_INLINE_FORCED float CMath::ceil( float f ) {
	return ceilf( f );
}

SMF_INLINE_FORCED float CMath::rint( float f ) {
	return floorf( f + 0.5f );
}

SMF_INLINE_FORCED float CMath::frac( float f ) {
	return f - floorf( f );
}

SMF_INLINE_FORCED int CMath::ftoi( float f ) {
	return (int) f;
}

SMF_INLINE_FORCED int CMath::ftoiFast( float f ) {
#if defined(_MSC_VER)
	int i;
	__asm fld		f
	__asm fistp		i		// use default rouding mode (round nearest)
	return i;
#elif 0						// round chop (C/C++ standard)
	int i, s, e, m, shift;
	i = *reinterpret_cast<int *>(&f);
	s = i >> IEEE_FLT_SIGN_BIT;
	e = ( ( i >> IEEE_FLT_MANTISSA_BITS ) & ( ( 1 << IEEE_FLT_EXPONENT_BITS ) - 1 ) ) - IEEE_FLT_EXPONENT_BIAS;
	m = ( i & ( ( 1 << IEEE_FLT_MANTISSA_BITS ) - 1 ) ) | ( 1 << IEEE_FLT_MANTISSA_BITS );
	shift = e - IEEE_FLT_MANTISSA_BITS;
	return ( ( ( ( m >> -shift ) | ( m << shift ) ) & ~( e >> 31 ) ) ^ s ) - s;
//#elif defined( __i386__ ) || defined (__i686__)
#elif 0
	int i = 0;
	__asm__ __volatile__ (
						  "fld %1\n" \
						  "fistp %0\n" \
						  : "=m" (i) \
						  : "m" (f) );
	return i;
#else
	return (int) f;
#endif
}

/*
========================
idMath::ftoi8
========================
*/
SMF_INLINE char CMath::ftoi8( float f ) {
#ifdef ID_WIN_X86_SSE_INTRIN
	sf_m128 x = _mm_load_ss( &f );
	x = _mm_max_ss( x, SIMD_SP_min_char );
	x = _mm_min_ss( x, SIMD_SP_max_char );
	return static_cast<char>( _mm_cvttss_si32( x ) );
#else
	// The converted result is clamped to the range [-128,127].
	int i = C_FLOAT_TO_INT( f );
	if ( i < -128 ) {
		return -128;
	} else if ( i > 127 ) {
		return 127;
	}
	return static_cast<char>( i );
#endif
}

/*
========================
idMath::ftoi16
========================
*/
SMF_INLINE short CMath::ftoi16( float f ) {
#ifdef ID_WIN_X86_SSE_INTRIN
	sf_m128 x = _mm_load_ss( &f );
	x = _mm_max_ss( x, SIMD_SP_min_short );
	x = _mm_min_ss( x, SIMD_SP_max_short );
	return static_cast<short>( _mm_cvttss_si32( x ) );
#else
	// The converted result is clamped to the range [-32768,32767].
	int i = C_FLOAT_TO_INT( f );
	if ( i < -32768 ) {
		return -32768;
	} else if ( i > 32767 ) {
		return 32767;
	}
	return static_cast<short>( i );
#endif
}

/*
========================
idMath::ftoui16
========================
*/
SMF_INLINE unsigned short CMath::ftoui16( float f ) {
	// TO DO - SSE ??

	// The converted result is clamped to the range [-32768,32767].
	int i = C_FLOAT_TO_INT( f );
	if ( i < 0 ) {
		return 0;
	} else if ( i > 65535 ) {
		return 65535;
	}
	return static_cast<unsigned short>( i );
}

/*
========================
idMath::ftob
========================
*/
SMF_INLINE sf_u8 CMath::ftob( float f ) {
#ifdef ID_WIN_X86_SSE_INTRIN
	// If a converted result is negative the value (0) is returned and if the
	// converted result is larger than the maximum byte the value (255) is returned.
	sf_m128 x = _mm_load_ss( &f );
	x = _mm_max_ss( x, SIMD_SP_zero );
	x = _mm_min_ss( x, SIMD_SP_255 );
	return static_cast<byte>( _mm_cvttss_si32( x ) );
#else
	// The converted result is clamped to the range [0,255].
	int i = C_FLOAT_TO_INT( f );
	if ( i < 0 ) {
		return 0;
	} else if ( i > 255 ) {
		return 255;
	}
	return static_cast<sf_u8>( i );
#endif
}


SMF_INLINE_FORCED unsigned long CMath::ftol( float f ) {
	return (unsigned long) f;
}

SMF_INLINE_FORCED unsigned long CMath::ftolFast( float f ) {
#if defined(_MSC_VER)
	// FIXME: this overflows on 31bits still .. same as ftoiFast
	unsigned long i;
	__asm fld		f
	__asm fistp		i		// use default rouding mode (round nearest)
	return i;
#elif 0						// round chop (C/C++ standard)
	int i, s, e, m, shift;
	i = *reinterpret_cast<int *>(&f);
	s = i >> IEEE_FLT_SIGN_BIT;
	e = ( ( i >> IEEE_FLT_MANTISSA_BITS ) & ( ( 1 << IEEE_FLT_EXPONENT_BITS ) - 1 ) ) - IEEE_FLT_EXPONENT_BIAS;
	m = ( i & ( ( 1 << IEEE_FLT_MANTISSA_BITS ) - 1 ) ) | ( 1 << IEEE_FLT_MANTISSA_BITS );
	shift = e - IEEE_FLT_MANTISSA_BITS;
	return ( ( ( ( m >> -shift ) | ( m << shift ) ) & ~( e >> 31 ) ) ^ s ) - s;
//#elif defined( __i386__ ) || defined (__i686__)
#elif 0
	// for some reason, on gcc I need to make sure i == 0 before performing a fistp
	int i = 0;
	__asm__ __volatile__ (
						  "fld %1\n" \
						  "fistp %0\n" \
						  : "=m" (i) \
						  : "m" (f) );
	return i;
#else
	return (unsigned long) f;
#endif
}


SMF_INLINE_FORCED int CMath::ceilInt(float x)
{
	return (int)ceilf(x);
}


SMF_INLINE_FORCED int CMath::floorInt(float x)
{
	return (int)floorf(x);
}

SMF_INLINE_FORCED float CMath::round(float x)
{
	return floor(x+0.5f);
}
SMF_INLINE_FORCED float CMath::round(float f,float prec)
{
	if (!isValid<float>(f)) return f;
    return (float)(static_cast<long int> ((f)*(pow(10,prec))+0.5)) / pow(10,prec);
}
SMF_INLINE_FORCED int CMath::roundInt(float x)
{
	return (int)round(x);
}

SMF_INLINE_FORCED float CMath::sign(float x)
{
	return x >= 0.f ? 1.f : -1.f;
}

SMF_INLINE_FORCED float CMath::signOrZero(float x, float epsilon)
{
	return abs(x) <= epsilon ? 0.f : sign(x);
}


SMF_INLINE_FORCED signed char CMath::clampChar( int i ) {
	if ( i < -128 ) {
		return -128;
	}
	if ( i > 127 ) {
		return 127;
	}
	return i;
}

SMF_INLINE_FORCED signed short CMath::clampShort( int i ) {
	if ( i < -32768 ) {
		return -32768;
	}
	if ( i > 32767 ) {
		return 32767;
	}
	return i;
}

SMF_INLINE_FORCED int CMath::clampInt( int min, int max, int value ) {
	if ( value < min ) {
		return min;
	}
	if ( value > max ) {
		return max;
	}
	return value;
}

SMF_INLINE_FORCED float CMath::clampFloat( float min, float max, float value ) {
	if ( value < min ) {
		return min;
	}
	if ( value > max ) {
		return max;
	}
	return value;
}

SMF_INLINE_FORCED float CMath::angleNormalize360( float angle ) {
	if ( ( angle >= 360.0f ) || ( angle < 0.0f ) ) {
		angle -= std::floor( angle / 360.0f ) * 360.0f;
	}
	return angle;
}

SMF_INLINE_FORCED float CMath::angleNormalize180( float angle ) {
	angle = angleNormalize360( angle );
	if ( angle > 180.0f ) {
		angle -= 360.0f;
	}
	return angle;
}

SMF_INLINE_FORCED float CMath::angleDelta( float angle1, float angle2 ) {
	return angleNormalize180( angle1 - angle2 );
}

SMF_INLINE_FORCED int CMath::floatHash( const float *array, const int numFloats ) {
	int i, hash = 0;
	const int *ptr;

	ptr = reinterpret_cast<const int *>( array );
	for ( i = 0; i < numFloats; i++ ) {
		hash ^= ptr[i];
	}
	return hash;
}
} //end MATH
} //end SMF
#endif /* !__MATH_MATH_H__ */
