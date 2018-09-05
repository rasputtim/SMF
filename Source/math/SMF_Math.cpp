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

#include "SMF_Config.h"
#include "math/SMF_Math.h"
#include "exceptions/all.h"
#include "sys/SMF_System.h"
#include <stdlib.h>
namespace SMF{
namespace MATH{

sf_s32 ReinterpretAsU32(float f)
{
	FloatIntReinterpret fi;
	fi.f = f;
	return fi.i;
}

sf_u64 ReinterpretAsU64(double d)
{
	DoubleU64Reinterpret di;
	di.d = d;
	return di.i;
};

/// Converts the bit pattern specified by the given integer to a floating point (this is a binary conversion, not numeral!).
float ReinterpretAsFloat(sf_s32 i)
{
	FloatIntReinterpret fi;
	fi.i = i;
	return fi.f;
};

double ReinterpretAsDouble(sf_u64 i)
{
	DoubleU64Reinterpret di;
	di.i = i;
	return di.d;
};


bool CMath::MATH_AUTOMATIC_SSE = true;

const int SMALLEST_NON_DENORMAL					= 1<<IEEE_FLT_MANTISSA_BITS;

const float	CMath::PI				= 3.14159265358979323846f;
const float	CMath::TWO_PI			= 2.0f * PI;
const float	CMath::HALF_PI			= 0.5f * PI;
const float	CMath::ONEFOURTH_PI	= 0.25f * PI;
const float CMath::E				= 2.71828182845904523536f;
const float CMath::SQRT_TWO		= 1.41421356237309504880f;
const float CMath::SQRT_THREE		= 1.73205080756887729352f;
const float	CMath::SQRT_1OVER2		= 0.70710678118654752440f;
const float	CMath::SQRT_1OVER3		= 0.57735026918962576450f;
const float	CMath::M_DEG2RAD		= PI / 180.0f;
const float	CMath::M_RAD2DEG		= 180.0f / PI;
const float	CMath::M_SEC2MS		= 1000.0f;
const float	CMath::M_MS2SEC		= 0.001f;

#if defined(_MSC_VER) || defined(EMSCRIPTEN)
#define _FLOAT_NAN ((float)std::numeric_limits<float>::quiet_NaN())
#define _INFINITY_FLOAT ((float)std::numeric_limits<float>::infinity())
#else
#define _FLOAT_NAN ((float)NAN)
#define _INFINITY_FLOAT ((float)INFINITY_FLOAT)
#endif

const float	CMath::INFINITY_FLOAT		= _INFINITY_FLOAT;
const float	CMath::NEG_INFINITY_FLOAT	= -_INFINITY_FLOAT;
const float	CMath::NAN_FLOAT			= _FLOAT_NAN;

const float	CMath::EPSILON_High      = 1.0E-16;
const float CMath::EPSILON_Medium    = 1.0E-12;
const float CMath::EPSILON_Low       = 1.0E-08;
const float CMath::EPSILON_SuperLow  = 1.0E-04;
const float CMath::FLT_EPSILON		=  EPSILON_Low;
const float	CMath::ZERO              = 0.0;
const float CMath::FLT_SMALLEST_NON_DENORMAL	= * reinterpret_cast< const float * >( & SMALLEST_NON_DENORMAL );	// 1.1754944e-038f

bool		CMath::initialized		= false;
sf_u32		CMath::iSqrt[SQRT_TABLE_SIZE];		// inverse square root lookup table

const float NOT_NECESSARILY_USED CMath::inf =          CMath::INFINITY_FLOAT;
const float NOT_NECESSARILY_USED CMath::negInf =       -CMath::INFINITY_FLOAT;
const float NOT_NECESSARILY_USED CMath::nan =          CMath::NAN_FLOAT;


/*
===============
CMath::init
===============
*/
void CMath::init() {
    union _flint fi, fo;

    for ( int i = 0; i < SQRT_TABLE_SIZE; i++ ) {
        fi.i	 = ((EXP_BIAS-1) << EXP_POS) | (i << LOOKUP_POS);
        fo.f	 = (float)( 1.0 / std::sqrt( fi.f ) );
        iSqrt[i] = ((sf_u32)(((fo.i + (1<<(SEED_POS-2))) >> SEED_POS) & 0xFF))<<SEED_POS;
    }

	iSqrt[SQRT_TABLE_SIZE / 2] = ((sf_u32)(0xFF))<<(SEED_POS);

	initialized = true;
}
bool	CMath::isInitialized(){
	return initialized;
}
/*
================
CMath::floatToBits
================
*/
int CMath::floatToBits( float f, int exponentBits, int mantissaBits ) {
	int i, sign, exponent, mantissa, value;

	SMF_ASSERT( exponentBits >= 2 && exponentBits <= 8 );
	SMF_ASSERT( mantissaBits >= 2 && mantissaBits <= 23 );

	int maxBits = ( ( ( 1 << ( exponentBits - 1 ) ) - 1 ) << mantissaBits ) | ( ( 1 << mantissaBits ) - 1 );
	int minBits = ( ( ( 1 <<   exponentBits       ) - 2 ) << mantissaBits ) | 1;

	float max = bitsToFloat( maxBits, exponentBits, mantissaBits );
	float min = bitsToFloat( minBits, exponentBits, mantissaBits );

	if ( f >= 0.0f ) {
		if ( f >= max ) {
			return maxBits;
		} else if ( f <= min ) {
			return minBits;
		}
	} else {
		if ( f <= -max ) {
			return ( maxBits | ( 1 << ( exponentBits + mantissaBits ) ) );
		} else if ( f >= -min ) {
			return ( minBits | ( 1 << ( exponentBits + mantissaBits ) ) );
		}
	}

	exponentBits--;
	i = *reinterpret_cast<int *>(&f);
	sign = ( i >> IEEE_FLT_SIGN_BIT ) & 1;
	exponent = ( ( i >> IEEE_FLT_MANTISSA_BITS ) & ( ( 1 << IEEE_FLT_EXPONENT_BITS ) - 1 ) ) - IEEE_FLT_EXPONENT_BIAS;
	mantissa = i & ( ( 1 << IEEE_FLT_MANTISSA_BITS ) - 1 );
	value = sign << ( 1 + exponentBits + mantissaBits );
	value |= ( ( INTSIGNBITSET( exponent ) << exponentBits ) | ( std::abs( exponent ) & ( ( 1 << exponentBits ) - 1 ) ) ) << mantissaBits;
	value |= mantissa >> ( IEEE_FLT_MANTISSA_BITS - mantissaBits );
	return value;
}

/*
================
CMath::bitsToFloat
================
*/
float CMath::bitsToFloat( int i, int exponentBits, int mantissaBits ) {
	static int exponentSign[2] = { 1, -1 };
	int sign, exponent, mantissa, value;

	SMF_ASSERT( exponentBits >= 2 && exponentBits <= 8 );
	SMF_ASSERT( mantissaBits >= 2 && mantissaBits <= 23 );

	exponentBits--;
	sign = i >> ( 1 + exponentBits + mantissaBits );
	exponent = ( ( i >> mantissaBits ) & ( ( 1 << exponentBits ) - 1 ) ) * exponentSign[( i >> ( exponentBits + mantissaBits ) ) & 1];
	mantissa = ( i & ( ( 1 << mantissaBits ) - 1 ) ) << ( IEEE_FLT_MANTISSA_BITS - mantissaBits );
	value = sign << IEEE_FLT_SIGN_BIT | ( exponent + IEEE_FLT_EXPONENT_BIAS ) << IEEE_FLT_MANTISSA_BITS | mantissa;
	return *reinterpret_cast<float *>(&value);
}

/** Compares the two values for equality, allowing the given amount of absolute error. */
bool CMath::equalsAbs(float a, float b, float epsilon)
{
	return abs(a-b) < epsilon;
}
bool CMath::equals(float a, float b){
return nearZero(a-b);
}

bool CMath::nearZero(float a){
   return a < FLT_EPSILON && a > -FLT_EPSILON;

}

//===IEEE================================

static int lookup_string (const char * p, int * precision, int * rounding, int * exception_mask)
{
  if (strcmp(p,"single-precision") == 0)
    {
      *precision = SMF_IEEE_SINGLE_PRECISION ;
    }
  else if (strcmp(p,"double-precision") == 0)
    {
      *precision = SMF_IEEE_DOUBLE_PRECISION ;
    }
  else if (strcmp(p,"extended-precision") == 0)
    {
      *precision = SMF_IEEE_EXTENDED_PRECISION ;
    }
  else if (strcmp(p,"round-to-nearest") == 0)
    {
      *rounding = SMF_IEEE_ROUND_TO_NEAREST ;
    }
  else if (strcmp(p,"round-down") == 0)
    {
      *rounding = SMF_IEEE_ROUND_DOWN ;
    }
  else if (strcmp(p,"round-up") == 0)
    {
      *rounding = SMF_IEEE_ROUND_UP ;
    }
  else if (strcmp(p,"round-to-zero") == 0)
    {
      *rounding = SMF_IEEE_ROUND_TO_ZERO ;
    }
  else if (strcmp(p,"mask-all") == 0)
    {
      *exception_mask = SMF_IEEE_MASK_ALL ;
    }
  else if (strcmp(p,"mask-invalid") == 0)
    {
      *exception_mask = SMF_IEEE_MASK_INVALID ;
    }
  else if (strcmp(p,"mask-denormalized") == 0)
    {
      *exception_mask = SMF_IEEE_MASK_DENORMALIZED ;
    }
  else if (strcmp(p,"mask-division-by-zero") == 0)
    {
      *exception_mask = SMF_IEEE_MASK_DIVISION_BY_ZERO ;
    }
  else if (strcmp(p,"mask-overflow") == 0)
    {
      *exception_mask = SMF_IEEE_MASK_OVERFLOW ;
    }
  else if (strcmp(p,"mask-underflow") == 0)
    {
      *exception_mask = SMF_IEEE_MASK_UNDERFLOW ;
    }
  else if (strcmp(p,"trap-inexact") == 0)
    {
      *exception_mask = SMF_IEEE_TRAP_INEXACT ;
    }
  else if (strcmp(p,"trap-common") == 0)
    {
      return 0 ;
    }
  else
    {
      return 1 ;
    }

  return 0 ;
}

int ieee_read_mode_string (const char * description,
                           int * precision,
                           int * rounding,
                           int * exception_mask)
{
  char * start ;
  char * end;
  char * p;

  int precision_count = 0 ;
  int rounding_count = 0 ;
  int exception_count = 0 ;

  start = (char *) malloc(strlen(description) + 1) ;

  if (start == 0)
    {
		Debug::debug(Debug::error,__FUNCTION__) << "no memory to parse mode string"<<endl ;
    }

  strcpy (start, description) ;

  p = start ;

  *precision = 0 ;
  *rounding = 0 ;
  *exception_mask = 0 ;

  do {
    int status ;
    int new_precision, new_rounding, new_exception ;

    end = strchr (p,',') ;

    if (end)
      {
        *end = '\0' ;
        do
          {
            end++ ;  /* skip over trailing whitespace */
          }
        while (*end == ' ' || *end == ',') ;
      }

    new_precision = 0 ;
    new_rounding = 0 ;
    new_exception = 0 ;

    status = lookup_string (p, &new_precision, &new_rounding, &new_exception) ;

    if (status)
      Debug::debug(Debug::error,__FUNCTION__) <<"unrecognized SMF_IEEE_MODE string.\nValid settings are:\n\n"
                 "  single-precision double-precision extended-precision\n"
                 "  round-to-nearest round-down round-up round-to-zero\n"
                 "  mask-invalid mask-denormalized mask-division-by-zero\n"
                 "  mask-overflow mask-underflow mask-all\n"
                 "  trap-common trap-inexact\n"
                 "\n"
                 "separated by commas. "
                 "(e.g. SMF_IEEE_MODE=\"round-down,mask-underflow\")"<<
                 endl ;

    if (new_precision)
      {
        *precision = new_precision ;
        precision_count ++ ;
        if (precision_count > 1)
          Debug::debug(Debug::error,__FUNCTION__) <<"attempted to set IEEE precision twice" <<endl;
      }

    if (new_rounding)
      {
        *rounding = new_rounding ;
        rounding_count ++ ;
        if (rounding_count > 1)
          Debug::debug(Debug::error,__FUNCTION__) <<"attempted to set IEEE rounding mode twice" <<endl ;
      }

    if (new_exception)
      {
        *exception_mask |= new_exception ;
        exception_count ++ ;
      }

    p = end ;

  } while (end && *p != '\0') ;

  free(start) ;

  return SMF_OK ;
}
void ieee_env_setup (void)
{
  const char * p = getenv("SMF_IEEE_MODE") ;

  int precision = 0, rounding = 0, exception_mask = 0 ;

  int comma = 0 ;

  if (p == 0)  /* SMF_IEEE_MODE environment variable is not set */
    return ;

  if (*p == '\0') /* SMF_IEEE_MODE environment variable is empty */
    return ;

  MATH::ieee_read_mode_string (p, &precision, &rounding, &exception_mask) ;

  System::Sys_ieee_set_FPU_mode (precision, rounding, exception_mask) ;

  Debug::debug(Debug::error,__FUNCTION__) << "SMF_IEEE_MODE=\"" << endl ;

  /* Print string with a preceeding comma if the list has already begun */

#define PRINTC(x) do {if(comma) Debug::debug(Debug::error,__FUNCTION__) <<","; Debug::debug(Debug::error,__FUNCTION__) << x; comma++ ;} while(0)

  switch (precision)
    {
    case SMF_IEEE_SINGLE_PRECISION:
      PRINTC("single-precision") ;
      break ;
    case SMF_IEEE_DOUBLE_PRECISION:
      PRINTC("double-precision") ;
      break ;
    case SMF_IEEE_EXTENDED_PRECISION:
      PRINTC("extended-precision") ;
      break ;
    }

  switch (rounding)
    {
    case SMF_IEEE_ROUND_TO_NEAREST:
      PRINTC("round-to-nearest") ;
      break ;
    case SMF_IEEE_ROUND_DOWN:
      PRINTC("round-down") ;
      break ;
    case SMF_IEEE_ROUND_UP:
      PRINTC("round-up") ;
      break ;
    case SMF_IEEE_ROUND_TO_ZERO:
      PRINTC("round-to-zero") ;
      break ;
    }

  if ((exception_mask & SMF_IEEE_MASK_ALL) == SMF_IEEE_MASK_ALL)
    {
      PRINTC("mask-all") ;
    }
  else if ((exception_mask & SMF_IEEE_MASK_ALL) == 0)
    {
      PRINTC("trap-common") ;
    }
  else
    {
      if (exception_mask & SMF_IEEE_MASK_INVALID)
        PRINTC("mask-invalid") ;

      if (exception_mask & SMF_IEEE_MASK_DENORMALIZED)
        PRINTC("mask-denormalized") ;

      if (exception_mask & SMF_IEEE_MASK_DIVISION_BY_ZERO)
        PRINTC("mask-division-by-zero") ;

      if (exception_mask & SMF_IEEE_MASK_OVERFLOW)
        PRINTC("mask-overflow") ;

      if (exception_mask & SMF_IEEE_MASK_UNDERFLOW)
        PRINTC("mask-underflow") ;
    }

  if (exception_mask & SMF_IEEE_TRAP_INEXACT)
    PRINTC("trap-inexact") ;

  Debug::debug(Debug::error,__FUNCTION__) <<"\\"<< endl ;
}
/* A table of sign characters, 0=positive, 1=negative. We print a space
   instead of a unary + sign for compatibility with bc */

static char signs[2]={' ','-'} ;
/* A table of character representations of nybbles */

static char nybble[16][5]={ /* include space for the \0 */
  "0000", "0001", "0010", "0011",
  "0100", "0101", "0110", "0111",
  "1000", "1001", "1010", "1011",
  "1100", "1101", "1110", "1111"
}  ;

static void sprint_nybble(int i, char *s)
{
  char *c ;
  c=nybble[i & 0x0f ];
  *s=c[0] ;  *(s+1)=c[1] ;  *(s+2)=c[2] ;  *(s+3)=c[3] ;
}

static void
sprint_byte(int i, char *s)
{
  char *c ;
  c=nybble[(i & 0xf0)>>4];
  *s=c[0] ;  *(s+1)=c[1] ;  *(s+2)=c[2] ;  *(s+3)=c[3] ;
  c=nybble[i & 0x0f];
  *(s+4)=c[0] ;  *(s+5)=c[1] ;  *(s+6)=c[2] ;  *(s+7)=c[3] ;
}
static void
make_float_bigendian (float * x)
{
  union {
    float f;
    unsigned char b[4];
  } u,v;

  u.f = *x ;

  v.b[0]=u.b[3] ;
  v.b[1]=u.b[2] ;
  v.b[2]=u.b[1] ;
  v.b[3]=u.b[0] ;

  *x=v.f ;
}

static void
make_double_bigendian (double * x)
{
  union {
    double d;
    unsigned char b[8];
  } u,v;

  u.d = *x ;

  v.b[0]=u.b[7] ;
  v.b[1]=u.b[6] ;
  v.b[2]=u.b[5] ;
  v.b[3]=u.b[4] ;
  v.b[4]=u.b[3] ;
  v.b[5]=u.b[2] ;
  v.b[6]=u.b[1] ;
  v.b[7]=u.b[0] ;

  *x=v.d ;
}

void ieee_fprintf_float (FILE * stream, const float * x) {
  ieee_float_rep_s r ;
  ieee_float_to_rep(x, &r) ;

  switch (r.type)
    {
    case SMF_IEEE_TYPE_NAN:
      fprintf(stream, "NaN") ;
      break ;
    case SMF_IEEE_TYPE_INF:
      fprintf(stream, "%cInf", signs[r.sign]) ;
      break ;
    case SMF_IEEE_TYPE_NORMAL:
      fprintf(stream, "%c1.%s*2^%d", signs[r.sign], r.mantissa, r.exponent) ;
      break ;
    case SMF_IEEE_TYPE_DENORMAL:
      fprintf(stream, "%c0.%s*2^%d", signs[r.sign], r.mantissa, r.exponent + 1) ;
      break ;
    case SMF_IEEE_TYPE_ZERO:
      fprintf(stream, "%c0", signs[r.sign]) ;
      break ;
    default:
      fprintf(stream, "[non-standard IEEE float]") ;
    }
}

void ieee_printf_float (const float * x)
{
  ieee_fprintf_float (stdout,x);
}

void ieee_fprintf_double (FILE * stream, const double * x) {
  ieee_double_rep_s r ;
  ieee_double_to_rep (x, &r) ;

  switch (r.type)
    {
    case SMF_IEEE_TYPE_NAN:
      fprintf(stream, "NaN") ;
      break ;
    case SMF_IEEE_TYPE_INF:
      fprintf(stream, "%cInf", signs[r.sign]) ;
      break ;
    case SMF_IEEE_TYPE_NORMAL:
      fprintf(stream, "%c1.%s*2^%d", signs[r.sign], r.mantissa, r.exponent) ;
      break ;
    case SMF_IEEE_TYPE_DENORMAL:
      fprintf(stream, "%c0.%s*2^%d", signs[r.sign], r.mantissa, r.exponent + 1) ;
      break ;
    case SMF_IEEE_TYPE_ZERO:
      fprintf(stream, "%c0", signs[r.sign]) ;
      break ;
    default:
      fprintf(stream, "[non-standard IEEE double]") ;
    }
}

void ieee_printf_double (const double * x)
{
  ieee_fprintf_double (stdout,x);
}

static int little_endian_p (void) {
  /* Are we little or big endian?  From Harbison & Steele.  */
  union
  {
    long l;
    char c[sizeof (long)];
  } u;
  u.l = 1;
  return (u.c[sizeof (long) - 1] == 1);
}

static int determine_ieee_type (int non_zero, int exponent, int max_exponent)
{
  if (exponent == max_exponent)
    {
      if (non_zero)
        {
          return SMF_IEEE_TYPE_NAN ;
        }
      else
        {
          return SMF_IEEE_TYPE_INF ;
        }
    }
  else if (exponent == 0)
    {
      if (non_zero)
        {
          return SMF_IEEE_TYPE_DENORMAL ;
        }
      else
        {
          return SMF_IEEE_TYPE_ZERO ;
        }
    }
  else
    {
      return SMF_IEEE_TYPE_NORMAL ;
    }
}

void ieee_float_to_rep (const float * x, ieee_float_rep_s * r)
{
  int e, non_zero;

  union {
    float f;
    struct  {
      unsigned char byte[4] ;
    } ieee ;
  } u;

  u.f = *x ;

  if (little_endian_p())
    make_float_bigendian(&(u.f)) ;

  /* note that r->sign is signed, u.ieee.byte is unsigned */

  if (u.ieee.byte[3]>>7)
    {
      r->sign = 1 ;
    }
  else
    {
      r->sign = 0 ;
    }

  e = (u.ieee.byte[3] & 0x7f) << 1 | (u.ieee.byte[2] & 0x80)>>7 ;

  r->exponent = e - 127 ;

  sprint_byte((u.ieee.byte[2] & 0x7f) << 1,r->mantissa) ;
  sprint_byte(u.ieee.byte[1],r->mantissa + 7) ;
  sprint_byte(u.ieee.byte[0],r->mantissa + 15) ;

  r->mantissa[23] = '\0' ;

  non_zero = u.ieee.byte[0] || u.ieee.byte[1] || (u.ieee.byte[2] & 0x7f);

  r->type = determine_ieee_type (non_zero, e, 255) ;
}

void ieee_double_to_rep (const double * x, ieee_double_rep_s * r)
{

  int e, non_zero;

  union
  {
    double d;
    struct  {
      unsigned char byte[8];
    } ieee ;
  } u;

  u.d= *x ;

  if (little_endian_p())
    make_double_bigendian(&(u.d)) ;

  /* note that r->sign is signed, u.ieee.byte is unsigned */

  if (u.ieee.byte[7]>>7)
    {
      r->sign = 1 ;
    }
  else
    {
      r->sign = 0 ;
    }


  e =(u.ieee.byte[7] & 0x7f)<<4 ^ (u.ieee.byte[6] & 0xf0)>>4 ;

  r->exponent = e - 1023 ;

  sprint_nybble(u.ieee.byte[6],r->mantissa) ;
  sprint_byte(u.ieee.byte[5],r->mantissa + 4) ;
  sprint_byte(u.ieee.byte[4],r->mantissa + 12) ;
  sprint_byte(u.ieee.byte[3],r->mantissa + 20) ;
  sprint_byte(u.ieee.byte[2],r->mantissa + 28) ;
  sprint_byte(u.ieee.byte[1],r->mantissa + 36) ;
  sprint_byte(u.ieee.byte[0],r->mantissa + 44) ;

  r->mantissa[52] = '\0' ;

  non_zero = (u.ieee.byte[0] || u.ieee.byte[1] || u.ieee.byte[2]
              || u.ieee.byte[3] || u.ieee.byte[4] || u.ieee.byte[5]
              || (u.ieee.byte[6] & 0x0f)) ;

  r->type = determine_ieee_type (non_zero, e, 2047) ;
}


} //END MATH
} //end SMF
