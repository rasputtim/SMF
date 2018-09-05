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

#include "math/SMF_Complex.h"
#include "util/SMF_StringUtils.h"

//#pragma hdrstop
namespace SMF {
namespace MATH{
CComplex complex_origin( 0.0f, 0.0f );


const CComplex CComplex::nan =CComplex(CMath::NAN_FLOAT, CMath::NAN_FLOAT);
/*
=============
CComplex::toString
=============
*/
const char *CComplex::toString( int precision ) const {
	return Util::CMyString::floatArraytoString( toFloatPtr(), getDimension(), precision );
}
/**********************************************************************
 * Elementary complex functions
 **********************************************************************/



#if 0
float _abs (CComplex z)
{                               /* return |z| */
  return hypot (z.r, z.i);
}

/* z=1/a */
CComplex _inverse (CComplex a)
{                               
  float s = 1.0 / (a).abs();

  CComplex z;
  z.set( (a.r * s) * s, -(a.i * s) * s);
  return z;
}

CComplex _mul_real (CComplex a, float x)
{                               /* z=a*x */
  CComplex z;
  z.set( x * a.r, x * a.i);
  return z;
}
#endif
CComplex _mul_imag (CComplex a, float y)
{                               /* z=a*iy */
  CComplex z;
  z.set( -y * a.i, y * a.r);
  return z;
}

CComplex _div (CComplex a, CComplex b)
{                               /* z=a/b */
  float ar = a.r, ai = a.i;
  float br = b.r, bi = b.i;

  float s = 1.0 / (b).abs();

  float sbr = s * br;
  float sbi = s * bi;

  float zr = (ar * sbr + ai * sbi) * s;
  float zi = (ai * sbr - ar * sbi) * s;

  CComplex z;
  z.set( zr, zi);
  return z;
}
/**********************************************************************
 * Properties of complex numbers
 **********************************************************************/

float _arg (CComplex z)
{                               /* return arg(z),  -pi < arg(z) <= +pi */
  float x = z.r;
  float y = z.i;

  if (x == 0.0 && y == 0.0)
    {
      return 0;
    }

  return CMath::atan(y, x);
}

float _logabs (CComplex z)
{                               /* return CMath::log|z| */
  float xabs = CMath::fabs (z.r);
  float yabs = CMath::fabs (z.i);
  float max, u;

  if (xabs >= yabs)
    {
      max = xabs;
      u = yabs / xabs;
    }
  else
    {
      max = yabs;
      u = xabs / yabs;
    }

  /* Handle underflow when u is close to 0 */

  return CMath::log (max) + 0.5 * CMath::log1p (u * u);
}




CComplex CComplex::sqrt (CComplex a)
{                               /* z=CMath::sqrt(a) */
  CComplex z;

  if (a.r == 0.0 && a.i == 0.0)
    {
	      z.set(0, 0);
    }
  else
    {
      float x = CMath::fabs (a.r);
      float y = CMath::fabs (a.i);
      float w;

      if (x >= y)
        {
          float t = y / x;
          w = CMath::sqrt (x) * CMath::sqrt (0.5 * (1.0 + CMath::sqrt (1.0 + t * t)));
        }
      else
        {
          float t = x / y;
          w = CMath::sqrt (y) * CMath::sqrt (0.5 * (t + CMath::sqrt (1.0 + t * t)));
        }

      if (a.r >= 0.0)
        {
          float ai = a.i;
          z.set( w, ai / (2.0 * w));
        }
      else
        {
          float ai = a.i;
          float vi = (ai >= 0) ? w : -w;
          z.set( ai / (2.0 * vi), vi);
        }
    }

  return z;
}

CComplex CComplex::sqrt_real (float x)
{                               /* z=CMath::sqrt(x) */
  CComplex z;

  if (x >= 0)
    {
      z.set( CMath::sqrt (x), 0.0);
    }
  else
    {
      z.set( 0.0, CMath::sqrt (-x));
    }

  return z;
}

CComplex CComplex::exp (CComplex a)
{                               /* z=CMath::exp(a) */
  float rho = CMath::exp (a.r);
  float theta = a.i;

  CComplex z;
  z.set( rho * CMath::cos (theta), rho * CMath::sin (theta));
  return z;
}

CComplex CComplex::pow (CComplex a, CComplex b)
{                               /* z=a^b */
  CComplex z;

  if (a.r == 0 && a.i == 0.0)
    {
      z.set( 0.0, 0.0);
    }
  else
    {
      float logr = _logabs (a);
      float theta = _arg (a);

      float br = b.r, bi = b.i;

      float rho = CMath::exp (logr * br - bi * theta);
      float beta = theta * br + bi * logr;

      z.set( rho * CMath::cos (beta), rho * CMath::sin (beta));
    }

  return z;
}

CComplex CComplex::pow_real (CComplex a, float b)
{                               /* z=a^b */
  CComplex z;

  if (a.r == 0 && a.i == 0)
    {
      z.set( 0, 0);
    }
  else
    {
      float logr = _logabs (a);
      float theta = _arg (a);
      float rho = CMath::exp (logr * b);
      float beta = theta * b;
      z.set( rho * CMath::cos (beta), rho * CMath::sin (beta));
    }

  return z;
}

CComplex CComplex::log (CComplex a)
{                               /* z=CMath::log(a) */
  float logr = _logabs (a);
  float theta = _arg (a);

  CComplex z;
  z.set( logr, theta);
  return z;
}

CComplex CComplex::log10 (CComplex a)
{                               /* z = log10(a) */
  return CComplex::log (a) * (1 / CMath::log (10.));
}

CComplex CComplex::log_b (CComplex a, CComplex b)
{
  return _div (CComplex::log (a), CComplex::log (b));
}

/***********************************************************************
 * Complex trigonometric functions
 ***********************************************************************/

CComplex CComplex::sin (CComplex a)
{                               /* z = CMath::sin(a) */
  float RR = a.r, I = a.i;

  CComplex z;

  if (I == 0.0) 
    {
      /* avoid returing negative zero (-0.0) for the imaginary part  */

      z.set( CMath::sin (RR), 0.0);  
    } 
  else 
    {
      z.set( CMath::sin(RR) * std::cosh (I), CMath::cos(RR) * std::sinh (I));
    }

  return z;
}

CComplex CComplex::cos (CComplex a)
{                               /* z = CMath::cos(a) */
  float R = a.r, I = a.i;

  CComplex z;

  if (I == 0.0) 
    {
      /* avoid returing negative zero (-0.0) for the imaginary part  */

      z.set( CMath::cos (R), 0.0);  
    } 
  else 
    {
      z.set( CMath::cos (R) * std::cosh (I), CMath::sin (R) * std::sinh (-I));
    }

  return z;
}

CComplex CComplex::tan (CComplex a)
{                               /* z = tan(a) */
  float R = a.r, I = a.i;

  CComplex z;

  if (CMath::fabs (I) < 1)
    {
      float D = CMath::pow (CMath::cos (R), 2.0) + CMath::pow (std::sinh (I), 2.0);

      z.set( 0.5 * CMath::sin (2 * R) / D, 0.5 * std::sinh (2 * I) / D);
    }
  else
    {
      float u = CMath::exp (-I);
      float C = 2 * u / (1 - CMath::pow (u, 2.0));
      float D = 1 + CMath::pow (CMath::cos (R), 2.0) * CMath::pow (C, 2.0);

      float S = CMath::pow (C, 2.0);
      float T = 1.0 / std::tanh (I);

      z.set( 0.5 * CMath::sin (2 * R) * S / D, T / D);
    }

  return z;
}

CComplex CComplex::sec (CComplex a)
{                               /* z = sec(a) */
  CComplex z = CComplex::cos (a);
  return (z).inverse();
}

CComplex CComplex::csc (CComplex a)
{                               /* z = csc(a) */
  CComplex z = CComplex::sin (a);
  return (z).inverse();
}


CComplex CComplex::cot (CComplex a)
{                               /* z = cot(a) */
  CComplex z = CComplex::tan (a);
  return (z).inverse();
}

/**********************************************************************
 * Inverse Complex Trigonometric Functions
 **********************************************************************/

CComplex CComplex::arcsin (CComplex a)
{                               /* z = arcsin(a) */
  float R = a.r, I = a.i;
  CComplex z;

  if (I == 0)
    {
      z = CComplex::arcsin_real (R);
    }
  else
    {
      float x = CMath::fabs (R), y = CMath::fabs (I);
      float r = hypot (x + 1, y), s = hypot (x - 1, y);
      float A = 0.5 * (r + s);
      float B = x / A;
      float y2 = y * y;

      float real, imag;

      const float A_crossover = 1.5, B_crossover = 0.6417;

      if (B <= B_crossover)
        {
          real = CMath::asin (B);
        }
      else
        {
          if (x <= 1)
            {
              float D = 0.5 * (A + x) * (y2 / (r + x + 1) + (s + (1 - x)));
              real = CMath::atan (x / CMath::sqrt (D));
            }
          else
            {
              float Apx = A + x;
              float D = 0.5 * (Apx / (r + x + 1) + Apx / (s + (x - 1)));
              real = CMath::atan (x / (y * CMath::sqrt (D)));
            }
        }

      if (A <= A_crossover)
        {
          float Am1;

          if (x < 1)
            {
              Am1 = 0.5 * (y2 / (r + (x + 1)) + y2 / (s + (1 - x)));
            }
          else
            {
              Am1 = 0.5 * (y2 / (r + (x + 1)) + (s + (x - 1)));
            }

          imag = CMath::log1p (Am1 + CMath::sqrt (Am1 * (A + 1)));
        }
      else
        {
          imag = CMath::log (A + CMath::sqrt (A * A - 1));
        }

      z.set( (R >= 0) ? real : -real, (I >= 0) ? imag : -imag);
    }

  return z;
}

CComplex CComplex::arcsin_real (float a)
{                               /* z = arcsin(a) */
  CComplex z;

  if (CMath::fabs (a) <= 1.0)
    {
      z.set( CMath::asin (a), 0.0);
    }
  else
    {
      if (a < 0.0)
        {
          z.set( -M_PI_2, CMath::acosh (-a));
        }
      else
        {
          z.set( M_PI_2, -CMath::acosh (a));
        }
    }

  return z;
}

CComplex CComplex::arccos (CComplex a)
{                               /* z = arccos(a) */
  float R = a.r, I = a.i;
  CComplex z;

  if (I == 0)
    {
      z = CComplex::arccos_real (R);
    }
  else
    {
      float x = CMath::fabs (R), y = CMath::fabs (I);
      float r = hypot (x + 1, y), s = hypot (x - 1, y);
      float A = 0.5 * (r + s);
      float B = x / A;
      float y2 = y * y;

      float real, imag;

      const float A_crossover = 1.5, B_crossover = 0.6417;

      if (B <= B_crossover)
        {
          real = CMath::acos (B);
        }
      else
        {
          if (x <= 1)
            {
              float D = 0.5 * (A + x) * (y2 / (r + x + 1) + (s + (1 - x)));
              real = CMath::atan (CMath::sqrt (D) / x);
            }
          else
            {
              float Apx = A + x;
              float D = 0.5 * (Apx / (r + x + 1) + Apx / (s + (x - 1)));
              real = CMath::atan ((y * CMath::sqrt (D)) / x);
            }
        }

      if (A <= A_crossover)
        {
          float Am1;

          if (x < 1)
            {
              Am1 = 0.5 * (y2 / (r + (x + 1)) + y2 / (s + (1 - x)));
            }
          else
            {
              Am1 = 0.5 * (y2 / (r + (x + 1)) + (s + (x - 1)));
            }

          imag = CMath::log1p (Am1 + CMath::sqrt (Am1 * (A + 1)));
        }
      else
        {
          imag = CMath::log (A + CMath::sqrt (A * A - 1));
        }

      z.set( (R >= 0) ? real : CMath::PI - real, (I >= 0) ? -imag : imag);
    }

  return z;
}

CComplex CComplex::arccos_real (float a)
{                               /* z = arccos(a) */
  CComplex z;

  if (CMath::fabs (a) <= 1.0)
    {
      z.set( CMath::acos (a), 0);
    }
  else
    {
      if (a < 0.0)
        {
          z.set( CMath::PI, -CMath::acosh (-a));
        }
      else
        {
          z.set( 0, CMath::acosh (a));
        }
    }

  return z;
}

CComplex CComplex::arctan (CComplex a)
{                               /* z = arctan(a) */
  float R = a.r, I = a.i;
  CComplex z;

  if (I == 0)
    {
      z.set( CMath::atan (R), 0);
    }
  else
    {
      /* FIXME: This is a naive implementation which does not fully
         take into account cancellation errors, overflow, underflow
         etc.  It would benefit from the Hull et al treatment. */

      float r = hypot (R, I);

      float imag;

      float u = 2 * I / (1 + r * r);

      /* FIXME: the following cross-over should be optimized but 0.1
         seems to work ok */

      if (CMath::fabs (u) < 0.1)
        {
          imag = 0.25 * (CMath::log1p (u) - CMath::log1p (-u));
        }
      else
        {
          float A = hypot (R, I + 1);
          float B = hypot (R, I - 1);
          imag = 0.5 * CMath::log (A / B);
        }

      if (R == 0)
        {
          if (I > 1)
            {
              z.set( M_PI_2, imag);
            }
          else if (I < -1)
            {
              z.set( -M_PI_2, imag);
            }
          else
            {
              z.set( 0, imag);
            };
        }
      else
        {
			z.set( 0.5 * CMath::atan(2 * R, ((1 + r) * (1 - r))), imag);
        }
    }

  return z;
}

CComplex CComplex::arcsec (CComplex a)
{                               /* z = arcsec(a) */
  CComplex z = (a).inverse();
  return CComplex::arccos (z);
}

CComplex CComplex::arcsec_real (float a)
{                               /* z = arcsec(a) */
  CComplex z;

  if (a <= -1.0 || a >= 1.0)
    {
      z.set( CMath::acos (1 / a), 0.0);
    }
  else
    {
      if (a >= 0.0)
        {
          z.set( 0, CMath::acosh (1 / a));
        }
      else
        {
          z.set( CMath::PI, -CMath::acosh (-1 / a));
        }
    }

  return z;
}

CComplex CComplex::arccsc (CComplex a)
{                               /* z = arccsc(a) */
  CComplex z = (a).inverse();
  return CComplex::arcsin (z);
}

CComplex CComplex::arccsc_real (float a)
{                               /* z = arccsc(a) */
  CComplex z;

  if (a <= -1.0 || a >= 1.0)
    {
      z.set( CMath::asin (1 / a), 0.0);
    }
  else
    {
      if (a >= 0.0)
        {
          z.set( M_PI_2, -CMath::acosh (1 / a));
        }
      else
        {
          z.set( -M_PI_2, CMath::acosh (-1 / a));
        }
    }

  return z;
}

CComplex CComplex::arccot (CComplex a)
{                               /* z = arccot(a) */
  CComplex z;

  if (a.r == 0.0 && a.i == 0.0)
    {
      z.set( M_PI_2, 0);
    }
  else
    {
      z = (a).inverse();
      z = CComplex::arctan (z);
    }

  return z;
}

/**********************************************************************
 * Complex Hyperbolic Functions
 **********************************************************************/

CComplex CComplex::sinh (CComplex a)
{                               /* z = sinh(a) */
  float R = a.r, I = a.i;

  CComplex z;
  z.set( std::sinh (R) * CMath::cos (I), std::cosh (R) * CMath::sin (I));
  return z;
}

CComplex CComplex::cosh (CComplex a)
{                               /* z = cosh(a) */
  float R = a.r, I = a.i;

  CComplex z;
  z.set( std::cosh (R) * CMath::cos (I), std::sinh (R) * CMath::sin (I));
  return z;
}

CComplex CComplex::tanh (CComplex a)
{                               /* z = tanh(a) */
  float R = a.r, I = a.i;

  CComplex z;

  if (CMath::fabs(R) < 1.0) 
    {
      float D = CMath::pow (CMath::cos (I), 2.0) + CMath::pow (std::sinh (R), 2.0);
      
      z.set( std::sinh (R) * std::cosh (R) / D, 0.5 * CMath::sin (2 * I) / D);
    }
  else
    {
      float D = CMath::pow (CMath::cos (I), 2.0) + CMath::pow (std::sinh (R), 2.0);
      float F = 1 + CMath::pow (CMath::cos (I) / std::sinh (R), 2.0);

      z.set( 1.0 / (std::tanh (R) * F), 0.5 * CMath::sin (2 * I) / D);
    }

  return z;
}

CComplex CComplex::sech (CComplex a)
{                               /* z = sech(a) */
  CComplex z = CComplex::cosh (a);
  return (z).inverse();
}

CComplex CComplex::csch (CComplex a)
{                               /* z = csch(a) */
  CComplex z = CComplex::sinh (a);
  return (z).inverse();
}

CComplex CComplex::coth (CComplex a)
{                               /* z = coth(a) */
  CComplex z = CComplex::tanh (a);
  return (z).inverse();
}

/**********************************************************************
 * Inverse Complex Hyperbolic Functions
 **********************************************************************/

CComplex CComplex::arcsinh (CComplex a)
{                               /* z = arcsinh(a) */
  CComplex z = _mul_imag(a, 1.0);
  z = CComplex::arcsin (z);
  z = _mul_imag (z, -1.0);
  return z;
}

CComplex CComplex::arccosh (CComplex a)
{                               /* z = arccosh(a) */
  CComplex z = CComplex::arccos (a);
  z = _mul_imag (z, z.i > 0 ? -1.0 : 1.0);
  return z;
}

CComplex CComplex::arccosh_real (float a)
{                               /* z = arccosh(a) */
  CComplex z;

  if (a >= 1)
    {
      z.set( CMath::acosh (a), 0);
    }
  else
    {
      if (a >= -1.0)
        {
          z.set( 0, CMath::acos (a));
        }
      else
        {
          z.set( CMath::acosh (-a), CMath::PI);
        }
    }

  return z;
}

CComplex CComplex::arctanh (CComplex a)
{                               /* z = arctanh(a) */
  if (a.i == 0.0)
    {
      return CComplex::arctanh_real (a.r);
    }
  else
    {
      CComplex z = _mul_imag(a, 1.0);
      z = CComplex::arctan (z);
      z = _mul_imag (z, -1.0);
      return z;
    }
}

CComplex CComplex::arctanh_real (float a)
{                               /* z = arctanh(a) */
  CComplex z;

  if (a > -1.0 && a < 1.0)
    {
      z.set( CMath::atanh (a), 0);
    }
  else
    {
      z.set( CMath::atanh (1 / a), (a < 0) ? M_PI_2 : -M_PI_2);
    }

  return z;
}

CComplex CComplex::arcsech (CComplex a)
{                               /* z = arcsech(a); */
  CComplex t = (a).inverse();
  return CComplex::arccosh (t);
}

CComplex CComplex::arccsch (CComplex a)
{                               /* z = arccsch(a) */
  CComplex t = (a).inverse();
  return CComplex::arcsinh (t);
}

CComplex CComplex::arccoth (CComplex a)
{                               /* z = arccoth(a) */
  CComplex t = (a).inverse();
  return CComplex::arctanh (t);
}

} //end MATH
} //end SMF
