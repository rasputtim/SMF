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

#include "math/SMF_Odife.h"
#include "math/SMF_Vector.h"
#include "math/SMF_Matriz.h"
#include "math/SMF_EulerAngles.h"
#include "util/SMF_StringUtils.h"

//#pragma hdrstop

namespace SMF{
namespace MATH{

static void central_deriv (const CNumDEqu * f, float x, float h,
               float *result, float *abserr_round, float *abserr_trunc)
{
  /* Compute the derivative using the 5-point rule (x-h, x-h/2, x,
     x+h/2, x+h). Note that the central point is not used.  

     Compute the error using the difference between the 5-point and
     the 3-point rule (x-h,x,x+h). Again the central point is not
     used. */

  float fm1 = SMF_FN_EVAL (f, x - h);
  float fp1 = SMF_FN_EVAL (f, x + h);

  float fmh = SMF_FN_EVAL (f, x - h / 2);
  float fph = SMF_FN_EVAL (f, x + h / 2);

  float r3 = 0.5 * (fp1 - fm1);
  float r5 = (4.0 / 3.0) * (fph - fmh) - (1.0 / 3.0) * r3;

  float e3 = (CMath::fabs (fp1) + CMath::fabs (fm1)) * CMath::FLT_EPSILON;
  float e5 = 2.0 * (CMath::fabs (fph) + CMath::fabs (fmh)) * CMath::FLT_EPSILON + e3;

  float dy = MAX (CMath::fabs (r3), CMath::fabs (r5)) * CMath::fabs (x) * CMath::FLT_EPSILON;

  /* The truncation error in the r5 approximation itself is O(h^4).
     However, for safety, we estimate the error from r5-r3, which is
     O(h^2).  By scaling h we will minimise this estimated error, not
     the actual truncation error in r5. */

  *result = r5 / h;
  *abserr_trunc = CMath::fabs ((r5 - r3) / h); /* Estimated truncation error O(h^2) */
  *abserr_round = CMath::fabs (e5 / h) + dy;   /* Rounding error (cancellations) */
}

int CNumDEqu::deriv_central ( float x, float h,
                   float *result, float *abserr)
{
  //const derive_function_cb * f=derive;
  float r_0, round, trunc, error;
  central_deriv (this, x, h, &r_0, &round, &trunc);
  error = round + trunc;

  if (round < trunc && (round > 0 && trunc > 0))
    {
      float r_opt, round_opt, trunc_opt, error_opt;

      /* Compute an optimised stepsize to minimize the total error,
         using the scaling of the truncation error (O(h^2)) and
         rounding error (O(1/h)). */

      float h_opt = h * CMath::pow (round / (2.0 * trunc), 1.0 / 3.0);
      central_deriv (this, x, h_opt, &r_opt, &round_opt, &trunc_opt);
      error_opt = round_opt + trunc_opt;

      /* Check that the new error is smaller, and that the new derivative 
         is consistent with the error bounds of the original estimate. */

      if (error_opt < error && CMath::fabs (r_opt - r_0) < 4.0 * error)
        {
          r_0 = r_opt;
          error = error_opt;
        }
    }

  *result = r_0;
  *abserr = error;

  return SMF_OK;
}


static void forward_deriv (const CNumDEqu * f, float x, float h,
               float *result, float *abserr_round, float *abserr_trunc)
{
  /* Compute the derivative using the 4-point rule (x+h/4, x+h/2,
     x+3h/4, x+h).

     Compute the error using the difference between the 4-point and
     the 2-point rule (x+h/2,x+h).  */

  float f1 = SMF_FN_EVAL (f, x + h / 4.0);
  float f2 = SMF_FN_EVAL (f, x + h / 2.0);
  float f3 = SMF_FN_EVAL (f, x + (3.0 / 4.0) * h);
  float f4 = SMF_FN_EVAL (f, x + h);

  float r2 = 2.0*(f4 - f2);
  float r4 = (22.0 / 3.0) * (f4 - f3) - (62.0 / 3.0) * (f3 - f2) +
    (52.0 / 3.0) * (f2 - f1);

  /* Estimate the rounding error for r4 */

  float e4 = 2 * 20.67 * (CMath::fabs (f4) + CMath::fabs (f3) + CMath::fabs (f2) + CMath::fabs (f1)) * CMath::FLT_EPSILON;

  float dy = MAX (CMath::fabs (r2), CMath::fabs (r4)) * CMath::fabs (x) * CMath::FLT_EPSILON;

  /* The truncation error in the r4 approximation itself is O(h^3).
     However, for safety, we estimate the error from r4-r2, which is
     O(h).  By scaling h we will minimise this estimated error, not
     the actual truncation error in r4. */

  *result = r4 / h;
  *abserr_trunc = CMath::fabs ((r4 - r2) / h); /* Estimated truncation error O(h) */
  *abserr_round = CMath::fabs (e4 / h) + dy;
}

int CNumDEqu::deriv_forward ( float x, float h,
                   float *result, float *abserr)
{
	
  float r_0, round, trunc, error;
  forward_deriv (this, x, h, &r_0, &round, &trunc);
  error = round + trunc;

  if (round < trunc && (round > 0 && trunc > 0))
    {
      float r_opt, round_opt, trunc_opt, error_opt;

      /* Compute an optimised stepsize to minimize the total error,
         using the scaling of the estimated truncation error (O(h)) and
         rounding error (O(1/h)). */

      float h_opt = h * CMath::pow (round / (trunc), 1.0 / 2.0);
      forward_deriv (this, x, h_opt, &r_opt, &round_opt, &trunc_opt);
      error_opt = round_opt + trunc_opt;

      /* Check that the new error is smaller, and that the new derivative 
         is consistent with the error bounds of the original estimate. */

      if (error_opt < error && CMath::fabs (r_opt - r_0) < 4.0 * error)
        {
          r_0 = r_opt;
          error = error_opt;
        }
    }

  *result = r_0;
  *abserr = error;

  return SMF_OK;
}

int CNumDEqu::deriv_backward ( float x, float h, float *result, float *abserr)
{
	
  return deriv_forward ( x, -h, result, abserr);
}



float f (float x, void * params)
{
  return CMath::pow (x, 1.5);
}

int CNumDEqu::testDerivative(void)
{
  //derive_function_cb F;
  float result, abserr;
  CNumDEqu F;
  F.function = &f;
  F.params = 0;

  Debug::debug(Debug::math,__FUNCTION__) <<"f(x) = x^(3/2)"<<endl;
  //CNumDEqu func(&F);
  F.deriv_central ( 2.0, 1e-8, &result, &abserr);
  Debug::debug(Debug::math,__FUNCTION__) <<"x = 2.0\n"<<endl;
  Debug::debug(Debug::math,__FUNCTION__) <<"f'(x) = "<< result << " +/- " <<  abserr <<endl;
  Debug::debug(Debug::math,__FUNCTION__) <<"exact = " <<  1.5 * sqrt(2.0) << endl;

  F.deriv_forward ( 0.0, 1e-8, &result, &abserr);
  Debug::debug(Debug::math,__FUNCTION__) <<"x = 0.0"<<endl;
  Debug::debug(Debug::math,__FUNCTION__) <<"f'(x) = "<< result << " +/- "<< abserr <<endl;
  Debug::debug(Debug::math,__FUNCTION__) <<"exact = "<< 0.0 <<endl;

  return 0;
}


//===============================================================
//
//	CODEqu_Euler
//
//===============================================================

/*
=============
CODEqu_Euler::CODEqu_Euler
=============
*/
CODEqu_Euler::CODEqu_Euler( const int dim, deriveFunction_t dr, const void *ud ) {
	dimension = dim;
	derivatives = new float[dim];
	derive = dr;
	userData = ud;
}

/*
=============
CODEqu_Euler::~CODEqu_Euler
=============
*/
CODEqu_Euler::~CODEqu_Euler() {
	delete[] derivatives;
}

/*
=============
CODEqu_Euler::evaluate
=============
*/
float CODEqu_Euler::evaluate( const float *state, float *newState, float t0, float t1 ) {
	float delta;
	int i;

	derive( t0, userData, state, derivatives );
	delta = t1 - t0;
	for ( i = 0; i < dimension; i++ ) {
		newState[i] = state[i] + delta * derivatives[i];
	}
	return delta;
}

//===============================================================
//
//	CODEqu_Midpoint
//
//===============================================================

/*
=============
CODEqu_Midpoint::CODEqu_Midpoint
=============
*/
CODEqu_Midpoint::CODEqu_Midpoint( const int dim, deriveFunction_t dr, const void *ud ) {
	dimension = dim;
	tmpState = new float[dim];
	derivatives = new float[dim];
	derive = dr;
	userData = ud;
}

/*
=============
CODEqu_Midpoint::~CODEqu_Midpoint
=============
*/
CODEqu_Midpoint::~CODEqu_Midpoint() {
	delete tmpState;
	delete derivatives;
}

/*
=============
CODEqu_Midpoint::~evaluate
=============
*/
float CODEqu_Midpoint::evaluate( const float *state, float *newState, float t0, float t1 ) {
	float delta, halfDelta;
    int i;

	delta = t1 - t0;
	halfDelta = delta * 0.5;
    // first step
	derive( t0, userData, state, derivatives );
	for ( i = 0; i < dimension; i++ ) {
		tmpState[i] = state[i] + halfDelta * derivatives[i];
	}
    // second step
	derive( t0 + halfDelta, userData, tmpState, derivatives );

	for ( i = 0; i < dimension; i++ ) {
		newState[i] = state[i] + delta * derivatives[i];
	}
	return delta;
}

//===============================================================
//
//	CODEqu_RK4
//
//===============================================================

/*
=============
CODEqu_RK4::CODEqu_RK4
=============
*/
CODEqu_RK4::CODEqu_RK4( const int dim, deriveFunction_t dr, const void *ud ) {
	dimension = dim;
	derive = dr;
	userData = ud;
	tmpState = new float[dim];
	d1 = new float[dim];
	d2 = new float[dim];
	d3 = new float[dim];
	d4 = new float[dim];
}

/*
=============
CODEqu_RK4::~CODEqu_RK4
=============
*/
CODEqu_RK4::~CODEqu_RK4() {
	delete tmpState;
	delete d1;
	delete d2;
	delete d3;
	delete d4;
}

/*
=============
CODEqu_RK4::evaluate
=============
*/
float CODEqu_RK4::evaluate( const float *state, float *newState, float t0, float t1 ) {
	float delta, halfDelta, sixthDelta;
	int i;

	delta = t1 - t0;
	halfDelta = delta * 0.5;
	// first step
	derive( t0, userData, state, d1 );
	for ( i = 0; i < dimension; i++ ) {
		tmpState[i] = state[i] + halfDelta * d1[i];
	}
	// second step
	derive( t0 + halfDelta, userData, tmpState, d2 );
	for ( i = 0; i < dimension; i++ ) {
		tmpState[i] = state[i] + halfDelta * d2[i];
	}
	// third step
	derive( t0 + halfDelta, userData, tmpState, d3 );
	for ( i = 0; i < dimension; i++ ) {
		tmpState[i] = state[i] + delta * d3[i];
	}
	// fourth step
	derive( t0 + delta, userData, tmpState, d4 );

	sixthDelta = delta * (1.0/6.0);
	for ( i = 0; i < dimension; i++ ) {
		newState[i] = state[i] + sixthDelta * (d1[i] + 2.0 * (d2[i] + d3[i]) + d4[i]);
	}
	return delta;
}

//===============================================================
//
//	CODEqu_RK4Adaptive
//
//===============================================================

/*
=============
CODEqu_RK4Adaptive::CODEqu_RK4Adaptive
=============
*/
CODEqu_RK4Adaptive::CODEqu_RK4Adaptive( const int dim, deriveFunction_t dr, const void *ud ) {
	dimension = dim;
	derive = dr;
	userData = ud;
	maxError = 0.01f;
	tmpState = new float[dim];
	d1 = new float[dim];
	d1half = new float [dim];
	d2 = new float[dim];
	d3 = new float[dim];
	d4 = new float[dim];
}

/*
=============
CODEqu_RK4Adaptive::~CODEqu_RK4Adaptive
=============
*/
CODEqu_RK4Adaptive::~CODEqu_RK4Adaptive() {
	delete tmpState;
	delete d1;
	delete d1half;
	delete d2;
	delete d3;
	delete d4;
}

/*
=============
CODEqu_RK4Adaptive::setMaxError
=============
*/
void CODEqu_RK4Adaptive::setMaxError( const float err ) {
	if ( err > 0.0f ) {
		maxError = err;
	}
}

/*
=============
CODEqu_RK4Adaptive::evaluate
=============
*/
float CODEqu_RK4Adaptive::evaluate( const float *state, float *newState, float t0, float t1 ) {
	float delta, halfDelta, fourthDelta, sixthDelta;
	float error, max;
	int i, n;

	delta = t1 - t0;

	for ( n = 0; n < 4; n++ ) {

		halfDelta = delta * 0.5;
		fourthDelta = delta * 0.25;

		// first step of first half delta
		derive( t0, userData, state, d1 );
		for ( i = 0; i < dimension; i++ ) {
			tmpState[i] = state[i] + fourthDelta * d1[i];
		}
		// second step of first half delta
		derive( t0 + fourthDelta, userData, tmpState, d2 );
		for ( i = 0; i < dimension; i++ ) {
			tmpState[i] = state[i] + fourthDelta * d2[i];
		}
		// third step of first half delta
		derive( t0 + fourthDelta, userData, tmpState, d3 );
		for ( i = 0; i < dimension; i++ ) {
			tmpState[i] = state[i] + halfDelta * d3[i];
		}
		// fourth step of first half delta
		derive( t0 + halfDelta, userData, tmpState, d4 );

		sixthDelta = halfDelta * (1.0/6.0);
		for ( i = 0; i < dimension; i++ ) {
			tmpState[i] = state[i] + sixthDelta * (d1[i] + 2.0 * (d2[i] + d3[i]) + d4[i]);
		}

		// first step of second half delta
		derive( t0 + halfDelta, userData, tmpState, d1half );
		for ( i = 0; i < dimension; i++ ) {
			tmpState[i] = state[i] + fourthDelta * d1half[i];
		}
		// second step of second half delta
		derive( t0 + halfDelta + fourthDelta, userData, tmpState, d2 );
		for ( i = 0; i < dimension; i++ ) {
			tmpState[i] = state[i] + fourthDelta * d2[i];
		}
		// third step of second half delta
		derive( t0 + halfDelta + fourthDelta, userData, tmpState, d3 );
		for ( i = 0; i < dimension; i++ ) {
			tmpState[i] = state[i] + halfDelta * d3[i];
		}
		// fourth step of second half delta
		derive( t0 + delta, userData, tmpState, d4 );

		sixthDelta = halfDelta * (1.0/6.0);
		for ( i = 0; i < dimension; i++ ) {
			newState[i] = state[i] + sixthDelta * (d1[i] + 2.0 * (d2[i] + d3[i]) + d4[i]);
		}

		// first step of full delta
		for ( i = 0; i < dimension; i++ ) {
			tmpState[i] = state[i] + halfDelta * d1[i];
		}
		// second step of full delta
		derive( t0 + halfDelta, userData, tmpState, d2 );
		for ( i = 0; i < dimension; i++ ) {
			tmpState[i] = state[i] + halfDelta * d2[i];
		}
		// third step of full delta
		derive( t0 + halfDelta, userData, tmpState, d3 );
		for ( i = 0; i < dimension; i++ ) {
			tmpState[i] = state[i] + delta * d3[i];
		}
		// fourth step of full delta
		derive( t0 + delta, userData, tmpState, d4 );

		sixthDelta = delta * (1.0/6.0);
		for ( i = 0; i < dimension; i++ ) {
			tmpState[i] = state[i] + sixthDelta * (d1[i] + 2.0 * (d2[i] + d3[i]) + d4[i]);
		}

		// get max estimated error
        max = 0.0;
		for ( i = 0; i < dimension; i++ ) {
			error = CMath::fabs( (newState[i] - tmpState[i]) / (delta * d1[i] + 1e-10) );
			if ( error > max ) {
				max = error;
			}
        }
		error = max / maxError;

        if ( error <= 1.0f ) {
			return delta * 4.0;
		}
		if ( delta <= 1e-7 ) {
			return delta;
		}
		delta *= 0.25;
	}
	return delta;
}
} //end MATH
} //end SMF
