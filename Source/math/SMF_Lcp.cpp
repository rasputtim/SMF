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

#include "math/SMF_Lincp.h"
#include "math/SMF_Vector.h"
#include "math/SMF_Matriz.h"
#include "structures/SMF_List.h"
#include "util/SMF_StringUtils.h"
#include "util/SMF_Debug.h"
//#pragma hdrstop
namespace SMF {
namespace MATH{
//static CSysVar lcp_showFailures( "lcp_showFailures", "0", SYSVAR_SYSTEM | SYSVAR_BOOL, "show LCP solver failures" );

const float LCP_BOUND_EPSILON			= 1e-5f;
const float LCP_ACCEL_EPSILON			= 1e-5f;
const float LCP_DELTA_ACCEL_EPSILON		= 1e-9f;
const float LCP_DELTA_FORCE_EPSILON		= 1e-9f;

//#define IGNORE_UNSATISFIABLE_VARIABLES

//===============================================================
//                                                        M
//  CLinCP_Square                                         MrE
//                                                        E
//===============================================================
/**
 * 
 *
 * \ingroup SMF_Math
 *
 * \if pt_br
 * \brief	Ajuda na solução do problema de \ref CLinCP \ref
 *
 * \elseif us_en
 * \brief 	Helper to solve LCP problem
 *
 * \endif
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class CLinCP_Square : public CLinCP {
public:
	virtual bool	solve( const CMatXD &o_m, CVecXD &o_x, const CVecXD &o_b, const CVecXD &o_lo, const CVecXD &o_hi, const int *o_boxIndex );

private:
	CMatXD			m;					// original matrix
	CVecXD			b;					// right hand side
	CVecXD			lo, hi;				// low and high bounds
	CVecXD			f, a;				// force and acceleration
	CVecXD			delta_f, delta_a;	// delta force and delta acceleration
	CMatXD			clamped;			// LU factored sub matrix for clamped variables
	CVecXD			diagonal;			// reciprocal of diagonal of U of the LU factored sub matrix for clamped variables
	int				numUnbounded;		// number of unbounded variables
	int				numClamped;			// number of clamped variables
	float **		rowPtrs;			// pointers to the rows of m
	int *			boxIndex;			// box index
	int *			side;				// tells if a variable is at the low boundary = -1, high boundary = 1 or inbetween = 0
	int *			permuted;			// index to keep track of the permutation
	bool			padded;				// set to true if the rows of the initial matrix are 16 byte padded

private:
	bool			FactorClamped();
	void			solveClamped( CVecXD &x, const float *b );
	void			Swap( int i, int j );
	void			AddClamped( int r );
	void			RemoveClamped( int r );
	void			CalcForceDelta( int d, float dir );
	void			CalcAccelDelta( int d );
	void			ChangeForce( int d, float step );
	void			ChangeAccel( int d, float step );
	void			GetMaxStep( int d, float dir, float &maxStep, int &limit, int &limitSide ) const;
};

/*
============
CLinCP_Square::FactorClamped
============
*/
bool CLinCP_Square::FactorClamped() {
	int i, j, k;
	float s, d;

	for ( i = 0; i < numClamped; i++ ) {
		memcpy( clamped[i], rowPtrs[i], numClamped * sizeof( float ) );
	}

	for ( i = 0; i < numClamped; i++ ) {

		s = CMath::fabs( clamped[i][i] );

		if ( s == 0.0f ) {
			return false;
		}

		diagonal[i] = d = 1.0f / clamped[i][i];
		for ( j = i + 1; j < numClamped; j++ ) {
			clamped[j][i] *= d;
		}

		for ( j = i + 1; j < numClamped; j++ ) {
			d = clamped[j][i];
			for ( k = i + 1; k < numClamped; k++ ) {
				clamped[j][k] -= d * clamped[i][k];
			}
		}
	}

	return true;
}

/*
============
CLinCP_Square::solveClamped
============
*/
void CLinCP_Square::solveClamped( CVecXD &x, const float *b ) {
	int i, j;
	float sum;

	// solve L
	for ( i = 0; i < numClamped; i++ ) {
		sum = b[i];
		for ( j = 0; j < i; j++ ) {
			sum -= clamped[i][j] * x[j];
		}
		x[i] = sum;
	}

	// solve U
	for ( i = numClamped - 1; i >= 0; i-- ) {
		sum = x[i];
		for ( j = i + 1; j < numClamped; j++ ) {
			sum -= clamped[i][j] * x[j];
		}
		x[i] = sum * diagonal[i];
	}
}

/*
============
CLinCP_Square::Swap
============
*/
void CLinCP_Square::Swap( int i, int j ) {

	if ( i == j ) {
		return;
	}

	swapElements( rowPtrs[i], rowPtrs[j] );
	m.swapColumns( i, j );
	b.swapElements( i, j );
	lo.swapElements( i, j );
	hi.swapElements( i, j );
	a.swapElements( i, j );
	f.swapElements( i, j );
	if ( boxIndex ) {
		swapElements( boxIndex[i], boxIndex[j] );
	}
	swapElements( side[i], side[j] );
	swapElements( permuted[i], permuted[j] );
}

/*
============
CLinCP_Square::AddClamped
============
*/
void CLinCP_Square::AddClamped( int r ) {
	int i, j;
	float sum;

	SMF_ASSERT( r >= numClamped );

	// add a row at the bottom and a column at the right of the factored
	// matrix for the clamped variables

	Swap( numClamped, r );

	// add row to L
	for ( i = 0; i < numClamped; i++ ) {
		sum = rowPtrs[numClamped][i];
		for ( j = 0; j < i; j++ ) {
			sum -= clamped[numClamped][j] * clamped[j][i];
		}
		clamped[numClamped][i] = sum * diagonal[i];
	}

	// add column to U
	for ( i = 0; i <= numClamped; i++ ) {
		sum = rowPtrs[i][numClamped];
		for ( j = 0; j < i; j++ ) {
			sum -= clamped[i][j] * clamped[j][numClamped];
		}
		clamped[i][numClamped] = sum;
	}

	diagonal[numClamped] = 1.0f / clamped[numClamped][numClamped];

	numClamped++;
}

/*
============
CLinCP_Square::RemoveClamped
============
*/
void CLinCP_Square::RemoveClamped( int r ) {
	int i, j;
	float *y0, *y1, *z0, *z1;
	double diag, beta0, beta1, p0, p1, q0, q1, d;

	SMF_ASSERT( r < numClamped );

	numClamped--;

	// no need to swap and update the factored matrix when the last row and column are removed
	if ( r == numClamped ) {
		return;
	}
        
	y0 = (float *) _allocafloat16( numClamped * sizeof( float ) );
	z0 = (float *) _allocafloat16( numClamped * sizeof( float ) );
	y1 = (float *) _allocafloat16( numClamped * sizeof( float ) );
	z1 = (float *) _allocafloat16( numClamped * sizeof( float ) );

	// the row/column need to be subtracted from the factorization
	for ( i = 0; i < numClamped; i++ ) {
		y0[i] = -rowPtrs[i][r];
	}

	memset( y1, 0, numClamped * sizeof( float ) );
	y1[r] = 1.0f;

	memset( z0, 0, numClamped * sizeof( float ) );
	z0[r] = 1.0f;

	for ( i = 0; i < numClamped; i++ ) {
		z1[i] = -rowPtrs[r][i];
	}

	// swap the to be removed row/column with the last row/column
	Swap( r, numClamped );

	// the swapped last row/column need to be added to the factorization
	for ( i = 0; i < numClamped; i++ ) {
		y0[i] += rowPtrs[i][r];
	}

	for ( i = 0; i < numClamped; i++ ) {
		z1[i] += rowPtrs[r][i];
	}
	z1[r] = 0.0f;

	// update the beginning of the to be updated row and column
	for ( i = 0; i < r; i++ ) {
		p0 = y0[i];
		beta1 = z1[i] * diagonal[i];

		clamped[i][r] += p0;
		for ( j = i+1; j < numClamped; j++ ) {
			z1[j] -= beta1 * clamped[i][j];
		}
		for ( j = i+1; j < numClamped; j++ ) {
			y0[j] -= p0 * clamped[j][i];
		}
		clamped[r][i] += beta1;
	}

	// update the lower right corner starting at r,r
	for ( i = r; i < numClamped; i++ ) {
		diag = clamped[i][i];

		p0 = y0[i];
		p1 = z0[i];
		diag += p0 * p1;

		if ( diag == 0.0f ) {
			Debug::debug(Debug::math,__FUNCTION__)<< "CLinCP_Square::RemoveClamped: updating factorization failed\n" <<endl;
			return;
		}

		beta0 = p1 / diag;

		q0 = y1[i];
		q1 = z1[i];
		diag += q0 * q1;

		if ( diag == 0.0f ) {
			Debug::debug(Debug::math,__FUNCTION__)<< "CLinCP_Square::RemoveClamped: updating factorization failed\n" <<endl;
			return;
		}

		d = 1.0f / diag;
		beta1 = q1 * d;

		clamped[i][i] = diag;
		diagonal[i] = d;

		for ( j = i+1; j < numClamped; j++ ) {

			d = clamped[i][j];

			d += p0 * z0[j];
			z0[j] -= beta0 * d;

			d += q0 * z1[j];
			z1[j] -= beta1 * d;

			clamped[i][j] = d;
		}

		for ( j = i+1; j < numClamped; j++ ) {

			d = clamped[j][i];

			y0[j] -= p0 * d;
			d += beta0 * y0[j];

			y1[j] -= q0 * d;
			d += beta1 * y1[j];

			clamped[j][i] = d;
		}
	}
	return;
}

/*
============
CLinCP_Square::CalcForceDelta

  modifies this->delta_f
============
*/
SMF_INLINE_FORCED void CLinCP_Square::CalcForceDelta( int d, float dir ) {
	int i;
	float *ptr;

	delta_f[d] = dir;

	if ( numClamped == 0 ) {
		return;
	}

	// get column d of matrix
	ptr = (float *) _allocafloat16( numClamped * sizeof( float ) );
	for ( i = 0; i < numClamped; i++ ) {
		ptr[i] = rowPtrs[i][d];
	}

	// solve force delta
	solveClamped( delta_f, ptr );

	// flip force delta based on direction
	if ( dir > 0.0f ) {
		ptr = delta_f.toFloatPtr();
		for ( i = 0; i < numClamped; i++ ) {
			ptr[i] = - ptr[i];
		}
	}
}

/*
============
CLinCP_Square::CalcAccelDelta

  modifies this->delta_a and uses this->delta_f
============
*/
SMF_INLINE_FORCED void CLinCP_Square::CalcAccelDelta( int d ) {
	int j;
	float dot;

	// only the not clamped variables, including the current variable, can have a change in acceleration
	for ( j = numClamped; j <= d; j++ ) {
		// only the clamped variables and the current variable have a force delta unequal zero
		SIMDProcessor->dot( dot, rowPtrs[j], delta_f.toFloatPtr(), numClamped );
		delta_a[j] = dot + rowPtrs[j][d] * delta_f[d];
	}
}

/*
============
CLinCP_Square::ChangeForce

  modifies this->f and uses this->delta_f
============
*/
SMF_INLINE_FORCED void CLinCP_Square::ChangeForce( int d, float step ) {
	// only the clamped variables and current variable have a force delta unequal zero
	SIMDProcessor->mulAdd( f.toFloatPtr(), step, delta_f.toFloatPtr(), numClamped );
	f[d] += step * delta_f[d];
}

/*
============
CLinCP_Square::ChangeAccel

  modifies this->a and uses this->delta_a
============
*/
SMF_INLINE_FORCED void CLinCP_Square::ChangeAccel( int d, float step ) {
	// only the not clamped variables, including the current variable, can have an acceleration unequal zero
	SIMDProcessor->mulAdd( a.toFloatPtr() + numClamped, step, delta_a.toFloatPtr() + numClamped, d - numClamped + 1 );
}

/*
============
CLinCP_Square::GetMaxStep
============
*/
void CLinCP_Square::GetMaxStep( int d, float dir, float &maxStep, int &limit, int &limitSide ) const {
	int i;
	float s;

	// default to a full step for the current variable
	if ( CMath::fabs( delta_a[d] ) > LCP_DELTA_ACCEL_EPSILON ) {
		maxStep = -a[d] / delta_a[d];
	} else {
		maxStep = 0.0f;
	}
	limit = d;
	limitSide = 0;

	// test the current variable
	if ( dir < 0.0f ) {
		if ( lo[d] != -CMath::INFINITY_FLOAT ) {
			s = ( lo[d] - f[d] ) / dir;
			if ( s < maxStep ) {
				maxStep = s;
				limitSide = -1;
			}
		}
	} else {
		if ( hi[d] != CMath::INFINITY_FLOAT ) {
			s = ( hi[d] - f[d] ) / dir;
			if ( s < maxStep ) {
				maxStep = s;
				limitSide = 1;
			}
		}
	}

	// test the clamped bounded variables
	for ( i = numUnbounded; i < numClamped; i++ ) {
		if ( delta_f[i] < -LCP_DELTA_FORCE_EPSILON ) {
			// if there is a low boundary
			if ( lo[i] != -CMath::INFINITY_FLOAT ) {
				s = ( lo[i] - f[i] ) / delta_f[i];
				if ( s < maxStep ) {
					maxStep = s;
					limit = i;
					limitSide = -1;
				}
			}
		} else if ( delta_f[i] > LCP_DELTA_FORCE_EPSILON ) {
			// if there is a high boundary
			if ( hi[i] != CMath::INFINITY_FLOAT ) {
				s = ( hi[i] - f[i] ) / delta_f[i];
				if ( s < maxStep ) {
					maxStep = s;
					limit = i;
					limitSide = 1;
				}
			}
		}
	}

	// test the not clamped bounded variables
	for ( i = numClamped; i < d; i++ ) {
		if ( side[i] == -1 ) {
			if ( delta_a[i] >= -LCP_DELTA_ACCEL_EPSILON ) {
				continue;
			}
		} else if ( side[i] == 1 ) {
			if ( delta_a[i] <= LCP_DELTA_ACCEL_EPSILON ) {
				continue;
			}
		} else {
			continue;
		}
		// ignore variables for which the force is not allowed to take any substantial value
		if ( lo[i] >= -LCP_BOUND_EPSILON && hi[i] <= LCP_BOUND_EPSILON ) {
			continue;
		}
		s = -a[i] / delta_a[i];
		if ( s < maxStep ) {
			maxStep = s;
			limit = i;
			limitSide = 0;
		}
	}
}

/*
============
CLinCP_Square::solve
============
*/
bool CLinCP_Square::solve( const CMatXD &o_m, CVecXD &o_x, const CVecXD &o_b, const CVecXD &o_lo, const CVecXD &o_hi, const int *o_boxIndex ) {
	int i, j, n, limit, limitSide, boxStartIndex;
	float dir, maxStep, dot, s;
	char *failed;

	// true when the matrix rows are 16 byte padded
	padded = ((o_m.getNumRows()+3)&~3) == o_m.getNumColumns();

	SMF_ASSERT( padded || o_m.getNumRows() == o_m.getNumColumns() );
	SMF_ASSERT( o_x.getSize() == o_m.getNumRows() );
	SMF_ASSERT( o_b.getSize() == o_m.getNumRows() );
	SMF_ASSERT( o_lo.getSize() == o_m.getNumRows() );
	SMF_ASSERT( o_hi.getSize() == o_m.getNumRows() );

	// allocate memory for permuted input
	f.setData( o_m.getNumRows(), VECX_ALLOCAFLOAT( o_m.getNumRows() ) );
	a.setData( o_b.getSize(), VECX_ALLOCAFLOAT( o_b.getSize() ) );
	b.setData( o_b.getSize(), VECX_ALLOCAFLOAT( o_b.getSize() ) );
	lo.setData( o_lo.getSize(), VECX_ALLOCAFLOAT( o_lo.getSize() ) );
	hi.setData( o_hi.getSize(), VECX_ALLOCAFLOAT( o_hi.getSize() ) );
	if ( o_boxIndex ) {
		boxIndex = (int *)_allocafloat16( o_x.getSize() * sizeof( int ) );
		memcpy( boxIndex, o_boxIndex, o_x.getSize() * sizeof( int ) );
	} else {
		boxIndex = NULL;
	}

	// we override the const on o_m here but on exit the matrix is unchanged
	m.setData( o_m.getNumRows(), o_m.getNumColumns(), const_cast<float *>(o_m[0]) );
	f.toZero();
	a.toZero();
	b = o_b;
	lo = o_lo;
	hi = o_hi;

	// pointers to the rows of m
	rowPtrs = (float **) _allocafloat16( m.getNumRows() * sizeof( float * ) );
	for ( i = 0; i < m.getNumRows(); i++ ) {
		rowPtrs[i] = m[i];
	}

	// tells if a variable is at the low boundary, high boundary or inbetween
	side = (int *) _allocafloat16( m.getNumRows() * sizeof( int ) );

	// index to keep track of the permutation
	permuted = (int *) _allocafloat16( m.getNumRows() * sizeof( int ) );
	for ( i = 0; i < m.getNumRows(); i++ ) {
		permuted[i] = i;
	}

	// permute input so all unbounded variables come first
	numUnbounded = 0;
	for ( i = 0; i < m.getNumRows(); i++ ) {
		if ( lo[i] == -CMath::INFINITY_FLOAT && hi[i] == CMath::INFINITY_FLOAT ) {
			if ( numUnbounded != i ) {
				Swap( numUnbounded, i );
			}
			numUnbounded++;
		}
	}

	// permute input so all variables using the boxIndex come last
	boxStartIndex = m.getNumRows();
	if ( boxIndex ) {
		for ( i = m.getNumRows() - 1; i >= numUnbounded; i-- ) {
			if ( boxIndex[i] >= 0 && ( lo[i] != -CMath::INFINITY_FLOAT || hi[i] != CMath::INFINITY_FLOAT ) ) {
				boxStartIndex--;
				if ( boxStartIndex != i ) {
					Swap( boxStartIndex, i );
				}
			}
		}
	}

	// sub matrix for factorization 
	clamped.setData( m.getNumRows(), m.getNumColumns(), MATX_ALLOCAFLOAT( m.getNumRows() * m.getNumColumns() ) );
	diagonal.setData( m.getNumRows(), VECX_ALLOCAFLOAT( m.getNumRows() ) );

	// all unbounded variables are clamped
	numClamped = numUnbounded;

	// if there are unbounded variables
	if ( numUnbounded ) {

		// factor and solve for unbounded variables
		if ( !FactorClamped() ) {
			Debug::debug(Debug::math,__FUNCTION__)<< "CLinCP_Square::solve: unbounded factorization failed\n" <<endl;
			return false;
		}
		solveClamped( f, b.toFloatPtr() );

		// if there are no bounded variables we are done
		if ( numUnbounded == m.getNumRows() ) {
			o_x = f;	// the vector is not permuted
			return true;
		}
	}

#ifdef IGNORE_UNSATISFIABLE_VARIABLES
	int numIgnored = 0;
#endif

	// allocate for delta force and delta acceleration
	delta_f.setData( m.getNumRows(), VECX_ALLOCAFLOAT( m.getNumRows() ) );
	delta_a.setData( m.getNumRows(), VECX_ALLOCAFLOAT( m.getNumRows() ) );

	// solve for bounded variables
	failed = NULL;
	for ( i = numUnbounded; i < m.getNumRows(); i++ ) {

		// once we hit the box start index we can initialize the low and high boundaries of the variables using the box index
		if ( i == boxStartIndex ) {
			for ( j = 0; j < boxStartIndex; j++ ) {
				o_x[permuted[j]] = f[j];
			}
			for ( j = boxStartIndex; j < m.getNumRows(); j++ ) {
				s = o_x[boxIndex[j]];
				if ( lo[j] != -CMath::INFINITY_FLOAT ) {
					lo[j] = - CMath::fabs( lo[j] * s );
				}
				if ( hi[j] != CMath::INFINITY_FLOAT ) {
					hi[j] = CMath::fabs( hi[j] * s );
				}
			}
		}

		// calculate acceleration for current variable
		SIMDProcessor->dot( dot, rowPtrs[i], f.toFloatPtr(), i );
		a[i] = dot - b[i];

		// if already at the low boundary
		if ( lo[i] >= -LCP_BOUND_EPSILON && a[i] >= -LCP_ACCEL_EPSILON ) {
			side[i] = -1;
			continue;
		}

		// if already at the high boundary
		if ( hi[i] <= LCP_BOUND_EPSILON && a[i] <= LCP_ACCEL_EPSILON ) {
			side[i] = 1;
			continue;
		}

		// if inside the clamped region
		if ( CMath::fabs( a[i] ) <= LCP_ACCEL_EPSILON ) {
			side[i] = 0;
			AddClamped( i );
			continue;
		}

		// drive the current variable into a valid region
		for ( n = 0; n < maxIterations; n++ ) {

			// direction to move
			if ( a[i] <= 0.0f ) {
				dir = 1.0f;
			} else {
				dir = -1.0f;
			}

			// calculate force delta
			CalcForceDelta( i, dir );

			// calculate acceleration delta: delta_a = m * delta_f;
			CalcAccelDelta( i );

			// maximum step we can take
			GetMaxStep( i, dir, maxStep, limit, limitSide );

			if ( maxStep <= 0.0f ) {
#ifdef IGNORE_UNSATISFIABLE_VARIABLES
				// ignore the current variable completely
				lo[i] = hi[i] = 0.0f;
				f[i] = 0.0f;
				side[i] = -1;
				numIgnored++;
#else
				failed = va( "invalid step size %.4f", maxStep );
#endif
				break;
			}

			// change force
			ChangeForce( i, maxStep );

			// change acceleration
			ChangeAccel( i, maxStep );

			// clamp/unclamp the variable that limited this step
			side[limit] = limitSide;
			switch( limitSide ) {
				case 0: {
					a[limit] = 0.0f;
					AddClamped( limit );
					break;
				}
				case -1: {
					f[limit] = lo[limit];
					if ( limit != i ) {
						RemoveClamped( limit );
					}
					break;
				}
				case 1: {
					f[limit] = hi[limit];
					if ( limit != i ) {
						RemoveClamped( limit );
					}
					break;
				}
			}

			// if the current variable limited the step we can continue with the next variable
			if ( limit == i ) {
				break;
			}
		}

		if ( n >= maxIterations ) {
			failed = va( "max iterations %d", maxIterations );
			break;
		}

		if ( failed ) {
			break;
		}
	}

#ifdef IGNORE_UNSATISFIABLE_VARIABLES
	if ( numIgnored ) {
		if ( lcp_showFailures.GetBool() ) {
			Debug::debug(Debug::math,__FUNCTION__)<< "CLCP_Symmetric::solve: "<< numIgnored <<" of " << m.getNumRows() - numUnbounded<<" bounded variables ignored\n" <<endl
		}
	}
#endif

	// if failed clear remaining forces
	if ( failed ) {
		//s if ( lcp_showFailures.GetBool() ) {
			Debug::debug(Debug::math,__FUNCTION__)<< "CLinCP_Square::solve: "<< failed<<" ( " <<m.getNumRows() - i << " of " << m.getNumRows() - numUnbounded <<" bounded variables ignored)"<<endl;
		//}
		for ( j = i; j < m.getNumRows(); j++ ) {
			f[j] = 0.0f;
		}
	}

#if defined(_DEBUG) && 0
	if ( !failed ) {
		// test whether or not the solution satisfies the complementarity conditions
		for ( i = 0; i < m.getNumRows(); i++ ) {
			a[i] = -b[i];
			for ( j = 0; j < m.getNumRows(); j++ ) {
				a[i] += rowPtrs[i][j] * f[j];
			}

			if ( f[i] == lo[i] ) {
				if ( lo[i] != hi[i] && a[i] < -LCP_ACCEL_EPSILON ) {
					int bah1 = 1;
				}
			} else if ( f[i] == hi[i] ) {
				if ( lo[i] != hi[i] && a[i] > LCP_ACCEL_EPSILON ) {
					int bah2 = 1;
				}
			} else if ( f[i] < lo[i] || f[i] > hi[i] || CMath::fabs( a[i] ) > 1.0f ) {
				int bah3 = 1;
			}
		}
	}
#endif

	// unpermute result
	for ( i = 0; i < f.getSize(); i++ ) {
		o_x[permuted[i]] = f[i];
	}

	// unpermute original matrix
	for ( i = 0; i < m.getNumRows(); i++ ) {
		for ( j = 0; j < m.getNumRows(); j++ ) {
			if ( permuted[j] == i ) {
				break;
			}
		}
		if ( i != j ) {
			m.swapColumns( i, j );
			swapElements( permuted[i], permuted[j] );
		}
	}

	return true;
}


//===============================================================
//                                                        M
//  CLCP_Symmetric                                      MrE
//                                                        E
//===============================================================
/**
 * 
 *
 * \ingroup SMF_Math
 *
 * \if pt_br
 * \brief	Ajuda na solução do problema de \ref CLinCP \ref
 *
 * \elseif us_en
 * \brief 	Helper to solve LCP problem
 *
 * \endif
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class CLCP_Symmetric : public CLinCP {
public:
	virtual bool	solve( const CMatXD &o_m, CVecXD &o_x, const CVecXD &o_b, const CVecXD &o_lo, const CVecXD &o_hi, const int *o_boxIndex );

private:
	CMatXD			m;					// original matrix
	CVecXD			b;					// right hand side
	CVecXD			lo, hi;				// low and high bounds
	CVecXD			f, a;				// force and acceleration
	CVecXD			delta_f, delta_a;	// delta force and delta acceleration
	CMatXD			clamped;			// LDLt factored sub matrix for clamped variables
	CVecXD			diagonal;			// reciprocal of diagonal of LDLt factored sub matrix for clamped variables
	CVecXD			solveCache1;		// intermediate result cached in solveClamped
	CVecXD			solveCache2;		// "
	int				numUnbounded;		// number of unbounded variables
	int				numClamped;			// number of clamped variables
	int				clampedChangeStart;	// lowest row/column changed in the clamped matrix during an iteration
	float **		rowPtrs;			// pointers to the rows of m
	int *			boxIndex;			// box index
	int *			side;				// tells if a variable is at the low boundary = -1, high boundary = 1 or inbetween = 0
	int *			permuted;			// index to keep track of the permutation
	bool			padded;				// set to true if the rows of the initial matrix are 16 byte padded

private:
	bool			FactorClamped();
	void			solveClamped( CVecXD &x, const float *b );
	void			Swap( int i, int j );
	void			AddClamped( int r, bool usesolveCache );
	void			RemoveClamped( int r );
	void			CalcForceDelta( int d, float dir );
	void			CalcAccelDelta( int d );
	void			ChangeForce( int d, float step );
	void			ChangeAccel( int d, float step );
	void			GetMaxStep( int d, float dir, float &maxStep, int &limit, int &limitSide ) const;
};

/*
============
CLCP_Symmetric::FactorClamped
============
*/
bool CLCP_Symmetric::FactorClamped() {

	clampedChangeStart = 0;

	for ( int i = 0; i < numClamped; i++ ) {
		memcpy( clamped[i], rowPtrs[i], numClamped * sizeof( float ) );
	}
	return SIMDProcessor->matX_LDLTFactor( clamped, diagonal, numClamped );
}

/*
============
CLCP_Symmetric::solveClamped
============
*/
void CLCP_Symmetric::solveClamped( CVecXD &x, const float *b ) {

	// solve L
	SIMDProcessor->matX_LowerTriangularsolve( clamped, solveCache1.toFloatPtr(), b, numClamped, clampedChangeStart );

	// solve D
	SIMDProcessor->mul( solveCache2.toFloatPtr(), solveCache1.toFloatPtr(), diagonal.toFloatPtr(), numClamped );

	// solve Lt
	SIMDProcessor->matX_LowerTriangularsolveTranspose( clamped, x.toFloatPtr(), solveCache2.toFloatPtr(), numClamped );

	clampedChangeStart = numClamped;
}

/*
============
CLCP_Symmetric::Swap
============
*/
void CLCP_Symmetric::Swap( int i, int j ) {

	if ( i == j ) {
		return;
	}

	swapElements( rowPtrs[i], rowPtrs[j] );
	m.swapColumns( i, j );
	b.swapElements( i, j );
	lo.swapElements( i, j );
	hi.swapElements( i, j );
	a.swapElements( i, j );
	f.swapElements( i, j );
	if ( boxIndex ) {
		swapElements( boxIndex[i], boxIndex[j] );
	}
	swapElements( side[i], side[j] );
	swapElements( permuted[i], permuted[j] );
}

/*
============
CLCP_Symmetric::AddClamped
============
*/
void CLCP_Symmetric::AddClamped( int r, bool usesolveCache ) {
	float d, dot;

	SMF_ASSERT( r >= numClamped );

	if ( numClamped < clampedChangeStart ) {
		clampedChangeStart = numClamped;
	}

	// add a row at the bottom and a column at the right of the factored
	// matrix for the clamped variables

	Swap( numClamped, r );

	// solve for v in L * v = rowPtr[numClamped]
	if ( usesolveCache ) {

		// the lower triangular solve was cached in solveClamped called by CalcForceDelta
		memcpy( clamped[numClamped], solveCache2.toFloatPtr(), numClamped * sizeof( float ) );
		// calculate row dot product
		SIMDProcessor->dot( dot, solveCache2.toFloatPtr(), solveCache1.toFloatPtr(), numClamped );

	} else {

		float *v = (float *) _allocafloat16( numClamped * sizeof( float ) );

		SIMDProcessor->matX_LowerTriangularsolve( clamped, v, rowPtrs[numClamped], numClamped );
		// add bottom row to L
		SIMDProcessor->mul( clamped[numClamped], v, diagonal.toFloatPtr(), numClamped );
		// calculate row dot product
		SIMDProcessor->dot( dot, clamped[numClamped], v, numClamped );
	}

	// update diagonal[numClamped]
	d = rowPtrs[numClamped][numClamped] - dot;

	if ( d == 0.0f ) {
		Debug::debug(Debug::math,__FUNCTION__)<< "CLCP_Symmetric::AddClamped: updating factorization failed\n" <<endl;
		numClamped++;
		return;
	}

	clamped[numClamped][numClamped] = d;
	diagonal[numClamped] = 1.0f / d;

	numClamped++;
}

/*
============
CLCP_Symmetric::RemoveClamped
============
*/
void CLCP_Symmetric::RemoveClamped( int r ) {
	int i, j, n;
	float *addSub, *original, *v, *ptr, *v1, *v2, dot;
	double sum, diag, newDiag, invNewDiag, p1, p2, alpha1, alpha2, beta1, beta2;

	SMF_ASSERT( r < numClamped );

	if ( r < clampedChangeStart ) {
		clampedChangeStart = r;
	}

	numClamped--;

	// no need to swap and update the factored matrix when the last row and column are removed
	if ( r == numClamped ) {
		return;
	}

	// swap the to be removed row/column with the last row/column
	Swap( r, numClamped );

	// update the factored matrix
	addSub = (float *) _allocafloat16( numClamped * sizeof( float ) );

	if ( r == 0 ) {

		if ( numClamped == 1 ) {
			diag = rowPtrs[0][0];
			if ( diag == 0.0f ) {
				Debug::debug(Debug::math,__FUNCTION__)<< "CLCP_Symmetric::RemoveClamped: updating factorization failed\n" <<endl;
				return;
			}
			clamped[0][0] = diag;
			diagonal[0] = 1.0f / diag;
			return;
		}

		// calculate the row/column to be added to the lower right sub matrix starting at (r, r)
		original = rowPtrs[numClamped];
		ptr = rowPtrs[r];
		addSub[0] = ptr[0] - original[numClamped];
		for ( i = 1; i < numClamped; i++ ) {
			addSub[i] = ptr[i] - original[i];
		}

	} else {

		v = (float *) _allocafloat16( numClamped * sizeof( float ) );

		// solve for v in L * v = rowPtr[r]
		SIMDProcessor->matX_LowerTriangularsolve( clamped, v, rowPtrs[r], r );

		// update removed row
		SIMDProcessor->mul( clamped[r], v, diagonal.toFloatPtr(), r );

		// if the last row/column of the matrix is updated
		if ( r == numClamped - 1 ) {
			// only calculate new diagonal
			SIMDProcessor->dot( dot, clamped[r], v, r );
			diag = rowPtrs[r][r] - dot;
			if ( diag == 0.0f ) {
				Debug::debug(Debug::math,__FUNCTION__)<< "CLCP_Symmetric::RemoveClamped: updating factorization failed\n" <<endl;
				return;
			}
			clamped[r][r] = diag;
			diagonal[r] = 1.0f / diag;
			return;
		}

		// calculate the row/column to be added to the lower right sub matrix starting at (r, r)
		for ( i = 0; i < r; i++ ) {
			v[i] = clamped[r][i] * clamped[i][i];
		}
		for ( i = r; i < numClamped; i++ ) {
			if ( i == r ) {
				sum = clamped[r][r];
			} else {
				sum = clamped[r][r] * clamped[i][r];
			}
			ptr = clamped[i];
			for ( j = 0; j < r; j++ ) {
				sum += ptr[j] * v[j];
			}
			addSub[i] = rowPtrs[r][i] - sum;
		}
	}

	// add row/column to the lower right sub matrix starting at (r, r)

	v1 = (float *) _allocafloat16( numClamped * sizeof( float ) );
	v2 = (float *) _allocafloat16( numClamped * sizeof( float ) );

	diag = CMath::SQRT_1OVER2;
	v1[r] = ( 0.5f * addSub[r] + 1.0f ) * diag;
	v2[r] = ( 0.5f * addSub[r] - 1.0f ) * diag;
	for ( i = r+1; i < numClamped; i++ ) {
		v1[i] = v2[i] = addSub[i] * diag;
	}

	alpha1 = 1.0f;
	alpha2 = -1.0f;

	// simultaneous update/downdate of the sub matrix starting at (r, r)
	n = clamped.getNumColumns();
	for ( i = r; i < numClamped; i++ ) {

		diag = clamped[i][i];
		p1 = v1[i];
		newDiag = diag + alpha1 * p1 * p1;

		if ( newDiag == 0.0f ) {
			Debug::debug(Debug::math,__FUNCTION__)<< "CLCP_Symmetric::RemoveClamped: updating factorization failed\n" <<endl;
			return;
		}

		alpha1 /= newDiag;
		beta1 = p1 * alpha1;
		alpha1 *= diag;

		diag = newDiag;
		p2 = v2[i];
		newDiag = diag + alpha2 * p2 * p2;

		if ( newDiag == 0.0f ) {
			Debug::debug(Debug::math,__FUNCTION__)<< "CLCP_Symmetric::RemoveClamped: updating factorization failed\n" <<endl;
			return;
		}

		clamped[i][i] = newDiag;
		diagonal[i] = invNewDiag = 1.0f / newDiag;

		alpha2 *= invNewDiag;
		beta2 = p2 * alpha2;
		alpha2 *= diag;

		// update column below diagonal (i,i)
		ptr = clamped.toFloatPtr() + i;

		for ( j = i+1; j < numClamped - 1; j += 2 ) {

			float sum0 = ptr[(j+0)*n];
			float sum1 = ptr[(j+1)*n];

			v1[j+0] -= p1 * sum0;
			v1[j+1] -= p1 * sum1;

			sum0 += beta1 * v1[j+0];
			sum1 += beta1 * v1[j+1];

			v2[j+0] -= p2 * sum0;
			v2[j+1] -= p2 * sum1;

			sum0 += beta2 * v2[j+0];
			sum1 += beta2 * v2[j+1];

			ptr[(j+0)*n] = sum0;
			ptr[(j+1)*n] = sum1;
		}

		for ( ; j < numClamped; j++ ) {

			sum = ptr[j*n];

			v1[j] -= p1 * sum;
			sum += beta1 * v1[j];

			v2[j] -= p2 * sum;
			sum += beta2 * v2[j];

			ptr[j*n] = sum;
		}
	}
}

/*
============
CLCP_Symmetric::CalcForceDelta

  modifies this->delta_f
============
*/
SMF_INLINE_FORCED void CLCP_Symmetric::CalcForceDelta( int d, float dir ) {
	int i;
	float *ptr;

	delta_f[d] = dir;

	if ( numClamped == 0 ) {
		return;
	}

	// solve force delta
	solveClamped( delta_f, rowPtrs[d] );

	// flip force delta based on direction
	if ( dir > 0.0f ) {
		ptr = delta_f.toFloatPtr();
		for ( i = 0; i < numClamped; i++ ) {
			ptr[i] = - ptr[i];
		}
	}
}

/*
============
CLCP_Symmetric::CalcAccelDelta

  modifies this->delta_a and uses this->delta_f
============
*/
SMF_INLINE_FORCED void CLCP_Symmetric::CalcAccelDelta( int d ) {
	int j;
	float dot;

	// only the not clamped variables, including the current variable, can have a change in acceleration
	for ( j = numClamped; j <= d; j++ ) {
		// only the clamped variables and the current variable have a force delta unequal zero
		SIMDProcessor->dot( dot, rowPtrs[j], delta_f.toFloatPtr(), numClamped );
		delta_a[j] = dot + rowPtrs[j][d] * delta_f[d];
	}
}

/*
============
CLCP_Symmetric::ChangeForce

  modifies this->f and uses this->delta_f
============
*/
SMF_INLINE_FORCED void CLCP_Symmetric::ChangeForce( int d, float step ) {
	// only the clamped variables and current variable have a force delta unequal zero
	SIMDProcessor->mulAdd( f.toFloatPtr(), step, delta_f.toFloatPtr(), numClamped );
	f[d] += step * delta_f[d];
}

/*
============
CLCP_Symmetric::ChangeAccel

  modifies this->a and uses this->delta_a
============
*/
SMF_INLINE_FORCED void CLCP_Symmetric::ChangeAccel( int d, float step ) {
	// only the not clamped variables, including the current variable, can have an acceleration unequal zero
	SIMDProcessor->mulAdd( a.toFloatPtr() + numClamped, step, delta_a.toFloatPtr() + numClamped, d - numClamped + 1 );
}

/*
============
CLCP_Symmetric::GetMaxStep
============
*/
void CLCP_Symmetric::GetMaxStep( int d, float dir, float &maxStep, int &limit, int &limitSide ) const {
	int i;
	float s;

	// default to a full step for the current variable
	if ( CMath::fabs( delta_a[d] ) > LCP_DELTA_ACCEL_EPSILON ) {
		maxStep = -a[d] / delta_a[d];
	} else {
		maxStep = 0.0f;
	}
	limit = d;
	limitSide = 0;

	// test the current variable
	if ( dir < 0.0f ) {
		if ( lo[d] != -CMath::INFINITY_FLOAT ) {
			s = ( lo[d] - f[d] ) / dir;
			if ( s < maxStep ) {
				maxStep = s;
				limitSide = -1;
			}
		}
	} else {
		if ( hi[d] != CMath::INFINITY_FLOAT ) {
			s = ( hi[d] - f[d] ) / dir;
			if ( s < maxStep ) {
				maxStep = s;
				limitSide = 1;
			}
		}
	}

	// test the clamped bounded variables
	for ( i = numUnbounded; i < numClamped; i++ ) {
		if ( delta_f[i] < -LCP_DELTA_FORCE_EPSILON ) {
			// if there is a low boundary
			if ( lo[i] != -CMath::INFINITY_FLOAT ) {
				s = ( lo[i] - f[i] ) / delta_f[i];
				if ( s < maxStep ) {
					maxStep = s;
					limit = i;
					limitSide = -1;
				}
			}
		} else if ( delta_f[i] > LCP_DELTA_FORCE_EPSILON ) {
			// if there is a high boundary
			if ( hi[i] != CMath::INFINITY_FLOAT ) {
				s = ( hi[i] - f[i] ) / delta_f[i];
				if ( s < maxStep ) {
					maxStep = s;
					limit = i;
					limitSide = 1;
				}
			}
		}
	}

	// test the not clamped bounded variables
	for ( i = numClamped; i < d; i++ ) {
		if ( side[i] == -1 ) {
			if ( delta_a[i] >= -LCP_DELTA_ACCEL_EPSILON ) {
				continue;
			}
		} else if ( side[i] == 1 ) {
			if ( delta_a[i] <= LCP_DELTA_ACCEL_EPSILON ) {
				continue;
			}
		} else {
			continue;
		}
		// ignore variables for which the force is not allowed to take any substantial value
		if ( lo[i] >= -LCP_BOUND_EPSILON && hi[i] <= LCP_BOUND_EPSILON ) {
			continue;
		}
		s = -a[i] / delta_a[i];
		if ( s < maxStep ) {
			maxStep = s;
			limit = i;
			limitSide = 0;
		}
	}
}

/*
============
CLCP_Symmetric::solve
============
*/
bool CLCP_Symmetric::solve( const CMatXD &o_m, CVecXD &o_x, const CVecXD &o_b, const CVecXD &o_lo, const CVecXD &o_hi, const int *o_boxIndex ) {
	int i, j, n, limit, limitSide, boxStartIndex;
	float dir, maxStep, dot, s;
	char *failed;

	// true when the matrix rows are 16 byte padded
	padded = ((o_m.getNumRows()+3)&~3) == o_m.getNumColumns();

	SMF_ASSERT( padded || o_m.getNumRows() == o_m.getNumColumns() );
	SMF_ASSERT( o_x.getSize() == o_m.getNumRows() );
	SMF_ASSERT( o_b.getSize() == o_m.getNumRows() );
	SMF_ASSERT( o_lo.getSize() == o_m.getNumRows() );
	SMF_ASSERT( o_hi.getSize() == o_m.getNumRows() );

	// allocate memory for permuted input
	f.setData( o_m.getNumRows(), VECX_ALLOCAFLOAT( o_m.getNumRows() ) );
	a.setData( o_b.getSize(), VECX_ALLOCAFLOAT( o_b.getSize() ) );
	b.setData( o_b.getSize(), VECX_ALLOCAFLOAT( o_b.getSize() ) );
	lo.setData( o_lo.getSize(), VECX_ALLOCAFLOAT( o_lo.getSize() ) );
	hi.setData( o_hi.getSize(), VECX_ALLOCAFLOAT( o_hi.getSize() ) );
	if ( o_boxIndex ) {
		boxIndex = (int *)_allocafloat16( o_x.getSize() * sizeof( int ) );
		memcpy( boxIndex, o_boxIndex, o_x.getSize() * sizeof( int ) );
	} else {
		boxIndex = NULL;
	}

	// we override the const on o_m here but on exit the matrix is unchanged
	m.setData( o_m.getNumRows(), o_m.getNumColumns(), const_cast<float *>(o_m[0]) );
	f.toZero();
	a.toZero();
	b = o_b;
	lo = o_lo;
	hi = o_hi;

	// pointers to the rows of m
	rowPtrs = (float **) _allocafloat16( m.getNumRows() * sizeof( float * ) );
	for ( i = 0; i < m.getNumRows(); i++ ) {
		rowPtrs[i] = m[i];
	}

	// tells if a variable is at the low boundary, high boundary or inbetween
	side = (int *) _allocafloat16( m.getNumRows() * sizeof( int ) );

	// index to keep track of the permutation
	permuted = (int *) _allocafloat16( m.getNumRows() * sizeof( int ) );
	for ( i = 0; i < m.getNumRows(); i++ ) {
		permuted[i] = i;
	}

	// permute input so all unbounded variables come first
	numUnbounded = 0;
	for ( i = 0; i < m.getNumRows(); i++ ) {
		if ( lo[i] == -CMath::INFINITY_FLOAT && hi[i] == CMath::INFINITY_FLOAT ) {
			if ( numUnbounded != i ) {
				Swap( numUnbounded, i );
			}
			numUnbounded++;
		}
	}

	// permute input so all variables using the boxIndex come last
	boxStartIndex = m.getNumRows();
	if ( boxIndex ) {
		for ( i = m.getNumRows() - 1; i >= numUnbounded; i-- ) {
			if ( boxIndex[i] >= 0 && ( lo[i] != -CMath::INFINITY_FLOAT || hi[i] != CMath::INFINITY_FLOAT ) ) {
				boxStartIndex--;
				if ( boxStartIndex != i ) {
					Swap( boxStartIndex, i );
				}
			}
		}
	}

	// sub matrix for factorization 
	clamped.setData( m.getNumRows(), m.getNumColumns(), MATX_ALLOCAFLOAT( m.getNumRows() * m.getNumColumns() ) );
	diagonal.setData( m.getNumRows(), VECX_ALLOCAFLOAT( m.getNumRows() ) );
	solveCache1.setData( m.getNumRows(), VECX_ALLOCAFLOAT( m.getNumRows() ) );
	solveCache2.setData( m.getNumRows(), VECX_ALLOCAFLOAT( m.getNumRows() ) );

	// all unbounded variables are clamped
	numClamped = numUnbounded;

	// if there are unbounded variables
	if ( numUnbounded ) {

		// factor and solve for unbounded variables
		if ( !FactorClamped() ) {
			Debug::debug(Debug::math,__FUNCTION__)<< "CLCP_Symmetric::solve: unbounded factorization failed\n" <<endl;
			return false;
		}
		solveClamped( f, b.toFloatPtr() );

		// if there are no bounded variables we are done
		if ( numUnbounded == m.getNumRows() ) {
			o_x = f;	// the vector is not permuted
			return true;
		}
	}

#ifdef IGNORE_UNSATISFIABLE_VARIABLES
	int numIgnored = 0;
#endif

	// allocate for delta force and delta acceleration
	delta_f.setData( m.getNumRows(), VECX_ALLOCAFLOAT( m.getNumRows() ) );
	delta_a.setData( m.getNumRows(), VECX_ALLOCAFLOAT( m.getNumRows() ) );

	// solve for bounded variables
	failed = NULL;
	for ( i = numUnbounded; i < m.getNumRows(); i++ ) {

		clampedChangeStart = 0;

		// once we hit the box start index we can initialize the low and high boundaries of the variables using the box index
		if ( i == boxStartIndex ) {
			for ( j = 0; j < boxStartIndex; j++ ) {
				o_x[permuted[j]] = f[j];
			}
			for ( j = boxStartIndex; j < m.getNumRows(); j++ ) {
				s = o_x[boxIndex[j]];
				if ( lo[j] != -CMath::INFINITY_FLOAT ) {
					lo[j] = - CMath::fabs( lo[j] * s );
				}
				if ( hi[j] != CMath::INFINITY_FLOAT ) {
					hi[j] = CMath::fabs( hi[j] * s );
				}
			}
		}

		// calculate acceleration for current variable
		SIMDProcessor->dot( dot, rowPtrs[i], f.toFloatPtr(), i );
		a[i] = dot - b[i];

		// if already at the low boundary
		if ( lo[i] >= -LCP_BOUND_EPSILON && a[i] >= -LCP_ACCEL_EPSILON ) {
			side[i] = -1;
			continue;
		}

		// if already at the high boundary
		if ( hi[i] <= LCP_BOUND_EPSILON && a[i] <= LCP_ACCEL_EPSILON ) {
			side[i] = 1;
			continue;
		}

		// if inside the clamped region
		if ( CMath::fabs( a[i] ) <= LCP_ACCEL_EPSILON ) {
			side[i] = 0;
			AddClamped( i, false );
			continue;
		}

		// drive the current variable into a valid region
		for ( n = 0; n < maxIterations; n++ ) {

			// direction to move
			if ( a[i] <= 0.0f ) {
				dir = 1.0f;
			} else {
				dir = -1.0f;
			}

			// calculate force delta
			CalcForceDelta( i, dir );

			// calculate acceleration delta: delta_a = m * delta_f;
			CalcAccelDelta( i );

			// maximum step we can take
			GetMaxStep( i, dir, maxStep, limit, limitSide );

			if ( maxStep <= 0.0f ) {
#ifdef IGNORE_UNSATISFIABLE_VARIABLES
				// ignore the current variable completely
				lo[i] = hi[i] = 0.0f;
				f[i] = 0.0f;
				side[i] = -1;
				numIgnored++;
#else
				failed = va( "invalid step size %.4f", maxStep );
#endif
				break;
			}

			// change force
			ChangeForce( i, maxStep );

			// change acceleration
			ChangeAccel( i, maxStep );

			// clamp/unclamp the variable that limited this step
			side[limit] = limitSide;
			switch( limitSide ) {
				case 0: {
					a[limit] = 0.0f;
					AddClamped( limit, ( limit == i ) );
					break;
				}
				case -1: {
					f[limit] = lo[limit];
					if ( limit != i ) {
						RemoveClamped( limit );
					}
					break;
				}
				case 1: {
					f[limit] = hi[limit];
					if ( limit != i ) {
						RemoveClamped( limit );
					}
					break;
				}
			}

			// if the current variable limited the step we can continue with the next variable
			if ( limit == i ) {
				break;
			}
		}

		if ( n >= maxIterations ) {
			failed = va( "max iterations %d", maxIterations );
			break;
		}

		if ( failed ) {
			break;
		}
	}

#ifdef IGNORE_UNSATISFIABLE_VARIABLES
	if ( numIgnored ) {
		if ( lcp_showFailures.GetBool() ) {
			Debug::debug(Debug::math,__FUNCTION__)<< "CLCP_Symmetric::solve: "<< numIgnored<< " of "<< m.getNumRows() - numUnbounded<<" bounded variables ignored"  <<endl;
		}
	}
#endif

	// if failed clear remaining forces
	if ( failed ) {
		//s if ( lcp_showFailures.GetBool() ) {
			Debug::debug(Debug::math,__FUNCTION__)<< "CLCP_Symmetric::solve: "<< failed<< " ( "<< m.getNumRows() - i<< " of "<<  m.getNumRows() - numUnbounded<<" bounded variables ignored)"<<endl;
		//}
		for ( j = i; j < m.getNumRows(); j++ ) {
			f[j] = 0.0f;
		}
	}

#if defined(_DEBUG) && 0
	if ( !failed ) {
		// test whether or not the solution satisfies the complementarity conditions
		for ( i = 0; i < m.getNumRows(); i++ ) {
			a[i] = -b[i];
			for ( j = 0; j < m.getNumRows(); j++ ) {
				a[i] += rowPtrs[i][j] * f[j];
			}

			if ( f[i] == lo[i] ) {
				if ( lo[i] != hi[i] && a[i] < -LCP_ACCEL_EPSILON ) {
					int bah1 = 1;
				}
			} else if ( f[i] == hi[i] ) {
				if ( lo[i] != hi[i] && a[i] > LCP_ACCEL_EPSILON ) {
					int bah2 = 1;
				}
			} else if ( f[i] < lo[i] || f[i] > hi[i] || CMath::fabs( a[i] ) > 1.0f ) {
				int bah3 = 1;
			}
		}
	}
#endif

	// unpermute result
	for ( i = 0; i < f.getSize(); i++ ) {
		o_x[permuted[i]] = f[i];
	}

	// unpermute original matrix
	for ( i = 0; i < m.getNumRows(); i++ ) {
		for ( j = 0; j < m.getNumRows(); j++ ) {
			if ( permuted[j] == i ) {
				break;
			}
		}
		if ( i != j ) {
			m.swapColumns( i, j );
			swapElements( permuted[i], permuted[j] );
		}
	}

	return true;
}


//===============================================================
//
//	CLinCP
//
//===============================================================

/*
============
CLinCP::AllocSquare
============
*/
CLinCP *CLinCP::AllocSquare() {
	CLinCP *lcp = new CLinCP_Square;
	lcp->setMaxIterations( 32 );
	return lcp;
}

/*
============
CLinCP::AllocSymmetric
============
*/
CLinCP *CLinCP::AllocSymmetric() {
	CLinCP *lcp = new CLCP_Symmetric;
	lcp->setMaxIterations( 32 );
	return lcp;
}

/*
============
CLinCP::~CLinCP
============
*/
CLinCP::~CLinCP() {
}

/*
============
CLinCP::setMaxIterations
============
*/
void CLinCP::setMaxIterations( int max ) {
	maxIterations = max;
}

/*
============
CLinCP::getMaxIterations
============
*/
int CLinCP::getMaxIterations() {
	return maxIterations;
}

} //end MATH
} //end SMF