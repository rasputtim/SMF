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

#ifndef _SMF__MATH_LCP_H__
#define _SMF__MATH_LCP_H__
#include "../SMF_Config.h"
namespace SMF{
namespace MATH{
	class CMatXD;
	class CVecXD;

/**
 * \class 	CLinCP
 *
 * \ingroup SMF_Math
 *
 * \if pt_br
 * \brief	Solucionador do problema de complementaridade Linear Misto (LCP) em Caixas 
 * \see http://en.wikipedia.org/wiki/Mixed_linear_complementarity_problem , http://en.wikipedia.org/wiki/Linear_complementarity_problem
 *  \p
 *  \p A é uma matriz de dimensões n*n e x,b, lo, hi são vetores de dimensão n.
    \p Solucione: Ax = b + t, onde t é um vetor de dimensão n, com condição 
	complementar: (x[i] - lo[i]) * (x[i] - hi[i]) * t[i] = 0
	\p de forma que para cada 0 <=i < n, uma das seguintes se aplica:
	\p 1. lo[i] < x[i] < hi[i], t[i] == 0
    \p 2. x[i] == lo[i], t[i] >= 0
    \p 3. x[i] == hi[i], t[i] <= 0
	\p variáveis parcialmente limitadas, ou ilimitadas podem ter lo[i] e/ou hi[i] 
	CMath::INFITITY_FLOAT negativadas/positivadas respectivamente.
    \p Se boxIndex != NULL e boxIndex[i] != -1 então

    \p lo[i] = - fabs( lo[i] * x[boxIndex[i]] )
    \p hi[i] = fabs( hi[i] * x[boxIndex[i]] )
	\p boxIndex[boxIndex[i]] deve ser -1
    \p Antes de calcular qualquer dos limites x[i[ com boxIndex[i] != -1 o solucionador
  calcula todos os limites x[i] e todos x[i] com boxIndex[i] == -1.

 * \elseif us_en
 * \brief 	  Box Constrained Mixed Linear Complementarity Problem solver

  A is a matrix of dimension n*n and x, b, lo, hi are vectors of dimension n

  solve: Ax = b + t, where t is a vector of dimension n, with
  complementarity condition: (x[i] - lo[i]) * (x[i] - hi[i]) * t[i] = 0
  such that for each 0 <= i < n one of the following holds:

    1. lo[i] < x[i] < hi[i], t[i] == 0
    2. x[i] == lo[i], t[i] >= 0
    3. x[i] == hi[i], t[i] <= 0

  Partly bounded or unbounded variables can have lo[i] and/or hi[i]
  set to negative/positive CMath::INFITITY_FLOAT respectively.

  If boxIndex != NULL and boxIndex[i] != -1 then

    lo[i] = - fabs( lo[i] * x[boxIndex[i]] )
    hi[i] = fabs( hi[i] * x[boxIndex[i]] )
	boxIndex[boxIndex[i]] must be -1
  
  Before calculating any of the bounded x[i] with boxIndex[i] != -1 the
  solver calculates all unbounded x[i] and all x[i] with boxIndex[i] == -1.

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
class SMF_API  CLinCP {
public:
	static CLinCP *	AllocSquare();		// A must be a square matrix
	static CLinCP *	AllocSymmetric();		// A must be a symmetric matrix

	virtual			~CLinCP();

	virtual bool	solve( const CMatXD &A, CVecXD &x, const CVecXD &b, const CVecXD &lo, const CVecXD &hi, const int *boxIndex = NULL ) = 0;
	virtual void	setMaxIterations( int max );
	virtual int		getMaxIterations();

protected:
	int				maxIterations;
};

}  //end MATH
}   //end SMF
#endif /* !__MATH_LCP_H__ */
