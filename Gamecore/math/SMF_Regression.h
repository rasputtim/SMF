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
#include "../sys/SMF_System.h"
#ifndef _SMF__MATH_LIN_REGRESSION_H__
#define _SMF__MATH_LIN_REGRESSION_H__

namespace SMF {
namespace MATH{


/**
 * \class CLinRegression
 *
 * \ingroup SMF_Math
 *
\see	http://pt.wikipedia.org/wiki/Regress%C3%A3o_linear

\brief	Linear Regression.  
This chapter describes routines for performing least squares fits to experimental data using linear combinations of functions. 

The data may be weighted or unweighted, i.e. with known or unknown errors. For weighted data the functions compute the best fit parameters and their associated covariance matrix. For unweighted data the covariance matrix is estimated from the scatter of the points, giving a variance-covariance matrix.

The functions are divided into separate versions for simple one- or two-parameter regression and multiple-parameter fits.
Overview

Least-squares fits are found by minimizing \chi^2 (chi-squared), the weighted sum of squared residuals over n experimental datapoints (x_i, y_i) for the model Y(c,x),

\chi^2 = \sum_i w_i (y_i - Y(c, x_i))^2

The p parameters of the model are c = {c_0, c_1, …}. The weight factors w_i are given by w_i = 1/\sigma_i^2, where \sigma_i is the experimental error on the data-point y_i. The errors are assumed to be Gaussian and uncorrelated. For unweighted data the chi-squared sum is computed without any weight factors.

The fitting routines return the best-fit parameters c and their p \times p covariance matrix. The covariance matrix measures the statistical errors on the best-fit parameters resulting from the errors on the data, \sigma_i, and is defined as C_{ab} = <\delta c_a \delta c_b> where < > denotes an average over the Gaussian error distributions of the underlying datapoints.

The covariance matrix is calculated by error propagation from the data errors \sigma_i. The change in a fitted parameter \delta c_a caused by a small change in the data \delta y_i is given by

\delta c_a = \sum_i (dc_a/dy_i) \delta y_i

allowing the covariance matrix to be written in terms of the errors on the data,

C_{ab} = \sum_{i,j} (dc_a/dy_i) (dc_b/dy_j) <\delta y_i \delta y_j>

For uncorrelated data the fluctuations of the underlying datapoints satisfy <\delta y_i \delta y_j> = \sigma_i^2 \delta_{ij}, giving a corresponding parameter covariance matrix of

C_{ab} = \sum_i (1/w_i) (dc_a/dy_i) (dc_b/dy_i) 

When computing the covariance matrix for unweighted data, i.e. data with unknown errors, the weight factors w_i in this sum are replaced by the single estimate w = 1/\sigma^2, where \sigma^2 is the computed variance of the residuals about the best-fit model, \sigma^2 = \sum (y_i - Y(c,x_i))^2 / (n-p). This is referred to as the variance-covariance matrix.

The standard deviations of the best-fit parameters are given by the square root of the corresponding diagonal elements of the covariance matrix, \sigma_{c_a} = \sqrt{C_{aa}}. The correlation coefficient of the fit parameters c_a and c_b is given by \rho_{ab} = C_{ab} / \sqrt{C_{aa} C_{bb}}. 

* \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CLinRegression 
{
public: 
/*
\brief The functions described in this section can be used to perform least-squares fits to a straight line model, Y(c,x) = c_0 + c_1 x.
\see \see http://mathworld.wolfram.com/LeastSquaresFitting.html  */

/*
\brief This function computes the best-fit linear regression coefficients (c0,c1) of the model Y = c_0 + c_1 X for the dataset (x, y), two vectors of length n with strides xstride and ystride. The errors on y are assumed unknown so the variance-covariance matrix for the parameters (c0, c1) is estimated from the scatter of the points around the best-fit line and returned via the parameters (cov00, cov01, cov11). 
The sum of squares of the residuals from the best-fit line is returned in sumsq. 
*/
static int fit_linear (const float * x, const size_t xstride,
                    const float * y, const size_t ystride,
                    const size_t n,
                    float * c0, float * c1, 
                    float * cov00, float * cov01, float * cov11, 
                    float * sumsq);

/**
    This function computes the best-fit linear regression coefficients (c0,c1) of the model Y = c_0 + c_1 X for the weighted dataset (x, y), two vectors of length n with strides xstride and ystride. 
	The vector w, of length n and stride wstride, specifies the weight of each datapoint. The weight is the reciprocal of the variance for each datapoint in y.

    The covariance matrix for the parameters (c0, c1) is computed using the weights and returned via the parameters (cov00, cov01, cov11). The weighted sum of squares of the residuals from the best-fit line, \chi^2, is returned in chisq. 
*/
static int fitWlinear (const float * x, const size_t xstride,
                     const float * w, const size_t wstride,
                     const float * y, const size_t ystride,
                     const size_t n,
                     float * c0, float * c1, 
                     float * cov00, float * cov01, float * cov11, 
                     float * chisq);
/**
This function uses the best-fit linear regression coefficients c0, c1 and their covariance cov00, cov01, cov11 to compute the fitted function y and its standard deviation y_err for the model Y = c_0 + c_1 X at the point x. 
\param [out] y
\param [out] y_err
**/
static int fitLinearEst (const float x, 
                    const float c0, const float c1, 
                    const float c00, const float c01, const float c11,
                    float *y, float *y_err);

/**
\brief The functions described in this section can be used to perform least-squares fits to a straight line model without a constant term, Y = c_1 X.
**/

/**
\brief     This function computes the best-fit linear regression coefficient c1 of the model Y = c_1 X for the datasets (x, y), two vectors of length n with strides xstride and ystride. The errors on y are assumed unknown so the variance of the parameter c1 is estimated from the scatter of the points around the best-fit line and returned via the parameter cov11. The sum of squares of the residuals from the best-fit line is returned in sumsq. 


**/
static int fitMul (const float * x, const size_t xstride,
                 const float * y, const size_t ystride,
                 const size_t n,
                 float * c1, 
                 float * cov11, 
                 float * sumsq);
/**
\brief     This function computes the best-fit linear regression coefficient c1 of the model Y = c_1 X for the weighted datasets (x, y), 
 two vectors of length n with strides xstride and ystride. The vector w, of length n and stride wstride, specifies the weight of each datapoint. The weight is the reciprocal of the variance for each datapoint in y.

    The variance of the parameter c1 is computed using the weights and returned via the parameter cov11. The weighted sum of squares of the residuals from the best-fit line, \chi^2, is returned in chisq. 

**/
static int fitWmul (const float * x, const size_t xstride,
                  const float * w, const size_t wstride,
                  const float * y, const size_t ystride,
                  const size_t n,
                  float * c1, 
                  float * cov11, 
                  float * sumsq);

/**
\brief This function uses the best-fit linear regression coefficient c1 and its covariance cov11 to compute the fitted function y and its standard deviation y_err for the model Y = c_1 X at the point x. 
\param [out] y
\param [out] y_err
**/
static int fitMul_est (const float x, const float c1, const float c11,
                 float *y, float *y_err);

static void TestRegression();
/**
The following program computes a least squares straight-line fit to a simple dataset, and outputs the best-fit line and its associated one standard-deviation error bars.
The result:
# best fit: Y = -106.6 + 0.06 X
# covariance matrix:
# [ 39602, -19.9
#   -19.9, 0.01]
# chisq = 0.8
and data of course

**/
static int testRegression2 (void);


};





} //end MATH
} //end SMF
#endif /* !__MATH_INTERPOLATE_H__ */
