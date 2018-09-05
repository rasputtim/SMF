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

#include "math/SMF_Regression.h"
#include "math/SMF_Vector.h"
#include "math/SMF_Matriz.h"
#include "structures/SMF_List.h"
#include "util/SMF_StringUtils.h"
#include "util/SMF_Debug.h"
#include "util/SMF_TestLib.h"
#include "sys/SMF_System.h"

//#pragma hdrstop
namespace SMF {
namespace MATH{
/* Fit the data (x_i, y_i) to the linear relationship 

   Y = c0 + c1 x

   returning, 

   c0, c1  --  coefficients
   cov00, cov01, cov11  --  variance-covariance matrix of c0 and c1,
   sumsq   --   sum of squares of residuals 

   This fit can be used in the case where the errors for the data are
   uknown, but assumed equal for all points. The resulting
   variance-covariance matrix estimates the error in the coefficients
   from the observed variance of the points around the best fit line.
*/

int CLinRegression::fit_linear (const float *x, const size_t xstride,
                const float *y, const size_t ystride,
                const size_t n,
                float *c0, float *c1,
                float *cov_00, float *cov_01, float *cov_11, float *sumsq)
{
  float m_x = 0, m_y = 0, m_dx2 = 0, m_dxdy = 0;

  size_t i;

  for (i = 0; i < n; i++)
    {
      m_x += (x[i * xstride] - m_x) / (i + 1.0);
      m_y += (y[i * ystride] - m_y) / (i + 1.0);
    }

  for (i = 0; i < n; i++)
    {
      const float dx = x[i * xstride] - m_x;
      const float dy = y[i * ystride] - m_y;

      m_dx2 += (dx * dx - m_dx2) / (i + 1.0);
      m_dxdy += (dx * dy - m_dxdy) / (i + 1.0);
    }

  /* In terms of y = a + b x */

  {
    float s2 = 0, d2 = 0;
    float b = m_dxdy / m_dx2;
    float a = m_y - m_x * b;

    *c0 = a;
    *c1 = b;

    /* Compute chi^2 = \sum (y_i - (a + b * x_i))^2 */

    for (i = 0; i < n; i++)
      {
        const float dx = x[i * xstride] - m_x;
        const float dy = y[i * ystride] - m_y;
        const float d = dy - b * dx;
        d2 += d * d;
      }

    s2 = d2 / (n - 2.0);        /* chisq per degree of freedom */

    *cov_00 = s2 * (1.0 / n) * (1 + m_x * m_x / m_dx2);
    *cov_11 = s2 * 1.0 / (n * m_dx2);

    *cov_01 = s2 * (-m_x) / (n * m_dx2);

    *sumsq = d2;
  }

  return SMF_OK;
}


/* Fit the weighted data (x_i, w_i, y_i) to the linear relationship 

   Y = c0 + c1 x

   returning, 

   c0, c1  --  coefficients
   s0, s1  --  the standard deviations of c0 and c1,
   r       --  the correlation coefficient between c0 and c1,
   chisq   --  weighted sum of squares of residuals */

int CLinRegression::fitWlinear (const float *x, const size_t xstride,
                 const float *w, const size_t wstride,
                 const float *y, const size_t ystride,
                 const size_t n,
                 float *c0, float *c1,
                 float *cov_00, float *cov_01, float *cov_11,
                 float *chisq)
{

  /* compute the weighted means and weighted deviations from the means */

  /* wm denotes a "weighted mean", wm(f) = (sum_i w_i f_i) / (sum_i w_i) */

  float W = 0, wm_x = 0, wm_y = 0, wm_dx2 = 0, wm_dxdy = 0;

  size_t i;

  for (i = 0; i < n; i++)
    {
      const float wi = w[i * wstride];

      if (wi > 0)
        {
          W += wi;
          wm_x += (x[i * xstride] - wm_x) * (wi / W);
          wm_y += (y[i * ystride] - wm_y) * (wi / W);
        }
    }

  W = 0;                        /* reset the total weight */

  for (i = 0; i < n; i++)
    {
      const float wi = w[i * wstride];

      if (wi > 0)
        {
          const float dx = x[i * xstride] - wm_x;
          const float dy = y[i * ystride] - wm_y;

          W += wi;
          wm_dx2 += (dx * dx - wm_dx2) * (wi / W);
          wm_dxdy += (dx * dy - wm_dxdy) * (wi / W);
        }
    }

  /* In terms of y = a + b x */

  {
    float d2 = 0;
    float b = wm_dxdy / wm_dx2;
    float a = wm_y - wm_x * b;

    *c0 = a;
    *c1 = b;

    *cov_00 = (1 / W) * (1 + wm_x * wm_x / wm_dx2);
    *cov_11 = 1 / (W * wm_dx2);

    *cov_01 = -wm_x / (W * wm_dx2);

    /* Compute chi^2 = \sum w_i (y_i - (a + b * x_i))^2 */

    for (i = 0; i < n; i++)
      {
        const float wi = w[i * wstride];

        if (wi > 0)
          {
            const float dx = x[i * xstride] - wm_x;
            const float dy = y[i * ystride] - wm_y;
            const float d = dy - b * dx;
            d2 += wi * d * d;
          }
      }

    *chisq = d2;
  }

  return SMF_OK;
}



int CLinRegression::fitLinearEst (const float x,
                    const float c0, const float c1,
                    const float c00, const float c01, const float c11,
                    float *y, float *y_err)
{
  *y = c0 + c1 * x;
  *y_err = CMath::sqrt(c00 + x * (2 * c01 + c11 * x));
  return SMF_OK;
}


int CLinRegression::fitMul (const float *x, const size_t xstride,
             const float *y, const size_t ystride,
             const size_t n, 
             float *c1, float *cov_11, float *sumsq)
{
  float m_x = 0, m_y = 0, m_dx2 = 0, m_dxdy = 0;

  size_t i;

  for (i = 0; i < n; i++)
    {
      m_x += (x[i * xstride] - m_x) / (i + 1.0);
      m_y += (y[i * ystride] - m_y) / (i + 1.0);
    }

  for (i = 0; i < n; i++)
    {
      const float dx = x[i * xstride] - m_x;
      const float dy = y[i * ystride] - m_y;

      m_dx2 += (dx * dx - m_dx2) / (i + 1.0);
      m_dxdy += (dx * dy - m_dxdy) / (i + 1.0);
    }

  /* In terms of y =  b x */

  {
    float s2 = 0, d2 = 0;
    float b = (m_x * m_y + m_dxdy) / (m_x * m_x + m_dx2);

    *c1 = b;

    /* Compute chi^2 = \sum (y_i -  b * x_i)^2 */

    for (i = 0; i < n; i++)
      {
        const float dx = x[i * xstride] - m_x;
        const float dy = y[i * ystride] - m_y;
        const float d = (m_y - b * m_x) + dy - b * dx;
        d2 += d * d;
      }

    s2 = d2 / (n - 1.0);        /* chisq per degree of freedom */

    *cov_11 = s2 * 1.0 / (n * (m_x * m_x + m_dx2));

    *sumsq = d2;
  }

  return SMF_OK;
}


int CLinRegression::fitWmul (const float *x, const size_t xstride,
              const float *w, const size_t wstride,
              const float *y, const size_t ystride,
              const size_t n, 
              float *c1, float *cov_11, float *chisq)
{

  /* compute the weighted means and weighted deviations from the means */

  /* wm denotes a "weighted mean", wm(f) = (sum_i w_i f_i) / (sum_i w_i) */

  float W = 0, wm_x = 0, wm_y = 0, wm_dx2 = 0, wm_dxdy = 0;

  size_t i;

  for (i = 0; i < n; i++)
    {
      const float wi = w[i * wstride];

      if (wi > 0)
        {
          W += wi;
          wm_x += (x[i * xstride] - wm_x) * (wi / W);
          wm_y += (y[i * ystride] - wm_y) * (wi / W);
        }
    }

  W = 0;                        /* reset the total weight */

  for (i = 0; i < n; i++)
    {
      const float wi = w[i * wstride];

      if (wi > 0)
        {
          const float dx = x[i * xstride] - wm_x;
          const float dy = y[i * ystride] - wm_y;

          W += wi;
          wm_dx2 += (dx * dx - wm_dx2) * (wi / W);
          wm_dxdy += (dx * dy - wm_dxdy) * (wi / W);
        }
    }

  /* In terms of y = b x */

  {
    float d2 = 0;
    float b = (wm_x * wm_y + wm_dxdy) / (wm_x * wm_x + wm_dx2);

    *c1 = b;

    *cov_11 = 1 / (W * (wm_x * wm_x + wm_dx2));

    /* Compute chi^2 = \sum w_i (y_i - b * x_i)^2 */

    for (i = 0; i < n; i++)
      {
        const float wi = w[i * wstride];

        if (wi > 0)
          {
            const float dx = x[i * xstride] - wm_x;
            const float dy = y[i * ystride] - wm_y;
            const float d = (wm_y - b * wm_x) + (dy - b * dx);
            d2 += wi * d * d;
          }
      }

    *chisq = d2;
  }

  return SMF_OK;
}

int CLinRegression::fitMul_est (const float x, 
                 const float c1, const float c11, 
                 float *y, float *y_err)
{
  *y = c1 * x;
  *y_err = CMath::sqrt(c11) * fabs (x);
  return SMF_OK;
}


void CLinRegression::TestRegression(){
	size_t norris_n = 36;

float norris_x[] = { 0.2, 337.4, 118.2, 884.6, 10.1, 226.5, 666.3, 996.3,
                      448.6, 777.0, 558.2, 0.4, 0.6, 775.5, 666.9, 338.0, 
                      447.5, 11.6, 556.0, 228.1, 995.8, 887.6, 120.2, 0.3, 
                      0.3, 556.8, 339.1, 887.2, 999.0, 779.0, 11.1, 118.3,
                      229.2, 669.1, 448.9, 0.5 } ;

float norris_y[] = { 0.1, 338.8, 118.1, 888.0, 9.2, 228.1, 668.5, 998.5,
                      449.1, 778.9, 559.2, 0.3, 0.1, 778.1, 668.8, 339.3, 
                      448.9, 10.8, 557.7, 228.3, 998.0, 888.8, 119.6, 0.3, 
                      0.6, 557.6, 339.3, 888.0, 998.5, 778.9, 10.2, 117.6,
                      228.9, 668.4, 449.2, 0.2};

size_t noint1_n = 11;
float noint1_x[] = { 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70 };
float noint1_y[] = { 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140};

size_t noint2_n = 3;
float noint2_x[] = { 4, 5, 6 } ;
float noint2_y[] = { 3, 4, 4 } ;


 float x[1000], y[1000], w[1000];

  size_t xstride = 2, wstride = 3, ystride = 5;
  size_t i;

  for (i = 0; i < norris_n; i++) 
    {
      x[i*xstride] = norris_x[i];
      w[i*wstride] = 1.0;
      y[i*ystride] = norris_y[i];
    }

  MATH::ieee_env_setup();

  {
    float c0, c1, cov00, cov01, cov11, sumsq;
       
    float expected_c0 = -0.262323073774029;
    float expected_c1 =  1.00211681802045; 
    float expected_cov00 = CMath::pow(0.232818234301152, 2.0);
    float expected_cov01 = -7.74327536339570e-05;  /* computed from octave */
    float expected_cov11 = CMath::pow(0.429796848199937E-03, 2.0);
    float expected_sumsq = 26.6173985294224;
    
    fit_linear (x, xstride, y, ystride, norris_n, 
                    &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
    
    /* fitWlinear (x, xstride, w, wstride, y, ystride, norris_n, 
                     &c0, &c1, &cov00, &cov01, &cov11, &sumsq); */
  
    testlib_test_rel (c0, expected_c0, 1e-10, "norris fit_linear c0") ;
    testlib_test_rel (c1, expected_c1, 1e-10, "norris fit_linear c1") ;
    testlib_test_rel (cov00, expected_cov00, 1e-10, "norris fit_linear cov00") ;
    testlib_test_rel (cov01, expected_cov01, 1e-10, "norris fit_linear cov01") ;
    testlib_test_rel (cov11, expected_cov11, 1e-10, "norris fit_linear cov11") ;
    testlib_test_rel (sumsq, expected_sumsq, 1e-10, "norris fit_linear sumsq") ;
  }

  {
    float c0, c1, cov00, cov01, cov11, sumsq;
       
    float expected_c0 = -0.262323073774029;
    float expected_c1 =  1.00211681802045; 
    float expected_cov00 = 6.92384428759429e-02;  /* computed from octave */
    float expected_cov01 = -9.89095016390515e-05; /* computed from octave */
    float expected_cov11 = 2.35960747164148e-07;  /* computed from octave */
    float expected_sumsq = 26.6173985294224;
    
    fitWlinear (x, xstride, w, wstride, y, ystride, norris_n, 
                     &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
  
    testlib_test_rel (c0, expected_c0, 1e-10, "norris fitWlinear c0") ;
    testlib_test_rel (c1, expected_c1, 1e-10, "norris fitWlinear c1") ;
    testlib_test_rel (cov00, expected_cov00, 1e-10, "norris fitWlinear cov00") ;
    testlib_test_rel (cov01, expected_cov01, 1e-10, "norris fitWlinear cov01") ;
    testlib_test_rel (cov11, expected_cov11, 1e-10, "norris fitWlinear cov11") ;
    testlib_test_rel (sumsq, expected_sumsq, 1e-10, "norris fitWlinear sumsq") ;
  }

  for (i = 0; i < noint1_n; i++) 
    {
      x[i*xstride] = noint1_x[i];
      w[i*wstride] = 1.0;
      y[i*ystride] = noint1_y[i];
    }

  {
    float c1, cov11, sumsq;
       
    float expected_c1 = 2.07438016528926; 
    float expected_cov11 = CMath::pow(0.165289256198347E-01, 2.0);  
    float expected_sumsq = 127.272727272727;
    
    fitMul (x, xstride, y, ystride, noint1_n, &c1, &cov11, &sumsq);
  
    testlib_test_rel (c1, expected_c1, 1e-10, "noint1 fitMul c1") ;
    testlib_test_rel (cov11, expected_cov11, 1e-10, "noint1 fitMul cov11") ;
    testlib_test_rel (sumsq, expected_sumsq, 1e-10, "noint1 fitMul sumsq") ;
  }

  {
    float c1, cov11, sumsq;
       
    float expected_c1 = 2.07438016528926; 
    float expected_cov11 = 2.14661371686165e-05; /* computed from octave */
    float expected_sumsq = 127.272727272727;
    
    fitWmul (x, xstride, w, wstride, y, ystride, noint1_n, &c1, &cov11, &sumsq);

    testlib_test_rel (c1, expected_c1, 1e-10, "noint1 fitWmul c1") ;
    testlib_test_rel (cov11, expected_cov11, 1e-10, "noint1 fitWmul cov11") ;
    testlib_test_rel (sumsq, expected_sumsq, 1e-10, "noint1 fitWmul sumsq") ;
  }


  for (i = 0; i < noint2_n; i++) 
    {
      x[i*xstride] = noint2_x[i];
      w[i*wstride] = 1.0;
      y[i*ystride] = noint2_y[i];
    }

  {
    float c1, cov11, sumsq;
       
    float expected_c1 = 0.727272727272727; 
    float expected_cov11 = CMath::pow(0.420827318078432E-01, 2.0);  
    float expected_sumsq = 0.272727272727273;
    
    fitMul (x, xstride, y, ystride, noint2_n, &c1, &cov11, &sumsq);
  
    testlib_test_rel (c1, expected_c1, 1e-10, "noint2 fitMul c1") ;
    testlib_test_rel (cov11, expected_cov11, 1e-10, "noint2 fitMul cov11") ;
    testlib_test_rel (sumsq, expected_sumsq, 1e-10, "noint2 fitMul sumsq") ;
  }

  {
    float c1, cov11, sumsq;
       
    float expected_c1 = 0.727272727272727; 
    float expected_cov11 = 1.29870129870130e-02 ; /* computed from octave */
    float expected_sumsq = 0.272727272727273;
    
    fitWmul (x, xstride, w, wstride, y, ystride, noint2_n, &c1, &cov11, &sumsq);

    testlib_test_rel (c1, expected_c1, 1e-10, "noint2 fitWmul c1") ;
    testlib_test_rel (cov11, expected_cov11, 1e-10, "noint2 fitWmul cov11") ;
    testlib_test_rel (sumsq, expected_sumsq, 1e-10, "noint2 fitWmul sumsq") ;
  }

  /* now summarize the results */

  exit (testlib_test_summary ());










}


int CLinRegression::testRegression2 (void)
{
  int i, n = 4;
  float x[4] = { 1970, 1980, 1990, 2000 };
  float y[4] = {   12,   11,   14,   13 };
  float w[4] = {  0.1,  0.2,  0.3,  0.4 };

  float c0, c1, cov00, cov01, cov11, chisq;

  fitWlinear (x, 1, w, 1, y, 1, n, 
                   &c0, &c1, &cov00, &cov01, &cov11, 
                   &chisq);

  printf ("# best fit: Y = %g + %g X\n", c0, c1);
  printf ("# covariance matrix:\n");
  printf ("# [ %g, %g\n#   %g, %g]\n", 
          cov00, cov01, cov01, cov11);
  printf ("# chisq = %g\n", chisq);

  for (i = 0; i < n; i++)
    printf ("data: %g %g %g\n", 
                   x[i], y[i], 1/sqrt(w[i]));

  printf ("\n");

  for (i = -30; i < 130; i++)
    {
      float xf = x[0] + (i/100.0) * (x[n-1] - x[0]);
      float yf, yf_err;

      fitLinearEst (xf, 
                          c0, c1, 
                          cov00, cov01, cov11, 
                          &yf, &yf_err);

      printf ("fit: %g %g\n", xf, yf);
      printf ("hi : %g %g\n", xf, yf + yf_err);
      printf ("lo : %g %g\n", xf, yf - yf_err);
    }
  return 0;
}



} //end MATH
} //end SMF