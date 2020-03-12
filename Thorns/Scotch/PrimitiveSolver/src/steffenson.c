/* roots/steffenson.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Reid Priedhorsky, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* steffenson.c -- steffenson root finding algorithm 

   This is Newton's method with an Aitken "delta-squared"
   acceleration of the iterates. This can improve the convergence on
   multiple roots where the ordinary Newton algorithm is slow.

   x[i+1] = x[i] - f(x[i]) / f'(x[i])

   x_accelerated[i] = x[i] - (x[i+1] - x[i])**2 / (x[i+2] - 2*x[i+1] + x[i])

   We can only use the accelerated estimate after three iterations,
   and use the unaccelerated value until then.

 */

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "steffenson.h"

typedef struct
  {
    double f, df;
    double x;
    double x_1;
    double x_2;
    int count;
  }
steffenson_state_t;

static int steffenson_init (void * vstate, gsl_function_fdf * fdf, double * root);
static int steffenson_iterate (void * vstate, gsl_function_fdf * fdf, double * root);
static const char *gsl_errtxt[] = {"GSL_SUCCESS", "GSL_FAILURE", "GSL_CONTINUE", "GSL_EDOM", "GSL_ERANGE", "GSL_EFAULT", "GSL_EINVAL", "GSL_EFAILED", "GSL_EFACTOR", "GSL_ESANITY", "GSL_ENOMEM", "GSL_EBADFUNC", "GSL_ERUNAWAY", "GSL_EMAXITER", "GSL_EZERODIV", "GSL_EBADTOL", "GSL_ETOL", "GSL_EUNDRFLW", "GSL_EOVRFLW", "GSL_ELOSS", "GSL_EROUND", "GSL_EBADLEN", "GSL_ENOTSQR", "GSL_ESING", "GSL_EDIVERGE", "GSL_EUNSUP", "GSL_EUNIMPL", "GSL_ECACHE", "GSL_ETABLE", "GSL_ENOPROG", "GSL_ENOPROGJ", "GSL_ETOLF", "GSL_ETOLX", "GSL_ETOLG", "GSL_EOF"};

static int
steffenson_init (void * vstate, gsl_function_fdf * fdf, double * root)
{
  steffenson_state_t * state = (steffenson_state_t *) vstate;

  const double x = *root ;

  state->f = GSL_FN_FDF_EVAL_F (fdf, x);
  state->df = GSL_FN_FDF_EVAL_DF (fdf, x) ;

  state->x = x;
  state->x_1 = 0.0;
  state->x_2 = 0.0;

  state->count = 1;

  return GSL_SUCCESS;

}

static int
steffenson_iterate (void * vstate, gsl_function_fdf * fdf, double * root)
{
  steffenson_state_t * state = (steffenson_state_t *) vstate;
  
  double x_new, f_new, df_new;

  double x_1 = state->x_1 ;
  double x = state->x ;

  if (state->df == 0.0)
    {
      CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"(%s) derivative is zero",gsl_errtxt[GSL_EZERODIV]);
      return GSL_EZERODIV;
    }

  x_new = x - (state->f / state->df);
  
  GSL_FN_FDF_EVAL_F_DF(fdf, x_new, &f_new, &df_new);

  state->x_2 = x_1 ;
  state->x_1 = x ;
  state->x = x_new;

  state->f = f_new ;
  state->df = df_new ;

  if (!isfinite (f_new))
    {
      CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"(%s) function value is not finite",gsl_errtxt[GSL_EBADFUNC]);
      return GSL_EBADFUNC;
    }

  if (state->count < 3)
    {
      *root = x_new ;
      state->count++ ;
    }
  else 
    {
      double u = (x - x_1) ;
      double v = (x_new - 2 * x + x_1);

      if (v == 0)
        *root = x_new;  /* avoid division by zero */
      else
        *root = x_1 - u * u / v ;  /* accelerated value */
    }

  if (!isfinite (df_new))
    {
      CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"(%s) derivative value is not finite",gsl_errtxt[GSL_EBADFUNC]);
      return GSL_EBADFUNC;
    }
      
  return GSL_SUCCESS;
}

int
prims_gsl_root_test_delta (double x1, double x0, double epsabs, double epsrel)
{
  const double tolerance = epsabs + epsrel * fabs(x1)  ;

  if (epsrel < 0.0)
  {
    CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"(%s) relative tolerance is negative",gsl_errtxt[GSL_EBADTOL]);
    return GSL_EBADTOL;
  }
  
  if (epsabs < 0.0)
  {
    CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"(%s) absolute tolerance is negative",gsl_errtxt[GSL_EBADTOL]);
    return GSL_EBADTOL;
  }
  
  if (fabs(x1 - x0) < tolerance || x1 == x0)
    return GSL_SUCCESS;
  
  return GSL_CONTINUE ;
}

int
prims_gsl_root_test_interval (double x_lower, double x_upper, double epsabs, double epsrel)
{
  const double abs_lower = fabs(x_lower) ;
  const double abs_upper = fabs(x_upper) ;

  double min_abs, tolerance;

  if (epsrel < 0.0)
  {
    CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"(%s) relative tolerance is negative",gsl_errtxt[GSL_EBADTOL]);
    return GSL_EBADTOL;
  }
  
  if (epsabs < 0.0)
  {
    CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"(%s) absolute tolerance is negative",gsl_errtxt[GSL_EBADTOL]);
    return GSL_EBADTOL;
  }

/* this interferes with how I use the routine :-)
  if (x_lower > x_upper)
  {
    CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,"(%s) lower bound larger than upper bound",gsl_errtxt[GSL_EINVAL]);
    return GSL_EINVAL;
  }
*/

  if ((x_lower > 0.0 && x_upper > 0.0) || (x_lower < 0.0 && x_upper < 0.0))
    {
      min_abs = abs_lower < abs_upper ? abs_lower : abs_upper;
    }
  else
    {
      min_abs = 0;
    }

  tolerance = epsabs + epsrel * min_abs  ;

  if (fabs(x_upper - x_lower) < tolerance)
    return GSL_SUCCESS;

  return GSL_CONTINUE ;
}

double steffenson_solve(gsl_function_fdf * fdf, double x0, const int max_iter,
  const double epsabs, const double epsrel, int *ierr)
{
  steffenson_state_t state;
  int iter, status;
  double x, xold;
  
  steffenson_init(&state, fdf, &x0);

  iter = 0;xold = x0;
  do
   {
     iter++;
     if((status = steffenson_iterate (&state, fdf, &x)) != GSL_SUCCESS)
       break;
     status = prims_gsl_root_test_delta (x, xold, epsabs, epsrel);
     if(status != GSL_SUCCESS && status != GSL_CONTINUE)
       break;
     xold = x;
   } while (status == GSL_CONTINUE && iter < max_iter);

  if(iter == max_iter) // could not find solution in time, warn
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING, "Steffenson solver could not find solution after %d iterations", max_iter);

  *ierr = status;

  return x;
}
