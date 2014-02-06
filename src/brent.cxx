/*
    eMF Copyright 2008 Sung-Hun Bae

    This file is part of eMF.

    eMF is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    eMF is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with eMF.  If not, see <http://www.gnu.org/licenses/>
*/

/* 
  one dimensional optimization without derivative 
  Brent's Algorithm
*/

#include "emf.h"
#include <gsl/gsl_min.h>
#include <gsl/gsl_blas.h>

#define BRENT_STOP_BRACKET  20
#define BRENT_STOP_ITER	    1000000

using namespace std;

void brent (ALLDATA &A, double &chisq)
{
  /* interval [lower...upper], initial guess */
  extern bool opt_verb;
  double lower,upper,guess,conv;
  int i,status,iter;

  lower=A.lb[A.fpar];
  upper=A.ub[A.fpar];
  guess=A.p[A.fpar];
  conv =A.conv[A.fpar];

  const gsl_min_fminimizer_type *T = gsl_min_fminimizer_brent;
  gsl_min_fminimizer *s = gsl_min_fminimizer_alloc (T);
  gsl_function F = {&brent_f, &A};

  /* turn off default GSL error handler */
  gsl_set_error_handler_off();

  /* bracket the minimum point */
  iter = 0;
  do {
    status = gsl_min_fminimizer_set (s, &F, guess, lower, upper);
    if (status == GSL_EINVAL) {
      iter++;
      upper *= 1.2;
      lower /= 1.2;
      A.lb[A.fpar] = lower;
      A.ub[A.fpar] = upper;
      if (opt_verb) 
	printf ("%5d [%10.7f, %10.7f]\n",iter,lower,upper);
      }
    } while (status == GSL_EINVAL && iter < BRENT_STOP_BRACKET);

  /* main minimization */
  iter = 0;
  do {
    iter++;
    status = gsl_min_fminimizer_iterate (s);
    if (status == GSL_EBADFUNC)
      terminate("Brent: singular point encountered at Inf or NaN\n");
    if (status == GSL_FAILURE)
      terminate("Brent: could not improve the approximation\n"); 
    chisq = gsl_min_fminimizer_f_minimum (s);
    guess = gsl_min_fminimizer_x_minimum (s);
    lower = gsl_min_fminimizer_x_lower (s);
    upper = gsl_min_fminimizer_x_upper (s);
    status = gsl_min_test_interval (lower,upper,conv,0.0);
    if (status == GSL_SUCCESS) 
      A.p[A.fpar] = guess;
    if (opt_verb) 
      printf ("%5d [%10.7f, %10.7f] %10.7f %10.7f f= %g\n",
        iter,lower,upper,guess,upper-lower,chisq);
    } while (status == GSL_CONTINUE && iter < BRENT_STOP_ITER);

  gsl_min_fminimizer_free (s);

  /* back to GSL default error handler */
  gsl_set_error_handler(NULL);
}

double brent_f (double x, void * par)
{
  ALLDATA *A = (ALLDATA *)par;
  A->p[A->fpar] = x;
  return chi2(*A, A->func);
}
