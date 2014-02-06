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

#include "emf.h"
#include <gsl/gsl_multimin.h>

using namespace std;

#define SIMPLEX_STOP_ITER 5000
#define SIMPLEX_STOP_SIZE 1e-4 

void simplex (ALLDATA &A, double &chisq) 
{
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function m;

  int np=0;
  for (int i=0;i<NP;i++)
    if (A.is[i]) np++;

  /* Starting point */
  x = gsl_vector_alloc (np);
  for (int i=0,j=0;i<NP; i++)
    if (A.is[i]) gsl_vector_set (x, j++, A.p[i]);

  /* Set initial step sizes to 0.1 */
  ss = gsl_vector_alloc (np);
  gsl_vector_set_all (ss, 0.1);

  /* Initialize method and iterate */
  m.n = np;
  m.f = &simplex_f;
  m.params = (void *)&A;

  s = gsl_multimin_fminimizer_alloc (T, np);
  gsl_multimin_fminimizer_set (s, &m, x, ss);

  unsigned int iter = 0;
  double size;
  int status;
  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
    if (status) break;
    size = gsl_multimin_fminimizer_size (s);
    status = gsl_multimin_test_size (size,SIMPLEX_STOP_SIZE);
    if (status == GSL_SUCCESS) { /* reserved */ }
#ifdef DEBUG
      printf ("%5d [",iter);
      for (i=0,j=0;i<NP;i++)
	if (A.is[i]) printf ("%10.7f ",gsl_vector_get (s->x, j++));
      printf("] f= %g sz=%.3f\n",s->fval,size);
#endif
    } while (status == GSL_CONTINUE && iter < SIMPLEX_STOP_ITER);

  chisq = s->fval;

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
}

double simplex_f (const gsl_vector *p,void *par)
{
  ALLDATA *A = (ALLDATA *)par;
  int i,j;
  /* overide parameters by given values */
  for (i=0,j=0;i<NP;i++)
    if (A->is[i]) A->p[i] = gsl_vector_get(p,j++);
  return chi2 (*A, A->func);
}
