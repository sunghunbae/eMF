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
  Levenberg-Marquardt Nonlinear Least Square Fitting 
  levmar : internal dynamics & diffusion fit

  local fitting:
  void levmar (ALLDATA &A,double &chisq, gsl_vector *e)
    e: error from covariance matrix

  global fitting:
  void levmar (ALLDATA &A,double &chisq)
  
  grand global fitting:
  void levmarw (ALLDATA &A, double &chisq)
*/

#include "emf.h"
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>
#define LEVMAR_STOP_GRAD    1e-12
#define LEVMAR_STOP_DELTA   1e-12
#define LEVMAR_STOP_ITER    1000000

using namespace std;

void static print_gsl_matrix (gsl_matrix *m)
{
  for (int i=0; i<m->size1; i++) {
    for (int j=0; j<m->size2; j++)
      printf("%12.4e ",gsl_matrix_get(m,i,j));
    printf("\n");
  }
}


void levmar (ALLDATA &A,double &chisq, gsl_vector *e)
{
  int i,j,np,status;
  unsigned int iter = 0;
  gsl_vector_set_zero (e);

  for (np=0,i=0;i<NP;i++)
    if (A.is[i]) np++;

  gsl_vector * p = gsl_vector_calloc(np); /* parameters */
  gsl_matrix * covar = gsl_matrix_alloc (np,np); 
  gsl_multifit_function_fdf f;
  gsl_multifit_fdfsolver *s;
  const gsl_multifit_fdfsolver_type *T;

  for (i=0,j=0;i<NP;i++)
    if (A.is[i]) gsl_vector_set(p,j++,A.p[i]);

  f.f   = &chi_R1R2NOE_f;
  f.df  = &chi_R1R2NOE_df;
  f.fdf = &chi_R1R2NOE_fdf;
  f.n   = A.NF*3;
  f.p   = np;
  f.params = &A;

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc (T, A.NF*3, np);
  gsl_multifit_fdfsolver_set (s, &f, p);

  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate (s);
    if (status) break;
    status = gsl_multifit_test_delta (s->dx, s->x, 
      LEVMAR_STOP_GRAD, LEVMAR_STOP_DELTA);
    } while (status == GSL_CONTINUE && iter < LEVMAR_STOP_ITER);
	
  chisq = SQR(gsl_blas_dnrm2(s->f));

  /* fitting error from covariance matrix */
  gsl_multifit_covar (s->J, 0.0, covar);

  /* debug */
  //printf("Jacobian:\n"); print_gsl_matrix(s->J);
  //printf("Covar:\n"); print_gsl_matrix(covar);

  for (j=0,i=0;i<NP;i++)
    if (A.is[i]) {
      gsl_vector_set (e, i, sqrt(gsl_matrix_get(covar,j,j)));
      j++;
      }
  gsl_matrix_free (covar);

  gsl_vector_free (p);
  gsl_multifit_fdfsolver_free(s);
}

void levmar (ALLDATA &A, double &chisq)
{
  int i,j,np,nr,status;
  unsigned int iter = 0;

  for (nr=0,i=0;i<A.NR;i++)
    if (A.flag[i]) nr++;
  for (np=0,i=0;i<NP;i++) 
    if (A.is[i]) np++;

  gsl_vector * p = gsl_vector_calloc(np);
  gsl_multifit_function_fdf f;
  gsl_multifit_fdfsolver *s;
  const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;

  for (i=0,j=0;i<NP;i++)
    if (A.is[i]) gsl_vector_set(p,j++,A.p[i]);

  if (A.func == GX2_R1R2NOE_BM) {
    f.f   = &chi_gbm_f;
    f.df  = &chi_gbm_df;
    f.fdf = &chi_gbm_fdf;
    f.n   = 3*A.NF*nr;
    f.p   = np;
    f.params = &A;
    s = gsl_multifit_fdfsolver_alloc (T, 3*A.NF*nr, np);
    }

  if (A.func == GX2_R2_OVER_R1) {
    f.f   = &chi_g21_f;
    f.df  = &chi_g21_df;
    f.fdf = &chi_g21_fdf;
    f.n   = A.NF*nr;
    f.p   = np;
    f.params = &A;
    s = gsl_multifit_fdfsolver_alloc (T, A.NF*nr, np);
    }

  gsl_multifit_fdfsolver_set (s, &f, p);

  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate (s);
    if (status) break;
    status = gsl_multifit_test_delta (s->dx, s->x,
      LEVMAR_STOP_GRAD,LEVMAR_STOP_DELTA);
    } while (status == GSL_CONTINUE && iter < LEVMAR_STOP_ITER);

  chisq = SQR(gsl_blas_dnrm2(s->f));

  gsl_multifit_fdfsolver_free(s);
  gsl_vector_free (p);
}

void levmarw (ALLDATA &A, double &chisq)
{
  int i,j,r,nr,np,nl,ng,status;
  unsigned int iter = 0;

  /* count global parameters */
  for (ng=0,i=0;i<NP;i++) {
    if ((A.attr[i] & P_ACTIVE) && ((A.attr[i] & P_LOCAL)==0) &&
      ((A.attr[i] & P_FIXED) == 0)) ng++;
    }

  /* count data points and local parameters */
  nr = nl = 0;
  for (r=0;r<A.NR;r++) {
    if (A.flag[r]) {
      nr++;
      nl += set_is_local (A,A.best[r]);
      }
    }

  gsl_vector * p = gsl_vector_calloc(ng+nl);
  gsl_multifit_function_fdf f;
  gsl_multifit_fdfsolver *s;
  const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;

  /* setup initial parameter vector :global */
  j=0;
  for (i=0;i<NP;i++) {
    if ((A.attr[i] & P_ACTIVE) && ((A.attr[i] & P_LOCAL)==0) && 
      ((A.attr[i] & P_FIXED) == 0)) 
      gsl_vector_set(p,j++,A.p[i]);
    }

  /* setup initial parameter vector :local */
  for (r=0;r<A.NR;r++) {
    if (A.flag[r]) {
      set_is_local (A,A.best[r]);
      for (i=0;i<NP;i++)
	if (A.is[i]) 
	  gsl_vector_set(p,j++,gsl_matrix_get(A.vP[i],r,A.best[r]));
      }// flag
    }// r

  if (A.func == GX2_R1R2NOE_BM) {
    f.f   = &chi_w_f;
    f.df  = &chi_w_df;
    f.fdf = &chi_w_fdf;
    f.n   = 3*A.NF*nr;
    f.p   = (ng+nl);
    f.params = &A;
    s = gsl_multifit_fdfsolver_alloc (T, 3*A.NF*nr, (ng+nl));
    }

  gsl_multifit_fdfsolver_set (s, &f, p);

  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate (s);
    if (status) break;
    status = gsl_multifit_test_delta (s->dx, s->x,
      LEVMAR_STOP_GRAD,LEVMAR_STOP_DELTA);
    } while (status == GSL_CONTINUE && iter < LEVMAR_STOP_ITER);

  chisq = SQR(gsl_blas_dnrm2(s->f));

  gsl_multifit_fdfsolver_free(s);
  gsl_vector_free (p);
}
