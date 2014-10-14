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
  multidimensional minimization with derivatives
  using conjugate gradient algorithm
*/

#include "emf.h"
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_deriv.h>

#define CONJGR_STEP_SIZE      	(1e-1)
#define CONJGR_STOP_ITER       	(400)
#define CONJGR_STOP_DELTA      	(1e-4)
#define CONJGR_STOP_GRAD      	(1e-9)
#define CONJGR_FTOL      	(1e-9)
#define N_DERIVATIVE_h		(1e-8)

void conjgr (ALLDATA &A, double &chisq) 
{
  const gsl_multimin_fdfminimizer_type * T;
  gsl_multimin_fdfminimizer * s;
  gsl_vector *x;

  int np=0;
  for (int i=0;i<NP;i++)
    if (A.is[i]) np++;

  /* Starting point */
  x = gsl_vector_alloc (np);
  for (int i=0,j=0;i<NP;i++)
    if (A.is[i]) gsl_vector_set (x,j++,A.p[i]);

  gsl_multimin_function_fdf m;
  m.f	    = &conjgr_f;
  m.df	    = &conjgr_df;   /* numerical derivative: &conjgr_ndf */
  m.fdf	    = &conjgr_fdf;  /* numerical derivative: &conjgr_nfdf */
  m.n	    = np;
  m.params  = (void *)&A;

  /* Fletcher-Reeves */
  T = gsl_multimin_fdfminimizer_conjugate_fr;
  /* Polak-Ribiere */
  //T = gsl_multimin_fdfminimizer_conjugate_pr;
  /* Broyden-Fletcher-Goldfarb-Shanno (BFGS) */
  //T = gsl_multimin_fdfminimizer_vector_bfgs;
  //T = gsl_multimin_fdfminimizer_vector_bfgs2; /* ??? */

  s = gsl_multimin_fdfminimizer_alloc (T, np);
  gsl_multimin_fdfminimizer_set (s,&m,x,CONJGR_STEP_SIZE,CONJGR_FTOL); 

  unsigned int iter = 0;
  double pf=0.0;
  int status;
  do {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate (s);
    if (status) break;
    status = gsl_multimin_test_gradient (s->gradient,CONJGR_STOP_GRAD);
    if (status == GSL_SUCCESS) {/* reserved */ break;}
    if (fabs(pf-s->f) < CONJGR_STOP_DELTA) break;
      else pf = s->f;
    } while (status == GSL_CONTINUE && iter < CONJGR_STOP_ITER);
  chisq = s->f;
  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);
}

double conjgr_f (const gsl_vector *p,void *par) 
{
  ALLDATA *A = (ALLDATA *)par;
  double chisq;
  int nr,r,i,j;

  /* overide parameters by given values */
  for (i=0,j=0;i<NP;i++)
    if (A->is[i]) A->p[i] = gsl_vector_get(p,j++);

  for (nr=0,r=0;r<A->NR;r++) if (A->flag[r]) nr++;

  if (A->func == GX2_R1R2NOE_BM) {
    gsl_vector *c = gsl_vector_alloc(3*A->NF*nr);
    chi_gbm_f (p,par,c);
    chisq = SQR(gsl_blas_dnrm2(c));
    gsl_vector_free(c);
    }

  if (A->func == GX2_R1R2NOE_FM) {
    gsl_vector *c = gsl_vector_alloc(3*A->NF*nr);
    chi_gfm_f (p,par,c);
    chisq = SQR(gsl_blas_dnrm2(c));
    gsl_vector_free(c);
    }

  if (A->func == GX2_R2_OVER_R1) {
    gsl_vector *c = gsl_vector_alloc(A->NF*nr);
    chi_g21_f (p,par,c);
    chisq = SQR(gsl_blas_dnrm2(c));
    gsl_vector_free(c);
    }

  return chisq;
}

void conjgr_df (const gsl_vector *p,void *par,gsl_vector *g) 
{
  ALLDATA *A = (ALLDATA *)par;
  int i,j,r,f,nr,np;
  double gradient;

  /* overide parameters by given values */
  for (i=0,j=0;i<NP;i++)
    if (A->is[i]) A->p[i] = gsl_vector_get(p,j++);

  for (i=0,np=0;i<NP;i++) if (A->is[i]) np++;
  for (nr=0,r=0;r<A->NR;r++) if (A->flag[r]) nr++;

  if (A->func == GX2_R1R2NOE_BM) {
    gsl_vector *c = gsl_vector_alloc(3*A->NF*nr);
    gsl_matrix *J = gsl_matrix_alloc(3*A->NF*nr,np);
    chi_gbm_fdf(p,par,c,J);
    /* convert jacobian matrix to gradient vector */
    for (j=0;j<np;j++) {
      gradient = 0.0;
      for (i=0,r=0;r<A->NR;r++) {
	if (A->flag[r]) {
	  for (f=0;f<A->NF;f++) {
	    gradient += 
	      2*gsl_vector_get(c,3*f*i+0)*gsl_matrix_get(J,3*f*i+0,j)+
	      2*gsl_vector_get(c,3*f*i+1)*gsl_matrix_get(J,3*f*i+1,j)+
	      2*gsl_vector_get(c,3*f*i+2)*gsl_matrix_get(J,3*f*i+2,j);
	    }// f
	  i++;
	  }// if
	}// r
      gsl_vector_set (g,j,gradient);
      }// j
    gsl_vector_free (c);
    gsl_matrix_free (J);	
    }

  if (A->func == GX2_R2_OVER_R1) {
    gsl_vector *c = gsl_vector_alloc(A->NF*nr);
    gsl_matrix *J = gsl_matrix_alloc(A->NF*nr,np);
    chi_g21_fdf(p,par,c,J);
    /* convert jacobian matrix to gradient vector */
    for (j=0;j<np;j++) {
      gradient = 0.0;
      for (i=0,r=0;r<A->NR;r++) {
	if (A->flag[r]) {
	  for (f=0;f<A->NF;f++)
	    gradient += 2*gsl_vector_get(c,f*i)*gsl_matrix_get(J,f*i,j);
	  i++;
	  }// if
	}// r
      gsl_vector_set (g,j,gradient);
      }// j
    gsl_vector_free (c);
    gsl_matrix_free (J);	
    }
}

void conjgr_fdf (const gsl_vector *p,void *par,double *chisq,gsl_vector *g) 
{
  ALLDATA *A = (ALLDATA *)par;
  int i,j,r,f,nr,np;
  double gradient;

  /* overide parameters by given values */
  for (i=0,j=0;i<NP;i++)
    if (A->is[i]) A->p[i] = gsl_vector_get(p,j++);

  for (i=0,np=0;i<NP;i++) if (A->is[i]) np++;
  for (nr=0,r=0;r<A->NR;r++) if (A->flag[r]) nr++;


  if (A->func == GX2_R1R2NOE_BM) {
    gsl_vector *c = gsl_vector_alloc(3*A->NF*nr);
    gsl_matrix *J = gsl_matrix_alloc(3*A->NF*nr,np);
    chi_gbm_fdf(p,par,c,J);
    *chisq = SQR(gsl_blas_dnrm2(c));
    /* convert jacobian matrix to gradient vector */
    for (j=0;j<np;j++) {
      gradient = 0.0;
      for (i=0,r=0;r<A->NR;r++) {
	if (A->flag[r]) {
	  for (f=0;f<A->NF;f++) {
	    gradient += 
	      2*gsl_vector_get(c,3*f*i+0)*gsl_matrix_get(J,3*f*i+0,j)+
	      2*gsl_vector_get(c,3*f*i+1)*gsl_matrix_get(J,3*f*i+1,j)+
	      2*gsl_vector_get(c,3*f*i+2)*gsl_matrix_get(J,3*f*i+2,j);
/*
	    gradient += gsl_matrix_get(J,3*f*i+0,j)+
	      gsl_matrix_get(J,3*f*i+1,j)+
	      gsl_matrix_get(J,3*f*i+2,j);
*/
	    }// f
	  i++;
	  }// if
	}// r
      gsl_vector_set (g,j,gradient);
      }// j
    gsl_vector_free (c);
    gsl_matrix_free (J);	
    }

  if (A->func == GX2_R2_OVER_R1) {
    gsl_vector *c = gsl_vector_alloc(A->NF*nr);
    gsl_matrix *J = gsl_matrix_alloc(A->NF*nr,np);
    chi_g21_fdf(p,par,c,J);
    *chisq = SQR(gsl_blas_dnrm2(c));
    /* convert jacobian matrix to gradient vector */
    for (j=0;j<np;j++) {
      gradient = 0.0;
      for (i=0,r=0;r<A->NR;r++) {
	if (A->flag[r]) {
	  for (f=0;f<A->NF;f++)
	    gradient += 2*gsl_vector_get(c,f*i+0)*gsl_matrix_get(J,f*i+0,j);
/*
	    gradient += gsl_matrix_get(J,f*i+0,j);
*/
	  i++;
	  }// if
	}// r
      gsl_vector_set (g,j,gradient);
      }// j
    gsl_vector_free (c);
    gsl_matrix_free (J);	
    }
}

double conjgr_nf (double x, void *par) 
{
  ALLDATA *A = (ALLDATA *)par;
  A->p[A->fpar] = x;
  return chi2 (*A, A->func);
}

void conjgr_ndf (const gsl_vector *p,void *par,gsl_vector *g) {
  ALLDATA *A = (ALLDATA *)par;
  double fret,abserr;
  int i,j;

  gsl_function F = {&conjgr_nf, par};
  for (i=0,j=0;i<NP;i++)
    if (A->is[i]) {
      A->p[i] = gsl_vector_get (p,j);
      A->fpar = i;
      /* x-h, x-h/2, x, x+h/2, x+h */
      gsl_deriv_central (&F,A->p[i], N_DERIVATIVE_h,&fret,&abserr);
      gsl_vector_set (g,j,fret);
      j++;
      }
}

void conjgr_nfdf (const gsl_vector *p,void *par,double *chisq,gsl_vector *g)
{
  *chisq = conjgr_f (p, par);
  conjgr_ndf (p,par,g);
}
