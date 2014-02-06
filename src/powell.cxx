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
#include <gsl/gsl_min.h>

#define POWELL_ITER_MAX		200
#define POWELL_FTOL		(1.0e-6)
#define POWELL_TINY		(1.0e-20)

using namespace std;

/******************************************************/
/* multi-dimensional minimization without derivatives */
/******************************************************/

void powell (ALLDATA &A, double &chisq) {
size_t iter = 0;
int i,j,ibig,status;
double fret,fp,fptt,del,t; 

const int n = NP;
gsl_matrix * xi  = gsl_matrix_calloc (n,n);

double *p   = new double [n];
double *pt  = new double [n];
double *ptt = new double [n];
double *xit = new double [n];

/* initialize directions (unit vectors) */
for (i=0;i<n;i++)
	if (A.is[i])
		gsl_matrix_set (xi,i,i,1.0);

/* save initial point */
for (i=0;i<n;i++) p[i] = pt[i] = A.p[i];

/* evaluate function */
fret = powell_f (p, &A);

/* iteration */
for (iter=0;;++iter) {
	fp = fret;
	ibig = 0;
	del = 0.0;
	for (i=0; i<n; i++) {
		if (A.is[i]) {
			for (j=0; j<n;j++) 
				xit[j] = gsl_matrix_get(xi,j,i);
			fptt = fret;

	 		powell_linmin (xit, &A, fret);

			if (fptt-fret > del) {
				del=fptt-fret;
				ibig=i;
				}
			}
		}

	/* tolerance test */
	if (2.0*(fp-fret) <= POWELL_FTOL*(fabs(fp)+fabs(fret))+POWELL_TINY) {
		gsl_matrix_free (xi);
		delete [] (p);
		delete [] (pt);
		delete [] (ptt);
		delete [] (xit);
		chisq = fret;
		return;
		}

	if (iter == POWELL_ITER_MAX)
		terminate("powell exceeds maximum iterations");

	for (j=0;j<n;j++) {
		ptt[j]=2.0*p[j]-pt[j];
		xit[j]=p[j]-pt[j];
		pt[j]=p[j];
		}
	
	fptt = powell_f(ptt, &A);

	if (fptt < fp) {
		t=2.0*(fp-2.0*fret+fptt)*SQR(fp-fret-del)-del*SQR(fp-fptt);
		if (t < 0.0) {

			powell_linmin (xit, &A, fret);

			for (j=0;j<n;j++) {
				gsl_matrix_set(xi,j,ibig,
					gsl_matrix_get(xi,j,n-1));
				gsl_matrix_set(xi,j,n-1,xit[j]);
				}
			}// t < 0
		}// fptt < fp
	}// iter
}

/****************************************************/
/* powell_linmin: minimization along a given vector */
/****************************************************/

void powell_linmin (double * xit, void * params, double &fret)
{
ALLDATA *A = (ALLDATA *)params;
size_t iter = 0;
int i,status;
double lower,upper,guess,conv;

/* check if xit vector is zero */
for (i=0; i< NP; i++) if (xit[i] != 0.0) break;
if (i == NP) return;

for (i=0; i<NP; i++) {
	A->u[i] = A->p[i];
	A->v[i] = xit[i];
	}

lower = 0.0;
upper = 1.0;
guess = 0.5;
conv  = 1.0E-3;

/** BEGIN Brent routine *******************************************/

/* turn off default GSL error handler (reporting error and abort) */
gsl_set_error_handler_off();

const gsl_min_fminimizer_type *T = gsl_min_fminimizer_brent;
gsl_min_fminimizer *s = gsl_min_fminimizer_alloc (T);
gsl_function F = {&powell_vector, A};

/* adjust boundaries */
	status = gsl_min_fminimizer_set (s, &F, guess, lower, upper);
/*
do {
	status = gsl_min_fminimizer_set (s, &F, guess, lower, upper);
	if (status == GSL_EINVAL) {
		upper *= 1.2;
		lower /= 1.2;
		}
	} while (status == GSL_EINVAL);
*/

do {
	iter++;
	status = gsl_min_fminimizer_iterate (s);
	if (status == GSL_EBADFUNC)
		terminate("Powell/Brent: singular point encountered at Inf or NaN\n");
	if (status == GSL_FAILURE)
		terminate("Powell/Brent: could not improve the current best approximation or bounding interval\n"); 
	fret = gsl_min_fminimizer_f_minimum (s);
	guess = gsl_min_fminimizer_x_minimum (s);
	lower = gsl_min_fminimizer_x_lower (s);
	upper = gsl_min_fminimizer_x_upper (s);
	status = gsl_min_test_interval (lower,upper,conv,0.0);
	}
while (status == GSL_CONTINUE && iter < POWELL_ITER_MAX);

gsl_min_fminimizer_free (s); /* free memory */
gsl_set_error_handler(NULL); /* back to GSL default error handler */

/** END Brent routine *********************************************/

/* update position and direction */
for (i=0;i<NP;i++) {
	xit[i] *= guess;
	A->p[i] =  A->u[i] + xit[i];
	}
}
		
/******************************************************/
/* powell_vector: artificial one-dimensional function */
/******************************************************/

double powell_vector (double x, void * params)
{
  ALLDATA *A = (ALLDATA *)params;
  for (int i=0; i< NP; i++)
    A->p[i] = A->u[i] + x * A->v[i];
  return chi2 (*A, A->func);
}

/*************************************/
/* powell_f: evaluate function value */
/*************************************/

double powell_f (double * p, void * params) {
  ALLDATA *A = (ALLDATA *)params;
  for (int i=0; i<NP; i++) 
    A->p[i] = p[i];
  return chi2 (*A, A->func);
}
