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
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_siman.h>

using namespace std;

/* set up parameters for this simulated annealing run */
/* how many points do we try before stepping */
#define N_TRIES 200
/* how many iterations for each T? */
#define ITERS_FIXED_T 200
/* max step size in random walk */
#define STEP_SIZE 0.1
/* Boltzmann constant */
#define K 1.0
/* initial temperature */
#define T_INITIAL 1.0e-2
/* damping factor for temperature */
#define MU_T 1.003
#define T_MIN 1.0e-6

double E1(void *xp)
{
  ALLDATA *A = (ALLDATA *)xp;
  return chi2 (*A,A->func);
}

/* metric */
double M1(void *xp, void *yp)
{
  ALLDATA *A = (ALLDATA *) xp;
  ALLDATA *B = (ALLDATA *) yp;
  int i;
  double distance;
  for (i=0;i<NP;i++)
    if (A->is[i]) distance += fabs(A->p[i] -B->p[i]);
  return distance;
}

/* step */
void S1(const gsl_rng * r, void *xp, double step_size)
{
  ALLDATA *A = (ALLDATA *) xp;
  int i;
  for (i=0;i<NP;i++)
    if (A->is[i])
      A->p[i] += step_size*(2*gsl_rng_uniform(r) -1);
}

/* print */
void P1(void *xp)
{
  ALLDATA *A = (ALLDATA *) xp;
  int i;
  for (i=0; i<NP; i++)
    if (A->is[i])
      printf (" %s %12g",A->pid[i],A->p[i]);
}

gsl_siman_params_t params = {N_TRIES, ITERS_FIXED_T, STEP_SIZE,
  K, T_INITIAL, MU_T, T_MIN};

void anneal (ALLDATA &A, double &chisq)
{
  extern gsl_rng * rng;

  gsl_siman_solve(rng,(void *)&A,E1,S1,M1,P1,NULL,NULL,NULL,sizeof(A),params);

}
