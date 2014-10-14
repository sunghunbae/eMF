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
  Grid Search

  definition of grid points
  lb ... ub, step
  1st grid point  : lb + step
  last grid point : ub - step
  total number of points: (grds - 1)
  lb and ub are excluded from grid evaluation
*/

#include "emf.h"

using namespace std;

void ax_grid_search (ALLDATA &A) 
{
  double x2,bic,aic,x2_,bic_,aic_;
  size_t r,m,n;
  int p = A.grds[_Dr_]+1;
  int q = A.grds[_tc_]+1;
  gsl_matrix *g = gsl_matrix_alloc (p*q,8);
  gsl_vector_view s;

  printf("******************** BEGIN GLOBAL AXIAL GRID SEARCH ********************\n");
  n=0;
  A.p[_Dr_] = A.lb[_Dr_];
  while (p > 0) {
    A.p[_tc_] = A.lb[_tc_];
    q = A.grds[_tc_]+1;
    while (q > 0) {
      x2 = x2_ = 0.0;
      bic = bic_ = 0.0;
      aic = aic_ = 0.0;
      for (r=0;r<A.NR;r++) { 
	select_all (A);
	if(A.flag[r]) {
	  A.r = r;
	  for (m=1;m<=NM;m++) 
	    fitmodel(m,0,A);
	  select_model(A);
	  x2  += gsl_matrix_get(A.x2, r,A.best[r]);
	  bic += gsl_matrix_get(A.bic,r,A.best[r]);
	  aic += gsl_matrix_get(A.aic,r,A.best[r]);
	  }// flag
	select_cluster(1,A);
	if (A.flag[r]) {
	  x2_  += gsl_matrix_get(A.x2, r,A.best[r]);
	  bic_ += gsl_matrix_get(A.bic,r,A.best[r]);
	  aic_ += gsl_matrix_get(A.aic,r,A.best[r]);
	  }// flag
	}// r
      printf("AXIAL_GRID_ALL> tc= %5.2f Dr= %4.2f ",A.p[_tc_],A.p[_Dr_]);
      printf("X2= %8.3f BIC= %8.3f AIC= %8.3f\n",x2,bic,aic);
      printf("AXIAL_GRID_SUB> tc= %5.2f Dr= %4.2f ",A.p[_tc_],A.p[_Dr_]);
      printf("X2= %8.3f BIC= %8.3f AIC= %8.3f\n",x2_,bic_,aic_);
      gsl_matrix_set (g,n,0,A.p[_tc_]);
      gsl_matrix_set (g,n,1,A.p[_Dr_]);
      gsl_matrix_set (g,n,2,x2);
      gsl_matrix_set (g,n,3,bic);
      gsl_matrix_set (g,n,4,aic);
      gsl_matrix_set (g,n,5,x2_);
      gsl_matrix_set (g,n,6,bic_);
      gsl_matrix_set (g,n,7,aic_);
      n++;
      A.p[_tc_] += A.step[_tc_];
      q--;
      }
    A.p[_Dr_] += A.step[_Dr_];
    p--;
    }
  printf("******************** END GLOBAL AXIAL GRID SEARCH ********************\n");

  /* best X2 */
  s = gsl_matrix_column (g,2); 
  n = gsl_vector_min_index (&s.vector);
  printf("BEST_ALL(X2)    tc= %5.2f Dr= %4.2f X2= %8.3f BIC= %8.3f AIC= %8.3f\n",
    gsl_matrix_get(g,n,0), gsl_matrix_get(g,n,1), gsl_matrix_get(g,n,2),
    gsl_matrix_get(g,n,3), gsl_matrix_get(g,n,4));
  s = gsl_matrix_column (g,5); 
  n = gsl_vector_min_index (&s.vector);
  printf("BEST_SUB(X2)    tc= %5.2f Dr= %4.2f X2= %8.3f BIC= %8.3f AIC= %8.3f\n",
    gsl_matrix_get(g,n,0), gsl_matrix_get(g,n,1), gsl_matrix_get(g,n,5),
    gsl_matrix_get(g,n,6), gsl_matrix_get(g,n,7));

  /* best BIC */
  s = gsl_matrix_column (g,3); 
  n = gsl_vector_min_index (&s.vector);
  printf("BEST_ALL(BIC)   tc= %5.2f Dr= %4.2f X2= %8.3f BIC= %8.3f AIC= %8.3f\n",
    gsl_matrix_get(g,n,0), gsl_matrix_get(g,n,1), gsl_matrix_get(g,n,2),
    gsl_matrix_get(g,n,3), gsl_matrix_get(g,n,4));
  s = gsl_matrix_column (g,6); 
  n = gsl_vector_min_index (&s.vector);
  printf("BEST_SUB(BIC)   tc= %5.2f Dr= %4.2f X2= %8.3f BIC= %8.3f AIC= %8.3f\n",
    gsl_matrix_get(g,n,0), gsl_matrix_get(g,n,1), gsl_matrix_get(g,n,5),
    gsl_matrix_get(g,n,6), gsl_matrix_get(g,n,7));
  A.p[_tc_] = gsl_matrix_get (g,n,0);
  A.p[_Dr_] = gsl_matrix_get (g,n,1);

  /* best AIC */
  s = gsl_matrix_column (g,4); 
  n = gsl_vector_min_index (&s.vector);
  printf("BEST_ALL(AIC)   tc= %5.2f Dr= %4.2f X2= %8.3f BIC= %8.3f AIC= %8.3f\n",
    gsl_matrix_get(g,n,0), gsl_matrix_get(g,n,1), gsl_matrix_get(g,n,2),
    gsl_matrix_get(g,n,3), gsl_matrix_get(g,n,4));
  s = gsl_matrix_column (g,7); 
  n = gsl_vector_min_index (&s.vector);
  printf("BEST_SUB(AIC)   tc= %5.2f Dr= %4.2f X2= %8.3f BIC= %8.3f AIC= %8.3f\n",
    gsl_matrix_get(g,n,0), gsl_matrix_get(g,n,1), gsl_matrix_get(g,n,5),
    gsl_matrix_get(g,n,6), gsl_matrix_get(g,n,7));

  gsl_matrix_free(g);
}

void grid_search (ALLDATA &A, int func, double &chisq, bool output)
{
  double unit;

  grid_index_end = 0;
  grid_min_x2 = -1.0;
  for (int i=0;i<NP;i++) grid_min_par[i]=0.0;
  recursive_s (A, func, 0);
  for (int i=0;i<NP;i++) A.p[i] = grid_min_par[i];
  chisq = grid_min_x2;

  if (output) {
    printf("MINIMIZE> grid search, target= %s\n",A.fid[func]);
    printf("MINIMIZE> X2 %.2e ",chisq);
    for (int i=0;i<NP;i++)
      if (A.is[i]) {
	if (i == _phi_ || i == _theta_) unit = 180.0/M_PI;	
	else unit = 1.0;
	printf("%s %.3f  ",A.pid[i],unit*A.p[i]);
	}
    printf("\n\n");
    }
}

/* recurseve N-dimensional grid search */

void recursive_s (ALLDATA &A, int func, int grid_index)
{

  if (grid_index == 0)
    for(size_t c=0;c<NP;c++) 
      if (A.is[c]) grid_index_end = c;

  if (grid_index == NP) return;
  else { // grid_index = 0 ... (NP-1)
    if (A.is[grid_index]) {
      for (size_t j=1;j<A.grds[grid_index];j++) {
	A.p[grid_index] = A.lb[grid_index] + A.step[grid_index] * j;
	recursive_s (A, func, grid_index+1);
	if (grid_index == grid_index_end) {
	  double chisq = chi2 (A, func);
	  if ((grid_min_x2 == -1.0 || grid_min_x2 > chisq) && chisq >= 0) {
	    grid_min_x2 = chisq;
	    for (size_t k=0;k<NP;k++) 
	      grid_min_par[k] = A.p[k];
	    } // if
	  } // if grid_index == grid_index_end
	} // for j
    return;
    } // if A.is[grid_index]
    else 
      recursive_s (A, func, grid_index+1);
    }; // else grid_index = 0 ... (NP-1)
}
