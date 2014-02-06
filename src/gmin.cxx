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

using namespace std;

bool global_minimize (ALLDATA &A, int func, double &chisq,int prefered)
{
  int i,dim=0;
  double conv, prev[NP];
  bool CONTINUE = false;
  A.func = func;

  for (i=0;i<NP;i++) {
    prev[i] = A.p[i]; /* for convergence test */
    if ( ((A.attr[i] & P_LOCAL) == 0) && (A.attr[i] & P_ACTIVE) &&
      ((A.attr[i] & P_FIXED) == 0) ) {
      dim++; 	
      A.is[i] = true;
      A.fpar = i;
      }
    else
      A.is[i] = false;
    }

  if (dim == 0) exit(1);

  printf("MINIMIZE> ");

  if (prefered == M_AUTO) {
    if (dim == 1) {
      printf("Brent\n"); 
      brent (A, chisq);
      }
    if (dim >  1) {
      printf("conjugate gradient");
      conjgr (A, chisq);
      //printf("Powell"); 
      //powell (A, chisq);
      //printf("Levenberg-Marquardt"); 
      //levmarw (A, chisq);
      //printf("simulated annealing"); 
      //anneal (A, chisq);
      //printf("downhill simplex"); 
      //simplex (A, chisq);
      //printf("grid search"); 
      //grid (A, chisq);
      }
    }// M_AUTO
  else 
    {
    if (prefered == M_LEVMAR) {
      printf("Levenberg-Marquardt"); 
      levmar (A, chisq);
      }
    if (prefered == M_LEVMARW) {
      printf("Levenberg-Marquardt"); 
      levmarw (A, chisq);
      }
    if (prefered == M_CONJGR) {
      printf("conjugate gradient");
      conjgr (A, chisq);
      }
    if (prefered == M_SIMPLEX) {
      printf("downhill simplex"); 
      simplex (A, chisq);
      }
    if (prefered == M_POWELL) {
      printf("Powell"); 
      powell (A, chisq);
      }
    if (prefered == M_ANNEAL) {
      printf("simulated annealing"); 
      anneal (A, chisq);
      }
    if (prefered == M_BRENT) {
      if (dim > 1) exit(1);
      else {
	printf("Brent"); 
	brent (A, chisq);
	}
      }

    }// prefered

  printf(", target= %s\n",A.fid[func]);
  printf("MINIMIZE> X2 %.2e ",chisq);
  for (int i=0;i<NP;i++) {
    if ((A.attr[i] & P_ACTIVE) && ((A.attr[i] & P_FIXED) == 0) &&
      ((A.attr[i] & P_LOCAL) == 0)) {
      if (i == _phi_ || i == _theta_ ) conv = 180.0/M_PI;
      else conv = 1.0;
      printf("%s %.3f",A.pid[i],conv*A.p[i]);
      if (fabs(A.p[i]-prev[i]) < A.conv[i]) printf("* ");
      else {
	CONTINUE = true;
	printf("  ");
	}
      }
    }
  printf("\n\n");
  if (CONTINUE) return false;
  else return true;
}
