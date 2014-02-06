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

void initialize (ALLDATA &A,int NM,int NF,int NR,int NK,int MC)
{
  int i;
  A.X	= gsl_matrix_calloc (NR,NF*3);
  A.Y	= gsl_matrix_calloc (NR,NF*3);
  A.SIG = gsl_matrix_calloc (NR,NF*3);
  A.TAU = gsl_matrix_calloc (NR,NK);
  A.WT  = gsl_matrix_calloc (NR,NK);
  A.dof = gsl_matrix_calloc (NR,NM+1);
  A.aic = gsl_matrix_calloc (NR,NM+1);
  A.bic = gsl_matrix_calloc (NR,NM+1);
  A.x2  = gsl_matrix_calloc (NR,NM+1);
  gsl_matrix_set_all (A.SIG,1);
  gsl_matrix_set_all (A.x2,-1);
  for (int i=0;i<NP;i++) 
    {
    A.vP[i] = gsl_matrix_calloc (NR,NM+1);
    A.aP[i] = gsl_matrix_calloc (NR,NM+1);
    A.dP[i] = gsl_matrix_alloc (NR,NM+1);
    gsl_matrix_set_all (A.dP[i],-1.0);
    }

  /* parameter */
  A.is  	= new bool [NP]; 
  A.p   	= new double [NP];
  A.lb  	= new double [NP];
  A.lk  	= new double [NP];
  A.ub  	= new double [NP];
  A.uk  	= new double [NP];
  A.grds 	= new int [NP];
  A.attr	= new char [NP];
  A.step 	= new double [NP];
  A.conv 	= new double [NP];
  A.u   	= new double [NP];
  A.v   	= new double [NP];
  /* anisotropic diffusion */
  A.t	= new double [5];
  A.x 	= new double [NR];
  A.y 	= new double [NR];
  A.z 	= new double [NR];
  A.ivec= new bool   [NR];
  A.alpha = new double [NR];
  /* info */
  A.num = new int [NR];
  A.clst  = new int [NR];
  A.best  = new int [NR];
  A.flag  = new bool [NR];
  A.NF 	= NF; 
  A.NR 	= NR; 
  A.NK 	= NK;

  /* Monte Carlso simulation X2 percentile (internal use) */
  A.percentile = new double [20];

  for (i=0;i<NR;i++) {
    A.num[i] = A.clst[i] = A.best[i] = 0;
    A.alpha[i] = A.x[i] = A.y[i] = A.z[i] = 0.0;
    A.flag[i] = A.ivec[i] = false;
    }
  for (i=0;i<NP;i++) {
    A.is[i] = false;
    A.grds[i] = 0;
    A.p[i]=A.lb[i]=A.ub[i]=A.lk[i]=A.uk[i]= 0.0;
    A.step[i]=1.0e-3;
    A.conv[i]=1.0e-5;
    A.attr[i] = P_RESET;
    }

  /* default diffusion parameter */
  A.p[_scf_] = 1.0;
  for (i=0;i<5;i++) A.t[i] = 0.0;

  strcpy(A.pid[_S2s_],"S2s");
  strcpy(A.pid[_S2f_],"S2f");
  strcpy(A.pid[_te_],"te");
  strcpy(A.pid[_Rex_],"Rex");
  strcpy(A.pid[_tc_],"tc");
  strcpy(A.pid[_tb_],"tb");
  strcpy(A.pid[_c_],"c");
  strcpy(A.pid[_scf_],"scf");
  strcpy(A.pid[_Dxx_],"Dxx");
  strcpy(A.pid[_Dyy_],"Dyy");
  strcpy(A.pid[_Dzz_],"Dzz");
  strcpy(A.pid[_Dr_],"Dr");
  strcpy(A.pid[_phi_],"phi");
  strcpy(A.pid[_theta_],"theta");
  strcpy(A.pid[_psi_],"psi");
  strcpy(A.fid[LX2_R1R2NOE],"local X2 of {R1,R2,NOE}");
  strcpy(A.fid[GX2_R2_OVER_R1],"global X2 of R2/R1");
  strcpy(A.fid[GX2_R1R2NOE_FM],"global X2 of {R1,R2,NOE} at float model");
  strcpy(A.fid[GX2_R1R2NOE_BM],"global X2 of {R1,R2,NOE} at best model");
  strcpy(A.fid[GX2_LIKELIHOOD],"global X2 of likelihood");
}

void freeing (ALLDATA &A)
{
  /* free memory */	
  gsl_matrix_free(A.X);
  gsl_matrix_free(A.Y);
  gsl_matrix_free(A.SIG);
  gsl_matrix_free(A.TAU);
  gsl_matrix_free(A.WT);
  gsl_matrix_free(A.dof);
  gsl_matrix_free(A.aic);
  gsl_matrix_free(A.bic);
  gsl_matrix_free(A.x2);
  for (int i=0;i<NP;i++) 
    {
    gsl_matrix_free(A.vP[i]);
    gsl_matrix_free(A.dP[i]);
    gsl_matrix_free(A.aP[i]);
    }
delete [] (A.is);
delete [] (A.p);
delete [] (A.lb);
delete [] (A.ub);
delete [] (A.grds);
delete [] (A.step);
delete [] (A.conv);
delete [] (A.attr);
delete [] (A.u);
delete [] (A.v);
delete [] (A.x);
delete [] (A.y);
delete [] (A.z);
delete [] (A.t);
delete [] (A.alpha);
delete [] (A.num);
delete [] (A.clst);
delete [] (A.best);
delete [] (A.flag);
delete [] (A.percentile);
}

void select_all (ALLDATA &A)
{
  const extern int D;
  int r;

  if (D == _GLOBAL_AXIAL_ || D == _GLOBAL_ANISOTROPIC_) {
    for (r=0;r<A.NR;r++)
      if (A.ivec[r]) 
	A.flag[r] = true;
    }
  else {
    for (r=0;r<A.NR;r++) 
      A.flag[r] = true;
    }
}

void select_window (ALLDATA &A,int pos,int size)
{
  const extern int D;
  int r;
  int start = ((pos-size) < 0 ? 0: pos-size);
  int end   = ((pos+size) < A.NR ? pos+size: A.NR-1);

  /* reset */
  for (r=0;r<A.NR;r++) A.flag[r] = false;

  printf(_INFO_ "selecting residues from %d to %d\n",A.num[start],A.num[end]);
  if (D == _GLOBAL_AXIAL_ || D == _GLOBAL_ANISOTROPIC_) {
    for (r=start;r<=end;r++) if (A.ivec[r]) A.flag[r] = true;
    }
  else {
    for (r=start;r<=end;r++) A.flag[r] = true;
    }
}


void select_cluster (int n, ALLDATA &A)
{
  for (int r=0;r<A.NR;r++)
    if(A.clst[r]==n) A.flag[r] = true;
    else A.flag[r] = false;
}

void select_optimizer (ALLDATA &A, int cluster, double s2)
{
  double S2;
  int r,b;

  /* 1. select cluster */
  if (cluster != -1) 
    for (r=0;r<A.NR;r++) {
      if (A.clst[r] == cluster) A.flag[r] = true;
      else A.flag[r] = false;
      }

  /* 2. select residues whose best model has order parameter larger than s2 */
  if (s2 > 0.0) {
    if (cluster != -1) {
      for (r=0;r<A.NR;r++) {
	if (A.flag[r]) {
	  b = A.best[r];
	  if (b!=0) {
	    S2 = gsl_matrix_get (A.vP[_S2f_],r,b)
	      *gsl_matrix_get(A.vP[_S2s_],r,b);
	    if (S2 < s2) A.flag[r] = false;
	    }
	  else  
	    A.flag[r] = false;
	  }// if flag
	}// r
      }// cluster
    else {
      for (r=0;r<A.NR;r++) {
	b = A.best[r];
	if (b!=0) {
	  S2 = gsl_matrix_get (A.vP[_S2f_],r,b)
	    *gsl_matrix_get(A.vP[_S2s_],r,b);
	  if (S2 < s2) A.flag[r] = false;
	  }
	else  
	  A.flag[r] = false;
	}//r
      }// else
    }// s2 != 0
}

void select_estimator (ALLDATA &A, int cluster, double noe, double r2stdev)
{
  double S2[A.NF],S1[A.NF],SS2[A.NF],SS1[A.NF];
  int r,f,n;

  /* 1. select cluster */
  if (cluster != -1) { 
    for (r=0;r<A.NR;r++) {
      if (A.clst[r] == cluster) 
	A.flag[r] *= true;
      else 
	A.flag[r] *= false;
      }// r
    }// if cluster is given
}


void setup_attribute (ALLDATA &A) 
{
  const extern int D;

  A.attr[_S2s_] |= (P_ACTIVE | P_LOCAL);
  A.attr[_S2f_] |= (P_ACTIVE | P_LOCAL);
  A.attr[_te_]  |= (P_ACTIVE | P_LOCAL);
  A.attr[_Rex_] |= (P_ACTIVE | P_LOCAL);

  if (D == _GLOBAL_ISOTROPIC_ )
    A.attr[_tc_] |= ((P_ACTIVE | P_ROTDIF) & ~P_LOCAL);
  if (D == _GLOBAL_AXIAL_) {
    A.attr[_tc_]    |= ((P_ACTIVE | P_ROTDIF) & (~P_LOCAL));
    A.attr[_Dr_]    |= ((P_ACTIVE | P_ROTDIF) & (~P_LOCAL));
    A.attr[_phi_]   |= ((P_ACTIVE | P_ROTDIF) & (~P_LOCAL));
    A.attr[_theta_] |= ((P_ACTIVE | P_ROTDIF) & (~P_LOCAL));
    }

  if (D == _WINDOW_ISOTROPIC_ )
    A.attr[_tc_] |= ((P_ACTIVE | P_ROTDIF) & ~P_LOCAL);

  if (D == _LOCAL_ISOTROPIC_ )
    A.attr[_tc_] |= (P_ACTIVE | P_ROTDIF | P_LOCAL);

  if (D == _GLOBAL_ANISOTROPIC_) {
    A.attr[_Dxx_]   |= ((P_ACTIVE | P_ROTDIF) & (~P_LOCAL));
    A.attr[_Dyy_]   |= ((P_ACTIVE | P_ROTDIF) & (~P_LOCAL));
    A.attr[_Dzz_]   |= ((P_ACTIVE | P_ROTDIF) & (~P_LOCAL));
    A.attr[_phi_]   |= ((P_ACTIVE | P_ROTDIF) & (~P_LOCAL));
    A.attr[_theta_] |= ((P_ACTIVE | P_ROTDIF) & (~P_LOCAL));
    A.attr[_psi_]   |= ((P_ACTIVE | P_ROTDIF) & (~P_LOCAL));
    }

  if (D == _GLOBAL_DISTRIBUTION_ )
    A.attr[_scf_] |= ((P_ACTIVE | P_ROTDIF) & (~P_LOCAL));

  if (D == _GLOBAL_BIMODAL_) {
    A.attr[_tc_] |= ((P_ACTIVE | P_ROTDIF) & ~P_LOCAL);
    A.attr[_tb_] |= (P_ACTIVE | P_ROTDIF | P_LOCAL);
    A.attr[_c_]  |= (P_ACTIVE | P_ROTDIF | P_LOCAL);
/*
    A.attr[_tb_] |= ((P_ACTIVE | P_ROTDIF) & ~P_LOCAL);
    A.attr[_c_]  |= ((P_ACTIVE | P_ROTDIF) & ~P_LOCAL);
*/
    }

  if (D == _LOCAL_BIMODAL_) {
    A.attr[_tc_] |= (P_ACTIVE | P_ROTDIF | P_LOCAL);
    A.attr[_tb_] |= (P_ACTIVE | P_ROTDIF | P_LOCAL);
    A.attr[_c_]  |= (P_ACTIVE | P_ROTDIF | P_LOCAL);
    }
  if (D == _LOCAL_COLE_COLE_ || D == _LOCAL_LORENTZIAN_) {
    A.attr[_tc_] |= (P_ACTIVE | P_ROTDIF | P_LOCAL);
    A.attr[_c_]  |= (P_ACTIVE | P_ROTDIF | P_LOCAL);
    }
}

bool diffusion_converged (ALLDATA &A,double *prev)
{
  for (int i=0;i<NP;i++) {
    if ((A.attr[i] & P_ACTIVE) && (A.attr[i] & P_LOCAL) == 0) {
      printf("CONVERGENCE>");
      printf("%s= %.2e (conv. limit = %.2e)\n",
	A.pid[i],A.p[i]-prev[i],A.conv[i]);
      if (fabs(A.p[i]-prev[i]) > A.conv[i]) return false;
      }
    }
  return true;
}
