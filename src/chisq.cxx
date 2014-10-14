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
  GSL functions for Levenberg-Marquardt Non-linear Least Square Fit
  called within levmar and levmar functions

  >>>>> local, {R1,R2,NOE}

  chi_R1R2NOE_f  : chi vector
  chi_R1R2NOE_df : Jacobian matrix
  chi_R1R2NOE_fdf: chi vector & Jacobian matrix

  chi[0] - (R1(p)-R1)/sig_R1 at field 1
  chi[1] - (R2(p)-R2)/sig_R2 at field 1
  chi[2] - (NOE(p)-NOE)/sig-NOE at field 1
  chi[3] - (R1(p)-R1)/sig_R1 at field 2
  chi[4] - (R2(p)-R2)/sig_R2 at field 2
  chi[5] - (NOE(p)-NOE)/sig_NOE at field 2
  ......
  Jacob[0][*] - dR1(p)/dp at field 1
  Jacob[1][*] - dR2(p)/dp at field 1
  Jacob[2][*] - dNOE(p)/dp at field 1
  Jacob[3][*] - dR1(p)/dp at field 2
  Jacob[4][*] - dR2(p)/dp at field 2
  Jacob[5][*] - dNOE(p)/dp at field 2
  ......

  >>>>> global, {R1,R2,NOE}, internal dynamics from best model

  chi_gbm_f  : chi vector
  chi_gbm_df : Jacobian matrix
  chi_gbm_fdf: chi vector & Jacobian matrix

  chi[0] - (R1(p)-R1)/sig_R1 at field 1, resid 1
  chi[1] - (R2(p)-R2)/sig_R2 at field 1, resid 1
  chi[2] - (NOE(p)-NOE)/sig-NOE at field 1, resid 1
  chi[3] - (R1(p)-R1)/sig_R1 at field 2, resid 1
  chi[4] - (R2(p)-R2)/sig_R2 at field 2, resid 1
  chi[5] - (NOE(p)-NOE)/sig_NOE at field 2, resid 1
  chi[6] - (R1(p)-R1)/sig_R1 at field 1, resid 2
  chi[7] - (R2(p)-R2)/sig_R2 at field 1, resid 2
  ......
  Jacob[0][*] - dR1(p)/dp at field 1, resid 1
  Jacob[1][*] - dR2(p)/dp at field 1, resid 1
  Jacob[2][*] - dNOE(p)/dp at field 1, resid 1
  Jacob[3][*] - dR1(p)/dp at field 2, resid 1
  Jacob[4][*] - dR2(p)/dp at field 2, resid 1
  Jacob[5][*] - dNOE(p)/dp at field 2, resid 1
  Jacob[6][*] - dR1(p)/dp at field 1, resid 1
  Jacob[7][*] - dR2(p)/dp at field 1, resid 1
  ......


  >>>>> global, {R2/R1}, internal dyanmics assumed

  chi_g21_f  : chi vector
  chi_g21_df : Jacobian matrix
  chi_g21_fdf: chi vector & Jacobian matrix

  e.g.
  chi[0] - (R2/R1(p)-R2/R1)/sig at field 1, resid 1
  chi[1] - (R2/R1(p)-R2/R1)/sig at field 2, resid 1
  chi[2] - (R2/R1(p)-R2/R1)/sig at field 1, resid 2
  chi[3] - (R2/R1(p)-R2/R1)/sig at field 2, resid 2
  ......
  Jacob[0][*] - d(R2(p)/R1(p))/dp at field 1, resid 1
  Jacob[1][*] - d(R2(p)/R1(p))/dp at field 2, resid 1
  Jacob[2][*] - d(R2(p)/R1(p))/dp at field 1, resid 2
  Jacob[3][*] - d(R2(p)/R1(p))/dp at field 2, resid 2
  ......

*/

#include "emf.h"

int chi_R1R2NOE_f (const gsl_vector *p, void *par, gsl_vector *c)
{
  ALLDATA *A = (ALLDATA *)par;
  double R1,R2,NOE,fv=0.0;
  int i,j,f;

  /* overide parameters by given values */
  for (i=0,j=0;i<NP;i++) 
    if (A->is[i]) A->p[i] = gsl_vector_get(p,j++);

  /* impose square-well restrain on each chi-square */
  fv = boundary_penalty (A);

  gsl_vector * relax = gsl_vector_calloc (A->NF*3);
  R1R2NOE (par,relax,NULL);

  for (f=0;f<A->NF;f++) {
    R1 = gsl_vector_get (relax,3*f+0);
    R2 = gsl_vector_get (relax,3*f+1);
    NOE = gsl_vector_get (relax,3*f+2);
    gsl_vector_set(c,3*f+0, fv+(R1-gsl_matrix_get(A->Y,A->r,3*f+0))
      /gsl_matrix_get(A->SIG,A->r,3*f+0));
    gsl_vector_set(c,3*f+1, fv+(R2-gsl_matrix_get(A->Y,A->r,3*f+1))
      /gsl_matrix_get(A->SIG,A->r,3*f+1));
    gsl_vector_set(c,3*f+2, fv+(NOE-gsl_matrix_get(A->Y,A->r,3*f+2))
      /gsl_matrix_get(A->SIG,A->r,3*f+2));
    } // field
  gsl_vector_free (relax);
  return GSL_SUCCESS;
}


int chi_R1R2NOE_df (const gsl_vector *p, void *par, gsl_matrix *J)
{
  ALLDATA *A = (ALLDATA *)par;
  double dR1j,dR2j,dNOEj;
  int f,h,i,j;

  /* overide parameters by given values */
  for (i=0,j=0;i<NP;i++) 
    if (A->is[i]) A->p[i] = gsl_vector_get(p,j++);

  gsl_vector * relax = gsl_vector_calloc (A->NF*3);
  gsl_matrix * jacob = gsl_matrix_calloc (A->NF*3,NP);
  R1R2NOE (par,relax,jacob);

  for (f=0;f<A->NF; f++) {
    for(i=0,j=0;i<NP;i++) {
      if (A->is[i]) {
	dR1j = gsl_matrix_get(jacob,3*f+0,i);
	dR2j = gsl_matrix_get(jacob,3*f+1,i);
	dNOEj = gsl_matrix_get(jacob,3*f+2,i);
	gsl_matrix_set(J,3*f+0,j,dR1j/gsl_matrix_get(A->SIG,A->r,3*f+0));
	gsl_matrix_set(J,3*f+1,j,dR2j/gsl_matrix_get(A->SIG,A->r,3*f+1));
	gsl_matrix_set(J,3*f+2,j,dNOEj/gsl_matrix_get(A->SIG,A->r,3*f+2));
	j++;
	} // A->is
      } // i
    }// f

  gsl_vector_free (relax);
  gsl_matrix_free (jacob);	
  return GSL_SUCCESS;
}

int chi_R1R2NOE_fdf (const gsl_vector *p,void *par,gsl_vector *c,gsl_matrix *J)
{
  ALLDATA *A = (ALLDATA *)par;
  double R1,R2,NOE,dR1j,dR2j,dNOEj,fv=0.0;
  int i,j,f;

  /* overide parameters by given values */
  for (i=0,j=0;i<NP;i++) 
    if (A->is[i]) A->p[i] = gsl_vector_get(p,j++);

  /* impose square-well restrain on each chi-square */
  fv = boundary_penalty (A);

  gsl_vector * relax = gsl_vector_calloc (A->NF*3);
  gsl_matrix * jacob = gsl_matrix_calloc (A->NF*3,NP);
  R1R2NOE (par,relax,jacob);

  for (f=0;f<A->NF;f++) {
    /* function evaluation */
    R1 = gsl_vector_get (relax,3*f+0);
    R2 = gsl_vector_get (relax,3*f+1);
    NOE = gsl_vector_get (relax,3*f+2);
    gsl_vector_set(c,3*f+0, fv+(R1-gsl_matrix_get(A->Y,A->r,3*f+0))
      /gsl_matrix_get(A->SIG,A->r,3*f+0));
    gsl_vector_set(c,3*f+1, fv+(R2-gsl_matrix_get(A->Y,A->r,3*f+1))
      /gsl_matrix_get(A->SIG,A->r,3*f+1));
    gsl_vector_set(c,3*f+2, fv+(NOE-gsl_matrix_get(A->Y,A->r,3*f+2))
      /gsl_matrix_get(A->SIG,A->r,3*f+2));

    /* derivative */
    for(i=0,j=0;i<NP;i++) {
      if (A->is[i]) {
	dR1j = gsl_matrix_get(jacob,3*f+0,i);
	dR2j = gsl_matrix_get(jacob,3*f+1,i);
	dNOEj = gsl_matrix_get(jacob,3*f+2,i);
	gsl_matrix_set(J,3*f+0,j,dR1j/
	  gsl_matrix_get(A->SIG,A->r,3*f+0));
	gsl_matrix_set(J,3*f+1,j,dR2j/
	  gsl_matrix_get(A->SIG,A->r,3*f+1));
	gsl_matrix_set(J,3*f+2,j,dNOEj/
	  gsl_matrix_get(A->SIG,A->r,3*f+2));
	j++;
	} // A->is
      } // i
    }// field
  gsl_vector_free (relax);
  gsl_matrix_free (jacob);	
  return GSL_SUCCESS;
}

int chi_gbm_f (const gsl_vector *p,void *par,gsl_vector *c) 
{
  ALLDATA *A = (ALLDATA *)par;
  gsl_vector_view chi;
  int r,i;

  for (i=0,r=0;r< A->NR; r++) {
    if (A->flag[r]) {
      A->r = r;
      fitmodel (A->best[r],0,*A);
      chi = gsl_vector_subvector (c,3*A->NF*i,3*A->NF);
      chi_R1R2NOE_f (p,par,&chi.vector);
      i++;
      } // flag
    } // r
  return GSL_SUCCESS;
}
	
int chi_gbm_df (const gsl_vector *p,void *par,gsl_matrix *J)
{
  ALLDATA *A = (ALLDATA *)par;
  gsl_matrix_view j;
  int r,i,np;

  for (np=0,i=0;i<NP;i++) 
    if (A->is[i]) np++;
  for (i=0,r=0;r< A->NR; r++) {
    if (A->flag[r]) {
      A->r = r;
      fitmodel (A->best[r],0,*A);
      j = gsl_matrix_submatrix (J,3*A->NF*i,0,3*A->NF,np);
      chi_R1R2NOE_df (p, par, &j.matrix);
      i++;
      } // flag
    } // r
  return GSL_SUCCESS;
}

int chi_gbm_fdf (const gsl_vector *p,void *par,gsl_vector *c,gsl_matrix *J)
{
  ALLDATA *A = (ALLDATA *)par;
  gsl_vector_view chi;
  gsl_matrix_view jac;
  int r,i,np;

  for (np=0,i=0;i<NP;i++) 
    if (A->is[i]) np++;
  for (i=0,r=0;r< A->NR; r++) {
    if (A->flag[r]) {
      A->r = r;
      fitmodel (A->best[r],0,*A);
      chi = gsl_vector_subvector (c,3*A->NF*i,3*A->NF);
      jac = gsl_matrix_submatrix (J,3*A->NF*i,0,3*A->NF,np);
      chi_R1R2NOE_fdf (p,par,&chi.vector,&jac.matrix);
      i++;
      } // flag
    } // r
  return GSL_SUCCESS;
}

int chi_g21_f (const gsl_vector *p,void *par,gsl_vector *c) 
{
  ALLDATA *A = (ALLDATA *)par;
  double R1,R2,fR1,fR2,dR1,dR2,sig,chi;
  int r,f,i,np;

  /* overide parameters by given values */
  for (np=0,i=0;i<NP;i++)
    if (A->is[i]) A->p[i] = gsl_vector_get(p,np++);

  gsl_vector *relax = gsl_vector_alloc(3*A->NF);

  for (i=0,r=0;r< A->NR; r++) {
    if (A->flag[r]) {
      A->r = r;
      R1R2NOE (A,relax,NULL);
      for (f=0;f<A->NF;f++) {
	fR1 = gsl_vector_get(relax,3*f+0); // fit. R1
	fR2 = gsl_vector_get(relax,3*f+1); // fit. R2
	R1 = gsl_matrix_get(A->Y,r,f*3+0);// exp. R1
	R2 = gsl_matrix_get(A->Y,r,f*3+1);// exp. R2
	dR1 = gsl_matrix_get(A->SIG,r,f*3+0);
	dR2 = gsl_matrix_get(A->SIG,r,f*3+1);
	sig = (R2/R1)*sqrt(SQR(dR1/R1)+SQR(dR2/R2));
	chi = (fR2/fR1-R2/R1)/sig;
	gsl_vector_set(c,A->NF*i+f,chi);
	}// f
      i++;
      } // flag
    } // r
  gsl_vector_free (relax);
  return GSL_SUCCESS;
}

int chi_g21_df (const gsl_vector *p,void *par,gsl_matrix *J)
{
  ALLDATA *A = (ALLDATA *)par;
  double R1,R2,fR1,fR2,dR1,dR2,sig,dy;
  int r,f,i,j,k,np;

  /* overide parameters by given values */
  for (np=0,i=0;i<NP;i++)
    if (A->is[i]) A->p[i] = gsl_vector_get(p,np++);

  gsl_vector *relax = gsl_vector_alloc (3*A->NF);
  gsl_matrix *jacob = gsl_matrix_alloc (3*A->NF,NP);

  for (i=0,r=0;r< A->NR; r++) {
    if (A->flag[r]) {
      A->r = r;
      R1R2NOE (A,relax,jacob);
      for (f=0;f<A->NF;f++) {
	fR1 = gsl_vector_get(relax,3*f+0); // fit. R1
	fR2 = gsl_vector_get(relax,3*f+1); // fit. R2
	R1 = gsl_matrix_get(A->Y,r,3*f+0);// exp. R2
	R2 = gsl_matrix_get(A->Y,r,3*f+1);// exp. R2
	dR1 = gsl_matrix_get(A->SIG,r,3*f+0);
	dR2 = gsl_matrix_get(A->SIG,r,3*f+1);
	sig = (R2/R1)*sqrt(SQR(dR1/R1)+SQR(dR2/R2));
	for (j=0,k=0;k<NP;k++) {
	  if (A->is[k]) {
	    dy = gsl_matrix_get(jacob,3*f+2,k)/R1
	      -R2*gsl_matrix_get(jacob,3*f+1,k)/SQR(R1);
	    gsl_matrix_set(J,A->NF*i+f,j,dy/sig);
	    j++;
	    } // if
	  }//k
	}// f
      i++;
      }// flag
    }// r

  gsl_vector_free (relax);
  gsl_matrix_free (jacob);
  return GSL_SUCCESS;
}

int chi_g21_fdf (const gsl_vector *p,void *par,gsl_vector *c,gsl_matrix *J)
{
  ALLDATA *A = (ALLDATA *)par;
  double R1,R2,fR1,fR2,dR1,dR2,sig,chi,dy;
  int r,f,i,j,k,np;

  /* overide parameters by given values */
  for (np=0,i=0;i<NP;i++)
    if (A->is[i]) A->p[i] = gsl_vector_get(p,np++);

  gsl_vector *relax = gsl_vector_alloc (3*A->NF);
  gsl_matrix *jacob = gsl_matrix_alloc (3*A->NF,NP);

  for (i=0,r=0;r< A->NR; r++) {
    if (A->flag[r]) {
      A->r = r;
      R1R2NOE (A,relax,jacob);
      for (f=0;f<A->NF;f++) {
	fR1 = gsl_vector_get(relax,3*f+0); // fit. R1
	fR2 = gsl_vector_get(relax,3*f+1); // fit. R2
	R1 = gsl_matrix_get(A->Y,r,3*f+0);// exp. R1
	R2 = gsl_matrix_get(A->Y,r,3*f+1);// exp. R2
	dR1 = gsl_matrix_get(A->SIG,r,3*f+0);
	dR2 = gsl_matrix_get(A->SIG,r,3*f+1);
	sig = (R2/R1)*sqrt(SQR(dR1/R1)+SQR(dR2/R2));
	chi = (fR2/fR1-R2/R1)/sig;
	gsl_vector_set(c,A->NF*i+f,chi);
	for (j=0,k=0;k<NP;k++) {
	  if (A->is[k]) {
	    dy = gsl_matrix_get(jacob,3*f+2,k)/R1
	      -R2*gsl_matrix_get(jacob,3*f+1,k)/SQR(R1);
	    gsl_matrix_set(J,A->NF*i+f,j,dy/sig);
	    j++;
	    }// if
	  }// k
	}// f
      i++;
      }// flag
    }// r

  gsl_vector_free (relax);
  gsl_matrix_free (jacob);
  return GSL_SUCCESS;
}

int chi_gfm_f (const gsl_vector *p,void *par,gsl_vector *c) 
{
  ALLDATA *A = (ALLDATA *)par;
  gsl_vector_view chi;
  int r,i,j,m;

  for (i=0,r=0;r< A->NR; r++) {
    if (A->flag[r]) {
      A->r = r;
      for (m=1;m<=NM;m++) fitmodel (m,0,*A);
      select_model(*A);
      for (j=0;j<NP;j++)
	A->p[j] = gsl_matrix_get(A->vP[j],r,A->best[r]);
      chi = gsl_vector_subvector (c,3*A->NF*i,3*A->NF);
      chi_R1R2NOE_f (p,par,&chi.vector);
      i++;
      } // flag
    } // r
  return GSL_SUCCESS;
}

double chi2 (ALLDATA &A, int func)
{
  double chisq,R1,R2,dR1,dR2,sig2,fR1,fR2,dy;
  int i,j,k,f,np;

  chisq = 0.0;
  if (func == LX2_R1R2NOE ) {
    /* for grid search */
    for (np=0,i=0;i<NP;i++) if (A.is[i]) np++;
    gsl_vector * p = gsl_vector_calloc(np);
    gsl_vector * chi = gsl_vector_calloc(A.NF*3);
    for (i=0,j=0;i<NP;i++) 
      if (A.is[i]) gsl_vector_set (p,j++,A.p[i]);
    chi_R1R2NOE_f (p,&A,chi);
    chisq = SQR(gsl_blas_dnrm2(chi));
    gsl_vector_free (p);
    gsl_vector_free (chi);
    }

  if (func == GX2_R2_OVER_R1 ) {
    /* for estimation of diffusion tensor */
    /* global chi-square of R2/R1 */
    gsl_vector * relax = gsl_vector_calloc (A.NF*3);
    double fv;
    for (A.r=0;A.r<A.NR;A.r++) {
      R1R2NOE (&A,relax,NULL);
      for (f=0;f<A.NF;f++) {
	/* only if both R1 and R2 are available */
	if (gsl_matrix_get(A.X, A.r, f*3+0) > 0.0 &&
	  gsl_matrix_get(A.X, A.r, f*3+1) > 0.0 && A.flag[A.r]) {
	  R1 = gsl_matrix_get(A.Y,A.r,f*3+0);
	  R2 = gsl_matrix_get(A.Y,A.r,f*3+1);
	  dR1 = gsl_matrix_get(A.SIG,A.r,f*3+0);
	  dR2 = gsl_matrix_get(A.SIG,A.r,f*3+1);
	  sig2 = SQR(R2/R1)*(SQR(dR1/R1)+SQR(dR2/R2));
	  fR1 = gsl_vector_get(relax,3*f+0); // fitted R1
	  fR2 = gsl_vector_get(relax,3*f+1); // fitted R2
	  dy=R2/R1-fR2/fR1;
	  chisq += dy*dy/sig2;
	  /* impose square-well restrain on each chi-square */
	  chisq += boundary_penalty (&A);
	  } // conditions
	} // f
      } // A.r
    gsl_vector_free (relax);
  }

  if (func == GX2_R2_OVER_R1_SIMPLE ) {
    /* for estimation of diffusion tensor */

    gsl_vector * relax = gsl_vector_calloc (A.NF*3);
    const double c0 = gamma_x/(5.0*gamma_h);
    const double c1 = 7.0*SQR(0.921/0.87);
    const double c2 = 13.0/2.0*SQR(0.955/0.87);
    double fv;
    for (A.r=0;A.r<A.NR;A.r++) {
      R1R2NOE (&A,relax,NULL);
      for (f=0;f<A.NF;f++) {
	/* only if R1, R2, NOE are available */
	if (gsl_matrix_get(A.X, A.r, f*3+0) > 0.0 &&
	    gsl_matrix_get(A.X, A.r, f*3+1) > 0.0 && 
	    gsl_matrix_get(A.X, A.r, f*3+2) > 0.0 && 
            A.flag[A.r]) 
        {
          double MHz = gsl_matrix_get(A.X,A.r,3*f+0);
          double wh = 1.0E-3*2.0*M_PI*MHz; // 1e+9 (rad/s)
          double wx = wh*gamma_x/gamma_h; // 1e+9 (rad/s)
          double jwx,jw0;
          Jw (0.,&A,jw0,NULL);
          Jw (wx,&A,jwx,NULL);

	  double R1 = gsl_matrix_get(A.Y,A.r,f*3+0);
	  double R2 = gsl_matrix_get(A.Y,A.r,f*3+1);
          double NOE = gsl_matrix_get(A.Y,A.r,f*3+2);
	  double dR1 = gsl_matrix_get(A.SIG,A.r,f*3+0);
	  double dR2 = gsl_matrix_get(A.SIG,A.r,f*3+1);
          double dNOE = gsl_matrix_get(A.SIG,A.r,f*3+2);
	  double HF = -c0*(1.0-NOE)*R1;
          double dHF = -c0*((-dNOE)*R1+(1.0-NOE)*dR1);
	  double R1p = R1-c1*HF;
	  double R2p = R2-c2*HF;
          double dR1p = dR1-c1*dHF;
          double dR2p = dR2-c2*dHF;
          double rho = (4.0/3.0)*(R1p)/(2.0*R2p-R1p);
	  double sig2 = SQR(rho)*(SQR((2.0*dR2p-dR1p)/(2.0*R2p-R1p))+SQR(dR1p/R1p));

	  dy=rho-jwx/jw0;
	  chisq += dy*dy/sig2;

	  /* impose square-well restrain on each chi-square */
	  chisq += boundary_penalty (&A);
	  } // conditions
	} // f
      } // A.r
    gsl_vector_free (relax);
    }

  if (func == GX2_R1R2NOE_BM ) {
    /* global chi-square of best model */
    for (A.r=0;A.r<A.NR;A.r++) {
      if (A.flag[A.r] && A.best[A.r]!=0) { 
	fitmodel(A.best[A.r],0,A);
	chisq += gsl_matrix_get(A.x2,A.r,A.best[A.r]);
	}// flag 
      }//A.r
    }

  if (func == GX2_R1R2NOE_FM ) {
    /* global chi-square of floating model */
    for (A.r=0;A.r<A.NR;A.r++) {
      if (A.flag[A.r]) { 
	for (i=1;i<=NM;i++) fitmodel(i,0,A);
	select_model(A);
	chisq += gsl_matrix_get(A.x2,A.r,A.best[A.r]);
	}// flag 
      }//A.r
    }

  if (func == GX2_LIKELIHOOD ) {
    for (np=0,i=0;i<NP;i++) if (A.is[i]) np++;
    gsl_vector * p = gsl_vector_calloc(np);
    gsl_vector *chi = gsl_vector_calloc (A.NF*3);
    for (i=0,j=0;i<NP;i++)
      if (A.is[i]) gsl_vector_set (p,j++,A.p[i]);
    chisq = 1.0;
    for (A.r=0;A.r<A.NR;A.r++) {
      if (A.flag[A.r]) {
	fitmodel(A.best[A.r],0,A);
	A.p[_S2f_] = gsl_matrix_get(A.vP[_S2f_],A.r,A.best[A.r]);
	A.p[_S2s_] = gsl_matrix_get(A.vP[_S2s_],A.r,A.best[A.r]);
	A.p[_te_] = gsl_matrix_get(A.vP[_te_],A.r,A.best[A.r]);
	A.p[_Rex_] = gsl_matrix_get(A.vP[_Rex_],A.r,A.best[A.r]);
	chi_R1R2NOE_f (p, &A, chi);
	for (k=0;k<A.NF*3;k++) {
	  chisq *= exp(-0.5*SQR(gsl_vector_get(chi,k)))
	    /sqrt(2*M_PI)/gsl_matrix_get(A.SIG,A.r,k);
	  } // k
	}// flag
      }// A.r
    gsl_vector_free (p);
    gsl_vector_free (chi);
    chisq *= -1; // for minimization
    }
  return chisq;
}

int chi_w_f (const gsl_vector *p, void *par, gsl_vector *c)
{
  ALLDATA *A = (ALLDATA *)par;
  double R1,R2,NOE,fv;
  int i,j,f,r,q;
  gsl_vector_view chi;
  gsl_vector * relax = gsl_vector_alloc (A->NF*3);
  
  gsl_vector_set_zero(c);

  /* overide parameters :global */
  j=0;
  for (i=0;i<NP;i++) {
    if ((A->attr[i] & P_ACTIVE) && ((A->attr[i] & P_LOCAL)==0) &&
      ((A->attr[i] & P_FIXED) == 0))
      A->p[i] = gsl_vector_get(p,j++);
    }
  /* calculate chi */
  for (q=0,r=0;r<A->NR;r++) {
    if (A->flag[r]) {
      A->r = r;
      set_is_local (*A,A->best[r]);
      /* overide parameters :local */
      for (i=0;i<NP;i++) 
	if (A->is[i]) A->p[i] = gsl_vector_get(p,j++);
      /* impose square-well restrain on each chi-square */
      fv = boundary_penalty(A);
      gsl_vector_set_zero (relax);
      R1R2NOE (par,relax,NULL);
      chi = gsl_vector_subvector (c,3*A->NF*q,3*A->NF);
      for (f=0;f<A->NF;f++) {
	R1 = gsl_vector_get (relax,3*f+0);
	R2 = gsl_vector_get (relax,3*f+1);
	NOE = gsl_vector_get (relax,3*f+2);
	gsl_vector_set(&chi.vector,3*f+0, fv+
	  (R1-gsl_matrix_get(A->Y,r,3*f+0))/gsl_matrix_get(A->SIG,r,3*f+0));
	gsl_vector_set(&chi.vector,3*f+1, fv+
	  (R2-gsl_matrix_get(A->Y,r,3*f+1))/gsl_matrix_get(A->SIG,r,3*f+1));
	gsl_vector_set(&chi.vector,3*f+2, fv+
	  (NOE-gsl_matrix_get(A->Y,r,3*f+2))/gsl_matrix_get(A->SIG,r,3*f+2));
	} // field
      q++;
      }// flag
    }// r
  gsl_vector_free (relax);
  return GSL_SUCCESS;
}

int chi_w_df (const gsl_vector *p, void *par, gsl_matrix *J)
{
  ALLDATA *A = (ALLDATA *)par;
  double dR1j,dR2j,dNOEj,z;
  int i,j,k,f,r,q,np,pos;
  gsl_vector * relax = gsl_vector_alloc (A->NF*3);
  gsl_matrix * jacob = gsl_matrix_alloc (A->NF*3,NP);

  gsl_matrix_set_zero(J);

  /* overide parameters :global */
  j=0;
  for (i=0;i<NP;i++) {
    if ((A->attr[i] & P_ACTIVE) && ((A->attr[i] & P_LOCAL)==0) &&
      ((A->attr[i] & P_FIXED) == 0)) {
      A->p[i] = gsl_vector_get(p,j++);
      }
    }
  
  /* calculate jacobian */
  pos = q = 0;
  for (r=0;r<A->NR;r++) {
    if (A->flag[r]) {
      A->r = r;
      np = set_is_local (*A,A->best[r]);
      /* overide parameters :local */
      for (i=0;i<NP;i++)
        if (A->is[i]) A->p[i] = gsl_vector_get(p,j++);

      gsl_vector_set_zero (relax);
      gsl_matrix_set_zero (jacob);
      R1R2NOE (par,relax,jacob);

      /* jacobian matrix :global contribution */
      for(i=0,j=0;i<NP;i++) {
	if ((A->attr[i] & P_ACTIVE) && ((A->attr[i] & P_LOCAL)==0) &&
	  ((A->attr[i] & P_FIXED) == 0)) {
	  for (f=0;f<A->NF;f++) {
	    dR1j = gsl_matrix_get(jacob,3*f+0,i);
	    dR2j = gsl_matrix_get(jacob,3*f+1,i);
	    dNOEj = gsl_matrix_get(jacob,3*f+2,i);
	    z = gsl_matrix_get(J,3*f*q+0,j);
	    gsl_matrix_set(J,3*f*q+0,j,z+dR1j/gsl_matrix_get(A->SIG,r,3*f+0));
	    z = gsl_matrix_get(J,3*f*q+1,j);
	    gsl_matrix_set(J,3*f*q+1,j,z+dR2j/gsl_matrix_get(A->SIG,r,3*f+1));
	    z = gsl_matrix_get(J,3*f*q+2,j);
	    gsl_matrix_set(J,3*f*q+2,j,z+dNOEj/gsl_matrix_get(A->SIG,r,3*f+2));
	    }// f
	  j++;
	  }// global
	}// i

      /* jacobian matrix :local */
      for(i=0,k=j+pos;i<NP;i++) {
	if (A->is[i]) {
	  for (f=0;f<A->NF;f++) {
	    dR1j = gsl_matrix_get(jacob,3*f+0,i);
	    dR2j = gsl_matrix_get(jacob,3*f+1,i);
	    dNOEj = gsl_matrix_get(jacob,3*f+2,i);
	    gsl_matrix_set(J,3*f*q+0,k,dR1j/gsl_matrix_get(A->SIG,r,3*f+0));
	    gsl_matrix_set(J,3*f*q+1,k,dR2j/gsl_matrix_get(A->SIG,r,3*f+1));
	    gsl_matrix_set(J,3*f*q+2,k,dNOEj/gsl_matrix_get(A->SIG,r,3*f+2));
	    } // f
	  k++;
	  } // A->is
	} // i

      pos += np;
      q++;
      }// flag
    }// r
  gsl_vector_free (relax);
  gsl_matrix_free (jacob);
  return GSL_SUCCESS;
}

int chi_w_fdf (const gsl_vector *p, void *par, gsl_vector *c,gsl_matrix *J)
{
  ALLDATA *A = (ALLDATA *)par;
  double R1,R2,NOE,dR1j,dR2j,dNOEj,fv,z;
  int i,j,k,f,r,q,np,pos;
  gsl_vector_view chi;
  gsl_vector * relax = gsl_vector_alloc (A->NF*3);
  gsl_matrix * jacob = gsl_matrix_alloc (A->NF*3,NP);

  gsl_vector_set_zero(c);
  gsl_matrix_set_zero(J);

  /* overide parameters :global */
  j=0;
  for (i=0;i<NP;i++) {
    if ((A->attr[i] & P_ACTIVE) && ((A->attr[i] & P_LOCAL)==0) &&
      ((A->attr[i] & P_FIXED) == 0)) {
      A->p[i] = gsl_vector_get(p,j++);
      }
    }
  
  /* calculate chi and jacobian */
  pos = q = 0;
  for (r=0;r<A->NR;r++) {
    if (A->flag[r]) {
      A->r = r;
      np = set_is_local (*A,A->best[r]);
      /* overide parameters :local */
      for (i=0;i<NP;i++)
        if (A->is[i]) A->p[i] = gsl_vector_get(p,j++);
      /* impose square-well restrain on each chi-square */
      fv = boundary_penalty(A);
      gsl_vector_set_zero (relax);
      gsl_matrix_set_zero (jacob);
      R1R2NOE (par,relax,jacob);

      chi = gsl_vector_subvector (c,3*A->NF*q,3*A->NF);

      /* chi vector */
      for (f=0;f<A->NF;f++) {
        R1 = gsl_vector_get (relax,3*f+0);
        R2 = gsl_vector_get (relax,3*f+1);
        NOE = gsl_vector_get (relax,3*f+2);
        gsl_vector_set(&chi.vector,3*f+0, fv+
          (R1-gsl_matrix_get(A->Y,r,3*f+0))/gsl_matrix_get(A->SIG,r,3*f+0));
        gsl_vector_set(&chi.vector,3*f+1, fv+
          (R2-gsl_matrix_get(A->Y,r,3*f+1))/gsl_matrix_get(A->SIG,r,3*f+1));
        gsl_vector_set(&chi.vector,3*f+2, fv+
          (NOE-gsl_matrix_get(A->Y,r,3*f+2))/gsl_matrix_get(A->SIG,r,3*f+2));
        } // field

      /* jacobian matrix :global contribution */
      for(i=0,j=0;i<NP;i++) {
	if ((A->attr[i] & P_ACTIVE) && ((A->attr[i] & P_LOCAL)==0) &&
	  ((A->attr[i] & P_FIXED) == 0)) {
	  for (f=0;f<A->NF;f++) {
	    dR1j = gsl_matrix_get(jacob,3*f+0,i);
	    dR2j = gsl_matrix_get(jacob,3*f+1,i);
	    dNOEj = gsl_matrix_get(jacob,3*f+2,i);
	    z = gsl_matrix_get(J,3*f*q+0,j);
	    gsl_matrix_set(J,3*f*q+0,j,z+dR1j/gsl_matrix_get(A->SIG,r,3*f+0));
	    z = gsl_matrix_get(J,3*f*q+1,j);
	    gsl_matrix_set(J,3*f*q+1,j,z+dR2j/gsl_matrix_get(A->SIG,r,3*f+1));
	    z = gsl_matrix_get(J,3*f*q+2,j);
	    gsl_matrix_set(J,3*f*q+2,j,z+dNOEj/gsl_matrix_get(A->SIG,r,3*f+2));
	    }// f
	  j++;
	  }// global
	}// i

      /* jacobian matrix :local */
      for(i=0,k=j+pos;i<NP;i++) {
	if (A->is[i]) {
	  for (f=0;f<A->NF;f++) {
	    dR1j = gsl_matrix_get(jacob,3*f+0,i);
	    dR2j = gsl_matrix_get(jacob,3*f+1,i);
	    dNOEj = gsl_matrix_get(jacob,3*f+2,i);
	    gsl_matrix_set(J,3*f*q+0,k,dR1j/gsl_matrix_get(A->SIG,r,3*f+0));
	    gsl_matrix_set(J,3*f*q+1,k,dR2j/gsl_matrix_get(A->SIG,r,3*f+1));
	    gsl_matrix_set(J,3*f*q+2,k,dNOEj/gsl_matrix_get(A->SIG,r,3*f+2));
	    } // f
	  k++;
	  } // A->is
	} // i

      pos += np;
      q++;
      }// flag
    }// r
  gsl_vector_free (relax);
  gsl_matrix_free (jacob);
  return GSL_SUCCESS;
}

double boundary_penalty (void *par)
{
  ALLDATA *A = (ALLDATA *)par;
  double fv=0.0;
  int i;
  /* impose square-well restrain on each chi-square */
  for (i=0;i<NP;i++) {
    if (A->is[i] && A->p[i] < A->lb[i])
      fv += A->lk[i]*SQR(A->p[i] - A->lb[i]);
    if (A->is[i] && A->p[i] > A->ub[i])
      fv += A->uk[i]*SQR(A->p[i] - A->ub[i]);
    }
  if ((D & _BIMODAL_) && (A->p[_tc_] > A->p[_tb_]))
    fv += 1.0e+6;
  return fv;
}
