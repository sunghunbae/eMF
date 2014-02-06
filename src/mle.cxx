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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

/* Monte Carlo Integration Type */
#define MC_PLAIN                0
#define MC_MISER                1
#define MC_VEGAS                2

/* Maximum Likelihood Estimation */
#define MLE_S2f                 0
#define MLE_S2s                 1
#define MLE_te                  2
#define MLE_Rex                 3
#define MLE_DIM                 4

using namespace std;

/***
Monte Carlo Simulation of Parameter Space using Metropolis Algorithm 
Michael Andrec and James H. Prestegard, JMR 130, 217-232 (1998)
Michael Andrec, Gaetano T. Montelione, and Ronald M. Levy, JMR 139, 408-421 (1999)
***/

void marginal_density (size_t calls,ALLDATA &A,double &res,double &err) {
gsl_monte_function G = { &likelihood, MLE_DIM, &A };
extern gsl_rng * rng;

double xl[MLE_DIM],xu[MLE_DIM];

/* set integration limit */
xl[MLE_S2f] = A.lb[_S2f_];
xl[MLE_S2s] = A.lb[_S2s_];
xl[MLE_te]  = A.lb[_te_];
xl[MLE_Rex] = A.lb[_Rex_];
xu[MLE_S2f] = A.ub[_S2f_];
xu[MLE_S2s] = A.ub[_S2s_];
xu[MLE_te]  = A.ub[_te_];
xu[MLE_Rex] = A.ub[_Rex_];

/* PLAIN */
/*
{
gsl_monte_plain_state *s = gsl_monte_plain_alloc (MLE_DIM);
gsl_monte_plain_integrate (&G,xl,xu,MLE_DIM,calls,rng,s,&res,&err);
gsl_monte_plain_free (s);
}
*/

/* MISER */
/*
{
gsl_monte_miser_state *s = gsl_monte_miser_alloc (MLE_DIM);
gsl_monte_miser_integrate (&G,xl,xu,MLE_DIM,calls,rng,s,&res,&err);
gsl_monte_miser_free (s);
}
*/

/* VEGAS */
/*
{
gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (MLE_DIM);
gsl_monte_vegas_integrate (&G,xl,xu,MLE_DIM,calls/5,rng,s,&res,&err);
do {
	gsl_monte_vegas_integrate (&G,xl,xu,MLE_DIM,calls/5,rng,s,&res,&err);
	} while (fabs (s->chisq - 1.0) > 0.5);
gsl_monte_vegas_free (s);
}
*/

/* Metropolis Algorithm */
{
double Q[MLE_DIM][calls];
double p[MLE_DIM],q[MLE_DIM];
double covar=0.1,U,prob_p,prob_q;
int i,j;
res = 0.0;
err = 0.0;// not estimated

// pre-run
for (j=0;j<MLE_DIM;j++)
	Q[j][0] = p[j] = gsl_ran_flat (rng,A.lb[j],A.ub[j]);
for (i=0;i<calls/5;i++) {
        /* Monte Carlo */
        for (j=0;j<MLE_DIM;j++) {
		q[j]=Q[j][i];
                p[j]=q[j]+gsl_ran_gaussian(rng,covar*(A.ub[j]-A.lb[j]));
		}
	prob_p = likelihood (p,MLE_DIM,&A);
	prob_q = likelihood (q,MLE_DIM,&A);
	U = gsl_ran_flat (rng,0,1);
	if (U <= prob_p/prob_q)
		for (j=0;j<MLE_DIM;j++) Q[j][i+1]=p[j];
	else for (j=0;j<MLE_DIM;j++) Q[j][i+1]=q[j];
        }// calls

for (j=0;j<MLE_DIM;j++) Q[j][0] = Q[j][i];

for (i=0;i<calls;i++) {
        /* Monte Carlo */
        for (j=0;j<MLE_DIM;j++) {
		q[j]=Q[j][i];
                p[j]=q[j]+gsl_ran_gaussian(rng,covar*(A.ub[j]-A.lb[j]));
		}
	prob_p = likelihood (p,MLE_DIM,&A);
	prob_q = likelihood (q,MLE_DIM,&A);
	U = gsl_ran_flat (rng,0,1);
	if (U <= prob_p/prob_q) {
		for (j=0;j<MLE_DIM;j++) Q[j][i+1]=p[j];
		res += prob_p;
		}
	else {
		for (j=0;j<MLE_DIM;j++) Q[j][i+1]=q[j];
		res += prob_q;
		}
        }// calls
}

}

/* Bayesian Statistics - calculate likelihood */
double likelihood (double *x, size_t dim, void *data) {
ALLDATA *A = (ALLDATA *)data;
gsl_vector * chi = gsl_vector_calloc(A->NF*3);
double prob = 1.0;
int k;
MLE_chi_R1R2NOE (x,A,chi);
for (k=0;k<(A->NF)*3;k++) {
	prob *= exp(-0.5*SQR(gsl_vector_get(chi,k)))
		/sqrt(2*M_PI)/gsl_matrix_get(A->SIG,A->r,k);
	}
gsl_vector_free (chi);
return prob;
}

/* Evaluate R1,R2,NOE and return chi-square */
void MLE_chi_R1R2NOE (double *p, void *data, gsl_vector *chi)
{
const extern double gamma_h,gamma_x,r_xh,csa_x;
ALLDATA *A = (ALLDATA *)data;
double wh,wx,d,fd_4,fd_8,c,fc,fc_6,fg;
double y0,y1,y2,y3,y4;
double MHz,R1,R2,NOE,Rex,fRex;
int i,j,f;

// field independent constants
d = 1.0E+3*gamma_h*gamma_x*mu0_h_8pi2/(r_xh*r_xh*r_xh);
fd_4 = d*d/4.0;
fd_8 = d*d/8.0;
fg = (gamma_h/gamma_x)*d*d/4.0;
for (f=0;f<A->NF;f++) {
	/* common for R1,R2,NOE if MHz is same */
	MHz = gsl_matrix_get(A->X,A->r,3*f+0); // unit E+6 (Hz or 1/s)
	wh= 1.0E-3*2.0*M_PI*MHz; // unit E+9 (rad/s)
	wx = wh*gamma_x/gamma_h; // unit E+9 (rad/s)
	fc = 1.0E+6*wx*wx*csa_x*csa_x/3.0; 
	MLE_Jw (0.0,  p,data,y0);
	MLE_Jw (wx,   p,data,y1);
	MLE_Jw (wh-wx,p,data,y2);
	MLE_Jw (wh,   p,data,y3);
	MLE_Jw (wh+wx,p,data,y4);

	/* R1 (1/s) */
	R1 = fd_4*(y2+3.0*y1+6.0*y4)+fc*y1;
	gsl_vector_set (chi,3*f+0, (R1 - gsl_matrix_get(A->Y,A->r,3*f+0)) 
		/ gsl_matrix_get (A->SIG,A->r,3*f+0) );

	/* R2 (1/s) */
	if (MHz != gsl_matrix_get(A->X,A->r,3*f+1)) {
		MHz = gsl_matrix_get(A->X,A->r,3*f+1);
		wh= 1.0E-3*2.0*M_PI*MHz; // unit E+9 (rad/s)
		wx = wh*gamma_x/gamma_h;
		c = 1.0E+3*wx*csa_x/sqrt(3.0);
		MLE_Jw (0.0,  p,data,y0);
		MLE_Jw (wx,   p,data,y1);
		MLE_Jw (wh-wx,p,data,y2);
		MLE_Jw (wh,   p,data,y3);
		MLE_Jw (wh+wx,p,data,y4);
		} // skip Jw evaluation if MHz is same

	fc_6 = fc/6.0;

	// quadratically scaled relative to the first R2 field
	fRex = SQR(MHz/gsl_matrix_get(A->X,A->r,1));

	R2 = fd_8*(4.0*y0+y2+3.0*y1+6.0*y3+6.0*y4)+fc_6*(4.0*y0+3.0*y1)
        	+fRex*p[MLE_Rex]; //Rex is quadratically scaled

	gsl_vector_set (chi,3*f+1, (R2 - gsl_matrix_get(A->Y,A->r,3*f+1)) 
		/ gsl_matrix_get (A->SIG,A->r,3*f+1) );

	/* NOE */
	if (MHz != gsl_matrix_get(A->X,A->r,3*f+2)) {
		MHz = gsl_matrix_get(A->X,A->r,3*f+2);
		wh= 1.0E-3*2.0*M_PI*MHz; // unit E+9 (rad/s)
		wx = wh*gamma_x/gamma_h;
		c = 1.0E+3*wx*csa_x/sqrt(3.0);
		MLE_Jw (0.0,  p,data,y0);
		MLE_Jw (wx,   p,data,y1);
		MLE_Jw (wh-wx,p,data,y2);
		MLE_Jw (wh,   p,data,y3);
		MLE_Jw (wh+wx,p,data,y4);
		} // skip Jw evaluation if MHz is same

	NOE = 1.0+fg*(6.0*y4-y2)/R1;

	gsl_vector_set (chi, 3*f+2, (NOE - gsl_matrix_get(A->Y,A->r,3*f+2)) 
		/ gsl_matrix_get (A->SIG,A->r,3*f+2) );

	}// field
}

void MLE_Jw (const double w,double *p,void *data, double &y)
/* w : 1E+9 rad/s */
/* y = spectral density */
{
const extern int D;
ALLDATA *A = (ALLDATA *)data;
double tk,wt,t,fac,sum_o,sum_i,sum_w;
int k;

double S2f = p[MLE_S2f] = 1.0;
double S2s = p[MLE_S2s];
double te  = p[MLE_te];

sum_o=sum_i=sum_w=0.0;

for (k=0;k<A->NK;k++) {
	if (D == _GLOBAL_DISTRIBUTION_) {
		tk = A->p[_scf_]*gsl_matrix_get(A->TAU,A->r,k);// (ns)
		wt = gsl_matrix_get(A->WT,A->r,k);
		}
	if (D == _GLOBAL_ISOTROPIC_ || D == _LOCAL_ISOTROPIC_) {
		tk = A->p[_tc_]; // (ns)
		wt = 1.0;
		}
	if (D == _GLOBAL_AXIAL_ || D == _GLOBAL_ANISOTROPIC_) {
		tk = A->t[k];// (ns)
		wt = gsl_matrix_get(A->WT,A->r,k);
		}
	if (tk > 0.0) {
		sum_w += wt;
		sum_o += wt * tk/(1.0+SQR(w*tk));// (ns)
		t  = te*tk/(te+tk);// (ns)
		sum_i += wt * t/(1.0+SQR(w*t));// (ns)
		} // tk > 0
	} // k

fac = 1.0E-9*(2.0/5.0)/sum_w;
y = fac*S2f*(S2s*sum_o + (1-S2s)*sum_i);
}
