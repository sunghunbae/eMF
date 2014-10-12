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

// get randomization number

void fitmodel (const int model, const int MC, ALLDATA &A)
{
  extern gsl_rng * rng;

  const extern double critx2;
  int i,j,k,mdata,mfit;
  bool valid, store_is[NP];
  double chisq, chisq_gr,store_p[NP];
  gsl_vector_view v;

  /* store is status */
  for (i=0;i<NP ;i++) store_is[i] = A.is[i];

  /* number of fitted parameters */
  mfit = set_is_local (A,model);

  /* number of data points */
  for(mdata=0,i=0;i<A.NF*3;i++)
    if (gsl_matrix_get(A.X,A.r,i) > 0) mdata++;

  /* degree of freedom */
  if (mdata-mfit < 0) return;
  gsl_matrix_set (A.dof,A.r,model,mdata-mfit);

  /* initialize parameter with grid search output */
  /* do grid search for new starting point */
  // grid_search (A, LX2_R1R2NOE, chisq_gr, false);

  /* randomization and simplex */

  /* error estimation from covariance matrix */
  gsl_vector *e = gsl_vector_alloc (NP);

  A.func = LX2_R1R2NOE;
  double chisq_min = 1e+9;
  size_t N = 100*mfit;
  while (N > 0) 
  {
    for(i=0;i<NP;i++) {
      if (A.is[i])
	A.p[i] = A.lb[i]+(A.ub[i]-A.lb[i])*gsl_rng_uniform(rng);
      }
    simplex (A,chisq_gr); 
    if (chisq_gr < chisq_min) {
      chisq_min = chisq_gr;
      for (i=0;i<NP;i++) store_p[i] = A.p[i];
      }
    N--;
  }
  chisq_gr = chisq_min;
  for (i=0;i<NP;i++) A.p[i] = store_p[i];

  /* store grid search results */
  for (i=0;i<NP;i++) store_p[i] = A.p[i];

  /* Levenberg-Marquardt Nonlinear Least Square Fitting */
  levmar (A, chisq, e);

  /* if LM gives worse chi-square than grid search */
  if (chisq > chisq_gr) {
    chisq = chisq_gr;
    for (i=0;i<NP;i++) A.p[i] = store_p[i];
    }

  gsl_matrix_set (A.x2, A.r,model,chisq);
  gsl_matrix_set (A.bic,A.r,model,chisq + mfit*log((double)mdata));
  gsl_matrix_set (A.aic,A.r,model,chisq + mfit*2);

  for (i=0;i<NP;i++)
    gsl_matrix_set (A.vP[i],A.r,model,A.p[i]);

  /* Monte Carlo Simulations */ 
  if (MC > 0) {
    /* store experimental data */
    double exp[A.NF*3], syn[A.NF*3];
    for (i=0;i<A.NF*3;i++)
      exp[i] = gsl_matrix_get(A.Y,A.r,i);
    gsl_vector *fitted = gsl_vector_alloc (A.NF*3);
    gsl_vector *simfit = gsl_vector_alloc (A.NF*3);
    R1R2NOE (&A,fitted,NULL);

    double sim;

    /* memory allocation for simulated data */
    gsl_matrix *mcpar = gsl_matrix_calloc(MC,NP);
    double *  mcx2 = new double [MC];
    extern int MC_trim;
    double simx2,s,ave,var,ep;
    int n, mc_trial = 0, mc_success = 0;

    do {
      /* add randomized noise to experimental data */
      /*
      for(j=0;j<A.NF;j++) {
	if (gsl_matrix_get(A.X,A.r,j*3+0) > 0) {
	  sim = gsl_ran_gaussian(rng, gsl_matrix_get(A.SIG,A.r,j*3+0));
	  gsl_matrix_set(A.Y,A.r,j*3+0,exp[j*3+0]+sim); 
	  }
	if (gsl_matrix_get(A.X,A.r,j*3+1) > 0) {
	  sim = gsl_ran_gaussian(rng, gsl_matrix_get(A.SIG,A.r,j*3+1));
	  gsl_matrix_set(A.Y,A.r,j*3+1,exp[j*3+1]+sim);
	  }
	if (gsl_matrix_get(A.X,A.r,j*3+2) > 0) {
	  sim = gsl_ran_gaussian(rng, gsl_matrix_get(A.SIG,A.r,j*3+2));
	  gsl_matrix_set(A.Y,A.r,j*3+2,exp[j*3+2]+sim);
	  }
	}
      */
      /* add randomized noise to fitted data */
      for(j=0;j<A.NF;j++) {
	if (gsl_matrix_get(A.X,A.r,j*3+0) > 0) {
	  sim = gsl_ran_gaussian(rng, gsl_matrix_get(A.SIG,A.r,j*3+0));
	  gsl_matrix_set(A.Y,A.r,j*3+0,gsl_vector_get(fitted,j*3+0)+sim); 
	  }
	if (gsl_matrix_get(A.X,A.r,j*3+1) > 0) {
	  sim = gsl_ran_gaussian(rng, gsl_matrix_get(A.SIG,A.r,j*3+1));
	  gsl_matrix_set(A.Y,A.r,j*3+1,gsl_vector_get(fitted,j*3+1)+sim);
	  }
	if (gsl_matrix_get(A.X,A.r,j*3+2) > 0) {
	  sim = gsl_ran_gaussian(rng, gsl_matrix_get(A.SIG,A.r,j*3+2));
	  gsl_matrix_set(A.Y,A.r,j*3+2,gsl_vector_get(fitted,j*3+2)+sim);
	  }
	}

      /* guess initial parameters from grid search */
      // grid_search (A, LX2_R1R2NOE, chisq_gr, false);

  A.func = LX2_R1R2NOE;
  double chisq_min = 1e+9;
  N = mfit*100;
  while (N > 0) {
    for(i=0;i<NP;i++) {
      if (A.is[i])
	A.p[i] = A.lb[i]+(A.ub[i]-A.lb[i])*gsl_rng_uniform(rng);
      }
    simplex (A,chisq_gr); 
    if (chisq_gr < chisq_min) {
      chisq_min = chisq_gr;
      for (i=0;i<NP;i++) store_p[i] = A.p[i];
      }
    N--;
    }
  chisq_gr = chisq_min;
  for (i=0;i<NP;i++) A.p[i] = store_p[i];

      /* store grid search results */
  for (i=0;i<NP;i++) store_p[i] = A.p[i];

      /* Levenberg-Marquardt Nonlinear Least Square Fitting */
      levmar (A,chisq,e);

      /* if LM gives worse chi-square than grid search */
      if (chisq > chisq_gr) {
	chisq = chisq_gr;
	for (i=0;i<NP;i++) A.p[i] = store_p[i];
	}

      /* check if fitted parameter is valid */
      valid = true;
      if (A.p[_S2s_] <  0 || A.p[_S2s_] >  1) valid = false;
      if (A.p[_S2f_] <  0 || A.p[_S2f_] >  1) valid = false;
      if (A.p[_te_]  <  0) valid = false;
      if (A.p[_Rex_] <  0) valid = false;
      if (valid) {
	for (i=0;i<NP;i++) 
	  gsl_matrix_set(mcpar,mc_success,i,A.p[i]);
	R1R2NOE (&A,simfit,NULL);
	for(simx2=0.0,j=0;j<A.NF;j++) {
	  if (gsl_matrix_get(A.X,A.r,j*3+0) > 0)
	    simx2 += SQR((gsl_vector_get(simfit,j*3+0)
		    -gsl_vector_get(fitted,j*3+0))/
		    gsl_matrix_get(A.SIG,A.r,j*3+0));
	  if (gsl_matrix_get(A.X,A.r,j*3+1) > 0)
	    simx2 += SQR((gsl_vector_get(simfit,j*3+1)
		    -gsl_vector_get(fitted,j*3+1))/
		    gsl_matrix_get(A.SIG,A.r,j*3+1));
	  if (gsl_matrix_get(A.X,A.r,j*3+2) > 0)
	    simx2 += SQR((gsl_vector_get(simfit,j*3+2)
		    -gsl_vector_get(fitted,j*3+2))/
		    gsl_matrix_get(A.SIG,A.r,j*3+2));
	  }
	mcx2[mc_success] = simx2;
	mc_success++;
	}
      mc_trial++;
      if (mc_trial > MC*100) exit(1);
      } while (mc_success < MC);
	
    /* restore experimental data */
    for (i=0;i<A.NF*3;i++)
      gsl_matrix_set(A.Y,A.r,i,exp[i]);
		
    /* simulated chi-square statistics */
    gsl_sort (mcx2,1,mc_success);
    for (i=0;i<20;i++) 
      A.percentile[i] = 
	gsl_stats_quantile_from_sorted_data (mcx2,1,mc_success,
	(double)(i+1)/20.0);

    // estimate error as standard error of simulated values
    for (i=0;i<NP;i++) {
      v = gsl_matrix_column (mcpar,i);
      gsl_sort_vector (&v.vector);
      // average 
      s=0.0;
      for (n=0,j=MC_trim;j<(mc_success-MC_trim);j++) {
	s += gsl_vector_get(&v.vector,j);
	n++;
	}
      ave=s/n;
      // standard deviation
      ep=var=0.0;
      for (j=MC_trim;j<(mc_success-MC_trim);j++) {
	s = gsl_vector_get(&v.vector,j)-ave;
	ep += s;
	var += s*s;
	}
      var=(var-ep*ep/n)/n;
      gsl_matrix_set (A.aP[i],A.r,model,ave);	    // mean
      gsl_matrix_set (A.dP[i],A.r,model,sqrt(var)); // stdev
      // confidence interval: gsl_cdf_tdist_Pinv (1.0-critx2,mc_success-1)
      // standard error = standard deviation / sqrt(n)
      //gsl_matrix_set (A.dP[i],A.r,model,
      //gsl_cdf_tdist_Pinv (1.0-critx2,mc_success-1)*sqrt(var/n));
      }
    gsl_matrix_free (mcpar);
    gsl_vector_free (fitted);
    gsl_vector_free (simfit);
    delete [] (mcx2);

    } // error estimated by standard deviation of MC fitting

  else {
    for (i=0;i<NP;i++)
      gsl_matrix_set (A.dP[i],A.r,model,gsl_vector_get(e,i));
    } // error estimated by covariance matrix if MC == 0

  gsl_vector_free (e);

  /* restore is status */
  for (i=0;i<NP;i++) A.is[i] = store_is[i];
}


void select_model(ALLDATA &A)
{
  extern double critx2;
  extern int NM,criterion;
  const int r = A.r;
  int best_model=0;// reset
  int second_best_model = 0; // reset
  bool valid;

  for(int m=1;m<=NM;m++) {
    /* consider only fitted models */
    if (gsl_matrix_get(A.x2,r,m) >= 0) {

      /* eliminate unrealistic models */
      valid = true;
      if (gsl_matrix_get(A.vP[_S2s_],r,m) <  0) valid = false;
      if (gsl_matrix_get(A.vP[_S2s_],r,m) >  1) valid = false;
      if (gsl_matrix_get(A.vP[_S2f_],r,m) <  0) valid = false;
      if (gsl_matrix_get(A.vP[_S2f_],r,m) >  1) valid = false;
      if (gsl_matrix_get(A.vP[_te_] ,r,m) <  0) valid = false;
      if (gsl_matrix_get(A.vP[_Rex_],r,m) <  0) valid = false;

      if (criterion == _BIC_) 
      {
        if (valid) 
        {
	  if (best_model == 0) best_model=m;
	  else if (gsl_matrix_get(A.bic,r,m) < gsl_matrix_get(A.bic,r,best_model)) best_model=m;
	} 
        else
        {
	  if (second_best_model == 0) second_best_model=m;
	  else if (gsl_matrix_get(A.bic,r,m) < gsl_matrix_get(A.bic,r,second_best_model)) second_best_model=m;
        }
      } // BIC

      if (criterion == _AIC_) 
      {
        if (valid)
        {
	  if (best_model == 0) best_model=m;
	  else if (gsl_matrix_get(A.aic,r,m) < gsl_matrix_get(A.aic,r,best_model)) best_model=m;
	} 
        else
        {
	  if (second_best_model == 0) second_best_model=m;
	  else if (gsl_matrix_get(A.aic,r,m) < gsl_matrix_get(A.aic,r,second_best_model)) second_best_model=m;
        }
      } // AIC

      }// for fitted models only
    } // model

  /* if all models are unrealistic, select the best BIC or AIC model */
  if (best_model != 0)
    A.best[r] = best_model;
  else
    A.best[r] = second_best_model;
}

int set_is_local (ALLDATA &A, int model)
{
/*
  Description of MODELS

  MODEL  0: Fix S2f=1, te=0, Rex=0, S2s=1
  MODEL  1: Fix S2f=1, te=0, Rex=0              Fit S2s
  MODEL  2: Fix S2f=1, Rex=0                    Fit S2s, te
  MODEL  3: Fix S2f=1, te=0                     Fit S2s, Rex
  MODEL  4: Fix S2f=1                           Fit S2s, te, Rex
  MODEL  5: Fix Rex=0                           Fit S2f, S2s, te(=ts)

  CAUTION: 
  active local diffusion parameters are fitted simultaneously unless fixed
*/

  int i,np;

  /* Include ONLY active Local parameters unless fixed */
  for (i=0;i<NP ;i++) {
    if ((A.attr[i] & P_LOCAL) && (A.attr[i] & P_ACTIVE) && 
      (A.attr[i] & P_FIXED) == 0)
      A.is[i] = true;
    else 	
      A.is[i] = false;
    }

  switch (model) {
    case 0:
      A.is[_S2f_]=false;A.p[_S2f_]=1.0;
      A.is[_S2s_]=false;A.p[_S2s_]=1.0;
      A.is[_te_]=false;A.p[_te_]=0.0;
      A.is[_Rex_]=false;A.p[_Rex_]=0.0;
      break;
    case 1:
      A.is[_S2f_]=false;A.p[_S2f_]=1.0;
      A.is[_S2s_]=true;
      A.is[_te_]=false;A.p[_te_]=0.0;
      A.is[_Rex_]=false;A.p[_Rex_]=0.0;
      break;
    case 2:
      A.is[_S2f_]=false;A.p[_S2f_]=1.0;
      A.is[_S2s_]=true;
      A.is[_te_]=true;
      A.is[_Rex_]=false;A.p[_Rex_]=0.0;
      break;
    case 3:
      A.is[_S2f_]=false;A.p[_S2f_]=1.0;
      A.is[_S2s_]=true;
      A.is[_te_]=false;A.p[_te_]=0.0;
      A.is[_Rex_]=true;
      break;
    case 4:
      A.is[_S2f_]=false;A.p[_S2f_]=1.0;
      A.is[_S2s_]=true;
      A.is[_te_]=true;
      A.is[_Rex_]=true;
      break;
    case 5:
      A.is[_S2f_]=true;
      A.is[_S2s_]=true;
      A.is[_te_]=true;
      A.is[_Rex_]=false;A.p[_Rex_]=0.0;
      break;
    }// switch

  for (np=0,i=0;i<NP;i++)
    if (A.is[i]) np++;

  return np;
}
