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

#ifndef _EMF_H_
#define _EMF_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <ctime>
#include <cmath>
#include <complex>
#include <list>
#include <map>
#include <set>
#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics.h>

using namespace std;

#define _INFO_		"# "
#define _WARNING_	"WARNING "
#define _ERROR_		">>>>> "
#define _ITEM_WIDTH_	30
#define _INDENT_	5
#define _MAXSTR_	256

#define mu0_h_8pi2		1.0545887 // unit E-41

/* rotational diffusion tensor type */

#define _NONE_			0
#define	_GLOBAL_		1	
#define	_LOCAL_			2
#define _ISOTROPIC_		4
#define _AXIAL_			8
#define _ANISOTROPIC_		16
#define _DISTRIBUTION_		32
#define _BIMODAL_		64
#define	_COLE_COLE_		128
#define	_LORENTZIAN_		256
#define _WINDOW_		512
#define _GLOBAL_ISOTROPIC_	5
#define _GLOBAL_AXIAL_		9
#define _GLOBAL_ANISOTROPIC_	17
#define _GLOBAL_DISTRIBUTION_	33
#define _GLOBAL_BIMODAL_	65
#define _LOCAL_ISOTROPIC_	6
#define _LOCAL_DISTRIBUTION_	34
#define _LOCAL_BIMODAL_		66
#define _LOCAL_COLE_COLE_	130
#define _LOCAL_LORENTZIAN_	258
#define _WINDOW_ISOTROPIC_	516

/* Parameters */
#define _S2f_	  0
#define _S2s_	  1
#define _te_	  2
#define _Rex_	  3
#define _tc_	  4
#define _Dr_	  5 
#define _scf_	  6
#define _Dxx_	  7
#define _Dyy_	  8
#define _Dzz_	  9
#define _phi_	  10
#define _theta_	  11
#define _psi_	  12
#define _tb_	  13
#define _c_	  14
#define NP	  15

/* parameter status (bitwise flag) */
#define P_RESET	  0
#define P_ACTIVE  1	/* unset= inactive */
#define P_FIXED	  2	/* unset= to be fitted */
#define P_LOCAL	  4	/* unset= global */
#define P_ROTDIF  8	/* unset= internal dynamics */

/* target function type */
#define LX2_R1R2NOE	1
#define GX2_R2_OVER_R1 	2
#define GX2_R2_OVER_R1_SIMPLE 	3
#define GX2_R1R2NOE_FM 	4
#define GX2_R1R2NOE_BM 	5
#define GX2_LIKELIHOOD	6

/* global minimization method */
#define M_AUTO	  0	
#define M_BRENT	  1
#define M_SIMPLEX 2
#define M_POWELL  3
#define M_CONJGR  4
#define M_LEVMAR  5
#define M_LEVMARW 6
#define M_ANNEAL  7

/* model selection criterion */
#define _BIC_	1
#define _AIC_	2
#define _STU_	3

typedef struct {
  /* info */
  int MA,NF,NR,NK;
  bool idtk,*flag; // NR
  int *num,*clst,*best; // NR

  /* function arguments (internal use) */
  int func,r,fpar;
  double *u,*v; // MA

  /* experimental data */
  gsl_matrix *X,*Y,*SIG;// NR x (NF*3)
  gsl_matrix *TAU,*WT;// NR x NK

  /* parameters */
  double *p,*lb,*ub,*lk,*uk,*step,*conv; // MA
  int *grds;
  char *attr;
  bool *is; // MA
  char pid[NP][80];
  char fid[NP][80];

  /* diffusion */
  bool *ivec;
  double *t;// max.5
  double *x,*y,*z,*alpha; // NR

  /* model info. NR x (NM+1) */
  gsl_matrix *x2,*dof,*aic,*bic;
  gsl_matrix *vP[NP],*dP[NP],*aP[NP];

  /* Monte Carlo simulation (internal use) */
  double *percentile;
  } ALLDATA;

/* func.c */
void  R1R2NOE (void *data,gsl_vector *,gsl_matrix *);
void  Jw_LSe  (const double,void *,double &,double *);
void  Jw_CC   (const double,void *,double &,double *);
void  Jw_Lz   (const double,void *,double &,double *);

/* chisq.c */
int chi_R1R2NOE_f (const gsl_vector *,void *,gsl_vector *);
int chi_R1R2NOE_df (const gsl_vector *,void *,gsl_matrix *);
int chi_R1R2NOE_fdf (const gsl_vector *,void *,gsl_vector *,gsl_matrix *);
int chi_gbm_f (const gsl_vector *,void *,gsl_vector *); 
int chi_gbm_df (const gsl_vector *,void *,gsl_matrix *);
int chi_gbm_fdf (const gsl_vector *,void *,gsl_vector *,gsl_matrix *);
int chi_g21_f (const gsl_vector *,void *,gsl_vector *); 
int chi_g21_df (const gsl_vector *,void *,gsl_matrix *);
int chi_g21_fdf (const gsl_vector *,void *,gsl_vector *,gsl_matrix *);
int chi_gfm_f (const gsl_vector *,void *,gsl_vector *); 
int chi_w_f (const gsl_vector *,void *,gsl_vector *); 
int chi_w_df (const gsl_vector *,void *,gsl_matrix *);
int chi_w_fdf (const gsl_vector *,void *,gsl_vector *,gsl_matrix *);
double	chi2 (ALLDATA &, int);
double	boundary_penalty (void *);

/* anisotropic.c */
void  anisotropic (ALLDATA &);
int   read_pdb (const char *,ALLDATA &);

/* rjmap.c */
void  rjmap (ALLDATA &, bool, FILE *);

/* gridsearch.c */
void  ax_grid_search (ALLDATA &);
void  grid_search (ALLDATA &, int, double &,bool);
void  recursive_s (ALLDATA &, int, int);

/* fitmodel.c */
void  fitmodel (const int,const int,ALLDATA &);
void  select_model (ALLDATA &);
int   set_is_local (ALLDATA &, int);

/* levmar.c */
void  levmar (ALLDATA &, double &, gsl_vector *);
void  levmar (ALLDATA &, double &);
void  levmarw(ALLDATA &, double &);

/* gop.c */
bool  global_minimize (ALLDATA &,int,double &,int);

/* brent.c */
void 	brent (ALLDATA &, double &);
double 	brent_f (double, void *);

/* simplex.c */
void   	simplex (ALLDATA &, double &);
double 	simplex_f (const gsl_vector *, void *);

/* conjgr.c */
void  conjgr (ALLDATA &, double &);
double	conjgr_f (const gsl_vector *,void *);
double	conjgr_nf (double,void *);
void  conjgr_df (const gsl_vector *,void *,gsl_vector *);
void  conjgr_fdf (const gsl_vector *,void *,double *,gsl_vector *);
void  conjgr_ndf (const gsl_vector *,void *,gsl_vector *);
void  conjgr_nfdf (const gsl_vector *,void *,double *,gsl_vector *);

/* powell.c */
void    powell (ALLDATA &, double &);
double  powell_f (double *, void *);
void    powell_linmin (double *, void *, double &);
double  powell_vector (double, void *);

/* anneal.c */
void	anneal (ALLDATA &, double &);

/* mle.c */
void  marginal_density (size_t, ALLDATA &,double &,double &);
double	likelihood (double *,size_t,void *);
void  MLE_chi_R1R2NOE (double *p,void *data,gsl_vector *chi);
void  MLE_Jw (const double w, double *p,void *data,double &y);

/* file.c */
void  chk_data (const char *,int &,int &,int &,bool);
void  read_data (const char *,ALLDATA &,bool);
void  read_config (const char *,ALLDATA &);
void  ltrim (string &);
void  rtrim (string &);
void  trim  (string &);
void  parse (const char *,vector <string> &);

/* data.c */
void  initialize (ALLDATA &,int,int,int,int,int);
void  freeing (ALLDATA &);
void  select_all (ALLDATA &);
void  select_window (ALLDATA &,int,int);
void  select_optimizer (ALLDATA &,int,double);
void  select_estimator (ALLDATA &,int,double,double);
void  select_cluster (int, ALLDATA &);
void  setup_attribute (ALLDATA &);
bool  diffusion_converged (ALLDATA &, double *);

/* print.c */
void  print_help ();
void  terminate (const string);
void  show_data (ALLDATA &);
void  show_text (const string);
void  print_selected (ALLDATA &);
void  print_iteration (int,ALLDATA &);
void  show_parameter (ALLDATA &);
void  show_models (ALLDATA &, char *);
void  show_info (const string);

inline double SQR (double a) {return a*a;}
inline double deg2rad (double deg) { return deg * (M_PI / 180.0); }
inline double rad2deg (double rad) { return rad * (180.0 / M_PI); }

#endif /* _EMF_H_ */
