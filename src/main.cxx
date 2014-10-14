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
#include <getopt.h>

using namespace std;

/* default constants */
double 	gamma_h	= 26.7522;
double 	gamma_x	= -2.7126;
double 	r_xh	= 1.02;
double 	csa_x	= -160.0;

/* default spectral density function */
void (* Jw) (const double w,void *data,double &y,double *dyda) = &Jw_LSe; 

/* default diffusion tensor */
int	D = _NONE_;

set <int> cluster_set,residue_set;
set <double> field_set;

/* adjust experimental error */
map <double,double> min_err1,min_err2,min_err3;
map <double,double> scl_err1,scl_err2,scl_err3;
map <double,double>::iterator e,e0;

/* PDB */
string 	pdbfile,atomx,atomh;

/* Maximum Likelihood Estimation of diffusion tensor */
int 	MLE_MC = 100000;
int 	MLE_cluster=-1;

/* diffusion tensor estimation by R2/R1 */
bool 	estimate = false;
int     estimateFunc = GX2_R2_OVER_R1;
double 	EstNOECut = 0.0;
double 	EstSTDCut = 0.0;
int 	EstCluster=-1;

/* axial grid search */
bool	axial_grid_search = false;

/* optimization of diffusion tensor */
bool 	optimize = false;
int 	OptMaxIter = 500;
int	OptMethod  = 0;
int 	OptCluster = -1;
double 	OptS2 = 0.0;

/* model selection & goodness of fit */
int 	NM = 5;
int 	criterion = _BIC_;
double 	critx2 = 0.05;

/* error estimation by Monte Carlso simulation */
int 	MC	= 0;
int 	MC_trim = 0;

/* for recursive grid search routine */
size_t 	grid_index_end;
double 	grid_min_par[NP];
double 	grid_min_x2,g_PROB;

/* internal use: seed for random number generation */
unsigned long int seed = time(0);
gsl_rng * rng;

/* long output */
bool opt_verb = false;

int main(int argc, char *argv[])
{
  /* signature */
  printf("  eMF 1.1 by Sung-Hun Bae 2008-2014\n\n");

  const gsl_rng_type * T = gsl_rng_ranlux389;
  rng = gsl_rng_alloc (T);
  gsl_rng_set (rng, seed);

  /* checking elapsed time */
  time_t start,finish;
  time(&start);

  int i;

  char *DATA_FILE = NULL;
  char *CONFIG_FILE = NULL;
  char *XMGR_PREFIX = NULL;

  int NF;// number of magnetic fields
  int NR;// number of residues
  int NK;// number of tau points
  int r,f,m;

  /* commands arguments */

  bool do_rjmap = false;
  bool do_test = false;
  bool idtk;

  int c;
  opterr = 0;
  while ((c = getopt(argc, argv, "hvsmtLIANDc:x:"))!= -1)
    switch (c) {
      case 'h': print_help(); break;
      case 'v': opt_verb = true; break;
      case 's': do_rjmap = true; break;
      case 't': do_test = true; break;
      case 'L': D = _LOCAL_ISOTROPIC_; break;
      case 'I': D = _GLOBAL_ISOTROPIC_; break;
      case 'A': D = _GLOBAL_AXIAL_; break;
      case 'N': D = _GLOBAL_ANISOTROPIC_; break;
      case 'D': D = _GLOBAL_DISTRIBUTION_; break;
      case 'c': CONFIG_FILE = optarg; break;
      case 'x': XMGR_PREFIX = optarg; break;
      case '?': print_help(); break;
      }// switch

  DATA_FILE = argv[optind]; 
  if (DATA_FILE == NULL) print_help();

  /* need tau distribution file ? */
  if (D == _GLOBAL_DISTRIBUTION_) idtk = true;
  else idtk = false;

  /* check data file & Allocate memory */
  chk_data (DATA_FILE,NF,NR,NK,idtk); 
  if (!idtk) {
    if (D & _BIMODAL_ ) NK = 2;
    else if (D & _AXIAL_ ) NK = 3;
    else if (D & _ANISOTROPIC_) NK = 5;
    else NK = 1;
    }

  /* allocate memory */
  ALLDATA A;
  initialize (A,NM,NF,NR,NK,MC);

  /* read configuration file */
  show_info("configuration");
  if (CONFIG_FILE == NULL) {
    printf("emf.conf\n");
    read_config ("./emf.conf",A);
    }
  else {
    printf("%s\n",CONFIG_FILE);
    read_config (CONFIG_FILE,A);
    }

  show_info("  gamma H");printf("%g\n",gamma_h);
  show_info("  gamma X");printf("%g\n",gamma_x);
  show_info("  CSA X (ppm)");printf("%g\n",csa_x);
  show_info("  bond length (A)");printf("%g\n",r_xh);

  show_info("  estimation of diffusion");
  if (estimate) {
    printf("Yes\n");
    show_info("   hnNOE cutoff");printf("%g\n",EstNOECut);
    show_info("   R2 stdev cutoff");printf("%g\n",EstSTDCut);
    show_info("   cluster");printf("%d\n",EstCluster);
    }
  else printf("No\n");

  show_info("  optimization of diffusion");
  if (optimize) {
    printf("Yes\n");
    show_info("    max. iteration");printf("%d\n",OptMaxIter);
    show_info("    cluster");printf("%d\n",OptCluster);
    show_info("    S2 cutoff");printf("%g\n",OptS2);
    }
  else printf("No\n");

  show_info("  model selection criterion");
  if (criterion == _BIC_) printf("BIC\n");
  if (criterion == _AIC_) printf("AIC\n");

  show_info("  Monte Carlo simulation");
  if (MC > 0) {
    printf("Yes\n");
    show_info("    number of simulations");printf("%d\n",MC);
    show_info("    number of trims");printf("%d\n",MC_trim);
    }
  else printf("No\n");

  /* read DATA file */
  read_data(DATA_FILE,A,idtk);

  /* check if errors are adjusted */
  bool error_adjusted = false;
  for (e0 = e = min_err1.begin();e!=min_err1.end();e++)
    if (e->second != 0) error_adjusted = true;
  for (e0 = e = scl_err1.begin();e!=scl_err1.end();e++)
    if (e->second != 0) error_adjusted = true;
  for (e0 = e = min_err2.begin();e!=min_err2.end();e++)
    if (e->second != 0) error_adjusted = true;
  for (e0 = e = scl_err2.begin();e!=scl_err2.end();e++)
    if (e->second != 0) error_adjusted = true;
  for (e0 = e = min_err3.begin();e!=min_err3.end();e++)
    if (e->second != 0) error_adjusted = true;
  for (e0 = e = scl_err3.begin();e!=scl_err3.end();e++)
    if (e->second != 0) error_adjusted = true;
  show_info("  error adjusted");
  if (error_adjusted) {
    printf("Yes\n");
    show_info("    R1 error minimum (%)");
    for (e0 = e = min_err1.begin();e!=min_err1.end();e++) {
      if (e != e0) printf(", %g",100*e->second);
      else printf("%g",100*e->second);
      }
    printf("\n");

    show_info("    R1 error scale");
    for (e0 = e = scl_err1.begin();e!=scl_err1.end();e++) {
      if (e != e0) printf(", %g",e->second);
      else printf("%g",e->second);
      }
    printf("\n");

    show_info("    R2 error minimum (%)");
    for (e0 = e = min_err2.begin();e!=min_err2.end();e++) {
      if (e != e0) printf(", %g",100*e->second);
      else printf("%g",100*e->second);
      }
    printf("\n");

    show_info("    R2 error scale");
    for (e0 = e = scl_err2.begin();e!=scl_err2.end();e++) {
      if (e != e0) printf(", %g",e->second);
      else printf("%g",e->second);
      }
    printf("\n");

    show_info("    hnNOE error minimum");
    for (e0 = e = min_err3.begin();e!=min_err3.end();e++) {
      if (e != e0) printf(", %g",e->second);
      else printf("%g",e->second);
      }
    printf("\n");

    show_info("    hnNOE error scale");
    for (e0 = e = scl_err3.begin();e!=scl_err3.end();e++) {
      if (e!= e0) printf(", %g",e->second);
      else printf("%g",e->second);
      }
    printf("\n");
    }

    else printf("No\n");

  /* read PDB file */
  if (D == _GLOBAL_AXIAL_ || D == _GLOBAL_ANISOTROPIC_) {
    if (pdbfile.empty()) {
      printf(_ERROR_ "pdb file is required but not defined\n");
      exit(1);
      }	
    else {
      int nvec = read_pdb(pdbfile.c_str(),A);
      show_info("PDB coordinate file");printf("%s\n",pdbfile.c_str());
      show_info("  number of vectors");printf("%d\n",nvec);
      show_info("  atom H");printf("%s\n",atomh.c_str());
      show_info("  atom X");printf("%s\n",atomx.c_str());
      if (nvec == 0) {
        printf(_ERROR_ "bond vectors are not defined\n");
        exit(1);
        }
      }
    }

  /* diffusion tensor */
  show_info("diffusion tensor");
  if (D == _NONE_)		  printf("not defined\n");
  if (D == _GLOBAL_ISOTROPIC_)	  printf("global isotropic\n");
  if (D == _GLOBAL_AXIAL_)	  printf("global axial\n");
  if (D == _GLOBAL_BIMODAL_)	  printf("global bimodal\n");
  if (D == _GLOBAL_ANISOTROPIC_)  printf("global anisotropic\n");
  if (D == _GLOBAL_DISTRIBUTION_) printf("global distribution\n");
  if (D == _LOCAL_ISOTROPIC_)	  printf("local isotropic\n");
  if (D == _LOCAL_BIMODAL_)	  printf("local bimodal\n");
  if (D == _LOCAL_COLE_COLE_)	  printf("local Cole-Cole\n");
  if (D == _LOCAL_LORENTZIAN_)	  printf("local Lorentzian\n");

  /* set up parameter attribute */
  setup_attribute (A);

  /* spectral density function type */
  show_info("spectral density function");
  if (Jw == &Jw_LSe) printf("Lipari-Szabo extended\n"); 

  /* break */
  printf("\n");
  show_parameter(A);
  printf("\n");

  /* show data */
  if(opt_verb) {
    show_data(A);
    printf(_INFO_ "\n");
    }

  if (do_test) {
    /* reserved for test */
    Jw = &Jw_LSe;
/*
    D = _LOCAL_ISOTROPIC_;
    printf(">>>>>> generating test data\n");
    //A.p[_c_] = 1.0;
    A.p[_tc_] = 10.;
    A.p[_S2s_] = 0.97;
    A.p[_S2f_] = 1.0;
    A.p[_te_] = 2.048;
    A.p[_Rex_] = 0.149;
    gsl_vector * relax = gsl_vector_calloc (3*A.NF);
    R1R2NOE((void *)&A,relax,NULL);
    for (int i=0;i<A.NF;i++) {
      printf("%f %f %f %f\n",gsl_matrix_get(A.X,0,3*i+0),
	gsl_vector_get(relax,3*i+0),
	gsl_vector_get(relax,3*i+1),
	gsl_vector_get(relax,3*i+2));
      }
    gsl_vector_free (relax);
*/
/*
    D = _GLOBAL_ISOTROPIC_;
    setup_attribute (A);
    printf(">>>>> testing chi_w_f\n");
    A.p[_S2s_] = 0.67;
    A.p[_S2f_] = 1.0;
    A.p[_te_] = 1.048;
    A.p[_Rex_] = 2.149;
    A.flag[0] = true;
    A.best[0] = 4;
    A.is[_tc_] = true;
    A.is[_S2s_] = true;
    A.is[_te_] = true;
    A.is[_Rex_] = true;
    gsl_vector *chi = gsl_vector_alloc(3*A.NF);
    gsl_vector *p = gsl_vector_alloc (4);
    gsl_vector_set (p,0,10);
    gsl_vector_set (p,1,0.97);
    gsl_vector_set (p,2,2.048);
    gsl_vector_set (p,3,0.149);
    chi_w_f(p,&A,chi);
    printf("chisq = %g\n",gsl_blas_dnrm2(chi));
    gsl_matrix *J = gsl_matrix_alloc (3*A.NF,4);
    chi_w_fdf(p,&A,chi,J);
    printf("chisq = %g\n",gsl_blas_dnrm2(chi));
    for (int i=0;i<3*A.NF;i++) {
      printf("jacob: ");
      for(int j=0;j<4;j++) 
	printf("%9.5e ",gsl_matrix_get(J,i,j));
      printf("\n");
      }
    gsl_matrix_free(J);
    gsl_vector_free(p);
    gsl_vector_free(chi);
*/
    D = _GLOBAL_AXIAL_;
    setup_attribute (A);
    A.p[_S2s_] = 0.8;
    A.p[_S2f_] = 1.0;
    A.p[_te_]  = 0.05;
    A.p[_Rex_] = 0.0;
    A.p[_Dr_]  = 1.9;
    A.p[_tc_]  = 6.4;
    A.p[_theta_] = 0.0;
    A.p[_phi_]   = 0.0;
    gsl_vector * relax = gsl_vector_calloc (3*A.NF);
    for (int r=0;r<A.NR;r++) {
      if (A.ivec[r]) {
	A.r = r;
	R1R2NOE((void *)&A,relax,NULL);
	for (int i=0;i<A.NF;i++) {
	  printf("%3d 1 %.3f %f 0.001 %f 0.001 %f 0.001\n",A.num[r],
	    gsl_matrix_get(A.X,0,3*i+0),
	    gsl_vector_get(relax,3*i+0),
	    gsl_vector_get(relax,3*i+1),
	    gsl_vector_get(relax,3*i+2));
	  }
	}
      }

    gsl_vector_free (relax);

    exit(1);
    }

  /* FIT GLOBAL MODELS options: -I, -A, -N, -D */

  // if diffusion is set to global
  if (D & _GLOBAL_) {
    //gsl_vector *global = gsl_vector_alloc (NP);
    double chisq, chisq_gr;
    int iter=0;

    /* 
      estimate rotational diffusion tensor:
	selected residues are assumed to be rigid
	initial diffusion parameters will be overwritten
	fit rotational diffusion tensor using R2/R1. 
    */

    if (estimate) {
      printf(_INFO_ "estimating diffusion tensor\n");
      select_all (A);
      select_estimator (A,EstCluster,EstNOECut,EstSTDCut);
      printf(_INFO_ "residues for estimation\n");
      print_selected(A);
      A.p[_te_] =0.01;
      A.p[_Rex_]=0.0;
      A.p[_S2f_]=1.0;
      A.p[_S2s_]=0.8;

      for (i=0;i<NP;i++) {
	if ((A.attr[i] & P_LOCAL) == 0 && (A.attr[i] & P_ACTIVE) && 
	  (A.attr[i] & P_FIXED) == 0)
	  A.is[i] = true;
	else
	  A.is[i] = false;
	}

      if (estimateFunc == GX2_R2_OVER_R1 ) 
      {
        grid_search (A, GX2_R2_OVER_R1, chisq_gr, true);

        global_minimize(A, GX2_R2_OVER_R1, chisq, M_SIMPLEX);

        gsl_vector * relax = gsl_vector_calloc (A.NF*3);
        for (A.r=0;A.r<A.NR;A.r++) 
        {
          R1R2NOE (&A,relax,NULL);
          for (f=0;f<A.NF;f++) 
          {
	    /* only if both R1 and R2 are available */
	    if (gsl_matrix_get(A.X, A.r, f*3+0) > 0.0 &&
	        gsl_matrix_get(A.X, A.r, f*3+1) > 0.0 && A.flag[A.r]) 
            {
              double MHz = gsl_matrix_get(A.X,A.r,f*3+0);
	      double R1 = gsl_matrix_get(A.Y,A.r,f*3+0);
	      double R2 = gsl_matrix_get(A.Y,A.r,f*3+1);
	      double dR1 = gsl_matrix_get(A.SIG,A.r,f*3+0);
	      double dR2 = gsl_matrix_get(A.SIG,A.r,f*3+1);
	      double sig2 = SQR(R2/R1)*(SQR(dR1/R1)+SQR(dR2/R2));
	      double fR1 = gsl_vector_get(relax,3*f+0); // fitted R1
	      double fR2 = gsl_vector_get(relax,3*f+1); // fitted R2
	      double dy=R2/R1-fR2/fR1;
	      double chisq = dy*dy/sig2;
              printf(_INFO_ "RESIDUE %4d %7.3f MHz exp= %7.3f cal= %7.3f dy= %7.3f chisq= %7.3f\n",
                A.num[A.r],MHz,R2/R1,fR2/fR1,dy,chisq);
	    }
	  }
        }
        gsl_vector_free (relax);
      }
    

      if (estimateFunc == GX2_R2_OVER_R1_SIMPLE ) 
      {
        grid_search (A, GX2_R2_OVER_R1_SIMPLE, chisq_gr, true);
        global_minimize(A, GX2_R2_OVER_R1_SIMPLE, chisq, M_SIMPLEX);

        const double c0 = gamma_x/(5.0*gamma_h);
        const double c1 = 7.0*SQR(0.921/0.87);
        const double c2 = 13.0/2.0*SQR(0.955/0.87);
        for (A.r=0;A.r<A.NR;A.r++) 
        {
          for (f=0;f<A.NF;f++) 
          {
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

	      double dy=rho-jwx/jw0;
	      double chisq = dy*dy/sig2;
              printf(_INFO_ "RESIDUE %4d %7.3f MHz exp= %7.3f cal= %7.3f dy= %7.3f chisq= %7.3f\n",
                A.num[A.r],MHz,rho,jwx/jw0,dy,chisq);
	    } // conditions
	  } // f
        } // A.r
      }

      } // estimate

      /* 
	optimize rotational diffusion tensor:
	with given diffusion tensor, fit models and
	find the best model. during model fitting,
	only internal dynamic variables are changed 
      */
	
      while (optimize && iter < OptMaxIter) {

	print_iteration(iter,A);

	if (iter == 0) 
	  select_optimizer (A,OptCluster,-1);
	else
	  select_optimizer (A,OptCluster,OptS2);

	for (chisq=0.0,r=0;r<NR;r++) 
        { 
	  if(A.flag[r]) 
          {
	    A.r = r;
	    for (m=1;m<=NM;m++) 
	      fitmodel(m,0,A);
	    select_model(A);
	    chisq += gsl_matrix_get(A.x2,r,A.best[r]);
	  }// flag
	}// r

	select_optimizer (A,OptCluster,OptS2);
	print_selected(A);

	if (global_minimize (A, GX2_R1R2NOE_BM, chisq, M_LEVMAR)) break;

        // store estimated parameters including global diffusion parameters
        //for(i=0;i<NP;i++)
        //  gsl_vector_set(global,i,A.p[i]);

	iter++;

	}// optimize

        if (optimize & (iter == OptMaxIter))
        {
	    printf(_ERROR_ "not converged at max. iteration\n");
	    exit(1);
        }

	/*
	  internal parameters are fitted 
	  with given optimized diffusion tensor
	  Monte Carlo simulation for error estimation
	  if MC is not chosen, covariance matrix is used instead
	*/
	
	select_all (A); 

	for (r=0;r<NR;r++) 
        { 
	  A.r = r;
	  /* if best model is not determined yet */
	  if (A.best[r] == 0) 
          {
	    for (m=1;m<=NM;m++) 
            {
/*
              for (i=0; i<NP; i++)
              {
                  gsl_matrix_set (A.vP[i],r,m, gsl_vector_get(global,i));
              }
*/
              fitmodel(m,0,A);
            }
	    select_model(A);
	  }

	  /* if MC simulations are to be run */ 
	  if (A.best[r] != 0 && MC > 0) 
          {
	    fitmodel(A.best[r],MC,A);
	  }

	  show_models(A, XMGR_PREFIX);
	  printf("\n");
	}// r

        //gsl_vector_free(global);

	} // -I,-A,-N,-D

  /* FIT LOCAL */
  if (D & _LOCAL_) {
    for (r=0;r<NR;r++) {
      A.r = r;
      for (m=1;m<=NM;m++) fitmodel(m,0,A);
      select_model(A);
      if (MC > 0 && A.best[r]!=0) fitmodel(A.best[r],MC,A);
      show_models(A, XMGR_PREFIX);
      printf("\n");
      }// r
    } // -L


  /* 
    REDUCED SPECTRAL DENSITY MAPPING 
    1. direct simple reduced spectral density mapping
    2. plot continuous spectral density function with 
      reduced spectral density points in xmgr format.
      it should be called after fitting.
  */

  if (do_rjmap) {
    printf(_INFO_ "REDUCED SPECTRAL DENSITY\n");
    printf(_INFO_ "\n");
    for (r=0;r<NR;r++) {
      A.r = r;
      rjmap (A, (r == 0 ? true: false), NULL);
      }
    }// do_rjmap

  /* EXIT */

  gsl_rng_free (rng);

  freeing (A); // freeing memory
  time(&finish);
  int hh=0,mm=0,ss=0,elapsed=finish-start;
  hh = (int)floor((double)(elapsed/3600));
  mm = (int)floor((double)((elapsed/60)-(hh*60)));
  ss = (elapsed - (hh*3600)-(mm*60));
  show_info("started");printf("%s",ctime(&start));
  show_info("finished");printf("%s",ctime(&finish));
  show_info("time elapsed");printf("%02d:%02d:%02d\n",hh,mm,ss);
}
