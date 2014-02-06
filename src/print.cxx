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
#include <gsl/gsl_cdf.h>

#define XMGR_PLOT_W_PTS	500
#define XMGR_PLOT_W_MAX	0.8   /* w: 1e+9 Hz */

using namespace std;

void print_help()
{
  /* usage */
  printf("  Usage: emf [commands] Data_File\n");
  printf("  -h                display this help and exit\n");
  printf("  -v                print long output\n");
  printf("  -L                diffusion model = Local\n");
  printf("  -I                diffusion model = Isotropic\n");
  printf("  -A                diffusion model = Axially symmetric\n");
  printf("  -N                diffusion model = aNisotropic\n");
  printf("  -D                diffusion model = Distributed\n");
  printf("  -x [xmgr_prefix]  xmgr output for calc. spectral densities\n");
  printf("  -c [config_file]  read configuration (default: emf.conf)\n");
  printf("  -s                reduced spectral density mapping\n");
  printf("\n");

  /* data file format */
  printf("  Data_File format (text file, # = comment)\n");
  printf("   column 1-3 : [resid #] [cluster #] [magnetic field (MHz)]\n");
  printf("   column 4-5 : [R1 (1/s)] [R1 error (1/s)]\n");
  printf("   column 6-7 : [R2 (1/s)] [R2 error (1/s)]\n");
  printf("   column 8-9 : [hetNOE] [hetNOE error]\n");
  printf("   column 10  : tc distribution file name (only for -D)\n");
  printf("\n");
  exit(1);
}

void terminate(const string error)
{
cout << endl;
cout << _WARNING_ << "ERROR " << error << endl;
cout << _WARNING_ << "terminated" << endl;
exit(1);
}

void show_data (ALLDATA &A)
{
  int r,f;
  for(r=0;r<A.NR;r++) {
    printf(_INFO_ "RESIDUE %4d CLUSTER %1d\n",A.num[r],A.clst[r]);
    for(f=0;f<A.NF;f++) {
      if (gsl_matrix_get(A.X,r,f*3+0) > 0)
	printf(_INFO_ "R1  %7.3f %7.3f %7.3f\n",
	  gsl_matrix_get(A.X,r,f*3+0),
	  gsl_matrix_get(A.Y,r,f*3+0),
	  gsl_matrix_get(A.SIG,r,f*3+0));
      if (gsl_matrix_get(A.X,r,f*3+1) > 0)
	printf(_INFO_ "R2  %7.3f %7.3f %7.3f\n",
	  gsl_matrix_get(A.X,r,f*3+1),
	  gsl_matrix_get(A.Y,r,f*3+1),
	  gsl_matrix_get(A.SIG,r,f*3+1));
      if (gsl_matrix_get(A.X,r,f*3+2) > 0)
	printf(_INFO_ "NOE %7.3f %7.3f %7.3f\n",
	  gsl_matrix_get(A.X,r,f*3+2),
	  gsl_matrix_get(A.Y,r,f*3+2),
	  gsl_matrix_get(A.SIG,r,f*3+2));
      } // f
    } // r
}

void print_selected(ALLDATA &A)
{
  int c,r;
  for (c=0,r=0;r<A.NR;r++) {
    if(A.flag[r]) {
      if (c%10==0) printf(_INFO_);
      if (A.best[r] == 0) printf("%4d  ",A.num[r]);
      if (A.best[r] != 0) printf("%4d:%1d ",A.num[r],A.best[r]);
      if (c%10==9) printf("\n");
	c++;
	}// 10 residues per line
      }// flag
  if (c%10!=0) printf("\n");
}

void print_iteration(int iter,ALLDATA &A)
{
  int i;
  double unit;
  printf(_INFO_ "iteration %2d ( ",iter);
  for(i=0;i<NP;i++) {
    if ((A.attr[i] & P_ACTIVE) && (A.attr[i] & P_ROTDIF) && 
      ((A.attr[i] & P_LOCAL) == 0)) {
      if (i == _phi_ || i == _theta_) unit = 180.0/M_PI;
      else unit = 1.0;
      printf("%s %g ",A.pid[i],unit*A.p[i]);
      }
    }
  printf(")\n");
}

void show_parameter (ALLDATA &A)
{
  double unit;
  printf(_INFO_ "   %-5s %8s [ %7s ... %7s | %3s, %7s ] %s\n",
    "Param","Value","Lower","Upper","GRD","Step","Convergence");
  for (int i=0;i<NP;i++) {
    if (A.attr[i] & P_ACTIVE) {
      printf(_INFO_);
      if (A.attr[i] & P_LOCAL) 
	printf(" "); 
      else 
	printf("G");
      if (A.attr[i] & P_FIXED) 
	printf("F"); 
      else 
	printf(" ");
      if (i == _phi_ || i == _theta_ || i == _psi_ ) unit = 180.0/M_PI;
      else unit = 1.0;
      printf(" %-5s %8.03f [ %7.03f ... %7.03f | %3d, %7.1e ] %7.1e\n",
	A.pid[i],unit*A.p[i],unit*A.lb[i],unit*A.ub[i],
	A.grds[i],unit*A.step[i],unit*A.conv[i]);
      } // if 
    } // for
}

void show_models (ALLDATA &A, char *XMGR_PREFIX)
{
  extern void (* Jw) (const double w,void *data,double &y,double *dyda);
  const extern double critx2;
  const extern int D,MC,NM;
  const int r=A.r, b=A.best[r];
  double MHz,unit;
  int i,f,m;

  /* diffusion parameters */
  for (i=0;i<NP;i++) 
    if (A.attr[i] & P_ROTDIF) {
      if (i == _phi_ || i == _theta_) unit = 180.0/M_PI;
      else unit = 1.0;
      printf("RESIDUE %-5d %5s %8.03f\n",A.num[r],A.pid[i],
	unit*gsl_matrix_get(A.vP[i],r,b));
      }
  /* additional information for axial */
  if (D & _AXIAL_) 
    printf("RESIDUE %-5d alpha %8.03f\n",A.num[r],A.alpha[r]*180.0/M_PI);

      
  /* model parameters */
  for (m=1;m<=NM;m++) {

    if (gsl_matrix_get(A.x2,r,m) >=0 ) {
    /* only for fit models */

      /* best model */
      if (m==b) printf(">"); else printf(" ");
      /* chi-square */
      if (gsl_matrix_get(A.dof,r,m) > 0) {
	if (gsl_cdf_chisq_Q (gsl_matrix_get(A.x2,r,m),
	  gsl_matrix_get(A.dof,r,m)) >= critx2)
	  printf("+M%1d ",m);/* accepted */
	else
	  printf("-M%1d ",m);/* rejected */
	}// chi-square
      else    
	printf("*M%1d ",m); /* dof = 0 */

      printf("S2 %5.3f %5.3f ",
	gsl_matrix_get(A.vP[_S2s_],r,m)*gsl_matrix_get(A.vP[_S2f_],r,m),
	sqrt(
	SQR(gsl_matrix_get(A.dP[_S2s_],r,m)*gsl_matrix_get(A.vP[_S2f_],r,m))+ 
	SQR(gsl_matrix_get(A.dP[_S2f_],r,m)*gsl_matrix_get(A.vP[_S2s_],r,m))
	));
      printf("S2s %5.3f %5.3f ",
	gsl_matrix_get(A.vP[_S2s_],r,m), gsl_matrix_get(A.dP[_S2s_],r,m));
      printf("S2f %5.3f %5.3f ",
	gsl_matrix_get(A.vP[_S2f_],r,m), gsl_matrix_get(A.dP[_S2f_],r,m));
      printf("te %6.3f %6.3f ",
	gsl_matrix_get(A.vP[_te_],r,m), gsl_matrix_get(A.dP[_te_],r,m));
      printf("Rex %5.2f %5.2f ",
	gsl_matrix_get(A.vP[_Rex_],r,m), gsl_matrix_get(A.dP[_Rex_],r,m));

      // model statistics
      printf("X2 %7.3f dof %1d BIC %7.3f AIC %7.3f\n",
	gsl_matrix_get(A.x2, r,m),
	(int)(gsl_matrix_get(A.dof,r,m)),
	gsl_matrix_get(A.bic,r,m),
	gsl_matrix_get(A.aic,r,m));
      
      } // if only fit models
    } // model

  /* sim average */
  if (MC > 0) {
    printf("# SIM M%1d avg.   ",b);
    printf("S2 %5.3f %5.3f ",
      gsl_matrix_get(A.aP[_S2s_],r,b)*gsl_matrix_get(A.aP[_S2f_],r,b),
	sqrt(
	SQR(gsl_matrix_get(A.dP[_S2s_],r,b)*gsl_matrix_get(A.vP[_S2f_],r,b))+ 
	SQR(gsl_matrix_get(A.dP[_S2f_],r,b)*gsl_matrix_get(A.vP[_S2s_],r,b))
	));
    printf("S2s %5.3f %5.3f ", 
      gsl_matrix_get(A.aP[_S2s_],r,b),gsl_matrix_get(A.dP[_S2s_],r,b));
    printf("S2f %5.3f %5.3f ",
      gsl_matrix_get(A.aP[_S2f_],r,b),gsl_matrix_get(A.dP[_S2f_],r,b));
    printf("te %6.3f %6.3f ",
      gsl_matrix_get(A.aP[_te_],r,b),gsl_matrix_get(A.dP[_te_],r,b));
    printf("Rex %5.2f %5.2f\n",
      gsl_matrix_get(A.aP[_Rex_],r,b),gsl_matrix_get(A.dP[_Rex_],r,b));
    }

  /* experimental & calculated relaxation data */

  gsl_vector *fit = gsl_vector_alloc (A.NF*3);
  gsl_vector *sim = gsl_vector_alloc (A.NF*3);

  for(i=0;i<NP;i++) A.p[i] = gsl_matrix_get(A.vP[i],r,b);
  R1R2NOE(&A,fit,NULL);

  if (MC > 0) {
    for(i=0;i<NP;i++) A.p[i] = gsl_matrix_get(A.aP[i],r,b);
    R1R2NOE(&A,sim,NULL);
    }

  for(f=0;f<A.NF;f++) {
    /* experimental */
    MHz = GSL_MAX (gsl_matrix_get(A.X,r,3*f),gsl_matrix_get(A.X,r,3*f+1));
    MHz = GSL_MAX (MHz, gsl_matrix_get(A.X,r,3*f+2));
    if (MHz > 0) {
      printf("# EXP    %7.2f ",MHz);
      if (gsl_matrix_get(A.X,r,3*f+0)>0)
	printf("%7.4f %7.4f ",
	  gsl_matrix_get(A.Y,r,3*f+0),gsl_matrix_get(A.SIG,r,3*f+0));
      else
	printf("%7s %7s "," "," ");
      if (gsl_matrix_get(A.X,r,3*f+1)>0)
	printf("%7.4f %7.4f ",
	  gsl_matrix_get(A.Y,r,3*f+1),gsl_matrix_get(A.SIG,r,3*f+1));
      else
	printf("%7s %7s "," "," ");

      if (gsl_matrix_get(A.X,r,3*f+2)>0)
	printf("%7.4f %7.4f\n",
	  gsl_matrix_get(A.Y,r,3*f+2),gsl_matrix_get(A.SIG,r,3*f+2));
      else
	printf("%7s %7s\n"," "," ");

      /* fit */
      printf("# FIT M%1d %7.2f ",b,MHz);
      if (gsl_matrix_get(A.X,r,3*f+0)>0)
	printf("%7.4f %7s ",gsl_vector_get(fit,3*f+0)," ");
      else
	printf("%7s %7s "," "," ");
      if (gsl_matrix_get(A.X,r,3*f+1)>0)
	printf("%7.4f %7s ",gsl_vector_get(fit,3*f+1)," ");
      else
	printf("%7s %7s "," "," ");
      if (gsl_matrix_get(A.X,r,3*f+2)>0)
	printf("%7.4f %7s\n",gsl_vector_get(fit,3*f+2)," ");
      else
	printf("%7s %7s\n"," "," ");

      if (MC > 0) {
	/* sim */
	printf("# SIM M%1d %7.2f ",b,MHz);
	if (gsl_matrix_get(A.X,r,3*f+0)>0)
	  printf("%7.4f %7s ",gsl_vector_get(sim,3*f+0)," ");
	else
	  printf("%7s %7s "," "," ");
	if (gsl_matrix_get(A.X,r,3*f+1)>0)
	  printf("%7.4f %7s ",gsl_vector_get(sim,3*f+1)," ");
	else
	  printf("%7s %7s "," "," ");
	if (gsl_matrix_get(A.X,r,3*f+2)>0)
	  printf("%7.4f %7s\n",gsl_vector_get(sim,3*f+2)," ");
	else
	  printf("%7s %7s\n"," "," ");
	}

      }
    }
  gsl_vector_free (fit);
  gsl_vector_free (sim);

  if (MC > 0) {
    /* sim X2 percentile */
    printf("# SIM M%1d  X2 percentile\n",b);
    for(i=0;i<4;i++) {
      printf(_INFO_ "    %4.2f %6.3f   %4.2f %6.3f   ",
      (double)(4*0+i+1)/20.0,A.percentile[4*0+i],
      (double)(4*1+i+1)/20.0,A.percentile[4*1+i]);
      printf("%4.2f %6.3f   %4.2f %6.3f   %4.2f %6.3f\n",
      (double)(4*2+i+1)/20.0,A.percentile[4*2+i],
      (double)(4*3+i+1)/20.0,A.percentile[4*3+i],
      (double)(4*4+i+1)/20.0,A.percentile[4*4+i]);
      }
    }

  /* spectral density function --> XMGR plot */
  if (XMGR_PREFIX) {
    char filename[256];
    sprintf(filename,"%s_%04d.xmgr",XMGR_PREFIX,A.num[r]);
    FILE *xmgr = fopen(filename,"w");
    if (xmgr == NULL) {
      printf(_ERROR_ "cannot open %s\n",filename);
      exit(1);
      }

    rjmap (A, false, xmgr);

    double y;
    double w = 0;
    double wmax = 2.0*M_PI*XMGR_PLOT_W_MAX;
    double wstep = wmax/XMGR_PLOT_W_PTS;
    fprintf(xmgr,"# RESIDUE: %d fitted spectral density function\n",A.num[r]);
    fprintf(xmgr,"# column 1:    w  (1e+9 rad/s)\n");
    fprintf(xmgr,"# column 2:  J(w) (ns/rad)\n");
    fprintf(xmgr,"@type xy\n");
    while (w < wmax) {
      Jw (w, (void *)&A, y, NULL);
      fprintf (xmgr,"%g %g\n",w,y*1e+9);
      w += wstep;
      }
    fprintf (xmgr,"&\n");
    fclose (xmgr);
    }// XMGR_PREFIX
}

void show_info(const string s) 
{
  string margin(_ITEM_WIDTH_-s.length(),' ');
  cout << _INFO_ << s << margin << "= ";
}
