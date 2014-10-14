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
  reduced spectral density mapping according to 
  Farrow et al. (1995) J. Biol. NMR 6, 153-162  
*/

#include "emf.h"

using namespace std;

void rjmap (ALLDATA &A, bool header, FILE *xmgr)
{
  double R1,dR1,R2,dR2,NOE,dNOE;
  double g,v,d,c,fd,fc,B0,wh,wx,k1,k2,Jh,Jn,J0,dJh,dJn,dJ0;
  int f;

  g= gamma_x/gamma_h;
  v= 2.0*M_PI*1e+6;
  k1= SQR(0.870/0.921);
  k2= SQR(0.870/0.955);
  d = 1.0E+3*gamma_h*gamma_x*mu0_h_8pi2/(r_xh*r_xh*r_xh);
  fd  = SQR(d)/4.0;

  if (xmgr) {
    fprintf(xmgr,"# RESIDUE %d reduced spectral densities\n",A.num[A.r]);
    fprintf(xmgr,"# column 1:    w  (1e+9 rad/s)\n");
    fprintf(xmgr,"# column 2:  J(w) (ns/rad)\n");
    fprintf(xmgr,"# column 3: dJ(w) (ns/rad)\n");
    } // xmgr
  else {
    if (header) {
      /* print out header information */ 
      for (f=0;f<A.NF;f++) {
	printf(_INFO_ "B0= %6.2f MHz, ",gsl_matrix_get(A.X,A.r,f*3));
      	printf("0.87H freq.= %6.2f MHz",gsl_matrix_get(A.X,A.r,f*3)*0.87);
	printf(" or %10.3e rad/s\n",gsl_matrix_get(A.X,A.r,f*3)*0.87*v);
	printf(_INFO_ "B0= %6.2f MHz, ",gsl_matrix_get(A.X,A.r,f*3));
	printf("    N freq.= %6.2f MHz",gsl_matrix_get(A.X,A.r,f*3)*g);
	printf(" or %10.3e rad/s\n",gsl_matrix_get(A.X,A.r,f*3)*g*v);
	}
      printf(_INFO_ "J : (ns/rad)\n");
      printf(_INFO_ "\n");
      printf(_INFO_ "RESID   J(0)     dJ(0)    J(N)     dJ(N)  ");
      printf("J(0.87H) dJ(0.87H) B0(MHz)\n");
      printf(_INFO_ "----- -------- -------- -------- -------- ");
      printf("-------- -------- --------\n");
      } // header
    } // else

  for (f=0;f<A.NF;f++) {
    B0 = gsl_matrix_get(A.X,A.r,f*3+0);
    wh = 1.0E-3*2.0*M_PI*B0;
    wx = wh*gamma_x/gamma_h;
    c = 1.0E+3*wx*csa_x/sqrt(3.0);
    fc = SQR(c);
    R1 = gsl_matrix_get(A.Y,A.r,f*3+0);
    R2 = gsl_matrix_get(A.Y,A.r,f*3+1);
    NOE = gsl_matrix_get(A.Y,A.r,f*3+2);
    dR1 = gsl_matrix_get(A.SIG,A.r,f*3+0);
    dR2 = gsl_matrix_get(A.SIG,A.r,f*3+1);
    dNOE= gsl_matrix_get(A.SIG,A.r,f*3+2);

    /* J(0.870wH) */
    Jh= 1/(5*fd)*g*(NOE-1)*R1;
    dJh= sqrt(SQR(1/(5*fd)*g*(1-NOE)*dR1)+SQR(-1/(5*fd)*g*R1*dNOE));

    /* J(wN) */
    Jn= (R1-7*fd*k1*Jh)/((3*fd)+fc);
    dJn= sqrt(SQR(dR1)+SQR(-7*fd*k1*dJh))/((3*fd)+fc);

    /* J(0) */
    J0= (R2-(3*fd/2+fc/2)*Jn-(13*fd/2)*k2*Jh)/(2*fd+2*fc/3);
    dJ0= sqrt(SQR(dR2)+SQR((3*fd/2+fc/2)*dJn)
      +SQR((13*fd/2)*k2*dJh))/(2*fd+2*fc/3);

    if (xmgr) {
      fprintf(xmgr,"@type xydy\n");
      fprintf(xmgr,"# %6.2f MHz\n",B0);
      fprintf(xmgr,"%6.2f %4.2e %4.2e\n",0.,J0*1e+9,dJ0*1e+9);
      fprintf(xmgr,"%6.2f %4.2e %4.2e\n",-wx,Jn*1e+9,dJn*1e+9);
      fprintf(xmgr,"%6.2f %4.2e %4.2e\n",wh,Jh*1e+9,dJh*1e+9);
      fprintf(xmgr,"&\n");
      } // xmgr
    else {
      /* print out */
      printf(" %5d ",A.num[A.r]);
      printf("%4.2e %4.2e ",J0*1e+9,dJ0*1e+9);
      printf("%4.2e %4.2e ",Jn*1e+9,dJn*1e+9);
      printf("%4.2e %4.2e ",Jh*1e+9,dJh*1e+9);
      printf(" %6.2f\n",gsl_matrix_get(A.X,A.r,f*3));
      } // else
    }// f
}
