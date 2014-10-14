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

/* 	
  Constants in MKS unit

  Included in the gsl_const_mks.h
  GSL_CONST_MKS_VACUUM_PERMEABILITY (1.25663706144e-6)  (kg m / A^2 s^2)
  GSL_CONST_MKS_PLANCKS_CONSTANT_H (6.62606876e-34)  (kg m^2 / s)
  GSL_CONST_MKS_PLANCKS_CONSTANT_HBAR (1.05457159642e-34) (kg m^2 / s) 

  mu0     : 4*M_PI*1e-7;        # N m^-2 or kg m^-1 s^-2 or T m A^-1
  h	  : 6.62606876e-34;   	# J*s or kg m^2 s^-1
  gamma_h : 26.7522e+7;       	# 1H  rad T^-1 s^-1, T = kg s^-2 A^-1
  gamma_c :  6.7283e+7;       	# 13C rad T^-1 s^-1, or kg^-1 s A
  gamma_n : -2.7126e+7;       	# 15N rad T^-1 s^-1
  rxh     : 1.02e-10;         	# 1H-15N m
  rxh     : 1.09e-10;         	# 1H-13C m
  csa     : -172e-6;          	# 15N, amide
  csa     : 25e-6;            	# 13C, methine alpha-carbon
  d       : gamma_h*gamma_x*r_xh^-3*mu0*h/(8*M_PI^2);
  c       : wx*csa/sqrt(3.0);

  NOTE: for unit analysis, A(ampere) ~ m
  CAUTION: internally, correlation times are treated in ns unit
	and angular frequency is treated in 1e+9 rad/sec.
  */

void R1R2NOE (void *data, gsl_vector *relax, gsl_matrix *jacob)
{
  /* 
  relax[0] - R1 at field 1
  relax[1] - R2 at field 1
  relax[2] - NOE at field 1
  relax[3] - R1 at field 2
  relax[4] - R2 at field 2
  relax[5] - NOE at field 2
  ...
  jacob[0][*] - dR1/da at field 1
  jacob[1][*] - dR2/da at field 1
  jacob[2][*] - dNOE/da at field 1
  jacob[3][*] - dR1/da at field 2
  jacob[4][*] - dR2/da at field 2
  jacob[5][*] - dNOE/da at field 2
  ...
 */

  ALLDATA *A = (ALLDATA *)data;
  double dy0da[NP],dy1da[NP],dy2da[NP],dy3da[NP],dy4da[NP];
  double wh,wx,d,fd_4,fd_8,c,fc,fc_6,fg;
  double y0,y1,y2,y3,y4;
  double MHz,R1,R1i,R2,NOE,Rex,fRex;
  double dR1j,dR2j,dNOEj;
  int i,j,f;

  /* field independent constants */
  d = 1.0E+3*gamma_h*gamma_x*mu0_h_8pi2/(r_xh*r_xh*r_xh);
  fd_4 = d*d/4.0;
  fd_8 = d*d/8.0;
  fg = (gamma_h/gamma_x)*d*d/4.0;

  for (f=0;f<A->NF;f++) {
    for (i=0;i<3;i++) 
      if (gsl_matrix_get(A->X,A->r,3*f+i) > 0) 
	break;
    MHz = gsl_matrix_get(A->X,A->r,3*f+i);
    wh = 1.0E-3*2.0*M_PI*MHz; // 1e+9 (rad/s)
    wx = wh*gamma_x/gamma_h; // 1e+9 (rad/s)
    fc = 1.0E+6*wx*wx*csa_x*csa_x/3.0;
    fc_6 = fc/6.0;
    if (jacob != NULL) {
      Jw (0.0,  data,y0,dy0da);
      Jw (wx,   data,y1,dy1da);
      Jw (wh-wx,data,y2,dy2da);
      Jw (wh,   data,y3,dy3da);
      Jw (wh+wx,data,y4,dy4da);
      }
    else {
      Jw (0.0,  data,y0,NULL);
      Jw (wx,   data,y1,NULL);
      Jw (wh-wx,data,y2,NULL);
      Jw (wh,   data,y3,NULL);
      Jw (wh+wx,data,y4,NULL);
      }

    /* R1 */
    if (gsl_matrix_get(A->X,A->r,3*f) > 0) { /* if R1 is available */
      R1 = fd_4*(y2+3.0*y1+6.0*y4)+fc*y1; /* (1/s) */
      gsl_vector_set (relax,3*f,R1);
      /* R1 jacobian */
      if (jacob != NULL) {
	for(j=0;j<NP;j++) {
	  dR1j = fd_4*(dy2da[j]+3.0*dy1da[j]+6.0*dy4da[j])+fc*dy1da[j];
	  gsl_matrix_set(jacob,3*f+0,j,dR1j);
	  }
	}
      }

    /* R2 */
    if (gsl_matrix_get(A->X,A->r,3*f+1) > 0) { /* if R2 is available */
      R2 = fd_8*(4.0*y0+y2+3.0*y1+6.0*y3+6.0*y4)+fc_6*(4.0*y0+3.0*y1);
      /* quadratically scaled relative to the first R2 field */
      fRex = SQR(gsl_matrix_get(A->X,A->r,3*f+1)/gsl_matrix_get(A->X,A->r,1));
      R2 += fRex*(A->p[_Rex_]);

      gsl_vector_set (relax,3*f+1,R2);

      /* R2 jacobian */
      if (jacob != NULL) {
	for(j=0;j<NP;j++) {
	  if (j == _Rex_) 
	    dR2j = fRex;
	  else
	    dR2j = fd_8*(4.0*dy0da[j]+dy2da[j]+3.0*dy1da[j]+6.0*dy3da[j]
	      +6.0*dy4da[j]) +fc_6*(4.0*dy0da[j]+3.0*dy1da[j]);
	  gsl_matrix_set(jacob,3*f+1,j,dR2j);
	  } // for
	} // if jacob
      }// if R2 is available

    /* NOE */
    if (gsl_matrix_get(A->X,A->r,3*f+2) > 0) { /* if NOE is available */
      NOE = 1.0+fg*(6.0*y4-y2)/(fd_4*(y2+3.0*y1+6.0*y4)+fc*y1);
      gsl_vector_set(relax,3*f+2,NOE);
      /* its jacobian if requested */
      if (jacob != NULL) {
	for (j=0;j<NP;j++) {
	  dR1j = fd_4*(dy2da[j]+3.0*dy1da[j]+6.0*dy4da[j])+fc*dy1da[j];
	  dNOEj = fg*((6.0*dy4da[j]-dy2da[j])*R1i
	    -(6.0*y4-y2)*SQR(R1i)*dR1j);
	  gsl_matrix_set(jacob,3*f+2,j,dNOEj);
	  }
	}
      } // if NOE is available
    }// f
}

void Jw_LSe (const double w,void *data,double &y,double *dyd)
{
  /* 

  w   : lamoir frequency, (1e+9 rad/s)
  data: relaxation parameters
  y   : resulting spectral density
  dyd : resulting partial derivatives

  NOTE: A->r should be set for axial, anisotropic, or distribution 

  Extended Lipari-Szabo formalism 

  Lipari and Szabo (1982) JACS 104, 4546-4559
  Clore et al. (1990) JACS 112, 4989-4991

  J(w) = (2/5)[S2f*S2s*tc/(1+w*w*tc*tc)+
    (1-S2f)*tf'/(1+w*w*tf'*tf')+(S2f-S2)*ts'/(1+w*w*ts'*ts')]
    where,1/tf'=1/tc+1/tf and 1/ts'=1/tc+1/ts

  when tf << ts < tc, simply:
  J(w) = (2/5)[S2f*S2s*tc/(1+w*w*tc*tc)+(S2f-S2)*ts'/(1+w*w*ts'*ts')]
    = (2/5)*S2f*[S2s*tc/(1+w*w*tc*tc)+(1-S2s)*ts'/(1+w*w*ts'*ts')]

  */

  ALLDATA *A = (ALLDATA *)data;
  double fac,p0,p1,p2,wt,tk,t;
  double sum_0,sum_1,sum_2,sum_3,sum_4,sum_w;
  double sum_5,sum_6,sum_7,sum_8,sum_9,sum_A;
  double sum_B,sum_C,sum_D,sum_E;
  int i,k;

  double S2f = A->p[_S2f_];
  double S2s = A->p[_S2s_];
  double te  = A->p[_te_];//(ns) 

  /* for axial */
  size_t r;
  double phi,theta,Dper,Dpar;
  double a,cos_a,cos2_a,sin2_a,cosin_a;
  double dyda,dadphi,dadtheta,dtkdDr,dtkdtc,dwtda;

  /* for bimodal */
  double dtkdtb,dwtdc;

  sum_0=sum_1=sum_2=sum_3=sum_4=sum_w=0.0;
  sum_5=sum_6=sum_7=sum_8=sum_9=sum_A=0.0;
  sum_B=sum_C=sum_D=sum_E=0.0;

  if (D & _AXIAL_) {
    r = A->r;
    phi	  = A->p[_phi_];  /* unit: radian */
    theta = A->p[_theta_];/* unit: radian */
    /* projection onto the pricipal z-axis of tensor */
    cos_a =  A->x[r]*sin(theta)*sin(phi)
	    -A->y[r]*sin(theta)*cos(phi)
	    +A->z[r]*cos(theta);
    A->alpha[r] = a = acos(cos_a);	/* unit: radian */
    cos2_a = cos_a*cos_a;
    sin2_a = 1.0-cos2_a;
    cosin_a= cos_a*sin(a);
    Dper = 1.0/(2*A->p[_tc_]*(2.0+A->p[_Dr_]));
    Dpar = A->p[_Dr_]*Dper;
    if (dyd != NULL) {
      dadphi	= (fabs(a) < 1e-20 ? 0.0 : 
	sin(theta)*(A->x[r]*cos(phi)+A->y[r]*sin(phi))/(-sin(a)));
      dadtheta	= (fabs(a) < 1e-20 ? theta :
	(cos(theta)*(A->x[r]*sin(phi)-A->y[r]*cos(phi))-sin(theta)*(A->z[r]))
	/(-sin(a)));
      }
    }

  for (k=0;k<A->NK;k++) {
    /* diffusion parameter setup */
    if (D & _ISOTROPIC_) {
      tk = A->p[_tc_]; // (ns)
      wt = 1.0;
      }
    if (D & _DISTRIBUTION_ ) {
      tk = A->p[_scf_]*gsl_matrix_get(A->TAU,A->r,k);// (ns)
      wt = gsl_matrix_get(A->WT,A->r,k);
      }
    if (D == _GLOBAL_ANISOTROPIC_) {
      tk = A->t[k];// (ns)
      wt = gsl_matrix_get(A->WT,A->r,k);
      }
    if (D & _AXIAL_) {
      switch (k) {
	case 0:
	  tk = 1.0/(6.0*Dper);
	  wt = SQR(1.5*cos2_a-0.5);
	  dtkdtc = (2.0+A->p[_Dr_])/3.0;
	  dtkdDr = (A->p[_tc_])/3.0;
	  dwtda = -6*(1.5*cos2_a-0.5)*cosin_a;
	  break;
	case 1:
	  tk = 1.0/(5.0*Dper + 1.0*Dpar);
	  wt = 3.0*sin2_a*cos2_a;
	  dtkdtc = (2.0+A->p[_Dr_])/(2.5+0.5*A->p[_Dr_]);
	  dtkdDr = (1.5*A->p[_tc_])/SQR(2.5+0.5*A->p[_Dr_]);
	  dwtda = 6*(cos2_a-sin2_a)*cosin_a;
	  break;
	case 2:
	  tk = 1.0/(2.0*Dper + 4.0*Dpar);
	  wt = 0.75*SQR(sin2_a);
	  dtkdtc = (2.0+A->p[_Dr_])/(1.0+2.0*A->p[_Dr_]);
	  dtkdDr = -(3.0*A->p[_tc_])/SQR(1.0+2.0*A->p[_Dr_]);
	  dwtda = 3*sin2_a*cosin_a;
	  break;
	  } 
	} // if axial

  if (D & _BIMODAL_) {
    switch (k) {
      case 0:
	tk = A->p[_tc_];
	wt = A->p[_c_];
	dtkdtc = 1.0;
	dtkdtb = 0.0;
	dwtdc = 1.0;
	break;
      case 1:
	tk = A->p[_tb_];
	wt = (1.0-A->p[_c_]);
	dtkdtc = 0.0;
	dtkdtb = 1.0;
	dwtdc = -1.0;
	break;
	} 
      } // if bimodal

  if (tk > 0.0 && wt > 0.0) {
    sum_w += wt;
    p0 = 1.0/(1.0+SQR(w*tk)); // unit-less
    p1 = tk*p0;// (ns)
    sum_0 += wt*p1;// (ns)
    if (dyd != NULL) { 
      sum_2 += wt*(p0-2.0*SQR(w*p1));
      sum_5 += p1*dwtdc;
      sum_7 += wt*(p0-2.0*SQR(w*p1))*dtkdtb;
      sum_9 += wt*(p0-2.0*SQR(w*p1))*dtkdtc;
      sum_B += wt*(p0-2.0*SQR(w*p1))*dtkdDr;
      sum_D += p1*dwtda;
      } // if partial derivatives are requested
    if(te > 0.0) {
      t  = te*tk/(te+tk);// (ns)
      p2 = t/(1.0+SQR(w*t));// (ns)
      sum_1 += wt*p2;// (ns)
      if (dyd != NULL) {
	sum_3 += wt*(p2/tk - p2/(te+tk) -2.0*SQR(w*p2*t/tk));
	sum_4 += wt*(p2/te - p2/(te+tk) -2.0*SQR(w*p2*t/te));
	sum_6 += p2*dwtdc;
	sum_8 += wt*(p2/tk - p2/(te+tk) -2.0*SQR(w*p2*t/tk))*dtkdtb;
	sum_A += wt*(p2/tk - p2/(te+tk) -2.0*SQR(w*p2*t/tk))*dtkdtc;
	sum_C += wt*(p2/tk - p2/(te+tk) -2.0*SQR(w*p2*t/tk))*dtkdDr;
	sum_E += p2*dwtda;
	} // if partial derivatives are requested
      } // te > 0
    } // tk > 0 && wt > 0

  } // k

  /* sum of weight should be 1 */
  fac = 1.0e-9*(2.0/5.0)*sum_w;
  y = fac*S2f*(S2s*sum_0 + (1.0-S2s)*sum_1);// (s/rad)

  if(dyd != NULL) { 
    for(i=0;i<NP;i++) dyd[i] = 0.0;
    dyd[_S2f_] = fac*(S2s*sum_0 + (1.0-S2s)*sum_1);
    dyd[_S2s_] = fac*S2f*(sum_0 - sum_1);
    dyd[_te_]  = fac*S2f*(1.0-S2s)*sum_4;// te (ns)
    if (D & _ISOTROPIC_) dyd[_tc_] = fac*S2f*(S2s*sum_2 + (1.0-S2s)*sum_3);
    if (D & _AXIAL_) { 
      dyd[_tc_]	    = fac*S2f*(S2s*sum_9 + (1.0-S2s)*sum_A);
      dyd[_Dr_]	    = fac*S2f*(S2s*sum_B + (1.0-S2s)*sum_C);
      dyda	    = fac*S2f*(S2s*sum_D + (1.0-S2s)*sum_E);
      dyd[_phi_]    = dyda*dadphi;
      dyd[_theta_]  = dyda*dadtheta;
      }
    if (D & _BIMODAL_) {
      dyd[_tc_]  = fac*S2f*(S2s*sum_9 + (1.0-S2s)*sum_A);
      dyd[_tb_]  = fac*S2f*(S2s*sum_7 + (1.0-S2s)*sum_8);
      dyd[_c_]   = fac*S2f*(S2s*sum_5 + (1.0-S2s)*sum_6);
      }
    } // if partial derivatives are requested
}
