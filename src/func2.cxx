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
#include <gsl/gsl_deriv.h>
#define CC_DERIVATIVE_h	(1.e-9)

using namespace std;

double Jw_CC_f (double x, void *par)
{
  ALLDATA *A = (ALLDATA *)par;
  double te,S2,t0,e;
  double t,w,f1,f2,f3;
  w  = A->u[0];
  t0 = (A->fpar == _tc_ ? x : A->p[_tc_]); 	// ns
  S2 = A->p[_S2s_];
  te = (A->fpar == _te_ ? x : A->p[_te_]); 	// ns
  e  = (A->fpar == _c_ ? x : A->p[_c_]);
  t = te*t0/(te + t0); // ns
  f1 = M_PI/2.*(1.-e);
  f2 = e*log(w*t0);
  f3 = log(w*t);
  return (1./5.)*1.e-9/w*(S2*cos(f1)/(cosh(f2)+sin(f1))+(1.0-S2)/cosh(f3));
}

void Jw_CC (const double x,void *data,double &y,double *dyda)
{
/*
  x : lamoir frequency, (1e+9 rad/s)
  data : relaxation parameters
  y : resulting spectral density
  dyda : resulting partial derivatives

  Cole-Cole spectral density function  
  Based on the Classical Lipari-Szabo Formalism (MODEL 2 & 4)
  e_e assumed to be 1

  Buevich and Baum (1999) J. Am. Chem. Soc. 121:8671-8672
  Buevich, Shinde, Inouye and Baum (2001) J. Biomolecular NMR 20:233-249
*/

  ALLDATA *A = (ALLDATA *)data;
  double te,S2,t0,e,t,w,wt0;
  double s,c,ch2,sh2,ch3,sh3,g0,g1,g2,f0,f1,f2,f3;
  double fret, abserr;
  int i;

  /* Jw_CC not defined at w = 0 */
  /* w = 0 -> 1 rad/s */
  w = ( fabs(x) < 1.e-9 ? 1.e-9 : fabs(x));

  t0 = A->p[_tc_]; 	// ns
  S2 = A->p[_S2s_];
  te = A->p[_te_]; 	// ns
  e  = A->p[_c_]; 	// width

  t = te*t0/(te + t0); // ns

  g0 = M_PI/2.;
  g1 = log(w*t0);
  f0 = (1./5.)*1.e-9/w;
  f1 = g0*(1.-e);
  f2 = g1*e;
  f3 = log(w*t);
  s = sin(f1);
  c = cos(f1);
  ch2 = cosh(f2); 
  sh2 = sinh(f2); 
  ch3 = cosh(f3);
  sh3 = sinh(f3);

  y = f0*(S2*c/(ch2+s)+(1.0-S2)/ch3);

  if (dyda != NULL) {
    dyda[_tc_] = -e*sh2*c*S2/(t0*SQR(s+ch2))
      -(te-t)*w*sh3*(1.-S2)/(t0*te*w*SQR(ch3));
    dyda[_te_] = - (t0-t)*w*sh3*(1.-S2)/(t0*te*w*SQR(ch3));
    dyda[_S2s_] = (c/(ch2+s)-1./ch3); 
    dyda[_S2f_] = 0.0;
    dyda[_c_] = g0*S2*s/(s+ch2)-c*S2*(g1*sh2-g0*c)/SQR(s+ch2);
    for (i=0;i<NP;i++) 
      dyda[i] *= f0;

    /*
    gsl_function F = {&Jw_CC_f, data};
    A->u[0] = w;
    A->fpar = _tc_;
    gsl_deriv_central (&F,A->p[_tc_],CC_DERIVATIVE_h,&fret,&abserr);
    dyda[_tc_] = fret;
    A->fpar = _te_;
    gsl_deriv_central (&F,A->p[_te_],CC_DERIVATIVE_h,&fret,&abserr);
    dyda[_te_] = fret;
    A->fpar = _c_;
    gsl_deriv_central (&F,A->p[_c_],CC_DERIVATIVE_h,&fret,&abserr);
    dyda[_c_] = fret;
    */


  }
}


/*
  Lorentzian spectral density function 
  Based on the Classical Lipari-Szabo Formalism (MODEL 2 & 4)

  Ochsenbein, Neumann, Guittet and Heijenoort (2000)
  Protein Science 11:957-964
  Modified by Sung-Hun Bae

  J(w) = S2*JLD(w)+(2/5)*(1-S2)*t/(1+w^2*t^2)
  t^-1 = tm^-1 + te^-1
*/

void Jw_Lz (const double x,void *data,double &y,double *dyda)
{
  ALLDATA *A = (ALLDATA *)data;
  double e,t,te,tm,tmax,S2,w,fa,cK,cS,w2,cSw2,cDi,ca,cb,cLn,cc,cd,ct;
  int i;

  w = x*1.0e+9;	// (rad/s)
  tm  = A->p[_tc_]*1.0e-9;	// (s)
  S2  = A->p[_S2f_]; 		// S2f
  te  = A->p[_te_]*1.0e-9;	// (s)
  e   = A->p[_c_];			// gwidth
	
  tmax = 100.0e-9;

  t = te*tm/(te + tm);
  fa = 1+SQR(w*t);

  cK = 1.0/(atan((tmax-tm)/e)+atan(tm/e));
  cS = SQR(e)+SQR(tm);
  w2 = SQR(w);
  cSw2 =  cS*w2;
  cDi = 1.0/(1.0+(-2.0*(SQR(e)-SQR(tm))+SQR(cS)*SQR(w))*SQR(w));
  ca = SQR(e)+SQR(tmax-tm);
  cb = cS*(1.0+w2*SQR(tmax));
  cLn = log(ca/cb);
  cc = SQR(e)*(1.0/cS+1.0/ca);
  cd = (tmax-tm)/ca + SQR(tm)/cS;
  ct = atan(tmax*w);
	
  y = 0.2*S2*cDi*(2.0*tm*(1.0+cSw2)-4.0*cK*e*w*tm*atan(w*tmax)
    +cK*e*(1.0-cSw2)*cLn) +0.4*(1-S2)*t/fa;

  if (dyda != NULL) {	
    for (i=0;i<NP;i++) dyda[i] = 0;

  /* tc */
  dyda[_tc_] = -(4.0*tm*w2*(cSw2+1.0))* ( e*(1.0-cSw2)*cLn*cK
    -4.0*e*tm*w*ct*cK +2.0*tm*(1.0+cSw2))*S2*0.2*SQR(cDi)
    +( -cc*(1.0-cSw2)*cLn*SQR(cK) -2.0*e*tm*w2*cLn*cK
    -4.0*e*w*ct*cK +4.0*tm*cc*w*ct*SQR(cK) +e*(1.0-cSw2)*cb*(
    -2.0*(tmax-tm)/cb -ca*(2.0*tm*SQR(tmax)*w2+2.0*tm)/SQR(cb)
    )*cK/ca +2.0*(cSw2+1.0)+4.0*SQR(tm)*w2)*S2*0.2*cDi
    +0.4*(1.0-S2)*(te-t)/((te+tm)*fa) 
    -0.8*(1.0-S2)*pow(t,3)*w2*(1.0/tm-1.0/(te+tm))/SQR(fa);

  /* S2f */
  dyda[_S2f_] = 0.2*cDi*(2.0*tm*(1.0+cSw2)-4.0*cK*e*w*tm*ct
    +cK*e*(1.0-cSw2)*cLn) -0.4*t/fa;

  /* S2s */
  dyda[_S2s_] = 0.0;

  /* te */
  dyda[_te_] = 0.4*(1.0-S2)*((tm-t)/((te+tm)*fa)
    -2.0*SQR(t)*(1.0-fa)/SQR(te*fa));

  /* e (width of distribution) */
  dyda[_c_] = ( (1.0-cSw2)*cLn*cK +e*cd*(1.0-cSw2)*cLn*SQR(cK)
    -2.0*SQR(e)*w2*cLn*cK -4.0*tm*w*ct*cK -4.0*e*tm*cd*w*ct*SQR(cK)
    +2.0*SQR(e)*(1.0-cSw2)*(1.0/ca-(SQR(tmax)*w2+1.0)/cb)*cK
    +4.0*e*tm*w2)*S2*0.2*cDi -(4.0*e*w2*(cSw2-1.0))* ( e*(1.0-cSw2)*cLn*cK
    -4.0*e*tm*w*ct*cK +2.0*tm*(1.0+cSw2))*S2*0.2*SQR(cDi);
  }
}
