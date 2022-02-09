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
  GSL_CONST_MKS_VACUUM_PERMEABILITY   (1.25663706144e-6)  (kg m / A^2 s^2)
  GSL_CONST_MKS_PLANCKS_CONSTANT_H    (6.62606876e-34)    (kg m^2 / s)
  GSL_CONST_MKS_PLANCKS_CONSTANT_HBAR (1.05457159642e-34) (kg m^2 / s) 

  mu0     : 4*M_PI*1e-7;        # N m^-2 or kg m^-1 s^-2 or T m A^-1
  h	      : 6.62606876e-34;   	# J*s or kg m^2 s^-1
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


void R1R2NOE(void *data, gsl_vector *relax, gsl_matrix *jacob)
{
    /* 
    relax[0] - R1  at field 1
    relax[1] - R2  at field 1
    relax[2] - NOE at field 1
    relax[3] - R1  at field 2
    relax[4] - R2  at field 2
    relax[5] - NOE at field 2
    ...
    jacob[0][*] - dR1/da  at field 1
    jacob[1][*] - dR2/da  at field 1
    jacob[2][*] - dNOE/da at field 1
    jacob[3][*] - dR1/da  at field 2
    jacob[4][*] - dR2/da  at field 2
    jacob[5][*] - dNOE/da at field 2
    ...
    */

    ALLDATA *A = (ALLDATA *)data;
    double dy0da[NP], dy1da[NP], dy2da[NP], dy3da[NP], dy4da[NP];
    double wh, wx, d, fd_4, fd_8, c, fc, fc_6, fg;
    double y0, y1, y2, y3, y4;
    double MHz, R1, R1i, R2, NOE, Rex, fRex;
    double dR1j, dR2j, dNOEj;
    int i, j, f;

    /* field independent constants */
    d = 1.0E+3 * gamma_h * gamma_x * mu0_h_8pi2 / (r_xh * r_xh * r_xh);
    fd_4 = d * d / 4.0;
    fd_8 = d * d / 8.0;
    fg = (gamma_h / gamma_x) * d * d / 4.0;

    for (f = 0; f < A->NF; f++)
    {
        for (i = 0; i < 3; i++)
        {
            if (gsl_matrix_get(A->X, A->r, 3 * f + i) > 0)
                break;
        }
        MHz = gsl_matrix_get(A->X, A->r, 3 * f + i);
        wh = 1.0E-3 * 2.0 * M_PI * MHz; // 1e+9 (rad/s)
        wx = wh * gamma_x / gamma_h;    // 1e+9 (rad/s)
        fc = 1.0E+6 * wx * wx * csa_x * csa_x / 3.0;
        fc_6 = fc / 6.0;
        if (jacob != NULL)
        {
            Jw(0.0,     data, y0, dy0da);
            Jw(wx,      data, y1, dy1da);
            Jw(wh - wx, data, y2, dy2da);
            Jw(wh,      data, y3, dy3da);
            Jw(wh + wx, data, y4, dy4da);
        }
        else
        {
            Jw(0.0,     data, y0, NULL);
            Jw(wx,      data, y1, NULL);
            Jw(wh - wx, data, y2, NULL);
            Jw(wh,      data, y3, NULL);
            Jw(wh + wx, data, y4, NULL);
        }

        /* R1 */
        if (gsl_matrix_get(A->X, A->r, 3 * f) > 0) /* if R1 is available */
        {                                                   
            R1 = fd_4 * (y2 + 3.0 * y1 + 6.0 * y4) + fc * y1; /* (1/s) */
            gsl_vector_set(relax, 3 * f, R1);
            /* R1 jacobian */
            if (jacob != NULL)
            {
                for (j = 0; j < NP; j++)
                {
                    dR1j = fd_4 * (dy2da[j] + 3.0 * dy1da[j] + 6.0 * dy4da[j]) + fc * dy1da[j];
                    gsl_matrix_set(jacob, 3 * f + 0, j, dR1j);
                }
            }
        }

        /* R2 */
        if (gsl_matrix_get(A->X, A->r, 3 * f + 1) > 0) /* if R2 is available */
        { 
            R2 = fd_8 * (4.0 * y0 + y2 + 3.0 * y1 + 6.0 * y3 + 6.0 * y4) + fc_6 * (4.0 * y0 + 3.0 * y1);
            /* quadratically scaled relative to the first R2 field */
            fRex = SQR(gsl_matrix_get(A->X, A->r, 3 * f + 1) / gsl_matrix_get(A->X, A->r, 1));
            R2 += fRex * (A->p[_Rex_]);

            gsl_vector_set(relax, 3 * f + 1, R2);

            /* R2 jacobian */
            if (jacob != NULL)
            {
                for (j = 0; j < NP; j++)
                {
                    if (j == _Rex_)
                    {
                        dR2j = fRex;
                    }
                    else
                    {
                        dR2j = fd_8 * (4.0 * dy0da[j] + dy2da[j] + 3.0 * dy1da[j] + 6.0 * dy3da[j] + 6.0 * dy4da[j]) +
                            fc_6 * (4.0 * dy0da[j] + 3.0 * dy1da[j]);
                    }
                    gsl_matrix_set(jacob, 3 * f + 1, j, dR2j);
                } // for
            }   // if jacob
        }     // if R2 is available

        /* NOE */
        if (gsl_matrix_get(A->X, A->r, 3 * f + 2) > 0) /* if NOE is available */
        { 
            NOE = 1.0 + fg * (6.0 * y4 - y2) / (fd_4 * (y2 + 3.0 * y1 + 6.0 * y4) + fc * y1);
            gsl_vector_set(relax, 3 * f + 2, NOE);
            /* its jacobian if requested */
            if (jacob != NULL)
            {
                for (j = 0; j < NP; j++)
                {
                    dR1j = fd_4 * (dy2da[j] + 3.0 * dy1da[j] + 6.0 * dy4da[j]) + fc * dy1da[j];
                    dNOEj = fg * ((6.0 * dy4da[j] - dy2da[j]) * R1i - (6.0 * y4 - y2) * SQR(R1i) * dR1j);
                    gsl_matrix_set(jacob, 3 * f + 2, j, dNOEj);
                }
            }
        } // if NOE is available
    }   // f
}


void Jw_LSe(const double w, void *data, double &fval, double *dyd)
{
    /* 

    Extended Lipari-Szabo formalism
    Lipari and Szabo (1982) JACS 104, 4546-4559
    Clore et al. (1990) JACS 112, 4989-4991 

    J(w) = (2/5)[S2f*S2s*tc/(1+w*w*tc*tc)+
    (1-S2f)*tf'/(1+w*w*tf'*tf')+(S2f-S2)*ts'/(1+w*w*ts'*ts')]
    where,1/tf'=1/tc+1/tf and 1/ts'=1/tc+1/ts

    when tf << ts < tc, simply:
    J(w) = (2/5)[S2f*S2s*tc/(1+w*w*tc*tc)+(S2f-S2)*ts'/(1+w*w*ts'*ts')]
    = (2/5)*S2f*[S2s*tc/(1+w*w*tc*tc)+(1-S2s)*ts'/(1+w*w*ts'*ts')]

    w     : lamoir frequency, (1e+9 rad/s)
    data  : relaxation parameters
    y     : resulting spectral density
    dyd   : resulting partial derivatives

    NOTE: A->r (residue index) should be set for axial, anisotropic, or distribution


    Anisotropic diffusion tensor (Euler angles)
    ===========================================

    according to (phi, theta, psi) convention
    https://mathworld.wolfram.com/EulerAngles.html
    
    A = BCD
    
    1. rotation by phi about z-axis using D matrix
    2. rotation by theta [0, pi] about x-axis using C matrix
    3. rotation by psi about z-axis using B matrix

    A = | a11 a12 a13 |
        | a21 a22 a23 |
        | a31 a32 a33 |

    a_(11)	=	cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi)
    a_(12)	=	cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi)
    a_(13)	=	sin(psi)*sin(theta)
    a_(21)	=	-sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi)
    a_(22)	=	-sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi)
    a_(23)	=	cos(psi)*sin(theta)
    a_(31)	=	sin(theta)*sin(phi)
    a_(32)	=	-sin(theta)*cos(phi)
    a_(33)	=	cos(theta)
    

    Axially symmetric diffusion tensor
    ==================================

    There are two possible solutions: protate and oblate.
    prolate (z-axis is the longest,  Dr >1)
    oblate  (z-axis is the shortest, Dr <1)
    
    approximation of axially symmetric diffusion tensor
    Dzz     --> Dpar
    Dxx=Dyy --> Dper
    Dr = Dpar/Dper = Dzz/Dxx
    
    Dper = 1.0 / (2 * tc * (2.0 + Dr ))

    For isotropic (Dr=1), Dper = Dpar = 1/(6*tc)
    For prolate (Dr>1), Dper < Dpar
    For oblate (Dr<1), Dper > Dpar
    
    A = BCD
    1. rotation by phi about z-axis using D matrix
    2. rotation by theta [0, pi] about x-axis using C matrix
    
    A = | a11 a12 a13 |
        | a21 a22 a23 |
        | a31 a32 a33 |

    a_(11)	=	cos(phi)
    a_(12)	=	sin(phi)
    a_(13)	=	0
    a_(21)	=	-cos(theta)*sin(phi)
    a_(22)	=	+cos(theta)*cos(phi)
    a_(23)	=	sin(theta)
    a_(31)	=	sin(theta)*sin(phi)
    a_(32)	=	-sin(theta)*cos(phi)
    a_(33)	=	cos(theta)

    */

    ALLDATA *A = (ALLDATA *)data;

    double fac, p0, p1, p2, wt, tk, t;
    double sum_wt  = 0.0;
    double sum_0 = 0.0;
    double sum_1 = 0.0;
    double deriv_2 = 0.0;
    double deriv_3 = 0.0;
    double deriv_4 = 0.0;
    double deriv_5 = 0.0;
    double deriv_6 = 0.0;
    double deriv_7 = 0.0;
    double deriv_8 = 0.0;
    double deriv_9 = 0.0;
    double deriv_A = 0.0;
    double deriv_B = 0.0;
    double deriv_C = 0.0;
    double deriv_D = 0.0;
    double deriv_E = 0.0;
    double deriv_F = 0.0;
    double deriv_G = 0.0; /* anisotropic */
    double deriv_H = 0.0; /* anisotropic */
    double deriv_I = 0.0; /* anisotropic */
    double deriv_J = 0.0; /* anisotropic */
    double deriv_K = 0.0; /* anisotropic */
    double deriv_L = 0.0; /* anisotropic */
    double deriv_M = 0.0; /* anisotropic */
    double deriv_N = 0.0; /* anisotropic */
    double deriv_O = 0.0; /* anisotropic */
    double deriv_P = 0.0; /* anisotropic */
    double deriv_Q = 0.0; /* anisotropic */
    double deriv_R = 0.0; /* anisotropic */
    
    double S2f      = A->p[_S2f_];      /* internal motion, unit: ns */
    double S2s      = A->p[_S2s_];      /* internal motion, unit: ns */
    double te       = A->p[_te_];       /* internal motion, unit: ns */
    
    double tc       = A->p[_tc_];       /* isotropic unit: ns */
    double Dr       = A->p[_Dr_];       /* axial */
    double phi      = A->p[_phi_];      /* axial & anisotropic, unit: radian */
    double theta    = A->p[_theta_];    /* axial & anisotropic, unit: radian */
    double psi      = A->p[_psi_];      /* anisotropic, unit: radian */
    double Dxx      = A->p[_Dxx_];      /* anisotropic, unit: 1/ns */
    double Dyy      = A->p[_Dyy_];      /* anisotropic, unit: 1/ns */
    double Dzz      = A->p[_Dzz_];      /* anisotropic, unit: 1/ns */

    /* axial */
    double Dper, Dpar;              
    double a, cos_a, cos2_a, sin2_a, cosin_a;
    double dyda, dadphi, dadtheta, dtkdDr, dtkdtc, dwtda;
    double dtkdDper, dtkdDpar, dDperdtc, dDperdDr, dDpardtc, dDpardDr;

    /* anisotropic */
    double Diso, L2, f1, xi, yi, zi, dx, dy, dz;
    double dtkdDxx, dtkdDyy, dtkdDzz;
    double dxidphi, dyidphi, dzidphi;
    double dxidtheta, dyidtheta, dzidtheta;
    double dxidpsi, dyidpsi, dzidpsi;
    double dwtdxi, dwtdyi, dwtdzi;
    double dwtdphi, dwtdtheta, dwtdpsi;

    /* for bimodal */
    double dtkdtb, dwtdc;

    int i, k;
    double x, y, z;     /* unit bond vector */
    
    /* preprocessing */
    if (D & _ISOTROPIC_)
    {
        // do nothing 
    }

    else if (D & _AXIAL_)
    {
        x = A->x[A->r];
        y = A->y[A->r];
        z = A->z[A->r];
        Dper = 1.0 / (2 * tc * (2.0 + Dr ));
        Dpar = Dr * Dper;
        dDperdtc = -1./(2.*SQR(tc)*(Dr + 2.0));
        dDperdDr = -1./(2.*tc*SQR(Dr + 2.0));
        dDpardtc = -Dr/(2.*SQR(tc)*(Dr + 2.0));
        dDpardDr = 1./(tc*SQR(Dr + 2.0));

        // a_(11)	=	cos(phi)
        // a_(12)	=	sin(phi)
        // a_(13)	=	0
        // a_(21)	=	-cos(theta)*sin(phi)
        // a_(22)	=	+cos(theta)*cos(phi)
        // a_(23)	=	sin(theta)
        // a_(31)	=	sin(theta)*sin(phi)
        // a_(32)	=	-sin(theta)*cos(phi)
        // a_(33)	=	cos(theta)
        // angle relative to the diffusion tensor z-axis 
        // or projection on to the tensor z-axis:
        xi =  cos(phi) * x + sin(phi) * y;
        yi = -cos(theta) * sin(phi) * x + cos(theta) * cos(phi) * y + sin(theta) * z;
        zi =  sin(theta) * sin(phi) * x  -sin(theta) * cos(phi) * y + cos(theta) * z;
        A->xi[A->r] = xi;
        A->yi[A->r] = yi;
        A->zi[A->r] = zi;
        A->alpha[A->r] = a = acos(zi); /* unit: radian */
        cos_a =  sin(theta) * sin(phi) * x  -sin(theta) * cos(phi) * y + cos(theta) * z;    
        cos2_a = cos_a * cos_a;
        sin2_a = 1.0 - cos2_a;
        cosin_a = cos_a * sin(a);
        // dxidphi   = -x*sin(phi) + y*cos(phi);
        // dxidtheta = 0.;
        // dyidphi   = -x*cos(phi)*cos(theta) - y*sin(phi)*cos(theta);
        // dyidtheta = (x*sin(phi)*sin(theta) - y*sin(theta)*cos(phi) + z*cos(theta));
        // dzidphi   = (x*cos(phi) + y*sin(phi))*sin(theta);
        // dzidtheta = (x*sin(phi) - y*cos(phi))*cos(theta) - z*sin(theta);
        dadphi   = (x*cos(phi) + y*sin(phi))*sin(theta);
        dadtheta = (x*sin(phi) - y*cos(phi))*cos(theta) - z*sin(theta);
    }

    else if (D & _ANISOTROPIC_)
    {
        x = A->x[A->r];
        y = A->y[A->r];
        z = A->z[A->r];
        Diso = (Dxx + Dyy + Dzz) / 3.0;
        L2 = (Dxx * Dyy + Dxx * Dzz + Dyy * Dzz) / 3.0;
        f1 = sqrt(Diso * Diso - L2);
        if (fabs(f1) < 1.0e-20) 
        {
            dx = 0.0;
            dy = 0.0;
            dz = 0.0;
        }
        else 
        {
            dx = (Dxx - Diso) / f1;
            dy = (Dyy - Diso) / f1;
            dz = (Dzz - Diso) / f1;
        }
        // direction cosines (xi*xi + yi*yi + zi*zi = 1)
        xi = (cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi)) * x +
             (cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi)) * y +
             (sin(psi)*sin(theta)) * z;
        yi = (-sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi)) * x +
             (-sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi)) * y +
             ( cos(psi)*sin(theta)) * z;
        zi = ( sin(theta)*sin(phi)) * x +
             (-sin(theta)*cos(phi)) * y +
             ( cos(theta)) * z;
        A->xi[A->r] = xi;
        A->yi[A->r] = yi;
        A->zi[A->r] = zi;
        // printf("direction cosine length= %.5f\n", (SQR(xi)+ SQR(yi) + SQR(zi)));
        dxidphi   = -x*sin(phi)*cos(psi) - x*sin(psi)*cos(phi)*cos(theta) - y*sin(phi)*sin(psi)*cos(theta) + y*cos(phi)*cos(psi);
        dxidtheta = (x*sin(phi)*sin(theta) - y*sin(theta)*cos(phi) + z*cos(theta))*sin(psi);
        dxidpsi   = -x*sin(phi)*cos(psi)*cos(theta) - x*sin(psi)*cos(phi) - y*sin(phi)*sin(psi) + y*cos(phi)*cos(psi)*cos(theta) + z*sin(theta)*cos(psi);
        dyidphi   =  x*sin(phi)*sin(psi) - x*cos(phi)*cos(psi)*cos(theta) - y*sin(phi)*cos(psi)*cos(theta) - y*sin(psi)*cos(phi);
        dyidtheta = (x*sin(phi)*sin(theta) - y*sin(theta)*cos(phi) + z*cos(theta))*cos(psi);
        dyidpsi   =  x*sin(phi)*sin(psi)*cos(theta) - x*cos(phi)*cos(psi) - y*sin(phi)*cos(psi) - y*sin(psi)*cos(phi)*cos(theta) - z*sin(psi)*sin(theta);
        dzidphi   = (x*cos(phi) + y*sin(phi))*sin(theta);
        dzidtheta =  x*sin(phi)*cos(theta) - y*cos(phi)*cos(theta) - z*sin(theta);
        dzidpsi   = 0.0;
    }

    /* 
        NK=1 for isotropic
        NK=3 for axially symmetric
        NK=5 for anisotropic
    */
    for (k = 0; k < A->NK; k++)
    {
        if (D & _ISOTROPIC_)
        {
            tk = tc;
            wt = 1.0;
        }

        else if (D & _AXIAL_)
        {
            switch (k)
            {
                case 0:
                    tk = 1.0 / (6.0 * Dper);
                    // dtkdDper = -1/(6*SQR(Dper));
                    // dtkdtc = dtkdDper * dDperdtc;
                    // dtkdDr = dtkdDper * dDperdDr;
                    // wt  = 3.0*SQR(yi)*SQR(zi);
                    // wt += 3.0*SQR(xi)*SQR(zi);
                    // dwtdphi    = 6.0 * ( yi * SQR(zi) * dyidphi    + SQR(yi) * zi * dzidphi   );
                    // dwtdtheta  = 6.0 * ( yi * SQR(zi) * dyidtheta  + SQR(yi) * zi * dzidtheta );
                    // dwtdphi   += 6.0 * ( xi * SQR(zi) * dxidphi    + SQR(xi) * zi * dzidphi   );
                    // dwtdtheta += 6.0 * ( xi * SQR(zi) * dxidtheta  + SQR(xi) * zi * dzidtheta );

                    wt = SQR(1.5 * cos2_a - 0.5);
                    dtkdtc = (2.0 + Dr) / 3.0;
                    dtkdDr = tc / 3.0;
                    dwtda = -6. * (1.5 * cos2_a - 0.5) * cosin_a;
                    break;
                case 1:
                    tk = 1.0 / (5.0 * Dper + 1.0 * Dpar);
                    // dtkdDper = -5./SQR(Dpar + 5.*Dper);
                    // dtkdDpar = -1./SQR(Dpar + 5.*Dper);
                    // dtkdtc = dtkdDper * dDperdtc + dtkdDpar * dDpardtc;
                    // dtkdDr = dtkdDper * dDperdDr + dtkdDpar * dDpardDr;
                    // wt  = 3.0*SQR(xi)*SQR(yi);
                    // wt += 0.25*(3.0*(SQR(xi)*SQR(xi) + SQR(yi)*SQR(yi) +SQR(zi)*SQR(zi))-1.0) - (1.0/12.0)*(
                    //         dx*(3.0*SQR(xi)*SQR(xi)+6.0*SQR(yi)*SQR(zi)-1.0) +
                    //         dy*(3.0*SQR(yi)*SQR(yi)+6.0*SQR(zi)*SQR(xi)-1.0) +
                    //         dz*(3.0*SQR(zi)*SQR(zi)+6.0*SQR(xi)*SQR(yi)-1.0));
                    // dwtdphi   = 6.0 * ( xi * SQR(yi) * dxidphi    + SQR(xi) * yi * dyidphi   );
                    // dwtdtheta = 6.0 * ( xi * SQR(yi) * dxidtheta  + SQR(xi) * yi * dyidtheta );
                    // dwtdxi = -xi*(dx*SQR(xi) + dy*SQR(zi) + dz*SQR(yi) - 3.0*SQR(xi));
                    // dwtdyi = -yi*(dx*SQR(zi) + dy*SQR(yi) + dz*SQR(xi) - 3.0*SQR(yi));
                    // dwtdzi = -zi*(dx*SQR(yi) + dy*SQR(xi) + dz*SQR(zi) - 3.0*SQR(zi));
                    // dwtdphi   += dwtdxi*dxidphi   + dwtdyi*dyidphi   + dwtdzi*dzidphi;
                    // dwtdtheta += dwtdxi*dxidtheta + dwtdyi*dyidtheta + dwtdzi*dzidtheta;

                    wt = 3.0 * sin2_a * cos2_a;
                    dtkdtc = (2.0 + Dr) / (2.5 + 0.5 * Dr);
                    dtkdDr = (1.5 * tc) / SQR(2.5 + 0.5 * Dr);
                    dwtda = 6. * (cos2_a - sin2_a) * cosin_a;
                    break;
                case 2:
                    tk = 1.0 / (2.0 * Dper + 4.0 * Dpar);
                    // dtkdDper = -1./(2.*SQR(2.*Dpar + Dper));
                    // dtkdDpar = -1./SQR(2.*Dpar + Dper);
                    // dtkdtc = dtkdDper * dDperdtc + dtkdDpar * dDpardtc;
                    // dtkdDr = dtkdDper * dDperdDr + dtkdDpar * dDpardDr;
                    // wt = 0.25*(3.0*(SQR(xi)*SQR(xi) + SQR(yi)*SQR(yi) +SQR(zi)*SQR(zi))-1.0) + (1.0/12.0)*(
                    //         dx*(3.0*SQR(xi)*SQR(xi)+6.0*SQR(yi)*SQR(zi)-1.0) +
                    //         dy*(3.0*SQR(yi)*SQR(yi)+6.0*SQR(zi)*SQR(xi)-1.0) +
                    //         dz*(3.0*SQR(zi)*SQR(zi)+6.0*SQR(xi)*SQR(yi)-1.0));
                    // dwtdxi = xi*(dx*SQR(xi) + dy*SQR(zi) + dz*SQR(yi) + 3.0*SQR(xi));
                    // dwtdyi = yi*(dx*SQR(zi) + dy*SQR(yi) + dz*SQR(xi) + 3.0*SQR(yi));
                    // dwtdzi = zi*(dx*SQR(yi) + dy*SQR(xi) + dz*SQR(zi) + 3.0*SQR(zi));
                    // dwtdphi   = dwtdxi*dxidphi   + dwtdyi*dyidphi   + dwtdzi*dzidphi;
                    // dwtdtheta = dwtdxi*dxidtheta + dwtdyi*dyidtheta + dwtdzi*dzidtheta;

                    wt = 0.75 * SQR(sin2_a);
                    dtkdtc = (2.0 + Dr) / (1.0 + 2.0 * Dr);                  
                    dtkdDr = -(3.0 * tc) / SQR(1.0 + 2.0 * Dr); 
                    dwtda = 3. * sin2_a * cosin_a;
                    break;
            }
        }

        else if (D & _ANISOTROPIC_)
        {
            switch (k)
            {
                case 0:
                    tk = 1.0 / (4 * Dxx + Dyy + Dzz);
                    wt = 3.0*SQR(yi)*SQR(zi);
                    dtkdDxx= -4.0 * SQR(tk);
                    dtkdDyy= -1.0 * SQR(tk);
                    dtkdDzz= -1.0 * SQR(tk);
                    dwtdphi   = 6.0 * ( yi * SQR(zi) * dyidphi    + SQR(yi) * zi * dzidphi   );
                    dwtdtheta = 6.0 * ( yi * SQR(zi) * dyidtheta  + SQR(yi) * zi * dzidtheta );
                    dwtdpsi   = 6.0 * ( yi * SQR(zi) * dyidpsi    + SQR(yi) * zi * dzidpsi   );
                    break;
                case 1:
                    tk = 1.0 / (Dxx + 4 * Dyy + Dzz);
                    wt = 3.0*SQR(xi)*SQR(zi);
                    dtkdDxx = -1.0 * SQR(tk);
                    dtkdDyy = -4.0 * SQR(tk);
                    dtkdDzz = -1.0 * SQR(tk);
                    dwtdphi   = 6.0 * ( xi * SQR(zi) * dxidphi    + SQR(xi) * zi * dzidphi   );
                    dwtdtheta = 6.0 * ( xi * SQR(zi) * dxidtheta  + SQR(xi) * zi * dzidtheta );
                    dwtdpsi   = 6.0 * ( xi * SQR(zi) * dxidpsi    + SQR(xi) * zi * dzidpsi   );
                    break;
                case 2:
                    tk = 1.0 / (Dxx + Dyy + 4 * Dzz);
                    wt = 3.0*SQR(xi)*SQR(yi);
                    dtkdDxx = -1.0 * SQR(tk);
                    dtkdDyy = -1.0 * SQR(tk);
                    dtkdDzz = -4.0 * SQR(tk);
                    dwtdphi   = 6.0 * ( xi * SQR(yi) * dxidphi    + SQR(xi) * yi * dyidphi   );
                    dwtdtheta = 6.0 * ( xi * SQR(yi) * dxidtheta  + SQR(xi) * yi * dyidtheta );
                    dwtdpsi   = 6.0 * ( xi * SQR(yi) * dxidpsi    + SQR(xi) * yi * dyidpsi   );
                    break;
                case 3:
                    tk = 1.0 / (6 * Diso + 6 * f1);
                    wt = 0.25*(3.0*(SQR(xi)*SQR(xi) + SQR(yi)*SQR(yi) +SQR(zi)*SQR(zi))-1.0) - (1.0/12.0)*(
                            dx*(3.0*SQR(xi)*SQR(xi)+6.0*SQR(yi)*SQR(zi)-1.0) +
                            dy*(3.0*SQR(yi)*SQR(yi)+6.0*SQR(zi)*SQR(xi)-1.0) +
                            dz*(3.0*SQR(zi)*SQR(zi)+6.0*SQR(xi)*SQR(yi)-1.0));
                    dtkdDxx = ((-2.* Dxx +     Dyy +     Dzz)/(3.*f1) - 2)/SQR(6*Diso + 6*f1);
                    dtkdDyy = ((     Dxx - 2.* Dyy +     Dzz)/(3.*f1) - 2)/SQR(6*Diso + 6*f1);
                    dtkdDzz = ((     Dxx +     Dyy - 2.* Dzz)/(3.*f1) - 2)/SQR(6*Diso + 6*f1);
                    dwtdxi = -xi*(dx*SQR(xi) + dy*SQR(zi) + dz*SQR(yi) - 3.0*SQR(xi));
                    dwtdyi = -yi*(dx*SQR(zi) + dy*SQR(yi) + dz*SQR(xi) - 3.0*SQR(yi));
                    dwtdzi = -zi*(dx*SQR(yi) + dy*SQR(xi) + dz*SQR(zi) - 3.0*SQR(zi));
                    dwtdphi   = dwtdxi*dxidphi   + dwtdyi*dyidphi   + dwtdzi*dzidphi;
                    dwtdtheta = dwtdxi*dxidtheta + dwtdyi*dyidtheta + dwtdzi*dzidtheta;
                    dwtdpsi   = dwtdxi*dxidpsi   + dwtdyi*dyidpsi   + dwtdzi*dzidpsi;
                    break;
                case 4:
                    tk = 1.0 / (6 * Diso - 6 * f1);
                    wt = 0.25*(3.0*(SQR(xi)*SQR(xi) + SQR(yi)*SQR(yi) +SQR(zi)*SQR(zi))-1.0) + (1.0/12.0)*(
                            dx*(3.0*SQR(xi)*SQR(xi)+6.0*SQR(yi)*SQR(zi)-1.0) +
                            dy*(3.0*SQR(yi)*SQR(yi)+6.0*SQR(zi)*SQR(xi)-1.0) +
                            dz*(3.0*SQR(zi)*SQR(zi)+6.0*SQR(xi)*SQR(yi)-1.0));
                    dtkdDxx = (( 2.* Dxx -     Dyy -     Dzz)/(3.*f1) - 2)/SQR(6*Diso - 6*f1);
                    dtkdDyy = ((    -Dxx + 2.* Dyy -     Dzz)/(3.*f1) - 2)/SQR(6*Diso - 6*f1);
                    dtkdDzz = ((    -Dxx -     Dyy + 2.* Dzz)/(3.*f1) - 2)/SQR(6*Diso - 6*f1);
                    dwtdxi = xi*(dx*SQR(xi) + dy*SQR(zi) + dz*SQR(yi) + 3.0*SQR(xi));
                    dwtdyi = yi*(dx*SQR(zi) + dy*SQR(yi) + dz*SQR(xi) + 3.0*SQR(yi));
                    dwtdzi = zi*(dx*SQR(yi) + dy*SQR(xi) + dz*SQR(zi) + 3.0*SQR(zi));
                    dwtdphi   = dwtdxi*dxidphi   + dwtdyi*dyidphi   + dwtdzi*dzidphi;
                    dwtdtheta = dwtdxi*dxidtheta + dwtdyi*dyidtheta + dwtdzi*dzidtheta;
                    dwtdpsi   = dwtdxi*dxidpsi   + dwtdyi*dyidpsi   + dwtdzi*dzidpsi;
                    break;
            }
        }
    
        else if (D & _BIMODAL_)
        {
            switch (k)
            {
                case 0:
                    tk = tc;
                    wt = A->p[_c_];
                    dtkdtc = 1.0;
                    dtkdtb = 0.0;
                    dwtdc = 1.0;
                    break;
                case 1:
                    tk = A->p[_tb_];
                    wt = (1.0 - A->p[_c_]);
                    dtkdtc = 0.0;
                    dtkdtb = 1.0;
                    dwtdc = -1.0;
                    break;
            }
        }

        else if (D & _DISTRIBUTION_)
        {
            tk = A->p[_scf_] * gsl_matrix_get(A->TAU, A->r, k); // (ns)
            wt = gsl_matrix_get(A->WT, A->r, k);
        }

        if (tk > 0.0 && wt > 0.0)
        {
            sum_wt += wt;
            p0 = 1.0 / (1.0 + SQR(w * tk));     // unit-less
            p1 = tk * p0;                       // (ns)
            sum_0 += wt * p1;                   // (ns)
            if (dyd != NULL)
            {
                deriv_2 += wt * (p0 - 2.0 * SQR(w * p1));
                deriv_5 += p1 * dwtdc;
                deriv_7 += wt * (p0 - 2.0 * SQR(w * p1)) * dtkdtb;
                deriv_9 += wt * (p0 - 2.0 * SQR(w * p1)) * dtkdtc;
                deriv_B += wt * (p0 - 2.0 * SQR(w * p1)) * dtkdDr;
                deriv_D += p1 * dwtda;
                deriv_G += wt * (p0 - 2.0 * SQR(w * p1)) * dtkdDxx;
                deriv_H += wt * (p0 - 2.0 * SQR(w * p1)) * dtkdDyy;
                deriv_I += wt * (p0 - 2.0 * SQR(w * p1)) * dtkdDzz;
                deriv_J += p1 * dwtdphi;
                deriv_K += p1 * dwtdtheta;
                deriv_L += p1 * dwtdpsi;
            } // if partial derivatives are requested
            if (te > 0.0)
            {
                t = te * tk / (te + tk);        // (ns)
                p2 = t / (1.0 + SQR(w * t));    // (ns)
                sum_1 += wt * p2;               // (ns)
                if (dyd != NULL)
                {
                    deriv_3 += wt * (p2 / tk - p2 / (te + tk) - 2.0 * SQR(w * p2 * t / tk));
                    deriv_4 += wt * (p2 / te - p2 / (te + tk) - 2.0 * SQR(w * p2 * t / te));
                    deriv_6 += p2 * dwtdc;
                    deriv_8 += wt * (p2 / tk - p2 / (te + tk) - 2.0 * SQR(w * p2 * t / tk)) * dtkdtb;
                    deriv_A += wt * (p2 / tk - p2 / (te + tk) - 2.0 * SQR(w * p2 * t / tk)) * dtkdtc;
                    deriv_C += wt * (p2 / tk - p2 / (te + tk) - 2.0 * SQR(w * p2 * t / tk)) * dtkdDr;
                    deriv_E += p2 * dwtda;
                    deriv_M += wt * (p2 / tk - p2 / (te + tk) - 2.0 * SQR(w * p2 * t / tk)) * dtkdDxx;
                    deriv_N += wt * (p2 / tk - p2 / (te + tk) - 2.0 * SQR(w * p2 * t / tk)) * dtkdDyy;
                    deriv_O += wt * (p2 / tk - p2 / (te + tk) - 2.0 * SQR(w * p2 * t / tk)) * dtkdDzz;
                    deriv_P += p2 * dwtdphi;
                    deriv_Q += p2 * dwtdtheta;
                    deriv_R += p2 * dwtdpsi;
                } // if partial derivatives are requested
            } // te > 0
        }  // tk > 0 && wt > 0
    } // k

    /* sum of weight should be 1 */
    fac = 1.0e-9 * (2.0 / 5.0) * sum_wt;
    fval = fac * S2f * (S2s * sum_0 + (1.0 - S2s) * sum_1); // (s/rad)

    /* derivative */
    if (dyd != NULL)
    {

        for (i = 0; i < NP; i++)
        {
            dyd[i] = 0.0;
        }

        dyd[_S2f_] = fac * (S2s * sum_0 + (1.0 - S2s) * sum_1);
        dyd[_S2s_] = fac * S2f * (sum_0 - sum_1);
        dyd[_te_]  = fac * S2f * (1.0 - S2s) * deriv_4; // te (ns)

        if (D & _ISOTROPIC_)
        {
            dyd[_tc_] = fac * S2f * (S2s * deriv_2 + (1.0 - S2s) * deriv_3);
        }

        else if (D & _AXIAL_)
        {
            dyd[_tc_] = fac * S2f * (S2s * deriv_9 + (1.0 - S2s) * deriv_A);
            dyd[_Dr_] = fac * S2f * (S2s * deriv_B + (1.0 - S2s) * deriv_C);
            dyda = fac * S2f * (S2s * deriv_D + (1.0 - S2s) * deriv_E);
            // if (fabs(a) < 1.0e-20) 
            // {
            //     dadphi = 0.0;
            //     dadtheta = theta;
            // }
            // else 
            // {
            //     dadphi   = sin(theta) * (x * cos(phi) + y * sin(phi)) / (-sin(a));
            //     dadtheta = cos(theta) * (x * sin(phi) - y * cos(phi)) - sin(theta) * z / (-sin(a));
            // }
            dyd[_phi_]   = dyda * dadphi;
            dyd[_theta_] = dyda * dadtheta;
            // dyd[_phi_]   = fac * S2f * (S2s * deriv_J + (1.0 - S2s) * deriv_P);
            // dyd[_theta_] = fac * S2f * (S2s * deriv_K + (1.0 - S2s) * deriv_Q);
        }

        else if (D & _ANISOTROPIC_)
        {
            dyd[_Dxx_]   = fac * S2f * (S2s * deriv_G + (1.0 - S2s) * deriv_M);
            dyd[_Dyy_]   = fac * S2f * (S2s * deriv_H + (1.0 - S2s) * deriv_N);
            dyd[_Dzz_]   = fac * S2f * (S2s * deriv_I + (1.0 - S2s) * deriv_O);
            dyd[_phi_]   = fac * S2f * (S2s * deriv_J + (1.0 - S2s) * deriv_P);
            dyd[_theta_] = fac * S2f * (S2s * deriv_K + (1.0 - S2s) * deriv_Q);
            dyd[_psi_]   = fac * S2f * (S2s * deriv_L + (1.0 - S2s) * deriv_R);
        }

        else if (D & _BIMODAL_)
        {
            dyd[_tc_] = fac * S2f * (S2s * deriv_9 + (1.0 - S2s) * deriv_A);
            dyd[_tb_] = fac * S2f * (S2s * deriv_7 + (1.0 - S2s) * deriv_8);
            dyd[_c_]  = fac * S2f * (S2s * deriv_5 + (1.0 - S2s) * deriv_6);
        }
    } // if partial derivatives are requested
}


void test_anisotropic ()
{
    int NF = 2;// number of magnetic fields
    int NR = 6;// number of residues
    int NK = 5;// number of tau points

    cout << "Testing anisotropic fitting" << endl;
    ALLDATA A;
    Jw = &Jw_LSe;
    initialize (A, NM, NF, NR, NK, MC);
    D = _GLOBAL_ANISOTROPIC_;
    setup_attribute (A);
    
    // internal dynamics
    A.p[_S2s_] = 1.00;
    A.p[_S2f_] = 0.90;
    A.p[_te_ ] = 0.05;
    A.p[_Rex_] = 0.55;
    
    // diffusion (Dxx <= Dyy <= Dzz)
    A.p[_Dxx_] = 0.009;
    A.p[_Dyy_] = 0.017;
    A.p[_Dzz_] = 0.025;
    A.p[_phi_  ] = (M_PI/180.) * 32.0;
    A.p[_theta_] = (M_PI/180.) * 46.0;
    A.p[_psi_  ] = (M_PI/180.) * 40.0;

    gsl_vector * relax = gsl_vector_calloc (3*A.NF);
    gsl_vector * error = gsl_vector_calloc (NP);

    double MHz, R1, R2, NOE, chisq;
    for (int r=0; r < NR; r++)
    {
        for (int f=0; f < NF; f++)
        {
            if (f == 0) 
            { 
                MHz = 500.38; 
            }
            else 
            {
                MHz = 600.13;
            }
            gsl_matrix_set(A.X,   r, 3*f+0, MHz);
            gsl_matrix_set(A.Y,   r, 3*f+0 , 1.0); // R1
            gsl_matrix_set(A.X,   r, 3*f+1, MHz);
            gsl_matrix_set(A.Y,   r, 3*f+1 , 1.0); // R2
            gsl_matrix_set(A.X,   r, 3*f+2, MHz);
            gsl_matrix_set(A.Y,   r, 3*f+2 , 1.0); // NOE
            gsl_matrix_set(A.SIG, r, 3*f+0 , 0.001); // dR1
            gsl_matrix_set(A.SIG, r, 3*f+1 , 0.001); // dR2
            gsl_matrix_set(A.SIG, r, 3*f+2 , 0.001); // dNOE
        }
        A.r = r;
        double phi   = gsl_ran_flat(rng, 0, 2.*M_PI);
        double theta = gsl_ran_flat(rng, 0, M_PI);
        A.x[r] = cos(phi)*sin(theta);
        A.y[r] = sin(phi)*sin(theta);
        A.z[r] = cos(theta);
        A.ivec[r] = true;

        R1R2NOE((void *)&A, relax, NULL);
        for (int f=0; f < NF; f++) 
        {
            MHz = gsl_matrix_get(A.X, r, 3*f+0);
            R1  = gsl_vector_get(relax, 3*f+0);
            R2  = gsl_vector_get(relax, 3*f+1);
            NOE = gsl_vector_get(relax, 3*f+2);
            printf("    %6.2f %8.3f %8.3f %8.3f     XH= %7.4f %7.4f %7.4f\n", 
                MHz, R1, R2, NOE, A.x[r], A.y[r], A.z[r]);
            gsl_matrix_set(A.Y,   r, 3*f+0 , R1);
            gsl_matrix_set(A.Y,   r, 3*f+1 , R2);
            gsl_matrix_set(A.Y,   r, 3*f+2 , NOE);
        }
        
    }

    printf("    Truth   Dxx %7.4f Dyy %7.4f Dzz %7.4f phi %9.4f theta %9.4f psi %9.4f\n",
        A.p[_Dxx_], A.p[_Dyy_], A.p[_Dzz_], A.p[_phi_]*(180./M_PI), A.p[_theta_]*(180./M_PI), A.p[_psi_]*(180./M_PI));
    
    for (int p=0; p < NP; p++)
    {
        A.is[p] = false;
    }
    
    A.is[_Dxx_] = true;
    A.is[_Dyy_] = true;
    A.is[_Dzz_] = true;
    // A.is[_phi_] = true;
    // A.is[_theta_] = true;
    // A.is[_psi_] = true;

    A.p[_Dxx_] = 0.002;
    A.p[_Dyy_] = 0.004;
    A.p[_Dzz_] = 0.006;
    // A.p[_phi_  ] = (M_PI/180.) * 1.0;
    // A.p[_theta_] = (M_PI/180.) * 1.0;
    // A.p[_psi_  ] = (M_PI/180.) * 1.0;
    
    printf("    Seed    Dxx %7.4f Dyy %7.4f Dzz %7.4f phi %9.4f theta %9.4f psi %9.4f\n",
        A.p[_Dxx_], A.p[_Dyy_], A.p[_Dzz_], A.p[_phi_]*(180./M_PI), A.p[_theta_]*(180./M_PI), A.p[_psi_]*(180./M_PI));
    levmar (A, chisq, error);
    printf("    Final   Dxx %7.4f Dyy %7.4f Dzz %7.4f phi %9.4f theta %9.4f psi %9.4f\n\n",
        A.p[_Dxx_], A.p[_Dyy_], A.p[_Dzz_], A.p[_phi_]*(180./M_PI), A.p[_theta_]*(180./M_PI), A.p[_psi_]*(180./M_PI));
    printf("            chisq   %g\n", chisq);
    printf("            Dxx     %7.4f +/- %7.5f\n", A.p[_Dxx_], gsl_vector_get(error, _Dxx_));
    printf("            Dyy     %7.4f +/- %7.5f\n", A.p[_Dyy_], gsl_vector_get(error, _Dyy_));
    printf("            Dzz     %7.4f +/- %7.5f\n", A.p[_Dzz_], gsl_vector_get(error, _Dzz_));
    printf("            phi     %7.4f +/- %7.5f\n", (180./M_PI)*A.p[_phi_],     (180./M_PI)*gsl_vector_get(error, _phi_));
    printf("            theta   %7.4f +/- %7.5f\n", (180./M_PI)*A.p[_theta_],   (180./M_PI)*gsl_vector_get(error, _theta_));
    printf("            psi     %7.4f +/- %7.5f\n", (180./M_PI)*A.p[_psi_],     (180./M_PI)*gsl_vector_get(error, _psi_));
    printf("\n");

    gsl_vector_free(relax);
    gsl_vector_free(error);
    freeing (A); // freeing memory
}


void test_axial ()
{
    int NF = 2;// number of magnetic fields
    int NR = 6;// number of residues
    int NK = 3;// number of tau points

    cout << "Testing axial fitting" << endl;
    ALLDATA A;
    Jw = &Jw_LSe;
    initialize (A, NM, NF, NR, NK, MC);
    D = _GLOBAL_AXIAL_;
    setup_attribute (A);
    
    // internal dynamics
    A.p[_S2s_] = 1.00;
    A.p[_S2f_] = 0.90;
    A.p[_te_ ] = 0.05;
    A.p[_Rex_] = 0.55;
    
    // diffusion
    A.p[_tc_ ] = 9.65;
    A.p[_Dr_ ] = 1.35;
    A.p[_phi_  ] = (M_PI/180.) * 32.0;
    A.p[_theta_] = (M_PI/180.) * 46.0;
 
    gsl_vector * relax = gsl_vector_calloc (3*A.NF);
    gsl_vector * error = gsl_vector_calloc (NP);

    double MHz, R1, R2, NOE, chisq;
    
    for (int r=0; r < NR; r++)
    {
        for (int f=0; f < NF; f++)
        {
            if (f == 0) 
            { 
                MHz = 500.38; 
            }
            else 
            {
                MHz = 600.13;
            }
            gsl_matrix_set(A.X,   r, 3*f+0, MHz);
            gsl_matrix_set(A.Y,   r, 3*f+0 , 1.0); // R1
            gsl_matrix_set(A.X,   r, 3*f+1, MHz);
            gsl_matrix_set(A.Y,   r, 3*f+1 , 1.0); // R2
            gsl_matrix_set(A.X,   r, 3*f+2, MHz);
            gsl_matrix_set(A.Y,   r, 3*f+2 , 1.0); // NOE
            gsl_matrix_set(A.SIG, r, 3*f+0 , 0.001); // dR1
            gsl_matrix_set(A.SIG, r, 3*f+1 , 0.001); // dR2
            gsl_matrix_set(A.SIG, r, 3*f+2 , 0.001); // dNOE
        }
        A.r = r;
        double phi   = gsl_ran_flat(rng, 0, 2.*M_PI);
        double theta = gsl_ran_flat(rng, 0, M_PI);
        A.x[r] = cos(phi)*sin(theta);
        A.y[r] = sin(phi)*sin(theta);
        A.z[r] = cos(theta);
        A.ivec[r] = true;

        R1R2NOE((void *)&A, relax, NULL);
        for (int f=0; f < NF; f++) 
        {
            MHz = gsl_matrix_get(A.X, r, 3*f+0);
            R1  = gsl_vector_get(relax, 3*f+0);
            R2  = gsl_vector_get(relax, 3*f+1);
            NOE = gsl_vector_get(relax, 3*f+2);
            printf("    %6.2f %8.3f %8.3f %8.3f     XH= %7.4f %7.4f %7.4f\n", 
                MHz, R1, R2, NOE, A.x[r], A.y[r], A.z[r]);
            gsl_matrix_set(A.Y,   r, 3*f+0 , R1);
            gsl_matrix_set(A.Y,   r, 3*f+1 , R2);
            gsl_matrix_set(A.Y,   r, 3*f+2 , NOE);
        }
        
    }

    printf("    Truth   tc %7.4f Dr %7.4f phi %9.4f theta %9.4f\n",
        A.p[_tc_], A.p[_Dr_], A.p[_phi_]*(180./M_PI), A.p[_theta_]*(180./M_PI));

    for (int p=0; p < NP; p++)
    {
        A.is[p] = false;
    }
    A.is[_tc_] = true;
    A.is[_Dr_] = true;
    // A.is[_phi_] = true;
    // A.is[_theta_] = true;
    A.p[_tc_]  = 4.4;
    A.p[_Dr_]  = 1.0;
    // A.p[_phi_  ] = (M_PI/180.) * 10.0;
    // A.p[_theta_] = (M_PI/180.) * 10.0;
    
    printf("    Seed    tc %7.4f Dr %7.4f phi %9.4f theta %9.4f\n",
        A.p[_tc_], A.p[_Dr_], A.p[_phi_]*(180./M_PI), A.p[_theta_]*(180./M_PI));
    levmar (A, chisq, error);
    printf("    Final   tc %7.4f Dr %7.4f phi %9.4f theta %9.4f\n\n",
        A.p[_tc_], A.p[_Dr_], A.p[_phi_]*(180./M_PI), A.p[_theta_]*(180./M_PI));
    printf("            chisq   %g\n", chisq);
    printf("            tc      %7.4f +/- %7.5f\n", A.p[_tc_], gsl_vector_get(error, _tc_));
    printf("            Dr      %7.4f +/- %7.5f\n", A.p[_Dr_], gsl_vector_get(error, _Dr_));
    printf("            phi     %7.4f +/- %7.5f\n", (180./M_PI)*A.p[_phi_],     (180./M_PI)*gsl_vector_get(error, _phi_));
    printf("            theta   %7.4f +/- %7.5f\n", (180./M_PI)*A.p[_theta_],   (180./M_PI)*gsl_vector_get(error, _theta_));
    printf("\n");

    gsl_vector_free(relax);
    gsl_vector_free(error);
    freeing (A); // freeing memory
}


void test_isotropic ()
{
    int NF = 2;// number of magnetic fields
    int NR = 1;// number of residues
    int NK = 1;// number of tau points

    cout << "Testing isotropic fitting" << endl;
    ALLDATA A;

    Jw = &Jw_LSe;

    initialize (A, NM, NF, NR, NK, MC);
    
    D = _GLOBAL_ISOTROPIC_;
    setup_attribute (A);
    
    // internal dynamics
    A.p[_S2s_] = 1.00;
    A.p[_S2f_] = 0.90;
    A.p[_te_ ] = 0.05;
    A.p[_Rex_] = 0.55;
    // diffusion
    A.p[_tc_]  = 9.65;

    gsl_vector * relax = gsl_vector_calloc (3*A.NF);
    gsl_vector * error = gsl_vector_calloc (NP);

    double MHz, R1, R2, NOE, chisq;
    printf("    Correct tc = %7.4f\n", A.p[_tc_]);

    for (int r=0; r < NR; r++)
    {
        for (int f=0; f < NF; f++)
        {
            if (f == 0) 
            { 
                MHz = 500.38; 
            }
            else 
            {
                MHz = 600.13;
            }
            gsl_matrix_set(A.X,   r, 3*f+0, MHz);
            gsl_matrix_set(A.Y,   r, 3*f+0 , 1.0); // R1
            gsl_matrix_set(A.X,   r, 3*f+1, MHz);
            gsl_matrix_set(A.Y,   r, 3*f+1 , 1.0); // R2
            gsl_matrix_set(A.X,   r, 3*f+2, MHz);
            gsl_matrix_set(A.Y,   r, 3*f+2 , 1.0); // NOE
            gsl_matrix_set(A.SIG, r, 3*f+0 , 0.01); // dR1
            gsl_matrix_set(A.SIG, r, 3*f+1 , 0.01); // dR2
            gsl_matrix_set(A.SIG, r, 3*f+2 , 0.01); // dNOE
        }
        A.r = r;
        R1R2NOE((void *)&A, relax, NULL);
        for (int f=0; f < NF; f++) 
        {
            MHz = gsl_matrix_get(A.X, r, 3*f+0);
            R1  = gsl_vector_get(relax, 3*f+0);
            R2  = gsl_vector_get(relax, 3*f+1);
            NOE = gsl_vector_get(relax, 3*f+2);
            printf("    %6.2f %8.3f %8.3f %8.3f\n", MHz, R1, R2, NOE);
            gsl_matrix_set(A.Y,   r, 3*f+0 , R1);
            gsl_matrix_set(A.Y,   r, 3*f+1 , R2);
            gsl_matrix_set(A.Y,   r, 3*f+2 , NOE);
        }
    }

    for (int p=0; p < NP; p++)
    {
        A.is[p] = false;
    }
    A.is[_tc_] = true;
    A.p[_tc_]  = 1.4;
    
    printf("    Initial tc = %7.4f\n", A.p[_tc_]);
    levmar (A, chisq, error);
    printf("    Final   tc = %7.4f +/- %7.6f\n\n", A.p[_tc_], gsl_vector_get(error, _tc_));

    gsl_vector_free(relax);
    gsl_vector_free(error);
    freeing (A); // freeing memory
}