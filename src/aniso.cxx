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

/*
	Euler angles (phi, theta, psi), "x-convention"

	A = BCD

	D = rotation by an angle phi about the z-axis
	C = rotation by an angle theta [0...PI] about the x-axis
	B = rotation by an angle psi about the z-axis (again)
	
	D
	| cos(phi) sin(phi) 0 |
	|-sin(phi) cos(phi) 0 |
	| 0        0        1 |

	C
	| 1        0        0      |
	| 0  cos(theta) sin(theta) |
	| 0 -sin(theta) cos(theta) |

	B
	| cos(psi) sin(psi) 0 |
	|-sin(psi) cos(psi) 0 |
	| 0        0        1 |

	A, overall rotation matrix

	| a11 a12 a13 |
	| a21 a22 a23 |
	| a31 a32 a33 |

	cartesian axis: 
	X-axis (1,0,0)
	Y-axis (0,1,0)
	Z-axis (0,0,1)

	principal axis: 
	X-axis (px[0],px[1],px[2])
	Y-axis (py[0],py[1],py[2])
	Z-axis (pz[0],pz[1],pz[2])

	rotate cartesian axis frame to principal axis frame

	px[0] =  cos(psi)*cos(theta)-cos(theta)*sin(phi)*sin(psi); // a11
	px[1] =  cos(psi)*sin(theta)+cos(theta)*cos(phi)*sin(psi); // a12
	px[2] =  sin(psi)*sin(theta); // a13
	py[0] = -sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi); // a21
	py[1] = -sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi); // a22
	py[2] = -cos(psi)*sin(theta); // a23
	pz[0] =  sin(theta)*sin(phi); // a31
	pz[1] = -sin(theta)*cos(phi); // a32
	pz[2] =  cos(theta); // a33

*/

/* 
 * fill-in vector-dependent weight for each correlatin time 
 * should be called whenever Dxx,Dyy,Dzz,phi,theta,psi or tm,Dr are changed
 *
 */

void anisotropic (ALLDATA &A)
{
/* fully anisotropic diffusion tensor */
if (D == _GLOBAL_ANISOTROPIC_) {
	double Dxx = A.p[_Dxx_];
	double Dyy = A.p[_Dyy_];
	double Dzz = A.p[_Dzz_];
	double Diso = (Dxx+Dyy+Dzz)/3.0;
	double L2 = (Dxx*Dyy+Dxx*Dzz+Dyy*Dzz)/3.0;
	double f1 = sqrt(Diso*Diso - L2);
	double dx,dy,dz;
	f1==0 ? dx=0 : dx = (Dxx-Diso)/f1;
	f1==0 ? dy=0 : dy = (Dyy-Diso)/f1;
	f1==0 ? dz=0 : dz = (Dzz-Diso)/f1;
	double phi = deg2rad(A.p[_phi_]);
	double theta = deg2rad(A.p[_theta_]);
	double psi = deg2rad(A.p[_psi_]);
	double x,y,z;
	double xi,yi,zi; // direction cosine
	double f2,f3,f4,f5,f6,f7,f8,f9;
	double px[3],py[3],pz[3]; 

	A.t[0] = 1.0/(4*Dxx + Dyy + Dzz);
	A.t[1] = 1.0/(Dxx + 4*Dyy + Dzz);
	A.t[2] = 1.0/(Dxx + Dyy + 4*Dzz);
	A.t[3] = 1.0/(6*Diso + 6*f1);
	A.t[4] = 1.0/(6*Diso - 6*f1);
	}
/* axially symmetric diffusion tensor (Dxx = Dyy) */
if (D == _GLOBAL_AXIAL_) {
	double phi = deg2rad (A.p[_phi_]);
	double theta = deg2rad (A.p[_theta_]);
	double x,y,z,cosa,cos2a,sin2a;

	A.p[_Dxx_] = A.p[_Dyy_] = 0.5/A.p[_tc_]/(2.0+A.p[_Dr_]);
	A.p[_Dzz_] = A.p[_Dxx_]*A.p[_Dr_];

	A.t[0] = 1.0/(6.0*A.p[_Dxx_]);
	A.t[1] = 1.0/(5.0*A.p[_Dxx_] + 1.0*A.p[_Dzz_]);
	A.t[2] = 1.0/(2.0*A.p[_Dxx_] + 4.0*A.p[_Dzz_]);

	for (int r=0;r<A.NR;r++) {
		if (A.flag[r]) {
			/* projection onto the pricipal z-axis */
			cosa = A.x[r]*sin(theta)*sin(phi)
				-A.y[r]*sin(theta)*cos(phi)
				+A.z[r]*cos(theta);
			cos2a = cosa*cosa;
			sin2a = 1.0-cos2a;
			gsl_matrix_set(A.WT,r,0,SQR(1.5*cos2a-0.5));
			gsl_matrix_set(A.WT,r,1,3.0*sin2a*cos2a);
			gsl_matrix_set(A.WT,r,2,0.75*sin2a*sin2a);
			A.alpha[r] = rad2deg(acos(cosa));
			} // residue loop
		}
	}// axially symmetric


/* fully anisotropic diffusion tensor */
if (D == _GLOBAL_ANISOTROPIC_) {
	double Dxx = A.p[_Dxx_];
	double Dyy = A.p[_Dyy_];
	double Dzz = A.p[_Dzz_];
	double Diso = (Dxx+Dyy+Dzz)/3.0;
	double L2 = (Dxx*Dyy+Dxx*Dzz+Dyy*Dzz)/3.0;
	double f1 = sqrt(Diso*Diso - L2);
	double dx,dy,dz;
	f1==0 ? dx=0 : dx = (Dxx-Diso)/f1;
	f1==0 ? dy=0 : dy = (Dyy-Diso)/f1;
	f1==0 ? dz=0 : dz = (Dzz-Diso)/f1;
	double phi = deg2rad(A.p[_phi_]);
	double theta = deg2rad(A.p[_theta_]);
	double psi = deg2rad(A.p[_psi_]);
	double x,y,z;
	double xi,yi,zi; // direction cosine
	double f2,f3,f4,f5,f6,f7,f8,f9;
	double px[3],py[3],pz[3]; 

	A.t[0] = 1.0/(4*Dxx + Dyy + Dzz);
	A.t[1] = 1.0/(Dxx + 4*Dyy + Dzz);
	A.t[2] = 1.0/(Dxx + Dyy + 4*Dzz);
	A.t[3] = 1.0/(6*Diso + 6*f1);
	A.t[4] = 1.0/(6*Diso - 6*f1);

	px[0] =  cos(psi)*cos(theta)-cos(theta)*sin(phi)*sin(psi); // a11
	px[1] =  cos(psi)*sin(theta)+cos(theta)*cos(phi)*sin(psi); // a12
	px[2] =  sin(psi)*sin(theta); // a13
	py[0] = -sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi); // a21
	py[1] = -sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi); // a22
	py[2] = -cos(psi)*sin(theta); // a23
	pz[0] =  sin(theta)*sin(phi); // a31
	pz[1] = -sin(theta)*cos(phi); // a32
	pz[2] =  cos(theta); // a33

	for (int r=0;r<A.NR;r++) {
		/* XH-bond vector (unit vector) */
		x = A.x[r];
		y = A.y[r];
		z = A.z[r];
		/* direction cosine relative to principal axis frame */
		xi = (x*px[0]+y*px[1]+z*px[2]);
		yi = (x*py[0]+y*py[1]+z*py[2]);
		zi = (x*pz[0]+y*pz[1]+z*pz[2]);
		f2 = xi*xi;
		f3 = yi*yi;
		f4 = zi*zi;
		f5 = (1.0/4.0)*(3*(f2*f2 + f3*f3 + f4*f4) - 1);
		f6 = 3*f2*f2 + 6*f3*f4 - 1;
		f7 = 3*f3*f3 + 6*f2*f4 - 1;
		f8 = 3*f4*f4 + 6*f2*f3 - 1;
		f9 = (1.0/12.0)*(dx*f6 + dy*f7 + dz*f8);

		gsl_matrix_set(A.WT,r,0,3*f3*f4);
		gsl_matrix_set(A.WT,r,1,3*f2*f4);
		gsl_matrix_set(A.WT,r,2,3*f2*f3);
		gsl_matrix_set(A.WT,r,3,f5-f9);
		gsl_matrix_set(A.WT,r,4,f5+f9);
		} // NR
	}
}


int read_pdb (const char* filename,ALLDATA &A)
{
  int c,i,NR=A.NR,res;
  char line[_MAXSTR_];
  string s,name,xyz,resSeq,resName;
  double x,y,z,l;
  double Xx[NR],Xy[NR],Xz[NR];
  double Hx[NR],Hy[NR],Hz[NR];
  bool isX[NR],isH[NR];

  ifstream pdb_file(filename);
  if (pdb_file.fail()) {
    stringstream error;
    error << "cannot open pdb file, "<<filename;
    terminate(error.str());
    }

  while (! pdb_file.eof()) {
    pdb_file.getline(line, sizeof(line));
    if (strlen(line) > 0) {
      s = line;
      if(s.find("ATOM") != string::npos) {
	// PDB format version 2.3
	name = s.substr(12,4);
	resName = s.substr(17,3);
	resSeq = s.substr(22,4);
	xyz = s.substr(30,24);
        trim(name);
	if (name == atomx) {
	  sscanf(resSeq.c_str(),"%d",&res);
	  for (i=0;i<NR;i++) 
	    if(A.num[i] == res) {
	      sscanf(xyz.c_str(),"%lf %lf %lf",&Xx[i],&Xy[i],&Xz[i]);
	      isX[i] = true;
	      }
	  } // X 
	if (name == atomh) {
	  sscanf(resSeq.c_str(),"%d",&res);
	  for (i=0;i<NR;i++) 
	    if(A.num[i] == res) {
	      sscanf(xyz.c_str(),"%lf %lf %lf",&Hx[i],&Hy[i],&Hz[i]);
	      isH[i] = true;
	      }
	  } // H
	} // ATOM
      } // line
    } // EOF

  // calculate XH-bond vector (unit-vector)
  for(c=0,i=0;i<NR;i++)
    if(isX[i] && isH[i]) {
      c++;
      A.x[i] = Xx[i] - Hx[i];
      A.y[i] = Xy[i] - Hy[i];
      A.z[i] = Xz[i] - Hz[i];
      l = sqrt (A.x[i]*A.x[i] + A.y[i]*A.y[i] + A.z[i]*A.z[i]);
      A.x[i] /= l;
      A.y[i] /= l;
      A.z[i] /= l;
      A.ivec[i] = true;
      }
  return c;
}
