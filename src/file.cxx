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
  Data File Format

  lines begining with # will be regarded as comments
  remark na or NA for unavailable data

  Column 1 : residue number
  Column 2 : cluster
  Column 3 : 1H Larmoir frequency (MHz).
  Column 4 : R1 (1/s)
  Column 5 : R1 uncertainty (1/s)
  Column 6 : R2 (1/s)
  Column 7 : R2 uncertainty (1/s)
  Column 8 : NOE
  Column 9 : NOE uncertainty 
  Column 10: tau distribution file name (STRING/OPTIONAL)

  Stored Data Structure:
  x[0] = w0,  y[0] = R1   at w0, sig[0] = R1 error at w0
  x[1] = w0,  y[1] = R2   at w0, sig[1] = R2 error at w0
  x[2] = w0,  y[2] = NOE  at w0, sig[2] = NOE error at w0
  x[3] = w1,  y[3] = R1   at w1, sig[3] = R1 error at w1
  x[4] = w1,  y[4] = R2   at w1, sig[4] = R2 error at w1
  x[5] = w1,  y[5] = NOE  at w1, sig[5] = NOE error at w1
  ...         ...                 ...
  x[n] = 0.0 for empty data
*/

#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>
#include <sstream>
#include <iterator>

#include "emf.h"
#define	_COMMENT_ '#'
#define	_UNKNOWN_ '?'

using namespace std;

// trim from start
void ltrim(std::string &s) 
{
  s.erase(
    s.begin(),
    std::find_if(s.begin(),
    s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
}

// trim from end
void rtrim(std::string &s) 
{
  s.erase(
    std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(),
    s.end());
}

// trim from both ends
void trim(std::string &s) 
{
  rtrim(s);
  ltrim(s);
}

void parse (const char *line,vector <string> &tokens)
{
  tokens.clear();
  string s(line), item;
  trim(s);
  istringstream iss(s);
  copy(istream_iterator<string>(iss),
    istream_iterator<string>(),
    back_inserter<vector<string> >(tokens));
  vector<string>::iterator it = tokens.begin();
  while (it!= tokens.end())
  {
    if ((*it).substr(0,1) == "#")
    {
      tokens.erase(it,tokens.end());
      break;
    }
    it++;
  }
}


/* Check ALLDATAFILE */
void chk_data (const char *file,int &NF,int &NR,int &NK,bool idtk)
{
  int k,l,cols;
  int res,clu,NC;
  double MHz,R1,dR1,R2,dR2,NOE,dNOE;
  char line[_MAXSTR_],dtk[_MAXSTR_];
  char R1_[_MAXSTR_],dR1_[_MAXSTR_],R2_[_MAXSTR_],dR2_[_MAXSTR_];
  char NOE_[_MAXSTR_],dNOE_[_MAXSTR_];
  set <double>::iterator f,f0;
  set <int>::iterator c,c0;

  FILE *data_file,*dtk_file;
  NF=NR=NK=NC=0; // initialize
  l=0;

  if ((data_file=fopen(file,"r"))==NULL) {
    printf(_ERROR_ "cannot open data file, %s\n",file);
    exit(1);
    }

  while (fgets(line,sizeof(line),data_file)) {
    l++;// inc. line counter
    cols = sscanf(line,"%d %d %lf %s %s %s %s %s %s %s",
      &res,&clu,&MHz,R1_,dR1_,R2_,dR2_,NOE_,dNOE_,dtk);

    if(cols == -1 || line[0] == _COMMENT_) continue; /* skip comments */

    if (R1_[0] != _UNKNOWN_) R1 = atof(R1_); else R1 = 0.0;
    if (dR1_[0]!= _UNKNOWN_) dR1 = atof(dR1_); else dR1 = 1;
    if (R2_[0] != _UNKNOWN_) R2 = atof(R2_); else R2 = 0.0;
    if (dR2_[0]!= _UNKNOWN_) dR2 = atof(dR2_); else dR2 = 1;
    if (NOE_[0]!= _UNKNOWN_) NOE = atof(NOE_); else NOE = 0.0;
    if (dNOE_[0]!=_UNKNOWN_) dNOE = atof(dNOE_); else dNOE = 1;

    /* check data file format */
    if (!(cols == 9 || cols == 10)) {
      printf(_ERROR_ "invalid format at %s, line %d\n",file,l);
      printf(_ERROR_ "%s\n",line);
      exit(1);
      }
    /* check magnetic field */
    if (MHz < 0.0  || MHz > 10000.0) {
      printf(_ERROR_ "invalid larmoir frequency at %s, line %d\n",file,l);
      printf(_ERROR_ "%s\n",line);
      exit(1);
      }
    /* check R1,R2,hnNOE data */
    if (R1 < 0.0 || R2 < 0.0 || NOE > 1.0) {
      printf(_ERROR_ "invalid data at %s, line %d\n",file,l);
      printf(_ERROR_ "%s\n",line);
      exit(1);
      }
    /* check R1,R2,hnNOE uncertainty */
    if (dR1 <= 0.0 || dR2 <= 0.0 || dNOE <= 0.0) {
      printf(_ERROR_ "invalid uncertainty at %s, line %d\n",file,l);
      printf(_ERROR_ "%s\n",line);
      exit(1);
      }
    /* check tau distribution file */
    if(cols==9 && idtk) {
      printf(_ERROR_ "requires tau distribution file\n");
      exit(1);
      }
    /* check distribution file */
    if(cols==10 && idtk) {
      if((dtk_file = fopen(dtk,"r"))==NULL) {
	printf(_ERROR_ "cannot open tau distribution file, %s\n",dtk);
	exit(1);
	}
      else {
	k=0;
	while (fgets(line,sizeof(line),dtk_file)) {
	  if(strlen(line) && line[0] != '#') k++; 
	  }
	if (k > NK) NK=k; 
	}//else
      fclose(dtk_file);
      }// idtk

    cluster_set.insert(clu);// count number of clusters 
    residue_set.insert(res);// count number of residues
    field_set.insert(MHz);// count number of fields
    min_err1[MHz] = 0.0;
    min_err2[MHz] = 0.0;
    min_err3[MHz] = 0.0;
    scl_err1[MHz] = 1.0;
    scl_err2[MHz] = 1.0;
    scl_err3[MHz] = 1.0;
    } // while

  fclose (data_file);

  NF = field_set.size();
  NR = residue_set.size();
  NC = cluster_set.size();

  show_info("data");printf("%s\n",file);
  show_info("  number of lines");printf("%d\n",l);
  show_info("  number of clusters");printf("%d (",NC);
  c=c0=cluster_set.begin(); 
  while (c!=cluster_set.end()) {
    if(c!=c0) printf(", ");
    printf("%d",(*c));
    c++;
    }
  printf(")\n");
  show_info("  number of residues");printf("%d\n",NR);
  show_info("  number of magnetic fields");printf("%d (",NF);
  f=f0=field_set.begin(); 
  while (f!=field_set.end()) {
    if(f!=f0) printf(", ");
    printf("%g",(*f));
    f++;
    }
  printf(")\n");
  if (idtk) {
    show_info("number of distribution points");
    printf("%d\n",NK);
    }
}

/* Read ALLDATAFILE  */
void read_data (const char *name,ALLDATA &A,bool idtk)
{
  int res,clu; 
  double MHz,R1,dR1,R2,dR2,NOE,dNOE,tau,wt;
  map <double,double>::iterator mp;
  set <int> residues;
  char dtk[_MAXSTR_],line[_MAXSTR_];
  char R1_[_MAXSTR_],dR1_[_MAXSTR_],R2_[_MAXSTR_],dR2_[_MAXSTR_];
  char NOE_[_MAXSTR_],dNOE_[_MAXSTR_];
  int r,f=0,k,cols;

  r=-1;//residue index

  FILE *data_file,*dtk_file;

  data_file = fopen(name,"r"); 
  while (fgets(line,sizeof(line),data_file)) {
    cols = sscanf(line,"%d %d %lf %s %s %s %s %s %s %s",
	&res,&clu,&MHz,R1_,dR1_,R2_,dR2_,NOE_,dNOE_,dtk);

    if(cols == -1 || line[0] == '#') continue;

    /* store data */
    if (residues.find(res) == residues.end()) {// new residue
      residues.insert(res);
      r++; // inc. residue index. r=[0...NR-1]
      f=0; // reset field index
      A.num[r]=res;
      A.clst[r]=clu;
      R1 = R2 = NOE = 0.0;
      dR1 = dR2 = dNOE = 1e-3;
      if (R1_[0] != _UNKNOWN_) R1 = atof(R1_);
      if (dR1_[0] != _UNKNOWN_) dR1 = atof(dR1_); 
      if (R1_[0] == _UNKNOWN_)
	gsl_matrix_set(A.X,r,0,0);
      else
	gsl_matrix_set(A.X,r,0,MHz);
	gsl_matrix_set(A.Y,r,0,R1);
	mp = scl_err1.find(MHz);
	gsl_matrix_set(A.SIG,r,0,dR1*(mp->second));
	mp = min_err1.find(MHz);
	gsl_matrix_set(A.SIG,r,0,max(dR1,(mp->second)*R1));

      if (R2_[0] != _UNKNOWN_) R2 = atof(R2_);
      if (dR2_[0] != _UNKNOWN_) dR2 = atof(dR2_);
      if (R2_[0] == _UNKNOWN_)
	gsl_matrix_set(A.X,r,1,0);
      else
	gsl_matrix_set(A.X,r,1,MHz);
	gsl_matrix_set(A.Y,r,1,R2);
	mp = scl_err2.find(MHz);
	gsl_matrix_set(A.SIG,r,1,dR2*(mp->second));
	mp = min_err2.find(MHz);
	gsl_matrix_set(A.SIG,r,1,max(dR2,(mp->second)*R2));

      if (NOE_[0] != _UNKNOWN_) NOE = atof(NOE_);
      if (dNOE_[0] != _UNKNOWN_) dNOE = atof(dNOE_);
      if (NOE_[0] == _UNKNOWN_)
	gsl_matrix_set(A.X,r,2,0);
      else
	gsl_matrix_set(A.X,r,2,MHz);
	gsl_matrix_set(A.Y,r,2,NOE);
	mp = scl_err3.find(MHz);
	gsl_matrix_set(A.SIG,r,2,dNOE*(mp->second));
	mp = min_err3.find(MHz);
	gsl_matrix_set(A.SIG,r,2,max(dNOE,(mp->second)));

      if(idtk) { /* Read tau distribution file */
	k = 0;
	dtk_file = fopen(dtk,"r");
	while (fgets(line,sizeof(line),dtk_file)) {
	  if(strlen(line)==0 || line[0] == '#') continue;
	  sscanf(line,"%lf %lf",&tau,&wt);
	  gsl_matrix_set(A.TAU,r,k,tau);
	  gsl_matrix_set(A.WT,r,k,wt);
	  k++;
	  }
	fclose (dtk_file);
	} // idtk
      }// for new residue

    else
      {// for existing residue
      f++; // inc. field index
      R1 = R2 = NOE = 0.0;
      dR1 = dR2 = dNOE = 1e-3;
      if (R1_[0] != _UNKNOWN_) R1 = atof(R1_);
      if (dR1_[0] != _UNKNOWN_) dR1 = atof(dR1_); 
      if (R1_[0] == _UNKNOWN_)
	gsl_matrix_set(A.X,r,3*f+0,0);
      else
	gsl_matrix_set(A.X,r,3*f+0,MHz);
	gsl_matrix_set(A.Y,r,3*f+0,R1);
	mp = scl_err1.find(MHz);
	gsl_matrix_set(A.SIG,r,3*f+0,dR1*(mp->second));
	mp = min_err1.find(MHz);
	gsl_matrix_set(A.SIG,r,3*f+0,max(dR1,(mp->second)*R1));

      if (R2_[0] != _UNKNOWN_) R2 = atof(R2_);
      if (dR2_[0] != _UNKNOWN_) dR2 = atof(dR2_);
      if (R2_[0] == _UNKNOWN_)
	gsl_matrix_set(A.X,r,3*f+1,0);
      else
	gsl_matrix_set(A.X,r,3*f+1,MHz);
	gsl_matrix_set(A.Y,r,3*f+1,R2);
	mp = scl_err2.find(MHz);
	gsl_matrix_set(A.SIG,r,3*f+1,dR2*(mp->second));
	mp = min_err2.find(MHz);
	gsl_matrix_set(A.SIG,r,3*f+1,max(dR2,(mp->second)*R2));

      if (NOE_[0] != _UNKNOWN_) NOE = atof(NOE_);
      if (dNOE_[0] != _UNKNOWN_) dNOE = atof(dNOE_);
      if (NOE_[0] == _UNKNOWN_)
	gsl_matrix_set(A.X,r,3*f+2,0);
      else
	gsl_matrix_set(A.X,r,3*f+2,MHz);
	gsl_matrix_set(A.Y,r,3*f+2,NOE);
	mp = scl_err3.find(MHz);
	gsl_matrix_set(A.SIG,r,3*f+2,dNOE*(mp->second));
	mp = min_err3.find(MHz);
	gsl_matrix_set(A.SIG,r,3*f+2,max(dNOE,(mp->second)));

      }// existing residue
    }// while
  fclose(data_file);
}

void read_config (const char *name, ALLDATA &A)
{
  map <double,double>::iterator it;
  double MHz,err1,err2,err3;

  char line[_MAXSTR_];
  vector <string> c;//parsed columns
  int i,j,l=0;
  FILE *config_file;

  if ((config_file = fopen(name,"r")) == NULL) {
    printf(_ERROR_ "cannot open config file\n");
    exit(1);
    }

  while (fgets(line,sizeof(line),config_file)) {

    l++;

    if(strlen(line)==0 || line[0] == '#') continue; /* skip comments */
    parse(line,c);//parse line into c (vector of string)
    if (c.size() == 0) continue;

    /* const [gamma_h ###] [gamma_x ###] [csa_x ###] [r_xh ###] */
    if(c[0]=="const") {
      for (i=1;i<c.size();i+=2) {
	if(c[i]=="gamma_h") sscanf(c[i+1].c_str(),"%lf",&gamma_h);
	else if(c[i]=="gamma_x") sscanf(c[i+1].c_str(),"%lf",&gamma_x);
	else if(c[i]=="csa_x") sscanf(c[i+1].c_str(),"%lf",&csa_x);
	else if(c[i]=="r_xh") sscanf(c[i+1].c_str(),"%lf",&r_xh);
	else {
	  printf(_ERROR_ "syntax error at config file, line %d\n",l);
	  exit(1);
	  }
	}
      } // const

    /* coor [file name] [-h atom name] [-x atom name] */
    if(c[0]=="coor") {
      pdbfile = c[1];
      atomh = "";
      atomx = "";
      for (i=2;i<c.size();i+=2) {
	if (c[i]=="-h") atomh=c[i+1];
	else if (c[i]=="-x") atomx=c[i+1];
	else {
	  printf(_ERROR_ "syntax error at config file, line %d\n",l);
	  exit(1);
	  }
	}
      }// coor

    /* mle [-cluster #] [-mc #] */
    if(c[0]=="mle") {
      printf(_INFO_);
      printf("maximum likelihood estimation ( ");
      for (i=1;i<c.size();i+=2) {
	if (c[i]=="-cluster") {
	  sscanf(c[i+1].c_str(),"%d",&MLE_cluster);
	  printf("cluster %d ",MLE_cluster);
	  }
	else if (c[i]=="-mc") {
	  sscanf(c[i+1].c_str(),"%d",&MLE_MC);
	  printf(" & MC %d ",MLE_MC);
	  }
	else {
	  printf(_ERROR_ "syntax error at config file, line %d\n",l);
	  exit(1);
	  }
	}
      printf(")\n");
      }// mle

    /* estimate [-cluster # ] [-noe #] [-r2stdev #] */
    if(c[0]=="estimate") {
      estimate = true;
      for (i=1;i<c.size();i+=2) {
	if (c[i] == "-cluster")
	  sscanf(c[i+1].c_str(),"%d",&EstCluster);
	else if (c[i] == "-noe")
	  sscanf(c[i+1].c_str(),"%lf",&EstNOECut);
	else if (c[i] == "-r2stdev")
	  sscanf(c[i+1].c_str(),"%lf",&EstSTDCut);
	else {
	  printf(_ERROR_ "syntax error at config file, line %d\n",l);
	  exit(1);
	  }
	}
      } // estimate

    /* axial_grid_search */
    if(c[0]=="axial_grid_search") {
      axial_grid_search = true;
      }

    /* optimize [-cluster #] [-s2 #] [-maxiter #] [-method #]*/
    if(c[0]=="optimize") {
      optimize = true;
      for (i=1;i<c.size();i+=2) {
	if(c[i]=="-cluster")
	  sscanf(c[i+1].c_str(),"%d",&OptCluster);
	else if(c[i]=="-s2")
	  sscanf(c[i+1].c_str(),"%lf",&OptS2);
	else if(c[i]=="-maxiter")
	  sscanf (c[i+1].c_str(),"%d",&OptMaxIter);
	else if(c[i]=="-method")
	  sscanf (c[i+1].c_str(),"%d",&OptMethod);
	else {
	  printf(_ERROR_ "syntax error at config file, line %d\n",l);
	  exit(1);
	  }
	}
      }// optimize

    /* fitmodel [-select AIC|BIC] [-critx2 #] [-mc #] [-mctrim #] */
    if(c[0]=="fitmodel") {
      for (i=1;i<c.size();i+=2) {
	if (c[i]=="-select") {
	  if (c[i+1]=="BIC") criterion = _BIC_;
	  if (c[i+1]=="AIC") criterion = _AIC_;
	  }
	else if (c[i]=="-mc")
	  sscanf(c[i+1].c_str(),"%d",&MC);
	else if (c[i]=="-mctrim")
	  sscanf(c[i+1].c_str(),"%d",&MC_trim);
	else if (c[i]=="-critx2")
	  sscanf(c[i+1].c_str(),"%lf",&critx2);
	else {
	  printf(_ERROR_ "syntax error at config file, line %d\n",l);
	  exit(1);
	  }
	}
      }// fitmodel

    /* error [-field #] [-min # # #] [-scale # # #] */
    if(c[0]=="error") {
      for (i=1;i<c.size();i++) {
	if (c[i]=="-field") {
	  if (c[i+1] == "*") MHz = -1.0;
	  else sscanf(c[i+1].c_str(),"%lf",&MHz);
	  i++;
	  }
	else if (c[i]=="-min") {
	  sscanf(c[i+1].c_str(),"%lf",&err1);
	  sscanf(c[i+2].c_str(),"%lf",&err2);
	  sscanf(c[i+3].c_str(),"%lf",&err3);
	  if (MHz == -1) {
	    for (it=min_err1.begin();it != min_err1.end();it++)
	      it->second = err1;
	    for (it=min_err2.begin();it != min_err2.end();it++)
	      it->second = err2;
	    for (it=min_err3.begin();it != min_err3.end();it++)
	      it->second = err3;
	    }
	  else {
	    min_err1[MHz]=err1;
	    min_err2[MHz]=err2;
	    min_err3[MHz]=err3;
	    }
	  i+=3;
	  }
	else if (c[i]=="-scale") {
	  sscanf(c[i+1].c_str(),"%lf",&err1);
	  sscanf(c[i+2].c_str(),"%lf",&err2);
	  sscanf(c[i+3].c_str(),"%lf",&err3);
	  if (MHz == -1) {
	    for (it=scl_err1.begin();it != scl_err1.end();it++)
	      it->second = err1;
	    for (it=scl_err2.begin();it != scl_err2.end();it++)
	      it->second = err2;
	    for (it=scl_err3.begin();it != scl_err3.end();it++)
	      it->second = err3;
	    }
	  else {
	    scl_err1[MHz]=err1;
	    scl_err2[MHz]=err2;
	    scl_err3[MHz]=err3;
	    }
	  i+=3;
	  }
	else {
	  printf(_ERROR_ "syntax error at config file, line %d\n",l);
	  exit(1);
	  }
	}
      }// error

    /* par [name][-lb #][-ub #][-gr #][-cv #][-lk #][-uk #][-fx #][-in #] */
    if(c[0]=="par") {
      for (i=0;i<NP;i++) 
	if (strcmp(c[1].c_str(),A.pid[i])==0) break;
      for (j=2;j<c.size();j+=2) {
	if (c[j] == "-lb") sscanf(c[j+1].c_str(),"%lf", &A.lb[i]);
	else if (c[j] == "-ub") sscanf(c[j+1].c_str(),"%lf", &A.ub[i]);
	else if (c[j] == "-gr") {
	  sscanf(c[j+1].c_str(),"%d", &A.grds[i]);
	  A.step[i]=(A.ub[i]-A.lb[i])/A.grds[i];
	  }
	else if (c[j] == "-cv") sscanf(c[j+1].c_str(),"%lf",&A.conv[i]);
	else if (c[j] == "-lk") sscanf(c[j+1].c_str(),"%lf",&A.lk[i]);
	else if (c[j] == "-uk") sscanf(c[j+1].c_str(),"%lf",&A.uk[i]);
	else if (c[j] == "-in") sscanf(c[j+1].c_str(),"%lf",&A.p[i]);
	else if (c[j] == "-fx") {
	  sscanf(c[j+1].c_str(),"%lf",&A.p[i]);
	  A.attr[i] |= P_FIXED;
	  }
	else {
	  printf(_ERROR_ "syntax error at config file, line %d\n",l);
	  printf(_ERROR_ "%s\n",line);
	  exit(1);
	  }
	} // j
      } // par

    }// while
  /* convert degree to radian */
  A.lb[_phi_]     *= M_PI/180.0;
  A.ub[_phi_]     *= M_PI/180.0;
  A.step[_phi_]   *= M_PI/180.0;
  A.conv[_phi_]   *= M_PI/180.0;
  A.p[_phi_]      *= M_PI/180.0;
  A.lb[_theta_]   *= M_PI/180.0;
  A.ub[_theta_]   *= M_PI/180.0;
  A.step[_theta_] *= M_PI/180.0;
  A.conv[_theta_] *= M_PI/180.0;
  A.p[_theta_]    *= M_PI/180.0;
  A.lb[_psi_]     *= M_PI/180.0;
  A.ub[_psi_]     *= M_PI/180.0;
  A.step[_psi_]   *= M_PI/180.0;
  A.conv[_psi_]   *= M_PI/180.0;
  A.p[_psi_]      *= M_PI/180.0;
  fclose(config_file);
}

