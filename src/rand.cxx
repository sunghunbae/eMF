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
  Random Search

  definition of grid points
  lb ... ub, step
  1st grid point  : lb + step
  last grid point : ub - step
  total number of points: (grds - 1)
  lb and ub are excluded from grid evaluation
*/

#include "emf.h"

using namespace std;

void grid_search (ALLDATA &A, int func, double &chisq, bool output)
{
  double unit;

  for (int i=0;i<MA;i++)
    if (A.is[i])
  for (int i=0;i<NP;i++) grid_min_par[i]=0.0;
  recursive_s (A, func, 0);
  for (int i=0;i<NP;i++) A.p[i] = grid_min_par[i];
  chisq = grid_min_x2;
}
