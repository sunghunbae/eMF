What is eMF ?
=============

eMF stands for easy Model Free. Why do I have to write another model free program ?
It was originally written to implement distributed correlation time approaches such as 
Cole-Cole and Lorentzian spectral density functions in the model free analysis of 
unfolded or disordered proteins. During my time with Drs. Peter E. Wright and Jane H. Dyson
at The Scripps Research Institute in La Jolla, CA, I studied the ps-ns dynamics of 
wild type and mutant prion proteins. Prion protein has a very long stretch of disordered 
amino acid sequence at the N-terminus which significantly affects observed NMR relaxation 
rates of the whole protein including the C-terminal structured domain. So, prion protein 
should be described by an ensemble of vastly different overall structures due to 
the disordered N-terminal domain.
eMF was born from an attempt to address this problem by applying the 'distributed' correlation time.
Although it was found later that mean rotational diffusion tensor determined from conventional
methods should be good enough for the prion case, eMF remained as my choice of model free analysis tool.

### Features

- separate config and data files
- reduced spectral mapping
- estimation and optimization of diffusion tensor (isotropic, axially symmetric, asymmetric, local, distributed)
- fit extended Lipari-Szabo models with selection options of BIC and AIC
- xmgrace plot of select model spectral density function

How to install ?
================

You need the following two packages installed in your system:

- CMake
- GNU Scientific Library (GSL)

```
$ mkdir build
$ cd build
$ make 
$ make install
```

By default, an executable binary file, eMF, will be installed in ```~/bin``` directory.
You may change the destination directory defined the ```CMakeLists.txt``` file :

```
SET(CMAKE_INSTALL_PREFIX ~/bin )
```

Please make sure that eMF is accessible in your working directory 
by adding ```~/bin``` or your destination directory in the 'path'.

How to prepare input files ?
============================

You need to prepare config and data files.
These files are plain/text file and please use the following examples files
as templates.

### Config file

```
#
# eMF configuration
#

# CONSTANTS
# const [gamma_h ###] [gamma_x ###] [csa_x ###] [r_xh ###] 
const gamma_h   26.7519 # gamma for H
const gamma_x   -2.71   # gamma for X
const r_xh       1.02   # bond length (A)
const csa_x   -172.0    # chemical shift anisotropy (ppm)

# PARAMETERS
# par [name][-lb #][-ub #][-gr #][-cv #][-lk #][-uk #][-fx #][-in #]
#      name= S2f S2s te Rex tc Dr a scf Dxx Dyy Dzz phi theta psi
# -lb : lower bound
# -ub : upper bound
# -gr : grid points
# -cv : convergence limit
# -lk : lower bound harmonic chi-square penalty
#       chisq += lk * (v - lb)^2
# -uk : upper bound harmonic chi-square penalty
#       chisq += uk * (v - ub)^2
# -fx : fix value
# -in : initial value

# internal dynamics

par S2f   -lb 0.0 -ub   1.0 -gr 20
par S2s   -lb 0.0 -ub   1.0 -gr 20
par te    -lb 0.0 -ub   3.0 -gr 30
par Rex   -lb 0.0 -ub  10.0 -gr 10

# diffusion = axially symmetric

par tc    -lb 7.5 -ub   9.5 -gr 20 -cv 1e-3   #-fx  8.630
par Dr    -lb 1.0 -ub   1.3 -gr 10 -cv 0.01   #-fx  1.195
par phi   -lb 0   -ub 360   -gr 24 -cv 0.1     -fx -2.385
par theta -lb 0   -ub 180   -gr 12 -cv 0.1     -fx  0.000

# PDB COORDINATES (NECESSARY FOR NON-ISOTROPIC MODELS)
# coor [file name] [-h atom name] [-x atom name]
coor 1rx7.rot.pdb -h H -x N

# ADJUST EXPERIMENTAL ERRORS
# error [-field #] [-min # # #] [-scale # # #]
# -field #     : specifiy field or * for all data
# -min # # #   : set minimum proportioanl errors for R1,R2,NOE 
# -scale # # # : scale errors for R1,R2,NOE           

#error -field * -min 0.05 0.05 0.05
# i.e if given errors are less than 5% of measured values,
#     set errors as 5% of measured values for all fields

# ESTIMATION OF DIFFUSION PARAMETER(S) by fitting R2/R1
# estimate [-cluster # ]

estimate -cluster 1

# OPTIMIZATION OF DIFFUSION PARAMETER(S) by fitting R1,R2,NOE
# optimize [-cluster #] [-s2 #] [-maxiter #]
#
# -cluster # : use cluster #
# -s2 #      : use residues with order parameter larger than s2
# -maxiter # : maximum iteration

optimize -cluster 1 -s2 0.7
# i.e. use cluster 1 AND S2 > 0.7 for optimization

# FITTING AND SELECTING MODELS
#
# fitmodel [-select AIC|BIC] [-mc #] [-mctrim #]
#
# -select AIC : use Akke information criterion
# -select BIC : use Bayesian information criterion

fitmodel -select BIC
```

### Data file

- column 1 : residue number
- column 2 : cluster number, used for inclusion or exclusion of data
- column 3 : B0, field strength in MHz
- column 4-5 : R1 relaxation rate and error (1/s)
- column 6-7 : R2 relaxation rate and error (1/s)
- column 8-9 : {1H}-X NOE and error
- column 10  : (only for distributed) tc distribution file name

```
  7 1 500.38  1.6265  0.0691 10.7276  0.4488  0.7760  0.0400
  7 1 600.13  1.2629  0.0436 11.1164  0.7537  0.8070  0.0880
  8 1 500.38  1.6935  0.0442 10.8300  0.3911  0.7730  0.0240
  8 1 600.13  1.3243  0.0559 12.1164  0.3416  0.7990  0.0220
  9 0 500.38  1.6324  0.0658 15.3615  0.6703  0.7930  0.0260
  9 0 600.13  1.2494  0.0476 16.6599  0.5612  0.7710  0.0660
 10 1 500.38  1.7236  0.0398 10.5304  0.5007  0.7900  0.0400
 10 1 600.13  1.3388  0.0534 10.7676  0.4242  0.7830  0.0420
 11 1 500.38  1.5871  0.0856 10.7627  0.4155  0.8040  0.0700
 11 1 600.13  1.2146  0.0245 10.0442  0.2884  0.7720  0.1080
```

### Estimation and optimization of diffusion tensor

Generally, diffusion tensor should be determined before fitting the models.

+===================+=====================+============+
| diffusion tensor  + parameter(s)        + eMF symbol +
+===================+=====================+============+
+ isotropic         + tc                  + I          +
+ axially symmetric + tc, Dr              + A          +
+ anisotropic       + tc, Dr, phi, theta  + N          +
+ distributed       + (datafile)          + D          +
+===================+=====================+============+
