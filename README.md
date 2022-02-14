What is eMF?
=============

eMF stands for easy Model Free. Why did I have to write another model free program?
It was originally written in order to implement distributed correlation time approaches 
such as Cole-Cole and Lorentzian spectral density functions for the model free analysis of 
unfolded or disordered proteins. 

During my time with Drs. Peter E. Wright and Jane H. Dyson at The Scripps Research Institute 
in La Jolla, CA, I studied the ps-ns dynamics of wild type and mutant prion proteins. 
Prion protein has a very long stretch of disordered amino acid sequence at the N-terminus 
which significantly alters observed NMR relaxation rates of the whole protein 
including the C-terminal structured domain. Apparently, it should be described by 
an ensemble of vastly different overall structures due to the disordered N-terminal domain.

eMF was born out of an attempt to resolve this problem by applying the 'distributed' correlation time.
It was written to be easy and straightforward to use and to put together recent model free protocols.

### Updated

- Compatible with GSL 2.x and GSL 1.x
- Using CMake default FindGSL (requires cmake >= 3.2)

### Features

- separate config and data files
- reduced spectral mapping
- estimation and optimization of diffusion tensor (isotropic, axially symmetric, asymmetric, local, distributed)
- fit extended Lipari-Szabo models with selection options of BIC and AIC
- xmgrace plot of select model spectral density function

### Citation

Please cite the following paper if you use eMF.

- Prion Proteins with Pathogenic and Protective Mutations Show Similar Structure and Dynamics.
Sung-Hun Bae, Giuseppe Legname, Ana Serban, Stanley B. Prusiner, Peter E. Wright, and H. Jane Dyson
Biochemistry 2009 48 (34), 8120-8128 http://pubs.acs.org/doi/abs/10.1021/bi900923b

### Related

- http://github.com/sunghunbae/BE2


How to install?
================

### Dependencies

You need the following two packages for proper installation:

- CMake
- GNU Scientific Library (GSL)

### Make

```
$ mkdir build
$ cd build
$ cmake ..
$ make 
$ make install
```

By default, an executable binary file(```eMF```) and 
a perl script(```eMF-digest```) will be installed in ```~/bin``` directory.
You may change the destination directory defined the ```CMakeLists.txt``` file :

```
SET(CMAKE_INSTALL_PREFIX ~/bin )
```

Please make sure that ```eMF``` and ```eMF-digest``` 
is accessible from your working directory 
by adding ```~/bin``` or your destination directory in the 'path'.

How to prepare input files?
============================

You need to prepare config and data files.
These files are plain/text file and please use the following examples files
as templates.

### Config file

Default config filename is ```emf.conf``` but you may use 
different config filename by ```-c [filename]```.

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

# diffusion tensor: isotropic (-I)
# par tc    -lb 7.5 -ub   9.5 -gr 20 -cv 1e-3

# diffusion tensor: axially symmetric (-A) 
# There are two possible solutions: protate and oblate.
#     prolate (z-axis is the largest,  Dr >1)
#     oblate  (z-axis is the smallest, Dr <1)
#
#     Dpar = Dzz       (parallel)
#     Dper = Dxx = Dyy (perpendicular)
#     Dr   = Dpar/Dper = Dzz/Dxx or Dzz/Dyy
#
#     Dper = 1.0 / (2 * tc * (2.0 + Dr ))
#     At the isotropic diffusion limit (Dr=1), Dper = Dpar = 1/(6*(tc))
#
# valid range of diffusion tensor angles
# phi   : 0-360 degree
# theta : 0-180 degree
    
par tc      -lb 7.5 -ub   9.5 -gr 20 -cv 1e-3
par Dr      -lb 1.0 -ub   1.3 -gr 10 -cv 0.01
par phi     -lb 0   -ub 360   -gr 24 -cv 0.1     -fx -2.385
par theta   -lb 0   -ub 180   -gr 12 -cv 0.1     -fx  0.000

# diffusion tensor: anisotropic (-N)
# valid range of diffusion tensor angles
# phi   : 0-360 degree
# theta : 0-180 degree
# psi   : 0-360 degree
#
# par Dxx   -lb 0.008 -ub 0.03  -gr 20  -cv 1e-3 
# par Dyy   -lb 0.008 -ub 0.03  -gr 20  -cv 1e-3
# par Dzz   -lb 0.008 -ub 0.03  -gr 20  -cv 1e-3
# par phi   -lb 0.0   -ub 360   -gr 24  -cv 1e-3  -in 173.128
# par theta -lb 0.0   -ub 180   -gr 12  -cv 1e-3  -in 75.123
# par psi   -lb 0.0   -ub 360   -gr 24  -cv 1e-3  -in 82.351

# PDB COORDINATES (NECESSARY FOR NON-ISOTROPIC MODELS)
# coor [file name] [-h atom name] [-x atom name]
coor 1rx7.rot.pdb -h H -x N

# ADJUST EXPERIMENTAL ERRORS
# error [-field #] [-min # # #] [-scale # # #]
# -field #     : specifiy field or * for all data
# -min # # #   : set minimum proportioanl errors for R1,R2,NOE 
# -scale # # # : scale errors for R1,R2,NOE           

# error -field * -min 0.05 0.05 0.05
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
# i.e. use cluster 1 AND S2 > 0.7 for optimization of diffusion tensor

# FITTING AND SELECTING MODELS
#
# fitmodel [-select AIC|BIC] [-mc #] [-mctrim #]
#
# -select AIC : use Akke information criterion
# -select BIC : use Bayesian information criterion

fitmodel -select BIC
```

### Data file
```
# column 1   : residue number
# column 2   : cluster number, used for inclusion or exclusion of data
# column 3   : B0, field strength in MHz
# column 4-5 : R1 relaxation rate and error (1/s)
# column 6-7 : R2 relaxation rate and error (1/s)
# column 8-9 : {1H}-X NOE and error
# column 10  : (only for distributed) tc distribution file name

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

How to run?
===========

### Help
```
$ eMF -h
```

### Reduced spectral mapping
Reduced spectral mapping only needs ```const``` definitions in the config file.
```
$ eMF -s dhfr.dat
```
Output will be like:
```
# REDUCED SPECTRAL DENSITY
# 
# B0= 500.38 MHz, 0.87H freq.= 435.33 MHz or  2.735e+09 rad/s
# B0= 500.38 MHz,     N freq.= -50.69 MHz or -3.185e+08 rad/s
# B0= 600.13 MHz, 0.87H freq.= 522.11 MHz or  3.281e+09 rad/s
# B0= 600.13 MHz,     N freq.= -60.79 MHz or -3.820e+08 rad/s
# J : (ns/rad)
# 
# RESID   J(0)     dJ(0)    J(N)     dJ(N)  J(0.87H) dJ(0.87H) B0(MHz)
# ----- -------- -------- -------- -------- -------- -------- --------
     2 2.68e+00 6.42e-02 3.52e-01 8.02e-03 6.62e-03 9.53e-04  500.38
     2 2.73e+00 9.90e-02 2.56e-01 1.07e-02 4.53e-03 6.82e-04  600.13
     3 3.15e+00 9.46e-02 3.43e-01 9.31e-03 5.27e-03 7.13e-04  500.38
     3 2.98e+00 1.38e-01 2.46e-01 1.03e-02 4.03e-03 9.79e-04  600.13
```

### Estimating and optimizing diffusion tensor

| diffusion tensor  | parameter(s)                    | eMF command |
|-------------------|---------------------------------|-------------|
| isotropic         | tc                              | I           |
| axially symmetric | tc, Dr, phi, theta              | A           |
| anisotropic       | Dxx, Dyy, Dzz, phi, theta, psi  | N           |
| distributed       | (datafile)                      | D           |

Diffusion tensor should be determined before fitting the models.
In general, determination of diffusion tensor is an iterative process
and takes several trials. So, you can disable the ```fitmodel``` 
at the last line in the config file 
until you find a suitable diffusion tensor:

```
#fitmodel -select BIC
```

Run 
```
$ eMF -A dhfr.dat 
```

### Fitting models

Once the diffusion tensor is determined, you can disable 
estimation and optimization and enable fitmodel:
```
#estimate -cluster 1
```
```
#optimize -cluster 1 -s2 0.7
```
```
fitmodel -select BIC
```
Run
```
$ eMF -A dhfr.dat
```

For xmgrace outputs (```resid_0007.xmgr, resid_0008.xmgr, ...```):
```
$ eMF -x resid -A dhfr.dat
```
Output will be like:
```
  eMF 1.0 by Sung-Hun Bae (2008)

# data                          = dhfr.dat
#   number of lines             = 10
#   number of clusters          = 2 (0, 1)
#   number of residues          = 5
#   number of magnetic fields   = 2 (500.38, 600.13)
# configuration                 = emf.conf
#   gamma H                     = 26.7519
#   gamma X                     = -2.71
#   CSA X (ppm)                 = -172
#   bond length (A)             = 1.02
#   estimation of diffusion     = No
#   optimization of diffusion   = No
#   model selection criterion   = BIC
#   Monte Carlo simulation      = No
#   error adjusted              = Yes
#     R1 error minimum (%)      = 0, 0
#     R1 error scale            = 1, 1
#     R2 error minimum (%)      = 0, 0
#     R2 error scale            = 1, 1
#     hnNOE error minimum       = 0, 0
#     hnNOE error scale         = 1, 1
# PDB coordinate file           = 1rx7.rot.pdb
#   number of vectors           = 5
#   atom H                      = H
#   atom X                      = N
# diffusion tensor              = global axial
# spectral density function     = Lipari-Szabo extended

#    Param    Value [   Lower ...   Upper | GRD,    Step ] Convergence
#    S2f      0.000 [   0.000 ...   1.000 |  20, 5.0e-02 ] 1.0e-05
#    S2s      0.000 [   0.000 ...   1.000 |  20, 5.0e-02 ] 1.0e-05
#    te       0.000 [   0.000 ...   3.000 |  30, 1.0e-01 ] 1.0e-05
#    Rex      0.000 [   0.000 ...  10.000 |  10, 1.0e+00 ] 1.0e-05
# GF tc       8.630 [   7.500 ...   9.500 |  20, 1.0e-01 ] 1.0e-03
# GF Dr       1.195 [   1.000 ...   1.300 |  10, 3.0e-02 ] 1.0e-02
# GF phi     -2.385 [   0.000 ... 360.000 |  24, 1.5e+01 ] 1.0e-01
# GF theta    0.000 [   0.000 ... 180.000 |  12, 1.5e+01 ] 1.0e-01

RESIDUE 7        tc    8.630
RESIDUE 7        Dr    1.195
RESIDUE 7       phi   -2.385
RESIDUE 7     theta    0.000
RESIDUE 7     alpha   34.971
>+M1 S2 0.828 0.018 S2s 0.828 0.018 S2f 1.000 0.000 te  0.000  0.000 Rex 0.00 0.00 X2 1.099 dof 5 BIC 2.890 AIC 3.099
 +M2 S2 0.828 0.018 S2s 0.828 0.018 S2f 1.000 0.000 te -1.095  0.000 Rex 0.00 0.00 X2 1.099 dof 4 BIC 4.682 AIC 5.099
 +M3 S2 0.821 0.022 S2s 0.821 0.022 S2f 1.000 0.000 te  0.000  0.000 Rex 0.22 0.42 X2 0.818 dof 4 BIC 4.402 AIC 4.818
 +M4 S2 0.815 0.140 S2s 0.815 0.140 S2f 1.000 0.000 te  0.009  0.217 Rex 0.28 1.41 X2 0.277 dof 3 BIC 5.652 AIC 6.277
 +M5 S2 0.824 3.155 S2s 0.985 2.659 S2f 0.837 2.270 te  0.131 29.834 Rex 0.00 0.00 X2 0.699 dof 3 BIC 6.074 AIC 6.699
# EXP     500.38  1.6265  0.0691 10.7276  0.4488  0.7760  0.0400
# FIT M1  500.38  1.6351         10.4611          0.8057        
# EXP     600.13  1.2629  0.0436 11.1164  0.7537  0.8070  0.0880
# FIT M1  600.13  1.2784         11.1412          0.8271        
```

### Digesting selected models

Selected models are marked by ```>``` in the eMF output. 
Parameters of these models can be conveniently extracted from
the eMF output file by ```eMF-digest```:

```
$ eMF-digest dhfr-test.out
```

Output will be like:
```
# Resid  S2.val  S2.err S2s.val S2s.err S2f.val S2f.err  te.val  te.err Rex.val Rex.err      X2     dof     BIC     AIC   Model
      7   0.828   0.018   0.828   0.018   1.000   0.000   0.000   0.000   0.000   0.000   1.099   5.000   2.890   3.099       1
      8   0.836   0.110   0.836   0.110   1.000   0.000   0.011   0.216   0.660   0.970   0.552   3.000   5.927   6.552       4
      9   0.816   0.022   0.816   0.022   1.000   0.000   0.000   0.000   4.270   0.410   3.180   4.000   6.763   7.180       3
     10   0.842   0.014   0.842   0.014   1.000   0.000   0.000   0.000   0.000   0.000   1.616   5.000   3.407   3.616       1
     11   0.781   0.011   0.781   0.011   1.000   0.000   0.000   0.000   0.000   0.000   8.499   5.000  10.291  10.499       1
```
