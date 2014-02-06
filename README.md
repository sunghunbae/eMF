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

### features

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

> By default, an executable binary file, eMF, will be installed in ```~/bin``` directory.
> You may change the destination directory defined the ```CMakeLists.txt``` file :
> 
```
SET(CMAKE_INSTALL_PREFIX ~/bin )
```

Please make sure that eMF is accessible in your working directory 
by adding ```~/bin``` or your destination directory in the 'path'.

How to use ?
============
