# Hydrus-1D-for-Linux

This repository is based on [another one that adapted vanilla Hydrus-1D source code to linux-based OS](https://github.com/AgriHarmony/HYDRUS-1-D-gfortran). Folder `src` contains sources from this repository and a makefile, adjusted to a bit newer version of source code delivered by Professor Jirka Simunek.

Base source code is based on Compaq Visual Fortran compiler that used Microsoft libraries, which are not present in modern gfortran compiler. Therefore it has been a little bit adapted by repository mentioned above.

The sources can be compiled with `gfortran` and `make`. In order to do so, one must call `make` command inside `src/` folder.
