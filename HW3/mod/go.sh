#!/bin/bash
rm -f a.out
rm -f *.mod
DBGFLAG="-ffree-line-length-none -g -fbacktrace -fcheck=all -pedantic -Wall -Wextra -W -Wno-unused-function -fopenmp"
OPTFLAG="-ffree-line-length-none -O3 -ftree-vectorize -ffast-math -funroll-loops -fomit-frame-pointer -pipe"
gfortran $DBGFLAG elemod.f90 msh.f90 lstiffmod.f90 gstiffmod.f90 solvemod.f90 main.f90
rm -f *.mod *.o
./a.out