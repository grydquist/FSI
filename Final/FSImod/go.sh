#!/bin/bash
rm -f a.out
rm -f d.txt
#rm -f *.mod
rm -f *i90
DBGFLAG="-ffree-line-length-none -g -fbacktrace -fcheck=all -pedantic -Wall -Wextra -W -Wno-unused-function -fopenmp"
OPTFLAG="-ffree-line-length-none -O3 -ftree-vectorize -ffast-math -funroll-loops -fomit-frame-pointer -pipe"
gfortran $DBGFLAG solmod.f90 tint.f90 elemod.f90 msh.f90 lresidelast.f90 solvemod.f90 elastmain.f90
#rm -f *.mod *.o
./a.out