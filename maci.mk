# $Id: maci.mk,v 1.7 2009/09/21 12:59:29 mjlaine Exp $
#
# defines for Intel Mac
#
CC=gcc
F90=gfortran
F77=gfortran
AR=ar
RANLIB=ranlib
RM=rm -f
MV=mv -f
ZIP=zip
MEXEND=mexmac
MEX=mex
ARFLAGS=
O=o

ifeq ($(debug),yes)
CCFLAGS=-g
else
CCFLAGS=-O2
endif
ifeq  ($(matlab_mat),yes)
CCFLAGS += -I/Applications/MATLAB73/extern/include
endif

ifeq ($(f90compiler),gfortran) 
F90=gfortran
F77=gfortran
F90FLAGS=-DGFORTRAN -frecord-marker=4 -fconvert=little-endian
F77FLAGS=-DGFORTRAN -frecord-marker=4 -fconvert=little-endian
ifeq ($(debug),yes)
F90FLAGS += -g -Wall -Wconversion -fbounds-check -Wno-unused-parameter
F77FLAGS += -g -Wall -Wconversion -fbounds-check -Wno-unused-parameter
else
F90FLAGS +=  -O3 -ftree-vectorize -mtune=native \
            -fexternal-blas
F77FLAGS +=  -O3 -ftree-vectorize -mtune=native \
            -fexternal-blas
endif
endif
