#
# defines for Mac
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
CCFLAGS += -I/usr/local/matlab/extern/include
endif

ifeq ($(f90compiler),gfortran) 
F90=gfortran
F77=gfortran
F90FLAGS=-DGFORTRAN -frecord-marker=4 -fconvert=little-endian
F77FLAGS=-DGFORTRAN -frecord-marker=4 -fconvert=little-endian
ifeq ($(debug),yes)
F90FLAGS += -g -fbounds-check
F77FLAGS += -g -fbounds-check
else
F90FLAGS += -mtune=G4 -mcpu=G4 -O2 -maltivec -ftree-vectorize
F77FLAGS += -mtune=G4 -mcpu=G4 -O2 -maltivec -ftree-vectorize
endif
endif
