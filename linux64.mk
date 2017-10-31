# $Id: linux64.mk,v 1.7 2012/11/22 11:40:26 mjlaine Exp $
#
# defines for Linux 64 bit
#
CC=gcc
F90=f90
F77=f90
AR=ar
RANLIB=ranlib
RM=rm -f
MV=mv -f
ZIP=zip
MEXEND=mexglna64
MEX=mex
MEXFLAGS=-fortran
ARFLAGS=
O=o

ifeq ($(debug),yes)
CCFLAGS=-g
MEXFLAGS+= -g
else
CCFLAGS=-O2
endif
ifeq  ($(matlab_mat),yes)
CCFLAGS += -I/usr/local/matlab/extern/include
endif

ifeq ($(f90compiler),gfortran) 
F90=gfortran
F77=gfortran
F90FLAGS= -DGFORTRAN -frecord-marker=4 -fconvert=little-endian
F77FLAGS= -DGFORTRAN -frecord-marker=4 -fconvert=little-endian
ifeq ($(debug),yes)
F90FLAGS += -g -Wall -fbounds-check 
F77FLAGS += -g -Wall -fbounds-check
else
F90FLAGS += -O2 -mtune=native
F77FLAGS += -O2 -mtune=native
endif
endif

### Inter fortran
ifeq ($(f90compiler),intel) 
F90=ifort
F77=ifort
CC=icc
F90FLAGS=-fpp -dps -vms -w90 -cm -Vaxlib
ifeq ($(debug),yes)
F90FLAGS += -g  -W1
else
F90FLAGS += -O2 -ip -funroll-loops
endif
F77FLAGS=$(F90FLAGS) -FI -w90 -w95
endif

### PGI
ifeq ($(f90compiler),pgi)
F90=pgf90
F77=pgf90
F90FLAGS=-DPGI
ifeq ($(debug),yes)
F90FLAGS += -g 
else
F90FLAGS += -O2
endif
F77FLAGS=$(F90FLAGS)
endif

### Pathscale
ifeq ($(f90compiler),pathscale)
F90=pathf90
F77=pathf90
CC=pathcc
F90FLAGS=-DPATHSCALE
ifeq ($(debug),yes)
F90FLAGS += -g 
else
F90FLAGS += -O2
endif
F77FLAGS=$(F90FLAGS)
endif


### Sun
ifeq ($(f90compiler),sun) 
F90=f95
F77=$(F90)
CC=gcc
## 
F90FLAGS= 
F77FLAGS= -f77
ifeq ($(debug),yes)
F90FLAGS += -g
F77FLAGS += -g
else
F90FLAGS += -O2
F77FLAGS += -O2
endif
endif


## Cray
ifeq ($(f90compiler),ftn) 
F90=ftn
F77=$(F90)
CC=cc
F90FLAGS= -target=linux
F77FLAGS= -target=linux
ifeq ($(debug),yes)
F90FLAGS += -g
F77FLAGS += -g
else
F90FLAGS += -O
F77FLAGS += -O
endif
endif
