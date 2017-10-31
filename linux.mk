# $Id: linux.mk,v 1.9 2012/11/21 11:37:39 mjlaine Exp $
#
# defines for Linux
#
CC=gcc
F90=f90
F77=f90
AR=ar
RANLIB=ranlib
RM=rm -f
MV=mv -f
ZIP=zip
MEXEND=mexglx
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
F90FLAGS=-mtune=native -DGFORTRAN -frecord-marker=4 -fconvert=little-endian
F77FLAGS=-mtune=native -DGFORTRAN -frecord-marker=4 -fconvert=little-endian
ifeq ($(debug),yes)
F90FLAGS += -g -Wall -fbounds-check 
F77FLAGS += -g -Wall -fbounds-check
else
F90FLAGS += -O2
F77FLAGS += -O2
endif
endif

ifeq ($(f90compiler),vast) 
# Vast f90 compiler
F90=f90vast
F77=f77
F90FLAGS=-m486
F77FLAGS=-m486
ifeq ($(debug),yes)
F90FLAGS += -g
F77FLAGS += -g
else
F90FLAGS += -O2
F77FLAGS += -O2
endif
endif


### Inter fortran
ifeq ($(f90compiler),intel) 
F90=ifc
F77=ifc
F90=ifort
F77=ifort
CC=icc
F90FLAGS=-tpp6 -axWMK -fpp -dps -vms -w90 -cm -Vaxlib
# P4
F90FLAGS=-tpp7 -fpp -dps -vms -w90 -cm -Vaxlib
# Xeon
F90FLAGS=-mtune=pentium4 -ssp  -parallel -tpp7 -fpp -dps -vms -w90 -cm -Vaxlib
ifeq ($(sys),Linuxia64)
# Itanium
# F90FLAGS=-mtune=itanium -parallel -fpp -dps -vms -w90 -cm -Vaxlib
F90FLAGS= -parallel -fpp -dps -vms -w90 -cm -Vaxlib
endif
ifeq ($(debug),yes)
F90FLAGS += -g  -W1
else
F90FLAGS += -O3 -ip -funroll-loops
endif
F77FLAGS=$(F90FLAGS) -FI -w90 -w95
endif

### absoft Pro Fortran
ifeq ($(f90compiler),absoft) 
F90=f90
F77=$(F90)
CC=gcc
## compatibility flags with gcc and g77
F90FLAGS= -YEXT_NAMES=LCS -s -B108 -YCFRL=1
F77FLAGS= -f -s -B108 -N90
F90FLAGS= -ffree
F77FLAGS= 
ifeq ($(debug),yes)
F90FLAGS += -g -m1
F77FLAGS += -g -m1
else
F90FLAGS += -O -B101
F77FLAGS += -O -B101
endif
F90FLAGS += -DABSOFT
endif


### Sun
ifeq ($(f90compiler),sun) 
F90=f95
F77=$(F90)
CC=gcc
MEXEND=mexsol
## 
F90FLAGS= 
F77FLAGS= -f77
ifeq ($(debug),yes)
F90FLAGS += -g
F77FLAGS += -g
else
F90FLAGS += -O4
F77FLAGS += -O4
endif
endif
