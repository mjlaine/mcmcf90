# $Id: aix.mk,v 1.1 2010/09/21 10:08:26 mjlaine Exp $
#
# defines for IBM AIX (cineca)
#
CC=xlc
F90=xlf90
F77=xlf
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

# -qsuffix=f=f90:cpp=F90
F90FLAGS += -WF,-DAIX 
F77FLAGS += -WF,-DAIX

ifeq ($(debug),yes)
F90FLAGS += -g
F77FLAGS += -g
else
F90FLAGS += -O
F77FLAGS += -O
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
