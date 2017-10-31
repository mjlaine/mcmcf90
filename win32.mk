# $Id: win32.mk,v 1.8 2015/01/29 11:31:37 mjlaine Exp $
#
# defines for VISUAL FORTRAN
#
CC=cl.exe
F90=df.exe
F77=df.exe
AR=lib.exe
RM=del
MV=move.exe
ZIP=zip.exe
MEXEND=dll
MEX=mex

ifeq ($(debug),yes) 
# debug
CCFLAGS=/nologo /MLd /W3 /Gm /GX /ZI /Od /I "c:\matlabr11\extern\include" \
         /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /FD /GZ /c 
F90FLAGS=/check:bounds /compile_only /debug:full /fpp /nologo \
        /traceback /warn:argument_checking /warn:nofileopt 

F77FLAGS=$(F90FLAGS)
ARFLAGS=
else
CCFLAGS=/nologo /ML /W3 /GX /O2 /I "c:\matlabr11\extern\include" \
        /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /c 
F90FLAGS=/compile_only /fpp /nologo /warn:nofileopt /tune:p5 \
         /optimize:5 /math_library:fast
endif

O=obj

# gfortran with win32
ifeq ($(f90compiler),gfortran) 
O=o
AR=ar
RANLIB=ranlib
RM=rm -f
MV=mv -f
CC=gcc
F90=gfortran
F77=gfortran
F90FLAGS=-DGFORTRAN -U_WIN32 -DWIN32
F77FLAGS=-DGFORTRAN -U_WIN32 -DWIN32
CCFLAGS=
ifeq ($(debug),yes)
F90FLAGS += -g -Wall -fbounds-check
F77FLAGS += -g -Wall -fbounds-check
else
F90FLAGS += -O2
F77FLAGS += -O2
endif
endif
