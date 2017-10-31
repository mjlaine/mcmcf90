# $Id: Makefile,v 1.58 2017/08/10 07:18:33 mjlaine Exp $
## MCMC f90 library  Makefile
## Uses (and needs) GNU make
##
## make newversion vesion.h
## updates the version.h and the library minor version number
##
## Marko Laine 2001 <marko.laine@fmi.fi>
## Copyrights licensed under a MIT License.
## See the accompanying LICENSE.txt file for terms.

MCMC_MAJOR_VERSION=1
MCMC_MINOR_VERSION=2
MCMC_VERSION=$(MCMC_MAJOR_VERSION).$(MCMC_MINOR_VERSION)

# Try to guess the system
sys := $(shell uname)$(shell uname -m)
# sys := aix
# sys := win32

# default is Linux64 (or win32) with gfortran
system=linux64
system=win32
f90compiler=gfortran

ifeq ($(sys),Linuxi686)
 system=linux
 f90compiler=gfortran
endif
ifeq ($(sys),Linuxx86_64)
 system=linux64
# f90compiler=gfortran
# f90compiler=intel
endif
ifeq ($(sys),Linuxia64)
 system=linux
# f90compiler=intel
endif
ifeq ($(sys),Darwini386)
 system=maci
 f90compiler=gfortran
endif
ifeq ($(sys),Darwinx86_64)
 system=maci
 f90compiler=gfortran
endif
ifeq ($(sys),aix)
 system=aix
 f90compiler=xlf90
endif

# produce debug code?
debug=no
#debug=yes

### no user modification below this, but check *.mk files ###

# read system specific compiler flags from a separate file
ifeq ($(system),win32)
  include win32.mk
endif
ifeq ($(system),mac)
  include mac.mk
endif
ifeq ($(system),maci)
  include maci.mk
endif
ifeq ($(system),linux)
  include linux.mk
endif
ifeq ($(system),linux64)
  include linux64.mk
endif
ifeq ($(system),aix)
  include aix.mk
endif

# main source files
LIBSRC1 = mcmc_main.F90 mcmc.F90 mcmcinit.F90 matutils.F90 mcmcprec.F90 \
 mcmcrand.F90 matfiles.F90 priorfun.f90 ssfunction0.f90 checkbounds0.f90 \
 mcmcrun1.F90 signalqq.c dchud.f dchdd.f initialize.F90 dump.F90 \
 ssfunction_er0.f90

# files included in mcmc.F90
LIBSRC2 = MCMC_DRAM.F90 MCMC_aux.F90 MCMC_init.F90 MCMC_userfun.F90 \
          MCMC_adapt.F90  MCMC_run.F90 MCMC_run1.F90 \
          MCMC_run_scam.F90 MCMC_run_ram.F90 MCMC_run_er.F90 MCMC_run1_er.F90 \
          MCMC_signal_handler.F90 MCMC_dump.F90

# for the src zip file
ALLSRC =  $(LIBSRC1) $(LIBSRC2)

LIBOBJ := $(patsubst %.f,%.$(O),$(filter %.f,$(LIBSRC1))) \
          $(patsubst %.F,%.$(O),$(filter %.F,$(LIBSRC1))) \
          $(patsubst %.f90,%.$(O),$(filter %.f90,$(LIBSRC1))) \
          $(patsubst %.F90,%.$(O),$(filter %.F90,$(LIBSRC1))) \
          $(patsubst %.c,%.$(O),$(filter %.c,$(LIBSRC1)))

# Test files
TSTFILES=testcases/Makefile testcases/mcmcrun.f90 \
	testcases/mcmcinit.nml testcases/data.dat testcases/mcmcpar.dat \
        testcases/mcmccov.dat testcases/mcmcsigma2.dat testcases/README.txt

# Matlab files
MFILES = 

# source code and other files needed
SRC = $(ALLSRC) lapack_inc.h external_inc.h  \
      version.h .version mkversion.sh mcmcinit.nml \
      Makefile linux.mk linux64.mk win32.mk mac.mk maci.mk aix.mk \
      $(MFILES) Readme.txt Readme_code.txt INSTALL.txt LICENSE.txt \
      $(TSTFILES)

.SUFFIXES: .a .$(O) .l$(O) .mod .f90 .F90 .F .f

all: libmcmcrun.a

zip: mcmcf90src.zip

mcmcf90src.zip: $(SRC)
	$(ZIP) -9 mcmcf90src.zip $(SRC)
clean:
	$(RM) *.$(O) *.l$(O) core *~ *.mod
realclean: clean
	$(RM) *.o *.obj *.lo *.il *.lobj *.a *.lib *.mod

standalone: libmcmcrun.a

libmcmcrun.a: $(LIBOBJ)
	$(AR) ruv libmcmcrun.a $(LIBOBJ)
	$(RANLIB) libmcmcrun.a

ifeq ($(system),aix)
defines_aix = $(defines:-D%=-WF,-D%)
.f90.$(O):
	$(F90) $(F90FLAGS) $(defines_aix) -c $<
.F90.$(O):
	$(F90) $(F90FLAGS) $(defines_aix) -c $<
.f.$(O):
	$(F77) $(F77FLAGS) $(defines_aix) -c $<
.F.$(O):
	$(F77) $(F77FLAGS) $(defines_aix) -c $<
.c.$(O):
	$(CC) $(CCFLAGS) $(defines_aix) -c $<
else

.f90.$(O):
	$(F90) $(F90FLAGS) $(defines) -c $<
.F90.$(O):
	$(F90) $(F90FLAGS) $(defines) -c $<
.f.$(O):
	$(F77) $(F77FLAGS) $(defines) -c $<
.F.$(O):
	$(F77) $(F77FLAGS) $(defines) -c $<
.c.$(O):
	$(CC) $(CCFLAGS) $(defines) -c $<
endif

# Cancel .mod.o rule
%.o : %.mod
# 

newversion:
	@. mkversion.sh > .ver
	@mv -f .ver .version

version.h: .version
	@/bin/echo -n "character(len=*), parameter :: " > .ver
	@/bin/echo -n Mcmc_Code_Version = \"$(MCMC_VERSION).`cat .version` >> .ver
	@/bin/echo ' -- '`date +%Y-%m-%d`'"' >> .ver
	@/bin/echo "character(len=*), parameter :: Mcmc_Lib_Compile_Date = __DATE__" >> .ver
	@mv -f .ver $@

# dependencies
mcmc_main.$(O): mcmc.$(O)
mdstmcmc.$(O): matutils.$(O) mcmc.$(O)
matutils.$(O): mcmcprec.$(O) mcmcrand.$(O) matfiles.$(O)
mcmc.$(O): mcmcinit.$(O) matutils.$(O) mcmcprec.$(O) mcmcrand.$(O) version.h $(LIBSRC2) mcmcrun1.$(O)
mdstexpanneal.$(O): mcmc.$(O) matutils.$(O) version.h
priorfun.$(O): mcmc.$(O)
mcmcrun1.$(O): matutils.$(O)
mcmcrand.$(O): mcmcprec.$(O)
mcmcinit.$(O): mcmcprec.$(O)
