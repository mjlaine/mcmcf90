Readme for mcmcf90 source code

Short explanation of the source files included.

mcmc_main.F90: 'main' computational program as a subroutine, to be
                called by the user.

ssfunction0.f90: dummy ssfunction to be replaced by the user at the
                 link time.

priorfun.f90: the default prior function, that might be replaced by the
              used at link time.

modules:

mcmcinit.mod

 mcmcinit.F90: mcmcinit.nml namelist definitions

mcmcmod.mod

 mcmc.F90: main mcmc computational module, includes the following files:

 MCMC_init.F90
 MCMC_run.F90
 MCMC_run_er.F90
 MCMC_run1.F90
 MCMC_run1_er.F90
 MCMC_run_scam.F90
 MCMC_run_ram.F90
 MCMC_adapt.F90
 MCMC_DRAM.F90
 MCMC_aux.F90
 MCMC_dump.F90
 MCMC_userfun.F90
 MCMC_signal_handler.F90

matutils.mod

 matutils.F90: matrix utility routines not directly related to mcmc

mcmcrand.mod

 mcmcrand.F90: random number generation

matfiles.mod

 matfiles.F90: code for writing mat v4 files

mcmcprec.mod

 mcmcprec.F90: variable precision and other constants


Copyrights licensed under a MIT License.
See the accompanying LICENSE.txt file for terms.

Two files dchdd.f and dchud.f are copied from the LINPACK library,
http://www.netlib.org/linpack/.
