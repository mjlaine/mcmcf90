## Makefile for mcmcrun project

PROJECT=mcmcrun

SRC=mcmcrun.F90

FC = gfortran
FFLAGS = -O3 -mtune=native -frecord-marker=4 -fexternal-blas
#FFLAGS = -g -Wall -fbounds-check -frecord-marker=4 -fexternal-blas
# In Linux, we need lapack and blas
LIBS = -lmcmcrun -llapack -lblas
# In Mac blas and lapack are provided in Accelerate framework
LIBS = -lmcmcrun -framework Accelerate

# cray
#FC = pgf90
#FFLAGS = -g
#LIBS = -lmcmcrun -lacml

# cray with intel
#FC = ftn
# Intel
#FFLAGS = -g -mkl
#LIBS = -lmcmcrun
# tai pgi
#FFLAGS = -g
#LIBS = -lmcmcrun -lacml


# Need to know where libmcmcrun.a and the *.mod files are
FFLAGS += -I.. -L..

all: $(PROJECT)

mcmcrun: $(SRC)
	$(FC) $(FFLAGS) -o $(PROJECT) $(SRC) $(LIBS)

run: mcmcrun
	./mcmcrun

clean:
	rm -f *.o *~ core $(PROJECT)

realclean:
	rm -f *.o *~ core *seed.dat chain.dat chain.mat \
        sschain.dat sschain.mat s2chain.dat s2chain.mat $(PROJECT) \
	mcmccovf.dat mcmcparf.dat mcmcsigma2f.dat mcmcinit_final.nml \
        mcmcmean.dat
	rm -rf *.dSYM
