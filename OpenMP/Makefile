# ===
# Compilers

CC=icpc
MPICC=mpiicpc

# ===
# Flags

CFLAGS=-std=c++14 -g 
OPTFLAGS=-O3 -no-prec-div -xcore-avx2 -ipo -qopt-report=5 -qopt-report-phase=vec
CFLAGS+=$(OPTFLAGS)

MPI_CFLAGS=$(CFLAGS)
OMP_CFLAGS=$(CFLAGS) -openmp

# ===
# Rules

LTM_omp: driverLTM_openmp.cpp RngStream.cpp RngStreamSupp.cpp funclass.cpp
	$(CC) $(OMP_CFLAGS) $^ -o $@

# omp: ranOpenMP.o RngStream.o
# 	$(CC) $(OMP_CFLAGS) -o $@ $^

# ranOpenMP.o: ranOpenMP.cpp
# 	$(CC) $(OMP_CFLAGS) -c $<

# RngStream.o: RngStream.cpp
# 	$(CC) $(OMP_CFLAGS) -c $<

.PHONY: clean

clean:
	rm -f *.o *.optrpt