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

mpi: driverLTM_mpi.cpp RngStream.cpp RngStreamSupp.cpp funclass.cpp
	$(MPICC) $(MPI_CFLAGS) $^ -o $@

# mpi: ranMPI.o RngStream.o RngStreamSupp.o
# 	$(MPICC) $(MPI_CFLAGS) -o $@ $^

# ranMPI.o: ranMPI.cpp
# 	$(MPICC) $(MPI_CFLAGS) -c $<

# RngStream.o: RngStream.cpp
# 	$(MPICC) $(MPI_CFLAGS) -c $<

# RngStreamSupp.o: RngStreamSupp.cpp
# 	$(MPICC) $(MPI_CFLAGS) -c $<

.PHONY: clean

clean:
	rm -f *.o *.optrpt