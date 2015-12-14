# ===
# Compilers

CC=icc
MPICC=mpiicc

# ===
# Flags

CFLAGS=-std=c99 -g 
OPTFLAGS=-O3 -no-prec-div -xcore-avx2 -ipo \
	 -qopt-report=5 -qopt-report-phase=vec
CFLAGS+=$(OPTFLAGS)

MPI_CFLAGS=$(CFLAGS)
OMP_CFLAGS=$(CFLAGS) -openmp

# ===
# Rules

sim: sim.c RngStream.c
	$(MPICC) $(MPI_CFLAGS) $^ -o $@ -lm

foo: foo.c RngStream.c
	$(MPICC) $(MPI_CFLAGS) $^ -o $@ -lm

bar: bar.c RngStream.c
	$(MPICC) $(MPI_CFLAGS) $^ -o $@

.PHONY: clean

clean:
	rm -f bar *.o *.optrpt *~