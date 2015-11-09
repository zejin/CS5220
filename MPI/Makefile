# ===
# Compilers

CC=gcc
MPICC=mpicc

# ===
# Flags

CFLAGS=-std=c99 -g 
OPTFLAGS=-O3 -march=native
CFLAGS+=$(OPTFLAGS)

MPI_CFLAGS=$(CFLAGS)
OMP_CFLAGS=$(CFLAGS) -fopenmp

# ===
# Rules

bar: bar.c RngStream.c
	$(MPICC) $(MPI_CFLAGS) $^ -o $@

.PHONY: clean

clean:
	rm -f bar *.o *.optrpt *~