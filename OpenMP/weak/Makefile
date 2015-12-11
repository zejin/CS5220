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

sim: sim.c RngStream.c
	$(CC) $(OMP_CFLAGS) $^ -o $@ -lm

foo: foo.c RngStream.c
	$(CC) $(OMP_CFLAGS) $^ -o $@ -lm

bar: bar.c RngStream.c
	$(CC) $(OMP_CFLAGS) $^ -o $@

.PHONY: clean

clean:
	rm -f bar *.o *.optrpt *~