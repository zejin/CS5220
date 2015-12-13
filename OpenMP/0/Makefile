# ===
# Compilers

CC=icc
MPICC=mpiicc

# ===
# Flags

CFLAGS=-std=c99 -g 
OPTFLAGS=-ipo -qopt-report=5 -qopt-report-phase=vec
CFLAGS+=$(OPTFLAGS)

MPI_CFLAGS=$(CFLAGS)
OMP_CFLAGS=$(CFLAGS) -openmp

# ===
# Rules

sim: sim.c RngStream.c
	$(CC) $(OMP_CFLAGS) $^ -o $@

foo: foo.c RngStream.c
	$(CC) $(OMP_CFLAGS) $^ -o $@

bar: bar.c RngStream.c
	$(CC) $(OMP_CFLAGS) $^ -o $@

.PHONY: clean

clean:
	rm -f bar *.o *.optrpt *~