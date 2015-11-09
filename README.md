# Using RngStreams for parallel random number generation in C/C++ and R

Simulation studies are in many respects the ideal application for parallel processing.
They are “naturally parallel” in the sense that computations may generally be conducted
without the need for inter-processor communication or modification of the base algorithm. 
Thus, near linear speedup can be realized in many cases that occur in practice.

However, in order for a division of labor across multiple processors to be productive, 
the random number streams being used by each processor must behave as if they are **independent** in a probabilistic sense.

Although the relationship between streams will likely be much complex in situations that arise in practice, 
the message remains the same: **inter-stream dependence** can sabotage a parallel random number generation scheme.

## Methods

Hill (2010) reviews available methods for parallel random number generation, including

1. **Random spacing**: workers are initialized to randomly spaced positions on the period of the same generator by assigning different, randomly generated seeds
2. **Sequence splitting**: dividing a sequence into non overlapping contiguous blocks
3. **Cycle division**: a more involved technique for sequence splitting where the period of a generator is deterministically divided into segments
4. **Parametrization**: parameters associated with a generator are varied to produce different streams

## Generators

1. Mersenne Twister (MT) generator (**R**’s default random uniform generator) with random spacing
2. Combined Multiple recursive generator (CMRG) with cycle division (**RngStreams** package) 
3. Multiplicative lagged Fibonacci generator with parametrization (**SPRNG** package)

#### Mersenne Twister (MT) generator with random spacing 
The random spacing method carries a risk of overlap between random streams on different processors with
an unlucky sampling of seeds, although the risk is small when using generators with large periods.

MT provides one example of a long-period generator with the period 2^19937 − 1 (Matsumoto and Nishimura 1998). 
This feature along with its widespread availability has made MT a popular choice for parallel generation via random spacing.

#### Combined Multiple recursive generator (CMRG) with cycle division
A method such as cycle division comes with a guarantee that the streams will not collide, 
whereas random spacing only renders such events unlikely.

The work of Matsumoto et al. (2007) shows that many modern random number generators exhibit correlations if their initial states are chosen using another linear generator with a similar modulus. 

Issues that arisewhen using linear congruential generators for parallel random number
generation have been well documented. Random spacing and sequence splitting
may result in random number streams that exhibit undesirable dependence properties.

Many of the problems with congruential generators can be avoided by using combined
multiple recursive generators (CMRG) with cycle division. 
A MRG is roughly equivalent to a multiple recursive generator with a period that can be as large as the
product of the periods of the generators used in the combination (L’Ecuyer and Tezuka 1991).

With a good choice for the parameters (L’Ecuyer 1999), these combined generators have also
fared well when subjected to statistical tests as in L’Ecuyer and Simard (2007). 
As a case in point, the **MRG32k3a** generator that provides the backbone generator for the
**RngStreams** package developed by L’Ecuyer et al. (2001) is known to have good
theoretical and statistical properties. It is cited by Lemieux (2009), along with MT, as
a generator that can be “safely used”.

Of course, cycle division becomes an option only if it can be done efficiently. For
multiple recursive generators L’Ecuyer (1990) has shown that there are computable
matrices that can be used to advance the state of the generator to any specified point
in its associated random number stream. This feature is what makes the combination
of generators easy to use in a cycle division paradigm. The **RngStreams** package
gives an implementation of this approach in the context of its **MRG32k3a** generator.

The **RngStreams** package is quite compact with all its [source code](http://www.iro.umontreal.ca/~lecuyer/myftp/streams00/) available.

#### Multiplicative lagged Fibonacci generator with parametrization

Another noteworthy package that provides functionality for parallel random number
generation is the **SPRNG** package of Mascagni and Srinivasan (2000). In contrast to
**RngStreams**, it produces different random number sequences for the processes
via parametrization of a single type of generator.

For example, one of the package’s
featured generators is a multiplicative lagged Fibonacci generator whose states fall into 2^1007 different equivalence classes of generators (each of period 281) that provide
the different streams.

L’Ecuyer et al. (2001) in reference to the **SPRNG** package
state that it is “not supported by the same theoretical analysis for the quality and
independence of the different streams” as RngStreams. There are also practical
difficulties that arise in using **SPRNG** due to its ties to **MPI** and a rather complex
implementation that necessitates the creation and linking of a compiled library.

## Goals

1. Understand the **MRG32k3a** generator in **RngStreams** package 
2. Employ **RngStreams** in C++ programs parallelized via **OpenMP** and **MPI**
3. Employ **RngStreams** in R programs parallelized through the [parallel](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf) and 
[rstream](https://cran.r-project.org/web/packages/rstream/rstream.pdf) packages
4. Employ **.Call** in R/C interfaces parallelized via **OpenMP**
5. Implement a Monte Carlo integration
6. Implement a test of independence

## References

1. A. T. Karl, R. Eubank, J. Milovanovic, M. Reiser, and D. Young (2014). Using RngStreams for parallel random number
generation in C++ and R. *Computational Statistics*, 29(5), 1301-1320.
