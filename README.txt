       birkhoff-decomposition       
====================================
-- description --
collection of birkhoff decomposition algorithms:
    - birkhoff.c (classic greedy decomposition)
    - binary_bottleneck.c (bottlenck matching decomposition based on binary search)
    - beta_solstice.c (modified Solstice with variable weight-scaling factor)
    - kpack.c (novel k-Pack algorithm, computes up to k concurrent bottleneck matchings and takes best)

-- input format --
bistochastic matrices in Matrix Market format (*.mtx files)
https://math.nist.gov/MatrixMarket/formats.html

-- output format --
matrix size and permutation count on first line
factors and permutations are printed on subsequent lines (one line each)
factor and permutation pairs are separated by empty lines
permutations are 0-indexed column indices sorted by row

-- compiling / building --
./builder.sh

-- run example --
./build/algorithms/birkhoff < examples/064_4_12_70_skewed.mtx
(all programs read from stdin and write to stdout)

-- clean up --
./cleaner.sh

-- references --
Solstice by Liu et al.
https://www.cs.cmu.edu/~dga/papers/solstice-conext15.pdf

bottleneck matching algorithms in Assignment Problems (172 ff.) by Burkard, R.E., Dell'Amico, M. and Martello, S.
https://books.google.de/books?id=nHIzbApLOr0C

MatchMaker C-lib by K. Kaya, J. Langguth, I. Panagiotas and B. Uçar
https://gitlab.inria.fr/bora-ucar/matchmaker

quadsort by scandum
https://github.com/scandum/quadsort

Matrix Market I/O C-lib by the National Institue of Standards and Technology (NIST)
https://math.nist.gov/MatrixMarket/mmio-c.html 

perfect-matching-enumeration
https://github.com/SteezyE/perfect-matching-enumeration
