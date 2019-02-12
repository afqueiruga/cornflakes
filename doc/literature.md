Tools:

PETSc: For MPI. Built in C. Has all the libraries. Traditional.
Chapel: PGAS language. Only does grid-array partitions.
UPC: NO!

RAJA: https://github.com/LLNL/RAJA does loop abstractions. Define a kernel, and it translates it to CUDA/pthreads/vectorized/etc
Kokkos: https://github.com/kokkos/kokkos similar to RAJA

Idea: Popcorn should target RAJA or Kokkos, or another languayge. cornflakes:assemble needs to be aware of this. Requires rewrite into C++. This abstractions is hidden _inside_ assemble; hypergraphs/dofmaps do not need to know about GPUs.