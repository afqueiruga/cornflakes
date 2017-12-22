Cornflakes
==========

2015-2017 Alejandro F Queiruga

Lawrence Berkeley Lab

Intro
-----

Cornflakes is a simple library for assembling
vectors and matrices from the kernels created by the accompanying
Domain Specific Language popcorn. In a sense, cornflakes
is the runtime to the popcorn DSL.
It's meant for super-quick prototyping for small problems
allowing for interactivity and inspection of the matrices,
etc.


Computational science codes are often broken down into repetitive applications of a **kernel** to a small chunk of data, e.g., a finite difference stencil applied nine nodes. Cornflakes abstracts the organization of the individual operators into a **hypergraph**. Each **hyperedge** of the graph represents one instance of applying the kernel calculation to the data.
Every **hypervertex** in the graph represents a location with associated data., e.g., the fields values on the finite difference grid points. The edge is the set of vertices containing data required to calculate the kernel. The value of the hypervertex is used to index into a global vector using a **DofMap** for each of the arguments of the kernels; e.g. vertex 2000 maps to degrees of freedom `{6000,6001,6002}` in the global vector. A large class of computations can be described by this abstraction, such as
- Finite Element meshes - the edge is an element
- Finite difference grids - the edge is a stencil
- Particle methods - the edge is a neighbor list

Or, in pictures,  
![hypergraphs](doc/figures/hypergraphs.pdf)  
 The result of each edge computation is a contribution to the global system whose local results are put together in the "big A" assembly:  
![bigA](
  http://latex.codecogs.com/gif.latex?\dpi{120}&space;\large&space;\mathbf{K}=\underset{\mathtt{edge}\in\mathcal{H}}&space;{\operatorname{\raisebox{-5pt}{\mbox{&space;\Huge&space;\textsf{\textbf{A}}}}}}&space;\mathbf{k}\left%28u\left[\mathtt{edge}\right]\right%29
  )  
Cornflakes is designed to be run in parallel. (But, like many things,
I haven't gotten to that part yet. I might just rewrite it
all in Julia first.)  The hypergraph of the problem is used to distribute the computations in parallel by assigning edges to processors using the overlap of vertices to minimize the communication. The parallel vectors can then be distributed figuring out which vertices ended up where.

Philosophy
----------

1. The user must retain control of main.
2. Agnostic to how the kernel was written. It just has to have
the kernel struct. The first DSL was already replaced (kerngen).
3. The code generation should run pre-compile time, not just-in-time.
It gets expensive. Some of my kernels take 5 minutes to run. The
production code should be given already compiled shared object files.

The main idiom is that a single calculation  is at a minimum two files:  

1. method\_pop.py: the popcron file specifying one or more kernels
2. simulation.py: the actual calculation that imports the husk,
uses cornflakes to describe a hypergraph and assemble a result by
applying the kernels to the hypergraph, and uses another library
to solve the matrices of your choice.

Motivations
-----------

Cornflakes is heavily inspired by FEniCS. (If you just need to do FEM, use FEniCS!
Cornflakes isn't ready for you.) My dissertation involved extensive modifications
to FEniCS in both Dolphin and FFC to express the multiphysics, nonlinear contacts
between elements. I managed to express the integrals in the DSL but saw the limitations
to express anything further. With dolphin, I had to copy-and-paste a lot of code, and
felt that I could make it more general.
I was further working on coupling Discrete Elements to FEniCS for particle-laden flow,
but I realized there must be a better way.

Cornflakes and popcorn were first written to perform the studies in
"Numerical experiments on the convergence properties of state-based peridynamic laws...",
where I needed a new DSL system to describe, automatically derive, and differentiate
general Peridynamic calculations.

So far I've cornflakes and popcorn to program:

1. Finite Elements
1. Peridynamics
1. Bonded Particles
1. Reproducing Kernel Particle Method
1. Fully coupled Peridynamics + 2D FEM + 1D dynamically mesh FEM all in cornflakes
1. FEM coupled to a linked finite difference code


Similarities to other packages
------------------------------

The design of cornflakes is similar to the following libraries:

1. [FEniCS](fenicsproject.org) - difficult to generalize past finite elements
1. [MOOSE](http://mooseframework.org) - No control over main; no DSL
1. [PyOp2](https://github.com/OP2/PyOP2) - Assembler only. It seems great, but it was developed in parallel to cornflakes! It's only been used for FFC kernels.


A distinction to TensorFlow is what graph entities represent. Each Hypergraph represents
kernel calculations that _can happen in parallel_, with the graph connectivity representing
how small bits of data overlap. The purpose of the cornflakes description is to facilitate
communication-minimizing graph partitioning across HPC systems. E.g., the use case is to
distribute tens of millions of calculations representing elements of a mesh across an HPC system,
which will have to be evaluated millions of times. In TensorFlow, the graph entities are chained
together in a Directed Acyclic Graph, which describes _the interdependencies for the ordering_ for
the calculations. The concepts are not mutually exclusive. (On another goal, I'm currently
trying to write a TensorFlow node for calling cornflakes assemblies.)



Example
-------


Requirements
------------

Popcorn in path, scipy, numpy. Output is to vtk files.

Installing
--------

1. mkdir build
2. cd build
3. ccmake ..
4. set CMAKE_INSTALL_PREFIX to somewhere local
5. Set USE_LIS or USE_PETSC
6. make
7. make test
8. make install

### OSX Specific:

Things on OSX are tricky. You'll need to know where to link things:

1. set CMAKE_C_COMPILER and CMAKE_CXX_COMPILER gcc-mp-5
2. set SWIG_EXECUTABLE:
/opt/local/bin/swig
3. Set EXTRA_INCLUDES:
port contents numpy | grep include
4. Set PYTHON_INCLUDE_DIR
port contents python27 | grep include
5. Set PYTHON_LIBRARY
port contents python27 | grep libpython2.7
