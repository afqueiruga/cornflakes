Philosophy
==========

These are the key points to the underlying philosophy dictating the architecture of this library.

This program is performming operations on large chunks of data to produce residual vectors, load vectors, stiffness matrices, potential scalars, output fields, etc. The output will be used to solve an ODE, PDE, Minimization problem by some other solver. This library will _not_ provide a matrix solver, minimizer, timestepper, etc. This library will make it easy to link to such a solver and facilitate assembling whatever it requires.

The fundamental unit of the program is the _kernel_ calculation. This represents the smallest unit of computational work to be performed. It is a common paradigm in scientific computing. (The nomenclature is taken from GPGPU progarmming.) Examples of kernels are the finite difference stencil, the local-element calculation in FEM, an image convolution filter, a pair-force interaction, etc.

Most scientific numerical discretizations and computations can be reprsented by a _hypergraph_. The _vertices_ of the graph represent a location for degrees of freedom (DOFs) to be associated with. Each kernel calculation takes in a few of these vertices to produce a small output. The set of these few vertices is a single _hyperedge_. Every hyperedge in a hypergraph represents one call of the kernel calculation on the data associated with the vertices.

Data is not tied to a vertex. Multiple different fields of data can be associated with a vertex of different characters. States of data are just pointers to long chunks of memory filled with floating point numbers. They can be swapping around willy-nilly on the whim of the user. Nothing more, nothing less. The association between vertices in the global graph and indices to actual data is performed by a Dofmap object. There is one Dofmap per type of field. I.e., there is a Dofmap for associating vertices with x,y,z, and another for p,rho, alpha, etc. An example of this usage is:

V1 ---> DofmapVector ---> [dof1, dof2, dof3]
xdata[ [dofs] ] ---> [ x1, x2, x3 ]
ydata[ [dofs] ] ---> [ y1, y2, y3 ]
V1 ---> DofmapTemperature ---> [ dof4 ]
Tdata[ [dof4] ] ---> [ T1 ]
V4 ---> DofmapPressure ---> [ dof1 ]
Pdata[ [dof1] ] ---> [ P1 ]

Hyperedges _cannot_ be labeled or sorted by the user. Hyperedges _cannot_ have data associated with them. To associate data with a hyperedge, create a new vertex. Helper graphing routines will be provided. Hypergraphs _can_ have hyperedges with multiple lengths (but the kernel must handle variable length input). The vertices in a hyperedge _can_ be ordered and the order will be preserved.

When partitioning, each processor gets a disjoint subset of the hyperedges in the total graph. Vertices are distributed as needed, overlapping ghost vertices are handled, and new Dofmaps are made.

The mapping from the output kernel calculation to the output large-residual-matrix-global-vector-thingy is handled by the same Vertex->Dofmap->index data structure. The output Dofmap will depend on vertices of the hyperedge, just like the input-data fetching. If the user does not have vertices that fit the neccessity (e.g., trying to calculate a length on a pair-bond and save the data to the bond), THE USER MUST ADD A VERTEX THAT REPRESENTS THE TOPOLOGY OF THE OUTPUT. _THE USER IS NOT ALLOWED TO ORDER THE HYPEREDGES_! This information must be stored on a vertex, proving the de facto ordering information. Deal with it.
Dofmap objects can be recycled to save on memory and reduce lines of code.

There can be different types and sizes of hyperedges associated with different computational kernels. They can be mixed-and matched into one large hypergraph. For simplicity, each kernel is associated with _one_ hypergraph, and it is up to the user to juggle different hypergraphs with one common vertex-naming scheme. For partitioning, a union of all the kernel-hypergraphs is taken, and then the processors re-divy up the subsets. 

For the most part, every piece of data is associated with a Vertex. The Big Exception is _global data_. It is allowed and handled as a fringe case when the kernel says the input field is associated with vertex position "-1" in the hyperedge. We could define a "global" vertex, but that'll just mess up the partitioning.
