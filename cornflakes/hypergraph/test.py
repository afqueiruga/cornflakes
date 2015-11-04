from pyHypergraph import Hypergraph
import numpy as np

h = Hypergraph(2)
h.Push_Edge(np.array([1,2],dtype=np.intc))

import mylibrary
import particle_placers
h = Hypergraph()
X = particle_placers.init_grid(11,11,[0.0,0.0],[10.0,0.0],[0.0,10.0],0.0)
#np.array([[ 1.0,2.0], [2.0,2.0], [10.0,1.0]])
mylibrary.Build_Particle_Graph(h.hg, X, 1.5)


pks=mylibrary.cvar.particle_kernel_strct

from dofmap import make_dofmap

dofmap = make_dofmap(h, pks, X.shape[1])

R = np.zeros(X.shape,dtype=np.double)
vel = X.copy()



mylibrary.assemble_vector_np(np.array([1,2,3],dtype=np.double), pks,h.hg, [X,vel])


import graphio
graphio.write_graph("foo.vtk",h,X, {'x':X,'v':vel})

