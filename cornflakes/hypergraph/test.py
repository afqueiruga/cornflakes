from pyHypergraph import Hypergraph
import numpy as np

h = Hypergraph(2)
h.Push_Edge(np.array([1,2],dtype=np.intc))

import mylibrary
import particle_placers
h = Hypergraph()
x = particle_placers.init_grid(10,10,[0.0,0.0],[12.0,2.0],[0.0,10.0],0.5)
#np.array([[ 1.0,2.0], [2.0,2.0], [10.0,1.0]])
mylibrary.Build_Particle_Graph(h.hg, x, 1.5)

v = h.view()

import graphio
graphio.write_graph("foo.vtk",h,x)
