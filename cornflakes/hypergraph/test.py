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


def make_dofmap(h, ke, Nnode):
    # Make the dofmap for where the kernel output will be
    # assembled into. It will look like
    # G1 ... E1 E3 ... N1 N2 N3
    # for now...
    # TODO: Do I need to give edges global indices??
    
    v = h.view()
    outlen = 0
    dofs = mylibrary.dofArray_frompointer(ke.outp)
    
    for i in xrange(ke.noutp):
        if dofs[i].loc == mylibrary.LOC_NODE:
            outlen+=dofs[i].len*h.hg.l_edge
        else:
            outlen+=dofs[i].len
    
    outmap = np.zeros((v.shape[0], outlen),dtype=np.intc)
    # Loop over the edges
    for i,e in enumerate(v):
        # Loop over the dofs
        off = 0
        for d in xrange(ke.noutp):
            if dofs[d].loc == mylibrary.LOC_NODE:
                for j,n in enumerate(e):
                    outmap[i,off+dofs[d].len*j:dofs[d].len*(j+1)] = range(n*dofs[d].len,(n+1)*dofs[d].len)
                off+= len(e)*dofs[d].len
            else:
                pass
    return outmap


pks=mylibrary.cvar.particle_kernel_strct
dofmap = make_dofmap(h, pks, X.shape[1])

R = np.zeros(X.shape,dtype=np.double)
vel = X.copy()



mylibrary.assemble_vector_np(np.array([1,2,3],dtype=np.double), pks,h.hg, [X,vel])


import graphio
graphio.write_graph("foo.vtk",h,X, {'x':X,'v':vel})

