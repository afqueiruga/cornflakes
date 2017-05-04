import cornflakes_library as cflib
import numpy as np


# interface class
#class Dofmap():
#    def __init__(self):
#        self.dm = cflib.dofmap_t()
#        # NO VTABLE SET
#    def Get(self,V):
#        return cflib.Dofmap_Get_np(self.dm,V)
#    def Max_Len(self):
#        return cflib.Dofmap_Max_Len(self.dm)
#    def Get_List(self,Vs):
#        return cflib.Dofmap_Get_List_np(self.dm, Vs)
#    def __del__(self):
#        cflib.Dofmap_Destroy(self.dm)
#class Dofmap_Strided(Dofmap):
#    def __init__(self, stride, offset=0):
#        self.dm = cflib.dofmap_t()
#        cflib.Dofmap_Strided(self.dm, stride,offset)
#class Dofmap_Tabled(Dofmap):
#    def __init__(self, stride, table, offset=0):
#        self.dm = cflib.dofmap_t()
#        cflib.Dofmap_Tabled(self.dm, table.reshape((table.size/stride,stride)), int(offset) )
#class Dofmap_From_Vertices(Dofmap):
#    def __init__(self, stride, vertices, offset=0):
#        self.dm = cflib.dofmap_t()
#        start = int(vertices.min())
#        stop =  int(vertices.max())+1
#        table = np.zeros(stop-start, dtype=np.intc)
#        for i,l in enumerate(vertices):
#            table[ l-start ] = i + offset
#        cflib.Dofmap_Tabled(self.dm, table.reshape((table.size/stride,stride)), -start )

Dofmap = cflib.Dofmap
Dofmap_Strided = cflib.Dofmap_Strided
Dofmap_Tabled = cflib.new_Dofmap_Tabled
Dofmap_From_Vertices = cflib.Dofmap_From_Vertices

# DEPRECATED
class Dofmap_dep():
    def __init__(self):
        self.vertex_dofs = []
        self.edge_dofs = []
        self.global_dofs = []

        self.R_order = []
    def __init__(self, h, ke, Nnode):
        dofs = cflib.dofArray_frompointer(ke.outp)
        
        nvertdof = 0
        nedgedof = 0
        nglobaldof = 0
        
        for i in xrange(ke.noutp):
            if dofs[i].loc == cflib.LOC_NODE:
                nvertdof += dofs[i].len
            elif dofs[i].loc == cflib.LOC_EDGE:
                nedgedof += dofs[i].len
            else:
                nglobaldof += dofs[i].len
        
        self.vertex_dofs = np.zeros( (Nnode, nvertdof),dtype=np.intc)
        for i in xrange(self.vertex_dofs.size):
            self.vertex_dofs.ravel()[i] = i
            
        self.edge_dofs = np.zeros( (h.hg.n_edge, nedgedof), dtype=np.intc)
        self.global_dofs = np.zeros( (nglobaldof), dtype=np.intc)

# DEPRECATED
def make_dofmap(h, ke, Nnode):
    # Make the dofmap for where the kernel output will be
    # assembled into. It will look like
    # G1 ... E1 E3 ... N1 N2 N3
    # for now...
    # TODO: Do I need to give edges global indices??
    
    v = h.view()
    outlen = 0
    dofs = cflib.dofArray_frompointer(ke.outp)
    
    for i in xrange(ke.noutp):
        if dofs[i].loc == cflib.LOC_NODE:
            outlen+=dofs[i].len*h.hg.l_edge
        else:
            outlen+=dofs[i].len
    
    outmap = np.zeros((v.shape[0], outlen),dtype=np.intc)
    # Loop over the edges
    for i,e in enumerate(v):
        # Loop over the dofs
        off = 0
        for d in xrange(ke.noutp):
            if dofs[d].loc == cflib.LOC_NODE:
                for j,n in enumerate(e):
                    outmap[i,off+dofs[d].len*j:dofs[d].len*(j+1)] = range(n*dofs[d].len,(n+1)*dofs[d].len)
                off+= len(e)*dofs[d].len
            else:
                pass
    return outmap

