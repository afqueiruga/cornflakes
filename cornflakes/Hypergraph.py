import cornflakes_library as cflib
import numpy as np

class Hypergraph():

    def __init__(self, alloc_init=1):
        self.hg = cflib.hypergraph_t()
        cflib.Hypergraph_Alloc(self.hg, alloc_init)

    def __del__(self):
        cflib.Hypergraph_Destroy(self.hg)

    def Push_Edge(self, edge):
        cflib.Hypergraph_Push_Edge(self.hg,np.array(edge,dtype=np.intc))
    #def Get_Edge(self, i):
    #    return hgswig.Hypergraph_Get_Edge_np(self.hg,i)
    def view(self):
        return [
            cflib.Hypergraph_Get_Edge_View_np(self.hg,i)
            for i in xrange(self.hg.n_types)
        ]
