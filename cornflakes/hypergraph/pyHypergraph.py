import mylibrary as hgswig
import numpy as np

class Hypergraph():

    def __init__(self, l_edge=None, alloc_init=10):
        self.hg = hgswig.hypergraph_t()
        if l_edge:
            hgswig.Hypergraph_Alloc(self.hg, l_edge, alloc_init)

    def __del__(self):
        hgswig.Hypergraph_Destroy(self.hg)

    def Push_Edge(self, edge):
        hgswig.Hypergraph_Push_Edge_np(self.hg,np.array(edge,dtype=np.intc))
    def Get_Edge(self, i):
        return hgswig.Hypergraph_Get_Edge_np(self.hg,i)
    def view(self):
        return hgswig.Hypergraph_Get_View_np(self.hg)
        
