import numpy as np
from . import cornflakes_library as cflib


class Hypergraph():

    def __init__(self, alloc_init=1, fname=None):
        self.hg = cflib.hypergraph_t()
        if alloc_init:
            cflib.Hypergraph_Alloc(self.hg, alloc_init)
        if fname is not None:
            with open(fname,'r') as f:
                for l in f:
                    self.Push_Edge([int(_) for _ in l.split() ])
    def __del__(self):
        cflib.Hypergraph_Destroy(self.hg)

    def Push_Edge(self, edge):
        cflib.Hypergraph_Push_Edge(self.hg,np.array(edge,dtype=np.intc))
    #def Get_Edge(self, i):
    #    return hgswig.Hypergraph_Get_Edge_np(self.hg,i)
    def view(self):
        return [
            cflib.Hypergraph_Get_Edge_View_np(self.hg,i)
            for i in range(self.hg.n_types)
        ]
    def __iter__(self):
        for i in range(self.hg.n_types):
            ev = cflib.Hypergraph_Get_Edge_View_np(self.hg,i)
            for e in ev:
                yield e
    def __len__(self):
        # TODO: Optimize this
        return sum([ _.shape[0] for _ in self.view()])

    def Add_Edge_Vertex(self,offset=0):
        # TODO I hate this method. Cornflakes graphs should be immutable
        hgnew = cflib.hypergraph_t()
        cflib.Add_Edge_Vertex(hgnew, self.hg, offset)
        cflib.Hypergraph_Destroy(self.hg)
        self.hg = hgnew

    def Write_To_File(self,fname):
        with open(fname,'w') as f:
            for e in self:
                f.write(" ".join([str(_) for _ in e])+"\n")
