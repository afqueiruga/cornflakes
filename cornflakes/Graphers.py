#
# Wrappers and Python routines for the C routines found in graphers.c
#

import cornflakes_library as cflib
from Hypergraph import Hypergraph
import numpy as np

def Build_Pair_Graph(x, cutoff):
    H = Hypergraph(0)
    cflib.Build_Pair_Graph( H.hg, x, cutoff)
    return H
def Build_Proximity_Graph_Uniform(x, cutoff):
    H = Hypergraph(0)
    cflib.Build_Proximity_Graph_Uniform(H.hg, x, cutoff)
    return H

def Build_Proximity_Graph_Variable(x, r):
    H = Hypergraph(0)
    cflib.Build_Proximity_Graph_Variable_np(H.hg, x,r)
    return H

def Build_Proximity_Graph_Given_Length(x, N_desired, cutoff):
    H = Hypergraph(0)
    #print x.shape
    r = np.empty( (x.shape[0],), dtype=x.dtype)
    #print r.shape
    cflib.Build_Proximity_Graph_Given_Length_np( H.hg, x, N_desired, cutoff, r)
    return H,r

def Tie_Cells_And_Particles(*args):
    raise NotImplementedError()
