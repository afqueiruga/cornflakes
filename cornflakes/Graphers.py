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
    cflib.Build_Proximity_Graph_Uniform(H.hg, X, cutoff)
    return H

def Build_Proximity_Graph_Variable(X, r):
    H = Hypergraph(0)
    cflib.Build_Proximity_Graph_Variable_np(H.hg, X,r)
    return H

def Tie_Cells_And_Particles(*args):
    raise NotImplementedError()
