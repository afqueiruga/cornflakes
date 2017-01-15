#
# Wrappers and Python routines for the C routines found in graphers.c
#

import cornflakes_library as cflib
from Hypergraph import Hypergraph
import numpy as np

def Build_Pair_Graph(x, cutoff, y=None):
    H = Hypergraph(0)
    if y is not None:
        cflib.Build_Pair_Graph_2Sets( H.hg, x, y, cutoff)
    else:
        cflib.Build_Pair_Graph( H.hg, x, cutoff)
    return H
def Build_Proximity_Graph_Uniform(x, cutoff, y=None):
    H = Hypergraph(0)
    if y is not None:
        cflib.Build_Proximity_Graph_2Sets_Uniform(H.hg, x, y, cutoff)
    else:
        cflib.Build_Proximity_Graph_Uniform(H.hg, x, cutoff)
    return H

def Build_Proximity_Graph_Variable(x, r, y=None):
    H = Hypergraph(0)
    if y is not None:
        cflib.Build_Proximity_Graph_2Sets_Variable_np(H.hg, x,y,r)
    else:
        cflib.Build_Proximity_Graph_Variable_np(H.hg, x,r)
    return H

def Build_Proximity_Graph_Given_Length(x, N_desired, cutoff, y=None):
    H = Hypergraph(0)
    if y is not None:
        r = np.empty( (y.shape[0],), dtype=x.dtype)
        cflib.Build_Proximity_Graph_2Sets_Given_Length_np( H.hg, x, y, N_desired,cutoff,r)
    else:
        r = np.empty( (x.shape[0],), dtype=x.dtype)
        cflib.Build_Proximity_Graph_Given_Length_np( H.hg, x, N_desired, cutoff, r)
    return H,r

def Tie_Cells_And_Particles(*args):
    raise NotImplementedError()

