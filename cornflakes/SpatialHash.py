import cornflakes_library as cflib
import numpy as np

class SpatialHash():

    def __init__(self, X):
        self.sh = cflib.spatialhash_t()
        cflib.Build_New_Hash(self.sh, X.shape[0], X.shape[1], X)

    def __del__(self):
        cflib.SpatialHash_destroy(sh)
