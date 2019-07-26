import unittest
import os
import cornflakes as cf

dm = cf.Dofmap_Strided(3,20)
assert dm.Max_Len() == 3, "wrong max"
assert dm.U.strided.stride == 3, "wrong stride"
assert dm.U.strided.offset == 20, "wrong offset"

n = dm.Get_np(20)
assert n.size == 3, "wrong array size"
assert n[0] == 80, "wrong calculated value!"
assert n[1] == 81, "wrong calculated value!"
assert n[2] == 82, "wrong calculated value!"
