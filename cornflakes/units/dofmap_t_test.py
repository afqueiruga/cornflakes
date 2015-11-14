from cornflakes import cflib

dm = cflib.dofmap_t()
cflib.Dofmap_Strided(dm,3,20)
assert cflib.Dofmap_Max_Len(dm) == 3, "wrong max"
assert dm.U.strided.stride == 3, "wrong stride"
assert dm.U.strided.offset == 20, "wrong offset"

n = cflib.Dofmap_Get_np(dm, 20)
assert n.size == 3, "wrong array size"
assert n[0] == 80, "wrong calculated value!"
assert n[1] == 81, "wrong calculated value!"
assert n[2] == 82, "wrong calculated value!"
