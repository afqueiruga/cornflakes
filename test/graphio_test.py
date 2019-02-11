import unittest
import os
import cornflakes as cf

abspath = os.path.dirname(os.path.abspath(__file__))

class TestLoadFiles(unittest.TestCase):
    def test_load_gmsh(self):
        H = cf.GraphIO.Load_gmsh(abspath+"/good_gmsh_v2.msh",gdim=3)
    def test_load_gmsh_wrong_version(self):
        with self.assertRaises(RuntimeError):
            cf.GraphIO.Load_gmsh(abspath+"/good_gmsh_v4.msh",gdim=3)
    def test_load_gmsh_no_file(self):
        with self.assertRaises(Exception):
            cf.GraphIO.Load_gmsh(abspath+"/not_a_file.msh",gdim=2)
        
if __name__=="__main__":
    unittest.main()
