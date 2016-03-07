from cornflakes import SpatialHash, cflib, ParticlePlacers as PP

def hash_test_A():
    X = init_grid(10,10, 0.0,(1.0,0.0),(0.0,1.0))
    SH = SpatialHash(X)

    
    
if __name__=="__main__":
    hash_test_A()
    print "All hash tests passed"
