from cornflakes import Hypergraph, cflib

def graph_test_A():
    H = Hypergraph(1)
    assert H.hg.n_alloc==1, "Allocation messed up"
    assert H.hg.n_types==0, "hg should be empty"
    
    H.Push_Edge([1,2,3])
    assert H.hg.n_types==1, "hg should contain a new type"
    print H.view()
    H.Push_Edge([3,4,5])
    assert H.hg.n_types==1, "hg should have the same number of tyes"
    print H.view()
    H.Push_Edge([3,4])
    assert H.hg.n_types==2, "hg should contain a new type"
    print H.view()
    H.Push_Edge([4,5])
    H.Push_Edge([1,2,3,4])
    print H.view()
    assert H.hg.n_types==3, "Wrong number of types"
    assert H.hg.n_alloc>=H.hg.n_types, "Insufficient memory allocated"
    
    v = H.view()
    v[0][0,0] = 10
    assert H.view()[0][0,0]==10, "Views should be mutable"
    print v
    print "Graph Test A Passed"

    
if __name__=="__main__":
    graph_test_A()
    print "All graph tests passed"
