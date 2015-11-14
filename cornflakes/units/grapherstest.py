from cornflakes import Hypergraph, cflib, write_graph, PP


H = Hypergraph(5)
X = PP.init_grid(10,10,[0,0],[1,0],[0,1],0.0)
cflib.Build_Adjacency_Graph_Uniform(H.hg, X, 0.15)
print H.view()
