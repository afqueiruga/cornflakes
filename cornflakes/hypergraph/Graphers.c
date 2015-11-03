#include "Graphers.h"

void Build_Particle_Graph(hypergraph_t * hg, int Npart, int dim, real_t * x, real_t cutoff) {
  spatialhash_t sh;
  Hypergraph_Alloc(hg,2, Npart);
  
  
  Build_New_Hash(&sh, Npart,dim,x, cutoff);

  void action(int a, int b) {
    int v[2] = {a,b};
    printf("%d,%d\n",a,b);
    Hypergraph_Push_Edge(hg,v);
  }
  SpatialHash_Scanall(&sh,x,action);
  
  SpatialHash_destroy(&sh);
}
