#include "Graphers.h"

#include "math.h"
double dist(int dim, real_t * x, real_t * y) {
  int i;
  double ac = 0.0;
  for(i=0;i<dim;i++) ac += (x[i]-y[i])*(x[i]-y[i]);
  return sqrt(ac);
}

void Build_Particle_Graph(hypergraph_t * hg, int Npart, int dim, real_t * x, real_t cutoff) {
  spatialhash_t sh;
  Hypergraph_Alloc(hg,2, Npart);
  
  Build_New_Hash(&sh, Npart,dim,x, cutoff);

  void action(int a, int b) {
    if(dist(dim, x+dim*a,x+dim*b)<=cutoff) {
      int v[2] = {a,b};
      Hypergraph_Push_Edge(hg,v);
    }
  }
  SpatialHash_Scanall(&sh,x,action);
  
  SpatialHash_destroy(&sh);
}
