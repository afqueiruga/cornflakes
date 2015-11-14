#include "graphers.h"

#include "math.h"
double dist(int dim, real_t * x, real_t * y) {
  int i;
  double ac = 0.0;
  for(i=0;i<dim;i++) ac += (x[i]-y[i])*(x[i]-y[i]);
  return sqrt(ac);
}

void Build_Pair_Graph(hypergraph_t * hg, int Npart, int dim, real_t * x, real_t cutoff) {
  spatialhash_t sh;
  //printf("alloc hg\n");
  Hypergraph_Alloc(hg,1); //2, Npart);
  //printf("alloc hash\n");
  Build_New_Hash(&sh, Npart,dim,x, cutoff);
  //printf("scan\n");
  void action(int a, int b) {
    if(dist(dim, x+dim*a,x+dim*b)<=cutoff) {
      int v[2] = {a,b};
      Hypergraph_Push_Edge(hg,2,v);
    }
  }
  SpatialHash_Scanall(&sh,x,action);
  //printf("destroy hash\n");
  SpatialHash_destroy(&sh);
}

void Build_Adjacency_Graph_Uniform(hypergraph_t * hg, int Npart, int dim, real_t * x, real_t cutoff) {
  spatialhash_t sh;
  int A;
  //printf("alloc hg\n");
  Hypergraph_Alloc(hg,1); //2, Npart);
  //printf("alloc hash\n");
  Build_New_Hash(&sh, Npart,dim,x, cutoff);
  //printf("scan\n");
  int list[30]; // TODO: Auto allocate this list
  int Nlist=0;
  void action(int FOO, int b) {
    if(dist(dim, x+dim*A,x+dim*b)<=cutoff) {
      list[Nlist] = b;
      Nlist++;
    }
  }

  for(A=0; A<Npart; A++) {
    list[0] = A;
    Nlist = 1;
    SpatialHash_ScanPt(&sh, x+dim*A, action);
    Hypergraph_Push_Edge(hg,Nlist,list);
  }
  //SpatialHash_Scanall(&sh,x,action);
  //printf("destroy hash\n");
  SpatialHash_destroy(&sh);
}
