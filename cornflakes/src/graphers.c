#include "graphers.h"


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
  int list[30]; // TODO: Auto allocate this list BADBADBAD
  int Nlist=0;
  void action(int FOO, int b) {
    if( b!=A &&  dist(dim, x+dim*A,x+dim*b)<=cutoff) {
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


void Add_Edge_Vertex(hypergraph_t * hgnew, hypergraph_t * hgold, int offset) {
  int i,j,k;
  Hypergraph_Alloc(hgnew, hgold->n_types);

  int edgenum = offset;
  for(i=0;i<hgold->n_types;i++) {
    hyperedges_t * he = hgold->he + i;
    hypervertex_t edge[ he->l_edge + 1];
    hypervertex_t * e;
    for(j=0;j<he->n_edge; j++) {
      e = Hyperedges_Get_Edge(he, j);
      for(k=0;k<he->l_edge;k++) {
	edge[k] = e[k];
      }
      edge[he->l_edge] = edgenum;
      Hypergraph_Push_Edge(hgnew,he->l_edge+1,edge);
      edgenum++;
    }
  }
}
