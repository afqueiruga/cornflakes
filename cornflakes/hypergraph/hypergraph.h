#ifndef __H_HYPERGRAPH__
#define __H_HYPERGRAPH__

typedef struct hypergraph_t {
  //int n_vertex; //How many vertices are referenced? (Or just the large number+1) The graph don't care
  int n_edge; //Number of edges
  int n_e_alloc; //How much memory has been allocated
  int l_edge; //How hyper is the edge? Number of verts per edge
  int * edges; //Pointer to the edges
} hypergraph_t;

/*
  Allocate a new hypergraph
*/
void Hypergraph_Alloc(hypergraph_t * hg, int l_edge, int alloc_init);
void Hypergraph_Push_Edge(hypergraph_t * hg, int * verts);
int * Hypergraph_Get_Edge(hypergraph_t * hg, int i);
void Hypergraph_Destroy(hypergraph_t * hg);
#endif
