#include "hypergraph.h"
#include "malloc.h"

/*
  Allocate a new hypergraph
*/
void Hypergraph_Alloc(hypergraph_t * hg, int l_edge, int alloc_init) {
  hg->n_edge = 0;
  hg->l_edge = l_edge;
  hg->n_e_alloc = alloc_init;
  hg->edges = (int*)malloc( sizeof(int)*l_edge* alloc_init );
}
/*
 * Free it.
 */
void Hypergraph_Destroy(hypergraph_t * hg) {
  free(hg->edges);
}

int * Hypergraph_Get_Edge(hypergraph_t * hg, int i) {
  if(i<0 || i>=hg->n_edge) {
    return NULL;
  }
  return hg->edges + hg->l_edge*i;
}


void Hypergraph_Push_Edge(hypergraph_t * hg, int * verts) {
  int i;
  // Do some automatic memory allocation
  if(hg->n_edge == hg->n_e_alloc) {
    hg->edges=realloc(hg->edges,2*hg->n_e_alloc*hg->l_edge*sizeof(int) );
    hg->n_e_alloc *= 2;
  }
  // Append the edge
  for(i=0;i<hg->l_edge;i++) {
    hg->edges[hg->n_edge*hg->l_edge + i] = verts[i];
  }
  hg->n_edge += 1;
}


