#include "hypergraph.h"
//#include "malloc.h"
#include "stdlib.h"

void Hypergraph_Alloc(hypergraph_t * hg, int n_alloc) {
  hg->n_types = 0;
  hg->n_alloc = n_alloc;
  hg->he = (hyperedges_t*)malloc(sizeof(hyperedges_t)*n_alloc );
}

void Hypergraph_Destroy(hypergraph_t * hg) {
  int i;
  for(i=0;i<hg->n_types;i++) {
    Hyperedges_Destroy(hg->he+i);
  }
  free(hg->he);
}
void Hypergraph_Push_Edge(hypergraph_t * hg, int l_edge, hypervertex_t* verts) {
  int i;
  // See if one exists
  for(i=0; i<hg->n_types; i++) {
    if(hg->he[i].l_edge == l_edge) {
      Hyperedges_Push_Edge(hg->he+i, verts);
      return;
    }
  }
  // Make a new collection for this type and append it
  // Need more room?
  if(hg->n_types == hg->n_alloc) {
    hg->he = realloc( hg->he, 2*hg->n_alloc*sizeof(hyperedges_t) );
    hg->n_alloc *= 2;
    
  }
  // Allocate and push
  Hyperedges_Alloc(hg->he+hg->n_types, l_edge, 10);
  Hyperedges_Push_Edge(hg->he+hg->n_types, verts);
  hg->n_types++;
}

/*
 *  Allocate a new hyperedge collection
 */
void Hyperedges_Alloc(hyperedges_t * he, int l_edge, int alloc_init) {
  he->n_edge = 0;
  he->l_edge = l_edge;
  he->n_e_alloc = alloc_init;
  he->edges = (int*)malloc( sizeof(int)*l_edge* alloc_init );
}
/*
 * Free it.
 */
void Hyperedges_Destroy(hyperedges_t * he) {
  free(he->edges);
}

hypervertex_t * Hyperedges_Get_Edge(hyperedges_t * he, int i) {
  if(i<0 || i>=he->n_edge) {
    return NULL;
  }
  return he->edges + he->l_edge*i;
}


void Hyperedges_Push_Edge(hyperedges_t * he, hypervertex_t * verts) {
  int i;
  // Do some automatic memory allocation
  if(he->n_edge == he->n_e_alloc) {
    he->edges=realloc(he->edges,2*he->n_e_alloc*he->l_edge*sizeof(int) );
    he->n_e_alloc *= 2;
  }
  // Append the edge
  for(i=0;i<he->l_edge;i++) {
    he->edges[he->n_edge*he->l_edge + i] = verts[i];
  }
  he->n_edge += 1;
}


