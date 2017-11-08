#include "filter.h"

#include <stdio.h>

// TODO: Update this to new call signature
void filter(kernel_t * ke, hypergraph_t * hg,
	    dofmap_t ** dofmaps, cfdata_t ** data,
	    hypergraph_t *htrue, hypergraph_t *hfalse)
{
  int i,j, hex,hx;
  hyperedges_t * he;
  
  if(htrue)  Hypergraph_Alloc(htrue, 1);
  if(hfalse) Hypergraph_Alloc(hfalse,1);
  /* Loop over the graph sets */
  for(he = hg->he; he < hg->he+hg->n_types ; he++) {
    /* Allocate the local vectors for this size edge */
    int len_ker_in = kernel_inps_len(ke, he->l_edge);
    int len_ker_out = kernel_outps_len(ke, he->l_edge);
    if(len_ker_out != 1) printf("This is not valid kernel for fitlering!\n");
    
    real_t ker_in[ len_ker_in];
    real_t ker_out[len_ker_out];
    
    /* Loop over the edges */
    hypervertex_t * edge;
    for(hex=0; hex<he->n_edge; hex++) {
      edge = Hyperedges_Get_Edge(he, hex);
      /* Collect the data */
      // collect(ker_in, ke, edge,he->l_edge, dofmaps,data); // TODO: Optimize by moving some overheard outside of loop
      /* Calculate the kernel */
      for(i=0;i<len_ker_out;i++) ker_out[i] = 0.0;
      ke->eval(he->l_edge, ker_in, ker_out);
      /* Add to the graph */
      if(ker_out[0]) {
	if(htrue)  Hypergraph_Push_Edge(htrue,he->l_edge,edge);
      } else {
	if(hfalse) Hypergraph_Push_Edge(hfalse,he->l_edge,edge);
      }
    }
  }

  
}

