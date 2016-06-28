#include "sparsity_pattern.h"

#include <stdlib.h>

void Sparsity_Init(sparsity_t * self, int N, int nrow_guess) {
  self->nnz = 0;
  self->N = N;
  self->Ibuild = malloc(N*sizeof(indexset_t*));
  for(int i=0; i<N; i++) {
    self->Ibuild[i] = malloc(sizeof(indexset_t));
    IndexSet_New(self->Ibuild[i], nrow_guess);
  }
}
void Sparsity_Add_NNZ(sparsity_t * self, index_t i, index_t j) {
  if(IndexSet_Insert(self->Ibuild[i],j)) {
    self->nnz++;
  }
}
void Sparsity_Fill(sparsity_t * self,
		   int onum, kernel_t * ke, hypergraph_t * hg, dofmap_t ** dms)
{
  int i,j,k,m, hex,hx;
  hyperedges_t * he;
  /* Loop over the graph sets */
  for(he = hg->he; he < hg->he+hg->n_types ; he++) {
    /* Loop over the edges */
    hypervertex_t * edge;
    for(hex=0; hex<he->n_edge; hex++) {
      edge = Hyperedges_Get_Edge(he, hex);
      
      /* Collect the dofs for this output*/
      int nselect, dim;
      hypervertex_t select[he->l_edge];
      int nalldofs = kernel_outp_ndof(ke,ke->outp+onum, he->l_edge);
      int alldofs[nalldofs];
      int iter=0;
      for(m=0; m<ke->outp[onum].nmap; m++) {
	dofmap_t * dmap;
	int maxlen;
	int mnum = ke->outp[onum].map_nums[m];
	//printf("t %d mnum %d\n",t,mnum);
	k_map_t kmap = ke->maps[ mnum ];
	kmap(edge,he->l_edge, select,&nselect, &dim);
      
	dmap = dms[mnum];
	maxlen = Dofmap_Max_Len(dmap);
	int dofs[maxlen], ndof;
	hypervertex_t V;
	for(j=0; j<nselect; j++) { 
	  V = select[j];
	  Dofmap_Get(dmap, V, dofs,&ndof);
	  for(k=0;k<ndof;k++) {
	    alldofs[iter+k] = dofs[k];
	  }
	  iter+=ndof;
	}
      }
      
      /* Push them into the pattern */
      for(i=0; i<nalldofs; i++ ) {
	for(j=0; j<nalldofs; j++) {
	  Sparsity_Add_NNZ(self, i,j);
	}
      }
    } //  end hex loop
  } // end he loop
}

void Sparsity_Destroy(sparsity_t * self) {
  for(int i=0; i<self->N; i++) {
    IndexSet_Destroy(self->Ibuild[i]);
    free(self->Ibuild[i]);
  }
  free(self->Ibuild);
}
