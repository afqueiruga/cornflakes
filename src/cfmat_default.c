#include "cfmat_default.h"

#include <stdlib.h>

#define data(x) CFMat_Default_Data(x)
real_t * CFMat_Default_Place(cfmat_t * self,
			      int n, int * dofs, real_t * ker_out) {
  int i,j;

    // Fill this block
    for(i=0;i<n;i++) {
      for(j=0;j<n;j++) {
	data(self)->IIiter[n*i + j ] = dofs[i];
	data(self)->JJiter[n*i + j ] = dofs[j];
	data(self)->Viter [n*i + j ] = ker_out[ n*i + j];
      }
    }
    // Advance our iterators.
    data(self)->IIiter += n*n;
    data(self)->JJiter += n*n;
    data(self)->Viter += n*n;
    return ker_out + n*n;

}
void CFMat_Default_Set_Value(cfmat_t * self,
			     int i, int j, real_t v) {
  *data(self)->IIiter = i;
  *data(self)->JJiter = j;
  *data(self)->Viter = v;
  data(self)->IIiter += 1;
  data(self)->JJiter += 1;
  data(self)->Viter += 1;
}
void CFMat_Default_Destroy(cfmat_t * self) {
  if(self->own) {
    free(data(self)->V);
    free(data(self)->II);
    free(data(self)->JJ);
  }
  free(data(self));
}
void CFMat_Default_Wipe(cfmat_t * self) {
  int i;
  /* for(i=0;i<self->N;i++) { */
    /* data(self)->V[i]=0.0; */
  /* } */
  data(self)->Viter = data(self)->V;
  data(self)->IIiter = data(self)->II;
  data(self)->JJiter = data(self)->JJ;
}
void CFMat_Default_Finalize(cfmat_t * self) {
  //pass
}

const _CFMAT_VTABLE_t CFMat_Default_vtable = {
  .Place = &CFMat_Default_Place,
  .Set_Value = &CFMat_Default_Set_Value,
  .Destroy = &CFMat_Default_Destroy,
  .Wipe = &CFMat_Default_Wipe,
  .Finalize = &CFMat_Default_Finalize
};

void CFMat_Default_New(cfmat_t * self, int onum,
			kernel_t * ke, hypergraph_t * hg, int ndof) {
  int matsize, oplen;
  int i;

  self->vtable = &CFMat_Default_vtable;
  self->own = 1;
  self->N = ndof;
  
  self->data = malloc(sizeof(struct CFMat_Default_data_t));
  
  matsize = 0;
  for( i=0; i<hg->n_types; i++ ) {
    oplen = kernel_outp_len(ke,ke->outp+onum,hg->he[i].l_edge);
    matsize += hg->he[i].n_edge * oplen;
  }
  data(self)->V = malloc( matsize * sizeof(real_t) );
  data(self)->II = malloc( matsize * sizeof(int) );
  data(self)->JJ = malloc( matsize * sizeof(int) );
  data(self)->nalloc = matsize;
  // For good measure:
  CFMat_Default_Wipe(self); 
}
void CFMat_Default_From_Array(cfmat_t * self, int ndof,
			       real_t * V, int * II, int * JJ) {
  self->vtable = &CFMat_Default_vtable;
  self->data = malloc(sizeof(struct CFMat_Default_data_t));

  self->N = ndof;

  self->own = 0;
  
  data(self)->V = V;
  data(self)->II = II;
  data(self)->JJ = JJ;
  
  CFMat_Default_Wipe(self);
}
