#include "cfmat_default.h"

#include <stdlib.h>

#define data(x) CFMat_Default_Data(x)
real_t * CFMat_Default_Place(cfmat_t * self,
			     int ln, int * ldofs,int rn, int * rdofs, real_t * ker_out) {
  int i,j;

    // Fill this block
    for(i=0;i<ln;i++) {
      for(j=0;j<rn;j++) {
	data(self)->IIiter[rn*i + j ] = ldofs[i];
	data(self)->JJiter[rn*i + j ] = rdofs[j];
	data(self)->Viter [rn*i + j ] = ker_out[ rn*i + j];
      }
    }
    // Advance our iterators.
    data(self)->IIiter += rn*ln;
    data(self)->JJiter += rn*ln;
    data(self)->Viter += rn*ln;
    return ker_out + rn*ln;

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
  for(i=0;i<self->N;i++) {
    data(self)->V[i]=0.0;
  }
  data(self)->Viter = data(self)->V;
  data(self)->IIiter = data(self)->II;
  data(self)->JJiter = data(self)->JJ;
}
void CFMat_Default_Finalize(cfmat_t * self) {
  //pass
}

const _CFMAT_VTABLE_t CFMat_Default_vtable = {
  .Place = &CFMat_Default_Place,
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
