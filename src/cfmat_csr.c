#include "cfmat_csr.h"

#include <stdlib.h>
#include <stdlib.h>

#define data(x) CFMat_CSR_Data(x)
real_t * CFMat_CSR_Place(cfmat_t * self,
			      int n, int * dofs, real_t * ker_out) {
  int i,j;

    // Fill this block
    for(i=0;i<n;i++) {
      for(j=0;j<n;j++) {
      }
    }
    
    return ker_out + n*n;
}
void CFMat_CSR_Set_Value(cfmat_t * self,
			     int i, int j, real_t v) {

}
void CFMat_CSR_Destroy(cfmat_t * self) {
  if(self->own) {
    free(data(self)->I);
    if(data(self)->J) free(data(self)->J);
    if(data(self)->V) free(data(self)->V);
  }
  free(data(self));
}
void CFMat_CSR_Wipe(cfmat_t * self) {
  int i;
  for(i=0;i<data(self)->nnz;i++) {
    /* data(self)->V[i]=0.0; */
  }
}
void CFMat_CSR_Finalize(cfmat_t * self) {
  //pass
}


const _CFMAT_VTABLE_t CFMat_CSR_vtable = {
  .Place = &CFMat_CSR_Place,
  .Set_Value = &CFMat_CSR_Set_Value,
  .Destroy = &CFMat_CSR_Destroy,
  .Wipe = &CFMat_CSR_Wipe,
  .Finalize = &CFMat_CSR_Finalize
};


void CFMat_CSR_New(cfmat_t * self, int N) {
  self->vtable = &CFMat_CSR_vtable;
  self->N = N;
  self->own = 1;
  self->data = malloc(sizeof(CFMat_CSR_data_t));
  data(self)-> I = malloc(sizeof(int)*N);
  data(self)->nnz = 0;
  data(self)->J = NULL;
  data(self)->V = NULL;
}
void CFMat_CSR_Add_Sparsity(cfmat_t * self, int onum,
			    kernel_t * ke, hypergraph_t * hg, dmap_t ** dms) {
  
}
void CFMat_CSR_Set_Sparsity(cfmat_t * self) {
  data(self)->J = malloc(data(self)->nnz * sizeof(int) );
  data(self)->V = malloc(data(self)->nnz * sizeof(real_t) );
}
