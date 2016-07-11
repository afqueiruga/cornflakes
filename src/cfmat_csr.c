#include "cfmat_csr.h"

#include <stdlib.h>
#include <stdlib.h>

#include "cfmat_default.h"

#define data(x) CFMat_CSR_Data(x)
void CFMat_CSR_Finalize_Sparsity(cfmat_t * self) {
  Sparsity_Make_CSR(&self->sparse,&data(self)->IA, &data(self)->JA);
  data(self)->nnz = self->sparse.nnz;
  data(self)->V = calloc(data(self)->nnz, sizeof(real_t));
  Sparsity_Destroy(&self->sparse);
}
void CFMat_CSR_Set_Value(cfmat_t * self,
			     int i, int j, real_t v) {
  /* Contract: i,j is in the sparsity structure. */
  int a = data(self)->IA[i], b = data(self)->IA[i+1];
  if(data(self)->JA[a] == j) {
    data(self)->V[a] += v;
  } else {
    do {
      if( j > data(self)->JA[(a+b)/2] ) {
	a = (a+b)/2;
      } else {
	b = (a+b)/2;
      }
    } while( b-a > 1 );
    data(self)->V[b] += v;
  }
}
real_t * CFMat_CSR_Place(cfmat_t * self,
			      int n, int * dofs, real_t * ker_out) {
  int i,j;
  for(i=0;i<n;i++) {
    for(j=0;j<n;j++) {
      CFMat_CSR_Set_Value(self, dofs[i],dofs[j], ker_out[n*i+j]);
    }
  }
  return ker_out + n*n;
}
void CFMat_CSR_Wipe(cfmat_t * self) {
  int i;
  for(i=0;i<data(self)->nnz;i++) {
    data(self)->V[i]=0.0;
  }
}
void CFMat_CSR_Finalize(cfmat_t * self) {
  //pass
}

void CFMat_CSR_Destroy(cfmat_t * self) {
  if(self->own) {
    if(data(self)->IA) free(data(self)->IA);
    if(data(self)->JA) free(data(self)->JA);
    if(data(self)->V) free(data(self)->V);
  }
  free(data(self));
}


const _CFMAT_VTABLE_t CFMat_CSR_vtable = {
  .Add_Sparsity = CFMat_Default_Add_Sparsity,
  .Finalize_Sparsity = CFMat_CSR_Finalize_Sparsity,
  .Place = CFMat_CSR_Place,
  .Set_Value = CFMat_CSR_Set_Value,
  .Destroy = CFMat_CSR_Destroy,
  .Wipe = CFMat_CSR_Wipe,
  .Finalize = CFMat_CSR_Finalize
};


void CFMat_CSR_New(cfmat_t * self, int N) {
  self->vtable = &CFMat_CSR_vtable;
  self->N = N;
  self->own = 1;
  self->data = malloc(sizeof(CFMat_CSR_data_t));
  Sparsity_Init(&self->sparse,self->N,20);

  data(self)->IA = NULL;
  data(self)->nnz = 0;
  data(self)->JA = NULL;
  data(self)->V = NULL;
}
