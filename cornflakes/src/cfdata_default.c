#include "cfdata_default.h"

#include <stdlib.h>

#define data(x) CFData_Default_Data(x)
/* Member methods */
real_t * CFData_Default_Place(cfdata_t * self,
			      int n, int * dofs, real_t * ker_out) {
  int i;
  for(i=0;i<n;i++) {
    data(self)[dofs[i]] += ker_out[i];
  }
  return ker_out + n;
}
void CFData_Default_Destroy(cfdata_t * self) {
  if(self->own) {
    free(data(self));
  }
}
void CFData_Default_Wipe(cfdata_t * self) {
  int i;
  for(i=0;i<self->N;i++) {
    data(self)[i]=0.0;
  }
}
void CFData_Default_Finalize(cfdata_t * self) {
  //pass
}
void CFData_Default_Get_Values(cfdata_t * self, int ndof,int *dofs, real_t * vals) {
  int k;
  for(k=0;k<ndof;k++) vals[k] = ((real_t*)self->data)[ dofs[k] ];
}

/* vtable */
const _CFDATA_VTABLE_t cfdata_default_vtable = {
  .Get_Values = CFData_Default_Get_Values,
  .Place = &CFData_Default_Place,
  .Destroy = &CFData_Default_Destroy,
  .Wipe = &CFData_Default_Wipe,
  .Finalize = &CFData_Default_Finalize
};

/* Constructors */
void CFData_Default_New_From_Ptr(cfdata_t * self, int N, real_t * payload) {
  self->vtable = &cfdata_default_vtable;
  self->data = payload;
  self->own = 0;
  self->N = N;
}
void CFData_Default_New(cfdata_t * self, int N) {
  self->vtable = &cfdata_default_vtable;
  self->data = malloc( N*sizeof(real_t) );
  self->own = 1;
  self->N = N;
}

