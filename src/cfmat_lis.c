#ifdef USE_LIS

#include "cfmat_lis.h"
#include <lis.h>

#define data(x) CFMat_LIS_Data(x)
real_t * CFMat_LIS_Place(cfmat_t * self,
			   int ln, int * ldofs, int rn, int * rdofs, real_t * ker_out) {
  for(int i=0;i<ln;i++) {
    for(int j=0;j<rn;j++) {
      lis_matrix_set_value(LIS_ADD_VALUE, ldofs[i],rdofs[j], ker_out[rn*i+j], data(self));
    }
  }
  return ker_out + n*n;
}
void CFMat_LIS_Destroy(cfmat_t * self) {
  if(self->own) {
    lis_matrix_unset(data(self));
    lis_matrix_destroy(data(self));
  }
}
void CFMat_LIS_Wipe(cfmat_t * self) {
  // TODO
  // I Don't think LIS Does this???
  // What about LIS Matrix Scale?
}
void CFMat_LIS_Finalize(cfmat_t * self) {
  lis_matrix_assemble(data(self));
}

const _CFMAT_VTABLE_t CFMat_LIS_vtable = {
  .Place = &CFMat_LIS_Place,
  .Destroy = &CFMat_LIS_Destroy,
  .Wipe = &CFMat_LIS_Wipe,
  .Finalize = &CFMat_LIS_Finalize
};
void CFMat_LIS_New(cfmat_t * self, int N) {
  self->vtable = &CFMat_LIS_vtable;
  self->N = N;
  self->own = 1;
  // TODO
  lis_matrix_create(/*MPI_COMM_WORLD*/0, (LIS_MATRIX*)&self->data);
  lis_matrix_set_size(data(self), 0,N);
  lis_matrix_set_type(data(self), LIS_MATRIX_CSR);
}
void CFMat_LIS_New_From_Ptr(cfmat_t * self, int N, void* payload) {
  self->vtable = &CFMat_LIS_vtable;
  self->N = N;
  self->own = 0;
  self->data = payload;
}

#endif
