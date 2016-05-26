#ifdef USE_LIS

#include "cfmat_lis.h"
#include <lis.h>

#define data(x) CFMat_LIS_Data(x)
real_t * CFMat_LIS_Place(cfmat_t * self,
			   int n, int * dofs, real_t * ker_out) {
  for(int i=0;i<n;i++) {
      for(int j=0;j<n;j++) {
	lis_matrix_set_value(LIS_ADD_VALUE, dofs[i],dofs[j], ker_out[n*i+j], data(self));
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
}
void CFMat_LIS_Finalize(cfmat_t * self) {
  // TODO
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
}
void CFMat_LIS_New_From_Ptr(cfmat_t * self, int N, void* payload) {
  self->vtable = &CFMat_LIS_vtable;
  self->N = N;
  self->own = 0;
  self->data = payload;
}

#endif
