#ifdef USE_LIS

#include "cfmat_lis_2.h"
#include "cfmat_default.h"
#include <lis.h>
#include <stdlib.h>

#define data(x) CFMat_preLIS_Overhead(x)
real_t * CFMat_preLIS_Place(cfmat_t * self,
			   int n, int * dofs, real_t * ker_out) {
  return CFMat_Place(&data(self)->triplet, n,dofs,ker_out);
}
void CFMat_preLIS_Set_Value(cfmat_t * self,int i, int j, real_t v) {
  CFMat_Set_Value(&data(self)->triplet, i,j,v);
  /* lis_matrix_set_value(LIS_ADD_VALUE, i,j,v, data(self)); */
}
void CFMat_preLIS_Destroy(cfmat_t * self) {
  if(self->own) {
    lis_matrix_unset(data(self)->K);
    lis_matrix_destroy(data(self)->K);
  }
  CFMat_Destroy(&data(self)->triplet);
  free(data(self));
}
void CFMat_preLIS_Wipe(cfmat_t * self) {
  CFMat_Wipe(&data(self)->triplet);
}
void CFMat_preLIS_Finalize(cfmat_t * self) {
  /* Now place in all of the triplet data into K */
  CFMat_Default_data_t* tripdat = CFMat_Default_Data(&data(self)->triplet);
  int nnz = tripdat->IIiter - tripdat->II;
  lis_matrix_set_coo(nnz, tripdat->II,tripdat->JJ,tripdat->V, data(self)->K);

  lis_matrix_assemble(data(self)->K);
}

const _CFMAT_VTABLE_t CFMat_preLIS_vtable = {
  .Place = &CFMat_preLIS_Place,
  .Set_Value = CFMat_preLIS_Set_Value,
  .Destroy = &CFMat_preLIS_Destroy,
  .Wipe = &CFMat_preLIS_Wipe,
  .Finalize = &CFMat_preLIS_Finalize
};
void CFMat_preLIS_New(cfmat_t * self, int N, int onum,
		      kernel_t * ke, hypergraph_t * hg, int ndof ) {
  self->vtable = &CFMat_preLIS_vtable;
  self->N = N;
  self->own = 1;
  self->data = malloc(sizeof(CFMat_preLIS_data_t));
  CFMat_Default_New(&data(self)->triplet, onum,ke,hg,ndof);
  lis_matrix_create(/*MPI_COMM_WORLD*/0, (LIS_MATRIX*)&data(self)->K);

  lis_matrix_set_size(data(self)->K, 0,N);
  lis_matrix_set_type(data(self)->K, LIS_MATRIX_COO);
  /* HOW TO PREALLOC */
}
void CFMat_preLIS_New_From_Ptr(cfmat_t * self, int N, void* payload) {
  self->vtable = &CFMat_preLIS_vtable;
  self->N = N;
  self->own = 0;
  self->data = malloc(sizeof(CFMat_preLIS_data_t));
  // TODO ERROR:
  CFMat_Default_New(&data(self)->triplet, 0, NULL,NULL,N);

  data(self)->K = payload;
}

#endif
