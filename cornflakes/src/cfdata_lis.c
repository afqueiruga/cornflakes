#ifdef USE_LIS
#include "cfdata_lis.h"


#define data(x) CFData_LIS_Data(x)

real_t * CFData_LIS_Place(cfdata_t * self,
			    int n, int * dofs, real_t * ker_out) {
  lis_vector_set_values(LIS_ADD_VALUE, n,dofs, ker_out, data(self));
  return ker_out + n;
}
void CFData_LIS_Destroy(cfdata_t * self) {
  if(self->own) lis_vector_destroy(data(self)); 
}
void CFData_LIS_Wipe(cfdata_t * self) {
  lis_vector_set_all( 0.0, data(self));
}
void CFData_LIS_Finalize(cfdata_t * self) {
  // Does it not need to?
}

void CFData_LIS_Get_Values(cfdata_t * self, int ndof, int *dofs, real_t * vals)
{
  int i;
  for(i=0;i<ndof;i++) {
    lis_vector_get_value((LIS_VECTOR)self->data, dofs[i], vals+i);
  }
}

const _CFDATA_VTABLE_t cfdata_lis_vtable = {
  .Get_Values = CFData_LIS_Get_Values,
  .Place = &CFData_LIS_Place,
  .Destroy = &CFData_LIS_Destroy,
  .Wipe = &CFData_LIS_Wipe,
  .Finalize = &CFData_LIS_Finalize
};

void CFData_LIS_New(cfdata_t * self, int N) {
  self->vtable = &cfdata_lis_vtable;
  self->N = N;
  self->own = 1;
  // TODO: MPI
  lis_vector_create(/*MPI_COMM_WORLD*/ 0,(LIS_VECTOR*)&self->data);
  lis_vector_set_size(data(self), 0,N);
}
void CFData_LIS_New_From_Ptr(cfdata_t * self, int N, LIS_VECTOR lvec) {
  self->vtable = &cfdata_lis_vtable;
  self->N = N;
  self->own = 0;
  self->data = lvec;
}

#endif
