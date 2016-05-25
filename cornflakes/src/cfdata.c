#include "cfdata.h"

void CFData_Get_Values(cfdata_t * self, int ndof,int * dofs, real_t * vals) {
  self->vtable->Get_Values(self,ndof,dofs,vals);
}
real_t * CFData_Place(cfdata_t * self, int n, int * dofs, real_t * vals) {
  return self->vtable->Place(self,n,dofs,vals);
}
void CFData_Wipe(cfdata_t * self) {
  self->vtable->Wipe(self);
}
void CFData_Finalize(cfdata_t * self) {
  self->vtable->Finalize(self);
}
void CFData_Destroy(cfdata_t * self) {
  self->vtable->Destroy(self);
}



#if 0
#ifdef USE_PETSC

void CFData_PETSc_Get_Values(cfdata_t * self, int ndof,int *dofs, real_t * vals)
{
  VecGetValues((Vec)self->data, ndof,dofs,  vals);
}

const _CFDATA_VTABLE_t cfdata_petsc_vtable = {
  .Get_Values = CFData_PETSc_Get_Values
};
void CFData_PETSc_New(cfdata_t * self, Vec pvec) {
  self->vtable = &cfdata_petsc_vtable;

  self->data = pvec;
}
#endif

#ifdef USE_LIS
#include <lis.h>
void CFData_LIS_Get_Values(cfdata_t * self, int ndof, int *dofs, real_t * vals)
{
  int i;
  for(i=0;i<ndof;i++) {
    lis_vector_get_value((LIS_VECTOR)self->data, dofs[i], vals+i);
  }
}
const _CFDATA_VTABLE_t cfdata_lis_vtable = {
  .Get_Values = CFData_LIS_Get_Values
};
void CFData_LIS_New(cfdata_t * self, LIS_VECTOR lvec) {
  self->vtable = &cfdata_lis_vtable;
  self->data = lvec;
}
#endif
#endif
