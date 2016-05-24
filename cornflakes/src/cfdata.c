#include "cfdata.h"

void CFData_Get(cfdata_t * self, int ndof,int * dofs, real_t * vals)
{
  self->vtable->Get_Values(self,ndof,dofs,vals);
}


void CFData_Default_Get_Values(cfdata_t * self, int ndof,int *dofs, real_t * vals)
{
  int k;
  for(k=0;k<ndof;k++) vals[k] = ((real_t*)self->data)[ dofs[k] ];
}
const _CFDATA_VTABLE_t cfdata_default_vtable = {
  .Get_Values = CFData_Default_Get_Values
};
void CFData_Default_New(cfdata_t * self, real_t * payload) {
  self->vtable = &cfdata_default_vtable;
  self->data = payload;
}

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
