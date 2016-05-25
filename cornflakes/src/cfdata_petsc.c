#ifdef USE_PETSC
#include "cfdata_petsc.h"

#define data(x) CFData_PETSc_Data(x)
/* Member methods */
real_t * CFData_PETSc_Place(cfdata_t * self,
			    int n, int * dofs, real_t * ker_out) {
  VecSetValues(data(self)->R, n,dofs, ker_out, ADD_VALUES);
  return ker_out + n;
}
void CFData_PETSc_Destroy(cfdata_t * self) {
  if(self->own) VecDestroy(&(data(self)->R));
}
void CFData_PETSc_Wipe(cfdata_t * self) {
  VecSet(data(self)->R,0.0);
}
void CFData_PETSc_Finalize(cfdata_t * self) {
  VecAssemblyBegin(data(self)->R);
  VecAssemblyEnd(data(self)->R);
}
void CFData_PETSc_Get_Values(cfdata_t * self, int ndof,int *dofs, real_t * vals) {
  VecGetValues((Vec)self->data, ndof,dofs,  vals);
}
/* vtable */
const _CFDATA_VTABLE_t cfdata_petsc_vtable = {
  .Get_Values = CFData_PETSc_Get_Values,
  .Place = &CFData_PETSc_Place,
  .Destroy = &CFData_PETSc_Destroy,
  .Wipe = &CFData_PETSc_Wipe,
  .Finalize = &CFData_PETSc_Finalize
};
/* Constructors */
void CFData_PETSc_New_From_Ptr(cfdata_t * self, int N, Vec payload) {
  self->vtable = cfdata_petsc_vtable;
  self->own = 0;
  self->N = N;
  self->data = payload;
}
void CFData_PETSc_New(cfdata_t * self, int N) {
  self->vtable = cfdata_petsc_vtable;
  self->own = 1;
  self->N = N;
  VecCreate(MPI_COMM_WORLD,&(data(self)->R));
  VecSetSizes(data(self)->R, PETSC_DECIDE, ndof);
}
#endif
