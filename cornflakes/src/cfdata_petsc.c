#ifdef USE_PETSC

#include "cfdata_petsc.h"
#include "petscvec.h"
#define data(x) CFData_PETSc_Data(x)
/* Member methods */
real_t * CFData_PETSc_Place(cfdata_t * self,
			    int n, int * dofs, real_t * ker_out) {
  VecSetValues(data(self), n,dofs, ker_out, ADD_VALUES);
  return ker_out + n;
}
void CFData_PETSc_Destroy(cfdata_t * self) {
  if(self->own) VecDestroy((Vec*)&(self->data));
}
void CFData_PETSc_Wipe(cfdata_t * self) {
  VecSet(data(self),0.0);
}
void CFData_PETSc_Finalize(cfdata_t * self) {
  VecAssemblyBegin(data(self));
  VecAssemblyEnd(data(self));
}
void CFData_PETSc_Get_Values(cfdata_t * self, int ndof,int *dofs, real_t * vals) {
  VecGetValues((Vec)self->data, ndof,dofs,  vals);
}
void CFData_PETSc_Get_Ptr(cfdata_t * self, real_t **ptr) {
  VecGetArray(data(self), ptr);
}
void CFData_PETSc_Release_Ptr(cfdata_t * self, real_t **ptr) {
  VecRestoreArray(data(self),ptr);
}


/* vtable */
const _CFDATA_VTABLE_t cfdata_petsc_vtable = {
  .Get_Values = CFData_PETSc_Get_Values,
  .Place = &CFData_PETSc_Place,
  .Destroy = &CFData_PETSc_Destroy,
  .Wipe = &CFData_PETSc_Wipe,
  .Finalize = &CFData_PETSc_Finalize,
  .Get_Ptr = &CFData_PETSc_Get_Ptr,
  .Release_Ptr = &CFData_PETSc_Release_Ptr
};


/* Constructors */
void CFData_PETSc_New_From_Ptr(cfdata_t * self, int N, Vec payload) {
  self->vtable = &cfdata_petsc_vtable;
  self->own = 0;
  self->N = N;
  self->data = payload;
}
void CFData_PETSc_New(cfdata_t * self, int N) {
  self->vtable = &cfdata_petsc_vtable;
  self->own = 1;
  self->N = N;
  VecCreate(MPI_COMM_WORLD,(Vec*)&self->data);
  VecSetSizes(data(self), PETSC_DECIDE, N);
}
void CFData_PETSc_New_Full(cfdata_t * self, int N, MPI_Comm comm, Vec like) {
  self->vtable = &cfdata_petsc_vtable;
  self->own = 1;
  self->N = N;
  if(like) {
    VecDuplicate(like,(Vec*)&self->data);
  } else {
    VecCreate(comm,(Vec*)&self->data);
    VecSetSizes(data(self), PETSC_DECIDE, N);
  }
}
#endif
