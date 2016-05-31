#ifdef USE_PETSC

#include "cfmat_petsc.h"
#include <petscmat.h>

#define data(x) CFMat_PETSc_Data(x)
real_t * CFMat_PETSc_Place(cfmat_t * self,
			   int n, int * dofs, real_t * ker_out) {
  MatSetValues(data(self), n,dofs, n,dofs,  ker_out, ADD_VALUES);
  return ker_out + n*n;
}
void CFMat_PETSc_Destroy(cfmat_t * self) {
  if(self->own) MatDestroy((Mat*)&(self->data));
}
void CFMat_PETSc_Wipe(cfmat_t * self) {
  MatZeroEntries(data(self));

}
void CFMat_PETSc_Finalize(cfmat_t * self) {
    MatAssemblyBegin(data(self),MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(data(self),MAT_FINAL_ASSEMBLY);
}

const _CFMAT_VTABLE_t CFMat_PETSc_vtable = {
  .Place = &CFMat_PETSc_Place,
  .Destroy = &CFMat_PETSc_Destroy,
  .Wipe = &CFMat_PETSc_Wipe,
  .Finalize = &CFMat_PETSc_Finalize
};
void CFMat_PETSc_New(cfmat_t * self, int N) {
  self->vtable = &CFMat_PETSc_vtable;
  self->N = N;
  self->own = 1;

  MatCreate(MPI_COMM_WORLD,(Mat*)&(self->data));
  MatSetSizes(data(self), PETSC_DECIDE,PETSC_DECIDE, N,N);
  //MatSetType(data(self)->K, MATMPIAIJ);
  MatSetFromOptions(data(self));
  MatSetUp(data(self));
}
void CFMat_PETSc_New_From_Ptr(cfmat_t * self, int N, Mat payload) {
  self->vtable = &CFMat_PETSc_vtable;
  self->N = N;
  self->own = 0;
  self->data = payload;
}

#endif
