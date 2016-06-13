#include <assert.h>
#include <stdlib.h>

#include <petsc.h>

#include "cornflakes.h"
#include "cfdata_petsc.h"
#include "cfmat_petsc.h"
int main(int argc, char **argv) {
  PetscInitialize(&argc,&argv,NULL,"Test for CFmat");
  int BCs[10] = { 0, 3,3, 5, 6, 7, 8 };
  int NBC = 7;
  int N = 20;
  indexmap_t imap;
  IndexMap_New(&imap, 0,20, BCs,NBC);

  cfmat_t K, KBC;
  cfdata_t u, R, RBC;
  CFMat_PETSc_New(&K,N - NBC);
  CFData_PETSc_New(&R,N - NBC);
  CFData_PETSc_New(&u,N);
  VecSet(CFData_PETSc_Data(&u), 1.0);
  CFMat_BC_New(&KBC, &K,&R,&u, &imap);

  int dofs[3] = {1,2,5};
  printf("%d\n",IndexMap_Get(&imap, 1));
  real_t vals[9] = {1.0,2.0,3.0,
		    4.0,5.0,6.0,
		    7.0,8.0,9.0};
  CFMat_Place(&KBC, 3,dofs,vals);
  CFMat_Finalize(&KBC);
  MatView(CFMat_PETSc_Data(&K),PETSC_VIEWER_STDOUT_SELF);
  VecView(CFData_PETSc_Data(&R),PETSC_VIEWER_STDOUT_SELF);
  CFMat_Destroy(&K);
  CFData_Destroy(&R);
  
  PetscFinalize();
  return 0;
}
