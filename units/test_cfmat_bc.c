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
  
  CFMat_BC_New(&KBC, &K,&R,&u, &imap);
  
  CFMat_Destroy(&K);
  CFData_Destroy(&R);
  
  PetscFinalize();
  return 0;
}
