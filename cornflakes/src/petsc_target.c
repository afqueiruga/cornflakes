#include "petsc_target.h"

#define data(x) (Target_PETSc_Data(x))
void Target_PETSc_Destroy(target_t * self) {
  switch(self.rank) {
    case 0:
    case 1:
      VecDestroy(&(data(self)->R));
      break;
    case 2:
      MatDestroy(&(data(self)->K));
      break;
    default:
      //
      break;
  }
}
void Target_PETSc_Wipe(target_t * self) {

}
real_t * Target_PETSc_Place(target_t * self,
			      int n, int * dofs, real_t * ker_out)
{

}
const _TARGET_VTABLE_t Table_PETSc_vtable = {
  .Place = &Target_PETSc_Place,
  .Destroy = &Target_PETSc_Destroy,
  .Wipe = &Target_PETSc_Wipe
};


 

void Target_PETSc_New(target_t * self, int onum,
		      kernel_t * ke, hypergraph_t * hg, int ndof) {

}
void Target_PETSc_From_Obj(target_t * self, int rank, void * obj) {

}
