#ifdef USE_PETSC
#include "petsc_target.h"

#define data(x) (Target_PETSc_Data(x))
void Target_PETSc_Destroy(target_t * self) {
  if(data(self)->own) {
    switch(self->rank) {
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
  free(self->data);
}
void Target_PETSc_Wipe(target_t * self) {
    switch(self->rank) {
    case 0:
    case 1:
      VecSet(data(self)->R,0.0);
      break;
    case 2:
      MatZeroEntries(data(self)->K);
      break;
    default:
      break;
    }
}
real_t * Target_PETSc_Place(target_t * self,
			      int n, int * dofs, real_t * ker_out)
{
  switch(self->rank) {
  case 2:
    MatSetValues(data(self)->K, n,dofs, n,dofs,  ker_out, ADD_VALUES);
    return ker_out + n*n;
  case 1:
    VecSetValues(data(self)->R, n,dofs, ker_out, ADD_VALUES);
    return ker_out + n;
  default:
    VecSetValue(data(self)->R, 0, ker_out[0], ADD_VALUES);
    return ker_out + 1;
  }
}
void Target_PETSc_Finalize(target_t * self) {
  switch(self->rank) {
  case 2:
    MatAssemblyBegin(data(self)->K,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(data(self)->K,MAT_FINAL_ASSEMBLY);
    break;
  case 1:
  default:
    VecAssemblyBegin(data(self)->R);
    VecAssemblyEnd(data(self)->R);
    
  }
}
const _TARGET_VTABLE_t Table_PETSc_vtable = {
  .Place = &Target_PETSc_Place,
  .Destroy = &Target_PETSc_Destroy,
  .Wipe = &Target_PETSc_Wipe,
  .Finalize = &Target_PETSc_Finalize
};


 

void Target_PETSc_New(target_t * self, int onum,
		      kernel_t * ke, hypergraph_t * hg, int ndof,
		      MPI_Comm comm, Vec like) {
  int i, matsize, oplen;
  self->rank = ke->outp[onum].rank;
  self->data = malloc(sizeof(struct Target_PETSc_data_t));
  self->vtable = &Table_PETSc_vtable;
  data(self)->own = 1;
  switch(self->rank) {
    case 0:
      VecCreate(comm,&(data(self)->R));
      VecSetSizes(data(self)->R, PETSC_DECIDE, 1);
      break;
    case 1:
      if(like) {
	VecDuplicate(like,&(data(self)->R));
      } else {
	VecCreate(comm,&(data(self)->R));
	VecSetSizes(data(self)->R, PETSC_DECIDE, ndof);
      }
      break;
    case 2:
      matsize = 0;
      for( i=0; i<hg->n_types; i++ ) {
	oplen = kernel_outp_len(ke,ke->outp+onum,hg->he[i].l_edge);
	matsize += hg->he[i].n_edge * oplen;
      }
      // Fill in with matsize nnz
      MatCreate(comm,&(data(self)->K));
      MatSetSizes(data(self)->K, PETSC_DECIDE,PETSC_DECIDE, ndof,ndof);
      //MatSetType(data(self)->K, MATMPIAIJ);
      MatSetFromOptions(data(self)->K);
      MatSetUp(data(self)->K);
      break;
    default:
      printf("Error: PETSc backend can't do rank > 2\n");
    }
}
void Target_PETSc_From_Obj(target_t * self, int rank, void * obj) {
  self->rank = rank;
  data(self)->own = 0;
  switch(self->rank) {
  case 0:
  case 1:
    data(self)->R = (Vec)obj;
    break;
  case 2:
    data(self)->K = (Mat)obj;
    break;
  default:
    printf("Error: PETSc backend can't do rank > 2\n");
  }
}
#endif
