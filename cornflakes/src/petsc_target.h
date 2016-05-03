#ifdef USE_PETSC
#ifndef __PETSC_TARGET_H
#define __PETSC_TARGET_H

#include "target.h"
#include <petscvec.h>
#include <petscmat.h>

void Target_PETSc_New(target_t * self, int onum,
		      kernel_t * ke, hypergraph_t * hg, int ndof,
		      MPI_Comm comm, Vec like);
void Target_PETSc_From_Obj(target_t * self, int rank, void * obj);
typedef struct Target_PETSc_data_t {
  int own;
  Vec R;
  Mat K;
} Target_PETSc_data_t;
#define Target_PETSc_Data(x) ((Target_PETSc_data_t*)((x)->data))

#endif
#endif
