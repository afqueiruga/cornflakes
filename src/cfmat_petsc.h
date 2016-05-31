#ifdef USE_PETSC
#ifndef __CFMAT_PETSC_H
#define __CFMAT_PETSC_H


#include "cfmat.h"
#include <petsc.h>

#define CFMat_PETSc_Data(x) ((Mat)(x)->data)
void CFMat_PETSc_New(cfmat_t * self, int N);
void CFMat_PETSc_New_From_Ptr(cfmat_t * self, int rank, Mat payload);


#endif
#endif
