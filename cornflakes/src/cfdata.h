#ifndef __FIELD_H
#define __FIELD_H

#include "kernel.h"

typedef struct _CFDATA_VTABLE_t _CFDATA_VTABLE_t;
typedef struct cfdata_t {
  int own;
  int N;
  void * data;
  const _CFDATA_VTABLE_t * vtable;
} cfdata_t;
struct _CFDATA_VTABLE_t {
  real_t * (*Place)(cfdata_t * self, int n, int * dofs, real_t * vals);
  void (*Get_Values)(cfdata_t * self, int ndof,int * dofs, real_t * vals);
  void (*Wipe)(cfdata_t * self);
  void (*Finalize)(cfdata_t * self);
  void (*Destroy)(cfdata_t * self);
};

void CFData_Get_Values(cfdata_t * self, int ndof,int * dofs, real_t * vals);
real_t * CFData_Place(cfdata_t * self, int n, int * dofs, real_t * vals);
void CFData_Wipe(cfdata_t * self);
void CFData_Finalize(cfdata_t * self);
void CFData_Destroy(cfdata_t * self);
#if 0
#ifdef USE_PETSC
#include <petscvec.h>
void CFData_PETSc_New(cfdata_t * self, Vec pvec);
#endif
#ifdef USE_LIS
#include <lis.h>
void CFData_LIS_New(cfdata_t * self, LIS_VECTOR lvec);
#endif
#endif
#endif
