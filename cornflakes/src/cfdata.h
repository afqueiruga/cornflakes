#ifndef __FIELD_H
#define __FIELD_H

#include "kernel.h"

typedef struct _CFDATA_VTABLE_t _CFDATA_VTABLE_t;
typedef struct cfdata_t {
  void * data;
  const _CFDATA_VTABLE_t * vtable;
} cfdata_t;
struct _CFDATA_VTABLE_t {
  void (*Get_Values)(cfdata_t * self, int ndof,int * dofs, real_t * vals);
};

void CFData_Get(cfdata_t * self, int ndof,int * dofs, real_t * vals);
void CFData_Default_New(cfdata_t * self, real_t * payload);

#ifdef USE_PETSC
#include <petscvec.h>
void CFData_PETSc_New(cfdata_t * self, Vec pvec);
#endif
#ifdef USE_LIS
#include <lis.h>
void CFData_LIS_New(cfdata_t * self, LIS_VECTOR lvec);
#endif

#endif
