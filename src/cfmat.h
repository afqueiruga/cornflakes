#ifndef __CFMAT_H
#define __CFMAT_H

#include "kernel.h"
#include "hypergraph.h"
#include "sparsity_pattern.h"

typedef struct _CFMAT_VTABLE_t _CFMAT_VTABLE_t;
typedef struct cfmat_t {
  int N;
  int own;
  sparsity_t sparse;
  const _CFMAT_VTABLE_t * vtable;
  void * data;
} cfmat_t;
struct _CFMAT_VTABLE_t {
  void (*Add_Sparsity)(cfmat_t * self, int n, int *dofs);
  void (*Finalize_Sparsity)(cfmat_t * self);
  real_t * (*Place)(cfmat_t * self, int n, int * dofs, real_t * vals);
  void (*Set_Value)(cfmat_t * self,int i, int j, real_t v);
  void (*Destroy)(cfmat_t * self);
  void (*Wipe)(cfmat_t * self);
  void (*Finalize)(cfmat_t * self);
};

void CFMat_Add_Sparsity(cfmat_t * self, int n, int *dofs);
void CFMat_Finalize_Sparsity(cfmat_t * self);
real_t * CFMat_Place(cfmat_t * self,int n, int * dofs, real_t * vals);
void CFMat_Set_Value(cfmat_t * self,int i, int j, real_t v);
void CFMat_Destroy(cfmat_t * self);
void CFMat_Wipe(cfmat_t * self);
void CFMat_Finalize(cfmat_t * self);

#endif
