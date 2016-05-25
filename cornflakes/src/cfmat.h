#ifndef __CFMAT_H
#define __CFMAT_H

#include "kernel.h"
#include "hypergraph.h"

typedef struct _CFMAT_VTABLE_t _CFMAT_VTABLE_t;
typedef struct cfmat_t {
  int rank;
  int N;
  void * data;
  const _CFMAT_VTABLE_t * vtable;
} cfmat_t;
struct _TARGET_VTABLE_t {
  real_t * (*Place)(target_t * self, int n, int * dofs, real_t * vals);
  void (*Destroy)(target_t * self);
  void (*Wipe)(target_t * self);
  void (*Finalize)(target_t * self);
};

real_t * CFMat_Place(cfmat_t * selfint n, int * dofs, real_t * vals);
void CFMat_Destroy(target_t * self);
void CFMat_Wipe(target_t * self);
void CFMat_Finalize(target_t * self);

#endif
