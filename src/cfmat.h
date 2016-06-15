#ifndef __CFMAT_H
#define __CFMAT_H

#include "kernel.h"
#include "hypergraph.h"

typedef struct _CFMAT_VTABLE_t _CFMAT_VTABLE_t;
typedef struct cfmat_t {
  int N;
  int own;
  void * data;
  const _CFMAT_VTABLE_t * vtable;
} cfmat_t;
struct _CFMAT_VTABLE_t {
  real_t * (*Place)(cfmat_t * self, int ln, int * ldofs,int rn, int * rdofs, real_t * vals);
  void (*Destroy)(cfmat_t * self);
  void (*Wipe)(cfmat_t * self);
  void (*Finalize)(cfmat_t * self);
};

real_t * CFMat_Place(cfmat_t * self, int ln, int * ldofs,int rn, int * rdofs, real_t * vals);
void CFMat_Destroy(cfmat_t * self);
void CFMat_Wipe(cfmat_t * self);
void CFMat_Finalize(cfmat_t * self);

#endif
