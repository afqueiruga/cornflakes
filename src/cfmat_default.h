#ifndef __CFMAT_DEFAULT_H
#define __CFMAT_DEFAULT_H

#include "cfmat.h"

typedef struct CFMat_Default_data_t {
  int nalloc;
  real_t * V;
  int * II;
  int * JJ;

  real_t * Viter;
  int * IIiter;
  int * JJiter;
} CFMat_Default_data_t;
#define CFMat_Default_Data(x) ((CFMat_Default_data_t*)((x)->data))

void CFMat_Default_New(cfmat_t * self, int onum,
			kernel_t * ke, hypergraph_t * hg, int ndof);
void CFMat_Default_From_Array(cfmat_t * self, int ndof,
			      real_t * V, int * II, int * JJ);

// Some of his members can be accessed by other implementations.
void CFMat_Default_Add_Sparsity(cfmat_t * self, int n, int *dofs);
void CFMat_Default_Finalize_Sparsity(cfmat_t * self);
void CFMat_Default_Finalize(cfmat_t * self);
#endif
