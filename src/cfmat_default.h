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
#endif
