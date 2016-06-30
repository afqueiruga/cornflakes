#ifndef __CFMAT_CSR_H
#define __CFMAT_CSR_H

#include "cfmat.h"

typedef struct CFMat_CSR_data_t {
  int nnz;
  int * IA; // Length N+1
  int * JA; // Length nnz
  real_t * V; // Length nnz

} CFMat_CSR_data_t;
#define CFMat_CSR_Data(x) ((CFMat_CSR_data_t*)((x)->data))

void CFMat_CSR_New(cfmat_t * self, int N);
/* void CFMat_CSR_Add_Sparsity(cfmat_t * self, int onum, */
			    /* kernel_t * ke, hypergraph_t * hg, dmap_t ** dms); */
/* void CFMat_CSR_Set_Sparsity(cfmat_t * self); */
/* void CFMat_CSR_From_Array(cfmat_t * self, int ndof, */
			      /* real_t * V, int * II, int * JJ); */
#endif
