#ifdef USE_LIS
#ifndef __CFMAT_LIS_2_H
#define __CFMAT_LIS_2_H

#include "cfmat.h"
#include <lis.h>

typedef struct CFMat_preLIS_data_t {
  LIS_MATRIX K;
  cfmat_t triplet;
} CFMat_preLIS_data_t;
#define CFMat_preLIS_Overhead(x) ((CFMat_preLIS_data_t*)((x)->data))
#define CFMat_preLIS_Data(x) ((LIS_MATRIX)(CFMat_preLIS_Overhead(x)->K))

void CFMat_preLIS_New(cfmat_t * self, int N, int onum,
		      kernel_t * ke, hypergraph_t * hg, int ndof);
void CFMat_preLIS_New_From_Ptr(cfmat_t * self, int N, void* payload);

#endif
#endif
