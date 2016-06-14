#ifdef USE_LIS
#ifndef __CFMAT_LIS_H
#define __CFMAT_LIS_H

#include "cfmat.h"
#include <lis.h>
#define CFMat_LIS_Data(x) ((LIS_MATRIX)((x)->data))
void CFMat_LIS_New(cfmat_t * self, int N);
void CFMat_LIS_New_From_Ptr(cfmat_t * self, int N, void* payload);

#endif
#endif
