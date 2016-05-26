#ifdef USE_LIS
#ifndef __CFDATA_LIS_H
#define __CFDATA_LIS_H

#include "cfdata.h"

#include <lis.h>

#define CFData_LIS_Data(x) ((LIS_VECTOR)((x)->data))
void CFData_LIS_New_From_Ptr(cfdata_t * self, int N, LIS_VECTOR lvec);
void CFData_LIS_New(cfdata_t * self, int N);
#endif
#endif
