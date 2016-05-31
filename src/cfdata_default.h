#ifndef __CFDATA_DEFAULT_H
#define __CFDATA_DEFAULT_H

#include "cfdata.h"

#define CFData_Default_Data(x) ( (real_t*)(x)->data )
void CFData_Default_New_From_Ptr(cfdata_t * self, int N, real_t * payload);
void CFData_Default_New(cfdata_t * self, int N);

#endif
