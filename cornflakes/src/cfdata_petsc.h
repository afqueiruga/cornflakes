#ifdef USE_PETSC
#ifndef __CFDATA_PETSC_H
#define __CFDATA_PETSC_H

#include "cfdata.h"
#include <petsc.h>

#define CFData_PETSc_Data(x) ( (Vec)(x)->data )
void CFData_PETSc_New_From_Ptr(cfdata_t * self, int N, Vec payload);
void CFData_PETSc_New(cfdata_t * self, int N);


#endif
#endif
