#ifndef __CFDATA_BC_H
#define __CFDATA_BC_H

#include "cfdata.h"
#include "indexmap.h"

/* 
 A lot of routines in the interface for this class are unimplemented!
 It is meant to be only a pseudo-target, not a data container! Manipulate 
 the pointer to R instead (or, R should be pointing to another CFData
 you are manipulating directly.)
 */
typedef struct cfdata_bc_data {
  cfdata_t * R;
  indexmap_t * map;
} cfdata_bc_data_t;
#define CFData_BC_Data(x) ( (cfdata_bc_data_t*)(x)->data )
void CFData_BC_New(cfdata_t * self,
		   cfdata_t * R, indexmap_t * map);
#endif
