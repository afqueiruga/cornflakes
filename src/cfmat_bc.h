#ifndef __CFMAT_BC_H
#define __CFMAT_BC_H

#include "cfmat.h"
#include "cfdata.h"
#include "indexmap.h"
/* 
   This is a wrapper for CFMat that applies BCs to a matrix
   It's really weird because adding entries to the matrix 
   also effects R by R' = R - K'u'
 */
typedef struct cfmat_bc_data {
  cfmat_t * K;
  cfdata_t * R; // Gets subtracted from
  cfdata_t * u; // Contains BCs
  indexmap_t * map;
} cfmat_bc_data_t;
#define CFMat_BC_Data(x) ((cfmat_bc_data_t*)(x)->data)
void CFMat_BC_New(cfmat_t * self,
		  cfmat_t * K, cfdata_t * R, cfdata_t * u,
		  indexmap_t * map);

#endif
