#ifndef __INDEXMAP_H
#define __INDEXMAP_H

#include "cfdata.h"
/*
  This class does mapping from indices to BC-ified indicies to help 
  CFMat_BC and CF_Data BC.
 */
typedef int index_t;
typedef struct IndexMap {
  int start,end; // Of the inputs
  int Nsys; // Of the output
  int * map;
} indexmap_t;

void IndexMap_New(indexmap_t * self,
		  index_t istart, index_t iend,
				   index_t *BCs, index_t NBC);
index_t IndexMap_Get(indexmap_t * self, index_t i);
void IndexMap_Destroy(indexmap_t * self);
void IndexMap_Set_Values(indexmap_t * self, cfdata_t * cfbc,
			 cfdata_t * orig);
void IndexMap_Get_Values(indexmap_t * self, cfdata_t * cfbc,
			 cfdata_t * orig);
void IndexMap_Push(indexmap_t * self, cfdata_t * cfbc, cfdata_t * orig);
void IndexMap_Pull(indexmap_t * self, cfdata_t * cfbc, cfdata_t * orig);
#endif
