#ifndef __INDEXMAP_H
#define __INDEXMAP_H

/*
  This class does mapping from indices to BC-ified indicies to help 
  CFMat_BC and CF_Data BC.
 */
typedef int index_t ;
typedef struct indexmap_str {
  int start,end;
  int * map;
} indexmap_t;
void IndexMap_New(indexmap_t * self,
		  index_t istart, index_t iend,
				   index_t *BCs, index_t NBC);
index_t IndexMap_Get(indexmap_t * self, index_t i);
void IndexMap_Destroy(indexmap_t * self);

#endif
