#include "indexmap.h"

#include <stdlib.h>
#include <assert.h>

int compare(const void *a, const void *b) {
  if( *(index_t*)a > *(index_t*)b ) {
    return -1;
  } else if( *(index_t*)a > *(index_t*)b ) {
    return 1;
  } else {
    return 0;
  }
}
void IndexMap_New(indexmap_t * self,
		  index_t istart, index_t iend,
		  index_t *BCs, index_t NBC) {
  qsort(BCs,NBC, sizeof(index_t),compare);
  self->map = malloc( sizeof(index_t) * (iend-istart) );
  self->end = iend;
  self->start = istart;
  // This algorithm has an error if there are duplicates! I should detect them
  index_t b=istart, i=0;
  index_t N_duplicates; // This is to double-check the algorithm
  // i indexes the matrix, b indexes the BCs
  for(i=istart; i < iend; i++) {
    if(BCs[b] == i) {
      // Mark a BC
      self->map[i-istart] = -1;
      // Detect duplicates: they'll mess up this algorithm
      while( i!=BCs[b] ) {
	b++;
	N_duplicates++;
      }
    } else {
      self->map[i-istart] = i-istart-b;
      b++;
    }
  }
  // Check that it all added up
  assert( b==iend-istart-NBC +N_duplicates );
}
index_t IndexMap_Get(indexmap_t * self, index_t i) {
  if( i>=self->start && i<self->end ) {
    return self->map[i];
  } else {
    return -1;
  }
}
void IndexMap_Destroy(indexmap_t * self) {
  free(self->map);
}
