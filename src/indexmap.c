#include "indexmap.h"

#include <stdlib.h>
#include <assert.h>

int compare(const void *a, const void *b) {
  if( *(index_t*)a > *(index_t*)b ) {
    return 1;
  } else if( *(index_t*)a < *(index_t*)b ) {
    return -1;
  } else {
    return 0;
  }
}
void IndexMap_New(indexmap_t * self,
		  index_t istart, index_t iend,
		  index_t *BCs, index_t NBC) {
  // Sort the BCs to use O(n) in the next loops
  // This routine's complexity: O(nbc log nbc) + O(nidx)
  qsort(BCs,NBC, sizeof(index_t),compare);
  // Need to scroll BCs up incase some of them are out of my ownership
  
  self->map = malloc( sizeof(index_t) * (iend-istart) );
  self->end = iend;
  self->start = istart;
  // This algorithm has an error if there are duplicates! I should detect them
  index_t b=0, i=0;
  index_t N_duplicates=0; // This is to double-check the algorithm
  // i indexes the matrix, b indexes the BCs
  for(i=istart; i < iend; i++) {
    if(b < NBC && BCs[b] == i) {
      // Mark a BC
      self->map[i-istart] = -1;
      b++;
      // Detect duplicates: they'll mess up this algorithm
      while( b<NBC && i==BCs[b] ) {
	N_duplicates++;
	b++;
      }
    } else {
      self->map[i-istart] = i-b;
      //b++;
    }
  }
  /*
  // Check that it all added up
  printf("%d %d %d\n",b,NBC,N_duplicates);
  for(int i=istart; i<iend; i++) {
    printf(" %d ",IndexMap_Get(self, i));
  }
  printf("\n");
  assert( self->map[iend-istart-1]==iend-istart-NBC +N_duplicates-1 );
  */
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
