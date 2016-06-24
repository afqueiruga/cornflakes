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
  self->Nsys = 0;
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
      self->map[i-istart] = self->Nsys;
      self->Nsys++; 
    }
  }
  // The error is if the _last_ DOF was a BC!!!!
  //self->Nsys = self->map[iend-istart-1] + 1;
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



void IndexMap_Push(indexmap_t * self, cfdata_t * cfbc,
		    cfdata_t * orig) {
  index_t i_bc,i_o;
  real_t v;
  for(i_o=0; i_o< orig->N; i_o++) {
    i_bc = IndexMap_Get(self,i_o);
    if(i_bc>=0) {
      CFData_Get_Values(cfbc, 1,&i_bc,&v);
      CFData_Place(orig, 1,&i_o,&v);
    }
  }
  CFData_Finalize(orig);
  //  CFData_Place(orig, data(self)->imap->Nsys,data(self)->map->map,
}
void IndexMap_Pull(indexmap_t * self, cfdata_t * cfbc,
		    cfdata_t * orig) {
  index_t i_bc,i_o;
  real_t v;
  /* for(i_bc=0;i_bc<data(self)->map->Nsys;i_bc++) { */
  for(i_o=0; i_o< orig->N; i_o++) {
    /* i_o = data(self)->map->map[i_bc]; */
    i_bc = IndexMap_Get(self,i_o);
    if(i_bc>=0) {
      CFData_Get_Values(orig, 1,&i_o,&v);
      CFData_Place(cfbc, 1,&i_bc,&v);
    }
  }
  CFData_Finalize(orig);
}
