#include "indexset.h"

#include <stdlib.h>

void IndexSet_New(indexset_t * self, int nalloc) {
  self->n = 0;
  if(nalloc<1) nalloc=32; // Don't be stupid.
  self->nalloc = nalloc;
  self->table = malloc(sizeof(index_t)*nalloc);
}
int IndexSet_Insert(indexset_t * self, index_t i) {
  // Realloc if needed:
  if(self->nalloc == self->n) {
    self->table = realloc(self->table, sizeof(index_t)*2*self->nalloc);
    self->nalloc *= 2;
  }
  // Fringe case: It's empty, and the binary search will mess up
  if(self->n == 0) {
    self->table[0] = i;
    self->n = 1;
    return 1;
  }
  int a,b;
  // Fringe case: It is at the beginning
  if(self->table[0]==i) {
    return 0;
  }
  // Fringe case: It needs to go at the beginning
  if(self->table[0] > i) {
    b = 0;
  } else {
    // Search through the list
    a=0, b=self->n;
    do {
      if(i>self->table[(a+b)/2]) {
	a = (a+b)/2;
      } else {
	b = (a+b)/2;
      }
    } while(b - a > 1);
    // i should be at a.
    if(b < self->n && self->table[b] == i) {
      return 0;
    }
  }
  
  // insert i into the map.
  self->n += 1;
  for(int idx = self->n-1; idx>=b; idx--) {
    self->table[idx+1] = self->table[idx];
  }
  self->table[b] = i;
  return 1;
}
void IndexSet_Destroy(indexset_t * self) {
  free(self->table);
}

