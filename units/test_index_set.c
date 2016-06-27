#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "cornflakes.h"

#define ADDANDCHECK(I) {\
  assert(IndexSet_Insert(&iset, I) == 1);\
  for(int i=0; i<iset.n; i++) printf("%d ",iset.table[i]); printf("\n"); \
  assert(IndexSet_Insert(&iset, I) == 0);\
  }
  
int main(int argc, char **argv) {
  indexset_t iset;
  
  IndexSet_New(&iset, 20);
  ADDANDCHECK(1);
  ADDANDCHECK(2);
  ADDANDCHECK(8);
  ADDANDCHECK(9);
  ADDANDCHECK(3);
  ADDANDCHECK(13);
  ADDANDCHECK(0);
  IndexSet_Destroy(&iset);
}
