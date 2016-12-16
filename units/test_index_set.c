#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "cornflakes.h"

#define ADDANDCHECK(I) {\
  assert(IndexSet_Insert(&iset, I) == 1);\
  printf("nalloc=%d : ",iset.nalloc); \
  for(int i=0; i<iset.n; i++) printf("%d ",iset.table[i]); \
  printf("\n");						   \
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
  ADDANDCHECK(14);
  ADDANDCHECK(4);
  ADDANDCHECK(5);
  ADDANDCHECK(6);
  ADDANDCHECK(7);
  ADDANDCHECK(11);
  ADDANDCHECK(12);
  ADDANDCHECK(21);
  ADDANDCHECK(22);
  ADDANDCHECK(23);
  ADDANDCHECK(25);
  ADDANDCHECK(24);
  ADDANDCHECK(15);
  ADDANDCHECK(16);
  ADDANDCHECK(17);
  ADDANDCHECK(18);
  for(int i=40; i<100; i++) {
    ADDANDCHECK(i);
  }
  IndexSet_Destroy(&iset);
}
