#include <assert.h>
#include <stdlib.h>

#include "cornflakes.h"

int main(int argc, char **argv) {
  int BCs[10] = { 0, 3, 5, 6, 7, 8 };
  int NBC = 6;
  
  indexmap_t imap;
  IndexMap_New(&imap, 0,20, BCs,NBC);

  for(int i=0; i<NBC; i++ ) {
    assert( IndexMap_Get(&imap, BCs[i]) == -1 );
  }
  assert( IndexMap_Get(&imap, -1) == -1 );
  assert( IndexMap_Get(&imap, 20) == -1 );

  assert( IndexMap_Get(&imap,  1) ==  0 );
  assert( IndexMap_Get(&imap, 19) == 13 );

  for(int i=0; i<20; i++) {
    printf(" %d ",IndexMap_Get(&imap, i));
  }
  printf("\n");

  return 0;
}
