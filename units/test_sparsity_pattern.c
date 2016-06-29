#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "cornflakes.h"

int main(int argc, char **argv) {
  sparsity_t sp;
  Sparsity_Init(&sp, 20, 4);

  Sparsity_Add_NNZ(&sp, 1,1);
  Sparsity_Add_NNZ(&sp, 1,2);
  Sparsity_Add_NNZ(&sp, 1,19);

  Sparsity_Add_NNZ(&sp, 7,19);
  Sparsity_Print(&sp);

  index_t *II, *JS;
  Sparsity_Make_CSR(&sp, &II,&JS);

  printf("II: ");
  for(int i=0; i<sp.N; i++) {
    printf("%d ",II[i]);
  }
  printf("\n");
  printf("JJ: ");
  for(int i=0; i<sp.nnz; i++) {
    printf("%d ",JS[i]);
  }
  printf("\n");
  
  free(II);
  free(JS);
  Sparsity_Destroy(&sp);
}
