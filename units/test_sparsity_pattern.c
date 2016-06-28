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

  Sparsity_Destroy(&sp);
}
