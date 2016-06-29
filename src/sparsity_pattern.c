#include "sparsity_pattern.h"

#include <stdlib.h>

void Sparsity_Init(sparsity_t * self, int N, int nrow_guess) {
  self->nnz = 0;
  self->N = N;
  self->Ibuild = malloc(N*sizeof(indexset_t*));
  for(int i=0; i<N; i++) {
    self->Ibuild[i] = malloc(sizeof(indexset_t));
    IndexSet_New(self->Ibuild[i], nrow_guess);
  }
}
void Sparsity_Add_NNZ(sparsity_t * self, index_t i, index_t j) {
  if(IndexSet_Insert(self->Ibuild[i],j)) {
    self->nnz++;
  }
}


void Sparsity_Print(sparsity_t * self) {
  for(int i=0; i<self->N; i++) {
    printf("%d :: ", i);
    for(int j=0; j<self->Ibuild[i]->n; j++) {
      printf("%d ",self->Ibuild[i]->table[j]);
    }
    printf("\n");
  }
}

void Sparsity_Destroy(sparsity_t * self) {
  for(int i=0; i<self->N; i++) {
    IndexSet_Destroy(self->Ibuild[i]);
    free(self->Ibuild[i]);
  }
  free(self->Ibuild);
}
