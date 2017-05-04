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

void Sparsity_Make_CSR(sparsity_t * self, index_t ** II, index_t **JS) {
  *II = malloc((self->N+1)*sizeof(index_t));
  *JS = malloc(self->nnz*sizeof(index_t));
  int itr = 0;
  (*II)[0] = 0;
  for(int i=0;i<self->N;i++) {
    int nrow = self->Ibuild[i]->n;
    (*II)[i+1] = itr+nrow;
    for(int j=0;j<nrow;j++) {
      (*JS)[itr + j] = self->Ibuild[i]->table[j];
    }
    itr += nrow;
  }
  /* assert( itr == self->nnz ); */
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
  //printf("Destroying sparsity of size %d with %d nnz\n",self->N, self->nnz);
  for(int i=0; i<self->N; i++) {
    IndexSet_Destroy(self->Ibuild[i]);
    free(self->Ibuild[i]);
  }
  free(self->Ibuild);
}
