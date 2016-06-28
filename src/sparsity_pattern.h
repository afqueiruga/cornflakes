#ifndef __SPARSITY_PATTERN_H
#define __SPARSITY_PATTERN_H

#include "indexset.h"

#include "kernel.h"
#include "dofmap.h"
#include "hypergraph.h"

typedef struct sparsity_str {
  indexset_t ** Ibuild;
  int N;
  int nnz;
  void * arena;
} sparsity_t;
  /* Analysis 
     We need to minimize memory reallocations and copies.
     We know N, but we don't know nnz, or nnz per row.
     nnz per row is small, O(1) e.g. 9, 21, etc = O(row)
     Need to insert J into I in O(1) time. 
     - Unsorted check is O(row)
     - Insert in order is O(log(row)) + O(row) for an array
     - Insert with linked lists needs mallocs... unless I arena it...
        It would be nnz * ( row  + malloc )
        Pro: No copies, Insert in order,
	Con: lots of mallocs, linear search time
      Best option I think
       [] ~~>x
       [] ~~> [  fixed-size payload  ][] ~~> [          ][x]
       [] ~~> [  of unsorted indices ][] ~~> [          ][x]
       And then sort them after
  */
void Sparsity_Init(sparsity_t * self, int N, int nrow_guess);
void Sparsity_Add_NNZ(sparsity_t * self, index_t i, index_t j);
void Sparsity_Fill(sparsity_t * self,
		   int onum, kernel_t * ke, hypergraph_t * hg, dofmap_t ** dms);
void Sparsity_Print(sparsity_t * self);
void Sparsity_Destroy(sparsity_t * self);
#endif
