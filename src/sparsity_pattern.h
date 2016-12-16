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
void Sparsity_Make_CSR(sparsity_t * self, index_t ** II, index_t **JS);
void Sparsity_Print(sparsity_t * self);
void Sparsity_Destroy(sparsity_t * self);
#endif

/*
options:
A)
1.make sparsity_t with indexmap passed in
2. fill is a method of sparsity_t
2.assign to cfmats individually
+ no extra members on cfmat_t
- have to call fill_sparsity multiple times for each target mat
- more complicated to set bcs

B) 
1. cfmats own a sparsity_t, with no knowledge of bc
2. a global fill_sparsity routine that operates on a list of targets
3. add_spasity is a method on cfmat
+ do all of the targets at once
+ easier for user to call it on kernels during set up
+ easier for the python wrapper (i.e. me as a user)
+ simpler sparsity_t
- more methods on cfmat
+ cfmats need a finalize_sparsity method anyways
+ cfmat_bc would pass through corrected indices

which makes more sense with cfmat_bc?
which is easier?
*/
