#ifndef __ASSEMBLE_H
#define __ASSEMBLE_H

#include "kernel.h"
#include "hypergraph.h"
#include "dofmap.h"

typedef struct assemble_target_t {
  int rank;
  real_t * V;
  int * II;
  int * JJ;
} assemble_target_t;

void assemble_targets(int ntarget, assemble_target_t * att,
		      kernel_t * ke, hyperedges_t * hg,
		      int * outmap, real_t ** data);


/* DEPRECATED: USE TARGETS INSTEAD */
void assemble_vector(real_t * R, kernel_t * ke,hyperedges_t * hg, int * outmap, real_t ** data);
void assemble_matrix(int * II, int * JJ, real_t * KK,
		     kernel_t * ke, hyperedges_t * hg,
		     int * outmap, real_t ** data);
void assemble_vector_matrix(real_t * R,
			    int * II, int * JJ, real_t * KK,
			    kernel_t * ke, hyperedges_t * hg,
			    int * outmap, real_t ** data);

#endif
