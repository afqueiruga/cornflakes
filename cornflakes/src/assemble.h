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

  real_t * Viter;
  int * IIiter;
  int * JJiter;
} assemble_target_t;


void setup_targets(kernel_t * ke, assemble_target_t * att, hypergraph_t * hg, int ndof);
void destroy_targets(kernel_t * ke, assemble_target_t * att);

void assemble_targets(kernel_t * ke, hypergraph_t * hg,
		      dofmap_t ** dofmaps, real_t ** data,
		      assemble_target_t * att);
/* More deprecation */
void assemble_targets_dep(int ntarget, assemble_target_t * att,
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
