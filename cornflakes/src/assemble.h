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


#endif
