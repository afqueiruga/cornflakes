#ifndef __ASSEMBLE_H
#define __ASSEMBLE_H

#include "kernel.h"
#include "hypergraph.h"
#include "dofmap.h"
#include "target.h"

void collect(real_t * ker_in, kernel_t * ke, hypervertex_t* edge, int l_edge,
	     dofmap_t ** dms, real_t ** data);

void assemble(kernel_t * ke, hypergraph_t * hg,
		      dofmap_t ** dofmaps, real_t ** data,
		      target_t * att);


/*
typedef struct assemble_target_t {
  int rank;
  real_t * V;
  int * II;
  int * JJ;

  real_t * Viter;
  int * IIiter;
  int * JJiter;
} assemble_target_t;
*/

//void setup_targets(kernel_t * ketarget_t * att, hypergraph_t * hg, int ndof);
//void destroy_targets(kernel_t * ke, target_t * att);

#endif
