#ifndef __ASSEMBLE_H
#define __ASSEMBLE_H

#include "kernel.h"
#include "hypergraph.h"
#include "dofmap.h"
#include "target.h"
#include "cfdata.h"

void collect(real_t * ker_in, kernel_t * ke, hypervertex_t* edge, int l_edge,
	     dofmap_t ** dms, cfdata_t ** data);
void place_targets(target_t * att,
		   kernel_t * ke,
		   real_t * ker_out, int len_ker_out,
		   dofmap_t ** dms,
		   hypervertex_t * edge, int l_edge);
void assemble(kernel_t * ke, hypergraph_t * hg,
	      dofmap_t ** dofmaps, cfdata_t ** data,
	      target_t * att);


void assemble2(kernel_t * ke, hypergraph_t * hg,
               cfdata_t ** data, dofmap_t ** idofmaps, // These are lined up
               target_t * att, dofmap_t ** odofmaps); // These are also lined up
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
