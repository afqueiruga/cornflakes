#ifndef __ASSEMBLE_H
#define __ASSEMBLE_H

#include "kernel.h"
#include "hypergraph.h"
#include "dofmap.h"
#include "cfdata.h"
#include "cfmat.h"

void collect(kernel_t * ke, hypervertex_t* edge, int l_edge,
			 cfdata_t ** data,	dofmap_t ** dms,
			 real_t * ker_in);
void assemble(kernel_t * ke, hypergraph_t * hg,
			  cfdata_t ** data, dofmap_t ** idofmaps, // These are lined up
			  void * targets, dofmap_t ** odofmaps); // These are also lined up

#endif
