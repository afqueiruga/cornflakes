#ifndef __FILL_SPARSITY_H
#define __FILL_SPARSITY_H

#include "kernel.h"
#include "hypergraph.h"
#include "dofmap.h"
#include "target.h"

void fill_sparsity2(kernel_t * ke, hypergraph_t * hg,
				   void * targets,	dofmap_t ** odofmaps);
#endif
