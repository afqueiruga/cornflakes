#ifndef __FILTER_H
#define __FILTER_H

#include "kernel.h"
#include "hypergraph.h"
#include "dofmap.h"
#include "assemble.h"

void filter(kernel_t * ke, hypergraph_t * hg,
	    dofmap_t ** dofmaps, cfdata_t ** data,
	    hypergraph_t *filtered);

#endif
