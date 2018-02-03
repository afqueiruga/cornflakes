#ifndef __FILTER_H
#define __FILTER_H

#include "kernel.h"
#include "hypergraph.h"
#include "dofmap.h"
#include "assemble.h"

void filter(kernel_t * ke, hypergraph_t * hg,
			cfdata_t ** data, dofmap_t ** idofmaps,
			hypergraph_t *htrue, hypergraph_t *hfalse);

#endif
