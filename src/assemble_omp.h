#ifndef __ASSEMBLE_OMP_H
#define __ASSEMBLE_OMP_H

#include "assemble.h"

void assemble_omp(kernel_t * ke, hypergraph_t * hg,
		  dofmap_t ** dofmaps, cfdata_t ** data,
		  target_t * att);

#endif
