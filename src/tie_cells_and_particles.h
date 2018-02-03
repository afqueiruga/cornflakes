#ifndef __TIE_CELLS_AND_PARTICLES_H
#define __TIE_CELLS_AND_PARTICLES_H

#include "assemble.h"
#include "util.h"

void Tie_Cells_and_Particles(hypergraph_t * hgnew,
			     hypergraph_t * mesh,
			     kernel_t * ke_circum,
			     kernel_t * ke_centroid,
			     kernel_t * ke_inside,
			     dofmap_t ** dofmaps,
			     cfdata_t ** data,
			     int Npart, int dim, real_t * x,
			     hypervertex_t * PV);

#endif
