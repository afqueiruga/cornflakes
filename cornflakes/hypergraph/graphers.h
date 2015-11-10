#ifndef __GRAPHERS_H
#define __GRAPHERS_H

#include "SpatialHash.h"
#include "hypergraph.h"

void Build_Particle_Graph(hypergraph_t * hg, int Npart, int dim, real_t * x, real_t cutoff);

#endif
