#ifndef __GRAPHERS_H
#define __GRAPHERS_H

#include "spatialhash.h"
#include "hypergraph.h"

void Build_Pair_Graph(hypergraph_t * hg, int Npart, int dim, real_t * x, real_t cutoff);
void Build_Adjacency_Graph_Uniform(hypergraph_t * hg, int Npart, int dim, real_t * x, real_t cutoff);
void Add_Edge_Vertex(hypergraph_t * hgnew, hypergraph_t * hgold, int offset);

#endif
