#ifndef __GRAPHERS_H
#define __GRAPHERS_H

#include "assemble.h"
#include "spatialhash.h"
#include "hypergraph.h"
#include "util.h"

void Build_Pair_Graph(hypergraph_t * hg, int Npart, int dim, real_t * x, real_t cutoff);
void Build_Adjacency_Graph_Uniform(hypergraph_t * hg, int Npart, int dim, real_t * x, real_t cutoff);
void Add_Edge_Vertex(hypergraph_t * hgnew, hypergraph_t * hgold, int offset);

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