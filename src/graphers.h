#ifndef __GRAPHERS_H
#define __GRAPHERS_H

#include "spatialhash.h"
#include "hypergraph.h"
#include "util.h"

void Build_Pair_Graph(hypergraph_t * hg, int Npart, int dim, real_t * x, real_t cutoff);
void Build_Proximity_Graph_Uniform(hypergraph_t * hg, int Npart, int dim, real_t * x, real_t cutoff);
void Build_Proximity_Graph_Variable(hypergraph_t * hg, int Npart, int dim, real_t * x, real_t * r);
void Build_Proximity_Graph_Given_Length(hypergraph_t * hg,
					int Npart, int dim, real_t * x,
					int N_desired, real_t cutoff,
					real_t * r);

void Build_Pair_Graph_2Sets(hypergraph_t * hg,
			    int Npart, int dim, real_t * x,
			    int Nparty, int dimy, real_t * y,
			    real_t cutoff);
void Build_Proximity_Graph_2Sets_Uniform(hypergraph_t * hg,
					 int Npart, int dim, real_t * x,
					 int Nparty, int dimy, real_t * y,
					 real_t cutoff);
void Build_Proximity_Graph_2Sets_Variable
    (hypergraph_t * hg,
     int Npart, int dim, real_t * x,
     int Nparty, int dimy, real_t * y,
     real_t * r);
void Build_Proximity_Graph_2Sets_Given_Length
    (hypergraph_t * hg,
     int Npart, int dim, real_t * x,
     int Nparty, int dimy, real_t * y,
     int N_desired, real_t cutoff,
     real_t * r); //r is an output!

void Add_Edge_Vertex(hypergraph_t * hgnew, hypergraph_t * hgold, int offset);


#endif
