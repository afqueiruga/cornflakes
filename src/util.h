#ifndef __UTIL_H
#define __UTIL_H

#include "spatialhash.h"
#include "hypergraph.h"


real_t dist(int dim, real_t * x, real_t * y);

void Interpolate(real_t * uold, real_t * Xold, int Nold,
		 real_t * unew, real_t * Xnew, int Nnew,
		 int udim, int xdim, real_t rad);

void load_gmsh(real_t ** x, int * N, int gdim,
	       hypergraph_t ** hg,
	       char * fname, ...);
void write_vtk(real_t * x, int gdim, int N, hypergraph_t * hg,
	       char * names, real_t ** data, int * l_data, int Ndata,
	       char * cnames, real_t **cdata, int * l_cdata, int Ncdata,

	      char * fname, ... );

#endif
