#ifndef __ASSEMBLE_H
#define __ASSEMBLE_H

#include "kernel.h"
#include "hypergraph.h"
  

void assemble_vector(real_t * R, kernel_t * ke,hypergraph_t * hg, int * outmap, real_t ** data);

#endif
