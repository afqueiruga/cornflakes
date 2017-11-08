#ifndef __ASSEMBLE_H
#define __ASSEMBLE_H

#include "kernel.h"
#include "hypergraph.h"
#include "dofmap.h"
#include "target.h"
#include "cfdata.h"

void assemble2(kernel_t * ke, hypergraph_t * hg,
               cfdata_t ** data, dofmap_t ** idofmaps, // These are lined up
               void * targets, dofmap_t ** odofmaps); // These are also lined up

#endif
