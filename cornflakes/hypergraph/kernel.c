#include "kernel.h"

void particle_kernel_calc(real_t* in,real_t* out) {
  
}

kernel_t particle_kernel_strct = {
  .ninp = 2,
  .inp = { {LOC_NODE,2},
	   {LOC_NODE,2},
	   {LOC_NODE,0},
	   {LOC_NODE,0},
	   {LOC_NODE,0} },
  .noutp = 1,
  .outp = { {LOC_NODE,2},
	    {LOC_NODE,0},
	    {LOC_NODE,0},
	    {LOC_NODE,0},
	    {LOC_NODE,0} },
  .assem = particle_kernel_calc
};
