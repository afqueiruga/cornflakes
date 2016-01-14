#ifndef __UTIL_H
#define __UTIL_H

#include "spatialhash.h"

real_t dist(int dim, real_t * x, real_t * y);

void Interpolate(real_t * uold, real_t * Xold, int Nold,
		 real_t * unew, real_t * Xnew, int Nnew,
		 int udim, int xdim, real_t rad);

#endif
