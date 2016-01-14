#include "util.h"


#include "math.h"
real_t dist(int dim, real_t * x, real_t * y) {
  int i;
  double ac = 0.0;
  for(i=0;i<dim;i++) ac += (x[i]-y[i])*(x[i]-y[i]);
  return sqrt(ac);
}

real_t smoother(real_t r, real_t rad) {
  return 1.0;
}
void Interpolate(real_t * uold, real_t * Xold, int Nold,
		 real_t * unew, real_t * Xnew, int Nnew,
		 int udim, int xdim, real_t rad)
{
  spatialhash_t sh;
  Build_New_Hash(&sh, Nold,xdim,Xold, rad);

  int A, i;
  real_t weight;
  
  void action(int FOO, int b) {
    real_t w = smoother( dist(xdim, Xnew+xdim*A,Xold+xdim*b), rad);
    weight += w;
    for(i=0;i<udim;i++) unew[udim*A+i] += uold[udim*b+i] * w;
  }

  for(A=0; A<Nnew; A++) {
    weight = 0.0;
    SpatialHash_ScanPt(&sh, Xnew+xdim*A, action);
    //printf("%lf\n",weight);
    for(i=0; i<udim; i++) unew[udim*A+i] /= weight;
  }

  SpatialHash_destroy(&sh);
}
