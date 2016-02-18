#include "util.h"

#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

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



void write_vtk(real_t * x, int gdim, int N, hypergraph_t * hg,
	       char * names, real_t ** data, int Ndata,
	       char * fname, ... )
{
  va_list args;
  va_start( args, fname );
  char buf[1024];
  vsprintf(buf,fname,args);
  FILE * fh = fopen(buf,"w");
  if(fh==NULL) {
    printf("ERROR: UNABLE TO OPEN OUTPUT FILE\n");
    exit(-1);
  }

  int A,i;
  
  fprintf(fh, "# vtk DataFile Version 2.0\nGraph connectivity\nASCII\nDATASET UNSTRUCTURED_GRID\n");

  // Write points
  fprintf(fh,"POINTS %d double\n",N);
  for(A=0; A<N; A++) {
    if(gdim==2)
      fprintf(fh,"%e %e 0\n",x[gdim*A+0],x[gdim*A+1]);
    else
      fprintf(fh,"%e %e %e\n",x[gdim*A+0],x[gdim*A+1],x[gdim*A+2]);
  }

  // Write edges
  hyperedges_t * he = hg->he+0;
  fprintf(fh,"\nCELLS %d %d\n",he->n_edge,he->n_edge*(he->l_edge+1));
  for(A=0;A<he->n_edge;A++) {
    fprintf(fh,"%d ",he->l_edge);
    for(i=0;i<he->l_edge;i++) fprintf(fh,"%d ",he->edges[he->l_edge*A+i]);
    fprintf(fh,"\n");
  }
  fprintf(fh,"\nCELL_TYPES %d\n",he->n_edge);
  int hetype = 0;
  switch(he->l_edge) { // From the manual
  case 1: hetype=1; break;
  case 2: hetype=3; break;
  case 3: hetype=5; break;
  case 4: hetype=9; break;
    //  case 4: hetype=10; break;
  case 8: hetype=12; break;
  default: hetype=1;
  }
  for(A=0;A<he->n_edge;A++) {
    fprintf(fh,"%d\n",hetype);
  }
  fprintf(fh,"\n\n");

  fclose(fh);
  // Write point data
}

