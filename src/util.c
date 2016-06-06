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

void load_gmsh(real_t ** x, int * N, int gdim,
	       hypergraph_t ** hg,
	       char * fname, ...)
{
  va_list args;
  va_start( args, fname );
  char buf[1024];
  vsprintf(buf,fname,args);
  FILE * fh = fopen(buf,"r");
  if(fh==NULL) {
    printf("ERROR: UNABLE TO OPEN INPUT FILE %s\n",buf);
    exit(-1);
  }
  
  /* Open a gmsh msh file and return the vertices and the loaded hypergraphs.
     Each group is given its own hypergraph */
  int nn, ne, i,k,A;
  /* Read in the header */
  fgets( buf,1024,fh);
  printf("%s",buf);
  fgets(buf,1024,fh);
  printf("%s",buf);
  fgets(buf,1024,fh);
  printf("%s",buf);

  /* The nodes block */
  fgets(buf,1024,fh);
  printf("%s",buf);
  fgets(buf,1024,fh);
  printf("%s",buf);

  sscanf(buf,"%d\n",&nn);
  printf("There are %d nodes\n",nn);
  *N = nn;
  *x = malloc(sizeof(real_t)*gdim* *N);
  #define IX(A,i) (gdim*(A)+(i))
  for(A=0;A<nn;A++) {
    fgets(buf,1024,fh);
    switch(gdim) {
    case 1:
      sscanf(buf,"%*d %lf %*f %*f\n", &(*x)[IX(A,0)] );
      break;
    case 2:
      sscanf(buf,"%*d %lf %lf %*f\n", &(*x)[IX(A,0)], &(*x)[IX(A,1)] );
      break;
    case 3:
      sscanf(buf,"%*d %lf %lf %lf\n", &(*x)[IX(A,0)], &(*x)[IX(A,1)], &(*x)[IX(A,2)] );
      break;
    }
  }
  fgets(buf,1024,fh);
  /* End nodes block */

  /* The element block */
  Hypergraph_Alloc(hg[0],1);
  
  fgets(buf,1024,fh);
  printf("%s",buf);
  fgets(buf,1024,fh);
  printf("%s",buf);

  sscanf(buf,"%d\n",&ne);
  printf("There are %d elements\n",ne);
  for(A=0; A<ne; A++) {
    int eid, etype,ntag, groupnum;
    //fgets(buf,1024, fh);
    fscanf(fh,"%d %d %d ", &eid,&etype,&ntag);
    fscanf(fh,"%d ",&groupnum);
    for(i=0;i<ntag-1;i++) fscanf(fh,"%*d ");
    int ledge = etype;//Hack: Need a lookup table
    hypervertex_t egg[ledge]; 
    for(i=0;i<ledge;i++) {
      fscanf(fh, "%d ",egg+i);
      egg[i]--;
    }
    fscanf(fh,"\n");
    Hypergraph_Push_Edge(hg[0],ledge,egg);
  }
  fgets(buf,1024,fh);

  /* End element block */
  
  fclose(fh);
}

void write_vtk(real_t * x, int gdim, int N, hypergraph_t * hg,
	       char * names, real_t ** data, int * l_data, int Ndata,
	       char * cnames, real_t **cdata, int * l_cdata, int Ncdata,
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

  int A,i,dnum;
  
  fprintf(fh, "# vtk DataFile Version 2.0\nGraph connectivity\nASCII\nDATASET UNSTRUCTURED_GRID\n");

  /* Write points */
  fprintf(fh,"POINTS %d double\n",N);
  for(A=0; A<N; A++) {
    if(gdim==2)
      fprintf(fh,"%e %e 0\n",x[gdim*A+0],x[gdim*A+1]);
    else
      fprintf(fh,"%e %e %e\n",x[gdim*A+0],x[gdim*A+1],x[gdim*A+2]);
  }

  /* Write edges */
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
  case 4: hetype=9; break; // TODO: Toggle based on gdim
    //  case 4: hetype=10; break;
  case 8: hetype=12; break;
  default: hetype=1;
  }
  for(A=0;A<he->n_edge;A++) {
    fprintf(fh,"%d\n",hetype);
  }
  fprintf(fh,"\n\n");

  /* Write point data */
  // TODO: Using gdim as the stride instead of l_data should be a bug.
  // HOWEVER: The vectors will match gdim (l_data == gdim for a vector)
  if(Ndata>0) fprintf(fh,"POINT_DATA %d\n",N);
  for(dnum=0; dnum<Ndata; dnum++) {
    switch(l_data[dnum]) {
    case 3: // A 3D vector
      fprintf(fh,"VECTORS %c double\n",names[dnum]);
      for(A=0;A<N;A++) {
	fprintf(fh,"%e %e %e\n",data[dnum][gdim*A+0],data[dnum][gdim*A+1],data[dnum][gdim*A+2]);
      }
      break;
    case 2: // A 2D vector
      fprintf(fh,"VECTORS %c double\n",names[dnum]);
      for(A=0;A<N;A++) {
	fprintf(fh,"%e %e 0\n",data[dnum][gdim*A+0],data[dnum][gdim*A+1]);
      }
      break;
    case 9: // A 3D tensor
    case 4: // A 2D tensor
      // TODO: fill in
    case 1: // A scalar
    default: // The default case is scalar field of the first component
      fprintf(fh,"SCALARS %c double\nLOOKUP_TABLE default\n",names[dnum]);
      for(A=0;A<N;A++) {
	fprintf(fh,"%e\n", data[dnum][l_data[dnum]*A]);
      }
    } // end switch
  } // end for dnum
  
  /* Write cell data */
  int Nelem = he->n_edge;
  if(Ncdata>0) fprintf(fh,"CELL_DATA %d\n",Nelem);
  for(dnum=0; dnum<Ncdata; dnum++) {
    switch(l_cdata[dnum]) {
    case 3: // A 3D vector
      fprintf(fh,"VECTORS %c double\n",cnames[dnum]);
      for(A=0;A<Nelem;A++) {
	fprintf(fh,"%e %e %e\n",cdata[dnum][gdim*A+0],cdata[dnum][gdim*A+1],cdata[dnum][gdim*A+2]);
      }
      break;
    case 2: // A 2D vector
      fprintf(fh,"VECTORS %c double\n",cnames[dnum]);
      for(A=0;A<Nelem;A++) {
	fprintf(fh,"%e %e 0\n",cdata[dnum][gdim*A+0],cdata[dnum][gdim*A+1]);
      }
      break;
    case 9: // A 3D tensor
    case 4: // A 2D tensor
      // TODO: fill in
    case 1: // A scalar
    default: // The default case is scalar field of the first component
      fprintf(fh,"SCALARS %c double\nLOOKUP_TABLE default\n",cnames[dnum]);
      for(A=0;A<Nelem;A++) {
	fprintf(fh,"%e\n", cdata[dnum][l_cdata[dnum]*A]);
      }
    } // end switch
  } // end for dnum
  
  fclose(fh);
}

