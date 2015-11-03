#include "SpatialHash.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#ifndef SQ
#define SQ(x) ((x)*(x))
#define MIN(A,B) ((A)>(B)?(B):(A))
#define MAX(A,B) ((A)>(B)?(A):(B))
#endif



void SpatialHash_init(spatialhash_t * sh, int Npart, int dim,
		      real_t *start, real_t * end,real_t *h)
{
  int i;
  
  sh->dim = dim;
  sh->Npart = Npart;
  

  int celltot = 1;
  for(i=0;i<dim;i++) {
    sh->start[i] = (start)[i];
    sh->h[i] = (h)[i];
    sh->Ncell[i] = ceil( ( (end)[i]-(start)[i] )/ (h)[i] );
    celltot *= sh->Ncell[i];
  }
  printf("ncells: %d\n",celltot);
  sh->cells = (int*)malloc(sizeof(int)*celltot);

  sh->list = (int*)malloc(sizeof(int)*Npart);
  
  for(i=0; i<celltot; i++) {
    sh->cells[i] = -1;
  }
  for(i=0; i<Npart; i++) {
    sh->list[i] = -1;
  }
}

void SpatialHash_destroy(spatialhash_t * sh) {
  int i;
  free(sh->cells);
  free(sh->list);
  sh->Npart=0;
  
  for(i=0;i<sh->dim;i++) {
    sh->Ncell[i]=0;
  }
}

void SpatialHash_Push(spatialhash_t * sh, int A, real_t * x) {
#define HASH(p,i) ((int)( ((p)-sh->start[i])/sh->h[i] + 0.5 ))
#define CC(H) ( (H)[0]+(sh->dim>1?sh->Ncell[0]*(H)[1]:0)+(sh->dim>2?sh->Ncell[1]*sh->Ncell[0]*((H)[2]):0) )
  int i;
  int hs[sh->dim];
  for(i=0;i<sh->dim;i++) hs[i]=HASH((x)[i],i);
  printf("hs: %d,%d\n",hs[0],hs[1]);
  sh->list[A] = sh->cells[CC(hs)];
  sh->cells[CC(hs)] = A;
  printf("cc: %d,%d\n",CC(hs),sh->cells[CC(hs)]);
#undef HASH
#undef CC
}

void SpatialHash_print(spatialhash_t * sh) {
  int i,j;
  for(i=0;i<sh->Ncell[1];i++) {
    for(j=0;j<sh->Ncell[0];j++) {
      printf("%d ",sh->cells[i*sh->Ncell[0]+j]);
      
    }
    printf("\n");
  }
}

void Build_New_Hash(spatialhash_t * sh,  int Npart,int dim, real_t * x, real_t binsize)
{

  // #define CC(a,b,c) ( this->Ncell[1]*this->Ncell[2]*(a) + this->Ncell[1]*(b) + (c) )

  int A,i;
  /* Determine the limits */
  real_t start[dim], end[dim],h[dim];
  for(i=0; i<dim; i++) {
    start[i] = 1e10;
    end[i] = -1e10;
    h[i] = binsize;
  }
  for(A=0; A<Npart ;A++) {
    for(i=0; i<dim; i++) {
      start[i] = MIN(start[i],x[A*dim + i]);
      end[i] = MAX(end[i],x[A*dim + i]);
    }
  }
  for(i=0;i<Npart;i++) {
    start[i] -= 1.0*binsize;
    end[i] += 2.0*binsize;
  }
  /* Initialize the hash */
  SpatialHash_init(sh, Npart,dim, start,end,h);
  
  /* Build the hash */
  for(A=0;A<Npart;A++) {
    printf("%d",A);
    SpatialHash_Push(sh, A, x+dim*A);
  }
  for(i=0;i<dim;i++) printf("%lf, %lf, %lf\n",start[i],end[i],h[i]);
  SpatialHash_print(sh);
}

void SpatialHash_Scanall(spatialhash_t * sh, real_t * x,void act(int,int) )
{
  #define HASH(p,i) ((int)( ((p)-sh->start[i])/sh->h[i] + 0.5 ))

  int A,i, hs[sh->dim];
  // Loop accross it to make contact pairs
  for(A=0;A<sh->Npart;A++) {
    // std::cout << "d" << a <<std::endl;
    for(i=0;i<sh->dim;i++) { hs[i]=HASH((x)[A*sh->dim+i],i); }
    printf("%d : %d,%d\n",A,hs[0],hs[1]);
    // We scan in this pattern on the grid
    //    ooo   ooo   xxx
    //    ooo   oXx   xxx
    //    ooo   xxx   xxx
    // to save work.
    int off[sh->dim];
    for(i=0;i<sh->dim;i++) off[i]=0;
    
    SpatialHash_Scan_Cell(sh,act, A, hs, off, 1); // The center
    off[0] = 1; SpatialHash_Scan_Cell(sh,act, A, hs, off, 0);
    if(sh->dim > 1) { // idk why this would be used in just 1D... but RIGOR!
      off[1]=1; 
      off[0] = -1; SpatialHash_Scan_Cell(sh,act, A, hs, off, 0);
      off[0] =  0; SpatialHash_Scan_Cell(sh,act, A, hs, off, 0);
      off[0] =  1; SpatialHash_Scan_Cell(sh,act, A, hs, off, 0);
      // Do the next flat if 3D
      if(sh->dim > 2) {
	off[2]=1;
	for(i=-1;i<2;i++) {
	  off[1]=i;
	  off[0] = -1; SpatialHash_Scan_Cell(sh,act, A, hs, off, 0);
	  off[0] =  0; SpatialHash_Scan_Cell(sh,act, A, hs, off, 0);
	  off[0] =  1; SpatialHash_Scan_Cell(sh,act, A, hs, off, 0);
	}
      }
    }
  } // end loop over entries
  #undef HASH
}

void SpatialHash_Scan_Cell(spatialhash_t * sh, void act(int,int),
			   int a,
			   int * hs, int * off, int filter)
{
#define CC(H) ( (H)[0]+(sh->dim>1?sh->Ncell[0]*(H)[1]:0)+(sh->dim>2?sh->Ncell[0]*sh->Ncell[1]*((H)[2]):0) )

  int i, iter;
  int c[sh->dim];
  for(i=0;i<sh->dim;i++) c[i]=hs[i]+off[i];
  for(i=0;i<sh->dim;i++) printf("%d,",c[i]); printf("\n");
  for(i=0;i<sh->dim;i++) if(c[i]<0||c[i]>=sh->Ncell[i]) return;
  
  iter = sh->cells[CC(c)];
  printf("  %d,%d\n",CC(c),iter); 
  while(iter>-1) {
    printf("  %d\n",iter);
    //if( !filter || a>iter ) {
      if(iter!=a) {
	act(a,iter);
      }
      //}
    iter=sh->list[iter];
  }
#undef CC
}
