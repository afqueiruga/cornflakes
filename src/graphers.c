#include "graphers.h"

#include "stdlib.h"
#include "math.h"

void Build_Pair_Graph(hypergraph_t * hg, int Npart, int dim, real_t * x, real_t cutoff) {
  spatialhash_t sh;
  Hypergraph_Alloc(hg,1);
  Build_New_Hash(&sh, Npart,dim,x, cutoff);
  void action(int a, int b) {
    if(dist(dim, x+dim*a,x+dim*b)<=cutoff) {
      int v[2] = {a,b};
      Hypergraph_Push_Edge(hg,2,v);
    }
  }
  SpatialHash_Scanall(&sh,x,action);
  SpatialHash_destroy(&sh);
}
void Build_Pair_Graph_2Sets(hypergraph_t * hg,
			    int Npart, int dim, real_t * x,
			    int Nparty, int dimy, real_t * y,
			    real_t cutoff) {
  spatialhash_t sh;
  Hypergraph_Alloc(hg,1);
  Build_New_Hash(&sh, Npart,dim,x, cutoff);
  void action(int a, int b) {
    if(dist(dim, y+dimy*a,x+dim*b)<=cutoff) {
      int v[2] = {a,b};
      Hypergraph_Push_Edge(hg,2,v);
    }
  }
  for(int A=0; A<Nparty; A++) {
    SpatialHash_ScanPt(&sh, y+dimy*A, action);
  }
  SpatialHash_destroy(&sh);
}

void Build_Proximity_Graph_Uniform(hypergraph_t * hg, int Npart, int dim, real_t * x, real_t cutoff) {
  spatialhash_t sh;
  int A;
  Hypergraph_Alloc(hg,1);
  Build_New_Hash(&sh, Npart,dim,x, cutoff);

  int listbuf = 20;
  int * list = malloc(sizeof(int)*listbuf); 
  int Nlist=0;
  void action(int FOO, int b) {
    if( b!=A &&  dist(dim, x+dim*A,x+dim*b)<=cutoff) {
      if (Nlist >= listbuf) {
	listbuf *= 2;
	list = realloc(list,sizeof(int)*listbuf);
      }
      list[Nlist] = b;
      Nlist++;
    }
  }

  for(A=0; A<Npart; A++) {
    list[0] = A;
    Nlist = 1;
    SpatialHash_ScanPt(&sh, x+dim*A, action);
    Hypergraph_Push_Edge(hg,Nlist,list);
  }
  free(list);
  SpatialHash_destroy(&sh);
}
void Build_Proximity_Graph_2Sets_Uniform(hypergraph_t * hg,
					 int Npart, int dim, real_t * x,
					 int Nparty, int dimy, real_t * y,
					 real_t cutoff) {
  spatialhash_t sh;
  int A;
  Hypergraph_Alloc(hg,1);
  Build_New_Hash(&sh, Npart,dim,x, cutoff);

  int listbuf = 20;
  int * list = malloc(sizeof(int)*listbuf); 
  int Nlist=0;
  void action(int FOO, int b) {
    if( dist(dim, y+dim*A,x+dim*b)<=cutoff) {
      if (Nlist >= listbuf) {
	listbuf *= 2;
	list = realloc(list,sizeof(int)*listbuf);
      }
      list[Nlist] = b;
      Nlist++;
    }
  }

  for(A=0; A<Nparty; A++) {
    list[0] = A;
    Nlist = 1;
    SpatialHash_ScanPt(&sh, y+dimy*A, action);
    Hypergraph_Push_Edge(hg,Nlist,list);
  }
  free(list);
  SpatialHash_destroy(&sh);
}



void Build_Proximity_Graph_Variable(hypergraph_t * hg,
				    int Npart, int dim, real_t * x,
				    real_t * r)
{
  spatialhash_t sh;
  int A;
  // Determine the maximum radius to use in the hash table
  real_t cutoff = 0.0;
  for(A=0;A<Npart;A++) if( r[A] > cutoff ) cutoff = r[A] ;
  cutoff *= 2.0001; // Just a wee bit extra for roundoff, lolz.
  
  Hypergraph_Alloc(hg,1);

  Build_New_Hash(&sh, Npart,dim,x, cutoff);

  int listbuf = 20;
  int * list = malloc(sizeof(int)*listbuf); 
  int Nlist=0;
  void action(int FOO, int b) {
    if( b!=A &&  dist(dim, x+dim*A,x+dim*b)<= r[A] ) {
      if (Nlist >= listbuf) {
	listbuf *= 2;
	list = realloc(list,sizeof(int)*listbuf);
      }
      list[Nlist] = b;
      Nlist++;
    }
  }

  for(A=0; A<Npart; A++) {
    list[0] = A;
    Nlist = 1;
    SpatialHash_ScanPt(&sh, x+dim*A, action);
    Hypergraph_Push_Edge(hg,Nlist,list);
  }
  free(list);
  SpatialHash_destroy(&sh);
}

void Build_Proximity_Graph_2Sets_Variable(hypergraph_t * hg,
					  int Npart, int dim, real_t * x,
					  int Nparty,int dimy,real_t * y,
					  real_t * r)
{
  spatialhash_t sh;
  int A;
  // Determine the maximum radius to use in the hash table
  real_t cutoff = 0.0;
  for(A=0;A<Nparty;A++) if( r[A] > cutoff ) cutoff = r[A] ;
  cutoff *= 2.0001; // Just a wee bit extra for roundoff, lolz.
  
  Hypergraph_Alloc(hg,1);

  Build_New_Hash(&sh, Npart,dim,x, cutoff);

  int listbuf = 20;
  int * list = malloc(sizeof(int)*listbuf); 
  int Nlist=0;
  void action(int FOO, int b) {
    if( dist(dim, y+dimy*A,x+dim*b)<= r[A] ) {
      if (Nlist >= listbuf) {
	listbuf *= 2;
	list = realloc(list,sizeof(int)*listbuf);
      }
      list[Nlist] = b;
      Nlist++;
    }
  }

  for(A=0; A<Nparty; A++) {
    list[0] = A;
    Nlist = 1;
    SpatialHash_ScanPt(&sh, y+dimy*A, action);
    Hypergraph_Push_Edge(hg,Nlist,list);
  }
  free(list);
  SpatialHash_destroy(&sh);
}


void Build_Proximity_Graph_Given_Length(hypergraph_t * hg,
					int Npart, int dim, real_t * x,
					int N_desired, real_t cutoff,
					real_t * r) //r is an output!
{
  spatialhash_t sh;
  int A;
  Hypergraph_Alloc(hg,1);
  Build_New_Hash(&sh, Npart,dim,x, cutoff);

  real_t        dists[N_desired+3];
  hypervertex_t list[N_desired+3];
  int Nlist; // How many we've put into the list so far
  
  void action(int FOO, int b) {
    // Calculate the distances
    if(b==A) return;
    real_t rad = dist(dim, x+dim*A,x+dim*b);
    // If the list is full, he might not fit
    if(Nlist == N_desired+2) {
      if( rad > dists[N_desired+1] )
	return;
    }
    // The list is empty, so I can't really search it!
    if(Nlist==1) {
      dists[1] = rad;
      list[1] = b;
      Nlist=2;
      return;
    }
    // Where does he go?
    int left=0,right=Nlist-1;
    // Find the spot:
    if( rad < dists[ 1 ] ) {
      right = 0;
    } else {
      do {
	if( rad > dists[ (left+right)/2 +1 ]) {
	  left = (left+right)/2;
	} else {
	  right = (left+right)/2;
	}
      } while(right - left > 1);
    }
    //Truncate
    if(Nlist<N_desired+2) Nlist++;
    // Shift the list and insert me
    for(int idx = Nlist-1; idx>=right; idx--) {
      list[idx+1 +1] = list[idx +1];
      dists[idx+1 +1] = dists[idx +1];
    }
    list[ right +1] = b;
    dists[right +1] = rad;
  }
  
  // Loop through all of the particles
  for(A=0; A<Npart; A++) {
    list[0] = A;
    Nlist = 1;
    SpatialHash_ScanPt(&sh, x+dim*A, action);
    Hypergraph_Push_Edge(hg,(N_desired > Nlist-1 ? Nlist : N_desired+1),list);

    r[A] = (dists[Nlist-1]+dists[Nlist-2])/2.0 + 1.0e-10;
  }

  SpatialHash_destroy(&sh);
}



void Build_Proximity_Graph_2Sets_Given_Length
(
 hypergraph_t * hg,
 int Npart, int dim, real_t * x,
 int Nparty,  int dimy, real_t * y,
 int N_desired, real_t cutoff,
 real_t * r) //r is an output!
{
  spatialhash_t sh;
  int A;
  Hypergraph_Alloc(hg,1);
  Build_New_Hash(&sh, Npart,dim,x, cutoff);

  real_t        dists[N_desired+3];
  hypervertex_t list[N_desired+3];
  int Nlist; // How many we've put into the list so far
  
  void action(int FOO, int b) {
    // Calculate the distances
    real_t rad = dist(dim, y+dim*A,x+dim*b);
    // If the list is full, he might not fit
    if(Nlist == N_desired+2) {
      if( rad > dists[N_desired+1] )
	return;
    }
    // The list is empty, so I can't really search it!
    if(Nlist==1) {
      dists[1] = rad;
      list[1] = b;
      Nlist=2;
      return;
    }
    // Where does he go?
    int left=0,right=Nlist-1;
    // Find the spot:
    if( rad < dists[ 1 ] ) {
      right = 0;
    } else {
      do {
	if( rad > dists[ (left+right)/2 +1 ]) {
	  left = (left+right)/2;
	} else {
	  right = (left+right)/2;
	}
      } while(right - left > 1);
    }
    //Truncate
    if(Nlist<N_desired+2) Nlist++;
    // Shift the list and insert me
    for(int idx = Nlist-1; idx>=right; idx--) {
      list[idx+1 +1] = list[idx +1];
      dists[idx+1 +1] = dists[idx +1];
    }
    list[ right +1] = b;
    dists[right +1] = rad;
  }
  
  // Loop through all of the particles
  for(A=0; A<Nparty; A++) {
    list[0] = A;
    Nlist = 1;
    SpatialHash_ScanPt(&sh, y+dim*A, action);
    Hypergraph_Push_Edge(hg,(N_desired > Nlist-1 ? Nlist : N_desired+1),list);

    r[A] = (dists[Nlist-1]+dists[Nlist-2])/2.0 + 1.0e-10;
  }

  SpatialHash_destroy(&sh);
}



void Add_Edge_Vertex(hypergraph_t * hgnew, hypergraph_t * hgold, int offset) {
  int i,j,k;
  Hypergraph_Alloc(hgnew, hgold->n_types);

  int edgenum = offset;
  for(i=0;i<hgold->n_types;i++) {
    hyperedges_t * he = hgold->he + i;
    hypervertex_t edge[ he->l_edge + 1];
    hypervertex_t * e;
    for(j=0;j<he->n_edge; j++) {
      e = Hyperedges_Get_Edge(he, j);
      for(k=0;k<he->l_edge;k++) {
	edge[k] = e[k];
      }
      edge[he->l_edge] = edgenum;
      Hypergraph_Push_Edge(hgnew,he->l_edge+1,edge);
      edgenum++;
    }
  }
}


