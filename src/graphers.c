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



/*
 * This routine assigns the particles in x to a cell in mesh
 * If the particles isn't inside of a mesh, it would be in the output.
 * It returns a graph that looks like
 * [ N1 ... don't care... Nn ] P0 P1 P2 ... Pn
 * where it is agnostic to the edge of the graph mesh, and just copies it over.
 * 
 * All it needs is a mesh, (doesn't really need to be a mesh) and 3 kerenls that 
 * gives the circumradius (or an upper bound), gives the centroid (or something 
 * close enough), and tests wether or not a list of points is inside of it.
 * They all need to have the same dofmap call signatures.
 */
#include "cfdata_default.h"
void Tie_Cells_and_Particles(hypergraph_t * hgnew,
			     hypergraph_t * mesh,
			     kernel_t * ke_circum,
			     kernel_t * ke_centroid,
			     kernel_t * ke_inside,
			     dofmap_t ** dofmaps,
			     cfdata_t ** data,
			     int Npart, int dim, real_t * x,
			     hypervertex_t * PV)
{
  int i, j;
  int n_cells = mesh->he[0].n_edge; // TODO: BUG. Hypergraph needs a total_edge_count method

  /* Setup storage for the results */
  target_t att[1];
  cfdata_t hs;
  CFData_Default_New(&hs, n_cells);
  Target_New_From_Ptr(att+0, 1,&hs);
  //setup_targets(ke_circum,att, mesh,n_cells);
  
  /* Calculate the maximum circum radius */
  assemble(ke_circum,mesh, dofmaps,data, att);
  real_t hmax = 0.0;
  for(int i=0; i<n_cells; i++) {
    if(CFData_Default_Data(&hs)[i] > hmax) {
      hmax = CFData_Default_Data(&hs)[i];
    }
  }
  hmax = sqrt(hmax);
  //printf("hmax: %lf\n",hmax);
  /* Build a hash using this circumradius and fill it with the particles */
  spatialhash_t sh;
  Build_New_Hash(&sh, Npart,dim,x, hmax); // There's a crash here
  int * found = calloc( sizeof(int),Npart);
  
  /* Now loop against the cells */
  hypervertex_t * edge;
  hyperedges_t * he;
  int hex;
  int listbuf = 20;
  int * list = calloc(sizeof(int),listbuf); 
  for(he=mesh->he; he<mesh->he+mesh->n_types ; he++) {
    int l_edge = he->l_edge;
    int len_ker_centroid = kernel_inps_len(ke_centroid, l_edge); 
    real_t ker_in[ len_ker_centroid ];
    real_t centroid[dim]; // Kernel doesn't have the choice... the output has to look like this, per the contract.
    
    for(hex=0; hex<he->n_edge; hex++) {
      edge = Hyperedges_Get_Edge(he, hex);

      /* Calculate the centroid of the cell */
      collect(ker_in,ke_centroid,edge,l_edge, dofmaps,data);
      for(i=0;i<dim;i++) centroid[i]=0.0;
      ke_centroid->eval(l_edge, ker_in, centroid);
      //printf("center: %lf %lf  ",centroid[0],centroid[1]);
      /* build a list of all the particles _near_ a cell */
      int Nlist=0;
      void action(int FOO, int b) {
	if (Nlist >= listbuf) {
	  listbuf *= 2;
	  list = realloc(list,sizeof(int)*listbuf);
	}
	list[Nlist] = b;
	Nlist++;
      }
      SpatialHash_ScanPt(&sh, centroid, action);
      //printf("nlist: %d\n", Nlist); 
      /* And prune the list by performing the check */
      int len_ker_inside = kernel_inps_len(ke_inside, l_edge+Nlist);
      real_t ker_in_inside[len_ker_inside + dim*Nlist];
      // I just do it myself...
      for(i=0;i<len_ker_centroid;i++) ker_in_inside[i] = ker_in[i];
      for( i=0; i<Nlist;i++) {
	for( j=0; j<dim; j++) {
	  ker_in_inside[len_ker_centroid + dim*i + j] = x[ dim*list[i] + j];
	}
      }
      /* for(i=0;i<len_ker_centroid+dim*Nlist;i++) printf("%lf ",ker_in_inside[i]); printf("\n"); */
      real_t testeval[Nlist];
      for(i=0;i<Nlist;i++) testeval[i]=0.0;
      ke_inside->eval(l_edge+Nlist, ker_in_inside, testeval);
      //for(i=0;i<Nlist;i++) printf("%lf ", testeval[i]);printf("\n"); 
      
      /* Push the result to the graph */
      int Naccept=0;
      hypervertex_t newedge[l_edge+Nlist];
      for(i=0;i<l_edge;i++) newedge[i] = edge[i];
      for(i=0;i<Nlist;i++) {
	if (testeval[i]>0 && found[list[i]]==0) {
	  newedge[l_edge+Naccept] = (PV? PV[list[i]] : list[i]);
	  found[list[i]]=1;
	  Naccept++;
	}
      }

      //printf("Naccept %d\n",Naccept);
      if(Naccept > 0) {
	Hypergraph_Push_Edge(hgnew,  l_edge+Naccept, newedge);
      }
    } // end edge loop
  } // end subgraph loop
  free(list);

  /* Clean up */
  SpatialHash_destroy(&sh);
  Target_Destroy(att+0);
  free(found);
}
