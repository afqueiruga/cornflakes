#include "graphers.h"

#include "stdlib.h"
#include "math.h"

void Build_Pair_Graph(hypergraph_t * hg, int Npart, int dim, real_t * x, real_t cutoff) {
  spatialhash_t sh;
  //printf("alloc hg\n");
  Hypergraph_Alloc(hg,1); //2, Npart);
  //printf("alloc hash\n");
  Build_New_Hash(&sh, Npart,dim,x, cutoff);
  //printf("scan\n");
  void action(int a, int b) {
    if(dist(dim, x+dim*a,x+dim*b)<=cutoff) {
      int v[2] = {a,b};
      Hypergraph_Push_Edge(hg,2,v);
    }
  }
  SpatialHash_Scanall(&sh,x,action);
  //printf("destroy hash\n");
  SpatialHash_destroy(&sh);
}

void Build_Adjacency_Graph_Uniform(hypergraph_t * hg, int Npart, int dim, real_t * x, real_t cutoff) {
  spatialhash_t sh;
  int A;
  //printf("alloc hg\n");
  Hypergraph_Alloc(hg,1); //2, Npart);
  //printf("alloc hash\n");
  Build_New_Hash(&sh, Npart,dim,x, cutoff);
  //printf("scan\n");

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
  //SpatialHash_Scanall(&sh,x,action);
  //printf("destroy hash\n");
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
void Tie_Cells_and_Particles(hypergraph_t * hgnew,
			     hypergraph_t * mesh,
			     kernel_t * ke_circum,
			     kernel_t * ke_centroid,
			     kernel_t * ke_inside,
			     dofmap_t ** dofmaps,
			     real_t ** data,
			     int Npart, int dim, real_t * x,
			     hypervertex_t * PV)
{
  int i, j;
  int n_cells = mesh->he[0].n_edge; // TODO: BUG. Hypergraph needs a total_edge_count method

  /* Setup storage for the results */
  assemble_target_t att[1];
  setup_targets(ke_circum,att, mesh,n_cells);
  
  /* Calculate the maximum circum radius */
  assemble_targets(ke_circum,mesh, dofmaps,data, att);
  real_t hmax = 0.0;
  for(int i=0; i<n_cells; i++) {
    if(att[0].V[i] > hmax) hmax = att[0].V[i];
  }
  hmax = sqrt(hmax);
  //printf("hmax: %lf\n",hmax);
  /* Build a hash using this circumradius and fill it with the particles */
  spatialhash_t sh;
  Build_New_Hash(&sh, Npart,dim,x, hmax);
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
  destroy_targets(ke_circum,att);
  free(found);
}
