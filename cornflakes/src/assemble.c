#include "assemble.h"

#include <stdio.h>

void collect(real_t * ker_in, kernel_t * ke, hypervertex_t* edge, int l_edge,
	     dofmap_t ** dms, real_t ** data)
{
  hypervertex_t V;
  int i, j, k, fnum, maxlen;
  dofmap_t * dmap;
  real_t * datum;
  real_t * ker_in_iter = ker_in;
  // TODO: Fringe case of global and variable length data.
  for(i=0;i<ke->ninp;i++) {
    fnum = ke->inp[i].field_number;
    dmap = dms[i];
    datum = data[i];
    maxlen = Dofmap_Max_Len(dmap);
    int dofs[maxlen];
    int ndof;
    for(j=0;j<l_edge;j++) {
      V = edge[j];
      Dofmap_Get(dmap, V, dofs,&ndof);
      for(k=0;k<ndof;k++) ker_in_iter[k] = datum[ dofs[k] ];
      ker_in_iter += ndof;
    }
  }

}


real_t * push_target(assemble_target_t * att, int len_loc_out, int hx, int * outmap, real_t * loc_out)
{
  
  int i,j;
  switch(att->rank) {
  case 1:
    for(i=0;i<len_loc_out;i++) {
      att->V[ outmap[i] ] += loc_out[i];
    }
    return loc_out + len_loc_out;
    break;
  case 2:
    for(i=0;i<len_loc_out;i++) {
      for(j=0;j<len_loc_out;j++) {
	att->II[ len_loc_out*len_loc_out*hx+len_loc_out*i + j] = outmap[ i ];
	att->JJ[ len_loc_out*len_loc_out*hx+len_loc_out*i + j] = outmap[ j ];
	att->V [ len_loc_out*len_loc_out*hx+len_loc_out*i + j] = loc_out[len_loc_out*i + j];
      }
    }
    return loc_out + len_loc_out*len_loc_out;
    break;
  default:
    att->V[0] += loc_out[0];
    return loc_out + 1;
    break;
  }
  
}

void place_targets(assemble_targets_t * att, real_t * ker_out, int len_ker_out,
		   hypervertex_t * edge, int l_edge)
{
  
}

void assemble_targets(kernel_t * ke, hypergraph_t * hg,
		      dofmap_t ** dofmaps, real_t ** data,
		      int * outmap, // TODO KILL
		      assemble_target_t * att)
{
  int i,j, hex,hx;
  hyperedges_t * he;
  /* Loop over the graph sets */
  for(he = hg->he; he < hg->he+hg->n_types ; he++) {
    /* Allocate the loca vectors for this size edge */
    int len_ker_in = kernel_inp_len(ke, he->l_edge);
    int len_ker_out = kernel_outp_len(ke, he->l_edge); // RIGHT NOW SAME DOFMAP FOR ALL OF THEM
    real_t ker_in[ len_ker_in];
    real_t ker_out[len_ker_out];

    /* Loop over the edges */
    hypervertex_t * edge;
    for(hex=0; hex<he->n_edge; hex++) {
      edge = Hyperedges_Get_Edge(he, hex);

      /* Collect the data */
      collect(ker_in, ke, edge,he->l_edge, dofmaps,data); // TODO: Optimize by moving some overheard outside of loop

      /* Calculate the kernel */
      ke->eval(he->l_edge, ker_in, ker_out);

      /* Push the data */
      
    }
  }

  
}




void assemble_targets_dep(int ntarget, assemble_target_t * att,
			  kernel_t * ke, hyperedges_t * hg,
			  int * outmap, real_t ** data);
{
  #if 0
  int i,j,hx;
  /* II,JJ,KK better have size len_loc_out^2 * he->n_edge */
  /*
   * Allocate local vectors
   */
  int len_loc_in = kernel_inp_len(ke, he->l_edge);
  int len_loc_out = kernel_outp_len(ke, he->l_edge);
  real_t loc_in[len_loc_in];
  int out_alloc = 0; // THE KERNEL SHOULD ENCODE THIS INFORMATION....
  for(i=0;i<ntarget;i++) {
    switch(att[i].rank) {
    case 0: out_alloc += 1; break;
    case 1: out_alloc += len_loc_out; break;
    case 2: out_alloc += len_loc_out*len_loc_out; break; // FOR NOW: ONLY SQUARE MATRICES
    }
  }
  real_t loc_out[out_alloc];

  /*
   * Loop over the edges in the graph
   */
  hypervertex_t * edge;
  for(hx=0; hx<he->n_edge; hx++) {
    edge = Hyperedges_Get_Edge(he,hx);
    //for(i=0;i<he->l_edge;i++) printf("%d ",edge[i]); printf("\n");

    /* Pull data */
    collect(loc_in, ke, hx,edge,he->l_edge, data);

    //for(i=0;i<len_loc_in;i++) printf("%lf ",loc_in[i]); printf("\n");
    /* Call the kernel */
    ke->assem( loc_in, loc_out );
    //for(i=0;i<len_loc_out;i++) printf("%lf ",loc_out[i]); printf("\n");
    real_t * iter_out = loc_out;
    for(i=0;i<ntarget;i++) {
      iter_out = push_target(att+i,len_loc_out, hx,outmap+len_loc_out*hx, iter_out);
    }
  }
  #endif
}




/************************************
 ***         DEPRECATED:          ***
 ************************************/


void assemble_vector(real_t * R, kernel_t * ke,hyperedges_t * he, int * outmap, real_t ** data)
{
#if 0
  int hx,i,j,A;
  /*
   * R will be accumalated onto!
   */

  /*
   * Allocate local vectors
   */
  int len_loc_in = kernel_inp_len(ke, he->l_edge);
  int len_loc_out = kernel_outp_len(ke, he->l_edge);
  real_t loc_in[len_loc_in]; 
  real_t loc_out[len_loc_out];

  //printf("local sizes are %d %d\n",len_loc_in,len_loc_out);

  //printf("%lf %lf\n", data[0][10],data[1][10]);
  /*
   * Loop over the edges in the graph
   */
  int * edge;
  for(hx=0; hx<he->n_edge; hx++) {
    edge = Hyperedges_Get_Edge(he,hx);
    //for(i=0;i<he->l_edge;i++) printf("%d ",edge[i]); printf("\n");

    /* Pull data */
    collect(loc_in, ke, hx,edge,he->l_edge, data);

    //for(i=0;i<len_loc_in;i++) printf("%lf ",loc_in[i]); printf("\n");
    /* Call the kernel */
    ke->assem( loc_in, loc_out );
    
    /* Push the data */
    for(i=0;i<len_loc_out;i++) {
      R[ outmap[len_loc_out*hx + i] ] += loc_out[i];
    }
  }
  #endif
}


void assemble_matrix(int * II, int * JJ, real_t * KK,
		     kernel_t * ke, hyperedges_t * he,
		     int * outmap, real_t ** data)
{
  #if 0
  int i,j,hx;
  /* II,JJ,KK better have size len_loc_out^2 * he->n_edge */
  /*
   * Allocate local vectors
   */
  int len_loc_in = kernel_inp_len(ke, he->l_edge);
  int len_loc_out = kernel_outp_len(ke, he->l_edge);
  real_t loc_in[len_loc_in]; 
  // Don't need this guy, we just plop it into KK
  //real_t loc_out[len_loc_out*len_loc_out]; // FOR NOW: ONLY SQUARE MATRICES

  /*
   * Loop over the edges in the graph
   */
  hypervertex_t * edge;
  
  for(hx=0; hx<he->n_edge; hx++) {
    edge = Hyperedges_Get_Edge(he,hx);
    //for(i=0;i<he->l_edge;i++) printf("%d ",edge[i]); printf("\n");

    /* Pull data */
    collect(loc_in, ke, hx,edge,he->l_edge, data);

    //for(i=0;i<len_loc_in;i++) printf("%lf ",loc_in[i]); printf("\n");
    /* Call the kernel */
    ke->assem( loc_in, KK+len_loc_out*len_loc_out*hx );
    
    /* Push the II,JJ indices */
    for(i=0;i<len_loc_out;i++) {
      for(j=0;j<len_loc_out;j++) {
	II[ len_loc_out*len_loc_out*hx + len_loc_out*i + j] = outmap[len_loc_out*hx + i ];
	JJ[ len_loc_out*len_loc_out*hx + len_loc_out*i + j] = outmap[len_loc_out*hx + j ];
      }
    }
  }

 #endif 
}


void assemble_vector_matrix(real_t * R,
			    int * II, int * JJ, real_t * KK,
			    kernel_t * ke, hyperedges_t * he,
			    int * outmap, real_t ** data)
{
  #if 0
    int i,j,hx;
  /* II,JJ,KK better have size len_loc_out^2 * he->n_edge */
  /*
   * Allocate local vectors
   */
  int len_loc_in = kernel_inp_len(ke, he->l_edge);
  int len_loc_out = kernel_outp_len(ke, he->l_edge);
  real_t loc_in[len_loc_in]; 
  real_t loc_out[len_loc_out + len_loc_out*len_loc_out]; // FOR NOW: ONLY SQUARE MATRICES

  /*
   * Loop over the edges in the graph
   */
  hypervertex_t * edge;
  
  for(hx=0; hx<he->n_edge; hx++) {
    edge = Hyperedges_Get_Edge(he,hx);
    //for(i=0;i<he->l_edge;i++) printf("%d ",edge[i]); printf("\n");

    /* Pull data */
    collect(loc_in, ke, hx,edge,he->l_edge, data);

    //for(i=0;i<len_loc_in;i++) printf("%lf ",loc_in[i]); printf("\n");
    /* Call the kernel */
    ke->assem( loc_in, loc_out );
    
    /* Push the II,JJ indices */
    for(i=0;i<len_loc_out;i++) {
      for(j=0;j<len_loc_out;j++) {
	II[ len_loc_out*len_loc_out*hx + len_loc_out*i + j] = outmap[len_loc_out*hx + i ];
	JJ[ len_loc_out*len_loc_out*hx + len_loc_out*i + j] = outmap[len_loc_out*hx + j ];
	KK[len_loc_out*len_loc_out*hx + len_loc_out*i + j] = loc_out[len_loc_out + len_loc_out*i + j];
      }
      R[outmap[len_loc_out*hx + i]] += loc_out[i];
    }
    
  }
#endif
}

